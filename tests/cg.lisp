(in-package #:coalton-linalg/tests)

;;; Conjugate gradients. Three kinds of evidence, none of which is the library's own word: the
;;; residual of every returned vector is formed here in exact rational arithmetic; the error
;;; of every iterate, measured exactly in the norm the method minimizes, is shown not to rise;
;;; and the whole iteration is recomputed here and compared bit for bit.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "CG" "COALTON-LINALG/CONJUGATE-GRADIENT")
  (sb-ext:add-package-local-nickname "CG-CALLER" "COALTON-LINALG/TESTS/CG-CALLERS"))

(defun two-norm (vector)
  "The 2-norm as the library forms it: a sum of squares in index order, then a root."
  (let ((sum 0d0))
    (declare (type double-float sum))
    (loop :for v :across vector :do (incf sum (* v v)))
    (sqrt sum)))

(defun exact-residual-squared (dense b x)
  "The squared 2-norm of b - A x, exactly, from the doubles of A, b and x."
  (let ((n (length b)))
    (loop :for i :below n
          :sum (expt (- (rational (aref b i))
                        (loop :for j :below n :sum (* (rational (coerce (aref dense i j) 'double-float)) (rational (aref x j)))))
                     2))))

(defun exact-energy-error (dense solution x)
  "(x - x*)^T A (x - x*), exactly: the squared error in the norm conjugate gradients minimize."
  (let* ((n (length x)) (e (make-array n)))
    (dotimes (i n) (setf (aref e i) (- (rational (aref x i)) (aref solution i))))
    (loop :for i :below n
          :sum (* (aref e i) (loop :for j :below n :sum (* (rational (coerce (aref dense i j) 'double-float)) (aref e j)))))))

(defun reference-cg (matrix b start precondition limit tolerance)
  "The iteration of docs/conjugate-gradient.md recomputed here with the library's product: the iterates, the running
norms and the number of steps."
  (let* ((n (length b)) (x (copy-seq start)) (q (zeros n)) (r (zeros n)) (z (zeros n)) (p (zeros n))
         (iterates (list (copy-seq x))) (norms '()) (rho 0d0) (fresh t) (k 0) (sum 0d0))
    (declare (type double-float rho sum))
    (flet ((form ()
             (csc:matvec! matrix x q)
             (setf sum 0d0)
             (dotimes (i n)
               (let ((v (- (aref b i) (aref q i))))
                 (setf (aref r i) v)
                 (incf sum (* v v)))))
           (dot (a c)
             (let ((s 0d0))
               (declare (type double-float s))
               (dotimes (i n) (incf s (* (aref a i) (aref c i))))
               s)))
      (form)
      (push (sqrt sum) norms)
      (loop :while (< k limit) :do
        (when fresh
          (when (<= (sqrt sum) tolerance) (return))
          (funcall precondition r z)
          (setf rho (dot r z))
          (replace p z)
          (setf fresh nil))
        (csc:matvec! matrix p q)
        (let ((alpha (/ rho (dot p q))))
          (setf sum 0d0)
          (dotimes (i n)
            (setf (aref x i) (+ (aref x i) (* alpha (aref p i))))
            (let ((v (- (aref r i) (* alpha (aref q i)))))
              (setf (aref r i) v)
              (incf sum (* v v))))
          (incf k)
          (push (sqrt sum) norms)
          (push (copy-seq x) iterates)
          (cond ((<= (sqrt sum) tolerance)
                 (form)
                 (setf fresh t))
                (t
                 (funcall precondition r z)
                 (let* ((next (dot r z)) (beta (/ next rho)))
                   (setf rho next)
                   (dotimes (i n) (setf (aref p i) (+ (aref z i) (* beta (aref p i)))))))))))
    (values (nreverse iterates) (nreverse norms) k)))

(defun check-cg (label dense matrix b start tolerance &key precondition reference-precondition (limit 200))
  "Solve and hold the outcome against exact arithmetic and against the recomputed iteration. Returns the steps taken."
  (let* ((n (length b)) (workspace (cg:make-cg n limit)) (x (copy-seq start))
         (apply (rf-caller:matrix-operator matrix))
         (code (if precondition
                   (cg-caller:pcg-code workspace apply precondition b x tolerance)
                   (cg-caller:cg-code workspace apply b x tolerance)))
         (steps (cg:step-count workspace)))
    (assert (= 0 code) () "~A: outcome ~D after ~D steps, running norms ~S." label code steps
            (subseq (cg:recurrence-norms workspace) 0 (1+ steps)))
    (assert (and (= n (cg:dimension workspace)) (= limit (cg:step-limit workspace))))
    ;; The reported norm is that of the residual formed from the operator, and is within the
    ;; tolerance; the exact residual of the vector returned agrees with it to rounding.
    (let* ((reported (cg:residual-norm workspace))
           (exact-square (exact-residual-squared dense b x))
           (rounding (* 64 n double-float-epsilon
                        (loop :for i :below n
                              :maximize (+ (abs (aref b i))
                                           (loop :for j :below n :sum (abs (* (aref dense i j) (aref x j)))))))))
      (assert (<= 0d0 reported tolerance) () "~A: converged with a reported residual of ~S against ~S." label reported
              tolerance)
      (assert (<= (expt (max 0 (- (rational reported) (rational rounding))) 2)
                  exact-square
                  (expt (+ (rational reported) (rational rounding)) 2))
              () "~A: reported residual ~S, exact squared residual ~S." label reported exact-square))
    ;; The iteration recomputed: every running norm and the returned vector, bit for bit.
    (multiple-value-bind (iterates norms reference-steps)
        (reference-cg matrix b start (or reference-precondition (lambda (r z) (replace z r))) limit tolerance)
      (assert (= steps reference-steps) () "~A: ~D steps, recomputed ~D." label steps reference-steps)
      (assert (same-bits-p (subseq (cg:recurrence-norms workspace) 0 (1+ steps))
                           (coerce norms '(simple-array double-float (*))))
              () "~A: running norms ~S, recomputed ~S." label (subseq (cg:recurrence-norms workspace) 0 (1+ steps)) norms)
      (assert (same-bits-p x (car (last iterates))) () "~A: the returned vector is not the recomputed one." label)
      ;; What the method is: each iterate is at least as near the solution as the one before,
      ;; in the norm of the operator, exactly measured. A part in 2^30 is allowed for rounding.
      (let ((solution (exact-solution dense b)))
        (loop :for (earlier later) :on (mapcar (lambda (iterate) (exact-energy-error dense solution iterate)) iterates)
              :while later
              :do (assert (<= later (* earlier (+ 1 (expt 2 -30)))) ()
                          "~A: the error in the operator's norm rose from ~E to ~E." label (float earlier 1d0)
                          (float later 1d0)))))
    steps))

(defun run-cg-convergence-tests ()
  (let ((state (sb-ext:seed-random-state 20261206)) (cases 0) (most 0) (preconditioned-total 0) (plain-total 0))
    (dolist (n '(1 2 3 5 8 12))
      (dotimes (repeat 4)
        (multiple-value-bind (dense matrix symbolic b) (test-system state n 0)
          (let* ((factor (ldlt:symbolic-factor symbolic))
                 (tolerance (* (scale-float 1d0 -40) (max 1d0 (two-norm b))))
                 (diagonal (make-array n :element-type 'double-float))
                 (start (zeros n)) (elsewhere (zeros n)))
            (dotimes (i n)
              (setf (aref diagonal i) (coerce (aref dense i i) 'double-float))
              (setf (aref elsewhere i) (coerce (- (random 9 state) 4) 'double-float)))
            (let ((plain (check-cg (format nil "plain ~D #~D" n repeat) dense matrix b start tolerance)))
              ;; In exact arithmetic the method ends in at most n steps; rounding may add a few.
              (assert (<= plain (+ n 3)) () "Order ~D took ~D steps." n plain)
              (setf most (max most plain))
              (incf plain-total plain))
            (check-cg (format nil "from elsewhere ~D #~D" n repeat) dense matrix b elsewhere tolerance)
            ;; Preconditioned by the diagonal, and by a factor of the matrix itself, which ends it at once.
            (incf preconditioned-total
                  (check-cg (format nil "diagonal ~D #~D" n repeat) dense matrix b start tolerance
                            :precondition (cg-caller:diagonal-preconditioner diagonal)
                            :reference-precondition (lambda (r z)
                                                      (dotimes (i n) (setf (aref z i) (/ (aref r i) (aref diagonal i)))))))
            (let ((exact-steps (check-cg (format nil "factor ~D #~D" n repeat) dense matrix b start tolerance
                                         :precondition (rf-caller:factor-solver factor)
                                         :reference-precondition (lambda (r z) (ldlt:solve! factor r z)))))
              (assert (<= exact-steps 2) () "Preconditioned by a factor of the matrix: ~D steps." exact-steps))
            (incf cases 4)))))
    (format t "CG-CONVERGENCE-PASS cases=~D most-plain-steps=~D plain-steps=~D diagonal-steps=~D~%" cases most plain-total
            preconditioned-total)))

(defun run-cg-outcome-tests ()
  (let* ((nan (not-a-number))
         (dense (rows-matrix '((4 1 0) (1 3 1) (0 1 2)))) (matrix (full-csc dense))
         (apply (rf-caller:matrix-operator matrix))
         (b (doubles 1d0 2d0 3d0)))
    (flet ((run (workspace operator rhs x tolerance) (cg-caller:cg-code workspace operator rhs x tolerance)))
      ;; Before any call.
      (let ((workspace (cg:make-cg 3 10)))
        (assert (and (= 0 (cg:step-count workspace)) (= 0 (cg:restart-count workspace))
                     (= -1d0 (cg:residual-norm workspace))))
        ;; Sizes and tolerances are refused before anything is read or written.
        (let ((x (doubles 9d0 9d0 9d0)))
          (assert (= 5 (run workspace apply (doubles 1d0 2d0) x 1d-9)))
          (assert (= 5 (run workspace apply b (doubles 9d0 9d0) 1d-9)))
          (assert (= 6 (run workspace apply b x -1d0)))
          (assert (= 6 (run workspace apply b x nan)))
          (assert (same-bits-p x (doubles 9d0 9d0 9d0)))
          (assert (= -1d0 (cg:residual-norm workspace))))
        ;; A right-hand side or a start that is not finite.
        (assert (= 3 (run workspace apply (doubles 1d0 nan 3d0) (zeros 3) 1d-9)))
        ;; The library's own product refuses such a vector by a typed error, so this operator checks nothing.
        (assert (= 3 (run workspace (cg-caller:scaling-operator 2d0) b
                          (doubles 0d0 sb-ext:double-float-positive-infinity 0d0) 1d-9)))
        (assert (= -1d0 (cg:residual-norm workspace)))
        ;; A sound call, then a refused one: the earlier report is withdrawn.
        (assert (= 0 (run workspace apply b (zeros 3) 1d-12)))
        (assert (and (plusp (cg:step-count workspace)) (<= 0d0 (cg:residual-norm workspace))))
        (assert (= 6 (run workspace apply b (zeros 3) -1d0)))
        (assert (and (= 0 (cg:step-count workspace)) (= -1d0 (cg:residual-norm workspace)))))
      ;; A limit of zero judges the start as it stands and leaves it untouched.
      (let ((workspace (cg:make-cg 3 0)) (x (doubles 1d0 1d0 1d0)))
        (assert (= 1 (run workspace apply b x 1d-12)))
        (assert (same-bits-p x (doubles 1d0 1d0 1d0)))
        (assert (= 0 (cg:step-count workspace)))
        (assert (= (cg:residual-norm workspace) (two-norm (doubles -4d0 -3d0 0d0))))
        (assert (= 0 (run workspace apply b x 10d0)))
        ;; With no step to take, a right-hand side that is not a number is still not a residual.
        (assert (= 3 (run workspace apply (doubles 1d0 nan 3d0) x 1d-12)))
        (assert (= -1d0 (cg:residual-norm workspace))))
      ;; A limit met: one step taken, and the report is the residual of what is returned.
      (let ((workspace (cg:make-cg 3 1)) (x (zeros 3)))
        (assert (= 1 (run workspace apply b x 1d-12)))
        (assert (= 1 (cg:step-count workspace)))
        (let ((exact (sqrt (float (exact-residual-squared dense b x) 1d0))))
          (assert (<= (abs (- (cg:residual-norm workspace) exact)) 1d-13) () "At the limit: reported ~S, exact ~S."
                  (cg:residual-norm workspace) exact)
          (assert (> exact 1d-3))))
      ;; A tolerance of zero is allowed and is met where the arithmetic is exact: the identity in one step.
      (let ((workspace (cg:make-cg 2 5)) (x (zeros 2))
            (identity-operator (rf-caller:matrix-operator (csc-from-entries 2 2 '((0 0 1) (1 1 1))))))
        (assert (= 0 (run workspace identity-operator (doubles 3d0 -4d0) x 0d0)))
        (assert (and (= 1 (cg:step-count workspace)) (= 0d0 (cg:residual-norm workspace))))
        (assert (same-bits-p x (doubles 3d0 -4d0)))
        ;; The same solve with one step allowed: the step that converges is the last one.
        (let ((tight (cg:make-cg 2 1)) (y (zeros 2)))
          (assert (= 0 (run tight identity-operator (doubles 3d0 -4d0) y 0d0)))
          (assert (= 1 (cg:step-count tight))))
        ;; A start that already solves the system takes no step.
        (assert (= 0 (run workspace identity-operator (doubles 3d0 -4d0) x 0d0)))
        (assert (= 0 (cg:step-count workspace)))
        ;; A zero right-hand side from a zero start is solved as it stands.
        (assert (= 0 (run workspace identity-operator (zeros 2) (zeros 2) 0d0))))
      ;; Operators that are not positive definite: a direction of zero curvature, one of negative
      ;; curvature, and the zero operator. The start is returned where the first direction fails.
      (dolist (entries '(((0 0 1) (1 1 -1)) ((0 0 1) (1 1 -2)) ()))
        (let ((workspace (cg:make-cg 2 5)) (x (zeros 2))
              (operator (rf-caller:matrix-operator (csc-from-entries 2 2 entries))))
          (assert (= 2 (run workspace operator (doubles 1d0 1d0) x 1d-9)) () "Entries ~S were taken as definite." entries)
          (assert (same-bits-p x (zeros 2)))
          (assert (and (= 0 (cg:step-count workspace)) (= -1d0 (cg:residual-norm workspace))))))
      ;; An operator whose product overflows, in the first direction's image.
      (let ((workspace (cg:make-cg 1 5)))
        (assert (= 3 (run workspace (cg-caller:scaling-operator 1d200) (doubles 1d200) (zeros 1) 1d-9)))
        (assert (= -1d0 (cg:residual-norm workspace))))
      ;; A curvature that is not a number: the first direction's image is spoiled. It is not
      ;; reported as an operator that is not positive.
      (let ((workspace (cg:make-cg 3 10)) (calls (indices 0)))
        (assert (= 3 (run workspace (cg-caller:drifting-operator apply 1 nan calls) b (zeros 3) 1d-12))))
      ;; A step length that overflows, from an operator of 1e-320 under a right-hand side of
      ;; 1e10: reported before the vector is touched.
      (let ((workspace (cg:make-cg 1 5)) (x (zeros 1)))
        (assert (= 3 (run workspace (cg-caller:scaling-operator 1d-320) (doubles 1d10) x 1d-9)))
        (assert (same-bits-p x (zeros 1)) () "A step length that is not finite was applied: ~S." x))
      ;; A running residual whose norm overflows while the curvature and the step length are
      ;; finite: the operator is diag(1e-200, 1e200) and the right-hand side (1e100, 1e-100).
      (let ((workspace (cg:make-cg 2 5)))
        (assert (= 3 (run workspace (cg-caller:diagonal-operator (doubles 1d-200 1d200)) (doubles 1d100 1d-100)
                          (zeros 2) 1d-9)))
        (assert (and (= 0 (cg:step-count workspace)) (= -1d0 (cg:residual-norm workspace)))))
      ;; Preconditioners: one that is not positive, one that refuses at once, one that refuses later.
      (let ((workspace (cg:make-cg 3 10)) (x (zeros 3)) (calls (indices 0)))
        (assert (= 2 (cg-caller:pcg-code workspace apply (cg-caller:scaled-preconditioner -1d0) b x 1d-12)))
        (assert (same-bits-p x (zeros 3)))
        (assert (= 2 (cg-caller:pcg-code workspace apply (cg-caller:scaled-preconditioner 0d0) b x 1d-12)))
        (assert (= 4 (cg-caller:pcg-code workspace apply
                                         (rf-caller:refusing-solver (cg-caller:scaled-preconditioner 1d0) 0 calls)
                                         b x 1d-12)))
        (assert (and (same-bits-p x (zeros 3)) (= 0 (cg:step-count workspace))))
        (setf (aref calls 0) 0)
        (assert (= 4 (cg-caller:pcg-code workspace apply
                                         (rf-caller:refusing-solver (cg-caller:scaled-preconditioner 1d0) 1 calls)
                                         b x 1d-12)))
        (assert (and (= 1 (cg:step-count workspace)) (= -1d0 (cg:residual-norm workspace))))
        ;; One that is positive on its first call and not on its second.
        (setf (aref calls 0) 0)
        (assert (= 2 (cg-caller:pcg-code workspace apply (cg-caller:switching-preconditioner calls) b (zeros 3) 1d-12)))
        (assert (= 1 (cg:step-count workspace)))
        ;; A positive multiple of the identity changes nothing but the scale of the directions.
        (let ((y (zeros 3)))
          (assert (= 0 (cg-caller:pcg-code workspace apply (cg-caller:scaled-preconditioner 0.25d0) b y 1d-12)))
          (assert (<= (cg:residual-norm workspace) 1d-12))))
      ;; One wrong product. The running residual then meets the tolerance while the residual of
      ;; the vector does not; the call must not report convergence on the running one, and
      ;; must go on from the residual it forms until that one meets the tolerance.
      (let* ((state (sb-ext:seed-random-state 77)))
        (multiple-value-bind (dense matrix symbolic b) (test-system state 8 0)
          (declare (ignore symbolic))
          (let* ((honest (cg:make-cg 8 200)) (drifting (cg:make-cg 8 200))
                 (operator (rf-caller:matrix-operator matrix))
                 (calls (indices 0))
                 (tolerance (* (scale-float 1d0 -40) (two-norm b)))
                 (x (zeros 8)) (y (zeros 8)))
            (assert (= 0 (run honest operator b x tolerance)))
            (assert (= 1 (cg:restart-count honest)))
            ;; Product 0 forms the first residual; product 2 is the second step's.
            (assert (= 0 (run drifting (cg-caller:drifting-operator operator 2 1.001d0 calls) b y tolerance)))
            (assert (>= (cg:restart-count drifting) 2) () "The wrong product was not noticed: ~D replacements."
                    (cg:restart-count drifting))
            (assert (> (cg:step-count drifting) (cg:step-count honest)))
            (assert (<= (sqrt (float (exact-residual-squared dense b y) 1d0)) (* tolerance (+ 1d0 1d-3))) ()
                    "After a wrong product the call reported convergence with a residual of ~E against ~E."
                    (sqrt (float (exact-residual-squared dense b y) 1d0)) tolerance)
            ;; Stopped at the limit two steps after the wrong product: the running residual and
            ;; the residual of the vector differ by a thousandth, and the report is the latter.
            (let ((short (cg:make-cg 8 3)) (w (zeros 8)))
              (setf (aref calls 0) 0)
              (assert (= 1 (run short (cg-caller:drifting-operator operator 2 1.001d0 calls) b w tolerance)))
              (let ((exact (sqrt (float (exact-residual-squared dense b w) 1d0)))
                    (running (aref (cg:recurrence-norms short) 3)))
                (assert (> (abs (- running exact)) (* 1d-6 exact)) ()
                        "The limit case's premise failed: running ~S, exact ~S." running exact)
                (assert (<= (abs (- (cg:residual-norm short) exact)) (* 1d-9 exact)) ()
                        "At the limit after a wrong product: reported ~S, exact ~S." (cg:residual-norm short)
                        exact))))))))
  (format t "CG-OUTCOMES-PASS~%"))

(defparameter *cg-public-control*
  "(declare inspect (coalton-linalg/conjugate-gradient:CgWorkspace -> UFix))
   (define (inspect workspace)
     (+ (coalton-linalg/conjugate-gradient:dimension workspace) (coalton-linalg/conjugate-gradient:step-count workspace)))
   (declare probe (coalton-linalg/conjugate-gradient:CgWorkspace -> UFix))
   (define (probe workspace) (inspect workspace))")

(defparameter *cg-encapsulation-probes*
  '(("cg-residual-access" .
     "(declare probe (coalton-linalg/conjugate-gradient:CgWorkspace -> arr:LispArray F64))
      (define (probe workspace) (.residual workspace))")
    ("cg-unwrapping" .
     "(declare probe (coalton-linalg/conjugate-gradient:CgWorkspace -> UFix))
      (define (probe workspace)
        (match workspace ((coalton-linalg/conjugate-gradient:MkCgWorkspace held) (.size held))))")
    ("cg-status-as-solve-status" .
     "(declare probe (coalton-linalg/conjugate-gradient:CgStatus -> ldlt:SolveStatus))
      (define (probe status) status)")))

(defun run-cg-lifetime-tests ()
  (let ((state (sb-ext:seed-random-state 78)))
    (multiple-value-bind (dense matrix symbolic b) (test-system state 8 1/4)
      (declare (ignore dense))
      (let* ((workspace (cg:make-cg 8 100)) (operator (rf-caller:matrix-operator matrix))
             (shifted (rf-caller:factor-solver (ldlt:symbolic-factor symbolic)))
             (tolerance (* (scale-float 1d0 -40) (two-norm b)))
             (start (zeros 8)) (x (zeros 8)) (first-answer (zeros 8)))
        ;; One workspace, call after call, plain and preconditioned: the same answers bit for bit.
        (assert (= 0 (cg-caller:cg-code workspace operator b first-answer tolerance)))
        (assert (= 0 (cg-caller:pcg-code workspace operator shifted b x tolerance)))
        (replace x start)
        (assert (= 0 (cg-caller:cg-code workspace operator b x tolerance)))
        (assert (same-bits-p x first-answer) () "The same call on a used workspace gave another vector.")
        ;; Steady state: compiled callers allocate nothing.
        (allocation-check "cg-loop"
                          (lambda () (assert (= 0 (cg-caller:cg-loop workspace operator b start x tolerance 1)))))
        (allocation-check "pcg-loop"
                          (lambda ()
                            (assert (= 0 (cg-caller:pcg-loop workspace operator shifted b start x tolerance 1))))))))
  (multiple-value-bind (accepted diagnostics)
      (coalton-linalg/tests/encapsulation::compile-probe "cg-public-control" *cg-public-control*)
    (assert accepted () "The valid conjugate gradient external caller failed:~%~A" diagnostics))
  (dolist (probe *cg-encapsulation-probes*)
    (assert (not (coalton-linalg/tests/encapsulation::compile-probe (car probe) (cdr probe))) ()
            "External caller ~A compiled despite required rejection." (car probe)))
  (dolist (name '("MKCGWORKSPACE" "CGSTORAGE"))
    (multiple-value-bind (symbol visibility) (find-symbol name "COALTON-LINALG/CONJUGATE-GRADIENT")
      (assert (and symbol (eq visibility ':internal)) ()
              "Storage constructor ~A must exist and remain internal; visibility is ~S." name visibility)))
  (format t "CG-LIFETIME-PASS rejected=~D public-control=1~%" (length *cg-encapsulation-probes*)))

(defun check-iterative-shared-rhs ()
  (let ((times-two (cg-caller:scaling-operator 2d0)) (half (cg-caller:scaled-preconditioner 0.5d0)))
    (dolist (method '(:cg :pcg :refine))
      (let ((v (doubles 2d0)))
        (assert (zerop (ecase method
                         (:cg (cg-caller:cg-code (cg:make-cg 1 2) times-two v v 0.5d0))
                         (:pcg (cg-caller:pcg-code (cg:make-cg 1 2) times-two half v v 0.5d0))
                         (:refine (rf-caller:refine-code (rf:make-refiner 1 2) times-two half v v 0.5d0 1d0)))))
        (assert (= 1d0 (aref v 0)))))))

(defun check-cg-tiny-residual ()
  (let ((op (cg-caller:scaling-operator 1d0)))
    (dolist (b (list (scale-float 1d0 -530) (scale-float 1d0 -600) least-positive-double-float))
      (let ((w (cg:make-cg 1 0)) (x (doubles 0d0)))
        (assert (= 1 (cg-caller:cg-code w op (doubles b) x 0d0)))
        (assert (= b (cg:residual-norm w))))
      (let ((w (cg:make-cg 1 2)) (x (doubles 0d0)))
        (assert (= 0 (cg-caller:cg-code w op (doubles b) x 0d0)))
        ;; Do not round the square to F64 before testing the tolerance.
        (assert (zerop (exact-residual-squared (rows-matrix '((1))) (doubles b) x)))
        (assert (= 1 (cg:step-count w))))))
  ;; Exercise the recurrence as well as a one-step solve on tiny data.
  (let* ((s (scale-float 1d0 -600)) (tol (scale-float s -45)) (b (doubles s (* 2 s))) (x (zeros 2))
         (w (cg:make-cg 2 4)) (op (cg-caller:diagonal-operator (doubles 1d0 2d0))))
    (assert (zerop (cg-caller:cg-code w op b x tol)))
    (assert (>= (cg:step-count w) 2))
    (assert (<= (exact-residual-squared (rows-matrix '((1 0) (0 2))) b x) (expt (rational tol) 2))))
  ;; Unlike the identity, the products in these step lengths have different bands.
  (let ((b (doubles (scale-float 1d0 -600))) (w (cg:make-cg 1 2)) (start (zeros 1)) (x (zeros 1)))
    (dolist (power '(300 -300))
      (let ((op (cg-caller:scaling-operator (scale-float 1d0 power))))
        (fill x 0d0)
        (assert (zerop (cg-caller:cg-code w op b x 0d0)))
        (assert (= (aref x 0) (scale-float (aref b 0) (- power))))
        (allocation-check "cg-tiny-products"
                          (lambda () (assert (zerop (cg-caller:cg-loop w op b start x 0d0 1)))))))))

(defun check-pcg-last-preconditioner ()
  (let* ((w (cg:make-cg 2 1)) (x (zeros 2)) (calls (indices 0))
         (pre (rf-caller:refusing-solver (cg-caller:scaled-preconditioner 1d0) 1 calls)))
    (assert (= 1 (cg-caller:pcg-code w (cg-caller:diagonal-operator (doubles 1d0 3d0))
                                      pre (doubles 1d0 1d0) x 0d0)))
    (assert (= 1 (aref calls 0)))
    (assert (equalp x (doubles 0.5d0 0.5d0)))
    (assert (= (sqrt 0.5d0) (cg:residual-norm w)))))

(defun run-cg-tests ()
  (check-iterative-shared-rhs)
  (check-cg-tiny-residual)
  (check-pcg-last-preconditioner)
  (run-cg-convergence-tests)
  (run-cg-outcome-tests)
  (run-cg-lifetime-tests)
  (format t "CG-TESTS-PASS~%"))
