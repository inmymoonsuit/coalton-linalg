(in-package #:coalton-linalg/tests)

;;; Damped least squares on the normal equations. The reference is the step in exact rational
;;; arithmetic, -(J^T J + mu D)^-1 J^T r from the stored doubles, by a dense elimination written
;;; for the tests.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "DM" "COALTON-LINALG/SPARSE/DAMPED")
  (sb-ext:add-package-local-nickname "DM-CALLER" "COALTON-LINALG/TESTS/DAMPED-CALLERS"))

(defun damped-scaled-matrix (m n entries power)
  "The stored matrix with every entry times 2^POWER."
  (csc-from-entries m n (mapcar (lambda (entry)
                                  (list (first entry) (second entry)
                                        (scale-float (coerce (third entry) 'double-float) power)))
                                entries)))

(defun damped-column-squares (m n entries)
  "The exact squared length of each column, from the stored doubles."
  (let ((a (dr-dense m n entries)))
    (loop :for j :below n :collect (loop :for i :below m :sum (expt (aref a i j) 2)))))

(defun damped-exact-step (m n entries damping mu r)
  "-(J^T J + mu D)^-1 J^T r in rational arithmetic. DAMPING is a list of rationals, one a column."
  (let* ((a (dr-dense m n entries))
         (normal (build-matrix n (lambda (p q)
                                   (+ (loop :for i :below m :sum (* (aref a i p) (aref a i q)))
                                      (if (= p q) (* (rational mu) (nth p damping)) 0)))))
         (inverse (rational-inverse normal))
         (gradient (loop :for p :below n :collect (loop :for i :below m :sum (* (aref a i p) (rational (aref r i)))))))
    (coerce (loop :for p :below n :collect (- (loop :for q :below n :sum (* (aref inverse p q) (nth q gradient)))))
            'vector)))

(defun damped-error (step exact)
  "The largest difference from the exact step, over the larger of 1 and the exact step's largest entry."
  (/ (loop :for value :across step :for reference :across exact :maximize (abs (- (rational value) reference)))
     (max 1 (loop :for reference :across exact :maximize (abs reference)))))

(defun run-damped-exact-tests ()
  (let ((state (sb-ext:seed-random-state 20261122)) (cases 0) (worst 0) (given 0) (columns 0))
    (loop :for (m n) :in '((3 2) (4 4) (5 3) (6 6) (8 5) (4 6) (2 5) (7 7)) :do
      (dotimes (repeat 6)
        (let* ((entries (seeded-entries state m n 55 :hub t))
               (matrix (csc-from-entries m n entries))
               (plan (dm:analyze matrix))
               (r (coerce (loop :repeat m :collect (coerce (- (random 11 state) 5) 'double-float))
                          '(simple-array double-float (*))))
               (d (coerce (loop :repeat n :collect (scale-float (coerce (1+ (random 7 state)) 'double-float)
                                                               (- (random 5 state) 2)))
                          '(simple-array double-float (*))))
               (step (zeros n)))
          (assert (= 0 (dm-caller:matrix-code plan matrix)))
          ;; A given diagonal.
          (assert (= 0 (dm-caller:damping-code plan d)))
          (dolist (mu (list 0.125d0 1d0 8d0 (scale-float 1d0 -10)))
            (assert (= 0 (dm-caller:factor-code plan mu)) () "~Dx~D #~D: factor at mu ~S refused." m n repeat mu)
            (assert (= 0 (dm-caller:step-code plan r step)))
            (let ((error (damped-error step (damped-exact-step m n entries (map 'list #'rational d) mu r))))
              (incf cases) (incf given)
              (setf worst (max worst error))
              (assert (<= error (expt 2 -30)) () "~Dx~D #~D, mu ~S, a given diagonal: the step is off by ~E." m n
                      repeat mu (float error 1d0))))
          ;; The columns' own squared lengths, plus a floor.
          (let ((floor (scale-float 1d0 -3)))
            (assert (= 0 (dm-caller:column-damping-code plan floor)))
            (dolist (mu (list 0.125d0 1d0 8d0))
              (assert (= 0 (dm-caller:factor-code plan mu)))
              (assert (= 0 (dm-caller:step-code plan r step)))
              (let ((error (damped-error step (damped-exact-step m n entries
                                                                 (mapcar (lambda (square) (+ square (rational floor)))
                                                                         (damped-column-squares m n entries))
                                                                 mu r))))
                (incf cases) (incf columns)
                (setf worst (max worst error))
                (assert (<= error (expt 2 -30)) () "~Dx~D #~D, mu ~S, the columns' own diagonal: the step is off by ~E."
                        m n repeat mu (float error 1d0))))))))
    ;; No damping at all, on a matrix of full column rank: the least squares step itself.
    (let* ((entries '((0 0 2) (1 0 1) (1 1 3) (2 1 1) (2 2 2) (3 0 1) (3 2 1)))
           (plan (dm:analyze (csc-from-entries 4 3 entries))) (step (zeros 3)) (r (doubles 1d0 -2d0 3d0 5d0)))
      (assert (= 0 (dm-caller:matrix-code plan (csc-from-entries 4 3 entries))))
      (assert (= 0 (dm-caller:damping-code plan (doubles 1d0 1d0 1d0))))
      (assert (= 0 (dm-caller:factor-code plan 0d0)))
      (assert (= 0 (dm-caller:step-code plan r step)))
      (let ((error (damped-error step (damped-exact-step 4 3 entries '(1 1 1) 0 r))))
        (incf cases)
        (assert (<= error (expt 2 -40)) () "Without damping the step is off by ~E." (float error 1d0))))
    (format t "DAMPED-EXACT-PASS cases=~D given=~D columns=~D worst=~,3E~%" cases given columns (float worst 1d0))))

(defun run-damped-scale-tests ()
  "The step of a matrix and a residual at other scales is the step at scale 1 times the residual's power over
the matrix's, bit for bit: nothing is squared before it is divided."
  (let* ((entries '((0 0 4) (0 1 1) (1 1 3) (1 2 -2) (2 0 1) (2 2 5) (3 1 2) (3 2 1)))
         (m 4) (n 3) (r (doubles 1d0 -2d0 3d0 5d0)) (d (doubles 2d0 0.5d0 3d0)) (cases 0))
    (flet ((step-at (matrix-power residual-power mode)
             (let* ((matrix (damped-scaled-matrix m n entries matrix-power)) (plan (dm:analyze matrix)) (step (zeros n)))
               (assert (= 0 (dm-caller:matrix-code plan matrix)))
               (ecase mode
                 (:given (assert (= 0 (dm-caller:damping-code
                                       plan (map '(simple-array double-float (*))
                                                 (lambda (value) (scale-float value (* 2 matrix-power))) d)))))
                 (:columns (assert (= 0 (dm-caller:column-damping-code plan 0d0)))))
               (assert (= 0 (dm-caller:factor-code plan 0.25d0)))
               (assert (= 0 (dm-caller:step-code plan (map '(simple-array double-float (*))
                                                           (lambda (value) (scale-float value residual-power)) r)
                                                 step))
                       () "At 2^~D and 2^~D the step was refused." matrix-power residual-power)
               step)))
      (dolist (mode '(:given :columns))
        (let ((reference (step-at 0 0 mode)))
          (loop :for (matrix-power residual-power)
                  :in (if (eq mode ':given)
                          '((300 0) (-300 0) (0 600) (0 -600) (400 -400) (-400 400) (500 500) (-500 -500))
                          '((600 0) (-600 0) (0 600) (0 -600) (900 -100) (-500 500) (1000 1000) (-1000 -1000)
                            (1020 0) (-1020 0)))
                :do (let ((step (step-at matrix-power residual-power mode))
                          (expected (map '(simple-array double-float (*))
                                         (lambda (value) (scale-float value (- residual-power matrix-power))) reference)))
                      (incf cases)
                      (assert (same-bits-p step expected) ()
                              "~A, matrix at 2^~D and residual at 2^~D: the step is ~S and the scaled step ~S." mode
                              matrix-power residual-power step expected))))))
    (format t "DAMPED-SCALE-PASS cases=~D~%" cases)))

(defun run-damped-outcome-tests ()
  (let* ((entries '((0 0 4) (0 1 1) (1 1 3) (1 2 -2) (2 0 1) (2 2 5) (3 1 2) (3 2 1)))
         (matrix (csc-from-entries 4 3 entries))
         (plan (dm:analyze matrix))
         (r (doubles 1d0 -2d0 3d0 5d0)) (out (doubles 9d0 9d0 9d0)) (untouched (doubles 9d0 9d0 9d0))
         (nan (not-a-number)) (infinity sb-ext:double-float-positive-infinity))
    (assert (and (= 4 (dm:row-count plan)) (= 3 (dm:column-count plan))))
    ;; Before anything is set.
    (assert (not (or (dm:matrix-set? plan) (dm:damping-set? plan) (dm:factor-ready? plan))))
    (assert (= 4 (dm-caller:factor-code plan 1d0)))
    (assert (= 1 (dm-caller:step-code plan r out)))
    (assert (= 3 (dm-caller:column-damping-code plan 1d0)))
    ;; A damping before a matrix is taken and kept; still nothing to factor.
    (assert (= 0 (dm-caller:damping-code plan (doubles 1d0 2d0 3d0))))
    (assert (and (dm:damping-set? plan) (= 4 (dm-caller:factor-code plan 1d0))))
    ;; Dampings that cannot be used leave the standing one.
    (assert (= 1 (dm-caller:damping-code plan (doubles 1d0 2d0))))
    (dolist (bad (list 0d0 -1d0 nan infinity))
      (assert (= 2 (dm-caller:damping-code plan (doubles 1d0 bad 3d0))) () "A damping entry of ~S was taken." bad))
    (assert (= 0 (dm-caller:matrix-code plan matrix)))
    (dolist (bad (list -1d0 nan infinity))
      (assert (= 2 (dm-caller:column-damping-code plan bad)) () "A floor of ~S was taken." bad))
    ;; The standing damping is still (1, 2, 3): the step says so.
    (assert (= 0 (dm-caller:factor-code plan 1d0)))
    (assert (dm:factor-ready? plan))
    (let ((step (zeros 3)))
      (assert (= 0 (dm-caller:step-code plan r step)))
      (assert (<= (damped-error step (damped-exact-step 4 3 entries '(1 2 3) 1 r)) (expt 2 -40))))
    ;; A mu that cannot be used leaves no factor.
    (dolist (bad (list -1d0 nan infinity))
      (assert (= 0 (dm-caller:factor-code plan 1d0)))
      (assert (= 3 (dm-caller:factor-code plan bad)) () "A mu of ~S was taken." bad)
      (assert (not (dm:factor-ready? plan)))
      (assert (= 1 (dm-caller:step-code plan r out))))
    ;; Steps that cannot be made leave the output.
    (assert (= 0 (dm-caller:factor-code plan 1d0)))
    (assert (= 2 (dm-caller:step-code plan (doubles 1d0 2d0 3d0) out)))
    (assert (= 2 (dm-caller:step-code plan r (doubles 9d0 9d0))))
    (assert (= 3 (dm-caller:step-code plan (doubles 1d0 nan 3d0 5d0) out)))
    (assert (= 3 (dm-caller:step-code plan (doubles 1d0 2d0 infinity 5d0) out)))
    (assert (same-bits-p out untouched) () "A refused step wrote ~S." out)
    ;; Setting the matrix or the damping withdraws the factor.
    (assert (dm:factor-ready? plan))
    (assert (= 0 (dm-caller:matrix-code plan matrix)))
    (assert (and (not (dm:factor-ready? plan)) (= 1 (dm-caller:step-code plan r out))))
    (assert (= 0 (dm-caller:factor-code plan 1d0)))
    (assert (= 0 (dm-caller:damping-code plan (doubles 3d0 2d0 1d0))))
    (assert (not (dm:factor-ready? plan)))
    ;; A matrix of another pattern or shape is a typed error and leaves none set. A stored
    ;; matrix holds only finite values, so there is no such matrix to refuse.
    (assert (= 91 (dm-caller:matrix-code plan (csc-from-entries 4 3 (cons '(3 0 1) entries)))))
    (assert (not (dm:matrix-set? plan)))
    (assert (= 4 (dm-caller:factor-code plan 1d0)))
    (assert (= 0 (dm-caller:matrix-code plan matrix)))
    (assert (= 92 (dm-caller:matrix-code plan (csc-from-entries 5 3 entries))))
    (assert (= 92 (dm-caller:matrix-code plan (csc-from-entries 4 4 entries))))
    ;; The columns' own damping is formed again from each matrix: twice the matrix, four times
    ;; the damping, and with a residual twice as large the same step.
    (let ((first (zeros 3)) (second (zeros 3)))
      (assert (= 0 (dm-caller:matrix-code plan matrix)))
      (assert (= 0 (dm-caller:column-damping-code plan 0d0)))
      (assert (= 0 (dm-caller:factor-code plan 0.5d0)))
      (assert (= 0 (dm-caller:step-code plan r first)))
      (assert (= 0 (dm-caller:matrix-code plan (damped-scaled-matrix 4 3 entries 1))))
      (assert (= 0 (dm-caller:factor-code plan 0.5d0)))
      (assert (= 0 (dm-caller:step-code plan (doubles 2d0 -4d0 6d0 10d0) second)))
      (assert (same-bits-p first second) () "The columns' damping did not follow the matrix: ~S then ~S." first second))
    ;; A matrix without full column rank: with no damping the factorization refuses a pivot and
    ;; names the column; with any it succeeds.
    (let* ((deficient '((0 0 1) (0 1 2) (1 0 2) (1 1 4) (2 0 3) (2 1 6) (2 2 1)))
           (singular (csc-from-entries 3 3 deficient)) (plan (dm:analyze singular)) (step (zeros 3)))
      (assert (= 0 (dm-caller:matrix-code plan singular)))
      (assert (= 0 (dm-caller:damping-code plan (doubles 1d0 1d0 1d0))))
      (assert (= 1 (dm-caller:factor-code plan 0d0)) () "A matrix of rank 2 was factored undamped: ~D."
              (dm-caller:factor-code plan 0d0))
      (assert (and (not (dm:factor-ready? plan)) (<= 0 (dm:refused-column plan) 1)) ()
              "The refused column is ~D." (dm:refused-column plan))
      (assert (= 0 (dm-caller:factor-code plan (scale-float 1d0 -20))))
      (assert (= -1 (dm:refused-column plan)))
      (assert (= 0 (dm-caller:step-code plan (doubles 1d0 2d0 3d0) step))))
    ;; A column of zeros. Its squared length is zero: a floor of zero is refused; the smallest
    ;; normal number as a floor is taken, the factorization succeeds at a small mu, and that
    ;; column's entry of the step is exactly zero.
    (let* ((holed '((0 0 2) (1 0 1) (0 2 3) (2 2 1) (1 1 0)))
           (matrix (csc-from-entries 3 3 holed)) (plan (dm:analyze matrix)) (step (doubles 9d0 9d0 9d0)))
      (assert (= 0 (dm-caller:matrix-code plan matrix)))
      (assert (= 2 (dm-caller:column-damping-code plan 0d0)))
      (assert (not (dm:damping-set? plan)))
      (assert (= 0 (dm-caller:column-damping-code plan least-positive-normalized-double-float)))
      (assert (= 0 (dm-caller:factor-code plan 1d-3)) () "A column of zeros under a floor refused: ~D."
              (dm-caller:factor-code plan 1d-3))
      (assert (= 0 (dm-caller:step-code plan (doubles 1d0 2d0 3d0) step)))
      (assert (= 0d0 (aref step 1)) () "The step along a column of zeros is ~S." (aref step 1))
      ;; The same matrix 2^600 times larger: the floor, divided with it, would vanish, and is
      ;; held at the least normal number instead, so the pivot is still positive.
      (let ((large (damped-scaled-matrix 3 3 holed 600)))
        (assert (= 0 (dm-caller:matrix-code plan large)))
        (assert (= 0 (dm-caller:factor-code plan 1d-3)))
        (assert (= 0 (dm-caller:step-code plan (doubles 1d0 2d0 3d0) step)))
        (assert (= 0d0 (aref step 1)))))
    ;; A step that cannot be held: a matrix of 1e-200 under a residual of 1e200.
    (let* ((tiny (csc-from-entries 1 1 '((0 0 1d-200)))) (plan (dm:analyze tiny)) (step (doubles 9d0)))
      (assert (= 0 (dm-caller:matrix-code plan tiny)))
      (assert (= 0 (dm-caller:column-damping-code plan 0d0)))
      (assert (= 0 (dm-caller:factor-code plan 1d-6)))
      (assert (= 4 (dm-caller:step-code plan (doubles 1d200) step)))
      (assert (= 9d0 (aref step 0)) () "A step that cannot be held wrote ~S." (aref step 0))
      (assert (= 0 (dm-caller:step-code plan (doubles 1d-200) step)))
      (assert (<= (abs (+ (aref step 0) (/ 1d0 (+ 1d0 1d-6)))) 1d-15)))
    ;; A damping far beyond the matrix: with mu zero it is not used at all, and with mu positive
    ;; the damped diagonal cannot be held.
    (let* ((small (csc-from-entries 1 1 (list (list 0 0 (scale-float 1d0 -600))))) (plan (dm:analyze small))
           (step (zeros 1)))
      (assert (= 0 (dm-caller:matrix-code plan small)))
      (assert (= 0 (dm-caller:damping-code plan (doubles 1d300))))
      (assert (= 0 (dm-caller:factor-code plan 0d0)))
      (assert (= 0 (dm-caller:step-code plan (doubles (scale-float 1d0 -600)) step)))
      (assert (= -1d0 (aref step 0)))
      (assert (= 2 (dm-caller:factor-code plan 1d0))))
    ;; A step multiplied back by two powers whose own product cannot be held, 2^520 for the
    ;; matrix and 2^520 for the residual, though the step can: the residual is all but
    ;; orthogonal to the one column, and the step is -2^988 exactly.
    (let* ((column (csc-from-entries 2 1 (list (list 0 0 (scale-float 1d0 -520)) (list 1 0 (scale-float 1d0 -520)))))
           (plan (dm:analyze column)) (step (zeros 1)))
      (assert (= 0 (dm-caller:matrix-code plan column)))
      (assert (= 0 (dm-caller:column-damping-code plan 0d0)))
      (assert (= 0 (dm-caller:factor-code plan 1d0)))
      (assert (= 0 (dm-caller:step-code plan (doubles (scale-float 1d0 520)
                                                      (scale-float (- (scale-float 1d0 -50) 1d0) 520))
                                        step)))
      (assert (= (- (scale-float 1d0 988)) (aref step 0)) () "The step multiplied back in two stages is ~S." (aref step 0)))
    ;; A damping divided by a power whose square cannot be held, 2^1076 for a matrix of entries
    ;; near 2^-540, though the divided damping can: the given entries are the smallest numbers
    ;; the format has and count for as much as the matrix. Against the exact step.
    (let* ((small-entries (mapcar (lambda (entry)
                                    (list (first entry) (second entry) (scale-float (coerce (third entry) 'double-float) -540)))
                                  entries))
           (small (csc-from-entries 4 3 small-entries)) (plan (dm:analyze small)) (step (zeros 3))
           (d (doubles (scale-float 1d0 -1073) (scale-float 1d0 -1074) (scale-float 1d0 -1072)))
           (residual (map '(simple-array double-float (*)) (lambda (value) (scale-float value -540)) r)))
      (assert (= 0 (dm-caller:matrix-code plan small)))
      (assert (= 0 (dm-caller:damping-code plan d)))
      (assert (= 0 (dm-caller:factor-code plan 1d0)) () "A damping of the smallest numbers refused: ~D."
              (dm-caller:factor-code plan 1d0))
      (assert (= 0 (dm-caller:step-code plan residual step)))
      (let ((error (damped-error step (damped-exact-step 4 3 small-entries (map 'list #'rational d) 1 residual))))
        (assert (<= error (expt 2 -40)) () "With a damping of the smallest numbers the step is off by ~E."
                (float error 1d0))))
    ;; A square matrix: the output may be the residual's own array.
    (let* ((square '((0 0 4) (1 0 1) (1 1 3) (2 1 1) (2 2 2) (0 2 1)))
           (matrix (csc-from-entries 3 3 square)) (plan (dm:analyze matrix))
           (apart (zeros 3)) (shared (doubles 1d0 -2d0 3d0)))
      (assert (= 0 (dm-caller:matrix-code plan matrix)))
      (assert (= 0 (dm-caller:column-damping-code plan 0d0)))
      (assert (= 0 (dm-caller:factor-code plan 0.25d0)))
      (assert (= 0 (dm-caller:step-code plan (doubles 1d0 -2d0 3d0) apart)))
      (assert (= 0 (dm-caller:step-code plan shared shared)))
      (assert (same-bits-p apart shared)))
    ;; A budget too small for the product or the factor is the analyses' typed error.
    (assert (= 98 (dm-caller:analysis-code matrix 2)))
    (assert (= 0 (dm-caller:analysis-code matrix 100)))
    ;; Steady state: a compiled caller allocates nothing, over several mu on one analysis.
    (let ((plan (dm:analyze matrix)) (step (zeros 3)) (dampings (doubles 1d-3 1d-2 1d-1 1d0)))
      (assert (= 0 (dm-caller:matrix-code plan matrix)))
      (assert (= 0 (dm-caller:column-damping-code plan least-positive-normalized-double-float)))
      (allocation-check "damped-loop" (lambda () (assert (= 0 (dm-caller:damped-loop plan matrix dampings r step 1))))))
    (format t "DAMPED-OUTCOMES-PASS~%")))

(defparameter *damped-public-control*
  "(declare probe (coalton-linalg/sparse/damped:DampedPlan -> UFix))
   (define (probe plan)
     (+ (coalton-linalg/sparse/damped:row-count plan) (coalton-linalg/sparse/damped:column-count plan)))")

(defparameter *damped-encapsulation-probes*
  '(("damped-plan-entries-access" .
     "(declare probe (coalton-linalg/sparse/damped:DampedPlan -> arr:LispArray F64))
      (define (probe plan) (.entries plan))")
    ("damped-plan-symbolic-access" .
     "(declare probe (coalton-linalg/sparse/damped:DampedPlan -> ldlt:Symbolic))
      (define (probe plan) (.symbolic plan))")
    ("damped-plan-unwrapping" .
     "(declare probe (coalton-linalg/sparse/damped:DampedPlan -> UFix))
      (define (probe plan)
        (match plan ((coalton-linalg/sparse/damped:MkDampedPlan storage) (.rows storage))))")
    ("damping-status-as-damped-status" .
     "(declare probe (coalton-linalg/sparse/damped:DampingStatus -> coalton-linalg/sparse/damped:DampedStatus))
      (define (probe status) status)")
    ("damped-status-as-step-status" .
     "(declare probe (coalton-linalg/sparse/damped:DampedStatus -> coalton-linalg/sparse/damped:StepStatus))
      (define (probe status) status)")))

(defun run-damped-encapsulation-tests ()
  "External Coalton callers reach the plan only through its public operations, and the status types stay apart."
  (multiple-value-bind (accepted diagnostics)
      (coalton-linalg/tests/encapsulation::compile-probe "damped-public-control" *damped-public-control*)
    (assert accepted () "The valid damped external caller failed:~%~A" diagnostics))
  (dolist (probe *damped-encapsulation-probes*)
    (assert (not (coalton-linalg/tests/encapsulation::compile-probe (car probe) (cdr probe))) ()
            "External caller ~A compiled despite required rejection." (car probe)))
  (dolist (name '("MKDAMPEDPLAN" "DAMPEDSTORAGE"))
    (multiple-value-bind (symbol visibility) (find-symbol name "COALTON-LINALG/SPARSE/DAMPED")
      (assert (and symbol (eq visibility ':internal)) ()
              "Storage constructor ~A must exist and remain internal; visibility is ~S." name visibility)))
  (format t "DAMPED-ENCAPSULATION-PASS rejected=~D public-control=1~%" (length *damped-encapsulation-probes*)))

(defun run-damped-tests ()
  (run-damped-exact-tests)
  (run-damped-scale-tests)
  (run-damped-outcome-tests)
  (run-damped-encapsulation-tests)
  (format t "DAMPED-TESTS-PASS~%"))
