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
    ;; The same of a floor: 1e-12 under a matrix at 2^-700 is not a number the format holds in the
    ;; divided matrix's units. The floor is taken, since it stays for later matrices; the factor is
    ;; where it is refused, and only with mu positive. Under the next matrix, at ordinary scale, the
    ;; same floor stands and is used.
    (let* ((small (csc-from-entries 1 1 (list (list 0 0 (scale-float 1d0 -700)))))
           (ordinary (csc-from-entries 1 1 (list (list 0 0 1d0))))
           (plan (dm:analyze small)) (step (doubles 9d0)))
      (assert (= 0 (dm-caller:matrix-code plan small)))
      (assert (= 0 (dm-caller:column-damping-code plan 1d-12)) () "A floor the scale takes out of range was not taken.")
      (let ((code (dm-caller:factor-code plan 1d0)))
        (assert (= 2 code) () "A floor of 1e-12 under a matrix at 2^-700: the factor's outcome is ~D." code))
      (assert (not (dm:factor-ready? plan)))
      (let ((code (dm-caller:step-code plan (doubles (scale-float 1d0 -700)) step)))
        (assert (and (= 1 code) (= 9d0 (aref step 0))) () "After that refusal the step's outcome is ~D and it wrote ~S."
                code (aref step 0)))
      (assert (= 0 (dm-caller:factor-code plan 0d0)))
      (assert (= 0 (dm-caller:step-code plan (doubles (scale-float 1d0 -700)) step)))
      (assert (= -1d0 (aref step 0)) () "With mu zero the step is ~S." (aref step 0))
      (assert (= 0 (dm-caller:matrix-code plan ordinary)))
      (let ((code (dm-caller:factor-code plan 1d0)))
        (assert (= 0 code) () "The same floor under a matrix at ordinary scale: the factor's outcome is ~D." code))
      (assert (= 0 (dm-caller:step-code plan (doubles 1d0) step)))
      (assert (<= (abs (+ (aref step 0) (/ 1d0 (+ 2d0 1d-12)))) 1d-15) () "The step with the floor is ~S." (aref step 0)))
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

(defun damped-diagonal-case (label diagonal kind damping mu residual)
  "A diagonal matrix, where the step is known exactly: entry j is -J_jj r_j / (J_jj^2 + mu D_j), D_j the entry
given or the column's square plus the floor. Asserts each outcome and the step to a relative 2^-48, an entry
whose exact value is zero being zero. Returns the plan."
  (let* ((n (length diagonal))
         (matrix (csc-from-entries n n (loop :for j :below n :for value :in diagonal :collect (list j j value))))
         (plan (dm:analyze matrix)) (out (zeros n)))
    (assert (= 0 (dm-caller:matrix-code plan matrix)) () "~A: the matrix was refused." label)
    (let ((code (ecase kind
                  (:given (dm-caller:damping-code plan (coerce damping '(simple-array double-float (*)))))
                  (:columns (dm-caller:column-damping-code plan damping)))))
      (assert (= 0 code) () "~A: the damping's outcome is ~D." label code))
    (let ((code (dm-caller:factor-code plan mu)))
      (assert (= 0 code) () "~A: the factor's outcome is ~D, refused column ~D." label code (dm:refused-column plan)))
    (let ((code (dm-caller:step-code plan (coerce residual '(simple-array double-float (*))) out)))
      (assert (= 0 code) () "~A: the step's outcome is ~D." label code))
    (loop :for j :from 0 :for value :in diagonal :for r :in residual
          :do (let* ((a (rational value))
                     (d (ecase kind (:given (rational (elt damping j))) (:columns (+ (* a a) (rational damping)))))
                     (denominator (+ (* a a) (* (rational mu) d)))
                     (exact (if (zerop denominator) 0 (- (/ (* a (rational r)) denominator))))
                     (found (rational (aref out j))))
                (assert (if (zerop exact) (zerop found) (<= (abs (- found exact)) (* (abs exact) (expt 2 -48)))) ()
                        "~A: entry ~D of the step is ~S, and exactly it is ~S." label j (aref out j) (float exact 1d0))))
    plan))

(defun run-damped-range-tests ()
  "Sizes inside one problem that one power of two for the whole matrix could not hold: a mu that is large
because the damping given is small beside the matrix squared, and columns far apart. Every step against its
exact value."
  (flet ((two (power) (scale-float 1d0 power)))
    ;; A damping entry far below the matrix squared, under a mu that makes it count.
    (damped-diagonal-case "J = 2^520, D = 1, mu = 2^1023" (list (two 520)) :given (list 1d0) (two 1023) (list (two 520)))
    (damped-diagonal-case "J = 2^512, D = 1, mu = 2^1000" (list (two 512)) :given (list 1d0) (two 1000) (list (two 512)))
    (damped-diagonal-case "J = diag(1, 2^-450), D = (1, 2^-1030), mu = 2^100" (list 1d0 (two -450)) :given
                          (list 1d0 (two -1030)) (two 100) (list 1d0 (two -450)))
    (damped-diagonal-case "J = 2^-520, D = 2^-1060, mu = 2^20" (list (two -520)) :given (list (two -1060)) (two 20)
                          (list 3d0))
    ;; The same with numbers that are not powers of two, so that a damping entry carried below the
    ;; normal range on its way would lose digits and not only be moved.
    (damped-diagonal-case "J = 1.5 2^520, D = 1.2345678901234567, mu = 0.75 2^1023" (list (* 1.5d0 (two 520))) :given
                          (list 1.2345678901234567d0) (* 0.75d0 (two 1023)) (list (* 1.1d0 (two 520))))
    (damped-diagonal-case "J = diag(1.3, 1.7 2^-450), D = (0.9, 1.2345678901234567 2^-1000), mu = 1.1 2^100"
                          (list 1.3d0 (* 1.7d0 (two -450))) :given (list 0.9d0 (* 1.2345678901234567d0 (two -1000)))
                          (* 1.1d0 (two 100)) (list 1d0 (* 1.9d0 (two -450))))
    ;; A damping far above the column's square under a mu small enough to make the two comparable.
    (damped-diagonal-case "J = 1.5 2^-700, D = 1.3, mu = 1.1 2^-1040" (list (* 1.5d0 (two -700))) :given (list 1.3d0)
                          (* 1.1d0 (two -1040)) (list 1d0))
    ;; The same of a floor under the columns' own damping.
    (damped-diagonal-case "J = diag(1, 2^-512), a floor of 2^-1030, mu = 1" (list 1d0 (two -512)) :columns (two -1030) 1d0
                          (list 1d0 (two -512)))
    (damped-diagonal-case "J = diag(1, 2^-512), a floor of 2^-1030, mu = 2^40" (list 1d0 (two -512)) :columns (two -1030)
                          (two 40) (list 1d0 (two -512)))
    ;; Columns far apart: neither is a column of zeros, and the matrix has full column rank.
    (damped-diagonal-case "J = diag(1, 2^-540), the columns' own damping, mu = 0" (list 1d0 (two -540)) :columns 0d0 0d0
                          (list 1d0 (two -540)))
    (damped-diagonal-case "J = diag(1, 2^-540), the columns' own damping, mu = 1" (list 1d0 (two -540)) :columns 0d0 1d0
                          (list 1d0 (two -540)))
    (damped-diagonal-case "J = diag(2^400, 2^-500), the columns' own damping, mu = 2^-3" (list (two 400) (two -500))
                          :columns 0d0 0.125d0 (list (two 100) (two -100)))
    (damped-diagonal-case "J = diag(2^500, 1, 2^-500), D = (2^1000, 1, 2^-1000), mu = 4"
                          (list (two 500) 1d0 (two -500)) :given (list (two 1000) 1d0 (two -1000)) 4d0
                          (list 1d0 1d0 1d0))
    (damped-diagonal-case "J = diag(1, 2^-540), D = (1, 2^-1060), mu = 0" (list 1d0 (two -540)) :given
                          (list 1d0 (two -1060)) 0d0 (list 1d0 (two -540))))
  ;; A step entry that is an ordinary number, and far below 1 in the units of the solve: the
  ;; residual's small entry under a large damping. One power of two for the residual put the
  ;; solve's unknown there below every number, and the entry came back as zero under StepSuccess.
  (flet ((two (power) (scale-float 1d0 power)))
    (damped-diagonal-case "the identity, D = (1, 2^700), mu = 1, r = (2^200, 2^-200)" (list 1d0 1d0) :given
                          (list 1d0 (two 700)) 1d0 (list (two 200) (two -200)))
    (damped-diagonal-case "diag(3, 5 2^40), D = (2^-300, 1.3 2^650), mu = 2^20, r = (1.1 2^300, 1.7 2^-300)"
                          (list 3d0 (* 5d0 (two 40))) :given (list (two -300) (* 1.3d0 (two 650))) (two 20)
                          (list (* 1.1d0 (two 300)) (* 1.7d0 (two -300))))
    (damped-diagonal-case "diag(1.5, 1.25, 1.75), D = (2^900, 1, 2^-900), mu = 1, r = (2^-100, 1, 2^-100)"
                          (list 1.5d0 1.25d0 1.75d0) :given (list (two 900) 1d0 (two -900)) 1d0
                          (list (two -100) 1d0 (two -100)))
    (damped-diagonal-case "diag(1, 2^-300), the columns' own damping with a floor of 2^-100, mu = 2^300, r = (2^100, 2^-200)"
                          (list 1d0 (two -300)) :columns (two -100) (two 300) (list (two 100) (two -200)))
    ;; The same where the matrix is not diagonal: the second column meets only the rows whose
    ;; residual is small. Against the step in rational arithmetic, entry by entry.
    (let* ((entries '((0 0 1) (1 0 1) (2 1 1) (3 1 2))) (matrix (csc-from-entries 4 2 entries)) (plan (dm:analyze matrix))
           (r (doubles (two 200) (* 3d0 (two 200)) (two -200) (* 5d0 (two -200)))) (out (zeros 2))
           (exact (damped-exact-step 4 2 entries (list 1 (expt 2 700)) 1 r)))
      (assert (= 0 (dm-caller:matrix-code plan matrix)))
      (assert (= 0 (dm-caller:damping-code plan (doubles 1d0 (two 700)))))
      (assert (= 0 (dm-caller:factor-code plan 1d0)))
      (assert (= 0 (dm-caller:step-code plan r out)))
      (dotimes (j 2)
        (assert (<= (abs (- (rational (aref out j)) (aref exact j))) (* (abs (aref exact j)) (expt 2 -48))) ()
                "Two columns on separate rows: entry ~D of the step is ~S, and exactly it is ~S." j (aref out j)
                (float (aref exact j) 1d0)))))
  ;; What one power of two for the residual, and one for the solve's right-hand side, cannot
  ;; carry is refused by name, with the output as it was, and not returned as a step with an entry
  ;; missing: a residual whose entries are more than 2^1022 apart, and a gradient that is so in the
  ;; units of the solve.
  (flet ((two (power) (scale-float 1d0 power))
         (attempt (diagonal damping mu residual)
           (let* ((n (length diagonal))
                  (matrix (csc-from-entries n n (loop :for j :below n :for value :in diagonal :collect (list j j value))))
                  (plan (dm:analyze matrix)) (out (make-array n :element-type 'double-float :initial-element 9d0)))
             (assert (= 0 (dm-caller:matrix-code plan matrix)))
             (assert (= 0 (dm-caller:damping-code plan (coerce damping '(simple-array double-float (*))))))
             (assert (= 0 (dm-caller:factor-code plan mu)))
             (values (dm-caller:step-code plan (coerce residual '(simple-array double-float (*))) out) out))))
    (dolist (case (list (list "a residual spanning 2^1100" (list 1d0 1d0) (list 1d0 1d0) 0d0 (list (two 600) (two -500)))
                        (list "a residual spanning 2^1030 beside a zero" (list 1d0 1d0 1d0) (list 1d0 1d0 1d0) 1d0
                              (list (two 500) 0d0 (* 1.2345678901234567d0 (two -530))))
                        (list "a gradient spanning 2^1100 in the solve's units" (list 1d0 1d0) (list 1d0 (two 1000)) 1d0
                              (list 1d0 (two -600)))
                        ;; Every division here is exact, the small entry being on the grid below the
                        ;; normal range; the solve's entry for it is below that range too, and the
                        ;; step's own entry, 3 2^-532, is not.
                        (list "an entry of the solve below the normal range" (list 1d0 1d0 1d0) (list 1d0 1d0 1d0) 1d0
                              (list (two 500) 0d0 (* 3d0 (two -531))))))
      (destructuring-bind (label diagonal damping mu residual) case
        (multiple-value-bind (code out) (attempt diagonal damping mu residual)
          (assert (and (= 5 code) (every (lambda (value) (= value 9d0)) out)) ()
                  "~A: the step's outcome is ~D and the output ~S." label code out))))
    ;; Numbers below the normal range that were given so, and that no division of this package's
    ;; made so, are not refused: the step is exact.
    (damped-diagonal-case "the identity, D = (1, 1), mu = 1, r = (1, 3 2^-1070)" (list 1d0 1d0) :given (list 1d0 1d0) 1d0
                          (list 1d0 (* 3d0 (two -1070))))
    ;; Just inside what is carried: a residual spanning 2^1000 under no damping to speak of.
    (damped-diagonal-case "the identity, D = (1, 1), mu = 2^-10, r = (2^500, 1.25 2^-500)" (list 1d0 1d0) :given
                          (list 1d0 1d0) (two -10) (list (two 500) (* 1.25d0 (two -500)))))
  ;; refused-column names a pivot that was refused, and is -1 after every other refusal.
  (let* ((small (csc-from-entries 1 1 (list (list 0 0 (scale-float 1d0 -600))))) (plan (dm:analyze small)))
    (assert (= 0 (dm-caller:matrix-code plan small)))
    (assert (= 0 (dm-caller:damping-code plan (doubles 1d300))))
    (assert (= 2 (dm-caller:factor-code plan 1d0)))
    (assert (= -1 (dm:refused-column plan)) () "After a refusal that was not of a pivot, refused-column is ~D."
            (dm:refused-column plan)))
  (let* ((matrix (csc-from-entries 2 2 (list (list 0 0 1d0) (list 1 1 0d0)))) (plan (dm:analyze matrix)))
    (assert (= 0 (dm-caller:matrix-code plan matrix)))
    (assert (= 0 (dm-caller:damping-code plan (doubles 1d0 1d0))))
    (assert (and (= 1 (dm-caller:factor-code plan 0d0)) (= 1 (dm:refused-column plan))))
    (assert (and (= 3 (dm-caller:factor-code plan -1d0)) (= -1 (dm:refused-column plan))) ()
            "After a mu that cannot be used, refused-column is ~D." (dm:refused-column plan))
    (assert (and (= 1 (dm-caller:factor-code plan 0d0)) (= 1 (dm:refused-column plan))))
    (assert (and (= 0 (dm-caller:factor-code plan 1d0)) (= -1 (dm:refused-column plan)))))
  ;; A column of zeros. Under the columns' own damping with a floor of zero its damping entry is zero:
  ;; the floor is refused where it is given, and where the column of zeros comes with a later matrix
  ;; the factor says so by name, at every mu, and names the column. With mu zero it is the pivot.
  (let* ((identity (csc-from-entries 2 2 (list (list 0 0 1d0) (list 1 1 1d0))))
         (zero-column (csc-from-entries 2 2 (list (list 0 0 1d0) (list 1 1 0d0))))
         (plan (dm:analyze identity)) (out (doubles 9d0 9d0)))
    (assert (= 0 (dm-caller:matrix-code plan zero-column)))
    (assert (= 2 (dm-caller:column-damping-code plan 0d0)))
    (assert (= 0 (dm-caller:matrix-code plan identity)))
    (assert (= 0 (dm-caller:column-damping-code plan 0d0)))
    (assert (= 0 (dm-caller:matrix-code plan zero-column)))
    (assert (dm:damping-set? plan))
    (dolist (mu (list 1d0 1d100 1d300 (scale-float 1d0 -1000)))
      (let ((code (dm-caller:factor-code plan mu)))
        (assert (and (= 5 code) (= 1 (dm:refused-column plan)) (not (dm:factor-ready? plan))) ()
                "A column of zeros under its own damping, mu = ~S: outcome ~D, refused column ~D." mu code
                (dm:refused-column plan)))
      (assert (and (= 1 (dm-caller:step-code plan (doubles 1d0 1d0) out)) (= 9d0 (aref out 1)))))
    (assert (and (= 1 (dm-caller:factor-code plan 0d0)) (= 1 (dm:refused-column plan))))
    ;; The next matrix has no such column, and the same damping serves.
    (assert (= 0 (dm-caller:matrix-code plan identity)))
    (assert (and (= 0 (dm-caller:factor-code plan 1d0)) (= -1 (dm:refused-column plan)))))
  ;; With a damping entry that is positive, any positive pivot gives the column of zeros a step
  ;; entry of exactly zero, and its pivot is taken to be 1: neither a floor nor a mu near the top
  ;; or the bottom of the range can be refused for it, and the other entry is exact.
  (dolist (case (list (list :columns 1d-12 1d0) (list :columns 1d300 1d0) (list :columns 1d-300 1d-300)
                      (list :given 1d300 1d300) (list :given 1d-300 1d-300) (list :given 1d0 8d0)))
    (destructuring-bind (kind damping mu) case
      (let* ((matrix (csc-from-entries 2 2 (list (list 0 0 3d0) (list 1 1 0d0)))) (plan (dm:analyze matrix))
             (out (doubles 9d0 9d0)))
        (assert (= 0 (dm-caller:matrix-code plan matrix)))
        (assert (= 0 (ecase kind
                       (:columns (dm-caller:column-damping-code plan damping))
                       (:given (dm-caller:damping-code plan (doubles 2d0 damping))))))
        (let ((code (dm-caller:factor-code plan mu)))
          (assert (= 0 code) () "A column of zeros, ~A damping ~S, mu ~S: outcome ~D." kind damping mu code))
        (assert (= 0 (dm-caller:step-code plan (doubles 1d0 5d0) out)))
        (let* ((first-damping (ecase kind (:columns (+ 9 (rational damping))) (:given 2)))
               (exact (- (/ 3 (+ 9 (* (rational mu) first-damping))))))
          (assert (and (zerop (aref out 1)) (<= (abs (- (rational (aref out 0)) exact)) (* (abs exact) (expt 2 -48)))) ()
                  "A column of zeros, ~A damping ~S, mu ~S: the step is ~S, and its first entry is exactly ~S." kind damping
                  mu out (float exact 1d0))))))
  (format t "DAMPED-RANGE-PASS~%"))

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
  (run-damped-range-tests)
  (run-damped-encapsulation-tests)
  (format t "DAMPED-TESTS-PASS~%"))
