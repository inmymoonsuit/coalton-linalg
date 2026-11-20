(in-package #:coalton-linalg/tests)

;;; Iterative refinement: the iteration is recomputed independently with the same building
;;; blocks and must agree bit for bit; the result is compared with an exact rational solution.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "RF" "COALTON-LINALG/REFINEMENT")
  (sb-ext:add-package-local-nickname "RF-CALLER" "COALTON-LINALG/TESTS/REFINEMENT-CALLERS"))

(defun full-csc (dense)
  "Every entry of a dense square matrix of numbers as a checked CSC matrix of doubles."
  (let ((n (array-dimension dense 0)))
    (csc-from-entries n n (loop :for j :below n :append
                            (loop :for i :below n :collect (list i j (aref dense i j)))))))

(defun dense-lower-pattern (dense)
  "Lower CSC arrays of a dense symmetric matrix, every entry of the lower triangle stored."
  (let ((n (array-dimension dense 0)) (pointers (list 0)) (rows '()) (values '()))
    (dotimes (j n)
      (loop :for i :from j :below n :do
        (push i rows)
        (push (coerce (aref dense i j) 'double-float) values))
      (push (length rows) pointers))
    (values n (apply #'indices (reverse pointers)) (apply #'indices (reverse rows)) (apply #'doubles (reverse values)))))

(defun infinity-norm (vector)
  (loop :for v :across vector :maximize (abs v)))

(defun exact-solution (dense b)
  "A^-1 b in rational arithmetic, from the double values of A and b."
  (let* ((n (array-dimension dense 0))
         (inverse (rational-inverse (build-matrix n (lambda (i j) (rational (coerce (aref dense i j) 'double-float)))))))
    (coerce (loop :for i :below n :collect (loop :for j :below n :sum (* (aref inverse i j) (rational (aref b j)))))
            'vector)))

(defun reference-refinement (matrix solve b start steps)
  "The iterates and residual norms of STEPS corrections, recomputed here with the library's own product:
x <- x + solve(b - A x). Returns the list of iterates and the list of norms, oldest first."
  (let* ((n (length b)) (x (copy-seq start)) (product (zeros n)) (residual (zeros n)) (correction (zeros n))
         (iterates '()) (norms '()))
    (dotimes (k (1+ steps))
      (csc:matvec! matrix x product)
      (dotimes (i n) (setf (aref residual i) (- (aref b i) (aref product i))))
      (push (copy-seq x) iterates)
      (push (infinity-norm residual) norms)
      (when (< k steps)
        (funcall solve residual correction)
        (dotimes (i n) (setf (aref x i) (+ (aref x i) (aref correction i))))))
    (values (nreverse iterates) (nreverse norms))))

(defun test-system (state n shift)
  "A symmetric positive definite matrix of integers, its operator, a solving factor of the matrix plus SHIFT
times the identity as the approximate solve, and a right-hand side."
  (let* ((basis (make-array (list (+ n 2) n))) (dense (make-array (list n n))))
    (dotimes (i (+ n 2)) (dotimes (j n) (setf (aref basis i j) (- (random 7 state) 3))))
    (dotimes (j n) (setf (aref basis j j) (+ 4 (random 3 state))))
    (dotimes (i n)
      (dotimes (j n)
        ;; Twenty on the diagonal keeps the smallest eigenvalue away from the shifts used below.
        (setf (aref dense i j) (+ (if (= i j) 20 0)
                                  (loop :for k :below (+ n 2) :sum (* (aref basis k i) (aref basis k j)))))))
    (multiple-value-bind (order pointers rows values) (dense-lower-pattern dense)
      (let* ((symbolic (ldlt:analyze (ldlt:pattern order pointers rows)))
             (matrix (full-csc dense))
             (b (make-array n :element-type 'double-float)))
        (dotimes (i n) (setf (aref b i) (coerce (- (random 19 state) 9) 'double-float)))
        ;; factor! subtracts its shift, so a negative one adds to the diagonal.
        (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess)
                    (ldlt:factor! symbolic values (coerce (- shift) 'double-float))))
        (values dense matrix symbolic b)))))

(defun run-refinement-iteration-tests ()
  "The history is the iteration's, the solution is the best iterate, and the best iterate solves the system."
  (let ((state (sb-ext:seed-random-state 20261202)) (cases 0) (most-steps 0))
    (dolist (n '(1 2 3 5 8))
      (dolist (shift '(0 1/64 1/4 2))
        (dotimes (repeat 3)
          (multiple-value-bind (dense matrix symbolic b) (test-system state n shift)
            (let* ((factor (ldlt:symbolic-factor symbolic))
                   (apply (rf-caller:matrix-operator matrix)) (solve (rf-caller:factor-solver factor))
                   (refiner (rf:make-refiner n 60)) (x (zeros n))
                   (tolerance (* 64 double-float-epsilon (max 1d0 (infinity-norm b))))
                   (code (rf-caller:refine-code refiner apply solve b x tolerance 1d0))
                   (steps (rf:step-count refiner)))
              (incf cases)
              (setf most-steps (max most-steps steps))
              (assert (= 0 code) () "n=~D shift ~A: refinement ended with code ~D after ~D steps, norms ~S." n shift
                      code steps (subseq (rf:residual-norms refiner) 0 (rf:residual-count refiner)))
              (assert (= (rf:residual-count refiner) (1+ steps)))
              (assert (= (rf:best-step refiner) steps) () "The converged iterate is not the returned one.")
              (assert (= 60 (rf:step-limit refiner)))
              (assert (= n (rf:dimension refiner)))
              ;; The same iteration recomputed here: every norm and the returned vector, bit for bit.
              (multiple-value-bind (iterates norms)
                  (reference-refinement matrix (lambda (r out) (ldlt:solve! factor r out)) b (zeros n) steps)
                (assert (same-bits-p (subseq (rf:residual-norms refiner) 0 (1+ steps))
                                     (coerce norms '(simple-array double-float (*))))
                        () "n=~D shift ~A: recorded norms ~S, recomputed ~S." n shift
                        (subseq (rf:residual-norms refiner) 0 (1+ steps)) norms)
                (assert (same-bits-p x (car (last iterates))) () "The returned vector is not the last iterate.")
                ;; Each step before the last brought the norm down, or the stall test would have stopped it.
                (loop :for (earlier later) :on norms :while later
                      :do (assert (< later earlier) () "Norms ~S did not fall." norms)))
              (assert (<= (aref (rf:residual-norms refiner) steps) tolerance))
              ;; The returned vector against the exact solution: the exact residual is the recorded
              ;; one to rounding, and the error is at most the inverse's norm times it.
              (let* ((exact (exact-solution dense b))
                     (scale (loop :for v :across exact :maximize (abs v)))
                     (error (loop :for i :below n :maximize (abs (- (rational (aref x i)) (aref exact i))))))
                (assert (<= error (* (expt 2 -30) (max scale 1))) ()
                        "n=~D shift ~A: the refined solution is off by ~E." n shift (float error 1d0))))))))
    (format t "REFINEMENT-ITERATION-PASS cases=~D most-steps=~D~%" cases most-steps)))

(defun run-refinement-outcome-tests ()
  (let* ((state (sb-ext:seed-random-state 20261203)) (n 5) (calls (indices 0)))
    (multiple-value-bind (dense matrix symbolic b) (test-system state n 1/4)
      (declare (ignore dense))
      (let* ((factor (ldlt:symbolic-factor symbolic)) (apply (rf-caller:matrix-operator matrix))
             (solve (rf-caller:factor-solver factor)) (plain (lambda (r out) (ldlt:solve! factor r out)))
             (tight (* 64 double-float-epsilon (infinity-norm b))))
        (flet ((scripted (&rest multipliers)
                 (setf (aref calls 0) 0)
                 (rf-caller:scripted-solver solve (apply #'doubles multipliers) calls)))
          ;; A correction of zero: the residual does not fall, the call stalls, and the starting
          ;; vector comes back untouched as the best iterate.
          (let ((refiner (rf:make-refiner n 10)) (x (doubles 1d0 2d0 3d0 4d0 5d0)))
            (assert (= 1 (rf-caller:refine-code refiner apply (scripted 0d0) b x tight 1d0)))
            (assert (equalp x (doubles 1d0 2d0 3d0 4d0 5d0)))
            (assert (and (= 1 (rf:step-count refiner)) (= 0 (rf:best-step refiner)) (= 2 (rf:residual-count refiner))))
            (assert (= (aref (rf:residual-norms refiner) 0) (aref (rf:residual-norms refiner) 1))))
          ;; A good step and then a bad one, with the stall test off: the limit is reached and the
          ;; vector returned is the first iterate, not the last.
          (let ((refiner (rf:make-refiner n 2)) (x (zeros n)))
            (assert (= 2 (rf-caller:refine-code refiner apply (scripted 1d0 -3d0) b x 0d0 1d300)))
            (assert (and (= 2 (rf:step-count refiner)) (= 1 (rf:best-step refiner)) (= 3 (rf:residual-count refiner))))
            (multiple-value-bind (iterates norms) (reference-refinement matrix plain b (zeros n) 1)
              (assert (same-bits-p x (second iterates)) () "The best iterate was not returned: ~S." x)
              (assert (= (second norms) (aref (rf:residual-norms refiner) 1)))
              (assert (> (aref (rf:residual-norms refiner) 2) (aref (rf:residual-norms refiner) 1)))))
          ;; The same bad second step with the stall test on stops there, and reports a stall.
          (let ((refiner (rf:make-refiner n 9)) (x (zeros n)))
            (assert (= 1 (rf-caller:refine-code refiner apply (scripted 1d0 -3d0) b x 0d0 1d0)))
            (assert (and (= 2 (rf:step-count refiner)) (= 1 (rf:best-step refiner))))
            (assert (same-bits-p x (second (reference-refinement matrix plain b (zeros n) 1))) ()
                    "A stalled call returned its last iterate, not its best: ~S." x))
          ;; The stall factor is the caller's: this solve brings the norm down by less than half
          ;; a step, which a factor of 1 accepts and a factor of 1/2 does not.
          (let ((slow (scripted 0.25d0)) (refiner (rf:make-refiner n 4)) (x (zeros n)))
            (assert (= 2 (rf-caller:refine-code refiner apply slow b x 0d0 1d0)))
            (assert (= 4 (rf:step-count refiner)))
            (let ((norms (rf:residual-norms refiner)))
              (assert (> (aref norms 1) (* 0.5d0 (aref norms 0))) () "The slow solve was not slow: ~S." norms))
            (setf (aref calls 0) 0)
            (fill x 0d0)
            (assert (= 1 (rf-caller:refine-code refiner apply slow b x 0d0 0.5d0)))
            (assert (and (= 1 (rf:step-count refiner)) (= 1 (rf:best-step refiner))) ()
                    "A stalled step that still improved was not kept."))
          ;; No steps allowed: the starting vector is judged as it stands.
          (let ((refiner (rf:make-refiner n 0)) (x (zeros n)))
            (assert (= 2 (rf-caller:refine-code refiner apply solve b x tight 1d0)))
            (assert (and (= 0 (rf:step-count refiner)) (= 1 (rf:residual-count refiner))))
            (assert (= 0 (rf-caller:refine-code refiner apply solve b x (* 2 (infinity-norm b)) 1d0))))
          ;; A solve that refuses on its second call: the first step's iterate stands.
          (let ((refiner (rf:make-refiner n 9)) (x (zeros n)))
            (setf (aref calls 0) 0)
            (assert (= 4 (rf-caller:refine-code refiner apply (rf-caller:refusing-solver solve 1 calls) b x 0d0 1d300)))
            (assert (and (= 1 (rf:step-count refiner)) (= 1 (rf:best-step refiner))))
            (assert (same-bits-p x (second (reference-refinement matrix plain b (zeros n) 1)))))
          ;; A right-hand side that is not finite, and a correction that is not: both are outcomes,
          ;; and the starting vector comes back.
          (let ((refiner (rf:make-refiner n 9)) (x (doubles 1d0 2d0 3d0 4d0 5d0))
                (bad (copy-seq b)))
            (setf (aref bad 2) (not-a-number))
            (assert (= 3 (rf-caller:refine-code refiner apply solve bad x tight 1d0)))
            (assert (equalp x (doubles 1d0 2d0 3d0 4d0 5d0)))
            (assert (= 3 (rf-caller:refine-code refiner apply (scripted sb-ext:double-float-positive-infinity)
                                                b x tight 1d0)))
            (assert (equalp x (doubles 1d0 2d0 3d0 4d0 5d0)) () "A non-finite step left ~S." x))
          ;; Arguments.
          (let ((refiner (rf:make-refiner n 3)) (x (doubles 9d0 9d0 9d0 9d0 9d0)))
            (assert (= 5 (rf-caller:refine-code refiner apply solve (doubles 1d0) x tight 1d0)))
            (assert (= 5 (rf-caller:refine-code refiner apply solve b (doubles 1d0) tight 1d0)))
            (assert (= 6 (rf-caller:refine-code refiner apply solve b x -1d0 1d0)))
            (assert (= 6 (rf-caller:refine-code refiner apply solve b x (not-a-number) 1d0)))
            (assert (= 6 (rf-caller:refine-code refiner apply solve b x tight 0d0)))
            (assert (= 6 (rf-caller:refine-code refiner apply solve b x tight (not-a-number))))
            (assert (equalp x (doubles 9d0 9d0 9d0 9d0 9d0)))
            (assert (and (= 0 (rf:step-count refiner)) (= 0 (rf:residual-count refiner))))))))
    ;; An operator that never reads one entry cannot hide a non-finite value there: the matrix
    ;; [[1, 0], [0, 0]] has no entry in its second column, and the correction's second entry is infinite.
    (let* ((matrix (csc-from-entries 2 2 '((0 0 1))))
           (refiner (rf:make-refiner 2 3)) (x (zeros 2)))
      (assert (= 3 (rf-caller:refine-code refiner (rf-caller:matrix-operator matrix)
                                          (rf-caller:poison-solver 1 sb-ext:double-float-positive-infinity)
                                          (doubles 1d0 0d0) x 0d0 1d300)))
      (assert (equalp x (doubles 0d0 0d0)) () "A non-finite entry the operator ignores was returned: ~S." x))
    ;; No unknowns.
    (let ((refiner (rf:make-refiner 0 3)) (empty (csc-from-entries 0 0 '())))
      (assert (= 0 (rf-caller:refine-code refiner (rf-caller:matrix-operator empty)
                                          (rf-caller:poison-solver 0 0d0) (zeros 0) (zeros 0) 0d0 1d0))))
    (format t "REFINEMENT-OUTCOMES-PASS~%")))

(defun run-refinement-use-tests ()
  "The two uses the contract describes: a proof factor refined toward the unshifted matrix, and a
regularized augmented factor refined toward the unregularized one."
  (let ((state (sb-ext:seed-random-state 20261204)) (proof-steps 0) (saddle-steps 0) (cases 0))
    (loop :for (m n) :in '((2 4) (3 5) (4 6) (5 9)) :do
      (dotimes (repeat 4)
        (let* ((entries (loop :for i :below m :append
                          (loop :for j :below n
                                ;; A diagonal leading block, so that J J^T is at least 9 I.
                                :when (or (= i j) (and (>= j m) (< (random 100 state) 55)))
                                  :collect (list i j (if (= i j) (+ 3 (random 3 state)) (- (random 5 state) 2))))))
               (source (csc-from-entries m n entries))
               (j-exact (let ((a (make-array (list m n) :initial-element 0)))
                          (dolist (entry entries a) (setf (aref a (first entry) (second entry)) (third entry)))))
               (gram (build-matrix m (lambda (a b) (loop :for k :below n :sum (* (aref j-exact a k) (aref j-exact b k))))))
               (b (make-array m :element-type 'double-float)))
          (dotimes (i m) (setf (aref b i) (coerce (- (random 9 state) 4) 'double-float)))
          (incf cases)
          ;; J J^T y = b from the Cholesky-form factor of J J^T - shift I, refined against
          ;; J (J^T y) applied as two products.
          (multiple-value-bind (order pointers rows values) (dense-lower-pattern gram)
            (let* ((symbolic (ldlt:analyze (ldlt:pattern order pointers rows))) (factor (ldlt:symbolic-factor symbolic))
                   (refiner (rf:make-refiner m 30)) (y (zeros m)) (shift 1/8)
                   (tolerance (* 256 double-float-epsilon (max 1d0 (infinity-norm b)))))
              (assert (= 0 (ch-caller:observe-cholesky symbolic values (coerce shift 'double-float))))
              (assert (= 0 (rf-caller:refine-code refiner (rf-caller:gram-operator source (zeros n))
                                                  (rf-caller:cholesky-solver factor) b y tolerance 1d0))
                      () "The proof factor did not refine to J J^T: norms ~S."
                      (subseq (rf:residual-norms refiner) 0 (rf:residual-count refiner)))
              (setf proof-steps (max proof-steps (rf:step-count refiner)))
              (assert (plusp (rf:step-count refiner)) () "The shifted solve needed no correction.")
              (let* ((exact (exact-solution gram b))
                     (error (loop :for i :below m :maximize (abs (- (rational (aref y i)) (aref exact i))))))
                (assert (<= error (* (expt 2 -30) (max 1 (loop :for v :across exact :maximize (abs v))))) ()
                        "The refined proof solve is off by ~E." (float error 1d0)))))
          ;; The minimum-norm x with J x = b from the augmented factor at a visible regularization,
          ;; refined against [[I, J^T], [J, 0]].
          (let* ((plan (aug:analyze source)) (size (+ n m))
                 (symbolic (ldlt:analyze-quasidefinite
                            (ldlt:pattern size (aug:column-pointers plan) (aug:row-indices plan)) n))
                 (factor (ldlt:symbolic-factor symbolic))
                 (augmented (zeros (aug:entry-count plan)))
                 (rhs (zeros size)) (z (zeros size)) (refiner (rf:make-refiner size 40))
                 (tolerance (* 256 double-float-epsilon (max 1d0 (infinity-norm b)))))
            (dotimes (i m) (setf (aref rhs (+ n i)) (aref b i)))
            (assert (= 0 (qd-caller:raw-error plan (csc:entries source) 1d0 (scale-float 1d0 -6) augmented)))
            (assert (= 0 (qd-caller:observe-signed symbolic augmented 0d0)))
            (assert (= 0 (rf-caller:refine-code refiner
                                                (rf-caller:saddle-operator source (zeros n) (zeros m) (zeros n) (zeros m))
                                                (rf-caller:factor-solver factor) rhs z tolerance 1d0))
                    () "The regularized factor did not refine to the unregularized system: norms ~S."
                    (subseq (rf:residual-norms refiner) 0 (rf:residual-count refiner)))
            (setf saddle-steps (max saddle-steps (rf:step-count refiner)))
            (let* ((exact (exact-minimum-norm (rect m n (lambda (i j) (rational (coerce (aref j-exact i j) 'double-float))))
                                              (loop :for v :across b :collect (rational v))))
                   (error (loop :for j :below n :for value :in exact
                                :maximize (abs (- (rational (aref z j)) value)))))
              (assert (<= error (* (expt 2 -30) (max 1 (loop :for v :in exact :maximize (abs v))))) ()
                      "The refined minimum-norm solution is off by ~E." (float error 1d0)))))))
    (format t "REFINEMENT-USES-PASS cases=~D proof-steps<=~D saddle-steps<=~D~%" cases proof-steps saddle-steps)))

(defun run-refinement-tests ()
  (run-refinement-iteration-tests)
  (run-refinement-outcome-tests)
  (run-refinement-use-tests)
  ;; Warmed repeated work through a compiled Coalton caller allocates nothing: the operator and
  ;; the solve are made once, as a caller would make them.
  (let ((state (sb-ext:seed-random-state 5)) (n 6))
    (multiple-value-bind (dense matrix symbolic b) (test-system state n 1/16)
      (declare (ignore dense))
      (let* ((apply (rf-caller:matrix-operator matrix))
             (solve (rf-caller:factor-solver (ldlt:symbolic-factor symbolic)))
             (refiner (rf:make-refiner n 30)) (start (zeros n)) (x (zeros n))
             (tolerance (* 64 double-float-epsilon (infinity-norm b))))
        (assert (= 0 (rf-caller:refine-loop refiner apply solve b start x tolerance 1d0 3)) ()
                "The allocation case does not converge.")
        (allocation-check "refinement" (lambda () (rf-caller:refine-loop refiner apply solve b start x tolerance 1d0 1))))))
  (format t "REFINEMENT-TESTS-PASS~%"))
