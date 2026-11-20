(in-package #:coalton-linalg/tests)

;;; The inverse 1-norm estimate against exact rational arithmetic on the stored F64 values.
;;; The references here share no code with the library: a dense rational Gauss-Jordan inverse.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "CONDITION-CALLER" "COALTON-LINALG/TESTS/CONDITION-CALLERS"))

(defun condition-success () (coalton:coalton coalton-linalg/sparse/spd-ldlt:ConditionSuccess))
(defun condition-invalid () (coalton:coalton coalton-linalg/sparse/spd-ldlt:ConditionInvalidFactor))
(defun condition-nonfinite () (coalton:coalton coalton-linalg/sparse/spd-ldlt:ConditionNonFiniteResult))
(defun factor-success () (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess))

(defun dense-lower (rows)
  "Lower CSC of a symmetric matrix given as a list of row lists, every lower entry stored."
  (let ((n (length rows)) (pointers (list 0)) (row-indices '()) (entries '()))
    (dotimes (j n)
      (loop for i from j below n do
        (push i row-indices)
        (push (coerce (nth j (nth i rows)) 'double-float) entries))
      (push (length row-indices) pointers))
    (values n (apply #'indices (reverse pointers)) (apply #'indices (reverse row-indices))
            (apply #'doubles (reverse entries)))))

(defun stored-rational (rows shift)
  "The matrix the library was given, minus shift on the diagonal, as exact rationals of the stored doubles."
  (let* ((n (length rows)) (matrix (make-array (list n n))))
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref matrix i j)
              (- (rational (coerce (nth (min i j) (nth (max i j) rows)) 'double-float))
                 (if (= i j) (rational shift) 0)))))
    matrix))

(defun rational-inverse (matrix)
  (let* ((n (array-dimension matrix 0)) (work (make-array (list n (* 2 n)) :initial-element 0)))
    (dotimes (i n)
      (dotimes (j n) (setf (aref work i j) (aref matrix i j)))
      (setf (aref work i (+ n i)) 1))
    (dotimes (column n)
      (let ((pivot (loop for r from column below n unless (zerop (aref work r column)) return r)))
        (unless pivot (error "The reference matrix is exactly singular."))
        (unless (= pivot column)
          (dotimes (j (* 2 n)) (rotatef (aref work pivot j) (aref work column j))))
        (let ((scale (aref work column column)))
          (dotimes (j (* 2 n)) (setf (aref work column j) (/ (aref work column j) scale))))
        (dotimes (r n)
          (unless (or (= r column) (zerop (aref work r column)))
            (let ((factor (aref work r column)))
              (dotimes (j (* 2 n))
                (decf (aref work r j) (* factor (aref work column j)))))))))
    (let ((inverse (make-array (list n n))))
      (dotimes (i n) (dotimes (j n) (setf (aref inverse i j) (aref work i (+ n j)))))
      inverse)))

(defun rational-norm1 (matrix)
  (let ((n (array-dimension matrix 0)))
    (loop for j below n maximize (loop for i below n sum (abs (aref matrix i j))))))

(defun estimate-for (rows shift)
  "Factor and estimate through the public interface. Returns estimate, norm, solves, and both exact norms."
  (multiple-value-bind (n pointers row-indices entries) (dense-lower rows)
    (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers row-indices)))
           (factor (ldlt:symbolic-factor symbolic)))
      (assert (eq (factor-success) (ldlt:factor! symbolic entries shift)))
      (assert (eq (condition-success) (ldlt:estimate-inverse-norm1! factor)))
      (let ((exact (stored-rational rows shift)))
        (values (ldlt:inverse-norm1-estimate factor) (ldlt:shifted-norm1 symbolic entries shift)
                (ldlt:estimate-solve-count factor)
                (rational-norm1 (rational-inverse exact)) (rational-norm1 exact))))))

(defun check-lower-bound (label rows shift slack floor)
  "The estimate never exceeds the exact inverse norm by more than SLACK, and reaches FLOOR of it."
  (multiple-value-bind (estimate norm solves exact-inverse exact-norm) (estimate-for rows shift)
    (let ((ratio (/ (rational estimate) exact-inverse)))
      (assert (<= ratio (+ 1 slack)) () "~A: estimate ~S exceeds the exact inverse norm ~S (ratio ~F)."
              label estimate (float exact-inverse 1d0) (float ratio 1d0))
      (assert (>= ratio floor) () "~A: estimate ~S is ~F of the exact inverse norm ~S, under the pinned ~F."
              label estimate (float ratio 1d0) (float exact-inverse 1d0) (float floor 1d0))
      (assert (= (rational norm) exact-norm) () "~A: shifted-norm1 ~S differs from the exact ~S."
              label norm (float exact-norm 1d0))
      (assert (<= 1 solves (coalton:coalton ldlt:maximum-estimate-solves)) () "~A made ~D solves." label solves)
      ;; An underestimate of the inverse norm makes the reciprocal an overestimate, never the reverse.
      (let ((reciprocal (ldlt:reciprocal-condition1 norm estimate)))
        (assert (>= (* (rational reciprocal) (+ 1 slack)) (/ 1 (* exact-norm exact-inverse))) ()
                "~A: reciprocal condition ~S is under the exact value." label reciprocal))
      (float ratio 1d0))))

(defun seeded-spd (state n)
  "B^T B + I for a small integer B: exactly representable and positive definite."
  (let ((b (loop repeat n collect (loop repeat n collect (- (random 7 state) 3)))))
    (loop for i below n collect
      (loop for j below n collect
        (+ (if (= i j) 1 0) (loop for row in b sum (* (nth i row) (nth j row))))))))

(defun run-condition-tests ()
  ;; Diagonal matrices: the ascent reaches the smallest entry's vertex, so the estimate is exact.
  (dolist (diagonal '((1 2 4 1/8) (8) (1/4 1/4 1/4) (1024 1 1/1024 2 3)))
    (let* ((n (length diagonal))
           (rows (loop for i below n collect (loop for j below n collect (if (= i j) (nth i diagonal) 0)))))
      (multiple-value-bind (estimate norm solves exact-inverse exact-norm) (estimate-for rows 0d0)
        (declare (ignore solves))
        (assert (= (rational estimate) exact-inverse) () "Diagonal estimate ~S is not exact." estimate)
        (assert (= (rational norm) exact-norm)))))
  ;; Two nearly parallel rows: the uniform start is almost a least-amplified direction, and only
  ;; the alternating probe sees the large part of the inverse.
  (let ((near '()))
    (dolist (power '(10 20 30))
      (let* ((e (expt 2 (- power))) (rows (list (list 1 (- 1 e)) (list (- 1 e) 1))))
        (push (check-lower-bound (format nil "near-parallel 2^-~D" power) rows 0d0 1/1000000 999/1000)
              near)))
    ;; A singular integer Gram matrix moved off singularity by a power of two. The computed factor
    ;; is then a perturbed matrix, so the comparison allows the relative error its conditioning admits.
    (dolist (power '(10 20 30))
      (let ((e (expt 2 (- power))))
        (push (check-lower-bound (format nil "near-singular 2^-~D" power)
                                 (list (list (+ 3 e) 1 4) (list 1 (+ 1 e) 2) (list 4 2 (+ 6 e)))
                                 0d0 1/1000 1/2)
              near)))
    (format t "CONDITION-NEAR-SINGULAR-PASS cases=~D ratios=~{~,6F~^ ~}~%" (length near) (reverse near)))
  ;; Fixed small families with exact references. The floors are pinned observations for these
  ;; matrices, not a theorem about the estimator.
  (let ((ratios '()))
    (push (check-lower-bound "tridiagonal 6"
                             (loop for i below 6 collect
                               (loop for j below 6 collect (cond ((= i j) 2) ((= 1 (abs (- i j))) -1) (t 0))))
                             0d0 1/1000000000 1/2)
          ratios)
    (push (check-lower-bound "shifted 3x3" '((4 -2 0) (-2 5 1) (0 1 3)) 1d0 1/1000000000 1/2) ratios)
    (let ((state (sb-ext:seed-random-state 20261119)) (cases 0))
      (loop for n from 2 to 7 do
        (dotimes (repeat 8)
          (incf cases)
          (push (check-lower-bound (format nil "seeded n=~D #~D" n repeat) (seeded-spd state n) 0d0
                                   1/1000000000 1/2)
                ratios)))
      (format t "CONDITION-EXACT-PASS cases=~D smallest-ratio=~,6F largest-ratio=~,6F~%"
              (+ cases 2) (reduce #'min ratios) (reduce #'max ratios))))
  ;; Scaling by a power of two scales every solve exactly, so the estimate scales exactly.
  (let ((rows '((4 -2 0) (-2 5 1) (0 1 3))))
    (multiple-value-bind (base) (estimate-for rows 0d0)
      (dolist (power '(-200 -3 7 200))
        (let ((scaled (loop for row in rows collect (loop for v in row collect (* v (expt 2 power))))))
          (multiple-value-bind (estimate) (estimate-for scaled 0d0)
            (assert (= estimate (scale-float base (- power))) ()
                    "Scaling by 2^~D gave ~S; expected ~S exactly." power estimate (scale-float base (- power)))))))
    ;; A shift applied by factor! is the same matrix as one stored with the shift already taken off.
    (multiple-value-bind (shifted) (estimate-for rows 1d0)
      (multiple-value-bind (stored) (estimate-for '((3 -2 0) (-2 4 1) (0 1 2)) 0d0)
        (assert (= shifted stored) () "Shifted estimate ~S differs from the stored-shift ~S." shifted stored))))
  (format t "CONDITION-IDENTITIES-PASS~%")
  ;; Lifetime: no estimate before a factor, none after a refusal, none surviving a refactor.
  (multiple-value-bind (n pointers row-indices entries) (dense-lower '((4 -2 0) (-2 5 1) (0 1 3)))
    (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers row-indices)))
           (factor (ldlt:symbolic-factor symbolic)))
      (assert (eq (condition-invalid) (ldlt:estimate-inverse-norm1! factor)))
      (assert (= -1d0 (ldlt:inverse-norm1-estimate factor)))
      (assert (= 0 (ldlt:estimate-solve-count factor)))
      (assert (eq (factor-success) (ldlt:factor! symbolic entries 0d0)))
      (assert (= -1d0 (ldlt:inverse-norm1-estimate factor)))
      (assert (eq (condition-success) (ldlt:estimate-inverse-norm1! factor)))
      (let ((first (ldlt:inverse-norm1-estimate factor)))
        (assert (> first 0d0))
        (assert (eq (condition-success) (ldlt:estimate-inverse-norm1! factor)))
        (assert (= first (ldlt:inverse-norm1-estimate factor)) () "A repeated estimate changed."))
      (assert (eq (factor-success) (ldlt:factor! symbolic entries 0d0)))
      (assert (= -1d0 (ldlt:inverse-norm1-estimate factor)) () "An estimate survived a refactor.")
      (assert (not (eq (factor-success) (ldlt:factor! symbolic entries 100d0))))
      (assert (eq (condition-invalid) (ldlt:estimate-inverse-norm1! factor)))
      (assert (= -1d0 (ldlt:inverse-norm1-estimate factor)))
      ;; shifted-norm1 refusals, through a compiled caller so the typed error is matched.
      (assert (= 0 (condition-caller:norm-error symbolic entries 0d0)))
      (assert (= 2 (condition-caller:norm-error symbolic (doubles 1d0) 0d0)))
      (assert (= 5 (condition-caller:norm-error symbolic entries sb-ext:double-float-positive-infinity)))
      (let ((bad (copy-seq entries)))
        (setf (aref bad 2) (- sb-ext:double-float-positive-infinity sb-ext:double-float-positive-infinity))
        (assert (= 5 (condition-caller:norm-error symbolic bad 0d0))))
      (let ((huge (copy-seq entries)))
        (fill huge most-positive-double-float)
        (assert (= 6 (condition-caller:norm-error symbolic huge 0d0))))))
  ;; A pivot so small that its reciprocal overflows: the factor stands, the estimate is refused.
  (let* ((symbolic (ldlt:analyze (ldlt:pattern 1 (indices 0 1) (indices 0))))
         (factor (ldlt:symbolic-factor symbolic)))
    (assert (eq (factor-success) (ldlt:factor! symbolic (doubles (scale-float 1d0 -1060)) 0d0)))
    (assert (eq (condition-nonfinite) (ldlt:estimate-inverse-norm1! factor)))
    (assert (= -1d0 (ldlt:inverse-norm1-estimate factor)))
    (assert (= 1 (ldlt:estimate-solve-count factor))))
  ;; The empty matrix: nothing to solve, and the reciprocal convention is 1.
  (let* ((symbolic (ldlt:analyze (ldlt:pattern 0 (indices 0) (indices))))
         (factor (ldlt:symbolic-factor symbolic)))
    (assert (eq (factor-success) (ldlt:factor! symbolic (doubles) 0d0)))
    (assert (eq (condition-success) (ldlt:estimate-inverse-norm1! factor)))
    (assert (= 0d0 (ldlt:inverse-norm1-estimate factor)))
    (assert (= 0d0 (ldlt:shifted-norm1 symbolic (doubles) 0d0)))
    (assert (= 1d0 (ldlt:reciprocal-condition1 0d0 0d0))))
  (assert (= 2d0 (ldlt:reciprocal-condition1 2d0 0.25d0)))
  (assert (= 0d0 (ldlt:reciprocal-condition1 2d0 -1d0)))
  (assert (= 0d0 (ldlt:reciprocal-condition1 0d0 5d0)))
  (assert (= 0d0 (ldlt:reciprocal-condition1 sb-ext:double-float-positive-infinity 1d0)))
  (assert (= 0d0 (ldlt:reciprocal-condition1 most-positive-double-float most-positive-double-float)))
  (format t "CONDITION-LIFETIME-PASS~%")
  ;; Warmed repeated work through a compiled Coalton caller allocates nothing.
  (multiple-value-bind (n pointers row-indices entries)
      (dense-lower (loop for i below 12 collect
                     (loop for j below 12 collect (cond ((= i j) 4) ((= 1 (abs (- i j))) -1) (t 0)))))
    (let ((symbolic (ldlt:analyze (ldlt:pattern n pointers row-indices))))
      (assert (= 0 (condition-caller:estimate-loop symbolic entries 0d0 3)))
      (allocation-check "condition-estimate"
                        (lambda () (condition-caller:estimate-loop symbolic entries 0.5d0 1)))))
  (format t "CONDITION-TESTS-PASS~%"))
