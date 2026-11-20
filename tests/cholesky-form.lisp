(in-package #:coalton-linalg/tests)

;;; The Cholesky-form numeric entry against exact rational arithmetic on the stored F64 values.
;;; Its purpose is a proof, so the test is the componentwise bound such a proof rests on: for
;;; every entry, the exact difference between the matrix and G G^T against the number of
;;; products that entry received.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "CH-CALLER" "COALTON-LINALG/TESTS/CHOLESKY-CALLERS"))

(defvar *ch-worst-entry* 0)
(defvar *ch-entries* 0)

(defun symmetric-lower-pattern (matrix &key keep)
  "Lower CSC of a dense symmetric matrix, storing the diagonal, every nonzero, and every position KEEP names."
  (let ((n (array-dimension matrix 0)) (pointers (list 0)) (rows '()) (values '()))
    (dotimes (j n)
      (loop :for i :from j :below n
            :when (or (= i j) (/= 0 (aref matrix i j)) (and keep (funcall keep i j)))
              :do (push i rows)
                 (push (coerce (aref matrix i j) 'double-float) values))
      (push (length rows) pointers))
    (values n (apply #'indices (reverse pointers)) (apply #'indices (reverse rows))
            (apply #'doubles (reverse values)))))

(defun cholesky-dense (factor n)
  "G as an exact rational lower triangular matrix in elimination order, and its structure as booleans."
  (let ((lp (ldlt:factor-columns factor)) (lr (ldlt:factor-rows factor)) (lv (ldlt:factor-values factor))
        (d (ldlt:factor-diagonal factor))
        (g (make-array (list n n) :initial-element 0))
        (stored (make-array (list n n) :initial-element nil)))
    (dotimes (k n)
      (setf (aref g k k) (rational (aref d k)) (aref stored k k) t)
      (loop :for q :from (aref lp k) :below (aref lp (1+ k)) :do
        (setf (aref g (aref lr q) k) (rational (aref lv q)) (aref stored (aref lr q) k) t)))
    (values g stored)))

(defun check-cholesky (label matrix &key (shift 0) keep)
  "Factor in Cholesky form through the public interface and check every entry's bound exactly."
  (multiple-value-bind (n pointers rows values) (symmetric-lower-pattern matrix :keep keep)
    (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers rows))) (factor (ldlt:symbolic-factor symbolic))
           (other (ldlt:analyze (ldlt:pattern n pointers rows)))
           (exact (build-matrix n (lambda (i j) (- (rational (coerce (aref matrix i j) 'double-float))
                                                   (if (= i j) (rational (coerce shift 'double-float)) 0))))))
      (assert (= 0 (ch-caller:observe-cholesky symbolic values (coerce shift 'double-float))) ()
              "~A: Cholesky-form factorization refused." label)
      (assert (and (ldlt:cholesky-ready? factor) (not (ldlt:factor-ready? factor))))
      (assert (= -1 (ldlt:refused-column symbolic)))
      ;; The proof factor is not a solving factor.
      (let ((out (make-array n :element-type 'double-float :initial-element 9d0)))
        (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveInvalidFactor)
                    (ldlt:solve! factor (make-array n :element-type 'double-float :initial-element 1d0) out)))
        (assert (every (lambda (v) (= v 9d0)) out))
        (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:ConditionInvalidFactor)
                    (ldlt:estimate-inverse-norm1! factor))))
      ;; The same order and structure as the solving factorization of the same pattern.
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess)
                  (ldlt:factor! other values (coerce shift 'double-float))))
      (let ((solving (ldlt:symbolic-factor other)))
        (assert (and (equalp (ldlt:factor-permutation factor) (ldlt:factor-permutation solving))
                     (equalp (ldlt:factor-columns factor) (ldlt:factor-columns solving))
                     (equalp (ldlt:factor-rows factor) (ldlt:factor-rows solving)))
                () "~A: the two factorizations differ in order or structure." label)
        ;; G is L sqrt(D) to rounding: a sanity relation, not the bound.
        (dotimes (k n)
          (assert (< (abs (- (expt (aref (ldlt:factor-diagonal factor) k) 2) (aref (ldlt:factor-diagonal solving) k)))
                     (* 1d-12 (aref (ldlt:factor-diagonal solving) k))))))
      (multiple-value-bind (g stored) (cholesky-dense factor n)
        (let ((perm (ldlt:factor-permutation factor)) (lr (ldlt:factor-rows factor)))
          (dotimes (k n) (assert (plusp (aref g k k))))
          ;; t'_j, the count a proof's bound uses: entries of G in row j left of the diagonal,
          ;; which a caller tallies from the row indices of the stored strict lower triangle.
          (let ((row-counts (make-array n :initial-element 0)))
            (loop :for i :across lr :do (incf (aref row-counts i)))
            (dotimes (j n)
              (assert (= (aref row-counts j) (loop :for k :below j :count (aref stored j k))))))
          (dotimes (j n)
            (loop :for i :from j :below n :do
              (let* ((products (loop :for k :below j :count (and (aref stored i k) (aref stored j k))))
                     (sum (loop :for k :to j :sum (* (aref g i k) (aref g j k))))
                     (magnitude (loop :for k :to j :sum (abs (* (aref g i k) (aref g j k)))))
                     (residual (abs (- (aref exact (aref perm i) (aref perm j)) sum)))
                     (allowed (* (+ products 2) *unit-roundoff* magnitude)))
                (when (aref stored i j)
                  (incf *ch-entries*)
                  (assert (<= residual (* allowed (+ 1 (expt 2 -40)))) ()
                          "~A: entry ~D,~D is off by ~E after ~D products; (products + 2) u |G||G^T| allows ~E."
                          label i j (float residual 1d0) products (float allowed 1d0))
                  (when (plusp magnitude)
                    (setf *ch-worst-entry*
                          (max *ch-worst-entry* (/ residual (* (+ products 2) *unit-roundoff* magnitude))))))
                ;; Outside the structure G is zero, and so must the matrix be.
                (unless (aref stored i j)
                  (assert (zerop (aref exact (aref perm i) (aref perm j))) ()
                          "~A: a nonzero entry lies outside the factor's structure." label)))))))
      (values factor symbolic))))

(defun seeded-sparse-definite (state n density)
  "A symmetric integer matrix with a random pattern and a strictly dominant positive diagonal."
  (let ((off (make-array (list n n) :initial-element 0)))
    (dotimes (i n)
      (dotimes (j i)
        (when (< (random 100 state) density)
          (let ((v (- (random 7 state) 3)))
            (setf (aref off i j) v (aref off j i) v)))))
    (build-matrix n (lambda (i j)
                      (if (= i j) (+ 1 (random 3 state) (loop :for k :below n :sum (abs (aref off i k)))) (aref off i j))))))

(defun run-cholesky-form-tests ()
  (setf *ch-worst-entry* 0 *ch-entries* 0)
  (let ((cases 0) (state (sb-ext:seed-random-state 20261127)))
    (flet ((check (label matrix &rest keys)
             (incf cases)
             (apply #'check-cholesky label (if (listp matrix) (rows-matrix matrix) matrix) keys)))
      ;; [[4,2],[2,5]] = G G^T with G = [[2,0],[1,2]], exactly.
      (let ((factor (check "exact" '((4 2) (2 5)))))
        (assert (equalp (ldlt:factor-diagonal factor) (doubles 2d0 2d0)))
        (assert (equalp (ldlt:factor-values factor) (doubles 1d0))))
      (check "empty" '())
      (check "one" '((9)))
      (check "definite" '((4 -2 0) (-2 5 1) (0 1 3)))
      (check "shifted" '((4 -2 0) (-2 5 1) (0 1 3)) :shift 1)
      (check "tridiagonal" (build-matrix 9 (lambda (i j) (cond ((= i j) 2) ((= 1 (abs (- i j))) -1) (t 0)))))
      ;; An arrow: one dense row and column, which the order leaves until last.
      (check "arrow" (build-matrix 8 (lambda (i j) (cond ((= i j) 9) ((or (= i 0) (= j 0)) 1) (t 0)))))
      ;; A stored zero is part of the structure and takes part as an exact zero.
      (check "stored zero" '((3 0 1) (0 2 0) (1 0 4)) :keep (lambda (i j) (and (= i 1) (= j 0))))
      (loop :for n :from 2 :to 10 :do
        (dolist (density '(20 50 100))
          (dotimes (repeat 3)
            (check (format nil "seeded n=~D density=~D #~D" n density repeat)
                   (seeded-sparse-definite state n density)))))
      ;; Dense positive definite matrices with irrational factors.
      (loop :for n :from 2 :to 7 :do
        (dotimes (repeat 3)
          (check (format nil "gram n=~D #~D" n repeat) (rows-matrix (seeded-spd state n)))))
      ;; A shift taken off the diagonal by the call is the matrix stored with it already taken off.
      (multiple-value-bind (n pointers rows values) (symmetric-lower-pattern (rows-matrix '((7 -2 0) (-2 8 1) (0 1 6))))
        (let ((shifted (ldlt:analyze (ldlt:pattern n pointers rows))))
          (assert (= 0 (ch-caller:observe-cholesky shifted values 3d0)))
          (multiple-value-bind (factor) (check-cholesky "stored shift" (rows-matrix '((4 -2 0) (-2 5 1) (0 1 3))))
            (assert (and (equalp (ldlt:factor-values factor) (ldlt:factor-values (ldlt:symbolic-factor shifted)))
                         (equalp (ldlt:factor-diagonal factor)
                                 (ldlt:factor-diagonal (ldlt:symbolic-factor shifted))))
                    () "A shift is not a stored shift.")))))
    (format t "CHOLESKY-FORM-BOUND-PASS cases=~D entries=~D worst-entry=~,3F of (products + 2) u |G||G^T|~%"
            cases *ch-entries* (float *ch-worst-entry* 1d0)))
  ;; Refusals: a pivot that is not positive, in a named column, and nothing left standing.
  (loop :for (rows shift column) :in '((((1 2) (2 1)) 0 1) (((1 1) (1 1)) 0 1) (((4 0) (0 9)) 4 0) (((-1)) 0 0))
        :do (multiple-value-bind (n pointers row-indices values) (symmetric-lower-pattern (rows-matrix rows))
             (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers row-indices)))
                    (factor (ldlt:symbolic-factor symbolic)))
               (assert (= 1 (ch-caller:observe-cholesky symbolic values (coerce shift 'double-float))) ()
                       "~S at shift ~D was not refused as non-positive." rows shift)
               (assert (= (aref (ldlt:factor-permutation factor) column) (ldlt:refused-column symbolic)))
               (assert (and (not (ldlt:cholesky-ready? factor)) (not (ldlt:factor-ready? factor)))))))
  (multiple-value-bind (n pointers rows values) (symmetric-lower-pattern (rows-matrix '((4 -2 0) (-2 5 1) (0 1 3))))
    (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers rows))) (factor (ldlt:symbolic-factor symbolic))
           (infinity sb-ext:double-float-positive-infinity)
           (out (doubles 9d0 9d0 9d0)))
      (assert (= 3 (ch-caller:observe-cholesky symbolic (doubles 1d0) 0d0)))
      (assert (= 4 (ch-caller:observe-cholesky symbolic values infinity)))
      (assert (= 2 (ch-caller:observe-cholesky symbolic (let ((bad (copy-seq values))) (setf (aref bad 1) infinity) bad)
                                                0d0)))
      ;; Each call replaces whatever stood: a solving factor is withdrawn by the proof call and
      ;; restored by the next solving call.
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess) (ldlt:factor! symbolic values 0d0)))
      (assert (and (ldlt:factor-ready? factor) (not (ldlt:cholesky-ready? factor))))
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:ConditionSuccess)
                  (ldlt:estimate-inverse-norm1! factor)))
      ;; A refused proof call withdraws a standing solving factor too.
      (assert (= 4 (ch-caller:observe-cholesky symbolic values infinity)))
      (assert (and (not (ldlt:factor-ready? factor)) (not (ldlt:cholesky-ready? factor))
                   (= -1d0 (ldlt:inverse-norm1-estimate factor)))
              () "A refused proof call left something standing.")
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess) (ldlt:factor! symbolic values 0d0)))
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:ConditionSuccess)
                  (ldlt:estimate-inverse-norm1! factor)))
      (assert (= 0 (ch-caller:observe-cholesky symbolic values 0d0)))
      (assert (and (ldlt:cholesky-ready? factor) (not (ldlt:factor-ready? factor))
                   (= -1d0 (ldlt:inverse-norm1-estimate factor)))
              () "The proof call left a solving factor or its estimate standing.")
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveInvalidFactor)
                  (ldlt:solve! factor (doubles 2d0 4d0 4d0) out)))
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess) (ldlt:factor! symbolic values 0d0)))
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess)
                  (ldlt:solve! factor (doubles 2d0 4d0 4d0) out)))
      (assert (close-vector-p out (doubles 1d0 1d0 1d0)))
      ;; Warmed repeated work through a compiled Coalton caller allocates nothing.
      (assert (= 0 (ch-caller:cholesky-loop symbolic values 0.5d0 3)))
      (allocation-check "cholesky-form" (lambda () (ch-caller:cholesky-loop symbolic values 0.5d0 1)))))
  (format t "CHOLESKY-FORM-TESTS-PASS~%"))
