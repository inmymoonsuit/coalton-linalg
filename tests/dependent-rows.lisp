(in-package #:coalton-linalg/tests)

;;; The pivot-skipping factorization and the dependent-rows module against exact rational
;;; arithmetic on the stored F64 values. The references share no code with the library: the
;;; exact pivots come from an elimination done here in the library's order, ranks and
;;; combinations from exact elimination and inversion, and every bound is compared with the
;;; exact quantity it bounds.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "ROWS" "COALTON-LINALG/SPARSE/DEPENDENT-ROWS")
  (sb-ext:add-package-local-nickname "DR-CALLER" "COALTON-LINALG/TESTS/DEPENDENT-ROWS-CALLERS"))

(defvar *sk-worst-entry* 0)
(defvar *sk-worst-solve* 0)
(defvar *sk-worst-null* 0)
(defvar *dr-worst-bound* 0)
(defvar *dr-worst-gap* 0)
(defvar *dr-worst-combination* 0)

;;; ---- The pivot-skipping factorization -------------------------------------------------

(defun exact-skipping-pivots (exact perm floor)
  "The exact pivots of P EXACT P^T when every index whose pivot is at most FLOOR is left out, and the list of
positions left out."
  (let* ((n (array-dimension exact 0))
         (work (build-matrix n (lambda (i j) (aref exact (aref perm i) (aref perm j)))))
         (pivots (make-array n)) (skipped '()))
    (dotimes (k n)
      (let ((pivot (aref work k k)))
        (setf (aref pivots k) pivot)
        (if (<= pivot floor)
            (push k skipped)
            (loop :for i :from (1+ k) :below n :do
              (let ((factor (/ (aref work i k) pivot)))
                (loop :for j :from (1+ k) :below n :do
                  (decf (aref work i j) (* factor (aref work k j)))))))))
    (values pivots (nreverse skipped))))

(defun check-skipping (label matrix &key (floor 1d-9) (shift 0) keep (state (sb-ext:seed-random-state 5)))
  "Factor a dense symmetric matrix of doubles, leaving out pivots at or below FLOOR, and hold every statement of
docs/skipping.md against exact arithmetic. Returns the number of indices left out."
  (multiple-value-bind (n pointers row-indices values) (symmetric-lower-pattern matrix :keep keep)
    (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers row-indices)))
           (factor (ldlt:symbolic-factor symbolic))
           (exact (stored-shifted matrix shift)))
      (assert (= 0 (ldlt:skipped-count factor)) () "~A: a skipped count before any factor." label)
      ;; A plain factor of the matrix made definite first: every pivot and every entry of the
      ;; storage then holds a number the skipping call has to replace.
      (expected-status (ldlt:factor! symbolic values (- (coerce (+ 1 (matrix-row-sum-norm exact)) 'double-float)))
                       (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess))
      (assert (notany #'zerop (ldlt:factor-diagonal factor)))
      (assert (= 0 (dr-caller:observe-skipping symbolic values (coerce shift 'double-float) floor)) ()
              "~A: the skipping factorization refused." label)
      (assert (and (ldlt:skipping-ready? factor) (not (ldlt:factor-ready? factor))
                   (not (ldlt:cholesky-ready? factor)))
              () "~A: the factor's state is wrong." label)
      (assert (= -1 (ldlt:refused-column symbolic)))
      (assert (and (= -1d0 (ldlt:factor-growth factor)) (= -1d0 (ldlt:factor-backward-bound factor))))
      (let* ((perm (ldlt:factor-permutation factor)) (lp (ldlt:factor-columns factor))
             (lr (ldlt:factor-rows factor)) (lv (ldlt:factor-values factor)) (d (ldlt:factor-diagonal factor))
             (kept (make-array n)) (stored (make-array (list n n) :initial-element nil))
             (row-count (make-array n :initial-element 0)) (widest 0))
        (dotimes (k n) (setf (aref kept k) (/= 0d0 (aref d k))))
        ;; The premise of every comparison below: no exact pivot lies between zero and far above
        ;; the floor, so the computed decisions are the exact ones.
        (multiple-value-bind (pivots skipped) (exact-skipping-pivots exact perm (rational floor))
          (dotimes (k n)
            (assert (or (<= (aref pivots k) (rational floor))
                        (> (aref pivots k) (* 1000 (max (rational floor) (expt 10 -12)))))
                    () "~A: the exact pivot ~E is too near the floor for this roster." label
                    (float (aref pivots k) 1d0))
            (assert (eq (aref kept k) (not (member k skipped))) ()
                    "~A: position ~D kept is ~S, and the exact pivot is ~E." label k (aref kept k)
                    (float (aref pivots k) 1d0)))
          (assert (= (length skipped) (ldlt:skipped-count factor)))
          (assert (equal (list (- n (length skipped)) 0) (pivot-counts factor)) ()
                  "~A: pivot counts ~S with ~D left out of ~D." label (pivot-counts factor) (length skipped) n)
          ;; A skipped index: a zero pivot and an empty column.
          (dotimes (k n)
            (setf (aref stored k k) t)
            (setf widest (max widest (- (aref lp (1+ k)) (aref lp k))))
            (loop :for q :from (aref lp k) :below (aref lp (1+ k)) :do
              (setf (aref stored (aref lr q) k) t)
              (incf (aref row-count (aref lr q)))
              (unless (aref kept k)
                (assert (= 0d0 (aref lv q)) () "~A: skipped column ~D holds ~S." label k (aref lv q)))))
          (dotimes (i n) (setf widest (max widest (aref row-count i))))
          (multiple-value-bind (product bound) (signed-factor-product factor n)
            ;; The kept principal submatrix, entry by entry, under the rounding analysis of the kernel.
            (dotimes (j n)
              (loop :for i :from j :below n :do
                (when (and (aref kept i) (aref kept j))
                  (let* ((products (loop :for k :below j
                                         :count (and (aref kept k) (aref stored i k) (aref stored j k))))
                         (residual (abs (- (aref exact (aref perm i) (aref perm j)) (aref product i j))))
                         (allowed (* (+ products 2) *unit-roundoff* (aref bound i j))))
                    (assert (<= residual (* allowed (+ 1 (expt 2 -40)))) ()
                            "~A: kept entry ~D,~D is off by ~E; ~E is allowed." label i j (float residual 1d0)
                            (float allowed 1d0))
                    (when (plusp allowed) (setf *sk-worst-entry* (max *sk-worst-entry* (/ residual allowed))))))))
            ;; The solve: zero at every skipped index, each kept row of the exact residual within
            ;; (3 w + 5) u |L||D||L^T| |x|, and nothing read from the skipped entries of b.
            (let ((b (zeros n)) (x (make-array n :element-type 'double-float :initial-element 7d0)))
              (dotimes (i n) (setf (aref b i) (coerce (- (random 9 state) 4) 'double-float)))
              (assert (= 0 (dr-caller:observe-skipping-solve factor b x)) () "~A: the solve refused." label)
              (dotimes (k n)
                (if (aref kept k)
                    (let* ((row (aref perm k))
                           (residual (abs (- (loop :for l :below n
                                                   :when (aref kept l)
                                                     :sum (* (aref exact row (aref perm l))
                                                             (rational (aref x (aref perm l)))))
                                             (rational (aref b row)))))
                           (allowed (* (+ (* 3 widest) 5) *unit-roundoff*
                                       (loop :for l :below n
                                             :sum (* (aref bound k l) (abs (rational (aref x (aref perm l)))))))))
                      (assert (<= residual (* allowed (+ 1 (expt 2 -40)))) ()
                              "~A: solve residual ~E in kept row ~D exceeds ~E." label (float residual 1d0) row
                              (float allowed 1d0))
                      (when (plusp allowed) (setf *sk-worst-solve* (max *sk-worst-solve* (/ residual allowed)))))
                    (assert (= 0d0 (aref x (aref perm k))) () "~A: the solution at skipped index ~D is ~S." label
                            (aref perm k) (aref x (aref perm k)))))
              (let ((other (copy-seq b)) (y (zeros n)))
                (dotimes (k n) (unless (aref kept k) (setf (aref other (aref perm k)) 123d0)))
                (assert (= 0 (dr-caller:observe-skipping-solve factor other y)))
                (assert (same-bits-p x y) () "~A: the solve read a skipped entry of the right-hand side." label)
                (assert (= 0 (dr-caller:observe-skipping-solve factor other other)))
                (assert (same-bits-p x other) () "~A: the solve in place differs." label)))
            ;; Each skipped index's vector: one at the index, zero at the other skipped indices and
            ;; at the kept ones eliminated after it, and M v at the kept ones before it within
            ;; (2 w + 4) u |L||D||L^T| |v|, the row of L the index left behind included.
            (dolist (k skipped)
              (let ((v (make-array n :element-type 'double-float :initial-element 7d0)))
                (assert (= 0 (dr-caller:observe-skipped-null symbolic (aref perm k) v)) ()
                        "~A: the vector of skipped index ~D was refused." label (aref perm k))
                (assert (= 1d0 (aref v (aref perm k))))
                (dotimes (j n)
                  (when (and (/= j k) (or (> j k) (not (aref kept j))))
                    (assert (= 0d0 (aref v (aref perm j))) ()
                            "~A: the vector of ~D is ~S at ~D, which is skipped or later." label (aref perm k)
                            (aref v (aref perm j)) (aref perm j))))
                (dotimes (j k)
                  (when (aref kept j)
                    (let* ((residual (abs (loop :for l :to k
                                                :sum (* (aref exact (aref perm j) (aref perm l))
                                                        (rational (aref v (aref perm l)))))))
                           (allowed (* (+ (* 2 widest) 4) *unit-roundoff*
                                       (loop :for l :to k
                                             :sum (* (aref bound j l) (abs (rational (aref v (aref perm l)))))))))
                      (assert (<= residual (* allowed (+ 1 (expt 2 -40)))) ()
                              "~A: M v for skipped ~D is ~E at kept ~D; ~E is allowed." label (aref perm k)
                              (float residual 1d0) (aref perm j) (float allowed 1d0))
                      (when (plusp allowed) (setf *sk-worst-null* (max *sk-worst-null* (/ residual allowed))))))))))
          ;; What does not accept this factor.
          (let ((b (zeros n)) (x (zeros n)))
            (assert (= 1 (dr-caller:observe-plain-solve factor b x)) () "~A: solve! accepted a skipping factor." label)
            (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveInvalidFactor)
                        (ldlt:solve-cholesky! factor b x)))
            (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:ConditionInvalidFactor)
                        (ldlt:estimate-inverse-norm1! factor))))
          (length skipped))))))

(defun seeded-semidefinite (state n rank density)
  "B B^T for an integer n-by-rank B: positive semidefinite, of rank at most RANK, with zero rows now and then."
  (let ((b (rect n rank (lambda (i j) (declare (ignore i j))
                          (if (< (random 100 state) density) (- (random 7 state) 3) 0)))))
    (build-matrix n (lambda (i j) (loop :for k :below rank :sum (* (aref b i k) (aref b j k)))))))

(defun run-skipping-factor-tests ()
  (setf *sk-worst-entry* 0 *sk-worst-solve* 0 *sk-worst-null* 0)
  (let ((state (sb-ext:seed-random-state 20261201)) (cases 0) (with-skips 0) (without 0) (most 0))
    (flet ((try (label matrix &rest keys)
             (let ((skipped (apply #'check-skipping label matrix :state state keys)))
               (incf cases)
               (if (plusp skipped) (incf with-skips) (incf without))
               (setf most (max most skipped))
               skipped)))
      (assert (= 0 (try "one by one" (rows-matrix '((4))))))
      (assert (= 1 (try "one zero" (rows-matrix '((0))))))
      (assert (= 1 (try "repeated index" (rows-matrix '((1 1) (1 1))))))
      (assert (= 2 (try "zero matrix" (rows-matrix '((0 0) (0 0))))))
      (assert (= 1 (try "sum of two" (rows-matrix '((2 0 2) (0 3 3) (2 3 5))))))
      (assert (= 0 (try "definite, nothing left out" (rows-matrix '((4 1 0) (1 3 1) (0 1 2))))))
      ;; The floor is inclusive: a pivot equal to it is left out, one just above is kept.
      (assert (= 1 (try "pivot at the floor" (rows-matrix '((4096 0 0) (0 1 0) (0 0 9216))) :floor 1d0)))
      (assert (= 0 (try "pivot above the floor" (rows-matrix '((4096 0 0) (0 2048 0) (0 0 9216))) :floor 1d0)))
      ;; The shift is taken before the pivot is judged.
      (assert (= 1 (try "shift makes a zero pivot" (rows-matrix '((4096 0 0) (0 1 0) (0 0 9216))) :shift 1 :floor 0d0)))
      ;; A negative pivot is at or below any floor: left out, not refused.
      (assert (= 1 (try "negative pivot" (rows-matrix '((4 0 0) (0 -5 0) (0 0 9))) :floor 0d0)))
      (loop :for n :from 2 :to 7 :do
        (loop :for rank :from 1 :to n :do
          (dotimes (repeat 3)
            (let ((matrix (seeded-semidefinite state n rank (if (evenp repeat) 70 45))))
              (try (format nil "semidefinite ~D rank<=~D #~D" n rank repeat) matrix
                   :keep (and (= repeat 2) (lambda (i j) (declare (ignore i j)) t))))))))
    (assert (and (> with-skips 40) (> without 15) (>= most 4)) ()
            "The skipping roster left indices out in ~D cases and none in ~D, ~D at most." with-skips without most)
    ;; With nothing at or below the floor the factor is the one factor! computes, bit for bit.
    (multiple-value-bind (n pointers row-indices values)
        (symmetric-lower-pattern (rows-matrix '((9 2 1 0) (2 8 0 3) (1 0 7 2) (0 3 2 6))))
      (let* ((plain (ldlt:analyze (ldlt:pattern n pointers row-indices)))
             (skipping (ldlt:analyze (ldlt:pattern n pointers row-indices))))
        (expected-status (ldlt:factor! plain values 0.5d0) (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess))
        (assert (= 0 (dr-caller:observe-skipping skipping values 0.5d0 0d0)))
        (assert (same-bits-p (ldlt:factor-diagonal (ldlt:symbolic-factor plain))
                             (ldlt:factor-diagonal (ldlt:symbolic-factor skipping))))
        (assert (same-bits-p (ldlt:factor-values (ldlt:symbolic-factor plain))
                             (ldlt:factor-values (ldlt:symbolic-factor skipping))))))
    (format t "SKIPPING-FACTOR-PASS cases=~D with-skips=~D most=~D worst-entry=~,3F worst-solve=~,3F worst-null=~,3F~%"
            cases with-skips most (float *sk-worst-entry* 1d0) (float *sk-worst-solve* 1d0)
            (float *sk-worst-null* 1d0))))

(defun run-skipping-refusal-tests ()
  (multiple-value-bind (n pointers row-indices values)
      (symmetric-lower-pattern (rows-matrix '((2 0 2) (0 3 3) (2 3 5))))
    (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers row-indices)))
           (factor (ldlt:symbolic-factor symbolic))
           (b (doubles 1d0 2d0 3d0)) (out (doubles 9d0 9d0 9d0)) (untouched (doubles 9d0 9d0 9d0))
           (nan (not-a-number)) (infinity sb-ext:double-float-positive-infinity))
      ;; Before any factor, and on the other kinds of factor.
      (assert (= 1 (dr-caller:observe-skipping-solve factor b out)))
      (assert (= 1 (dr-caller:observe-skipped-null symbolic 2 out)))
      (expected-status (ldlt:factor! symbolic values -1d0) (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess))
      (assert (= 1 (dr-caller:observe-skipping-solve factor b out)))
      (assert (= 1 (dr-caller:observe-skipped-null symbolic 2 out)))
      (assert (= 0 (ldlt:skipped-count factor)))
      (assert (same-bits-p out untouched))
      ;; A refused call leaves no factor of any kind, whatever stood before.
      (flet ((refused (code floor shift &optional (entries values))
               (expected-status (ldlt:factor! symbolic values -1d0)
                                (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess))
               (assert (= code (dr-caller:observe-skipping symbolic entries shift floor)) ()
                       "Floor ~S, shift ~S: expected refusal ~D." floor shift code)
               (assert (not (or (ldlt:skipping-ready? factor) (ldlt:factor-ready? factor))))
               (assert (= 0 (ldlt:skipped-count factor)))
               (assert (equal '(0 0) (pivot-counts factor)))))
        (refused 5 -1d0 0d0)
        ;; With the shift the matrix is definite: every pivot is positive, and only the check
        ;; of the floor itself can refuse.
        (refused 5 -1d0 -1d0)
        (refused 5 (- (scale-float 1d0 -1074)) 0d0)
        (refused 5 nan 0d0)
        (refused 5 infinity 0d0)
        (refused 4 0d0 nan)
        (refused 4 0d0 infinity)
        (refused 3 0d0 0d0 (doubles 1d0 2d0))
        (let ((poisoned (copy-seq values)))
          (setf (aref poisoned (1- (length poisoned))) nan)
          (refused 2 0d0 0d0 poisoned)
          (assert (<= 0 (ldlt:refused-column symbolic)) () "A non-finite value named no column.")))
      ;; On a standing skipping factor: sizes, a right-hand side that is not finite, indices.
      (assert (= 0 (dr-caller:observe-skipping symbolic values 0d0 1d-9)))
      (assert (= 1 (ldlt:skipped-count factor)))
      (assert (= 2 (dr-caller:observe-skipping-solve factor (doubles 1d0 2d0) out)))
      (assert (= 2 (dr-caller:observe-skipping-solve factor b (doubles 9d0 9d0))))
      (assert (= 3 (dr-caller:observe-skipping-solve factor (doubles 1d0 nan 3d0) out)))
      (assert (= 3 (dr-caller:observe-skipping-solve factor (doubles 1d0 2d0 infinity) out)))
      (assert (same-bits-p out untouched) () "A refused skipping solve wrote its output.")
      (let ((skipped (position 0d0 (ldlt:factor-diagonal factor)))
            (perm (ldlt:factor-permutation factor)))
        (dotimes (k 3)
          (assert (= (if (= k skipped) 0 2) (dr-caller:observe-skipped-null symbolic (aref perm k) (zeros 3))) ()
                  "Index ~D: the vector's outcome is wrong." (aref perm k))))
      (assert (= 2 (dr-caller:observe-skipped-null symbolic 3 out)))
      (assert (= 2 (dr-caller:observe-skipped-null symbolic 4000000 out)))
      (assert (= 3 (dr-caller:observe-skipped-null symbolic 0 (zeros 2))))
      (assert (same-bits-p out untouched) () "A refused null-vector call wrote its output.")))
  ;; A vector that overflows. The pivots are 1e-320, 1e280 and 1e306, the last index is left out
  ;; with a finite negative pivot, and the substitution for its vector multiplies 1e300 by 5e13.
  (let ((matrix (make-array '(4 4) :initial-element 0d0)))
    (flet ((put (i j v) (setf (aref matrix i j) v (aref matrix j i) v)))
      (put 0 0 1d-320) (put 1 0 1d-20) (put 1 1 2d280) (put 2 1 1d293) (put 2 2 2d306) (put 3 2 5d306))
    (multiple-value-bind (n pointers row-indices values) (symmetric-lower-pattern matrix)
      (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers row-indices)))
             (factor (ldlt:symbolic-factor symbolic)) (out (doubles 9d0 9d0 9d0 9d0)))
        (assert (= 0 (dr-caller:observe-skipping symbolic values 0d0 0d0)))
        (assert (equalp #(0 1 2 3) (ldlt:factor-permutation factor)) ()
                "The overflow case relies on the natural order; the order is ~S." (ldlt:factor-permutation factor))
        (assert (and (= 1 (ldlt:skipped-count factor)) (= 0d0 (aref (ldlt:factor-diagonal factor) 3))) ()
                "The overflow case's premise failed: pivots ~S." (ldlt:factor-diagonal factor))
        (assert (= 4 (dr-caller:observe-skipped-null symbolic 3 out)) ()
                "A vector that overflows was reported as ~S." out)
        (assert (same-bits-p out (doubles 9d0 9d0 9d0 9d0)) () "A non-finite vector was written."))))
  (format t "SKIPPING-REFUSALS-PASS~%"))

;;; ---- Dependent rows ---------------------------------------------------------------------

(defun dr-dense (m n entries)
  "The stored matrix as exact rationals of its doubles."
  (let ((a (rect m n (lambda (i j) (declare (ignore i j)) 0))))
    (dolist (entry entries)
      (setf (aref a (first entry) (second entry)) (rational (coerce (third entry) 'double-float))))
    a))

(defun dr-rank (a)
  (length (exact-basis-columns a)))

(defun dr-rows (a rows)
  (rect (length rows) (array-dimension a 1) (lambda (i j) (aref a (nth i rows) j))))

(defun dr-product-squared (a v)
  "The squared 2-norm of A^T v, exactly; V is a vector of doubles or rationals."
  (loop :for j :below (array-dimension a 1)
        :sum (expt (loop :for i :below (array-dimension a 0) :sum (* (aref a i j) (rational (aref v i)))) 2)))

(defun dr-gram-shifted (a shift)
  "A A^T - SHIFT I, exactly."
  (let ((m (array-dimension a 0)) (n (array-dimension a 1)))
    (build-matrix m (lambda (p q)
                      (- (loop :for j :below n :sum (* (aref a p j) (aref a q j))) (if (= p q) shift 0))))))

(defun dr-combination (a kept row)
  "The exact y on the KEPT rows nearest ROW in the least squares sense: (A_S A_S^T)^-1 A_S a_row."
  (let* ((s (dr-rows a kept)) (count (length kept)) (n (array-dimension a 1))
         (inverse (rational-inverse (dr-gram-shifted s 0)))
         (right (loop :for p :below count :collect (loop :for j :below n :sum (* (aref s p j) (aref a row j))))))
    (loop :for p :below count :collect (loop :for q :below count :sum (* (aref inverse p q) (nth q right))))))

(defun dr-distance-squared (a kept row)
  "The exact squared distance from ROW to the span of the KEPT rows."
  (let ((v (make-array (array-dimension a 0) :initial-element 0)))
    (setf (aref v row) 1)
    (loop :for i :in kept :for y :in (dr-combination a kept row) :do (setf (aref v i) (- y)))
    (dr-product-squared a v)))

(defun dr-left-out (analysis)
  (coerce (subseq (rows:left-out-rows analysis) 0 (rows:left-out-count analysis)) 'list))

(defun check-rows (label m n entries threshold &key exact-rank-expected (combination-tolerance (expt 2 -36)))
  "Find the dependent rows and hold every statement of docs/dependent-rows.md that does not depend on how near the
threshold a row lies. With EXACT-RANK-EXPECTED the kept rows must also be a basis and the combinations the exact ones."
  (let* ((matrix (csc-from-entries m n entries)) (a (dr-dense m n entries))
         (analysis (rows:analyze matrix)) (out (make-array m :element-type 'double-float :initial-element 7d0)))
    (assert (and (= m (rows:row-count analysis)) (= n (rows:column-count analysis))))
    (assert (and (not (rows:rows-found? analysis)) (= 0 (rows:kept-count analysis))
                 (= 0 (rows:left-out-count analysis)) (= -1d0 (rows:gap-bound analysis))))
    (assert (= 0 (dr-caller:observe-rows analysis matrix threshold)) () "~A: the search refused." label)
    (assert (rows:rows-found? analysis))
    (let* ((left (dr-left-out analysis))
           (kept (loop :for i :below m :unless (member i left) :collect i)) (sum 0))
      (assert (= (length left) (length (remove-duplicates left))))
      (assert (= (rows:kept-count analysis) (length kept)))
      (assert (= (+ (rows:kept-count analysis) (rows:left-out-count analysis)) m))
      (when exact-rank-expected
        (assert (= (dr-rank a) (length kept)) () "~A: ~D rows kept, and the exact rank is ~D." label (length kept)
                (dr-rank a))
        (assert (= (length kept) (dr-rank (dr-rows a kept))) () "~A: the kept rows ~S are not independent." label kept))
      ;; The kept values: the source's, with zeros for the rows left out.
      (let ((values (make-array (csc:entry-count matrix) :element-type 'double-float :initial-element 7d0))
            (ri (csc:row-indices matrix)) (source (csc:entries matrix)))
        (assert (= 1 (dr-caller:kept-values-code analysis values)))
        (dotimes (q (length values))
          (assert (= (aref values q) (if (member (aref ri q) left) 0d0 (aref source q))) ()
                  "~A: kept value ~D is ~S." label q (aref values q))))
      (assert (= (if left -1d0 0d0) (rows:gap-bound analysis)) () "~A: a gap bound before the vectors: ~S." label
              (rows:gap-bound analysis))
      (loop :for position :from 0 :for row :in left :do
        (assert (= -1d0 (rows:residual-bound analysis position)))
        (assert (= 0 (dr-caller:observe-null analysis position out)) () "~A: the vector of row ~D was refused." label row)
        (assert (= 1d0 (aref out row)))
        (dolist (other left)
          (unless (= other row)
            (assert (= 0d0 (aref out other)) () "~A: the vector of row ~D is ~S at left-out row ~D." label row
                    (aref out other) other)))
        ;; The bound holds for the exact product of the stored values and the vector as written.
        (let ((bound (rational (rows:residual-bound analysis position))) (exact (dr-product-squared a out)))
          (assert (<= exact (expt bound 2)) () "~A: row ~D: the residual is ~E and the bound ~E." label row
                  (sqrt (float exact 1d0)) (float bound 1d0))
          (when (plusp bound) (setf *dr-worst-bound* (max *dr-worst-bound* (/ exact (expt bound 2)))))
          (incf sum exact))
        (when exact-rank-expected
          (loop :for i :in kept :for y :in (dr-combination a kept row) :do
            (let ((off (abs (+ (rational (aref out i)) y))))
              (setf *dr-worst-combination* (max *dr-worst-combination* (/ off (max 1 (abs y)))))
              (assert (<= off (* combination-tolerance (max 1 (abs y)))) ()
                      "~A: row ~D on kept row ~D: multiplier ~S, exact ~E." label row i (- (aref out i))
                      (float y 1d0))))))
      ;; The gap bound: at least the root of the sum of the exact squares, and so at least
      ;; singular value r + 1, which the exact inertia of A A^T - g^2 I confirms.
      (let ((gap (rational (rows:gap-bound analysis))))
        (assert (>= gap 0) () "~A: no gap bound after every vector." label)
        (assert (<= sum (expt gap 2)) () "~A: the gap bound ~E is below the exact ~E." label (float gap 1d0)
                (sqrt (float sum 1d0)))
        (when (plusp gap) (setf *dr-worst-gap* (max *dr-worst-gap* (/ sum (expt gap 2)))))
        (when (and left (plusp m))
          (assert (<= (exact-inertia (dr-gram-shifted a (expt gap 2))) (length kept)) ()
                  "~A: more than ~D singular values exceed the gap bound ~E." label (length kept) (float gap 1d0))))
      (values analysis matrix a left kept))))

(defun with-combinations (state m n density extra)
  "Entries of an integer matrix whose last EXTRA rows are small integer combinations of the M rows before them."
  (let* ((base (seeded-entries state m n density))
         (a (dr-dense m n base)) (entries (copy-list base)))
    (dotimes (e extra)
      (let ((weights (loop :repeat m :collect (- (random 5 state) 2))))
        (dotimes (j n)
          (let ((value (loop :for i :below m :for w :in weights :sum (* w (aref a i j)))))
            (unless (zerop value) (push (list (+ m e) j value) entries))))))
    (reverse entries)))

(defun run-dependent-rows-exact-tests ()
  (setf *dr-worst-bound* 0 *dr-worst-gap* 0 *dr-worst-combination* 0)
  (let ((state (sb-ext:seed-random-state 20261202)) (cases 0) (with-left 0) (left-total 0))
    (flet ((try (label m n entries &optional (threshold 1d-4))
             (multiple-value-bind (analysis matrix a left) (check-rows label m n entries threshold :exact-rank-expected t)
               (declare (ignore analysis matrix a))
               (incf cases)
               (when left (incf with-left) (incf left-total (length left)))
               (length left))))
      (assert (= 0 (try "one row" 1 2 '((0 0 3) (0 1 4)))))
      (assert (= 1 (try "repeated row" 2 3 '((0 0 1) (0 1 2) (1 0 1) (1 1 2)))))
      (assert (= 1 (try "sum of rows" 3 4 '((0 0 1) (0 1 1) (1 2 1) (1 3 1) (2 0 1) (2 1 1) (2 2 1) (2 3 1)))))
      (assert (= 1 (try "more rows than columns" 3 2 '((0 0 1) (1 1 1) (2 0 1) (2 1 1)))))
      (assert (= 2 (try "two of four" 4 3 '((0 0 1) (0 1 2) (1 1 1) (1 2 1) (2 0 1) (2 1 3) (2 2 1) (3 0 2) (3 1 4)))))
      ;; A row with no stored entry, a row of stored zeros, and no columns at all.
      (assert (= 1 (try "empty row" 3 2 '((0 0 1) (2 1 1)))))
      (assert (= 1 (try "row of stored zeros" 2 2 '((0 0 1) (0 1 2) (1 0 0) (1 1 0)))))
      (assert (= 2 (try "no columns" 2 0 '())))
      (assert (= 0 (try "no rows" 0 3 '())))
      ;; A dense column, as a hub gives.
      (assert (= 2 (try "hub with a repeated spoke and a sum of two" 5 4
                        '((0 0 1) (0 1 1) (1 0 1) (1 2 1) (2 0 1) (2 3 1) (3 0 1) (3 1 1) (4 0 2) (4 1 1) (4 2 1)))))
      (loop :for (m n) :in '((2 3) (3 3) (3 5) (4 4) (4 6) (5 5) (5 8) (6 4) (7 9)) :do
        (dotimes (extra 4)
          (dotimes (repeat 3)
            (try (format nil "combinations ~Dx~D +~D #~D" m n extra repeat) (+ m extra) n
                 (with-combinations state m n (if (evenp repeat) 65 40) extra))))))
    (assert (and (> with-left 70) (> left-total 150)) ()
            "The exact roster left rows out in ~D cases, ~D rows in all." with-left left-total)
    (format t "DEPENDENT-ROWS-EXACT-PASS cases=~D with-left-out=~D rows-left-out=~D worst-bound=~,3F worst-gap=~,3F worst-multiplier=~,2E~%"
            cases with-left left-total (float *dr-worst-bound* 1d0) (float *dr-worst-gap* 1d0)
            (float *dr-worst-combination* 1d0))))

(defun run-dependent-rows-near-tests ()
  "Rows at a known distance from the others, with a threshold on each side of it, down to the method's floor."
  (let ((cases 0))
    ;; Row 1 is row 0 plus delta in a column of its own: its distance from row 0 is delta.
    (loop :for exponent :from -4 :downto -24 :by 2 :do
      (let* ((delta (scale-float 1d0 exponent))
             (entries (list '(0 0 1) '(0 1 1) '(1 0 1) '(1 1 1) (list 1 2 delta) '(2 3 5))))
        (multiple-value-bind (analysis matrix a left) (check-rows "near, kept" 3 4 entries (* 0.5d0 delta))
          (declare (ignore analysis matrix a))
          (assert (null left) () "Distance 2^~D with a threshold of half of it: rows ~S were left out." exponent left))
        (multiple-value-bind (analysis matrix a left kept) (check-rows "near, left out" 3 4 entries (* 2d0 delta))
          (declare (ignore matrix))
          (assert (= 1 (length left)) () "Distance 2^~D with a threshold of twice it: rows ~S were left out."
                  exponent left)
          ;; The bound reports the distance of the row left out, whichever of the two it is: at
          ;; least it, and within a part in a million of it.
          (let ((bound (rows:residual-bound analysis 0))
                (distance (sqrt (float (dr-distance-squared a kept (first left)) 1d0))))
            (assert (<= (* delta 0.99d0) distance (* delta 1.01d0)))
            (assert (<= (* distance (- 1d0 1d-12)) bound (* distance (+ 1d0 1d-6))) ()
                    "Distance 2^~D: exact ~S, bound ~S." exponent distance bound)
            (assert (<= bound (rows:gap-bound analysis) (* bound (+ 1d0 1d-12))))))
        (incf cases 2)))
    ;; Two rows near a third, both left out, each at its own distance: the gap bound is the
    ;; root of the sum of the squares, and neither alone.
    (let* ((delta (scale-float 1d0 -12))
           (entries (list '(0 0 1) '(0 1 1) '(1 0 1) '(1 1 1) (list 1 2 delta) '(2 0 1) '(2 1 1) (list 2 3 delta))))
      (multiple-value-bind (analysis matrix a left) (check-rows "two near rows" 3 4 entries (* 4d0 delta))
        (declare (ignore matrix a))
        (assert (= 2 (length left)) () "Two near rows: left out ~S." left)
        (let ((larger (max (rows:residual-bound analysis 0) (rows:residual-bound analysis 1))))
          (assert (> (rows:gap-bound analysis) (* 1.2d0 larger)) ()
                  "The gap bound ~S does not carry both residual bounds, the larger being ~S."
                  (rows:gap-bound analysis) larger))
        (incf cases)))
    ;; Below the floor of the method. The distance is 2^-30 and the threshold 2^-32, so the row
    ;; is independent at the caller's scale; in the product A A^T the distance is squared and
    ;; lost, and the row is left out. The bound, formed from A itself, reports the distance
    ;; that was lost, above the threshold: the caller can see what happened.
    (let* ((delta (scale-float 1d0 -30)) (threshold (scale-float 1d0 -32))
           (entries (list '(0 0 1) '(0 1 1) '(1 0 1) '(1 1 1) (list 1 2 delta))))
      (multiple-value-bind (analysis matrix a left kept) (check-rows "below the floor" 2 3 entries threshold)
        (declare (ignore matrix))
        (assert (= 1 (length left)) () "The floor case no longer shows the floor: rows left out ~S." left)
        (let ((distance (sqrt (float (dr-distance-squared a kept (first left)) 1d0))))
          (assert (< threshold (* 0.99d0 delta) distance))
          (assert (<= (* distance (- 1d0 1d-12)) (rows:residual-bound analysis 0) (* distance (+ 1d0 1d-6))) ()
                  "The bound ~S does not report the lost distance ~S." (rows:residual-bound analysis 0) distance))
        (incf cases)))
    ;; Refinement. Rows 0 and 1 differ by about 1e-5 of their size, so the kept pair is
    ;; conditioned near 1e5 and the factor of the squared problem gives the combination to five
    ;; or six digits. Row 2 is 3 row 0 - 2 row 1, rounded entry by entry, and is the nearest
    ;; thing to dependent that these values allow. The factor's own vector is measured here on
    ;; a separate analysis; null-vector! must be far nearer the exact least squares combination.
    (let* ((e 1d-5) (threshold 1d-7)
           (r0 (list 0.1d0 0.7d0 0.3d0))
           (r1 (mapcar (lambda (base direction) (+ base (* e direction))) r0 (list 0.3d0 -0.2d0 0.9d0)))
           (r2 (mapcar (lambda (p q) (- (* 3 p) (* 2 q))) r0 r1))
           (entries (append (loop :for i :from 0 :for row :in (list r0 r1 r2)
                                  :append (loop :for j :from 0 :for value :in row :collect (list i j value)))
                            '((3 3 1) (3 4 2)))))
      (multiple-value-bind (analysis matrix a left kept) (check-rows "refined" 4 5 entries threshold)
        (assert (and (= 1 (length left)) (< (first left) 3)) () "The refinement case left out ~S." left)
        (let* ((row (first left)) (exact (dr-combination a kept row))
               (plan (gram:analyze-aat matrix))
               (product (zeros (gram:entry-count plan)))
               (symbolic (ldlt:analyze (ldlt:pattern 4 (gram:column-pointers plan) (gram:row-indices plan))))
               (own (zeros 4)) (refined (zeros 4)))
          (gram:refill-csc! plan matrix product)
          (assert (= 0 (dr-caller:observe-skipping symbolic product 0d0 (* threshold threshold))))
          (assert (= 0 (dr-caller:observe-skipped-null symbolic row own)) ()
                  "The separate analysis did not leave out row ~D." row)
          (assert (= 0 (dr-caller:observe-null analysis 0 refined)))
          (flet ((off (vector)
                   (loop :for i :in kept :for y :in exact
                         :maximize (/ (abs (+ (rational (aref vector i)) y)) (max 1 (abs y))))))
            (assert (> (off own) (expt 2 -26)) ()
                    "The refinement case's premise failed: the factor's own vector is already within ~E."
                    (float (off own) 1d0))
            (assert (< (off refined) (expt 2 -32)) ()
                    "The refined combination is off by ~E; the factor's own was off by ~E."
                    (float (off refined) 1d0) (float (off own) 1d0))
            (format t "DEPENDENT-ROWS-REFINEMENT unrefined=~,2E refined=~,2E~%" (float (off own) 1d0)
                    (float (off refined) 1d0))))
        (incf cases)))
    ;; Refinement where one correction is not enough and the row left out is at a real
    ;; distance. The first two rows differ by 3e-7 of their size, so each correction gains only
    ;; about three digits; the third is their combination moved 6e-8 off the span, in a column
    ;; of its own, and is left out at a threshold of 1e-7. The row left out has a product with
    ;; its own residual that never falls, the square of its distance. If the refinement's
    ;; residual kept that entry, the first correction would look like a stall and the
    ;; combination would stop at six digits.
    (let* ((e 3d-7) (threshold 1d-7)
           (r0 (list 0.1d0 0.7d0 0.3d0))
           (r1 (mapcar (lambda (base direction) (+ base (* e direction))) r0 (list 0.3d0 -0.2d0 0.9d0)))
           (r2 (mapcar (lambda (p q) (- (* 3 p) (* 2 q))) r0 r1))
           (entries (append (loop :for i :from 0 :for row :in (list r0 r1 r2)
                                  :append (loop :for j :from 0 :for value :in row :collect (list i j value)))
                            '((2 5 6d-8) (3 3 1) (3 4 2)))))
      (multiple-value-bind (analysis matrix a left kept) (check-rows "refined, at a distance" 4 6 entries threshold)
        (assert (equal '(2) left) () "The distant refinement case left out ~S." left)
        (let* ((exact (dr-combination a kept 2))
               (plan (gram:analyze-aat matrix))
               (product (zeros (gram:entry-count plan)))
               (symbolic (ldlt:analyze (ldlt:pattern 4 (gram:column-pointers plan) (gram:row-indices plan))))
               (own (zeros 4)) (refined (zeros 4)))
          (gram:refill-csc! plan matrix product)
          (assert (= 0 (dr-caller:observe-skipping symbolic product 0d0 (* threshold threshold))))
          (assert (= 0 (dr-caller:observe-skipped-null symbolic 2 own)))
          (assert (= 0 (dr-caller:observe-null analysis 0 refined)))
          (flet ((off (vector)
                   (loop :for i :in kept :for y :in exact
                         :maximize (/ (abs (+ (rational (aref vector i)) y)) (max 1 (abs y))))))
            (assert (> (off own) (expt 2 -16)) ()
                    "The distant refinement case's premise failed: the factor's own vector is within ~E."
                    (float (off own) 1d0))
            (assert (< (off refined) (expt 2 -26)) ()
                    "The refined combination at a distance is off by ~E; the factor's own was off by ~E."
                    (float (off refined) 1d0) (float (off own) 1d0))
            (format t "DEPENDENT-ROWS-REFINEMENT-DISTANT unrefined=~,2E refined=~,2E~%" (float (off own) 1d0)
                    (float (off refined) 1d0)))))
      (incf cases))
    (format t "DEPENDENT-ROWS-NEAR-PASS cases=~D~%" cases)))

(defun certify-kept-rows (analysis matrix n top-scale bottom-scale)
  "The proof of docs/dependent-rows.md that every singular value of the kept rows exceeds sqrt(TOP-SCALE
BOTTOM-SCALE): the certificate of docs/quasidefinite.md on the plan of the whole matrix, with zeros written for the
rows left out, asking for as many negative pivots as rows kept. NIL when it does not prove it."
  (let* ((plan (aug:analyze matrix))
         (symbolic (ldlt:analyze-quasidefinite
                    (ldlt:pattern (aug:dimension plan) (aug:column-pointers plan) (aug:row-indices plan)) n))
         (factor (ldlt:symbolic-factor symbolic))
         (values (zeros (csc:entry-count matrix)))
         (out (zeros (aug:entry-count plan)))
         (kept (rows:kept-count analysis)) (left (rows:left-out-count analysis))
         (margin (* top-scale (scale-float 1d0 -44))))
    (assert (= 1 (dr-caller:kept-values-code analysis values)))
    (loop :repeat 4 :do
      (let ((top (+ top-scale margin)) (trailing (+ bottom-scale margin)))
        (unless (= 0 (qd-caller:raw-error plan values top (- trailing) out))
          (return-from certify-kept-rows nil))
        (if (= 0 (qd-caller:observe-signed symbolic out 0d0))
            (let ((needed (* (ldlt:factor-backward-bound factor) (+ 1d0 (scale-float 1d0 -50)))))
              (when (and (equal (list (+ n left) kept) (pivot-counts factor))
                         (>= (- top top-scale) needed) (>= (- trailing bottom-scale) needed))
                (return-from certify-kept-rows t))
              (setf margin (max (* 2 margin) (* 2 needed))))
            (setf margin (* 2 margin)))))
    nil))

(defun run-dependent-rows-certificate-tests ()
  "The kept rows are proved independent on the whole matrix's own plan, and never falsely."
  (let ((state (sb-ext:seed-random-state 20261203)) (cases 0) (issued 0) (issued-with-left 0) (withheld-false 0)
        (withheld-true 0))
    (flet ((try (label m n entries top bottom &optional (threshold 1d-4))
             (multiple-value-bind (analysis matrix a left kept) (check-rows label m n entries threshold)
               (let ((certified (certify-kept-rows analysis matrix n top bottom))
                     (true (or (null kept)
                               (= (length kept)
                                  (exact-inertia (dr-gram-shifted (dr-rows a kept) (* (rational top) (rational bottom))))))))
                 (incf cases)
                 (assert (or (not certified) true) ()
                         "~A: certified that the kept rows' sigma^2 exceeds ~E, which is false." label (* top bottom))
                 (cond ((and certified left) (incf issued) (incf issued-with-left))
                       (certified (incf issued))
                       (true (incf withheld-true))
                       (t (incf withheld-false)))
                 certified))))
      ;; A repeated row: the whole matrix has a zero singular value and is never certified;
      ;; the kept row, of squared norm 5, is.
      (assert (try "repeated row" 2 3 '((0 0 1) (0 1 2) (1 0 1) (1 1 2)) 2d0 1d0))
      (assert (not (try "repeated row, threshold above" 2 3 '((0 0 1) (0 1 2) (1 0 1) (1 1 2)) 2d0 3d0)))
      (assert (try "sum of rows" 3 4 '((0 0 1) (0 1 1) (1 2 1) (1 3 1) (2 0 1) (2 1 1) (2 2 1) (2 3 1)) 1d0 0.5d0))
      (assert (try "empty row" 3 2 '((0 0 1) (2 1 1)) 1d0 0.5d0))
      (loop :for (m n) :in '((2 3) (3 4) (3 6) (4 5) (5 7)) :do
        (dotimes (extra 3)
          (dotimes (repeat 4)
            (let* ((entries (with-combinations state m n 60 extra))
                   (top (coerce (max 1 (reduce #'max entries :key (lambda (e) (abs (third e))) :initial-value 0))
                                'double-float)))
              (dolist (bottom (list (scale-float 1d0 -40) 1d-6 1d-2 0.3d0 3d0 40d0))
                (try (format nil "combinations ~Dx~D +~D #~D bottom=~E" m n extra repeat bottom) (+ m extra) n entries
                     top (/ bottom top))))))))
    (assert (and (> issued 100) (> issued-with-left 50) (> withheld-false 40)) ()
            "The certificate roster issued ~D, ~D of them with rows left out, and rightly withheld ~D." issued
            issued-with-left withheld-false)
    ;; The case that fools the search. Row k of this matrix is column k of a Kahan matrix, at
    ;; distance 0.8^k from the rows before it, so with a threshold of 0.05 every row is kept;
    ;; the smallest singular value is far below the threshold all the same. The search cannot
    ;; see that, and the certificate at the threshold must be withheld.
    (let* ((size 10) (c 0.6d0) (s 0.8d0) (threshold 0.05d0)
           (entries (loop :for i :below size :append
                      (loop :for j :to i
                            :collect (list i j (* (expt s j) (if (= i j) 1d0 (- c))))))))
      (multiple-value-bind (analysis matrix a left) (check-rows "Kahan columns" size size entries threshold)
        (assert (null left) () "The fooling case no longer fools: rows ~S were left out." left)
        (assert (< (exact-inertia (dr-gram-shifted a (expt (rational threshold) 2))) size) ()
                "The fooling case's premise failed: every singular value exceeds the threshold.")
        (assert (not (certify-kept-rows analysis matrix size 1d0 (* threshold threshold))) ()
                "A kept set with a singular value below the threshold was certified at it.")))
    (format t "DEPENDENT-ROWS-CERTIFICATE-PASS cases=~D issued=~D issued-with-rows-left-out=~D withheld-though-true=~D withheld-false=~D~%"
            cases issued issued-with-left withheld-true withheld-false)))

(defun run-product-bound-tests ()
  "transpose-product-bound against the exact norm of the exact product."
  (let ((state (sb-ext:seed-random-state 20261204)) (cases 0) (worst 0) (rounding-mattered 0))
    (flet ((try (label m n entries v)
             (let* ((matrix (csc-from-entries m n entries)) (a (dr-dense m n entries))
                    (vector (apply #'doubles v))
                    (bound (rows:transpose-product-bound matrix vector))
                    (exact (dr-product-squared a vector))
                    ;; What a plain computation in the library's order gives.
                    (plain (loop :for j :below n
                                 :sum (expt (rational
                                             (let ((sum 0d0))
                                               (dotimes (i m)
                                                 (when (find-if (lambda (e) (and (= i (first e)) (= j (second e)))) entries)
                                                   (incf sum (* (coerce (aref a i j) 'double-float) (aref vector i)))))
                                               sum))
                                            2))))
               (incf cases)
               (assert (<= exact (expt (rational bound) 2)) () "~A: the norm is ~E and the bound ~E." label
                       (sqrt (float exact 1d0)) bound)
               (when (< plain exact) (incf rounding-mattered))
               (when (plusp bound) (setf worst (max worst (/ exact (expt (rational bound) 2)))))
               (values bound exact plain))))
      ;; The plain computation gives zero and the exact norm is one: only the rounding term saves the bound.
      (multiple-value-bind (bound exact plain) (try "absorbed" 3 1 '((0 0 1d16) (1 0 1) (2 0 -1d16)) '(1d0 1d0 1d0))
        (assert (and (= plain 0) (= exact 1) (>= bound 1d0)) () "The absorption case: plain ~S exact ~S bound ~S."
                plain exact bound))
      ;; A product that underflows to zero: the exact norm is 2^-1200.
      (multiple-value-bind (bound exact plain)
          (try "underflow" 1 1 (list (list 0 0 (scale-float 1d0 -600))) (list (scale-float 1d0 -600)))
        (assert (and (= plain 0) (plusp exact) (plusp bound))))
      ;; No rounding at all: the bound is within its own margin of the exact norm.
      (multiple-value-bind (bound exact) (try "exact" 2 2 '((0 0 3) (1 1 4)) '(1d0 1d0))
        (assert (= exact 25))
        (assert (<= 5d0 bound (* 5d0 (+ 1d0 1d-13))) () "A product without rounding is bounded by ~S." bound))
      ;; The same product with its factors at the ends of the range, where a square of either
      ;; cannot be held: the bound is still within its margin of the norm, 5 times the scales.
      (loop :for (label matrix-power vector-power) :in '(("large" 600 0) ("small" -600 0) ("large vector" 0 600)
                                                          ("small vector" 0 -600) ("opposite" 900 -900)
                                                          ("both small" -500 -500) ("both large" 500 500))
            :do (let ((norm (scale-float 5d0 (+ matrix-power vector-power))))
                  (multiple-value-bind (bound exact)
                      (try label 2 2 (list (list 0 0 (scale-float 3d0 matrix-power)) (list 1 1 (scale-float 4d0 matrix-power)))
                           (list (scale-float 1d0 vector-power) (scale-float 1d0 vector-power)))
                    (assert (= exact (expt (rational norm) 2)))
                    (assert (<= norm bound (* norm (+ 1d0 1d-13))) () "~A: a norm of ~S is bounded by ~S." label norm bound))))
      (try "empty" 2 0 '() '(1d0 2d0))
      ;; One entry of 1 and sixty of 1e-8, each in its own column: every square after the first
      ;; is below half a unit in the last place of the running sum and is rounded away, so the
      ;; computed sum of squares falls short of the exact one by sixty parts in 10^16.
      (multiple-value-bind (bound exact plain)
          (try "sum rounds down" 1 61 (cons '(0 0 1) (loop :for j :from 1 :to 60 :collect (list 0 j 1d-8))) '(1d0))
        (declare (ignore bound))
        (assert (and (= plain exact) (> exact (+ 1 (* 50 (expt 2 -53))))) ()
                "The rounding-down case's premise failed: exact ~S." exact))
      (loop :for (m n) :in '((1 1) (2 3) (3 2) (4 4) (5 3) (6 7)) :do
        (dotimes (repeat 12)
          (let ((entries (loop :for i :below m :append
                           (loop :for j :below n
                                 :when (< (random 100 state) 70)
                                   :collect (list i j (* (- (random 2001 state) 1000)
                                                         (expt 10d0 (- (random 9 state) 4)))))))
                (v (loop :repeat m :collect (/ (- (random 2001 state) 1000) 7d0))))
            (try (format nil "seeded ~Dx~D #~D" m n repeat) m n entries v)))))
    (assert (> rounding-mattered 10) () "Only ~D cases had a plain norm below the exact one." rounding-mattered)
    ;; A product too large to hold: infinity, never a number that is not one. Two such products
    ;; of opposite sign in one column have an exact sum of zero, and the bound is finite: what it
    ;; allows for is the rounding of the two, sixteen digits below their size.
    (let ((infinity sb-ext:double-float-positive-infinity))
      (assert (= infinity (rows:transpose-product-bound (csc-from-entries 1 1 '((0 0 1d300))) (doubles 1d10))))
      (let ((bound (rows:transpose-product-bound (csc-from-entries 2 1 '((0 0 1d300) (1 0 -1d300))) (doubles 1d10 1d10))))
        (assert (< 0d0 bound 1d296) () "Two products of 1e310 that cancel are bounded by ~S." bound))
      ;; A vector that is not finite gives infinity.
      (assert (= infinity (rows:transpose-product-bound (csc-from-entries 1 1 '((0 0 1))) (doubles infinity)))))
    (assert (= 92 (dr-caller:bound-code (csc-from-entries 2 2 '((0 0 1) (1 1 1))) (doubles 1d0))))
    (format t "PRODUCT-BOUND-PASS cases=~D rounding-mattered=~D worst=~,3F~%" cases rounding-mattered
            (float worst 1d0))))

(defun run-dependent-rows-lifetime-tests ()
  (let* ((entries '((0 0 1) (0 1 2) (1 1 1) (1 2 1) (2 0 1) (2 1 3) (2 2 1) (3 0 2) (3 1 4)))
         (matrix (csc-from-entries 4 3 entries))
         (analysis (rows:analyze matrix))
         (out (doubles 9d0 9d0 9d0 9d0)) (untouched (doubles 9d0 9d0 9d0 9d0))
         (values (zeros 9)) (nan (not-a-number)))
    ;; Before any search.
    (assert (= 1 (dr-caller:observe-null analysis 0 out)))
    (assert (= 0 (dr-caller:kept-values-code analysis values)))
    (assert (= -1d0 (rows:residual-bound analysis 0)))
    ;; A threshold that cannot be used leaves no set, whatever stood before.
    (dolist (threshold (list -1d0 nan sb-ext:double-float-positive-infinity))
      (assert (= 0 (dr-caller:observe-rows analysis matrix 1d-6)))
      (assert (rows:rows-found? analysis))
      (assert (= 2 (dr-caller:observe-rows analysis matrix threshold)) () "Threshold ~S was accepted." threshold)
      (assert (and (not (rows:rows-found? analysis)) (= 0 (rows:kept-count analysis))
                   (= 0 (rows:left-out-count analysis)) (= -1d0 (rows:gap-bound analysis))))
      (assert (= 1 (dr-caller:observe-null analysis 0 out)))
      (assert (= 0 (dr-caller:kept-values-code analysis values))))
    ;; The search is made on the matrix divided by the power of two at its largest entry, 4 here,
    ;; and compares squares. A positive threshold so small beside that entry that its square is
    ;; not a normal number is refused, and the smallest whose square is one is taken. A threshold
    ;; beyond every row, 1e200, whose own square is not finite, leaves every row out.
    (assert (= 0 (dr-caller:observe-rows analysis matrix 1d-6)))
    (assert (= 2 (dr-caller:observe-rows analysis matrix 1d-160)) () "A threshold whose square underflows was accepted.")
    (assert (and (not (rows:rows-found? analysis)) (= 0 (rows:left-out-count analysis))))
    (assert (= 0 (dr-caller:observe-rows analysis matrix 1d-153)) () "A threshold whose square is normal was refused.")
    (assert (= 0 (dr-caller:observe-rows analysis matrix 1d200)) () "A threshold beyond every row was refused.")
    (assert (= 4 (rows:left-out-count analysis)) () "A threshold beyond every row left out ~D of 4."
            (rows:left-out-count analysis))
    ;; A matrix whose stored values are all zero has no largest entry for a threshold to be small
    ;; beside: every threshold is beyond every row, and every row is left out, as for [1e-300].
    (dolist (case (list (list 1 1 '((0 0 0d0))) (list 2 3 '((0 0 0d0) (1 2 0d0))) (list 2 2 '())))
      (destructuring-bind (m n zero-entries) case
        (let* ((zero-matrix (csc-from-entries m n zero-entries)) (zero-analysis (rows:analyze zero-matrix)))
          (dolist (threshold (list 0d0 1d-300 1d-160 1d-6 1d200))
            (let ((code (dr-caller:observe-rows zero-analysis zero-matrix threshold)))
              (assert (and (= 0 code) (rows:rows-found? zero-analysis) (= m (rows:left-out-count zero-analysis))
                           (= 0 (rows:kept-count zero-analysis)))
                      () "A ~D by ~D matrix of zeros under a threshold of ~S: outcome ~D, ~D rows left out." m n threshold
                      code (rows:left-out-count zero-analysis)))))))
    (let* ((tiny (csc-from-entries 1 1 (list (list 0 0 1d-300)))) (tiny-analysis (rows:analyze tiny)))
      (assert (and (= 0 (dr-caller:observe-rows tiny-analysis tiny 1d-160)) (= 1 (rows:left-out-count tiny-analysis)))))
    ;; The same rows at other scales are searched as the rows themselves: under a threshold of
    ;; zero and under one scaled with them, the rows left out, the vector of the first and the
    ;; values handed back are those of the matrix as given, the values scaled.
    (flet ((scaled (power)
             (csc-from-entries 4 3 (mapcar (lambda (entry)
                                             (list (first entry) (second entry)
                                                   (scale-float (coerce (third entry) 'double-float) power)))
                                           entries)))
           (answer (source threshold power)
             (let ((kept (zeros 9)) (vector (zeros 4)))
               (assert (= 0 (dr-caller:observe-rows analysis source threshold)) ()
                       "At 2^~D the search under ~S gave ~D." power threshold
                       (dr-caller:observe-rows analysis source threshold))
               (assert (= 1 (dr-caller:kept-values-code analysis kept)))
               (when (plusp (rows:left-out-count analysis)) (assert (= 0 (dr-caller:observe-null analysis 0 vector))))
               (list (rows:left-out-count analysis)
                     (subseq (rows:left-out-rows analysis) 0 (rows:left-out-count analysis))
                     vector (map 'vector (lambda (value) (scale-float value (- power))) kept)))))
      (let ((exact (answer matrix 0d0 0)) (near (answer matrix 1d-6 0)))
        (assert (plusp (first near)) () "The fixture leaves no row out.")
        (dolist (power '(-1000 -600 600 1000))
          (assert (equalp exact (answer (scaled power) 0d0 power)) () "At 2^~D the search under zero differs." power)
          (assert (equalp near (answer (scaled power) (scale-float 1d-6 power) power)) ()
                  "At 2^~D the search under a scaled threshold differs." power)))
      ;; Rows 2^600 times smaller under a threshold that is not: every row is shorter than it.
      (assert (= 0 (dr-caller:observe-rows analysis (scaled -600) 1d-150)))
      (assert (= 4 (rows:left-out-count analysis)) () "Rows shorter than the threshold were kept: ~D left out."
              (rows:left-out-count analysis)))
    ;; One row 2^600 times smaller than the others: its square cannot be held beside theirs.
    ;; Under a threshold of zero that is refused by name, since the row would be left out as
    ;; dependent; under a threshold whose square can be held the row is shorter than the
    ;; threshold, and is rightly left out.
    (let ((mixed (csc-from-entries 4 3 (mapcar (lambda (entry)
                                                 (list (first entry) (second entry)
                                                       (if (= 3 (first entry))
                                                           (scale-float (coerce (third entry) 'double-float) -600)
                                                           (third entry))))
                                               entries))))
      (assert (= 0 (dr-caller:observe-rows analysis matrix 1d-6)))
      (assert (= 3 (dr-caller:observe-rows analysis mixed 0d0)) ()
              "A row whose square underflows was searched under a threshold of zero: code ~D, ~D left out."
              (dr-caller:observe-rows analysis mixed 0d0) (rows:left-out-count analysis))
      (assert (and (not (rows:rows-found? analysis)) (= 0 (rows:left-out-count analysis))))
      (assert (= 0 (dr-caller:observe-rows analysis mixed 1d-100)))
      (assert (find 3 (subseq (rows:left-out-rows analysis) 0 (rows:left-out-count analysis))) ()
              "The row shorter than the threshold was kept."))
    ;; A matrix of another pattern or shape is a typed error, and leaves no set.
    (assert (= 0 (dr-caller:observe-rows analysis matrix 1d-6)))
    (assert (= 91 (dr-caller:observe-rows analysis (csc-from-entries 4 3 (cons '(3 2 1) entries)) 1d-6)))
    (assert (not (rows:rows-found? analysis)))
    (assert (= 92 (dr-caller:observe-rows analysis (csc-from-entries 4 4 entries) 1d-6)))
    (assert (= 92 (dr-caller:observe-rows analysis (csc-from-entries 5 3 entries) 1d-6)))
    (assert (not (rows:rows-found? analysis)))
    ;; On a standing set.
    (assert (= 0 (dr-caller:observe-rows analysis matrix 1d-6)))
    (assert (= 2 (rows:left-out-count analysis)))
    (assert (= 3 (dr-caller:observe-null analysis 0 (zeros 3))))
    (assert (= 2 (dr-caller:observe-null analysis 2 out)))
    (assert (= 2 (dr-caller:observe-null analysis 4000000 out)))
    (assert (same-bits-p out untouched) () "A refused null-vector! wrote its output.")
    (assert (= 92 (dr-caller:kept-values-code analysis (zeros 8))))
    (assert (= -1d0 (rows:residual-bound analysis 2)))
    ;; The gap bound needs every vector, and a new search withdraws the bounds.
    (assert (= -1d0 (rows:gap-bound analysis)))
    (assert (= 0 (dr-caller:observe-null analysis 1 out)))
    (assert (and (= -1d0 (rows:residual-bound analysis 0)) (<= 0d0 (rows:residual-bound analysis 1))))
    (assert (= -1d0 (rows:gap-bound analysis)))
    (assert (= 0 (dr-caller:observe-null analysis 0 out)))
    (assert (<= 0d0 (rows:gap-bound analysis)))
    (assert (= 0 (dr-caller:observe-rows analysis matrix 1d-6)))
    (assert (and (= -1d0 (rows:residual-bound analysis 0)) (= -1d0 (rows:residual-bound analysis 1))
                 (= -1d0 (rows:gap-bound analysis))))
    ;; The same analysis on other values of the same pattern, and back: the same answers, bit for bit.
    (let ((first-vector (zeros 4)) (again (zeros 4))
          (independent (csc-from-entries 4 3 '((0 0 1) (0 1 2) (1 1 1) (1 2 1) (2 0 1) (2 1 3) (2 2 5) (3 0 2) (3 1 7)))))
      (assert (= 0 (dr-caller:observe-null analysis 0 first-vector)))
      (assert (= 0 (dr-caller:observe-rows analysis independent 1d-6)))
      (assert (= 1 (rows:left-out-count analysis)))
      (assert (= 0 (dr-caller:observe-rows analysis matrix 1d-6)))
      (assert (= 2 (rows:left-out-count analysis)))
      (assert (= 0 (dr-caller:observe-null analysis 0 again)))
      (assert (same-bits-p first-vector again) () "The same values gave another vector after other values."))
    ;; A row of 1e-160 beside one of 1e150, under a threshold of zero: its square cannot be held
    ;; beside the other's. Before that was refused for itself, as an underflow, this case reached
    ;; the factorization, whose first multiplier overflowed.
    (let* ((tiny (csc-from-entries 2 1 '((0 0 1d-160) (1 0 1d150))))
           (small (rows:analyze tiny)))
      (assert (= 3 (dr-caller:observe-rows small tiny 0d0)) ()
              "The row too small to square gave ~D with ~D rows left out." (dr-caller:observe-rows small tiny 0d0)
              (rows:left-out-count small))
      (assert (not (rows:rows-found? small))))
    ;; A budget too small for the product or the factor is the analyses' typed error.
    (assert (= 98 (dr-caller:analysis-code matrix 3)))
    (assert (= 0 (dr-caller:analysis-code matrix 100)))
    ;; Steady state: a compiled caller allocates nothing.
    (let ((kept (zeros 9)))
      (allocation-check "dependent-rows-loop"
                        (lambda () (assert (= 0 (dr-caller:rows-loop analysis matrix 1d-6 out kept 1))))))
    (multiple-value-bind (n pointers row-indices entries)
        (symmetric-lower-pattern (rows-matrix '((2 0 2) (0 3 3) (2 3 5))))
      (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers row-indices)))
             (rhs (doubles 1d0 2d0 3d0)) (vector (zeros 3)))
        (assert (= 0 (dr-caller:observe-skipping symbolic entries 0d0 1d-9)))
        (let ((skipped (aref (ldlt:factor-permutation (ldlt:symbolic-factor symbolic))
                             (position 0d0 (ldlt:factor-diagonal (ldlt:symbolic-factor symbolic))))))
          (allocation-check "skipping-loop"
                            (lambda ()
                              (assert (= 0 (dr-caller:skipping-loop symbolic entries 1d-9 skipped rhs vector 1)))))))))
  (format t "DEPENDENT-ROWS-LIFETIME-PASS~%"))

(defparameter *rows-public-control*
  "(declare inspect (coalton-linalg/sparse/dependent-rows:RowAnalysis -> UFix))
   (define (inspect analysis)
     (+ (coalton-linalg/sparse/dependent-rows:kept-count analysis)
        (coalton-linalg/sparse/dependent-rows:left-out-count analysis)))
   (declare probe (coalton-linalg/sparse/dependent-rows:RowAnalysis -> UFix))
   (define (probe analysis) (inspect analysis))")

(defparameter *rows-encapsulation-probes*
  '(("row-analysis-entries-access" .
     "(declare probe (coalton-linalg/sparse/dependent-rows:RowAnalysis -> arr:LispArray F64))
      (define (probe analysis) (.entries analysis))")
    ("row-analysis-symbolic-access" .
     "(declare probe (coalton-linalg/sparse/dependent-rows:RowAnalysis -> ldlt:Symbolic))
      (define (probe analysis) (.symbolic analysis))")
    ("row-analysis-unwrapping" .
     "(declare probe (coalton-linalg/sparse/dependent-rows:RowAnalysis -> UFix))
      (define (probe analysis)
        (match analysis ((coalton-linalg/sparse/dependent-rows:MkRowAnalysis storage) (.rows storage))))")
    ("skip-status-as-factor-status" .
     "(declare probe (ldlt:SkipFactorStatus -> ldlt:FactorStatus))
      (define (probe status) status)")
    ("null-vector-status-as-solve-status" .
     "(declare probe (ldlt:NullVectorStatus -> ldlt:SolveStatus))
      (define (probe status) status)")
    ("row-status-as-null-status" .
     "(declare probe (coalton-linalg/sparse/dependent-rows:RowStatus -> coalton-linalg/sparse/dependent-rows:NullStatus))
      (define (probe status) status)")))

(defun run-dependent-rows-encapsulation-tests ()
  "External Coalton callers reach the analysis only through its public operations, and the status types stay apart."
  (multiple-value-bind (accepted diagnostics)
      (coalton-linalg/tests/encapsulation::compile-probe "rows-public-control" *rows-public-control*)
    (assert accepted () "The valid dependent-rows external caller failed:~%~A" diagnostics))
  (dolist (probe *rows-encapsulation-probes*)
    (assert (not (coalton-linalg/tests/encapsulation::compile-probe (car probe) (cdr probe))) ()
            "External caller ~A compiled despite required rejection." (car probe)))
  (dolist (name '("MKROWANALYSIS" "ROWANALYSISSTORAGE"))
    (multiple-value-bind (symbol visibility) (find-symbol name "COALTON-LINALG/SPARSE/DEPENDENT-ROWS")
      (assert (and symbol (eq visibility ':internal)) ()
              "Storage constructor ~A must exist and remain internal; visibility is ~S." name visibility)))
  (format t "DEPENDENT-ROWS-ENCAPSULATION-PASS rejected=~D public-control=1~%"
          (length *rows-encapsulation-probes*)))

(defun run-dependent-rows-tests ()
  (run-skipping-factor-tests)
  (run-skipping-refusal-tests)
  (run-dependent-rows-exact-tests)
  (run-dependent-rows-near-tests)
  (run-dependent-rows-certificate-tests)
  (run-product-bound-tests)
  (run-dependent-rows-lifetime-tests)
  (run-dependent-rows-encapsulation-tests)
  (format t "DEPENDENT-ROWS-TESTS-PASS~%"))
