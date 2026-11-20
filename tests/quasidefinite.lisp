(in-package #:coalton-linalg/tests)

;;; The signed-pivot sparse factorization and the augmented-system plan against exact rational
;;; arithmetic on the stored F64 values. The references share no code with the library: the
;;; augmented matrix is assembled densely from the source entries, the exported factor is
;;; multiplied out exactly, and the inertia comes from an exact symmetric elimination.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "AUG" "COALTON-LINALG/SPARSE/AUGMENTED")
  (sb-ext:add-package-local-nickname "QD-CALLER" "COALTON-LINALG/TESTS/QUASIDEFINITE-CALLERS"))

(defvar *qd-worst-reconstruction* 0)
(defvar *qd-worst-residual* 0)
(defvar *qd-largest-growth* 0)
(defvar *qd-latest-growth* 0)
(defvar *qd-worst-backward* 0)

(defun exact-backward-bound (factor n longest bound)
  "The contract's formula on the exported factor, in exact arithmetic: gamma(w + 2) trace(|L||D||L^T|) plus the
underflow term eta N (w + 2) (2 + 2 max|L|). BOUND is |L||D||L^T| and LONGEST the longest row of L."
  (let* ((reach (+ longest 2))
         (relative (/ (* reach *unit-roundoff*) (- 1 (* reach *unit-roundoff*))))
         (trace (loop :for i :below n :sum (aref bound i i)))
         (steepest (reduce #'max (ldlt:factor-values factor) :key (lambda (v) (abs (rational v))) :initial-value 0)))
    (+ (* relative trace)
       (* (rational least-positive-double-float) (* n reach (+ 2 (* 2 steepest)))))))

(defun exact-frobenius-squared (exact perm product n)
  "The squared Frobenius norm of P K P^T - L D L^T, exactly."
  (loop :for i :below n :sum
    (loop :for j :below n :sum (expt (- (aref exact (aref perm i) (aref perm j)) (aref product i j)) 2))))

(defun csc-from-entries (m n entries)
  "A checked CSC matrix from a list of (row column value); a stored zero stays structural."
  (let ((count (length entries)))
    (csc:from-coo m n
                  (make-array count :element-type '(unsigned-byte 62) :initial-contents (mapcar #'first entries))
                  (make-array count :element-type '(unsigned-byte 62) :initial-contents (mapcar #'second entries))
                  (make-array count :element-type 'double-float
                                    :initial-contents (mapcar (lambda (e) (coerce (third e) 'double-float)) entries))
                  (coalton:coalton coalton-linalg/sparse/csc:RejectDuplicates))))

(defun seeded-entries (state m n density &key hub)
  "Stored entries of an m-by-n integer matrix; HUB also stores the whole first column."
  (loop :for i :below m :append
    (loop :for j :below n
          :when (or (and hub (= j 0)) (< (random 100 state) density))
            :collect (list i j (let ((v (- (random 7 state) 3))) (if (zerop v) 1 v))))))

(defun augmented-exact (m n entries top bottom)
  "[[top I, A^T], [A, -bottom I]] as an exact rational matrix."
  (let ((k (make-array (list (+ n m) (+ n m)) :initial-element 0)))
    (dotimes (j n) (setf (aref k j j) (rational top)))
    (dotimes (i m) (setf (aref k (+ n i) (+ n i)) (- (rational bottom))))
    (dolist (entry entries)
      (destructuring-bind (i j v) entry
        (setf (aref k (+ n i) j) (rational (coerce v 'double-float))
              (aref k j (+ n i)) (rational (coerce v 'double-float)))))
    k))

(defun lower-pattern-of (matrix)
  "Every lower entry of a dense symmetric matrix as a lower CSC pattern with its values."
  (let ((n (array-dimension matrix 0)) (pointers (list 0)) (rows '()) (values '()))
    (dotimes (j n)
      (loop :for i :from j :below n :do
        (push i rows)
        (push (coerce (aref matrix i j) 'double-float) values))
      (push (length rows) pointers))
    (values n (apply #'indices (reverse pointers)) (apply #'indices (reverse rows))
            (apply #'doubles (reverse values)))))

(defun signed-factor-product (factor n)
  "L D L^T and |L| |D| |L|^T of the exported factor as exact rational matrices, in elimination order."
  (let ((lp (ldlt:factor-columns factor)) (lr (ldlt:factor-rows factor)) (lv (ldlt:factor-values factor))
        (d (ldlt:factor-diagonal factor))
        (l (build-matrix n (lambda (i j) (if (= i j) 1 0)))))
    (dotimes (k n)
      (loop :for q :from (aref lp k) :below (aref lp (1+ k)) :do
        (setf (aref l (aref lr q) k) (rational (aref lv q)))))
    (let ((block (build-matrix n (lambda (i j) (if (= i j) (rational (aref d i)) 0)))))
      (values (matrix-product (matrix-product l block) (matrix-transpose l))
              (matrix-product (matrix-product (matrix-map #'abs l) (matrix-map #'abs block))
                              (matrix-transpose (matrix-map #'abs l)))))))

(defun pivot-counts (factor)
  (list (ldlt:positive-pivot-count factor) (ldlt:negative-pivot-count factor)))

(defun check-signed (label exact n pointers row-indices values
                     &key (shift 0) leading (state (sb-ext:seed-random-state 11)))
  "Factor with signed pivots through the public interface; EXACT is the rational matrix that is factored.
With LEADING the analysis holds the indices from LEADING on, as analyze-quasidefinite does."
  (let* ((symbolic (if leading
                       (ldlt:analyze-quasidefinite (ldlt:pattern n pointers row-indices) leading)
                       (ldlt:analyze (ldlt:pattern n pointers row-indices))))
         (factor (ldlt:symbolic-factor symbolic)))
    (assert (= -1d0 (ldlt:factor-growth factor)) () "~A: a growth was reported before any factor." label)
    (assert (= -1d0 (ldlt:factor-backward-bound factor)) () "~A: a backward bound was reported before any factor." label)
    (assert (= 0 (qd-caller:observe-signed symbolic values (coerce shift 'double-float))) ()
            "~A: signed factorization refused." label)
    (assert (ldlt:factor-ready? factor))
    (assert (= -1 (ldlt:refused-column symbolic)))
    (let ((perm (ldlt:factor-permutation factor)) (lp (ldlt:factor-columns factor)) (lr (ldlt:factor-rows factor))
          (stored (make-array (list n n) :initial-element nil))
          (row-count (make-array n :initial-element 0)) (widest 0))
      (dotimes (k n)
        (setf (aref stored k k) t)
        (setf widest (max widest (- (aref lp (1+ k)) (aref lp k))))
        (loop :for q :from (aref lp k) :below (aref lp (1+ k)) :do
          (setf (aref stored (aref lr q) k) t)
          (incf (aref row-count (aref lr q)))))
      (dotimes (i n) (setf widest (max widest (aref row-count i))))
      (multiple-value-bind (product bound) (signed-factor-product factor n)
        ;; Each entry: the exact difference from L D L^T against (products + 2) u |L||D||L^T|,
        ;; the rounding analysis of this kernel, whose products multiply three numbers.
        (dotimes (j n)
          (loop :for i :from j :below n :do
            (let* ((products (loop :for k :below j :count (and (aref stored i k) (aref stored j k))))
                   (residual (abs (- (aref exact (aref perm i) (aref perm j)) (aref product i j))))
                   (allowed (* (+ products 2) *unit-roundoff* (aref bound i j))))
              (assert (<= residual (* allowed (+ 1 (expt 2 -40)))) ()
                      "~A: entry ~D,~D is off by ~E after ~D products; (products + 2) u |L||D||L^T| allows ~E."
                      label i j (float residual 1d0) products (float allowed 1d0))
              (when (plusp allowed)
                (setf *qd-worst-reconstruction* (max *qd-worst-reconstruction* (/ residual allowed)))))))
        ;; The reported growth is the exact one: the largest entry of |L||D||L^T| over the largest
        ;; magnitude of the matrix factored.
        (let* ((top (matrix-max exact))
               (growth (if (plusp top) (/ (matrix-max bound) top) 0)))
          (setf *qd-latest-growth* growth)
          (setf *qd-largest-growth* (max *qd-largest-growth* growth))
          (assert (<= (abs (- (rational (ldlt:factor-growth factor)) growth)) (* (expt 2 -40) (+ 1 growth))) ()
                  "~A: reported growth ~S, exact ~E." label (ldlt:factor-growth factor) (float growth 1d0)))
        ;; The backward bound holds for the exact Frobenius norm of P K P^T - L D L^T, and is the
        ;; contract's formula on the exported factor, rounded up and by no more than rounding.
        (let* ((reported (rational (ldlt:factor-backward-bound factor)))
               (off (exact-frobenius-squared exact perm product n))
               (formula (exact-backward-bound factor n (reduce #'max row-count :initial-value 0) bound)))
          (assert (<= off (expt reported 2)) ()
                  "~A: the backward error is ~E and the reported bound ~E." label (sqrt (float off 1d0))
                  (float reported 1d0))
          (assert (<= formula reported (+ (* formula (+ 1 (expt 2 -30))) (* 8 (rational least-positive-double-float))))
                  () "~A: reported bound ~E, formula ~E." label (float reported 1d0) (float formula 1d0))
          (setf *qd-worst-backward* (max *qd-worst-backward* (/ off (expt reported 2)))))
        (assert (= n (reduce #'+ (pivot-counts factor))) () "~A: pivot counts ~S do not sum to ~D." label
                (pivot-counts factor) n)
        (multiple-value-bind (positive negative zero) (exact-inertia exact)
          (when (zerop zero)
            (assert (equal (list positive negative) (pivot-counts factor)) ()
                    "~A: pivot counts ~S differ from the exact inertia ~S." label (pivot-counts factor)
                    (list positive negative))))
        ;; A solve: each row of the exact residual against (3 w + 5) u |L||D||L^T| |x|, w being the
        ;; longest row or column of L, which is the bound the contract states.
        (let ((b (make-array n :element-type 'double-float))
              (x (make-array n :element-type 'double-float :initial-element 0d0)))
          (dotimes (i n) (setf (aref b i) (coerce (- (random 9 state) 4) 'double-float)))
          (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess) (ldlt:solve! factor b x))
                  () "~A: solve refused." label)
          (dotimes (k n)
            (let* ((row (aref perm k))
                   (residual (abs (- (loop :for j :below n :sum (* (aref exact row j) (rational (aref x j))))
                                     (rational (aref b row)))))
                   (allowed (* (+ (* 3 widest) 5) *unit-roundoff*
                               (loop :for l :below n :sum (* (aref bound k l) (abs (rational (aref x (aref perm l)))))))))
              (assert (<= residual (* allowed (+ 1 (expt 2 -40)))) ()
                      "~A: solve residual ~E in row ~D exceeds ~E." label (float residual 1d0) row (float allowed 1d0))
              (when (plusp allowed)
                (setf *qd-worst-residual* (max *qd-worst-residual* (/ residual allowed)))))))))
    (values factor symbolic)))

(defun check-held-order (label factor m n entries)
  "Every source row that stores an entry is eliminated after at least one of the columns it stores an entry in."
  (let ((position (make-array (+ n m))) (perm (ldlt:factor-permutation factor)))
    (dotimes (k (+ n m)) (setf (aref position (aref perm k)) k))
    (dotimes (i m)
      (let ((columns (loop :for entry :in entries :when (= i (first entry)) :collect (second entry))))
        (when columns
          (assert (some (lambda (j) (< (aref position j) (aref position (+ n i)))) columns) ()
                  "~A: row ~D is eliminated before every column it touches." label i))))))

(defun full-row-rank-p (m n entries)
  "Whether the source matrix has exactly m independent rows, by exact elimination."
  (let ((transposed (make-array (list n m) :initial-element 0)))
    (dolist (entry entries)
      (setf (aref transposed (second entry) (first entry)) (rational (coerce (third entry) 'double-float))))
    (= m (length (exact-basis-columns transposed)))))

(defun check-augmented (label m n entries top bottom &key held)
  "Lay the augmented matrix out through the plan, compare it entry by entry with the dense assembly, and factor it."
  (let* ((matrix (csc-from-entries m n entries))
         (plan (aug:analyze matrix)) (size (+ n m))
         (exact (augmented-exact m n entries top bottom))
         (out (make-array (aug:entry-count plan) :element-type 'double-float :initial-element 9d0)))
    (assert (= size (aug:dimension plan)))
    (assert (= n (aug:leading-dimension plan)))
    (assert (= (+ size (length entries)) (aug:entry-count plan)))
    (assert (= 0 (qd-caller:matrix-error plan matrix (coerce top 'double-float) (coerce bottom 'double-float) out)))
    ;; Every stored position holds the exact entry, and every exact nonzero is stored.
    (let ((stored (make-array (list size size) :initial-element nil))
          (cp (aug:column-pointers plan)) (ri (aug:row-indices plan)))
      (dotimes (j size)
        (loop :for q :from (aref cp j) :below (aref cp (1+ j)) :do
          (let ((i (aref ri q)))
            (assert (>= i j) () "~A: entry above the diagonal at ~D,~D." label i j)
            (assert (not (aref stored i j)) () "~A: entry ~D,~D stored twice." label i j)
            (setf (aref stored i j) t)
            (assert (= (rational (aref out q)) (aref exact i j)) ()
                    "~A: plan value ~S at ~D,~D, exact ~S." label (aref out q) i j (aref exact i j)))))
      (dotimes (j size)
        (loop :for i :from j :below size :do
          (assert (or (aref stored i j) (zerop (aref exact i j))) ()
                  "~A: nonzero entry ~D,~D is missing from the plan." label i j)))
      (assert (= (count t (make-array (* size size) :displaced-to stored)) (aug:entry-count plan))))
    (multiple-value-bind (factor symbolic)
        (check-signed label exact size (aug:column-pointers plan) (aug:row-indices plan) out
                      :leading (and held n))
      (when held (check-held-order label factor m n entries))
      ;; A quasi-definite matrix has exactly as many positive pivots as leading unknowns, in any order.
      (when (and (plusp top) (plusp bottom))
        (assert (equal (list n m) (pivot-counts factor)) ()
                "~A: pivot counts ~S, expected ~S." label (pivot-counts factor) (list n m)))
      (values factor symbolic plan))))

(defun run-quasidefinite-exact-tests ()
  (setf *qd-worst-reconstruction* 0 *qd-worst-residual* 0 *qd-largest-growth* 0 *qd-worst-backward* 0)
  (let ((cases 0) (state (sb-ext:seed-random-state 20261126)))
    ;; Augmented systems of tall, square and wide sources, sparse and with a hub column, from
    ;; a well-separated regularization down to 2^-30.
    ;; Each system is analysed twice: by plain minimum degree, and with every row held until one
    ;; of its columns is eliminated. Both must reproduce the matrix; only the held order keeps
    ;; the factor small when the trailing block is small.
    (let ((plain-growth 0) (held-growth 0))
      (loop :for (m n) :in '((1 1) (1 3) (2 2) (3 2) (2 5) (4 4) (3 6) (6 3) (5 5) (4 7) (7 4) (6 8))
            :do (dolist (hub '(nil t))
                 (let* ((entries (seeded-entries state m n 45 :hub hub))
                        (independent (full-row-rank-p m n entries)))
                   ;; Dependent rows leave pivots of the size of the trailing block after
                   ;; cancellation, so a block at rounding level is kept for independent rows.
                   (dolist (bottom (append (list 1 (expt 2 -10) (expt 2 -20) (expt 2 -30))
                                           (and independent (list (expt 2 -50)))))
                     (let ((label (format nil "augmented ~Dx~D hub=~A bottom=2^~D" m n hub (round (log bottom 2)))))
                       (incf cases 2)
                       (check-augmented label m n entries 1 bottom)
                       (when (= bottom (expt 2 -30)) (setf plain-growth (max plain-growth *qd-latest-growth*)))
                       (check-augmented (concatenate 'string label " held") m n entries 1 bottom :held t)
                       (setf held-growth (max held-growth *qd-latest-growth*)))))))
      ;; Pinned observations on this roster, not theorems: plain minimum degree reaches a growth
      ;; near the reciprocal of the trailing block, and the held order stays small at every size of it.
      (assert (> plain-growth (expt 2 24)) () "Plain minimum degree grew only ~E at 2^-30." (float plain-growth 1d0))
      (assert (< held-growth (expt 2 12)) () "The held order grew ~E." (float held-growth 1d0))
      (format t "QUASIDEFINITE-GROWTH plain-at-2^-30=~,3E held-at-every-size=~,3E~%"
              (float plain-growth 1d0) (float held-growth 1d0)))
    ;; Dependent rows do not matter: the trailing block keeps every pivot away from zero.
    (dolist (bottom (list 1 (expt 2 -20)))
      (incf cases 2)
      (check-augmented "repeated rows" 3 4 '((0 0 1) (0 1 1) (1 2 1) (1 3 1) (2 0 1) (2 1 1)) 1 bottom)
      (check-augmented "repeated rows held" 3 4 '((0 0 1) (0 1 1) (1 2 1) (1 3 1) (2 0 1) (2 1 1)) 1 bottom
                       :held t))
    ;; Empty shapes, a stored zero, and other diagonal scalars.
    (incf cases 6)
    (check-augmented "empty row held" 2 2 '((0 0 1) (0 1 2)) 1 1/2 :held t)
    (check-augmented "no rows" 0 3 '() 1 1)
    (check-augmented "no columns" 3 0 '() 1 1)
    (check-augmented "nothing" 0 0 '() 1 1)
    (check-augmented "stored zero" 2 2 '((0 0 2) (1 0 0) (1 1 3)) 1 1/2)
    (check-augmented "other scalars" 2 3 '((0 0 1) (0 2 -2) (1 1 3)) 4 1/8)
    (format t "QUASIDEFINITE-AUGMENTED-PASS cases=~D worst-reconstruction=~,3F worst-solve=~,3F ~
               largest-growth=~,3E (fractions of the entrywise bounds; growth is max |L||D||L|^T over max |K|)~%"
            cases (float *qd-worst-reconstruction* 1d0) (float *qd-worst-residual* 1d0)
            (float *qd-largest-growth* 1d0))
    ;; Symmetric matrices with a strictly dominant diagonal of mixed signs: no pivot can vanish
    ;; in any order, and the inertia is the signs of the diagonal.
    (let ((dominant 0))
      (loop :for n :from 1 :to 8 :do
        (dotimes (repeat 6)
          (let* ((signs (loop :repeat n :collect (if (zerop (random 2 state)) 1 -1)))
                 (off (let ((a (make-array (list n n) :initial-element 0)))
                        (dotimes (i n)
                          (dotimes (j i)
                            (when (< (random 100 state) 50)
                              (let ((v (- (random 7 state) 3)))
                                (setf (aref a i j) v (aref a j i) v)))))
                        a))
                 (matrix (build-matrix n (lambda (i j)
                                           (if (= i j)
                                               (* (nth i signs) (1+ (loop :for k :below n :sum (abs (aref off i k)))))
                                               (aref off i j))))))
            (incf dominant)
            (multiple-value-bind (size pointers rows values) (lower-pattern-of matrix)
              (let ((factor (check-signed (format nil "dominant n=~D #~D" n repeat) matrix size pointers rows values)))
                (assert (equal (list (count 1 signs) (count -1 signs)) (pivot-counts factor)) ()
                        "Dominant inertia ~S for signs ~S." (pivot-counts factor) signs))))))
      (format t "QUASIDEFINITE-DOMINANT-PASS cases=~D~%" dominant))))

(defun copied-sparse-factor (factor)
  (list (copy-seq (ldlt:factor-permutation factor)) (copy-seq (ldlt:factor-values factor))
        (copy-seq (ldlt:factor-diagonal factor))))

(defun run-quasidefinite-identity-tests ()
  (let ((success (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess))
        (non-positive (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive)))
    ;; On a positive definite matrix the two entries share every bit, and all pivots count positive.
    (multiple-value-bind (n pointers rows values) (lower-pattern-of (rows-matrix '((4 -2 0) (-2 5 1) (0 1 3))))
      (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers rows))) (factor (ldlt:symbolic-factor symbolic)))
        (assert (eq success (ldlt:factor! symbolic values 0d0)))
        (assert (equal '(3 0) (pivot-counts factor)) () "factor! counts ~S." (pivot-counts factor))
        ;; Growth is measured by the signed entry only.
        (assert (= -1d0 (ldlt:factor-growth factor)) () "factor! reported a growth of ~S." (ldlt:factor-growth factor))
        (assert (= -1d0 (ldlt:factor-backward-bound factor)) () "factor! reported a backward bound.")
        (let ((positive (copied-sparse-factor factor)))
          (assert (= 0 (qd-caller:observe-signed symbolic values 0d0)))
          (assert (equal '(3 0) (pivot-counts factor)))
          (assert (equalp positive (copied-sparse-factor factor)) () "The signed entry changed a definite factor.")
          ;; The growth and the backward bound describe a signed factor only: each of the other
          ;; two entries withdraws them.
          (assert (and (plusp (ldlt:factor-growth factor)) (plusp (ldlt:factor-backward-bound factor))))
          (assert (eq success (ldlt:factor! symbolic values 0d0)))
          (assert (and (= -1d0 (ldlt:factor-growth factor)) (= -1d0 (ldlt:factor-backward-bound factor))) ()
                  "factor! kept a signed factor's growth or backward bound.")
          (assert (= 0 (qd-caller:observe-signed symbolic values 0d0)))
          (assert (eq success (ldlt:factor-cholesky! symbolic values 0d0)))
          (assert (and (= -1d0 (ldlt:factor-growth factor)) (= -1d0 (ldlt:factor-backward-bound factor))) ()
                  "factor-cholesky! kept a signed factor's growth or backward bound.")
          (assert (eq success (ldlt:factor! symbolic values 0d0))))
        ;; Negated, every pivot is negative: the positive definite entry refuses, the signed one does not.
        (let ((negated (map '(simple-array double-float (*)) #'- values)))
          (assert (eq non-positive (ldlt:factor! symbolic negated 0d0)))
          (assert (equal '(0 0) (pivot-counts factor)) () "A refused factor! kept counts ~S." (pivot-counts factor))
          (assert (= 0 (qd-caller:observe-signed symbolic negated 0d0)))
          (assert (equal '(0 3) (pivot-counts factor))))
        ;; A shift applied by the call is the matrix stored with the shift already subtracted.
        (assert (= 0 (qd-caller:observe-signed symbolic values 7d0)))
        (assert (equal '(0 3) (pivot-counts factor)) () "Shifted counts ~S." (pivot-counts factor))
        (let ((shifted (copied-sparse-factor factor)))
          (multiple-value-bind (n2 p2 r2 stored) (lower-pattern-of (rows-matrix '((-3 -2 0) (-2 -2 1) (0 1 -4))))
            (let* ((other (ldlt:analyze (ldlt:pattern n2 p2 r2))))
              (assert (= 0 (qd-caller:observe-signed other stored 0d0)))
              (assert (equalp shifted (copied-sparse-factor (ldlt:symbolic-factor other))) ()
                      "A shift is not a stored shift.")))))))
  ;; Growth is measured against the matrix that is factored: with a shift, the diagonal counts
  ;; as shifted. Here the stored diagonal is 100 and the factored one is 2.
  (multiple-value-bind (n pointers rows values) (lower-pattern-of (rows-matrix '((100 1 0) (1 100 1) (0 1 100))))
    (check-signed "shifted growth" (rows-matrix '((2 1 0) (1 2 1) (0 1 2))) n pointers rows values :shift 98))
  ;; Growth belongs to one factorization: refactoring other values on the same analysis gives
  ;; the growth a fresh analysis gives for those values.
  ;; One column read by three rows: plain minimum degree takes two of the rows before the
  ;; column, so a small trailing block gives a large growth and a block of 1 a small one.
  (let* ((entries '((0 0 1) (1 0 2) (2 0 3)))
         (matrix (csc-from-entries 3 1 entries)) (plan (aug:analyze matrix))
         (pattern (ldlt:pattern (aug:dimension plan) (aug:column-pointers plan) (aug:row-indices plan)))
         (reused (ldlt:analyze pattern)) (fresh (ldlt:analyze pattern))
         (small (make-array (aug:entry-count plan) :element-type 'double-float :initial-element 0d0))
         (large (make-array (aug:entry-count plan) :element-type 'double-float :initial-element 0d0)))
    (assert (= 0 (qd-caller:raw-error plan (csc:entries matrix) 1d0 (scale-float 1d0 -30) small)))
    (assert (= 0 (qd-caller:raw-error plan (csc:entries matrix) 1d0 1d0 large)))
    (assert (= 0 (qd-caller:observe-signed reused small 0d0)))
    (let ((first (ldlt:factor-growth (ldlt:symbolic-factor reused))))
      (assert (= 0 (qd-caller:observe-signed reused large 0d0)))
      (assert (= 0 (qd-caller:observe-signed fresh large 0d0)))
      (assert (> first (* 1000 (ldlt:factor-growth (ldlt:symbolic-factor fresh)))) ()
              "The two growths do not differ enough to tell them apart: ~S and ~S." first
              (ldlt:factor-growth (ldlt:symbolic-factor fresh)))
      (assert (= (ldlt:factor-growth (ldlt:symbolic-factor reused)) (ldlt:factor-growth (ldlt:symbolic-factor fresh)))
              () "A refactorization reported ~S where a fresh one reports ~S."
              (ldlt:factor-growth (ldlt:symbolic-factor reused)) (ldlt:factor-growth (ldlt:symbolic-factor fresh)))))
  ;; The inverse-norm estimate works from a signed factor as it does from a definite one.
  (multiple-value-bind (factor symbolic)
      (check-augmented "estimate" 2 3 '((0 0 1) (0 1 2) (1 1 1) (1 2 -1)) 1 1/4)
    (declare (ignore symbolic))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:ConditionSuccess)
                (ldlt:estimate-inverse-norm1! factor)))
    (let* ((exact (augmented-exact 2 3 '((0 0 1) (0 1 2) (1 1 1) (1 2 -1)) 1 1/4))
           (truth (rational-norm1 (rational-inverse exact)))
           (ratio (/ (rational (ldlt:inverse-norm1-estimate factor)) truth)))
      (assert (<= 1/3 ratio (+ 1 (expt 2 -40))) () "Estimate over the exact inverse norm is ~F." (float ratio 1d0))))
  (format t "QUASIDEFINITE-IDENTITIES-PASS~%"))

(defun run-quasidefinite-lifetime-tests ()
  ;; A zero pivot is the one pivot the signed entry refuses, and it names the column.
  (loop :for (rows column) :in '((((0 1) (1 0)) 0) (((1 1) (1 1)) 1) (((1 2) (2 4)) 1)) :do
    (multiple-value-bind (n pointers row-indices values) (lower-pattern-of (rows-matrix rows))
      (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers row-indices))) (factor (ldlt:symbolic-factor symbolic))
             (out (doubles 9d0 9d0)))
        (assert (= 1 (qd-caller:observe-signed symbolic values 0d0)) () "~S was not refused as a zero pivot." rows)
        (assert (= (aref (ldlt:factor-permutation factor) column) (ldlt:refused-column symbolic)) ()
                "~S refused at ~D." rows (ldlt:refused-column symbolic))
        (assert (not (ldlt:factor-ready? factor)))
        (assert (= -1d0 (ldlt:factor-growth factor)) () "A refused signed factorization kept a growth.")
        (assert (= -1d0 (ldlt:factor-backward-bound factor)) () "A refused signed factorization kept a backward bound.")
        (assert (equal '(0 0) (pivot-counts factor)) () "A refused signed factorization kept counts.")
        (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveInvalidFactor)
                    (ldlt:solve! factor (doubles 1d0 1d0) out)))
        (assert (every (lambda (v) (= v 9d0)) out)))))
  ;; Every other refusal withdraws a standing factor, its counts and its estimate.
  (multiple-value-bind (n pointers rows values) (lower-pattern-of (rows-matrix '((1 2 0) (2 -3 1) (0 1 5))))
    (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers rows))) (factor (ldlt:symbolic-factor symbolic))
           (infinity sb-ext:double-float-positive-infinity))
      (loop :for (attempt expected) :in (list (list (lambda () (qd-caller:observe-signed symbolic (doubles 1d0) 0d0)) 3)
                                            (list (lambda () (qd-caller:observe-signed symbolic values infinity)) 4)
                                            (list (lambda () (qd-caller:observe-signed symbolic values (not-a-number))) 4)
                                            (list (lambda () (let ((bad (copy-seq values)))
                                                               (setf (aref bad 1) infinity)
                                                               (qd-caller:observe-signed symbolic bad 0d0)))
                                                  2)
                                            ;; A finite entry whose square overflows in the update.
                                            (list (lambda () (let ((bad (copy-seq values)))
                                                               (setf (aref bad 1) 1d200)
                                                               (qd-caller:observe-signed symbolic bad 0d0)))
                                                  2))
            :do (assert (= 0 (qd-caller:observe-signed symbolic values 0d0)))
               (assert (equal '(2 1) (pivot-counts factor)) () "Counts ~S before a refusal." (pivot-counts factor))
               (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:ConditionSuccess)
                           (ldlt:estimate-inverse-norm1! factor)))
               (assert (= expected (funcall attempt)))
               (assert (and (not (ldlt:factor-ready? factor)) (equal '(0 0) (pivot-counts factor))
                            (= -1d0 (ldlt:inverse-norm1-estimate factor))
                            (= -1d0 (ldlt:factor-growth factor))
                            (= -1d0 (ldlt:factor-backward-bound factor)))
                       () "A refusal with code ~D left something standing." expected))))
  ;; The held analysis takes any leading dimension up to the size; at either end nothing
  ;; distinguishes the vertices and the order is plain minimum degree.
  (multiple-value-bind (n pointers rows values) (lower-pattern-of (rows-matrix '((1 2 0 1) (2 -3 1 0) (0 1 5 2) (1 0 2 -4))))
    (declare (ignore values))
    (let ((plain (ldlt:factor-permutation (ldlt:symbolic-factor (ldlt:analyze (ldlt:pattern n pointers rows))))))
      (dolist (leading '(0 4))
        (assert (equalp plain (ldlt:factor-permutation
                               (ldlt:symbolic-factor
                                (ldlt:analyze-quasidefinite (ldlt:pattern n pointers rows) leading))))
                () "A leading dimension of ~D changed the order." leading))
      (assert (= 0 (qd-caller:held-analysis-error (ldlt:pattern n pointers rows) 2)))
      (assert (= 2 (qd-caller:held-analysis-error (ldlt:pattern n pointers rows) 5)))))
  ;; Dependent rows with a trailing block at rounding level: the last row's pivot is what is
  ;; left of that block after cancellation. The call either refuses a zero pivot or completes;
  ;; when it completes with counts other than columns and rows, the counts are the warning.
  (let* ((entries '((0 0 1) (0 1 1) (1 0 1) (1 1 1) (2 0 3) (2 1 3)))
         (plan (aug:analyze (csc-from-entries 3 2 entries)))
         (symbolic (ldlt:analyze-quasidefinite
                    (ldlt:pattern (aug:dimension plan) (aug:column-pointers plan) (aug:row-indices plan)) 2))
         (out (make-array (aug:entry-count plan) :element-type 'double-float :initial-element 0d0))
         (outcomes '()))
    (dolist (power '(-10 -40 -50 -52 -60))
      (assert (= 0 (qd-caller:raw-error plan (csc:entries (csc-from-entries 3 2 entries)) 1d0
                                        (scale-float 1d0 power) out)))
      (let ((code (qd-caller:observe-signed symbolic out 0d0)))
        (assert (member code '(0 1)) () "Dependent rows at 2^~D gave code ~D." power code)
        (push (list power code (pivot-counts (ldlt:symbolic-factor symbolic))) outcomes)))
    (assert (equal '(0 (2 3)) (rest (assoc -10 outcomes))) () "Dependent rows at 2^-10: ~S." (assoc -10 outcomes))
    (format t "QUASIDEFINITE-DEPENDENT power, code, counts: ~S~%" (reverse outcomes)))
  ;; A hub: every row stores an entry in column 0 and in one column of its own. Held, each row
  ;; follows its own column and the hub comes at the end, so the factor has no fill at all. Taking every
  ;; column first, which the held rule must not degenerate into, would join all the rows.
  (let* ((m 12)
         (entries (loop :for i :below m :append (list (list i 0 1) (list i (1+ i) 2))))
         (plan (aug:analyze (csc-from-entries m (1+ m) entries)))
         (pattern (ldlt:pattern (aug:dimension plan) (aug:column-pointers plan) (aug:row-indices plan)))
         (held (ldlt:symbolic-factor (ldlt:analyze-quasidefinite pattern (aug:leading-dimension plan)))))
    (assert (= (aug:entry-count plan) (+ (aug:dimension plan) (length (ldlt:factor-rows held)))) ()
            "The held order filled the hub: ~D factor entries for ~D matrix entries."
            (+ (aug:dimension plan) (length (ldlt:factor-rows held))) (aug:entry-count plan))
    ;; The hub and the last row are left with one neighbour each and tie; either may come last.
    (assert (>= (position 0 (ldlt:factor-permutation held)) (- (aug:dimension plan) 2)) ()
            "The hub column is not eliminated at the end: ~S." (ldlt:factor-permutation held)))
  ;; The plan refuses before it writes: lengths, aliasing and non-finite inputs leave the output alone.
  (let* ((entries '((0 0 1) (0 2 -2) (1 1 3))) (matrix (csc-from-entries 2 3 entries))
         (plan (aug:analyze matrix)) (values (csc:entries matrix))
         (out (make-array (aug:entry-count plan) :element-type 'double-float :initial-element 9d0))
         (infinity sb-ext:double-float-positive-infinity))
    (flet ((untouched () (every (lambda (v) (= v 9d0)) out)))
      (assert (= 2 (qd-caller:raw-error plan (doubles 1d0 2d0) 1d0 1d0 out)))
      (assert (= 2 (qd-caller:raw-error plan values 1d0 1d0 (make-array 3 :element-type 'double-float))))
      (assert (= 5 (qd-caller:raw-error plan values infinity 1d0 out)))
      (assert (= 5 (qd-caller:raw-error plan values 1d0 (not-a-number) out)))
      (assert (= 5 (qd-caller:raw-error plan (doubles 1d0 infinity 3d0) 1d0 1d0 out)))
      (assert (untouched) () "A refused refill wrote ~S." out)
      ;; One array for both: the lengths must first agree, which needs a source with no leading or trailing block.
      (let* ((empty (csc-from-entries 0 0 '())) (empty-plan (aug:analyze empty)) (shared (doubles)))
        (assert (= 9 (qd-caller:raw-error empty-plan shared 1d0 1d0 shared))))
      ;; The checked entry compares the whole source pattern.
      (assert (= 2 (qd-caller:matrix-error plan (csc-from-entries 3 3 entries) 1d0 1d0 out)))
      (assert (= 1 (qd-caller:matrix-error plan (csc-from-entries 2 3 '((0 0 1) (0 1 -2) (1 1 3))) 1d0 1d0 out)))
      (assert (= 1 (qd-caller:matrix-error plan (csc-from-entries 2 3 '((0 0 1) (1 1 3))) 1d0 1d0 out)))
      ;; The same column pointers with one entry in another row: only the row indices differ.
      (assert (= 1 (qd-caller:matrix-error plan (csc-from-entries 2 3 '((1 0 1) (0 2 -2) (1 1 3))) 1d0 1d0 out)))
      (assert (untouched))
      (assert (= 0 (qd-caller:raw-error plan values 2d0 0.5d0 out)))
      (assert (equalp out (doubles 2d0 1d0  2d0 3d0  2d0 -2d0  -0.5d0 -0.5d0)) () "Refill wrote ~S." out))
    ;; The budget counts both diagonals and every source entry, and is refused before allocation.
    (assert (= 0 (qd-caller:analysis-error matrix 8 8)))
    (assert (= 8 (qd-caller:analysis-error matrix 7 8)))
    (assert (= 8 (qd-caller:analysis-error matrix 4 8)))
    (assert (= 8 (qd-caller:analysis-error matrix 0 8)))
    (assert (= 7 (qd-caller:analysis-error
                  (csc:from-csc (1+ (coalton:coalton aug:maximum-analysis-size)) 1 (indices 0 0) (indices) (doubles))
                  10 0))))
  (format t "QUASIDEFINITE-LIFETIME-PASS~%"))

(defparameter *augmented-public-control*
  "(declare inspect (coalton-linalg/sparse/augmented:AugmentedPlan -> UFix))
   (define (inspect plan)
     (+ (coalton-linalg/sparse/augmented:dimension plan) (coalton-linalg/sparse/augmented:entry-count plan)))
   (declare probe (coalton-linalg/sparse/augmented:AugmentedPlan -> UFix))
   (define (probe plan) (inspect plan))")

(defparameter *augmented-encapsulation-probes*
  '(("augmented-construction" .
     "(declare probe (Void -> coalton-linalg/sparse/augmented:AugmentedPlan))
      (define (probe)
        (coalton-linalg/sparse/augmented:AugmentedPlan
          0 0 (arr:make 1 (the UFix 0)) (arr:make 0 (the UFix 0)) 0
          (arr:make 1 (the UFix 0)) (arr:make 0 (the UFix 0))))")
    ("augmented-source-access" .
     "(declare probe (coalton-linalg/sparse/augmented:AugmentedPlan -> arr:LispArray UFix))
      (define (probe plan) (.source-columns plan))")
    ("signed-status-confusion" .
     "(declare probe (ldlt:SignedFactorStatus -> ldlt:FactorStatus))
      (define (probe status) status)")))

(defun run-quasidefinite-encapsulation-tests ()
  "External Coalton callers reach the plan only through its public operations, and the two status types stay apart."
  (multiple-value-bind (accepted diagnostics)
      (coalton-linalg/tests/encapsulation::compile-probe "augmented-public-control" *augmented-public-control*)
    (assert accepted () "The valid augmented external caller failed:~%~A" diagnostics))
  (dolist (probe *augmented-encapsulation-probes*)
    (assert (not (coalton-linalg/tests/encapsulation::compile-probe (car probe) (cdr probe))) ()
            "External caller ~A compiled despite required rejection." (car probe)))
  (dolist (name '("MKAUGMENTEDPLAN" "AUGMENTEDPLANSTORAGE"))
    (multiple-value-bind (symbol visibility) (find-symbol name "COALTON-LINALG/SPARSE/AUGMENTED")
      (assert (and symbol (eq visibility ':internal)) ()
              "Storage constructor ~A must exist and remain internal; visibility is ~S." name visibility)))
  (format t "QUASIDEFINITE-ENCAPSULATION-PASS rejected=~D public-control=1~%"
          (length *augmented-encapsulation-probes*)))

(defun run-quasidefinite-underflow-tests ()
  "Near the bottom of the range the relative bound alone is false, and the reported bound still holds."
  ;; [[a, b], [b, c]] with a = 2^-1060, b = (2^14 + 1) eta, c = 5 2^-1060. The multiplier is
  ;; 1 + 2^-14 exactly and its product with b is (2^14 + 2 + 2^-14) eta, which rounds to a
  ;; multiple of eta: the second pivot is off by 2^-14 eta, far more than gamma times anything here.
  (let* ((eta least-positive-double-float)
         (a (scale-float 1d0 -1060)) (b (* 16385d0 eta)) (c (scale-float 5d0 -1060))
         (exact (rows-matrix (list (list (rational a) (rational b)) (list (rational b) (rational c)))))
         (symbolic (ldlt:analyze (ldlt:pattern 2 (indices 0 2 3) (indices 0 1 1))))
         (factor (ldlt:symbolic-factor symbolic)))
    (assert (= 0 (qd-caller:observe-signed symbolic (doubles a b c) 0d0)))
    (assert (equal '(2 0) (pivot-counts factor)))
    (assert (equalp #(0 1) (ldlt:factor-permutation factor)))
    (multiple-value-bind (product bound) (signed-factor-product factor 2)
      (let ((off (abs (- (aref exact 1 1) (aref product 1 1))))
            (relative-only (* 3 *unit-roundoff* (aref bound 1 1)))
            (reported (rational (ldlt:factor-backward-bound factor))))
        (assert (= off (* (expt 2 -14) (rational eta))) () "The planned underflow did not occur: off by ~E."
                (float off 1d0))
        (assert (> off (* 1000 relative-only)) () "The relative bound was meant to fail here.")
        (assert (<= (exact-frobenius-squared exact (ldlt:factor-permutation factor) product 2) (expt reported 2))
                () "The reported bound ~E does not cover an underflow of ~E." (float reported 1d0) (float off 1d0))
        ;; The term is the stated one, a few dozen eta here and not a blanket allowance.
        (let ((formula (exact-backward-bound factor 2 1 bound)))
          (assert (<= formula reported (+ (* formula (+ 1 (expt 2 -30))) (* 8 (rational eta)))) ()
                  "The reported bound is ~F eta and the formula ~F eta." (float (/ reported (rational eta)) 1d0)
                  (float (/ formula (rational eta)) 1d0))))))
  (format t "QUASIDEFINITE-UNDERFLOW-PASS~%"))

(defun exactly-above-p (m n entries product)
  "Whether every squared singular value of the source exceeds PRODUCT, by exact inertia of J J^T - PRODUCT I."
  (let ((j (make-array (list m n) :initial-element 0)))
    (dolist (entry entries)
      (setf (aref j (first entry) (second entry)) (rational (coerce (third entry) 'double-float))))
    (multiple-value-bind (positive negative zero)
        (exact-inertia (build-matrix m (lambda (a b)
                                         (- (loop :for k :below n :sum (* (aref j a k) (aref j b k)))
                                            (if (= a b) product 0)))))
      (declare (ignore negative zero))
      (= positive m))))

(defun certify-singular-values (m n entries top-scale bottom-scale &key first-margin)
  "The worked use of docs/quasidefinite.md. T when the factorization proves that every singular value of the
source exceeds sqrt(TOP-SCALE BOTTOM-SCALE); NIL when it does not prove it, which is not a proof of the contrary.
FIRST-MARGIN replaces the first guess at the margin, which is otherwise a 2^-44 part of TOP-SCALE."
  (let* ((matrix (csc-from-entries m n entries)) (plan (aug:analyze matrix))
         (symbolic (ldlt:analyze-quasidefinite
                    (ldlt:pattern (aug:dimension plan) (aug:column-pointers plan) (aug:row-indices plan)) n))
         (factor (ldlt:symbolic-factor symbolic))
         (out (make-array (aug:entry-count plan) :element-type 'double-float :initial-element 0d0))
         (margin (or first-margin (* top-scale (scale-float 1d0 -44)))))
    (loop :repeat 4 :do
      (let ((top (+ top-scale margin)) (trailing (+ bottom-scale margin)))
        ;; refill! writes minus its last scalar on the trailing diagonal; the proof wants plus.
        (unless (= 0 (qd-caller:raw-error plan (csc:entries matrix) top (- trailing) out))
          (return-from certify-singular-values nil))
        (if (= 0 (qd-caller:observe-signed symbolic out 0d0))
            (let ((needed (* (ldlt:factor-backward-bound factor) (+ 1d0 (scale-float 1d0 -50)))))
              (when (and (equal (list n m) (pivot-counts factor))
                         (>= (- top top-scale) needed) (>= (- trailing bottom-scale) needed))
                (return-from certify-singular-values t))
              (setf margin (max (* 2 margin) (* 2 needed))))
            (setf margin (* 2 margin)))))
    nil))

(defun run-quasidefinite-certificate-tests ()
  "A certificate is never issued for a false statement, and is issued where the margin is plain."
  (let ((state (sb-ext:seed-random-state 20261130)) (issued 0) (withheld-true 0) (withheld-false 0) (cases 0))
    (flet ((try (label m n entries top bottom)
             (let ((certified (certify-singular-values m n entries top bottom))
                   (true (exactly-above-p m n entries (* (rational top) (rational bottom)))))
               (incf cases)
               (assert (or (not certified) true) ()
                       "~A: certified that sigma^2 exceeds ~E, which is false." label (* top bottom))
               (cond (certified (incf issued)) (true (incf withheld-true)) (t (incf withheld-false)))
               certified)))
      ;; Plain margins: [I | B] has every singular value at least 1.
      (assert (try "identity block" 3 5 '((0 0 1) (1 1 1) (2 2 1) (0 3 2) (1 4 -1) (2 4 3)) 3d0 (/ 0.25d0 3d0)))
      (assert (try "one row" 1 2 '((0 0 3) (0 1 4)) 4d0 1d0))
      (assert (not (try "one row, threshold above" 1 2 '((0 0 3) (0 1 4)) 4d0 7d0)))
      ;; A threshold exactly at sigma^2 = 25 is not exceeded.
      (assert (not (try "one row, threshold at" 1 2 '((0 0 3) (0 1 4)) 5d0 5d0)))
      ;; The signs alone are not a proof. For A = [1 1 5] the squared singular value is 27; with
      ;; t = 3 and s = 9 the statement that it exceeds s t = 27 is false, and yet the computed
      ;; pivots of [[3 I, A^T], [A, 9]] are three positive and one negative, as if it were true:
      ;; the last pivot comes out as a small negative number where the exact one is zero.
      (let* ((entries '((0 0 1) (0 1 1) (0 2 5)))
             (matrix (csc-from-entries 1 3 entries)) (plan (aug:analyze matrix))
             (symbolic (ldlt:analyze-quasidefinite
                        (ldlt:pattern (aug:dimension plan) (aug:column-pointers plan) (aug:row-indices plan)) 3))
             (out (make-array (aug:entry-count plan) :element-type 'double-float :initial-element 0d0)))
        (assert (not (exactly-above-p 1 3 entries 27)))
        (assert (= 0 (qd-caller:raw-error plan (csc:entries matrix) 3d0 -9d0 out)))
        (assert (= 0 (qd-caller:observe-signed symbolic out 0d0)))
        (assert (equal '(3 1) (pivot-counts (ldlt:symbolic-factor symbolic))) ()
                "The trap is not one: the counts without a margin are ~S."
                (pivot-counts (ldlt:symbolic-factor symbolic)))
        (assert (not (try "signs alone would certify" 1 3 entries 3d0 9d0)))
        ;; A first margin below the rounding of 3 and 9 writes the same matrix, so the counts are
        ;; (3, 1) again: only the comparison of the computed margins with the bound refuses it.
        (assert (= 3d0 (+ 3d0 (scale-float 1d0 -70))))
        (assert (not (certify-singular-values 1 3 entries 3d0 9d0 :first-margin (scale-float 1d0 -70))) ()
                "A margin lost to rounding certified a false statement."))
      ;; Dependent rows have a zero singular value: nothing positive is ever certified.
      (dolist (bottom (list (scale-float 1d0 -60) (scale-float 1d0 -30) 1d-3 1d0))
        (assert (not (try "repeated row" 2 3 '((0 0 1) (0 1 2) (1 0 1) (1 1 2)) 2d0 bottom)))
        (assert (not (try "sum of rows" 3 4 '((0 0 1) (0 1 1) (1 2 1) (1 3 1) (2 0 1) (2 1 1) (2 2 1) (2 3 1))
                          1d0 bottom)))
        (assert (not (try "more rows than columns" 3 2 '((0 0 1) (1 1 1) (2 0 1) (2 1 1)) 1d0 bottom))))
      ;; Nearly dependent: sigma_min^2 is about 2^-41; a threshold far below it is provable, one above is false.
      (let ((entries (list '(0 0 1) '(0 1 1) '(1 0 1) '(1 1 1) (list 1 2 (scale-float 1d0 -20)))))
        (assert (try "nearly dependent, far below" 2 3 entries 1d0 (scale-float 1d0 -60)))
        (assert (not (try "nearly dependent, above" 2 3 entries 1d0 (scale-float 1d0 -30)))))
      ;; Seeded sources at thresholds from negligible to larger than any singular value.
      (loop :for (m n) :in '((1 3) (2 2) (2 4) (3 3) (3 5) (4 4) (4 6) (5 7) (6 6)) :do
        (dotimes (repeat 6)
          (let* ((entries (seeded-entries state m n 55 :hub (evenp repeat)))
                 (top (coerce (max 1 (reduce #'max entries :key (lambda (e) (abs (third e))) :initial-value 0))
                              'double-float)))
            (dolist (bottom (list (scale-float 1d0 -40) 1d-6 1d-2 0.3d0 3d0 40d0))
              (try (format nil "seeded ~Dx~D #~D bottom=~E" m n repeat bottom) m n entries top (/ bottom top)))))))
    ;; The roster must show all three outcomes, or it shows nothing.
    (assert (and (> issued 150) (> withheld-false 100)) ()
            "The certificate roster issued ~D and rightly withheld ~D." issued withheld-false)
    (format t "QUASIDEFINITE-CERTIFICATE-PASS cases=~D issued=~D withheld-though-true=~D withheld-false=~D~%"
            cases issued withheld-true withheld-false)))

(defun run-quasidefinite-tests ()
  (run-quasidefinite-exact-tests)
  (run-quasidefinite-identity-tests)
  (format t "QUASIDEFINITE-BACKWARD-BOUND largest backward error over its reported bound=~,3E~%"
          (sqrt (float *qd-worst-backward* 1d0)))
  (run-quasidefinite-underflow-tests)
  (run-quasidefinite-certificate-tests)
  (run-quasidefinite-lifetime-tests)
  (run-quasidefinite-encapsulation-tests)
  ;; Warmed repeated work through a compiled Coalton caller allocates nothing: a refill, a signed
  ;; factorization, the counts, and two solves against the same factor.
  (let* ((state (sb-ext:seed-random-state 91))
         (matrix (csc-from-entries 8 12 (seeded-entries state 8 12 35 :hub t)))
         (plan (aug:analyze matrix))
         (symbolic (ldlt:analyze-quasidefinite
                    (ldlt:pattern (aug:dimension plan) (aug:column-pointers plan) (aug:row-indices plan))
                    (aug:leading-dimension plan)))
         (augmented (make-array (aug:entry-count plan) :element-type 'double-float :initial-element 0d0))
         (rhs (make-array 20 :element-type 'double-float :initial-element 1d0))
         (out (make-array 20 :element-type 'double-float :initial-element 0d0)))
    (assert (= 0 (qd-caller:step-loop plan symbolic (csc:entries matrix) 1d0 1d-6 augmented rhs out 3)) ()
            "The allocation case does not run.")
    (assert (equal '(12 8) (pivot-counts (ldlt:symbolic-factor symbolic))))
    (allocation-check "quasidefinite-step"
                      (lambda () (qd-caller:step-loop plan symbolic (csc:entries matrix) 1d0 1d-6 augmented rhs out 1))))
  (format t "QUASIDEFINITE-TESTS-PASS~%"))
