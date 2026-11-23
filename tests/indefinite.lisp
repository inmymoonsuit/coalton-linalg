(in-package #:coalton-linalg/tests)

;;; Dense symmetric indefinite LDL^T against exact rational arithmetic on the stored F64 values.
;;; The references share no code with the library: an exact reconstruction of the exported
;;; factor, an exact inertia by symmetric elimination with its own pivot rule, and an exact
;;; rational inverse.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "IL" "COALTON-LINALG/DENSE/INDEFINITE-LDLT")
  (sb-ext:add-package-local-nickname "IL-CALLER" "COALTON-LINALG/TESTS/INDEFINITE-CALLERS"))

(defparameter *unit-roundoff* (expt 2 -53))
(defvar *il-ones* 0)
(defvar *il-twos* 0)
(defvar *il-worst-residual* 0)
(defvar *il-worst-solve* 0)
(defvar *il-singular-skipped* 0)

(defun rows-matrix (rows)
  (let* ((n (length rows)) (m (make-array (list n n))))
    (dotimes (i n) (dotimes (j n) (setf (aref m i j) (nth j (nth i rows)))))
    m))

(defun build-matrix (n fn)
  (let ((m (make-array (list n n))))
    (dotimes (i n) (dotimes (j n) (setf (aref m i j) (funcall fn i j))))
    m))

(defun not-a-number ()
  (- sb-ext:double-float-positive-infinity sb-ext:double-float-positive-infinity))

(defun column-major (matrix &key junk-upper)
  "The library's input layout. With JUNK-UPPER the strict upper triangle holds NaN, which must be ignored."
  (let* ((n (array-dimension matrix 0))
         (a (make-array (* n n) :element-type 'double-float :initial-element 0d0)))
    (dotimes (j n)
      (dotimes (i n)
        (setf (aref a (+ i (* j n)))
              (if (and junk-upper (< i j)) (not-a-number) (coerce (aref matrix i j) 'double-float)))))
    a))

(defun stored-shifted (matrix shift)
  "Exact rationals of the doubles the library reads from the lower triangle, less the shift on the diagonal."
  (let ((n (array-dimension matrix 0)) (s (rational (coerce shift 'double-float))))
    (build-matrix n (lambda (i j)
                      (- (rational (coerce (aref matrix (max i j) (min i j)) 'double-float))
                         (if (= i j) s 0))))))

(defun matrix-product (a b)
  (let ((n (array-dimension a 0)))
    (build-matrix n (lambda (i j) (loop :for k :below n :sum (* (aref a i k) (aref b k j)))))))

(defun matrix-map (fn a)
  (build-matrix (array-dimension a 0) (lambda (i j) (funcall fn (aref a i j)))))

(defun matrix-transpose (a)
  (build-matrix (array-dimension a 0) (lambda (i j) (aref a j i))))

(defun matrix-max (a)
  (let ((n (array-dimension a 0)))
    (loop :for i :below n :maximize (loop :for j :below n :maximize (abs (aref a i j))) :into worst
          :finally (return (or worst 0)))))

(defun matrix-row-sum-norm (a)
  (let ((n (array-dimension a 0)))
    (loop :for i :below n :maximize (loop :for j :below n :sum (abs (aref a i j))) :into worst
          :finally (return (or worst 0)))))

(defun exact-inertia (matrix)
  "Positive, negative and zero eigenvalue counts by exact symmetric elimination (Sylvester's law)."
  (let* ((n (array-dimension matrix 0))
         (a (build-matrix n (lambda (i j) (aref matrix i j))))
         (live (loop :for i :below n :collect i)) (positive 0) (negative 0) (zero 0))
    (loop :while live :do
      (let ((i (find-if (lambda (i) (/= 0 (aref a i i))) live)))
        (cond
          (i
           (if (plusp (aref a i i)) (incf positive) (incf negative))
           (setf live (remove i live))
           (dolist (p live)
             (dolist (q live)
               (decf (aref a p q) (/ (* (aref a p i) (aref a i q)) (aref a i i))))))
          ;; Every live diagonal is zero: an off-diagonal b gives the block [[0,b],[b,0]].
          (t
           (let ((pair (loop :named found :for p :in live :do
                         (loop :for q :in live
                               :when (and (/= p q) (/= 0 (aref a p q)))
                                 :do (return-from found (cons p q))))))
             (cond
               (pair
                (let* ((i (car pair)) (j (cdr pair)) (b (aref a i j)))
                  (incf positive) (incf negative)
                  (setf live (remove i (remove j live)))
                  (dolist (p live)
                    (dolist (q live)
                      (decf (aref a p q)
                            (/ (+ (* (aref a p i) (aref a j q)) (* (aref a p j) (aref a i q))) b))))))
               (t (incf zero (length live)) (setf live nil))))))))
    (values positive negative zero)))

(defun block-inertia (factor)
  "Eigenvalue signs of the exported block diagonal, each 2-by-2 block by its exact trace and determinant."
  (let ((n (il:dimension factor)) (d (il:factor-diagonal factor)) (e (il:factor-subdiagonal factor))
        (blocks (il:factor-block-sizes factor)) (positive 0) (negative 0) (zero 0))
    (loop :with k := 0 :while (< k n) :do
      (if (= 2 (aref blocks k))
          (let* ((a (rational (aref d k))) (b (rational (aref e k))) (c (rational (aref d (1+ k))))
                 (determinant (- (* a c) (* b b))) (trace (+ a c)))
            (cond ((minusp determinant) (incf positive) (incf negative))
                  ((plusp determinant) (if (plusp trace) (incf positive 2) (incf negative 2)))
                  ((plusp trace) (incf positive) (incf zero))
                  ((minusp trace) (incf negative) (incf zero))
                  (t (incf zero 2)))
            (incf k 2))
          (let ((v (aref d k)))
            (cond ((plusp v) (incf positive)) ((minusp v) (incf negative)) (t (incf zero)))
            (incf k))))
    (values positive negative zero)))

(defun check-factor-structure (label factor)
  "Layout promised by the contract: a permutation, well-formed blocks, unit lower L, bounded multipliers."
  (let* ((n (il:dimension factor)) (perm (il:factor-permutation factor)) (l (il:factor-lower factor))
         (d (il:factor-diagonal factor)) (e (il:factor-subdiagonal factor))
         (blocks (il:factor-block-sizes factor))
         (threshold (rational (coalton:coalton il:pivot-threshold)))
         (limit (* (/ 1 (- 1 threshold)) (+ 1 (expt 2 -40)))))
    (assert (equal (sort (coerce perm 'list) #'<) (loop :for i :below n :collect i)) ()
            "~A: ~S is not a permutation." label perm)
    (loop :with k := 0 :while (< k n) :do
      (ecase (aref blocks k)
        ((1) (assert (= 0d0 (aref e k)) () "~A: subdiagonal ~S at a 1-by-1 pivot ~D." label (aref e k) k)
         (incf *il-ones*) (incf k))
        ((2) (assert (and (< (1+ k) n) (= 0 (aref blocks (1+ k))) (= 0d0 (aref e (1+ k)))) ()
                     "~A: malformed 2-by-2 block at ~D: ~S ~S." label k blocks e)
         (let ((a (abs (rational (aref d k)))) (b (abs (rational (aref e k))))
               (c (abs (rational (aref d (1+ k))))))
           (assert (and (plusp b) (< a (* threshold b)) (< c (* threshold b))) ()
                   "~A: block at ~D has diagonals ~S ~S against off-diagonal ~S." label k a c b))
         (assert (= 0d0 (aref l (+ (1+ k) (* k n)))) () "~A: L is not zero inside the block at ~D." label k)
         (incf *il-twos*) (incf k 2))))
    (dotimes (j n)
      (dotimes (i n)
        (let ((v (aref l (+ i (* j n)))))
          (cond ((= i j) (assert (= 1d0 v) () "~A: L diagonal ~S at ~D." label v i))
                ((< i j) (assert (= 0d0 v) () "~A: L upper entry ~S at ~D,~D." label v i j))
                (t (assert (<= (abs (rational v)) limit) ()
                           "~A: multiplier ~S at ~D,~D exceeds ~F." label v i j (float limit 1d0)))))))))

(defun exported-product (factor)
  "L D L^T and |L| |D| |L|^T of the exported factor as exact rational matrices, in pivoted order."
  (let* ((n (il:dimension factor)) (lower (il:factor-lower factor))
         (d (il:factor-diagonal factor)) (e (il:factor-subdiagonal factor))
         (blocks (il:factor-block-sizes factor))
         (l (build-matrix n (lambda (i j) (rational (aref lower (+ i (* j n)))))))
         (block (build-matrix n (lambda (i j)
                                  (cond ((= i j) (rational (aref d i)))
                                        ((and (= i (1+ j)) (= 2 (aref blocks j))) (rational (aref e j)))
                                        ((and (= j (1+ i)) (= 2 (aref blocks i))) (rational (aref e i)))
                                        (t 0))))))
    (values (matrix-product (matrix-product l block) (matrix-transpose l))
            (matrix-product (matrix-product (matrix-map #'abs l) (matrix-map #'abs block))
                            (matrix-transpose (matrix-map #'abs l))))))

(defun check-indefinite (label matrix shift &key (state (sb-ext:seed-random-state 7)))
  "Factor through the public interface and compare everything the contract states with exact arithmetic."
  (let* ((n (array-dimension matrix 0)) (factor (il:make-factor n))
         (exact (stored-shifted matrix shift)))
    (assert (= 0 (il-caller:observe-factor factor (column-major matrix :junk-upper t)
                                           (coerce shift 'double-float)))
            () "~A: factorization refused." label)
    (assert (il:factor-ready? factor))
    (check-factor-structure label factor)
    (let ((perm (il:factor-permutation factor)))
      (multiple-value-bind (product bound) (exported-product factor)
        (let ((worst (loop :for i :below n :maximize
                       (loop :for j :below n :maximize
                         (abs (- (aref exact (aref perm i) (aref perm j)) (aref product i j)))) :into w
                           :finally (return (or w 0))))
              (scale (* *unit-roundoff* (+ (matrix-max exact) (matrix-max bound)))))
          (assert (<= worst (* 4 (1+ n) scale)) ()
                  "~A: reconstruction error ~E exceeds ~E." label (float worst 1d0) (float (* 4 (1+ n) scale) 1d0))
          (when (plusp scale) (setf *il-worst-residual* (max *il-worst-residual* (/ worst scale)))))
        (multiple-value-bind (positive negative zero) (exact-inertia exact)
          (multiple-value-bind (block-positive block-negative block-zero) (block-inertia factor)
            (assert (equal (list block-positive block-negative block-zero)
                           (list (il:positive-count factor) (il:negative-count factor) (il:zero-count factor)))
                    () "~A: reported inertia ~S differs from the block diagonal's ~S." label
                    (list (il:positive-count factor) (il:negative-count factor) (il:zero-count factor))
                    (list block-positive block-negative block-zero)))
          (if (plusp zero)
              (incf *il-singular-skipped*)
              (let ((b (make-array n :element-type 'double-float))
                    (x (make-array n :element-type 'double-float :initial-element 0d0)))
                (assert (equal (list positive negative 0)
                               (list (il:positive-count factor) (il:negative-count factor) (il:zero-count factor)))
                        () "~A: inertia ~S differs from the exact ~S." label
                        (list (il:positive-count factor) (il:negative-count factor) (il:zero-count factor))
                        (list positive negative 0))
                (dotimes (i n) (setf (aref b i) (coerce (- (random 9 state) 4) 'double-float)))
                (assert (= 0 (il-caller:observe-solve factor b x)) () "~A: solve refused." label)
                (let* ((residual (loop :for i :below n :maximize
                                   (abs (- (loop :for j :below n :sum (* (aref exact i j) (rational (aref x j))))
                                           (rational (aref b i)))) :into w
                                       :finally (return (or w 0))))
                       (scale (* *unit-roundoff*
                                 (+ (* (matrix-row-sum-norm bound)
                                       (loop :for v :across x :maximize (abs (rational v)) :into w
                                             :finally (return (or w 0))))
                                    (loop :for v :across b :maximize (abs (rational v)) :into w
                                          :finally (return (or w 0)))))))
                  (assert (<= residual (* 8 (1+ n) scale)) ()
                          "~A: solve residual ~E exceeds ~E." label (float residual 1d0)
                          (float (* 8 (1+ n) scale) 1d0))
                  (when (plusp scale) (setf *il-worst-solve* (max *il-worst-solve* (/ residual scale))))))))))
    factor))

(defun seeded-symmetric (state n zero-diagonal)
  (let ((m (make-array (list n n) :initial-element 0)))
    (dotimes (i n)
      (dotimes (j (1+ i))
        (let ((v (if (and zero-diagonal (= i j)) 0 (- (random 7 state) 3))))
          (setf (aref m i j) v (aref m j i) v))))
    m))

(defun inertia-of (rows &optional (shift 0))
  (let* ((matrix (rows-matrix rows)) (factor (il:make-factor (length rows))))
    (assert (= 0 (il-caller:observe-factor factor (column-major matrix) (coerce shift 'double-float))))
    (values (list (il:positive-count factor) (il:negative-count factor) (il:zero-count factor)) factor)))

(defun copied-factor (factor)
  (list (copy-seq (il:factor-permutation factor)) (copy-seq (il:factor-block-sizes factor))
        (copy-seq (il:factor-lower factor)) (copy-seq (il:factor-diagonal factor))
        (copy-seq (il:factor-subdiagonal factor))))

(defun run-indefinite-exact-tests ()
  (setf *il-ones* 0 *il-twos* 0 *il-worst-residual* 0 *il-worst-solve* 0 *il-singular-skipped* 0)
  (let ((cases 0))
    (flet ((check (label rows shift)
             (incf cases)
             (check-indefinite label (if (listp rows) (rows-matrix rows) rows) shift)))
      ;; Named cases.
      (check "empty" '() 0)
      (check "one positive" '((3)) 0)
      (check "one negative" '((-3)) 0)
      (check "one zero" '((0)) 0)
      (check "exchange" '((0 1) (1 0)) 0)
      (check "positive diagonal, indefinite" '((1 2) (2 1)) 0)
      (check "definite" '((4 -2 0) (-2 5 1) (0 1 3)) 0)
      (check "negative definite" '((-4 2 0) (2 -5 -1) (0 -1 -3)) 0)
      (check "zero matrix" '((0 0 0) (0 0 0) (0 0 0)) 0)
      (check "rank one" '((1 1) (1 1)) 0)
      (check "tiny diagonals" (list (list (expt 2 -70) 1) (list 1 (expt 2 -70))) 0)
      ;; A chain of growing off-diagonals makes the search walk before it finds its pair.
      (check "walking search" (build-matrix 7 (lambda (i j) (if (= 1 (abs (- i j))) (1+ (max i j)) 0))) 0)
      (check "walking search, shifted" (build-matrix 7 (lambda (i j) (if (= 1 (abs (- i j))) (1+ (max i j)) 0))) 1/2)
      ;; Seeded integer matrices; a zero diagonal forces 2-by-2 pivots. One half is not an
      ;; algebraic integer, so the shifted matrices are never exactly singular.
      (let ((state (sb-ext:seed-random-state 20261120)))
        (loop :for n :from 1 :to 10 :do
          (dotimes (repeat 8)
            (dolist (zero-diagonal '(nil t))
              (let ((matrix (seeded-symmetric state n zero-diagonal)))
                (check (format nil "seeded n=~D #~D zero-diagonal=~A" n repeat zero-diagonal) matrix 0)
                (check (format nil "seeded n=~D #~D zero-diagonal=~A shifted" n repeat zero-diagonal)
                       matrix 1/2)))))
        ;; Congruences S diag(+-1) S^T with S unit lower triangular: the inertia is the signature's.
        (loop :for n :from 2 :to 9 :do
          (dotimes (repeat 4)
            (let* ((signs (loop :repeat n :collect (if (zerop (random 2 state)) 1 -1)))
                   (s (build-matrix n (lambda (i j) (cond ((= i j) 1) ((> i j) (- (random 5 state) 2)) (t 0)))))
                   (matrix (matrix-product
                            (matrix-product s (build-matrix n (lambda (i j) (if (= i j) (nth i signs) 0))))
                            (matrix-transpose s)))
                   (factor (check (format nil "congruence n=~D #~D" n repeat) matrix 0)))
              (assert (equal (list (count 1 signs) (count -1 signs) 0)
                             (list (il:positive-count factor) (il:negative-count factor) (il:zero-count factor)))
                      () "Congruence inertia ~S for signature ~S."
                      (list (il:positive-count factor) (il:negative-count factor)) signs))))))
    (assert (plusp *il-twos*) () "No 2-by-2 pivot was exercised.")
    (format t "INDEFINITE-EXACT-PASS cases=~D pivots-1x1=~D pivots-2x2=~D exactly-singular=~D ~
               worst-reconstruction=~,3F worst-solve=~,3F (units of roundoff times scale)~%"
            cases *il-ones* *il-twos* *il-singular-skipped*
            (float *il-worst-residual* 1d0) (float *il-worst-solve* 1d0))))

(defun run-indefinite-inertia-tests ()
  ;; The threshold is pinned here, so the structure checks do not take their bound from the library alone.
  (assert (< (abs (- (coalton:coalton il:pivot-threshold) (/ (+ 1 (sqrt 17d0)) 8))) 1d-15) ()
          "The pivot threshold is ~S." (coalton:coalton il:pivot-threshold))
  ;; Counting pivot blocks or diagonal signs is wrong; each of these would expose it.
  (assert (equal '(1 1 0) (inertia-of '((0 1) (1 0)))))
  (assert (equal '(1 1 0) (inertia-of '((1 2) (2 1)))))
  (assert (equal '(1 1 0) (inertia-of '((-1 2) (2 -1)))))
  (multiple-value-bind (inertia factor) (inertia-of '((1 2) (2 1)))
    (declare (ignore inertia))
    (assert (equalp #(2 0) (il:factor-block-sizes factor)) () "The indefinite positive-diagonal case took no 2-by-2 pivot.")
    (assert (and (plusp (aref (il:factor-diagonal factor) 0)) (plusp (aref (il:factor-diagonal factor) 1)))))
  ;; Completed singular factors: exact arithmetic throughout, so the zero count is exact.
  (assert (equal '(0 0 3) (inertia-of '((0 0 0) (0 0 0) (0 0 0)))))
  (assert (equal '(1 0 1) (inertia-of '((1 1) (1 1)))))
  (assert (equal '(1 1 1) (inertia-of '((0 1 0) (1 0 2) (0 2 0)))))
  ;; Eigenvalues of the second-difference matrix are 2 - 2 cos(k pi / 7); a shift counts those below it.
  (let* ((rows (loop :for i :below 6 :collect
                 (loop :for j :below 6 :collect (cond ((= i j) 2) ((= 1 (abs (- i j))) -1) (t 0)))))
         (spectrum (loop :for k :from 1 :to 6 :collect (- 2 (* 2 (cos (/ (* k pi) 7)))))))
    (dolist (shift '(0 1/4 1 2 3 7/2 4))
      (let ((below (count-if (lambda (v) (< v shift)) spectrum)))
        (assert (equal (list (- 6 below) below 0) (inertia-of rows shift)) ()
                "Shift ~A: inertia ~S, expected ~D below." shift (inertia-of rows shift) below))))
  ;; The augmented matrix [[I, J^T], [J, 0]] of an m-by-n J of rank r has inertia (n, r, m - r).
  (let ((state (sb-ext:seed-random-state 4242)) (full 0))
    (loop :for (m n) :in '((1 2) (2 3) (2 5) (3 4) (3 7) (4 6)) :do
      (dotimes (repeat 4)
        (let* ((j (let ((a (make-array (list m n))))
                    (dotimes (p m) (dotimes (q n) (setf (aref a p q) (- (random 5 state) 2))))
                    a))
               (size (+ n m))
               (k (build-matrix size (lambda (p q)
                                       (cond ((and (< p n) (< q n)) (if (= p q) 1 0))
                                             ((and (>= p n) (< q n)) (aref j (- p n) q))
                                             ((and (< p n) (>= q n)) (aref j (- q n) p))
                                             (t 0)))))
               (gram (build-matrix m (lambda (p q) (loop :for c :below n :sum (* (aref j p c) (aref j q c)))))))
          (multiple-value-bind (rank-positive rank-negative rank-zero) (exact-inertia gram)
            (declare (ignore rank-negative))
            (when (zerop rank-zero)
              (incf full)
              (let ((factor (check-indefinite (format nil "augmented ~Dx~D #~D" m n repeat) k 0)))
                (assert (equal (list n rank-positive 0)
                               (list (il:positive-count factor) (il:negative-count factor) (il:zero-count factor))))
                ;; The first n unknowns of K [x; y] = [0; b] are the minimum-norm solution J^T (J J^T)^-1 b.
                (let* ((b (loop :repeat m :collect (- (random 7 state) 3)))
                       (rhs (make-array size :element-type 'double-float :initial-element 0d0))
                       (out (make-array size :element-type 'double-float :initial-element 0d0))
                       (inverse (rational-inverse gram))
                       (multipliers (loop :for p :below m :collect
                                      (loop :for q :below m :sum (* (aref inverse p q) (nth q b)))))
                       (expected (loop :for c :below n :collect
                                   (loop :for p :below m :sum (* (aref j p c) (nth p multipliers))))))
                  (dotimes (p m) (setf (aref rhs (+ n p)) (coerce (nth p b) 'double-float)))
                  (assert (= 0 (il-caller:observe-solve factor rhs out)))
                  (loop :for c :below n :do
                    (assert (<= (abs (- (rational (aref out c)) (nth c expected)))
                                (* (expt 2 -40) (+ 1 (abs (nth c expected)))))
                            () "Minimum-norm component ~D is ~S, expected ~S." c (aref out c)
                            (float (nth c expected) 1d0))))))))))
    (assert (>= full 12) () "Only ~D augmented systems had full row rank." full)
    (format t "INDEFINITE-AUGMENTED-PASS full-row-rank=~D of 24~%" full))
  ;; Dependent rows: the augmented matrix is singular, the factor completes, the solve says so.
  (multiple-value-bind (inertia factor)
      (inertia-of '((1 0 0 1 1) (0 1 0 0 0) (0 0 1 0 0) (1 0 0 0 0) (1 0 0 0 0)))
    (assert (equal '(3 1 1) inertia) () "Dependent-row augmented inertia ~S." inertia)
    (let ((out (make-array 5 :element-type 'double-float :initial-element 9d0)))
      (assert (= 3 (il-caller:observe-solve factor (doubles 0d0 0d0 0d0 1d0 1d0) out)))
      (assert (every (lambda (v) (= v 9d0)) out) () "A singular solve wrote its output.")))
  (format t "INDEFINITE-INERTIA-PASS~%"))

(defun run-indefinite-identity-tests ()
  ;; Scaling by a power of two leaves the pivots, the permutation and L unchanged and scales D
  ;; exactly, including where an off-diagonal squared would overflow or underflow.
  (let* ((matrix (rows-matrix '((0 3 1 -2) (3 0 2 1) (1 2 1 0) (-2 1 0 -1))))
         (base (il:make-factor 4)) (rhs (doubles 1d0 -2d0 3d0 1d0))
         (base-x (make-array 4 :element-type 'double-float :initial-element 0d0)))
    (assert (= 0 (il-caller:observe-factor base (column-major matrix) 0d0)))
    (assert (find 2 (il:factor-block-sizes base)) () "The scaling case has no 2-by-2 pivot.")
    (assert (= 0 (il-caller:observe-solve base rhs base-x)))
    (dolist (power '(-600 -31 29 600))
      (let ((scaled (il:make-factor 4)) (x (make-array 4 :element-type 'double-float :initial-element 0d0)))
        (assert (= 0 (il-caller:observe-factor
                      scaled (column-major (matrix-map (lambda (v) (* v (expt 2 power))) matrix)) 0d0)))
        (assert (equalp (il:factor-permutation base) (il:factor-permutation scaled)))
        (assert (equalp (il:factor-block-sizes base) (il:factor-block-sizes scaled)))
        (assert (equalp (il:factor-lower base) (il:factor-lower scaled)) () "Scaling by 2^~D changed L." power)
        (dotimes (i 4)
          (assert (= (aref (il:factor-diagonal scaled) i) (scale-float (aref (il:factor-diagonal base) i) power)))
          (assert (= (aref (il:factor-subdiagonal scaled) i)
                     (scale-float (aref (il:factor-subdiagonal base) i) power))))
        (assert (= 0 (il-caller:observe-solve scaled rhs x)))
        (dotimes (i 4)
          (assert (= (aref x i) (scale-float (aref base-x i) (- power))) ()
                  "Scaling by 2^~D gave x[~D] = ~S." power i (aref x i)))))
    ;; A shift applied by factor! is the matrix stored with the shift already subtracted.
    (let ((shifted (il:make-factor 4)) (stored (il:make-factor 4)))
      (assert (= 0 (il-caller:observe-factor shifted (column-major matrix) 2d0)))
      (assert (= 0 (il-caller:observe-factor
                    stored (column-major (build-matrix 4 (lambda (i j) (- (aref matrix i j) (if (= i j) 2 0)))))
                    0d0)))
      (assert (equalp (copied-factor shifted) (copied-factor stored)) () "A shift is not a stored shift."))
    ;; The strict upper triangle is not read, and a repeated call repeats every bit.
    (let ((clean (il:make-factor 4)) (junk (il:make-factor 4)))
      (assert (= 0 (il-caller:observe-factor clean (column-major matrix) 0d0)))
      (assert (= 0 (il-caller:observe-factor junk (column-major matrix :junk-upper t) 0d0)))
      (assert (equalp (copied-factor clean) (copied-factor junk)) () "The upper triangle was read.")
      (let ((first (copied-factor clean)))
        (assert (= 0 (il-caller:observe-factor clean (column-major matrix) 0d0)))
        (assert (equalp first (copied-factor clean)) () "A repeated factorization changed."))))
  (format t "INDEFINITE-IDENTITIES-PASS~%"))

(defun run-indefinite-lifetime-tests ()
  (let* ((matrix (rows-matrix '((4 -2 0) (-2 5 1) (0 1 3))))
         (a (column-major matrix)) (factor (il:make-factor 3))
         (rhs (doubles 2d0 4d0 4d0)) (out (make-array 3 :element-type 'double-float :initial-element 9d0))
         (infinity sb-ext:double-float-positive-infinity))
    (flet ((untouched () (every (lambda (v) (= v 9d0)) out))
           (cleared ()
             (and (not (il:factor-ready? factor)) (= 0 (il:positive-count factor))
                  (= 0 (il:negative-count factor)) (= 0 (il:zero-count factor)))))
      ;; Nothing stands before the first factorization.
      (assert (cleared))
      (assert (= 1 (il-caller:observe-solve factor rhs out)))
      (assert (untouched))
      (assert (= 0 (il-caller:observe-factor factor a 0d0)))
      (assert (equal '(3 0 0) (list (il:positive-count factor) (il:negative-count factor) (il:zero-count factor))))
      ;; Solve refusals leave the output alone; the factor stays ready.
      (assert (= 2 (il-caller:observe-solve factor (doubles 1d0 2d0) out)))
      (assert (= 2 (il-caller:observe-solve factor rhs (make-array 2 :element-type 'double-float))))
      (assert (= 4 (il-caller:observe-solve factor (doubles 1d0 infinity 3d0) out)))
      (assert (= 4 (il-caller:observe-solve factor (doubles 1d0 (not-a-number) 3d0) out)))
      (assert (untouched))
      (assert (il:factor-ready? factor))
      (assert (= 0 (il-caller:observe-solve factor rhs out)))
      (assert (close-vector-p out (doubles 1d0 1d0 1d0)))
      ;; The right-hand side and the output may be one array.
      (let ((shared (doubles 2d0 4d0 4d0)))
        (assert (= 0 (il-caller:observe-solve factor shared shared)))
        (assert (close-vector-p shared (doubles 1d0 1d0 1d0))))
      ;; Every refused factorization clears the factor, whatever stood before.
      (fill out 9d0)
      (dolist (attempt (list (lambda () (il-caller:observe-factor factor (doubles 1d0 2d0 3d0) 0d0))
                             (lambda () (il-caller:observe-factor factor a infinity))
                             (lambda () (il-caller:observe-factor factor a (not-a-number)))
                             (lambda ()
                               (let ((bad (copy-seq a)))
                                 (setf (aref bad 5) infinity)
                                 (il-caller:observe-factor factor bad 0d0)))
                             (lambda ()
                               (let ((bad (copy-seq a)))
                                 (setf (aref bad 0) (not-a-number))
                                 (il-caller:observe-factor factor bad 0d0)))))
        (assert (= 0 (il-caller:observe-factor factor a 0d0)))
        (assert (/= 0 (funcall attempt)))
        (assert (cleared) () "A refused factorization left a factor standing.")
        (assert (= 1 (il-caller:observe-solve factor rhs out)))
        (assert (untouched)))
      (assert (= 2 (il-caller:observe-factor factor (doubles 1d0 2d0 3d0) 0d0)))
      (assert (= 3 (il-caller:observe-factor factor a infinity)))
      (assert (= 1 (il-caller:observe-factor factor (let ((bad (copy-seq a))) (setf (aref bad 5) infinity) bad) 0d0)))))
  ;; Overflow while eliminating is a refusal, not a factor holding infinities.
  (let* ((big (/ most-positive-double-float 1.2d0))
         (factor (il:make-factor 2)))
    (assert (= 1 (il-caller:observe-factor factor (doubles big big big (- big)) 0d0)))
    (assert (not (il:factor-ready? factor)))
    ;; The first pivot was counted before the overflow; a refusal reports no inertia.
    (assert (equal '(0 0 0) (list (il:positive-count factor) (il:negative-count factor) (il:zero-count factor)))
            () "A refused factorization kept a partial inertia."))
  ;; A shift that overflows a diagonal is refused as non-finite arithmetic.
  (let ((factor (il:make-factor 1)))
    (assert (= 1 (il-caller:observe-factor factor (doubles most-positive-double-float) (- most-positive-double-float))))
    (assert (not (il:factor-ready? factor))))
  ;; A finite factor whose solution overflows.
  (let ((factor (il:make-factor 1)) (out (doubles 9d0)))
    (assert (= 0 (il-caller:observe-factor factor (doubles (scale-float 1d0 -1000)) 0d0)))
    (assert (= 5 (il-caller:observe-solve factor (doubles (scale-float 1d0 100)) out)))
    (assert (= 9d0 (aref out 0))))
  ;; The dimension bound is checked before any allocation that depends on it.
  (assert (= 0 (il-caller:allocation-error 3)))
  (assert (= 3 (il-caller:allocation-error (1+ (coalton:coalton il:maximum-dimension)))))
  (assert (= 3 (il-caller:allocation-error most-positive-fixnum)))
  ;; Within the dimension bound, storage is refused by an entry budget before it is allocated.
  (assert (= 16777216 (coalton:coalton il:default-entry-limit)))
  (assert (= 4 (il-caller:allocation-error 4097)) () "The default budget admitted 4097 * 4097 entries.")
  (assert (= 4 (il-caller:allocation-error (coalton:coalton il:maximum-dimension))))
  (assert (= 0 (il-caller:budget-error 3 9)))
  (assert (= 9 (il-caller:budget-error 3 8)) () "A budget of 8 admitted a 3-by-3 factor.")
  (assert (= 0 (il-caller:budget-error 0 0)))
  (assert (= 1 (il-caller:budget-error 1 0)))
  (format t "INDEFINITE-LIFETIME-PASS~%"))

(defun check-large-two-pivot ()
  (let* ((h most-positive-double-float) (a (/ h 4)) (f (il:make-factor 2)) (x (doubles 0d0 0d0)))
    (assert (= 0 (il-caller:observe-factor f (doubles a h 0d0 (- a)) 0d0)))
    (assert (= 0 (il-caller:observe-solve f (doubles (/ h 2) 0d0) x)))
    (assert (close-vector-p x (doubles (float 2/17 1d0) (float 8/17 1d0)))))
  ;; Dividing the RHS by the off-diagonal first can overflow even when its
  ;; product with the inverse block is representable.
  (let* ((h (/ most-positive-double-float 8)) (rhs (* 7 h)) (f (il:make-factor 2)) (x (doubles 0d0 0d0)))
    (assert (= 0 (il-caller:observe-factor f (doubles 0.375d0 0.75d0 0d0 -0.375d0) 0d0)))
    (assert (= 0 (il-caller:observe-solve f (doubles rhs 0d0) x)))
    (assert (close-vector-p x (doubles (float (* 8/15 (rational rhs)) 1d0)
                                      (float (* 16/15 (rational rhs)) 1d0))))
    (assert (= 0 (il-caller:observe-factor f (doubles 0.625d0 1.25d0 0d0 -0.625d0) 0d0)))
    (assert (= 0 (il-caller:observe-solve f (doubles rhs rhs) x)))
    (assert (close-vector-p x (doubles (float (* 24/25 (rational rhs)) 1d0)
                                      (float (* 8/25 (rational rhs)) 1d0)))))
  ;; The third row exercises the same block inverse while constructing L.
  (let* ((h most-positive-double-float) (a (/ h 4)) (f (il:make-factor 3)))
    (assert (= 0 (il-caller:observe-factor f (doubles a h (/ h 2) 0d0 (- a) 0d0 0d0 0d0 a) 0d0)))
    (assert (close-vector-p (doubles (aref (il:factor-lower f) 2) (aref (il:factor-lower f) 5))
                            (doubles (float 2/17 1d0) (float 8/17 1d0))))))

(defun run-indefinite-tests ()
  (check-large-two-pivot)
  (run-indefinite-exact-tests)
  (run-indefinite-inertia-tests)
  (run-indefinite-identity-tests)
  (run-indefinite-lifetime-tests)
  (run-indefinite-encapsulation-tests)
  ;; Warmed repeated work through a compiled Coalton caller allocates nothing.
  (let* ((state (sb-ext:seed-random-state 99))
         (matrix (seeded-symmetric state 12 t))
         (a (column-major matrix)) (factor (il:make-factor 12))
         (rhs (make-array 12 :element-type 'double-float :initial-element 1d0))
         (out (make-array 12 :element-type 'double-float :initial-element 0d0)))
    (assert (= 0 (il-caller:factor-solve-loop factor a rhs out 3)) () "The allocation case does not factor and solve.")
    (assert (find 2 (il:factor-block-sizes factor)) () "The allocation case has no 2-by-2 pivot.")
    (allocation-check "indefinite-factor-solve"
                      (lambda () (il-caller:factor-solve-loop factor a rhs out 1))))
  (format t "INDEFINITE-TESTS-PASS~%"))

(defparameter *indefinite-public-control*
  "(declare inspect (coalton-linalg/dense/indefinite-ldlt:Factor -> UFix))
   (define (inspect factor)
     (if (coalton-linalg/dense/indefinite-ldlt:factor-ready? factor)
         1
         (coalton-linalg/dense/indefinite-ldlt:dimension factor)))
   (declare probe (Void -> UFix))
   (define (probe)
     (inspect (coalton-linalg/dense/indefinite-ldlt:make-factor 3)))")

(defparameter *indefinite-encapsulation-probes*
  '(("dense-factor-construction" .
     "(declare probe (Void -> coalton-linalg/dense/indefinite-ldlt:Factor))
      (define (probe)
        (coalton-linalg/dense/indefinite-ldlt:Factor
          0 (arr:make 0 0.0d0) (arr:make 0 0.0d0) (arr:make 0 0.0d0)
          (arr:make 0 (the UFix 0)) (arr:make 0 (the UFix 0)) (arr:make 0 0.0d0)
          (arr:make 3 (the UFix 0)) (arr:make 4 (the UFix 0))))")
    ("dense-factor-destructuring" .
     "(declare probe (coalton-linalg/dense/indefinite-ldlt:Factor -> UFix))
      (define (probe factor)
        (match factor
          ((coalton-linalg/dense/indefinite-ldlt:Factor _ _ _ _ _ _ _ _ state) (arr:length state))))")
    ("dense-factor-state-access" .
     "(declare probe (coalton-linalg/dense/indefinite-ldlt:Factor -> UFix))
      (define (probe factor) (arr:aref (.state factor) 0))")
    ("dense-factor-scratch-access" .
     "(declare probe (coalton-linalg/dense/indefinite-ldlt:Factor -> arr:LispArray F64))
      (define (probe factor) (.scratch factor))")
    ("dense-status-confusion" .
     "(declare probe (ldlt:FactorStatus -> coalton-linalg/dense/indefinite-ldlt:FactorStatus))
      (define (probe status) status)")))

(defun run-indefinite-encapsulation-tests ()
  "External Coalton callers can use the dense factor only through its public operations."
  (multiple-value-bind (accepted diagnostics)
      (coalton-linalg/tests/encapsulation::compile-probe "dense-public-control" *indefinite-public-control*
                                                         :execute t)
    (assert accepted () "The valid dense external caller failed:~%~A" diagnostics))
  (dolist (probe *indefinite-encapsulation-probes*)
    (assert (not (coalton-linalg/tests/encapsulation::compile-probe (car probe) (cdr probe))) ()
            "External caller ~A compiled despite required rejection." (car probe)))
  (dolist (name '("MKFACTOR" "FACTORSTORAGE"))
    (multiple-value-bind (symbol visibility) (find-symbol name "COALTON-LINALG/DENSE/INDEFINITE-LDLT")
      (assert (and symbol (eq visibility ':internal)) ()
              "Storage constructor ~A must exist and remain internal; visibility is ~S." name visibility)))
  (format t "INDEFINITE-ENCAPSULATION-PASS rejected=~D public-control=1~%"
          (length *indefinite-encapsulation-probes*)))
