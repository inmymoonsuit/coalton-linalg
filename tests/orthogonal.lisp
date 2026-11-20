(in-package #:coalton-linalg/tests)

;;; Pivoted Householder QR and the complete orthogonal decomposition against exact rational
;;; arithmetic on the stored F64 values. The references share no code with the library: Q is
;;; read back column by column and multiplied out exactly, and the minimum-norm least-squares
;;; solution comes from an exact rank factorization.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "OG" "COALTON-LINALG/DENSE/ORTHOGONAL")
  (sb-ext:add-package-local-nickname "OG-CALLER" "COALTON-LINALG/TESTS/ORTHOGONAL-CALLERS"))

(defvar *og-worst-orthogonality* 0)
(defvar *og-worst-reconstruction* 0)
(defvar *og-worst-solution* 0)
(defvar *og-rank-deficient* 0)

(defun rect (m n fn)
  (let ((a (make-array (list m n))))
    (dotimes (i m) (dotimes (j n) (setf (aref a i j) (funcall fn i j))))
    a))

(defun rect-of-rows (rows)
  (rect (length rows) (length (first rows)) (lambda (i j) (nth j (nth i rows)))))

(defun rect-product (a b)
  (let ((inner (array-dimension a 1)))
    (rect (array-dimension a 0) (array-dimension b 1)
          (lambda (i j) (loop for k below inner sum (* (aref a i k) (aref b k j)))))))

(defun rect-transpose (a)
  (rect (array-dimension a 1) (array-dimension a 0) (lambda (i j) (aref a j i))))

(defun rect-max (a)
  (let ((worst 0))
    (dotimes (i (array-dimension a 0))
      (dotimes (j (array-dimension a 1))
        (setf worst (max worst (abs (aref a i j))))))
    worst))

(defun rect-column-major (a)
  (let* ((m (array-dimension a 0)) (n (array-dimension a 1))
         (out (make-array (* m n) :element-type 'double-float :initial-element 0d0)))
    (dotimes (j n) (dotimes (i m) (setf (aref out (+ i (* j m))) (coerce (aref a i j) 'double-float))))
    out))

(defun rect-stored (a)
  "Exact rationals of the doubles the library reads."
  (rect (array-dimension a 0) (array-dimension a 1)
        (lambda (i j) (rational (coerce (aref a i j) 'double-float)))))

(defun zero-vector (n) (make-array n :element-type 'double-float :initial-element 0d0))

(defun exact-basis-columns (a)
  "Indices of a maximal independent set of columns, by exact elimination."
  (let* ((m (array-dimension a 0)) (n (array-dimension a 1))
         (work (rect m n (lambda (i j) (aref a i j)))) (row 0) (basis '()))
    (dotimes (j n)
      (when (< row m)
        (let ((p (loop for i from row below m unless (zerop (aref work i j)) return i)))
          (when p
            (dotimes (c n) (rotatef (aref work p c) (aref work row c)))
            (loop for i from (1+ row) below m do
              (let ((factor (/ (aref work i j) (aref work row j))))
                (unless (zerop factor)
                  (loop for c from j below n do (decf (aref work i c) (* factor (aref work row c)))))))
            (push j basis)
            (incf row)))))
    (nreverse basis)))

(defun exact-minimum-norm (a b)
  "The minimum-norm least-squares solution, exactly, as C^T (C C^T)^-1 (B^T B)^-1 B^T b with A = B C."
  (let* ((m (array-dimension a 0)) (n (array-dimension a 1))
         (basis (exact-basis-columns a)) (r (length basis)))
    (if (zerop r)
        (make-list n :initial-element 0)
        (let* ((bm (rect m r (lambda (i k) (aref a i (nth k basis)))))
               (left (rect-product (rational-inverse (rect-product (rect-transpose bm) bm)) (rect-transpose bm)))
               (c (rect-product left a))
               (coefficients (rect-product left (rect m 1 (lambda (i j) (declare (ignore j)) (nth i b)))))
               (x (rect-product (rect-product (rect-transpose c)
                                              (rational-inverse (rect-product c (rect-transpose c))))
                                coefficients)))
          (loop for i below n collect (aref x i 0))))))

(defun og-q-matrix (decomposition)
  "Q as the exact rationals of the floats apply-q! returns for each unit vector."
  (let* ((m (og:row-count decomposition)) (q (make-array (list m m))))
    (dotimes (j m)
      (let ((v (zero-vector m)))
        (setf (aref v j) 1d0)
        (assert (= 0 (og-caller:observe-apply-q decomposition v)))
        (dotimes (i m) (setf (aref q i j) (rational (aref v i))))))
    q))

(defun og-r-matrix (decomposition)
  (let ((m (og:row-count decomposition)) (work (og:factor-work decomposition)))
    (rect m (og:column-count decomposition)
          (lambda (i j) (if (<= i j) (rational (aref work (+ i (* j m)))) 0)))))

(defun largest-column-norm (exact)
  "An upper bound, as a rational, on the largest Euclidean column norm."
  (let ((worst 0))
    (dotimes (j (array-dimension exact 1))
      (setf worst (max worst (loop for i below (array-dimension exact 0) sum (expt (aref exact i j) 2)))))
    (* (rational (sqrt (coerce worst 'double-float))) (+ 1 (expt 2 -40)))))

(defun check-qr-accuracy (label decomposition exact)
  "Orthogonality of Q, A P = Q R, Q^T against Q, and the ordering the pivoting promises."
  (let* ((m (og:row-count decomposition)) (n (og:column-count decomposition))
         (steps (og:step-count decomposition)) (perm (og:factor-permutation decomposition))
         (q (og-q-matrix decomposition)) (r (og-r-matrix decomposition))
         (scale (* *unit-roundoff* (largest-column-norm exact))))
    (assert (= steps (min m n)))
    (assert (equal (sort (coerce perm 'list) #'<) (loop for j below n collect j)) ()
            "~A: ~S is not a permutation." label perm)
    (let* ((gram (rect-product (rect-transpose q) q))
           (defect (rect-max (rect m m (lambda (i j) (- (aref gram i j) (if (= i j) 1 0)))))))
      (assert (<= defect (* 8 (1+ m) *unit-roundoff*)) ()
              "~A: Q^T Q differs from the identity by ~E." label (float defect 1d0))
      (setf *og-worst-orthogonality* (max *og-worst-orthogonality* (/ defect *unit-roundoff*))))
    (let* ((product (rect-product q r))
           (error (rect-max (rect m n (lambda (i j) (- (aref product i j) (aref exact i (aref perm j))))))))
      (assert (<= error (* 8 (+ m n 1) scale)) ()
              "~A: A P differs from Q R by ~E, allowed ~E." label (float error 1d0)
              (float (* 8 (+ m n 1) scale) 1d0))
      (when (plusp scale) (setf *og-worst-reconstruction* (max *og-worst-reconstruction* (/ error scale)))))
    ;; apply-qt! is the transpose of what apply-q! applies.
    (dotimes (i m)
      (let ((v (zero-vector m)))
        (setf (aref v i) 1d0)
        (assert (= 0 (og-caller:observe-apply-qt decomposition v)))
        (dotimes (j m)
          (assert (<= (abs (- (rational (aref v j)) (aref q i j))) (* 8 (1+ m) *unit-roundoff*)) ()
                  "~A: Q^T e_~D entry ~D is ~S, Q has ~S." label i j (aref v j) (float (aref q i j) 1d0)))))
    ;; The diagonal of R does not grow, and each entry dominates the rest of its row block.
    (dotimes (k steps)
      (assert (= (abs (aref r k k)) (rational (og:pivot-magnitude decomposition k))))
      (when (< (1+ k) steps)
        (assert (>= (* (abs (aref r k k)) (+ 1 (expt 2 -40))) (abs (aref r (1+ k) (1+ k)))) ()
                "~A: |r| grows at ~D." label k))
      (loop for j from (1+ k) below n do
        (assert (>= (* (expt (aref r k k) 2) (+ 1 (expt 2 -40)))
                    (loop for i from k to (min j (1- steps)) sum (expt (aref r i j) 2)))
                () "~A: column ~D outweighs the pivot at ~D." label j k)))
    (assert (= 0d0 (og:pivot-magnitude decomposition steps)))))

(defun check-orthogonal (label matrix b &key (relative 1d-9))
  "Factor, pick the rank by a relative threshold, and compare the solution and nullspace with exact arithmetic."
  (let* ((m (array-dimension matrix 0)) (n (array-dimension matrix 1))
         (decomposition (og:make-decomposition m n)) (exact (rect-stored matrix))
         (exact-b (mapcar (lambda (v) (rational (coerce v 'double-float))) b))
         (exact-rank (length (exact-basis-columns exact))))
    (assert (= 0 (og-caller:observe-factor decomposition (rect-column-major matrix))) ()
            "~A: factorization refused." label)
    (check-qr-accuracy label decomposition exact)
    (let ((rank (og:rank-above decomposition (* relative (og:pivot-magnitude decomposition 0))))
          (x (zero-vector n))
          (rhs (make-array m :element-type 'double-float
                             :initial-contents (mapcar (lambda (v) (coerce v 'double-float)) b))))
      (assert (= rank exact-rank) () "~A: rank ~D by threshold, exact rank ~D." label rank exact-rank)
      (when (< rank (min m n)) (incf *og-rank-deficient*))
      (assert (= 0 (og-caller:observe-rank decomposition rank)))
      (assert (and (og:rank-set? decomposition) (= rank (og:selected-rank decomposition))))
      (assert (= 0 (og-caller:observe-solve decomposition rhs x)) () "~A: solve refused." label)
      (let* ((expected (exact-minimum-norm exact exact-b))
             (size (+ 1 (loop for v in expected maximize (abs v) into w finally (return (or w 0)))))
             (error (loop for i below n maximize (abs (- (rational (aref x i)) (nth i expected))) into w
                          finally (return (or w 0))))
             (residual-squared
               (loop for i below m sum
                 (expt (- (nth i exact-b) (loop for j below n sum (* (aref exact i j) (nth j expected)))) 2)))
             (b-size (+ 1 (loop for v in exact-b maximize (abs v) into w finally (return (or w 0))))))
        (assert (<= error (* (expt 2 -36) size)) ()
                "~A: solution error ~E against the exact minimum-norm solution of size ~E." label
                (float error 1d0) (float size 1d0))
        (setf *og-worst-solution* (max *og-worst-solution* (/ error size *unit-roundoff*)))
        (assert (<= (abs (- (og:solve-residual-norm decomposition) (sqrt (coerce residual-squared 'double-float))))
                    (* 1d-10 (float b-size 1d0)))
                () "~A: residual norm ~S, exact ~S." label (og:solve-residual-norm decomposition)
                (sqrt (coerce residual-squared 'double-float))))
      ;; The selected nullspace: orthonormal, annihilated by A, and orthogonal to the solution.
      (let ((basis (loop for index below (- n rank) collect
                     (let ((z (zero-vector n)))
                       (assert (= 0 (og-caller:observe-nullspace decomposition index z)))
                       (map 'list #'rational z))))
            (tolerance (expt 2 -40)))
        (assert (= 4 (og-caller:observe-nullspace decomposition (- n rank) (zero-vector n))))
        (loop for z in basis for index from 0 do
          (dotimes (i m)
            (assert (<= (abs (loop for j below n sum (* (aref exact i j) (nth j z))))
                        (* tolerance (+ 1 (rect-max exact))))
                    () "~A: A times nullspace vector ~D is not small in row ~D." label index i))
          (assert (<= (abs (loop for j below n sum (* (rational (aref x j)) (nth j z))))
                      (* tolerance (+ 1 (loop for v across x maximize (abs (rational v))))))
                  () "~A: the solution is not orthogonal to nullspace vector ~D." label index)
          (loop for other in basis for other-index from 0 do
            (assert (<= (abs (- (loop for j below n sum (* (nth j z) (nth j other)))
                                (if (= index other-index) 1 0)))
                        tolerance)
                    () "~A: nullspace vectors ~D and ~D are not orthonormal." label index other-index)))))
    decomposition))

(defun seeded-product (state m n rank)
  "An m-by-n integer matrix of rank at most RANK: a product of two integer factors, or zero."
  (if (zerop rank)
      (rect m n (lambda (i j) (declare (ignore i j)) 0))
      (rect-product (rect m rank (lambda (i j) (declare (ignore i j)) (- (random 5 state) 2)))
                    (rect rank n (lambda (i j) (declare (ignore i j)) (- (random 5 state) 2))))))

(defun run-orthogonal-exact-tests ()
  (setf *og-worst-orthogonality* 0 *og-worst-reconstruction* 0 *og-worst-solution* 0 *og-rank-deficient* 0)
  (let ((cases 0))
    (flet ((check (label matrix b)
             (incf cases)
             (check-orthogonal label (if (listp matrix) (rect-of-rows matrix) matrix) b)))
      (check "line fit" '((1 0) (1 1) (1 2)) '(1 3 5))
      (check "line fit, inconsistent" '((1 0) (1 1) (1 2)) '(1 3 4))
      (check "one constraint" '((1 1)) '(2))
      (check "rank one, consistent" '((1 2) (2 4)) '(1 2))
      (check "rank one, inconsistent" '((1 2) (2 4)) '(1 0))
      (check "zero matrix" '((0 0 0) (0 0 0)) '(3 4))
      (check "one by one" '((5)) '(10))
      (check "one by one zero" '((0)) '(7))
      (check "no rows" (rect 0 3 (lambda (i j) (declare (ignore i j)) 0)) '())
      (check "no columns" (rect 3 0 (lambda (i j) (declare (ignore i j)) 0)) '(1 2 3))
      (check "nothing" (rect 0 0 (lambda (i j) (declare (ignore i j)) 0)) '())
      (check "zero column between" '((1 0 2) (3 0 4) (5 0 7)) '(1 1 1))
      (check "repeated columns" '((1 1 2) (2 2 1) (3 3 0) (4 4 5)) '(1 0 2 1))
      ;; The first pivot column needs no reflection; the columns after it must still be measured
      ;; over the rows that remain.
      (check "first column already reduced" '((4 3 0) (0 1 2) (0 0 1)) '(1 2 3))
      ;; Seeded integer matrices of every rank from zero to full, tall, square and wide.
      (let ((state (sb-ext:seed-random-state 20261121)))
        (loop for (m n) in '((1 1) (2 1) (1 3) (2 2) (3 2) (2 4) (3 3) (5 3) (3 6) (4 4) (6 4) (4 7)
                             (5 5) (7 5) (5 8) (6 6))
              do (loop for rank from 0 to (min m n) do
                   (dotimes (repeat 3)
                     (check (format nil "seeded ~Dx~D rank<=~D #~D" m n rank repeat)
                            (seeded-product state m n rank)
                            (loop repeat m collect (- (random 9 state) 4))))))))
    (assert (plusp *og-rank-deficient*) () "No rank-deficient case was exercised.")
    (format t "ORTHOGONAL-EXACT-PASS cases=~D rank-deficient=~D worst-orthogonality=~,3F ~
               worst-reconstruction=~,3F worst-solution=~,1F (units of roundoff times scale)~%"
            cases *og-rank-deficient* (float *og-worst-orthogonality* 1d0)
            (float *og-worst-reconstruction* 1d0) (float *og-worst-solution* 1d0))))

(defun run-orthogonal-rank-tests ()
  ;; Thresholds are absolute and strict: an entry equal to the threshold is not above it.
  (let ((decomposition (og:make-decomposition 3 3))
        (x (zero-vector 3)))
    (assert (= 0 (og:rank-above decomposition 0d0)) () "A rank was reported before any factorization.")
    (assert (= 0 (og-caller:observe-factor decomposition (doubles 1d0 0d0 0d0  0d0 4d0 0d0  0d0 0d0 2d0))))
    (assert (equalp (list 4d0 2d0 1d0) (loop for k below 3 collect (og:pivot-magnitude decomposition k))))
    (assert (equalp #(1 2 0) (og:factor-permutation decomposition)))
    (loop for (threshold expected) in '((-1d0 3) (0d0 3) (0.999d0 3) (1d0 2) (1.5d0 2) (2d0 1) (4d0 0) (9d0 0))
          do (assert (= expected (og:rank-above decomposition threshold)) ()
                     "Threshold ~S gave rank ~D." threshold (og:rank-above decomposition threshold)))
    ;; A chosen rank below the exact one solves the truncated problem: the smallest pivot's
    ;; column is left out and its right-hand side entry becomes residual.
    (assert (= 0 (og-caller:observe-rank decomposition 2)))
    (assert (= 0 (og-caller:observe-solve decomposition (doubles 3d0 8d0 6d0) x)))
    (assert (equalp x (doubles 0d0 2d0 3d0)) () "Truncated solution ~S." x)
    (assert (= 3d0 (og:solve-residual-norm decomposition)))
    (assert (= 0 (og-caller:observe-rank decomposition 3)))
    (assert (= 0 (og-caller:observe-solve decomposition (doubles 3d0 8d0 6d0) x)))
    (assert (equalp x (doubles 3d0 2d0 3d0)))
    (assert (= 0d0 (og:solve-residual-norm decomposition)))
    (assert (= 2 (og-caller:observe-rank decomposition 4)))
    (assert (not (og:rank-set? decomposition)) () "A refused rank left the earlier one standing."))
  ;; A rank that reaches an exactly zero pivot is refused rather than divided by.
  (let ((decomposition (og:make-decomposition 2 3)))
    (assert (= 0 (og-caller:observe-factor decomposition (doubles 1d0 2d0  2d0 4d0  0d0 0d0))))
    (assert (= 0 (og-caller:observe-rank decomposition 1)))
    (assert (= 0 (og-caller:observe-factor decomposition (doubles 0d0 0d0  0d0 0d0  0d0 0d0))))
    (assert (= 3 (og-caller:observe-rank decomposition 1)))
    (assert (= 0 (og-caller:observe-rank decomposition 0))))
  ;; The Kahan matrix: its remaining columns tie in norm at every step, so column pivoting
  ;; leaves it alone, and the last diagonal entry of R is far above the smallest singular
  ;; value. A diagonal perturbation of a few rounding units, decreasing along the diagonal,
  ;; keeps rounding from breaking the ties the other way. The bound on the smallest singular
  ;; value comes from the exact inverse of the stored matrix.
  (let* ((n 30) (c (cos 1.2d0)) (s (sin 1.2d0))
         (matrix (rect n n (lambda (i j)
                             (cond ((= i j) (+ (expt s i) (* 25 (- n i) double-float-epsilon)))
                                   ((< i j) (* (- c) (expt s i)))
                                   (t 0d0)))))
         (exact (rect-stored matrix)) (decomposition (og:make-decomposition n n)))
    (assert (= 0 (og-caller:observe-factor decomposition (rect-column-major matrix))))
    (check-qr-accuracy "Kahan" decomposition exact)
    (let* ((inverse (rational-inverse exact))
           (largest (loop for j below n maximize (loop for i below n sum (expt (aref inverse i j) 2))))
           (singular-bound (/ 1 (sqrt (coerce largest 'double-float))))
           (last-pivot (og:pivot-magnitude decomposition (1- n)))
           (ratio (/ last-pivot singular-bound)))
      (assert (equalp (og:factor-permutation decomposition) (coerce (loop for j below n collect j) 'vector)) ()
              "Kahan: columns were exchanged: ~S." (og:factor-permutation decomposition))
      (assert (= n (og:rank-above decomposition (* 0.5d0 last-pivot))))
      (assert (> ratio 1000d0) () "Kahan: last pivot ~S over the singular value bound ~S is only ~S."
              last-pivot singular-bound ratio)
      (format t "ORTHOGONAL-KAHAN n=~D last-pivot=~,3E smallest-singular-value<=~,3E ratio>=~,1F~%"
              n last-pivot singular-bound ratio)))
  (format t "ORTHOGONAL-RANK-PASS~%"))

(defun copied-decomposition (decomposition)
  (list (copy-seq (og:factor-permutation decomposition)) (copy-seq (og:factor-work decomposition))
        (copy-seq (og:factor-tau decomposition))))

(defun run-orthogonal-identity-tests ()
  ;; Scaling by a power of two scales R exactly and leaves the reflectors and the pivoting
  ;; unchanged, including where an unscaled sum of squares would overflow or underflow.
  (let* ((matrix (rect-of-rows '((3 1 -2) (1 0 2) (-2 2 1) (1 1 1))))
         (base (og:make-decomposition 4 3)) (rhs (doubles 1d0 -2d0 3d0 1d0)) (base-x (zero-vector 3)))
    (assert (= 0 (og-caller:observe-factor base (rect-column-major matrix))))
    (assert (= 0 (og-caller:observe-rank base 2)))
    (assert (= 0 (og-caller:observe-solve base rhs base-x)))
    (dolist (power '(-600 -31 29 600))
      (let ((scaled (og:make-decomposition 4 3)) (x (zero-vector 3)))
        (assert (= 0 (og-caller:observe-factor
                      scaled (rect-column-major
                              (rect 4 3 (lambda (i j) (* (aref matrix i j) (expt 2 power))))))))
        (assert (equalp (og:factor-permutation base) (og:factor-permutation scaled)))
        (assert (equalp (og:factor-tau base) (og:factor-tau scaled)) () "Scaling by 2^~D changed tau." power)
        (dotimes (j 3)
          (dotimes (i 4)
            (let ((b (aref (og:factor-work base) (+ i (* j 4))))
                  (v (aref (og:factor-work scaled) (+ i (* j 4)))))
              ;; On and above the diagonal is R, which scales; below it are reflector tails, which do not.
              (assert (= v (if (<= i j) (scale-float b power) b)) ()
                      "Scaling by 2^~D gave work[~D,~D] = ~S from ~S." power i j v b))))
        (assert (= 0 (og-caller:observe-rank scaled 2)))
        (assert (= 0 (og-caller:observe-solve scaled rhs x)))
        (dotimes (i 3)
          (assert (= (aref x i) (scale-float (aref base-x i) (- power))) ()
                  "Scaling by 2^~D gave x[~D] = ~S." power i (aref x i)))))
    ;; A repeated call repeats every bit, and a new rank does not need a new factorization.
    (let ((first (copied-decomposition base)) (x (zero-vector 3)))
      (assert (= 0 (og-caller:observe-factor base (rect-column-major matrix))))
      (assert (equalp first (copied-decomposition base)) () "A repeated factorization changed.")
      (assert (= 0 (og-caller:observe-rank base 3)))
      (assert (= 0 (og-caller:observe-rank base 2)))
      (assert (equalp first (copied-decomposition base)) () "Selecting a rank changed the factorization.")
      (assert (= 0 (og-caller:observe-solve base rhs x)))
      (assert (equalp x base-x) () "Returning to an earlier rank changed its solution.")))
  ;; Q^T undoes Q.
  (let ((decomposition (og:make-decomposition 3 2)) (v (doubles 1d0 2d0 3d0)))
    (assert (= 0 (og-caller:observe-factor decomposition (doubles 1d0 1d0 1d0  0d0 1d0 2d0))))
    (assert (= 0 (og-caller:observe-apply-q decomposition v)))
    (assert (= 0 (og-caller:observe-apply-qt decomposition v)))
    (assert (close-vector-p v (doubles 1d0 2d0 3d0))))
  (format t "ORTHOGONAL-IDENTITIES-PASS~%"))

(defun run-orthogonal-lifetime-tests ()
  (let* ((decomposition (og:make-decomposition 3 2))
         (a (doubles 1d0 1d0 1d0  0d0 1d0 2d0)) (rhs (doubles 1d0 3d0 5d0))
         (x (doubles 9d0 9d0)) (infinity sb-ext:double-float-positive-infinity))
    (flet ((untouched () (every (lambda (v) (= v 9d0)) x)))
      ;; Nothing stands before the first factorization.
      (assert (not (og:factor-ready? decomposition)))
      (assert (= 1 (og-caller:observe-rank decomposition 1)))
      (assert (= 1 (og-caller:observe-solve decomposition rhs x)))
      (assert (= 1 (og-caller:observe-nullspace decomposition 0 (zero-vector 2))))
      (assert (= 1 (og-caller:observe-apply-q decomposition (zero-vector 3))))
      (assert (= 1 (og-caller:observe-apply-qt decomposition (zero-vector 3))))
      (assert (= -1d0 (og:solve-residual-norm decomposition)))
      ;; A factor without a rank does not solve.
      (assert (= 0 (og-caller:observe-factor decomposition a)))
      (assert (and (og:factor-ready? decomposition) (not (og:rank-set? decomposition))))
      (assert (= 2 (og-caller:observe-solve decomposition rhs x)))
      (assert (= 2 (og-caller:observe-nullspace decomposition 0 (zero-vector 2))))
      (assert (untouched))
      (assert (= 0 (og-caller:observe-rank decomposition 2)))
      ;; Solve refusals leave the output alone and withdraw the residual.
      (assert (= 3 (og-caller:observe-solve decomposition (doubles 1d0 2d0) x)))
      (assert (= 3 (og-caller:observe-solve decomposition rhs (zero-vector 3))))
      (assert (= 5 (og-caller:observe-solve decomposition (doubles 1d0 infinity 3d0) x)))
      (assert (= 5 (og-caller:observe-solve decomposition (doubles 1d0 (not-a-number) 3d0) x)))
      (assert (untouched))
      (assert (= 0 (og-caller:observe-solve decomposition rhs x)))
      (assert (close-vector-p x (doubles 1d0 2d0)))
      (assert (<= 0d0 (og:solve-residual-norm decomposition) 1d-14))
      (assert (= 5 (og-caller:observe-solve decomposition (doubles 1d0 infinity 3d0) x)))
      (assert (= -1d0 (og:solve-residual-norm decomposition)) () "A refused solve kept an earlier residual.")
      ;; Full column rank has no nullspace vector to give.
      (assert (= 4 (og-caller:observe-nullspace decomposition 0 (zero-vector 2))))
      (assert (= 0 (og-caller:observe-rank decomposition 1)))
      (assert (= 3 (og-caller:observe-nullspace decomposition 0 (zero-vector 3))))
      (assert (= 0 (og-caller:observe-nullspace decomposition 0 (zero-vector 2))))
      (assert (= 4 (og-caller:observe-nullspace decomposition 1 (zero-vector 2))))
      ;; Reflector applications check their vector.
      (assert (= 3 (og-caller:observe-apply-q decomposition (zero-vector 2))))
      (assert (= 5 (og-caller:observe-apply-qt decomposition (doubles 1d0 infinity 3d0))))
      ;; Every factorization call withdraws the rank; a refused one withdraws the factor.
      (assert (= 0 (og-caller:observe-factor decomposition a)))
      (assert (not (og:rank-set? decomposition)) () "A new factorization kept the earlier rank.")
      (fill x 9d0)
      (dolist (attempt (list (lambda () (og-caller:observe-factor decomposition (doubles 1d0 2d0)))
                             (lambda () (let ((bad (copy-seq a))) (setf (aref bad 4) infinity)
                                          (og-caller:observe-factor decomposition bad)))
                             (lambda () (let ((bad (copy-seq a))) (setf (aref bad 0) (not-a-number))
                                          (og-caller:observe-factor decomposition bad)))))
        (assert (= 0 (og-caller:observe-factor decomposition a)))
        (assert (= 0 (og-caller:observe-rank decomposition 2)))
        (assert (/= 0 (funcall attempt)))
        (assert (and (not (og:factor-ready? decomposition)) (not (og:rank-set? decomposition))
                     (= 0 (og:selected-rank decomposition)) (= 0d0 (og:pivot-magnitude decomposition 0))
                     (= 0 (og:rank-above decomposition 0d0)))
                () "A refused factorization left something standing.")
        (assert (= 1 (og-caller:observe-solve decomposition rhs x)))
        (assert (untouched)))
      (assert (= 2 (og-caller:observe-factor decomposition (doubles 1d0 2d0))))
      (assert (= 1 (og-caller:observe-factor decomposition (let ((bad (copy-seq a))) (setf (aref bad 4) infinity) bad))))))
  ;; A column whose norm overflows is a refusal, not an infinite R.
  (let ((decomposition (og:make-decomposition 2 1)) (big (/ most-positive-double-float 1.2d0)))
    (assert (= 1 (og-caller:observe-factor decomposition (doubles big big))))
    (assert (not (og:factor-ready? decomposition))))
  ;; A finite decomposition whose solution overflows.
  (let ((decomposition (og:make-decomposition 1 1)) (x (doubles 9d0)))
    (assert (= 0 (og-caller:observe-factor decomposition (doubles (scale-float 1d0 -1000)))))
    (assert (= 0 (og-caller:observe-rank decomposition 1)))
    (assert (= 6 (og-caller:observe-solve decomposition (doubles (scale-float 1d0 100)) x)))
    (assert (= 9d0 (aref x 0))))
  ;; Storage is refused by name before it is allocated.
  (let ((limit (coalton:coalton og:maximum-dimension)) (budget (coalton:coalton og:default-entry-limit)))
    (assert (= 0 (og-caller:storage-error 3 2 6)))
    (assert (= 4 (og-caller:storage-error 3 2 5)))
    (assert (= 0 (og-caller:storage-error 0 5 0)))
    (assert (= 3 (og-caller:storage-error (1+ limit) 1 budget)))
    (assert (= 3 (og-caller:storage-error 1 (1+ limit) budget)))
    (assert (= 4 (og-caller:storage-error limit limit budget)))
    (assert (= 4 (og-caller:storage-error 4097 4096 budget))))
  (format t "ORTHOGONAL-LIFETIME-PASS~%"))

(defparameter *orthogonal-public-control*
  "(declare inspect (coalton-linalg/dense/orthogonal:Decomposition -> UFix))
   (define (inspect decomposition)
     (if (coalton-linalg/dense/orthogonal:factor-ready? decomposition)
         1
         (coalton-linalg/dense/orthogonal:step-count decomposition)))
   (declare probe (Void -> UFix))
   (define (probe)
     (inspect (coalton-linalg/dense/orthogonal:make-decomposition 3 2)))")

(defparameter *orthogonal-encapsulation-probes*
  '(("orthogonal-construction" .
     "(declare probe (Void -> coalton-linalg/dense/orthogonal:Decomposition))
      (define (probe)
        (coalton-linalg/dense/orthogonal:Decomposition
          0 0 0 (arr:make 0 0.0d0) (arr:make 0 0.0d0) (arr:make 0 (the UFix 0)) (arr:make 0 0.0d0)
          (arr:make 0 0.0d0) (arr:make 0 0.0d0) (arr:make 0 0.0d0) (arr:make 3 (the UFix 0))
          (arr:make 1 0.0d0)))")
    ("orthogonal-destructuring" .
     "(declare probe (coalton-linalg/dense/orthogonal:Decomposition -> UFix))
      (define (probe decomposition)
        (match decomposition
          ((coalton-linalg/dense/orthogonal:Decomposition _ _ steps _ _ _ _ _ _ _ _ _) steps)))")
    ("orthogonal-state-access" .
     "(declare probe (coalton-linalg/dense/orthogonal:Decomposition -> UFix))
      (define (probe decomposition) (arr:aref (.state decomposition) 0))")
    ("orthogonal-trapezoid-access" .
     "(declare probe (coalton-linalg/dense/orthogonal:Decomposition -> arr:LispArray F64))
      (define (probe decomposition) (.trapezoid decomposition))")
    ("orthogonal-status-confusion" .
     "(declare probe (coalton-linalg/dense/orthogonal:RankStatus -> coalton-linalg/dense/orthogonal:SolveStatus))
      (define (probe status) status)")))

(defun run-orthogonal-encapsulation-tests ()
  "External Coalton callers can use the decomposition only through its public operations."
  (multiple-value-bind (accepted diagnostics)
      (coalton-linalg/tests/encapsulation::compile-probe "orthogonal-public-control" *orthogonal-public-control*
                                                         :execute t)
    (assert accepted () "The valid orthogonal external caller failed:~%~A" diagnostics))
  (dolist (probe *orthogonal-encapsulation-probes*)
    (assert (not (coalton-linalg/tests/encapsulation::compile-probe (car probe) (cdr probe))) ()
            "External caller ~A compiled despite required rejection." (car probe)))
  (dolist (name '("MKDECOMPOSITION" "DECOMPOSITIONSTORAGE"))
    (multiple-value-bind (symbol visibility) (find-symbol name "COALTON-LINALG/DENSE/ORTHOGONAL")
      (assert (and symbol (eq visibility :internal)) ()
              "Storage constructor ~A must exist and remain internal; visibility is ~S." name visibility)))
  (format t "ORTHOGONAL-ENCAPSULATION-PASS rejected=~D public-control=1~%"
          (length *orthogonal-encapsulation-probes*)))

(defun run-orthogonal-tests ()
  (run-orthogonal-exact-tests)
  (run-orthogonal-rank-tests)
  (run-orthogonal-identity-tests)
  (run-orthogonal-lifetime-tests)
  (run-orthogonal-encapsulation-tests)
  ;; Warmed repeated work through a compiled Coalton caller allocates nothing.
  (let* ((state (sb-ext:seed-random-state 77))
         (matrix (seeded-product state 9 12 6))
         (a (rect-column-major matrix)) (decomposition (og:make-decomposition 9 12))
         (rhs (make-array 9 :element-type 'double-float :initial-element 1d0))
         (x (zero-vector 12)) (basis (zero-vector 12)))
    (assert (= 0 (og-caller:pipeline-loop decomposition a rhs x basis 1d-9 3)) ()
            "The allocation case does not run.")
    (assert (= 6 (og:selected-rank decomposition)) () "The allocation case has rank ~D."
            (og:selected-rank decomposition))
    (allocation-check "orthogonal-pipeline"
                      (lambda () (og-caller:pipeline-loop decomposition a rhs x basis 1d-9 1))))
  (format t "ORTHOGONAL-TESTS-PASS~%"))
