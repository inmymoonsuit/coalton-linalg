(in-package #:coalton-linalg/tests)

;;; The sparse LU factorization with partial pivoting against exact rational arithmetic on the
;;; stored F64 values. The reference shares no code with the library: the exported factors are
;;; laid out densely and multiplied out exactly, every entry of P A Q - L U is compared with the
;;; rounding analysis of the kernel, and the pivot rule is read off the multipliers.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "LU" "COALTON-LINALG/SPARSE/LU")
  (sb-ext:add-package-local-nickname "LU-CALLER" "COALTON-LINALG/TESTS/LU-CALLERS"))

(defvar *lu-worst-entry* 0)
(defvar *lu-worst-backward* 0)
(defvar *lu-worst-solve* 0)
(defvar *lu-worst-transpose* 0)
(defvar *lu-largest-growth* 0)

(defun gamma (count)
  "count u / (1 - count u), exactly."
  (/ (* count *unit-roundoff*) (- 1 (* count *unit-roundoff*))))

(defun lu-dense-factors (handle n)
  "L with its unit diagonal and U as exact rational matrices in pivot order, the number of entries of U above the
diagonal in each column, and the pivot step of each row of the matrix. Asserts the structure the contract states."
  (let* ((rows (lu:row-permutation handle))
         (lp (lu:lower-pointers handle)) (li (lu:lower-rows handle)) (lx (lu:lower-values handle))
         (up (lu:upper-pointers handle)) (ui (lu:upper-rows handle)) (ux (lu:upper-values handle))
         (diagonal (lu:upper-diagonal handle))
         (step (make-array n :initial-element nil))
         (l (build-matrix n (lambda (i j) (if (= i j) 1 0))))
         (u (build-matrix n (lambda (i j) (declare (ignore i j)) 0)))
         (above (make-array n :initial-element 0)))
    (dotimes (k n)
      (assert (and (< (aref rows k) n) (null (aref step (aref rows k)))) () "Row ~D is pivoted twice." (aref rows k))
      (setf (aref step (aref rows k)) k))
    (assert (and (= 0 (aref lp 0)) (= 0 (aref up 0))))
    (assert (= (aref lp n) (lu:lower-entry-count handle)))
    (assert (= (aref up n) (lu:upper-entry-count handle)))
    (dotimes (k n)
      (setf (aref u k k) (rational (aref diagonal k)))
      (assert (/= 0 (aref diagonal k)) () "Pivot ~D is zero." k)
      (loop :for p :from (aref lp k) :below (aref lp (1+ k)) :do
        (let ((i (aref step (aref li p))))
          (assert (> i k) () "Column ~D of L holds row ~D, pivoted at step ~D." k (aref li p) i)
          (assert (zerop (aref l i k)) () "Column ~D of L holds row ~D twice." k (aref li p))
          ;; A stored zero is told from an absent entry by this marker, removed below.
          (setf (aref l i k) (if (zerop (aref lx p)) :zero (rational (aref lx p))))))
      (setf (aref above k) (- (aref up (1+ k)) (aref up k)))
      (loop :for p :from (aref up k) :below (aref up (1+ k)) :do
        (let ((i (aref ui p)))
          (assert (< i k) () "Column ~D of U holds step ~D." k i)
          (assert (zerop (aref u i k)) () "Column ~D of U holds step ~D twice." k i)
          (setf (aref u i k) (if (zerop (aref ux p)) :zero (rational (aref ux p)))))))
    (dotimes (i n)
      (dotimes (j n)
        (when (eq (aref l i j) :zero) (setf (aref l i j) 0))
        (when (eq (aref u i j) :zero) (setf (aref u i j) 0))))
    (values l u above step)))

(defun check-lu (label n entries &key (threshold 1d0) order (limit 100000) (state (sb-ext:seed-random-state 3)))
  "Factor and hold every statement of docs/lu.md against exact arithmetic. Returns :singular when the matrix was
refused as singular, after checking that it is; the handle otherwise."
  (let* ((matrix (csc-from-entries n n entries)) (a (dr-dense n n entries))
         (handle (if order
                     (lu:make-lu-with-order n (apply #'indices order) limit)
                     (lu:make-lu matrix)))
         (code (lu-caller:observe-factor handle matrix threshold)))
    (assert (and (= n (lu:dimension handle)) (not (and (/= code 0) (lu:lu-ready? handle)))))
    (when (= code 1)
      (assert (< (dr-rank a) n) () "~A: refused as singular, and the exact rank is ~D of ~D." label (dr-rank a) n)
      (assert (< (lu:singular-column handle) n))
      (assert (and (= -1d0 (lu:lu-growth handle)) (= -1d0 (lu:lu-backward-bound handle))
                   (= 0 (lu:lower-entry-count handle)) (= 0 (lu:upper-entry-count handle))))
      (return-from check-lu :singular))
    (assert (= code 0) () "~A: factorization outcome ~D." label code)
    (assert (and (lu:lu-ready? handle) (= n (lu:singular-column handle))))
    (multiple-value-bind (l u above step) (lu-dense-factors handle n)
      (let* ((rows (lu:row-permutation handle)) (columns (lu:column-permutation handle))
             (product (matrix-product l u))
             (bound (matrix-product (matrix-map #'abs l) (matrix-map #'abs u)))
             (off 0) (total 0) (longest (reduce #'max above :initial-value 0)))
        (assert (equalp (sort (copy-seq columns) #'<) (coerce (loop :for i :below n :collect i) 'vector)))
        (when order (assert (equalp columns (coerce order 'vector))))
        ;; Each entry of P A Q - L U against gamma(products + 1) |L| |U|: the products that entry
        ;; received are at most the entries of U above the diagonal in its column, and an entry
        ;; of L was divided once more.
        (dotimes (i n)
          (dotimes (j n)
            (let* ((residual (abs (- (aref a (aref rows i) (aref columns j)) (aref product i j))))
                   (allowed (* (gamma (1+ (aref above j))) (aref bound i j))))
              (assert (<= residual allowed) () "~A: entry ~D,~D is off by ~E; ~E is allowed." label i j
                      (float residual 1d0) (float allowed 1d0))
              (when (plusp allowed) (setf *lu-worst-entry* (max *lu-worst-entry* (/ residual allowed))))
              (incf off (expt residual 2))
              (incf total (aref bound i j)))))
        ;; The pivot rule, read off the multipliers. No multiplier exceeds one over the
        ;; threshold; with a threshold of one a multiplier of magnitude exactly one is a tie,
        ;; which went to the diagonal of the matrix or to the smaller row.
        (dotimes (k n)
          (loop :for i :from (1+ k) :below n :do
            (assert (<= (abs (aref l i k)) (* (/ 1 (rational threshold)) (+ 1 (expt 2 -50)))) ()
                    "~A: multiplier ~E in column ~D exceeds one over the threshold." label (float (aref l i k) 1d0) k)
            (when (and (= threshold 1d0) (= 1 (abs (aref l i k))))
              (assert (or (= (aref rows k) (aref columns k)) (< (aref rows k) (aref rows i))) ()
                      "~A: column ~D: row ~D tied with the pivot row ~D and is the smaller." label k (aref rows i)
                      (aref rows k)))))
        (assert (every #'integerp step))
        ;; Growth: the largest magnitude in U over the largest in the matrix.
        (let* ((top (rect-max a)) (growth (if (plusp top) (/ (matrix-max u) top) 0)))
          (setf *lu-largest-growth* (max *lu-largest-growth* growth))
          (assert (<= (abs (- (rational (lu:lu-growth handle)) growth)) (* (expt 2 -40) (+ 1 growth))) ()
                  "~A: reported growth ~S, exact ~E." label (lu:lu-growth handle) (float growth 1d0)))
        ;; The backward bound: at least the exact Frobenius norm, and the contract's formula on
        ;; the exported factors, rounded up and by no more than rounding.
        (let* ((reported (rational (lu:lu-backward-bound handle)))
               (stored (+ n (lu:lower-entry-count handle) (lu:upper-entry-count handle)))
               (formula (+ (* (gamma (1+ longest)) total)
                           (* (rational least-positive-double-float) stored (+ longest 2 (matrix-max u))))))
          (assert (<= off (expt reported 2)) () "~A: the backward error is ~E and the reported bound ~E." label
                  (sqrt (float off 1d0)) (float reported 1d0))
          (assert (<= formula reported (+ (* formula (+ 1 (expt 2 -30))) (* 8 (rational least-positive-double-float))))
                  () "~A: reported bound ~E, formula ~E." label (float reported 1d0) (float formula 1d0))
          (when (plusp reported) (setf *lu-worst-backward* (max *lu-worst-backward* (/ off (expt reported 2))))))
        ;; Solves, each row of the exact residual against gamma(3n) |L| |U| |x|, and its
        ;; transposed form; and the same bits in place.
        (let ((b (zeros n)) (x (make-array n :element-type 'double-float :initial-element 7d0))
              (y (make-array n :element-type 'double-float :initial-element 7d0)))
          (dotimes (i n) (setf (aref b i) (coerce (- (random 9 state) 4) 'double-float)))
          (assert (= 0 (lu-caller:observe-solve handle b x)) () "~A: solve refused." label)
          (dotimes (i n)
            (let* ((row (aref rows i))
                   (residual (abs (- (loop :for j :below n :sum (* (aref a row j) (rational (aref x j))))
                                     (rational (aref b row)))))
                   (allowed (* (gamma (* 3 n))
                               (loop :for j :below n
                                     :sum (* (aref bound i j) (abs (rational (aref x (aref columns j)))))))))
              (assert (<= residual allowed) () "~A: solve residual ~E in row ~D exceeds ~E." label
                      (float residual 1d0) row (float allowed 1d0))
              (when (plusp allowed) (setf *lu-worst-solve* (max *lu-worst-solve* (/ residual allowed))))))
          (assert (= 0 (lu-caller:observe-transpose-solve handle b y)) () "~A: transposed solve refused." label)
          (dotimes (j n)
            (let* ((column (aref columns j))
                   (residual (abs (- (loop :for i :below n :sum (* (aref a i column) (rational (aref y i))))
                                     (rational (aref b column)))))
                   (allowed (* (gamma (* 3 n))
                               (loop :for i :below n
                                     :sum (* (aref bound i j) (abs (rational (aref y (aref rows i)))))))))
              (assert (<= residual allowed) () "~A: transposed solve residual ~E in column ~D exceeds ~E." label
                      (float residual 1d0) column (float allowed 1d0))
              (when (plusp allowed) (setf *lu-worst-transpose* (max *lu-worst-transpose* (/ residual allowed))))))
          (let ((in-place (copy-seq b)))
            (assert (= 0 (lu-caller:observe-solve handle in-place in-place)))
            (assert (same-bits-p x in-place) () "~A: the solve in place differs." label)
            (setf in-place (copy-seq b))
            (assert (= 0 (lu-caller:observe-transpose-solve handle in-place in-place)))
            (assert (same-bits-p y in-place) () "~A: the transposed solve in place differs." label)))))
    handle))

(defun seeded-square (state n density &key (diagonal t) hub)
  "Stored entries of an integer n-by-n matrix; DIAGONAL stores every diagonal entry, HUB the whole first row."
  (loop :for i :below n :append
    (loop :for j :below n
          :when (or (and diagonal (= i j)) (and hub (= i 0)) (< (random 100 state) density))
            :collect (list i j (let ((v (- (random 9 state) 4))) (if (zerop v) 2 v))))))

(defun run-lu-exact-tests ()
  (setf *lu-worst-entry* 0 *lu-worst-backward* 0 *lu-worst-solve* 0 *lu-worst-transpose* 0 *lu-largest-growth* 0)
  (let ((state (sb-ext:seed-random-state 20261205)) (cases 0) (singular 0) (exchanged 0))
    (flet ((try (label n entries &rest keys)
             (let ((handle (apply #'check-lu label n entries :state state keys)))
               (incf cases)
               (cond ((eq handle :singular) (incf singular))
                     ((notevery #'= (lu:row-permutation handle) (lu:column-permutation handle)) (incf exchanged)))
               handle)))
      (try "one by one" 1 '((0 0 -3)))
      (try "identity" 3 '((0 0 1) (1 1 1) (2 2 1)))
      ;; No entry on the diagonal anywhere: every pivot is an exchange.
      (try "cyclic permutation" 4 '((1 0 2) (2 1 3) (3 2 5) (0 3 7)))
      (try "upper triangular" 3 '((0 0 2) (0 1 1) (0 2 4) (1 1 3) (1 2 5) (2 2 6)))
      (try "lower triangular" 3 '((0 0 2) (1 0 1) (2 0 4) (1 1 3) (2 1 5) (2 2 6)))
      (try "stored zeros" 3 '((0 0 2) (1 0 0) (2 0 4) (1 1 3) (2 1 0) (2 2 6) (0 2 0)))
      ;; In the natural order the first column holds 1 on the diagonal and 2 below it. At a
      ;; threshold of one half the diagonal is taken; at six tenths the larger entry is.
      (let ((entries '((0 0 1) (1 0 2) (0 1 3) (1 1 1))))
        (assert (equalp #(0 1) (lu:row-permutation (try "diagonal within the threshold" 2 entries
                                                        :order '(0 1) :threshold 0.5d0))))
        (assert (equalp #(1 0) (lu:row-permutation (try "diagonal outside the threshold" 2 entries
                                                        :order '(0 1) :threshold 0.6d0))))
        (assert (equalp #(1 0) (lu:row-permutation (try "strict partial pivoting" 2 entries :order '(0 1))))))
      ;; A tie between two rows neither of which is the diagonal's goes to the smaller row; a
      ;; tie with the diagonal goes to the diagonal.
      (assert (equalp #(1 0 2) (lu:row-permutation (try "tie to the smaller row" 3
                                                        '((1 0 4) (2 0 -4) (0 1 1) (1 1 1) (2 2 1) (0 2 3))
                                                        :order '(0 1 2)))))
      (assert (equalp #(1 0) (lu:row-permutation (try "tie to the diagonal" 2 '((0 1 -4) (1 1 4) (0 0 1) (1 0 2))
                                                      :order '(1 0)))))
      ;; Dense, large enough that the factors grow past their first allocation.
      (let ((dense (loop :for i :below 8 :append
                     (loop :for j :below 8 :collect (list i j (+ (* 3 (mod (* (1+ i) (+ 2 j)) 7)) (if (= i j) 20 0)))))))
        (assert (> (lu:upper-entry-count (try "dense, growing storage" 8 dense)) 16)))
      (loop :for n :from 2 :to 8 :do
        (dotimes (repeat 8)
          (let ((entries (seeded-square state n (nth (mod repeat 3) '(25 45 75))
                                        :diagonal (/= 0 (mod repeat 4)) :hub (= 5 repeat))))
            (try (format nil "seeded ~D #~D" n repeat) n entries)
            (try (format nil "seeded ~D #~D, natural order" n repeat) n entries
                 :order (loop :for j :below n :collect j))
            (try (format nil "seeded ~D #~D, reversed, threshold 0.1" n repeat) n entries
                 :order (loop :for j :from (1- n) :downto 0 :collect j) :threshold 0.1d0)))))
    (assert (and (> cases 150) (> exchanged 60) (> singular 3)) ()
            "The LU roster has ~D cases, ~D with a row exchange and ~D singular." cases exchanged singular)
    (format t "LU-EXACT-PASS cases=~D with-exchanges=~D singular=~D worst-entry=~,3F worst-backward=~,3F worst-solve=~,3F worst-transpose=~,3F largest-growth=~,2F~%"
            cases exchanged singular (float *lu-worst-entry* 1d0) (float *lu-worst-backward* 1d0)
            (float *lu-worst-solve* 1d0) (float *lu-worst-transpose* 1d0) (float *lu-largest-growth* 1d0))))

(defun run-lu-outcome-tests ()
  (let ((nan (not-a-number)) (infinity sb-ext:double-float-positive-infinity))
    ;; Singular matrices, and the column that showed it.
    (flet ((singular (label n entries order column)
             (let* ((matrix (csc-from-entries n n entries))
                    (handle (lu:make-lu-with-order n (apply #'indices order) 1000)))
               (assert (= 1 (lu-caller:observe-factor handle matrix 1d0)) () "~A was not refused as singular." label)
               (assert (= column (lu:singular-column handle)) () "~A: singular column ~D, expected ~D." label
                       (lu:singular-column handle) column)
               (assert (not (lu:lu-ready? handle)))
               (assert (= 1 (lu-caller:observe-solve handle (zeros n) (zeros n))))
               (assert (= 1 (lu-caller:observe-transpose-solve handle (zeros n) (zeros n)))))))
      (singular "A column with no entry" 3 '((0 0 1) (1 0 2) (2 2 3)) '(0 1 2) 1)
      (singular "A column of stored zeros" 2 '((0 0 0) (1 0 0) (0 1 1) (1 1 1)) '(0 1) 0)
      (singular "A repeated column" 2 '((0 0 1) (1 0 1) (0 1 1) (1 1 1)) '(0 1) 1)
      (singular "A repeated column, the other order" 2 '((0 0 1) (1 0 1) (0 1 1) (1 1 1)) '(1 0) 0)
      (singular "A row with no entry" 3 '((0 0 1) (0 1 2) (2 1 1) (0 2 5) (2 2 3)) '(0 1 2) 2))
    ;; An entry that overflows as the second column is formed: 1e308 + 1e308.
    (let* ((matrix (csc-from-entries 2 2 '((0 0 1) (1 0 -1) (0 1 1d308) (1 1 1d308))))
           (handle (lu:make-lu-with-order 2 (indices 0 1) 1000)))
      (assert (= 2 (lu-caller:observe-factor handle matrix 1d0)))
      (assert (not (lu:lu-ready? handle)))
      ;; The same handle then factors a sound matrix as a fresh handle does, bit for bit: the
      ;; refused call left nothing behind in the workspace.
      (let* ((sound (csc-from-entries 2 2 '((0 0 3) (1 0 1) (0 1 2) (1 1 5))))
             (fresh (lu:make-lu-with-order 2 (indices 0 1) 1000)))
        (assert (= 0 (lu-caller:observe-factor handle sound 1d0)))
        (assert (= 0 (lu-caller:observe-factor fresh sound 1d0)))
        (assert (same-bits-p (lu:upper-diagonal handle) (lu:upper-diagonal fresh)))
        (assert (same-bits-p (subseq (lu:upper-values handle) 0 1) (subseq (lu:upper-values fresh) 0 1)))
        (assert (same-bits-p (subseq (lu:lower-values handle) 0 1) (subseq (lu:lower-values fresh) 0 1)))))
    ;; An entry of U that overflows while every candidate for the pivot stays finite: the third
    ;; column holds 1e308 in the two rows already pivoted, the second of which takes it twice.
    (let* ((matrix (csc-from-entries 3 3 '((0 0 1) (1 0 -1) (1 1 1) (0 2 1d308) (1 2 1d308) (2 2 1))))
           (handle (lu:make-lu-with-order 3 (indices 0 1 2) 1000)))
      (assert (= 2 (lu-caller:observe-factor handle matrix 1d0)) () "An infinite entry of U was accepted.")
      (assert (not (lu:lu-ready? handle))))
    ;; The budget counts the entries of both factors off their diagonals: a dense matrix of
    ;; order four has twelve.
    (let* ((entries (loop :for i :below 4 :append (loop :for j :below 4 :collect (list i j (+ 1 (* i 4) j (if (= i j) 30 0))))))
           (matrix (csc-from-entries 4 4 entries)))
      (loop :for limit :from 0 :to 11 :do
        (let ((handle (lu:make-lu-with-order 4 (indices 0 1 2 3) limit)))
          (assert (= limit (lu:entry-limit handle)))
          (assert (= 3 (lu-caller:observe-factor handle matrix 1d0)) () "A budget of ~D held twelve entries." limit)
          (assert (not (lu:lu-ready? handle)))))
      (let ((handle (lu:make-lu-with-order 4 (indices 0 1 2 3) 12)))
        (assert (= 0 (lu-caller:observe-factor handle matrix 1d0)))
        (assert (= 12 (+ (lu:lower-entry-count handle) (lu:upper-entry-count handle))))
        ;; A handle that ran out factors a smaller matrix afterwards, and this one again.
        (let ((tight (lu:make-lu-with-order 4 (indices 0 1 2 3) 11))
              (diagonal (csc-from-entries 4 4 '((0 0 1) (1 1 2) (2 2 3) (3 3 4)))))
          (assert (= 3 (lu-caller:observe-factor tight matrix 1d0)))
          (assert (= 0 (lu-caller:observe-factor tight diagonal 1d0)))
          (assert (= 0 (+ (lu:lower-entry-count tight) (lu:upper-entry-count tight)))))
        ;; Wrong shapes and thresholds leave no factors, whatever stood before.
        (dolist (threshold (list 0d0 -1d0 1.5d0 nan infinity))
          (assert (= 0 (lu-caller:observe-factor handle matrix 1d0)))
          (assert (= 5 (lu-caller:observe-factor handle matrix threshold)) () "Threshold ~S was accepted." threshold)
          (assert (and (not (lu:lu-ready? handle)) (= -1d0 (lu:lu-backward-bound handle))
                       (= 0 (lu:lower-entry-count handle)))))
        (assert (= 0 (lu-caller:observe-factor handle matrix (scale-float 1d0 -1074))))
        (dolist (other (list (csc-from-entries 3 3 '((0 0 1) (1 1 1) (2 2 1)))
                             (csc-from-entries 4 5 '((0 0 1)))
                             (csc-from-entries 5 4 '((0 0 1)))))
          (assert (= 0 (lu-caller:observe-factor handle matrix 1d0)))
          (assert (= 4 (lu-caller:observe-factor handle other 1d0)))
          (assert (not (lu:lu-ready? handle))))
        ;; Solves on standing factors: sizes and right-hand sides that are not finite.
        (assert (= 0 (lu-caller:observe-factor handle matrix 1d0)))
        (let ((out (doubles 9d0 9d0 9d0 9d0)) (untouched (doubles 9d0 9d0 9d0 9d0)) (b (doubles 1d0 2d0 3d0 4d0)))
          (dolist (solve (list #'lu-caller:observe-solve #'lu-caller:observe-transpose-solve))
            (assert (= 2 (funcall solve handle (doubles 1d0 2d0 3d0) out)))
            (assert (= 2 (funcall solve handle b (zeros 5))))
            (assert (= 3 (funcall solve handle (doubles 1d0 nan 3d0 4d0) out)))
            (assert (= 3 (funcall solve handle (doubles 1d0 2d0 3d0 infinity) out)))
            (assert (same-bits-p out untouched) () "A refused solve wrote its output.")))))
    ;; A solution that overflows: a pivot of 1e-320 under a right-hand side of 1.
    (let* ((matrix (csc-from-entries 1 1 '((0 0 1d-320)))) (handle (lu:make-lu matrix)) (out (doubles 9d0)))
      (assert (= 0 (lu-caller:observe-factor handle matrix 1d0)))
      (assert (= 4 (lu-caller:observe-solve handle (doubles 1d0) out)))
      (assert (= 4 (lu-caller:observe-transpose-solve handle (doubles 1d0) out)))
      (assert (= 9d0 (aref out 0))))
    ;; Typed errors of construction.
    (assert (= 92 (lu-caller:make-code (csc-from-entries 2 3 '((0 0 1) (1 1 1))))))
    (assert (= 0 (lu-caller:make-code (csc-from-entries 0 0 '()))))
    (assert (= 93 (lu-caller:make-ordered-code 3 (indices 0 1) 100)))
    (assert (= 93 (lu-caller:make-ordered-code 3 (indices 0 1 2 2) 100)))
    (assert (= 93 (lu-caller:make-ordered-code 3 (indices 0 1 3) 100)))
    (assert (= 93 (lu-caller:make-ordered-code 3 (indices 0 1 1) 100)))
    (assert (= 0 (lu-caller:make-ordered-code 3 (indices 2 0 1) 100)))
    (assert (= 98 (lu-caller:order-code (csc-from-entries 3 3 '((0 0 1) (0 1 1) (0 2 1) (1 1 1) (2 2 1))) 4)))
    (assert (= 0 (lu-caller:order-code (csc-from-entries 3 3 '((0 0 1) (0 1 1) (0 2 1) (1 1 1) (2 2 1))) 100))))
  ;; Underflow. Every entry of this matrix is below 2^-1039, so the relative term of the bound,
  ;; which is a multiple of the sum of |L| |U|, is below 2^-1090. The one product of the
  ;; factorization is 2^-10 times 2^-1066, a quarter of the smallest positive number, and
  ;; rounds to zero: the error is 2^-1076, and the bound holds through its underflow term alone.
  (let* ((least least-positive-double-float)
         (matrix (csc-from-entries 2 2 (list (list 0 0 (scale-float 1d0 -1040)) (list 1 0 (scale-float 1d0 -1050))
                                             (list 0 1 (scale-float 1d0 -1066)) (list 1 1 least))))
         (handle (lu:make-lu-with-order 2 (indices 0 1) 100)))
    (assert (= 0 (lu-caller:observe-factor handle matrix 1d0)))
    (assert (equalp #(0 1) (lu:row-permutation handle)))
    (multiple-value-bind (l u) (lu-dense-factors handle 2)
      (let* ((a (dr-dense 2 2 (list (list 0 0 (scale-float 1d0 -1040)) (list 1 0 (scale-float 1d0 -1050))
                                    (list 0 1 (scale-float 1d0 -1066)) (list 1 1 least))))
             (product (matrix-product l u))
             (total (loop :for i :below 2 :sum (loop :for j :below 2
                                                     :sum (aref (matrix-product (matrix-map #'abs l) (matrix-map #'abs u)) i j))))
             (error (sqrt (float (loop :for i :below 2 :sum (loop :for j :below 2 :sum (expt (- (aref a i j) (aref product i j)) 2)))
                                 1d0)))
             (exact (loop :for i :below 2 :sum (loop :for j :below 2 :sum (expt (- (aref a i j) (aref product i j)) 2)))))
        (declare (ignorable error))
        (assert (= exact (expt 2 -2152)) () "The underflow case's premise failed: the squared error is ~S." exact)
        (assert (> exact (expt (* (gamma 2) total) 2)) ()
                "The underflow case's premise failed: the relative term alone covers the error.")
        (assert (>= (expt (rational (lu:lu-backward-bound handle)) 2) exact) ()
                "The backward bound ~S is below the error that underflow made." (lu:lu-backward-bound handle)))))
  (format t "LU-OUTCOMES-PASS~%"))

(defparameter *lu-public-control*
  "(declare inspect (coalton-linalg/sparse/lu:Lu -> UFix))
   (define (inspect handle)
     (+ (coalton-linalg/sparse/lu:dimension handle) (coalton-linalg/sparse/lu:lower-entry-count handle)))
   (declare probe (coalton-linalg/sparse/lu:Lu -> UFix))
   (define (probe handle) (inspect handle))")

(defparameter *lu-encapsulation-probes*
  '(("lu-work-access" .
     "(declare probe (coalton-linalg/sparse/lu:Lu -> arr:LispArray F64))
      (define (probe handle) (.x handle))")
    ("lu-unwrapping" .
     "(declare probe (coalton-linalg/sparse/lu:Lu -> UFix))
      (define (probe handle)
        (match handle ((coalton-linalg/sparse/lu:MkLu storage) (.size storage))))")
    ("lu-status-as-solve-status" .
     "(declare probe (coalton-linalg/sparse/lu:LuStatus -> coalton-linalg/sparse/lu:LuSolveStatus))
      (define (probe status) status)")
    ("lu-solve-status-as-ldlt-solve-status" .
     "(declare probe (coalton-linalg/sparse/lu:LuSolveStatus -> ldlt:SolveStatus))
      (define (probe status) status)")))

(defun run-lu-lifetime-tests ()
  ;; One handle, several matrices of its order, of other patterns too: each as a fresh handle gives it.
  (let* ((order (indices 2 0 3 1)) (handle (lu:make-lu-with-order 4 order 1000))
         (first (csc-from-entries 4 4 '((0 0 4) (1 0 1) (1 1 5) (2 1 2) (2 2 6) (3 2 1) (3 3 7) (0 3 2))))
         (second (csc-from-entries 4 4 '((0 1 3) (1 0 2) (2 3 5) (3 2 4) (0 0 1) (3 3 1) (2 0 6)))))
    (dolist (matrix (list first second first))
      (let ((fresh (lu:make-lu-with-order 4 order 1000)))
        (assert (= 0 (lu-caller:observe-factor handle matrix 1d0)))
        (assert (= 0 (lu-caller:observe-factor fresh matrix 1d0)))
        (assert (equalp (lu:row-permutation handle) (lu:row-permutation fresh)))
        (assert (same-bits-p (lu:upper-diagonal handle) (lu:upper-diagonal fresh)))
        (assert (= (lu:lower-entry-count handle) (lu:lower-entry-count fresh)))
        (assert (same-bits-p (subseq (lu:lower-values handle) 0 (lu:lower-entry-count handle))
                             (subseq (lu:lower-values fresh) 0 (lu:lower-entry-count fresh))))
        (assert (same-bits-p (subseq (lu:upper-values handle) 0 (lu:upper-entry-count handle))
                             (subseq (lu:upper-values fresh) 0 (lu:upper-entry-count fresh)))))))
  ;; Steady state: a compiled caller allocates nothing once the factors have grown to their size.
  (let* ((entries (loop :for i :below 8 :append
                    (loop :for j :below 8 :collect (list i j (+ (* 3 (mod (* (1+ i) (+ 2 j)) 7)) (if (= i j) 20 0))))))
         (matrix (csc-from-entries 8 8 entries)) (handle (lu:make-lu matrix))
         (rhs (doubles 1d0 2d0 3d0 4d0 5d0 6d0 7d0 8d0)) (out (zeros 8)))
    (allocation-check "lu-factor-solve-loop"
                      (lambda () (assert (= 0 (lu-caller:factor-solve-loop handle matrix 1d0 rhs out 1))))))
  (multiple-value-bind (accepted diagnostics)
      (coalton-linalg/tests/encapsulation::compile-probe "lu-public-control" *lu-public-control*)
    (assert accepted () "The valid LU external caller failed:~%~A" diagnostics))
  (dolist (probe *lu-encapsulation-probes*)
    (assert (not (coalton-linalg/tests/encapsulation::compile-probe (car probe) (cdr probe))) ()
            "External caller ~A compiled despite required rejection." (car probe)))
  (dolist (name '("MKLU" "LUSTORAGE"))
    (multiple-value-bind (symbol visibility) (find-symbol name "COALTON-LINALG/SPARSE/LU")
      (assert (and symbol (eq visibility ':internal)) ()
              "Storage constructor ~A must exist and remain internal; visibility is ~S." name visibility)))
  (format t "LU-LIFETIME-PASS rejected=~D public-control=1~%" (length *lu-encapsulation-probes*)))

(defun run-lu-tests ()
  (run-lu-exact-tests)
  (run-lu-outcome-tests)
  (run-lu-lifetime-tests)
  (format t "LU-TESTS-PASS~%"))
