(in-package #:coalton-linalg/tests)

;;; The scale-sweep roster. One question is asked of every exported call that takes or returns
;;; floating-point numbers, or returns a status: what does it do with the same problem at
;;; another scale? Multiplying every input of a problem by a power of two changes no
;;; significand, so each output has one right value at each scale: the output at scale 1
;;; times a known power of the factor. A call is run across the whole range of the format and
;;; each run is put in one class:
;;;
;;;   :exact    every output is that value, bit for bit;
;;;   :close    every output is within the entry's tolerance of it, measured against the largest
;;;             magnitude in its array;
;;;   :looser   an upper bound came out above its scaled value, which a bound may;
;;;   :lost     the call succeeded and an output's right value is too small to be a normal
;;;             number: what came back is zero or has lost digits;
;;;   :refused  the call refused, and some output's right value cannot be held at this scale;
;;;   :limit    the call refused although every output's right value can be held: the scale is
;;;             outside what the routine's own arithmetic reaches;
;;;   :no-input an input of the problem is not a normal number at this scale;
;;;   :wrong    the call succeeded and an output is neither its right value nor lost.
;;;
;;; :wrong fails wherever it appears. The other classes are declared for each entry as the
;;; ranges of the sweep they hold over, and the roster fails when a run leaves its declared
;;; class in either direction: a range that shrinks is a regression, and one that grows is
;;; an improvement the entry has to be changed to record.
;;;
;;; The members are read from the sources, as the aliasing roster's are, and each entry names
;;; the calls it answers for. A name is accepted only if the entry's own code calls it.

(defparameter *sweep-exponents*
  '(-1020 -1000 -800 -600 -520 -512 -510 -400 -300 -200 -100 -50 0 50 100 200 300 400 510 512 520 600 800 1000 1020)
  "The powers of two a problem is multiplied by. All even, so that a square root of the factor is a power of two.")

(define-condition sweep-premise (error)
  ((text :initarg :text :reader sweep-premise-text))
  (:report (lambda (condition stream) (write-string (sweep-premise-text condition) stream))))

(defun sweep-finite-p (value)
  (or (not (floatp value)) (and (= value value) (not (sb-ext:float-infinity-p value)))))

(defun sweep-number (value k)
  "VALUE times 2^K. An input the format cannot hold as a normal number is no input: SWEEP-PREMISE is signalled."
  (let ((scaled (scale-float (coerce value 'double-float) k)))
    (unless (and (sweep-finite-p scaled)
                 (or (zerop value) (>= (abs scaled) least-positive-normalized-double-float)))
      (error 'sweep-premise :text (format nil "~S times 2^~D is not a normal number." value k)))
    scaled))

(defun sweep-array (array k)
  (map '(simple-array double-float (*)) (lambda (value) (sweep-number value k)) array))

(defun sweep-rows (rows)
  "A dense matrix from its rows, of any shape."
  (let* ((m (length rows)) (n (length (first rows))) (a (make-array (list m n))))
    (dotimes (i m) (dotimes (j n) (setf (aref a i j) (nth j (nth i rows)))))
    a))

(defun sweep-stored (dense k)
  "The nonzero entries of a dense matrix, each times 2^K, as a checked CSC matrix."
  (let ((m (array-dimension dense 0)) (n (array-dimension dense 1)))
    (csc-from-entries m n (loop :for j :below n :append
                            (loop :for i :below m
                                  :unless (zerop (aref dense i j))
                                    :collect (list i j (sweep-number (aref dense i j) k)))))))

(defun sweep-column-major (dense k)
  (let* ((m (array-dimension dense 0)) (n (array-dimension dense 1))
         (a (make-array (* m n) :element-type 'double-float)))
    (dotimes (j n) (dotimes (i m) (setf (aref a (+ i (* j m))) (sweep-number (aref dense i j) k))))
    a))

(defun sweep-lower (dense k)
  "The lower triangle of a symmetric dense matrix as a pattern's arrays, the values times 2^K."
  (multiple-value-bind (n pointers rows values) (symmetric-lower-pattern dense)
    (values n pointers rows (sweep-array values k))))

(defmacro sweep-need (form constructor)
  "Run FORM, a call that returns a status, and refuse the run unless the status is CONSTRUCTOR."
  (let ((status (gensym "STATUS")))
    `(let ((,status ,form))
       (unless (eq ,status (coalton:coalton ,constructor))
         (error "~(~A~) returned ~A" ',(first form) ,status)))))

(defun sweep-out (label value power &key (kind :law) bound)
  "One output of a run: its right value at scale 2^k is its value at scale 1 times 2^(POWER k). KIND :law for a
value, :bound for an upper bound, which may come out larger, and :within for a quantity with no value of its own
that must lie in [0, BOUND 2^(POWER k)]."
  (list label (if (arrayp value) (copy-seq value) value) power kind bound))

;;; ---- The class of one run ------------------------------------------------------------------

(defparameter *sweep-order* '(:exact :close :looser :lost :wrong))

(defun sweep-in-range-p (expected)
  (or (zerop expected)
      (<= (rational least-positive-normalized-double-float) (abs expected) (rational most-positive-double-float))))

(defun sweep-element (actual reference largest shift kind bound tolerance)
  "The class of one number against its value at scale 1. LARGEST is the largest magnitude at scale 1 in its array."
  (cond
    ((eq kind :within)
     (if (and (sweep-finite-p actual) (<= 0 (rational actual) (* (rational bound) (expt 2 shift)))) :exact :wrong))
    ((not (sweep-finite-p actual)) (if (and (eq kind :bound) (= actual actual)) :looser :wrong))
    (t
     (let ((expected (* (rational reference) (expt 2 shift))) (value (rational actual)))
       (cond ((= value expected) :exact)
             ;; A bound whose scaled value the format cannot hold is still a bound if it is above it.
             ((and (eq kind :bound) (not (sweep-in-range-p expected))) (if (> value expected) :looser :wrong))
             ((not (sweep-in-range-p expected)) :lost)
             ((<= (abs (- value expected)) (* tolerance largest (expt 2 shift))) :close)
             ((and (eq kind :bound) (> value expected)) :looser)
             (t :wrong))))))

(defun sweep-values (value)
  (if (arrayp value) (coerce value 'list) (list value)))

(defun sweep-holdable-p (reference k)
  "Whether every value of the outputs at scale 1 has a scaled value the format can hold at 2^K."
  (loop :for (nil value power kind) :in reference
        :always (or (not (eq kind :law))
                    (loop :for r :in (sweep-values value)
                          :always (and (sweep-finite-p r) (sweep-in-range-p (* (rational r) (expt 2 (* power k)))))))))

(defun sweep-class (run reference k tolerance)
  "The class of the run at scale 2^K against the outputs at scale 1, and for a run that is not :exact what decided
it: the output and the number, or the refusal's own words."
  (let ((outputs (handler-case (funcall run k)
                   (sweep-premise (condition)
                     (return-from sweep-class (values :no-input (princ-to-string condition))))
                   (error (condition)
                     (return-from sweep-class
                       (values (if (sweep-holdable-p reference k) :limit :refused)
                               (substitute #\Space #\Newline (princ-to-string condition)))))))
        (class :exact) (decided nil))
    (assert (= (length outputs) (length reference)) ()
            "At 2^~D the run has ~D outputs and at scale 1 it has ~D." k (length outputs) (length reference))
    (loop :for (label value power kind bound) :in outputs
          :for (nil reference-value) :in reference
          :do (let* ((shift (* power k)) (actual (sweep-values value)) (expected (sweep-values reference-value))
                     (largest (loop :for r :in expected :maximize (if (sweep-finite-p r) (abs (rational r)) 0))))
                (assert (integerp shift) () "~A: 2^(~S * ~D) is not a power of two." label power k)
                (flet ((take (found detail)
                         (when (> (position found *sweep-order*) (position class *sweep-order*))
                           (setf class found decided detail))))
                  (if (/= (length actual) (length expected))
                      (take :wrong (list label :length (length actual) :at-scale-1 (length expected)))
                      (loop :for a :in actual :for r :in expected :for position :from 0
                            :do (take (sweep-element a r largest shift kind bound tolerance)
                                      (list label position a :at-scale-1 r)))))))
    (values class decided)))

(defun sweep-segments (classes)
  "A list of (k class), in order, as the ranges each class holds over: (first-k last-k class)."
  (let ((segments '()))
    (loop :for (k class) :in classes
          :do (if (and segments (eq (third (first segments)) class))
                  (setf (second (first segments)) k)
                  (push (list k k class) segments)))
    (nreverse segments)))

;;; ---- The members ---------------------------------------------------------------------------

(defmacro sweep-entry (name covers tolerance declared run)
  "An entry: its name, the calls it answers for, the tolerance of :close, the declared ranges, and the run, a
function of the exponent that returns the outputs. The run's source is kept, to check the calls against."
  `(list ,name ',covers ,tolerance ',declared ,run ',run))

(defun sweep-arrow () (sweep-rows '((4 1 1 1) (1 3 0 0) (1 0 3 0) (1 0 0 3))))
(defun sweep-unsymmetric () (sweep-rows '((0 2 0 1) (3 0 1 0) (0 1 4 0) (1 0 2 5))))
;;; The two indefinite matrices begin with a two-by-two pivot, which each entry asserts, and put it
;;; at the ends of the sweep. In the first the block is (6, 14; 14, -6): at 2^1020 its off-diagonal
;;; entry times the block's dimensionless determinant, -1.18, is too large to hold. In the second it
;;; is (1/4, 7/16; 7/16, -1/4): at 2^-1020 the products that eliminate with it are subnormal. Those
;;; are the two magnitudes the block's arithmetic is moved away from.
(defun sweep-indefinite () (sweep-rows '((6 1 14) (1 2 3) (14 3 -6))))
(defun sweep-indefinite-small () (sweep-rows '((1/4 1/4 7/16) (1/4 1/2 5/16) (7/16 5/16 -1/4))))

(defun sweep-block-first (factor)
  "The premise of the two: the factorization opened with a two-by-two pivot."
  (unless (= 2 (aref (il:factor-block-sizes factor) 0))
    (error 'sweep-premise :text (format nil "The first pivot is not a two-by-two block: ~S."
                                        (il:factor-block-sizes factor)))))

(defun sweep-tall () (sweep-rows '((1 5 2) (4 1 7) (3 2 1) (2 6 3))))
;;; Diagonal entries of opposite sign whose difference, at the top of the sweep, is too large to hold.
(defun sweep-symmetric () (sweep-rows '((12 1 2) (1 -11 1) (2 1 5))))
(defun sweep-dependent () (sweep-rows '((1 2 0) (0 1 3) (2 0 1) (1 3 3))))
(defun sweep-repeated () (sweep-rows '((2 2 0) (2 2 0) (0 0 3))))
(defun sweep-twice () (sweep-rows '((1 2 0) (0 1 3) (2 0 1) (1 2 0))))
(defun sweep-wide () (sweep-rows '((1 2 0 3) (0 1 4 1))))

(defun sweep-prefix (array count) (subseq array 0 count))

(defun sweep-roster ()
  (let ((b4 (doubles 1d0 -2d0 3d0 5d0)) (b3 (doubles 1d0 -2d0 3d0)) (small (scale-float 1d0 -40)))
    (list
     (sweep-entry "norms of a vector" (numeric:norm2 numeric:norm-inf) 1d-15
         ((-1020 1020 :exact))
       (lambda (k)
         (let ((v (sweep-array (doubles 3d0 -4d0 12d0) k)))
           (list (sweep-out "norm2" (numeric:norm2 v) 1) (sweep-out "norm-inf" (numeric:norm-inf v) 1)))))

     (sweep-entry "a count as a number" (coalton-linalg/internal/float:count->float) 0d0
         ((-1020 1020 :exact))
       (lambda (k)
         (declare (ignore k))
         (list (sweep-out "the count" (coalton-linalg/internal/float:count->float 12345) 0))))

     (sweep-entry "the power of two at a magnitude" (coalton-linalg/internal/float:binade!) 0d0
         ((-1020 1020 :exact))
       (lambda (k)
         (let ((pair (zeros 2)))
           (coalton-linalg/internal/float:binade! (sweep-number 5d0 k) pair 0)
           (list (sweep-out "the power" (aref pair 0) 1) (sweep-out "its reciprocal" (aref pair 1) -1)))))

     (sweep-entry "a stored matrix and its products, the matrix scaled"
         (csc:from-csc csc:from-coo csc:entries csc:matvec csc:transpose-matvec csc:matvec! csc:transpose-matvec!)
         1d-15
         ((-1020 1000 :exact) (1020 1020 :refused))
       (lambda (k)
         (let* ((copied (csc:from-csc 3 2 (indices 0 2 3) (indices 0 2 1) (sweep-array (doubles 5d0 7d0 3d0) k)))
                (summed (csc:from-coo 2 2 (indices 0 0 1) (indices 0 0 1) (sweep-array (doubles 1.5d0 2.5d0 3d0) k)
                                      (coalton:coalton coalton-linalg/sparse/csc:SumDuplicates)))
                (matrix (sweep-stored (sweep-tall) k)) (x (doubles 1d0 -2d0 3d0)) (y (doubles 1d0 2d0 -1d0 3d0))
                (product (zeros 4)) (transposed (zeros 3)))
           (csc:matvec! matrix x product)
           (csc:transpose-matvec! matrix y transposed)
           (list (sweep-out "from-csc" (csc:entries copied) 1) (sweep-out "from-coo" (csc:entries summed) 1)
                 (sweep-out "matvec" (csc:matvec matrix x) 1) (sweep-out "matvec!" product 1)
                 (sweep-out "transpose-matvec" (csc:transpose-matvec matrix y) 1)
                 (sweep-out "transpose-matvec!" transposed 1)))))

     (sweep-entry "the products of a stored matrix, the vector scaled"
         (csc:matvec csc:transpose-matvec csc:matvec! csc:transpose-matvec!) 1d-15
         ((-1020 1000 :exact) (1020 1020 :refused))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-tall) 0)) (x (sweep-array (doubles 1d0 -2d0 3d0) k))
                (y (sweep-array (doubles 1d0 2d0 -1d0 3d0) k)) (product (zeros 4)) (transposed (zeros 3)))
           (csc:matvec! matrix x product)
           (csc:transpose-matvec! matrix y transposed)
           (list (sweep-out "matvec" (csc:matvec matrix x) 1) (sweep-out "matvec!" product 1)
                 (sweep-out "transpose-matvec" (csc:transpose-matvec matrix y) 1)
                 (sweep-out "transpose-matvec!" transposed 1)))))

     (sweep-entry "the Gram matrix" (gram:refill! gram:refill-csc!) 1d-15
         ((-1020 -600 :lost) (-520 400 :exact) (510 1020 :refused))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-tall) k)) (plan (gram:analyze-ata matrix))
                (raw (zeros (gram:entry-count plan))) (checked (zeros (gram:entry-count plan))))
           (gram:refill! plan (csc:entries matrix) raw)
           (gram:refill-csc! plan matrix checked)
           (list (sweep-out "refill!" raw 2) (sweep-out "refill-csc!" checked 2)))))

     (sweep-entry "the augmented matrix" (aug:refill! aug:refill-csc!) 1d-15
         ((-1020 1020 :exact))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-wide) k)) (plan (aug:analyze matrix))
                (raw (zeros (aug:entry-count plan))) (checked (zeros (aug:entry-count plan))))
           (aug:refill! plan (csc:entries matrix) (sweep-number 1d0 k) (sweep-number 0.25d0 k) raw)
           (aug:refill-csc! plan matrix (sweep-number 1d0 k) (sweep-number 0.25d0 k) checked)
           (list (sweep-out "refill!" raw 1) (sweep-out "refill-csc!" checked 1)))))

     (sweep-entry "scaling by rows and columns"
         (sc:row-maxima! sc:column-maxima! sc:equilibrate! sc:scale! sc:scale-vector!) 1d-15
         ((-1020 1020 :exact))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-tall) k)) (values (csc:entries matrix))
                (row-maxima (zeros 4)) (column-maxima (zeros 3)) (rows (zeros 4)) (columns (zeros 3))
                (scaled (zeros (length values))) (vector (zeros 3)))
           (sc:row-maxima! matrix values row-maxima)
           (sc:column-maxima! matrix values column-maxima)
           (sc:equilibrate! matrix values rows columns)
           (sc:scale! matrix values rows columns scaled)
           (sc:scale-vector! (doubles 0.5d0 2d0 4d0) (sweep-array b3 k) vector)
           (list (sweep-out "row-maxima!" row-maxima 1) (sweep-out "column-maxima!" column-maxima 1)
                 (sweep-out "the matrix scaled by equilibrate!" scaled 0) (sweep-out "scale-vector!" vector 1)))))

     (sweep-entry "scaling a symmetric matrix" (sc:equilibrate-symmetric! sc:scale-symmetric!) 1d-15
         ((-1020 1020 :exact))
       (lambda (k)
         (let* ((dense (sweep-arrow))
                (matrix (csc-from-entries 4 4 (loop :for j :below 4 :append
                                                (loop :for i :from j :below 4
                                                      :unless (zerop (aref dense i j))
                                                        :collect (list i j (sweep-number (aref dense i j) k))))))
                (values (csc:entries matrix)) (scales (zeros 4)) (scaled (zeros (length values))))
           (sweep-need (sc:equilibrate-symmetric! matrix values scales) coalton-linalg/sparse/scaling:DiagonalScaled)
           (sc:scale-symmetric! matrix values scales scaled)
           (list (sweep-out "the matrix scaled" scaled 0)))))

     (sweep-entry "the sparse factorization with positive pivots, its solve and its condition"
         (ldlt:factor! ldlt:solve! ldlt:factor-values ldlt:factor-diagonal ldlt:estimate-inverse-norm1!
          ldlt:inverse-norm1-estimate ldlt:shifted-norm1 ldlt:reciprocal-condition1)
         1d-15
         ((-1020 1000 :exact) (1020 1020 :close))
       (lambda (k)
         (multiple-value-bind (n pointers rows values) (sweep-lower (sweep-arrow) k)
           (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers rows))) (factor (ldlt:symbolic-factor symbolic))
                  (shift (sweep-number 0.5d0 k)) (x (zeros 4)))
             (sweep-need (ldlt:factor! symbolic values shift) coalton-linalg/sparse/spd-ldlt:FactorSuccess)
             (sweep-need (ldlt:solve! factor (sweep-array b4 k) x) coalton-linalg/sparse/spd-ldlt:SolveSuccess)
             (sweep-need (ldlt:estimate-inverse-norm1! factor) coalton-linalg/sparse/spd-ldlt:ConditionSuccess)
             (let ((norm (ldlt:shifted-norm1 symbolic values shift)) (inverse (ldlt:inverse-norm1-estimate factor)))
               (list (sweep-out "L" (ldlt:factor-values factor) 0) (sweep-out "D" (ldlt:factor-diagonal factor) 1)
                     (sweep-out "x" x 0) (sweep-out "the norm" norm 1) (sweep-out "the inverse's norm" inverse -1)
                     (sweep-out "the reciprocal condition" (ldlt:reciprocal-condition1 norm inverse) 0)))))))

     (sweep-entry "the Cholesky form and its solve"
         (ldlt:factor-cholesky! ldlt:solve-cholesky! ldlt:factor-values ldlt:factor-diagonal) 1d-15
         ((-1020 1020 :exact))
       (lambda (k)
         (multiple-value-bind (n pointers rows values) (sweep-lower (sweep-arrow) k)
           (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers rows))) (factor (ldlt:symbolic-factor symbolic))
                  (x (zeros 4)))
             (sweep-need (ldlt:factor-cholesky! symbolic values (sweep-number 0.5d0 k))
                         coalton-linalg/sparse/spd-ldlt:FactorSuccess)
             (sweep-need (ldlt:solve-cholesky! factor (sweep-array b4 k) x) coalton-linalg/sparse/spd-ldlt:SolveSuccess)
             (list (sweep-out "G below its diagonal" (ldlt:factor-values factor) 1/2)
                   (sweep-out "G's diagonal" (ldlt:factor-diagonal factor) 1/2) (sweep-out "x" x 0))))))

     (sweep-entry "the factorization with pivots of either sign, its growth and its bound"
         (ldlt:factor-signed! ldlt:solve! ldlt:factor-values ldlt:factor-diagonal ldlt:factor-growth
          ldlt:factor-backward-bound)
         1d-15
         ((-1020 -1000 :looser) (-800 1000 :exact) (1020 1020 :refused))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-wide) k)) (plan (aug:analyze matrix))
                (symbolic (ldlt:analyze-quasidefinite
                           (ldlt:pattern 6 (aug:column-pointers plan) (aug:row-indices plan)) 4))
                (factor (ldlt:symbolic-factor symbolic)) (values (zeros (aug:entry-count plan))) (x (zeros 6)))
           (aug:refill! plan (csc:entries matrix) (sweep-number 1d0 k) (sweep-number 0.25d0 k) values)
           (sweep-need (ldlt:factor-signed! symbolic values 0d0) coalton-linalg/sparse/spd-ldlt:SignedSuccess)
           (sweep-need (ldlt:solve! factor (sweep-array (doubles 1d0 -2d0 3d0 5d0 2d0 -1d0) k) x)
                       coalton-linalg/sparse/spd-ldlt:SolveSuccess)
           (list (sweep-out "L" (ldlt:factor-values factor) 0) (sweep-out "D" (ldlt:factor-diagonal factor) 1)
                 (sweep-out "x" x 0) (sweep-out "the growth" (ldlt:factor-growth factor) 0)
                 (sweep-out "the backward bound" (ldlt:factor-backward-bound factor) 1 :kind :bound)))))

     (sweep-entry "the factorization that leaves pivots out, its solve and its null vector"
         (ldlt:factor-skipping! ldlt:solve-skipping! ldlt:skipped-null-vector! ldlt:factor-values ldlt:factor-diagonal)
         1d-15
         ((-1020 1020 :exact))
       (lambda (k)
         (multiple-value-bind (n pointers rows values) (sweep-lower (sweep-repeated) k)
           (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers rows))) (factor (ldlt:symbolic-factor symbolic))
                  (x (zeros 3)) (v (zeros 3)))
             (sweep-need (ldlt:factor-skipping! symbolic values 0d0 (sweep-number 0.5d0 k))
                         coalton-linalg/sparse/spd-ldlt:SkipSuccess)
             (unless (= 1 (ldlt:skipped-count factor))
               (error 'sweep-premise :text (format nil "~D pivots were left out, not 1." (ldlt:skipped-count factor))))
             (sweep-need (ldlt:solve-skipping! factor (sweep-array b3 k) x) coalton-linalg/sparse/spd-ldlt:SolveSuccess)
             (let ((index (aref (ldlt:factor-permutation factor) (position 0d0 (ldlt:factor-diagonal factor)))))
               (sweep-need (ldlt:skipped-null-vector! symbolic index v) coalton-linalg/sparse/spd-ldlt:NullVectorSuccess))
             (list (sweep-out "L" (ldlt:factor-values factor) 0) (sweep-out "D" (ldlt:factor-diagonal factor) 1)
                   (sweep-out "x" x 0) (sweep-out "the null vector" v 0))))))

     (sweep-entry "the dense symmetric factorization of small entries"
         (il:factor! il:solve! il:factor-lower il:factor-diagonal il:factor-subdiagonal) 1d-15
         ((-1020 1020 :exact))
       (lambda (k)
         (let ((factor (il:make-factor 3)) (x (zeros 3)))
           (sweep-need (il:factor! factor (sweep-column-major (sweep-indefinite-small) k) 0d0)
                       coalton-linalg/dense/indefinite-ldlt:FactorSuccess)
           (sweep-block-first factor)
           (sweep-need (il:solve! factor (sweep-array b3 k) x) coalton-linalg/dense/indefinite-ldlt:SolveSuccess)
           (list (sweep-out "L" (il:factor-lower factor) 0) (sweep-out "D" (il:factor-diagonal factor) 1)
                 (sweep-out "D below its diagonal" (il:factor-subdiagonal factor) 1) (sweep-out "x" x 0)
                 (sweep-out "the blocks" (il:factor-block-sizes factor) 0)))))

     (sweep-entry "the dense symmetric factorization and its solve"
         (il:factor! il:solve! il:factor-lower il:factor-diagonal il:factor-subdiagonal) 1d-15
         ((-1020 -1020 :close) (-1000 1020 :exact))
       (lambda (k)
         (let ((factor (il:make-factor 3)) (x (zeros 3)))
           (sweep-need (il:factor! factor (sweep-column-major (sweep-indefinite) k) 0d0)
                       coalton-linalg/dense/indefinite-ldlt:FactorSuccess)
           (sweep-block-first factor)
           (sweep-need (il:solve! factor (sweep-array b3 k) x) coalton-linalg/dense/indefinite-ldlt:SolveSuccess)
           (list (sweep-out "L" (il:factor-lower factor) 0) (sweep-out "D" (il:factor-diagonal factor) 1)
                 (sweep-out "D below its diagonal" (il:factor-subdiagonal factor) 1) (sweep-out "x" x 0)
                 (sweep-out "the blocks" (il:factor-block-sizes factor) 0)))))

     (sweep-entry "the orthogonal factorization, its ranks, solves and vectors"
         (og:factor! og:pivot-magnitude og:rank-above og:set-rank! og:solve! og:solve-residual-norm
          og:nullspace-vector! og:apply-q! og:apply-qt! og:factor-work og:factor-tau)
         1d-13
         ((-1020 -1020 :close) (-1000 1020 :exact))
       (lambda (k)
         (let ((decomposition (og:make-decomposition 4 3)) (x (zeros 3)) (null (zeros 3))
               (q (doubles 1d0 2d0 -1d0 3d0)) (qt (doubles 1d0 2d0 -1d0 3d0)))
           (sweep-need (og:factor! decomposition (sweep-column-major (sweep-tall) k))
                       coalton-linalg/dense/orthogonal:FactorSuccess)
           (let* ((work (copy-seq (og:factor-work decomposition))) (tau (copy-seq (og:factor-tau decomposition)))
                  (triangle (loop :for j :below 3 :append (loop :for i :to j :collect (aref work (+ i (* j 4))))))
                  (tails (loop :for j :below 3 :append
                           (loop :for i :from (1+ j) :below 4 :collect (aref work (+ i (* j 4))))))
                  (pivots (loop :for i :below 3 :collect (og:pivot-magnitude decomposition i)))
                  (rank (og:rank-above decomposition (sweep-number 0.5d0 k))))
             (sweep-need (og:apply-q! decomposition q) coalton-linalg/dense/orthogonal:SolveSuccess)
             (sweep-need (og:apply-qt! decomposition qt) coalton-linalg/dense/orthogonal:SolveSuccess)
             (sweep-need (og:set-rank! decomposition 3) coalton-linalg/dense/orthogonal:RankSuccess)
             (sweep-need (og:solve! decomposition (sweep-array b4 k) x) coalton-linalg/dense/orthogonal:SolveSuccess)
             (let ((residual (og:solve-residual-norm decomposition)))
               (sweep-need (og:set-rank! decomposition 2) coalton-linalg/dense/orthogonal:RankSuccess)
               (sweep-need (og:nullspace-vector! decomposition 0 null) coalton-linalg/dense/orthogonal:SolveSuccess)
               (list (sweep-out "R" (coerce triangle 'vector) 1) (sweep-out "the reflectors' tails" (coerce tails 'vector) 0)
                     (sweep-out "tau" tau 0) (sweep-out "the pivots' magnitudes" (coerce pivots 'vector) 1)
                     (sweep-out "the rank above a threshold" rank 0) (sweep-out "Q v" q 0) (sweep-out "Q^T v" qt 0)
                     (sweep-out "x" x 0) (sweep-out "the residual" residual 1)
                     (sweep-out "a nullspace vector" null 0)))))))

     (sweep-entry "the symmetric eigenvalues" (jc:eigen! jc:eigen-with-sweep-limit! jc:eigenvalues jc:eigenvectors)
         1d-13
         ((-1020 1020 :exact))
       (lambda (k)
         (let ((plain (jc:make-symmetric-eigen 3)) (limited (jc:make-symmetric-eigen 3)))
           (sweep-need (jc:eigen! plain (sweep-column-major (sweep-symmetric) k)) coalton-linalg/dense/jacobi:JacobiSuccess)
           (sweep-need (jc:eigen-with-sweep-limit! limited (sweep-column-major (sweep-symmetric) k) 30)
                       coalton-linalg/dense/jacobi:JacobiSuccess)
           (list (sweep-out "the eigenvalues" (jc:eigenvalues plain) 1) (sweep-out "the eigenvectors" (jc:eigenvectors plain) 0)
                 (sweep-out "the eigenvalues, with a limit" (jc:eigenvalues limited) 1)
                 (sweep-out "the eigenvectors, with a limit" (jc:eigenvectors limited) 0)))))

     (sweep-entry "the singular values"
         (jc:svd! jc:svd-with-sweep-limit! jc:singular-values jc:left-vectors jc:right-vectors jc:rank-above)
         1d-13
         ((-1020 1020 :exact))
       (lambda (k)
         (let ((plain (jc:make-svd 4 3)) (limited (jc:make-svd 4 3)))
           (sweep-need (jc:svd! plain (sweep-column-major (sweep-tall) k)) coalton-linalg/dense/jacobi:JacobiSuccess)
           (sweep-need (jc:svd-with-sweep-limit! limited (sweep-column-major (sweep-tall) k) 30)
                       coalton-linalg/dense/jacobi:JacobiSuccess)
           (list (sweep-out "the singular values" (jc:singular-values plain) 1)
                 (sweep-out "the left vectors" (jc:left-vectors plain) 0)
                 (sweep-out "the right vectors" (jc:right-vectors plain) 0)
                 (sweep-out "the rank above a threshold" (jc:rank-above plain (sweep-number 0.5d0 k)) 0)
                 (sweep-out "the singular values, with a limit" (jc:singular-values limited) 1)))))

     (sweep-entry "the sparse LU factorization, its solves, its growth and its bound"
         (lu:factor! lu:solve! lu:solve-transpose! lu:lower-values lu:upper-values lu:upper-diagonal lu:lu-growth
          lu:lu-backward-bound)
         1d-15
         ((-1020 -1000 :looser) (-800 1000 :exact) (1020 1020 :looser))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-unsymmetric) k)) (handle (lu:make-lu matrix)) (x (zeros 4)) (y (zeros 4)))
           (sweep-need (lu:factor! handle matrix 1d0) coalton-linalg/sparse/lu:LuSuccess)
           (sweep-need (lu:solve! handle (sweep-array b4 k) x) coalton-linalg/sparse/lu:LuSolveSuccess)
           (sweep-need (lu:solve-transpose! handle (sweep-array b4 k) y) coalton-linalg/sparse/lu:LuSolveSuccess)
           (list (sweep-out "L" (sweep-prefix (lu:lower-values handle) (lu:lower-entry-count handle)) 0)
                 (sweep-out "U" (sweep-prefix (lu:upper-values handle) (lu:upper-entry-count handle)) 1)
                 (sweep-out "U's diagonal" (lu:upper-diagonal handle) 1) (sweep-out "x" x 0)
                 (sweep-out "the transposed x" y 0) (sweep-out "the growth" (lu:lu-growth handle) 0)
                 (sweep-out "the backward bound" (lu:lu-backward-bound handle) 1 :kind :bound)))))

     (sweep-entry "conjugate gradients, the matrix and the right-hand side scaled together"
         (cg:cg! cg:residual-norm cg:recurrence-norms) 1d-9
         ((-1020 -1000 :no-input) (-800 1020 :exact))
       (lambda (k)
         (let ((apply (rf-caller:matrix-operator (sweep-stored (sweep-arrow) k))) (workspace (cg:make-cg 4 50))
               (x (zeros 4)))
           (sweep-need (cg:cg! workspace apply (sweep-array b4 k) x (sweep-number small k))
                       coalton-linalg/conjugate-gradient:CgConverged)
           (list (sweep-out "x" x 0) (sweep-out "the steps" (cg:step-count workspace) 0)
                 (sweep-out "the residual" (cg:residual-norm workspace) 1 :kind :within :bound small)
                 (sweep-out "the first running norm" (aref (cg:recurrence-norms workspace) 0) 1)))))

     (sweep-entry "conjugate gradients, the right-hand side scaled"
         (cg:cg! cg:residual-norm cg:recurrence-norms) 1d-9
         ((-1020 -1000 :no-input) (-800 1020 :exact))
       (lambda (k)
         (let ((apply (rf-caller:matrix-operator (sweep-stored (sweep-arrow) 0))) (workspace (cg:make-cg 4 50))
               (x (zeros 4)))
           (sweep-need (cg:cg! workspace apply (sweep-array b4 k) x (sweep-number small k))
                       coalton-linalg/conjugate-gradient:CgConverged)
           (list (sweep-out "x" x 1) (sweep-out "the steps" (cg:step-count workspace) 0)
                 (sweep-out "the residual" (cg:residual-norm workspace) 1 :kind :within :bound small)
                 (sweep-out "the first running norm" (aref (cg:recurrence-norms workspace) 0) 1)))))

     (sweep-entry "preconditioned conjugate gradients, the matrix, its preconditioner and the right-hand side scaled"
         (cg:pcg! cg:residual-norm cg:recurrence-norms) 1d-9
         ((-1020 -1000 :no-input) (-800 1000 :exact) (1020 1020 :close))
       (lambda (k)
         (let ((apply (rf-caller:matrix-operator (sweep-stored (sweep-arrow) k))) (workspace (cg:make-cg 4 50))
               (x (zeros 4)))
           (sweep-need (cg:pcg! workspace apply
                                (cg-caller:diagonal-preconditioner (sweep-array (doubles 4d0 3d0 3d0 3d0) k))
                                (sweep-array b4 k) x (sweep-number small k))
                       coalton-linalg/conjugate-gradient:CgConverged)
           (list (sweep-out "x" x 0) (sweep-out "the steps" (cg:step-count workspace) 0)
                 (sweep-out "the residual" (cg:residual-norm workspace) 1 :kind :within :bound small)
                 (sweep-out "the first running norm" (aref (cg:recurrence-norms workspace) 0) 1)))))

     (sweep-entry "the damped step under the columns' own damping, the matrix scaled"
         (dm:set-column-damping! dm:factor! dm:step!) 1d-15
         ((-1020 1020 :exact))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-tall) k)) (plan (dm:analyze matrix)) (step (zeros 3)))
           (dm:set-matrix! plan matrix)
           (sweep-need (dm:set-column-damping! plan 0d0) coalton-linalg/sparse/damped:DampingSet)
           (sweep-need (dm:factor! plan 0.25d0) coalton-linalg/sparse/damped:DampedSuccess)
           (sweep-need (dm:step! plan b4 step) coalton-linalg/sparse/damped:StepSuccess)
           (list (sweep-out "the step" step -1)))))

     (sweep-entry "the damped step under a given damping, the matrix and the residual scaled"
         (dm:set-damping! dm:factor! dm:step!) 1d-15
         ((-1020 -512 :no-input) (-510 510 :exact) (512 1020 :no-input))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-tall) k)) (plan (dm:analyze matrix)) (step (zeros 3)))
           (dm:set-matrix! plan matrix)
           (sweep-need (dm:set-damping! plan (sweep-array (doubles 2d0 0.5d0 3d0) (* 2 k)))
                       coalton-linalg/sparse/damped:DampingSet)
           (sweep-need (dm:factor! plan 0.25d0) coalton-linalg/sparse/damped:DampedSuccess)
           (sweep-need (dm:step! plan (sweep-array b4 k) step) coalton-linalg/sparse/damped:StepSuccess)
           (list (sweep-out "the step" step 0)))))

     (sweep-entry "the damped step, the residual scaled" (dm:factor! dm:step!) 1d-15
         ((-1020 1020 :exact))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-tall) 0)) (plan (dm:analyze matrix)) (step (zeros 3)))
           (dm:set-matrix! plan matrix)
           (sweep-need (dm:set-column-damping! plan 0d0) coalton-linalg/sparse/damped:DampingSet)
           (sweep-need (dm:factor! plan 0.25d0) coalton-linalg/sparse/damped:DampedSuccess)
           (sweep-need (dm:step! plan (sweep-array b4 k) step) coalton-linalg/sparse/damped:StepSuccess)
           (list (sweep-out "the step" step 1)))))

     (sweep-entry "kept rows, their vectors and bounds, proved again, and the certificate"
         (rd:find-kept-rows! rd:reprove! rd:residual-bound rd:gap-bound rd:multiplier-values rd:certify!
          rd:proved-lower-bound)
         1d-15
         ((-1020 -1020 :no-input) (-1000 -1000 :looser) (-800 1020 :exact))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-dependent) k)) (kept (rd:analyze matrix))
                (floor (sweep-number (scale-float 1d0 -20) k)))
           (sweep-need (rd:find-kept-rows! kept matrix floor 1d-12) coalton-linalg/sparse/rank-deficient:FindFound)
           (when (and (zerop k) (/= 1 (rd:left-out-count kept)))
             (error 'sweep-premise :text (format nil "~D rows were left out, not 1." (rd:left-out-count kept))))
           (let ((found-bound (rd:residual-bound kept 0)))
             (sweep-need (rd:reprove! kept matrix) coalton-linalg/sparse/rank-deficient:Reproved)
             (sweep-need (rd:certify! kept floor floor) coalton-linalg/sparse/rank-deficient:Certified)
             (list (sweep-out "how many rows are left out" (rd:left-out-count kept) 0)
                   (sweep-out "the vector" (sweep-prefix (rd:multiplier-values kept) (rd:multiplier-count kept 0)) 0)
                   (sweep-out "the bound as found" found-bound 1 :kind :bound)
                   (sweep-out "the bound proved again" (rd:residual-bound kept 0) 1 :kind :bound)
                   (sweep-out "the gap bound" (rd:gap-bound kept) 1 :kind :bound)
                   (sweep-out "the lower bound proved" (rd:proved-lower-bound kept) 1))))))

     (sweep-entry "the step on certified kept rows, the matrix and the residual scaled together"
         (rd:step! rd:step-solve-bound rd:step-bound rd:solve-growth rd:solve-backward-bound rd:predicted-rate)
         1d-15
         ((-1020 -1020 :no-input) (-1000 1020 :exact))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-dependent) k)) (kept (rd:analyze matrix))
                (floor (sweep-number (scale-float 1d0 -20) k)) (step (zeros 3)))
           (sweep-need (rd:find-kept-rows! kept matrix floor 1d-12) coalton-linalg/sparse/rank-deficient:FindFound)
           (sweep-need (rd:certify! kept floor floor) coalton-linalg/sparse/rank-deficient:Certified)
           (sweep-need (rd:step! kept (sweep-array b4 k) step) coalton-linalg/sparse/rank-deficient:StepTaken)
           (list (sweep-out "the step" step 0)
                 (sweep-out "the solve's bound" (rd:step-solve-bound kept) 0 :kind :bound)
                 (sweep-out "the whole bound" (rd:step-bound kept) 0 :kind :bound)
                 (sweep-out "the solve factor's growth" (rd:solve-growth kept) 0)
                 (sweep-out "its backward bound" (rd:solve-backward-bound kept) 0 :kind :bound)
                 (sweep-out "the predicted rate" (rd:predicted-rate kept) 0 :kind :bound)))))

     (sweep-entry "the step on certified kept rows, the residual scaled" (rd:step! rd:step-solve-bound rd:step-bound)
         1d-15
         ((-1020 -1000 :looser) (-800 1020 :exact))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-dependent) 0)) (kept (rd:analyze matrix)) (step (zeros 3)))
           (sweep-need (rd:find-kept-rows! kept matrix 1d-6 1d-12) coalton-linalg/sparse/rank-deficient:FindFound)
           (sweep-need (rd:certify! kept 1d-6 1d-6) coalton-linalg/sparse/rank-deficient:Certified)
           (sweep-need (rd:step! kept (sweep-array b4 k) step) coalton-linalg/sparse/rank-deficient:StepTaken)
           (list (sweep-out "the step" step 1)
                 (sweep-out "the solve's bound" (rd:step-solve-bound kept) 1 :kind :bound)
                 (sweep-out "the whole bound" (rd:step-bound kept) 1 :kind :bound)))))

     (sweep-entry "bounds on the largest singular value" (rd:largest-singular-bounds!) 1d-15
         ((-1020 1020 :exact))
       (lambda (k)
         (let ((out (zeros 2)))
           (rd:largest-singular-bounds! (sweep-stored (sweep-tall) k) 20 (zeros 3) (zeros 4) out)
           (list (sweep-out "the lower bound" (aref out 0) 1) (sweep-out "the upper bound" (aref out 1) 1 :kind :bound)))))

     (sweep-entry "refinement" (rf:refine! rf:residual-norms) 1d-9
         ((-1020 -1000 :no-input) (-800 1020 :exact))
       (lambda (k)
         (multiple-value-bind (n pointers rows values) (sweep-lower (sweep-arrow) k)
           (let* ((symbolic (ldlt:analyze (ldlt:pattern n pointers rows))) (factor (ldlt:symbolic-factor symbolic))
                  (refiner (rf:make-refiner 4 20)) (x (zeros 4)))
             ;; The solver is the factor of the matrix less a twentieth of the identity, so that
             ;; corrections are needed.
             (sweep-need (ldlt:factor! symbolic values (sweep-number 0.05d0 k)) coalton-linalg/sparse/spd-ldlt:FactorSuccess)
             (sweep-need (rf:refine! refiner (rf-caller:matrix-operator (sweep-stored (sweep-arrow) k))
                                     (rf-caller:factor-solver factor) (sweep-array b4 k) x (sweep-number small k) 1d0)
                         coalton-linalg/refinement:RefineConverged)
             (unless (plusp (rf:step-count refiner))
               (error 'sweep-premise :text "The refinement needed no correction."))
             (list (sweep-out "x" x 0) (sweep-out "the steps" (rf:step-count refiner) 0)
                   (sweep-out "the first residual" (aref (rf:residual-norms refiner) 0) 1)
                   (sweep-out "the last residual"
                              (aref (rf:residual-norms refiner) (1- (rf:residual-count refiner))) 1
                              :kind :within :bound small))))))

     (sweep-entry "the rows that depend on the rows before them"
         (rows:find-dependent-rows! rows:kept-values! rows:null-vector! rows:residual-bound rows:gap-bound
          rows:transpose-product-bound)
         1d-9
         ((-1020 -1020 :no-input) (-1000 -1000 :looser) (-800 1000 :exact) (1020 1020 :looser))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-dependent) k)) (analysis (rows:analyze matrix))
                (kept (zeros (length (csc:entries matrix)))) (v (zeros 4)))
           (sweep-need (rows:find-dependent-rows! analysis matrix (sweep-number (scale-float 1d0 -20) k))
                       coalton-linalg/sparse/dependent-rows:RowsFound)
           (when (and (zerop k) (/= 1 (rows:left-out-count analysis)))
             (error 'sweep-premise :text (format nil "~D rows were left out, not 1." (rows:left-out-count analysis))))
           (unless (rows:kept-values! analysis kept) (error "No kept set stands."))
           (when (plusp (rows:left-out-count analysis))
             (sweep-need (rows:null-vector! analysis 0 v) coalton-linalg/sparse/dependent-rows:NullFound))
           (list (sweep-out "how many rows are left out" (rows:left-out-count analysis) 0)
                 (sweep-out "the row left out" (aref (rows:left-out-rows analysis) 0) 0)
                 (sweep-out "the kept values" kept 1) (sweep-out "the null vector" v 0)
                 (sweep-out "the residual bound" (rows:residual-bound analysis 0) 1 :kind :bound)
                 (sweep-out "the gap bound" (rows:gap-bound analysis) 1 :kind :bound)
                 (sweep-out "a bound on a transposed product"
                            (rows:transpose-product-bound matrix (doubles 1d0 2d0 -1d0 3d0)) 1 :kind :bound)))))

     (sweep-entry "a row that repeats another, under a threshold of zero"
         (rows:find-dependent-rows! rows:kept-values! rows:null-vector!) 1d-9
         ((-1020 1020 :exact))
       (lambda (k)
         (let* ((matrix (sweep-stored (sweep-twice) k)) (analysis (rows:analyze matrix))
                (kept (zeros (length (csc:entries matrix)))) (v (zeros 4)))
           (sweep-need (rows:find-dependent-rows! analysis matrix 0d0) coalton-linalg/sparse/dependent-rows:RowsFound)
           (when (and (zerop k) (/= 1 (rows:left-out-count analysis)))
             (error 'sweep-premise :text (format nil "~D rows were left out, not 1." (rows:left-out-count analysis))))
           (unless (rows:kept-values! analysis kept) (error "No kept set stands."))
           (when (plusp (rows:left-out-count analysis))
             (sweep-need (rows:null-vector! analysis 0 v) coalton-linalg/sparse/dependent-rows:NullFound))
           (list (sweep-out "how many rows are left out" (rows:left-out-count analysis) 0)
                 (sweep-out "the row left out" (aref (rows:left-out-rows analysis) 0) 0)
                 (sweep-out "the kept values" kept 1) (sweep-out "the null vector" v 0))))))))

;;; ---- The roster against the sources, and its run -------------------------------------------

(defun sweep-member-p (type-text)
  "A member of the family: a function whose type mentions a floating-point number, or whose result is a status."
  (and (search "->" type-text)
       (or (search "F64" type-text)
           (let ((end (search "Status)" type-text :from-end t))) (and end (= end (- (length type-text) 7)))))))

(defun sweep-key (symbol)
  (format nil "~(~A~):~(~A~)" (package-name (symbol-package symbol)) (symbol-name symbol)))

(defun sweep-mentions-p (form symbol)
  (cond ((eq form symbol) t)
        ((consp form) (or (sweep-mentions-p (car form) symbol) (sweep-mentions-p (cdr form) symbol)))
        (t nil)))

(defun sweep-check-entry (entry)
  "Run one entry across the sweep. Returns its ranges; prints them for an entry still to be declared."
  (destructuring-bind (name covers tolerance declared run source) entry
    (dolist (symbol covers)
      (assert (sweep-mentions-p source symbol) () "~A answers for ~S and does not call it." name symbol))
    (let* ((reference (funcall run 0)) (decisions '())
           (classes (loop :for k :in *sweep-exponents*
                          :collect (multiple-value-bind (class decided) (sweep-class run reference k tolerance)
                                     (unless (eq class :exact) (push (list k class decided) decisions))
                                     (list k class))))
           (segments (sweep-segments classes)))
      (assert (eq :exact (second (assoc 0 classes))) () "~A: two runs at scale 1 differ." name)
      (cond
        ((eq declared :measure)
         (format t "SWEEP-MEASURED ~S~%    ~S~%" name segments)
         (dolist (decision (reverse decisions)) (format t "      ~S~%" decision)))
        (t
         (assert (not (find :wrong classes :key #'second)) ()
                 "~A: a run succeeds with an output that is not its scaled value: ~S." name
                 (find :wrong decisions :key #'second))
         (assert (equal segments declared) ()
                 "~A: across the sweep the runs are~%    ~S~%and the roster declares~%    ~S" name segments declared)))
      segments)))

(defun run-scale-sweep-tests ()
  (let* ((roster (sweep-roster))
         (inventory (source-inventory #'sweep-member-p))
         (answered (remove-duplicates (loop :for entry :in roster :append (mapcar #'sweep-key (second entry)))
                                      :test #'string=))
         (measured 0) (runs 0))
    (assert (> (length inventory) 60) () "The inventory found only ~D calls." (length inventory))
    (multiple-value-bind (unanswered unknown) (alias-missing inventory answered)
      (assert (null unanswered) () "Calls with floating-point arguments or results that no entry of the sweep answers for: ~S."
              unanswered)
      (assert (null unknown) () "The sweep answers for names that are not such calls: ~S." unknown))
    (assert (equal '("coalton-linalg/new:call!") (alias-missing (cons "coalton-linalg/new:call!" inventory) answered)))
    (dolist (entry roster)
      (when (eq (fourth entry) :measure) (incf measured))
      (sweep-check-entry entry)
      (incf runs (length *sweep-exponents*)))
    ;; Control: an entry that answers for a call its run does not make is refused.
    (assert (handler-case
                (progn (sweep-check-entry
                        (sweep-entry "a control" (numeric:norm2 numeric:norm-inf) 1d-15 ((-1020 1020 :exact))
                          (lambda (k) (list (sweep-out "norm2" (numeric:norm2 (sweep-array (doubles 3d0 4d0) k)) 1)))))
                       nil)
              (error () t))
            () "An entry answered for a call it does not make.")
    ;; Controls for the classes: a value off its scaled value, a value the format cannot hold, a
    ;; refusal, and a bound that came out larger are each put in their class.
    (flet ((classes (run) (mapcar #'second (loop :with reference := (funcall run 0)
                                                 :for k :in '(0 100 1000)
                                                 :collect (list k (sweep-class run reference k 1d-12))))))
      (assert (equal '(:exact :exact :exact) (classes (lambda (k) (list (sweep-out "v" (sweep-number 3d0 k) 1))))))
      (assert (equal '(:exact :wrong :wrong)
                     (classes (lambda (k) (list (sweep-out "v" (if (zerop k) 3d0 (* 1.5d0 (sweep-number 3d0 k))) 1))))))
      (assert (equal '(:exact :close :close)
                     (classes (lambda (k) (list (sweep-out "v" (if (zerop k) 3d0 (* (+ 1d0 1d-14) (sweep-number 3d0 k))) 1))))))
      (assert (equal '(:exact :exact :lost) (classes (lambda (k) (list (sweep-out "v" (if (= k 1000) 0d0 (sweep-number 3d0 (* 2 k))) 2))))))
      (assert (equal '(:exact :exact :wrong)
                     (classes (lambda (k) (list (sweep-out "v" (if (= k 1000) sb-ext:double-float-positive-infinity
                                                                 (sweep-number 3d0 k))
                                                           1))))))
      (assert (equal '(:exact :limit :limit) (classes (lambda (k) (if (zerop k) (list (sweep-out "v" 3d0 1)) (error "no"))))))
      (assert (equal '(:exact :limit :refused) (classes (lambda (k) (if (zerop k) (list (sweep-out "v" 3d0 2)) (error "no"))))))
      (assert (equal '(:exact :looser :looser)
                     (classes (lambda (k) (list (sweep-out "v" (if (zerop k) 3d0 sb-ext:double-float-positive-infinity) 2
                                                           :kind :bound))))))
      (assert (equal '(:exact :close :wrong)
                     (classes (lambda (k) (list (sweep-out "v" (if (zerop k) (doubles 3d0 1d-20)
                                                                 (doubles (sweep-number 3d0 k)
                                                                          (if (= k 100) 0d0 (sweep-number 1d0 k))))
                                                           1))))))
      (assert (equal '(:exact :looser :looser)
                     (classes (lambda (k) (list (sweep-out "v" (if (zerop k) 3d0 (* 4d0 (sweep-number 3d0 k))) 1 :kind :bound))))))
      (assert (equal '(:exact :wrong :wrong)
                     (classes (lambda (k) (list (sweep-out "v" (if (zerop k) 3d0 (* 0.5d0 (sweep-number 3d0 k))) 1 :kind :bound))))))
      (assert (equal '(:exact :exact :no-input) (classes (lambda (k) (list (sweep-out "v" (sweep-number 1d10 (* 2 k)) 2)))))))
    (format t "SCALE-SWEEP-~A entries=~D calls=~D runs=~D~%" (if (zerop measured) "PASS" "MEASURED") (length roster)
            (length answered) runs)))
