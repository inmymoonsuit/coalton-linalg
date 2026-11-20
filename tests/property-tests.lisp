(in-package #:coalton-linalg/tests)

(defun check-small-factor (matrix factor &optional (tolerance 1d-12))
  "Independently multiply the exported factors and compare every matrix entry."
  (let* ((n (array-dimension matrix 0))
         (permutation (coalton-linalg/sparse/spd-ldlt:factor-permutation factor))
         (columns (coalton-linalg/sparse/spd-ldlt:factor-columns factor))
         (rows (coalton-linalg/sparse/spd-ldlt:factor-rows factor))
         (values (coalton-linalg/sparse/spd-ldlt:factor-values factor))
         (diagonal (coalton-linalg/sparse/spd-ldlt:factor-diagonal factor))
         (lower (make-array (list n n) :element-type 'double-float :initial-element 0d0))
         (largest 0d0)
         (error-max 0d0))
    (dotimes (i n)
      (assert (and (> (aref diagonal i) 0d0)
                   (not (sb-ext:float-infinity-p (aref diagonal i)))))
      (setf (aref lower i i) 1d0)
      (loop :for q :from (aref columns i) :below (aref columns (1+ i)) :do
        (setf (aref lower (aref rows q) i) (aref values q))))
    (dotimes (i n)
      (dotimes (j n)
        (let* ((expected (aref matrix (aref permutation i) (aref permutation j)))
               (actual (loop :for k :below n
                             :sum (* (aref lower i k) (aref diagonal k) (aref lower j k))
                               :of-type double-float)))
          (setf largest (max largest (abs expected))
                error-max (max error-max (abs (- actual expected)))))))
    (assert (<= error-max (* tolerance (max least-positive-normalized-double-float largest))) ()
            "Independent small factor product differs by ~E." error-max)))

(defun repeated-refactor-properties ()
  (let* ((symbolic
           (coalton-linalg/sparse/spd-ldlt:analyze
             (coalton-linalg/sparse/spd-ldlt:pattern
               4 (indices 0 2 4 6 7) (indices 0 1 1 2 2 3 3))))
         (factor (coalton-linalg/sparse/spd-ldlt:symbolic-factor symbolic))
         (values (doubles 0d0 0d0 0d0 0d0 0d0 0d0 0d0))
         (answer (doubles 1d0 -2d0 3d0 -4d0))
         (rhs (doubles 0d0 0d0 0d0 0d0))
         (out (doubles 0d0 0d0 0d0 0d0))
         (matrix (make-array '(4 4) :element-type 'double-float :initial-element 0d0)))
    (dotimes (step 96)
      (let* ((d0 (+ 4d0 (mod step 7)))
             (d1 (+ 5d0 (mod step 5)))
             (d2 (+ 6d0 (mod step 3)))
             (d3 (+ 7d0 (mod step 11)))
             (e0 (if (evenp step) 0.25d0 -0.5d0))
             (e1 (if (zerop (mod step 3)) 0d0 0.75d0))
             (e2 (if (evenp step) -0.125d0 0.125d0))
             (shift (if (evenp step) 0d0 0.125d0)))
        (setf (aref values 0) d0 (aref values 1) e0 (aref values 2) d1
              (aref values 3) e1 (aref values 4) d2 (aref values 5) e2
              (aref values 6) d3)
        (setf (aref matrix 0 0) (- d0 shift)
              (aref matrix 1 1) (- d1 shift)
              (aref matrix 2 2) (- d2 shift)
              (aref matrix 3 3) (- d3 shift)
              (aref matrix 0 1) e0 (aref matrix 1 0) e0
              (aref matrix 1 2) e1 (aref matrix 2 1) e1
              (aref matrix 2 3) e2 (aref matrix 3 2) e2)
        (dotimes (i 4)
          (setf (aref rhs i)
                (loop :for j :below 4 :sum (* (aref matrix i j) (aref answer j)) :of-type double-float)))
        (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess) (coalton-linalg/sparse/spd-ldlt:factor! symbolic values shift)))
        (check-small-factor matrix factor)
        (fill out 42d0)
        (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess) (coalton-linalg/sparse/spd-ldlt:solve! factor rhs out)))
        (assert (close-vector-p out answer 2d-13))
        (let ((aliased (copy-seq rhs)))
          (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess) (coalton-linalg/sparse/spd-ldlt:solve! factor aliased aliased)))
          (assert (close-vector-p aliased answer 2d-13)))
        (when (zerop (mod step 8))
          (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive) (coalton-linalg/sparse/spd-ldlt:factor! symbolic values 1000d0)))
          (fill out 42d0)
          (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveInvalidFactor) (coalton-linalg/sparse/spd-ldlt:solve! factor rhs out)))
          (assert (every (lambda (x) (= x 42d0)) out)))))
    (format t "REFACTOR-PROPERTIES-PASS cases=96~%")))

(defun explicit-floating-boundary-tests ()
  "Test binary input/shift boundaries without calling a floating pivot an exact rank certificate."
  (let* ((symbolic
           (coalton-linalg/sparse/spd-ldlt:analyze
             (coalton-linalg/sparse/spd-ldlt:pattern 2 (indices 0 2 3) (indices 0 1 1))))
         (rounded (+ 1d0 (scale-float 1d0 -54)))
         (values (doubles 1d0 1d0 rounded)))
    (assert (= rounded 1d0))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive) (coalton-linalg/sparse/spd-ldlt:factor! symbolic values 0d0))))
  (let* ((symbolic
           (coalton-linalg/sparse/spd-ldlt:analyze
             (coalton-linalg/sparse/spd-ldlt:pattern 1 (indices 0 1) (indices 0))))
         (values (doubles 1d0)))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess) (coalton-linalg/sparse/spd-ldlt:factor!
                      symbolic values (- 1d0 (scale-float 1d0 -53)))))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive) (coalton-linalg/sparse/spd-ldlt:factor! symbolic values 1d0)))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive) (coalton-linalg/sparse/spd-ldlt:factor!
                     symbolic values (+ 1d0 (scale-float 1d0 -52)))))
    (setf (aref values 0) least-positive-double-float)
    (let* ((result (coalton-linalg/sparse/spd-ldlt:factor! symbolic values 0d0))
           (zero (doubles 0d0)))
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess) result))
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess) (coalton-linalg/sparse/spd-ldlt:solve!
                      (coalton-linalg/sparse/spd-ldlt:symbolic-factor symbolic) zero zero)))))
  (format t "EXPLICIT-FLOATING-BOUNDARIES-PASS~%"))

(defun rounded-gram-observations ()
  "Bound the behavior on nearly singular rounded inputs; either finite positive pivots or refusal is valid."
  (let ((accepted 0) (refused 0))
    (dotimes (seed 24)
      (let* ((u (doubles (+ 0.1d0 (* seed 0.01d0)) 0.3d0 0d0))
             (v (doubles 0.2d0 (+ 0.4d0 (* seed 0.01d0)) 0d0))
             (matrix (make-array '(3 3) :element-type 'double-float))
             (values (doubles 0d0 0d0 0d0 0d0 0d0 0d0))
             (symbolic
               (coalton-linalg/sparse/spd-ldlt:analyze
                 (coalton-linalg/sparse/spd-ldlt:pattern
                   3 (indices 0 3 5 6) (indices 0 1 2 1 2 2)))))
        (setf (aref u 2) (+ (aref u 0) (aref u 1))
              (aref v 2) (+ (aref v 0) (aref v 1)))
        (dotimes (i 3)
          (dotimes (j 3)
            (setf (aref matrix i j)
                  (+ (* (aref u i) (aref u j)) (* (aref v i) (aref v j))))))
        (loop :for (row col) :in '((0 0) (1 0) (2 0) (1 1) (2 1) (2 2))
              :for index :from 0 :do
                (setf (aref values index) (aref matrix row col)))
        (let ((status (coalton-linalg/sparse/spd-ldlt:factor! symbolic values 0d0)))
          (cond
            ((eq status (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess))
             (incf accepted)
             (check-small-factor matrix (coalton-linalg/sparse/spd-ldlt:symbolic-factor symbolic)))
            ((eq status (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive))
             (incf refused)
             (assert (<= 0 (coalton-linalg/sparse/spd-ldlt:refused-column symbolic) 2))
             (assert (not (coalton-linalg/sparse/spd-ldlt:factor-ready?
                            (coalton-linalg/sparse/spd-ldlt:symbolic-factor symbolic)))))
            (t (error "Unexpected boundary-observation status ~S." status))))))
    (format t "ROUNDED-GRAM-OBSERVATIONS-PASS accepted=~D refused=~D~%" accepted refused)))


(defun exact-three-by-three-determinant (matrix)
  "Compute a 3x3 determinant from the exact rational values of its stored F64 entries."
  (assert (equal (array-dimensions matrix) '(3 3)))
  (flet ((q (row column) (rational (aref matrix row column))))
    (+ (* (q 0 0) (- (* (q 1 1) (q 2 2)) (* (q 1 2) (q 2 1))))
       (- (* (q 0 1) (- (* (q 1 0) (q 2 2)) (* (q 1 2) (q 2 0)))))
       (* (q 0 2) (- (* (q 1 0) (q 2 1)) (* (q 1 1) (q 2 0)))))))

(defun exact-singular-rounding-observation ()
  "Prove singularity of serialized F64 input, then observe its positive computed pivots."
  (let ((data
          (with-open-file
              (stream (asdf:system-relative-pathname
                       "coalton-linalg" "tests/exact-singular-rounding.sexp")
                      :direction :input)
            (assert (< (file-length stream) 4096))
            (let* ((*read-eval* nil)
                   (eof (gensym "EOF"))
                   (value (read stream nil eof)))
              (assert (not (eq value eof)))
              (assert (eq (read stream nil eof) eof))
              value))))
    (destructuring-bind
        (&key size columns rows values null-vector control-last-diagonal) data
      (assert (= size 3))
      (assert (equal columns '(0 3 5 6)))
      (assert (equal rows '(0 1 2 1 2 2)))
      (assert (and (= (length values) 6) (every (lambda (x) (typep x 'double-float)) values)))
      (assert (and (= (length null-vector) 3) (every #'integerp null-vector)
                   (some (lambda (x) (not (zerop x))) null-vector)))
      (assert (typep control-last-diagonal 'double-float))
      (let* ((pointers (apply #'indices columns))
             (row-indices (apply #'indices rows))
             (entries (apply #'doubles values))
             (matrix (make-array '(3 3) :element-type 'double-float :initial-element 0d0))
             (symbolic
               (coalton-linalg/sparse/spd-ldlt:analyze
                 (coalton-linalg/sparse/spd-ldlt:pattern size pointers row-indices)))
             (factor (coalton-linalg/sparse/spd-ldlt:symbolic-factor symbolic)))
        (dotimes (column size)
          (loop :for q :from (aref pointers column) :below (aref pointers (1+ column))
                :for row := (aref row-indices q)
                :do (setf (aref matrix row column) (aref entries q)
                         (aref matrix column row) (aref entries q))))
        (flet ((kernel-residual ()
                 (loop :for row :below size
                       :collect (loop :for column :below size
                                     :for x :in null-vector
                                     :sum (* (rational (aref matrix row column)) x)))))
          (assert (= (rational (aref matrix 0 0)) 3))
          (assert (= (- (* (rational (aref matrix 0 0)) (rational (aref matrix 1 1)))
                        (* (rational (aref matrix 0 1)) (rational (aref matrix 1 0))))
                     2))
          (assert (zerop (exact-three-by-three-determinant matrix)))
          (assert (every #'zerop (kernel-residual)))
          ;; A nonzero 2x2 minor and zero 3x3 determinant prove rank exactly two.
          (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess)
                      (coalton-linalg/sparse/spd-ldlt:factor! symbolic entries 0d0)))
          (assert (= -1 (coalton-linalg/sparse/spd-ldlt:refused-column symbolic)))
          (assert (equalp (coalton-linalg/sparse/spd-ldlt:factor-permutation factor) #(0 1 2)))
          (check-small-factor matrix factor)
          (let* ((diagonal (coalton-linalg/sparse/spd-ldlt:factor-diagonal factor))
                 (minimum-pivot (reduce #'min diagonal))
                 (matrix-scale (reduce #'max entries :key #'abs)))
            (assert (< 0d0 minimum-pivot (* 8d0 double-float-epsilon matrix-scale)))
            (format t "EXACT-SINGULAR-ROUNDING rank=2 determinant=0 min-pivot=~E~%" minimum-pivot))
          (setf (aref entries 5) control-last-diagonal
                (aref matrix 2 2) control-last-diagonal)
          (assert (= (exact-three-by-three-determinant matrix) 2))
          (assert (not (every #'zerop (kernel-residual))))
          (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess)
                      (coalton-linalg/sparse/spd-ldlt:factor! symbolic entries 0d0)))
          (assert (= -1 (coalton-linalg/sparse/spd-ldlt:refused-column symbolic)))
          (check-small-factor matrix factor)
          (let ((out (doubles 0d0 0d0 0d0)))
            (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess)
                        (coalton-linalg/sparse/spd-ldlt:solve! factor (doubles 13d0 5d0 21d0) out)))
            (assert (close-vector-p out (doubles 1d0 -2d0 3d0) 2d-13)))
          (format t "EXACT-SINGULAR-CONTROL-PASS rank=3 determinant=2~%"))))))

(defun run-property-tests ()
  (repeated-refactor-properties)
  (explicit-floating-boundary-tests)
  (rounded-gram-observations)
  (exact-singular-rounding-observation)
  (format t "NUMERICAL-PROPERTIES-PASS~%"))
