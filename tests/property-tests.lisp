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
      (loop for q from (aref columns i) below (aref columns (1+ i)) do
        (setf (aref lower (aref rows q) i) (aref values q))))
    (dotimes (i n)
      (dotimes (j n)
        (let* ((expected (aref matrix (aref permutation i) (aref permutation j)))
               (actual (loop for k below n
                             sum (* (aref lower i k) (aref diagonal k) (aref lower j k))
                               of-type double-float)))
          (setf largest (max largest (abs expected))
                error-max (max error-max (abs (- actual expected)))))))
    (assert (<= error-max (* tolerance (max least-positive-normalized-double-float largest))) ()
            "Independent small factor product differs by ~E." error-max)))

(defun repeated-refactor-properties ()
  (let* ((symbolic
           (coalton-linalg/sparse/spd-ldlt:analyse
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
                (loop for j below 4 sum (* (aref matrix i j) (aref answer j)) of-type double-float)))
        (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess) (coalton-linalg/sparse/spd-ldlt:factor-status
                        (coalton-linalg/sparse/spd-ldlt:factor! symbolic values shift))))
        (check-small-factor matrix factor)
        (fill out 42d0)
        (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess) (coalton-linalg/sparse/spd-ldlt:solve! factor rhs out)))
        (assert (close-vector-p out answer 2d-13))
        (let ((aliased (copy-seq rhs)))
          (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess) (coalton-linalg/sparse/spd-ldlt:solve! factor aliased aliased)))
          (assert (close-vector-p aliased answer 2d-13)))
        (when (zerop (mod step 8))
          (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive) (coalton-linalg/sparse/spd-ldlt:factor-status
                         (coalton-linalg/sparse/spd-ldlt:factor! symbolic values 1000d0))))
          (fill out 42d0)
          (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveInvalidFactor) (coalton-linalg/sparse/spd-ldlt:solve! factor rhs out)))
          (assert (every (lambda (x) (= x 42d0)) out)))))
    (format t "REFACTOR-PROPERTIES-PASS cases=96~%")))

(defun explicit-floating-boundary-tests ()
  "Test binary input/shift boundaries without calling a floating pivot an exact rank certificate."
  (let* ((symbolic
           (coalton-linalg/sparse/spd-ldlt:analyse
             (coalton-linalg/sparse/spd-ldlt:pattern 2 (indices 0 2 3) (indices 0 1 1))))
         (rounded (+ 1d0 (scale-float 1d0 -54)))
         (values (doubles 1d0 1d0 rounded)))
    (assert (= rounded 1d0))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive) (coalton-linalg/sparse/spd-ldlt:factor-status
                   (coalton-linalg/sparse/spd-ldlt:factor! symbolic values 0d0)))))
  (let* ((symbolic
           (coalton-linalg/sparse/spd-ldlt:analyse
             (coalton-linalg/sparse/spd-ldlt:pattern 1 (indices 0 1) (indices 0))))
         (values (doubles 1d0)))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess) (coalton-linalg/sparse/spd-ldlt:factor-status
                    (coalton-linalg/sparse/spd-ldlt:factor!
                      symbolic values (- 1d0 (scale-float 1d0 -53))))))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive) (coalton-linalg/sparse/spd-ldlt:factor-status
                   (coalton-linalg/sparse/spd-ldlt:factor! symbolic values 1d0))))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive) (coalton-linalg/sparse/spd-ldlt:factor-status
                   (coalton-linalg/sparse/spd-ldlt:factor!
                     symbolic values (+ 1d0 (scale-float 1d0 -52))))))
    (setf (aref values 0) least-positive-double-float)
    (let* ((result (coalton-linalg/sparse/spd-ldlt:factor! symbolic values 0d0))
           (zero (doubles 0d0)))
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess) (coalton-linalg/sparse/spd-ldlt:factor-status result)))
      (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess) (coalton-linalg/sparse/spd-ldlt:solve!
                      (coalton-linalg/sparse/spd-ldlt:result-factor result) zero zero)))))
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
               (coalton-linalg/sparse/spd-ldlt:analyse
                 (coalton-linalg/sparse/spd-ldlt:pattern
                   3 (indices 0 3 5 6) (indices 0 1 2 1 2 2)))))
        (setf (aref u 2) (+ (aref u 0) (aref u 1))
              (aref v 2) (+ (aref v 0) (aref v 1)))
        (dotimes (i 3)
          (dotimes (j 3)
            (setf (aref matrix i j)
                  (+ (* (aref u i) (aref u j)) (* (aref v i) (aref v j))))))
        (loop for (row col) in '((0 0) (1 0) (2 0) (1 1) (2 1) (2 2))
              for index from 0 do
                (setf (aref values index) (aref matrix row col)))
        (let* ((result (coalton-linalg/sparse/spd-ldlt:factor! symbolic values 0d0))
               (status (coalton-linalg/sparse/spd-ldlt:factor-status result)))
          (cond
            ((eq status (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess))
             (incf accepted)
             (check-small-factor matrix (coalton-linalg/sparse/spd-ldlt:result-factor result)))
            ((eq status (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorNonPositive))
             (incf refused)
             (assert (<= 0 (coalton-linalg/sparse/spd-ldlt:result-column result) 2))
             (assert (not (coalton-linalg/sparse/spd-ldlt:factor-ready?
                            (coalton-linalg/sparse/spd-ldlt:symbolic-factor symbolic)))))
            (t (error "Unexpected boundary-observation status ~S." status))))))
    (format t "ROUNDED-GRAM-OBSERVATIONS-PASS accepted=~D refused=~D~%" accepted refused)))

(defun run-property-tests ()
  (repeated-refactor-properties)
  (explicit-floating-boundary-tests)
  (rounded-gram-observations)
  (format t "NUMERICAL-PROPERTIES-PASS~%"))
