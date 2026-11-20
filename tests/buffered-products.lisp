(in-package #:coalton-linalg/tests)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "CSC" "COALTON-LINALG/SPARSE/CSC")
  (sb-ext:add-package-local-nickname "BUFFERED" "COALTON-LINALG/TESTS/BUFFERED-PRODUCTS"))

(defun buffered-output (size)
  "Use nonfinite old contents so success must overwrite rather than accumulate."
  (let ((output (make-array size :element-type 'double-float
                           :initial-element (sb-kernel:make-double-float #x7ff80000 37))))
    (when (> size 1) (setf (aref output 1) -0d0))
    output))

(defun buffered-observe (label matrix input output transposed expected
                         &key (preserve-output t))
  "Classify a real typed call and preserve all inputs; preflight refusals also preserve output bits."
  (let* ((input-before (copy-seq input))
         (output-before (copy-seq output))
         (pointers-before (copy-seq (csc:column-pointers matrix)))
         (rows-before (copy-seq (csc:row-indices matrix)))
         (entries-before (copy-seq (csc:entries matrix)))
         (actual (funcall (if transposed #'buffered:transposed-error #'buffered:forward-error)
                          matrix input output)))
    (assert (= expected actual) () "~A expected typed code ~D, got ~D." label expected actual)
    (assert (same-bits-p input input-before) () "~A changed the input vector." label)
    (assert (equalp pointers-before (csc:column-pointers matrix)) () "~A changed pointers." label)
    (assert (equalp rows-before (csc:row-indices matrix)) () "~A changed row indices." label)
    (assert (same-bits-p entries-before (csc:entries matrix)) () "~A changed matrix entries." label)
    (when preserve-output
      (assert (same-bits-p output output-before) () "~A changed output on preflight refusal." label))
    actual))

(defun buffered-success (label matrix input output transposed expected)
  (buffered-observe label matrix input output transposed 0 :preserve-output nil)
  (assert (same-bits-p output expected) () "~A produced ~S, expected ~S." label output expected))

(defun buffered-rectangle ()
  "A = [[2,0,-1],[0,3,4]], including an explicitly stored zero."
  (csc:from-csc 2 3 (indices 0 2 3 5) (indices 0 1 1 0 1)
                (doubles 2d0 0d0 3d0 -1d0 4d0)))

(defun buffered-success-tests ()
  (let ((matrix (buffered-rectangle)))
    (buffered-success "rectangular Ax overwrite" matrix (doubles 1d0 2d0 3d0)
                      (buffered-output 2) nil (doubles -1d0 18d0))
    (buffered-success "rectangular transpose Ax overwrite" matrix (doubles -2d0 5d0)
                      (buffered-output 3) t (doubles -4d0 15d0 22d0)))
  (dolist (shape '((0 0) (0 3) (4 0) (2 3)))
    (destructuring-bind (m n) shape
      (let ((matrix (csc:from-csc m n
                      (make-array (1+ n) :element-type '(unsigned-byte 62) :initial-element 0)
                      (indices) (doubles))))
        (dolist (transposed '(nil t))
          (let* ((input-size (if transposed m n))
                 (output-size (if transposed n m))
                 (input (make-array input-size :element-type 'double-float :initial-element 1d0))
                 (output (buffered-output output-size))
                 (expected (make-array output-size :element-type 'double-float :initial-element 0d0)))
            (buffered-success (format nil "empty/zero shape ~Dx~D transpose=~S" m n transposed)
                              matrix input output transposed expected))))))
  (let* ((matrix (csc:from-csc 2 2 (indices 0 1 2) (indices 0 1) (doubles 2d0 3d0)))
         (input (csc:entries matrix)))
    (dolist (transposed '(nil t))
      (buffered-success "matrix entries are a permitted read alias"
                        matrix input (buffered-output 2) transposed (doubles 4d0 9d0))
      ;; Equality of contents is not identity or overlapping storage.
      (let ((other-input (doubles 2d0 3d0))
            (output (doubles 2d0 3d0)))
        (assert (not (eq other-input output)))
        (assert (not (eq (csc:entries matrix) output)))
        (buffered-success "equal-valued distinct arrays" matrix other-input output
                          transposed (doubles 4d0 9d0)))))
  (format t "BUFFERED-SUCCESS-PASS rectangular=2 empty-or-zero=8 read-alias=2 distinct-equal=2~%"))

(defun buffered-preflight-tests ()
  (let ((matrix (buffered-rectangle)))
    (buffered-observe "Ax wrong input length" matrix (doubles 1d0 2d0) (buffered-output 2) nil 2)
    (buffered-observe "Ax wrong output length" matrix (doubles 1d0 2d0 3d0) (buffered-output 3) nil 2)
    (buffered-observe "transpose Ax wrong input length" matrix (doubles 1d0) (buffered-output 3) t 2)
    (buffered-observe "transpose Ax wrong output length" matrix (doubles 1d0 2d0) (buffered-output 2) t 2))
  (let ((matrix (csc:from-csc 2 2 (indices 0 1 2) (indices 0 1) (doubles 2d0 3d0))))
    (dolist (transposed '(nil t))
      (let ((input (doubles 1d0 -0d0)))
        (buffered-observe "input/output identity" matrix input input transposed 9))
      (buffered-observe "output/matrix-entry identity" matrix (doubles 1d0 2d0)
                        (csc:entries matrix) transposed 9)
      (buffered-observe "output aliases both read inputs" matrix (csc:entries matrix)
                        (csc:entries matrix) transposed 9)
      ;; Alias rejection precedes the full finite-input scan.
      (let ((input (buffered-output 2)))
        (buffered-observe "alias before nonfinite-input validation" matrix input input transposed 9))
      ;; Dimensions precede aliasing, while neither branch may write anything.
      (let ((input (buffered-output 1)))
        (buffered-observe "dimensions before alias validation" matrix input input transposed 2))))
  (let ((matrix (csc:from-csc 0 0 (indices 0) (indices) (doubles))))
    (dolist (transposed '(nil t))
      (let ((input (doubles)))
        (buffered-observe "empty input/output identity" matrix input input transposed 9))
      (buffered-observe "empty output/matrix-entry identity" matrix (doubles)
                        (csc:entries matrix) transposed 9)))
  (let ((matrix (csc:from-csc 2 3 (indices 0 1 1 1) (indices 0) (doubles 1d0))))
    (dolist (transposed '(nil t))
      (dolist (bad (list sb-ext:double-float-positive-infinity
                        sb-ext:double-float-negative-infinity
                        (sb-kernel:make-double-float #x7ff80000 41)))
        (dolist (used '(nil t))
          (let* ((size (if transposed 2 3))
                 (input (make-array size :element-type 'double-float :initial-element 1d0)))
            (setf (aref input (if used 0 (1- size))) bad)
            (buffered-observe "used or unused nonfinite input" matrix input
                              (buffered-output (if transposed 3 2)) transposed 5))))))
  ;; Even an operation producing no output validates every input component.
  (dolist (transposed '(nil t))
    (let ((matrix (if transposed
                      (csc:from-csc 2 0 (indices 0) (indices) (doubles))
                      (csc:from-csc 0 2 (indices 0 0 0) (indices) (doubles)))))
      (buffered-observe "nonfinite input with zero-length output"
                        matrix (doubles 1d0 sb-ext:double-float-positive-infinity)
                        (doubles) transposed 5)))
  (format t "BUFFERED-PREFLIGHT-PASS dimensions=6 aliases=12 nonfinite=14~%"))

(defun buffered-overflow-tests ()
  (dolist (transposed '(nil t))
    (let ((matrix (if transposed
                      (csc:from-csc 2 1 (indices 0 2) (indices 0 1)
                                    (doubles most-positive-double-float most-positive-double-float))
                      (csc:from-csc 1 2 (indices 0 1 2) (indices 0 0)
                                    (doubles most-positive-double-float most-positive-double-float)))))
      (buffered-observe "sum overflow" matrix (doubles 1d0 1d0) (buffered-output 1)
                        transposed 6 :preserve-output nil)
      (buffered-observe "product overflow" matrix (doubles 2d0 0d0) (buffered-output 1)
                        transposed 6 :preserve-output nil))
    (let* ((values (doubles most-positive-double-float most-positive-double-float
                            (- most-positive-double-float)))
           (matrix (if transposed
                       (csc:from-csc 3 1 (indices 0 3) (indices 0 1 2) values)
                       (csc:from-csc 1 3 (indices 0 1 2 3) (indices 0 0 0) values))))
      (buffered-observe "later cancellation does not recover overflow"
                        matrix (doubles 1d0 1d0 1d0) (buffered-output 1)
                        transposed 6 :preserve-output nil)))
  (format t "BUFFERED-ARITHMETIC-REFUSAL-PASS cases=6 output-state=unspecified~%"))

(defun buffered-allocation-case (label matrix input output transposed expected)
  "Measure 10000 successful calls in a real Coalton loop after 10000 warmup calls."
  (let ((caller (if transposed #'buffered:repeat-transposed! #'buffered:repeat-forward!))
        (input-before (copy-seq input))
        (entries-before (copy-seq (csc:entries matrix))))
    (sb-ext:gc :full t)
    (let ((before-warmup (sb-ext:get-bytes-consed)))
      (funcall caller matrix input output 10000)
      (let ((after-warmup (sb-ext:get-bytes-consed)))
        (funcall caller matrix input output 10000)
        (let ((after-measurement (sb-ext:get-bytes-consed)))
          (format t "BUFFERED-ALLOCATION ~A warmup-calls=10000 warmup-bytes=~D calls=10000 bytes=~D~%"
                  label (- after-warmup before-warmup) (- after-measurement after-warmup))
          (assert (zerop (- after-measurement after-warmup)) ()
                  "~A allocated ~D bytes during 10000 successful Coalton calls."
                  label (- after-measurement after-warmup)))))
    (assert (same-bits-p input input-before))
    (assert (same-bits-p (csc:entries matrix) entries-before))
    (assert (same-bits-p output expected) () "~A lost the computed output." label)))

(defun buffered-allocation-tests ()
  (let ((matrix (buffered-rectangle)))
    (buffered-allocation-case "Ax-rectangular" matrix (doubles 1d0 2d0 3d0)
                              (buffered-output 2) nil (doubles -1d0 18d0))
    (buffered-allocation-case "transpose-Ax-rectangular" matrix (doubles -2d0 5d0)
                              (buffered-output 3) t (doubles -4d0 15d0 22d0)))
  (let ((matrix (csc:from-csc 0 3 (indices 0 0 0 0) (indices) (doubles))))
    (buffered-allocation-case "Ax-empty-output" matrix (doubles 1d0 2d0 3d0)
                              (doubles) nil (doubles))
    (buffered-allocation-case "transpose-Ax-empty-input" matrix (doubles)
                              (buffered-output 3) t (doubles 0d0 0d0 0d0)))
  (let* ((matrix (csc:from-csc 2 2 (indices 0 1 2) (indices 0 1) (doubles 2d0 3d0)))
         (input (csc:entries matrix)))
    (buffered-allocation-case "Ax-matrix-entry-input" matrix input
                              (buffered-output 2) nil (doubles 4d0 9d0))
    (buffered-allocation-case "transpose-Ax-matrix-entry-input" matrix input
                              (buffered-output 2) t (doubles 4d0 9d0))))

(defun run-buffered-product-tests ()
  (buffered-success-tests)
  (buffered-preflight-tests)
  (buffered-overflow-tests)
  (buffered-allocation-tests)
  (format t "BUFFERED-PRODUCTS-PASS~%")
  t)