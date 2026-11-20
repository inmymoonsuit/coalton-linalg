;;;; Public-interface analysis measurements; no numerical factorization is run.
(defpackage #:coalton-linalg/sparse-analysis-tests
  (:use #:cl)
  (:export #:run-case #:run-directory))
(in-package #:coalton-linalg/sparse-analysis-tests)

(defun read-pattern (path)
  "Read bounded integer-only CSC data; parsing is outside the measured analysis."
  (with-open-file (stream path :direction :input)
    (let ((*read-eval* nil)
          (eof (gensym "EOF")))
      (flet ((next-index ()
               (let ((value (read stream nil eof)))
                 (unless (typep value '(integer 0 #.(1- (expt 2 62))))
                   (error "Expected a nonnegative UFix in ~A, got ~S." path value))
                 value)))
        (let* ((n (next-index))
               (nnz (next-index)))
          (unless (and (<= n 1000000) (<= n nnz 10000000))
            (error "Analysis test input exceeds the reader's dimension/entry limits: ~A." path))
          (let ((columns (make-array (1+ n) :element-type '(unsigned-byte 62)))
                (rows (make-array nnz :element-type '(unsigned-byte 62))))
            (dotimes (i (length columns)) (setf (aref columns i) (next-index)))
            (dotimes (i nnz) (setf (aref rows i) (next-index)))
            (unless (eq (read stream nil eof) eof)
              (error "Trailing analysis input in ~A." path))
            (values (coalton-linalg/sparse/spd-ldlt:Pattern n columns rows) n nnz)))))))

(defun case-name (path)
  "Remove the dedicated .pattern.txt suffix from a generated case filename."
  (let* ((filename (file-namestring path))
         (suffix ".pattern.txt")
         (end (- (length filename) (length suffix))))
    (unless (and (> end 0) (string= suffix filename :start2 end))
      (error "Expected NAME.pattern.txt, got ~A." filename))
    (subseq filename 0 end)))

(defun output-path (directory name suffix)
  (merge-pathnames (concatenate 'string name suffix)
                   (uiop:ensure-directory-pathname directory)))

(defun write-indices (path vector)
  (with-open-file (stream path :direction :output :if-exists :supersede
                               :if-does-not-exist :create)
    (loop for value across vector do (format stream "~D~%" value))))

(defun run-case (input-path output-directory)
  "Analyze one pattern and export its public structure plus scoped timing/allocation.
Elapsed time and SBCL's program-wide bytes-consed delta cover only analyze.
They do not measure peak memory. External process measurements need their own
baseline and scope. One call is measured; compilation, parsing and output are excluded."
  (let ((name (case-name input-path)))
    (ensure-directories-exist (output-path output-directory name ".analysis.meta.txt"))
    (dolist (suffix '(".perm.out.txt" ".columns.out.txt" ".rows.out.txt" ".analysis.meta.txt"))
      (let ((old (output-path output-directory name suffix)))
        (when (probe-file old) (delete-file old))))
    (multiple-value-bind (pattern n input-nnz) (read-pattern input-path)
      (sb-ext:gc :full t)
      (format t "SPARSE-ANALYSIS-BEGIN ~A n=~D input-nnz=~D~%" name n input-nnz)
      (finish-output)
      (let* ((before-bytes (sb-ext:get-bytes-consed))
             (before-time (get-internal-real-time))
             (symbolic (coalton-linalg/sparse/spd-ldlt:analyze pattern))
             (elapsed (- (get-internal-real-time) before-time))
             (bytes (- (sb-ext:get-bytes-consed) before-bytes))
             (factor (coalton-linalg/sparse/spd-ldlt:symbolic-factor symbolic))
             (permutation (coalton-linalg/sparse/spd-ldlt:factor-permutation factor))
             (columns (coalton-linalg/sparse/spd-ldlt:factor-columns factor))
             (rows (coalton-linalg/sparse/spd-ldlt:factor-rows factor)))
        (assert (not (coalton-linalg/sparse/spd-ldlt:factor-ready? factor)))
        (assert (= -1 (coalton-linalg/sparse/spd-ldlt:refused-column symbolic)))
        (format t "SPARSE-ANALYSIS-END ~A ticks=~D bytes-consed=~D strict-l-nnz=~D~%"
                name elapsed bytes (length rows))
        (finish-output)
        (write-indices (output-path output-directory name ".perm.out.txt") permutation)
        (write-indices (output-path output-directory name ".columns.out.txt") columns)
        (write-indices (output-path output-directory name ".rows.out.txt") rows)
        (with-open-file (stream (output-path output-directory name ".analysis.meta.txt")
                               :direction :output :if-exists :supersede
                               :if-does-not-exist :create)
          (format stream "format: 1~%status: success~%size: ~D~%input_nnz: ~D~%strict_l_nnz: ~D~%"
                  n input-nnz (length rows))
          (format stream "analysis_ticks: ~D~%internal_time_units_per_second: ~D~%"
                  elapsed internal-time-units-per-second)
          (format stream "analysis_seconds: ~,9F~%analysis_bytes_consed: ~D~%"
                  (coerce (/ elapsed internal-time-units-per-second) 'double-float) bytes)
          (format stream "measurement_scope: analyze only; bytes-consed is program-wide allocation, not peak memory~%"))
        (format t "SPARSE-ANALYSIS-PASS ~A~%" name)
        (list :name name :size n :input-nnz input-nnz :strict-l-nnz (length rows)
              :analysis-ticks elapsed :analysis-bytes-consed bytes)))))

(defun run-directory (input-directory output-directory)
  "Run generated pattern files in deterministic filename order, serially."
  (let ((paths (sort (directory (merge-pathnames "*.pattern.txt"
                                                (uiop:ensure-directory-pathname input-directory)))
                     #'string< :key #'file-namestring)))
    (unless paths (error "No generated analysis patterns found in ~A." input-directory))
    (dolist (path paths) (run-case path output-directory))
    (format t "SPARSE-ANALYSIS-DIRECTORY-PASS cases=~D~%" (length paths))))