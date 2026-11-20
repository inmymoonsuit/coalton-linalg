(in-package #:coalton-linalg/tools)

;;; Optional artifact driver for fixed-pattern Gram products. For every <name>.J.txt in a folder
;;; it analyzes the orientation named by <name>.meta.txt, refills twice, and writes files an
;;; independent reviewer can recompute from. See docs/gram.md for the calling contract and
;;; its Artifact driver section for the file protocol.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "GRAM" "COALTON-LINALG/SPARSE/GRAM")
  (sb-ext:add-package-local-nickname "CSC" "COALTON-LINALG/SPARSE/CSC"))

(defun read-rectangular (directory name)
  "Read <name>.J.txt: `rows columns nnz`, then zero-based `row column value` triplets in any order."
  (with-open-file (in (fixture-path directory name ".J.txt"))
    (let* ((m (read-numeric in)) (n (read-numeric in)) (nz (read-numeric in)))
      (unless (and (integerp m) (integerp n) (integerp nz) (>= m 0) (>= n 0) (>= nz 0))
        (error "Invalid rectangular header: ~A" name))
      (let ((rows (uarray nz)) (columns (uarray nz)) (vals (darray nz)))
        (dotimes (q nz)
          (let ((i (read-numeric in)) (j (read-numeric in)) (v (read-numeric in)))
            (unless (and (integerp i) (integerp j) (<= 0 i) (<= 0 j) (realp v))
              (error "Malformed rectangular entry ~A/~D." name q))
            (setf (aref rows q) i (aref columns q) j (aref vals q) (coerce v 'double-float))))
        (when (read-numeric in) (error "Trailing rectangular entries: ~A" name))
        (values m n rows columns vals)))))

(defun gram-orientation (directory name)
  "The product <name>.meta.txt names on its `gram:` line: :ata for `Jt J`, :aat for `J Jt`."
  (let ((path (fixture-path directory name ".meta.txt")) (prefix "gram: "))
    (unless (probe-file path) (error "No meta file names a Gram orientation: ~A" name))
    (with-open-file (in path)
      (loop :for line := (read-line in nil nil) :while line
            :when (and (<= (length prefix) (length line)) (string= prefix line :end2 (length prefix)))
              :do (let ((rest (subseq line (length prefix))))
                   (cond ((and (<= 4 (length rest)) (string= "Jt J" rest :end2 4)) (return ':ata))
                         ((and (<= 4 (length rest)) (string= "J Jt" rest :end2 4)) (return ':aat))
                         (t (error "Unknown Gram orientation in ~A: ~A" name rest))))
            :finally (error "No `gram:` line in the meta file: ~A" name)))))

(defparameter *constructor-spellings*
  '("InvalidStructure" "InvalidDimensions" "InvalidPermutation" "InvalidAliasing" "DuplicateEntry"
    "NonFiniteInput" "NonFiniteResult" "AnalysisLimitExceeded" "AnalysisStorageExceeded")
  "LinalgError constructors as the source spells them; the Lisp reader upcases the class names.")

(defun refusal-name (condition)
  "The LinalgError constructor a typed library refusal carries, spelled as in the source."
  (let* ((full (symbol-name (class-name (class-of condition))))
         (slash (position #\/ full :from-end t))
         (bare (if slash (subseq full (1+ slash)) full)))
    (or (find bare *constructor-spellings* :test #'string-equal) bare)))

(defun write-gram-file (path plan vals)
  (let ((cp (gram:column-pointers plan)) (ri (gram:row-indices plan)))
    (with-open-file (out path :direction :output :if-exists :supersede :if-does-not-exist :create)
      (format out "~D ~D~%" (gram:dimension plan) (length ri))
      (dotimes (j (gram:dimension plan))
        (loop :for q :from (aref cp j) :below (aref cp (1+ j)) :do
          (format out "~D ~D ~A~%" (aref ri q) j (f64-string (aref vals q))))))))

(defun run-gram-fixtures (input output &key (repeats 20) (allocation-repeats 1000))
  "Write each system's lower Gram product, a rescaled second refill, and one measurement line."
  (ensure-directories-exist (uiop:ensure-directory-pathname output))
  (let ((paths (sort (directory (merge-pathnames "*.J.txt" (uiop:ensure-directory-pathname input)))
                     #'string< :key #'namestring))
        (failures 0) (count 0)
        (linalg-error (find-symbol "LINALGERROR" "COALTON-LINALG/ERRORS")))
    (with-open-file (table (merge-pathnames "gram-table.txt" (uiop:ensure-directory-pathname output))
                           :direction :output :if-exists :supersede :if-does-not-exist :create)
      (format table "name orientation n nnzG analyse_ms refill_ms refill_bytes~%")
      (dolist (path paths)
        (let* ((file (file-namestring path))
               (name (subseq file 0 (- (length file) (length ".J.txt")))))
          (incf count)
          (dolist (suffix '(".G.out.txt" ".G2.out.txt" ".refused.out.txt"))
            (let ((stale (fixture-path output name suffix)))
              (when (probe-file stale) (delete-file stale))))
          (handler-case
              (multiple-value-bind (m n rows columns vals) (read-rectangular input name)
                (let* ((orientation (gram-orientation input name))
                       (matrix (csc:from-coo m n rows columns vals
                                             (coalton:coalton coalton-linalg/sparse/csc:RejectDuplicates)))
                       (analyze (if (eq orientation ':ata) #'gram:analyze-ata #'gram:analyze-aat))
                       (plan (funcall analyze matrix))
                       (source (csc:entries matrix))
                       (first-out (darray (gram:entry-count plan)))
                       (second-out (darray (gram:entry-count plan)))
                       (scaled (map '(simple-array double-float (*)) (lambda (v) (* -0.5d0 v)) source))
                       (tag (if (eq orientation ':ata) "JtJ" "JJt")))
                  (gram:refill! plan scaled second-out)
                  (gram:refill! plan source first-out)
                  ;; Scaling every input by -1/2 is exact in binary, so each product is exactly a quarter.
                  (unless (every (lambda (a b) (= b (* 0.25d0 a))) first-out second-out)
                    (incf failures)
                    (format t "GRAM-SCALING-MISMATCH ~A~%" name))
                  (write-gram-file (fixture-path output name ".G.out.txt") plan first-out)
                  (write-gram-file (fixture-path output name ".G2.out.txt") plan second-out)
                  (let* ((analyse-ms (best-ms (lambda () (funcall analyze matrix)) repeats))
                         (refill-thunk (lambda () (gram:refill! plan source first-out)))
                         (refill-ms (best-ms refill-thunk repeats))
                         (bytes (bytes-for refill-thunk (max 1000 allocation-repeats)))
                         (line (format nil "~A ~A ~D ~D ~,6F ~,6F ~D" name tag (gram:dimension plan)
                                       (gram:entry-count plan) analyse-ms refill-ms bytes)))
                    (unless (zerop bytes) (incf failures))
                    (write-line line) (write-line line table))))
            (error (condition)
              (cond
                ((and linalg-error (typep condition linalg-error))
                 (let ((line (format nil "~A refused ~A" name (refusal-name condition))))
                   (with-open-file (out (fixture-path output name ".refused.out.txt")
                                        :direction :output :if-exists :supersede :if-does-not-exist :create)
                     (write-line (refusal-name condition) out))
                   (write-line line) (write-line line table)))
                (t
                 (incf failures)
                 (format t "GRAM-DRIVER-ERROR ~A ~A~%" name condition))))))))
    (if (zerop failures)
        (format t "GRAM-FIXTURES-PASS ~D~%" count)
        (error "GRAM-FIXTURES-FAIL ~D of ~D" failures count))))
