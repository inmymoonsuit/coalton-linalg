(in-package #:coalton-linalg/tests)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "CSC" "COALTON-LINALG/SPARSE/CSC")
  (sb-ext:add-package-local-nickname "NUMERIC" "COALTON-LINALG/NUMERIC")
  (sb-ext:add-package-local-nickname "LDLT" "COALTON-LINALG/SPARSE/SPD-LDLT"))

(defun primitive-doubles (&rest values)
  (make-array (length values) :element-type 'double-float
              :initial-contents (mapcar (lambda (value) (coerce value 'double-float)) values)))

(defun primitive-indices (&rest values)
  (make-array (length values) :element-type '(unsigned-byte 62) :initial-contents values))

(defun primitive-refuses (label thunk)
  (assert (handler-case (progn (funcall thunk) nil) (error () t)) ()
          "~A should signal an error." label))

(defun primitive-finite-p (value)
  (and (typep value 'double-float)
       (not (sb-ext:float-nan-p value))
       (not (sb-ext:float-infinity-p value))))

(defun primitive-close (label actual expected &optional (relative 2d-14))
  "Use a relative comparison without an absolute floor that would hide tiny zeros."
  (assert (primitive-finite-p actual) () "~A returned nonfinite ~S." label actual)
  (assert (if (zerop expected)
              (zerop actual)
              (<= (abs (- actual expected)) (* relative (abs expected))))
          () "~A: expected ~S, got ~S." label expected actual))

(defun primitive-vector (label actual expected &optional (relative 2d-14))
  (assert (= (length actual) (length expected)) () "~A length differs." label)
  (loop for a across actual for e across expected
        do (primitive-close label a e relative)))

(defun primitive-dense (m n entries)
  "Build a dense mathematical reference from explicitly supplied (row column value) entries."
  (let ((dense (make-array (list m n) :element-type 'double-float :initial-element 0d0)))
    (dolist (entry entries dense)
      (destructuring-bind (row column value) entry
        (setf (aref dense row column) value)))))

(defun primitive-product (dense vector transpose-p)
  "Accumulate exact rational products, then round once; small test cases are exactly representable."
  (let* ((m (array-dimension dense 0))
         (n (array-dimension dense 1))
         (result (make-array (if transpose-p n m) :element-type 'double-float)))
    (dotimes (i (length result) result)
      (let ((sum 0))
        (dotimes (j (if transpose-p m n))
          (incf sum (* (rational (if transpose-p (aref dense j i) (aref dense i j)))
                       (rational (aref vector j)))))
        (setf (aref result i) (coerce sum 'double-float))))))

(defun primitive-dense-permute (dense row-order column-order)
  (let ((answer (make-array (list (length row-order) (length column-order))
                            :element-type 'double-float)))
    (dotimes (i (length row-order) answer)
      (dotimes (j (length column-order))
        (setf (aref answer i j) (aref dense (aref row-order i) (aref column-order j)))))))

(defun primitive-sort-entries (entries)
  (sort (copy-list entries)
        (lambda (a b) (or (< (second a) (second b))
                         (and (= (second a) (second b)) (< (first a) (first b)))))))

(defun primitive-matrix (label matrix dense entries)
  "Check every dense coordinate and every stored entry, retaining structural zeros."
  (let* ((m (array-dimension dense 0))
         (n (array-dimension dense 1))
         (pointers (csc:column-pointers matrix))
         (rows (csc:row-indices matrix))
         (values (csc:entries matrix))
         (ordered (primitive-sort-entries entries))
         (actual (make-array (list m n) :element-type 'double-float :initial-element 0d0)))
    (assert (= m (csc:row-count matrix)) () "~A row count differs." label)
    (assert (= n (csc:column-count matrix)) () "~A column count differs." label)
    (assert (= (length ordered) (csc:entry-count matrix) (length rows) (length values))
            () "~A stored entry count differs." label)
    (assert (= (1+ n) (length pointers)) () "~A pointer length differs." label)
    (assert (= 0 (aref pointers 0)) () "~A starts at nonzero pointer." label)
    (assert (= (length rows) (aref pointers n)) () "~A final pointer differs." label)
    (let ((position 0))
      (dotimes (column n)
        (assert (<= (aref pointers column) (aref pointers (1+ column)) (length rows)))
        (loop with previous = -1
              for q from (aref pointers column) below (aref pointers (1+ column))
              for row = (aref rows q)
              do (assert (< previous row m) () "~A is not canonical CSC." label)
                 (setf previous row)
                 (destructuring-bind (er ec ev) (nth position ordered)
                   (assert (= row er) () "~A stored row differs." label)
                   (assert (= column ec) () "~A stored column differs." label)
                   (primitive-close label (aref values q) ev 0d0))
                 (setf (aref actual row column) (aref values q))
                 (incf position)))
      (assert (= position (length ordered))))
    (dotimes (i m)
      (dotimes (j n)
        (primitive-close label (aref actual i j) (aref dense i j) 0d0)))
    matrix))

(defun primitive-from-entries (m n entries policy)
  (csc:from-coo m n
                (apply #'primitive-indices (mapcar #'first entries))
                (apply #'primitive-indices (mapcar #'second entries))
                (apply #'primitive-doubles (mapcar #'third entries))
                policy))

(defun primitive-fresh-storage (label original derived)
  (assert (not (eq (csc:column-pointers original) (csc:column-pointers derived))) ()
          "~A shares column storage." label)
  (assert (not (eq (csc:row-indices original) (csc:row-indices derived))) ()
          "~A shares row storage." label)
  (assert (not (eq (csc:entries original) (csc:entries derived))) ()
          "~A shares value storage." label))

(defun primitive-norm-tests ()
  (dolist (value (list 0d0 -0d0 3d0 -9d0 least-positive-double-float
                       most-positive-double-float))
    (assert (numeric:finite? value)))
  (let ((nan (sb-kernel:make-double-float #x7ff80000 0)))
    (dolist (value (list nan sb-ext:double-float-positive-infinity
                         sb-ext:double-float-negative-infinity))
      (assert (not (numeric:finite? value)))
      (primitive-refuses "nonfinite norm2"
                         (lambda () (numeric:norm2 (primitive-doubles 0d0 value))))
      (primitive-refuses "nonfinite norm-inf"
                         (lambda () (numeric:norm-inf (primitive-doubles 0d0 value))))))
  (dolist (case (list (list (primitive-doubles) 0d0 0d0)
                     (list (primitive-doubles 0d0 -0d0 0d0) 0d0 0d0)
                     (list (primitive-doubles 3d0 -4d0) 5d0 4d0)
                     (list (primitive-doubles 3d200 -4d200) 5d200 4d200)
                     (list (primitive-doubles 3d-200 -4d-200) 5d-200 4d-200)
                     (list (primitive-doubles 1d200 -2d-200 0d0) 1d200 1d200)
                     (list (primitive-doubles least-positive-double-float)
                           least-positive-double-float least-positive-double-float)
                     (list (primitive-doubles (- most-positive-double-float))
                           most-positive-double-float most-positive-double-float)))
    (destructuring-bind (vector euclidean maximum) case
      (let ((before (copy-seq vector)))
        (primitive-close "norm2" (numeric:norm2 vector) euclidean)
        (primitive-close "norm-inf" (numeric:norm-inf vector) maximum 0d0)
        (assert (same-bits-p vector before))
        (primitive-close "norm2 reversed" (numeric:norm2 (reverse vector)) euclidean))))
  (primitive-refuses "unrepresentable norm2"
                     (lambda () (numeric:norm2
                                 (primitive-doubles most-positive-double-float
                                                    most-positive-double-float))))
  (format t "PRIMITIVE-NORMS-PASS~%"))

(defun primitive-generated-case ()
  "A deterministic rectangular roster with stored zeros and uneven merge lengths."
  (list 7 11
        (loop for row below 7 append
              (loop for column below 11
                    when (< (mod (+ (* 11 row) (* 7 column)) 5) 2)
                    collect (list row column
                                  (coerce (- (mod (+ row (* 3 column)) 5) 2) 'double-float))))))

(defun primitive-construction-tests ()
  (let* ((pointers (primitive-indices 0 2 2 5 6))
         (rows (primitive-indices 0 2 0 1 2 1))
         (values (primitive-doubles 2d0 -1d0 0d0 3d0 4d0 -2d0))
         (entries '((0 0 2d0) (2 0 -1d0) (0 2 0d0) (1 2 3d0) (2 2 4d0) (1 3 -2d0)))
         (dense (primitive-dense 3 4 entries))
         (matrix (csc:from-csc 3 4 pointers rows values)))
    (assert (not (eq pointers (csc:column-pointers matrix))))
    (assert (not (eq rows (csc:row-indices matrix))))
    (assert (not (eq values (csc:entries matrix))))
    (assert (equalp pointers (primitive-indices 0 2 2 5 6)))
    (assert (equalp rows (primitive-indices 0 2 0 1 2 1)))
    (assert (equalp values (primitive-doubles 2d0 -1d0 0d0 3d0 4d0 -2d0)))
    (fill pointers 0) (fill rows 0) (fill values 99d0)
    (primitive-matrix "owned CSC" matrix dense entries)
    (primitive-matrix "unsorted unique COO"
                      (primitive-from-entries 3 4 (reverse entries)
                                              (coalton:coalton csc:RejectDuplicates))
                      dense entries)
    (let* ((coo-rows (primitive-indices 1 0 1))
           (coo-columns (primitive-indices 1 0 1))
           (coo-values (primitive-doubles 3d0 7d0 -3d0))
           (sum (csc:from-coo 2 2 coo-rows coo-columns coo-values
                              (coalton:coalton csc:SumDuplicates))))
      (assert (equalp coo-values (primitive-doubles 3d0 7d0 -3d0)))
      (assert (equalp coo-rows (primitive-indices 1 0 1)))
      (assert (equalp coo-columns (primitive-indices 1 0 1)))
      (fill coo-rows 0) (fill coo-columns 0) (fill coo-values 99d0)
      (primitive-matrix "duplicate cancellation retains zero" sum
                        (primitive-dense 2 2 '((0 0 7d0) (1 1 0d0)))
                        '((0 0 7d0) (1 1 0d0)))))
  (let* ((large (scale-float 1d0 54))
         (policy (coalton:coalton csc:SumDuplicates))
         (lost (primitive-from-entries 1 1 (list (list 0 0 large)
                                                 (list 0 0 1d0) (list 0 0 (- large))) policy))
         (kept (primitive-from-entries 1 1 (list (list 0 0 large)
                                                 (list 0 0 (- large)) (list 0 0 1d0)) policy)))
    (primitive-matrix "duplicate input order rounds away one" lost
                      (primitive-dense 1 1 '((0 0 0d0))) '((0 0 0d0)))
    (primitive-matrix "duplicate input order retains one" kept
                      (primitive-dense 1 1 '((0 0 1d0))) '((0 0 1d0))))
  (flet ((bad-csc (m n pointers rows values)
           (primitive-refuses "malformed CSC"
                              (lambda () (csc:from-csc m n pointers rows values)))))
    (bad-csc 2 2 (primitive-indices 0 0) (primitive-indices) (primitive-doubles))
    (bad-csc 2 2 (primitive-indices 1 1 1) (primitive-indices 0) (primitive-doubles 1d0))
    (bad-csc 2 2 (primitive-indices 0 2 1) (primitive-indices 0) (primitive-doubles 1d0))
    (bad-csc 2 2 (primitive-indices 0 0 0) (primitive-indices 0) (primitive-doubles 1d0))
    (bad-csc 2 1 (primitive-indices 0 1) (primitive-indices 2) (primitive-doubles 1d0))
    (bad-csc 2 1 (primitive-indices 0 2) (primitive-indices 1 0) (primitive-doubles 1d0 2d0))
    (bad-csc 2 1 (primitive-indices 0 2) (primitive-indices 0 0) (primitive-doubles 1d0 2d0))
    (bad-csc 2 1 (primitive-indices 0 1) (primitive-indices 0) (primitive-doubles))
    (bad-csc 0 1 (primitive-indices 0 1) (primitive-indices 0) (primitive-doubles 1d0))
    (dolist (bad (list sb-ext:double-float-positive-infinity
                       sb-ext:double-float-negative-infinity
                       (sb-kernel:make-double-float #x7ff80000 0)))
      (bad-csc 1 1 (primitive-indices 0 1) (primitive-indices 0) (primitive-doubles bad))))
  (dolist (policy (list (coalton:coalton csc:RejectDuplicates)
                        (coalton:coalton csc:SumDuplicates)))
    (flet ((bad-coo (m n rows columns values)
             (primitive-refuses "malformed COO"
                                (lambda () (csc:from-coo m n rows columns values policy)))))
      (bad-coo 2 2 (primitive-indices 0) (primitive-indices) (primitive-doubles 1d0))
      (bad-coo 2 2 (primitive-indices 0) (primitive-indices 0) (primitive-doubles))
      (bad-coo 2 2 (primitive-indices 2) (primitive-indices 0) (primitive-doubles 1d0))
      (bad-coo 2 2 (primitive-indices 0) (primitive-indices 2) (primitive-doubles 1d0))
      (bad-coo 0 0 (primitive-indices 0) (primitive-indices 0) (primitive-doubles 1d0))
      (dolist (bad (list sb-ext:double-float-positive-infinity
                         sb-ext:double-float-negative-infinity
                         (sb-kernel:make-double-float #x7ff80000 0)))
        (bad-coo 1 1 (primitive-indices 0) (primitive-indices 0) (primitive-doubles bad)))))
  (primitive-refuses "duplicate zero rejected"
                     (lambda () (primitive-from-entries 1 1 '((0 0 0d0) (0 0 0d0))
                                                        (coalton:coalton csc:RejectDuplicates))))
  (primitive-refuses "duplicate sum intermediate overflow"
                     (lambda () (primitive-from-entries
                                 1 1 (list (list 0 0 most-positive-double-float)
                                           (list 0 0 most-positive-double-float)
                                           (list 0 0 (- most-positive-double-float)))
                                 (coalton:coalton csc:SumDuplicates))))
  (destructuring-bind (m n entries) (primitive-generated-case)
    (let* ((raw (loop for (row column value) in entries append
                     (list (list row column 7d0) (list row column value) (list row column -7d0))))
           (scrambled (append (nthcdr 17 raw) (subseq raw 0 17)))
           (matrix (primitive-from-entries m n scrambled (coalton:coalton csc:SumDuplicates))))
      (primitive-matrix "many interleaved duplicates" matrix (primitive-dense m n entries) entries)
      (primitive-refuses "many interleaved duplicates rejected"
                         (lambda () (primitive-from-entries m n scrambled
                                                            (coalton:coalton csc:RejectDuplicates))))))
  (format t "PRIMITIVE-CONSTRUCTION-PASS~%"))

(defun primitive-transform-tests ()
  (dolist (case (append '((0 0 ()) (0 3 ()) (4 0 ()) (2 3 ())
                 (1 1 ((0 0 0d0)))
                 (3 4 ((0 0 2d0) (2 0 -1d0) (0 2 0d0)
                       (1 2 3d0) (2 2 4d0) (1 3 -2d0))))
                       (list (primitive-generated-case))))
    (destructuring-bind (m n entries) case
      (let* ((matrix (primitive-from-entries m n (reverse entries)
                                             (coalton:coalton csc:RejectDuplicates)))
             (dense (primitive-dense m n entries))
             (x (apply #'primitive-doubles (loop for j below n collect (- (* 2 j) 3))))
             (y (apply #'primitive-doubles (loop for i below m collect (+ i 2))))
             (x-before (copy-seq x))
             (y-before (copy-seq y))
             (ax (csc:matvec matrix x))
             (aty (csc:transpose-matvec matrix y))
             (transposed (csc:transpose matrix))
             (transposed-entries (mapcar (lambda (entry)
                                          (list (second entry) (first entry) (third entry)))
                                        entries))
             (rows (apply #'primitive-indices
                          (if (= m 3) '(2 0 1) (reverse (loop for i below m collect i)))))
             (columns (apply #'primitive-indices
                             (if (= n 4) '(2 0 3 1) (reverse (loop for j below n collect j)))))
             (rows-before (copy-seq rows))
             (columns-before (copy-seq columns))
             (permuted (csc:permute matrix rows columns))
             (permuted-entries (mapcar (lambda (entry)
                                        (list (position (first entry) rows)
                                              (position (second entry) columns) (third entry)))
                                      entries)))
        (primitive-matrix "constructed matrix" matrix dense entries)
        (let ((roundtrip (csc:from-csc m n (csc:column-pointers matrix)
                                        (csc:row-indices matrix) (csc:entries matrix))))
          (primitive-matrix "CSC roundtrip including empty dimensions" roundtrip dense entries)
          (primitive-fresh-storage "CSC roundtrip" matrix roundtrip))
        (primitive-vector "rectangular Ax" ax (primitive-product dense x nil) 0d0)
        (primitive-vector "rectangular transpose Ax" aty (primitive-product dense y t) 0d0)
        (assert (same-bits-p x x-before))
        (assert (same-bits-p y y-before))
        (assert (not (eq x ax)))
        (assert (not (eq y aty)))
        (assert (not (eq ax (csc:matvec matrix x))))
        (primitive-matrix "transpose" transposed (primitive-dense n m transposed-entries)
                          transposed-entries)
        (primitive-fresh-storage "transpose" matrix transposed)
        (primitive-matrix "double transpose" (csc:transpose transposed) dense entries)
        (primitive-vector "transpose matvec identity" (csc:matvec transposed y) aty 0d0)
        (primitive-vector "transpose transpose-matvec identity"
                          (csc:transpose-matvec transposed x) ax 0d0)
        (primitive-matrix "new to old permutation" permuted
                          (primitive-dense-permute dense rows columns) permuted-entries)
        (primitive-fresh-storage "permutation" matrix permuted)
        (assert (equalp rows rows-before))
        (assert (equalp columns columns-before))
        (fill rows 0) (fill columns 0)
        (primitive-matrix "permutation owns its maps" permuted
                          (primitive-dense-permute dense rows-before columns-before) permuted-entries)
        (primitive-matrix "transformations preserve source" matrix dense entries)
        (primitive-refuses "Ax wrong length"
                           (lambda () (csc:matvec matrix (make-array (1+ n)
                                            :element-type 'double-float :initial-element 0d0))))
        (primitive-refuses "transpose Ax wrong length"
                           (lambda () (csc:transpose-matvec matrix (make-array (1+ m)
                                            :element-type 'double-float :initial-element 0d0)))))))
  (let ((matrix (primitive-from-entries 3 4 '((0 0 2d0) (2 2 3d0))
                                        (coalton:coalton csc:RejectDuplicates))))
    (dolist (rows (list (primitive-indices 0 1) (primitive-indices 0 1 1)
                        (primitive-indices 0 1 3)))
      (primitive-refuses "invalid row permutation"
                         (lambda () (csc:permute matrix rows (primitive-indices 0 1 2 3)))))
    (dolist (columns (list (primitive-indices 0 1 2) (primitive-indices 0 1 2 2)
                           (primitive-indices 0 1 2 4)))
      (primitive-refuses "invalid column permutation"
                         (lambda () (csc:permute matrix (primitive-indices 0 1 2) columns)))))
  (let* ((matrix (primitive-from-entries 2 2 '((0 0 2d0) (1 1 3d0))
                                         (coalton:coalton csc:RejectDuplicates)))
         (input (csc:entries matrix))
         (before (copy-seq input)))
    (primitive-vector "matrix values used as Ax input" (csc:matvec matrix input)
                      (primitive-doubles 4d0 9d0) 0d0)
    (primitive-vector "matrix values used as transpose Ax input" (csc:transpose-matvec matrix input)
                      (primitive-doubles 4d0 9d0) 0d0)
    (assert (same-bits-p input before)))
  (format t "PRIMITIVE-TRANSFORMS-PASS~%"))

(defun primitive-matvec-refusal-tests ()
  (let ((matrix (primitive-from-entries 2 3 '((0 0 1d0))
                                        (coalton:coalton csc:RejectDuplicates))))
    (dolist (bad (list sb-ext:double-float-positive-infinity
                       sb-ext:double-float-negative-infinity
                       (sb-kernel:make-double-float #x7ff80000 0)))
      (primitive-refuses "unreferenced nonfinite Ax input"
                         (lambda () (csc:matvec matrix (primitive-doubles 1d0 0d0 bad))))
      (primitive-refuses "unreferenced nonfinite transpose Ax input"
                         (lambda () (csc:transpose-matvec matrix (primitive-doubles 1d0 bad))))))
  (dolist (transpose-p '(nil t))
    (let ((matrix (primitive-from-entries
                   (if transpose-p 2 1) (if transpose-p 1 2)
                   (if transpose-p
                       (list (list 0 0 most-positive-double-float)
                             (list 1 0 most-positive-double-float))
                       (list (list 0 0 most-positive-double-float)
                             (list 0 1 most-positive-double-float)))
                   (coalton:coalton csc:RejectDuplicates))))
      (primitive-refuses "matvec sum overflow"
                         (lambda () (if transpose-p
                                        (csc:transpose-matvec matrix (primitive-doubles 1d0 1d0))
                                        (csc:matvec matrix (primitive-doubles 1d0 1d0)))))
      (primitive-refuses "matvec product overflow"
                         (lambda () (if transpose-p
                                        (csc:transpose-matvec matrix (primitive-doubles 2d0 0d0))
                                        (csc:matvec matrix (primitive-doubles 2d0 0d0)))))))
  (format t "PRIMITIVE-MATVEC-REFUSALS-PASS~%"))

(defun primitive-ldlt-bridge-tests ()
  (let* ((matrix (primitive-from-entries 2 2 '((0 0 4d0) (1 0 1d0) (1 1 3d0))
                                         (coalton:coalton csc:RejectDuplicates)))
         (symbolic (ldlt:analyse-csc matrix))
         (again (ldlt:analyse-csc matrix))
         (result (ldlt:factor! symbolic (csc:entries matrix) 0d0))
         (out (primitive-doubles 99d0 99d0)))
    (assert (not (eq (ldlt:symbolic-factor symbolic) (ldlt:symbolic-factor again))))
    (assert (eq (coalton:coalton ldlt:FactorSuccess) (ldlt:factor-status result)))
    (assert (eq (coalton:coalton ldlt:SolveSuccess) (ldlt:solve! (ldlt:result-factor result) (primitive-doubles 6d0 7d0) out)))
    (primitive-vector "analyse-csc solve" out (primitive-doubles 1d0 2d0))
    (assert (eq (coalton:coalton ldlt:FactorSuccess) (ldlt:factor-status (ldlt:factor! again (csc:entries matrix) 0d0))))
    (assert (eq (coalton:coalton ldlt:FactorNonPositive) (ldlt:factor-status (ldlt:factor! again (csc:entries matrix) 10d0))))
    (assert (ldlt:factor-ready? (ldlt:symbolic-factor symbolic)))
    (primitive-matrix "analysis preserves matrix" matrix
                      (primitive-dense 2 2 '((0 0 4d0) (1 0 1d0) (1 1 3d0)))
                      '((0 0 4d0) (1 0 1d0) (1 1 3d0))))
  (let* ((empty (primitive-from-entries 0 0 nil (coalton:coalton csc:RejectDuplicates)))
         (symbolic (ldlt:analyse-csc empty)))
    (assert (eq (coalton:coalton ldlt:FactorSuccess) (ldlt:factor-status (ldlt:factor! symbolic (primitive-doubles) 0d0)))))
  (dolist (case '((2 3 ((0 0 1d0) (1 1 1d0)))
                 (2 2 ((0 0 2d0) (0 1 1d0) (1 1 2d0)))
                 (2 2 ((0 0 2d0) (1 0 1d0)))
                 (2 2 ((1 0 1d0) (1 1 2d0)))))
    (destructuring-bind (m n entries) case
      (primitive-refuses "analyse-csc shape or lower diagonal structure"
                         (lambda () (ldlt:analyse-csc
                                     (primitive-from-entries m n entries
                                                            (coalton:coalton csc:RejectDuplicates)))))))
  (let* ((indefinite (primitive-from-entries 1 1 '((0 0 -1d0))
                                             (coalton:coalton csc:RejectDuplicates)))
         (symbolic (ldlt:analyse-csc indefinite)))
    (assert (eq (coalton:coalton ldlt:FactorNonPositive) (ldlt:factor-status (ldlt:factor! symbolic (csc:entries indefinite) 0d0)))))
  (format t "PRIMITIVE-LDLT-BRIDGE-PASS~%"))

(defun primitive-oracle-controls ()
  (primitive-refuses "norm oracle planted zero"
                     (lambda () (primitive-close "control" 0d0 5d-200)))
  (primitive-refuses "dense product oracle planted wrong component"
                     (lambda () (primitive-vector "control" (primitive-doubles 2d0 11d0 7d0)
                                                  (primitive-doubles 2d0 12d0 7d0) 0d0)))
  (let ((matrix (primitive-from-entries 2 2 '((0 0 1d0))
                                        (coalton:coalton csc:RejectDuplicates))))
    (primitive-refuses "stored-zero oracle planted omission"
                       (lambda () (primitive-matrix "control" matrix
                                                    (primitive-dense 2 2 '((0 0 1d0) (1 1 0d0)))
                                                    '((0 0 1d0) (1 1 0d0))))))
  (let* ((entries '((0 0 1d0) (1 1 2d0) (2 2 3d0)))
         (matrix (primitive-from-entries 3 3 entries (coalton:coalton csc:RejectDuplicates)))
         (wrong (csc:permute matrix (primitive-indices 1 2 0) (primitive-indices 0 1 2))))
    (primitive-refuses "permutation oracle planted inverse direction"
                       (lambda () (primitive-matrix
                                   "control" wrong
                                   (primitive-dense 3 3 '((0 2 3d0) (1 0 1d0) (2 1 2d0)))
                                   '((0 2 3d0) (1 0 1d0) (2 1 2d0))))))
  (format t "PRIMITIVE-ORACLE-CONTROLS-PASS count=4~%"))

(defun run-primitive-tests ()
  "Exercise the public primitive interface against dense formulas and explicit contracts."
  (primitive-norm-tests)
  (primitive-construction-tests)
  (primitive-transform-tests)
  (primitive-matvec-refusal-tests)
  (primitive-ldlt-bridge-tests)
  (primitive-oracle-controls)
  (format t "PRIMITIVE-TESTS-PASS~%")
  t)
