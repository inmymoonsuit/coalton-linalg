(in-package #:coalton-linalg/tests)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "GRAM" "COALTON-LINALG/SPARSE/GRAM"))

(define-condition gram-property-mismatch (simple-error) ())

(defun gram-require (test message &rest arguments)
  (unless test
    (error 'gram-property-mismatch
           :format-control message :format-arguments arguments)))

(defun gram-reference-pattern (m n entries ata-p)
  "Intersect explicit coordinate sets, independently of production incidence/group traversal."
  (let* ((dimension (if ata-p n m))
         (inner (if ata-p m n))
         (present (make-array (list m n) :element-type 'bit :initial-element 0))
         (pointers (make-array (1+ dimension) :element-type '(unsigned-byte 62)))
         (rows nil)
         (count 0))
    (dolist (entry entries)
      (setf (aref present (first entry) (second entry)) 1))
    (dotimes (column dimension)
      (setf (aref pointers column) count)
      (loop for row from column below dimension
            when (or (= row column)
                     (loop for k below inner
                           thereis (if ata-p
                                       (and (= 1 (aref present k row))
                                            (= 1 (aref present k column)))
                                       (and (= 1 (aref present row k))
                                            (= 1 (aref present column k))))))
            do (push row rows)
               (incf count)))
    (setf (aref pointers dimension) count)
    (values dimension pointers
            (make-array count :element-type '(unsigned-byte 62)
                        :initial-contents (nreverse rows)))))

(defun gram-reference-values (m n entries ata-p pointers rows)
  "Use exact rational products of the supplied F64 entries, then round once per stored value."
  (let* ((dense (primitive-dense m n entries))
         (inner (if ata-p m n))
         (dimension (if ata-p n m))
         (answer (make-array (length rows) :element-type 'double-float)))
    (dotimes (column dimension answer)
      (loop for position from (aref pointers column) below (aref pointers (1+ column))
            for row = (aref rows position)
            do (let ((total 0))
                 (dotimes (k inner)
                   (incf total
                         (if ata-p
                             (* (rational (aref dense k row))
                                (rational (aref dense k column)))
                             (* (rational (aref dense row k))
                                (rational (aref dense column k))))))
                 (setf (aref answer position) (coerce total 'double-float)))))))

(defun gram-check-pattern-data (label dimension pointers rows count expected-dimension expected-pointers expected-rows)
  "Check canonical lower CSC including all diagonals and every structural intersection."
  (gram-require (= dimension expected-dimension) "~A: Gram dimension differs." label)
  (gram-require (= (length pointers) (1+ dimension)) "~A: pointer length differs." label)
  (gram-require (= (aref pointers 0) 0) "~A: first pointer is nonzero." label)
  (gram-require (= (aref pointers dimension) (length rows) count)
                "~A: terminal pointer or entry count differs." label)
  (dotimes (column dimension)
    (let ((begin (aref pointers column))
          (end (aref pointers (1+ column))))
      (gram-require (<= 0 begin end (length rows)) "~A: pointers are out of range." label)
      (gram-require (< begin end) "~A: a diagonal is missing." label)
      (gram-require (= (aref rows begin) column) "~A: a column does not start at its diagonal." label)
      (loop for position from begin below end
            for row = (aref rows position)
            do (gram-require (<= column row (1- dimension)) "~A: row is outside the lower triangle." label)
               (when (> position begin)
                 (gram-require (< (aref rows (1- position)) row) "~A: rows are not strictly increasing." label)))))
  (gram-require (equalp pointers expected-pointers) "~A: pointer pattern differs from coordinate intersections." label)
  (gram-require (equalp rows expected-rows) "~A: row pattern differs from coordinate intersections." label))

(defun gram-check-plan (label plan dimension pointers rows)
  (gram-check-pattern-data label (gram:dimension plan) (gram:column-pointers plan)
                           (gram:row-indices plan) (gram:entry-count plan)
                           dimension pointers rows))

(defun gram-check-refill (label plan input output matrix expected operation)
  "Observe caller output and preserve source bits and the plan's borrowed structural arrays."
  (let ((input-before (copy-seq input))
        (plan-pointers (copy-seq (gram:column-pointers plan)))
        (plan-rows (copy-seq (gram:row-indices plan)))
        (source-pointers (copy-seq (csc:column-pointers matrix)))
        (source-rows (copy-seq (csc:row-indices matrix)))
        (source-values (copy-seq (csc:entries matrix))))
    (funcall operation)
    (gram-require (= (length output) (length expected)) "~A: output length differs." label)
    (dotimes (i (length output))
      (gram-require (and (primitive-finite-p (aref output i))
                        (= (aref output i) (aref expected i)))
                    "~A: value ~D expected ~S, got ~S." label i (aref expected i) (aref output i)))
    (gram-require (same-bits-p input input-before) "~A: input values changed." label)
    (gram-require (equalp plan-pointers (gram:column-pointers plan)) "~A: plan pointers changed." label)
    (gram-require (equalp plan-rows (gram:row-indices plan)) "~A: plan rows changed." label)
    (gram-require (equalp source-pointers (csc:column-pointers matrix)) "~A: source pointers changed." label)
    (gram-require (equalp source-rows (csc:row-indices matrix)) "~A: source rows changed." label)
    (gram-require (same-bits-p source-values (csc:entries matrix)) "~A: source values changed." label))
  output)

(defun gram-case-values (entries pass seed)
  "Keep coordinates fixed while replacing values, including a complete zero generation."
  (loop for (row column value) in entries
        for index from 0
        collect (list row column
                      (case pass
                        (0 value)
                        (1 0d0)
                        (2 (/ (coerce (- (mod (+ seed (* 7 index)) 17) 8) 'double-float) 8d0))))))

(defun gram-explicit-cases ()
  (list
   (list "empty" 0 0 nil)
   (list "no-rows" 0 5 nil)
   (list "no-columns" 4 0 nil)
   (list "no-entries" 3 4 nil)
   (list "stored-zeros" 3 3 '((0 0 0d0) (0 2 -0d0) (2 1 0d0)))
   (list "diagonal" 3 3 '((0 0 2d0) (1 1 -3d0) (2 2 0.5d0)))
   (list "rectangle" 2 3 '((0 0 2d0) (1 0 0d0) (1 1 3d0) (0 2 -1d0) (1 2 4d0)))
   (list "disjoint-columns" 5 3 '((0 0 1d0) (4 0 2d0) (1 1 3d0) (3 2 -1d0)))
   (list "exact-cancellation" 2 2 '((0 0 1d0) (1 0 1d0) (0 1 1d0) (1 1 -1d0)))
   (list "wide-gaps" 7 5 '((0 0 2d0) (6 0 -1d0) (0 2 3d0) (4 2 0d0) (6 4 0.25d0)))
   (list "dense-ones" 4 3 (loop for row below 4 append (loop for column below 3 collect (list row column 1d0))))
   (list "exact-dependent" 3 3 '((0 0 1d0) (0 1 1d0) (0 2 2d0)
                                (1 0 1d0) (1 2 1d0) (2 0 1d0) (2 2 1d0)))))

(defun gram-seeded-case (index seed)
  "Small dyadic coefficients keep each Gram product and partial sum exactly representable."
  (let* ((m (1+ (mod (* 7 index) 13)))
         (n (1+ (mod (+ (* 5 index) 2) 15)))
         (family (mod index 5))
         (exponent (nth (mod (floor index 5) 3) '(-160 0 160)))
         (state (mod (+ seed (* 7919 index)) (ash 1 32)))
         (entries nil))
    (labels ((next (bound)
               (setf state (mod (+ (* 1664525 state) 1013904223) (ash 1 32)))
               (mod state bound)))
      (dotimes (row m)
        (dotimes (column n)
          (when (case family
                  (0 t)
                  (1 (<= (abs (- row column)) 1))
                  (2 (zerop (mod (+ row column) 3)))
                  (3 (and (not (zerop (mod row 4))) (not (zerop (mod column 3))) (< (next 5) 3)))
                  (4 (< (next 7) 2)))
            (push (list row column
                        (scale-float (/ (coerce (- (next 33) 16) 'double-float) 8d0) exponent))
                  entries)))))
    (list (format nil "seeded-~2,'0D" index) m n (nreverse entries))))

(defun gram-run-exact-case (case seed)
  "Reuse both plans and their outputs for three value generations through both refill entry points."
  (destructuring-bind (name m n entries) case
    (let* ((matrix (primitive-from-entries m n (reverse entries) (coalton:coalton csc:RejectDuplicates)))
           (original-values (copy-seq (csc:entries matrix)))
           (original-pointers (copy-seq (csc:column-pointers matrix)))
           (original-rows (copy-seq (csc:row-indices matrix))))
      (dolist (ata-p '(t nil))
        (multiple-value-bind (dimension pointers rows) (gram-reference-pattern m n entries ata-p)
          (let* ((plan (if ata-p (gram:analyze-ata matrix) (gram:analyze-aat matrix)))
                 (limited (if ata-p
                              (gram:analyze-ata-with-limit matrix (length rows))
                              (gram:analyze-aat-with-limit matrix (length rows))))
                 (output (make-array (length rows) :element-type 'double-float
                                     :initial-element (sb-kernel:make-double-float #x7ff80000 53)))
                 (label (format nil "~A/~A" name (if ata-p "ATA" "AAT"))))
            (gram-check-plan label plan dimension pointers rows)
            (gram-check-plan (concatenate 'string label "/exact-budget") limited dimension pointers rows)
            (dotimes (pass 3)
              (let* ((current-entries (gram-case-values entries pass seed))
                     (input (make-array (length entries) :element-type 'double-float
                                        :initial-contents (mapcar #'third (primitive-sort-entries current-entries))))
                     (current-matrix (primitive-from-entries
                                      m n (reverse current-entries) (coalton:coalton csc:RejectDuplicates)))
                     (expected (gram-reference-values m n current-entries ata-p pointers rows)))
                (gram-check-refill (format nil "~A/pass~D/raw" label pass)
                                   plan input output current-matrix expected
                                   (lambda () (gram:refill! plan input output)))
                (gram-check-refill (format nil "~A/pass~D/checked" label pass)
                                   plan input output current-matrix expected
                                   (lambda () (gram:refill-csc! plan current-matrix output)))
                (gram-check-plan label plan dimension pointers rows)))
            (gram-require (same-bits-p original-values (csc:entries matrix))
                          "~A: analysis source values changed during reuse." label)
            (gram-require (equalp original-pointers (csc:column-pointers matrix))
                          "~A: analysis source pointers changed during reuse." label)
            (gram-require (equalp original-rows (csc:row-indices matrix))
                          "~A: analysis source rows changed during reuse." label)))))))

(defun gram-rounded-source (ata-p)
  "Arrange an off-diagonal sum as 2^54,1,-2^54 in increasing contracted index."
  (let* ((large (scale-float 1d0 54))
         (entries (if ata-p
                      (list '(0 0 1d0) '(1 0 1d0) '(2 0 1d0)
                            (list 0 1 large) '(1 1 1d0) (list 2 1 (- large)))
                      (list '(0 0 1d0) '(0 1 1d0) '(0 2 1d0)
                            (list 1 0 large) '(1 1 1d0) (list 1 2 (- large)))))
         (m (if ata-p 3 2))
         (n (if ata-p 2 3)))
    (values (primitive-from-entries m n entries (coalton:coalton csc:RejectDuplicates))
            m n entries)))

(defun gram-rounding-tests ()
  (dolist (ata-p '(t nil))
    (multiple-value-bind (matrix m n entries) (gram-rounded-source ata-p)
      (let* ((plan (if ata-p (gram:analyze-ata matrix) (gram:analyze-aat matrix)))
             (input (copy-seq (csc:entries matrix)))
             (output (primitive-doubles 913d0 913d0 913d0))
             (pointers (primitive-indices 0 2 3))
             (rows (primitive-indices 0 1 1))
             (expected (primitive-doubles 3d0 0d0 (scale-float 1d0 109)))
             (exact (gram-reference-values m n entries ata-p pointers rows)))
        (gram-require (= 1d0 (aref exact 1)) "Rounded-order reference lost its exact unit off-diagonal.")
        (gram-check-plan "rounded Gram pattern" plan 2 pointers rows)
        (gram-check-refill "rounded Gram raw refill" plan input output matrix expected
                           (lambda () (gram:refill! plan input output)))
        (gram-check-refill "rounded Gram checked refill" plan input output matrix expected
                           (lambda () (gram:refill-csc! plan matrix output))))))
  (format t "GRAM-ROUNDED-ORDER-PASS orientations=2 refills=4 exact-offdiagonal=1 computed=0~%"))

(defun gram-control-rejected (label operation)
  "Count only a mismatch detected by this oracle, not arbitrary runtime/API failures."
  (let ((rejected nil))
    (handler-case (funcall operation)
      (gram-property-mismatch () (setf rejected t)))
    (assert rejected () "Gram oracle missed planted defect ~A." label)))

(defun gram-oracle-controls ()
  (let* ((entries '((0 0 0d0) (0 1 0d0)))
         (matrix (primitive-from-entries 1 2 entries (coalton:coalton csc:RejectDuplicates)))
         (plan (gram:analyze-ata matrix))
         (expected-pointers (primitive-indices 0 2 3))
         (expected-rows (primitive-indices 0 1 1)))
    (gram-check-plan "structural control baseline" plan 2 expected-pointers expected-rows)
    (dolist (control (list (list "missing mandatory diagonal" (primitive-indices 0 1 1) (primitive-indices 0))
                          (list "pruned structural zero" (primitive-indices 0 1 2) (primitive-indices 0 1))
                          (list "upper instead of lower triangle" (primitive-indices 0 1 3) (primitive-indices 0 0 1))))
      (destructuring-bind (name pointers rows) control
        (gram-control-rejected name
          (lambda () (gram-check-pattern-data name 2 pointers rows (length rows)
                                              2 expected-pointers expected-rows))))))
  (let* ((entries '((0 0 1d0) (1 0 3d0) (0 1 2d0) (1 1 4d0)))
         (matrix (primitive-from-entries 2 2 entries (coalton:coalton csc:RejectDuplicates)))
         (plan (gram:analyze-ata matrix))
         (other (gram:analyze-aat matrix))
         (expected (primitive-doubles 10d0 14d0 20d0)))
    (dolist (name '("missing product" "wrong orientation" "retained stale output" "source mutation"))
      (let ((input (copy-seq (csc:entries matrix)))
            (output (primitive-doubles 913d0 913d0 913d0)))
        (gram-check-refill "numeric control baseline" plan input output matrix expected
                           (lambda () (gram:refill! plan input output)))
        (gram-control-rejected name
          (lambda ()
            (gram-check-refill name plan input output matrix expected
              (lambda ()
                (gram:refill! (if (string= name "wrong orientation") other plan) input output)
                (cond ((string= name "missing product") (decf (aref output 0) 1d0))
                      ((string= name "retained stale output") (incf (aref output 0) 913d0))
                      ((string= name "source mutation") (incf (aref input 0) 1d0))))))))))
  (multiple-value-bind (matrix) (gram-rounded-source t)
    (let* ((plan (gram:analyze-ata matrix))
           (input (copy-seq (csc:entries matrix)))
           (output (primitive-doubles 913d0 913d0 913d0))
           (expected (primitive-doubles 3d0 0d0 (scale-float 1d0 109))))
      (gram-check-refill "rounding control baseline" plan input output matrix expected
                         (lambda () (gram:refill! plan input output)))
      (gram-control-rejected "reordered accumulation"
        (lambda ()
          (gram-check-refill "reordered accumulation" plan input output matrix expected
            (lambda ()
              (gram:refill! plan input output)
              (setf (aref output 1) (+ (+ (aref input 3) (aref input 5)) (aref input 4)))))))))
  (format t "GRAM-ORACLE-CONTROLS-PASS rejected=8~%"))

(defun run-gram-property-tests ()
  "Independent structural intersections, rational products and ordered-F64 controls."
  (let* ((seed 20261119)
         (cases (append (gram-explicit-cases)
                        (loop for index below 64 collect (gram-seeded-case index seed)))))
    (loop for case in cases
          for ordinal from 0
          do (gram-run-exact-case case (+ seed (* 104729 ordinal))))
    (format t "GRAM-EXACT-PROPERTIES-PASS seed=~D cases=~D orientations=~D refills=~D~%"
            seed (length cases) (* 2 (length cases)) (* 12 (length cases))))
  (gram-rounding-tests)
  (gram-oracle-controls)
  t)
