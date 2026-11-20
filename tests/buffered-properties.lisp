(in-package #:coalton-linalg/tests)

(define-condition buffered-property-mismatch (simple-error) ())

(defun buffered-require (test message &rest arguments)
  (unless test
    (error 'buffered-property-mismatch
           :format-control message :format-arguments arguments)))

(defun buffered-public-product (transpose-p matrix input output)
  (if transpose-p
      (csc:transpose-matvec! matrix input output)
      (csc:matvec! matrix input output)))

(defun buffered-check-call (label matrix input output transpose-p expected
                            &optional (kernel #'buffered-public-product))
  "Check a public product against an independent answer and preserve all source bits."
  (let ((input-before (copy-seq input))
        (pointers-before (copy-seq (csc:column-pointers matrix)))
        (rows-before (copy-seq (csc:row-indices matrix)))
        (values-before (copy-seq (csc:entries matrix))))
    (funcall kernel transpose-p matrix input output)
    (buffered-require (= (length output) (length expected))
                      "~A: output length differs." label)
    (dotimes (i (length output))
      (buffered-require (and (primitive-finite-p (aref output i))
                            (= (aref output i) (aref expected i)))
                        "~A: output ~D expected ~S, got ~S."
                        label i (aref expected i) (aref output i)))
    (buffered-require (same-bits-p input input-before)
                      "~A: input changed." label)
    (buffered-require (equalp (csc:column-pointers matrix) pointers-before)
                      "~A: matrix pointers changed." label)
    (buffered-require (equalp (csc:row-indices matrix) rows-before)
                      "~A: matrix row indices changed." label)
    (buffered-require (same-bits-p (csc:entries matrix) values-before)
                      "~A: matrix values changed." label)
    (buffered-require
     (same-bits-p output (if transpose-p
                            (csc:transpose-matvec matrix input)
                            (csc:matvec matrix input)))
     "~A: allocating and buffered arithmetic differ in bits." label))
  output)

(defun buffered-rng (seed)
  "Return a local, repeatable integer generator, independent of CL random state."
  (let ((state (mod seed (ash 1 32))))
    (lambda (bound)
      (setf state (mod (+ (* 1664525 state) 1013904223) (ash 1 32)))
      (mod state bound))))

(defun buffered-explicit-cases ()
  (list
   (list "empty" 0 0 nil)
   (list "no-rows" 0 7 nil)
   (list "no-columns" 6 0 nil)
   (list "no-entries" 4 5 nil)
   (list "stored-zeros" 3 4 '((0 0 0d0) (2 2 -0d0) (1 3 0d0)))
   (list "wide-gaps" 4 7
         '((0 0 2d0) (3 0 -1d0) (1 1 0d0) (0 2 3d0)
           (1 2 -2d0) (3 4 0.125d0) (0 6 -4d0) (1 6 5d0)))
   (list "tall-gaps" 7 3
         '((0 0 -3d0) (2 0 1d0) (6 0 0d0)
           (1 2 0.25d0) (4 2 7d0) (6 2 -2d0)))
   (list "exact-cancellation" 3 5
         (list (list 0 0 (scale-float 1d0 40))
               (list 0 1 (- (scale-float 1d0 40)))
               '(0 2 0.375d0) '(1 1 3d0) '(1 3 -3d0)
               '(2 0 -0.5d0) '(2 4 0.5d0)))))

(defun buffered-seeded-case (index seed)
  "Small dyadic families keep every product and partial sum exactly representable."
  (let* ((next (buffered-rng (+ seed (* 7919 index))))
         (m (1+ (mod (* 7 index) 19)))
         (n (1+ (mod (+ (* 11 index) 5) 19)))
         (family (mod index 6))
         (exponent (nth (mod (floor index 6) 3) '(-200 0 200)))
         (entries nil))
    (dotimes (row m)
      (dotimes (column n)
        (when (case family
                ((0) t)
                ((1) (<= (abs (- row column)) 1))
                ((2) (zerop (mod (+ row column) 3)))
                ((3) (and (not (zerop (mod row 4)))
                          (not (zerop (mod column 5)))
                          (< (funcall next 5) 3)))
                ((4) (= row column))
                ((5) (< (funcall next 7) 2)))
          (push (list row column
                      (scale-float (/ (coerce (- (funcall next 33) 16) 'double-float)
                                      8d0)
                                   exponent))
                entries))))
    (list (format nil "seeded-~3,'0D" index) m n (nreverse entries))))

(defun buffered-run-exact-case (case seed)
  "Check two input vectors per orientation, reusing the same prefilled output."
  (destructuring-bind (name m n entries) case
    (let* ((matrix (primitive-from-entries
                    m n (reverse entries) (coalton:coalton csc:RejectDuplicates)))
           (dense (primitive-dense m n entries))
           (next (buffered-rng seed)))
      (dolist (transpose-p '(nil t))
        (let* ((input-size (if transpose-p m n))
               (output (make-array (if transpose-p n m)
                                   :element-type 'double-float :initial-element 913d0)))
          (dotimes (pass 2)
            (let ((input (make-array input-size :element-type 'double-float)))
              (dotimes (i input-size)
                (setf (aref input i)
                      (if (zerop pass)
                          (/ (coerce (- (funcall next 25) 12) 'double-float) 4d0)
                          1d0)))
              (buffered-check-call
               (format nil "~A/~A/pass~D" name (if transpose-p "At" "A") pass)
               matrix input output transpose-p
               (primitive-product dense input transpose-p)))))))))

(defun buffered-rounding-case (transpose-p)
  "The ordered terms 2^54, 1, -2^54 round to zero; their exact sum is one."
  (let* ((large (scale-float 1d0 54))
         (entries (if transpose-p
                      (list (list 0 0 large) '(1 0 1d0) (list 2 0 (- large)))
                      (list (list 0 0 large) '(0 1 1d0) (list 0 2 (- large))))))
    (primitive-from-entries (if transpose-p 3 1) (if transpose-p 1 3)
                            entries (coalton:coalton csc:RejectDuplicates))))

(defun buffered-rounding-tests ()
  (dolist (transpose-p '(nil t))
    (let ((matrix (buffered-rounding-case transpose-p)))
      (buffered-check-call
       (if transpose-p "rounded CSC order At" "rounded CSC order A")
       matrix (primitive-doubles 1d0 1d0 1d0) (primitive-doubles 913d0)
       transpose-p (primitive-doubles 0d0))))
  (format t "BUFFERED-ROUNDED-ORDER-PASS cases=2~%"))

(defun buffered-control-rejected (label thunk)
  "Accept only a mismatch from this oracle, never an unrelated API or runtime error."
  (let ((rejected nil))
    (handler-case (funcall thunk)
      (buffered-property-mismatch () (setf rejected t)))
    (assert rejected () "Buffered oracle missed planted defect ~A." label)))

(defun buffered-oracle-controls ()
  (let ((entries '((0 0 2d0) (1 0 -1d0) (0 1 3d0) (1 1 4d0))))
    (dolist
        (control
         (list
          (cons "missing contribution"
                (lambda (transpose-p matrix input output)
                  (buffered-public-product transpose-p matrix input output)
                  (decf (aref output 0) 2d0)))
          (cons "wrong orientation"
                (lambda (transpose-p matrix input output)
                  (buffered-public-product (not transpose-p) matrix input output)))
          (cons "retained stale output"
                (lambda (transpose-p matrix input output)
                  (buffered-public-product transpose-p matrix input output)
                  (incf (aref output 0) 913d0)))
          (cons "unwritten empty row"
                (lambda (transpose-p matrix input output)
                  (buffered-public-product transpose-p matrix input output)
                  (setf (aref output 2) 913d0)))
          (cons "input mutation"
                (lambda (transpose-p matrix input output)
                  (buffered-public-product transpose-p matrix input output)
                  (incf (aref input 0) 1d0)))))
      (let* ((matrix (primitive-from-entries
                      3 3 entries (coalton:coalton csc:RejectDuplicates)))
             (dense (primitive-dense 3 3 entries))
             (input (primitive-doubles 1d0 2d0 3d0))
             (output (primitive-doubles 913d0 913d0 913d0))
             (expected (primitive-product dense input nil)))
        (buffered-check-call "control baseline" matrix input output nil expected)
        (buffered-control-rejected
         (car control)
         (lambda ()
           (buffered-check-call (car control) matrix input output nil expected
                                (cdr control)))))))
  (let ((matrix (buffered-rounding-case nil))
        (input (primitive-doubles 1d0 1d0 1d0))
        (output (primitive-doubles 913d0))
        (expected (primitive-doubles 0d0)))
    (buffered-check-call "rounding control baseline" matrix input output nil expected)
    (buffered-control-rejected
     "reordered accumulation"
     (lambda ()
       (buffered-check-call
        "reordered accumulation" matrix input output nil expected
        (lambda (transpose-p a x out)
          (declare (ignore transpose-p))
          (let* ((values (csc:entries a))
                 (first (* (aref values 0) (aref x 0)))
                 (middle (* (aref values 1) (aref x 1)))
                 (last (* (aref values 2) (aref x 2))))
            (setf (aref out 0) (+ (+ first last) middle))))))))
  (format t "BUFFERED-ORACLE-CONTROLS-PASS rejected=6~%"))

(defun run-buffered-property-tests ()
  "Independent exact products plus separate ordinary-F64 accumulation controls."
  (let* ((seed 20261119)
         (cases (append (buffered-explicit-cases)
                        (loop :for index :below 96
                              :collect (buffered-seeded-case index seed)))))
    (loop :for case :in cases
          :for ordinal :from 0
          :do (buffered-run-exact-case case (+ seed (* 104729 ordinal))))
    (format t "BUFFERED-EXACT-PROPERTIES-PASS seed=~D cases=~D products=~D~%"
            seed (length cases) (* 4 (length cases))))
  (buffered-rounding-tests)
  (buffered-oracle-controls)
  t)
