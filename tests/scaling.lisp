(in-package #:coalton-linalg/tests)

;;; Scaling by powers of two: every claim of docs/scaling.md checked exactly, on the stored values.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "SC" "COALTON-LINALG/SPARSE/SCALING")
  (sb-ext:add-package-local-nickname "SC-CALLER" "COALTON-LINALG/TESTS/SCALING-CALLERS"))

(defun power-of-two-p (x)
  (and (plusp x) (= 0.5d0 (decode-float x))))

(defun zeros (n) (make-array n :element-type 'double-float :initial-element 0d0))

(defun stored-entries (matrix values)
  "Every stored entry as (row column value) in the matrix's CSC order, with VALUES in place of its own."
  (let ((cp (csc:column-pointers matrix)) (ri (csc:row-indices matrix)) (found '()))
    (dotimes (j (csc:column-count matrix))
      (loop for q from (aref cp j) below (aref cp (1+ j)) do
        (push (list (aref ri q) j (aref values q)) found)))
    (nreverse found)))

(defun largest-by (entries count key)
  "The largest magnitude among the entries of each row (KEY #'first) or column (KEY #'second), exactly."
  (let ((largest (make-array count :initial-element 0)))
    (dolist (entry entries)
      (setf (aref largest (funcall key entry))
            (max (aref largest (funcall key entry)) (abs (rational (third entry))))))
    largest))

(defvar *scaling-cases* 0)

(defun check-equilibrated (label m n entries)
  "Equilibrate, scale and inspect one rectangular matrix through the public interface, checking each claim exactly."
  (let* ((matrix (csc-from-entries m n entries)) (values (csc:entries matrix))
         (rows (zeros m)) (columns (zeros n)) (scaled (zeros (length values)))
         (before (copy-seq values)))
    (incf *scaling-cases*)
    (assert (= 0 (sc-caller:equilibrate-error matrix values rows columns)) () "~A: equilibrate! refused." label)
    (assert (same-bits-p values before) () "~A: equilibrate! changed the values." label)
    (assert (and (every #'power-of-two-p rows) (every #'power-of-two-p columns)) ()
            "~A: a scale factor is not a power of two: ~S ~S." label rows columns)
    (assert (= 0 (sc-caller:scale-error matrix values rows columns scaled)) () "~A: scale! refused." label)
    (let* ((original (stored-entries matrix values)) (result (stored-entries matrix scaled))
           (row-before (largest-by original m #'first)) (column-before (largest-by original n #'second))
           (row-after (largest-by result m #'first)) (column-after (largest-by result n #'second)))
      ;; Each scaled value is the exact product: powers of two change no digit.
      (loop for (i j v) in original for (nil nil s) in result do
        (assert (= (rational s) (* (rational (aref rows i)) (rational v) (rational (aref columns j)))) ()
                "~A: entry ~D,~D was scaled to ~S, not the exact product." label i j s))
      ;; A row or column with no nonzero keeps a factor of 1; every other has its largest
      ;; scaled magnitude in [1/2, 1).
      (dotimes (i m)
        (if (zerop (aref row-before i))
            (assert (= 1d0 (aref rows i)) () "~A: empty row ~D was scaled by ~S." label i (aref rows i))
            (assert (and (<= 1/2 (aref row-after i)) (< (aref row-after i) 1)) ()
                    "~A: row ~D has largest scaled magnitude ~E." label i (float (aref row-after i) 1d0))))
      (dotimes (j n)
        (if (zerop (aref column-before j))
            (assert (= 1d0 (aref columns j)) () "~A: empty column ~D was scaled by ~S." label j (aref columns j))
            (assert (and (<= 1/2 (aref column-after j)) (< (aref column-after j) 1)) ()
                    "~A: column ~D has largest scaled magnitude ~E." label j (float (aref column-after j) 1d0))))
      ;; The maxima the library reports are the exact ones.
      (let ((reported-rows (zeros m)) (reported-columns (zeros n)))
        (assert (= 0 (sc-caller:maxima-error matrix scaled reported-rows t)))
        (assert (= 0 (sc-caller:maxima-error matrix scaled reported-columns nil)))
        (assert (and (every (lambda (a b) (= (rational a) b)) reported-rows row-after)
                     (every (lambda (a b) (= (rational a) b)) reported-columns column-after))
                () "~A: reported maxima ~S ~S." label reported-rows reported-columns))
      ;; An equilibrated matrix needs no further scaling.
      (let ((again-rows (zeros m)) (again-columns (zeros n)))
        (assert (= 0 (sc-caller:equilibrate-error matrix scaled again-rows again-columns)))
        (assert (and (every (lambda (r) (= r 1d0)) again-rows) (every (lambda (c) (= c 1d0)) again-columns)) ()
                "~A: a second equilibration gave ~S ~S." label again-rows again-columns)))
    (values rows columns scaled)))

(defun ranged-entries (state m n percent spread)
  "Stored entries whose magnitudes range over 2^-SPREAD to 2^SPREAD, with a few stored zeros."
  (loop for i below m append
    (loop for j below n
          when (< (random 100 state) percent)
            collect (list i j (if (zerop (random 12 state))
                                  0d0
                                  (* (if (zerop (random 2 state)) 1 -1)
                                     (scale-float (+ 1d0 (random 1d0 state)) (- (random (1+ (* 2 spread)) state) spread))))))))

(defun lower-csc (n dense)
  "The lower triangle of a dense symmetric matrix of doubles as a checked CSC matrix, every diagonal stored."
  (csc-from-entries n n (loop for j below n append (loop for i from j below n collect (list i j (aref dense i j))))))

(defun check-symmetric-scaled (label n dense)
  "Scale a symmetric positive definite matrix by its diagonal and check each claim exactly."
  (let* ((matrix (lower-csc n dense)) (values (csc:entries matrix))
         (scales (zeros n)) (scaled (zeros (length values))))
    (incf *scaling-cases*)
    (assert (= 0 (sc-caller:symmetric-code matrix values scales)) () "~A: the diagonal was refused." label)
    (assert (every #'power-of-two-p scales) () "~A: scale factors ~S." label scales)
    (assert (= 0 (sc-caller:symmetric-scale-error matrix values scales scaled)))
    (loop for (i j v) in (stored-entries matrix values) for (nil nil s) in (stored-entries matrix scaled) do
      (assert (= (rational s) (* (rational (aref scales i)) (rational v) (rational (aref scales j)))) ()
              "~A: entry ~D,~D was scaled to ~S, not the exact product." label i j s)
      (if (= i j)
          (assert (and (<= 1/2 (rational s)) (< (rational s) 2)) () "~A: scaled diagonal ~D is ~S." label i s)
          ;; For a positive semidefinite matrix an entry is at most the geometric mean of its two diagonals.
          (assert (< (abs (rational s)) 2) () "~A: scaled entry ~D,~D is ~S." label i j s)))
    (values scales scaled)))

(defun run-scaling-exact-tests ()
  (let ((state (sb-ext:seed-random-state 20261201)))
    (check-equilibrated "one entry" 1 1 '((0 0 3)))
    (check-equilibrated "identity" 3 3 '((0 0 1) (1 1 1) (2 2 1)))
    (check-equilibrated "already in range" 2 2 '((0 0 0.75d0) (0 1 0.5d0) (1 1 0.9d0)))
    ;; Lengths beside angles: one row a million times the other, one column a thousandth.
    (check-equilibrated "mixed units" 2 3 '((0 0 2d6) (0 1 3d6) (0 2 4d3) (1 0 1.5d0) (1 2 2d-3)))
    ;; The column pass has to lift a column whose entries are small in every row.
    (multiple-value-bind (rows columns)
        (check-equilibrated "a small column" 2 2 (list '(0 0 1) (list 0 1 (scale-float 1d0 -40))
                                                       '(1 0 1) (list 1 1 (scale-float 1d0 -37))))
      (assert (equalp rows (doubles 0.5d0 0.5d0)))
      (assert (equalp columns (doubles 1d0 (scale-float 1d0 37)))))
    (check-equilibrated "stored zeros" 3 3 '((0 0 0) (1 0 4) (1 1 0) (2 2 0) (0 2 5)))
    (check-equilibrated "empty row and column" 3 3 '((0 0 2) (2 0 8)))
    (check-equilibrated "no entries" 2 3 '())
    (check-equilibrated "wide" 1 4 (list '(0 0 1d200) '(0 1 -3) (list 0 2 (scale-float 1d0 -200)) '(0 3 7d100)))
    (check-equilibrated "tall" 4 1 (list '(0 0 1d200) '(1 0 -3) (list 2 0 (scale-float 1d0 -200)) '(3 0 7d100)))
    (check-equilibrated "extremes" 2 2 (list (list 0 0 most-positive-double-float) '(0 1 1)
                                             (list 1 1 least-positive-normalized-double-float)))
    (loop for (m n) in '((2 2) (3 3) (3 5) (5 3) (4 4) (6 6) (5 8) (8 5) (8 8)) do
      (dolist (percent '(30 60 100))
        (dolist (spread '(3 40 300))
          (dotimes (repeat 2)
            (check-equilibrated (format nil "seeded ~Dx~D ~D% spread ~D #~D" m n percent spread repeat)
                                m n (ranged-entries state m n percent spread))))))
    ;; Symmetric positive definite: Gram matrices of columns that differ by many orders.
    (loop for n in '(1 2 3 4 6 8) do
      (dolist (spread '(0 10 100))
        (dotimes (repeat 4)
          (let* ((b (make-array (list (+ n 2) n)))
                 (dense (make-array (list n n) :element-type 'double-float)))
            (dotimes (j n)
              (let ((column-scale (scale-float 1d0 (- (random (1+ (* 2 spread)) state) spread))))
                (dotimes (i (+ n 2))
                  (setf (aref b i j) (* column-scale (- (random 9 state) 4))))
                ;; A column of zeros would give a zero diagonal; keep every column present.
                (setf (aref b j j) (* column-scale (1+ (random 4 state))))))
            (dotimes (i n)
              (dotimes (j n)
                (setf (aref dense i j) (loop for k below (+ n 2) sum (* (aref b k i) (aref b k j))))))
            (check-symmetric-scaled (format nil "gram n=~D spread ~D #~D" n spread repeat) n dense)))))
    (multiple-value-bind (scales scaled)
        (check-symmetric-scaled "named" 2 (make-array '(2 2) :element-type 'double-float
                                                             :initial-contents '((4d0 2d0) (2d0 64d0))))
      ;; 4 = 0.5 2^3 and 64 = 0.5 2^7: factors 2^-1 and 2^-3, scaled diagonal 1 and 1, off-diagonal 2 / 16.
      (assert (equalp scales (doubles 0.5d0 0.125d0)))
      (assert (equalp scaled (doubles 1d0 0.125d0 1d0)))))
  (format t "SCALING-EXACT-PASS cases=~D~%" *scaling-cases*))

(defun run-scaling-edge-tests ()
  ;; A row whose largest magnitude is subnormal cannot be brought to 1/2: its factor is the largest power of two.
  (let* ((tiny (* 16 least-positive-double-float))
         (matrix (csc-from-entries 2 2 (list (list 0 0 tiny) '(1 1 1))))
         (rows (zeros 2)) (columns (zeros 2)) (scaled (zeros 2)))
    (assert (= 0 (sc-caller:equilibrate-error matrix (csc:entries matrix) rows columns)))
    (assert (= (aref rows 0) (scale-float 1d0 1023)) () "A subnormal row was scaled by ~S." (aref rows 0))
    (assert (= 0 (sc-caller:scale-error matrix (csc:entries matrix) rows columns scaled)))
    ;; The column pass finishes what the row pass could not.
    (assert (and (<= 0.5d0 (aref scaled 0)) (< (aref scaled 0) 1d0)) () "The subnormal entry became ~S." (aref scaled 0)))
  ;; A subnormal diagonal is likewise beyond reach: its factor is the largest the method has, 2^511.
  (let* ((matrix (csc-from-entries 2 2 (list (list 0 0 (scale-float 1d0 -1040)) '(1 1 1))))
         (scales (zeros 2)) (scaled (zeros 2)))
    (assert (= 0 (sc-caller:symmetric-code matrix (csc:entries matrix) scales)))
    (assert (equalp scales (doubles (scale-float 1d0 511) 1d0)) () "Scales ~S." scales)
    (assert (= 0 (sc-caller:symmetric-scale-error matrix (csc:entries matrix) scales scaled)))
    (assert (equalp scaled (doubles (scale-float 1d0 -18) 1d0)) () "Scaled ~S." scaled))
  ;; An entry far below its row's largest underflows when scaled: it is rounded, not refused, and the
  ;; pattern keeps it.
  (let* ((matrix (csc-from-entries 1 2 (list (list 0 0 (scale-float 1d0 1000)) (list 0 1 (scale-float 1d0 -1000)))))
         (rows (zeros 1)) (columns (zeros 2)) (scaled (zeros 2)))
    (assert (= 0 (sc-caller:equilibrate-error matrix (csc:entries matrix) rows columns)))
    (assert (= 0 (sc-caller:scale-error matrix (csc:entries matrix) rows columns scaled)))
    (assert (= 0.5d0 (aref scaled 0)))
    (assert (zerop (aref scaled 1)) () "An entry 2^-2001 of its row's size became ~S." (aref scaled 1)))
  ;; Refusals name their kind and leave the outputs alone.
  (let* ((matrix (csc-from-entries 2 3 '((0 0 1) (0 2 -2) (1 1 3)))) (values (csc:entries matrix))
         (rows (doubles 7d0 7d0)) (columns (doubles 7d0 7d0 7d0)) (out (doubles 7d0 7d0 7d0))
         (bad (copy-seq values)))
    (flet ((untouched ()
             (and (equalp rows (doubles 7d0 7d0)) (equalp columns (doubles 7d0 7d0 7d0))
                  (equalp out (doubles 7d0 7d0 7d0)))))
      (assert (= 2 (sc-caller:equilibrate-error matrix (doubles 1d0 2d0) rows columns)))
      (assert (= 2 (sc-caller:equilibrate-error matrix values (doubles 1d0) columns)))
      (assert (= 2 (sc-caller:equilibrate-error matrix values rows (doubles 1d0 1d0))))
      ;; Aliasing can only be met where the lengths agree: a square matrix with as many entries as rows.
      (let* ((square (csc-from-entries 2 2 '((0 0 1) (1 1 2)))) (stored (csc:entries square))
             (both (doubles 7d0 7d0)) (other (doubles 7d0 7d0)))
        (assert (= 9 (sc-caller:equilibrate-error square stored stored other)))
        (assert (= 9 (sc-caller:equilibrate-error square stored other stored)))
        (assert (= 9 (sc-caller:equilibrate-error square stored both both)))
        (assert (= 9 (sc-caller:maxima-error square stored stored t)))
        (assert (= 9 (sc-caller:maxima-error square stored stored nil)))
        (assert (and (equalp both (doubles 7d0 7d0)) (equalp other (doubles 7d0 7d0))
                     (equalp stored (doubles 1d0 2d0)))))
      (setf (aref bad 1) sb-ext:double-float-positive-infinity)
      (assert (= 5 (sc-caller:equilibrate-error matrix bad rows columns)))
      (assert (= 5 (sc-caller:scale-error matrix bad (doubles 1d0 1d0) (doubles 1d0 1d0 1d0) out)))
      (assert (= 5 (sc-caller:maxima-error matrix bad rows t)))
      (assert (= 5 (sc-caller:maxima-error matrix bad columns nil)))
      (assert (= 2 (sc-caller:maxima-error matrix values columns t)))
      (assert (= 2 (sc-caller:maxima-error matrix values rows nil)))
      (assert (= 2 (sc-caller:scale-error matrix values (doubles 1d0) (doubles 1d0 1d0 1d0) out)))
      (assert (= 2 (sc-caller:scale-error matrix values (doubles 1d0 1d0) (doubles 1d0 1d0) out)))
      (assert (= 2 (sc-caller:scale-error matrix values (doubles 1d0 1d0) (doubles 1d0 1d0 1d0) (doubles 0d0))))
      (assert (= 5 (sc-caller:scale-error matrix values (doubles 1d0 (not-a-number)) (doubles 1d0 1d0 1d0) out)))
      (assert (= 5 (sc-caller:scale-error matrix values (doubles 1d0 1d0)
                                          (doubles 1d0 sb-ext:double-float-negative-infinity 1d0) out)))
      (assert (= 9 (sc-caller:scale-error matrix values (doubles 1d0 1d0) (doubles 1d0 1d0 1d0) values)))
      (let ((three (doubles 1d0 1d0 1d0)))
        (assert (= 9 (sc-caller:scale-error matrix values (doubles 1d0 1d0) three three))))
      (assert (untouched) () "A refused call wrote ~S ~S ~S." rows columns out)
      ;; Factors the caller chose can overflow; that is an arithmetic refusal.
      (assert (= 6 (sc-caller:scale-error matrix values (doubles 1d300 1d0) (doubles 1d300 1d0 1d0) out)))))
  ;; Vectors: the output may be the input, never the scales.
  (let ((scales (doubles 2d0 0.5d0 4d0)) (x (doubles 3d0 -8d0 0.25d0)) (out (doubles 7d0 7d0 7d0)))
    (assert (= 0 (sc-caller:vector-error scales x out)))
    (assert (equalp out (doubles 6d0 -4d0 1d0)))
    (assert (= 0 (sc-caller:vector-error scales x x)))
    (assert (equalp x (doubles 6d0 -4d0 1d0)))
    (assert (= 9 (sc-caller:vector-error scales x scales)))
    (assert (= 2 (sc-caller:vector-error scales (doubles 1d0) out)))
    (assert (= 2 (sc-caller:vector-error scales x (doubles 1d0))))
    (assert (= 5 (sc-caller:vector-error scales (doubles 1d0 (not-a-number) 1d0) out)))
    (assert (= 5 (sc-caller:vector-error (doubles 1d0 sb-ext:double-float-positive-infinity 1d0) x out)))
    (assert (= 6 (sc-caller:vector-error (doubles 1d300 1d0 1d0) (doubles 1d300 1d0 1d0) out)))
    (assert (equalp scales (doubles 2d0 0.5d0 4d0))))
  ;; Symmetric: the structure is checked, a diagonal that is not positive is an outcome, and
  ;; neither writes the scales.
  (let ((scales (doubles 7d0 7d0)) (out (doubles 7d0 7d0 7d0)))
    (flet ((code (m n entries) (let ((matrix (csc-from-entries m n entries)))
                                 (sc-caller:symmetric-code matrix (csc:entries matrix) scales))))
      (assert (= 10 (code 2 2 '((0 0 4) (1 0 1) (1 1 0)))))
      (assert (= 10 (code 2 2 '((0 0 -4) (1 0 1) (1 1 9)))))
      (assert (= 1 (code 2 2 '((0 0 4) (1 0 1)))))
      (assert (= 1 (code 2 2 '((0 0 4) (0 1 1) (1 1 9)))))
      (assert (= 1 (code 2 2 '((1 0 1) (1 1 9))))))
    (let ((wide (csc-from-entries 2 3 '((0 0 1) (1 1 1)))))
      (assert (= 2 (sc-caller:symmetric-code wide (csc:entries wide) scales))))
    (let* ((matrix (csc-from-entries 2 2 '((0 0 4) (1 0 1) (1 1 9)))) (values (csc:entries matrix)))
      (assert (= 2 (sc-caller:symmetric-code matrix values (doubles 1d0))))
      (let* ((diagonal (csc-from-entries 2 2 '((0 0 4) (1 1 9)))) (stored (csc:entries diagonal)))
        (assert (= 9 (sc-caller:symmetric-code diagonal stored stored)))
        (assert (equalp stored (doubles 4d0 9d0))))
      (assert (= 5 (sc-caller:symmetric-code matrix (doubles 4d0 (not-a-number) 9d0) scales)))
      (assert (= 2 (sc-caller:symmetric-scale-error matrix values (doubles 1d0) out)))
      (assert (= 2 (sc-caller:symmetric-scale-error matrix values (doubles 1d0 1d0) (doubles 0d0))))
      (assert (= 9 (sc-caller:symmetric-scale-error matrix values (doubles 1d0 1d0) values)))
      (assert (= 5 (sc-caller:symmetric-scale-error matrix values (doubles 1d0 (not-a-number)) out)))
      (assert (= 6 (sc-caller:symmetric-scale-error matrix values (doubles 1d200 1d200) out))))
    (assert (equalp scales (doubles 7d0 7d0)) () "A refused symmetric call wrote the scales: ~S." scales))
  (format t "SCALING-EDGE-PASS~%"))

(defun run-scaling-tests ()
  (setf *scaling-cases* 0)
  (run-scaling-exact-tests)
  (run-scaling-edge-tests)
  ;; Warmed repeated work through a compiled Coalton caller allocates nothing.
  (let* ((state (sb-ext:seed-random-state 7))
         (matrix (csc-from-entries 6 9 (ranged-entries state 6 9 70 30)))
         (lower (lower-csc 3 (make-array '(3 3) :element-type 'double-float
                                                :initial-contents '((4d0 1d0 0d0) (1d0 900d0 2d0) (0d0 2d0 1d-6)))))
         (rows (zeros 6)) (columns (zeros 9)) (scaled (zeros (csc:entry-count matrix)))
         (x (make-array 6 :element-type 'double-float :initial-element 1d0)) (y (zeros 6))
         (d (zeros 3)) (lower-scaled (zeros (csc:entry-count lower))))
    (assert (= 0 (sc-caller:scaling-loop matrix rows columns scaled x y lower d lower-scaled 3)) ()
            "The allocation case does not run.")
    (allocation-check "scaling-loop"
                      (lambda () (sc-caller:scaling-loop matrix rows columns scaled x y lower d lower-scaled 1))))
  (format t "SCALING-TESTS-PASS~%"))
