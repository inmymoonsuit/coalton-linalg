(in-package #:coalton-linalg/tests)

;;; The analysis against an independent elimination written with explicit neighbour sets:
;;; the permutation, every column of the factor's structure and the entry budget, with and
;;; without held indices. The reference stores every fill edge and looks for the least key
;;; by a full scan; it shares no representation with the library's ordering.

(defun reference-elimination (n pointers rows leading)
  "Least degree, then least index; indices from LEADING on are held until an eliminated index below LEADING
releases one held neighbour, of least degree and then least index. Returns the order and each step's neighbours."
  (let ((adjacent (make-array n :initial-element '())) (alive (make-array n :initial-element t))
        (held (make-array n)) (order '()) (columns '()))
    (dotimes (j n)
      (loop for q from (1+ (aref pointers j)) below (aref pointers (1+ j)) do
        (let ((i (aref rows q)))
          (push j (aref adjacent i))
          (push i (aref adjacent j)))))
    (dotimes (i n) (setf (aref held i) (>= i leading)))
    (dotimes (step n)
      (let ((best nil) (best-key nil))
        (dotimes (i n)
          (when (aref alive i)
            (let ((key (+ (length (aref adjacent i)) (if (aref held i) (1+ n) 0))))
              (when (or (null best) (< key best-key)) (setf best i best-key key)))))
        (let ((around (sort (copy-list (aref adjacent best)) #'<)))
          (push best order)
          (push around columns)
          (setf (aref alive best) nil)
          (dolist (i around)
            (setf (aref adjacent i) (remove best (aref adjacent i)))
            (dolist (j around) (unless (= i j) (pushnew j (aref adjacent i)))))
          (when (< best leading)
            (let ((chosen nil))
              (dolist (i around)
                (when (and (aref held i)
                           (or (null chosen) (< (length (aref adjacent i)) (length (aref adjacent chosen)))))
                  (setf chosen i)))
              (when chosen (setf (aref held chosen) nil)))))))
    (values (nreverse order) (nreverse columns))))

(defun pattern-from-edges (n edges)
  "Lower CSC pattern arrays, diagonal first in every column, from a list of (i j) with i different from j."
  (let ((below (make-array n :initial-element '())))
    (dolist (edge edges)
      (pushnew (max (first edge) (second edge)) (aref below (min (first edge) (second edge)))))
    (let ((pointers (list 0)) (rows '()))
      (dotimes (j n)
        (push j rows)
        (dolist (i (sort (copy-list (aref below j)) #'<)) (push i rows))
        (push (length rows) pointers))
      (values (apply #'indices (reverse pointers)) (apply #'indices (reverse rows))))))

(defvar *ordering-cases* 0)
(defvar *ordering-entries* 0)

(defun check-ordering (label n edges)
  "Compare the library's analysis with the reference at several numbers of held indices."
  (multiple-value-bind (pointers rows) (pattern-from-edges n edges)
    (dolist (leading (remove-duplicates (list n 0 (floor n 3) (floor n 2) (max 0 (1- n)))))
      (multiple-value-bind (order columns) (reference-elimination n pointers rows leading)
        (let* ((pattern (ldlt:pattern n pointers rows))
               (factor (ldlt:symbolic-factor (if (= leading n)
                                                 (ldlt:analyze pattern)
                                                 (ldlt:analyze-quasidefinite pattern leading))))
               (perm (ldlt:factor-permutation factor)) (lp (ldlt:factor-columns factor))
               (lr (ldlt:factor-rows factor)) (position (make-array n))
               (entries (+ n (reduce #'+ columns :key #'length))))
          (incf *ordering-cases*)
          (incf *ordering-entries* entries)
          (assert (equal order (coerce perm 'list)) ()
                  "~A, ~D held from ~D: order ~S, reference ~S." label (- n leading) leading perm order)
          (dotimes (k n) (setf (aref position (aref perm k)) k))
          (assert (= n (1- (length lp))))
          (loop for k below n for column in columns do
            (let ((expected (sort (mapcar (lambda (i) (aref position i)) column) #'<))
                  (actual (coerce (subseq lr (aref lp k) (aref lp (1+ k))) 'list)))
              (assert (equal expected actual) ()
                      "~A, held from ~D: column ~D holds rows ~S, reference ~S." label leading k actual expected)))
          ;; The budget is the factor's entries with its diagonal: one fewer is refused and names them.
          (when (= leading n)
            (flet ((limit-error (limit)
                     (uiop:symbol-call :coalton-linalg/tests/typed-status :analysis-pattern-limit-error
                                       n pointers rows limit entries)))
              (assert (= 0 (limit-error entries)) () "~A: a budget of exactly ~D entries was refused." label entries)
              (when (plusp entries)
                (assert (= 8 (limit-error (1- entries))) ()
                        "~A: a budget of ~D for ~D entries gave code ~D." label (1- entries) entries
                        (limit-error (1- entries)))))))))))

(defun seeded-edges (state n percent)
  (loop for i below n append
    (loop for j below i when (< (random 100 state) percent) collect (list i j))))

(defun run-ordering-tests ()
  (setf *ordering-cases* 0 *ordering-entries* 0)
  (let ((state (sb-ext:seed-random-state 20261128)))
    ;; Shapes with a known character: nothing to do, no fill, all fill, one hub, long fronts.
    (check-ordering "one vertex" 1 '())
    (check-ordering "no edges" 5 '())
    (check-ordering "one edge" 2 '((1 0)))
    (check-ordering "path" 9 (loop for i from 1 below 9 collect (list i (1- i))))
    (check-ordering "cycle" 8 (cons '(7 0) (loop for i from 1 below 8 collect (list i (1- i)))))
    (check-ordering "star, hub first" 9 (loop for i from 1 below 9 collect (list i 0)))
    (check-ordering "star, hub last" 9 (loop for i below 8 collect (list 8 i)))
    (check-ordering "complete" 7 (loop for i below 7 append (loop for j below i collect (list i j))))
    (check-ordering "complete less one edge" 7
                    (remove '(6 0) (loop for i below 7 append (loop for j below i collect (list i j))) :test #'equal))
    (check-ordering "two components" 10
                    (append (loop for i from 1 below 5 collect (list i 0))
                            (loop for i from 6 below 10 append (loop for j from 5 below i collect (list i j)))))
    (dolist (side '(3 4 5 6))
      (check-ordering (format nil "grid ~D" side) (* side side)
                      (loop for r below side append
                        (loop for c below side
                              append (append (and (< (1+ c) side) (list (list (+ (* r side) c 1) (+ (* r side) c))))
                                             (and (< (1+ r) side) (list (list (+ (* (1+ r) side) c) (+ (* r side) c)))))))))
    ;; Two cliques joined through one vertex, and a clique with a tail: elements that absorb others.
    (check-ordering "bow tie" 9
                    (append (loop for i below 5 append (loop for j below i collect (list i j)))
                            (loop for i from 4 below 9 append (loop for j from 4 below i collect (list i j)))))
    (check-ordering "clique with a tail" 10
                    (append (loop for i below 6 append (loop for j below i collect (list i j)))
                            (loop for i from 6 below 10 collect (list i (1- i)))))
    ;; A wheel: a hub adjacent to every rim vertex, which leaves the hub in many small elements.
    (check-ordering "wheel" 41 (append (loop for i from 1 to 40 collect (list i 0))
                                       (loop for i from 1 to 40 collect (list (1+ (mod i 40)) i))))
    ;; Two hubs over a path: both are in every element a path vertex leaves.
    (check-ordering "two hubs" 30 (append (loop for i from 2 below 30 append (list (list i 0) (list i 1)))
                                          (loop for i from 3 below 30 collect (list i (1- i)))))
    ;; An augmented pattern: columns 0..5, rows 6..9, a hub column and two rows sharing a column.
    (check-ordering "augmented" 10 '((6 0) (7 0) (8 0) (9 0) (6 1) (7 2) (8 3) (9 3) (9 4)))
    ;; Seeded graphs from sparse to nearly complete.
    (dolist (n '(2 3 4 5 6 8 10 13 17 24 32 45))
      (dolist (percent '(4 10 20 35 60 90))
        (dotimes (repeat (if (<= n 10) 4 2))
          (check-ordering (format nil "seeded n=~D ~D% #~D" n percent repeat) n (seeded-edges state n percent)))))
    ;; Larger sparse graphs, where most vertices are light and elements overlap without nesting.
    (dolist (n '(80 130))
      (dolist (percent '(2 4 7))
        (check-ordering (format nil "seeded n=~D ~D%" n percent) n (seeded-edges state n percent)))))
  (format t "ORDERING-REFERENCE-PASS analyses=~D factor-entries=~D~%" *ordering-cases* *ordering-entries*))
