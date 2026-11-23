(in-package #:coalton-linalg/tests)

;;; The aliasing roster. One question is asked of every exported call that takes two or more
;;; arrays: what happens when two of them are one array? The answer the library gives is that
;;; such a call either does exactly what it does with two distinct arrays of equal contents,
;;; or refuses before it writes anything.
;;;
;;; The members are not listed by hand. They are read from the sources: every exported
;;; function whose declared type has two or more arrays among its own parameters. A call
;;; added to the library joins the inventory by being declared, and this file fails until it
;;; has an entry here or a stated reason for having none.

;;; ---- The inventory, from the sources ---------------------------------------------------

(defun alias-whitespace-p (character)
  (member character '(#\Space #\Newline #\Return #\Tab)))

(defun alias-balanced (text start)
  "The parenthesized form that opens at START."
  (let ((depth 0))
    (loop :for i :from start :below (length text) :do
      (case (char text i)
        ((#\() (incf depth))
        ((#\)) (decf depth) (when (zerop depth) (return-from alias-balanced (subseq text start (1+ i)))))))
    (error "Unbalanced form at ~D." start)))

(defun alias-parameters (type-text)
  "The parameters of a function type '(A * B * (C * D -> E) -> F)', split where it is one level deep. NIL for a
type that is not a function's."
  (let ((inner (subseq type-text 1 (1- (length type-text)))) (depth 0) (parts '()) (current '()))
    (flet ((finish () (push (string-trim " " (coerce (reverse current) 'string)) parts) (setf current '())))
      (loop :for i :from 0 :below (length inner) :do
        (let ((c (char inner i)))
          (case c ((#\() (incf depth)) ((#\)) (decf depth)))
          (cond
            ((and (zerop depth) (char= c #\-) (< (1+ i) (length inner)) (char= (char inner (1+ i)) #\>))
             (finish)
             (return-from alias-parameters (nreverse parts)))
            ((and (zerop depth) (char= c #\*) (plusp i) (alias-whitespace-p (char inner (1- i)))
                  (< (1+ i) (length inner)) (alias-whitespace-p (char inner (1+ i))))
             (finish))
            (t (push c current))))))
    nil))

(defun alias-exported-names (text)
  (let ((start (search "(:export" text)) (end (search "(cl:in-package" text)) (names '()))
    (when (and start end)
      (loop :with i := start
            :for at := (search "#:" text :start2 i :end2 end)
            :while at
            :do (let ((stop (position-if (lambda (c) (or (alias-whitespace-p c) (member c '(#\( #\))))) text
                                         :start (+ at 2))))
                  (push (subseq text (+ at 2) stop) names)
                  (setf i stop))))
    names))

(defun source-inventory (member-p)
  "Every exported name of the library whose declared type MEMBER-P accepts, as 'package:name' strings, read from the
sources. MEMBER-P is given the type as written, on one line."
  (let ((found '())
        (root (asdf:system-source-directory "coalton-linalg")))
    (dolist (path (directory (merge-pathnames "src/**/*.ct" root)))
      (let* ((text (substitute #\Space #\Return (uiop:read-file-string path :external-format :utf-8)))
             (package-start (+ (search "(cl:defpackage #:" text) (length "(cl:defpackage #:")))
             (package (subseq text package-start (position-if #'alias-whitespace-p text :start package-start)))
             (exported (alias-exported-names text)))
        (loop :with i := 0
              :for at := (search "(declare " text :start2 i)
              :while at
              :do (let* ((name-start (+ at (length "(declare ")))
                         (name-end (position-if #'alias-whitespace-p text :start name-start))
                         (name (subseq text name-start name-end))
                         (form (alias-balanced text at))
                         (type-text (string-trim " " (substitute #\Space #\Newline
                                                                 (subseq form (- name-end at) (1- (length form)))))))
                    (setf i (+ at (length form)))
                    (when (and (member name exported :test #'string=) (plusp (length type-text))
                               (char= (char type-text 0) #\() (funcall member-p type-text))
                      (push (format nil "~A:~A" package name) found))))))
    (sort found #'string<)))

(defun alias-inventory ()
  "The exported functions with two or more arrays among their own parameters."
  (source-inventory (lambda (type-text)
                      (<= 2 (count-if (lambda (parameter) (eql 0 (search "arr:LispArray" parameter)))
                                      (alias-parameters type-text))))))

;;; ---- The verdict on one call -------------------------------------------------------------

(defun alias-same-p (a b)
  (and (= (length a) (length b))
       (if (typep a '(simple-array double-float (*))) (same-bits-p a b) (equalp a b))))

(defun alias-run (call arrays)
  "The call's outcome as a code, 0 for success, and whatever further arrays it returned. The library's typed
error is a refusal. Any other condition is a defect and is not caught."
  (handler-case (multiple-value-bind (code extras) (funcall call arrays)
                  (values (if (integerp code) code 0) extras))
    (coalton-linalg/errors:linalgerror () (values 99 nil))))

(defun alias-pairs (arrays)
  "Every pair of positions holding arrays of one element type and one length: the pairs that can be one array."
  (loop :for i :from 0 :below (length arrays) :append
    (loop :for j :from (1+ i) :below (length arrays)
          :when (and (equal (array-element-type (nth i arrays)) (array-element-type (nth j arrays)))
                     (= (length (nth i arrays)) (length (nth j arrays))))
            :collect (list i j))))

(defun alias-verdict (label build call i j)
  "Call once with positions I and J holding two arrays of equal contents, and once with them holding one array.
:REFUSED when the second call refuses and leaves every array as it was; :SAME when it does what the first did.
Anything else is an error."
  (let* ((distinct (funcall build)) (shared (funcall build)))
    (setf (nth j distinct) (copy-seq (nth i distinct)))
    (setf (nth j shared) (nth i shared))
    (let ((before (mapcar #'copy-seq distinct)) (shared-before (mapcar #'copy-seq shared)))
      (multiple-value-bind (reference-code reference-extras) (alias-run call distinct)
        ;; The premise: with two arrays the call is one the library accepts. Without it a mistake in
        ;; an entry would read as a refusal.
        (assert (= 0 reference-code) ()
                "~A: the call with arguments ~D and ~D two equal arrays gives ~D; the entry is wrong." label i j
                reference-code)
        (multiple-value-bind (code extras) (alias-run call shared)
          (cond
            ((/= code 0)
             (loop :for array :in shared :for was :in shared-before :for position :from 0 :do
               (assert (alias-same-p array was) ()
                       "~A: with arguments ~D and ~D one array the call refused, code ~D, after writing argument ~D."
                       label i j code position))
             :refused)
            (t
             (let ((wrote-i (not (alias-same-p (nth i distinct) (nth i before))))
                   (wrote-j (not (alias-same-p (nth j distinct) (nth j before)))))
               (assert (not (and wrote-i wrote-j)) ()
                       "~A: arguments ~D and ~D are both written, and the call accepted them as one array." label i j)
               (loop :for position :from 0 :for array :in shared :for reference :in distinct :do
                 (let ((expected (cond ((and (= position i) wrote-j) (nth j distinct))
                                       ((and (= position j) wrote-i) (nth i distinct))
                                       (t reference))))
                   (assert (alias-same-p array expected) ()
                           "~A: with arguments ~D and ~D one array, argument ~D ends as ~S; with two equal arrays, ~S."
                           label i j position array expected)))
               (assert (and (= (length extras) (length reference-extras))
                            (every #'alias-same-p extras reference-extras))
                       () "~A: with arguments ~D and ~D one array the result differs." label i j))
             :same)))))))

;;; ---- The members -------------------------------------------------------------------------

;;; The matrices are chosen so that every factorization here reorders or pivots. A solve that read its
;;; right-hand side after it had begun to write its result would then read the wrong entries; with a
;;; matrix taken in its own order it could read the right ones by accident.

(defun alias-arrow ()
  "Symmetric positive definite, the first row and column full: the ordering moves that column to the end."
  (rows-matrix '((4 1 1 1) (1 3 0 0) (1 0 3 0) (1 0 0 3))))

(defun alias-unsymmetric ()
  "Nonsingular with a zero on the diagonal in two columns, so that rows are exchanged."
  (rows-matrix '((0 2 0 1) (3 0 1 0) (0 1 4 0) (1 0 2 5))))

(defun alias-indefinite ()
  "Symmetric and indefinite with two zeros on the diagonal, so that the pivots are not taken in order."
  (rows-matrix '((0 1 2) (1 0 3) (2 3 1))))

(defun alias-columns ()
  "Nonsingular, the largest column not the first, so that columns are exchanged."
  (rows-matrix '((1 5 2) (0 1 7) (3 0 1))))

(defun alias-stored (dense)
  "The nonzero entries of a dense square matrix as a checked CSC matrix."
  (let ((n (array-dimension dense 0)))
    (csc-from-entries n n (loop :for j :below n :append
                            (loop :for i :below n
                                  :unless (zerop (aref dense i j)) :collect (list i j (aref dense i j)))))))

(defun alias-diagonal ()
  "Three rows, three columns, three entries: every value, row and column array has one length."
  (csc-from-entries 3 3 '((0 0 2) (1 1 4) (2 2 8))))

(defun alias-reordered (label permutation)
  "The premise of a fixture: its factorization is not taken in the matrix's own order."
  (assert (loop :for k :from 0 :for entry :across permutation :thereis (/= k entry)) ()
          "~A is factored in its own order, ~S: the fixture no longer reorders." label permutation))

(defun alias-arrow-symbolic ()
  "The analysis of the arrow matrix and its values. The premise, that its order is not the matrix's own, is
asserted here."
  (multiple-value-bind (order pointers rows values) (symmetric-lower-pattern (alias-arrow))
    (let ((symbolic (ldlt:analyze (ldlt:pattern order pointers rows))))
      (alias-reordered "The arrow matrix" (ldlt:factor-permutation (ldlt:symbolic-factor symbolic)))
      (values symbolic values))))

(defun alias-success (status expected)
  (if (eq status expected) 0 1))

(defun alias-roster ()
  "One entry for each member: its key, how to make its arrays afresh, how to call it on them, and the outcome
expected of each pair that can be one array. An entry of (key :exempt reason) has no test."
  (let ((definite (alias-stored (alias-arrow))) (unsymmetric (alias-stored (alias-unsymmetric)))
        (diagonal (alias-diagonal))
        (success (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess)))
    (flet ((right-hand-side () (list (doubles 1d0 -2d0 3d0 5d0) (doubles 0d0 0d0 0d0 0d0)))
           (short-right-hand-side () (list (doubles 1d0 -2d0 3d0) (doubles 0d0 0d0 0d0)))
           (diagonal-values () (copy-seq (csc:entries diagonal)))
           (matrix-arrays (matrix)
             (list (csc:column-pointers matrix) (csc:row-indices matrix) (csc:entries matrix))))
      (list
       (list "coalton-linalg/conjugate-gradient:cg!" #'right-hand-side
             (lambda (a) (cg-caller:cg-code (cg:make-cg 4 50) (rf-caller:matrix-operator definite)
                                            (first a) (second a) 1d-12))
             '(((0 1) :same)))
       (list "coalton-linalg/conjugate-gradient:pcg!" #'right-hand-side
             (lambda (a) (cg-caller:pcg-code (cg:make-cg 4 50) (rf-caller:matrix-operator definite)
                                             (cg-caller:scaled-preconditioner 0.5d0) (first a) (second a) 1d-12))
             '(((0 1) :same)))
       (list "coalton-linalg/lsqr:lsqr!" #'right-hand-side
             (lambda (a) (ls-caller:lsqr-code (lsqr:make-lsqr 4 4 50) (rf-caller:matrix-operator unsymmetric)
                                              (ls-caller:transpose-operator unsymmetric) (first a) (second a)
                                              0d0 1d-9 0d0))
             '(((0 1) :same)))
       (list "coalton-linalg/refinement:refine!" #'right-hand-side
             (lambda (a)
               (multiple-value-bind (symbolic values) (alias-arrow-symbolic)
                 (expected-status (ldlt:factor! symbolic values 0d0)
                                  (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess))
                 (rf-caller:refine-code (rf:make-refiner 4 20) (rf-caller:matrix-operator definite)
                                        (rf-caller:factor-solver (ldlt:symbolic-factor symbolic))
                                        (first a) (second a) 1d-12 1d0)))
             '(((0 1) :same)))
       (list "coalton-linalg/dense/indefinite-ldlt:solve!" #'short-right-hand-side
             (lambda (a)
               (let ((factor (il:make-factor 3)))
                 (assert (= 0 (il-caller:observe-factor factor (column-major (alias-indefinite)) 0d0)))
                 (alias-reordered "The indefinite matrix" (il:factor-permutation factor))
                 (il-caller:observe-solve factor (first a) (second a))))
             '(((0 1) :same)))
       (list "coalton-linalg/dense/orthogonal:solve!" #'short-right-hand-side
             (lambda (a)
               (let ((decomposition (og:make-decomposition 3 3)))
                 (assert (= 0 (og-caller:observe-factor decomposition (column-major (alias-columns)))))
                 (alias-reordered "The matrix of unequal columns" (og:factor-permutation decomposition))
                 (assert (= 0 (og-caller:observe-rank decomposition 3)))
                 (og-caller:observe-solve decomposition (first a) (second a))))
             '(((0 1) :same)))
       (list "coalton-linalg/sparse/spd-ldlt:solve!" #'right-hand-side
             (lambda (a)
               (multiple-value-bind (symbolic values) (alias-arrow-symbolic)
                 (ldlt:factor! symbolic values 0d0)
                 (alias-success (ldlt:solve! (ldlt:symbolic-factor symbolic) (first a) (second a)) success)))
             '(((0 1) :same)))
       (list "coalton-linalg/sparse/spd-ldlt:solve-cholesky!" #'right-hand-side
             (lambda (a)
               (multiple-value-bind (symbolic values) (alias-arrow-symbolic)
                 (ldlt:factor-cholesky! symbolic values 0d0)
                 (alias-success (ldlt:solve-cholesky! (ldlt:symbolic-factor symbolic) (first a) (second a)) success)))
             '(((0 1) :same)))
       (list "coalton-linalg/sparse/spd-ldlt:solve-skipping!" #'right-hand-side
             (lambda (a)
               (multiple-value-bind (symbolic values) (alias-arrow-symbolic)
                 (assert (= 0 (dr-caller:observe-skipping symbolic values 0d0 0d0)))
                 (dr-caller:observe-skipping-solve (ldlt:symbolic-factor symbolic) (first a) (second a))))
             '(((0 1) :same)))
       (list "coalton-linalg/sparse/lu:solve!" #'right-hand-side
             (lambda (a)
               (let ((handle (lu:make-lu unsymmetric)))
                 (assert (= 0 (lu-caller:observe-factor handle unsymmetric 1d0)))
                 (alias-reordered "The unsymmetric matrix, by rows," (lu:row-permutation handle))
                 (lu-caller:observe-solve handle (first a) (second a))))
             '(((0 1) :same)))
       (list "coalton-linalg/sparse/lu:solve-transpose!" #'right-hand-side
             (lambda (a)
               (let ((handle (lu:make-lu unsymmetric)))
                 (assert (= 0 (lu-caller:observe-factor handle unsymmetric 1d0)))
                 (alias-reordered "The unsymmetric matrix, by rows," (lu:row-permutation handle))
                 (lu-caller:observe-transpose-solve handle (first a) (second a))))
             '(((0 1) :same)))
       ;; A square matrix of order two, so that the two work vectors and the output have one length.
       (list "coalton-linalg/sparse/rank-deficient:largest-singular-bounds!"
             (lambda () (list (doubles 0d0 0d0) (doubles 0d0 0d0) (doubles 0d0 0d0)))
             (lambda (a) (rd-caller:bounds-code (csc-from-entries 2 2 '((0 0 3) (1 0 1) (1 1 4))) 10
                                                (first a) (second a) (third a)))
             '(((0 1) :refused) ((0 2) :refused) ((1 2) :refused)))
       ;; Three rows on three columns, the third the sum of the others, so that the residual and the step
       ;; have one length and a row is left out.
       (list "coalton-linalg/sparse/rank-deficient:step!" #'short-right-hand-side
             (lambda (a)
               (let* ((matrix (csc-from-entries 3 3 '((0 0 1) (0 1 5) (0 2 2) (1 1 1) (1 2 7) (2 0 1) (2 1 6) (2 2 9))))
                      (kept (rd:analyze matrix)))
                 (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
                 (assert (= 1 (rd:left-out-count kept)) () "The premise: one row is left out; ~D are." (rd:left-out-count kept))
                 (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
                 (rd-caller:step-code kept (first a) (second a))))
             '(((0 1) :same)))
       (list "coalton-linalg/sparse/damped:step!" #'short-right-hand-side
             (lambda (a)
               (let* ((matrix (alias-stored (alias-columns))) (plan (dm:analyze matrix)))
                 (dm:set-matrix! plan matrix)
                 (assert (= 0 (dm-caller:column-damping-code plan 0d0)))
                 (assert (= 0 (dm-caller:factor-code plan 0.25d0)))
                 (dm-caller:step-code plan (first a) (second a))))
             '(((0 1) :same)))
       (list "coalton-linalg/sparse/csc:matvec!" #'right-hand-side
             (lambda (a) (csc:matvec! definite (first a) (second a)) 0)
             '(((0 1) :refused)))
       (list "coalton-linalg/sparse/csc:transpose-matvec!" #'right-hand-side
             (lambda (a) (csc:transpose-matvec! definite (first a) (second a)) 0)
             '(((0 1) :refused)))
       ;; Three rows, one column, two entries: the column pointers and the row indices can then be
       ;; one array, [0, 2].
       (list "coalton-linalg/sparse/csc:from-csc"
             (lambda () (list (indices 0 2) (indices 0 1) (doubles 5d0 7d0)))
             (lambda (a) (values 0 (matrix-arrays (csc:from-csc 3 1 (first a) (second a) (third a)))))
             '(((0 1) :same)))
       (list "coalton-linalg/sparse/csc:from-coo"
             (lambda () (list (indices 0 1 2) (indices 2 1 0) (doubles 5d0 7d0 9d0)))
             (lambda (a) (values 0 (matrix-arrays
                                    (csc:from-coo 3 3 (first a) (second a) (third a)
                                                  (coalton:coalton coalton-linalg/sparse/csc:RejectDuplicates)))))
             '(((0 1) :same)))
       (list "coalton-linalg/sparse/csc:permute"
             (lambda () (list (indices 2 0 3 1) (indices 1 3 2 0)))
             (lambda (a) (values 0 (matrix-arrays (csc:permute definite (first a) (second a)))))
             '(((0 1) :same)))
       (list "coalton-linalg/sparse/gram:refill!"
             (lambda () (list (diagonal-values) (doubles 0d0 0d0 0d0)))
             (lambda (a) (gram-caller:raw-error (gram:analyze-ata diagonal) (first a) (second a)))
             '(((0 1) :refused)))
       ;; The output is longer than the input by the order of the augmented matrix, which is at
       ;; least 1 wherever there is an entry: no pair has one length, and the roster checks that.
       (list "coalton-linalg/sparse/augmented:refill!"
             (lambda () (list (diagonal-values) (make-array 9 :element-type 'double-float :initial-element 0d0)))
             (lambda (a) (qd-caller:raw-error (aug:analyze diagonal) (first a) 1d0 1d0 (second a)))
             '())
       (list "coalton-linalg/sparse/scaling:row-maxima!"
             (lambda () (list (diagonal-values) (doubles 0d0 0d0 0d0)))
             (lambda (a) (sc-caller:maxima-error diagonal (first a) (second a) t))
             '(((0 1) :refused)))
       (list "coalton-linalg/sparse/scaling:column-maxima!"
             (lambda () (list (diagonal-values) (doubles 0d0 0d0 0d0)))
             (lambda (a) (sc-caller:maxima-error diagonal (first a) (second a) nil))
             '(((0 1) :refused)))
       (list "coalton-linalg/sparse/scaling:equilibrate!"
             (lambda () (list (diagonal-values) (doubles 0d0 0d0 0d0) (doubles 0d0 0d0 0d0)))
             (lambda (a) (sc-caller:equilibrate-error diagonal (first a) (second a) (third a)))
             '(((0 1) :refused) ((0 2) :refused) ((1 2) :refused)))
       (list "coalton-linalg/sparse/scaling:scale!"
             (lambda () (list (diagonal-values) (doubles 0.5d0 2d0 4d0) (doubles 8d0 0.25d0 1d0) (doubles 0d0 0d0 0d0)))
             (lambda (a) (sc-caller:scale-error diagonal (first a) (second a) (third a) (fourth a)))
             '(((0 1) :same) ((0 2) :same) ((0 3) :refused) ((1 2) :same) ((1 3) :refused) ((2 3) :refused)))
       (list "coalton-linalg/sparse/scaling:scale-vector!"
             (lambda () (list (doubles 0.5d0 2d0 4d0) (doubles 3d0 -8d0 0.25d0) (doubles 0d0 0d0 0d0)))
             (lambda (a) (sc-caller:vector-error (first a) (second a) (third a)))
             '(((0 1) :same) ((0 2) :refused) ((1 2) :same)))
       (list "coalton-linalg/sparse/scaling:equilibrate-symmetric!"
             (lambda () (list (diagonal-values) (doubles 0d0 0d0 0d0)))
             (lambda (a) (sc-caller:symmetric-code diagonal (first a) (second a)))
             '(((0 1) :refused)))
       (list "coalton-linalg/sparse/scaling:scale-symmetric!"
             (lambda () (list (diagonal-values) (doubles 0.5d0 2d0 4d0) (doubles 0d0 0d0 0d0)))
             (lambda (a) (sc-caller:symmetric-scale-error diagonal (first a) (second a) (third a)))
             '(((0 1) :same) ((0 2) :refused) ((1 2) :refused)))
       ;; Not public: the package is the library's own, and its one caller gives it the two arrays
       ;; of a pattern it has just validated as a pattern, which one array cannot be both of.
       (list "coalton-linalg/sparse/ordering:minimum-degree" :exempt
             "an internal package; its only caller passes the two arrays of a validated pattern")))))

(defun alias-missing (inventory keys)
  "Members of the inventory the roster does not name, and names in the roster that are not members."
  (values (set-difference inventory keys :test #'string=) (set-difference keys inventory :test #'string=)))

(defun alias-check-entry (entry)
  "Run every pair of one entry and compare with what it declares. Returns the outcomes as a list of (pair outcome)."
  (destructuring-bind (key build call expected) entry
    (assert (= 0 (alias-run call (funcall build))) () "~A: the entry's own call is refused." key)
    (let* ((pairs (alias-pairs (funcall build)))
           (outcomes (loop :for (i j) :in pairs :collect (list (list i j) (alias-verdict key build call i j)))))
      (assert (equal outcomes expected) ()
              "~A: the pairs that can be one array behave as ~S; the roster declares ~S." key outcomes expected)
      outcomes)))

(defun run-aliasing-tests ()
  (let* ((roster (alias-roster)) (inventory (alias-inventory)) (keys (mapcar #'first roster))
         (same 0) (refused 0) (exempt 0) (without-pair 0))
    ;; The roster names every member the sources have, and nothing else.
    (assert (> (length inventory) 20) () "The inventory found only ~D calls: ~S." (length inventory) inventory)
    (multiple-value-bind (unanswered unknown) (alias-missing inventory keys)
      (assert (null unanswered) () "Calls that take two or more arrays and have no entry in the aliasing roster: ~S."
              unanswered)
      (assert (null unknown) () "Entries of the aliasing roster that name no such call: ~S." unknown))
    ;; Control: a call the sources would add is reported until it is answered for.
    (assert (equal '("coalton-linalg/new:call!") (alias-missing (cons "coalton-linalg/new:call!" inventory) keys)))
    (dolist (entry roster)
      (cond
        ((eq (second entry) :exempt) (incf exempt))
        (t
         (let ((outcomes (alias-check-entry entry)))
           (when (null outcomes) (incf without-pair))
           (incf same (count :same outcomes :key #'second))
           (incf refused (count :refused outcomes :key #'second))))))
    ;; Controls for what a refusal is: the library's typed error, and nothing else that can be signalled.
    (assert (= 99 (alias-run (lambda (arrays)
                               (declare (ignore arrays))
                               (rows:analyze-with-limit (csc-from-entries 2 2 '((0 0 1) (1 0 1) (1 1 1))) 1))
                             '())))
    (assert (handler-case (progn (alias-run (lambda (arrays) (declare (ignore arrays)) (error "an index out of range")) '())
                                 nil)
              (coalton-linalg/errors:linalgerror () nil)
              (error () t))
            () "A condition that is no refusal was taken for one.")
    ;; Controls for the verdict itself: a call that reads what it has already overwritten, and
    ;; a refusal that comes after a write, are both reported.
    (flet ((reported-p (call)
             (handler-case (progn (alias-verdict "control" (lambda () (list (doubles 1d0 2d0 3d0) (doubles 0d0 0d0 0d0)))
                                                 call 0 1)
                                  nil)
               (error () t))))
      (assert (reported-p (lambda (a)
                            ;; out[i] = in[i - 1] + in[i]: right with two arrays, wrong with one.
                            (let ((in (first a)) (out (second a)))
                              (dotimes (i 3) (setf (aref out i) (+ (aref in i) (if (plusp i) (aref in (1- i)) 0d0)))))
                            0))
              () "The verdict accepted a call that reads what it has overwritten.")
      (assert (reported-p (lambda (a)
                            (setf (aref (second a) 0) 42d0)
                            (if (eq (first a) (second a)) 9 0)))
              () "The verdict accepted a refusal that came after a write.")
      (assert (not (reported-p (lambda (a) (replace (second a) (first a)) 0)))
              () "The verdict reported a call that is right with one array."))
    (format t "ALIASING-ROSTER-PASS calls=~D exempt=~D without-a-pair=~D pairs-same=~D pairs-refused=~D~%"
            (length roster) exempt without-pair same refused)))
