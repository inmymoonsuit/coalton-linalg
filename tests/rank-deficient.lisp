(in-package #:coalton-linalg/tests)

;;; Kept rows, their vectors kept from call to call, and the certificate. Every statement is
;;; checked against exact rational arithmetic on the stored doubles: a bound against the exact
;;; norm of the exact product, a certificate against the exact inertia of the kept rows' Gram
;;; matrix less the square of what it claims.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "RD" "COALTON-LINALG/SPARSE/RANK-DEFICIENT")
  (sb-ext:add-package-local-nickname "RD-CALLER" "COALTON-LINALG/TESTS/RANK-DEFICIENT-CALLERS"))

(defun rd-redundant-entries (state bases columns extras)
  "The stored entries of a matrix of small integers: BASES seeded rows, and EXTRAS rows each a combination of one
or two of them with small integer weights, the rows then shuffled. Returns the entries and the row count."
  (let ((rows '()))
    (dotimes (b bases)
      (let ((row (make-array columns :initial-element 0)))
        (dotimes (j columns)
          (when (< (random 100 state) 60)
            (setf (aref row j) (let ((v (- (random 7 state) 3))) (if (zerop v) 2 v)))))
        (when (every #'zerop row) (setf (aref row (random columns state)) 1))
        (push row rows)))
    (let ((base (coerce (reverse rows) 'vector)))
      (dotimes (e extras)
        (let ((row (make-array columns :initial-element 0)))
          (dotimes (term (1+ (random 2 state)))
            (let ((weight (nth (random 4 state) '(-2 -1 1 2))) (source (aref base (random bases state))))
              (dotimes (j columns) (incf (aref row j) (* weight (aref source j))))))
          (when (every #'zerop row) (replace row (aref base (random bases state))))
          (push row rows))))
    (let ((all (coerce rows 'vector)))
      (loop :for i :from (1- (length all)) :downto 1
            :do (rotatef (aref all i) (aref all (random (1+ i) state))))
      (values (loop :for i :from 0 :for row :across all
                    :append (loop :for j :below columns :unless (zerop (aref row j)) :collect (list i j (aref row j))))
              (length all)))))

(defun rd-vector (kept position m)
  "The kept vector at a position as a dense vector: 1 at its own row, its values at kept rows."
  (let ((v (zeros m)) (start (rd:multiplier-start kept position)) (rows (rd:multiplier-rows kept))
        (values (rd:multiplier-values kept)))
    (setf (aref v (aref (rd:left-out-rows kept) position)) 1d0)
    (dotimes (e (rd:multiplier-count kept position) v)
      (setf (aref v (aref rows (+ start e))) (aref values (+ start e))))))

(defun rd-check-bounds (label kept matrix a m)
  "Every standing bound holds for the exact product with the vector as kept, is no looser than the product bound
over the whole matrix, and the gap bound holds for their sum. Returns the sum of the exact squares."
  (let ((total 0))
    (dotimes (p (rd:left-out-count kept))
      (let* ((v (rd-vector kept p m)) (exact (dr-product-squared a v)) (bound (rd:residual-bound kept p)))
        (incf total exact)
        (assert (and (plusp bound) (<= exact (expt (rational bound) 2))) ()
                "~A, row ~D: the product's norm is ~E and its bound ~E." label p (sqrt (float exact 1d0)) bound)
        (let ((whole (rows:transpose-product-bound matrix v)))
          (assert (<= bound (* whole (+ 1d0 1d-12))) ()
                  "~A, row ~D: the bound over the rows touched, ~E, exceeds the one over the whole matrix, ~E." label p
                  bound whole))))
    (let ((gap (rd:gap-bound kept)))
      (assert (<= total (expt (rational gap) 2)) () "~A: the gap bound ~E is below the root of the sum, ~E." label gap
              (sqrt (float total 1d0))))
    total))

(defun rd-kept-rows (kept m)
  (let ((out (coerce (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)) 'list)))
    (loop :for i :below m :unless (member i out) :collect i)))

(defun rd-proved-p (a kept m bound)
  "Whether every singular value of the kept rows exceeds BOUND, exactly: their Gram matrix less its square is
positive definite."
  (let* ((kept-rows (rd-kept-rows kept m)) (r (length kept-rows)))
    (multiple-value-bind (positive negative zero) (exact-inertia (dr-gram-shifted (dr-rows a kept-rows) (expt (rational bound) 2)))
      (declare (ignore negative zero))
      (= positive r))))

(defun run-rank-deficient-exact-tests ()
  (let ((state (sb-ext:seed-random-state 20261123)) (cases 0) (vectors 0) (certified 0) (refused 0) (kinds (list 0 0 0 0 0 0)))
    (loop :for (bases columns extras) :in '((2 4 1) (3 5 2) (4 6 3) (5 5 4) (6 9 5) (4 4 2) (7 10 6))
          :do (dotimes (repeat 6)
                (multiple-value-bind (entries m) (rd-redundant-entries state bases columns extras)
                  (let* ((label (format nil "~D+~D rows by ~D #~D" bases extras columns repeat))
                         (matrix (csc-from-entries m columns entries)) (a (dr-dense m columns entries))
                         (kept (rd:analyze matrix)) (rank (dr-rank a)))
                    (incf cases)
                    (assert (not (rd:kept-found? kept)))
                    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)) () "~A: the search gave ~D." label
                            (rd-caller:find-code kept matrix 1d-6 1d-12))
                    (assert (and (rd:kept-found? kept) (= rank (rd:kept-count kept))
                                 (= (- m rank) (rd:left-out-count kept)))
                            () "~A: kept ~D rows of exact rank ~D." label (rd:kept-count kept) rank)
                    (incf vectors (rd:left-out-count kept))
                    (rd-check-bounds label kept matrix a m)
                    ;; The combinations are the exact ones to rounding: the gap is tiny beside the matrix.
                    (assert (<= (rd:gap-bound kept) 1d-10) () "~A: the gap bound is ~E." label (rd:gap-bound kept))
                    ;; A vector cut harder is another vector, and its bound is still a bound.
                    (let ((before (loop :for p :below (rd:left-out-count kept) :sum (rd:multiplier-count kept p))))
                      (dolist (drop '(1d-3 0.5d0 0.99d0))
                        (assert (= 0 (rd-caller:find-code kept matrix 1d-6 drop)))
                        (assert (<= (loop :for p :below (rd:left-out-count kept) :sum (rd:multiplier-count kept p)) before))
                        (rd-check-bounds (format nil "~A, drop ~S" label drop) kept matrix a m)))
                    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
                    ;; The certificate, over a ladder of floors: what it issues is true exactly, and
                    ;; below the bound it proves first it is always issued.
                    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)) () "~A: no certificate at 1e-9: ~D." label
                            (rd-caller:certify-code kept 1d-9 1d-9))
                    (let ((first (rd:proved-lower-bound kept)))
                      (assert (and (rd:certified? kept) (>= first 1d-9) (rd-proved-p a kept m first)) ()
                              "~A: the bound proved, ~E, is not true of the kept rows." label first)
                      (loop :for power :from -9 :to 3
                            :do (let* ((floor (expt 10d0 power)) (code (rd-caller:certify-code kept 1d-9 floor)))
                                  (incf (nth code kinds))
                                  (cond
                                    ((= code 0)
                                     (incf certified)
                                     (assert (and (>= (rd:proved-lower-bound kept) floor)
                                                  (rd-proved-p a kept m (rd:proved-lower-bound kept)))
                                             () "~A: certified at ~E with ~E, which is false." label floor
                                             (rd:proved-lower-bound kept)))
                                    (t
                                     (incf refused)
                                     (assert (member code '(2 3)) () "~A at ~E: refusal ~D." label floor code)
                                     (assert (and (not (rd:certified? kept)) (= -1d0 (rd:proved-lower-bound kept))))
                                     (assert (> floor (* 0.2d0 first)) ()
                                             "~A: refused at ~E, though ~E was proved before." label floor first))))))
                    ;; A lower floor the gap bound does not come under.
                    (when (plusp (rd:left-out-count kept))
                      (assert (= 1 (rd-caller:certify-code kept (* 0.5d0 (rd:gap-bound kept)) 1d-9))))))))
    (assert (and (> certified 100) (> refused 100)) () "The floors were ~D certified and ~D refused." certified refused)
    (format t "RANK-DEFICIENT-EXACT-PASS cases=~D vectors=~D certified=~D refused=~D by-status=~S~%" cases vectors certified
            refused kinds)))

(defun run-rank-deficient-reprove-tests ()
  "Bounds proved again for other values are bounds for those values and those vectors, whatever the values are."
  (let ((state (sb-ext:seed-random-state 99)) (cases 0) (outdated 0))
    (loop :for (bases columns extras) :in '((3 5 2) (5 7 4) (6 6 3))
          :do (dotimes (repeat 5)
                (multiple-value-bind (entries m) (rd-redundant-entries state bases columns extras)
                  (let* ((matrix (csc-from-entries m columns entries)) (a (dr-dense m columns entries))
                         (kept (rd:analyze matrix)))
                    (assert (= 1 (rd-caller:reprove-code kept matrix)))
                    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
                    (let ((bounds (loop :for p :below (rd:left-out-count kept) :collect (rd:residual-bound kept p)))
                          (gap (rd:gap-bound kept)) (starts (loop :for p :to (rd:left-out-count kept) :collect (rd:multiplier-start kept p)))
                          (values (copy-seq (rd:multiplier-values kept))))
                      ;; The same values: the same bounds, bit for bit, and the vectors untouched.
                      (assert (= 0 (rd-caller:reprove-code kept matrix)))
                      (assert (and (equal bounds (loop :for p :below (rd:left-out-count kept) :collect (rd:residual-bound kept p)))
                                   (= gap (rd:gap-bound kept))))
                      ;; Every value moved by a part in a million: the vectors are no longer the
                      ;; combinations, the bounds say by how much, and they hold exactly.
                      (let* ((moved (mapcar (lambda (entry)
                                              (list (first entry) (second entry)
                                                    (* (third entry) (+ 1d0 (* 1d-6 (- (random 2d0 state) 1d0))))))
                                            entries))
                             (other (csc-from-entries m columns moved)) (b (dr-dense m columns moved)))
                        (incf cases)
                        (assert (= 0 (rd-caller:reprove-code kept other)))
                        (assert (equal starts (loop :for p :to (rd:left-out-count kept) :collect (rd:multiplier-start kept p))))
                        (assert (same-bits-p values (rd:multiplier-values kept)) () "A reproof changed a vector.")
                        (rd-check-bounds "moved values" kept other b m)
                        (when (plusp (rd:left-out-count kept))
                          ;; Out of date on purpose: the certificate that held is refused on the gap.
                          (assert (> (rd:gap-bound kept) (* 1d3 gap)))
                          (assert (= 1 (rd-caller:certify-code kept 1d-9 1d-9)) () "Vectors a part in a million out of date were certified at 1e-9.")
                          (incf outdated)
                          ;; And at a floor above what they now leave, the certificate is issued and true.
                          (let ((floor (* 4d0 (rd:gap-bound kept))))
                            (when (= 0 (rd-caller:certify-code kept floor floor))
                              (assert (rd-proved-p b kept m (rd:proved-lower-bound kept))))))
                        ;; Back to the first values: the first bounds.
                        (assert (= 0 (rd-caller:reprove-code kept matrix)))
                        (assert (= gap (rd:gap-bound kept)))
                        (rd-check-bounds "first values again" kept matrix a m)))))))
    (format t "RANK-DEFICIENT-REPROVE-PASS cases=~D out-of-date-refused=~D~%" cases outdated)))

(defun run-rank-deficient-scale-tests ()
  "The same matrix at another scale gives the same rows, the same vectors bit for bit, and bounds and a proved
lower bound that are the first ones times the scale."
  (let* ((entries '((0 0 1) (0 1 2) (1 1 1) (1 2 3) (2 0 2) (2 2 1) (3 0 1) (3 1 3) (3 2 3) (4 0 3) (4 1 2) (4 2 4)))
         (m 5) (n 3) (cases 0))
    (flet ((at (power)
             (let* ((matrix (csc-from-entries m n (mapcar (lambda (entry)
                                                            (list (first entry) (second entry)
                                                                  (scale-float (coerce (third entry) 'double-float) power)))
                                                          entries)))
                    (kept (rd:analyze matrix)) (floor (scale-float 1d-6 power)))
               (assert (= 0 (rd-caller:find-code kept matrix floor 1d-12)) () "At 2^~D the search gave ~D." power
                       (rd-caller:find-code kept matrix floor 1d-12))
               (assert (= 0 (rd-caller:certify-code kept floor floor)) () "At 2^~D the certificate gave ~D." power
                       (rd-caller:certify-code kept floor floor))
               (list (rd:left-out-count kept) (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept))
                     (copy-seq (rd:multiplier-rows kept)) (copy-seq (rd:multiplier-values kept))
                     (loop :for p :below (rd:left-out-count kept) :collect (rd:residual-bound kept p))
                     (rd:gap-bound kept) (rd:proved-lower-bound kept)))))
      (let ((reference (at 0)))
        (assert (= 2 (first reference)))
        (dolist (power '(-900 -600 -100 100 600 900))
          (let ((scaled (at power)))
            (incf cases)
            (assert (and (equalp (subseq reference 0 3) (subseq scaled 0 3)) (same-bits-p (fourth reference) (fourth scaled)))
                    () "At 2^~D the rows or the vectors differ." power)
            (assert (and (equal (mapcar (lambda (b) (scale-float b power)) (fifth reference)) (fifth scaled))
                         (= (scale-float (sixth reference) power) (sixth scaled))
                         (= (scale-float (seventh reference) power) (seventh scaled)))
                    () "At 2^~D the bounds are ~S, ~S and ~S; scaled from scale 1 they are ~S, ~S and ~S." power
                    (fifth scaled) (sixth scaled) (seventh scaled)
                    (mapcar (lambda (b) (scale-float b power)) (fifth reference)) (scale-float (sixth reference) power)
                    (scale-float (seventh reference) power))))))
    (format t "RANK-DEFICIENT-SCALE-PASS cases=~D~%" cases)))

(defun run-singular-bounds-tests ()
  "The two bounds on the largest singular value hold exactly: the square of the upper less J^T J has no negative
eigenvalue, and J^T J less the square of the lower has one that is not negative."
  (let ((state (sb-ext:seed-random-state 5)) (cases 0) (widest 0d0))
    (loop :for (m n) :in '((1 1) (2 3) (4 4) (5 3) (6 8) (9 9))
          :do (dotimes (repeat 6)
                (let* ((entries (seeded-entries state m n 55 :hub t)) (matrix (csc-from-entries m n entries))
                       (a (dr-dense m n entries)) (out (zeros 2))
                       (gram (build-matrix n (lambda (p q) (loop :for i :below m :sum (* (aref a i p) (aref a i q)))))))
                  (assert (= 0 (rd-caller:bounds-code matrix 30 (zeros n) (zeros m) out)))
                  (let ((lower (aref out 0)) (upper (aref out 1)))
                    (incf cases)
                    (assert (and (plusp lower) (<= lower upper)) () "The bounds are ~S and ~S." lower upper)
                    (setf widest (max widest (/ upper lower)))
                    (multiple-value-bind (positive negative zero)
                        (exact-inertia (build-matrix n (lambda (p q) (- (if (= p q) (expt (rational upper) 2) 0) (aref gram p q)))))
                      (declare (ignore positive zero))
                      (assert (= 0 negative) () "~D by ~D: the upper bound ~S is below the largest singular value." m n upper))
                    (multiple-value-bind (positive negative zero)
                        (exact-inertia (build-matrix n (lambda (p q) (- (aref gram p q) (if (= p q) (expt (rational lower) 2) 0)))))
                      (declare (ignore negative))
                      (assert (plusp (+ positive zero)) () "~D by ~D: the lower bound ~S is above the largest singular value." m n
                              lower))))))
    ;; No step at all still gives two bounds; a matrix of zeros gives a lower bound of zero.
    (let ((out (zeros 2)) (matrix (csc-from-entries 2 2 '((0 0 3) (1 1 4)))))
      (assert (= 0 (rd-caller:bounds-code matrix 0 (zeros 2) (zeros 2) out)))
      (assert (<= (aref out 0) 4d0 (aref out 1)))
      (assert (= 0 (rd-caller:bounds-code matrix 30 (zeros 2) (zeros 2) out)))
      (assert (<= (* 4d0 (- 1d0 1d-9)) (aref out 0) 4d0 (aref out 1)) () "A diagonal matrix's bounds are ~S." out)
      (assert (= 0 (rd-caller:bounds-code (csc-from-entries 2 2 '((0 0 0))) 5 (zeros 2) (zeros 2) out)))
      (assert (and (= 0d0 (aref out 0)) (<= 0d0 (aref out 1) 1d-300)))
      ;; Lengths and shared arrays are typed errors.
      (assert (= 92 (rd-caller:bounds-code matrix 5 (zeros 3) (zeros 2) out)))
      (assert (= 92 (rd-caller:bounds-code matrix 5 (zeros 2) (zeros 2) (zeros 1))))
      (let ((shared (zeros 2)))
        (assert (= 93 (rd-caller:bounds-code matrix 5 shared shared out)))
        (assert (= 93 (rd-caller:bounds-code matrix 5 shared (zeros 2) shared)))
        (assert (= 93 (rd-caller:bounds-code matrix 5 (zeros 2) shared shared)))))
    (assert (<= widest 4d0) () "The widest pair of bounds is a factor ~S apart." widest)
    (format t "SINGULAR-BOUNDS-PASS cases=~D widest=~,3F~%" cases widest)))

(defun run-rank-deficient-outcome-tests ()
  (let* ((entries '((0 0 1) (0 1 2) (1 1 1) (1 2 3) (2 0 2) (2 2 1) (3 0 1) (3 1 3) (3 2 3) (4 0 3) (4 1 2) (4 2 4)))
         (matrix (csc-from-entries 5 3 entries)) (kept (rd:analyze matrix))
         (nan (not-a-number)) (infinity sb-ext:double-float-positive-infinity))
    (assert (and (= 5 (rd:row-count kept)) (= 3 (rd:column-count kept))))
    ;; Before any search.
    (assert (and (not (rd:kept-found? kept)) (= 0 (rd:kept-count kept)) (= 0 (rd:left-out-count kept))
                 (= -1d0 (rd:gap-bound kept)) (= -1d0 (rd:residual-bound kept 0)) (= -1d0 (rd:proved-lower-bound kept))
                 (= 0 (rd:multiplier-count kept 0)) (= 0 (rd:multiplier-start kept 0))))
    (assert (= 1 (rd-caller:reprove-code kept matrix)))
    (assert (= 5 (rd-caller:certify-code kept 1d-9 1d-9)))
    ;; A drop tolerance that is not between 0 and 1, and the search's own refusals, leave nothing.
    (dolist (drop (list 0d0 1d0 -1d-3 2d0 nan))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
      (assert (= 4 (rd-caller:find-code kept matrix 1d-6 drop)) () "A drop tolerance of ~S was taken." drop)
      (assert (and (not (rd:kept-found? kept)) (= -1d0 (rd:gap-bound kept)) (= 5 (rd-caller:certify-code kept 1d-9 1d-9))
                   (= 0 (rd:kept-count kept)) (= 0 (rd:left-out-count kept)) (= 0 (rd:multiplier-count kept 0)))))
    (dolist (threshold (list -1d0 nan infinity))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
      (assert (= 3 (rd-caller:find-code kept matrix threshold 1d-12)))
      (assert (not (rd:kept-found? kept))))
    (let ((mixed (csc-from-entries 5 3 (mapcar (lambda (entry)
                                                 (list (first entry) (second entry)
                                                       (if (= 4 (first entry))
                                                           (scale-float (coerce (third entry) 'double-float) -600)
                                                           (third entry))))
                                               entries))))
      (assert (= 2 (rd-caller:find-code kept mixed 0d0 1d-12)))
      (assert (not (rd:kept-found? kept))))
    ;; A stand: two rows left out of five on three columns.
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
    (assert (and (= 3 (rd:kept-count kept)) (= 2 (rd:left-out-count kept)) (plusp (rd:gap-bound kept))
                 (= -1d0 (rd:residual-bound kept 2)) (= 0 (rd:multiplier-count kept 2))))
    ;; Floors that cannot be used.
    (dolist (pair (list (list 0d0 1d-9) (list -1d0 1d-9) (list 1d-9 infinity) (list nan 1d-9) (list 1d-9 nan) (list 1d-3 1d-9)))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= 4 (rd-caller:certify-code kept (first pair) (second pair))) () "Floors ~S were taken." pair)
      (assert (and (not (rd:certified? kept)) (= -1d0 (rd:proved-lower-bound kept)))))
    ;; A floor above every singular value of the kept rows, by a factorization accurate there.
    (assert (= 2 (rd-caller:certify-code kept 1d-9 1d3)))
    ;; A certificate does not survive a reproof or a search.
    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
    (assert (rd:certified? kept))
    (assert (= 0 (rd-caller:reprove-code kept matrix)))
    (assert (and (not (rd:certified? kept)) (= -1d0 (rd:proved-lower-bound kept))))
    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
    (assert (and (not (rd:certified? kept)) (= -1d0 (rd:proved-lower-bound kept))))
    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
    (assert (= 4 (rd-caller:find-code kept matrix 1d-6 2d0)))
    (assert (and (not (rd:certified? kept)) (= -1d0 (rd:proved-lower-bound kept))))
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
    ;; A matrix of another pattern or shape is a typed error from a reproof and from a search.
    (assert (= 91 (rd-caller:reprove-code kept (csc-from-entries 5 3 (cons '(0 2 1) entries)))))
    (assert (= 92 (rd-caller:reprove-code kept (csc-from-entries 6 3 entries))))
    (assert (= 91 (rd-caller:find-code kept (csc-from-entries 5 3 (cons '(0 2 1) entries)) 1d-6 1d-12)))
    ;; Two kept rows 4e-8 apart: their smallest singular value is about 2e-8. A floor of 1e-6 is
    ;; above it, and the certificate says the rows were not proved above the floor. A floor of
    ;; 1e-8 and below is below it, and true, but the proof's own error at such a floor is larger
    ;; than the floor: that is said by a name of its own, and is neither a certificate nor a
    ;; claim about the rows. Rows 1e-6 apart are resolved, and what is proved of them is exactly true.
    (flet ((pair (gap) (list '(0 0 1) '(0 1 1) '(1 0 1) (list 1 1 (+ 1d0 gap)))))
      (let* ((close (csc-from-entries 2 2 (pair 4d-8))) (kept (rd:analyze close)))
        (assert (= 0 (rd-caller:find-code kept close 0d0 1d-12)))
        (assert (= 2 (rd:kept-count kept)) () "The premise: rows 4e-8 apart are both kept; ~D were." (rd:kept-count kept))
        (assert (= 2 (rd-caller:certify-code kept 1d-300 1d-6)))
        (dolist (floor '(1d-8 1d-9 1d-10 1d-12))
          (assert (= 3 (rd-caller:certify-code kept 1d-300 floor)) () "Rows 4e-8 apart at a floor of ~E: ~D." floor
                  (rd-caller:certify-code kept 1d-300 floor)))
        (assert (and (not (rd:certified? kept)) (= -1d0 (rd:proved-lower-bound kept)))))
      (let* ((apart (csc-from-entries 2 2 (pair 1d-6))) (kept (rd:analyze apart)))
        (assert (= 0 (rd-caller:find-code kept apart 0d0 1d-12)))
        (assert (= 0 (rd-caller:certify-code kept 1d-300 1d-10)))
        (assert (and (> (rd:proved-lower-bound kept) 1d-7)
                     (rd-proved-p (dr-dense 2 2 (pair 1d-6)) kept 2 (rd:proved-lower-bound kept)))
                () "Rows 1e-6 apart: the bound proved is ~S." (rd:proved-lower-bound kept))))
    ;; No row left out: a gap bound of zero, and a certificate on the kept rows alone.
    (let* ((full (csc-from-entries 2 2 '((0 0 2) (1 0 1) (1 1 3)))) (kept (rd:analyze full)))
      (assert (= 0 (rd-caller:find-code kept full 1d-6 1d-12)))
      (assert (and (= 0 (rd:left-out-count kept)) (= 0d0 (rd:gap-bound kept))))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (rd-proved-p (dr-dense 2 2 '((0 0 2) (1 0 1) (1 1 3))) kept 2 (rd:proved-lower-bound kept))))
    ;; One row on one column. Its singular value is 5; the ladder of bounds tried starts at half
    ;; of one, so a floor of 3 is proved only by trying the floor itself, and is then what is
    ;; proved; a floor of 6 is above the row.
    (let* ((one (csc-from-entries 1 1 '((0 0 5)))) (kept (rd:analyze one)))
      (assert (= 0 (rd-caller:find-code kept one 1d-6 1d-12)))
      (assert (and (= 1 (rd:kept-count kept)) (= 0d0 (rd:gap-bound kept)) (= 0 (rd-caller:certify-code kept 1d-9 1d-9))
                   (<= 1d-9 (rd:proved-lower-bound kept) 5d0))
              () "One entry of 5: kept ~D, proved ~S." (rd:kept-count kept) (rd:proved-lower-bound kept))
      (assert (and (= 0 (rd-caller:certify-code kept 1d-9 3d0)) (= 3d0 (rd:proved-lower-bound kept)))
              () "One entry of 5 at a floor of 3: proved ~S." (rd:proved-lower-bound kept))
      (assert (= 2 (rd-caller:certify-code kept 1d-9 6d0))))
    ;; The drop is relative to each vector's own largest entry at a kept row. Row 2 is row 0
    ;; over 1024, and row 3 is that plus row 1. Cut at a half, row 2's vector keeps its one entry
    ;; of 1/1024, and row 3's keeps the entry of 1 and loses the one of 1/1024; what is lost is
    ;; row 0 over 1024, whose length is the root of 5, and the bound says so.
    (let* ((entries '((0 0 1024) (0 1 2048) (1 1 1) (1 2 1) (2 0 1) (2 1 2) (3 0 1) (3 1 3) (3 2 1)))
           (matrix (csc-from-entries 4 3 entries)) (a (dr-dense 4 3 entries)) (kept (rd:analyze matrix)))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
      (assert (and (equalp #(2 3) (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)))
                   (equal '(1 2) (loop :for p :below 2 :collect (rd:multiplier-count kept p))))
              () "The premise: rows 2 and 3 are left out with one and two entries; rows ~S with ~S were."
              (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept))
              (loop :for p :below (rd:left-out-count kept) :collect (rd:multiplier-count kept p)))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 0.5d0)))
      (assert (and (equal '(1 1) (loop :for p :below 2 :collect (rd:multiplier-count kept p)))
                   (= (/ -1d0 1024) (aref (rd:multiplier-values kept) (rd:multiplier-start kept 0)))
                   (= 1 (aref (rd:multiplier-rows kept) (rd:multiplier-start kept 1)))
                   (= -1d0 (aref (rd:multiplier-values kept) (rd:multiplier-start kept 1))))
              () "Cut at a half the vectors hold ~S entries, ~S at rows ~S."
              (loop :for p :below 2 :collect (rd:multiplier-count kept p)) (subseq (rd:multiplier-values kept) 0 2)
              (subseq (rd:multiplier-rows kept) 0 2))
      (rd-check-bounds "a relative drop" kept matrix a 4)
      (assert (and (<= (rd:residual-bound kept 0) 1d-12) (<= (sqrt 5d0) (rd:residual-bound kept 1) (* (sqrt 5d0) (+ 1d0 1d-12))))
              () "Cut at a half the bounds are ~S and ~S." (rd:residual-bound kept 0) (rd:residual-bound kept 1)))
    ;; The proof's margin. Two kept rows, a rotation times diag(1, 1/16), so that their smallest
    ;; singular value is 1/16 to rounding, above it or below; and 1/16 is the first bound the
    ;; ladder tries. A factorization at exactly that bound cannot tell the two apart, so what is
    ;; proved must be a bound the margin covers, and must be true.
    (let ((below 0))
      (dolist (triple '((3 4 5) (5 12 13) (8 15 17) (7 24 25) (20 21 29) (12 35 37) (9 40 41) (28 45 53) (11 60 61) (16 63 65)))
        (dolist (swap '(nil t))
          (let* ((c (/ (coerce (if swap (second triple) (first triple)) 'double-float) (third triple)))
                 (s (/ (coerce (if swap (first triple) (second triple)) 'double-float) (third triple)))
                 (entries (list (list 0 0 c) (list 0 1 (- (* s 0.0625d0))) (list 1 0 s) (list 1 1 (* c 0.0625d0))))
                 (matrix (csc-from-entries 2 2 entries)) (a (dr-dense 2 2 entries)) (kept (rd:analyze matrix)))
            (assert (= 0 (rd-caller:find-code kept matrix 1d-9 1d-12)))
            (assert (and (= 2 (rd:kept-count kept)) (= 0 (rd-caller:certify-code kept 1d-9 1d-9))))
            (unless (rd-proved-p a kept 2 0.0625d0) (incf below))
            (assert (rd-proved-p a kept 2 (rd:proved-lower-bound kept)) ()
                    "Rows with a smallest singular value at 1/16 to rounding, ~S ~S: ~S was proved and is false." triple swap
                    (rd:proved-lower-bound kept)))))
      (assert (>= below 4) () "The premise: some of the twenty have their singular value at or below 1/16; ~D do." below))
    ;; A row with no entry and a row of stored zeros beside real rows: both are left out with a
    ;; vector of one entry, their own, and bounds that are true of a product that is exactly zero.
    (let* ((entries '((0 0 1) (0 1 2) (2 0 0) (2 1 0) (3 0 2) (3 1 4) (4 1 5)))
           (matrix (csc-from-entries 5 2 entries)) (a (dr-dense 5 2 entries)) (kept (rd:analyze matrix)))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
      (assert (and (= 2 (rd:kept-count kept)) (equalp #(1 2 3) (subseq (rd:left-out-rows kept) 0 3))
                   (equal '(0 0 1) (loop :for p :below 3 :collect (rd:multiplier-count kept p))))
              () "Rows left out ~S with ~S entries." (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept))
              (loop :for p :below (rd:left-out-count kept) :collect (rd:multiplier-count kept p)))
      (rd-check-bounds "an empty row and a row of zeros" kept matrix a 5)
      (assert (and (<= (rd:residual-bound kept 0) 1d-300) (<= (rd:residual-bound kept 1) 1d-150)))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (rd-proved-p a kept 5 (rd:proved-lower-bound kept))))
    ;; A matrix of stored zeros: no row is kept, the statement about the kept rows is empty, and
    ;; the certificate is the gap bound's alone. Proved again for a matrix that is not zero, the
    ;; vectors bound nothing below its entries, and the certificate is refused on the gap.
    (let* ((zero (csc-from-entries 2 2 '((0 0 0) (1 1 0)))) (kept (rd:analyze zero)))
      (assert (= 0 (rd-caller:find-code kept zero 1d-6 1d-12)))
      (assert (and (= 0 (rd:kept-count kept)) (= 2 (rd:left-out-count kept)) (<= 0d0 (rd:gap-bound kept) 1d-150)))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= 0 (rd-caller:reprove-code kept (csc-from-entries 2 2 '((0 0 3) (1 1 4))))))
      (assert (<= 5d0 (rd:gap-bound kept) 5.000001d0) () "The gap bound of diag(3, 4) with nothing kept is ~S." (rd:gap-bound kept))
      (assert (= 1 (rd-caller:certify-code kept 1d-9 1d-9))))
    ;; A bound that cannot be held. A row 2^-510 of another is kept and the other left out with a
    ;; multiplier of 2^510; against values where the kept row is of ordinary size the product's
    ;; square is beyond the format. The vectors stand, the bounds and any certificate do not, and
    ;; a reproof against values that can be bounded restores them.
    (let* ((n 8)
           (small (append (loop :for j :below n :collect (list 0 j (scale-float 1d0 -510)))
                          (loop :for j :below n :collect (list 1 j 1d0))))
           (large (append (loop :for j :below n :collect (list 0 j 1.9d0)) (loop :for j :below n :collect (list 1 j 1d0))))
           (matrix (csc-from-entries 2 n small)) (kept (rd:analyze matrix)))
      (assert (= 0 (rd-caller:find-code kept matrix 0d0 1d-12)))
      (assert (and (= 1 (rd:left-out-count kept)) (= 1 (aref (rd:left-out-rows kept) 0))
                   (= (- (scale-float 1d0 510)) (aref (rd:multiplier-values kept) 0)))
              () "The premise: row 1 is left out with a multiplier of 2^510; the vector holds ~S."
              (aref (rd:multiplier-values kept) 0))
      (let ((gap (rd:gap-bound kept)))
        (assert (= 2 (rd-caller:reprove-code kept (csc-from-entries 2 n large))))
        (assert (and (rd:kept-found? kept) (not (rd:certified? kept)) (= -1d0 (rd:gap-bound kept))
                     (= -1d0 (rd:residual-bound kept 0)) (= 1 (rd:multiplier-count kept 0))
                     (= 5 (rd-caller:certify-code kept 1d-9 1d-9))))
        (assert (= 0 (rd-caller:reprove-code kept matrix)))
        (assert (= gap (rd:gap-bound kept)))))
    ;; A budget too small is the analyses' typed error.
    (assert (= 98 (rd-caller:analysis-code matrix 3)))
    (assert (= 0 (rd-caller:analysis-code matrix 200)))
    ;; Steady state: proving the bounds again and certifying allocate nothing.
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
    (allocation-check "kept-rows-loop" (lambda () (assert (= 0 (rd-caller:kept-loop kept matrix 1d-9 1d-9 1)))))
    (format t "RANK-DEFICIENT-OUTCOMES-PASS~%")))

(defparameter *rank-deficient-public-control*
  "(declare probe (coalton-linalg/sparse/rank-deficient:KeptRows -> UFix))
   (define (probe kept)
     (+ (coalton-linalg/sparse/rank-deficient:kept-count kept) (coalton-linalg/sparse/rank-deficient:left-out-count kept)))")

(defparameter *rank-deficient-encapsulation-probes*
  '(("kept-rows-entries-access" .
     "(declare probe (coalton-linalg/sparse/rank-deficient:KeptRows -> arr:LispArray F64))
      (define (probe kept) (.entries kept))")
    ("kept-rows-unwrapping" .
     "(declare probe (coalton-linalg/sparse/rank-deficient:KeptRows -> UFix))
      (define (probe kept)
        (match kept ((coalton-linalg/sparse/rank-deficient:MkKeptRows held) (.rows held))))")
    ("find-status-as-reprove-status" .
     "(declare probe (coalton-linalg/sparse/rank-deficient:FindStatus -> coalton-linalg/sparse/rank-deficient:ReproveStatus))
      (define (probe status) status)")
    ("reprove-status-as-certificate-status" .
     "(declare probe (coalton-linalg/sparse/rank-deficient:ReproveStatus -> coalton-linalg/sparse/rank-deficient:CertificateStatus))
      (define (probe status) status)")))

(defun run-rank-deficient-encapsulation-tests ()
  "External Coalton callers reach the handle only through its public operations, and the status types stay apart."
  (multiple-value-bind (accepted diagnostics)
      (coalton-linalg/tests/encapsulation::compile-probe "rank-deficient-public-control" *rank-deficient-public-control*)
    (assert accepted () "The valid kept-rows external caller failed:~%~A" diagnostics))
  (dolist (probe *rank-deficient-encapsulation-probes*)
    (assert (not (coalton-linalg/tests/encapsulation::compile-probe (car probe) (cdr probe))) ()
            "External caller ~A compiled despite required rejection." (car probe)))
  (dolist (name '("MKKEPTROWS" "KEPTSTORAGE"))
    (multiple-value-bind (symbol visibility) (find-symbol name "COALTON-LINALG/SPARSE/RANK-DEFICIENT")
      (assert (and symbol (eq visibility ':internal)) ()
              "Storage constructor ~A must exist and remain internal; visibility is ~S." name visibility)))
  (format t "RANK-DEFICIENT-ENCAPSULATION-PASS rejected=~D public-control=1~%" (length *rank-deficient-encapsulation-probes*)))

;;; ---- The step ------------------------------------------------------------------------------

(defun rd-b-matrix (a kept m)
  "B exactly: the stored matrix with each row left out replaced by its kept combination of the kept rows, the
multipliers being the kept vector's stored doubles with the sign changed."
  (let* ((n (array-dimension a 1)) (b (rect m n (lambda (i j) (aref a i j)))) (rows (rd:multiplier-rows kept))
         (values (rd:multiplier-values kept)))
    (dotimes (p (rd:left-out-count kept) b)
      (let ((own (aref (rd:left-out-rows kept) p)) (start (rd:multiplier-start kept p)))
        (dotimes (j n)
          (setf (aref b own j)
                (loop :for e :below (rd:multiplier-count kept p)
                      :sum (* (- (rational (aref values (+ start e)))) (aref a (aref rows (+ start e)) j)))))))))

(defun rd-distance-squared (step exact)
  "The exact squared distance of the step written, a vector of doubles, from minus an exact solution."
  (loop :for x :across step :for e :in exact :sum (expt (+ (rational x) e) 2)))

(defun rd-check-step (label kept a m residual step)
  "Both bounds of the standing step hold exactly: the solve's against the exact step for B, and, where the stored
matrix has exactly the kept count for its rank, the whole one against its exact minimum-norm least squares
step. Returns the two exact squared distances; the second is nil where the rank is not the kept count."
  (let* ((b-values (coerce residual 'list)) (rhs (mapcar #'rational b-values))
         (for-b (rd-distance-squared step (exact-minimum-norm (rd-b-matrix a kept m) rhs)))
         (solve-bound (rd:step-solve-bound kept)) (whole-bound (rd:step-bound kept))
         (ranked (= (dr-rank a) (rd:kept-count kept)))
         (for-a (and ranked (rd-distance-squared step (exact-minimum-norm a rhs)))))
    (assert (and (rd:step-taken? kept) (<= 0d0 solve-bound whole-bound)) () "~A: the bounds are ~S and ~S." label
            solve-bound whole-bound)
    (assert (<= for-b (expt (rational solve-bound) 2)) ()
            "~A: the step is ~E from the exact step for B, and the solve's bound is ~E." label (sqrt (float for-b 1d0))
            solve-bound)
    (when ranked
      (assert (<= for-a (expt (rational whole-bound) 2)) ()
              "~A: the step is ~E from the exact step, and its bound is ~E." label (sqrt (float for-a 1d0)) whole-bound))
    (values for-b for-a)))

(defun rd-seeded-vector (state count)
  (let ((v (zeros count)))
    (dotimes (i count v) (setf (aref v i) (coerce (- (random 13 state) 6) 'double-float)))))

(defun run-step-exact-tests ()
  (let ((state (sb-ext:seed-random-state 20261124)) (cases 0) (stale 0) (stale-refused 0) (moved 0)
        (tally (make-list (1+ rd:correction-limit) :initial-element 0)) (widest 0d0))
    (loop :for (bases columns extras) :in '((2 4 1) (3 5 2) (4 6 3) (5 5 4) (6 9 5) (4 4 2) (7 10 6) (3 6 0))
          :do (dotimes (repeat 5)
                (multiple-value-bind (entries m) (rd-redundant-entries state bases columns extras)
                  (let* ((label (format nil "step, ~D+~D rows by ~D #~D" bases extras columns repeat))
                         (matrix (csc-from-entries m columns entries)) (a (dr-dense m columns entries))
                         (kept (rd:analyze matrix)) (out (zeros columns)))
                    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
                    (assert (= 1 (rd-caller:step-code kept (zeros m) out)) () "~A: a step was taken with no certificate." label)
                    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
                    ;; Two residuals: one anywhere, one in the range of the matrix.
                    (dolist (residual (list (rd-seeded-vector state m)
                                            (let ((x (rd-seeded-vector state columns)) (r (zeros m)))
                                              (dolist (entry entries r)
                                                (incf (aref r (first entry)) (* (third entry) (aref x (second entry))))))))
                      (incf cases)
                      (assert (= 0 (rd-caller:step-code kept residual out)) () "~A: the step gave ~D." label
                              (rd-caller:step-code kept residual out))
                      (rd-check-step label kept a m residual out)
                      (incf (nth (rd:step-corrections kept) tally))
                      ;; Fresh vectors: the whole bound is rounding beside the step and the residual.
                      (let ((size (+ (sqrt (loop :for x :across out :sum (* x x))) (sqrt (loop :for x :across residual :sum (* x x))))))
                        (setf widest (max widest (/ (rd:step-bound kept) (max size 1d-300))))
                        (assert (<= (rd:step-bound kept) (* 1d-9 size)) () "~A: the bound is ~E beside a size of ~E." label
                                (rd:step-bound kept) size)))
                    ;; The kept vectors out of date on purpose, and still under the floor: every row left
                    ;; out is multiplied by 1 + 2^-12, which keeps the rank and the row space and makes
                    ;; each kept combination wrong by that much. The step for B is then not the step.
                    (when (plusp (rd:left-out-count kept))
                      (let* ((left (coerce (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)) 'list))
                             (later (mapcar (lambda (entry)
                                              (list (first entry) (second entry)
                                                    (if (member (first entry) left)
                                                        (* (third entry) (+ 1d0 (scale-float 1d0 -12)))
                                                        (third entry))))
                                            entries))
                             (other (csc-from-entries m columns later)) (b (dr-dense m columns later))
                             (anywhere (rd-seeded-vector state m))
                             (in-range (let ((x (rd-seeded-vector state columns)) (r (zeros m)))
                                         (dolist (entry later r)
                                           (incf (aref r (first entry)) (* (third entry) (aref x (second entry))))))))
                        (assert (= 0 (rd-caller:reprove-code kept other)))
                        (assert (= 1 (rd-caller:step-code kept anywhere out)))
                        (let ((floor (* 2d0 (rd:gap-bound kept))))
                          (cond
                            ((/= 0 (rd-caller:certify-code kept floor floor)) (incf stale-refused))
                            (t
                             (dolist (residual (list anywhere in-range))
                               (let ((code (rd-caller:step-code kept residual out)))
                                 (assert (member code '(0 5)) () "~A, out of date: the step gave ~D." label code)
                                 (cond
                                   ((= code 5)
                                    (incf stale-refused)
                                    (assert (>= (* 4.001d0 (rd:gap-bound kept)) (rd:proved-lower-bound kept)) ()
                                            "~A: bound not reached with a gap of ~E under a proved ~E." label
                                            (rd:gap-bound kept) (rd:proved-lower-bound kept)))
                                   (t
                                    (incf stale)
                                    (multiple-value-bind (for-b for-a)
                                        (rd-check-step (format nil "~A, out of date" label) kept b m residual out)
                                      (when (> for-a (* 4 (max for-b (expt (rational (rd:step-solve-bound kept)) 2))))
                                        (incf moved))))))))))))))))
    (assert (and (> stale 20) (> moved 10)) ()
            "The premise: out-of-date vectors gave ~D steps, of which ~D were further from the exact step than from B's." stale moved)
    (format t "STEP-EXACT-PASS cases=~D out-of-date=~D of-them-moved=~D out-of-date-refused=~D by-corrections=~S widest-fresh-bound=~,2E~%"
            cases stale moved stale-refused tally widest)))

(defun run-step-scale-tests ()
  "The matrix and the residual each at its own scale give the step at scale 1 times the residual's power over the
matrix's, bit for bit, and bounds that are the first ones times the same."
  (let* ((entries '((0 0 1) (0 1 2) (1 1 1) (1 2 3) (2 0 2) (2 2 1) (3 0 1) (3 1 3) (3 2 3) (4 0 3) (4 1 2) (4 2 4)))
         (m 5) (n 3) (residual (doubles 1d0 -2d0 3d0 5d0 -4d0)) (cases 0))
    (flet ((at (matrix-power residual-power)
             (let* ((matrix (csc-from-entries m n (mapcar (lambda (entry)
                                                            (list (first entry) (second entry)
                                                                  (scale-float (coerce (third entry) 'double-float) matrix-power)))
                                                          entries)))
                    (kept (rd:analyze matrix)) (floor (scale-float 1d-6 matrix-power)) (out (zeros n))
                    (r (map '(simple-array double-float (*)) (lambda (x) (scale-float x residual-power)) residual)))
               (assert (= 0 (rd-caller:find-code kept matrix floor 1d-12)))
               (assert (= 0 (rd-caller:certify-code kept floor floor)))
               (assert (= 0 (rd-caller:step-code kept r out)) () "At 2^~D and 2^~D the step gave ~D." matrix-power residual-power
                       (rd-caller:step-code kept r out))
               (list out (rd:step-solve-bound kept) (rd:step-bound kept) (rd:step-corrections kept)))))
      (let ((reference (at 0 0)))
        (dolist (powers '((-900 0) (900 0) (0 -900) (0 900) (500 -500) (-500 500) (-400 -400) (300 300) (-1000 -1000)))
          (let* ((scaled (at (first powers) (second powers))) (shift (- (second powers) (first powers))))
            (incf cases)
            (assert (same-bits-p (map '(simple-array double-float (*)) (lambda (x) (scale-float x shift)) (first reference))
                                 (first scaled))
                    () "At 2^~D and 2^~D the step is ~S, and scaled from scale 1 it is ~S." (first powers) (second powers)
                    (first scaled) (map 'vector (lambda (x) (scale-float x shift)) (first reference)))
            ;; A bound below the normal range is rounded upward, so it is compared from both sides.
            (flet ((scaled-bound-p (first-bound bound)
                     (let ((exact (* (rational first-bound) (expt 2 shift))))
                       (if (>= exact least-positive-normalized-double-float)
                           (= exact (rational bound))
                           (<= exact (rational bound) (+ exact (* 2 (rational least-positive-double-float))))))))
              (assert (and (scaled-bound-p (second reference) (second scaled)) (scaled-bound-p (third reference) (third scaled))
                           (= (fourth reference) (fourth scaled)))
                      () "At 2^~D and 2^~D the bounds are ~S and ~S; at scale 1 they are ~S and ~S." (first powers)
                      (second powers) (second scaled) (third scaled) (second reference) (third reference)))))))
    (format t "STEP-SCALE-PASS cases=~D~%" cases)))

(defun run-step-outcome-tests ()
  (let* ((entries '((0 0 1) (0 1 2) (1 1 1) (1 2 3) (2 0 2) (2 2 1) (3 0 1) (3 1 3) (3 2 3) (4 0 3) (4 1 2) (4 2 4)))
         (matrix (csc-from-entries 5 3 entries)) (a (dr-dense 5 3 entries)) (kept (rd:analyze matrix))
         (residual (doubles 1d0 -2d0 3d0 5d0 -4d0)) (out (doubles 7d0 7d0 7d0))
         (nan (not-a-number)) (infinity sb-ext:double-float-positive-infinity))
    (flet ((untouched () (assert (equalp out (doubles 7d0 7d0 7d0)) () "A refused step wrote ~S." out))
           (nothing () (assert (and (not (rd:step-taken? kept)) (= -1d0 (rd:step-bound kept)) (= -1d0 (rd:step-solve-bound kept))
                                    (= 0 (rd:step-corrections kept))))))
      (assert (= 8 rd:correction-limit))
      ;; No certificate: before a search, after one, after a refused certificate, after a reproof.
      (assert (= 1 (rd-caller:step-code kept residual out)))
      (nothing)
      (assert (and (= -1d0 (rd:solve-growth kept)) (= -1d0 (rd:solve-backward-bound kept)) (= -1d0 (rd:predicted-rate kept))))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
      (assert (= 1 (rd-caller:step-code kept residual out)))
      (assert (= 2 (rd-caller:certify-code kept 1d-9 1d3)))
      (assert (= 1 (rd-caller:step-code kept residual out)))
      (untouched)
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= -1d0 (rd:solve-growth kept)))
      ;; Sizes and residuals that cannot be used leave everything as it was.
      (assert (= 2 (rd-caller:step-code kept (zeros 4) out)))
      (assert (= 2 (rd-caller:step-code kept (zeros 6) out)))
      (assert (= 2 (rd-caller:step-code kept residual (zeros 2))))
      (assert (= 2 (rd-caller:step-code kept residual (zeros 4))))
      (dolist (bad (list nan infinity (- infinity)))
        (let ((r (copy-seq residual)))
          (setf (aref r 3) bad)
          (assert (= 3 (rd-caller:step-code kept r out)) () "A residual holding ~S was taken." bad)))
      (untouched)
      (nothing)
      ;; A step, and what it leaves readable.
      (assert (= 0 (rd-caller:step-code kept residual out)))
      (rd-check-step "the outcome matrix" kept a 5 residual out)
      (assert (and (plusp (rd:solve-growth kept)) (plusp (rd:solve-backward-bound kept)) (plusp (rd:predicted-rate kept))
                   (<= (rd:step-corrections kept) rd:correction-limit))
              () "Growth ~S, backward bound ~S, predicted rate ~S, ~D corrections." (rd:solve-growth kept)
              (rd:solve-backward-bound kept) (rd:predicted-rate kept) (rd:step-corrections kept))
      ;; A second step on the same certificate uses the factor the first left, and is the same step.
      (let ((first (copy-seq out)) (bound (rd:step-bound kept)) (growth (rd:solve-growth kept)))
        (fill out 7d0)
        (assert (= 0 (rd-caller:step-code kept residual out)))
        (assert (and (same-bits-p first out) (= bound (rd:step-bound kept)) (= growth (rd:solve-growth kept))))
        ;; A refused step withdraws the standing one and writes nothing.
        (fill out 7d0)
        (assert (= 2 (rd-caller:step-code kept (zeros 4) out)))
        (untouched)
        (nothing)
        (assert (= growth (rd:solve-growth kept)))
        ;; A residual of zeros: the step is zero.
        (assert (= 0 (rd-caller:step-code kept (zeros 5) out)))
        (assert (and (every #'zerop out) (<= 0d0 (rd:step-bound kept) 1d-150)) () "A residual of zeros gave ~S within ~S." out
                (rd:step-bound kept))
        ;; A step does not survive a certificate, a reproof or a search, and neither does its factor.
        (assert (= 0 (rd-caller:step-code kept residual out)))
        (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
        (nothing)
        (assert (= -1d0 (rd:solve-growth kept)))
        (assert (= 0 (rd-caller:step-code kept residual out)))
        (assert (and (same-bits-p first out) (rd:step-taken? kept)))
        (assert (= 0 (rd-caller:reprove-code kept matrix)))
        (nothing)
        (assert (and (= -1d0 (rd:solve-growth kept)) (= 1 (rd-caller:step-code kept residual out))))
        (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
        (assert (= 0 (rd-caller:step-code kept residual out)))
        (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
        (nothing)
        (assert (= -1d0 (rd:solve-backward-bound kept)))))
    ;; The sign and the size, by hand: two equal rows [1; 1] and the residual [1, 3]. The least
    ;; squares answer is 2, the mean, and the step is minus it; the kept row alone would give -1.
    (let* ((pair '((0 0 1) (1 0 1))) (matrix (csc-from-entries 2 1 pair)) (kept (rd:analyze matrix)) (out (zeros 1)))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= 0 (rd-caller:step-code kept (doubles 1d0 3d0) out)))
      (assert (<= (abs (+ (aref out 0) 2d0)) 1d-14) () "Two equal rows against [1, 3]: the step is ~S." (aref out 0))
      (rd-check-step "two equal rows" kept (dr-dense 2 1 pair) 2 (doubles 1d0 3d0) out))
    ;; No row left out, on more columns than rows: the minimum-norm step, with nothing to solve
    ;; for the rows left out.
    (let* ((wide '((0 0 1) (0 1 2) (0 3 3) (1 1 1) (1 2 4) (1 3 1))) (matrix (csc-from-entries 2 4 wide))
           (kept (rd:analyze matrix)) (out (zeros 4)) (residual (doubles 3d0 -5d0)))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
      (assert (and (= 0 (rd:left-out-count kept)) (= 0 (rd-caller:certify-code kept 1d-9 1d-9))))
      (assert (= 0 (rd-caller:step-code kept residual out)))
      (rd-check-step "no row left out" kept (dr-dense 2 4 wide) 2 residual out)
      (assert (<= (rd:step-solve-bound kept) (rd:step-bound kept) (* (rd:step-solve-bound kept) (+ 1d0 1d-12))) ()
              "With no row left out the two bounds differ: ~S and ~S." (rd:step-solve-bound kept) (rd:step-bound kept)))
    ;; Rows with no entry and rows of stored zeros beside real ones; one entry; a matrix of zeros,
    ;; whose step is zero.
    (let* ((entries '((0 0 1) (0 1 2) (2 0 0) (2 1 0) (3 0 2) (3 1 4) (4 1 5)))
           (matrix (csc-from-entries 5 2 entries)) (kept (rd:analyze matrix)) (out (zeros 2))
           (residual (doubles 1d0 2d0 3d0 4d0 5d0)))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= 0 (rd-caller:step-code kept residual out)))
      (rd-check-step "an empty row and a row of zeros" kept (dr-dense 5 2 entries) 5 residual out))
    (let* ((one (csc-from-entries 1 1 '((0 0 5)))) (kept (rd:analyze one)) (out (zeros 1)))
      (assert (= 0 (rd-caller:find-code kept one 1d-6 1d-12)))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= 0 (rd-caller:step-code kept (doubles 10d0) out)))
      (assert (<= (abs (+ (aref out 0) 2d0)) 1d-14) () "One entry of 5 against 10: the step is ~S." (aref out 0)))
    ;; Steps that are not numbers the format holds, on one entry: 1 over 3, 1 over 7, 2 over 3 and
    ;; the like. The step written is a rounded number and its residual can compute to exactly
    ;; zero; the bound must still be above the distance to the exact quotient, which is not zero.
    (let ((nonzero 0))
      (dolist (pair '((3 1) (7 1) (3 2) (11 5) (13 1) (9 7) (49 1) (6 5)))
        (let* ((entries (list (list 0 0 (first pair)))) (one (csc-from-entries 1 1 entries)) (kept (rd:analyze one))
               (out (zeros 1)) (residual (doubles (coerce (second pair) 'double-float))))
          (assert (= 0 (rd-caller:find-code kept one 1d-6 1d-12)))
          (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
          (assert (= 0 (rd-caller:step-code kept residual out)))
          (multiple-value-bind (for-b for-a) (rd-check-step (format nil "one entry, ~S" pair) kept (dr-dense 1 1 entries) 1 residual out)
            (declare (ignore for-b))
            (when (plusp for-a) (incf nonzero)))))
      (assert (>= nonzero 6) () "The premise: the quotients are not numbers the format holds; ~D of 8 steps were inexact." nonzero))
    (let* ((zero (csc-from-entries 2 2 '((0 0 0) (1 1 0)))) (kept (rd:analyze zero)) (out (doubles 7d0 7d0)))
      (assert (= 0 (rd-caller:find-code kept zero 1d-6 1d-12)))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= 0 (rd-caller:step-code kept (doubles 3d0 4d0) out)) () "A matrix of zeros: the step gave ~D."
              (rd-caller:step-code kept (doubles 3d0 4d0) out))
      (assert (and (every #'zerop out) (<= 0d0 (rd:step-bound kept) 1d-100)) () "A matrix of zeros gave the step ~S within ~S." out
              (rd:step-bound kept)))
    ;; Bound not reached. One kept row [1, 0] and a second that was equal to it when the vector was
    ;; found and is 1.3 times it now: the gap bound is 0.3 and the kept row is proved above 0.35,
    ;; so the certificate holds between them, and four times the gap is not below what is proved.
    (let* ((kept-matrix (csc-from-entries 2 2 '((0 0 1) (1 0 1)))) (kept (rd:analyze kept-matrix)) (out (doubles 7d0 7d0)))
      (assert (= 0 (rd-caller:find-code kept kept-matrix 1d-6 1d-12)))
      (assert (= 0 (rd-caller:reprove-code kept (csc-from-entries 2 2 '((0 0 1) (1 0 1.3d0))))))
      (assert (<= 0.3d0 (rd:gap-bound kept) 0.3001d0) () "The premise: a gap bound of 0.3; it is ~S." (rd:gap-bound kept))
      (assert (= 0 (rd-caller:certify-code kept 0.35d0 0.35d0)))
      (assert (= 5 (rd-caller:step-code kept (doubles 1d0 1d0) out)))
      (assert (and (equalp out (doubles 7d0 7d0)) (not (rd:step-taken? kept)) (= -1d0 (rd:step-bound kept))))
      ;; The same rows 1.01 apart: the step is taken, and its bound covers the exact step.
      (let ((near '((0 0 1) (1 0 1.01d0))))
        (assert (= 0 (rd-caller:reprove-code kept (csc-from-entries 2 2 near))))
        (assert (= 0 (rd-caller:certify-code kept 0.05d0 0.05d0)))
        (assert (= 0 (rd-caller:step-code kept (doubles 1d0 1d0) out)))
        (multiple-value-bind (for-b for-a) (rd-check-step "rows 1.01 apart" kept (dr-dense 2 2 near) 2 (doubles 1d0 1d0) out)
          (assert (> for-a (* 100 (max for-b (expt (rational (rd:step-solve-bound kept)) 2)))) ()
                  "The premise: the step for B is not the exact step; the squared distances are ~E and ~E."
                  (float for-b 1d0) (float for-a 1d0)))
        ;; A residual the kept combination maps to nothing, [1, -1]: the step for B is zero and the
        ;; exact step is not, and what covers the difference is the residual outside the range.
        (assert (= 0 (rd-caller:step-code kept (doubles 1d0 -1d0) out)))
        (assert (<= (abs (aref out 0)) 1d-15) () "The premise: the step for B against [1, -1] is zero; it is ~S." (aref out 0))
        (multiple-value-bind (for-b for-a) (rd-check-step "rows 1.01 apart, a residual outside the range" kept
                                                          (dr-dense 2 2 near) 2 (doubles 1d0 -1d0) out)
          (declare (ignore for-b))
          (assert (> for-a 1d-6) () "The premise: the exact step is not zero; its squared distance is ~E." (float for-a 1d0)))))
    ;; Two rows a part in a million apart, kept both: the solve is ill conditioned, and its bound
    ;; still covers the exact step.
    (let* ((close (list (list 0 0 1) (list 0 1 1) (list 1 0 1) (list 1 1 1.000001d0)))
           (matrix (csc-from-entries 2 2 close)) (kept (rd:analyze matrix))
           (out (zeros 2)) (residual (doubles 1d0 2d0)))
      (assert (= 0 (rd-caller:find-code kept matrix 0d0 1d-12)))
      (assert (and (= 2 (rd:kept-count kept)) (= 0 (rd-caller:certify-code kept 1d-300 1d-10))))
      (assert (= 0 (rd-caller:step-code kept residual out)))
      (multiple-value-bind (for-b for-a) (rd-check-step "rows a part in a million apart" kept (dr-dense 2 2 close) 2 residual out)
        (declare (ignore for-b))
        (assert (and (> (abs (aref out 0)) 1d5) (<= (rd:step-bound kept) (* 1d-3 (abs (aref out 0))))) ()
                "Rows a part in a million apart: the step is ~S within ~S, at a squared distance of ~E." out
                (rd:step-bound kept) (float for-a 1d0))))
    ;; A large multiplier: a kept row 2^-20 the size of the row left out, which is 2^20 times it
    ;; plus another. The system of the rows left out then loses digits to cancellation; the
    ;; bound is wider for it, and still a bound.
    (dolist (power '(10 16 20 21))
      (let* ((small (scale-float 1d0 (- power)))
             (uneven (list (list 0 0 small) (list 0 1 small) (list 1 1 1) (list 1 2 1) (list 2 0 1) (list 2 1 2) (list 2 2 1)))
             (matrix (csc-from-entries 3 3 uneven)) (kept (rd:analyze matrix)) (out (zeros 3)))
        (assert (= 0 (rd-caller:find-code kept matrix 1d-9 1d-12)))
        (assert (and (= 1 (rd:left-out-count kept)) (= 2 (aref (rd:left-out-rows kept) 0))
                     (= (- (scale-float 1d0 power)) (aref (rd:multiplier-values kept) 0)))
                () "The premise: row 2 is left out with a multiplier of 2^~D; rows ~S, vector ~S." power
                (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)) (subseq (rd:multiplier-values kept) 0 2))
        (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)) () "The premise: the uneven rows are certified at 2^~D; ~D." power
                (rd-caller:certify-code kept 1d-9 1d-9))
        (dolist (residual (list (doubles 1d0 2d0 -3d0) (doubles 0.3d0 -0.7d0 0.9d0) (doubles 1d0 0d0 0d0) (doubles 0d0 0d0 1d0)))
          (assert (= 0 (rd-caller:step-code kept residual out)))
          (rd-check-step (format nil "a multiplier of 2^~D" power) kept (dr-dense 3 3 uneven) 3 residual out))
        (format t "STEP-LARGE-MULTIPLIER 2^~D: solve bound ~,2E and whole bound ~,2E of the step's norm~%" power
                (/ (rd:step-solve-bound kept) (sqrt (loop :for x :across out :sum (* x x))))
                (/ (rd:step-bound kept) (sqrt (loop :for x :across out :sum (* x x)))))))
    ;; A step that cannot be held: a residual of 2^1000 against a matrix of 2^-1000.
    (let* ((tiny (csc-from-entries 2 2 (list (list 0 0 (scale-float 1d0 -1000)) (list 1 1 (scale-float 3d0 -1000)))))
           (kept (rd:analyze tiny)) (out (doubles 7d0 7d0)) (floor (scale-float 1d-6 -1000)))
      (assert (= 0 (rd-caller:find-code kept tiny floor 1d-12)))
      (assert (= 0 (rd-caller:certify-code kept floor floor)))
      (assert (= 6 (rd-caller:step-code kept (doubles (scale-float 1d0 1000) 1d0) out)))
      (assert (and (equalp out (doubles 7d0 7d0)) (not (rd:step-taken? kept))))
      (assert (= 0 (rd-caller:step-code kept (doubles (scale-float 1d0 20) 1d0) out)))
      (assert (<= (abs (+ 1d0 (scale-float (aref out 0) -1020))) 1d-14) () "2^20 over 2^-1000: the step is ~S." (aref out 0)))
    ;; Steady state: a reproof, a certificate and two steps allocate nothing.
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
    (allocation-check "kept-rows-step-loop"
                      (lambda () (assert (= 0 (rd-caller:step-loop kept matrix 1d-9 1d-9 residual out 1)))))
    (format t "STEP-OUTCOMES-PASS~%")))

(defun run-step-tests ()
  (run-step-exact-tests)
  (run-step-scale-tests)
  (run-step-outcome-tests)
  (format t "STEP-TESTS-PASS~%"))

;;; ---- The refresh of kept vectors -----------------------------------------------------------

(defun rd-vectors (kept m)
  "Every kept vector as a dense vector, in order."
  (loop :for p :below (rd:left-out-count kept) :collect (rd-vector kept p m)))

(defun rd-scaled-rows (entries rows factor)
  "ENTRIES with every entry of the ROWS given multiplied by FACTOR."
  (mapcar (lambda (entry)
            (list (first entry) (second entry) (if (member (first entry) rows) (* (third entry) factor) (third entry))))
          entries))

(defun run-refresh-exact-tests ()
  "Rows left out are multiplied by 1 + 2^-12, which keeps the rank and the row space and makes each kept combination
wrong by that much. A refresh brings every bound back under its share of the target, proved against the new values
exactly; the certificate that was refused is issued; and the step is the exact step again, to a bound of rounding
size."
  (let ((state (sb-ext:seed-random-state 20261125)) (cases 0) (replaced 0) (partial 0) (everywhere 0))
    (loop :for (bases columns extras) :in '((2 4 1) (3 5 2) (4 6 3) (5 5 4) (6 9 5) (4 4 2) (7 10 6))
          :do (dotimes (repeat 4)
                (multiple-value-bind (entries m) (rd-redundant-entries state bases columns extras)
                  (let* ((label (format nil "refresh, ~D+~D rows by ~D #~D" bases extras columns repeat))
                         (matrix (csc-from-entries m columns entries)) (kept (rd:analyze matrix)) (out (zeros columns)))
                    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
                    (let* ((d (rd:left-out-count kept))
                           (left (coerce (subseq (rd:left-out-rows kept) 0 d) 'list))
                           (kept-rows (rd-kept-rows kept m))
                           (factor (+ 1d0 (scale-float 1d0 -12)))
                           (later (rd-scaled-rows entries left factor))
                           (other (csc-from-entries m columns later)) (b (dr-dense m columns later))
                           (first-vectors (rd-vectors kept m)))
                      ;; With no certificate ever issued there is no shift to solve with: the bounds are
                      ;; proved, as a reproof would, and nothing is replaced.
                      (assert (= 3 (rd-caller:refresh-code kept other 1d-9)) () "~A: a refresh with no shift gave ~D." label
                              (rd-caller:refresh-code kept other 1d-9))
                      (assert (and (= 0 (rd:refreshed-count kept)) (> (rd:gap-bound kept) 1d-6)))
                      (rd-check-bounds (format nil "~A, no shift" label) kept other b m)
                      (assert (= 0 (rd-caller:reprove-code kept matrix)))
                      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
                      ;; The same values: nothing is out of date, and nothing changes.
                      (assert (= 0 (rd-caller:refresh-code kept matrix 1d-9)))
                      (assert (and (= 0 (rd:refreshed-count kept)) (every #'same-bits-p first-vectors (rd-vectors kept m))
                                   (not (rd:certified? kept))))
                      ;; The later values, with a certificate and a step standing for the first.
                      (incf cases)
                      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
                      (assert (= 0 (rd-caller:step-code kept (rd-seeded-vector state m) out)))
                      (assert (and (rd:certified? kept) (rd:step-taken? kept)))
                      (assert (= 0 (rd-caller:refresh-code kept other 1d-9)) () "~A: the refresh gave ~D." label
                              (rd-caller:refresh-code kept other 1d-9))
                      (assert (= d (rd:refreshed-count kept)) () "~A: ~D of ~D vectors were replaced." label
                              (rd:refreshed-count kept) d)
                      (incf replaced d)
                      (assert (and (equal kept-rows (rd-kept-rows kept m)) (not (rd:certified? kept)) (not (rd:step-taken? kept))
                                   (= -1d0 (rd:proved-lower-bound kept))))
                      (rd-check-bounds label kept other b m)
                      (assert (<= (rd:gap-bound kept) 1d-10) () "~A: the gap bound after the refresh is ~E." label
                              (rd:gap-bound kept))
                      ;; Each new vector is the old one with its multipliers scaled, to rounding.
                      (loop :for old :in first-vectors :for new :in (rd-vectors kept m) :for own :in left
                            :do (dotimes (i m)
                                  (let ((expected (if (= i own) 1d0 (* factor (aref old i)))))
                                    (assert (<= (abs (- (aref new i) expected)) (* 1d-9 (max 1d0 (abs expected)))) ()
                                            "~A: a refreshed vector holds ~S at row ~D, not ~S." label (aref new i) i expected))))
                      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
                      (let ((residual (rd-seeded-vector state m)))
                        (assert (= 0 (rd-caller:step-code kept residual out)))
                        (rd-check-step label kept b m residual out)
                        (let ((size (+ (sqrt (loop :for x :across out :sum (* x x))) (sqrt (loop :for x :across residual :sum (* x x))))))
                          (assert (<= (rd:step-bound kept) (* 1d-9 size)) () "~A: after the refresh the step's bound is ~E beside ~E."
                                  label (rd:step-bound kept) size)))
                      ;; Only the first row left out moved: only its vector is replaced.
                      (when (> d 1)
                        (let* ((moved (rd-scaled-rows entries (list (first left)) factor))
                               (one-moved (csc-from-entries m columns moved)) (c (dr-dense m columns moved)))
                          (assert (= 0 (rd-caller:refresh-code kept matrix 1d-9)))
                          (assert (= d (rd:refreshed-count kept)))
                          (let ((before (rd-vectors kept m)))
                            (assert (= 0 (rd-caller:refresh-code kept one-moved 1d-9)))
                            (assert (= 1 (rd:refreshed-count kept)) () "~A: one row moved and ~D vectors were replaced." label
                                    (rd:refreshed-count kept))
                            (incf partial)
                            (loop :for old :in (rest before) :for new :in (rest (rd-vectors kept m))
                                  :do (assert (same-bits-p old new) () "~A: a vector that was not out of date changed." label))
                            (rd-check-bounds (format nil "~A, one row moved" label) kept one-moved c m)
                            (assert (<= (rd:gap-bound kept) 1d-10)))))
                      ;; Every row multiplied by a factor of its own, 1.125 to 1.5: the kept rows change
                      ;; as well, the span and the rank do not, and each multiplier changes by the ratio
                      ;; of two factors.
                      (let* ((every (mapcar (lambda (entry)
                                              (list (first entry) (second entry)
                                                    (* (third entry) (+ 1d0 (/ (+ 1 (mod (first entry) 4)) 8d0)))))
                                            entries))
                             (all-moved (csc-from-entries m columns every)) (e (dr-dense m columns every))
                             (residual (rd-seeded-vector state m)))
                        (assert (= 0 (rd-caller:refresh-code kept all-moved 1d-9)))
                        (when (> d 0) (incf everywhere (rd:refreshed-count kept)))
                        (rd-check-bounds (format nil "~A, every row moved" label) kept all-moved e m)
                        (assert (<= (rd:gap-bound kept) 1d-10) () "~A, every row moved: the gap bound is ~E." label
                                (rd:gap-bound kept))
                        (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
                        (assert (= 0 (rd-caller:step-code kept residual out)))
                        (rd-check-step (format nil "~A, every row moved" label) kept e m residual out)
                        (let ((size (+ (sqrt (loop :for x :across out :sum (* x x))) (sqrt (loop :for x :across residual :sum (* x x))))))
                          (assert (<= (rd:step-bound kept) (* 1d-9 size)) () "~A, every row moved: the step's bound is ~E beside ~E."
                                  label (rd:step-bound kept) size))))))))
    (assert (and (> replaced 50) (> partial 10) (> everywhere 50)) ()
            "The premise: ~D vectors replaced, ~D matrices with one row moved, ~D replaced with every row moved." replaced partial
            everywhere)
    (format t "REFRESH-EXACT-PASS cases=~D replaced=~D one-row-moved=~D replaced-with-every-row-moved=~D~%" cases replaced partial
            everywhere)))

(defun run-refresh-outcome-tests ()
  (let* ((entries '((0 0 1) (0 1 2) (1 1 1) (1 2 3) (2 0 2) (2 2 1) (3 0 1) (3 1 3) (3 2 3) (4 0 3) (4 1 2) (4 2 4)))
         (matrix (csc-from-entries 5 3 entries)) (kept (rd:analyze matrix))
         (nan (not-a-number)) (infinity sb-ext:double-float-positive-infinity))
    ;; Before a search; targets that cannot be used, which change nothing; another pattern or shape.
    (assert (and (= 1 (rd-caller:refresh-code kept matrix 1d-9)) (= 0 (rd:refreshed-count kept))))
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
    (dolist (target (list 0d0 -1d-9 nan infinity))
      (assert (= 2 (rd-caller:refresh-code kept matrix target)) () "A target of ~S was taken." target)
      (assert (rd:certified? kept)))
    (assert (= 91 (rd-caller:refresh-code kept (csc-from-entries 5 3 (cons '(0 2 1) entries)) 1d-9)))
    (assert (= 92 (rd-caller:refresh-code kept (csc-from-entries 6 3 entries) 1d-9)))
    ;; After a refused search nothing is kept, there is nothing to refresh, and the count of an
    ;; earlier refresh is not read.
    (let ((left (coerce (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)) 'list)))
      (assert (= 0 (rd-caller:refresh-code kept (csc-from-entries 5 3 (rd-scaled-rows entries left 1.25d0)) 1d-9)))
      (assert (= 2 (rd:refreshed-count kept))))
    (assert (= 4 (rd-caller:find-code kept matrix 1d-6 2d0)))
    (assert (= 0 (rd:refreshed-count kept)))
    (assert (= 1 (rd-caller:refresh-code kept matrix 1d-9)))
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
    ;; A vector that changes its rows. Row 2 holds a stored zero in the second column and is row 0;
    ;; with a 1 there it is row 0 plus row 1, and its vector gains an entry; with the zero again
    ;; it loses it. The system of the rows left out is analyzed again each time.
    (flet ((corner (x) (list '(0 0 1) '(1 1 1) '(2 0 1) (list 2 1 x))))
      (let* ((plain (csc-from-entries 3 2 (corner 0))) (kept (rd:analyze plain)) (out (zeros 2)))
        (assert (= 0 (rd-caller:find-code kept plain 1d-6 1d-12)))
        (assert (and (= 1 (rd:left-out-count kept)) (= 2 (aref (rd:left-out-rows kept) 0)) (= 1 (rd:multiplier-count kept 0)))
                () "The premise: row 2 is left out with one entry; rows ~S, ~D entries."
                (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)) (rd:multiplier-count kept 0))
        (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
        (dolist (x '(1 0 3 0))
          (let ((entries (corner x)))
            (assert (= 0 (rd-caller:refresh-code kept (csc-from-entries 3 2 entries) 1d-9)))
            (assert (and (= 1 (rd:refreshed-count kept)) (= (if (zerop x) 1 2) (rd:multiplier-count kept 0))) ()
                    "With ~D in the corner: ~D replaced, ~D entries." x (rd:refreshed-count kept) (rd:multiplier-count kept 0))
            (rd-check-bounds (format nil "a corner of ~D" x) kept (csc-from-entries 3 2 entries) (dr-dense 3 2 entries) 3)
            (assert (<= (rd:gap-bound kept) 1d-12))
            (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
            (assert (= 0 (rd-caller:step-code kept (doubles 1d0 2d0 4d0) out)))
            (rd-check-step (format nil "a corner of ~D" x) kept (dr-dense 3 2 entries) 3 (doubles 1d0 2d0 4d0) out)))))
    ;; A vector that keeps its length and changes its row: row 2 is row 0, then row 1, while row 3
    ;; stays row 0. The system of the rows left out has as many entries and another pattern: the
    ;; two rows left out share a kept row, and then do not.
    (flet ((turned (a b) (list '(0 0 1) '(1 1 1) (list 2 0 a) (list 2 1 b) '(3 0 1))))
      (let* ((plain (csc-from-entries 4 2 (turned 1 0))) (kept (rd:analyze plain)) (out (zeros 2))
             (residual (doubles 1d0 2d0 4d0 -3d0)))
        (assert (= 0 (rd-caller:find-code kept plain 1d-6 1d-12)))
        (assert (and (equalp #(2 3) (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)))
                     (= 1 (rd:multiplier-count kept 0)) (= 1 (rd:multiplier-count kept 1))
                     (= 0 (aref (rd:multiplier-rows kept) 0)) (= 0 (aref (rd:multiplier-rows kept) 1)))
                () "The premise: rows 2 and 3 are left out, each as row 0; rows ~S at ~S."
                (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)) (subseq (rd:multiplier-rows kept) 0 2))
        (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
        (dolist (pair '((0 1) (1 0) (0 2)))
          (let ((entries (turned (first pair) (second pair))))
            (assert (= 0 (rd-caller:refresh-code kept (csc-from-entries 4 2 entries) 1d-9)))
            (assert (and (= 1 (rd:refreshed-count kept)) (= 1 (rd:multiplier-count kept 0))
                         (= (if (zerop (first pair)) 1 0) (aref (rd:multiplier-rows kept) 0)))
                    () "Row 2 as ~S: ~D replaced, ~D entries, at row ~D." pair (rd:refreshed-count kept)
                    (rd:multiplier-count kept 0) (aref (rd:multiplier-rows kept) 0))
            (rd-check-bounds (format nil "row 2 as ~S" pair) kept (csc-from-entries 4 2 entries) (dr-dense 4 2 entries) 4)
            (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
            (assert (= 0 (rd-caller:step-code kept residual out)))
            (rd-check-step (format nil "row 2 as ~S" pair) kept (dr-dense 4 2 entries) 4 residual out)
            (assert (<= (rd:step-bound kept) 1d-9) () "Row 2 as ~S: the step's bound is ~E." pair (rd:step-bound kept))))))
    ;; A row that no longer depends on the kept rows. Row 2 holds a stored zero in a column no kept
    ;; row reaches; with a half there it is half a unit from their span. The refresh leaves the
    ;; best combination, its bound is that distance, and the certificate refuses on the gap: the
    ;; rank has changed, and the caller searches again.
    (flet ((spur (y) (list '(0 0 1) '(1 1 1) '(2 0 1) (list 2 2 y))))
      (let* ((plain (csc-from-entries 3 3 (spur 0))) (kept (rd:analyze plain)))
        (assert (= 0 (rd-caller:find-code kept plain 1d-6 1d-12)))
        (assert (and (= 1 (rd:left-out-count kept)) (= 2 (aref (rd:left-out-rows kept) 0))))
        (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
        (let ((entries (spur 0.5d0)))
          (assert (= 0 (rd-caller:refresh-code kept (csc-from-entries 3 3 entries) 1d-9)))
          (rd-check-bounds "a row half a unit from the span" kept (csc-from-entries 3 3 entries) (dr-dense 3 3 entries) 3)
          (assert (<= 0.5d0 (rd:gap-bound kept) 0.500001d0) () "A row half a unit from the span: the gap bound is ~S."
                  (rd:gap-bound kept))
          (assert (= 1 (rd-caller:certify-code kept 1d-9 1d-9))))))
    ;; The same matrix at other scales: the same vectors bit for bit, and a gap bound times the scale.
    (flet ((at (power)
             (let* ((scaled (lambda (list) (mapcar (lambda (entry)
                                                     (list (first entry) (second entry)
                                                           (scale-float (coerce (third entry) 'double-float) power)))
                                                   list)))
                    (first-matrix (csc-from-entries 5 3 (funcall scaled entries))) (kept (rd:analyze first-matrix))
                    (floor (scale-float 1d-6 power)))
               (assert (= 0 (rd-caller:find-code kept first-matrix floor 1d-12)))
               (assert (= 0 (rd-caller:certify-code kept floor floor)))
               (let* ((left (coerce (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)) 'list))
                      (later (rd-scaled-rows (funcall scaled entries) left 1.25d0)))
                 (assert (= 0 (rd-caller:refresh-code kept (csc-from-entries 5 3 later) floor)))
                 (assert (= 2 (rd:refreshed-count kept)))
                 (list (copy-seq (rd:multiplier-rows kept)) (copy-seq (rd:multiplier-values kept)) (rd:gap-bound kept))))))
      (let ((reference (at 0)))
        (dolist (power '(-900 -300 300 900))
          (let ((scaled (at power)))
            (assert (and (equalp (first reference) (first scaled)) (same-bits-p (second reference) (second scaled))
                         (= (scale-float (third reference) power) (third scaled)))
                    () "At 2^~D the refreshed vectors are ~S with a gap bound of ~S; at scale 1 ~S and ~S." power (second scaled)
                    (third scaled) (second reference) (third reference))))))
    ;; Steady state: a refresh that replaces two vectors, a certificate and a step allocate nothing.
    (let* ((left (progn (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
                        (coerce (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)) 'list)))
           (other (csc-from-entries 5 3 (rd-scaled-rows entries left 1.25d0))) (out (zeros 3))
           (residual (doubles 1d0 -2d0 3d0 5d0 -4d0)))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (allocation-check "kept-rows-refresh-loop"
                        (lambda () (assert (= 0 (rd-caller:refresh-loop kept matrix other 1d-9 residual out 1))))))
    (format t "REFRESH-OUTCOMES-PASS~%")))

(defun run-refresh-tests ()
  (run-refresh-exact-tests)
  (run-refresh-outcome-tests)
  (format t "REFRESH-TESTS-PASS~%"))

(defun run-rank-deficient-tests ()
  (run-rank-deficient-exact-tests)
  (run-rank-deficient-reprove-tests)
  (run-rank-deficient-scale-tests)
  (run-singular-bounds-tests)
  (run-rank-deficient-outcome-tests)
  (run-step-tests)
  (run-refresh-tests)
  (run-rank-deficient-encapsulation-tests)
  (format t "RANK-DEFICIENT-TESTS-PASS~%"))
