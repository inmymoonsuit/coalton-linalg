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
    (let ((gap (rd:gap-bound kept)))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= 91 (rd-caller:reprove-code kept (csc-from-entries 5 3 (cons '(0 2 1) entries)))))
      (assert (= 92 (rd-caller:reprove-code kept (csc-from-entries 6 3 entries))))
      ;; The values held were not touched, and what was proved of them stands.
      (assert (and (<= 0d0 gap) (= gap (rd:gap-bound kept)) (rd:certified? kept)) ()
              "After a reproof's typed error the gap bound is ~S, it was ~S, and a certificate stands: ~S."
              (rd:gap-bound kept) gap (rd:certified? kept)))
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

(defun run-step-subnormal-tests ()
  "A step below the normal range is rounded entry by entry when it is written, by up to half the smallest
positive number each, and the two bounds are of the step as written: the rounding of the step's own entries
is in them. With five entries at ties the distance is 1.118 of that number, above a bound that held only the
rounding of its own product."
  (let ((cases 0) (worst 0))
    (dolist (order '(1 2 5 9))
      (dolist (extra '(-15 -14 -3))
        (let* ((big (scale-float 1d0 1000))
               (matrix (csc-from-entries order order (loop :for i :below order :collect (list i i big))))
               (kept (rd:analyze matrix)) (floor (* big 1d-6))
               ;; (1 + 2^extra) 2^-60 over 2^1000: 2^-1060 and a part below the subnormal grid.
               (value (* (+ 1d0 (scale-float 1d0 extra)) (scale-float 1d0 -60)))
               (r (make-array order :element-type 'double-float :initial-element value)) (out (zeros order)))
          (assert (= 0 (rd-caller:find-code kept matrix floor 1d-15)))
          (assert (and (= order (rd:kept-count kept)) (= 0 (rd-caller:certify-code kept floor floor))) ()
                  "The premise: every row of 2^1000 I is kept and certified; ~D kept." (rd:kept-count kept))
          (assert (= 0 (rd-caller:step-code kept r out)))
          (assert (every (lambda (x) (< 0 (abs x) least-positive-normalized-double-float)) out) ()
                  "The premise: the step is below the normal range and not zero; it is ~S." out)
          (let* ((exact (- (/ (rational value) (rational big))))
                 (distance-squared (loop :for i :below order :sum (expt (- (rational (aref out i)) exact) 2)))
                 (eta (rational least-positive-double-float)))
            (incf cases)
            (setf worst (max worst (/ distance-squared (expt eta 2))))
            (dolist (bound (list (rd:step-solve-bound kept) (rd:step-bound kept)))
              (assert (<= distance-squared (expt (rational bound) 2)) ()
                      "Order ~D, 2^~D: the step is ~,4F of the smallest positive number from the exact one, and a bound is ~D of it."
                      order extra (sqrt (float (/ distance-squared (expt eta 2)) 1d0)) (/ (rational bound) eta))
              ;; And not more than the rounding of the bound's own product and of the step's entries.
              (assert (<= (rational bound) (* (+ 1 (ceiling (1+ order) 2)) eta)) ()
                      "Order ~D: a bound of ~D of the smallest positive number." order (/ (rational bound) eta)))))))
    (format t "STEP-SUBNORMAL-PASS cases=~D worst-distance=~,4F~%" cases (sqrt (float worst 1d0)))))

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
      ;; It is known from the gap bound and the bound proved alone, and nothing was factored for it.
      (assert (= -1d0 (rd:solve-growth kept)) () "A step refused on the certificate's two numbers made a factor, of growth ~S."
              (rd:solve-growth kept))
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
  (run-step-subnormal-tests)
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
  (let ((state (sb-ext:seed-random-state 20261125)) (cases 0) (replaced 0) (partial 0) (everywhere 0) (iterations 0))
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
                      ;; Before a refresh the later values leave every vector out of date, and the bounds
                      ;; say by how much. These vectors keep their rows, so a refresh needs no certificate
                      ;; to have been issued: it is made on each vector's own rows.
                      (assert (= 0 (rd-caller:reprove-code kept other)))
                      (assert (> (rd:gap-bound kept) 1d-6))
                      (rd-check-bounds (format nil "~A, before a refresh" label) kept other b m)
                      (assert (= 0 (rd-caller:refresh-code kept other 1d-9)) () "~A: a refresh before any certificate gave ~D."
                              label (rd-caller:refresh-code kept other 1d-9))
                      (assert (and (= d (rd:refreshed-count kept)) (<= (rd:gap-bound kept) 1d-10)))
                      (rd-check-bounds (format nil "~A, refreshed before any certificate" label) kept other b m)
                      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
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
                        ;; The factor of the rows left out's system is the search's, of the vectors as
                        ;; they were. A change of a part in four thousand: it preconditions this step's
                        ;; solves well, and is not made again.
                        (assert (and (= 1 (rd:system-factor-count kept)) (plusp (rd:system-iterations kept)))
                                () "~A: ~D iterations and ~D factorizations after a refresh." label (rd:system-iterations kept)
                                (rd:system-factor-count kept))
                        (incf iterations (rd:system-iterations kept))
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
    (assert (plusp iterations) () "The premise: some step after a refresh was preconditioned; ~D iterations in all." iterations)
    (format t "REFRESH-EXACT-PASS cases=~D replaced=~D one-row-moved=~D replaced-with-every-row-moved=~D preconditioned-iterations=~D~%"
            cases replaced partial everywhere iterations)))

(defun rd-growing-entries (state bases columns extras wide)
  "A matrix with no zero entry: BASES seeded rows, then EXTRAS rows each a multiple of one of them, or with WIDE each
the sum of small multiples of two or three. The two forms have one pattern and one row space."
  (let ((rows '()))
    (dotimes (b bases)
      (push (let ((row (make-array columns)))
              (dotimes (j columns row) (setf (aref row j) (let ((v (- (random 9 state) 4))) (if (zerop v) 5 v)))))
            rows))
    (let ((base (coerce (reverse rows) 'vector)))
      (dotimes (e extras)
        (let ((row (make-array columns :initial-element 0)) (first-source (random bases state)))
          (dotimes (j columns) (setf (aref row j) (* 2 (aref (aref base first-source) j))))
          (when wide
            (dotimes (more (1+ (random 2 state)))
              (let ((source (mod (+ first-source 1 more) bases)))
                (dotimes (j columns) (incf (aref row j) (* 4 (aref (aref base source) j))))))
            ;; A row that came out with a zero is put back as the plain multiple, so that the pattern stays whole.
            (when (some #'zerop row)
              (dotimes (j columns) (setf (aref row j) (* 2 (aref (aref base first-source) j))))))
          (push row rows))))
    (let ((all (coerce (reverse rows) 'vector)))
      (loop :for i :from 0 :for row :across all
            :append (loop :for j :below columns :collect (list i j (aref row j)))))))

(defun run-refresh-wider-tests ()
  "A vector that must take rows it did not hold: rows left out that are multiples of one kept row become sums of
two or three. The rows a vector holds no longer serve, the kept rows' solve finds the others, and the bounds, the
certificate and the step are exact again. With no certificate issued the solve has no shift, and says so."
  (let ((cases 0) (widened 0) (blocked 0))
    (loop :for (bases columns extras) :in '((3 5 2) (4 6 3) (5 8 4) (4 7 2))
          :do (dotimes (repeat 5)
                (let* ((seed (+ 77 (* 13 cases)))
                       (narrow (rd-growing-entries (sb-ext:seed-random-state seed) bases columns extras nil))
                       (wide (rd-growing-entries (sb-ext:seed-random-state seed) bases columns extras t))
                       (m (+ bases extras)) (label (format nil "wider, ~D+~D rows by ~D #~D" bases extras columns repeat))
                       (matrix (csc-from-entries m columns narrow)) (later (csc-from-entries m columns wide))
                       (b (dr-dense m columns wide)) (kept (rd:analyze matrix)) (out (zeros columns))
                       (state (sb-ext:seed-random-state (+ 5 seed))))
                  (incf cases)
                  (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
                  (assert (= (dr-rank b) (rd:kept-count kept)) () "~A: the premise, that the later matrix has the rank kept." label)
                  (let* ((d (rd:left-out-count kept)) (before (rd:multiplier-start kept d)) (kept-rows (rd-kept-rows kept m))
                         (first-vectors (rd-vectors kept m)))
                    ;; No certificate has been issued: the vectors that need the kept rows' solve are left.
                    (let ((code (rd-caller:refresh-code kept later 1d-9)))
                      (assert (member code '(0 3)) () "~A: a refresh with no shift gave ~D." label code)
                      (rd-check-bounds (format nil "~A, no shift" label) kept later b m)
                      (when (= code 3)
                        (incf blocked)
                        (assert (and (< (rd:refreshed-count kept) d) (> (rd:gap-bound kept) 1d-6)) ()
                                "~A: no shift, yet ~D of ~D replaced and a gap of ~E." label (rd:refreshed-count kept) d
                                (rd:gap-bound kept))
                        (assert (= 1 (rd-caller:certify-code kept 1d-9 1d-9)))))
                    ;; With one: every vector is brought under its share, whatever rows it needs.
                    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
                    (assert (every #'same-bits-p first-vectors (rd-vectors kept m)))
                    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
                    (assert (= 0 (rd-caller:refresh-code kept later 1d-9)) () "~A: the refresh gave ~D." label
                            (rd-caller:refresh-code kept later 1d-9))
                    (assert (equal kept-rows (rd-kept-rows kept m)))
                    (rd-check-bounds label kept later b m)
                    (assert (<= (rd:gap-bound kept) 1d-10) () "~A: the gap bound after the refresh is ~E." label (rd:gap-bound kept))
                    (when (> (rd:multiplier-start kept d) before) (incf widened))
                    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
                    (let ((residual (rd-seeded-vector state m)))
                      (assert (= 0 (rd-caller:step-code kept residual out)))
                      (rd-check-step label kept b m residual out)
                      (let ((size (+ (sqrt (loop :for x :across out :sum (* x x))) (sqrt (loop :for x :across residual :sum (* x x))))))
                        (assert (<= (rd:step-bound kept) (* 1d-9 size)) () "~A: the step's bound is ~E beside ~E." label
                                (rd:step-bound kept) size)))))))
    (assert (and (> widened 10) (> blocked 10)) ()
            "The premise: ~D of ~D cases widened a vector, and ~D were blocked with no shift." widened cases blocked)
    (format t "REFRESH-WIDER-PASS cases=~D widened=~D blocked-without-a-shift=~D~%" cases widened blocked)))

(defun run-refresh-restart-tests ()
  "A factor too far out of date to precondition with. Twelve rows left out, each a multiple of a kept row, are
multiplied by twelve different powers of two: the system of the rows left out then has a dozen eigenvalues the old
factor does not gather, the preconditioned iteration does not finish in its eight steps, and the factor is made
again. The step is the exact step either way."
  (let ((restarted 0) (cases 0))
    (dotimes (repeat 4)
      (let* ((entries (rd-growing-entries (sb-ext:seed-random-state (+ 300 repeat)) 3 5 12 nil)) (m 15) (columns 5)
             (matrix (csc-from-entries m columns entries)) (kept (rd:analyze matrix)) (out (zeros columns))
             (state (sb-ext:seed-random-state (+ 900 repeat))) (label (format nil "restart #~D" repeat)))
        (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
        (assert (= 12 (rd:left-out-count kept)) () "~A: the premise, twelve rows left out; ~D are." label (rd:left-out-count kept))
        (assert (and (= 0 (rd-caller:certify-code kept 1d-9 1d-9)) (= 1 (rd:system-factor-count kept))))
        (let* ((left (coerce (subseq (rd:left-out-rows kept) 0 12) 'list))
               (later (mapcar (lambda (entry)
                                (let ((at (position (first entry) left)))
                                  (list (first entry) (second entry)
                                        (if at (* (third entry) (scale-float 1d0 (1+ at))) (third entry)))))
                              entries))
               (other (csc-from-entries m columns later)) (b (dr-dense m columns later)) (residual (rd-seeded-vector state m)))
          (incf cases)
          (assert (= 0 (rd-caller:refresh-code kept other 1d-9)))
          (assert (and (= 12 (rd:refreshed-count kept)) (= 1 (rd:system-factor-count kept))) ()
                  "~A: ~D replaced, ~D factorizations after the refresh." label (rd:refreshed-count kept)
                  (rd:system-factor-count kept))
          (rd-check-bounds label kept other b m)
          (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
          (assert (= 0 (rd-caller:step-code kept residual out)))
          (rd-check-step label kept b m residual out)
          (let ((size (+ (sqrt (loop :for x :across out :sum (* x x))) (sqrt (loop :for x :across residual :sum (* x x))))))
            (assert (<= (rd:step-bound kept) (* 1d-9 size)) () "~A: the step's bound is ~E beside ~E." label (rd:step-bound kept) size))
          (when (= 2 (rd:system-factor-count kept))
            (incf restarted)
            (assert (<= rd:system-iteration-limit (rd:system-iterations kept)) ()
                    "~A: the factor was made again after ~D iterations." label (rd:system-iterations kept))
            ;; The next step finds the factor current: no iteration, and no factorization.
            (assert (= 0 (rd-caller:step-code kept residual out)))
            (assert (and (= 0 (rd:system-iterations kept)) (= 2 (rd:system-factor-count kept))))))))
    (assert (= 8 rd:system-iteration-limit))
    (assert (>= restarted 3) () "The premise: the factor was made again in ~D of ~D cases." restarted cases)
    (format t "REFRESH-RESTART-PASS cases=~D factor-made-again=~D~%" cases restarted)))

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
    ;; The fraction of its share a replaced vector is corrected to: from 0, which is as far as
    ;; the corrections go, to 1; nothing else.
    (assert (= 0.0625d0 rd:default-refresh-fraction))
    (dolist (fraction (list -0.5d0 1.5d0 nan infinity))
      (assert (= 2 (rd-caller:refresh-to-code kept matrix 1d-9 fraction)) () "A fraction of ~S was taken." fraction)
      (assert (rd:certified? kept)))
    (dolist (fraction (list 0d0 0.0625d0 1d0))
      (assert (= 0 (rd-caller:refresh-to-code kept matrix 1d-9 fraction))))
    ;; A fraction of zero still takes a vector on its own rows where that is at rounding: the rows
    ;; left out a quarter larger, no certificate issued, both vectors replaced.
    (let* ((fresh (rd:analyze matrix))
           (left (progn (assert (= 0 (rd-caller:find-code fresh matrix 1d-6 1d-12)))
                        (coerce (subseq (rd:left-out-rows fresh) 0 (rd:left-out-count fresh)) 'list)))
           (later (rd-scaled-rows entries left 1.25d0)))
      (assert (= 0 (rd-caller:refresh-to-code fresh (csc-from-entries 5 3 later) 1d-9 0d0)))
      (assert (and (= 2 (rd:refreshed-count fresh)) (<= (rd:gap-bound fresh) 1d-13)) ()
              "Corrected as far as it goes: ~D replaced, a gap bound of ~S." (rd:refreshed-count fresh) (rd:gap-bound fresh))
      (rd-check-bounds "a fraction of zero" fresh (csc-from-entries 5 3 later) (dr-dense 5 3 later) 5))
    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
    ;; The step's whole bound in its terms: none negative, and the whole no more than their sum.
    (let ((out (zeros 3)) (residual (doubles 1d0 -2d0 3d0 5d0 -4d0)))
      (assert (and (= -1d0 (rd:step-norm-term kept)) (= -1d0 (rd:step-outside-term kept)) (= 0 (rd:system-iterations kept))
                   (= 1 (rd:system-factor-count kept))))
      (assert (= 0 (rd-caller:step-code kept residual out)))
      (let ((solve (rd:step-solve-bound kept)) (in-norm (rd:step-norm-term kept)) (outside (rd:step-outside-term kept))
            (whole (rd:step-bound kept)))
        (assert (and (plusp solve) (plusp in-norm) (plusp outside) (<= whole (* (+ solve in-norm outside) (+ 1d0 1d-12)))
                     (>= whole (* 0.99d0 (+ solve in-norm outside))))
                () "The whole bound ~S and its terms ~S, ~S and ~S." whole solve in-norm outside))
      (assert (= 0 (rd:system-iterations kept)))
      (assert (= 2 (rd-caller:step-code kept (zeros 4) out)))
      (assert (and (= -1d0 (rd:step-norm-term kept)) (= -1d0 (rd:step-outside-term kept)) (= 0 (rd:system-iterations kept)))))
    (let ((gap (rd:gap-bound kept)) (certified (rd:certified? kept)))
      (assert (= 91 (rd-caller:refresh-code kept (csc-from-entries 5 3 (cons '(0 2 1) entries)) 1d-9)))
      (assert (= 92 (rd-caller:refresh-code kept (csc-from-entries 6 3 entries) 1d-9)))
      (assert (and (<= 0d0 gap) (= gap (rd:gap-bound kept)) (eq certified (rd:certified? kept))) ()
              "After a refresh's typed error the gap bound is ~S and it was ~S." (rd:gap-bound kept) gap))
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
        ;; The vector needs a row it does not hold, and no certificate has been issued: there is
        ;; no shift to solve with, the bound is proved as a reproof would, and nothing is replaced.
        (assert (= 3 (rd-caller:refresh-code kept (csc-from-entries 3 2 (corner 1)) 1d-9)))
        (assert (and (= 0 (rd:refreshed-count kept)) (= 1 (rd:multiplier-count kept 0)) (<= 1d0 (rd:gap-bound kept) 1.000001d0))
                () "With no shift: ~D replaced, a gap bound of ~S." (rd:refreshed-count kept) (rd:gap-bound kept))
        (assert (= 0 (rd-caller:reprove-code kept plain)))
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

(defun run-refresh-error-tests ()
  "A refresh that has exchanged its vectors and then cannot analyze the system of the rows left out, for
want of entries: nothing proved of the old vectors may stand beside the new ones. Sixteen rows in a chain
are kept and sixteen copies left out, each with a vector of one entry; changed to (1, 1.25) every copy has
left the span and its refreshed vector holds all sixteen kept rows, so the system is full."
  (flet ((chain (second)
           (csc-from-entries 32 17 (append (loop :for i :below 16 :append (list (list i i 1d0) (list i (1+ i) 1d0)))
                                           (loop :for j :below 16
                                                 :append (list (list (+ 16 j) j 1d0) (list (+ 16 j) (1+ j) second)))))))
    (let* ((first (chain 1d0)) (second (chain 1.25d0)) (kept (rd:analyze-with-limit first 120)))
      (assert (= 0 (rd-caller:find-code kept first 1d-6 1d-12)))
      (assert (and (= 16 (rd:left-out-count kept)) (= 1 (rd:multiplier-count kept 0))))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-3)))
      (let ((code (rd-caller:refresh-code kept second 1d-9)))
        (assert (= 98 code) () "The premise: the refresh ends in the storage error; its outcome is ~D." code))
      (assert (= 16 (rd:multiplier-count kept 0)) ()
              "The premise: the vectors were exchanged before the error; the first holds ~D rows."
              (rd:multiplier-count kept 0))
      (assert (and (= -1d0 (rd:gap-bound kept)) (= -1d0 (rd:residual-bound kept 0)) (not (rd:certified? kept))) ()
              "After the error a gap bound of ~S and a residual bound of ~S stand beside the new vectors."
              (rd:gap-bound kept) (rd:residual-bound kept 0))
      (assert (= 5 (rd-caller:certify-code kept 1d-9 1d-3)))
      (assert (= 1 (rd-caller:support-code kept 1d-6)))
      (let ((out (zeros 17)))
        (assert (= 1 (rd-caller:step-code kept (make-array 32 :element-type 'double-float :initial-element 1d0) out))))
      ;; The handle is not lost: proving the bounds again, for the vectors now held, restores it.
      (assert (= 0 (rd-caller:reprove-code kept second)))
      (assert (< 0d0 (rd:gap-bound kept) 0.3d0) () "Proved again, the gap bound is ~S." (rd:gap-bound kept))))
  (format t "REFRESH-ERROR-PASS~%"))

(defun run-refresh-tests ()
  (run-refresh-exact-tests)
  (run-refresh-wider-tests)
  (run-refresh-restart-tests)
  (run-refresh-outcome-tests)
  (run-refresh-error-tests)
  (format t "REFRESH-TESTS-PASS~%"))

;;; ---- Supports, dependent sets, and the right null-space basis ------------------------------

(defun rd-exact-multipliers (a kept m)
  "For each row left out, in order: its exact multipliers against the kept rows, w with a_t = w J_S, as a list in
the kept rows' order. The kept rows are independent, so w = (J_S J_S^T)^-1 J_S a_t."
  (let* ((kept-rows (rd-kept-rows kept m)) (r (length kept-rows)) (n (array-dimension a 1)))
    (if (zerop r)
        (loop :repeat (rd:left-out-count kept) :collect '())
        (let* ((js (dr-rows a kept-rows)) (inverse (rational-inverse (rect-product js (rect-transpose js)))))
          (loop :for p :below (rd:left-out-count kept)
                :collect (let* ((own (aref (rd:left-out-rows kept) p))
                                (b (loop :for k :below r :collect (loop :for j :below n :sum (* (aref js k j) (aref a own j))))))
                           (loop :for k :below r :collect (loop :for l :below r :for x :in b :sum (* (aref inverse k l) x)))))))))

(defun rd-exact-sets (kept m exact)
  "The dependent sets from exact multipliers: rows joined wherever a multiplier is not zero. Returns the sets that
hold a row left out, each a sorted list of rows, sorted; a row in none of them is in no dependency."
  (let ((parent (make-array m)) (kept-rows (rd-kept-rows kept m)))
    (dotimes (i m) (setf (aref parent i) i))
    (labels ((root (i) (if (= (aref parent i) i) i (setf (aref parent i) (root (aref parent i))))))
      (loop :for p :from 0 :for w :in exact
            :do (loop :for k :in kept-rows :for x :in w
                      :unless (zerop x) :do (setf (aref parent (root k)) (root (aref (rd:left-out-rows kept) p)))))
      (let ((roots (remove-duplicates (loop :for p :below (rd:left-out-count kept) :collect (root (aref (rd:left-out-rows kept) p))))))
        (sort (loop :for top :in roots :collect (loop :for i :below m :when (= (root i) top) :collect i))
              #'< :key #'first)))))

(defun rd-labelled-sets (kept m)
  "The same from the labels the library gives: the rows of each label below the row count."
  (let ((tags (rd:dependent-sets kept)))
    (sort (loop :for tag :below (rd:dependent-set-count kept)
                :collect (loop :for i :below m :when (= (aref tags i) tag) :collect i))
          #'< :key #'first)))

(defun run-support-exact-tests ()
  "The sets the library labels are the sets exact multipliers give, the stored multipliers are within the bound it
states of the exact ones, and a row in no dependency is labelled so."
  (let ((state (sb-ext:seed-random-state 20261126)) (cases 0) (sets 0) (alone 0) (several 0))
    (loop :for (bases columns extras) :in '((2 4 1) (3 5 2) (4 6 3) (5 5 4) (6 9 5) (4 4 2) (7 10 6) (8 12 3) (3 6 0))
          :do (dotimes (repeat 5)
                (multiple-value-bind (entries m) (rd-redundant-entries state bases columns extras)
                  (let* ((label (format nil "supports, ~D+~D rows by ~D #~D" bases extras columns repeat))
                         (matrix (csc-from-entries m columns entries)) (a (dr-dense m columns entries)) (kept (rd:analyze matrix)))
                    (incf cases)
                    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
                    (assert (= 1 (rd-caller:support-code kept 1d-6)) () "~A: classed with no certificate." label)
                    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
                    (assert (= 0 (rd-caller:support-code kept 1d-6)) () "~A: the classing gave ~D, vector ~D entry ~D." label
                            (rd-caller:support-code kept 1d-6) (rd:undecided-vector kept) (rd:undecided-entry kept))
                    (assert (and (rd:supports-classed? kept) (= m (rd:undecided-vector kept))))
                    (let* ((exact (rd-exact-multipliers a kept m)) (kept-rows (rd-kept-rows kept m))
                           (expected (rd-exact-sets kept m exact)) (found (rd-labelled-sets kept m)))
                      ;; The bound on each vector's multipliers, exactly.
                      (loop :for p :from 0 :for w :in exact
                            :do (let* ((v (rd-vector kept p m)) (bound (rd:multiplier-error-bound kept p))
                                       (distance (loop :for k :in kept-rows :for x :in w :sum (expt (+ (rational (aref v k)) x) 2))))
                                  (assert (and (plusp bound) (<= distance (expt (rational bound) 2))) ()
                                          "~A, vector ~D: its multipliers are ~E from the exact ones, and the bound is ~E." label p
                                          (sqrt (float distance 1d0)) bound)))
                      ;; Every stored entry is a real multiplier here, and is classed in.
                      (dotimes (e (rd:multiplier-start kept (rd:left-out-count kept)))
                        (assert (= 1 (aref (rd:support-flags kept) e))))
                      (assert (equal expected found) () "~A: the sets are ~S, and exactly they are ~S." label found expected)
                      (incf sets (length found))
                      (when (> (length found) 1) (incf several))
                      ;; A row in no set carries the row count.
                      (dotimes (i m)
                        (let ((inside (some (lambda (set) (member i set)) expected)))
                          (unless inside (incf alone))
                          (assert (if inside (< (aref (rd:dependent-sets kept) i) (rd:dependent-set-count kept))
                                      (= m (aref (rd:dependent-sets kept) i)))
                                  () "~A: row ~D is labelled ~D." label i (aref (rd:dependent-sets kept) i)))))))))
    (assert (and (> sets 30) (> alone 20) (> several 3)) ()
            "The premise: ~D sets in all, ~D rows in no dependency, ~D matrices with more than one set." sets alone several)
    (format t "SUPPORT-EXACT-PASS cases=~D sets=~D rows-in-no-dependency=~D with-several-sets=~D~%" cases sets alone several)))

(defun run-support-outcome-tests ()
  (let* ((entries '((0 0 1) (0 1 2) (1 1 1) (1 2 3) (2 0 2) (2 2 1) (3 0 1) (3 1 3) (3 2 3) (4 0 3) (4 1 2) (4 2 4)))
         (matrix (csc-from-entries 5 3 entries)) (kept (rd:analyze matrix)) (nan (not-a-number)))
    ;; Before a search, before a certificate, and cuts that cannot be used.
    (assert (and (= 1 (rd-caller:support-code kept 1d-6)) (not (rd:supports-classed? kept)) (= 0 (rd:dependent-set-count kept))
                 (= -1d0 (rd:multiplier-error-bound kept 0))))
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
    (assert (= 1 (rd-caller:support-code kept 1d-6)))
    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
    (dolist (cut (list 0d0 1d0 -1d-6 2d0 nan))
      (assert (= 0 (rd-caller:support-code kept 1d-6)))
      (assert (= 2 (rd-caller:support-code kept cut)) () "A cut of ~S was taken." cut)
      (assert (and (not (rd:supports-classed? kept)) (= 0 (rd:dependent-set-count kept)) (= -1d0 (rd:multiplier-error-bound kept 0)))))
    ;; Classed supports do not survive a certificate, a reproof, a refresh or a search.
    (dolist (undo (list (lambda () (rd-caller:certify-code kept 1d-9 1d-9)) (lambda () (rd-caller:reprove-code kept matrix))
                        (lambda () (rd-caller:refresh-code kept matrix 1d-9)) (lambda () (rd-caller:find-code kept matrix 1d-6 1d-12))))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= 0 (rd-caller:support-code kept 1d-6)))
      (assert (and (rd:supports-classed? kept) (= 1 (rd:dependent-set-count kept)) (plusp (rd:multiplier-error-bound kept 0))
                   (= -1d0 (rd:multiplier-error-bound kept 2))))
      (assert (= 0 (funcall undo)))
      (assert (and (not (rd:supports-classed? kept)) (= 0 (rd:dependent-set-count kept)))))
    ;; A multiplier at the cut, above it and below it. Row 2 is row 0 plus a small multiple of row 1.
    (flet ((mixed (small) (list '(0 0 1) '(1 1 1) '(2 0 1) (list 2 1 small))))
      (dolist (trial (list (list 1d-3 0 '(0 0 0) 1) (list 1d-6 3 nil nil) (list 1d-9 0 '(0 3 0) 0)))
        (destructuring-bind (small code tags flag) trial
          (let* ((entries (mixed small)) (matrix (csc-from-entries 3 2 entries)) (kept (rd:analyze matrix)))
            (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-15)))
            (assert (and (= 1 (rd:left-out-count kept)) (= 2 (aref (rd:left-out-rows kept) 0)) (= 2 (rd:multiplier-count kept 0)))
                    () "The premise: row 2 is left out with two entries; rows ~S, ~D entries."
                    (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)) (rd:multiplier-count kept 0))
            (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
            (assert (= code (rd-caller:support-code kept 1d-6)) () "A multiplier of ~S against a cut of 1e-6 gave ~D." small
                    (rd-caller:support-code kept 1d-6))
            (cond
              ((= code 3)
               (assert (and (= 0 (rd:undecided-vector kept)) (= 1 (rd:undecided-entry kept)) (not (rd:supports-classed? kept)))
                       () "Undecided at vector ~D, entry ~D." (rd:undecided-vector kept) (rd:undecided-entry kept)))
              (t
               (assert (and (equalp (coerce tags 'vector) (rd:dependent-sets kept)) (= 1 (rd:dependent-set-count kept))
                            (= 1 (aref (rd:support-flags kept) 0)) (= flag (aref (rd:support-flags kept) 1)))
                       () "A multiplier of ~S: labels ~S, flags ~S." small (rd:dependent-sets kept)
                       (subseq (rd:support-flags kept) 0 2))))))))
    ;; The cut is against the vector's largest entry with its own 1 counted. Row 2 is a thousandth
    ;; of row 0 plus a hundred-millionth of row 1: beside its own 1 the second multiplier is below
    ;; a millionth, and out; beside the first multiplier alone it would be above, and in.
    (let* ((entries (list (list 0 0 1) (list 1 1 1) (list 2 0 1d-3) (list 2 1 1d-8)))
           (matrix (csc-from-entries 3 2 entries)) (kept (rd:analyze matrix)))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-15)))
      (assert (and (= 1 (rd:left-out-count kept)) (= 2 (aref (rd:left-out-rows kept) 0)) (= 2 (rd:multiplier-count kept 0))))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= 0 (rd-caller:support-code kept 1d-6)))
      (assert (and (equalp #(0 3 0) (rd:dependent-sets kept)) (= 1 (aref (rd:support-flags kept) 0))
                   (= 0 (aref (rd:support-flags kept) 1)))
              () "Small multipliers: labels ~S, flags ~S." (rd:dependent-sets kept) (subseq (rd:support-flags kept) 0 2)))
    ;; Kept rows a part in a million apart: the bound on the multipliers is the vector's residual
    ;; over the bound proved on the kept rows, here near 1e-8. The vector holds both kept rows and
    ;; both multipliers are near 1, so a cut below that bound classes them as one above it does:
    ;; there is no kept row left to class out. (This test once held a refusal here, for a kept row
    ;; the vector does not hold; there is none. run-support-held-tests has one that exists.)
    (let* ((close (list (list 0 0 1) (list 0 1 1) (list 1 0 1) (list 1 1 1.000001d0) (list 2 0 2) (list 2 1 2.000001d0)))
           (matrix (csc-from-entries 3 2 close)) (kept (rd:analyze matrix)))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-7 1d-12)))
      (assert (and (= 2 (rd:kept-count kept)) (= 0 (rd-caller:certify-code kept 1d-12 1d-10))) ()
              "The premise: two rows a part in a million apart are kept and certified; ~D kept." (rd:kept-count kept))
      (assert (= 0 (rd-caller:support-code kept 1d-6)))
      (assert (and (equalp #(0 0 0) (rd:dependent-sets kept)) (<= 1d-12 (rd:multiplier-error-bound kept 0) 1d-6)) ()
              "Rows a part in a million apart: labels ~S, a bound of ~S." (rd:dependent-sets kept) (rd:multiplier-error-bound kept 0))
      (assert (= 2 (rd:multiplier-count kept 0)) () "The premise: the vector holds both kept rows; it holds ~D."
              (rd:multiplier-count kept 0))
      (dolist (cut '(1d-12 1d-15))
        (let ((code (rd-caller:support-code kept cut)))
          (assert (and (= 0 code) (equalp #(0 0 0) (rd:dependent-sets kept)) (rd:supports-classed? kept)) ()
                  "A cut of ~S gave ~D and labels ~S." cut code (rd:dependent-sets kept)))))
    ;; A row of zeros and a row with no entry are each a dependent set of one; a row no vector
    ;; holds is in none.
    (let* ((entries '((0 0 1) (0 1 2) (2 0 0) (2 1 0) (3 0 2) (3 1 4) (4 1 5)))
           (matrix (csc-from-entries 5 2 entries)) (kept (rd:analyze matrix)))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
      (assert (= 0 (rd-caller:support-code kept 1d-6)))
      (assert (and (= 3 (rd:dependent-set-count kept)) (equalp #(2 0 1 2 5) (rd:dependent-sets kept))) ()
              "An empty row and a row of zeros: ~D sets, labels ~S." (rd:dependent-set-count kept) (rd:dependent-sets kept)))
    ;; The same matrix at other scales: the same labels and the same bound.
    (flet ((at (power)
             (let* ((matrix (csc-from-entries 5 3 (mapcar (lambda (entry)
                                                            (list (first entry) (second entry)
                                                                  (scale-float (coerce (third entry) 'double-float) power)))
                                                          entries)))
                    (kept (rd:analyze matrix)) (floor (scale-float 1d-6 power)))
               (assert (= 0 (rd-caller:find-code kept matrix floor 1d-12)))
               (assert (= 0 (rd-caller:certify-code kept floor floor)))
               (assert (= 0 (rd-caller:support-code kept 1d-6)))
               (list (copy-seq (rd:dependent-sets kept)) (rd:multiplier-error-bound kept 0) (rd:multiplier-error-bound kept 1)))))
      (let ((reference (at 0)))
        (dolist (power '(-900 -300 300 900))
          (assert (equalp reference (at power)) () "At 2^~D the labels and bounds are ~S; at scale 1 ~S." power (at power) reference))))
    ;; Steady state: a reproof, a certificate and a classing allocate nothing.
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
    (assert (= 0 (rd-caller:support-loop kept matrix 1d-9 1d-6 1)))
    (allocation-check "kept-rows-support-loop" (lambda () (assert (= 0 (rd-caller:support-loop kept matrix 1d-9 1d-6 1)))))
    (format t "SUPPORT-OUTCOMES-PASS~%")))

(defun run-null-basis-tests ()
  "The handle on the transpose is a basis of the right null space. For a matrix whose columns depend on one
another, the kept vectors of its transpose are one vector for each column left out, with 1 at that column and 0
at every other column left out, so they are independent; there are as many as the column count less the rank;
and the bound of each holds for the exact product of the matrix with it."
  (let ((state (sb-ext:seed-random-state 20261127)) (cases 0) (vectors 0))
    (loop :for (bases rows extras) :in '((2 4 1) (3 5 2) (5 7 4) (6 6 3) (4 9 2))
          :do (dotimes (repeat 5)
                ;; A matrix of ROWS rows whose columns are BASES seeded vectors and EXTRAS combinations of them.
                (multiple-value-bind (transposed n) (rd-redundant-entries state bases rows extras)
                  (let* ((entries (mapcar (lambda (entry) (list (second entry) (first entry) (third entry))) transposed))
                         (j (dr-dense rows n entries)) (rank (dr-rank j))
                         (transpose (csc-from-entries n rows transposed)) (kept (rd:analyze transpose)))
                    (incf cases)
                    (assert (= 0 (rd-caller:find-code kept transpose 1d-6 1d-12)))
                    (assert (and (= rank (rd:kept-count kept)) (= (- n rank) (rd:left-out-count kept))) ()
                            "~D columns of rank ~D: ~D basis vectors." n rank (rd:left-out-count kept))
                    (assert (= 0 (rd-caller:certify-code kept 1d-9 1d-9)))
                    (let ((left (coerce (subseq (rd:left-out-rows kept) 0 (rd:left-out-count kept)) 'list)))
                      (dotimes (p (rd:left-out-count kept))
                        (let* ((z (rd-vector kept p n)) (bound (rd:residual-bound kept p))
                               ;; J z, exactly, by the matrix itself.
                               (product (loop :for i :below rows
                                              :sum (expt (loop :for c :below n :sum (* (aref j i c) (rational (aref z c)))) 2))))
                          (incf vectors)
                          (assert (<= product (expt (rational bound) 2)) ()
                                  "A basis vector's product has norm ~E and its bound is ~E." (sqrt (float product 1d0)) bound)
                          (loop :for column :in left :for q :from 0
                                :do (assert (= (if (= q p) 1d0 0d0) (aref z column)) ()
                                            "Basis vector ~D holds ~S at the column of vector ~D." p (aref z column) q)))))))))
    (format t "NULL-BASIS-PASS cases=~D vectors=~D~%" cases vectors)))

(defun run-support-tests ()
  (run-support-exact-tests)
  (run-support-outcome-tests)
  (run-null-basis-tests)
  (format t "SUPPORT-TESTS-PASS~%"))

;;; ---- Cases from the independent review of 2026-11-23 -----------------------------------------

(defun run-certify-floor-tests ()
  "A floor that is a singular value of the kept rows exactly, where the first trial meets a pivot of zero,
and a floor so far above the rows that its product with the scale is not a number: both are rows not above
the floor, and neither is a floor below what the method resolves."
  (let* ((matrix (csc-from-entries 1 1 (list (list 0 0 5d0)))) (kept (rd:analyze matrix)))
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-15)))
    (assert (= 0 (rd-caller:certify-code kept 1d-9 3d0)))
    (let ((code (rd-caller:certify-code kept 1d-9 5d0)))
      (assert (= 2 code) () "The one row [5] asked at a floor of 5: outcome ~D." code))
    (assert (= 2 (rd-caller:certify-code kept 1d-9 6d0))))
  (let* ((matrix (csc-from-entries 2 2 (list (list 0 0 5d0)))) (fresh (rd:analyze matrix)))
    ;; The same on a handle with no bound remembered, and with a row of zeros left out beside it.
    (assert (= 0 (rd-caller:find-code fresh matrix 1d-6 1d-15)))
    (let ((code (rd-caller:certify-code fresh 1d-9 5d0)))
      (assert (= 2 code) () "[5] beside an empty row at a floor of 5: outcome ~D." code)))
  (dolist (order '(1 2 4))
    (let* ((matrix (csc-from-entries order order (loop :for i :below order :collect (list i i 1d0))))
           (kept (rd:analyze matrix)))
      (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-15)))
      (let ((code (rd-caller:certify-code kept 1d-9 1d0)))
        (assert (= 2 code) () "The identity of order ~D asked at a floor of 1: outcome ~D." order code))
      (assert (= 0 (rd-caller:certify-code kept 1d-9 0.5d0)))
      ;; No row is left out, the gap bound is zero, and zero is below every lower floor there is.
      (let ((code (rd-caller:certify-code kept least-positive-double-float 0.5d0)))
        (assert (and (= 0d0 (rd:gap-bound kept)) (= 0 code)) ()
                "No row left out and a lower floor of the smallest positive number: outcome ~D, gap bound ~S." code
                (rd:gap-bound kept)))))
  (let* ((matrix (csc-from-entries 1 1 (list (list 0 0 (scale-float 1d0 -600))))) (kept (rd:analyze matrix)))
    (assert (= 0 (rd-caller:find-code kept matrix (scale-float 1d0 -620) 1d-15)))
    (let ((code (rd-caller:certify-code kept 1d0 (scale-float 1d0 500))))
      (assert (= 2 code) () "A row of 2^-600 asked at a floor of 2^500: outcome ~D." code))
    ;; A floor below the row is certified as before.
    (assert (= 0 (rd-caller:certify-code kept (scale-float 1d0 -700) (scale-float 1d0 -650)))))
  ;; A bound proved below the normal range is not rounded up: 0.3 proved of [[1.5]], and then of the same
  ;; row at 2^-1021, where 0.3 times that power lies halfway between two numbers.
  (let* ((matrix (csc-from-entries 1 1 (list (list 0 0 1.5d0))))
         (tiny (csc-from-entries 1 1 (list (list 0 0 (scale-float 1.5d0 -1021))))) (kept (rd:analyze matrix)))
    (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-15)))
    (assert (and (= 0 (rd-caller:certify-code kept 1d-9 0.3d0)) (= 0.3d0 (rd:proved-lower-bound kept))) ()
            "The premise: 0.3 is proved of [[1.5]]; the bound is ~S." (rd:proved-lower-bound kept))
    (assert (= 0 (rd-caller:reprove-code kept tiny)))
    (assert (= 0 (rd-caller:certify-code kept (scale-float 1d0 -1060) (scale-float 0.25d0 -1021))))
    (let ((proved (rational (rd:proved-lower-bound kept))) (exact (* (rational 0.3d0) (expt 2 -1021))))
      (assert (/= exact (rational (float exact 1d0))) () "The premise: the exact bound is not a double.")
      (assert (<= (* (rational 0.25d0) (expt 2 -1021)) proved exact) ()
              "A bound of ~D of the smallest positive number, where the bound proved is ~D of it."
              (/ proved (rational least-positive-double-float)) (/ exact (rational least-positive-double-float)))))
  (format t "CERTIFY-FLOOR-PASS~%"))

(defun run-subnormal-singular-bounds-tests ()
  "The lower bound on the largest singular value stays a lower bound where it is below the normal range."
  (let ((eta least-positive-double-float) (cases 0))
    (dolist (case (list (list 3 1 '((0 0) (1 0) (2 0)) 3)
                        (list 3 2 '((0 0) (1 1) (2 0) (2 1)) 3)
                        (list 2 1 '((0 0) (1 0)) 2)
                        (list 5 1 '((0 0) (1 0) (2 0) (3 0) (4 0)) 5)))
      (destructuring-bind (m n places squared) case
        (dolist (multiple '(1 3 7))
          (dolist (steps '(0 30))
            (let* ((value (* multiple eta))
                   (matrix (csc-from-entries m n (mapcar (lambda (place) (list (first place) (second place) value)) places)))
                   (out (zeros 2)))
              (assert (= 0 (rd-caller:bounds-code matrix steps (zeros n) (zeros m) out)))
              (incf cases)
              ;; For one column the largest singular value is the column's length exactly; for the
              ;; other matrix it is at least the longest column's.
              (when squared
                (assert (<= (expt (rational (aref out 0)) 2) (* squared (expt (rational value) 2))) ()
                        "~D by ~D of ~D eta, ~D steps: a lower bound of ~D eta above the length, whose square is ~D eta squared."
                        m n multiple steps (/ (rational (aref out 0)) (rational eta)) (* squared multiple multiple)))
              (assert (>= (expt (rational (aref out 1)) 2) (* (or squared 2) (expt (rational value) 2))) ()
                      "~D by ~D of ~D eta: an upper bound of ~D eta." m n multiple (/ (rational (aref out 1)) (rational eta)))
              (assert (<= (aref out 0) (aref out 1))))))))
    (format t "SUBNORMAL-SINGULAR-BOUNDS-PASS cases=~D~%" cases)))

(defun run-support-held-tests ()
  "What cannot be classed out is a kept row the vector does not hold. A vector that holds every kept row,
and a matrix with no kept row, have none."
  ;; Every kept row held: [1; 1] searched, then proved again for [1; 1.01], so that the vector's bound,
  ;; 0.01 over the 0.125 proved of the kept row, is above the cut. Its one multiplier is 1, and in.
  (let* ((first (csc-from-entries 2 1 (list (list 0 0 1d0) (list 1 0 1d0))))
         (second (csc-from-entries 2 1 (list (list 0 0 1d0) (list 1 0 1.01d0))))
         (kept (rd:analyze first)))
    (assert (= 0 (rd-caller:find-code kept first 1d-6 1d-15)))
    (assert (and (= 1 (rd:kept-count kept)) (= 1 (rd:multiplier-count kept 0))))
    (assert (= 0 (rd-caller:reprove-code kept second)))
    (assert (= 0 (rd-caller:certify-code kept 0.05d0 0.05d0)))
    (let ((code (rd-caller:support-code kept 0.05d0)))
      (assert (and (= 0 code) (equalp #(0 0) (rd:dependent-sets kept)) (= 1 (rd:dependent-set-count kept))) ()
              "A vector that holds its one kept row: outcome ~D, labels ~S, undecided at ~D, ~D." code
              (rd:dependent-sets kept) (rd:undecided-vector kept) (rd:undecided-entry kept))
      (assert (>= (rd:multiplier-error-bound kept 0) 0.05d0) ()
              "The premise: the multipliers' bound, ~S, is not below the cut." (rd:multiplier-error-bound kept 0))))
  ;; No kept row at all: two rows of stored zeros, proved again for diag(3, 4) and certified at a
  ;; floor above them. Each row left out is a set of one.
  (let* ((zeros-matrix (csc-from-entries 2 2 (list (list 0 0 0d0) (list 1 1 0d0))))
         (diagonal (csc-from-entries 2 2 (list (list 0 0 3d0) (list 1 1 4d0))))
         (kept (rd:analyze zeros-matrix)))
    (assert (= 0 (rd-caller:find-code kept zeros-matrix 1d-6 1d-15)))
    (assert (and (= 0 (rd:kept-count kept)) (= 0 (rd-caller:reprove-code kept diagonal))))
    (assert (= 0 (rd-caller:certify-code kept 10d0 10d0)))
    (let ((code (rd-caller:support-code kept 1d-6)))
      (assert (and (= 0 code) (= 2 (rd:dependent-set-count kept)) (equalp #(0 1) (rd:dependent-sets kept))) ()
              "No kept row: outcome ~D, ~D sets, labels ~S." code (rd:dependent-set-count kept) (rd:dependent-sets kept))))
  ;; A kept row the vector does not hold, and a bound above the cut: that is the refusal, and the
  ;; entry named is the vector's count. Row 2 is row 0 at the search and holds row 0 alone; proved
  ;; again with 0.01 of row 1 in it, its exact multiplier there is 0.01 and its bound 0.08.
  (let* ((first (csc-from-entries 3 2 (list (list 0 0 1d0) (list 1 1 1d0) (list 2 0 1d0) (list 2 1 0d0))))
         (second (csc-from-entries 3 2 (list (list 0 0 1d0) (list 1 1 1d0) (list 2 0 1d0) (list 2 1 0.01d0))))
         (kept (rd:analyze first)))
    (assert (= 0 (rd-caller:find-code kept first 1d-6 1d-15)))
    (assert (and (= 2 (rd:kept-count kept)) (= 1 (rd:multiplier-count kept 0))) ()
            "The premise: two rows kept and a vector that holds one; ~D kept, ~D held." (rd:kept-count kept)
            (rd:multiplier-count kept 0))
    (assert (= 0 (rd-caller:reprove-code kept second)))
    (assert (= 0 (rd-caller:certify-code kept 0.05d0 0.05d0)))
    (assert (= 3 (rd-caller:support-code kept 0.05d0)))
    (assert (and (= 0 (rd:undecided-vector kept)) (= 1 (rd:undecided-entry kept)) (not (rd:supports-classed? kept))) ()
            "Undecided at vector ~D, entry ~D." (rd:undecided-vector kept) (rd:undecided-entry kept))
    ;; A cut above the bound classes it out, and row 1 is in no dependency for these kept rows.
    (assert (= 0 (rd-caller:support-code kept 0.5d0)))
    (assert (equalp #(0 3 0) (rd:dependent-sets kept)) () "Labels ~S." (rd:dependent-sets kept)))
  (format t "SUPPORT-HELD-PASS~%"))

(defun run-support-order-tests ()
  "The sets are those of the kept rows at hand. The same three rows, one circuit, in two orders: kept as
(1,0) and (0,1) the third row's multiplier at the second is 2^-8, below a cut of 2^-4, and that row is in no
dependency; kept as (1,0) and (1,2^-8) both multipliers are 2^8, and the three are one set. Each classing is
proved; they are not the same sets."
  (let ((c (scale-float 1d0 -8)))
    (flet ((labels-of (entries m n upper)
             (let* ((matrix (csc-from-entries m n entries)) (kept (rd:analyze matrix)))
               (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-15)))
               (assert (= 0 (rd-caller:certify-code kept 1d-9 upper)))
               (assert (= 0 (rd-caller:support-code kept 0.0625d0)))
               (copy-seq (rd:dependent-sets kept)))))
      (let ((one (labels-of (list (list 0 0 1d0) (list 1 1 1d0) (list 2 0 1d0) (list 2 1 c)) 3 2 1d-3))
            (other (labels-of (list (list 0 0 1d0) (list 1 0 1d0) (list 1 1 c) (list 2 1 1d0)) 3 2 1d-4)))
        (assert (and (equalp #(0 3 0) one) (equalp #(0 0 0) other)) ()
                "The same rows in two orders: labels ~S and ~S." one other))
      (let ((one (labels-of (list (list 0 0 1d0) (list 1 0 c)) 2 1 1d-3))
            (other (labels-of (list (list 0 0 c) (list 1 0 1d0)) 2 1 1d-4)))
        (assert (and (equalp #(2 0) one) (equalp #(0 0) other)) ()
                "A row and 2^-8 of it, in two orders: labels ~S and ~S." one other))))
  (format t "SUPPORT-ORDER-PASS~%"))

(defun run-support-distance-tests ()
  "The bound on a vector's multipliers where it decides something: vectors found for one matrix and proved
again for another, so that the stored multipliers are hundredths from the exact ones, the least squares
multipliers of the row against the kept rows. The exact distance is under the bound on every one."
  (let ((cases 0) (largest 0))
    (dolist (case (list (list 2 1 '((0 0 1) (1 0 1)) '((0 0 1) (1 0 1.01d0)) 0.05d0)
                        (list 3 2 '((0 0 1) (1 1 1) (2 0 1) (2 1 1)) '((0 0 1) (1 1 1) (2 0 1.02d0) (2 1 0.97d0)) 0.2d0)
                        (list 3 2 '((0 0 2) (1 1 3) (2 0 2) (2 1 -3)) '((0 0 2) (1 1 3) (2 0 2.05d0) (2 1 -3.1d0)) 0.5d0)
                        (list 4 2 '((0 0 1) (1 1 1) (2 0 1) (2 1 1) (3 0 1) (3 1 -1))
                              '((0 0 1) (1 1 1) (2 0 1.03d0) (2 1 0.98d0) (3 0 0.96d0) (3 1 -1.01d0)) 0.2d0)))
      (destructuring-bind (m n first second floor) case
        (let* ((matrix (csc-from-entries m n first)) (other (csc-from-entries m n second)) (kept (rd:analyze matrix))
               (a (dr-dense m n second)))
          (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
          (assert (= 0 (rd-caller:reprove-code kept other)))
          (assert (= 0 (rd-caller:certify-code kept floor floor)) () "~D by ~D: not certified at ~S; the gap bound is ~S." m n
                  floor (rd:gap-bound kept))
          (assert (= 0 (rd-caller:support-code kept 0.25d0)))
          (let ((kept-rows (rd-kept-rows kept m)))
            (loop :for p :from 0 :for w :in (rd-exact-multipliers a kept m)
                  :do (let* ((v (rd-vector kept p m)) (bound (rd:multiplier-error-bound kept p))
                             (distance (loop :for k :in kept-rows :for x :in w :sum (expt (+ (rational (aref v k)) x) 2))))
                        (incf cases)
                        (setf largest (max largest (/ distance (expt (rational bound) 2))))
                        (assert (>= distance (expt 1/1000 2)) ()
                                "The premise: the stored multipliers are not the exact ones; they are ~E apart."
                                (sqrt (float distance 1d0)))
                        (assert (<= distance (expt (rational bound) 2)) ()
                                "~D by ~D, vector ~D: its multipliers are ~E from the exact ones, and the bound is ~E." m n p
                                (sqrt (float distance 1d0)) bound)))))))
    (format t "SUPPORT-DISTANCE-PASS vectors=~D largest-distance-over-bound=~,3F~%" cases (sqrt (float largest 1d0)))))

(defun rd-rational-root (x digits)
  "The square root of a non-negative rational, from below, to within 10^-DIGITS over its denominator."
  (/ (isqrt (* (numerator x) (denominator x) (expt 10 (* 2 digits)))) (* (denominator x) (expt 10 digits))))

(defun rd-truncated-solution (a m r)
  "For a matrix of two columns whose two singular values differ: J_1^+ r, J_1 the matrix truncated to its larger
singular value, as (G - l2 I) J^T r / (l1 (l1 - l2)) with G = J^T J and l1 > l2 its eigenvalues. The root in
them is taken to eighty digits, so this is the exact vector to far below a rounding of a double. Second value: l2."
  (let* ((g11 (loop :for i :below m :sum (expt (aref a i 0) 2))) (g22 (loop :for i :below m :sum (expt (aref a i 1) 2)))
         (g12 (loop :for i :below m :sum (* (aref a i 0) (aref a i 1))))
         (trace (+ g11 g22)) (root (rd-rational-root (- (expt trace 2) (* 4 (- (* g11 g22) (expt g12 2)))) 80))
         (l1 (/ (+ trace root) 2)) (l2 (/ (- trace root) 2))
         (t1 (loop :for i :below m :for x :in r :sum (* (aref a i 0) x)))
         (t2 (loop :for i :below m :for x :in r :sum (* (aref a i 1) x)))
         (scale (* l1 (- l1 l2))))
    (values (list (/ (+ (* (- g11 l2) t1) (* g12 t2)) scale) (/ (+ (* g12 t1) (* (- g22 l2) t2)) scale)) l2)))

(defun run-step-truncation-tests ()
  "The whole bound where something is truncated. One row is kept and the others were equal to it when the
vectors were found; proved again they are a hundredth from it and a thousandth out of its line, so the stored
matrix has rank two and the step is that of its truncation to rank one, which is not its own step and not the
step for B. The step written is within step-bound of the truncated matrix's step."
  (let ((cases 0) (largest 0))
    (dolist (case (list (list 2 '((0 0 1) (1 0 1) (1 1 0)) '((0 0 1) (1 0 1.01d0) (1 1 0.001d0)))
                        (list 2 '((0 0 1) (1 0 1) (1 1 0)) '((0 0 1) (1 0 0.98d0) (1 1 0.004d0)))
                        (list 3 '((0 0 1) (1 0 1) (1 1 0) (2 0 1) (2 1 0))
                              '((0 0 1) (1 0 1.01d0) (1 1 0.001d0) (2 0 0.99d0) (2 1 -0.002d0)))
                        (list 3 '((0 0 2) (1 0 2) (1 1 0) (2 0 2) (2 1 0))
                              '((0 0 2) (1 0 2.03d0) (1 1 0.01d0) (2 0 1.96d0) (2 1 0.005d0)))))
      (destructuring-bind (m first second) case
        (let* ((matrix (csc-from-entries m 2 first)) (other (csc-from-entries m 2 second)) (kept (rd:analyze matrix))
               (a (dr-dense m 2 second)))
          (assert (= 0 (rd-caller:find-code kept matrix 1d-6 1d-12)))
          (assert (= 1 (rd:kept-count kept)))
          (assert (= 0 (rd-caller:reprove-code kept other)))
          (assert (= 0 (rd-caller:certify-code kept 0.5d0 0.5d0)) () "Not certified; the gap bound is ~S." (rd:gap-bound kept))
          (dolist (residual (list (make-array m :element-type 'double-float :initial-element 1d0)
                                  (coerce (loop :for i :below m :collect (- (* 0.75d0 (1+ i)) 1.25d0))
                                          '(simple-array double-float (*)))))
            (let ((out (zeros 2)) (rhs (map 'list #'rational residual)))
              (assert (= 0 (rd-caller:step-code kept residual out)))
              ;; The rate predicted is of the iteration matrix's norm, the root of two times 2^-10 where
              ;; the factorization's own error adds nothing, and not of its two blocks, 2^-10 each.
              (assert (>= (rd:predicted-rate kept) (* 1.41421356d0 (expt 2d0 -10))) ()
                      "A predicted rate of ~S, under the root of two times 2^-10." (rd:predicted-rate kept))
              (multiple-value-bind (truncated l2) (rd-truncated-solution a m rhs)
                (let* ((distance (rd-distance-squared out truncated)) (bound (rational (rd:step-bound kept)))
                       (own (rd-distance-squared out (exact-minimum-norm a rhs))))
                  (incf cases)
                  (setf largest (max largest (/ distance (expt bound 2))))
                  (assert (and (= 2 (dr-rank a)) (plusp l2)) () "The premise: the stored matrix has rank two.")
                  ;; The matrix's own step is another vector: what is bounded is the truncation's.
                  (assert (> own 1/100) ()
                          "The premise: the step is not the stored matrix's own; it is ~E from that."
                          (sqrt (float own 1d0)))
                  (assert (<= distance (+ (expt bound 2) (expt 10 -60))) ()
                          "~D rows: the step is ~E from the truncated matrix's step, and step-bound is ~E." m
                          (sqrt (float distance 1d0)) (float bound 1d0)))))))))
    (format t "STEP-TRUNCATION-PASS cases=~D largest-distance-over-bound=~,3F~%" cases (sqrt (float largest 1d0)))))

(defun run-rank-deficient-tests ()
  (run-rank-deficient-exact-tests)
  (run-rank-deficient-reprove-tests)
  (run-rank-deficient-scale-tests)
  (run-singular-bounds-tests)
  (run-rank-deficient-outcome-tests)
  (run-step-tests)
  (run-refresh-tests)
  (run-support-tests)
  (run-certify-floor-tests)
  (run-subnormal-singular-bounds-tests)
  (run-support-held-tests)
  (run-support-order-tests)
  (run-support-distance-tests)
  (run-step-truncation-tests)
  (run-rank-deficient-encapsulation-tests)
  (format t "RANK-DEFICIENT-TESTS-PASS~%"))
