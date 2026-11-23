(in-package #:coalton-linalg/tests)

;;; Least squares on an operator. Four kinds of evidence, none of which is the library's own
;;; word: the residual and the gradient of every returned vector are formed here in exact
;;; rational arithmetic; the returned vector is held against the exact minimizer, of minimum
;;; norm where there is no damping; what the method is said to do from one iterate to the next
;;; is measured exactly; and the whole iteration is recomputed here and compared bit for bit.

(eval-when (:compile-toplevel :load-toplevel :execute)
  (sb-ext:add-package-local-nickname "LSQR" "COALTON-LINALG/LSQR")
  (sb-ext:add-package-local-nickname "LS-CALLER" "COALTON-LINALG/TESTS/LSQR-CALLERS"))

(defun ls-binade (magnitude)
  "The power of two p with p <= magnitude < 2 p, and 1 for zero: what the library divides by."
  (if (zerop magnitude) 1d0 (scale-float 1d0 (1- (nth-value 1 (decode-float magnitude))))))

(defun ls-norm (vector sum)
  "The 2-norm from a sum of squares in index order, measured again with the entries lifted where the sum is tiny."
  (declare (type double-float sum))
  (if (< sum least-positive-normalized-double-float)
      (let ((total 0d0) (lift (scale-float 1d0 600)))
        (declare (type double-float total))
        (loop :for value :across vector :do (let ((v (* value lift))) (incf total (* v v))))
        (/ (sqrt total) lift))
      (if (and (= sum sum) (< sum sb-ext:double-float-positive-infinity))
          (sqrt sum)
          (let ((total 0d0) (lift (scale-float 1d0 600)))
            (declare (type double-float total))
            (loop :for value :across vector :do (let ((v (/ value lift))) (incf total (* v v))))
            (* (sqrt total) lift)))))

(defun ls-twice-scaled (value first second)
  "VALUE times two powers of two, by their product where that is a normal number and one after the other where not."
  (let ((together (* first second)))
    (if (and (>= together least-positive-normalized-double-float) (< together sb-ext:double-float-positive-infinity))
        (* value together)
        (* (* value first) second))))

(defun ls-hypotenuse (a b)
  (let* ((x (abs a)) (y (abs b)) (large (if (> x y) x y)) (small (if (> x y) y x)))
    (if (= large 0d0)
        0d0
        (let ((q (/ small large))) (* large (sqrt (+ 1d0 (* q q))))))))

(defun ls-scaled-copy (vector factor)
  (map '(simple-array double-float (*)) (lambda (value) (* value factor)) vector))

(defun reference-lsqr (matrix m n b damping residual-tolerance gradient-tolerance limit)
  "The iteration of docs/lsqr.md recomputed here with the library's two products. Returns the outcome's code, the
vector, the steps, the confirmations, the two histories, the two norms last formed, the operator norm estimate and
every iterate in the caller's units, oldest first."
  (let* ((rhs (zeros m)) (u (zeros m)) (p (zeros m)) (v (zeros n)) (w (zeros n)) (q (zeros n)) (x (zeros n))
         (residuals '()) (gradients '()) (iterates '()) (sum 0d0)
         (b-up (ls-binade (loop :for value :across b :maximize (abs value)))) (b-down (/ 1d0 b-up)))
    (declare (type double-float sum))
    (dotimes (i m)
      (let ((value (* (aref b i) b-down)))
        (setf (aref rhs i) (aref b i) (aref u i) value)
        (incf sum (* value value))))
    (let ((beta-one (ls-norm u sum)))
      (when (= beta-one 0d0)
        (return-from reference-lsqr (values 0 x 0 0 (list 0d0) (list 0d0) 0d0 0d0 -1d0 (list (copy-seq x)))))
      (dotimes (i m) (setf (aref u i) (/ (aref u i) beta-one)))
      (csc:transpose-matvec! matrix u q)
      (let* ((a-up (ls-binade (loop :for value :across q :maximize (abs value)))) (a-down (/ 1d0 a-up)))
        (setf sum 0d0)
        (dotimes (j n)
          (let ((value (* (aref q j) a-down)))
            (setf (aref v j) value)
            (incf sum (* value value))))
        (let* ((alpha-one (ls-norm v sum)) (damp (* damping a-down)) (back (* b-up a-down))
               (alpha alpha-one) (rho-bar alpha-one) (phi-bar beta-one) (damped 0d0) (frobenius 0d0)
               (estimate (* alpha-one a-up))
               (residual-trigger residual-tolerance) (gradient-trigger gradient-tolerance)
               (ended (= alpha-one 0d0)) (confirmed nil) (k 0) (confirmations 0)
               (formed-r 0d0) (formed-g 0d0) (divided-r 0d0) (divided-g 0d0))
          (declare (type double-float alpha rho-bar phi-bar damped frobenius))
          (unless (= alpha-one 0d0)
            (dotimes (j n)
              (let ((value (/ (aref v j) alpha-one)))
                (setf (aref v j) value (aref w j) value))))
          (push (* beta-one b-up) residuals)
          (push (ls-twice-scaled (* alpha-one beta-one) a-up b-up) gradients)
          (push (ls-scaled-copy x back) iterates)
          (labels ((confirm ()
                     (csc:matvec! matrix x p)
                     (let ((largest 0d0))
                       (dotimes (i m)
                         (let ((value (- (aref rhs i) (* (aref p i) back))))
                           (setf (aref p i) value)
                           (when (> (abs value) largest) (setf largest (abs value)))))
                       (let* ((r-up (ls-binade largest)) (r-down (/ 1d0 r-up)) (image-top 0d0) (pulled-top 0d0))
                         (setf sum 0d0)
                         (dotimes (i m)
                           (let ((value (* (aref p i) r-down)))
                             (setf (aref p i) value)
                             (incf sum (* value value))))
                         (setf divided-r (ls-norm p sum))
                         (csc:transpose-matvec! matrix p q)
                         (dotimes (j n)
                           (let ((image (* (aref q j) a-down)) (pulled (abs (* damp (* damp (aref x j))))))
                             (setf (aref q j) image)
                             (when (> (abs image) image-top) (setf image-top (abs image)))
                             (when (> pulled pulled-top) (setf pulled-top pulled))))
                         (let ((image-led (>= (ls-twice-scaled image-top r-up b-down) pulled-top)))
                           (setf sum 0d0)
                           (dotimes (j n)
                             (let* ((pulled (* damp (* damp (aref x j))))
                                    (value (cond ((= pulled 0d0)
                                                  (if image-led (aref q j) (ls-twice-scaled (aref q j) r-up b-down)))
                                                 (image-led (- (aref q j) (ls-twice-scaled pulled b-up r-down)))
                                                 (t (- (ls-twice-scaled (aref q j) r-up b-down) pulled)))))
                               (setf (aref q j) value)
                               (incf sum (* value value))))
                           (setf divided-g (ls-norm q sum))
                           (setf formed-r (* divided-r r-up)
                                 formed-g (ls-twice-scaled divided-g a-up (if image-led r-up b-up))))))
                     (incf confirmations))
                   (within (reported divided tolerance)
                     (and (<= reported tolerance) (or (> tolerance 0d0) (= divided 0d0))))
                   (tightened (tolerance running formed)
                     (let ((trigger (* tolerance (/ running formed)))) (if (> trigger 0d0) trigger -1d0)))
                   (done (code)
                     (return-from reference-lsqr
                       (values code (ls-scaled-copy x back) k confirmations (nreverse residuals) (nreverse gradients)
                               formed-r formed-g estimate (nreverse iterates)))))
            (loop
              (let* ((running-r (first residuals)) (running-g (first gradients))
                     (hit-r (<= running-r residual-trigger)) (hit-g (<= running-g gradient-trigger)))
                (when (or hit-r hit-g)
                  (confirm)
                  (setf confirmed t)
                  (when (within formed-r divided-r residual-tolerance) (done 0))
                  (when (within formed-g divided-g gradient-tolerance) (done 1))
                  (when hit-r (setf residual-trigger (tightened residual-tolerance running-r formed-r)))
                  (when hit-g (setf gradient-trigger (tightened gradient-tolerance running-g formed-g)))))
              (when (or ended (= k limit)) (return))
              (setf confirmed nil)
              (csc:matvec! matrix v p)
              (let ((alpha-now alpha) (beta 0d0) (alpha-new 0d0))
                (setf sum 0d0)
                (dotimes (i m)
                  (let ((value (- (* (aref p i) a-down) (* alpha-now (aref u i)))))
                    (setf (aref u i) value)
                    (incf sum (* value value))))
                (setf beta (ls-norm u sum))
                (unless (= beta 0d0)
                  (dotimes (i m) (setf (aref u i) (/ (aref u i) beta)))
                  (csc:transpose-matvec! matrix u q)
                  (setf sum 0d0)
                  (dotimes (j n)
                    (let ((value (- (* (aref q j) a-down) (* beta (aref v j)))))
                      (setf (aref v j) value)
                      (incf sum (* value value))))
                  (setf alpha-new (ls-norm v sum))
                  (unless (= alpha-new 0d0)
                    (dotimes (j n) (setf (aref v j) (/ (aref v j) alpha-new)))))
                (let* ((rho-tilde (ls-hypotenuse rho-bar damp)) (rho (ls-hypotenuse rho-tilde beta))
                       (psi (* (/ damp rho-tilde) phi-bar)) (phi-tilde (* (/ rho-bar rho-tilde) phi-bar))
                       (c (/ rho-tilde rho)) (sine (/ beta rho)) (theta (* sine alpha-new)) (phi (* c phi-tilde))
                       (forward (/ phi rho)) (carry (/ theta rho)))
                  (setf rho-bar (- 0d0 (* c alpha-new)) phi-bar (* sine phi-tilde))
                  (setf damped (+ damped (* psi psi)))
                  (setf frobenius (+ frobenius (+ (* alpha-now alpha-now) (* beta beta))))
                  (setf sum 0d0)
                  (dotimes (j n)
                    (let ((value (+ (aref x j) (* forward (aref w j)))))
                      (setf (aref x j) value)
                      (incf sum (* value value))
                      (setf (aref w j) (- (aref v j) (* carry (aref w j))))))
                  (setf alpha alpha-new)
                  (incf k)
                  (setf estimate (* (sqrt frobenius) a-up))
                  (let* ((pull (* damp (cond ((= damp 0d0) 0d0)
                                             ((< sum sb-ext:double-float-positive-infinity) (sqrt sum))
                                             (t (ls-norm x sum)))))
                         (left-over (- (+ (* phi-bar phi-bar) damped) (* pull pull))))
                    (push (* (cond ((= damp 0d0) (abs phi-bar)) ((> left-over 0d0) (sqrt left-over)) (t 0d0)) b-up)
                          residuals)
                    (push (ls-twice-scaled (* (* alpha-new beta) (/ (abs phi) rho)) a-up b-up) gradients))
                  (push (ls-scaled-copy x back) iterates)
                  (when (or (= beta 0d0) (= alpha-new 0d0) (= rho-bar 0d0)) (setf ended t)))))
            (unless confirmed (confirm))
            (cond ((within formed-r divided-r residual-tolerance) (done 0))
                  ((within formed-g divided-g gradient-tolerance) (done 1))
                  (ended (done 2))
                  (t (done 3)))))))))

;;; Exact quantities, from the doubles of the matrix, the right-hand side and a vector.

(defun ls-stored (a)
  (rect (array-dimension a 0) (array-dimension a 1) (lambda (i j) (rational (coerce (aref a i j) 'double-float)))))

(defun ls-csc (a)
  "The nonzero entries of a matrix of numbers as a checked CSC matrix of doubles."
  (let ((m (array-dimension a 0)) (n (array-dimension a 1)))
    (csc-from-entries m n (loop :for j :below n :append
                            (loop :for i :below m :unless (zerop (aref a i j)) :collect (list i j (aref a i j)))))))

(defun ls-exact-residual (a b x)
  "b - A x as a list of rationals."
  (loop :for i :below (array-dimension a 0)
        :collect (- (rational (aref b i))
                    (loop :for j :below (array-dimension a 1) :sum (* (aref a i j) (rational (aref x j)))))))

(defun ls-exact-gradient (a residual x damping)
  "A^T r - damping^2 x as a list of rationals."
  (loop :for j :below (array-dimension a 1)
        :collect (- (loop :for i :below (array-dimension a 0) :for r :in residual :sum (* (aref a i j) r))
                    (* (expt (rational damping) 2) (rational (aref x j))))))

(defun ls-squares (values) (loop :for value :in values :sum (* value value)))

(defun ls-root (square) (sqrt (coerce square 'double-float)))

(defun ls-exact-minimizer (a b damping)
  "The minimizer of |b - A x|^2 + damping^2 |x|^2 as a list of rationals; of minimum norm where damping is zero."
  (let ((m (array-dimension a 0)) (n (array-dimension a 1)) (rhs (map 'list #'rational b)))
    (if (zerop damping)
        (exact-minimum-norm a rhs)
        (let* ((shift (expt (rational damping) 2))
               (inverse (rational-inverse
                         (build-matrix n (lambda (i j)
                                           (+ (if (= i j) shift 0)
                                              (loop :for k :below m :sum (* (aref a k i) (aref a k j))))))))
               (product (loop :for j :below n :collect (loop :for i :below m :for r :in rhs :sum (* (aref a i j) r)))))
          (loop :for i :below n :collect (loop :for j :below n :for g :in product :sum (* (aref inverse i j) g)))))))

(defun ls-null-part-squared (a x)
  "The squared norm of the part of x outside the range of A^T: x less the minimum-norm solution of A z = A x."
  (let* ((m (array-dimension a 0)) (n (array-dimension a 1))
         (image (loop :for i :below m :collect (loop :for j :below n :sum (* (aref a i j) (rational (aref x j))))))
         (inside (exact-minimum-norm a image)))
    (loop :for j :below n :for z :in inside :sum (expt (- (rational (aref x j)) z) 2))))

(defun ls-case (state m n rank &key zero-row zero-column consistent)
  "A matrix of small integers of rank at most RANK, as a product, and a right-hand side of small integers, in its
range when CONSISTENT."
  (let* ((left (rect m rank (lambda (i k) (declare (ignore i k)) (- (random 7 state) 3))))
         (right (rect rank n (lambda (k j) (declare (ignore k j)) (- (random 7 state) 3))))
         (a (rect-product left right)))
    (when zero-row (dotimes (j n) (setf (aref a zero-row j) 0)))
    (when zero-column (dotimes (i m) (setf (aref a i zero-column) 0)))
    (values a
            (if consistent
                (let ((from (loop :repeat n :collect (- (random 7 state) 3))))
                  (apply #'doubles (loop :for i :below m
                                         :collect (coerce (loop :for j :below n :for z :in from :sum (* (aref a i j) z))
                                                          'double-float))))
                (apply #'doubles (loop :repeat m :collect (coerce (- (random 9 state) 4) 'double-float)))))))

(defun ls-frobenius (a)
  (ls-root (loop :for i :below (array-dimension a 0)
                 :sum (loop :for j :below (array-dimension a 1) :sum (expt (aref a i j) 2)))))

(defun check-lsqr (label a b damping residual-tolerance gradient-tolerance &key (limit 200) rank)
  "Solve, and hold the outcome against exact arithmetic and against the recomputed iteration. Returns the steps, the
exact error of the returned vector over 1 + the minimizer's norm, the exact null-space part over the vector's norm,
and the iterates."
  (let* ((m (array-dimension a 0)) (n (array-dimension a 1)) (matrix (ls-csc a)) (exact (ls-stored a))
         (workspace (lsqr:make-lsqr m n limit)) (x (make-array n :element-type 'double-float :initial-element 9d0))
         (apply (rf-caller:matrix-operator matrix)) (apply-transpose (ls-caller:transpose-operator matrix))
         (code (ls-caller:lsqr-code workspace apply apply-transpose b x damping residual-tolerance gradient-tolerance))
         (steps (lsqr:step-count workspace))
         (frobenius (ls-frobenius a)) (b-norm (two-norm b)) (x-norm (two-norm x)))
    (assert (<= code 1) () "~A: outcome ~D after ~D steps; residual ~S, gradient ~S." label code steps
            (lsqr:residual-norm workspace) (lsqr:gradient-norm workspace))
    (assert (and (= m (lsqr:rows workspace)) (= n (lsqr:columns workspace)) (= limit (lsqr:step-limit workspace))))
    ;; The norms reported are those of the vector returned: the exact ones agree to rounding.
    (let* ((residual (ls-exact-residual exact b x))
           (gradient (ls-exact-gradient exact residual x damping))
           (exact-r (ls-root (ls-squares residual))) (exact-g (ls-root (ls-squares gradient)))
           (slack-r (* (scale-float 1d0 -40) (+ b-norm (* frobenius x-norm))))
           (slack-g (* (scale-float 1d0 -40) (+ (* frobenius (+ b-norm (* frobenius x-norm))) (* damping damping x-norm))))
           (reported-r (lsqr:residual-norm workspace)) (reported-g (lsqr:gradient-norm workspace)))
      (assert (<= (abs (- reported-r exact-r)) slack-r) () "~A: reported residual ~S, exact ~S." label reported-r exact-r)
      (assert (<= (abs (- reported-g exact-g)) slack-g) () "~A: reported gradient ~S, exact ~S." label reported-g exact-g)
      ;; The outcome says which tolerance was met, of the norms as formed.
      (if (= code 0)
          (assert (<= reported-r residual-tolerance) () "~A: residual met at ~S against ~S." label reported-r
                  residual-tolerance)
          (assert (and (<= reported-g gradient-tolerance) (> reported-r residual-tolerance)) ()
                  "~A: gradient met at ~S against ~S with residual ~S against ~S." label reported-g gradient-tolerance
                  reported-r residual-tolerance)))
    ;; The iteration recomputed: the outcome, the counts, both histories, both reports and the vector.
    (multiple-value-bind (reference-code reference-x reference-steps reference-confirmations residuals gradients
                          formed-r formed-g estimate iterates)
        (reference-lsqr matrix m n b damping residual-tolerance gradient-tolerance limit)
      (assert (and (= code reference-code) (= steps reference-steps)
                   (= (lsqr:confirmation-count workspace) reference-confirmations))
              () "~A: outcome ~D in ~D steps and ~D confirmations; recomputed ~D, ~D and ~D." label code steps
              (lsqr:confirmation-count workspace) reference-code reference-steps reference-confirmations)
      (assert (same-bits-p x reference-x) () "~A: the returned vector ~S is not the recomputed one ~S." label x
              reference-x)
      (assert (same-bits-p (subseq (lsqr:recurrence-residual-norms workspace) 0 (1+ steps))
                           (coerce residuals '(simple-array double-float (*))))
              () "~A: running residuals ~S, recomputed ~S." label
              (subseq (lsqr:recurrence-residual-norms workspace) 0 (1+ steps)) residuals)
      (assert (same-bits-p (subseq (lsqr:recurrence-gradient-norms workspace) 0 (1+ steps))
                           (coerce gradients '(simple-array double-float (*))))
              () "~A: running gradients ~S, recomputed ~S." label
              (subseq (lsqr:recurrence-gradient-norms workspace) 0 (1+ steps)) gradients)
      (assert (same-bits-p (doubles (lsqr:residual-norm workspace) (lsqr:gradient-norm workspace)
                                    (lsqr:operator-norm-estimate workspace))
                           (doubles formed-r formed-g estimate))
              () "~A: reports ~S ~S ~S, recomputed ~S ~S ~S." label (lsqr:residual-norm workspace)
              (lsqr:gradient-norm workspace) (lsqr:operator-norm-estimate workspace) formed-r formed-g estimate)
      ;; The estimate of the operator's norm is from below, to rounding, for as many steps as the sequences
      ;; would have run in exact arithmetic. Past that it goes on growing, and is not held to anything here.
      (when (and rank (<= steps rank))
        (assert (<= (lsqr:operator-norm-estimate workspace) (* frobenius (+ 1d0 (scale-float 1d0 -40)))) ()
                "~A: an operator norm estimate of ~S for a Frobenius norm of ~S after ~D steps." label
                (lsqr:operator-norm-estimate workspace) frobenius steps))
      ;; Against the exact minimizer, and for the part of the vector the operator cannot see.
      (let* ((minimizer (ls-exact-minimizer exact b damping))
             (minimizer-norm (ls-root (ls-squares minimizer)))
             (distance (ls-root (loop :for j :below n :for z :in minimizer :sum (expt (- (rational (aref x j)) z) 2))))
             (outside (if (zerop damping) (ls-root (ls-null-part-squared exact x)) 0d0)))
        (values steps (/ distance (+ 1d0 minimizer-norm)) (if (zerop x-norm) 0d0 (/ outside x-norm)) iterates
                minimizer)))))

(defun check-lsqr-monotone (label a b damping iterates minimizer)
  "What the method is, measured exactly on every iterate: the damped objective does not rise, the norm does not fall
and the distance to the minimizer does not rise. A part in 2^30 is allowed for rounding, and a distance already at
rounding level is not compared."
  (let* ((exact (ls-stored a)) (n (array-dimension a 1)) (shift (expt (rational damping) 2))
         (floor (* (expt 2 -70) (+ 1 (ls-squares minimizer)))) (slack (+ 1 (expt 2 -30))))
    (loop :for (earlier later) :on (mapcar (lambda (x)
                                             (let ((length (loop :for value :across x :sum (expt (rational value) 2))))
                                               (list (+ (ls-squares (ls-exact-residual exact b x)) (* shift length))
                                                     length
                                                     (loop :for j :below n :for z :in minimizer
                                                           :sum (expt (- (rational (aref x j)) z) 2)))))
                                           iterates)
          :while later
          :do (assert (<= (first later) (* (first earlier) slack)) () "~A: the damped objective rose from ~E to ~E."
                      label (float (first earlier) 1d0) (float (first later) 1d0))
              (assert (>= (* (second later) slack) (second earlier)) () "~A: the norm fell from ~E to ~E." label
                      (float (second earlier) 1d0) (float (second later) 1d0))
              (when (> (third earlier) floor)
                (assert (<= (third later) (* (third earlier) slack)) () "~A: the distance rose from ~E to ~E." label
                        (float (third earlier) 1d0) (float (third later) 1d0)))
              ;; The residual alone: the objective less the damped term. It does not rise either, since
              ;; the objective does not rise and the norm does not fall.
              (let ((before (- (first earlier) (* shift (second earlier)))) (after (- (first later) (* shift (second later))))
                    (resolved (* (expt 2 -70) (+ 1 (ls-squares (mapcar #'rational (coerce b 'list)))))))
                (when (> before resolved)
                  (assert (<= after (+ (* before slack) resolved)) () "~A: the residual rose from ~E to ~E." label
                          (float before 1d0) (float after 1d0)))))))

(defparameter *lsqr-shapes*
  '((6 4 4) (4 4 4) (4 6 4) (7 5 3) (5 7 3) (5 5 2) (3 1 1) (1 3 1) (8 6 6) (6 6 1)
    (6 5 3 :zero-row 2 :zero-column 1) (5 6 3 :zero-row 0 :zero-column 5)))

(defun run-lsqr-exact-tests ()
  (let ((state (sb-ext:seed-random-state 20261123)) (cases 0) (most 0) (worst-distance 0d0) (worst-outside 0d0)
        (residual-met 0) (iterates-measured 0) (short-calls 0))
    (dolist (shape *lsqr-shapes*)
      (destructuring-bind (m n rank &rest options) shape
        (dotimes (repeat 2)
          (dolist (damping '(0d0 0.25d0 1d0 8d0))
            (dolist (consistent (if (zerop damping) '(nil t) '(nil)))
              (multiple-value-bind (a b) (apply #'ls-case state m n rank :consistent consistent options)
                (let* ((label (format nil "~D by ~D rank ~D damping ~A~:[~; consistent~] #~D" m n rank damping consistent
                                      repeat))
                       (b-norm (two-norm b)) (frobenius (ls-frobenius a))
                       ;; A consistent system is asked for its residual; the others for their gradient.
                       (residual-tolerance (if consistent (* (scale-float 1d0 -30) b-norm) 0d0))
                       (gradient-tolerance (if consistent 0d0 (* (scale-float 1d0 -30) frobenius b-norm))))
                  (multiple-value-bind (steps distance outside iterates minimizer)
                      (check-lsqr label a b damping residual-tolerance gradient-tolerance :rank rank)
                    ;; In exact arithmetic the sequences end within as many steps as the rank; rounding may add some.
                    (assert (<= steps (+ rank 4)) () "~A took ~D steps." label steps)
                    (assert (<= distance (scale-float 1d0 -30)) () "~A: the vector is ~E from the minimizer." label distance)
                    (assert (<= outside (scale-float 1d0 -40)) () "~A: ~E of the vector is outside the range of the transpose."
                            label outside)
                    (check-lsqr-monotone label a b damping iterates minimizer)
                    (incf iterates-measured (length iterates))
                    ;; Each iterate is what a call with that step limit returns.
                    (when (zerop (mod cases 7))
                      (let ((matrix (ls-csc a)))
                        (loop :for k :from 0 :below steps :for iterate :in iterates :do
                          (let ((workspace (lsqr:make-lsqr m n k)) (x (zeros n)))
                            (ls-caller:lsqr-code workspace (rf-caller:matrix-operator matrix)
                                                 (ls-caller:transpose-operator matrix) b x damping 0d0 0d0)
                            (assert (same-bits-p x iterate) () "~A: a call limited to ~D steps gave ~S, iterate ~D is ~S."
                                    label k x k iterate)
                            (incf short-calls)))))
                    (when consistent (incf residual-met))
                    (setf most (max most steps) worst-distance (max worst-distance distance)
                          worst-outside (max worst-outside outside))
                    (incf cases)))))))))
    (format t "LSQR-EXACT-PASS cases=~D consistent=~D most-steps=~D iterates-measured=~D short-calls=~D worst-distance=~,3E worst-outside=~,3E~%"
            cases residual-met most iterates-measured short-calls worst-distance worst-outside)))

(defun ls-agree (label a b damping residual-tolerance gradient-tolerance limit)
  "A call on a stored matrix that must be the recomputed iteration: the outcome, the counts, the reports and the vector,
bit for bit. Returns the code, the vector and the workspace."
  (let* ((m (array-dimension a 0)) (n (array-dimension a 1)) (matrix (ls-csc a))
         (workspace (lsqr:make-lsqr m n limit)) (x (make-array n :element-type 'double-float :initial-element 9d0))
         (code (ls-caller:lsqr-code workspace (rf-caller:matrix-operator matrix) (ls-caller:transpose-operator matrix)
                                    b x damping residual-tolerance gradient-tolerance)))
    (multiple-value-bind (reference-code reference-x steps confirmations residuals gradients formed-r formed-g estimate)
        (reference-lsqr matrix m n b damping residual-tolerance gradient-tolerance limit)
      (declare (ignore residuals gradients))
      (assert (and (= code reference-code) (= steps (lsqr:step-count workspace))
                   (= confirmations (lsqr:confirmation-count workspace)) (same-bits-p x reference-x)
                   (same-bits-p (doubles (lsqr:residual-norm workspace) (lsqr:gradient-norm workspace)
                                         (lsqr:operator-norm-estimate workspace))
                                (doubles formed-r formed-g estimate)))
              () "~A: outcome ~D, ~D steps, ~D confirmations, ~S, reports ~S ~S ~S; recomputed ~D, ~D, ~D, ~S, ~S ~S ~S."
              label code (lsqr:step-count workspace) (lsqr:confirmation-count workspace) x (lsqr:residual-norm workspace)
              (lsqr:gradient-norm workspace) (lsqr:operator-norm-estimate workspace) reference-code steps confirmations
              reference-x formed-r formed-g estimate))
    (values code x workspace)))

(defun run-lsqr-outcome-tests ()
  (let* ((nan (not-a-number)) (infinity sb-ext:double-float-positive-infinity)
         (a (rect-of-rows '((4 1 0) (1 3 1) (0 1 2) (1 0 1)))) (matrix (ls-csc a)) (exact (ls-stored a))
         (forward (rf-caller:matrix-operator matrix)) (backward (ls-caller:transpose-operator matrix))
         (b (doubles 1d0 2d0 3d0 4d0)) (nines (doubles 9d0 9d0 9d0))
         ;; A^T b is (10, 10, 12).
         (b-norm (two-norm b)) (start-gradient (two-norm (doubles 10d0 10d0 12d0)))
         (drifting-confirmations 0) (drifting-steps 0)
         (small-damping-steps 0) (small-damping-confirmations 0) (small-damping-zeros 0)
         (slight-steps 0) (slight-confirmations 0) (slight-below 0))
    (flet ((run (workspace apply apply-transpose rhs x damping residual-tolerance gradient-tolerance)
             (ls-caller:lsqr-code workspace apply apply-transpose rhs x damping residual-tolerance gradient-tolerance))
           (unreported-p (workspace)
             (and (= -1d0 (lsqr:residual-norm workspace)) (= -1d0 (lsqr:gradient-norm workspace))))
           (near (value expected) (<= (abs (- value expected)) (* 1d-14 (abs expected)))))
      ;; Before any call.
      (let ((workspace (lsqr:make-lsqr 4 3 20)) (x (copy-seq nines)) (calls (indices 0)))
        (assert (and (= 0 (lsqr:step-count workspace)) (= 0 (lsqr:confirmation-count workspace)) (unreported-p workspace)
                     (= -1d0 (lsqr:operator-norm-estimate workspace))))
        ;; Sizes, tolerances and the damping are refused before anything is read or written.
        (assert (= 6 (run workspace forward backward (doubles 1d0 2d0 3d0) x 0d0 0d0 0d0)))
        (assert (= 6 (run workspace forward backward b (doubles 9d0 9d0) 0d0 0d0 0d0)))
        (dolist (bad (list -1d0 nan))
          (assert (= 7 (run workspace forward backward b x 0d0 bad 0d0)))
          (assert (= 7 (run workspace forward backward b x 0d0 0d0 bad))))
        (dolist (bad (list -1d0 nan infinity))
          (assert (= 8 (run workspace forward backward b x bad 0d0 0d0)) () "A damping of ~S was taken." bad))
        ;; A right-hand side that is not finite, and a first image that is not: x is still untouched.
        (assert (= 4 (run workspace forward backward (doubles 1d0 nan 3d0 4d0) x 0d0 0d0 0d0)))
        (assert (= 4 (run workspace forward backward (doubles 1d0 infinity 3d0 4d0) x 0d0 0d0 0d0)))
        (assert (= 4 (run workspace forward (cg-caller:drifting-operator backward 0 nan calls) b x 0d0 0d0 0d0)))
        (assert (and (same-bits-p x nines) (unreported-p workspace) (= 0 (lsqr:step-count workspace))) ()
                "A refused call left ~S, ~D steps, residual ~S." x (lsqr:step-count workspace)
                (lsqr:residual-norm workspace))
        ;; A sound call, then a refused one: the earlier report is withdrawn.
        (assert (= 1 (run workspace forward backward b x 0d0 0d0 1d-9)))
        (assert (and (plusp (lsqr:step-count workspace)) (plusp (lsqr:confirmation-count workspace))
                     (<= 0d0 (lsqr:gradient-norm workspace) 1d-9) (plusp (lsqr:residual-norm workspace))
                     (plusp (lsqr:operator-norm-estimate workspace))))
        (assert (= 7 (run workspace forward backward b x 0d0 -1d0 0d0)))
        (assert (and (= 0 (lsqr:step-count workspace)) (= 0 (lsqr:confirmation-count workspace)) (unreported-p workspace)
                     (= -1d0 (lsqr:operator-norm-estimate workspace)))))
      ;; Sizes the format cannot hold together: x untouched.
      (let ((workspace (lsqr:make-lsqr 1 1 5)) (x (doubles 9d0)))
        (flet ((scaled (factor) (cg-caller:scaling-operator factor)))
          (assert (= 5 (run workspace (scaled 1d-300) (scaled 1d-300) (doubles 1d300) x 0d0 0d0 0d0)))
          (assert (= 5 (run workspace (scaled 1d300) (scaled 1d300) (doubles 1d-300) x 0d0 0d0 0d0)))
          (assert (= 5 (run workspace (scaled 1d-300) (scaled 1d-300) (doubles 1d0) x 1d300 0d0 0d0)))
          (assert (and (= 9d0 (aref x 0)) (unreported-p workspace)))
          ;; An operator and a right-hand side of 1e200 are solved, and exactly.
          (assert (= 0 (run workspace (scaled 1d200) (scaled 1d200) (doubles 1d200) x 0d0 0d0 0d0)))
          (assert (and (= 1d0 (aref x 0)) (= 0d0 (lsqr:residual-norm workspace)) (= 1 (lsqr:step-count workspace))) ()
                  "An operator and a right-hand side of 1e200 gave ~S with residual ~S." x (lsqr:residual-norm workspace))
          ;; A gradient that is not zero and too small to hold does not meet a tolerance of zero, and
          ;; meets any other.
          (let ((none (lsqr:make-lsqr 1 1 0)) (tiny (scale-float 1d0 -600)))
            (assert (= 3 (run none (scaled tiny) (scaled tiny) (doubles tiny) x 0d0 0d0 0d0)))
            (assert (and (= 0d0 (lsqr:gradient-norm none)) (= tiny (lsqr:residual-norm none)) (= 0d0 (aref x 0))))
            (assert (= 1 (run none (scaled tiny) (scaled tiny) (doubles tiny) x 0d0 0d0 least-positive-double-float))))))
      ;; A zero right-hand side: zero is the answer and no product is made.
      (let ((workspace (lsqr:make-lsqr 4 3 20)) (x (copy-seq nines)) (calls (indices 0)))
        (assert (= 0 (run workspace (cg-caller:drifting-operator forward 0 nan calls)
                          (cg-caller:drifting-operator backward 0 nan calls) (zeros 4) x 0d0 0d0 0d0)))
        (assert (and (same-bits-p x (zeros 3)) (= 0 (aref calls 0)) (= 0d0 (lsqr:residual-norm workspace))
                     (= 0d0 (lsqr:gradient-norm workspace)) (= 0 (lsqr:step-count workspace))
                     (= 0 (lsqr:confirmation-count workspace)) (= -1d0 (lsqr:operator-norm-estimate workspace)))))
      ;; A right-hand side the transpose sends to zero: zero is a minimizer, found with no step.
      (multiple-value-bind (code x workspace) (ls-agree "orthogonal" (rect-of-rows '((1) (0))) (doubles 0d0 3d0) 0d0 0d0 0d0 5)
        (assert (and (= 1 code) (= 0d0 (aref x 0)) (= 0 (lsqr:step-count workspace)) (= 1 (lsqr:confirmation-count workspace))
                     (= 3d0 (lsqr:residual-norm workspace)) (= 0d0 (lsqr:gradient-norm workspace)))))
      ;; The sequences end and the arithmetic is exact: one by one, with and without a damping.
      (multiple-value-bind (code x workspace) (ls-agree "one by one" (rect-of-rows '((1))) (doubles 3d0) 0d0 0d0 0d0 5)
        (assert (and (= 0 code) (= 3d0 (aref x 0)) (= 1 (lsqr:step-count workspace)) (= 0d0 (lsqr:residual-norm workspace)))))
      (multiple-value-bind (code x workspace) (ls-agree "one by one, damped" (rect-of-rows '((1))) (doubles 3d0) 1d0 0d0 0d0 5)
        (assert (and (member code '(1 2)) (near (aref x 0) 1.5d0) (= 1 (lsqr:step-count workspace))
                     (<= (lsqr:gradient-norm workspace) 1d-15))
                () "One by one with a damping of 1: outcome ~D, ~S, gradient ~S." code x (lsqr:gradient-norm workspace)))
      ;; The sequences end, by a transpose that returns zero once, at a vector that is not a minimizer:
      ;; that is said, and the report is of the vector returned.
      (let ((workspace (lsqr:make-lsqr 4 3 20)) (x (copy-seq nines)) (calls (indices 0)))
        (assert (= 2 (run workspace forward (cg-caller:drifting-operator backward 0 0d0 calls) b x 0d0 0d0 0d0)))
        (assert (and (same-bits-p x (zeros 3)) (= 0 (lsqr:step-count workspace)) (= 1 (lsqr:confirmation-count workspace))
                     (near (lsqr:residual-norm workspace) b-norm) (near (lsqr:gradient-norm workspace) start-gradient))
                () "After a transpose that returned zero once: ~S, residual ~S, gradient ~S." x
                (lsqr:residual-norm workspace) (lsqr:gradient-norm workspace)))
      ;; A limit of zero judges the zero vector: neither tolerance met, then each met.
      (multiple-value-bind (code x workspace) (ls-agree "no step" a b 0d0 0d0 0d0 0)
        (assert (and (= 3 code) (same-bits-p x (zeros 3)) (= 0 (lsqr:step-count workspace))
                     (= 1 (lsqr:confirmation-count workspace)) (near (lsqr:residual-norm workspace) b-norm)
                     (near (lsqr:gradient-norm workspace) start-gradient))))
      (assert (= 0 (ls-agree "no step, residual met" a b 0d0 10d0 0d0 0)))
      (assert (= 1 (ls-agree "no step, gradient met" a b 0d0 0d0 20d0 0)))
      (assert (= 1 (ls-agree "any gradient" a b 0d0 0d0 infinity 20)))
      ;; One step allowed: the report is the residual and the gradient of what is returned.
      (multiple-value-bind (code x workspace) (ls-agree "one step" a b 0.5d0 0d0 0d0 1)
        (let* ((residual (ls-exact-residual exact b x)) (gradient (ls-exact-gradient exact residual x 0.5d0)))
          (assert (and (= 3 code) (= 1 (lsqr:step-count workspace)) (= 1 (lsqr:confirmation-count workspace))
                       (near (lsqr:residual-norm workspace) (ls-root (ls-squares residual)))
                       (near (lsqr:gradient-norm workspace) (ls-root (ls-squares gradient)))
                       (> (lsqr:gradient-norm workspace) 1d-3))
                  () "At a limit of one: outcome ~D, residual ~S, gradient ~S." code (lsqr:residual-norm workspace)
                  (lsqr:gradient-norm workspace))
          ;; An image that is not a number in the second step: what is returned is the first iterate.
          (let ((spoiled (lsqr:make-lsqr 4 3 20)) (y (copy-seq nines)) (calls (indices 0)))
            (assert (= 4 (run spoiled (cg-caller:drifting-operator forward 1 nan calls) backward b y 0.5d0 0d0 0d0)))
            (assert (and (same-bits-p y x) (= 1 (lsqr:step-count spoiled)) (unreported-p spoiled)) ()
                    "After a spoiled second step: ~S, the first iterate is ~S." y x))
          ;; In the first step: zero.
          (let ((spoiled (lsqr:make-lsqr 4 3 20)) (y (copy-seq nines)) (calls (indices 0)))
            (assert (= 4 (run spoiled (cg-caller:drifting-operator forward 0 nan calls) backward b y 0.5d0 0d0 0d0)))
            (assert (and (same-bits-p y (zeros 3)) (= 0 (lsqr:step-count spoiled)) (unreported-p spoiled))))))
      ;; A damping far above the operator: up to 2^500 times its size the problem is solved, the answer
      ;; being A^T b over the damping's square to a part in 2^900 or so; above that the iterate could
      ;; not be held in the iteration's units, and the call says so with x untouched.
      (let* ((damping (scale-float 1d0 480))
             (minimizer (ls-exact-minimizer exact b damping)))
        (multiple-value-bind (code x workspace) (ls-agree "a large damping" a b damping 0d0 (scale-float 1d0 -30) 20)
          (assert (and (= 1 code) (<= (lsqr:step-count workspace) 2)
                       (loop :for j :below 3 :for z :in minimizer
                             :always (<= (abs (- (rational (aref x j)) z)) (* (expt 2 -40) (abs z)))))
                  () "A damping of 2^480: outcome ~D in ~D steps, ~S against ~S." code (lsqr:step-count workspace) x
                  (mapcar (lambda (z) (float z 1d0)) minimizer))))
      (let ((workspace (lsqr:make-lsqr 4 3 20)) (x (copy-seq nines)))
        (assert (= 5 (run workspace forward backward b x (scale-float 1d0 520) 0d0 0d0)))
        (assert (and (same-bits-p x nines) (unreported-p workspace))))
      ;; Vectors far below unit size, whose squares cannot be held. The operator is diag(1, 2^-200) and
      ;; the right-hand side (1, 2^-600): the second left vector has the one entry -2^-600 before it is
      ;; normalized, and after one step the gradient is (0, 2^-800). Their norms are those numbers, not
      ;; the zero their squares round to: the gradient does not meet a tolerance of zero, and the second
      ;; step is taken and solves the system.
      (let ((graded (rect-of-rows (list (list 1d0 0d0) (list 0d0 (scale-float 1d0 -200)))))
            (small (doubles 1d0 (scale-float 1d0 -600))))
        (multiple-value-bind (code x workspace) (ls-agree "a tiny gradient" graded small 0d0 0d0 0d0 1)
          (assert (and (= 3 code) (= 1 (lsqr:step-count workspace)) (= (scale-float 1d0 -800) (lsqr:gradient-norm workspace))
                       (= (scale-float 1d0 -600) (lsqr:residual-norm workspace))
                       (same-bits-p x (doubles 1d0 (scale-float 1d0 -800))))
                  () "After one step on the graded operator: outcome ~D, ~S, residual ~S, gradient ~S." code x
                  (lsqr:residual-norm workspace) (lsqr:gradient-norm workspace)))
        (multiple-value-bind (code x workspace) (ls-agree "a tiny left vector" graded small 0d0 0d0 0d0 5)
          (assert (and (= 0 code) (<= 2 (lsqr:step-count workspace) 5) (= 0d0 (lsqr:residual-norm workspace))
                       (same-bits-p x (doubles 1d0 (scale-float 1d0 -400))))
                  () "On the graded operator: outcome ~D in ~D steps, ~S, residual ~S." code (lsqr:step-count workspace) x
                  (lsqr:residual-norm workspace))))
      ;; One wrong product under a residual tolerance, on a system that has a solution: the same three
      ;; things are asked of the residual's estimate.
      (let ((state (sb-ext:seed-random-state 81)))
        (multiple-value-bind (square b6) (ls-case state 6 6 6 :consistent t)
          (let* ((stored (ls-csc square)) (exact6 (ls-stored square))
                 (apply (rf-caller:matrix-operator stored)) (apply-transpose (ls-caller:transpose-operator stored))
                 (tolerance (* (scale-float 1d0 -30) (two-norm b6)))
                 (honest (lsqr:make-lsqr 6 6 40)) (drifting (lsqr:make-lsqr 6 6 40))
                 (x (zeros 6)) (y (zeros 6)) (calls (indices 0)))
            (assert (= 0 (run honest apply apply-transpose b6 x 0d0 tolerance 0d0)))
            (assert (= 1 (lsqr:confirmation-count honest)) () "The honest call took ~D confirmations."
                    (lsqr:confirmation-count honest))
            (let* ((code (run drifting (cg-caller:drifting-operator apply 1 1.001d0 calls) apply-transpose b6 y 0d0
                              tolerance 0d0))
                   (exact-residual (ls-root (ls-squares (ls-exact-residual exact6 b6 y))))
                   (confirmations (lsqr:confirmation-count drifting)) (steps (lsqr:step-count drifting)))
              (assert (> exact-residual (* 1000 tolerance)) () "The wrong product's premise failed: residual ~E against ~E."
                      exact-residual tolerance)
              (assert (member code '(2 3)) () "After a wrong product the call reported outcome ~D with a residual of ~E against ~E."
                      code exact-residual tolerance)
              (assert (<= (abs (- (lsqr:residual-norm drifting) exact-residual)) (* 1d-11 exact-residual)) ()
                      "After a wrong product: reported residual ~S, exact ~S." (lsqr:residual-norm drifting) exact-residual)
              (assert (and (>= confirmations 2) (<= (* 3 confirmations) steps)) ()
                      "One wrong product under a residual tolerance cost ~D confirmations in ~D steps." confirmations
                      steps)))))
      ;; The sequences end inside the iteration. The operator is one column, (1, k), and the right-hand
      ;; side (1, 0): the second right vector is exactly zero, after one step. With tolerances of zero
      ;; the call ends there, whether the gradient it forms is exactly zero or a rounding away from it,
      ;; and for some k it is a rounding away.
      (let ((exhausted 0))
        (dolist (k '(1 2 3 5 7))
          (multiple-value-bind (code x workspace)
              (ls-agree (format nil "one column, (1, ~D)" k) (rect-of-rows (list (list 1) (list k))) (doubles 1d0 0d0) 0d0
                        0d0 0d0 10)
            (assert (and (member code '(1 2)) (= 1 (lsqr:step-count workspace)) (near (aref x 0) (/ 1d0 (+ 1 (* k k))))
                         (<= (lsqr:gradient-norm workspace) 1d-15))
                    () "One column (1, ~D): outcome ~D in ~D steps, ~S, gradient ~S." k code (lsqr:step-count workspace) x
                    (lsqr:gradient-norm workspace))
            (when (= code 2) (incf exhausted))))
        (assert (>= exhausted 2) () "Of five one-column problems ~D ended at a gradient that is not zero." exhausted))
      ;; An estimate that is exactly zero step after step while what it estimates is not. The operator
      ;; is diag(1, 2^-200), the right-hand side (1, 2^-600) and the damping 2^-700: every term of the
      ;; residual's estimate is a square too small to hold. With tolerances of zero the first such
      ;; estimate buys a confirmation, which shows the residual is not zero, and the estimate is then
      ;; switched off: three estimates of zero cost one confirmation, not three. The second confirmation
      ;; is the gradient's, at the third step, and finds the system solved.
      (multiple-value-bind (code x workspace)
          (ls-agree "an estimate of zero" (rect-of-rows (list (list 1d0 0d0) (list 0d0 (scale-float 1d0 -200))))
                    (doubles 1d0 (scale-float 1d0 -600)) (scale-float 1d0 -700) 0d0 0d0 4)
        (declare (ignore x))
        (setf small-damping-steps (lsqr:step-count workspace)
              small-damping-confirmations (lsqr:confirmation-count workspace)
              small-damping-zeros (count 0d0 (subseq (lsqr:recurrence-residual-norms workspace) 1
                                                     (1+ (lsqr:step-count workspace)))))
        (assert (and (= 0 code) (= 3 small-damping-steps) (= 3 small-damping-zeros) (= 2 small-damping-confirmations)
                     (= 0d0 (lsqr:residual-norm workspace)))
                () "With an estimate of zero: outcome ~D, ~D steps, ~D estimates of zero, ~D confirmations, residual ~S."
                code small-damping-steps small-damping-zeros small-damping-confirmations (lsqr:residual-norm workspace)))
      ;; One wrong product. The running estimates then meet the tolerance while the gradient of the
      ;; vector does not. The call must not report a tolerance met on an estimate; it must say how
      ;; far it got, of the vector returned; and a wrong estimate must not buy a confirmation a step.
      (let ((state (sb-ext:seed-random-state 79)))
        (multiple-value-bind (wide b8) (ls-case state 8 6 6)
          (let* ((stored (ls-csc wide)) (exact8 (ls-stored wide))
                 (apply (rf-caller:matrix-operator stored)) (apply-transpose (ls-caller:transpose-operator stored))
                 (tolerance (* (scale-float 1d0 -30) (ls-frobenius wide) (two-norm b8)))
                 (honest (lsqr:make-lsqr 8 6 40)) (drifting (lsqr:make-lsqr 8 6 40))
                 (x (zeros 6)) (y (zeros 6)) (calls (indices 0)))
            (assert (= 1 (run honest apply apply-transpose b8 x 0d0 0d0 tolerance)))
            (assert (= 1 (lsqr:confirmation-count honest)) () "The honest call took ~D confirmations."
                    (lsqr:confirmation-count honest))
            ;; Product 1 of the operator is the second step's.
            (let* ((code (run drifting (cg-caller:drifting-operator apply 1 1.001d0 calls) apply-transpose b8 y 0d0 0d0
                              tolerance))
                   (residual (ls-exact-residual exact8 b8 y))
                   (exact-gradient (ls-root (ls-squares (ls-exact-gradient exact8 residual y 0d0)))))
              (setf drifting-confirmations (lsqr:confirmation-count drifting) drifting-steps (lsqr:step-count drifting))
              (assert (> exact-gradient (* 1000 tolerance)) () "The wrong product's premise failed: gradient ~E against ~E."
                      exact-gradient tolerance)
              (assert (member code '(2 3)) () "After a wrong product the call reported outcome ~D with a gradient of ~E against ~E."
                      code exact-gradient tolerance)
              ;; The gradient is a small difference of large products here, so it is held to less than `near`.
              (assert (<= (abs (- (lsqr:gradient-norm drifting) exact-gradient)) (* 1d-11 exact-gradient)) ()
                      "After a wrong product: reported ~S, exact ~S."
                      (lsqr:gradient-norm drifting) exact-gradient)
              ;; A product wrong by less, a part in 2^14. The estimates then settle far below the tolerance
              ;; and stay there step after step, while the gradient stays above it: the case in which an
              ;; estimate that went on triggering as before would buy a confirmation at every one of them.
              (let* ((slight (lsqr:make-lsqr 8 6 40)) (z (zeros 6)) (slight-calls (indices 0))
                     (slight-code (run slight (cg-caller:drifting-operator apply 1 (+ 1d0 (scale-float 1d0 -14)) slight-calls)
                                       apply-transpose b8 z 0d0 0d0 tolerance))
                     (slight-gradient
                       (ls-root (ls-squares (ls-exact-gradient exact8 (ls-exact-residual exact8 b8 z) z 0d0)))))
                (setf slight-steps (lsqr:step-count slight) slight-confirmations (lsqr:confirmation-count slight)
                      slight-below (count-if (lambda (value) (<= value tolerance))
                                             (subseq (lsqr:recurrence-gradient-norms slight) 0 (1+ slight-steps))))
                (assert (and (member slight-code '(2 3)) (> slight-gradient (* 1000 tolerance))) ()
                        "After a slightly wrong product: outcome ~D with a gradient of ~E against ~E." slight-code
                        slight-gradient tolerance)
                (assert (and (>= slight-below 20) (<= slight-confirmations 5)) ()
                        "A slightly wrong product: ~D estimates within the tolerance cost ~D confirmations in ~D steps."
                        slight-below slight-confirmations slight-steps))
              (assert (and (>= drifting-confirmations 2) (<= (* 3 drifting-confirmations) drifting-steps)) ()
                      "One wrong product cost ~D confirmations in ~D steps." drifting-confirmations drifting-steps))))))
    (format t "LSQR-OUTCOMES-PASS wrong-product-steps=~D wrong-product-confirmations=~D zero-estimate-steps=~D zero-estimate-confirmations=~D zero-estimates=~D slightly-wrong-steps=~D slightly-wrong-estimates-within=~D slightly-wrong-confirmations=~D~%"
            drifting-steps drifting-confirmations small-damping-steps small-damping-confirmations small-damping-zeros
            slight-steps slight-below slight-confirmations)))

(defun run-lsqr-scale-tests ()
  "The same problem at other powers of two: the right-hand side, the operator with its damping, and both."
  (let ((state (sb-ext:seed-random-state 20261124)) (runs 0))
    (dolist (shape '((6 4 4 0d0) (4 6 4 0.25d0) (7 5 3 0d0) (5 5 5 1d0)))
      (destructuring-bind (m n rank damping) shape
        (multiple-value-bind (a b) (ls-case state m n rank)
          (let* ((matrix (ls-csc a)) (forward (rf-caller:matrix-operator matrix))
                 (backward (ls-caller:transpose-operator matrix))
                 (gradient-tolerance (* (scale-float 1d0 -30) (ls-frobenius a) (two-norm b)))
                 (base (lsqr:make-lsqr m n 50)) (x (zeros n))
                 (code (ls-caller:lsqr-code base forward backward b x damping 0d0 gradient-tolerance))
                 (steps (lsqr:step-count base)))
            (assert (= 1 code))
            (dolist (powers '((900 0) (-900 0) (0 400) (0 -400) (500 -400) (-500 400) (300 300) (-300 -300)))
              (destructuring-bind (of-b of-a) powers
                (let* ((sb (scale-float 1d0 of-b)) (sa (scale-float 1d0 of-a))
                       (workspace (lsqr:make-lsqr m n 50)) (y (zeros n))
                       (scaled-code (ls-caller:lsqr-code workspace (ls-caller:scaled-operator forward sa)
                                                         (ls-caller:scaled-operator backward sa) (ls-scaled-copy b sb) y
                                                         (* damping sa) 0d0 (* (* gradient-tolerance sa) sb)))
                       (label (format nil "~D by ~D, right-hand side by 2^~D and operator by 2^~D" m n of-b of-a)))
                  (flet ((history (workspace reader factor)
                           (ls-scaled-copy (subseq (funcall reader workspace) 0 (1+ steps)) factor)))
                    (assert (and (= code scaled-code) (= steps (lsqr:step-count workspace))
                                 (= (lsqr:confirmation-count base) (lsqr:confirmation-count workspace)))
                            () "~A: outcome ~D in ~D steps against ~D in ~D." label scaled-code
                            (lsqr:step-count workspace) code steps)
                    (assert (same-bits-p y (ls-scaled-copy x (scale-float 1d0 (- of-b of-a)))) ()
                            "~A: ~S is not ~S scaled." label y x)
                    (assert (same-bits-p (doubles (lsqr:residual-norm workspace) (lsqr:gradient-norm workspace)
                                                  (lsqr:operator-norm-estimate workspace))
                                         (doubles (* (lsqr:residual-norm base) sb) (* (* (lsqr:gradient-norm base) sa) sb)
                                                  (* (lsqr:operator-norm-estimate base) sa)))
                            () "~A: reports ~S ~S ~S against ~S ~S ~S." label (lsqr:residual-norm workspace)
                            (lsqr:gradient-norm workspace) (lsqr:operator-norm-estimate workspace)
                            (lsqr:residual-norm base) (lsqr:gradient-norm base) (lsqr:operator-norm-estimate base))
                    (assert (and (same-bits-p (history workspace #'lsqr:recurrence-residual-norms 1d0)
                                              (history base #'lsqr:recurrence-residual-norms sb))
                                 (same-bits-p (history workspace #'lsqr:recurrence-gradient-norms 1d0)
                                              (ls-scaled-copy (history base #'lsqr:recurrence-gradient-norms sa) sb)))
                            () "~A: the histories are not the scaled ones." label))
                  (incf runs))))))))
    (format t "LSQR-SCALE-PASS runs=~D~%" runs)))

(defparameter *lsqr-public-control*
  "(declare inspect (coalton-linalg/lsqr:LsqrWorkspace -> UFix))
   (define (inspect workspace)
     (+ (coalton-linalg/lsqr:rows workspace) (coalton-linalg/lsqr:step-count workspace)))
   (declare probe (coalton-linalg/lsqr:LsqrWorkspace -> UFix))
   (define (probe workspace) (inspect workspace))")

(defparameter *lsqr-encapsulation-probes*
  '(("lsqr-direction-access" .
     "(declare probe (coalton-linalg/lsqr:LsqrWorkspace -> arr:LispArray F64))
      (define (probe workspace) (.direction workspace))")
    ("lsqr-unwrapping" .
     "(declare probe (coalton-linalg/lsqr:LsqrWorkspace -> UFix))
      (define (probe workspace)
        (match workspace ((coalton-linalg/lsqr:MkLsqrWorkspace held) (.limit held))))")
    ("lsqr-status-as-cg-status" .
     "(declare probe (coalton-linalg/lsqr:LsqrStatus -> coalton-linalg/conjugate-gradient:CgStatus))
      (define (probe status) status)")))

(defun ls-range-run (rows b damping residual-tolerance gradient-tolerance limit)
  "Solve with the matrix whose rows are given, and measure the vector returned exactly. Returns the code, the
vector, the workspace, and the exact squared norms of its residual and gradient."
  (let* ((m (length rows)) (n (length (first rows)))
         (a (make-array (list m n) :initial-contents rows)) (matrix (ls-csc a)) (exact (ls-stored a))
         (rhs (coerce b '(simple-array double-float (*))))
         (workspace (lsqr:make-lsqr m n limit)) (x (make-array n :element-type 'double-float :initial-element 9d0))
         (code (ls-caller:lsqr-code workspace (rf-caller:matrix-operator matrix) (ls-caller:transpose-operator matrix)
                                    rhs x damping residual-tolerance gradient-tolerance))
         (finite (every (lambda (value) (< (abs value) sb-ext:double-float-positive-infinity)) x))
         (residual (and finite (ls-exact-residual exact rhs x))))
    (values code x workspace (and finite (ls-squares residual))
            (and finite (ls-squares (ls-exact-gradient exact residual x damping))))))

(defun ls-reported-p (reported exact-square)
  "Whether a norm reported is the exact one to a part in 2^40, compared as squares so that none is lost."
  (let ((square (expt (rational reported) 2)))
    (<= (abs (- square exact-square)) (* exact-square (expt 2 -39)))))

(defun run-lsqr-range-tests ()
  "Sizes inside one problem that one power of two for the right-hand side, and one for the operator taken from
its first image, do not carry. What a status says, and what is reported, is of the vector returned."
  (flet ((two (power) (scale-float 1d0 power)))
    ;; An entry of the right-hand side more than 2^1022 below its largest is not seen by the iteration.
    ;; It is seen when the residual is formed: neither tolerance is met, and both reports are exact.
    (dolist (b (list (list (two 100) (two -980)) (list 1d308 1d-17)))
      (multiple-value-bind (code x workspace residual gradient) (ls-range-run '((1d0 0d0) (0d0 1d0)) b 0d0 0d0 0d0 10)
        (assert (and (member code '(2 3)) (= (first b) (aref x 0)) (= 0d0 (aref x 1))) ()
                "b = ~S: outcome ~D and x = ~S." b code x)
        (assert (and (plusp residual) (ls-reported-p (lsqr:residual-norm workspace) residual)
                     (ls-reported-p (lsqr:gradient-norm workspace) gradient))
                () "b = ~S: a residual of ~S and a gradient of ~S are reported; exactly their squares are ~S and ~S." b
                (lsqr:residual-norm workspace) (lsqr:gradient-norm workspace) (float residual 1d0) (float gradient 1d0))
        ;; A second call on that residual, which the contract gives for an answer to improve, finds the rest.
        (let ((rest-b (list (- (first b) (aref x 0)) (- (second b) (aref x 1)))))
          (multiple-value-bind (rest-code rest-x) (ls-range-run '((1d0 0d0) (0d0 1d0)) rest-b 0d0 0d0 0d0 10)
            (assert (and (= 0 rest-code) (= 0d0 (aref rest-x 0)) (= (second b) (aref rest-x 1))) ()
                    "b = ~S: the second call's outcome is ~D and its vector ~S." b rest-code rest-x)))))
    ;; A column far below the other: the answer's second entry is 1/s, and a sum of its squares
    ;; that nothing reads without damping is above the largest number.
    (dolist (power '(-511 -512 -520 -600))
      (let ((s (two power)))
        (multiple-value-bind (code x workspace residual)
            (ls-range-run (list (list 1d0 0d0) (list 1d0 0d0) (list 0d0 s) (list 0d0 s)) '(1d0 1d0 1d0 1d0) 0d0 1d-12 0d0 20)
          (assert (and (= 0 code) (<= (abs (- (aref x 0) 1d0)) 1d-12) (<= (abs (- (* (aref x 1) s) 1d0)) 1d-12)
                       (<= residual (expt (rational 1d-12) 2)) (<= (lsqr:residual-norm workspace) 1d-12))
                  () "A column of 2^~D: outcome ~D, x = ~S, a residual of ~S." power code x (lsqr:residual-norm workspace)))))
    ;; The same under a damping far below the small column: the iterate's norm, which the damped
    ;; residual's estimate reads, has a square above the largest number.
    (let ((s (two -520)))
      (multiple-value-bind (code x workspace residual)
          (ls-range-run (list (list 1d0 0d0) (list 1d0 0d0) (list 0d0 s) (list 0d0 s)) '(1d0 1d0 1d0 1d0) (two -600) 1d-12 0d0 20)
        (assert (and (= 0 code) (<= (abs (- (aref x 0) 1d0)) 1d-12) (<= (abs (- (* (aref x 1) s) 1d0)) 1d-12)
                     (<= residual (expt (rational 1d-12) 2)))
                () "A column of 2^-520 under a damping of 2^-600: outcome ~D, x = ~S, a residual of ~S." code x
                (lsqr:residual-norm workspace))))
    ;; A row far above the other: the second left vector has an entry of L before it is made a unit.
    (dolist (power '(511 512 600))
      (let ((big (two power)))
        (multiple-value-bind (code x workspace residual)
            (ls-range-run (list (list big big 0d0 0d0) (list 0d0 0d0 1d0 1d0)) (list (/ (two 400) big) (two 400)) 0d0
                          (two 360) 0d0 20)
          (declare (ignorable x))
          (assert (and (= 0 code) (<= residual (expt (rational (two 360)) 2)) (<= (lsqr:residual-norm workspace) (two 360))) ()
                  "A row of 2^~D: outcome ~D, a residual of ~S." power code (lsqr:residual-norm workspace)))))
    ;; A right-hand side all but orthogonal to the range: the transpose's first image is 2^-512 of
    ;; the operator, and the second left vector 2^512 before it is made a unit.
    (multiple-value-bind (code x workspace residual gradient)
        (ls-range-run '((1d0) (0d0)) (list (two -512) 1d0) 0d0 0d0 (two -572) 10)
      (assert (and (= 1 code) (= (two -512) (aref x 0)) (= 1 residual) (= 0 gradient)
                   (= 1d0 (lsqr:residual-norm workspace)) (= 0d0 (lsqr:gradient-norm workspace)))
              () "b = (2^-512, 1) on the first unit vector: outcome ~D, x = ~S, reports ~S and ~S." code x
              (lsqr:residual-norm workspace) (lsqr:gradient-norm workspace)))
    ;; The two scale factors of the gradient are of opposite size, and their product is 1: the
    ;; gradient of the first iterate is 2^-601, as it is for the same problem at scale 1, and is
    ;; above a tolerance of 2^-700.
    (flet ((at (matrix-power)
             (let ((entry (two matrix-power)) (b (make-list 4 :initial-element (two (- matrix-power)))))
               (multiple-value-bind (code x workspace residual gradient)
                   (ls-range-run (list (list entry) (list entry) (list entry) (list (- entry))) b
                                 (two (- matrix-power 300)) 0d0 (two -700) 1)
                 (declare (ignore residual))
                 (assert (and (= 3 code) (= (expt 2 -1202) gradient) (= (two -601) (lsqr:gradient-norm workspace))
                              (= (two (1- (* -2 matrix-power))) (aref x 0)))
                         () "The operator at 2^~D: outcome ~D, x = ~S, a gradient of ~S whose exact square is ~S."
                         matrix-power code x (lsqr:gradient-norm workspace) (float gradient 1d0))
                 (aref (lsqr:recurrence-gradient-norms workspace) 1)))))
      (let ((plain (at 0)) (low (at -500)) (high (at 400)))
        (assert (and (plusp plain) (= plain low high)) ()
                "The gradient's estimate after one step is ~S at scale 1, ~S and ~S elsewhere." plain low high)))
    ;; The gradient's running estimates late in a run are far below 1, and with the operator at
    ;; 2^-500 and the right-hand side at 2^500 one factor alone takes them below every number. They
    ;; are the estimates at scale 1, every one, and so are the confirmations they trigger.
    (flet ((history (matrix-power)
             (let ((rows (mapcar (lambda (row) (mapcar (lambda (value) (* value (two matrix-power))) row))
                                 '((4d0 1d0 0d0) (1d0 3d0 1d0) (0d0 1d0 2d0) (1d0 0d0 1d0))))
                   (b (mapcar (lambda (value) (* value (two (- matrix-power)))) '(1d0 2d0 3d0 4d0))))
               (multiple-value-bind (code x workspace) (ls-range-run rows b 0d0 0d0 0d0 40)
                 (declare (ignore x))
                 (list code (lsqr:step-count workspace) (lsqr:confirmation-count workspace)
                       (copy-seq (subseq (lsqr:recurrence-gradient-norms workspace) 0 (1+ (lsqr:step-count workspace)))))))))
      (let ((plain (history 0)) (low (history -500)) (high (history 400)))
        (assert (and (equalp plain low) (equalp plain high) (some #'plusp (fourth plain))
                     (< 0d0 (reduce #'min (remove 0d0 (fourth plain))) (two -530)))
                () "The gradient's estimates at three scales of one problem: ~S, ~S and ~S." plain low high)))
    ;; A right-hand side the transpose sends to zero is solved by zero whatever the damping, at any scale.
    (dolist (power '(0 500 -500))
      (multiple-value-bind (code x workspace residual gradient)
          (ls-range-run (list (list (two power)) (list 0d0)) (list 0d0 (* 3d0 (two power))) (* 8d0 (two power)) 0d0 0d0 10)
        (declare (ignore residual))
        (assert (and (= 1 code) (= 0d0 (aref x 0)) (= 0 gradient) (= (* 3d0 (two power)) (lsqr:residual-norm workspace))) ()
                "A right-hand side orthogonal to the range at 2^~D: outcome ~D, x = ~S, a residual of ~S." power code x
                (lsqr:residual-norm workspace))))
    ;; The operator norm estimate is not a bound once a singular value has converged, which a
    ;; dominant one does within a few steps: here it is above the Frobenius norm before the sixth
    ;; step of forty.
    (let* ((n 40) (values (loop :for i :below n :collect (if (= i 0) 1000d0 (+ 1d0 (/ (float (1- i) 1d0) 38d0)))))
           (rows (loop :for i :below n :collect (loop :for j :below n :collect (if (= i j) (elt values i) 0d0))))
           (b (loop :for i :below n :collect (+ 1d0 (* 0.37d0 (sin (float (* 7 (1+ i)) 1d0))))))
           (frobenius (sqrt (loop :for value :in values :sum (* value value)))))
      (multiple-value-bind (code x workspace) (ls-range-run rows b 0d0 0d0 0d0 6)
        (declare (ignore x))
        (assert (and (= 3 code) (= 6 (lsqr:step-count workspace)) (> (lsqr:operator-norm-estimate workspace) (* 1.2d0 frobenius)))
                () "One dominant singular value: after ~D steps the estimate is ~S and the Frobenius norm ~S."
                (lsqr:step-count workspace) (lsqr:operator-norm-estimate workspace) frobenius)))
    ;; Both tolerances at once: the status names one that the exact norm meets, the residual's first.
    (dolist (case (list (list '((4d0 1d0 0d0) (1d0 3d0 1d0) (0d0 1d0 2d0) (1d0 0d0 1d0)) '(1d0 2d0 3d0 4d0) 0d0 1d-9 1d-9 1)
                        (list '((4d0 1d0 0d0) (1d0 3d0 1d0) (0d0 1d0 2d0)) '(5d0 5d0 3d0) 0d0 1d-9 1d-9 0)
                        (list '((4d0 1d0 0d0) (1d0 3d0 1d0) (0d0 1d0 2d0) (1d0 0d0 1d0)) '(1d0 2d0 3d0 4d0) 0.5d0 1d-9 1d-9 1)))
      (destructuring-bind (rows b damping residual-tolerance gradient-tolerance expected) case
        (multiple-value-bind (code x workspace residual gradient)
            (ls-range-run rows b damping residual-tolerance gradient-tolerance 60)
          (declare (ignore x workspace))
          (assert (and (= expected code)
                       (if (= code 0)
                           (<= residual (expt (rational residual-tolerance) 2))
                           (and (> residual (expt (rational residual-tolerance) 2))
                                (<= gradient (expt (rational gradient-tolerance) 2)))))
                  () "Both tolerances, damping ~S: outcome ~D with exact squares ~S and ~S." damping code
                  (float residual 1d0) (float gradient 1d0))))))
  (format t "LSQR-RANGE-PASS~%"))

(defun run-lsqr-lifetime-tests ()
  (let ((state (sb-ext:seed-random-state 80)))
    (multiple-value-bind (a b) (ls-case state 5 5 5)
      (let* ((matrix (ls-csc a)) (forward (rf-caller:matrix-operator matrix))
             (backward (ls-caller:transpose-operator matrix))
             (tolerance (* (scale-float 1d0 -30) (ls-frobenius a) (two-norm b)))
             (workspace (lsqr:make-lsqr 5 5 60)) (first-answer (zeros 5)) (x (zeros 5)) (shared (copy-seq b)))
        ;; One workspace, call after call, damped and not: the same answers bit for bit.
        (assert (= 1 (ls-caller:lsqr-code workspace forward backward b first-answer 0d0 0d0 tolerance)))
        (assert (= 1 (ls-caller:lsqr-code workspace forward backward b x 2d0 0d0 tolerance)))
        (assert (= 1 (ls-caller:lsqr-code workspace forward backward b x 0d0 0d0 tolerance)))
        (assert (same-bits-p x first-answer) () "The same call on a used workspace gave another vector.")
        ;; The right-hand side and the solution as one array.
        (assert (= 1 (ls-caller:lsqr-code workspace forward backward shared shared 0d0 0d0 tolerance)))
        (assert (same-bits-p shared first-answer) () "With one array for both, ~S; with two, ~S." shared first-answer)
        ;; Steady state: a compiled caller allocates nothing.
        (allocation-check "lsqr-loop"
                          (lambda ()
                            (assert (= 0 (ls-caller:lsqr-loop workspace forward backward b x 0.5d0 0d0 tolerance 1))))))))
  (multiple-value-bind (accepted diagnostics)
      (coalton-linalg/tests/encapsulation::compile-probe "lsqr-public-control" *lsqr-public-control*)
    (assert accepted () "The valid least squares external caller failed:~%~A" diagnostics))
  (dolist (probe *lsqr-encapsulation-probes*)
    (assert (not (coalton-linalg/tests/encapsulation::compile-probe (car probe) (cdr probe))) ()
            "External caller ~A compiled despite required rejection." (car probe)))
  (dolist (name '("MKLSQRWORKSPACE" "LSQRSTORAGE"))
    (multiple-value-bind (symbol visibility) (find-symbol name "COALTON-LINALG/LSQR")
      (assert (and symbol (eq visibility ':internal)) ()
              "Storage constructor ~A must exist and remain internal; visibility is ~S." name visibility)))
  (format t "LSQR-LIFETIME-PASS rejected=~D public-control=1~%" (length *lsqr-encapsulation-probes*)))

(defun run-lsqr-tests ()
  (run-lsqr-exact-tests)
  (run-lsqr-outcome-tests)
  (run-lsqr-scale-tests)
  (run-lsqr-range-tests)
  (run-lsqr-lifetime-tests)
  (format t "LSQR-TESTS-PASS~%"))
