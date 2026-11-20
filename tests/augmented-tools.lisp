(in-package #:coalton-linalg/tools)

;;; Regressions for the augmented artifact driver: one small system whose held and plain
;;; factors, growth and solution are exact, a zero pivot named with its column, a typed
;;; refusal file, a system that names no run, and removal of stale artifacts.

(defun expect-signed (output name index first-line numbers last-line)
  "NUMBERS are the refused column, the two counts and the growth. The backward bound follows them: -1 with
no factor, and otherwise positive and tiny, since these factors are exact."
  (let* ((path (fixture-path output name (format nil ".~D.signed.out.txt" index)))
         (lines (file-lines path)) (written (file-numbers path)))
    (unless (and (= 7 (length lines)) (string= first-line (first lines)) (string= last-line (car (last lines)))
                 (= 5 (length written)) (equalp numbers (subseq written 0 4))
                 (if (string= last-line "solve none")
                     (= -1d0 (fifth written))
                     (< 0d0 (fifth written) 1d-12)))
      (error "Augmented driver wrote ~S for run ~D of ~A." lines index name))))

(defun run-augmented-tools-tests ()
  (let* ((base (asdf:system-source-directory "coalton-linalg"))
         (dir (merge-pathnames "build/augmented-runner-input/" base))
         (output (merge-pathnames "build/augmented-runner-output/" base)))
    (ensure-directories-exist dir)
    (dolist (stale (directory (merge-pathnames "*.*" dir))) (delete-file stale))
    ;; A = [[1], [2]], so K = [[top, 1, 2], [1, -bottom, 0], [2, 0, -bottom]] and with b = (6, 6)
    ;; and top = bottom = 1 the solution is (3, -3, 0) exactly, in either order.
    (write-text (merge-pathnames "tall.J.txt" dir) (format nil "2 1 2~%0 0 1~%1 0 2~%"))
    (write-text (merge-pathnames "tall.b.txt" dir) (format nil "6~%6~%"))
    (write-text (merge-pathnames "tall.meta.txt" dir)
                (format nil "name: tall~%augmented: 1 1 held~%augmented: 1 1 plain~%augmented: 1 0 held~%"))
    (write-text (merge-pathnames "twice.J.txt" dir) (format nil "1 1 2~%0 0 1~%0 0 1~%"))
    (write-text (merge-pathnames "twice.b.txt" dir) (format nil "1~%"))
    (write-text (merge-pathnames "twice.meta.txt" dir) (format nil "augmented: 1 1 held~%"))
    ;; A system of some other driver's: no `augmented:` line, so no run and a line that says so.
    (write-text (merge-pathnames "quiet.J.txt" dir) (format nil "1 1 1~%0 0 1~%"))
    (write-text (merge-pathnames "quiet.meta.txt" dir) (format nil "gram: Jt J~%"))
    (run-augmented-fixtures dir output :repeats 1 :allocation-repeats 1000)
    ;; Held: the column first, then each row. Pivots 1, -2, -3; |L||D||L^T| has diagonal 1, 3, 9
    ;; and the largest magnitude factored is 2, so growth is 4.5.
    (expect-numbers (fixture-path output "tall" ".0.K.out.txt")
                    '(3 5 0 0 1d0 1 0 1d0 2 0 2d0 1 1 -1d0 2 2 -1d0))
    (expect-numbers (fixture-path output "tall" ".0.perm.out.txt") '(0 1 2))
    (expect-numbers (fixture-path output "tall" ".0.D.out.txt") '(1d0 -2d0 -3d0))
    (expect-numbers (fixture-path output "tall" ".0.L.out.txt") '(3 3 1 0 1d0 2 0 2d0 2 1 1d0))
    (expect-numbers (fixture-path output "tall" ".0.x.out.txt") '(3d0 -3d0 0d0))
    (expect-signed output "tall" 0 "status SignedSuccess" '(-1 1 2 4.5d0) "solve SolveSuccess")
    ;; Plain minimum degree takes a row first, then the column, then the other row: pivots
    ;; -1, 2, -3, the diagonal of |L||D||L^T| is 1, 3, 5 and growth is 5/2.
    (expect-numbers (fixture-path output "tall" ".1.perm.out.txt") '(1 0 2))
    (expect-numbers (fixture-path output "tall" ".1.D.out.txt") '(-1d0 2d0 -3d0))
    (expect-numbers (fixture-path output "tall" ".1.L.out.txt") '(3 2 1 0 -1d0 2 1 1d0))
    (expect-numbers (fixture-path output "tall" ".1.x.out.txt") '(3d0 -3d0 0d0))
    (expect-signed output "tall" 1 "status SignedSuccess" '(-1 1 2 2.5d0) "solve SolveSuccess")
    ;; With a zero trailing block the two rows are dependent after the column: an exact zero pivot.
    (expect-numbers (fixture-path output "tall" ".2.K.out.txt")
                    '(3 5 0 0 1d0 1 0 1d0 2 0 2d0 1 1 0d0 2 2 0d0))
    (expect-signed output "tall" 2 "status SignedZeroPivot" '(2 0 0 -1d0) "solve none")
    (dolist (suffix '(".2.perm.out.txt" ".2.D.out.txt" ".2.L.out.txt" ".2.x.out.txt"))
      (assert (not (probe-file (fixture-path output "tall" suffix)))))
    (unless (equal '("DuplicateEntry") (file-lines (fixture-path output "twice" ".refused.out.txt")))
      (error "Augmented driver did not name the duplicate-entry refusal."))
    (assert (not (probe-file (fixture-path output "twice" ".0.K.out.txt"))))
    (assert (not (probe-file (fixture-path output "quiet" ".0.K.out.txt"))))
    (let ((lines (file-lines (merge-pathnames "augmented-table.txt" output))))
      (unless (and (= 6 (length lines)) (find "quiet skipped" lines :test #'string=)
                   (find "twice refused DuplicateEntry" lines :test #'string=))
        (error "The augmented table holds ~S." lines)))
    ;; A rerun with fewer runs removes the artifacts of the runs that no longer exist, and a
    ;; corrected system replaces its refusal.
    (write-text (merge-pathnames "tall.meta.txt" dir) (format nil "augmented: 1 1 plain~%"))
    (write-text (merge-pathnames "twice.J.txt" dir) (format nil "1 1 1~%0 0 3~%"))
    (run-augmented-fixtures dir output :repeats 1 :allocation-repeats 1000)
    (expect-numbers (fixture-path output "tall" ".0.perm.out.txt") '(1 0 2))
    (dolist (suffix '(".1.K.out.txt" ".1.x.out.txt" ".1.signed.out.txt" ".2.K.out.txt" ".2.signed.out.txt"))
      (assert (not (probe-file (fixture-path output "tall" suffix)))))
    (assert (not (probe-file (fixture-path output "twice" ".refused.out.txt"))))
    ;; [[1, 3], [3, -1]] with right-hand side (0, 1): x = (0.3, -0.1).
    (let ((x (file-numbers (fixture-path output "twice" ".0.x.out.txt"))))
      (unless (and (= 2 (length x)) (< (abs (- (first x) 0.3d0)) 1d-15) (< (abs (+ (second x) 0.1d0)) 1d-15))
        (error "Augmented driver wrote ~S for the corrected system." x)))
    ;; A malformed run line is a driver error, never a silent skip.
    (write-text (merge-pathnames "tall.meta.txt" dir) (format nil "augmented: 1 held~%"))
    (unless (handler-case (progn (run-augmented-fixtures dir output :repeats 1 :allocation-repeats 1000) nil)
              (error () t))
      (error "Augmented driver accepted a run line with no trailing diagonal."))
    ;; So is a system that names a run and has no right-hand side.
    (write-text (merge-pathnames "tall.meta.txt" dir) (format nil "augmented: 1 1 held~%"))
    (delete-file (merge-pathnames "tall.b.txt" dir))
    (unless (handler-case (progn (run-augmented-fixtures dir output :repeats 1 :allocation-repeats 1000) nil)
              (error () t))
      (error "Augmented driver accepted a system with no right-hand side."))
    (format t "AUGMENTED-RUNNER-TESTS-PASS~%")))
