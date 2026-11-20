(in-package #:coalton-linalg/tools)

;;; Regressions for the dense indefinite artifact driver: an exact factor file, shifts read from
;;; the meta file, a singular solve named instead of written, a named factor refusal, and
;;; removal of an earlier run's indexed artifacts.

(defun file-lines (path)
  (with-open-file (in path)
    (loop :for line := (read-line in nil nil) :while line :collect line)))

(defun file-numbers (path)
  "Every token of a driver output file read as a number; the `counts` word is skipped."
  (with-open-file (in path)
    (loop :for token := (read-numeric in) :while token
          :unless (symbolp token) :collect token)))

(defun expect-numbers (path expected)
  (let ((actual (file-numbers path)))
    (unless (equalp actual expected)
      (error "Indefinite driver wrote ~S in ~A; expected ~S." actual (file-namestring path) expected))))

(defun run-indefinite-tools-tests ()
  (let* ((base (asdf:system-source-directory "coalton-linalg"))
         (dir (merge-pathnames "build/indefinite-runner-input/" base))
         (output (merge-pathnames "build/indefinite-runner-output/" base)))
    (ensure-directories-exist dir)
    (dolist (stale (directory (merge-pathnames "*.*" dir))) (delete-file stale))
    ;; [[0,1],[1,0]] with no meta file: one run at shift 0, one 2-by-2 block, x = [5,3].
    (write-text (merge-pathnames "exchange.S.txt" dir) (format nil "2~%0~%1~%1~%0~%"))
    (write-text (merge-pathnames "exchange.b.txt" dir) (format nil "3~%5~%"))
    ;; [[1,1],[1,1]] at shifts 0 and 2: singular, then [[-1,1],[1,-1]] which is singular too;
    ;; at shift 1 it is the exchange matrix. The strict upper entry is junk and must not be read.
    (write-text (merge-pathnames "ones.S.txt" dir) (format nil "2~%1~%1~%77~%1~%"))
    (write-text (merge-pathnames "ones.b.txt" dir) (format nil "1~%2~%"))
    (write-text (merge-pathnames "ones.meta.txt" dir) (format nil "name: ones~%shift: 0~%shift: 1~%shift: 2~%"))
    ;; Elimination overflows: a named factor refusal and no factor file.
    (write-text (merge-pathnames "huge.S.txt" dir)
                (format nil "2~%1.4d308~%1.4d308~%0~%-1.4d308~%"))
    (write-text (merge-pathnames "huge.b.txt" dir) (format nil "1~%1~%"))
    (write-text (merge-pathnames "empty.S.txt" dir) (format nil "0~%"))
    (write-text (merge-pathnames "empty.b.txt" dir) "")
    (run-indefinite-fixtures dir output :repeats 1 :allocation-repeats 10)
    ;; n, permutation, block sizes, diagonal, subdiagonal, L column-major, counts.
    (expect-numbers (fixture-path output "exchange" ".0.factor.out.txt")
                    '(2 0 1 2 0 0d0 0d0 1d0 0d0 1d0 0d0 0d0 1d0 1 1 0))
    (expect-numbers (fixture-path output "exchange" ".0.x.out.txt") '(5d0 3d0))
    (expect-numbers (fixture-path output "ones" ".0.factor.out.txt")
                    '(2 0 1 1 1 1d0 0d0 0d0 0d0 1d0 1d0 0d0 1d0 1 0 1))
    (unless (equal '("SolveSingular") (file-lines (fixture-path output "ones" ".0.solve.out.txt")))
      (error "Indefinite driver did not name the singular solve."))
    (assert (not (probe-file (fixture-path output "ones" ".0.x.out.txt"))))
    (expect-numbers (fixture-path output "ones" ".1.factor.out.txt")
                    '(2 0 1 2 0 0d0 0d0 1d0 0d0 1d0 0d0 0d0 1d0 1 1 0))
    (expect-numbers (fixture-path output "ones" ".1.x.out.txt") '(2d0 1d0))
    (unless (equal '("SolveSingular") (file-lines (fixture-path output "ones" ".2.solve.out.txt")))
      (error "Indefinite driver did not name the second singular solve."))
    (unless (equal '("FactorNonFinite") (file-lines (fixture-path output "huge" ".0.refused.out.txt")))
      (error "Indefinite driver did not name the factor refusal."))
    (assert (not (probe-file (fixture-path output "huge" ".0.factor.out.txt"))))
    (expect-numbers (fixture-path output "empty" ".0.factor.out.txt") '(0 0 0 0))
    ;; One table line a run: four systems, six runs, after the header.
    (let ((lines (file-lines (merge-pathnames "indefinite-table.txt" output))))
      (unless (= 7 (length lines)) (error "Indefinite table has ~D lines; expected 7." (length lines))))
    ;; A rerun with fewer shifts removes the artifacts of the runs that no longer exist.
    (write-text (merge-pathnames "ones.meta.txt" dir) (format nil "shift: 1~%"))
    (run-indefinite-fixtures dir output :repeats 1 :allocation-repeats 10)
    (expect-numbers (fixture-path output "ones" ".0.x.out.txt") '(2d0 1d0))
    (assert (not (probe-file (fixture-path output "ones" ".0.solve.out.txt"))))
    (assert (not (probe-file (fixture-path output "ones" ".1.factor.out.txt"))))
    (assert (not (probe-file (fixture-path output "ones" ".2.solve.out.txt"))))
    ;; A malformed system is a driver error, never a silent skip.
    (write-text (merge-pathnames "short.S.txt" dir) (format nil "2~%1~%2~%3~%"))
    (write-text (merge-pathnames "short.b.txt" dir) (format nil "1~%1~%"))
    (unless (handler-case (progn (run-indefinite-fixtures dir output :repeats 1 :allocation-repeats 10) nil)
              (error () t))
      (error "Indefinite driver accepted a truncated matrix."))
    (delete-file (merge-pathnames "short.S.txt" dir))
    (format t "INDEFINITE-RUNNER-TESTS-PASS~%")))
