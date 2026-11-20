(in-package #:coalton-linalg/tools)

;;; Regressions for the orthogonal artifact driver: an exact QR file, runs by threshold and by
;;; forced rank, a truncated solve with its nullspace vector, named rank and factor refusals,
;;; and removal of an earlier run's indexed artifacts.

(defun run-orthogonal-tools-tests ()
  (let* ((base (asdf:system-source-directory "coalton-linalg"))
         (dir (merge-pathnames "build/orthogonal-runner-input/" base))
         (output (merge-pathnames "build/orthogonal-runner-output/" base)))
    (ensure-directories-exist dir)
    (dolist (stale (directory (merge-pathnames "*.*" dir))) (delete-file stale))
    ;; diag(1, 2): the second column is the pivot, Q exchanges the two rows, R = diag(2, 1).
    (write-text (merge-pathnames "diagonal.M.txt" dir) (format nil "2 2~%1~%0~%0~%2~%"))
    (write-text (merge-pathnames "diagonal.b.txt" dir) (format nil "3~%8~%"))
    (write-text (merge-pathnames "diagonal.meta.txt" dir)
                (format nil "name: diagonal~%threshold: 1e-10~%rank: 1~%rank: 3~%"))
    ;; A zero row vector with no meta file: one run at the default threshold, rank 0.
    (write-text (merge-pathnames "nothing.M.txt" dir) (format nil "1 2~%0~%0~%"))
    (write-text (merge-pathnames "nothing.b.txt" dir) (format nil "5~%"))
    (write-text (merge-pathnames "forced.M.txt" dir) (format nil "1 2~%0~%0~%"))
    (write-text (merge-pathnames "forced.b.txt" dir) (format nil "5~%"))
    (write-text (merge-pathnames "forced.meta.txt" dir) (format nil "rank: 1~%"))
    ;; A column whose norm overflows: refused at the factorization, no QR file.
    (write-text (merge-pathnames "huge.M.txt" dir) (format nil "2 1~%1.5d308~%1.5d308~%"))
    (write-text (merge-pathnames "huge.b.txt" dir) (format nil "1~%1~%"))
    (run-orthogonal-fixtures dir output :repeats 1 :allocation-repeats 10)
    ;; m, n, permutation, R column-major (2 by 2), thin Q column-major (2 by 2).
    (expect-numbers (fixture-path output "diagonal" ".qr.out.txt")
                    '(2 2 1 0 2d0 0d0 0d0 1d0 0d0 1d0 1d0 0d0))
    (expect-numbers (fixture-path output "diagonal" ".0.x.out.txt") '(3d0 4d0))
    (expect-numbers (fixture-path output "diagonal" ".0.null.out.txt") '())
    (expect-numbers (fixture-path output "diagonal" ".0.run.out.txt") '(2 0d0))
    ;; Forced rank 1 keeps the pivot column only: x = [0, 4], residual 3, nullspace e0.
    (expect-numbers (fixture-path output "diagonal" ".1.x.out.txt") '(0d0 4d0))
    (expect-numbers (fixture-path output "diagonal" ".1.null.out.txt") '(1d0 0d0))
    (expect-numbers (fixture-path output "diagonal" ".1.run.out.txt") '(1 3d0))
    (unless (equal '("RankTooLarge") (file-lines (fixture-path output "diagonal" ".2.refused.out.txt")))
      (error "Orthogonal driver did not name the rank refusal."))
    (assert (not (probe-file (fixture-path output "diagonal" ".2.x.out.txt"))))
    (expect-numbers (fixture-path output "nothing" ".qr.out.txt") '(1 2 0 1 0d0 0d0 1d0))
    (expect-numbers (fixture-path output "nothing" ".0.x.out.txt") '(0d0 0d0))
    (expect-numbers (fixture-path output "nothing" ".0.null.out.txt") '(1d0 0d0 0d0 1d0))
    (expect-numbers (fixture-path output "nothing" ".0.run.out.txt") '(0 5d0))
    (unless (equal '("RankSingular") (file-lines (fixture-path output "forced" ".0.refused.out.txt")))
      (error "Orthogonal driver did not name the singular rank."))
    (unless (equal '("FactorNonFinite") (file-lines (fixture-path output "huge" ".refused.out.txt")))
      (error "Orthogonal driver did not name the factor refusal."))
    (assert (not (probe-file (fixture-path output "huge" ".qr.out.txt"))))
    ;; One table line a run, or one a refused factorization, after the header.
    (let ((lines (file-lines (merge-pathnames "orthogonal-table.txt" output))))
      (unless (= 7 (length lines)) (error "Orthogonal table has ~D lines; expected 7." (length lines))))
    ;; A rerun with fewer runs removes the artifacts of the runs that no longer exist.
    (write-text (merge-pathnames "diagonal.meta.txt" dir) (format nil "rank: 1~%"))
    (run-orthogonal-fixtures dir output :repeats 1 :allocation-repeats 10)
    (expect-numbers (fixture-path output "diagonal" ".0.x.out.txt") '(0d0 4d0))
    (assert (not (probe-file (fixture-path output "diagonal" ".1.x.out.txt"))))
    (assert (not (probe-file (fixture-path output "diagonal" ".1.run.out.txt"))))
    (assert (not (probe-file (fixture-path output "diagonal" ".2.refused.out.txt"))))
    ;; A malformed system is a driver error, never a silent skip.
    (write-text (merge-pathnames "short.M.txt" dir) (format nil "2 2~%1~%2~%3~%"))
    (write-text (merge-pathnames "short.b.txt" dir) (format nil "1~%1~%"))
    (unless (handler-case (progn (run-orthogonal-fixtures dir output :repeats 1 :allocation-repeats 10) nil)
              (error () t))
      (error "Orthogonal driver accepted a truncated matrix."))
    (delete-file (merge-pathnames "short.M.txt" dir))
    (format t "ORTHOGONAL-RUNNER-TESTS-PASS~%")))
