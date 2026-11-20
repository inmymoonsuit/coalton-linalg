(in-package #:coalton-linalg/tools)

;;; Regressions for the Jacobi artifact driver: exact eigenvalue and singular value files, runs
;;; by sweep limit with a named non-convergence, strict rank lines by threshold, a named
;;; overflow refusal, a name used for both kinds refused, and removal of stale artifacts.

(defun run-jacobi-tools-tests ()
  (let* ((base (asdf:system-source-directory "coalton-linalg"))
         (dir (merge-pathnames "build/jacobi-runner-input/" base))
         (output (merge-pathnames "build/jacobi-runner-output/" base)))
    (ensure-directories-exist dir)
    (dolist (stale (directory (merge-pathnames "*.*" dir))) (delete-file stale))
    ;; diag(3, 1): already diagonal, so one sweep; the values come back ascending with their vectors.
    (write-text (merge-pathnames "plain.S.txt" dir) (format nil "2~%3~%0~%0~%1~%"))
    ;; [[2, 1], [1, 2]] needs a rotation: a limit of 0 or 1 cannot finish, a limit of 5 can.
    (write-text (merge-pathnames "coupled.S.txt" dir) (format nil "2~%2~%1~%77~%2~%"))
    (write-text (merge-pathnames "coupled.meta.txt" dir)
                (format nil "name: coupled~%sweep limit: 0~%sweep limit: 1~%sweep limit: 5~%"))
    ;; diag(3, 4) with two thresholds: 0.5 of the largest keeps both, 0.75 is exactly 3 and keeps one.
    (write-text (merge-pathnames "square.M.txt" dir) (format nil "2 2~%3~%0~%0~%4~%"))
    (write-text (merge-pathnames "square.meta.txt" dir) (format nil "threshold: 0.5~%threshold: 0.75~%"))
    (write-text (merge-pathnames "huge.M.txt" dir) (format nil "2 1~%1.5d308~%1.5d308~%"))
    (run-jacobi-fixtures dir output :repeats 1 :allocation-repeats 1)
    ;; n, values, vectors column-major, sweeps.
    (expect-numbers (fixture-path output "plain" ".0.eigen.out.txt") '(2 1d0 3d0 0d0 1d0 1d0 0d0 1))
    (dolist (index '(0 1))
      (unless (equal '("JacobiNotConverged")
                     (file-lines (fixture-path output "coupled" (format nil ".~D.refused.out.txt" index))))
        (error "Jacobi driver did not name the sweep-limit refusal of run ~D." index))
      (assert (not (probe-file (fixture-path output "coupled" (format nil ".~D.eigen.out.txt" index))))))
    (let ((numbers (file-numbers (fixture-path output "coupled" ".2.eigen.out.txt"))))
      (unless (and (= 8 (length numbers)) (= 2 (first numbers))
                   (< (abs (- (second numbers) 1d0)) 1d-15) (< (abs (- (third numbers) 3d0)) 1d-15)
                   (= 2 (car (last numbers))))
        (error "Jacobi driver wrote ~S for the coupled matrix." numbers)))
    ;; m, n, values, left vectors, right vectors, sweeps, then one rank a threshold.
    (expect-numbers (fixture-path output "square" ".0.svd.out.txt")
                    '(2 2 4d0 3d0 0d0 1d0 1d0 0d0 0d0 1d0 1d0 0d0 1 2 1))
    (unless (equal '("JacobiNonFinite") (file-lines (fixture-path output "huge" ".0.refused.out.txt")))
      (error "Jacobi driver did not name the overflow refusal."))
    (unless (= 5 (length (file-lines (merge-pathnames "eigen-table.txt" output))))
      (error "The eigenvalue table does not have a header and four runs."))
    (unless (= 3 (length (file-lines (merge-pathnames "svd-table.txt" output))))
      (error "The singular value table does not have a header and two runs."))
    ;; A rerun with fewer runs removes the artifacts of the runs that no longer exist.
    (write-text (merge-pathnames "coupled.meta.txt" dir) (format nil "sweep limit: 5~%"))
    (run-jacobi-fixtures dir output :repeats 1 :allocation-repeats 1)
    (assert (probe-file (fixture-path output "coupled" ".0.eigen.out.txt")))
    (assert (not (probe-file (fixture-path output "coupled" ".0.refused.out.txt"))))
    (assert (not (probe-file (fixture-path output "coupled" ".1.refused.out.txt"))))
    (assert (not (probe-file (fixture-path output "coupled" ".2.eigen.out.txt"))))
    ;; One name for both kinds would write the same refusal file twice: a driver error.
    (write-text (merge-pathnames "plain.M.txt" dir) (format nil "1 1~%1~%"))
    (unless (handler-case (progn (run-jacobi-fixtures dir output :repeats 1 :allocation-repeats 1) nil)
              (error () t))
      (error "Jacobi driver accepted one name for a symmetric and a rectangular system."))
    (delete-file (merge-pathnames "plain.M.txt" dir))
    ;; A malformed system is a driver error, never a silent skip.
    (write-text (merge-pathnames "short.S.txt" dir) (format nil "2~%1~%2~%3~%"))
    (unless (handler-case (progn (run-jacobi-fixtures dir output :repeats 1 :allocation-repeats 1) nil)
              (error () t))
      (error "Jacobi driver accepted a truncated matrix."))
    (delete-file (merge-pathnames "short.S.txt" dir))
    (format t "JACOBI-RUNNER-TESTS-PASS~%")))
