(in-package #:coalton-linalg/tests)

(defun run-runtime-tests ()
  (let* ((s (coalton-linalg/sparse/spd-ldlt:analyse
             (coalton-linalg/sparse/spd-ldlt:pattern 2 (indices 0 2 3) (indices 0 1 1))))
         (a (doubles 4d0 2d0 3d0)))
    (sb-ext:gc :full t)
    (coalton-linalg/tests/callers:dynamic-shifts s a 10000)
    (let* ((before (sb-ext:get-bytes-consed))
           (status (coalton-linalg/tests/callers:dynamic-shifts s a 10000))
           (bytes (- (sb-ext:get-bytes-consed) before)))
      (format t "DYNAMIC-SHIFT status=~D calls=10000 bytes=~D~%" status bytes)
      (assert (zerop status))
      (assert (zerop bytes))))
  (let* ((a (doubles))
         (s (coalton-linalg/sparse/spd-ldlt:analyse
              (coalton-linalg/sparse/spd-ldlt:pattern 0 (indices 0) (indices))))
         (r (coalton-linalg/sparse/spd-ldlt:factor! s a 0d0)))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:FactorSuccess) (coalton-linalg/sparse/spd-ldlt:factor-status r)))
    (assert (eq (coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess) (coalton-linalg/sparse/spd-ldlt:solve!
                    (coalton-linalg/sparse/spd-ldlt:result-factor r) a a)))
    (format t "EMPTY-SYSTEM-PASS~%")))