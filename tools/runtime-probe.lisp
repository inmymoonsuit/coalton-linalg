(load (merge-pathnames "../scripts/bootstrap.lisp" *load-truename*))
(load (merge-pathnames "tests/interface.lisp" *library-root*))
(in-package #:coalton-linalg/tests)

(defun probe-idle ()
  (let ((end (+ (get-internal-real-time)
                (ceiling internal-time-units-per-second 20))))
    (loop while (< (get-internal-real-time) end)))
  0)

(defun probe-batch (thunk count)
  (declare (type function thunk) (type fixnum count))
  (let ((started (get-internal-real-time))
        (before (sb-ext:get-bytes-consed)))
    (dotimes (i count) (funcall thunk))
    (values (- (sb-ext:get-bytes-consed) before)
            (- (get-internal-real-time) started))))

(defun run-allocation-probe ()
  (let* ((mode (or (uiop:getenv "COALTON_LINALG_RUNTIME_PROBE") "factor-first"))
         (paused (string= mode "paused"))
         (idle-first (string= mode "idle-first"))
         (s (coalton-linalg/sparse/spd-ldlt:analyse
             (coalton-linalg/sparse/spd-ldlt:pattern
              3 (indices 0 3 4 5) (indices 0 1 2 1 2))))
         (a (doubles 5d0 1d0 1d0 2d0 2d0))
         (factor-thunk (lambda () (coalton-linalg/sparse/spd-ldlt:factor! s a 0d0)))
         (control-thunk (lambda () (setf *allocation-control* (cons 1 2))))
         (bytes (make-array 7 :element-type '(unsigned-byte 64) :initial-element 0))
         (ticks (make-array 7 :element-type '(unsigned-byte 64) :initial-element 0)))
    (unless (member mode '("factor-first" "idle-first" "paused") :test #'string=)
      (error "Unknown allocation probe mode ~S." mode))
    ;; This diagnostic changes only this disposable process's runtime thread.
    (when paused (sb-impl::finalizer-thread-stop))
    (unwind-protect
         (progn
           (sb-ext:gc :full t)
           ;; No warmup and no I/O occur between GC and the first measurement.
           (flet ((record (slot thunk count)
                    (multiple-value-bind (allocated elapsed) (probe-batch thunk count)
                      (setf (aref bytes slot) allocated (aref ticks slot) elapsed))))
             (if idle-first
                 (progn (record 0 #'probe-idle 1) (record 1 factor-thunk 10000))
                 (progn (record 0 factor-thunk 10000) (record 1 #'probe-idle 1)))
             (record 2 factor-thunk 10000)
             (record 3 #'probe-idle 1)
             (record 4 factor-thunk 10000)
             (record 5 control-thunk 10000)
             (record 6 factor-thunk 10000)))
      (when paused (sb-impl::finalizer-thread-start)))
    (format t "RUNTIME-PROBE mode=~A ticks-per-second=~D~%" mode internal-time-units-per-second)
    (dotimes (i 7)
      (format t "slot=~D operation=~A bytes=~D ticks=~D~%"
              i (case i (0 (if idle-first "idle" "factor"))
                        (1 (if idle-first "factor" "idle"))
                        (3 "idle") (5 "allocating-control") (otherwise "factor"))
              (aref bytes i) (aref ticks i)))
    (assert (plusp (aref bytes 5)) () "Allocating control was not detected.")
    (format t "RUNTIME-PROBE-CONTROL-PASS~%")))

(run-allocation-probe)
