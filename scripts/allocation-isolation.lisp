(defun call-with-allocation-isolation (thunk)
  "Run bounded checks in this disposable process with its finalizer paused."
  (let ((stop (find-symbol "FINALIZER-THREAD-STOP" "SB-IMPL"))
        (start (find-symbol "FINALIZER-THREAD-START" "SB-IMPL")))
    (unless (and stop start (fboundp stop) (fboundp start))
      (error "The pinned SBCL finalizer hooks are unavailable."))
    (format t "ALLOCATION-ISOLATION finalizer-paused~%")
    (unwind-protect
         (progn (funcall stop) (funcall thunk))
      (funcall start)
      (format t "ALLOCATION-ISOLATION finalizer-restored~%"))))

(defun check-isolation-restoration ()
  "Plant a failure inside the measurement policy and require thread restoration."
  (let ((caught nil))
    (handler-case
        (call-with-allocation-isolation (lambda () (error "PLANTED-ISOLATION-FAILURE")))
      (simple-error (condition)
        (unless (string= (princ-to-string condition) "PLANTED-ISOLATION-FAILURE")
          (error condition))
        (setf caught t)))
    (assert caught)
    (let ((name (find-symbol "*FINALIZER-THREAD*" "SB-IMPL")))
      (unless (and name (boundp name) (symbol-value name)
                   (sb-thread:thread-alive-p (symbol-value name)))
        (error "Finalizer thread was not restored after the planted failure.")))
    (format t "ISOLATION-RESTORATION-CONTROL-PASS~%")))