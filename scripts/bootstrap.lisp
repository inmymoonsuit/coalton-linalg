(require :asdf)
(unless (and (member :sbcl *features*) (= most-positive-fixnum (1- (expt 2 62))))
  (error "This initial library targets 64-bit SBCL."))

(defparameter *library-root*
  (uiop:pathname-parent-directory-pathname
    (uiop:pathname-directory-pathname *load-truename*)))

(load (merge-pathnames "scripts/register-dependencies.lisp" *library-root*))

(let ((safety (uiop:getenv "COALTON_LINALG_SAFETY")))
  (unless (member safety '("0" "1") :test #'equal)
    (error "COALTON_LINALG_SAFETY must be 0 or 1."))
  (when (equal safety "0") (pushnew :coalton-linalg-safety-zero *features*)))

(asdf:load-system "coalton")
(asdf:load-asd (merge-pathnames "coalton-linalg.asd" *library-root*))
;; Safety is a reader feature: never reuse a library FASL from the other policy.
;; A warning from the compiler while the library is built is a failure here. The build would
;; otherwise go on: a function inlined from another package was once compiled with a variable
;; left unbound, reported only as a warning, and returned wrong numbers. The compiler handles
;; its own warnings, so they are recognized in what it prints.
(let* ((printed (make-string-output-stream))
       (*error-output* (make-broadcast-stream *error-output* printed)))
  (asdf:load-system "coalton-linalg" :force t)
  (when (search "caught WARNING" (get-output-stream-string printed))
    (error "The compiler warned while building coalton-linalg; see its output above.")))
(format t "TOOLCHAIN SBCL=~A COALTON_ENV=~A safety=~A~%"
        (lisp-implementation-version) (uiop:getenv "COALTON_ENV")
        (uiop:getenv "COALTON_LINALG_SAFETY"))