;;;; Register the pinned, already installed sources in a fresh process.
;;;; The check runner supplies absolute roots; see docs/setup.md.
(require :asdf)

(defpackage #:coalton-linalg/setup
  (:use #:cl))
(in-package #:coalton-linalg/setup)

(defun required-directory (variable)
  "Resolve an absolute directory supplied by the standalone check runner."
  (let ((value (uiop:getenv variable)))
    (unless (and value (plusp (length value)))
      (error "~A must name an installed directory." variable))
    (let ((path (uiop:ensure-directory-pathname
                 (uiop:parse-native-namestring value))))
      (unless (and (uiop:absolute-pathname-p path)
                   (uiop:directory-exists-p path))
        (error "~A must name an existing absolute directory: ~A" variable value))
      (truename path))))

(defun require-file (root relative)
  "Refuse incomplete installations before Quicklisp can initialize a dist."
  (or (probe-file (merge-pathnames relative root))
      (error "Missing installed dependency file: ~A"
             (merge-pathnames relative root))))

(defun check-output-mapping (source-root cache-root namespace)
  "Check that a nested output retains its relative path in the named cache scope."
  (let* ((probe (merge-pathnames "cache-probe/deep/file.fasl" source-root))
         (actual (asdf:apply-output-translations probe))
         (expected (merge-pathnames
                    (format nil "*/~A/cache-probe/deep/file.fasl" namespace)
                    cache-root)))
    (unless (pathname-match-p actual expected)
      (error "ASDF cache mapping for ~A escaped its relative namespace: ~A"
             namespace actual))))
(defun register-dependencies (library-root)
  "Isolate dependency discovery and compiled output for this disposable process.
No downloading or Quicklisp quickloading is performed. Optional local Quicklisp
startup hooks are refused. The ordinary library ASDF system does not call this."
  (when (or (find-package "QUICKLISP") (find-package "COALTON"))
    (error "Dependency registration requires a fresh Lisp process without Quicklisp or Coalton."))
  (unless (equal (uiop:getenv "COALTON_ENV") "release")
    (error "Set COALTON_ENV=release before starting this process."))
  (unless (asdf:version-satisfies (asdf:asdf-version) "3.3.1")
    (error "The standalone recipe requires the supported SBCL's ASDF 3.3.1 or newer."))
  (let* ((coalton-root (required-directory "COALTON_LINALG_COALTON_ROOT"))
         (quicklisp-root (required-directory "COALTON_LINALG_QUICKLISP_ROOT"))
         (cache-root (required-directory "XDG_CACHE_HOME"))
         (setup (require-file quicklisp-root "setup.lisp"))
         (search-functions (copy-list asdf:*system-definition-search-functions*)))
    (require-file coalton-root "coalton.asd")
    (dolist (relative '("quicklisp/quicklisp.asd"
                        "dists/quicklisp/distinfo.txt"
                        "dists/quicklisp/releases.txt"
                        "dists/quicklisp/systems.txt"))
      (require-file quicklisp-root relative))
    (when (or (probe-file (merge-pathnames "asdf-config/init.lisp" quicklisp-root))
              (directory (merge-pathnames "local-init/*.lisp" quicklisp-root))
              (directory (merge-pathnames "local-init/*.cl" quicklisp-root)))
      (error "Use an isolated Quicklisp root without asdf-config/init.lisp or local-init hooks."))
    (setf (get :coalton-config :compiler-mode) "release"
          (get :coalton-config :perform-specialization) t
          (get :coalton-config :emit-type-annotations) t
          asdf:*central-registry* nil)
    (asdf:initialize-source-registry
     `(:source-registry
       (:directory ,library-root)
       (:tree ,coalton-root)
       (:directory ,(merge-pathnames "quicklisp/" quicklisp-root))
       (:tree ,(merge-pathnames "dists/quicklisp/software/" quicklisp-root))
       :ignore-inherited-configuration))
    (let ((mappings (list (cons coalton-root "coalton")
                          (cons quicklisp-root "quicklisp")
                          (cons library-root "library"))))
      (asdf:initialize-output-translations
       `(:output-translations
         ,@(loop for (root . namespace) in mappings
                 collect `(,(namestring root)
                           (,cache-root :implementation ,namespace)))
         (t (,cache-root :implementation "other"))
         :ignore-inherited-configuration))
      (dolist (mapping mappings)
        (check-output-mapping (car mapping) cache-root (cdr mapping))))
    (unwind-protect
         (load setup :verbose nil :print nil)
      (setf asdf:*system-definition-search-functions* search-functions
            asdf:*central-registry* nil))
    (set (or (find-symbol "*LOCAL-PROJECT-DIRECTORIES*" "QUICKLISP-CLIENT")
             (error "The installed Quicklisp client lacks its local-project setting."))
         nil)
    (format t "DEPENDENCIES registered explicit Coalton/Quicklisp roots; inherited ASDF configuration disabled.~%")))

(register-dependencies
 (uiop:pathname-parent-directory-pathname
  (uiop:pathname-directory-pathname *load-truename*)))