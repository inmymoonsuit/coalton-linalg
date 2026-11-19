(asdf:defsystem "coalton-linalg"
  :description "Numerical linear algebra in Coalton; currently sparse SPD LDL factorization."
  :version "0.1.0"
  :license "MIT"
  :defsystem-depends-on ("coalton-asdf")
  :depends-on ("coalton")
  :in-order-to ((asdf:test-op (asdf:test-op "coalton-linalg/tests")))
  :around-compile
  (lambda (compile)
    (let ((sb-ext:*derive-function-types* t)
          (sb-ext:*block-compile-default* :specified))
      (funcall compile)))
  :components ((:module "src" :serial t
                :components ((:ct-file "errors") (:ct-file "numeric")
                             (:module "sparse" :serial t
                              :components ((:ct-file "csc") (:ct-file "spd-ldlt")))))))

(asdf:defsystem "coalton-linalg/examples"
  :description "Runnable public-interface examples."
  :depends-on ("coalton-linalg")
  :components ((:module "examples" :components ((:ct-file "spd-solve") (:ct-file "sparse-operations")))))

(asdf:defsystem "coalton-linalg/tests"
  :description "Interface, Coalton caller, and allocation regression tests."

  :depends-on ("coalton-linalg" "coalton-linalg/examples")
  :serial t
  :components ((:module "tests" :serial t
                :components ((:file "interface") (:ct-file "callers")
                             (:file "runtime") (:static-file "exact-singular-rounding.sexp") (:file "property-tests")
                             (:file "primitives") (:file "encapsulation") (:ct-file "typed-errors") (:ct-file "typed-status") (:file "suite"))))
  :perform (asdf:test-op (operation system)
             (declare (ignore operation system))
             (uiop:symbol-call :coalton-linalg/tests :run-tests)))

(asdf:defsystem "coalton-linalg/tools"
  :description "Optional fixture reader, artifact writer, and benchmark driver."
  :depends-on ("coalton-linalg")
  :components ((:module "tools" :components ((:file "fixtures")))))

(asdf:defsystem "coalton-linalg/tools-tests"
  :description "Fixture parser and artifact lifecycle regressions."
  :depends-on ("coalton-linalg/tools")
  :components ((:module "tests" :components ((:file "tools"))))
  :perform (asdf:test-op (operation system)
             (declare (ignore operation system))
             (uiop:symbol-call :coalton-linalg/tools :run-tools-tests)))