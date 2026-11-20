(asdf:defsystem "coalton-linalg"
  :description "Numerical linear algebra in Coalton: sparse CSC operations, Gram products and SPD LDL factorization."
  :version "0.9.1"
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
                              :components ((:ct-file "csc") (:ct-file "limits") (:ct-file "ordering")
                                           (:ct-file "gram") (:ct-file "scaling") (:ct-file "augmented")
                                           (:ct-file "spd-ldlt")))
                             (:module "dense" :serial t
                              :components ((:ct-file "limits") (:ct-file "indefinite-ldlt")
                                           (:ct-file "orthogonal") (:ct-file "jacobi")))))))

(asdf:defsystem "coalton-linalg/examples"
  :description "Runnable public-interface examples."
  :depends-on ("coalton-linalg")
  :components ((:module "examples" :components ((:ct-file "spd-solve") (:ct-file "sparse-operations")
                                                (:ct-file "normal-equations") (:ct-file "minimum-norm") (:ct-file "dependent-constraints") (:ct-file "singular-rank")
                                                (:ct-file "augmented-minimum-norm") (:ct-file "scaled-solve")))))

(asdf:defsystem "coalton-linalg/tests"
  :description "Interface, Coalton caller, and allocation regression tests."

  :depends-on ("coalton-linalg" "coalton-linalg/examples")
  :serial t
  :components ((:module "tests" :serial t
                :components ((:file "interface") (:ct-file "callers")
                             (:file "runtime") (:static-file "exact-singular-rounding.sexp") (:file "property-tests")
                             (:file "primitives") (:file "encapsulation") (:ct-file "typed-errors") (:ct-file "typed-status")
                             (:ct-file "buffered-callers") (:file "buffered-products")
                             (:file "buffered-properties")
                             (:ct-file "gram-callers") (:file "gram-runtime")
                             (:file "gram-properties")
                             (:ct-file "condition-callers") (:file "condition")
                             (:ct-file "indefinite-callers") (:file "indefinite")
                             (:ct-file "orthogonal-callers") (:file "orthogonal")
                             (:ct-file "jacobi-callers") (:file "jacobi")
                             (:ct-file "cholesky-callers") (:file "cholesky-form")
                             (:ct-file "quasidefinite-callers") (:file "quasidefinite") (:file "ordering")
                             (:ct-file "scaling-callers") (:file "scaling")
                             (:file "suite"))))
  :perform (asdf:test-op (operation system)
             (declare (ignore operation system))
             (uiop:symbol-call :coalton-linalg/tests :run-tests)))

(asdf:defsystem "coalton-linalg/tools"
  :description "Optional fixture reader, artifact writer, and benchmark driver."
  :depends-on ("coalton-linalg")
  :components ((:module "tools" :serial t :components ((:file "fixtures") (:file "gram-fixtures")
                                                       (:file "indefinite-fixtures") (:file "orthogonal-fixtures")
                                                       (:file "jacobi-fixtures") (:file "augmented-fixtures")))))

(asdf:defsystem "coalton-linalg/tools-tests"
  :description "Fixture parser and artifact lifecycle regressions."
  :depends-on ("coalton-linalg/tools")
  :components ((:module "tests" :serial t :components ((:file "tools") (:file "gram-tools") (:file "indefinite-tools")
                                                       (:file "orthogonal-tools") (:file "jacobi-tools")
                                                       (:file "augmented-tools"))))
  :perform (asdf:test-op (operation system)
             (declare (ignore operation system))
             (uiop:symbol-call :coalton-linalg/tools :run-tools-tests)
             (uiop:symbol-call :coalton-linalg/tools :run-gram-tools-tests)
             (uiop:symbol-call :coalton-linalg/tools :run-indefinite-tools-tests)
             (uiop:symbol-call :coalton-linalg/tools :run-orthogonal-tools-tests)
             (uiop:symbol-call :coalton-linalg/tools :run-jacobi-tools-tests)
             (uiop:symbol-call :coalton-linalg/tools :run-augmented-tools-tests)))