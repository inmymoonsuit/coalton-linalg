(load (merge-pathnames "bootstrap.lisp" *load-truename*))
(asdf:load-system "coalton-linalg/tools")
(load (merge-pathnames "allocation-isolation.lisp" *load-truename*))
(call-with-allocation-isolation
  (lambda ()
    (uiop:symbol-call :coalton-linalg/tools :run-jacobi-fixtures
      (uiop:getenv "COALTON_LINALG_INPUT") (uiop:getenv "COALTON_LINALG_OUTPUT")
      :repeats (parse-integer (uiop:getenv "COALTON_LINALG_REPEATS"))
      :allocation-repeats (parse-integer (uiop:getenv "COALTON_LINALG_ALLOCATION_REPEATS")))))
