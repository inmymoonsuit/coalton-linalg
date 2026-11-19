(defpackage #:coalton-linalg/tests/encapsulation
  (:use #:cl)
  (:export #:run-encapsulation-tests #:run-status-type-tests #:run-typed-status-tests))
(in-package #:coalton-linalg/tests/encapsulation)

(defun compile-probe (name body &key execute)
  "Compile a real external Coalton caller; preserve diagnostics for accepted and rejected probes."
  (let* ((directory (asdf:system-relative-pathname "coalton-linalg" "build/interface-probes/"))
         (source (merge-pathnames (format nil "~A.lisp" name) directory))
         (binary (merge-pathnames (format nil "~A.fasl" name) directory))
         (log (merge-pathnames (format nil "~A.log" name) directory))
         (package-name (format nil "COALTON-LINALG/EXTERNAL-PROBE/~A" (string-upcase name)))
         (accepted nil)
         (diagnostics
           (with-output-to-string (stream)
             (ensure-directories-exist source)
             (with-open-file (out source :direction :output :if-exists :supersede)
               (format out "(cl:defpackage ~S~%  (:use #:coalton #:coalton-prelude)~%  (:local-nicknames (#:ldlt #:coalton-linalg/sparse/spd-ldlt)~%                    (#:arr #:coalton/lisparray)))~%(cl:in-package ~S)~%(coalton-toplevel~%~A)~%"
                       package-name package-name body))
             (let ((*standard-output* stream)
                   (*error-output* stream)
                   (*trace-output* stream)
                   (*compile-verbose* nil)
                   (*compile-print* nil)
                   (*package* (find-package :cl-user)))
               (handler-case
                   (multiple-value-bind (output warnings failure)
                       (compile-file source :output-file binary)
                     (declare (ignore warnings))
                     (setf accepted (and output (not failure)))
                     (when (and accepted execute)
                       (load output :verbose nil :print nil)
                       (funcall (symbol-function (find-symbol "PROBE" package-name)))))
                 (error (condition)
                   (setf accepted nil)
                   (format stream "~&PROBE-CONDITION ~A: ~A~%"
                           (type-of condition) condition)))))))
    (with-open-file (out log :direction :output :if-exists :supersede)
      (write-string diagnostics out))
    (format t "EXTERNAL-CALLER ~A ~A diagnostics=~A~%"
            name (if accepted "ACCEPTED" "REJECTED") log)
    (values accepted diagnostics)))

(defparameter *public-ownership-control*
  "(declare probe (Void -> ldlt:NumericFactor))
  (define (probe)
    (ldlt:symbolic-factor
      (ldlt:analyse
        (ldlt:Pattern 0 (arr:make 1 (the UFix 0)) (arr:make 0 (the UFix 0))))))")

(defparameter *encapsulation-probes*
  '(("factor-construction" .
     "(declare probe (Void -> ldlt:NumericFactor))
      (define (probe)
        (ldlt:NumericFactor 0 (arr:make 0 (the UFix 0))
                             (arr:make 1 (the UFix 0))
                             (arr:make 0 (the UFix 0))
                             (arr:make 0 0.0d0)
                             (arr:make 0 0.0d0)
                             (arr:make 0 0.0d0)
                             (arr:make 1 (the UFix 0))))")
    ("factor-destructuring" .
     "(declare probe (ldlt:NumericFactor -> UFix))
      (define (probe factor)
        (match factor
          ((ldlt:NumericFactor _ _ _ _ _ _ _ state) (arr:length state))))")
    ("symbolic-construction" .
     "(declare probe (ldlt:NumericFactor * ldlt:FactorResult -> ldlt:Symbolic))
      (define (probe factor result)
        (ldlt:Symbolic factor 0
                       (arr:make 1 (the UFix 0))
                       (arr:make 0 (the UFix 0))
                       (arr:make 0 (the UFix 0))
                       (arr:make 1 (the UFix 0))
                       (arr:make 0 (the UFix 0))
                       (arr:make 0 (the UFix 0))
                       result (arr:make 0 result) (arr:make 0 result)))")
    ("symbolic-destructuring" .
     "(declare probe (ldlt:Symbolic -> ldlt:NumericFactor))
      (define (probe symbolic)
        (match symbolic
          ((ldlt:Symbolic factor _ _ _ _ _ _ _ _ _ _) factor)))")
    ("factor-state-access" .
     "(declare probe (ldlt:NumericFactor -> UFix))
      (define (probe factor) (arr:aref (.state factor) 0))")
    ("factor-scratch-access" .
     "(declare probe (ldlt:NumericFactor -> arr:LispArray F64))
      (define (probe factor) (.scratch factor))")
    ("symbolic-factor-access" .
     "(declare probe (ldlt:Symbolic -> ldlt:NumericFactor))
      (define (probe symbolic) (.numeric symbolic))")
    ("symbolic-result-access" .
     "(declare probe (ldlt:Symbolic -> ldlt:FactorResult))
      (define (probe symbolic) (.success symbolic))")))

(defun require-accepted-control (name body &key execute)
  (multiple-value-bind (accepted diagnostics) (compile-probe name body :execute execute)
    (assert accepted () "Valid external caller ~A failed:~%~A" name diagnostics)))

(defun require-rejected-probes (probes)
  "Aggregate all unexpected acceptances so the red run records the complete exposed family."
  (let ((accepted nil))
    (dolist (probe probes)
      (when (compile-probe (car probe) (cdr probe))
        (push (car probe) accepted)))
    (when accepted
      (error "External callers compiled despite required rejection: ~{~A~^, ~}."
             (nreverse accepted)))))

(defun require-internal-storage-constructors ()
  "Check the real package export registry as well as external compilation behavior."
  (dolist (name '("MKNUMERICFACTOR" "MKSYMBOLIC" "NUMERICFACTORSTORAGE" "SYMBOLICSTORAGE"))
    (multiple-value-bind (symbol visibility)
        (find-symbol name "COALTON-LINALG/SPARSE/SPD-LDLT")
      (assert (and symbol (eq visibility :internal)) ()
              "Storage constructor ~A must exist and remain internal; visibility is ~S."
              name visibility)))
  (format t "PRIVATE-CONSTRUCTOR-EXPORTS-PASS count=4~%"))

(defun run-encapsulation-tests ()
  "The pre-fix public constructors make this regression red; ordinary public creation must work."
  (require-accepted-control "ownership-public-control" *public-ownership-control* :execute t)
  (require-rejected-probes *encapsulation-probes*)
  (require-internal-storage-constructors)
  (format t "ENCAPSULATION-PASS cases=8 public-control=1~%")
  t)

(defparameter *public-status-control*
  "(declare probe (ldlt:FactorStatus * ldlt:SolveStatus -> UFix))
   (define (probe factor-status solve-status)
     (+ (match factor-status
          ((ldlt:FactorSuccess) 1)
          ((ldlt:FactorNonPositive) 2)
          ((ldlt:FactorNonFinite) 3)
          ((ldlt:FactorBadValues) 4)
          ((ldlt:FactorBadShift) 5))
        (match solve-status
          ((ldlt:SolveSuccess) 10)
          ((ldlt:SolveInvalidFactor) 20)
          ((ldlt:SolveBadSize) 30)
          ((ldlt:SolveNonFiniteRhs) 40)
          ((ldlt:SolveNonFiniteResult) 50))))")

(defparameter *status-type-probes*
  '(("factor-status-as-solve-status" .
     "(declare probe (ldlt:FactorStatus -> ldlt:SolveStatus))
      (define (probe status) status)")
    ("solve-status-as-factor-status" .
     "(declare probe (ldlt:SolveStatus -> ldlt:FactorStatus))
      (define (probe status) status)")
    ("solve-return-as-factor-status" .
     "(declare probe (ldlt:NumericFactor * arr:LispArray F64 * arr:LispArray F64 -> ldlt:FactorStatus))
      (define (probe factor rhs output) (ldlt:solve! factor rhs output))")
    ("factor-return-as-solve-status" .
     "(declare probe (ldlt:Symbolic * arr:LispArray F64 * F64 -> ldlt:SolveStatus))
      (define (probe symbolic entries shift)
        (ldlt:factor-status (ldlt:factor! symbolic entries shift)))")))

(defun run-status-type-tests ()
  "Compile a valid exhaustive match before rejecting cross-family outcome assignments."
  (require-accepted-control "status-public-control" *public-status-control*)
  (require-rejected-probes *status-type-probes*)
  (format t "STATUS-TYPE-SEPARATION-PASS cases=4 public-control=1~%")
  t)

(defun run-typed-status-tests ()
  "Require all ten actual factor/solve outcomes to reach their external Coalton match branches."
  (let ((actual (uiop:symbol-call :coalton-linalg/tests/typed-status :exercise-status-branches)))
    (assert (equalp actual #(0 1 2 3 4 0 1 2 3 4)) ()
            "Typed status branches differ: ~S." actual))
  (format t "TYPED-STATUS-BRANCHES-PASS factor=5 solve=5~%")
  t)
