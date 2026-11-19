(defpackage #:coalton-linalg/toolchain
  (:use #:cl)
  (:export #:write-metadata))

(in-package #:coalton-linalg/toolchain)

(defun write-metadata (output)
  "Record loaded systems and compiler settings; do not load or compile systems."
  (with-open-file (stream output :direction :output :if-exists :supersede)
    (format stream "SBCL~C~A~%" #\Tab (lisp-implementation-version))
    (format stream "ASDF~C~A~%" #\Tab (asdf:asdf-version))
    (flet ((setting (key value)
             (format stream "ENV~C~A~C~A~%" #\Tab key #\Tab value)))
      (setting "COALTON_ENV" (uiop:getenv "COALTON_ENV"))
      (setting "COALTON_RELEASE" (if (coalton-impl/settings:coalton-release-p) "true" "false"))
      (setting "COALTON_SPECIALIZATION"
               (if coalton-impl/settings:*coalton-disable-specialization* "false" "true"))
      (setting "COALTON_TYPE_ANNOTATIONS"
               (if coalton-impl/settings:*emit-type-annotations* "true" "false")))
    (dolist (name (sort (copy-list (asdf:already-loaded-systems)) #'string<))
      (let* ((system (asdf:find-system name nil))
             (source (and system (asdf:system-source-file system))))
        (format stream "SYSTEM~C~A~C~A~C~A~%" #\Tab name #\Tab
                (or (and system (asdf:component-version system)) "") #\Tab
                (if source (namestring source) "<built-in>")))))
  output)
