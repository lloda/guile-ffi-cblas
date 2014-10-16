; run this with guile -s test-me.scm

(add-to-load-path (string-append (getcwd) "/mod"))
(load "test/test-ffi-cblas.scm")
