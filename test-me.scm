; run this with guile -s test-me.scm
; or GUILE_FFI_CBLAS_LIBCBLAS_NAME=libblas GUILE_FFI_CBLAS_LIBCBLAS_PATH=/usr/lib $GUILE -s test-me.scm

(add-to-load-path (string-append (getcwd) "/mod"))
(load "test/test-ffi-cblas.scm")
;; (load "test/test-ffi-blis.scm")
