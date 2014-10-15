
; Access CBLAS through Guile's FFI.
; (c) Daniel Llorens - 2014

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(define-module (ffi cblas))
(import (system foreign) (srfi srfi-1))

; @TODO As an alternative go through installation.
(define libcblas (dynamic-link (let ((lpath (getenv "GUILE_FFI_CBLAS_LIBCBLAS_PATH")))
                                 (if (and lpath (not (string=? lpath "")))
                                   (string-append lpath file-name-separator-string "libcblas")
                                   "libcblas"))))

;; Consider http://wiki.call-cc.org/eggref/4/blas#usage
;; The three levels are thus:

;; 1) ([original library name] ...) ~ (unsafe-xxx! ...) the result of
;; pointer->procedure.

;; 2) (name! ...) ~ (xxx! ...) requires compatible objects (no copies), expects
;; pre-sized arrays for the return. Increments are implicit in the array
;; object. But this means that slicing should be easier than it is in plain
;; Guile.

;; 3) (name ...) ~ (xxx ...) converts arguments as far as possible. Returns new
;; typed arrays.

; double cblas_ddot (const int N, const double *X, const int incX, const double *Y, const int incY);
(define cblas_ddot
  (pointer->procedure double (dynamic-func "cblas_ddot" libcblas)
                      (list int '* int '* int)))

(define (ddot A B)
  (unless (typed-array? A 'f64) (throw 'bad-type (array-type A)))
  (unless (typed-array? B 'f64) (throw 'bad-type (array-type B)))
  (unless (= 1 (array-rank A)) (throw 'bad-rank (array-rank A)))
  (unless (= 1 (array-rank B)) (throw 'bad-rank (array-rank B)))
  (unless (= (array-length A) (array-length B))
    (throw 'bad-sizes (array-length A) (array-length B)))
  (unless (= 0 (shared-array-offset A) (shared-array-offset B))
    (throw 'bad-sizes (array-length A) (array-length B)))
  (cblas_ddot (array-length A)
              (bytevector->pointer (shared-array-root A) (* (shared-array-offset A) (sizeof double)))
              (first (shared-array-increments A))
              (bytevector->pointer (shared-array-root B) (* (shared-array-offset B) (sizeof double)))
              (first (shared-array-increments B))))

(export cblas_ddot)
(export ddot)
