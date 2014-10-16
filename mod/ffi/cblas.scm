
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

(define (check-array A rank type)
  (unless (= rank (array-rank A)) (throw 'bad-rank (array-rank A)))
  (unless (typed-array? A type) (throw 'bad-type type (array-type A))))

(define (check-2-arrays A B rank type)
  (check-array A rank type)
  (check-array B rank type)
  (unless (= (array-length A) (array-length B))
    (throw 'bad-sizes (array-length A) (array-length B)))
  (unless (= 0 (caar (array-shape A)) (caar (array-shape B)))
    (throw 'bad-base-indices (array-length A) (array-length B))))

(define (pointer-to-first A type)
  (bytevector->pointer (shared-array-root A) (* (shared-array-offset A) (sizeof type))))

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
  (check-2-arrays A B 1 'f64)
  (cblas_ddot (array-length A)
              (pointer-to-first A double) (first (shared-array-increments A))
              (pointer-to-first B double) (first (shared-array-increments B))))

(export cblas_ddot)
(export ddot)
