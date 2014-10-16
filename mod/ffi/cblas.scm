
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

(define (pointer-to-first A typesize)
  (bytevector->pointer (shared-array-root A) (* (shared-array-offset A) typesize)))

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

(define-syntax define-dot
  (syntax-rules ()
    ((_ name srfi4-type type cblas-name cblas-name-string)
     (begin
       (define cblas-name (pointer->procedure type (dynamic-func cblas-name-string libcblas)
                                              (list int '* int '* int)))
       (define (name A B)
         (check-2-arrays A B 1 srfi4-type)
         (cblas-name (array-length A)
                     (pointer-to-first A (sizeof type)) (first (shared-array-increments A))
                     (pointer-to-first B (sizeof type)) (first (shared-array-increments B))))))))

; double cblas_ddot (const int N, const double *X, const int incX, const double *Y, const int incY)
(define-dot ddot 'f64 double cblas_ddot "cblas_ddot")
(export cblas_ddot ddot)

; float cblas_sdot (const int N, const float *X, const int incX, const float *Y, const int incY)
(define-dot sdot 'f32 float cblas_sdot "cblas_sdot")
(export cblas_sdot sdot)

(define-syntax define-dot-complex
  (syntax-rules ()
    ((_ name srfi4-type typesize cblas-name cblas-name-string)
     (begin
       (define cblas-name (pointer->procedure void (dynamic-func cblas-name-string libcblas)
                                              (list int '* int '* int '*)))
       (define (name A B)
         (check-2-arrays A B 1 srfi4-type)
         (let ((C (make-typed-array srfi4-type *unspecified*)))
           (cblas-name (array-length A)
                       (pointer-to-first A typesize) (first (shared-array-increments A))
                       (pointer-to-first B typesize) (first (shared-array-increments B))
                       (pointer-to-first C typesize))
           (array-ref C)))))))

; void cblas_zdotu_sub (const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu)
(define-dot-complex zdotu 'c64 (* 2 (sizeof double)) cblas_zdotu_sub "cblas_zdotu_sub")
(export zdotu cblas_zdotu_sub)

; void cblas_zdotc_sub (const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc)
(define-dot-complex zdotc 'c64 (* 2 (sizeof double)) cblas_zdotc_sub "cblas_zdotc_sub")
(export zdotc cblas_zdotc_sub)

; void cblas_cdotu_sub (const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu)
(define-dot-complex cdotu 'c32 (* 2 (sizeof float)) cblas_cdotu_sub "cblas_cdotu_sub")
(export cdotu cblas_cdotu_sub)

; void cblas_cdotc_sub (const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc)
(define-dot-complex cdotc 'c32 (* 2 (sizeof float)) cblas_cdotc_sub "cblas_cdotc_sub")
(export cdotc cblas_cdotc_sub)
