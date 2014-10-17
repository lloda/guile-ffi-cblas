
; Access CBLAS through Guile's FFI.
; (c) Daniel Llorens - 2014

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(define-module (ffi cblas))
(import (system foreign) (srfi srfi-1) (srfi srfi-11))

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

(define (srfi4-type-size srfi4-type)
  (case srfi4-type
    ((f32) 4)
    ((f64 c32) 8)
    ((c64) 16)
    (else (throw 'bad-array-type srfi4-type))))

(define (srfi4-type->type srfi4-type)
  (case srfi4-type
    ((f32) float)
    ((f64) double)
    ((c32 c64) '*)
    (else (throw 'no-ffi-type-for-type srfi4-type))))

(define (pointer-to-first A)
  (bytevector->pointer (shared-array-root A)
                       (* (shared-array-offset A) (srfi4-type-size (array-type A)))))

(define (scalar->cblas-arg srfi4-type a)
  (case srfi4-type
    ((f32 f64) a)
    ((c32 c64) (pointer-to-first (make-typed-array srfi4-type a)))
    (else (throw 'bad-array-type srfi4-type))))

(define (stride A i)
  (list-ref (shared-array-increments A) i))

(define (dim A i)
  (list-ref (array-dimensions A) i))

;; Consider http://wiki.call-cc.org/eggref/4/blas#usage
;; The three levels would be:

;; 1) ([original library name] ...) ~ (unsafe-xxx! ...) the result of
;; pointer->procedure.

;; 2) (name! ...) ~ (xxx! ...) requires compatible objects (no copies), expects
;; pre-sized arrays for the return. Increments are implicit in the array
;; object. But this means that slicing should be easier than it is in plain
;; Guile.

;; 3) (name ...) ~ (xxx ...) converts arguments as far as possible. Returns new
;; typed arrays.

; -----------------------------
; sum_i(a_i * b_i): sdot ddot cdotu cdotc zdotu zdotc
; -----------------------------

(define-syntax define-dot-real
  (syntax-rules ()
    ((_ name srfi4-type cblas-name cblas-name-string)
     (begin
       (define cblas-name (pointer->procedure (srfi4-type->type srfi4-type)
                                              (dynamic-func cblas-name-string libcblas)
                                              (list int '* int '* int)))
       (define (name A B)
         (check-2-arrays A B 1 srfi4-type)
         (cblas-name (array-length A)
                     (pointer-to-first A) (stride A 0)
                     (pointer-to-first B) (stride B 0)))))))

; double cblas_ddot (const int N, const double *X, const int incX, const double *Y, const int incY)
(define-dot-real ddot 'f64 cblas_ddot "cblas_ddot")
(export cblas_ddot ddot)

; float cblas_sdot (const int N, const float *X, const int incX, const float *Y, const int incY)
(define-dot-real sdot 'f32 cblas_sdot "cblas_sdot")
(export cblas_sdot sdot)

(define-syntax define-dot-complex
  (syntax-rules ()
    ((_ name srfi4-type cblas-name cblas-name-string)
     (begin
       (define cblas-name (pointer->procedure void
                                              (dynamic-func cblas-name-string libcblas)
                                              (list int '* int '* int '*)))
       (define (name A B)
         (check-2-arrays A B 1 srfi4-type)
         (let ((C (make-typed-array srfi4-type *unspecified*)))
           (cblas-name (array-length A)
                       (pointer-to-first A) (stride A 0)
                       (pointer-to-first B) (stride B 0)
                       (pointer-to-first C))
           (array-ref C)))))))

; void cblas_zdotu_sub (const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu)
(define-dot-complex zdotu 'c64 cblas_zdotu_sub "cblas_zdotu_sub")
(export zdotu cblas_zdotu_sub)

; void cblas_zdotc_sub (const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc)
(define-dot-complex zdotc 'c64 cblas_zdotc_sub "cblas_zdotc_sub")
(export zdotc cblas_zdotc_sub)

; void cblas_cdotu_sub (const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu)
(define-dot-complex cdotu 'c32 cblas_cdotu_sub "cblas_cdotu_sub")
(export cdotu cblas_cdotu_sub)

; void cblas_cdotc_sub (const int N, const void *X, const int incX, const void *Y, const int incY, void *dotc)
(define-dot-complex cdotc 'c32 cblas_cdotc_sub "cblas_cdotc_sub")
(export cdotc cblas_cdotc_sub)

; -----------------------------
; a*x + y -> y: saxpy daxpy caxpy zaxpy
; -----------------------------

; @TODO pointer-to-this-value support in the ffi, for old C decls that take double * for complex.
(define-syntax define-axpy
  (syntax-rules ()
    ((_ name srfi4-type cblas-name cblas-name-string)
     (begin
       (define cblas-name (pointer->procedure void
                                              (dynamic-func cblas-name-string libcblas)
                                              (list int (srfi4-type->type srfi4-type) '* int '* int)))
       (define (name a X Y)
         (check-2-arrays X Y 1 srfi4-type)
         (cblas-name (array-length X) (scalar->cblas-arg srfi4-type a)
                     (pointer-to-first X) (stride X 0)
                     (pointer-to-first Y) (stride Y 0)))))))

; void cblas_daxpy (const int N, const double alpha, const double *X, const int incX, double *Y, const int incY)
(define-axpy daxpy! 'f64 cblas_daxpy "cblas_daxpy")
(export cblas_daxpy daxpy!)

; void cblas_saxpy (const int N, const float alpha, const float *X, const int incX, float *Y, const int incY)
(define-axpy saxpy! 'f32 cblas_saxpy "cblas_saxpy")
(export cblas_saxpy saxpy!)

; void cblas_zaxpy (const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY)
(define-axpy zaxpy! 'c64 cblas_zaxpy "cblas_zaxpy")
(export cblas_zaxpy zaxpy!)

; void cblas_caxpy (const int N, const void *alpha, const void *X, const int incX, void *Y, const int incY)
(define-axpy caxpy! 'c32 cblas_caxpy "cblas_caxpy")
(export cblas_caxpy caxpy!)

; -----------------------------
; alpha*sum_j(A_{ij}*X_j) + beta*Y_i -> Y_i: sgemv dgemv cgemv zgemv
; -----------------------------

; CBLAS_ORDER
(define CblasRowMajor 101)
(define CblasColMajor 102)

; CBLAS_TRANSPOSE
(define CblasNoTrans 111)
(define CblasTrans 112)
(define CblasConjTrans 113)

; CBLAS_UPLO
(define CblasUpper 121)
(define CblasLower 122)

; CBLAS_DIAG
(define CblasNonUnit 131)
(define CblasUnit 132)

; CBLAS_SIDE
(define CblasLeft 141)
(define CblasRight 142)

(define (lead-and-order A)
  (let ((A-strides (shared-array-increments A)))
    (cond ((= 1 (first A-strides))
           (values (second A-strides) CblasColMajor))
          ((= 1 (second A-strides))
           (values (first A-strides) CblasRowMajor))
          (else (throw 'unsupported-stride-for-cblas-matrix)))))

(define (check-arrays-AXY A X Y type)
  (check-array A 2 type)
  (check-array X 1 type)
  (check-array Y 1 type)
  (unless (= (dim A 1) (array-length X))
    (throw 'bad-size-X (dim A 1) (array-length X)))
  (unless (= (dim A 0) (array-length Y))
    (throw 'bad-size-X (dim A 0) (array-length Y))))

(define-syntax define-gemv
  (syntax-rules ()
    ((_ name srfi4-type cblas-name cblas-name-string)
     (begin
       (define cblas-name (pointer->procedure
                           void
                           (dynamic-func cblas-name-string libcblas)
                           (list int int int int (srfi4-type->type srfi4-type) '* int
                                 '* int (srfi4-type->type srfi4-type) '* int)))
       (define (name alpha A X beta Y)
         (check-arrays-AXY A X Y srfi4-type)
         (let-values (((A-lead A-order) (lead-and-order A)))
           (cblas-name A-order CblasNoTrans
                       (dim A 0) (dim A 1) (scalar->cblas-arg srfi4-type alpha)
                       (pointer-to-first A) A-lead
                       (pointer-to-first X) (stride X 0) (scalar->cblas-arg srfi4-type beta)
                       (pointer-to-first Y) (stride Y 0))))))))

; void cblas_dgemv (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
;                   const int M, const int N, const double alpha, const double *A, const int lda,
;                   const double *X, const int incX, const double beta, double *Y, const int incY)
(define-gemv dgemv! 'f64 cblas_dgemv "cblas_dgemv")
(export cblas_dgemv dgemv!)

; void cblas_sgemv (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
;                   const int M, const int N, const float alpha, const float *A, const int lda,
;                   const float *X, const int incX, const float beta, float *Y, const int incY)
(define-gemv sgemv! 'f32 cblas_sgemv "cblas_sgemv")
(export cblas_sgemv sgemv!)

; void cblas_zgemv (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
;                   const int M, const int N, const void * alpha, const void *A, const int lda,
;                   const void *X, const int incX, const void * beta, void *Y, const int incY)
(define-gemv zgemv! 'c64 cblas_zgemv "cblas_zgemv")
(export cblas_zgemv zgemv!)

; void cblas_cgemv (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
;                   const int M, const int N, const void * alpha, const void *A, const int lda,
;                   const void *X, const int incX, const void * beta, void *Y, const int incY)
(define-gemv cgemv! 'c32 cblas_cgemv "cblas_cgemv")
(export cblas_cgemv cgemv!)
