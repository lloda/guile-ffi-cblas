; -*- mode: scheme; coding: utf-8 -*-
; FFI CBLAS module.

; (c) Daniel Llorens - 2014-2015, 2017, 2019
; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU Lesser General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

;;; Commentary:
;; Access CBLAS (http://www.netlib.org/blas/#_cblas) through Guile's FFI.
;;; Code:

(define-module (ffi cblas))
(import (system foreign) (srfi srfi-1) (srfi srfi-11) (ffi blis arrays)
          (ice-9 match))

; TODO As an alternative go through installation.
(eval-when (compile load eval)
  (define libcblas (dynamic-link (let ((lpath (getenv "GUILE_FFI_CBLAS_LIBPATH"))
                                       (lname (or (getenv "GUILE_FFI_CBLAS_LIBNAME") "libcblas")))
                                   (if (and lpath (not (string=? lpath "")))
                                     (string-append lpath file-name-separator-string lname)
                                     lname))))
  (define have-rotg?
    (catch #t
      (lambda () (dynamic-func "cblas_crotg" libcblas))
      (lambda args #f))))


; -----------------------------
; wrapper utilities
; -----------------------------

(define (srfi4-type->type srfi4-type)
  (case srfi4-type
    ((f32) float)
    ((f64) double)
    ((c32 c64) '*)
    (else (throw 'no-ffi-type-for-type srfi4-type))))

(define (srfi4-type->real srfi4-type)
  (case srfi4-type
    ((f32 c32) 'f32)
    ((f64 c64) 'f64)
    (else (throw 'no-real-type-for-type srfi4-type))))

(define (srfi4-type->real-type srfi4-type)
  (case srfi4-type
    ((f32 c32) float)
    ((f64 c64) double)
    (else (throw 'no-ffi-type-for-real-type srfi4-type))))

; BUG CBLAS expects different when the inc is negative (expects pointer to size*inc element, so first in memory; not to the logically first element).
(define (pointer-to-first A)
  (bytevector->pointer (shared-array-root A)
                       (* (shared-array-offset A) (srfi4-type-size (array-type A)))))

(define (scalar->arg srfi4-type a)
  (case srfi4-type
    ((f32 f64) a)
    ((c32 c64) (bytevector->pointer (make-typed-array srfi4-type a 1) 0))
    (else (throw 'bad-array-type srfi4-type))))


; -----------------------------
; CBLAS flags
; -----------------------------

; CBLAS_ORDER
(define CblasRowMajor 101)
(define CblasColMajor 102)

; CBLAS_TRANSPOSE
(define CblasNoTrans 111)
(define CblasTrans 112)
(define CblasConjTrans 113)
(define AtlasConj 114) ; not standard

; CBLAS_UPLO
(define CblasUpper 121)
(define CblasLower 122)

; CBLAS_DIAG
(define CblasNonUnit 131)
(define CblasUnit 132)

; CBLAS_SIDE
(define CblasLeft 141)
(define CblasRight 142)

(export CblasRowMajor CblasColMajor
        CblasNoTrans CblasTrans CblasConjTrans AtlasConj
        CblasUpper CblasLower
        CblasNonUnit CblasUnit
        CblasLeft CblasRight)

(define (fliptr CBLAS_TRANSPOSE)
  (cond
   ((= CBLAS_TRANSPOSE CblasNoTrans) CblasTrans)
   ((= CBLAS_TRANSPOSE AtlasConj) CblasConjTrans)
   ((= CBLAS_TRANSPOSE CblasTrans) CblasNoTrans)
   ((= CBLAS_TRANSPOSE CblasConjTrans) AtlasConj)
   (else (throw 'bad-CBLAS_TRANSPOSE-1 CBLAS_TRANSPOSE))))

(define (tr? CBLAS_TRANSPOSE)
  (cond
   ((= CBLAS_TRANSPOSE CblasNoTrans) #f)
   ((= CBLAS_TRANSPOSE AtlasConj) #f)
   ((= CBLAS_TRANSPOSE CblasTrans) #t)
   ((= CBLAS_TRANSPOSE CblasConjTrans) #t)
   (else (throw 'bad-CBLAS_TRANSPOSE-2 CBLAS_TRANSPOSE))))


; -----------------------------
; reference & legend
; -----------------------------

;; Consider http://wiki.call-cc.org/eggref/4/blas#usage
;; The three variants per binding would be:

;; 1) ([original library name] ...) ~ (unsafe-xxx! ...) the result of
;; pointer->procedure.

;; 2) (name! ...) ~ (xxx! ...) requires compatible objects (no copies), expects
;; pre-sized arrays for the return. Increments are implicit in the array
;; object. But this means that slicing should be easier than it is in plain
;; Guile.

;; 3) (name ...) ~ (xxx ...) converts arguments as far as possible. Returns new
;; typed arrays.

#|
LEVEL 1

    Single and Double

    *   SROTG - setup Givens rotation
        SROTMG - setup modified Givens rotation
        SROT - apply Givens rotation
        SROTM - apply modified Givens rotation
    *   SSWAP - swap x and y
    *   SSCAL - x = a*x
    *   SCOPY - copy x into y
    *   SAXPY - y = a*x + y
    *   SDOT - dot product
        SDSDOT - dot product with extended precision accumulation
    *   SNRM2 - Euclidean norm
    *   SCNRM2- Euclidean norm
    *   SASUM - sum of absolute values
    *   ISAMAX - index of max abs value

    Complex and Double Complex

    *   CROTG - setup Givens rotation
        CSROT - apply Givens rotation
    *   CSWAP - swap x and y
    *   CSCAL - x = a*x
    *   CSSCAL - x = a*x
    *   CCOPY - copy x into y
    *   CAXPY - y = a*x + y
    *   CDOTU - dot product
    *   CDOTC - dot product, conjugating the first vector
    *   SCASUM - sum of absolute values
    *   ICAMAX - index of max abs value

LEVEL 2

    Single and Double

    *   SGEMV - matrix vector multiply
        SGBMV - banded matrix vector multiply
        SSYMV - symmetric matrix vector multiply
        SSBMV - symmetric banded matrix vector multiply
        SSPMV - symmetric packed matrix vector multiply
        STRMV - triangular matrix vector multiply
        STBMV - triangular banded matrix vector multiply
        STPMV - triangular packed matrix vector multiply
        STRSV - solving triangular matrix problems
        STBSV - solving triangular banded matrix problems
        STPSV - solving triangular packed matrix problems
    *   SGER - performs the rank 1 operation A := alpha*x*y' + A
        SSYR - performs the symmetric rank 1 operation A := alpha*x*x' + A
        SSPR - symmetric packed rank 1 operation A := alpha*x*x' + A
        SSYR2 - performs the symmetric rank 2 operation, A := alpha*x*y' + alpha*y*x' + A
        SSPR2 - performs the symmetric packed rank 2 operation, A := alpha*x*y' + alpha*y*x' + A

    Complex and Double Complex

    *   CGEMV - matrix vector multiply
        CGBMV - banded matrix vector multiply
        CHEMV - hermitian matrix vector multiply
        CHBMV - hermitian banded matrix vector multiply
        CHPMV - hermitian packed matrix vector multiply
        CTRMV - triangular matrix vector multiply
        CTBMV - triangular banded matrix vector multiply
        CTPMV - triangular packed matrix vector multiply
        CTRSV - solving triangular matrix problems
        CTBSV - solving triangular banded matrix problems
        CTPSV - solving triangular packed matrix problems
    *   CGERU - performs the rank 1 operation A := alpha*x*y' + A
    *   CGERC - performs the rank 1 operation A := alpha*x*conjg( y' ) + A
        CHER - hermitian rank 1 operation A := alpha*x*conjg(x') + A
        CHPR - hermitian packed rank 1 operation A := alpha*x*conjg( x' ) + A
        CHER2 - hermitian rank 2 operation
        CHPR2 - hermitian packed rank 2 operation

LEVEL 3

    Single and Double

    *   SGEMM - matrix matrix multiply
        SSYMM - symmetric matrix matrix multiply
        SSYRK - symmetric rank-k update to a matrix
        SSYR2K - symmetric rank-2k update to a matrix
        STRMM - triangular matrix matrix multiply
        STRSM - solving triangular matrix with multiple right hand sides

    Complex and Double Complex

    *   CGEMM - matrix matrix multiply
        CSYMM - symmetric matrix matrix multiply
        CHEMM - hermitian matrix matrix multiply
        CSYRK - symmetric rank-k update to a matrix
        CHERK - hermitian rank-k update to a matrix
        CSYR2K - symmetric rank-2k update to a matrix
        CHER2K - hermitian rank-2k update to a matrix
        CTRMM - triangular matrix matrix multiply
        CTRSM - solving triangular matrix with multiple right hand sides
|#


; -----------------------------
; a b -> (values c s): srotg drotg crotg zrotg
; -----------------------------

; TODO pointer-to-this-value support in the ffi, for old C decls that take double * for complex.
; int cblas_srotg ( FIXME )
(define-syntax define-rotg
  (lambda (x)
    (syntax-case x ()
      ((_ type_ cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name)))
                     (ctype #'(quote type_))
                     (stype #'(srfi4-type->real (syntax->datum #'type_)))
                     (docstring (string-append (symbol->string (syntax->datum #'cblas-name))
                                               " a b -> (values c s)")))
         (if have-rotg?
           #'(begin
               (define cblas-name
                 (pointer->procedure
                  void (dynamic-func cblas-name-string libcblas)
                  (list '* '* '* '*)))
               (define (name a b)
                 docstring
                 (let* ((a (make-typed-array ctype a))
                        (b (make-typed-array ctype b))
                        (c (make-typed-array stype *unspecified*))
                        (s (make-typed-array ctype *unspecified*)))
                   (cblas-name (pointer-to-first a)
                               (pointer-to-first b)
                               (pointer-to-first c)
                               (pointer-to-first s))
                   (values (array-ref c) (array-ref s)))))
           #'(if #f #f)))))))

(define-sdcz rotg cblas_?rotg ?rotg)


; -----------------------------
; sum_i(a_i * b_i): sdot ddot cdotu cdotc zdotu zdotc
; float cblas_sdot (const int N, const float *X, const int incX, const float *Y, const int incY)
; void cblas_cdotu_sub (const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu)
; -----------------------------

(define-syntax define-dot-real
  (lambda (x)
    (syntax-case x ()
      ((_  srfi4-type name cblas-name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name
               (pointer->procedure
                (srfi4-type->type srfi4-type) (dynamic-func cblas-name-string libcblas)
                (list int '* int '* int)))
             (define (name A B)
               (check-2-arrays A B 1 srfi4-type)
               (cblas-name (array-length A)
                           (pointer-to-first A) (stride A 0)
                           (pointer-to-first B) (stride B 0)))))))))

(define-dot-real 'f32 sdot cblas_sdot)
(define-dot-real 'f64 ddot cblas_ddot)

(export cblas_sdot cblas_ddot sdot ddot)

(define-syntax define-dot-complex
  (lambda (x)
    (syntax-case x ()
      ((_ srfi4-type cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name
               (pointer->procedure
                void (dynamic-func cblas-name-string libcblas)
                (list int '* int '* int '*)))
             (define (name A B)
               (check-2-arrays A B 1 srfi4-type)
               (let ((C (make-typed-array srfi4-type *unspecified*)))
                 (cblas-name (array-length A)
                             (pointer-to-first A) (stride A 0)
                             (pointer-to-first B) (stride B 0)
                             (pointer-to-first C))
                 (array-ref C)))))))))

(define-dot-complex 'c32 cblas_cdotu_sub cdotu)
(define-dot-complex 'c32 cblas_cdotc_sub cdotc)
(define-dot-complex 'c64 cblas_zdotu_sub zdotu)
(define-dot-complex 'c64 cblas_zdotc_sub zdotc)

(export cdotu cdotc zdotu zdotc
        cblas_cdotu_sub cblas_cdotc_sub cblas_zdotu_sub cblas_zdotc_sub)


; -----------------------------
; x -> y: scopy dcopy ccopy zcopy
; void cblas_scopy (const int N, const float *X, const int incX, float *Y, const int incY)
; -----------------------------

; TODO pointer-to-this-value support in the ffi, for old C decls that take double * for complex.
(define-syntax define-copy
  (lambda (x)
    (syntax-case x ()
      ((_ srfi4-type cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name)))
                     (type #'(quote srfi4-type)))
         #'(begin
             (define cblas-name
               (pointer->procedure
                void (dynamic-func cblas-name-string libcblas)
                (list int '* int '* int)))
             (define (name X Y)
               (check-2-arrays X Y 1 type)
               (cblas-name (array-length X)
                           (pointer-to-first X) (stride X 0)
                           (pointer-to-first Y) (stride Y 0)))))))))

(define-sdcz copy cblas_?copy ?copy!)


; -----------------------------
; x -> y: sswap dswap cswap zswap
; void cblas_sswap (const int N, const float *X, const int incX, float *Y, const int incY)
; -----------------------------

; TODO pointer-to-this-value support in the ffi, for old C decls that take double * for complex.
(define-syntax define-swap
  (lambda (x)
    (syntax-case x ()
      ((_ srfi4-type cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name)))
                     (type #'(quote srfi4-type)))
         #'(begin
             (define cblas-name
               (pointer->procedure
                void (dynamic-func cblas-name-string libcblas)
                (list int '* int '* int)))
             (define (name X Y)
               (check-2-arrays X Y 1 type)
               (cblas-name (array-length X)
                           (pointer-to-first X) (stride X 0)
                           (pointer-to-first Y) (stride Y 0)))))))))

(define-sdcz swap cblas_?swap ?swap!)


; -----------------------------
; a*x + y -> y: saxpy daxpy caxpy zaxpy
; void cblas_saxpy (const int N, const float alpha, const float *X, const int incX, float *Y, const int incY)
; -----------------------------

; TODO pointer-to-this-value support in the ffi, for old C decls that take double * for complex.
(define-syntax define-axpy
  (lambda (x)
    (syntax-case x ()
      ((_ srfi4-type cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name)))
                     (type #'(quote srfi4-type)))
         #'(begin
             (define cblas-name
               (pointer->procedure
                void (dynamic-func cblas-name-string libcblas)
                (list int (srfi4-type->type type) '* int '* int)))
             (define (name a X Y)
               (check-2-arrays X Y 1 type)
               (cblas-name (array-length X) (scalar->arg type a)
                           (pointer-to-first X) (stride X 0)
                           (pointer-to-first Y) (stride Y 0)))))))))

(define-sdcz axpy cblas_?axpy ?axpy!)


; -----------------------------
; alpha * X_i -> X_i: sscal cscal dscal zscal csscal zdscal
; void cblas_sscal (const int N, const float alpha, const float *X, const int incX)
; -----------------------------

(define-syntax define-scal
  (lambda (x)
    (syntax-case x ()
      ((_ ctype_ stype_ cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name)))
                     (ctype #'(quote ctype_))
                     (stype #'(quote stype_)))
         #'(begin
             (define cblas-name
               (pointer->procedure
                void (dynamic-func cblas-name-string libcblas)
                (list int (srfi4-type->type stype) '* int)))
             (define (name alpha X)
               (check-array X 1 ctype)
               (cblas-name (array-length X) (scalar->arg stype alpha)
                           (pointer-to-first X) (stride X 0))))))
      ((_ ctype cblas-name name)
       #'(define-scal ctype ctype cblas-name name)))))

(define-sdcz scal cblas_?scal ?scal!)

(define-scal c32 f32 cblas_csscal csscal!)
(define-scal c64 f64 cblas_zdscal zdscal!)
(export cblas_csscal cblas_zdscal csscal! zdscal!)


; -----------------------------
; sqrt(sum_i(conj(X_i)*X_i)): snrm2 dnrm2 cnrm2 znrm2
; sum('absolute value'(X_i))): sasum dasum casum zasum
; float cblas_snrm2 (const int N, const float *X, const int incX)
; -----------------------------

; 'absolute value' is |Re|+|Im| in the complex case; cf LawsonEtAl1979,
; 'Basic linear algebra subprograms for Fortran usage, p. 311.

(define-syntax define-nrm2/asum
  (lambda (x)
    (syntax-case x ()
      ((_ srfi4-type cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name
               (pointer->procedure
                (srfi4-type->real-type srfi4-type) (dynamic-func cblas-name-string libcblas)
                (list int '* int)))
             (define (name X)
               (check-array X 1 srfi4-type)
               (cblas-name (array-length X) (pointer-to-first X) (stride X 0)))))))))

(define-nrm2/asum 'f32 cblas_snrm2  snrm2)
(define-nrm2/asum 'f64 cblas_dnrm2  dnrm2)
(define-nrm2/asum 'c32 cblas_scnrm2 cnrm2)
(define-nrm2/asum 'c64 cblas_dznrm2 znrm2)
(export cblas_snrm2 cblas_dnrm2 cblas_scnrm2 cblas_dznrm2 snrm2 dnrm2 cnrm2 znrm2)

; float cblas_sasum2 (const int N, const float *X, const int incX)
(define-nrm2/asum 'f32 cblas_sasum  sasum)
(define-nrm2/asum 'f64 cblas_dasum  dasum)
(define-nrm2/asum 'c32 cblas_scasum casum)
(define-nrm2/asum 'c64 cblas_dzasum zasum)
(export cblas_sasum cblas_dasum cblas_scasum cblas_dzasum sasum dasum casum zasum)


; -----------------------------
; i | max_j('absolute value'(X_j)) = X_i: i?amax
; int cblas_isamax (const int N, const float *X, const int incX)
; -----------------------------

; 'absolute value' is |Re|+|Im| in the complex case; cf LawsonEtAl1979,
; 'Basic linear algebra subprograms for Fortran usage, p. 311.

(define-syntax define-iamax
  (lambda (x)
    (syntax-case x ()
      ((_ srfi4-type cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name)))
                     (type #'(quote srfi4-type)))
         #'(begin
             (define cblas-name
               (pointer->procedure
                int (dynamic-func cblas-name-string libcblas)
                (list int '* int)))
             (define (name X)
               (check-array X 1 type)
               (cblas-name (array-length X) (pointer-to-first X) (stride X 0)))))))))

(define-sdcz iamax cblas_i?amax i?amax)


; -----------------------------
; alpha*x_i*(maybe conj)(y_j) + A_{i, j} -> A_{i, j}: sger dger cgeru cgerc zgeru cgerc
; void cblas_sger (const enum CBLAS_ORDER Order, const int M, const int N, const float alpha, const float *X,
;                  const int incX, const float *Y, const int incY, float *A, const int lda)
; -----------------------------

(define (check-arrays-AXY A X Y type)
  (check-array A 2 type)
  (check-array X 1 type)
  (check-array Y 1 type)
  (unless (= (dim A 1) (array-length X))
    (throw 'bad-size-X (dim A 1) (array-length X)))
  (unless (= (dim A 0) (array-length Y))
    (throw 'bad-size-X (dim A 0) (array-length Y))))

(define-syntax define-ger
  (lambda (x)
    (syntax-case x ()
      ((_ srfi4-type cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name
               (pointer->procedure
                void (dynamic-func cblas-name-string libcblas)
                (list int int int (srfi4-type->type srfi4-type) '* int '* int '* int)))
             (define (name alpha X Y A)
               (check-arrays-AXY A Y X srfi4-type)
               (let-values (((A-lead A-order) (lead-and-order A)))
                 (cblas-name A-order
                             (dim A 0) (dim A 1) (scalar->arg srfi4-type alpha)
                             (pointer-to-first X) (stride X 0)
                             (pointer-to-first Y) (stride Y 0)
                             (pointer-to-first A) A-lead)))))))))

(define-ger 'f32 cblas_sger  sger!)
(define-ger 'f64 cblas_dger  dger!)
(define-ger 'c32 cblas_cgeru cgeru!)
(define-ger 'c32 cblas_cgerc cgerc!)
(define-ger 'c64 cblas_zgeru zgeru!)
(define-ger 'c64 cblas_zgerc zgerc!)
(export cblas_sger cblas_dger cblas_cgeru cblas_cgerc cblas_zgeru cblas_zgerc
        sger! dger! cgeru! cgerc! zgeru! zgerc!)


; -----------------------------
; alpha*sum_j(A_{ij} * X_j) + beta*Y_i -> Y_i: sgemv dgemv cgemv zgemv
; void cblas_sgemv (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
;                   const int M, const int N, const float alpha, const float *A, const int lda,
;                   const float *X, const int incX, const float beta, float *Y, const int incY)
; -----------------------------

(define (lead-and-order A)
  (let ((A-strides (shared-array-increments A)))
    (cond ((= 1 (first A-strides))
           (values (second A-strides) CblasColMajor))
          ((= 1 (second A-strides))
           (values (first A-strides) CblasRowMajor))
          (else (throw 'unsupported-stride-for-cblas-matrix)))))

(define-syntax define-gemv
  (lambda (x)
    (syntax-case x ()
      ((_ srfi4-type cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name)))
                     (type #'(quote srfi4-type)))
         #'(begin
             (define cblas-name
               (pointer->procedure
                void (dynamic-func cblas-name-string libcblas)
                (list int int int int (srfi4-type->type type) '* int
                      '* int (srfi4-type->type type) '* int)))
             (define (name alpha A TransA X beta Y)
               (let ((M (dim A 0))
                     (N (dim A 1)))
                 (unless (= M (array-length Y)) (throw 'mismatched-AY M (array-length Y)))
                 (unless (= N (array-length X)) (throw 'mismatched-AX N (array-length X)))
                 (let-values (((A-lead A-order) (lead-and-order A)))
                   (cblas-name A-order TransA
                               M N (scalar->arg type alpha)
                               (pointer-to-first A) A-lead
                               (pointer-to-first X) (stride X 0) (scalar->arg type beta)
                               (pointer-to-first Y) (stride Y 0)))))))))))

(define-sdcz gemv cblas_?gemv ?gemv!)


; -----------------------------
; alpha * sum_k(A_{ik}*B_{kj}) + beta * C_{ij} -> C_{ij}: sgemm dgemm cgemm zgemm
; void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
;                  const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
;                  const int K, const float alpha, const float *A,
;                  const int lda, const float *B, const int ldb,
;                  const float beta, float *C, const int ldc)
; -----------------------------

(define-syntax define-gemm
  (lambda (x)
    (syntax-case x ()
      ((_ srfi4-type cblas-name name)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name)))
                     (type #'(quote srfi4-type)))
         #'(begin
             (define cblas-name
               (pointer->procedure
                void (dynamic-func cblas-name-string libcblas)
                (list int int int int int int
                      (srfi4-type->type type) '* int '* int
                      (srfi4-type->type type) '* int)))
             (define (name alpha A TransA B TransB beta C)
               (check-array A 2 type)
               (check-array B 2 type)
               (check-array C 2 type)
               (let-values (((A-lead A-order) (lead-and-order A))
                            ((B-lead B-order) (lead-and-order B))
                            ((C-lead C-order) (lead-and-order C)))
                 (let ((M (dim C 0))
                       (N (dim C 1))
                       (K (dim A (if (tr? TransA) 0 1))))
                   (unless (= M (dim A (if (tr? TransA) 1 0))) (throw 'mismatched-CA))
                   (unless (= N (dim B (if (tr? TransB) 0 1))) (throw 'mismatched-CB))
                   (unless (= K (dim B (if (tr? TransB) 1 0))) (throw 'mismatched-AB))
                   (let ((TransA (if (eqv? C-order A-order) TransA (fliptr TransA)))
                         (TransB (if (eqv? C-order B-order) TransB (fliptr TransB))))
                     (cblas-name C-order TransA TransB M N K
                                 (scalar->arg type alpha)
                                 (pointer-to-first A) A-lead
                                 (pointer-to-first B) B-lead
                                 (scalar->arg type beta)
                                 (pointer-to-first C) C-lead)))))))))))

(define-sdcz gemm cblas_?gemm ?gemm!)
