
; Access CBLAS through Guile's FFI.
; (c) Daniel Llorens - 2014-2015, 2017, 2019

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(define-module (ffi cblas))
(import (system foreign) (srfi srfi-1) (srfi srfi-11) (ffi arrays))

; TODO As an alternative go through installation.
(define libcblas (dynamic-link (let ((lpath (getenv "GUILE_FFI_CBLAS_LIBCBLAS_PATH"))
                                     (lname (or (getenv "GUILE_FFI_CBLAS_LIBCBLAS_NAME") "libcblas")))
                                 (if (and lpath (not (string=? lpath "")))
                                   (string-append lpath file-name-separator-string lname)
                                   lname))))


; -----------------------------
; wrapper utilities
; -----------------------------

(define (srfi4-type->type srfi4-type)
  (case srfi4-type
    ((f32) float)
    ((f64) double)
    ((c32 c64) '*)
    (else (throw 'no-ffi-type-for-type srfi4-type))))

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
    ((c32 c64) (pointer-to-first (make-typed-array srfi4-type a)))
    (else (throw 'bad-array-type srfi4-type))))


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
(define-syntax define-rotg
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name ctype stype)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name)))
                     (docstring (string-append (symbol->string (syntax->datum #'cblas-name))
                                               " a b -> (values c s)")))
         #'(begin
             (define cblas-name (pointer->procedure void
                                                    (dynamic-func cblas-name-string libcblas)
                                                    (list '* '* '* '*)))
             (define (name a b)
               docstring
               (let ((a (make-typed-array ctype a))
                     (b (make-typed-array ctype b))
                     (c (make-typed-array stype *unspecified*))
                     (s (make-typed-array ctype *unspecified*)))
                 (cblas-name (pointer-to-first a)
                             (pointer-to-first b)
                             (pointer-to-first c)
                             (pointer-to-first s))
                 (values (array-ref c) (array-ref s))))))))))

; int cblas_isamax (const int N, const float *X, const int incX)
(define-rotg srotg cblas_srotg 'f32 'f32)
(define-rotg drotg cblas_drotg 'f64 'f64)
(define-rotg crotg cblas_crotg 'c32 'f32)
(define-rotg zrotg cblas_zrotg 'c64 'f64)

(export cblas_srotg cblas_drotg cblas_crotg cblas_zrotg
        srotg drotg crotg zrotg)


; -----------------------------
; x -> i: isamax idamax icamax izamax
; -----------------------------

; TODO pointer-to-this-value support in the ffi, for old C decls that take double * for complex.
(define-syntax define-iamax
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure int
                                                    (dynamic-func cblas-name-string libcblas)
                                                    (list int '* int)))
             (define (name X)
               (check-array X 1 srfi4-type)
               (cblas-name (array-length X)
                           (pointer-to-first X) (stride X 0)))))))))

; int cblas_isamax (const int N, const float *X, const int incX)
(define-iamax isamax cblas_isamax 'f32)
(define-iamax idamax cblas_idamax 'f64)
(define-iamax icamax cblas_icamax 'c32)
(define-iamax izamax cblas_izamax 'c64)

(export cblas_isamax cblas_idamax cblas_icamax cblas_izamax
        isamax idamax icamax izamax)


; -----------------------------
; sum_i(a_i * b_i): sdot ddot cdotu cdotc zdotu zdotc
; -----------------------------

(define-syntax define-dot-real
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure (srfi4-type->type srfi4-type)
                                                    (dynamic-func cblas-name-string libcblas)
                                                    (list int '* int '* int)))
             (define (name A B)
               (check-2-arrays A B 1 srfi4-type)
               (cblas-name (array-length A)
                           (pointer-to-first A) (stride A 0)
                           (pointer-to-first B) (stride B 0)))))))))

; float cblas_sdot (const int N, const float *X, const int incX, const float *Y, const int incY)
(define-dot-real sdot cblas_sdot 'f32)
(define-dot-real ddot cblas_ddot 'f64)

(export cblas_sdot cblas_ddot
        sdot ddot)

(define-syntax define-dot-complex
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
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
                 (array-ref C)))))))))

; void cblas_cdotu_sub (const int N, const void *X, const int incX, const void *Y, const int incY, void *dotu)
(define-dot-complex cdotu cblas_cdotu_sub 'c32)
(define-dot-complex cdotc cblas_cdotc_sub 'c32)
(define-dot-complex zdotu cblas_zdotu_sub 'c64)
(define-dot-complex zdotc cblas_zdotc_sub 'c64)

(export cdotu cdotc zdotu zdotc
        cblas_cdotu_sub cblas_cdotc_sub cblas_zdotu_sub cblas_zdotc_sub)


; -----------------------------
; x -> y: scopy dcopy ccopy zcopy
; -----------------------------

; TODO pointer-to-this-value support in the ffi, for old C decls that take double * for complex.
(define-syntax define-copy
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure void
                                                    (dynamic-func cblas-name-string libcblas)
                                                    (list int '* int '* int)))
             (define (name X Y)
               (check-2-arrays X Y 1 srfi4-type)
               (cblas-name (array-length X)
                           (pointer-to-first X) (stride X 0)
                           (pointer-to-first Y) (stride Y 0)))))))))

; void cblas_scopy (const int N, const float *X, const int incX, float *Y, const int incY)
(define-copy scopy! cblas_scopy 'f32)
(define-copy dcopy! cblas_dcopy 'f64)
(define-copy ccopy! cblas_ccopy 'c32)
(define-copy zcopy! cblas_zcopy 'c64)

(export cblas_scopy cblas_dcopy cblas_ccopy cblas_zcopy
        scopy! dcopy! ccopy! zcopy!)


; -----------------------------
; x -> y: sswap dswap cswap zswap
; -----------------------------

; TODO pointer-to-this-value support in the ffi, for old C decls that take double * for complex.
(define-syntax define-swap
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure void
                                                    (dynamic-func cblas-name-string libcblas)
                                                    (list int '* int '* int)))
             (define (name X Y)
               (check-2-arrays X Y 1 srfi4-type)
               (cblas-name (array-length X)
                           (pointer-to-first X) (stride X 0)
                           (pointer-to-first Y) (stride Y 0)))))))))

; void cblas_sswap (const int N, const float *X, const int incX, float *Y, const int incY)
(define-swap sswap! cblas_sswap 'f32)
(define-swap dswap! cblas_dswap 'f64)
(define-swap cswap! cblas_cswap 'c32)
(define-swap zswap! cblas_zswap 'c64)

(export cblas_sswap cblas_dswap cblas_cswap cblas_zswap
        sswap! dswap! cswap! zswap!)


; -----------------------------
; a*x + y -> y: saxpy daxpy caxpy zaxpy
; -----------------------------

; TODO pointer-to-this-value support in the ffi, for old C decls that take double * for complex.
(define-syntax define-axpy
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure void
                                                    (dynamic-func cblas-name-string libcblas)
                                                    (list int (srfi4-type->type srfi4-type) '* int '* int)))
             (define (name a X Y)
               (check-2-arrays X Y 1 srfi4-type)
               (cblas-name (array-length X) (scalar->arg srfi4-type a)
                           (pointer-to-first X) (stride X 0)
                           (pointer-to-first Y) (stride Y 0)))))))))

; void cblas_saxpy (const int N, const float alpha, const float *X, const int incX, float *Y, const int incY)
(define-axpy saxpy! cblas_saxpy 'f32)
(define-axpy daxpy! cblas_daxpy 'f64)
(define-axpy caxpy! cblas_caxpy 'c32)
(define-axpy zaxpy! cblas_zaxpy 'c64)

(export cblas_saxpy cblas_daxpy cblas_caxpy cblas_zaxpy
        saxpy! daxpy! caxpy! zaxpy!)


; -----------------------------
; alpha * X_i -> X_i: sscal cscal dscal zscal csscal zdscal
; -----------------------------

(define-syntax define-scal
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name srfi4-type scalar-srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure
                                 void
                                 (dynamic-func cblas-name-string libcblas)
                                 (list int (srfi4-type->type scalar-srfi4-type) '* int)))
             (define (name alpha X)
               (check-array X 1 srfi4-type)
               (cblas-name (array-length X) (scalar->arg scalar-srfi4-type alpha)
                           (pointer-to-first X) (stride X 0)))))))))

; void cblas_sscal (const int N, const float alpha, const float *X, const int incX)
(define-scal sscal! cblas_sscal 'f32 'f32)
(define-scal dscal! cblas_dscal 'f64 'f64)
(define-scal cscal! cblas_cscal 'c32 'c32)
(define-scal zscal! cblas_zscal 'c64 'c64)
(define-scal csscal! cblas_csscal 'c32 'f32)
(define-scal zdscal! cblas_zdscal 'c64 'f64)

(export cblas_sscal cblas_dscal cblas_cscal cblas_zscal cblas_csscal cblas_zdscal
        sscal! dscal! cscal! zscal! csscal! zdscal!)


; -----------------------------
; sqrt(sum_i(conj(X_i)*X_i)): snrm2 dnrm2 cnrm2 znrm2
; sum('absolute value'(X_i))): sasum dasum casum zasum
; 'absolute value' is |Re|+|Im| in the complex case; cf LawsonEtAl1979,
; 'Basic linear algebra subprograms for Fortran usage, p. 311.
; -----------------------------

(define-syntax define-nrm2/asum
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure
                                 (srfi4-type->real-type srfi4-type)
                                 (dynamic-func cblas-name-string libcblas)
                                 (list int '* int)))
             (define (name X)
               (check-array X 1 srfi4-type)
               (cblas-name (array-length X) (pointer-to-first X) (stride X 0)))))))))

; float cblas_snrm2 (const int N, const float *X, const int incX)
(define-nrm2/asum snrm2 cblas_snrm2 'f32)
(define-nrm2/asum dnrm2 cblas_dnrm2 'f64)
(define-nrm2/asum cnrm2 cblas_scnrm2 'c32)
(define-nrm2/asum znrm2 cblas_dznrm2 'c64)

(export cblas_snrm2 cblas_dnrm2 cblas_scnrm2 cblas_dznrm2
        snrm2 dnrm2 cnrm2 znrm2)

; float cblas_sasum2 (const int N, const float *X, const int incX)
(define-nrm2/asum sasum cblas_sasum 'f32)
(define-nrm2/asum dasum cblas_dasum 'f64)
(define-nrm2/asum casum cblas_scasum 'c32)
(define-nrm2/asum zasum cblas_dzasum 'c64)

(export cblas_sasum cblas_dasum cblas_scasum cblas_dzasum
        sasum dasum casum zasum)


; -----------------------------
; i | max_j('absolute value'(X_j)) = X_i: isamax idamax icamax izamax
; 'absolute value' is |Re|+|Im| in the complex case; cf LawsonEtAl1979,
; 'Basic linear algebra subprograms for Fortran usage, p. 311.
; -----------------------------

(define-syntax define-iamax
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure
                                 int
                                 (dynamic-func cblas-name-string libcblas)
                                 (list int '* int)))
             (define (name X)
               (check-array X 1 srfi4-type)
               (cblas-name (array-length X) (pointer-to-first X) (stride X 0)))))))))

; int cblas_isamax (const int N, const float *X, const int incX)
(define-iamax isamax cblas_isamax 'f32)
(define-iamax idamax cblas_idamax 'f64)
(define-iamax icamax cblas_icamax 'c32)
(define-iamax izamax cblas_izamax 'c64)

(export cblas_isamax cblas_idamax cblas_icamax cblas_izamax
        isamax idamax icamax izamax)


; -----------------------------
; alpha*x_i*(maybe conj)(y_j) + A_{i, j} -> A_{i, j}: sger dger cgeru cgerc zgeru cgerc
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
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure
                                 void
                                 (dynamic-func cblas-name-string libcblas)
                                 (list int int int (srfi4-type->type srfi4-type) '* int '* int '* int)))
             (define (name alpha X Y A)
               (check-arrays-AXY A Y X srfi4-type)
               (let-values (((A-lead A-order) (lead-and-order A)))
                 (cblas-name A-order
                             (dim A 0) (dim A 1) (scalar->arg srfi4-type alpha)
                             (pointer-to-first X) (stride X 0)
                             (pointer-to-first Y) (stride Y 0)
                             (pointer-to-first A) A-lead)))))))))

; void cblas_sger (const enum CBLAS_ORDER Order, const int M, const int N, const float alpha, const float *X,
;                  const int incX, const float *Y, const int incY, float *A, const int lda)
(define-ger sger! cblas_sger 'f32)
(define-ger dger! cblas_dger 'f64)
(define-ger cgeru! cblas_cgeru 'c32)
(define-ger cgerc! cblas_cgerc 'c32)
(define-ger zgeru! cblas_zgeru 'c64)
(define-ger zgerc! cblas_zgerc 'c64)

(export cblas_sger cblas_dger cblas_cgeru cblas_cgerc cblas_zgeru cblas_zgerc
        sger! dger! cgeru! cgerc! zgeru! zgerc!)


; -----------------------------
; alpha*sum_j(A_{ij} * X_j) + beta*Y_i -> Y_i: sgemv dgemv cgemv zgemv
; -----------------------------

; CBLAS_ORDER
(define CblasRowMajor 101)
(define CblasColMajor 102)
(export CblasRowMajor CblasColMajor)

; CBLAS_TRANSPOSE
(define CblasNoTrans 111)
(define CblasTrans 112)
(define CblasConjTrans 113)
(define AtlasConj 114)
(export CblasNoTrans CblasTrans CblasConjTrans AtlasConj)

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

; CBLAS_UPLO
(define CblasUpper 121)
(define CblasLower 122)
(export CblasUpper CblasLower)

; CBLAS_DIAG
(define CblasNonUnit 131)
(define CblasUnit 132)
(export CblasNonUnit CblasUnit)

; CBLAS_SIDE
(define CblasLeft 141)
(define CblasRight 142)
(export CblasLeft CblasRight)

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
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure
                                 void
                                 (dynamic-func cblas-name-string libcblas)
                                 (list int int int int (srfi4-type->type srfi4-type) '* int
                                       '* int (srfi4-type->type srfi4-type) '* int)))
             (define (name alpha A TransA X beta Y)
               (let ((M (dim A 0))
                     (N (dim A 1)))
                 (unless (= M (array-length Y)) (throw 'mismatched-AY M (array-length Y)))
                 (unless (= N (array-length X)) (throw 'mismatched-AX N (array-length X)))
                 (let-values (((A-lead A-order) (lead-and-order A)))
                   (cblas-name A-order TransA
                               M N (scalar->arg srfi4-type alpha)
                               (pointer-to-first A) A-lead
                               (pointer-to-first X) (stride X 0) (scalar->arg srfi4-type beta)
                               (pointer-to-first Y) (stride Y 0)))))))))))

; void cblas_sgemv (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
;                   const int M, const int N, const float alpha, const float *A, const int lda,
;                   const float *X, const int incX, const float beta, float *Y, const int incY)
(define-gemv sgemv! cblas_sgemv 'f32)
(define-gemv dgemv! cblas_dgemv 'f64)
(define-gemv zgemv! cblas_zgemv 'c64)
(define-gemv cgemv! cblas_cgemv 'c32)

(export cblas_dgemv cblas_sgemv cblas_zgemv cblas_cgemv
        dgemv! sgemv! zgemv! cgemv!)


; -----------------------------
; alpha * sum_k(A_{ik}*B_{kj}) + beta * C_{ij} -> C_{ij}: sgemm dgemm cgemm zgemm
; -----------------------------

(define-syntax define-gemm
  (lambda (x)
    (syntax-case x ()
      ((_ name cblas-name srfi4-type)
       (with-syntax ((cblas-name-string (symbol->string (syntax->datum #'cblas-name))))
         #'(begin
             (define cblas-name (pointer->procedure
                                 void
                                 (dynamic-func cblas-name-string libcblas)
                                 (list int int int int int int
                                       (srfi4-type->type srfi4-type) '* int '* int
                                       (srfi4-type->type srfi4-type) '* int)))
             (define (name alpha A TransA B TransB beta C)
               (check-array A 2 srfi4-type)
               (check-array B 2 srfi4-type)
               (check-array C 2 srfi4-type)
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
                                 (scalar->arg srfi4-type alpha)
                                 (pointer-to-first A) A-lead
                                 (pointer-to-first B) B-lead
                                 (scalar->arg srfi4-type beta)
                                 (pointer-to-first C) C-lead)))))))))))

; void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
;                  const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
;                  const int K, const float alpha, const float *A,
;                  const int lda, const float *B, const int ldb,
;                  const float beta, float *C, const int ldc)
(define-gemm sgemm! cblas_sgemm 'f32)
(define-gemm dgemm! cblas_dgemm 'f64)
(define-gemm cgemm! cblas_cgemm 'c32)
(define-gemm zgemm! cblas_zgemm 'c64)

(export cblas_sgemm cblas_dgemm cblas_cgemm cblas_zgemm
        sgemm! dgemm! cgemm! zgemm!)
