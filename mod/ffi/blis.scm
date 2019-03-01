
; Access BLIS through Guile's FFI.
; (c) Daniel Llorens - 2014-2015, 2019

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(define-module (ffi blis))
(import (system foreign) (srfi srfi-1) (srfi srfi-11) (ffi arrays))

; TODO As an alternative go through installation.
(define libblis (dynamic-link (let ((lpath (getenv "GUILE_FFI_CBLAS_LIBBLIS_PATH"))
                                     (lname (or (getenv "GUILE_FFI_CBLAS_LIBBLIS_NAME") "libblis")))
                                 (if (and lpath (not (string=? lpath "")))
                                   (string-append lpath file-name-separator-string lname)
                                   lname))))
(dynamic-call "bli_init" libblis)

(define dim_t int64)
(define inc_t int64)
(define trans_t int)
(define conj_t int)


; -----------------------------
; wrapper utilities
; -----------------------------

(define (pointer-to-first A)
  (bytevector->pointer (shared-array-root A)
                       (* (shared-array-offset A) (srfi4-type-size (array-type A)))))

(define (scalar->arg srfi4-type a)
  (pointer-to-first (make-typed-array srfi4-type a)))


; -----------------------------
; BLIS flags
; -----------------------------

; https://github.com/flame/blis/blob/master/docs/BLISTypedAPI.md

(define BLIS_NO_TRANSPOSE 0)
(define BLIS_TRANSPOSE 8)
(define BLIS_CONJ_NO_TRANSPOSE 16)
(define BLIS_CONJ_TRANSPOSE 24)

(define (fliptr t)
  (cond
   ((= t BLIS_NO_TRANSPOSE) BLIS_TRANSPOSE)
   ((= t BLIS_TRANSPOSE) BLIS_NO_TRANSPOSE)
   ((= t BLIS_CONJ_NO_TRANSPOSE) BLIS_CONJ_TRANSPOSE)
   ((= t BLIS_CONJ_TRANSPOSE) BLIS_CONJ_NO_TRANSPOSE)
   (else (throw 'bad-transpose-1 t))))

(define (tr? t)
  (cond
   ((= t BLIS_NO_TRANSPOSE) #f)
   ((= t BLIS_TRANSPOSE) #t)
   ((= t BLIS_CONJ_NO_TRANSPOSE) #f)
   ((= t BLIS_CONJ_TRANSPOSE) #t)
   (else (throw 'bad-transpose-2 t))))

(export BLIS_NO_TRANSPOSE BLIS_TRANSPOSE BLIS_CONJ_NO_TRANSPOSE BLIS_CONJ_TRANSPOSE tr? fliptr)

(define BLIS_NO_CONJUGATE 0)
(define BLIS_CONJUGATE (ash 1 4))

(export BLIS_NO_CONJUGATE BLIS_CONJUGATE)


; -----------------------------
; level-1v: addv amaxv axpyv axpbyv copyv *dotv dotxv invertv scal2v scalv setv subv swapv xpbyv
; -----------------------------

(define-syntax define-dotv
  (lambda (x)
    (syntax-case x ()
      ((_ name blis-name type)
       (with-syntax ((blis-name-string (symbol->string (syntax->datum #'blis-name)))
                     (docstring (format #f "~a conjx [conj_t] conjy [conj_t] x [#~a(…)] y [#~a(…)]  -> rho"
                                        (symbol->string (syntax->datum #'name))
                                        (syntax->datum #'type) (syntax->datum #'type))))
         #'(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func blis-name-string libblis)
                                (list conj_t conj_t dim_t '* inc_t '* inc_t '*)))
             (define (name conjX conjY X Y)
               docstring
               (check-2-arrays X Y 1 (quote type))
               (let ((rho (make-typed-array (quote type) 0)))
                 (blis-name conjX conjY (array-length X)
                            (pointer-to-first X) (stride X 0)
                            (pointer-to-first Y) (stride Y 0)
                            (pointer-to-first rho))
                 (array-ref rho)))))))))

#|
void bli_?dotv
     (
       conj_t  conjx,
       conj_t  conjy,
       dim_t   n,
       ctype*  x, inc_t incx,
       ctype*  y, inc_t incy,
       ctype*  rho
       );
|#

(define-dotv sdotv bli_sdotv f32)
(define-dotv ddotv bli_ddotv f64)
(define-dotv cdotv bli_cdotv c32)
(define-dotv zdotv bli_zdotv c64)

(export bli_sdotv bli_ddotv bli_cdotv bli_zdotv
        sdotv ddotv cdotv zdotv)


; -----------------------------
; gemm: alpha * sum_k(A_{ik}*B_{kj}) + beta * C_{ij} -> C_{ij}
; -----------------------------

(define-syntax define-gemm
  (lambda (x)
    (syntax-case x ()
      ((_ name! name blis-name srfi4-type)
       (with-syntax ((blis-name-string (symbol->string (syntax->datum #'blis-name))))
         #'(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func blis-name-string libblis)
                                (list trans_t trans_t dim_t dim_t dim_t
                                      '* '* inc_t inc_t '* inc_t inc_t
                                      '* '* inc_t inc_t)))
             (define (name! alpha A transA B transB beta C)
               (check-array A 2 srfi4-type)
               (check-array B 2 srfi4-type)
               (check-array C 2 srfi4-type)
               (let ((M (dim C 0))
                     (N (dim C 1))
                     (K (dim A (if (tr? transA) 0 1))))
                 (unless (= M (dim A (if (tr? transA) 1 0))) (throw 'mismatched-CA))
                 (unless (= N (dim B (if (tr? transB) 0 1))) (throw 'mismatched-CB))
                 (unless (= K (dim B (if (tr? transB) 1 0))) (throw 'mismatched-AB))
                 (blis-name transA transB M N K
                            (scalar->arg srfi4-type alpha)
                            (pointer-to-first A) (stride A 0) (stride A 1)
                            (pointer-to-first B) (stride B 0) (stride B 1)
                            (scalar->arg srfi4-type beta)
                            (pointer-to-first C) (stride C 0) (stride C 1))))
             (define (name alpha A transA B transB)
               (let ((C (make-typed-array srfi4-type *unspecified*
                                          (dim A (if (tr? transA) 1 0))
                                          (dim B (if (tr? transB) 0 1)))))
                 (name! alpha A transA B transB 0. C)
                 C))))))))

;; void bli_?gemm( trans_t transa,
;;                 trans_t transb,
;;                 dim_t   m,
;;                 dim_t   n,
;;                 dim_t   k,
;;                 ctype*  alpha,
;;                 ctype*  a, inc_t rsa, inc_t csa,
;;                 ctype*  b, inc_t rsb, inc_t csb,
;;                 ctype*  beta,
;;                 ctype*  c, inc_t rsc, inc_t csc )

(define-gemm sgemm! sgemm bli_sgemm 'f32)
(define-gemm dgemm! dgemm bli_dgemm 'f64)
(define-gemm cgemm! cgemm bli_cgemm 'c32)
(define-gemm zgemm! zgemm bli_zgemm 'c64)

(export bli_sgemm bli_dgemm bli_cgemm bli_zgemm
        sgemm! dgemm! cgemm! zgemm!
        sgemm dgemm cgemm zgemm)


; -----------------------------
; gemv: alpha*sum_j(A_{ij} * X_j) + beta*Y_i -> Y_i
; -----------------------------

(define-syntax define-gemv
  (lambda (x)
    (syntax-case x ()
      ((_ name! name blis-name srfi4-type)
       (with-syntax ((blis-name-string (symbol->string (syntax->datum #'blis-name))))
         #'(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func blis-name-string libblis)
                                (list trans_t conj_t dim_t dim_t '* '* inc_t inc_t
                                      '* inc_t '* '* inc_t)))
             (define (name! alpha A transA X conjX beta Y)
               (check-array A 2 srfi4-type)
               (check-array X 1 srfi4-type)
               (check-array Y 1 srfi4-type)
               (let ((M (array-length Y))
                     (N (array-length X)))
                 (unless (= M (dim A (if (tr? transA) 1 0))) (throw 'mismatched-YA))
                 (unless (= N (dim A (if (tr? transA) 0 1))) (throw 'mismatched-XA))
                 (blis-name transA conjX M N
                            (scalar->arg srfi4-type alpha)
                            (pointer-to-first A) (stride A 0) (stride A 1)
                            (pointer-to-first X) (stride X 0)
                            (scalar->arg srfi4-type beta)
                            (pointer-to-first Y) (stride Y 0))))
             (define (name alpha A transA X conjX)
               (let ((Y (make-typed-array srfi4-type *unspecified*
                                          (dim A (if (tr? transA) 1 0)))))
                 (name! alpha A transA X conjX 0. Y)
                 Y))))))))

;; void bli_?gemv( trans_t transa,
;;                 conj_t  conjx,
;;                 dim_t   m,
;;                 dim_t   n,
;;                 ctype*  alpha,
;;                 ctype*  a, inc_t rsa, inc_t csa,
;;                 ctype*  x, inc_t incx,
;;                 ctype*  beta,
;;                 ctype*  y, inc_t incy );

(define-gemv sgemv! sgemv bli_sgemv 'f32)
(define-gemv dgemv! dgemv bli_dgemv 'f64)
(define-gemv cgemv! cgemv bli_cgemv 'c32)
(define-gemv zgemv! zgemv bli_zgemv 'c64)

(export bli_sgemv bli_dgemv bli_cgemv bli_zgemv
        sgemv! dgemv! cgemv! zgemv!
        sgemv dgemv cgemv zgemv)


; -----------------------------
; ger: alpha*x_i*y_j + A_{i, j} -> A_{i, j}
; -----------------------------

(define-syntax define-ger
  (lambda (x)
    (syntax-case x ()
      ((_ name! name blis-name srfi4-type)
       (with-syntax ((blis-name-string (symbol->string (syntax->datum #'blis-name))))
         #'(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func blis-name-string libblis)
                                (list conj_t conj_t dim_t dim_t '* '* inc_t '* inc_t
                                      '* inc_t inc_t)))
             (define (name! alpha X conjX Y conjY A)
               (check-array A 2 srfi4-type)
               (check-array X 1 srfi4-type)
               (check-array Y 1 srfi4-type)
               (let ((M (array-length X))
                     (N (array-length Y)))
                 (unless (= M (dim A 0)) (throw 'mismatched-XA))
                 (unless (= N (dim A 1)) (throw 'mismatched-YA))
                 (blis-name conjX conjY (array-length X) (array-length Y)
                            (scalar->arg srfi4-type alpha)
                            (pointer-to-first X) (stride X 0)
                            (pointer-to-first Y) (stride Y 0)
                            (pointer-to-first A) (stride A 0) (stride A 1))))
             (define (name alpha X conjX Y conjY)
               (let ((A (make-typed-array srfi4-type *unspecified*
                                          (array-length X) (array-length Y))))
                 (name! alpha X conjX Y conjY A)
                 A))))))))

;; void bli_?ger( conj_t  conjx,
;;                conj_t  conjy,
;;                dim_t   m,
;;                dim_t   n,
;;                ctype*  alpha,
;;                ctype*  x, inc_t incx,
;;                ctype*  y, inc_t incy,
;;                ctype*  a, inc_t rsa, inc_t csa );

(define-ger sger! sger bli_sger 'f32)
(define-ger dger! dger bli_dger 'f64)
(define-ger cger! cger bli_cger 'c32)
(define-ger zger! zger bli_zger 'c64)

(export bli_sger bli_dger bli_cger bli_zger
        sger! dger! cger! zger!
        sger dger cger zger)
