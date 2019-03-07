
; Access BLIS through Guile's FFI.
; (c) Daniel Llorens - 2014-2015, 2019

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(define-module (ffi blis))
(import (system foreign) (srfi srfi-1) (srfi srfi-11) (ffi arrays) (ice-9 match) (srfi srfi-26))

; TODO As an alternative go through installation.
(define libblis (dynamic-link (let ((lpath (getenv "GUILE_FFI_CBLAS_LIBBLIS_PATH"))
                                     (lname (or (getenv "GUILE_FFI_CBLAS_LIBBLIS_NAME") "libblis")))
                                 (if (and lpath (not (string=? lpath "")))
                                   (string-append lpath file-name-separator-string lname)
                                   lname))))
(dynamic-call "bli_init" libblis)

(define gint_t int64)
(define dim_t gint_t)
(define inc_t gint_t)
(define doff_t gint_t)
(define trans_t int)
(define conj_t int)
(define rank0_t '*)
(define rank1_t '*)
(define rank2_t '*)


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
; https://github.com/flame/blis/blob/master/frame/include/bli_type_defs.h

(define BLIS_NO_TRANSPOSE 0)
(define BLIS_TRANSPOSE 8)
(define BLIS_CONJ_NO_TRANSPOSE 16)
(define BLIS_CONJ_TRANSPOSE 24)

(define BLIS_NO_CONJUGATE 0)
(define BLIS_CONJUGATE (ash 1 4))

(define BLIS_NONUNIT_DIAG 0)
(define BLIS_UNIT_DIAG (ash 1 8))

(define BLIS_ZEROS 0)
(define BLIS_LOWER (logior (ash 1 7) (ash 1 6)))
(define BLIS_UPPER (logior (ash 1 5) (ash 1 6)))
(define BLIS_DENSE (ash 1 5))

(define BLIS_LEFT 0)
(define BLIS_RIGHT 1)

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

(export BLIS_NO_TRANSPOSE BLIS_TRANSPOSE BLIS_CONJ_NO_TRANSPOSE BLIS_CONJ_TRANSPOSE tr? fliptr
        BLIS_NO_CONJUGATE BLIS_CONJUGATE
        BLIS_NONUNIT_DIAG BLIS_UNIT_DIAG
        BLIS_ZEROS BLIS_LOWER BLIS_UPPER BLIS_DENSE
        BLIS_LEFT BLIS_RIGHT)


; -----------------------------
; level-1v: addv amaxv *axpyv *axpbyv copyv *dotv dotxv invertv scal2v scalv setv subv swapv xpbyv
; -----------------------------

#|
y := beta * y + alpha * conjx(x)

void bli_?axpbyv
     (
       conj_t  conjx,
       dim_t   n,
       ctype*  alpha,
       ctype*  x, inc_t incx,
       ctype*  beta,
       ctype*  y, inc_t incy
     )
|#

(define-syntax define-axpbyv
  (lambda (x)
    (syntax-case x ()
      ((_ type blis-name name!)
       (with-syntax ((type #'(quote type)))
         #`(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func #,(symbol->string (syntax->datum #'blis-name)) libblis)
                                (list conj_t dim_t rank0_t rank1_t inc_t rank0_t rank1_t inc_t)))
             (define (name! conjX alpha X beta Y)
               #,(let ((t (syntax->datum #'type_)))
                   (format #f "(~a conjx [conj_t] alpha [~a] x [#~a(…)] beta [~a] y [#~a(…)])\n\n~a\n"
                           (symbol->string (syntax->datum #'name!)) t t t t
                           "y := beta * y + alpha * conjx(x)"))
               (check-2-arrays X Y 1 type)
               (blis-name conjX (array-length X)
                          (scalar->arg type alpha)
                          (pointer-to-first X) (stride X 0)
                          (scalar->arg type beta)
                          (pointer-to-first Y) (stride Y 0)))))))))

(define-sdcz axpbyv bli_?axpbyv ?axpbyv!)
(define-auto (axpbyv! conjX alpha X beta Y) X ?axpbyv!)

#|
y := y + alpha * conjx(x)

void bli_?axpyv
     (
       conj_t  conjx,
       dim_t   n,
       ctype*  alpha,
       ctype*  x, inc_t incx,
       ctype*  y, inc_t incy
     )
|#

(define-syntax define-axpyv
  (lambda (x)
    (syntax-case x ()
      ((_ type_ blis-name name!)
       (with-syntax ((type #'(quote type_)))
         #`(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func (symbol->string (syntax->datum #'blis-name)) libblis)
                                (list conj_t dim_t rank0_t rank1_t inc_t rank1_t inc_t)))
             (define (name! conjX alpha X Y)
               #,(let ((t (syntax->datum #'type_)))
                   (format #f "(~a conjx [conj_t] alpha [~a] X [#~a(…)] Y [#~a(…)])\n\n~a\n"
                           (symbol->string (syntax->datum #'name!)) t t t
                           "y := y + alpha * conjx(x)"))
               (check-2-arrays X Y 1 type)
               (blis-name conjX (array-length X)
                          (scalar->arg type alpha)
                          (pointer-to-first X) (stride X 0)
                          (pointer-to-first Y) (stride Y 0)))))))))

(define-sdcz axpyv bli_?axpyv ?axpyv!)
(define-auto (axpyv! conjX alpha X Y) X ?axpyv!)

#|
rho := conjx(x)^T * conjy(y)

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

(define-syntax define-dotv
  (lambda (x)
    (syntax-case x ()
      ((_ type_ blis-name name)
       (with-syntax ((type #'(quote type_)))
         #`(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func #,(symbol->string (syntax->datum #'blis-name)) libblis)
                                (list conj_t conj_t dim_t rank1_t inc_t rank1_t inc_t rank0_t)))
             (define (name conjX conjY X Y)
               #,(let ((t (syntax->datum #'type_)))
                   (format #f "(~a conjx [conj_t] conjy [conj_t] x [#~a(…)] y [#~a(…)])\n\t-> rho [~a]\n\n~a\n"
                           (symbol->string (syntax->datum #'name)) t t t
                           "rho := conjx(x)^T * conjy(y)"))
               (check-2-arrays X Y 1 type)
               (let ((rho (make-typed-array type 0)))
                 (blis-name conjX conjY (array-length X)
                            (pointer-to-first X) (stride X 0)
                            (pointer-to-first Y) (stride Y 0)
                            (pointer-to-first rho))
                 (array-ref rho)))))))))

(define-sdcz dotv bli_?dotv ?dotv)
(define-auto (dotv conjX conjY X Y) X ?dotv)


; -----------------------------
; gemm: alpha * sum_k(A_{ik}*B_{kj}) + beta * C_{ij} -> C_{ij}
; -----------------------------

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

(define-syntax define-gemm
  (lambda (x)
    (syntax-case x ()
      ((_ type blis-name name! name)
       (with-syntax ((type #'(quote type)))
         #`(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func #,(symbol->string (syntax->datum #'blis-name)) libblis)
                                (list trans_t trans_t dim_t dim_t dim_t
                                      rank0_t rank2_t inc_t inc_t
                                      rank2_t inc_t inc_t
                                      rank0_t rank2_t inc_t inc_t)))
             (define (name! transA transB alpha A B beta C)
               (check-array A 2 type)
               (check-array B 2 type)
               (check-array C 2 type)
               (let ((M (dim C 0))
                     (N (dim C 1))
                     (K (dim A (if (tr? transA) 0 1))))
                 (unless (= M (dim A (if (tr? transA) 1 0))) (throw 'mismatched-CA))
                 (unless (= N (dim B (if (tr? transB) 0 1))) (throw 'mismatched-CB))
                 (unless (= K (dim B (if (tr? transB) 1 0))) (throw 'mismatched-AB))
                 (blis-name transA transB M N K
                            (scalar->arg type alpha)
                            (pointer-to-first A) (stride A 0) (stride A 1)
                            (pointer-to-first B) (stride B 0) (stride B 1)
                            (scalar->arg type beta)
                            (pointer-to-first C) (stride C 0) (stride C 1))))
             (define (name transA transB alpha A B)
               (let ((C (make-typed-array type *unspecified*
                                          (dim A (if (tr? transA) 1 0))
                                          (dim B (if (tr? transB) 0 1)))))
                 (name! transA transB alpha A B 0. C)
                 C))))))))

(define-sdcz gemm bli_?gemm ?gemm! ?gemm)
(define-auto (gemm! transA transB alpha A B beta C) A ?gemm!)


; -----------------------------
; gemv: alpha*sum_j(A_{ij} * X_j) + beta*Y_i -> Y_i
; -----------------------------

;; void bli_?gemv( trans_t transa,
;;                 conj_t  conjx,
;;                 dim_t   m,
;;                 dim_t   n,
;;                 ctype*  alpha,
;;                 ctype*  a, inc_t rsa, inc_t csa,
;;                 ctype*  x, inc_t incx,
;;                 ctype*  beta,
;;                 ctype*  y, inc_t incy );

(define-syntax define-gemv
  (lambda (x)
    (syntax-case x ()
      ((_ type blis-name name! name)
       (with-syntax ((type #'(quote type)))
         #`(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func #,(symbol->string (syntax->datum #'blis-name)) libblis)
                                (list trans_t conj_t dim_t dim_t
                                      rank0_t rank2_t inc_t inc_t
                                      rank1_t inc_t
                                      rank0_t rank1_t inc_t)))
             (define (name! transA conjX alpha A X beta Y)
               (check-array A 2 type)
               (check-array X 1 type)
               (check-array Y 1 type)
               (let ((M (array-length Y))
                     (N (array-length X)))
                 (unless (= M (dim A (if (tr? transA) 1 0))) (throw 'mismatched-YA))
                 (unless (= N (dim A (if (tr? transA) 0 1))) (throw 'mismatched-XA))
                 (blis-name transA conjX M N
                            (scalar->arg type alpha)
                            (pointer-to-first A) (stride A 0) (stride A 1)
                            (pointer-to-first X) (stride X 0)
                            (scalar->arg type beta)
                            (pointer-to-first Y) (stride Y 0))))
             (define (name transA conjX alpha A X)
               (let ((Y (make-typed-array type *unspecified*
                                          (dim A (if (tr? transA) 1 0)))))
                 (name! transA conjX alpha A X 0 Y)
                 Y))))))))

(define-sdcz gemv bli_?gemv ?gemv! ?gemv)
(define-auto (gemv! transA conjX alpha A X beta Y) A ?gemv!)


; -----------------------------
; ger: alpha*x_i*y_j + A_{i, j} -> A_{i, j}
; -----------------------------

;; void bli_?ger( conj_t  conjx,
;;                conj_t  conjy,
;;                dim_t   m,
;;                dim_t   n,
;;                ctype*  alpha,
;;                ctype*  x, inc_t incx,
;;                ctype*  y, inc_t incy,
;;                ctype*  a, inc_t rsa, inc_t csa );

(define-syntax define-ger
  (lambda (x)
    (syntax-case x ()
      ((_ type blis-name name! name)
       (with-syntax ((type #'(quote type)))
         #`(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func #,(symbol->string (syntax->datum #'blis-name)) libblis)
                                (list conj_t conj_t dim_t dim_t
                                      rank0_t rank1_t inc_t rank1_t inc_t
                                      rank2_t inc_t inc_t)))
             (define (name! conjX conjY alpha X Y A)
               (check-array A 2 type)
               (check-array X 1 type)
               (check-array Y 1 type)
               (let ((M (array-length X))
                     (N (array-length Y)))
                 (unless (= M (dim A 0)) (throw 'mismatched-XA))
                 (unless (= N (dim A 1)) (throw 'mismatched-YA))
                 (blis-name conjX conjY (array-length X) (array-length Y)
                            (scalar->arg type alpha)
                            (pointer-to-first X) (stride X 0)
                            (pointer-to-first Y) (stride Y 0)
                            (pointer-to-first A) (stride A 0) (stride A 1))))
             (define (name conjX conjY alpha X Y)
               (let ((A (make-typed-array type 0
                                          (array-length X) (array-length Y))))
                 (name! conjX conjY alpha X Y A)
                 A))))))))

(define-sdcz ger bli_?ger ?ger! ?ger)
(define-auto (ger! conjX conjY alpha X Y A) X ?ger!)
