
; Access BLIS through Guile's FFI.
; (c) Daniel Llorens - 2014-2015, 2019

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(define-module (ffi blis))
(import (system foreign) (srfi srfi-1) (srfi srfi-11) (ffi arrays) (ice-9 match))

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

(eval-when (expand load eval)
  (define (level stx-name t)
    (let* ((s (symbol->string (syntax->datum stx-name)))
           (i (string-index s #\?))
           (fmt (string-replace s "~a" i (+ i 1))))
      (datum->syntax stx-name (string->symbol (format #f fmt t))))))

; FIXME ellipsis n0 ...
; FIXME export from here.
(define-syntax define-sdcz
  (lambda (x)
    (syntax-case x ()
      ((_ definer n0 n1)
       #`(begin
           (definer f32 #,(level #'n0 's) #,(level #'n1 's))
           (definer f64 #,(level #'n0 'd) #,(level #'n1 'd))
           (definer c32 #,(level #'n0 'c) #,(level #'n1 'c))
           (definer c64 #,(level #'n0 'z) #,(level #'n1 'z))))
      ((_ definer n0 n1 n2)
       #`(begin
           (definer f32 #,(level #'n0 's) #,(level #'n1 's) #,(level #'n2 's))
           (definer f64 #,(level #'n0 'd) #,(level #'n1 'd) #,(level #'n2 'd))
           (definer c32 #,(level #'n0 'c) #,(level #'n1 'c) #,(level #'n2 'c))
           (definer c64 #,(level #'n0 'z) #,(level #'n1 'z) #,(level #'n2 'z)))))))


; -----------------------------
; BLIS flags
; -----------------------------

; https://github.com/flame/blis/blob/master/docs/BLISTypedAPI.md

(define BLIS_NO_TRANSPOSE 0)
(define BLIS_TRANSPOSE 8)
(define BLIS_CONJ_NO_TRANSPOSE 16)
(define BLIS_CONJ_TRANSPOSE 24)

(define BLIS_NO_CONJUGATE 0)
(define BLIS_CONJUGATE (ash 1 4))

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
        BLIS_NO_CONJUGATE BLIS_CONJUGATE)


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
       (with-syntax
           ((blis-name-string (symbol->string (syntax->datum #'blis-name)))
            (docstring (format #f "(~a conjx [conj_t] alpha [~a] x [#~a(…)] beta [~a] y [#~a(…)])\n\n~a\n"
                               (symbol->string (syntax->datum #'name!))
                               (syntax->datum #'type) (syntax->datum #'type)
                               (syntax->datum #'type) (syntax->datum #'type)
                               "y := beta * y + alpha * conjx(x)")))
         #'(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func blis-name-string libblis)
                                (list conj_t dim_t '* '* inc_t '* '* inc_t)))
             (define (name! conjX alpha X beta Y)
               docstring
               (check-2-arrays X Y 1 (quote type))
               (blis-name conjX (array-length X)
                          (scalar->arg (quote type) alpha)
                          (pointer-to-first X) (stride X 0)
                          (scalar->arg (quote type) beta)
                          (pointer-to-first Y) (stride Y 0)))))))))

(define-sdcz define-axpbyv bli_?axpbyv ?axpbyv!)

(define (axpbyv! conjX alpha X beta Y)
  "
(axpbyv! conjX [conj_t] alpha [type] X [#type(…)] beta [type] Y [#type(…)])

y := beta * y + alpha * conjx(x)

See also: saxpbyv! daxpbyv! caxpbyv! zaxpbyv! axpyv!"

  ((match (array-type X)
     (('f32 saxpbyv!) ('f64 daxpbyv!) ('c32 caxpbyv!) ('c64 saxpbyv!)))
   conjX alpha X beta Y))

(export bli_saxpbyv bli_daxpbyv bli_caxpbyv bli_zaxpbyv
        saxpbyv! daxpbyv! caxpbyv! zaxpbyv!
        axpbyv!)

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
      ((_ type blis-name name!)
       (with-syntax
           ((blis-name-string (symbol->string (syntax->datum #'blis-name)))
            (docstring (format #f "(~a conjx [conj_t] alpha [~a] x [#~a(…)] y [#~a(…)])\n\n~a\n"
                               (symbol->string (syntax->datum #'name!))
                               (syntax->datum #'type) (syntax->datum #'type)
                               (syntax->datum #'type)
                               "y := y + alpha * conjx(x)")))
         #'(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func blis-name-string libblis)
                                (list conj_t dim_t '* '* inc_t '* inc_t)))
             (define (name! conjX alpha X Y)
               docstring
               (check-2-arrays X Y 1 (quote type))
               (blis-name conjX (array-length X)
                            (scalar->arg (quote type) alpha)
                            (pointer-to-first X) (stride X 0)
                            (pointer-to-first Y) (stride Y 0)))
             (export blis-name name!)))))))

(define-sdcz define-axpyv bli_?axpyv ?axpyv!)

(define (axpyv! conjX alpha X Y)
  "
(axpyv! conjX [conj_t] alpha [type] X [#type()] Y [#type()])

y := y + alpha * conjx(x)

See also: saxpyv! daxpyv! caxpyv! zaxpyv! axpbyv!"

  ((match (array-type X)
     ('f32 saxpyv!) ('f64 daxpyv!) ('c32 caxpyv!) ('c64 zaxpyv!))
   conjX alpha X Y))

(export axpyv!)

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
      ((_ type blis-name name)
       (with-syntax
           ((blis-name-string (symbol->string (syntax->datum #'blis-name)))
            (docstring (format #f "(~a conjx [conj_t] conjy [conj_t] x [#~a(…)] y [#~a(…)])\n\t-> rho [~a]\n\n~a\n"
                               (symbol->string (syntax->datum #'name))
                               (syntax->datum #'type) (syntax->datum #'type)
                               (syntax->datum #'type)
                               "rho := conjx(x)^T * conjy(y)")))
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
                 (array-ref rho)))
             (export blis-name name)))))))

(define-sdcz define-dotv bli_?dotv ?dotv)

(define (dotv conjX conjY X Y)
  "
(dotv conjX [conj_t] conjY [conj_t] X [#type(…)] Y [#type(…)]\n\t-> rho [type])

rho := conjx(x)^T * conjy(y)

See also: sdotv ddotv cdotv vdotv"

  ((match (array-type X)
     ('f32 sdotv) ('f64 ddotv) ('c32 cdotv) ('c64 zdotv))
   conjX conjY X Y))

(export dotv)


; -----------------------------
; gemm: alpha * sum_k(A_{ik}*B_{kj}) + beta * C_{ij} -> C_{ij}
; -----------------------------

(define-syntax define-gemm
  (lambda (x)
    (syntax-case x ()
      ((_ type blis-name name! name)
       (with-syntax
           ((blis-name-string (symbol->string (syntax->datum #'blis-name))))
         #'(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func blis-name-string libblis)
                                (list trans_t trans_t dim_t dim_t dim_t
                                      '* '* inc_t inc_t '* inc_t inc_t
                                      '* '* inc_t inc_t)))
             (define (name! alpha A transA B transB beta C)
               (check-array A 2 (quote type))
               (check-array B 2 (quote type))
               (check-array C 2 (quote type))
               (let ((M (dim C 0))
                     (N (dim C 1))
                     (K (dim A (if (tr? transA) 0 1))))
                 (unless (= M (dim A (if (tr? transA) 1 0))) (throw 'mismatched-CA))
                 (unless (= N (dim B (if (tr? transB) 0 1))) (throw 'mismatched-CB))
                 (unless (= K (dim B (if (tr? transB) 1 0))) (throw 'mismatched-AB))
                 (blis-name transA transB M N K
                            (scalar->arg (quote type) alpha)
                            (pointer-to-first A) (stride A 0) (stride A 1)
                            (pointer-to-first B) (stride B 0) (stride B 1)
                            (scalar->arg (quote type) beta)
                            (pointer-to-first C) (stride C 0) (stride C 1))))
             (define (name alpha A transA B transB)
               (let ((C (make-typed-array (quote type) *unspecified*
                                          (dim A (if (tr? transA) 1 0))
                                          (dim B (if (tr? transB) 0 1)))))
                 (name! alpha A transA B transB 0. C)
                 C))
             (export blis-name name! name)))))))

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

(define-sdcz define-gemm bli_?gemm ?gemm! ?gemm)


; -----------------------------
; gemv: alpha*sum_j(A_{ij} * X_j) + beta*Y_i -> Y_i
; -----------------------------

(define-syntax define-gemv
  (lambda (x)
    (syntax-case x ()
      ((_ type blis-name name! name)
       (with-syntax
           ((blis-name-string (symbol->string (syntax->datum #'blis-name))))
         #'(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func blis-name-string libblis)
                                (list trans_t conj_t dim_t dim_t '* '* inc_t inc_t
                                      '* inc_t '* '* inc_t)))
             (define (name! alpha A transA X conjX beta Y)
               (check-array A 2 (quote type))
               (check-array X 1 (quote type))
               (check-array Y 1 (quote type))
               (let ((M (array-length Y))
                     (N (array-length X)))
                 (unless (= M (dim A (if (tr? transA) 1 0))) (throw 'mismatched-YA))
                 (unless (= N (dim A (if (tr? transA) 0 1))) (throw 'mismatched-XA))
                 (blis-name transA conjX M N
                            (scalar->arg (quote type) alpha)
                            (pointer-to-first A) (stride A 0) (stride A 1)
                            (pointer-to-first X) (stride X 0)
                            (scalar->arg (quote type) beta)
                            (pointer-to-first Y) (stride Y 0))))
             (define (name alpha A transA X conjX)
               (let ((Y (make-typed-array (quote type) *unspecified*
                                          (dim A (if (tr? transA) 1 0)))))
                 (name! alpha A transA X conjX 0. Y)
                 Y))
             (export blis-name name! name)))))))

;; void bli_?gemv( trans_t transa,
;;                 conj_t  conjx,
;;                 dim_t   m,
;;                 dim_t   n,
;;                 ctype*  alpha,
;;                 ctype*  a, inc_t rsa, inc_t csa,
;;                 ctype*  x, inc_t incx,
;;                 ctype*  beta,
;;                 ctype*  y, inc_t incy );

(define-sdcz define-gemv bli_?gemv ?gemv! ?gemv)


; -----------------------------
; ger: alpha*x_i*y_j + A_{i, j} -> A_{i, j}
; -----------------------------

(define-syntax define-ger
  (lambda (x)
    (syntax-case x ()
      ((_ type blis-name name! name)
       (with-syntax
           ((blis-name-string (symbol->string (syntax->datum #'blis-name))))
         #'(begin
             (define blis-name (pointer->procedure
                                void (dynamic-func blis-name-string libblis)
                                (list conj_t conj_t dim_t dim_t '* '* inc_t '* inc_t
                                      '* inc_t inc_t)))
             (define (name! alpha X conjX Y conjY A)
               (check-array A 2 (quote type))
               (check-array X 1 (quote type))
               (check-array Y 1 (quote type))
               (let ((M (array-length X))
                     (N (array-length Y)))
                 (unless (= M (dim A 0)) (throw 'mismatched-XA))
                 (unless (= N (dim A 1)) (throw 'mismatched-YA))
                 (blis-name conjX conjY (array-length X) (array-length Y)
                            (scalar->arg (quote type) alpha)
                            (pointer-to-first X) (stride X 0)
                            (pointer-to-first Y) (stride Y 0)
                            (pointer-to-first A) (stride A 0) (stride A 1))))
             (define (name alpha X conjX Y conjY)
               (let ((A (make-typed-array (quote type) *unspecified*
                                          (array-length X) (array-length Y))))
                 (name! alpha X conjX Y conjY A)
                 A))
             (export blis-name name! name)))))))

;; void bli_?ger( conj_t  conjx,
;;                conj_t  conjy,
;;                dim_t   m,
;;                dim_t   n,
;;                ctype*  alpha,
;;                ctype*  x, inc_t incx,
;;                ctype*  y, inc_t incy,
;;                ctype*  a, inc_t rsa, inc_t csa );

(define-sdcz define-ger bli_?ger ?ger! ?ger)
