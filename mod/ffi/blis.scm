
; Access BLIS through Guile's FFI.
; (c) Daniel Llorens - 2014-2015

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(define-module (ffi blis))
(import (system foreign) (srfi srfi-1) (srfi srfi-11))

; https://code.google.com/p/blis/wiki/BLISAPIQuickReference

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

(define (pointer-to-first A)
  (bytevector->pointer (shared-array-root A)
                       (* (shared-array-offset A) (srfi4-type-size (array-type A)))))

(define (scalar->blis-arg srfi4-type a)
  (pointer-to-first (make-typed-array srfi4-type a)))

(define (stride A i)
  (list-ref (shared-array-increments A) i))

(define (dim A i)
  (list-ref (array-dimensions A) i))


; -----------------------------
; BLIS flags
; -----------------------------

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
; alpha * sum_k(A_{ik}*B_{kj}) + beta * C_{ij} -> C_{ij}
; -----------------------------

(define-syntax define-gemm
  (lambda (x)
    (syntax-case x ()
      ((_ name! name blis-name srfi4-type)
       (with-syntax ((blis-name-string (symbol->string (syntax->datum #'blis-name))))
         #'(begin
             (define blis-name (pointer->procedure
                                void
                                (dynamic-func blis-name-string libblis)
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
                            (scalar->blis-arg srfi4-type alpha)
                            (pointer-to-first A) (stride A 0) (stride A 1)
                            (pointer-to-first B) (stride B 0) (stride B 1)
                            (scalar->blis-arg srfi4-type beta)
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
; alpha*sum_j(A_{ij} * X_j) + beta*Y_i -> Y_i
; -----------------------------

(define-syntax define-gemv
  (lambda (x)
    (syntax-case x ()
      ((_ name! name blis-name srfi4-type)
       (with-syntax ((blis-name-string (symbol->string (syntax->datum #'blis-name))))
         #'(begin
             (define blis-name (pointer->procedure
                                void
                                (dynamic-func blis-name-string libblis)
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
                            (scalar->blis-arg srfi4-type alpha)
                            (pointer-to-first A) (stride A 0) (stride A 1)
                            (pointer-to-first X) (stride X 0)
                            (scalar->blis-arg srfi4-type beta)
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
; alpha*x_i*y_j + A_{i, j} -> A_{i, j}
; -----------------------------

(define-syntax define-ger
  (lambda (x)
    (syntax-case x ()
      ((_ name! name blis-name srfi4-type)
       (with-syntax ((blis-name-string (symbol->string (syntax->datum #'blis-name))))
         #'(begin
             (define blis-name (pointer->procedure
                                void
                                (dynamic-func blis-name-string libblis)
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
                            (scalar->blis-arg srfi4-type alpha)
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
