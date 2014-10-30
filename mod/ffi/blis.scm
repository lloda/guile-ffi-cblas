
; Access BLIS through Guile's FFI.
; (c) Daniel Llorens - 2014

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(define-module (ffi blis))
(import (system foreign) (srfi srfi-1) (srfi srfi-11))

; https://code.google.com/p/blis/wiki/BLISAPIQuickReference

; @TODO As an alternative go through installation.
(define libblis (dynamic-link (let ((lpath (getenv "GUILE_FFI_BLIS_LIBBLIS_PATH")))
                                (if (and lpath (not (string=? lpath "")))
                                  (string-append lpath file-name-separator-string "libblis")
                                  "libblis"))))
(dynamic-call "bli_init" libblis)

(define dim_t int64)
(define inc_t int64)

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
; alpha * sum_k(A_{ik}*B_{kj}) + beta * C_{ij} -> C_{ij}: sgemm dgemm cgemm zgemm
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

(define-syntax define-gemm
  (lambda (x)
    (syntax-case x ()
      ((_ name! name blis-name srfi4-type)
       (with-syntax ((blis-name-string (symbol->string (syntax->datum (syntax blis-name)))))
         (syntax
          (begin
            (define blis-name (pointer->procedure
                               void
                               (dynamic-func blis-name-string libblis)
                               (list int int dim_t dim_t dim_t
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
                C)))))))))

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

(export bli_sgemm bli_dgemm bli_cgemm bli_zgemm)
(export sgemm! dgemm! cgemm! zgemm!)
(export sgemm dgemm cgemm zgemm)
