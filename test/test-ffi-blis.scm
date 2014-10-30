
; Tests for (ffi blis).
; (c) Daniel Llorens - 2014

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(import (ffi blis) (srfi srfi-64) (srfi srfi-1) (ice-9 match))

; ---------------------------------
; Various sorts of arrays.
; ---------------------------------

(define (make-v-compact type)
  (make-typed-array type *unspecified* 10))

(define (make-v-strided type)
  (make-shared-array (make-typed-array type *unspecified* 20) (lambda (i) (list (* i 2))) 10))

(define (make-v-offset type)
  (make-shared-array (make-typed-array type *unspecified* 12) (lambda (i) (list (+ i 2))) 10))

(define (fill-A1! A)
  (case (array-type A)
    ((f32 f64) (array-index-map! A (lambda (i) (- i 4))))
    ((c32 c64) (array-index-map! A (lambda (i) (make-rectangular (+ i 1) (- i 5)))))
    (else (throw 'bad-array-type (array-type A))))
  A)

(define (fill-B1! A)
  (case (array-type A)
    ((f32 f64) (array-index-map! A (lambda (i) (+ i 3))))
    ((c32 c64) (array-index-map! A (lambda (i) (make-rectangular (+ i 2) (- 5 i)))))
    (else (throw 'bad-array-type (array-type A))))
  A)

(define (make-M-c-order type)
  (make-typed-array type *unspecified* 10 10))

(define (make-M-fortran-order type)
  (transpose-array (make-typed-array type *unspecified* 10 10) 1 0))

(define (make-M-strided type)
  (make-shared-array (make-typed-array type *unspecified* 20 10) (lambda (i j) (list (* i 2) j)) 10 10))

(define (make-M-strided-both type)
  (make-shared-array (make-typed-array type *unspecified* 20 20) (lambda (i j) (list (* i 2) (* j 2))) 10 10))

(define (make-M-offset type)
  (make-shared-array (make-typed-array type *unspecified* 12 13) (lambda (i j) (list (+ i 2) (+ j 3))) 10 10))

(define (fill-A2! A)
  (case (array-type A)
    ((f32 f64) (array-index-map!
                A (lambda (i j) (+ 4 (* i 1) (* j j 2)))))
    ((c32 c64) (array-index-map!
                A (lambda (i j) (make-rectangular (+ 4 (* i 1) (* j j 2)) (+ -3 (* i i 1) (* j 2) -4)))))
    (else (throw 'bad-array-type (array-type A))))
  A)

; ---------------------------------
; Utilities
; ---------------------------------

(define (conj a) (make-rectangular (real-part a) (- (imag-part a))))

(define (array-copy A)
  (let ((B (apply make-typed-array (array-type A) *unspecified* (array-dimensions A))))
    (array-copy! A B)
    B))

; ---------------------------------
; Test types
; ---------------------------------

(define* (test-approximate-array tag test expected err)
  (test-begin tag)
  (array-for-each (lambda (test expected) (test-approximate test expected err))
                  test expected)
  (test-end tag))

; ---------------------------------
; @TODO sgemm dgemm cgemm zgemm
; ---------------------------------

(define (apply-transpose-flag A TransA)
  (cond ((equal? TransA BLIS_NO_TRANSPOSE) A)
        ((equal? TransA BLIS_TRANSPOSE) (transpose-array A 1 0))
        ((equal? TransA BLIS_CONJ_NO_TRANSPOSE) (let ((B (array-copy A))) (array-map! B conj A) B))
        ((equal? TransA BLIS_CONJ_TRANSPOSE) (let ((B (array-copy A))) (array-map! B conj A) (transpose-array B 1 0)))
        (else (throw 'bad-transpose-flag TransA))))

; alpha * sum_k(A_{ik}*B_{kj}) + beta * C_{ij} -> C_{ij}
(define (ref-gemm! alpha A TransA B TransB beta C)
  (let* ((A (apply-transpose-flag A TransA))
         (B (apply-transpose-flag B TransB))
         (M (first (array-dimensions C)))
         (N (second (array-dimensions C)))
         (K (first (array-dimensions B))))
     (do ((i 0 (+ i 1))) ((= i M))
       (do ((j 0 (+ j 1))) ((= j N))
         (array-set! C (* beta (array-ref C i j)) i j)
         (do ((k 0 (+ k 1))) ((= k K))
           (array-set! C (+ (array-ref C i j) (* alpha (array-ref A i k) (array-ref B k j))) i j))))))

(define (test-gemm tag gemm! alpha A TransA B TransB beta C)
  (let ((C1 (array-copy C))
        (C2 (array-copy C))
        (AA (array-copy A))
        (BB (array-copy B)))
    (gemm! alpha A TransA B TransB beta C1)
    (ref-gemm! alpha A TransA B TransB beta C2)
    ;; (test-approximate-array tag C1 C2 1e-15) ; @TODO  as a single test.
    (test-begin tag)
    (test-equal C1 C2)
    (test-end tag)))

(map
 (match-lambda
     ((srfi4-type gemm!)
; some extra tests with non-square matrices.
      (let ((A (fill-A2! (make-typed-array srfi4-type *unspecified* 4 3)))
            (B (fill-A2! (make-typed-array srfi4-type *unspecified* 3 5)))
            (C (fill-A2! (make-typed-array srfi4-type *unspecified* 4 5))))
        (test-gemm "gemm-1" gemm! 1. A BLIS_NO_TRANSPOSE B BLIS_NO_TRANSPOSE 1. C)
        (test-gemm "gemm-2" gemm! 1. A BLIS_TRANSPOSE C BLIS_NO_TRANSPOSE 1. B)
        (test-gemm "gemm-3" gemm! 1. C BLIS_NO_TRANSPOSE B BLIS_TRANSPOSE 1. A))
      (let ((A (fill-A2! (transpose-array (make-typed-array 'f64 *unspecified* 4 3) 1 0)))
            (B (fill-A2! (transpose-array (make-typed-array 'f64 *unspecified* 3 5) 1 0)))
            (C (fill-A2! (transpose-array (make-typed-array 'f64 *unspecified* 4 5) 1 0))))
        (test-gemm "gemm-4" dgemm! 1. A BLIS_TRANSPOSE B BLIS_TRANSPOSE 1. (transpose-array C 1 0))
        (test-gemm "gemm-5" dgemm! 1. A BLIS_NO_TRANSPOSE C BLIS_TRANSPOSE 1. (transpose-array B 1 0))
        (test-gemm "gemm-6" dgemm! 1. C BLIS_TRANSPOSE B BLIS_NO_TRANSPOSE 1. (transpose-array A 1 0)))
      (for-each
       (lambda (make-A)
         (for-each
          (lambda (make-B)
            (for-each
             (lambda (make-C)
               (for-each
                (lambda (TransA)
                  (for-each
                   (lambda (TransB)
                     (test-gemm (format #f "gemm:~a:~a:~a:~a:~a:~a" srfi4-type (procedure-name make-A)
                                        (procedure-name make-B) (procedure-name make-C)
                                        TransA TransB)
                                gemm! 3. (fill-A2! (make-A srfi4-type)) TransA
                                (fill-A2! (make-B srfi4-type)) TransB
                                2. (fill-A2! (make-C srfi4-type))))
                   (list BLIS_TRANSPOSE BLIS_NO_TRANSPOSE BLIS_CONJ_NO_TRANSPOSE BLIS_CONJ_TRANSPOSE)))
                (list BLIS_TRANSPOSE BLIS_NO_TRANSPOSE BLIS_CONJ_NO_TRANSPOSE BLIS_CONJ_TRANSPOSE)))
             (list make-M-c-order make-M-fortran-order make-M-offset make-M-strided make-M-strided-both)))
          (list make-M-c-order make-M-fortran-order make-M-offset make-M-strided make-M-strided-both)))
       (list make-M-c-order make-M-fortran-order make-M-offset make-M-strided make-M-strided-both))))
 `((f32 ,sgemm!)
   (f64 ,dgemm!)
   (c32 ,cgemm!)
   (c64 ,zgemm!)))

(unless (zero? (test-runner-fail-count (test-runner-current)))
  (error "FAILED test-ffi-cblas.csm"))
