
; Tests for (ffi blis).
; (c) Daniel Llorens - 2014-2015, 2019

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(import (ffi blis) (srfi srfi-64) (srfi srfi-1) (ice-9 match) (srfi srfi-26))
(include "common.scm")

(set! test-log-to-file #f)

(define (apply-transpose-flag A flag)
  (cond ((or (equal? flag BLIS_NO_TRANSPOSE) (equal? flag BLIS_NO_CONJUGATE)) A)
        ((equal? flag BLIS_TRANSPOSE) (transpose-array A 1 0))
        ((or (equal? flag BLIS_CONJ_NO_TRANSPOSE) (equal? flag BLIS_CONJUGATE))
         (let ((B (array-copy A))) (array-map! B conj A) B))
        ((equal? flag BLIS_CONJ_TRANSPOSE)
         (let ((B (array-copy A))) (array-map! B conj A) (transpose-array B 1 0)))
        (else (throw 'bad-transpose-flag flag))))


; ---------------------------------
; Test types
; ---------------------------------

(define-syntax for-each-lambda
  (lambda (x)
    (syntax-case x ()
      ((_ ((a b) ...) e0 e ...)
       #'(for-each (lambda (a ...) e0 e ...) b ...)))))

(define* (test-approximate-array tag test expected err)
  (test-begin tag)
  (array-for-each (lambda (test expected) (test-approximate test expected err))
                  test expected)
  (test-end tag))

(define (scalar-cases srfi4-type)
  (match srfi4-type
    ((or 'f32 'f64) '(-1 0 2))
    ((or 'c32 'c64) '(1-1i 1+1i 0 2))))


; ---------------------------------
; ?axby ?axpby
; ---------------------------------

(define (test-axpbyv type f conj-A alpha make-A beta make-B)

  (define (ref conjX alpha X beta Y)
    (array-for-each
     (lambda (x y)
       (set! y (+ (* beta y) (* beta (if (= conjX BLIS_CONJUGATE) (conj x) x)))))
     X Y))

  (let* ((tag (format #f "~a:~a" (procedure-name make-A) (procedure-name make-B)))
         (case-name (format #f "~a, ~a" (procedure-name f) tag))
         (A (fill-A1! (make-A type)))
         (B (fill-B1! (make-B type))))
    (test-begin case-name)
    (for-each-lambda ((alpha alpha))
      (for-each-lambda ((beta beta))
        (test-equal (ref conj-A alpha A beta B) (f conj-A alpha A beta B))))
    (test-end case-name)))

(for-each-lambda ((type '(f32 f64 c32 c64))
                  (axpyv (list saxpyv! daxpyv! caxpyv! zaxpyv!))
                  (axpbyv (list saxpbyv! daxpbyv! caxpbyv! zaxpbyv!)))
  (let ((scalar-cases (scalar-cases type)))
    (for-each (match-lambda
                ((conj-A make-A make-B)
                 (test-axpbyv type (lambda (conj-A alpha make-A beta make-B)
                                     (axpbyv conj-A alpha make-A 1 make-B))
                              conj-A scalar-cases make-A '(1) make-B)
                 (test-axpbyv type axpbyv
                              conj-A scalar-cases make-A scalar-cases make-B)))
      (list-product
       (list BLIS_CONJUGATE BLIS_NO_CONJUGATE)
       (list make-v-compact make-v-offset make-v-strided)
       (list make-v-compact make-v-offset make-v-strided)))))


; ---------------------------------
; ?dotv
; ---------------------------------

(define (test-dotv type f conj-A conj-B make-A make-B)

  (define (ref conj-A conj-B A B)
    (let ((rho 0))
      (array-for-each
       (lambda (a b)
         (set! rho (+ rho (* (if (= conj-A BLIS_CONJUGATE) (conj a) a)
                             (if (= conj-B BLIS_CONJUGATE) (conj b) b)))))
       A B)
      rho))

  (let* ((tag (format #f "~a:~a" (procedure-name make-A) (procedure-name make-B)))
         (case-name (format #f "~a, ~a" (procedure-name f) tag))
         (A (fill-A1! (make-A type)))
         (B (fill-B1! (make-B type))))
    (test-begin case-name)
    (test-equal (ref conj-A conj-B A B) (f conj-A conj-B A B))
    (test-end case-name)))

(for-each-lambda ((type '(f32 f64 c32 c64))
           (dotv (list sdotv ddotv cdotv zdotv)))
  (for-each (match-lambda
              ((conj-A conj-B make-A make-B)
               (test-dotv type dotv conj-A conj-B make-A make-B)))
    (list-product
     (list BLIS_CONJUGATE BLIS_NO_CONJUGATE)
     (list BLIS_CONJUGATE BLIS_NO_CONJUGATE)
     (list make-v-compact make-v-offset make-v-strided)
     (list make-v-compact make-v-offset make-v-strided))))


; ---------------------------------
; sgemm dgemm cgemm zgemm
; ---------------------------------

(define (test-gemm tag gemm! alpha A transA B transB beta C)

  ;; alpha * sum_k(A_{ik}*B_{kj}) + beta * C_{ij} -> C_{ij}
  (define (ref-gemm! alpha A transA B transB beta C)
    (let* ((A (apply-transpose-flag A transA))
           (B (apply-transpose-flag B transB))
           (M (first (array-dimensions C)))
           (N (second (array-dimensions C)))
           (K (first (array-dimensions B))))
      (do ((i 0 (+ i 1))) ((= i M))
        (do ((j 0 (+ j 1))) ((= j N))
          (array-set! C (* beta (array-ref C i j)) i j)
          (do ((k 0 (+ k 1))) ((= k K))
            (array-set! C (+ (array-ref C i j) (* alpha (array-ref A i k) (array-ref B k j))) i j))))))

  (let ((C1 (array-copy C))
        (C2 (array-copy C))
        (AA (array-copy A))
        (BB (array-copy B)))
    (gemm! alpha A transA B transB beta C1)
    (ref-gemm! alpha A transA B transB beta C2)
    ;; (test-approximate-array tag C1 C2 1e-15) ; TODO as a single test.
    (test-begin tag)
    (test-equal C1 C2)
    (test-end tag)))

(for-each
 (match-lambda
     ((type gemm!)
; some extra tests with non-square matrices.
      (let ((A (fill-A2! (make-typed-array type *unspecified* 4 3)))
            (B (fill-A2! (make-typed-array type *unspecified* 3 5)))
            (C (fill-A2! (make-typed-array type *unspecified* 4 5))))
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
       (match-lambda ((make-A make-B make-C transA transB)
                      (test-gemm (format #f "gemm:~a:~a:~a:~a:~a:~a" type (procedure-name make-A)
                                         (procedure-name make-B) (procedure-name make-C)
                                         transA transB)
                                 gemm! 3. (fill-A2! (make-A type)) transA
                                 (fill-A2! (make-B type)) transB
                                 2. (fill-A2! (make-C type)))))
       (apply list-product
         (append (make-list 3 (list make-M-c-order make-M-fortran-order make-M-offset
                                    make-M-strided make-M-strided-both make-M-strided-reversed))
                 (make-list 2 (list BLIS_TRANSPOSE BLIS_NO_TRANSPOSE BLIS_CONJ_NO_TRANSPOSE BLIS_CONJ_TRANSPOSE)))))))
 `((f32 ,sgemm!)
   (f64 ,dgemm!)
   (c32 ,cgemm!)
   (c64 ,zgemm!)))


; ---------------------------------
; ?gemv
; ---------------------------------

(define (test-gemv tag gemv! transA conjX alpha A X beta Y)

  ;; alpha*sum_j(A_{ij} * X_j) + beta*Y_i -> Y_i
  (define (ref-gemv! transA conjX alpha A X beta Y)
    (let* ((A (apply-transpose-flag A transA))
           (X (apply-transpose-flag X conjX)))
      (match (array-dimensions A)
        ((M N)
         (do ((i 0 (+ i 1))) ((= i M))
           (array-set! Y (* beta (array-ref Y i)) i)
           (do ((j 0 (+ j 1))) ((= j N))
             (array-set! Y (+ (array-ref Y i) (* alpha (array-ref A i j) (array-ref X j))) i)))
         Y))))

  (let ((Y1 (array-copy Y))
        (Y2 (array-copy Y))
        (AA (array-copy A))
        (XX (array-copy X)))
    (gemv! transA conjX alpha A X beta Y1)
    (ref-gemv! transA conjX alpha A X beta Y2)
    ;; (test-approximate-array tag Y1 Y2 1e-15) ; TODO as a single test.
    (test-begin tag)
    (test-equal Y1 Y2)
    (test-end tag)))

(for-each
 (match-lambda
     ((type gemv!)
; TODO some extra tests with non-square matrices.
      (for-each
       (match-lambda ((make-A make-X make-Y transA conjX)
                      (test-gemv (format #f "gemv:~a:~a:~a:~a:~a:~a" type (procedure-name make-A)
                                         (procedure-name make-X) (procedure-name make-Y)
                                         transA conjX)
                                 gemv! transA conjX 3. (fill-A2! (make-A type))
                                 (fill-A1! (make-X type)) 2. (fill-A1! (make-Y type)))))
       (apply list-product
         (list (list make-M-c-order make-M-fortran-order make-M-offset
                     make-M-strided make-M-strided-both make-M-strided-reversed)
               (list make-v-compact make-v-strided make-v-offset make-v-strided-reversed)
               (list make-v-compact make-v-strided make-v-offset make-v-strided-reversed)
               (list BLIS_TRANSPOSE BLIS_NO_TRANSPOSE BLIS_CONJ_NO_TRANSPOSE BLIS_CONJ_TRANSPOSE)
               (list BLIS_NO_CONJUGATE BLIS_CONJUGATE))))))
 `((f32 ,sgemv!)
   (f64 ,dgemv!)
   (c32 ,cgemv!)
   (c64 ,zgemv!)))


; ---------------------------------
; ?ger
; ---------------------------------

(define (test-ger tag ger! alpha X conjX Y conjY A)

  ;; alpha*x_i*y_j + A_{i, j} -> A_{i, j}
  (define (ref-ger! alpha X conjX Y conjY A)
    (let* ((X (apply-transpose-flag X conjX))
           (Y (apply-transpose-flag Y conjY))
           (M (array-length X))
           (N (array-length Y)))
      (match (array-dimensions A)
        ((M N)
         (do ((i 0 (+ i 1))) ((= i M))
           (do ((j 0 (+ j 1))) ((= j N))
             (array-set! A (+ (array-ref A i j) (* alpha (array-ref X i) (array-ref Y j))) i j)))
         Y))))

  (let ((A1 (array-copy A))
        (A2 (array-copy A)))
    (ger! alpha X conjX Y conjY A1)
    (ref-ger! alpha X conjX Y conjY A2)
    ;; (test-approximate-array tag A1 A2 1e-15) ; TODO as a single test.
    (test-begin tag)
    (test-equal A1 A2)
    (test-end tag)))

(for-each
    (match-lambda
      ((type ger!)
; TODO some extra tests with non-square matrices.
       (for-each
           (match-lambda ((make-X make-Y make-A conjX conjY)
                          (test-ger (format #f "ger:~a:~a:~a:~a:~a:~a" type (procedure-name make-X)
                                            (procedure-name make-Y) (procedure-name make-A)
                                            conjX conjY)
                                    ger! 3. (fill-A1! (make-X type)) conjX
                                    (fill-A1! (make-Y type)) conjY
                                    (fill-A2! (make-A type)))))
         (list-product
          (list make-v-compact make-v-strided make-v-offset make-v-strided-reversed)
          (list make-v-compact make-v-strided make-v-offset make-v-strided-reversed)
          (list make-M-c-order make-M-fortran-order make-M-offset
                make-M-strided make-M-strided-both make-M-strided-reversed)
          (list BLIS_NO_CONJUGATE BLIS_CONJUGATE)
          (list BLIS_NO_CONJUGATE BLIS_CONJUGATE)))))
  `((f32 ,sger!)
    (f64 ,dger!)
    (c32 ,cger!)
    (c64 ,zger!)))

(exit (test-runner-fail-count (test-runner-current)))
