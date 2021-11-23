
; Test helpers.
; (c) Daniel Llorens - 2014, 2020

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU Lesser General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(import (srfi srfi-1) (srfi srfi-11))

; ---------------------------------
; Various sorts of arrays.
; ---------------------------------

(define (make-v-compact type)
  (make-typed-array type *unspecified* 10))

(define (make-v-strided type)
  (make-shared-array (make-typed-array type *unspecified* 20) (lambda (i) (list (* i 2))) 10))

(define (make-v-offset type)
  (make-shared-array (make-typed-array type *unspecified* 12) (lambda (i) (list (+ i 2))) 10))

(define (make-v-strided-reversed type)
  (make-shared-array (make-typed-array type *unspecified* 20) (lambda (i) (list (- 19 (* i 2)))) 10))

(define (make-v-z type)
  (make-shared-array (make-typed-array type *unspecified*) (lambda (i) (list)) 10))

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

(define (make-M-strided-reversed type)
  (make-shared-array (make-typed-array type *unspecified* 20 20) (lambda (i j) (list (- 19 (* i 2)) (- 19 (* j 2)))) 10 10))

(define (make-M-offset type)
  (make-shared-array (make-typed-array type *unspecified* 12 13) (lambda (i j) (list (+ i 2) (+ j 3))) 10 10))

(define (make-M-z0 type)
  (make-shared-array (make-typed-array type *unspecified* 12) (lambda (i j) (list (+ 2 j))) 10 10))

(define (make-M-z1 type)
  (make-shared-array (make-typed-array type *unspecified* 12) (lambda (i j) (list (+ 2 i))) 10 10))

(define (make-M-z00 type)
  (make-shared-array (make-typed-array type *unspecified*) (lambda (i j) (list)) 10 10))

(define (make-M-overlap type)
  (make-shared-array (make-typed-array type *unspecified* 20) (lambda (i j) (list (+ i j))) 10 10))

(define (make-M-overlap-reversed type)
  (make-shared-array (make-typed-array type *unspecified* 20) (lambda (i j) (list (- 18 i j))) 10 10))

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

(define (array-map op A)
  (let ((B (apply make-typed-array (array-type A) *unspecified* (array-dimensions A))))
    (array-map! B op A)
    B))

(define (list-product . rest)
  "make a list of all lists with 1st element from 1st arg, 2nd element from

   (list-product '(1 2)) => ((1) (2))
   (list-product '(1 2) '(3 4)) => ((1 3) (1 4) (2 3) (2 4))
   (list-product '(1 2) '(3 4) '(5 6)) => ((1 3 5) (1 3 6) etc.)"
  (cond
    ((null? rest)
      '())
    ((null? (cdr rest))
      (map list (car rest)))
    (else
      (let ((list-product-rest (apply list-product (cdr rest))))
        (append-map!
          (lambda (p)
            (map
              (lambda (v) (cons p v))
              list-product-rest))
          (car rest))))))
