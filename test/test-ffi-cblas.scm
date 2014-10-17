
; Tests for (ffi cblas).
; (c) Daniel Llorens - 2014

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(import (ffi cblas) (srfi srfi-64) (srfi srfi-1) (ice-9 match))

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
; Test types
; ---------------------------------

(define* (test-approximate-array test expected err)
  (test-begin "approximate array")
  (array-for-each (lambda (test expected) (test-approximate test expected err))
                  test expected)
  (test-end "approximate array"))

; ---------------------------------
; sdot ddot cdotc cdotu zdotc zdotu
; ---------------------------------

(define (conj a) (make-rectangular (real-part a) (- (imag-part a))))
(define dotc-step (lambda (a b c) (+ c (* (conj a) b))))
(define dotu-step (lambda (a b c) (+ c (* a b))))

(define (test-dot srfi4-type dot dot-step make-A make-B)
  (let* ((tag (format #f "~a:~a" (procedure-name make-A) (procedure-name make-B)))
         (case-name (format #f "~a, ~a" (procedure-name dot) tag))
         (A (fill-A1! (make-A srfi4-type)))
         (B (fill-B1! (make-B srfi4-type))))
    (test-begin case-name)
    (test-equal (dot A B) (fold dot-step 0. (array->list A) (array->list B)))
    (test-end case-name)))

(map
 (match-lambda
     ((srfi4-type dot dot-step)
      (for-each
       (lambda (make-X)
         (for-each
          (lambda (make-Y)
            (test-dot srfi4-type dot dot-step make-X make-Y))
          (list make-v-compact make-v-offset make-v-strided)))
       (list make-v-compact make-v-offset make-v-strided))))
 `((f64 ,ddot ,dotu-step)
   (f32 ,sdot ,dotu-step)
   (c64 ,zdotc ,dotc-step)
   (c64 ,zdotu ,dotu-step)
   (c32 ,cdotu ,dotu-step)
   (c32 ,cdotu ,dotu-step)))

; ---------------------------------
; saxpy daxpy caxpy zaxpy
; ---------------------------------

(define (test-axpy srfi4-type axpy! make-A make-B)
  (let* ((tag (format #f "~a:~a" (procedure-name make-A) (procedure-name make-B)))
         (case-name (format #f "~a, ~a" (procedure-name axpy!) tag))
         (A (fill-A1! (make-A srfi4-type)))
         (B (fill-B1! (make-B srfi4-type))))
    (let ((Alist (array->list A))
          (Blist (array->list B)))
      (test-begin case-name)
      (axpy! 3 A B)
      (test-equal B (list->typed-array srfi4-type 1 (map (lambda (a b) (+ (* 3 a) b)) Alist Blist)))
      (test-equal A (list->typed-array srfi4-type 1 Alist))
      (axpy! 1.9 A B)
      (test-approximate-array
       B (list->typed-array srfi4-type 1 (map (lambda (a b) (+ (* a (+ 3 1.9)) b)) Alist Blist)) 1e-14)
      (test-equal A (list->typed-array srfi4-type 1 Alist))
      (test-end case-name))))

(map
 (match-lambda
     ((srfi4-type axpy!)
      (for-each
       (lambda (make-X)
         (for-each
          (lambda (make-Y)
            (test-axpy srfi4-type axpy! make-X make-Y))
          (list make-v-compact make-v-offset make-v-strided)))
       (list make-v-compact make-v-offset make-v-strided))))
 `((f64 ,daxpy!)
   (f32 ,saxpy!)
   (c64 ,zaxpy!)
   (c32 ,caxpy!)))

; ---------------------------------
; sgemv dgemv cgemv zgemv
; ---------------------------------

(define (array-copy A)
  (let ((B (apply make-typed-array (array-type A) *unspecified* (array-dimensions A))))
    (array-copy! A B)
    B))

(define (ref-gemv! alpha A X beta Y)
  (match (array-dimensions A)
    ((M N)
     (do ((i 0 (+ i 1))) ((= i M))
       (array-set! Y (* beta (array-ref Y i)) i)
       (do ((j 0 (+ j 1))) ((= j N))
         (array-set! Y (+ (array-ref Y i) (* alpha (array-ref A i j) (array-ref X j))) i)))
     Y)))

(define (test-gemv srfi4-type gemv! make-A make-X make-Y)
  (let* ((tag (format #f "~a:~a:~a" (procedure-name make-A) (procedure-name make-X) (procedure-name make-Y)))
         (case-name (format #f "~a, ~a" (procedure-name gemv!) tag))
         (A (fill-A2! (make-A srfi4-type)))
         (X (fill-A1! (make-X srfi4-type)))
         (Y (fill-B1! (make-Y srfi4-type))))
    (let ((A1 (array-copy A))
          (X1 (array-copy X)))
      (test-begin case-name)
      (let ((Y1 (array-copy Y))
            (Y2 (array-copy Y)))
        (gemv! 2. A X 3. Y1)
        (ref-gemv! 2. A X 3. Y2)
        (test-equal Y1 Y2)
        (test-equal A A1)
        (test-equal X X1)
        (test-end case-name)))))

(map
 (match-lambda
     ((srfi4-type gemv!)
      (for-each
       (lambda (make-A)
         (for-each
          (lambda (make-X)
            (for-each
             (lambda (make-Y)
               (test-gemv srfi4-type gemv! make-A make-X make-Y))
             (list make-v-compact make-v-offset make-v-strided)))
          (list make-v-compact make-v-offset make-v-strided)))
       (list make-M-c-order make-M-fortran-order make-M-offset make-M-strided))))
 `((f64 ,dgemv!)
   (f32 ,sgemv!)
   (c64 ,zgemv!)
   (c32 ,cgemv!)))

(unless (zero? (test-runner-fail-count (test-runner-current)))
  (error "FAILED test-ffi-cblas.csm"))
