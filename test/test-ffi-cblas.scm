
; Tests for (ffi cblas).
; (c) Daniel Llorens - 2014

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(import (ffi cblas) (srfi srfi-64) (srfi srfi-1) (ice-9 match))

; Various sorts of arrays.

(define (make-A-compact type)
  (make-typed-array type *unspecified* 10))

(define (make-A-strided type)
  (make-shared-array (make-typed-array type *unspecified* 20) (lambda (i) (list (* i 2))) 10))

(define (make-A-offset type)
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

; Basic tests.

(define (conj a) (make-rectangular (real-part a) (- (imag-part a))))
(define dotc-step (lambda (a b c) (+ c (* (conj a) b))))
(define dotu-step (lambda (a b c) (+ c (* a b))))

(define (test-dot srfi4-type dot dot-step tag make-A make-B)
  (let ((case-name (format #f "~a, ~a" (procedure-name dot) tag))
        (A (fill-A1! (make-A srfi4-type)))
        (B (fill-B1! (make-B srfi4-type))))
    (test-begin case-name)
    (test-equal (dot A B) (fold dot-step 0. (array->list A) (array->list B)))
    (test-end case-name)))

(map
 (match-lambda
     ((srfi4-type dot dot-step)
      (test-dot srfi4-type dot dot-step "compact-compact" make-A-compact make-A-compact)
      (test-dot srfi4-type dot dot-step "compact-compact" make-A-compact make-A-compact)
      (test-dot srfi4-type dot dot-step "compact-strided" make-A-compact make-A-strided)
      (test-dot srfi4-type dot dot-step "compact-offset" make-A-compact make-A-offset)
      (test-dot srfi4-type dot dot-step "strided-compact" make-A-strided make-A-compact)
      (test-dot srfi4-type dot dot-step "strided-strided" make-A-strided make-A-strided)
      (test-dot srfi4-type dot dot-step "strided-offset" make-A-strided make-A-offset)
      (test-dot srfi4-type dot dot-step "offset-compact" make-A-offset make-A-compact)
      (test-dot srfi4-type dot dot-step "offset-strided" make-A-offset make-A-strided)
      (test-dot srfi4-type dot dot-step "offset-offset" make-A-offset make-A-offset)))
 `((f64 ,ddot ,dotu-step)
   (f32 ,sdot ,dotu-step)
   (c64 ,zdotc ,dotc-step)
   (c64 ,zdotu ,dotu-step)
   (c32 ,cdotu ,dotu-step)
   (c32 ,cdotu ,dotu-step)))

(define* (test-approximate-array test expected err)
  (test-begin "approximate array")
  (array-for-each (lambda (test expected) (test-approximate test expected err))
                  test expected)
  (test-end "approximate array"))

(define (test-axpy srfi4-type axpy! tag make-A make-B)
  (let ((case-name (format #f "~a, ~a" (procedure-name axpy!) tag))
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
      (test-axpy srfi4-type axpy! "compact-compact" make-A-compact make-A-compact)
      (test-axpy srfi4-type axpy! "compact-compact" make-A-compact make-A-compact)
      (test-axpy srfi4-type axpy! "compact-strided" make-A-compact make-A-strided)
      (test-axpy srfi4-type axpy! "compact-offset" make-A-compact make-A-offset)
      (test-axpy srfi4-type axpy! "strided-compact" make-A-strided make-A-compact)
      (test-axpy srfi4-type axpy! "strided-strided" make-A-strided make-A-strided)
      (test-axpy srfi4-type axpy! "strided-offset" make-A-strided make-A-offset)
      (test-axpy srfi4-type axpy! "offset-compact" make-A-offset make-A-compact)
      (test-axpy srfi4-type axpy! "offset-strided" make-A-offset make-A-strided)
      (test-axpy srfi4-type axpy! "offset-offset" make-A-offset make-A-offset)))
 `((f64 ,daxpy!)
   (f32 ,saxpy!)
   (c64 ,zaxpy!)
   (c32 ,caxpy!)))

(unless (zero? (test-runner-fail-count (test-runner-current)))
  (error "FAILED test-ffi-cblas.csm"))
