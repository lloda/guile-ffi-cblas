
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

; Basic tests.

(define (test-dot srfi4-type dot tag make-A make-B)
  (let ((case-name (format #f "~a, ~a" (procedure-name dot) tag))
        (A (make-A srfi4-type))
        (B (make-B srfi4-type)))
    (array-index-map! A values)
    (array-index-map! B values)
    (test-begin case-name)
    (test-equal (dot A B) (fold (lambda (a c) (+ (* a a) c)) 0. (iota 10)))
    (test-end case-name)))

(map
 (match-lambda
     ((srfi4-type dot)
      (test-dot srfi4-type dot "compact-compact" make-A-compact make-A-compact)
      (test-dot srfi4-type dot "compact-compact" make-A-compact make-A-compact)
      (test-dot srfi4-type dot "compact-strided" make-A-compact make-A-strided)
      (test-dot srfi4-type dot "compact-offset" make-A-compact make-A-offset)
      (test-dot srfi4-type dot "strided-compact" make-A-strided make-A-compact)
      (test-dot srfi4-type dot "strided-strided" make-A-strided make-A-strided)
      (test-dot srfi4-type dot "strided-offset" make-A-strided make-A-offset)
      (test-dot srfi4-type dot "offset-compact" make-A-offset make-A-compact)
      (test-dot srfi4-type dot "offset-strided" make-A-offset make-A-strided)
      (test-dot srfi4-type dot "offset-offset" make-A-offset make-A-offset)))
 `((f64 ,ddot)
   (f32 ,sdot)))

(define (conj a) (make-rectangular (real-part a) (- (imag-part a))))
(define dotc-step (lambda (a b c) (+ c (* (conj a) b))))
(define dotu-step (lambda (a b c) (+ c (* a b))))

(define (test-dot-complex srfi4-type dot dot-step tag make-A make-B)
  (let ((case-name (format #f "~a, ~a" (procedure-name dot) tag))
        (A (make-A srfi4-type))
        (B (make-B srfi4-type)))
    (array-index-map! A (lambda (i) (make-rectangular (+ i 1) (- i 5))))
    (array-index-map! B (lambda (i) (make-rectangular (+ i 2) (- 5 i))))
    (test-begin case-name)
    (test-equal (dot A B) (fold dot-step 0. (array->list A) (array->list B)))
    (test-end case-name)))

(map
 (match-lambda
     ((srfi4-type dot dot-step)
      (test-dot-complex srfi4-type dot dot-step "compact-compact" make-A-compact make-A-compact)
      (test-dot-complex srfi4-type dot dot-step "compact-compact" make-A-compact make-A-compact)
      (test-dot-complex srfi4-type dot dot-step "compact-strided" make-A-compact make-A-strided)
      (test-dot-complex srfi4-type dot dot-step "compact-offset" make-A-compact make-A-offset)
      (test-dot-complex srfi4-type dot dot-step "strided-compact" make-A-strided make-A-compact)
      (test-dot-complex srfi4-type dot dot-step "strided-strided" make-A-strided make-A-strided)
      (test-dot-complex srfi4-type dot dot-step "strided-offset" make-A-strided make-A-offset)
      (test-dot-complex srfi4-type dot dot-step "offset-compact" make-A-offset make-A-compact)
      (test-dot-complex srfi4-type dot dot-step "offset-strided" make-A-offset make-A-strided)
      (test-dot-complex srfi4-type dot dot-step "offset-offset" make-A-offset make-A-offset)))
 `((c64 ,zdotc ,dotc-step)
   (c64 ,zdotu ,dotu-step)
   (c32 ,cdotu ,dotu-step)
   (c32 ,cdotu ,dotu-step)))

(unless (zero? (test-runner-fail-count (test-runner-current)))
  (error "FAILED test-ffi-cblas.csm"))
