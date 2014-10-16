
; Tests for (ffi cblas).
; (c) Daniel Llorens - 2014

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(import (ffi cblas) (srfi srfi-64) (srfi srfi-1))

; Various sorts of arrays.

(define (make-A-compact type)
  (let ((A (make-typed-array type *unspecified* 10)))
    (array-index-map! A values)
    A))

(define (make-A-strided type)
  (let ((A (make-shared-array (make-typed-array type *unspecified* 20) (lambda (i) (list (* i 2))) 10)))
    (array-index-map! A values)
    A))

(define (make-A-offset type)
  (let ((A (make-shared-array (make-typed-array type *unspecified* 12) (lambda (i) (list (+ i 2))) 10)))
    (array-index-map! A values)
    A))

; Basic tests.

(define (test-ddot tag make-A make-B)
  (let ((case-name (format #f "ddot, ~a" tag))
        (A (make-A 'f64))
        (B (make-B 'f64)))
    (test-begin case-name)
    (test-equal (ddot A B) (fold (lambda (a c) (+ (* a a) c)) 0. (iota 10)))
    (test-end case-name)))

(test-ddot "compact-compact" make-A-compact make-A-compact)
(test-ddot "compact-compact" make-A-compact make-A-compact)
(test-ddot "compact-strided" make-A-compact make-A-strided)
(test-ddot "compact-offset" make-A-compact make-A-offset)
(test-ddot "strided-compact" make-A-strided make-A-compact)
(test-ddot "strided-strided" make-A-strided make-A-strided)
(test-ddot "strided-offset" make-A-strided make-A-offset)
(test-ddot "offset-compact" make-A-offset make-A-compact)
(test-ddot "offset-strided" make-A-offset make-A-strided)
(test-ddot "offset-offset" make-A-offset make-A-offset)

(unless (zero? (test-runner-fail-count (test-runner-current)))
  (error "FAILED test-ffi-cblas.csm"))
