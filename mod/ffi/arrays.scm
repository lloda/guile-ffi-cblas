
; Common functions for CBLAS & BLIS bindings.
; (c) Daniel Llorens - 2014-2015, 2017, 2019

; This library is free software; you can redistribute it and/or modify it under
; the terms of the GNU General Public License as published by the Free
; Software Foundation; either version 3 of the License, or (at your option) any
; later version.

(define-module (ffi arrays)
  #:export (srfi4-type-size
            check-array check-2-arrays
            stride dim
            syntax->list define-sdcz))

(import (system foreign) (srfi srfi-1) (srfi srfi-11) (srfi srfi-26))

(define (stride A i)
  (list-ref (shared-array-increments A) i))

(define (dim A i)
  (list-ref (array-dimensions A) i))

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

; https://www.scheme.com/csug8/syntax.html §11.3
(define syntax->list
  (lambda (ls)
    (syntax-case ls ()
      (() '())
      ((x . r) (cons #'x (syntax->list #'r))))))

(eval-when (expand load eval)
  (define (subst-qmark stx-name t)
    (let* ((s (symbol->string (syntax->datum stx-name)))
           (i (string-index s #\?))
           (fmt (string-replace s "~a" i (+ i 1))))
      (datum->syntax stx-name (string->symbol (format #f fmt t))))))

(define-syntax define-sdcz
  (lambda (x)
    (syntax-case x ()
      ((_ definer n ...)
       (cons #'begin
             (append-map
              (lambda (tag t)
                (let ((fun (map (cut subst-qmark <> t) (syntax->list #'(n ...)))))
; #`(quote #,(datum->syntax x tag)) to write out a symbol, but assembling docstrings seems harder (?)
                  (list (cons* #'definer (datum->syntax x tag) fun)
                        (cons* #'export fun))))
              '(f32 f64 c32 c64)
              '(s d c z)))))))
