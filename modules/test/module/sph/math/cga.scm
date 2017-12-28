(define-test-module (test module sph math cga)
  (import
    (sph math cga))

  (define-test (blade-ids a)
    (let* ((dim (first a)) (ids (blade-ids (make-list dim 1)))) (= (blade-count dim) (length ids))))

  (define-test (blade-id->string) "also tests the sorting of blade-ids"
    (let*
      ( (expected
          (list "s" "e1"
            "e2" "e3" "e4" "e12" "e13" "e23" "e14" "e24" "e34" "e123" "e124" "e134" "e234" "e1234"))
        (strings (map blade-id->string (blade-ids (make-list 4 1)))))
      (equal? expected strings)))

  (test-execute-procedures-lambda blade-id->string (blade-count 1 2 2 4 3 8 4 16 5 32 6 64)
    (blade-ids 1 #t 2 #t 3 #t 4 #t 5 #t 6 #t) (mv-null? ((unquote (mv-new))) #t)
    (blade-geometric-product (unquote (list (blade-new 1 3) (blade-new 2 4))) (3 . 12))
    (mv-simplify
      ( (unquote
          (mv-new (blade-new 1 3) (blade-new 2 4) (blade-new 2 3) (blade-new 2 3) (blade-new 1 4))))
      ((1 . 7) (2 . 10)))
    (mv-geometric-product
      (unquote
        (list (mv-new (blade-new 1 3) (blade-new 2 4)) (mv-new (blade-new 1 5) (blade-new 2 9))))
      #t)))
