(define-test-module (test module sph math cga)
  (import
    (sph math cga))

  (define-test (basis-ids a)
    (let* ((dim (first a)) (ids (basis-ids dim))) (= (basis-count dim) (length ids))))

  (define-test (basis-id->string) "also tests the sorting of basis-ids"
    (let*
      ( (expected
          (list "s" "e1"
            "e2" "e3"
            "e4" "e1e2"
            "e1e3" "e2e3" "e1e4" "e2e4" "e3e4" "e1e2e3" "e1e2e4" "e1e3e4" "e2e3e4" "e1e2e3e4"))
        (strings (map basis-id->string (basis-ids 4))))
      (equal? expected strings)))

  ; mv-norm, mv-op
  ;(define p1 (vec-new 0.96 0.0 0.02))
  ;(define p2 (vec-new -0.19 0.86 0.16))
  ;(debug-log (mv-norm (mv-op p1 p2)))
  ;op: (0.8256 0.15739999999999998 -0.0172)
  ;op norm: 0.8406461562393538

  (test-execute-procedures-lambda basis-id->string (basis-count 1 2 2 4 3 8 4 16 5 32 6 64)
    (basis-ids 1 #t 2 #t 3 #t 4 #t 5 #t 6 #t) (mv-null? ((unquote (mv-new))) #t)
    (basis-ip (unquote (list (basis-new 1 1) (basis-new 1 1))) (0 . 1))
    (basis-op (unquote (list (basis-new 1 1) (basis-new 1 1))) (0 . 0))
    (basis-gp (unquote (list (basis-new 1 3) (basis-new 2 4))) (3 . 12)
      (unquote (list (basis-new 1 1) (basis-new 1 1))) (0 . 1)
      (unquote (list (basis-new 1 1) (basis-new 2 1))) (3 . 1))
    (mv-simplify
      ( (unquote
          (mv-new (basis-new 1 3) (basis-new 2 4) (basis-new 2 3) (basis-new 2 3) (basis-new 1 4))))
      ((1 . 7) (2 . 10)))
    #;(mv-gp
      (unquote
        (list (mv-new (basis-new 1 3) (basis-new 2 4)) (mv-new (basis-new 1 5) (basis-new 2 9))))
      #t)))
