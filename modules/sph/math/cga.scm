(library (sph math cga)
  (export)
  (import
    (guile)
    (rnrs arithmetic bitwise)
    (rnrs sorting)
    (sph)
    (sph hashtable)
    (sph list)
    (sph set))

  (define sph-math-cga-description
    "blade: (integer:id . number:scale-factor)
     blade-id
       each basis blade is a bit and outer products are xor combinations")

  ;-- generic

  (define* (bit-range count #:optional (start 0))
    "integer [integer] -> (integer ...)
     return count next bits in a list, that is, integers with only one bit set at one position.
     example
       (bit-range 3) -> (000 001 010)"
    (let (bit-last (bitwise-arithmetic-shift-left 1 (+ start (- count 1))))
      (pair start
        (let loop ((a (+ 1 start)))
          (if (< a bit-last) (pair a (loop (bitwise-arithmetic-shift-left a 1))) null)))))

  (define (bit-count a)
    "integer -> integer
     count all 1 bits in \"a\""
    ; repeatedly shift bits right and check the first bit
    (let loop ((bits a) (result 0))
      (if (zero? bits) result
        (loop (bitwise-arithmetic-shift-right bits 1)
          (if (= 1 (bitwise-and 1 bits)) (+ 1 result) result)))))

  (define (factorial a) (if (<= a 1) 1 (* a (factorial (- a 1)))))
  ;
  ;-- basis blade representation
  ;
  (define blade-id-grade bit-count)

  (define (blade-id-sign a b) "integer:blade-id integer:blade-id -> 1/-1"
    ; start after grade-one bit
    (let loop ((sum 0) (n (bitwise-arithmetic-shift-right a 1)))
      (if (zero? n) (if (zero? (bitwise-and 1 sum)) 1 -1)
        (loop (+ sum (blade-id-grade (bitwise-and n b))) (bitwise-arithmetic-shift-right n 1)))))

  (define* (blade-id-geometric-product a b #:optional metric)
    "integer:id integer:id [(integer ...)] -> blade
     create a new basis blade that represents the geometric product of two basis blades"
    (let
      (result
        (blade-new (bitwise-xor a b) (* (blade-id-sign a b) (blade-scalar a) (blade-scalar b))))
      (if metric
        ; rescale according to metric
        (let loop
          ( (meet (bitwise-and (blade-id a) (blade-id b))) (scale (blade-scalar result))
            (metric metric))
          (if (zero? meet) (blade-new (blade-id result) scale)
            (loop (bitwise-arithmetic-shift-right meet 1)
              (if (zero? (bitwise-and 1 meet)) scale (* scale (first metric))) (tail metric))))
        result)))

  (define (blade-id-outer-product a b)
    "integer integer -> blade
     creates a new blade that represents the result of the outer product
     for the given basis blades"
    ; if any overlap then zero
    (if (zero? (bitwise-and a b)) (blade-id-geometric-product a b) (blade-new 0 0)))

  (define (blade-id->string a)
    "integer ->string
     return a string name for a basis blade id"
    (let loop ((a a) (grade 0) (numbers null))
      (if (zero? a) (if (> grade 0) (apply string-append "e" (reverse numbers)) "s")
        (loop (bitwise-arithmetic-shift-right a 1) (+ grade 1)
          (if (= 1 (bitwise-and 1 a)) (pair (number->string (+ grade 1)) numbers) numbers)))))

  (define (blade-ids-sort a) "list -> list"
    ; sort by grade and value
    (list-sort (l (a b) (or (< (blade-id-grade a) (blade-id-grade b)) (< a b))) a))

  (define (blade-ids metric)
    "(integer ...) -> (pair:blade ...)
     create basis blade ids of all grades for the current metric"
    (blade-ids-sort (map-integers (expt 2 (length metric)) identity)))

  (define (blade-equal? a b)
    "note: float scalar comparison can be inexact and lead to unexpected results"
    (and (= (blade-id a) (blade-id b)) (= (blade-scalar a) (blade-scalar b))))

  (define (blade-inner-product-f a grade-a grade-b type)
    "type: left-contraction/right-contraction/hestenes"
    (let (a-grade (blade-id-grade (blade-id a)))
      (case type
        (left-contraction
          (if (or (> grade-a grade-b) (not (= a-grade (- grade-b grade-a)))) (blade-new 0 0) a))
        (right-contraction
          (if (or (< grade-a grade-b) (not (= a-grade (- grade-a grade-b)))) (blade-new 0 0) a))
        (hestenes
          (if (or (zero? grade-a) (zero? grade-b)) (blade-new 0 0)
            (if (= a-grade (abs (- grade-a grade-b))) a (blade-new 0 0)))))))

  (define (blade-grade-inversion a)
    "blade -> blade
     grade inversion"
    (blade-new (blade-id a) (expt -1 (blade-id-grade (blade-id a)))))

  (define (blade-reverse a) "blade -> blade"
    (let (g (blade-id-grade (blade-id a))) (blade-new (blade-id a) (expt -1 (/ (* g (- g 1)) 2)))))

  (define (blade-conjugate a)
    "blade -> blade
     clifford conjugate"
    (let (g (blade-id-grade (blade-id a))) (blade-new (blade-id a) (expt -1 (/ (* g (+ g 1)) 2)))))

  (define (blade-new id s) "integer:bits integer:scalar:scale-factor -> pair:blade" (pair id s))
  (define blade-id first)
  (define blade-scalar second)

  (define* (blade-count dimensions #:optional grade)
    "integer [integer] -> integer
     return the number of basis blades of a given grade.
     when grade is not given then return the total number of basis blades"
    (if grade (/ (factorial dimensions) (factorial (- dimensions grade)) (factorial grade))
      (expt 2 dimensions)))

  ;-- other
  ;
  (define dim 4)
  (define b (blade-ids (make-list dim 1)))
  (if (not (= (blade-count dim) (length b))) (display-line "error: unexpected blade count"))

  (each (l (a) (display-line (list (number->string a 2) (blade-id-grade a) (blade-id->string a))))
    b)

  (define (space-new metric) "(integer:1/-1 ...) ->"))
