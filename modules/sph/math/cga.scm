(library (sph math cga)
  (export
    blade-conjugate
    blade-count
    blade-equal?
    blade-gp
    blade-grade
    blade-grade-inversion
    blade-id
    blade-id->string
    blade-id-from-string
    blade-id-grade
    blade-id-less?
    blade-id-sign
    blade-ids
    blade-ip
    blade-new
    blade-op
    blade-reverse
    blade-scalar
    mv-dual
    mv-gp
    mv-grade-inversion
    mv-ip
    mv-new
    mv-null?
    mv-op
    mv-reverse
    mv-scalar
    mv-scalar-product
    mv-simplify
    mv-subtract
    mv-sum
    space-new
    sph-math-cga-description)
  (import
    (rnrs arithmetic bitwise)
    (sph)
    (only (guile) string-split string-join)
    (only (rnrs sorting) list-sort)
    (only (sph alist) alist-q)
    (only (sph list one) group)
    (only (sph list)
      produce
      map-map
      map-integers
      list-sort-with-accessor)
    (only (sph string) string-equal?))

  (define sph-math-cga-description
    "mostly non-optimising implementation of geometric algebra elements.
     blade: (integer:id . number:scale-factor)
     multivector: (blade ...)
     blade-id
       each basis blade is identified by a bit
       outer products between basis blades are bitwise xor combinations of the ids")

  ; outline
  ;   generic
  ;   blade-id
  ;   blade
  ;   multivector
  ;   other

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
  (define blade-id-grade bit-count)

  (define (blade-id-sign a b)
    "integer:blade-id integer:blade-id -> 1/-1
     sign change for geometric product of arguments"
    ; start after grade-one bit
    (let loop ((sum 0) (n (bitwise-arithmetic-shift-right a 1)))
      (if (zero? n) (if (zero? (bitwise-and 1 sum)) 1 -1)
        (loop (+ sum (blade-id-grade (bitwise-and n b))) (bitwise-arithmetic-shift-right n 1)))))

  (define (blade-id->string a)
    "integer -> string
     return a string name for a basis blade id"
    (let loop ((a a) (grade 0) (numbers null))
      (if (zero? a) (if (> grade 0) (string-join (reverse numbers) "e" (q prefix)) "s")
        (loop (bitwise-arithmetic-shift-right a 1) (+ grade 1)
          (if (= 1 (bitwise-and 1 a)) (pair (number->string (+ grade 1)) numbers) numbers)))))

  (define (blade-id-from-string a) "string -> integer"
    (if (string-equal? "s" a) 0
      (fold (l (a result) (bitwise-xor result (expt 2 (- a 1)))) 0
        (map string->number (tail (string-split a #\e))))))

  (define (blade-id-less? a b)
    "integer integer -> boolean
     compare blade ids for sorting. considers grade primarily and id value secondarily"
    (or (< (blade-id-grade a) (blade-id-grade b)) (< a b)))

  (define (blade-ids dimensions)
    "(integer ...) -> (pair:blade ...)
     create basis blade ids for all grades for the given metric signature.
     the result list will have (expt 2 dimensions) elements"
    (list-sort blade-id-less? (map-integers (expt 2 dimensions) identity)))

  (define (blade-equal? a b)
    "note: float scalar comparison can be inexact and lead to unexpected results"
    (and (= (blade-id a) (blade-id b)) (= (blade-scalar a) (blade-scalar b))))

  (define* (blade-gp a b #:optional metric)
    "blade blade [(integer ...)] -> blade
     create a new basis blade that represents the geometric product of two basis blades"
    (let*
      ( (a-id (blade-id a)) (b-id (blade-id b)) (id (bitwise-xor a-id b-id))
        (scale (* (blade-id-sign a-id b-id) (blade-scalar a) (blade-scalar b))))
      (if metric
        ; rescale according to metric
        (let loop ((meet (bitwise-and a-id b-id)) (scale scale) (metric metric))
          (if (zero? meet) (blade-new id scale)
            (loop (bitwise-arithmetic-shift-right meet 1)
              (if (zero? (bitwise-and 1 meet)) scale (* scale (first metric))) (tail metric))))
        (blade-new id scale))))

  (define (blade-op a b)
    "blade -> blade
     creates a new blade that represents the result of the outer product
     for the given basis blades"
    ; if any overlap then zero
    (if (zero? (bitwise-and (blade-id a) (blade-id b))) (blade-gp a b) (blade-new 0 0)))

  (define (blade-ip-f a grade-a grade-b type)
    "blade integer integer symbol -> blade
     convert a geometric product to an inner product.
     type: left/right/h, where left/right is for contraction"
    (let (a-grade (blade-id-grade (blade-id a)))
      (case type
        ( (left)
          (if (or (> grade-a grade-b) (not (= a-grade (- grade-b grade-a)))) (blade-new 0 0) a))
        ( (right)
          (if (or (< grade-a grade-b) (not (= a-grade (- grade-a grade-b)))) (blade-new 0 0) a))
        ( (h hestenes)
          (if (or (zero? grade-a) (zero? grade-b)) (blade-new 0 0)
            (if (= a-grade (abs (- grade-a grade-b))) a (blade-new 0 0)))))))

  (define* (blade-ip a b #:optional metric type)
    "blade blade [list symbol] -> blade
     create a basis blade that represents the inner product of the arguments"
    (blade-ip-f (blade-gp a b metric) (blade-grade a) (blade-grade b) (or type (q left))))

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
  (define blade-scalar tail)

  (define* (blade-count dimensions #:optional grade)
    "integer [integer] -> integer
     return the number of basis blades of a given grade.
     when grade is not given then return the total number of basis blades"
    (if grade (/ (factorial dimensions) (factorial (- dimensions grade)) (factorial grade))
      (expt 2 dimensions)))

  (define (blade-grade a) (blade-id-grade (blade-id a)))
  (define (mv-new . a) "blade ... -> list:multivector" a)

  (define (mv-simplify a)
    "list:multivector -> list:multivector
     remove blades with zero scale, unite blades with the same grade"
    (let loop ((b (list-sort-with-accessor < blade-id a)))
      (if (null? b) b
        (let (a (first b))
          (if (zero? (blade-scalar a))
            ; skip zero entries
            (loop (tail b))
            ; unite entries with equal grade.
            ; process entries to the right and merge the first with the current one
            (let (b (loop (tail b)))
              (if (or (null? b) (not (= (blade-id a) (blade-id (first b))))) (pair a b)
                (pair (blade-new (blade-id a) (+ (blade-scalar a) (blade-scalar (first b))))
                  (tail b)))))))))

  (define (mv-op a b) (mv-simplify (produce blade-op a b)))

  (define* (mv-ip a b #:optional metric type)
    (mv-simplify (produce (l (a b) (blade-ip a b metric type)) a b)))

  (define mv-gp
    (let (scalar-gp (l (a b) (map (l (a) (blade-new (blade-id a) (* b (blade-scalar a)))) a)))
      (l* (a b #:optional metric) "mv/scala mv/scalar list -> mv"
        (mv-simplify
          (if (pair? a) (if (pair? b) (produce blade-gp a b) (scalar-gp a b)) (scalar-gp b a))))))

  (define (mv-sum . a) "mv/scalar ... -> mv"
    (mv-simplify (apply append (map (l (a) (if (pair? a) a (blade-new 0 a))) a))))

  (define (mv-subtract a . b) "mv/scalar ... -> mv"
    (mv-simplify
      (fold
        (l (b result) (if (pair? b) (append (mv-gp b -1) result) (pair (blade-new 0 (- a)) result)))
        null b)))

  (define (mv-null? a)
    "mv -> boolean
     null if empty or all blades have a scale of 0"
    (null? (mv-simplify a)))

  (define (mv-reverse a) (map blade-reverse a))
  (define (mv-grade-inversion a) (map blade-grade-inversion a))

  (define (mv-dual a metric) "mv list ->"
    (mv-ip (mv-new (blade-new (- (bitwise-arithmetic-shift-left 1 (length metric)) 1) 1)) metric))

  (define (mv-scalar a)
    "mv -> scalar
     returns the sum of all scalar values in multivector (not the pseudoscalar)"
    (apply + (filter (l (a) (zero? (blade-id a))) a)))

  (define* (mv-scalar-product a b #:optional metric) "list:mv list:mv [list:metric]-> mv"
    (mv-scalar (mv-ip a b metric)))

  (define (mv-type base-names)
    "string ... -> procedure:{scalar ... -> mv}
     returns a procedure that when called with scalars returns a multivector with the
     bases specified by base-names"
    (let (bases (map blade-id-from-string base-names)) (l a (map blade-new bases a))))

  (define (mv-types-new config) "((string ...) ...) -> (procedure ...)" (map mv-type config))

  (define* (space-new metric #:key types)
    "(integer ...) [#:types ((basis-blade-name ...) ...)]
     example
       (define s (space-new (list 1 1 1 1)))
       (define vec-new (list-ref (alist-ref-q s mv-types) 0))"
    (let*
      ( (basis-blades (blade-ids (length metric)))
        (blades-per-grade
          (map-map blade-id->string (tail (map tail (group basis-blades blade-id-grade))))))
      (alist-q basis-blades basis-blades
        mv-types (mv-types-new (append blades-per-grade (or types null)))))))
