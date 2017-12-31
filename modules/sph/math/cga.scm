(library (sph math cga)
  (export
    basis-conjugate
    basis-count
    basis-equal?
    basis-gp
    basis-grade
    basis-grade-inversion
    basis-id
    basis-id->string
    basis-id->symbol
    basis-id-from-string
    basis-id-from-symbol
    basis-id-grade
    basis-id-less?
    basis-id-sign
    basis-ids
    basis-ip
    basis-new
    basis-op
    basis-reverse
    basis-scalar
    mv-coordinates
    mv-dual
    mv-gp
    mv-grade-inversion
    mv-ip
    mv-new
    mv-norm
    mv-null?
    mv-op
    mv-reverse
    mv-scalar
    mv-simplify
    mv-sp
    mv-subtract
    mv-sum
    space-basis
    space-new
    space-type
    sph-math-cga-description)
  (import
    (rnrs arithmetic bitwise)
    (sph)
    (only (guile)
      string-split
      string-join
      compose)
    (only (rnrs sorting) list-sort)
    (only (sph alist)
      alist-q
      alist-ref
      alist-ref-q)
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
     basis-id
       each basis blade is identified by a bit.
       outer products between basis blades are bitwise xor combinations of the ids")

  ; outline
  ;   generic
  ;   basis-id
  ;   basis
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
  (define basis-id-grade bit-count)

  (define (basis-id-sign a b)
    "integer:basis-id integer:basis-id -> 1/-1
     sign change for geometric product of arguments"
    ; start after grade-one bit
    (let loop ((sum 0) (n (bitwise-arithmetic-shift-right a 1)))
      (if (zero? n) (if (zero? (bitwise-and 1 sum)) 1 -1)
        (loop (+ sum (basis-id-grade (bitwise-and n b))) (bitwise-arithmetic-shift-right n 1)))))

  (define (basis-id->string a)
    "integer -> string
     return a string name for a basis blade id"
    (let loop ((a a) (grade 0) (numbers null))
      (if (zero? a) (if (> grade 0) (string-join (reverse numbers) "e" (q prefix)) "s")
        (loop (bitwise-arithmetic-shift-right a 1) (+ grade 1)
          (if (= 1 (bitwise-and 1 a)) (pair (number->string (+ grade 1)) numbers) numbers)))))

  (define (basis-id->symbol a)
    "integer -> symbol
     like basis-id->string but returns a symbol"
    (string->symbol (basis-id->string a)))

  (define (basis-id-from-symbol a)
    "symbol -> integer
     like basis-id-from-string but takes a symbol"
    (basis-id-from-string (symbol->string a)))

  (define (basis-id-from-string a) "string -> integer"
    (if (string-equal? "s" a) 0
      (fold (l (a result) (bitwise-xor result (expt 2 (- a 1)))) 0
        (map string->number (tail (string-split a #\e))))))

  (define (basis-id-less? a b)
    "integer integer -> boolean
     compare blade ids for sorting. considers grade primarily and id value secondarily"
    (or (< (basis-id-grade a) (basis-id-grade b)) (< a b)))

  (define (basis-ids dimensions)
    "(integer ...) -> (pair:blade ...)
     create basis blade ids for all grades for the given metric signature.
     the result list will have (expt 2 dimensions) elements"
    (list-sort basis-id-less? (map-integers (expt 2 dimensions) identity)))

  (define (basis-equal? a b)
    "note: float scalar comparison can be inexact and lead to unexpected results"
    (and (= (basis-id a) (basis-id b)) (= (basis-scalar a) (basis-scalar b))))

  (define* (basis-gp a b #:optional metric)
    "blade blade [(integer ...)] -> blade
     create a new basis blade that represents the geometric product of two basis blades"
    (let*
      ( (a-id (basis-id a)) (b-id (basis-id b)) (id (bitwise-xor a-id b-id))
        (scale (* (basis-id-sign a-id b-id) (basis-scalar a) (basis-scalar b))))
      (if metric
        ; rescale according to metric
        (let loop ((meet (bitwise-and a-id b-id)) (scale scale) (metric metric))
          (if (zero? meet) (basis-new id scale)
            (loop (bitwise-arithmetic-shift-right meet 1)
              (if (zero? (bitwise-and 1 meet)) scale (* scale (first metric))) (tail metric))))
        (basis-new id scale))))

  (define (basis-op a b)
    "blade -> blade
     creates a new blade that represents the result of the outer product
     for the given basis blades"
    ; if any overlap then zero
    (if (zero? (bitwise-and (basis-id a) (basis-id b))) (basis-gp a b) (basis-new 0 0)))

  (define (basis-ip-f a grade-a grade-b type)
    "blade integer integer symbol -> blade
     convert a geometric product to an inner product.
     type: left/right/h, where left/right is for contraction"
    (let (a-grade (basis-id-grade (basis-id a)))
      (case type
        ( (left)
          (if (or (> grade-a grade-b) (not (= a-grade (- grade-b grade-a)))) (basis-new 0 0) a))
        ( (right)
          (if (or (< grade-a grade-b) (not (= a-grade (- grade-a grade-b)))) (basis-new 0 0) a))
        ( (h hestenes)
          (if (or (zero? grade-a) (zero? grade-b)) (basis-new 0 0)
            (if (= a-grade (abs (- grade-a grade-b))) a (basis-new 0 0)))))))

  (define* (basis-ip a b #:optional metric type)
    "blade blade [list symbol] -> blade
     create a basis blade that represents the inner product of the arguments"
    (basis-ip-f (basis-gp a b metric) (basis-grade a) (basis-grade b) (or type (q left))))

  (define (basis-grade-inversion a)
    "blade -> blade
     grade inversion"
    (basis-new (basis-id a) (expt -1 (basis-id-grade (basis-id a)))))

  (define (basis-reverse a) "blade -> blade"
    (let (g (basis-id-grade (basis-id a))) (basis-new (basis-id a) (expt -1 (/ (* g (- g 1)) 2)))))

  (define (basis-conjugate a)
    "blade -> blade
     clifford conjugate"
    (let (g (basis-id-grade (basis-id a))) (basis-new (basis-id a) (expt -1 (/ (* g (+ g 1)) 2)))))

  (define (basis-new id s) "integer:bits integer:scalar:scale-factor -> pair:blade" (pair id s))
  (define basis-id first)
  (define basis-scalar tail)

  (define* (basis-count dimensions #:optional grade)
    "integer [integer] -> integer
     return the number of basis blades of a given grade.
     when grade is not given then return the total number of basis blades"
    (if grade (/ (factorial dimensions) (factorial (- dimensions grade)) (factorial grade))
      (expt 2 dimensions)))

  (define (basis-grade a) (basis-id-grade (basis-id a)))
  (define (mv-new . a) "blade ... -> list:multivector" a)

  (define (mv-simplify a)
    "list:multivector -> list:multivector
     remove blades with zero scale, unite blades with the same grade"
    (let loop ((b (list-sort-with-accessor < basis-id a)))
      (if (null? b) b
        (let (a (first b))
          (if (zero? (basis-scalar a))
            ; skip zero entries
            (loop (tail b))
            ; unite entries with equal grade.
            ; process entries to the right and merge the first with the current one
            (let (b (loop (tail b)))
              (if (or (null? b) (not (= (basis-id a) (basis-id (first b))))) (pair a b)
                (pair (basis-new (basis-id a) (+ (basis-scalar a) (basis-scalar (first b))))
                  (tail b)))))))))

  (define (mv-op a b)
    ; "produce" creates a cartesian product
    (mv-simplify (produce basis-op a b)))

  (define* (mv-ip a b #:optional metric type)
    (mv-simplify (produce (l (a b) (basis-ip a b metric type)) a b)))

  (define mv-gp
    (let (scalar-gp (l (a b) (map (l (a) (basis-new (basis-id a) (* b (basis-scalar a)))) a)))
      (l* (a b #:optional metric) "mv/scala mv/scalar list -> mv"
        (mv-simplify
          (if (pair? a) (if (pair? b) (produce basis-gp a b) (scalar-gp a b)) (scalar-gp b a))))))

  (define (mv-sum . a) "mv/scalar ... -> mv"
    (mv-simplify (apply append (map (l (a) (if (pair? a) a (basis-new 0 a))) a))))

  (define (mv-subtract a . b) "mv/scalar ... -> mv"
    (mv-simplify
      (fold
        (l (b result) (if (pair? b) (append (mv-gp b -1) result) (pair (basis-new 0 (- a)) result)))
        null b)))

  (define (mv-null? a)
    "mv -> boolean
     null if empty or all blades have a scale of 0"
    (null? (mv-simplify a)))

  (define (mv-reverse a) (map basis-reverse a))
  (define (mv-grade-inversion a) (map basis-grade-inversion a))

  (define (mv-dual a metric) "mv list ->"
    (mv-ip (mv-new (basis-new (- (bitwise-arithmetic-shift-left 1 (length metric)) 1) 1)) metric))

  (define (mv-scalar a)
    "mv -> scalar
     returns the sum of all scalar values in multivector (not the pseudoscalar)"
    (apply + (filter (l (a) (zero? (basis-id a))) a)))

  (define* (mv-sp a b #:optional metric) "list:mv list:mv [list:metric]-> mv"
    (mv-scalar (mv-ip a b metric)))

  (define (mv-coordinates a) (map tail a))

  (define (mv-norm a)
    "mv -> number
     magnitude"
    ; sum the squares of the scalars and take the square root
    (sqrt (fold (l (a result) (let (a (basis-scalar a)) (+ result (* a a)))) 0 a)))

  (define (mv-type base-names)
    "symbol/integer ... -> procedure:{scalar ... -> mv}
     returns a procedure that when called with scalars returns a multivector with the
     bases specified by base-names"
    (let (bases (map (l (a) (if (symbol? a) (basis-id-from-symbol a) a)) base-names))
      (l a (map basis-new bases a))))

  (define (mv-types-new config) "((any:key symbol ...) ...) -> (procedure ...)"
    (map (l (a) (pair (first a) (mv-type (tail a)))) config))

  (define (space-type s key)
    "alist any value ... -> procedure:{number ... -> mv}
     return a procedure for creating a multivector with a specific combination of bases.
     custom types can be created using the \"#:types\" argument when calling \"space-new\".
     the default types are multivectors for each grade of basis blades.
     the keys for selecting default types is the grade number, eg
     (space-type-get s 2) gets the default bivector type"
    (alist-ref (alist-ref-q s types) key))

  (define (space-basis s name . scalars)
    "alist symbol number ... -> basis/false
     create a basis blade object corresponding to the given space object.
     examples
       (space-basis-new s (q e1) 1)
       (space-basis-new s (q e1e2) 2.5)"
    (and-let* ((id (alist-ref (alist-ref-q s bases) name))) (basis-new id scalars)))

  (define* (space-new metric #:key types)
    "(integer ...) [#:types ((custom-type-name basis-name-or-id ...) ...)] -> alist
     example
       (define s (space-new (list 1 1 1 1)))
       (define vec-new (space-type s 1))"
    (let*
      ( (names-and-ids (map (l (a) (pair (basis-id->symbol a) a)) (basis-ids (length metric))))
        (blades-per-grade
          (group
            ; ignore 0-blade
            (tail names-and-ids) (compose basis-id-grade tail)))
        (ids-per-grade (map (l (a) (pair (first a) (map tail (tail a)))) blades-per-grade)))
      (alist-q bases names-and-ids types (mv-types-new (append ids-per-grade (or types null)))))))
