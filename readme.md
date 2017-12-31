experimental code for implementing geometric algebras.
see also [learning about conformal geometric algebra 2017](http://sph.mn/c/view/13e).

not thoroughly validated and largely didactic at this point.

current features
* basis blade and multivector objects
* geometric/inner/outer product for basis blades and multivectors
* purely functional
* no macros and performance optimisations yet

# license
gpl3+

# dependencies
* [guile](https://www.gnu.org/software/guile/guile.html) >= 2.2
* [sph-lib](https://github.com/sph-mn/sph-lib)

# installation
copy or link the path ``sph/math/cga.scm`` into a path that is listed in the environment variable $GUILE_LOAD_PATH

# module name
```
(sph math cga)
```

# module exports
```
basis-id-grade
basis-id-sign
basis-id->string
basis-id->symbol
basis-id-from-symbol
basis-id-from-string
basis-id-less?
basis-ids
basis-equal?
basis-gp
basis-op
basis-ip-f
basis-ip
basis-grade-inversion
basis-reverse
basis-conjugate
basis-new
basis-id
basis-scalar
basis-count
basis-grade
mv-new
mv-simplify
mv-op
mv-ip
mv-gp
mv-sum
mv-subtract
mv-null?
mv-reverse
mv-grade-inversion
mv-dual
mv-scalar
mv-sp
mv-type
space-type
space-base
space-new
```

# example usage
```
(import (sph math cga))

(mv-gp
  (mv-new (basis-new 1 3) (basis-new 2 4))
  (mv-new (basis-new 1 5) (basis-new 2 9)))

; default multivector types
(define 2d-metric (list 1 1))
(define s (space-new 2d-metric))
(define vec-new (space-type s 1))
(define biv-new (space-type s 2))
(define tri-new (space-type s 3))
(mv-op (biv-new 2) (biv-new 3))

; custom types
(define s (space-new (list 1 1) #:types (q ((custom s e2 e1e3)))))
(define custom-new (space-type s (q custom)))
(custom-new 1 2 3)
```
