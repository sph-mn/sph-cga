experimental code for implementing geometric algebras.
see also [learning about conformal geometric algebra 2017](http://sph.mn/c/view/13e).

not thoroughly validated and largely didactic at this point.

current features
* basis blade and multivector objects
* geometric/inner/outer product for basis blades and multivectors
* purely functional
* k-blade pre-generation for the desired number of dimensions
* no macros and optimisations yet

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
sph-math-cga-description
bit-range
bit-count
factorial
blade-id-grade
blade-id-sign
blade-id->string
blade-id-from-string
blade-id-less?
blade-ids
blade-equal?
blade-gp
blade-op
blade-ip-f
blade-ip
blade-grade-inversion
blade-reverse
blade-conjugate
blade-new
blade-id
blade-scalar
blade-count
blade-grade
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
mv-scalar-product
mv-type
mv-types-new
space-new
```

# example usage
```
(import (sph math cga))

(mv-gp
  (mv-new (blade-new 1 3) (blade-new 2 4))
  (mv-new (blade-new 1 5) (blade-new 2 9)))
```
