experimental code for implementing geometric algebras.
see also [learning about conformal geometric algebra 2017](http://sph.mn/c/view/13e)

not thoroughly validated and largely didactic at this point.

current features
* basis blades and multivectors
* k-blade pre-generation for the desired number of dimensions
* some basic operations on basis blades and multivectors
* no macros and optimisations yet

# license
gpl3+

dependencies
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
blade-conjugate
blade-count
blade-equal?
blade-geometric-product
blade-grade
blade-grade-inversion
blade-id
blade-id->string
blade-id-grade
blade-id-less?
blade-id-sign
blade-ids
blade-inner-product
blade-new
blade-outer-product
blade-reverse
blade-scalar
mv-dual
mv-geometric-product
mv-grade-inversion
mv-inner-product
mv-new
mv-null?
mv-outer-product
mv-reverse
mv-scalar
mv-scalar-product
mv-simplify
mv-subtract
mv-sum
```

# example usage
```
(import (sph math cga))

(mv-geometric-product
  (mv-new (blade-new 1 3) (blade-new 2 4))
  (mv-new (blade-new 1 5) (blade-new 2 9)))
```
