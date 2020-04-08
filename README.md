
### guile-ffi-cblas

This is a set of Guile FFI bindings for two libraries of linear algebra
subprograms, **CBLAS** and **BLIS**. They provide operations such as vector dot
product, matrix-vector product, matrix-matrix product, and so on.

The bindings for either library are independent of each other; you do not need to
have **BLIS** installed to use the **CBLAS** bindings or viceversa. I am
packaging them together because they share a fair amount of code.

**CBLAS** (or **BLAS**) is by far the more popular library<sup id="a1">[1](#f1)</sup>,
and it probably has the fastest implementation on your system. However **BLIS**
is much more regular and nicer to program for.  These bindings are for the **BLIS**'
‘typed’ API.<sup id="a2">[2](#f2)</sup>

Both sets of bindings work in the same way. Upon importing `(ffi cblas)` (or
`(ffi blis)`), **CBLAS** (or **BLIS**) will be loaded from the default dynamic
library path (see ‘Installation notes’ below). There are up to three bindings
for each function, here using `ZGEMM` as an example:

- `cblas_zgemm` (raw binding): the raw C function by `pointer->procedure`. Don't
  use this if you aren't familiar with Guile's FFI.

- `zgemm!` (typed binding): takes array arguments of type `'c64` and operates by
  effect, without making copies. All the arguments must be properly sized. The
  return value is unspecified.

- `zgemm` (functional binding): takes array arguments of compatible types and
  returns a newly constructed array. The arguments will be converted as
  necessary, which may result in copies.  The returned array will be of `'c64`
  type.

In principle, for the last two bindings, you don't need to care whether your
array is row-major or column-major or what the strides are. The bindings will
extract the required strides from the array arguments. However, since
**CBLAS** doesn't support arbitrary strides (e.g. it only supports a column
stride for matrix arguments, assuming column-major order), some array arguments
will cause the typed binding to fail, or result in extra copies with the
functional binding. **BLIS** doesn't have this problem.

If the function doesn't return an array (e.g. `cdot`) then we only provide
two bindings (e.g. `cblas_cdot` and `cdot`).

The bindings also provide type generic versions of the functions (e.g. `dotv`
for **BLIS** `sdotv ddotv cdotv zdotv`). These simply call one of the
typed variants according to the type of the first array argument.

Enter `,help (ffi cblas)` (or `,help (ffi blis)`) at the Guile REPL to list all
the bindings available.

Note that this package is a work in progress and that there are bugs. For
example, negative strides require specific handling for **CBLAS** and are not
supported yet. **BLIS** doesn't mess up negative strides, so it avoids
this problem.

<b id="f1">¹</b> See [http://www.netlib.org/blas/](http://www.netlib.org/blas/) or [https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms). [↩](#a1)

<b id="f2">²</b> See [(https://github.com/flame/blis/blob/master/docs/BLISTypedAPI.md)](https://github.com/flame/blis/blob/master/docs/BLISTypedAPI.md). [↩](#a2)

### Installation notes

`guile-ffi-cblas` uses `dynamic-link` to load the dynamic libraries for
**CBLAS**/**BLIS**. To do this, the names of the respective library files must
be known. The default names `libcblas`/`libblis` can be configured with the
environment variables `GUILE_FFI_CBLAS_LIBNAME`/`GUILE_FFI_BLIS_LIBNAME`. The
name for **BLIS** is a safe bet but the name for **CBLAS** is sadly variable;
possible alternative names are `libblas` or `libatlas`.

If your **CBLAS**/**BLIS** libraries are not installed in the default dynamic
library search path, you can configure specific paths for `guile-ffi-cblas` with
the environment variables
`GUILE_FFI_CBLAS_LIBPATH`/`GUILE_FFI_BLIS_LIBPATH`. There are other variables
that control where `dynamic-link` searches for libraries (`LTDL_LIBRARY_PATH`,
`LD_LIBRARY_PATH`) and you may prefer to set those instead. See also
[https://notabug.org/ZelphirKaltstahl/guile-ml#using-guile-ffi-cblas](https://notabug.org/ZelphirKaltstahl/guile-ml#using-guile-ffi-cblas) for
notes on using `guile-ffi-cblas` on Guix.

### Running the tests

The tests use SRFI-64.

```
$GUILE -L mod -s test/test-ffi-cblas.scm
$GUILE -L mod -s test/test-ffi-blis.scm
```

Depending on your installation (see above) you might need

```
GUILE_FFI_CBLAS_LIBNAME=libotherblas \
GUILE_FFI_CBLAS_LIBPATH=/custom/path/lib \
$GUILE ... etc.
```

and similarly for `BLIS`.

### Coverage

---

#### CBLAS level 1

* srotg drotg crotg zrotg <sup id="a3">[3](#f3)</sup>
* sscal dscal cscal zscal csscal zdscal
* sswap dswap cswap zswap
* scopy dcopy ccopy zcopy
* saxpy daxpy caxpy zaxpy
* sdot ddot cdotu zdotu cdotc zdotc
* snrm2 dnrm2 scnrm2 dznrm2
* sasum dasum scasum dzasum
* isamax idamax icamax izamax

<b id="f3">³</b> Some distributions of `libcblas` do not provide these. `guile-ffi-cblas` will still work, just without these bindings. [↩](#a3)

#### CBLAS level 2

* sgemv dgemv cgemv zgemv
* sger dger cgeru zgeru cgerc zgerc

#### CBLAS level 3

* sgemm dgemm cgemm zgemm

---

#### BLIS level 1

* sdaxpy ddaxpy cdaxpy zdaxpy
* sdaxpby ddaxpby cdaxpby zdaxpby
* sdotv ddotv cdotv zdotv

#### BLIS level 2

* sgemv dgemv cgemv zgemv
* sger dger cger zger

#### BLIS level 3

* sgemm dgemm cgemm zgemm
