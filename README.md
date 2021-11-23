
### guile-ffi-cblas

This is a set of Guile FFI bindings for **CBLAS**, the library of linear algebra
subprograms<sup id="a1">[1](#f1)</sup>.

To use the bindings, import `(ffi cblas)`. **CBLAS** will be loaded from the
default dynamic library path (see ‘Installation notes’ below). There are up to
three bindings for each function, here using `ZGEMM` as an example:

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
functional binding.

If the function doesn't return an array (e.g. `sdot`) then we only provide
two bindings (e.g. `cblas_sdot` and `sdot`).

Enter `,help (ffi cblas)` at the Guile REPL to list all the bindings
available<sup id="a2">[2](#f2)</sup>.

Note that this package is a work in progress and that there are bugs. For
example, negative strides require specific handling for **CBLAS** and are not
supported yet.

### Installation notes

`guile-ffi-cblas` uses `dynamic-link` to load the dynamic libraries for
**CBLAS**. To do this, the names of the respective library files must be
known. The default name `libcblas` can be configured with the environment
variable `GUILE_FFI_CBLAS_LIBNAME`. Some alternative names that I've seen are
`libblas` or `libatlas`.

If your **CBLAS** libraries are not installed in the default dynamic library
search path, you can configure specific paths for `guile-ffi-cblas` with the
environment variable `GUILE_FFI_CBLAS_LIBPATH`. There are other variables that
control where `dynamic-link` searches for libraries (`LTDL_LIBRARY_PATH`,
`LD_LIBRARY_PATH`) and you may prefer to set those instead. See also
[https://notabug.org/ZelphirKaltstahl/guile-ml#using-guile-ffi-cblas](https://notabug.org/ZelphirKaltstahl/guile-ml#using-guile-ffi-cblas)
for notes on using `guile-ffi-cblas` on Guix<sup id="a3">[3](#f3)</sup>.

### Running the tests

The tests use SRFI-64.

```
$GUILE -L mod -s test/test-ffi-cblas.scm
```

Depending on your installation (see above) you might need

```
GUILE_FFI_CBLAS_LIBNAME=libotherblas \
GUILE_FFI_CBLAS_LIBPATH=/custom/path/lib \
$GUILE ... etc.
```

### Coverage

#### CBLAS level 1

* `srotg` `drotg` `crotg` `zrotg` <sup id="a4">[4](#f4)</sup>
* `sscal` `dscal` `cscal` `zscal` `csscal` `zdscal`
* `sswap` `dswap` `cswap` `zswap`
* `scopy` `dcopy` `ccopy` `zcopy`
* `saxpy` `daxpy` `caxpy` `zaxpy`
* `sdot` `ddot` `cdotu` `zdotu` `cdotc` `zdotc`
* `snrm2` `dnrm2` `scnrm2` `dznrm2`
* `sasum` `dasum` `scasum` `dzasum`
* `isamax` `idamax` `icamax` `izamax`

#### CBLAS level 2

* `sgemv` `dgemv` `cgemv` `zgemv`
* `sger` `dger` `cgeru` `zgeru` `cgerc` `zgerc`

#### CBLAS level 3

* `sgemm` `dgemm` `cgemm` `zgemm`

***

<b id="f1">¹</b> See [http://www.netlib.org/blas/](http://www.netlib.org/blas/) or [https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms). [↩](#a1)

<b id="f2">²</b> Unfortunately this triggers a bug in the current
version. [↩](#a2)

<b id="f3">³</b> A previous version of this library also included **BLIS**
bindings, but now I have moved those to a separate library
(`guile-ffi-blis`). [↩](#a3) .

<b id="f4">⁴</b> Some distributions of `libcblas` do not provide
these. `guile-ffi-cblas` will still work, just without these bindings.
