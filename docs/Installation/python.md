# Python Support

In the `python` folder we also supply a `f2py` and `meson` based system that will enable use of (some of) the functionality in python. This can be installed from the `python` folder using

`pip install -e . --no-build-isolation`

or via a conda environment file like:
```
name: snakeTesting
channels:
  - conda-forge
dependencies:
  - conda-forge::python=3.12.9
  - conda-forge::fpm            # Fortran package manager
  - conda-forge::meson-python   # For building fortran for use in python
  - conda-forge::ninja
  - conda-forge::pip
  - conda-forge::numpy
  - pip:
       - ../../libs/FLUE_Lat/python
```
where the repository has been cloned inside the `../../libs` directory (this can also be done with `git submodule`).

Available routines are:

IO Routines:

`write[ILDG/OQCD](filename, U, NS, NT)`
`read[ILDG/OQCD/CSSM](filename, NS, NT)`

plaquette routines:

`plaq(U)/splaq(U)/tplaq(U)` returns the average of all plaquettes, space-space plaquettes and space-time plaquettes respectively.

'gauge' routines:

`stoutSmearLinks(U, rho, nweeps)` returns a new copy of the gaugefield which has been stout-linked smeared (spatial links).

The IO routines return a numpy array of type double precision complex, of shape `(NT, NS, NS, NS, 4, 3, 3)`. Note this is not the same as the Fortran, the translation is done in `c_wrapper/FLUE_c.f90`. This is so that it matches i.e. [lyncs_io](https://github.com/Lyncs-API/lyncs.io).

This `f2py` and `meson` build approach is inspired by [fortran_meson_py](https://github.com/SalvadorBrandolin/fortran_meson_py/) by Salvador Brandolin with thanks to the [Fortran-Lang Discourse](https://fortran-lang.discourse.group/t/packaging-a-fpm-project-with-python-bindings-a-little-guide-and-insights-from-our-experience/8495/9)
