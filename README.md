[![Fortitude](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/PlasmaFAIR/fortitude/main/docs/assets/badge/v0.json)](https://github.com/PlasmaFAIR/fortitude)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://pre-commit.com/)
# <img width="30" height="30" alt="FLUE Favicon" src="https://github.com/user-attachments/assets/407cb74b-5fdf-4ab4-b46c-79cee3e57c85" />LUE
<p align="center">
  <img width="3845" height="1027" alt="FLUE banner" src="https://github.com/user-attachments/assets/afea2fee-3099-41f6-b1d0-1634a8d89e6e" />
</p>

Code to calculate various gluonic quantities in lattice QCD. These include quantities such as various Wilson loops or $F_{\mu\nu}$. See the first table below. Primarily this supports $SU(3)$, but there is a limited amount of $SU(2)$ support.

**Note you may need to set the stacksize to larger (unlimited is easiest) to use this code**


| Method                               | Implementation Status | Notes                                                                      |
|--------------------------------------|-----------------------|----------------------------------------------------------------------------|
| Generic 'path' based Wilson line     | &check;               | Returns the multiplication along the path from a starting point            |
| Plaquette (spatial, temporal, total) | &check;               |                                                                            |
| Polyakov Loop                        | &check;               | hard coded only                                                            |
| Clover Fmunu                         | &check; &cross;       | Function exists but is not currently used anywhere                         |
| 5-loop improved Fmunu                | &check;               | hep-lat/0203008                                                            |
| B^2 Calculation                      | &check;               | Currently using 5-loop Fmunu                                               |
| Wmunu                                | &check; &cross;       | Deprecated functions to do nMuxnNu Wilson loops. Does not use generic path |
| Gluon Propagator                      | &cross;               | Work in progress - POSTPONED.                                              |
| uzerobar			                   | &check;	           | Calculates u_0 using Landau gauge and plaquette definitions		        |
| stout-link smearing                  | &check;               | Stout smear (Morningstar & Peardon, [10.1103/PhysRevD.69.054501](https://doi.org/10.1103/PhysRevD.69.054501)) gauge links. Spatial link smearing only                                                                      |
| $SU(3)$ & $SU(2)$ heatbath           | &check;                | Heatbath algorithm to generate quenched $SU(3)$ or $SU(2)$. Naive implementation (no overrelaxation) but can use Wilson or tree-level Symanzik (Luescher-Weisz action). No improved action for $SU(2)$.                   |

Note that $SU(2)$ supports only the `Generic path based Wilson line` and the `generic path plaquette`. $SU(2)$ modules and routines have SU2 in the name.

---
## Gauge Formats $(SU3)$

The supported formats are `openqcd` and an ILDG-like binary format. This is the ILDG binary data format in big endian. The `cssm` format stores some metadata, and then only the first two rows of the $SU(3)$ matrix. The `cssm` format assumes big endian.

| Format   | Read    | Write   |
|----------|---------|---------|
| openqcd  | &check; | &check; |
| ildg-bin | &check; | &cross; |
| cssm     | &check; | &cross; |


The internal representation of a $SU(3)$ gauge field is `(3, 3, 4, NT, NS, NS, NS)`.

---

---
## Gauge Formats $(SU2)$

The `cssm` format stores some metadata, the real part of the whole $SU(2)$ matrix, the imaginary part of the whole $SU(2)$ matrix, then some more metadata. The `cssm` format assumes big endian.

The HKLS (Hands-Kim-Lawlor-Skullerud) supports the compact format used by [su2hmc](https://doi.org/10.5281/zenodo.12910604).

The NRQ2CD format is `(ip5d2,2,2,4,0:1)` where `ip5d2` is $NX \times NY \times NZ \times NT / 2` where 0 means odd sites and 1 is even sites with NX iterating fastest.

| Format   | Read    | Write   |
|----------|---------|---------|
| HKLS     | &check; | &check; |
| cssm     | &cross; | &check; |
| NRQ2CD   | &check; | &check; |

The internal representation of a $SU(3)$ gauge field is `(2, 2, 4, NT, NS, NS, NS)`.

----
# Run

Run a program using a command like the below. Here magnetic is the name of the app, we set the compiler to ifx, add the -qopenmp flag, specify the 'release' default flags of fpm, and us the mag.toml input file.

fpm run magnetic --compiler ifx --flag "-qopenmp" --profile release -- mag.toml

If using the intel compiler (ifx) you can enable openMP style parallelism with do-concurrent using the "-qopenmp" flag as well as ensuring that `LOCALITYSUPPORT` is defined in the `fpm.toml` file. gfortran<15 does not fully support the 2018 Fortran standard and so you need to remove the `LOCALITYSUPPORT` macro. Parallelisation may still work with gfortran using the flag ` -ftree-parallelize-loops=N` but is less likely to perform well. Alternatively use the 'OMP' macro to activate openMP parallelism instead. The openMP parallelism is limited to `genPlaquette` measurements.
```
[preprocess]
[preprocess.cpp]
suffixes = ["F90"]
macros=['SETGITHASH=Yes', 'LOCALITYSUPPORT=1']
```

| Program            | purpose					                                    	                                               | args		                     |
|--------------------|-------------------------------------------------------------------------------------------------------|-----------------------------|
| magnetic           | Calculates the magnetic portion of Fmunu                                                             | mag.toml		                 |
| ILDG_to_OQCD       | Convert a gauge field ILDG-bin in big endian to openqcd                                               | inputFile outputFile NT NS  |
| CSSM_to_OQCD       | Convert a gauge field cssm to openqcd                                                                 | inputFile outputFile NT NS  |
| OQCD_to_ILDG       | Convert a gauge field openqcd to ILDG-bin in big endian                                               | inputFile outputFile NT NS  |
| SU2_HKLS_to_CSSM   | Convert a $SU(2)$ gaugefield in HKLS to cssm                                                          | inputFile outputFile NT NS  |
| SU2_HKLS_to_NRQ2CD | Convert a $SU(2)$ gaugefield in HKLS to NRQ2CD                                                        | inputFile outputFile NT NS  |
| OQCD_stoutSmear    | Stout smear (spatial) an openqcd gauge field and print average unsmeared and smeared plaquette values | inputFile rho nSweeps NT NS |
| SU3_heatbath       | Generate $SU(3)$ gaugefields and measure plaquette. Hard-coded parameters in this file                |                             |
| SU2_heatbath       | Generate $SU(2)$ gaugefields and measure plaquette. Hard-coded parameters in this file                |                             |

----
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

## Tests
Some of the functionality has tests in the `test` folder using the `testdata` reference data. Run the tests using `fpm test`

---
### AI Assistance
The tests were predominatly generated using Microsoft Copilot (GPT 5.5). Refinements to the heatbath code also used Microsoft Copilot.
