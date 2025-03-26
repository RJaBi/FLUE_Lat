# FLUE

Code to calculate various gluonic quantities in lattice QCD. These include quantities such as various wilson loops or $F_{\mu\nu}$. See the first table below.



| Method                               | Implementation Status | Notes                                                                      |
|--------------------------------------|-----------------------|----------------------------------------------------------------------------|
| Generic 'path' based wilson line     | &check;               | Returns the multiplication along the path from a starting point            |
| Plaquette (spatial, temporal, total) | &check;               | both generic path and hard coded versions                                  |
| Polyakov Loop                        | &check;               | hard coded only                                                            |
| Clover Fmunu                         | &check; &cross;       | Function exists but is not currently used anywhere                         |
| 5-loop improved Fmunu                | &check;               | hep-lat/0203008                                                            |
| B^2 Calculation                      | &check;               | Currently using 5-loop Fmunu                                               |
| W_munu                               | &check; &cross;       | Deprecated functions to do nMuxnNu wilson loops. Does not use generic path |
| Gluon Propagtor                      | &cross;               | Work in progress.                                                          |
| uzerobar			       | &check;	       | Calculates u_0 using Landau gauge and plaquette definitions		    |


---
## Gauge IO

The supported formats are `openqcd` and an ILDG-like binary format. This is the ILDG binary data format, but in big-endian.

---

Based upon [fortran_meson_py](https://github.com/SalvadorBrandolin/fortran_meson_py) with thanks to the [Fortran-Lang Discourse](https://fortran-lang.discourse.group/t/packaging-a-fpm-project-with-python-bindings-a-little-guide-and-insights-from-our-experience/8495/9)

----
#Run

Run a program using a command like the below. Here magnetic is the name of the app, we set the compiler to ifx, add the -qopenmp flag, specify the 'release' default flags of fpm, and us the mag.toml input file.

fpm run magnetic --compiler ifx --flag "-qopenmp" --profile release -- mag.toml

If using the intel compiler (ifx) you can enable openMP style parallelism with do-concurrent using the "-qopenmp" flag as well as ensuring that `LOCALITYSUPPORT` is defined in the `fpm.toml` file. gfortran does not fully support the 2018 Fortran standard and so you need to remove the `LOCAILTYSUPPORT` macro. Parallelisation may still work with gfortran using the flag ` -ftree-parallelize-loops=N` but is less likely to perform well. Alternatively use the 'OMP' macro to activate openMP parallelism instead.
```
[preprocess]
[preprocess.cpp]
suffixes = ["F90"]
macros=['SETGITHASH=Yes', 'LOCALITYSUPPORT=1']
```
