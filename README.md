# FLUE

Code to calculate the gluon propagator. Based upon [fortran_meson_py](https://github.com/SalvadorBrandolin/fortran_meson_py) with thanks to the [Fortran-Lang Discourse](https://fortran-lang.discourse.group/t/packaging-a-fpm-project-with-python-bindings-a-little-guide-and-insights-from-our-experience/8495/9)


app/magnetic
* Calculate F_{ij}^2
app/uzerobar
* Calculates u_0 from Landau fixed configs and the plaquettes


----
#Run

Run a program using a command like the below. Here magnetic is the name of the app, we set the compiler to ifx, add the -qopenmp flag, specify the 'release' default flags of fpm, and us the mag.toml input file.

fpm run magnetic --compiler ifx --flag "-qopenmp" --profile release -- mag.toml

If using the intel compiler (ifx) you can enable openMP style parallelism using the "-qopenmp" flag as well as ensuring that `LOCALITYSUPPORT` is defined in the `fpm.toml` file. gfortran does not fully support the 2018 Fortran standard and so you need to remove the `LOCAILTYSUPPORT` macro. Parallelisation may still work with gfortran using the flag ` -ftree-parallelize-loops=N` but is less likely to perform well.
```
[preprocess]
[preprocess.cpp]
suffixes = ["F90"]
macros=['SETGITHASH=Yes', 'LOCALITYSUPPORT=1']
```
