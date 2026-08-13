# Fortran Compilation

This package is designed to be compiled using the [fortran-package-manager](https://fpm.fortran-lang.org/), i.e.

```fpm build```

If using the intel compiler (ifx) you can enable openMP style parallelism with do-concurrent using the "-qopenmp" flag as well as ensuring that `LOCALITYSUPPORT` is defined in the `fpm.toml` file. gfortran<15 does not fully support the 2018 Fortran standard and so you need to remove the `LOCALITYSUPPORT` macro. Parallelisation may still work with gfortran using the flag ` -ftree-parallelize-loops=N` but is less likely to perform well. 


The preprocessor options are set in the `[preprocess]` section of the `fpm.toml` file.

```
[preprocess]
[preprocess.cpp]
suffixes = ["F90"]
macros=['SETGITHASH=Yes', 'LOCALITYSUPPORT=1']
```


If you have `SETGITHASH=Yes`, then `src/version.F90` will look for a file `GITHASH.txt` in `src` and put the contents into a variable so that the git commit number can be printed at run-time. This `GITHASH.txt` file may be generated using the `src/generateGitHash.sh` script.


