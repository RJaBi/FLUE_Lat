# OQCD_stoutSmear

Reads an openqcd format gaugefield. Performs some (spatial) [stout-link smearing](https://doi.org/10.1103/PhysRevD.69.054501) and writes out the unsmeared (bare) and smeared average plaquette, spatial plaquette and temporal plaquette


Run using

```
fpm run OQCD_stoutSmear -- inputFile rho nSmear NT NS
```


Where `rho` is the smearing strength and `nSmear` is the number of smearing sweeps. Typical numbers are `rho=0.14` and `nSmear=1`.
