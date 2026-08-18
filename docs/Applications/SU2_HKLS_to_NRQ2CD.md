# SU2_HKLS_to_NRQ2CD

Reads a $SU(2)$ Hands-Kim-Lawlor-Skullerud (HKLS) format gaugefield. Writes it to NRQ2CD format.


HKLS format is used by [su2hmc](https://doi.org/10.5281/zenodo.12910604). The NRQ2CD format is used by [Seyong Kim](https://orcid.org/0000-0002-2102-7398)'s non-relativistic QC$_2$D [code](https://doi.org/10.1016/j.physletb.2012.04.002).


Run using


```
fpm run SU2_HKLS_to_NRQ2CD -- inputFile outputFile NT NS [endianess=big|little(default)]
```


Note that here we can specify the endianess of the input gaugefield. The default is little endian as used by the $C$ version of the `su2hmc` code.
