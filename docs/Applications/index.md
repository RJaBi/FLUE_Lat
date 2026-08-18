# Applications

FLUE contains a number of applications. These are shown below. They may be run using `fpm`:


```
fpm run [Application Name] -- [Arg1] [Arg2] [Arg3]...
```


| Program            | purpose					                                    	                                     | args		                     |
|--------------------|-------------------------------------------------------------------------------------------------------|-----------------------------|
| magnetic           | Calculates the magnetic portion of Fmunu                                                              | mag.toml		                 |
| ILDG_to_OQCD       | Convert a gauge field ILDG-bin in big endian to openqcd                                               | inputFile outputFile NT NS  |
| UNIT_to_OQCD       | Write a unit gaugefield (i.e. colour matrix is Real 3x3 Identity) to openqcd                          | inputFile outputFile NT NS  |
| CSSM_to_OQCD       | Convert a gauge field cssm to openqcd                                                                 | inputFile outputFile NT NS  |
| OQCD_to_ILDG       | Convert a gauge field openqcd to ILDG-bin in big endian                                               | inputFile outputFile NT NS  |
| SU2_HKLS_to_CSSM   | Convert a $SU(2)$ gaugefield in HKLS to cssm                                                          | inputFile outputFile NT NS  |
| SU2_HKLS_to_NRQ2CD | Convert a $SU(2)$ gaugefield in HKLS to NRQ2CD                                                        | inputFile outputFile NT NS  |
| OQCD_stoutSmear    | Stout smear (spatial) an openqcd gauge field and print average unsmeared and smeared plaquette values | inputFile rho nSweeps NT NS |
| SU3_heatbath       | Generate $SU(3)$ gaugefields and measure plaquette. Hard-coded parameters in this file                |                             |
| SU2_heatbath       | Generate $SU(2)$ gaugefields and measure plaquette. Hard-coded parameters in this file                |                             |
