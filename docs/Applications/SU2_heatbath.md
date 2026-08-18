# SU2_heatbath

Run's $SU(2)$ heatbath using the standard Cabibbo-Marinari method for Wilson gauge. Does not save the configurations, just measures the plaquette after each (thermallised) trajectory.


All parameters are set using parameters in the `.f90` file. Notably this is *not* production worthy code as it uses a fixed seed.


Run using


```
fpm run SU2_heatbath
```
