# SU3_heatbath

Run's $SU(3)$ heatbath using the standard Cabibbo-Marinari method for Wilson/Iwasaki/Symanzik gauge with optional anisotropy. Does not save the configurations, just measures the plaquette after each (thermallised) trajectory.


All parameters are set using parameters in the `.f90` file. Notably this is *not* production worthy code as it uses a fixed seed.


Run using


```
fpm run SU3_heatbath
```



Options for actionTag are 'Wilson', 'Symanzik' and 'Iwasaki'. Note that Symanzik and Iwasaki still use temporal rectangles even for anisotropic actions.
