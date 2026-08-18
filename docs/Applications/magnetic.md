# Magnetic

This program calculates the chromomagnetic portion of the field strength tensor, i.e.

$B^2 = -\sum_{i=1}^{3} \mathrm{Tr}\left(F_{jk}F_{jk}\right)$

using a five-loop clover discretisation. The result is averaged over all lattice sites. The real trace of the three spatial field strength tensors $F_{23},F_{21},F_{12}$ are summed.


The code uses a toml as input file and [tomlf](https://github.com/toml-f/toml-f) to read it. The toml file must be of form

```toml

fixNum = 2
fixLabels = ['G2L-8', 'G2L-36']

[fix.G2L-8]
gaugePath='conf/Gen2L/8x32/'
gaugeFormat='openqcd'
cfgList='conf/Gen2L/G2l_8x32.list'
NT=8
NS=32

[fix.G2L-36]
gaugePath='conf/Gen2L/36x32/'
gaugeFormat='openqcd'
cfgList='conf/Gen2L/G2l_36x32.list'
NT=36
NS=32
```

where the `gaugePath` is either the full or relative path to the directory containing the gaugefields, `gaugeFormat` is either `openqcd` or `cssmILDG` (for ILDG binary) and `cfgList` contains a list of the filenames of the gaugefield files, each on a separate line.


The application will read and produce $B^2$ for each gaugefield as well as the total plaquette, spatial and temporal plaquettes. It will perform a 1-remove jackknife analysis using [FJSample](https://github.com/RJaBi/FJsample) and print for each 'label' in `fixLabels`. It will work for different numbers of fixLabels, and do each set sequentially.


Run using 

```
fpm run magnetic -- mytoml.toml
```
