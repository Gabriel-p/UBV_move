UBV-Move
=============

[![DOI](https://zenodo.org/badge/14873360.svg)](https://zenodo.org/badge/latestdoi/14873360)

Moves stars in a UBV color-color diagram to their intrinsic position over a
given ZAMS. Uses an interpolated [Schmidt-Kaler (1982)][1] ZAMS and the standard
extinction law:

    E(U-B) = 0.72*E(B-V) + 0.05*E(B-V) ** 2

The code stores all possible solutions for each observed star, ie: all `E(B-V)`
values that move the star over the ZAMS according to the extinction law. Each
extinction solution positions the star at a given distance, with fixed values
for its absolute magnitude and spectral type.

See the [Wiki][2] for details on how to install and use.

________________________________

Please cite as:

```
@misc{perren_2017_ubvmove,
  author       = {Perren, G. I.},
  title        = {UBV_move: version 1.1},
  month        = mar,
  year         = 2017,
  doi          = {10.5281/zenodo.375672},
  url          = {https://github.com/Gabriel-p/UBV_move}
}
```


[1]: http://www.fcaglp.unlp.edu.ar/~egiorgi/cumulos/herramientas/tracks/zams.txt
[2]: https://github.com/Gabriel-p/UBV_move/wiki