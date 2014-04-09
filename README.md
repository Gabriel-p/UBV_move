UBV-Move
=============

Moves stars in a UBV color-color diagram to their intrinsic position over a given ZAMS.
We use an interpolated Schmidt-Kaler ZAMS and the usual extinction law.

The code stores all possible solutions for each observed star, ie: all extinction values
that move the star over the ZAMS according to the reddening line. Each solution positions
the star at a given distance. The density map shows *all* the possible solutions for
*all* the observed stars and interprets the maximum value as the most probable extinction
and distance associated to the (putative) cluster.

Dirty code, definitely needs improvement. Works nonetheless.

### Data input

The code expects an input photometric data file named `data_input.dat` to exist in the
same folder where the code is located. The format of the columns in the file should
follow the order:

    ID x  y  mag  e_mag bv  e_bv ub  e_ub

### Running

From terminal the code can be executed with the command:

    python UBV_move.py

Alternatively you can also pass a name for the cluster and a maximum value of `E(B-V)` extinction to be used  with:

    python UBV_move.py NGC001 3.5

### Output

The code outputs a text file called `clust_name_stars_out.dat` (where `clust_name`)
is the name assigned when calling the code) with all the extinction solutions found for each star along with their resolved distances, absolute magnitudes and spectral types.

Two figures are also created: `out1.png`, see *Fig. 1*, which shows the observed CMD, the new positions of stars over the ZAMS (only for those stars with a unique solution) and
`out2.png` which shows a density map of all the extinction and distance values found for
all the observed stars in the CMD in various bin sizes, see *Fig. 2*.

![Output 1](/out1.png "Example first output image")
**Fig. 1** *Left*: Observed cluster. *Center*: Corrected stars placed over main sequence.
*Right*: Only stars with distance and extinction values close to the maximum density
value in the density map (see below)

![Output 2](/out2.png "Example second output image")
**Fig. 2** Density map of extinction vs distances values showing all possible solutions
for all the observed stars.