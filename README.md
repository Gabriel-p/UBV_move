UBV-Move
=============

Moves stars in a UBV color-color diagram to their intrinsic position over a given ZAMS.
Uses an interpolated Schmidt-Kaler ZAMS and the standard extinction law:

    E(U-B) = 0.72*E(B-V) + 0.05*E(B-V) ** 2

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

    ID x  y  mag  e_mag  bv  e_bv  ub  e_ub

Any extra columns will be ignored.

### Running

From terminal the code can be executed with the command:

    python UBV_move.py

Alternatively you can also pass a name for the cluster and a maximum value of `E(B-V)`
extinction to be used  with:

    python UBV_move.py NGC001 3.5

If no name or maximum extinction value are passed the code will assume "cluster" and 4.
respectively by default.

### Output

The code outputs a text file called `clust_name_stars_out.dat` (where `clust_name`
is the name assigned when calling the code) with all the solutions found for each star
which includes extinction, distances, absolute magnitudes and spectral types.

Two figures are also created, the first one is called `clust_name_dens_map.png`,
see **Fig. 1**, and shows a density map of all the extinction and distance values found
for all the observed stars in the CMD in various bin sizes. The maximum density value is
displayed and assumed as the most likely extinction and distance values for the cluster.

![Output 1](/out1.png "Example density map.")
**Fig. 1** Density map of extinction vs distances values showing all possible solutions
for all the observed stars.

The second figure is called `clust_name_CMD.png`, see **Fig. 2**, and shows the observed
CMD (left), the new positions of stars over the ZAMS (center, only for those stars with
a unique solution) and a CMD showing only those observed stars with at least one solution
including extinction and distance values within a given range of the most likely extinction
and distance values assumed for the cluster (right).

![Output 2](/out2.png "Example CMDs.")
**Fig. 2** *Left*: Observed cluster. *Center*: Corrected stars placed over main sequence.
*Right*: Only stars with distance and extinction values close to the maximum density
value in the density map (see below)

