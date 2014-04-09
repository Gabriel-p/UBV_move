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
same folder where the code is located. The format of the columsn in the file should
follow the order:

    ID x  y  mag  e_mag bv  e_bv ub  e_ub

## Running

From terminal the code can be executed with the command:

    python UBV_move.py

Alternatively you can also pass a name for the cluster and a maximum value of `E(B-V)` extinction to be used  with:

    python UBV_move.py NGC001 3.5

## Output

![Output 1](/out1.png "Example first output image")
**Fig. 1** *Left*: Observed cluster. *Center*: Corrected stars placed over main sequence.
*Right*: Only stars with distance and extinction values close to the maximum density
value in the density map (see below)

![Output 2](/out2.png "Example second output image")
**Fig. 2** Density map of extinction vs distances values showing all possible solutions
for all the observed stars.