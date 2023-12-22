# Read me for collion.cpp

a gas in collisional ionization equilibrium

This program varies the gas kinetic temperature in an environment where collisional 
ionization will dominate.  This is very similar to conditions in the solar corona. 
The gas temperature ranges between 1e3 K and 1e8 K in 0.1 dex steps.  The hydrogen 
density is 1.0 cm-3.  The predicted ionization fractions are placed in the file collion.txt.
This file is formatted similar to Jordan 1969, MNRAS, 142, 501.

The script split-collion.pl will spit collion.txt into separate files for each element.
