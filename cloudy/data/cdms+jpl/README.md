## Notes

This directory contains a selection of molecular line lists that were
calculated with the CALPGM program suite written by Herbert Pickett and which
are hosted on the Cologne Database for Molecular Spectroscopy (CDMS) and the
Jet Propulsion Laboratory (JPL). These files are being redistributed with
permission.

Files with names c???0??.cat come from JPL.
Files with names c???5??.cat come from CDMS.

* * *

Interfaces to the respective databases can be found here:

* [CDMS](https://cdms.astro.uni-koeln.de)
* [JPL](http://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.html)

------------------

Pages giving the partition function for each molecule can be found here:

https://cdms.astro.uni-koeln.de/classic/entries/partition_function.html
http://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.cat

<!--
N/A anymore

These are needed to convert the line intensity given in the CALPGM files into
Einstein-A coefficients. The conversion formulas are here:

http://www.astro.uni-koeln.de/cdms/catalog#equations
-->

------------------

The CALPGM program suite written by Herbert Pickett can be found here:

https://cdms.astro.uni-koeln.de/classic/pickett

------------------

These are the papers describing the databases. Each publication using results
from either database should cite the appropriate papers.

CDMS:

H. S. P. Müller, S. Thorwirth, D. A. Roth, and G. Winnewisser; Astron.
Astrophys. 370 (2001) L49 – L52.

H. S. P. Müller, F. Schlöder, J. Stutzki, and G. Winnewisser, THE COLOGNE
DATABASE FOR MOLECULAR SPECTROSCOPY, CDMS: A USEFUL TOOL FOR ASTRONOMERS AND
SPECTROSCOPISTS, J. Mol. Struct. 742, 215–227 (2005).

JPL:

H. M. Pickett, R. L. Poynter, E. A. Cohen, M. L. Delitsky, J. C. Pearson, and
H. S. P. Muller, "Submillimeter, Millimeter, and Microwave Spectral Line
Catalog," J. Quant. Spectrosc. & Rad. Transfer 60, 883-890 (1998).

------------------

The script get_database_files.pl contained in this directory will download all
the necessary data files (including the partition function pages) and will
create the file masterlist. Please EDIT THIS SCRIPT to add or remove molecules
to the masterlist. The script will need the program wget to be installed. To
execute simply type:

make squeaky-clean
./get_database_files.pl

This only needs to be done when there have been updates, or when you want to add
new molecules or isotopologues.

In order to convert the CALPGM format files into LAMDA compatible files that can
be read by Cloudy, simply type:

make [ -j <n> ]

In order to remove the converted database files (and the binary) simply type

make clean

