## Notes

This is a collection of scripts that are useful to the Cloudy project.
The scripts themselves have documentation.

adf042stout - convert an ADF04 file to Stout format

--------------

The following scripts derive a publication quality table from the information embedded
in the Stout data files (first two) and a table describing individual species models:

db-ref-bib2json.pl - gather references from Stout species in masterlist, save in data/stout/refs.json
db-ref-json2tex.pl - convert refs.json data to Latex table in current directory
db-species-tex.pl - make table of species and origin of their data

Please read the comments at the top of each script for usage instructions.
These require that the following packages be installed, as explained in db-ref-bib2json.pl 
http://search.cpan.org/dist/Astro-ADS
http://search.cpan.org/~makamaka/JSON-2.90/
http://search.cpan.org/~ambs/Text-BibTeX-0.71/

--------------

generate_md5sums - create sums to check that certain files have not been altered

test_properties check whether the files in the 
Cloudy repo have the correct file properties (such as svn:eol-style). 

test_repos - driver for test_properties

HS-ratio-make.pl - create table of Case B H-like lines to 
compare predictions with Hummer & Storey.  Used to create files
HS-species.dat in tsuite auto.  Marios Chatzikos.
