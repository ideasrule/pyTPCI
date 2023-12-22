## Notes

LineLabels.in produces a report of the species and lines predicted
by the code.  The input file contains its own documentation.

This directory also contains PDF versions of the Quick Start Guide and
Parts 1 and 2 of Hazy, the documentation. This Quick Start Guide
gives an overview.  Part 1 gives a complete list of all commands. 
Part 2 describes the output that is generated, shows how to extract 
observed quantities from the predictions, and explains how to call
Cloudy as part of a larger program. 

Part 3, which is badly out of date at the time of this writing, describes
the physics behind the simulation.  The past few years have seen many
expansions of the code's capabilities, especially in infrared emission,
molecular physics, dynamics and advective flows, and time dependent simulations.
Unfortunately, Part 3 of Hazy has not had a high enough
priority to be updated with available resources.  I have not given up and
do intend to update this section eventually.  For now, an ADS search on
"Ferland, G" will find the most recent papers that describe advances in the
physics.

The LaTeX source needed to build the documentation is included in folders
below docs/latex.  The Perl script latex/CompileAll.pl will compile all
four parts of the documentation, leaving PDF files in the appropriate
folders.  To only compile one volume do 
```
CompileAll.pl hazy1
```

To check on compile-time errors do
```
grep underfined */*.log in the latex directory.
```

To erase all intermediate files produced by latex and bibtex
(but not the PDFs), do
```
cleanAll.pl
```
in the latex directory.  This may be needed if major changes to the document
structure are made so that pdflatex becomes confused.  To only clean one
volume do
```
cleanAll.pl hazy1
```
