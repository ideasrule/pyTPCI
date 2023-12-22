## Notes

UTA data sets from several papers, cited in \*.dat files:

* UTA\_Gu06.dat
* UTA\_Kisielius.dat

------------------


### NOTATION FOR THE OPACITY PROJECT DATA, nrb00\*

These sets of inner shell transitions are from Nigel Badnell's
[web site](http://amdpp.phys.strath.ac.uk/tamoc/DATA/PE).

[Badnell et al. (2005)](https://adsabs.harvard.edu/abs/2005MNRAS.360..458B)
describe the methods.
They were generated as part of the Seaton/Opacity Project and follow their
convention.
The names designate the electron target because the R-matrix
radiative data is generated from an electron collision problem,
taken down to bound states. 

The names begin with ```nrb00```, followed by the isoelectronic sequence of
the next higher ion, followed by the element.
The number after the element is related to the charge, but think of the number 
as Roman since target electron charge = Roman recombined.
So, in this convention, data for O VI, which is Li-sequence, would be He since
that is the parent ion, followed by o6.  So,
```O VI == he_o6ic1-2.dat```.


Data within the file:

IRSL are the (non-energy-order) levels in
```
   INDX  IRSL   CODE                 S L   WJ        WNR
```
The first column here is lower level, second upper level.
```
  IRSL  IRSL  N   L      DEL(RYD)       B(SEC)      R(SEC)      A(SEC):       1
```
B is Einstein up, R is down and A is Auger.
see [here](https://www.adas.ac.uk/man/appxa-38.pdf)
for the rest of the file spec.

The full BibTeX record of that paper is:
```
@ARTICLE{2005MNRAS.360..458B,
   author = {{Badnell}, N.~R. and {Bautista}, M.~A. and {Butler}, K. and 
	{Delahaye}, F. and {Mendoza}, C. and {Palmeri}, P. and {Zeippen}, C.~J. and 
	{Seaton}, M.~J.},
    title = "{Updated opacities from the Opacity Project}",
  journal = {\mnras},
   eprint = {arXiv:astro-ph/0410744},
 keywords = {atomic processes, radiative transfer, stars: interiors},
     year = 2005,
    month = jun,
   volume = 360,
    pages = {458-464},
      doi = {10.1111/j.1365-2966.2005.08991.x},
   adsurl = {http://adsabs.harvard.edu/abs/2005MNRAS.360..458B},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
