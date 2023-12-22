## Notes

The files in this directory specify standard chemical compositions and,
optionally, grains.
One of these files will be used if the command
```
abundances "filename.abn"
```
appears.

These ```.abn``` files may also specify grains, although grains are not
included by default.
```default.abn``` specifies our default composition.
```default-reference.abn``` is a copy of that file.

The ```\*.dpl``` files specify depletion factors that account for the elements
that are used to build grains.
These may be parsed by the
```
metals deplete "filename.dpl"
```
or
```
metals deplete Jenkins 2009 Fstar=0.5 "filename.dpl"
```
commands.
The default file, described below, will be used if no file appears between the quotes.

* * *

### Specifying abundances

The default is the "solar" composition in ```default.abn```.
It is used if no others are specified.
```default-backup.abn``` is a backup of this file, used to reestablish default
if you change it.
This is a fairly old set of "solar" abundances that we maintain for backwards
compatibility.

To change the default abundance set, copy one of the abundance files (e.g.,
```default.abn``` or ```solar84.abn```) in your working directory, preferrably
with a custom name, edit it at will, and use it with the command
```
abundances "mycustom.abn"
```

#### Available files

* solar84.abn - abundances used by "old solar" option, used in versions 84-94
of the code.

* Other .abn files - see the comments within the file for more details.
The file name gives a good indication of its contents

* * *

### Default depletion files

The command
```
metals deplete
```
by default will read the ISM_CloudyClassic.dpl files.
The format should be clear.
If the element is not specified, then its depletion is set to unity.

The command
```
metals deplete Jankins2009 Fstar
```
reads the ISM_Jenkins09_Tab4.dpl file instead.
If this file is updated or replaced then its current format must be maintained.

Chamani Gunasekera created the spreadsheet from Table 4 of
[Jenkins et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJ...700.1299J).
S does not have data in Table 4, so Chamani used the information provided in
Section 9.


* * *

### Format of the \*.abn files

Comments begin with #.

The file parsing ends with a line containing a field of stars, as
```
***************
```
Everything after the field of stars is ignored, so remaining lines can be used
to document the abundances, in addition to comment lines starting with #.

The grains command can be included.
The abundance files recognize all grain command options.
This makes it possible to include grains in your mix of gas and dust.

Names of elements are followed by the abundance by number relative to hydrogen.

Hydrogen may be in specified, or may not be.
If the abundance of hydrogen is given and is not equal to unity the code will
rescale the abundances by the entered value for hydrogen. 
If hydrogen is not specified it is assumed to have an abundance of 1.

All elements heavier than hydrogen are turned off before the ```.abn``` files
are read.
An element is turned on if it appears in the abn file.
If an element is not included in the ```.abn``` file, it will not be included
in the calculation.

* * *

Original coding by Joshua Schlueter, NSF REU 2012 Summer
