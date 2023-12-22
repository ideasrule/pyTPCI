Read me for Cloudy data files
=============================

last reviewed November 3, 2022

* * *

Overview
--------

This documents the data files included in the Cloudy data distribution.
The instructions for setting up the code are on the
[web site](https://www.nublado.org).

Three types of files are present in this (top) directory:

* **_\*.ini_** These are files that are used to set up Cloudy.
You can change these files to alter the default behavior of the code.
These can, for instance, change the continuum resolution or add new entries
into the main emission-line output.

* **_\*.dat_** These are the foundation atomic/molecular data files that are
needed for the code to operate.
Do not change these files.

* **_\*.in_** These are input scripts that are used to compile various helper
files such as stellar atmospheres or types of grains. 

* * *


Database atomic/molecular models
--------------------------------

Cloudy draws the majority of its atomic and molecular data from external
databases, namely:

* the _Chianti_ database (v10.0.1),
[Del Zanna et al. 2021 ApJ, 909, 38](https://ui.adsabs.harvard.edu/abs/2021ApJ...909...38D)

* the _Stout_ database,
[Lykins et al. 2015 ApJ, 807, 118](https://ui.adsabs.harvard.edu/abs/2015ApJ...807..118L)

* the LAMDA database,
[van der Tak et al. 2020, Atoms, 8, 15](https://ui.adsabs.harvard.edu/abs/2020Atoms...8...15V)

See the respective papers for details.

The commands
```
database [chianti|stout|lamda]
```
are provided to manipulate each database en masse, while the command
```
species "Fe+" [dataset="Tayal18"]
```
is provided for refined control of individual model ions.
Notice that the keyword ```dataset``` may be used to select between various
datasets, as in the Fe+ example above.


* * *

Compiling ancillary files
-------------------------

The download includes all the files you will need to get Cloudy running.
The download does not include the stellar continua file that are needed for
the table stars command to work.
The web site describes how to set up these continua and this only needs to be
done if you want to use the table stars command.

The script ```make_data.sh``` may be used to construct ancillary files used
by Cloudy.
Some of them, relating to grains, are discussed below.

### (Possibly) compiling the size-distributed grains

It is easy to create new grain tipes.
Compiled opacities are already included in the data files (the \*.opc files),
so nothing need be done if you are happy with the default setup.
They would need to be recompiled if you change the energy grid of the code,
or wish to use a different grain refractive index or size distribution.
Distributed grains are new in C96 and were added in collaboration with
Peter van Hoof and Peter G. Martin.
These use optical properties for each grain material (the \*.rfi files) and
size distribution files (the \*.szd files) to create grain opacities
(the \*.opc files) that are a function of grain material and size.
The code can then do a more realistic simulation of the grain emission and
physics.


### The grain properties files, \*.szd, \*.rfi, \*.opc

These files specify the size distribution (\*.szd) and refractive indices
(\*.rfi) for the new treatment of grain physics that is used with the
```compile grains``` command.
The \*.opc files give the actual opacities used by the code.

For details, see Appendix A in Hazy 1.


* * *

User-defined files that change code behavior
--------------------------------------------

### The \*.ini files

These are a series of files that add commands to the input stream when you run
a model.
They are added by including an **init** command that names one of the following
files.

* **c84.ini** - makes code behave more like version 84  
* **fast.ini** - this includes a command that disables elements to make the code 
run faster, at the expense of a less accurate simulation  
* **honly.ini** - hydrogen only init file  
* **hheonly.ini** - init file for H, He only  
* **ism.ini** - turns off level 2 lines and only includes prominent elements
for depleted mixture.
_**NB**_ This does not add grains to the mix even though many elements are
strongly depleted.
This is not consistent, and so grains should be added separately.  
* **pdr\_leiden.ini** and **pdr\_leiden\_hack.ini** are used to compute the
PDR models given in Roellig et al.

### FeII bands in the output

The data file _FeII\_bands.ini_ is used to specify a series of FeII bands that
are entered into the main emission line output.
These bands are described further in the dat file and in the part of Hazy where
FeII is discussed.


### Continuum bands in the output

The data file _continuum\_bands.ini_ is used to define a series of wavelength
bands.
Each band has a lower and upper wavelength and the code will integrate all
emission over these bands.
The total luminosity or intensity is entered into the main emission line output.


### Emission line lists for cdGetLineList

These are the default lists of emission lines that can be read into arrays of
emission lines by calling cdGetLineList.
This is useful when the code is being called as a subroutine of other larger
codes, as a way to obtain a list of emission lines whose intensities are to
be extracted from the calculation.
This process is described in the section of a later part of Hazy on calling
Cloudy as a subroutine.
These files are meant to be changed by the user. 

The files have names LineList\*.dat and the end of the filename indicates the
purpose.
At present the lists are the following:

* **LineList\_BLR.dat** - lines seen in the spectra of BLR of AGN  
* **LineList\_NLR.dat** - lines seen in the spectra of NLR of AGN  
* **LineList\_HII.dat** - lines for HII Regions  
* **LineList\_strong.dat** - a minimarl list of emission lines  
* **LineList\_PDR.dat** - a list of PDR lines  
* **LineList\_PDR\_H2.dat** - a PDR line list with many H2 lines from the large
molecule  
* **LineList\_HeH.dat** - a set of H and He emission lines  
* **LineList\_He\_like.dat** - lines for the He-like iso sequence

Some line wavelengths may change over time as the accuracy of energy levels
improve.
The ```table lines "name.dat"``` command is used to check that that the lines
in the file name.dat are still correct.
Any line list file that is included here in the data director should also be
included in an one of the test cases in the test suite.


### The resolution of the continuum mesh

The file _continuum\_mesh.ini_ contains data the defines the continuum mesh
used during a calculation.
It is possible to make the continuum have any resolution at all by changing
this file.
The continuum resolution largely sets the execution time so be careful.
The file contains documentation of its use.

* * *

Hummer and Storey Case B data files
-----------------------------------

The files HS\_e1b.dat through HS\_e8b.dat, and HS\_e1a.dat through HS\_e8a.dat
are from
[Storey P.J., Hummer D.G. 1995, MNRAS 272, 41](https://ui.adsabs.harvard.edu/abs/1995MNRAS.272...41S).

* * *

Mewe files
----------

The files mewe\_fluor.dat and mewe\_nelectron.dat are tables 2 and 3 of the
atomic data from 
[Kaastra, J.S., & Mewe, R., 1993, A&AS, 97, 443](https://ui.adsabs.harvard.edu/abs/1993A%26AS...97..443K).

The file mewe\_gbar.dat  is the Mewe data files for g-bar collision strengths.

* * *

Molecular hydrogen files
------------------------

Files related to the structure of the H2 molecule are contained within the
_h2/_ directory.
File names reflect their contents, so that files:

* *energy_?.dat* contain the energy structure of the molecule in various
configurations

* *transprob_?.dat* contain transition probabilities

* *coll_rates_?.dat* contain collision strengths for impact by a variety of
projectiles striking H2

* *dissprob_?.dat* contain dissociation probabilities for various electronic
configurations of H2  
 
* *hminus_deposit.dat* distribution functions for populating levels of H2
following formation from H\-

* * *

UTA data file
-------------

All UTA data are stored in the _UTA/_ directory.
There are two files with names UTA\*.dat, the following:

* UTA\_Gu06.dat - from Gu et al. ApJ, 641, 1227.  These are now the default for Fe0 through Fe+12.  
* UTA\_Kisielius.dat - from Kisielius et al. MNRAS 344, 696, used for Fe+13 - Fe+15.

All other files (303 in total) are taken from the Opacity Project.
For details, see README.md therein.

* * *


Level 1 and Level 2 data files
------------------------------

* level1.dat, the main set of level 1 lines, this file is designed to be edited by humans.  
* level2.dat, the level 2 lines, do not edit this file!

* * *


Miscellaneous directories
-------------------------

Other directories exist with useful information:

* _abundances/_ defines elemental and isotopic abundances used by Cloudy

* _SED/_ contains the Spectral Energy Distributions of various sources 

* * *


Visit [www.nublado.org](https://www.nublado.org) for more details and the latest updates.

Good luck!  
Gary J. Ferland
