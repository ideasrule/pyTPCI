# Cloudy

Cloudy is an _ab initio_ spectral synthesis code designed to model a wide range
of interstellar "clouds", from H II regions and planetary nebulae, to Active
Galactic Nuclei, and the hot intracluster medium that permeates galaxy clusters.

Cloudy has been in continuous development since 1978, led by Gary Ferland, and
in close collaboration with a number of scientists -- see the
[list of contributors](others.txt).


# Version

The current version of Cloudy is C23, released in 2023.
A summary of what is new is available
[here](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/NewC23).

If you used Cloudy in your research, please cite our most recent
[release paper](https://ui.adsabs.harvard.edu/abs/2023RMxAA..59..327C)

## Brief History

Cloudy recently migrated to a pure git version control system from a
subversion (SVN) system (with limited support for git).
Cloudy had been on a SVN repository for about 15 years, which is still 
maintained as a read-only reference at
[https://trac.nublado.org](https://trac.nublado.org).

The migration was done on 2020 Dec 2 at r14364.
Only the trunk and a few actively maintained branches were migrated.

Previous releases of Cloudy are still available on the SVN site,
as well as tarballs on our
[release folder](https://data.nublado.org/cloudy_releases).


# Directory structure

There are seven directories below this one containing the:
1. atomic, molecular, grains data, as well as SEDs (```data/```);
1. documentation (```docs/```);
1. doxygen setup files (```doxygen/```);
1. a unit test library (```library/```);
1. some helpful scripts (```scripts/```);
1. the source (```source/```);
1. and test suite (```tsuite/```).
The test suite directory, tsuite, has a number of directories below it,
each exercising different aspects of Cloudy.

These directories contain all files needed to build and execute Cloudy.
Each directory has a readme file giving more information on its contents.
It is important to maintain this directory structure when the download is
opened on your computer.


# Building Cloudy

Instructions for building the code on various platforms are available on
[the wiki](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/CompileCode).
Makefiles are provided for most popular compilers (see ```source/sys_*```).


# Documentation

The ```docs``` directory contains Hazy, Cloudy's documentation.


# API

Cloudy's API is described with [Doxygen](https://doxygen.nl).
A precompiled version is available
[online](https://data.nublado.org/doxygen/c22.00).

If you wish to produce a local instance, please follow the instructions in
the ```doxygen/``` directory.


# Contact us

See the project's [website](https://nublado.org) for new versions, bug fixes,
etc.

If you have any questions, please post them on the
[Cloudy user group](https://cloudyastrophysics.groups.io/g/Main/topics).
