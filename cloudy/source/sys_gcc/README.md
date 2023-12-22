## Notes

This directory contains the files needed to generate
an O3 optimized executable with g++

===========================

if you receive the error
> #error "This g++ version cannot compile Cloudy and must not be used!"
it is because you have version 2.96 or 3.4 of g++.  These versions cannot
produce a valid executable.  Update to a more recent version of the compiler.

===========================

To change Makefile so that the default build has other properties:

Edit the third line of the Makefile in sys_gcc to read
```
<TAB>$(MAKE) -f ../Makefile SRCDIR=.. $(MAKECMDGOALS) options
```
where "options" contains the options you want to pass on to make, e.g.,
```
-j 4 EXTRA="-DFOO -DBAR"
```
