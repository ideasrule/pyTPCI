## Notes

This directory contains the files needed to generate an O3 optimized
executable with g++, using a shared library.  This means that the
executable is small, but needs to be able to access libcloudy.so at
run time.  The libcloudy.so file may also be useful for access by
third-party packages (e.g. perl or python).

To build do
```
source sourceme.txt
```
or enter
```
make
```
at the command prompt.

clean - to remove files that were created by make, do
make clean

The Makefile also builds a python link module using swig.  If swig is
installed, and apart from platform dependencies,

```
% cd gcc_shared
% make -f ../Makefile SRCDIR=.. _cloudy.so
[stuff]
% python
>>> import cloudy
>>> cloudy.cdInit()
>>> cloudy.cdRead("title A model")
4000
[etc...]
% make -f ../Makefile SRCDIR=.. Cloudy.so
[stuff]
% perl
use Cloudy;
print Cloudy::cdInit();
print Cloudy::cdRead("title A model"),"\n";
^D
4000
```
should work.  

Small finesses will allow this to work with Ruby, Tcl, etc. as
advertised on the Swig website.
 
* * *

If you receive the error
> #error "This g++ version cannot compile Cloudy and must not be used!"
it is because you have version 2.96 or 3.4 of g++.  These versions cannot
produce a valid executable.  Update to a more recent version of the compiler.
