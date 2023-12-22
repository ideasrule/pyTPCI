## Notes

Read me file for the sample programs

These are a series of programs that show how Cloudy can be
used as a subroutine of other, larger, programs.  There are
two ways to use these programs.  You can create a library, compile
the new main program, and link the main program to this library. 
Or you can compile the program and then link it to all of the
object files created when Cloudy itself was compiled. 

Each sample program includes the header files cddefines.h and cddrive.h. 
These are included in the Cloudy source directory.  You will
need to alter these lines so that they point to where the files live
on your machine or you will need to include the path to these files
on your compile line. 

The perl script run_programs.pl can be set up to use any of the compilers
supported by source/sys_xxx.  It reads a list of program directories
from the file run_programs.dat, builds each by linking against the
object files in source directory, and runs the programs.
See the instructions at the top of that script to find out more.

The shell script 
```
complink.sh directory [gcc/clang] 
```
will use gcc to build the cpp program within directory by default.
If a compiler is provided uses the given compiler and assumes that
cloudy was build in the sys subdirectory.

* template.cpp  -  This is a template you can use to write your own
program that calls Cloudy as a subroutine.  This provides a
template C++ main program that includes the needed exception handlers.

* collion  -  This is an example of a gas in collisional ionization
equilibrium.  There is no ionizing radiation.  This is something
like the solar corona.

* collcool - cooling in collisional ionization

* comp4 - This checks that the code is properly initialized when
it sets up.  Two different models are computed twice and the
output goes to file1.txt (the first pair) and file2.txt (the
second).  These files should be exactly the same.

* hazy_coolingcurve - This is an example of a gas in collisional
ionization equilibrium in which the temperature is varied and
the heating and cooling derived. 

* hazy_kmt - This redoes the "S-curve" calculation originally
done by Krolick, McKee, & Tarter (1981; ApJ, 249, 422), 

* hizlte - a high metallicity gas in thermodynamic equilibrium.

* mpi - This is my main routine for running large grids on our HP
cluster.  It uses MPI to place a grid point on each processor.

* vary_nete - Both the electron density and temperature are varied
over a broad range and several [O III] lines are predicted. 
This shows how to make diagnostic diagrams where line ratios
are used to deduce temperature and density.

* varyn - This is a series of calculations in which the density
is varied and the intensities of some [OII] lines relative to
[OII] 3727 are printed. 

Visit [https://nublado.org](https://nublado.org)  for details and
latest updates. 
