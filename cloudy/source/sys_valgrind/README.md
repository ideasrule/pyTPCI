## Notes

This directory contains the files needed to generate
an O3 optimized executable with g++ that can be used with valgrind

it includes the option -DNOINIT

* * *

Run valgrind with option --leak-check=full as in
```
valgrind --leak-check=full ./cloudy.exe < test.in > test.out
```

The default setup in this directory will use the g++ compiler.
To use a different compiler, add a line like the following to 
the file Makefile.conf in this directory:
```
CXX = icc
```
This example would use the Intel compiler.

* * *

If you receive the error
> #error "This g++ version cannot compile Cloudy and must not be used!"
it is because you have version 2.96 or 3.4 of g++.  These versions cannot
produce a valid executable.  Update to a more recent version of the compiler.

