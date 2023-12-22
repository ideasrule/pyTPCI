## Notes

This generates an executable compiled with g++ with
gprof profiling enabled.
To build enter
```
make
```
at the command prompt.

This builds a gcc executable with
coverage analysis switched on, which might be some help here.  Takes a
long while to compile, but seems to run not too much more slowly than
the standard version.  

The runs generate a load of files with ".gcda" extension , which build
up stats for all subsequent runs until they are deleted (so one could
look at coverage for the whole test suite).  

```make reset``` will remove all these statistics files.

```make analyze``` runs "gcov" itself, sending the coverage statistics for
all the routines calls (a lot of text...) to "coverage", and
generating ".gcov" files in the subdirectory which give line-by-line
execution counts.  It uses this output to print some basic coverage
statistics.

The following are useful for analysing the output from gprof
>   grep "^ *-: *[1-9]" *.gcov | wc
>   grep "^ *#" *.gcov | wc
>   grep "^ *[1-9]" *.gcov | wc
