## Notes

This generates an MPI parallelized executable with mpiCC (based on g++)
for running a large grid.

To build enter
```
make
```
at the command prompt, to build debug do
```
make debug
```
add "-j n" where n is the number of threads to use in the build.

This will not work on a Mac, since it uses mpiCC which is misinterpreted as
mpicc on the Mac file system.  Macs use ++ rather than CC to indicate C++.
Try replacing mpiCC with mpic++

to use the MacPorts openmpi mv Makefile_conf.Mac to Makefile.conf

