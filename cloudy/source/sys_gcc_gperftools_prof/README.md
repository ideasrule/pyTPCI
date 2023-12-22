## Notes

This generates an executable compiled with g++ with
gperftools profiling enabled.
To build enter
```
make
```
at the command prompt.

```
gcc [...] -o myprogram -lprofiler
CPUPROFILE=/tmp/profile ./myprogram
```
See the Google Performance Tools
[repository](https://github.com/gperftools/gperftools).
