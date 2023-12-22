## Notes

Set up to run through Clang static analyzer.
```
  make clean
  scan-build -o scan make
```
-- the make can be run in parallel with the -j <n> argument as usual.

This can be used on systems where scan-build and clang++ are
installed in the same location.  If they are not, as they are no
in Lexington, then use the  
```
--use-analyzer
```
option to specify the location of clang++.  On our machines
the entry
```
scan-build -o scan --use-analyzer /usr/local/bin/clang++ make 
```
will create the output.
