## Notes

This generates an O3 optimized executable with icc
To build enter
```
make
```
at the command prompt.

sys_icc_mac:

on Mac with xcode 3.2.2 and icc 11.1 add
```
-use-asm
```
to both compile and link steps in Makefile.conf.

