## Notes

This generates a g++ executable with bounds-checking enabled, requires g++ 4.8.0 or later.

to build enter
```
make [ -j <n> ]
```
at the command prompt

the three tests included here should crash due to array bounds exceeded
```
crashBoundsStack.in
crashBoundsStatic.in
crashBoundsVector.in
```
llvm 3.2 supports this, but you need the compiler-rt 
runtime libraries for the program to link. You can use sys_gccBounds with 
```
make -j <n> CXX=clang++ 
```
to produce a binary with bounds checking enabled. 

If, when using clang++, the error message produced gives the error
locations as a list of hexadecimal pointers, and there is a message
'Trying to symbolize code, but external symbolizer is not
initialized!', you may need to re-run the executable as e.g.
```
ASAN_SYMBOLIZER_PATH=$(which llvm-symbolizer-3.4) ../cloudy.exe
```
in order to get a proper stack trace.
