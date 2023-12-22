# Read me for template.cpp 

A template for writing your own program to call Cloudy as a subroutine

OVERVIEW:
---------
This program is a template you can use to start writing your own
program that calls Cloudy as a subroutine. You can freely modify
the code inside the try block, but should normally not touch any
of the catch blocks. They are designed to catch any exceptions
that Cloudy may throw. If they were not caught, your code would
crash and your output would be lost!

To demonstrate how you can define an input script for Cloudy,
a single call to cdRead is included which exercises the smoke test.
You should replace that with a series of calls to cdRead to define
your own input script. The call to cdDrive causes the script to be
executed. Since you typically want to calculate a series of models
in the program (e.g. some sort of special grid) a simple loop is
included in the program to demonstrate how this can be done.
You can of course change the loop in any way that suits you,
or even use multiple nested loops.


COMPILATION:
------------
To compile the template into an executable, run the script
```
./complink-template.sh
```
It points to ```../../source``` for the necessary
header files (\*.h) and the static library (libcloudy.a).
This is all you need to do, if you compile Cloudy in source/ with g++.

If you want to build the template with a compiler other than g++,
e.g., clang++, you must also build Cloudy with the same compiler.

If you build Cloudy from one of the compiler-specific subdirectories
(for clang++, source/sys_llvm), you need to make sure that the -L
(path to library) option points to the appropriate directory, while
the SYS_CONFIG macro points to the cloudyconfig.h file therein, for
instance:
```
 clang++ template.cpp \
	-DSYS_CONFIG=\"/path/to/cloudy/source/sys_llvm/cloudyconfig.h\" \
	-I ../../source \
	-L ../../source/sys_llvm -lcloudy \
	-o template.exe
```

This is VERY important.

The header file cloudyconfig.h is produced prior to compilation
of any source files, it holds information on the capabilities of
the compiler, and fine-tunes the building process.  Without the
SYS_CONFIG macro, the compiler will pick up the, if any,
cloudyconfig.h file found in source/ (as directed by -I option),
which may have been produced for a different compiler.  This
could cause mysterious compilation errors.
