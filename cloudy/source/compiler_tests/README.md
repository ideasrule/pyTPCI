## Notes

This test suite verifies that illegal constructs are rejected by the compiler.

To run the test suite:
```
  ./run_tests.pl [ <compiler-name> ]
```
The default compiler name is g++. Compiler output will be in run_tests.log.
This file should be checked to see if it contains the expected error messages
(the perl script only checks if the compiler returns an error code). All tests
are expected to fail.
