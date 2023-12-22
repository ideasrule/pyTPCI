# Cloudy tests readme file

Read me file for the standard tests directory
=============================================

**For more information visit the Cloudy web site,
[www.nublado.org](http://www.nublado.org)**

* * *

The standard tests
------------------

The standard test cases that are computed every night in Lexington are the set
of files with names ending in ".in".
Each contains the commands needed to compute a particular model.
When they are computed they produce output files with the same name but ending
in ".out".
Additional files are created too - these mostly results of all the monitor
commands, ending in ".asr", and overviews of the model results ".ovr".

The file [doc\_tsuite.htm](doc_tsuite.htm) gives a list of all test cases, the
commands contained in each, and a description of its purpose.
This file is automatically created by _doc\_tsuite.pl_

The purpose of each test case is given in the documentation that follows the
input commands.
Cloudy stops reading input commands if the end of file or a blank line is
encountered.
Each ```\*.in``` file begins with the commands, followed by the blank line to
tell the code to stop reading, then followed by a description of the purpose
of the test.
This description is totally ignored - the command parser stops when it
encounters the first empty line.

The first part of the name of many test scripts indicate its purpose.
For instance, all PDR models start with "pdr\_" and all quasar BLR models
start with "blr\_". 

The input scripts follow a particular order of commands although this is not
necessary.
The file template\_scripts.txt describes the various groups of commands and
gives the overall format for an input script.
New scripts should follow this format.

* * *

Asserts - reliability in the face of complexity
-----------------------------------------------

Cloudy is designed to be autonomous and self-aware.
The code uses extensive self-checking to insure that the results are valid.
This  philosophy is described in Ferland 2001, ASP Conference Series, Vol 247,
_Spectroscopic Challenges of Photoionized Plasmas_, G Ferland & D Savin,
editors ([astro-ph/0210161](http://xxx.lanl.gov/ftp/astro-ph/papers/0210/0210161.pdf)). 

Asserts provide the ability to automatically validate complex results.
There are two types of asserts here - the first are a set of commands that are
included in the input files and tell the code what answer to expect, and the
second are C++ language macros that are part of the source and confirm the
internal decisions made by the code.
To avoid confusion, the first type of asserts have been renamed to monitors
since C10.

All of the files in the test suite include **monitor** commands.
The **monitor** command is described in the **Miscellaneous Commands** section
of **Hazy I**, and provides the infrastructure needed for complete automatic
testing of the code.
This command tells the code what answer to expect.
If the computed results do not agree then the code prints a standard comment
and exits with an error condition.
These **monitor** commands have nothing to do with the simulation itself, and
would not be included in an actual calculation.
You should ignore them, or even delete all of them (a Perl script,
tests\_remove\_asserts.pl, is provided as part of the test suite to do this).

The source code also includes many C++ **assert** macros that are designed to
validate the internal decisions made by the code.
The **assert** macro only exists if the NDEBUG macro is not set when the code
is compiled.
(If the NDEBUG macro is set then **assert** macros within the source are
ignored by the compiler.)
The test cases should be run at least once with the **assert** macro active,
that is, do not include a compiler option to define the NDEBUG macro.
In most compilers the NDEBUG macro is set when compiler optimization is set to
a high level.
In practice this means that the entire code should be compiled with only low
optimization and the tests computed to validate the platform.
Then recompile with higher optimization for production runs, and recompute the
test cases.

* * *

The "run" command
-----------------

When executed as a stand-alone program the code expects to read commands from
standard input and write results to standard output.
I compute single models by defining a shell script called "**run**`".`
It is described in more detail in the
[Cloudy wiki](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/RunningC22).

* * *

The Perl scripts
----------------

The test suite includes a series of Perl scripts that compute test cases and
then check for errors.
These provide an automatic way to validate the code.  

### Notes on Perl

Perl comments start with a sharp sign "#" and end with the end of line.
Perl variable names begin with a dollar sign "$".
A Perl script is executed by typing the name of the script, as in
**`run_parallel.pl`** or **`./run_parallel.pl`**

### run\_parallel.pl

This script was written by Peter van Hoof and will run the test suite using a
number of processors.
The beginning of the script says how to use it.

### rerun\_parallel.pl

This script is very similar to run\_parallel.pl, with two important differences.
First, it does not clean up the test suite directory before running, and second,
by default it only runs those sims that failed in an earlier run.
This script is mainly intended for developers.

### checkall.pl

This script searches for problems in all of the test cases (the \*.out files)
in the current directory.
It first looks for botched monitors and warnings.
These indicate very serious problems.
A list of any models with these problems is placed in the file **serious.txt**.
Next it whether the string PROBLEM was printed in any of the models, and writes
a list of these to the file **minor.txt**.
A few of these can occur in a normal series of models and they do not,
by themselves, indicate a serious problem.
Finally it looks for all models that did not end.. 

### doc\_tsuite.pl

This script creates two files that document the test suite.
The file **doc\_tsuite.html** is an html representation of the set of tests,
while **doc\_tsuite.txt** is a comma delimited table that can be incorporated
into a word processor.

* * *

If you find a problem
---------------------

Cloudy should run on all platforms without errors.
Botched monitors or outright crashes should never happen.
I cannot fix it if I do not know it is broken.
Please let me know of any problems.
My email address is [gary@pa.uky.edu](mailto:gary@cloud9.pa.uky.edu)

* * *

Visit [https://nublado.org](https://nublado.org) for details and latest updates.

Good luck,  
Gary J. Ferland

Last edited: Nov 3, 2022
