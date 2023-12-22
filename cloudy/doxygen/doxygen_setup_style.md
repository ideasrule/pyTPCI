# Doxygen Style

## Setup

The setup of the Doxygen output is contained in file `Doxygen`.
To change the preset options, run
``` 
$ doxywizard Doxyfile
``` 
where `doxywizard` is a GUI program that, at least on Ubuntu, is a standalone
program and needs to be installed in addition to `doxygen`.

The Expert tab has many options.
The Doxygen setup file ignores most of these options, and creates only HTML
output in a new subdirectory `html`, created in this directory.


To run `doxygen`, one can use the `Run` tab on `doxywizard`, or run the command
```
$ cd doxygen
$ doxygen Doxyfile
```
that is, run it from this directory.


## Source Code Markup

Within the codebase, doxygen markup is fully contained in the headers.
All doxygen special comment blocks are using the markup for
[C-like languages](https://www.doxygen.nl/manual/docblocks.html#specialblock).
Some basic rules:
 - Special comment blocks are indicated by ``/**``.
 - All comment blocks are immediately before what they describe, with the
   exception of some struct members, where the comment immediately follows ("/**<" syntax).
 - Commands are indicated with the with the syntax `\command`.
 - Within special comment blocks, the \, @, &, $, #, <, >, % characters must be
   escaped using a preceeding '\'.
 - LaTeX formulas can be entered using the `\f$` sentinel, e.g.,
```
\f$ - same as in line $  - opposite is another \f$
```

The following commands are in use:
 - `\file <filename> (description)` - description for a file, will appear in
   the output above any descriptions of items contained in the file.
 - `\verbatim`  - Doxygen outputs text enclosed in verbatim/endverbatim tags
   exactly (ie, preserving whitespace and newlines)
 - `\endverbatim`
 - `\param [in|out|in,out] <parameter name> (description)` - describe a
   parameter of a function.  Parameter name is the name of the variable and 
   does not include type.  The '[in|out|in,out]' is optional.
 - `\post (description)` - describe the post conditions for a function
 - `\return (description)` - describe what the function returns (in output
   "Returns "+description is printed)
 - `\author` Joe Blow


The ideal declaration should look like the following:
```
 /**
  routine_name This is the long description of what routine_name does, appears after 
	the function name
  \brief this is routine brief description - only 1 line long
  \param iz  a description of what parameter 1 is for
  \param [in] in  a description of what parameter 2 is for
  \param [out] *out description
  \author Joe Blow
  \return explain the return value
 */ 
double routine_name(long int iz, 
  long int in ,
  double *out );
```

The full description of these commands and all others is under the
[Special Commands](https://www.doxygen.nl/manual/commands.html) section of the
online reference.
