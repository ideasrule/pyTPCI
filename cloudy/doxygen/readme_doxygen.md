# Readme for Doxygen

This directory contains the file needed to create Doxygen documentation for the
Cloudy source.  You must have `doxygen`, `latex`, and `graphviz` installed for
this to work.

To create the Doxygen documentation run the commands
```
$ cd doxygen 
$ doxygen Doxyfile
```
(i.e., run `doxygen` from this directory).
The documentation will be stored in the `html` directory that is created.
Open the `html/index.html` file to view the documentation.


## Requirements

### Doxygen

Doxygen is a source code documentation system that is widely used in open
source projects.
You will need a copy of the `doxygen` executable on your system to create the
documentation.

It is available [on its website](https://www.doxygen.nl), along with its
[manual](https://www.doxygen.nl/manual/index.html).
The full description of its commands is under the
[Special Commands](https://www.doxygen.nl/manual/commands.html) section.

### Graphviz

Doxygen must be able to find the `graphviz`.
This is used to create equations from embedded LaTex.
You may download `graphviz` from [its website](http://www.graphviz.org).

If you receive the error message
```
> sh: dot: command not found
> Problems running dot. Check your installation!
```
this means that `doxygen` cannot find `graphviz`.


## Files

This directory includes the setup file `Doxyfile` that is needed to run
`doxygen`.
It was created with the gui that is lauchned with the command 
```
doxywizard Doxyfile
```
this is used to set the parameters for the generated output.

The document file `doxygen_setup_style.txt` in this directory contains some
notes on how Cloudy uses `doxygen`.



## Questions

If you run into problems, visit the
[Cloudy forum](https://cloudyastrophysics.groups.io), where you can search
among previous questions for help, or post a new question.

Good luck,<br>
Gary Ferland<br>
[https://www.nublado.org](https://www.nublado.org)
