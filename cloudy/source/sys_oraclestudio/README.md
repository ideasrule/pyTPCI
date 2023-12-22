## Notes

This directory contains the files needed to generate
an optimized executable with the Sun Studio compiler.

The -Yl,/usr/bin flag is needed to work around an incompatibility between the
(modified) linker supplied by Sun and the object file format in recent Linux
distributions. It is only needed on 64-bit systems. See the bug report
[here](http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6533627)
for further details.

------------------------

To build enter
```
make
``` 
at the command prompt.
