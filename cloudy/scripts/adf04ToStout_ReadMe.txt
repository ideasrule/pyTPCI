Adf04ToStout
============

This script converts ADF04 files from Open ADAS to Stout format.

Syntax:
    adf042stout.py <name of adf04file>
    
Original version 2013 by Matt Lykins

The ADF04 file is designed to be read with Fortran fixed format.
This presents a few challenges.

The file may, or may not, have a last number with no temperature.
If it is present there will be one more "collision strength" than
temperature.  This last number may be negative, in which case
it will run into the previous number.

The collision strengths themselves are in exponential format but
with the "e" suppressed.  So 1.2e-2 is written as 1.2-2.
