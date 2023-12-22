## Notes

```
# The test suites use the following TABLE SED commands to import stellar SEDs
# without requiring the user to install the corresponding grid:
#
# table SED "star_mihalas_46000.dat" 
# table SED "star_kurucz_35000.dat" 
# table SED "star_kurucz_39600.dat" 
#
# These files were originally created using the SAVE TRANSMITTED CONTINUUM
# command and then converted into SED tables to avoid the need to recompile them
# every time the Cloudy frequency mesh changes. Below is the script that was
# used to create the original transmitted continuum output.
#
# uncomment each of the following pairs to create transmitted radiation field
# needed for the test suite to function properly
# 
# table star kurucz 39600
# save transmitted continuum "star_kurucz_39600.dat"
# 
# table star mihalas 46000
# save transmitted continuum "star_mihalas_46000.dat"
# 
table star kurucz 35000
save transmitted continuum "star_kurucz_35000.dat"
# 
set dr -10
hden -5
init "honly.ini"
ionization parameter -2
stop zone 1
constant temperature 4
```
