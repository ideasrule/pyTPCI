#!/usr/bin/perl -w
#
# this script looks for problems in the Cloudy test suite.
#
# execute this script from within the directory where the *.out files live.
# The command is: ./checkall.pl
#
# see the counterpart of this script in the auto directory for a more detailed
# description of how it works.
#
# Peter van Hoof

system "../auto/checkall.pl";
