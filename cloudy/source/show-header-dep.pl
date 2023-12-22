#!/usr/bin/perl
#################################################################################
#  show-header-dep.pl::
#     List all the headers included directly by the header file in question.  All
#  header files are assessed if none issued on the command line.
#  Chatzikos, Sep 23, 2013
#################################################################################

use strict;
use warnings;


use lib "./";
use headers;


&headers::header_dependencies();
&headers::prtHdrDep( @ARGV );
