#!/usr/bin/perl
#################################################################################
#  show-header-alldep.pl::
#     Display a list of all the headers invoked in the call-tree of a header.  If
#  no header files are requested on the command line, the script operates on all
#  the headers present in the current directory.
#  Chatzikos, Sep 23, 2013
#################################################################################
use strict;
use warnings;

use lib "./";
use headers;

&headers::header_dependencies();
&headers::prtHdrAlldep( @ARGV );
