#!/usr/bin/perl
#################################################################################
#  show-header-tree.pl::
#     Display the tree of header calls for the entire project or the header files
#  provided on the command line.
#
#  Chatzikos, Sep 23, 2013
#################################################################################

use strict;
use warnings;


use lib "./";
use headers;

&headers::header_dependencies();
&headers::prtHdrBranches( @ARGV );
