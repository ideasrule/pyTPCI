#!/usr/bin/perl
#################################################################################
# show-header-req.pl::
#    Find the header dependencies for a minimal compilation of each header file.
# 'cddefines.h' are used by default, but not reported.
# Chatzikos, Sep 18, 2013
#################################################################################

use strict;
use warnings;

use lib "./";
use headers;


my $ncpus = 1;
if( @ARGV )
{
	$ncpus = int($ARGV[0]);
	die "Usage:\t $0 <ncpus>\n"
	 if( $ncpus <= 0 );
}

print	&headers::header_requirements( $ncpus );
