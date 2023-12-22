#!/usr/bin/perl

use strict;
use warnings;

if( $#ARGV != 2 )
{
    die "usage $0 <filename> <prefix> <delimiter>\n";
}

my $fnam = shift;
my $prefix = shift;
my $delim = shift;

open( Inp, "<$fnam" ) || die "file $fnam not found!\n";
my $onam = $prefix . "00";
open( Out, ">$onam" );
my $nfile = 1;
while( <Inp> )
{
    print Out "$_";
    if( /$delim/ )
    {
	close( Out );
	++$onam;
	open( Out, ">$onam" );
	++$nfile;
    }
}
close( Inp );
close( Out );
print "$nfile\n";
