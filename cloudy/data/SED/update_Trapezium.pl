#!/usr/bin/perl
#
# Convert the transmitted continuum specified below into the SED format, so that
# it can be used with the 'table' and 'table SED' commands.  Note that the
# intensity of the original continuum is in nuFnu units, and that must be
# specified in the first line of data in the SED file.  Note also that the
# continuum contains bins of zero intensity, which must be removed, or the table
# command will fail.
#
# Chatzikos, Marios	2015-Oct-29
#
use strict;
use warnings;

#
# Read Trapezium Transmitted continuum
#
my $TrapeziumTran = "Trapezium/Trapezium_WMbasic.tran";
open FILE, "< $TrapeziumTran"
  or die "Could not open: $TrapeziumTran\n";
my @SED = <FILE>;
close FILE;

#
# Get rid of header, marked by the second line that matches /^#$/
#
my $ncomment_lines = 0;
while( $ncomment_lines != 2 )
{
	++$ncomment_lines
	  if( $SED[0] =~ m/^#$/ );
	shift( @SED );
}

#
# Write output file
#
my $TrapeziumSED = "Trapezium.sed";
open FILE, "> $TrapeziumSED"
  or die "Could not open: $TrapeziumSED\n";

print FILE "# energy\t nuFnu\n";
my $nlines = 0;
foreach my $line ( @SED )
{
	++$nlines;
	my( $energy, $nuFnu ) = split( /\s+/, $line );
	my $isZero = 0;
	if( $nuFnu == 0. )
	{
		$nuFnu = 1e-30;
		$isZero = 1;
	}
	printf FILE "%.5e\t%.3e", $energy, $nuFnu;
	print FILE "\t nuFnu"	if( $nlines == 1 );
	print FILE "\n";
	last	if( $isZero );
}

close FILE;

print "Written file: $TrapeziumSED\n";
print "Dismissed ". (@SED - $nlines) ." lines with zero intensity\n";
