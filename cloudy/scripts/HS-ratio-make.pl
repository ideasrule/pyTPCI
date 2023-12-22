#!/usr/bin/perl

use strict;
use warnings;

print "Enter input file name:\t";
chomp( my $file = <STDIN> );

open FILE, "< $file"
  or die "Error: Could not open:\t $file\n";
my @contents = <FILE>;
close FILE
   or warn "Warning: Could not close:\t $file\n";

print "Enter line label:\t";
chomp( my $emline = <STDIN> );

$file = "HSrat.dat";
open FILE, "> $file"
  or die "Error: Could not open:\t $file\n";
foreach my $oldline ( @contents )
{
	next	if( $oldline !~ m/^Ca B/i );
	my $newline = $oldline;
	$newline =~ s/^Ca B/$emline/i;
	print FILE $newline;
	print FILE $oldline;
	print FILE "#\n";
}
close FILE
   or warn "Warning: Could not close:\t $file\n";
print "Done.  Created:\t". $file ."\n";
