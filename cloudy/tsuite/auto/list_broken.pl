#!/usr/bin/perl

# This is a little helper script for rerun_parallel.pl that outputs all the
# sims that had problems according to what is reported in the file
# serious.txt. This will include any type of problem like botched monitors,
# aborts, failures, and crashes.
#
# Peter van Hoof

use strict;
use warnings;

my $list = "";

open( my $Foo, "<serious.txt" );
while( <$Foo> ) {
	my @word = split( / +/ );
	my $input = $word[0];
	if( $input =~ /\.out:$/ ) {
		$input =~ s/\.out:/.in/;
#         prevent input scripts from being added twice
		if( $list !~ / $input / ) {
			$list .= " $input ";
		}
	}
}
close( $Foo );

# remove leading and trailing blank...
$list =~ s/^ //;
$list =~ s/ $//;

print "$list\n";
