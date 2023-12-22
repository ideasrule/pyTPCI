#!/usr/bin/perl

# This is a little helper script for run_parallel.pl that outputs the order in which the
# input scripts should be run. The idea is to put the CPU intensive jobs upfront and the
# quick ones last so that you get good load balancing on the CPUs. The order is set up
# in a rather heuristic way, and needs to be updated manually once in a while to reflect
# changes in the way Cloudy runs and/or new test suite scripts. It has no parameters and
# prints the list of jobs to stdout, separated by spaces.
#
# If present, it reads the contents of the file skip.dat. This should contain a list of
# input file names that you want to be skipped. There may be multiple script names on one
# line, separated by spaces, or the names may be on multiple lines. Wildcards are not
# allowed. It is OK if this file is absent, in which case all jobs will be executed.
#
# Peter van Hoof

$skip = "";
if( -r "skip.dat" ) {
	open Foo, "<skip.dat";
	while( <Foo> ) {
		chomp;
		$skip .= " $_ ";
	}
	close Foo;
}

$list = "";

# this funny sequence is to assure that the heavy jobs are upfront...
while( defined( $input = glob("time*.in orion*.in h2*orion*.in *xdr*.in m17*.in *extreme*.in " .
							  "d*.in feii*Z20*.in feii*.in *HTT*.in *leiden_v*.in *.in") ) ) {
#     prevent input scripts from being added twice
	if( ! ( $list =~ / $input / ) && ! ( $skip =~ / $input / ) ) {
		$list .= " $input ";
	}
}

# remove leading and trailing blank...
$list =~ s/^ //;
$list =~ s/ $//;

print "$list\n";
