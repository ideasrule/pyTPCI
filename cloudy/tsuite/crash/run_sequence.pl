#!/usr/bin/perl

# This is a little helper script for run_parallel.pl that outputs the order in which the
# input scripts should be run. The idea is to put the CPU intensive jobs upfront and the
# quick ones last so that you get good load balancing on the CPUs. The order is set up
# in a rather heuristic way, and needs to be updated manually once in a while to reflect
# changes in the way Cloudy runs and/or new test suite scripts. It has no parameters and
# prints the list of jobs to stdout, separated by spaces.
#
# Peter van Hoof

$list = "";

# this funny sequence is to assure that the heavy jobs are upfront...
while( defined( $input = glob("limit*.in *.in") ) ) {
#     prevent input scripts from being added twice
    if( ! ( $list =~ / $input / ) && ! ( $skip =~ / $input / ) ) {
	$list .= " $input ";
    }
}

# remove leading and trailing blank...
$list =~ s/^ //;
$list =~ s/ $//;

print "$list\n";
