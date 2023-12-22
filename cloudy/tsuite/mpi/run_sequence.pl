#!/usr/bin/perl

# This is a little helper script for run_all.pl that outputs the order in
# which the input scripts should be run. It has no parameters and prints the
# list of jobs to stdout, separated by spaces.
#
# Peter van Hoof

$list = "";

while( defined( $input = glob("*.in") ) ) {
    $list .= " $input ";
}

# remove leading and trailing blank...
$list =~ s/^ //;
$list =~ s/ $//;

print "$list\n";
