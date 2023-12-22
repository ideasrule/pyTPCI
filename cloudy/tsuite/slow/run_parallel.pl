#!/usr/bin/perl

# This is a little perl script that runs the slow test suite in parallel on
# suitable machines or clusters. See the file auto/run_parallel.pl for a full
# description of how this script should be used.
#
# Peter van Hoof

use strict;
use warnings;

do "../run_service.pm";

# this sets the path to various executables, sets up the list of CPU slots,
# and determines the number of CPUs to be used
my $nproc = RunService::initialize( ".." );

# clean up the test suite
system "./clean_tsuite.pl";

# create list of all input scripts in correct order for load balancing
my $jobList1;
chomp( $jobList1 = `./run_sequence.pl` );
# create list of input scripts to be executed
my $jobList2;
if( $#ARGV >= 2 ) {
	$jobList2 = RunService::create_list( @ARGV );
}
else {
	$jobList2 = $jobList1;
}

my $jobList = RunService::shuffle_list( $jobList1, $jobList2 );

if( $jobList ne "" ) {
	print "\nNow I will run the test suite using $nproc processors\n\n";

#     this does the heavy lifting...
	RunService::run_jobs( $jobList );

#     now do the checking
	system "./checkall.pl ; ./wl_checkall.pl";
}
