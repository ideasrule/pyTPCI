#!/usr/bin/perl

# This is a perl script that reruns failed sims in the test suite in parallel
# on suitable machines or clusters. This script assumes that you initially ran
# the test suite with run_parallel.pl and discovered problems that way. After
# fixing the code, you can rerun the failed sims with this script. It looks
# for failed scripts in the file serious.txt and reruns those in the same
# sequence as run_parallel.pl would have done. This implies that any sim with
# botched asserts or crashes will be run again.
#
# The calling sequence is the same as for run_parallel.pl. See that file for
# further details.
#
# Peter van Hoof

use strict;
use warnings;

do "../run_service.pm";

# this sets the path to various executables, sets up the list of CPU slots,
# and determines the number of CPUs to be used
my $nproc = RunService::initialize( ".." );

# create list of all input scripts in correct order for load balancing
my $jobList1;
chomp( $jobList1 = `./run_sequence.pl` );
# create list of input scripts with problems
my $jobList2;
if( $#ARGV >= 2 ) {
	$jobList2 = RunService::create_list( @ARGV );
}
else {
	chomp( $jobList2 = `./list_broken.pl` );
}

my $jobList = RunService::shuffle_list( $jobList1, $jobList2 );

if( $jobList ne "" ) {
	print "\nNow I will run the test suite using $nproc processors\n\n";

#     this does the heavy lifting...
	RunService::run_jobs( $jobList );

	if( $nproc > 1 && $jobList =~ /\bfunc_trans_read.*\.in\b/ ) {
#         just in case func_trans_read was run before func_trans_punch was complete
		my @jobs=();
		foreach my $job ( split ( / +/, $jobList ) ) {
			push @jobs, $job if ( $job =~ /\bfunc_trans_read.*\.in\b/ );
		}
		RunService::run_jobs( join (' ', @jobs ) );
	}

#     now do the checking
	system "./checkall.pl ; ./wl_checkall.pl";
}
