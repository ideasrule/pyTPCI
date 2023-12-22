#!/usr/bin/perl

# This is a perl script that reruns failed sims in the auto and slow test
# suite in parallel on suitable machines or clusters. See the file
# auto/rerun_parallel.pl for a full description of how this script should be
# used.
#
# Peter van Hoof

use strict;
use warnings;

do "./run_service.pm";

# this sets the path to various executables, sets up the list of CPU slots,
# and determines the number of CPUs to be used
my $nproc = RunService::initialize( "." );

# create list of all input scripts in correct order for load balancing
chomp( my $jobListAuto = `cd auto ; ./run_sequence.pl` );
chomp( my $jobListSlow = `cd slow ; ./run_sequence.pl` );

my $jobList1 = "";
foreach my $input ( split( / +/, $jobListSlow ) ) {
	$jobList1 .= "slow/$input ";
}
foreach my $input ( split( / +/, $jobListAuto ) ) {
	$jobList1 .= "auto/$input ";
}
# create list of input scripts to be executed
my $jobList2 = "";
if( $#ARGV >= 2 ) {
	$jobList2 = RunService::create_list( @ARGV );
}
else {
#     create list of failed scripts
	chomp( my $jobFailAuto = `cd auto ; ./list_broken.pl` );
	chomp( my $jobFailSlow = `cd slow ; ./list_broken.pl` );

	foreach my $input ( split( / +/, $jobFailSlow ) ) {
		$jobList2 .= "slow/$input ";
	}
	foreach my $input ( split( / +/, $jobFailAuto ) ) {
		$jobList2 .= "auto/$input ";
	}
}

my $jobList = RunService::shuffle_list( $jobList1, $jobList2 );

if( $jobList ne "" ) {
	print "\nNow I will run the test suite using $nproc processors\n\n";

#     this does the heavy lifting...
	RunService::run_jobs( $jobList );

	if( $nproc > 1 && $jobList =~ m$\bauto/func_trans_read.*\.in\b$ ) {
#         just in case func_trans_read was run before func_trans_punch was complete
		my @jobs=();
		foreach my $job ( split ( / +/, $jobList ) ) {
			push @jobs, $job if ( $job =~ /\bfunc_trans_read.*\.in\b/ );
		}
		RunService::run_jobs( join (' ', @jobs ) );
	}

#     now do the checking
	print "\n\n\n";
	print "========================================\n";
	print "      Checking auto test suite          \n";
	print "========================================\n";
	system "cd auto ; ./checkall.pl ; ./wl_checkall.pl";
	print "\n\n\n";
	print "========================================\n";
	print "      Checking slow test suite          \n";
	print "========================================\n";
	system "cd slow ; ./checkall.pl ; ./wl_checkall.pl";
}
