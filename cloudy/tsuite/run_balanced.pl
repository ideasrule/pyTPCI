#!/usr/bin/perl

# This is a little perl script that runs both the auto and slow test suite in
# parallel on suitable machines or clusters. See the file auto/run_parallel.pl
# for a full description of how this script should be used.
#
# Peter van Hoof
#==============================================================================
# This script is an adaptation of Peter's that keeps job load even across
# processors, based on timings from previous runs.  Timings made by perl are
# to be preferred over Cloudy's, due to complications in non-sequential grid
# calculations.
#
# The script accepts a maximum exec time to discreminate against slow sims.
# This is given as the third argument.  Additional arguments are interpreted
# as sims, and are the only ones executed.
#
# Chatzikos, Dec 17, 2013	First Version
# 	     Jan 21, 2014	Maximum time cutoff

use strict;
use warnings;

use RunServiceBal;
use Timings;

# order matters for reproducing Peter's default order
my @dirs = qw/ slow auto /;

sub check_suite
{
	my $dir = shift;
	print "\n\n\n";
	print "========================================\n";
	print "      Checking $dir test suite          \n";
	print "========================================\n";
	system "cd $dir; ./checkall.pl";
}


# this sets the path to various executables, sets up the list of CPU slots,
# and determines the number of CPUs to be used
my ($nproc, $maxTime, @usims) = &RunServiceBal::initialize( "." );

die "Not-a-number: maxTime=\"$maxTime\"\n" if( not &Timings::isNumber( $maxTime ) );
die "Negative number: maxTime= $maxTime\n" if( $maxTime < 0. );


my $jobList = &Timings::getOrder( $nproc, $maxTime, \@dirs, @usims );

die "Empty job list\n"
	if ($jobList eq "");

#	print	"JOBLIST:\t \"$jobList\"\n";


# clean up the test suite
# do after load balancing, in case the *.out files
# are needed
system "./clean_tsuite.pl > /dev/null";

print "\nNow I will run the test suite using $nproc processors\n\n";



# this does the heavy lifting...
my ($tslot, %timings) = &RunServiceBal::run_jobs( $jobList );


# now do the checking
foreach my $dir ( sort @dirs ) { &check_suite( $dir ); }


&Timings::report_slot_timings	($tslot, %timings);
&Timings::record_timings	($nproc, \@dirs, \@usims, %timings);
