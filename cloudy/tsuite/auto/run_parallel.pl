#!/usr/bin/perl

# This is a perl script that runs the auto test suite in parallel on suitable
# machines or clusters. There are actually 3 incarnations of this script: one
# for running the auto test suite (tsuite/auto/run_parallel.pl, this script),
# one for running the slow test suite (tsuite/slow/run_parallel.pl), and one
# for running both test suites combined (tsuite/run_parallel.pl). Some
# examples of the calling sequence for all of these scripts are:
#
# run_parallel.pl
# run_parallel.pl sys_icc
# run_parallel.pl . 8
# run_parallel.pl /path/to/cloudy.exe
# run_parallel.pl "valgrind --leak-check=full cloudy.exe" 8
# run_parallel.pl sys_icc dlx
# run_parallel.pl sys_sunstudio static
#
# where the first parameter is the path to the Cloudy executable, or a more
# complex command between double quotes. If the parameter is completely
# omitted (or a "." is entered), it will default to "../../source/cloudy.exe".
# You can also supply the name of one of the sys_xxx directories, in which
# case the path will default to "../../source/sys_xxx/cloudy.exe". In all
# other cases cloudy.exe should either be in your search path, or the full
# pathname should be used. The second number is the number of processors to
# use. It can be any number > 0. This parameter is optional and will default
# to the number of processors on your computer (this number is only correctly
# determined under Linux and Windows, you will need to set it manually on
# other systems, most notably Mac).
#
# Optionally you can supply additional parameters after the first two that
# were shown above. If present, these will be interpreted as a list of input
# scripts to be executed. The default is to execute all scripts that are
# present. This example will execute only the PDR scripts:
#
# run_parallel.pl . 8 pdr*.in
#
# There is also limited support for doing parallel runs on distributed memory
# clusters. You can do this by supplying "bcx", "dlx" or "static" as the second
# parameter. The option "static" is appropriate if you know beforehand which
# machines you will be using (i.e. when there is no batch system). You need to
# supply a file ".cloudy_hosts.txt" in your home directory that contains the
# following:
#
#   <machine> [cpu=<cpucount>] [user=<userid>]
#   <machine> [cpu=<cpucount>] [user=<userid>]
#    ...
#
# <machine> should be the host name (ssh should be able to DNS resolve this
# name), and cpu=<cpucount> is the number threads you want to run on that
# machine (optional, the default value is 1). If you have a different user ID
# on that machine, you can supply that as well with user=<userid>. You should
# be able to do a password-less login with ssh to these machines for this
# script to work. Supply only one machine per line. You can embed comments
# with a '#' in the first column. This format is deliberately chosen to be
# identical to the one used by MPI (at least LAM MPI).
#
# If there is a batch system on your cluster, you will generally not know
# beforehand which machines will be reserved for you. You may be able to get
# the necessary information by other means, such as environment variables. If
# you supply "bcx" as the second parameter that is exactly what is being done.
# The routine init_slots_bcx() in tsuite/run_service.pm is then called, which
# is specifically set up for our IBM xSeries cluster called bcx in Lexington.
# It runs Moab+Slurm as a scheduling system. If you have a similar cluster, it
# may work for you as well, but most likely you will have to adapt the
# routine. Please contact your local system administrator or helpdesk if you
# need assistence. Every cluster is different so we cannot help you with
# that...
#
# In general using more CPUs will reduce your waiting time (although you may
# spend more time in the input queue if you have a batch system). However, if
# you use more than N CPUs there will be no further reduction in waiting time
# because it will be set by the longest running single job. The value of N
# depends on the test suite you are running. For the auto test suite N = 12,
# for the slow test suite N = 8, and for the combined test suites N = 8-12.
# These numbers may vary though, depending on the platform you are using and
# changes in the code and/or test suites. They appear to be reasonably stable
# though...
#
# NB NB - The script will first delete any remaining output from a previous
# run of the test suite by calling the script clean_tsuite.pl.
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
