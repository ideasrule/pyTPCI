# This file contains routines that are common to all the run_parallel.pl
# scripts. It should not be called directly (there isn't even any to code to
# execute if you tried...). See the file auto/run_parallel.pl for a full
# description of how the script should be used.
#
# Peter van Hoof

use strict;
use warnings;
use Math::BigInt;

package RunService;

# global variables, these are all defined in initialize()
# except for @slots and @pids, they are not changed afterwards
my $base_dir; # path to the directory from which run_parallel.pl was called
my $exe;      # path to the Cloudy executable
my $runexe_r; # path to the script running a single instance of Cloudy, uses -r syntax
my @slots;    # array of host names associated with each processor slot
my @pids;     # array of process IDs associated with each processor slot
my @users;    # array of user IDs to use on each processor slot
my @procmask; # processor mask for taskset

#
# Set up the global variables needed by run_jobs() and subsidiaries.
#

sub initialize {
	chomp( $base_dir = `pwd` );

	# Set $pin = 1 to set task affinity of one process per physical core
	# Would be better to control this with an input flag...
	my $pin = 0;
	if (getpwuid $< eq 'cloudtests') {
		$pin = 1;
	}

	if( $#ARGV == -1 ) {
		$ARGV[0] = ".";
	}
	if( $#ARGV == 0 ) {
		$ARGV[1] = "";
	}

#     find the Cloudy executable
#     don't check for existence of $exe! it may be a complex command...
#     e.g. we should be able to type: run_parallel.pl "valgrind cloudy.exe"
	if( -d "$_[0]/../source/$ARGV[0]" ) {
#         don't use $_[0] here, this path should always be relative to the tsuite itself!
		$exe = "../../source/$ARGV[0]/cloudy.exe";
	}
	else {
		$exe = "$ARGV[0]";
	}

	$runexe_r = "$base_dir/$_[0]/auto/run_single_r.pl";

	my $nproc;
#     test HOSTNAME in case we forget to use the correct parameter on the dlx...
#     using "<#CPUs>" instead of "dlx" still runs on the dlx, but would load all
#     jobs on the first node, while leaving the remaining allocated nodes unused...
#     there actually are two login nodes; like this the test below works on both.
	if( $ARGV[1] eq "dlx" || ( exists $ENV{"HOSTNAME"} && $ENV{"HOSTNAME"} =~ /^dlxlogin2/ ) ) {
		init_slots_dlx();
		$nproc = @slots;
	}
	elsif( $ARGV[1] eq "static" ) {
		init_slots_static();
		$nproc = @slots;
	}
	else {
		if( "$ARGV[1]" eq "" ) {
			$nproc = get_ncpus();
		}
		else {
			$nproc = $ARGV[1];
		}
		init_slots_local( $nproc , $pin );
	}
	return $nproc;
}


#
# This processes each of the Cloudy input scripts provided in jobList, which is the sole parameter to
# this routine. It is a long string containing a space separated list of all the input script names.
#

sub run_jobs {
	my $jobList = shift;
	my $np_used = 0;

	foreach my $input ( split( / +/, $jobList ) ) {
		chomp( $input );
		my $s = wait_for_slot( \$np_used );
		my $command = "";
		if( "$slots[$s]" ne "localhost" ) {
			$command .= "ssh -x ";
			if( "$users[$s]" ne "" ) {
				$command .= "$users[$s]@";
			}
			$command .= "$slots[$s] ";
		}
		elsif ( $procmask[$s] ne 0x0 ) {
			$command .= "taskset " . $procmask[$s]->as_hex() . " ";
		}
		$command .= "$runexe_r $base_dir \"$exe\" $input";
		if( my $pid = fork() ) {
#             this is the parent process: do the bookkeeping and continue to next job
			print "started $input on host $slots[$s] (slot $s, pid $pid)\n";
			$pids[$s] = $pid;
			++$np_used;
		}
		else {
#             this is the child process: execute Cloudy and exit
			system "$command";
			exit 0;
		}
	}

#     no more jobs, so wait for everything to finish...
	drain_slots( \$np_used );
}


#
# This processes each of the Cloudy input scripts provided in jobList in MPI mode.
#

sub run_jobs_mpi {
	my $jobList = shift;
	my $nproc = shift;

	foreach my $input ( split( / +/, $jobList ) ) {
		chomp( $input );
		system "$runexe_r $base_dir \"mpirun -np $nproc $exe\" $input";
	}
}


#
# This routine sets up the list of available CPU slots and the host names associated with them.
# It is intended to be used on a distrbuted cluster like our Dell cluster called dlx.
# Please contact your system administrator to find out how to modify this for your local cluster.
#

sub init_slots_dlx {
#     check if we are running under the batch system
	if( ! exists $ENV{"SLURM_NODELIST"} ) {
		die "environment variable SLURM_NODELIST not set, did you use the batch system?\n";
	}
#     get the list of nodes assigned to us by slurm
#     it will look something like: "cnode335" or "cnode[195-209,227-233,238-245,276,316]"
	my $nodelist = $ENV{"SLURM_NODELIST"};
#     after these manipulations it will be a list of numbers that can be parsed by perl
	$nodelist =~ s/cnode//;
	$nodelist =~ s/-/../g;
	$nodelist =~ s/[\[\]]//g;
#     for the 2nd example, $nodelist would now be "195..209,227..233,238..245,276,316"
#     this will parse $nodelist and store the node numbers one by one in the array @nodes
#     the @nodes array will hold the NUMBERS of the nodes assigned to us, not the full names
	my @nodes = ( eval "$nodelist" );
#     now set up the array of slots, each slot corresponds to a single core
#     for each of these cores we will store the host name in the @slots array
#     the @pids array will hold the process ID of the job running in that slot
	my $n = 0;
	for( my $i=0; $i <= $#nodes; ++$i ) {
		my $nodestr = sprintf "%3.3d", $nodes[$i];
#         there are 16 cores per node on the dlx, so we create 16 slots per host
		for( 1 .. 16 ) {
#             store the name of the host; on the dlx they are "cnode001" .. "cnode256"
			$slots[$n] = "cnode$nodestr";
			$pids[$n] = 0;
			$procmask[$n] = Math::BigInt->new(0x0);
			$users[$n++] = "";
		}
	}
}


#
# This routine sets up the CPU slots using a file that specifies the machines to use
#

sub init_slots_static {
#     name of the file specifying the static cluster
	my $hostList = $ENV{"HOME"} . "/.cloudy_hosts.txt";
	my $n = 0;
	if( ! open( FOO, "<$hostList" ) ) {
		die "host file $hostList not found!\n";
	}
	while( <FOO> ) {
		if( ! ( /^#/ || /^\s*$/ ) ) {
			chomp;
			s/^\s*//;
			my @fields = split( /\s+/ );
			my $host = shift @fields;
			my $ncore = 1;
			my $user = "";
			my $next = shift @fields;
			while( $next ne "" ) {
				if( $next =~ /cpu=/ ) {
					$next =~ s/cpu=//;
					$ncore = $next;
				}
				if( $next =~ /user=/ ) {
					$next =~ s/user=//;
					$user = $next;
				}
				$next = shift @fields;
			}
#             Now set up the array of slots, each slot corresponds to a single core.
#             For each of these cores we will store the host name in the @slots array.
#             The @pids array will hold the process ID of the job running in that slot.
			for( 1 .. $ncore ) {
				$slots[$n] = "$host";
				$pids[$n] = 0;
				$procmask[$n] = Math::BigInt->new(0x0);
				$users[$n++] = "$user";
			}
		}
	}
	close( FOO );
}

#
# Subroutine to parse Linux /proc/cpuinfo into an array
#

sub parse_cpuinfo {
	my @info;
	my $i = {};
	open(INFO,'/proc/cpuinfo') or die ("could not open /proc/cpuinfo");
	while(<INFO>) {
		chomp;
		if (/^\s*$/)
		{
			push @info, $i;
			$i = {};
		} else {
			my ($k, $v) = split(/\s+:\s+/);
			$i->{$k} = $v;
		}
	}
	if (scalar keys %$i) {
		push @info, $i;
	}
	close(INFO);
	return @info;
}

#
# Subroutine to create per-core taskset masks from cpu_info data
#
#   $procmask[$i] = 0xbits;
# with least significant bit corresponding to first core.
# Typically with hyperthreading, first half of logical processor
# numbers cycle through the physical cores, next half cover the
# additional hyperthreads in same order.  So
#   $procmask[$i] = Math::BigInt->new(0x1) << $i;
# with $nproc = number of physical processors would populate first
# thread on each core
#   $bigone = Math::BigInt->new(0x1);
#   $procmask[$i] = ($bigone << $i) | ($bigone << ($nproc+$i));
# would allow tasks to swap between threads on each core

sub create_masks {
	my @info = @_;
	my %threads;
	foreach my $i (@info) {
		my $proc = "$i->{'physical id'}:$i->{'core id'}";
		if (!exists $threads{$proc}) {
			$threads{$proc} = ();
		}
		push @{$threads{$proc}}, $i;
	}
	my @masks;
	my $bigone = Math::BigInt->new(0x1);
	for my $proc (sort keys %threads)
	{
		my $mask = Math::BigInt->new(0x0);
		for my $t (@{$threads{$proc}}) {
			$mask |= $bigone << $t->{'processor'};
		}
		push @masks, $mask;
	}
	return @masks;
}

#
# This routine sets up the CPU slots for running on a local multi-core machine
#

sub init_slots_local {
	my $nproc = shift;
	my $pin = shift;
	my (@info, @masks);
	if ($pin)
	{
		@info = parse_cpuinfo();
		@masks = create_masks(@info);
		foreach my $m (@masks) {
			print $m->as_hex(),"\n";
		}
	}
	for( my $i=0; $i < $nproc; ++$i ) {
		$slots[$i] = "localhost";
		$pids[$i] = 0;
		if ($pin and $nproc <= $#masks+1) {
			$procmask[$i] = $masks[$i];
		} else {
			$procmask[$i] = Math::BigInt->new(0x0);
		}
		$users[$i] = "";
	}
}


#
# wait until a slot is available and return the slot number to the caller
#

sub wait_for_slot {
	my $np_used = shift;
#     first check if we already have an empty slot
#     this only happens in the initial phase...
	for( my $i=0; $i <= $#pids; ++$i ) {
		if( $pids[$i] == 0 ) {
			return $i;
		}
	}
#     none found, we have to wait for any job to finish
	my $pid = wait();
#     now find the slot occupied by this job and free it
	for( my $i=0; $i <= $#pids; ++$i ) {
		if( $pids[$i] == $pid ) {
			$pids[$i] = 0;
			--$$np_used;
			return $i;
		}
	}
}


#
# wait until all jobs have finished
#

sub drain_slots {
	my $np_used = shift;
	while( $$np_used > 0 ) {
		my $pid = wait();
		for( my $i=0; $i <= $#pids; ++$i ) {
			if( $pids[$i] == $pid ) {
				$pids[$i] = 0;
				--$$np_used;
			}
		}
	}
}


#
# determine the number of CPUs on this host (doesn't work on all systems)
#

sub get_ncpus {
	my $ncpus = 1;
	if( exists $ENV{'NUMBER_OF_PROCESSORS'} )
	{
#         the following should work under Windows.
		$ncpus = $ENV{'NUMBER_OF_PROCESSORS'};
	}
	else
	{
#         this branch is for non-Windows
		if( -r "/proc/cpuinfo" )
		{
#             this branch is for Linux
			$ncpus = 0;
			open( FOO, "</proc/cpuinfo" );
			while( <FOO> )
			{
				if( /^processor/ )
				{
					++$ncpus;
				}
			}
			close( FOO );
		}
		else
		{
#	          this branch is for Mac and BSD variants
			my $sysctl;
			if( -x "/usr/sbin/sysctl" )
			{
#		          this is for Mac and older BSD variants
				$sysctl = "/usr/sbin/sysctl";
			}
			elsif( -x "/sbin/sysctl" )
			{
#		          this is for newer BSD variants
				$sysctl = "/sbin/sysctl";
			}
			if( $sysctl ne "" )
			{
				chomp( my $res = `$sysctl hw.ncpu` );
				my @words = split( / +/, $res );
				$ncpus = $words[1];
			}
		}
	}
	return $ncpus;
}

sub create_list {
	my (@argv) = @_;
	my $jobList = "";
	for( my $i=2; $i <= $#argv; ++$i ) {
		$jobList .= " $argv[$i] ";
	}
#     remove leading and trailing blank...
	$jobList =~ s/^ //;
	$jobList =~ s/ $//;
	return $jobList;
}

sub shuffle_list {
	my $jobList1 = shift;
	my $jobList2 = shift;
#     this reshuffles $jobList2 into the order of $jobList1...
#     in an attempt to load balance this job list.
	my $jobList = "";
	while( defined( my $input = glob($jobList1) ) ) {
		if( $input !~ /^grid00000/ && $jobList2 =~ /\b\Q$input\E\b/ ) {
			$jobList .= " $input ";
		}
	}
#     remove leading and trailing blank...
	$jobList =~ s/^ //;
	$jobList =~ s/ $//;
	return $jobList;
}

1;
