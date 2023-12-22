# Process the script timings obtained from Perl, or directly from the Cloudy
# output, and save them in files in the current directory, or alternatively,
# in auto and slow.
#
# Chatzikos, 2013-Dec-18	First Version
# 	     2014-Jan-21	Use Peter's function to sort a list of
# 	     			command-line sims to an 'optimal' order

use strict;
use warnings;

use RunServiceBal;



package   Timings;
require   Exporter;

our @ISA       = qw/Exporter/;
our @EXPORT    = qw/ getOrder report_slot_timings record_timings /;
our @EXPORT_OK = qw/ get_min_max get_avg_std isNumber /;
our $VERSION   = 2013.12.18;


my $timings_fname_base = "timings";
my @dirs = qw/ slow auto /;



#
#
# Note that use of bare blocks below ensures that 'my' global variables
# of different procedures (see TAGS above) cannot interfere with each other.
#
#





#################################################################################
#                              IO GENERAL                                      #
#################################################################################
sub get_contents
{
	my ($file) = @_;
	open FILE, "< $file"	or die "Error: Could not open:\t $file\n";
	my @contents = <FILE>;
	close FILE 		or warn "Warning: Could not close:\t $file\n";
	return	@contents;
}


sub get_fields
{
	my $line = shift;
	chomp( $line );
	my @fields = split(/\s+/, $line);
	shift( @fields )	while( @fields and $fields[0] eq "" );
	return	@fields;
}


sub isNumber
{
	my ($input) = @_;
	my $isno = 0;

	return $isno	if (not defined($input) or $input =~ m/^\s*$/);

	# Accepts any number: integer, real, positive, negative, and even
	# numbers in scientific/engineering notation.
	# First pattern has an optional real part, while the second has an optional integer part
	# Combining this into one statement would lead to the patterns "+" or "-", matching as numeric
	#
	#
	if( $input =~ /^([\+-]?\d+((\.\d+)?([D-Ed-e][\+-]?\d+)?)?|[\+-]?(\d+)?(\.\d+)?([D-Ed-e][\+-]?\d+)?)$/ )
	{
		$isno = 1;
	}

	return $isno;
}





#################################################################################
#				FILE & DIR NAMES				#
#################################################################################
sub get_dir_fname
{
	my ($dir, $prefix) = @_;
	return	$dir ."/". $prefix ."_". $timings_fname_base .".txt";
}





#################################################################################
#			READ TIMINGS & OPTIMAL SEQUENCE				#
#################################################################################
sub sort_timings_alpha
{
	my %timings = @_;
	return	sort (keys %timings);
}


sub sort_timings_num
{
	my %timings = @_;
	#	foreach my $sim ( keys %timings )
	#	{
	#		print "$sim\t $timings{$sim}\n";
	#	}
	return	sort { $timings{$b} <=> $timings{$a} } ( keys %timings );
}


sub read_timings_from_cloudy_output
{
	my ($maxTime, $nproc, $dir) = @_;

	my @outfiles = glob ( "$dir/*.out" );
	if( @outfiles == 0 )
	{
		#	warn "\t*** Warning: No *.out files in:\t $dir\n";
		return;
	}

	my %timings;
	my $untimed = 0;

	foreach my $file (@outfiles)
	{
		my $ifile = $file;
		$ifile =~ s/\.out/\.in/;

		my @contents = &get_contents( $file );

		my $time = 0;
		my @lines = grep { m/Cloudy ends:/ } @contents;
		foreach my $line ( @lines )
		{
			chomp($line);
			my @els = split (/ExecTime\(s\)/, $line);
			last if( @els == 0 );
			$time += $els[1];
		}

		$time /= $nproc
			if( @lines > 1 );

		$timings{$ifile} = $time
			if( $maxTime == 0 or ($maxTime > 0 and $time <= $maxTime ) );

		if( $timings{$ifile} == 0 )
		{
			#	%timings = ();
			#	last;
			print "*** Warning: No timing in Cloudy output: $file\n";
			$untimed++;
		}
	}

	print "\t $untimed\t files did not have timing info in $dir/\n"
		if ($untimed > 0);

	return	%timings;
}


sub read_timings_report
{
	my ($maxTime, $outfile, $dir) = @_;

	my @contents = &get_contents( $outfile );
	return	()	if( not @contents );

	my %timings;
	foreach my $line ( @contents )
	{
		# Skip comment and empty lines
		chomp($line);
		next	if( $line =~ m/^#/ or $line =~ m/^\s*$/ ); 

		my ($file, $time) = split (/\s+/, $line);

		$file =~ s/\.out/\.in/;
		$file = "$dir/$file"	if( defined($dir) and $file !~ m#/# );

		$timings{$file} = $time
			if( $maxTime == 0 or ($maxTime > 0 and $time <= $maxTime ) );
	}

	return	%timings;
}


sub read_master_timings
{
	my $maxTime = shift;
	my @masters = reverse sort glob( "$timings_fname_base*" );
	return	if( not @masters );

	# This guarantees that 'timings-r9000.txt' will be
	# picked up if 'timings-r9000-2sims.txt' also exists
	#
	my $master_fname = shift( @masters );
	return	if( $master_fname =~ m/sims.txt$/ );

	return	join( " ", &sort_timings_num( &read_timings_report( $maxTime, $master_fname ) ) );
}


sub read_timings_from_file
{
	my ($maxTime, $dir, $prefix) = @_;

	my $file = &get_dir_fname($dir, $prefix);
	return	()	if (not -e $file);

	return	&read_timings_report( $maxTime, $file, $dir );
}


sub read_timings
{
	my ($maxTime, $nproc, $dir) = @_;

	# First, try to read the perl time estimates.
	#
	my %timings = &read_timings_from_file($maxTime, $dir, "perl");

	# If not, try to get the timings directly
	# from the Cloudy output files.
	#
	%timings = &read_timings_from_cloudy_output($maxTime, $nproc, $dir)
		if( not %timings );

	#	my @in = glob("$dir/*.in");
	#	if( scalar(keys %timings) != scalar(@in) )
	#	{
	#		warn "Some scripts have no timing info.  Setting times to 0\n";
	#		foreach my $script ( @in )
	#		{
	#			$timings{$script} = 0
	#				if( not exists $timings{$script} );
	#		}
	#	}

	return %timings;
}


sub getOrder
{
	my ($nproc, $maxTime, $dirs, @usims) = @_;

	# First, look for master file in present dir
	#
	my $jobOrder = &read_master_timings( $maxTime );

	# If not found, search for output in directories
	#
	if( not defined( $jobOrder ) )
	{
		my %timings;
		foreach my $dir ( @$dirs )
		{
			%timings = ( %timings, &read_timings( $maxTime, $nproc, $dir ) );
		}

		if( %timings )
		{
			$jobOrder = join(" ", &sort_timings_num( %timings ) );
		}
		else
		{
			# If nothing else found, go with Peter's default order
			#
			foreach my $dir ( @$dirs )
			{
				my $list = `cd $dir; ./run_sequence.pl`;
				my @fields = &get_fields( $list );
				for( my $i = 0; $i < @fields; $i++ )
				{
					$fields[$i] = "$dir/". $fields[$i];
				}
				$list = join(" ", @fields);
				$jobOrder .= " $list ";
			}
			substr($jobOrder,  0, 1, "");
			substr($jobOrder, -1, 1, "");
		}
	}

	if( @usims )
	{
		# Use command-line issued sims only, sorted in the 'optimal' order
		my $jobOrder2 = join( " ", @usims );
		$jobOrder = &RunServiceBal::shuffle_list( $jobOrder, $jobOrder2 );
	}
	elsif( $maxTime > 0 )
	{
		print "I will run the ". scalar(my @a = split(/\s/, $jobOrder)) ." sims with t <= $maxTime seconds\n";
	}
	else
	{
		# Make sure we haven't skipped any sims
		foreach my $dir ( @$dirs )
		{
			my @sims = glob( "$dir/*.in" );
			foreach my $sim ( @sims )
			{
				next	if( $sim =~ m/grid\d+/ );
				$jobOrder .= " $sim"
				if( $jobOrder !~ $sim );
			}
		}
	}

	return	$jobOrder;
}





#################################################################################
#			TIMING PROCESSING & OUTPUT				#
#################################################################################
sub get_rev
{
	my @dirs = @_;
	my @out;
	foreach my $dir ( @dirs ) { @out = ( @out, glob("$dir/*.out") ); }
	my $rev = "";
	if( @out )
	{
		my $out = shift( @out );	@out = ();
		my @contents = &get_contents( $out );
		my $line = shift( @contents );	@contents = ();
		my @fields = &get_fields( $line );
		$rev = "-". substr( $fields[2], 0, length($fields[2])-1 );
	}
	return	$rev;
}


sub get_master_fname
{
	my ($ncpus, @usims) = @_;
	my $nsims = "";
	$nsims = "-". @usims ."sims"	if( @usims );
	return	$timings_fname_base . &get_rev( @dirs ) ."-n$ncpus". $nsims .".txt";
}


sub write_file
{
	my ($fname, $timings, $dir) = @_;

	if( not open FILE,     "> $fname" )
	{
		warn	"Could not open timings output file:\t $fname\n";
		return	0;
	}

	printf FILE	"#%-40s\t%s\n", " Sim", "Time [s]";

	foreach my $file ( &sort_timings_alpha( %$timings ) )
	{
		next	if( defined($dir) and not -e "$dir/$file" );
		printf FILE "%-40s\t%8.2f\n", $file, $$timings{$file};
	}

	close FILE;

	print "Done. Created:\t $fname\n";

	return	1;
}


sub write_master_timings
{
	my ($ncpus, $usims, %timings) = @_;
	return	&write_file( &get_master_fname( $ncpus, @$usims ), \%timings );
}


sub write_timings_per_dir
{
	my ($dir, $prefix, %timings) = @_;
	return	&write_file( &get_dir_fname($dir, $prefix), \%timings, $dir );
}


sub record_timings
{
	my( $ncpus, $dirs, $usims, %timings ) = @_;

	if( %timings )
	{
		# Save timings in dir-specific files, if master failed
		if( not &write_master_timings( $ncpus, $usims, %timings ) )
		{
			foreach my $dir ( @$dirs )
			{
				&write_timings_per_dir( $dir, "perl", %timings );
			}
		}
	}

	return;
}


sub report_slot_timings
{
	my( $tslot, %timings ) = @_;

	print "\n" x3 ."Timings per slot:\n";
	for( my $i = 0; $i < @$tslot; $i++ )
	{
		printf "\tCPU: %2d\t %g\n", $i+1, $$tslot[$i];
	}
	my ($avg, $std) = &get_avg_std( @$tslot );
	my ($min, $max) = &get_min_max( @$tslot );
	print	"CPU Timings: avg=$avg, std=$std, Dtmax=". ($max-$min) ."\n";
	print	"Texec: ". $max ." s\t (";
	printf "%d h  ",   (int($max /3600)      );
	printf "%d m  ",   ((int($max /60)  % 60));
	printf "%d s)\n\n",  ((($max % 3600) % 60));

	return;
}





#################################################################################
#				STATISTICS					#
#################################################################################
sub get_avg_std
{
	my (@array) = @_;

	my $nelem = scalar(@array);
	my ($avg, $std) = ( 0, 0 );
	foreach my $val (@array)
	{
		$avg += $val;
		$std += $val * $val;
	}
	$avg /= $nelem;
	$std = sqrt( $std/($nelem) - $avg * $avg );

	return	($avg, $std);
}


sub get_min_max
{
	my (@array) = @_;
	@array = sort { $a <=> $b } @array;
	return	($array[0], $array[-1]);
}

1;
