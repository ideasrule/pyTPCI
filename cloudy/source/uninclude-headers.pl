#!/usr/bin/perl
#################################################################################
#  uninclude-headers.pl::
#
#     For each .cpp file in the current directory, this script extracts a list of
#  the included header files, except for those issued from within #if--#endif
#  blocks.  It then determines whether each of these headers is indeed required
#  for the compilation of the .cpp file by omitting it from the body of the file
#  and looking for compilation errors.
#     The scipt automatically removes headers already invoked by a header file
#  that must be included.  It also rearranges the order of #include statements,
#  if there is a pair of headers one of which invokes a file the other requires.
#  It removes duplicate declarations, and tries to bring most #include statements
#  at the beginning of the file.
#     There are two known limitations at the moment.  First, #include statements
#  from within #if-blocks cannot be moved to the beginning of the file.  More
#  importantly, headers included from within #define--#undef blocks can be
#  erroneously removed.  Case in point, physconst.cpp.
#     Exercise caution when running the script and double check what has been
#  done, even if 'make' completes successfully.  YOU MAY HAVE TO EDIT A FILE OR
#  TWO BY HAND!
#
#  Chatzikos, Sep 23, 2013		First Version
#
#################################################################################

use strict;
use warnings;


use threads;
use threads::shared;


my (@header_output_files, @nheaders_files, @nheaders_head, @nheaders_fail) : shared;
my (@output_files, @nfiles, @nhead, @nfail, @nmiss, @nmult) : shared;


my $compile_flags = '-I ../library/UnitTest++/src';
my $redundant_cpp = "redundant-headers-in-cpp.list";

my @dontAnalyze	  = qw/ cdstd.h vectorize*.h /;


my $cpp_extension = ".cpp";
my $dry_run = 0;

use lib "./";
use headers;





#################################################################################
#				USER INPUT					#
#################################################################################
sub printUsage
{
	die "$0 [-h] [-d] [<ncpus>] [<source>]\n"
	 .	"where:\n"
	 .	"\t      -h:\t print this message and exit\n"
	 .	"\t      -d:\t dry run, produce a summary, but change no files\n"
	 .	"\t <ncpus>:\t number of CPUs to use [OPTIONAL]\n"
	 .	"\t<source>:\t list of source (.cpp) files to process [OPTIONAL]\n"
	 .	"\t         \t must always follow <ncpus>\n";
}

sub remove_tmp_files
{
	my @test_cpp = glob "test-*.cpp";

	foreach my $this_cpp ( @test_cpp )
	{
		if( $this_cpp =~ m/^test-\d+\.cpp$/ )
		{
			print "Removing $this_cpp...\n";
			unlink $this_cpp;

			my $objfile = $this_cpp;
			$objfile =~ s/\.cpp$/.o/;
			unlink $objfile;
		}
	}
}

sub ConfirmInput
{
	my( $ncpus, $cpp ) = @_;

	if( defined( $ncpus ) and $ncpus <= 0 )
	{
		die "Error: Number of CPUs must be positive; got $ncpus\n";
	}

	if( @$cpp )
	{
		foreach my $this_cpp ( @$cpp )
		{
			#
			# This is impossible when called with command-line arguments
			#
			die "Error: Undefined source file\n"
			 if( not defined( $this_cpp ) or $this_cpp eq "" );

			die "Error: Source file ('$this_cpp') does not exist or empty\n"
			 if( not -e $this_cpp or not -s $this_cpp );
		}
	}
}

sub Init
{
	my @argv = @_;

	&printUsage()	if( $argv[0] eq "-h" );
	if( $argv[0] eq "-d" )
	{
		shift( @argv );
		$dry_run = 1;
	}

	my $ncpus = int( shift( @argv ) );
	my @cpp = @argv;


	&ConfirmInput( $ncpus, \@cpp );

	&remove_tmp_files();
	$ncpus = 1	if( not defined( $ncpus ) );
	@cpp = glob "*.cpp"	if( not @cpp );

	if( 0 )
	{
		print "ncpus=\t $ncpus\n";
		print "cpp:\t @cpp\n";
		die;
	}

	return	( $ncpus, @cpp );
}

#################################################################################
#				FILE PROCESSING					#
#################################################################################
sub copy_all_but_include
{
	my ($tmp, $excl, @contents) = @_;

	open TMP, "> $tmp"
	  or die "copy_all_but_include: Error: Could not open:\t $tmp\n";

	my $exclude = join( " ", @$excl );
	foreach my $line ( @contents )
	{
		my $header = &headers::parse_include( $line );
		next	if( $line =~ m/^#\s*include/ and $exclude =~ m/\b\Q$header\E\b/ );
		print TMP $line;
	}

	close TMP
	   or warn "copy_all_but_include: Warning: Could not close:\t $tmp\n";

	return;
}



sub fix_cpp
{
	my ($cpp, $contents, $headers) = @_;

	open CPP, "> $cpp"
	  or die "fix_cpp: Error: Could not open:\t $cpp\n";

	my $isinif = 0;
	foreach my $line ( @$contents)
	{
		if( $line !~ m/^#\s*include/			or
		    $line =~ m/^#\s*if/  			or
		    $line =~ m/^#\s*endif/ 			or
		    $isinif > 0 )
		{
			$isinif++	if( $line =~ m/^#\s*if/ );
			$isinif--	if( $line =~ m/^#\s*endif/ );
			print CPP $line;
		}
		elsif( $isinif == 0 )
		{
			my $head = &headers::parse_include( $line );
			if( not exists $$headers{$head}		or
			    $$headers{$head}{comp} == 0 )
			{
				print	CPP $line;
			}
		}
	}

	close CPP
	   or warn "fix_cpp: Warning: Could not close:\t $cpp\n";

	return;
}



sub cpp_dep
{
	my ($thread, @cpp) = @_;

	my $redundant = "$redundant_cpp-$thread";
	open REDUN, "> $redundant"
	  or die "Error: Could not open:\t $redundant\n";

	my $tmp = "test-$thread$cpp_extension";
	print "thread:\t $thread\t #cpp:\t". @cpp ."\n"		if( 0 );

	my ($nfiles, $nheads, $nfails, $nmuldc, $nmissdc) = (0, 0, 0, 0, 0);
	$nfiles[$thread] = $nhead [$thread] = $nfail [$thread] =
	$nmult [$thread] = $nmiss [$thread] = 0;
	foreach my $cpp ( @cpp )
	{
		print "$cpp\n"	if( 0 );

		my $header_for_this_source = $cpp;
		$header_for_this_source =~ s/\.cpp$/.h/;

		my @noAnalysis = ( $header_for_this_source, @dontAnalyze );

		my @contents = &headers::get_contents( $cpp );
		my @headers  = &headers::get_includes_from_cpp( \@noAnalysis, @contents );

		my( $nmult, $order, $mult_headers, %headers );
		($nmult, $mult_headers, @contents)
				    = &headers::uniq_head_contents( \@headers, @contents );
		      %headers      = &headers::header_calls   ( @headers );
		($order, @contents) = &headers::reorder_headers( $cpp, \%headers, @contents );

		my $nhead = scalar( @$order );
		foreach my $head ( @$order ) { $headers{$head}{comp} = -1; }

		foreach my $head ( @$order )
		{
			print "\tHEADER:\t $head\n"	if( 0 );

			my @exclude = ( $head );
			if( $headers{$head}{called} )
			{
				my @hdep = @{ $headers{$head}{called} };

				foreach my $h ( @hdep )
				{
					if( $headers{$h}{comp} == 0 )
					{
						$headers{$head}{comp} = 1;
						last;
					}
					else
					{
						push( @exclude, $h );
					}
				}
			}
			next	if( $headers{$head}{comp} == 1 );

			&copy_all_but_include( $tmp, \@exclude, @contents );
			my $out = &headers::do_compile( $tmp, $compile_flags );
			$headers{$head}{comp} = ( (length($out) > 0) ? 0 : 1 );

			if( $headers{$head}{comp} == 0 and $headers{$head}{calls} )
			{
				foreach my $ohead ( @{ $headers{$head}{calls} } )
				{
					$headers{$ohead}{comp} = 1;
				}
			}

			unlink $tmp
			    or warn "Could not delete:\t $tmp\n";
		}

		my $nredund = 0;
		foreach my $head ( @$order ) { $nredund += $headers{$head}{comp}; }

		#	print "nredund = $nredund\n";
		#	foreach my $head ( @$order )
		#	{
		#		print $head ."\t". $headers{$head}{comp} ."\n";
		#		if( $headers{$head}{called} )
		#		{
		#			my @a = @{ $headers{$head}{called} };
		#			print "\tcalled by:\t@a\n";
		#		}
		#		if( $headers{$head}{calls} )
		#		{
		#			my @a = @{ $headers{$head}{calls} };
		#			print "\tcalls    :\t@a\n";
		#		}
		#	}


		if( $nredund )
		{
			&fix_cpp( $cpp, \@contents, \%headers )
				if( not $dry_run );

			print REDUN "$cpp:\t $nhead headers, $nredund redundant, $nmult *multiple\n";
			print REDUN "=" x (length($cpp)+1) ."\n";
			foreach my $head ( keys %headers )
			{
				print REDUN "$head\n"
					if( $headers{$head}{comp} );
			}
			print REDUN "\n";

			if( $mult_headers )
			{
				foreach my $head ( @$mult_headers )
				{
					print REDUN "*$head\n";
				}
				print REDUN "\n";
			}

			++$nfiles;
			$nfails += $nredund;
		}
		$nheads += $nhead;
		$nmuldc += $nmult;
	}
	$nfiles[$thread] += $nfiles;
	$nhead [$thread] += $nheads;
	$nfail [$thread] += $nfails;
	$nmult [$thread] += $nmuldc;

	close REDUN, "> $redundant"
	   or warn "Warning: Could not close:\t $redundant\n";

	$output_files[$thread] = $redundant;

	return;
}



sub concat_thread_output
{
	my ($filename, @output_files) = @_;
	my $out = join(" ", @output_files);
	system("cat $out > $filename");
	system("rm $out");
}




#################################################################################
#										#
#				MAIN PROGRAM					#
#										#
#################################################################################

my( $ncpus, @cpp ) = &Init( @ARGV );

if( not -e "cloudyconfig.h" )
{
	die "Please run 'make' prior to running the script\n";
}

print "The script produces a summary in:\t $redundant_cpp\n";


&headers::ResolveHdr( $ncpus );

#	@cpp = qw/ mole_reactions.cpp /;
&headers::do_parallel( $ncpus, \@cpp, \&cpp_dep );
#	&cpp_dep(0, @cpp);
#	die;


#  Reports...
#
my ($nfiles, $nheads, $nfails, $nmult) = (0, 0, 0, 0);
for( my $i = 0; $i < @nfiles; $i++ )
{
	$nfiles += $nfiles[$i];
	$nheads += $nhead [$i];
	$nfails += $nfail [$i];
	$nmult  += $nmult [$i];
}
printf "%3i of %3i files contained redundant headers\n", $nfiles, scalar(@cpp);
printf "headers included: %4i, redundant: %4i, duplicates: %4i\n", $nheads, $nfails, $nmult;

&concat_thread_output( $redundant_cpp, @output_files );
