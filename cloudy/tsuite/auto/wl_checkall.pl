#!/usr/bin/env perl
#
#  wl_checkall.pl:
# 	Check all output files and write all unique wavelength misses per sim to
#  file.  Totals (including unique lines across all sims) are printed at stdout.
#
#  Chatzikos, Marios			2013-Jul-25
#  Chatzikos, Marios			2021-Mar-17
#  	- Exported functions to module WL.pm
#  	- Refactored all calls into sub main()
#

use strict;
use warnings;

use lib ".";
use WL;

my $outfile = "wvlng.txt";



sub main
{
	open OUTFILE, "> $outfile"
	  or die "Error: Could not open: $outfile\n";
	
	
	my @AllFixes;
	my ($nlines, $nfiles, $nuniq) = (0, 0, 0);
	foreach my $output ( glob "*.out" )
	{
		my $script = $output;
		$script =~ s/\.in$/.out/;

		my @fix_list = &WL::get_findline_fails( $script, $output, 0 );
		if( @fix_list )
		{
			for( my $i = 0; $i < scalar(@fix_list); $i++)
			{
				print OUTFILE "$output: ". $fix_list[$i]{error};
			}
	
			$nlines += @fix_list;
			$nfiles++;
	
			for (my $i = 0; $i < scalar(@fix_list); $i++)
			{
				my $found = 0;
				for( my $ifix = 0; $ifix < scalar(@AllFixes); $ifix++ )
				{
					if( $AllFixes[$ifix]{label} eq $fix_list[$i]{label} and
					    &WL::do_wl_match($AllFixes[$ifix]{oldwl}, $fix_list[$i]{oldwl} ) and
					    &WL::do_wl_match($AllFixes[$ifix]{newwl}, $fix_list[$i]{newwl} ) )
					{
						$found++;
						last;
					}
				}
	
				if( not $found )
				{
					my %fix;
					$fix{label} = $fix_list[$i]{label};
					$fix{oldwl} = $fix_list[$i]{oldwl};
					$fix{newwl} = $fix_list[$i]{newwl};
					push( @AllFixes, { %fix } );
					$nuniq++;
				}
			}
		}
	}
	
	close OUTFILE
	   or warn "Warning: Could not close: $outfile\n";
	
	
	print	"\nLooking for unmatched wavelengths:";
	if( -s $outfile )
	{
		print	"\nWARNING! In $nfiles sims:"
		  .	" $nlines unidentified wavelengths"
		  .	" ($nuniq unique across all sims).\n";
		print	  "WARNING! See file:  $outfile\n";
	}
	else
	{
		print	" good, none found.\n";
	}
	
	print	"\n";
}

&main();
