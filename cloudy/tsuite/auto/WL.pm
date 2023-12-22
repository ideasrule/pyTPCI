#!/usr/bin/env perl
#
#  Provide functionality common to wl_checkall.pl and fix_input_scripts_wl.pl
#
#  Chatzikos, Marios			2021-Mar-17
# 	- Initial version; refactored from code in these scripts to address
# 	  a bug that came up from the 2020-Aug-25 fix to the fix script not
# 	  being fully ported to the check script.

package WL;

use strict;
use warnings;

$| = 1;



#################################################################################
#										#
#				FIELD PROCESSING				#
#										#
#################################################################################
#
#
#  Conversion factor for wavelength comparison.  This is required to match the
#  intended wavelength reported by Cloudy with the one that appears in the script
#  itself.
#  The factor is meant to be applied to the first of the input wavelengths.
#
sub convert_wl
{
	my ($wl_unit1, $wl_unit2) = @_;


	my $conv_fact = 1;
	if( $wl_unit2 ne $wl_unit1 )
	{
		if( $wl_unit2 eq "A" or $wl_unit2 eq "")
		{
			$conv_fact = 1e4	if( $wl_unit1 eq "m" );
			$conv_fact = 1e6	if( $wl_unit1 eq "c" );
		}
		elsif( $wl_unit2 eq "m" )
		{
			$conv_fact = 1e-4	if( $wl_unit1 eq "A" or $wl_unit1 eq "" );
			$conv_fact = 1e2	if( $wl_unit1 eq "c" );
		}
		elsif( $wl_unit2 eq "c" )
		{
			$conv_fact = 1e-6	if( $wl_unit1 eq "A" or $wl_unit1 eq "" );
			$conv_fact = 1e-2	if( $wl_unit1 eq "m" );
		}
	}
	#	print "$wl_unit2\t  VS  \t$unit_oldwl\t $conv_fact:\t $scrpt_wl\t VS ".	($oldwl * $conv_fact)	."\n";

	return	$conv_fact;
}



#
#  Process wavelength to report the numerical value, the unit, and the number of
#  significant figures found.
#
sub wl_proc
{
	my ($wl) = @_;

	my $wl_unit = "";
	$wl_unit = substr($wl, length($wl)-1, 1, "")	if ($wl =~ m/[a-zA-Z]$/);

	my $idot = index( $wl, "." );
	my $nsigfig = 0;
	if( $idot >= 0 )
	{
		$nsigfig = length(substr($wl, $idot+1));
	}

	return	($wl, $wl_unit, $nsigfig);
}



#
#  Given two wavelengths (including units), produce a string of the first
#  wavelength that best matches the format of the second wavelength.  This is
#  used to pattern match the wavelength reported by Cloudy against the one
#  that exists in the script (likely of lower precision).
#
sub wl_match
{
	my ($wl1, $wl2) = @_;

	my ($wl1str, $wl1_unit, $wl1_nsigfig) = &wl_proc($wl1);
	my ($wl2str, $wl2_unit, $wl2_nsigfig) = &wl_proc($wl2);
	my $conv_fact = &convert_wl($wl1_unit, $wl2_unit);

	return	sprintf("%.*f", $wl2_nsigfig, $wl1str * $conv_fact);
}



#
#  Report if two wavelengths match after taking into account differences in
#  units and precision.
#
sub do_wl_match
{
	my ($wl1, $wl2) = @_;

	my (undef, $wl2_unit) = &wl_proc($wl2);
	my $match_wl = &wl_match($wl1, $wl2);

	my $match = 0;
	$match = 1	if( $wl2 eq "$match_wl$wl2_unit" );
	#	print "$wl1 \t $wl2\t $match_wl\t $match\n";

	return	($match);
}



#
#  Confirm that the new wl is within Cloudy's standard error (1e-4).
#
sub wldiff_error
{
	my ($oldwl, $newwl) = @_;

	my ($wl1, $wl1_unit) = &wl_proc($oldwl);
	my ($wl2, $wl2_unit) = &wl_proc($newwl);
	my $conv_fact = &convert_wl($wl1_unit, $wl2_unit);

	$wl2 *= $conv_fact;

	my $tol = int( log($wl1) / log(10.0) );
	$tol = 0.5 * 10.0**($tol-4.0);

	my $error = abs($wl1 - $conv_fact*$wl2);

	return	($tol, $error);
}



#################################################################################
#										#
#			ARRAY & FILE MANIPULATION				#
#										#
#################################################################################

sub get_file_contents
{
	my ($filename) = @_;

	open FILE, "< $filename" or die  "Could not open:\t $filename\n";
	my @contents = <FILE>;
	close FILE		or warn "Could not close:\t $filename\n";

	return	@contents;
}

sub write_file_contents
{
	my ($fname, @newcontents) = @_;

	$fname .= ".new";

	open  NEWFILE, "> $fname" or die  "Could not open tmp sim:\t $fname\n";
	print NEWFILE @newcontents;
	close NEWFILE		 or print "Could not close tmp sim:\t $fname\n";

	return	$fname;
}



#################################################################################
#										#
#				SCRIPT PROCESSING				#
#										#
#################################################################################

#
#  Collect the unique findline failures from the main Cloudy output file.
#
sub get_findline_fails
{
	my ($script, $output, $verbose) = @_;

	my @contents = &get_file_contents( $output );
	my @fix_list;

	# Uncomment for old reporting scheme of 3 lines of test per mismatch.
	#	my $FindLine_Fail_Mesg = "PROBLEM findline did not find line with label";
	#	my $FindLine_Suggestion = "The closest with correct label was";
	my $FindLine_Fail_Mesg = "WARNING: no exact matching lines found for";
	my $FindLine_Incorrect_Label = "WARNING: Line with incorrect label found";
	my $FindLine_Suggestion = "Taking best match as";

	for( my $i = 0; $i < scalar(@contents); $i++ )
	{
		my $line = $contents[$i];

		if( $line =~ $FindLine_Fail_Mesg )
		{
			chomp($line);
			$line =~ s/\.$//;

			my %fix;

			my @words = split(/\"/, $line);
			$fix{error} = $line;
			$fix{label} = $words[1];
			$fix{oldwl} = $words[2];
			$fix{oldwl} =~ s/\s+//g;
			#	@words      = split(/\s+/, $words[2]);
			#	$fix{oldwl} = $words[3];
			#	$fix{oldwl} =~ s/\.$//;


			++$i;
			while( $contents[$i] =~ $FindLine_Incorrect_Label )
			{
				++$i;
			}
			$line = $contents[$i];

			if( $line !~ $FindLine_Suggestion )
			{
				print	"Script: $script\t WARNING!  No closest match found! "
				  .	"Ignore line:\t" . $fix{label} ." @".$fix{oldwl} ."\n";
				next;
			}

			@words	    = split(/\"/, $line);
			@words	    = split(/\s+/, $words[1]);
			$fix{newwl} = $words[$#words];


			# Needed for error checking...
			$fix{found} = 0;
			$fix{matches} = [];


			#	print "label= ". $fix{label} ."\t oldwl = ". $fix{oldwl} ."\t ". $fix{newwl} ."\n";
			
			my $is_listed = 0;
			for( my $i = 0; $i < scalar(@fix_list); $i++ )
			{
				if( $fix_list[$i]{label} eq $fix{label}  and
				    &do_wl_match($fix_list[$i]{oldwl}, $fix{oldwl}) and
				    &do_wl_match($fix_list[$i]{newwl}, $fix{newwl}) )
				{
					$is_listed++;
					last;
				}
			}



			if( not $is_listed )
			{
				my ($tol, $error) = &wldiff_error( $fix{oldwl}, $fix{newwl} );
				if( $error > $tol and $verbose )
				{
					print	"\tScript: $script\t WARNING!\t Large WL diff:\t";
					printf	"\tlabel= %s\t oldwl= %s\t newwl= %s\t error= %.4f vs tol= %.4f\n",
						$fix{label}, $fix{oldwl}, $fix{newwl}, $error, $tol;
					#	print	"\tWARNING! Double check line this line.\n";
				}
				$fix{error} .= "\t closest wavelength:\t $fix{newwl}\n";
				push( @fix_list, { %fix } );
			}
		}
	}

	if( 0 )
	{
		for my $f ( @fix_list )
		{
			for ( sort keys %$f )
			{
				print "$_:\t$$f{$_}\n";
			}
			print "-" x 20 ."\n";
		}
		die;
	}

	return	@fix_list;
}



1;
