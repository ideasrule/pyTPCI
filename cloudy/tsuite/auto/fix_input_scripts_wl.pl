#!/usr/bin/env perl
#
#  fix_input_scripts_wl.pl:
#
#	Invocation:
#		fix_input_scripts_wl.pl <list of .in files>
#  or		fix_input_scripts_wl.pl
#  the second version operates on all .in files present in the current directory.
#
#	The script assumes that Cloudy was run on each input file, so that for
#  each .in file there is an accompanying .out file, e.g., test.in and test.out.
#
#	The script searches the output (.out) files for failures of the Cloudy
#  findline() function, and compiles a list of suggested updates.  These changes
#  are applied to a number of commands and any lists of lines that may appear in
#  the input script.  Files invoked from the script that may contain commands or
#  linelists are also parsed and corrected if needed.
#
#  	The user is prompted to confirm or reject each update, through the options
#  Y,y,N,n.  Capital letters repeat the operation (acceptance or rejection) on
#  subsequent instances of the same change without further prompting.
#
#	A limited number of Cloudy commands is currently supported.  To process
#  another command, one must do two things:
#
#  1) Define two routines of the form "&match_any_cmd()" and "&match_this_cmd()"
#     (see below for examples) that evaluate true when the script line matches
#     any command of that type, or the command with the intended change has been
#     performed already (in a previous run).
#
#  2) Append the @cmds array below with references to these routines.
#
#  Chatzikos, Marios			2013-Jul-25
#  Chatzikos, Marios			2015-Jan-25
#	- Confirmation now needed for each unique change, and applied to all
#	  subsequent instances of that change.
#	- Fixed logical errors with input script updates -- previously done
#	  only when script versbose (=1).
#  Chatzikos, Marios                    2015-Jan-29
#	- Response to confirmation prompts now distinguish between lowercase and
#	  capital letters.  The latter apply to all prompts for the same change.
#	  The former only to the current one.
#  Chatzikos, Marios			2015-Oct-01
#  	- Updated to the new transition mismatch reporting scheme.  It uses two
#  	  lines of output (instead of 3), and different messages.  See function
#  	  get_findline_fails() for details.
#  Chatzikos, Marios			2020-Aug-25
#  	- Updated to search for the wavelength suggestion, instead of assuming it
#  	  is always 3 lines after the emitted warning.
#  Chatzikos, Marios			2021-Mar-17
#  	- Exported functions to module WL.pm
#  	- Updated to report proposed wavelength for lines not encountered in
#  	  a script and its datafiles.  These lines are typically hardcoded in
#  	  Cloudy, e.g., lines involving the Halpha, or Hbeta wavelengths.
#  Chatzikos, Marios			2021-Apr-09
#  	- Enabled processing of 'set blend' command.
#  Chatzikos, Marios			2022-Jun-08
#  	- Enabled processing of data/blends.ini by default.
#  	- Enabled processing of PROBLEMs reported by findline.
#  	  (These are followed by an abort; files have to be updated piecemeal.)
#

use strict;
use warnings;

use lib ".";
use WL;


$| = 1;


my $end_of_lines = qr/^end\s+(of\s+|)line/;
my $datadir = "../../data";
my $blends = "$datadir/blends.ini";

#
#  These are references to functions aimed to match certain commands.  The array
#  is populated below, after the functions are defined.
#
my @cmds;


#
#  Switch on to get detailed report
#
my $verbose = 0;


#
#  Global fix list, used with storing, and later reusing, change confirmations
#
my @Global_FixList;




#################################################################################
#                                                                               #
#                               SET GLOBAL VARIABLES                            #
#                                                                               #
#################################################################################

my $current_filename;
my ($comment, $datestr);

sub set_date
{
	# get the current date
	my @DD = localtime;
	my @abbr = qw( jan feb mar apr may jun jul aug sep oct nov dec );
	$DD[5] %= 100;
	$datestr = sprintf( "%2.2d %s %2.2d", $DD[5], $abbr[$DD[4]], $DD[3] );
}

sub set_comment
{
	print "please give a short text to be appended to the new comments in the input scripts: ";
	$comment = <STDIN>;
	chomp( $comment );
	print   "\n";
}

sub set_current_filename
{
	$current_filename = shift;
}

sub get_global_update
{
	my( $label, $oldwl ) = @_;

	my $update = undef;
	foreach my $hashref ( @Global_FixList )
	{
		next	if( $$hashref{label} ne $label );
		next	if( $$hashref{oldwl} ne $oldwl );
		$update = $$hashref{update};
		#	print "$label, $oldwl, $update\n";
	}

	return	$update;
}

sub set_global_update
{
	my( $label, $oldwl, $update ) = @_;

	my %hash;
	$hash{label} = $label;
	$hash{oldwl} = $oldwl;
	$hash{update} = $update;

	push( @Global_FixList, { %hash } );

	return;
}




#################################################################################
#										#
#				USER INTERACTION				#
#										#
#################################################################################
sub confirm_change
{
	my( $valid_response, $response );
	while( not defined( $valid_response ) )
	{
		print "Confirm change? [Y/y/N/n]  ";
		$response = <STDIN>;
		chomp( $response );
		if( $response =~ m/^y$/i )
		{
			$response = 'y';
			$valid_response = 1;
		}
		elsif( $response =~ m/^n$/i )
		{
			$response = 'n';
			$valid_response = 1;
		}
		else
		{
			print "Unacceptable answer: '$response'.  Please try again...\n";
		}
	}

	return $response;
}




#################################################################################
#										#
#				FIELD PROCESSING				#
#										#
#################################################################################

#
#  Subroutines used to identify script lines of certain commands (*_any_* subs)
#  and commands of a given transition (*_this_*).
#  These are mainly used as input to &process_line().
#

#
#  Routines to match commands....
#
sub match_this
{
	my ($script_line, $label, $wvlng) = @_;
	my $match = 0;
	$match = 1
		if( $script_line =~ /$label/i	and
		    $script_line =~ /$wvlng/i );
	return	$match;
}


sub match_any_monitor_line
{
	my ($script_line) = @_;
	my $match = 0;
	$match = 1
		if( $script_line =~ m/^moni[a-zA-Z]*/i	and
		    $script_line =~ m/line\s+/i );
	return	$match;
}


sub match_this_monitor_line
{
	my ($script_line, $label, $wvlng) = @_;
	my $match = 0;
	$match = 1
		if( &match_any_monitor_line($script_line)	and
		    &match_this($script_line, $label, $wvlng) );
	return	$match;
}



sub match_any_normalize
{
	my ($script_line) = @_;
	my $match = 0;
	$match = 1
		if( $script_line =~ m/^norm[a-zA-Z]*/ );
	return $match;
}



sub match_this_normalize
{
	my ($script_line, $label, $wvlng) = @_;
	my $match = 0;
	$match = 1
		if( &match_any_normalize($script_line)	and
		    &match_this($script_line, $label, $wvlng) );
	return $match;
}



sub match_any_stopline
{
	my ($script_line) = @_;
	my $match = 0;
	$match = 1
		if( $script_line =~ m/^stop\s+line/ );
	return $match;
}



sub match_this_stopline
{
	my ($script_line, $label, $wvlng) = @_;
	my $match = 0;
	$match = 1
		if( &match_any_stopline($script_line) and
		    &match_this($script_line, $label, $wvlng) );
	return $match;
}



#
#  Routines to match line lists....
#

sub match_no_comment_line
{
	my ($script_line) = @_;
	my $match = 0;
	$match = 1
		if( $script_line !~ m/^#/ );
	return	$match;
}



sub match_listed_line
{
	my ($script_line, $label, $wvlng) = @_;
	my $match = 0;
	$match = 1
		if( &match_no_comment_line($script_line) and
		    &match_this($script_line, $label, $wvlng) );
	return	$match;
}



sub match_any_linelist_cmd
{
	my ($script_line) = @_;
	my $match = 0;
	$match = 1
		if( $script_line =~ m/^(opti|prin|save|set blend)/ );
	return	$match;
}



sub match_end_of_lines
{
	my ($script_line) = @_;
	my $match = 0;
	$match = 1
		if( $script_line =~ m/$end_of_lines/ );
	return $match;
}



#
#  Populate the array of function references that can regexp-match a set of
#  commands.  To expand the number of commands implemented, add suitable
#  functions above, and enter references to them at the end of the array.
#
#  NOTE:  These all should refer to commands.  Lists of lines are treated
#  separately.
#
$cmds[0]{any}  = \&match_any_monitor_line;
$cmds[0]{this} = \&match_this_monitor_line;
$cmds[1]{any}  = \&match_any_normalize;
$cmds[1]{this} = \&match_this_normalize;
$cmds[2]{any}  = \&match_any_stopline;
$cmds[2]{this} = \&match_this_stopline;



#
#  The parser may read emission lines in two ways:  As arguments to commands, in
#  which case the species is given between double quotations somewhere in the line,
#  typically followed by the wavelength.  Or as part of a list of emission lines,
#  in which case, the species is given without quotations in columns 1-4 followed
#  by the wavelength.
#
#  Lists can appear in independent files, or in the input script terminated by an
#  "end (of) lines" command.
#
#  The 3 subs below deal with properly reading and reporting the species and line
#  wavelength.
#

#
#  Emulate the parser, and use the first number in the command (possibly ending in
#  a single character unit) as the wavelength of the transition.
#
sub get_quoted_emline
{
	my ($script_line) = @_;

	my @words = split(/\"/, $script_line);

	my $scrpt_label = $words[1];
	$scrpt_label = "\Q$scrpt_label\E";

	my $scrpt_wl;
	@words = split( /\s+/, $script_line );
	foreach my $word ( @words )
	{
		if ($word =~ m/^\d+(\.\d*|)([a-zA-Z]|)$/ )
		{
			$scrpt_wl = $word;
			last;
		}
	}

	return	($scrpt_label, $scrpt_wl);
}



# Parse lines of the form
#	"H  1  1216A come what may"
sub get_listed_emline
{
	my( $script_line ) = @_;

	my ($label) = substr($script_line, 0, 4, " " x 4);
	$label = "\Q$label\E";
	my (undef, $oldwl, $rest) = split(/\s+/, $script_line, 3);

	#	print "script_line:\t \"$script_line\"\n";
	#	print "label= \"$label\"\t oldwl= $oldwl\n";

	return	($label, $oldwl, $rest);
}



sub get_emission_line
{
	my ($script_line, $quoted) = @_;

	if( defined($quoted) and $quoted )
	{
		return	&get_quoted_emline($script_line);
	}
	else
	{
		return	&get_listed_emline($script_line);
	}
}



#################################################################################
#										#
#				LINE PROCESSING					#
#										#
#################################################################################

#
#  Form the new command with the correct wavelength information, and prepend it
#  by a suitable comment line detailing the change.
#
sub fix_line
{
	my ($is_llst, $script_line, $scrpt_wl, $newwl, $update) = @_;

	my $sentinel = "##";
	   $sentinel = "#"	if( $is_llst );

	my $comment_line = "$sentinel >>chng $datestr, wl from $scrpt_wl to $newwl, $comment\n";

	my ($prefix, $suffix) = split(/$scrpt_wl/, $script_line);
	my $new_line = $prefix . $newwl;

	$suffix =~ s/^[a-zA-z]//	if( $suffix =~ m/^[a-zA-Z]/ );
	$new_line .= $suffix;

	print "\n";
	print "$current_filename: Current line:\t$script_line";
	print "$current_filename: Updated line:\t$new_line";

	$update = &confirm_change()
		if( not defined( $update ) );

	if( $update eq 'y' )
	{
		$new_line = $comment_line . $new_line;
	}
	else
	{
		$new_line = $script_line;
	}

	#	print "$script_line";
	#	print "$new_line\n\n";

	return  ( $new_line, $update );
}



#
#  Process an input line and, if needed, replace the old wavelength.  If a match
#  is found, store the matching line for the sake of error checking at the end.
#
#  The subroutine reference $match_this_line checks against the possibility that
#  the change was done in a previous run.  One of the &match_this_*() functions
#  above should be used for that.
#
#  The parameter $is_lab_quoted distinguishes between input in a command line vs
#  in a list of lines.
#
sub process_line
{
	my ($is_llst, $script_line, $match_this_line, $is_lab_quoted, @fix_list) = @_;

	#	print "SCRIPT LINE:\t $script_line";
	my ($scrpt_label, $scrpt_wl) = &get_emission_line( $script_line, $is_lab_quoted );
	#	print "SCRIPT LINE:\t quoted= $is_lab_quoted\t label= \"$scrpt_label\"\t oldwl= \"$scrpt_wl\"\n";

	my ($index, $fix_made) = (-1, 0);
	my ($newwl, $match_old) = (0, 0);
	for( my $i = 0; $i < scalar(@fix_list); $i++ )
	{
		my $label = "\Q$fix_list[$i]{label}\E";
		#	print "\tlabel=\"$label\"\t scrpt_label=\"$scrpt_label\"\n";
		next
			if( "\U$label\E" ne "\U$scrpt_label\E" ); 

		$newwl = $fix_list[$i]{newwl};

		# Change done in a previous run
		if( &{ $match_this_line } ($script_line, $label, $newwl) )
		{
			$fix_list[$i]{found}++;
			push ( @{ $fix_list[$i]{matches} }, $script_line );
			last;
		}

		my $oldwl = $fix_list[$i]{oldwl};
		#	print "\t i= $i\t label = \"$label\"\t oldwl = $oldwl\t newwl = $newwl\n";

		if( &WL::do_wl_match($oldwl, $scrpt_wl) )
		{
			$index = $i;
			$fix_list[$i]{found}++;
			push ( @{ $fix_list[$i]{matches} }, $script_line );
			#	print "$script_line\t label= $label\t oldwl= $oldwl\t index= $index\t newwl= $newwl\n";
			#	print "=>\tslab = $scrpt_label\t swl = $scrpt_wl\t --VS-- \t lab = $label\t oldwl = $oldwl\t newwl = $newwl\t WITH\t match_old = $match_old\n\n";
			last;
		}
	}

	if( $index != -1 )
	{
		my $update = &get_global_update(
				$fix_list[ $index ]{label},
				$fix_list[ $index ]{oldwl} );

		my $new_update;
		( $script_line, $new_update ) =
			&fix_line($is_llst, $script_line, $scrpt_wl, $newwl,
					$update );

		&set_global_update(
				$fix_list[ $index ]{label},
				$fix_list[ $index ]{oldwl},
				$new_update )
				if( not defined( $update ) );

		$fix_made = ( $new_update eq 'y' ? 1 : 0 );
	}
	#	print "FINAL SCRIPT LINE:\t $script_line\n";

	return	($script_line, $fix_made, @fix_list);
}





#################################################################################
#										#
#			ARRAY & FILE MANIPULATION				#
#										#
#################################################################################

sub show_array_of_fixes
{
	my (@fixes) = @_;

	for( my $i = 0; $i < scalar(@fixes); $i++)
	{
		printf "\ti=%3d\t label=\"%s\"\t oldwl=\"%s\"\t newwl=\"%s\"\n",
			$i + 1,
			$fixes[$i]{label},
			$fixes[$i]{oldwl},
			$fixes[$i]{newwl};
	}
	print	"\n";

	return;
}

sub get_file_contents
{
	my ($filename) = @_;

	my @contents = &WL::get_file_contents( $filename );

	&set_current_filename( $filename );

	return  @contents;
}


#
# Get filenames that appear between quotes as part of a command.  Some commands
# expect an input filename after the output filename.  The $offset parameter
# allows picking out the intended one.
#
sub get_aux_files
{
	my ($string, $offset, @contents) = @_;

	my @files = grep(/^$string/, @contents);

	my (@local_files, @remdir_files);
	foreach my $file (@files)
	{
		my @words = split(/\"/, $file);

		$file = undef;
		$file = $words[$offset]
			if( defined($words[$offset]) );

		next
			if( not defined( $file ) );

		if( -s $file )
		{
			push( @local_files, $file );
		}
		else
		{
			push( @remdir_files, "$datadir/$file" )
				if( -s "$datadir/$file" );
		}
	}

	return	(\@local_files, \@remdir_files);
}



sub is_hash_entered
{
	my ($hval, $key, @ahash) = @_;

	my $is_there = 0;

	for( my $i = 0; $i < scalar(@ahash); $i++)
	{
		if( $hval eq $ahash[$i]{$key} )
		{
			$is_there = 1;
			last;
		}
	}

	return	$is_there;
}



sub join_arrays
{
	my ($arr1, $arr2) = @_;
	my @arr1 = @{ $arr1 };
	my @arr2 = @{ $arr2 };

	for( my $i = 0; $i < scalar(@arr2); $i++)
	{
		if( not grep( $arr2[$i], @arr1 ) )
		{
			push( @arr1, $arr2[$i] );
		}
	}

	return	(@arr1);
}





#################################################################################
#										#
#				SCRIPT PROCESSING				#
#										#
#################################################################################

#
# Parse findline error message to extract the label and wavelength.
#
sub parse_findline_error
{
	my( $line ) = @_;

	#
	# Error message looks like:
	#  'PROBLEM findline did not find line with label (between quotes) \
	#   "O  4" and wavelength 1397.20A.'
	#
	my @words = split( /\"/, $line );
	my $label = $words[1];
	@words = split( /\s+/, $words[-1] );
	my $wl = $words[-1];
	$wl =~ s/\.$//;

	return ( $label, $wl );
}


#
# Parse findline warning message to extract the label and wavelength.
#
sub parse_findline_warning
{
	my( $line ) = @_;

	#
	# Warning message looks like:
	#  'WARNING: no exact matching lines found for: "O  4  1397.20A"'
	#
	my @words = split( /\"/, $line );
	my $label = $words[1];
	my $wl = $words[-1];
	$wl =~ s/\s+//g;

	return ( $label, $wl );
}


#
#  Process findline error messages.
#
sub process_findline_messages
{
	my( $script, $contents, $msg_Incorrect_Label, $msg_Suggestion,
		$parse_message, $fix_list ) = @_;

	my %fix;

	my $line = shift( @$contents );

	my( $label, $wl ) = &{ $parse_message } ( $line );
	$fix{label} = $label;
	$fix{oldwl} = $wl;

	while( defined( $contents )
		and @$contents
		and $$contents[0] =~ "\Q$msg_Incorrect_Label\E" )
	{
		shift @$contents;
	}
	$line = shift @$contents;

	if( $line !~ $msg_Suggestion )
	{
		print	"Script: $script\t WARNING!  No closest match found! "
		  .	"Ignore line:\t" . $fix{label} ." @".$fix{oldwl} ."\n";
		return;
	}

	my @words   = split(/\"/, $line);
	@words	    = split(/\s+/, $words[1]);
	$fix{newwl} = $words[-1];


	# Needed for error checking...
	$fix{found} = 0;
	$fix{matches} = [];

	#	print "label= ". $fix{label} ."\t oldwl = ". $fix{oldwl} ."\t ". $fix{newwl} ."\n";
	
	my $is_listed = 0;
	for( my $i = 0; $i < scalar(@$fix_list); $i++ )
	{
		if( $$fix_list[$i]{label} eq $fix{label}  and
		    &WL::do_wl_match($$fix_list[$i]{oldwl}, $fix{oldwl}) and
		    &WL::do_wl_match($$fix_list[$i]{newwl}, $fix{newwl}) )
		{
			$is_listed++;
			last;
		}
	}

	if( not $is_listed )
	{
		my ($tol, $error) = &WL::wldiff_error( $fix{oldwl}, $fix{newwl} );
		if( $error > $tol and $verbose )
		{
			print	"\tScript: $script\t WARNING!\t Large WL diff:\t";
			printf	"\tlabel= %s\t oldwl= %s\t newwl= %s\t error= %.4f vs tol= %.4f\n",
				$fix{label}, $fix{oldwl}, $fix{newwl}, $error, $tol;
			#	print	"\tWARNING! Double check line this line.\n";
		}
		push( @$fix_list, { %fix } );
	}
}


#
#  Collect the unique findline failures from the main Cloudy output file.
#
sub get_findline_fails
{
	my ($script, $output) = @_;

	my @contents = &get_file_contents( $output );
	my @fix_list;

	my $FindLine_Error_Mesg = "PROBLEM findline did not find line with label";
	my $FindLine_Error_Incorrect_Label = "The closest line \(any label\) was";
	my $FindLine_Error_Suggestion = "The closest with correct label was";

	my $FindLine_Warn_Mesg = "WARNING: no exact matching lines found for";
	my $FindLine_Warn_Incorrect_Label = "WARNING: Line with incorrect label found";
	my $FindLine_Warn_Suggestion = "Taking best match as";

	while( @contents )
	{
		if( $contents[0] =~ $FindLine_Error_Mesg )
		{
			&process_findline_messages( $script,
						\@contents,
						$FindLine_Error_Incorrect_Label,
						$FindLine_Error_Suggestion,
						\&parse_findline_error,
						\@fix_list );
		}
		elsif( $contents[0] =~ $FindLine_Warn_Mesg )
		{
			&process_findline_messages( $script,
						\@contents,
						$FindLine_Warn_Incorrect_Label,
						$FindLine_Warn_Suggestion,
						\&parse_findline_warning,
						\@fix_list );
		}
		else
		{
			shift( @contents );
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



#
#  Problems with disambiguation:
#  In some cases, "H  1" lines and their "Ca B" equivalents have
#  different recommended wavelength fixes.
#  Bring these to the user's attention, and change the "Ca B" new
#  wavelength to be in agreement with the "H  1" one.
#
sub process_H1_caseb
{
	my (@fix_list) = @_;

	my $ncab = 0;
	for( my $i = 0; $i < scalar(@fix_list); $i++ )
	{
		next
			if( $fix_list[$i]{label} !~ m/^Ca (A|B)$/i );

		for( my $j = 0; $j < scalar(@fix_list); $j++ )
		{
			if( $fix_list[$j]{label} =~ m/^H  1$/i	and
				&WL::do_wl_match($fix_list[$j]{oldwl}, $fix_list[$i]{oldwl})	and
			    not &WL::do_wl_match($fix_list[$j]{newwl}, $fix_list[$i]{newwl}) )
			{
				print	"\t  => Setting new wl for:  "
				  .	"\"". $fix_list[$i]{label} ."\" @". $fix_list[$i]{oldwl} .","
				  .	"\t from: "             . $fix_list[$i]{newwl}
				  .	"\t   to: "             . $fix_list[$j]{newwl}
				  .	"\t to match with \""   . $fix_list[$j]{label} ."\"\n";
				$fix_list[$i]{newwl} = $fix_list[$j]{newwl};
				$ncab = 1;
				last;
			}
		}
	}

	print "\n"
		if( $ncab > 0 );
	
	return  @fix_list;
}



#
#  Loop over and fix an array of lines.  The pattern matching sub reference
#  $match_any_line should be one of the &match_any_*() functions above.
#
sub fix_contents
{
	my ($is_llst, $contents, $fix_list, $match_any_line, $match_this_line, $is_lab_quoted) = @_;

	my $nfixes = 0;
	my @newcontents;
	foreach my $script_line ( @$contents )
	{
		if( &{ $match_any_line } ($script_line) )
		{
			# A line label in a list may be quoted
			my $is_lab_quoted_this = $is_lab_quoted;
			if( $script_line =~ m/^"/ )
			{
				$is_lab_quoted_this = 1;
			}

			my $fix = 0;
			($script_line, $fix, @$fix_list) = &process_line( $is_llst, $script_line, $match_this_line, $is_lab_quoted_this, @$fix_list );
			$nfixes += $fix;
		}

		push( @newcontents, $script_line );
	}

	return	($nfixes, \@newcontents, $fix_list);
}



#
#  Invoke the general &fix_contents() routine to carry out the wavelength
#  updates for the command specified.
#
sub fix_commands
{
	my ($is_llst, $content_ref, $fix_list_ref, $match_any_cmd, $match_this_cmd) = @_;

	my $is_lab_quoted = 1;

	return	&fix_contents($is_llst, $content_ref, $fix_list_ref, $match_any_cmd, $match_this_cmd, $is_lab_quoted);
}



#
#  Fix list of lines...
#
sub fix_lines_listed
{
	my ($is_llst, $content_ref, $fix_list_ref) = @_;

	my $is_lab_quoted = 0;

	return &fix_contents($is_llst, $content_ref, $fix_list_ref, \&match_no_comment_line, \&match_listed_line, $is_lab_quoted);
}



#
#  Fix groups of lines in an input script, terminated by 'end (of) lines'.
#
sub fix_script_linelist
{
	my ($is_llst, $contents_ref, $fix_list_ref) = @_;
	my @contents = @{ $contents_ref };


	my $start = 0;
	my $nfixes_tot = 0;
	my @newcontents;
	for( my $i = 0; $i < scalar(@contents); $i++ )
	{
		if( &match_end_of_lines( $contents[ $i ] ) )
		{
			my $end = --$i;
			do { $end--; } while( $end >= 0 and not &match_any_linelist_cmd( $contents[$end] ) );
			my @tmp = @contents[ ++$end .. $i ];
			my $nfixes = 0;
			($nfixes, $contents_ref, $fix_list_ref) = &fix_lines_listed($is_llst, \@tmp, $fix_list_ref);
			$nfixes_tot += $nfixes;
			push( @newcontents, @contents[ $start..($end-1)] );
			push( @newcontents, @{ $contents_ref } );
			push( @newcontents, $contents[++$i]    );
			$start = $i+1;
			next;
		}
	}

	push( @newcontents, @contents[ $start..$#contents ] );

	return	($nfixes_tot, \@newcontents, $fix_list_ref);
}



#
#  Process a script of commands, that may also contain lists of lines terminated
#  by the "end of lines" sentinel, to update wavelengths to the values suggested
#  by Cloudy.  The commands processed are listed in the @cms array of hashes, see
#  above.
#
sub fix_script
{
	my (undef, $contents_ref, $fix_list_ref) = @_;

	my $is_llst = 0;
	my $nfixes_tot = 0;

	#  Commands
	#
	for( my $icmd = 0; $icmd < scalar(@cmds); $icmd++ )
	{
		my %cmd = %{ $cmds[$icmd] };
		if( grep { &{ $cmd{any} } ( $_ ) }  @{ $contents_ref } )
		{
			my $nfixes = 0;
			($nfixes, $contents_ref, $fix_list_ref) = &fix_commands( $is_llst, $contents_ref, $fix_list_ref, $cmd{any}, $cmd{this} );
			$nfixes_tot += $nfixes;
		}
	}

	#  Lists of lines
	#
	if( grep { &match_end_of_lines( $_ ) } @{ $contents_ref } )
	{
		my $nfixes = 0;
		($nfixes, $contents_ref, $fix_list_ref) = &fix_script_linelist( $is_llst, $contents_ref, $fix_list_ref );
		$nfixes_tot += $nfixes;
	}


	return	($nfixes_tot, $contents_ref, $fix_list_ref);
}



#
#  Apply a fix on a list of files.  The kind of fix is determined by the
#  $file_fix_sub sub-ref.  This is useful for treating external files, such as
#  those found in 'init' and 'save line' commands.  The parameter $is_llst
#  determines the comment sentinel: '##' for scripts, and '#' for line lists.
#
sub fix_filelist
{
	my ($is_llst, $llist_ref, $fix_list_ref, $file_fix_sub) = @_;
	my @llist    = @{ $llist_ref    };

	my (@llists, $content_ref, %fix);
	foreach my $llist ( @llist )
	{
		next
			if( &is_hash_entered($llist, "file", @llists) );

		my @contents = &get_file_contents( $llist );

		my $nfixes_tot = 0;
		($nfixes_tot, $content_ref, $fix_list_ref) = &{ $file_fix_sub }( $is_llst, \@contents, $fix_list_ref );
		@contents = @{ $content_ref };

		if( $nfixes_tot > 0 )
		{
			$fix{file} = $llist;
			$fix{newfile} = &WL::write_file_contents( $llist, @contents );

			push( @llists, { %fix } );
		}
	}

	return	(\@llists, $fix_list_ref);
}



#
#  Final reporting...
#
sub check_fixed_wl
{
	my ($fix_list, $silent) = @_;

	return	if not defined $fix_list;

	$silent = 0	if not defined( $silent );

	my $all_found = 0;
	for (my $i = 0; $i < @$fix_list; $i++)
	{
		my $found = $$fix_list[$i]{found};
		my $newwl = $$fix_list[$i]{newwl};

		if( $found == 1 )
		{
			$all_found++;
			splice(@$fix_list, $i, 1);
			$i--;
			next;
		}
		elsif( $found == 0 )
		{
			print	"\t No matches for line:\t"
			  .	"\"". $$fix_list[$i]{label} ."\" ". $$fix_list[$i]{oldwl}
			  .	"\t\t=> ". $$fix_list[$i]{newwl} ."\n"
			  if not $silent;
		}
		elsif( $found > 1 )
		{
			if( 0 )
			{
				print "\t WARNING:\t Multiple matches for line:\t"
				  .	"\"". $$fix_list[$i]{label} ."\" @". $$fix_list[$i]{oldwl} ."\n";

				my @matches = @{ $$fix_list[$i]{matches} };
				for (my $i = 0; $i < scalar(@matches); $i++)
				{
					print "\t\t =>\t ". $matches[$i]
						if( $matches[$i] !~ m/$newwl/ );
				}
				print "\n";
			}

			$all_found++;
			splice(@$fix_list, $i, 1);
			$i--;
		}
	}

	return;
}


#
#  Warn about unmatched wavelengths during update.  These may be due to updates
#  done to linelists or init files through a different script.
#
sub warn_unmatched_emlines
{
	my( $fix_list, $script, $newscript, $local_inits, $remdir_inits, $llists, $remdir_llists ) = @_;

	return
		if( not @$fix_list );

	if( $verbose )
	{
		print "\n\t WARNING: Could not match ". @$fix_list ." lines.\n";
		print "\t Inspect scripts and init/linelists:\n";
		print "\t     SCRIPT:\t $script \t  VS  \t $newscript\n";

		for( my $i = 0; $i < @$local_inits; $i++ )
		{
			print "\t      INIT:\t ". $$local_inits[$i]{file} ."\n";
		}
		for( my $i = 0; $i < @$remdir_inits; $i++ )
		{
			print "\t DATA-INIT:\t ". $$remdir_inits[$i]{file} ."\n";
		}
		for( my $i = 0; $i < @$llists; $i++ )
		{
			print "\t      LLIST:\t ". $$llists[$i]{file} ."\n";
		}
		for( my $i = 0; $i < @$remdir_llists; $i++ )
		{
			print "\t DATA-LLIST:\t ". $$remdir_llists[$i]{file} ."\n";
		}
	}
	else
	{
		print "\t Script: $script\t unmatched wavelengths:\t ". @$fix_list ."\n";
	}
}



#
#  Perform file updates where applicable
#
sub update_files
{
	my( $script, $newscript, $local_inits, $remdir_inits, $llists, $remdir_llists ) = @_;

	rename "$newscript", "$script"
		if( defined( $newscript ) );
	for( my $i = 0; $i < @$local_inits; $i++ )
	{
		#	print "\t init: ". $$local_inits[$i]{file} ."\n";
		rename	$$local_inits[$i]{newfile}, $$local_inits[$i]{file};
	}
	for( my $i = 0; $i < @$remdir_inits; $i++ )
	{
		#	print "\t remdir_inits: ". $$remdir_inits[$i]{file} ."\n";
		rename	$$remdir_inits[$i]{newfile}, $$remdir_inits[$i]{file};
	}
	for( my $i = 0; $i < @$llists; $i++ )
	{
		#	print "\t llist: ". $$llists[$i]{file} ."\n";
		rename	$$llists[$i]{newfile}, $$llists[$i]{file};
	}
	for( my $i = 0; $i < @$remdir_llists; $i++ )
	{
		#	print "\t remdir_llists: ". $$remdir_llists[$i]{file} ."\n";
		rename	$$remdir_llists[$i]{newfile}, $$remdir_llists[$i]{file};
	}
}





#################################################################################
#										#
#				MAIN PROGRAM					#
#										#
#################################################################################

&set_date();
&set_comment();

my @scripts = @ARGV;
@scripts = glob "*.in"	if( not @scripts );

foreach my $script ( @scripts )
{
	if( $script !~ m/\.in$/ )
	{
		warn	"WARNING: Not an input script:\t $script\n";
		next;
	}

	my $output = $script;
	$output =~ s/\.in$/.out/;

	next
		if( not -s $output );

	my @fix_list = &get_findline_fails( $script, $output );
	@fix_list = &process_H1_caseb( @fix_list );
	&show_array_of_fixes( @fix_list )
		if( $verbose );


	my @contents = &get_file_contents( $script );
	my ($nfixes_tot, $contents, $fix_list_ref) = &fix_script( undef, \@contents, \@fix_list );
	my $newscript;
	$newscript = &WL::write_file_contents( $script, @$contents )
		if( $nfixes_tot > 0 );


	#
	#  Process auxiliary files
	#

	#
	#  External files that contain commands
	#
	my ($local_ref, $remdir_ref) = &get_aux_files("init", 1, @$contents);

	($local_ref, $fix_list_ref) = &fix_filelist( 0, $local_ref, $fix_list_ref, \&fix_script );
	my @local_inits = @{ $local_ref };

	($local_ref, $fix_list_ref) = &fix_filelist( 0, $remdir_ref, $fix_list_ref, \&fix_script );
	my @remdir_inits = @{ $local_ref };


	#
	#  External files that contain lists of lines
	#
	($local_ref, $remdir_ref) = &get_aux_files(qr/table\s+line/, 1, @$contents);
	my @local_tline  = @{ $local_ref  };
	my @remdir_tline = @{ $remdir_ref };

	($local_ref, $remdir_ref) = &get_aux_files(qr/save\s+line/, 3, @$contents);
	my @local_llist  = @{ $local_ref  };
	my @remdir_llist = @{ $remdir_ref };

	# Join lists of linelists
	#
	@local_llist  = &join_arrays(\@local_llist , \@local_tline );
	@remdir_llist = &join_arrays(\@remdir_llist, \@remdir_tline);
	@local_tline = @remdir_tline = ();


	my (@llists, @remdir_llists);
	($local_ref, $fix_list_ref) = &fix_filelist( 1, \@local_llist, $fix_list_ref, \&fix_lines_listed );
	@llists = @{ $local_ref };

	($local_ref, $fix_list_ref) = &fix_filelist( 1, \@remdir_llist, $fix_list_ref, \&fix_lines_listed );
	@remdir_llists = @{ $local_ref };

	&check_fixed_wl( $fix_list_ref, 1 );

	#
	# data/blends.ini
	#   This file is used by default -- no command is used to specify it.
	#   Must check if there are still wavelengths to be fixed.
	#
	if( @{ $fix_list_ref } )
	{
		my @blends = ( $blends );
		($local_ref, $fix_list_ref) =
			&fix_filelist( 0, \@blends, $fix_list_ref, \&fix_lines_listed );
		push( @remdir_llists, @{ $local_ref } );
	}

	#
	#  Done
	#
	&check_fixed_wl( $fix_list_ref );

	#
	#  Warn of unmatched wavelengths, but don't prevent file updates.
	#  These emission lines may be contained in an init or linelist file
	#  that has been already updated through a different script.
	#
	&warn_unmatched_emlines( $fix_list_ref, $script, $newscript,
				 \@local_inits, \@remdir_inits, \@llists, \@remdir_llists );

	&update_files( $script, $newscript,
				 \@local_inits, \@remdir_inits, \@llists, \@remdir_llists );

#	die;
}
