#!/usr/bin/perl -w
#
# this script looks for problems in the Cloudy test suite it looks for botched
# monitors, warnings and for simulations that did not end. all serious errors
# will be reported in serious.txt all lines with PROBLEM will be in minor.txt
# sims that are close to thowing a monitor (but OK) are in close.txt, those
# with DEBUG print statements or TestCode calls are in debug.txt, and those
# with WARNING print statements are in warnings.txt
#
# execute this script from within the directory where the *.out files live.
# The command is: ./checkall.pl
#
# Peter van Hoof

# this was produced by the optimizers and is not really part of the test suite
# its presence will fool this script....
if( -r "optimal.in" )
{
	unlink "optimal.in";
}

print( "\nThis script checks for problems in the output from the Cloudy test suite.\n");
print( "  It checks for three things in the simulations:\n");
print( "  1) Botched monitors or W- messages:  These are serious problems.\n");
print( "  2) Whether the PROBLEM, WARNING, or DEBUG string were printed.\n");
print( "  3) Crashes, in which the code did not end - these are also serious problems.\n");

$nMod = 0;
$nOK = 0;
$nSkip = 0;
$nFail = 0;
$nCrash = 0;
$nBotch = 0;
$nWarn = 0;
$nNoAssert = 0;
$nMinor = 0;
$nDebug = 0;
$nWarnMsg = 0;
$testgrid = "";

open( CRASHED, ">crashed.txt" );
open( SERIOUS, ">serious.txt" );
open( BOTCH, ">serious.asr" );
open( PERFORMANCE, ">performance.asr" );
open( MINOR, ">minor.txt" );
open( CLOSE, ">close.txt" );
open( DEBUG, ">debug.txt" );
open( SKIP, ">skip.txt" );
open( WARNING, ">warnings.txt" );

# loop over all input scripts
while( defined( $input = glob("*.in") ) )
{
	# count the number of simulations
	++$nMod;

	# form name of output file
	$output = $input;
	$output =~ s/\.in/\.out/gi;

	if( -r $output )
	{
		$checkEnd = 0;
		$delimiter = 0;
		$monitors = 0;
		$nWarnMessage = 0;
		$lgCrash = 0;
		$lgFailure = 0;
		$lgBotch = 0;
		$lgAtmosNotFound = 0;
		open( OUTP, $output );
		while( <OUTP> )
		{
# find "DISASTER", "Botched" and "W-" strings in all the output files and write
# corresponding line into serious.txt file.
			if( /DISASTER/ || /Botched/ || /W-/ )
			{
				print SERIOUS "$output: $_";
				if( /W-/ )
				{
					++$nWarnMessage;
				}
			}
			if( /botch/ )
			{
				if( /itrz/ || /zone/ || /iter/ )
				{
					print PERFORMANCE "$output: $_";
				}
				else
				{
					print BOTCH "$output: $_";
				}
				$lgBotch = 1;
			}
# now look for string PROBLEM, which indicates an internal problem
# during the calculation - this is not a serious problem
			if( /PROBLEM/ && ! /DISASTER/ )
			{
				print MINOR "$output: $_";
				++$nMinor;
			}
			if( /ChkMonitor ----/ )
			{
				print CLOSE "$output: $_";
			}
			if( /DEBUG/ || /Test code/ )
			{
				print DEBUG "$output: $_";
				++$nDebug;
			}
# check output files to see whether they contain any monitors output
			if( /ChkMonitor/ )
			{
				++$monitors;
			}
# determine if the program crashed.  
			if( /Cloudy exited OK\]/ || /something went wrong\]/ )
			{
				++$checkEnd;
			}
			if( /something bad has happened./ )
			{
				$lgCrash = 1;
			}
# this is the last line, so lgBotch and nWarnMessage are already determined.
			if( /something went wrong\]/ && !$lgBotch && $nWarnMessage == 0 )
			{
				$lgFailure = 1;
			}
			if( /Error: stellar atmosphere file not found./ )
			{
				$lgAtmosNotFound = 1;
			}
# this counts the number of sims in a grid
			if( /GRID_DELIMIT/ )
			{
				++$delimiter;
			}
			if( /warning/i )
			{
				print WARNING "$output: $_";
				++$nWarnMsg;
			}
		}
		close( OUTP );

# non-grid runs have no delimiter, so pretend there was a single delimiter
		if( $delimiter == 0 )
		{
			$delimiter = 1;
		}

# special hack to check whether the output from a repeated model reproduces exactly...
		$script = "test_grid.pl";
		if( -x $script && $output eq "func_test_cycle.out" )
		{
			$testgrid .= "\nNow I will check whether the output in $output repeats exactly:\n";
			$testgrid1 = `./$script $output xx 4`;
			if( $testgrid1 =~ /PROBLEM/ )
			{
				print SERIOUS "PROBLEM: the output in $output did not repeat exactly\n";
				$lgFailure = 1;
			}
			$testgrid .= $testgrid1;
		}
		if( -x $script && $output eq "func_test_grid.out" )
		{
			$testgrid .= "\nNow I will check whether the output in $output repeats exactly:\n";
			$testgrid2 = `./$script $output xy 1`;
			if( $testgrid2 =~ /PROBLEM/ )
			{
				print SERIOUS "PROBLEM: the output in $output did not repeat exactly\n";
				$lgFailure = 1;
			}
			$testgrid .= $testgrid2;
		}
		if( -x $script && $output eq "func_test_cycle2.out" )
		{
			$testgrid .= "\nNow I will check whether the output in $output repeats exactly:\n";
			$testgrid3 = `./$script $output xz 48`;
			if( $testgrid3 =~ /PROBLEM/ )
			{
				print SERIOUS "PROBLEM: the output in $output did not repeat exactly\n";
				$lgFailure = 1;
			}
			$testgrid .= $testgrid3;
		}

		if( $lgAtmosNotFound )
		{
			++$nSkip;
			print SKIP "$input: stellar atmosphere file is not installed\n";
		}
		elsif( $checkEnd != $delimiter || $lgCrash )
		{
			++$nCrash;
			print SERIOUS "$output: crashed\n";
			print CRASHED "$output: crashed\n";
			$lgSeriousProblems = 1;
		}
		elsif( $lgFailure )
		{
			++$nFail;
			print SERIOUS "$output: returned with FAILURE status\n";
			$lgSeriousProblems = 1;
		}
		elsif( $monitors == 0 )
		{
			++$nNoAssert;
			print SERIOUS "$output: no monitors found\n";
			print CRASHED "$output: no monitors found\n";
			$lgSeriousProblems = 1;
		}
		elsif( $nWarnMessage > 0 )
		{
			++$nWarn;
		}
		elsif( $lgBotch )
		{
			++$nBotch;
		}
		else
		{
			++$nOK;
		}
	}
	else
	{
		++$nSkip;
		print SKIP "$input: no corresponding output was found\n";
	}
}

close( CRASHED );
close( SERIOUS );
close( BOTCH );
close( PERFORMANCE );
close( MINOR );
close( CLOSE );
close( DEBUG );
close( SKIP );
close( WARNING );

if( $nMod != $nOK+$nSkip+$nCrash+$nFail+$nWarn+$nBotch+$nNoAssert )
{
	print STDERR "\nINTERNAL ERROR: the number of sims don't add up!!!!!\n\n";
}

print( "\nLooking for crashes, botched results, and W- messages:");
# if this file has non-zero length, we detected a problem
if( -s "serious.txt" )
{
	print "\nWARNING! crashes, failures, botches, or W- messages were found. This is a serious problem.\n";
	if( $nCrash > 0 )
	{
		print "Cloudy crashed in $nCrash sims.\n";
	}
	print "Check the file serious.txt for names, serious.asr for botched monitors.\n";
	system "grep Cloudy < serious.txt";
	$lgSeriousProblems = 1;
}
else
{
	print " good, none found.\n";
	$lgSeriousProblems = 0;
}

print( "\nLooking for PROBLEM string in the output files:");
if( -s "minor.txt" )
{
	print "\nPROBLEM string was found.  Check minor.txt for a list of names.\n";
	print "This can occur a few times and is not a serious concern.\n";
	print "There are $nMinor lines of text in the file minor.txt.\n";
	$lgMinorProblems = 1;
}
else
{
	print " good, none found.\n";
	$lgMinorProblems = 0;
}

print( "\nLooking for WARNING string in the output files:");
if( -s "warnings.txt" )
{
	print "\nWARNING messages were found.  Check warnings.txt for names.\n";
	print "There are $nWarnMsg lines of text in the file warnings.txt.\n";
	$lgWarnMsg = 1;
}
else
{
	print " good, none found.\n";
	$lgWarnMsg = 0;
}
	
if( -s "close.txt" )
{
	print "\nSome sims are close to a monitor.  Check close.txt for a list of names.\n";
	print "This is for information only, not a problem.\n";
}

print( "\nLooking for DEBUG string or TestCode call:");
if( -s "debug.txt" )
{
	print "\nDEBUG string or TestCode calls were found.  Check debug.txt for names.\n";
	print "There are $nDebug lines of text in the file debug.txt.\n";
	$lgDebugPrint = 1;
}
else
{
	print " good, none found.\n";
	$lgDebugPrint = 0;
}
	
print "$testgrid";

# check whether any sims were skipped
print "\nLooking for skipped sims:";
if( $nSkip > 0 )
{
	print " $nSkip sims were skipped - a list is in skip.txt\n" ;
}
else
{
	print " good: no sims were skipped.\n";
}

printf( "\n%i sims are present", $nMod );
if( $nSkip > 0 )
{
	printf( ", %i were skipped", $nSkip );
}
if( $nCrash > 0 )
{
	printf( ", %i crashed", $nCrash );
}
if( $nFail > 0 )
{
	printf( ", %i failed", $nFail );
}
if( $nNoAssert > 0 )
{
	printf( ", %i have no monitors", $nNoAssert );
}
if( $nWarn > 0 )
{
	printf( ", %i have warnings", $nWarn );
}
if( $nBotch > 0 )
{
	printf( ", %i botched", $nBotch );
}
if( $nOK > 0 )
{
	printf( ", %i are OK", $nOK );
}
print ".\n";

if( $lgMinorProblems )
{
	print "Minor problems were found, this is normal, list is in minor.txt\n";
}

if( $lgDebugPrint )
{
	print "DEBUG prints or TestCode calls were found, list is in debug.txt\n";
}

if( $lgWarnMsg )
{
	print "WARNING messages were found, list is in warnings.txt\n";
}

if( $lgSeriousProblems )
{
	print "Serious problems were found. See serious.txt for a summary\n";
}
else
{
	if( $nOK > 0 )
	{
		print "The test suite was successful; you have a valid executable.\n";
	}
}
