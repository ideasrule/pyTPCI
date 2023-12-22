#!/usr/bin/perl -w

# print file name for all files corresponding to *.in that do not 
# start with a sharp sign, the default first character for header
# information in a punch file.  Run this to confirm that headers
# are produced properly in punch output

# loop over all input scripts
foreach $input (<*.in>)
{
	# form name of output files
	$outname = $input;
	$outname =~ s/\.in$//;
	foreach $output (<$outname.*>)
	{
		open IN, $output or die "Couldn't open $output";
		$line = <IN>;
		close IN;
		if ($line) {
				print "$output\n" if not $line =~ /^#/;
		}
  }
}

