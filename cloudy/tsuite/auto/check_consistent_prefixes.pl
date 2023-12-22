#!/usr/bin/perl -w

# check each save command in each *.in file to see that the prefix is consistent.

# loop over all input scripts
foreach $input (<*.in>)
{
	# form name of output files
	$outname = $input;
	$outname =~ s/\.in$//;
	system "grep -H '^save ' $input | grep -v '\"$outname'";
}

