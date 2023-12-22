#!/usr/bin/perl -w
#
# look for models that did not end.
# these will be recomputed
#

# check all output files to see whether they contain "Cloudy ends:"
# string and write corresponding line into checkend.txt file.
# if size of file $o_dir/checkend.txt is zero, then program crashed,
# and we want to rerun it.  

$nMod = 0;
while ( defined( $input = glob("*.in") ) )
{
	$output = $input;
	$output =~ s/\.in/.out/gi;
	#system "grep 'Cloudy ends'  $output > checkend.txt";
	system "grep 'ChkMonitor'  $output > checkend.txt";

#	if zero length, code did  not end
	if ( -z "checkend.txt" )
	{
		++$nMod;

		system "vs $output";
	}
}
printf( "\n %i models had crashed, and were viewed.\n ",$nMod);

