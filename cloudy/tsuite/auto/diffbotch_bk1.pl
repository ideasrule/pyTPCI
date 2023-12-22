#!/usr/bin/perl -w


# view all models with botched asserts  

$nMod = 0;
while(defined($output= glob("*.out")) )
{
	#system "grep 'Cloudy ends'  $output >checkend.txt";
        system "grep Botched  $output >checkend.txt";
        system "grep 'W-'  $output >>checkend.txt";

#	if non-zero length, model had botched asserts
#       -s option returns number of bytes if file exists and has finite size
	if(-s "checkend.txt")
	{
           ++$nMod;
           # previous file is named .bak
           $previous = $output;
           $previous =~ s/\.out//gi;
           $bak = "$previous".".bk1";

#         actually view the output
          system "vsdiff $bak $output";
	}
}
printf( "\n %i models had botched asserts, and were viewed.\n ",$nMod);

