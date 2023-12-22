#!/usr/bin/perl -w


# view all models with botched asserts  

$nMod = 0;
while(defined($output= glob("*.out")) )
{
        system "grep \"ABORT ABORT\"  $output >checkend.txt";

#	if non-zero length, model had abort
#       -s option returns number of bytes if file exists and has finite size
	if(-s "checkend.txt")
	{
           ++$nMod;

#         actually view the output
          system "vs $output";
	}
}
printf( "\n %i models had aborts and were viewed.\n ",$nMod);

