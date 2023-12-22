#!/usr/bin/perl -w


# edit all input models that do not match a certain string  

$nMod = 0;
while(defined($input= glob("*.in")) )
{
        system "grep 'punch dr'  $input >checkend.txt";

#	if non-zero length, model had botched asserts
#       -s option returns number of bytes if file exists and has finite size
	if( (-s "checkend.txt")==0 )
	{
           ++$nMod;

#         actually view the output
          system "gvim $input";
	}
}
printf( "\n %i models had the string, and were edited.\n ",$nMod);

