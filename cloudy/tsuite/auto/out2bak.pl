#!/usr/bin/perl
# copy all outpout *.out to bak
  
# bring in the perl module that includes copy
use File::Copy;

$nMod = 0;
while(defined($output= glob("*.out")) )
{
           # previous file is named .bak
           $previous = $output;
           $previous =~ s/\.out//gi;
           $bak = "$previous".".bak";
           # simply copy the source here
           copy( $output  , $bak );
          ++$nMod
}
printf( "\n %i out files copied to bak.\n ",$nMod);

