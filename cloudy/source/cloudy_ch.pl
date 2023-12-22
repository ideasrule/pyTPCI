#!/usr/bin/perl
# Program to swap the position of indexes of a two-dimentional array
#abund.xIonFracs[b][a], i.e., abund.xIonFracs[b][a]<->abund.xIonFracs[a][b] 

#Store the changes initially into a temporary directory.
$ndir='cloudy_ch/';
if(!-e $ndir)
{
  system "mkdir $ndir";
}
else
{
 system "rm -f $ndir/*.*";
}

#Reading the files and swapping the indexes of array abund.xIonFracs
while(defined($input=glob("*.c *.h")))
{ 
   if((system "grep 'abund.xIonFracs' $input")==0)
   {
     $input=~s/\n//;
     $out=$input;
     $out=~s/$input/$ndir$input/;
     open(OUTFILE,">$out");
     open(INFILE,"$input");
     while(<INFILE>)
     {
       if($_=~/abund.xIonFracs\[\S*\]\[\S*\]/)
       {
         $_=~s/abund.xIonFracs\[(\w*|\w*\.\w*\[\w*\]\W\w*|\w*\W\w*)\]\[(\w*|\w*\.\w*\[\w*\]\W\w*|\w*\W\w*)\]/abund.xIonFracs\[$2\]\[$1\]/gi;
         print OUTFILE "$_"; 
       }

       else
       {
	print OUTFILE "$_"; 
       }
     }
 }
#closing both the input and output files.
close(OUTFILE);
close(INFILE);
}

# Copying the files of temporary directory into original files and
# then removing the temporary directory 
#system " cd cloudy_ch; cp *.* ..; rm *.* .*";
#system "rmdir cloudy_ch";

