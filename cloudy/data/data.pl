#!/usr/bin/perl

#This program takes input  data files "*.dat" .
#The output it generates is a series of references which is stored in 
#'finaldata.txt' .

#'tempdat.tmp' is an intermediate field
$tdfile='tempdat.tmp';

#finaldata.txt is the final formatted file
$fdata='finaldata.txt';

#opening the temporary file 'tempdat.tmp' and assigning a handle to it. 
open(TDFILE,">$tdfile"); 

#Opening and reading the .c and .h files	 	 
while(defined($indata=glob("*.dat")))
{ 
  #assigning a handle to the file as it is opened
  open(DATAF,"$indata");
  while(<DATAF>)
  {	
    if($_=~/>>refer/)
    {
      #printing the output in a "filename(tab)reference" format 
      $_=~s/\W*//;
      $_=~s/\W*$//;	
      print TDFILE "$indata\t$_\n";
    }
   }
}

#close the two files
close(DATAF);
close(TDFILE);

open(TDFILE,"$tdfile");
open(FDATA,">$fdata");

while(<TDFILE>)
{
  $_=~s/\n$//;
  if($_=~/refercon/)
  {
    $_=~s/.*\s*.*refercon\s/ /;
    print FDATA "$_";      #writing lines to 'final.txt'
  }
  else
  {
    $_=~s/\t*refer(\t|\s*)/\t/;
    print FDATA "\n$_";
  }
}

#closing the two files
close(FDATA);
close(TDFILE);
unlink($tdfile);	#delete the temporary file 'tempdat.tmp

print "The file created is:\n";
print "$fdata.\n\n";












