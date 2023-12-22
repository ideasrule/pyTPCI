#!/usr/bin/perl -w

#*****************************************************************************#
#*  This progarm reads all the input "*.in" removes the 'monitor's within    *#
#*  the file and saves it in another directory 'tsuite'                      *#
#*                                                                           *#
#*  Program written by Geetashree Chakravorty, Graduate Student              *#
#*  Computer Science Department, University of Kentucky                      *#
#*  in the year 2001 for Dr. Gary Ferland                                    *#
#*****************************************************************************#

#assign directory name to variable
$ndir='tsuite/';

#Checks if directory exists. If 'yes' then delete all files in directory else 
#create a new directory
if(!-e $ndir)
{
  system "mkdir $ndir";
}
else
{
  system "rm -f $ndir/*.in";
}

#start of source program
#open *.in files
while(defined($input=glob("*.in")))
{ 
#Reading in the input file names that are to be executed

  $input=~s/\n//; 	  #removing any newline characters
  $out=$input;	  #giving name of output file same as the input file.
  $out=~s/$input/$ndir$input/;
  open(OUTFILE,">$out"); 
  open(INFILE, "$input");
  while(<INFILE>)
  {
    if ($_!~/monitor/)     #if not 'monitor'
    {
      print OUTFILE "$_";  # print to files in /tsuite
    } 
  }
   close(OUTFILE);
   close(INFILE);
}

print "The formatted files are in the directory $ndir\n";
#End of program




