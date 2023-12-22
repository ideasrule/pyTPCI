#!/usr/bin/perl
#This program takes input from a number of ".cpp" and ".h" files
#The output it generates is a series of references which is stored in 
#'doc_atomic_data_refer.txt' and a series of old references stored in 'doc_atomic_data_refer_old.txt'

#'tempfile.tmp' is an intermediate fiel
$tfile='tempfile.tmp';

#'atomicdata.list' is the final formatted output file
$atdat='doc_atomic_data_refer.txt';

#'oldsatomicdata.list' is the file which contains old references
$olatdat='doc_atomic_data_refer_old.txt';

#opening the temporary file 'tempfile.tmp' and assigning a handle to it. 
open(TFILE,">$tfile"); 

#Opening and reading the .c and .h files	 	 
while(defined($infiles=glob("*.cpp *.h")))
{ 
  #assigning a handle to the file as it is opened
  open(FNAME,"$infiles");
  while(<FNAME>)
  {	
    if($_=~/>>refer/)
    {
      #printing the output in a "filename(tab)reference" format 
      print TFILE "$infiles\t$_";
    }
   }
}

# get all references in the data files
$in_dir = "c:/projects/cloudy/trunk/data/";
#Opening and reading the .dat files in the data dir
while(defined($infiles=glob("$in_dir*.dat")))
{ 
  #assigning a handle to the file as it is opened
  open(FNAME,"$infiles");
  while(<FNAME>)
  {	
    if($_=~/>>refer/)
    {
      #printing the output in a "filename(tab)reference" format 
      print TFILE "$infiles\t$_";
    }
  }
}
#Opening and reading refractive index files in the data dir
while(defined($infiles=glob("$in_dir*.rfi")))
{ 
  #assigning a handle to the file as it is opened
  open(FNAME,"$infiles");
  while(<FNAME>)
  {	
    if($_=~/>>refer/)
    {
      #printing the output in a "filename(tab)reference" format 
      print TFILE "$infiles\t$_";
    }
  }
}

# get all references in the test case files
$in_dir = "c:/projects/cloudy/trunk/tsuite/auto/";
#Opening and reading the .dat files in the data dir
while(defined($infiles=glob("$in_dir*.in")))
{ 
  #assigning a handle to the file as it is opened
  open(FNAME,"$infiles");
  while(<FNAME>)
  {	
    if($_=~/>>refer/)
    {
      #printing the output in a "filename(tab)reference" format 
      print TFILE "$infiles\t$_";
    }
  }
}

#close the two files
close(FNAME);
close(TFILE);

#Opening the 'tempfile.tmp','atomicdata.list' and 'oldatomicdata.list'
#Takes input from 'tempfile.tmp' and writes the 
#new references to 'atomicdata.list' and old references to 'oldatomicdata.list'
open(TFILE,"$tfile");
open(ATDAT,">$atdat");
open(OLDAT,">$olatdat");

while(<TFILE>)
{	
  $_=~s/\n$//;
  if($_=~/refercon/)
  {
    $_=~s/.*\s*.*\W*>>refercon\s/ /;
    $_=~s/\W*$//;
    print ATDAT "$_";      #writing lines to 'atomicdata.list'
  }

  elsif($_=~/referold/)
  {
    $_=~s/\t*\W*\>\>referold(\t|\s*)/\t/;
    $_=~s/\W*$//;
    print OLDAT "\n$_";	   #writing lines to 'oldatomicdata.list'
  }

  else
  {
  $_=~s/\t*\W*\>\>refer(\t|\s*)/\t/;
  $_=~s/\W*$//;
  print ATDAT "\n$_";      #writing lines to 'atomicdata.list'
  }
}
close(TFILE);
close(ATDAT);
close(OLDAT);
unlink($tfile);    #deletes the 'tempfile.tmp'
 

print "\nI created files \n";
print "$atdat containing new references\n";
print "$olatdat containing old references.\n\n";









