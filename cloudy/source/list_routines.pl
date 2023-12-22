#!/usr/bin/perl

#****************************************************************************#
#* This program takes input from all the '*.c' files under the directory    *#
#* '/source/' and gives the filename, subroutines used and description of   *#
#* the subroutines into an output data file 'sortsrc.list' and html file    *#
#* sorttable.html.							    *#
#*  								            *#
#*  Program written by Geetashree Chakravorty, Graduate Student       	    *#
#*  Computer Science Department, University of Kentucky		      	    *#
#*  in the year 2002 for Dr. Gary Ferland. Modified: August,2003            *#
#****************************************************************************#

#Set the temporary file 'tempsrc.tmp' to a variable
$tsrc='tempsrc.tmp';

#Set the final file 'finalsrc.tmp' to a variable
$final='finalsrc.tmp';

#Setting name of 'list_routines.txt' to variable
$stsrc='list_routines.txt';

#Setting name of 'list_routines.html' to variable
$sttab='list_routines.html';

#open 'tempsrc.tmp' file
open(TSRC,">$tsrc"); 

#Opening and reading the .c files under Directory '/source/'	 	 
while(defined($infile=glob("*.cpp")))
{ 
  #assigning a handle to the file as it is opened
  open(FNAME,"$infile");
  
  while(<FNAME>)
  {
    # Removes the first two lines with copyright information present in the source c files
    if($_=~/^\/\* This file is part of Cloudy/ || $_=~/license.txt/)
    {
      $_=~s/.*//;
    }
    elsif($_=~/\#/)
    {
      last;
    }
    else
    {       #printing the output in a "filename(tab)reference" format 
      $_=~s/\/\*\s?//;          # Remove /*
      $_=~s/\t/ /g;
      $_=~s/^(\w+)(\s|\,)/$1\t/gi;	# setting a tab between the subroutines and the rest of the line
      $_=~s/\t\,/\t/; 		# Replace '(tab),' by (tab)
      $_=~s/\*\///;		# Remove */
      print TSRC "$infile\t$_";
   }
  }
}

#close the two files
close(FNAME);
close(TSRC);

#opening the temporary and final files and formatting the temporary 
#file to have continuing lines into 1 line and setting the result to
#final file 'finalsrc.tmp'
open(FSRC,">$final"); 
open(TSRC,"$tsrc"); 

while(<TSRC>)
{
  $_=~s/\n$//;
  if($_=~/^\w+\.\S\s*$/)
  {
    $_=~s/^\w+\.\S\s*$//;  #Removing the lines which has only the C file name
  }
  else
  {
    if($_=~/\*/)
    {
      $_=~s/.*\s*\*/ /;
      print FSRC "$_";
    }
    else
    {
      print FSRC "\n$_";
    }
  }
}

#closing both the opened files.
close(TSRC);
close(FSRC);
unlink($tsrc);

#Getting the sorted data in 'finalsrc.tmp' into a new file 'sortsrc.list'
system "sort -f -k 1 -k 2 <$final> $stsrc";
unlink($final);

#Creating tables in HTML with the data in 'sortsrc.list' in tables
open(SSRC,"$stsrc"); # opens input file 'sortsrc.list'
open(STAB,">$sttab"); #opens output HTML file 'sorttable.html'

print STAB "<TABLE border=\"1\" cellpadding=\"5\" cellspacing=\"2\" colspan=\"5\" scope=\"colgroup\"
	summary=\"This is a table of the list of C programs used in '/source/' 
	and the various subroutines used within the programs\">
	<CAPTION> Table from data in 'sortsrc.list'</CAPTION>
    	<TR>
		<TH scope=\"col\">Filename</TH>
		<TH scope=\"col\">Subroutine</TH>
		<TH scope=\"col\">Description</TH> 
	</TR>";
while(<SSRC>)
{
	if($_=~/^\n$/)
	{
	  $_=~s/^\n$//;
        }
	else{
	$_=~s/\t/\<\/TD\>\<TD\>/g;  #formatting the data into rows and columns
	print STAB "<TR><TD>$_</TD></TR>";}
}
print STAB "</TABLE>";

#closing both the final files 
close(SSRC);
close(STAB);

print "The output file is: $stsrc.\n";
print "The HTML file is: $sttab.\n\n"; 



