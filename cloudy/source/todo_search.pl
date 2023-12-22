#!/usr/bin/perl
#****************************************************************************#
#* The program uses the C source files as the input scans it and retrieves  *#
#* all instances of TODO and the comments following them		    *#
#*		 Then it formats the arguments of the functions into a      *#
#* table having "label for ion", "wavelength", "type of entry for line"     *#
#* and "preceeding comments" ( in order) as the columns.                    *#
#*                                                                          *#
#*  Program written by Geetashree Chakravorty, Graduate Student             *#
#*  Computer Science Department, University of Kentucky                     *#
#*  in the year 2003 for Dr. Gary Ferland.                                  *#
#****************************************************************************#

use Cwd;        # To find the working directory path

$currpath = getcwd;
$tempout = "outputtodo.txt"; # Temporary output file
$finalout = "todo_table.txt";    # Final output file

#Setting name of 'tabledtodo.html' to variable
$finalhtm = "todo_table.htm";

$sourcedir = "../source/";   # Directory having the source C files. This needs to be customized
$prev ="";

# Move to current path
if ( !chdir ($currpath))
{
   print (" invalid directory\n");
   exit;
}

# Opening the temporary file
open(TEMPFILE, "> $tempout ");
# this line is needed so that the "file" utility correctly identifies the output as ascii text.
print TEMPFILE "#######################################################\n";
while (defined ($sourcefile = glob ("$sourcedir*.h $sourcedir*.cpp"))){
  open (SFILE,"< $sourcefile");
  $flag = 0;
  $countline = 0;
  
  # If source files are in a different directory
  if($sourcefile =~/\//){  
  $sourcefile =~s/(.*\/)*(\w*\.\w)/$2/;	#Extracting the source file name
  }
  while($line = <SFILE>){
    $line =~s/^\s*//;
    $countline ++;
    
    # If line has TODO
    if($line =~/TODO/i && $flag == 0){
      # Retrieving the priority 
      $line_no = $countline;
      $priority = $line;
      $priority =~s/(.*TODO\s*)(\d*)(\s.*)/$2/i;
      $priority =~s/\n$//;

      $comment = $line;
      if( $priority ne '' )
      {
      	$comment =~s/(.*TODO\s+$priority)(.*)/$2/i;
      }
      else
      {
      	$comment =~s/(.*TODO\s+)(.*)/$2/i;
      }
      $comment =~s/^\s*-\s+//;
      $comment =~s/\*\///;
      $comment =~s/\/\*\*//;
      $comment =~s/\t*//gi;
      $prev = join '',$prev,$comment;  
      
      # Checks if TODO comments extend more than one line
      if($line !~/\*\// and $line !~ m/\/\// ){
        $flag = 1;
      }
      else{ # Print into output temporary file
        print TEMPFILE "$priority\t$sourcefile\t$line_no\t$prev";
	$prev="";      
      }
    }
    
    # if TODO comment is more than one line
    elsif($flag == 1){
      $line =~s/^(\s*|\t)\*//;
      $line =~s/\t/ /gi;
      $prev = join '',$prev,$line;
      if($line =~/.*\*\// || $line =~/\//){
        $flag = 0;
	$comment = $prev;
	$comment =~s/\t*//gi;
	$comment =~s/\n//;
	$comment =~s/\*\///;
	$comment =~s/\///;
	# Print into output temporary file
	print TEMPFILE "$priority\t$sourcefile\t$line_no\t$comment";
	$prev="";
      }   
      else{
        $prev =~s/\n$/ /;
      }
    }
  } 
  close(SFILE);  
}
close(TEMPFILE);

# Sorting the tabled data from the temporary file to the final output file
system "sort -f -k 1 <$tempout> $finalout";
unlink($tempout);	# Delete temporary file

#opens output HTML file 'tabledtodo.html'
open(FHTM,">$finalhtm"); 
print FHTM "<HTML><BODY>
<TABLE border=\"1\" cellpadding=\"5\" cellspacing=\"2\" colspan=\"5\" scope=\"colgroup\">
        <CAPTION> This is a table of the list of TODO used in CLOUDY source files
        and the comments following them. \n Table from data in 'tabledtodo.txt'</CAPTION>
        <TR>
                <TH scope=\"col\">Priority</TH>
		<TH scope=\"col\">Filename</TH>
                <TH scope=\"col\">Line number</TH>
                <TH scope=\"col\">Comment</TH>
        </TR>";

# Copying the data from tabledtodo.txt to tabledtodo.html
open(FINAL, "< $finalout");
while(<FINAL>){
	next if $_ =~ m/^#+$/;
	$_=~s/\t/\<\/TD\>\<TD\>/g;  #formatting the data into rows and columns
	print FHTM "<TR><TD>$_</TD></TR>"
}
print FHTM "</TABLE></BODY></HTML>\n";
close(FINAL);
close(FHTM);
print "\ncreated output files ", $finalout ," and ", $finalhtm , "\n\n";

