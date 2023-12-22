#!/usr/bin/perl

#  This program takes  the '*.in' files, one at a time, 
#  creates a HTML file, 'doc_tsuite.htm' and prints the inputs 
#  and description into the HTML file. The names of the input 
#  files and the tasks they perform are printed in a text file 
#  'doc_tsuite.txt'.

#  Program written by Geetashree Chakravorty, Graduate Student
#  Computer Science Department, University of Kentucky
#  in the year 2001 for Dr. Gary Ferland


#Assigment of variables to the two output files
# Variable assignment for the output HTML file
$htests='doc_tsuite.htm';  
# Variable assignment of the output Text file
$tests='doc_tsuite.txt';    
$tfile='tempfile.tmp';

#open the output HTML file
open(THTML,">$htests");

#open the output text file
open(TTEXT,">$tests"); 

#two variables $bool1 and $bool2 has been assigned as 
#flags to determine whether to print the lines with 'monitor' and '##'
#initially $bool1 and $bool2 has been set to '0' which means the 
#printed result will not print the lines having 'monitor' and '##'. 

#To have 'monitor' and '##' back on we have to remove the comment sign 
#before the next two print commands and comment the default values set 
#and then set $bool1=<> and $bool2=<>.

# the conditions are as follow:
#	$bool1	$bool2
#	  0	  0   --->no monitor and ## (Default)
#	  0	  1   --->no monitor
#	  1       0   --->no ##
#	  1       1   --->all included 	
	 	 
#print "Do you want to have 'monitor' on (Enter 0 for No and 1 for Yes): \n";
$bool1=0;
#print "Do you want the '##' on (Enter 0 for No and 1 for Yes): \n";
$bool2=0;

print THTML "<html>\n
<head>\n
<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">\n
   <meta name=\"GENERATOR\" content=\"Mozilla/4.77 [en] (X11; U; Linux 2.4.3-12 i686) [Netscape]\">\n
   <title>Cloudy test suite</title>\n

</head>
<body>";
print THTML "<h4>  This HTML file was created by the program doc_tsuite.pl, doc_tsuite.txt contains a tab delimited list of the files. </h4>\n";

#start of source program
#open *.in files
while(defined($input=glob("*.in")))
{ 
    open(TEST,"$input");   #TEST is the handle to the input .in files

    #$flag is a flag whih is set to 1 when it encounters a blank space
    $flag=0;

    #prints the name of the file as header
    print THTML "<hr><kbd><h2>$input</h2>  ";
    #reading in each line of file, formatting it and printing the formated result in 'text2.html'

    while(<TEST>)
    {
	if($_=~/^title/i)
	{
	    $rtitle=$_;
	    $rtitle=~s/title\s//;  # Removes 'title' from line  
	    print THTML "</kbd><h2><I>$rtitle</I></h2><kbd>"; 
	    chomp( $rtitle );
	    $pinput = $input;
	    # no need for all the .in
	    $pinput =~s/\.in//gi;
	    # change _ to \_ to keep latex happy
	    $pinput =~s/_/\\_/gi;
	    $rtitle =~s/_/\\_/gi;
	    $rtitle =~s/\^/\\^/gi;
	    $rtitle =~s/&/\\&/gi;
	    print TTEXT "$pinput & $rtitle \\\\ \n";
	}
	if($flag==0)
	{
	    if($_!~/^\n*$/)  
	    {
		# printing files in a line to line basis
		if ($bool1==0 && $bool2==0)
		{	
		    if($_!~/##/)    # if not comments
		    {
			if($_!~/monitor/)   # if not monitors
			{
			    print THTML "<br>$_";
			}
		    }
		}
		elsif($bool1==0 && $bool2==1)
		{
		    if($_!~ /monitor/)    # if not monitors
		    {
			print THTML "<br>$_";
		    }
		}	
		elsif($bool1==1 && $bool2==0)
		{	
		    if($_!~/##/)       # if not comments 
		    {
			print THTML "<br>$_";
		    }
		}
		elsif($bool1==1 && $bool2==1)
		{	
		    print THTML "<br>$_";
		}
		else
		{
		    print "Error in setting values.\n"
		}
	    }
	    else
	    {
		print THTML "</kbd>";	
		$flag=1;
	    }
	}
	if($flag==1)
	{
	    if($_=~/^\n$/)
	    {
		print THTML "<p>";	# prints description in a paragraph format.
	    }
	    if($_=~/^(\n*Checks|\-)/)
	    {
		print THTML "<br>$_";   # prints 'Checks'
	    }
	    else
	    {
		print THTML "$_";
	    }
	}
    }
    if($flag==0)
    {
	print THTML "</kbd>";
    }
    print THTML "</p>";
    close(TEST);
}

print THTML "<hr>\n
</body>\n
</html>\n";

#closing both 'doc_tsuite.txt' and 'doc_tsuite.html'
close(TTEXT);
close(THTML);

print "\n$tests contains a list of the file names and titles.\n";
print "$htests has the formatted contents of each test case.\n";
















