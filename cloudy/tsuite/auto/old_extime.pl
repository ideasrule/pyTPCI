#!/usr/bin/perl -w

# this script is obsolete - use the ir_extime.pl instead 

#**********************************************************************#
#*  Perl script to check the current Execution time of a model with   *#
#*  it's previous execution time. If difference is above a threshold  *#
#*  sends a mail stating so. It also creates a logs/ directory and    *#
#*  keeps a log of the execution time of the models, every time the   *#
#*  models are executed.                                              *#
#*  								      *#
#*  Program written by Geetashree Chakravorty, Graduate Student       *#
#*  Computer Science Department, University of Kentucky		      *#
#*  in the year 2002 for Dr. Gary Ferland			      *#
#**********************************************************************#


# Path set to use Perl mail utility SENDMAIL
#$sendmailpath="/usr/sbin/sendmail";
#open(SENDMAIL, "|$sendmailpath -t");  
 
#$from="gary\@pa.uky.edu";
#$to="gary\@pa.uky.edu"; 
#$subject="Incease in execution time";

# Common directory
$curr_dir="c://projects//cloudy//trunk//tsuite";

# Directory to store all the log files
$lgdir='/oldlogs/';

if( !chdir( $curr_dir ) )
{
   printf(" invalid directory\n");
   exit;
}

# If logs/ directory not present then create one
if(!-e $lgdir)
{
  system "mkdir $lgdir";
}

if( !chdir( "$curr_dir/auto" ) )
{
   printf(" invalid directory\n");
   exit;
}

$tlog='extime.log';
$efile='mailinc.log';
$thold=5;

open(LGFILE,">$tlog");
open(ELOG,">$efile");

print LGFILE "Records of execution times of models if .bak file exists";
print LGFILE "\n********************************************************"; 
print ELOG "List of models with increase in execution time\n\n";
print ELOG "Model \t %increase in execution time\n"; 

# Scans the .out files
while (defined($ofile=glob("*.out")))
{
   $bfile=$ofile;
   $bfile=~s/\.out/\.bak/;    # Substitutes '.out' by '.bak' in $bfile 

# Creates a '.log' file having same name as .out. Saves .log file in 
# directory /logs
   open(OUTFILE,"$ofile");
   $lfile=$ofile;
   $lfile=~s/\.out/\.log/;    
   $lfile=~s/$lfile/$curr_dir\/$lgdir\/$lfile/;

# Opening the .log files in append mode
   open(LOGFILE,">>$lfile");

# get current time and date
  ($sec,$min,$hour,$mday,$mon,$year) = localtime(time);
   $mon=$mon+1; #Since month returned is 1 less than actual
   $year=~s/^\d//;
   print LOGFILE "Date: $mon/$mday/$year\tTime: $hour:$min:$sec";

# Checks and computes if .bak file exists
   if(-e $bfile)
   {
      open(BKFILE,"$bfile");
      $mo=$ofile;
      $mo=~s/\.out//;        # removes '.out'
      print LGFILE "\nModel         :  $mo\n";

# Gets execution time from .bak file. Extracts the time by ignoring 
# everything in the '.out' files except the time after the string
# 'ExecTime' 
      while ($bk = <BKFILE>)
      {
         if ($bk=~/ExecTime/) 
         {
            $prev=$bk;
            $prev=~s/.*ExecTime//gi;  # gives only the execution time
            print LGFILE "Previous time : $prev";
          }
      }

# Gets execution time from .out file
      while ($new = <OUTFILE>)
      {                                                    
         if ($new=~/ExecTime/)
         {
            $curr=$new;
            $curr=~s/.*ExecTime//gi;  # gives only the execution time
            print LGFILE "New time      : $curr";
         }
      }

# Calculates difference between the execution times 
      $diff=$curr - $prev;
      print LGFILE "Difference    : $diff\n";

# Checking if present execution time is greater than previous execution time
      if ($diff !~/^-/) 
      {
         $ptage=($diff / $prev) * 100;
         if ($ptage > $thold) 
         {
            print LGFILE "% increase in time beyond threshold for model $mo  : $ptage \n";
            print ELOG " $mo\t   $ptage\n"; 
         }
      }
    }
# Prints the execution time to the .log file
    print LOGFILE "\t Exectime: $curr";
    close(LOGFILE);
}

 close(BKFILE);
 close(OUTFILE);
 close(LGFILE);
 close(ELOG);
 
# Checks if there is any model which took more time than the threshold.
# If yes than send a mail.
if (-e $efile){

# For Windows machines
   system("c:\\u\\blat\\blat.exe $efile -t gary\@pa.uky.edu -s \"Increase in exe time found in automatic testing\" " );

# For Unix machines
#   open(ELOG,"$efile");
#   print SENDMAIL "Subject: $subject\n";
#   print SENDMAIL "From: $from\n";
#   print SENDMAIL "To: $to\n\n";
#   while (<ELOG>){
#         print SENDMAIL "$_";
#      }
#   close(ELOG);

}
# End of program 
