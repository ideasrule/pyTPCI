#!/usr/bin/perl
# program to find out the C++ files using the header(.h) files.

# this will be list of headers
$hfiles='headers.txt';

# List of header files  and C++ files using them  
$result='listfiles.list'; 

# List of unused header files
$unused='headers_unused.txt';

print "A list of header files will be placed in ", $hfiles , "\n";
print "All source files using each header will be placed in ", $result , "\n";
print "Unused header files will be placed in ", $unused , "\n";
print "This will take a while...\n\n";

#Start of main
#system "ls *.h > $hfiles";
open(ioHEADERS,">$hfiles"); #ioHEADERS-Output file
open(LFILE,">$result");
open(UFILE,">$unused");

while(defined($header=glob("*.h")))
{
  print ioHEADERS "$header\n";
}
close(ioHEADERS);

open(IHFILE,"$hfiles"); #IHFILE-Input file
while(<IHFILE>)
{
  $header=$_;
  print LFILE "$header";
  $header=~s/\n//;
  $nuses = 0;
  while(defined($input=glob("*.h *.cpp")))  #Scanning through the C++ files
  {
     next if $header eq $input;
     open(CFILE,"$input");
     while(<CFILE>)
     {
       if($_=~/$header/)
       {
         print LFILE "\t$input";
         $nuses++;
       }     
     }
  }
  print LFILE "\n\n"; 
  print UFILE "$header\n"
    if $nuses == 0;
}
close(IHFILE);
close(LFILE);
close(UFILE);
#End of program
