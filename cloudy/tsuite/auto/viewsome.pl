#!/usr/bin/perl -w
#This progarm reads in the names of all the input "*.in" files 
#in runsome.dat copies its models and saves them as "*.out"

#Assigning a variable to 'runsome.dat'
$rsome='viewsome.dat';

#opening 'runsome.dat'
open(RSOM,"$rsome");

#The command to be executed
$exe = "/home59/home37/gary/cloudy/trunk/alpha/cloudy_alpha.exe";
# this is the path on my XP box
#$exe = "c:/Projects/Cloudy/trunk/debug/trunk.exe";
$exe = "c:/projects/Cloudy/trunk/release/trunk.exe";

print "\n";

#Reading in the input file names that are to be executed
while($r=<RSOM>)
{ 
  system "vs $r";   #copy the contents of the .in files to its counterpart .out files
}
close(RSOM);
 

