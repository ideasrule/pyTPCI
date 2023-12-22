#!/usr/bin/perl
# copy all outpout *.out to bak
  

# bring in the perl module that includes copy
$oldname="casea";
$newname="h_casea";

system( "sed -es/$oldname/$newname/ < $newname.in > test.txt ") ;
rename "test.txt" $newname.in ;

