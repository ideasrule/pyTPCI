#!/usr/bin/perl
  
# renames test suite files from $oldname to $newname
# also renames punch file names
# and log files

$oldname="limit_lte_hhe";
$newname="limit_lte_hhe_ste";

$nMod = 0;
while(defined($input= glob("$oldname.*")) )
{
           # previous file is named .bak
           $new = $input;
           $new =~ s/$oldname/$newname/gi;
           print "old: " , $input , " to " ,$new, "\n" ;
           rename $input , $new ;
          ++$nMod
}
printf( "\n %i files were renamed.\n ",$nMod);

# change the name of the log files - these live in a different dir
rename "../itrlogs/$oldname.log" , "../itrlogs/$newname.log";

# change name of file within the main input script
system( "sed -es/$oldname/$newname/ < $newname.in > test.txt ") ;
rename "$newname.in" , "$oldname.sav" ;
rename "test.txt" , "$newname.in" ;
print "names within ", $newname , " changed\n";


