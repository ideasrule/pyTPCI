#!/usr/bin/perl
$command = $ARGV[0];
$oldnam = $ARGV[1];
$newnam = $ARGV[2];
$remnam = $ARGV[3];
$outfile = $ARGV[4];
$temp_file = "$outfile.tmp";

# make sure that double quotes are escaped when passed to the shell
$command =~ s/"/\\"/g;

system( "$command > $temp_file" ) && &error_exit();
# this produces the dependencies for the object file
system( "cat $temp_file | sed \"s|$oldnam|$newnam|\" > $outfile" ) && &error_exit();
# this creates rules without recipes or prerequisites for each of the dependencies
# this helps if the dependency file is out of date and contains a dependency
# that no longer exists: this will now be ignored giving make a fighting chance...
# see http://scottmcpeak.com/autodepend/autodepend.html for more details.
system( "cat $temp_file | sed -e \'s/.*://\' -e \'s/\\\\\$//\' | fmt -1 | grep -v $remnam | sed -e 's/^ *//' -e 's/\$/:/' >> $outfile" ) && &error_exit();
unlink $temp_file;

sub error_exit {
    unlink $temp_file, $outfile;
    exit 1;
}
