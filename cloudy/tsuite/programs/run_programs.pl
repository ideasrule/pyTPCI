#!/usr/bin/perl
# This progarm reads in the names of all the program files 
# in the file run_programs.dat, compiles and links them, then runs the tests
#
# Some examples of the calling sequence for this script are:
#
# run_programs.pl
# run_programs.pl sys_icc
#
# in the first form it will use the binary in the source directory itself, in
# the second form you can supply the name of one of the sys_xxx subdirectories.

# location of source directory
$DirSource  = "../../source";

if( ! -r "$DirSource/Makefile" )
{
	die "PROBLEM could not find file $DirSource/Makefile\n";
}

# this is where the object files should be
$DirObject  = "$DirSource/$ARGV[0]";

if( ! -d "$DirObject" )
{
	die "PROBLEM could not find directory $DirObject\n";
}

# now determine the compiler and flags to use by polling make
$cxx = `cd $DirObject; make echo-cxx | egrep -v '^make' | tail -n 1`;
chomp( $cxx );
$cxxflags = `cd $DirObject; make echo-cxxflags | egrep -v '^make' | tail -n 1`;
chomp( $cxxflags );
# escape double quotes
$cxxflags =~ s/"/\\"/g;

# this file contains the names of the directories containing the programs
$rsome = 'run_programs.dat';
open( RSOM, "<$rsome" ) or die "PROBLEM could not open $rsome\n";

# read the program names that are to be executed, each is in a directory of the same name
while($r=<RSOM>)
{
	# check is beginning of line (^) is not a word
	if( $r !~ /^\W/ )
	{
		# remove trailing newline
		chomp( $r );
		chdir( $r ) or die "PROBLEM could not change to $r\n";
		$command = "$cxx $cxxflags $r.cpp -o $r.exe -I../$DirSource -L../$DirObject -lcloudy";
		# print "$command\n";
		print "compiling $r.cpp\n";
		$res1 = system( "$command" );
		if( $res1 != 0 )
		{
			print "PROBLEM compilation of $r.cpp failed.\n";
		}
		else
		{
			print "executing $r.exe\n";
			$res2 = system( "./$r.exe" );
			if( $res2 != 0 )
			{
				$res2 >>= 8;
				print "PROBLEM executing $r.exe failed, return code: $res2\n";
			}
		}
		chdir( ".." ) or die "PROBLEM could not return home\n";
	}
}

close(RSOM);


 

