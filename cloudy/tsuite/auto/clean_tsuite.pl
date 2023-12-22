#!/usr/bin/perl

# This perl script cleans out the test suite after a test run, it will
# only leave files that are part of the official distribution and the
# cloudy executable (with a name ending in .exe) if that was present!
#
# Peter van Hoof

while( defined( $input = glob("*") ) ) {
	@ll = split( /\./, "$input" );
	if( $#ll != 1
	    or ( $ll[1] ne "in"
		 and $ll[1] ne "pl"
		 and $ll[1] ne "pm"
		 and $ll[1] ne "htm"
		 and $ll[1] ne "md"
		 and $ll[1] ne "jpg"
		 and $ll[1] ne "dat"
		 and $ll[1] ne "txt"
		 and $ll[1] ne "exe"
		 and $ll[1] ne "pdf"
		 and $ll[1] ne "vsz" ) ) {
		print "deleting $input...\n";
		unlink "$input";
	}
	elsif( $ll[0] =~ /^grid[0-9]{9}/ ) {
		print "deleting $input...\n";
		unlink "$input";
	}
}
unlink "mpi_batch_jobs.txt";
unlink "checkend.txt";
unlink "close.txt";
unlink "crashed.txt";
unlink "debug.txt";
unlink "minor.txt";
unlink "serious.txt";
unlink "skip.txt";
unlink "wvlng.txt";
unlink "warnings.txt";
unlink "optimal.in";
