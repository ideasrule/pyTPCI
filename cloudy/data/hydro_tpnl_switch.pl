#!/usr/bin/env perl
#
# Switch to the dataset of choice for the hydrogen Einstein A values.
#
# The default files are:
# 	hydro_tpn_n100.dat
#	hydro_tpnl_n100.dat
# and include values up to principal quantum number n=100.
# The user has the option of using the dataset that extends up to n=2000:
# 	hydro_tpn_n2000.dat.xz
#	hydro_tpnl_n2000.dat.xz
# Either is accessed through the symbolic links
# 	hydro_tpn.dat
#	hydro_tpnl.dat.
# This script points the links to the dataset of choice, decompressing the
# zipped files, if necessary.
#
# The only input on the command-line is the number 100 or 2000.
# 
# Created: Oct 19, 2020
# Author: Chatzikos, Marios
#
# Updated: Oct 21, 2020
# Author: Chatzikos, Marios
# Comment: Switched from bzip2 compression to xz, for a gain of a factor of 2
#	   in file size.  Now the zipped hydro_tpnl_n2000.dat file takes up
#	   about 40 MB (90 MB with bzip2, 450 MB uncompressed).
#

use strict;
use warnings;


our @links = ( "hydro_tpn.dat", "hydro_tpnl.dat" );

our @n100_files = ( "hydro_tpn_n100.dat", "hydro_tpnl_n100.dat" );
our @n2000_files = ( "hydro_tpn_n2000.dat", "hydro_tpnl_n2000.dat" );



sub read_command_line
{
	my $usage = "Usage: $0 <n>\n"
		  . "where:\n"
		  . "\t<n>: maximum principal quantum number; 100 or 2000\n";

	die $usage	if not @ARGV;

	my $n = shift @ARGV;

	die "Only 100 and 2000 accepted.  Given number not supported: '$n'\n"
		if( $n ne "100" and $n ne "2000" );

	return $n;
}

sub unzip_n2000
{
	for my $file ( @n2000_files )
	{
		my $xz = $file . ".xz";
		if( -e $xz )
		{
			print "Unzipping $xz...\t";

			die "Error: $xz did not successfully unzip\n"
			 if( system( "unxz $nz" ) != 0 );
			die "Error: $file not found!\n"
			 if( not -e $file );

			print "done!\n";
		}
	}
}

sub reset_links
{
	my $n = shift;

	for my $link ( @links )
	{
		unlink( $link ) if( -l $link );

		my $file = $link;
		$file =~ s/\.dat$/_n$n.dat/;

		die "Error: Target file $file does not exist\n"
		 if( not -e $file );

		warn "Warning: Failed to link $file to $link\n"
		  if( system( "ln -s $file $link" ) != 0 );
	}
}

sub main
{
	my $n = &read_command_line();
	&unzip_n2000()	if( $n eq "2000" );
	&reset_links( $n );
	return;
}

&main();
