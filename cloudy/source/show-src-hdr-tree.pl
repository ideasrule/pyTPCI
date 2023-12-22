#!/usr/bin/perl
#################################################################################
#  show-src-hdr-tree.pl::
#     Display the tree of header calls for all the source files in the directory,
#  or those specified in the command line.
#
#  Chatzikos, Dec 22, 2016
#################################################################################

use strict;
use warnings;


use lib "./";
use headers;


#################################################################################
#                               USER INPUT                                      #
#################################################################################
sub printUsage
{
	die "$0 [-h] [<source>]\n"
	 .      "where:\n"
	 .      "\t      -h:\t print this message and exit\n"
	 .      "\t<source>:\t list of source (.cpp) files to process [OPTIONAL]\n"
	 .      "\t         \t if not given, operate on *.cpp\n";
}

sub ConfirmInput
{
	my( $cpp ) = @_;

	if( @$cpp )
	{
		foreach my $this_cpp ( @$cpp )
		{
			#
			# This is impossible when called with command-line arguments
			#
		       die "Error: Undefined source file\n"
			 if( not defined( $this_cpp ) or $this_cpp eq "" );

			die "Error: Source file ('$this_cpp') does not exist or empty\n"
			 if( not -e $this_cpp or not -s $this_cpp );
		}
	}
}

sub Init
{
	my @argv = @_;

	my @cpp = @argv;

	&printUsage()   if( $cpp[0] eq "-h" );

	&ConfirmInput( \@cpp );

	@cpp = glob "*.cpp"     if( not @cpp );

	if( 0 )
	{
		print "cpp:\t @cpp\n";
		die;
	}

	return  ( @cpp );
}



#################################################################################
#                               SOURCE FILE PROCESSING				#
#################################################################################
sub process_source
{
	my $this_cpp = shift;

	die "Error: Undefined source file\n"
	 if( not defined( $this_cpp ) or $this_cpp eq "" );

	die "Error: Source file ($this_cpp) does not exist / empty\n"
	 if( not -e $this_cpp or not -s $this_cpp );

	my @contents = &headers::get_contents( $this_cpp );
	my @headers = &headers::get_includes_from_cpp( undef, @contents );

	&headers::prtSrcTree( $this_cpp, @headers );

	print "\n";
}


#################################################################################
#                               MAIN PROGRAM					#
#################################################################################
my( @cpp ) = &Init( @ARGV );

&headers::header_dependencies();

foreach my $this_cpp ( @cpp )
{
	&process_source( $this_cpp );
}
