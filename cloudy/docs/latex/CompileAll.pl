#!/usr/bin/perl -w
# compile Quick Start Guide and all three parts of Hazy with pdflatex
# do multiple times to resolve bibtex entries
#
# 2015-Oct-14	allow individual directories on the command line
# 2015-Oct-21	automate number of recompilations
#
sub need_another_run
{
	my $dir = shift;
	open FILE, "<${dir}.log";
	my @lines = <FILE>;
	close FILE;
	return 1
		if( grep( /Rerun to get cross-references right/, @lines )
		 or grep( /Rerun to get citations correct/, @lines ) );
}
sub run_tex
{
	my $dir = shift;
	substr($dir, -1, 1, "")	if( $dir =~ m/\/$/ );
	chdir( "$dir" ) or die "invalid directory $dir\n";
	system( "pdflatex $dir" ) == 0 or die;
	system( "bibtex $dir" ) == 0 or die;
	my $counter = 0;
	do
	{
		system( "pdflatex $dir" ) == 0 or die;
		++$counter;
	} while( &need_another_run( $dir ) and $counter < 5 );
	chdir( ".." );
}
sub run_tex_once
{
	my $dir = shift;
	substr($dir, -1, 1, "")	if( $dir =~ m/\/$/ );
	chdir( "$dir" ) or die "invalid directory $dir\n";
	system( "pdflatex $dir" ) == 0 or die;
	chdir( ".." );
}

my @dirs_def = qw/ hazy1 hazy2 hazy3 QuickStart /;

my @dirs;
push(@dirs, @ARGV)	if( @ARGV );
push(@dirs, @dirs_def)	if( not @dirs );

foreach my $dir ( @dirs ) { &run_tex( $dir ); }

#
# Fix cross-references between Hazy 1 & Hazy 2
#
if( join( ' ', @dirs ) =~ m/hazy1/ and join( ' ', @dirs ) =~ m/hazy2/ )
{
	&run_tex_once( 'hazy1' );
	&run_tex_once( 'hazy2' );
}
