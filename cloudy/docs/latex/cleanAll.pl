#!/usr/bin/perl -w
# Remove temporary files from the specified directories.
# By default operate on all parts of Hazy and the QuickStart
# guide.
#
# 2015-Oct-27
#
sub clean_dir
{
	my $dir = shift;
	foreach my $file ( glob "${dir}.*" )
	{
		unlink "$file" if( $file !~ /\.tex$/ );
	}
	foreach my $file ( glob "*.aux" )
	{
		unlink "$file";
	}
}
sub run_clean
{
	my $dir = shift;
	substr($dir, -1, 1, "") if( $dir =~ m/\/$/ );
	chdir( "$dir" ) or die "invalid directory $dir\n";
	&clean_dir( $dir );
	chdir( ".." );
}

my @dirs_def = qw/ hazy1 hazy2 hazy3 QuickStart /;

my @dirs;
push(@dirs, @ARGV)      if( @ARGV );
push(@dirs, @dirs_def)  if( not @dirs );

foreach my $dir ( @dirs ) { &run_clean( $dir ); }
