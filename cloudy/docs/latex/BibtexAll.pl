#!/usr/bin/perl -w
# compile all three parts of Hazy with pdflatex
#
chdir( "hazy1" ) or die "invalid directory hazy1\n";
system( "bibtex hazy1" );
#
chdir( "../hazy2" ) or die "invalid directory hazy2\n";
system( "bibtex hazy2" );
#
chdir( "../hazy3" ) or die "invalid directory hazy3\n";
system( "bibtex hazy3" );
#
chdir( "../QuickStart" ) or die "invalid directory hazy3\n";
system( "bibtex QuickStart" );
