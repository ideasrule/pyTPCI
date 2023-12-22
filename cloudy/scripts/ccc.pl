#!/usr/bin/perl

use warnings;
use strict;

use constant false => 0;
use constant true  => 1;

if( $#ARGV != 0 ) {
    die "usage: \"$0 <script-name>\"";
}

my $oldnam = $ARGV[0];
my $newnam = $oldnam . ".new";

open( my $ioIN, "<", $oldnam ) or die "could not open $oldnam";
open( my $ioOUT, ">", $newnam ) or die "could not open $newnam";

print "correcting $oldnam -> $newnam\n";

my $protect = false;
while( <$ioIN> ) {
    my $newline = $_;
    if( /^opti/i && /line/i ) {
	$protect = true;
    }
    if( /^save/i && /line/i && !/xspe/i && !/near/i && !/species/i ) {
	if( /emis/i ) {
	    $protect = true;
	}
	if( / rt /i ) {
	    $protect = true;
	}
	if( /zone/i && /cumu/i ) {
	    $protect = true;
	}
	if( /opti/i && /some/i ) {
	    $protect = true;
	}
    }
    if( /^prin/i && /line/i && / sum/i ) {
	$protect = true;
    }
    if( /^end/i ) {
	$protect = false;
    }
    my $test = $_;
    $test =~ s/##.*//;
    # c or C needs to be followed by whitespace
    # don't replace when inside a line list (could be a carbon species)
    if( /^c\s/i  && !$protect ) {
	$newline =~ s/^c/#/i;
    }
    # the rest below all create hidden comments, so replace with '##'
    # protect the special EOF marker '***' starting in column 1...
    elsif( /^\*/ && !/^\*\*\*/ ) {
	$newline =~ s/^\*/##/;
    }
    # the following comment characters may appear anywhere on the line
    # only replace the first instance of any of these
    # using $test avoids replacing '##' -> '###'
    elsif( $test =~ /[#%;]|\/\// ) {
	$newline =~ s/[#%;]|\/\//##/;
    }
    print $ioOUT $newline;
}

close($ioIN);
close($ioOUT);
