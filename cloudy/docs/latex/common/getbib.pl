#!/usr/bin/perl
#
# Usage: perl getbib.pl
#

# -- replaces lines in the form @ABSTRACT{ADSCODE.....} with the full BiBTeX
# reference, automatically downloaded from the NASA ADS

use LWP::Simple;

# http://manas.tungare.name/projects/isbn-to-bibtex/?isbn=0201835959

$old = 'bibliography.bib';
$new = 'bibliographynew.bib';
open IN, $old;
open OUT, ">$new";
while (<IN>)
{
    if (/^\@\w*\{([^\}]{19})\}/) {
	$ref = $1;
	$resp = get "http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=$ref&data_type=BIBTEX&db_key=AST&nocookieset=1";
	print ".";
	if ($resp eq "") {
	    print "Did not find $ref\n";
	    print OUT $_;
	} else {
	    foreach $l (split("\n",$resp)) {
		next if $l =~ /^\s*$/;
		next if $l =~ /^(Retrieved|Query)/;
		print OUT "$l\n";
	    }
	}
    } else {
	print OUT $_;
    }
}
print "\n";
close IN;
close OUT;
unlink $old;
rename $new, $old;
