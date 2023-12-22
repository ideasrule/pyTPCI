#!/usr/bin/perl
$code = "";
open IN, "refer.tex";
open OUT, ">bibliography.bib";
while ($f=<IN>) {
    if ($f =~ /^\\bibitem\[[^\]]*\]\{(\w*)\}/) {
	$lab = $1;
    } elsif ($f =~ /%([^:]*):\s*(.*)/) {
	$src = $1;
	$code = $2;
	if ($src ne "ADS") {
	    $code = $src.':'.$code;
	}	
    } elsif ($f =~ /^\s*$/) {
	if ($code) {
	    print OUT "\@ARTICLE{$lab,\n crossref = \"$code\"\n}\n";
	    print OUT "\@ARTICLE{$code}\n";
	    $code = "";
	}
    }
    print OUT $f;
}
close OUT;
close IN;
