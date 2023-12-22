#!/usr/bin/perl
open IN, "refer1.tex";
open OUT, ">refer2.tex";
while ($f=<IN>) {
		if ($f =~ /^\\bibitem/) {
				$r = $c = "";
				while ($l=<IN>) {
						if ($l =~ /^%/) {
								$c = $l;
								next;
						}
						$r .= $l;
						last if $l =~ /^\s*$/;
				}
				if ($r =~ /^([\w ]*),[^,]*,\s*[\w ]*,[^,]*,\s*[\w ]*,[^,]*,?.*(\d{4,4}\w?),/ms) {
						$n = $1.' et al.';
						$y = $2;
				} elsif ($r =~ 
								 /^([\w ]*),[^,]*,\s*([\w ]*),.*\\&\s*([\w ]*),[^,]*,?.*(\d{4,4}\w?),/ms) 
				{
						$n = $1.', '.$2.' \& '.$3;
						$y = $4;
				} elsif ($r =~ 
								 /^([\w ]*),[^,]*,\s+\\&\s*([\w ]*),[^,]*,?.*(\d{4,4}\w?),/ms) 
				{
						$n = $1.' \& '.$2;
						$y = $3;
				} else {
						$r =~ /^([\w ]*).*(\d{4,4}\w?),/ms;
						$n = $1;
						$y = $2;
				}
				$f =~ s/bibitem/bibitem[$n($y)]/;
				print OUT $f,$c,$r;
		}
		else
		{
				print OUT $f;
		}
}
close OUT;
close IN;
