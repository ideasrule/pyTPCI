#!/usr/bin/perl
open IN, 'refer.tex';
foreach (<IN>) {
		if (/\\bibitem\[([^\]]*)\]\{(.*)\}/) {
				if (exists($code{$1})) {
						$code{$1} = ""; # duplicate
						print "Duplicate $1\n";
				} else {
						$code{$1} = $2;
				}
		}
}
close IN;

sub findcode {
		my ($f) = @_;
		my ($l) = $f;
		$l =~ s/\s+/ /mg;
		$l =~ s/(\d{4,4}\w?)$/\($1\)/;
		$l =~ s/\s*\(/\(/;
		return $code{$l} if exists $code{$l};
		return "";
}
sub process_t {
		my ($f) = @_;
		my ($c) = findcode($f);
		return '\\citet{'.$c.'}' if ($c);
		return $f;
}
sub process_p {
		my ($f) = @_;
		my ($c) = findcode($f);
		return '\\citealp{'.$c.'}' if ($c);
		return $f;
}

foreach $f (<*.tex>) 
{
		next if $f eq 'refer.tex';
		$new = $f.'.new';
		$old = $f.'.old';
		open IN, $f;
		open OUT, '>',$f.'.new';
		$p = "";
		foreach $l (<IN>)
		{
				if ($l =~ /^\s*$/)
				{
						#$p =~ s/((\w+,\s+)+\d{4,4})/$1 !/mg;
						# Aaaa (1999)
						# Aaaa \& Aaaa (1999)
						$p =~ s/((((van\s+)?[A-Z]\w*,\s+)?(van\s+)?[A-Z]\w*,?\s+(\\&|and)\s+)?(van\s+)?[A-Z]\w*\s+\(\d{4,4}\w?\))/process_t($1)/mge;
						# Aaaa et al. 1999
						$p =~ s/((van\s+)?[A-Z]\w*,?\s+et\s+al\.\s*\(\d{4,4}\w?\))/process_t($1)/mge;
						# Aaaa \& Aaaa 1999
						$p =~ s/((((van\s+)?[A-Z]\w*,\s+)?(van\s+)?[A-Z]\w*,?\s+(\\&|and)\s+)?(van\s+)?[A-Z]\w*\s+\d{4,4}\w?)/process_p($1)/mge;
						# Aaaa et al. 1999
						$p =~ s/((van\s+)?[A-Z]\w*,?\s+et\s+al\.\s*\d{4,4}\w?)/process_p($1)/mge;
						print OUT $p,$l;
						$p="";
				}
				else
				{
						$p .= $l;
				}
		}
		print OUT $p,$l;
		close OUT;
		close IN;
		rename $f, $old;
		rename $new, $f;
}
