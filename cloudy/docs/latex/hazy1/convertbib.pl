#!/usr/bin/perl

$refdone = 0;
$l="";
open IN, "refer.tex";
open OUT, ">refer1.tex";
foreach $n (<IN>) {
		if ($n =~ /\\chapter\{REFERENCES\}/) {
				print OUT "\\bibliographystyle{plainnat}\n\\begin{thebibliography}{}\n";
		} elsif ($n =~ /^\s*$/) {
				printmeta($l);
				print OUT $l;
				print OUT $n;
				$l = "";
		} else {
				$l .= $n;
		}
}
if ($l) {
		printmeta($l);
		print OUT $l;		
}
print OUT "\\end{thebibliography}\n";
close OUT;
close IN;

sub predot {
		my ($w,$l) = @_;
		return ('.'x($l-length($w))).$w;
}
sub posdot {
		my ($w,$l) = @_;
		return $w.('.'x($l-length($w)));
}
sub printmeta {
		my ($line) = @_;
		if ($line =~ /(\w*).*(\d{4,4}\w?),/ms) {
				$code = $1.$2;
				print "Duplicate $code\n" if exists $refs{$code};
				$refs{$code} = 1;
				print OUT "\\bibitem\{${1}$2\}\n";
		}
		if ($line =~ /(\w).*\D(\d+)\w?,\s*(.+),\s*(\d+),\s*(L)?(\d+)/ms) {
				my $i = $1;
				my $y = $2;
				my $j = $3;					 
				my $v = $4;
				my $p = $6;
				my $l = $5;
				$j =~ s/\\&/&/;
				$j =~ s/Astrophysical Letters/ApL/;
				$j =~ s/Phys\s*Rev\s*/PhRv/;
				$j =~ s/Geochim\.\s*Cosmochim\.\s*Acta/GeCoA/;
				$j =~ s/Space\s*Science\s*Reviews/SSRv/;
				$j =~ s/Space\s*Sci\.?\s*Rev\.?$/SSRv/;
				$j =~ s/At\s*Dat\s*Nuc\s*Dat\s*Tab/ADNDT/;
				$j =~ s/Science/Sci/;
				$j =~ s/Rev\.?\s*Mod\.?\s*Phys\.?/RvMP/;
				$j =~ s/Z\.?\s*Physik/ZPhy/;
				$j =~ s/Bull\.?\s*Astr\.?\s*Inst.\.?\s*Netherlands/BAN/;
				$j =~ s/J\.?\s*Comp\.?\s*Appl\.?\s*Math/JCoAM/;
				$j =~ s/American\s*Journal\s*of\s*Physics/AmJPh/;
				$j =~ s/Phys Scrip T/PhST/;
				$j =~ s/Physics Reports/PhR/;
				$j =~ s/Solar Physics/SoPh/;
				$j =~ s/Rev Mexicana/RMxAA/;
				$j =~ s/Nature/Natur/;
				print OUT '%ADS: ',predot($y,4).posdot($j,5).predot($v,4).
						predot($l,1).predot($p,4).$i."\n" if not $n =~ /\%ADS/;
		}
}
