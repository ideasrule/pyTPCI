#!/usr/bin/perl
#
# Usage: perl addref.pl
#
$old = 'refer.bib';
$new = 'refernew.bib';
open IN, $old;
open OUT, ">$new";
%adstags = ( 
						 'A&A' => 'A&A..', 
						 'A&AS' => 'A&AS.',
						 'ADNDT' => 'ADNDT',
             'American Journal of Physics' => 'AmJPh',
						 'Ann. Phys.' => 'AnPhy',
						 'ApJ' => 'ApJ..',
						 'ApJS' => 'ApJS.',
						 'Ap&SS' => 'Ap&SS',
						 'ARAA' => 'ARA&A',
						 'ARA&A' => 'ARA&A',
						 'Astrophysical Letters' => 'ApL..',
						 'At Dat Nuc Dat Tab' => 'ADNDT',
						 'BAAS' => 'BAAS.',
						 'Bull. Astr. Inst. Netherlands' => 'BAN..',
						 'Geochim. Cosmochim. Acta' => 'GeCoA',
						 'J. Comp. Appl. Math' => 'JCoAM',
						 'J Phys B' => 'JPhB.',
						 'J. Phys. B.' => 'JPhB.',
						 'J.Phys. B' => 'JPhB.',
						 'J.Phys.B.' => 'JPhB.',
						 'MNRAS' => 'MNRAS',
						 'Nature' => 'Natur',
						 'Observatory' => 'Obs..',
						 'PASP' => 'PASP.',
						 'Phys Rev' => 'PhRv.',
						 'PhysRevA' => 'PhRvA',
						 'Phys Rev A' => 'PhRvA',
						 'Phys Scrip' => 'PhyS.',
						 'Rep. Prog. in Physics' => 'RPPh.',
						 'Revista Brasileira de Fisica' => 'RBrFi',
						 'Rev Mexicana' => 'RMxAA',
						 'Rev Mod Phys' => 'RvMP.',
						 'Rev. Mod. Phys' => 'RvMP.',
						 'RMxAC' => 'RMxAC',
						 'Science' => 'Sci..',
						 'Solar Physics' => 'SoPh.',
						 'Space Science Review' => 'SSRv.',
						 'Space Science Reviews' => 'SSRv.',
						 'Space Sci. Rev.' => 'SSRv.',
						 'Z. Physic' => 'ZPhy.'
						 );

while (<IN>)
{
		next if /^@/ or /^$/;
		$ADSCODE="ABSENT";
		if (/(\w).*(\d{4,4})\w?,\s*([^\d]+)\s*,\s*(\w?\d+)\s*,\s*(\w?\d+)/) { #,(\w+).*,\s*(\d+),\s*(\d+)/) {
#				print "$1 $2 $3 $4 $5\n";
				$vl = $4;
				$vl = ('.' x (4-length($vl))).$vl;
				$pg = $5;
				$pg = ('.' x (5-length($pg))).$pg;
				if (exists $adstags{$3}) {
						$ADSCODE = "$2$adstags{$3}$vl$pg$1";
				}
		} else {
#				print $_;
		}
		print OUT $_;
		print OUT "\@ARTICLE{$ADSCODE}\n\n";
}
close IN;
close OUT;
unlink $old;
rename $new, $old;
exit;
