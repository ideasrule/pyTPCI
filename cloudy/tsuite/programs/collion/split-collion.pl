#!/usr/bin/perl
#######################################################################################
#
#  split-collion.pl::
#
#	This script may be used to transform the contents of 'collion.txt' to a
#  format suitable for plotting (eg, with Gnuplot).
#
#  	It splits 'tsuite/programs/collion.txt' into several files, one per element.
#  Each file contains the ionization fractions of all possible stages; for elements
#  heavier than Silicon, the fractions are brought together into a single table.
#
#	To use, create a link to 'collion.txt' in a directory, if not the present,
#  and run the script with no arguments.  This will generate files for all available
#  elements.  To produce output for a select number of elements, issue their symbols
#  (eg, 'h li fe') on the command line (case-insensitive).  The output files are
#  tab-delimited and suitable for use with Gnuplot or Veusz.
#
#	A Gnuplot script is also generated that produces one plot for each requested
#  element, which displays the temperature variation of the ionization fractions for
#  all ionization stages.  The output may be tweaked by modifying the User Options
#  below. Both PDF and PS can be generated, but that may depend on your system.
#
#	Note that the script has been tested only under Linux.  If problems occur
#  under MacOS X, or Windows, please let us know at cloudyastrophysics.groups.io
#
#
#  Chatzikos,	2013-Oct-14		First Version
#		2013-Oct-19		Added plots with Gnuplot.
#					Changed command line to accept chemical
#					symbols, instead of names.
#		2013-Oct-20		Output text files are now tab-delimited.
#					Introduced PS output, and parameterized
#					the script to allow easy choice of output.
#
#######################################################################################

use	strict;
use	warnings;


##############################
#	User Options
#
my $write_Gnuplot_file	= 1;
my $terminal_PDF	= "pdfcairo";	# Gnuplot PDF terminal
my $use_PDF		= 0;		# else, PS
my $use_PS_keys_right	= 1;		# else, keys go below plot
my $run_Gnuplot		= 1;
my $gnuplot		= "gnuplot";
my $convert_PS_to_PDF	= 0;		# in case of PS files
my $ps2pdf		= "ps2pdf";	# converter
my $keep_PS		= 1;

my $linewidth		= 3;
my $font		= "Helvetica, 15";	# type, size
#
# Some more parameters are defined below 
##############################


my %elements =
(
	H 	=>	{ Z =>  1, name => "Hydrogen"	},
	He	=>	{ Z =>  2, name => "Helium"	},
	Li	=>	{ Z =>  3, name => "Lithium"	},
	Be	=>	{ Z =>  4, name => "Beryllium"	},
	B 	=>	{ Z =>  5, name => "Boron"	},
	C 	=>	{ Z =>  6, name => "Carbon"	},
	N 	=>	{ Z =>  7, name => "Nitrogen"	},
	O 	=>	{ Z =>  8, name => "Oxygen"	},
	F 	=>	{ Z =>  9, name => "Fluorine"	},
	Ne	=>	{ Z => 10, name => "Neon"	},
	Na	=>	{ Z => 11, name => "Sodium"	},
	Mg	=>	{ Z => 12, name => "Magnesium"	},
	Al	=>	{ Z => 13, name => "Aluminium"	},
	Si	=>	{ Z => 14, name => "Silicon"	},
	P 	=>	{ Z => 15, name => "Phosphorus"	},
	S 	=>	{ Z => 16, name => "Sulphur"	},
	Cl	=>	{ Z => 17, name => "Chlorine"	},
	Ar	=>	{ Z => 18, name => "Argon"	},
	K 	=>	{ Z => 19, name => "Potassium"	},
	Ca	=>	{ Z => 20, name => "Calcium"	},
	Sc	=>	{ Z => 21, name => "Scandium"	},
	Ti	=>	{ Z => 22, name => "Titanium"	},
	V 	=>	{ Z => 23, name => "Vanadium"	},
	Cr	=>	{ Z => 24, name => "Chromium"	},
	Mn	=>	{ Z => 25, name => "Manganese"	},
	Fe	=>	{ Z => 26, name => "Iron"	},
	Co	=>	{ Z => 27, name => "Cobalt"	},
	Ni	=>	{ Z => 28, name => "Nickel"	},
	Cu	=>	{ Z => 29, name => "Copper"	},
	Zn	=>	{ Z => 30, name => "Zinc"	},
);



my @elements_req;
foreach my $elm ( @ARGV ) { push( @elements_req, $elements{ucfirst( lc( $elm ) )}{name} ); }
my $elements_req = join("\t", @elements_req);


my $inpfile = "collion.txt";
open FILE, "< $inpfile"		or die "Could not open:\t $inpfile\n";
my @contents = <FILE>;
close FILE			or warn "Could not close:\t $inpfile\n";



my @Elements;
foreach my $line ( @contents )
{
	next	if( $line !~ m/Element/ );
	my @words = split(/\s+/, $line);
	shift(@words)	while( $words[0] ne "Element" );
	shift(@words);
	shift(@words);
	next if( length($elements_req) and $elements_req !~ m/($words[0]\t|$words[0]$)/ );
	next if( join(' ', @Elements) =~ $words[0] );
	foreach my $elm ( keys %elements )
	{
		if( $elements{$elm}{name} eq $words[0] )
		{
			push( @Elements, $elm )
				if( join(' ', @Elements) !~ m/($elm\s|$elm$)/ );
			last;
		}
	}
}



foreach my $elm ( @Elements )
{
	my @output;
	my $i;
	for( $i = 0; $i < @contents; $i++ )
	{
		last	if( $contents[$i] =~ $elements{$elm}{name} );
	}

	my $j = 0;
	my $nelm = 0;
	for( ; $i < @contents; $i++ )
	{
		my $line = $contents[$i];
		chomp( $line );
		next	if( $line =~ m/^\s*$/ ); 
		last	if( $line =~ m/Element/ and
			    $line !~ $elements{$elm}{name} );
		if( $line =~ m/Element/ )
		{
			$nelm++;
			$j = 0;
			next;
		}

		my @fields = split(/\s+/, $line);
		while( $fields[0] =~ m/^\s*$/ )
		{
			shift( @fields );
		}

		if( $fields[0] eq "Te" )
		{
			for( my $i = 1; $i < @fields; $i++ )
			{
				if( $fields[$i] == 1 )
				{
					$fields[$i] = "$elm";
				}
				else
				{
					$fields[$i] = "$elm+". ($fields[$i]-1);
				}
			}
		}

		if( $nelm == 1 )
		{
			$output[$j] = join("\t", @fields);
		}
		else
		{
			shift( @fields );
			$output[$j] .= "\t". join("\t", @fields);
		}
		$j++;
	}

	for( my $k = 0; $k < @output; $k++ )
	{
		$output[$k] = "#". $output[$k]
			if( $output[$k] =~ m/Te/ );
		$output[$k] .= "\n";
	}

	my $file = "collion-". $elements{$elm}{name} .".txt";
	open FILE, "> $file"	or die "Could not open:\t $file\n";
	print FILE @output;
	close FILE		or die "Could not close:\t $file\n";

	$elements{$elm}{data} = $file;
}

exit
	unless( $write_Gnuplot_file );



my $plot = "MakeIonFrac.plot";
open GPLOT, "> $plot"	or die "Could not open:\t $plot\n";

# Additional parameters
my( $height, $width, $width_key_col, $height_key_col, $key_tmargin, $key_lmargin, $plot_bmargin );

if( $use_PDF or $use_PS_keys_right )
{
	$height = 15;
	$width  = 15;
	$width_key_col = 4;
}
else
{
	$height = 6;
	$width  = 10;
	$height_key_col = 2.0;
	$plot_bmargin   = $height_key_col / ($height+$height_key_col);
	$key_lmargin	= 0.07;
	$key_tmargin	= $plot_bmargin * 0.6;
}
print GPLOT "set xlabel \"log T [K]\"\n";
print GPLOT "set ylabel \"log f_{ion}\"\n";
print GPLOT "set xtics scale 1.5, 1\n";
print GPLOT "set mxtics 5\n";
print GPLOT "set ytics scale 1.5, 1\n";
print GPLOT "set mytics 5\n";
if( $use_PDF or $use_PS_keys_right )
{
	print GPLOT "set key outside right top vertical enhanced\n";
}
else
{
	print GPLOT "set key left top horizontal enhanced at screen $key_lmargin,$key_tmargin maxcols 7\n";
}
print GPLOT "\n";



my @psfiles;
foreach my $elm ( @Elements )
{
	my $ncol = 1;
	$ncol++
		if( $elements{$elm}{Z} >= 24 );
	if( $use_PDF )
	{
		print GPLOT "set terminal $terminal_PDF size ". int($width+$ncol*$width_key_col) ."cm,". int($height) ." color enhanced font \"$font\" lw $linewidth dashed\n";
	}
	elsif( $use_PS_keys_right )
	{
		print GPLOT "set terminal postscript size ". int($width+$ncol*$width_key_col) ."cm,". int($height) ."cm color enhanced font \"$font\" lw $linewidth dashed\n";
		print GPLOT "set rmargin at screen ". ($width / ($width+$ncol*$width_key_col)) ."\n";
	}
	else
	{
		print GPLOT "set terminal postscript size ". int($width) ."in,". int($height+$height_key_col) ."in color enhanced font \"$font\" lw $linewidth dashed\n";
		print GPLOT "set bmargin at screen $plot_bmargin\n";
		print GPLOT "set rmargin at screen 0.98\n";
	}

	my $psfile = $elm . &file_suffix();
	push( @psfiles, $psfile );

	print GPLOT "set output \"$psfile\"\n";
	print GPLOT "set title \"". $elements{$elm}{name} ."\"\n";
	print GPLOT "plot \t\"". $elements{$elm}{data} ."\" u 1:2 w l ". &get_lt( 0 ) ." t \"". $elm ."^{+0}\",\t\\\n";
	for( my $i = 1; $i <= $elements{$elm}{Z}; $i++ )
	{
		print GPLOT "\t\"\" u 1:". ($i+2) ." w l ". &get_lt( $i ) ." t \"". $elm ."^{+$i}\"";
		print GPLOT ",\t \\"
			if( $i != $elements{$elm}{Z} );
		print GPLOT "\n";
	}
	print GPLOT "\n";
}

close GPLOT	or warn "Could not close:\t $plot\n";

exit
	if( not $run_Gnuplot );


system( "$gnuplot $plot" );

exit
	if( $use_PDF or not $convert_PS_to_PDF );


foreach my $ps ( @psfiles )
{
	my $pdf = $ps;
	$pdf =~ s/\.ps$/.pdf/;
	system( "$ps2pdf $ps $pdf" );
	unlink $ps
		unless( $keep_PS );
}



sub file_suffix
{
	return	".pdf"	if( $use_PDF );
	return	".ps";
}



sub get_lt
{
	my ($ion) = @_;
	my $lt = int( $ion / 7 ) + 1;
	return	"lc $ion lt $lt";
}
