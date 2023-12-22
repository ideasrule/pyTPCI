#!/usr/bin/perl
#
# Produce a file that contains a TeX table suitable for publication. The
# table lists all the atomic species used in Cloudy along with the database
# that holds the relevant data (Stout, Chianti, or none).  Because this
# depends on how Cloudy is called, the script operates on the output of the
# 'save species labels' command, i.e., it makes no effort to resolve the
# default species sources, but leaves that to Cloudy itself.
#
# An example is the file
#	tsuite/auto/func_lines.spclab
# which may be produced by running the sim
#	tsuite/auto/func_lines.in
#
# By default, the script produces a standard table (via the table environment),
# but it can also format the table with the deluxetable environment, defined
# in deluxetable.sty.  The latter is included by default when using the aastex
# document class, but it can also be used standalone, obtained from
# 	http://casa.colorado.edu/~danforth/comp/tex/deluxetable.sty
# An up-to-date link may be found here:
#	https://aas.org/faq/there-separate-deluxetable-style-file-i-can-use-outside-aastex
# To produce a deluxetable, the script must be called with the -d flag.
#
# In total, the table constists of 465 atomic species, so it cannot fit in a
# single page.  For each species, the species symbol and the employed database
# are reported.  The rows are set in order primarily of increasing Z, and then
# of increasing ionization stage.  A more compact output may be produced by
# folding several pairs of 'species-db' columns in the same table.  The number
# of column pairs per table is controlled by the -np option.  In that case, the
# species are arranged column-first, so the last row of one column is followed
# by the first row of the following column, and so on, even across table breaks.
#
# The number of rows per table may be set with the -nr option.  In regular table
# mode, this determines the number of rows printed per table.  However, because
# deluxetable does its own calculations for breaking the table across pages, it
# most likely will end up listing a different number of rows per table than the
# user might have guessed.  This may be lead to species discontinuities in each
# column, past the user-specified row number.  E.g., in a column, O^{6+} may be
# followed by Mg^{4+}, because O^{7+} has already been placed in the first row
# of the next column.  Some experimentation is advised when folding several
# column pairs with deluxetable.
#
# The table may appear 'busy' when several column-pairs are printed per table.
# In deluxetable, they are separated by vertical lines.  Although the default
# regular table environment permits vertical lines, these are not shown with
# the MNRAS style (mn2e).  For that reason, extra blank columns are inserted
# between successive 'species-db' column pairs.  Their width is controlled by
# the -bc option.  A zero-length choice, e.g., '0cm', eliminates the extra
# columns.
#
# The default of the regular table environment is to print a table in half the
# page width.  Because the table may overflow into the other half, or outside
# of the page, the -a flag is provided to use the table* environment, instead,
# which prints only one table per page width.  (It may still overflow for a high
# number of folded column pairs.)  On the other hand, the default deluxetable
# mode is to print the table at its natural width.  The -a flag may be used to
# force the table to take up the entire page width, but usually that does not
# produce appealing results.
#
# The appearance of the table also depends on the font-size used, controlled
# with the -f option.
#
# Finally, the -e flag may be used to produce a minimal TeX document, and a
# script to compile it.  This is offered as an aid to fine-tune the table
# before using it in a larger TeX document.
#
# Usage examples:
#
# 1. Create a table that takes up the entire page width (table*), holding 4
#    pairs of 'species-db' columns across the page, each table containing 35
#    rows of data per page, in small font, and generate a TeX file and a
#    script to test the output table:
#
#	 ./db-species-tex.pl ../tsuite/auto/func_lines.spclab -e -a \
#	 	-np=4 -f=s -nr=35
#
#
# 2. Create a deluxetable holding 4 pairs of 'species-db' columns across the
#    page, each table containing 38 rows of data per page, in small font, and
#    generate a TeX file and a script to test the output table:
# 
# 	./db-species-tex.pl ../tsuite/auto/func_lines.spclab -e -d \
# 		-np=4 -f=s -nr=38
#
#
# Chatzikos, 2016-Jan-28
#

use strict;
use warnings;

use File::Basename;
use Getopt::Long;

use lib ".";
use BiblioToTeX;

my @save_argv = @ARGV;


#################################################################################
#			COMMAND LINE INTERFACE					#
#################################################################################
my $numColPairsDef = 1;
my $numRowsDef = 40;
my $fontSizeDef = 'r';
my $blankColumnSizeDef = '0.2cm';
sub printUsage
{
	my $prog = &File::Basename::basename( $0 );
	print	"Usage:\n"
	 .	"\t$prog [-h] [-d] [-a] [-np=<np> -nr=<nr>] [-bc=<bc>] [-f=<f>] [-e] sp_file\n"
	 .	"where:\n"
	 .	"\t-h      :\t show this message and exit\n"
	 .	"\t-d      :\t use the deluxetable format [OPTIONAL]\n"
	 .	"\t-a      :\t print table across the entire page [OPTIONAL]\n"
	 .	"\t-np=<np>:\t number of 'species-db' column pairs per row"
	 .		" [OPTIONAL; DEFAULT: $numColPairsDef]\n";
	&BiblioToTeX::print_numRows_help( $numRowsDef );
	print	"\t-bc=<bc>:\t size of blank columns between sets of 'species-db' columns\n"
	 .	"\t         \t in regular (not deluxe) format [OPTIONAL; DEFAULT: '$blankColumnSizeDef']\n"
	 .	"\t         \t use '0cm' to eliminate the blank columns altogether\n";
	&BiblioToTeX::print_fontSize_help( $fontSizeDef );
	&BiblioToTeX::print_exampleTeX_help();
	print	"\tsp_file :\t file output of command 'save species labels'\n"
	 .	"\t         \t run the script ../tsuite/auto/func_lines.in to get the default state of Cloudy\n"
	 .	"\n\tNOTE: Some experimentation will be required to get the deluxe table\n"
	 .	  "\t      to like right on paper, if fiddling with -np, -nr, or -f.\n"
	 .	  "\t      See script comments.";
	die	"\n";
}

sub getInput
{
	my( $showUsage, $deluxeTable, $acrossPage, $numColPairs, $numRows, $fontSize,
		$blankColumnSize, $exampleTeX );

	&Getopt::Long::GetOptions
		(
			"h"	=>	\$showUsage,
			"d"	=>	\$deluxeTable,
			"a"	=>	\$acrossPage,
			"np=i"	=>	\$numColPairs,
			"nr=i"	=>	\$numRows,
			"f=s"	=>	\$fontSize,
			"bc=s"	=>	\$blankColumnSize,
			"e"	=>	\$exampleTeX,
		)
	or die "Error in command-line arguments\n";

	&printUsage()	if( defined( $showUsage ) );

	if( not defined( $numColPairs ) )
	{
		$numColPairs = $numColPairsDef;
	}
	elsif( $numColPairs <= 0 )
	{
		die "Error: Invalid number of 'species-db' column pairs per table:\t $numColPairs\n";
	}

	$numRows = &BiblioToTeX::check_numRows( $numRows, $numRowsDef );

	$fontSize = &BiblioToTeX::check_fontSize( $fontSize, $fontSizeDef );

	if( not defined( $deluxeTable ) )
	{
		if( not defined( $blankColumnSize ) )
		{
			$blankColumnSize = $blankColumnSizeDef;
		}
		else
		{
			die "Error: Blank column size expected in format "
			 ,	"'0.5cm', got:\t$blankColumnSize\n"
				if( $blankColumnSize !~
					m/^\d*(\.\d+|)(pt|cm|mm|in|ex|em)$/ );
		}
		$blankColumnSize = '0cm'
			if( $blankColumnSize =~ m/^0[a-z][a-z]$/i );
	}

	my $species_list = shift( @ARGV );
	&printUsage()
		if( not defined( $species_list ) );
	die "Error: File does not exist:\t $species_list\n"
		if( not -e $species_list );
	die "Error: Empty file:\t $species_list\n"
		if( not -s $species_list );

	return	( $deluxeTable, $acrossPage, $numColPairs, $numRows,
			$fontSize, $blankColumnSize, $exampleTeX,
			$species_list );
}


#################################################################################
#				INPUT / OUTPUT					#
#################################################################################
sub read_species_list
{
	my $species_list = shift;
	my $contents = &BiblioToTeX::read_contents( $species_list );
	my %species;
	foreach my $line ( @$contents )
	{
		next	if( $line =~ m/^#/ );
		chomp( $line );
		my( $species, $db ) = split( /\t+/, $line );
		if( $line =~ m/lamda/i )
		{
			$species{$species}{molecule} = $species;
			$species{$species}{db} = $db;
		}
		elsif( $species !~ m/(CRP|CRPHOT|PHOTON|e-|grn)/ )
		{
			if( $species =~ m/(-|\*)/ )
			{
				$species{$species}{molecule} = $species;
				next;
			}
			my $plus_sign;
			$plus_sign = 1	if( $species =~ m/\+/ );
			my ( $element, $ion ) = split( '\+', $species );
			if( not defined( &BiblioToTeX::isAtomicElement( $element ) ) )
			{
				$species{$species}{molecule} = $species;
				next;
			}
			if( not defined( $plus_sign ) )
			{
				if( not defined( $ion ) )
				{
					$ion = 0;
				}
				else
				{
					die	"Error: No plus sign, but ion defined for: $line\n";
				}
			}
			else
			{
				if( not defined( $ion ) or $ion eq "" )
				{
					$ion = 1;
				}
				#	print "$line\t". ( $ion ) ."\n";
			}
			#	print	"$line\t =>\t $species\t $ion\t $db\n";
			$species{$species}{element} = $element;
			$species{$species}{ion} = $ion;
			$species{$species}{db} = $db;
		}
	}

	return	{ %species };
}

# Put all elements and species per element in order
# of increasing Z, and increasing ionization stage
#
sub order_species
{
	my $species = shift;

	my @all_species;
	foreach my $element ( @BiblioToTeX::elements_ordered )
	{
		my @these_species;
		foreach my $sp ( sort keys %$species )
		{
			next	if( defined( $$species{$sp}{molecule} ) );
			if( $$species{$sp}{element} =~ m/^$element$/i )
			{
				push( @these_species, $$species{$sp} );
			}
		}

		next	if( not @these_species );

		my @order =
			sort
			{
				$these_species[$a]{ion} <=> $these_species[$b]{ion}
			} (0 .. $#these_species );

		my @sorted_species;
		foreach my $ispecies ( @order )
		{
			my $species = $these_species[ $ispecies ];
			last
				if( $$species{ion} ==
					$BiblioToTeX::elements{
						ucfirst(lc($$species{element}))
					} );
			#	print	ucfirst(lc($$species{element})) ."\t".
			#		$$species{ion} ."\t". 
			#		$BiblioToTeX::elements{
			#			ucfirst(lc($$species{element}))
			#		} ."\n";
			push( @sorted_species, $species );
		}

		push( @all_species, @sorted_species );
	}

	foreach my $sp ( sort keys %$species )
	{
		next	if( not exists( $$species{$sp}{molecule} ) );
		push( @all_species, $$species{$sp} );
	}

	return	\@all_species;
}

sub print_TeX_table_header
{
	my( $deluxeTable, $acrossPage, $numColPairs_perTable, $fontSize,
		$blankColumnSize, $FILE, $TeX, $ntable_breaks ) = @_;

	my $font_cmd = &BiblioToTeX::getFontSize( $fontSize );

	if( defined( $deluxeTable ) )
	{
		&BiblioToTeX::print_deluxe_instructions( $FILE, $TeX );
		print $FILE '\renewcommand{\baselinestretch}{1.0}' ."\n";
		print $FILE '\begin{deluxetable}{';
		for( my $i = 0; $i < $numColPairs_perTable; $i++ )
		{
			print $FILE 'll';
			print $FILE '|'
				if( $i < $numColPairs_perTable-1 );
		}
		print $FILE "}\n";
		print $FILE "\t". '\tablewidth{0pt}' ."\n"
			if( not defined( $acrossPage ) );
		print $FILE "\t". "\\tabletypesize{$font_cmd}\n";
		print $FILE "\t". '\tablecolumns{'. (2*$numColPairs_perTable) .'}' ."\n";
		print $FILE "\t". '\tablecaption' ."\n";
		print $FILE "\t". "{\n";
		print $FILE "\t\t". "Species data bases.\n";
		print $FILE "\t\t". '\label{table:species-db}' ."\n";
		print $FILE "\t". "}\n";
		print $FILE "\t". '\tablehead' ."\n";
		print $FILE "\t". "{\n";
		for( my $i = 0; $i < $numColPairs_perTable; $i++ )
		{
			print $FILE "\t\t". '\colhead{Species}' ."\t&\n";
			print $FILE "\t\t". '\colhead{Database}';
			print $FILE "\t&"
				if( $i < $numColPairs_perTable-1 );
			print $FILE "\n";
		}
		print $FILE "\t". "}\n";
		print $FILE "\t". '\startdata' ."\n";
	}
	else
	{
		print $FILE '\begin{table';
		print $FILE "*"	if( $acrossPage );
		print $FILE "}\n";
		print $FILE "\t". "{$font_cmd" ."\n";
		print $FILE "\t". '\addtocounter{table}{-1}' ."\n"
			if( $ntable_breaks > 0 );
		print $FILE "\t". '\centering' ."\n";
		print $FILE "\t". '\caption'	."\n";
		print $FILE "\t". "{\n";
		if( $ntable_breaks > 0 )
		{
			print $FILE "\t\t". '{\it (Cont.{})}' ."\n";
		}
		else
		{
			print $FILE "\t\t". '\label{table:species-db}' ."\n";
			print $FILE "\t\t". "Species data bases.\n";
		}
		print $FILE "\t". "}\n";

		print $FILE "\t". '\begin{tabular}{|';
		for( my $i = 0; $i < $numColPairs_perTable; $i++ )
		{
			print $FILE 'l|l|';
			if( $blankColumnSize ne '0cm' )
			{
				print $FILE "p{$blankColumnSize}|"
					if( $i < $numColPairs_perTable-1 );

			}
			else
			{
				print $FILE '|'
					if( $i < $numColPairs_perTable-1 );
			}
		}
		print $FILE "}\n";
		print $FILE "\t\t". '\hline' ."\n";
		print $FILE "\t\t". '\hline' ."\n";
		for( my $i = 0; $i < $numColPairs_perTable; $i++ )
		{
			print $FILE "\t\t". 'Species' ."\t&\n";
			print $FILE "\t\t". 'Database';
			print $FILE "\t&"
				if( $blankColumnSize ne '0cm' and
					$i < $numColPairs_perTable-1 );
			print $FILE "\n";
			print $FILE "\t\t\t&\n"
				if( $i < $numColPairs_perTable-1 );
		}
		print $FILE "\t\t". '\\\\' ."\n";
		print $FILE "\t\t". '\hline' ."\n";
	}
}

sub print_TeX_table_footer
{
	my( $deluxeTable, $acrossPage, $FILE ) = @_;

	if( defined( $deluxeTable ) )
	{
		print $FILE "\t". '\enddata' ."\n";
		print $FILE '\end{deluxetable}' ."\n";
	}
	else
	{
		print $FILE "\t". '\end{tabular}' ."\n";
		print $FILE "\t". '}' ."\n";
		print $FILE '\end{table';
		print $FILE "*"	if( $acrossPage );
		print $FILE "}\n";
	}
}

sub get_array_chunk
{
	my( $array, $nelements ) = @_;

	my @chunk;
	if( @$array > $nelements )
	{
		@chunk = splice( @$array, 0, $nelements );
	}
	else
	{
		#	my $nsp = @$array;
		#	print "$nsp\n";
		#	print $$array[$nsp-1]{ion} ."\n";
		@chunk = splice( @$array, 0, @$array );
	}
	return	\@chunk;
}

sub print_one_species_in_TeX
{
	my( $deluxeTable, $FILE, $species ) = @_;

	#	print "$$species{element}\t$$species{ion}\n";

	# Column 1:	Species
	#
	if( defined( $$species{molecule} ) )
	{
		my $label = "";
		my $superscript_started = 0;
		my $subscript_started = 0;
		for( my $i = 0; $i < length( $$species{molecule} ); $i++ )
		{
			my $this_character = substr( $$species{molecule}, $i, 1 );
			if( $this_character eq '^' ) 
			{
				$label .= "\$^\{";
				$superscript_started = 1;
			}
			elsif( $this_character =~ m/^\d$/ )
			{
				if( $superscript_started == 0 and
					$subscript_started == 0 )
				{
					$label .= "\$_\{";
					$subscript_started = 1;
				}
				$label .= $this_character;
			}
			elsif( $this_character =~ m/^[a-zA-Z]$/ )
			{
				if( $superscript_started or $subscript_started )
				{
					$label .= "\}\$";
					$superscript_started = 0;
					$subscript_started = 0;
				}
				$label .= $this_character;
			}
			elsif( $this_character =~ m/^(\+|-|\*)/ )
			{
				if( $superscript_started or $subscript_started )
				{
					$label .= "\}\$";
					$superscript_started = 0;
					$subscript_started = 0;
				}
				$this_character = '\star'	if( $this_character eq '*' );
				$label .= "\$^\{". $this_character . "\}\$";
			}
			else
			{
				warn "Warning: Unmapped character: '$this_character'\n";
			}
		}
		if( $subscript_started or $subscript_started )
		{
			$label .= "\}\$";
			$superscript_started = 0;
			$subscript_started = 0;
		}
		printf $FILE	"\t\t\t%s",
			$label;
	}
	else
	{
		printf	$FILE	"\t\t\t%s",
			ucfirst( lc( $$species{element} ) );

		if( $$species{ion} == 1 )
		{
			print	$FILE	'$^{+}$';
		}
		elsif( $$species{ion} > 1 )
		{
			printf	$FILE
				"\$^{%d+}\$",
				$$species{ion};
		}
	}
	print	$FILE	"\t&\n";

	print	$FILE	"\t\t\t";
	if( not defined( $$species{db} ) or $$species{db} eq "" )
	{
		if( defined( $deluxeTable ) )
		{
			print	$FILE	'\nodata';
		}
		else
		{
			print	$FILE	'\ldots';
		}
	}
	else
	{
		print	$FILE	$$species{db};
	}
}

sub print_table_cont
{
	my( $FILE, $deluxeTable, $numColPairs_perTable, $numRows_perPage,
			$blankColumnSize, $species_chunk ) = @_;

	my $numRows = ( @$species_chunk >=
				$numRows_perPage * $numColPairs_perTable )
			? $numRows_perPage
			: int( @$species_chunk / $numColPairs_perTable ) +
				(( @$species_chunk % $numColPairs_perTable > 0 )
					? 1 : 0 );

	for( my $irow = 0; $irow < $numRows; $irow++ )
	{
		#	print "irow = $irow\t numRows = $numRows\n";
		for( my $icol = 0; $icol < $numColPairs_perTable; $icol ++ )
		{
			my $index = $irow + $icol * $numRows;
			if( $index < @$species_chunk )
			{
				#	print	"\ticol = $icol"
				#	  .	"\t ncols = $numColPairs_perTable\n";
				&print_one_species_in_TeX( $deluxeTable, $FILE,
						$$species_chunk[ $index ] );
				print	$FILE	"\t&"
					if( not defined( $deluxeTable ) and
						$blankColumnSize ne '0cm' and
						$icol < $numColPairs_perTable-1 );
				print	$FILE	"\t&"
					if( defined( $deluxeTable ) and
						$icol < $numColPairs_perTable-1 );
				print	$FILE	"\n";
				print	$FILE	"\t\t\t\t&\n"
					if( not defined( $deluxeTable ) and
						$icol < $numColPairs_perTable-1 );
			}
			else
			{
				print	$FILE	"\t\t\t\t&\n"
					if( not defined( $deluxeTable ) and
						$icol == $numColPairs_perTable-1 );
			}
		}
		print	$FILE	"\t\t\\\\\n";
		print	$FILE	"\t\t\\hline\n";
	}
}

sub print_TeX_table
{
	my( $deluxeTable, $acrossPage, $numColPairs_perTable, $numRows_perPage,
		$fontSize, $blankColumnSize, $species ) = @_;

	my $ordered_species = &order_species( $species );

	my $TeX = "species-db.tex";
	my $FILE;
	open $FILE, "> $TeX"
	  or die "Error: Could not open file: $TeX\n";

	&BiblioToTeX::print_comment_command_line( $FILE, @save_argv );

	my $ntable_breaks = 0;

	&print_TeX_table_header( $deluxeTable, $acrossPage,
				$numColPairs_perTable, $fontSize,
				$blankColumnSize, $FILE, $TeX,
				$ntable_breaks );

	do
	{
		my $species_chunk =
			&get_array_chunk( $ordered_species,
				$numColPairs_perTable * $numRows_perPage );

		if( not defined( $deluxeTable ) )
		{
			if( $ntable_breaks > 0 )
			{
				&print_TeX_table_footer( undef, $acrossPage, $FILE );
				&print_TeX_table_header( undef, $acrossPage,
					$numColPairs_perTable, $fontSize,
					$blankColumnSize, $FILE, $TeX,
					$ntable_breaks );
			}
			$ntable_breaks++;
		}

		&print_table_cont( $FILE, $deluxeTable,
			$numColPairs_perTable, $numRows_perPage,
			$blankColumnSize, $species_chunk );
	}
	while( @$ordered_species );

	&print_TeX_table_footer( $deluxeTable, $acrossPage, $FILE );

	close $FILE
	   or warn "Warning: Could not close file: $TeX\n";

	print "Created:\t $TeX\n";

	return	$TeX;
}


#################################################################################
#				MAIN PROGRAM					#
#################################################################################
my( $deluxeTable, $acrossPage, $numColPairs, $numRows, $fontSize,
	$blankColumnSize, $exampleTeX, $species_list ) = &getInput();

&BiblioToTeX::order_elements();

my $species = &read_species_list( $species_list );
my $tableTeX = &print_TeX_table( $deluxeTable, $acrossPage, $numColPairs,
			$numRows, $fontSize, $blankColumnSize, $species );
&BiblioToTeX::create_exampleTeX_Script( $deluxeTable, $tableTeX )
	if( defined( $exampleTeX ) );
