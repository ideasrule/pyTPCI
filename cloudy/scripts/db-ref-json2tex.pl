#!/usr/bin/perl
#
# SUMMARY:
#
# Generate a TeX file suitable for publication holding the references for the
# data available in the data base of choice (Stout).
#
# DESCRIPTION:
#
# By default, the script produces a standard table (via the table* environment),
# but it can also format the table with the deluxetable environment, defined
# in deluxetable.sty.  The latter is included by default when using the aastex
# document class, but it can also be used standalone, obtained from
#       http://casa.colorado.edu/~danforth/comp/tex/deluxetable.sty
# An up-to-date link may be found here:
#       https://aas.org/faq/there-separate-deluxetable-style-file-i-can-use-outside-aastex
# To produce a deluxetable, the script must be called with the -d flag.
#
# The table holds four columns that identify the species, as well as the list of
# references for the energy levels, transition probabilities, collision
# strengths separately.  Species that are included in the default masterlist
# are distinguished with a double dagger on their right.
#
# The table is sorted in order primarily of increasing Z, and then of increasing
# ionization stage.  In regular table* mode, the number of rows may be
# controlled with the -nr option.  Species are not broken up across table breaks
# so this may be lead to some variation in table length from one table to the
# next.  This does not apply when deluxetable is used, because that environment
# figures out its own table breaks.  In that case, species may end up broken up
# across tables.
#
# By default, each species is followed by a horizontal rule for distinction.  In
# deluxetable a rule may appear as the top entry to the table data, just under
# the horizontal line that bounds the column headings.  Some fine-tuning of the
# table by hand may be needed to produce visually appealing results.
#
# The font size across all tables may be controlled with the -f option.
#
# Finally, the -e flag may be used to produce a minimal TeX document, and a
# script to compile it.  This is offered as an aid to fine-tune the table
# before using it in a larger TeX document.
#
# Usage examples:
#
# 1. Convert the references file (data/stout/refs.json) to a TeX table that
#    contains up to 35 rows of data per page, in small font, and generate a
#    TeX file and a script to test the output table:
#    
# 	./db-ref-json2tex.pl -nr=35 -f=s -e
#
# 2. As above, but format the output with deluxetable, and print it at the
#    default scriptsize font size:
#
# 	 ./db-ref-json2tex.pl -nr=35 -e -d
#
# Chatzikos, 2016-Jan-28
# Chatzikos, 2022-Apr-12
# 	Bugfix: Count TeX rows only if data ref exists
# Chatzikos, 2023-Aug-17
# 	Update to new structure of references, containing datasets instead of
# 	references within files.  Alternate datasets are shown in the Species
# 	column, prepended by 'or'.
# Chatzikos, 2023-Aug-18
# 	Fix bug with number of lines per table.
# Chatzikos, 2023-Oct-10
# 	Update to Perl 5.34.
#

use warnings;
use strict;

use Cwd;
use File::Basename;
use Getopt::Long;
use JSON;

use lib ".";
use BiblioToTeX;

# Immediately flush output
#
$| = 1;

my @save_ARGV = @ARGV;


#################################################################################
#			COMMAND-LINE INTERFACE					#
#################################################################################
my $numRowsDef = 40;
my $fontSizeDef = 'r';
sub printUsage
{
	my $program = &File::Basename::basename( $0 );
	print	"Usage:\n"
	 .	"   \t$program [-h] [-d] [-f=<f>] [-db=<db>] [-e]\n"
	 .	"where:\n"
	 .	"\t-h      :\t show this message and exit\n" 
	 .	"\t-d      :\t use the deluxetable format [OPTIONAL]\n";
	&BiblioToTeX::print_numRows_help( $numRowsDef );
	&BiblioToTeX::print_fontSize_help( $fontSizeDef );
	print	"\t-db=<db>:\t data base, one of: stout, chianti, lamda, all [OPTIONAL]\n";
	&BiblioToTeX::print_exampleTeX_help();
	die	"\n";
}

sub getInput
{
	my( $showUsage, $deluxeTable, $database, $numRows, $fontSize,
		$makeExample );

	&Getopt::Long::GetOptions
		(
			"h"	=>	\$showUsage,
			"d"	=>	\$deluxeTable,
			"nr=i"	=>	\$numRows,
			"f=s"	=>	\$fontSize,
			"db=s"	=>	\$database,
			"e"	=>	\$makeExample,
		)
	or die "Error in command-line arguments\n";

	&printUsage()	if( defined( $showUsage ) );

	my $db = &BiblioToTeX::check_database( $database );
	$numRows = &BiblioToTeX::check_numRows( $numRows, $numRowsDef );
	$fontSize = &BiblioToTeX::check_fontSize( $fontSize, $fontSizeDef );

	return	( $deluxeTable, $numRows, $fontSize, $db, $makeExample );
}


#################################################################################
#					TEX					#
#################################################################################
sub create_TeX_file
{
	my( $db ) = @_;

	my $TeX = "$db-refs.tex";
	my $FILE;
	open $FILE, "> $TeX"
	  or die "Error: Could not open file: $TeX\n";

	return	( $TeX, $FILE );
}

sub print_TeX_table_header
{
	my( $deluxeTable, $fontSize, $db, $FILE, $ntable_breaks ) = @_;

	my $font_cmd = &BiblioToTeX::getFontSize( $fontSize );

	if( defined( $deluxeTable ) )
	{
		print $FILE '\renewcommand{\baselinestretch}{1.0}' ."\n";
		print $FILE '\begin{deluxetable}{llll}' ."\n";
		print $FILE "\t". '\tablewidth{0pt}' ."\n";
		print $FILE "\t". "\\tabletypesize{$font_cmd}\n";
		print $FILE "\t". '\tablecolumns{4}' ."\n";
		print $FILE "\t". '\tablecaption' ."\n";
		print $FILE "\t". "{\n";
		print $FILE "\t\t". "References for the ".
				$BiblioToTeX::db_title{$db} ." database.\n";
		print $FILE "\t\t". 'Species marked with $\ddagger$ are listed '
		  .		"in Stout.ini, the default masterlist.\n";
		print $FILE "\t\t". "Alternate species are shown in the first"
		  .		" column, prepended by `or'.\n";
		print $FILE "\t\t". '\label{table:'. $db .'-refs}' ."\n";
		print $FILE "\t". "}\n";
		print $FILE "\t". '\tablehead' ."\n";
		print $FILE "\t". "{\n";
		print $FILE "\t\t". '\colhead{Species}' ."\t&\n";
		print $FILE "\t\t". '\colhead{Energy}' ."\t&\n";
		print $FILE "\t\t". '\colhead{Transition}' ."\t&\n";
		print $FILE "\t\t". '\colhead{Collision}' ."\n";
		print $FILE "\t". "}\n";
		print $FILE "\t". '\startdata' ."\n";
	}
	else
	{
		print $FILE '\begin{table*}' ."\n";
		print $FILE "\t". "{$font_cmd\n";
		print $FILE "\t". '\addtocounter{table}{-1}' ."\n"
			if( $ntable_breaks > 0 );
		print $FILE "\t". '\centering' ."\n";
		print $FILE "\t". '\caption'	."\n";
		print $FILE "\t". "{\n";
		if( $ntable_breaks == 0 )
		{
			print $FILE "\t". "\tReferences for ".
				$BiblioToTeX::db_title{$db} ." database.\n";
			print $FILE "\t\t". 'Species marked with $\ddagger$ are listed '
			  .		"in Stout.ini, the default masterlist.\n";
			print $FILE "\t\t". "Alternate species are shown in the first"
			  .		" column, prepended by `or'.\n";
		}
		else
		{
			print $FILE "\t\t". '{\it (Cont.{})}' ."\n"
		}
		print $FILE "\t". "}\n";
		print $FILE "\t". '\label{table:'. $db .'-refs}' ."\n"
			if( $ntable_breaks == 0 );

		print $FILE "\t". '\begin{tabular}{|l|p{5cm}|p{5cm}|p{5cm}|}' ."\n";
		print $FILE "\t\t". '\hline' ."\n";
		print $FILE "\t\t". '\hline' ."\n";
		print $FILE "\t\t\t". 'Species &' ."\n";
		print $FILE "\t\t\t". 'Energy &' ."\n";
		print $FILE "\t\t\t". 'Transition &' ."\n";
		print $FILE "\t\t\t". 'Collision' ."\n";
		print $FILE "\t\t". '\\\\' ."\n";
		print $FILE "\t\t". '\hline' ."\n";
	}
}

sub print_TeX_table_footer
{
	my( $deluxeTable, $FILE ) = @_;

	if( defined( $deluxeTable ) )
	{
		print $FILE "\t". '\enddata' ."\n";
		print $FILE '\end{deluxetable}' ."\n";
	}
	else
	{
		print $FILE "\t". '\end{tabular}' ."\n";
		print $FILE "\t". '}' ."\n";
		print $FILE '\end{table*}' ."\n";
	}
}

sub get_nrows_in_TeX_table
{
	my( $species ) = @_;

	my %nrows;
	for my $dataset ( keys %{ $$species{ref} } )
	{
		my $ref = $$species{ref}{$dataset};
		$nrows{$dataset} = {};

		$nrows{$dataset}{energy} = @{ $$ref{energy} }
			if( exists $$ref{energy} );
		$nrows{$dataset}{trans} = @{ $$ref{trans} }
			if( exists $$ref{trans} );
		$nrows{$dataset}{coll} = @{ $$ref{coll} }
			if( exists $$ref{coll} );
	}

	my $nrows = 0;
	for my $dataset ( keys %nrows )
	{
		next if( not %{ $nrows{$dataset} } );
		$nrows{$dataset}{max} =
			&List::Util::max( values %{ $nrows{$dataset} } );
		$nrows += $nrows{$dataset}{max};
	}
	print "$$species{element}\t$$species{ion}\t $nrows\n" if 0;

	return	( $nrows, { %nrows } );
}

sub fill_species_column
{
	my( $FILE, $species, $dataset, $deluxeTable ) = @_;

	print	$FILE	"\t\t\t";
	if( defined( $species ) )
	{
		printf $FILE "%s", ucfirst( $$species{element} );

		if( $$species{ion} == 2 )
		{
			print $FILE '$^{+}$';
		}
		elsif( $$species{ion} > 2 )
		{
			printf $FILE "\$^{%d+}\$", $$species{ion}-1;
		}
		print $FILE "\t". '$\ddagger$'
		  if( $$species{list} eq "default" );
	}
	elsif( defined( $dataset ) )
	{
		print $FILE "or $dataset:";
	}
	else
	{
		if( defined( $deluxeTable ) )
		{
			print $FILE '\nodata';
		}
	}
	print $FILE "\t&\n";
}

sub print_species_multiple_TeX_rows
{
	my( $deluxeTable, $FILE, $species, $nrows_per_dset_file ) = @_;

	#	print "$$species{element}\t$$species{ion}\n";

	my @file_order = qw/ energy trans coll /;
	my $dataset_order = &BiblioToTeX::set_dataset_order( $species );
#	print "$$species{element}\t$$species{ion}\t@$dataset_order\n";

	for my $dataset ( @$dataset_order )
	{
		if( $dataset eq $BiblioToTeX::default_dataset )
		{
			&fill_species_column( $FILE, $species );
		}
		else
		{
			&fill_species_column( $FILE, undef, $dataset, $deluxeTable );
		}

		my $ref = $$species{ref}{$dataset};

		for( my $i = 0; $i < $$nrows_per_dset_file{$dataset}{max}; $i++ )
		{
			if( $i > 0 )
			{
				&fill_species_column( $FILE, undef, undef, $deluxeTable );
			}

			foreach my $file ( @file_order )
			{
				my $citation;
				print	$FILE	"\t\t\t";
				if( exists $$nrows_per_dset_file{$dataset}{$file}
				    and defined( $$nrows_per_dset_file{$dataset}{$file} )
				    and $i < $$nrows_per_dset_file{$dataset}{$file} )
				{
					$citation = $$ref{$file}[ $i ]{bibcode};
					if( defined( $citation ) )
					{
						my $crossref =
							&BiblioToTeX::get_crossref_abbrv( $citation );
						if( not defined( $crossref ) )
						{
							$crossref = $citation;
						}
						$citation = "\\citet{". $crossref ."}";
					}
					else
					{
						$citation = $$ref{$file}[ $i ]{name};
						if( $citation =~ m/(\w+ NIST|NIST \w+)/ and
						    $citation !~ m/ \d\d\d\d-\d\d-\d/ )
						{
							#print "$citation\t =>\t";
							$citation = 'NIST';
							#print "$citation\n";
						}
					}
					print $FILE $citation
						if( defined( $citation ) );
				}
				elsif( $i == 0 and $file eq "coll" )
				{
					print $FILE 'baseline';
				}
				else
				{
					if( defined( $deluxeTable ) )
					{
						print $FILE '\nodata';
					}
				}
				print $FILE "\t&"
					if( $file ne $file_order[-1] );
				print $FILE "\n";
			}
			print $FILE "\t\t". '\\\\' ."\n";
		}
	}
	print	$FILE	"\t\t". '\hline' ."\n";
}

sub print_TeX_table_Stout
{
	my( $deluxeTable, $numRows_perTable, $fontSize, $db, $FILE, $species ) = @_;

	my $numRows_perTable_current = 0;
	my $ntable_breaks = 0;

	&print_TeX_table_header( $deluxeTable, $fontSize, $db, $FILE,
					$ntable_breaks );

	foreach my $element ( @BiblioToTeX::elements_ordered )
	{
		my @these_species;
		foreach my $sp ( sort keys %$species )
		{
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

		foreach my $ispecies ( @order )
		{
			my $species = $these_species[ $ispecies ];
			my( $nrows_for_species, $nrows_per_dset_file ) =
					&get_nrows_in_TeX_table( $species );

			next	if( $nrows_for_species == 0 );

			if( not defined( $deluxeTable )
			    and $numRows_perTable_current + $nrows_for_species
		    		> $numRows_perTable )
			{
				$ntable_breaks++;
				$numRows_perTable_current = $nrows_for_species;
				&print_TeX_table_footer( $deluxeTable, $FILE );
				&print_TeX_table_header( $deluxeTable, $fontSize,
						$db, $FILE, $ntable_breaks );
			}
			else
			{
				$numRows_perTable_current += $nrows_for_species;
			}

			&print_species_multiple_TeX_rows( $deluxeTable, $FILE,
				$species, $nrows_per_dset_file );
		}
	}

	&print_TeX_table_footer( $deluxeTable, $FILE );

	return;
}

sub print_TeX_table
{
	my( $deluxeTable, $numRows, $fontSize, $db, $species ) = @_;

	my( $TeX, $FILE ) = &create_TeX_file( $db );

	&BiblioToTeX::print_comment_command_line( $FILE, @save_ARGV );

	if( $db eq "stout" )
	{
		&print_TeX_table_Stout( $deluxeTable, $numRows, $fontSize,
					$db, $FILE, $species );
	}

	close $FILE
	   or warn "Warning: Could not close file: $TeX\n";

	print "Created:\t $TeX\n";

	return	$TeX;
}

sub main
{
	my( $deluxeTable, $numRows, $fontSize, $db, $exampleTeX ) = &getInput();

	&BiblioToTeX::set_globals();
	&BiblioToTeX::order_elements();

	&BiblioToTeX::load_cloudy_bibliography();

	foreach my $db ( @$db )
	{
		my $curdir = &Cwd::cwd();

		my $dbdir = $BiblioToTeX::data_dir ."/$db";
		chdir $dbdir
		   or die "Error: Could not chdir to db: $dbdir\n";

		my $species = &BiblioToTeX::load_json_or_die();

		chdir $curdir
		   or die "Error: Could not chdir back to: $curdir\n";

		my $tableTeX =
			&print_TeX_table( $deluxeTable, $numRows, $fontSize, $db,
						$species );
		&BiblioToTeX::create_exampleTeX_Script( $deluxeTable, $tableTeX, 1 )
			if( defined( $exampleTeX ) );
	}
}

&main();
