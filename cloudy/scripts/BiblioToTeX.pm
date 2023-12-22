# This package provides some functionality common to the scripts that treat the
# gathering of references from the atomic data base files, and the generation of
# TeX tables relevant to atomic data.
#
# Chatzikos, 2016-Jan-28
# Chatzikos, 2022-Apr-12
# 	Bugfix: In custom_to_json, protect against missing 'ref' key.
# Chatzikos, 2022-Apr-13
# 	Rename prep_bibtex to define_month_macros, and disable its calling
# 	-- current installation issues warnings like:
# 	   warning: overriding existing definition of macro "sep"
# Chatzikos, 2023-Aug-17
#	Change structure of JSON file: make 'ref' contain datasets, and move
#	its previous contents into the 'default' dataset.
#	Bugfix: Species masterlist was overwritten with value in JSON file.
# Chatzikos, 2023-Oct-10
#	Update to Perl 5.34.
#
use warnings;
use strict;

package BiblioToTeX;

require Exporter;
our @ISA = qw/ Exporter /;

use Cwd;
use File::Basename;
use Text::BibTeX;
use URI::Escape;
use JSON;
use Const::Fast;

# Immediately flush output
#
$| = 1;

#
# Global data
#
our $root_dir = "../";
our $data_dir = "data/";
our $bibliography;

const our @all_db => qw/ stout chianti lamda /;
const our %db_title =>
(
	stout	=>	"Stout",
	chianti	=>	"Chianti",
	lamda	=>	"LAMDA",
);


#################################################################################
#					CHEMISTRY				#
#################################################################################
our @elements_ordered;

our %elements =
(
	"H"	=>	 1,
	"He"	=>	 2,
	"Li"	=>	 3,
	"Be"	=>	 4,
	"B"	=>	 5,
	"C"	=>	 6,
	"N"	=>	 7,
	"O"	=>	 8,
	"F"	=>	 9,
	"Ne"	=>	10,
	"Na"	=>	11,
	"Mg"	=>	12,
	"Al"	=>	13,
	"Si"	=>	14,
	"P"	=>	15,
	"S"	=>	16,
	"Cl"	=>	17,
	"Ar"	=>	18,
	"K"	=>	19,
	"Ca"	=>	20,
	"Sc"	=>	21,
	"Ti"	=>	22,
	"V"	=>	23,
	"Cr"	=>	24,
	"Mn"	=>	25,
	"Fe"	=>	26,
	"Co"	=>	27,
	"Ni"	=>	28,
	"Cu"	=>	29,
	"Zn"	=>	30,
);

sub order_elements
{
	@elements_ordered = sort { $elements{$a} <=> $elements{$b} } keys( %elements );
}

sub isAtomicElement
{
	my $species = shift;
	my $isElement;
	foreach my $element ( @elements_ordered )
	{
		$isElement = 1
			if( $element =~ m/^$species$/i );
	}
	return	$isElement;
}


#################################################################################
#				   COMMAND-LINE					#
#################################################################################
sub check_database
{
	my $database = shift;

	my @db;
	if( not defined( $database ) )
	{
		push( @db, 'stout' );
	}
	else
	{
		if( $database eq 'stout' )
		{
			push( @db, 'stout' );
		}
		elsif( $database eq 'chianti' )
		{
			push( @db, 'chianti' );
		}
		elsif( $database eq 'lamda' )
		{
			push( @db, 'lamda' );
		}
		elsif( $database eq 'all' )
		{
			push( @db, @all_db );
		}
		else
		{
			die "Error: Unrecognized database option:\t $database\n";
		}
	}

	return	\@db;
}

sub get_string_command_line
{
	my @argv = @_;
	return	$0 ." ". join( " ", @argv );
}

#################################################################################
#					UTILITY					#
#################################################################################
sub set_globals
{
	my $curdir = &Cwd::cwd();
	$root_dir = &File::Basename::dirname( $curdir );
	$data_dir = $root_dir ."/". $data_dir;
	$bibliography = $root_dir . "/docs/latex/common/bibliography2.bib";
}

sub read_contents
{
	my $file = shift;
	open FILE, "< $file " or die "Error: Could not open file: $file\n";
	my @contents = <FILE>;
	close FILE or warn "Warning: Could not close file: $file\n";
	return	\@contents;
}

sub is_in_array
{
	my( $val, $array ) = @_;

	return 0
		if( not defined( $array ) );

	for ( @$array )
	{
		return 1 if $val eq $_;
	}

	return 0;
}


#################################################################################
#					STORAGE					#
#################################################################################
const my $JSON_filename => "refs.json";

#
# Default dataset identifier, used in JSON file.
# Other datasets are identified by the dataset string in their filenames;
# e.g., 'Tayal18' is the identifier for the fe_2_Tayal18.* set of files
#
const our $default_dataset => "default";  # used in JSON file

sub load_json
{
	my( $all_species ) = @_;

	if( -s $JSON_filename )
	{
		my $contents = join( '',
			@{ &BiblioToTeX::read_contents( $JSON_filename ) } );
		my $species_stored = &JSON::from_json( $contents );
		if( defined( $all_species ) )
		{
			foreach my $this_species ( sort keys %$species_stored )
			{
				#
				# Preserve the masterlist, just read in!
				#
				my $list = $$all_species{$this_species}{list};
				$$all_species{$this_species} =
					$$species_stored{$this_species};
				$$all_species{$this_species}{list} = $list;
			}
		}
		else
		{
			$all_species = $species_stored;
		}
	}
}

sub load_json_or_die
{
	die "Error: JSON file does not exist:\t$JSON_filename\n"
		if( not -e $JSON_filename );
	die "Error: JSON file is empty:\t$JSON_filename\n"
		if( not -s $JSON_filename );

	return	&load_json();
}

sub set_dataset_order
{
	my ( $species_data ) = @_;

	my @dataset_order = ( $default_dataset );

	for my $dataset ( sort keys %{ $$species_data{ref} } )
	{
		push( @dataset_order, $dataset )
			if( not &is_in_array( $dataset, \@dataset_order ) );
	}

	return \@dataset_order;
}

sub custom_to_json
{
	my( $hash ) = @_;

	my $tab = " " x 1;

	my @order_subhash = qw/element ion list/;

	&order_elements()
		if( not @elements_ordered );

	my @species;
	foreach my $element ( @elements_ordered )
	{
		for( my $ion = 1; $ion <= $elements{$element}; $ion++ )
		{
			my $species = lc( $element ) ."_". $ion;
			next
				if( not exists $$hash{$species} );
			push( @species, $species );
		}
	}

	my $string = "{\n";
	for( my $ispecies = 0; $ispecies < @species; $ispecies++ )
	{
		my $ntabs = 1;

		my $species = $species[ $ispecies ];
		$string .= ($tab x $ntabs) .'"'. $species .'" : {'."\n";

		$ntabs++;
		foreach my $subhash_elm ( @order_subhash )
		{
			$string .= ($tab x $ntabs) . '"'. $subhash_elm .'" : "'.
					$$hash{$species}{$subhash_elm} ."\",\n";
		}

		if( not exists $$hash{$species}{ref} )
		{
			# drop last comma
			$string = substr( $string, 0, -2 ) ."\n";
		}
		else
		{
			$string .= ($tab x $ntabs) .'"ref" : {'."\n";

			my $dataset_order = &set_dataset_order( $$hash{$species} );

			$ntabs++;
			for my $dataset ( @$dataset_order )
			{
				$ntabs++;
				$string .= ($tab x $ntabs) ."\"$dataset\" : {"."\n";

				my @datatypes = sort keys %{ $$hash{$species}{ref}{$dataset} };
				for my $datatype ( @datatypes )
				{
					$ntabs++;
					$string .= ($tab x $ntabs) ."\"$datatype\" : [";
					my $nrefs = @{ $$hash{$species}{ref}{$dataset}{$datatype} };
					if( $nrefs > 0 )
					{
						$string .= "\n";
						for( my $iref = 0; $iref < $nrefs; $iref++ )
						{
							my $ref = $$hash{$species}{ref}{$dataset}{$datatype}[ $iref ];
							$ntabs++;
							$string .= ($tab x $ntabs) ."{\n";
							$ntabs++;
							my @keys = sort keys %$ref;
							for( my $ikey = 0; $ikey < @keys; $ikey++ )
							{
								my $key = $keys[ $ikey ];
								$string .= ($tab x $ntabs) ."\"$key\" : \"$$ref{$key}\"";
								$string .= ","
									if( $ikey != @keys-1 );
								$string .= "\n";
							}
							$ntabs--;
							$string .= ($tab x $ntabs) ."}";
							$string .= ","
								if( $iref != $nrefs-1 );
							$string .= "\n";
							$ntabs--;
						}
						$string .= ($tab x $ntabs);
					}
					$string .=  "]";
					$string .= ","
						if( $datatype ne $datatypes[-1] );
					$string .= "\n";
					$ntabs--;
				}

				$string .= ($tab x $ntabs) ."}";
				$string .= ","
					if( $dataset ne $$dataset_order[-1] );
				$string .= "\n";
				$ntabs--;
			}

			$string .= ($tab x $ntabs) ."}\n";
			$ntabs--;
		}

		$string .= ($tab x $ntabs) ."}";
		$string .= ","
			if( $ispecies != @species-1 );
		$string .= "\n";
	}

	$string .= "}\n";
	return	$string;
}

sub store_json
{
	my( $species_data ) = @_;

	my $pp_species = &custom_to_json( $species_data );
	#	die $pp_species;
	#	my $pp_species = &JSON::to_json( $species_data, {pretty => 1} );

	open FILE, "> $JSON_filename"
	  or die "Error: Could not open:\t $JSON_filename\n";

	print FILE $pp_species;

	close FILE
	   or warn "Warning: Could not close:\t $JSON_filename\n";

	return;
}


#################################################################################
#					BIBTEX					#
#################################################################################
our @bibtex_basic_fields = qw/author title journal year volume pages /;

sub print_bibentry
{
	my( $entry, $keys ) = @_;

	my @keys;
	if( not defined( $keys ) )
	{
		@keys = sort keys %$entry;
	}
	else
	{
		@keys = @$keys;
	}

	printf "%10s\t=>\t%s\n", "type", $$entry{type}	if( exists( $$entry{type} ) );
	printf "%10s\t=>\t%s\n", "key", $$entry{key}	if( exists( $$entry{key} ) );
	foreach my $key ( sort keys %$entry )
	{
		next	if( $key eq 'type' or $key eq 'key' );
		next	if( not defined( $$entry{$key} ) );
		printf "%10s\t=>\t", $key;
		if( ref( $$entry{$key} ) eq "ARRAY" )
		{
			print  join( "; ", @{$$entry{$key}} ) ."\n";
		}
		else
		{
			print $$entry{$key} ."\n";
		}
	}
	print "\n";
}

sub parse_bibtex_entry
{
	my $entry = shift;

	my %this_entry;
	$this_entry{type} = $entry->type;
	$this_entry{key} = $entry->key;
	foreach my $field ( @bibtex_basic_fields )
	{
		$this_entry{$field} = $entry->get( $field );
	}
	$this_entry{crossref} = $entry->get( 'crossref' );
	if( not defined( $this_entry{author} ) )
	{
		$this_entry{author} = "";
		#	print "UNDEF:\n";
		#	&print_bibentry( { %this_entry } );
	}
	else
	{
		my @authors = $entry->split( 'author' );
		$this_entry{authors} = \@authors;
	}
	#	&print_bibentry( { %this_entry } );

	return	{ %this_entry };
}

sub convert_bibtex_to_hash
{
	my $bibdef = shift;

	my $entry = Text::BibTeX::Entry->new();
	$entry->parse_s( join( "\n", @$bibdef ) );
	#	&BiblioToTeX::print_bibentry( $entry );
	$$entry{year} = $entry->get( 'year' );
	$$entry{journal} = $entry->get( 'journal' );
	$$entry{volume} = $entry->get( 'volume' );
	$$entry{pages} = $entry->get( 'pages' );
	$$entry{title} = $entry->get( 'title' );

	my @authors = $entry->split( 'author' );
	$$entry{authors} = \@authors;

	return	$entry;
}

sub insert_bibtex_field
{
	my( $bibdef, $field, $value ) = @_;

	my $closing_brace = pop( @$bibdef );
	my $last_line = pop( @$bibdef );
	chomp( $last_line );
	$last_line .= ",\n"
		if( $last_line !~ m/, *$/ );
	push( @$bibdef, $last_line );

	my $line = sprintf( "%9s = {%s}\n", $field, $value );
	push( @$bibdef, $line );
	push( @$bibdef, $closing_brace );
}


#################################################################################
#				BIBLIOGRAPHY					#
#################################################################################
our %crossref;
our %biblio_order;

sub get_crossref_abbrv
{
	my $adscode = shift;

	# Convert HTML code to ASCII
	$adscode =~ s/%26/&/;

	my $humancode;
	foreach my $key ( %crossref )
	{
		$humancode = $crossref{$key}
			if( $key eq $adscode and
				defined( $crossref{$key} ) );
	}
	return	$humancode;
}

sub add_crossref_abbrv
{
	my( $humancode, $adscode ) = @_;
	$crossref{$adscode} = $humancode;
}

sub sort_biblio
{
	my( $biblio ) = @_;
	#	for( my $i=0; $i<3; $i++)
	#	{
	#		&print_bibentry( $$biblio[ $i ] );
	#	}
	#	die;

	my @order =
	sort
	{
		$$biblio[$a]{year} <=> $$biblio[$b]{year}
				  ||
		$$biblio[$a]{author} cmp $$biblio[$b]{author}
	}
	( 0 .. (@$biblio-1) );

	#	for( my $i=0; $i< @$biblio; $i++)
	#	{
	#		&print_bibentry( $$biblio[ $order[ $i ] ] );
	#	}
	#	die;

	foreach my $i ( @order )
	{
		my $year = $$biblio[ $i ]{year};
		#	print "i = $i\t year = $year\n";
		push( @{ $biblio_order{$year} }, $$biblio[ $i ] );
	}

	#	foreach my $year ( sort { $a <=> $b } keys( %biblio_order ) )
	#	{
	#		print "year = $year\n";
	#		foreach my $href ( @{ $biblio_order{$year} } )
	#		{
	#			&print_bibentry( $href );
	#		}
	#		print "\n";
	#	}
	#	die;
}

sub define_month_macros
{
	foreach my $month ( qw/jan feb mar apr may jun jul aug sep oct nov dec/ )
	{
		&Text::BibTeX::add_macro_text( $month, $month );
	}
}

sub load_cloudy_bibliography
{
	&define_month_macros()
		if 0;	# Enable if you get warnings about
			# unknown macros for months

	my $bibfile = Text::BibTeX::File->new( $bibliography ) ;

	my @biblio;
	while( my $entry = Text::BibTeX::Entry->new( $bibfile ) )
	{
		next	unless $entry->parse_ok;
		if( defined( $entry->get( 'crossref' ) ) )
		{
			&add_crossref_abbrv( $entry->key, $entry->get( 'crossref' ) );
		}
		else
		{
			push( @biblio, &parse_bibtex_entry( $entry ) );
		}
	}
	$bibfile->close;

	print "Acquired ". @biblio ." unique references and ".
		( keys %crossref ) ." cross references\n\n";

	&sort_biblio( \@biblio );
}

sub is_bibcode_in_Cloudy_bib
{
	my( $bibcode, $year ) = @_;

	die "Error: Undefined bibcode\n"
		if( not defined( $bibcode ) or $bibcode eq "" );
	die "Error: Undefined year\n"
		if( not defined( $year ) or $year eq "" );

	# Sanitize bibcode in case it comes from a URL.
	$bibcode = &URI::Escape::uri_unescape( $bibcode );

	my $found;
	foreach my $href ( @{ $biblio_order{$year} } )
	{
		#	print "\t bibcode: $bibcode\t $$href{key}\n";
		if( $$href{key} eq $bibcode )
		{
			$found = 1;
			last;
		}
	}
	return	$found;
}

sub drop_brackets
{
	my( $string, $ending ) = @_;
	return	""
		if( not defined( $string ) );
	$string =~ s/^\{//;
	if( defined( $ending ) )
	{
		$string =~ s/\}$//;
	}
	else
	{
		$string =~ s/\}//;
	}
	return	$string;
}

sub resolve_author_name
{
	my $author = shift;
	$author = &drop_brackets( $author );
	my( $lastname, $givennames ) = split( /,/, $author );
	$lastname =~ s/^ *//g;
	$lastname =~ s/ *$//g;
	if( defined( $givennames ) )
	{
		$givennames =~ s/^ *//g;
		$givennames =~ s/ *$//g;
	}
	return	( $lastname, $givennames );
}

sub get_bibcode_from_Cloudy_bib
{
	my $bt = shift;

	my @authors;
	foreach my $auth ( @{ $$bt{authors} } )
	{
		my( $lastname, $givennames ) = &resolve_author_name( $auth );
		last
			if( not defined( $lastname ) );
		$lastname =~ s/\W//g;
		#	print "$lastname\n"; 
		push( @authors, $lastname );
	}

	my( $bibcode, $title );
	foreach my $href ( @{ $biblio_order{$$bt{year}} } )
	{
		#	print "\t". $$href{author} ."\n";
		#	foreach my $auth ( @{ $$href{authors} } )
		#	{
		#		print "\t\t$auth\n";
		#	}
		#	print "\n";

		#	&print_bibentry( $href );

		next	if( not keys %{ $href} );

		next if( not defined $$href{authors} );
		next if( not exists( $$bt{authors_etal} ) and
				@{ $$href{authors} } < @authors );

		$bibcode = $$href{key};
		$title = $$href{title};
		for( my $i = 0; $i < @authors; $i++ )
		{
			last if( $authors[ $i ] =~ m/et al/ );
			#	print "'$$href{authors}[$i]'\t '$authors[$i]'\n";
			my $bib_auth = $$href{authors}[ $i ];
			$bib_auth =~ s/\W//g;
			if( $bib_auth !~ $authors[ $i ] )
			{
				#	print "$bib_auth\t $authors[$i]\n";
				$bibcode = $title = undef;
				last;
			}
		}

		if( defined( $bibcode ) and exists( $$bt{volume} ) )
		{
			#	print "Volume: $$bt{volume}\t ref Volume: $$href{volume}\n";
			$bibcode = $title = undef
				if( not exists( $$href{volume} ) or
					not defined( $$href{volume} ) or
					$$href{volume} != $$bt{volume} );
		}

		if( defined( $bibcode ) and exists( $$bt{pages} ) )
		{
			#	print "Pages : $$bt{pages}\t ref Pages: $$href{pages}\n";
			$bibcode = $title = undef
				if( not exists( $$href{pages} ) or
					not defined( $$href{pages} ) or
					$$href{pages} !~ $$bt{pages} );
		}

		last
			if( defined( $bibcode ) );
	}

	if( defined( $bibcode ) )
	{
		$title = &drop_brackets( $title, 1 );
		print	"Reference found in Cloudy biblio:\n"
		  .	"\tbibcode :\t $bibcode\n"
		  .	"\ttitle   :\t \"$title\"\n";
	}
	else
	{
		print	"Reference not found in Cloudy biblio\n";
	}

	return	$bibcode;
}

sub update_biblio_hash
{
	my $bibdef = shift;

	my $bt = &convert_bibtex_to_hash( $bibdef );
	#	print $$bt{year} ."\n";
	push( @{ $biblio_order{$$bt{year}} }, $bt );
}


#################################################################################
#					TEX					#
#################################################################################
sub print_comment_command_line
{
	my( $FILE, @argv ) = @_;

	print $FILE "%\n";
	print $FILE "% This file was created with the command:\n";
	print $FILE "%\t". &get_string_command_line( @argv ) . "\n";
	print $FILE "%\n";
}

sub print_deluxe_instructions
{
	my( $FILE, $TeXfile ) = @_;
	print $FILE "%\n"
		.   "% To use you must include one of these lines in the preamble of your .tex:\n"
		.   "% \t\\documentclass{aastex}\n"
		.   "% or:\n"
		.   "% \t\\usepackage{deluxetable}\n"
		.   "% You may download deluxetable.sty fron, e.g.:\n"
		.   "%\t http://casa.colorado.edu/~danforth/comp/tex/deluxetable.sty\n"
		.   "%\n";
	print	"\t See comments at the top of $TeXfile for\n"
	  .	"\t instructions on how to use the deluxe table.\n";
}

sub print_fontSize_help
{
	my $fontSizeDef = shift;

	print	"\t-f=<f>  :\t font size, one of: 'n' (normalsize), 's' (small),\n"
	 .	"\t         \t 'f' (footnotesize), 'r' (scriptsize) [OPTIONAL; DEFAULT: $fontSizeDef]\n";
}

sub check_fontSize
{
	my( $fontSize, $fontSizeDef ) = @_;

	if( not defined( $fontSize ) )
	{
		$fontSize = $fontSizeDef;
	}
	elsif( $fontSize !~ m/(n|s|f|r)/ )
	{
		die "Error: Invalid font size:\t $fontSize\n";
	}

	return	$fontSize;
}

sub getFontSize
{
	my $fontsize = shift;
	my $font_cmd;
	if( $fontsize eq 'n' )
	{
		$font_cmd = '\normalsize';
	}
	elsif( $fontsize eq 's' )
	{
		$font_cmd = '\small';
	}
	elsif( $fontsize eq 'f' )
	{
		$font_cmd = '\footnotesize';
	}
	elsif( $fontsize eq 'r' )
	{
		$font_cmd = '\scriptsize';
	}
	return	$font_cmd;
}

sub print_numRows_help
{
	my $numRowsDef = shift;
	print	"\t-nr=<nr>:\t number of rows per printed table [OPTIONAL; DEFAULT: $numRowsDef]\n";
}

sub check_numRows
{
	my( $numRows, $numRowsDef ) = @_;

	if( not defined( $numRows ) )
	{
		$numRows = $numRowsDef;
	}
	elsif( $numRows <= 0 )
	{
		die "Error: Invalid number of rows per page:\t $numRows\n";
	}

	return	$numRows;
}

sub print_exampleTeX_help
{
	print	"\t-e      :\t create an example TeX file to test the TeX output,\n";
	print	"\t         \t a script to run pdflatex, and a script to clean up afterwards\n";
}

sub create_run_script
{
	my( $fname, $hasBiblio ) = @_;

	my $script_name = "mktable-$fname.sh";
	open FILE, "> $script_name"
	  or die "Error: Could not open:\t $script_name\n";
	print FILE "pdflatex $fname\n";
	if( defined( $hasBiblio ) )
	{
		print FILE '[ $? == 0 ] && '. "bibtex   $fname\n";
		print FILE '[ $? == 0 ] && '. "pdflatex $fname\n";
		print FILE '[ $? == 0 ] && '. "pdflatex $fname\n";
		print FILE '[ $? == 0 ] && '. "pdflatex $fname\n";
	}
	else
	{
		print FILE '[ $? == 0 ] && '. "pdflatex $fname\n";
	}
	close FILE
	   or die "Error: Could not close:\t $script_name\n";

	print "Created example script:\t $script_name\n";

	return	$script_name;
}

sub print_instructions
{
	my( $FILE, $tableTeX, $fname, $script_name, $comment ) = @_;

	print $FILE "\n";
	print $FILE "# "	if( defined( $comment ) );
	print $FILE "NB NB: to remove intermediate pdflatex files, run:\t\t '. $script_name'\n";
	print $FILE "# "	if( defined( $comment ) );
	print $FILE "NB NB: to also remove ${fname}.tex, and the scripts, run: '. $script_name all'\n";
	print $FILE "# "	if( defined( $comment ) );
	print $FILE "NB NB: ${fname}.pdf, and $tableTeX are not removed by this script\n\n";
}

sub create_cleanup_script
{
	my( $tableTeX, $fname, $hasBiblio, $run_script ) = @_;

	my $script_name = "cleanup-$fname.sh";
	open FILE, "> $script_name"
	  or die "Error: Could not open:\t $script_name\n";

	&print_instructions( \*FILE, $tableTeX, $fname, $script_name, 1 );

	my @suffix = qw/ log aux /;
	push( @suffix, qw/ blg bbl / )	if( $hasBiblio );
	foreach my $suffix ( @suffix )
	{
		my $file = $fname .".". $suffix;
		print FILE "[ -e $file ] && rm $file\n";
	}

	print FILE 'if [ "$1" == "all" ]; then' ."\n";
	print FILE "\trm ". $fname .".tex\n";
	print FILE "\trm $run_script\n";
	print FILE "\trm $script_name\n";	# commit suicide
	print FILE 'fi' ."\n";

	close FILE
	   or die "Error: Could not close:\t $script_name\n";

	print "Created example script:\t $script_name\n";

	&print_instructions( \*STDOUT, $tableTeX, $fname, $script_name );

	return	$script_name;
}

sub create_exampleTeX_Script
{
	my( $deluxeTable, $tableTeX, $hasBiblio ) = @_;

	my $tableTeX_basename = $tableTeX;
	$tableTeX_basename =~ s/\.tex$//;

	my $fname = "$tableTeX_basename-list.tex";
	   $fname = "$tableTeX_basename-deluxe-list.tex"
	   	if( defined( $deluxeTable ) );
	open FILE, "> $fname"
	  or die "Error: Could not open:\t $fname\n";
	if( defined( $deluxeTable ) )
	{
		print	FILE	'\documentclass{aastex}' ."\n";
		if( defined( $hasBiblio ) )
		{
			print	FILE	"\\usepackage{natbib}\n";
			print	FILE	"\\bibliographystyle{apj}\n";
		}
	}
	else
	{
		print	FILE	'\documentclass[usenatbib,onecolumn]{mn2e}' ."\n";
		if( defined( $hasBiblio ) )
		{
			print	FILE	"\\bibliographystyle{mn2e}\n";
		}
	}
	print	FILE	'\usepackage{../docs/latex/common/aas_macros}' ."\n";
	print	FILE	'\begin{document}'	."\n";
	print	FILE	"\\input{$tableTeX_basename}\n";
	if( defined( $hasBiblio ) )
	{
		print	FILE	"\\clearpage\n";
		print	FILE	"\\bibliography{../docs/latex/common/bibliography2}\n";
	}
	print	FILE	'\end{document}'	."\n";
	close FILE
	   or die "Error: Could not close:\t $fname\n";

	print "\nCreated example TeX:\t $fname\n";
	$fname =~ s/\.tex$//;

	my $run_script = &create_run_script( $fname, $hasBiblio );
	&create_cleanup_script( $tableTeX, $fname, $hasBiblio, $run_script );
}

1;
