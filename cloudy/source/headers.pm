#################################################################################
# headers.pm:
#    Construct header dependency tree; determine required headers for each header
# to compile cleanly; and provide facilities to check source (.cpp) files for the
# header requirements.
#
# Chatzikos, Sep 23			First Version 
#################################################################################

use strict;
use warnings;

package headers;                                                       
require Exporter;                                                        

use threads;
use threads::shared;


our @ISA        = qw/Exporter/;                                          
our @EXPORT     = qw/ ResolveHdr /;
our @EXPORT_OK  = qw/	do_parallel do_compile header_requirements
			get_contents parse_include get_includes_from_cpp uniq_head_contents
			header_dependencies ResolveHdr getHdrBranch getHdrAlldep getHdrDep
			reorder_headers header_calls
			prtHdrTree prtHdrAlldep prtHdrDep prtHdrBranches /;
our $VERSION    = 2013.0920;



my %AllHeaders;

my (@output) : shared;

my $cpp_extension = ".cpp";
my $compile_flags = "";



#################################################################################
#                                PARALLELISM                                    #
#################################################################################
sub get_thread_lims
{
	my( $ncpus, $files ) = @_;

	my( $files_per_thread, $nleft) = (1, 0);
	if( @$files <= $ncpus )
	{
		$ncpus = @$files;
	}
	else
	{
		$files_per_thread = int( @$files / ($ncpus) );
		$nleft = @$files % $ncpus;
	}
	print "files/thread = $files_per_thread \t nleft= $nleft\n"
		if( 0 );

	my( @ibeg, @iend );
	$ibeg[ 0 ] = 0;
	$iend[ 0 ] = $files_per_thread-1;
	$iend[ 0 ] ++	if( $nleft );
	for( my $ithread = 1; $ithread < $ncpus; $ithread++ )
	{
		my $this_files_per_thread = $files_per_thread;
		   $this_files_per_thread++
			if( $ithread < $nleft );
		$ibeg[ $ithread ] = $iend[ $ithread-1 ] + 1;
		$iend[ $ithread ] = $iend[ $ithread-1 ] + $this_files_per_thread;
		$iend[ $ithread ] = @$files - 1
			if( $iend[ $ithread ] >= @$files );
	}

	return	(\@ibeg, \@iend);
}



sub do_parallel
{
	my ($ncpus, $files, $run_thread) = @_;

	print "files: #". @$files ."\t @$files\n"	if( 0 );

	my ($ibeg, $iend) = &get_thread_lims( $ncpus, $files );

	my @thread;
	for (my $i=0; $i < @$ibeg; $i++)
	{
		my @files_loc = @$files[ $$ibeg[$i] .. $$iend[$i] ];
		print "$i:\t loc_files:\t ".
			$$ibeg[$i] ."--". $$iend[$i] ."\t @files_loc\n"
			if( 0 );
		$thread[$i] = threads->new($run_thread, $i, @files_loc);
	}

	for (my $i=0; $i < @$ibeg; $i++)
	{
		$thread[$i]->join;
	}
	return;
}





#################################################################################
#                                  COMPILATION                                  #
#################################################################################
sub do_compile
{
	my ($cpp, $compile_flags) = @_;
	my $obj = $cpp;
	$obj =~ s/\.cpp$/.o/;

	my $command = "g++ -c -Wall -Wmissing-declarations $cpp $compile_flags -o $obj 2>&1"; 
	print "$command\n"	if( 0 );

	my $out = `$command`;
	print "$cpp:\t$out\n"	if( 0 and $out ne "" and $out ne "\n" );

	if( -f $obj )
	{
		unlink $obj
		    or warn "do_compile:\t Could not delete:\t $obj\n";
	}

	return	$out;
}




sub make_test_cpp
{
	my ($test_cpp, $header, @Reqrd) = @_;

	open FILE, "> $test_cpp"	or die  "make_test_cpp: Error: Could not open:\t $test_cpp\n";
	print FILE "#include \"cddefines.h\"\n";
	foreach my $reqr ( @Reqrd )
	{
		print FILE "#include \"$reqr\"\n";
	}
	print FILE "#include \"$header\"\n";
	close FILE			or warn "make_test_cpp: Warning: Could not open:\t $test_cpp\n";

	return;
}



sub decipher_err
{
	my ( $error ) = @_;
	my @recm;

	my @err = split(/\n/, $error);
	my @undefnd = grep( /declared/, @err );

	foreach my $line ( @undefnd )
	{
		# Emilimate unicode characters for
		# left and right single quotation marks
		#
		$line =~ s/\xE2\x80\x98/\'/;
		$line =~ s/\xE2\x80\x99/\'/;
		my (undef, $var) = split(/\'/, $line);
		my @matches = `grep \"class\\s\\+$var\" *.h`;
		@matches = `grep \"struct\\s\\+$var\" *.h`
			if( @matches == 0 );
		foreach my $match ( @matches )
		{
			my ($head) = split(/:/, $match);
			push( @recm, $head );
		}
	}

	return	&unique( @recm );
}



sub header_req
{
	my ($thread, @heads) = @_;
	my $cpp = "test-$thread$cpp_extension";

	$output[$thread] = "";
	foreach my $header ( @heads ) 
	{
		print "$header\n"
			if( 0 );

		my @Reqrd;
		next
			if( $header eq "physconst_template.h" );

		my $error;

		do {
			&make_test_cpp( $cpp, $header, @Reqrd );
			$error = &do_compile( $cpp, $compile_flags );
			my @recm = &decipher_err( $error );
			push( @Reqrd, @recm );
		} while( length( $error ) > 0 );
	
		if( @Reqrd )
		{
			$output[$thread] .= sprintf( "%22s:\t %s\n", $header, "@Reqrd" );
		}

		unlink $cpp
		    or die "header_req: Error: Could not delete:\t $cpp\n";
	}

	return;
}



sub header_requirements
{
	my( $ncpus ) = @_;
	my @hd = glob "*.h";
	&do_parallel( $ncpus, \@hd, \&header_req );
	my $out = "";
	for(my $i = 0; $i < @output; $i++ ) { $out .= $output[$i]; }
	return	$out;
}





#################################################################################
#                           HEADER TREE CONSTRUCTION                            #
#################################################################################
sub get_contents
{
	my ($file) = @_;
	open FILE, "< $file"	or die "get_contents: Could not open file:\t $file\n";
	my @contents = <FILE>;
	close FILE		or warn "get_contents: Could not close file:\t $file\n";
	return	@contents;
}



sub unique
{
	my @headers = @_;

	my @unique;
	foreach my $header ( @headers )
	{
		my $all_head = join(" ", @unique);
		if( $all_head !~ m/\b\Q$header\E\b/ )
		{
			push( @unique, $header );
		}
	}

	return	@unique;
}



sub parse_include
{
	my ($line) = @_;

	my $header;
	if( $line =~ m/\"/ )
	{
		my @words = split(/\"/, $line);
		$header = $words[1];
	}
	elsif( $line =~ m/\</ )
	{
		my @words;
		@words = split(/\</, $line);
		@words = split(/\>/, $words[1]);
		$header = $words[0];
	}
	return	$header;
}



sub get_includes_from_headers
{
	my ( @contents) = @_;

	my %headers;
	my $isinif = 0;		# account for the guard
	foreach my $line ( @contents )
	{
		if( $isinif == 1 and $line =~ m/^#\s*include/ )
		{
			my $header = &parse_include( $line );
			$headers{$header} = -1;
		}
		else
		{
			$isinif++	if( $line =~ m/^#\s*if/ );
			$isinif--	if( $line =~ m/^#\s*endif/ );
		}
	}

	return	%headers;
}



sub get_includes_from_cpp
{
	my ($skipHdr, @contents) = @_;

	my @headers;
	my $isinif  = 0;		# account for guards
	foreach my $line ( @contents )
	{
		if( $isinif == 0 and $line =~ m/^#\s*include/ )
		{
			my $header = &parse_include( $line );
			push( @headers, $header );
		}
		else
		{
			$isinif++	if( $line =~ m/^#\s*if/ );
			$isinif--	if( $line =~ m/^#\s*endif/ );
		}
	}

	@headers = &unique( @headers );
	foreach my $skip ( @$skipHdr )
	{
		my $i;
		for( $i = 0; $i < @headers; $i++ )
		{
			last
				if( $headers[$i] =~ $skip );
		}
		if( $i != @headers )
		{
			splice( @headers, $i, 1 );
		}
	}

	return	&unique( @headers );
}



sub uniq_head_contents
{
	my ($h, @contents) = @_;
	my @headers = @$h;

	my %headers;
	foreach my $h ( @headers ) { $headers{$h} = 0; }

	my @uniq_contents;
	my $isinif  = 0;		# account for guards
	foreach my $line ( @contents )
	{
		if( $isinif == 0 and $line =~ m/^#\s*include/ )
		{
			my $found = 0;
			foreach my $head ( @headers )
			{
				if( $line =~ m/\b\Q$head\E\b/ )
				{
					++$headers{$head};
					push( @uniq_contents, $line )
						if( $headers{$head} == 1 );
					$found = 1;
					last;
				}
			}
			push( @uniq_contents, $line )
				if( $found == 0 );
		}
		else
		{
			$isinif++	if( $line =~ m/^#\s*if/ );
			$isinif--	if( $line =~ m/^#\s*endif/ );
			push( @uniq_contents, $line );
		}
	}

	my $nmult = 0;
	my @mult_headers;
	foreach my $h ( @headers )
	{
		my $n = $headers{$h} -1;
		push( @mult_headers, $h )
			if( $n > 0 );
		$nmult += $n; 
	}

	return	( $nmult, \@mult_headers, @uniq_contents );
}



sub get_include_statm
{
	my ($h, @contents) = @_;
	my @headers = @$h;

	my %lines;
	foreach my $h ( @headers ) { $lines{$h} = 0; }

	my @uniq_contents;
	my $isinif  = 0;		# account for guards
	foreach my $line ( @contents )
	{
		if( $isinif == 0 and $line =~ m/^#\s*include/ )
		{
			foreach my $head ( @headers )
			{
				if( $line =~ m/\b\Q$head\E\b/ )
				{
					$lines{$head} = $line;
					last;
				}
			}
		}
		else
		{
			$isinif++	if( $line =~ m/^#\s*if/ );
			$isinif--	if( $line =~ m/^#\s*endif/ );
		}
	}

	return	%lines;
}



sub get_last_include_line
{
	my ($cpp, @contents) = @_;

	my $i;
	for( $i = 0; $i < @contents; $i++ )
	{
		last	if( $contents[$i] =~ m/^(static|int|double|realnum|void|namespace|t_)/i );
	}


	my $iline = ( $i >= @contents ? @contents-1 : $i );
	for( ; $iline >= 0; $iline-- )
	{
		last	if( $contents[$iline] =~ m/^#\s*(undef|endif|include)/ );
	}


	return	$iline;
}



sub gather_headers_top
{
	my ($cpp, $last_headers, $inc, @contents) = @_;
	my %includes = %$inc;

	my $last = join( " ", @$last_headers );
	my $iline = &get_last_include_line( $cpp, @contents );

	my @newcontents;

	my $isinif  = 0;
	for( my $i = 0; $i <= $iline; $i++ )
	{
		my $line = $contents[$i];
		if( $isinif == 0 and $line =~ m/^#\s*include/ )
		{
			my $header = &parse_include( $line );
			next	if( $last =~ m/\b\Q$header\E\b/ );
			if( exists $includes{$header} )
			{
				push( @newcontents, $includes{$header} );
			}
			else
			{
				push( @newcontents, $line );
			}
		}
		else
		{
			$isinif++	if( $line =~ m/^#\s*if/ );
			$isinif--	if( $line =~ m/^#\s*endif/ );
			push( @newcontents, $line );
		}
	}

	for( my $i = 0; $i < @$last_headers; $i++ )
	{
		push( @newcontents, $includes{$$last_headers[$i]} );
	}

	for( my $i = $iline+1; $i < @contents; $i++ )
	{
		my $line = $contents[$i];
		push( @newcontents, $line );
	}

	return	@newcontents;
}



sub get_header_req
{
	my ($ncpus) = @_;

	my $output = &header_requirements( $ncpus );
	my @output = split(/\n/, $output);
	foreach my $hreq ( @output )
	{
		my ($header, $req) = split(/\:/, $hreq);
		$header =~ s/\s//g;
		my @req = split(/\s+/, $req);
		shift( @req );
		push( @{ $AllHeaders{$header}{Reqr} }, @req );
	}
	return;
}



sub concat_alldep
{
	my (%calls) = @_;

	return	()
		if( not $calls{Dep} );

	my @alldep;
	push( @alldep, @{ $calls{Dep} } );

	foreach my $dep ( keys %calls )
	{
		next	if( $dep =~ m/(AllDep|Dep|Reqr)/ );
		if( $calls{$dep} )
		{
			my @this_dep = &concat_alldep( %{ $calls{$dep} } );
			push( @alldep, @this_dep )
				if( @this_dep );
		}
	}

	return	@alldep;
}



sub header_dependencies
{
	foreach my $header ( glob "*.h" )
	{
		my @contents = &get_contents( $header );
		my %headers = &get_includes_from_headers( @contents );
		$AllHeaders{$header}{Dep} = [ keys %headers ];
	}


	my @order = sort{ scalar($AllHeaders{$a}{Dep}) <=> scalar($AllHeaders{$b}{Dep})}	keys %AllHeaders;
	#	print @order."\n";

	foreach my $header ( @order )
	{
		#	print $header."\n";
		foreach my $called ( @{ $AllHeaders{$header}{Dep} } )
		{
			if( $called =~ m/\.h$/ and $AllHeaders{$called}{Dep} )
			{
				$AllHeaders{$header}{$called} = $AllHeaders{$called};
			}
			else
			{
				$AllHeaders{$header}{$called} = ();
			}
		}
	}
	#	foreach my $header ( @order ) { print "$header\n"; }
	#	@order = keys %AllHeaders; print @order."\n";
	#	die;


	foreach my $header ( @order )
	{
		my @alldep = &concat_alldep( %{ $AllHeaders{$header} } );
		@alldep = &unique( @alldep );
		push( @{ $AllHeaders{$header}{AllDep} }, @alldep );
	}

	return;
}





#################################################################################
#                            ACCESS TREE BRANCHES                               #
#################################################################################
sub getHdrBranch
{
	my ($header) = @_;
	my %this = ();
	if( $AllHeaders{$header} )
	{
		%this = $AllHeaders{$header};
	}
	return	%this;
}



sub getHdrAlldep
{
	my ($header) = @_;
	my @this = ();
	if( $AllHeaders{$header} )
	{
		@this = $AllHeaders{$header}{AllDep};
	}
	return	@this;
}



sub getHdrDep
{
	my ($header) = @_;
	my @this = ();
	if( $AllHeaders{$header} )
	{
		@this = $AllHeaders{$header}{Dep};
	}
	return	@this;
}



sub reorder_headers
{
	my ($cpp, $headers, @contents) = @_;
	my (@put_first, @put_last);

	my %headers = %$headers;
	my @headers = keys %$headers;

	foreach my $head ( @headers )
	{
		if( not $AllHeaders{$head}{Reqr} )
		{
			push( @put_first, $head );
		}
		else
		{
			push( @put_last, $head );
		}
	}

	my @add;
	foreach my $head ( @put_last )
	{
		foreach my $req ( @{ $AllHeaders{$head}{Reqr} } )
		{
			foreach my $h ( @headers )
			{
				last
					if( not $AllHeaders{$h}{AllDep} );
				foreach my $ih ( @{ $AllHeaders{$h}{AllDep} } )
				{
					if( $ih eq $req )
					{
						push( @add, $h );
						last;
					}
				}
			}
		}
	}

	unshift( @put_last, @add );
	@put_last = &unique( @put_last );

	my @reorder = ( @put_first, @put_last );
	@reorder = &unique( @reorder );

	my %incl_statm = &get_include_statm( \@headers, @contents );
	@contents = &gather_headers_top( $cpp, \@put_last, \%incl_statm, @contents );

	return	( \@reorder, @contents );
}



sub header_calls 
{
	my @headers = @_;

	my %included;
	foreach my $header ( @headers )
	{
		$included{$header}{calls}  = ();
		$included{$header}{called} = ();
	}


	foreach my $header ( @headers )
	{
		my @calls = ();
		if( $AllHeaders{$header}{AllDep} )
		{
			my @a;
			push( @a, @{ $AllHeaders{$header}{AllDep} } );
			foreach my $hdep ( @a )
			{
				foreach my $hincl ( @headers )
				{
					next	if( $hincl eq $header );
					if( $hincl eq $hdep )
					{
						push( @calls, $hdep );
						last;
					}
				}
			}
			push( @{ $included{$header}{calls} }, @calls )
				if( scalar( @calls ) );
		}
	}


	my (@indep, @dep);
	foreach my $header ( @headers )
	{
		my @depnd;
		foreach my $hdep ( @headers )
		{
			next	if( $header eq $hdep or not $AllHeaders{$hdep}{AllDep} );
			my $adep = join( " ", @{ $AllHeaders{$hdep}{AllDep} } );
			if( $adep =~ m/\b\Q$header\E\b/ )
			{
				push( @depnd, $hdep );
				push( @{ $included{$header}{called} }, $hdep );
			}
		}
	}
	@indep = &unique( @indep );
	@dep   = &unique( @dep   );

	#	foreach my $dep (@dep)
	#	{
	#		print	"$dep\n";
	#		my @a = @{ $included{$dep}{called} };
	#		foreach my $adep (@a) { print "\t$adep"; }
	#		print "\n";
	#	}

	return	%included;
}





#################################################################################
#                                 PRINTING                                      #
#################################################################################
sub getNdep
{
	my @deps = @_;
	my $ndep = 0;
	foreach my $dep ( @deps )
	{
		next	if( $dep =~ m/(AllDep|Dep|Reqr)/ );
		$ndep++;
	}
	return	$ndep;
}



my @vert_bars;
my $ibar_min = 0;
my $connector = "  |_____";
sub prt_tree_lines
{
	my ($nvert) = @_;
	for( my $i = $ibar_min; $i < $nvert; $i++ ) 
	{
		print	"  |"
			if( $vert_bars[$i] > 0 );
		print	"\t";
	}
	print $connector;
}



sub prtSubHdr
{
	my ($ntab, %called) = @_;

	$vert_bars[$ntab] = &getNdep( keys %called ) -1;
	#	print $ntab ."\t" . $vert_bars[$ntab] ."\n";

	foreach my $callee ( keys %called )
	{
		next	if( $callee =~ m/(AllDep|Dep|Reqr)/ );
		&prt_tree_lines( $ntab );
		print "$callee\n";
		if( $called{$callee} )
		{
			my %this_dep = %{ $called{$callee} };
			&prtSubHdr( $ntab+1, %this_dep );
		}
		--$vert_bars[$ntab];
	}

	return;
}



sub prtHdrBranches
{
	my (@headers) = @_;

	if( not @headers )
	{
		&prtHdrTree();
	}
	else
	{
		$vert_bars[0] = 0;
		$ibar_min = 1;
		for( my $i = 0; $i < @headers; $i++ )
		{
			my $header = $headers[$i];
			print	"$header\n";
			&prtSubHdr( 1, %{ $AllHeaders{$header} } );
			print	"#" x 25 ."\n"
				if( $i != @headers-1 );
		}
	}

	return;
}



sub prtTree
{
	my @order = @_;

	$vert_bars[0] = 1;
	$ibar_min = 0;
	foreach my $header ( @order )
	{
		print	"$connector$header\n";
		if( $AllHeaders{$header}{Dep} )
		{
			#	print &getNdep( keys %{ $AllHeaders{$header} } ) ."\n";
			&prtSubHdr( 1, %{ $AllHeaders{$header} } );
		}
	}
	return;
}



sub prtHdrTree
{
	my @order = sort keys %AllHeaders;

	print "Cloudy\n";
	&prtTree( @order );
}



sub prtSrcTree
{
	my( $this_cpp, @headers ) = @_;

	print "$this_cpp\n";
	&prtTree( @headers );
}



sub prtHdrAlldep
{
	my @heads = @_;

	@heads = keys %AllHeaders
		if( not @heads );
	@heads = sort @heads;

	foreach my $header ( @heads )
	{
		next	if( not $AllHeaders{$header}{AllDep} );
		if( @{ $AllHeaders{$header}{AllDep} } )
		{
			my @this = @{ $AllHeaders{$header}{AllDep} };
			my $deps = join(" ", @this);
			printf "%22s:\t %s\n", $header, $deps;
		}
	}
}



sub prtHdrDep
{
	my @heads = @_;

	@heads = keys %AllHeaders
		if( not @heads );
	@heads = sort @heads;

	foreach my $header ( @heads )
	{
		next	if( not defined( $AllHeaders{$header}{Dep} ) );
		if( scalar( @{ $AllHeaders{$header}{Dep} } ) )
		{
			my @this = @{ $AllHeaders{$header}{Dep} };
			my $deps = join(" ", @this);
			printf "%22s:\t %s\n", $header, $deps;
		}
	}
}



sub ResolveHdr
{
	my ($ncpus) = @_;
	&header_dependencies();
	&get_header_req( $ncpus );
}

1;
