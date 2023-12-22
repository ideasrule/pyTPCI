#!/usr/bin/perl
$ARGV[0] =~ s/^\.\///;
@ll = split( /\./, "$ARGV[0]" );
if( $#ll > 0 ) {
    $ext = "$ll[$#ll]";
}
else {
    $ext = "";
}
# filenames with an '@' cause svn errors... add trailing '@' as a workaround
if( $ARGV[0] =~ /@/ ) {
    $arg = $ARGV[0] . "@";
}
else {
    $arg = $ARGV[0];
}
if( ! -d $ARGV[0] ) {
    $out = `svn proplist -v $arg 2>&1`;
    $filetype = `file -b $ARGV[0]`;
    $ascii = ( $filetype =~ /text/ || $filetype =~ /program/ || $filetype =~ /Palm OS/ );
    $script = ( $filetype =~ /ASCII text executable/ );
    $binary = ( $filetype =~ /PDF document/ || $filetype =~ /Excel/ || $filetype =~ /image/i ||
		$filetype =~ /swap file/ || $filetype =~ /CDF/ || $filetype =~ /Composite Document/ ||
		$filetype =~ /PowerPoint/ );
    $symlink = ( $filetype =~ /symbolic link/ );
    $empty = ( $filetype =~ /empty/ );
}
else {
    $out = "directory";
}
if( $out !~ /not under version control/ && $out !~ /not a working copy/ && $out !~ /directory/ ) {
#     the test for $script needs to come first since they also have $ascii set...
    if( $script ) {
	if( $out !~ /svn:eol-style\n    native/ || $out !~ /svn:executable/ ) {
	    &print_out( "executable script", $out );
	}

	&test_eol( $ARGV[0] );
    }
    elsif( $ascii ) {
	if( $out !~ /svn:eol-style\n    native/ || $out =~ /svn:executable/ ) {
	    &print_out( "plain text file", $out );
	}

	&test_eol( $ARGV[0] );
    }
    elsif( $binary ) {
	if( $out !~ /svn:mime-type\n    application\/octet-stream/ || $out =~ /svn:executable/ ) {
	    &print_out( "binary file", $out );
	}
    }
    elsif( $symlink ) {
	if( $out !~ /svn:special/ ) {
	    &print_out( "symbolic link", $out );
	}
    }
    elsif( $empty ) {
	print "---------------------------\n";
	print "$ARGV[0]:\n File is empty.\n"
    }
    else {
	print "---------------------------\n";
	print "WARNING -- file $ARGV[0]: file type not recognized.\n" 
   }
}

sub print_out {
    open FOO, ">>fix_properties";
    print "---------------------------\n";
    $_[1] =~ s/Properties on \'$ARGV[0]\':\n//;
    print "$ARGV[0]: detected file type - $_[0]\n";
    if( $_[0] eq "plain text file" ) {
	if( $_[1] !~ /svn:eol-style\n    native/ ) {
	    print "  property svn:eol-style not set.\n";
	    print FOO "svn propset svn:eol-style native $arg\n";
	}
	if( $_[1] =~ /svn:executable/ ) {
	    print "  property svn:executable should not be set.\n";
	    print FOO "svn propdel svn:executable $arg\n";
	}
    }
    elsif( $_[0] eq "executable script" ) {
	if( $_[1] !~ /svn:eol-style\n    native/ ) {
	    print "  property svn:eol-style not set.\n";
	    print FOO "svn propset svn:eol-style native $arg\n";
	}
	if( $_[1] !~ /svn:executable/ ) {
	    print "  property svn:executable not set.\n";
	    print FOO "svn propset svn:executable \"*\" $arg\n";
	}
    }
    elsif( $_[0] eq "binary file" ) {
	if( $_[1] !~ /svn:mime-type\n    application\/octet-stream/ ) {
	    print "  property svn:mime-type not set.\n";
	    print FOO "svn propset svn:mime-type application/octet-stream $arg\n";
	}
	if( $_[1] =~ /svn:executable/ ) {
	    print "  property svn:executable should not be set.\n";
	    print FOO "svn propdel svn:executable $arg\n";
	}
    }
    elsif( $_[0] eq "symbolic link" ) {
	if( $_[1] !~ /svn:special/ ) {
	    print "  property svn:special not set.\n";
	    print FOO "svn propset svn:special \"*\" $arg\n";
	}
    }
    else {
	print "INTERNAL ERROR -- file type $_[0] not recognized.\n";
    }
    if( $_[1] ne "" ) {
	print "\nThese properties are currently set on \'$ARGV[0]\':\n";
	print $_[1];
    }
    close FOO;
}

sub test_eol {
    open FOO, ">>fix_properties";
    if( ! -z $_[0] ) {
	open( foo, "<$_[0]" );
	sysseek( foo,-1,2 );
	$last = getc( foo );
	close( foo );
	if( "$last" ne "\n" ) {
	    print "$_[0]:\n No newline at end of file\n";
	    print FOO "echo \"\" >> $ARGV[0]\n";
	}
    }
    close FOO;
}
