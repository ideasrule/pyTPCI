#!/usr/bin/perl
$update = 1;
foreach $f (<*.h>)
{    
    $new = $f.'.new';
    $old = $f.'.old';
    open FILE, $f;
    open NEW, ">$new" if $update;
    $guard = $f.'_';
    $guard =~ tr/a-z\./A-Z_/;
    $oldguard = '_'+$guard;
    $state = 0;
    $chg = 0;
    foreach $l (<FILE>)
    {
	$o = $l;
	if ($state == 0 && $l =~ /\#ifndef\s+(\w+)/) {
	    # print "Found $f $guard\n$l";
	    $1 == $guard or $1 == $oldguard or
		die "Current guard does not match in $f\n$l";
	    $l = "#ifndef $guard\n";
	    $state = 1;
	} elsif ($state == 1 && $l =~ /\#define\s+(\w+)/) {
	    # print "Defined $f $guard\n";
	    $1 == $guard or $1 == $oldguard or
		die "Current guard does not match in $f\n$l";
	    $l = "#define $guard\n";
	    $state = 2;
	} elsif ($state == 2 && $l =~ /\#endif\s+(\/\*\s*\w+\s*\*\/)/) {
	    if ($1 == $oldguard or $1 == $guard) {
		# print "Closed $l -- $f $guard\n";
		$state = 3;
		$l = "#endif /* $guard */\n";
	    }
	} elsif ($state == 3 && ! $l =~ /^\n/ ) {
	    die "Found trailing matter or unlabelled endif in $f\n$l";
	}
	print NEW $l if $update;
	$chg = 1 if not $chg and $o ne $l;
    }
    close FILE;
    close NEW;
    if ($update) {	
	if ($chg) {
	    rename $f, $old;
	    rename $new, $f;
	    print "$f changed\n";
	} else {
	    unlink $new;
	}
    }
}
