#!/usr/bin/perl
# program to determine a list of all parser option keys
foreach $f (<*.cpp>) {
    open FILE,$f or die "Couldn't open $f";
    while (<FILE>) {
	if (/p\.(nMatch|GetParam|GetRange)/) {
	    @l = split(/(?<!\\)\"/);
	    for ($i=1; $i<$#l; $i += 2) {
		$k{$l[$i]} = 1;
	    }
	}
    }
    close FILE;
}
foreach $word (sort keys %k)
{
    print "\"$word\"\n";
}
print "\n";
print <<EOT
Note: there may be accidental additions for strings appearing on the same 
line as command keys: it doesn't seem worth the complexity of coding up 
these cases at present. 
EOT
