#!/usr/bin/perl

# This script tests for identical output among individual sims in "grid repeat" or "grid cycle" runs.

$out = $ARGV[0];
$prefix = $ARGV[1];
$stride = $ARGV[2];
$delimiter = "GRID_DELIMIT -- grid";

# grep "Cloudy was called" and grab the number of times
$num_sims = qx/grep "$delimiter" $out | wc -l/;
chomp($num_sims);

# delete any old split files
while( defined( $input = glob("$prefix??") ) )
{
    unlink "$input";
}

# get rid of "duplicate reaction", "reverse reaction" warnings and "Using STOUT model" messages
# strip the "xxxxxxxxx" from "GRID_DELIMIT -- gridxxxxxxxxx"
system("egrep -v 'Warning: duplicate species|Warning: duplicate reaction|Warning! No reverse reaction|Using STOUT model' $out | sed 's/grid[0-9]\\{9\\}/grid/' > tmpfile");
# split output file and count number of resultant files
# should be number of sims plus 2 for header and footer
$command = "./split.pl tmpfile $prefix GRID_DELIMIT";
my $num_split = qx/$command/;
unlink "tmpfile";

# now substract the 2
$num_split -= 2;

# make sure $num_sims and $num_split agree
if( $num_sims != $num_split )
{
    die "PROBLEM These should agree: ".$num_sims." and ".$num_split."\n";
}

unlink $prefix . "00";
# generate filename for footer
$footer = $prefix . sprintf( "%02i",$num_sims+1 );
unlink "$footer";

$res = 0;
for( $i=1; $i <= $stride; ++$i )
{
    # generate the filename
    $sim1 = $prefix . sprintf( "%02i", $i );
    for( $j=$i+$stride; $j <= $num_sims; $j += $stride )
    {
	$sim2 = $prefix . sprintf( "%02i", $j );
	# Compare everything to the last sim.
	# This minimizes output if first is different from all others which are themselves identical.
	# That is the most common pathology.
	$res += system( "diff -q $sim1 $sim2" );
    }
}
if( $res == 0 )
{
    print "all sims are equal\n";
}
else
{
    print "PROBLEM: some sims did not repeat exactly\n";
}
