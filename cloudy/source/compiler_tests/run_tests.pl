#!/usr/bin/perl
$logfile = "run_tests.log";
if( $#ARGV == 0 ) {
    $CXX = $ARGV[0];
}
else {
    $CXX = "g++";
}
print "running the test suite using the $CXX compiler\n";
print "output will be in $logfile\n\n";
$test = 0;
$xfail = 0;
$pass = 0;
unlink $logfile;
system( "touch $logfile" );
while( defined( $input = glob("test*.cpp") ) ) {
    $res = system( "$CXX -O0 -c -I.. $input >> $logfile 2>&1" );
    if( ! $res ) {
	print "\nunexpected pass on $input\n";
	$pass += 1;
    }
    else {
	$xfail += 1;
    }
    $test += 1;
    print ".";
}
print "\n\n";

print "number of tests executed: $test\n";
if( $xfail > 0 ) {
    print "number of expected failures: $xfail\n";
}
if( $pass > 0 ) {
    print "number of unexpected passes: $pass\n";
}
if( $test == $xfail ) {
    print "you have a valid implementation!\n";
}
