#!/usr/bin/perl

# This script sends an email to the address in ~/.forward, if any.

use Cwd;
$curdir = getcwd;

$output = "A script has finished. It lives here: " . $curdir . ".  Time is " . `date`;
chomp($output);

$forwardfile="$ENV{HOME}/.forward";
print "Checking if $forwardfile exists....\n";
if ( (-e $forwardfile) && (-s $forwardfile) )
{
    print "Forwarding address exists and is non-empty.\n";
    open (ff, $forwardfile) || die ("Could not open file. <br> $!");
    $email = <ff>;
    chomp($email);
    print "Trying to send mail to " . $email . ".\n";
    # now send mail    
    system "echo $output | mail -s \"run_parallel.pl has finished.\" $email";
    print "Did it work?\n";
}
else
{
    print "Could not find an email address.\n";
}

