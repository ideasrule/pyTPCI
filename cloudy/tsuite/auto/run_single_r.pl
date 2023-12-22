#!/usr/bin/perl

# This is a little helper script for run_parallel.pl that runs a single instance of Cloudy
# The syntax is:
#
# run_single.pl "pwd" "command" "path_to_input_script"
#
# pwd: path to where run_parallel.pl was invoked
# command: the command that should be executed
#      if a relative path is used, it should be relative to where the input script resides!
# path_to_input_script: path to the input script relative to pwd
#
# Peter van Hoof

$command = $ARGV[1];
$path_script = $ARGV[2];

$script = $path_script;
$script =~ s|^.*/||;

chdir "$ARGV[0]";

if( $path_script =~ m|/| ) {
	$wd = $path_script;
	$wd =~ s|/[^/]*$||;
	chdir "$wd";
}

$output2 = $script;
$output2 =~ s/\.in/.err/gi;

$script =~ s/\.in//gi;

system "$command -r $script 2> $output2";
