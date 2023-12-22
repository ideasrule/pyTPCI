#!/usr/bin/env perl
#
# Run script to remove tests/ directory from Git history.
# 
# The script is expected to operate on a branch other than master.
# It automatically figures out the range of revisions covered by
# the branch, and will operate on them.
#
# NB NB  The first rev reported is the last rev on master prior to
#        branch creation.  This is necessary b/c revision ranges
#        are not lower-bound inclusive under Git.
#
# User confirmation is required at every step of the way.
#
# Git commands are suggested so that the user can independently
# confirm the right operations will be carried out -- best done
# in another terminal.
#
# Chatzikos, Marios			2021-May-25
#


my $tests_dir = "tests/";

sub check_with_user
{
	print "Is this OK? (y/n)\t";
	my $resp = <STDIN>;
	chomp( $resp );
	if( $resp !~ m/^y$/i )
	{
		die "*** Script aborted by user\n";
	}
}

sub get_branch_name
{
	my $branch = `git rev-parse --abbrev-ref HEAD`;
	if( $branch eq "" )
	{
		die "Error: Could not get current branch name\n";
	}
	chomp( $branch );
	if( $branch eq "master" )
	{
		die "Error: Currently on '$branch' (see 'git branch').\n"
		 .  "       Switch to a branch to run script (see 'git checkout').\n";
	}
	print "Operating on branch: $branch\n";
	&check_with_user();

	return $branch;
}

sub get_last_master_commit
{
	my( $branch ) = @_;

	my $output = `git log master..$branch --oneline | wc -l`;
	if( $output eq "" )
	{
		die "Error: Could not get number of commits on branch\n";
	}
	chomp( $output );

	my $nrevs = $output + 1;
	$output = `git log $branch --oneline | head -n $nrevs | tail -n 1`;
	if( $output eq "" )
	{
		die "Error: Could not get last commit before branching\n";
	}
	chomp( $output );
	my( $rev, $message ) = split( /\s+/, $output );

	print "\n";
	print "Will run script for revs: $rev..HEAD (see 'git log --oneline $rev..HEAD').\n";
	&check_with_user();

	--$nrevs;

	return( $nrevs, $rev );
}

sub rm_tests_directory
{
	my( $branch, $rev ) = @_;
	print "\n";
	print "About to delete directory 'tests/' from git $branch history...\n";
	&check_with_user();
	system( "git filter-branch --tree-filter 'rm -rf $tests_dir' $rev..HEAD" ) == 0
	  or die "\n$tests_dir could not be removed from git history\n";
	print "Done!  Removed directory $tests_dir from git history!\n";
}

sub main
{
	print "Notice: You are about to remove directory"
	  .	" $tests_dir from git history\n";
	print "        If anything goes wrong, it's best to start over from a\n"
          .   "        fresh clone of the repository\n\n";

	my $branch = &get_branch_name();
	my( $nrevs, $rev ) = &get_last_master_commit( $branch );
	&rm_tests_directory( $branch, $rev );

	print
	  "\nTip: "
	  .  "Consider dropping empty commits with:"
	  .  " 'git rebase -i HEAD~$nrevs'\n";
}

&main();
