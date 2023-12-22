#!/usr/bin/perl

# this is the file produced by checkall.pl containing all the botched monitors
$assert_file = "serious.asr";
# this is the file produced by checkall.pl containing all the botched performance monitors
$performance_file = "performance.asr";
# this is the file produced by checkall.pl containing all close botches
$close_file = "close.txt";

if( ! -r $assert_file )
{
    die "$assert_file not found, bailing out\n";
}

if( ! -r $performance_file )
{
    die "$performance_file not found, bailing out\n";
}

if( ! -r $close_file )
{
    die "$close_file not found, bailing out\n";
}

print "please give a short text to be appended to the new comments in the input scripts: ";
$comment = <STDIN>;
chomp( $comment );

# get the current date
@DD = localtime;
my @abbr = qw( jan feb mar apr may jun jul aug sep oct nov dec );
$DD[5] %= 100;
$datestr = sprintf( "%2.2d %s %2.2d", $DD[5], $abbr[$DD[4]], $DD[3] );

print "\ndeleting backup scripts (*.bak, if present)...\n";
while( defined( $backupfile = glob("*.bak") ) )
{
    unlink "$backupfile";
}

$answer = "";
$oldscript = "";
$finalmessage = "\nthe following files have been modified:\n";
$finalmessage2 = "\nTHE FOLLOWING FILES HAD FAILURES:\n";
$finalmessage3 = "\nthese botches were not processed:\n";

open( foo, "<$assert_file" );
&process_file;
close( foo );

open( foo, "<$performance_file" );
&process_file;
close( foo );

open( foo, "<$close_file" );
&process_file;
close( foo );

&print_final_messages;

sub process_file
{
  INPUTLINE:
    while( <foo> )
    {
	chomp;
	$lineorg = $_;
	$line = $_;
	print "\n$line\n";
      GETANSWER:
	while( ! ( $answer =~ /^a.*/i ) )
	{
	    print "\nWould you like to fix this botch [(y)es, (n)o, ,(a)ll, (q)uit]: ";
	    $answer = <STDIN>;
	    chomp( $answer );
	    if( ( $answer eq "" ) || ( $answer =~ /^y.*/i ) || ( $answer =~ /^a.*/i ) )
	    {
		last GETANSWER;
	    }
	    elsif( $answer =~ /^n.*/i )
	    {
		next INPUTLINE;
	    }
	    elsif( $answer =~ /^q.*/i )
	    {
		&clean_up;
		&print_final_messages;
		exit 1;
	    }
	    else
	    {
		print "response not understood: $answer\n";
	    }
	}
#         first extract the output file name and form an input file name from that...
	$newscript = $line;
	$newscript =~ s/:.*//;
	$newscript =~ s/\.out/.in/;
	$newscripttmp = $newscript . ".tmp";
	$newscriptbak = $newscript . ".bak";
#         do we have to open a different input script?
	if( $oldscript ne $newscript )
	{
#             finish copying the previous script
	    &clean_up;
#             then open the next script ($newscript) and its modified version ($newscripttmp)
	    open( old, "<$newscript" );
	    open( new, ">$newscripttmp" );
	    $oldscript = $newscript;
	    $oldscripttmp = $newscripttmp;
	    $oldscriptbak = $newscriptbak;
	}
#         now search the line for the line label and the old and new monitored value
	if( $line =~ /botch>>/ )
	{
	    $line =~ s/.*botch>>//;
	}
	elsif( $line =~ /ChkMonitor -+/ )
	{
	    $line =~ s/.*ChkMonitor -+ *//;
	}
	$lower_limit = 0;
	$upper_limit = 0;
	@fields = split( / +/, $line );
      FIELD:
	for( $i=0; $i <= $#fields; ++$i )
	{
	    if( $fields[$i] eq "=" || $fields[$i] eq "<" || $fields[$i] eq ">" )
	    {
		$oldval = $fields[$i+1];
		$newval = $fields[$i-1];
		$top = $i-1;
		$lower_limit = $fields[$i] eq ">";
		$upper_limit = $fields[$i] eq "<";
		last FIELD;
	    }
	}
#         build the regular expression to match the monitor command in the input script
	$regexp = "^moni.*";
	for( $j=0; $j < $top; ++$j )
	{
#             this does a lot of magic to convert the monitor output back to an input command
	    &do_field_heuristics( $fields[$j], $line );
	    if( $fields[$j] ne "" )
	    {
		$regexp .= "$fields[$j].*";
	    }
	}
	$oldvaltrim = $oldval;
	&do_field_heuristics( $oldvaltrim, $line );
	$regexp .= "${oldvaltrim}[\\.0]*";
	if( $lower_limit || $upper_limit )
	{
	    $oldval = sprintf( "%.4g", $oldval );
	    if( $upper_limit )
	    {
#                if an upper limit is botched we need to increase the new value
		$newval = 1.05*$newval;
		$newval = sprintf( "%.1f", $newval );
		print "Please enter the new upper limit [$newval]: ";
	    }
	    if( $lower_limit )
	    {
		$newval = 0.95*$newval;
		$newval = sprintf( "%.1f", $newval );
		print "Please enter the new lower limit [$newval]: ";
	    }
	    $answer2 = <STDIN>;
	    chomp( $answer2 );
	    if( $answer2 ne "" )
	    {
		$newval = $answer2;
	    }
	}
	else
	{
#             round both the old and new monitored value
	    $oldval = sprintf( "%.4g", $oldval );
	    $newval = sprintf( "%.4g", $newval );
	}
#         now copy the script until we have found the correct monitor command
	$done = 0;
      SCRIPTLINE:
	while( <old> )
	{
	    if( /$regexp/i )
	    {
#                 the correct monitor command has been found; insert comment
		print new "## >>chng $datestr, from $oldval to $newval, $comment\n";
#                 ... and modify the monitor command itself
		$oldvalmatch = "\\b${oldvaltrim}[\\.0]*\\b";
		if( $oldvalmatch =~ /^\\b-/ )
		{
#                     can't use '\b' at the start if ${oldvaltrim} starts with '-'...
		    $oldvalmatch =~ s/^\\b-/-\\b/;
		}
		s/$oldvalmatch/$newval/;
#                 protect '+' sign in case exponential notation was used
		$newval =~ s/\+/\\+/;
#                 now double-check if the change has been done
		if( /$newval/ )
		{
		    $done = 1;
		}
		print new "$_";
#                 stop copying now and determine first if we need to change another monitor command
		last SCRIPTLINE;
	    }
#             this line is not the monitor command we were looking for -> do a verbatim copy
	    print new "$_";
	}
	if( ! $done )
	{
	    print "\nPROBLEM the following line could not be processed:\n";
	    print "$lineorg\n";
	    print "the regular expression was \"$regexp\", the input script was $newscript\n";
	    print "please fix this monitor manually, or repair this script\n";
	    if( ! ( $finalmessage2 =~ /$oldscript/ ) )
	    {
		$finalmessage2 .= "$oldscript\n";
	    }
	    $finalmessage3 .= "$lineorg\n";
#             the failure means that we are now at the end of the file, so clean up and start over
	    &clean_up;
	    open( old, "<$oldscript" );
	    open( new, ">$oldscripttmp" );
	}
    }
    &clean_up;
}

sub clean_up
{
#      copy remainder of the current script, it is not complete yet...
    if( $oldscript ne "" )
    {
	while( <old> )
	{
	    print new "$_";
	}
	close( old );
	close( new );
#         don't overwrite an existing backup script. this e.g. happens when a sim
#         has both full and close botches. when the close botches have been processed
#         we want to keep the first backup so that the user can see all the changes
	if( ! -r "$oldscriptbak" )
	{
	    rename "$oldscript", "$oldscriptbak";
	}
#         this will clobber $oldscript if it still exists
	rename "$oldscripttmp", "$oldscript";
	if( ! ( $finalmessage =~ /$oldscript/ ) )
	{
	    $finalmessage .= "$oldscript\n";
	}
    }
}

sub do_field_heuristics
{
#     modify the monitor output so that the input command can be recovered
#     remove trailing A from wavelength and provide less accurate alternatives for wavelength
    if( $_[0] =~ /^[+-]?\d*\.?\d*[Amc]$/ )
    {
	$num = $_[0];
#         remove unit so that we have only the number
	chop( $num );
#         determine precision with which the number was printed
	$pos = index($num, ".");
	$prec = ( $pos >= 0 ) ? length($num) - $pos - 1 : 0;
#         print number with two decimal places less since that is still used very often in scripts
	$prec = ( $prec >= 2 ) ? $prec-2 : 0;
	$num2 = sprintf( "%.*f", $prec, $num );
#	  remove trailing zeros, and possibly a trailing dot
	$num3 = $num;
	$num3 =~ s/0+$//;
	$num3 =~ s/\.$//;
	$unit = substr($_[0], -1);
#         remove 'A' since it is never used in input scripts, but put other units back
	if( $unit eq "A" )
	{
	    $unit = "";
	}
#         use curly brackets here to prevent them from being excaped
#         they will be replaced with parentheses further below
#         this construct allows either the precise or less precise
	$_[0] = "{$num|$num2|$num3}$unit";
    }
#     remove trailing zeros (even if it is the only digit in the field)
    if( $_[0] =~ /0+$/ )
    {
	$_[0] =~ s/0+$//;
    }
#     remove trailing dot from a number, it may not be present in the input
    if( $_[0] =~ /^[+-]?\d+\.$/ )
    {
	$_[0] =~ s/\.$//;
    }
#     protect '+', '.', '[', ']', '(', ')', and '^', they have special meaning in regular expressions
    $_[0] =~ s/\+/\\+/;
    $_[0] =~ s/\./\\./g;
    $_[0] =~ s/\[/\\[/g;
    $_[0] =~ s/\]/\\]/g;
    $_[0] =~ s/\(/\\(/g;
    $_[0] =~ s/\)/\\)/g;
    $_[0] =~ s/\^/\\^/g;
#     replace curly brackets with parentheses
    $_[0] =~ s/{/(/;
    $_[0] =~ s/}/)/;
#     now do some straight substitutions...
    if( $_[0] =~ /d Fe/ )
    {
	$_[0] = "feii +departure";
    }
    if( $_[0] =~ /d He/ )
    {
	$_[0] = "he-like.+departure";
    }
    if( $_[0] =~ /dHyd/ )
    {
	$_[0] = "h-like.+departure";
    }
    if( $_[0] =~ /GPot/ )
    {
	$_[0] = "grain +potential";
    }
    if( $_[0] =~ /GTem/ )
    {
	$_[0] = "grain +temperature";
    }
    if( $_[0] =~ /MapH/ )
    {
	$_[0] = "map +heating";
    }
    if( $_[0] =~ /MapC/ )
    {
	$_[0] = "map +cooling";
    }
    if( $_[0] =~ /sion/ )
    {
	$_[0] = "csupra";
    }
    if( $_[1] =~ /^H2/ && $_[0] =~ /^[1-9]0[1-9]$/ )
    {
	$_[0] =~ s/([1-9])0([1-9])/$1 +$2/;
    }
}

sub print_final_messages
{
    if( $finalmessage =~ /\.in\n$/ )
    {
	print "$finalmessage\n";
    }
    if( $finalmessage2 =~ /\.in\n$/ )
    {
	print "$finalmessage2\n";
    }
    if( ! ( $finalmessage3 =~ /processed:\n$/ ) )
    {
	print "$finalmessage3\n";
    }
}
