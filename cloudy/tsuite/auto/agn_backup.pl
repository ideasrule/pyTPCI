#!/usr/bin/perl

$out_dir = "c:/projects/cloudy/trunk/tsuite/auto/";
$agn_dir= "c:/projects/";

# move to the cloudy current directory
if( !chdir( $agn_dir ) )
{
	printf(" invalid directory for agn\n");
	printf(" was ==%s==\n",$agn_dir );
	if( $lgLOGGING ) 
	{
		printf(ioLOG " invalid directory for agn\n");
		printf(ioLOG " was ==%s==\n",$agn_dir );
	}
	exit(1);
}
# system("ls $agn_dir");
# make backup copy of agn word docs
system( 
"tar -cvf \"$out_dir\"\"agn.tar\" agn/*.doc agn/*/*.doc ");

