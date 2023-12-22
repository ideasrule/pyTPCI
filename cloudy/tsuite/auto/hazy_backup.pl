#!/usr/bin/perl

$out_dir = "c:/projects/cloudy/trunk/tsuite/auto/";
$hazy_dir= "c:/projects/";

# move to the cloudy current directory
if( !chdir( $hazy_dir ) )
{
	printf(" invalid directory for hazy\n");
	printf(" was ==%s==\n",$hazy_dir );
	if( $lgLOGGING ) 
	{
		printf(ioLOG " invalid directory for hazy\n");
		printf(ioLOG " was ==%s==\n",$hazy_dir );
	}
	exit(1);
}
# system("ls $hazy_dir");
# make backup copy of hazy
system( 
"tar -cvf \"$out_dir\"\"hazy.tar\" hazy/*.doc hazy/*.jnb hazy/*.spw ");

