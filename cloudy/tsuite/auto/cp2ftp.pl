#!/usr/bin/perl -w
# all statements that start with a # are comments in perl
# -w turns on some warnings in above


# this is the perl script that is automatically run every night to verify cloudy
# it ends with sending an email message announcing success or failure - this needs
# to be modifies of deleted if you are not me

# the variable $exe must be set below to point to a valid cloudy executable
# for the OS in use

# bring in the perl module that includes copy
use File::Copy;

# input and output directories
$source_dir = "c:/projects/cloudy/trunk/source/";
$data_dir = "c:/projects/cloudy/trunk/data/";
$tests_dir = "c:/projects/cloudy/trunk/tsuite/auto/";
$ftp_dir = "c:/Inetpub/ftproot/pub/cloudy/";
$hazy_dir = "c:/projects/hazy/";


# ===========================================================================
# data
# remove all files in ftp data subdir
system( 
"rm \"$ftp_dir\"\"data/\"*.*");

# move to the data directory
if( !chdir( $data_dir ) )
{
       printf(" could not move to data directory\n");
       printf(" was ==%s==\n",$data_dir );
       exit(1);
}

# copy all dat files
while( defined( $input = glob("*.dat") ) )
{
   copy( $input , "$ftp_dir"."data/"."$input" );
}  
# copy all htm scripts
while( defined( $input = glob("*.htm") ) )
{
   copy( $input , "$ftp_dir"."data/"."$input" );
}  
# copy all jpg scripts
while( defined( $input = glob("*.jpg") ) )
{
   copy( $input , "$ftp_dir"."data/"."$input" );
}  
# copy all in scripts
while( defined( $input = glob("*.in") ) )
{
   copy( $input , "$ftp_dir"."data/"."$input" );
}  
# copy all ini scripts
while( defined( $input = glob("*.ini") ) )
{
   copy( $input , "$ftp_dir"."data/"."$input" );
}  
# copy all rfi scripts
while( defined( $input = glob("*.rfi") ) )
{
   copy( $input , "$ftp_dir"."data/"."$input" );
}  
# copy all szd scripts
while( defined( $input = glob("*.szd") ) )
{
   copy( $input , "$ftp_dir"."data/"."$input" );
}  
print( "data copy complete\n");

# ===========================================================================
# hazy
# remove all files in ftp hazy subdir
system( 
"rm \"$ftp_dir\"\"hazy/\"*.*");

# move to the main hazy directory
if( !chdir( $hazy_dir ) )
{
       printf(" could not move to hazy directory\n");
       printf(" was ==%s==\n",$hazy_dir );
       exit(1);
}

while( defined( $input = glob("*.pdf") ) )
{
   copy( $input , "$ftp_dir"."hazy/"."$input" );
}  
print( "Hazy copy complete\n");

# ===========================================================================
# tests
# remove all files in ftp tsuite subdir
system( 
"rm \"$ftp_dir\"\"tests/\"*.*");

# move to the tsuite/auto directory
if( !chdir( $tests_dir ) )
{
       printf(" could not move to tests directory\n");
       printf(" was ==%s==\n",$tests_dir );
       exit(1);
}

# copy all input files
while( defined( $input = glob("*.in") ) )
{
   copy( $input , "$ftp_dir"."tests/"."$input" );
}  
# copy all data files
while( defined( $input = glob("*.dat") ) )
{
   copy( $input , "$ftp_dir"."tests/"."$input" );
}  
# copy all perl scripts
while( defined( $input = glob("*.pl") ) )
{
   copy( $input , "$ftp_dir"."tests/"."$input" );
}  
# copy all html scripts
while( defined( $input = glob("*.htm") ) )
{
   copy( $input , "$ftp_dir"."tests/"."$input" );
}  
# copy all jpg scripts
while( defined( $input = glob("*.jpg") ) )
{
   copy( $input , "$ftp_dir"."tests/"."$input" );
}  
print( "tsuite copy complete\n");

# ===========================================================================
# source
# remove all files in ftp source subdir
system( 
"rm \"$ftp_dir\"\"source/\"*.*");

# move to the main source directory
if( !chdir( $source_dir ) )
{
       printf(" could not move to source directory\n");
       printf(" was ==%s==\n",$source_dir );
       exit(1);
}

# copy all files
while( defined( $input = glob("*") ) )
{
   copy( $input , "$ftp_dir"."source/"."$input" );
}  
print( "source copy complete\n");

exit(0);


