#!/usr/bin/perl -w
# bring in the perl module that includes copy
use File::Copy;

# --------------------------------------------------------------------------
# this block - set up paths to all targets on various machines 

# this is where the original source files live on my pc
$original_programs_dir = "c:/projects/cloudy/current/tsuite/programs";

# the root dir where the unix files will live on my pc
$local_unix_programs_dir = "c:/projects/cloudy/current/tsuite/programs/unix";

# the win dir within the unix set
$local_win_programs_dir = "c:/projects/cloudy/current/tsuite/programs/windows"; 

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# this block - clean out the target directories so we get fresh copies
# move to the old unix source directory, where we will delete all files
#
# dir that are cleaned out:
#  unix source
#  unix tsuite
#
# --------
#
# move to pc unix source dir
if( !chdir( "$local_unix_programs_dir" ) )
{
       printf(" invalid pc unix target directory for programs\n");
       exit(1);
}
print(" moved to pc unix target source directory, about to delete all files\n");
# remove the unix source files
unlink <*> ; 

# move to win source dir
if( !chdir( "$local_win_programs_dir" ) )
{
       printf(" invalid win target directory for programs\n");
       exit(1);
}
print(" moved to win target source directory, about to delete all files\n");
# remove the win source files
unlink <*> ; 

# end block - all target directories cleaned out the so we get fresh copies

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
#
# this block, copy source to various places
#
# move to where the original source files live 
if( !chdir( "$original_programs_dir" ) )
{
       printf(" invalid original target directory for source\n");
       exit(1);
}

# first copy the gif and jpeg files without changes 
copy( "hazy_coolingcurve.gif" , $local_unix_programs_dir );
copy( "hazy_coolingcurve.gif" , $local_win_programs_dir );
copy( "hazy_kmt.gif" , $local_unix_programs_dir );
copy( "hazy_kmt.gif" , $local_win_programs_dir );
copy( "hazy_hizlte.gif" , $local_unix_programs_dir );
copy( "hazy_hizlte.gif" , $local_win_programs_dir );
copy( "varyn.gif" , $local_unix_programs_dir );
copy( "varyn.gif" , $local_win_programs_dir );
copy( "hizlte.gif" , $local_unix_programs_dir );
copy( "hizlte.gif" , $local_win_programs_dir );
copy( "vary_nete.jpg" , $local_unix_programs_dir );
copy( "vary_nete.jpg" , $local_win_programs_dir );
copy( "clouds.jpg" , $local_unix_programs_dir );
copy( "clouds.jpg" , $local_win_programs_dir );

# now copy readme file to each
$output = "$local_unix_programs_dir"."/"."readme_programs.htm";
system "dos2unix < \"readme_programs.htm\"  > $output ";
copy( "readme_programs.htm" , $local_win_programs_dir );

# now copy each source to the win and unix dirs
$output = "$local_unix_programs_dir"."/"."hazy_coolingcurve.c";
copy( "hazy_coolingcurve/hazy_coolingcurve.c" , $local_win_programs_dir );
system "dos2unix < \"hazy_coolingcurve/hazy_coolingcurve.c\"  > $output ";

# now copy each source to the win and unix dirs
$output = "$local_unix_programs_dir"."/"."vary_nete.c";
copy( "vary_nete/vary_nete.c" , $local_win_programs_dir );
system "dos2unix < \"hazy_coolingcurve/hazy_coolingcurve.c\"  > $output ";

# now copy each source to the win and unix dirs
$output = "$local_unix_programs_dir"."/"."mpi.c";
copy( "mpi/mpi.c" , $local_win_programs_dir );
system "dos2unix < \"mpi/mpi.c\"  > $output ";


# now zip up the files we created
if( !chdir( "$original_programs_dir" ) )
{
       printf(" invalid unix home directory for source\n");
       exit(1);
}
print(" moved to programs home directory\n");

# remove old gzipd file if it exists
if( -e "programs.tar.gz" )
{
   unlink("programs.tar.gz" );
}

# remove old gzipd file if it exists
if( -e "programs.zip" )
{
   unlink("programs.zip" );
}

# ----------------------------------------------------------
#   tar the files
#
#  unix source
#
# tar pc unix source
if( !chdir( "$local_unix_programs_dir" ) )
{
       printf(" invalid directory for local unix programs\n");
       exit(1);
}
$tarfile = "../"."programs.tar";
system "tar -cf $tarfile * ";
print( "about to gzip $tarfile\n" );
system "gzip $tarfile ";

# how go to win home and zip the win source dir
#
# move to programs home directory
if( !chdir( "$original_programs_dir" ) )
{
       printf(" invalid directory for win programs\n");
       exit(1);
}

# win source
system "wzzip programs.zip $local_win_programs_dir ";

print("\n=========================\n");

