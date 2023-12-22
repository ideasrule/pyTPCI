#!/usr/bin/perl
# this script checks whether the installed g++ supports certain features we need
$res = "";
$os = `uname -s`;
$arch = "";
if( $os =~ /Darwin/ )
{
    my $tmp = `uname -m`;
    chomp( $tmp );
    $arch = "arm64" if( $tmp =~ m/arm64/ );
}
$version = `$ARGV[0] --version 2> /dev/null`;
# remove comments in parentheses, they may or may not contain spaces
$version =~ s/\(.*?\)//g;
@v1 = split( /\s+/, $version );
$compiler = $v1[0];
# this is needed for CC
$compiler =~ s/://;
# g++ may identify itself as e.g. g++-7
$compiler =~ s/^g\+\+.*/g++/;
# at this point $compiler will hold "g++", "clang", "icc", or "CC" for supported compilers
if( $compiler eq "g++" ) {
    @v2 = split( /\./, $v1[1] );
#     $v2[0] is major revision, $v2[1] is minor revision, $v2[2] is microversion
#     precompiling headers is supported from g++ 3.4.0 onwards
#     it's broken on FreeBSD 11.0-RELEASE, so disable it on all BSD variants to be safe
#     it is also broken when using the openMPI mpiCC wrapper, so test for that as well
    if( $os !~ /bsd/i && $ARGV[0] !~ /mpi/ ) {
	$res .= "precompile ";
    }
#     on modern FreeBSD installations, LLVM is the standard compiler
#     if you use GNU g++, you need different runtime support libraries
    if( $os =~ /FreeBSD/ ) {
	$libdir = "/usr/local/lib/gcc";
	if( $v2[0] <= 4 ) {
	    $libdir .= "$v2[0]$v2[1]";
	}
	else {
	    $libdir .= "$v2[0]";
	}
	if( -d $libdir ) {
	    $res .= "rpath=$libdir ";
	}
    }
#     vectorization is supported from g++ 4.0.0 onwards
    $res .= "vectorize ";
#     ignore AVX capabilities on Darwin since assembler doesn't support it
    if( $os !~ /Darwin/ ) {
#         -march=native is supported from g++ 4.2.0 onwards, however we only
#         use it from 4.6.0 onwards to avoid sub-optimal AVX compilations...
	$res .= "native ";
    }
}
# test if the -Wl,-export-dynamic flag is supported
open FOO, ">tmp_cloudyconfig.cpp";
print FOO <<'EOF';
int main()
{
    return 0;
}
EOF
close FOO;
if( "$^O" ne "cygwin" ) {
    $ret = system( "$ARGV[0] -Wl,-export-dynamic -o tmp_cloudyconfig.out tmp_cloudyconfig.cpp > /dev/null 2>&1" );
    if( $ret == 0 ) {
	$res .= "dynamic ";
    }
}
unlink "tmp_cloudyconfig.out";
if( $arch ne "arm64" )
{
    # test if the -no-pie and -fno-pie flags are supported
    if( $compiler eq "clang" || $os =~ /Darwin/ ) {
        $npflags = "-fno-pie";
    }
    else  {
        $npflags = "-no-pie -fno-pie";
    }
    $ret = system( "$ARGV[0] $npflags -o tmp_cloudyconfig.out tmp_cloudyconfig.cpp > /dev/null 2>&1" );
    if( $ret == 0 ) {
        $res .= "$npflags ";
    }
}
unlink "tmp_cloudyconfig.cpp";
unlink "tmp_cloudyconfig.out";
# test the necessary headers and libraries are installed to demangle routine names
open FOO, ">tmp_cloudyconfig.cpp";
print FOO <<'EOG';
extern "C" int cplus_demangle_noret(const char*, char*, unsigned long);
#include <iostream>
int main()
{
    char buf[1024];
    int p = cplus_demangle_noret( "__1cMParseCrashDo6FrnGParser__v_", buf, sizeof(buf) );
    if( p == 0 )
        std::cout << buf << "\n";
    else
        std::cout << "invalid mangled name!\n";
    return 0;
}
EOG
close FOO;
$ret = system( "$ARGV[0] -o tmp_cloudyconfig.out tmp_cloudyconfig.cpp -ldemangle > /dev/null 2>&1" );
$string = "";
if( $ret == 0 )
{
    $string = `./tmp_cloudyconfig.out`;
    chomp( $string );
}
if( "$string" eq "ParseCrashDo(Parser&)" )
{
    $res .= "demangle_sun ";
}
unlink "tmp_cloudyconfig.cpp";
unlink "tmp_cloudyconfig.out";
# remove trailing spaces
$res =~ s/ +$//;
if( $res ne "" ) {
    print "$res\n";
}
