#!/bin/sh

# This is run as ./complink.sh arg1 [arg2]
# with no second argument will use gcc and source/
# arg2='gcc' or arg2='clang'
# if arg2 is given then it will compile in source/sys_gcc/ or source/sys_llvm
# directories respectively

dir=""
comp="g++"

if [ "$#" -eq "2" ] && [ "$2" != "gcc" ] && [ "$2" != "clang" ]; then
    echo "Sorry, I don't understand this compiler name. Your option was:" $2
    echo "Accepted options are 'clang' and 'gcc'"
    exit 1
fi

if [ "$2" = "gcc" ]; then
    dir="sys_gcc"
    comp="g++"
elif [ "$2" = "clang" ]; then
    dir="sys_llvm"
    comp="clang++"
fi

$comp $1/$1.cpp -DSYS_CONFIG=\"`pwd`/../../source/$dir/cloudyconfig.h\" -I ../../source  -I ../../source/$dir ../../source/$dir/libcloudy.a -o $1/$1.exe
