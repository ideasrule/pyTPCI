#!/bin/sh
g++ template.cpp -DSYS_CONFIG=\"`pwd`/../../source/cloudyconfig.h\" -I ../../source -L ../../source -lcloudy -o template.exe
# clang++ template.cpp -DSYS_CONFIG=\"`pwd`/../../source/sys_llvm/cloudyconfig.h\" -I ../../source -L ../../source/sys_llvm -lcloudy -o template.exe
