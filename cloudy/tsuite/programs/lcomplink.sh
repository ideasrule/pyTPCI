#!/bin/sh
clang++ $1/$1.cpp -I ../../source  ../../source/sys_llvm/libcloudy.a -o $1/$1.exe
