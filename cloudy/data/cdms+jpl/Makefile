# Makefile for converting the CDMS / JPL data files

# on modern Mac Xcode and FreeBSD releases, g++ is disappearing in favor of clang++
# so make sure we have a fallback in case g++ is missing...
DEFCXX := ${shell which g++ 2> /dev/null}
ifeq ($(DEFCXX),)
  DEFCXX := ${shell which clang++ 2> /dev/null}
endif

# use 'override' here to make sure we do not inherit the settings from the main compilation
# these could be for an incompatible compiler...
override CXX = $(DEFCXX)
override OPT = -O3 -ftrapping-math -fno-math-errno
# -W is the old name for -Wextra, g++ 3.3 only supports the former
override CXXFLAGS = -ansi ${OPT} -Wall -W -g
override EXTRA =

.PHONY: all clean squeaky-clean

data := ${notdir ${wildcard *.cat}}
processed := ${data:.cat=.dat}

all: ${processed}

convert_calpgm: convert_calpgm.cpp
	${CXX} ${CXXFLAGS} ${EXTRA} -o convert_calpgm convert_calpgm.cpp

%.dat: %.cat convert_calpgm
	./convert_calpgm $<

clean:
	rm -f convert_calpgm
	rm -rf convert_calpgm.dSYM
	rm -f c??????.dat

squeaky-clean: clean
	rm -f partition_function_cdms.html
	rm -f partition_function_jpl.html
	rm -f c??????.cat
	rm -f masterlist
