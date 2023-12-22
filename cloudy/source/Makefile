# Makefile for compiling Cloudy on a UNIX platform
# This works only with GNU make!

# on modern Mac Xcode and FreeBSD releases, g++ is disappearing in favor of clang++
# so make sure we have a fallback in case g++ is missing...
DEFCXX := ${shell which g++ 2> /dev/null}
ifeq ($(DEFCXX),)
  DEFCXX := ${shell which clang++ 2> /dev/null}
endif

STD = -std=c++11
OPT = -O3 -ftrapping-math -fno-math-errno -fasynchronous-unwind-tables -Wno-deprecated-declarations
DEBUGOPT = -finline -fasynchronous-unwind-tables
CXX = $(DEFCXX)
# -W is the old name for -Wextra, g++ 3.3 only supports the former
CXXFLAGS = $(STD) $(OPT) -Wall -W -g
CXXFLAGSEXC = -fnon-call-exceptions
LDFLAGS = $(OPT) -Wall -g
LDLIBS =
NO_TRAP_FLAGS =
# this is intended for supplying extra command line parameters, e.g. extra -D flags
EXTRA = 
SRCDIR = .
OBJEXT = o
LIBEXT = a

.PHONY: all default data debug testclean clean distclean debug-test test valgrind-test depend echo-cxx echo-cxxflags

# Default rule comes first -- having this implement 'make all' is standard
all: default

# include a configuration file if one exists (or carry on without
# worrying if it doesn't)
-include Makefile.conf

main := maincl.cpp vh128sum.cpp
fullsource := ${notdir ${wildcard ${SRCDIR}/*.cpp}}
source := ${filter-out Test%,${fullsource}}
testsource := ${filter Test%,${fullsource}}
libsource := ${filter-out ${main},${source}}
includes := ${notdir ${wildcard ${SRCDIR}/*.h}}
fullobjects := ${notdir ${wildcard *.${OBJEXT}}}
fulldepends := ${fullobjects:.${OBJEXT}=.d}
libobjects := ${libsource:.cpp=.${OBJEXT}}
testobjects := ${testsource:.cpp=.${OBJEXT}}

ifneq (${SRCDIR},.)
  vpath %.cpp ${SRCDIR}
  vpath %.h   ${SRCDIR}
endif

# this is the compiler used to create the dependencies
CXXDEPEND = $(DEFCXX)

# check compiler capabilities
CAPABILITIES := $(shell $(SRCDIR)/capabilities.pl $(CXX))
PRECOMPILE := $(filter precompile,$(CAPABILITIES))
VECTORIZE := $(filter vectorize,$(CAPABILITIES))
NOPIE := $(filter -no-pie,$(CAPABILITIES))
FNOPIE := $(filter -fno-pie,$(CAPABILITIES))
NATIVE := $(filter native,$(CAPABILITIES))
DYNAMIC := $(filter dynamic,$(CAPABILITIES))
DEMANGLE := $(filter demangle%,$(CAPABILITIES))
RPATH := $(filter rpath=%,$(CAPABILITIES))

# DEP_GOALS will be an empty string for targets that don't need dependencies
NODEP_GOALS = clean distclean data
DEP_GOALS = ${filter-out ${NODEP_GOALS},${MAKECMDGOALS}}

CXXFLAGS += $(NOPIE) $(FNOPIE)
LDFLAGS += $(NOPIE) $(FNOPIE)

# the default target does need dependencies
ifeq (${MAKECMDGOALS},)
  DEP_GOALS = need_dependencies
endif

ifeq ($(VECTORIZE),vectorize)
  OPT += -ftree-vectorize
endif

ifeq ($(NATIVE),native)
  OPT += -march=native
endif

ifeq ($(DYNAMIC),dynamic)
  LDFLAGS += -Wl,-export-dynamic
endif

ifeq ($(DEMANGLE),demangle_sun)
  LDLIBS += -ldemangle
endif

ifneq ($(RPATH),)
  LDFLAGS += -Wl,-$(RPATH)
endif

DEBUG_GOALS = $(filter debug%,$(MAKECMDGOALS))
ifneq ($(DEBUG_GOALS),)
  OPT = $(DEBUGOPT)
endif

CXXDEFINES = ${filter -D%,${CXXFLAGS}}

# Precompiling the cddefines.h header speeds up compilation dramatically with g++
# versions from 3.4.0 onwards -- unfortunately, only one header can be precompiled.
# See http://gcc.gnu.org/onlinedocs/gcc-4.2.1/gcc/Precompiled-Headers.html

ifeq (${SRCDIR},.)
  presource := cddefines.h
else
  presource := ${SRCDIR}/cddefines.h
endif
substname := ${presource}
preobjects :=
predepends :=
ifeq ($(PRECOMPILE),precompile)
  substname := cddefines.h.gch
  preobjects := ${substname}
  predepends := ${preobjects:.gch=.d}
endif

RANLIB := ${shell which ranlib}
ifneq ($(RANLIB),)
  RANLIB = ranlib
else
  RANLIB = ar ts
endif

PWD := ${shell pwd}

ifneq ($(word 2,$(PWD)),)
  $(error Error: an installation path containing whitespace is invalid. See https://nublado.org/wiki/DownLoad)
endif

CLDCONFIG = $(PWD)/cloudyconfig.h

CIP = $(PWD)/${SRCDIR}/..
ifneq (${CLOUDY_LAPACK_PATH},)
  CXXFLAGS += -DLAPACK -I${CLOUDY_LAPACK_PATH}/include
  LDLIBS += -L${CLOUDY_LAPACK_PATH}/lib -llapack
endif

ifeq (${MAKECMDGOALS},valgrind-test)
  CXXFLAGS += -DVALGRIND
endif

# add in user-defined extra flags
CXXFLAGS += $(EXTRA)
LDFLAGS += $(EXTRA)

CXXFLAGS += -DSYS_CONFIG=\"$(CLDCONFIG)\"

# this is a collection of compiler flags that enable FPE traps
TRAP_FPE_FILTER = -fnonstd -ftrap=common -Ktrap=fp

# we dont need the -DREVISION here, disable FPE traps for unit tests
CXXUNITTESTFLAGS := $(filter-out ${TRAP_FPE_FILTER},${CXXFLAGS}) $(NO_TRAP_FLAGS)

ifneq ($(wildcard ./gitversion.sh),)
  REVISION := ${shell ./gitversion.sh }
else
  REVISION := ${shell ../gitversion.sh }
endif
CXXFLAGS += -DREVISION=\"$(REVISION)\"

# this is a filter to erase optimization flags from CXXFLAGS
OPTFILTER = -O -O1 -O2 -O3 -O4 -O5 -Os -xO -xO1 -xO2 -xO3 -xO4 -xO5 /O /O1 /O2 /O3 -fnonstd -fsimple=1 -fast -xvector%
CFGFILTER = -DSYS_CONFIG%

CXXFLAGSNOOPT := ${filter-out ${OPTFILTER},${CXXFLAGS}}
CXXFLAGSNOCFG := ${filter-out ${CFGFILTER},${CXXFLAGS}}
CXXFLAGSNOOPTNOCFG := ${filter-out ${CFGFILTER},${CXXFLAGSNOOPT}}

CXXFLAGSVECT = $(CXXFLAGS)
# this will pick up both g++ and clang++
ifeq ($(findstring g++,$(CXX)),g++)
  CXXFLAGSVECT += -ffast-math
endif

# Makefile.targets can override or add to the DEFAULT variable.
DEFAULT = cloudy.exe data vh128sum.exe

# include additional targets from file if one exists (or carry on without
# worrying if it doesn't)
-include Makefile.targets

debug: default

default: $(DEFAULT)

cloudy.exe: ${preobjects} maincl.${OBJEXT} libcloudy.${LIBEXT}
	${CXX} ${LDFLAGS} -o cloudy.exe maincl.${OBJEXT} -L. -lcloudy ${LDLIBS}

libcloudy.${LIBEXT}: ${libobjects}
	ar cr libcloudy.${LIBEXT} $^
	${RANLIB} libcloudy.${LIBEXT}

data:
	cd ${SRCDIR}/../data/cdms+jpl && $(MAKE)

vh128sum.exe: vh128sum.${OBJEXT} vectorhash.${OBJEXT}
	${CXX} ${LDFLAGS} -o vh128sum.exe vh128sum.${OBJEXT} vectorhash.${OBJEXT}

testclean:
	rm -f Test*.${OBJEXT}
	rm -f Test*.d
	rm -f runtests.exe
	rm -rf runtests.exe.dSYM
	rm -f temp.??????

clean: testclean
	rm -f *.${OBJEXT}
	rm -f *.d
	rm -f *.h.gch
	rm -rf SunWS_cache
	rm -rf cloudy.exe.dSYM
	rm -rf vh128sum.exe.dSYM
	rm -f libcloudy.${LIBEXT}
	rm -f $(DEFAULT)
	rm -rf tmp_cloudyconfig.*
	rm -f cloudyconfig.h

distclean: clean
	rm -f ${SRCDIR}/Makefile.conf
	rm -rf lib
	rm -rf include
	cd ${SRCDIR}/../library/UnitTest++ && $(MAKE) clean
	cd ${SRCDIR}/../data/cdms+jpl && $(MAKE) clean

# Ensure version update if *any* source files change
version.${OBJEXT}: ${main} ${libsource} ${includes}

DEPCMD := ${CXXDEPEND} ${STD} -MM -MG -DMM ${CXXDEFINES}

cpu.${OBJEXT}: ${SRCDIR}/cpu.cpp ${CLDCONFIG} ${preobjects}
	${CXX} ${CXXFLAGS} ${CXXFLAGSEXC} -DCLOUDY_ROOT=\"${CIP}\" -c -o $@ $<
	${SRCDIR}/make_depend.pl "${DEPCMD} $<" ${presource} ${substname} ${CLDCONFIG} ${@:.${OBJEXT}=.d}

# an explicit -O0 is needed for some compilers like icc
parse_crashdo.${OBJEXT}: ${SRCDIR}/parse_crashdo.cpp ${CLDCONFIG} ${preobjects}
	${CXX} ${CXXFLAGSNOOPT} ${DEBUGOPT} -c -o $@ $<
	${SRCDIR}/make_depend.pl "${DEPCMD} $<" ${presource} ${substname} ${CLDCONFIG} ${@:.${OBJEXT}=.d}

# vectorize_reduce.cpp needs more agressive optimization to enable vectorization of the loops in this file
vectorize_reduce.${OBJEXT}: ${SRCDIR}/vectorize_reduce.cpp ${CLDCONFIG} ${preobjects}
	${CXX} ${CXXFLAGSVECT} -c -o $@ $<
	${SRCDIR}/make_depend.pl "${DEPCMD} $<" ${presource} ${substname} ${CLDCONFIG} ${@:.${OBJEXT}=.d}

Test%.${OBJEXT}: Test%.cpp lib/UnitTest++/README ${CLDCONFIG} ${preobjects}
	${CXX} ${CXXUNITTESTFLAGS} -Ilib/UnitTest++/src -I${SRCDIR} -c -o $@ $<
	${SRCDIR}/make_depend.pl "${DEPCMD} $<" ${presource} ${substname} ${CLDCONFIG} ${@:.${OBJEXT}=.d}

%.${OBJEXT}: %.cpp ${CLDCONFIG} ${preobjects}
	${CXX} ${CXXFLAGS} -c -o $@ $<
	${SRCDIR}/make_depend.pl "${DEPCMD} $<" ${presource} ${substname} ${CLDCONFIG} ${@:.${OBJEXT}=.d}

%.h.gch: %.h ${CLDCONFIG}
	${CXX} ${CXXFLAGS} -o $@ $<
	${SRCDIR}/make_depend.pl "${DEPCMD} $<" ${@:.h.gch=.o} $@ ${CLDCONFIG} ${@:.gch=.d}

${CLDCONFIG}: ${SRCDIR}/configure.sh
	@echo "Creating cloudyconfig.h"
	@${SRCDIR}/configure.sh ${CXX} ${CXXFLAGSNOCFG} -I${SRCDIR}

debug-test: test

test: runtests.exe
	./runtests.exe

valgrind-test: runtests.exe
	valgrind --leak-check=full --track-fds=yes ./runtests.exe

TESTLIBS := -L. -lcloudy -Llib/UnitTest++ -lUnitTest++

runtests.exe: lib/UnitTest++/libUnitTest++.${LIBEXT} libcloudy.${LIBEXT} ${testobjects}
	${CXX} ${CXXFLAGS} -${OBJEXT} $@ ${filter-out %.${LIBEXT}, $^} ${LDFLAGS} ${TESTLIBS} ${LDLIBS}

lib/UnitTest++/libUnitTest++.${LIBEXT}: lib/UnitTest++/README
	cd lib/UnitTest++; \
	$(MAKE) CC=${CXX} CCFLAGS="-O0 ${CXXFLAGSNOOPTNOCFG}" LDFLAGS="${LDFLAGS}" LDLIBS="${LDLIBS}"

# each sys_xxx dir has its own copy of UnitTest++ as it is compiler dependent.
# this also assures that if we do "make distclean" in another sys_xxx, it does
# not break the installation here...
# make the README file the target as using the directory lib/UnitTest++ is not a good idea
# directories get their timestamp changed whenever a file inside is changed...
lib/UnitTest++/README:
	mkdir -p lib
	cp -a ${SRCDIR}/../library/UnitTest++ lib
	cd lib/UnitTest++ && $(MAKE) clean

depend:
	@echo "make depend is no longer needed."

echo-cxx:
	@echo "$(CXX)"

echo-cxxflags:
	@echo "$(CXXFLAGS)"

# only include dependencies when we really need them, this prevents unnecessary builds of the dependencies
# the way we create dependency files is discussed in http://scottmcpeak.com/autodepend/autodepend.html
ifneq (${DEP_GOALS},)
-include ${predepends} ${fulldepends}
endif
