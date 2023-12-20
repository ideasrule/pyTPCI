/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <fenv.h>
#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include "cddefines.h"
#include "ran.h"

int main ()
{
	ioQQQ = stdout;
	// disable FP exceptions, they would lead to spurious crashes in the tests
	fesetenv(FE_DFL_ENV);
	// initialize random number generator here since it is used in multiple units
	ran.init(0x0efae2bd7f0a3fddULL, cpu.i().nRANK());
	// initialize search path for open_data()
	cpu.i().initPath();
	// we do not want a backtrace since we are deliberately causing problems...
	cpu.i().disableBacktrace();
	return UnitTest::RunAllTests();
}
