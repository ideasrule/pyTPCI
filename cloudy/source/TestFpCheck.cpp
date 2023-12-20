/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"

namespace {
	TEST(FpBoundRealnumOK)
	{
		realnum lo=1.,hi=2.,x=1.5;
		CHECK(fp_bound(lo,x,hi));
	}
	TEST(FpBoundRealnumFail)
	{
		realnum lo=1.,hi=2.,x=2.5;
		CHECK(!fp_bound(lo,x,hi));
	}
	TEST(FpBoundDoubleOK)
	{
		double lo=1.,hi=2.,x=1.5;
		CHECK(fp_bound(lo,x,hi));
	}
	TEST(FpBoundDoubleFail)
	{
		double lo=1.,hi=2.,x=2.5;
		CHECK(!fp_bound(lo,x,hi));
	}
	TEST(FpBoundRealnumTolOK)
	{
		realnum lo=1.,hi=2.,x=1.5;
		CHECK(fp_bound_tol(lo,x,hi,(realnum)1.0e-3));
	}
	TEST(FpBoundRealnumTolFail)
	{
		realnum lo=1.,hi=2.,x=2.5;
		CHECK(!fp_bound_tol(lo,x,hi,(realnum)1.0e-3));
	}
	TEST(FpBoundRealnumTolWide)
	{
		realnum lo=1.,hi=2.,x=2.5;
		CHECK(fp_bound_tol(lo,x,hi,(realnum)1.0));
	}
	TEST(FpBoundDoubleTolOK)
	{
		double lo=1.,hi=2.,x=1.5;
		CHECK(fp_bound_tol(lo,x,hi,1e-3));
	}
	TEST(FpBoundDoubleTolFail)
	{
		double lo=1.,hi=2.,x=2.5;
		CHECK(!fp_bound_tol(lo,x,hi,1e-3));
	}
	TEST(FpBoundDoubleTolWide)
	{
		double lo=1.,hi=2.,x=2.5;
		CHECK(fp_bound_tol(lo,x,hi,1.0));
	}
}
