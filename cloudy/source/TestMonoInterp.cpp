/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "monointerp.h"

#undef NUMVALS
#define NUMVALS(a) sizeof(a)/sizeof((a)[0])

namespace {
	static double xvals1[] = {0.0,1.0}, 
		yvals1[] = {0.0,1.0},
		xvals2[] = {0.0,1.0,2.0,3.0,4.0}, 
		yvals2[] = {0.0,1.0,4.0,9.0,16.0},
		xvals3[] = {0.0,1.0, 2.0,3.0,4.0,5.0,6.0}, 
		yvals3[] = {0.0,1.0,-1.0,2.0,1.0,1.0,0.0};
	struct MonoInterpFixture
	{
		Monointerp m1, m2, m3;

		MonoInterpFixture () : 
			m1(xvals1,yvals1,NUMVALS(xvals1)),
			m2(xvals2,yvals2,NUMVALS(xvals2)),
			m3(xvals3,yvals3,NUMVALS(xvals3))
			{}
		~MonoInterpFixture () {}
	};
	TEST_FIXTURE(MonoInterpFixture,TestMonoInterpLimits)
	{
		CHECK( fp_equal(0.0,m1(0.0)) );
		CHECK( fp_equal(1.0,m1(1.0)) );
	}
	TEST_FIXTURE(MonoInterpFixture,TestMonoInterpInterp)
	{
		CHECK( fp_equal(0.5,m1(0.5)) );
		CHECK( fp_equal(0.333,m1(0.333)) );
		CHECK( fp_equal(0.99,m1(0.99)) );
	}
	TEST_FIXTURE(MonoInterpFixture,TestMonoInterpExtrap)
	{
		CHECK( fp_equal(0.,m1(-0.5)) );
		CHECK( fp_equal(1.,m1(1.5)) );
	}
	TEST_FIXTURE(MonoInterpFixture,TestMonoInterpQuad)
	{
		for (unsigned int i=0; i<NUMVALS(xvals2); ++i)
			CHECK( fp_equal(yvals2[i],m2(xvals2[i])) );
		CHECK_CLOSE(1.5*1.5,m2(1.5),5e-2);
		CHECK_CLOSE(1.1*1.1,m2(1.1),5e-2);
		CHECK_CLOSE(1.9*1.9,m2(1.9),5e-2);
	}
	TEST_FIXTURE(MonoInterpFixture,TestMonoInterpQuadExtrap)
	{
		CHECK( fp_equal(yvals2[0],m2(-0.5)) );
		CHECK( fp_equal(yvals2[NUMVALS(yvals2)-1],m2(NUMVALS(yvals2)-0.5)) );
	}
	TEST_FIXTURE(MonoInterpFixture,TestMonoInterpJump)
	{
		for (unsigned int i=0; i<NUMVALS(xvals3); ++i)
			CHECK( fp_equal(yvals3[i],m3(xvals3[i])) );
	}
	TEST_FIXTURE(MonoInterpFixture,TestMonoInterpJumpExtrap)
	{
		CHECK( fp_equal(yvals3[0],m3(-0.5)) );
		CHECK( fp_equal(yvals3[NUMVALS(yvals3)-1],m3(NUMVALS(yvals3)-0.5)) );
	}
	TEST_FIXTURE(MonoInterpFixture,TestMonoInterpJumpMonotonic)
	{
		const double eps=1e-5;
		for (unsigned int i=0; i<NUMVALS(xvals3)-1; ++i)
		{
			CHECK((yvals3[i]-m3(xvals3[i]+eps))*(yvals3[i+1]-m3(xvals3[i]+eps))<=0);
			CHECK((yvals3[i]-m3(xvals3[i+1]-eps))
			      *(yvals3[i+1]-m3(xvals3[i+1]-eps))<=0);
		}
	}
}
