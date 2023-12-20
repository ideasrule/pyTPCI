/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "vectorize.h"

//
// This file should ONLY contain simple routines that implement vectorizable loops.
// These loops could then also be targets for parallelization, e.g. with openMP.
// It is imperative to check that the loops are written up in such a way that the
// compiler allows the vectorization to take place. With g++ this can be done with
//   g++ [ other options ] -ftree-vectorizer-verbose=1 <somefile.cpp>
// This will generate a report about loop vectorization on stderr.
//
// Note that for the reduction loops more aggressive optimization options are needed
// in order for them to be vectorized. Most importantly, the compiler assumes that
// -fno-trapping-math is in effect in this file. This can cause problems with
// conditional statements, e.g. like this one
//   y = ( x > 0. ) ? log(x) : -70.;
// These sort of conditionals should be avoided here.
//

// calculate Sum( a[i] )
double reduce_a(const double* a, long ilo, long ihi)
{
	DEBUG_ENTRY( "reduce_a()" );

	double sum = 0.;
	for( long i=ilo; i < ihi; ++i )
		sum += a[i];
	return sum;
}

// calculate Sum( a[i] )
sys_float reduce_a(const sys_float* a, long ilo, long ihi)
{
	DEBUG_ENTRY( "reduce_a()" );

	sys_float sum = 0.f;
	for( long i=ilo; i < ihi; ++i )
		sum += a[i];
	return sum;
}

// calculate Sum( a[i]*b[i] )
double reduce_ab(const double* a, const double* b, long ilo, long ihi)
{
	DEBUG_ENTRY( "reduce_ab()" );

	double sum = 0.;
	for( long i=ilo; i < ihi; ++i )
		sum += a[i]*b[i];
	return sum;
}

// calculate Sum( a[i]*double(b[i]) )
double reduce_ab(const double* a, const sys_float* b, long ilo, long ihi)
{
	DEBUG_ENTRY( "reduce_ab()" );

	double sum = 0.;
	for( long i=ilo; i < ihi; ++i )
		sum += a[i]*double(b[i]);
	return sum;
}

// calculate Sum( a[i]*b[i] )
sys_float reduce_ab(const sys_float* a, const sys_float* b, long ilo, long ihi)
{
	DEBUG_ENTRY( "reduce_ab()" );

	sys_float sum = 0.f;
	for( long i=ilo; i < ihi; ++i )
		sum += a[i]*b[i];
	return sum;
}

// calculate Sum( a[i]*b[i]*c[i] )
double reduce_abc(const double* a, const double* b, const double* c, long ilo, long ihi)
{
	DEBUG_ENTRY( "reduce_abc()" );

	double sum = 0.;
	for( long i=ilo; i < ihi; ++i )
		sum += a[i]*b[i]*c[i];
	return sum;
}

// calculate Sum( a[i]*b[i]*double(c[i]) )
double reduce_abc(const double* a, const double* b, const sys_float* c, long ilo, long ihi)
{
	DEBUG_ENTRY( "reduce_abc()" );

	double sum = 0.;
	for( long i=ilo; i < ihi; ++i )
		sum += a[i]*b[i]*double(c[i]);
	return sum;
}

// calculate Sum( a[i]*double(b[i])*double(c[i]) )
double reduce_abc(const double* a, const sys_float* b, const sys_float* c, long ilo, long ihi)
{
	DEBUG_ENTRY( "reduce_abc()" );

	double sum = 0.;
	for( long i=ilo; i < ihi; ++i )
		sum += a[i]*double(b[i])*double(c[i]);
	return sum;
}

// calculate Sum( a[i]*b[i]*c[i] )
sys_float reduce_abc(const sys_float* a, const sys_float* b, const sys_float* c, long ilo, long ihi)
{
	DEBUG_ENTRY( "reduce_abc()" );

	sys_float sum = 0.f;
	for( long i=ilo; i < ihi; ++i )
		sum += a[i]*b[i]*c[i];
	return sum;
}

// calculate Sum( a[i]*b[i] ), sum_a = Sum( a[i] )
double reduce_ab_a(const double* a, const double* b, long ilo, long ihi, double* sum_a)
{
	DEBUG_ENTRY( "reduce_ab_a()" );

	double sum1 = 0.;
	double sum2 = 0.;
	for( long i=ilo; i < ihi; ++i )
	{
		double one = a[i];
		sum1 += one;
		sum2 += one*b[i];
	}
	*sum_a = sum1;
	return sum2;
}

// calculate Sum( double(a[i])*b[i] ), sum_a = Sum( double(a[i]) )
double reduce_ab_a(const sys_float* a, const double* b, long ilo, long ihi, double* sum_a)
{
	DEBUG_ENTRY( "reduce_ab_a()" );

	double sum1 = 0.;
	double sum2 = 0.;
	for( long i=ilo; i < ihi; ++i )
	{
		double one = double(a[i]);
		sum1 += one;
		sum2 += one*b[i];
	}
	*sum_a = sum1;
	return sum2;
}

// calculate Sum( a[i]*double(b[i]) ), sum_a = Sum( a[i] )
double reduce_ab_a(const double* a, const sys_float* b, long ilo, long ihi, double* sum_a)
{
	DEBUG_ENTRY( "reduce_ab_a()" );

	double sum1 = 0.;
	double sum2 = 0.;
	for( long i=ilo; i < ihi; ++i )
	{
		double one = a[i];
		sum1 += one;
		sum2 += one*double(b[i]);
	}
	*sum_a = sum1;
	return sum2;
}

// calculate Sum( a[i]*b[i] ), sum_a = Sum( a[i] )
sys_float reduce_ab_a(const sys_float* a, const sys_float* b, long ilo, long ihi, sys_float* sum_a)
{
	DEBUG_ENTRY( "reduce_ab_a()" );

	sys_float sum1 = 0.f;
	sys_float sum2 = 0.f;
	for( long i=ilo; i < ihi; ++i )
	{
		sys_float one = a[i];
		sum1 += one;
		sum2 += one*b[i];
	}
	*sum_a = sum1;
	return sum2;
}

// calculate Sum( a[i]*b[i]*c[i] ), sum_ab = Sum( a[i]*b[i] )
double reduce_abc_ab(const double* a, const double* b, const double* c, long ilo, long ihi, double* sum_ab)
{
	DEBUG_ENTRY( "reduce_abc_ab()" );

	double sum2 = 0.;
	double sum3 = 0.;
	for( long i=ilo; i < ihi; ++i )
	{
		double one = a[i]*b[i];
		sum2 += one;
		sum3 += one*c[i];
	}
	*sum_ab = sum2;
	return sum3;
}
