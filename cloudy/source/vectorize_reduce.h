/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_REDUCE_H
#define VECTORIZE_REDUCE_H

// NB NB -- do not include this file directly, include vectorize.h

// calculate Sum( a[i] )
double reduce_a(const double* a, long ilo, long ihi);
// calculate Sum( a[i] )
sys_float reduce_a(const sys_float* a, long ilo, long ihi);


// calculate Sum( a[i]*b[i] )
double reduce_ab(const double* a, const double* b, long ilo, long ihi);
// calculate Sum( a[i]*double(b[i]) )
double reduce_ab(const double* a, const sys_float* b, long ilo, long ihi);
// calculate Sum( a[i]*b[i] )
sys_float reduce_ab(const sys_float* a, const sys_float* b, long ilo, long ihi);


// calculate Sum( a[i]*b[i]*c[i] )
double reduce_abc(const double* a, const double* b, const double* c, long ilo, long ihi);
// calculate Sum( a[i]*b[i]*c[i] )
double reduce_abc(const double* a, const double* b, const sys_float* c, long ilo, long ihi);
// calculate Sum( a[i]*b[i]*c[i] )
double reduce_abc(const double* a, const sys_float* b, const sys_float* c, long ilo, long ihi);
// calculate Sum( a[i]*b[i]*c[i] )
sys_float reduce_abc(const sys_float* a, const sys_float* b, const sys_float* c, long ilo, long ihi);


// calculate Sum( a[i]*b[i] ), sum_a = Sum( a[i] )
double reduce_ab_a(const double* a, const double* b, long ilo, long ihi, double* sum_a);
// calculate Sum( double(a[i])*b[i] ), sum_a = Sum( double(a[i]) )
double reduce_ab_a(const sys_float* a, const double* b, long ilo, long ihi, double* sum_a);
// calculate Sum( a[i]*double(b[i]) ), sum_a = Sum( a[i] )
double reduce_ab_a(const double* a, const sys_float* b, long ilo, long ihi, double* sum_a);
// calculate Sum( a[i]*b[i] ), sum_a = Sum( a[i] )
sys_float reduce_ab_a(const sys_float* a, const sys_float* b, long ilo, long ihi, sys_float* sum_a);


// calculate Sum( a[i]*b[i]*c[i] ), sum_ab = Sum( a[i]*b[i] )
double reduce_abc_ab(const double* a, const double* b, const double* c, long ilo, long ihi, double* sum_ab);

#endif
