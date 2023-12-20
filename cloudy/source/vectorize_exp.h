/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_EXP_H
#define VECTORIZE_EXP_H

// NB NB -- do not include this file directly, include vectorize.h

// calculate y[i] = exp(x[i]), for i=nlo; i < nhi;
void vexp(const double x[], double y[], long nlo, long nhi);

// calculate y[i] = exp10(x[i]), for i=nlo; i < nhi;
void vexp10(const double x[], double y[], long nlo, long nhi);

// calculate y[i] = expm1(x[i]), for i=nlo; i < nhi;
void vexpm1(const double x[], double y[], long nlo, long nhi);

// calculate y[i] = expf(x[i]), for i=nlo; i < nhi;
void vexp(const sys_float x[], sys_float y[], long nlo, long nhi);

// calculate y[i] = exp10f(x[i]), for i=nlo; i < nhi;
void vexp10(const sys_float x[], sys_float y[], long nlo, long nhi);

// calculate y[i] = expm1f(x[i]), for i=nlo; i < nhi;
void vexpm1(const sys_float x[], sys_float y[], long nlo, long nhi);

void vexp(double *y, double x0, double x1, double x2, double x3);
void vexp10(double *y, double x0, double x1, double x2, double x3);
void vexpm1(double *y, double x0, double x1, double x2, double x3);
void vexp(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
void vexp10(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
void vexpm1(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
void vexp(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3);
void vexp10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3);
void vexpm1(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3);
void vexp(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	  sys_float x6, sys_float x7);
void vexp10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5,
	    sys_float x6, sys_float x7);
void vexpm1(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5,
	    sys_float x6, sys_float x7);
void vexp(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	  sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	  sys_float x13, sys_float x14, sys_float x15);
void vexp10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	    sys_float x13, sys_float x14, sys_float x15);
void vexpm1(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	    sys_float x13, sys_float x14, sys_float x15);

#endif
