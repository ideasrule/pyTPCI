/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_SQRT_H
#define VECTORIZE_SQRT_H

// NB NB -- do not include this file directly, include vectorize.h

// calculate y[i] = sqrt(x[i]), for i=nlo; i < nhi;
void vsqrt(const double x[], double y[], long nlo, long nhi);

// calculate y[i] = hypot(x1[i], x2[i]), for i=nlo; i < nhi;
void vhypot(const double x1[], const double x2[], double y[], long nlo, long nhi);

// calculate y[i] = sqrtf(x[i]), for i=nlo; i < nhi;
void vsqrt(const sys_float x[], sys_float y[], long nlo, long nhi);

// calculate y[i] = hypotf(x1[i], x2[i]), for i=nlo; i < nhi;
void vhypot(const sys_float x1[], const sys_float x2[], sys_float y[], long nlo, long nhi);

void vsqrt(double *y, double x0, double x1, double x2, double x3);
void vhypot(double *y, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);
void vsqrt(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
void vsqrt(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3);
void vhypot(sys_float *y, sys_float x0, sys_float y0, sys_float x1, sys_float y1, sys_float x2, sys_float y2, 
	    sys_float x3, sys_float y3);
void vsqrt(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	   sys_float x6, sys_float x7);
void vhypot(sys_float *y, sys_float x0, sys_float y0, sys_float x1, sys_float y1, sys_float x2, sys_float y2, 
	    sys_float x3, sys_float y3, sys_float x4, sys_float y4, sys_float x5, sys_float y5, sys_float x6,
	    sys_float y6, sys_float x7, sys_float y7);
void vsqrt(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	   sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	   sys_float x13, sys_float x14, sys_float x15);

#endif
