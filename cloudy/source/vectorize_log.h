/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_LOG_H
#define VECTORIZE_LOG_H

// NB NB -- do not include this file directly, include vectorize.h

// calculate y[i] = log(x[i]), for i=nlo; i < nhi;
void vlog(const double x[], double y[], long nlo, long nhi);

// calculate y[i] = log10(x[i]), for i=nlo; i < nhi;
void vlog10(const double x[], double y[], long nlo, long nhi);

// calculate y[i] = log1p(x[i]), for i=nlo; i < nhi;
void vlog1p(const double x[], double y[], long nlo, long nhi);

// calculate y[i] = logf(x[i]), for i=nlo; i < nhi;
void vlog(const sys_float x[], sys_float y[], long nlo, long nhi);

// calculate y[i] = log10f(x[i]), for i=nlo; i < nhi;
void vlog10(const sys_float x[], sys_float y[], long nlo, long nhi);

// calculate y[i] = log1pf(x[i]), for i=nlo; i < nhi;
void vlog1p(const sys_float x[], sys_float y[], long nlo, long nhi);

void vlog(double *y, double x0, double x1, double x2, double x3);
void vlog10(double *y, double x0, double x1, double x2, double x3);
void vlog1p(double *y, double x0, double x1, double x2, double x3);
void vlog(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
void vlog10(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
void vlog1p(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
void vlog(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3);
void vlog10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3);
void vlog1p(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3);
void vlog(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	  sys_float x6, sys_float x7);
void vlog10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5,
	    sys_float x6, sys_float x7);
void vlog1p(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5,
	    sys_float x6, sys_float x7);
void vlog(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	  sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	  sys_float x13, sys_float x14, sys_float x15);
void vlog10(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	    sys_float x13, sys_float x14, sys_float x15);
void vlog1p(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	    sys_float x13, sys_float x14, sys_float x15);

#endif
