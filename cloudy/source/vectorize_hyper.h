/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORIZE_HYPERBOLIC_H
#define VECTORIZE_HYPERBOLIC_H

// calculate y[i] = asinh(x[i]), for i=nlo; i < nhi;
void vasinh(const double x[], double y[], long nlo, long nhi);

// fast version of vasinh with restricted domain [0,sqrt(DBL_MAX)]
void vfast_asinh(const double x[], double y[], long nlo, long nhi);

// calculate y[i] = asinhf(x[i]), for i=nlo; i < nhi;
void vasinh(const sys_float x[], sys_float y[], long nlo, long nhi);

// fast version of vasinh with restricted domain [0,sqrt(FLT_MAX)]
void vfast_asinh(const sys_float x[], sys_float y[], long nlo, long nhi);

void vasinh(double *y, double x0, double x1, double x2, double x3);
void vfast_asinh(double *y, double x0, double x1, double x2, double x3);
void vasinh(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
void vfast_asinh(double *y, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
void vasinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3);
void vfast_asinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3);
void vasinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7);
void vfast_asinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
		 sys_float x6, sys_float x7);
void vasinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
	    sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
	    sys_float x13, sys_float x14, sys_float x15);
void vfast_asinh(sys_float *y, sys_float x0, sys_float x1, sys_float x2, sys_float x3, sys_float x4, sys_float x5, 
		 sys_float x6, sys_float x7, sys_float x8, sys_float x9, sys_float x10, sys_float x11, sys_float x12,
		 sys_float x13, sys_float x14, sys_float x15);

#endif
