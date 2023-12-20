/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "thirdparty_quadpack.h"
#include "physconst.h"

namespace {

	sys_float qngf(sys_float x)
	{
		sys_float ret_val;
		
		ret_val = exp(x) / (x * x + 1.f);
		return ret_val;
	}
	
	TEST(TestQng)
	{
		sys_float a, b;
		long ier;
		long neval;
		sys_float epsabs, abserr, epsrel, result;
		
		a = 0.f;
		b = 1.f;
		epsabs = 0.f;
		epsrel = .001f;
		qng_(E_fp_fp(qngf), 
			  &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &
			  ier);
		if (0) printf("QNG: Result=%.8g abserr=%.8g neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 1.2707243f, 1e-5f ) );
		CHECK( fp_equal_tol( abserr, 7.574107e-6f, 1e-9f ) );
		CHECK_EQUAL( 21, neval );
		CHECK_EQUAL( 0, ier );
	}

	sys_float qagf(sys_float x)
	{
		sys_float ret_val;

		ret_val = 2.f / (sin(x * 10.f * sys_float(PI)) + 2.f);
		return ret_val;
	}

	TEST(TestQag)
	{
		sys_float a, b;
		long ier, key;
		long last, lenw;
		sys_float work[400];
		long neval, limit, iwork[100];
		sys_float epsabs, abserr, epsrel, result;

		a = 0.f;
		b = 1.f;
		epsabs = 0.f;
		epsrel = .001f;
		key = 6;
		limit = 100;
		lenw = limit << 2;
		qag_(E_fp_fp(qagf), &a, &b, &epsabs, &epsrel, &key, &result, &abserr, &
			  neval, &ier, &limit, &lenw, &last, iwork, work);
		if (0) printf("QAG: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 1.154701f, 1e-5f ) );
		CHECK( abserr < 0.0008f );
		CHECK_EQUAL( 183, neval );
		CHECK_EQUAL( 0, ier );
	} 

	sys_float qagsf(sys_float x)
	{
		sys_float ret_val;

		ret_val = 0.f;
		if (x > 0.f) {
			ret_val = 1.f / sqrt(x);
		}
		return ret_val;
	} /* qagsf_ */

	TEST(TestQags)
	{
		sys_float a, b;
		long ier;
		long last, lenw;
		sys_float work[400];
		long neval, limit, iwork[100];
		sys_float epsabs, abserr, epsrel, result;

		a = 0.f;
		b = 1.f;
		epsabs = 0.f;
		epsrel = .001f;
		limit = 100;
		lenw = limit << 2;
		qags_(E_fp_fp(&qagsf), 
				&a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &
				limit, &lenw, &last, iwork, work);
		if (0) printf("QAGS: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 2.f, 1e-5f ) );
		CHECK( fp_equal_tol( abserr, 5.841255e-6f, 1e-8f ) );
		CHECK_EQUAL( 231, neval );
		CHECK_EQUAL( 0, ier );
	}

	sys_float qagpf(sys_float x)
	{
		sys_float ret_val, r__1, r__2;
		double d__1, d__2;
		const double c_b13 = -.25;
		const double c_b14 = -.55;

		ret_val = 0.f;
		if (x != 1.f/7.f && x != 2.f/3.f) {
			d__1 = (r__1 = x - 1.f/7.f, fabs(r__1));
			d__2 = (r__2 = x - 2.f/3.f, fabs(r__2));
			ret_val = pow(d__1, c_b13) * pow(d__2, c_b14);
		}
		return ret_val;
	}

	TEST(TestQagp)
	{
		sys_float a, b;
		long ier;
		long last, lenw;
		sys_float work[404];
		long npts2;
		long neval, leniw, limit, iwork[204];
		sys_float epsabs, abserr, epsrel, points[4], result;

		a = 0.f;
		b = 1.f;
		npts2 = 4;
		points[0] = 1.f/7.f;
		points[1] = 2.f/3.f;
		epsabs = 0.0;
		epsrel = 1e-3;
		limit = 100;
		leniw = (limit << 1) + npts2;
		lenw = (limit << 2) + npts2;
		qagp_(E_fp_fp(qagpf), 
				&a, &b, &npts2, points, &epsabs, &epsrel, &result, &
				abserr, &neval, &ier, &leniw, &lenw, &last, iwork, work);
		if (0) printf("QAGP: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 4.253580f, 1e-5f ) );
		CHECK( abserr <  8.5e-4f );
		CHECK_EQUAL( neval , 861 );
		CHECK_EQUAL( ier, 0 );
	}

	sys_float qagif(sys_float x)
	{
		sys_float ret_val;

		ret_val = 0.f;
		if (x > 0.f) {
			ret_val = sqrt(x) * log(x) / ((x + 1.f) * (x + 2.f));
		}
		return ret_val;
	} /* qagif_ */

	TEST(TestQagi)
	{
		long inf, ier;
		sys_float boun;
		long last, lenw;
		sys_float work[400];
		long neval, limit, iwork[100];
		sys_float epsabs, abserr, epsrel, result;

		boun = 0.f;
		inf = 1;
		epsabs = 0.f;
		epsrel = .001f;
		limit = 100;
		lenw = limit << 2;
		qagi_(E_fp_fp(qagif), 
				&boun, &inf, &epsabs, &epsrel, &result, &abserr, &
				neval, &ier, &limit, &lenw, &last, iwork, work);
		if (0) printf("QAGI: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 3.079564f, 1e-4f ) );
		CHECK( abserr < 9.e-4f );
		CHECK_EQUAL( 315, neval );
		CHECK_EQUAL( 0, ier );
	} /* tqagi_ */

	sys_float qawof(sys_float x)
	{
		sys_float ret_val;

		ret_val = 0.f;
		if (x > 0.f) {
			ret_val = exp(-(x)) * log(x);
		}
		return ret_val;
	} /* qawof_ */

	TEST(TestQawo)
	{
		sys_float a, b;
		long ier, last, lenw;
		sys_float work[925];
		long maxp1;
		sys_float omega;
		long neval;
		long leniw, limit, iwork[200];
		sys_float epsabs, abserr;
		long integr;
		sys_float epsrel, result;

		a = 0.f;
		b = 1.f;
		omega = 10.f;
		integr = 1;
		epsabs = 0.f;
		epsrel = .001f;
		limit = 100;
		leniw = limit << 1;
		maxp1 = 21;
		lenw = (limit << 2) + maxp1 * 25;
		qawo_(E_fp_fp(qawof), 
				&a, &b, &omega, &integr, &epsabs, &epsrel, &result, &
				abserr, &neval, &ier, &leniw, &maxp1, &lenw, &last, iwork, work);
		if (0) printf("QAWO: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, -1.776392e-1f, 1e-5f ) );
		CHECK( abserr < 2.e-7f );
		CHECK_EQUAL( neval , 255 );
		CHECK_EQUAL( ier, 0 );
	} /* tqawo_ */


	sys_float qawff(sys_float x)
	{
		sys_float ret_val=0.f;

		if (x > 0.f) {
			ret_val = sin(x * 50.f) / (x * sqrt(x));
		}
		return ret_val;
	} /* qawff_ */

	TEST(TestQawf)
	{
		sys_float a;
		long ier, lst;
		long lenw;
		sys_float work[1025];
		long maxp1;
		sys_float omega;
		long neval, leniw, limit, iwork[250];
		sys_float epsabs, abserr;
		long integr, limlst;
		sys_float result;

		a = 0.f;
		omega = 8.f;
		integr = 2;
		epsabs = .001f;
		limlst = 50;
		limit = 100;
		leniw = (limit << 1) + limlst;
		maxp1 = 21;
		lenw = (leniw << 1) + maxp1 * 25;
		qawf_(E_fp_fp(qawff), &a, &omega, &integr, &epsabs, &result, &abserr, &
				neval, &ier, &limlst, &lst, &leniw, &maxp1, &lenw, iwork, work);
		if (0) printf("QAWF: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 1.422553f, 1e-5f ) );
		CHECK( abserr < 8e-4f );
		CHECK_EQUAL( neval , 7110 );
		CHECK_EQUAL( ier, 0 );
	}

	sys_float qawsf(sys_float x)
	{
		sys_float ret_val;

		ret_val = sin(x * 10.f);
		return ret_val;
	} /* qawsf_ */

	TEST(TestQaws)
	{
		sys_float a, b;
		long ier;
		sys_float alfa, beta;
		long last, lenw;
		sys_float work[400];
		long neval, limit;
		long iwork[100];
		sys_float epsabs, abserr;
		long integr;
		sys_float epsrel, result;

		a = 0.f;
		b = 1.f;
		alfa = -.5f;
		beta = -.5f;
		integr = 1;
		epsabs = 0.f;
		epsrel = .001f;
		limit = 100;
		lenw = limit << 2;
		qaws_(E_fp_fp(qawsf), 
				&a, &b, &alfa, &beta, &integr, &epsabs, &epsrel, &
				result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
		if (0) printf("QAWS: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 5.350192e-1f, 1e-5f ) );
		CHECK( abserr < 3e-7f );
		CHECK_EQUAL( neval , 50 );
		CHECK_EQUAL( ier, 0 );
	} /* tqaws_ */

	sys_float qawcf(sys_float x)
	{
		sys_float ret_val;

		ret_val = 1.f / (x * x + 1e-4f);
		return ret_val;
	} /* qawcf_ */

	TEST(TestQawc)
	{
		sys_float a, b, c__;
		long ier;
		long last, lenw;
		sys_float work[400];
		long neval, limit, iwork[100];
		sys_float epsabs, abserr, epsrel, result;

		a = -1.f;
		b = 1.f;
		c__ = .5f;
		epsabs = 0.f;
		epsrel = .001f;
		limit = 100;
		lenw = limit << 2;
		qawc_(E_fp_fp(qawcf), 
				&a, &b, &c__, &epsabs, &epsrel, &result, &abserr, &
				neval, &ier, &limit, &lenw, &last, iwork, work);
		if (0) printf("QAWC: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, -6.284617e2f, 1e-3f ) );
		CHECK( abserr < 1.7e-1f );
		CHECK_EQUAL( neval , 225 );
		CHECK_EQUAL( ier, 0 );
	} /* tqawc_ */

	double dqngf(double x)
	{
		double ret_val;

		ret_val = exp(x) / (x * x + 1.);
		return ret_val;
	} /* dqngf_ */

	TEST(TestDQng)
	{
		double a, b;
		long ier;
		long neval;
		double epsabs, abserr, epsrel, result;

		a = 0.;
		b = 1.;
		epsabs = 0.;
		epsrel = .001;
		dqng_(D_fp_fp(dqngf), 
				&a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier);
		if (0) printf("DQNG: Result=%.8g abserr=%.8g neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 1.2707241, 1e-5) );
		CHECK( fp_equal_tol( abserr, 1.4107872e-14, 1e-15) );
		CHECK_EQUAL( neval , 21 );
		CHECK_EQUAL( ier, 0 );
	} /* tdqngf_ */

	double dqagf(double x)
	{
		double ret_val;

		ret_val = 2. / (sin(x * 10.0 * PI) + 2.0);
		return ret_val;
	} /* dqagf_ */

	TEST(TestDQag)
	{
		double a, b;
		long ier, key;
		long last, lenw;
		double work[400];
		long neval, limit, iwork[100];
		double epsabs, abserr, epsrel, result;

		a = 0.;
		b = 1.;
		epsabs = 0.;
		epsrel = .001;
		key = 6;
		limit = 100;
		lenw = limit << 2;
		dqag_(D_fp_fp(dqagf), 
				&a, &b, &epsabs, &epsrel, &key, &result, &abserr, &
				neval, &ier, &limit, &lenw, &last, iwork, work);
		if (0) printf("DQAG: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 1.154701, 1e-5) );
		CHECK( fp_equal_tol( abserr, 7.609303e-4, 1e-7) );
		CHECK_EQUAL( neval , 183 );
		CHECK_EQUAL( ier, 0 );
	} /* tdqagf_ */

	double dqagsf(double x)
	{
		double ret_val;

		ret_val = 0.;
		if (x > 0.) {
			ret_val = 1. / sqrt(x);
		}
		return ret_val;
	} /* dqagsf_ */

	TEST(TestDQags)
	{
		double a, b;
		long ier;
		long last, lenw;
		double work[400];
		long neval, limit, iwork[100];
		double epsabs, abserr, epsrel, result;

		a = 0.;
		b = 1.;
		epsabs = 0.;
		epsrel = .001;
		limit = 100;
		lenw = limit << 2;
		dqags_(D_fp_fp(&dqagsf), 
				 &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &
				 limit, &lenw, &last, iwork, work);
		if (0) printf("DQAGS: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 2., 1e-5) );
		CHECK( abserr < 6.e-15 );
		CHECK_EQUAL( 231, neval );
		CHECK_EQUAL( 0, ier );
	} /* tdqags_ */

	double dqagpf(double x)
	{
		double ret_val, r__1, r__2;
		double d__1, d__2;
		const double c_b13 = -.25;
		const double c_b14 = -.55;

		ret_val = 0.f;
		if (x != 1./7. && x != 2./3.) {
			d__1 = (r__1 = x - 1./7., fabs(r__1));
			d__2 = (r__2 = x - 2./3., fabs(r__2));
			ret_val = pow(d__1, c_b13) * pow(d__2, c_b14);
		}
		return ret_val;
	} /* dqagpf_ */

	TEST(TestDQagp)
	{
		double a, b;
		long ier;
		long last, lenw;
		double work[404];
		long npts2;
		long neval, leniw, limit, iwork[204];
		double epsabs, abserr, epsrel, points[4], result;

		a = 0.;
		b = 1.;
		npts2 = 4;
		points[0] = 1./7.;
		points[1] = 2./3.;
		epsabs = 0.0;
		epsrel = 1e-3;
		limit = 100;
		leniw = (limit << 1) + npts2;
		lenw = (limit << 2) + npts2;
		dqagp_(D_fp_fp(dqagpf), 
				 &a, &b, &npts2, points, &epsabs, &epsrel, &result, &
				 abserr, &neval, &ier, &leniw, &lenw, &last, iwork, work);
		if (0) printf("DQAGP: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 4.253688, 1e-5) );
		CHECK( fp_equal_tol( abserr, 8.758089e-4, 1e-8) );
		CHECK_EQUAL( 861, neval );
		CHECK_EQUAL( 0, ier );
	} /* tdqagp_ */

	double dqagif(double x)
	{
		double ret_val;

		ret_val = 0.;
		if (x > 0.) {
			ret_val = sqrt(x) * log(x) / ((x + 1.) * (x + 2.));
		}
		return ret_val;
	} /* dqagif_ */

	TEST(TestDQagi)
	{
		long inf, ier;
		double boun;
		long last, lenw;
		double work[400];
		long neval, limit, iwork[100];
		double epsabs, abserr, epsrel, result;

		boun = 0.;
		inf = 1;
		epsabs = 0.;
		epsrel = .001;
		limit = 100;
		lenw = limit << 2;
		dqagi_(D_fp_fp(dqagif), 
				 &boun, &inf, &epsabs, &epsrel, &result, &abserr, &
				 neval, &ier, &limit, &lenw, &last, iwork, work);
		if (0) printf("DQAGI: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 3.079555, 1e-5) );
		CHECK( fp_equal_tol( abserr, 7.775459e-4, 1e-8) );
		CHECK_EQUAL( 315, neval );
		CHECK_EQUAL( 0, ier );
	} /* tdqagi_ */

	double dqawof(double x)
	{
		double ret_val;

		ret_val = 0.;
		if (x > 0.) {
			ret_val = exp(-(x)) * log(x);
		}
		return ret_val;
	} /* dqawof_ */

	TEST(TestDQawo)
	{
		double a, b;
		long ier, last, lenw;
		double work[925];
		long maxp1;
		double omega;
		long neval;
		long leniw, limit, iwork[200];
		double epsabs, abserr;
		long integr;
		double epsrel, result;

		a = 0.;
		b = 1.;
		omega = 10.;
		integr = 1;
		epsabs = 0.;
		epsrel = .001;
		limit = 100;
		leniw = limit << 1;
		maxp1 = 21;
		lenw = (limit << 2) + maxp1 * 25;
		dqawo_(D_fp_fp(dqawof), 
				 &a, &b, &omega, &integr, &epsabs, &epsrel, &result, &
				 abserr, &neval, &ier, &leniw, &maxp1, &lenw, &last, iwork, work);
		if (0) printf("DQAWO: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, -1.776392e-1, 1e-5) );
		CHECK( fp_equal_tol( abserr, 2.816840e-8, 1e-10) );
		CHECK_EQUAL( 255 , neval );
		CHECK_EQUAL( 0 , ier );
	} /* tdqawo_ */

	double dqawff(double x)
	{
		double ret_val=0.0;

		if (x > 0.f) {
			ret_val = sin(x * 50.) / (x * sqrt(x));
		}
		return ret_val;
	} /* dqawff_ */

	TEST(TestDQawf)
	{
		double a;
		long ier, lst;
		long lenw;
		double work[1025];
		long maxp1;
		double omega;
		long neval, leniw, limit, iwork[250];
		double epsabs, abserr;
		long integr, limlst;
		double result;

		a = 0.;
		omega = 8.;
		integr = 2;
		epsabs = .001;
		limlst = 50;
		limit = 100;
		leniw = (limit << 1) + limlst;
		maxp1 = 21;
		lenw = (leniw << 1) + maxp1 * 25;
		dqawf_(D_fp_fp(dqawff), 
				 &a, &omega, &integr, &epsabs, &result, &abserr, &
				 neval, &ier, &limlst, &lst, &leniw, &maxp1, &lenw, iwork, work);
		if (0) printf("DQAWF: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 1.422552, 1e-5) );
		CHECK( fp_equal_tol( abserr, 8.448439e-4, 1e-7) );
		CHECK_EQUAL( 7080, neval );
		CHECK_EQUAL( 0, ier );
	} /* tdqawf_ */

	double dqawsf(double x)
	{
		double ret_val;

		ret_val = sin(x * 10.);
		return ret_val;
	} /* dqawsf_ */

	TEST(TestDQaws)
	{
		double a, b;
		long ier;
		double alfa, beta;
		long last, lenw;
		double work[400];
		long neval, limit;
		long iwork[100];
		double epsabs, abserr;
		long integr;
		double epsrel, result;

		a = 0.;
		b = 1.;
		alfa = -.5;
		beta = -.5;
		integr = 1;
		epsabs = 0.;
		epsrel = .001;
		limit = 100;
		lenw = limit << 2;
		dqaws_(D_fp_fp(dqawsf), 
				 &a, &b, &alfa, &beta, &integr, &epsabs, &epsrel, &
				 result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
		if (0) printf("DQAWS: Result=%e abserr=%e neval=%ld ier=%ld\n",
				 result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, 5.350191e-1, 1e-5) );
		CHECK( fp_equal_tol( abserr, 1.936312e-12, 1e-14) );
		CHECK_EQUAL( 50, neval );
		CHECK_EQUAL( 0, ier );
	} /* tdqaws_ */

	double dqawcf(double x)
	{
		double ret_val;

		ret_val = 1. / (x * x + 1e-4);
		return ret_val;
	} /* dqawcf_ */

	TEST(TestDQawc)
	{
		double a, b, c__;
		long ier;
		long last, lenw;
		double work[400];
		long neval, limit, iwork[100];
		double epsabs, abserr, epsrel, result;

		a = -1.;
		b = 1.;
		c__ = .5;
		epsabs = 0.;
		epsrel = .001;
		limit = 100;
		lenw = limit << 2;
		dqawc_(D_fp_fp(dqawcf), 
				 &a, &b, &c__, &epsabs, &epsrel, &result, &abserr, &
				 neval, &ier, &limit, &lenw, &last, iwork, work);
		if (0) printf("DQAWC: Result=%e abserr=%e neval=%ld ier=%ld\n",
						  result,abserr,neval,ier);
		CHECK( fp_equal_tol( result, -6.284617e2, 1e-3) );
		CHECK( fp_equal_tol( abserr, 1.62477e-1, 1e-3) );
		CHECK( neval <= 255 );
		CHECK_EQUAL( 0, ier );
	} /* tdqawc_ */

}
