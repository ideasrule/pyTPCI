/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "ran.h"
#include "thirdparty.h"

namespace {

	// constants below were calculated with xmaxima 5.17.0 with fpprec set to 30

	TEST(TestFactorial)
	{
		CHECK( fp_equal( factorial(0), 1. ) );
		CHECK( fp_equal( factorial(1), 1. ) );
		CHECK( fp_equal( factorial(2), 2. ) );
		CHECK( fp_equal( factorial(10), 3628800. ) );
		CHECK( fp_equal( factorial(170), 7.257415615307998967396728211129316e306 ) );
		// Validate against standard library function
		CHECK( fp_equal( factorial(170), tgamma(171.), 40 ) );
	}

	TEST(TestLogFactorial)
	{
		CHECK( fp_equal( lfactorial(0), 0. ) );
		CHECK( fp_equal( lfactorial(1), 0. ) );
		CHECK( fp_equal( lfactorial(2), 3.01029995663981195213738894725e-1, 10 ) );
		CHECK( fp_equal( lfactorial(10), 6.55976303287679375117476123996e0, 10 ) );
		CHECK( fp_equal( lfactorial(170), 3.06860781994828320192380273111e2, 10 ) );
		CHECK( fp_equal( lfactorial(1700), 4.75547688279135071178908638597e3, 10 ) );	
		// Validate against standard library function
		CHECK( fp_equal( lfactorial(1700), lgamma(1701.)/log(10.), 10 ) );
	}

	TEST(TestCDGamma)
	{
		complex<double> y = cdgamma( complex<double>(1.,0.) );
		CHECK( fp_equal( y.real(), 1., 10 ) && fp_equal( y.imag(), 0. ) );
		y = cdgamma( complex<double>(2.,0.) );
		CHECK( fp_equal( y.real(), 1., 10 ) && fp_equal( y.imag(), 0. ) );
		y = cdgamma( complex<double>(11.,0.) );
		CHECK( fp_equal( y.real(), 3628800., 50 ) && fp_equal( y.imag(), 0. ) );	
		y = cdgamma( complex<double>(-0.5,0.) );
		CHECK( fp_equal( y.real(), -3.544907701811032054596334966682277e0, 10 ) &&
		       fp_equal( y.imag(), 0. ) );
		y = cdgamma( complex<double>(0.,1.) );
		CHECK( fp_equal( y.real(), -1.549498283018106851249551304838863e-1, 10 ) &&
		       fp_equal( y.imag(), -4.980156681183560427136911174621973e-1, 10 ) );
		y = cdgamma( complex<double>(-1.,-2.) );
		CHECK( fp_equal( y.real(), -3.23612885501927256868232032760969e-2, 30 ) &&
		       fp_equal( y.imag(), -1.122942423463261735043406872030743e-2, 30 ) );		       
	}

	// constants below were taken from Abramowitz & Stegun, Handbook of Mathematical Functions

	TEST(TestBesselJ0)
	{
		CHECK( fp_equal( bessel_j0(0.), 1., 10 ) );
		CHECK( fp_equal_tol( bessel_j0(1.), 0.765197686557967, 1.e-15 ) );
		CHECK( fp_equal_tol( bessel_j0(2.9), -0.224311545791968, 1.e-15 ) );
		CHECK( fp_equal_tol( bessel_j0(4.8), -0.240425327291183, 1.e-15 ) );
		CHECK( fp_equal_tol( bessel_j0(8.3), 0.096006100895010, 1.e-15 ) );
		CHECK( fp_equal_tol( bessel_j0(17.4), -0.118955856336348, 1.e-15 ) );
	}

	TEST(TestBesselJ1)
	{
		CHECK( fp_equal( bessel_j1(0.), 0., 10 ) );
		CHECK( fp_equal( bessel_j1(1.e-30), 5.e-31, 10 ) );
		CHECK( fp_equal_tol( bessel_j1(1.e-15), 5.e-16, 1.e-27 ) );
		CHECK( fp_equal_tol( bessel_j1(1.), 0.4400505857, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_j1(2.9), 0.3754274818, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_j1(4.8), -0.2984998581, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_j1(8.3), 0.2657393020, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_j1(17.4), -0.1532161760, 1.e-10 ) );
	}

	TEST(TestBesselJn)
	{
		CHECK( fp_equal( bessel_jn(2,0.), 0., 10 ) );
		CHECK( fp_equal( bessel_jn(2,1.e-30), 1.25e-61, 10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,2.e-15), 5.e-31, 1.e-42 ) );
		CHECK( fp_equal_tol( bessel_jn(2,2.e-10), 5.e-21, 1.e-32 ) );
		CHECK( fp_equal_tol( bessel_jn(2,0.099999999999), 0.0012489587, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,0.100000000001), 0.0012489587, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,1.), 0.1149034849, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,2.9), 0.4832270505, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,4.8), 0.1160503864, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,8.3), -0.0319725341, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_jn(2,17.4), 0.1013448016, 1.e-10 ) );

		CHECK( fp_equal( bessel_jn(8,0.), 0., 10 ) );
		CHECK( fp_equal( bessel_jn(8,2.e-20), 1.e-160/40320., 10 ) );
		CHECK( fp_equal_tol( bessel_jn(8,2.e-15), 1.e-120/40320., 1.e-136 ) );
		CHECK( fp_equal_tol( bessel_jn(8,1.), 9.4223e-8, 1.e-12 ) );
		CHECK( fp_equal_tol( bessel_jn(8,2.8), 2.9367e-4, 1.e-8 ) );
		CHECK( fp_equal_tol( bessel_jn(8,4.8), 1.4079e-2, 1.e-6 ) );
		CHECK( fp_equal_tol( bessel_jn(8,8.2), 0.24257, 1.e-5 ) );

		CHECK( fp_equal_tol( bessel_jn(20,1.), 3.873503009e-25, 1.e-34 ) );
		CHECK( fp_equal_tol( bessel_jn(50,1.), 2.906004948e-80, 1.e-89 ) );
		CHECK( fp_equal_tol( bessel_jn(100,1.), 8.431828790e-189, 1.e-198 ) );
		CHECK( fp_equal_tol( bessel_jn(100,100.), 9.636667330e-2, 1.e-11 ) );
	}

	TEST(TestBesselY0)
	{
		CHECK( fp_equal_tol( bessel_y0(1.e-30), -44.0499402279, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y0(1.), 0.0882569642, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y0(2.9), 0.4079117692, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y0(4.8), -0.2723037945, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y0(8.3), 0.2595149638, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y0(17.4), -0.1497391883, 1.e-10 ) );
	}

	TEST(TestBesselY1)
	{
		CHECK( fp_equal_tol( bessel_y1(1.e-30), -6.36619772368e29, 1.e18 ) );
		CHECK( fp_equal_tol( bessel_y1(1.), -0.7812128213, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y1(2.9), 0.2959400546, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y1(4.8), 0.2135651673, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y1(8.3), -0.0805975035, 1.e-10 ) );
		CHECK( fp_equal_tol( bessel_y1(17.4), 0.1147053859, 1.e-10 ) );
	}

	TEST(TestBesselYn)
	{
		CHECK( fp_equal_tol( bessel_yn(2,1.e-30), -1.27323954474e60, 1.e49 ) );
		CHECK( fp_equal_tol( bessel_yn(2,1.), -1.65068261, 1.e-8 ) );
		CHECK( fp_equal_tol( bessel_yn(2,2.9), -0.20381518, 1.e-8 ) );
		CHECK( fp_equal_tol( bessel_yn(2,4.8), 0.36128928, 1.e-8 ) );
		CHECK( fp_equal_tol( bessel_yn(2,8.3), -0.27893605, 1.e-8 ) );
		CHECK( fp_equal_tol( bessel_yn(2,17.4), 0.16292372, 1.e-8 ) );

		CHECK( fp_equal_tol( bessel_yn(8,2.e-20), -1.60428182637e163, 1.e152 ) );
		CHECK( fp_equal_tol( bessel_yn(8,1.), -4.2567e5, 1.e1 ) );
		CHECK( fp_equal_tol( bessel_yn(8,2.8), -1.4486e2, 1.e-2 ) );
		CHECK( fp_equal_tol( bessel_yn(8,4.8), -3.5855, 1.e-4 ) );
		CHECK( fp_equal_tol( bessel_yn(8,8.2), -0.35049, 1.e-5 ) );

		CHECK( fp_equal_tol( bessel_yn(20,1.), -4.113970315e22, 1.e13 ) );
		CHECK( fp_equal_tol( bessel_yn(50,1.), -2.191142813e77, 1.e68 ) );
		CHECK( fp_equal_tol( bessel_yn(100,1.), -3.775287810e185, 1.e176 ) );
		CHECK( fp_equal_tol( bessel_yn(100,100.), -1.669214114e-1, 1.e-10 ) );
	}

	// constants below were calculated using http://keisan.casio.com/exec/system/1180573473

	TEST(TestBesselI0)
	{
		CHECK( fp_equal( bessel_i0(0.), 1. ) );
		CHECK( fp_equal( bessel_i0(1.), 1.266065877752008335598 ) );
		CHECK( fp_equal( bessel_i0(2.9), 4.50274866132627436631 ) );
		CHECK( fp_equal( bessel_i0(4.8), 22.79367799310579796012 ) );
		CHECK( fp_equal( bessel_i0(8.3), 566.255055846491838162, 10 ) );
		CHECK( fp_equal( bessel_i0(17.4), 3471961.712717317753354, 10 ) );

		CHECK( fp_equal( bessel_i0(-1.), 1.266065877752008335598 ) );
		CHECK( fp_equal( bessel_i0(-17.4), 3471961.712717317753354, 10 ) );

		CHECK( fp_equal( bessel_i0_scaled(1.), exp(-1.)*bessel_i0(1.) ) );
		CHECK( fp_equal( bessel_i0_scaled(17.4), exp(-17.4)*bessel_i0(17.4) ) );

		CHECK( fp_equal( bessel_i0_scaled(-1.), exp(-1.)*bessel_i0(-1.) ) );
		CHECK( fp_equal( bessel_i0_scaled(-17.4), exp(-17.4)*bessel_i0(-17.4) ) );
	}

	TEST(TestBesselI1)
	{
		CHECK( fp_equal( bessel_i1(0.), 0. ) );
		CHECK( fp_equal( bessel_i1(1.), 0.5651591039924850272077 ) );
		CHECK( fp_equal( bessel_i1(2.9), 3.6126072124369077367 ) );
		CHECK( fp_equal( bessel_i1(4.8), 20.25283460023855998949 ) );
		CHECK( fp_equal( bessel_i1(8.3), 530.959765452190423905, 7 ) );
		CHECK( fp_equal( bessel_i1(17.4), 3370668.4094771099123, 10 ) );

		CHECK( fp_equal( bessel_i1(-1.), -0.5651591039924850272077 ) );
		CHECK( fp_equal( bessel_i1(-17.4), -3370668.4094771099123, 10 ) );

		CHECK( fp_equal( bessel_i1_scaled(1.), exp(-1.)*bessel_i1(1.) ) );
		CHECK( fp_equal( bessel_i1_scaled(17.4), exp(-17.4)*bessel_i1(17.4) ) );

		CHECK( fp_equal( bessel_i1_scaled(-1.), exp(-1.)*bessel_i1(-1.) ) );
		CHECK( fp_equal( bessel_i1_scaled(-17.4), exp(-17.4)*bessel_i1(-17.4) ) );
	}

	TEST(TestBesselI0I1)
	{
		double i0val, i1val;
		bessel_i0_i1(0.3, &i0val, &i1val);
		CHECK( fp_equal( i0val, 1.02262687935159699112 ) );
		CHECK( fp_equal( i1val, 0.1516938400035927803287 ) );
		bessel_i0_i1_scaled(0.3, &i0val, &i1val);
		CHECK( fp_equal( i0val, 1.02262687935159699112*exp(-0.3) ) );
		CHECK( fp_equal( i1val, 0.1516938400035927803287*exp(-0.3) ) );

		bessel_i0_i1(-0.3, &i0val, &i1val);
		CHECK( fp_equal( i0val, 1.02262687935159699112 ) );
		CHECK( fp_equal( i1val, -0.1516938400035927803287 ) );
		bessel_i0_i1_scaled(-0.3, &i0val, &i1val);
		CHECK( fp_equal( i0val, 1.02262687935159699112*exp(-0.3) ) );
		CHECK( fp_equal( i1val, -0.1516938400035927803287*exp(-0.3) ) );

		bessel_i0_i1(22., &i0val, &i1val);
		CHECK( fp_equal( i0val, 306692993.6403647456135 ) );
		CHECK( fp_equal( i1val, 299639606.8773789797263 ) );
		bessel_i0_i1_scaled(22., &i0val, &i1val);
		CHECK( fp_equal( i0val, 306692993.6403647456135*exp(-22.) ) );
		CHECK( fp_equal( i1val, 299639606.8773789797263*exp(-22.) ) );

		bessel_i0_i1(-22., &i0val, &i1val);
		CHECK( fp_equal( i0val, 306692993.6403647456135 ) );
		CHECK( fp_equal( i1val, -299639606.8773789797263 ) );
		bessel_i0_i1_scaled(-22., &i0val, &i1val);
		CHECK( fp_equal( i0val, 306692993.6403647456135*exp(-22.) ) );
		CHECK( fp_equal( i1val, -299639606.8773789797263*exp(-22.) ) );
	}

	TEST(TestBesselK0)
	{
		CHECK_THROW( bessel_k0(0.), domain_error );
		CHECK_THROW( bessel_k0(-1.), domain_error );
		CHECK( fp_equal( bessel_k0(2.e-30), 68.50033712491983765993 ) );
		CHECK( fp_equal( bessel_k0(1.), 0.4210244382407083333356 ) );
		CHECK( fp_equal( bessel_k0(2.9), 0.0390062345662234241006 ) );
		CHECK( fp_equal( bessel_k0(4.8), 0.00459724631672465789944 ) );
		CHECK( fp_equal( bessel_k0(8.3), 1.065830600438036175396e-4, 5 ) );
		CHECK( fp_equal( bessel_k0(17.4), 8.279919509749720001621e-9, 10 ) );

		CHECK( fp_equal( bessel_k0_scaled(1.), exp(1.)*bessel_k0(1.) ) );
		CHECK( fp_equal( bessel_k0_scaled(17.4), exp(17.4)*bessel_k0(17.4) ) );
	}

	TEST(TestBesselK1)
	{
		CHECK_THROW( bessel_k1(0.), domain_error );
		CHECK_THROW( bessel_k1(-1.), domain_error );
		CHECK( fp_equal( bessel_k1(2.e-30), 5.e29 ) );
		CHECK( fp_equal( bessel_k1(1.), 0.6019072301972345747375 ) );
		CHECK( fp_equal( bessel_k1(2.9), 0.0452864232983614435612 ) );
		CHECK( fp_equal( bessel_k1(4.8), 0.00505517644405629981585 ) );
		CHECK( fp_equal( bessel_k1(8.3), 1.128300939464441907577e-4 ) );
		CHECK( fp_equal( bessel_k1(17.4), 8.514610381504642118954e-9, 10 ) );

		CHECK( fp_equal( bessel_k1_scaled(1.), exp(1.)*bessel_k1(1.) ) );
		CHECK( fp_equal( bessel_k1_scaled(17.4), exp(17.4)*bessel_k1(17.4) ) );
	}

	TEST(TestBesselK0K1)
	{
		double k0val, k1val;
		CHECK_THROW( bessel_k0_k1(0., &k0val, &k1val), domain_error );
		CHECK_THROW( bessel_k0_k1(-1., &k0val, &k1val), domain_error );
		bessel_k0_k1(0.3, &k0val, &k1val);
		CHECK( fp_equal( k0val, 1.372460060544297376645 ) );
		CHECK( fp_equal( k1val, 3.055992033457324978851 ) );
		bessel_k0_k1_scaled(0.3, &k0val, &k1val);
		CHECK( fp_equal( k0val, 1.372460060544297376645*exp(0.3) ) );
		CHECK( fp_equal( k1val, 3.055992033457324978851*exp(0.3) ) );

		bessel_k0_k1(22., &k0val, &k1val);
		CHECK( fp_equal( k0val, 7.412351614084865368946e-11 ) );
		CHECK( fp_equal( k1val, 7.578981163485331089804e-11 ) );
		bessel_k0_k1_scaled(22., &k0val, &k1val);
		CHECK( fp_equal( k0val, 7.412351614084865368946e-11*exp(22.) ) );
		CHECK( fp_equal( k1val, 7.578981163485331089804e-11*exp(22.) ) );
	}

	// constants below were taken from Abramowitz & Stegun, Handbook of Mathematical Functions

	TEST(TestEllpk)
	{
		CHECK( fp_equal( ellpk(1.), PI/2., 10 ) );
		CHECK( fp_equal_tol( ellpk(0.86), 1.630575548881754, 1.e-15 ) );
		CHECK( fp_equal_tol( ellpk(0.56), 1.806327559107699, 1.e-15 ) );
		CHECK( fp_equal_tol( ellpk(0.36), 1.995302777664729, 1.e-15 ) );
		CHECK( fp_equal_tol( ellpk(0.13), 2.455338028321380, 1.e-15 ) );
		CHECK( fp_equal_tol( ellpk(0.01), 3.695637362989875, 1.e-15 ) );
	}

	// constants below were calculated with xmaxima 5.34.1 with fpprec set to 30

	TEST(TestIgam)
	{
		CHECK( fp_equal_tol( igam(0.3,0.2), 0.65750672426972171, 1.e-15 ) );
		CHECK( fp_equal_tol( igam(0.46,0.8), 0.81316729720991554, 1.e-15 ) );
		CHECK( fp_equal_tol( igam(0.67,2.2), 0.94307235949930435, 1.e-15 ) );
		CHECK( fp_equal_tol( igam(1.2,4.5), 0.98302397308764477, 1.e-15 ) );
		CHECK( fp_equal_tol( igam(1.,10.), 0.99995460007023752, 1.e-15 ) );
		CHECK( fp_equal_tol( igam(0.3,1e-5), 0.03523536061556257, 1.e-15 ) );
	}

	TEST(TestIgamc)
	{
		CHECK( fp_equal_tol( igamc(0.3,0.2), 1.-0.65750672426972171, 1.e-15 ) );
		CHECK( fp_equal_tol( igamc(0.46,0.8), 1.-0.81316729720991554, 1.e-15 ) );
		CHECK( fp_equal_tol( igamc(0.67,2.2), 1.-0.94307235949930435, 1.e-15 ) );
		CHECK( fp_equal_tol( igamc(1.2,4.5), 1.-0.98302397308764477, 1.e-15 ) );
		CHECK( fp_equal_tol( igamc(1.,10.), 1.-0.99995460007023752, 1.e-15 ) );
		CHECK( fp_equal_tol( igamc(0.3,1e-5), 1.-0.03523536061556257, 1.e-15 ) );
	}

	TEST(TestIgamc_scaled)
	{
		CHECK( fp_equal_tol( igamc_scaled(0.3,0.2), 0.41832223162827345394, 1.e-15 ) );
		CHECK( fp_equal_tol( igamc_scaled(0.46,0.8), 0.41580382684020182212, 1.e-15 ) );
		CHECK( fp_equal_tol( igamc_scaled(0.67,2.2), 0.51377272400971086105, 1.e-15 ) );
		CHECK( fp_equal_tol( igamc_scaled(1.2,4.5), 1.52813324353067303945, 1.5e-15 ) );
		CHECK( fp_equal_tol( igamc_scaled(1.,10.), 1., 1.e-15 ) );
		CHECK( fp_equal_tol( igamc_scaled(0.3,1e5), 1.057055858519714446e-4 , 1.e-15 ) );
	}

	// constants below were taken from Abramowitz & Stegun, Handbook of Mathematical Functions

	TEST(TestExpn)
	{
		CHECK( fp_equal_tol( expn(1,0.25), 0.25*0.9408157528 - log(0.25) - EULER, 1.e-10 ) );
		CHECK( fp_equal_tol( expn(1,0.79), 0.316277004, 1.e-9 ) );
		CHECK( fp_equal_tol( expn(1,1.97), 0.050976988, 1.e-9 ) );

		CHECK( fp_equal_tol( expn(2,0.25), 0.8643037 + log(0.25)/4., 1.e-7 ) );
		CHECK( fp_equal_tol( expn(2,0.79), 0.2039860, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(2,1.97), 0.0390322, 1.e-7 ) );

		CHECK( fp_equal( expn(3,0.), 0.5, 10 ) );
		CHECK( fp_equal_tol( expn(3,0.15), 0.3822761, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(3,0.79), 0.1463479, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(3,1.97), 0.0312817, 1.e-7 ) );

		CHECK( fp_equal( expn(4,0.), 1./3., 10 ) );
		CHECK( fp_equal_tol( expn(4,0.15), 0.2677889, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(4,0.79), 0.1127433, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(4,1.97), 0.0259440, 1.e-7 ) );

		CHECK( fp_equal( expn(10,0.), 1./9., 10 ) );
		CHECK( fp_equal_tol( expn(10,0.15), 0.0938786, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(10,0.79), 0.0459453, 1.e-7 ) );
		CHECK( fp_equal_tol( expn(10,1.97), 0.0124964, 1.e-7 ) );
	}

	TEST(TestE1)
	{
		CHECK( fp_equal_tol( e1(0.25), 0.25*0.9408157528 - log(0.25) - EULER, 2.e-7 ) );
		CHECK( fp_equal_tol( e1(0.79), 0.316277004, 2.e-7 ) );
		CHECK( fp_equal_tol( e1(1.97), 0.050976988, 1.4e-9 ) );

		CHECK( fp_equal_tol( e1_scaled(0.25), (0.25*0.9408157528 - log(0.25) - EULER)*exp(0.25), 2.6e-7 ) );
		CHECK( fp_equal_tol( e1_scaled(0.79), 0.316277004*exp(0.79), 4.4e-7 ) );
		CHECK( fp_equal_tol( e1_scaled(1.97), 0.050976988*exp(1.97), 1.e-8 ) );
	}

	TEST(TestE2)
	{
		CHECK( fp_equal_tol( e2(0.25), 0.8643037 + log(0.25)/4., 1.e-7 ) );
		CHECK( fp_equal_tol( e2(0.79), 0.2039860, 1.e-7 ) );
		CHECK( fp_equal_tol( e2(1.97), 0.0390322, 1.e-7 ) );
	}

	TEST(TestErf)
	{
		/* constants calculated with xmaxima 5.22.1 */
		CHECK( fp_equal_tol( erf(0.), 0., 1.e-22 ) );
		// erf(x) loses some precision around 1.e-10, but should still be plenty good...
		CHECK( fp_equal_tol( erf(1.e-10), 1.1283791671081724525e-10, 3.e-21 ) );
		CHECK( fp_equal_tol( erf(1.e-5), 1.1283791670579000307e-5, 1.e-17 ) );
		CHECK( fp_equal_tol( erf(0.1), 1.1246291601828489259e-1, 1.e-13 ) );
		CHECK( fp_equal_tol( erf(0.5), 5.2049987781304653768e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erf(1.), 8.4270079294971486934e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erf(2.), 9.9532226501895273416e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erf(10.), 1.0, 1.e-12 ) );
		CHECK( fp_equal( erf(-1.), -erf(1.) ) );
	}

	TEST(TestErfc)
	{
		/* constants calculated with xmaxima 5.22.1 */
		CHECK( fp_equal_tol( erfc(-1.), 1.8427007929497148693e0, 1.e-12 ) );
		CHECK( fp_equal_tol( erfc(0.), 1., 1.e-12 ) );
		CHECK( fp_equal_tol( erfc(1.e-5), 9.99988716208329421e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erfc(0.1), 8.8753708398171510741e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erfc(0.5), 4.7950012218695346232e-1, 1.e-12 ) );
		CHECK( fp_equal_tol( erfc(1.), 1.5729920705028513066e-1, 1.e-13 ) );
		CHECK( fp_equal_tol( erfc(2.), 4.6777349810472658396e-3, 1.e-15 ) );
		CHECK( fp_equal_tol( erfc(10.), 2.088487583762544757e-45, 1.e-57 ) );
	}

	TEST(TestErfce)
	{
		/* constants taken from Finn G.D. & Mugglestone D., 1965, MNRAS 129, 221 */
		/* this is the voigt function H(a,0) normalized according to 9-44 of Mihalas */
		CHECK( fp_equal_tol( erfce(0.), 1., 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.01), 9.88815e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.02), 9.77826e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.05), 9.45990e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.10), 8.96456e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.20), 8.09019e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(0.55), 5.90927e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( erfce(1.00), 4.27583e-1, 1.e-6 ) );
	}

	TEST(TestVoigtH0)
	{
		/* constants taken from Finn G.D. & Mugglestone D., 1965, MNRAS 129, 221 */
		/* this is the voigt function H(a,0) normalized according to 9-44 of Mihalas */
		CHECK( fp_equal_tol( VoigtH0(0.), 1., 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.01), 9.88815e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.02), 9.77826e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.05), 9.45990e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.10), 8.96456e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.20), 8.09019e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(0.55), 5.90927e-1, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtH0(1.00), 4.27583e-1, 1.e-6 ) );
	}

	TEST(TestVoigtU0)
	{
		/* constants taken from Finn G.D. & Mugglestone D., 1965, MNRAS 129, 221 */
		/* this is the voigt function U(a,0) normalized according to 9-45 of Mihalas */
		CHECK( fp_equal_tol( VoigtU0(0.), 1./SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.01), 9.88815e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.02), 9.77826e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.05), 9.45990e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.10), 8.96456e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.20), 8.09019e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(0.55), 5.90927e-1/SQRTPI, 1.e-6 ) );
		CHECK( fp_equal_tol( VoigtU0(1.00), 4.27583e-1/SQRTPI, 1.e-6 ) );
	}

	TEST(TestGegen)
	{
		/* constants calculated with xmaxima 5.30.0 */
		CHECK( fp_equal_tol( gegenbauer(0,1,1.), 1., 1.e-12 ) );
		CHECK( fp_equal_tol( gegenbauer(0,1,0.5), 1., 1.e-12 ) );
		CHECK( fp_equal_tol( gegenbauer(1,1,0.5), 1., 1.e-12 ) );	
		CHECK( fp_equal_tol( gegenbauer(2,2,0.55), 1.63, 1.e-12 ) );
		CHECK( fp_equal_tol( gegenbauer(2,0.2,0.5), -0.08, 1.e-12 ) );
		CHECK( fp_equal_tol( gegenbauer(2,3,0.553), 4.339416, 1.e-12 ) );
		CHECK( fp_equal_tol( gegenbauer(3,3,0.553), 0.25699016, 1.e-12 ) );
		CHECK( fp_equal_tol( gegenbauer(3,3,-0.553), -0.25699016, 1.e-12 ) );
		// test case from Levrie & Piessens, http://nalag.cs.kuleuven.be/papers/paul/report74/report74.ps.gz, corrected using xmaxima, xmaxima quoted error is 3e-13
		CHECK( fp_equal_tol( gegenbauer(40,0.5,0.999), 0.333404687422222, 1.e-9 ) );
		// test case from http://nalag.cs.kuleuven.be/papers/paul/report74/report74.ps.gz, corrected using xmaxima, xmaxima quoted error is 3e-13
		CHECK( fp_equal_tol( gegenbauer(40,0.5,0.999), 0.333404687422222, 1.e-9 ) );
		{
			UltraGen u(1,1.);
			CHECK( fp_equal_tol( u.val(), 1., 1.e-12 ) );
		}		
		{
			UltraGen u(1,0.5);
			CHECK( fp_equal_tol( u.val(), 1., 1.e-12 ) );
		}
		{
			UltraGen u(2,0.5);
			u.step();
			CHECK( fp_equal_tol( u.val(), 1., 1.e-12 ) );
		}
		{
			UltraGen u(4,0.55);
			u.step(); u.step();
			CHECK( fp_equal_tol( u.val(), 1.63, 1.e-12 ) );
		}
		{
			UltraGen u(5,0.553);
			u.step(); u.step();
			CHECK( fp_equal_tol( u.val(), 4.339416, 1.e-12 ) );
		}
		{
			UltraGen u(6,0.553);
			u.step(); u.step(); u.step();
			CHECK( fp_equal_tol( u.val(), 0.25699016, 1.e-12 ) );
		}
		CHECK( fp_equal_tol( gegenbauer(3,3,-0.553), -0.25699016, 1.e-12 ) );
		
	}

	TEST(TestSixj)
	{
		long a=1,b=2,c=3,s=a+b+c;
		double cc = ipow(-1,s)/sqrt(double((2*b+1)*(2*c+1)));
		CHECK( fp_equal_tol( sjs(2*a,2*b,2*c,0,2*c,2*b), cc, 1.e-15 ) );

		// test failure modes
		// this set of j's fails only Triangle2( j1, j2, j3 )
		CHECK( sjs( 0, 6, 2, 4, 2, 2 ) == 0. );
		// this set of j's fails only Triangle2( j4, j5, j3 )
		CHECK( sjs( 2, 2, 0, 0, 4, 2 ) == 0. );
		// this set of j's fails only Triangle2( j1, j5, j6 )
		CHECK( sjs( 2, 4, 6, 6, 2, 6 ) == 0. );
		// this set of j's fails only Triangle2( j2, j4, j6 )
		CHECK( sjs( 2, 6, 4, 0, 4, 2 ) == 0. );

		// test sjs() for small arguments by comparison with SixJFull()
		for( int j1=0; j1 < 5; ++j1 )
		for( int j2=0; j2 < 5; ++j2 )
		for( int j3=0; j3 < 5; ++j3 )
		for( int j4=0; j4 < 5; ++j4 )
		for( int j5=0; j5 < 5; ++j5 )
		for( int j6=0; j6 < 5; ++j6 )
		{
			double s1 = sjs( j1, j2, j3, j4, j5, j6 );
			if( !Triangle2( j1, j2, j3 ) || !Triangle2( j4, j5, j3 ) ||
			    !Triangle2( j1, j5, j6 ) || !Triangle2( j2, j4, j6 ) )
			{
				CHECK( s1 == 0. );
			}
			else
			{
				double s2 = SixJFull( j1, j2, j3, j4, j5, j6 );
				CHECK( fp_equal_tol( s1, s2, 1.e-15 ) );
			}
		}

		// test a few handpicked 6j symbols with large arguments
		// results were taken from
		// http://www-stone.ch.cam.ac.uk/cgi-bin/wigner.cgi
		double s1 = sjs( 42, 38, 34, 46, 36, 32 );
		double s2 = -146218147./1778945129108.*sqrt(20735./57.);
		CHECK( fp_equal_tol( s1, s2, 1.e-15 ) );
		s1 = sjs( 46, 44, 42, 40, 38, 36 );
		s2 = -16185716849919./4093621176253.*sqrt(201./135605470.);
		CHECK( fp_equal_tol( s1, s2, 1.e-15 ) );
		s1 = sjs( 43, 39, 40, 12, 48, 49 );
		s2 = 1111./11600540.*sqrt(2982639./595.);
		CHECK( fp_equal_tol( s1, s2, 1.e-15 ) );

		// test sjs() for random arguments by comparison with SixJFull()
		for( int i=0; i < 1000; ++i )
		{
			int j1, j2, j3, j4, j5, j6;
			const int JMAX = 110;
			do
			{
				j1 = ran.u16()%JMAX;
				j2 = ran.u16()%JMAX;
				j3 = ran.u16()%JMAX;
				j4 = ran.u16()%JMAX;
				j5 = ran.u16()%JMAX;
				j6 = ran.u16()%JMAX;
			}
			while( !Triangle2( j1, j2, j3 ) || !Triangle2( j4, j5, j3 ) ||
			       !Triangle2( j1, j5, j6 ) || !Triangle2( j2, j4, j6 ) );

			s1 = sjs( j1, j2, j3, j4, j5, j6 );
			s2 = SixJFull( j1, j2, j3, j4, j5, j6 );
			// for very large arguments, cancellation errors start growing
			// so we need to increase the tolerance to 3e-12 here...
			CHECK( fp_equal_tol( s1, s2, 3.e-12 ) );
		}
	}
	TEST(TestRecSixj)
	{
		const int NPT=256;
		long a=1,b=2,c=3,s=a+b+c;
		double cc = pow(-1,s)/sqrt((2*b+1)*(2*c+1));
		double sixcof[NPT];
		long ier;
		double l1min, l1max, lmatch;
		rec6j(sixcof, b, c, 0, c, b, &l1min, &l1max, &lmatch, NPT, &ier);
		CHECK( fp_equal_tol( sixcof[a-(int)l1min], cc, 1.e-15 ) );

		// test failure modes
		// this set of j's fails only Triangle2( j1, j2, j3 )
		rec6j(sixcof, 3, 1, 2, 1, 1, &l1min, &l1max, &lmatch, NPT, &ier);
		CHECK( 0 < l1min );
		// this set of j's fails only Triangle2( j4, j5, j3 )
		rec6j(sixcof, 1, 0, 0, 2, 1, &l1min, &l1max, &lmatch, NPT, &ier);
		CHECK( ier == -1 );
		// this set of j's fails only Triangle2( j1, j5, j6 )
		rec6j(sixcof, 1, 3, 3, 1, 3, &l1min, &l1max, &lmatch, NPT, &ier);
		CHECK( 1 < l1min );
		// this set of j's fails only Triangle2( j2, j4, j6 )
		rec6j(sixcof, 3, 2, 0, 2, 1, &l1min, &l1max, &lmatch, NPT, &ier);
		CHECK( ier == -1 );

		// test rec6j() for small arguments by comparison with SixJFull()
		for( int j2=0; j2 < 5; ++j2 )
		for( int j3=0; j3 < 5; ++j3 )
		for( int j4=0; j4 < 5; ++j4 )
		for( int j5=0; j5 < 5; ++j5 )
		for( int j6=0; j6 < 5; ++j6 )
		{
			rec6j(sixcof, 0.5*j2, 0.5*j3, 0.5*j4, 0.5*j5, 0.5*j6, &l1min, &l1max, &lmatch, NPT, &ier);
			for( int j1=0; j1 < 5; ++j1 )
			{
				double pos = (0.5*j1-l1min);
				if( !Triangle2( j1, j2, j3 ) || !Triangle2( j4, j5, j3 ) ||
					 !Triangle2( j1, j5, j6 ) || !Triangle2( j2, j4, j6 ) )
				{
					CHECK( ier == -1 || pos != int(pos) || 
							 0.5*j1 < l1min || 0.5*j1 > l1max );
				}
				else
				{
					CHECK( ier >= 0 );
					double s1 = sixcof[int(pos)];
					double s2 = SixJFull( j1, j2, j3, j4, j5, j6 );
					CHECK( fp_equal_tol( s1, s2, 1.e-15 ) );
				}
			}
		}

		// test a few handpicked 6j symbols with large arguments
		// results were taken from
		// http://www-stone.ch.cam.ac.uk/cgi-bin/wigner.cgi
		rec6j(sixcof, 19, 17, 23, 18, 16, &l1min, &l1max, &lmatch, NPT, &ier);
		double s1 = sixcof[int(21-l1min)];
		double s2 = -146218147./1778945129108.*sqrt(20735./57.);
		CHECK( fp_equal_tol( s1, s2, 1.e-15 ) );
		rec6j(sixcof, 22, 21, 20, 19, 18, &l1min, &l1max, &lmatch, NPT, &ier);
		s1 = sixcof[int(23-l1min)];
		s2 = -16185716849919./4093621176253.*sqrt(201./135605470.);
		CHECK( fp_equal_tol( s1, s2, 1.e-15 ) );
		rec6j(sixcof, 19.5, 20, 6, 24, 24.5, &l1min, &l1max, &lmatch, NPT, &ier);
		s1 = sixcof[int(21.5-l1min)];
		s2 = 1111./11600540.*sqrt(2982639./595.);
		CHECK( fp_equal_tol( s1, s2, 1.e-15 ) );

		// test rec6j() for random arguments by comparison with SixJFull()
		for( int i=0; i < 1000; ++i )
		{
			int j1, j2, j3, j4, j5, j6;
			const int JMAX = 110;
			do
			{
				j1 = ran.u16()%JMAX;
				j2 = ran.u16()%JMAX;
				j3 = ran.u16()%JMAX;
				j4 = ran.u16()%JMAX;
				j5 = ran.u16()%JMAX;
				j6 = ran.u16()%JMAX;
			}
			while( !Triangle2( j1, j2, j3 ) || !Triangle2( j4, j5, j3 ) ||
			       !Triangle2( j1, j5, j6 ) || !Triangle2( j2, j4, j6 ) );

			rec6j(sixcof, 0.5*j2, 0.5*j3, 0.5*j4, 0.5*j5, 0.5*j6, 
					&l1min, &l1max, &lmatch, NPT, &ier);
			s1 = sixcof[int(0.5*j1-l1min)];
			s2 = SixJFull( j1, j2, j3, j4, j5, j6 );
			// for very large arguments, cancellation errors start growing
			// for SixJFull, so we need to increase the tolerance to 3e-11 here...
			CHECK( fp_equal_tol( s1, s2, 3.e-11 ) );
		}

		// test rec6j for problematic arguments
		// comparison values were calculated using a 512-bit version of sjs()
		rec6j(sixcof, 85., 84., 89.5, 89.5, 89.5, &l1min, &l1max, &lmatch, NPT, &ier);
		s1 = sixcof[int(86.-l1min)];
		CHECK( fp_equal_tol( s1, 3.3988213869175631e-08, 4.e-18 ) );
		// check endpoints of recursion
		CHECK( fp_equal( l1min, 1. ) );
		CHECK( fp_equal( l1max, 169. ) );
		s1 = sixcof[int(1.-l1min)];
		CHECK( fp_equal_tol( s1, 0.0035632864755963242, 4.e-18 ) );
		s1 = sixcof[int(169.-l1min)];
		CHECK( fp_equal_tol( s1, 1.785014472974611e-17, 2.e-31 ) );
	}

	TEST(TestVoigtH)
	{
		// check that the Voigt profile returned by VoigtH() is properly normalized
		const int NP = 200;
		realnum v[NP], a, y[NP];
		// for a > 0.1, VoigtH() calls humlik(), for smaller values it
		// calls FastVoigtH().
		//
		// humlik() is set up for a relative precision of 1e-4 over its
		// entire range, but looks to be more precise in practice (at
		// least at v=0)
		// FastVoigtH() is set up for a rel. precision of 2.5e-3 over its
		// entire range, but should be better than 1e-4 for a < 0.0235
		a = 0.0002_r;
		for( int i=0; i < 9; ++i )
		{
			// test both humlik() and FastVoigtH()
			for( int i=0; i < NP; ++i )
				v[i] = realnum(i)*max(a,1_r)/5_r;
			VoigtH( a, v, y, NP );
			realnum integral = 0_r;
			// We need the integral from -infinity to +infinity, which is simply
			// two times the integral from 0 to +infinity. Hence we omit the
			// division by 2 in the trapezoidal rule
			for( int i=1; i < NP; ++i )
				integral += (y[i]+y[i-1])*(v[i]-v[i-1]);
			// add in the integral over the Lorentz wings assuming H(a,v) = c/v^2
			integral += 2_r*v[NP-1]*y[NP-1];
			// VoigtH() calculates H(a,v), so integral should be sqrt(pi)
			CHECK( fp_equal_tol( integral, realnum(SQRTPI), 1.e-4_r ) );
			// also check the central value...
			realnum h0 = realnum(VoigtH0(a));
			CHECK( fp_equal_tol( y[0], h0, 1.e-4_r*h0 ) );
			a *= 10_r;
		}
		// now do some spot checks
		// constants taken from
		// >>refer	Zaghloul M.R. & Ali A.N. 2011, ACM Transactions on Mathematical Software, 38, 15
		// available from http://arxiv.org/abs/1106.0151
		v[0] = 5.76_r;
		VoigtH( 1.e-20_r, v, y, 1 );
		// this constant comes from page 21 ("Present algorithm").
		CHECK( fp_equal_tol( y[0], 3.900779639194698e-015_r, 4.e-19_r ) );
		v[0] = 6.3e-2_r;
		v[1] = 6.3e+0_r;
		v[2] = 6.3e+2_r;
		VoigtH( 1.e-20_r, v, y, 3 );
		// these constants come from Table 2 of the same paper ("Present algorithm").
		CHECK( fp_equal_tol( y[0], 9.960388660702479e-001_r, 1.e-4_r ) );
		CHECK( fp_equal_tol( y[1], 5.792460778844116e-018_r, 6.e-22_r ) );
		CHECK( fp_equal_tol( y[2], 1.421495882582395e-026_r, 1.4e-30_r ) );
		VoigtH( 1.e-14_r, v, y, 3 );
		CHECK( fp_equal_tol( y[0], 9.960388660702366e-001_r, 1.e-4_r ) );
		CHECK( fp_equal_tol( y[1], 1.536857621303163e-016_r, 1.5e-20_r ) );
		CHECK( fp_equal_tol( y[2], 1.421495882582395e-020_r, 1.4e-24_r ) );
		VoigtH( 1.e-12_r, v, y, 3 );
		CHECK( fp_equal_tol( y[0], 9.960388660691284e-001_r, 1.e-4_r ) );
		CHECK( fp_equal_tol( y[1], 1.479513723737753e-014_r, 1.5e-18_r ) );
		CHECK( fp_equal_tol( y[2], 1.421495882582395e-018_r, 1.4e-22_r ) );
		VoigtH( 1.e-10_r, v, y, 3 );
		CHECK( fp_equal_tol( y[0], 9.960388659583033e-001_r, 1.e-4_r ) );
		CHECK( fp_equal_tol( y[1], 1.478940284762099e-012_r, 1.5e-16_r ) );
		CHECK( fp_equal_tol( y[2], 1.421495882582395e-016_r, 1.4e-20_r ) );
		VoigtH( 1.e-6_r, v, y, 3 );
		CHECK( fp_equal_tol( y[0], 9.960377466254801e-001_r, 1.e-4_r ) );
		CHECK( fp_equal_tol( y[1], 1.478934493028404e-008_r, 1.5e-12_r ) );
		CHECK( fp_equal_tol( y[2], 1.421495882582395e-012_r, 1.4e-16_r ) );
		VoigtH( 1.e-2_r, v, y, 3 );
		CHECK( fp_equal_tol( y[0], 9.849424862549039e-001_r, 1.e-4_r ) );
		CHECK( fp_equal_tol( y[1], 1.478930389133934e-004_r, 1.5e-8_r ) );
		CHECK( fp_equal_tol( y[2], 1.421495882224242e-008_r, 1.4e-12_r ) );
		VoigtH( 1.e+1_r, v, y, 3 );
		CHECK( fp_equal_tol( y[0], 5.613881832823886e-002_r, 6.e-6_r ) );
		CHECK( fp_equal_tol( y[1], 4.040671157393835e-002_r, 4.e-6_r ) );
		CHECK( fp_equal_tol( y[2], 1.421137820009847e-005_r, 1.4e-9_r ) );
		VoigtH( 1.2e+1_r, v, y, 3 );
		CHECK( fp_equal_tol( y[0], 4.685295149211637e-002_r, 5.e-6_r ) );
		CHECK( fp_equal_tol( y[1], 3.684277239564798e-002_r, 4.e-6_r ) );
		CHECK( fp_equal_tol( y[2], 1.705176395541707e-005_r, 1.7e-9_r ) );
		VoigtH( 1.5e+1_r, v, y, 3 );
		CHECK( fp_equal_tol( y[0], 3.752895161491574e-002_r, 4.e-6_r ) );
		CHECK( fp_equal_tol( y[1], 3.194834330452605e-002_r, 3.e-6_r ) );
		CHECK( fp_equal_tol( y[2], 2.131035743074598e-005_r, 2.e-9_r ) );
		VoigtH( 2.e+2_r, v, y, 3 );
		CHECK( fp_equal_tol( y[0], 2.820912377324508e-003_r, 3.e-7_r ) );
		CHECK( fp_equal_tol( y[1], 2.818116555672206e-003_r, 3.e-7_r ) );
		CHECK( fp_equal_tol( y[2], 2.582702147491469e-004_r, 3.e-8_r ) );
		VoigtH( 1.e+5_r, v, y, 3 );
		CHECK( fp_equal_tol( y[0], 5.641895835193228e-006_r, 6.e-10_r ) );
		CHECK( fp_equal_tol( y[1], 5.641895812802746e-006_r, 6.e-10_r ) );
		CHECK( fp_equal_tol( y[2], 5.641671917237128e-006_r, 6.e-10_r ) );
		v[0] = 1.e+0_r;
		VoigtH( 1.e-20_r, v, y, 1 );
		CHECK( fp_equal_tol( y[0], 3.678794411714423e-001_r, 4.e-5_r ) );
		v[0] = 5.5e+0_r;
		VoigtH( 1.e-14_r, v, y, 1 );
		CHECK( fp_equal_tol( y[0], 7.307386729528773e-014_r, 7.e-18_r ) );
		v[0] = 3.9e+4_r;
		VoigtH( 1.e+0_r, v, y, 1 );
		CHECK( fp_equal_tol( y[0], 3.709333226385423e-010_r, 4.e-14_r ) );
		v[0] = 1.e+0_r;
		VoigtH( 2.8e+4_r, v, y, 1 );
		CHECK( fp_equal_tol( y[0], 2.014962794529686e-005_r, 2.e-9_r ) );
	}

	TEST(TestVoigtU)
	{
		// check that the Voigt profile returned by VoigtU() is properly normalized
		const int NP = 200;
		realnum v[NP], a, y[NP];
		// for a > 0.1, VoigtU() calls humlik(), for smaller values it
		// calls FastVoigtH() and divides by sqrt(pi).
		//
		// humlik() is set up for a relative precision of 1e-4 over its
		// entire range, but looks to be more precise in practice (at
		// least at v=0)
		// FastVoigtH() is set up for a rel. precision of 2.5e-3 over its
		// entire range, but should be better than 1e-4 for a < 0.0235
		a = 0.0002_r;
		for( int i=0; i < 9; ++i )
		{
			// test both humlik() and FastVoigtH()
			for( int i=0; i < NP; ++i )
				v[i] = realnum(i)*max(a,1_r)/5_r;
			VoigtU( a, v, y, NP );
			realnum integral = 0_r;
			// We need the integral from -infinity to +infinity, which is simply
			// two times the integral from 0 to +infinity. Hence we omit the
			// division by 2 in the trapezoidal rule
			for( int i=1; i < NP; ++i )
				integral += (y[i]+y[i-1])*(v[i]-v[i-1]);
			// add in the integral over the Lorentz wings assuming U(a,v) = c/v^2
			integral += 2_r*v[NP-1]*y[NP-1];
			// VoigtU() calculates U(a,v), so integral should be 1
			CHECK( fp_equal_tol( integral, 1_r, 1.e-4_r ) );
			// also check the central value...
			CHECK( fp_equal_tol( y[0], realnum(VoigtU0(a)), 1.e-4_r ) );
			a *= 10_r;
		}
	}

	TEST(TestMD5string)
	{
		string test;
		// md5sum of an empty file...
		CHECK( MD5string( test ) == "d41d8cd98f00b204e9800998ecf8427e" );
		CHECK( MD5string( test ).length() == NMD5 );
		// check if padding is done correctly
		// an extra block of padding needs to be added when length%64 == 56
		test = "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 55 );
		CHECK( MD5string( test ) == "426ec4ac35ad38d125f6efb39da03098" );
		test = "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 56 );
		CHECK( MD5string( test ) == "d03607b2c89adc0c4abf5a0f1d9e40c9" );
		test = "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 57 );
		CHECK( MD5string( test ) == "bac1b47748411cb6eee0cae3befb8377" );
		string test64 = "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		test = test64 + "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 64+55 );
		CHECK( MD5string( test ) == "10d49aad1fc69976376fbe7c8c5ed118" );
		test = test64 + "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 64+56 );
		CHECK( MD5string( test ) == "61ec7da14576f3b585038c6d72cd5bd5" );
		test = test64 + "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh";
		CHECK( test.length() == 64+57 );
		CHECK( MD5string( test ) == "f17a0475a26d0930e2a35bb320c10e0d" );
		// check that leading zeros are printed correctly
		test = "ghijklmn";
		CHECK( MD5string( test ) == "0256b9cea63bc1f97b8c5aea92c24a98" );
	}

	TEST(TestTrieNode)
	{
		trieNode* root = new trieNode;
		insertToken(root, "NO");
		insertToken(root, "NOINIT");
		insertToken(root, "NONSENSE");
		insertToken(root, "SENSE");
		CHECK( findUniqueLen(root, "NO") == 2 );
		CHECK( findUniqueLen(root, "NOINIT") == 3 );
		CHECK( findUniqueLen(root, "NONSENSE") == 3 );
		CHECK( findUniqueLen(root, "SENSE") == 1 );
		// only ASCII characters are allowed
		string nonASCII(1,(char)128);
		CHECK_THROW( insertToken(root,nonASCII), bad_assert );
		delete root;
	}

	TEST(TestLevenshteinDistance)
	{
		// test all possible combinations of up to two instances of insertion, deletion, or substitution
		CHECK( LevenshteinDistance("test", "test") == 0 );
		CHECK( LevenshteinDistance("test", "test1") == 1 );
		CHECK( LevenshteinDistance("test", "t1est2") == 2 );
		CHECK( LevenshteinDistance("test", "tst") == 1 );
		CHECK( LevenshteinDistance("test", "ts") == 2 );
		CHECK( LevenshteinDistance("test", "tfst") == 1 );
		CHECK( LevenshteinDistance("test", "tasd") == 2 );
		CHECK( LevenshteinDistance("test", "tstr") == 2 );
		CHECK( LevenshteinDistance("test", "tsd") == 2 );
		CHECK( LevenshteinDistance("test", "teysd") == 2 );
	}

}
