/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*GrainMakeDiffuse main routine for generating the grain diffuse emission, called by RT_diffuse */
#include "cddefines.h"
#include "rfield.h"
#include "phycon.h"
#include "dense.h"
#include "hmi.h"
#include "thermal.h"
#include "trace.h"
#include "iterations.h"
#include "vectorize.h"
#include "grainvar.h"
#include "grains.h"

inline double no_atoms(size_t nd)
{
	return gv.bin[nd].AvVol*gv.bin[nd].dustp[0]/ATOMIC_MASS_UNIT/gv.bin[nd].atomWeight;
}

/* NB NB -- the sequence below has been carefully chosen and should NEVER be
 *          altered unless you really know what you are doing !! */
/* >>chng 03 jan 16, introduced QH_THIGH_FAIL and started using splint_safe and spldrv_safe
 *                   throughout the code; this solves the parispn_a031_sil.in bug, PvH */
/* >>chng 03 jan 16, rescheduled QH_STEP_FAIL as non-fatal; it can be triggered at very low temps, PvH */
typedef enum {
	/* the following are OK */
	/* 0        1               2              3    */
	QH_OK, QH_ANALYTIC, QH_ANALYTIC_RELAX, QH_DELTA, 
	/* the following are mild errors we already recovered from */
	/*     4              5              6             7        */
	QH_NEGRATE_FAIL, QH_LOOP_FAIL, QH_ARRAY_FAIL, QH_THIGH_FAIL,
	/* any of these errors will prompt qheat to do another iteration */
	/*  8          9             10             11       */
	QH_RETRY, QH_CONV_FAIL, QH_BOUND_FAIL, QH_DELTA_FAIL,
	/* any error larger than QH_NO_REBIN will cause GetProbDistr_LowLimit to return
	 * before even attempting to rebin the results; we may be able to recover though... */
	/*   12          13            14            15       */
	QH_NO_REBIN, QH_LOW_FAIL, QH_HIGH_FAIL, QH_STEP_FAIL,
	/* any case larger than QH_FATAL is truly pathological; there is no hope of recovery */
	/* 16          17           18             19      */
	QH_FATAL, QH_WIDE_FAIL, QH_NBIN_FAIL, QH_REBIN_FAIL
} QH_Code;

/*================================================================================*/
/* definitions relating to quantum heating routines */

/* this is the minimum number of bins for quantum heating calculation to be valid */
static const long NQMIN = 10L;

/* this is the lowest value for dPdlnT that should be included in the modeling */
static const double PROB_CUTOFF_LO = 1.e-15;
static const double PROB_CUTOFF_HI = 1.e-20;
static const double SAFETY = 1.e+8;

/* during the first NSTARTUP steps, use special algorithm to calculate stepsize */
static const long NSTARTUP = 5L;

/* if the average number of multiple events is above this number
 * don't try full quantum heating treatment. */
static const double MAX_EVENTS = 150.;

/* maximum number of tries for quantum heating routine */
/* >>chng 02 jan 30, changed LOOP_MAX from 10 -> 20, PvH */
static const long LOOP_MAX = 20L;

/* if all else fails, divide temp estimate by this factor */
static const double DEF_FAC = 3.;

/* total probability for all bins should not deviate more than this from 1. */
static const double PROB_TOL = 0.02;

/* after NQTEST steps, make an estimate if prob distr will converge in NQGRID steps */
/* >>chng 02 jan 30, change NQTEST from 1000 -> 500, PvH */
static const long NQTEST = 500L;

/* if the ratio fwhm/Umax is lower than this number
 * don't try full quantum heating treatment. */
static const double FWHM_RATIO = 0.1;
/* if the ratio fwhm/Umax is lower than this number
 * don't even try analytic approximation, use delta function instead */
static const double FWHM_RATIO2 = 0.007;

/* maximum number of steps for integrating quantum heating probability distribution */
static const long MAX_LOOP = 2*NQGRID;

/* this is the tolerance used while integrating the temperature distribution of the grains */
static const double QHEAT_TOL = 5.e-3;

/* maximum number of tries before we declare that probability distribution simply won't fit */
static const long WIDE_FAIL_MAX = 3;

/* multipliers for PROB_TOL used in GetProbDistr_HighLimit */
static const double STRICT = 1.;
static const double RELAX = 15.;

/* when rebinning quantum heating results, make ln(qtemp[i]/qtemp[i-1]) == QT_RATIO */
/* this constant determines the accuracy with which the Wien tail of the grain emission is
 * calculated; if x = h*nu/k*T_gr, d = QT_RATIO-1., and p(T) = p0 * T_gr^C, then the
 * relative accuracy of the flux in the Wien tail is
 *   rel. acc. = d^2/24*fabs(C^2 + (2x-1)*C + x^2 -2*x)
 * A typical value of x would be x = 12 (this corresponds to Bnu ~ 1e-2*Bnu(peak)), and C
 * varies around -4 near the highest temperatures, so QT_RATIO = 1.074 should converge the
 * spectrum to 1% at x = 12 and better closer to the peak of the spectrum. Since C can
 * vary, the actual error should be between 0 and 2% at x=12. */
static const double QT_RATIO = 1.074;


/*================================================================================*/
/* global variables */

/* these data define the enthalpy function for silicates
 * derived from:
 * >>refer	grain	physics	Guhathakurta & Draine, 1989, ApJ, 345, 230
 * coefficients converted to rydberg energy units, per unit atom
 * assuming a density of 3.3 g/cm^3 and pure MgSiFeO4.
 * this is not right, but the result is correct because number
 * of atoms will be calculated using the same assumption. */

/* this is the specific density of silicate in g/cm^3 */
static const double DEN_SIL = 3.30;

/* these are the mean molecular weights per atom for MgSiFeO4 and SiO2 in amu */
static const double MW_SIL = 24.6051;
/*static const double MW_SIO2 = 20.0283;*/

static const double tlim[5]={0.,50.,150.,500.,DBL_MAX};
static const double ppower[4]={2.00,1.30,0.68,0.00};
static const double cval[4]={
	1.40e3/DEN_SIL*ATOMIC_MASS_UNIT*MW_SIL/EN1RYD,
	2.20e4/DEN_SIL*ATOMIC_MASS_UNIT*MW_SIL/EN1RYD,
	4.80e5/DEN_SIL*ATOMIC_MASS_UNIT*MW_SIL/EN1RYD,
	3.41e7/DEN_SIL*ATOMIC_MASS_UNIT*MW_SIL/EN1RYD};


/* initialize phiTilde */
STATIC void qheat_init(size_t,/*@out@*/vector<double>&,/*@out@*/double*);
/* worker routine, this implements the algorithm of Guhathakurtha & Draine */
STATIC void GetProbDistr_LowLimit(size_t,double,double,double,/*@in@*/const vector<double>&,
				  /*@in@*/const vector<double>&,/*@out@*/vector<double>&,
				  /*@out@*/vector<double>&,/*@out@*/vector<double>&,
				  /*@out@*/long*,/*@out@*/double*,long*,/*@out@*/QH_Code*);
/* try two consecutive integration steps using stepsize "step/2." (yielding p[k]),
 * and also one double integration step using stepsize "step" (yielding p2k). */
STATIC double TryDoubleStep(vector<double>&,vector<double>&,vector<double>&,vector<double>&,
			    vector<double>&,const vector<double>&,const vector<double>&,double,
			    /*@out@*/double*,double,long,size_t,/*@out@*/bool*);
/* calculate logarithmic integral from (x1,y1) to (x2,y2) */
STATIC double log_integral(double,double,double,double,double,double,double,double);
/* scan the probability distribution for valid range */
STATIC void ScanProbDistr(const vector<double>&,const vector<double>&,long,double,long,double,
			  /*@out@*/long*,/*@out@*/long*,/*@out@*/long*,long*,QH_Code*);
/* rebin the quantum heating results to speed up RT_diffuse */
STATIC long RebinQHeatResults(size_t,long,long,vector<double>&,vector<double>&,vector<double>&,
			      vector<double>&,vector<double>&,vector<double>&,vector<double>&,QH_Code*);
/* calculate approximate probability distribution in high intensity limit */
STATIC void GetProbDistr_HighLimit(long,double,double,double,/*@out@*/vector<double>&,/*@out@*/vector<double>&,
				   /*@out@*/vector<double>&,/*@out@*/double*,
				   /*@out@*/long*,/*@out@*/double*,/*@out@*/QH_Code*);
/* derivative of the enthalpy function dU/dT */
STATIC double uderiv(double,size_t);
/* enthalpy function */
STATIC double ufunct(double,size_t,/*@out@*/bool*);
/* inverse enthalpy function */
STATIC double inv_ufunct(double,size_t,/*@out@*/bool*);
/* helper function for calculating enthalpy, uses Debye theory */
STATIC double DebyeDeriv(double,long);

/* >>chng 01 oct 29, introduced gv.bin[nd].cnv_H_pGR, cnv_GR_pH, etc. PvH */


#ifndef __AVX__

/* this assures 1e-4 relative precision in the evaluation of exp(x)-1 below */
static const realnum LIM1 = 2.e-4f;
static const realnum LIM2 = sqrtf(6.e-4f);
static const realnum LIM3 = cbrtf(24.e-4f);

// fast version of expm1() with reduced precision, it implicitly assumes x >= 0!
inline realnum fast_expm1(realnum x)
{
	if( x < LIM1 )
		return x;
	else if( x < LIM2 )
		return (x/2.f + 1.f)*x;
	else if( x < LIM3 )
		return ((x/6.f + 0.5f)*x +1.f)*x;
	else
		return expf(x) - 1.f;
}

#endif

STATIC long GrainMakeDiffuseSingle(double Tgrain, double fracpop, avx_ptr<realnum>& flux, long nflux)
{
	avx_ptr<realnum> arg(nflux), val(nflux);
	const realnum x = log(0.999f*FLT_MAX);

	double hokT = TE1RYD/Tgrain;
	// NB NB -- the loops in this routine consume lots of CPU time, they should be well optimized!!
	for( long i=0; i < nflux; i++ )
	{
		arg[i] = hokT*rfield.anu(i);
		if( arg[i] > x )
		{
			nflux = i;
			break;
		}
	}
#ifdef __AVX__
	vexpm1( arg.ptr0(), val.ptr0(), 0, nflux );
	for( long i=0; i < nflux; i++ )
		flux[i] += fracpop/val[i];
#else
	for( long i=0; i < nflux; i++ )
		flux[i] += fracpop/fast_expm1(arg[i]);
#endif
	return nflux;
}

/* main routine for generating the grain diffuse emission, called by RT_diffuse */
void GrainMakeDiffuse()
{
	DEBUG_ENTRY( "GrainMakeDiffuse()" );

	const double factor = 2.*PI4*pow2(FR1RYD/SPEEDLIGHT)*FR1RYD;

	/* save grain emission per unit volume */
	for( long i=0; i < rfield.nflux; i++ )
	{
		/* used in RT_diffuse to derive total emission */
		gv.GrainEmission[i] = 0.;
		gv.SilicateEmission[i] = 0.;
		gv.GraphiteEmission[i] = 0.;
	}

	vector<double> qtemp(NQGRID);
	vector<double> qprob(NQGRID);
	avx_ptr<realnum> flux(rfield.nflux);

	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		/* this local copy is necessary to keep lint happy */
		bool lgLocalQHeat = gv.bin[nd].lgQHeat;
		/* >>chng 04 nov 09, do not evaluate quantum heating if abundance is negligible, PvH
		 * this prevents PAH's deep inside molecular regions from failing if GrnVryDepth is used */
		/* >>chng 04 dec 31, introduced separate thresholds near I-front and in molecular region, PvH */
		realnum threshold = ( dense.xIonDense[ipHYDROGEN][0]+dense.xIonDense[ipHYDROGEN][1] > hmi.H2_total ) ?
			gv.dstAbundThresholdNear : gv.dstAbundThresholdFar;
		long qnbin=-200;

		if( lgLocalQHeat && gv.bin[nd].dstAbund >= threshold )
		{
			qheat(qtemp,qprob,&qnbin,nd);

			if( gv.bin[nd].lgUseQHeat )
			{
				ASSERT( qnbin > 0 );
			}
		}
		else
		{
			/* >> chng 04 dec 31, repaired bug lgUseQHeat not set when abundance below threshold, PvH */
			gv.bin[nd].lgUseQHeat = false;
		}

		long loopmax = 0;
		memset( flux.data(), 0, size_t(rfield.nflux*sizeof(flux[0])) );

		if( lgLocalQHeat && gv.bin[nd].lgUseQHeat )
		{
			for( long j=0; j < qnbin; j++ )
			{
				long maxi = GrainMakeDiffuseSingle(qtemp[j], qprob[j], flux, rfield.nflux);
				loopmax = max(loopmax, maxi);
			}
		}
		else
		{
			loopmax = GrainMakeDiffuseSingle(gv.bin[nd].tedust, 1., flux, rfield.nflux);
		}

		realnum* gt_ptr;
		/* unit emission for each different grain type */
		strg_type scase = gv.which_strg[gv.bin[nd].matType];
		switch( scase )
		{
		case STRG_SIL:
			gt_ptr = get_ptr(gv.SilicateEmission);
			break;
		case STRG_CAR:
			gt_ptr = get_ptr(gv.GraphiteEmission);
			break;
		default:
			TotalInsanity();
		}

		double fac = factor*gv.bin[nd].cnv_H_pCM3;
		// use two separate loops so that they can be vectorized
		for( long i=0; i < loopmax; i++ )
			flux[i] *= realnum(fac*gv.bin[nd].dstab1[i]*rfield.anu2(i)*rfield.widflx(i));
		for( long i=0; i < loopmax; i++ )
		{
			/* remember local emission -- these are zeroed out on each zone 
			 * above, and now incremented so is unit emission from this zone */
			gv.GrainEmission[i] += flux[i];
			gt_ptr[i] += flux[i];
		}
	}

#	ifndef NDEBUG
	/*********************************************************************************
	 *
	 * Following are three checks on energy and charge conservation by the grain code.
	 * Their primary function is as an internal consistency check, so that coding
	 * errors get caught as early as possible. This is why the program terminates as
	 * soon as any one of the checks fails.
	 *
	 * NB NB - despite appearances, these checks do NOT guarantee overall energy
	 *         conservation in the Cloudy model to the asserted tolerance, see note 1B !
	 *
	 * Note 1: there are two sources for energy imbalance in the grain code (see A & B).
	 *   A: Interpolation in dstems. The code calculates how much energy the grains
	 *      emit in thermal radiation (gv.bin[nd].GrainHeat), and converts that into
	 *      an (average) grain temperature by reverse interpolation in dstems. If
	 *      quantum heating is not used, that temperature is used directly to generate
	 *      the local diffuse emission. Hence the finite resolution of the dstems grid
	 *      can lead to small errors in flux. This is tested in Check 1. The maximum
	 *      error of interpolating in dstems scales with NDEMS^-3. The same problem
	 *      can also occur when quantum heating is used, however, the fact that many
	 *      bins are used will probably lead to a cancellation effect of the errors.
	 *   B: RT_OTS_Update gets called AFTER grain() has completed, so the grain heating
	 *      was calculated using a different version of the OTS fields than the one
	 *      that gets fed into the RT routines (where the actual attenuation of the
	 *      radiation fields by the grain opacity is done). This can lead to an energy
	 *      imbalance, depending on how accurate the convergence of the OTS fields is.
	 *      This is outside the control of the grain code and is therefore NOT checked.
	 *      Rather, the grain code remembers the contribution from the old OTS fields
	 *      (through gv.bin[nd].BolFlux) and uses that in Check 3. In most models the
	 *      difference will be well below 0.1%, but in AGN type models where OTS continua
	 *      are important, the energy imbalance can be of the order of 0.5% of the grain
	 *      heating (status nov 2001). On 04 jan 25 the initialization of phiTilde has
	 *      been moved to qheat, implying that phiTilde now uses the updated version of
	 *      the OTS fields. The total amount of radiated energy however is still based
	 *      on gv.bin[nd].GrainHeat which uses the old version of the OTS fields.
	 *   C: Energy conservation for collisional processes is guaranteed by adding in
	 *      (very small) correction terms. These corrections are needed to cover up
	 *      small imperfection in the theory, and cannot be avoided without making the
	 *      already very complex theory even more complex.
	 *   D: Photo-electric heating and collisional cooling can have an important effect
	 *      on the total heating balance of the gas. Both processes depend strongly on
	 *      the grain charge, so assuring proper charge balance is important as well.
	 *      This is tested in Check 2.
	 *
	 * Note 2: for quantum heating it is important to resolve the Maxwell distribution
	 *   of the electrons and ions far enough into the tail that the total amount of
	 *   energy contained in the remaining tail is negligible. If this is not the case,
	 *   the assert at the beginning of the qheat() routine will fail. This is because
	 *   the code filling up the phiTilde array in GrainCollHeating() assumes a value for
	 *   the average particle energy based on a Maxwell distribution going to infinity.
	 *   If the maximum energy used is too low, the assumed average energy would be
	 *   incorrect.
	 *
	 *********************************************************************************/

	bool lgNoTdustFailures = true;
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		if( !gv.bin[nd].lgTdustConverged )
		{
			lgNoTdustFailures = false;
			break;
		}
	}

	/* CHECK 1: does the grain thermal emission conserve energy ? */
	double BolFlux = 0.;
	for( long i=0; i < rfield.nflux; i++ )
	{
		BolFlux += gv.GrainEmission[i]*rfield.anu(i)*EN1RYD;
	}
	double Comparison1 = 0.;
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		if( gv.bin[nd].tedust < gv.bin[nd].Tsublimat )
			Comparison1 += CONSERV_TOL*gv.bin[nd].GrainHeat;
		else
			/* for high temperatures the interpolation in dstems
			 * is less accurate, so we have to be more lenient */
			Comparison1 += 10.*CONSERV_TOL*gv.bin[nd].GrainHeat;
	}

	/* >>chng 04 mar 11, add constant grain temperature to pass assert */
	/* >>chng 04 jun 01, deleted test for constant grain temperature, PvH */
	ASSERT( fabs(BolFlux-gv.GrainHeatSum) < Comparison1 );

	/* CHECK 2: assert charging balance */
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		double ave = 0.5*(gv.bin[nd].RateDn+gv.bin[nd].RateUp);
		ASSERT( fabs(gv.bin[nd].RateDn-gv.bin[nd].RateUp) < CONSERV_TOL*ave );
	}

	if( lgNoTdustFailures && gv.lgDHetOn && gv.lgDColOn && thermal.ConstGrainTemp == 0. )
	{
		/* CHECK 3: calculate the total energy donated to grains, must be balanced by
		 * the energy emitted in thermal radiation plus various forms of gas heating */
		Comparison1 = 0.;
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			Comparison1 += gv.bin[nd].BolFlux;
		}
		/* add in collisional heating of grains by plasma (if positive) */
		Comparison1 += MAX2(gv.GasCoolColl,0.);
		/* add in net amount of chemical energy donated by recombining ions and molecule formation */
		Comparison1 += gv.GrainHeatChem;

		/*              thermal emis        PE effect          gas heating by coll    thermionic emis */
		double Comparison2 = gv.GrainHeatSum+thermal.heating(0,13)+thermal.heating(0,14)+thermal.heating(0,25);

		/* >>chng 06 jun 02, add test on gv.GrainHeatScaleFactor so that assert not thrown
		 * when set grain heat command is used */
		ASSERT( gv.GrainHeatScaleFactor != 1.f || gv.lgBakesPAH_heat ||
			fabs(Comparison1-Comparison2)/Comparison2 < CONSERV_TOL );
	}
#	endif
	return;
}


/****************************************************************************
 *
 * qheat: driver routine for non-equilibrium heating of grains
 *
 * This routine calls GetProbDistr_LowLimit, GetProbDistr_HighLimit
 * (which do the actual non-equilibrium calculations), and does the
 * subsequent quality control.
 *
 * Written by Peter van Hoof (UK, CITA, QUB).
 *
 ****************************************************************************/

/* this is the new version of the quantum heating code, used starting Cloudy 96 beta 3 */

void qheat(/*@out@*/ vector<double>& qtemp, /* qtemp[NQGRID] */
	   /*@out@*/ vector<double>& qprob, /* qprob[NQGRID] */
	   /*@out@*/ long int *qnbin,
	   size_t nd)
{
	bool lgBoundErr,
	  lgDelta,
	  lgNegRate,
	  lgOK,
	  lgTried;
	long int i,
	  nWideFail;
	QH_Code ErrorCode,
	  ErrorCode2,
	  ErrorStart;
	double c0,
	  c1,
	  c2,
	  check,
	  DefFac,
	  deriv,
	  fwhm,
	  FwhmRatio,
	  integral,
	  minBracket,
	  maxBracket,
	  new_tmin,
	  NumEvents,
	  rel_tol,
	  Tmax,
	  tol,
	  Umax,
	  xx,
	  y;

	DEBUG_ENTRY( "qheat()" );

	/* sanity checks */
	ASSERT( gv.bin[nd].lgQHeat );
	ASSERT( nd < gv.bin.size() );

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "\n >>>>qheat called for grain %s\n", gv.bin[nd].chDstLab );
	}

	/* >>chng 01 aug 22, allocate space */
	/* phiTilde is continuum corrected for photo-electric effect, in events/H/s/cell, default depl */
	vector<double> phiTilde(rfield.nflux_with_check);
	vector<double> Phi(rfield.nflux_with_check);
	vector<double> PhiDrv(rfield.nflux_with_check);
	vector<double> dPdlnT(NQGRID);

	qheat_init( nd, phiTilde, &check );

	check += gv.bin[nd].GrainHeatColl-gv.bin[nd].GrainCoolTherm;

	xx = integral = 0.;
	c0 = c1 = c2 = 0.;
	lgNegRate = false;
	tol = 1.;

	/* >>chng 01 nov 29, rfield.nflux -> gv.qnflux, PvH */
	/* >>chng 03 jan 26, gv.qnflux -> gv.bin[nd].qnflux, PvH */
	for( i=gv.bin[nd].qnflux-1; i >= 0; i-- )
	{
		/* >>chng 97 jul 17, to summed continuum */
		/* >>chng 00 mar 30, to phiTilde, to keep track of photo-electric effect and collisions, by PvH */
		/* >>chng 01 oct 10, use trapezoidal rule for integrating Phi, reverse direction of integration. */
		/* >>chng 01 oct 30, change normalization of Phi, PhiDrv from <unit>/cm^3 -> <unit>/grain, PvH */
		/* phiTilde has units events/H/s, PhiDrv[i] has units events/grain/s/Ryd */
		/* there are minus signs here because we are integrating from infinity downwards */
		y = -phiTilde[i]*gv.bin[nd].cnv_H_pGR;
		PhiDrv[i] = y/rfield.widflx(i);
		xx -= y;
		/* Phi[i] is integral from exactly rfield.anumin(i) to infinity to second-order precision, PvH */
		/* Phi[i] has units events/grain/s */
		Phi[i] = xx;

#		ifndef NDEBUG
		/* trapezoidal rule is not needed for integral, this is also second-order correct */
		integral += phiTilde[i]*gv.bin[nd].cnv_H_pCM3*rfield.anu(i)*EN1RYD;
#		endif

		/* c<n> has units Ryd^(n+1)/grain/s */
		c0 += Phi[i]*rfield.widflx(i);
		c1 += Phi[i]*rfield.anu(i)*rfield.widflx(i);
		c2 += Phi[i]*rfield.anu2(i)*rfield.widflx(i);

		lgNegRate = lgNegRate || ( phiTilde[i] < 0. );
	}

	/* sanity check */
	ASSERT( fabs(check-integral)/check <= CONSERV_TOL );

#	if 0
	{
		char fnam[50];
		FILE *file;

		sprintf(fnam,"Lambda_%2.2ld.asc",nd);
		file = open_data(fnam,"w");
		for( i=0; i < NDEMS; ++i )
			fprintf(file,"%e %e %e\n",
				exp(gv.dsttmp[i]),
				ufunct(exp(gv.dsttmp[i]),nd,&lgBoundErr),
				exp(gv.bin[nd].dstems[i])*gv.bin[nd].cnv_H_pGR/EN1RYD);
		fclose(file);

		sprintf(fnam,"Phi_%2.2ld.asc",nd);
		file = open_data(fnam,"w");
		for( i=0; i < gv.bin[nd].qnflux; ++i )
			fprintf(file,"%e %e\n", rfield.anu(i),Phi[i]);
		fclose(file);
	}
#	endif

	/* Tmax is where p(U) will peak (at least in high intensity limit) */
	Tmax = gv.bin[nd].tedust;
	/* grain enthalpy at peak of p(U), in Ryd */
	Umax = ufunct(Tmax,nd,&lgBoundErr);
	ASSERT( !lgBoundErr ); /* this should never happen */
	/* y is dln(Lambda)/dlnT at Tmax */
	spldrv_safe(gv.dsttmp,gv.bin[nd].dstems,gv.bin[nd].dstslp2,NDEMS,log(Tmax),&y,&lgBoundErr);
	ASSERT( !lgBoundErr ); /* this should never happen */
	/* deriv is dLambda/dU at Umax, in 1/grain/s */
	deriv = y*c0/(uderiv(Tmax,nd)*Tmax);
	/* estimated FWHM of probability distribution, in Ryd */
	fwhm = sqrt(8.*LN_TWO*c1/deriv);

	NumEvents = pow2(fwhm)*c0/(4.*LN_TWO*c2);
	FwhmRatio = fwhm/Umax;

	/* >>chng 01 nov 15, change ( NumEvents > MAX_EVENTS2 ) --> ( FwhmRatio < FWHM_RATIO2 ), PvH */
	lgDelta = ( FwhmRatio < FWHM_RATIO2 );
	/* high intensity case is always OK since we will use equilibrium treatment */
	lgOK = lgDelta;

	ErrorStart = QH_OK;

	if( lgDelta ) 
	{
		/* in this case we ignore negative rates, equilibrium treatment is good enough */
		lgNegRate = false;
		ErrorStart = MAX2(ErrorStart,QH_DELTA);
	}

	if( lgNegRate )
		ErrorStart = MAX2(ErrorStart,QH_NEGRATE_FAIL);

	ErrorCode = ErrorStart;

	if( trace.lgTrace && trace.lgDustBug )
	{
		double Rate2 = 0.;
		for( int nz=0; nz < gv.bin[nd].nChrg; nz++ )
			Rate2 += gv.bin[nd].chrg(nz).FracPop*gv.bin[nd].chrg(nz).HeatingRate2;

		fprintf( ioQQQ, "   grain heating: %.4e, integral %.4e, total rate %.4e lgNegRate %c\n",
			 gv.bin[nd].GrainHeat,integral,Phi[0],TorF(lgNegRate));
		fprintf( ioQQQ, "   av grain temp %.4e av grain enthalpy (Ryd) %.4e\n",
			 gv.bin[nd].tedust,Umax);
		fprintf( ioQQQ, "   fwhm^2/(4ln2*c2/c0): %.4e fwhm (Ryd) %.4e fwhm/Umax %.4e\n",
			 NumEvents,fwhm,FwhmRatio );
		fprintf( ioQQQ, "   HeatingRate1 %.4e HeatingRate2 %.4e lgQHTooWide %c\n",
			 gv.bin[nd].HeatingRate1*gv.bin[nd].cnv_H_pCM3, Rate2*gv.bin[nd].cnv_H_pCM3,
			 TorF(gv.bin[nd].lgQHTooWide) );
	}

	/* these two variables will bracket qtmin, they should only be needed during the initial search phase */
	minBracket = GRAIN_TMIN;
	maxBracket = gv.bin[nd].tedust;

	/* >>chng 02 jan 30, introduced lgTried to avoid running GetProbDistr_HighLimit over and over..., PvH */
	lgTried = false;
	/* >>chng 02 aug 06, introduced QH_WIDE_FAIL and nWideFail, PvH */
	nWideFail = 0;
	/* >>chng 03 jan 27, introduced DefFac to increase factor for repeated LOW_FAIL's, PvH */
	DefFac = DEF_FAC;
	/* >>chng 04 nov 10, introduced rel_tol to increase precision in case of repeated CONV_FAIL's, PvH */
	rel_tol = 1.;

	/* if average number of multiple photon events is too high, lgOK is set to true */
	/* >>chng 02 aug 12, added gv.bin[nd].lgQHTooWide to prevent unnecessary looping here.
	 * In general the number of integration steps needed to integrate the probability distribution
	 * will increase monotonically with depth into the cloud. Hence, once the distribution becomes
	 * too wide to fit into NQGRID steps (which only happens for extremely cold grains in deeply
	 * shielded conditions) there is no hope of ever converging GetProbDistr_LowLimit and the code
	 * will waste a lot of CPU time establishing this for every zone again. So once the distribution
	 * becomes too wide we immediately skip to the analytic approximation to save time, PvH */
	for( i=0; i < LOOP_MAX && !lgOK && !gv.bin[nd].lgQHTooWide; i++ )
	{
		if( gv.bin[nd].qtmin >= gv.bin[nd].tedust )
		{
			/* >>chng 02 jul 31, was gv.bin[nd].qtmin = 0.7*gv.bin[nd].tedust */
			/* >>chng 03 nov 10, changed Umax/exp(+... to Umax*exp(-... to avoid overflow, PvH */
			double Ulo = Umax*exp(-sqrt(-log(PROB_CUTOFF_LO)/(4.*LN_TWO))*fwhm/Umax);
			double MinEnth = exp(gv.bin[nd].DustEnth[0]);
			Ulo = MAX2(Ulo,MinEnth);
			gv.bin[nd].qtmin = inv_ufunct(Ulo,nd,&lgBoundErr);
			ASSERT( !lgBoundErr ); /* this should never happen */
			/* >>chng 02 jul 30, added this test; prevents problems with ASSERT below, PvH */
			if( gv.bin[nd].qtmin <= minBracket || gv.bin[nd].qtmin >= maxBracket )
				gv.bin[nd].qtmin = sqrt(minBracket*maxBracket);
		}
		gv.bin[nd].qtmin = MAX2(gv.bin[nd].qtmin,GRAIN_TMIN);

		ASSERT( minBracket <= gv.bin[nd].qtmin && gv.bin[nd].qtmin <= maxBracket );

		ErrorCode = ErrorStart;

		/* >>chng 01 nov 15, add ( FwhmRatio >= FWHM_RATIO ), PvH */
		if( FwhmRatio >= FWHM_RATIO && NumEvents <= MAX_EVENTS )
		{
			GetProbDistr_LowLimit(nd,rel_tol,Umax,fwhm,Phi,PhiDrv,qtemp,qprob,dPdlnT,qnbin,
					      &new_tmin,&nWideFail,&ErrorCode);

			/* >>chng 02 jan 07, various changes to improve convergence for very cold grains, PvH */
			if( ErrorCode == QH_DELTA_FAIL && fwhm < Umax && !lgTried )
			{
				double dummy;

				/* this situation can mean two things: either the photon rate is so high that
				 * the code needs too many steps to integrate the probability distribution,
				 * or alternatively, tmin is far too low and the code needs too many steps
				 * to overcome the rising side of the probability distribution.
				 * So we call GetProbDistr_HighLimit first to determine if the former is the
				 * case; if that fails then the latter must be true and we reset QH_DELTA_FAIL */
				ErrorCode = MAX2(ErrorStart,QH_ANALYTIC);
				/* use dummy to avoid losing estimate for new_tmin from GetProbDistr_LowLimit */
				/* >>chng 02 aug 06, introduced STRICT and RELAX, PvH */
				GetProbDistr_HighLimit(nd,STRICT,Umax,fwhm,qtemp,qprob,dPdlnT,&tol,qnbin,&dummy,
						       &ErrorCode);

				if( ErrorCode >= QH_RETRY )
				{
					ErrorCode = QH_DELTA_FAIL;
					lgTried = true;
				}
			}

			/* >>chng 02 aug 07 major rewrite of the logic below */
			if( ErrorCode < QH_NO_REBIN )
			{
				if( new_tmin < minBracket || new_tmin > maxBracket )
					++nWideFail;

				if( nWideFail < WIDE_FAIL_MAX )
				{
					if( new_tmin <= minBracket )
						new_tmin = sqrt(gv.bin[nd].qtmin*minBracket);
					if( new_tmin >= maxBracket )
						new_tmin = sqrt(gv.bin[nd].qtmin*maxBracket);
				}
				else
				{
					ErrorCode = MAX2(ErrorCode,QH_WIDE_FAIL);
				}

				if( ErrorCode == QH_CONV_FAIL )
				{
					rel_tol *= 0.9;
				}
			}
			else if( ErrorCode == QH_LOW_FAIL )
			{
				double help1 = gv.bin[nd].qtmin*sqrt(DefFac);
				double help2 = sqrt(gv.bin[nd].qtmin*maxBracket);
				minBracket = gv.bin[nd].qtmin;
				new_tmin = MIN2(help1,help2);
				/* increase factor in case we get repeated LOW_FAIL's */
				DefFac += 1.5;
			}
			else if( ErrorCode == QH_HIGH_FAIL )
			{
				double help = sqrt(gv.bin[nd].qtmin*minBracket);
				maxBracket = gv.bin[nd].qtmin;
				new_tmin = MAX2(gv.bin[nd].qtmin/DEF_FAC,help);
			}
			else
			{
				new_tmin = sqrt(minBracket*maxBracket);
			}
		}
		else
		{
			GetProbDistr_HighLimit(nd,STRICT,Umax,fwhm,qtemp,qprob,dPdlnT,&tol,qnbin,&new_tmin,
					       &ErrorCode);
		}

		gv.bin[nd].qtmin = new_tmin;

		lgOK = ( ErrorCode < QH_RETRY );

		if( ErrorCode >= QH_FATAL )
			break;

		if( ErrorCode != QH_LOW_FAIL )
			DefFac = DEF_FAC;

		if( trace.lgTrace && trace.lgDustBug ) 
		{
			fprintf( ioQQQ, "   GetProbDistr returns code %d\n", ErrorCode );
			if( !lgOK )
			{
				fprintf( ioQQQ, " >>>>qheat trying another iteration, qtmin bracket %.4e %.4e",
					 minBracket,maxBracket );
				fprintf( ioQQQ, " nWideFail %ld\n", nWideFail );
			}
		}
	}

	if( ErrorCode == QH_WIDE_FAIL )
		gv.bin[nd].lgQHTooWide = true;

	/* >>chng 03 jan 17, added test for !lgDelta, PvH */
	/* if( gv.bin[nd].lgQHTooWide ) */
	if( gv.bin[nd].lgQHTooWide && !lgDelta )
		ErrorCode = MAX2(ErrorCode,QH_WIDE_FAIL);

/* 	if( ErrorCode >= QH_RETRY ) */
/* 		printf( "      *** PROBLEM  loop not converged, errorcode %d\n",ErrorCode ); */

	/* The quantum heating code tends to run into trouble when it goes deep into the neutral zone,
	 * especially if the original spectrum was very hard, as is the case in high excitation PNe or AGN.
	 * You then get a bipartition in the spectrum where most of the photons have low energy, while
	 * there still are hard X-ray photons left. The fact that the average energy per photon is low
	 * forces the quantum code to take tiny little steps when integrating the probability distribution,
	 * while the fact that X-ray photons are still present implies that big temperature spikes still
	 * occur and hence the temperature distribution is very wide. Therefore the code needs a zillion
	 * steps to integrate the probability distribution and simply runs out of room. As a last resort
	 * try the analytic approximation with relaxed constraints used below. */
	/* >>chng 02 oct 03, vary Umax and fwhm to force fit with fwhm/Umax remaining constant */
	/* >>chng 03 jan 17, changed test so that last resort always gets tried when lgOK = lgDelta = false, PvH */
	/* if( !lgOK && FwhmRatio >= FWHM_RATIO && NumEvents <= MAX_EVENTS ) */
	if( !lgOK && !lgDelta )
	{
		double Umax2 = Umax*sqrt(tol);
		double fwhm2 = fwhm*sqrt(tol);

		for( i=0; i < LOOP_MAX; ++i )
		{
			double dummy;

			ErrorCode2 = MAX2(ErrorStart,QH_ANALYTIC);
			GetProbDistr_HighLimit(nd,RELAX,Umax2,fwhm2,qtemp,qprob,dPdlnT,&tol,qnbin,&dummy,
					       &ErrorCode2);

			lgOK = ( ErrorCode2 < QH_RETRY );
			if( lgOK )
			{
				gv.bin[nd].qtmin = dummy;
				ErrorCode = ErrorCode2;
				break;
			}
			else
			{
				Umax2 *= sqrt(tol);
				fwhm2 *= sqrt(tol);
			}
		}
	}

	if( nzone == 1 )
		gv.bin[nd].qtmin_zone1 = gv.bin[nd].qtmin;

	gv.bin[nd].lgUseQHeat = ( lgOK && !lgDelta );
	gv.bin[nd].lgEverQHeat = ( gv.bin[nd].lgEverQHeat || gv.bin[nd].lgUseQHeat );

	if( lgOK )
	{
		if( trace.lgTrace && trace.lgDustBug )
			fprintf( ioQQQ, " >>>>qheat converged with code: %d\n", ErrorCode );
	}
	else
	{
		*qnbin = 0;
		++gv.bin[nd].QHeatFailures;
		fprintf( ioQQQ, " PROBLEM  qheat did not converge grain %s in zone %ld, error code = %d\n",
			 gv.bin[nd].chDstLab,nzone,ErrorCode );		
	}

	if( gv.QHSaveFile != NULL && ( iterations.lgLastIt || !gv.lgQHPunLast ) ) 
	{
		fprintf( gv.QHSaveFile, "\nDust Temperature Distribution: grain %s zone %ld\n",
			 gv.bin[nd].chDstLab,nzone );

		fprintf( gv.QHSaveFile, "Equilibrium temperature: %.2f\n", gv.bin[nd].tedust );

		if( gv.bin[nd].lgUseQHeat ) 
		{
			/* >>chng 01 oct 09, remove qprob from output, it depends on step size, PvH */
			fprintf( gv.QHSaveFile, "Number of bins: %ld\n", *qnbin );
			fprintf( gv.QHSaveFile, "  Tgrain      dP/dlnT\n" );
			for( i=0; i < *qnbin; i++ ) 
			{
				fprintf( gv.QHSaveFile, "%.4e %.4e\n", qtemp[i],dPdlnT[i] );
			}
		}
		else 
		{
			fprintf( gv.QHSaveFile, "**** quantum heating was not used\n" );
		}
	}
	return;
}


/* initialize phiTilde */
STATIC void qheat_init(size_t nd,
		       /*@out@*/ vector<double>& phiTilde,  /* phiTilde[rfield.nflux_with_check] */
		       /*@out@*/ double *check)
{
	long i,
	  nz;
	double sum = 0.;

	/*@-redef@*/
	enum {DEBUG_LOC=false};
	/*@+redef@*/

	DEBUG_ENTRY( "qheat_init()" );

	/* sanity checks */
	ASSERT( gv.bin[nd].lgQHeat );
	ASSERT( nd < gv.bin.size() );

	*check = 0.;

	/* >>chng 01 nov 29, rfield.nflux -> gv.qnflux, PvH */
	/* >>chng 03 jan 26, gv.qnflux -> gv.bin[nd].qnflux, PvH */
	for( i=0; i < gv.bin[nd].qnflux; i++ )
	{
		phiTilde[i] = 0.;
	}

	/* fill in auxiliary array for quantum heating routine
	 * it reshuffles the input spectrum times the cross section to take
	 * the photo-electric effect into account. this prevents the quantum
	 * heating routine from having to calculate this effect over and over
	 * again; it can do a straight integration instead, making the code
	 * a lot simpler and faster. this initializes the array for non-ionizing
	 * energies, the reshuffling for higher energies is done in the next loop
	 * phiTilde has units events/H/s/cell at default depletion */

	double NegHeatingRate = 0.;

	for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
	{
		double check1 = 0.;
		ChargeBin& gptr = gv.bin[nd].chrg(nz);

		// rfield.nPositive may have increased since the last call to GrainDrive()
		// if so, arrays like gptr.fac1 would not be initialized up to nPositive
		long limit = min( rfield.nPositive, gptr.nfill );

		/* integrate over incident continuum for non-ionizing energies */
		for( i=0; i < min(gptr.ipThresInf,limit); i++ )
		{
			check1 += rfield.SummedCon[i]*gv.bin[nd].dstab1[i]*rfield.anu(i);
			phiTilde[i] += gptr.FracPop*rfield.SummedCon[i]*gv.bin[nd].dstab1[i];
		}

		/* >>chng 01 mar 02, use new expressions for grain cooling and absorption
		 * cross sections following the discussion with Joe Weingartner, PvH */
		for( i=gptr.ipThresInf; i < limit; i++ )
		{
			long ipLo2 = gptr.ipThresInfVal;
			double cs1 = ( i >= ipLo2 ) ? gv.bin[nd].dstab1[i]*gptr.yhat_primary[i] : 0.;

			check1 += rfield.SummedCon[i]*gptr.fac1[i];
			/* this accounts for the photons that are fully absorbed by grain */
			phiTilde[i] += gptr.FracPop*rfield.SummedCon[i]*MAX2(gv.bin[nd].dstab1[i]-cs1,0.);

			/* >>chng 01 oct 10, use bisection search to find ip. On C scale now */

			/* this accounts for photons that eject an electron from the valence band */
			if( cs1 > 0. )
			{
				/* we treat photo-ejection and all subsequent de-excitation cascades
				 * from the conduction/valence band as one simultaneous event */
				/* the condition cs1 > 0. assures i >= ipLo2 */
				/* ratio is number of ejected electrons per primary ionization */
				double ratio = ( gv.lgWD01 ) ? 1. : gptr.yhat[i]/gptr.yhat_primary[i];
				/* ehat is average energy of ejected electron at infinity */
				double ehat = gptr.ehat[i];
				double cool1, sign = 1.;
				realnum xx;

				if( gptr.DustZ <= -1 )
					cool1 = gptr.ThresSurf + gptr.PotSurf + ehat;
				else
					cool1 = gptr.ThresSurfVal + gptr.PotSurf + ehat;
				/* secondary electrons are assumed to have the same Elo and Ehi as the
				 * primary electrons that excite them. This neglects the threshold for
				 * the ejection of the secondary electron and can cause xx to become
				 * negative if Ehi is less than that threshold. To conserve energy we
				 * will simply assume a negative rate here. Since secondary electrons
				 * are generally not important this should have little impact on the
				 * overall temperature distribution */
				xx = rfield.anu(i) - (realnum)(ratio*cool1);
				if( xx < 0.f )
				{
					xx = -xx;
					sign = -1.;
				}
				long ipLo = rfield.ipointC( max(xx,rfield.emm()) );
				/* for grains in hard X-ray environments, the coarseness of the grid can
				 * lead to inaccuracies in the integral over phiTilde that would trip the
				 * sanity check in qheat(), here we correct for the energy mismatch */
				double corr = xx/rfield.anu(ipLo);
				phiTilde[ipLo] += sign*corr*gptr.FracPop*rfield.SummedCon[i]*cs1;
			}

			/* no need to account for photons that eject an electron from the conduction band */
			/* >>chng 01 dec 11, cool2 always equal to rfield.anu(i) -> no grain heating */
		}

		*check += gptr.FracPop*check1*EN1RYD*gv.bin[nd].cnv_H_pCM3;

		sum += gptr.FracPop*check1*EN1RYD*gv.bin[nd].cnv_H_pCM3;

		if( DEBUG_LOC )
		{
			double integral = 0.;
			for( i=0; i < gv.bin[nd].qnflux; i++ )
			{
				integral += phiTilde[i]*gv.bin[nd].cnv_H_pCM3*rfield.anu(i)*EN1RYD;
			}
			dprintf( ioQQQ, " integral test 1: integral %.6e %.6e\n", integral, sum );
		}

		/* add quantum heating due to recombination of electrons, subtract thermionic cooling */

		/* gptr.HeatingRate2 is net heating rate in erg/H/s at standard depl
		 * includes contributions for recombining electrons, autoionizing electrons
		 * subtracted by thermionic emissions here since it is inverse process
		 *
		 * NB - in extreme conditions this rate may be negative (if there
		 * is an intense radiation field leading to very hot grains, but no ionizing
		 * photons, hence very few free electrons). we assume that the photon rates
		 * are high enough under those circumstances to avoid phiTilde becoming negative,
		 * but we will check that in qheat1 anyway. */

		/* >>chng 03 nov 06, check for extremely low HeatingRate and save CPU time, pah_crash.in, PvH */
		if( gptr.HeatingRate2*gv.bin[nd].cnv_H_pCM3 > 0.05*CONSERV_TOL*gv.bin[nd].GrainHeat ) 
		{
			double Sum,ESum,DSum,E_av2,Corr;
			double fac = BOLTZMANN/EN1RYD*phycon.te;
			/* E0 is barrier that electron needs to overcome, zero for positive grains */
			/* >>chng 03 jan 23, added second term to correctly account for autoionizing states
			 *                   where ThresInfInc is negative, tested in small_grain.in, PvH */
			double E0 = -(MIN2(gptr.PotSurfInc,0.) + MIN2(gptr.ThresInfInc,0.));
			/* >>chng 01 mar 02, this should be energy gap between top electron and infinity, PvH */
			/* >>chng 01 nov 21, use correct barrier: ThresInf[nz] -> ThresInfInc[nz], PvH */
			/* >>chng 03 jan 23, replaced -E0 with MIN2(PotSurfInc[nz],0.), PvH */
			double Einf = gptr.ThresInfInc + MIN2(gptr.PotSurfInc,0.);
			/* this is average energy deposited by one event, in erg
			 * this value is derived from distribution assumed here, and may
			 * not be the same as HeatElectrons/(CollisionRateElectr*eta) !! */
			/* >>chng 01 nov 21, use correct barrier: ThresInf[nz] -> ThresInfInc[nz], PvH */
			/* >>chng 03 jan 23, replaced ThresInfInc[nz] with MAX2(ThresInfInc[nz],0.), PvH */
			double E_av = MAX2(gptr.ThresInfInc,0.)*EN1RYD + 2.*BOLTZMANN*phycon.te;
			/* this is rate in events/H/s at standard depletion */
			double rate = gptr.HeatingRate2/E_av;

			double ylo = -exp(-E0/fac);
			/* this is highest kinetic energy of electron that can be represented in phiTilde */
			/* >>chng 01 nov 29, rfield.nflux -> gv.qnflux, PvH */
			/* >>chng 03 jan 26, gv.qnflux -> gv.bin[nd].qnflux, PvH */
			double Ehi = rfield.anumax(gv.bin[nd].qnflux-1)-Einf;
			double yhi = ((E0-Ehi)/fac-1.)*exp(-Ehi/fac);
			/* renormalize rate so that integral over phiTilde*anu gives correct total energy */
			rate /= yhi-ylo;

			/* correct for fractional population of this charge state */
			rate *= gptr.FracPop;

			/* >>chng 03 jan 24, add code to correct for discretization errors, hotdust.in, PvH */
			vector<double> RateArr(gv.bin[nd].qnflux);
			Sum = ESum = DSum = 0.;

			/* >>chng 04 jan 21, replaced gv.bin[nd].qnflux -> gv.bin[nd].qnflux2, PvH */
			for( i=0; i < gv.bin[nd].qnflux2; i++ ) 
			{
				Ehi = rfield.anumax(i) - Einf;
				if( Ehi >= E0 ) 
				{
					/* Ehi is kinetic energy of electron at infinity */
					yhi = ((E0-Ehi)/fac-1.)*exp(-Ehi/fac);
					/* >>chng 01 mar 24, use MAX2 to protect against roundoff error, PvH */
					RateArr[i] = rate*MAX2(yhi-ylo,0.);
					Sum += RateArr[i];
					ESum += rfield.anu(i)*RateArr[i];
#					ifndef NDEBUG
					DSum += rfield.widflx(i)*RateArr[i];
#					endif
					ylo = yhi;
				}
				else
				{
					RateArr[i] = 0.;
				}
			}
			E_av2 = ESum/Sum*EN1RYD;
			ASSERT( fabs(E_av-E_av2) <= DSum/Sum*EN1RYD );
			Corr = E_av/E_av2;

			for( i=0; i < gv.bin[nd].qnflux2; i++ ) 
			{
				phiTilde[i] += RateArr[i]*Corr;
			}

			sum += gptr.FracPop*gptr.HeatingRate2*gv.bin[nd].cnv_H_pCM3;

			if( DEBUG_LOC )
			{
				double integral = 0.;
				for( i=0; i < gv.bin[nd].qnflux; i++ )
				{
					integral += phiTilde[i]*gv.bin[nd].cnv_H_pCM3*rfield.anu(i)*EN1RYD;
				}
				dprintf( ioQQQ, " integral test 2: integral %.6e %.6e\n", integral, sum );
			}
		}
		else
		{
			NegHeatingRate += gptr.FracPop*gptr.HeatingRate2*gv.bin[nd].cnv_H_pCM3;
		}
	}

	/* ============================================================================= */

	/* add quantum heating due to molecule/ion collisions */

	/* gv.bin[nd].HeatingRate1 is heating rate in erg/H/s at standard depl
	 * includes contributions from molecules/neutral atoms and recombining ions
	 *
	 * in fully ionized conditions electron heating rates will be much higher
	 * than ion and molecule rates since electrons are so much faster and grains
	 * tend to be positive. in non-ionized conditions the main contribution will
	 * come from neutral atoms and molecules, so it is appropriate to treat both
	 * the same. in fully ionized conditions we don't care since unimportant.
	 *
	 * NB - if grains are hotter than ambient gas, the heating rate may become negative.
	 * if photon rates are not high enough to prevent phiTilde from becoming negative,
	 * we will raise a flag while calculating the quantum heating in qheat1 */

	/* >>chng 03 nov 06, check for extremely low HeatingRate and save CPU time, PvH */
	if( gv.bin[nd].HeatingRate1*gv.bin[nd].cnv_H_pCM3 > 0.05*CONSERV_TOL*gv.bin[nd].GrainHeat )
	{
		/* limits for Taylor expansion of (1+x)*exp(-x) */
		/* these choices will assure only 6 digits precision */
		const double LIM2 = cbrt(3.e-6);
		const double LIM3 = powpq(8.e-6,1,4);
		/* if gas temperature is higher than grain temperature we will
		 * consider Maxwell-Boltzmann distribution of incoming particles
		 * and ignore distribution of outgoing particles, if grains
		 * are hotter than ambient gas, we use reverse treatment */
		double fac = BOLTZMANN/EN1RYD*MAX2(phycon.te,gv.bin[nd].tedust);
		/* this is average energy deposited/extracted by one event, in erg */
		double E_av = 2.*BOLTZMANN*MAX2(phycon.te,gv.bin[nd].tedust);
		/* this is rate in events/H/s at standard depletion */
		double rate = gv.bin[nd].HeatingRate1/E_av;

		double ylo = -1.;
		/* this is highest energy of incoming/outgoing particle that can be represented in phiTilde */
		/* >>chng 01 nov 29, rfield.nflux -> gv.qnflux, PvH */
		/* >>chng 03 jan 26, gv.qnflux -> gv.bin[nd].qnflux, PvH */
		double Ehi = rfield.anumax(gv.bin[nd].qnflux-1);
		double yhi = -(Ehi/fac+1.)*exp(-Ehi/fac);
		/* renormalize rate so that integral over phiTilde*anu gives correct total energy */
		rate /= yhi-ylo;

		for( i=0; i < gv.bin[nd].qnflux2; i++ ) 
		{
			/* Ehi is kinetic energy of incoming/outgoing particle
			 * we assume that Ehi-E0 is deposited/extracted from grain */
			/* Ehi = rfield.anumax(i); */
			double x = rfield.anumax(i)/fac;
			/* (1+x)*exp(-x) = 1 - 1/2*x^2 + 1/3*x^3 - 1/8*x^4 + O(x^5)
			 *               = 1 - Sum_n=2^infty (-x)^n/(n*(n-2)!)      */
			if( x > LIM3 )
				yhi = -(x+1.)*exp(-x);
			else if( x > LIM2 )
				yhi = -(((1./3.)*x - 0.5)*x*x + 1.);
			else
				yhi = -(1. - 0.5*x*x);

			/* >>chng 01 mar 24, use MAX2 to protect against roundoff error, PvH */
			phiTilde[i] += rate*MAX2(yhi-ylo,0.);
			ylo = yhi;
		}

		sum += gv.bin[nd].HeatingRate1*gv.bin[nd].cnv_H_pCM3;

		if( DEBUG_LOC )
		{
			double integral = 0.;
			for( i=0; i < gv.bin[nd].qnflux; i++ )
			{
				integral += phiTilde[i]*gv.bin[nd].cnv_H_pCM3*rfield.anu(i)*EN1RYD;
			}
			dprintf( ioQQQ, " integral test 3: integral %.6e %.6e\n", integral, sum );
		}
	}
	else
	{
		NegHeatingRate += gv.bin[nd].HeatingRate1*gv.bin[nd].cnv_H_pCM3;
	}

	/* here we account for the negative heating rates, we simply do that by scaling the entire
	 * phiTilde array down by a constant factor such that the total amount of energy is conserved
	 * This treatment assures that phiTilde never goes negative, which avoids problems further on */
	if( NegHeatingRate < 0. )
	{
		double scale_fac = (sum+NegHeatingRate)/sum;
		for( i=0; i < gv.bin[nd].qnflux; i++ )
			phiTilde[i] *= scale_fac;

		if( DEBUG_LOC )
		{
			sum += NegHeatingRate;

			double integral = 0.;
			for( i=0; i < gv.bin[nd].qnflux; i++ )
			{
				integral += phiTilde[i]*gv.bin[nd].cnv_H_pCM3*rfield.anu(i)*EN1RYD;
			}
			dprintf( ioQQQ, " integral test 4: integral %.6e %.6e\n", integral, sum );
		}
	}

	return;
}


/*******************************************************************************************
 *
 * GetProbDistr_LowLimit: main routine for calculating non-equilibrium heating of grains
 *
 * This routine implements the algorithm outlined in:
 * >>refer	grain	physics	Guhathakurtha & Draine, 1989, ApJ, 345, 230
 *
 * The original (fortran) version of the code was written by Kevin Volk.
 *
 * Heavily modified and adapted for new style grains by Peter van Hoof.
 *
 *******************************************************************************************/

STATIC void GetProbDistr_LowLimit(size_t nd,
				  double rel_tol,
				  double Umax,
				  double fwhm,
				  /*@in@*/ const vector<double>& Phi,    /* Phi[NQGRID]      */
				  /*@in@*/ const vector<double>& PhiDrv, /* PhiDrv[NQGRID]   */
				  /*@out@*/ vector<double>& qtemp,       /* qtemp[NQGRID]    */
				  /*@out@*/ vector<double>& qprob,       /* qprob[NQGRID]    */
				  /*@out@*/ vector<double>& dPdlnT,      /* dPdlnT[NQGRID]   */
				  /*@out@*/ long int *qnbin,
				  /*@out@*/ double *new_tmin,
				  long *nWideFail,
				  /*@out@*/ QH_Code *ErrorCode)
{
	bool lgAllNegSlope,
	  lgBoundErr;
	long int j,
	  k,
	  l,
	  nbad,
	  nbin,
	  nend=0,
	  nmax,
	  nok,
	  nstart=0,
	  nstart2=0;
	double dCool=0.,
	  dlnLdlnT,
	  dlnpdlnU,
	  fac = 0.,
	  maxVal,
	  NextStep,
	  qtmin1, 
	  RadCooling,
	  sum,
	  y;
	vector<double> delu(NQGRID);
	vector<double> Lambda(NQGRID);
	vector<double> p(NQGRID);
	vector<double> u1(NQGRID);


	DEBUG_ENTRY( "GetProbDistr_LowLimit()" );

	/* sanity checks */
	ASSERT( nd < gv.bin.size() );

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "   GetProbDistr_LowLimit called for grain %s\n", gv.bin[nd].chDstLab );
		fprintf( ioQQQ, "    got qtmin1 %.4e\n", gv.bin[nd].qtmin);
	}

	qtmin1 = gv.bin[nd].qtmin;
	qtemp[0] = qtmin1;
	/* u1 holds enthalpy in Ryd/grain */
	u1[0] = ufunct(qtemp[0],nd,&lgBoundErr);
	ASSERT( !lgBoundErr ); /* this should never happen */
	/* >>chng 00 mar 22, taken out factor 4, factored in hden and dstAbund
	 * interpolate in dstems array instead of integrating explicitly, by PvH */
	splint_safe(gv.dsttmp,gv.bin[nd].dstems,gv.bin[nd].dstslp2,NDEMS,log(qtemp[0]),&y,&lgBoundErr);
	ASSERT( !lgBoundErr ); /* this should never happen */
	/* Lambda holds the radiated energy for grains in this bin, in Ryd/s/grain */
	Lambda[0] = exp(y)*gv.bin[nd].cnv_H_pGR/EN1RYD;
	/* set up first step of integration */
	/* >>chng 01 nov 14, set to 2.*Lambda[0]/Phi[0] instead of u1[0],
	 * this choice assures that p[1] doesn't make a large jump from p[0], PvH */
	delu[0] = 2.*Lambda[0]/Phi[0];
	p[0] = PROB_CUTOFF_LO;
	dPdlnT[0] = p[0]*qtemp[0]*uderiv(qtemp[0],nd);
	RadCooling = 0.5*p[0]*Lambda[0]*delu[0];
	NextStep = 0.01*Lambda[0]/Phi[0];
	/* >>chng 03 nov 10, added extra safeguard against stepsize too small, PvH */
	if( NextStep < 5.*DBL_EPSILON*u1[0] )
	{
		*ErrorCode = MAX2(*ErrorCode,QH_STEP_FAIL);
		return;
	}

	nbad = 0;
	k = 0;

	*qnbin = 0;
	*new_tmin = qtmin1;
	lgAllNegSlope = true;
	maxVal = dPdlnT[0];
	nmax = 0;
	double p_max = p[0];

	/* this test neglects a negative contribution which is impossible to calculate, so it may
	 * fail to detect cases where the probability distribution starts dropping immediately,
	 * we will use a second test using the variable lgAllNegSlope below to catch those cases, PvH */
	spldrv_safe(gv.dsttmp,gv.bin[nd].dstems,gv.bin[nd].dstslp2,NDEMS,log(qtemp[0]),&dlnLdlnT,&lgBoundErr);
	ASSERT( !lgBoundErr ); /* this should never happen */
	dlnpdlnU = u1[0]*Phi[0]/Lambda[0] - dlnLdlnT*u1[0]/(qtemp[0]*uderiv(qtemp[0],nd));
	if( dlnpdlnU < 0. )
	{
		/* >>chng 03 nov 06, confirm this by integrating first step..., pah_crash.in, PvH */
		(void)TryDoubleStep(u1,delu,p,qtemp,Lambda,Phi,PhiDrv,NextStep,&dCool,p_max,k,nd,&lgBoundErr);
		dPdlnT[2] = p[2]*qtemp[2]*uderiv(qtemp[2],nd);

		if( dPdlnT[2] < dPdlnT[0] )
		{
			/* if dPdlnT starts falling immediately, 
			 * qtmin1 was too high and convergence is impossible */
			*ErrorCode = MAX2(*ErrorCode,QH_HIGH_FAIL);
			return;
		}
	}

	/* NB NB -- every break in this loop should set *ErrorCode (except for regular stop criterion) !! */
	for( l=0; l < MAX_LOOP; ++l )
	{
		double rerr = TryDoubleStep(u1,delu,p,qtemp,Lambda,Phi,PhiDrv,NextStep,&dCool,p_max,k,nd,&lgBoundErr);

		/* this happens if the grain temperature in qtemp becomes higher than GRAIN_TMAX
		 * nothing that TryDoubleStep returns can be trusted, so this check should be first */
		if( lgBoundErr )
		{
			nbad += 2;
			*ErrorCode = MAX2(*ErrorCode,QH_THIGH_FAIL);
			break;
		}

		/* estimate new stepsize */
		if( rerr > rel_tol*QHEAT_TOL )
		{
			nbad += 2;

			/* step is rejected, decrease stepsize and try again */
			NextStep *= sqrt(0.9*rel_tol*QHEAT_TOL/rerr);

			/* stepsize too small, this can happen at extreme low temperatures */
			if( NextStep < 5.*DBL_EPSILON*u1[k] )
			{
				*ErrorCode = MAX2(*ErrorCode,QH_STEP_FAIL);
				break;
			}

			continue;
		}
		else
		{
			/* step was OK, adjust stepsize */
			k += 2;

			p_max = max(p_max,p[k-1]);
			p_max = max(p_max,p[k]);

			/* >>chng 03 nov 10, safeguard against division by zero, PvH */
			NextStep *= MIN2(cbrt(0.9*rel_tol*QHEAT_TOL/MAX2(rerr,1.e-50)),4.);
			NextStep = MIN2(NextStep,Lambda[k]/Phi[0]);
		}

		dPdlnT[k-1] = p[k-1]*qtemp[k-1]*uderiv(qtemp[k-1],nd);
		dPdlnT[k] = p[k]*qtemp[k]*uderiv(qtemp[k],nd);

		lgAllNegSlope = lgAllNegSlope && ( dPdlnT[k] < dPdlnT[k-2] );

		if( dPdlnT[k-1] > maxVal )
		{
			maxVal = dPdlnT[k-1];
			nmax = k-1;
		}
		if( dPdlnT[k] > maxVal )
		{
			maxVal = dPdlnT[k];
			nmax = k;
		}

		RadCooling += dCool;

// 		if( nzone >= 24 && nd == 0 ) {
// 		printf(" k %ld T[k] %.6e U[k] %.6e p[k] %.6e dPdlnT[k] %.6e\n",k-1,qtemp[k-1],u1[k-1],p[k-1],dPdlnT[k-1]);
// 		printf(" k %ld T[k] %.6e U[k] %.6e p[k] %.6e dPdlnT[k] %.6e\n",k,qtemp[k],u1[k],p[k],dPdlnT[k]);
// 		}

		/* if qtmin is far too low, p[k] can become extremely large, exceeding
		 * even double precision range. the following check should prevent overflows */
		/* >>chng 01 nov 07, sqrt(DBL_MAX) -> sqrt(DBL_MAX/100.) so that sqrt(p[k]*p[k+1]) is safe */
		if( p[k] > sqrt(DBL_MAX/100.) ) 
		{
			*ErrorCode = MAX2(*ErrorCode,QH_LOW_FAIL);
			break;
		}

		/* this may catch a bug in the Compaq C compiler V6.3-025
		 * if this gets triggered, try compiling with -ieee */
		ASSERT( p[k] > 0. && dPdlnT[k] > 0. && RadCooling > 0. );

		/* do a check for negative slope and if there will be enough room to store results */
		if( k > 0 && k%NQTEST == 0 )
		{
			double wid, avStep, factor;
			/* >>chng 02 jul 31, do additional test for HIGH_FAIL,
			 * first test before main loop doesn't catch everything, PvH */
			if( lgAllNegSlope )
			{
				*ErrorCode = MAX2(*ErrorCode,QH_HIGH_FAIL);
				break;
			}

			/* this is a lower limit for the total width of the probability distr */
			/* >>chng 02 jan 30, corrected calculation of wid and avStep, PvH */
			wid = (sqrt(-log(PROB_CUTOFF_LO)/(4.*LN_TWO)) +
			       sqrt(-log(PROB_CUTOFF_HI)/(4.*LN_TWO)))*fwhm/Umax;
			avStep = log(u1[k]/u1[0])/(double)k;
			/* make an estimate for the number of steps needed */
			/* >>chng 02 jan 30, included factor 1.5 because stepsize increases near peak, PvH */
			/* >>chng 02 aug 06, changed 1.5 to sliding scale because of laqheat2.in test, PvH */
			factor = 1.1 + 3.9*(1.0 - sqrt((double)k/(double)NQGRID));
			if( wid/avStep > factor*(double)NQGRID )
			{
				*ErrorCode = MAX2(*ErrorCode,QH_ARRAY_FAIL);
				break;
			}
		}

		/* if we run out of room to store results, do regular break
		 * the code below will sort out if integration is valid or not */
		if( k >= NQGRID-2 )
		{
			*ErrorCode = MAX2(*ErrorCode,QH_ARRAY_FAIL);
			break;
		}

		/* force thermal equilibrium of the grains */
		fac = RadCooling*gv.bin[nd].cnv_GR_pCM3*EN1RYD/gv.bin[nd].GrainHeat;

		/* this is regular stop criterion */
		if( dPdlnT[k] < dPdlnT[k-1] && dPdlnT[k]/fac < PROB_CUTOFF_HI )
		{
			break;
		}
	}

	if( l == MAX_LOOP )
		*ErrorCode = MAX2(*ErrorCode,QH_LOOP_FAIL);

	nok = k;
	nbin = k+1;

	/* there are insufficient bins to attempt rebinning */
	if( *ErrorCode < QH_NO_REBIN && nbin < NQMIN ) 
		*ErrorCode = MAX2(*ErrorCode,QH_NBIN_FAIL);

	/* >>chng 02 aug 07, do some basic checks on the distribution first */
	if( *ErrorCode < QH_NO_REBIN )
		ScanProbDistr(u1,dPdlnT,nbin,maxVal,nmax,qtmin1,&nstart,&nstart2,&nend,nWideFail,ErrorCode);

	if( *ErrorCode >= QH_NO_REBIN )
	{
		return;
	}

	for( j=0; j < nbin; j++ )
	{
		p[j] /= fac;
		dPdlnT[j] /= fac;
	}
	RadCooling /= fac;

	/* >>chng 02 aug 08, moved RebinQHeatResults from here further down, this improves new_tmin estimate */
	*new_tmin = 0.;
	for( j=nstart; j < nbin; j++ ) 
	{
		if( dPdlnT[j] < PROB_CUTOFF_LO ) 
		{
			*new_tmin = qtemp[j];
		}
		else 
		{
			if( j == nstart )
			{
				/* if dPdlnT[nstart] is too high, but qtmin1 is already close to GRAIN_TMIN,
				 * then don't bother signaling a QH_BOUND_FAIL. grains below GRAIN_TMIN have the
				 * peak of their thermal emission beyond 3 meter, so they really are irrelevant
				 * since free-free emission from electrons will drown this grain emission... */
				if( dPdlnT[j] > SAFETY*PROB_CUTOFF_LO && qtmin1 > 1.2*GRAIN_TMIN )
					*ErrorCode = MAX2(*ErrorCode,QH_BOUND_FAIL);

				/* >>chng 02 aug 11, use nstart2 for more reliable extrapolation */
				if( dPdlnT[nstart2] < 0.999*dPdlnT[nstart2+NSTARTUP] ) 
				{
					/* >>chng 02 aug 09, new formula for extrapolating new_tmin, PvH */
					/* this assumes that at low temperatures the behaviour
					 * is as follows:   dPdlnT(T) = C1 * exp( -C2/T**3 ) */
					double T1 = qtemp[nstart2];
					double T2 = qtemp[nstart2+NSTARTUP];
					double delta_y = log(dPdlnT[nstart2+NSTARTUP]/dPdlnT[nstart2]);
					double c2 = delta_y/(1./pow3(T1)-1./pow3(T2));
					double help = c2/pow3(T1) + log(dPdlnT[nstart2]/PROB_CUTOFF_LO);
					*new_tmin = cbrt(c2/help);
				}

				/* >>chng 04 nov 09, in case of lower bound failure, assure qtmin is lowered, PvH */ 
				if( dPdlnT[j] > SAFETY*PROB_CUTOFF_LO && *new_tmin >= qtmin1 )
				{
					double delta_x = log(qtemp[nstart2+NSTARTUP]/qtemp[nstart2]);
					double delta_y = log(dPdlnT[nstart2+NSTARTUP]/dPdlnT[nstart2]);
					delta_x *= log(PROB_CUTOFF_LO/dPdlnT[nstart2])/delta_y;
					*new_tmin = qtemp[nstart2]*exp(delta_x);
					if( *new_tmin < qtmin1 )
						/* in general this estimate will be too low -> use geometric mean */
						*new_tmin = sqrt( *new_tmin*qtmin1 );
					else
						/* last resort... */
						*new_tmin = qtmin1/DEF_FAC;
				}
			}
			break;
		}
	}
	*new_tmin = MAX3(*new_tmin,qtmin1/DEF_FAC,GRAIN_TMIN);

	ASSERT( *new_tmin < gv.bin[nd].tedust );

	/* >>chng 02 jan 30, prevent excessive looping when prob distribution simply won't fit, PvH */
	if( dPdlnT[nbin-1] > SAFETY*PROB_CUTOFF_HI )
	{
		if( *ErrorCode == QH_ARRAY_FAIL || *ErrorCode == QH_LOOP_FAIL )
		{
			++(*nWideFail);

			if( *nWideFail < WIDE_FAIL_MAX )
			{
				/* this indicates that low end was OK, but we ran out of room
				 * to store the high end -> try GetProbDistr_HighLimit instead */
				*ErrorCode = MAX2(*ErrorCode,QH_DELTA_FAIL);
			}
			else
			{
				*ErrorCode = MAX2(*ErrorCode,QH_WIDE_FAIL);
			}
		}
		else
		{
			*ErrorCode = MAX2(*ErrorCode,QH_BOUND_FAIL);
		}
	}

	/* >>chng 01 may 11, rebin the quantum heating results
	 *
	 * for grains in intense radiation fields, the code needs high resolution for stability
	 * and therefore produces lots of small bins, even though the grains never make large
	 * excurions from the equilibrium temperature; adding in the resulting spectra in RT_diffuse
	 * takes up an excessive amount of CPU time where most CPU is spent on grains for which
	 * the quantum treatment matters least, and moreover on temperature bins with very low
	 * probability; rebinning the results on a coarser grid should help reduce the overhead */
	/* >>chng 02 aug 07, use nstart and nend while rebinning */

	nbin = RebinQHeatResults(nd,nstart,nend,p,qtemp,qprob,dPdlnT,u1,delu,Lambda,ErrorCode);

	/* >>chng 01 jul 13, add fail-safe for failure in RebinQHeatResults */
	if( nbin == 0 )
	{
		return;
	}

	*qnbin = nbin;

	sum = 0.;
	for( j=0; j < nbin; j++ )
	{
		sum += qprob[j];
	}

	/* the fact that the probability normalization fails may indicate that the distribution is
	 * too close to a delta function to resolve, but another possibility is that the radiation
	 * field is extremely diluted allowing a sizable fraction of the grains to cool below
	 * GRAIN_TMIN. In the latter case we don't raise QH_CONV_FAIL since these very cool grains
	 * only contribute at very long radio wavelengths (more than 1 meter) */
	if( fabs(sum-1.) > PROB_TOL && qtmin1 > 1.2*GRAIN_TMIN )
		*ErrorCode = MAX2(*ErrorCode,QH_CONV_FAIL);

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ,
			 "    zone %ld %s nbin %ld nok %ld nbad %ld total prob %.4e rel_tol %.3e new_tmin %.4e\n",
			 nzone,gv.bin[nd].chDstLab,nbin,nok,nbad,sum,rel_tol,*new_tmin );
	}
	return;
}


/* try two consecutive integration steps using stepsize "step/2." (yielding p[k]),
 * and also one double integration step using stepsize "step" (yielding p2k).
 * the difference fabs(p2k-p[k])/(3.*p[k]) can be used to estimate the relative
 * accuracy of p[k] and will be used to adapt the stepsize to an optimal value */
STATIC double TryDoubleStep(vector<double>& u1,
			    vector<double>& delu,
			    vector<double>& p,
			    vector<double>& qtemp,
			    vector<double>& Lambda,
			    const vector<double>& Phi,
			    const vector<double>& PhiDrv,
			    double step,
			    /*@out@*/ double *cooling,
			    double p_max,
			    long index,
			    size_t nd,
			    /*@out@*/ bool *lgBoundFail)
{
	long i,
	  j,
	  jlo,
	  k=-1000;
	double bval_jk,
	  cooling2,
	  p2k,
	  RelErrCool,
	  RelErrPk,
	  sum,
	  sum2 = -DBL_MAX,
	  trap1,
	  trap12 = -DBL_MAX,
	  trap2,
	  uhi,
	  ulo,
	  umin,
	  y;

	DEBUG_ENTRY( "TryDoubleStep()" );

	/* sanity checks */
	ASSERT( index >= 0 && index < NQGRID-2 && nd < gv.bin.size() && step > 0. );

	ulo = rfield.anumin(0);
	/* >>chng 01 nov 29, rfield.nflux -> gv.qnflux, PvH */
	/* >>chng 03 jan 26, gv.qnflux -> gv.bin[nd].qnflux, PvH */
	uhi = rfield.anumax(gv.bin[nd].qnflux-1);

	/* >>chng 01 nov 21, skip initial bins if they have very low probability */
	jlo = 0;
	while( p[jlo] < PROB_CUTOFF_LO*p_max )
		jlo++;

	for( i=1; i <= 2; i++ )
	{
		bool lgErr;
		long ipLo = 0;
		// using gv.bin[nd].qnflux is OK since gv.bin[nd].qnflux < gv.bin[nd].nflux_with_check
		long ipHi = gv.bin[nd].qnflux;
		k = index + i;
		delu[k] = step/2.;
		u1[k] = u1[k-1] + delu[k];
		qtemp[k] = inv_ufunct(u1[k],nd,lgBoundFail);
		splint_safe(gv.dsttmp,gv.bin[nd].dstems,gv.bin[nd].dstslp2,NDEMS,log(qtemp[k]),&y,&lgErr);
		*lgBoundFail = *lgBoundFail || lgErr;
		Lambda[k] = exp(y)*gv.bin[nd].cnv_H_pGR/EN1RYD;

		sum = sum2 = 0.;
		trap1 = trap2 = trap12 = 0.;

		/* this loop uses up a large fraction of the total CPU time, it should be well optimized */
		for( j=jlo; j < k; j++ )
		{
			umin = u1[k] - u1[j];

			if( umin >= uhi ) 
			{
				/* for increasing j, umin will decrease monotonically. If ( umin > uhi ) holds for
				 * the current value of j, it must have held for the previous value as well. Hence
				 * both trap1 and trap2 are zero at this point and we would only be adding zero
				 * to sum. Therefore we skip this step to save CPU time */
				continue;
			}
			else if( umin > ulo )
			{
				/* do a bisection search such that rfield.anumin(ipLo) <= umin < rfield.anumin(ipHi)
				 * explicit bisection search is faster, which is important here to save CPU time.
				 * on the first iteration ipLo equals 0 and the first while loop will be skipped;
				 * after that umin is monotonically decreasing, and ipHi is retained from the
				 * previous iteration since it is a valid upper limit; ipLo will equal ipHi-1 */
				long ipStep = 1;
				/* >>chng 03 feb 03 rjrw: hunt for lower bracket */
				while( rfield.anumin(ipLo) > umin )
				{
					ipHi = ipLo;
					ipLo -= ipStep;
					if( ipLo <= 0 ) 
					{
						ipLo = 0;
						break;
					}
					ipStep *= 2;
				}
				/* now do regular bisection search */
				while( ipHi-ipLo > 1 )
				{
					long ipMd = (ipLo+ipHi)/2;
					if( rfield.anumin(ipMd) > umin )
						ipHi = ipMd;
					else
						ipLo = ipMd;
				}
				/* Phi[i] is integral of PhiDrv from exactly rfield.anumin(i) to infinity */
				bval_jk = Phi[ipLo] + (umin - rfield.anumin(ipLo))*PhiDrv[ipLo];
			}
			else
			{
				bval_jk = Phi[0];
			}

			/* these two quantities are needed to take double step from index -> index+2 */
			trap12 = trap1;
			sum2 = sum;

			/* bval_jk*gv.bin[nd].cnv_CM3_pGR is the total excitation rate from j to k and
			 * higher due to photon absorptions and particle collisions, it already implements
			 * Eq. 2.17 of Guhathakurtha & Draine, in events/grain/s */
			/* >>chng 00 mar 27, factored in hden (in Phi), by PvH */
			/* >>chng 00 mar 29, add in contribution due to particle collisions, by PvH */
			/* >>chng 01 mar 30, moved multiplication of bval_jk with gv.bin[nd].cnv_CM3_pGR
			 *   out of loop, PvH */
			trap2 = p[j]*bval_jk;
			/* Trapezoidal rule, multiplication with factor 0.5 is done outside loop */
			sum += (trap1 + trap2)*delu[j];
			trap1 = trap2;
		}

		/* >>chng 00 mar 27, multiplied with delu, by PvH */
		/* >>chng 00 apr 05, taken out Lambda[0], improves convergence at low end dramatically!, by PvH */
		/* qprob[k] = sum*gv.bin[nd].cnv_CM3_pGR*delu[k]/(Lambda[k] - Lambda[0]); */
		/* this expression includes factor 0.5 from trapezoidal rule above */
		/* p[k] = 0.5*(sum + (trap1 + p[k]*Phi[0])*delu[k])/Lambda[k] */
		p[k] = (sum + trap1*delu[k])/(2.*Lambda[k] - Phi[0]*delu[k]);

		// total failure -> force next step to be smaller
		if( p[k] <= 0. )
			return 3.*QHEAT_TOL;
	}

	/* this is estimate for p[k] using one double step of size "step" */
	p2k = (sum2 + trap12*step)/(2.*Lambda[k] - Phi[0]*step);

	// total failure -> force next step to be smaller
	if( p2k <= 0. )
		return 3.*QHEAT_TOL;

	RelErrPk = fabs(p2k-p[k])/p[k];

	/* this is radiative cooling due to the two probability bins we just added */
	/* simple trapezoidal rule will not do here, RelErrCool would never converge */
	double z[8];
	vlog(z, u1[k-2], u1[k-1], u1[k], p[k-2]*Lambda[k-2], p[k-1]*Lambda[k-1], p[k]*Lambda[k], p2k*Lambda[k], 1.);
	*cooling = log_integral(u1[k-2],p[k-2]*Lambda[k-2],u1[k-1],p[k-1]*Lambda[k-1],z[0],z[3],z[1],z[4]);
	*cooling += log_integral(u1[k-1],p[k-1]*Lambda[k-1],u1[k],p[k]*Lambda[k],z[1],z[4],z[2],z[5]);

	/* same as cooling, but now for double step of size "step" */
	cooling2 = log_integral(u1[k-2],p[k-2]*Lambda[k-2],u1[k],p2k*Lambda[k],z[0],z[3],z[2],z[6]);

	/* p[0] is not reliable, so ignore convergence test on cooling on first step */
	RelErrCool = ( index > 0 ) ? fabs(cooling2-(*cooling))/(*cooling) : 0.;

//	dprintf( ioQQQ, " TryDoubleStep k %ld p[k-1] %.4e p[k] %.4e p2k %.4e\n",k,p[k-1],p[k],p2k );
	/* error scales as O(step^3), so this is relative accuracy of p[k] or cooling */
	return MAX2(RelErrPk,RelErrCool)/3.;
}


/* calculate logarithmic integral from (xx1,yy1) to (xx2,yy2) */
STATIC double log_integral(double xx1,
			   double yy1,
			   double xx2,
			   double yy2,
			   double log_xx1,
			   double log_yy1,
			   double log_xx2,
			   double log_yy2)
{
	DEBUG_ENTRY( "log_integral()" );

	double xx = log_xx2 - log_xx1;
	double eps = xx + log_yy2 - log_yy1;
	if( fabs(eps) < 1.e-4 )
	{
		return xx1*yy1*xx*(((eps/24. + 1./6.)*eps + 0.5)*eps + 1.);
	}
	else
	{
		return (xx2*yy2 - xx1*yy1)*xx/eps;
	}
}


/* scan the probability distribution for valid range */
STATIC void ScanProbDistr(const vector<double>& u1,      /* u1[nbin] */
			  const vector<double>& dPdlnT,  /* dPdlnT[nbin] */
			  long nbin,
			  double maxVal,
			  long nmax,
			  double qtmin1,
			  /*@out@*/long *nstart,
			  /*@out@*/long *nstart2,
			  /*@out@*/long *nend,
			  long *nWideFail,
			  QH_Code *ErrorCode)
{
	bool lgSetLo,
	  lgSetHi;
	long i;
	double deriv_max,
	  minVal;
	const double MIN_FAC_LO = 1.e4;
	const double MIN_FAC_HI = 1.e4;

	DEBUG_ENTRY( "ScanProbDistr()" );

	/* sanity checks */
	ASSERT( nbin > 0 && nmax >= 0 && nmax < nbin && maxVal > 0. );

	/* sometimes the probability distribution will start falling before settling on
	 * a rising slope. In such a case nstart will point to the first rising point,
	 * while nstart2 will point to the point with the steepest derivative on the
	 * rising slope. The code will use the distribution from nstart to nend as a
	 * valid probability distribution, but the will use the points near nstart2
	 * to extrapolate a new value for qtmin if needed */
	minVal = maxVal;
	*nstart = nmax;
	for( i=nmax; i >= 0; --i )
	{
		if( dPdlnT[i] < minVal )
		{
			*nstart = i;
			minVal = dPdlnT[i];
		}
	}
	deriv_max = 0.;
	*nstart2 = nmax;
	for( i=nmax; i > *nstart; --i )
	{
		double deriv = log(dPdlnT[i]/dPdlnT[i-1])/log(u1[i]/u1[i-1]);
		if( deriv > deriv_max )
		{
			*nstart2 = i-1;
			deriv_max = deriv;
		}
	}
	*nend = nbin-1;

	/* now do quality control; these checks are more stringent than the ones in GetProbDistr_LowLimit */
	lgSetLo = ( nmax >= *nend || maxVal/dPdlnT[*nend] < MIN_FAC_HI );
	/* >>chng 03 jan 22, prevent problems if both dPdlnT and its derivative are continuously rising,
	 *                   in which case both lgSetLo and lgSetHi are set and QH_WIDE_FAIL is triggered;
	 *                   this can happen if qtmin is far too low; triggered by pahtest.in, PvH */
	if( lgSetLo )
		/* use relaxed test if lgSetLo is already set */
		lgSetHi = ( nmax <= *nstart || maxVal/dPdlnT[*nstart] < MIN_FAC_LO );
	else
		lgSetHi = ( nmax <= *nstart2 || maxVal/dPdlnT[*nstart2] < MIN_FAC_LO );

	if( lgSetLo && lgSetHi )
	{
		++(*nWideFail);

		if( *nWideFail >= WIDE_FAIL_MAX )
			*ErrorCode = MAX2(*ErrorCode,QH_WIDE_FAIL);
	}

	if( lgSetLo )
		*ErrorCode = MAX2(*ErrorCode,QH_LOW_FAIL);

	/* if dPdlnT[nstart] is too high, but qtmin1 is already close to GRAIN_TMIN,
	 * then don't bother signaling a QH_HIGH_FAIL. grains below GRAIN_TMIN have the
	 * peak of their thermal emission beyond 3 meter, so they really are irrelevant
	 * since free-free emission from electrons will drown this grain emission... */
	if( lgSetHi && qtmin1 > 1.2*GRAIN_TMIN )
		*ErrorCode = MAX2(*ErrorCode,QH_HIGH_FAIL);

	/* there are insufficient bins to attempt rebinning */
	if( *ErrorCode < QH_NO_REBIN && (*nend - *nstart) < NQMIN ) 
		*ErrorCode = MAX2(*ErrorCode,QH_NBIN_FAIL);

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "    ScanProbDistr nstart %ld nstart2 %ld nend %ld nmax %ld maxVal %.3e",
			 *nstart,*nstart2,*nend,nmax,maxVal );
		fprintf( ioQQQ, " dPdlnT[nstart] %.3e dPdlnT[nstart2] %.3e dPdlnT[nend] %.3e code %d\n",
			 dPdlnT[*nstart],dPdlnT[*nstart2],dPdlnT[*nend],*ErrorCode );
	}

	if( *ErrorCode >= QH_NO_REBIN )
	{
		*nstart = -1;
		*nstart2 = -1;
		*nend = -2;
	}
	return;
}


/* rebin the quantum heating results to speed up RT_diffuse */
STATIC long RebinQHeatResults(size_t nd,
			      long nstart,
			      long nend,
			      vector<double>& p,
			      vector<double>& qtemp,
			      vector<double>& qprob,
			      vector<double>& dPdlnT,
			      vector<double>& u1,
			      vector<double>& delu,
			      vector<double>& Lambda,
			      QH_Code *ErrorCode)
{
	long i,
	  newnbin;
	double fac,
	  help,
	  mul_fac,
	  PP1,
	  PP2,
	  RadCooling,
	  T1,
	  T2,
	  Ucheck,
	  uu1,
	  uu2;

	DEBUG_ENTRY( "RebinQHeatResults()" );

	/* sanity checks */
	ASSERT( nd < gv.bin.size() );
	/* >>chng 02 aug 07, changed oldnbin -> nstart..nend */
	ASSERT( nstart >= 0 && nstart < nend && nend < NQGRID );

	/* leading entries may have become very small or zero -> skip */
	for( i=nstart; i <= nend && dPdlnT[i] < PROB_CUTOFF_LO; i++ ) {}

	/* >>chng 04 oct 17, add fail-safe to keep lint happy, but this should never happen... */
	if( i >= NQGRID )
	{
		*ErrorCode = MAX2(*ErrorCode,QH_REBIN_FAIL);
		return 0;
	}

	vector<double> new_delu(NQGRID);
	vector<double> new_dPdlnT(NQGRID);
	vector<double> new_Lambda(NQGRID);
	vector<double> new_p(NQGRID);
	vector<double> new_qprob(NQGRID);
	vector<double> new_qtemp(NQGRID);
	vector<double> new_u1(NQGRID);
	double z[8];

	newnbin = 0;

	T1 = qtemp[i];
	PP1 = p[i];
	uu1 = u1[i];

	/* >>chng 04 feb 01, change 2.*NQMIN -> 1.5*NQMIN, PvH */
	help = pow(qtemp[nend]/qtemp[i],1./(1.5*NQMIN));
	mul_fac = MIN2(QT_RATIO,help);

	Ucheck = u1[i];
	RadCooling = 0.;

	while( i < nend )
	{
		bool lgBoundErr;
		bool lgDone= false;
		double s0 = 0.;
		double s1 = 0.;
		double wid = 0.;
		double xx,y;

		T2 = T1*mul_fac;

		do 
		{
			double p1,p2,L1,L2,frac,slope;
			if( qtemp[i] <= T1 && T1 <= qtemp[i+1] )
			{
				/* >>chng 01 nov 15, copy uu2 into uu1 instead, PvH */
				/* uu1 = ufunct(T1,nd); */
				double xrlog = log(qtemp[i+1]/qtemp[i]);
				if( xrlog > 0. )
				{
					frac = log(T1/qtemp[i]);
					slope = log(p[i+1]/p[i])/xrlog;
					p1 = p[i]*exp(frac*slope);
					slope = log(Lambda[i+1]/Lambda[i])/xrlog;
					L1 = Lambda[i]*exp(frac*slope);
				}
				else
				{
					/* pathological case where slope is extremely steep (daniela.in) */
					p1 = sqrt(p[i]*p[i+1]);
					L1 = sqrt(Lambda[i]*Lambda[i+1]);
				}
			}
			else
			{
				/* >>chng 01 nov 15, copy uu2 into uu1 instead, PvH */
				/* uu1 = u1[i]; */
				p1 = p[i];
				L1 = Lambda[i];
			}
			if( qtemp[i] <= T2 && T2 <= qtemp[i+1] )
			{
				/* >>chng 02 apr 30, make sure this doesn't point beyond valid range, PvH */
				help = ufunct(T2,nd,&lgBoundErr);
				uu2 = MIN2(help,u1[i+1]);
				ASSERT( !lgBoundErr ); /* this should never be triggered */
				double xrlog = log(qtemp[i+1]/qtemp[i]);
				if( xrlog > 0. )
				{
					frac = log(T2/qtemp[i]);
					slope = log(p[i+1]/p[i])/xrlog;
					p2 = p[i]*exp(frac*slope);
					slope = log(Lambda[i+1]/Lambda[i])/xrlog;
					L2 = Lambda[i]*exp(frac*slope);
				}
				else
				{
					/* pathological case where slope is extremely steep */
					p2 = sqrt(p[i]*p[i+1]);
					L2 = sqrt(Lambda[i]*Lambda[i+1]);
				}
				lgDone = true;
			}
			else
			{
				uu2 = u1[i+1];
				p2 = p[i+1];
				L2 = Lambda[i+1];
				/* >>chng 01 nov 15, this caps the range in p(U) integrated in one bin
				 * it helps avoid spurious QH_BOUND_FAIL's when flank is very steep, PvH */
				if( MAX2(p2,PP1)/MIN2(p2,PP1) > sqrt(SAFETY) )
				{
					lgDone = true;
					T2 = qtemp[i+1];
				}
				++i;
			}
			PP2 = p2;
			wid += uu2 - uu1;
			/* sanity check */
			ASSERT( wid >= 0. );
			vlog(z,uu1,uu2,p1,p2,p1*L1,p2*L2,1.,1.);
			s0 += log_integral(uu1,p1,uu2,p2,z[0],z[2],z[1],z[3]);
			s1 += log_integral(uu1,p1*L1,uu2,p2*L2,z[0],z[4],z[1],z[5]);
			uu1 = uu2;

		} while( i < nend && ! lgDone );

		/* >>chng 01 nov 14, if T2 == qtemp[oldnbin-1], the code will try another iteration
		 * break here to avoid zero divide, the assert on Ucheck tests if we are really finished */
		/* >>chng 01 dec 04, change ( s0 == 0. ) to ( s0 <= 0. ), PvH */
		if( s0 <= 0. )
		{
			ASSERT( wid == 0. );
			break;
		}

		new_qprob[newnbin] = s0;
		new_Lambda[newnbin] = s1/s0;
		xx = log(new_Lambda[newnbin]*EN1RYD*gv.bin[nd].cnv_GR_pH);
		splint_safe(gv.bin[nd].dstems,gv.dsttmp,gv.bin[nd].dstslp,NDEMS,xx,&y,&lgBoundErr);
		ASSERT( !lgBoundErr ); /* this should never be triggered */
		new_qtemp[newnbin] = exp(y);
		new_u1[newnbin] = ufunct(new_qtemp[newnbin],nd,&lgBoundErr);
		ASSERT( !lgBoundErr ); /* this should never be triggered */
		new_delu[newnbin] = wid;
		new_p[newnbin] = new_qprob[newnbin]/new_delu[newnbin];
		new_dPdlnT[newnbin] = new_p[newnbin]*new_qtemp[newnbin]*uderiv(new_qtemp[newnbin],nd);

		Ucheck += wid;
		RadCooling += new_qprob[newnbin]*new_Lambda[newnbin];

		T1 = T2;
		PP1 = PP2;
		++newnbin;
	}

	/* >>chng 01 jul 13, add fail-safe */
	if( newnbin < NQMIN )
	{
		*ErrorCode = MAX2(*ErrorCode,QH_REBIN_FAIL);
		return 0;
	}

	fac = RadCooling*EN1RYD*gv.bin[nd].cnv_GR_pCM3/gv.bin[nd].GrainHeat;

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "     RebinQHeatResults found tol1 %.4e tol2 %.4e\n",
			 fabs(u1[nend]/Ucheck-1.), fabs(fac-1.) );
	}

	/* do quality control */
	/* >>chng 02 apr 30, tighten up check, PvH */
	ASSERT( fabs(u1[nend]/Ucheck-1.) < 10.*sqrt((double)(nend-nstart+newnbin))*DBL_EPSILON );

	if( fabs(fac-1.) > CONSERV_TOL )
		*ErrorCode = MAX2(*ErrorCode,QH_CONV_FAIL);

	for( i=0; i < newnbin; i++ )
	{
		/* renormalize the distribution to assure energy conservation */
		p[i] = new_p[i]/fac;
		qtemp[i] = new_qtemp[i];
		qprob[i] = new_qprob[i]/fac;
		dPdlnT[i] = new_dPdlnT[i]/fac;
		u1[i] = new_u1[i];
		delu[i] = new_delu[i];
		Lambda[i] = new_Lambda[i];

		/* sanity checks */
		ASSERT( qtemp[i] > 0. && qprob[i] > 0. );

/*  		printf(" rk %ld T[k] %.6e U[k] %.6e p[k] %.6e dPdlnT[k] %.6e\n",i,qtemp[i],u1[i],p[i],dPdlnT[i]); */
	}
	return newnbin;
}


/* calculate approximate probability distribution in high intensity limit */
STATIC void GetProbDistr_HighLimit(long nd,
				   double TolFac,
				   double Umax,
				   double fwhm,
				   /*@out@*/vector<double>& qtemp,
				   /*@out@*/vector<double>& qprob,
				   /*@out@*/vector<double>& dPdlnT,
				   /*@out@*/double *tol,
				   /*@out@*/long *qnbin,
				   /*@out@*/double *new_tmin,
				   /*@out@*/QH_Code *ErrorCode)
{
	bool lgBoundErr,
	  lgErr;
	long i,
	  nbin;
	double c1,
	  c2,
	  delu[NQGRID],
	  fac,
	  fac1,
	  fac2,
	  help1,
	  help2,
	  L1,
	  L2,
	  Lambda[NQGRID],
	  mul_fac,
	  p[NQGRID],
	  p1,
	  p2,
	  RadCooling,
	  sum,
	  T1,
	  T2,
	  Tlo,
	  Thi,
	  Ulo,
	  Uhi,
	  uu1,
	  uu2,
	  xx,
	  y;

	DEBUG_ENTRY( "GetProbDistr_HighLimit()" );

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "   GetProbDistr_HighLimit called for grain %s\n", gv.bin[nd].chDstLab );
	}

	c1 = sqrt(4.*LN_TWO/PI)/fwhm*exp(-pow2(fwhm/Umax)/(16.*LN_TWO));
	c2 = 4.*LN_TWO*pow2(Umax/fwhm);

	fac1 = fwhm/Umax*sqrt(-log(PROB_CUTOFF_LO)/(4.*LN_TWO));
	/* >>chng 03 nov 10, safeguard against underflow, PvH */
	help1 = Umax*exp(-fac1);
	help2 = exp(gv.bin[nd].DustEnth[0]);
	Ulo = MAX2(help1,help2);
	/* >>chng 03 jan 28, ignore lgBoundErr on lower boundary, low-T grains have negigible emission, PvH */
	Tlo = inv_ufunct(Ulo,nd,&lgBoundErr);

	fac2 = fwhm/Umax*sqrt(-log(PROB_CUTOFF_HI)/(4.*LN_TWO));
	/* >>chng 03 nov 10, safeguard against overflow, PvH */
	if( fac2 > log(DBL_MAX/10.) )
	{
		*ErrorCode = MAX2(*ErrorCode,QH_WIDE_FAIL);
		return;
	}
	Uhi = Umax*exp(fac2);
	Thi = inv_ufunct(Uhi,nd,&lgBoundErr);

	nbin = 0;

	T1 = Tlo;
	uu1 = ufunct(T1,nd,&lgErr);
	lgBoundErr = lgBoundErr || lgErr;
	help1 = log(uu1/Umax);
	p1 = c1*exp(-c2*pow2(help1));
	splint_safe(gv.dsttmp,gv.bin[nd].dstems,gv.bin[nd].dstslp2,NDEMS,log(T1),&y,&lgErr);
	lgBoundErr = lgBoundErr || lgErr;
	L1 = exp(y)*gv.bin[nd].cnv_H_pGR/EN1RYD;

	/* >>chng 03 nov 10, safeguard against underflow, PvH */
	if( uu1*p1*L1 < 1.e5*DBL_MIN )
	{
		*ErrorCode = MAX2(*ErrorCode,QH_WIDE_FAIL);
		return;
	}

	/* >>chng 04 feb 01, change 2.*NQMIN -> 1.2*NQMIN, PvH */
	help1 = pow(Thi/Tlo,1./(1.2*NQMIN));
	mul_fac = MIN2(QT_RATIO,help1);

	sum = 0.;
	RadCooling = 0.;

	double z[8];

	do 
	{
		double s0,s1,wid;

		T2 = T1*mul_fac;
		uu2 = ufunct(T2,nd,&lgErr);
		lgBoundErr = lgBoundErr || lgErr;
		help1 = log(uu2/Umax);
		p2 = c1*exp(-c2*pow2(help1));
		splint_safe(gv.dsttmp,gv.bin[nd].dstems,gv.bin[nd].dstslp2,NDEMS,log(T2),&y,&lgErr);
		lgBoundErr = lgBoundErr || lgErr;
		L2 = exp(y)*gv.bin[nd].cnv_H_pGR/EN1RYD;

		wid = uu2 - uu1;
		vlog(z,uu1,uu2,p1,p2,p1*L1,p2*L2,1.,1.);
		s0 = log_integral(uu1,p1,uu2,p2,z[0],z[2],z[1],z[3]);
		s1 = log_integral(uu1,p1*L1,uu2,p2*L2,z[0],z[4],z[1],z[5]);

		qprob[nbin] = s0;
		Lambda[nbin] = s1/s0;
		xx = log(Lambda[nbin]*EN1RYD*gv.bin[nd].cnv_GR_pH);
		splint_safe(gv.bin[nd].dstems,gv.dsttmp,gv.bin[nd].dstslp,NDEMS,xx,&y,&lgErr);
		lgBoundErr = lgBoundErr || lgErr;
		qtemp[nbin] = exp(y);
		delu[nbin] = wid;
		p[nbin] = qprob[nbin]/delu[nbin];
		dPdlnT[nbin] = p[nbin]*qtemp[nbin]*uderiv(qtemp[nbin],nd);

		sum += qprob[nbin];
		RadCooling += qprob[nbin]*Lambda[nbin];

		T1 = T2;
		uu1 = uu2;
		p1 = p2;
		L1 = L2;

		++nbin;

	} while( T2 < Thi && nbin < NQGRID );

	fac = RadCooling*EN1RYD*gv.bin[nd].cnv_GR_pCM3/gv.bin[nd].GrainHeat;

	for( i=0; i < nbin; ++i )
	{
		qprob[i] /= fac;
		dPdlnT[i] /= fac;
	}

	*tol = sum/fac;
	*qnbin = nbin;
	*new_tmin = qtemp[0];
	*ErrorCode = MAX2(*ErrorCode,QH_ANALYTIC);

	/* do quality control */
	if( TolFac > STRICT )
		*ErrorCode = MAX2(*ErrorCode,QH_ANALYTIC_RELAX);

	if( lgBoundErr )
		*ErrorCode = MAX2(*ErrorCode,QH_THIGH_FAIL);

	if( fabs(sum/fac-1.) > PROB_TOL )
		*ErrorCode = MAX2(*ErrorCode,QH_CONV_FAIL);

	if( dPdlnT[0] > SAFETY*PROB_CUTOFF_LO || dPdlnT[nbin-1] > SAFETY*PROB_CUTOFF_HI )
		*ErrorCode = MAX2(*ErrorCode,QH_BOUND_FAIL);

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "     GetProbDistr_HighLimit found tol1 %.4e tol2 %.4e\n",
			 fabs(sum-1.), fabs(sum/fac-1.) );
		fprintf( ioQQQ, "    zone %ld %s nbin %ld total prob %.4e new_tmin %.4e\n",
			 nzone,gv.bin[nd].chDstLab,nbin,sum/fac,*new_tmin );
	}
	return;
}


/* calculate derivative of the enthalpy function dU/dT (aka the specific heat) at a given temperature, in Ryd/K */
STATIC double uderiv(double temp, 
		     size_t nd)
{
	enth_type ecase;
	long int i,
	  j;
	double N_C,
	  N_H;
	double deriv = 0.,
	  hok[3] = {1275., 1670., 4359.},
	  numer,
	  dnumer,
	  denom,
	  ddenom,
	  x;


	DEBUG_ENTRY( "uderiv()" );

	if( temp <= 0. ) 
	{
		fprintf( ioQQQ, " uderiv called with non-positive temperature: %.6e\n" , temp );
		cdEXIT(EXIT_FAILURE);
	}
	ASSERT( nd < gv.bin.size() );

	const double CALORIE = 4.184e7; /* this is the thermochemical calorie in erg */

	ecase = gv.which_enth[gv.bin[nd].matType];
	switch( ecase )
	{
	case ENTH_CAR:
		numer = (4.15e-22/EN1RYD)*pow(temp,3.3);
		dnumer = (3.3*4.15e-22/EN1RYD)*pow(temp,2.3);
		denom = 1. + 6.51e-03*temp + 1.5e-06*temp*temp + 8.3e-07*pow(temp,2.3);
		ddenom = 6.51e-03 + 2.*1.5e-06*temp + 2.3*8.3e-07*pow(temp,1.3);
		/* dU/dT for pah/graphitic grains in Ryd/K, derived from:
		 * >>refer	grain	physics	Guhathakurta & Draine, 1989, ApJ, 345, 230 */
		deriv = (dnumer*denom - numer*ddenom)/pow2(denom);
		break;
	case ENTH_CAR2:
		/* dU/dT for graphite grains in Ryd/K, using eq 9 of */
		/* >>refer	grain	physics	Draine B.T., and Li A., 2001, ApJ, 551, 807 */
		deriv = (DebyeDeriv(temp/863.,2) + 2.*DebyeDeriv(temp/2504.,2))*BOLTZMANN/EN1RYD;
		break;
	case ENTH_SIL:
		/* dU/dT for silicate grains (and grey grains) in Ryd/K */
		/* limits of tlim set above, 0 and DBL_MAX, so always OK */
		/* >>refer	grain	physics	Guhathakurta & Draine, 1989, ApJ, 345, 230 */
		for( j = 0; j < 4; j++ )
		{
			if( temp > tlim[j] && temp <= tlim[j+1] )
			{
				deriv = cval[j]*pow(temp,ppower[j]);
				break;
			}
		}
		break;
	case ENTH_SIL2:
		/* dU/dT for silicate grains in Ryd/K, using eq 11 of */
		/* >>refer	grain	physics	Draine B.T., and Li A., 2001, ApJ, 551, 807 */
		deriv = (2.*DebyeDeriv(temp/500.,2) + DebyeDeriv(temp/1500.,3))*BOLTZMANN/EN1RYD;
		break;
	case ENTH_PAH:
		/* dU/dT for PAH grains in Ryd/K, using eq A.4 of */
		/* >>refer	grain	physics	Dwek E., Arendt R.G., Fixsen D.J. et al., 1997, ApJ, 475, 565 */
		/* this expression is only valid upto 2000K */
		x = log10(MIN2(temp,2000.));
		deriv = exp10(-21.26+3.1688*x-0.401894*pow2(x))/EN1RYD;
		break;
	case ENTH_PAH2:
		/* dU/dT for PAH grains in Ryd/K, approximately using eq 33 of */
		/* >>refer	grain	physics	Draine B.T., and Li A., 2001, ApJ, 551, 807 */
		/* N_C and N_H should actually be nint(N_C) and nint(N_H),
		 * but this can lead to FP overflow for very large grains */
		N_C = no_atoms(nd);
		if( N_C <= 25. )
		{
			N_H = 0.5*N_C;
		}
		else if( N_C <= 100. )
		{
			N_H = 2.5*sqrt(N_C);
		}
		else
		{
			N_H = 0.25*N_C;
		}
		deriv = 0.;
		for( i=0; i < 3; i++ )
		{
			double help1 = hok[i]/temp;
			if( help1 < 300. )
			{
				double help2 = exp(help1);
				double help3 = ( help1 < 1.e-7 ) ? help1*(1.+0.5*help1) : help2-1.;
				deriv += N_H/(N_C-2.)*pow2(help1)*help2/pow2(help3)*BOLTZMANN/EN1RYD;
			}
		}
		deriv += (DebyeDeriv(temp/863.,2) + 2.*DebyeDeriv(temp/2504.,2))*BOLTZMANN/EN1RYD;
		break;
	case ENTH_SIC:
		/* c_p for alpha-SiC in cal/(mol*K), derived from Table 4 in */
		/* >>refer	grain	physics	Chekhovskoy V. Ya., 1971, J. Chem. Thermodynamics, 3, 289 */
		if( temp < 300. )
			deriv = 13.25*(3.979984e-01*DebyeDeriv(temp/747.,3) + 6.020016e-01*DebyeDeriv(temp/1647.,3));
		else
			deriv = 13.250 - 2035./temp + 288.e5/pow2(temp)*exp(-5680./temp); /* Eq. 3 */
		// now convert to Ryd/K per atom
		deriv *= CALORIE/(EN1RYD*2.*AVOGADRO);
		break;
	default:
		fprintf( ioQQQ, " uderiv called with unknown type for enthalpy function: %d\n", ecase );
		cdEXIT(EXIT_FAILURE);
	}

	/* >>chng 00 mar 23, use formula 3.1 of Guhathakurtha & Draine, by PvH */
	/* >>chng 03 jan 17, use MAX2(..,1) to prevent crash for extremely small grains, PvH */
	deriv *= MAX2(no_atoms(nd)-2.,1.);

	if( deriv <= 0. ) 
	{
		fprintf( ioQQQ, " uderiv finds non-positive derivative: %.6e, what's up?\n" , deriv );
		cdEXIT(EXIT_FAILURE);
	}
	return deriv;
}


/* calculate the enthalpy of a grain at a given temperature, in Ryd */
STATIC double ufunct(double temp, 
		     size_t nd,
		     /*@out@*/ bool *lgBoundErr)
{
	double enthalpy,
	  y;

	DEBUG_ENTRY( "ufunct()" );

	if( temp <= 0. ) 
	{
		fprintf( ioQQQ, " ufunct called with non-positive temperature: %.6e\n" , temp );
		cdEXIT(EXIT_FAILURE);
	}
	ASSERT( nd < gv.bin.size() );

	/* >>chng 02 apr 22, interpolate in DustEnth array to get enthalpy, by PvH */
	splint_safe(gv.dsttmp,gv.bin[nd].DustEnth,gv.bin[nd].EnthSlp,NDEMS,log(temp),&y,lgBoundErr);
	enthalpy = exp(y);

	ASSERT( enthalpy > 0. );
	return enthalpy;
}


/* this is the inverse of ufunct: determine grain temperature as a function of enthalpy */
STATIC double inv_ufunct(double enthalpy, 
			 size_t nd,
			 /*@out@*/ bool *lgBoundErr)
{
	double temp,
	  y;

	DEBUG_ENTRY( "inv_ufunct()" );

	if( enthalpy <= 0. ) 
	{
		fprintf( ioQQQ, " inv_ufunct called with non-positive enthalpy: %.6e\n" , enthalpy );
		cdEXIT(EXIT_FAILURE);
	}
	ASSERT( nd < gv.bin.size() );

	/* >>chng 02 apr 22, interpolate in DustEnth array to get temperature, by PvH */
	splint_safe(gv.bin[nd].DustEnth,gv.dsttmp,gv.bin[nd].EnthSlp2,NDEMS,log(enthalpy),&y,lgBoundErr);
	temp = exp(y);

	ASSERT( temp > 0. );
	return temp;
}


/* initialize interpolation arrys for grain enthalpy */
void InitEnthalpy()
{
	DEBUG_ENTRY( "InitEnthalpy()" );

	double z[8];
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		double tdust2 = GRAIN_TMIN;
		double C_V2 = uderiv(tdust2,nd);
		/* at low temps, C_V = C*T^3 -> U = C*T^4/4 = C_V*T/4 */
		gv.bin[nd].DustEnth[0] = C_V2*tdust2/4.;
		double tdust1 = tdust2;
		double C_V1 = C_V2;

		for( long i=1; i < NDEMS; i++ )
		{
			double tmid,Cmid;
			tdust2 = exp(gv.dsttmp[i]);
			C_V2 = uderiv(tdust2,nd);
			tmid = sqrt(tdust1*tdust2);
			/* this ensures accuracy for silicate enthalpy */
			for( long j=1; j < 4; j++ )
			{
				if( tdust1 < tlim[j] && tlim[j] < tdust2 )
				{
					tmid = tlim[j];
					break;
				}
			}
			Cmid = uderiv(tmid,nd);
			vlog(z,tdust1,tmid,tdust2,C_V1,Cmid,C_V2,1.,1.);
			gv.bin[nd].DustEnth[i] = gv.bin[nd].DustEnth[i-1] +
				log_integral(tdust1,C_V1,tmid,Cmid,z[0],z[3],z[1],z[4]) +
				log_integral(tmid,Cmid,tdust2,C_V2,z[1],z[4],z[2],z[5]);
			tdust1 = tdust2;
			C_V1 = C_V2;
		}
	}

	/* conversion for logarithmic interpolation */
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		vlog(gv.bin[nd].DustEnth,gv.bin[nd].DustEnth,0,NDEMS);
		/* set up coefficients for splint */
		spline(gv.dsttmp,gv.bin[nd].DustEnth,NDEMS,2e31,2e31,gv.bin[nd].EnthSlp);
		spline(gv.bin[nd].DustEnth,gv.dsttmp,NDEMS,2e31,2e31,gv.bin[nd].EnthSlp2);
	}
	return;
}


/* helper function for calculating specific heat, uses Debye theory */
STATIC double DebyeDeriv(double x,
			 long n)
{
	long i,
	  nn;
	double res;

	DEBUG_ENTRY( "DebyeDeriv()" );

	ASSERT( x > 0. );
	ASSERT( n == 2 || n == 3 );

	if( x < 0.001 )
	{
		/* for general n this is Gamma(n+2)*zeta(n+1)*powi(x,n) */
		if( n == 2 )
		{
			res = 6.*1.202056903159594*pow2(x);
		}
		else if( n == 3 )
		{
			res = 24.*1.082323233711138*pow3(x);
		}
		else
			/* added to keep lint happy - note that assert above confimred that n is 2 or 3,
			 * but lint flagged possible flow without defn of res */
			TotalInsanity();
	}
	else
	{
		nn = 4*MAX2(4,2*(long)(0.05/x));
		vector<double> xx(nn);
		vector<double> rr(nn);
		vector<double> aa(nn);
		vector<double> ww(nn);
		gauss_legendre(nn,xx,aa);
		gauss_init(nn,0.,1.,xx,aa,rr,ww);

		res = 0.;
		for( i=0; i < nn; i++ )
		{
			double help1 = rr[i]/x;
			if( help1 < 300. )
			{
				double help2 = exp(help1);
				double help3 = ( help1 < 1.e-7 ) ? help1*(1.+0.5*help1) : help2-1.;
				res += ww[i]*powi(rr[i],n+1)*help2/pow2(help3);
			}
		}
		res /= pow2(x);
	}
	return (double)n*res;
}
