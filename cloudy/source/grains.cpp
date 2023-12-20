/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*grain main routine to converge grains thermal solution */
#include "cddefines.h"
#include "atmdat.h"
#include "atmdat_adfa.h"
#include "rfield.h"
#include "hmi.h"
#include "trace.h"
#include "conv.h"
#include "ionbal.h"
#include "thermal.h"
#include "phycon.h"
#include "doppvel.h"
#include "heavy.h"
#include "ipoint.h"
#include "elementnames.h"
#include "grainvar.h"
#include "grains.h"
#include "iso.h"
#include "mole.h"
#include "dense.h"
#include "vectorize.h"
#include "parser.h"

/* the next three defines are for debugging purposes only, uncomment to activate */
/*  #define WD_TEST2 1 */
/*  #define IGNORE_GRAIN_ION_COLLISIONS 1 */
/*  #define IGNORE_THERMIONIC 1 */

/* no parentheses around PTR needed since it needs to be an lvalue */
#define FREE_CHECK(PTR) { ASSERT( PTR != NULL ); free( PTR ); PTR = NULL; }
#define FREE_SAFE(PTR) { if( PTR != NULL ) free( PTR ); PTR = NULL; }

static const long MAGIC_AUGER_DATA = 20060126L;

static const bool INCL_TUNNEL = true;
static const bool NO_TUNNEL = false;

/*================================================================================*/
/* these are used for setting up grain emissivities in InitEmissivities() */

/* NTOP is number of bins for temps between GRAIN_TMID and GRAIN_TMAX */
static const long NTOP = NDEMS/5;

/*================================================================================*/
/* these are used when iterating the grain charge in GrainCharge() */
static const double TOLER = CONSERV_TOL/10.;
static const long BRACKET_MAX = 50L;

/* >>chng 06 feb 07, increased CT_LOOP_MAX (10 -> 25), T_LOOP_MAX (30 -> 50), pah.in, PvH */

/* maximum number of tries to converge charge/temperature in GrainChargeTemp() */
static const long CT_LOOP_MAX = 25L;

/* maximum number of tries to converge grain temperature in GrainChargeTemp() */
static const long T_LOOP_MAX = 50L;

/* these will become the new tolerance levels used throughout the code */
static double HEAT_TOLER = DBL_MAX;
static double HEAT_TOLER_BIN = DBL_MAX;
static double CHRG_TOLER = DBL_MAX;
/* static double CHRG_TOLER_BIN = DBL_MAX; */

/*================================================================================*/
/* miscellaneous grain physics */

/* a_0 thru a_2 constants for calculating IP_V and EA, in cm */
static const double AC0 = 3.e-9;
static const double AC1G = 4.e-8;
static const double AC2G = 7.e-8;

/* constants needed to calculate energy distribution of secondary electrons */
static const double ETILDE = 2.*SQRT2/EVRYD; /* sqrt(8) eV */
static const double INV_ETILDE = 1./ETILDE;

/* constant for thermionic emissions, 7.501e20 e/cm^2/s/K^2 */
static const double THERMCONST = PI4*ELECTRON_MASS*POW2(BOLTZMANN)/POW3(HPLANCK);

/* sticking probabilities */
static const double STICK_ELEC = 0.5;
static const double STICK_ION = 1.0;

/** evaluate e^2/a, the potential of one electron */
inline double one_elec(long nd)
{
	return ELEM_CHARGE/EVRYD/gv.bin[nd].Capacity;
}

/** convert grain potential in Ryd to charge in electrons */
inline double pot2chrg(double x,
		       long nd)
{
	return x/one_elec(nd) - 1.;
}

/** convert grain charge in electrons into potential in Ryd */
inline double chrg2pot(double x,
		       long nd)
{
	return (x+1.)*one_elec(nd);
}

/** mean pathlength travelled by electrons inside the grain, in cm (Eq. 11 of WDB06) */
inline double elec_esc_length(double e, // energy of electron in Ryd
			      long nd)
{
	// calculate escape length in cm
	if( e <= gv.bin[nd].le_thres )
		return 1.e-7;
	else
		return 3.e-6*gv.bin[nd].eec*powpq(e*EVRYD*1.e-3,3,2);
}
	
/* read data for electron energy spectrum of Auger electrons */
STATIC void ReadAugerData();
/* initialize the Auger data for grain bin nd between index ipBegin <= i < ipEnd */
STATIC void InitBinAugerData(size_t,long,long);
/* initialize grain emissivities */
STATIC void InitEmissivities();
/* PlanckIntegral compute total radiative cooling due to large grains */
STATIC double PlanckIntegral(double,size_t,long);
/* invalidate charge dependent data from previous iteration */
STATIC void NewChargeData(long);
/* GrnStdDpth returns the grain abundance as a function of depth into cloud */
STATIC double GrnStdDpth(long);
/* iterate grain charge and temperature */
STATIC void GrainChargeTemp();
/* GrainCharge compute grains charge */
STATIC void GrainCharge(size_t,/*@out@*/double*);
/* grain electron recombination rates for single charge state */
STATIC double GrainElecRecomb1(size_t,long,/*@out@*/double*,/*@out@*/double*);
/* grain electron emission rates for single charge state */
STATIC double GrainElecEmis1(size_t,long,/*@out@*/double*,/*@out@*/double*,/*@out@*/double*,/*@out@*/double*);
/* correction factors for grain charge screening (including image potential
 * to correct for polarization of the grain as charged particle approaches). */
STATIC void GrainScreen(long,size_t,long,double*,double*);
/* helper function for GrainScreen */
STATIC double ThetaNu(double);
/* update items that depend on grain potential */
STATIC void UpdatePot(size_t,long,long,/*@out@*/double[],/*@out@*/double[]);
/* calculate charge state populations */
STATIC void GetFracPop(size_t,long,/*@in@*/double[],/*@in@*/double[],/*@out@*/long*);
/* this routine updates all quantities that depend only on grain charge and radius */
STATIC void UpdatePot1(size_t,long,long,long);
/* this routine updates all quantities that depend on grain charge, radius and temperature */
STATIC void UpdatePot2(size_t,long);
/* Helper function to calculate primary and secondary yields and the average electron energy at infinity */
STATIC void Yfunc(long,long,const realnum[],const realnum[],const realnum[],double,const double[],const double[],
		  /*@out@*/realnum[],/*@out@*/realnum[],/*@out@*/realnum[],/*@out@*/realnum[],long,long);
STATIC void Yfunc(long,long,const realnum[],const realnum[],double,double,double,
		  /*@out@*/realnum[],/*@out@*/realnum[],/*@out@*/realnum[],/*@out@*/realnum[],long,long);
/* This calculates the y0 function for band electrons (Sect. 4.1.3/4.1.4 of WDB06) */
STATIC void y0b(size_t,long,/*@out@*/realnum[],long,long);
/* This calculates the y0 function for band electrons (Eq. 16 of WD01) */
STATIC void y0b01(size_t,long,/*@out@*/realnum[],long,long);
/* This calculates the y0 function for primary/secondary and Auger electrons (Eq. 9 of WDB06) */
STATIC double y0psa(size_t,long,long,double);
/* This calculates the y1 function for primary/secondary and Auger electrons (Eq. 6 of WDB06) */
STATIC double y1psa(size_t,long,double);
/* This calculates the y2 function for primary and Auger electrons (Eq. 8 of WDB06) */
inline void y2pa(double,const double[],long,/*@out@*/realnum[],/*@out@*/realnum[],long,long);
/* This calculates the y2 function for secondary electrons (Eq. 20-21 of WDB06) */
inline void y2s(double,const double[],long,const realnum[],/*@out@*/realnum[],/*@out@*/realnum[],long,long);
/* determine charge Z0 ion recombines to upon impact on grain */
STATIC void UpdateRecomZ0(size_t,long);
/* helper routine for UpdatePot */
STATIC void GetPotValues(size_t,long,/*@out@*/double*,/*@out@*/double*,/*@out@*/double*,
			 /*@out@*/double*,/*@out@*/double*,/*@out@*/double*,bool);
/* given grain nd in charge state nz, and incoming ion (nelem,ion),
 * detemine outgoing ion (nelem,Z0) and chemical energy ChEn released
 * ChemEn is net contribution of ion recombination to grain heating */
STATIC void GrainIonColl(size_t,long,long,long,const double[],const double[],/*@out@*/long*,
			 /*@out@*/realnum*,/*@out@*/realnum*);
/* initialize ion recombination rates on grain species nd */
STATIC void GrainChrgTransferRates(long);
/* this routine updates all grain quantities that depend on radius, except gv.dstab and gv.dstsc */
STATIC void GrainUpdateRadius1();
/* this routine adds all the grain opacities in gv.dstab and gv.dstsc */
STATIC void GrainUpdateRadius2();
/* GrainTemperature computes grains temperature, and gas cooling */
STATIC void GrainTemperature(size_t,/*@out@*/realnum*,/*@out@*/double*,/*@out@*/double*,
			     /*@out@*/double*);
/* helper routine for initializing quantities related to the photo-electric effect */
STATIC void PE_init(size_t,long,long,/*@out@*/double*,/*@out@*/double*,/*@out@*/double*,
		    /*@out@*/double*,/*@out@*/double*,/*@out@*/double*,/*@out@*/double*);
/* GrainCollHeating computes grains collisional heating cooling */
STATIC void GrainCollHeating(size_t,/*@out@*/realnum*,/*@out@*/realnum*);
/* GrnVryDpth user supplied function for the grain abundance as a function of depth into cloud */
STATIC double GrnVryDpth(size_t);


/* this routine is called by IterStart(), so anything that needs to be reset before each
 * iteration starts should be put here; typically variables that are integrated over radius */
void GrainStartIter()
{
	DEBUG_ENTRY( "GrainStartIter()" );

	if( gv.lgDustOn() && gv.lgGrainPhysicsOn )
	{
		gv.lgNegGrnDrg = false;
		gv.TotalDustHeat = 0.;
		gv.GrnElecDonateMax = 0.;
		gv.GrnElecHoldMax = 0.;
		gv.dphmax = 0.f;
		gv.dclmax = 0.f;

		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			gv.bin[nd].ZloSave = gv.bin[nd].chrg(0).DustZ;
			gv.bin[nd].qtmin = ( gv.bin[nd].qtmin_zone1 > 0. ) ?
				gv.bin[nd].qtmin_zone1 : DBL_MAX;
			gv.bin[nd].avdust = 0.;
			gv.bin[nd].avdpot = 0.;
			gv.bin[nd].avdft = 0.;
			gv.bin[nd].avDGRatio = 0.;
			gv.bin[nd].TeGrainMax = -1.f;
			gv.bin[nd].lgEverQHeat = false;
			gv.bin[nd].QHeatFailures = 0L;
			gv.bin[nd].lgQHTooWide = false;
			gv.bin[nd].lgPAHsInIonizedRegion = false;
			gv.bin[nd].nChrgOrg = gv.bin[nd].nChrg;
		}
	}
	return;
}


/* this routine is called by IterRestart(), so anything that needs to be
 * reset or saved after an iteration is finished should be put here */
void GrainRestartIter()
{
	DEBUG_ENTRY( "GrainRestartIter()" );

	if( gv.lgDustOn() && gv.lgGrainPhysicsOn )
	{
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			gv.bin[nd].lgIterStart = true;
			gv.bin[nd].nChrg = gv.bin[nd].nChrgOrg;
		}
	}
	return;
}


/* this routine is called by ParseSet() */
void SetNChrgStates(long nChrg)
{
	DEBUG_ENTRY( "SetNChrgStates()" );

	ASSERT( nChrg >= 2 && nChrg <= NCHU );
	gv.nChrgRequested = nChrg;
	return;
}


/*GrainsInit, called one time by opacitycreateall at initialization of calculation, 
 * called after commands have been parsed,
 * not after every iteration or every model */
void GrainsInit()
{
	long int i,
	  nelem;
	unsigned int ns;

	DEBUG_ENTRY( "GrainsInit()" );

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, " GrainsInit called.\n" );
	}

	gv.dstab.resize( rfield.nflux_with_check );
	gv.dstsc.resize( rfield.nflux_with_check );
	gv.GrainEmission.resize( rfield.nflux_with_check );
	gv.GraphiteEmission.resize( rfield.nflux_with_check );
	gv.SilicateEmission.resize( rfield.nflux_with_check );

	/* >>chng 02 jan 15, initialize to zero in case grains are not used, needed in IonIron(), PvH */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		gv.elmSumAbund[nelem] = 0.f;
	}

	for( i=0; i < rfield.nflux_with_check; i++ )
	{
		gv.dstab[i] = 0.;
		gv.dstsc[i] = 0.;
		/* >>chng 01 sep 12, moved next three initializations from GrainZero(), PvH */
		gv.GrainEmission[i] = 0.;
		gv.SilicateEmission[i] = 0.;
		gv.GraphiteEmission[i] = 0.;
	}

	if( !gv.lgDustOn() )
	{
		/* grains are not on, set all heating/cooling agents to zero */
		gv.GrainHeatInc = 0.;
		gv.GrainHeatDif = 0.;
		gv.GrainHeatLya = 0.;
		gv.GrainHeatCollSum = 0.;
		gv.GrainHeatSum = 0.;
		gv.GasCoolColl = 0.;
		thermal.setHeating(0,13,0.);
		thermal.setHeating(0,14,0.);
		thermal.setHeating(0,25,0.);

		if( trace.lgTrace && trace.lgDustBug )
		{
			fprintf( ioQQQ, " GrainsInit exits.\n" );
		}
		return;
	}

#ifdef WD_TEST2
	gv.lgWD01 = true;
#endif

	HEAT_TOLER = conv.HeatCoolRelErrorAllowed / 3.;
	HEAT_TOLER_BIN = HEAT_TOLER / sqrt((double)gv.bin.size());
	CHRG_TOLER = conv.EdenErrorAllowed / 3.;
	/* CHRG_TOLER_BIN = CHRG_TOLER / sqrt(gv.bin.size()); */

	ReadAugerData();

	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		double help,atoms,p_rad,ThresInf,ThresInfVal,Emin,d[5];
		long low1,low2,low3,lowm;

		/* sanity checks */
		ASSERT( gv.bin[nd].nChrg >= 2 && gv.bin[nd].nChrg <= NCHU );

		if( gv.bin[nd].DustWorkFcn < rfield.emm() || gv.bin[nd].DustWorkFcn > rfield.egamry() )
		{
			fprintf( ioQQQ, " Grain work function for %s has insane value: %.4e\n",
				 gv.bin[nd].chDstLab,gv.bin[nd].DustWorkFcn );
			cdEXIT(EXIT_FAILURE);
		}

		/* this is QHEAT ALL command */
		if( gv.lgQHeatAll )
		{
			gv.bin[nd].lgQHeat = true;
		}

		/* this is NO GRAIN QHEAT command, always takes precedence */
		if( !gv.lgQHeatOn ) 
		{
			gv.bin[nd].lgQHeat = false;
		}

		/* >>chng 04 jun 01, disable quantum heating when constant grain temperature is used, PvH */
		if( thermal.ConstGrainTemp > 0. )
		{
			gv.bin[nd].lgQHeat = false;
		}

		if (ENABLE_QUANTUM_HEATING)
		{
			gv.bin[nd].lgQHTooWide = false;
			gv.bin[nd].qtmin = DBL_MAX;
		}

		gv.bin[nd].lgIterStart = true;

		if( gv.bin[nd].nDustFunc>DF_STANDARD || gv.bin[nd].matType == MAT_PAH || 
			gv.bin[nd].matType == MAT_PAH2 )
			gv.lgAnyDustVary = true;

		/* grain abundance may depend on radius,
		 * invalidate for now; GrainUpdateRadius1() will set correct value */
		gv.bin[nd].dstAbund = -FLT_MAX;

		gv.bin[nd].GrnDpth = 1.f;

		gv.bin[nd].qtmin_zone1 = -1.;

		/* this is threshold in Ryd above which to use X-ray prescription for electron escape length */
		gv.bin[nd].le_thres = gv.lgWD01 ? FLT_MAX :
			(realnum)(powpq(pow((double)gv.bin[nd].dustp[0],0.85)/30.,2,3)*1.e3/EVRYD);

		/* >>chng 00 jun 19, this value is absolute lower limit for the grain
		 * potential, electrons cannot be bound for lower values..., PvH */
		zmin_type zcase = gv.which_zmin[gv.bin[nd].matType];
		switch( zcase )
		{
		case ZMIN_CAR:
			// this is Eq. 23a + 24 of WD01
			help = gv.bin[nd].AvRadius*1.e7;
			help = ceil(-(1.2*POW2(help)+3.9*help+0.2)/1.44);
			break;
		case ZMIN_SIL:
			// this is Eq. 23b + 24 of WD01
			help = gv.bin[nd].AvRadius*1.e7;
			help = ceil(-(0.7*POW2(help)+2.5*help+0.8)/1.44);
			break;
		default:
			fprintf( ioQQQ, " GrainsInit detected unknown Zmin type: %d\n" , zcase );
			cdEXIT(EXIT_FAILURE);
		}

		/* this is to assure that gv.bin[nd].LowestZg > LONG_MIN */
		ASSERT( help > (double)(LONG_MIN+1) );
		low1 = nint(help);

		/* >>chng 01 apr 20, iterate to get LowestZg such that the exponent in the thermionic
		 * rate never becomes positive; the value can be derived by equating ThresInf >= 0;
		 * the new expression for Emin (see GetPotValues) cannot be inverted analytically,
		 * hence it is necessary to iterate for LowestZg. this also automatically assures that
		 * the expressions for ThresInf and LowestZg are consistent with each other, PvH */
		low2 = low1;
		GetPotValues(nd,low2,&ThresInf,&d[0],&d[1],&d[2],&d[3],&d[4],INCL_TUNNEL);
		if( ThresInf < 0. )
		{
			low3 = 0;
			/* do a bisection search for the lowest charge such that
			 * ThresInf >= 0, the end result will eventually be in low3 */
			while( low3-low2 > 1 )
			{
				lowm = (low2+low3)/2;
				GetPotValues(nd,lowm,&ThresInf,&d[0],&d[1],&d[2],&d[3],&d[4],INCL_TUNNEL);
				if( ThresInf < 0. )
					low2 = lowm;
				else
					low3 = lowm;
			}
			low2 = low3;
		}

		/* the first term implements the minimum charge due to autoionization
		 * the second term assures that the exponent in the thermionic rate never
		 * becomes positive; the expression was derived by equating ThresInf >= 0 */
		gv.bin[nd].LowestZg = MAX2(low1,low2);
		gv.bin[nd].ZloSave = LONG_MIN;
		ASSERT( gv.bin[nd].sd.size() == 0 );
		
		ns = 0;
		gv.bin[nd].sd.emplace_back( ShellData() );

		/* this is data for valence band */
		gv.bin[nd].sd[ns].nelem = -1;
		gv.bin[nd].sd[ns].ns = -1;
		gv.bin[nd].sd[ns].ionPot = gv.bin[nd].DustWorkFcn;

		/* now add data for inner shell photoionization */
		for( nelem=ipLITHIUM; nelem < LIMELM && !gv.lgWD01; nelem++ )
		{
			if( gv.bin[nd].elmAbund[nelem] > 0. )
			{
				if( gv.AugerData[nelem].nSubShell == 0 )
				{
					fprintf( ioQQQ, " Grain Auger data are missing for element %s\n",
						 elementnames.chElementName[nelem] );
					fprintf( ioQQQ, " Please include the NO GRAIN X-RAY TREATMENT command "
						 "to disable the Auger treatment in grains.\n" );
					cdEXIT(EXIT_FAILURE);
				}

				for( unsigned int j=0; j < gv.AugerData[nelem].nSubShell; j++ )
				{
					++ns;
					gv.bin[nd].sd.emplace_back( ShellData() );

					gv.bin[nd].sd[ns].nelem = nelem;
					gv.bin[nd].sd[ns].ns = j;
					gv.bin[nd].sd[ns].ionPot = gv.AugerData[nelem].IonThres[j];
				}
			}
		}

		GetPotValues(nd,gv.bin[nd].LowestZg,&d[0],&ThresInfVal,&d[1],&d[2],&d[3],&Emin,INCL_TUNNEL);

		for( ns=0; ns < gv.bin[nd].sd.size(); ns++ )
		{
			long ipLo;
			double Ethres = ( ns == 0 ) ? ThresInfVal : gv.bin[nd].sd[ns].ionPot;
			ShellData *sptr = &gv.bin[nd].sd[ns];

			sptr->ipLo = rfield.ithreshC( Ethres );

			ipLo = sptr->ipLo;
			// allow room for adjustment of rfield.nPositive later on
			long len = rfield.nflux_with_check - ipLo;

			sptr->p.reserve( len );
			sptr->p.alloc( ipLo, rfield.nPositive );

			sptr->y01.reserve( len );
			sptr->y01.alloc( ipLo, rfield.nPositive );

			/* there are no Auger electrons from the band structure */
			if( ns > 0 )
			{
				sptr->nData = gv.AugerData[sptr->nelem].nData[sptr->ns];
				sptr->AvNr.resize( sptr->nData );
				sptr->Ener.resize( sptr->nData );
				sptr->y01A.resize( sptr->nData );

				for( long n=0; n < sptr->nData; n++ )
				{
					sptr->AvNr[n] = gv.AugerData[sptr->nelem].AvNumber[sptr->ns][n];
					sptr->Ener[n] = gv.AugerData[sptr->nelem].Energy[sptr->ns][n];

					sptr->y01A[n].reserve( len );
					sptr->y01A[n].alloc( ipLo, rfield.nPositive );
				}
			}
		}

		gv.bin[nd].y0b06.resize( rfield.nflux_with_check );

		InitBinAugerData( nd, 0, rfield.nPositive );

		gv.bin[nd].nfill = rfield.nPositive;

		/* >>chng 00 jul 13, new sticking probability for electrons */
		/* the second term is chance that electron passes through grain,
		 * 1-p_rad is chance that electron is ejected before grain settles
		 * see discussion in 
		 * >>refer	grain	physics	Weingartner & Draine, 2001, ApJS, 134, 263 */
		/** \todo xray - StickElec depends on Te ???? use elec_esc_length(1.5*kTe,nd) ???? */
		gv.bin[nd].StickElecPos = STICK_ELEC*(1. - exp(-gv.bin[nd].AvRadius/elec_esc_length(0.,nd)));
		atoms = gv.bin[nd].AvVol*gv.bin[nd].dustp[0]/ATOMIC_MASS_UNIT/gv.bin[nd].atomWeight;
		p_rad = 1./(1.+exp(20.-atoms));
		gv.bin[nd].StickElecNeg = gv.bin[nd].StickElecPos*p_rad;

		/* >>chng 02 feb 15, these quantities depend on radius and are normally set
		 * in GrainUpdateRadius1(), however, it is necessary to initialize them here
		 * as well so that they are valid the first time hmole is called. */
		gv.bin[nd].GrnDpth = (realnum)GrnStdDpth(nd);
		gv.bin[nd].dstAbund = (realnum)(gv.bin[nd].dstfactor*gv.GrainMetal*gv.bin[nd].GrnDpth);
		ASSERT( gv.bin[nd].dstAbund > 0.f );
		/* grain unit conversion, <unit>/H (default depl) -> <unit>/cm^3 (actual depl) */
		gv.bin[nd].cnv_H_pCM3 = dense.gas_phase[ipHYDROGEN]*gv.bin[nd].dstAbund;
		gv.bin[nd].cnv_CM3_pH = 1./gv.bin[nd].cnv_H_pCM3;
		/* grain unit conversion, <unit>/cm^3 (actual depl) -> <unit>/grain */
		gv.bin[nd].cnv_CM3_pGR = gv.bin[nd].cnv_H_pGR/gv.bin[nd].cnv_H_pCM3;
		gv.bin[nd].cnv_GR_pCM3 = 1./gv.bin[nd].cnv_CM3_pGR;
	}

	/* >>chng 02 dec 19, these quantities depend on radius and are normally set
	 * in GrainUpdateRadius1(), however, it is necessary to initialize them here
	 * as well so that they are valid for the initial printout in Cloudy, PvH */
	/* calculate the summed grain abundances, these are valid at the inner radius;
	 * these numbers depend on radius and are updated in GrainUpdateRadius1() */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		gv.elmSumAbund[nelem] = 0.f;
	}

	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		for( nelem=0; nelem < LIMELM; nelem++ )
		{
			gv.elmSumAbund[nelem] += gv.bin[nd].elmAbund[nelem]*(realnum)gv.bin[nd].cnv_H_pCM3;
		}
	}

	gv.nzone = -1;
	gv.GrnRecomTe = -1.;

	/* >>chng 01 nov 21, total grain opacities depend on charge and therefore on radius,
	 *                   invalidate for now; GrainUpdateRadius2() will set correct values */
	for( i=0; i < rfield.nflux_with_check; i++ )
	{
		/* these are total absorption and scattering cross sections,
		 * the latter should contain the asymmetry factor (1-g) */
		gv.dstab[i] = -DBL_MAX;
		gv.dstsc[i] = -DBL_MAX;
	}

	InitEmissivities();
	InitEnthalpy();

	if( trace.lgDustBug && trace.lgTrace )
	{
		fprintf( ioQQQ, "     There are %ld grain types turned on.\n", (unsigned long)gv.bin.size() );

		fprintf( ioQQQ, "     grain depletion factors, dstfactor*GrainMetal=" );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, "%10.2e", gv.bin[nd].dstfactor*gv.GrainMetal );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "     nChrg =" );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, " %ld", gv.bin[nd].nChrg );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "     lowest charge (e) =" );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, " %ld", gv.bin[nd].LowestZg );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "     nDustFunc flag for user requested custom depth dependence:" );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, "%2i", gv.bin[nd].nDustFunc );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "     Quantum heating flag:" );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, "%2c", TorF(gv.bin[nd].lgQHeat) );
		}
		fprintf( ioQQQ, "\n" );

		/* >>chng 01 nov 21, removed total abs and sct cross sections, they are invalid */
		fprintf( ioQQQ, "     NU(Ryd), Abs cross sec per proton\n" );

		fprintf( ioQQQ, "    Ryd   " );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, " %-12.12s", gv.bin[nd].chDstLab );
		}
		fprintf( ioQQQ, "\n" );

		for( i=0; i < rfield.nflux_with_check; i += 40 )
		{
			fprintf( ioQQQ, "%10.2e", rfield.anu(i) );
			for( size_t nd=0; nd < gv.bin.size(); nd++ )
			{
				fprintf( ioQQQ, " %10.2e  ", gv.bin[nd].dstab1[i] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, "     NU(Ryd), Sct cross sec per proton\n" );

		fprintf( ioQQQ, "    Ryd   " );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, " %-12.12s", gv.bin[nd].chDstLab );
		}
		fprintf( ioQQQ, "\n" );

		for( i=0; i < rfield.nflux_with_check; i += 40 )
		{
			fprintf( ioQQQ, "%10.2e", rfield.anu(i) );
			for( size_t nd=0; nd < gv.bin.size(); nd++ )
			{
				fprintf( ioQQQ, " %10.2e  ", gv.bin[nd].pure_sc1[i] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, "     NU(Ryd), Q abs\n" );

		fprintf( ioQQQ, "    Ryd   " );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, " %-12.12s", gv.bin[nd].chDstLab );
		}
		fprintf( ioQQQ, "\n" );

		for( i=0; i < rfield.nflux_with_check; i += 40 )
		{
			fprintf( ioQQQ, "%10.2e", rfield.anu(i) );
			for( size_t nd=0; nd < gv.bin.size(); nd++ )
			{
				fprintf( ioQQQ, " %10.2e  ", gv.bin[nd].dstab1[i]*4./gv.bin[nd].IntArea );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, "     NU(Ryd), Q sct\n" );

		fprintf( ioQQQ, "    Ryd   " );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, " %-12.12s", gv.bin[nd].chDstLab );
		}
		fprintf( ioQQQ, "\n" );

		for( i=0; i < rfield.nflux_with_check; i += 40 )
		{
			fprintf( ioQQQ, "%10.2e", rfield.anu(i) );
			for( size_t nd=0; nd < gv.bin.size(); nd++ )
			{
				fprintf( ioQQQ, " %10.2e  ", gv.bin[nd].pure_sc1[i]*4./gv.bin[nd].IntArea );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, "     NU(Ryd), asymmetry factor\n" );

		fprintf( ioQQQ, "    Ryd   " );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, " %-12.12s", gv.bin[nd].chDstLab );
		}
		fprintf( ioQQQ, "\n" );

		for( i=0; i < rfield.nflux_with_check; i += 40 )
		{
			fprintf( ioQQQ, "%10.2e", rfield.anu(i) );
			for( size_t nd=0; nd < gv.bin.size(); nd++ )
			{
				fprintf( ioQQQ, " %10.2e  ", gv.bin[nd].asym[i] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, " GrainsInit exits.\n" );
	}
	return;
}

/* read data for electron energy spectrum of Auger electrons */
STATIC void ReadAugerData()
{
	DEBUG_ENTRY( "ReadAugerData()" );

	DataParser d("auger_spectrum.dat", ES_NONE);
	d.getline();
	d.checkMagic(MAGIC_AUGER_DATA);

	while( true )
	{
		d.getline();
		long nelem;
		d.getToken(nelem);

		if( nelem < 0 )
			break;

		ASSERT( nelem < LIMELM );

		AEInfo *ptr = &gv.AugerData[nelem];

		d.getline();
		d.getToken(ptr->nSubShell);
		ASSERT( ptr->nSubShell > 0 );

		ptr->nData.resize( ptr->nSubShell );
		ptr->IonThres.resize( ptr->nSubShell );
		ptr->Energy.resize( ptr->nSubShell );
		ptr->AvNumber.resize( ptr->nSubShell );

		for( unsigned int ns=0; ns < ptr->nSubShell; ns++ )
		{
			unsigned int ss;

			d.getline();
			d.getToken(ss);
			ASSERT( ns == ss );

			d.getline();
			d.getToken(ptr->IonThres[ns]);
			ptr->IonThres[ns] /= EVRYD;

			d.getline();
			d.getToken(ptr->nData[ns]);
			ASSERT( ptr->nData[ns] > 0 );

			ptr->Energy[ns].resize( ptr->nData[ns] );
			ptr->AvNumber[ns].resize( ptr->nData[ns] );

			for( unsigned int n=0; n < ptr->nData[ns]; n++ )
			{
				d.getline();
				d.getToken(ptr->Energy[ns][n]);
				d.getToken(ptr->AvNumber[ns][n]);
				ptr->Energy[ns][n] /= EVRYD;
				ASSERT( ptr->Energy[ns][n] < ptr->IonThres[ns] );
			}
		}
	}
}

/** initialize the Auger data for grain bin nd between index ipBegin <= i < ipEnd */
STATIC void InitBinAugerData(size_t nd,
			     long ipBegin,
			     long ipEnd)
{
	DEBUG_ENTRY( "InitBinAugerData()" );

	long i, ipLo, nelem;
	unsigned int ns;

	flex_arr<realnum> temp( ipBegin, ipEnd );
	temp.zero();

	/* this converts gv.bin[nd].elmAbund[nelem] to particle density inside the grain */
	double norm = gv.bin[nd].cnv_H_pGR/gv.bin[nd].AvVol;

	/* this loop calculates the probability that photoionization occurs in a given shell */
	for( ns=0; ns < gv.bin[nd].sd.size(); ns++ )
	{
		ipLo = max( gv.bin[nd].sd[ns].ipLo, ipBegin );

		gv.bin[nd].sd[ns].p.realloc( ipEnd );

		/** \todo xray - Compton recoil still needs to be added here */

		for( i=ipLo; i < ipEnd; i++ )
		{
			long nel,nsh;
			double phot_ev,cs;

			phot_ev = rfield.anu(i)*EVRYD;

			if( ns == 0 )
			{
				/* this is the valence band, defined as the sum of any
				 * subshell not treated explicitly as an inner shell below */
				gv.bin[nd].sd[ns].p[i] = 0.;

				for( nelem=ipHYDROGEN; nelem < LIMELM && !gv.lgWD01; nelem++ )
				{
					if( gv.bin[nd].elmAbund[nelem] == 0. )
						continue;

					long nshmax = Heavy.nsShells[nelem][0];

					for( nsh = gv.AugerData[nelem].nSubShell; nsh < nshmax; nsh++ )
					{
						nel = nelem+1;
						cs = t_ADfA::Inst().phfit(nelem+1,nel,nsh+1,phot_ev);
						gv.bin[nd].sd[ns].p[i] +=
							(realnum)(norm*gv.bin[nd].elmAbund[nelem]*cs*1e-18);
					}
				}

				temp[i] += gv.bin[nd].sd[ns].p[i];
			}
			else
			{
				/* this is photoionization from inner shells */
				nelem = gv.bin[nd].sd[ns].nelem;
				nel = nelem+1;
				nsh = gv.bin[nd].sd[ns].ns;
				cs = t_ADfA::Inst().phfit(nelem+1,nel,nsh+1,phot_ev);
				gv.bin[nd].sd[ns].p[i] =
					(realnum)(norm*gv.bin[nd].elmAbund[nelem]*cs*1e-18);
				temp[i] += gv.bin[nd].sd[ns].p[i];
			}
		}
	}

	for( i=ipBegin; i < ipEnd && !gv.lgWD01; i++ )
	{
		/* this is Eq. 10 of WDB06 */
		if( rfield.anu(i) > 20./EVRYD )
			gv.bin[nd].inv_att_len[i] = temp[i];
	}

	for( ns=0; ns < gv.bin[nd].sd.size(); ns++ )
	{
		ipLo = max( gv.bin[nd].sd[ns].ipLo, ipBegin );
		/* renormalize so that sum of probabilities is 1 */
		for( i=ipLo; i < ipEnd; i++ )
		{
			if( temp[i] > 0. )
				gv.bin[nd].sd[ns].p[i] /= temp[i];
			else
				gv.bin[nd].sd[ns].p[i] = ( ns == 0 ) ? 1.f : 0.f;
		}
	}

	temp.clear();

	for( ns=0; ns < gv.bin[nd].sd.size(); ns++ )
	{
		long n;
		ShellData *sptr = &gv.bin[nd].sd[ns];

		ipLo = max( sptr->ipLo, ipBegin );

		/* initialize the yield for primary electrons */
		sptr->y01.realloc( ipEnd );

		for( i=ipLo; i < ipEnd; i++ )
		{
			double elec_en,yzero,yone;

			elec_en = MAX2(rfield.anu(i) - sptr->ionPot,0.);
			yzero = y0psa( nd, ns, i, elec_en );

			/* this is size-dependent geometrical yield enhancement
			 * defined in Weingartner & Draine, 2001; modified in WDB06 */
			yone = y1psa( nd, i, elec_en );

			if( ns == 0 )
			{
				gv.bin[nd].y0b06[i] = (realnum)yzero;
				sptr->y01[i] = (realnum)yone;
			}
			else
			{
				sptr->y01[i] = (realnum)(yzero*yone);
			}
		}

		/* there are no Auger electrons from the band structure */
		if( ns > 0 )
		{
			/* initialize the yield for Auger electrons */
			for( n=0; n < sptr->nData; n++ )
			{
				sptr->y01A[n].realloc( ipEnd );

				for( i=ipLo; i < ipEnd; i++ )
				{
					double yzero = sptr->AvNr[n] * y0psa( nd, ns, i, sptr->Ener[n] );

					/* this is size-dependent geometrical yield enhancement
					 * defined in Weingartner & Draine, 2001; modified in WDB06 */
					double yone = y1psa( nd, i, sptr->Ener[n] );

					sptr->y01A[n][i] = (realnum)(yzero*yone);
				}
			}
		}
	}
}

STATIC void InitEmissivities()
{
	double fac,
	  fac2,
	  tdust;
	long int i;

	DEBUG_ENTRY( "InitEmissivities()" );

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "  InitEmissivities starts\n" );
		fprintf( ioQQQ, "    ND    Tdust       Emis       BB Check   4pi*a^2*<Q>\n" );
	}

	ASSERT( NTOP >= 2 && NDEMS > 2*NTOP );
	fac = exp(log(GRAIN_TMID/GRAIN_TMIN)/(double)(NDEMS-NTOP));
	tdust = GRAIN_TMIN;
	for( i=0; i < NDEMS-NTOP; i++ )
	{
		gv.dsttmp[i] = log(tdust);
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			gv.bin[nd].dstems[i] = log(PlanckIntegral(tdust,nd,i));
		}
		tdust *= fac;
	}

	/* temperatures above GRAIN_TMID are unrealistic -> make grid gradually coarser */
	fac2 = exp(log(GRAIN_TMAX/GRAIN_TMID/powi(fac,NTOP-1))/(double)((NTOP-1)*NTOP/2));
	for( i=NDEMS-NTOP; i < NDEMS; i++ )
	{
		gv.dsttmp[i] = log(tdust);
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			gv.bin[nd].dstems[i] = log(PlanckIntegral(tdust,nd,i));
		}
		fac *= fac2;
		tdust *= fac;
	}

#ifndef NDEBUG
	/* sanity checks */
	double mul = 1.;
	ASSERT( fabs(gv.dsttmp[0] - log(GRAIN_TMIN)) < 10.*mul*DBL_EPSILON );
	mul = sqrt((double)(NDEMS-NTOP));
	ASSERT( fabs(gv.dsttmp[NDEMS-NTOP] - log(GRAIN_TMID)) < 10.*mul*DBL_EPSILON );
	mul = (double)NTOP + sqrt((double)NDEMS);
	ASSERT( fabs(gv.dsttmp[NDEMS-1] - log(GRAIN_TMAX)) < 10.*mul*DBL_EPSILON );
#endif

	/* now find slopes form spline fit */
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		/* set up coefficients for spline */
		spline(gv.bin[nd].dstems,gv.dsttmp,NDEMS,2e31,2e31,gv.bin[nd].dstslp);
		spline(gv.dsttmp,gv.bin[nd].dstems,NDEMS,2e31,2e31,gv.bin[nd].dstslp2);
	}

#	if 0
	/* test the dstems interpolation */
	nd = nint(fudge(0));
	ASSERT( nd >= 0 && nd < gv.bin.size() );
	for( i=0; i < 2000; i++ )
	{
		double x,y,z;
		z = exp10(-40. + (double)i/50.);
		splint(gv.bin[nd].dstems,gv.dsttmp,gv.bin[nd].dstslp,NDEMS,log(z),&y);
		if( exp(y) > GRAIN_TMIN && exp(y) < GRAIN_TMAX )
		{
			x = PlanckIntegral(exp(y),nd,3);
			printf(" input %.6e temp %.6e output %.6e rel. diff. %.6e\n",z,exp(y),x,(x-z)/z);
		}
	}
	cdEXIT(EXIT_SUCCESS);
#	endif
	return;
}


/* PlanckIntegral compute total radiative cooling due to grains */
STATIC double PlanckIntegral(double tdust, 
			     size_t nd, 
			     long int ip)
{
	DEBUG_ENTRY( "PlanckIntegral()" );

	/******************************************************************
	 *
	 * >>>chng 99 mar 12, this sub rewritten following Peter van Hoof
	 * comments.  Original coding was in single precision, and for
	 * very low temperature the exponential was set to zero.  As 
	 * a result Q was far too large for grain temperatures below 10K
	 *
	 ******************************************************************/

	/* Boltzmann factors for Planck integration */
	double TDustRyg = TE1RYD/tdust;
	double x = 0.999*log(DBL_MAX);
	long mini = 0, maxi = rfield.nflux_with_check;
	avx_ptr<double> arg(rfield.nflux_with_check), expval(rfield.nflux_with_check);
	for( long i=0; i < rfield.nflux_with_check; i++ )
	{
		/* this is hnu/kT for grain at this temp and photon energy */
		arg[i] = TDustRyg*rfield.anu(i);
		if( arg[i] < 0.9999e-5 )
			++mini;
		if( arg[i] > x+0.0001 )
		{
			maxi = i;
			break;
		}
	}
	vexp( arg.ptr0(), expval.ptr0(), mini, maxi );
		
	double integral1 = 0.;  /* integral(Planck) */
	double integral2 = 0.;  /* integral(Planck*abs_cs) */

	for( long i=0; i < maxi; i++ )
	{
		double ExpM1;
		/* want the number exp(hnu/kT) - 1, two expansions */
		if( arg[i] < 1.e-5 )
		{
			/* for small arg expand exp(hnu/kT) - 1 to second order */
			ExpM1 = arg[i]*(1. + arg[i]/2.);
		}
		else
		{
			/* for large arg, evaluate the full expansion */
			ExpM1 = expval[i] - 1.;
		}

		double Planck1 = PI4*2.*HPLANCK/POW2(SPEEDLIGHT)*POW2(FR1RYD)*POW2(FR1RYD)*
			rfield.anu3(i)/ExpM1*rfield.widflx(i);
		double Planck2 = Planck1*gv.bin[nd].dstab1[i];

		/* add integral over RJ tail, maybe useful for extreme low temps */
		if( i == 0 ) 
		{
			integral1 = Planck1/rfield.widflx(0)*rfield.anu(0)/3.;
			integral2 = Planck2/rfield.widflx(0)*rfield.anu(0)/5.;
		}
		/* if we are in the Wien tail - exit */
		if( Planck1/integral1 < DBL_EPSILON && Planck2/integral2 < DBL_EPSILON )
			break;

		integral1 += Planck1;
		integral2 += Planck2;
	}

	/* this is an option to print out every few steps, when 'trace grains' is set */
	if( trace.lgDustBug && trace.lgTrace && ip%10 == 0 )
	{
		fprintf( ioQQQ, "  %4ld %11.4e %11.4e %11.4e %11.4e\n", (unsigned long)nd, tdust, 
		  integral2, integral1/4./5.67051e-5/powi(tdust,4), integral2*4./integral1 );
	}

	ASSERT( integral2 > 0. );
	return integral2;
}


/* invalidate charge dependent data from previous iteration */
STATIC void NewChargeData(long nd)
{
	long nz;

	DEBUG_ENTRY( "NewChargeData()" );

	for( nz=0; nz < NCHS; nz++ )
	{
		gv.bin[nd].chrg(nz).RSum1 = -DBL_MAX;
		gv.bin[nd].chrg(nz).RSum2 = -DBL_MAX;
		gv.bin[nd].chrg(nz).ESum1a = -DBL_MAX;
		gv.bin[nd].chrg(nz).ESum1b = -DBL_MAX;
		gv.bin[nd].chrg(nz).ESum2 = -DBL_MAX;

		gv.bin[nd].chrg(nz).hots1 = -DBL_MAX;
		gv.bin[nd].chrg(nz).bolflux1 = -DBL_MAX;
		gv.bin[nd].chrg(nz).pe1 = -DBL_MAX;

		/** \todo	2	should any of the following 3 statements be removed? */
		gv.bin[nd].chrg(nz).ThermRate = -DBL_MAX;
		gv.bin[nd].chrg(nz).GrainHeat = -DBL_MAX;
		gv.bin[nd].chrg(nz).HeatingRate2 = -DBL_MAX;
	}

	if( !fp_equal(phycon.te,gv.GrnRecomTe) )
	{
		for( nz=0; nz < NCHS; nz++ )
		{
			memset( gv.bin[nd].chrg(nz).eta, 0, (LIMELM+2)*sizeof(double) );
			memset( gv.bin[nd].chrg(nz).xi, 0, (LIMELM+2)*sizeof(double) );
		}
	}

	if( nzone != gv.nzone )
	{
		for( nz=0; nz < NCHS; nz++ )
		{
			gv.bin[nd].chrg(nz).hcon1 = -DBL_MAX;
		}
	}
	return;
}


/* GrnStdDpth sets the standard behavior of the grain abundance as a function 
 * of depth into cloud - user-define code should go in GrnVryDpth */
STATIC double GrnStdDpth(long int nd)
{
	double GrnStdDpth_v;

	DEBUG_ENTRY( "GrnStdDpth()" );

	/* NB NB - this routine defines the standard depth dependence of the grain abundance,
	 * to implement user-defined behavior of the abundance (invoked with the "function"
	 * keyword on the command line), modify the routine GrnVryDpth at the end of this file,
	 * DO NOT MODIFY THIS ROUTINE */

	if( gv.bin[nd].nDustFunc == DF_STANDARD )
	{
		if( gv.bin[nd].matType == MAT_PAH || gv.bin[nd].matType == MAT_PAH2 )
		{
			/* default behavior for PAH's */
			if( gv.chPAH_abundance == "H" )
			{
				/* the scale factor is the hydrogen atomic fraction, small when gas is ionized
				 * or molecular and unity when atomic. This function is observed for PAHs
				 * across the Orion Bar, the PAHs are strong near the ionization front and
				 * weak in the ionized and molecular gas */
				/* >>chng 04 sep 28, propto atomic fraction */
				GrnStdDpth_v = dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN];
			}
			else if( gv.chPAH_abundance == "H,H2" )
			{
				/* the scale factor is the total of the hydrogen atomic and molecular fractions,
				 * small when gas is ionized and unity when atomic or molecular. This function is
				 * observed for PAHs across the Orion Bar, the PAHs are strong near the ionization
				 * front and weak in the ionized and molecular gas */
				/* >>chng 04 sep 28, propto atomic fraction */
				GrnStdDpth_v = (dense.xIonDense[ipHYDROGEN][0]+2*hmi.H2_total)/dense.gas_phase[ipHYDROGEN];
			}
			else if( gv.chPAH_abundance == "CON" )
			{
				/* constant abundance - unphysical, used only for testing */
				GrnStdDpth_v = 1.;
			}
			else
			{
				fprintf(ioQQQ,"Invalid argument to SET PAH: %s\n",gv.chPAH_abundance.c_str());
				TotalInsanity();
			}
		}
		else
		{
			/* default behavior for all other types of grains */
			GrnStdDpth_v = 1.;
		}
	}
	else if( gv.bin[nd].nDustFunc == DF_USER_FUNCTION )
	{
		GrnStdDpth_v = GrnVryDpth(nd);
	}
	else if( gv.bin[nd].nDustFunc == DF_SUBLIMATION )
	{
		// abundance depends on temperature relative to sublimation
		// "grain function sublimation" command
		GrnStdDpth_v = sexp( pow3( gv.bin[nd].tedust / gv.bin[nd].Tsublimat ) );
	}
	else
	{
		TotalInsanity();
	}

	GrnStdDpth_v = max(1.e-10,GrnStdDpth_v);

	return GrnStdDpth_v;
}


/* this is the main routine that drives the grain physics */
void GrainDrive()
{
	DEBUG_ENTRY( "GrainDrive()" );

	/* gv.lgGrainPhysicsOn set false with no grain physics command */
	if( gv.lgDustOn() && gv.lgGrainPhysicsOn )
	{
		static double tesave = -1.;
		static long int nzonesave = -1;

		/* option to only reevaluate grain physics if something has changed.  
		 * gv.lgReevaluate is set false with keyword no reevaluate on grains command 
		 * option to force constant reevaluation of grain physics - 
		 * by default is true 
		 * usually reevaluate grains at all times, but NO REEVALUATE will
		 * save some time but may affect stability */
		if( gv.lgReevaluate || conv.lgSearch || nzonesave != nzone || 
			/* need to reevaluate the grains when temp changes since */
			! fp_equal( phycon.te, tesave ) || 
			/* >>chng 03 dec 30, check that electrons locked in grains are not important,
			 * if they are, then reevaluate */
			 fabs(gv.TotalEden)/dense.eden > conv.EdenErrorAllowed/5. ||
			 /* >>chng 04 aug 06, always reevaluate when thermal effects of grains are important,
			  * first is collisional energy exchange with gas, second is grain photoionization */
			 (fabs( gv.GasCoolColl ) + fabs( thermal.heating(0,13) ))/SDIV(thermal.ctot)>0.1 )
		{
			nzonesave = nzone;
			tesave = phycon.te;

			if( trace.nTrConvg >= 5 )
			{
				fprintf( ioQQQ, "        grain5 calling GrainChargeTemp\n");
			}
			/* find dust charge and temperature - this must be called at least once per zone
			 * since grain abundances, set here, may change with depth */
			GrainChargeTemp();

			/* >>chng 04 jan 31, moved call to GrainDrift to ConvPresTempEdenIoniz(), PvH */
		}
	}
	else if( gv.lgDustOn() && !gv.lgGrainPhysicsOn )
	{
		/* very minimalistic treatment of grains; only extinction of continuum is considered
		 * however, the absorbed energy is not reradiated, so this creates thermal imbalance! */
		if( gv.nCalledGrainDrive == 0 )
		{
			long nelem, ion, ion_to;

			/* when not doing grain physics still want some exported quantities
			 * to be reasonable, grain temperature used for H2 formation */
			gv.GasHeatPhotoEl = 0.;
			for( size_t nd=0; nd < gv.bin.size(); nd++ )
			{
				long nz;

				/* this disables warnings about PAHs in the ionized region */
				gv.bin[nd].lgPAHsInIonizedRegion = false;

				for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
				{
					gv.bin[nd].chrg(nz).DustZ = nz;
					gv.bin[nd].chrg(nz).FracPop = ( nz == 0 ) ? 1. : 0.;
					gv.bin[nd].chrg(nz).nfill = 0;
					gv.bin[nd].chrg(nz).tedust = 100.f;
				}

				gv.bin[nd].AveDustZ = 0.;
				gv.bin[nd].dstpot = chrg2pot(0.,nd);

				gv.bin[nd].tedust = 100.f;
				gv.bin[nd].TeGrainMax = 100.;

				/* set all heating/cooling agents to zero */
				gv.bin[nd].BolFlux = 0.;
				gv.bin[nd].GrainCoolTherm = 0.;
				gv.bin[nd].GasHeatPhotoEl = 0.;
				gv.bin[nd].GrainHeat = 0.;
				gv.bin[nd].GrainHeatColl = 0.;
				gv.bin[nd].ChemEn = 0.;
				gv.bin[nd].ChemEnH2 = 0.;
				gv.bin[nd].thermionic = 0.;

				gv.bin[nd].lgUseQHeat = false;
				gv.bin[nd].lgEverQHeat = false;
				gv.bin[nd].QHeatFailures = 0;

				gv.bin[nd].DustDftVel = 0.;

				gv.bin[nd].avdust = gv.bin[nd].tedust;
				gv.bin[nd].avdft = 0.f;
				gv.bin[nd].avdpot = (realnum)(gv.bin[nd].dstpot*EVRYD);
				gv.bin[nd].avDGRatio = -1.f;

				/* >>chng 06 jul 21, add this here as well as in GrainTemperature so that can
				 * get fake heating when grain physics is turned off */
				if( 0 && gv.lgBakesPAH_heat )
				{
					/* this is a dirty hack to get BT94 PE heating rate
					 * for PAH's included, for Lorentz Center 2004 PDR meeting, PvH */
					/*>>>refer	PAH	heating	Bakes, E.L.O., & Tielens, A.G.G.M. 1994, ApJ, 427, 822 */
					/* >>chng 05 aug 12, change from +=, which added additional heating to what exists already,
					 * to simply = to set the heat, this equation gives total heating */
					gv.bin[nd].GasHeatPhotoEl = 1.e-24*hmi.UV_Cont_rel2_Habing_TH85_depth*
						dense.gas_phase[ipHYDROGEN]*(4.87e-2/(1.0+4e-3*pow((hmi.UV_Cont_rel2_Habing_TH85_depth*
						sqrt(phycon.te)/dense.eden),0.73)) + 3.65e-2*pow(phycon.te/1.e4,0.7)/
						(1.+2.e-4*(hmi.UV_Cont_rel2_Habing_TH85_depth*sqrt(phycon.te)/dense.eden)))/gv.bin.size() *
						gv.GrainHeatScaleFactor;
					gv.GasHeatPhotoEl += gv.bin[nd].GasHeatPhotoEl;
				}
			}

			gv.TotalEden = 0.;
			gv.GrnElecDonateMax = 0.f;
			gv.GrnElecHoldMax = 0.f;

			for( nelem=0; nelem < LIMELM; nelem++ )
			{
				for( ion=0; ion <= nelem+1; ion++ )
				{
					for( ion_to=0; ion_to <= nelem+1; ion_to++ )
					{
						gv.GrainChTrRate[nelem][ion][ion_to] = 0.f;
					}
				}
			}

			/* set all heating/cooling agents to zero */
			gv.GrainHeatInc = 0.;
			gv.GrainHeatDif = 0.;
			gv.GrainHeatLya = 0.;
			gv.GrainHeatCollSum = 0.;
			gv.GrainHeatSum = 0.;
			gv.GrainHeatChem = 0.;
			gv.GasCoolColl = 0.;
			gv.TotalDustHeat = 0.f;
			gv.dphmax = 0.f;
			gv.dclmax = 0.f;

			thermal.setHeating(0,13,0.);
			thermal.setHeating(0,14,0.);
			thermal.setHeating(0,25,0.);
		}

		if( gv.nCalledGrainDrive == 0 || gv.lgAnyDustVary )
		{
			GrainUpdateRadius1();
			GrainUpdateRadius2();
		}
	}

	++gv.nCalledGrainDrive;
	return;
}

/* iterate grain charge and temperature */
STATIC void GrainChargeTemp()
{
	long int i,
	  ion,
	  ion_to,
	  nelem,
	  nz;
	realnum dccool = FLT_MAX;
	double delta,
	  GasHeatNet,
	  hcon = DBL_MAX,
	  hla = DBL_MAX,
	  hots = DBL_MAX,
	  oldtemp,
	  oldTotalEden,
	  ratio,
	  ThermRatio;

	static long int oldZone = -1;
	static double oldTe = -DBL_MAX,
	  oldHeat = -DBL_MAX;

	DEBUG_ENTRY( "GrainChargeTemp()" );

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "\n GrainChargeTemp called lgSearch%2c\n\n", TorF(conv.lgSearch) );
	}

	oldTotalEden = gv.TotalEden;

	/* these will sum heating agents over grain populations */
	gv.GrainHeatInc = 0.;
	gv.GrainHeatDif = 0.;
	gv.GrainHeatLya = 0.;
	gv.GrainHeatCollSum = 0.;
	gv.GrainHeatSum = 0.;
	gv.GrainHeatChem = 0.;

	gv.GasCoolColl = 0.;
	gv.GasHeatPhotoEl = 0.;
	gv.GasHeatTherm = 0.;

	gv.TotalEden = 0.;

	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		for( ion=0; ion <= nelem+1; ion++ )
		{
			for( ion_to=0; ion_to <= nelem+1; ion_to++ )
			{
				gv.GrainChTrRate[nelem][ion][ion_to] = 0.f;
			}
		}
	}

	/* this sets dstAbund and conversion factors, but not gv.dstab and gv.dstsc! */
	GrainUpdateRadius1();

	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		double one;
		double ChTdBracketLo = 0., ChTdBracketHi = -DBL_MAX;
		long relax = ( conv.lgSearch ) ? 3 : 1;

		/* >>chng 02 nov 11, added test for the presence of PAHs in the ionized region, PvH */
		if( gv.bin[nd].matType == MAT_PAH || gv.bin[nd].matType == MAT_PAH2 )
		{
			if( dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN] > 0.50 )
			{
				gv.bin[nd].lgPAHsInIonizedRegion = true;
			}
		}

		/* >>chng 01 sep 13, dynamically allocate backup store, remove ncell dependence, PvH */
		/* allocate data inside loop to avoid accidental spillover to next iteration */
		/* >>chng 04 jan 18, no longer delete and reallocate data inside loop to speed up the code, PvH */
		NewChargeData(nd);

		if( trace.lgTrace && trace.lgDustBug )
		{
			fprintf( ioQQQ, " >>GrainChargeTemp starting grain %s\n",
				 gv.bin[nd].chDstLab );
		}

		delta = 2.*TOLER;
		/* >>chng 01 nov 29, relax max no. of iterations during initial search */
		for( i=0; i < relax*CT_LOOP_MAX && delta > TOLER; ++i )
		{
			string which;
			long j;
			double TdBracketLo = 0., TdBracketHi = -DBL_MAX;
			double ThresEst = 0.;
			oldtemp = gv.bin[nd].tedust;

			/* solve for charge using previous estimate for grain temp
			 * grain temp only influences thermionic emissions
			 * Thermratio is fraction thermionic emissions contribute
			 * to the total electron loss rate of the grain */
			GrainCharge(nd,&ThermRatio);

			ASSERT( gv.bin[nd].GrainHeat > 0. );
			ASSERT( gv.bin[nd].tedust >= GRAIN_TMIN && gv.bin[nd].tedust <= GRAIN_TMAX );

			/* >>chng 04 may 31, in conditions where collisions become an important
			 * heating/cooling source (e.g. gas that is predominantly heated by cosmic
			 * rays), the heating rate depends strongly on the assumed dust temperature.
			 * hence it is necessary to iterate for the dust temperature. PvH */
			gv.bin[nd].lgTdustConverged = false;
			for( j=0; j < relax*T_LOOP_MAX; ++j )
			{
				double oldTemp2 = gv.bin[nd].tedust;
				double oldHeat2 = gv.bin[nd].GrainHeat;
				double oldCool = gv.bin[nd].GrainGasCool;

				/* now solve grain temp using new value for grain potential */
				GrainTemperature(nd,&dccool,&hcon,&hots,&hla);

				gv.bin[nd].GrainGasCool = dccool;

				if( trace.lgTrace && trace.lgDustBug )
				{
					fprintf( ioQQQ, "  >>loop %ld BracketLo %.6e BracketHi %.6e",
						 j, TdBracketLo, TdBracketHi );
				}

				/* this test assures that convergence can only happen if GrainHeat > 0
				 * and therefore the value of tedust is guaranteed to be valid as well */
				/* >>chng 04 aug 05, test that gas cooling is converged as well,
				 * in deep PDRs gas cooling depends critically on grain temperature, PvH */
				if( fabs(gv.bin[nd].GrainHeat-oldHeat2) < HEAT_TOLER*gv.bin[nd].GrainHeat &&
				    fabs(gv.bin[nd].GrainGasCool-oldCool) < HEAT_TOLER_BIN*thermal.ctot )
				{
					gv.bin[nd].lgTdustConverged = true;
					if( trace.lgTrace && trace.lgDustBug )
						fprintf( ioQQQ, " converged\n" );
					break;
				}

				/* update the bracket for the solution */
				if( gv.bin[nd].tedust < oldTemp2 )
					TdBracketHi = oldTemp2;
				else
					TdBracketLo = oldTemp2;

				/* GrainTemperature yields a new estimate for tedust, and initially
				 * that estimate will be used. In most zones this will converge quickly.
				 * However, sometimes the solution will oscillate and converge very
				 * slowly. So, as soon as j >= 2 and the bracket is set up, we will
				 * force convergence by using a bisection search within the bracket */
				/** \todo	2	this algorithm might be more efficient with Brent */

				/* this test assures that TdBracketHi is initialized */
				if( TdBracketHi > TdBracketLo )
				{
					/* if j >= 2, the solution is converging too slowly
					 * so force convergence by doing a bisection search */
					if( ( j >= 2 && TdBracketLo > 0. ) ||
					    gv.bin[nd].tedust <= TdBracketLo ||
					    gv.bin[nd].tedust >= TdBracketHi )
					{
						gv.bin[nd].tedust = (realnum)(0.5*(TdBracketLo + TdBracketHi));
						if( trace.lgTrace && trace.lgDustBug )
							fprintf( ioQQQ, " bisection\n" );
					}
					else
					{
						if( trace.lgTrace && trace.lgDustBug )
							fprintf( ioQQQ, " iteration\n" );
					}
				}
				else
				{
					if( trace.lgTrace && trace.lgDustBug )
						fprintf( ioQQQ, " iteration\n" );
				}

				ASSERT( gv.bin[nd].tedust >= GRAIN_TMIN && gv.bin[nd].tedust <= GRAIN_TMAX );
			}

			if( gv.bin[nd].lgTdustConverged )
			{
				/* update the bracket for the solution */
				if( gv.bin[nd].tedust < oldtemp )
					ChTdBracketHi = oldtemp;
				else
					ChTdBracketLo = oldtemp;
			}
			else
			{
				bool lgBoundErr;
				double y, x = log(gv.bin[nd].tedust);
				/* make sure GrainHeat is consistent with value of tedust */
				splint_safe(gv.dsttmp,gv.bin[nd].dstems,gv.bin[nd].dstslp2,NDEMS,x,&y,&lgBoundErr);
				gv.bin[nd].GrainHeat = exp(y)*gv.bin[nd].cnv_H_pCM3;

				fprintf( ioQQQ," PROBLEM  temperature of grain species %s (Tg=%.3eK) not converged\n",
					 gv.bin[nd].chDstLab , gv.bin[nd].tedust );
				ConvFail("grai","");
			}

			ASSERT( gv.bin[nd].GrainHeat > 0. );
			ASSERT( gv.bin[nd].tedust >= GRAIN_TMIN && gv.bin[nd].tedust <= GRAIN_TMAX );

			/* delta estimates relative change in electron emission rate
			 * due to the update in the grain temperature, if it is small
			 * we won't bother to iterate (which is usually the case)
			 * the formula assumes that thermionic emission is the only
			 * process that depends on grain temperature */
			/** \todo	2	should collisional heating/cooling be included here? */
			ratio = gv.bin[nd].tedust/oldtemp;
			for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
			{
				ThresEst += gv.bin[nd].chrg(nz).FracPop*gv.bin[nd].chrg(nz).ThresInf;
			}
			delta = ThresEst*TE1RYD/gv.bin[nd].tedust*(ratio - 1.);
			/** \todo	2	use something like log(ThermRatio) + log(delta) ???? */
			delta = ( delta < 0.9*log(DBL_MAX) ) ?
				ThermRatio*fabs(POW2(ratio)*exp(delta)-1.) : DBL_MAX;

			/* >>chng 06 feb 07, bracket grain temperature to force convergence when oscillating, PvH */
			if( delta > TOLER )
			{
				if( trace.lgTrace && trace.lgDustBug )
					which = "iteration";

				/* The loop above yields a new estimate for tedust, and initially that
				 * estimate will be used. In most zones this will converge very quickly.
				 * However, sometimes the solution will oscillate and converge very
				 * slowly. So, as soon as i >= 2 and the bracket is set up, we will
				 * force convergence by using a bisection search within the bracket */
				/** \todo	2	this algorithm might be more efficient with Brent */

				/* this test assures that ChTdBracketHi is initialized */
				if( ChTdBracketHi > ChTdBracketLo )
				{
					/* if i >= 2, the solution is converging too slowly
					 * so force convergence by doing a bisection search */
					if( ( i >= 2 && ChTdBracketLo > 0. ) ||
					    gv.bin[nd].tedust <= ChTdBracketLo ||
					    gv.bin[nd].tedust >= ChTdBracketHi )
					{
						gv.bin[nd].tedust = (realnum)(0.5*(ChTdBracketLo + ChTdBracketHi));
						if( trace.lgTrace && trace.lgDustBug )
							which = "bisection";
					}
				}
			}

			if( trace.lgTrace && trace.lgDustBug )
			{
				fprintf( ioQQQ, " >>GrainChargeTemp finds delta=%.4e, ", delta );
				fprintf( ioQQQ, " old/new temp=%.5e %.5e, ", oldtemp, gv.bin[nd].tedust );
				if( delta > TOLER ) 
					fprintf( ioQQQ, "doing another %s\n", which.c_str() );
				else 
					fprintf( ioQQQ, "converged\n" );
			}
		}
		if( delta > TOLER )
		{
			fprintf( ioQQQ, " PROBLEM  charge/temperature not converged for %s zone %.2f\n",
				 gv.bin[nd].chDstLab , fnzone );
			ConvFail("grai","");
		}

		/* add in ion recombination rates on this grain species */
		/* ionbal.lgGrainIonRecom is 1 by default, set to 0 with
		 * no grain neutralization command */
		if( ionbal.lgGrainIonRecom )
			GrainChrgTransferRates(nd);

		/* >>chng 04 jan 31, moved call to UpdateRadius2 outside loop, PvH */

		/* following used to keep track of heating agents in printout
		 * no physics done with GrainHeatInc
		 * dust heating by incident continuum, and elec friction before ejection */
		gv.GrainHeatInc += hcon;
		/* remember total heating by diffuse fields, for printout (includes Lya) */
		gv.GrainHeatDif += hots;
		/* GrainHeatLya - total heating by LA in this zone, erg cm-3 s-1, only here
		 * for eventual printout, hots is total ots line heating */
		gv.GrainHeatLya += hla;

		/* this will be total collisional heating, for printing in lines */
		gv.GrainHeatCollSum += gv.bin[nd].GrainHeatColl;

		/* GrainHeatSum is total heating of all grain types in this zone,
		 * will be carried by total cooling, only used in lines to print tot heat
		 * printed as entry "GraT    0 " */
		gv.GrainHeatSum += gv.bin[nd].GrainHeat;

		/* net amount of chemical energy donated by recombining ions and molecule formation */
		gv.GrainHeatChem += gv.bin[nd].ChemEn + gv.bin[nd].ChemEnH2;

		/* dccool is gas cooling due to collisions with grains - negative if net heating 
		 * zero if NO GRAIN GAS COLLISIONAL EXCHANGE command included */
		gv.GasCoolColl += dccool;
		gv.GasHeatPhotoEl += gv.bin[nd].GasHeatPhotoEl;
		gv.GasHeatTherm += gv.bin[nd].thermionic;

		/* this is grain charge in e/cm^3, positive number means grain supplied free electrons */
		/* >>chng 01 mar 24, changed DustZ+1 to DustZ, PvH */
		one = 0.;
		for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
		{
			one += gv.bin[nd].chrg(nz).FracPop*(double)gv.bin[nd].chrg(nz).DustZ*
				gv.bin[nd].cnv_GR_pCM3;
		}
		/* electron density contributed by grains, cm-3 */
		gv.TotalEden += one;
		{
			/*@-redef@*/
			enum {DEBUG_LOC=false};
			/*@+redef@*/
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ," DEBUG grn chr nz\t%.2f\teden\t%.3e\tnd\t%li",
					fnzone,
					dense.eden,
					(unsigned long)nd);
				fprintf(ioQQQ,"\tne\t%.2e\tAveDustZ\t%.2e\t%.2e\t%.2e\t%.2e",
					one,
					gv.bin[nd].AveDustZ,
					gv.bin[nd].chrg(0).FracPop,(double)gv.bin[nd].chrg(0).DustZ,
					gv.bin[nd].cnv_GR_pCM3);
				fprintf(ioQQQ,"\n");
			}
		}

		if( trace.lgTrace && trace.lgDustBug )
		{
			fprintf(ioQQQ,"     %s Pot %.5e Thermal %.5e GasCoolColl %.5e" , 
				gv.bin[nd].chDstLab, gv.bin[nd].dstpot, gv.bin[nd].GrainHeat, dccool );
			fprintf(ioQQQ," GasPEHeat %.5e GasThermHeat %.5e ChemHeat %.5e\n\n" , 
				gv.bin[nd].GasHeatPhotoEl, gv.bin[nd].thermionic, gv.bin[nd].ChemEn );
		}
	}

	/* >>chng 04 aug 06, added test of convergence of the net gas heating/cooling, PvH */
	GasHeatNet = gv.GasHeatPhotoEl + gv.GasHeatTherm - gv.GasCoolColl;

	if( !fp_equal(phycon.te,gv.GrnRecomTe) )
	{
		oldZone = gv.nzone;
		oldTe = gv.GrnRecomTe;
		oldHeat = gv.GasHeatNet;
	}

	/* >>chng 04 aug 07, added estimate for heating derivative, PvH */
	if( nzone == oldZone && !fp_equal(phycon.te,oldTe) )
	{
		gv.dHeatdT = (GasHeatNet-oldHeat)/(phycon.te-oldTe);
	}

	/* >>chng 04 sep 15, add test for convergence of gv.TotalEden, PvH */
	if( nzone != gv.nzone || !fp_equal(phycon.te,gv.GrnRecomTe) ||
	    fabs(gv.GasHeatNet-GasHeatNet) > HEAT_TOLER*thermal.ctot ||
	    fabs(gv.TotalEden-oldTotalEden) > CHRG_TOLER*dense.eden )
	{
		/* >>chng 04 aug 07, add test whether eden on grain converged */
		/* flag that change in eden was too large */
		/*conv.lgConvEden = false;*/
		if( fabs(gv.TotalEden-oldTotalEden) > CHRG_TOLER*dense.eden )
		{
			conv.setConvIonizFail( "grn eden chg" , oldTotalEden, gv.TotalEden);
		}
		else if( fabs(gv.GasHeatNet-GasHeatNet) > HEAT_TOLER*thermal.ctot )
		{
			conv.setConvIonizFail( "grn het chg" , gv.GasHeatNet, GasHeatNet);
		}
		else if( !fp_equal(phycon.te,gv.GrnRecomTe) )
		{
			conv.setConvIonizFail( "grn ter chg" , gv.GrnRecomTe, phycon.te);
		}
		else if( nzone != gv.nzone )
		{
			conv.setConvIonizFail( "grn zon chg" , gv.nzone, nzone );
		}
		else
			TotalInsanity();
	}

	/* printf( "DEBUG GasHeatNet %.6e -> %.6e TotalEden %e -> %e conv.lgConvIoniz %c\n",
	   gv.GasHeatNet, GasHeatNet, gv.TotalEden, oldTotalEden, TorF(conv.lgConvIoniz) ); */
	/* printf( "DEBUG %.2f %e %e\n", fnzone, phycon.te, dense.eden ); */

	/* update total grain opacities in gv.dstab and gv.dstsc,
	 * they depend on grain charge and may depend on depth
	 * also add in the photo-dissociation cs in gv.dstab */
	GrainUpdateRadius2();

	gv.nzone = nzone;
	gv.GrnRecomTe = phycon.te;
	gv.GasHeatNet = GasHeatNet;

#ifdef WD_TEST2
	printf("wd test: proton fraction %.5e Total DustZ %.6f heating/cooling rate %.5e %.5e\n",
	       dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN],
	       gv.bin[0]->AveDustZ,gv.GasHeatPhotoEl/dense.gas_phase[ipHYDROGEN]/fudge(0),
	       gv.GasCoolColl/dense.gas_phase[ipHYDROGEN]/fudge(0));
#endif

	if( trace.lgTrace )
	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			fprintf( ioQQQ, "     %.2f Grain surface charge transfer rates\n", fnzone );

			for( nelem=0; nelem < LIMELM; nelem++ )
			{
				if( dense.lgElmtOn[nelem] )
				{
					printf( "      %s:", elementnames.chElementSym[nelem] );
					for( ion=dense.IonLow[nelem]; ion <= dense.IonHigh[nelem]; ++ion )
					{
						for( ion_to=0; ion_to <= nelem+1; ion_to++ )
						{
							if( gv.GrainChTrRate[nelem][ion][ion_to] > 0.f )
							{
								printf( "  %ld->%ld %.2e", ion, ion_to,
									gv.GrainChTrRate[nelem][ion][ion_to] );
							}
						}
					}
					printf( "\n" );
				}
			}
		}

		fprintf( ioQQQ, "     %.2f Grain contribution to electron density %.2e\n", 
			fnzone , gv.TotalEden );

		fprintf( ioQQQ, "     Grain electons: " );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			double sum = 0.;
			for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
			{
				sum += gv.bin[nd].chrg(nz).FracPop*(double)gv.bin[nd].chrg(nz).DustZ*
					gv.bin[nd].cnv_GR_pCM3;
			}
			fprintf( ioQQQ, " %.2e", sum );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "     Grain potentials:" );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, " %.2e", gv.bin[nd].dstpot );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "     Grain temperatures:" );
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			fprintf( ioQQQ, " %.2e", gv.bin[nd].tedust );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "     GrainCollCool: %.6e\n", gv.GasCoolColl );
	}

	/*if( nzone > 900) 
		fprintf(ioQQQ,"DEBUG cool\t%.2f\t%e\t%e\t%e\n",
		fnzone,
		phycon.te ,
		dense.eden,
		gv.GasCoolColl );*/
	return;
}


STATIC void GrainCharge(size_t nd,
			/*@out@*/double *ThermRatio) /* ratio of thermionic to total rate */
{
	bool lgBigError;
	long backup,
	  i,
	  loopMax,
	  newZlo,
	  nz,
	  power,
	  stride,
	  stride0,
	  Zlo;
	double crate,
	  csum1,
	  csum1a,
	  csum1b,
	  csum2,
	  csum3,
	  netloss0 = -DBL_MAX,
	  netloss1 = -DBL_MAX,
	  rate_up[NCHU],
	  rate_dn[NCHU],
	  step;

	DEBUG_ENTRY( "GrainCharge()" );

	/* find dust charge */
	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "    Charge loop, search %c,", TorF(conv.lgSearch) );
	}

	ASSERT( nd < gv.bin.size() );

	for( nz=0; nz < NCHS; nz++ )
	{
		gv.bin[nd].chrg(nz).FracPop = -DBL_MAX;
	}

	/* this algorithm determines the value of Zlo and the charge state populations
	 * in the n-charge state model as described in:
	 *
	 * >>refer	grain	physics	van Hoof P.A.M., Weingartner J.C., et al., 2004, MNRAS, 350, 1330
	 *
	 * the algorithm first uses the n charge states to bracket the solution by
	 * separating the charge states with a stride that is an integral power of
	 * n-1, e.g. like this for an n == 4 calculation:
	 *
	 *  +gain    +gain    -gain    -gain
	 *    |        |        |        |
	 *   -8        1        10       19
	 *
	 * for each of the charge states the total electron emission and recombination
	 * rates are calculated. the solution is considered properly bracketed if the
	 * sign of the net gain rate (emission-recombination) is different for the first
	 * and the last of the n charge states. then the algorithm searches for two
	 * adjacent charge states where the net gain rate changes sign and divides
	 * that space into n-1 equal parts, like here
	 *
	 *   net gain  +  +  +  -
	 *             |  |  |  |
	 *        Zg   1  4  7  10
	 *
	 * note that the first and the last charge state can be retained from the
	 * previous iteration and do not need to be recalculated (UpdatePot as well
	 * as GrainElecEmis1 and GrainElecRecomb1 have mechanisms for re-using
	 * previously calculated data, so GrainCharge doesn't need to worry about
	 * this). The dividing algorithm is repeated until the stride equals 1, like
	 * here
	 *
	 *        net gain   +---
	 *                   ||||
	 *             Zg    7  10
	 *
	 * finally, the bracket may need to be shifted in order for the criterion
	 * J1 x J2 <= 0 to be fulfilled (see the paper quoted above for a detailed
	 * explanation). in the example here one shift might be necessary:
	 *
	 *        net gain  ++--
	 *                  ||||
	 *             Zg   6  9
	 *
	 * for n == 2, the algorithm would have to be slightly different. In order
	 * to avoid complicating the code, we force the code to use n == 3 in the
	 * first two steps of the process, reverting back to n == 2 upon the last
	 * step. this should not produce any noticeable additional CPU overhead */

	backup = gv.bin[nd].nChrg;
	gv.bin[nd].nChrg = MAX2(gv.bin[nd].nChrg,3);

	stride0 = gv.bin[nd].nChrg-1;

	/* set up initial bracket for grain charge, will be checked below */
	if( gv.bin[nd].lgIterStart )
	{
		if( iteration == 1 )
		{
			// this is the very first time the grain charge is determined
			// so here we try to get a rough first estimate of the charge
			long ilo = rfield.ipointC(0.9);
			long ihi = rfield.nflux;
			// determine average photon energy above 0.9 Ryd to ionize grain
			double sum1, sum2 = reduce_ab_a( get_ptr(rfield.flux[0]), rfield.anuptr(), ilo, ihi, &sum1 );
			double anuav = safe_div( sum2, sum1, 0. );
			if( anuav > 1. )
			{
				// this is a hard radiation field -> assume grains are positively charged
				Zlo = 0;
				step = max(pot2chrg(anuav/2.,nd),1.);
			}
			else
			{
				// this is likely a PDR model -> assume grains are somewhat negative
				step = max(pot2chrg(0.2,nd),1.);
				Zlo = -long(nint(step/4.));
			}
			power = (int)(log(step)/log((double)stride0));
			power = MAX2(power,0);
			double xxx = powi((double)stride0,power);
			stride = nint(xxx);
		}
		else
		{
			// use the initial solution from the previous iteration
			stride = 1;
			Zlo = gv.bin[nd].ZloSave;
		}
	}
	else
	{
		/* the previous solution is the best choice here */
		stride = 1;
		Zlo = gv.bin[nd].chrg(0).DustZ;
	}
	Zlo = max(Zlo,gv.bin[nd].LowestZg);
	UpdatePot( nd, Zlo, stride, rate_up, rate_dn );

	/* check if the grain charge is correctly bracketed */
	for( i=0; i < BRACKET_MAX; i++ )
	{
		netloss0 = rate_up[0] - rate_dn[0];
		netloss1 = rate_up[gv.bin[nd].nChrg-1] - rate_dn[gv.bin[nd].nChrg-1];

		if( netloss0*netloss1 <= 0. )
			break;

		if( netloss1 > 0. )
			Zlo += (gv.bin[nd].nChrg-1)*stride;

		if( i > 0 )
			stride *= stride0;

		if( netloss1 < 0. )
			Zlo -= (gv.bin[nd].nChrg-1)*stride;

		Zlo = MAX2(Zlo,gv.bin[nd].LowestZg);
		UpdatePot( nd, Zlo, stride, rate_up, rate_dn );
	}

	if( netloss0*netloss1 > 0. ) {
		fprintf( ioQQQ, " insanity: could not bracket grain charge for %s\n", gv.bin[nd].chDstLab );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* home in on the charge */
	while( stride > 1 )
	{
		stride /= stride0;

		netloss0 = rate_up[0] - rate_dn[0];
		for( nz=0; nz < gv.bin[nd].nChrg-1; nz++ )
		{
			netloss1 = rate_up[nz+1] - rate_dn[nz+1];

			if( netloss0*netloss1 <= 0. )
			{
				Zlo = gv.bin[nd].chrg(nz).DustZ;
				break;
			}

			netloss0 = netloss1;
		}
		UpdatePot( nd, Zlo, stride, rate_up, rate_dn );
	}

	ASSERT( netloss0*netloss1 <= 0. );

	gv.bin[nd].nChrg = backup;

	/* >>chng 04 feb 15, relax upper limit on initial search when nChrg is much too large, PvH */ 
	loopMax = ( gv.bin[nd].lgIterStart ) ? 4*gv.bin[nd].nChrg : 2*gv.bin[nd].nChrg;

	lgBigError = true;
	for( i=0; i < loopMax; i++ )
	{
		GetFracPop( nd, Zlo, rate_up, rate_dn, &newZlo );

		if( newZlo == Zlo )
		{
			lgBigError = false;
			break;
		}

		Zlo = newZlo;
		UpdatePot( nd, Zlo, 1, rate_up, rate_dn );
	}

	if( lgBigError ) {
		fprintf( ioQQQ, " insanity: could not converge grain charge for %s\n", gv.bin[nd].chDstLab );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	gv.bin[nd].AveDustZ = 0.;
	crate = csum3 = 0.;
	for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
	{
		double d[4];

		gv.bin[nd].AveDustZ += gv.bin[nd].chrg(nz).FracPop*gv.bin[nd].chrg(nz).DustZ;

		crate += gv.bin[nd].chrg(nz).FracPop*GrainElecEmis1(nd,nz,&d[0],&d[1],&d[2],&d[3]);
		csum3 += gv.bin[nd].chrg(nz).FracPop*d[3];
	}
	gv.bin[nd].dstpot = chrg2pot(gv.bin[nd].AveDustZ,nd);
	*ThermRatio = ( crate > 0. ) ? csum3/crate : 0.;

	ASSERT( *ThermRatio >= 0. );

	gv.bin[nd].lgIterStart = false;

	if( trace.lgTrace && trace.lgDustBug )
	{
		double d[4];

		fprintf( ioQQQ, "\n" );

		crate = csum1a = csum1b = csum2 = csum3 = 0.;
		for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
		{
			crate += gv.bin[nd].chrg(nz).FracPop*
				GrainElecEmis1(nd,nz,&d[0],&d[1],&d[2],&d[3]);
			csum1a += gv.bin[nd].chrg(nz).FracPop*d[0];
			csum1b += gv.bin[nd].chrg(nz).FracPop*d[1];
			csum2 += gv.bin[nd].chrg(nz).FracPop*d[2];
			csum3 += gv.bin[nd].chrg(nz).FracPop*d[3];
		}

		fprintf( ioQQQ, "    ElecEm  rate1a=%.4e, rate1b=%.4e, ", csum1a, csum1b );
		fprintf( ioQQQ, "rate2=%.4e, rate3=%.4e, sum=%.4e\n", csum2, csum3, crate );
		if( crate > 0. ) 
		{
			fprintf( ioQQQ, "      rate1a/sum=%.4e, rate1b/sum=%.4e, ", csum1a/crate, csum1b/crate );
			fprintf( ioQQQ, "rate2/sum=%.4e, rate3/sum=%.4e\n", csum2/crate, csum3/crate );
		}

		crate = csum1 = csum2 = 0.;
		for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
		{
			crate += gv.bin[nd].chrg(nz).FracPop*GrainElecRecomb1(nd,nz,&d[0],&d[1]);
			csum1 += gv.bin[nd].chrg(nz).FracPop*d[0];
			csum2 += gv.bin[nd].chrg(nz).FracPop*d[1];
		}

		fprintf( ioQQQ, "    ElecRc  rate1=%.4e, rate2=%.4e, sum=%.4e\n", csum1, csum2, crate );
		if( crate > 0. ) 
		{
			fprintf( ioQQQ, "      rate1/sum=%.4e, rate2/sum=%.4e\n", csum1/crate, csum2/crate );
		}

		fprintf( ioQQQ, "    Charging rates:" );
		for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
		{
			fprintf( ioQQQ, "    Zg %ld up %.4e dn %.4e",
				 gv.bin[nd].chrg(nz).DustZ, rate_up[nz], rate_dn[nz] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "    FracPop: fnzone %.2f nd %ld AveZg %.5e",
			 fnzone, (unsigned long)nd, gv.bin[nd].AveDustZ );
		for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
		{
			fprintf( ioQQQ, "    Zg %ld %.5f", Zlo+nz, gv.bin[nd].chrg(nz).FracPop );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "  >Grain potential:%12.12s %.6fV", 
			 gv.bin[nd].chDstLab, gv.bin[nd].dstpot*EVRYD );
		for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
		{
			fprintf( ioQQQ, "    Thres[%ld]: %.4f V ThresVal[%ld]: %.4f V", 
				 gv.bin[nd].chrg(nz).DustZ, gv.bin[nd].chrg(nz).ThresInf*EVRYD,
				 gv.bin[nd].chrg(nz).DustZ, gv.bin[nd].chrg(nz).ThresInfVal*EVRYD );
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}


/* grain electron recombination rates for single charge state */
STATIC double GrainElecRecomb1(size_t nd,
			       long nz,
			       /*@out@*/ double *sum1,
			       /*@out@*/ double *sum2)
{
	long ion,
	  nelem;
	double eta,
	  rate,
	  Stick,
	  ve,
	  xi;

	DEBUG_ENTRY( "GrainElecRecomb1()" );

	ASSERT( nd < gv.bin.size() );
	ASSERT( nz >= 0 && nz < gv.bin[nd].nChrg );

	/* >>chng 01 may 31, try to find cached results first */
	/* >>chng 04 feb 11, moved cached data to ChargeBin, PvH */
	if( gv.bin[nd].chrg(nz).RSum1 >= 0. )
	{
		*sum1 = gv.bin[nd].chrg(nz).RSum1;
		*sum2 = gv.bin[nd].chrg(nz).RSum2;
		rate = *sum1 + *sum2;
		return rate;
	}

	/* -1 makes psi correct for impact by electrons */
	ion = -1;
	/* VE is mean (not RMS) electron velocity */
	/*ve = TePowers.sqrte*6.2124e5;*/
	ve = sqrt(8.*BOLTZMANN/PI/ELECTRON_MASS*phycon.te);

	Stick = ( gv.bin[nd].chrg(nz).DustZ <= -1 ) ? gv.bin[nd].StickElecNeg : gv.bin[nd].StickElecPos;
	/* >>chng 00 jul 19, replace classical results with results including image potential
	 * to correct for polarization of the grain as charged particle approaches. */
	GrainScreen(ion,nd,nz,&eta,&xi);
	/* this is grain surface recomb rate for electrons */
	*sum1 = ( gv.bin[nd].chrg(nz).DustZ > gv.bin[nd].LowestZg ) ? Stick*dense.eden*ve*eta : 0.;

	/* >>chng 00 jul 13, add in gain rate from atoms and ions, PvH */
	*sum2 = 0.;

#ifndef IGNORE_GRAIN_ION_COLLISIONS
	for( ion=0; ion <= LIMELM; ion++ )
	{
		double CollisionRateAll = 0.;

		for( nelem=MAX2(ion-1,0); nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] && dense.xIonDense[nelem][ion] > 0. &&
			    gv.bin[nd].chrg(nz).RecomZ0[nelem][ion] > ion )
			{
				/* this is rate with which charged ion strikes grain */
				CollisionRateAll += STICK_ION*dense.xIonDense[nelem][ion]*
					GetAveVelocity( dense.AtomicWeight[nelem] )*
					(double)(gv.bin[nd].chrg(nz).RecomZ0[nelem][ion]-ion);
			}
		}

		if( CollisionRateAll > 0. )
		{
			/* >>chng 00 jul 19, replace classical results with results
			 * including image potential to correct for polarization 
			 * of the grain as charged particle approaches. */
			GrainScreen(ion,nd,nz,&eta,&xi);
			*sum2 += CollisionRateAll*eta;
		}
	}
#endif

	rate = *sum1 + *sum2;

	/* >>chng 01 may 31, store results so that they may be used agian */
	gv.bin[nd].chrg(nz).RSum1 = *sum1;
	gv.bin[nd].chrg(nz).RSum2 = *sum2;

	ASSERT( *sum1 >= 0. && *sum2 >= 0. );
	return rate;
}


/* grain electron emission rates for single charge state */
STATIC double GrainElecEmis1(size_t nd,
			     long nz,
			     /*@out@*/ double *sum1a,
			     /*@out@*/ double *sum1b,
			     /*@out@*/ double *sum2,
			     /*@out@*/ double *sum3)
{
	DEBUG_ENTRY( "GrainElecEmis1()" );

	ASSERT( nd < gv.bin.size() );
	ASSERT( nz >= 0 && nz < gv.bin[nd].nChrg );

	/* >>chng 01 may 31, try to find cached results first */
	/* >>chng 04 feb 11, moved cached data to ChargeBin, PvH */
	if( gv.bin[nd].chrg(nz).ESum1a >= 0. )
	{
		*sum1a = gv.bin[nd].chrg(nz).ESum1a;
		*sum1b = gv.bin[nd].chrg(nz).ESum1b;
		*sum2 = gv.bin[nd].chrg(nz).ESum2;
		/* don't cache thermionic rates as they depend on grain temp */
		*sum3 = 4.*gv.bin[nd].chrg(nz).ThermRate;
		return *sum1a + *sum1b + *sum2 + *sum3;
	}

	/* this is the loss rate due to photo-electric effect (includes Auger and secondary electrons) */
	/* >>chng 01 dec 18, added code for modeling secondary electron emissions, PvH */
	/* this code does a crude correction for the Auger effect in grains,
	 * it is roughly valid for neutral and negative grains, but overestimates
	 * the effect for positively charged grains. Many of the Auger electrons have
	 * rather low energies and will not make it out of the potential well for
	 * high grain potentials typical of AGN conditions, see Table 4.1 & 4.2 of
	 * >>refer	grain	physics	Dwek E. & Smith R.K., 1996, ApJ, 459, 686 */
	/* >>chng 06 jan 31, this code has been completely rewritten following
	 * >>refer	grain	physics	Weingartner J.C., Draine B.T, Barr D.K., 2006, ApJ, 645, 1188 */

	ChargeBin& gptr = gv.bin[nd].chrg(nz);

	long ipLo = gptr.ipThresInfVal;
	long ipHi = rfield.nPositive;
#	ifdef WD_TEST2
	*sum1a = reduce_abc( get_ptr(gv.bin[nd].dstab1), get_ptr(rfield.flux[0]), gptr.yhat.ptr0(), ipLo, ipHi );
#	else
	*sum1a = reduce_abc( get_ptr(gv.bin[nd].dstab1), get_ptr(rfield.SummedCon), gptr.yhat.ptr0(), ipLo, ipHi );
#	endif
	/* normalize to rates per cm^2 of projected grain area */
	*sum1a /= gv.bin[nd].IntArea/4.;

	*sum1b = 0.;
	if( gptr.DustZ <= -1 )
	{
		ipLo = gptr.ipThresInf;
		/* >>chng 00 jul 17, use description of Weingartner & Draine, 2001 */
#		ifdef WD_TEST2
		*sum1b = reduce_ab( gptr.cs_pdt.ptr0(), get_ptr(rfield.flux[0]), ipLo, ipHi );
#		else
		*sum1b = reduce_ab( gptr.cs_pdt.ptr0(), get_ptr(rfield.SummedCon), ipLo, ipHi );
#		endif
		*sum1b /= gv.bin[nd].IntArea/4.;
	}

	/* >>chng 00 jun 19, add in loss rate due to recombinations with ions, PvH */
	*sum2 = 0.;
#	ifndef IGNORE_GRAIN_ION_COLLISIONS
	for( long ion=0; ion <= LIMELM; ion++ )
	{
		double CollisionRateAll = 0.;
		for( long nelem=MAX2(ion-1,0); nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] && dense.xIonDense[nelem][ion] > 0. &&
			    ion > gptr.RecomZ0[nelem][ion] )
			{
				/* this is rate with which charged ion strikes grain */
				CollisionRateAll += dense.xIonDense[nelem][ion]*
					GetAveVelocity( dense.AtomicWeight[nelem] )*
					(double)(ion-gptr.RecomZ0[nelem][ion]);
			}
		}
		CollisionRateAll *= STICK_ION;

		if( CollisionRateAll > 0. )
		{
			/* >>chng 00 jul 19, replace classical results with results
			 * including image potential to correct for polarization 
			 * of the grain as charged particle approaches. */
			double eta, xi;
			GrainScreen(ion,nd,nz,&eta,&xi);
			*sum2 += CollisionRateAll*eta;
		}
	}
#	endif

	/* >>chng 01 may 30, moved calculation of ThermRate to UpdatePot */
	/* >>chng 01 jan 19, multiply by 4 since thermionic emissions scale with total grain
	 * surface area while the above two processes scale with projected grain surface area, PvH */
	*sum3 = 4.*gv.bin[nd].chrg(nz).ThermRate;

	/** \todo - add ionizations due to cosmic rays */

	/* >>chng 01 may 31, store results so that they may be used agian */
	gptr.ESum1a = *sum1a;
	gptr.ESum1b = *sum1b;
	gptr.ESum2 = *sum2;

	ASSERT( *sum1a >= 0. && *sum1b >= 0. && *sum2 >= 0. && *sum3 >= 0. );
	return *sum1a + *sum1b + *sum2 + *sum3;
}


/* correction factors for grain charge screening (including image potential
 * to correct for polarization of the grain as charged particle approaches). */
STATIC void GrainScreen(long ion,
			size_t nd,
			long nz,
			/*@out@*/ double *eta,
			/*@out@*/ double *xi)
{
	DEBUG_ENTRY( "GrainScreen()" );

	/* >>chng 04 jan 31, start caching eta, xi results, PvH */
	/* add 1 to allow for electron charge ion = -1 */
	long ind = ion+1;

	ASSERT( ind >= 0 && ind < LIMELM+2 );

	if( gv.bin[nd].chrg(nz).eta[ind] > 0. )
	{
		*eta = gv.bin[nd].chrg(nz).eta[ind];
		*xi = gv.bin[nd].chrg(nz).xi[ind];
		return;
	}

	/* >>refer	grain	physics	Draine & Sutin, 1987, ApJ, 320, 803
	 * eta = J-tilde (eq. 3.3 thru 3.5), xi = Lambda-tilde/2. (eq. 3.8 thru 3.10) */
	if( ion == 0 ) 
	{
		*eta = 1.;
		*xi = 1.;
	}
	else 
	{
		/* >>chng 01 jan 03, assume that grain charge is distributed in two states just below and
		 * above the average charge, instead of the delta function we assume elsewhere. by averaging
		 * over the distribution we smooth out the discontinuities of the the Draine & Sutin expressions
		 * around nu == 0. they were the cause of temperature instabilities in globule.in. PvH */
		/* >>chng 01 may 07, revert back to single charge state since only integral charge states are
		 * fed into this routine now, making the two-charge state approximation obsolete, PvH */
		double nu = (double)gv.bin[nd].chrg(nz).DustZ/(double)ion;
		double tau = gv.bin[nd].Capacity*BOLTZMANN*phycon.te*1.e-7/POW2((double)ion*ELEM_CHARGE);
		if( nu < 0. )
		{
			*eta = (1. - nu/tau)*(1. + sqrt(2./(tau - 2.*nu)));
			*xi = (1. - nu/(2.*tau))*(1. + 1./sqrt(tau - nu));
		}
		else if( nu == 0. ) 
		{
			*eta = 1. + sqrt(PI/(2.*tau));
			*xi = 1. + 0.75*sqrt(PI/(2.*tau));
		}
		else 
		{
			double theta_nu = ThetaNu(nu);
			/* >>>chng 00 jul 27, avoid passing functions to macro so set to local temp var */
			double xxx = 1. + 1./sqrt(4.*tau+3.*nu);
			*eta = POW2(xxx)*exp(-theta_nu/tau);
#			ifdef WD_TEST2
			*xi = (1. + nu/(2.*tau))*(1. + 1./sqrt(3./(2.*tau)+3.*nu))*exp(-theta_nu/tau);
#			else
			/* >>chng 01 jan 24, use new expression for xi which only contains the excess
			 * energy above the potential barrier of the incoming particle (accurate to
			 * 2% or better), and then add in potential barrier separately, PvH */
			xxx = 0.25*powpq(nu/tau,3,4)/(powpq(nu/tau,3,4) + powpq((25.+3.*nu)/5.,3,4)) +
				(1. + 0.75*sqrt(PI/(2.*tau)))/(1. + sqrt(PI/(2.*tau)));
			*xi = (MIN2(xxx,1.) + theta_nu/(2.*tau))*(*eta);
#			endif
		}

		ASSERT( *eta >= 0. && *xi >= 0. );
	}

	gv.bin[nd].chrg(nz).eta[ind] = *eta;
	gv.bin[nd].chrg(nz).xi[ind] = *xi;

	return;
}


STATIC double ThetaNu(double nu)
{
	DEBUG_ENTRY( "ThetaNu()" );

	double theta_nu;
	if( nu > 0. )
	{
		double xxx;
		const double REL_TOLER = 10.*DBL_EPSILON;

		/* >>chng 01 jan 24, get first estimate for xi_nu and iteratively refine, PvH */
		double xi_nu = 1. + 1./sqrt(3.*nu);
		double xi_nu2 = POW2(xi_nu);
		do 
		{
			double old = xi_nu;
			/* >>chng 04 feb 01, use Newton-Raphson to speed up convergence, PvH */
                        /* xi_nu = sqrt(1.+sqrt((2.*POW2(xi_nu)-1.)/(nu*xi_nu))); */
			double fnu = 2.*xi_nu2 - 1. - nu*xi_nu*POW2(xi_nu2 - 1.);
			double dfdxi = 4.*xi_nu - nu*((5.*xi_nu2 - 6.)*xi_nu2 + 1.);
			xi_nu -= fnu/dfdxi;
			xi_nu2 = POW2(xi_nu);
			xxx = fabs(old-xi_nu);
		} while( xxx > REL_TOLER*xi_nu );

		theta_nu = nu/xi_nu - 1./(2.*xi_nu2*(xi_nu2-1.));
	}
	else
	{
		theta_nu = 0.;
	}
	return theta_nu;
}


/* update items that depend on grain potential */
STATIC void UpdatePot(size_t nd,
		      long Zlo,
		      long stride,
		      /*@out@*/ double rate_up[], /* rate_up[NCHU] */
		      /*@out@*/ double rate_dn[]) /* rate_dn[NCHU] */
{
	DEBUG_ENTRY( "UpdatePot()" );

	ASSERT( nd < gv.bin.size() );
	ASSERT( Zlo >= gv.bin[nd].LowestZg );
	ASSERT( stride >= 1 );

	/* >>chng 00 jul 17, use description of Weingartner & Draine, 2001 */
	/* >>chng 01 mar 21, assume that grain charge is distributed in two states just below and
	 * above the average charge. */
	/* >>chng 01 may 07, this routine now completely supports the hybrid grain
	 * charge model, and the average charge state is not used anywhere anymore, PvH */
	/* >>chng 01 may 30, reorganize code such that all relevant data can be copied
	 * when a valid set of data is available from a previous call, this saves CPU time, PvH */
	/* >>chng 04 jan 17, reorganized code to use pointers to the charge bins, PvH */

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, " %ld/%ld", Zlo, stride );
	}

	if( gv.bin[nd].nfill < rfield.nPositive )
	{
		InitBinAugerData( nd, gv.bin[nd].nfill, rfield.nPositive );
		gv.bin[nd].nfill = rfield.nPositive;
	}

	for( long nz=0; nz < gv.bin[nd].nChrg; nz++ )
	{
		long Zg = Zlo + nz*stride;

		/* check if charge state is already present */
		long ind = NCHS-1;
		for( long zz=0; zz < NCHS-1; zz++ )
		{
			if( gv.bin[nd].chrg(zz).DustZ == Zg )
			{
				ind = zz;
				break;
			}
		}

		/* in the code below the old charge bins are shifted to the back in order to assure
		 * that the most recently used ones are upfront; they are more likely to be reused */
		long ptr = gv.bin[nd].ichrg[ind];

		/* make room to swap in charge state */
		for( long zz=ind-1; zz >= nz; zz-- )
			gv.bin[nd].ichrg[zz+1] = gv.bin[nd].ichrg[zz];

		gv.bin[nd].ichrg[nz] = ptr;

		if( gv.bin[nd].chrg(nz).DustZ != Zg )
			UpdatePot1(nd,nz,Zg,0);
		else if( gv.bin[nd].chrg(nz).nfill < rfield.nPositive )
			UpdatePot1(nd,nz,Zg,gv.bin[nd].chrg(nz).nfill);

		UpdatePot2(nd,nz);

		double d[4];
		rate_up[nz] = GrainElecEmis1(nd,nz,&d[0],&d[1],&d[2],&d[3]);
		rate_dn[nz] = GrainElecRecomb1(nd,nz,&d[0],&d[1]);

		/* sanity checks */
		ASSERT( gv.bin[nd].chrg(nz).DustZ == Zg );
		ASSERT( gv.bin[nd].chrg(nz).nfill >= rfield.nPositive );
		ASSERT( rate_up[nz] >= 0. && rate_dn[nz] >= 0. );
	}

	/* determine highest energy to be considered by quantum heating routines.
	 * since the Boltzmann distribution is resolved, the upper limit has to be
	 * high enough that a negligible amount of energy is in the omitted tail */
	/* >>chng 03 jan 26, moved this code from GrainChargeTemp to UpdatePot
	 *                   since the new code depends on grain potential, HTT91.in, PvH */
	double BoltzFac = (-log(CONSERV_TOL) + 8.)*BOLTZMANN/EN1RYD;
	double HighEnergy = 0.;
	for( long nz=0; nz < gv.bin[nd].nChrg; nz++ )
	{
		/* >>chng 04 jan 21, changed phycon.te -> MAX2(phycon.te,gv.bin[nd].tedust), PvH */
		HighEnergy = MAX2(HighEnergy,
		  MAX2(gv.bin[nd].chrg(nz).ThresInfInc,0.) + BoltzFac*MAX2(phycon.te,gv.bin[nd].tedust));
	}
	HighEnergy = min(HighEnergy,rfield.anu(rfield.nflux));
	gv.bin[nd].qnflux2 = ipoint(HighEnergy);
	gv.bin[nd].qnflux = max(rfield.nPositive,gv.bin[nd].qnflux2);
	gv.bin[nd].qnflux = min(gv.bin[nd].qnflux,rfield.nflux);
	return;
}


/* calculate charge state populations */
STATIC void GetFracPop(size_t nd,
		       long Zlo,
		       /*@in@*/ double rate_up[], /* rate_up[NCHU] */
		       /*@in@*/ double rate_dn[], /* rate_dn[NCHU] */
		       /*@out@*/ long *newZlo)
{
	DEBUG_ENTRY( "GetFracPop()" );

	ASSERT( nd < gv.bin.size() );
	ASSERT( Zlo >= gv.bin[nd].LowestZg );

	bool lgRedo;
	double netloss[2], pop[2][NCHU-1];

	/* solve level populations for levels 0..nChrg-2 (i == 0) and
	 * levels 1..nChrg-1 (i == 1), and determine net electron loss rate
	 * for each of those subsystems. Next we demand that
	 *          netloss[1] <= 0 <= netloss[0]
	 * and determine FracPop by linearly adding the subsystems such that
	 *     0 == frac*netloss[0] + (1-frac)*netloss[1]
	 * this assures that all charge state populations are positive */
	do
	{
		for( long i=0; i < 2; i++ )
		{
			double sum = pop[i][0] = 1.;
			long nlim = gv.bin[nd].nChrg-1;
			for( long j=1; j < nlim; j++ )
			{
				long nz = i + j;
				if( rate_dn[nz] > 10.*rate_up[nz-1]/sqrt(DBL_MAX) )
				{
					pop[i][j] = pop[i][j-1]*rate_up[nz-1]/rate_dn[nz];
					sum += pop[i][j];
				}
				else
				{
					for( long k=0; k < j; k++ )
					{
						pop[i][k] = 0.;
					}
					pop[i][j] = 1.;
					sum = pop[i][j];
				}
				/* guard against overflow */
				if( pop[i][j] > sqrt(DBL_MAX) )
				{
					for( long k=0; k <= j; k++ )
					{
						pop[i][k] /= DBL_MAX/10.;
					}
					sum /= DBL_MAX/10.;
				}
			}
			netloss[i] = 0.;
			for( long j=0; j < nlim; j++ )
			{
				long nz = i + j;
				pop[i][j] /= sum;
				netloss[i] += pop[i][j]*(rate_up[nz] - rate_dn[nz]);
			}
		}

		/* ascertain that the choice of Zlo was correct, this is to ensure positive
		 * level populations and continuous emission and recombination rates */
		if( netloss[0]*netloss[1] > 0. )
			*newZlo = ( netloss[1] > 0. ) ? Zlo + 1 : Zlo - 1;
		else
			*newZlo = Zlo;

		/* >>chng 04 feb 15, add protection for roundoff error when nChrg is much too large;
		 * netloss[0/1] can be almost zero, but may have arbitrary sign due to roundoff error;
		 * this can lead to a spurious lowering of Zlo below LowestZg, orion_pdr10.in,
		 * since newZlo cannot be lowered, lower nChrg instead and recalculate, PvH */
		/* >>chng 04 feb 15, also lower nChrg if population of highest charge state is marginal */
		if( gv.bin[nd].nChrg > 2 &&
		    ( *newZlo < gv.bin[nd].LowestZg ||
		    ( *newZlo == Zlo && pop[1][gv.bin[nd].nChrg-2] < DBL_EPSILON ) ) )
		{
			gv.bin[nd].nChrg--;
			lgRedo = true;
		}
		else
		{
			lgRedo = false;
		}

#		if 0
		printf( " fnzone %.2f nd %ld Zlo %ld newZlo %ld netloss %.4e %.4e nChrg %ld lgRedo %d\n",
			fnzone, nd, Zlo, *newZlo, netloss[0], netloss[1], gv.bin[nd].nChrg, lgRedo );
#		endif
	}
	while( lgRedo );

	if( *newZlo < gv.bin[nd].LowestZg )
	{
		fprintf( ioQQQ, " could not converge charge state populations for %s\n", gv.bin[nd].chDstLab );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( *newZlo == Zlo )
	{
		/* frac1 == 1-frac0, but calculating it like this avoids cancellation errors */
		double frac0 = netloss[1]/(netloss[1]-netloss[0]);
		double frac1 = -netloss[0]/(netloss[1]-netloss[0]);

		gv.bin[nd].chrg(0).FracPop = frac0*pop[0][0];
		gv.bin[nd].chrg(gv.bin[nd].nChrg-1).FracPop = frac1*pop[1][gv.bin[nd].nChrg-2];
		for( long nz=1; nz < gv.bin[nd].nChrg-1; nz++ )
		{
			gv.bin[nd].chrg(nz).FracPop = frac0*pop[0][nz] + frac1*pop[1][nz-1];
		}

#		ifndef NDEBUG
		double test1 = 0., test2 = 0., test3 = 0.;
		for( long nz=0; nz < gv.bin[nd].nChrg; nz++ )
		{
			ASSERT( gv.bin[nd].chrg(nz).FracPop >= 0. );
			test1 += gv.bin[nd].chrg(nz).FracPop;
			test2 += gv.bin[nd].chrg(nz).FracPop*rate_up[nz];
			test3 += gv.bin[nd].chrg(nz).FracPop*(rate_up[nz]-rate_dn[nz]);
		}
		double x1 = fabs(test1-1.);
		double x2 = fabs( safe_div(test3, test2, 0.) );
		ASSERT( MAX2(x1,x2) < 10.*sqrt((double)gv.bin[nd].nChrg)*DBL_EPSILON );
#		endif
	}
	return;
}


/* this routine updates all quantities that depend only on grain charge and radius
 *
 * NB NB - All global data in grain.c and grainvar.h that are charge dependent should
 *         be calculated here or in UpdatePot2 (except gv.bin[nd].chrg(nz).FracPop
 *         which is special).
 *
 * NB NB - the code assumes that the data calculated here remain valid throughout
 *         the model, i.e. do NOT depend on grain temperature, etc.
 */
STATIC void UpdatePot1(size_t nd,
		       long nz,
		       long Zg,
		       long ipStart)
{
	DEBUG_ENTRY( "UpdatePot1()" );

	/* constants in the expression for the photodetachment cross section
	 * taken from 
	 * >>refer	grain	physics	Weingartner & Draine, ApJS, 2001, 134, 263 */
	const double INV_DELTA_E = EVRYD/3.;
	const double CS_PDT = 1.2e-17;

	/* >>chng 04 jan 23, replaced rfield.nflux with rfield.nflux_with_check so
	 *                   that data remains valid throughout the model, PvH */
	/* >>chng 04 jan 26, start using ipStart to solve the validity problem 
	 *      ipStart == 0 means that a full initialization needs to be done
	 *      ipStart > 0  means that rfield.nflux was updated and yhat(_primary), ehat,
	 *                   cs_pdt, fac1, and fac2 need to be reallocated, PvH */
	if( ipStart == 0 )
	{
		gv.bin[nd].chrg(nz).DustZ = Zg;

		/* invalidate eta and xi storage */
		memset( gv.bin[nd].chrg(nz).eta, 0, (LIMELM+2)*sizeof(double) );
		memset( gv.bin[nd].chrg(nz).xi, 0, (LIMELM+2)*sizeof(double) );

		GetPotValues(nd,Zg,&gv.bin[nd].chrg(nz).ThresInf,&gv.bin[nd].chrg(nz).ThresInfVal,
			     &gv.bin[nd].chrg(nz).ThresSurf,&gv.bin[nd].chrg(nz).ThresSurfVal,
			     &gv.bin[nd].chrg(nz).PotSurf,&gv.bin[nd].chrg(nz).Emin,INCL_TUNNEL);

		/* >>chng 01 may 09, do not use tunneling corrections for incoming electrons */
		/* >>chng 01 nov 25, add gv.bin[nd].chrg(nz).ThresInfInc, PvH */
		double d[2];
		GetPotValues(nd,Zg-1,&gv.bin[nd].chrg(nz).ThresInfInc,&d[0],&gv.bin[nd].chrg(nz).ThresSurfInc,
			     &d[1],&gv.bin[nd].chrg(nz).PotSurfInc,&gv.bin[nd].chrg(nz).EminInc,NO_TUNNEL);

		gv.bin[nd].chrg(nz).ipThresInfVal = rfield.ithreshC( gv.bin[nd].chrg(nz).ThresInfVal );
	}

	ASSERT( gv.bin[nd].chrg(nz).DustZ == Zg );

	long ipLo = gv.bin[nd].chrg(nz).ipThresInfVal;

	/* remember how far the yhat(_primary), ehat, cs_pdt, fac1, and fac2 arrays were filled in */
	gv.bin[nd].chrg(nz).nfill = rfield.nPositive;
	long nfill = gv.bin[nd].chrg(nz).nfill;

	/* >>chng 04 feb 07, only allocate arrays from ipLo to nfill to save memory, PvH */
	long len = rfield.nflux_with_check - ipLo;
	if( ipStart == 0 )
	{
		gv.bin[nd].chrg(nz).yhat.reserve( len );
		gv.bin[nd].chrg(nz).yhat_primary.reserve( len );
		gv.bin[nd].chrg(nz).ehat.reserve( len );
		gv.bin[nd].chrg(nz).yhat.alloc( ipLo, nfill );
		gv.bin[nd].chrg(nz).yhat_primary.alloc( ipLo, nfill );
		gv.bin[nd].chrg(nz).ehat.alloc( ipLo, nfill );
	}
	else
	{
		gv.bin[nd].chrg(nz).yhat.realloc( nfill );
		gv.bin[nd].chrg(nz).yhat_primary.realloc( nfill );
		gv.bin[nd].chrg(nz).ehat.realloc( nfill );
	}

	double GrainPot = chrg2pot(Zg,nd);
	const double *anu = rfield.anuptr();

	if( nfill > ipLo )
	{
		double DustWorkFcn = gv.bin[nd].DustWorkFcn;
		double Elo = -gv.bin[nd].chrg(nz).PotSurf;
		double ThresInfVal = gv.bin[nd].chrg(nz).ThresInfVal;
		double Emin = gv.bin[nd].chrg(nz).Emin;
		double Wcorr = ThresInfVal + Emin - GrainPot;

		/** \todo xray - secondaries from incident electrons still need to be added in */
		/** \todo xray - primary, secondary, auger electrons need to be added into suprathermals */

		ShellData *sptr = &gv.bin[nd].sd[0];

		ASSERT( sptr->ipLo <= ipLo );

		long ipBegin = max(ipLo,ipStart);
		avx_ptr<realnum> yzero(ipBegin, nfill);
		y0b(nd, nz, yzero.ptr0(), ipBegin, nfill);

		/* calculate yield for band structure */
		avx_ptr<double> Ehi(ipBegin, nfill), Ehi_band(ipBegin, nfill), Eel(ipBegin, nfill);

		for( long i=ipBegin; i < nfill; i++ )
		{
			Ehi[i] = anu[i] - ThresInfVal - Emin;
			Ehi_band[i] = Ehi[i];
			Eel[i] = anu[i] - DustWorkFcn;
		}

		avx_ptr<realnum> yp(ipBegin, nfill), ya(ipBegin, nfill), ys(ipBegin, nfill), eyp(ipBegin, nfill),
			eya(ipBegin, nfill), eys(ipBegin, nfill), Yp1(ipBegin, nfill), Ys1(ipBegin, nfill),
			Ehp(ipBegin, nfill), Ehs(ipBegin, nfill);
	     
		Yfunc(nd, nz, yzero.ptr0(), sptr->y01.ptr0(), sptr->p.ptr0(), Elo, Ehi.ptr0(), Eel.ptr0(),
		      Yp1.ptr0(), Ys1.ptr0(), Ehp.ptr0(), Ehs.ptr0(), ipBegin, nfill);

		if( nfill > ipBegin )
		{
			memset( ya.data(), 0, size_t((nfill-ipBegin)*sizeof(ya[ipBegin])) );
			memset( eya.data(), 0, size_t((nfill-ipBegin)*sizeof(eya[ipBegin])) );
		}

		// write this as two separate loops so that they get vectorized
		for( long i=ipBegin; i < nfill; i++ )
		{
			yp[i] = Yp1[i];
			eyp[i] = Yp1[i]*Ehp[i];
		}
		for( long i=ipBegin; i < nfill; i++ )
		{
			ys[i] = Ys1[i];
			eys[i] = Ys1[i]*Ehs[i];
		}

		/* add in yields for inner shell ionization */
		unsigned int nsmax = gv.bin[nd].sd.size();
		for( unsigned int ns=1; ns < nsmax; ns++ )
		{
			sptr = &gv.bin[nd].sd[ns];

			long ipBeginShell = max(ipBegin, sptr->ipLo);
			double ionPot = sptr->ionPot;
			for( long i=ipBeginShell; i < nfill; i++ )
			{
				Ehi[i] = Ehi_band[i] + Wcorr - ionPot;
				Eel[i] = anu[i] - ionPot;
			}

			Yfunc(nd, nz, sptr->y01.ptr0(), NULL, sptr->p.ptr0(), Elo, Ehi.ptr0(), Eel.ptr0(),
			      Yp1.ptr0(), Ys1.ptr0(), Ehp.ptr0(), Ehs.ptr0(), ipBeginShell, nfill);

			// write this as two separate loops so that they get vectorized
			for( long i=ipBeginShell; i < nfill; i++ )
			{
				yp[i] += Yp1[i];
				eyp[i] += Yp1[i]*Ehp[i];
			}
			for( long i=ipBeginShell; i < nfill; i++ )
			{
				ys[i] += Ys1[i];
				eys[i] += Ys1[i]*Ehs[i];
			}

			/* add in Auger yields */
			long nmax = sptr->nData;
			for( long n=0; n < nmax; n++ )
			{
				avx_ptr<realnum> max(ipBeginShell, nfill);
				double AvNr = sptr->AvNr[n];
				double Ener = sptr->Ener[n];
				double Ehic = Ener - GrainPot;
				double Eelc = Ener;
				for( long i=ipBeginShell; i < nfill; i++ )
					max[i] = AvNr*sptr->p[i];

				Yfunc(nd, nz, sptr->y01A[n].ptr0(), max.ptr0(), Elo, Ehic, Eelc,
				      Yp1.ptr0(), Ys1.ptr0(), Ehp.ptr0(), Ehs.ptr0(), ipBeginShell, nfill);

				// write this as two separate loops so that they get vectorized
				for( long i=ipBeginShell; i < nfill; i++ )
				{
					ya[i] += Yp1[i];
					eya[i] += Yp1[i]*Ehp[i];
				}
				for( long i=ipBeginShell; i < nfill; i++ )
				{
					ys[i] += Ys1[i];
					eys[i] += Ys1[i]*Ehs[i];
				}
			}
		}

		/* add up primary, Auger, and secondary yields */
		for( long i=ipBegin; i < nfill; i++ )
		{
			gv.bin[nd].chrg(nz).yhat[i] = yp[i] + ya[i] + ys[i];
			gv.bin[nd].chrg(nz).yhat_primary[i] = min(yp[i],1.f);
			// if gv.bin[nd].chrg(nz).ThresInfVal and rfield.anu(ipLo) are extremely
			// close, gv.bin[nd].chrg(nz).yhat[i] and ehat[i] can underflow to zero
			gv.bin[nd].chrg(nz).ehat[i] = ( gv.bin[nd].chrg(nz).yhat[i] > 0.f ) ?
				(eyp[i] + eya[i] + eys[i])/gv.bin[nd].chrg(nz).yhat[i] : 0.f;

			ASSERT( yp[i] <= 1.00001f );
			ASSERT( gv.bin[nd].chrg(nz).ehat[i] >= 0.f );
		}
	}

	if( ipStart == 0 )
	{
		/* >>chng 01 jan 08, ThresInf[nz] and ThresInfVal[nz] may become zero in
		 * initial stages of grain potential search, PvH */
		/* >>chng 01 oct 10, use bisection search to find ipThresInf, ipThresInfVal. On C scale now */
		gv.bin[nd].chrg(nz).ipThresInf = rfield.ithreshC( gv.bin[nd].chrg(nz).ThresInf );
	}

	ipLo = gv.bin[nd].chrg(nz).ipThresInf;

	len = rfield.nflux_with_check - ipLo;

	if( Zg <= -1 )
	{
		/* >>chng 04 feb 07, only allocate arrays from ipLo to nfill to save memory, PvH */
		if( ipStart == 0 )
		{
			gv.bin[nd].chrg(nz).cs_pdt.reserve( len );
			gv.bin[nd].chrg(nz).cs_pdt.alloc( ipLo, rfield.nflux );

			double c1 = -CS_PDT*(double)Zg;
			double ThresInf = gv.bin[nd].chrg(nz).ThresInf;
			double cnv_GR_pH = gv.bin[nd].cnv_GR_pH;

			for( long i=ipLo; i < rfield.nflux; i++ )
			{
				double x = (anu[i] - ThresInf)*INV_DELTA_E;
				double cs = c1*x/POW2(1.+(1./3.)*POW2(x));
				/* this cross section must be at default depletion for consistency
				 * with dstab1, hence no factor dstAbund should be applied here */
				gv.bin[nd].chrg(nz).cs_pdt[i] = MAX2(cs,0.)*cnv_GR_pH;
			}
		}
	}

	/* >>chng 04 feb 07, added fac1, fac2 to optimize loop for photoelectric heating, PvH */
	if( ipStart == 0 )
	{
		gv.bin[nd].chrg(nz).fac1.reserve( len );
		gv.bin[nd].chrg(nz).fac2.reserve( len );
		gv.bin[nd].chrg(nz).fac1.alloc( ipLo, nfill );
		gv.bin[nd].chrg(nz).fac2.alloc( ipLo, nfill );
	}
	else
	{
		gv.bin[nd].chrg(nz).fac1.realloc( nfill );
		gv.bin[nd].chrg(nz).fac2.realloc( nfill );
	}

	if( nfill > ipLo )
	{
		for( long i=max(ipLo,ipStart); i < nfill; i++ )
		{
			double cs1,cs2,cs_tot,cool1,cool2,ehat1,ehat2;

			/* >>chng 04 jan 25, created separate routine PE_init, PvH */
			PE_init(nd,nz,i,&cs1,&cs2,&cs_tot,&cool1,&cool2,&ehat1,&ehat2);

			gv.bin[nd].chrg(nz).fac1[i] = cs_tot*anu[i]-cs1*cool1-cs2*cool2;
			gv.bin[nd].chrg(nz).fac2[i] = cs1*ehat1 + cs2*ehat2;

			ASSERT( gv.bin[nd].chrg(nz).fac1[i] >= 0. && gv.bin[nd].chrg(nz).fac2[i] >= 0. );
		}
	}

	if( ipStart == 0 )
	{
		/* >>chng 00 jul 05, determine ionization stage Z0 the ion recombines to */
		/* >>chng 04 jan 20, use all stages here so that result remains valid throughout the model */
		UpdateRecomZ0(nd,nz);
	}

	/* invalidate the remaining fields */
	gv.bin[nd].chrg(nz).FracPop = -DBL_MAX;

	gv.bin[nd].chrg(nz).RSum1 = -DBL_MAX;
	gv.bin[nd].chrg(nz).RSum2 = -DBL_MAX;
	gv.bin[nd].chrg(nz).ESum1a = -DBL_MAX;
	gv.bin[nd].chrg(nz).ESum1b = -DBL_MAX;
	gv.bin[nd].chrg(nz).ESum2 = -DBL_MAX;

	gv.bin[nd].chrg(nz).tedust = 1.f;

	gv.bin[nd].chrg(nz).hcon1 = -DBL_MAX;
	gv.bin[nd].chrg(nz).hots1 = -DBL_MAX;
	gv.bin[nd].chrg(nz).bolflux1 = -DBL_MAX;
	gv.bin[nd].chrg(nz).pe1 = -DBL_MAX;

	gv.bin[nd].chrg(nz).BolFlux = -DBL_MAX;
	gv.bin[nd].chrg(nz).GrainHeat = -DBL_MAX;
	gv.bin[nd].chrg(nz).GrainHeatColl = -DBL_MAX;
	gv.bin[nd].chrg(nz).GasHeatPhotoEl = -DBL_MAX;
	gv.bin[nd].chrg(nz).GasHeatTherm = -DBL_MAX;
	gv.bin[nd].chrg(nz).GrainCoolTherm = -DBL_MAX;
	gv.bin[nd].chrg(nz).ChemEnIon = -DBL_MAX;
	gv.bin[nd].chrg(nz).ChemEnH2 = -DBL_MAX;

	gv.bin[nd].chrg(nz).HeatingRate2 = -DBL_MAX;

	/* sanity check */
	ASSERT( gv.bin[nd].chrg(nz).ipThresInf <= gv.bin[nd].chrg(nz).ipThresInfVal );
	return;
}


/* this routine updates all quantities that depend on grain charge, radius and temperature
 *
 * NB NB - All global data in grain.c and grainvar.h that are charge dependent should
 *         be calculated here or in UpdatePot1 (except gv.bin[nd].chrg(nz).FracPop
 *         which is special).
 *
 * NB NB - the code assumes that the data calculated here may vary throughout the model,
 *         e.g. because of a dependence on grain temperature
 */
STATIC void UpdatePot2(size_t nd,
		       long nz)
{
	DEBUG_ENTRY( "UpdatePot2()" );

	/* >>chng 00 jun 19, add in loss rate due to thermionic emission of electrons, PvH */
	double ThermExp = gv.bin[nd].chrg(nz).ThresInf*TE1RYD/gv.bin[nd].tedust;
	/* ThermExp is guaranteed to be >= 0. */
	gv.bin[nd].chrg(nz).ThermRate = THERMCONST*gv.bin[nd].ThermEff*POW2(gv.bin[nd].tedust)*exp(-ThermExp);
#	if defined( WD_TEST2 ) || defined( IGNORE_THERMIONIC )
	gv.bin[nd].chrg(nz).ThermRate = 0.;
#	endif
	return;
}


/* Helper function to calculate primary and secondary yields and the average electron energy at infinity */
STATIC void Yfunc(long nd,
		  long nz,
		  const realnum y0[],
		  const realnum y1[],
		  const realnum maxval[],
		  double Elo,
		  const double Ehi[],
		  const double Eel[],
		  /*@out@*/ realnum Yp[],
		  /*@out@*/ realnum Ys[],
		  /*@out@*/ realnum Ehp[],
		  /*@out@*/ realnum Ehs[],
		  long ilo,
		  long ihi)
{
	DEBUG_ENTRY( "Yfunc()" );

	if( ihi <= ilo )
		return;

	long Zg = gv.bin[nd].chrg(nz).DustZ;

	double eps;
	pe_type pcase = gv.which_pe[gv.bin[nd].matType];
	if( pcase == PE_CAR )
		eps = 117./EVRYD;
	else if( pcase == PE_SIL )
		eps = 155./EVRYD;
	else
	{
		fprintf( ioQQQ, " Yfunc: unknown type for PE effect: %d\n" , pcase );
		cdEXIT(EXIT_FAILURE);
	}

	avx_ptr<realnum> y2pr(ilo, ihi), y2sec(ilo, ihi), y01cap(ilo, ihi);

	y2pa( Elo, Ehi, Zg, Ehp, y2pr.ptr0(), ilo, ihi );
	y2s( Elo, Ehi, Zg, y2pr.ptr0(), Ehs, y2sec.ptr0(), ilo, ihi );

	if( y1 == NULL )
	{
		for( long i=ilo; i < ihi; ++i )
			y01cap[i] = min(y0[i],maxval[i]);
	}
	else
	{
		for( long i=ilo; i < ihi; ++i )
			y01cap[i] = min(y0[i]*y1[i],maxval[i]);
	}
		
	for( long i=ilo; i < ihi; ++i )
	{
		ASSERT( Ehi[i] >= Elo );

		if( y2pr[i] > 0.f )
		{
			Yp[i] = y2pr[i]*y01cap[i];

			/* this is Eq. 18 of WDB06 */
			/* Eel may be negative near threshold -> set yield to zero */
			realnum f3 = max(Eel[i],0.)/(eps*elec_esc_length(Eel[i],nd)*gv.bin[nd].eyc);
			Ys[i] = y2sec[i]*f3*y01cap[i];
		}
		else
		{
			Yp[i] = 0.f;
			Ys[i] = 0.f;
			Ehp[i] = 0.f;
			Ehs[i] = 0.f;
		}
	}
}


// overloaded version of Yfunc exploiting the fact that Ehi and Eel are independent of frequency
STATIC void Yfunc(long nd,
		  long nz,
		  const realnum y01[],
		  const realnum maxval[],
		  double Elo,
		  double Ehi,
		  double Eel,
		  /*@out@*/ realnum Yp[],
		  /*@out@*/ realnum Ys[],
		  /*@out@*/ realnum Ehp[],
		  /*@out@*/ realnum Ehs[],
		  long ilo,
		  long ihi)
{
	DEBUG_ENTRY( "Yfunc()" );

	if( ihi <= ilo )
		return;

	ASSERT( Ehi >= Elo );

	// since both Elo and Ehi are constant, there is no frequency dependence of the results of
	// y2pa() and y2s() -> we only evaluate one frequency and copy this to the rest of the array.

	avx_ptr<realnum> y2pr(ilo, ilo+1), y2sec(ilo, ilo+1);
	avx_ptr<double> Ehiloc(ilo, ilo+1);

	Ehiloc[ilo] = Ehi;
	long Zg = gv.bin[nd].chrg(nz).DustZ;

	y2pa( Elo, Ehiloc.ptr0(), Zg, Ehp, y2pr.ptr0(), ilo, ilo+1 );

	if( y2pr[ilo] > 0.f )
	{
		y2s( Elo, Ehiloc.ptr0(), Zg, y2pr.ptr0(), Ehs, y2sec.ptr0(), ilo, ilo+1 );

		for( long i=ilo+1; i < ihi; ++i )
		{
			Ehp[i] = Ehp[ilo];
			Ehs[i] = Ehs[ilo];
		}

		double eps;
		pe_type pcase = gv.which_pe[gv.bin[nd].matType];
		if( pcase == PE_CAR )
			eps = 117./EVRYD;
		else if( pcase == PE_SIL )
			eps = 155./EVRYD;
		else
		{
			fprintf( ioQQQ, " Yfunc: unknown type for PE effect: %d\n" , pcase );
			cdEXIT(EXIT_FAILURE);
		}

		/* this is Eq. 18 of WDB06 */
		/* Eel may be negative near threshold -> set yield to zero */
		realnum f3 = max(Eel,0.)/(eps*elec_esc_length(Eel,nd)*gv.bin[nd].eyc);

		for( long i=ilo; i < ihi; ++i )
		{
			realnum y01cap = min(y01[i],maxval[i]);
			Yp[i] = y2pr[ilo]*y01cap;
			Ys[i] = y2sec[ilo]*f3*y01cap;
		}
	}
	else
	{
		memset( Yp+ilo, 0, size_t((ihi-ilo)*sizeof(Yp[ilo])) );
		memset( Ys+ilo, 0, size_t((ihi-ilo)*sizeof(Ys[ilo])) );
		memset( Ehp+ilo, 0, size_t((ihi-ilo)*sizeof(Ehp[ilo])) );
		memset( Ehs+ilo, 0, size_t((ihi-ilo)*sizeof(Ehs[ilo])) );
	}
}


/* This calculates the y0 function for band electrons (Sect. 4.1.3/4.1.4 of WDB06) */
STATIC void y0b(size_t nd,
		long nz,
		/*@out@*/realnum yzero[],
		long ilo,
		long ihi)
{
	DEBUG_ENTRY( "y0b()" );

	if( ihi <= ilo )
		return;

	if( gv.lgWD01 )
		y0b01( nd, nz, yzero, ilo, ihi );
	else
	{
		realnum Ethres_low = 20./EVRYD;
		realnum Ethres_high = 50./EVRYD;
		long ip20 = rfield.ithreshC(Ethres_low);
		long ip50 = rfield.ithreshC(Ethres_high);
		long ipr1a = ilo;
		//long ipr1b = min(ip20,ihi);
		long ipr2a = max(ilo,ip20);
		long ipr2b = min(ip50,ihi);
		long ipr3a = max(ilo,ip50);
		long ipr3b = ihi;

		y0b01( nd, nz, yzero, ipr1a, ipr2b );
		avx_ptr<realnum> arg(ipr2a, ipr2b), val(ipr2a, ipr2b), val2(ipr2a, ipr2b);
		for( int i=ipr2a; i < ipr2b; i++ )
			arg[i] = realnum(rfield.anu(i))/Ethres_low;
		vlog( arg.ptr0(), val.ptr0(), ipr2a, ipr2b );
		for( int i=ipr2a; i < ipr2b; i++ )
			arg[i] = gv.bin[nd].y0b06[i]/yzero[i];
		vlog( arg.ptr0(), val2.ptr0(), ipr2a, ipr2b );
		for( int i=ipr2a; i < ipr2b; i++ )
			/* constant 1.09135666... is 1./log(50./20.) */
			arg[i] = val2[i]*val[i]*1.0913566679f;
		vexp( arg.ptr0(), val.ptr0(), ipr2a, ipr2b );
		for( int i=ipr2a; i < ipr2b; i++ )
			yzero[i] *= val[i];
		for( int i=ipr3a; i < ipr3b; i++ )
			yzero[i] = gv.bin[nd].y0b06[i];
	}
	// result may underflow to zero if anu(ilo) is very close to threshold
	ASSERT( yzero[ilo] >= 0.f );
	for( int i=ilo+1; i < ihi; i++ )
		ASSERT( yzero[i] > 0.f );
}


/* This calculates the y0 function for band electrons (Eq. 16 of WD01) */
STATIC void y0b01(size_t nd,
		  long nz,
		  /*@out@*/realnum yzero[],
		  long ilo,
		  long ihi)
{
	DEBUG_ENTRY( "y0b01()" );

	pe_type pcase = gv.which_pe[gv.bin[nd].matType];
	switch( pcase )
	{
	case PE_CAR:
		/* >>refer	grain	physics	Bakes & Tielens, 1994, ApJ, 427, 822 */
		for( int i=ilo; i < ihi; i++ )
		{
			double xv = max((rfield.anu(i)-gv.bin[nd].chrg(nz).ThresSurfVal)/gv.bin[nd].DustWorkFcn,0.);
			double xv2 = xv*xv;
			xv = xv2*xv2*xv;
			yzero[i] = realnum(xv/((1./9.e-3) + (3.7e-2/9.e-3)*xv));
		}
		break;
	case PE_SIL:
		/* >>refer	grain	physics	Weingartner & Draine, 2001 */
		for( int i=ilo; i < ihi; i++ )
		{
			double xv = max((rfield.anu(i)-gv.bin[nd].chrg(nz).ThresSurfVal)/gv.bin[nd].DustWorkFcn,0.);
			yzero[i] = realnum(xv/(2.+10.*xv));
		}
		break;
	default:
		fprintf( ioQQQ, " y0b01: unknown type for PE effect: %d\n" , pcase );
		cdEXIT(EXIT_FAILURE);
	}
}


/* This calculates the y0 function for primary/secondary and Auger electrons (Eq. 9 of WDB06) */
STATIC double y0psa(size_t nd,
		    long ns,    /* shell number */
		    long i,     /* incident photon energy is anu[i] */
		    double Eel) /* emitted electron energy */
{
	DEBUG_ENTRY( "y0psa()" );

	ASSERT( i >= gv.bin[nd].sd[ns].ipLo );

	/* this is l_e/l_a */
	double leola = elec_esc_length(Eel,nd)*gv.bin[nd].inv_att_len[i];

	ASSERT( leola > 0. );

	/* this is Eq. 9 of WDB06 */
	double yzero;
	if( leola < 1.e4 )
		yzero = gv.bin[nd].sd[ns].p[i]*leola*(1. - leola*log(1.+1./leola));
	else
	{
		double x = 1./leola;
		yzero = gv.bin[nd].sd[ns].p[i]*(((-1./5.*x+1./4.)*x-1./3.)*x+1./2.);
	}
	ASSERT( yzero > 0. );
	return yzero;
}


/* This calculates the y1 function for primary/secondary and Auger electrons (Eq. 6 of WDB06) */
STATIC double y1psa(size_t nd,
		    long i,     /* incident photon energy is anu[i] */
		    double Eel) /* emitted electron energy */
{
	DEBUG_ENTRY( "y1psa()" );

	double bf, beta = gv.bin[nd].AvRadius*gv.bin[nd].inv_att_len[i];
	if( beta > 1.e-4 ) 
		bf = pow2(beta) - 2.*beta + 2. - 2.*exp(-beta);
	else 
		bf = ((1./60.*beta - 1./12.)*beta + 1./3.)*pow3(beta);

	double af, alpha = beta + gv.bin[nd].AvRadius/elec_esc_length(Eel,nd);
	if( alpha > 1.e-4 ) 
		af = pow2(alpha) - 2.*alpha + 2. - 2.*exp(-alpha);
	else 
		af = ((1./60.*alpha - 1./12.)*alpha + 1./3.)*pow3(alpha);

	double yone = pow2(beta/alpha)*af/bf;

	ASSERT( yone > 0. );
	return yone;
}


/* This calculates the y2 function for primary and Auger electrons (Eq. 8 of WDB06) */
inline void y2pa(double Elo,
		 const double Ehi[],
		 long Zg,
		 /*@out@*/realnum Ehp[],
		 /*@out@*/realnum y2pr[],
		 long ilo,
		 long ihi)
{
	DEBUG_ENTRY( "y2pa()" );

	if( Zg > -1 )
	{
		for( long i=ilo; i < ihi; ++i )
		{
			if( Ehi[i] > 0. )
			{
				double x = Elo/Ehi[i];
				Ehp[i] = 0.5f*realnum(Ehi[i]*(1.-2.*x)/(1.-3.*x));
				// use Taylor expansion for small arguments to avoid blowing assert
				y2pr[i] = ( abs(x) > 1e-4 ) ?
					realnum((1.-3.*x)/pow3(1.-x)) : realnum(1. - (3. + 8.*x)*x*x);
				ASSERT( Ehp[i] > 0.f && Ehp[i] <= realnum(Ehi[i]) && y2pr[i] > 0.f && y2pr[i] <= 1.f );
			}
			else
			{
				Ehp[i] = 0.f;
				y2pr[i] = 0.f;
			}
		}
	}
	else
	{
		for( long i=ilo; i < ihi; ++i )
		{
			if( Ehi[i] > Elo )
			{
				Ehp[i] = 0.5f*realnum(Elo+Ehi[i]);
				y2pr[i] = 1.f;
				ASSERT( Ehp[i] >= realnum(Elo) && Ehp[i] <= realnum(Ehi[i]) );
			}
			else
			{
				Ehp[i] = 0.f;
				y2pr[i] = 0.f;
			}
		}
	}
}


/* This calculates the y2 function for secondary electrons (Eqs. 20-21 of WDB06) */
inline void y2s(double Elo,
		const double Ehi[],
		long Zg,
		const realnum y2pr[],
		/*@out@*/realnum Ehs[],
		/*@out@*/realnum y2sec[],
		long ilo,
		long ihi)
{
	DEBUG_ENTRY( "y2s()" );

	if( ihi <= ilo )
		return;

	avx_ptr<realnum> arg(ilo, ihi), val(ilo, ihi);
	avx_ptr<double> arg2(ilo, ihi), val2(ilo, ihi), alpha(ilo, ihi);
	avx_ptr<bool> lgVecFun(ilo, ihi);
	memset( arg.data(), 0, size_t((ihi-ilo)*sizeof(arg[ilo])) );
	memset( arg2.data(), 0, size_t((ihi-ilo)*sizeof(arg2[ilo])) );
	memset( lgVecFun.data(), 0, size_t((ihi-ilo)*sizeof(lgVecFun[ilo])) );
	long iveclo = -1, ivechi = -1;

	double yl = Elo*INV_ETILDE;
	double sR0 = (1. + yl*yl);
	double sqR0 = sqrt(sR0);

	if( Zg > -1 )
	{
		double asinh_yl = asinh(yl);
		avx_ptr<double> E0(ilo, ihi), N0(ilo, ihi);

		for( long i=ilo; i < ihi; ++i )
		{
			double yh = Ehi[i]*INV_ETILDE;
			double x = yh - yl;
			arg2[i] = 1. + x*x;
		}
		vsqrt(arg2.ptr0(), val2.ptr0(), ilo, ihi);

		for( long i=ilo; i < ihi; ++i )
		{
			if( !gv.lgWD01 && Ehi[i] > 0. && y2pr[i] > 0.f )
			{
				double yh = Ehi[i]*INV_ETILDE;
				double x = yh - yl;
				if( x < 0.01 )
				{
					// use series expansions to avoid cancellation error
					double x2 = x*x, x3 = x2*x, x4 = x3*x, x5 = x4*x;
					double yh2 = yh*yh, yh3 = yh2*yh, yh4 = yh3*yh, yh5 = yh4*yh;
					double help1 = 2.*x-yh;
					double help2 = (6.*x3-15.*yh*x2+12.*yh2*x-3.*yh3)/4.;
					double help3 = (22.*x5-95.*yh*x4+164.*yh2*x3-141.*yh3*x2+60.*yh4*x-10.*yh5)/16.;
					N0[i] = yh*(help1 - help2 + help3)/x2;

					help1 = (3.*x-yh)/3.;
					help2 = (15.*x3-25.*yh*x2+15.*yh2*x-3.*yh3)/20.;
					help3 = (1155.*x5-3325.*yh*x4+4305.*yh2*x3-2961.*yh3*x2+1050.*yh4*x-150.*yh5)/
						1680.;
					E0[i] = ETILDE*yh2*(help1 - help2 + help3)/x2;
				}
				else
				{
					double sqRh = val2[i];
					alpha[i] = sqRh/(sqRh - 1.);
					if( yh/sqR0 < 0.01 )
					{
						// use series expansions to avoid cancellation error
						double z = yh*(yh - 2.*yl)/sR0;
						N0[i] = ((((7./256.*z-5./128.)*z+1./16.)*z-1./8.)*z+1./2.)*z/(sqRh-1.);

						double yl2 = yl*yl, yl3 = yl2*yl, yl4 = yl3*yl;
						double help1 = yl/2.;
						double help2 = (2.*yl2-1.)/3.;
						double help3 = (6.*yl3-9.*yl)/8.;
						double help4 = (8.*yl4-24.*yl2+3.)/10.;
						double h = yh/sR0;
						E0[i] = -alpha[i]*Ehi[i]*(((help4*h+help3)*h+help2)*h+help1)*h/sqR0;
					}
					else
					{
						N0[i] = alpha[i]*(1./sqR0 - 1./sqRh);
						arg[i] = realnum(x);
						if( iveclo < 0 )
							iveclo = i;
						ivechi = i + 1;
						lgVecFun[i] = true;
						E0[i] = alpha[i]*ETILDE*(asinh_yl - yh/sqRh);
					}
				}
				ASSERT( N0[i] > 0. && N0[i] <= 1. );
			}
		}
		vfast_asinh(arg.ptr0(), val.ptr0(), iveclo, ivechi);

		for( long i=ilo; i < ihi; ++i )
		{
			if( !gv.lgWD01 && Ehi[i] > 0. && y2pr[i] > 0.f )
			{
				if( lgVecFun[i] )
					E0[i] += alpha[i]*ETILDE*val[i];
				Ehs[i] = realnum(E0[i]/N0[i]);
				ASSERT( Ehs[i] > 0.f && Ehs[i] <= realnum(Ehi[i]) );
				y2sec[i] = realnum(N0[i]);
			}
			else
			{
				Ehs[i] = 0.f;
				y2sec[i] = 0.f;
			}
		}
	}
	else
	{
		for( long i=ilo; i < ihi; ++i )
		{
			double yh = Ehi[i]*INV_ETILDE;
			double x = yh - yl;
			arg2[i] = 1. + x*x;
		}
		vsqrt(arg2.ptr0(), val2.ptr0(), ilo, ihi);

		for( long i=ilo; i < ihi; ++i )
		{
			if( !gv.lgWD01 && Ehi[i] > Elo && y2pr[i] > 0.f )
			{
				double yh = Ehi[i]*INV_ETILDE;
				double x = yh - yl;
				double x2 = x*x;
				if( x > 0.025 )
				{
					double sqRh = val2[i];
					alpha[i] = sqRh/(sqRh - 1.);
					arg[i] = realnum(x);
					if( iveclo < 0 )
						iveclo = i;
					ivechi = i + 1;
					lgVecFun[i] = true;
					Ehs[i] = realnum(alpha[i]*ETILDE*(yl - yh/sqRh));
				}
				else
				{
					// use series expansion to avoid cancellation error
					Ehs[i] = realnum(Ehi[i] - (Ehi[i]-Elo)*((-37./840.*x2 + 1./10.)*x2 + 1./3.));
				}
			}
		}
		vfast_asinh(arg.ptr0(), val.ptr0(), iveclo, ivechi);

		for( long i=ilo; i < ihi; ++i )
		{
			if( !gv.lgWD01 && Ehi[i] > Elo && y2pr[i] > 0.f )
			{
				if( lgVecFun[i] )
					Ehs[i] += realnum(alpha[i]*ETILDE)*val[i];
				ASSERT( Ehs[i] >= realnum(Elo) && Ehs[i] <= realnum(Ehi[i]) );
				y2sec[i] = 1.f;
			}
			else
			{
				Ehs[i] = 0.f;
				y2sec[i] = 0.f;
			}
		}
	}
}


STATIC void UpdateRecomZ0(size_t nd,
			  long nz)
{
	DEBUG_ENTRY( "UpdateRecomZ0()" );

	long Zg = gv.bin[nd].chrg(nz).DustZ;

	double d[5], phi_s_up[LIMELM+1], phi_s_dn[2];
	phi_s_up[0] = gv.bin[nd].chrg(nz).ThresSurf;
	for( long i=1; i <= LIMELM; i++ )
		GetPotValues(nd,Zg+i,&d[0],&d[1],&phi_s_up[i],&d[2],&d[3],&d[4],INCL_TUNNEL);
	phi_s_dn[0] = gv.bin[nd].chrg(nz).ThresSurfInc;
	GetPotValues(nd,Zg-2,&d[0],&d[1],&phi_s_dn[1],&d[2],&d[3],&d[4],NO_TUNNEL);

	/* >>chng 01 may 09, use GrainIonColl which properly tracks step-by-step charge changes */
	for( long nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			for( long ion=0; ion <= nelem+1; ion++ )
			{
				GrainIonColl(nd,nz,nelem,ion,phi_s_up,phi_s_dn,
					     &gv.bin[nd].chrg(nz).RecomZ0[nelem][ion],
					     &gv.bin[nd].chrg(nz).RecomEn[nelem][ion],
					     &gv.bin[nd].chrg(nz).ChemEn[nelem][ion]);
			}
		}
	}
	return;
}

STATIC void GetPotValues(size_t nd,
			 long Zg,
			 /*@out@*/ double *ThresInf,
			 /*@out@*/ double *ThresInfVal,
			 /*@out@*/ double *ThresSurf,
			 /*@out@*/ double *ThresSurfVal,
			 /*@out@*/ double *PotSurf,
			 /*@out@*/ double *Emin,
			 bool lgUseTunnelCorr)
{
	DEBUG_ENTRY( "GetPotValues()" );

	/* >>chng 01 may 07, this routine now completely supports the hybrid grain charge model,
	 * the change for this routine is that now it is only fed integral charge states; calculation
	 * of IP has also been changed in accordance with Weingartner & Draine, 2001, PvH */

	double dZg = (double)Zg;
	/* this is average grain potential in Rydberg */
	double dstpot = chrg2pot(dZg,nd);

	/* >>chng 01 mar 20, include O(a^-2) correction terms in ionization potential */
	/* these are correction terms for the ionization potential that are
	 * important for small grains. See Weingartner & Draine, 2001, Eq. 2 */
	double IP_v = gv.bin[nd].DustWorkFcn + dstpot - 0.5*one_elec(nd) +
		(dZg+2.)*AC0/gv.bin[nd].AvRadius*one_elec(nd);

	/* >>chng 01 mar 01, use new expresssion for ThresInfVal, ThresSurfVal following the discussion
	 * with Joe Weingartner. Also include the Schottky effect (see 
	 * >>refer	grain	physics	Spitzer, 1948, ApJ, 107, 6,
	 * >>refer	grain	physics	Draine & Sutin, 1987, ApJ, 320, 803), PvH */
	if( Zg <= -1 )
	{
		pot_type pcase = gv.which_pot[gv.bin[nd].matType];
		double IP;

		IP = gv.bin[nd].DustWorkFcn - gv.bin[nd].BandGap + dstpot - 0.5*one_elec(nd);
		switch( pcase )
		{
		case POT_CAR:
			IP -= AC1G/(gv.bin[nd].AvRadius+AC2G)*one_elec(nd);
			break;
		case POT_SIL:
			/* do nothing */
			break;
		default:
			fprintf( ioQQQ, " GetPotValues detected unknown type for ionization pot: %d\n" , pcase );
			cdEXIT(EXIT_FAILURE);
		}

		/* prevent valence electron from becoming less bound than attached electron; this
		 * can happen for very negative, large graphitic grains and is not physical, PvH */
		IP_v = MAX2(IP,IP_v);

		if( Zg < -1 )
		{
			/* >>chng 01 apr 20, use improved expression for tunneling effect, PvH */
			double help = fabs(dZg+1);
			/* this is the barrier height solely due to the Schottky effect */
			*Emin = -ThetaNu(help)*one_elec(nd);
			if( lgUseTunnelCorr )
			{
				/* this is the barrier height corrected for tunneling effects */
				*Emin *= 1. - 2.124e-4/(pow(gv.bin[nd].AvRadius,(realnum)0.45)*pow(help,0.26));
			}
		}
		else
		{
			*Emin = 0.;
		}

		*ThresInf = IP - *Emin;
		*ThresInfVal = IP_v - *Emin;
		*ThresSurf = *ThresInf;
		*ThresSurfVal = *ThresInfVal;
		*PotSurf = *Emin;
	}
	else
	{
		*ThresInf = IP_v;
		*ThresInfVal = IP_v;
		*ThresSurf = *ThresInf - dstpot;
		*ThresSurfVal = *ThresInfVal - dstpot;
		*PotSurf = dstpot;
		*Emin = 0.;
	}
	return;
}


/* given grain nd in charge state nz, and incoming ion (nelem,ion),
 * detemine outgoing ion (nelem,Z0) and chemical energy ChEn released
 * ChemEn is net contribution of ion recombination to grain heating */
STATIC void GrainIonColl(size_t nd,
			 long int nz,
			 long int nelem,
			 long int ion,
			 const double phi_s_up[], /* phi_s_up[LIMELM+1] */
			 const double phi_s_dn[], /* phi_s_dn[2] */
			 /*@out@*/long *Z0,
			 /*@out@*/realnum *ChEn,
			 /*@out@*/realnum *ChemEn)
{
	DEBUG_ENTRY( "GrainIonColl()" );

	long save = ion;

	if( ion > 0 && rfield.anu(Heavy.ipHeavy[nelem][ion-1]-1) > (realnum)phi_s_up[0] )
	{
		/* ion will get electron(s) */
		*ChEn = 0.f;
		*ChemEn = 0.f;
		double phi_s = phi_s_up[0];
		do 
		{
			*ChEn += rfield.anu(Heavy.ipHeavy[nelem][ion-1]-1) - (realnum)phi_s;
			*ChemEn += rfield.anu(Heavy.ipHeavy[nelem][ion-1]-1);
			/* this is a correction for the imperfections in the n-charge state model:
			 * since the transfer gets modeled as n single-electron transfers, instead of one
			 * n-electron transfer, a correction for the difference in binding energy is needed */
			*ChemEn -= (realnum)(phi_s - phi_s_up[0]);
			--ion;
			phi_s = phi_s_up[save-ion];
		} while( ion > 0 && rfield.anu(Heavy.ipHeavy[nelem][ion-1]-1) > (realnum)phi_s );

		*Z0 = ion;
	}
	else if( ion <= nelem && gv.bin[nd].chrg(nz).DustZ > gv.bin[nd].LowestZg &&
		 rfield.anu(Heavy.ipHeavy[nelem][ion]-1) < (realnum)phi_s_dn[0] )
	{
		/* grain will get electron(s) */
		*ChEn = 0.f;
		*ChemEn = 0.f;
		long Zg = gv.bin[nd].chrg(nz).DustZ;
		double phi_s = phi_s_dn[0];
		do 
		{
			*ChEn += (realnum)phi_s - rfield.anu(Heavy.ipHeavy[nelem][ion]-1);
			*ChemEn -= rfield.anu(Heavy.ipHeavy[nelem][ion]-1);
			/* this is a correction for the imperfections in the n-charge state model:
			 * since the transfer gets modeled as n single-electron transfers, instead of one
			 * n-electron transfer, a correction for the difference in binding energy is needed */
			*ChemEn += (realnum)(phi_s - phi_s_dn[0]);
			++ion;
			--Zg;

			double d[5];
			if( ion-save < 2 )
				phi_s = phi_s_dn[ion-save];
			else
				GetPotValues(nd,Zg-1,&d[0],&d[1],&phi_s,&d[2],&d[3],&d[4],NO_TUNNEL);

		} while( ion <= nelem && Zg > gv.bin[nd].LowestZg &&
			 rfield.anu(Heavy.ipHeavy[nelem][ion]-1) < (realnum)phi_s );
		*Z0 = ion;
	}
	else
	{
		/* nothing happens */
		*ChEn = 0.f;
		*ChemEn = 0.f;
		*Z0 = ion;
	}
/*  	printf(" GrainIonColl: nelem %ld ion %ld -> %ld, ChEn %.6f\n",nelem,save,*Z0,*ChEn); */
	return;
}


/* initialize grain-ion charge transfer rates on grain species nd */
STATIC void GrainChrgTransferRates(long nd)
{
	DEBUG_ENTRY( "GrainChrgTransferRates()" );

	double fac0 = STICK_ION*gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3;

#	ifndef IGNORE_GRAIN_ION_COLLISIONS

	for( long nz=0; nz < gv.bin[nd].nChrg; nz++ )
	{
		ChargeBin& gptr = gv.bin[nd].chrg(nz);
		double fac1 = gptr.FracPop*fac0;

		if( fac1 == 0. )
			continue;

		for( long ion=0; ion <= LIMELM; ion++ )
		{
			/* >>chng 00 jul 19, replace classical results with results including image potential
			 * to correct for polarization of the grain as charged particle approaches. */
			double eta, xi;
			GrainScreen(ion,nd,nz,&eta,&xi);

			double fac2 = eta*fac1;

			if( fac2 == 0. )
				continue;

			for( long nelem=MAX2(0,ion-1); nelem < LIMELM; nelem++ )
			{
				if( dense.lgElmtOn[nelem] && ion != gptr.RecomZ0[nelem][ion] )
				{
					gv.GrainChTrRate[nelem][ion][gptr.RecomZ0[nelem][ion]] +=
						(realnum)(fac2*GetAveVelocity( dense.AtomicWeight[nelem] )*atmdat.lgCTOn);
				}
			}
		}
	}
#	endif
	return;
}


/* this routine updates all grain quantities that depend on radius,
 * except gv.dstab and gv.dstsc which are updated in GrainUpdateRadius2() */
STATIC void GrainUpdateRadius1()
{
	DEBUG_ENTRY( "GrainUpdateRadius1()" );

	for( long nelem=0; nelem < LIMELM; nelem++ )
	{
		gv.elmSumAbund[nelem] = 0.f;
	}

	/* grain abundance may be a function of depth */
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		gv.bin[nd].GrnDpth = (realnum)GrnStdDpth(nd);
		gv.bin[nd].dstAbund = (realnum)(gv.bin[nd].dstfactor*gv.GrainMetal*gv.bin[nd].GrnDpth);
		ASSERT( gv.bin[nd].dstAbund > 0.f );

		/* grain unit conversion, <unit>/H (default depl) -> <unit>/cm^3 (actual depl) */
		gv.bin[nd].cnv_H_pCM3 = dense.gas_phase[ipHYDROGEN]*gv.bin[nd].dstAbund;
		gv.bin[nd].cnv_CM3_pH = 1./gv.bin[nd].cnv_H_pCM3;
		/* grain unit conversion, <unit>/cm^3 (actual depl) -> <unit>/grain */
		gv.bin[nd].cnv_CM3_pGR = gv.bin[nd].cnv_H_pGR/gv.bin[nd].cnv_H_pCM3;
		gv.bin[nd].cnv_GR_pCM3 = 1./gv.bin[nd].cnv_CM3_pGR;

		/* >>chng 01 dec 05, calculate the number density of each element locked in grains,
		 * summed over all grain bins. this number uses the actual depletion of the grains
		 * and is already multiplied with hden, units cm^-3. */
		for( long nelem=0; nelem < LIMELM; nelem++ )
		{
			gv.elmSumAbund[nelem] += gv.bin[nd].elmAbund[nelem]*(realnum)gv.bin[nd].cnv_H_pCM3;
		}
	}
	return;
}


/* this routine adds all the grain opacities in gv.dstab and gv.dstsc, this could not be
 * done in GrainUpdateRadius1 since charge and FracPop must be converged first */
STATIC void GrainUpdateRadius2()
{
	DEBUG_ENTRY( "GrainUpdateRadius2()" );

	memset( get_ptr(gv.dstab), 0, size_t(rfield.nflux_with_check*sizeof(gv.dstab[0])) );
	memset( get_ptr(gv.dstsc), 0, size_t(rfield.nflux_with_check*sizeof(gv.dstsc[0])) );

	/* >>chng 06 oct 05 rjrw, reorder loops */
	/* >>chng 11 dec 12 reorder loops so they can be vectorized, PvH */
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		double dstAbund = gv.bin[nd].dstAbund;

		/* >>chng 01 mar 26, from nupper to nflux */
		for( long i=0; i < rfield.nflux; i++ )
		{
			/* these are total absorption and scattering cross sections,
			 * the latter should contain the asymmetry factor (1-g) */
			/* this is effective area per proton, scaled by depletion
			 * dareff(nd) = darea(nd) * dstAbund(nd) */
			/* grain abundance may be a function of depth */
			/* >>chng 02 dec 30, separated scattering cross section and asymmetry factor, PvH */
			gv.dstab[i] += gv.bin[nd].dstab1[i]*dstAbund;
		}
		for( long i=0; i < rfield.nflux; i++ )
		{
			gv.dstsc[i] += gv.bin[nd].pure_sc1[i]*gv.bin[nd].asym[i]*dstAbund;
		}

		for( long nz=0; nz < gv.bin[nd].nChrg; nz++ )
		{
			ChargeBin& gptr = gv.bin[nd].chrg(nz);
			if( gptr.DustZ <= -1 )
			{
				double FracAbund = gptr.FracPop*dstAbund;

				for( long i=gptr.ipThresInf; i < rfield.nflux; i++ )
					gv.dstab[i] += FracAbund*gptr.cs_pdt[i];
			}
		}
	}

	for( long i=0; i < rfield.nflux; i++ )
	{
		/* this must be positive, zero in case of uncontrolled underflow */
		ASSERT( gv.dstab[i] > 0. && gv.dstsc[i] > 0. );
	}
	return;
}


/* GrainTemperature computes grains temperature, and gas cooling */
STATIC void GrainTemperature(size_t nd,
			     /*@out@*/ realnum *dccool,
			     /*@out@*/ double *hcon,
			     /*@out@*/ double *hots,
			     /*@out@*/ double *hla)
{
	DEBUG_ENTRY( "GrainTemperature()" );

	/* sanity checks */
	ASSERT( nd < gv.bin.size() );

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "    GrainTemperature starts for grain %s\n", gv.bin[nd].chDstLab );
	}

	/* >>chng 01 may 07, this routine now completely supports the hybrid grain
	 * charge model, and the average charge state is not used anywhere anymore, PvH */

	/* direct heating by incident continuum (all energies) */
	*hcon = 0.;
	/* heating by diffuse ots fields */
	*hots = 0.;
	/* heating by Ly alpha alone, for output only, is already included in hots */
	*hla = 0.;

	long ipLya = iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).ipCont() - 1;

	/* integrate over ionizing continuum; energy goes to dust and gas
	 * GasHeatPhotoEl is what heats the gas */
	gv.bin[nd].GasHeatPhotoEl = 0.;

	gv.bin[nd].GrainCoolTherm = 0.;
	gv.bin[nd].thermionic = 0.;

	realnum dcheat = 0.f;
	*dccool = 0.f;

	gv.bin[nd].BolFlux = 0.;

	/* >>chng 04 jan 25, moved initialization of phiTilde to qheat_init(), PvH */

	for( long nz=0; nz < gv.bin[nd].nChrg; nz++ )
	{
		ChargeBin& gptr = gv.bin[nd].chrg(nz);

		/* >>chng 04 may 31, introduced lgReEvaluate2 to save time when iterating Tdust, PvH */
		bool lgReEvaluate1 = gptr.hcon1 < 0.;
		bool lgReEvaluate2 = gptr.hots1 < 0.;

		long ip0 = 0;
		long ip1 = min(gptr.ipThresInf,rfield.nPositive);
		long ip2 = rfield.nPositive;
		double hcon1, hots1, pe1, bolflux1, hla1;
		if( lgReEvaluate1 )
		{
			hcon1 = reduce_ab( get_ptr(gv.bin[nd].dstab1_x_anu), get_ptr(rfield.flux[0]), ip0, ip1 ) +
				reduce_ab( gptr.fac1.ptr0(), get_ptr(rfield.flux[0]), ip1, ip2 );
			gptr.hcon1 = hcon1;
		}
		else
		{
			hcon1 = gptr.hcon1;
		}
		if( lgReEvaluate2 )
		{
			hots1 = reduce_ab( get_ptr(gv.bin[nd].dstab1_x_anu), get_ptr(rfield.SummedDif), ip0, ip1 ) +
				reduce_ab( gptr.fac1.ptr0(), get_ptr(rfield.SummedDif), ip1, ip2 );
#			ifdef WD_TEST2
			pe1 = reduce_ab( gptr.fac2.ptr0(), get_ptr(rfield.flux[0]), ip1, ip2 );
#			else
			pe1 = reduce_ab( gptr.fac2.ptr0(), get_ptr(rfield.SummedCon), ip1, ip2 );
#			endif
#			ifndef NDEBUG
			bolflux1 = reduce_ab( get_ptr(gv.bin[nd].dstab1_x_anu), get_ptr(rfield.SummedCon), ip0, ip2 );
			if( gptr.DustZ <= -1 )
				bolflux1 +=
					reduce_abc( gptr.cs_pdt.ptr0(), rfield.anuptr(), get_ptr(rfield.SummedCon), ip1, ip2 );
#			else
			bolflux1 = 0.;
#			endif
			gptr.hots1 = hots1;
			gptr.pe1 = pe1;
			gptr.bolflux1 = bolflux1;
		}
		else
		{
			hots1 = gptr.hots1;
			pe1 = gptr.pe1;
			bolflux1 = gptr.bolflux1;
		}

		/*  heating by Ly A on dust in this zone,
		 *  only used for printout; Ly-a is already in OTS fields */
		/* >>chng 00 apr 18, moved calculation of hla, by PvH */
		/* >>chng 04 feb 01, moved calculation of hla1 outside loop for optimization, PvH */
		if( ipLya < MIN2(gptr.ipThresInf,rfield.nPositive) )
		{
			hla1 = rfield.otslin[ipLya]*gv.bin[nd].dstab1_x_anu[ipLya];
		}
		else if( ipLya < rfield.nPositive )
		{
			/* >>chng 00 apr 18, include photo-electric effect, by PvH */
			hla1 = rfield.otslin[ipLya]*gptr.fac1[ipLya];
		}
		else
		{
			hla1 = 0.;
		}

		ASSERT( hcon1 >= 0. && hots1 >= 0. && hla1 >= 0. && bolflux1 >= 0. && pe1 >= 0. );

		*hcon += gptr.FracPop*hcon1;
		*hots += gptr.FracPop*hots1;
		*hla += gptr.FracPop*hla1;
		gv.bin[nd].BolFlux += gptr.FracPop*bolflux1;
		if( gv.lgDHetOn )
			gv.bin[nd].GasHeatPhotoEl += gptr.FracPop*pe1;

#		ifndef NDEBUG
		if( trace.lgTrace && trace.lgDustBug )
		{
			fprintf( ioQQQ, "    Zg %ld bolflux: %.4e\n", gptr.DustZ,
			  gptr.FracPop*bolflux1*EN1RYD*gv.bin[nd].cnv_H_pCM3 );
		}
#		endif

		/* add in thermionic emissions (thermal evaporation of electrons), it gives a cooling
		 * term for the grain. thermionic emissions will not be treated separately in quantum
		 * heating since they are only important when grains are heated to near-sublimation 
		 * temperatures; under those conditions quantum heating effects will never be important.
		 * in order to maintain energy balance they will be added to the ion contribution though */
		/* ThermRate is normalized per cm^2 of grain surface area, scales with total grain area */
		double rate = gptr.FracPop*gptr.ThermRate*gv.bin[nd].IntArea*gv.bin[nd].cnv_H_pCM3;
		/* >>chng 01 mar 02, PotSurf[nz] term was incorrectly taken into account, PvH */
		double EhatThermionic = 2.*BOLTZMANN*gv.bin[nd].tedust + MAX2(gptr.PotSurf*EN1RYD,0.);
		gv.bin[nd].GrainCoolTherm += rate * (EhatThermionic + gptr.ThresSurf*EN1RYD);
		gv.bin[nd].thermionic += rate * (EhatThermionic - gptr.PotSurf*EN1RYD);
	}

	/* norm is used to convert all heating rates to erg/cm^3/s */
	double norm = EN1RYD*gv.bin[nd].cnv_H_pCM3;

	/* hcon is radiative heating by incident radiation field */
	*hcon *= norm;

	/* hots is total heating of the grain by diffuse fields */
	*hots *= norm;

	/* heating by Ly alpha alone, for output only, is already included in hots */
	*hla *= norm;

	gv.bin[nd].BolFlux *= norm;

	/* heating by thermal collisions with gas does work
	 * DCHEAT is grain collisional heating by gas
	 * DCCOOL is gas cooling due to collisions with grains
	 * they are different since grain surface recombinations
	 * heat the grains, but do not cool the gas ! */
	/* >>chng 03 nov 06, moved call after renorm of BolFlux, so that GrainCollHeating can look at it, PvH */
	GrainCollHeating(nd,&dcheat,dccool);

	/* GasHeatPhotoEl is what heats the gas */
	gv.bin[nd].GasHeatPhotoEl *= norm;

	if( gv.lgBakesPAH_heat )
	{
		/* this is a dirty hack to get BT94 PE heating rate
		 * for PAH's included, for Lorentz Center 2004 PDR meeting, PvH */
		/*>>>refer	PAH	heating	Bakes, E.L.O., & Tielens, A.G.G.M. 1994, ApJ, 427, 822 */
		/* >>chng 05 aug 12, change from +=, which added additional heating to what exists already,
		 * to simply = to set the heat, this equation gives total heating */
		gv.bin[nd].GasHeatPhotoEl = 1.e-24*hmi.UV_Cont_rel2_Habing_TH85_depth*
			dense.gas_phase[ipHYDROGEN]*(4.87e-2/(1.0+4e-3*pow((hmi.UV_Cont_rel2_Habing_TH85_depth*
			/*>>chng 06 jul 21, use phycon.sqrte in next two lines */
			phycon.sqrte/dense.eden),0.73)) + 3.65e-2*pow(phycon.te/1.e4,0.7)/
			(1.+2.e-4*(hmi.UV_Cont_rel2_Habing_TH85_depth*phycon.sqrte/dense.eden)))/gv.bin.size();

	}

	/* >>chng 06 jun 01, add optional scale factor, set with command
	 * set grains heat, to rescale PE heating as per Allers et al. 2005 */
	gv.bin[nd].GasHeatPhotoEl *= gv.GrainHeatScaleFactor;

	/*  find power absorbed by dust and resulting temperature
	 *
	 * hcon is heating from incident continuum (all energies)
	 * hots is heating from ots continua and lines
	 * dcheat is net grain collisional and chemical heating by
	 *    particle collisions and recombinations
	 * GrainCoolTherm is grain cooling by thermionic emissions
	 *
	 * GrainHeat is net heating of this grain type,
	 *    to be balanced by radiative cooling */
	gv.bin[nd].GrainHeat = *hcon + *hots + dcheat - gv.bin[nd].GrainCoolTherm;

	/* remember collisional heating for this grain species */
	gv.bin[nd].GrainHeatColl = dcheat;

	/* >>chng 04 may 31, replace ASSERT of GrainHeat > 0. with if-statement and let
	 * GrainChargeTemp sort out the consquences of GrainHeat becoming negative, PvH */
	/* in case where the thermionic rates become very large,
	 * or collisional cooling dominates, this may become negative */
	if( gv.bin[nd].GrainHeat > 0. )
	{
		bool lgOutOfBounds;
		/*  now find temperature, GrainHeat is sum of total heating of grain
		 *  >>chng 97 jul 17, divide by abundance here */
		double y, x = log(MAX2(DBL_MIN,gv.bin[nd].GrainHeat*gv.bin[nd].cnv_CM3_pH));
		/* >>chng 96 apr 27, as per Peter van Hoof comment */
		splint_safe(gv.bin[nd].dstems,gv.dsttmp,gv.bin[nd].dstslp,NDEMS,x,&y,&lgOutOfBounds);
		gv.bin[nd].tedust = (realnum)exp(y);
	}
	else
	{
		gv.bin[nd].GrainHeat = -1.;
		gv.bin[nd].tedust = -1.;
	}

	if( thermal.ConstGrainTemp > 0. )
	{
		bool lgOutOfBounds;
		/* use temperature set with constant grain temperature command */
		gv.bin[nd].tedust = thermal.ConstGrainTemp;
		/* >>chng 04 jun 01, make sure GrainHeat is consistent with value of tedust, PvH */
		double y, x = log(gv.bin[nd].tedust);
		splint_safe(gv.dsttmp,gv.bin[nd].dstems,gv.bin[nd].dstslp2,NDEMS,x,&y,&lgOutOfBounds);
		gv.bin[nd].GrainHeat = exp(y)*gv.bin[nd].cnv_H_pCM3;
	}

	/*  save for later possible printout */
	gv.bin[nd].TeGrainMax = (realnum)MAX2(gv.bin[nd].TeGrainMax,gv.bin[nd].tedust);

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "  >GrainTemperature finds %s Tdst %.5e hcon %.4e ",
			 gv.bin[nd].chDstLab, gv.bin[nd].tedust, *hcon);
		fprintf( ioQQQ, "hots %.4e dcheat %.4e GrainCoolTherm %.4e\n", 
			 *hots, dcheat, gv.bin[nd].GrainCoolTherm );
	}
	return;
}


/* helper routine for initializing quantities related to the photo-electric effect */
STATIC void PE_init(size_t nd,
		    long nz,
		    long i,
		    /*@out@*/ double *cs1,
		    /*@out@*/ double *cs2,
		    /*@out@*/ double *cs_tot,
		    /*@out@*/ double *cool1,
		    /*@out@*/ double *cool2,
		    /*@out@*/ double *ehat1,
		    /*@out@*/ double *ehat2)
{
	DEBUG_ENTRY( "PE_init()" );

	ChargeBin& gptr = gv.bin[nd].chrg(nz);
	long ipLo1 = gptr.ipThresInfVal;
	long ipLo2 = gptr.ipThresInf;

	/* sanity checks */
	ASSERT( nd < gv.bin.size() );
	ASSERT( nz >= 0 && nz < gv.bin[nd].nChrg );
	ASSERT( i >= 0 && i < rfield.nPositive );

	/** \todo xray - add fluoresence in energy balance */

	/* contribution from valence band */
	if( i >= ipLo1 )
	{
		/* effective cross section for photo-ejection */
		*cs1 = gv.bin[nd].dstab1[i]*gptr.yhat[i];
		/* >>chng 00 jul 17, use description of Weingartner & Draine, 2001 */
		/* ehat1 is the average energy of the escaping electron at infinity */
		*ehat1 = gptr.ehat[i];

		/* >>chng 01 nov 27, changed de-excitation energy to conserve energy,
		 * this branch treats valence band ionizations, but for negative grains an
		 * electron from the conduction band will de-excite into the hole in the
		 * valence band, reducing the amount of potential energy lost. It is assumed
		 * that no photons are ejected in the process. PvH */
		/* >>chng 06 mar 19, reorganized this routine in the wake of the introduction
		 * of the WDB06 X-ray physics. The basic functionality is still the same, but
		 * the meaning is not. A single X-ray photon can eject multiple electrons from
		 * either the conduction band, valence band or an inner shell. In the WDB06
		 * approximation all these electrons are assumed to be ejected from a grain
		 * with the same charge. After the primary ejection, Auger cascades will fill
		 * up any inner shell holes, so energetically it is as if all these electrons
		 * come from the outermost band (either conduction or valence band, depending
		 * on the grain charge). Recombination will also be into the outermost band
		 * so that way energy conservation is assured. It is assumed that these Auger
		 * cascades are so fast that they can be treated as a single event as far as
		 * quantum heating is concerned. */

		/* cool1 is the amount by which photo-ejection cools the grain */
		if( gptr.DustZ <= -1 )
			*cool1 = gptr.ThresSurf + gptr.PotSurf + *ehat1;
		else
			*cool1 = gptr.ThresSurfVal + gptr.PotSurf + *ehat1;

		ASSERT( *ehat1 >= 0. && *cool1 >= gptr.PotSurf );
	}
	else
	{
		*cs1 = 0.;
		*ehat1 = 0.;
		*cool1 = 0.;
	}

	/* contribution from conduction band */
	if( gptr.DustZ <= -1 && i >= ipLo2 )
	{
		/* effective cross section for photo-detechment */
		*cs2 = gptr.cs_pdt[i];
		/* ehat2 is the average energy of the escaping electron at infinity */
		*ehat2 = rfield.anu(i) - gptr.ThresSurf - gptr.PotSurf;
		/* cool2 is the amount by which photo-detechment cools the grain */
		*cool2 = rfield.anu(i);

		ASSERT( *ehat2 >= 0. && *cool2 > 0. );
	}
	else
	{
		*cs2 = 0.;
		*ehat2 = 0.;
		*cool2 = 0.;
	}

	*cs_tot = gv.bin[nd].dstab1[i] + *cs2;		
	return;
}


/* GrainCollHeating compute grains collisional heating cooling */
STATIC void GrainCollHeating(size_t nd,
			     /*@out@*/ realnum *dcheat,
			     /*@out@*/ realnum *dccool)
{
	long int ion,
	  nelem,
	  nz;
	H2_type ipH2;
	double Accommodation,
	  CollisionRateElectr,      /* rate electrons strike grains */
	  CollisionRateMol,         /* rate molecules strike grains */
	  CollisionRateIon,         /* rate ions strike grains */
	  CoolTot,
	  CoolBounce,
	  CoolEmitted,
	  CoolElectrons,
	  CoolMolecules,
	  CoolPotential,
	  CoolPotentialGas,
	  eta,
	  HeatTot,
	  HeatBounce,
	  HeatCollisions,
	  HeatElectrons,
	  HeatIons,
	  HeatMolecules,
	  HeatRecombination, /* sum of abundances of ions times velocity times ionization potential times eta */
	  HeatChem,
	  HeatCor,
	  Stick,
	  ve,
	  WeightMol,
	  xi;

	/* energy deposited into grain by formation of a single H2 molecule, in eV,
	 * >>refer	grain	physics	Takahashi J., Uehara H., 2001, ApJ, 561, 843 */
	const double H2_FORMATION_GRAIN_HEATING[H2_TOP] = { 0.20, 0.4, 1.72 };

	DEBUG_ENTRY( "GrainCollHeating()" );


	/* >>chng 01 may 07, this routine now completely supports the hybrid grain
	 * charge model, and the average charge state is not used anywhere anymore, PvH */

	/* this subroutine evaluates the gas heating-cooling rate
	 * (erg cm^-3 s^-1) due to grain gas collisions.
	 * the net effect can be positive or negative,
	 * depending on whether the grains or gas are hotter
	 * the physics is described in 
	 * >>refer	grain	physics	Baldwin, Ferland, Martin et al., 1991, ApJ 374, 580 */

	HeatTot = 0.;
	CoolTot = 0.;

	HeatIons = 0.;

	gv.bin[nd].ChemEn = 0.;

	/* loop over the charge states */
	for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
	{
		ChargeBin& gptr = gv.bin[nd].chrg(nz);

		/* HEAT1 will be rate collisions heat the grain
		 * COOL1 will be rate collisions cool the gas kinetics */
		double Heat1 = 0.;
		double Cool1 = 0.;
		double ChemEn1 = 0.;

		/* ============================================================================= */
		/* heating/cooling due to neutrals and positive ions */

		/* loop over all stages of ionization */
		for( ion=0; ion <= LIMELM; ion++ )
		{
			/* this is heating of grains due to recombination energy of species,
			 * and assumes that every ion is fully neutralized upon striking the grain surface.
			 * all radiation produced in the recombination process is absorbed within the grain
			 *
			 * ion=0 are neutrals, ion=1 are single ions, etc
			 * each population is weighted by the AVERAGE velocity
			 * */
			CollisionRateIon = 0.;
			CoolPotential = 0.;
			CoolPotentialGas = 0.;
			HeatRecombination = 0.;
			HeatChem = 0.;

			/* >>chng 00 jul 19, replace classical results with results including image potential
			 * to correct for polarization of the grain as charged particle approaches. */
			GrainScreen(ion,nd,nz,&eta,&xi);

			for( nelem=MAX2(0,ion-1); nelem < LIMELM; nelem++ )
			{
				if( dense.lgElmtOn[nelem] && dense.xIonDense[nelem][ion] > 0. )
				{
					double CollisionRateOne;

					/* >>chng 00 apr 05, use correct accomodation coefficient, by PvH
					 * the coefficient is defined at the end of appendix A.10 of BFM
					 * assume ion sticking prob is unity */
#if defined( IGNORE_GRAIN_ION_COLLISIONS )
					Stick = 0.;
#elif defined( WD_TEST2 )
					Stick = ( ion == gptr.RecomZ0[nelem][ion] ) ?
						0. : STICK_ION;
#else
					Stick = ( ion == gptr.RecomZ0[nelem][ion] ) ?
						gv.bin[nd].AccomCoef[nelem] : STICK_ION;
#endif
					/* this is rate with which charged ion strikes grain */
					/* >>chng 00 may 02, this had left 2./SQRTPI off */
					/* >>chng 00 may 05, use average speed instead of 2./SQRTPI*Doppler, PvH */
					CollisionRateOne = Stick*dense.xIonDense[nelem][ion]*
						GetAveVelocity( dense.AtomicWeight[nelem] );
					CollisionRateIon += CollisionRateOne;
					/* >>chng 01 nov 26, use PotSurfInc when appropriate:
					 * the values for the surface potential used here make it
					 * consistent with the rest of the code and preserve energy.
					 * NOTE: For incoming particles one should use PotSurfInc with
					 * Schottky effect for positive ion, for outgoing particles
					 * one should use PotSurf for Zg+ion-Z_0-1 (-1 because PotSurf
					 * assumes electron going out), these corrections are small
					 * and will be neglected for now, PvH */
					if( ion >= gptr.RecomZ0[nelem][ion] )
					{
						CoolPotential += CollisionRateOne * (double)ion *
							gptr.PotSurf;
						CoolPotentialGas += CollisionRateOne *
							(double)gptr.RecomZ0[nelem][ion] *
							gptr.PotSurf;
					}
					else
					{
						CoolPotential += CollisionRateOne * (double)ion *
							gptr.PotSurfInc;
						CoolPotentialGas += CollisionRateOne *
							(double)gptr.RecomZ0[nelem][ion] *
							gptr.PotSurfInc;
					}
					/* this is sum of all energy liberated as ion recombines to Z0 in grain */
					/* >>chng 00 jul 05, subtract energy needed to get 
					 * electron out of grain potential well, PvH */
					/* >>chng 01 may 09, chemical energy now calculated in GrainIonColl, PvH */
					HeatRecombination += CollisionRateOne *
						gptr.RecomEn[nelem][ion];
					HeatChem += CollisionRateOne * gptr.ChemEn[nelem][ion];
				}
			}

			/* >>chng 00 may 01, Boltzmann factor had multiplied all of factor instead
			 * of only first and last term.  pvh */

			/* equation 29 from Balwin et al 91 */
			/* this is direct collision rate, 2kT * xi, first term in eq 29 */
			HeatCollisions = CollisionRateIon * 2.*BOLTZMANN*phycon.te*xi;
			/* this is change in energy due to charge acceleration within grain's potential 
			 * this is exactly balanced by deceleration of incoming electrons and accelaration
			 * of outgoing photo-electrons and thermionic emissions; all these terms should
			 * add up to zero (total charge of grain should remain constant) */
			CoolPotential *= eta*EN1RYD;
			CoolPotentialGas *= eta*EN1RYD;
			/* this is recombination energy released within grain */
			HeatRecombination *= eta*EN1RYD;
			HeatChem *= eta*EN1RYD;
			/* energy carried away by neutrals after recombination, so a cooling term */
			CoolEmitted = CollisionRateIon * 2.*BOLTZMANN*gv.bin[nd].tedust*eta;

			/* total GraC 0 in the emission line output */
			Heat1 += HeatCollisions - CoolPotential + HeatRecombination - CoolEmitted;

			/* rate kinetic energy lost from gas - gas cooling - eq 32 in BFM */
			/* this GrGC 0 in the main output */
			/* >>chng 00 may 05, reversed sign of gas cooling contribution */
			Cool1 += HeatCollisions - CoolEmitted - CoolPotentialGas;

			ChemEn1 += HeatChem;
		}

		/* remember grain heating by ion collisions for quantum heating treatment */
		HeatIons += gptr.FracPop*Heat1;

		if( trace.lgTrace && trace.lgDustBug )
		{
			fprintf( ioQQQ, "    Zg %ld ions heat/cool: %.4e %.4e\n", gptr.DustZ,
			  gptr.FracPop*Heat1*gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3,
			  gptr.FracPop*Cool1*gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3 );
		}

		/* ============================================================================= */
		/* heating/cooling due to electrons */

		ion = -1;
		Stick = ( gptr.DustZ <= -1 ) ? gv.bin[nd].StickElecNeg : gv.bin[nd].StickElecPos;
		/* VE is mean (not RMS) electron velocity */
		/*ve = TePowers.sqrte*6.2124e5;*/
		ve = sqrt(8.*BOLTZMANN/PI/ELECTRON_MASS*phycon.te);

		/* electron arrival rate - eqn 29 again */
		CollisionRateElectr = Stick*dense.eden*ve;

		/* >>chng 00 jul 19, replace classical results with results including image potential
		 * to correct for polarization of the grain as charged particle approaches. */
		GrainScreen(ion,nd,nz,&eta,&xi);

		if( gptr.DustZ > gv.bin[nd].LowestZg )
		{
			HeatCollisions = CollisionRateElectr*2.*BOLTZMANN*phycon.te*xi;
			/* this is change in energy due to charge acceleration within grain's potential 
			 * this term (perhaps) adds up to zero when summed over all charged particles */
			CoolPotential = CollisionRateElectr * (double)ion*gptr.PotSurfInc*eta*EN1RYD;
			/* >>chng 00 jul 05, this is term for energy released due to recombination, PvH */
			HeatRecombination = CollisionRateElectr * gptr.ThresSurfInc*eta*EN1RYD;
			HeatBounce = 0.;
			CoolBounce = 0.;
		}
		else
		{
			HeatCollisions = 0.;
			CoolPotential = 0.;
			HeatRecombination = 0.;
			/* >>chng 00 jul 05, add in terms for electrons that bounce off grain, PvH */
			/* >>chng 01 mar 09, remove these terms, their contribution is negligible, and replace
			 * them with similar terms that describe electrons that are captured by grains at Z_min,
			 * these electrons are not in a bound state and the grain will quickly autoionize, PvH */
			HeatBounce = CollisionRateElectr * 2.*BOLTZMANN*phycon.te*xi;
			/* >>chng 01 mar 14, replace (2kT_g - phi_g) term with -EA; for autoionizing states EA is
			 * usually higher than phi_g, so more energy is released back into the electron gas, PvH */ 
			CoolBounce = CollisionRateElectr *
				(-gptr.ThresSurfInc-gptr.PotSurfInc)*EN1RYD*eta;
			CoolBounce = MAX2(CoolBounce,0.);
		}

		/* >>chng 00 may 02, CoolPotential had not been included */
		/* >>chng 00 jul 05, HeatRecombination had not been included */
		HeatElectrons = HeatCollisions-CoolPotential+HeatRecombination+HeatBounce-CoolBounce;
		Heat1 += HeatElectrons;

		CoolElectrons = HeatCollisions+HeatBounce-CoolBounce;
		Cool1 += CoolElectrons;

		if( trace.lgTrace && trace.lgDustBug )
		{
			fprintf( ioQQQ, "    Zg %ld electrons heat/cool: %.4e %.4e\n", gptr.DustZ,
			  gptr.FracPop*HeatElectrons*gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3,
			  gptr.FracPop*CoolElectrons*gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3 );
		}

		/* add quantum heating due to recombination of electrons, subtract thermionic cooling */

		/* calculate net heating rate in erg/H/s at standard depl
		 * include contributions for recombining electrons, autoionizing electrons
		 * and subtract thermionic emissions here since it is inverse process
		 *
		 * NB - in extreme conditions this rate may become negative (if there
		 * is an intense radiation field leading to very hot grains, but no ionizing
		 * photons, hence very few free electrons). we assume that the photon rates
		 * are high enough under those circumstances to avoid phiTilde becoming negative,
		 * but we will check that in qheat1 anyway. */
		gptr.HeatingRate2 = HeatElectrons*gv.bin[nd].IntArea/4. -
			gv.bin[nd].GrainCoolTherm*gv.bin[nd].cnv_CM3_pH;

		/* >>chng 04 jan 25, moved inclusion into phitilde to qheat_init(), PvH */

		/* heating/cooling above is in erg/s/cm^2 -> multiply with projected grain area per cm^3 */
		/* GraC 0 is integral of dcheat, the total collisional heating of the grain */
		HeatTot += gptr.FracPop*Heat1;

		/* GrGC 0 total cooling of gas integrated */
		CoolTot += gptr.FracPop*Cool1;

		gv.bin[nd].ChemEn += gptr.FracPop*ChemEn1;
	}

	/* ============================================================================= */
	/* heating/cooling due to molecules */

	/* these rates do not depend on charge, hence they are outside of nz loop */

	/* sticking prob for H2 onto grain,
	 * estimated from accomodation coefficient defined at end of A.10 in BFM */
	WeightMol = 2.*dense.AtomicWeight[ipHYDROGEN];
	Accommodation = 2.*gv.bin[nd].atomWeight*WeightMol/POW2(gv.bin[nd].atomWeight+WeightMol);
	/* molecular hydrogen onto grains */
#ifndef IGNORE_GRAIN_ION_COLLISIONS
	/*CollisionRateMol = Accommodation*findspecies("H2")->den* */
	CollisionRateMol = Accommodation*hmi.H2_total*
		sqrt(8.*BOLTZMANN/PI/ATOMIC_MASS_UNIT/WeightMol*phycon.te);
	/* >>chng 03 feb 12, added grain heating by H2 formation on the surface, PvH 
	 * >>refer	grain	H2 heat	Takahashi & Uehara, ApJ, 561, 843 */
	ipH2 = gv.which_H2distr[gv.bin[nd].matType];
	/* this is rate in erg/cm^3/s */
	/* >>chng 04 may 26, changed dense.gas_phase[ipHYDROGEN] -> dense.xIonDense[ipHYDROGEN][0], PvH */
	gv.bin[nd].ChemEnH2 = gv.bin[nd].rate_h2_form_grains_used*dense.xIonDense[ipHYDROGEN][0]*
		H2_FORMATION_GRAIN_HEATING[ipH2]*EN1EV;
	/* convert to rate per cm^2 of projected grain surface area used here */
	gv.bin[nd].ChemEnH2 /=	gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3;
#else
	CollisionRateMol = 0.;
	gv.bin[nd].ChemEnH2 = 0.;
#endif

	/* now add in CO */
	WeightMol = dense.AtomicWeight[ipCARBON] + dense.AtomicWeight[ipOXYGEN];
	Accommodation = 2.*gv.bin[nd].atomWeight*WeightMol/POW2(gv.bin[nd].atomWeight+WeightMol);
#ifndef IGNORE_GRAIN_ION_COLLISIONS
	CollisionRateMol += Accommodation*findspecieslocal("CO")->den*
		sqrt(8.*BOLTZMANN/PI/ATOMIC_MASS_UNIT/WeightMol*phycon.te);
#else
	CollisionRateMol = 0.;
#endif

	/* xi and eta are unity for neutrals and so ignored */
	HeatCollisions = CollisionRateMol * 2.*BOLTZMANN*phycon.te;
	CoolEmitted = CollisionRateMol * 2.*BOLTZMANN*gv.bin[nd].tedust;

	HeatMolecules = HeatCollisions - CoolEmitted + gv.bin[nd].ChemEnH2;
	HeatTot += HeatMolecules;

	/* >>chng 00 may 05, reversed sign of gas cooling contribution */
	CoolMolecules = HeatCollisions - CoolEmitted;
	CoolTot += CoolMolecules;

	gv.bin[nd].RateUp = 0.;
	gv.bin[nd].RateDn = 0.;
	HeatCor = 0.;
	for( nz=0; nz < gv.bin[nd].nChrg; nz++ )
	{
		double d[4];
		double rate_dn = GrainElecRecomb1(nd,nz,&d[0],&d[1]);
		double rate_up = GrainElecEmis1(nd,nz,&d[0],&d[1],&d[2],&d[3]);

		gv.bin[nd].RateUp += gv.bin[nd].chrg(nz).FracPop*rate_up;
		gv.bin[nd].RateDn += gv.bin[nd].chrg(nz).FracPop*rate_dn;

		 /** \todo	2	a self-consistent treatment for the heating by Compton recoil should be used */
		HeatCor += (gv.bin[nd].chrg(nz).FracPop*rate_up*gv.bin[nd].chrg(nz).ThresSurf -
			    gv.bin[nd].chrg(nz).FracPop*rate_dn*gv.bin[nd].chrg(nz).ThresSurfInc +
			    gv.bin[nd].chrg(nz).FracPop*rate_up*gv.bin[nd].chrg(nz).PotSurf -
			    gv.bin[nd].chrg(nz).FracPop*rate_dn*gv.bin[nd].chrg(nz).PotSurfInc)*EN1RYD;
	}
	/* >>chng 01 nov 24, correct for imperfections in the n-charge state model,
	 * these corrections should add up to zero, but are actually small but non-zero, PvH */
	HeatTot += HeatCor;

	if( trace.lgTrace && trace.lgDustBug )
	{
		fprintf( ioQQQ, "    molecules heat/cool: %.4e %.4e heatcor: %.4e\n",
			 HeatMolecules*gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3,
			 CoolMolecules*gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3,
			 HeatCor*gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3 );
	}

	*dcheat = (realnum)(HeatTot*gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3);
	*dccool = ( gv.lgDColOn ) ? (realnum)(CoolTot*gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3) : 0.f;

	gv.bin[nd].ChemEn *= gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3;
	gv.bin[nd].ChemEnH2 *= gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3;

	/* add quantum heating due to molecule/ion collisions */

	/* calculate heating rate in erg/H/s at standard depl
	 * include contributions from molecules/neutral atoms and recombining ions
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
	/* >>chng 01 nov 26, add in HeatCor as well, otherwise energy imbalance will result, PvH */
	gv.bin[nd].HeatingRate1 = (HeatMolecules+HeatIons+HeatCor)*gv.bin[nd].IntArea/4.;

	/* >>chng 04 jan 25, moved inclusion into phiTilde to qheat_init(), PvH */
	return;
}


/* GrainDrift computes grains drift velocity */
void GrainDrift()
{
	long int i, 
	  loop; 
	double alam, 
	  corr, 
	  dmomen, 
	  fac, 
	  fdrag, 
	  g0, 
	  g2, 
	  phi2lm, 
	  psi, 
	  rdust, 
	  si, 
	  vdold, 
	  volmom;

	DEBUG_ENTRY( "GrainDrift()" );

	vector<realnum> help( rfield.nPositive );
	for( i=0; i < rfield.nPositive; i++ )
	{
		help[i] = (rfield.flux[0][i]+rfield.ConInterOut[i]+rfield.outlin[0][i]+rfield.outlin_noplot[i])*
			rfield.anu(i);
	}

	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		/* find momentum absorbed by grain */
		dmomen = 0.;
		for( i=0; i < rfield.nPositive; i++ )
		{
			/* >>chng 02 dec 30, separated scattering cross section and asymmetry factor, PvH */
			dmomen += help[i]*(gv.bin[nd].dstab1[i] + gv.bin[nd].pure_sc1[i]*gv.bin[nd].asym[i]);
		}
		ASSERT( dmomen >= 0. );
		dmomen *= EN1RYD*4./gv.bin[nd].IntArea;

		/* now find force on grain, and drift velocity */
		fac = 2*BOLTZMANN*phycon.te;

		/* now PSI defined by 
		 * >>refer	grain	physics	Draine and Salpeter 79 Ap.J. 231, 77 (1979) */
		psi = gv.bin[nd].dstpot*TE1RYD/phycon.te;
		if( psi > 0. )
		{
			rdust = 1.e-6;
			alam = log(20.702/rdust/psi*phycon.sqrte/dense.SqrtEden);
		}
		else
		{
			alam = 0.;
		}

		phi2lm = POW2(psi)*alam;
		corr = 2.;
		/* >>chng 04 jan 31, increased loop limit 10 -> 50, precision -> 0.001, PvH */
		for( loop = 0; loop < 50 && fabs(corr-1.) > 0.001; loop++ )
		{
			vdold = gv.bin[nd].DustDftVel;

			/* interactions with protons */
			si = gv.bin[nd].DustDftVel/phycon.sqrte*7.755e-5;
			g0 = 1.5045*si*sqrt(1.+0.4418*si*si);
			g2 = si/(1.329 + POW3(si));

			/* drag force due to protons, both linear and square in velocity
			 * equation 4 from D+S Ap.J. 231, p77. */
			fdrag = fac*dense.xIonDense[ipHYDROGEN][1]*(g0 + phi2lm*g2);

			/* drag force due to interactions with electrons */
			si = gv.bin[nd].DustDftVel/phycon.sqrte*1.816e-6;
			g0 = 1.5045*si*sqrt(1.+0.4418*si*si);
			g2 = si/(1.329 + POW3(si));
			fdrag += fac*dense.eden*(g0 + phi2lm*g2);

			/* drag force due to collisions with hydrogen and helium atoms */
			si = gv.bin[nd].DustDftVel/phycon.sqrte*7.755e-5;
			g0 = 1.5045*si*sqrt(1.+0.4418*si*si);
			fdrag += fac*(dense.xIonDense[ipHYDROGEN][0] + 1.1*dense.xIonDense[ipHELIUM][0])*g0;

			/* drag force due to interactions with helium ions */
			si = gv.bin[nd].DustDftVel/phycon.sqrte*1.551e-4;
			g0 = 1.5045*si*sqrt(1.+0.4418*si*si);
			g2 = si/(1.329 + POW3(si));
			fdrag += fac*dense.xIonDense[ipHELIUM][1]*(g0 + phi2lm*g2);

			/* this term does not work
			 *  2      HEIII*(G0+4.*PSI**2*(ALAM-0.693)*G2) )
			 * this is total momentum absorbed by dust per unit vol */
			volmom = dmomen/SPEEDLIGHT;

			if( fdrag > 0. )
			{
				corr = sqrt(volmom/fdrag);
				gv.bin[nd].DustDftVel = (realnum)(vdold*corr);
			}
			else
			{
				corr = 1.;
				gv.lgNegGrnDrg = true;
				gv.bin[nd].DustDftVel = 0.;
			}

			if( trace.lgTrace && trace.lgDustBug )
			{
				fprintf( ioQQQ, "     %2ld new drift velocity:%10.2e momentum absorbed:%10.2e\n", 
				  loop, gv.bin[nd].DustDftVel, volmom );
			}
		}
	}
	return;
}

/* GrnVryDpth sets the grain abundance as a function of depth into cloud 
 * this is intended as a playpen where the user can alter things at will 
 * standard, documented, code should go in GrnStdDpth */
STATIC double GrnVryDpth(

/* nd is the number of the grain bin. The values are listed in the Cloudy output,
 * under "Average Grain Properties", and can easily be obtained by doing a trial
 * run without varying the grain abundance and setting stop zone to 1 */

	size_t nd)
{
	DEBUG_ENTRY( "GrnVryDpth()" );

	if( nd >= gv.bin.size() )
	{
		fprintf( ioQQQ, "invalid parameter for GrnVryDpth\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* --- DEFINE THE USER-SUPPLIED GRAIN ABUNDANCE FUNCTION HERE --- */

	/* This is the code that gets activated by the keyword "function" on the command line */

	/* NB some quantities may still be undefined on the first call to this routine. */

	/* the scale factor is the hydrogen atomic fraction, small when gas is ionized or molecular and
	 * unity when atomic. This function is observed for PAHs across the Orion Bar, the PAHs are
	 * strong near the ionization front and weak in the ionized and molecular gas */
	double GrnVryDpth_v = dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN];

	/* This routine must return a scale factor >= 1.e-10 for the grain abundance at this position.
	 * See Section A.3 in Hazy for more details on how grain abundances are calculated. The function
	 * A_rel(r) mentioned there is this function times the multiplication factor set with the METALS
	 * command (usually 1) and the multiplication factor set with the GRAINS command */
	return max(1.e-10,GrnVryDpth_v);
}
