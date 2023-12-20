/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "iso.h"

#include "two_photon.h"
#include "freebound.h"
#include "parser.h"
#include "elementnames.h"
#include "atmdat.h"

t_isoCTRL iso_ctrl;

t_iso_sp iso_sp[NISO][LIMELM];

long int max_num_levels = 0;

void t_isoCTRL::zero()
{
	DEBUG_ENTRY( "t_isoCTRL::zero()" );

	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		/* option to disable continuum lowering */
		lgContinuumLoweringEnabled[ipISO] = true;

		/* flag set by compile he-like command, says to regenerate table of recombination coef */
		lgCompileRecomb[ipISO] = false;
		lgNoRecombInterp[ipISO] = false;

		/* how the gbar cs will be treated - set with atom he-like gbar command */
		/** \todo	2	change this to CS_new */
		lgCS_Vriens[ipISO] = false;
		lgCS_Lebedev[ipISO]= true;
		lgCS_Vrinceanu[ipISO] = false;
		lgCS_Fujim[ipISO] = false;
		lgCS_vrgm[ipISO] = false;
		lgCS_PSM20[ipISO] = true;
		lgCS_PS64[ipISO] = true;
		lgCS_PSClassic[ipISO] = false;
		lgCS_VOS12[ipISO] = false; // lgCS_Vrinceanu[ipISO] == true overrides
		lgCS_VOS12QM[ipISO] = false;
		/*Seaton  M. J. 1962, Proc. Phys. Soc. 79, 1105 treatment for l<=3 is default */
		lgCS_Seaton[ipISO] = true;
		lgCS_B72[ipISO] = false;
		lgCS_PSdeg[ipISO] = true;

		fixit("make this the default for ipH_LIKE if not too slow.");
		lgCS_Vrinceanu[ipH_LIKE] = false;
		lgCS_PS64[ipH_LIKE] = true;

		lgCS_therm_ave[ipISO] = false;
		lgCS_VOS_thermal[ipISO] = false;
		lgCS_None[ipISO] = false;
		/* when set try actually set to 1 or 2, depending on which fit is to be used,
		 * 1 is the broken power law fit */
		/* >>chng 02 dec 21, change to broken power law fit */
		nCS_new[ipISO] = 1;
		/* This flag says whether the density is high enough that helium is sufficiently l-mixed. */
		lgCritDensLMix[ipISO] = true;
		/* This is the flag saying whether to generate errors.  false means don't.	*/
		lgRandErrGen[ipISO] = false;
		/* this is the flag saying whether we should include excess recombination in the
		 * helike sequence.  Should only be off if testing effect of top off approximations. */
		lgTopoff[ipISO] = true;
		/* Dielectronic recombination for helike ions is on by default.	*/
		lgDielRecom[ipISO] = true;

		/* number of Lyman lines to include in opacities, this can be vastly larger
		 * than the number of actual levels in the model atom */
		nLyman[ipISO] = 100;
		nLyman_max[ipISO] = 100;
		nLyman_alloc[ipISO] = 100;

		/* controls whether l-mixing and collisional ionization included */
		lgColl_l_mixing[ipISO] = true;
		lgColl_excite[ipISO] = true;
		lgColl_ionize[ipISO] = true;
		lgLTE_levels[ipISO] = false;
		lgPrintNumberOfLevels = false;

		for (long nelem=0; nelem<LIMELM; ++nelem)
		{
			RRC_TeUsed[ipISO][nelem]=0.;
			/* Masing of levels is allowed by default */
			lgNoMaser[ipISO][nelem] = false;
		}
	}

	/* Dielectronic recombination forming hydrogen-like ions does not exist. */
	lgDielRecom[ipH_LIKE] = false;

	/* smallest transition probability allowed */
	SmallA = 1e-30f;

	/* reset with SET IND2 command, turns on/off induced two photon */
	lgInd2nu_On = false;

	/* hydrogen redistribution functions */
	ipLyaRedist[ipH_LIKE] = ipPRD;
	ipResoRedist[ipH_LIKE] = ipCRD;
	ipSubRedist[ipH_LIKE] = ipCRDW;

	/* this is the upper level for each Lya, which uses the special ipLY_A */
	nLyaLevel[ipH_LIKE] = ipH2p;
	nLyaLevel[ipHE_LIKE] = ipHe2p1P;

	/* he-like redistribution functions */
	ipLyaRedist[ipHE_LIKE] = ipPRD;
	ipResoRedist[ipHE_LIKE] = ipCRD;
	ipSubRedist[ipHE_LIKE] = ipCRDW;

	lgPessimisticErrors = false;

	/* do not average collision strengths - evaluate at kT 
	 * set true with command SET COLLISION STRENGHTS AVERAGE */
	lgCollStrenThermAver = false;

	/* option to save fine structure components on the line stack */
	lgKeepFS = false;
}

void t_iso_sp::Reset()
{
	// this is flag indicating which type of model atom to use 
	strcpy( chTypeAtomUsed , "none" );
	CaseBCheck = 0.;
	/* a first guess at the recombination coefficients */
	RadRec_caseB = 1e-13;
	lgLevelsLowered = false;
	lgLevelsEverLowered = false;
	lgMustReeval = false;
	lgPopsRescaled = false;
	/* error generation done yet? false means not done.	*/
	lgErrGenDone = false;
	for( vector<two_photon>::iterator it = TwoNu.begin(); it != TwoNu.end(); ++it )
		(*it).Reset();
	for( vector<freeBound>::iterator it = fb.begin(); it != fb.end(); ++it )
		(*it).Reset();
}

long t_iso_sp::QN2Index(QNPack ind)
{
	auto p = QNPack2Index.find(ind);
	if( p != QNPack2Index.end() )
		return p->second;
	else
	{
		// also match an unresolved term
		QNPack weak1 = ind | 0xffffULL;
		// and also match a collaped level
		QNPack weak3 = ind | 0xffffffffffffULL;
		for( long ip=0; ip < numLevels_max; ip++ )
		{
			QNPack ind2 = QN2ind(st[ip].n(), st[ip].l(), st[ip].S(), 2*st[ip].j()+1);
			if( ind2 == ind || ind2 == weak1 || ind2 == weak3 )
			{
				QNPack2Index[ind] = ip;
				return ip;
			}
		}
		return -1;
	}
}

void iso_init()
{
	DEBUG_ENTRY( "iso_init()" );

	/* variables for H-like sequence */
	/* default number of levels for hydrogen iso sequence */
	for( int nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		iso_sp[ipH_LIKE][nelem].n_HighestResolved_max = 5;
		iso_sp[ipH_LIKE][nelem].nCollapsed_max = 2;
	}

	/* add more collapsed levels for these abundant elements.  Tests
	 * show collapsed levels have very little impact on time */
	iso_sp[ipH_LIKE][ipCARBON].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipNITROGEN].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipOXYGEN].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipNEON].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipSILICON].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipMAGNESIUM].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipSULPHUR].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][ipIRON].nCollapsed_max = 5;
	iso_sp[ipH_LIKE][LIMELM-1].nCollapsed_max = 5;

	/* H and He are special cases since very high resolution, S/N
	 * optical / IR data are common
	 */
	iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max = 10;
	iso_sp[ipH_LIKE][ipHYDROGEN].nCollapsed_max = 15;

	iso_sp[ipH_LIKE][ipHELIUM].n_HighestResolved_max = 10;
	iso_sp[ipH_LIKE][ipHELIUM].nCollapsed_max = 15;

	/* variables for He-like sequence */
	/* "he-like" hydrogen (H-) is treated elsewhere */
	iso_sp[ipHE_LIKE][ipHYDROGEN].n_HighestResolved_max = -LONG_MAX;
	iso_sp[ipHE_LIKE][ipHYDROGEN].numLevels_max = -LONG_MAX;
	iso_sp[ipHE_LIKE][ipHYDROGEN].nCollapsed_max = -LONG_MAX;

	for( int nelem=ipHELIUM; nelem < LIMELM; ++nelem )
	{
		/* put at least three resolved and 1 collapsed in every element for he-like
		 * An n shell for He has twice the number of resolved levels due to singlet/triplet */
		iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max = 3;
		iso_sp[ipHE_LIKE][nelem].nCollapsed_max = 2;
	}

	/* And n=5 for these because they are most abundant */
	iso_sp[ipHE_LIKE][ipCARBON].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipCARBON].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipNITROGEN].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipNITROGEN].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipOXYGEN].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipOXYGEN].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipNEON].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipNEON].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipSILICON].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipSILICON].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipMAGNESIUM].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipMAGNESIUM].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipSULPHUR].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipSULPHUR].nCollapsed_max = 5;
	iso_sp[ipHE_LIKE][ipIRON].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][ipIRON].nCollapsed_max = 5;
	/* also set this, for exercising any possible issues with highest charge models */
	iso_sp[ipHE_LIKE][LIMELM-1].n_HighestResolved_max = 5;
	iso_sp[ipHE_LIKE][LIMELM-1].nCollapsed_max = 5;

	/* He I is very special case - we want to do a good job, lots of levels */
	iso_sp[ipHE_LIKE][ipHELIUM].n_HighestResolved_max = 6;
	iso_sp[ipHE_LIKE][ipHELIUM].nCollapsed_max = 20;

	for( int ipISO = ipH_LIKE; ipISO < NISO; ipISO++ )
	{
		for( int nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* set this to LONG_MAX, reduce to actual number later,
			 * then verify number of levels is not increased after initial coreload */
			iso_sp[ipISO][nelem].numLevels_alloc = LONG_MAX;
			iso_update_num_levels( ipISO, nelem );
		}
	}
}
long int iso_Max_Emitting_Level(long nelem, long ipISO, bool lgPrnIsoCollapsed)
{

	/* 17 02 18, gjf, set number of levels we want to print,
	 * first was default in <= C13, only print resolved levels,
	 * With C17 we report collapsed levels
	 * since their emission is accurate if the gas is dense enough to l-mix
	 * the l levels within the n shell.  Through much of development
	 * after c13 we reported to the very top of the model.  Tests
	 * show that the highest levels have "edge" effects introduced by
	 * the top of the model.  We do not report emission from the top four levels
	 * for this reason.
	 * lgPrnIsoCollapsed is init to true and set false with
	 * "print line iso collapsed off" command. */
	long int nLoop  = iso_sp[ipISO][nelem].numLevels_max - iso_sp[ipISO][nelem].nCollapsed_max;
	/* possible that nCollapsed_max < 4 */

	if( lgPrnIsoCollapsed )
	{
		if(ipISO==ipH_LIKE)
		{
			int nTooHigh = MIN2( 4, iso_sp[ipH_LIKE][nelem].nCollapsed_max );
			nLoop  = iso_sp[ipISO][nelem].numLevels_max-nTooHigh;
		}
		else if(ipISO==ipHE_LIKE)
			nLoop  = iso_sp[ipISO][nelem].numLevels_max;
	}

	return nLoop;
}

void iso_init_energies()
{
	DEBUG_ENTRY( "iso_init_energies()" );

	string sep = cpu.i().chDirSeparator();

	DataParser d;
	for( long ipISO=ipH_LIKE; ipISO < NISO; ++ipISO )
	{
		for( long nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			string elm = elementnames.chElementSym[ nelem ];
			if( elm.substr( 1, 1 ) == " " )
				elm = elm.substr( 0, 1 );
			elm[ 0 ] = tolower( elm[ 0 ] );

			string species = elm + '_' + to_string( nelem - ipISO + 1 );
			string fnam = "stout" + sep + elm + sep + species + sep + species + ".nrg.dat"; 
			if( false ) printf( "'%s'\n", fnam.c_str() );

			d.open(fnam, ES_NONE);
			d.getline();
			d.checkMagic(ENERGIESMAGIC);
			while( d.getline() )
			{
				long n, l, s, j;
				d.getToken(n);
				d.getToken(l);
				d.getToken(s);
				d.getToken(j);
				QNPack ind = QN2ind(n, l, s, 2*j+1);
				d.getToken(iso_sp[ipISO][nelem].Energy[ind]);
				if( ! d.lgEOL() )
				{
					double error;
					d.getToken( error );
				}
				d.checkEOL();
			}
			d.checkEOD();

			double IP = atmdat.getIonPot(nelem, nelem-ipISO);
			iso_sp[ipISO][nelem].IonPot = Energy(IP, "Ryd").WN();
		}
	}
}
