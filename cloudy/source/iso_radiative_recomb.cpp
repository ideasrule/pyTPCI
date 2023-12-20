/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_radiative_recomb find state specific creation and destruction rates for iso sequences */
/*iso_RRCoef_Te - evaluate radiative recombination coef at some temperature */
/*iso_recomb_check - called by SanityCheck to confirm that recombination coef are ok,
 * return value is relative error between new calculation of recom, and interp value */

#include "cddefines.h"
#include "atmdat_adfa.h"
#include "conv.h"
#include "cosmology.h"
#include "elementnames.h"
#include "helike_recom.h"
#include "hydrogenic.h"
#include "ionbal.h"
#include "iso.h"
#include "opacity.h"
#include "phycon.h"
#include "prt.h"
#include "save.h"
#include "thirdparty.h"
#include "trace.h"
#include "rt.h"
#include "freebound.h"
#include "dense.h"
#include "integrate.h"
#include "container_classes.h"

/* this will save log of radiative recombination rate coefficients at N_ISO_TE_RECOMB temperatures.
 * there will be NumLevRecomb[ipISO][nelem] levels saved in RRCoef[nelem][level][temp] */
static multi_arr<double,4> RRCoef/*[LIMELM][NumLevRecomb[ipISO][nelem]][N_ISO_TE_RECOMB]*/;
static multi_arr<long, 2> NumLevRecomb;
static multi_arr<double,3> TotalRecomb;	/*[ipISO][nelem][i]*/

/* the array of logs of temperatures at which RRCoef was defined */
static double TeRRCoef[N_ISO_TE_RECOMB];
static double kTRyd,global_EthRyd; 
static long int globalZ,globalISO;
static long int globalN, globalL, globalS;

STATIC double TempInterp( double* TempArray, double* ValueArray, long NumElements, double temp );
STATIC double iso_recomb_integrand(double EE);
STATIC void iso_put_recomb_error( long ipISO, long nelem );

STATIC double iso_radrecomb_from_cross_section(long ipISO, double temp, long nelem, long ipLo)
{
	double alpha,RecomIntegral=0.,b,E1,E2,step,OldRecomIntegral,TotChangeLastFive;
	double change[5] = {0.,0.,0.,0.,0.};

	DEBUG_ENTRY( "iso_radrecomb_from_cross_section()" );

	if( ipISO==ipH_LIKE && ipLo == 0 )
		return t_ADfA::Inst().H_rad_rec(nelem+1,ipLo,temp);

	global_EthRyd = iso_sp[ipISO][nelem].fb[ipLo].xIsoLevNIonRyd;

	/* Factors outside integral in Milne relation.	*/
	b = MILNE_CONST * iso_sp[ipISO][nelem].st[ipLo].g() * powpq(temp,-3,2);

	if( ipISO==ipH_LIKE )
		b /= 2.;
	else if( ipISO==ipHE_LIKE ) 	
		b /= 4.;
		
	/* kT in Rydbergs.	*/
	kTRyd = temp / TE1RYD;
	globalISO = ipISO;
	globalZ = nelem;
	globalN = N_(ipLo);
	globalL = L_(ipLo);
	globalS = S_(ipLo);

	/* Begin integration.	*/
	/* First define characteristic step */
	E1 = global_EthRyd;

	if( ipISO==ipH_LIKE )
		step = MIN2( 0.125*kTRyd, 0.5*E1 );
	else if( ipISO==ipHE_LIKE ) 	
		step = MIN2( 0.25*kTRyd, 0.5*E1 );
	else
		TotalInsanity();
	
	E2 = E1 + step;
	/* Perform initial integration, from threshold to threshold + step.	*/
	RecomIntegral = qg32( E1, E2, iso_recomb_integrand);
	/* Repeat the integration, adding each new result to the total, 
	 * except that the step size is doubled every time, since values away from 
	 * threshold tend to fall off more slowly.	*/
	do
	{
		OldRecomIntegral = RecomIntegral;
		E1 = E2;
		step *= 1.25;
		E2 = E1 + step;
		RecomIntegral += qg32( E1, E2, iso_recomb_integrand);
		change[4] = change[3];
		change[3] = change[2];
		change[2] = change[1];
		change[1] = change[0];
		change[0] = (RecomIntegral - OldRecomIntegral)/RecomIntegral;
		TotChangeLastFive = change[0] + change[1] + change[2] + change[3] + change[4];
	/* Continue integration until the upper limit exceeds 100kTRyd, an arbitrary
	 * point at which the integrand component exp(electron energy/kT) is very small,
	 * making neglible cross-sections at photon energies beyond that point,
	 * OR when the last five steps resulted in less than a 1 percent change.	*/
	} while ( ((E2-global_EthRyd) < 100.*kTRyd) && ( TotChangeLastFive > 0.0001) );

	/* Calculate recombination coefficient.	*/
	alpha = b * RecomIntegral;

	alpha = MAX2( alpha, SMALLDOUBLE );

	return alpha;
}

/*iso_recomb_integrand, used in Milne relation for iso sequences - the energy is photon Rydbergs.	*/
STATIC double iso_recomb_integrand(double ERyd)
{
	double x1, temp;

	/* Milne relation integrand	*/
	x1 = ERyd * ERyd * exp(-1.0 * ( ERyd - global_EthRyd ) / kTRyd);
	temp = iso_cross_section( ERyd , global_EthRyd, globalN, globalL, globalS, globalZ, globalISO );
	x1 *= temp;

	return x1;
}

double iso_cross_section( double EgammaRyd , double EthRyd, long n, long l, long S, long globalZ, long globalISO )
{
	double cross_section;
	DEBUG_ENTRY( "iso_cross_section()" );

	if( globalISO == ipH_LIKE )
		cross_section = H_cross_section( EgammaRyd , EthRyd, n, l, globalZ );
	else if( globalISO == ipHE_LIKE )
		cross_section =  He_cross_section( EgammaRyd , EthRyd, n, l, S, globalZ );	
	else
		TotalInsanity();

	return cross_section;
}

/*=======================================================*/
/* iso_radiative_recomb get rad recomb rate coefficients for iso sequences */
void iso_radiative_recomb(
						  long ipISO,
						  /* nelem on the c scale, He is 1 */
						  long nelem )
{
	long ipFirstCollapsed, LastN=0L, ThisN=1L, ipLevel; 
	double topoff, TotMinusExplicitResolved,
		TotRRThisN=0., TotRRLastN=0., Total_DR_Added=0.;
	double RecExplictLevels, TotalRadRecomb, RecCollapsed;
	static double TeUsed[NISO][LIMELM]={
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
		 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
		 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
		 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
		 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};

	DEBUG_ENTRY( "iso_radiative_recomb()" );

	iso_put_recomb_error( ipISO, nelem );

	/* evaluate recombination escape probability for all levels */

	/* define radiative recombination rates for all levels */ 
	/* this will be the sum of all levels explicitly included in the model atom */
	RecExplictLevels = 0.;
	
	/* number of resolved levels, so first collapsed level is [ipFirstCollapsed] */
	ipFirstCollapsed = iso_sp[ipISO][nelem].numLevels_local-iso_sp[ipISO][nelem].nCollapsed_local;

	ASSERT( ipFirstCollapsed == iso_sp[ipISO][nelem].numLevels_local - iso_sp[ipISO][nelem].nCollapsed_local );
	if( !fp_equal(phycon.te,TeUsed[ipISO][nelem]) || iso_sp[ipISO][nelem].lgMustReeval || !conv.nTotalIoniz || cosmology.lgDo )
	{
		TeUsed[ipISO][nelem] = phycon.te;

		for( ipLevel=0; ipLevel<ipFirstCollapsed; ++ipLevel )
		{
			/* this is radiative recombination rate coefficient */
			double RadRec;

			if( !iso_ctrl.lgNoRecombInterp[ipISO] )
			{
				RadRec = iso_RRCoef_Te( ipISO, nelem, phycon.te, ipLevel );
			}
			else
			{
				RadRec = iso_radrecomb_from_cross_section( ipISO, phycon.te, nelem, ipLevel);
			}
			ASSERT( RadRec > 0. );

			iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad] = RadRec;

			ASSERT( iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad] > 0. );

			RecExplictLevels += iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad];

			if( iso_ctrl.lgDielRecom[ipISO] )
			{
				iso_sp[ipISO][nelem].fb[ipLevel].DielecRecomb = iso_dielec_recomb_rate( ipISO, nelem, ipLevel )
						* ionbal.DR_Badnell_suppress_fact[nelem][nelem-ipISO];
				Total_DR_Added	+= iso_sp[ipISO][nelem].fb[ipLevel].DielecRecomb;
			}
		}

		/**************************************************/
		/***  Add on recombination to collapsed levels  ***/
		/**************************************************/
		RecCollapsed = 0.;
		for( ipLevel=ipFirstCollapsed; ipLevel<iso_sp[ipISO][nelem].numLevels_local; ++ipLevel )
		{
			/* use hydrogenic for collapsed levels */
			double RadRec = t_ADfA::Inst().H_rad_rec(nelem+1-ipISO, N_(ipLevel), phycon.te);

			/* this is radiative recombination rate coefficient */
			iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad] = RadRec;

			RecCollapsed += iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad];

			ASSERT( iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad] > 0. );

			if( iso_ctrl.lgDielRecom[ipISO] )
			{
				iso_sp[ipISO][nelem].fb[ipLevel].DielecRecomb = iso_dielec_recomb_rate( ipISO, nelem, ipLevel )
						* ionbal.DR_Badnell_suppress_fact[nelem][nelem-ipISO];
				Total_DR_Added	+= iso_sp[ipISO][nelem].fb[ipLevel].DielecRecomb;
			}
		}
		for( ipLevel=iso_sp[ipISO][nelem].numLevels_local; ipLevel<iso_sp[ipISO][nelem].numLevels_max;++ipLevel )
		{
			iso_sp[ipISO][nelem].fb[ipLevel].DielecRecomb = 0.;
		}
		/* >>chng 06 aug 17, from numLevels_max to numLevels_local */
		for( ipLevel = 0; ipLevel<iso_sp[ipISO][nelem].numLevels_local; ipLevel++ )
		{
			if( N_(ipLevel) == ThisN )
			{
				TotRRThisN += iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad];
			}
			else
			{
				ASSERT( iso_ctrl.lgRandErrGen[ipISO] || (TotRRThisN<TotRRLastN) || (ThisN<=2L) || (phycon.te>3E7) || (nelem!=ipHELIUM) || (ThisN==iso_sp[ipISO][nelem].n_HighestResolved_local+1) );
				LastN = ThisN;
				ThisN = N_(ipLevel);
				TotRRLastN = TotRRThisN;
				TotRRThisN = iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad];

				{
					/* Print the sum of recombination coefficients per n at current temp.	*/
					/*@-redef@*/
					enum {DEBUG_LOC=false};
					/*@+redef@*/
					static long RUNONCE = false;

					if( !RUNONCE && DEBUG_LOC )
					{
						static long FIRSTTIME = true;

						if( FIRSTTIME )
						{
							fprintf( ioQQQ,"Sum of Radiative Recombination at current iso, nelem, temp = %li %li %.2f\n", 
								ipISO, nelem, phycon.te);
							FIRSTTIME= false;
						}

						fprintf( ioQQQ,"%li\t%.2e\n",LastN,TotRRLastN );
					}
					RUNONCE = true;
				}
			}
		}

		/* Get total recombination into all levels, including those not explicitly considered. */
		if( iso_sp[ipISO][nelem].lgLevelsLowered )
		{
			TotalRadRecomb = 0.;
			for( ipLevel = 0; ipLevel<iso_sp[ipISO][nelem].numLevels_local; ipLevel++ )
				TotalRadRecomb += iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad];
		}
		else
		{
			// Not lowered, choose between fits for total for all levels
			if( iso_ctrl.lgNoRecombInterp[ipISO] )
			{
				/* We are not interpolating, must calculate total now, as sum of resolved and collapsed levels... */
				TotalRadRecomb = RecCollapsed + RecExplictLevels;
				
				/* Plus additional levels out to a predefined limit... */
				for( long nLo = N_(ipFirstCollapsed) + iso_sp[ipISO][nelem].nCollapsed_max; nLo < NHYDRO_MAX_LEVEL; nLo++ )
				{
					TotalRadRecomb += t_ADfA::Inst().H_rad_rec(nelem+1-ipISO, nLo, phycon.te);
				}
				/* Plus a bunch more for good measure */
				for( long nLo = NHYDRO_MAX_LEVEL; nLo<=SumUpToThisN; nLo++ )
				{
					TotalRadRecomb += Recomb_Seaton59( nelem+1-ipISO, phycon.te, nLo );
				}
			}
			else 
			{
				/* We are interpolating, and total was calculated here in iso_recomb_setup */
				TotalRadRecomb = iso_RRCoef_Te( ipISO, nelem, phycon.te, iso_sp[ipISO][nelem].numLevels_max - iso_sp[ipISO][nelem].nCollapsed_max );
			}
		}

		if(ipISO==0 && nelem==0 )
		{
			// insure rec coef will be always evaluated
			TeUsed[ipISO][nelem] = phycon.te*0.999;
		}

		/* If generating random error, apply one to total recombination */
		if( iso_ctrl.lgRandErrGen[ipISO] )
		{
			/* Error for total recombination */
			/* this has to be from iso.numLevels_max instead of iso.numLevels_local because
			 * the error factors for rrc are always stored at iso.numLevels_max, regardless of
			 * continuum lowering effects. */
			iso_put_error(ipISO,nelem,iso_sp[ipISO][nelem].numLevels_max,
				      iso_sp[ipISO][nelem].numLevels_max,IPRAD,0.0001f,0.0001f);
		}

		/* this is case B recombination, sum without the ground included */
		iso_sp[ipISO][nelem].RadRec_caseB = TotalRadRecomb - iso_sp[ipISO][nelem].fb[0].RadRecomb[ipRecRad];
		ASSERT( iso_sp[ipISO][nelem].RadRec_caseB > 0.);

		// Restore highest levels opacity stack
		for( unsigned i = 0; i < iso_sp[ipISO][nelem].HighestLevelOpacStack.size(); ++i )
		{
			long index = iso_sp[ipISO][nelem].numLevels_max-1;
			opac.OpacStack[iso_sp[ipISO][nelem].fb[index].ipOpac-1+i] = iso_sp[ipISO][nelem].HighestLevelOpacStack[i];
		}

		/* Add topoff (excess) recombination to top level.  This is only done if atom is full size,
		 * that is, no pressure lowered of the continuum the current conditions.  Radiative recombination
		 * to non-existent states does not occur, those would dominate the topoff */
		if( !iso_sp[ipISO][nelem].lgLevelsLowered )
		{
			/* at this point we have RecExplictLevels, the sum of radiative recombination 
			 * to all levels explicitly included in the model atom the total 
			 * recombination rate.  The difference is the amount of "topoff" we will need to do */
			TotMinusExplicitResolved = TotalRadRecomb - RecExplictLevels;

			topoff = TotMinusExplicitResolved - RecCollapsed;

			/* the t_ADfA::Inst().rad_rec fits are too small at high temperatures, so this atom is
			 * better than the topoff.  Only a problem for helium itself, at high temperatures.
			 * complain if the negative topoff is not for this case */
			if( topoff < 0. && (nelem!=ipHELIUM || phycon.te < 1e5 ) &&
				fabs( topoff/TotalRadRecomb ) > 0.01 )
			{
				fprintf(ioQQQ," PROBLEM  negative RR topoff for %li, rel error was %.2e, temp was %.2f\n",  
					nelem, topoff/TotalRadRecomb, phycon.te );
			}

			/* option to turn off topoff with atom xx-like topoff off, or if levels
			 * into our model, so topoff is not needed since we have a complete model */
			if( !iso_ctrl.lgTopoff[ipISO] || iso_sp[ipISO][nelem].lgLevelsLowered )
				topoff = 0;

			topoff = MAX2( 0., topoff );
			double scale_factor = 1. + topoff/iso_sp[ipISO][nelem].fb[iso_sp[ipISO][nelem].numLevels_max-1].RadRecomb[ipRecRad];
			ASSERT( scale_factor >= 1. );

			// Scale highest level opacities to be consistent with recombination topoff
			for( unsigned i = 0; i < iso_sp[ipISO][nelem].HighestLevelOpacStack.size(); ++i )
			{
				long index = iso_sp[ipISO][nelem].numLevels_max-1;
				opac.OpacStack[iso_sp[ipISO][nelem].fb[index].ipOpac-1+i] *= scale_factor;
			}
		
			/* We always have at least one collapsed level if continuum is not lowered.  Put topoff there.	*/
			iso_sp[ipISO][nelem].fb[iso_sp[ipISO][nelem].numLevels_max-1].RadRecomb[ipRecRad] += topoff;

			/* check for negative DR topoff, but only if Total_DR_Added is not negligible compared with TotalRadRecomb */
			if( Total_DR_Added > TotalRadRecomb/100. )
			{
				if( ipISO == ipHE_LIKE &&
					ionbal.DR_Badnell_rate_coef[nelem][nelem-ipISO] > 0. &&
					Total_DR_Added > 1.02 * ionbal.DR_Badnell_rate_coef[nelem][nelem-ipISO] )
				{
					fprintf(ioQQQ," PROBLEM  negative DR topoff for %li iso %li, tot1, tot2 = %.2e\t%.2e rel error was %.2e, temp was %.2e, eden was %.2e\n",
						nelem,
						ipISO,
						Total_DR_Added,
						ionbal.DR_Badnell_rate_coef[nelem][nelem-ipISO],
						Total_DR_Added / ionbal.DR_Badnell_rate_coef[nelem][nelem-ipISO] - 1.0,
						phycon.te,
						dense.eden );
				}
			}

			ASSERT( iso_sp[ipISO][nelem].numLevels_max == iso_sp[ipISO][nelem].numLevels_local );

			if( iso_ctrl.lgDielRecom[ipISO] && iso_ctrl.lgTopoff[ipISO] )
			{
				/* \todo 2 suppress this total rate for continuum lowering using factors from Jordan (1969). */
				/* put extra DR in top level */
				iso_sp[ipISO][nelem].fb[iso_sp[ipISO][nelem].numLevels_max-1].DielecRecomb +=
					MAX2( 0., ionbal.DR_Badnell_rate_coef[nelem][nelem-ipISO] - Total_DR_Added );
			}
		}
	}

	/**************************************************************/
	/***  Stuff escape probabilities, and effective rad recomb  ***/
	/**************************************************************/

	/* total effective radiative recombination, initialize to zero */
	iso_sp[ipISO][nelem].RadRec_effec = 0.;

	for( ipLevel=0; ipLevel<iso_sp[ipISO][nelem].numLevels_local; ++ipLevel )
	{
		/* option for case b conditions, kill ground state recombination */
		if( opac.lgCaseB && ipLevel==0 )
		{
			iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecEsc] = 1e-10;
			iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecNetEsc] = 1e-10;
		}
		else if( cosmology.lgDo && ipLevel==0 )
		{
			iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecEsc] = 0.;
			iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecNetEsc] = 0.;
		}
		else
		{
			iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecEsc] =
				RT_recom_effic(iso_sp[ipISO][nelem].fb[ipLevel].ipIsoLevNIonCon);

			/* net escape prob includes dest by background opacity */
			iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecNetEsc] =
				MIN2( 1.,
				iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecEsc] +
				(1.-iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecEsc])*
				iso_sp[ipISO][nelem].fb[ipLevel].ConOpacRatio );
		}

		ASSERT( iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecEsc] >= 0. );
		ASSERT( iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecNetEsc] >= 0. );
		ASSERT( iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecNetEsc] <= 1. );

		/* sum of all effective rad rec */
		iso_sp[ipISO][nelem].RadRec_effec += iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad]*
		  iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecNetEsc];
	}

	/* zero out escape probabilities of levels above numLevels_local */
	for( ipLevel=iso_sp[ipISO][nelem].numLevels_local; ipLevel<iso_sp[ipISO][nelem].numLevels_max; ++ipLevel )
	{
		/* this is escape probability */
		iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecEsc] = 0.;
		/* net escape prob includes dest by background opacity */
		iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecNetEsc] = 0.;
	}

	/* trace escape probabilities */
	if( trace.lgTrace && trace.lgIsoTraceFull[ipISO] && (nelem == trace.ipIsoTrace[ipISO]) )
	{
		fprintf( ioQQQ, "       iso_radiative_recomb trace ipISO=%3ld Z=%3ld\n", ipISO, nelem );

		/* print continuum escape probability */
		fprintf( ioQQQ, "       iso_radiative_recomb recomb effic" );
		for( ipLevel=0; ipLevel < iso_sp[ipISO][nelem].numLevels_local; ipLevel++ )
		{
			fprintf( ioQQQ,PrintEfmt("%10.3e", iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecEsc] ));
		}
		fprintf( ioQQQ, "\n" );

		/* net recombination efficiency factor, including background opacity*/
		fprintf( ioQQQ, "       iso_radiative_recomb recomb net effic" );
		for( ipLevel=0; ipLevel < iso_sp[ipISO][nelem].numLevels_local; ipLevel++ )
		{
			fprintf( ioQQQ,PrintEfmt("%10.3e", iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecNetEsc]) );
		}
		fprintf( ioQQQ, "\n" );

		/* inward continuous optical depths */
		fprintf( ioQQQ, "       iso_radiative_recomb in optic dep" );
		for( ipLevel=0; ipLevel < iso_sp[ipISO][nelem].numLevels_local; ipLevel++ )
		{	
			fprintf( ioQQQ,PrintEfmt("%10.3e",  opac.TauAbsGeo[0][iso_sp[ipISO][nelem].fb[ipLevel].ipIsoLevNIonCon-1] ));
		}
		fprintf( ioQQQ, "\n" );

		/* outward continuous optical depths*/
		fprintf( ioQQQ, "       iso_radiative_recomb out op depth" );
		for( ipLevel=0; ipLevel < iso_sp[ipISO][nelem].numLevels_local; ipLevel++ )
		{
			fprintf( ioQQQ,PrintEfmt("%10.3e",  opac.TauAbsGeo[1][iso_sp[ipISO][nelem].fb[ipLevel].ipIsoLevNIonCon-1] ));
		}
		fprintf( ioQQQ, "\n" );

		/* print radiative recombination coefficients */
		fprintf( ioQQQ, "       iso_radiative_recomb rad rec coef " );
		for( ipLevel=0; ipLevel < iso_sp[ipISO][nelem].numLevels_local; ipLevel++ )
		{
			fprintf( ioQQQ,PrintEfmt("%10.3e", iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad]) );
		}
		fprintf( ioQQQ, "\n" );
	}

	if( trace.lgTrace && (trace.lgHBug||trace.lgHeBug) )
	{
		/* print total recombination only one time per temperature so that this statement can be
		 * enabled and produce useful tables of recombination vs Te
		 */
		static double Tused = -1;
		if( Tused != phycon.te )
		{
			Tused = phycon.te;
			fprintf( ioQQQ, "     iso_radiative_recomb total rec coef Te %10.3e", phycon.te );
			fprintf( ioQQQ,PrintEfmt("%10.3e", iso_sp[ipISO][nelem].RadRec_effec ));
			fprintf( ioQQQ, " case A=" );
			fprintf( ioQQQ,PrintEfmt("%10.3e",
				iso_sp[ipISO][nelem].RadRec_caseB + iso_sp[ipISO][nelem].fb[ipH1s].RadRecomb[ipRecRad] ) );
			fprintf( ioQQQ, " case B=");
			fprintf( ioQQQ,PrintEfmt("%10.3e", iso_sp[ipISO][nelem].RadRec_caseB ));
			fprintf( ioQQQ, "\n" );
		}
	}

	/****************************/
	/***  begin sanity check  ***/
	/****************************/
	{
		bool lgOK = true;
		for( ipLevel=0; ipLevel < iso_sp[ipISO][nelem].numLevels_local; ipLevel++ )
		{
			if( iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad] <= 0. )
			{
				fprintf( ioQQQ, 
					" PROBLEM iso_radiative_recomb non-positive recombination coefficient for ipISO=%3ld Z=%3ld lev n=%3ld rec=%11.2e te=%11.2e\n", 
				  ipISO, nelem, ipLevel, iso_sp[ipISO][nelem].fb[ipLevel].RadRecomb[ipRecRad] , phycon.te);
					lgOK = false;
			}
		}
		/* bail if we found problems */
		if( !lgOK )
		{
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}
		/*end sanity check */
	}

	/* confirm that we have good rec coef at bottom and top of atom/ion */
	ASSERT( iso_sp[ipISO][nelem].fb[0].RadRecomb[ipRecRad] > 0. );
	ASSERT( iso_sp[ipISO][nelem].fb[iso_sp[ipISO][nelem].numLevels_local-1].RadRecomb[ipRecRad] > 0. );

	/* set true when save recombination coefficients command entered */
	if( save.lgioRecom )
	{
		/* this prints Z on physical, not C, scale */
		fprintf( save.ioRecom, "%s %s %2li ", 
			iso_ctrl.chISO[ipISO], elementnames.chElementSym[nelem], nelem+1 );
		fprintf( save.ioRecom,PrintEfmt("%9.2e", iso_sp[ipISO][nelem].RadRec_effec ));
		fprintf( save.ioRecom, "\n" );
	}

	return;
}

STATIC void iso_put_recomb_error( long ipISO, long nelem )
{
	/* optimistic and pessimistic errors for HeI recombination */
	static map<QNPack,realnum> He_errorOpti = {
		{QN2ind(1,0,1,1), 0.0000_r}, {QN2ind(2,0,3,3), 0.0009_r},
		{QN2ind(2,0,1,1), 0.0037_r}, {QN2ind(2,1,3,-1), 0.0003_r}, {QN2ind(2,1,1,3), 0.0018_r},
		{QN2ind(3,0,3,3), 0.0009_r}, {QN2ind(3,0,1,1), 0.0050_r}, {QN2ind(3,1,3,-1), 0.0007_r},
		{QN2ind(3,2,3,-1), 0.0003_r}, {QN2ind(3,2,1,5), 0.0001_r}, {QN2ind(3,1,1,3), 0.0007_r},
		{QN2ind(4,0,3,3), 0.0045_r}, {QN2ind(4,0,1,1), 0.0071_r}, {QN2ind(4,1,3,-1), 0.0005_r},
		{QN2ind(4,2,3,-1), 0.0005_r}, {QN2ind(4,2,1,5), 0.0004_r}, {QN2ind(4,3,3,-1), 0.0005_r},
		{QN2ind(4,3,1,7), 0.0004_r}, {QN2ind(4,1,1,3), 0.0009_r}, {QN2ind(5,0,3,3), 0.0045_r},
		{QN2ind(5,0,1,1), 0.0071_r}, {QN2ind(5,1,3,-1), 0.0005_r}, {QN2ind(5,2,3,-1), 0.0005_r},
		{QN2ind(5,2,1,5), 0.0004_r}, {QN2ind(5,3,3,-1), 0.0005_r}, {QN2ind(5,3,1,7), 0.0004_r},
		{QN2ind(5,4,3,-1), 0.0005_r}, {QN2ind(5,4,1,9), 0.0004_r}, {QN2ind(5,1,1,3), 0.0009_r}
	};
	static map<QNPack,realnum> He_errorPess = {
		{QN2ind(1,0,1,1), 0.0100_r}, {QN2ind(2,0,3,3), 0.0100_r},
		{QN2ind(2,0,1,1), 0.0060_r}, {QN2ind(2,1,3,-1), 0.0080_r}, {QN2ind(2,1,1,3), 0.0200_r},
		{QN2ind(3,0,3,3), 0.0200_r}, {QN2ind(3,0,1,1), 0.0200_r}, {QN2ind(3,1,3,-1), 0.0200_r},
		{QN2ind(3,2,3,-1), 0.0600_r}, {QN2ind(3,2,1,5), 0.0600_r}, {QN2ind(3,1,1,3), 0.0080_r},
		{QN2ind(4,0,3,3), 0.0200_r}, {QN2ind(4,0,1,1), 0.0200_r}, {QN2ind(4,1,3,-1), 0.0070_r},
		{QN2ind(4,2,3,-1), 0.0100_r}, {QN2ind(4,2,1,5), 0.0100_r}, {QN2ind(4,3,3,-1), 0.0020_r},
		{QN2ind(4,3,1,7), 0.0030_r}, {QN2ind(4,1,1,3), 0.0070_r}, {QN2ind(5,0,3,3), 0.0300_r},
		{QN2ind(5,0,1,1), 0.0300_r}, {QN2ind(5,1,3,-1), 0.0100_r}, {QN2ind(5,2,3,-1), 0.0200_r},
		{QN2ind(5,2,1,5), 0.0200_r}, {QN2ind(5,3,3,-1), 0.0200_r}, {QN2ind(5,3,1,7), 0.0200_r},
		{QN2ind(5,4,3,-1), 0.0010_r}, {QN2ind(5,4,1,9), 0.0004_r}, {QN2ind(5,1,1,3), 0.0090_r}
	};

	/* now put recombination errors into iso.Error[ipISO] array */
	for( long ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
	{
		if( ipISO==ipHE_LIKE && nelem==ipHELIUM )
		{
			long n, l, s;
			if( ipHi >= iso_sp[ipISO][nelem].numLevels_max - iso_sp[ipISO][nelem].nCollapsed_max )
			{
				n = 5; // treat collapsed levels as 5 3S
				l = 0;
				s = 3;
			}
			else
			{
				n = min(iso_sp[ipISO][nelem].st[ipHi].n(), 5);
				l = min(iso_sp[ipISO][nelem].st[ipHi].l(), 4);
				s = iso_sp[ipISO][nelem].st[ipHi].S();
			}
			QNPack ind = QN2ind(n, l ,s);
			auto p_opt = He_errorOpti.find(ind);
			auto p_pes = He_errorPess.find(ind);

			ASSERT( p_opt != He_errorOpti.end() && p_pes != He_errorPess.end() );

			iso_put_error(ipISO,nelem,iso_sp[ipISO][nelem].numLevels_max,ipHi,IPRAD,p_opt->second,p_pes->second);
		}
		else
			iso_put_error(ipISO,nelem,iso_sp[ipISO][nelem].numLevels_max,ipHi,IPRAD,0.1_r,0.1_r);
	}
}

void iso_radiative_recomb_effective( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_radiative_recomb_effective()" );

	/* Find effective recombination coefficients */
	for( long ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
	{
		iso_sp[ipISO][nelem].fb[ipHi].RadEffec = 0.;

		/* >>chng 06 aug 17, from numLevels_max to numLevels_local */
		for( long ipHigher=ipHi; ipHigher < iso_sp[ipISO][nelem].numLevels_local; ipHigher++ )
		{
			ASSERT( iso_sp[ipISO][nelem].CascadeProb[ipHigher][ipHi] >= 0. );
			ASSERT( iso_sp[ipISO][nelem].fb[ipHigher].RadRecomb[ipRecRad] >= 0. );

			iso_sp[ipISO][nelem].fb[ipHi].RadEffec += iso_sp[ipISO][nelem].CascadeProb[ipHigher][ipHi] *
				iso_sp[ipISO][nelem].fb[ipHigher].RadRecomb[ipRecRad];
		}
	}

	/**************************************************************/
	/***  option to print effective recombination coefficients  ***/
	/**************************************************************/
	{
		enum {DEBUG_LOC=false};

		if( DEBUG_LOC )
		{
			const int maxPrt=10;

			fprintf( ioQQQ,"Effective recombination, ipISO=%li, nelem=%li, Te = %e\n", ipISO, nelem, phycon.te );
			fprintf( ioQQQ, "N\tL\tS\tRadEffec\tLifetime\n" );

			for( long i=0; i<maxPrt; i++ )
			{
				fprintf( ioQQQ, "%li\t%li\t%li\t%e\t%e\n", N_(i), L_(i), S_(i),
					iso_sp[ipISO][nelem].fb[i].RadEffec,
					MAX2( 0., iso_sp[ipISO][nelem].st[i].lifetime() ) );
			}
			fprintf( ioQQQ, "\n" );
		}
	}

	/* If we have the variable set, find errors in rad rates */
	if( iso_ctrl.lgRandErrGen[ipISO] )
	{
		dprintf( ioQQQ, "ipHi\tipLo\tWL\tEmiss\tSigmaEmiss\tRadEffec\tSigRadEff\tBrRat\tSigBrRat\n" );

		/* >>chng 06 aug 17, from numLevels_max to numLevels_local */
		for( long ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
		{
			iso_sp[ipISO][nelem].fb[ipHi].SigmaRadEffec = 0.;

			/* >>chng 06 aug 17, from numLevels_max to numLevels_local */
			for( long ipHigher=ipHi; ipHigher < iso_sp[ipISO][nelem].numLevels_local; ipHigher++ )
			{
				ASSERT( iso_sp[ipISO][nelem].ex[iso_sp[ipISO][nelem].numLevels_max][ipHigher].Error[IPRAD] >= 0. );
				ASSERT( iso_sp[ipISO][nelem].ex[ipHigher][ipHi].SigmaCascadeProb >= 0. );

				/* iso.RadRecomb has to appear here because iso.Error is only relative error */ 
				iso_sp[ipISO][nelem].fb[ipHi].SigmaRadEffec += pow2( iso_sp[ipISO][nelem].ex[iso_sp[ipISO][nelem].numLevels_max][ipHigher].Error[IPRAD] *
					iso_sp[ipISO][nelem].CascadeProb[ipHigher][ipHi] * iso_sp[ipISO][nelem].fb[ipHigher].RadRecomb[ipRecRad]) +
					pow2( iso_sp[ipISO][nelem].ex[ipHigher][ipHi].SigmaCascadeProb * iso_sp[ipISO][nelem].fb[ipHigher].RadRecomb[ipRecRad]);
			}

			ASSERT( iso_sp[ipISO][nelem].fb[ipHi].SigmaRadEffec >= 0. );
			iso_sp[ipISO][nelem].fb[ipHi].SigmaRadEffec = sqrt( iso_sp[ipISO][nelem].fb[ipHi].SigmaRadEffec );

			for( long ipLo = 0; ipLo < ipHi; ipLo++ )
			{
				if( (( L_(ipLo) == L_(ipHi) + 1 ) || ( L_(ipHi) == L_(ipLo) + 1 )) )
				{	
					double EnergyInRydbergs = iso_sp[ipISO][nelem].fb[ipLo].xIsoLevNIonRyd - iso_sp[ipISO][nelem].fb[ipHi].xIsoLevNIonRyd;
					realnum wavelength = (realnum)(RYDLAM/MAX2(1E-8,EnergyInRydbergs));
					double emissivity = iso_sp[ipISO][nelem].fb[ipHi].RadEffec * iso_sp[ipISO][nelem].BranchRatio[ipHi][ipLo] * EN1RYD * EnergyInRydbergs;
					double sigma_emiss = 0., SigmaBranchRatio = 0.;

					if( ( emissivity > 2.E-29 ) && ( wavelength < 1.E6 ) && (N_(ipHi)<=5) )
					{
						SigmaBranchRatio = iso_sp[ipISO][nelem].BranchRatio[ipHi][ipLo] * sqrt(
							pow2( (double)iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[IPRAD] ) +
							pow2( iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot*iso_sp[ipISO][nelem].st[ipHi].lifetime() ) );

						sigma_emiss =  EN1RYD * EnergyInRydbergs * sqrt( 
							pow2( (double)iso_sp[ipISO][nelem].fb[ipHi].SigmaRadEffec * iso_sp[ipISO][nelem].BranchRatio[ipHi][ipLo] ) +
							pow2( SigmaBranchRatio * iso_sp[ipISO][nelem].fb[ipHi].RadEffec ) );

						/* \todo 2 make this a trace option */
						dprintf( ioQQQ, "%li\t%li\t", ipHi, ipLo );
						prt_wl( ioQQQ, wavelength );
						fprintf( ioQQQ, "\t%e\t%e\t%e\t%e\t%e\t%e\n", 
							emissivity,
							sigma_emiss,
							iso_sp[ipISO][nelem].fb[ipHi].RadEffec,
							iso_sp[ipISO][nelem].fb[ipHi].SigmaRadEffec,
							iso_sp[ipISO][nelem].BranchRatio[ipHi][ipLo],
							SigmaBranchRatio);
					}
				}
			}
		}
	}

	return;
}
/*iso_RRCoef_Te evaluated radiative recombination coef at some temperature */
double iso_RRCoef_Te( long ipISO, long nelem, double temp, long n )
{
	double rate = 0.;

	DEBUG_ENTRY( "iso_RRCoef_Te()" );

	ASSERT( !iso_ctrl.lgNoRecombInterp[ipISO] );

	/* if n is equal to the number of levels, return the total recomb, else recomb for given level.	*/
	if( n == iso_sp[ipISO][nelem].numLevels_max - iso_sp[ipISO][nelem].nCollapsed_max )
	{
		rate = TempInterp( TeRRCoef, &TotalRecomb[ipISO][nelem][0], N_ISO_TE_RECOMB, temp );
	}
	else
	{
		rate = TempInterp( TeRRCoef, &RRCoef[ipISO][nelem][n][0], N_ISO_TE_RECOMB, temp );
	}

	/* that was the log, now make linear */
	rate = exp10( rate );

	return rate;
}

/*iso_recomb_check called by SanityCheck to confirm that recombination coef are ok
 * return value is relative error between new calculation of recom, and interp value */
double iso_recomb_check( long ipISO, long nelem, long level, double temperature )
{
	double RecombRelError ,
		RecombInterp,
		RecombCalc;

	DEBUG_ENTRY( "iso_recomb_check()" );

	/* actually calculate the recombination coefficient from the Milne relation,
	 * normally only done due to compile he-like command */
	RecombCalc = iso_radrecomb_from_cross_section( ipISO, temperature , nelem , level );

	/* interpolate the recombination coefficient, this is the usual method */
	RecombInterp = iso_RRCoef_Te( ipISO, nelem, temperature, level );

	RecombRelError = ( RecombInterp - RecombCalc ) / MAX2( RecombInterp , RecombCalc );

	return RecombRelError;
}

/* allocate space needed for iso recombination tables */
void iso_recomb_alloc()
{
	DEBUG_ENTRY( "iso_recomb_alloc()" );

	/* The number of recombination coefficients to be read from file for each element.	*/
	NumLevRecomb.alloc(NISO, LIMELM);
	TotalRecomb.alloc(NISO, LIMELM, N_ISO_TE_RECOMB);

	RRCoef.reserve(NISO);
	for( long ipISO=0; ipISO<NISO; ipISO++ )
	{
		RRCoef.reserve(ipISO, LIMELM);
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			long int MaxLevels, maxN;

			if( nelem == ipISO )
				maxN = RREC_MAXN;
			else
				maxN = LIKE_RREC_MAXN( nelem );

			NumLevRecomb[ipISO][nelem] = iso_get_total_num_levels( ipISO, maxN, 0 );

			if( nelem == ipISO || dense.lgElmtOn[nelem] )
			{
				/* must always have at least NumLevRecomb[ipISO][nelem] levels since that is number 
				* that will be read in from he rec data file, but possible to require more */
				MaxLevels = MAX2( NumLevRecomb[ipISO][nelem] , iso_sp[ipISO][nelem].numLevels_max );

				/* always define this */
				/* >>chng 02 jan 24, RRCoef will be iso_sp[ipISO][nelem].numLevels_max levels, not iso.numLevels_max,
				* code will stop if more than this is requested */
				RRCoef.reserve(ipISO, nelem, MaxLevels);

				for( long ipLo=0; ipLo < MaxLevels;++ipLo )
				{
					RRCoef.reserve(ipISO, nelem, ipLo, N_ISO_TE_RECOMB);
				}
			}
		}
	}
	RRCoef.alloc();

	for(long i = 0; i < N_ISO_TE_RECOMB; i++)
	{
		/* this is the vector of temperatures */
		TeRRCoef[i] = 0.25*(i);
	}

	/* >>chng 06 jun 06, NP, assert thrown at T == 1e10 K, just bump the 
	 * high temperature end slightly.  */
	TeRRCoef[N_ISO_TE_RECOMB-1] += 0.01f;

	return;
}

void iso_recomb_auxiliary_free()
{
	DEBUG_ENTRY( "iso_recomb_auxiliary_free()" );

	NumLevRecomb.clear();

	return;
}

void iso_recomb_setup( long ipISO )
{
	double RadRecombReturn;
	long int i, i1, i2, i3, i4, i5;
	long int ipLo, nelem;

	const char* chFilename[NISO] = { "h_iso_recomb.dat", "he_iso_recomb.dat" };

	bool lgEOL;

	DEBUG_ENTRY( "iso_recomb_setup()" );

	/* if we are compiling the recombination data file, we must interpolate in temperature */
	if( iso_ctrl.lgCompileRecomb[ipISO] )
	{
		iso_ctrl.lgNoRecombInterp[ipISO] = false;
	}

	if( !iso_ctrl.lgNoRecombInterp[ipISO] )
	{
		/******************************************************************/
		/**  Establish radiative recombination rate coefficients - RRC	***/
		/******************************************************************/
		/* This flag says we are not compiling the data file	*/
		if( !iso_ctrl.lgCompileRecomb[ipISO] )
		{
			if( trace.lgTrace )
				fprintf( ioQQQ," iso_recomb_setup opening %s:", chFilename[ipISO] );

			/* Now try to read from file...*/
			FILE *ioDATA = open_data( chFilename[ipISO], "r", AS_OPTIONAL );
			if( ioDATA == NULL )
			{
				fprintf( ioQQQ, " Defaulting to on-the-fly computation," );
				fprintf( ioQQQ, " but the code runs much faster if you compile this file!\n" );
				for( nelem = ipISO; nelem < LIMELM; nelem++ )
				{
					if( dense.lgElmtOn[nelem] )
					{
						/* Zero out the recombination sum array.	*/
						for(i = 0; i < N_ISO_TE_RECOMB; i++)
						{
							TotalRecomb[ipISO][nelem][i] = 0.;
						}

						/* NumLevRecomb[ipISO][nelem] corresponds to n = 40 for H and He and 20 for ions, at present	*/
						/* There is no need to fill in values for collapsed levels, because we do not need to
						* interpolate for a given temperature, just calculate it directly with a hydrogenic routine.	*/
						for( ipLo=0; ipLo < iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max; ipLo++ )
						{
							/* loop over temperatures to produce array of recombination coefficients	*/
							for(i = 0; i < N_ISO_TE_RECOMB; i++)
							{
								/* Store log of recombination coefficients, in N_ISO_TE_RECOMB half dec steps */
								RadRecombReturn = iso_radrecomb_from_cross_section( ipISO, exp10( TeRRCoef[i] ) ,nelem,ipLo);
								TotalRecomb[ipISO][nelem][i] += RadRecombReturn;
								RRCoef[ipISO][nelem][ipLo][i] = log10(RadRecombReturn);
							}
						}
						for(i = 0; i < N_ISO_TE_RECOMB; i++)
						{
							for( i1 = iso_sp[ipISO][nelem].n_HighestResolved_max+1; i1< NHYDRO_MAX_LEVEL; i1++ )
							{
								TotalRecomb[ipISO][nelem][i] += t_ADfA::Inst().H_rad_rec(nelem+1-ipISO,i1, exp10(TeRRCoef[i]));
							}
							for( i1 = NHYDRO_MAX_LEVEL; i1<=SumUpToThisN; i1++ )
							{
								TotalRecomb[ipISO][nelem][i] += Recomb_Seaton59( nelem+1-ipISO, exp10(TeRRCoef[i]), i1 );
							}
							TotalRecomb[ipISO][nelem][i] = log10( TotalRecomb[ipISO][nelem][i] );
						}
					}
				}
			}
			/* Data file is present and readable...begin read.	*/
			else 
			{
				/* check that magic number is ok */
				string chLine;
				if( !read_whole_line( chLine, ioDATA ) )
				{
					fprintf( ioQQQ, " iso_recomb_setup could not read first line of %s.\n", chFilename[ipISO]);
					cdEXIT(EXIT_FAILURE);
				}
				i = 1;
				i1 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
				i2 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
				i3 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
				i4 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
				if( i1 !=RECOMBMAGIC || i2 !=NumLevRecomb[ipISO][ipISO] || i3 !=NumLevRecomb[ipISO][ipISO+1] || i4 !=N_ISO_TE_RECOMB )
				{
					fprintf( ioQQQ, 
						" iso_recomb_setup: the version of %s is not the current version.\n", chFilename[ipISO] );
					fprintf( ioQQQ, 
						" iso_recomb_setup: I expected to find the numbers  %i %li %li %i and got %li %li %li %li instead.\n" ,
						RECOMBMAGIC ,
						NumLevRecomb[ipISO][ipISO],
						NumLevRecomb[ipISO][ipISO+1],
						N_ISO_TE_RECOMB,
						i1 , i2 , i3, i4 );
					fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine.c_str() );
					fprintf( ioQQQ, 
						" iso_recomb_setup: please recompile the data file with the COMPile RECOmb COEFficient H-LIke [or HE-Like] command.\n" );
					cdEXIT(EXIT_FAILURE);
				}

				i5 = 1;
				/* now read in the data */
				for( nelem = ipISO; nelem < LIMELM; nelem++ )
				{
					for( ipLo=0; ipLo <= NumLevRecomb[ipISO][nelem]; ipLo++ )
					{
						i5++;
						/* get next line image */
						if( !read_whole_line( chLine, ioDATA ) )
						{
							fprintf( ioQQQ, " iso_recomb_setup could not read line %li of %s.\n", i5,
								chFilename[ipISO] );
							cdEXIT(EXIT_FAILURE);
						}
						/* each line starts with element and level number */
						i3 = 1;
						i1 = (long)FFmtRead(chLine.c_str(),&i3,chLine.length(),&lgEOL);
						i2 = (long)FFmtRead(chLine.c_str(),&i3,chLine.length(),&lgEOL);
						/* check that these numbers are correct */
						if( i1!=nelem || i2!=ipLo )
						{
							fprintf( ioQQQ, " iso_recomb_setup detected insanity in %s.\n", chFilename[ipISO]);
							fprintf( ioQQQ, 
								" iso_recomb_setup: please recompile the data file with the COMPile RECOmb COEFficient H-LIke [or HE-Like] command.\n" );
							cdEXIT(EXIT_FAILURE);
						}

						/* loop over temperatures to produce array of recombination coefficients	*/
						for(i = 0; i < N_ISO_TE_RECOMB; i++)
						{
							double ThisCoef = FFmtRead(chLine.c_str(),&i3,chLine.length(),&lgEOL);

							if( nelem == ipISO || dense.lgElmtOn[nelem] )
							{
								/* The last line for each element is the total recombination for each temp.	*/
								if( ipLo == NumLevRecomb[ipISO][nelem] )
								{
									TotalRecomb[ipISO][nelem][i] = ThisCoef;
								}
								else
									RRCoef[ipISO][nelem][ipLo][i] = ThisCoef;
							}

							if( lgEOL )
							{
								fprintf( ioQQQ, " iso_recomb_setup detected insanity in %s.\n", chFilename[ipISO]);
								fprintf( ioQQQ, 
									" iso_recomb_setup: please recompile the data file with the COMPile RECOmb COEFficient H-LIke [or HE-Like] command.\n" );
								cdEXIT(EXIT_FAILURE);
							}
						}
					}

					/* following loop only executed if we need more levels than are
						* stored in the recom coef data set
						* do not do collapsed levels since will use H-like recom there */
					if( nelem == ipISO || dense.lgElmtOn[nelem] ) 
					{
						for( ipLo=NumLevRecomb[ipISO][nelem]; ipLo < iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max; ipLo++ )
						{
								for(i = 0; i < N_ISO_TE_RECOMB; i++)
								{
									/* Store log of recombination coefficients, in N_ISO_TE_RECOMB half dec steps */
									RRCoef[ipISO][nelem][ipLo][i] = log10(iso_radrecomb_from_cross_section( ipISO, exp10( TeRRCoef[i] ) ,nelem,ipLo));
								}
							}
						}
				}

				/* check that ending magic number is ok */
				if( !read_whole_line( chLine, ioDATA ) )
				{
					fprintf( ioQQQ, " iso_recomb_setup could not read last line of %s.\n", chFilename[ipISO]);
					cdEXIT(EXIT_FAILURE);
				}
				i = 1;
				i1 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
				i2 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
				i3 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
				i4 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

				if( i1 !=RECOMBMAGIC || i2 !=NumLevRecomb[ipISO][ipISO] || i3 !=NumLevRecomb[ipISO][ipISO+1] || i4 !=N_ISO_TE_RECOMB )
				{
					fprintf( ioQQQ, 
						" iso_recomb_setup: the version of %s is not the current version.\n", chFilename[ipISO] );
					fprintf( ioQQQ, 
						" iso_recomb_setup: I expected to find the numbers  %i %li %li %i and got %li %li %li %li instead.\n" ,
						RECOMBMAGIC ,
						NumLevRecomb[ipISO][ipISO],
						NumLevRecomb[ipISO][ipISO+1],
						N_ISO_TE_RECOMB,
						i1 , i2 , i3, i4 );
					fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine.c_str() );
					fprintf( ioQQQ, 
						" iso_recomb_setup: please recompile the data file with the COMPile RECOmb COEFficient H-LIke [or HE-Like] command.\n" );
					cdEXIT(EXIT_FAILURE);
				}

				/* close the data file */
				fclose( ioDATA );
			}
		}
		/* We are compiling the he_iso_recomb.dat data file.	*/
		else if( iso_ctrl.lgCompileRecomb[ipISO] )
		{
			/* option to create table of recombination coefficients,
				* executed with the compile he-like command */
			FILE *ioRECOMB;

			ASSERT( !iso_ctrl.lgNoRecombInterp[ipISO] );

			/*RECOMBMAGIC the magic number for the table of recombination coefficients */
			/*NumLevRecomb[ipISO][nelem] the number of levels that will be done */
			/* create recombination coefficients  */
			ioRECOMB = open_data( chFilename[ipISO], "w" );
			fprintf(ioRECOMB,"%i\t%li\t%li\t%i\t%s isoelectronic sequence recomb data, created by COMPile RECOmb COEFficient H-LIke [or HE-Like] command, with %li %s levels, %li ion levels, and %i temperatures.\n",
				RECOMBMAGIC ,
				NumLevRecomb[ipISO][ipISO],
				NumLevRecomb[ipISO][ipISO+1],
				N_ISO_TE_RECOMB,
				iso_ctrl.chISO[ipISO],
				NumLevRecomb[ipISO][ipISO],
				elementnames.chElementSym[ipISO],
				NumLevRecomb[ipISO][ipISO+1],
				N_ISO_TE_RECOMB );

			for( nelem = ipISO; nelem < LIMELM; nelem++ )
			{
				/* this must pass since compile xx-like command reset numLevels to the macro */
				ASSERT( NumLevRecomb[ipISO][nelem] <= iso_sp[ipISO][nelem].numLevels_max );

				/* Zero out the recombination sum array.	*/
				for(i = 0; i < N_ISO_TE_RECOMB; i++)
				{
					TotalRecomb[ipISO][nelem][i] = 0.;
				}

				for( ipLo=ipHe1s1S; ipLo < NumLevRecomb[ipISO][nelem]; ipLo++ )
				{
					fprintf(ioRECOMB, "%li\t%li", nelem, ipLo );
					/* loop over temperatures to produce array of recombination coefficients	*/
					for(i = 0; i < N_ISO_TE_RECOMB; i++)
					{
						/* Store log of recombination coefficients, in N_ISO_TE_RECOMB half dec steps */
						RadRecombReturn = iso_radrecomb_from_cross_section( ipISO, exp10( TeRRCoef[i] ) ,nelem,ipLo);
						TotalRecomb[ipISO][nelem][i] += RadRecombReturn;
						RRCoef[ipISO][nelem][ipLo][i] = log10(RadRecombReturn);
						fprintf(ioRECOMB, "\t%f", RRCoef[ipISO][nelem][ipLo][i] );
					}
					fprintf(ioRECOMB, "\n" );
				}

				/* Store one additional line in XX_iso_recomb.dat that gives the total recombination,
				 * as computed by the sum so far, plus levels up to NHYDRO_MAX_LEVEL using Verner's fits,
				 * plus levels up to SumUpToThisN using Seaton 59, for each element and each temperature.	*/
				fprintf(ioRECOMB, "%li\t%li", nelem, NumLevRecomb[ipISO][nelem] );
				for(i = 0; i < N_ISO_TE_RECOMB; i++)
				{
					for( i1 = ( (nelem == ipISO) ? (RREC_MAXN + 1) : (LIKE_RREC_MAXN( nelem ) + 1) ); i1< NHYDRO_MAX_LEVEL; i1++ )
					{
						TotalRecomb[ipISO][nelem][i] += t_ADfA::Inst().H_rad_rec(nelem+1-ipISO,i1, exp10(TeRRCoef[i]));
					}
					for( i1 = NHYDRO_MAX_LEVEL; i1<=SumUpToThisN; i1++ )
					{
						TotalRecomb[ipISO][nelem][i] += Recomb_Seaton59( nelem+1-ipISO, exp10(TeRRCoef[i]), i1 );
					}
					fprintf(ioRECOMB, "\t%f", log10( TotalRecomb[ipISO][nelem][i] ) );
				}
				fprintf(ioRECOMB, "\n" );
			}
			/* end the file with the same information */
			fprintf(ioRECOMB,"%i\t%li\t%li\t%i\t%s isoelectronic sequence recomb data, created by COMPile RECOmb COEFficient [H-LIke/HE-Like] command, with %li %s levels, %li ion levels, and %i temperatures.\n",
				RECOMBMAGIC ,
				NumLevRecomb[ipISO][ipISO],
				NumLevRecomb[ipISO][ipISO+1],
				N_ISO_TE_RECOMB,
				iso_ctrl.chISO[ipISO],
				NumLevRecomb[ipISO][ipISO],
				elementnames.chElementSym[ipISO],
				NumLevRecomb[ipISO][ipISO+1],
				N_ISO_TE_RECOMB );

			fclose( ioRECOMB );

			fprintf( ioQQQ, "iso_recomb_setup: compilation complete, %s created.\n", chFilename[ipISO] );
			fprintf( ioQQQ, "The compilation is completed successfully.\n");
			cdEXIT(EXIT_SUCCESS);
		}
	}

	return;
}

double iso_dielec_recomb_rate( long ipISO, long nelem, long ipLo )
{
	double rate;
	long ipTe, i;
	double TeDRCoef[NUM_DR_TEMPS];
	const freeBound *fb = &iso_sp[ipISO][nelem].fb[ipLo];
	const double Te_over_Z1_Squared[NUM_DR_TEMPS] = {
		1.00000,	1.30103,	1.69897,	2.00000,	2.30103,	2.69897,	3.00000,
		3.30103,	3.69897,	4.00000,	4.30103,	4.69897,	5.00000,	5.30103,
		5.69897,	6.00000,	6.30103,	6.69897,	7.00000 };

	DEBUG_ENTRY( "iso_dielec_recomb_rate()" );

	/* currently only two iso sequences and only he-like is applicable. */
	ASSERT( ipISO == ipHE_LIKE );
	ASSERT( ipLo >= 0 );

	/* temperature grid is nelem^2 * constant temperature grid above. */
	for( i=0; i<NUM_DR_TEMPS; i++ )
	{
		TeDRCoef[i] = Te_over_Z1_Squared[i] + 2. * log10( (double) nelem );
	}

	if( ipLo == ipHe1s1S )
	{
		rate = 0.;
	}
	else if( ipLo<iso_sp[ipISO][nelem].numLevels_max )
	{
		if( phycon.alogte <= TeDRCoef[0] )
		{
			/* Take lowest tabulated value for low temperature end. */
			rate = fb->DielecRecombVsTemp[0];
		}
		else if( phycon.alogte >= TeDRCoef[NUM_DR_TEMPS-1] )
		{
			/* use T^-1.5 extrapolation at high temperatures. */
			rate = fb->DielecRecombVsTemp[NUM_DR_TEMPS-1] *
				exp10(  1.5* (TeDRCoef[NUM_DR_TEMPS-1] - phycon.alogte ) ) ;
		}
		else
		{
			/* find temperature in tabulated values.  */
			ipTe = hunt_bisect( TeDRCoef, NUM_DR_TEMPS, phycon.alogte );			

			ASSERT( (ipTe >=0) && (ipTe < NUM_DR_TEMPS-1)  );

			if( fb->DielecRecombVsTemp[ipTe+1] == 0. )
				rate = 0.;
			else if( fb->DielecRecombVsTemp[ipTe] == 0. )
				rate = fb->DielecRecombVsTemp[ipTe+1];
			else
			{
				/* interpolate between tabulated points */
				rate = log10( fb->DielecRecombVsTemp[ipTe]) +
					(phycon.alogte-TeDRCoef[ipTe])*
					(log10(fb->DielecRecombVsTemp[ipTe+1])-log10(fb->DielecRecombVsTemp[ipTe]))/
					(TeDRCoef[ipTe+1]-TeDRCoef[ipTe]);

				rate = exp10(  rate );
			}
		}
	}
	else 
	{
		rate = 0.;
	}

	ASSERT( rate >= 0. && rate < 1.0e-12 );

	return rate*iso_ctrl.lgDielRecom[ipISO];
}

/* TempInterp - interpolate on an array */
/** \todo	2	use a canned interpolation routine, no need for special one here */
STATIC double TempInterp( double* TempArray, double* ValueArray, long NumElements, double temp )
{
	static long int ipTe=-1;
	double rate = 0.;
	long i0;

	DEBUG_ENTRY( "TempInterp()" );

	double alogte = log10(temp);

	if( ipTe < 0 )
	{
		/* te totally unknown */
		if( ( alogte < TempArray[0] ) || 
			( alogte > TempArray[NumElements-1] ) )
		{
			fprintf(ioQQQ," TempInterp called with te out of bounds \n");
			cdEXIT(EXIT_FAILURE);
		}
		ipTe = hunt_bisect( TempArray, NumElements, alogte );			
	}
	else if( alogte < TempArray[ipTe] )
	{
		/* temp is too low, must also lower ipTe */
		ASSERT( alogte > TempArray[0] );
		/* decrement ipTe until it is correct */
		while( ( alogte < TempArray[ipTe] ) && ipTe > 0)
			--ipTe;
	}
	else if( alogte > TempArray[ipTe+1] )
	{
		/* temp is too high */
		ASSERT( alogte <= TempArray[NumElements-1] );
		/* increment ipTe until it is correct */
		while( ( alogte > TempArray[ipTe+1] ) && ipTe < NumElements-1)
			++ipTe;
	}

	ASSERT( (ipTe >=0) && (ipTe < NumElements-1) );

	/* ipTe should now be valid */
	ASSERT( ( alogte >= TempArray[ipTe] )
		&& ( alogte <= TempArray[ipTe+1] ) && ( ipTe < NumElements-1 ) );

	if( ValueArray[ipTe+1] == 0. && ValueArray[ipTe] == 0. )
	{
		rate = 0.;
	}
	else
	{
		/* Do a four-point interpolation */
		const int ORDER = 3; /* order of the fitting polynomial */
		i0 = max(min(ipTe-ORDER/2,NumElements-ORDER-1),0);
		rate = lagrange( &TempArray[i0], &ValueArray[i0], ORDER+1, alogte );
	}

	return rate;
}
