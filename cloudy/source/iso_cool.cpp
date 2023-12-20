/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_cool compute net cooling due to species in iso-sequences */
#include "cddefines.h"
#include "hydrogenic.h"
#include "elementnames.h"
#include "phycon.h"
#include "dense.h"
#include "thermal.h"
#include "cooling.h"
#include "iso.h"
#include "freebound.h"
#include "conv.h"
#include "rfield.h"
#include "opacity.h"
#include "vectorize.h"

STATIC void iso_rad_rec_cooling_discrete( const long ipISO, const long nelem);
STATIC double iso_rad_rec_cooling_approx( long ipISO, long nelem );
STATIC double iso_rad_rec_cooling_extra( long ipISO, long nelem, const double& ThinCoolingSum );

// set to true to enable debug print of contributors to collisional ionization cooling
const bool lgPrintIonizCooling = false;

void iso_cool(
	   /* iso sequence, 0 for hydrogenic */
	   long int ipISO , 
		/* nelem is charge -1, so 0 for H itself */
		long int nelem)
{
	long int ipbig; //, n;
	double biggest = 0.,
	  dCdT_all, 
	  edenIonAbund, 
	  CollisIonizatCoolingTotal, 
	  CollisIonizatHeatingTotal, 
	  dCollisIonizatCoolingTotalDT, 
	  HeatExcited,
	  heat_max,
	  CollisIonizatCooling, 
	  CollisIonizatHeating, 
	  CollisIonizatCoolingDT, 
	  hlone;

	valarray<double> CollisIonizatCoolingArr, 
	  CollisIonizatCoolingDTArr,
	  SavePhotoHeat,
	  SaveInducCool;

	long int nlo_heat_max  , nhi_heat_max;
	int coolnum = thermal.ncltot;

	/* place to put string labels for iso lines */
	char chLabel[NCOLNT_LAB_LEN+1];

	DEBUG_ENTRY( "iso_cool()" );

	t_iso_sp* sp = &iso_sp[ipISO][nelem];

	/* validate the incoming data */
	ASSERT( nelem >= ipISO );
	ASSERT( ipISO < NISO );
	ASSERT( nelem < LIMELM );
	/* local number of levels may be less than allocated number if continuum
	 * has been lowered due to high density */
	ASSERT( sp->numLevels_local <= sp->numLevels_max );

	if( dense.xIonDense[nelem][nelem-ipISO]<=0. ||
		!dense.lgElmtOn[nelem] )
	{
		/* all global variables must be zeroed here or set below */
		sp->coll_ion = 0.;
		sp->cLya_cool = 0.;
		sp->cLyrest_cool = 0.;
		sp->cBal_cool = 0.;
		sp->cRest_cool = 0.;
		sp->xLineTotCool = 0.;
		sp->RadRecCool = 0.;
		sp->FreeBnd_net_Cool_Rate = 0.;
		sp->dLTot = 0.;
		sp->RecomInducCool_Rate = 0.;
		// zero out line coolants
		for( long ipHi=1; ipHi < sp->numLevels_max; ++ipHi )
		{
			for( long ipLo=0; ipLo < ipHi; ++ipLo )
				CollisionZero( sp->trans(ipHi,ipLo).Coll() );
		}
		return;
	}

	/* create some space, these go to numLevels_local instead of _max
	 * since continuum may have been lowered by density */
	if( lgPrintIonizCooling )
	{
		CollisIonizatCoolingArr.resize( sp->numLevels_local );
		CollisIonizatCoolingDTArr.resize( sp->numLevels_local );
	}
	SavePhotoHeat.resize( sp->numLevels_local );
	SaveInducCool.resize( sp->numLevels_local );

	/***********************************************************************
	 *                                                                     *
	 * collisional ionization cooling, less three-body recombination  heat *
	 *                                                                     *
	 ***********************************************************************/

	/* will be net collisional ionization cooling, units erg/cm^3/s */
	CollisIonizatCoolingTotal = 0.;
	CollisIonizatHeatingTotal = 0.;
	dCollisIonizatCoolingTotalDT = 0.;

	// collisional ionization cooling and three body heating 
	// do not include the highest level since it has a fictitious collisional ionization
	// rate, to bring level into LTE with the continuum
	int lim = MIN2( sp->numLevels_max-1 , sp->numLevels_local );
	for( long n=0; n < lim; ++n )
	{
		/* total collisional ionization cooling less three body heating */
		CollisIonizatCooling =
			EN1RYD*sp->fb[n].xIsoLevNIonRyd*sp->fb[n].ColIoniz*dense.EdenHCorr*
			sp->st[n].Pop();
		CollisIonizatCoolingTotal += CollisIonizatCooling;
		CollisIonizatHeating =
			EN1RYD*sp->fb[n].xIsoLevNIonRyd*sp->fb[n].ColIoniz*dense.EdenHCorr*
			sp->fb[n].PopLTE* dense.eden*
			dense.xIonDense[nelem][nelem+1-ipISO];
		CollisIonizatHeatingTotal += CollisIonizatHeating;

		double CollisIonizatDiff = CollisIonizatCooling-CollisIonizatHeating;
		/* the derivative of the cooling */
		CollisIonizatCoolingDT = CollisIonizatDiff*
			(sp->fb[n].xIsoLevNIonRyd*TE1RYD/POW2(phycon.te)- thermal.halfte);


		dCollisIonizatCoolingTotalDT += CollisIonizatCoolingDT;
		// save values for debug printout
		if( lgPrintIonizCooling )
		{
			CollisIonizatCoolingArr[n] = CollisIonizatDiff;
			CollisIonizatCoolingDTArr[n] = CollisIonizatCoolingDT;
		}
	}

	double CollisIonizatNetCooling = 
		CollisIonizatCoolingTotal-CollisIonizatHeatingTotal;
	/* save net collisional ionization cooling less H-3 body heating
	 * for inclusion in printout */
	sp->coll_ion = CollisIonizatNetCooling;

	/* add this derivative to total */
	thermal.dCooldT += dCollisIonizatCoolingTotalDT;

	/* create a meaningful label */
	sprintf(chLabel , "IScion%2s%2s" , elementnames.chElementSym[ipISO] , 
		elementnames.chElementSym[nelem] );

	// Adding heating and cooling separately allows ionization solvers
	// to sense convergence close to LTE, but (at r5525) breaks thermal
	// equilibrium.
	if (0)
	{
		// New treatment -- separates cooling and heating as advertised
		/* dump the coolant onto the cooling stack */
		CoolAdd(chLabel,(realnum)nelem,CollisIonizatCoolingTotal);
		
		/* heating(0,3) is all three body heating, opposite of collisional 
		 * ionization cooling,
		 * would be unusual for this to be non-zero since would require excited
		 * state departure coefficients to be greater than unity */
		thermal.AddHeating(0,3, CollisIonizatHeatingTotal);
	}
	else
	{
		// Old treatment
		CoolAdd(chLabel,(realnum)nelem,MAX2(0.,CollisIonizatNetCooling));
		if (CollisIonizatNetCooling < 0)
			thermal.AddHeating(0,3,-CollisIonizatNetCooling);
	}

	/* debug block printing individual contributors to collisional ionization cooling */
	if( lgPrintIonizCooling && nelem==1 && ipISO==1 )
	{
		fprintf(ioQQQ,"DEBUG coll ioniz cool contributors:");
		for( long n=0; n < sp->numLevels_local; n++ )
		{
			if( CollisIonizatCoolingArr[n] / SDIV( CollisIonizatNetCooling ) > 0.1 )
				fprintf(ioQQQ," %2li %.1e",
					n,
					CollisIonizatCoolingArr[n]/ SDIV( CollisIonizatNetCooling ) );
		}
		fprintf(ioQQQ,"\n");
		fprintf(ioQQQ,"DEBUG coll ioniz derivcontributors:");
		for( long n=0; n < sp->numLevels_local; n++ )
		{
			if( CollisIonizatCoolingDTArr[n] / SDIV( dCollisIonizatCoolingTotalDT ) > 0.1 )
				fprintf(ioQQQ," %2li %.1e",
					n,
					CollisIonizatCoolingDTArr[n]/ SDIV( dCollisIonizatCoolingTotalDT ) );
		}
		fprintf(ioQQQ,"\n");
	}

	/***********************************************************************
	 *                                                                     *
	 * hydrogen recombination free-bound free bound cooling                *
	 *                                                                     *
	 ***********************************************************************/

	/* this is the product of the ion abundance times the electron density */
	edenIonAbund = dense.eden*dense.xIonDense[nelem][nelem+1-ipISO];

	sp->RadRecCool = 0.;
	if( dense.gas_phase[nelem] > 1e-3 * dense.xNucleiTotal )
	{
		iso_rad_rec_cooling_discrete( ipISO, nelem );
	}
	else
	{
		double ThinCoolingSum = iso_rad_rec_cooling_approx( ipISO, nelem );
		/* this is a cooling "topoff" term - significant for approach to STE */
		sp->RadRecCool = iso_rad_rec_cooling_extra( ipISO, nelem, ThinCoolingSum );
	}
	
	/* this is now total free-bound cooling */
	for( long n=0; n < sp->numLevels_local; n++ )
	{
		sp->RadRecCool += sp->fb[n].RadRecCoolCoef;
	}

	// now multiply by density squared to get real cooling, erg cm-3 s-1
	sp->RadRecCool *= edenIonAbund;

	{
		/* debug block for state specific recombination cooling */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC  )
		{
			if( nelem==ipISO )
			{
				for( long n=0; n < sp->numLevels_local; n++ )
				{
					fprintf(ioQQQ,"\t%.2f",sp->fb[n].RadRecCoolCoef/sp->RadRecCool);
				}
				fprintf(ioQQQ,"\n");
			}
		}
	}


	/***********************************************************************
	 *                                                                     
	 * heating  by photoionization of                                      
	 * excited states of all species                            
	 *                                                                     
	 ***********************************************************************/

	/* photoionization of excited levels */
	HeatExcited = 0.;
	ipbig = -1000;
	for( long n=1; n < sp->numLevels_local; ++n )
	{
		ASSERT( sp->fb[n].PhotoHeat >= 0. );
		ASSERT( sp->st[n].Pop() >= 0. );

		SavePhotoHeat[n] = sp->st[n].Pop()*
			sp->fb[n].PhotoHeat;
		HeatExcited += SavePhotoHeat[n];
		if( SavePhotoHeat[n] > biggest )
		{
			biggest = SavePhotoHeat[n];
			ipbig = (int)n;
		}
	}
	{
		/* debug block for heating due to photo of each n */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC  && ipISO==0 && nelem==0  && nzone > 700)
		{
			/* this was not done above */
			SavePhotoHeat[ipH1s] = sp->st[ipH1s].Pop()*
				sp->fb[ipH1s].PhotoHeat;
			fprintf(ioQQQ,"ipISO = %li nelem=%li ipbig=%li biggest=%g HeatExcited=%.2e ctot=%.2e\n",
				ipISO , nelem,
				ipbig , 
				biggest,
				HeatExcited ,
				thermal.ctot);
			/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
			for( long n=ipH1s; n< sp->numLevels_local; ++n )
			{
				fprintf(ioQQQ,"DEBUG phot heat%2li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n", 
					n,
					SavePhotoHeat[n]/HeatExcited,
					dense.xIonDense[nelem][nelem+1-ipISO],
					sp->st[n].Pop(),
					sp->fb[n].PhotoHeat,
					sp->fb[n].gamnc);
			}
		}
	}

	/* FreeBnd_net_Cool_Rate is net cooling due to recombination 
	 * RadRecCool is total radiative recombination cooling sum to all levels,
	 * with n>=2 photoionization heating subtracted */
	sp->FreeBnd_net_Cool_Rate = sp->RadRecCool - HeatExcited;
	/*fprintf(ioQQQ,"DEBUG Hn2\t%.3e\t%.3e\n",
		-sp->RadRecCool/SDIV(thermal.htot),
		HeatExcited/SDIV(thermal.htot));*/

	/* heating(0,1) is all excited state photoionization heating from ALL 
	 * species, this is set to zero in CoolEvaluate before loop where this 
	 * is called. */
	thermal.AddHeating(0,1, MAX2(0.,-sp->FreeBnd_net_Cool_Rate));

	/* net free-bound cooling minus free-free heating */
	/* create a meaningful label */
	sprintf(chLabel , "ISrcol%2s%2s" , elementnames.chElementSym[ipISO]  , 
		elementnames.chElementSym[nelem]);
	CoolAdd(chLabel, (realnum)nelem, MAX2(0.,sp->FreeBnd_net_Cool_Rate));

	/* if rec coef goes as T^-.8, but times T, so propto t^.2 */
	thermal.dCooldT += 0.2*sp->FreeBnd_net_Cool_Rate*phycon.teinv;

	/***********************************************************************
	 *                                                                     *
	 * induced recombination cooling                                       *
	 *                                                                     *
	 ***********************************************************************/

	sp->RecomInducCool_Rate = 0.;
	/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
	for( long n=0; n < sp->numLevels_local; ++n )
	{
		/* >>chng 02 jan 22, removed cinduc, replace with RecomInducCool */
		SaveInducCool[n] = sp->fb[n].RecomInducCool_Coef*sp->fb[n].PopLTE*edenIonAbund;
		sp->RecomInducCool_Rate += SaveInducCool[n];
		thermal.dCooldT += SaveInducCool[n]*
			(sp->fb[n].xIsoLevNIonRyd/phycon.te_ryd - 1.5)*phycon.teinv;
	}

	{
		/* print rec cool, induced rec cool, photo heat */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && ipISO==0 && nelem==5 /**/ )
		{
			fprintf(ioQQQ," ipISO=%li nelem=%li ctot = %.2e\n",
				ipISO,
				nelem,
				thermal.ctot);
			fprintf(ioQQQ,"sum\t%.2e\t%.2e\t%.2e\n", 
				HeatExcited,
				sp->RadRecCool,
				sp->RecomInducCool_Rate);
			fprintf(ioQQQ,"sum\tp ht\tr cl\ti cl\n");

			/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
			for( long n=0; n< sp->numLevels_local; ++n )
			{
				fprintf(ioQQQ,"%li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e  \n", 
					n,
					SavePhotoHeat[n],
					sp->fb[n].RadRecCoolCoef,
					SaveInducCool[n] ,
					sp->fb[n].RecomInducCool_Coef,
					sp->fb[n].PopLTE,
					sp->fb[n].RecomInducRate);
			}
			fprintf(ioQQQ," \n");
		}
	}
	/* create a meaningful label - induced rec cooling */
	sprintf(chLabel , "ISicol%2s%2s" , elementnames.chElementSym[ipISO]  , 
		elementnames.chElementSym[nelem]);
	/* induced rec cooling */
	CoolAdd(chLabel,(realnum)nelem,sp->RecomInducCool_Rate);

	/* find total collisional energy exchange due to bound-bound */
	sp->xLineTotCool = 0.;
	dCdT_all = 0.;
	heat_max = 0.;
	nlo_heat_max = -1;
	nhi_heat_max = -1;

	if (0)
	{
		long ipHi = 0;
		fprintf(ioQQQ,"%ld %ld %ld %g\n",nelem,ipISO,ipHi,
				  sp->st[ipHi].Pop()/sp->st[ipHi].g());
		for( ipHi=1; ipHi < sp->numLevels_local; ++ipHi )
		{
			fprintf(ioQQQ,"%ld %ld %ld %g\n",nelem,ipISO,ipHi,
					  sp->st[ipHi].Pop()/(sp->st[ipHi].Boltzmann()*
												 sp->st[ipHi].g()));
		}
		
	}

	/* loop does not include highest levels - their population may be
	 * affected by topoff */
	vector<double> Boltzmann(sp->numLevels_local-1, 1.0);
	avx_ptr<double> arg(1,sp->numLevels_local), bstep(1,sp->numLevels_local);
	for( long ipHi=1; ipHi < sp->numLevels_local; ++ipHi )
	{
		fixit("Wouldn't need to mask this out if levels were in order");
		arg[ipHi] = -max((sp->st[ipHi].energy().K()-sp->st[ipHi-1].energy().K()) / phycon.te, -38.);
	}
	vexp( arg.ptr0(), bstep.ptr0(), 1, sp->numLevels_local );

	ColliderDensities colld(colliders);

	for( long ipHi=1; ipHi < sp->numLevels_local; ++ipHi )
	{
		for (long ipLo = 0; ipLo < ipHi; ++ipLo )
			Boltzmann[ipLo] *= bstep[ipHi];

		for( long ipLo=0; ipLo < ipHi; ++ipLo )
		{
			/* collisional energy exchange between ipHi and ipLo - net cool */
			hlone = 
			  sp->trans(ipHi,ipLo).Coll().ColUL( colld ) *
			  (sp->st[ipLo].Pop()*
			  Boltzmann[ipLo]*
			  sp->st[ipHi].g()/sp->st[ipLo].g() - 
			  sp->st[ipHi].Pop())*
			  sp->trans(ipHi,ipLo).EnergyErg();

			if( hlone > 0. )
				sp->trans(ipHi,ipLo).Coll().cool() = hlone;
			else
				sp->trans(ipHi,ipLo).Coll().heat() = -1.*hlone;

			sp->xLineTotCool += hlone;

			/* next get derivative */
			if( hlone > 0. )
			{
				/* usual case, this line was a net coolant - for derivative 
				 * taking the exponential from ground gave better
				 * representation of effects */
				dCdT_all += hlone*
					(sp->trans(ipHi,ipH1s).EnergyK()*thermal.tsq1 - thermal.halfte);
			}
			else
			{
				/* this line heats the gas, remember which one it was,
				 * the following three vars never are used, but could be for
				 * debugging */
				if( hlone < heat_max )
				{
					heat_max = hlone;
					nlo_heat_max = ipLo;
					nhi_heat_max = ipHi;
				} 
				dCdT_all -= hlone*thermal.halfte;
			}
		}
	}
	{
		/* debug block announcing which line was most important */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			if( nelem==ipISO )
				fprintf(ioQQQ,"%li %li %.2e\n", nlo_heat_max, nhi_heat_max, heat_max );
		}
	}

	/* Lyman line collisional heating/cooling */
	/* Lya itself */
	sp->cLya_cool = 0.;
	/* lines higher than Lya */
	sp->cLyrest_cool = 0.;

	for( long ipHi=1; ipHi < sp->numLevels_local; ipHi++ )
	{
		hlone = sp->trans(ipHi,ipH1s).Coll().ColUL( colld ) *
		  (sp->st[0].Pop()*
		  sp->st[ipHi].Boltzmann()*
		  sp->st[ipHi].g()/sp->st[0].g() - 
		  sp->st[ipHi].Pop())* 
		  sp->trans(ipHi,0).EnergyErg();

		if( ipHi == iso_ctrl.nLyaLevel[ipISO] )
		{
			sp->cLya_cool = hlone;

		}
		else
		{
			/* sum energy in higher lyman lines */
			sp->cLyrest_cool += hlone;
		}
	}

	/* Balmer line collisional heating/cooling and derivative 
	 * only used for print out, incorrect if not H so don't calculate it */
	sp->cBal_cool = 0.;
	if (nelem == ipHYDROGEN)
	{
		double boltzH2s = bstep[2];
		double boltzH2p = 1.;
		for( long ipHi=3; ipHi < sp->numLevels_local; ipHi++ )
		{
			boltzH2s *= bstep[ipHi];
			boltzH2p *= bstep[ipHi];
			
			/* single line cooling */
			long ipLo = ipH2s;
			hlone = 
				sp->trans(ipHi,ipLo).Coll().ColUL( colld ) *
				(sp->st[ipLo].Pop()*
				 boltzH2s*
				 sp->st[ipHi].g()/sp->st[ipLo].g() - 
				 sp->st[ipHi].Pop())*
				sp->trans(ipHi,ipLo).EnergyErg();
			ASSERT( sp->st[ipHi].energy().Erg() >
					  sp->st[ipLo].energy().Erg()	);
			
			ipLo = ipH2p;
			hlone += 
				sp->trans(ipHi,ipLo).Coll().ColUL( colld ) *
				(sp->st[ipLo].Pop()*
				 boltzH2p*
				 sp->st[ipHi].g()/sp->st[ipLo].g() - 
				 sp->st[ipHi].Pop())*
				sp->trans(ipHi,ipLo).EnergyErg();
			ASSERT( sp->st[ipHi].energy().Erg() >
					  sp->st[ipLo].energy().Erg()	);
			
			/* this is only used to add to line array */
			sp->cBal_cool += hlone;
		}
	}

	/* all hydrogen lines from Paschen on up, but not Balmer 
	 * only used for printout, incorrect if not H so don't calculate it */
	sp->cRest_cool = 0.;
	if (nelem == ipHYDROGEN)
	{
		for (unsigned lev = 0; lev < Boltzmann.size(); ++lev)
			Boltzmann[lev] = 1.;
		
		for( long ipHi=4; ipHi < sp->numLevels_local; ++ipHi )
		{
			for ( long ipLo = 0; ipLo < ipHi; ++ipLo )
				Boltzmann[ipLo] *= bstep[ipHi];

			for( long ipLo=3; ipLo < ipHi; ++ipLo )
			{
				hlone = 
					sp->trans(ipHi,ipLo).Coll().ColUL( colld ) *
					(sp->st[ipLo].Pop()*
					 Boltzmann[ipLo]*
					 sp->st[ipHi].g()/sp->st[ipLo].g() - 
					 sp->st[ipHi].Pop())*
					sp->trans(ipHi,ipLo).EnergyErg();
				
				sp->cRest_cool += hlone;
			}
		}
	}

	/* add total line heating or cooling into stacks, derivatives */
	/* line energy exchange can be either heating or coolant
	 * must add this to total heating or cooling, as appropriate */
	/* create a meaningful label */
	sprintf(chLabel , "ISclin%2s%2s" , elementnames.chElementSym[ipISO] , 
		elementnames.chElementSym[nelem]);
	if( sp->xLineTotCool > 0. )
	{
		/* species is a net coolant label gives iso sequence, "wavelength" gives element */
		CoolAdd(chLabel,(realnum)nelem,sp->xLineTotCool);
		thermal.dCooldT += dCdT_all;
		sp->dLTot = 0.;
	}
	else
	{
		/* species is a net heat source, thermal.heating(0,23)was set to 0 in CoolEvaluate*/
		thermal.AddHeating(0,23,-sp->xLineTotCool);
		CoolAdd(chLabel,(realnum)nelem,0.);
		sp->dLTot = -dCdT_all;
	}

	{
		/* debug print for understanding contributors to HLineTotCool */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			if( nelem == 0 )
			{
				fprintf(ioQQQ,"%.2e la %.2f restly %.2f barest %.2f hrest %.2f\n",
					sp->xLineTotCool ,
					sp->cLya_cool/sp->xLineTotCool ,
					sp->cLyrest_cool/sp->xLineTotCool ,
					sp->cBal_cool/sp->xLineTotCool ,
					sp->cRest_cool/sp->xLineTotCool );
			}
		}
	}
	{
		/* this is an option to print out each rate affecting each level in strict ste,
		 * this is a check that the rates are indeed in detailed balance */
		enum {DEBUG_LOC=false};
		enum {LTEPOP=true};
		/* special debug print for gas near STE */
		if( DEBUG_LOC  && (nelem==1 || nelem==0) )
		{
			/* LTEPOP is option to assume STE for rates */
			if( LTEPOP )
			{
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for( long n=ipH1s; n<sp->numLevels_local; ++n )
				{
					fprintf(ioQQQ,"%li\t%li\t%g\t%g\t%g\t%g\tT=\t%g\t%g\t%g\t%g\n", nelem,n,
						sp->fb[n].gamnc *sp->fb[n].PopLTE, 
						/* induced recom, intergral times hlte */
						/*sp->fb[n].RadRecomb[ipRecRad]+sp->rinduc[n]  ,*/
						/* >>chng 02 jan 22, remove rinduc, replace with RecomInducRate */
						sp->fb[n].RadRecomb[ipRecRad]+ 
							sp->fb[n].RecomInducRate*sp->fb[n].PopLTE  ,
						/* induced rec */
						sp->fb[n].RecomInducRate*sp->fb[n].PopLTE  ,
						/* spontaneous recombination */
						sp->fb[n].RadRecomb[ipRecRad] ,
						/* heating, followed by two processes that must balance it */
						sp->fb[n].PhotoHeat*sp->fb[n].PopLTE, 
						sp->fb[n].RecomInducCool_Coef*sp->fb[n].PopLTE+sp->fb[n].RadRecCoolCoef,
						/* induced rec cooling, integral times hlte */
						sp->fb[n].RecomInducCool_Coef*sp->fb[n].PopLTE ,
						/* free-bound cooling per unit vol */
						sp->fb[n].RadRecCoolCoef );
				}
			}
			else
			{
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for( long n=ipH1s; n<sp->numLevels_local; ++n )
				{
					fprintf(ioQQQ,"%li\t%li\t%g\t%g\t%g\t%g\tT=\t%g\t%g\t%g\t%g\n", nelem,n,
						sp->fb[n].gamnc*sp->st[n].Pop(), 
						/* induced recom, intergral times hlte */
						sp->fb[n].RadRecomb[ipRecRad]*edenIonAbund+
							sp->fb[n].RecomInducRate*sp->fb[n].PopLTE *edenIonAbund ,
						sp->fb[n].RadRecomb[ipRecRad]*edenIonAbund ,
						sp->fb[n].RecomInducRate*sp->fb[n].PopLTE *edenIonAbund ,
						/* heating, followed by two processes that must balance it */
						SavePhotoHeat[n], 
						SaveInducCool[n]+sp->fb[n].RadRecCoolCoef*edenIonAbund ,
						/* induced rec cooling, integral times hlte */
						SaveInducCool[n] ,
						/* free-bound cooling per unit vol */
						sp->fb[n].RadRecCoolCoef*edenIonAbund );
				}
			}
		}
	}

	/* Add Coolant together */
	for( int i = coolnum; i < thermal.ncltot ; i++ )
		thermal.elementcool[nelem] += thermal.cooling[i];
	
	return;
}

STATIC double iso_rad_rec_cooling_approx( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_rad_rec_cooling_approx()" );
	
	t_iso_sp* sp = &iso_sp[ipISO][nelem];

	/* now do case b sum to compare with exact value below */
	double ThinCoolingSum = 0.;

	/* radiative recombination cooling for all excited states */
	for( long n=0; n < sp->numLevels_local; n++ )
	{
		double thin = 0.;

		if( n==0 && ipISO==ipH_LIKE )
		{
			/* do ground with special approximate fits to Ferland et al. '92 */
			thin = HydroRecCool(
				/* n is the prin quantum number on the physical scale */
				1 , 
				/* nelem is the charge on the C scale, 0 is hydrogen */
				nelem);
		}
		else
		{
			/* this is the cooling before correction for optical depths */
			 thin = sp->fb[n].RadRecomb[ipRecRad]*
				/* arg is the scaled temperature, T * n^2 / Z^2, 
				 * n is principal quantum number, Z is charge, 1 for H */
				HCoolRatio( 
				phycon.te * POW2( (double)sp->st[n].n() / (double)(nelem+1-ipISO) ))*
				/* convert results to energy per unit vol */
				phycon.te * BOLTZMANN;
		}

		/* the cooling, corrected for optical depth */
		sp->fb[n].RadRecCoolCoef = sp->fb[n].RadRecomb[ipRecNetEsc] * thin;

		/* keep track of case b sum for topoff below */
		if( n > 0 )
			ThinCoolingSum += thin;
	}

	return ThinCoolingSum;
}

STATIC double iso_rad_rec_cooling_extra( long ipISO, long nelem, const double& ThinCoolingSum )
{
	DEBUG_ENTRY( "iso_rad_rec_cooling_extra()" );
	
	t_iso_sp* sp = &iso_sp[ipISO][nelem];
	
	double RecCoolExtra = 0.;

	/* Case b sum of optically thin radiative recombination cooling. 
	 * add any remainder to the sum from above - high precision is needed 
	 * to get STE result to converge close to equilibrium - only done for
	 * H-like ions where exact result is known */
	if( ipISO == ipH_LIKE )
	{
		double ThinCoolingCaseB;
		/* these expressions are only valid for hydrogenic sequence */
		if( nelem == 0 )
		{
			/*expansion for hydrogen itself */
			ThinCoolingCaseB = (-25.859117 + 
			0.16229407*phycon.telogn[0] + 
			0.34912863*phycon.telogn[1] - 
			0.10615964*phycon.telogn[2])/(1. + 
			0.050866793*phycon.telogn[0] - 
			0.014118924*phycon.telogn[1] + 
			0.0044980897*phycon.telogn[2] + 
			6.0969594e-5*phycon.telogn[3]);
		}
		else
		{
			/* same expansion but for hydrogen ions */
			ThinCoolingCaseB = (-25.859117 + 
			0.16229407*(phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]) + 
			0.34912863*POW2(phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]) - 
			0.10615964*powi( (phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]),3) )/(1. + 
			0.050866793*(phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]) - 
			0.014118924*POW2(phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]) + 
			0.0044980897*powi( (phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]),3) + 
			6.0969594e-5*powi( (phycon.telogn[0]-phycon.sqlogz[nelem-ipISO]),4) );
		}

		/* now convert to linear cooling coefficient */
		ThinCoolingCaseB = POW3(1.+nelem-ipISO)*exp10(ThinCoolingCaseB)/(phycon.te/POW2(1.+nelem-ipISO) );

		/* this is the error, expect positive since do not include infinite number of levels */
		RecCoolExtra = ThinCoolingCaseB - ThinCoolingSum;
	}
	else
	{
		RecCoolExtra = 0.;
	}

	/* don't let the extra be negative - should be positive if H-like, negative for
	 * he-like only due to real difference in recombination coefficients */
	RecCoolExtra = MAX2(0., RecCoolExtra );

	return RecCoolExtra * sp->fb[sp->numLevels_local-1].RadRecomb[ipRecNetEsc];
}

STATIC void iso_rad_rec_cooling_discrete( const long ipISO, const long nelem )
{
	DEBUG_ENTRY( "iso_rad_rec_cooling_discrete()" );
	
	if( fp_equal( phycon.te, iso_ctrl.RRC_TeUsed[ipISO][nelem] ) && conv.nTotalIoniz )
		return;

	iso_ctrl.RRC_TeUsed[ipISO][nelem] = phycon.te;

	ASSERT( nelem >= ipISO );
	ASSERT( nelem < LIMELM );

	// recombination continua for all iso seq - 
	// if this stage of ionization exists 
	if( dense.IonHigh[nelem] >= nelem+1-ipISO  )
	{
		t_iso_sp* sp = &iso_sp[ipISO][nelem];

		// loop over all levels to include recombination diffuse continua,
		// pick highest energy continuum point that opacities extend to 
		long ipHi = rfield.nflux;
		// >>chng 06 aug 17, should go to numLevels_local instead of _max.
		avx_ptr<double> arg(sp->numLevels_local), val(sp->numLevels_local);
		for( long n=0; n < sp->numLevels_local; n++ )
		{
			long ipLo = sp->fb[n].ipIsoLevNIonCon-1;
			double thresh = sp->fb[n].xIsoLevNIonRyd;
			double widflx = rfield.anumax(ipLo) - thresh;
			arg[n] = -widflx/phycon.te_ryd;
		}
		vexpm1( arg.ptr0(), val.ptr0(), 0, sp->numLevels_local );
		for( long n=0; n < sp->numLevels_local; n++ )
		{
			// the number is (2 pi me k/h^2) ^ -3/2 * 8 pi/c^2 / ge - it includes
			// the stat weight of the free electron in the demominator 
			double gamma = 0.5*MILNE_CONST*sp->st[n].g()/iso_ctrl.stat_ion[ipISO]/phycon.te/phycon.sqrte;

			// loop over all recombination continua 
			// escaping part of recombinations are added to rfield.ConEmitLocal 
			// added to ConInterOut at end of routine
			long ipLo = sp->fb[n].ipIsoLevNIonCon-1;
			long offset = -sp->fb[n].ipIsoLevNIonCon + sp->fb[n].ipOpac;
			double thresh = sp->fb[n].xIsoLevNIonRyd;
			ASSERT( rfield.anumin(ipLo) <= thresh && thresh < rfield.anumax(ipLo) );
			// the first zone is special since it contains the threshold, we treat this outside the loop
			long nu = ipLo;
			double widflx = rfield.anumax(ipLo) - thresh;
			// fac is the fraction of the cell width that is above threshold
			double fac = widflx/rfield.widflx(ipLo);
			double bfac = -fac*phycon.te_ryd*val[n]/widflx;
			double efac = exp(-widflx/phycon.te_ryd);
			double energyAboveThresh = widflx/2.;
			double photon = bfac*rfield.widflx(nu)*opac.OpacStack[nu+offset]*rfield.anu2(nu);
			double RadRecCoolCoef = energyAboveThresh*photon;
			for( nu=ipLo+1; nu < ipHi; nu++ )
			{
				// efac = exp(-(rfield.anumin(nu)-thresh)/phycon.te_ryd)
				bfac = efac*rfield.ContBoltzHelp1[nu];
				efac *= rfield.ContBoltzHelp2[nu];
				energyAboveThresh = rfield.anu(nu) - thresh;
				/* gamma*photon is in photons cm^3 s^-1 per cell */
				/* multiplication with gamma is done outside this loop */
				photon = bfac*rfield.widflx(nu)*opac.OpacStack[nu+offset]*rfield.anu2(nu);
				/* multiplication with sp->fb[n].RadRecomb[ipRecNetEsc] is done outside this loop */
				double one = energyAboveThresh*photon;
				RadRecCoolCoef += one;
				// no point in wasting time on adding tiny numbers
				if( one < 1.e-20*RadRecCoolCoef )
					break;
			}

			// convert to erg cm-3 s-1
			sp->fb[n].RadRecCoolCoef = gamma*RadRecCoolCoef*sp->fb[n].RadRecomb[ipRecNetEsc]*EN1RYD;
		}

		// no RRC emission from levels that do not exist
		for( long n=sp->numLevels_local;n<sp->numLevels_max; n++ )
		{
			sp->fb[n].RadRecCoolCoef = 0.;
		}
	}

	return;
}
