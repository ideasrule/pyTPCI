/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "atmdat.h"
#include "phycon.h"
#include "taulines.h"
#include "atoms.h"
#include "rfield.h"
#include "conv.h"
#include "secondaries.h"
#include "thermal.h"
#include "cooling.h"
#include "ionbal.h"
#include "iso.h"
#include "mole.h"
#include "dense.h"
#include "lines_service.h"
#include "trace.h"
#include "doppvel.h"
#include "oxy.h"
#include "hydrogenic.h"
#include "continuum.h"
#include "vectorize.h"
#include "container_classes.h"

STATIC double LeidenCollRate(long, long, const TransitionProxy& ,double);
STATIC double StoutCollRate(long ipSpecies, long ipCollider, const TransitionProxy&, double ftemp);
STATIC double ChiantiCollRate(long ipSpecies, long ipCollider, const TransitionProxy&, double ftemp);
STATIC void setXtraRatesO1(const TransitionProxy& tr, double &xtraExRate, double &xtraDxRate);
STATIC void setXtraRatesCa2(const TransitionProxy& tr, double &xtraDxRate);
STATIC void setXtraRatesFe2(const TransitionProxy& tr, double &xtraExRate, double &xtraDxRate);

static const bool DEBUGSTATE = false;

static realnum dBaseAbund(long ipSpecies)
{
	realnum abund;
	/* first find current density (cm-3) of species */
	if( dBaseSpecies[ipSpecies].lgMolecular )
	{
		/** \todo	0	this pointer should be cached one time, and the species
		 * removed from the list if it is not computed */
		molezone *SpeciesCurrent = findspecieslocal(dBaseSpecies[ipSpecies].chLabel);
		if( !exists( SpeciesCurrent ) )
		{
			/* did not find the species - print warning for now */
			if( !conv.nTotalIoniz )
				fprintf(ioQQQ," PROBLEM dBase_solve did not find molecular species %li\n",ipSpecies);
		}
		abund = (realnum)SpeciesCurrent->den;
	}
	else
	{
		/* an atom or ion */
		ASSERT( dBaseStates[ipSpecies][0].nelem()<=LIMELM && dBaseStates[ipSpecies][0].IonStg()<=dBaseStates[ipSpecies][0].nelem()+1 );
		abund = dense.xIonDense[ dBaseStates[ipSpecies][0].nelem()-1 ][ dBaseStates[ipSpecies][0].IonStg()-1 ];
	}
	
	abund *= dBaseSpecies[ipSpecies].fracType * dBaseSpecies[ipSpecies].fracIsotopologue;
	return abund;
}

void dBaseTrim(void)
{
	DEBUG_ENTRY( "dBaseTrim()" );
	// initialization at start of each iteration
	if( conv.nTotalIoniz == 0)
	{
		for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
		{
			dBaseSpecies[ipSpecies].lgActive = true;
		}
	}

	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		realnum abund = dBaseAbund(ipSpecies);
		bool lgMakeInActive = (abund <= 1e-20 * dense.xNucleiTotal);
		if( lgMakeInActive && dBaseSpecies[ipSpecies].lgActive )
		{
			// zero out populations and intensities, if previously not set
			dBaseStates[ipSpecies][0].Pop() = 0.;
			for(long ipHi = 1; ipHi < dBaseSpecies[ipSpecies].numLevels_max; ipHi++ )
			{	
				dBaseStates[ipSpecies][ipHi].Pop() = 0.;
			}
			for (TransitionList::iterator tr=dBaseTrans[ipSpecies].begin(); 
				  tr != dBaseTrans[ipSpecies].end(); ++tr)
			{
				(*tr).Emis().xIntensity() = 0.;
				(*tr).Emis().xObsIntensity() = 0.;
				(*tr).Coll().col_str() = 0.;
				(*tr).Coll().cool() = 0.;
				(*tr).Coll().heat() = 0.;
			}
			dBaseSpecies[ipSpecies].lgActive = false;
		}

		if( !lgMakeInActive )
			dBaseSpecies[ipSpecies].lgActive = true;
	}
}

void dBaseUpdateCollCoeffs(void)
{
	DEBUG_ENTRY( "dBaseUpdateCollCoeffs()" );
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		if( !dBaseSpecies[ipSpecies].lgActive )
			continue;
		
		/*Setting all the collision strengths and collision rate to zero*/
		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			for( long k=0; k<ipNCOLLIDER; ++k )
				(*tr).Coll().rate_coef_ul_set()[k] = 0.f;
		}
		/* update the collision rates */
		/* molecule */
		if( dBaseSpecies[ipSpecies].database == "LAMDA" )
		{
			for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
				 tr != dBaseTrans[ipSpecies].end(); ++tr)
			{
				int ipHi = (*tr).ipHi();
				if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
					continue;
				for( long intCollNo=0; intCollNo<ipNCOLLIDER; intCollNo++)
				{
					/*using the collision rate coefficients directly*/
					(*tr).Coll().rate_coef_ul_set()[intCollNo] =
						LeidenCollRate(ipSpecies, intCollNo, *tr, phycon.te);
				}
				tr->Coll().is_gbar() = 0;
			}
		}
		/* Chianti */
		else if( dBaseSpecies[ipSpecies].database == "Chianti" )
		{
			for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
				 tr != dBaseTrans[ipSpecies].end(); ++tr)
			{
				if ((*tr).ipHi() >= dBaseSpecies[ipSpecies].numLevels_local)
					continue;
				for( long intCollNo=0; intCollNo<ipNCOLLIDER; intCollNo++)
				{
					(*tr).Coll().rate_coef_ul_set()[intCollNo] =
						ChiantiCollRate(ipSpecies, intCollNo, *tr, phycon.te);
				}
				tr->Coll().is_gbar() = 0;
			}
		}
		/* Stout */
		else if( dBaseSpecies[ipSpecies].database == "Stout" )
		{
			for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
				 tr != dBaseTrans[ipSpecies].end(); ++tr)
			{
				if ((*tr).ipHi() >= dBaseSpecies[ipSpecies].numLevels_local)
					continue;
				for( long intCollNo=0; intCollNo<ipNCOLLIDER; intCollNo++)
				{
					(*tr).Coll().rate_coef_ul_set()[intCollNo] =
							StoutCollRate(ipSpecies, intCollNo, *tr, phycon.te);
				}
				tr->Coll().is_gbar() = 0;
			}
		}
		else
			TotalInsanity();

		/* guess some missing data */
		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			const CollisionProxy &coll_temp = (*tr).Coll();

				/* make educated guesses for some missing data */
			if( dBaseSpecies[ipSpecies].lgMolecular )
			{
				/*The collision rate coefficients for helium should not be present and that for molecular hydrogen should be present*/
				if( AtmolCollRateCoeff[ipSpecies][ipATOM_HE].temps.size() == 0 &&
					AtmolCollRateCoeff[ipSpecies][ipH2].temps.size() != 0 )
				{
					coll_temp.rate_coef_ul_set()[ipATOM_HE] = 0.7f * coll_temp.rate_coef_ul()[ipH2];
				}

				/* Put in something for hydrogen collisions if not in database */
				if( AtmolCollRateCoeff[ipSpecies][ipATOM_H].temps.size() == 0 )
				{
					if( AtmolCollRateCoeff[ipSpecies][ipATOM_HE].temps.size() != 0 ) //He0
					{
						coll_temp.rate_coef_ul_set()[ipATOM_H] = 2.0f * coll_temp.rate_coef_ul()[ipATOM_HE];
					}
					else if( AtmolCollRateCoeff[ipSpecies][ipH2_ORTHO].temps.size() != 0 ) //ortho-H2
					{
						coll_temp.rate_coef_ul_set()[ipATOM_H] = 1.4f * coll_temp.rate_coef_ul()[ipH2_ORTHO];
					}
					else if( AtmolCollRateCoeff[ipSpecies][ipH2_PARA].temps.size() != 0 ) //para-H2
					{
						coll_temp.rate_coef_ul_set()[ipATOM_H] = 1.4f * coll_temp.rate_coef_ul()[ipH2_PARA];
					}
					else if( AtmolCollRateCoeff[ipSpecies][ipH2].temps.size() != 0 ) // total H2
					{
						coll_temp.rate_coef_ul_set()[ipATOM_H] = 1.4f * coll_temp.rate_coef_ul()[ipH2];
					}
					else
						coll_temp.rate_coef_ul_set()[ipATOM_H] = 1e-13f * (*(*tr).Lo()).g();
				}

				/* Put in something for proton collisions if not in database */
				if( AtmolCollRateCoeff[ipSpecies][ipPROTON].temps.size() == 0 )
				{
					if( AtmolCollRateCoeff[ipSpecies][ipHE_PLUS].temps.size() != 0 ) //He+
					{
						coll_temp.rate_coef_ul_set()[ipPROTON] = 2.0f * coll_temp.rate_coef_ul()[ipHE_PLUS];
					}
					else
						coll_temp.rate_coef_ul_set()[ipPROTON] = 1e-13f * (*(*tr).Lo()).g();

				}
				
#if	0
				/* if nothing else has been done, just put a small rate coefficient in */
				for( long intCollNo=0; intCollNo<ipNCOLLIDER; intCollNo++)
				{
					if( coll_temp.rate_coef_ul()[intCollNo] == 0. )
						coll_temp.rate_coef_ul_set()[intCollNo] = 1e-13;
				}
#endif
			}
			else
			{
				/* test for transitions without collision data */
				if( (*tr).Coll().rate_coef_ul_set()[ipELECTRON] == 0. )
				{
					if( atmdat.lgGbarOn && (*tr).Emis().gf() != 0.)
					{
						/* All transitions without collision data should use gbar if enabled */
						MakeCS(*tr);
					}
					else
					{
						//If gbar is off or no Aul, use the default collision strength value.
						coll_temp.col_str() = atmdat.collstrDefault;
					}
        				(*tr).Coll().is_gbar() = 1;

					coll_temp.rate_coef_ul_set()[ipELECTRON] = (COLL_CONST*coll_temp.col_str())/
						((*tr->Hi()).g()*phycon.sqrte);
				}

				//For some species col_str was showing up as zero in save line data when there was a non-zero coll rate
				if( tr->Coll().col_str() == -FLT_MAX && tr->Coll().rate_coef_ul()[ipELECTRON] > 0.0 )
				{
					tr->Coll().col_str() = (realnum)( tr->Coll().rate_coef_ul()[ipELECTRON] *
										((*tr->Hi()).g()*phycon.sqrte)/COLL_CONST);
				}
			}
		}
	}
}

/*Solving for the level populations*/

void dBase_solve()
{
	DEBUG_ENTRY( "dBase_solve()" );

	if( nSpecies==0 )
		return;

	static long maxNumLevels = 1;
	static vector<double> g, ex, pops, depart, source, sink;
	static multi_arr<double,2> AulEscp, AulDest, AulPump, CollRate;

	if( g.size() == 0 )
	{
		for( long ipSpecies=0; ipSpecies < nSpecies; ipSpecies++ )
			maxNumLevels = MAX2( maxNumLevels, dBaseSpecies[ipSpecies].numLevels_max );

		ASSERT( maxNumLevels > 0 );

		g.resize(maxNumLevels);
		ex.resize(maxNumLevels);
		pops.resize(maxNumLevels);
		depart.resize(maxNumLevels);
		source.resize(maxNumLevels);
		sink.resize(maxNumLevels);

		AulEscp.alloc(maxNumLevels, maxNumLevels);
		AulDest.alloc(maxNumLevels, maxNumLevels);
		AulPump.alloc(maxNumLevels, maxNumLevels);
		CollRate.alloc(maxNumLevels, maxNumLevels);
	}

	// zero all of these values
	vzero( g );
	vzero( ex );
	vzero( pops );
	vzero( depart );
	vzero( source );
	vzero( sink );
	AulEscp.zero();
	AulDest.zero();
	AulPump.zero();
	CollRate.zero();

	avx_ptr<double> arg(1,maxNumLevels), bstep(1,maxNumLevels);

	double totalHeating = 0.;
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		dBaseSpecies[ipSpecies].CoolTotal = 0.;

		if( !dBaseSpecies[ipSpecies].lgActive )
			continue;

#if 0
		//limit for now to small number of levels
		dBaseSpecies[ipSpecies].numLevels_local = MIN2( dBaseSpecies[ipSpecies].numLevels_local, 10 );
#endif

		// we always hit search phase first, reset number of levels
		if( conv.lgSearch )
			dBaseSpecies[ipSpecies].numLevels_local = dBaseSpecies[ipSpecies].numLevels_max;

		realnum abund = dBaseAbund(ipSpecies);

		const char *spName = dBaseSpecies[ipSpecies].chLabel;
		for( long ipLo = 0; ipLo < dBaseSpecies[ipSpecies].numLevels_local; ipLo++ )
		{
			/* statistical weights & Excitation Energies*/
			g[ipLo] = dBaseStates[ipSpecies][ipLo].g() ;
			// parts of the code assert that ground is at zero energy - this is
			// not true for the stored molecular data - so rescale to zero
			ex[ipLo] = dBaseStates[ipSpecies][ipLo].energy().WN() - 
					  dBaseStates[ipSpecies][0].energy().WN();
			/* zero some rates */	
			source[ipLo] = 0.;
			sink[ipLo] = 0.;
		}

		// non-zero was due to roundoff errors on 32-bit
		if( ex[0] <= dBaseStates[ipSpecies][0].energy().WN()* 10. *DBL_EPSILON )
			ex[0] = 0.;
		else
			TotalInsanity();

		for( long ipHi= 0; ipHi<dBaseSpecies[ipSpecies].numLevels_local; ipHi++)
		{
			for( long ipLo= 0; ipLo<dBaseSpecies[ipSpecies].numLevels_local; ipLo++)
			{
				AulEscp[ipHi][ipLo] = 0.;
			}
		}
		for( long ipHi= 0; ipHi<dBaseSpecies[ipSpecies].numLevels_local; ipHi++)
		{
			for( long ipLo= 0; ipLo<dBaseSpecies[ipSpecies].numLevels_local; ipLo++)
			{
				AulDest[ipHi][ipLo] = 0.;
			}
		}
		for( long ipLo= 0; ipLo<dBaseSpecies[ipSpecies].numLevels_local; ipLo++)
		{
			for( long ipHi= 0; ipHi<dBaseSpecies[ipSpecies].numLevels_local; ipHi++)
			{
				AulPump[ipLo][ipHi] = 0.;
			}
		}

		const bool isO1 = ( strcmp( spName,"O  1" ) == 0 ),
			isO3 = (strcmp(spName, "O  3") == 0),
			isCa2 = (strcmp(spName, "Ca 2") == 0),
			isN1 = (strcmp(spName, "N  1") == 0),
			isMg2 = (strcmp(spName, "Mg 2") == 0);

		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local || (*tr).ipCont() <= 0)
				continue;
			int ipLo = (*tr).ipLo();
			AulEscp[ipHi][ipLo] = (*tr).Emis().Aul()*(*tr).Emis().Pesc_total();
			AulDest[ipHi][ipLo] = (*tr).Emis().Aul()*(*tr).Emis().Pdest();
			AulPump[ipLo][ipHi] = (*tr).Emis().pump();
			AulPump[ipHi][ipLo] = AulPump[ipLo][ipHi]*g[ipLo]/g[ipHi];

			//Zero out xtra Ex and Dx rates before filling
			double xtraExRate = 0.;
			double xtraDxRate = 0.;

			//Set the extra Ex and Dx rates and add to AulPump
			if( isO1 )
			{
				static const double lo1_nrg = 0.;
				static const double lo2_nrg = 158.265;
				static const double lo3_nrg = 226.977;

				static const double hi1_nrg = 97488.378;
				static const double hi2_nrg = 97488.448;
				static const double hi3_nrg = 97488.538;

				if( ( fp_equal( tr->Lo()->energy().WN(),lo1_nrg ) ||
						fp_equal( tr->Lo()->energy().WN(),lo2_nrg ) ||
						fp_equal( tr->Lo()->energy().WN(),lo3_nrg ) ) &&
						( fp_equal( tr->Hi()->energy().WN(),hi1_nrg ) ||
						fp_equal( tr->Hi()->energy().WN(),hi2_nrg ) ||
						fp_equal( tr->Hi()->energy().WN(),hi3_nrg ) ) )
				{
					setXtraRatesO1(*tr, xtraExRate, xtraDxRate);
				}

				// add deexcitation rate from level 3 to all 3 below it
				if( tr->ipHi() == 3 )
					xtraDxRate += oxy.d6300/3;
			}
			else if( isO3 )
			{
				// add deexcitation rate from level 3 to all 3 below it
				if( tr->ipHi() == 3 )
					xtraDxRate += oxy.d5007r/3;
			}
			else if( isCa2 )
			{
				//Deexcitation of levels 2 through 5 to the ground state
				if( tr->ipLo() == 0 && tr->ipHi() < 5 )
					setXtraRatesCa2(*tr, xtraDxRate);
			}
			else if( isN1 )
			{
				//Deexcitation of levels 2 and 3 of N 1 to the ground state
				if( tr->ipLo() == 0 && (tr->ipHi() == 1 || tr->ipHi() == 2) )
					xtraDxRate += atoms.d5200r;
			}
			else if( isMg2 )
			{
				//Deexcitation of levels 4 and 5 of Mg 2
				if( tr->ipLo() == 0 && (tr->ipHi() == 4 || tr->ipHi() == 5) )
					xtraDxRate += atoms.rateMg2;
			}
			else if( strcmp(spName, "Fe 2") == 0 )
			{
				setXtraRatesFe2(*tr, xtraExRate, xtraDxRate);
			}

			AulPump[ipLo][ipHi] += xtraExRate;
			AulPump[ipHi][ipLo] += xtraDxRate;
		}

		/*Set all the heating and cooling to zero*/
		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			CollisionZero( (*tr).Coll() );
		}

		/*Updating the CollRate*/

		/*Set all the collision rates to zero*/
		for( long ipHi= 0; ipHi<dBaseSpecies[ipSpecies].numLevels_local; ipHi++)
		{
			for( long ipLo= 0; ipLo<dBaseSpecies[ipSpecies].numLevels_local; ipLo++)
			{
				CollRate[ipHi][ipLo] = 0.;
			}
		}

		ColliderDensities colld(colliders);
		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			int ipLo = (*tr).ipLo();
			CollRate[ipHi][ipLo] = (*tr).Coll().ColUL( colld );
		}

		for(int ipHi=1; ipHi<dBaseSpecies[ipSpecies].numLevels_local; ++ipHi)
			arg[ipHi] = -(ex[ipHi]-ex[ipHi-1])*T1CM / phycon.te;
		vexp( arg.ptr0(), bstep.ptr0(), 1, dBaseSpecies[ipSpecies].numLevels_local );
		for (int ipLo=0; ipLo<dBaseSpecies[ipSpecies].numLevels_local-1; ++ipLo)
		{
			double boltz_over_glo=1./ g[ipLo];
			for (int ipHi=ipLo+1; ipHi<dBaseSpecies[ipSpecies].numLevels_local; ++ipHi)
			{
				boltz_over_glo *= bstep[ipHi];
				CollRate[ipLo][ipHi] = CollRate[ipHi][ipLo] * g[ipHi] * boltz_over_glo;
			}
		}
			
		/* now add in excitations resulting from cosmic ray secondaries */
		// \todo 2 add branch to do forbidden transitions	
		// this g-bar only works for permitted lines

		/* get secondaries for all permitted lines by scaling LyA 
		 * excitation by ratio of cross section (oscillator strength/energy) 
		 * Born approximation or plane-wave approximation */
		double sfac =  secondaries.x12tot * 
			iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,0).EnergyWN() /
			iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,0).Emis().gf();

		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			if( (*tr).ipCont() > 0 )
			{
				int ipHi = (*tr).ipHi();
				if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
					continue;
				int ipLo = (*tr).ipLo();
				(*tr).Coll().rate_lu_nontherm_set() = sfac *
					((*tr).Emis().gf()/(*tr).EnergyWN());
					
				CollRate[ipLo][ipHi] += (*tr).Coll().rate_lu_nontherm();
				CollRate[ipHi][ipLo] += (*tr).Coll().rate_lu_nontherm() * g[ipLo] / g[ipHi];
			}
		}

		{
			// debug dump one line
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC && (nzone>110)/**/ )
			{

				//static bool runOnce = false;
				if( atmdat.ipSpecIon[ipIRON][1] == ipSpecies)
				{
					//runOnce = false;
					for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
						 tr != dBaseTrans[ipSpecies].end(); ++tr)
					{
						int ipLo = tr->ipLo();
						int ipHi = tr->ipHi();
						if( ipLo == 39 && ipHi == 139 )
						{
							fprintf(ioQQQ,"DUMPING LINE:%i:%i\n",ipLo+1,ipHi+1);
							DumpLine(*tr);
						}
						if( ipLo == 36 && ipHi == 133 )
						{
							fprintf(ioQQQ,"DUMPING LINE:%i:%i\n",ipLo+1,ipHi+1);
							DumpLine(*tr);
						}
					}
				}
			}
		}

		multi_arr<double,2> Cool(dBaseSpecies[ipSpecies].numLevels_local, dBaseSpecies[ipSpecies].numLevels_local);
		double grnd_excit = 0.0;
		double cooltl, coolder;
		int nNegPop;
		bool lgZeroPop, lgDeBug = false;

		/* solve the n-level atom */
		static Atom_LevelN atom_levelN;
		/* dBaseSpecies[ipSpecies].numLevels_local is the number of levels to compute*/ 

		atom_levelN(
			dBaseSpecies[ipSpecies].numLevels_local,
			/* ABUND is total abundance of species, used for nth equation
			 * if balance equations are homogeneous */
			abund, 
			/* g(dBaseSpecies[ipSpecies].numLevels_local) is statistical weight of levels */
			g, 
			/* EX(dBaseSpecies[ipSpecies].numLevels_local) is excitation potential of levels, deg K or wavenumbers
			 * 0 for lowest level, all are energy rel to ground NOT d(ENER) */
			ex, 
			/* this is 'K' for ex[] as Kelvin deg, is 'w' for wavenumbers */
			'w',
			/* populations [cm-3] of each level as deduced here */
			pops, 
			/* departure coefficient, derived below */
			depart,
			/* net transition rate, A * esc prob, s-1 */
			AulEscp, 
			/* AulDest[ihi][ilo] is destruction rate, trans from ihi to ilo, A * dest prob,
			 * asserts confirm that [ihi][ilo] is zero */
			AulDest, 
			/* AulPump[ilo][ihi] is pumping rate from lower to upper level (s^-1), (hi,lo) must be zero  */
			AulPump, 
			/* collision rates (s^-1), evaluated here and returned for cooling by calling function,
			 * unless following flag is true.  If true then calling function has already filled
			 * in these rates.  CollRate[ipSpecies][j] is rate from ipSpecies to j */
			CollRate,
			/* this is an additional creation rate from continuum, normally zero, units cm-3 s-1 */
			source,
			/* this is an additional destruction rate to continuum, normally zero, units s-1 */
			sink,
			/* total cooling and its derivative, set here but nothing done with it*/
			&cooltl, 
			&coolder, 
			/* string used to identify species in case of error */
			spName, 
			dBaseSpecies[ipSpecies].lgPrtMatrix,
			dBaseSpecies[ipSpecies].lgImgMatrix,
			/* nNegPop flag indicating what we have done
			 * positive if negative populations occurred
			 * zero if normal calculation done
			 * negative if too cold (for some atoms other routine will be called in this case) */
			&nNegPop,
			/* true if populations are zero, either due to zero abundance of very low temperature */
			&lgZeroPop,
			/* option to print debug information */
			lgDeBug,
			/* option to do the molecule in LTE */
			dBaseSpecies[ipSpecies].lgLTE,
			/* cooling per line */
			&Cool,
			NULL,
			/* excitation rate out of the ground state */
			&grnd_excit);


		if ( ! dBaseSpecies[ipSpecies].lgMolecular )
			ionbal.ExcitationGround
				[dBaseStates[ipSpecies][0].nelem()-1]
				[(*(*dBaseTrans[ipSpecies].begin()).Lo()).IonStg()-1]
			   += grnd_excit;


		if( nNegPop > 0 )
		{
			/* negative populations occurred */
			if( conv.lgSearch )
				fprintf(ioQQQ," dBase_solve fixup: negative pops set to zero during search phase, continuing.\n");
			else
			{
				fprintf(ioQQQ," PROBLEM in dBase_solve, atom_levelN returned negative population .\n");
				ContNegative();
				ShowMe();
				cdEXIT( EXIT_FAILURE );
			}
		}

		// highest levels may have no population
		while( (pops[dBaseSpecies[ipSpecies].numLevels_local-1]<=0 ) &&
			(dBaseSpecies[ipSpecies].numLevels_local > 1) )
				--dBaseSpecies[ipSpecies].numLevels_local;

		for( long j=0;j< dBaseSpecies[ipSpecies].numLevels_local; j++ )
		{
			dBaseStates[ipSpecies][j].Pop() = pops[j];
			dBaseStates[ipSpecies][j].DepartCoef() = depart[j];
			dBaseStates[ipSpecies][j].status() = LEVEL_ACTIVE;
		}
		for( long j=dBaseSpecies[ipSpecies].numLevels_local;
			j< dBaseSpecies[ipSpecies].numLevels_max; j++ )
		{
			dBaseStates[ipSpecies][j].Pop() = 0.;
			dBaseStates[ipSpecies][j].DepartCoef() = 0.;
			dBaseStates[ipSpecies][j].status() = LEVEL_INACTIVE;
		}

		/*Atmol  line*/
		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			int ipLo = (*tr).ipLo();
			(*tr).Coll().cool() = max(Cool[ipHi][ipLo],0.);
			(*tr).Coll().heat() = max(-Cool[ipHi][ipLo],0.);

			if ( (*tr).ipCont() > 0 )
			{
				/* population of lower level rel to ion, corrected for stim em */
				(*tr).Emis().PopOpc() = (*(*tr).Lo()).Pop() - (*(*tr).Hi()).Pop()*
					(*(*tr).Lo()).g()/(*(*tr).Hi()).g();

				set_xIntensity(*tr);
			
				/* it's possible for this sum to be zero.  Set ratio to zero in that case. */
				if( CollRate[ipLo][ipHi]+AulPump[ipLo][ipHi] > 0. )
				{
					(*tr).Emis().ColOvTot() = CollRate[ipLo][ipHi]/
						(CollRate[ipLo][ipHi]+AulPump[ipLo][ipHi]);
				}
				else
					(*tr).Emis().ColOvTot() = 0.;
				
				// this is only used for save line printout.  Maybe colliders may be involved, but
				// this simple approximation of a "collision strength" should be good enough for
				// the purposes of that printout.
				(*tr).Coll().col_str() = (realnum)( (*tr).Coll().rate_coef_ul()[ipELECTRON] *
					(g[ipHi]*phycon.sqrte)/COLL_CONST);
			}
		}

		//LyA destruction by Fe II
		if( ipSpecies == atmdat.ipSpecIon[ipIRON][1])
		{
			/* the hydrogen Lya destruction rate, then probability */
			hydro.dstfe2lya = 0.;
			/* count how many photons were removed-added */
			for( TransitionList::iterator tr = dBaseTrans[ipSpecies].begin(); tr != dBaseTrans[ipSpecies].end();++tr)
			{
				double exRate = 0.;
				double dxRate = 0.;
				setXtraRatesFe2(*tr,exRate,dxRate);
				hydro.dstfe2lya += (realnum)(
						tr->Lo()->Pop()*exRate -
						tr->Hi()->Pop()*dxRate);
			}

			/* the destruction prob comes from
			 * dest rate = n(2p) * A(21) * PDest */
			double pop = iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop();
			if( pop > SMALLFLOAT )
			{
				hydro.dstfe2lya /= (realnum)(pop * iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().Aul());
			}
			else
			{
				hydro.dstfe2lya = 0.;
			}
			/* NB - do not update hydro.dstfe2lya here if large FeII not on since
			 * done in FeII overlap */
		}

		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi < dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			EmLineZero( (*tr).Emis() );		
		}

		dBaseSpecies[ipSpecies].CoolTotal = cooltl;
		CoolAdd( dBaseSpecies[ipSpecies].chLabel, 0., max(0.,dBaseSpecies[ipSpecies].CoolTotal) );
		if( dBaseSpecies[ipSpecies].lgMolecular )
			thermal.elementcool[LIMELM] += max(0.,dBaseSpecies[ipSpecies].CoolTotal);
		else
			thermal.elementcool[dBaseStates[ipSpecies][0].nelem()-1] += max(0.,dBaseSpecies[ipSpecies].CoolTotal);
		totalHeating += max(0., -dBaseSpecies[ipSpecies].CoolTotal);
		thermal.dCooldT += coolder;

		/* option to print departure coefficients */
		{
			enum {DEBUG_LOC=false};

			if( DEBUG_LOC )
			{
				fprintf( ioQQQ, " Departure coefficients for species %li\n", ipSpecies );
				for( long j=0; j< dBaseSpecies[ipSpecies].numLevels_local; j++ )
				{
					fprintf( ioQQQ, " level %li \t Depar Coef %e\n", j, depart[j] );
				}
			}
		}
	}

	// total heating for all data base species
	thermal.setHeating(0,27,totalHeating);

	return;
}

/*Leiden*/
STATIC double LeidenCollRate(long ipSpecies, long ipCollider, const TransitionProxy& tr, double ftemp)
{
	DEBUG_ENTRY( "LeidenCollRate()" );
	double ret_collrate = InterpCollRate( AtmolCollRateCoeff[ipSpecies][ipCollider], tr.ipHi(), tr.ipLo(), ftemp);
	return ret_collrate;
}

/*STOUT*/
STATIC double StoutCollRate(long ipSpecies, long ipCollider, const TransitionProxy& tr, double ftemp)
{
	DEBUG_ENTRY( "StoutCollRate()" );

	double rate = 0.;
	int n = StoutCollData[ipSpecies].ntemps(tr.ipHi(),tr.ipLo(),ipCollider);
	if( n < 2 )
		return 0.;

	const double *x = StoutCollData[ipSpecies].temps(tr.ipHi(),tr.ipLo(),ipCollider);
	const double *y = StoutCollData[ipSpecies].collstrs(tr.ipHi(),tr.ipLo(),ipCollider);

	//If the temperature is above or below the temperature range, use the CS from the closest temperature.
	//Otherwise, do the linear interpolation.
	double fupsilon = 0.;
	if( ftemp < x[0] )
	{
		fupsilon = y[0];
	}
	else if( ftemp > x[n-1] )
	{
		fupsilon = y[n-1];
	}
	else
	{
		fupsilon = linint(&x[0],&y[0],n,ftemp);
	}

	ASSERT(fupsilon > 0);

	// deexcitation rate coefficient or collision strength?
	bool lgIsRate = StoutCollData[ipSpecies].lgIsRate(tr.ipHi(),tr.ipLo(),ipCollider);
	/* We can deal with deexcitation rate coefficients and collision strengths currently */
	if( lgIsRate )
	{
		rate = fupsilon;
		tr.Coll().col_str() = (rate * (*tr.Hi()).g()*sqrt(ftemp))/COLL_CONST;
	}
	else
	{
		/* convert the collision strength to a collision rate coefficient */
		/* This formula converting collision strength to collision rate coefficient works fine for the electrons*/
		/* For any other collider the mass would be different*/
		if( ipCollider == ipELECTRON )
		{
			rate = (COLL_CONST*fupsilon)/((*tr.Hi()).g()*sqrt(ftemp));
			tr.Coll().col_str() = fupsilon;
		}
		else
		{
			fprintf(ioQQQ,"PROBLEM: Stout data format does not support using collision strengths with "
					"non-electron colliders.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	return rate;
}

/*CHIANTI*/
STATIC double ChiantiCollRate(long ipSpecies, long ipCollider, const TransitionProxy& tr, double ftemp)
{
	DEBUG_ENTRY( "ChiantiCollRate()" );

	double rate = 0.;
	double fupsilon = CHIANTI_Upsilon(ipSpecies, ipCollider, tr.ipHi(), tr.ipLo(), ftemp);

	/* NB NB - if proton colliders, the upsilons returned here are actually already rate coefficients. */
	/* these are designated by a collider index and a transition type */
	if( ipCollider == ipPROTON )
	{
		rate = fupsilon;
	}
	else if( ipCollider == ipELECTRON )
	{
		/* convert the collision strength to a collision rate coefficient */
		/*This formula converting collision strength to collision rate coefficient works fine for the electrons*/
		/*For any other collider the mass would be different*/
		rate = (COLL_CONST*fupsilon)/((*tr.Hi()).g()*sqrt(ftemp));
	}
	else
		rate = 0.;

	return rate;
}

double CHIANTI_Upsilon(long ipSpecies, long ipCollider, long ipHi, long ipLo, double ftemp)
{
	double fdeltae,fscalingparam,fkte,fxt,fsups,fups;
	int intxs,inttype,intsplinepts;

	DEBUG_ENTRY( "CHIANTI_Upsilon()" );

	if( AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].collspline.size() == 0 )
	{
		return 0.;
	}

	intsplinepts = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].nSplinePts;
	inttype = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].intTranType;
	fdeltae = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].EnergyDiff;
	fscalingparam = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].ScalingParam;

	fkte = ftemp/fdeltae/1.57888e5;

	/*Way the temperature is scaled*/
	/*Burgess&Tully 1992:Paper gives only types 1 to 4*/
	/*Found that the expressions were the same for 5 & 6 from the associated routine DESCALE_ALL*/
	/*What about 7,8&9?*/
	if( inttype ==1 || inttype==4 )
	{
		fxt = 1-(log(fscalingparam)/(log(fkte+fscalingparam)));
	}
	else if(inttype  == 2 || inttype == 3||inttype == 5 || inttype == 6)
	{
		fxt = fkte/(fkte+fscalingparam);
	}
	else
		TotalInsanity();

	vector<double> xs( intsplinepts );

	/*Creating spline points array*/
	double* spl = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].collspline.data();
	for(intxs=0;intxs<intsplinepts;intxs++)
	{
		double coeff = (double)1/(intsplinepts-1);
		xs[intxs] = coeff*intxs;
		if(DEBUGSTATE)
		{
			printf("The xs and spl values are %f and %f \n",xs[intxs],spl[intxs]);
			getchar();
		}
	}

	double* xs2 = xs.data();

	const bool SPLINE_INTERP=false;
	if (! SPLINE_INTERP)
	{
		fsups = linint( xs2, spl, intsplinepts, fxt);
	}
	else
	{
		/*Finding the second derivative*/
		double *spl2 = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].SplineSecDer.data();
		
		if(DEBUGSTATE)
		{
			printf("\n");
			for(intxs=0;intxs<intsplinepts;intxs++)
			{
				printf("The %d value of 2nd derivative is %f \n",intxs+1,spl2[intxs]);
			}
		}
		
		/*Extracting out the value*/
		splint(xs2,spl,spl2,intsplinepts,fxt,&fsups);
	}

	/*Finding upsilon*/
	if(inttype == 1)
	{
		fups = fsups*log(fkte+EE);
	}
	else if(inttype == 2)
	{
		fups = fsups;
	}
	else if(inttype == 3)
	{
		fups = fsups/(fkte+1.0) ;
	}
	else if(inttype == 4)
	{
		fups = fsups*log(fkte+fscalingparam) ;
	}
	else if(inttype == 5)
	{
		fups = fsups/fkte ;
	}
	else if(inttype == 6)
	{
		fups = exp10(fsups) ;
	}
	else
	{
		TotalInsanity();
	}

	if( fups < 0. ) 
	{
		fprintf( ioQQQ," WARNING: Negative upsilon in species %s, collider %li, indices %4li %4li, Te = %e.\n",
				dBaseSpecies[ipSpecies].chLabel, ipCollider, ipHi, ipLo, ftemp );
		fups = 0.;
	}
	ASSERT(fups>=0);
	return(fups);
}

STATIC void setXtraRatesO1(const TransitionProxy& tr, double &xtraExRate, double &xtraDxRate)
{
	DEBUG_ENTRY("setXtraRatesO1()");

	double esab,
	  eslb,
	  esoi,
	  flb,
	  foi,
	  opaclb,
	  opacoi,
	  xlb,
	  xoi;
	double aoi = tr.Emis().Aul();
	double alb = iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH3p,ipH1s).Emis().Aul();

	fixit("ticket #78 refers");
	// The code below should be calculating the O I 1025 pumping by H Ly beta, as well as
	// the inverse process (this can become important in hydrogen-deficient environments).
	// It now uses the Elitzur & Netzer (1985, ApJ, 291, 464) theory, which is no longer
	// valid since the line overlap code prevents us from getting at the escape probability
	// of individual lines.

	/* A's from Pradhan; OI pump line; Ly beta, 8446 */

	// line overlap code makes this the escape probability of the combined lines
	esab = tr.Emis().Pesc_total();

	// these two are no longer correct, the line overlap code makes it impossible
	// to get at the escape probabilities of the individual lines...
	esoi = tr.Emis().Pesc_total();
	eslb = iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH3p,ipH1s).Emis().Pesc_total();

	/* all trace output turned on with "trace ly beta command' */
	if( trace.lgTr8446 && trace.lgTrace )
	{
		fprintf( ioQQQ,
			"       P8446 finds Lbeta, OI widths=%10.2e%10.2e and esc prob=%10.2e%10.2e esAB=%10.2e\n",
		  GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]), GetDopplerWidth(dense.AtomicWeight[ipOXYGEN]), eslb, esoi, esab );
	}

	/* find relative opacities for two lines */
	opacoi = 2.92e-9*dense.xIonDense[ipOXYGEN][0]*0.5556/GetDopplerWidth(dense.AtomicWeight[ipOXYGEN]);
	opaclb = 1.22e-8*iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop()/GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]);

	/* these are x sub a (OI) and x sub b (ly beta) defined in Elitz+Netz */
	xoi = opacoi/(opacoi + opaclb);
	xlb = opaclb/(opacoi + opaclb);

	/* find relative line-widths, assume same rest freq */
	foi = MIN2(GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]),GetDopplerWidth(dense.AtomicWeight[ipOXYGEN]))/GetDopplerWidth(dense.AtomicWeight[ipOXYGEN]);
	flb = MIN2(GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]),GetDopplerWidth(dense.AtomicWeight[ipOXYGEN]))/GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN])*
		MAX2(0.,1.- tr.Emis().Pesc_total());

	if( trace.lgTr8446 && trace.lgTrace )
	{
		fprintf( ioQQQ,
			"       P8446 opac Lb, OI=%10.2e%10.2e X Lb, OI=%10.2e%10.2e FLb, OI=%10.2e%10.2e\n",
		  opaclb, opacoi, xlb, xoi, flb, foi );
	}

	/* pumping of OI by L-beta - this goes into OI matrix as 1-5 rate
	 * lgInducProcess set false with no induced command, usually true */
	if( rfield.lgInducProcess )
	{
		xtraExRate += (realnum)((flb*alb*iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH3p].Pop()*
				xoi*(1. - esab)/dense.xIonDense[ipOXYGEN][0]));
		/* net decay rate from upper level */
		xtraDxRate += (realnum)(aoi*(1. - (1. - foi)*(1. - esoi) - xoi*(1. - esab)*foi));
	}
	else
	{
		xtraExRate += 0.;
		xtraDxRate += 0.;
	}
	return;
}

STATIC void setXtraRatesCa2(const TransitionProxy& tr, double &xtraDxRate )
{
	DEBUG_ENTRY( "setXtraRatesCa2()" );

	double hlgam,
	PhotoRate2,
	PhotoRate3,
	PhotoRate4,
	PhotoRate5;

	/* photoionization of evcited levels by Ly-alpha */
	hlgam = rfield.otslin[ iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).ipCont() -1];
	PhotoRate5 = 1.7e-18*hlgam;
	PhotoRate4 = 8.4e-19*hlgam;
	PhotoRate3 = 7.0e-18*hlgam;
	PhotoRate2 = 4.8e-18*hlgam;

	if( tr.ipHi() + 1 == 2)
	{
		xtraDxRate += PhotoRate2;
	}
	else if( tr.ipHi() + 1 == 3 )
	{
		xtraDxRate += PhotoRate3;
	}
	else if( tr.ipHi() + 1 == 4 )
	{
		xtraDxRate += PhotoRate4;
	}
	else if( tr.ipHi() + 1 == 5 )
	{
		xtraDxRate += PhotoRate5;
	}
	else
	{
		xtraDxRate += 0.;
	}
	return;
}
/*setXtraRatesFe2 find rate of Lya excitation of the FeII atom */
STATIC void setXtraRatesFe2(const TransitionProxy& tr, double &xtraExRate, double &xtraDxRate)
{
	DEBUG_ENTRY( "setXtraRatesFe2()" );

	if( ! hydro.lgLyaFeIIPumpOn )
		return;

	/* get rates FeII atom is pumped */

	/* elsewhere in this file the dest prob hydro.dstfe2lya is defined from
	 * quantites derived here, and the resulting populations */
	
	/*************trapeze form La profile:de,EnerLyaProf1,EnerLyaProf2,EnerLyaProf3,EnerLyaProf4*************************
	 * */
	/* width of Lya in cm^-1 */
	/* HLineWidth has units of cm/s, as was evaluated in PresTotCurrent */
	/* the factor is 1/2 of E(Lya, cm^-1_/c */
	double de = 1.37194e-06*hydro.HLineWidth*2.0/3.0;
	/* 82259 is energy of Lya in wavenumbers, so these are the form of the trapezoid */
	double EnerLyaProf1 = 82259.0 - de*2.0;
	double EnerLyaProf2 = 82259.0 - de;
	double EnerLyaProf3 = 82259.0 + de;
	double EnerLyaProf4 = 82259.0 + de*2.0;

	double PhotOccNumLyaCenter = 0.;

	/* find Lya photon occupation number */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop() > SMALLFLOAT )
	{
		/* This is the photon occupation number at the Lya line center */
		PhotOccNumLyaCenter =
			MAX2(0.,1.0-
			iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().Pesc_total())/
			(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop()/iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop()*3. - 1.0);
	}


	/* on first iteration optical depth in line is inward only, on later
	 * iterations is total optical depth */
	double taux = 0.;
	if( iteration == 1 )
	{
		taux = tr.Emis().TauIn();
	}
	else
	{
		taux = tr.Emis().TauTot();
	}

	/* the energy of the FeII line */
	double EnergyWN = tr.EnergyWN();

	if( EnergyWN >= EnerLyaProf1 && EnergyWN <= EnerLyaProf4  &&  taux > 1 )
	{
		/* this branch, line is within the Lya profile */

		/*
		 * Lya source function, at peak is PhotOccNumLyaCenter,
		 *
		 *     Prof2    Prof3
		 *       ----------
		 *      /          \
		 *     /            \
		 *    /              \
		 *  ======================
		 * Prof1              Prof4
		 *
		 */

		double PhotOccNum_at_nu = 0.,
			PumpRate = 0.;
		if( EnergyWN < EnerLyaProf2 )
		{
			/* linear interpolation on edge of trapazoid */
			PhotOccNum_at_nu = PhotOccNumLyaCenter*(EnergyWN - EnerLyaProf1)/ de;
		}
		else if( EnergyWN < EnerLyaProf3 )
		{
			/* this is the central plateau */
			PhotOccNum_at_nu = PhotOccNumLyaCenter;
		}
		else
		{
			/* linear interpolation on edge of trapazoid */
			PhotOccNum_at_nu = PhotOccNumLyaCenter*(EnerLyaProf4 - EnergyWN)/de;
		}

		/* at this point Lya source function at FeII line energy is defined, but
		 * we need to multiply by line width in Hz,
		 * >>refer	Fe2	pump	Netzer, H., Elitzur, M., & Ferland, G. J. 1985, ApJ, 299, 752-762*/

#if 0
		/** \todo 2 change this number to speed of light. */
		/* width of FeII line in Hz  */
		double FeIILineWidthHz = 1.e8/(EnergyWN*299792.5)*sqrt(log(taux))*GetDopplerWidth(dense.AtomicWeight[ipIRON]);

		/* final Lya pumping rate, s^-1*/
		PumpRate = FeIILineWidthHz * PhotOccNum_at_nu * tr.Emis().Aul()*
		  powi(82259.0f/EnergyWN,3);
#endif
		/* above must be bogus, use just occ num times A */
		PumpRate = tr.Emis().Aul()* PhotOccNum_at_nu;

		/* Lya pumping rate from ipHi to lower n */
		xtraDxRate += PumpRate;

		/* Lya pumping rate from n to upper ipHi */
		xtraExRate += PumpRate * (*tr.Hi()).g()/(*tr.Lo()).g();
	}

	return;
}
