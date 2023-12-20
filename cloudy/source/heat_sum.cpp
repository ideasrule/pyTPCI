/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*HeatSum evaluate heating and secondary ionization for current conditions */
/*HeatZero is called by ConvBase */
#include "cddefines.h"
#include "thermal.h"
#include "heavy.h"
#include "trace.h"
#include "secondaries.h"
#include "conv.h"
#include "called.h"
#include "coolheavy.h"
#include "iso.h"
#include "h2.h"
#include "hmi.h"
#include "ionbal.h"
#include "phycon.h"
#include "numderiv.h"
#include "grainvar.h"
#include "deuterium.h"
#include "mole.h"
#include "dense.h"
#include "cooling.h"
#include "rfield.h"

/* this is the faintest relative heating we will print */
static const double FAINT_HEAT = 0.02;

static const bool PRT_DERIV = false;

STATIC void ElectronFractions(realnum& ElectronFraction, realnum& xValenceDonorDensity)
{
	if (1)
	{
		long ipNoble[] = { 
			ipHELIUM, ipNEON, ipARGON, ipKRYPTON
			/*, ipXENON, ipRADON, ipUNUNOCTIUM */ 
		};
		
		long ipFullShell = -1, ipNextFull = 0;
		xValenceDonorDensity = 0.;
		for (long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem)
		{
			// H -> 1; He -> 2; Li -> 1 (nelem==2,ipFullShell==1), etc.
			xValenceDonorDensity += 
				(nelem-ipFullShell)*dense.gas_phase[nelem];
			if (nelem == ipNoble[ipNextFull])
			{
				ipFullShell = nelem;
				++ipNextFull;
			}
		}
		
		//fprintf(ioQQQ,"%g %g\n",dense.eden,xValenceDonorDensity);

		ElectronFraction = (realnum)(MAX2(MIN2(1.0,dense.eden/
															xValenceDonorDensity),1e-4));
	}
	else
	{
		/* >>chng 03 apr 29, move evaluation of xNeutralParticleDensity to here 
		 * from PresTotCurrent since only used for secondary ionization */
		/* this is total neutral particle density for collisional ionization 
		 * for very high ionization conditions it may be zero */
		/** total number of neutral colliders (atoms & molecules) */
		realnum xNeutralParticleDensity = 0.;
		for( long nelem=0; nelem < LIMELM; nelem++ )
		{
			xNeutralParticleDensity += dense.xIonDense[nelem][0];
		}
		
		xNeutralParticleDensity += deut.xIonDense[0];
		
		/* now add all the heavy molecules */
		for(long i=0; i < mole_global.num_calc; i++ )
		{
			/* add in if real molecule and not counted above */
			if(mole.species[i].location == NULL && mole_global.list[i]->charge == 0 && mole_global.list[i]->isIsotopicTotalSpecies())
				xNeutralParticleDensity += (realnum)mole.species[i].den;
		}
		
		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
				fprintf(ioQQQ," xIonDense xNeutralParticleDensity tot\t%.3e",xNeutralParticleDensity);
				for( long nelem=0; nelem < LIMELM; nelem++ )
				{
					fprintf(ioQQQ,"\t%.2e",dense.xIonDense[nelem][0]);
				}
				fprintf(ioQQQ,"\n");
			}
		}
		xValenceDonorDensity = (xNeutralParticleDensity+dense.EdenTrue);
		ElectronFraction = (realnum)(MAX2(MIN2(1.0,dense.EdenTrue/
															xValenceDonorDensity),1e-4));
	}
}

void SecIoniz( void )
{
	DEBUG_ENTRY( "SecIoniz()" );
	/* secondary ionization and excitation due to Compton scattering */
	double cosmic_ray_ionization_rate , 
		pair_production_ionization_rate , 
		fast_neutron_ionization_rate , SecExcitLyaRate;

	/* ionization and excitation rates from hydrogen itself */
	double SeconIoniz_iso[NISO] , 
		SeconExcit_iso[NISO];

	long ion,
		limit;

	double secmet ,
		smetla;
	long ipISO,
		ns;
	double secmetsave[LIMELM];

	realnum SaveOxygen1 , SaveCarbon1;

	/* routine to sum up total heating, and print agents if needed
	 * it also computes the true derivative, dH/dT */
	double Cosmic_ray_sec2prim;
	realnum sec2prim_par_1;
	realnum sec2prim_par_2;

	/*******************************************************************
	 *
	 * reevaluate the secondary ionization and excitation rates 
	 *
	 *******************************************************************/
	/* ElectronFraction is electron fraction 
	 * analytic fits to 
	 * >>>refer	sec	ioniz	Shull & Van Steenberg (1985; Ap.J. 298, 268).
	 * lgSecOFF turns off secondary ionizations, sets heating efficiency to 100% */

	/* at very low ionization - as per 
	 * >>>refer	sec	ioniz	Xu & McCray 1991, Ap.J. 375, 190.
	 * everything goes to asymptote that is not present in Shull and
	 * Van Steenberg - do this by not letting ElecFrac fall below 1e-4 */
	/* this uses the correct electron density, which will not equal the
	 * current electron density if we are not converged.  Using the 
	 * correct value aids convergence onto it */
	
	// ElectronFraction should be the fraction of electrons available
	// for collisions which are free, as opposed to bound in atoms,
	// ions or molecules, which seems like the most plausible
	// generalization of the definition of X in S&vS.

	realnum xValenceDonorDensity, ElectronFraction;

	ElectronFractions(ElectronFraction, xValenceDonorDensity);

	// xBoundValenceDensity must be consistent with ElectronFraction
	// for expressions for ionizations below to have bounded behaviour
	// in the high ElectronFraction limit, i.e. have a factor of (1-EF)
	realnum xBoundValenceDensity = (1.0f-ElectronFraction)*xValenceDonorDensity;

	if( secondaries.lgSecOFF )
	{
		secondaries.HeatEfficPrimary = 1.;
		secondaries.SecIon2PrimaryErg = 0.;
		secondaries.SecExcitLya2PrimaryErg = 0.;
		Cosmic_ray_sec2prim = 0.05f;
	}
	/* >>chng 03 apr 29, only evaluate one time per zone since drove oscillations
	 * in He+ - He0 ionization front in very high Z models, like hizqso */
	else 
	{

		/* this is heating efficiency for high-energy photo ejections.  It is the ratio
		 * of the final heat added to the energy of the photo electron.  Fully
		 * ionized, this is unity, and less than unity for neutral media.  
		 * It is a scale factor that multiplies the
		 * high energy heating rate */
		secondaries.HeatEfficPrimary = 0.9971f*(1.f - pow(1.f - pow(ElectronFraction,(realnum)0.2663f),(realnum)1.3163f));

		/* secondary ionizations per primary Rydberg - the number of secondary 
		 * ionizations produced for each Rydberg of high energy photo-electron 
		 * energy deposited.  It multiplies the high-energy heating rate.  
		 * It is zero for a fully ionized gas and goes to 0.3908 for very neutral gas */
		secondaries.SecIon2PrimaryErg = 0.3908f*pow(1.f - pow(ElectronFraction,(realnum)0.4092f),(realnum)1.7592f);
		/* by dividing by the energy of one Rydberg, converts heating rate 
		 * in ergs into secondary ionization rate */
		secondaries.SecIon2PrimaryErg /= ((realnum)EN1RYD);

		/* This is cosmic ray secondaries per primary (???), 
		 * derived to approximate the curve given in Shull and
		 * van Steenberg for cosmic rays at 35 eV */		
		sec2prim_par_1 = -(1.251f + 195.932f * ElectronFraction);
		sec2prim_par_2 = 1.f + 46.814f * ElectronFraction - 44.969f * 
			ElectronFraction * ElectronFraction;

		Cosmic_ray_sec2prim = (exp(sec2prim_par_1/SDIV( sec2prim_par_2)));

		/*  number of secondary excitations of L-alpha per erg of high
		 * energy light - actually all Lyman lines 
		 *
		 *  Note--This formula is derived for primary energies greater 
		 *  than 100 eV and is only an approximation.  This will
		 *  over predict the secondary ionization of L-alpha.  We cannot make
		 *  fitting functions for this equation at low energies like we did
		 *  for the heating efficiency and for the number of secondaries
		 *  because the Shull and van Steenberg paper did not publish a similar
		 *  curve for L-alpha */
		secondaries.SecExcitLya2PrimaryErg = 0.4766f*pow(1.f - pow(ElectronFraction,(realnum)0.2735f),(realnum)1.5221f)/((realnum)EN1RYD);

		if( (trace.lgTrace && trace.lgSecIon) )
		{
			fprintf( ioQQQ, 
				" csupra H0 %.2e  HeatSum eval sec ion effic, ElectronFraction = %.3e HeatEfficPrimary %.3e SecIon2PrimaryErg %.3e\n", 
				secondaries.csupra[ipHYDROGEN][0],
				ElectronFraction,
				secondaries.HeatEfficPrimary ,
				secondaries.SecIon2PrimaryErg );
		}
	}

	/* add in effects of high-energy opacity of CO using C and O atomic opacities
	 * k-shell and valence photo of C and O in CO is not explicitly counted elsewhere
	 * this trick roughly accounts for it*/
	SaveOxygen1 = dense.xIonDense[ipOXYGEN][0];
	SaveCarbon1 = dense.xIonDense[ipCARBON][0];

	/* atomic C and O will include CO during the heating sum loop */
	fixit("This breaks the invariant for total mass.");
			/*
			* In any case, is CO the only species which contributes?
			* -- might expect all other molecular species to do so,
			* i.e. the item added should perhaps be
			* dense.xMolecules[nelem]) -- code duplicated in
			* opacity_addtotal
			*/
	dense.xIonDense[ipOXYGEN][0] += (realnum) findspecieslocal("CO")->den;
	dense.xIonDense[ipCARBON][0] += (realnum) findspecieslocal("CO")->den;
	
	/* metals secondary ionization, Lya excitation */
	secmet = 0.;
	smetla = 0.;
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		SeconIoniz_iso[ipISO] = 0.;
		SeconExcit_iso[ipISO] = 0.;
	}

	/* this loop includes hydrogenic, which is special case due to presence
	 * of substantial excited state populations */
	for( long nelem=0; nelem<LIMELM; ++nelem)
	{
		secmetsave[nelem] = 0.;
		
		if( dense.lgElmtOn[nelem] )
		{
			/* sum heat over all stages of ionization that exist */
			/* first do the iso sequences,
			 * h-like and he-like are special because full atom always done, 
			 * can be substantial,  pops in excited states, with little in ground 
			 * (true near LTE), these are done in following loop */

			limit = MIN2( dense.IonHigh[nelem] , nelem+1-NISO );

			/* loop over all elements/ions done with as two-level systems */
			for( ion=dense.IonLow[nelem]; ion < limit; ion++ )
			{
				/* this will be heating for a single stage of ionization */
				double HeatingHi = 0.;
				for( ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
				{
					/* heating by various sub-shells */
					HeatingHi += ionbal.PhotoRate_Shell[nelem][ion][ns][2];
				}

				/* secondary ionization rate,  */
				secmetsave[nelem] += HeatingHi*secondaries.SecIon2PrimaryErg* dense.xIonDense[nelem][ion];

				/* LA excitation rate, =0 if ionized, s-1 cm-3 */
				smetla += HeatingHi*secondaries.SecExcitLya2PrimaryErg* dense.xIonDense[nelem][ion];
			}
			secmet += secmetsave[nelem];

			if( secondaries.SecIon2PrimaryErg > SMALLFLOAT && xBoundValenceDensity > SMALLFLOAT )
			{
				/* prevent crash in very high ionization conditions 
				 * where xBoundValenceDensity is zero */
				/* secondary ionization rate,  */
				
				/* this branch loop over all ions done with full iso sequence */
				limit = MAX2( limit, dense.IonLow[nelem] );
				for( ion=MAX2(0,limit); ion < dense.IonHigh[nelem]; ion++ )
				{
					/* this is the iso sequence */
					ipISO = nelem-ion;
					ASSERT(ipISO < NISO);
					/* this will be heating for a single stage of ionization */
					double HeatingHi = 0.;
					/* the outer shell contains the Compton recoil part */
					for( ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
					{
						/* heating, erg s-1 atom-1, by low energy, and then high energy, photons */
						HeatingHi += ionbal.PhotoRate_Shell[nelem][ion][ns][2];
					}
					
					SeconIoniz_iso[ipISO] += HeatingHi*secondaries.SecIon2PrimaryErg* 
						iso_sp[ipISO][nelem].st[0].Pop()/
						xBoundValenceDensity;
					
					/* LA excitation rate, =0 if ionized, units excitations s-1 */
					SeconExcit_iso[ipISO] += HeatingHi*secondaries.SecExcitLya2PrimaryErg* 
						iso_sp[ipISO][nelem].st[0].Pop()/
						xBoundValenceDensity;
					
					ASSERT( SeconIoniz_iso[ipISO]>=0. && 
							  SeconExcit_iso[ipISO]>=0.);
				}
			}
		}
	}
	if( trace.lgTrace && trace.lgSecIon )
	{
		double savemax=0.;
		long int ipsavemax=-1;
		for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		{
			if( secmetsave[nelem] > savemax )
			{
				savemax = secmetsave[nelem];
				ipsavemax = nelem;
			}
		}
		fprintf( ioQQQ, 
			"   HeatSum secmet largest contributor element %li frac of total %.3e, total %.3e\n", 
			ipsavemax,
			savemax/SDIV(secmet),
			secmet);
	}

	/* now reset the abundances */
	dense.xIonDense[ipOXYGEN][0] = SaveOxygen1;
	dense.xIonDense[ipCARBON][0] = SaveCarbon1;

	/* convert secmet to proper final units */
	/*fprintf(ioQQQ,"bugggg\t%li\t%.3e\t%.3e\t%.3e\n", nzone , 
		secmet , xBoundValenceDensity , secmet / xBoundValenceDensity );*/
	/* prevent crash when xBoundValenceDensity is zero */
	if( secondaries.SecIon2PrimaryErg > SMALLFLOAT && xBoundValenceDensity > SMALLFLOAT )
	{
		/* convert from s-1 cm-3 to s-1 - a true rate */
		secmet /= xBoundValenceDensity;
		smetla /= xBoundValenceDensity;
	}
	else
	{
		/* zero */
		secmet = 0.;
		smetla = 0.;
	}

	/* add on cosmic rays - most important in neutral or molecular gas, use
	 * difference between total H and H+ density as substitute for sum of
	 * H in atoms and all molecules, but only in limit where difference is
	 * significant */
#	if 0
	double hneut;
	if( (dense.gas_phase[ipHYDROGEN] - dense.xIonDense[ipHYDROGEN][1])/
		dense.gas_phase[ipHYDROGEN]<0.001 )
	{
		/* limit where most H is ionized - simply use atomic and H2 */
		hneut = dense.xIonDense[ipHYDROGEN][0] + 2.*(findspecieslocal("H2")->den+findspecieslocal("H2*")->den);
	}
	else
	{
		/* limit where neutral is significant, use different between total and H+ */
		hneut = dense.gas_phase[ipHYDROGEN] - dense.xIonDense[ipHYDROGEN][1];
	}
#	endif


	/*******************************************************************
	 *
	 * secondary ionization and excitation rates 
	 *
	 *******************************************************************/

	/* the terms cosmic_ray_ionization_rate & SecExcitLyaRate contain energies added in highen.  
	 * The only significant one is usually bound Compton heating except when 
	 * cosmic rays are present */

	if( secondaries.SecIon2PrimaryErg > SMALLFLOAT && xBoundValenceDensity > SMALLFLOAT )
	{
		/* add on ionization rate s-1 cm-3 due to pair production */
		/* next two terms account for secondary ionization and excitation due to pair productions.
		 * The code integrates over the radiation field where \gamma <-> e- e+ is possible -
		 * that is ionbal.PairProducPhotoRate.  For all SEDs of realistic sources (AGN, or x-ray binaries)
		 * the radiation field over this range is small compared with the FUV/x-ray.  Pair production will,
		 * as a result, only be important energetically in well-shielded regions where much of the SED
		 * is absorbed.  In these regions gas would be neutral or molecular and the pair would produce
		 * the same effects as a cosmic ray - heating, excitation, and ionization.  That is why the code is
		 * here. Given the shape of the SED these effects are the largest pair production should produce -
		 * even then it is negligible.
		 */

		/* the photon rates in ionbal.PairProducPhotoRate are the integral over the range that the
		 * do pair production using cross sections from Figure 41 of Experimental Nuclear Physics,
		 * Vol1, E.Segre, ed
		 * factor of four in DensityRatio accounts for additional nucleons in alpha particle */
		realnum DensityRatio = (dense.gas_phase[ipHYDROGEN] + 
			4.F*dense.gas_phase[ipHELIUM])/xBoundValenceDensity;

		pair_production_ionization_rate = 
			ionbal.PairProducPhotoRate[2]*secondaries.SecIon2PrimaryErg*DensityRatio;

		/* total secondary excitation rate cm-3 s-1 for Lya due to pair production*/
		SecExcitLyaRate = 
			ionbal.PairProducPhotoRate[2]*secondaries.SecExcitLya2PrimaryErg*1.3333*DensityRatio;

		/* ionization rate of fast neutrons */
		fast_neutron_ionization_rate = 
			ionbal.xNeutronHeatRate*secondaries.SecIon2PrimaryErg*DensityRatio;

		/* Secondary excitation rate due to neutrons */
		SecExcitLyaRate += 
			ionbal.xNeutronHeatRate*secondaries.SecExcitLya2PrimaryErg*1.3333*DensityRatio;

		/* cosmic ray Lya excitation rate */
		SecExcitLyaRate += 
			/* multiply by atomic H and He densities */
			ionbal.CosRayHeatNeutralParticles*secondaries.SecExcitLya2PrimaryErg*1.3333*DensityRatio;
	}
	else
	{
		/* zero */
		pair_production_ionization_rate = 0.;
		SecExcitLyaRate = 0.;
		fast_neutron_ionization_rate = 0.;
	}

	/* cosmic ray ionization */
	/* >>>chng 99 apr 29, term in PhotoRate was not present */
	cosmic_ray_ionization_rate = 
		/* this term is H^0 cosmic ray ionization, set in highen, did not multiply by
		 * collider density in highen, so do not divide by it here */
		/* >>chng 99 jun 29, added SecIon2PrimaryErg to cosmic ray rate, also multiply
		 * by number of secondaries per primary*/
		ionbal.CosRayIonRate*Cosmic_ray_sec2prim +
		/* this is the cosmic ray heating rate */
		ionbal.CosRayHeatNeutralParticles*secondaries.SecIon2PrimaryErg;

	/* find total suprathermal collisional ionization rate
	 * CSUPRA is H0 col ionize rate from non-thermal electrons (inverse sec)
	 * x12tot is LA excitation rate, units s-1
	 * SECCMP evaluated in HIGHEN, =ioniz rate for cosmic rays, sec bound */

	/* option to force secondary ionization rate, normally false */
	if( secondaries.lgCSetOn )
	{
		for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		{
			for( ion=0; ion<nelem+1; ++ion )
			{
				secondaries.csupra[nelem][ion] = secondaries.SetCsupra*secondaries.csupra_effic[nelem][ion];
			}
		}
		secondaries.x12tot = secondaries.SetCsupra;
	}
	else
	{
		double csupra = cosmic_ray_ionization_rate + 
			pair_production_ionization_rate +
			fast_neutron_ionization_rate +
			SeconIoniz_iso[ipH_LIKE] + 
			SeconIoniz_iso[ipHE_LIKE] + 
			secmet;

		/** \todo	2	find correct high-energy limit for these */
		/* now fill in ionization rates for all elements and ion stages */
		for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		{
			for( ion=0; ion<nelem+1; ++ion )
			{
				secondaries.csupra[nelem][ion] = (realnum)csupra*secondaries.csupra_effic[nelem][ion];
			}
		}

		fixit("why does x12tot include non Ly-a stuff here, and get used to scale other lines elsewhere?");
		/* secondary suprathermal excitation of Ly-alpha, excitations s-1 */
		secondaries.x12tot = SecExcitLyaRate + 
			SeconExcit_iso[ipH_LIKE] + 
			SeconExcit_iso[ipHE_LIKE] + 
			smetla;

	}

	return;
}

void HeatSum( void )
{
	/* use to dim some vectors */
	const int NDIM = 40;

	long int ion,
	  ipnt, 
	  ipsave[NDIM], 
	  j, 
	  jpsave[NDIM], 
	  limit;

	double HeatingLo ,
		HeatingHi;
	long ipISO,
		ns;
	static long int nhit=0, 
	  nzSave=0;
	double photo_heat_2lev_atoms,
		photo_heat_ISO_atoms ,
		photo_heat_UTA ,
		OtherHeat , 
		oldfac, 
		save[NDIM];
	static double oldheat=0., 
	  oldtemp=0.;

	realnum SaveOxygen1 , SaveCarbon1;

	/* routine to sum up total heating, and print agents if needed
	 * it also computes the true derivative, dH/dT */
	double Cosmic_ray_heat_eff;

	DEBUG_ENTRY( "HeatSum()" );

	if ( lgConvBaseHeatTest )
	{
		fixit("much or all of this secondary stuff should be updated much earlier: solvers in conv_base use outdated details");
		SecIoniz();
	}

	/*******************************************************************
	 *
	 * reevaluate the secondary ionization and excitation rates 
	 *
	 *******************************************************************/
	/* ElectronFraction is electron fraction 
	 * analytic fits to 
	 * >>>refer	sec	ioniz	Shull & Van Steenberg (1985; Ap.J. 298, 268).
	 * lgSecOFF turns off secondary ionizations, sets heating efficiency to 100% */

	/* at very low ionization - as per 
	 * >>>refer	sec	ioniz	Xu & McCray 1991, Ap.J. 375, 190.
	 * everything goes to asymptote that is not present in Shull and
	 * Van Steenberg - do this by not letting ElecFrac fall below 1e-4 */
	/* this uses the correct electron density, which will not equal the
	 * current electron density if we are not converged.  Using the 
	 * correct value aids convergence onto it */
	
	// ElectronFraction should be the fraction of electrons available
	// for collisions which are free, as opposed to bound in atoms,
	// ions or molecules, which seems like the most plausible
	// generalization of the definition of X in S&vS.

	realnum xValenceDonorDensity, ElectronFraction;

	ElectronFractions(ElectronFraction, xValenceDonorDensity);

	// xBoundValenceDensity must be consistent with ElectronFraction
	// for expressions for ionizations below to have bounded behaviour
	// in the high ElectronFraction limit, i.e. have a factor of (1-EF)
	realnum xBoundValenceDensity = (1.0f-ElectronFraction)*xValenceDonorDensity;

	if( secondaries.lgSecOFF )
	{
		Cosmic_ray_heat_eff = 0.95;
	}
	/* >>chng 03 apr 29, only evaluate one time per zone since drove oscillations
	 * in He+ - He0 ionization front in very high Z models, like hizqso */
	else 
	{
		/* cosmic ray heating, counts as non-ionizing heat since already result of cascade,
		 * this was 35 eV per secondary ionization */
		/* We want to use the heating efficiency that is applicable to a 35 eV
		 * primary electron, the current efficiency used is for 100eV cosmic ray
		 * Here we use the correct value as given in Wolfire et al. 1995 */

		/* In general the equation for the cosmic ray heating rate is:*
		 *															  *
		 *															  *
		 * CR_Heat_Rate = (density)*(cosmic ray ionization rate)*     *
		 *				  (energy of primary electron)*               *
		 *				  (Heating efficiency at that energy)		  *
		 *															  *
		 * The product of the last two terms gives the amount of heat *
		 * added by the primary electron, and is dependent upon the   *
		 * electron fraction.  We are using the same average primary  *
		 * electron energy as Wolfire et al. (1995).  We do not,	  *
		 * however, use their formula for the heating efficiency.     *
		 * This is because the coefficients (f1, f2, and f3) were     *
		 * originally intended for primary electron energies >100eV.  *
		 * Instead we derived a heating efficiency based on the curves*
		 * given in Shull and van Steenberg (1985).  We interpolated  *
		 * for an energy of 35 eV the heating efficiency for electron *
		 * fractions between 1e-4 and 1								  *
		 **************************************************************/

		/*printf("Here is ElectronFraction:\t%.3e\n", ElectronFraction);*/
		Cosmic_ray_heat_eff = - 8.189 - 18.394 * ElectronFraction - 6.608 * ElectronFraction * ElectronFraction * log(ElectronFraction)
			+ 8.322 * exp( ElectronFraction )  + 4.961 * sqrt(ElectronFraction);
	}

	/*******************************************************************
	 *
	 * get total heating from all species
	 *
	 *******************************************************************/

	/* get total heating */
	photo_heat_2lev_atoms = 0.;
	photo_heat_ISO_atoms = 0.;
	photo_heat_UTA = 0.;

	/* add in effects of high-energy opacity of CO using C and O atomic opacities
	 * k-shell and valence photo of C and O in CO is not explicitly counted elsewhere
	 * this trick roughly accounts for it*/
	SaveOxygen1 = dense.xIonDense[ipOXYGEN][0];
	SaveCarbon1 = dense.xIonDense[ipCARBON][0];

	/* atomic C and O will include CO during the heating sum loop */
	fixit("This breaks the invariant for total mass.");
			/*
			* In any case, is CO the only species which contributes?
			* -- might expect all other molecular species to do so,
			* i.e. the item added should perhaps be
			* dense.xMolecules[nelem]) -- code duplicated in
			* opacity_addtotal
			*/
	dense.xIonDense[ipOXYGEN][0] += (realnum) findspecieslocal("CO")->den;
	dense.xIonDense[ipCARBON][0] += (realnum) findspecieslocal("CO")->den;

	/* this will hold cooling due to metal collisional ionization */
	CoolHeavy.colmet = 0.;
	
	/* this loop includes hydrogenic, which is special case due to presence
	 * of substantial excited state populations */
	for( long nelem=0; nelem<LIMELM; ++nelem)
	{
		thermal.heavycollcool[nelem] = 0.;
		
		if( dense.lgElmtOn[nelem] )
		{
			/* sum heat over all stages of ionization that exist */
			/* first do the iso sequences,
			 * h-like and he-like are special because full atom always done, 
			 * can be substantial,  pops in excited states, with little in ground 
			 * (true near LTE), these are done in following loop */

			limit = MIN2( dense.IonHigh[nelem] , nelem+1-NISO );

			/* loop over all elements/ions done with as two-level systems */
			for( ion=dense.IonLow[nelem]; ion < limit; ion++ )
			{
				/* this will be heating for a single stage of ionization */
				HeatingLo = 0.;
				HeatingHi = 0.;
				for( ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
				{
					/* heating by various sub-shells */
					HeatingLo += ionbal.PhotoRate_Shell[nelem][ion][ns][1];
					HeatingHi += ionbal.PhotoRate_Shell[nelem][ion][ns][2];
				}

				/* total photoelectric heating, all shells, for this stage 
				 * of ionization */
				thermal.setHeating(nelem,ion, dense.xIonDense[nelem][ion]* 
					(HeatingLo + HeatingHi*secondaries.HeatEfficPrimary) );
				/* heating due to UTA ionization */
				photo_heat_UTA += dense.xIonDense[nelem][ion]* 
					ionbal.UTA_heat_rate[nelem][ion];

				/* add to total heating */
				photo_heat_2lev_atoms += thermal.heating(nelem,ion);
				/*if( nzone>290 && thermal.heating(nelem,ion)/thermal.htot>0.01 )
					fprintf(ioQQQ,"buggg\t%li %li %.3f\n", nelem,ion,thermal.heating(nelem,ion)/thermal.htot);*/

				/* Cooling due to collisional ionization of all elements by thermal electrons
				 * CollidRate[nelem][ion][1] is cooling, erg/s/atom, evaluated in ion_collis
				 *
				 * we count energy from the continuum, so bound states have negative energy and free
				 * electrons have ~kT energy.  During a collisional ionization the kinetic energy of
				 * the first free electron is converted into the ionization energy of the ion plus the
				 * kinetic energies of the two resulting free electrons,
				 * KE(1 electron before collision) = IP + KE(2 electrons after collision)
				 * The change in kinetic energy, the cooling, is IP. */
				CoolHeavy.colmet += dense.xIonDense[nelem][ion]*ionbal.CollIonRate_Ground[nelem][ion][1];
				thermal.heavycollcool[nelem] += dense.xIonDense[nelem][ion]*ionbal.CollIonRate_Ground[nelem][ion][1];

			}
	
			/* this branch loop over all ions done with full iso sequence */
			limit = MAX2( limit, dense.IonLow[nelem] );
			for( ion=MAX2(0,limit); ion < dense.IonHigh[nelem]; ion++ )
			{
				/* this is the iso sequence */
				ipISO = nelem-ion;
				/* this will be heating for a single stage of ionization */
				HeatingLo = 0.;
				HeatingHi = 0.;
				/* the outer shell contains the Compton recoil part */
				for( ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
				{
					/* heating, erg s-1 atom-1, by low energy, and then high energy, photons */
					HeatingLo += ionbal.PhotoRate_Shell[nelem][ion][ns][1];
					HeatingHi += ionbal.PhotoRate_Shell[nelem][ion][ns][2];
				}

				/* net heating, erg cm-3 s-1 */
				thermal.setHeating(nelem,ion, iso_sp[ipISO][nelem].st[0].Pop()*
					(HeatingLo + HeatingHi*secondaries.HeatEfficPrimary) );

				/* heating due to UTA ionization, erg cm-3 s-1 */
				photo_heat_UTA += dense.xIonDense[nelem][ion]* 
					ionbal.UTA_heat_rate[nelem][ion];

				/* add to total heating */
				photo_heat_ISO_atoms += thermal.heating(nelem,ion);
			}

			/* make sure stages of ionization with no abundances,
			 * also has no heating */
			for( ion=0; ion<dense.IonLow[nelem]; ion++ )
			{
				ASSERT( thermal.heating(nelem,ion) == 0. );
			}
			for( ion=dense.IonHigh[nelem]; ion<nelem+1; ion++ )
			{
				ASSERT( thermal.heating(nelem,ion) == 0. );
			}
		}
	}

	/* now reset the abundances */
	dense.xIonDense[ipOXYGEN][0] = SaveOxygen1;
	dense.xIonDense[ipCARBON][0] = SaveCarbon1;

	/* >>chng 01 dec 20, do full some over all secondaries */
	/* bound Compton recoil heating */
	/* >>chng 02 mar 28, save heating in this var rather than heating(0,18) 
	 * since now saved in photo heat 
	 * this is only used for a printout and in lines, not as heat since already counted*/
	ionbal.CompRecoilHeatLocal = 0.;
	for( long nelem=0; nelem<LIMELM; ++nelem )
	{
		for( ion=0; ion<nelem+1; ++ion )
		{
			ionbal.CompRecoilHeatLocal += 
				ionbal.CompRecoilHeatRate[nelem][ion]*secondaries.HeatEfficPrimary*dense.xIonDense[nelem][ion];
		}
	}
	/* >>chng 05 nov 26, include this term - bound ionization of H2 
	 * assume total cs is that of two separated H */
	ionbal.CompRecoilHeatLocal += 
		2.*ionbal.CompRecoilHeatRate[ipHYDROGEN][0]*secondaries.HeatEfficPrimary*hmi.H2_total;

	/* find heating due to charge transfer 
	 * >>chng 05 oct 29, move from here to conv base so that can be treated as cooling if 
	 * negative heating */
	thermal.setHeating(0,24,thermal.char_tran_heat);

	/* save compton heating rate into main heating save array, 
	 * no secondary ionizations from compton heating cooling */
	thermal.setHeating(0,19, rfield.cmheat);

	/* heating due to pair production */
	thermal.setHeating(0,21, 
		ionbal.PairProducPhotoRate[2]*secondaries.HeatEfficPrimary*(dense.gas_phase[ipHYDROGEN] + 4.F*dense.gas_phase[ipHELIUM]) );
	/* last term above is number of nuclei in helium */

	/* this is heating due to fast neutrons, assumed to secondary ionize */
	thermal.setHeating(0,20,
		ionbal.xNeutronHeatRate*secondaries.HeatEfficPrimary*dense.gas_phase[ipHYDROGEN] +
		/* turbulent heating, assumed to be a non-ionizing heat agent, added here */
		ionbal.ExtraHeatRate );

	/* UTA heating */
	thermal.setHeating(0,7, photo_heat_UTA);
	/*fprintf(ioQQQ,"DEBUG UTA heat %.3e\n", photo_heat_UTA/SDIV(thermal.htot ));*/

	// diatomic photoionization heating 
	thermal.setHeating(0,18,0.);
	{
		double heat = 0.;
		for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		{
			/*>>KEYWORD	H2 photoionization heating */
			heat += (*diatom)->Abund() *
				( (*diatom)->photo_heat_soft + (*diatom)->photo_heat_hard*secondaries.HeatEfficPrimary);
		}
		thermal.setHeating(0,18,heat);
	}
	/** \todo	1	add part of hard heat to secondaries */

	/* >>chng 05 nov 27, approximate heating due to H2+, H3+ photoionization
	 * assuming H0 rates */
	/*>>KEYWORD	H2+ photoionization heating; H3+ photoionization heating */
	thermal.setHeating(0,26, (findspecieslocal("H2+")->den+findspecieslocal("H3+")->den) *
		(ionbal.PhotoRate_Shell[ipHYDROGEN][0][0][1] + 
		ionbal.PhotoRate_Shell[ipHYDROGEN][0][0][2]*secondaries.HeatEfficPrimary) );

	/* add on cosmic rays - most important in neutral or molecular gas, use
	 * difference between total H and H+ density as substitute for sum of
	 * H in atoms and all molecules, but only in limit where difference is
	 * significant */
#	if 0
	double hneut;
	if( (dense.gas_phase[ipHYDROGEN] - dense.xIonDense[ipHYDROGEN][1])/
		dense.gas_phase[ipHYDROGEN]<0.001 )
	{
		/* limit where most H is ionized - simply use atomic and H2 */
		hneut = dense.xIonDense[ipHYDROGEN][0] + 2.*(findspecieslocal("H2")->den+findspecieslocal("H2*")->den);
	}
	else
	{
		/* limit where neutral is significant, use different between total and H+ */
		hneut = dense.gas_phase[ipHYDROGEN] - dense.xIonDense[ipHYDROGEN][1];
	}
#	endif

	/* cosmic ray heating */
	thermal.setHeating(1,6,
		ionbal.CosRayHeatNeutralParticles*
		xBoundValenceDensity * Cosmic_ray_heat_eff +
		/* cosmic ray heating of thermal electrons */
		/* >>refer	CR	physics	Ferland, G.J. & Mushotzky, R.F., 1984, ApJ, 286, 42 */
		ionbal.CosRayHeatThermalElectrons * dense.eden );

	/* now sum up all heating agents not included in sum over normal bound levels above */
	OtherHeat = 0.;
	for( long nelem=0; nelem<LIMELM; ++nelem)
	{
		/* we now have ionization solution for this element,sum heat over
		 * all stages of ionization that exist */
		/* >>>chng 99 may 08, following loop had started at nelem+3, and so missed [1][0],
		 * which is excited states of hydrogenic species.  increase this limit */
		for( long i=nelem+1; i < LIMELM; i++ )
		{
			OtherHeat += thermal.heating(nelem,i);
		}
	}

	thermal.htot = OtherHeat + photo_heat_2lev_atoms + photo_heat_ISO_atoms;

	/* following checks whether total heating is strange, if we are not in search phase */
	if( called.lgTalk && !conv.lgSearch )
	{
		/* print this warning if not constant temperature and neg heat */
		if( thermal.htot < 0. && !thermal.lgTemperatureConstant)
		{
			fprintf( ioQQQ, 
				" Total heating is <0; is this model collisionally ionized? zone is %li\n",
				nzone );
		}
		else if( thermal.htot == 0. )
		{
			fprintf( ioQQQ, 
				" Total heating is 0; is the density small? zone is %li\n",
				nzone);
		}
	}

	/* add on line heating to this array, heatl was evaluated in sumcool
	 * not zero, because of roundoff error */
	if( thermal.heatl/thermal.ctot < -1e-15 )
	{
		fprintf( ioQQQ, " HeatSum gets negative line heating,%10.2e%10.2e this is insane.\n", 
		  thermal.heatl, thermal.ctot );

		fprintf( ioQQQ, " this is zone%4ld\n", nzone );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/*******************************************************************
	 *
	 * get derivative of total heating
	 *
	 *******************************************************************/

	/* now get derivative of heating due to photoionization, 
	 * >>chng 96 jan 14
	 *>>>>>NB heating decreasing with increasing temp is negative dH/dT */
	thermal.dHeatdT = -0.7*(photo_heat_2lev_atoms+photo_heat_ISO_atoms+
		photo_heat_UTA)/phycon.te;
	/* >>chng 04 feb 28, add this correction factor - when gas totally neutral heating
	 * does not depend on temperature - when ionized depends on recom coefficient - this
	 * smoothly joins two limits */
	thermal.dHeatdT *= dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN];
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dhdt  0 %.3e  %.3e  %.3e \n", 
		thermal.dHeatdT,
		photo_heat_2lev_atoms,
		photo_heat_ISO_atoms);

	/* oldfac was factor used in old implementation
	 * any term using it should be rethought */
	oldfac = -0.7/phycon.te;

	/* carbon monoxide */
	thermal.dHeatdT += thermal.heating(0,9)*oldfac;

	/* Ly alpha collisional heating */
	thermal.dHeatdT += thermal.heating(0,10)*oldfac;

	/* line heating */
	thermal.dHeatdT += thermal.heating(0,22)*oldfac;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dhdt  1 %.3e \n", thermal.dHeatdT);

	/* free free heating, propto T^-0.5
	 * >>chng 96 oct 30, from heating(1,12) to CoolHeavy.brems_heat_total,
	 * do cooling separately assume CoolHeavy.brems_heat_total goes as T^-3/2
	 * dHTotDT = dHTotDT + heating(1,12) * (-0.5/te) */
	thermal.dHeatdT += CoolHeavy.brems_heat_total*(-1.5/phycon.te);

	/* >>chng 04 aug 07, use better estimate for heating derivative; needed in PDRs, PvH */
	/* this includes PE, thermionic, and collisional heating of the gas by the grains */
	thermal.dHeatdT += gv.dHeatdT;

	/* helium triplets heating */
	thermal.dHeatdT += thermal.heating(1,2)*oldfac;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dhdt  2 %.3e \n", thermal.dHeatdT);

	/* hydrogen molecule collisional deexcitation heating */
	/* >>chng 04 jan 25, add max to prevent cooling from entering here */
	/*thermal.dHeatdT += MAX2(0.,thermal.heating(0,8))*oldfac;*/
	if( hmi.HeatH2Dexc_used > 0. )
		thermal.dHeatdT += hmi.deriv_HeatH2Dexc_used;

	/* >>chng 96 nov 15, had been oldfac below, wrong sign
	 * H- H minus heating, goes up with temp since rad assoc does */
	thermal.dHeatdT += thermal.heating(0,15)*0.7/phycon.te;

	/* H2+ heating */
	thermal.dHeatdT += thermal.heating(0,16)*oldfac;

	/* Balmer continuum and all other excited states
	 * - T goes up so does pop and heating
	 * but this all screwed up by change in eden */
	thermal.dHeatdT += thermal.heating(0,1)*oldfac;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dhdt  3 %.3e \n", thermal.dHeatdT);

	/* all three body heating, opposite of coll ion cooling */
	thermal.dHeatdT += thermal.heating(0,3)*oldfac;

	/* bound electron recoil heating
	thermal.dHeatdT += ionbal.CompRecoilHeatLocal*oldfac; */

	/* Compton heating */
	thermal.dHeatdT += thermal.heating(0,19)*oldfac;

	/* extra heating source, usually zero */
	thermal.dHeatdT += thermal.heating(0,20)*oldfac;

	/* pair production */
	thermal.dHeatdT += thermal.heating(0,21)*oldfac;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dhdt  4 %.3e \n", thermal.dHeatdT);

	/* derivative of heating due to collisions of H lines, heating(1,24) */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem<LIMELM; ++nelem)
		{
			thermal.dHeatdT += iso_sp[ipISO][nelem].dLTot;
		}
	}

	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dhdt  5 %.3e \n", thermal.dHeatdT);

	/* possibility of getting empirical heating derivative
	 * normally false, set true with 'set numerical derivatives' command */
	if( NumDeriv.lgNumDeriv )
	{
		if( ((nzone > 2 && nzone == nzSave) && 
		     ! fp_equal( oldtemp, phycon.te )) && nhit > 4 )
		{
			/* hnit is number of tries on this zone - use to stop numerical problems
			 * do not evaluate numerical derivative until well into solution */
			double deriv = (oldheat - thermal.htot)/(oldtemp - phycon.te);
			thermal.dHeatdT = deriv;
		}

		/* this happens when new zone starts */
		if( nzone != nzSave )
		{
			nhit = 0;
		}
		nzSave = nzone;
		nhit += 1;
		oldheat = thermal.htot;
		oldtemp = phycon.te;
	}

	if( trace.lgTrace && trace.lgHeatBug )
	{
		ipnt = 0;
		/* this loops through the 2D array that contains all agents counted in htot */
		for( long i=0; i < LIMELM; i++ )
		{
			for( j=0; j < LIMELM; j++ )
			{
				if( thermal.heating(i,j)/thermal.htot > FAINT_HEAT )
				{
					ipsave[ipnt] = i;
					jpsave[ipnt] = j;
					save[ipnt] = thermal.heating(i,j)/thermal.htot;
					/* increment ipnt, but do not let it get too big */
					ipnt = MIN2((long)NDIM-1,ipnt+1);
				}
			}
		}

		/* now check for possible line heating not counted in 1,23
		 * there should not be any significant heating source here
		 * since they would not be counted in derivative correctly */
		for(long i=0; i < thermal.ncltot; i++ )
		{
			if( thermal.heatnt[i]/thermal.htot > FAINT_HEAT )
			{
				ipsave[ipnt] = -1;
				jpsave[ipnt] = i;
				save[ipnt] = thermal.heatnt[i]/thermal.htot;
				ipnt = MIN2((long)NDIM-1,ipnt+1);
			}
		}

		fprintf( ioQQQ, 
		  "    HeatSum HTOT %.3e Te:%.3e dH/dT%.2e other %.2e photo 2lev %.2e photo iso %.2e\n", 
		  thermal.htot, 
		  phycon.te, 
		  thermal.dHeatdT ,
		  /* total heating is sum of following three terms
		   * OtherHeat is a sum over all other heating agents */
		  OtherHeat , 
		  photo_heat_2lev_atoms, 
		  photo_heat_ISO_atoms);

		fprintf( ioQQQ, "  " );
		for(long i=0; i < ipnt; i++ )
		{
			/*lint -e644 these three are initialized above */
			fprintf( ioQQQ, "   [%ld][%ld]%6.3f",
				ipsave[i], 
				jpsave[i],
				save[i] );
			/*lint +e644 these three are initialized above */
			/* throw a cr every n numbers */
			if( !(i%8) && i>0 && i!=(ipnt-1) )
			{
				fprintf( ioQQQ, "\n  " );
			}
		}
		fprintf( ioQQQ, "\n" );
	}
	return;
}

/* =============================================================================*/
/* HeatZero is called by ConvBase */
void HeatZero( void )
{
	DEBUG_ENTRY( "HeatZero()" );

	for(long i=0; i < LIMELM; i++ )
	{
		for(long j=0; j < LIMELM; j++ )
		{
			thermal.setHeating(i,j,0.);
		}
	}
	return;
}
