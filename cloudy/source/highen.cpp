/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*highen do high energy radiation field - gas interaction, Compton scattering, etc */
#include "cddefines.h"
#include "highen.h"
#include "trace.h"
#include "heavy.h"
#include "radius.h"
#include "magnetic.h"
#include "hextra.h"
#include "thermal.h"
#include "ionbal.h"
#include "opacity.h"
#include "pressure.h"
#include "gammas.h"
#include "rfield.h"
#include "doppvel.h"
#include "dense.h"

void highen( void )
{
	long int i,
		ion,
		nelem;

	double CosRayDen, 
	  crsphi, 
	  heatin, 
	  sqthot;

	double hin;

	DEBUG_ENTRY( "highen()" );


	/**********************************************************************
	 *
	 *						COMPTON RECOIL IONIZATION
	 *
	 * bound electron scattering of >2.3 kev photons if neutral
	 * lgComptonOn usually t, f if "NO COMPTON EFFECT" command given
	 * lgCompRecoil usually t, f if "NO RECOIL IONIZATION" command given 
	 *
	 **********************************************************************/

	/* lgComptonOn turned false with no compton command,
	 * lgCompRecoil - no recoil ionization */
	if( rfield.lgComptonOn && ionbal.lgCompRecoil )
	{
		for( nelem=0; nelem<LIMELM; ++nelem )
		{
			for( ion=0; ion<nelem+1; ++ion )
			{
				double CompRecoilIonRate = 0.,
					CompRecoilHeatRate = 0.;
				if( dense.xIonDense[nelem][ion] > 0. )
				{
					/* recoil ionization starts at 194 Ryd = 2.6 keV */
					/* this is the ionization potential of the valence shell */
					/* >>chng 02 may 27, lower limit is now 1 beyond actual threshold
					 * since recoil energy at threshold was very small, sometimes negative */
					long limit = MIN2( rfield.nflux , rfield.nPositive );
					double valence_IP_Ryd = Heavy.Valence_IP_Ryd[nelem][ion];
					for( i=ionbal.ipCompRecoil[nelem][ion]; i < limit; ++i)
					{
						crsphi = opac.OpacStack[i-1+opac.iopcom] * rfield.SummedCon[i];

						/* direct hydrogen ionization due to compton scattering, 
						 * does not yet include secondaries,
						 * last term accounts for number of valence electrons that contribute */
						CompRecoilIonRate += crsphi;

						/* recoil energy in Rydbergs
						 * heating modified for suprathermal secondaries below; ANU2=ANU**2 */
						/* >>chng 02 mar 27, use general formula for recoil energy */
						/*energy = 2.66e-5*rfield.anu2(i) - 1.;*/
						double recoil_energy = rfield.anu2(i) * 
							( EN1RYD / 
							  ( ELECTRON_MASS * SPEEDLIGHT * SPEEDLIGHT) ) - valence_IP_Ryd;

						/* heating is in rydbergs because SecIon2PrimaryErg, SecExcitLya2PrimaryErg, HeatEfficPrimary in ryd */
						CompRecoilHeatRate += crsphi*recoil_energy;
					}
					/* net heating rate, per atom, convert ryd/sec/cm3 to ergs/sec/atom */
					CompRecoilHeatRate *= EN1RYD;

					/* this is number of electrons in valence shell of this ion */
					long nElec = ionbal.nCompRecoilElec[nelem-ion];

					ionbal.CompRecoilHeatRate[nelem][ion] = nElec*CompRecoilHeatRate;
					ionbal.CompRecoilIonRate[nelem][ion] = nElec*CompRecoilIonRate;

					ASSERT( ionbal.CompRecoilHeatRate[nelem][ion] >= 0.);
					ASSERT( ionbal.CompRecoilIonRate[nelem][ion] >= 0.);
				}
			}
		}
	}
	else
	{
		for( nelem=0; nelem<LIMELM; ++nelem )
		{
			for( ion=0; ion<nelem+1; ++ion )
			{
				ionbal.CompRecoilIonRate[nelem][ion] = 0.;
				ionbal.CompRecoilHeatRate[nelem][ion] = 0.;
			}
		}
	}

	/**********************************************************************
	 *
	 *                          COSMIC RAYS
	 *
	 * heating and ionization by cosmic rays, other relativistic particles
	 * CRYDEN=density (1/CM**3), neutral rate assumes 15ev total
	 * energy loss, 13.6 into ionization, 1.4 into heating
	 * units erg/sec/cm**3
	 * iff not specified, CRTEMP is 2.6E9K
	 *
	 **********************************************************************/

	if( hextra.cryden > 0. )
	{
		ASSERT( hextra.crtemp > 0. );
		/* this is current cosmic ray density, as determined from original density times
		 * possible dependence on radius */
		if( hextra.lg_CR_B_equipartition )
		{
			/* >>chng 06 jun 02, add this option
			 *  this is case where cr are in equipartition with magnetic field -
			 * set with COSMIC RAY EQUIPARTITION command */
			CosRayDen = hextra.background_density * 
				/* ratio of energy density in current B to typical galactic
				 * galactic background energy density of 1.8 eV cm-3 is from
				 *>>refer	cr	background	Webber, W.R. 1998, ApJ, 506, 329 */ 
				magnetic.energydensity / 
				(CR_EDEN_GAL_BACK_EV_CMM3/*1.8eV cm-3*/ * EN1EV/*erg eV-1*/ );
		}
		else
		{
			/* this is usual case, CR density may depend on radius, usually does not */
			CosRayDen = hextra.cryden*pow(radius.Radius/radius.rinner,(double)hextra.crpowr);
		}

		/* cosmic ray energy density rescaled by ratio to background ion rate and B field */
		hextra.cr_energydensity = CosRayDen/hextra.background_density * 
			(CR_EDEN_GAL_BACK_EV_CMM3/*1.8eV cm-3*/ * EN1EV/*erg eV-1*/ );

		/* related to current temperature, when thermal */
		sqthot = sqrt(hextra.crtemp);

		/* rate hot electrons heat cold ones, Balbus and McKee 1982
		 * units erg sec^-1 cm^-3,
		 * in sumheat we will multipy this rate by sum of neturals, but for this
		 * term we actually want eden, so mult by eden over sum of neut */
		ionbal.CosRayHeatThermalElectrons = 5.5e-14/sqthot*CosRayDen;

		/* ionization rate; Balbus and McKee */
		ionbal.CosRayIonRate = 1.22e-4/sqthot*
			log(2.72*pow(hextra.crtemp/1e8,0.097))*CosRayDen;

		/* option for thermal CRs, first is the usual (and default) relativistic case */
		if( hextra.crtemp > 2e9 )
		{
			/* usual circumstance; relativistic cosmic rays, 
			 * cosmic ray ionization rate s-1 cm-3; ext rel limit */
			ionbal.CosRayIonRate *= 3.;

		}
		else
		{
			/*  option for thermal cosmic rays */
			ionbal.CosRayIonRate *= 10.;
		}
		/* >>chng 04 jan 27, from 0.093 to 2.574 as per following */
		/* cr heating from Table 1 of
		 *>>refer	cr	heating	Wolfire et al.1995, ApJ, 443, 152
		 * For every ionization due to cosmic rays, ~35eV of heat is added 
		 * to the system.  This manifests itself in the ionbal.CosRayHeatNeutralParticles term
		 * by the 2.574*EN1RYD term, which is just the energy in ergs in 35 eV.
		 * Change made by Nick Abel and Gargi Shaw, 04 Jan 27.   In heatsum 
		 * we  multiply by the number of secondaries that occur */
		ionbal.CosRayHeatNeutralParticles = ionbal.CosRayIonRate*2.574*EN1RYD;

		if( trace.lgTrace )
		{
			fprintf( ioQQQ, "     highen: cosmic ray density;%10.2e CRion rate;%10.2e CR heat rate=;%10.2e CRtemp;%10.2e\n", 
			  CosRayDen, ionbal.CosRayIonRate, ionbal.CosRayHeatNeutralParticles, hextra.crtemp );
		}
	}
	else
	{
		ionbal.CosRayIonRate = 0.;
		ionbal.CosRayHeatNeutralParticles = 0.;
	}
	/* >>chng 06 may 23, Penning ionization
	ionbal.CosRayIonRate += 1e-9 * 
		iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe2s3S].Pop; */

	/*fprintf(ioQQQ,"DEBUG cr %.2f %.3e %.3e %.3e\n",
		fnzone,
		hextra.cryden ,
		ionbal.CosRayIonRate ,
		ionbal.CosRayHeatNeutralParticles );*/

	/**********************************************************************
	 *
	 * add on extra heating due to turbulence, goes into [1] of [x][0][11][0]
	 *
	 **********************************************************************/

	/* TurbHeat added with hextra command, DispScale added with turbulence dissipative */
	if( (hextra.TurbHeat+DoppVel.DispScale) != 0. )
	{
		/* turbulent heating only goes into the low-energy heat of this element */
		/* >>>>chng 00 apr 28, functional form of radius dependence had bee turrad/depth
		 * and so went to infinity at the illuminated face.  Changed to current form as
		 * per Ivan Hubeny comment */
		if( hextra.lgHextraDepth )
		{
			/* if turrad is >0 then vary heating with depth */
			ionbal.ExtraHeatRate = 
				hextra.TurbHeat*sexp(radius.depth /hextra.turrad);

			/* >>chng 00 aug 16, added option for heating from back side */
			if( hextra.turback != 0. )
			{
				ionbal.ExtraHeatRate += 
					hextra.TurbHeat*sexp((hextra.turback - radius.depth) /hextra.turrad);
			}
		}
		else if( hextra.lgHextraDensity )
		{
			/* depends on density */
			ionbal.ExtraHeatRate = 
				hextra.TurbHeat*dense.gas_phase[ipHYDROGEN]/hextra.HextraScaleDensity;
		}
		else if( hextra.lgHextraSS )
		{
			/* with SS disk model */
			ionbal.ExtraHeatRate = 
			  hextra.HextraSSalpha*pressure.PresGasCurr*sqrt(hextra.HextraSS_M*GRAV_CONST*
					powi((double)hextra.HextraSSradius,-3));
		}
		/* this is turbulence dissipate command */
		else if( DoppVel.DispScale > 0. )
		{
			double turb = DoppVel.TurbVel * sexp( radius.depth / DoppVel.DispScale );
			/* if turrad is >0 then vary heating with depth */
			/* >>chng 01 may 10, add extra factor of length over 1e13 cm */
			ionbal.ExtraHeatRate = 3.45e-28 / 2.82824 * turb * turb * turb *
					( dense.gas_phase[ipHYDROGEN] / 1e10 ) / (DoppVel.DispScale/1e13);
		}
		else
		{
			/* constant extra heating */
			ionbal.ExtraHeatRate = hextra.TurbHeat;
		}
	}

	else
	{
		ionbal.ExtraHeatRate = 0.;
	}

	/**********************************************************************
	 *
	 * option to add on fast neutron heating, goes into [0] & [2] of [x][0][11][0]
	 *
	 **********************************************************************/
	if( hextra.lgNeutrnHeatOn )
	{
		/* hextra.totneu is energy flux erg cm-2 s-1
		 * CrsSecNeutron is 4E-26 cm^-2, cross sec for stopping rel neutrons
		 * this is heating erg s-1 due to fast neutrons, assumed to secondary ionize */
		/* neutrons assumed to only secondary ionize */
		ionbal.xNeutronHeatRate = hextra.totneu*hextra.CrsSecNeutron;
	}
	else
	{
		ionbal.xNeutronHeatRate = 0.;
	}


	/**********************************************************************
	 *
	 * pair production in elec field of nucleus 
	 *
	 **********************************************************************/
	t_phoHeat photoHeat;
	ionbal.PairProducPhotoRate[0] = GammaK(opac.ippr,rfield.nflux,opac.ioppr,1.,&photoHeat);
	ionbal.PairProducPhotoRate[1] = photoHeat.HeatLowEnr;
	ionbal.PairProducPhotoRate[2] = photoHeat.HeatHiEnr;

	/**********************************************************************
	 *
	 * Compton energy exchange 
	 *
	 **********************************************************************/
	rfield.cmcool = 0.;
	rfield.cmheat = 0.;
	heatin = 0.;
	/* lgComptonOn usually t, turns off Compton */
	if( rfield.lgComptonOn )
	{
		for( i=0; i < rfield.nflux; i++ )
		{

			/* Compton cooling
			 * CSIGC is Tarter expression times ANU(I)*3.858E-25
			 * 6.338E-6 is k in inf mass Rydbergs, still needs factor of TE */
			rfield.comup[i] = (double)(rfield.flux[0][i]+rfield.ConInterOut[i]+
			  rfield.outlin[0][i]+ rfield.outlin_noplot[i])*rfield.csigc[i]*(dense.eden*4.e0/
			  TE1RYD*1e-15);

			rfield.cmcool += rfield.comup[i];

			/* Compton heating 
			 * CSIGH is Tarter expression times ANU(I)**2 * 3.858E-25
			 * CMHEAT is just spontaneous, HEATIN is just induced */
			rfield.comdn[i] = (double)(rfield.flux[0][i]+rfield.ConInterOut[i]+
			   rfield.outlin[0][i]+ rfield.outlin_noplot[i])*rfield.csigh[i]*dense.eden*1e-15;

			/* induced Compton heating */
			hin = (double)(rfield.flux[0][i]+rfield.ConInterOut[i]+rfield.outlin[0][i]+rfield.outlin_noplot[i])*
			  rfield.csigh[i]*rfield.OccNumbIncidCont[i]*dense.eden*1e-15;
			rfield.comdn[i] += hin;
			heatin += hin;

			/* following is total compton heating */
			rfield.cmheat += rfield.comdn[i];
		}

		/* remember how important induced compton heating is */
		if( rfield.cmheat > 0. )
		{
			rfield.cinrat = MAX2(rfield.cinrat,heatin/rfield.cmheat);
		}

		if( trace.lgTrace && trace.lgComBug )
		{
			fprintf( ioQQQ, "     HIGHEN: COOL num=%8.2e HEAT num=%8.2e\n", 
			  rfield.cmcool, rfield.cmheat );
		}
	}

	if( trace.lgTrace && trace.lgComBug )
	{
		fprintf( ioQQQ, 
			"     HIGHEN finds heating fracs= frac(compt)=%10.2e "
			" f(pair)%10.2e totHeat=%10.2e\n", 
		  rfield.cmheat/thermal.htot, 
		  thermal.heating(0,21)/thermal.htot,
		  thermal.htot	);
	}
	return;
}
