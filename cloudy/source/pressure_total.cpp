/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PresTotCurrent determine the gas and line radiation pressures for current conditions,
 * this sets values of pressure.PresTotlCurr, also calls tfidle */
#include "cddefines.h"
#include "taulines.h"
#include "opacity.h"
#include "hydrogenic.h"
#include "conv.h"
#include "iso.h"
#include "wind.h"
#include "magnetic.h"
#include "phycon.h"
#include "thermal.h"
#include "h2.h"
#include "dense.h"
#include "dynamics.h"
#include "trace.h"
#include "rt.h"
#include "pressure.h"
#include "radius.h"
#include "rfield.h"
#include "doppvel.h"
#include "rt_escprob.h"

/* this sets values of pressure.PresTotlCurr, also calls tfidle */
void PresTotCurrent()
{
	static long int
	  /* used in debug print statement to see which line adds most pressure */
	  ipLineTypePradMax=-1 , 
	  /* points to line causing radiation pressure */
	  ipLinePradMax=-LONG_MAX,
	  ip2=-LONG_MAX,
	  ip3=-LONG_MAX,
	  ip4=-LONG_MAX;

	/* the line radiation pressure variables that must be preserved since
	 * a particular line radiation pressure may not be evaluated if it is
	 * not important */
	static double Piso_seq[NISO]={0.,0.},
		PLevel2=0.,
		PHfs=0.,
		PCO=0.,
		P_H2=0.,
		P_dBase=0.;

	double 
		RadPres1, 
		TotalPressure_v, 
		pmx;

	DEBUG_ENTRY( "PresTotCurrent()" );

	if( !conv.nTotalIoniz )
	{
		/* resetting ipLinePradMax, necessary for 
		 * multiple cloudy calls in single coreload. */ 
		ipLinePradMax=-LONG_MAX;
		//pressure.PresTotlCurr = 0.;
	}

	/* PresTotCurrent - evaluate total pressure, dyne cm^-2
	 * and radiative acceleration */

	/* several loops over the emission lines for radiation pressure and
	 * radiative acceleration are expensive due to the number of lines and
	 * the memory they occupy.  Many equations of state do not include
	 * radiation pressure or radiative acceleration.  Only update these
	 * when included in EOS - when not included only evaluate once per zone.
	 * do it once per zone since we will still report these quantities.
	 * this flag says whether we must update all terms */
	bool lgMustEvaluate = false;

	/* this says we already have an evaluation in this zone so do not 
	 * evaluate terms known to be trivial, even when reevaluating major terms */
	bool lgMustEvaluateTrivial = false;
	/* if an individual agent is larger than this fraction of the total current
	 * radiation pressure then it is not trivial 
	 * conv.PressureErrorAllowed is relative error allowed in pressure */
	double TrivialLineRadiationPressure = conv.PressureErrorAllowed/10. * 
		pressure.pres_radiation_lines_curr;

	/* reevaluate during search phase when conditions are changing dramatically */
	if( conv.lgSearch )
	{
		lgMustEvaluate = true;
		lgMustEvaluateTrivial = true;
	}

	/* reevaluate if zone or iteration has changed */
	static long int nzoneEvaluated=0, iterationEvaluated=0;
	if( nzone!=nzoneEvaluated || iteration!=iterationEvaluated )
	{
		lgMustEvaluate = true;
		lgMustEvaluateTrivial = true;
		/* this flag says which source of radiation pressure was strongest */
		ipLineTypePradMax = -1;
	}

	/* constant pressure or dynamical sim - we must reevaluate terms 
	 * but do not need to reevaluate trivial contributors */
	if( (strcmp(dense.chDenseLaw,"WIND") == 0 ) ||
		 (strcmp(dense.chDenseLaw,"DYNA") == 0 ) ||
		(strcmp(dense.chDenseLaw,"CPRE") == 0 ) )
		lgMustEvaluate = true;

	if( lgMustEvaluate )
	{
		/* we are reevaluating quantities in this zone and iteration,
		 * remember that we did it */
		nzoneEvaluated = nzone;
		iterationEvaluated = iteration;
	}

	/* update density - temperature variables which may have changed */
	/* must call TempChange since ionization has changed, there are some
	* terms that affect collision rates (H0 term in electron collision) */
	TempChange(phycon.te , false);

	ASSERT(lgElemsConserved());

	/* evaluate the total ionization energy on a scale where neutral atoms
	 * correspond to an energy of zero, so the ions have a positive energy */
	phycon.EnergyIonization = 0.;
#if 0
	for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		for( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
		{
			/* lgMETALS is true, set false with "no advection metals" command */
			int kadvec = dynamics.lgMETALS;
			/* this gives the iso sequence for this ion - should it be included in the
			 * advected energy equation? lgISO is true, set false with 
			 * "no advection H-like" and "he-like" - for testing*/
			ipISO = nelem-ion;
			fixit(); /* should this be kadvec = kadvec && dynamics.lgISO[ipISO]; ? */
			if( ipISO >= 0 && ipISO<NISO )
				kadvec = dynamics.lgISO[ipISO];
			for( long i=1; i<=ion; ++i )
			{
				long int nelec = nelem+1 - i;
				/* this is the sum of all the energy needed to bring the atom up
				 * to the ion+1 level of ionization - at this point a positive number */
				phycon.EnergyIonization += dense.xIonDense[nelem][ion] * 
					t_ADfA::Inst().getEthresh(Heavy.nsShells[nelem][i-1],nelec+1,nelem+1)/EVRYD*kadvec;
			}
		}
	}
	/* convert to ergs/cm^3 */
	phycon.EnergyIonization *= EN1RYD;
#endif

	/** \todo	2	this is the total binding energy of the molecules, and 
	 * is negative, the energy need to get back to free atoms 
	 * never set and only appears in print statements */
	phycon.EnergyBinding = 0.;

	/* find total number of atoms and ions */
	SumDensities();

	/* the current gas pressure */
	pressure.PresGasCurr = dense.pden*phycon.te*BOLTZMANN;
	/*fprintf(ioQQQ,"DEBUG gassss %.2f %.4e %.4e %.4e \n", 
		fnzone, pressure.PresGasCurr , dense.pden , pressure.PresInteg );*/

	/* >>chng 01 nov 02, move to here from dynamics routines */
	/* >>chng 02 jun 18 (rjrw), fix to use local values */
	/* WJH 21 may 04, now recalculate wind v for the first zone too.
	 * This is necessary when we are forcing the sub or supersonic branch */
	if(!(wind.lgBallistic() || wind.lgStatic()))
	{
		wind.windv = DynaFlux(radius.depth)/(dense.xMassDensity);
	}

	/* this is the current ram pressure, at this location */
	pressure.PresRamCurr = dense.xMassDensity*POW2( wind.windv );

	/* this is the current turbulent pressure, not now added to equation of state 
	 * >>chng 06 mar 29, add Heiles_Troland_F / 6. term*/
	pressure.PresTurbCurr = dense.xMassDensity*POW2( DoppVel.TurbVel ) *
		DoppVel.Heiles_Troland_F / 6.;

	/** \todo	0	add this press term due to cosmic rays - hextra.cr_energydensity */

	/* radiative acceleration, lines and continuum */
	if( lgMustEvaluate )
	{
		/* continuous radiative acceleration */
		double rforce = 0.;
		double relec = 0.;
		for( long i=0; i < rfield.nflux; i++ )
		{
			rforce += (rfield.flux[0][i] + rfield.outlin[0][i] + rfield.outlin_noplot[i]+ rfield.ConInterOut[i])*
				rfield.anu(i)*(opac.opacity_abs[i] + opac.opacity_sct[i]);

			/* electron scattering acceleration - used only to derive force multiplier */
			relec += 
				(rfield.flux[0][i] + rfield.outlin[0][i] + rfield.outlin_noplot[i]+ 
				rfield.ConInterOut[i]) * 
				opac.OpacStack[i-1+opac.iopcom]*dense.eden*rfield.anu(i);
		}
		/* total continuum radiative acceleration */
		wind.AccelCont = (realnum)(rforce*EN1RYD/SPEEDLIGHT/dense.xMassDensity);

		wind.AccelElectron = (realnum)(relec*EN1RYD/SPEEDLIGHT/dense.xMassDensity);

		/* line acceleration; xMassDensity is gm per cc */
		wind.AccelLine = (realnum)(RT_line_driving()/SPEEDLIGHT/dense.xMassDensity);

		/* total acceleration cm s^-2 */
		wind.AccelTotalOutward = wind.AccelCont + wind.AccelLine;

		/* find wind.fmul - the force multiplier; wind.AccelElectron can be zero for very low 
		 * densities - fmul is of order unity - wind.AccelLine and wind.AccelCont 
		 * are both floats to will underflow long before wind.AccelElectron will - fmul is only used
		 * in output, not in any physics */
		if( wind.AccelElectron > SMALLFLOAT )
			wind.fmul = (realnum)( (wind.AccelLine + wind.AccelCont) / wind.AccelElectron);
		else
			wind.fmul = 0.;

		/* inward acceleration of gravity cm s^-2 */
		if (wind.comass == 0.)
		{
			wind.AccelGravity = 0.;
		}
		else
		{
			double reff = radius.Radius-radius.drad/2.;
			wind.AccelGravity = (realnum)(
				/* wind.comass - mass of star in solar masses */
				GRAV_CONST*wind.comass*SOLAR_MASS/POW2(reff)*
				/* wind.DiskRadius normally zero, set with disk option on wind command */
				(1.-wind.DiskRadius/reff) );
		}

#		if 0
		if( fudge(-1) )
			fprintf(ioQQQ,"DEBUG pressure_total updates AccelTotalOutward to %.2e grav %.2e\n",
			wind.AccelTotalOutward , wind.AccelGravity );
#		endif
	}

	/* must always evaluate H La width since used elsewhere */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().PopOpc() > 0. )
	{
		hydro.HLineWidth = (realnum)RT_LineWidth(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s),GetDopplerWidth(dense.AtomicWeight[ipHYDROGEN]));
	}
	else
		hydro.HLineWidth = 0.;


	/* find max rad pressure even if capped
	 * lgLineRadPresOn is turned off by NO RADIATION PRESSURE command */
	if( trace.lgTrace )
	{
		fprintf(ioQQQ,
			"     PresTotCurrent nzone %li iteration %li lgMustEvaluate:%c "
			"lgMustEvaluateTrivial:%c " 
			"pressure.lgLineRadPresOn:%c " 
			"rfield.lgDoLineTrans:%c \n", 
			nzone , iteration , TorF(lgMustEvaluate) , TorF(lgMustEvaluateTrivial),
			TorF(pressure.lgLineRadPresOn), TorF(rfield.lgDoLineTrans) );
	}

	if( lgMustEvaluate && pressure.lgLineRadPresOn && rfield.lgDoLineTrans )
	{
		/* RadPres is pressure due to lines, lgPres_radiation_ON turns off or on */
		pressure.pres_radiation_lines_curr = 0.;
		/* used to remember largest radiation pressure source */
		pmx = 0.;
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			if( lgMustEvaluateTrivial || Piso_seq[ipISO] > TrivialLineRadiationPressure )
			{
				Piso_seq[ipISO] = 0.;
				for( long nelem=ipISO; nelem < LIMELM; nelem++ )
				{
					/* does this ion stage exist? */
					if( dense.IonHigh[nelem] >= nelem + 1 - ipISO  )
					{
						/* do not include highest levels since maser can occur due to topoff,
						 * and pops are set to small number in this case */
						for( long ipHi=1; ipHi <iso_sp[ipISO][nelem].numLevels_local - iso_sp[ipISO][nelem].nCollapsed_local; ipHi++ )
						{
							for( long ipLo=0; ipLo < ipHi; ipLo++ )
							{
								if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() <= 0 )
									continue;

								ASSERT( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() > iso_ctrl.SmallA );

								if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().PopOpc() > SMALLFLOAT &&
									/* test that have not overrun optical depth scale */
									( (iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauTot() - 
										iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauIn()) > SMALLFLOAT ) &&
									iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pesc() > FLT_EPSILON*100. )
								{
									RadPres1 =  PressureRadiationLine( iso_sp[ipISO][nelem].trans(ipHi,ipLo), GetDopplerWidth(dense.AtomicWeight[nelem]) );
									
									if( RadPres1 > pmx )
									{
										ipLineTypePradMax = 2;
										pmx = RadPres1;
										ip4 = ipISO;
										ip3 = nelem;
										ip2 = ipHi;
										ipLinePradMax = ipLo;
									}
									Piso_seq[ipISO] += RadPres1;
									{
										/* option to print particulars of some line when called */
										enum {DEBUG_LOC=false};
										if( DEBUG_LOC && ipISO==ipH_LIKE && ipLo==3 && ipHi==5 && nzone > 260 )
										{
											fprintf(ioQQQ,
												"DEBUG prad1 \t%.2f\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t\n",
												fnzone,
												RadPres1,
												iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().PopOpc(),
												iso_sp[ipISO][nelem].st[ipLo].Pop(),
												iso_sp[ipISO][nelem].st[ipHi].Pop(),
												iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pesc());
										}
									}
								}
							}
						}
					}
				}
				ASSERT( Piso_seq[ipISO] >= 0. );
			}
			pressure.pres_radiation_lines_curr += Piso_seq[ipISO];
		}
		{
			/* option to print particulars of some line when called */
			enum {DEBUG_LOC=false};
#			if 0
			if( DEBUG_LOC /*&& iteration > 1*/ && nzone > 260 )
			{
				fprintf(ioQQQ,
					"DEBUG prad2 \t%li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
					nzone,
					pmx,
					iso_sp[ipISO][ip3].trans(ipLinePradMax,ip2).Emis().PopOpc(),
					iso_sp[ipISO][ip3].st[0].Pop(),
					iso_sp[ipISO][ip3].st[2].Pop(),
					iso_sp[ipISO][ip3].st[3].Pop(),
					iso_sp[ipISO][ip3].st[4].Pop(),
					iso_sp[ipISO][ip3].st[5].Pop(),
					iso_sp[ipISO][ip3].st[6].Pop());
			}
#			endif
			if( DEBUG_LOC /*&& iteration > 1 && nzone > 150 */)
			{
				fprintf(ioQQQ,
					"DEBUG prad3\t%.2e\t%li\t%li\t%li\t%li\t%.2e\t%.2e\t%.2e\n",
					pmx,
					ip4,
					ip3,
					ip2,
					ipLinePradMax,
					iso_sp[ip4][ip3].trans(ip2,ipLinePradMax).Emis().PopOpc(),
					iso_sp[ip4][ip3].st[ip2].Pop(),
					1.-iso_sp[ip4][ip3].trans(ip2,ipLinePradMax).Emis().Pesc() );
			}
		}

		if( lgMustEvaluateTrivial || PLevel2 > TrivialLineRadiationPressure )
		{
			/* level 2 lines */
			PLevel2 = 0.;
			for( long i=0; i < nWindLine; i++ )
			{
				if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
				{
					if( (*TauLine2[i].Hi()).Pop() > 1e-30 )
					{
						RadPres1 = PressureRadiationLine( TauLine2[i], GetDopplerWidth(dense.AtomicWeight[(*TauLine2[i].Hi()).nelem()-1]) );

						PLevel2 += RadPres1;
						if( RadPres1 > pmx )
						{
							ipLineTypePradMax = 4;
							pmx = RadPres1;
							ipLinePradMax = i;
						}
					}
				}
			}
			ASSERT( PLevel2 >= 0. );
		}
		pressure.pres_radiation_lines_curr += PLevel2;

		/* fine structure lines */
		if( lgMustEvaluateTrivial || PHfs > TrivialLineRadiationPressure )
		{
			PHfs = 0.;
			for( size_t i=0; i < HFLines.size(); i++ )
			{
				if( (*HFLines[i].Hi()).Pop() > 1e-30 )
				{
					RadPres1 = PressureRadiationLine( HFLines[i], GetDopplerWidth(dense.AtomicWeight[(*HFLines[i].Hi()).nelem()-1]) );

					PHfs += RadPres1;
					if( RadPres1 > pmx )
					{
						ipLineTypePradMax = 8;
						pmx = RadPres1;
						ipLinePradMax = i;
					}
				}
			}
			ASSERT( PHfs >= 0. );
		}
		pressure.pres_radiation_lines_curr += PHfs;

		/* radiation pressure due to H2 lines */
		if( lgMustEvaluateTrivial || P_H2 > TrivialLineRadiationPressure )
		{
			P_H2 = 0.;
			for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
			{
				P_H2 += (*diatom)->H2_RadPress();
				/* flag to remember H2 radiation pressure */
				if( P_H2 > pmx )
				{
					pmx = P_H2;
					ipLineTypePradMax = 9;
				}
				ASSERT( P_H2 >= 0. );
			}
		}
		pressure.pres_radiation_lines_curr += P_H2;

		/* do lines from third-party databases */
		if( lgMustEvaluateTrivial || P_dBase > TrivialLineRadiationPressure )
		{
			P_dBase = 0.;
			for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
			{
				if( dBaseSpecies[ipSpecies].lgActive )
				{
					realnum DopplerWidth = GetDopplerWidth( dBaseSpecies[ipSpecies].fmolweight );
					for (TransitionList::iterator tr=dBaseTrans[ipSpecies].begin(); 
						  tr != dBaseTrans[ipSpecies].end(); ++tr)
					{	
						int ipHi = (*tr).ipHi();
						if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
							continue;
						int ipLo = (*tr).ipLo();
						if( (*tr).ipCont() > 0 && (*(*tr).Hi()).Pop() > 1e-30 )
						{
							RadPres1 =  PressureRadiationLine( *tr, DopplerWidth );
							
							if( RadPres1 > pmx )
							{
								ipLineTypePradMax = 10;
								pmx = RadPres1;
								ip3 = ipSpecies;
								ip2 = ipHi;
								ipLinePradMax = ipLo;
							}
							P_dBase += RadPres1;
						}
					}
				}
			}
			ASSERT( P_dBase >= 0. );
		}
		pressure.pres_radiation_lines_curr += P_dBase;

	}
	else if( !pressure.lgLineRadPresOn || !rfield.lgDoLineTrans )
	{
		/* case where radiation pressure is not turned on */
		ipLinePradMax = -1;
		ipLineTypePradMax = 0;
	}

	fixit("all other line stacks need to be included here");
	// can we just sweep over line stack?  Is that ready yet?

	/* the ratio of radiation to total (gas + continuum + etc) pressure */
	pressure.pbeta = (realnum)(pressure.pres_radiation_lines_curr/SDIV(pressure.PresTotlCurr));

	/* this would be a major logical error */
	if( pressure.pres_radiation_lines_curr < 0. )
	{
		fprintf(ioQQQ,
			"PresTotCurrent: negative pressure, constituents follow %e %e %e %e \n",
		Piso_seq[ipH_LIKE],
		Piso_seq[ipHE_LIKE],
		PLevel2,
		PCO);
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* following test will never succeed, here to trick lint, ipLineTypePradMax is only used
	 * when needed for debug */
	if( trace.lgTrace && ipLineTypePradMax <0 )
	{
		fprintf(ioQQQ,
			"     PresTotCurrent, pressure.pbeta = %e, ipLineTypePradMax%li ipLinePradMax=%li \n", 
			pressure.pbeta,ipLineTypePradMax,ipLinePradMax );
	}

	/* this is the total internal energy of the gas, the amount of energy needed
	 * to produce the current level populations, relative to ground,
	 * only used for energy terms in advection */
	phycon.EnergyExcitation = 0.;
#if 0
	fixit(); /* the quantities phycon.EnergyExcitation, phycon.EnergyIonization, phycon.EnergyBinding
		  * are not used anywhere, except in print statements... */
	broken(); /* the code below contains serious bugs! It is supposed to implement loops
		   * over quantum states in order to evaluate the potential energy stored in
		   * excited states of all atoms, ions, and molecules. The code below however
		   * often implements loops over all combinations of upper and lower levels!
		   * This code needs to be rewritten once quantumstates are fully implemented. */
	ipLo = 0;
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.IonHigh[nelem] == nelem + 1 - ipISO )
			{
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				for( long ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
				{
					phycon.EnergyExcitation += 
						iso_sp[ipISO][nelem].st[ipHi].Pop() * 
						iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyErg() *
						/* last term is option to turn off energy term, no advection hlike, he-like */
						dynamics.lgISO[ipISO];
				}
			}
		}
	}

	if( dynamics.lgISO[ipH_LIKE] )
		/* internal energy of H2 */
		phycon.EnergyExcitation += H2_InterEnergy();

	/* this is option to turn off energy effects of advection, no advection metals */
	if( dynamics.lgMETALS )
	{
		for( long i=0; i < nWindLine; i++ )
		{
			if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
			{
				phycon.EnergyExcitation += 
					(*TauLine2[i].Hi()).Pop()* TauLine2[i].EnergyErg;
			}
		}
		for( long i=0; i < nHFLines; i++ )
		{
			phycon.EnergyExcitation += 
				(*HFLines[i].Hi()).Pop()* HFLines[i].EnergyErg;
		}

		/* internal energy of large FeII atom */
	}
#endif

	/* ==================================================================
	 * end internal energy of atoms and molecules */

	/* evaluate some parameters to do with magnetic field */
	Magnetic_evaluate();

	/*lint -e644 Piso_seq not init */
	if( trace.lgTrace && (pressure.PresTotlCurr > SMALLFLOAT) )
	{
		fprintf(ioQQQ,
			"     PresTotCurrent zn %.2f Ptot:%.5e Pgas:%.3e Prad:%.3e Pmag:%.3e Pram:%.3e "
			"gas parts; H:%.2e He:%.2e L2:%.2e CO:%.2e fs%.2e h2%.2e turb%.2e\n",
			fnzone,
			pressure.PresTotlCurr, 
			pressure.PresGasCurr/pressure.PresTotlCurr,
			pressure.pres_radiation_lines_curr*pressure.lgPres_radiation_ON/pressure.PresTotlCurr,
			magnetic.pressure/pressure.PresTotlCurr,
			pressure.PresRamCurr*pressure.lgPres_ram_ON/pressure.PresTotlCurr,
			Piso_seq[ipH_LIKE]/pressure.PresTotlCurr,
			Piso_seq[ipHE_LIKE]/pressure.PresTotlCurr,
			PLevel2/pressure.PresTotlCurr,
			PCO/pressure.PresTotlCurr,
			PHfs/pressure.PresTotlCurr,
			P_H2/pressure.PresTotlCurr,
			pressure.PresTurbCurr*DoppVel.lgTurb_pressure/pressure.PresTotlCurr);
		/*lint +e644 Piso_seq not initialized */
	}

	/* Conserved quantity in steady-state energy equation for dynamics:
	 * thermal energy, since recombination is treated as cooling
	 * (i.e. is loss of electron k.e. to emitted photon), so don't
	 * include
	 * ...phycon.EnergyExcitation + phycon.EnergyIonization + phycon.EnergyBinding
	 * */

	/* constant density case is also hypersonic case */
	if( dynamics.lgTimeDependentStatic )
	{
		/* this branch is time dependent single-zone */
		/* \todo	1	this has 3/2 on the PresGasCurr while true dynamics case below
		 * has 5/2 - so this is not really the enthalpy density - better
		 * would be to always use this term and add the extra PresGasCurr
		 * when the enthalpy is actually needed */
		phycon.EnthalpyDensity =  
			0.5*POW2(wind.windv)*dense.xMassDensity +	/* KE */
			3./2.*pressure.PresGasCurr +				/* thermal plus PdV work */
			magnetic.EnthalpyDensity;					/* magnetic terms */
			//pressure.RhoGravity * (radius.Radius/(radius.Radius+radius.drad)); /* gravity */
	}
	else
	{
		/* this branch either advective or constant pressure */
		/*fprintf(ioQQQ,"DEBUG enthalpy HIT2\n");*/
		/* this is usual dynamics case, or time-varying case where pressure
		 * is kept constant and PdV work is added to the cell */
		phycon.EnthalpyDensity =  
			0.5*POW2(wind.windv)*dense.xMassDensity +	/* KE */
			5./2.*pressure.PresGasCurr +				/* thermal plus PdV work */
			magnetic.EnthalpyDensity;					/* magnetic terms */
			//pressure.RhoGravity * (radius.Radius/(radius.Radius+radius.drad)); /* gravity */
	}

	/* stop model from running away on first iteration, when line optical
	 * depths are not defined correctly anyway.
	 * if( iter.le.itermx .and. RadPres.ge.GasPres ) then
	 * >>chng 97 jul 23, only cap radiation pressure on first iteration */
	if( iteration <= 1 && pressure.pres_radiation_lines_curr >= pressure.PresGasCurr )
	{
		/* stop RadPres from exceeding the gas pressure on first iteration */
		pressure.pres_radiation_lines_curr = 
			MIN2(pressure.pres_radiation_lines_curr,pressure.PresGasCurr);
		pressure.lgPradCap = true;
	}

	/* remember the globally most important line, in the entire model 
	 * test on nzone so we only do test after solution is stable */
	if( pressure.pbeta > pressure.RadBetaMax && nzone )
	{
		pressure.RadBetaMax = pressure.pbeta;
		pressure.ipPradMax_line = ipLinePradMax;
		pressure.ipPradMax_nzone = nzone;
		if( ipLineTypePradMax == 2 )
		{
			/* hydrogenic */
			pressure.chLineRadPres = "ISO   ";
			ASSERT( ip4 < NISO && ip3<LIMELM );
			ASSERT( ipLinePradMax>=0 && ip2>=0 && ip3>=0 && ip4>=0 );
			pressure.chLineRadPres += chLineLbl(iso_sp[ip4][ip3].trans(ip2,ipLinePradMax));
			{
				/* option to print particulars of some line when called */
				enum {DEBUG_LOC=false};
				/*lint -e644 Piso_seq not initialized */
				/* trace serious constituents in radiation pressure */
				if( DEBUG_LOC  )
				{
					fprintf(ioQQQ,"DEBUG iso prad\t%li\t%li\tISO,nelem=\t%li\t%li\tnu, no=\t%li\t%li\t%.4e\t%.4e\t%e\t%e\t%e\n",
					iteration, nzone, 
					ip4,ip3,ip2,ipLinePradMax,
					iso_sp[ip4][ip3].trans(ip2,ipLinePradMax).Emis().TauIn(),
					iso_sp[ip4][ip3].trans(ip2,ipLinePradMax).Emis().TauTot(),
					iso_sp[ip4][ip3].trans(ip2,ipLinePradMax).Emis().Pesc(),
					iso_sp[ip4][ip3].trans(ip2,ipLinePradMax).Emis().Pelec_esc(),
					iso_sp[ip4][ip3].trans(ip2,ipLinePradMax).Emis().Pdest());
					if( ip2==5 && ipLinePradMax==4 )
					{
						double width;
						fprintf(ioQQQ,"hit it\n");
						width = RT_LineWidth(iso_sp[ip4][ip3].trans(ip2,ipLinePradMax),GetDopplerWidth(dense.AtomicWeight[ip3]));
						fprintf(ioQQQ,"DEBUG width %e\n", width);
					}
				}
			}
		}
		else if( ipLineTypePradMax == 4 )
		{
			/* level 2 lines */
			ASSERT( ipLinePradMax>=0  );
			pressure.chLineRadPres = "Level2 ";
			pressure.chLineRadPres += chLineLbl(TauLine2[ipLinePradMax]);
		}
		else if( ipLineTypePradMax == 5 )
		{
			cdEXIT( EXIT_FAILURE );
		}
		else if( ipLineTypePradMax == 6 )
		{
			cdEXIT( EXIT_FAILURE );
		}
		else if( ipLineTypePradMax == 7 )
		{
			/* FeII lines */
			pressure.chLineRadPres = "Fe II ";
		}
		else if( ipLineTypePradMax == 8 )
		{
			/* hyperfine struct lines */
			pressure.chLineRadPres = "hyperf ";
			ASSERT( ipLinePradMax>=0  );
			pressure.chLineRadPres += chLineLbl(HFLines[ipLinePradMax]);
		}
		else if( ipLineTypePradMax == 9 )
		{
			/* large H2 lines */
			pressure.chLineRadPres = "large H2 ";
		}
		else if( ipLineTypePradMax == 10 )
		{
			/* database lines */
			pressure.chLineRadPres = "dBaseLin " ;
			pressure.chLineRadPres += dBaseSpecies[ip3].chLabel ;
		}
		else
		{
			fprintf(ioQQQ," PresTotCurrent ipLineTypePradMax set to %li, this is impossible.\n", ipLineTypePradMax);
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}
	}

	if( trace.lgTrace && pressure.pbeta > 0.5 )
	{
		fprintf(ioQQQ,
			"     PresTotCurrent Pbeta:%.2f due to %s\n",
			pressure.pbeta ,
		   pressure.chLineRadPres.c_str()
			);
	}

	/* >>chng 02 jun 27 - add in magnetic pressure */
	/* start to bring total pressure together */
	TotalPressure_v = pressure.PresGasCurr;

	/* these flags are set false by default since constant density is default,
	 * set true for constant pressure or dynamics */
	TotalPressure_v += pressure.PresRamCurr * pressure.lgPres_ram_ON;

	/* magnetic pressure, evaluated in magnetic.c - this can be negative for an ordered field! 
	 * option is on by default, turned off with constant density, or constant gas pressure, cases */
	/** \todo	0	code has variable magnetic energydensity and pressure, which are equal,
	 * as they must be - del one or the other */
	TotalPressure_v += magnetic.pressure * pressure.lgPres_magnetic_ON;

	/* turbulent pressure
	 * >>chng 06 mar 24, include this by default */
	TotalPressure_v += pressure.PresTurbCurr * DoppVel.lgTurb_pressure;

	/* radiation pressure only included in total eqn of state when this flag is
	 * true, set with constant pressure command */
	/* option to add in internal line radiation pressure */
	TotalPressure_v += pressure.pres_radiation_lines_curr * pressure.lgPres_radiation_ON;

	{
		enum{DEBUG_LOC=false};
		if( DEBUG_LOC && pressure.PresTotlCurr > SMALLFLOAT /*&& iteration > 1*/ )
		{
			fprintf(ioQQQ,"pressureee%li\t%.4e\t%.4e\t%.4e\t%.3f\t%.3f\t%.3f\n", 
				nzone,
				pressure.PresTotlError,
			    pressure.PresTotlCurr, 
				TotalPressure_v ,
				pressure.PresGasCurr/pressure.PresTotlCurr,
				pressure.pres_radiation_lines_curr/pressure.PresTotlCurr ,
				pressure.PresRamCurr/pressure.PresTotlCurr
				);
		}
	}

	if( TotalPressure_v< 0. )
	{
		ASSERT( magnetic.pressure < 0. );

		/* negative pressure due to ordered field overwhelms total pressure - punt */
		fprintf(ioQQQ," The negative pressure due to ordered magnetic field overwhelms the total outward pressure - please reconsider the geometry & field.\n");
		cdEXIT(EXIT_FAILURE);
	}

	ASSERT( TotalPressure_v > 0. );

	/* remember highest pressure anywhere */
	pressure.PresMax = MAX2(pressure.PresMax,(realnum)TotalPressure_v);

	/* this is what we came for - set the current pressure */
	pressure.PresTotlCurr = TotalPressure_v;

	return;
}
