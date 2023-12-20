/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolEvaluate main routine to call others, to evaluate total cooling */
#include "cddefines.h"
#include "cool_eval.h"
#include "taulines.h"
#include "wind.h"
#include "coolheavy.h"
#include "radius.h"
#include "conv.h"
#include "h2.h"
#include "rt.h"
#include "opacity.h"
#include "ionbal.h"
#include "trace.h"
#include "dynamics.h"
#include "grainvar.h"
#include "atmdat.h"
#include "atoms.h"
#include "called.h"
#include "hmi.h"
#include "numderiv.h"
#include "magnetic.h"
#include "phycon.h"
#include "hyperfine.h"
#include "iso.h"
#include "thermal.h"
#include "cooling.h"
#include "pressure.h"
#include "mole.h"
#include "rfield.h"
#include "doppvel.h"
#include "freebound.h"
#include "dense.h"
#include "atmdat_gaunt.h"
#include "iterations.h"
#include "vectorize.h"


/* H2 cooling according to
 * >>refer h2 cool Glover & Abel 2008, MNRAS, 388, 1627; Sections 2.3.1-2.3.6 */
STATIC double CoolH2_GA08 ( double Temp );

/*fndneg search cooling array to find negative values */
STATIC void fndneg(void);

/*fndstr search cooling stack to find strongest values */
STATIC void fndstr(double tot, double dc);

STATIC void CoolHyperfine (void);

STATIC double eeBremsCooling(double Te,    // electron temperature in K
			     double eden); // electron density in cm^-3

/* set true to debug derivative of heating and cooling */
static const bool PRT_DERIV = false;

void CoolEvaluate(double *tot)
{
	static long int nhit = 0, 
	  nzSave=0;

	static double oltcool=0., 
	  oldtemp=0.;

	DEBUG_ENTRY( "CoolEvaluate()" );

	/* returns tot, the total cooling,
	 * and dc, the derivative of the cooling */

	if( trace.lgTrace )
		fprintf( ioQQQ, "   COOLR TE:%.4e zone %li %li Cool:%.4e Heat:%.4e eden:%.4e edenTrue:%.4e\n", 
		phycon.te, 
		nzone, conv.nPres2Ioniz ,
		thermal.ctot , thermal.htot,dense.eden,dense.EdenTrue );

	/* must call TempChange since ionization has changed, there are some
	 * terms that affect collision rates (H0 term in electron collision) */
	TempChange(phycon.te , false);

	/* now zero out the cooling stack */
	CoolZero();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  0 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);
	if( gv.lgGrainPhysicsOn )
	{
		/* grain heating and cooling */
		/* grain recombination cooling, evaluated elsewhere
		* can either heat or cool the gas, do cooling here */
		CoolAdd("dust",0,MAX2(0.,gv.GasCoolColl));

		/* grain cooling proportional to temperature ^3/2 */
		thermal.dCooldT += MAX2(0.,gv.GasCoolColl)*3./(2.*phycon.te);

		/* these are the various heat agents from grains */
		/* options to force gas heating or cooling by grains to zero - for tests only ! */
		if( gv.lgDustOn() && gv.lgDHetOn )
		{
			/* rate dust heats gas by photoelectric effect */
			thermal.setHeating(0,13, gv.GasHeatPhotoEl);

			/* if grains hotter than gas then collisions with gas act
			* to heat the gas, add this in here
			* a symmetric statement appears in COOLR, where cooling is added on */
			thermal.setHeating(0,14, MAX2(0.,-gv.GasCoolColl) );

			/* this is gas heating due to thermionic emissions */
			thermal.setHeating(0,25, gv.GasHeatTherm );
		}
		else
		{
			thermal.setHeating(0,13,0.);
			thermal.setHeating(0,14,0.);
			thermal.setHeating(0,25,0.);
		}
	}
	else if( gv.lgBakesPAH_heat )
	{
		/* >>chng 06 jul 21, option to include Bakes PAH hack with grain physics off,
		 * needed to test dynamics models */
		thermal.setHeating(0,13,gv.GasHeatPhotoEl);
	}

	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  1 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* molecular molecules molecule cooling */
	if( mole_global.lgNoMole )
	{
		/* this branch - do not include molecules */
		hmi.hmicol = 0.;
		CoolHeavy.brems_cool_hminus = 0.;
		/* line cooling within simple H2 molecule - zero when big used */
		CoolHeavy.h2line = 0.;
		/*  H + H+ => H2+ cooling */
		CoolHeavy.H2PlsCool = 0.;
		CoolHeavy.HD = 0.;

		/* thermal.heating(0,8) is heating due to collisions within X of H2 */
		thermal.setHeating(0,8,0.);
		/* thermal.heating(0,15) is H minus heating*/
		thermal.setHeating(0,15,0.);
		/* thermal.heating(0,16) is H2+ heating */
		thermal.setHeating(0,16,0.);
		hmi.HeatH2Dish_used = 0.;
		hmi.HeatH2Dexc_used = 0.;
		hmi.deriv_HeatH2Dexc_used = 0.;
	}

	else
	{
		/* save various molecular heating/cooling agent */
		thermal.setHeating(0,15, hmi.hmihet);
		thermal.setHeating(0,16, hmi.h2plus_heat);
		/* now get heating from H2 molecule, either simple or from big one */
		if( h2.lgEnabled  && hmi.lgH2_Thermal_BigH2 )
		{
			if( h2.lgEvaluated )
			{
				/* these are explicitly from big H2 molecule,
				 * first is heating due to radiative pump of excited states, followed by
				 * radiative decay into continuum of X, followed by dissociation of molecule
				 * with kinetic energy, typically 0.25 - 0.5 eV per event */
				hmi.HeatH2Dish_used = h2.HeatDiss;
				hmi.HeatH2Dexc_used = h2.HeatDexc;
				if (0)
					fprintf(ioQQQ,"DEBUG big %.2f\t%.5e\t%.2e\t%.2e\t%.2e\n", 
							  fnzone , phycon.te, hmi.HeatH2Dexc_used,
							  hmi.H2_total, dense.gas_phase[ipHYDROGEN] );
				/* negative sign because right term is really deriv of heating,
				 * but will be used below as deriv of cooling */
				hmi.deriv_HeatH2Dexc_used = -h2.HeatDexc_deriv;
			}
			else
			{
				hmi.HeatH2Dish_used = 0;
				hmi.HeatH2Dexc_used = 0;
				hmi.deriv_HeatH2Dexc_used = 0;
			}
		}

		else if( hmi.chH2_small_model_type == 'T' )
		{
			/* TH85 dissociation heating */
			/* these come from approximations in TH85, see comments above */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_TH85;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_TH85;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_TH85;
		}
		else if( hmi.chH2_small_model_type == 'H' )
		{
			/* Burton et al. 1990 */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_BHT90;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_BHT90;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_BHT90;
		}
		else if( hmi.chH2_small_model_type == 'B' )
		{
			/* Bertoldi & Draine */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_BD96;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_BD96;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_BD96;
		}
		else if( hmi.chH2_small_model_type == 'E' )
		{
			/* this is the default when small H2 used */
			hmi.HeatH2Dish_used = hmi.HeatH2Dish_ELWERT;
			hmi.HeatH2Dexc_used = hmi.HeatH2Dexc_ELWERT;
			hmi.deriv_HeatH2Dexc_used = hmi.deriv_HeatH2Dexc_ELWERT;
		}
		else
			TotalInsanity();

		/* heating due to photodissociation heating */
		thermal.setHeating(0,17,hmi.HeatH2Dish_used);

		/* heating due to continuum photodissociation */
		thermal.setHeating(0,28,0.);
		{
			double heat = 0.;
			for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
			{
				if( (*diatom)->lgEnabled && mole_global.lgStancil )
				{
					heat += (*diatom)->Cont_Diss_Heat_Rate();
				}
			}
			thermal.setHeating(0,28,MAX2( 0., heat ));
			CoolAdd("H2cD",0,MAX2(0.,-heat));
		}

		/* heating (usually cooling in big H2) due to collisions within X */
		/* add to heating is net heating is positive */
		thermal.setHeating(0,8, MAX2(0.,hmi.HeatH2Dexc_used) +
				// HD heating, cooling is CoolHeavy.HD
				MAX2(0.,hd.HeatDexc) + MAX2(0.,hd.HeatDiss));

		/* add to cooling if net heating is negative */
		CoolAdd("H2cX",0,MAX2(0.,-hmi.HeatH2Dexc_used));
		/*fprintf(ioQQQ,"DEBUG coolh2\t%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
			fnzone, phycon.te, dense.eden, hmi.H2_total, thermal.ctot, -hmi.HeatH2Dexc_used );*/
		/* add to net derivative */
		/*thermal.dCooldT += MAX2(0.,-hmi.HeatH2Dexc_used)* ( 30172. * thermal.tsq1 - thermal.halfte );*/
		/* >>chng 04 jan 25, check sign to prevent cooling from entering here,
		 * also enter neg sign since going into cooling stack (bug), in heatsum
		 * same term adds to deriv of heating */
		if( hmi.HeatH2Dexc_used < 0. )
			thermal.dCooldT -= hmi.deriv_HeatH2Dexc_used;

		/*  H + H+ => H2+ cooling */
		CoolHeavy.H2PlsCool = (realnum)(MAX2((2.325*phycon.te-1875.)*1e-20,0.)*
		  dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipHYDROGEN][1]*1.66e-11);

		if( h2.lgEnabled )
		{
			/* this is simplified approximation to H2 rotation cooling,
			 * big molecule does this far better */
			CoolHeavy.h2line = 0.;
		}
		else
		{
			CoolHeavy.h2line = CoolH2_GA08( phycon.te ) * hmi.H2_total;

			{
				enum {DEBUG_LOC=false};
				if( DEBUG_LOC && nzone>187&& iteration > 1/**/)
				{
					fprintf(ioQQQ,"h2coolbug\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
						phycon.te, 
						CoolHeavy.h2line, 
						hmi.H2_total, 
						findspecieslocal("H-")->den, 
						hmi.HMinus_photo_rate );
				}
			}

		}

		if( hd.lgEnabled )
		{
			// heating was thermal.setHeating(0,8 above, with H2
			CoolHeavy.HD = -MIN2(0.,hd.HeatDexc) - MIN2(0.,hd.HeatDiss);
		}
		else
		{
			/* >>chng 22 aug 13, use Flower et al. (2000, MNRAS, 314, 753)
			 * HD cooling function
			 * http://ccp7.dur.ac.uk/cooling_by_HD/node3.html */
#if 0
			factor = 0.25*pow2( log10(double(dense.gas_phase[ipHYDROGEN])))+
					(0.283978*log10(double(dense.gas_phase[ipHYDROGEN]))-1.27333)*
					(sin(2.03657*log10(phycon.te)+4.63258))-2.08189*
					log10(double(dense.gas_phase[ipHYDROGEN]))+4.66288;

			CoolHeavy.HD = hmi.HD_total*exp10((0.5*log10(double(dense.gas_phase[ipHYDROGEN]))) +
					(-26.2982*pow(log10(phycon.te), -0.215807)) - sqrt(factor));
#endif

			double aa=-26.2982, bb=-0.215807, omeg=2.03657, phi=4.63258,
					c1=0.283978, c2=-1.27333, d1=-2.08189, d2=4.66288;

			// the sum of hydrogen nuclei in all forms
			double y = log10(dense.gas_phase[ipHYDROGEN]);
			double x = phycon.alogte;

			double w = 0.5 * y + aa * pow(x,bb)
				- sqrt( 0.25 * y*y
				+ (c1*y+c2) * sin(omeg*x+phi) + (d1*y+d2) );

			CoolHeavy.HD = hmi.HD_total * exp10(w);

		}
	}

	fixit("test and enable chemical heating by default");
#if 0
	double chemical_heating = mole.chem_heat();	
	thermal.setHeating(0,29, MAX2(0.,chemical_heating) );
	/* add to cooling if net heating is negative */
	CoolAdd("Chem",0,MAX2(0.,-chemical_heating));
#endif

	/* cooling due to charge transfer ionization / recombination */
	CoolAdd("CT C" , 0. , thermal.char_tran_cool );

	/*  H- FB; H + e -> H- + hnu */
	/*  H- FF is in with H ff */
	CoolAdd("H-fb",0,hmi.hmicol);

	/* >>chng 96 nov 15, fac of 2 in deriv to help convergence in very dense
	 * models where H- is important, this takes change in eden into
	 * partial account */
	thermal.dCooldT += 2.*hmi.hmicol*phycon.teinv;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  2 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	CoolAdd("H2ln",0,CoolHeavy.h2line);
	/* >>chng 00 oct 21, added coef of 3.5, sign had been wrong */
	/*thermal.dCooldT += CoolHeavy.h2line*phycon.teinv;*/
	/* >>chng 03 mar 17, change 3.5 to 15 as per behavior in primal.in */
	/*thermal.dCooldT += 3.5*CoolHeavy.h2line*phycon.teinv;*/
	/* >>chng 03 may 18, from 15 to 30 as per behavior in primal.in - overshoots happen */
	/*thermal.dCooldT += 15.*CoolHeavy.h2line*phycon.teinv;*/
	/*>>chng 03 oct 03, from 30 to 3, no more overshoots in primalin */
	/*thermal.dCooldT += 30.*CoolHeavy.h2line*phycon.teinv;*/
	thermal.dCooldT += 3.0*CoolHeavy.h2line*phycon.teinv;

	{
		/* problems with H2 cooling */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC /*&& nzone>300 && iteration > 1*/)
		{
			fprintf(ioQQQ,"CoolEvaluate debuggg\t%.2f\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",
			fnzone, 
			phycon.te, 
			hmi.H2_total , 
			CoolHeavy.h2line,
			findspecieslocal("H-")->den , 
			dense.eden);
		}
	}

	CoolAdd("HDro",0,CoolHeavy.HD);
	thermal.dCooldT += CoolHeavy.HD*phycon.teinv;

	CoolAdd("H2+ ",0,CoolHeavy.H2PlsCool);
	thermal.dCooldT += CoolHeavy.H2PlsCool*phycon.teinv;

	/* heating due to three-body, will be incremented in iso_cool*/
	thermal.setHeating(0,3,0.);
	/* heating due to hydrogen lines */
	thermal.setHeating(0,23,0.);
	/* heating due to photoionization of all excited states of hydrogen species */
	thermal.setHeating(0,1,0.);

	/* isoelectronic species cooling, mainly lines, and level ionization */
	for( long int ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long int nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always call iso_cool since we must zero those variables
			* that would have been set had the species been present */
			iso_cool( ipISO , nelem );
		}
	}

	/* free-free free free brems emission for all ions */
	/* highest frequency where we have non-zero Boltzmann factors */
	long limit = min( rfield.ipMaxBolt, rfield.nflux );

	CoolHeavy.brems_cool_h = 0.;
	CoolHeavy.brems_cool_hminus = 0.;
	CoolHeavy.brems_cool_he = 0.;
	CoolHeavy.brems_cool_metals = 0.;
	CoolHeavy.brems_heat_total = 0.;

	if( CoolHeavy.lgFreeOn )
	{
		ASSERT(rfield.ipEnergyBremsThin < rfield.nflux_with_check);
		ASSERT(limit < rfield.nflux_with_check);

		t_brems_den sum;
		t_gaunt::Inst().brems_sum_ions(sum);
		double bfac = dense.eden * FREE_FREE_EMIS / phycon.sqrte * EN1RYD;

		/* do hydrogen first, before main loop since want to break out as separate
		 * coolant, and what to add on H- brems */
		CoolHeavy.brems_cool_h = sum.den_Hp*bfac*t_gaunt::Inst().brems_cool( 1, phycon.te );
		CoolHeavy.brems_cool_hminus = sum.den_Hm*bfac*t_gaunt::Inst().brems_cool( -1, phycon.te );

		/* now do helium, both He+ and He++ */
		CoolHeavy.brems_cool_he = sum.den_Hep*bfac*t_gaunt::Inst().brems_cool( 1, phycon.te ) +
			sum.den_Hepp*bfac*t_gaunt::Inst().brems_cool( 2, phycon.te );

		/* heavy elements */
		CoolHeavy.brems_cool_metals = 0.;
		for( long ion=1; ion < LIMELM+1; ++ion )
			if( sum.den_ion[ion] > 0. )
				CoolHeavy.brems_cool_metals +=
					sum.den_ion[ion]*bfac*t_gaunt::Inst().brems_cool( ion, phycon.te );

		/* ipEnergyBremsThin is index to energy where gas becomes optically thin to brems,
		 * so this loop is over optically thin frequencies 
		 * do not include optically thick part as net emission since self absorbed */
		CoolHeavy.brems_heat_total = 0.;
		for( long int i=rfield.ipEnergyBremsThin; i < limit; i++ )
		{
			/* the total heating due to bremsstrahlung */
			CoolHeavy.brems_heat_total +=
				opac.FreeFreeOpacity[i]*rfield.flux[0][i]*rfield.anu(i)*EN1RYD;
		}

		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nzone>60 /*&& iteration > 1*/)
		{
			double sumfield = 0., sumtot=0., sum1=0., sum2=0.;
			for( long int i=rfield.ipEnergyBremsThin; i<limit;  i++ )
			{
				sumtot += opac.FreeFreeOpacity[i]*rfield.flux[0][i]*rfield.anu(i);
				sumfield += rfield.flux[0][i]*rfield.anu(i);
				sum1 += opac.FreeFreeOpacity[i]*rfield.flux[0][i]*rfield.anu(i);
				sum2 += opac.FreeFreeOpacity[i]*rfield.flux[0][i];
			}
			fprintf(ioQQQ,"DEBUG brems heat\t%.2f\t%.3e\t%.3e\t%.3e\t%e\t%.3e\t%.3e\n",
				fnzone,
				CoolHeavy.brems_heat_total,
				sumtot/SDIV(sumfield) ,
				sum1/SDIV(sum2),
				phycon.te , 
				t_gaunt::Inst().gauntff(1,phycon.te,rfield.anu(1218)),
				opac.FreeFreeOpacity[1218]);
		}
	}

	/* these two terms are both large, nearly canceling, near lte */
	CoolHeavy.brems_cool_net = 
		CoolHeavy.brems_cool_h + 
		CoolHeavy.brems_cool_he + 
		CoolHeavy.brems_cool_hminus + 
		CoolHeavy.brems_cool_metals - 
		CoolHeavy.brems_heat_total;
	/*fprintf(ioQQQ,"DEBUG brems\t%.2f\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
		fnzone,
		phycon.te,
		CoolHeavy.brems_cool_net,
		CoolHeavy.brems_cool_h ,
		CoolHeavy.brems_cool_he ,
		CoolHeavy.brems_cool_hminus,
		CoolHeavy.brems_cool_metals ,
		CoolHeavy.brems_heat_total);*/

	/* net free free brems cooling, count as cooling if positive */
	CoolAdd( "FF c" , 0, MAX2(0.,CoolHeavy.brems_cool_net) );

	/* now stuff into heating array if negative */
	thermal.setHeating(0,11, MAX2(0.,-CoolHeavy.brems_cool_net) );

	/* >>chng 96 oct 30, from HFFNet to just FreeFreeCool,
	 * since HeatSum picks up CoolHeavy.brems_heat_total */
	thermal.dCooldT += CoolHeavy.brems_cool_h*thermal.halfte;
	thermal.dCooldT += CoolHeavy.brems_cool_he*thermal.halfte;
	thermal.dCooldT += CoolHeavy.brems_cool_hminus*thermal.halfte;
	thermal.dCooldT += CoolHeavy.brems_cool_metals*thermal.halfte;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  3 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* >>chng 02 jun 21, net cooling already includes this */
	/* end of brems cooling */

	/* heavy element recombination cooling, do not count hydrogenic since
	 * already done above, also helium singlets have been done */
	/* >>chng 02 jul 21, put in charge dependent rec term */
	CoolHeavy.heavfb = 0.;
	for( long int nelem=ipLITHIUM; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* note that detailed iso seq atoms are done in iso_cool */
			long limit_lo = MAX2( 1 , dense.IonLow[nelem] );
			long limit_hi = MIN2( nelem-NISO+1, dense.IonHigh[nelem] );
			for( long int ion=limit_lo; ion<=limit_hi; ++ion )
			{
				/* factor of 0.9 is roughly correct for nebular conditions, see
				 * >>refer	H	rec cooling	LaMothe, J., & Ferland, G.J., 2001, PASP, 113, 165 */
				/* note that ionbal.RR_rate_coef_used is the rate coef, cm3 s-1, needs eden */
				/* >>chng 02 nov 07, move rec arrays around, this now has ONLY rad rec,
				 * previously had included grain rec and three body */
				/* recombination cooling for iso-seq atoms are done in iso_cool */
				double one = dense.xIonDense[nelem][ion] * ionbal.RR_rate_coef_used[nelem][ion-1]*
					dense.eden * phycon.te * BOLTZMANN;
				/*fprintf(ioQQQ,"debugggfb\t%li\t%li\t%.3e\t%.3e\t%.3e\n", nelem, ion, one, 
					dense.xIonDense[nelem][ion] , ionbal.RR_rate_coef_used[nelem][ion]);*/
				CoolHeavy.heavfb += one;
				thermal.elementcool[nelem] += one;
			}
		}
	}

	/*fprintf(ioQQQ,"debuggg hvFB\t%i\t%.2f\t%.2e\t%.2e\n",iteration, fnzone,CoolHeavy.heavfb, dense.eden);*/

	CoolAdd("hvFB",0,CoolHeavy.heavfb);
	thermal.dCooldT += CoolHeavy.heavfb*.113*phycon.teinv;
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  4 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	CoolHeavy.eebrm = eeBremsCooling(phycon.te,dense.eden);

	/* >>chng 97 mar 12, added deriv */
	thermal.dCooldT += CoolHeavy.eebrm*thermal.halfte;
	CoolAdd("eeff",0,CoolHeavy.eebrm);

	/* add advective heating and cooling */
	/* this is cooling due to loss of matter from this region */
	CoolAdd("adve",0,dynamics.Cool() );
	/* >>chng02 dec 04, rm factor of 8 in front of dCooldT */
	thermal.dCooldT += dynamics.dCooldT();
	/* local heating due to matter moving into this location */
	thermal.setHeating(1,5, dynamics.Heat() );
	thermal.dHeatdT += dynamics.dHeatdT;

	/* total Compton cooling */
	CoolHeavy.tccool = rfield.cmcool*phycon.te;
	CoolAdd("Comp",0,CoolHeavy.tccool);
	thermal.dCooldT += rfield.cmcool;

	/* option for "extra" cooling, expressed as power-law in temperature, these
	 * are set with the CEXTRA command */
	if( thermal.lgCExtraOn )
	{
		CoolHeavy.cextxx = 
			(realnum)(thermal.CoolExtra*pow(phycon.te/1e4,(double)thermal.cextpw));
	}
	else
	{
		CoolHeavy.cextxx = 0.;
	}
	CoolAdd("Extr",0,CoolHeavy.cextxx);

	/* cooling due to wind expansion, only for winds expansion cooling */
	if( wind.lgBallistic() )
	{
		realnum dDensityDT = -(realnum)(wind.AccelTotalOutward/wind.windv + 2.*wind.windv/
		  radius.Radius);
		CoolHeavy.expans = -2.5*pressure.PresGasCurr*dDensityDT;
	}
	else if( dynamics.lgTimeDependentStatic && 
			! dynamics.isInitialRelaxIteration( iteration ) )
	{
		realnum dens = scalingDensity();
		realnum dDensityDT = 
			(realnum)((dens-dynamics.Upstream_density)/
			(dynamics.timestep*0.5*(dens+dynamics.Upstream_density)));
		// pdV work term
		CoolHeavy.expans = -pressure.PresGasCurr*dDensityDT;
	}
	else
	{
		CoolHeavy.expans = 0.;
	}
	CoolAdd("Expn",0,CoolHeavy.expans);
	thermal.dCooldT += CoolHeavy.expans/phycon.te;

	/* cyclotron cooling */
	/* coef is 4/3 /8pi /c * vtherm(elec) */
	CoolHeavy.cyntrn = 4.5433e-25f*magnetic.pressure*PI8*dense.eden*phycon.te;
	CoolAdd("Cycl",0,CoolHeavy.cyntrn);
	thermal.dCooldT += CoolHeavy.cyntrn/phycon.te;

	/* heavy element collisional ionization
	 * derivative should be zero since increased coll ion rate
	 * decreases neutral fraction by proportional amount */
	CoolAdd("Hvin",0,CoolHeavy.colmet);


	double xIonDenseSave[LIMELM][LIMELM+1];
	if( atmdat.lgChiantiOn ||atmdat.lgStoutOn)
	{
		for( int nelem=0; nelem < LIMELM; nelem++ )
		{		
			for( int ion=0; ion<=nelem+1; ++ion )
			{
				xIonDenseSave[nelem][ion] = dense.xIonDense[nelem][ion];
				// zero abundance of species if we are using Chianti for this ion
				if( dense.lgIonChiantiOn[nelem][ion] || dense.lgIonStoutOn[nelem][ion] )
					dense.xIonDense[nelem][ion] = 0.;
			}
		}
	}

	(*TauDummy).Zero();
	(*(*TauDummy).Hi()).g() = 0.;
	(*(*TauDummy).Lo()).g() = 0.;
	(*(*TauDummy).Hi()).IonStg() = 0;
	(*(*TauDummy).Lo()).IonStg() = 0;
	(*(*TauDummy).Hi()).nelem() = 0;
	(*(*TauDummy).Lo()).nelem() = 0;
	(*TauDummy).Emis().Aul() = 0.;
	(*TauDummy).EnergyWN() = 0.;
	(*TauDummy).WLAng() = 0.;

	// reset abundances to original values, may have been set zero to protect against old cloudy lines
	if( atmdat.lgChiantiOn || atmdat.lgStoutOn)
	{
		// this clause, first reset abundances set to zero when Chianti included
		for( int nelem=0; nelem < LIMELM; nelem++ )
		{
			for( int ion=0; ion<=nelem+1; ++ion )
			{
				dense.xIonDense[nelem][ion] = xIonDenseSave[nelem][ion];
			}
		}
	}

	/* opacity project lines Dima Verner added with g-bar approximation */
	CoolDima();

	/* do external database lines */
	dBaseTrim();
	dBaseUpdateCollCoeffs();
 	dBase_solve();

	/* Hyperfine line cooling */
	CoolHyperfine();

	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  6 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* Print number of levels for each species */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			static bool lgMustPrintHeader=true;
			if( lgMustPrintHeader )
			{
				lgMustPrintHeader = false;
				printf("DEBUG Levels\t%s",dBaseSpecies[0].chLabel );
				for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
				{
					printf("\t%s",dBaseSpecies[ipSpecies].chLabel );
				}
				printf("\n" );
				printf("DEBUG Max\t%li" ,dBaseSpecies[0].numLevels_max );
				for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
				{
					printf( "\t%li" ,dBaseSpecies[ipSpecies].numLevels_max );
				}
				printf("\n");
			}
			printf("DEBUG Local\t%li" ,dBaseSpecies[0].numLevels_local );
			for( long ipSpecies=1; ipSpecies<nSpecies; ipSpecies++ )
			{
				printf("\t%li" ,dBaseSpecies[ipSpecies].numLevels_local );
			}
			printf("\n");
		}
	}

	/* now add up all the coolants */
	CoolSum(tot);
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT 14 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);

	/* negative cooling */
	if( *tot <= 0. )
	{
		fprintf( ioQQQ, " COOLR; cooling is <=0, this is impossible.\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* bad derivative */
	if( thermal.dCooldT == 0. )
	{
		fprintf( ioQQQ, " COOLR; cooling slope <=0, this is impossible.\n" );
		if( *tot > 0. && dense.gas_phase[ipHYDROGEN] < 1e-4 )
		{
			fprintf( ioQQQ, " Probably due to very low density.\n" );
		}
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( trace.lgTrace )
	{
		fndstr(*tot,thermal.dCooldT);
	}

	/* lgTSetOn true for constant temperature model */
	if( (((((!thermal.lgTemperatureConstant) && *tot < 0.) && called.lgTalk) && 
	  !conv.lgSearch) && thermal.lgCNegChk) && nzone > 0 )
	{
		fprintf( ioQQQ, 
			" NOTE Negative cooling, zone %4ld, =%10.2e coola=%10.2e CHION=%10.2e Te=%10.2e\n", 
		  nzone, 
		  *tot, 
		  iso_sp[ipH_LIKE][ipHYDROGEN].cLya_cool, 
		  iso_sp[ipH_LIKE][ipHYDROGEN].coll_ion, 
		  phycon.te );
		fndneg();
	}

	/* possibility of getting empirical cooling derivative
	 * normally false, set true with 'set numerical derivatives' command */
	if( NumDeriv.lgNumDeriv )
	{
		if( ((nzone > 2 && nzone == nzSave) && ! fp_equal( oldtemp, phycon.te )) && nhit > 4 )
		{
			/* hnit is number of tries on this zone - use to stop numerical problems
			 * do not evaluate numerical deriv until well into solution */
			thermal.dCooldT = (oltcool - *tot)/(oldtemp - phycon.te);
		}
		if( nzone != nzSave )
			nhit = 0;

		nzSave = nzone;
		nhit += 1;
		oltcool = *tot;
		oldtemp = phycon.te;
	}
	return;
}


/*******************************************************************************
 *
 * 		GLOVER & ABEL (2008) H2 LOW-DENSITY COOLING
 *
 ******************************************************************************/

/* Implementation of log10 of H2 cooling rate, given by
 * >>refer eqn(27)	Glover & Abel 2008, MNRAS, 388, 1627 */
STATIC double ga08_sum ( double Temp, double *coeff, int ncoeff )
{
	DEBUG_ENTRY( "ga08_sum()" );
	double logT = log10( Temp * 0.001 );
	double logTn = 1.0;
	double sum = 0.;
	for( int i = 0; i < ncoeff; i++ )
	{
		sum += coeff[i] * logTn;
		logTn *= logT;
	}

	return sum;
}

/* Cooling rate of ortho-H2 due to collisions with atomic H at T <= 100 K
 * >>refer eqn(28)	Glover & Abel 2008, MNRAS, 388, 1627 */
STATIC double ga08_oH2_H_b100 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_H_b100()" );
	const double DE = 852.5;	// energy in K
	return 5.09e-27 * pow( Temp * 0.001, 0.5) * dsexp( DE / Temp );
}

/* Cooling rate of ortho-H2 due to collisions with atomic H at 100 K < T < 1000K
 * >>refer eqn(27) & Table 1 second column	Glover & Abel 2008, MNRAS, 388, 1627 */
STATIC double ga08_oH2_H_b1000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_H_b1000()" );
	double coeff_oH2_H_cool[] = { -24.330855, 4.4404496, -4.0460989, -1.1390725, 9.8094223, 8.6273872 };
	int ncoeff = sizeof( coeff_oH2_H_cool ) / sizeof( coeff_oH2_H_cool[0] );
	return ga08_sum( Temp, coeff_oH2_H_cool, ncoeff );
}

/* Cooling rate of ortho-H2 due to collisions with atomic H at 1000 K <= T < 6000K
 * >>refer eqn(27) & Table 1 third column	Glover & Abel 2008, MNRAS, 388, 1627 */
STATIC double ga08_oH2_H_b6000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_H_b6000()" );
	/* >>refer h2 cool H Glover & Abel 2008, MNRAS, 388, 1627; Table 1, third column */
	double coeff_oH2_H_warm[] = { -24.329086, 4.6105087, -3.9505350, 12.363818, -32.403165, 48.853562, -38.542008, 12.066770 };
	int ncoeff = sizeof( coeff_oH2_H_warm ) / sizeof( coeff_oH2_H_warm[0] );
	return ga08_sum( Temp, coeff_oH2_H_warm, ncoeff );
}

/* Cooling rate of ortho-H2 due to collisions with atomic H at T > 6000K
 * Artificial tail to give the cooling rate a finite value across temperature range */
/* Modify the ortho-H2 cooling rate due to atomic H collisions for a small
 * temperature range above 100 K to obtain a continuous function across 100 K.
 * Returns log10 of correction. */
STATIC double ga08_oH2_H_stitch_100 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_H_stitch_100()" );
	const double Tbreak = 100.;
	static double Temp_c = -1.;
	if( Temp_c < 0. )
	{
		// _bb -> below break; _ab -> above break
		double crate_bb = log10( ga08_oH2_H_b100( Tbreak ) );
		double crate_ab = ga08_oH2_H_b1000( Tbreak );	// this is log10( rate )
		Temp_c = Tbreak / ( 1. + (crate_ab - crate_bb) / LOG10_E );
		if( 0 )
			printf("oH2 H, at 100K:\t %f (%.10e vs %.10e)\n", Temp_c, crate_bb, crate_ab);
	}

	double corr = 0.;
	if( Temp > Tbreak && Temp <= Temp_c )
		corr = LOG10_E * ( 1. - Temp / Temp_c );
	return corr;
}

/* Modify the ortho-H2 cooling rate due to atomic H collisions for a small
 * temperature range below 1000 K to obtain a continuous function across 1000 K.
 * Returns log10 of correction. */
STATIC double ga08_oH2_H_stitch_1000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_H_stitch_1000()" );
	const double Tbreak = 1000.;
	static double Temp_c = -1.;
	if( Temp_c < 0. )
	{
		// _bb -> below break; _ab -> above break
		double crate_bb = ga08_oH2_H_b1000( Tbreak );
		double crate_ab = ga08_oH2_H_b6000( Tbreak );	// this is log10( rate )
		Temp_c = Tbreak / ( 1. + (crate_ab - crate_bb) / LOG10_E );
		if( 0 )
			printf("oH2 H, at 1000K:\t %f (%.10e vs %.10e)\n", Temp_c, crate_bb, crate_ab);
	}

	double corr = 0.;
	if( Temp >= Temp_c && Temp_c < Tbreak )
		corr = LOG10_E * ( Temp / Temp_c - 1. );
	return corr;
}

/* Cooling rate of ortho-H2 due to collisions with atomic H at any T.
 * Extrapolate below 100 K. */
STATIC double ga08_oH2_H ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_H()" );
	static double cool_above_6000 = -1.;
	if( cool_above_6000 < 0. )
	{
		cool_above_6000 = pow( 10., ga08_oH2_H_b6000( 6000.0 ) );
	}

	double cool = 0.;
	if( Temp <= 100. )
	{
		cool = ga08_oH2_H_b100( Temp );
	}
	else if( Temp < 1000. )
	{
		cool = ga08_oH2_H_b1000( Temp );
		cool += ga08_oH2_H_stitch_100( Temp );
		cool += ga08_oH2_H_stitch_1000( Temp );
		cool = pow(10., cool);
	}
	else if( Temp < 6000. )
	{
		cool = pow( 10., ga08_oH2_H_b6000( Temp ) );
	}
	else
	{
		cool = cool_above_6000;
	}
	return cool;
}



/* Cooling rate of para-H2 due to collisions with atomic H at T <= 100 K
 * >>refer eqn(29)	Glover & Abel 2008, MNRAS, 388, 1627 */
STATIC double ga08_pH2_H_b100 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_H_b100()" );
	const double E20 = 509.85;	// energy in K
	return 8.16e-26 * pow(Temp * 0.001, 0.5) * dsexp( E20 / Temp );
}

/* Cooling rate of para-H2 due to collisions with atomic H at 100 K < T < 1000K
 * >>refer eqn(27) & Table 2 second column	Glover & Abel 2008, MNRAS, 388, 1627 */
STATIC double ga08_pH2_H_b1000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_H_b1000()" );
	double coeff_pH2_H_cool[] = { -24.216387, 3.3237480, -11.642384, -35.553366, -35.105689, -10.922078 };
	int ncoeff = sizeof( coeff_pH2_H_cool ) / sizeof( coeff_pH2_H_cool[0] );
	return ga08_sum( Temp, coeff_pH2_H_cool, ncoeff );
}

/* Cooling rate of para-H2 due to collisions with atomic H at 1000 K <= T < 6000K
 * >>refer eqn(27) & Table 2 third column	Glover & Abel 2008, MNRAS, 388, 1627 */
STATIC double ga08_pH2_H_b6000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_H_b6000()" );
	double coeff_pH2_H_warm[] = { -24.216387, 4.2046488, -1.3155285, -1.6552763, 4.1780102, -0.56949697, -3.3824407, 1.0904027 };
	int ncoeff = sizeof( coeff_pH2_H_warm ) / sizeof( coeff_pH2_H_warm[0] );
	return ga08_sum( Temp, coeff_pH2_H_warm, ncoeff );
}

/* Modify the para-H2 cooling rate due to atomic H collisions for a small
 * temperature range above 100 K to obtain a continuous function across 100 K.
 * Returns log10 of correction. */
STATIC double ga08_pH2_H_stitch_100 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_H_stitch_100()" );
	const double Tbreak = 100.;
	static double Temp_c = -1.;
	if( Temp_c < 0. )
	{
		// _bb -> below break; _ab -> above break
		double crate_bb = log10( ga08_pH2_H_b100( Tbreak ) );
		double crate_ab = ga08_pH2_H_b1000( Tbreak );	// this is log10( rate )
		Temp_c = Tbreak / ( 1. + (crate_ab - crate_bb) / LOG10_E );
		if( 0 )
			printf("pH2 H, at 100K:\t %f\n", Temp_c);
	}

	double corr = 0.;
	if( Temp > Tbreak && Temp <= Temp_c )
		corr = LOG10_E * ( 1. - Temp / Temp_c );
	return corr;
}

/* Cooling rate of para-H2 due to collisions with atomic H at any T.
 * Extrapolate below 100 K. */
STATIC double ga08_pH2_H ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_H()" );
	static double cool_above_6000 = -1.;
	if( cool_above_6000 < 0. )
	{
		cool_above_6000 = pow( 10., ga08_pH2_H_b6000( 6000.0 ) );
	}

	double cool = 0.;
	if( Temp <= 100. )
	{
		cool = ga08_pH2_H_b100( Temp );
	}
	else
	if( Temp < 1000. )
	{
		cool = ga08_pH2_H_b1000( Temp );
		cool +=	ga08_pH2_H_stitch_100( Temp );
		cool = pow(10., cool);
	}
	else if( Temp < 6000. )
	{
		cool = pow( 10., ga08_pH2_H_b6000( Temp ) );
	}
	else
	{
		cool = cool_above_6000;
	}
	return cool;
}



/* Cooling rate of para-H2 due to collisions with para-H2 at 100 K <= T <= 6000 K
 * >>refer h2 cool para-para	Glover & Abel 2008, MNRAS, 388, 1627; Table 3, second column */
STATIC double ga08_pH2_pH2_a100_b6000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_pH2_a100_b6000()" );
	double coeff_pH2_pH2[] = { -23.889798, 1.8550774, -0.55593388, 0.28429361, -0.20581113, 0.13112378 };
	int ncoeff = sizeof( coeff_pH2_pH2 ) / sizeof( coeff_pH2_pH2[0] );
	return ga08_sum( Temp, coeff_pH2_pH2, ncoeff );
}

/* Cooling rate of para-H2 due to collisions with para-H2 at any T.
 * Extrapolate below 100 K. */
STATIC double ga08_pH2_pH2 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_pH2()" );
	static double cool_above_6000 = -1.;
	if( cool_above_6000 < 0. )
	{
		cool_above_6000 = pow( 10., ga08_pH2_pH2_a100_b6000( 6000. ) );
	}

	double cool = 0.;
	if( Temp <= 6000. )
		cool = pow( 10., ga08_pH2_pH2_a100_b6000( Temp ) );
	else
		cool = cool_above_6000; 
	return cool;
}

/* Cooling rate of para-H2 due to collisions with para-H2 at 100 K <= T <= 6000 K
 * >>refer h2 cool para-ortho	Glover & Abel 2008, MNRAS, 388, 1627; Table 3, third column */
STATIC double ga08_pH2_oH2_a100_b6000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_oH2_a100_b6000()" );
	double coeff_pH2_oH2[] = { -23.748534, 1.76676480, -0.58634325, 0.31074159, -0.17455629, 0.18530758 };
	int ncoeff = sizeof( coeff_pH2_oH2 ) / sizeof( coeff_pH2_oH2[0] );
	return ga08_sum( Temp, coeff_pH2_oH2, ncoeff );
}

/* Cooling rate of para-H2 due to collisions with ortho-H2 at any T.
 * Extrapolate below 100 K. */
STATIC double ga08_pH2_oH2 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_oH2()" );
	static double cool_above_6000 = -1.;
	if( cool_above_6000 < 0. )
	{
		cool_above_6000 = pow( 10., ga08_pH2_oH2_a100_b6000( 6000.0 ) );
	}

	double cool = 0.;
	if( Temp <= 6000. )
		cool = pow( 10., ga08_pH2_oH2_a100_b6000( Temp ) );
	else
		cool = cool_above_6000;
	return cool;
}



/* Cooling rate of ortho-H2 due to collisions with para-H2 at 100 K <= T <= 6000 K
 * >>refer h2 cool ortho-para	Glover & Abel 2008, MNRAS, 388, 1627; Table 4, second column */
STATIC double ga08_oH2_pH2_a100_b6000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_pH2_a100_b6000()" );
	double coeff_oH2_pH2[] = { -24.126177, 2.3258217, -1.0082491, 0.54823768, -0.33679759, 0.20771406 };
	int ncoeff = sizeof( coeff_oH2_pH2 ) / sizeof( coeff_oH2_pH2[0] );
	return ga08_sum( Temp, coeff_oH2_pH2, ncoeff );
}

/* Cooling rate of ortho-H2 due to collisions with para-H2 at any T.
 * Extrapolate below 100 K. */
STATIC double ga08_oH2_pH2 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_pH2()" );
	static double cool_above_6000 = -1.;
	if( cool_above_6000 < 0. )
	{
		cool_above_6000 = pow( 10., ga08_oH2_pH2_a100_b6000( 6000. ) );
	}
	double cool = 0.;
	if( Temp <= 6000. )
		cool =  pow( 10., ga08_oH2_pH2_a100_b6000( Temp ) );
	else
		cool = cool_above_6000;
	return cool;
}



/* Cooling rate of ortho-H2 due to collisions with ortho-H2 at 100 K <= T <= 6000 K
 * >>refer h2 cool ortho-ortho	Glover & Abel 2008, MNRAS, 388, 1627; Table 4, third column */
STATIC double ga08_oH2_oH2_a100_b6000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_oH2_a100_b6000()" );
	double coeff_oH2_oH2[] = { -24.020047, 2.2687566, -1.0200304, 0.83561432, -0.40772247, 0.096025713 };
	int ncoeff = sizeof( coeff_oH2_oH2 ) / sizeof( coeff_oH2_oH2[0] );
	return ga08_sum( Temp, coeff_oH2_oH2, ncoeff );
}

/* Cooling rate of ortho-H2 due to collisions with ortho-H2 at any T.
 * Extrapolate below 100 K. */
STATIC double ga08_oH2_oH2 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_oH2()" );
	static double cool_above_6000 = -1.;
	if( cool_above_6000 < 0. )
	{
		cool_above_6000 = pow( 10., ga08_oH2_oH2_a100_b6000( 6000.0 ) );
	}

	double cool = 0.;
	if( Temp <= 6000. )
		cool = pow( 10., ga08_oH2_oH2_a100_b6000( Temp ) );
	else
		cool = cool_above_6000;
	return cool;
}



/* Cooling rate of para-H2 due to collisions with He at T <= 6000 K
 * >>refer h2 cool He-para	Glover & Abel 2008, MNRAS, 388, 1627; Table 5, second column */
STATIC double ga08_pH2_He_b6000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_He_b6000()" );
	double coeff_pH2_He[] = { -23.489029 , 1.8210825 , -0.59110559 , 0.42280623 , -0.30171138 , 0.12872839 };
	int ncoeff = sizeof( coeff_pH2_He ) / sizeof( coeff_pH2_He[0] );
	return ga08_sum( Temp, coeff_pH2_He, ncoeff );
}

/* Cooling rate of para-H2 due to collisions with He at any T */
STATIC double ga08_pH2_He ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_He()" );
	static double cool_above_6000 = -1.;
	if( cool_above_6000 < 0. )
	{
		cool_above_6000 = pow( 10., ga08_pH2_He_b6000( 6000.0 ) );
	}

	double cool = 0.;
	if( Temp <= 6000. )
		cool = pow( 10., ga08_pH2_He_b6000( Temp ) );
	else
		cool = cool_above_6000;
	return	cool;
}


/* Cooling rate of ortho-H2 due to collisions with He at T <= 6000 K
 * >>refer h2 cool He-ortho	Glover & Abel 2008, MNRAS, 388, 1627; Table 5, third column */
STATIC double ga08_oH2_He_b6000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_He_b6000()" );
	double coeff_oH2_He[] = { -23.7749 , 2.40654 , -1.23449 , 0.739874 , -0.258940 , 0.120573 };
	int ncoeff = sizeof( coeff_oH2_He ) / sizeof( coeff_oH2_He[0] );
	return ga08_sum( Temp, coeff_oH2_He, ncoeff );
}


/* Cooling rate of ortho-H2 due to collisions with He at any T */
STATIC double ga08_oH2_He ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_He()" );
	static double cool_above_6000 = -1.;
	if( cool_above_6000 < 0. )
	{
		cool_above_6000 = pow( 10., ga08_oH2_He_b6000( 6000.0 ) );
	}

	double cool = 0.;
	if( Temp <= 6000. )
		cool = pow( 10., ga08_oH2_He_b6000( Temp ) );
	else
		cool = cool_above_6000;
	return	cool;
}



/* Cooling rate of para-H2 due to collisions with He at 10 K <= T <= 10000 K
 * >>refer h2 cool p-para	Glover & Abel 2008, MNRAS, 388, 1627; Table 6, second column */
STATIC double ga08_pH2_p_a10_b1e4 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_p_a10_b1e4()" );
	double coeff_pH2_p[] = { -21.757160 , 1.3998367 , -0.37209530 , 0.061554519 , -0.37238286 , 0.23314157 };
	int ncoeff = sizeof( coeff_pH2_p ) / sizeof( coeff_pH2_p[0] );
	return ga08_sum( Temp, coeff_pH2_p, ncoeff );
}

/* Cooling rate of para-H2 due to collisions with p at any T */
STATIC double ga08_pH2_p ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_p()" );
	static double cool_above_10000 = -1.;
	if( cool_above_10000 < 0. )
	{
		cool_above_10000 = pow( 10., ga08_pH2_p_a10_b1e4( 1e4 ) );
	}

	double cool = 0.;
	if( Temp <= 1e4 )
		cool = pow( 10., ga08_pH2_p_a10_b1e4( Temp ) );
	else
		cool = cool_above_10000;
	return	cool;
}


/* Cooling rate of ortho-H2 due to collisions with p at 10 K <= T <= 10000 K
 * >>refer h2 cool p-ortho	Glover & Abel 2008, MNRAS, 388, 1627; Table 6, third column */
STATIC double ga08_oH2_p_a10_b1e4 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_p_a10_b1e4()" );
	double coeff_oH2_p[] = { -21.706641 , 1.3901283 , -0.34993699 , 0.075402398 , -0.23170723 , 0.068938876 };
	int ncoeff = sizeof( coeff_oH2_p ) / sizeof( coeff_oH2_p[0] );
	return ga08_sum( Temp, coeff_oH2_p, ncoeff );
}

/* Cooling rate of ortho-H2 due to collisions with p at any T */
STATIC double ga08_oH2_p ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_p()" );
	static double cool_above_10000 = -1.;
	if( cool_above_10000 < 0. )
	{
		cool_above_10000 = pow( 10., ga08_oH2_p_a10_b1e4( 1e4 ) );
	}

	double cool = 0.;
	if( Temp <= 1e4 )
		cool = pow( 10., ga08_oH2_p_a10_b1e4( Temp ) );
	else
		cool = cool_above_10000;
	return	cool;
}



/* Cooling rate of para-H2 due to collisions with e at 10 K <= T <= 1000 K
 * >>refer h2 cool e-para	Glover & Abel 2008, MNRAS, 388, 1627; Table 7, second column */
STATIC double ga08_pH2_e_a10_b1000 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_e_a10_b1000()" );
	double coeff_pH2_e_cool[] = { -22.817869 , 0.95653474 , 0.79283462 , 0.56811779 , 0.27895033 , 0.056049813 };
	int ncoeff = sizeof( coeff_pH2_e_cool ) / sizeof( coeff_pH2_e_cool[0] );
	return ga08_sum( Temp, coeff_pH2_e_cool, ncoeff );
}

/* Cooling rate of para-H2 due to collisions with e at 1000 K < T <= 10000 K
 * >>refer h2 cool e-para	Glover & Abel 2008, MNRAS, 388, 1627; Table 7, third column */
STATIC double ga08_pH2_e_a1000_b1e4 ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_e_a1000_b1e4()" );
	double coeff_pH2_e_warm[] = { -22.817869 , 0.66916141 , 7.1191428 , -11.176835 , 7.0467275 , -1.6471816 };
	int ncoeff = sizeof( coeff_pH2_e_warm ) / sizeof( coeff_pH2_e_warm[0] );
	return ga08_sum( Temp, coeff_pH2_e_warm, ncoeff );
}

/* Cooling rate of para-H2 due to collisions with e at any T */
STATIC double ga08_pH2_e ( double Temp )
{
	DEBUG_ENTRY( "ga08_pH2_e()" );
	static double cool_above_10000 = -1.;
	if( cool_above_10000 < 0. )
	{
		cool_above_10000 = pow( 10., ga08_pH2_e_a1000_b1e4( 1e4 ) );
	}

	const double E20 = 509.85;	// energy in K
	double cool = 0.;

	if( Temp <= 1000. )
		cool = pow( 10., ga08_pH2_e_a10_b1000( Temp ) );
	else if( Temp <= 1e4 )
		cool = pow( 10., ga08_pH2_e_a1000_b1e4( Temp ) );
	else
		cool = cool_above_10000;

	return dsexp( E20 / Temp ) * cool; 
}



/* Cooling rate of ortho-H2 due to collisions with e at 10 K <= T <= 10000 K
 * >>refer h2 cool e-para	Glover & Abel 2008, MNRAS, 388, 1627; Table 7, second column */
STATIC double ga08_oH2_e_a10_b1e4 ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_e_a10_b1e4()" );
	double coeff_oH2_e[] = { -21.703215 , 0.76059565 , 0.50644890 , 0.050371349 , -0.10372467 , -0.035709409 };
	int ncoeff = sizeof( coeff_oH2_e ) / sizeof( coeff_oH2_e[0] );
	return ga08_sum( Temp, coeff_oH2_e, ncoeff );
}

/* Cooling rate of ortho-H2 due to collisions with e at any T */
STATIC double ga08_oH2_e ( double Temp )
{
	DEBUG_ENTRY( "ga08_oH2_e()" );
	static double cool_above_10000 = -1.;
	if( cool_above_10000 < 0. )
	{
		cool_above_10000 = pow( 10., ga08_oH2_e_a10_b1e4( 1e4 ) );
	}

	double cool = 0.;
	if( Temp <= 1e4 )
		cool = pow( 10., ga08_oH2_e_a10_b1e4( Temp ) );
	else
		cool = cool_above_10000;

	const double E31 = 845.;
	return dsexp( E31 / Temp ) * cool;
}



/* H2 cooling per H2 molecule according to
 * >>refer h2 cool Glover & Abel 2008, MNRAS, 388, 1627; Sections 2.3.1-2.3.6 */
STATIC double CoolH2_GA08 ( double Temp )
{
	DEBUG_ENTRY( "CoolH2_GA08()" );

	double x_para = safe_div( h2.para_density, double( dense.gas_phase[ ipHYDROGEN ] ) );
	double x_ortho = safe_div( h2.ortho_density, double( dense.gas_phase[ ipHYDROGEN ] ) );
	double f_para = safe_div( h2.para_density, h2.para_density + h2.ortho_density );
	double f_ortho = 1. - f_para;

	/* Collisions with atomic H */
	double cool_H = dense.xIonDense[ipHYDROGEN][0] *
		( f_ortho * ga08_oH2_H( Temp )
		+ f_para * ga08_pH2_H( Temp ) );

	/* Collisions with H2 */
	double cool_H2 = hmi.H2_total *
		( x_para * x_para * ga08_pH2_pH2( Temp )
		+ x_para * x_ortho * ga08_pH2_oH2( Temp )
		+ x_ortho * x_para * ga08_oH2_pH2( Temp )
		+ x_ortho * x_ortho * ga08_oH2_oH2( Temp ) );

	/* Collisions with He */
	double cool_He = dense.xIonDense[ipHELIUM][0] *
		( f_para * ga08_pH2_He( Temp )
		+ f_ortho * ga08_oH2_He( Temp ) );

	/* Collisions with protons */
	double cool_p = dense.xIonDense[ipHYDROGEN][1] *
		( f_para * ga08_pH2_p( Temp )
		+ f_ortho * ga08_oH2_p( Temp ) );

	/* non-equilibrium cooling */
	//	/* >>refer h2 cool eq(35) J=1<->0	Glover & Abel 2008, MNRAS, 388, 1627; Table 6, third column */
	//	const double E10 = 170.5;
	//	cool_rate_per_H2 += 4.76e-24 * ( 9.*dsexp( E10 / Temp ) * x_para - x_ortho ) * dense.xIonDense[ipHYDROGEN][1];

	/* Collisions with electrons */
	double cool_e = dense.eden *
		( f_para * ga08_pH2_e( Temp )
		+ f_ortho * ga08_oH2_e( Temp ) );

	double cooling_low = cool_H + cool_H2 + cool_He + cool_p + cool_e;
	double cooling_high = h2.interpolate_LTE_Cooling( Temp );

	{
		enum { DEBUG_LOC = false };
		if( DEBUG_LOC && iterations.lgLastIt )
		{
			printf( "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
				Temp,
				dense.gas_phase[ ipHYDROGEN ],
				hmi.H2_total,
				cool_H,
				cool_H2,
				cool_He,
				cool_p ,
				cool_e ,
				cooling_low,
				cooling_high,
				cooling_high / (1. + cooling_high / cooling_low),
				cooling_high / (1. + cooling_high / cooling_low) * hmi.H2_total
				);
		}
	}

	// originally cooling_high / (1. + cooling_high / cooling_low)
	// See Eq. (39) therein
	//
	return safe_div( cooling_low * cooling_high, cooling_low + cooling_high );
}

STATIC void CoolHyperfine (void)
{
	DEBUG_ENTRY( "CoolHyperfine()" );
	static double	TeEvalCS = 0., TeEvalCS_21cm=0.;
	static double 	electron_rate_21cm,
		atomic_rate_21cm,
		proton_rate_21cm;


	for( int ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( int nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if ( ! dense.lgElmtOn[nelem] )
					continue;
			ionbal.ExcitationGround[nelem][nelem-ipISO] +=
					iso_sp[ipISO][nelem].fb[0].RateLevel2Cont;
		}
	}

	{
		enum {DEBUG_LOC=false};
		if ( DEBUG_LOC )
		{
			for (int nelem = 0; nelem < LIMELM; nelem++)
			{
				fprintf(ioQQQ, "ionbal.ExcitationGround: nelem = %d", nelem);
				for (int ion = 0; ion < nelem+1; ion++)
				{
					fprintf(ioQQQ,"\t%.6e", ionbal.ExcitationGround[nelem][ion]);
				}
				fprintf(ioQQQ, "\n");
			}
		}
	}

	/* evaluate H 21 cm spin changing collisions */
	if( !fp_equal(phycon.te,TeEvalCS_21cm) )
	{
		{
			/* this prints table of rates at points given in original data paper */
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
#				define N21CM_TE	16
				int n;
				double teval[N21CM_TE]={2.,5.,10.,20.,50.,100.,200.,500.,1000.,
					2000.,3000.,5000.,7000.,10000.,15000.,20000.};
				for( n = 0; n<N21CM_TE; ++n )
				{
					fprintf(
						ioQQQ,"DEBUG 21 cm deex Te=\t%.2e\tH0=\t%.2e\tp=\t%.2e\te=\t%.2e\n",
						teval[n], 
						H21cm_H_atom( teval[n] ),
						H21cm_proton( teval[n] ),
						H21cm_electron( teval[n] ) );
				}
				cdEXIT(EXIT_FAILURE);
#				undef N21CM_TE
			}
		}
		/*only evaluate T dependent part when Te changes, but update
		 * density part below since densities may constantly change */
		atomic_rate_21cm = H21cm_H_atom( phycon.te );
		proton_rate_21cm = H21cm_proton( phycon.te );
		electron_rate_21cm = H21cm_electron( phycon.te );
		TeEvalCS_21cm = phycon.te;
	}
	/* H 21 cm emission/population,
	* cs will be sum of e cs and H cs converted from rate */
	double cs = (electron_rate_21cm * dense.eden + 
		atomic_rate_21cm * dense.xIonDense[ipHYDROGEN][0] +
		proton_rate_21cm * dense.xIonDense[ipHYDROGEN][1] ) *
		3. / dense.cdsqte;
	PutCS(  cs , HFLines[0] );

	/* do level pops for H 21 cm which is a special case since Lya pumping in included */
	RT_line_one_escape( HFLines[0], true,0.f, GetDopplerWidth(dense.AtomicWeight[(*HFLines[0].Hi()).nelem()-1]) );
	H21_cm_pops();

	hyperfine.cooling_total = HFLines[0].Coll().cool();
	if( PRT_DERIV )
		fprintf(ioQQQ,"DEBUG dCdT  5 %.3e dHdT %.3e\n",thermal.dCooldT , thermal.dHeatdT);


	if( !fp_equal(phycon.te,TeEvalCS) )
	{
		/* do cross-sections for all except 21 cm  */
		for( size_t i=1; i < HFLines.size(); i++ )
		{
			cs = HyperfineCS( i );
			PutCS(  cs , HFLines[i] );
		}
		TeEvalCS = phycon.te;
	}


	/* now do level pops for all except 21 cm  */
	for( size_t i=1; i < HFLines.size(); i++ )
	{
		/* remember current gas-phase abundance of this isotope */
		double save = dense.xIonDense[(*HFLines[i].Hi()).nelem()-1][(*HFLines[i].Hi()).IonStg()-1];

		/* bail if no abundance */
		if( save<=0. ) 
			continue;

		RT_line_one_escape( HFLines[i], true,0.f, GetDopplerWidth(dense.AtomicWeight[(*HFLines[i].Hi()).nelem()-1]) );

		/* set gas-phase abundance to total times isotope ratio */
		dense.xIonDense[(*HFLines[i].Hi()).nelem()-1][(*HFLines[i].Hi()).IonStg()-1] *= 
			hyperfine.HFLabundance[i];

		/* use the collision strength generated above and find pops and cooling */
		atom_level2( HFLines[i], true );

		/* put the correct gas-phase abundance back in the array */
		dense.xIonDense[(*HFLines[i].Hi()).nelem()-1][(*HFLines[i].Hi()).IonStg()-1] = save;

		/* find total cooling due to hyperfine structure lines */
		hyperfine.cooling_total += HFLines[i].Coll().cool();
	}

	return;
}

static const double EPS = 0.01;

/*fndneg search cooling array to find negative values */
STATIC void fndneg(void)
{
	long int i;
	double trig;

	DEBUG_ENTRY( "fndneg()" );

	trig = fabs(thermal.htot*EPS);
	for( i=0; i < thermal.ncltot; i++ )
	{
		if( thermal.cooling[i] < 0. && fabs(thermal.cooling[i]) > trig )
		{
			fprintf( ioQQQ, " negative line=%s %.2f fraction of heating=%.3f\n", 
			  thermal.chClntLab[i], thermal.collam[i], thermal.cooling[i]/
			  thermal.htot );
		}

		if( thermal.heatnt[i] > trig )
		{
			fprintf( ioQQQ, " heating line=%s %.2f fraction of heating=%.3f\n", 
			  thermal.chClntLab[i], thermal.collam[i], thermal.heatnt[i]/
			  thermal.htot );
		}
	}
	return;
}

/*fndstr search cooling stack to find strongest values */
STATIC void fndstr(double tot, 
  double dc)
{
	char chStrngLab[NCOLNT_LAB_LEN+1];
	long int i;
	realnum wl;
	double str, 
	  strong;

	DEBUG_ENTRY( "fndstr()" );

	strong = 0.;
	wl = -FLT_MAX;
	for( i=0; i < thermal.ncltot; i++ )
	{
		if( fabs(thermal.cooling[i]) > strong )
		{
			/* this is the wavelength of the coolant, 0 for a continuum*/
			wl = thermal.collam[i];
			/* make sure labels are all valid*/
			/*>>chng 06 jun 06, bug fix, assert length was ==4, should be <=NCOLNT_LAB_LEN */
			ASSERT( strlen( thermal.chClntLab[i] ) <= NCOLNT_LAB_LEN );
			strcpy( chStrngLab, thermal.chClntLab[i] );
			strong = fabs(thermal.cooling[i]);
		}
	}

	str = strong;

	fprintf( ioQQQ, 
		"   fndstr cool: TE=%10.4e Ne=%10.4e C=%10.3e dC/dT=%10.2e ABS(%s %.1f)=%.2e nz=%ld\n", 
	  phycon.te, dense.eden, tot, dc, chStrngLab
	  , wl, str, nzone );

	/* option for extensive printout of lines */
	if( trace.lgCoolTr )
	{
		realnum ratio;

		/* flag all significant coolants, first zero out the array */
		coolpr(ioQQQ,thermal.chClntLab[0],1,0.,"ZERO");

		/* push all coolants onto the stack */
		for( i=0; i < thermal.ncltot; i++ )
		{
			/* usually positive, although can be neg for coolants that heats, 
			 * only do positive here */
			ratio = (realnum)(thermal.cooling[i]/thermal.ctot);
			if( ratio >= EPS )
			{
				/*>>chng 99 jan 27, only cal when ratio is significant */
				coolpr(ioQQQ,thermal.chClntLab[i],thermal.collam[i], ratio,"DOIT");
			}
		}

		/* complete the printout for positive coolants */
		coolpr(ioQQQ,"DONE",1,0.,"DONE");

		/* now do cooling that was counted as a heat source if significant */
		if( thermal.heating(0,22)/thermal.ctot > 0.05 )
		{
			fprintf( ioQQQ, 
				"     All coolant heat greater than%6.2f%% of the total will be printed.\n", 
			  EPS*100. );

			coolpr(ioQQQ,"ZERO",1,0.,"ZERO");
			for( i=0; i < thermal.ncltot; i++ )
			{
				ratio = (realnum)(thermal.heatnt[i]/thermal.ctot);
				if( fabs(ratio) >=EPS )
				{
					coolpr(ioQQQ,thermal.chClntLab[i],thermal.collam[i],
					  ratio,"DOIT");
				}
			}
			coolpr(ioQQQ,"DONE",1,0.,"DONE");
		}
	}
	return;
}


// Nozawa et al. (2009) Table 1: Coefficients a^I_ik (taken from Itoh et al. 2002a).
static const double aI[11][11] = { // a^I_ik = aI[i][k]
	{  3.15847e+00,-2.52430e+00, 4.04877e-01, 6.13466e-01, 6.28867e-01, 3.29441e-01,
	  -1.71486e-01,-3.68685e-01,-7.59200e-02, 1.60187e-01, 8.37729e-02 },
	{  2.46819e-02, 1.03924e-01, 1.98935e-01, 2.18843e-01, 1.20482e-01,-4.82390e-02,
	  -1.20811e-01,-4.46133e-04, 8.88749e-02, 2.50320e-02,-1.28900e-02 },
	{ -2.11118e-02,-8.53821e-02,-1.52444e-01,-1.45660e-01,-4.63705e-02, 8.16592e-02,
	   9.87296e-02,-3.24743e-02,-8.82637e-02,-7.52221e-03, 1.99419e-02 },
	{  1.24009e-02, 4.73623e-02, 7.51656e-02, 5.07201e-02,-2.25247e-02,-8.17151e-02,
	  -4.59297e-02, 5.05096e-02, 5.58818e-02,-9.11885e-03,-1.71348e-02 },
	{ -5.41633e-03,-1.91406e-02,-2.58034e-02,-2.23048e-03, 5.07325e-02, 5.94414e-02,
	  -2.11247e-02,-5.05387e-02, 9.20453e-03, 1.67321e-02,-3.47663e-03 },
	{  1.70070e-03, 5.39773e-03, 4.13361e-03,-1.14273e-02,-3.23280e-02,-2.19399e-02,
	   1.76310e-02, 2.23352e-02,-4.59817e-03,-8.24286e-03,-3.90032e-04 },
	{ -3.05111e-04,-7.26681e-04, 4.67015e-03, 1.24789e-02,-1.16976e-02,-1.13488e-02,
	   6.31446e-02, 1.33830e-02,-8.54735e-02,-6.47349e-03, 3.72266e-02 },
	{ -1.21721e-04,-7.47266e-04,-2.20675e-03,-2.74351e-03,-1.00402e-03,-2.38863e-03,
	  -2.28987e-03, 7.79323e-03, 7.98332e-03,-3.80435e-03,-4.25035e-03 },
	{  1.77611e-04, 8.73517e-04,-2.67582e-03,-4.57871e-03, 2.96622e-02, 1.89850e-02,
	  -8.84093e-02,-2.93629e-02, 1.02966e-01, 1.38957e-02,-4.22093e-02 },
	{ -2.05480e-05,-6.92284e-05, 2.95254e-05,-1.70374e-04,-5.43191e-04, 2.50978e-03,
	   4.45570e-03,-2.80083e-03,-5.68093e-03, 1.10618e-03, 2.33625e-03 },
	{ -3.58754e-05,-1.80305e-04, 1.40751e-03, 2.06757e-03,-1.23098e-02,-8.81767e-03,
	   3.46210e-02, 1.23727e-02,-4.04801e-02,-5.68689e-03, 1.66733e-02 }
};

// Nozawa et al. (2009) Table 2: Coefficients b^I_i (taken from Itoh et al. 2002a).
static const double bI[11] = { // b^I_i = bI[i]
	2.21564e+00, 1.83879e-01, -1.33575e-01, 5.89871e-02, -1.45904e-02, -7.10244e-04,
	2.80940e-03, -1.70485e-03, 5.26075e-04, 9.94159e-05, -1.06851e-04
};

// Nozawa et al. (2009) Table 3: Coefficients a^II_ij and b^II_ij.
static const double aII[9][3] = { // a^II_ij = aII[j][i]
	{ -3.7369800e+01,-9.3647000e+00, 9.2170000e-01 },
	{  3.8036590e+02, 9.5918600e+01,-1.3498800e+01 },
	{ -1.4898014e+03,-3.9701720e+02, 7.6453900e+01 },
	{  2.8614150e+03, 8.4293760e+02,-2.1783010e+02 },
	{ -2.3263704e+03,-9.0730760e+02, 3.2097530e+02 },
	{ -6.9161180e+02, 3.0688020e+02,-1.8806670e+02 },
	{  2.8537893e+03, 2.9129830e+02,-8.2416100e+01 },
	{ -2.0407952e+03,-2.9902530e+02, 1.6371910e+02 },
	{  4.9259810e+02, 7.6346100e+01,-6.0024800e+01 }
};

static const double bII[9][2] = { // b^II_ij = bII[j][i]
	{ -1.1628100e+01,-8.6991000e+00 },
	{  1.2560660e+02, 6.3383000e+01 },
	{ -5.3274890e+02,-1.2889390e+02 },
	{  1.1423873e+03,-1.3503120e+02 },
	{ -1.1568545e+03, 9.7758380e+02 },
	{  7.5010200e+01,-1.6499529e+03 },
	{  9.9681140e+02, 1.2586812e+03 },
	{ -8.8818950e+02,-4.0474610e+02 },
	{  2.5013860e+02, 2.7335400e+01 }
};

// Nozawa et al. (2009) Table 4: Coefficients c^II_ij.
static const double cII[7][5] = { // c^II_ij = cII[j][i-2]
	{ -5.7752000e+00, 3.0558600e+01,-5.4327200e+01, 3.6262500e+01,-8.4082000e+00 },
	{  4.6209700e+01,-2.4821770e+02, 4.5096760e+02,-3.1009720e+02, 7.4792500e+01 },
	{ -1.6072800e+02, 8.7419640e+02,-1.6165987e+03, 1.1380531e+03,-2.8295400e+02 },
	{  3.0500700e+02,-1.6769028e+03, 3.1481061e+03,-2.2608347e+03, 5.7639300e+02 },
	{ -3.2954200e+02, 1.8288677e+03,-3.4783930e+03, 2.5419361e+03,-6.6193900e+02 },
	{  1.9107700e+02,-1.0689366e+03, 2.0556693e+03,-1.5252058e+03, 4.0429300e+02 },
	{ -4.6271800e+01, 2.6056560e+02,-5.0567890e+02, 3.8008520e+02,-1.0223300e+02 }
};

// Gf1[i][j-2] = Gamma(i + j/8. + 1.)
static const double Gf1[3][5] = {
	{ 9.0640248e-01, 8.8891357e-01, 8.8622693e-01, 8.9657428e-01, 9.1906253e-01 },
	{ 1.1330031e+00, 1.2222562e+00, 1.3293404e+00, 1.4569332e+00, 1.6083594e+00 },
	{ 2.5492570e+00, 2.9028584e+00, 3.3233510e+00, 3.8244497e+00, 4.4229884e+00 }
};

// Gf1[i][j-2] = Gamma(i + j/8. + 1.)/(i + j/8. + 1.)
static const double Gf2[2][5] = {
	{ 7.2512198e-01, 6.4648260e-01, 5.9081795e-01, 5.5173802e-01, 5.2517859e-01 },
	{ 5.0355693e-01, 5.1463417e-01, 5.3173616e-01, 5.5502217e-01, 5.8485797e-01 }
};

// Nozawa et al. (2009) Table 5: Coefficients a^III_ij and b^III_ij.
static const double aIII[9][3] = { // a^III_ij = aIII[j][i]
	{  5.2163300e+01, 4.9713900e+01, 6.4751200e+01 },
	{ -2.5703130e+02,-1.8977460e+02,-2.1389560e+02 },
	{  4.4681610e+02, 2.7102980e+02, 1.7414320e+02 },
	{ -2.9305850e+02,-2.6978070e+02, 1.3650880e+02 },
	{  0.0000000e+00, 4.2048120e+02,-2.7148990e+02 },
	{  7.7047400e+01,-5.7662470e+02, 8.9321000e+01 },
	{ -2.3871800e+01, 4.3277900e+02, 5.8258400e+01 },
	{  0.0000000e+00,-1.6053650e+02,-4.6080700e+01 },
	{  1.9970000e-01, 2.3392500e+01, 8.7301000e+00 }
};

static const double bIII[9][2] = { // b^III_ij = bIII[j][i]
	{ -8.5862000e+00, 3.7643220e+02 },
	{  3.4134800e+01,-1.2233635e+03 },
	{ -1.1632870e+02, 6.2867870e+02 },
	{  2.9654510e+02, 2.2373946e+03 },
	{ -3.9342070e+02,-3.8288387e+03 },
	{  2.3754970e+02, 2.1217933e+03 },
	{ -3.0600000e+01,-5.5166700e+01 },
	{ -2.7617000e+01,-3.4943210e+02 },
	{  8.8453000e+00, 9.2205900e+01 }
};

STATIC double eeBremsCooling(double Te,   // electron temperature in K
			     double eden) // electron density in cm^-3
{
	DEBUG_ENTRY( "eeBremsCooling()" );

	// Evaluate e-e brems cooling; the fitting functions are described in
	// >>refer	ee	brems	Nozawa S., Takahashi K., Kohyama Y., Itoh N., 2009, A&A 499, 661
	// >>refer	ee	brems	Itoh N., Kawana Y., Nozawa S., 2002, Il Nuovo Cimento, 117B, 359

	double mec2 = ELECTRON_MASS*pow2(SPEEDLIGHT);
	double tau = BOLTZMANN*Te/mec2;
	// common factor, see Eq. 11 of Nozawa et al. (2009)
	double fac = mec2*pow2(eden)*SIGMA_THOMSON*SPEEDLIGHT*FINE_STRUCTURE*powpq(tau,3,2);
	double G = 0.;
	if( Te <= 5.754399e5 )
	{
		// powerlaw extrapolation to get reasonable behavior at low Te
		G = 1.680322*pow(Te/5.754399e5,0.116);
	}
	else if( Te <= 1.160451e7 )
	{
		// Eq. 28 of Nozawa et al. (2009)
		double Theta = (log10(tau) + 2.65)/1.35;
		// region I (50 eV < kTe <= 1 keV), Eq. 31 of Nozawa et al. (2009)
		double Thp = 1.;
		for( int i=0; i < 11; ++i )
		{
			G += bI[i]*Thp;
			Thp *= Theta;
		}
		G *= sqrt(8./(3.*PI));
	}
	else if( Te <= 3.481352e9 )
	{
		double A[3]={0.,0.,0.}, B[2]={0.,0.}, C[5]={0.,0.,0.,0.,0.};
		// region II (1 keV < kTe <= 300 keV), Eq. 39 of Nozawa et al. (2009)
		double tau8 = powpq(tau,1,8); // tau^(1/8)
		double t8p = 1.;
		for( int j=0; j < 9; ++j )
		{
			A[0] += aII[j][0]*t8p;
			A[1] += aII[j][1]*t8p;
			A[2] += aII[j][2]*t8p;
			B[0] += bII[j][0]*t8p;
			B[1] += bII[j][1]*t8p;
			t8p *= tau8;
		}
		double tau6 = powpq(tau,1,6); // tau^(1/6)
		double t6p = 1.;
		for( int j=0; j < 7; ++j )
		{
			C[0] += cII[j][0]*t6p;
			C[1] += cII[j][1]*t6p;
			C[2] += cII[j][2]*t6p;
			C[3] += cII[j][3]*t6p;
			C[4] += cII[j][4]*t6p;
			t6p *= tau6;
		}
		G = A[0] + A[1] + 2.*A[2] + B[0] + 0.5*B[1];
		for( int i=0; i < 3; ++i )
			for( int j=2; j < 7; ++j )
				G += Gf1[i][j-2]*A[i]*C[j-2];
		for( int i=0; i < 2; ++i )
			for( int j=2; j < 7; ++j )
				G += Gf2[i][j-2]*B[i]*C[j-2];
	}
	else
	{
		double A[3]={0.,0.,0.}, B[2]={0.,0.};
		// region III (300 keV < Te <= 7 Mev), Eq. 45 of Nozawa et al. (2009)
		double tau8 = powpq(tau,1,8); // tau^(1/8)
		double t8p = 1.;
		for( int j=0; j < 9; ++j )
		{
			A[0] += aIII[j][0]*t8p;
			A[1] += aIII[j][1]*t8p;
			A[2] += aIII[j][2]*t8p;
			B[0] += bIII[j][0]*t8p;
			B[1] += bIII[j][1]*t8p;
			t8p *= tau8;
		}
		G = A[0] + A[1] + 2.*A[2] + B[0] + 0.5*B[1];
	}
	return fac*G;
}

void eeBremsSpectrum(double Te,            // electron temperature in K
					 vector<realnum>& jnu, // 4pi*jnu in photons/cm^3/s/cell
					 vector<double>& knu)  // opacity in cm^-1
{
	DEBUG_ENTRY( "eeBremsSpectrum()" );

	// Calculate e-e brems spectrum; the fitting functions are described in
	// >>refer	ee	brems	Nozawa S., Takahashi K., Kohyama Y., Itoh N., 2009, A&A 499, 661
	// >>refer	ee	brems	Itoh N., Kawana Y., Nozawa S., 2002, Il Nuovo Cimento, 117B, 359

	double mec2 = ELECTRON_MASS*pow2(SPEEDLIGHT);
	// if Te is below 5.754399e5, use the shape for 5.754399e5 K
	double tau = BOLTZMANN*max(Te,5.754399e5)/mec2;
	// common factor, see Eq. 6 of Nozawa et al. (2009)
	// multiplication by ne^2 is done outside this routine to avoid needless reevaluations
	double fac = SIGMA_THOMSON*SPEEDLIGHT*FINE_STRUCTURE/sqrt(tau)*FR1RYD*HPLANCK/mec2;
	// adjust normalization if we are extrapolating below 50 eV
	// the extrapolation for the spectrum used below gives a Te^2 dependence,
	// eeBremsCooling() requires it to be Te^1.616, so here we correct the difference
	if( Te < 5.754399e5 )
		fac *= pow(Te/5.754399e5,-0.384);
	// fits are valid between 1e-4 <= hnu/kTe <= 10; convert these linits to Ryd here
	double Elo = max( 1.e-4*Te/TE1RYD, rfield.anu(0) );
	long klo = fp_equal(Elo,rfield.anu(0)) ? 0 : rfield.ithreshC(Elo);
	double Ehi = min( 10.*Te/TE1RYD, rfield.anu(rfield.nflux) );
	long khi = fp_equal(Ehi,rfield.anu(rfield.nflux)) ? rfield.nflux : rfield.ithreshC(Ehi);
	// normalization constant for Bnu
	double fac2 = PI8*pow3(FR1RYD)/pow2(SPEEDLIGHT);
	avx_ptr<double> arg(rfield.nflux), boltz(rfield.nflux);
	double* bfac;
	if( fp_equal( Te, phycon.te ) )
	{
		bfac = rfield.ContBoltz.data();
	}
	else
	{
		for( long k=klo; k < khi; ++k )
			arg[k] = rfield.anu(k)*TE1RYD/Te;
		vexp( arg.ptr0(), boltz.ptr0(), klo, khi );
		bfac = boltz.ptr0();
	}

	if( Te <= 1.160451e7 )
	{
		// region I (50 eV < kTe <= 1 keV), Eq. 26 of Nozawa et al. (2009)
		//
		// Below 50 eV (is 5.754399e5 K) we extrapolate the expressions assuming
		// that P_ee(x) is independent of Te. Test between Te = 1e6 and 1e7 K show
		// that this is reasonable. The normalization constant is then adjusted to
		// give the correct integral in accordance with what eeBremsCooling()
		// returns.
		//
		double Theta = (log10(tau) + 2.65)/1.35; // Eq. 28 of Nozawa et al. (2009)
		double Thp[11];
		Thp[0] = 1.;
		for( int i=1; i < 11; ++i )
			Thp[i] = Thp[i-1]*Theta;

		avx_ptr<double> x(klo,khi), y(klo,khi), z(klo,khi);
		double hokT = TE1RYD/Te;
		for( long k=klo; k < khi; ++k )
			x[k] = rfield.anu(k)*hokT;
		vlog10( x.ptr0(), y.ptr0(), klo, khi );
		vexpm1( x.ptr0(), z.ptr0(), klo, khi );
		for( long k=klo; k < khi; ++k )
		{
			double xx = (y[k] + 1.5)/2.5; // Eq. 29 of Nozawa et al. (2009)
			double xxp[11];
			xxp[0] = 1.;
			for( int i=1; i < 11; ++i )
				xxp[i] = xxp[i-1]*xx;

			double G = 0.;
			for( int i=0; i < 11; ++i )
				for( int j=0; j < 11; ++j )
					G += aI[i][j]*Thp[i]*xxp[j];
			G *= sqrt(8./(3.*PI));
			double j_nu = fac*bfac[k]/x[k]*G*rfield.widflx(k);
			// calculate opacity from Kirchhoff's law, use 4pi*Bnu in photons/cm^2/s/cell
			double Bnu = fac2*rfield.anu2(k)*rfield.widflx(k)/z[k];
			jnu[k] = realnum(j_nu);
			knu[k] = j_nu/Bnu;
		}
	}
	else if( Te <= 3.481352e9 )
	{
		// region II (1 keV < kTe <= 300 keV), Eq. 32 of Nozawa et al. (2009)
		double A[3]={0.,0.,0.}, B[2]={0.,0.}, C[5]={0.,0.,0.,0.,0.};
		double tau8 = powpq(tau,1,8); // tau^(1/8)
		double t8p = 1.;
		for( int j=0; j < 9; ++j )
		{
			A[0] += aII[j][0]*t8p;
			A[1] += aII[j][1]*t8p;
			A[2] += aII[j][2]*t8p;
			B[0] += bII[j][0]*t8p;
			B[1] += bII[j][1]*t8p;
			t8p *= tau8;
		}
		double tau6 = powpq(tau,1,6); // tau^(1/6)
		double t6p = 1.;
		for( int j=0; j < 7; ++j )
		{
			C[0] += cII[j][0]*t6p;
			C[1] += cII[j][1]*t6p;
			C[2] += cII[j][2]*t6p;
			C[3] += cII[j][3]*t6p;
			C[4] += cII[j][4]*t6p;
			t6p *= tau6;
		}

		for( long k=klo; k < khi; ++k )
		{
			double x = rfield.anu(k)*TE1RYD/Te;
			// use the identity -Ei(-x) = E1(x)
			double G = (A[2]*x + A[1])*x + A[0] + e1_scaled(x)*(B[1]*x + B[0]);
			double F = 1.;
			double x8 = powpq(x,1,8); // x^1/8
			double x8p = x8*x8;
			for( int i=2; i < 7; ++i )
			{
				F += C[i-2]*x8p;
				x8p *= x8;
			}
			double j_nu = fac*bfac[k]/x*G*F*rfield.widflx(k);
			// calculate opacity from Kirchhoff's law, use 4pi*Bnu in photons/cm^2/s/cell
			double Bnu = fac2*rfield.anu2(k)*rfield.widflx(k)/expm1(x);
			jnu[k] = realnum(j_nu);
			knu[k] = j_nu/Bnu;
		}
	}
	else
	{
		// region III (300 keV < Te <= 7 Mev), Eq. 40 of Nozawa et al. (2009)
		double A[3]={0.,0.,0.}, B[2]={0.,0.};
		double tau8 = powpq(tau,1,8); // tau^(1/8)
		double t8p = 1.;
		for( int j=0; j < 9; ++j )
		{
			A[0] += aIII[j][0]*t8p;
			A[1] += aIII[j][1]*t8p;
			A[2] += aIII[j][2]*t8p;
			B[0] += bIII[j][0]*t8p;
			B[1] += bIII[j][1]*t8p;
			t8p *= tau8;
		}

		for( long k=klo; k < khi; ++k )
		{
			double x = rfield.anu(k)*TE1RYD/Te;
			// use the identity -Ei(-x) = E1(x)
			double G = (A[2]*x + A[1])*x + A[0] + e1_scaled(x)*(B[1]*x + B[0]);
			double j_nu = fac*bfac[k]/x*G*rfield.widflx(k);
			// calculate opacity from Kirchhoff's law, use 4pi*Bnu in photons/cm^2/s/cell
			double Bnu = fac2*rfield.anu2(k)*rfield.widflx(k)/expm1(x);
			jnu[k] = realnum(j_nu);
			knu[k] = j_nu/Bnu;
		}
	}

	//
	// extrapolate the emissivity and opacity to cover the entire frequency range of Cloudy
	//
	// Itoh calculated the total cooling rate by integrating his fits from 0 to infinity.
	// This implies a different extrapolation than what we use below, and therefore the
	// integral over jnu[] is not exactly equal to what eeBremsCooling() reports above.
	// These two always agree within 0.03% though, with the largest discrepancy at 10 GK.
	//
	double x1l = rfield.anuln(klo);
	double x2l = rfield.anuln(klo+1);
	double y1 = jnu[klo]/rfield.widflx(klo);
	double y2 = jnu[klo+1]/rfield.widflx(klo+1);
	double z1 = knu[klo];
	double slope = log(y1/y2)/(x1l-x2l);
	// extrapolate downwards as jnu = C*nu^slope
	avx_ptr<double> efac(rfield.nflux);
	for( long k=0; k < klo; ++k )
	{
		x2l = rfield.anuln(k);
		arg[k] = slope*(x2l-x1l);
	}
	vexp( arg.ptr0(), efac.ptr0(), 0, klo );
	for( long k=0; k < klo; ++k )
	{
		jnu[k] = realnum(y1*rfield.widflx(k)*efac[k]);
		double af = rfield.anu(klo)/rfield.anu(k);
		knu[k] = z1*efac[k]*af;
	}
	x1l = rfield.anuln(khi-1);
	double x1 = rfield.anu(khi-1);
	x2l = rfield.anuln(khi-2);
	double x2 = rfield.anu(khi-2);
	y1 = jnu[khi-1]/rfield.widflx(khi-1);
	y2 = jnu[khi-2]/rfield.widflx(khi-2);
	z1 = knu[khi-1];
	slope = (log(y2/y1) + (x2-x1)*TE1RYD/Te)/(x2l-x1l);
	// extrapolate upwards as jnu = C*nu^slope*exp(-hnu/kTe)
	avx_ptr<double> arg2(rfield.nflux), efac2(rfield.nflux);
	for( long k=khi; k < rfield.nflux; ++k )
	{
		x2 = rfield.anu(k);
		x2l = rfield.anuln(k);
		arg[k] = (x1-x2)*TE1RYD/Te - slope*(x1l-x2l);
		arg2[k] = (slope-2.)*(x2l-x1l);
	}
	vexp( arg.ptr0(), efac.ptr0(), khi, rfield.nflux );
	vexp( arg2.ptr0(), efac2.ptr0(), khi, rfield.nflux );
	for( long k=khi; k < rfield.nflux; ++k )
	{
		jnu[k] = realnum(y1*rfield.widflx(k)*efac[k]);
		knu[k] = z1*efac2[k];
	}
}
