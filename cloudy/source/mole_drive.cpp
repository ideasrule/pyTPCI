/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*mole_drive - public routine, calls mole_solve to converge molecules */
#include "cddefines.h"
#include "dense.h"
#include "hmi.h"
#include "thermal.h"
#include "iso.h"
#include "ion_trim.h"
#include "phycon.h"
#include "radius.h"
#include "secondaries.h"
#include "timesc.h"
#include "trace.h"
#include "co.h"
#include "h2.h"
#include "mole_priv.h"
#include "mole.h"
#include "freebound.h"

/* this is the error tolerance for reporting non-convergence */
static const double MOLETOLER = 0.10;

STATIC void mole_effects(void);
STATIC void mole_h_rate_diagnostics(void);
STATIC void mole_update_limiting_reactants();

/*mole_drive main driver for heavy molecular equilibrium routines      */
void mole_drive(void)
{
	DEBUG_ENTRY( "mole_drive()" );

	mole_update_species_cache();  /* Update densities of species controlled outside the chemical network */

	mole_update_limiting_reactants();
 	
	mole_update_rks();
	
	double error = mole_solve();	
	
	bool lgConverged = (error < MOLETOLER); 

	if( !lgConverged )
	{
		// should something be done with this?
	}

	mole_ion_trim();

	mole_effects();
	
	return;
}

void mole_update_sources(void)
{
	DEBUG_ENTRY( "mole_update_sources()" );

	mole_update_species_cache();
	mole_eval_sources(mole_global.num_total);
}

STATIC void mole_effects()
{
	double 
		plte;

	DEBUG_ENTRY( "mole_effects()" );

	mole_eval_sources(mole_global.num_total);

	co.CODissHeat = (realnum)(mole.findrate("PHOTON,CO=>C,O")*1e-12);

	thermal.setHeating(0,9, co.CODissHeat);

	/* total H2 - all forms */
	hmi.H2_total = findspecieslocal("H2")->den + findspecieslocal("H2*")->den;
	hmi.HD_total = findspecieslocal("HD")->den + findspecieslocal("HD*")->den;
	/* first guess at ortho and para densities */
	h2.ortho_density = 0.75*hmi.H2_total;
	h2.para_density = 0.25*hmi.H2_total;
	{
		hmi.H2_total_f = (realnum)hmi.H2_total;
		h2.ortho_density_f = (realnum)h2.ortho_density;
		h2.para_density_f = (realnum)h2.para_density;
	}

	/* save rate H2 is destroyed units s-1 */
	/* >>chng 05 mar 18, TE, add terms - 
		total destruction rate is: dest_tot = n_H2g/n_H2tot * dest_H2g + n_H2s/n_H2tot * dest_H2s */
	/* as reactions that change H2s to H2g and vice versa are not counted destruction processes, the terms c[ipMH2g][ipMH2s] *
	   and c[ipMH2s][ipMH2g], which have a different sign than [ipMH2g][ipMH2g] and [ipMH2s][ipMH2s], have to be added	*/
	hmi.H2_rate_destroy = mole.sink_rate_tot("H2");

	/* establish local timescale for H2 molecule destruction */
	if(hmi.H2_rate_destroy > SMALLFLOAT )
	{
		/* units are now seconds */
		timesc.time_H2_Dest_here = 1./hmi.H2_rate_destroy;
	}
	else
	{
		timesc.time_H2_Dest_here = 0.;
	}
	
	/* local timescale for H2 formation 
	 * both grains and radiative attachment 
	 * >>chng 08 mar 15 rjrw -- ratio formation rate to density of _H2_ so comparable with destruction rate,
	 * use full rates not rate coefficients so have correct units */
	timesc.time_H2_Form_here = (mole.findrate("H,H,grn=>H2*,grn")+mole.findrate("H,H,grn=>H2,grn"));
	/* timescale is inverse of this rate */
	if( timesc.time_H2_Form_here > SMALLFLOAT )
	{
		/* units are now seconds */
		timesc.time_H2_Form_here =  hmi.H2_total/timesc.time_H2_Form_here;
	}
	else
	{
		timesc.time_H2_Form_here = 0.;
	}


	/* heating due to H2 dissociation */

	/* H2 photodissociation heating, eqn A9 of Tielens & Hollenbach 1985a */
	/*hmi.HeatH2Dish_TH85 = (float)(1.36e-23*findspecieslocal("H2")->den*h2esc*hmi.UV_Cont_rel2_Habing_TH85_depth);*/
	/* >>chng 04 feb 07, more general to express heating in terms of the assumed
	 * photo rates - the 0.25 was obtained by inverting A8 & A9 of TH85 to find that
	 * there are 0.25 eV per dissociative pumping, ie, 10% of total 
	 * this includes both H2g and H2s - TH85 say just ground but they include
	 * process for both H2 and H2s - as we did above - both must be in
	 * heating term */
	/* >>chng 05 mar 11, TE, old had used H2_Solomon_dissoc_rate_used, which was only
	 * H2g.  in regions where Solomon process is fast, H2s has a large population
	 * and the heating rate was underestimated. */
	/* >>chng 05 jun 23, 
	 * >>chng 05 dec 05, TE, modified to approximate the heating better for the
	 * new approximation */
	/* >>chng 00 nov 25, explicitly break out this heat agent */
	/* 2.6 eV of heat per deexcitation, consider difference
	 * between deexcitation (heat) and excitation (cooling) */
	/* >>chng 04 jan 27, code moved here and simplified */
	/* >>chng 05 jul 10, GS*/ 
	/*  average collisional rate for H2* to H2g calculated from big H2, GS*/
	
	/* TH85 dissociation heating - this is ALWAYS defined for reference
	 * may be output for comparison with other rates*/
	hmi.HeatH2Dish_TH85 = 0.25 * EN1EV *
		(hmi.H2_Solomon_dissoc_rate_used_H2g * findspecieslocal("H2")->den +
		 hmi.H2_Solomon_dissoc_rate_used_H2s * findspecieslocal("H2*")->den);
	
	/* TH85 deexcitation heating */
	hmi.HeatH2Dexc_TH85 = ((mole.findrate("H2*,H2=>H2,H2") + mole.findrate("H2*,H=>H2,H"))
												 - (mole.findrate("H2,H2=>H2*,H2") + mole.findrate("H2,H=>H2*,H"))) * 4.17e-12;
	/* this is derivative wrt temperature, only if counted as a cooling source */
	hmi.deriv_HeatH2Dexc_TH85 = (realnum)(MAX2(0.,-hmi.HeatH2Dexc_TH85)* ( 30172. * thermal.tsq1 - thermal.halfte ) );
	
	if( hmi.chH2_small_model_type == 'H' )
	{
		/* Burton et al. 1990 */
		hmi.HeatH2Dish_BHT90 = 0.25 * EN1EV *
			(hmi.H2_Solomon_dissoc_rate_used_H2g * findspecieslocal("H2")->den +
			 hmi.H2_Solomon_dissoc_rate_used_H2s * findspecieslocal("H2*")->den);
		
		/* Burton et al. 1990 heating */
		hmi.HeatH2Dexc_BHT90 = ((mole.findrate("H2*,H2=>H2,H2") + mole.findrate("H2*,H=>H2,H")) 
														- (mole.findrate("H2,H2=>H2*,H2") + mole.findrate("H2,H=>H2*,H"))) * 4.17e-12;
		/* this is derivative wrt temperature, only if counted as a cooling source */
		hmi.deriv_HeatH2Dexc_BHT90 = (realnum)(MAX2(0.,-hmi.HeatH2Dexc_BHT90)* ( 30172. * thermal.tsq1 - thermal.halfte ) );
	}
	else if( hmi.chH2_small_model_type == 'B')
	{
		/* Bertoldi & Draine */
		hmi.HeatH2Dish_BD96 = 0.25 * EN1EV *
			(hmi.H2_Solomon_dissoc_rate_used_H2g * findspecieslocal("H2")->den +
			 hmi.H2_Solomon_dissoc_rate_used_H2s * findspecieslocal("H2*")->den);
		/* Bertoldi & Draine heating, same as TH85 */
		hmi.HeatH2Dexc_BD96 = ((mole.findrate("H2*,H2=>H2,H2") + mole.findrate("H2*,H=>H2,H"))
													 - (mole.findrate("H2,H2=>H2*,H2") + mole.findrate("H2,H=>H2*,H"))) * 4.17e-12;
		/* this is derivative wrt temperature, only if counted as a cooling source */
		hmi.deriv_HeatH2Dexc_BD96 = (realnum)(MAX2(0.,-hmi.HeatH2Dexc_BD96)* ( 30172. * thermal.tsq1 - thermal.halfte ) );
	}
	else if(hmi.chH2_small_model_type == 'E')
	{
		/* heating due to dissociation of H2
		 * >>chng 05 oct 19, TE, define new approximation for the heating due to the destruction of H2
		 *	use this approximation for the specified cloud parameters, otherwise
		 * use limiting cases for 1 <= G0, G0 >= 1e7, n >= 1e7, n <= 1 */
		
		double log_density, 
			f1, f2,f3, f4, f5;
		static double log_G0_face = -1;
		static double k_f4;
		
		
		/* test for G0 
		 * this is a constant so only do it in zone 0 */
		if( !nzone )
		{
			if(hmi.UV_Cont_rel2_Draine_DB96_face <= 1.) 
			{ 
				log_G0_face = 0.;
			}
			else if(hmi.UV_Cont_rel2_Draine_DB96_face >= 1e7) 
			{ 
				log_G0_face = 7.;
			}
			else 
			{ 
				log_G0_face = log10(hmi.UV_Cont_rel2_Draine_DB96_face); 
			}
			/*>>chng 06 oct 24 TE change Go face for spherical geometry*/
			log_G0_face /= radius.r1r0sq;
		}

		/* test for density */
		if(dense.gas_phase[ipHYDROGEN] <= 1.) 
		{ 
			log_density = 0.; 
		}
		else if(dense.gas_phase[ipHYDROGEN] >= 1e7) 
		{ 
			log_density = 7.; 
		}
		else 
		{ 
			log_density = log10(dense.gas_phase[ipHYDROGEN]); 
		}
		
		f1 = 0.15 * log_density + 0.75;
		f2 = -0.5 * log_density + 10.;
		
		hmi.HeatH2Dish_ELWERT = 0.25 * EN1EV *  f1 * 
			(hmi.H2_Solomon_dissoc_rate_used_H2g * findspecieslocal("H2")->den +
			 hmi.H2_Solomon_dissoc_rate_used_H2s * findspecieslocal("H2*")->den ) + 
			f2 * secondaries.x12tot * EN1EV * hmi.H2_total;
		
		/*fprintf( ioQQQ, "f1: %.2e, f2: %.2e,heat Solomon: %.2e",f1,f2,hmi.HeatH2Dish_TH85);*/
		
		
		/* heating due to collisional deexcitation within X of H2
		 * >>chng 05 oct 19, TE, define new approximation for the heating due to the destruction of H2
		 *	use this approximation for the specified cloud parameters, otherwise
		 * use limiting cases for 1 <= G0, G0 >= 1e7, n >= 1e7, n <= 1 */
		
		/* set value of k_f4 by testing on value of G0 */
		if(hmi.UV_Cont_rel2_Draine_DB96_face <= 1.) 
		{ 
			log_G0_face = 0.;
		}
		else if(hmi.UV_Cont_rel2_Draine_DB96_face >= 1e7) 
		{ 
			log_G0_face = 7.;
		}
		else 
		{ 
			log_G0_face = log10(hmi.UV_Cont_rel2_Draine_DB96_face); 
		}
		/* 06 oct 24, TE introduce effects of spherical geometry */
		log_G0_face /= radius.r1r0sq;
		
		/* terms only dependent on G0_face */
		k_f4 = -0.25 * log_G0_face + 1.25;

		/* test for density */
		if(dense.gas_phase[ipHYDROGEN] <= 1.) 
		{ 
			log_density = 0.; 
			f4 = 0.; 
		}
		else if(dense.gas_phase[ipHYDROGEN] >= 1e7) 
		{ 
			log_density = 7.; 
			f4 = pow2(k_f4) * exp10(  2.2211 * log_density - 29.8506);
		}
		else 
		{ 
			log_density = log10(dense.gas_phase[ipHYDROGEN]); 
			f4 = pow2(k_f4) * exp10(  2.2211 * log_density - 29.8506);
		}
		
		f3 = MAX2(0.1, -4.75 * log_density + 24.25);
		f5 = MAX2(1.,0.95 * log_density - 1.45) * 0.2 * log_G0_face;
		
		hmi.HeatH2Dexc_ELWERT = ((mole.findrate("H2*,H2=>H2,H2") + mole.findrate("H2*,H=>H2,H"))
			- (mole.findrate("H2,H2=>H2*,H2") + mole.findrate("H2,H=>H2*,H"))) * 4.17e-12 * f3 + 
			f4 * (mole.species[findnuclide("H")->ipMl[0]].den/dense.gas_phase[ipHYDROGEN]) +
			f5 * secondaries.x12tot * EN1EV * hmi.H2_total;
		
		if(log_G0_face == 0.&& dense.gas_phase[ipHYDROGEN] > 1.) 
			hmi.HeatH2Dexc_ELWERT *= 0.001 / dense.gas_phase[ipHYDROGEN];
		
		/* >>chng 06 oct 24, TE introduce effects of spherical geometry */
		/*if(radius.depth/radius.rinner >= 1.0) */
		hmi.HeatH2Dexc_ELWERT /= POW2(radius.r1r0sq);
		
		/* this is derivative wrt temperature, only if counted as a cooling source */
		hmi.deriv_HeatH2Dexc_ELWERT = (realnum)(MAX2(0.,-hmi.HeatH2Dexc_ELWERT)* ( 30172. * thermal.tsq1 - thermal.halfte ) );
		
		/*fprintf( ioQQQ, "\tf3: %.2e, f4: %.2e, f5: %.2e, heat coll dissoc: %.2e\n",f3,f4,f5,hmi.HeatH2Dexc_TH85);*/
	}
	/* end Elwert branch for photo rates */
	else
		TotalInsanity();


	if( findspecieslocal("H-")->den > 0. && hmi.rel_pop_LTE_Hmin > 0. )
	{
		hmi.hmidep = (double)findspecieslocal("H-")->den/ SDIV( 
			((double)dense.xIonDense[ipHYDROGEN][0]*dense.eden*hmi.rel_pop_LTE_Hmin));
	}
	else
	{
		hmi.hmidep = 1.;
	}

	/* this will be net volume heating rate, photo heat - induc cool */
	hmi.hmihet = hmi.HMinus_photo_heat*findspecieslocal("H-")->den - hmi.HMinus_induc_rec_cooling;

	/* H2+  +  HNU  =>  H+  + H */
	hmi.h2plus_heat = hmi.h2plus_heatcoef * findspecieslocal("H2+")->den;

	/* departure coefficient for H2 defined rel to n(1s)**2
	 * (see Littes and Mihalas Solar Phys 93, 23) */
	plte = (double)dense.xIonDense[ipHYDROGEN][0] * hmi.rel_pop_LTE_H2g * (double)dense.xIonDense[ipHYDROGEN][0];
	if( plte > 0. )
	{
		hmi.h2dep = findspecieslocal("H2")->den/plte;
	}
	else
	{
		hmi.h2dep = 1.;
	}

	/* departure coefficient of H2+ defined rel to n(1s) n(p)
	 * sec den was HI before 85.34 */
	plte = (double)dense.xIonDense[ipHYDROGEN][0]*hmi.rel_pop_LTE_H2p*(double)dense.xIonDense[ipHYDROGEN][1];
	if( plte > 0. )
	{
		hmi.h2pdep = findspecieslocal("H2+")->den/plte;
	}
	else
	{
		hmi.h2pdep = 1.;
	}

	/* departure coefficient of H3+ defined rel to N(H2+) N(p) */
	if( hmi.rel_pop_LTE_H3p > 0. )
	{
		hmi.h3pdep = findspecieslocal("H3+")->den/hmi.rel_pop_LTE_H3p;
	}
	else
	{
		hmi.h3pdep = 1.;
	}

	mole_h_rate_diagnostics();

	return;
}

STATIC void mole_h_rate_diagnostics(void)
{
	int i;
	double 
		rate; 

	/* identify dominant H2 formation process */
	{
		/* following should be set true to identify H2 formation and destruction processes */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && (nzone>50) )
		{
			double createsum ,create_from_Hn2 , create_3body_Ho, create_h2p, 
				create_h3p, create_h3pe, create_grains, create_hminus;
			double destroysum, destroy_hm ,destroy_solomon ,destroy_2h ,destroy_hp,
				destroy_h,destroy_hp2,destroy_h3p;

			create_from_Hn2 = mole.findrate("H,H=>H2"); 			/* H(n=2) + H(n=1) -> H2 */
			create_3body_Ho = mole.findrate("H,H,H=>H2,H");
			create_h2p = mole.findrate("H,H2+=>H+,H2*");
			create_h3p = mole.findrate("H,H3+=>H2,H2+");
			create_h3pe = mole.findrate("e-,H3+=>H2,H"); 
			create_grains = mole.findrate("H,H,grn=>H2*,grn")+mole.findrate("H,H,grn=>H2,grn");
			create_hminus = (mole.findrate("H,H-=>H2,e-")+mole.findrate("H,H-=>H2*,e-"));

			createsum = create_from_Hn2 + create_3body_Ho + create_h2p +
				create_h3p + create_h3pe + create_grains + create_hminus;

			fprintf(ioQQQ,"H2 create zone\t%.2f \tsum\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				fnzone,
				createsum,
				create_hminus   / createsum,
				create_from_Hn2 / createsum,
				create_3body_Ho / createsum,
				create_h2p      / createsum,
				create_h3p      / createsum,
				create_h3pe     / createsum,
				create_grains   / createsum	);

			/* all h2 molecular hydrogen destruction processes */
			/* >>chng 04 jan 28, had wrong Boltzmann factor,
			 * caught by gargi Shaw */
			destroy_hm = (mole.findrate("H2,e-=>H,H-")+mole.findrate("H2*,e-=>H,H-"));
			/*destroy_hm2 = eh2hhm;*/
			destroy_solomon = hmi.H2_Solomon_dissoc_rate_used_H2g * findspecieslocal("H2")->den;
			destroy_2h = (mole.findrate("H2,e-=>H,H,e-")+mole.findrate("H2*,e-=>H,H,e-"));
			destroy_hp = mole.findrate("H2,H+=>H3+");
			destroy_h = mole.findrate("H,H2=>H,H,H");
			destroy_hp2 = mole.findrate("H2,H+=>H2+,H");
			destroy_h3p = mole.findrate("H2,H3+=>H2,H2+,H");
			destroysum = destroy_hm + /*destroy_hm2 +*/ destroy_solomon + destroy_2h + 
				destroy_hp+ destroy_h+ destroy_hp2+ destroy_h3p;

			fprintf(ioQQQ,"H2 destroy\t%.3f \t%.2e\tsum\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				fnzone,
				destroysum,
				destroy_hm / destroysum ,
				destroy_solomon / destroysum ,
				destroy_2h / destroysum ,
				destroy_hp / destroysum ,
				destroy_h / destroysum ,
				destroy_hp2 / destroysum ,
				destroy_h3p / destroysum );

		}
	}

	{
		/* following should be set true to identify H- formation and destruction processes */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && (nzone>140)/**/ )
		{
			double create_from_Ho,create_3body_Ho,create_batach,destroy_photo,
				destroy_coll_heavies,destroy_coll_electrons,destroy_Hattach,destroy_fhneut,
				destsum , createsum;

			create_from_Ho = mole.findrate("H,e-=>H-,PHOTON");
			create_3body_Ho = mole.findrate("H,e-,e-=>H-,e-");
			/* total formation is sum of g and s attachment */
			create_batach = (mole.findrate("H2,e-=>H,H-") + mole.findrate("H2*,e-=>H,H-"));
			destroy_photo = mole.findrate("H-,PHOTON=>H,e-");
			destroy_coll_heavies = 0.;			
			for(ChemNuclideList::iterator atom = nuclide_list.begin(); atom != nuclide_list.end(); ++atom )
			{
				if( !(*atom)->lgHasLinkedIon())
					continue;
				long nelem = (*atom)->el()->Z-1;
				if (dense.lgElmtOn[nelem] && nelem > ipHELIUM) 
				{
					char react[32];
					sprintf(react,"H-,%s=>H,%s",mole_global.list[(*atom)->ipMl[1]]->label.c_str(),(*atom)->label().c_str() );
					destroy_coll_heavies += mole.findrate(react);
				}
			}
			destroy_coll_electrons = mole.findrate("H-,e-=>H-,e-,e-");
			destroy_Hattach = (mole.findrate("H,H-=>H2,e-")+mole.findrate("H,H-=>H2*,e-"));
			destroy_fhneut = mole.findrate("H-,H+=>H,H");

			destsum = destroy_photo + destroy_coll_heavies + destroy_coll_electrons + 
				destroy_Hattach + destroy_fhneut;
			fprintf(ioQQQ,"H- destroy zone\t%.2f\tTe\t%.4e\tsum\t%.2e\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", 
				fnzone,
				phycon.te,
				destsum,
				destroy_photo/destsum , 
				destroy_coll_heavies/destsum,
				destroy_coll_electrons/destsum, 
				destroy_Hattach/destsum,
				destroy_fhneut/destsum );

			createsum = create_from_Ho+create_3body_Ho+create_batach;
			fprintf(ioQQQ,"H- create\t%.2f\tTe\t%.4e\tsum\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				fnzone,
				phycon.te,
				createsum,
				dense.eden,
				create_from_Ho/createsum,
				create_3body_Ho/createsum,
				create_batach/createsum);
		}
	}
	{
		/* following should be set true to print populations */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			/* these are the raw results */
			fprintf( ioQQQ, " MOLE raw; hi\t%.2e" , dense.xIonDense[ipHYDROGEN][0]);
			for( i=0; i < mole_global.num_calc; i++ )
			{
				fprintf( ioQQQ, "\t%-4.4s\t%.2e", mole_global.list[i]->label.c_str(), mole.species[i].den );
			}
			fprintf( ioQQQ, " \n" );
		}
	}

	if( trace.lgTrace && trace.lgTraceMole )
	{
		/* these are the raw results */
		fprintf( ioQQQ, " raw; " );
		for( i=0; i < mole_global.num_calc; i++ )
		{
			fprintf( ioQQQ, " %-4.4s:%.2e", mole_global.list[i]->label.c_str(), mole.species[i].den );
		}
		fprintf( ioQQQ, " \n" );

	}

	/* option to print rate H2 forms */
	/* trace.lgTraceMole is trace molecules option,
	 * punch htwo */
	if( (trace.lgTrace && trace.lgTraceMole) )
	{
		/** total H2 creation rate, cm-3 s-1 */
		double H2_rate_create = mole.source_rate_tot("H2") + mole.source_rate_tot("H2*");

		if( H2_rate_create > SMALLFLOAT )
		{
			fprintf( ioQQQ, " Create H2, rate=%10.2e grain;%5.3f hmin;%5.3f bhedis;%5.3f h2+;%5.3f hmi.radasc:%5.3f hmi.h3ph2p:%5.3f hmi.h3petc:%5.3f\n", 
				H2_rate_create, 
				(mole.findrate("H,H,grn=>H2*,grn")+mole.findrate("H,H,grn=>H2,grn"))/H2_rate_create, 
				(mole.findrate("H,H-=>H2,e-")+mole.findrate("H,H-=>H2*,e-"))/H2_rate_create, 
				mole.findrate("H,H,H=>H2,H")/H2_rate_create, 
				mole.findrate("H,H2+=>H+,H2*")/H2_rate_create, 
				mole.findrate("H,H=>H2")/H2_rate_create, 
				mole.findrate("H,H3+=>H2,H2+")/H2_rate_create, 
				mole.findrate("H2,H3+=>H2,H2+,H")/H2_rate_create );
		}
		else
		{
			fprintf( ioQQQ, " Create H2, rate=0\n" );
		}
	}

	/* this is H2+ */
	if( trace.lgTrace && trace.lgTraceMole )
	{
		rate = mole.findrate("H2,H+=>H2+,H") + mole.findrate("H,H+,e-=>H2+,e-") + 
		  mole.findrate("H,H3+=>H2,H2+") + mole.findrate("H2,H3+=>H2,H2+,H");
		if( rate > 1e-25 )
		{
			fprintf( ioQQQ, " Create H2+, rate=%10.2e hmi.rh2h2p;%5.3f b2pcin;%5.3f hmi.h3ph2p;%5.3f hmi.h3petc+;%5.3f\n", 
			  rate, mole.findrate("H2,H+=>H2+,H")/rate, 
							 mole.findrate("H,H+,e-=>H2+,e-")/rate, mole.findrate("H,H3+=>H2,H2+")/rate, mole.findrate("H2,H3+=>H2,H2+,H")/rate );
		}
		else
		{
			fprintf( ioQQQ, " Create H2+, rate=0\n" );
		}
	}

	if( trace.lgTrace && trace.lgTraceMole )
	{

		double destroy_coll_heavies = 0.;

		for(ChemNuclideList::iterator atom = nuclide_list.begin(); atom != nuclide_list.end(); ++atom )
		{
			if( !(*atom)->lgHasLinkedIon())
				continue;
			long nelem = (*atom)->el()->Z-1;
			if (dense.lgElmtOn[nelem] && nelem > ipHELIUM) 
			{
				char react[32];
				sprintf(react,"H-,%s=>H,%s",mole_global.list[(*atom)->ipMl[1]]->label.c_str(),(*atom)->label().c_str() );
				destroy_coll_heavies += mole.findrate(react);
			}
		}
		fprintf( ioQQQ, " MOLE, Dep Coef, H-:%10.2e H2:%10.2e H2+:%10.2e\n", 
		  hmi.hmidep, hmi.h2dep, hmi.h2pdep );
		fprintf( ioQQQ, "     H- creat: Rad atch%10.3e Induc%10.3e bHneut%10.2e 3bod%10.2e b=>H2%10.2e N(H-);%10.2e b(H-);%10.2e\n", 
						 mole.findrate("H,e-=>H-,PHOTON"), hmi.HMinus_induc_rec_rate*findspecieslocal("H-")->den, mole.findrate("H,H=>H-,H+"), 
						 mole.findrate("H,e-,e-=>H-,e-"), 
						 mole.findrate("H2,e-=>H,H-"), findspecieslocal("H-")->den, hmi.hmidep );

		fprintf( ioQQQ, "     H- destr: Photo;%10.3e mut neut%10.2e e- coll ion%10.2e =>H2%10.2e x-ray%10.2e p+H-%10.2e\n", 
						 mole.findrate("H-,PHOTON=>H,e-"), destroy_coll_heavies, 
						 mole.findrate("H-,e-=>H-,e-,e-"), 
						 (mole.findrate("H,H-=>H2,e-")+mole.findrate("H,H-=>H2*,e-")), iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].gamnc*findspecieslocal("H-")->den, 
						 mole.findrate("H-,H+=>H,H") );
		fprintf( ioQQQ, "     H- heating:%10.3e Ind cooling%10.2e Spon cooling%10.2e\n", 
						 hmi.hmihet, hmi.HMinus_induc_rec_cooling, hmi.hmicol );
	}

	/* identify creation and destruction processes for H2+ */
	if( trace.lgTrace && trace.lgTraceMole )
	{
		rate = findspecieslocal("H2+")->snk;
		if( rate != 0. )
		{
			rate *= findspecieslocal("H2+")->den;
			if( rate > 0. )
			{
				fprintf( ioQQQ, 
								 " Destroy H2+: rate=%10.2e e-;%5.3f phot;%5.3f hard gam;%5.3f H2col;%5.3f h2phhp;%5.3f pion;%5.3f bh2h2p:%5.3f\n", 
								 rate, mole.findrate("H2+,e-=>H,H+,e-")/rate, mole.findrate("H2+,PHOTON=>H,H+")/rate, 
								 mole.findrate("H2+,CRPHOT=>H,H+")/rate, mole.findrate("H2,H2+=>H,H3+")/rate, mole.findrate("H2+,H2=>H,H+,H2")/rate, 
								 mole.findrate("H2+,H+=>H,H+,H+")/rate, mole.findrate("H,H2+=>H+,H2*")/rate );

				fprintf( ioQQQ, 
								 " Create  H2+: rate=%.2e HII HI;%.3f Col H2;%.3f HII H2;%.3f HI HI;%.3f\n", 
								 rate, 
								 mole.findrate("H+,H=>H2+,PHOTON")/rate, 
								 mole.findrate("H2,CRPHOT=>H2+,e-")/rate, 
								 mole.findrate("H2,H+=>H2+,H")/rate, 
								 mole.findrate("H,H+,e-=>H2+,e-")/rate );
			}
			else
			{
				fprintf( ioQQQ, " Create  H2+: rate= is zero\n" );
			}
		}
	}

	{
		/* following should be set true to print populations */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"mole bugg\t%.3f\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
				fnzone,
				iso_sp[ipH_LIKE][ipHYDROGEN].fb[0].gamnc,
				hmi.HMinus_photo_rate,
				findspecieslocal("H2")->den , 
				findspecieslocal("H-")->den ,
				findspecieslocal("H+")->den);
		}
	}

	{
		/*@-redef@*/
		/* this debug print statement compares H2 formation through grn vs H- */
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && nzone>140 )
		{
			fprintf(ioQQQ," debuggggrn grn\t%.2f\t%.3e\t%.3e\tfrac\t%.3e\tH-\t%.3e\t%.3e\tfrac\t%.3e\t%.3e\t%.3e\t%.3e\n",
				fnzone ,
				mole.findrate("H,H,grn=>H2*,grn")+mole.findrate("H,H,grn=>H2,grn") , 
				hmi.H2_forms_grains+hmi.H2star_forms_grains ,
							mole.findrate("H,H,grn=>H2*,grn")/SDIV(mole.findrate("H,H,grn=>H2*,grn")+mole.findrate("H,H,grn=>H2,grn")),
							(mole.findrate("H,H-=>H2,e-")+mole.findrate("H,H-=>H2*,e-")),
				hmi.H2star_forms_hminus+hmi.H2_forms_hminus,
							frac_H2star_hminus(),
							(mole.findrate("H,H-=>H2,e-")+mole.findrate("H,H-=>H2*,e-")),findspecieslocal("H")->den,findspecieslocal("H-")->den
				);
		}
	}

	{
		/*@-redef@*/
		/* often the H- route is the most efficient formation mechanism for H2,
		 * will be through rate called mole_global.list[unresolved_element_list[ipHYDROGEN]->ipMl]->den*hmi.assoc_detach
		 * this debug print statement is to trace h2 oscillations */
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && nzone>140/*&& iteration > 1*/)
		{
			/* rapid increase in H2 density caused by rapid increase in hmi.rel_pop_LTE_H2g */
			fprintf(ioQQQ,
							"hmi.assoc_detach_backwards_grnd\t%.2f\t%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", 
							/* total forward rate */
							fnzone,
							phycon.te, 
							dense.eden,
							(mole.findrate("H,H-=>H2,e-")+mole.findrate("H,H-=>H2*,e-")),
							mole.findrate("H2,e-=>H,H-"), 
							mole.findrate("H2*,e-=>H,H-"), 
							mole.findrate("H,e-=>H-,PHOTON"),
							hmi.HMinus_induc_rec_rate*findspecieslocal("H-")->den,
					      mole.species[findnuclide("H")->ipMl[0]].den,
							mole.species[findnuclide("H")->ipMl[1]].den,
							findspecieslocal("H-")->den,
							hmi.H2_total,
							hmi.rel_pop_LTE_Hmin,
							hmi.rel_pop_LTE_H2g,
							hmi.rel_pop_LTE_H2s
				);
		}
	}

	if( (trace.lgTrace && trace.lgTraceMole) )
	{
		if( hmi.H2_rate_destroy > SMALLFLOAT )
		{
			fprintf( ioQQQ, 
			  " H2 destroy rate=%.2e DIS;%.3f bat;%.3f h2dis;%.3f photoionize_rate;%.3f h2h2p;%.3f E-h;%.3f hmi.h2hph3p;%.3f sec;%.3f\n", 
							 hmi.H2_rate_destroy, 
							 hmi.H2_Solomon_dissoc_rate_used_H2g / hmi.H2_rate_destroy, 
							 mole.findrate("H2,e-=>H,H-") / (hmi.H2_total*hmi.H2_rate_destroy), 
							 mole.findrate("H,H2=>H,H,H") / (hmi.H2_total*hmi.H2_rate_destroy), 
							 h2.photoionize_rate / hmi.H2_rate_destroy, 
							 mole.findrate("H2,H+=>H2+,H") / (hmi.H2_total* hmi.H2_rate_destroy), 
							 (mole.findrate("H2,e-=>H,H,e-")+mole.findrate("H2*,e-=>H,H,e-")) / (hmi.H2_total* hmi.H2_rate_destroy), 
							 mole.findrate("H2,H+=>H3+") / (hmi.H2_total* hmi.H2_rate_destroy) ,
							 secondaries.csupra[ipHYDROGEN][0]*2.02 / hmi.H2_rate_destroy
			  );
		}
		else
		{
			fprintf( ioQQQ, " Destroy H2: rate=0\n" );
		}
	}

	{
		/* following should be set true to print populations */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC )
		{
			if( DEBUG_LOC && (nzone > 570) ) 
			{
				fprintf(ioQQQ,"Temperature %g\n",phycon.te);
				fprintf(ioQQQ," Net mol ion rate [%g %g] %g\n",mole.source[ipHYDROGEN][1],mole.sink[ipHYDROGEN][1],
						  mole.source[ipHYDROGEN][1]-mole.sink[ipHYDROGEN][1]*mole.species[findnuclide("H")->ipMl[1]].den);
			}
		}
	}
}

STATIC void mole_update_limiting_reactants()
{
	DEBUG_ENTRY( "mole_update_limiting_reactants()" );

	/* largest fraction of atoms in molecules */
	for( long i=0; i<mole_global.num_calc; ++i )
	{
		mole.species[i].xFracLim = 0.;
		mole.species[i].atomLim = null_nuclide;
		for (molecule::nNucsMap::iterator atom = mole_global.list[i]->nNuclide.begin();
			  atom != mole_global.list[i]->nNuclide.end(); ++atom)
		{
			long nelem = atom->first->el()->Z-1;
			if( dense.lgElmtOn[nelem])
			{
				realnum densAtomInSpecies = (realnum)( mole.species[i].den * atom->second );
				realnum densAtom = dense.gas_phase[nelem] * atom->first->frac; 
				if( mole.species[i].xFracLim * densAtom < densAtomInSpecies )
				{
					mole.species[i].xFracLim = densAtomInSpecies/densAtom;
					mole.species[i].atomLim = atom->first.get();
				}
			}
		}
		//if( mole.species[i].atomLim != null_nuclide )
		//	fprintf( ioQQQ, "DEBUGGG mole_update_limiting_reactants %li\t%s\t%s\t%.12e\t%.12e\n", 
		//		i, mole_global.list[i]->label.c_str(), mole.species[i].atomLim->label().c_str(), mole.species[i].xFracLim, mole.species[i].den );
	}
	
	return;
}

