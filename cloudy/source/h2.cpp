/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "h2.h"
#include "hmi.h"
#include "atmdat.h"
	
vector<diatomics*> diatoms;

diatomics h2("h2", 4100., &hmi.H2_total, Yan_H2_CS);
diatomics hd("hd", 4100., &hmi.HD_total, Yan_H2_CS);

diatomics::diatomics( const string& a, const double& e_star, const double* const abund, double (*fun)(double) ):
	trans(a, &states), ENERGY_H2_STAR (e_star), dense_total(abund)
	{
		DEBUG_ENTRY( "diatomics::diatomics()" );
		fixit("should probably force path lower and label upper case.");
		path = a;
		label = a;
		{
			unsigned j;
			for( j = 0; j < a.size(); ++j )
				label[j] = toupper( label[j] );
			shortlabel = label;
			for( j = a.size(); j < 4; ++j )
				label += ' ';
		}
		/* option to turn off or on gbar collisions of the collision rate,
		* default is to have it on */
		/* turn h2.lgColl_gbar on/off with atom h2 gbar on off */
		lgColl_gbar = true;
		/* include collision rates that come from real calculations,
		* off with atom h2 collisions off command */
		lgColl_deexec_Calc = true;
		lgColl_dissoc_coll = true;
		lgEnabled = false;
		lgEvaluated = false;
		/* option to turn off H2 - grain collision & deexcitation,
		* atom h2 grain collision on/off */
		lgH2_grain_deexcitation = false;
		/* option to scramble collision data */
		lgH2_NOISE = false;
		lgH2_NOISECOSMIC = false; 
		/* option to turn off ortho-para collisions, command SPECIES H2 COLLISIONS ORTHO PARA OFF */
		lgH2_ortho_para_coll_on = true;
		/* which set of H2 - He collision data to use?  default is ORNL data set
		 * also Le Bourlot Meudon set available,
		 * set to ORNL with command
		 * atom H2 He collisions ORNL */
		lgH2_He_ORNL = true;
		/*>>chng 08 feb 27, GS, ORNL or Meudon H2 - H2 collision data 
		 * >>chng 09 may 11, make it the default */
		lgH2_ORH2_ORNL = true;
		lgH2_PAH2_ORNL = true;
		/* flag to force using LTE level populations */
		lgLTE = false;
		lgREAD_DATA = false;
		loop_h2_oscil = -1;
		HeatDiss = 0.;
		HeatDexc = 0.;
		HeatDexc_old = 0.;
		HeatDexc_deriv = 0.;
		HeatChangeOld = 0.;
		HeatChange = 0.;
		photo_heat_soft = 0.;
		photo_heat_hard = 0.;
		photodissoc_BigH2_H2s = 0.;
		photodissoc_BigH2_H2g = 0.;
		photoion_opacity_fun = fun;
		Solomon_dissoc_rate_s = 0.;
		Solomon_dissoc_rate_g = 0.;
		spon_diss_tot = 0.;
		rate_grain_op_conserve = 0.;
		rate_grain_J1_to_J0 = 0.;
		Average_A = 0.;
		Average_collH2_deexcit = 0.;
		Average_collH_deexcit = 0.;
		Average_collH2_excit = 0.;
		Average_collH_excit = 0.;
		Average_collH_dissoc_g = 0.;
		Average_collH_dissoc_s = 0.;
		Average_collH2_dissoc_g = 0.;
		Average_collH2_dissoc_s = 0.;

		/* these remember the largest and smallest factors needed to
		 * renormalize the H2 chemistry */
		renorm_max = 1.;
		renorm_min = 1.;
		/* ortho and para column densities */
		ortho_colden = 0.;
		para_colden = 0.;
		ortho_para_old = 0.;
		ortho_para_older = 0.;
		ortho_para_current = 0.;
		TeUsedBoltz = -1.;
		TeUsedColl = -1.;

		/* counters used by H2_itrzn to find number of calls of h2 per zone */
		nH2_pops  = 0;
		nH2_zone = 0;
		/** \todo	3	these should be const since cannot change, are flags */
		/* these are used to set trace levels of output by options on atom h2 trace command 
		 * first is minimum level of trace, keyword is FINAL */
		n_trace_final = 1;
		/* intermediate level of trace, info per iteration, key ITERATION */
		n_trace_iterations = 2;
		/* full trace, keyword is FULL */
		n_trace_full = 3;
		/* print matrices used in solving X */
		n_trace_matrix = 4;
		/* holds options for trace set with atom h2 command */
		nTRACE = false;
		/* this is number of electronic levels to include in the print and save output 
		* changed with the PRINT LINE H2 ELECTRONIC and PUNCH H2 commands 
		* by default only include X */
		nElecLevelOutput = 1;
		/* the number of electronic quantum states to include.
		 * To do both Lyman and Werner bands want nelec = 3,
		 * default is to do all bands included */
		n_elec_states = N_ELEC;
		nCall_this_zone = 0;
		/* the number of levels used in the matrix solution
		 * of the level populations - set with atom h2 matrix nlevel,
		 * >>chng 04 oct 05, make default 30 levels 
		 * >>chng 04 dec 23, make default 70 levels */
		nXLevelsMatrix = 70;
		ndim_allocated = 0;

		levelAsEval = -1;
		lgFirst = true;
		nzone_eval = -1;
		nzoneEval = -1;
		iteration_evaluated = -1;

		/* this is used to establish zone number for evaluation of number of levels in matrix */
		nzone_nlevel_set = -1;
		nzoneAsEval = -1;
		/* this is used to establish zone number for evaluation of number of levels in matrix */
		nzone_nlevel_set = 0;

		/* this is the smallest ratio of H2/H where we will bother with the large H2 molecule
		 * this value was chosen so that large mole used at very start of TH85 standard PDR,
		 * NB - this appears in headinfo and must be updated there if changed here */
		/* >>chng 03 jun 02, from 1e-6 to 1e-8 - in orion veil large H2 turned on half way
		 * across, and Solomon process was very fast since all lines optically thin.  correct
		 * result has some shielding, so reset to lower value so that H2 comes on sooner. */
		H2_to_H_limit = 1e-8;
		iterationAsEval = -1; 

		strcpy( chH2ColliderLabels[0] , "H0" );
		strcpy( chH2ColliderLabels[1] , "He" );
		strcpy( chH2ColliderLabels[2] , "H2 o" );
		strcpy( chH2ColliderLabels[3] , "H2 p" );
		strcpy( chH2ColliderLabels[4] , "H+" );
	}

void diatoms_init( void )
{
	DEBUG_ENTRY( "diatoms_init()" );

	diatoms.clear();
	diatoms.push_back( &h2 );
	diatoms.push_back( &hd );

	// molecular hydrogen, H2
	// H
	h2.coll_source[0].magic = 110416;
	/*>>refer	H2	H collision	Lique, F., 2015, MNRAS, 453, 810 */
	h2.coll_source[0].filename = "coll_rates_H_15.dat";
	// He
	h2.coll_source[1].magic = 110416;
	h2.coll_source[1].filename = "coll_rates_He_ORNL.dat";
	// H2 ortho
	h2.coll_source[2].magic = 110416;
	h2.coll_source[2].filename = "coll_rates_H2ortho_ORNL.dat";
	// H2 para	
	h2.coll_source[3].magic = 110416;
	h2.coll_source[3].filename = "coll_rates_H2para_ORNL.dat";
	// proton
	h2.coll_source[4].magic = 110416;
	h2.coll_source[4].filename = "coll_rates_Hp.dat";

	// HD
	// H
	hd.coll_source[0].magic = 110416;
	hd.coll_source[0].filename = "coll_rates_H_07.dat";
	// He
	hd.coll_source[1].magic = 110416;
	hd.coll_source[1].filename = "coll_rates_He_Flower.dat";
	// H2 ortho
	hd.coll_source[2].magic = 110416;
	hd.coll_source[2].filename = "coll_rates_H2ortho_Flower.dat";
	// H2 para	
	hd.coll_source[3].magic = 110416;
	hd.coll_source[3].filename = "coll_rates_H2para_Flower.dat";
	// proton
	hd.coll_source[4].magic = 110416;
	hd.coll_source[4].filename = "coll_rates_Hp.dat";

	return;
}

