/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HMI_H_
#define HMI_H_

#include "module.h"

/**hmirat computes radiative association rate for H- 
\param te
*/
double hmirat(double te);

/** hmi.h - parameters dealing with hydrogen molecules */
struct t_hmi : public module {
	const char *chName() const
	{
		return "hmi";
	}

	void zero();
	void comment(t_warnings&) {}

	/** the total H2 density [cm-3], NOT 2*n(H2), the sum of H2 and H2* */
	double H2_total;
	realnum H2_total_f; // single-precision version of above
	double HD_total;

	/** rate ground of H2 is destroyed */
	double H2_rate_destroy;

	/** hminus heating, free bound */
	double hmihet, 
	  hmitot, 
	  hmicol;

	/** mean cross section (cm-2) for H2 Lyman absorption */
	realnum H2Opacity;

	/** these are departure coef for H-, H2, H2+, and HeH,
	 * defined in hmole */
	double hmidep, 
	  h2dep, 
	  h2pdep, 
	  h3pdep;

	/** heating due to photo dissoc of H2+ */
	double h2plus_heatcoef, h2plus_heat;
	/** fraction of H2+ in excited states */
	double h2plus_exc_frac;

	/** H- photo dissoc rate */
	double HMinus_photo_rate;

	realnum 
	  /** the largest fraction of total heat anywhere in model */
	  HeatH2DexcMax,
	  /** the largest fraction of total cooling anywhere in model */
	  CoolH2DexcMax,
	  h2dfrc, 
	  h2dtot,
	  /** fraqction of cooling carried by H2 lines */
	  h2line_cool_frac;

	double HMinus_induc_rec_cooling,
	  HMinus_induc_rec_rate,
	  HMinus_photo_heat;

	long int iheh1, 
	  iheh2;

	/** UV flux relative to Habing value, used for some simple molecular photodissociation rates,
	 * as defined by Tielens & Hollenbach 1985 */
	realnum UV_Cont_rel2_Habing_TH85_face,
	  UV_Cont_rel2_Habing_TH85_depth,
	  /** the special version of g0 with adjustable bounds */
	  UV_Cont_rel2_Habing_spec_depth;

	/** UV flux relative to Habing value, used for some simple molecular photodissociation rates,
	 * as defined by Draine & Bertoldi 1996 -0 we try to do this the way they describe,
	 * since they say that this will agree with their large H2 molecule, first
	 * define field at the illuminated face, then get value at depth using their
	 * form of the extinction and shielding, rather than our exact calculation */
	realnum UV_Cont_rel2_Draine_DB96_face , 
		UV_Cont_rel2_Draine_DB96_depth;

	/** the Solomon process excitation, H2g -> H2*, rate from Tielens & Hollenbach 85 */
	double H2_H2g_to_H2s_rate_TH85;

	/** the Solomon process excitation, H2g -> H2*, rate from Burton et al. 1990 */
	double H2_H2g_to_H2s_rate_BHT90;
	
	/** the Solomon process excitation, H2g -> H2*, rate for the Bertodi & Draine model */
	double H2_H2g_to_H2s_rate_BD96;
	
	/** the Solomon process excitation, H2g -> H2*, rate for Elwert et al. model in prep.*/
	double H2_H2g_to_H2s_rate_ELWERT;
	
	/** the Solomon process excitation, H2g -> H2*, - actually used */
	double H2_H2g_to_H2s_rate_used;

	/** the Solomon process dissociate rate from Tielens & Hollenbach 85 */
	double H2_Solomon_dissoc_rate_used_H2g;
	double H2_Solomon_dissoc_rate_TH85_H2g;
	double H2_Solomon_dissoc_rate_BHT90_H2g;
	double H2_Solomon_dissoc_rate_BD96_H2g;
	double H2_Solomon_dissoc_rate_ELWERT_H2g;

	double H2_Solomon_dissoc_rate_used_H2s;
	double H2_Solomon_dissoc_rate_TH85_H2s;
	double H2_Solomon_dissoc_rate_BHT90_H2s;
	double H2_Solomon_dissoc_rate_BD96_H2s;
	double H2_Solomon_dissoc_rate_ELWERT_H2s;

	/** the Solomon process rate H2 dissociates into X continuum - actually used */
	/**double H2_Solomon_dissoc_rate_used;*/
	/** H2 + hnu => 2H from TH85 */
	/** H2 + hnu => 2H actually used */
	double H2_photodissoc_used_H2g;
	double H2_photodissoc_used_H2s;
	double H2_photodissoc_ELWERT_H2g;
	double H2_photodissoc_ELWERT_H2s;
	double H2_photodissoc_TH85;
	double H2_photodissoc_BHT90;


	/** continuum array index for H minus threshold  */
	long int iphmin;

	/** largest local fraction heating due to dissoc of H2+ */
	realnum h2pmax;

	/** binding energy for change in H2 population while on grain surface,
	 * set with "set h2 Tad " command */
	realnum Tad;

	double 

	  /** HeatH2Dish_used is heating due to H2 dissociation actually used*/
	  HeatH2Dish_used, 
	  HeatH2Dish_TH85, 
	  HeatH2Dish_BD96 ,
	  HeatH2Dish_BHT90, 
	  HeatH2Dish_ELWERT ,

	  /** HeatH2Dexc_used is heating due to collisional deexcitation of vib-excited 
	   * H2 actually used */
	  HeatH2Dexc_used,
	  HeatH2Dexc_TH85,
	  HeatH2Dexc_BD96,
	  HeatH2Dexc_BHT90,
	  HeatH2Dexc_ELWERT;

	/** these are derivative wrt temp for collisional processes within X */
	realnum 
		deriv_HeatH2Dexc_used,
		deriv_HeatH2Dexc_TH85 ,
		deriv_HeatH2Dexc_BD96 ,
		deriv_HeatH2Dexc_BHT90 ,
		deriv_HeatH2Dexc_ELWERT;

	/** these are the H- and grain formation rates, added above and below a
	 * certain energy (2.6 eV) for production of H2 or H2* in small network */
	double H2_forms_grains ,
		H2_forms_hminus,
		H2star_forms_grains,
		H2star_forms_hminus;

	/** say how to do thermal solution, if true (default) use results of large molecule,
	 * if false use TH85 approximations */
	bool lgH2_Thermal_BigH2,
	/** say how to do chemistry (formation and destruction), 
	 * if true (default) use results of large molecule,
	 * if false use TH85 approximations */
	  lgH2_Chemistry_BigH2;

	/** the set h2 small model command tells code says which of the small model H2
	 * to use.  Default is Elwert */
	char chH2_small_model_type;

	/** method used for grain formation pumping */
	char chGrainFormPump;

	/** the set h2 jura command tells code which treatment of H2 formation to use */
	char chJura;

	/** this is a scale factor to multiply the Jura rate, default is unity, changed
	 * with the set jura scale command */
	realnum ScaleJura;

	/** H2 formation rate as set with set h2 rate command units S^-1, actual depl */
	double rate_h2_form_grains_set;  

	/** this is set to zero, but to positive number with atom h2 fraction command
	 * this sets the H2 density by multiplying the hydrogen density to become the H2 density */
	double H2_frac_abund_set;

	/** Boltzmann factor for hmi */
	double exphmi,
	/** related to the LTE populations of H-, H2, and H2+
	 * each is a constant with temperature dependence, and
	 * needs to be multiplied by the densities of the separated
	 * components to become the LTE density.  
	 * following is n(H-) / [ n(e) n(H) ], units cm3 */
	  rel_pop_LTE_Hmin,
	/** related to the LTE population of H2s, following is  
	 * n(H2s) / [n(H) n(H) ], units cm3 */
	  rel_pop_LTE_H2s;
	/** LTE population for H2+, following is
	 * n(H2+) / [n(H) n(p) ], units cm3 */
	double rel_pop_LTE_H2p,
	/** related to the LTE population of H2 in ground, following is  
	 * n(H2) / [n(H) n(H) ], units cm3 */
	  rel_pop_LTE_H2g,
	  /** related to population of H3+ */
	  rel_pop_LTE_H3p;

	/** hack to kill effects of H2* in chemistry network "set leiden hack h2* off */
	bool lgLeiden_Keep_ipMH2s;
	bool lgLeidenCRHack;

};
extern t_hmi hmi;

#endif /* HMI_H_ */
