/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "hmi.h"
t_hmi hmi;

void t_hmi::zero()
{
	DEBUG_ENTRY( "t_hmi::zero()" );
	H2_total = 0.;
	H2_total_f = 0.f;
	HD_total = 0.;
	H2_frac_abund_set = 0.;
	hmihet = 0.;
	h2plus_exc_frac = 0.;
	h2plus_heatcoef = 0.;
	h2plus_heat = 0.;
	HeatH2Dish_used = 0.;
	HeatH2Dexc_used = 0.;
	HeatH2Dish_TH85 = 0.;
	HeatH2Dexc_TH85 = 0.;
	UV_Cont_rel2_Draine_DB96_face = 0.;
	UV_Cont_rel2_Draine_DB96_depth = 0.;
	UV_Cont_rel2_Habing_TH85_face = 0.;
	UV_Cont_rel2_Habing_TH85_depth = 0.;
	HeatH2DexcMax = 0.;
	CoolH2DexcMax = 0.;
	hmitot = 0.;
	H2Opacity = 0.;
	hmidep = 1.;
	h2dep = 1.;
	h2pdep = 1.;
	h3pdep = 1.;

	/* option to kill effects of H2 in CO chemistry - set with
	 * set Leiden hack h2* off */
	lgLeiden_Keep_ipMH2s = true;
	lgLeidenCRHack = true;

	/* this says which estimate of the rate of the Solomon process to use,
	 * default is Tielens & Hollenbach 1985a, changed with 
	 * set h2 Solomon command, options are TH85 and BD96,
	 * the second for the Bertoldi & Draine rates - they
	 * differ by 1 dex.   when large H2 turned on this is ignored */
	/* the Tielens & Hollenbach 1985 treatment */
	chH2_small_model_type = 'T';
	/* the improved H2 formalism given by 
	*>>refer	H2	dissoc	Burton, M.G., Hollenbach, D.J. & Tielens, A.G.G.M 
	>>refcon	1990, ApJ, 365, 620 */
	chH2_small_model_type = 'H';
	/* the Bertoldi & Draine 1996 treatment */
	/* >>chng 03 nov 15, change default to BD96 */
	chH2_small_model_type = 'B';
	/* >>chng 05 dec 08, use the Elwert et al. approximations as the default */
	chH2_small_model_type = 'E';

	/* set NaN */
	set_NaN( HeatH2Dish_used ); 
	set_NaN( HeatH2Dish_TH85 );
	set_NaN( HeatH2Dish_BD96 );
	set_NaN( HeatH2Dish_BHT90 );
	set_NaN( HeatH2Dish_ELWERT );

	/** HeatH2Dexc_used is heating due to collisional deexcitation of vib-excited 
	* H2 actually used */
	set_NaN( HeatH2Dexc_used );
	set_NaN( HeatH2Dexc_TH85 );
	set_NaN( HeatH2Dexc_BD96 );
	set_NaN( HeatH2Dexc_BHT90 );
	set_NaN( HeatH2Dexc_ELWERT );

	/** these are derivative wrt temp for collisional processes within X */
	set_NaN( deriv_HeatH2Dexc_used );
	set_NaN( deriv_HeatH2Dexc_TH85 );
	set_NaN( deriv_HeatH2Dexc_BD96 );
	set_NaN( deriv_HeatH2Dexc_BHT90 );
	set_NaN( deriv_HeatH2Dexc_ELWERT );

	set_NaN( H2_Solomon_dissoc_rate_used_H2g );
	set_NaN( H2_Solomon_dissoc_rate_TH85_H2g );
	set_NaN( H2_Solomon_dissoc_rate_BHT90_H2g );
	set_NaN( H2_Solomon_dissoc_rate_BD96_H2g );
	set_NaN( H2_Solomon_dissoc_rate_ELWERT_H2g );

	set_NaN( H2_Solomon_dissoc_rate_used_H2s );
	set_NaN( H2_Solomon_dissoc_rate_TH85_H2s );
	set_NaN( H2_Solomon_dissoc_rate_BHT90_H2s );
	set_NaN( H2_Solomon_dissoc_rate_BD96_H2s );
	set_NaN( H2_Solomon_dissoc_rate_ELWERT_H2s );

	/** the Solomon process rate H2 dissociates into X continuum - actually used */
	/**set_NaN( H2_Solomon_dissoc_rate_used );*/
	/** H2 + hnu => 2H from TH85 */
	/** H2 + hnu => 2H actually used */
	set_NaN( H2_photodissoc_used_H2g );
	set_NaN( H2_photodissoc_used_H2s );
	set_NaN( H2_photodissoc_ELWERT_H2g );
	set_NaN( H2_photodissoc_ELWERT_H2s );
	set_NaN( H2_photodissoc_TH85 );
	set_NaN( H2_photodissoc_BHT90 );

	/* default grain formation pumping - Takahashi 2001 */
	chGrainFormPump = 'T';

	/* set which approximation for Jura rate - Cazaux & Tielens 
	 * >>refer	H2	form	Cazaux, S., & Tielens, A.G.G.M., 2002, ApJ, 575, L29 */
	/*>> chng 14 aug 28, to Eley Rideal */
	chJura = 'E';

	/* scale factor to multiply Jura rate, set Jura rate command */
	ScaleJura = 1.f;

	/* binding energy for change in H2 population while on grain surface,
	 * set with "set h2 Tad" command */
	Tad = 800.;

	lgH2_Thermal_BigH2 = true;
	lgH2_Chemistry_BigH2 = true;
}
