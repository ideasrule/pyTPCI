/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef COLDEN_H_
#define COLDEN_H_

/* colden.h */

#include "module.h"

enum {
	/** total hydrogen column density, all forms, xxx + 2*H2 */
	ipCOL_HTOT,
	/** column density in electrons */
	ipCOL_elec,
	/** number of column densities to remember */
	NCOLD};

struct t_colden : public module {

	const char *chName() const
	{
		return "colden";
	}
	void zero();
	void comment(t_warnings&) {}

	/**save total column densities in various species for this and 
	 * previous iteration, to check whether
	 * it changed by too much (a bad sign)
	 * column densities, mostly h molecule related */
	realnum colden[NCOLD], 
		/** the previous iteration's coluumn density */
	  colden_old[NCOLD];

	/** integral of n(H2) / v(H2) */
	realnum coldenH2_ov_vel;

	/** integral of ne np over radius */
	double dlnenp;

	/** integral of ne n(He+) over radius */
	double dlnenHep;

	/** integral of ne n(He++) over radius */
	double dlnenHepp;

	/** integral of ne n(C+) over radius */
	double dlnenCp;

	/** integral of n(H0) / Tspin */
	double H0_ov_Tspin;

	/** integral of n(OH) / Tkin */
	double OH_ov_Tspin; 

	/** variables to do with mean mass per particle over model,
	 * called \<Mol\> in final print out, used to get Jeans mass */
	realnum rjnmin, 
	  ajmmin;
	realnum TotMassColl, 
	  tmas, 
	  wmas;

	/** column densities in the lower and upper level of 1s of H0 */
	double H0_21cm_upper;
	double H0_21cm_lower;

	};
extern t_colden colden;

#endif /* COLDEN_H_ */
