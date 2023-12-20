/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ITERATIONS_H_
#define ITERATIONS_H_

#include "module.h"

/**IterStart, set and save values of many variables at start of iteration */
void IterStart(void);

/**IterRestart, reset values of many variables at start of iteration */
void IterRestart(void);

/** close out this iteration */
void IterEnd(void);

/**iter_end_check called by Cloudy after each zone to determine whether iteration is complete
 * returns true if iteration is complete, false if not */
int iter_end_check(void);

struct t_iterations : public module {

	const char *chName() const
	{
		return "iterations";
	}
	void zero();
	void comment(t_warnings&) {}

	void alloc();

	/** these are the variables that control how many iterations are
	 * to be done, and number of the current iteration
	 * itermx is number of iterations to perform, set with iterate command
	 * upper limit is parameter variable ItrDim */
	long int itermx;

	/** amount of space that has been allocated for max iterations */
	long int iter_alloc;

	/**number of zones to print on each iteration*/
	vector<long int> IterPrnt/**[ITR DIM]*/;

	/** this is false on any but the last iteration
	 * set true in startr if iter > itermx */
	bool lgLastIt;

	/** flag indicating that another iteration is needed, includes
	 * various criteria, done in prt_comment.cpp for ALL models.
	 * This does check some line optical depths.
	 * comment included in final printout*/
	bool lgIterAgain;

	/* has the optical depth scale converged, or is another iteration needed?
	 * These checks are only done when ITERATE TO CONVERGENCE is requested */
	bool lgOpticalDepthonverged;

	/** next three implement set coverage command to limit iterations and zones */
	bool lgConverge_set;
	long int lim_zone;
	long int lim_iter;

	/** lgEndDflt true if still at default */
	bool lgEndDflt;

	/** default limit to number of zones */
	long int nEndDflt;

	/** nend[iteration] is limiting number of zones for that iteration */
	vector<long int> nend /**<[ITR DIM]*/;

	/** total physical thickness of modeled region, (NOT OUTER RADIUS)
	 * this can set set as a stopping criteria, but if
	 * not set is 1e30 before first iteration.  At end of each iteration, thickness
	 * is set to total depth from illuminated face to outer edge */
	vector<double> StopThickness/**[ITR DIM]*/;

	/** stopping radius for the model, set with STOP RADIUS command */
	vector<double> StopRadius/**[ITR DIM]*/;
	
};
extern t_iterations iterations;

#endif /* ITERATIONS_H_ */
