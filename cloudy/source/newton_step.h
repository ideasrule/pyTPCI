/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef NEWTON_STEP_H_
#define NEWTON_STEP_H_

class GroupMap;

/* mole_co_priv.h */
/** Take one Newton step
\param *nBad
\param *error
*/
bool newton_step(GroupMap &mole_map, const valarray<double> &oldmols, valarray<double> &newmols, realnum *eqerror, realnum *error, const long n, 
					  double *rlimit, double *rmax, 
					  valarray<double> &escale,
					  void (*jacobn)(GroupMap &mole_map, 
							const valarray<double> &b2vec,
							double * const ervals, double * const amat,
										  const bool lgJac, bool *lgConserved));

typedef void (*error_print_t)(long, long, 
															const valarray<double> &, 
															const valarray<double> &);

int32 solve_system(const valarray<double> &a, valarray<double> &b, 
									 long int n, error_print_t error_print);

#endif /* NEWTON_STEP_H_ */
