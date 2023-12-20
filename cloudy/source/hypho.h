/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HYPHO_H_
#define HYPHO_H_

/**hypho generate hydrogenic photoionization cross sections 
\param zed
\param n
\param lmin
\param lmax
\param en
\param ncell
\param anu[]
\param bfnu[]
*/
void hypho(double zed, 
  long int n, 
  long int lmin, 
  long int lmax, 
  double en, 
  long int ncell, 
  realnum anu[], 
  realnum bfnu[]);

#endif /* HYPHO_H_ */
