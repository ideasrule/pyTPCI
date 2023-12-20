/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef GAMMAS_H_
#define GAMMAS_H_

struct t_phoHeat;

/**\file gammas.h
 * gammas.h - all the routines to evaluate gamma functions <BR>
 * contains the following:<BR>
 *
 * GammaBn evaluate photoionization rate for single shell with induced recomb<BR>
 * GammaPrt special version of gamma function to print strong contributors <BR>
 * GammaK evaluate photoionization rate for single shell <BR>
 */

/**GammaBn evaluate photoionization rate for single shell with induced recomb 
\param n1
\param n2
\param ip
\param thresh
\param *ainduc
\param *rcool
\param *photoHeat
*/
double GammaBn(long int n1, 
  long int n2, 
  long int ip, 
  double thresh, 
  double *ainduc, 
  double *rcool,
  t_phoHeat *photoHeat);

/**GammaPrtShells for the element nelem and ion, print total photo rate, subshells,
 * and call GamaPrt for important subshells 
 \param nelem
 \param ion
 */
void GammaPrtShells( long nelem , long int );

/**GammaPrt special version of gamma function to print strong contributors 
\param n1
\param n2
\param ip
\param io io unit we will write to 
\param total
\param threshold
*/
void GammaPrt(long int n1, 
	  long int n2, 
	  long int ip, 
	  FILE * io, 
	  double total, 
	  double threshold);

/**GammaK evaluate photoionization rate for single shell */
double GammaK(long int n1, 
  long int n2, 
  long int ipOpac, 
  double yield1,
  t_phoHeat *photoHeat);

/** GammaPrtRate will print resulting rates for ion and element 
\param *ioFILE io unit we will write to
\param ion stage of ionization on C scale, 0 for atom
\param ipZ 0 for H, etc 
\param lgPRT true - then print photo sources for valence shell
*/
void GammaPrtRate(FILE * ioFILE,
	long int ion ,
	long int ipZ,
	bool lgPRT );

#endif /* GAMMAS_H_ */
