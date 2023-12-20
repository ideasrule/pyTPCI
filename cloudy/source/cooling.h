/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef COOLING_H_
#define COOLING_H_

const bool lgConvBaseHeatTest=false;

/**CoolZero set cooling and heating stack to zero */
void CoolZero(void);

/**CoolAdd add coolants to the cooling stack, called in evaluation of cooling function 
\param *chLabel
\param xlambda
\param cool
*/
void CoolAdd(
  const char *chLabel, 
  realnum xlambda, 
  double cool);

/**CoolSum  total cooling from all entries into cooling stack */
void CoolSum(double *total);

/**CoolEvaluate main routine to call others, to evaluate total cooling 
\param tot total cooling */
void CoolEvaluate(double *tot);

/**coolpr stores coolants before block printed, when printing cooling agents 
\param *io the label for the coolant
\param *chLabel
\param lambda the wavelength
\param ratio the ratio of this coolant, to total cooling, may be negative
\param *chJOB which job, either ZERO, DOIT, or DONE 
*/
void coolpr(
	FILE *io,
	const char *chLabel ,
	realnum lambda,
	double ratio,
	const char *chJOB );

/**HeatSum evaluate all heating agents to determine total heating for this zone,
 * called at end of ionize */
void HeatSum(void);

void SecIoniz( void );

/**HeatZero zeroes out the heating array, called at start of ionize*/
void HeatZero(void);

void CoolDima(void);

#endif /* COOLING_H_ */
