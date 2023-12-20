/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef RT_ESCPROB_H_
#define RT_ESCPROB_H_

class TransitionProxy;

class DestType
{
public:
	enum dest_t {
		ipINVALID = -1,
		/* core function for K2 destruction */
		ipDEST_K2 = 1,
		/* core function for complete redist destruction */
		ipDEST_INCOM,
		/* core function for simple destruction */
		ipDEST_SIMPL,
		/* Lyman alpha destruction */
		ipDEST_LYA
	};
	enum dest_t t;
	double dest;
DestType() : t(ipINVALID), dest(0.0) {}
};

/**esc_PRD_1side fundamental escape probability radiative transfer routine for incomplete redistribution 
\param tau
\param a
*/
double esc_PRD_1side(double tau, 
							double a);

/**esc_CRDwing_1side fundamental escape probability radiative transfer routine, for complete redistribution */
double esc_CRDwing_1side(double tau, 
								 double a );

/**RTesc_lya escape prob for hydrogen atom Lya, using Hummer and Kunasz results 
\param *esin
\param *dest
\param abund
\param t line structure
\param DopplerWidth
*/
double RTesc_lya(
	/* the inward escape probability */
	double *esin, 
	/* the destruction probility */
	double *dest, 
	/* abundance of the species */
	double abund, 
	const TransitionProxy& t,
	realnum DopplerWidth);

/**esc_PRD escape probability radiative transfer for incomplete redistribution 
\param tau
\param tout
\param damp
*/
double esc_PRD(
	double tau, 
	double tau_out, 
	double damp );

/**esc_CRDwing escape probability CRD with wings, for subordinate lines 
\param tau
\param tout
\param damp 
*/
double esc_CRDwing(
	double tau_in, 
	double tau_out, 
	double damp);

/**esc_CRDcore escape probability CRD with no wings, for subordinate lines 
\param tau
\param tout
*/
double esc_CRDcore(
	double tau_in, 
	double tau_out);

/**esca0k2 derive Hummer's K2 escape probability for Doppler core only 
\param taume
*/
double esca0k2(double taume);

/**escpcn continuum escape probability 
\param tau
\param hnukt
*/
double esccon(double tau, 
				  double hnukt);

/**RT_DestProb returns line destruction probability due to continuum opacity 
\param abund abundance of species
\param crsec its line absorption cross section
\param ipanu pointer to energy within continuum array, to get background opacity,
 this is on the f not c scale
\param widl line width
\param escp escape probability
\param nCore type of redistribution function
*/
void RT_DestProb(
	const TransitionProxy& t,
	/* line width */
	double widl, 
	/* type of redistribution function */
	const DestType& nCore);

/**RT_LineWidth compute line width (cm/sec), using optical depth array information 
\param t
\param DopplerWidth
*/
double RT_LineWidth(const TransitionProxy& t, realnum DopplerWidth);
#endif /* RT_ESCPROB_H_ */
