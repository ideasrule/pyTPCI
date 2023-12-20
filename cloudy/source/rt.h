/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef RT_H_
#define RT_H_

#include "module.h"

class TransitionProxy;

/* in following two the logical variable says whether to do the
 * escape probabilities too (true) or just the destruction probabilities (false) */

/**RT_line_one do line radiative transfer 
\param t line structure
\param pestrk Stark escape probability
\param DopplerWidth 
\param lgShield_this_zone this is option to not include line self shielding across this zone.
this can cause pump to depend on zone thickness, and leads to unstable
feedback in some models with the large H2 molecule, due to Solomon
process depending on zone thickness and level populations.
*/
void RT_line_one_escape(const TransitionProxy& t, 
		 bool lgShield_this_zone,
		 realnum pestrk,
		 realnum DopplerWidth);
void RT_line_one_fine(const TransitionProxy& t, 
		 bool lgShield_this_zone,
		 realnum pestrk,
		 realnum DopplerWidth);

// Reset the fine opacity array and velocity shift
void RT_fine_clear();

typedef void (*linefunc)(const TransitionProxy& t, 
								 bool lgShield_this_zone,
								 realnum pestrk,
								 realnum DopplerWidth);

/**MakeRT drive static or wind metal line radiative transfer,
 */
void RT_line_all( linefunc line_one );

void RT_line_all_escape( realnum* error ); 

/**rt_line_driving radiative acceleration due to line absorption of incident continuum
 * returns line radiative acceleration cm s-2 */
double RT_line_driving(void);


/**rt_continuum_shield_fcn computing continuum shielding due to single line 
\param t
*/
double RT_continuum_shield_fcn( const TransitionProxy &t, 
										  bool lgShieldThisZone, double dTau );

/**RT_diffuse fill in DiffCont array with diffuse emission for this zone */
void RT_diffuse(void);

/**RT_continuum attenuation of diffuse and beamed continua */
void RT_continuum(void);

/**RT_OTS compute diffuse fields due to helium atom, metals,
 * hydro done in HydroOTS */
void RT_OTS(void);

/**RT_OTS_AddLine add local destruction of lines to ots field 
\param ots
\param ip pointer on the f scale
*/
void RT_OTS_AddLine(double ots, 
  long int ip );

/**RTOTSUpdata sum flux, otscon, otslin, ConInterOut, outlin, 
 * to form SummeDif, SummedCon SummedOcc,  
 * int * is pointer to energy in continuum where this change happened,
 \param SumOTS
  \return sum of ots rates
 */
void RT_OTS_Update(double* SumOTS);

/** zero these things out, called in zero */
void RT_OTS_Zero( void );

/**RT_OTS_ChkSum sanity check confirms summed continua reflect contents of individuals */
void RT_OTS_ChkSum(
	long int ipPnt);

/**RT_line_one_tauinc increment optical depths for all heavy element lines, zone by zone 
\param t
\param mas_species
\param mas_ion
\param mas_hi
\param mas_lo
\param DopplerWidth
*/
void RT_line_one_tauinc(
	const TransitionProxy & t ,
	long int mas_species,
	long int mas_ion,
	long int mas_hi,
	long int mas_lo,
	realnum DopplerWidth);

/**RT_tau_init set initial outward optical depths at start of first iteration */
void RT_tau_init(void);

/**RT_line_one_tau_reset computes average of old and new optical depths for new scale at end of iter 
\param t
*/
class TransitionProxy;
void RT_line_one_tau_reset(
	const TransitionProxy & t );

/**RT_tau_reset update total optical depth scale, called after iteration is complete */
void RT_tau_reset(void);

/**RT_tau_inc increment optical depths once per zone, called after radius_increment */
void RT_tau_inc(void);

/**RT_OTS_PrtRate print ots arrays, called by ionize 
\param weak arg is weakest rate to print
\param chFlag flag, 'c' continuum, 'l' line, 'b' both
*/
void RT_OTS_PrtRate(
	  double weak ,
	  int chFlag );


#if 0
/**wrapper to call RT_LineWidth for the special case of Lya 
\param vth answer will be in whatever units vth is in, since only generates
	optical depth dependent scale factor
*/
double RT_LyaWidth(
			 double tauin, 
			 double tauout, 
			 double a, 
			 double vth);
#endif

/**rt_recom_effic generate escape probability function for continua, 
\param ip
*/
double RT_recom_effic(
	long int ip);

/**rt_stark compute stark broadening escape probabilities using Puetter formalism */
void RT_stark(void);

/** DEST0 is the smallest destruction probability to return
 * in high metallicity models */
/* #define DEST0 1e-8 */
#define DEST0 SMALLFLOAT

struct t_rt : public module {

	const char *chName() const
	{
		return "rt";
	}

	void zero();
	void comment(t_warnings&) {}

	/** wayin - escape probability in inward direction */
	realnum wayin, 

	/** wayout - escape probability in outward direction
	 * =1 when outward optical depths unknown (touton .false.) */
	  wayout; 

	/** fractin = wayin / (wayin+wayout) when outer defined, else zero */
	realnum fracin;

	/** 1 or two, set to double optical depth scale, set to 2 with double command
	* default is 1 */
	realnum DoubleTau;

	/** offset in continuum array for energy where x-ray opacity determined */
	long int ipxry;

	/** optical depth at this energy */
	realnum tauxry;

	/** option to turn off fine structure line optical depths */
	bool lgFstOn;

	/** which type of line continuum shielding function should be used?  */
	int nLineContShield;

	/** include electron scattering escape for lines? */
	bool lgElecScatEscape;

	/** allow photon loss by absorption by background opacities */
	bool lgAbsorLineEscape;

	/** include electron scattering escape for line intensity? */
	bool lgElecScatInten;

	/** dTauMase is smallest maser optical depth in atoms, set in
	 * RT_tau_inc for H, and in tauchn for heavy elements
	 * it is negative or zero */
	realnum dTauMase;

	/** set true in radius_next if maser ever sets zone thickness */
	bool lgMaserSetDR;

	/** flag set true in tauchn if maser cap ever hit, 
	 * causes comment to be printed in prtComments */
	bool lgMaserCapHit;

	/** these identify the species that had a major maser, for debugging */
	long int mas_species , mas_ion , mas_hi , mas_lo;

	/** flag saying that stark broadening is enabled, set false with no stark */
	bool lgStarkON;

};

extern t_rt rt;

double RT_EscLVG( double tau, double sigma );

/** these are all possible values of rt.nLineContShield,
 * first is default, these are set with set continuum shielding */
enum
{
	LINE_CONT_SHIELD_PESC = 1,
	LINE_CONT_SHIELD_FEDERMAN,
	LINE_CONT_SHIELD_FEDERMAN_BUG,
	LINE_CONT_SHIELD_FERLAND,
	LINE_CONT_SHIELD_RODGERS,
	LINE_CONT_SHIELD_INTEGRAL
};

void prt_trans_opc_debug( const char *LineGroup, const TransitionProxy &t );

#endif /* RT_H_ */
