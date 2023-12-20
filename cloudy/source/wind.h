/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef WIND_H_
#define WIND_H_

/** parameters for wind model; wind.h */
struct Wind {

	/** initial wind velocity cm s-1, will be negative for D-critical flow */
	realnum windv0;

	/** central object mass in solar units, to compute inward acceleration */
	realnum comass;

	/** current wind velocity cm s-1, code tests for wind solution by
	 * seeing whether this is positive, static solution if 0 */
	realnum windv;

	/** acceleration within local zone */
	realnum dvdr;

	/** test for static model */
	bool lgStatic(void) const
		{
			ASSERT( ( windv0 == 0.) == m_lgStatic);
			return m_lgStatic;
		}

	/** test for model in ballistic approximation */
	bool lgBallistic(void) const
		{
			return m_lgBallistic ;
		}

	/** mass flux at inner radius, n*v [cm-3 cm s-1 == cm-2 s-1], 
	 * set in parsecommands at inner boundary.  needs r^2 to be real mass flux 
	 * density set by conservation of this quantity in wind model */
	realnum emdot;

	/** *flag set if wind not ok, (below sonic point) */
	bool lgWindOK;

	/** variables for tracking average radiative acceleration - used to report
	 * average acceleration at end of calculation */
	realnum AccelAver, acldr;

	/** inward gravitational acceleration, cm s-2, a positive number */
	realnum AccelGravity;

	/** total outward radiative acceleration, lines and continuum, cm s-2 */
	realnum AccelTotalOutward;

	/** continuum radiative acceleration - cm s-2 */
	realnum AccelCont;

	/** electron scattering radiative acceleration - cm s-2 */
	realnum AccelElectron;

	/** line radiative acceleration - cm s-2 */
	realnum AccelLine;

	/** force multiplier - dimensionless ratio of total to electron scattering 
	 * acceleration */
	realnum fmul;

	/** largest acceleration, usually at illuminated face - cm s-2 */
	realnum AccelMax;

	/** flag set if wind solution had negative velocity */
	bool lgVelPos;

	/* option to do rotating disk, set true with DISK option on wind command */
	bool lgDisk;

	/* for rotating disk this is inner radius, else it is zero,
	 * divide by this to get centrifugal acceleration */
	double DiskRadius;

	Wind(void) : windv0(0.), m_lgStatic(true), m_lgBallistic(false) {}

	void setDefault(void)
		{
			m_lgStatic = false;
			m_lgBallistic = false;
		}

	void setStatic(void)
		{
			m_lgStatic = true;
			m_lgBallistic = false;
		}

	void setBallistic(void)
		{
			m_lgStatic = false;
			m_lgBallistic = true;
		}

	private:
	bool m_lgStatic;

	bool m_lgBallistic;

	};

extern Wind wind;


#endif /* WIND_H_ */
