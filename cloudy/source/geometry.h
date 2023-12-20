/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

/** geometry.h */

#include "module.h"

struct t_geometry : public module {
	const char *chName() const
	{
		return "geometry";
	}
	
	void zero();
	void comment(t_warnings&) {}

	/** true if geometry is plane parallel, set in comment */
	bool lgGeoPP;

	/** 1/cos incident continuum entering AngleIllumRadian
	 * from normal - this is 1./cos(normal), so > 1. */
	realnum DirectionalCosin;

	/** FillFac is filling factor, default 1,
	* others are parameters to make function of radius */
	realnum FillFac, 
	  filpow, 
	  fiscal;

	/** flag saying whether spherical geometry */
	bool lgSphere;

	/** default open geometry, covgeo is unity, covrt is zero.  
	 * closed geometry, covgeo is unity, and covrt is unity */

	/* covering factors, account for possible less than total coverageof 4\pi */ 
	/** linearly affects line luminosities but not intensities (since they are
	 * per unit area of cloud.
	 * When sphere command is entered with no optional covering factor this is
	 * set to 1 and covrt is set to covgeo.  The default is an open geometry,
	 * and covgeo is unity for this as well */
	realnum covgeo;

	/** the aperture covering factor; must be between 0 and 1:
	 * iEmissPower == 2 : in this case the variable is set equal to covgeo
	 * iEmissPower == 1 : the fraction of a large circle in the plane that is being
	 *                    observed that actually passes through nebular material.
	 * iEmissPower == 0 : is gas present only at the front or the back (D_0 = 1/2)
	 *                    or on both sides (D_0 = 1) of the central source.
	 * this is the variable D_alpha in the Hazy discussion of the APERTURE command */
	realnum covaper;

	/** radiative transfer covering factor accounts for diffuse radiation
	 * escaping from system in open geometry.  weakly affects intensities and
	 * luminosities.  For open geometry this is zero (half of diffuse fields
	 * escape from cloud and do not strike gas, for a closed geometry this is unity 
	 * also set to covgeo (def unity) if covering factor command is entered */
	realnum covrt;

	/** flag saying that spherical geometry is static, set with sphere static */
	bool lgStatic;

	/** this is the power of r that is used in the line flux integral
	 * 2 : default case - integrate the whole (spherical) nebula
	 * 1 : simulate a long slit over the central source
	 * 0 : simulate a pencil beam centered on the central source
	 * this is the variable alpha in the Hazy discussion of the APERTURE command */
	long int iEmissPower;

	/** size of the aperture:
	 * iEmissPower == 2 : in this case the variable is not used
	 * iEmissPower == 1 : slit width in arcsec
	 * iEmissPower == 0 : surface area of the pencil beam in arcsec^2 */
	realnum size;

	/** has the aperture size been set by the user? */
	bool lgSizeSet;

	/** flag saying that it is ok to not iterate when sphere static is set,
	 * set with (OK) option on sphere static, used for testing hydrogen atom */
	bool lgStaticNoIt;

	/** nprint is how many zones to print */
	long int  nprint;

	/* the largest number of zones needed in any iteration, used to allocate
	 * arrays that save source function */
	long int nend_max;

	/** lgZonSet set if stopping zone specified */
	bool lgZoneSet;

	/** lgZoneTrp set if stopped due to zone number */
	bool lgZoneTrp;

	};
extern t_geometry geometry;



#endif /* GEOMETRY_H_ */
