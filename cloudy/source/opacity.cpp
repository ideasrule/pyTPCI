/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "opacity.h"
t_opac opac;

void t_opac::zero()
{
	DEBUG_ENTRY( "t_opac::zero()" );
	/* do we want to save negative opacities */
	lgNegOpacIO = false;

	otsmin = 0.;

	/* this flag says to use the static opacities,
	 * only evaluate them at start of each new zone.
	 * when set false with 
	 * no static opacities
	 * command, always reevaluate them */
	lgOpacStatic =  true;

	/* set true in radinc if negative opacities ever occur */
	lgOpacNeg = false;

	/* can turn of scattering opacities for some geometries */
	lgScatON = true;

	/* variables having to do with compiling and/or using the
	 * ancillary file of stored opacities */
	lgCompileOpac = false;
	/* "no file opacity" command sets following var false, says not to use file */
	/** \todo	2	file opacities are disabled for now - reinstate this when arrays settle down */
	lgUseFileOpac = false;
	stimax[0] = 0.;
	stimax[1] = 0.;
}
