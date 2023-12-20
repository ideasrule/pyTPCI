/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MAGNETIC_H_
#define MAGNETIC_H_

/* magnetic.h */

#include "module.h"

/**ParseMagnet parse magnetic field command  
\param *chCard
*/
class Parser;
void ParseMagnet(Parser &p );

/**Magnetic_reinit - reinitialized magnetic field at start of new iteration */
void Magnetic_reinit(void);

/**Magnetic_evaluate evaluate some parameters to do with magnetic field */
void Magnetic_evaluate(void);

struct t_magnetic : public module {

	const char *chName() const
	{
		return "magnetic";
	}

	/** zero initialize magnetic field parameters */
	void zero(void);
	void comment(t_warnings&) {}

	/** this says mag fields are turned on */
	bool lgB;

	/** enthalpy density associated with field */
	double EnthalpyDensity;

	/** this is magnetic pressure at current location - this can be negative for ordered field */
	double pressure;

	/** energy density at current location, this is positive */
	double energydensity;

	};
extern t_magnetic magnetic;


#endif /* MAGNETIC_H_ */
