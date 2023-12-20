/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef NUMDERIV_H_
#define NUMDERIV_H_

#include "module.h"

/* numderiv.h */
struct t_NumDeriv : public module {
	const char *chName() const
	{
		return "NumDeriv";
	}
	void zero();
	void comment(t_warnings&) {}

	bool lgNumDeriv;
	};
extern t_NumDeriv NumDeriv;


#endif /* NUMDERIV_H_ */
