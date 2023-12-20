/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CO_H_
#define CO_H_

/* co.h */

#include "module.h"

struct t_co : public module {
	const char *chName() const
	{
		return "co";
	}
	
	void zero();
	void comment(t_warnings&) {}

	/** CODissHeat is CO Photodissociation heating */
	realnum CODissHeat, 
	  /**< largest fraction of total heating */
	  codfrc, 
	  /**< total heating integrated over cloud */
	  codtot;
};

extern t_co co;

#endif /* CO_H_ */
