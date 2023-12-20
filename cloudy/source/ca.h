/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CA_H_
#define CA_H_

#include "module.h"

/** ca2mly.h */
struct t_ca : public module {
	const char *chName() const
	{
		return "ca";
	}

	void zero();
	void comment(t_warnings&) {}

	/** amount of Lya removed by photoionization of excited CaII levels */
	realnum Ca2RmLya;

	/** cooling due to CaII */
	realnum Cak, 
	  Cah, 
	  Cax, 
	  Cay, 
	  Caz, 
	  Caf1, 
	  Caf2;
	double
	  Cair, 
	  c7306, 
	  Cakh;

	/** parameters dealing with dest of CaII by Lya */
	realnum dstCala, 
	  dCakh, 
	  dCaf12, 
	  Ca3d, 
	  Ca4p;

	/** summed pop of CaII excited states */
	realnum popca2ex;

	};
extern t_ca ca;


#endif /* CA_H_ */
