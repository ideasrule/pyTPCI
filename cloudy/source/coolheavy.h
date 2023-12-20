/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef COOLHEAVY_H_
#define COOLHEAVY_H_

#include "module.h"

/** coolheavy.h */
struct t_CoolHeavy : public module {
	const char *chName() const
	{
		return "CoolHeavy";
	}
	void zero();
	void comment(t_warnings&) {}

	double Fe231,
	  Fe232, 
	  Fe221,
	  eebrm, 
	  h2line, 
	  HD,
	  colmet, 
	  tccool, 
	  expans, 
	  cextxx,
	  cyntrn, 
	  heavfb;

	  /* Fe12 lines added 01 Aug
	  double fe13_1216 , fe13_3000 , fe13_2302; */

	  /** flag set false if no free-free heating command entered */
	  bool lgFreeOn;
	  double brems_cool_h, 
		 brems_cool_hminus,
		 brems_cool_he, 
		 brems_cool_metals, 
		 /** total brems heating, all opacity sources */
		 brems_heat_total,
		 /** net brems cooling, sum of cool minus heat */
		 brems_cool_net;


	/** cooling due to H + H+ => H2+ */
	realnum H2PlsCool;

	};
extern t_CoolHeavy CoolHeavy;

#endif /* COOLHEAVY_H_ */
