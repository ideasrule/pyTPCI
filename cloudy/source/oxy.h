/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef OXY_H_
#define OXY_H_

/** oxy.h */

#include "module.h"

struct t_oxy : public module {
	const char *chName() const
	{
		return "oxy";
	}
	void zero();
	void comment(t_warnings&) {}

	realnum poiii2, 
	  d5007r, 
	  d6300,
	  poiii3, 
	  poiii2Max, 
	  poiii3Max, 
	  r5007Max, 
	  r4363Max, 
	  poiexc, 
	  poimax;

	/** photoionization rate for OII, producing 1665 */
	realnum AugerO3;

	/** lines produced by relaxation following photoionization */
	realnum s3727, 
	  s7325;

	/** array indices for photoionization into excited state of oii */
	long int i2d, 
	  i2p;

	};
extern t_oxy oxy;


#endif /* OXY_H_ */
