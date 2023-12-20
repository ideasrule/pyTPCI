/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef HE_H_
#define HE_H_

#include "module.h"

/** he.h saves vars that are only for he */
struct t_he : public module {

	const char* chName() const
	{
		return "he";
	}

	void zero();
	void comment(t_warnings&) {}

	/** array index wihtin continuum array for Bowen 374A resonance line */
	long int ip660;

	/** these remember total He0 dest rate relative to exits from triplets */
	/** largest fraction of He0 dest due to ionization from 23S */
	double frac_he0dest_23S, 
	 /** fraction of that max due to photoioniztion */
	  frac_he0dest_23S_photo;
	/** zone where this max happened */
	long nzone;

	};
extern t_he he;


#endif /* HE_H_ */
