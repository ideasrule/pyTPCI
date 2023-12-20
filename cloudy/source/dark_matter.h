/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef DARK_MATTER_H_
#define DARK_MATTER_H_

#include "module.h"

struct t_dark_matter : public module
{
	const char* chName() const
	{
		return "dark";
	}
	void zero();
	void comment(t_warnings&) {}

	bool lgNFW_Set;
	double r_200;
	double r_s;

};
extern t_dark_matter dark;


#endif /* DARK_MATTER_H_ */
