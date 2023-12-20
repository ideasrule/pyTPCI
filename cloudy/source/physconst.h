/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PHYSCONST_H_
#define PHYSCONST_H_

void prt_phys_constants(FILE* io);

class t_physconst : public Singleton<t_physconst>
{
	friend class Singleton<t_physconst>;

	map<string,double> p_constant;
protected:
 	t_physconst();
public:
	void add_constant(const string& s, double v)
	{
		p_constant[s] = v;
	}
	void prt_constants(FILE* io) const
	{
		map<string,double>::const_iterator p;
		for( p = p_constant.begin(); p !=  p_constant.end(); ++p )
			fprintf( io, "%-20s %.15g\n", (*p).first.c_str(), (*p).second );
	}
	double get_constant(const string& s) const
	{
		map<string,double>::const_iterator p = p_constant.find(s);
		if( p == p_constant.end() )
			return 0.;
		else
			return (*p).second;
	}
};

#ifdef HAVE_CONSTEXPR
#define NEW_CONSTANT(NAME, VALUE) constexpr double NAME = (VALUE)
#else
#define NEW_CONSTANT(NAME, VALUE) const double NAME = (VALUE)
#endif
#include "physconst_template.h"
#undef NEW_CONSTANT

#endif /* PHYSCONST_H_ */
