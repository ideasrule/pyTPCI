/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PREDCONT_H_
#define PREDCONT_H_

/* predcont.h */

#include "energy.h"

/** array of energies and pointers for continuum points */
class t_PredCont : public Singleton<t_PredCont>
{
	friend class Singleton<t_PredCont>;
protected:
	t_PredCont();
private:
	vector<EnergyEntry> p_val;
	// the offset of the nFnu entries in the line stack.
	long p_offset;
public:
	long find(double energy, const char* unit = "Ryd") const;
	long add(double energy, const char* unit = "Ryd");
	EnergyEntry& operator[] ( size_t i )
	{
		return p_val[i];
	}
	size_t size() const
	{
		return p_val.size();
	}
	void set_offset(long offset)
	{
		p_offset = offset;
	}
	long offset() const
	{
		return p_offset;
	}
};

#endif /* PREDCONT_H_ */
