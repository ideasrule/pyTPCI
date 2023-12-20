/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef BROKE_H_
#define BROKE_H_

/** broke.h */

#include "module.h"

class FixitList: public Singleton<FixitList>
{
	friend class Singleton<FixitList>;
protected:
	FixitList() {}
public:
	vector<string> list;
};

struct t_broke : public module {
	const char *chName() const
	{
		return "broke";
	}
	void zero() {}
	void comment(t_warnings&);
	/** logical flag saying that the code is broken, set by calling broken(); 
	 * causes a warning to be printed at the end of the calculation.  prototype
	 * is in cddefines.h */
	bool lgBroke;
	/** flag set with call to fixit, says that code needs attention, but not
	 * broken.  only causes a caution. */
	bool lgFixit;
	bool lgPrintFixits;

	/** says that new code is in place that needs to be checked */
	bool lgCheckit;

	t_broke()
	{
		// the code is ok at startup; only init here so that code remains broken
		// in a grid if any single model announces broken
		lgBroke = false;
		lgFixit = false;
		lgCheckit = false;
		lgPrintFixits = false;
	}
};
extern t_broke broke;


#endif /* BROKE_H_ */
