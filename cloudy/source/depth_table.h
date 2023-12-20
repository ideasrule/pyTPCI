/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef DEPTH_TABLE_H_
#define DEPTH_TABLE_H_

class DepthTable
{
public:
        /*true if quantities are interpreted linearly*/
        bool lgLinear;
	/* lg is true if depth, false if radius to be used*/
	bool lgDepth;
	/**dist is log radius in cm, val is log value*/
	vector<double> dist;
	vector<double> val;

	/**number of values in above table */
	long int nvals;

	/**tabval, adapted from dense_tabden interpolate on table of points for density with dlaw table command, by K Volk 
		\param r0
		\param depth
	*/
	double tabval( double r0, double depth) const;
	void clear()
	{
		nvals = 0;
		dist.resize(0);
		val.resize(0);
	}
};

#endif // DEPTH_TABLE_H_
