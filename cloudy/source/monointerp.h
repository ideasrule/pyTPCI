/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MONOINTERP_H_
#define MONOINTERP_H_

// Usage:
// Monointerp m(xvals,yvals, npt); -- Constructor
// y_interp = m(x_interp);         -- Interpolate value

class Monointerp {
	const std::vector<double> m_x, m_y;
	std::vector<double> m_g;

public:
	// CONSTRUCTORS
	Monointerp ( const double x[], const double y[], long n);
	~Monointerp(void);
private:
	Monointerp ( const Monointerp& );             // Not implemented
	Monointerp& operator= ( const Monointerp& );  // Not implemented

	// MANIPULATORS -- None

	// ACCESSORS
public:
	double operator() ( double xval) const;

};
#endif
