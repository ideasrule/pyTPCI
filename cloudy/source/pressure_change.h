/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PRESSURE_CHANGE_H_
#define PRESSURE_CHANGE_H_

/**PressureChange evaluate the current pressure, and change needed to get it to PresTotlInit,
 * return value is true is density was changed, false if no changes were necessary
 * \param dP_chng_factor this is change factor, 1 at first, becomes smaller as oscillations occur
 */ 
class PresMode
{
public:
	int global, zone;
	void set();
};

class solverState
{
public:
	double dp, erp, press;
	int lastzone;
	explicit solverState() : dp(-1.), erp(-1.), press(0), lastzone(-1) {}
};

double pressureZone(const PresMode &presmode);
void PressureChange(double dP_chng_factor, const PresMode &presmode, solverState &st, bool& lgStable);
double zoneDensity();


#endif // PRESSURE_CHANGE_H_
