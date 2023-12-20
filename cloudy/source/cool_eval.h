/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef COOL_EVAL_H_
#define COOL_EVAL_H_

void eeBremsSpectrum(double Te,            // electron temperature in K
					 vector<realnum>& jnu, // 4pi*jnu in photons/cm^3/s/cell
					 vector<double>& knu); // opacity in cm^-1

#endif // COOL_EVAL_H_
