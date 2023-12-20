/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef FLUX_H_
#define FLUX_H_

#include "energy.h"

class Flux
{
	typedef enum {
		FU_NONE, FU_ERG_S, FU_W, FU_JY, FU_MJY, FU_MJY_SR, FU_CM2,
		FU_M2, FU_A, FU_NM, FU_MU, FU_HZ, FU_SR, FU_SQAS, FU_TOP
	} fu_flag;
	typedef bitset<FU_TOP> fu_bits;

	Energy p_energy;
	double p_flux;
	fu_bits p_userunits;

	fu_bits p_InternalFluxUnitNoCheck(const string& unit, size_t& len) const;
	fu_bits p_InternalFluxUnit(const string& unit) const;
	bool p_ValidFluxUnit(fu_bits) const;
	void p_set(Energy e, double value, fu_bits bits);
	double p_get(fu_bits bits) const;
public:
	// CONSTRUCTORS
	Flux()
	{
		set(0., 0.);
		p_userunits.reset();
	}
	Flux(Energy e, double flux)
	{
		set(e, flux);
	}
	Flux(Energy e, double flux, const string& unit)
	{
		set(e, flux, unit);
	}
	// MUTATORS
	void set(Energy e, double flux)
	{
		set( e, flux, "erg/s/cm2" );
	}
	void set(Energy e, double flux, const string& unit)
	{
		p_set( e, flux, p_InternalFluxUnit(unit) );
	}
	// ACCESSORS
	double get() const
	{
		return p_flux;
	}
	double get(const string& unit) const
	{
		return p_get( p_InternalFluxUnit(unit) );
	}
	Energy E() const
	{
		return p_energy;
	}
	string uu() const;
	friend inline bool ValidFluxUnit(const string& unit);
};

// parse standard flux unit from Cloudy input line
string StandardFluxUnit(const char*);

// check if standard unit is valid
inline bool ValidFluxUnit(const string& unit)
{
	Flux f;
	size_t p;
	Flux::fu_bits bits = f.p_InternalFluxUnitNoCheck(unit,p);
	return ( p == unit.length() && f.p_ValidFluxUnit(bits) );
}

#endif /* FLUX_H_ */
