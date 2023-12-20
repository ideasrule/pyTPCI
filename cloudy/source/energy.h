/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ENERGY_H_
#define ENERGY_H_

#include "physconst.h"

class Energy
{
	double m_energy;
public:
	// CONSTRUCTORS
	Energy () : m_energy(0.0) {}
	Energy (double energy) : m_energy(energy) {}
	Energy (double energy, const char *unit) : m_energy(0.0)
	{
		set(energy, unit);
	}
	// OPERATORS
	bool operator< (const Energy& e) const
	{
		return ( m_energy < e.m_energy );
	}
	// MUTATORS
	void set(double energy)
	{
		m_energy = energy;
	}
	void set(double energy, const char *unit);
	// ACCESSORS
	double get(const char * unit) const;
	double Ryd() const
	{
		return m_energy;
	}
	double K() const
	{
		return m_energy*TE1RYD;
	}
	double Erg() const
	{
		return m_energy*EN1RYD;
	}
	double WN() const
	{
		return m_energy*RYD_INF;
	}
	double eV() const
	{
		return m_energy*EVRYD;
	}
	double keV() const
	{
		return 1e-3*eV();
	}
	double MeV() const
	{
		return 1e-6*eV();
	}
	double Hz() const
	{
		return m_energy*FR1RYD;
	}
	double kHz() const
	{
		return Hz()*1e-3;
	}
	double MHz() const
	{
		return Hz()*1e-6;
	}
	double GHz() const
	{
		return Hz()*1e-9;
	}
	double Angstrom() const
	{
		return RYDLAM/m_energy;
	}
	double nm() const
	{
		return Angstrom()*1e-1;
	}
	double micron() const
	{
		return Angstrom()*1e-4;
	}
	double mm() const
	{
		return Angstrom()*1e-7;
	}
	double cm() const
	{
		return Angstrom()*1e-8;
	}
};

//
//! EnergyEntry: class for storing a continuum energy and its associated pointer as a pair.
//! This class is safe to construct even before the mesh is set up, as in that case calculating
//! the pointer is delayed until it is actually needed. The energy can be changed after
//! construction using the set() methods, but only if the mesh is already set up.
//
class EnergyEntry : public Energy
{
	long p_ip;
	void p_set_ip();
public:
	EnergyEntry() : Energy(0.)
	{
		p_ip = -1;
	}
	EnergyEntry(double energy) : Energy(energy)
	{
		p_ip = -1;
	}
	EnergyEntry(double energy, const char *unit) : Energy(energy,unit)
	{
		p_ip = -1;
	}
	void set(double energy, const char *unit)
	{
		Energy::set(energy,unit);
		p_set_ip();
	}
	void set(double energy)
	{
		Energy::set(energy);
		p_set_ip();
	}
	// pointer on C scale
	long ip_C()
	{
		// doing it this way assures that we can create an EnergyEntry before
		// the mesh is set up; ipoint() will check that mesh is actually set up
		if( p_ip < 0 )
			p_set_ip();
		return p_ip;
	}
	// pointer on fortran scale
	long ip_fortran()
	{
		return ip_C() + 1;
	}
};

const char *StandardEnergyUnitNoAbort(const char *chCard);

inline const char *StandardEnergyUnit(const char *chCard)
{
	const char *unit = StandardEnergyUnitNoAbort(chCard);
	if( unit == nullptr )
	{
		fprintf( ioQQQ, " No energy / wavelength unit was recognized on this line:\n %s\n\n", chCard );
		fprintf( ioQQQ, " See Hazy for details.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return unit;
}

void ConserveEnergy();

#endif /* ENERGY_H_ */

