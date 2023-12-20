/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef DEUTERIUM_H_
#define DEUTERIUM_H_

#include "module.h"

class t_deuterium : public module
{
public:
	t_deuterium()
	{
		lgElmtOn = false;
		gas_phase = 0.f;
		xIonDense[0] = 0.;
		xIonDense[1] = 0.;
		fractionation = 0.f;
		m_xMolecules = 0.f;
	}
	bool lgElmtOn;
	realnum gas_phase;
	double xIonDense[2];
	realnum fractionation;
private:
	/** xMolecules density locked in molecules, 
	 * this is included in gas_phase */
	realnum m_xMolecules;
public:
	realnum xMolecules(void)
	{
		return m_xMolecules;
	}
	void updateXMolecules();
	void zero();
	void comment(t_warnings&) {}
	const char* chName() const
	{
		return "deuterium";
	}
};

extern t_deuterium deut;

void ScaleDensitiesDeuterium( const realnum &factor );
void SetDeuteriumFractionation( const realnum &frac );
void SetGasPhaseDeuterium( const realnum &Hdensity );
void SetDeuteriumIonization( const double &xNeutral, const double &xIonized );
void InitDeuteriumIonization();

#endif /* DEUTERIUM_H_ */
