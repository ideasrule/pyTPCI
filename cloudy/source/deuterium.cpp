/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "deuterium.h"
#include "mole.h"

t_deuterium deut;

void t_deuterium::zero()
{
	xIonDense[0] = 0.;
	xIonDense[1] = 0.;
	m_xMolecules = 0.;
}

void t_deuterium::updateXMolecules()
{
	total_molecule_deut(m_xMolecules);
}

void ScaleDensitiesDeuterium( const realnum &factor )
{
	deut.gas_phase *= factor;
	deut.xIonDense[0] *= factor;
	deut.xIonDense[1] *= factor;
	return;
}

void InitDeuteriumIonization()
{
	deut.zero();
	deut.xIonDense[0] = deut.gas_phase;
	deut.xIonDense[1] = 0.;
	return;
}

void SetDeuteriumIonization( const double &xNeutral, const double &xIonized )
{
	if( !deut.lgElmtOn )
		return;

	realnum total = deut.gas_phase - deut.xMolecules();

	double neut = total * xNeutral/( xNeutral + xIonized );
	double ionz = total * xIonized/( xNeutral + xIonized );
	
	deut.xIonDense[0] = neut;
	deut.xIonDense[1] = ionz;

	//if( iteration >= 3 )
	//	fprintf( ioQQQ, "DEBUGGG SetDeuteriumIonization %li\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
	//	nzone, deut.gas_phase, deut.xMolecules(), total, deut.xIonDense[0], deut.xIonDense[1] );
	
	return;
}

void SetDeuteriumFractionation( const realnum &frac )
{
	if( !deut.lgElmtOn )
		return;
	deut.fractionation = frac;
	return;
}

void SetGasPhaseDeuterium( const realnum &Hdensity )
{
	if( !deut.lgElmtOn )
		return;
	deut.gas_phase = Hdensity * deut.fractionation;
}

