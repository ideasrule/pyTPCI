/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "collision.h"

#include "dense.h"
#include "h2.h"
#include "hmi.h"

ColliderList::ColliderList(const t_dense &d)
{
	DEBUG_ENTRY( "ColliderList::ColliderList()" );

	realnum e_amu = ELECTRON_MASS*AVOGADRO;
	// weight of 1 eV of binding energy in amu
	realnum eV_amu = EN1EV*AVOGADRO/pow2(SPEEDLIGHT);

	list.resize( ipNCOLLIDER );
	list[ipELECTRON].charge = -1;
	list[ipELECTRON].mass_amu = e_amu;

	// the masses derived below are for a standard mixture of isotopes, not for p or alpha!
	list[ipPROTON].charge = 1;
	list[ipPROTON].mass_amu = d.AtomicWeight[ipHYDROGEN] - e_amu + 13.6f*eV_amu;

	list[ipHE_PLUS].charge = 1;
	list[ipHE_PLUS].mass_amu = d.AtomicWeight[ipHELIUM] - e_amu + 24.6f*eV_amu;

	list[ipALPHA].charge = 2;
	list[ipALPHA].mass_amu = d.AtomicWeight[ipHELIUM] - 2.f*e_amu + 79.0f*eV_amu;

	list[ipATOM_H].charge = 0;
	list[ipATOM_H].mass_amu = d.AtomicWeight[ipHYDROGEN];

	list[ipATOM_HE].charge = 0;
	list[ipATOM_HE].mass_amu = d.AtomicWeight[ipHELIUM];

	list[ipH2_ORTHO].charge = 0;
	list[ipH2_ORTHO].mass_amu = 2.f*d.AtomicWeight[ipHYDROGEN] - 4.75f*eV_amu;

	list[ipH2_PARA].charge = 0;
	list[ipH2_PARA].mass_amu = 2.f*d.AtomicWeight[ipHYDROGEN] - 4.75f*eV_amu;

	list[ipH2].charge = 0;
	list[ipH2].mass_amu = 2.f*d.AtomicWeight[ipHYDROGEN] - 4.75f*eV_amu;
}

void ColliderList::init()
{
	DEBUG_ENTRY( "ColliderList::init()" );

	list[ipELECTRON].density = &(dense.EdenHCorr);
	list[ipPROTON].density = &(dense.xIonDense[ipHYDROGEN][1]);
	list[ipHE_PLUS].density = &(dense.xIonDense[ipHELIUM][1]);
	list[ipALPHA].density = &(dense.xIonDense[ipHELIUM][2]);
	list[ipATOM_H].density = &(dense.xIonDense[ipHYDROGEN][0]);
	list[ipATOM_HE].density = &(dense.xIonDense[ipHELIUM][0]);
	list[ipH2_ORTHO].density = &(h2.ortho_density);
	list[ipH2_PARA].density = &(h2.para_density);
	list[ipH2].density = &(hmi.H2_total);
}

/*CollisionJunk set all elements of transition struc to dangerous values */
void CollisionJunk( const CollisionProxy & t )
{

	DEBUG_ENTRY( "CollisionJunk()" );

	/* Coll->cooling and Coll->heating due to collisional excitation */
	t.cool() = -FLT_MAX;
	t.heat() = -FLT_MAX;

	/* collision strengths for transition */
	t.col_str() = -FLT_MAX;

	t.is_gbar() = -1;

	for( long i=0; i<ipNCOLLIDER; i++ )
		t.rate_coef_ul_set()[i] = 0.f;

	t.rate_lu_nontherm_set() = 0.f;

	return;
}

/*CollisionZero zeros out the structure */
void CollisionZero( const CollisionProxy & t )
{

	DEBUG_ENTRY( "CollisionZero()" );

	/* Coll->cooling and Coll->heating due to collisional excitation */
	t.cool() = 0.;
	t.heat() = 0.;

	return;
}

