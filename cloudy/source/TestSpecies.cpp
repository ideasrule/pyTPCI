/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "species.h"
#include "elementnames.h"

namespace {
	TEST(TestElemSymbolToIndex)
	{
		CHECK( elem_symbol_to_index( "H" ) == ipHYDROGEN );
		CHECK( elem_symbol_to_index( "H " ) == ipHYDROGEN );

		CHECK( elem_symbol_to_index( "He" ) == ipHELIUM );
		CHECK( elem_symbol_to_index( "HE" ) == ipHELIUM );

		CHECK( elem_symbol_to_index( "Hg" ) == -1 );
	}
	TEST(TestIsElementSym)
	{
		CHECK( isElementSym( "H" ) );
		CHECK( isElementSym( "H " ) );

		CHECK( isElementSym( "He" ) );
		CHECK( isElementSym( "HE" ) );

		CHECK( ! isElementSym( "Hg" ) );
	}
	TEST(TestIsAtomicIonValid)
	{
		CHECK( ! isAtomicIonValid( ipHYDROGEN, -1 ) );
		CHECK( isAtomicIonValid( ipHYDROGEN, 0 ) );
		CHECK( isAtomicIonValid( ipHYDROGEN, 1 ) );
		CHECK( ! isAtomicIonValid( ipHYDROGEN, 2 ) );

		CHECK( isAtomicIonValid( ipHELIUM, 0 ) );
		CHECK( isAtomicIonValid( ipHELIUM, 1 ) );
		CHECK( isAtomicIonValid( ipHELIUM, 2 ) );
		CHECK( ! isAtomicIonValid( ipHELIUM, 3 ) );

		CHECK( isAtomicIonValid( ipZINC, 0 ) );
		CHECK( isAtomicIonValid( ipZINC, 1 ) );
		CHECK( isAtomicIonValid( ipZINC, 30 ) );
		CHECK( ! isAtomicIonValid( ipZINC, 31 ) );
	}
	TEST(TestBareNucleus)
	{
		CHECK( ! isBareNucleus( ipHYDROGEN, -1 ) );
		CHECK( ! isBareNucleus( ipHYDROGEN, 0 ) );
		CHECK( isBareNucleus( ipHYDROGEN, ipHYDROGEN+1 ) );
		CHECK( ! isBareNucleus( ipHYDROGEN, ipHYDROGEN+2 ) );

		CHECK( ! isBareNucleus( ipZINC, 0 ) );
		CHECK( isBareNucleus( ipZINC, ipZINC+1 ) );
		CHECK( ! isBareNucleus( ipZINC, ipZINC+2 ) );
	}
	TEST(TestParseChemical)
	{
		string species;
		string element_or_molecule;
		long ion_stage;

		species = "";
		CHECK( ! parse_chemical( species, element_or_molecule, ion_stage ) );
		CHECK( element_or_molecule == "" );
		CHECK( ion_stage == -1 );

		species = "H";
		CHECK( parse_chemical( species, element_or_molecule, ion_stage ) );
		CHECK( element_or_molecule == "H" );
		CHECK( ion_stage == 0 );

		species = "H+";
		CHECK( parse_chemical( species, element_or_molecule, ion_stage ) );
		CHECK( element_or_molecule == "H" );
		CHECK( ion_stage == 1 );

		species = "H+1";
		CHECK( parse_chemical( species, element_or_molecule, ion_stage ) );
		CHECK( element_or_molecule == "H" );
		CHECK( ion_stage == 1 );

		species = "H+0";
		CHECK( parse_chemical( species, element_or_molecule, ion_stage ) );
		CHECK( element_or_molecule == "H" );
		CHECK( ion_stage == 0 );

		species = "HCl";
		CHECK( parse_chemical( species, element_or_molecule, ion_stage ) );
		CHECK( element_or_molecule == "HCl" );
		CHECK( ion_stage == 0 );

		species = "He+";
		CHECK( parse_chemical( species, element_or_molecule, ion_stage ) );
		CHECK( element_or_molecule == "He" );
		CHECK( ion_stage == 1 );

		species = "He+2";
		CHECK( parse_chemical( species, element_or_molecule, ion_stage ) );
		CHECK( element_or_molecule == "He" );
		CHECK( ion_stage == 2 );

		species = "Fe+24";
		CHECK( parse_chemical( species, element_or_molecule, ion_stage ) );
		CHECK( element_or_molecule == "Fe" );
		CHECK( ion_stage == 24 );

		species = "HCl+";
		CHECK( parse_chemical( species, element_or_molecule, ion_stage ) );
		CHECK( element_or_molecule == "HCl" );
		CHECK( ion_stage == 1 );
	}
	TEST(ChemicalToSpectral)
	{
		string spectrum;

		chemical_to_spectral( "", spectrum );
		CHECK( spectrum == "" );

		chemical_to_spectral( "H", spectrum );
		CHECK( spectrum == "H  1" );

		chemical_to_spectral( "H+0", spectrum );
		CHECK( spectrum == "H  1" );

		chemical_to_spectral( "HCl", spectrum );
		CHECK( spectrum == "HCl" );

		chemical_to_spectral( "He+", spectrum );
		CHECK( spectrum == "He 2" );

		chemical_to_spectral( "He+2", spectrum );
		CHECK( spectrum == "" );

		chemical_to_spectral( "Fe+24", spectrum );
		CHECK( spectrum == "Fe25" );

		chemical_to_spectral( "HCl+", spectrum );
		CHECK( spectrum == "HCl+" );
	}
}
