/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <fenv.h>
#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include "cddefines.h"
#include "prt.h"

namespace {
	TEST(PrtMatrixPruneSpeciesNoLevels)
	{
		prt.matrix.setSpecies( "Cu+13" );
		CHECK( prt.matrix.species == "Cu+13" );
	}
	TEST(PrtMatrixPruneSpeciesAllLevels)
	{
		prt.matrix.setSpecies( "Cu+13[:]" );
		CHECK( prt.matrix.species == "Cu+13" );
	}
	TEST(PrtMatrixPruneSpeciesFromLevels)
	{
		prt.matrix.setSpecies( "Cu+13[1:10]" );
		CHECK( prt.matrix.species == "Cu+13" );
	}

	TEST(PrtMatrixPruneSpeciesLevelsNoLevels)
	{
		prt.matrix.setSpecies( "Cu+13" );
		CHECK( prt.matrix.speciesLevels == "Cu+13[:]" );
	}
	TEST(PrtMatrixPruneSpeciesLevelsAllLevels)
	{
		prt.matrix.setSpecies( "Cu+13[:]" );
		CHECK( prt.matrix.speciesLevels == "Cu+13[:]" );
	}
	TEST(PrtMatrixPruneSpeciesLevelsFromSpecies)
	{
		prt.matrix.setSpecies( "Cu+13[1:10]" );
		CHECK( prt.matrix.speciesLevels == "Cu+13[1:10]" );
	}

	//
	// prt.matrix.resolveLevels() cannot be tested,
	// unless we read in atomic data
	//
}
