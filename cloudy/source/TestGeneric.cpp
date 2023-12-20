/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <fenv.h>
#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include "cddefines.h"
#include "generic_state.h"

namespace {
	TEST(extractLevelsAllBlank)
	{
		string species;
		long nLevelLo, nLevelHi;
		bool lgLevels;
		extractLevels( "", species, nLevelLo, nLevelHi, lgLevels );
		CHECK( species == "" );
		CHECK( nLevelLo == -1 && nLevelHi == -1 && ! lgLevels );
	}
	TEST(extractLevelsOnlySpecies)
	{
		string species;
		long nLevelLo, nLevelHi;
		bool lgLevels;
		extractLevels( "Cu+13", species, nLevelLo, nLevelHi, lgLevels );
		CHECK( species == "Cu+13" );
		CHECK( nLevelLo == -1 && nLevelHi == -1 && ! lgLevels );
	}
	TEST(extractLevelsOnlyLoLevel)
	{
		string species;
		long nLevelLo, nLevelHi;
		bool lgLevels;
		extractLevels( "[2:]", species, nLevelLo, nLevelHi, lgLevels );
		CHECK( species == "" );
		CHECK( nLevelLo == 2 && nLevelHi == 0 && lgLevels );
	}
	TEST(extractLevelsOnlyHiLevel)
	{
		string species;
		long nLevelLo, nLevelHi;
		bool lgLevels;
		extractLevels( "[:22]", species, nLevelLo, nLevelHi, lgLevels );
		CHECK( species == "" );
		CHECK( nLevelLo == 1 && nLevelHi == 22 && lgLevels );
	}
	TEST(extractLevelsOnlySpeciesLoLevel)
	{
		string species;
		long nLevelLo, nLevelHi;
		bool lgLevels;
		extractLevels( "H[2:]", species, nLevelLo, nLevelHi, lgLevels );
		CHECK( species == "H" );
		CHECK( nLevelLo == 2 && nLevelHi == 0 && lgLevels );
	}
	TEST(extractLevelsOnlySpeciesHiLevel)
	{
		string species;
		long nLevelLo, nLevelHi;
		bool lgLevels;
		extractLevels( "^13CO[:22]", species, nLevelLo, nLevelHi, lgLevels );
		CHECK( species == "^13CO" );
		CHECK( nLevelLo == 1 && nLevelHi == 22 && lgLevels );
	}
	TEST(extractLevelsAllGiven)
	{
		string species;
		long nLevelLo, nLevelHi;
		bool lgLevels;
		extractLevels( "CH3CH2OH[10:200]", species, nLevelLo, nLevelHi, lgLevels );
		CHECK( species == "CH3CH2OH" );
		CHECK( nLevelLo == 10 && nLevelHi == 200 && lgLevels );
	}
}
