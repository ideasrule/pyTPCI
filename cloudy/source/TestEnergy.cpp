/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "energy.h"

namespace {

	TEST(TestEnergyCTor)
	{
		Energy Ezero;
		CHECK( fp_equal( Ezero.Ryd(), 0. ) );
		Energy Eryd0(1.);
		CHECK( fp_equal( Eryd0.Ryd(), 1. ) );
		Energy Eryd1(1.,"Ryd");
		CHECK( fp_equal( Eryd1.Ryd(), 1. ) );
		Energy Eerg(1.,"erg");
		CHECK( fp_equal( Eerg.Ryd(), 1./EN1RYD ) );
		Energy EeV(1.,"eV");
		CHECK( fp_equal( EeV.Ryd(), 1./EVRYD ) );
		Energy EkeV(1.,"keV");
		CHECK( fp_equal( EkeV.Ryd(), 1.e3/EVRYD ) );
		Energy EMeV(1.,"MeV");
		CHECK( fp_equal( EMeV.Ryd(), 1.e6/EVRYD ) );
		Energy Ewavn(1.,"cm^-1");
		CHECK( fp_equal( Ewavn.Ryd(), 1./RYD_INF ) );
		Energy Ecm(1.,"cm");
		CHECK( fp_equal( Ecm.Ryd(), RYDLAM/1.e8 ) );
		Energy Emm(1.,"mm");
		CHECK( fp_equal( Emm.Ryd(), RYDLAM/1.e7 ) );
		Energy Eum(1.,"um");
		CHECK( fp_equal( Eum.Ryd(), RYDLAM/1.e4 ) );
		Energy Enm(1.,"nm");
		CHECK( fp_equal( Enm.Ryd(), RYDLAM/1.e1 ) );
		Energy EA(1.,"A");
		CHECK( fp_equal( EA.Ryd(), RYDLAM ) );
		Energy EHz(1.,"Hz");
		CHECK( fp_equal( EHz.Ryd(), 1./FR1RYD ) );
		Energy EkHz(1.,"kHz");
		CHECK( fp_equal( EkHz.Ryd(), 1.e3/FR1RYD ) );
		Energy EMHz(1.,"MHz");
		CHECK( fp_equal( EMHz.Ryd(), 1.e6/FR1RYD ) );
		Energy EGHz(1.,"GHz");
		CHECK( fp_equal( EGHz.Ryd(), 1.e9/FR1RYD ) );
		Energy EK(1.,"K");
		CHECK( fp_equal( EK.Ryd(), 1./TE1RYD ) );
	}

	TEST(TestEnergySet)
	{
		Energy E;
		E.set(10.);
		CHECK( fp_equal( E.Ryd(), 10. ) );
		E.set(10.,"Ryd");
		CHECK( fp_equal( E.Ryd(), 10. ) );
		E.set(10.,"erg");
		CHECK( fp_equal( E.Ryd(), 10./EN1RYD ) );
		E.set(10.,"eV");
		CHECK( fp_equal( E.Ryd(), 10./EVRYD ) );
		E.set(10.,"keV");
		CHECK( fp_equal( E.Ryd(), 1.e4/EVRYD ) );
		E.set(10.,"MeV");
		CHECK( fp_equal( E.Ryd(), 1.e7/EVRYD ) );
		E.set(10.,"cm^-1");
		CHECK( fp_equal( E.Ryd(), 10./RYD_INF ) );
		E.set(10.,"cm");
		CHECK( fp_equal( E.Ryd(), RYDLAM/1.e9 ) );
		E.set(10.,"mm");
		CHECK( fp_equal( E.Ryd(), RYDLAM/1.e8 ) );
		E.set(10.,"um");
		CHECK( fp_equal( E.Ryd(), RYDLAM/1.e5 ) );
		E.set(10.,"nm");
		CHECK( fp_equal( E.Ryd(), RYDLAM/1.e2 ) );
		E.set(10.,"A");
		CHECK( fp_equal( E.Ryd(), RYDLAM/10. ) );
		E.set(10.,"Hz");
		CHECK( fp_equal( E.Ryd(), 10./FR1RYD ) );
		E.set(10.,"kHz");
		CHECK( fp_equal( E.Ryd(), 1.e4/FR1RYD ) );
		E.set(10.,"MHz");
		CHECK( fp_equal( E.Ryd(), 1.e7/FR1RYD ) );
		E.set(10.,"GHz");
		CHECK( fp_equal( E.Ryd(), 1.e10/FR1RYD ) );
		E.set(10.,"K");
		CHECK( fp_equal( E.Ryd(), 10./TE1RYD ) );
	}

	TEST(TestEnergyGet)
	{
		Energy E( 1. );
		CHECK( fp_equal( E.get("Ryd"), 1. ) );
		CHECK( fp_equal( E.get("erg"), EN1RYD ) );
		CHECK( fp_equal( E.get("eV"), EVRYD ) );
		CHECK( fp_equal( E.get("keV"), EVRYD/1.e3 ) );
		CHECK( fp_equal( E.get("MeV"), EVRYD/1.e6 ) );
		CHECK( fp_equal( E.get("cm^-1"), RYD_INF ) );
		CHECK( fp_equal( E.get("cm"), RYDLAM/1.e8 ) );
		CHECK( fp_equal( E.get("mm"), RYDLAM/1.e7 ) );
		CHECK( fp_equal( E.get("um"), RYDLAM/1.e4 ) );
		CHECK( fp_equal( E.get("nm"), RYDLAM/1.e1 ) );
		CHECK( fp_equal( E.get("A"), RYDLAM ) );
		CHECK( fp_equal( E.get("Hz"), FR1RYD ) );
		CHECK( fp_equal( E.get("kHz"), FR1RYD/1.e3 ) );
		CHECK( fp_equal( E.get("MHz"), FR1RYD/1.e6 ) );
		CHECK( fp_equal( E.get("GHz"), FR1RYD/1.e9 ) );
		CHECK( fp_equal( E.get("K"), TE1RYD ) );
	}

	TEST(TestEnergyGet2)
	{
		Energy E( 1. );
		CHECK( fp_equal( E.Ryd(), 1. ) );
		CHECK( fp_equal( E.Erg(), EN1RYD ) );
		CHECK( fp_equal( E.eV(), EVRYD ) );
		CHECK( fp_equal( E.keV(), EVRYD/1.e3 ) );
		CHECK( fp_equal( E.MeV(), EVRYD/1.e6 ) );
		CHECK( fp_equal( E.WN(), RYD_INF ) );
		CHECK( fp_equal( E.cm(), RYDLAM/1.e8 ) );
		CHECK( fp_equal( E.mm(), RYDLAM/1.e7 ) );
		CHECK( fp_equal( E.micron(), RYDLAM/1.e4 ) );
		CHECK( fp_equal( E.nm(), RYDLAM/1.e1 ) );
		CHECK( fp_equal( E.Angstrom(), RYDLAM ) );
		CHECK( fp_equal( E.Hz(), FR1RYD ) );
		CHECK( fp_equal( E.kHz(), FR1RYD/1.e3 ) );
		CHECK( fp_equal( E.MHz(), FR1RYD/1.e6 ) );
		CHECK( fp_equal( E.GHz(), FR1RYD/1.e9 ) );
		CHECK( fp_equal( E.K(), TE1RYD ) );
	}

	TEST(TestEnergyUnitConversion)
	{
		Energy E( 10. );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" RYDBERG ") ), 10. ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" ERG ")), 10.*EN1RYD ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" EV ") ), 10.*EVRYD ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" KEV ") ), EVRYD/1.e2 ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" MEV ") ), EVRYD/1.e5 ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" WAVENUMBERS ") ), 10.*RYD_INF ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" CM ") ), RYDLAM/1.e9 ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" CENTIMETER ") ), RYDLAM/1.e9 ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" MM ") ), RYDLAM/1.e8 ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" MICRON ") ), RYDLAM/1.e5 ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" NM ") ), RYDLAM/1.e2 ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" ANGSTROM ") ), RYDLAM/10. ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" HZ ") ), 10.*FR1RYD ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" KHZ ") ), FR1RYD/1.e2 ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" MHZ ") ), FR1RYD/1.e5 ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" GHZ ") ), FR1RYD/1.e8 ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" K ")  ), 10.*TE1RYD ) );
		CHECK( fp_equal( E.get( StandardEnergyUnit(" KELVIN ") ), 10.*TE1RYD ) );
	}

	TEST(TestEnergyOperator)
	{
		Energy E10( 10. );
		Energy E11( 11. );
		CHECK( E10 < E11 );
		CHECK( !( E11 < E11 ) );
		CHECK( !( E11 < E10 ) );
	}

}
