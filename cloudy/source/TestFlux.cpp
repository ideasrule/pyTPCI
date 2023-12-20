/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "flux.h"
#include "parser.h"

namespace {

	TEST(TestFluxCTor)
	{
		Energy E(1.,"Ryd");
		Flux Fzero;
		CHECK( fp_equal( Fzero.get(), 0. ) );
		Flux Ferg_s_cm20(E,1.);
		CHECK( fp_equal( Ferg_s_cm20.get(), 1. ) );
		Flux Ferg_s_cm21(E,1.,"erg/s/cm2");
		CHECK( fp_equal( Ferg_s_cm21.get(), 1. ) );
		Flux Ferg_s_cm2_A(E,1.,"erg/s/cm2/A");
		CHECK( fp_equal( Ferg_s_cm2_A.get(), E.Angstrom() ) );
		Flux Ferg_s_cm2_micron(E,1.,"erg/s/cm2/micron");
		CHECK( fp_equal( Ferg_s_cm2_micron.get(), E.micron() ) );
		Flux Ferg_s_cm2_Hz(E,1.,"erg/s/cm2/Hz");
		CHECK( fp_equal( Ferg_s_cm2_Hz.get(), E.Hz() ) );
		Flux FW_m2(E,1.,"W/m2");
		CHECK( fp_equal( FW_m2.get(), 1.e3 ) );
		Flux FW_m2_A(E,1.,"W/m2/A");
		CHECK( fp_equal( FW_m2_A.get(), 1.e3*E.Angstrom() ) );
		Flux FW_m2_nm(E,1.,"W/m2/nm");
		CHECK( fp_equal( FW_m2_nm.get(), 1.e3*E.nm() ) );
		Flux FW_m2_micron(E,1.,"W/m2/micron");
		CHECK( fp_equal( FW_m2_micron.get(), 1.e3*E.micron() ) );
		Flux FW_m2_Hz(E,1.,"W/m2/Hz");
		CHECK( fp_equal( FW_m2_Hz.get(), 1.e3*E.Hz() ) );
		Flux FW_cm2(E,1.,"W/cm2");
		CHECK( fp_equal( FW_cm2.get(), 1.e7 ) );
		Flux FW_cm2_A(E,1.,"W/cm2/A");
		CHECK( fp_equal( FW_cm2_A.get(), 1.e7*E.Angstrom() ) );
		Flux FW_cm2_micron(E,1.,"W/cm2/micron");
		CHECK( fp_equal( FW_cm2_micron.get(), 1.e7*E.micron() ) );
		Flux FW_cm2_Hz(E,1.,"W/cm2/Hz");
		CHECK( fp_equal( FW_cm2_Hz.get(), 1.e7*E.Hz() ) );
		Flux FW_cm2_Hz_sr(E,1.,"W/cm2/Hz/sr");
		CHECK( fp_equal( FW_cm2_Hz_sr.get(), 1.e7*E.Hz()*PI4 ) );
		Flux FW_cm2_Hz_arcsec2(E,1.,"W/cm2/Hz/arcsec2");
		CHECK( fp_equal( FW_cm2_Hz_arcsec2.get(), 1.e7*E.Hz()*SQAS_SKY ) );
		Flux FJy(E,1.,"Jy");
		CHECK( fp_equal( FJy.get(), 1.e-23*E.Hz() ) );
		Flux FmJy(E,1.,"mJy");
		CHECK( fp_equal( FmJy.get(), 1.e-26*E.Hz() ) );
		Flux FMJy_sr(E,1.,"MJy/sr");
		CHECK( fp_equal( FMJy_sr.get(), 1.e-17*E.Hz()*PI4 ) );
	}

	TEST(TestFluxSet)
	{
		Energy E(1.,"Ryd");
		Flux F;
		F.set(E,10.);
		CHECK( fp_equal( F.get(), 10. ) );
		CHECK_EQUAL( "erg/s/cm2", F.uu() );

		F.set(E,10.,"erg/s/cm2");
		CHECK( fp_equal( F.get(), 10. ) );
		CHECK_EQUAL( "erg/s/cm2", F.uu() );

		F.set(E,10.,"erg/s/cm2/A");
		CHECK( fp_equal( F.get(), 10.*E.Angstrom() ) );
		CHECK_EQUAL( "erg/s/cm2/A", F.uu() );

		F.set(E,10.,"erg/s/cm2/micron");
		CHECK( fp_equal( F.get(), 10.*E.micron() ) );
		CHECK_EQUAL( "erg/s/cm2/micron", F.uu() );

		F.set(E,10.,"erg/s/cm2/Hz");
		CHECK( fp_equal( F.get(), 10.*E.Hz() ) );
		CHECK_EQUAL( "erg/s/cm2/Hz", F.uu() );

		F.set(E,10.,"W/m2");
		CHECK( fp_equal( F.get(), 1.e4 ) );
		CHECK_EQUAL( "W/m2", F.uu() );

		F.set(E,10.,"W/m2/A");
		CHECK( fp_equal( F.get(), 1.e4*E.Angstrom() ) );
		CHECK_EQUAL( "W/m2/A", F.uu() );

		F.set(E,10.,"W/m2/nm");
		CHECK( fp_equal( F.get(), 1.e4*E.nm() ) );
		CHECK_EQUAL( "W/m2/nm", F.uu() );

		F.set(E,10.,"W/m2/micron");
		CHECK( fp_equal( F.get(), 1.e4*E.micron() ) );
		CHECK_EQUAL( "W/m2/micron", F.uu() );

		F.set(E,10.,"W/m2/Hz");
		CHECK( fp_equal( F.get(), 1.e4*E.Hz() ) );
		CHECK_EQUAL( "W/m2/Hz", F.uu() );

		F.set(E,10.,"W/cm2");
		CHECK( fp_equal( F.get(), 1.e8 ) );
		CHECK_EQUAL( "W/cm2", F.uu() );

		F.set(E,10.,"W/cm2/A");
		CHECK( fp_equal( F.get(), 1.e8*E.Angstrom() ) );
		CHECK_EQUAL( "W/cm2/A", F.uu() );

		F.set(E,10.,"W/cm2/micron");
		CHECK( fp_equal( F.get(), 1.e8*E.micron() ) );
		CHECK_EQUAL( "W/cm2/micron", F.uu() );

		F.set(E,10.,"W/cm2/Hz");
		CHECK( fp_equal( F.get(), 1.e8*E.Hz() ) );
		CHECK_EQUAL( "W/cm2/Hz", F.uu() );

		F.set(E,10.,"W/cm2/Hz/sr");
		CHECK( fp_equal( F.get(), 1.e8*E.Hz()*PI4 ) );
		CHECK_EQUAL( "W/cm2/Hz/sr", F.uu() );

		F.set(E,10.,"W/cm2/Hz/arcsec2");
		CHECK( fp_equal( F.get(), 1.e8*E.Hz()*SQAS_SKY ) );
		CHECK_EQUAL( "W/cm2/Hz/arcsec2", F.uu() );

		F.set(E,10.,"Jy");
		CHECK( fp_equal( F.get(), 1.e-22*E.Hz() ) );
		CHECK_EQUAL( "Jy", F.uu() );

		F.set(E,10.,"mJy");
		CHECK( fp_equal( F.get(), 1.e-25*E.Hz() ) );
		CHECK_EQUAL( "mJy", F.uu() );

		F.set(E,10.,"MJy/sr");
		CHECK( fp_equal( F.get(), 1.e-16*E.Hz()*PI4 ) );
		CHECK_EQUAL( "MJy/sr", F.uu() );
	}

	TEST(TestFluxGet)
	{
		Energy E(1.,"Ryd");
		Flux F(E,1.);
		CHECK( fp_equal( F.get("erg/s/cm2"), 1. ) );
		CHECK( fp_equal( F.get("erg/s/cm2/A"), 1./E.Angstrom() ) );
		CHECK( fp_equal( F.get("erg/s/cm2/micron"), 1./E.micron() ) );
		CHECK( fp_equal( F.get("erg/s/cm2/Hz"), 1./E.Hz() ) );
		CHECK( fp_equal( F.get("W/m2"), 1.e-3 ) );
		CHECK( fp_equal( F.get("W/m2/A"), 1.e-3/E.Angstrom() ) );
		CHECK( fp_equal( F.get("W/m2/micron"), 1.e-3/E.micron() ) );
		CHECK( fp_equal( F.get("W/m2/Hz"), 1.e-3/E.Hz() ) );
		CHECK( fp_equal( F.get("W/cm2"), 1.e-7 ) );
		CHECK( fp_equal( F.get("W/cm2/A"), 1.e-7/E.Angstrom() ) );
		CHECK( fp_equal( F.get("W/cm2/nm"), 1.e-7/E.nm() ) );
		CHECK( fp_equal( F.get("W/cm2/micron"), 1.e-7/E.micron() ) );
		CHECK( fp_equal( F.get("W/cm2/Hz"), 1.e-7/E.Hz() ) );
		CHECK( fp_equal( F.get("W/cm2/Hz/sr"), 1.e-7/E.Hz()/PI4 ) );
		CHECK( fp_equal( F.get("W/cm2/Hz/arcsec2"), 1.e-7/E.Hz()/SQAS_SKY ) );
		CHECK( fp_equal( F.get("Jy"), 1.e23/E.Hz() ) );
		CHECK( fp_equal( F.get("mJy"), 1.e26/E.Hz() ) );
		CHECK( fp_equal( F.get("MJy/sr"), 1.e17/E.Hz()/PI4 ) );
	}

	TEST(TestFluxUnitConversion)
	{
		Energy E(1.,"Ryd");
		Flux F(E,10.);
		CHECK( fp_equal( F.get( StandardFluxUnit("ERG/S/SQCM") ), 10. ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("ERG/S/SQCM/A ") ), 10./F.E().Angstrom() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("ERG/S/SQCM/MICRON") ), 10./F.E().micron() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("ERG/S/SQCM/HZ") ), 10./F.E().Hz() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQM") ), 1.e-2 ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQM/A ") ), 1.e-2/F.E().Angstrom() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQM/MICRON") ), 1.e-2/F.E().micron() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQM/HZ") ), 1.e-2/F.E().Hz() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQCM") ), 1.e-6 ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQCM/A ") ), 1.e-6/F.E().Angstrom() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQCM/MICRON") ), 1.e-6/F.E().micron() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQCM/HZ") ), 1.e-6/F.E().Hz() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQCM/NM/SR") ), 1.e-6/F.E().nm()/PI4 ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQCM/HZ/SR") ), 1.e-6/F.E().Hz()/PI4 ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("W/SQCM/HZ/SQAS") ), 1.e-6/F.E().Hz()/SQAS_SKY ) );
		CHECK( fp_equal( F.get( StandardFluxUnit(" JY ") ), 1.e24/F.E().Hz() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("JANSKY") ), 1.e24/F.E().Hz() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit(" MJY") ), 1.e27/F.E().Hz() ) );
		CHECK( fp_equal( F.get( StandardFluxUnit("MJY/SR") ), 1.e18/F.E().Hz()/PI4 ) );
	}

	TEST(TestValidFluxUnit)
	{
		CHECK( ValidFluxUnit( "erg/s/cm2/arcsec2" ) );
		CHECK( !ValidFluxUnit( "erg/s/cm2/arcsec" ) );
		CHECK( !ValidFluxUnit( "erg/cm2/arcsec2" ) );
		CHECK( !ValidFluxUnit( "erg/s/cm2/cm2/arcsec2" ) );
		CHECK( !ValidFluxUnit( "erg/s/cm2/A/nm" ) );
		CHECK( !ValidFluxUnit( "erg/s/cm2/" ) );
		CHECK( !ValidFluxUnit( "MJy/arcsec2" ) );
		CHECK( !ValidFluxUnit( "mJy/sr" ) );
		CHECK( !ValidFluxUnit( "Jy/sr" ) );
	}

	// test some command lines that could possibly confuse the unit readers...
	TEST(TestEnergyFluxParser1)
	{
		char chCard[100];
		strcpy( chCard, "OPTIMIZE FLUX 0.03 RYD 24 ERG/S/SQCM/MICRON" );
		Parser p;
		p.setline(chCard);
		double energy = p.FFmtRead();
		double flux = p.FFmtRead();
		Energy E( energy, p.StandardEnergyUnit() );
		Flux F( E, flux, p.StandardFluxUnit() );
		CHECK( fp_equal( E.Ryd(), 0.03 ) );
		CHECK( fp_equal( F.get("erg/s/cm2/micron"), 24. ) );
	}

	TEST(TestEnergyFluxParser2)
	{
		char chCard[100];
		strcpy( chCard, "OPTIMIZE FLUX 12 MICRON 24 ERG/S/SQCM" );
		Parser p;
		p.setline( chCard );
		double energy = p.FFmtRead();
		double flux = p.FFmtRead();
		Energy E( energy, p.StandardEnergyUnit() );
		Flux F( E, flux, p.StandardFluxUnit() );
		CHECK( fp_equal( E.micron(), 12. ) );
		CHECK( fp_equal( F.get("erg/s/cm2"), 24. ) );
	}

	TEST(TestEnergyFluxParser3)
	{
		char chCard[100];
		strcpy( chCard, "OPTIMIZE FLUX 12E-14 ERG 24 W/SQCM" );
		Parser p;
		p.setline(chCard );
		double energy = p.FFmtRead();
		double flux = p.FFmtRead();
		Energy E( energy, p.StandardEnergyUnit() );
		Flux F( E, flux, p.StandardFluxUnit() );
		CHECK( fp_equal( E.Erg(), 12.e-14 ) );
		CHECK( fp_equal( F.get("W/cm2"), 24. ) );
	}

}
