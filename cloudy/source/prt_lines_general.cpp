/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_general put general information and energetics into line intensity stack */
#include "cddefines.h"
#include "taulines.h"
#include "coolheavy.h"
#include "thermal.h"
#include "continuum.h"
#include "geometry.h"
#include "dynamics.h"
#include "iso.h"
#include "rfield.h"
#include "trace.h"
#include "ionbal.h"
#include "radius.h"
#include "lines.h"
#include "dense.h"
#include "lines_service.h"

void lines_general(void)
{
	long int i, 
	  nelem, 
	  ipnt;

	double 
	  HeatMetal ,
	  ee511; 

	DEBUG_ENTRY( "lines_general()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_general called\n" );
	}

	i = StuffComment( "general properties" );
	linadd( 0., (realnum)i , "####", 'i',
		" start of general properties");

	/* this entry only works correctly if the APERTURE command is not in effect */
	if( geometry.iEmissPower == 2 )
	{
		linadd(continuum.totlsv/radius.dVeffAper,0,"Inci",'i',
		       "total luminosity in incident continuum");
		/* ipass is flag to indicate whether to only set up line array
		 * (ipass=0) or actually evaluate lines intensities (ipass=1) */
		if( LineSave.ipass > 0 )
		{
			continuum.totlsv = 0.;
		}
	}

	linadd(thermal.htot,0,"TotH",'i',
		"  total heating, all forms, information since individuals added later ");

	linadd(thermal.ctot,0,"TotC",'i',
		"  total cooling, all forms, information since individuals added later ");

	linadd(thermal.heating(0,0),0,"BFH1",'h',
		"  hydrogen photoionization heating, ground state only ");

	linadd(thermal.heating(0,1),0,"BFHx",'h',
		"  net hydrogen photoionization heating less rec cooling, all excited states normally zero, positive if excited states are net heating ");

	linadd(thermal.heating(0,22),0,"Line",'h',
		"  heating due to induced lines absorption of continuum ");
	if( thermal.htot > 0. )
	{
		if( thermal.heating(0,22)/thermal.htot > thermal.HeatLineMax )
		{
			thermal.HeatLineMax = (realnum)(thermal.heating(0,22)/thermal.htot);
		}
	}

	linadd(thermal.heating(1,0)+thermal.heating(1,1)+thermal.heating(1,2),0,"BFHe",'h',
	  "  total helium photoionization heating, all stages ");

	HeatMetal = 0.;
	/* some sums that will be printed in the stack */
	for( nelem=2; nelem<LIMELM; ++nelem)
	{
		/* we now have final solution for this element */
		for( i=dense.IonLow[nelem]; i < dense.IonHigh[nelem]; i++ )
		{
			ASSERT( i < LIMELM );
			/* total metal photo heating for LINES */
			HeatMetal += thermal.heating(nelem,i);
		}
	}

	linadd(HeatMetal,0,"TotM",'h',
		"  total heavy element photoionization heating, all stages ");

	linadd(thermal.heating(0,21),0,"pair",'h',
		"  heating due to pair production ");

	/* ipass is flag to indicate whether to only set up line array
	 * (ipass=0) or actually evaluate lines intensities (ipass=1) */
	if( LineSave.ipass > 0 )
	{
		/* this will be max local heating due to bound compton */
		ionbal.CompHeating_Max = MAX2( ionbal.CompHeating_Max , safe_div(ionbal.CompRecoilHeatLocal,thermal.htot,0.0));
	}
	else
	{
		ionbal.CompHeating_Max = 0.;
	}

	linadd(ionbal.CompRecoilHeatLocal,0,"Cbnd",'h',
		"  heating due to bound compton scattering ");

	linadd(rfield.cmheat,0,"ComH",'h',
		"  Compton heating ");

	linadd(CoolHeavy.tccool,0,"ComC",'c',
		"  total Compton cooling ");

	/* record max local heating due to advection */
	dynamics.HeatMax = MAX2( dynamics.HeatMax , safe_div(dynamics.Heat(),thermal.htot,0.) );
	/* record max local cooling due to advection */
	dynamics.CoolMax = MAX2( dynamics.CoolMax , safe_div(dynamics.Cool(),thermal.htot,0.) );

	linadd(dynamics.Cool()  , 0 , "advC" , 'i',
		"  cooling due to advection " );

	linadd(dynamics.Heat() , 0 , "advH" , 'i' ,
		"  heating due to advection ");

	linadd( thermal.char_tran_heat ,0,"CT H",'h',
		" heating due to charge transfer ");

	linadd( thermal.char_tran_cool ,0,"CT C",'c',
		" cooling due to charge transfer ");

	linadd(thermal.heating(1,6),0,"CR H",'h',
		" cosmic ray heating ");

	linadd(thermal.heating(0,20),0,"extH",'h',
		" extra heat added to this zone, from HEXTRA command ");

	linadd(CoolHeavy.cextxx,0,"extC",'c',
		" extra cooling added to this zone, from CEXTRA command ");

	// 511 keV annihilation line, counts as recombination line since
	// neither heating nor cooling, but does remove energy
	ee511 = (dense.gas_phase[ipHYDROGEN] + 4.*dense.gas_phase[ipHELIUM])*ionbal.PairProducPhotoRate[0]*2.*8.20e-7;
	PntForLine(2.427e-2,"e-e+",&ipnt);
	lindst(ee511,(realnum)2.427e-2,"e-e+",ipnt,'r',true,
		" 511keV annihilation line " );

	linadd(CoolHeavy.expans,0,"Expn",'c',
		"  expansion cooling, only non-zero for wind ");

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].RadRecCool,0,"H FB",'i',
		"  H radiative recombination cooling ");

	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].FreeBnd_net_Cool_Rate),0,"HFBc",'c',
		"  net free-bound cooling ");

	linadd(MAX2(0.,-iso_sp[ipH_LIKE][ipHYDROGEN].FreeBnd_net_Cool_Rate),0,"HFBh",'h',
		"  net free-bound heating ");

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].RecomInducCool_Rate,0,"Hind",'c',
		"  cooling due to induced rec of hydrogen ");

	linadd(CoolHeavy.cyntrn,0,"Cycn",'c',
		"  cyclotron cooling ");

	// cooling due to iso-sequence species
	for( size_t ipISO = ipH_LIKE; ipISO <= size_t(ipHE_LIKE); ipISO++ )
	{
		for( size_t nelem = ipISO; nelem < size_t(LIMELM); nelem++ )
		{
			string chLabel_base = chIonLbl( nelem+1, nelem+1-ipISO );
			// label may be too long for linadd
			size_t nchars = min( NCHLAB-1, chLabel_base.length() );
			char chLabel_cool[NCHLAB] = { 0 },
				chLabel_heat[NCHLAB] = { 0 };
			strncpy( chLabel_cool, chLabel_base.c_str(), NCHLAB-1 );
			strncpy( chLabel_heat, chLabel_base.c_str(), NCHLAB-1 );
			if( nchars <  NCHLAB-1 )
			{
				chLabel_cool[ nchars ] = 'c';
				chLabel_heat[ nchars ] = 'h';
			}
			// this is information, 'i', since individual lines
			// have been added as cooling or heating
			linadd(max(0., iso_sp[ipISO][nelem].xLineTotCool),0, chLabel_cool,'c',
				" net cooling due to iso-seq species");
			linadd(max(0., -iso_sp[ipISO][nelem].xLineTotCool),0, chLabel_heat,'h',
				" heating due to iso-seq species");
		}
	}

	// cooling due to database species
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		// label may be too long for linadd
		size_t nchars = min( NCHLAB-1, strlen( dBaseSpecies[ipSpecies].chLabel ) );
		char chLabel_cool[NCHLAB] = { 0 },
			chLabel_heat[NCHLAB] = { 0 };
		strncpy( chLabel_cool, dBaseSpecies[ipSpecies].chLabel, NCHLAB-1 );
		strncpy( chLabel_heat, dBaseSpecies[ipSpecies].chLabel, NCHLAB-1 );
		if( nchars <  NCHLAB-1 )
		{
			chLabel_cool[ nchars ] = 'c';
			chLabel_heat[ nchars ] = 'h';
		}
		// this is information, 'i', since individual lines
		// have been added as cooling or heating
		linadd(max(0., dBaseSpecies[ipSpecies].CoolTotal),0, chLabel_cool,'c',
			" net cooling due to database species");
		linadd(max(0., -dBaseSpecies[ipSpecies].CoolTotal),0, chLabel_heat,'h',
			" heating due to database species");
	}

	return;
}

