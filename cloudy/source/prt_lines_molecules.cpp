/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_molecules put energetics, H, and He lines into line intensity stack */
#include "cddefines.h"
#include "coolheavy.h"
#include "thermal.h"
#include "dense.h"
#include "hmi.h"
#include "phycon.h"
#include "h2.h"
#include "co.h"
#include "radius.h"
#include "lines.h"
#include "lines_service.h"
#include "mole.h"

void lines_molecules(void)
{
	long int i;

	DEBUG_ENTRY( "lines_molecules()" );

	/* molecules */
	i = StuffComment( "molecules" );
	linadd( 0., (realnum)i , "####", 'i',
		"  molecules");


	/* >>refer	H2	rot	Lepp, S., & Shull, J.M., 1983, ApJ, 270, 578-582
	 * roughly two microns */
	linadd(CoolHeavy.h2line,20000.,"H2 l",'c', 
		   "cooling due H2 rotation lines from simple model" );
	/* remember largest fraction of H2 cooling for possible comment */
	hmi.h2line_cool_frac = (realnum)MAX2( CoolHeavy.h2line/thermal.ctot , hmi.h2line_cool_frac );

	/* HD rotation cooling */
	linadd(CoolHeavy.HD,0,"HDro",'c',
		"HD rotation cooling");

	/* molecular hydrogen heating */
	hmi.h2dtot += (realnum)(hmi.HeatH2Dish_used*radius.dVeffAper);
	hmi.h2dfrc = (realnum)(MAX2(safe_div(hmi.HeatH2Dish_used,thermal.htot,0.0),hmi.h2dfrc));

	/* largest fraction of heating due to photo dissoc of H2+ */
	hmi.h2pmax = MAX2(hmi.h2pmax,(realnum)(safe_div(thermal.heating(0,16),thermal.htot,0.0)));

	linadd(hmi.HeatH2Dish_used,0,"H2dH",'h',
		"heating by H2 dissociation by photons and cosmic rays");

	/*remember largest fraction of heating due to H2 vib deexcitation */	
	hmi.HeatH2DexcMax = MAX2(safe_div((realnum)hmi.HeatH2Dexc_used,(realnum)thermal.htot,(realnum)0.0f),hmi.HeatH2DexcMax);

	/*remember largest fraction of cooling due to H2 cooling */
	hmi.CoolH2DexcMax = MAX2(safe_div(-(realnum)hmi.HeatH2Dexc_used,(realnum)thermal.htot,(realnum)0.0f),hmi.CoolH2DexcMax);

	linadd( MAX2(0.,hmi.HeatH2Dexc_used),0,"H2vH",'h',
		"heating by coll deexcit of vib-excited H2");

	linadd( MAX2(0.,-hmi.HeatH2Dexc_used) ,0,"H2vC",'c',
		" cooling by coll deexcit of vib-excited H2");

	/* line emission by vib-excited H2 */
	if( h2.lgEnabled )
	{

		linadd( 0. ,0,"H2 v",'i',
			"  when large molecule is turned on do not print this simple estimate  line emission by vib-excited H2 ");
	}
	else
	{
		linadd( findspecieslocal("H2*")->den*2e-7*4.17e-12,0,"H2 v",'i',
			" H2 vib-excited lines from Tielens & Hollenbach 1985");
	}

	/* add in explicit lines from the large H2 molecule 
	 * routine in mole_h2_io.c */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_LinesAdd();

	linadd(hmi.hmicol,0,"H-FB",'c',
		"	 neg H ion free-bound emission, H + e -> H- + hnu ");

	linadd(CoolHeavy.brems_cool_hminus,0,"H-FF",'i',
		" neg H ion free-free emission ");

	/* H-alpha produced by H- mutual neutralization */
	linadd(mole.findrate("H-,H+=>H,H")*3.032e-12,6562.80,"H-CT",'i',
		"  H-alpha produced by H- mutual neutralization ");

	/* remember total heating */
	hmi.hmitot += hmi.hmihet*radius.dVeffAper;

	linadd(MAX2(0.,hmi.hmihet),0,"H- H",'h',
		"  H- heating ");

	linadd(MAX2(0.,-hmi.hmihet),0,"H-Hc",'c',
		"  induced H- cooling ");

	linadd(CoolHeavy.H2PlsCool,0,"H2+ ",'c',
		"  H+ + H => H2+ + photon continuum cooling ");

	linadd(hmi.h2plus_heat,0,"H2+p",'h',
		"  H2+ photo dissoc heating ");

	linadd(MAX2(3.27e-12+phycon.te*BOLTZMANN,0.)*dense.xIonDense[ipHYDROGEN][1]*dense.xIonDense[ipHELIUM][0]*1e-20+
	  (1.76e-11+phycon.te*BOLTZMANN)*dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipHELIUM][1]*1e-16,0,"HEH+",'i'	,
	  "  HeH+ formation cooling ");

	/* carbon monoxide heating */
	co.codtot += co.CODissHeat*(realnum)radius.dVeffAper;
	co.codfrc = (realnum)MAX2(safe_div((double)co.CODissHeat,thermal.htot,0.0),co.codfrc);

	linadd(co.CODissHeat,0,"COdh",'h',
		"  carbon monoxide co photodissociation ");

	return;
}
