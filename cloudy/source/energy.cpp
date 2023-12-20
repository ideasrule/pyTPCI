/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "continuum.h"
#include "rfield.h"
#include "doppvel.h"
#include "geometry.h"
#include "hextra.h"
#include "ipoint.h"
#include "lines.h"
#include "opacity.h"
#include "radius.h"
#include "secondaries.h"
#include "struc.h"
#include "thermal.h"
#include "warnings.h"
#include "wind.h"
#include "dynamics.h"
#include "lines_service.h"

static const char *ENERGY_RYD    = "Ryd";
static const char *ENERGY_ERG    = "erg";
static const char *ENERGY_EV     = "eV";
static const char *ENERGY_KEV    = "keV";
static const char *ENERGY_MEV    = "MeV";
static const char *ENERGY_WN     = "cm^-1";
static const char *ENERGY_CM     = "cm";
static const char *ENERGY_MM     = "mm";
static const char *ENERGY_MICRON = "um";
static const char *ENERGY_NM     = "nm";
static const char *ENERGY_A      = "A";
static const char *ENERGY_HZ     = "Hz";
static const char *ENERGY_KHZ    = "kHz";
static const char *ENERGY_MHZ    = "MHz";
static const char *ENERGY_GHZ    = "GHz";
static const char *ENERGY_K      = "K";

inline bool isSameUnit(const char *unit, const char *ref)
{
	return strcmp(unit,ref) == 0;
}

const char *StandardEnergyUnitNoAbort(const char *chCard)
{
	DEBUG_ENTRY( "StandardEnergyUnit()" );

	// use " MIC" so that it cannot pick up on "W/SQM/MICRON" in flux units
	if( nMatch(" MIC",chCard) )
	{
		/* micron */
		return ENERGY_MICRON;
	}
	else if( nMatch(" EV ",chCard) )
	{
		/* eV */
		return ENERGY_EV;
	}
	else if( nMatch(" KEV",chCard) )
	{
		/* keV */
		return ENERGY_KEV;
	}
	else if( nMatch(" MEV",chCard) )
	{
		/* MeV */
		return ENERGY_MEV;
	}
	else if( nMatch("WAVE",chCard) )
	{
		/* wavenumbers */
		return ENERGY_WN;
	}
	else if( nMatch("CENT",chCard) || nMatch(" CM ",chCard) )
	{
		/* centimeter */
		return ENERGY_CM;
	}
	else if( nMatch(" MM ",chCard) )
	{
		/* millimeter */
		return ENERGY_MM;
	}
	else if( nMatch(" NM ",chCard) )
	{
		/* nanometer */
		return ENERGY_NM;
	}
	else if( nMatch("ANGS",chCard) )
	{
		/* angstrom */
		return ENERGY_A;
	}
	else if( nMatch(" HZ ",chCard) )
	{
		/* hertz */
		return ENERGY_HZ;
	}
	else if( nMatch(" KHZ",chCard) )
	{
		/* kilohertz */
		return ENERGY_KHZ;
	}
	else if( nMatch(" MHZ",chCard) )
	{
		/* megahertz */
		return ENERGY_MHZ;
	}
	else if( nMatch(" GHZ",chCard) )
	{
		/* gigahertz */
		return ENERGY_GHZ;
	}
	else if( nMatch("KELV",chCard) || nMatch(" K ",chCard) )
	{
		/* kelvin */
		return ENERGY_K;
	}
	else if( nMatch(" RYD",chCard) )
	{
		/* RydbERG must come before ERG to avoid false trip */
		return ENERGY_RYD;
	}
	// use "ERG " so that it cannot pick up on "ERG/S" in flux units
	else if( nMatch(" ERG ",chCard) )
	{
		/* erg */
		return ENERGY_ERG;
	}
	else
	{
		/* not recognized */
		return nullptr;
	}
}

double Energy::get(const char *unit) const
{
	DEBUG_ENTRY( "Energy::get()" );

	/* internal unit is Ryd */
	if( isSameUnit(unit,ENERGY_RYD) )
	{
		return Ryd();
	}
	else if( isSameUnit(unit,ENERGY_ERG) )
	{
		return Erg();
	}
	else if( isSameUnit(unit,ENERGY_EV) )
	{
		return eV();
	}
	else if( isSameUnit(unit,ENERGY_KEV) )
	{
		return keV();
	}
	else if( isSameUnit(unit,ENERGY_MEV) )
	{
		return MeV();
	}
	else if( isSameUnit(unit,ENERGY_WN) )
	{
		return WN();
	}
	else if( isSameUnit(unit,ENERGY_CM) )
	{
		return cm();
	}
	else if( isSameUnit(unit,ENERGY_MM) )
	{
		return mm();
	}
	else if( isSameUnit(unit,ENERGY_MICRON) )
	{
		return micron();
	}
	else if( isSameUnit(unit,ENERGY_NM) )
	{
		return nm();
	}
	else if( isSameUnit(unit,ENERGY_A) )
	{
		return Angstrom();
	}
	else if( isSameUnit(unit,ENERGY_HZ) )
	{
		return Hz();
	}
	else if( isSameUnit(unit,ENERGY_KHZ) )
	{
		return kHz();
	}
	else if( isSameUnit(unit,ENERGY_MHZ) )
	{
		return MHz();
	}
	else if( isSameUnit(unit,ENERGY_GHZ) )
	{
		return GHz();
	}
	else if( isSameUnit(unit,ENERGY_K) )
	{
		return K();
	}
	else
	{
		fprintf( ioQQQ, " insane units in Energy::get: \"%s\"\n", unit );
		cdEXIT(EXIT_FAILURE);
	}
}

void Energy::set(double value, const char *unit)
{
	DEBUG_ENTRY( "Energy::set()" );

	/* internal unit is Ryd */
	if( isSameUnit(unit,ENERGY_RYD) )
	{
		set(value);
	}
	else if( isSameUnit(unit,ENERGY_ERG) )
	{
		set(value/EN1RYD);
	}
	else if( isSameUnit(unit,ENERGY_MEV) )
	{
		set(1e6*value/EVRYD);
	}
	else if( isSameUnit(unit,ENERGY_KEV) )
	{
		set(1e3*value/EVRYD);
	}
	else if( isSameUnit(unit,ENERGY_EV) )
	{
		set(value/EVRYD);
	}
	else if( isSameUnit(unit,ENERGY_WN) )
	{
		set(value/RYD_INF);
	}
	else if( isSameUnit(unit,ENERGY_A) )
	{
		set(RYDLAM/value);
	}
	else if( isSameUnit(unit,ENERGY_NM) )
	{
		set(RYDLAM/(1.e1*value));
	}
	else if( isSameUnit(unit,ENERGY_MICRON) )
	{
		set(RYDLAM/(1.e4*value));
	}
	else if( isSameUnit(unit,ENERGY_MM) )
	{
		set(RYDLAM/(1.e7*value));
	}
	else if( isSameUnit(unit,ENERGY_CM) )
	{
		set(RYDLAM/(1.e8*value));
	}
	else if( isSameUnit(unit,ENERGY_HZ) )
	{
		set(value/FR1RYD);
	}
	else if( isSameUnit(unit,ENERGY_KHZ) )
	{
		set(1e3*value/FR1RYD);
	}
	else if( isSameUnit(unit,ENERGY_MHZ) )
	{
		set(1e6*value/FR1RYD);
	}
	else if( isSameUnit(unit,ENERGY_GHZ) )
	{
		set(1e9*value/FR1RYD);
	}
	else if( isSameUnit(unit,ENERGY_K) )
	{
		set(value/TE1RYD);
	}
	else
	{
		fprintf( ioQQQ, " insane units in Energy::set: \"%s\"\n", unit );
		cdEXIT(EXIT_FAILURE);
	}
}

void EnergyEntry::p_set_ip()
{
	DEBUG_ENTRY( "EnergyEntry::p_set_ip()" );

	double energy = Ryd();
	if( energy < rfield.emm() || energy > rfield.egamry() )
	{
		fprintf( ioQQQ, " The energy %g Ryd is not within the default Cloudy range\n", energy );
		cdEXIT(EXIT_FAILURE);
	}
	p_ip = ipoint(energy) - 1;
	ASSERT( p_ip >= 0 );
}

STATIC double sum_radiation( void );

/*badprt print out coolants if energy not conserved */
STATIC void badprt(
		/* total is total luminosity available in radiation */
		double total);

void ConserveEnergy()
{
	DEBUG_ENTRY( "ConserveEnergy()" );

	/* check whether energy is conserved
	 * following is outward continuum */
	double outin, outout, outtot;
	outsum(&outtot,&outin,&outout);
	/* sum_recom_lines is the sum of all recombination line energy */
	double sum_recom_lines = totlin('r');
	if( sum_recom_lines == 0. )
		bangin("  !Total recombination line energy is 0.");

	/* sum_coolants is the sum of all coolants */
	double sum_coolants = totlin('c');
	if( sum_coolants == 0. )
		bangin("  !Total cooling is zero.");

	double flow_energy;
	if( !(wind.lgBallistic() || wind.lgStatic()) )
	{
		/* approximate energy density coming out in wind
		 * should ideally include flow into grid and internal energies */
		flow_energy = (2.5*struc.GasPressure[0]+0.5*struc.DenMass[0]*wind.windv0*wind.windv0)
			*(-wind.windv0);
	}
	else
	{
		flow_energy = 0.;
	}

	if( 0 && geometry.lgSphere && (!thermal.lgTemperatureConstant) && (!dynamics.lgTimeDependentStatic) &&
		(hextra.cryden == 0.) && ((hextra.TurbHeat+DoppVel.DispScale) == 0.) && !secondaries.lgCSetOn )
	{
		double sum_rad = sum_radiation();

		if( fabs( 1. - sum_rad/(continuum.TotalLumin*POW2(radius.rinner)) ) > 0.01 )
			fprintf( ioQQQ, "PROBLEM DISASTER ConserveEnergy failed %li\t%li\t%e\t%e\t%e\n",
				iteration, nzone, sum_rad, continuum.TotalLumin*POW2(radius.rinner),
					sum_rad/(continuum.TotalLumin*POW2(radius.rinner)) );
	}

	/* TotalLumin is total incident photon luminosity, evaluated in setup
	 * sum_coolants is evaluated here, is sum of all coolants */
	/* this test is meaningless when the APERTURE command is in effect since
	 * sum_coolants and sum_recom_lines react to the APERTURE command, while
	 * continuum.TotalLumin does not */
	if( !dynamics.lgTimeDependentStatic && (sum_coolants + sum_recom_lines + flow_energy) > continuum.TotalLumin*geometry.covgeo &&
	    !thermal.lgTemperatureConstant && geometry.iEmissPower == 2 )
	{
		if( (hextra.cryden == 0.) && ((hextra.TurbHeat+DoppVel.DispScale) == 0.) && !secondaries.lgCSetOn )
		{
			ostringstream oss1;
			oss1 << " W-Radiated luminosity (cool+rec+wind=" << scientific << setprecision(2);
			oss1 << sum_coolants << "+" << sum_recom_lines << "+" << flow_energy;
			oss1 << ") is greater than that in incident radiation (TotalLumin=";
			oss1 << setw(8) << continuum.TotalLumin << ").  Power radiated was ";
			oss1 << setw(8) << thermal.power << ".";
			warnin(oss1.str());
			/* write same thing directly to output (above will be sorted later) */
			fprintf( ioQQQ, "\n\n PROBLEM DISASTER This calculation DID NOT CONSERVE ENERGY!  This is zone %li\n\n\n",
					nzone);

			if( !continuum.lgCheckEnergyEveryZone )
				fprintf( ioQQQ, "Rerun with *set check energy every zone* command to do"
					 " energy check for each zone.\n\n" );

			/* the case b command can cause this problem - say so if case b was set */
			if( opac.lgCaseB )
				fprintf( ioQQQ, "\n The CASE B command was entered - this can have unphysical effects,"
					 " including non-conservation of energy.\n Why was it needed?\n\n" );

			/* print out significant heating and cooling */
			badprt(continuum.TotalLumin);

			ostringstream oss2;
			oss2 <<  " W-Something is really wrong: the ratio of radiated to incident luminosity is ";
			oss2 << scientific << setprecision(2) << (sum_coolants + sum_recom_lines)/continuum.TotalLumin;
			oss2 << ".";
			warnin(oss2.str());

			/* this can be caused by the force te command */
			if( thermal.ConstTemp > 0. )
			{
				fprintf( ioQQQ, " This may have been caused by the FORCE TE command.\n" );
				fprintf( ioQQQ, " Remove it and run again.\n" );
			}
			else
			{
				ShowMe();
				throw cloudy_abort("the calculation did not conserve energy");
			}
		}
	}
}

STATIC double sum_radiation( void )
{
	double flxref = 0., flxatt = 0., conem = 0.;
	long nEmType = 0;
	double sum = 0.;
	
	const realnum *trans_coef_total=rfield.getCoarseTransCoef();
	for( long j=0; j<rfield.nflux; j++ )
	{
		/* the reflected continuum */
		flxref += (rfield.ConRefIncid[nEmType][j] + rfield.ConEmitReflec[nEmType][j] + rfield.reflin[nEmType][j]) * rfield.anu(j);
		ASSERT( flxref >= 0. );

		/* the attenuated incident continuum */
		flxatt += rfield.flux[nEmType][j]*radius.r1r0sq*trans_coef_total[j] * rfield.anu(j);
		ASSERT( flxatt >= 0. );

		/* the outward emitted continuum */
		conem += (rfield.ConEmitOut[nEmType][j] + rfield.outlin[nEmType][j] + rfield.outlin_noplot[j])*radius.r1r0sq*geometry.covgeo * rfield.anu(j);
		ASSERT( conem >= 0. );
	}

	sum = (flxref + flxatt + conem) * EN1RYD * POW2(radius.rinner);
	ASSERT( sum >= 0. );

	return sum;
}

/*badprt print out coolants if energy not conserved */
STATIC void badprt(
		/* total is total luminosity available in radiation */
		double total)
{
	/* following is smallest ratio to print */
	static const double RATIO = 0.02;
	char chInfo;
	long int i;
	realnum sum_coolants, 
	  sum_recom_lines;

	DEBUG_ENTRY( "badprt()" );

	fprintf( ioQQQ, " badprt: all entries with greater than%6.2f%% of incident continuum follow.\n", 
	  RATIO*100. );
	fprintf( ioQQQ, " badprt: Intensities are relative to total energy in incident continuum.\n" );

	/* now find sum of recombination lines */
	chInfo = 'r';
	sum_recom_lines = (realnum)totlin('r');
	fprintf( ioQQQ, 
		" Sum of energy in recombination lines is: %.2e ; relative to total incident radiation is: %.2e\n", 
		sum_recom_lines, 
		sum_recom_lines/MAX2(1e-30,total) );

	fprintf(ioQQQ," all strong information lines \n line  wl  ener/total\n");
	/* now print all strong lines */
	for( i=0; i < LineSave.nsum; i++ )
	{
		if( LineSave.lines[i].chSumTyp() == chInfo && LineSave.lines[i].SumLine(0)/total > RATIO )
		{
			fprintf( ioQQQ, " ");
			LineSave.lines[i].prt(ioQQQ);
			fprintf( ioQQQ, " %7.3f %c\n", LineSave.lines[i].SumLine(0)/total, chInfo );
		}
	}

	fprintf(ioQQQ," all strong cooling lines \n line  wl  ener/total\n");
	chInfo = 'c';
	sum_coolants = (realnum)totlin('c');
	fprintf( ioQQQ, " Sum of coolants (abs) = %.2e (rel)= %.2e\n", 
	  sum_coolants, sum_coolants/MAX2(1e-30,total) );
	for( i=0; i < LineSave.nsum; i++ )
	{
		if( LineSave.lines[i].chSumTyp() == chInfo && LineSave.lines[i].SumLine(0)/total > RATIO )
		{
			fprintf( ioQQQ, " ");
			LineSave.lines[i].prt(ioQQQ);
			fprintf( ioQQQ, " %7.3f %c\n", LineSave.lines[i].SumLine(0)/total, chInfo );
		}
	}

	fprintf(ioQQQ," all strong heating lines \n line  wl  ener/total\n");
	chInfo = 'h';
	fprintf( ioQQQ, " Sum of heat (abs) = %.2e (rel)= %.2e\n", 
	  thermal.power, thermal.power/MAX2(1e-30,total) );
	for( i=0; i < LineSave.nsum; i++ )
	{
		if( LineSave.lines[i].chSumTyp() == chInfo && LineSave.lines[i].SumLine(0)/total > RATIO )
		{
			fprintf( ioQQQ, " ");
			LineSave.lines[i].prt(ioQQQ);
			fprintf( ioQQQ, " %7.3f %c\n", LineSave.lines[i].SumLine(0)/total, chInfo );
		}
	}

	return;
}
