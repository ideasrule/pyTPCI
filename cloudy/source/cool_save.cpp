/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolSave save coolants */
#include "cddefines.h"
#include "thermal.h"
#include "dynamics.h"
#include "radius.h"
#include "conv.h"
#include "phycon.h"
#include "save.h"
#include "grainvar.h"
#include "hmi.h"
#include "coolheavy.h"

/*CoolSave save coolants */
void CoolSave(FILE * io, const char chJob[])
{
	long int i, 
	  ip, 
	  is;

	int nFail;

	double cset, 
		cool_total,
		heat_total;

	DEBUG_ENTRY( "CoolSave()" );

	/* cannot do one-time init since thermal.ncltot can change */
	vector<long>index(thermal.ncltot);
	vector<realnum>csav(thermal.ncltot), 
		sgnsav(thermal.ncltot);

	cool_total = thermal.ctot;
	heat_total = thermal.htot;

	/* >>chng 06 mar 17, comment out following block and replace with this 
	 * removing dynamics heating & cooling and report only physical
	 * heating and cooling 
	 * NB the heating and cooling as punched no longer need be
	 * equal for a converged model */
	cool_total -= dynamics.Cool();
	heat_total -= dynamics.Heat();
#	if 0
	if(dynamics.Cool > dynamics.Heat()) 
	{
		cool_total -= dynamics.Heat();
		heat_total -= dynamics.Heat();
	} 
	else
	{
		cool_total -= dynamics.Cool;
		heat_total -= dynamics.Cool;
	}
#	endif

	/* cset will be weakest cooling to consider
	 * WeakHeatCool set with 'set weakheatcool' command
	 * default is 0.05 */
	cset = cool_total*save.WeakHeatCool;

	/* first find all strong lines, both + and - sign */
	for( i=0; i < thermal.ncltot; i++ )
	{
		/* NB NB exclude advective cooling from sorting */
		if( strncmp( thermal.chClntLab[i], "adve", 4 ) == 0 )
		{
			csav[i] = 0;
			continue;
		}

		if( cool_total>SMALLFLOAT )
			csav[i] = (realnum)( safe_div( MAX2(thermal.cooling[i],thermal.heatnt[i]), cool_total, 0. ));
		else
			csav[i] = 0.;
		/* save sign to remember if heating or cooling line */
		if( thermal.heatnt[i] == 0. )
		{
			sgnsav[i] = 1.;
		}
		else
		{
			sgnsav[i] = -1.;
		}
	}

	ip = thermal.ncltot;

	/* order strongest to weakest */
	/* now sort by decreasing importance */
	/*spsort netlib routine to sort array returning sorted indices */
	spsort(
		  /* input array to be sorted */
		  &csav[0], 
		  /* number of values in x */
		  ip, 
		  /* permutation output array */
		  &index[0], 
		  /* flag saying what to do - 1 sorts into increasing order, not changing
		   * the original routine */
		  -1, 
		  /* error condition, should be 0 */
		  &nFail);

	/* warn if tcovergence failure occurred */
	if( !conv.lgConvTemp )
	{
		fprintf( io, "#>>>>  Temperature not converged.\n" );
	}
	else if( !conv.lgConvEden )
	{
		fprintf( io, "#>>>>  Electron density not converged.\n" );
	}
	else if( !conv.lgConvIoniz() )
	{
		fprintf( io, "#>>>>  Ionization not converged.\n" );
	}
	else if( !conv.lgConvPres )
	{
		fprintf( io, "#>>>>  Pressure not converged.\n" );
	}

	if( strcmp(chJob,"EACH") == 0 )
	{
		/* begin the print out with zone number, total heating and cooling */
		fprintf( io, "%.5e\t%.4e\t%.4e", 
				 radius.depth_mid_zone, 
				 phycon.te, 
				 cool_total );
		double debug_ctot = 0.;

		for( int i=0 ; i < LIMELM ; i++)
		{
			fprintf( io, "\t%.4e", thermal.elementcool[i]+thermal.heavycollcool[i]);
			debug_ctot += MAX2(0., thermal.elementcool[i]+thermal.heavycollcool[i]);
		}

		/* molecular cooling, except those who are listed separately */
		fprintf( io, "\t%.4e", thermal.elementcool[LIMELM]);
		debug_ctot += MAX2( 0., thermal.elementcool[LIMELM] );

		/* dust */
		fprintf( io, "\t%.4e", MAX2(0.,gv.GasCoolColl) );
		debug_ctot += MAX2(0.,gv.GasCoolColl);

		/* H2cX - H2 molecule cooling*/
		fprintf( io, "\t%.4e", MAX2(0.,-hmi.HeatH2Dexc_used) );
		debug_ctot += MAX2(0.,-hmi.HeatH2Dexc_used);

		/* CT C - charge transfering cooling*/ 
		fprintf( io, "\t%.4e", thermal.char_tran_cool );
		debug_ctot += MAX2( 0., thermal.char_tran_cool );

		/* H-fb - H + e -> H- + h*niu*/
		fprintf( io, "\t%.4e", hmi.hmicol );
		debug_ctot += MAX2( 0., hmi.hmicol );

		/* H2ln - line cooling within simple H2 molecule (rotation) */
		fprintf( io, "\t%.4e", CoolHeavy.h2line );
		debug_ctot += MAX2( 0., CoolHeavy.h2line );

		/* HDro - HD cooling*/
		fprintf( io, "\t%.4e", CoolHeavy.HD );
		debug_ctot += MAX2( 0., CoolHeavy.HD );

		/* H2+ - H + H+ -> H2+ cooling*/
		fprintf( io, "\t%.4e", CoolHeavy.H2PlsCool );
		debug_ctot += MAX2( 0., CoolHeavy.H2PlsCool );

		/* FFcm, splite into FF_H and FF_M */
		/*fprintf( io, "\t%.4e", MAX2(0.,CoolHeavy.brems_cool_net) );
		  debug_ctot += MAX2(0.,CoolHeavy.brems_cool_net); */

		/* Due to the calculation of bremsstrahlung heating, FF_H and FF_M result may not reliable when CoolHeavy.brems_cool_metals is not far more than bremsstrahlung heating */
		if( CoolHeavy.brems_cool_net > 0. )
		{
			/* FF_H, free-free cooling from H and He */
			fprintf( io, "\t%.4e", CoolHeavy.brems_cool_h + CoolHeavy.brems_cool_he + CoolHeavy.brems_cool_hminus );
			debug_ctot += MAX2( 0., CoolHeavy.brems_cool_h + CoolHeavy.brems_cool_he + CoolHeavy.brems_cool_hminus );
			/* FF_M, free-free cooling from metal */
			fprintf( io, "\t%.4e", CoolHeavy.brems_cool_metals - CoolHeavy.brems_heat_total );
			debug_ctot += MAX2( 0., CoolHeavy.brems_cool_metals - CoolHeavy.brems_heat_total );
		}
		else
		{
			double zero = 0.;
			fprintf( io, "\t%.4e", zero );
			fprintf( io, "\t%.4e", zero );
		}

		/* NB - heavy element recombination cooling (hvFB, CoolHeavy.heavfb) has been splited into each element*/

		/* eeff - electron electron bremsstrahlung cooling */
		fprintf( io, "\t%.4e", CoolHeavy.eebrm );
		debug_ctot += MAX2( 0., CoolHeavy.eebrm );

		/* adve - dynamics cooling, do not add this to debug_ctot because cool_total does not contain this term */
		fprintf( io, "\t%.4e", dynamics.Cool() );

		/* comp - Compton cooling*/
		fprintf( io, "\t%.4e", CoolHeavy.tccool );
		debug_ctot += MAX2( 0., CoolHeavy.tccool );

		/* Extr - extra cooling, set with CEXTRA command*/
		fprintf( io, "\t%.4e", CoolHeavy.cextxx );
		debug_ctot += MAX2( 0., CoolHeavy.cextxx );

		/* Expn - wind expansion cooling*/
		fprintf( io, "\t%.4e", CoolHeavy.expans );
		debug_ctot += MAX2( 0., CoolHeavy.expans );

		/* Cycl - cyclotron cooling*/
		fprintf( io, "\t%.4e", CoolHeavy.cyntrn );
		debug_ctot += MAX2( 0., CoolHeavy.cyntrn );

		/* Hvin - heavy element collisional ionization cooling (CoolHeavy.colmet), splited into each element */

		/* dima - report, but don't add Dima cooling;
		 * already included in elementcool through atom_level2() */
		fprintf( io, "\t%.4e", thermal.dima );

		fprintf( io, " \n" );

		{
			enum{ DEBUG_COOLING = false };
			if( DEBUG_COOLING )
			{
				fprintf( ioQQQ, "DEBUG COOLING:\t"
						"recomputed: %6e\t"
						"should be: %.6e\t"
						"fractional diff: %.4e\n",
						debug_ctot,
						cool_total,
						debug_ctot / cool_total - 1. );
			}
		}

		/* check if all coolants are added together */
		if( fabs( (debug_ctot - cool_total)/cool_total ) > 1e-10 )
		{
			fprintf( ioQQQ , "PROBLEM with the SAVE EACH COOLING output\n" );
			fprintf( ioQQQ , "PROBLEM One or more coolants have been lost, the sum of the reported cooling is %.4e\n", debug_ctot );
			fprintf( ioQQQ , "PROBLEM The total cooling is %.4e\n", cool_total );
			fprintf( ioQQQ , "PROBLEM The difference is %.4e\n", cool_total - debug_ctot );
			cdEXIT(EXIT_FAILURE);
		}
	}
	else if( strcmp(chJob,"COOL") == 0 )
	{
		/*>>chng 06 jun 06, change start of save to give same info as heating 
		 * as per comment by Yumihiko Tsuzuki */
		/* begin the print out with zone number, total heating and cooling */
		fprintf( io, "%.5e\t%.4e\t%.4e\t%.4e\t%.4e", 
				radius.depth_mid_zone, 
				phycon.te, 
				heat_total, 
				cool_total,
				dynamics.Cool() );

		/* now print the coolants 
		 * keep sign of coolant, for strong negative cooling 
		 * order is ion, wavelength, fraction of total */
		for( is=0; is < ip; is++ )
		{
			if(is > 4 && (thermal.cooling[index[is]] < cset && thermal.heatnt[index[is]] < cset))
				break;
			fprintf( io, "\t%s %.1f\t%.7f", 
					 thermal.chClntLab[index[is]], 
					 thermal.collam[index[is]], 
					 sign(csav[index[is]],sgnsav[index[is]]) );
		}
		fprintf( io, " \n" );
	}
	else
		TotalInsanity();

	return;
}
