/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveHeat save contributors to local heating, with save heat command, called by punch_do */
#include "cddefines.h"
#include "thermal.h"
#include "radius.h"
#include "conv.h"
#include "dense.h"
#include "taulines.h"
#include "phycon.h"
#include "elementnames.h"
#include "dynamics.h"
#include "save.h"
#include "lines_service.h"

/*SaveHeat save contributors to local heating, with save heat command, 
 * called by punch_do */
void SaveHeat(FILE* io)
{
	DEBUG_ENTRY( "SaveHeat()" );

	vector<realnum>SaveVal(LIMELM*LIMELM, -FLT_MAX);
	vector<long>ipsave(LIMELM*LIMELM, INT_MIN);
	vector<long>jpsave(LIMELM*LIMELM, INT_MIN);
	vector<long>ipOrdered(LIMELM*LIMELM);
	vector<string> chLabel(LIMELM*LIMELM);

	double cool_total = thermal.ctot;
	double heat_total = thermal.htot;

	/* >>chng 06 mar 17, comment out following block and replace with this 
	 * removing dynamics heating & cooling and report only physical
	 * heating and cooling 
	 * NB the heating and cooling as punched no longer need be
	 * equal for a converged model */
	cool_total -= dynamics.Cool();
	heat_total -= dynamics.Heat();
	if( false )
	{
		if(dynamics.Cool() > dynamics.Heat()) 
		{
			cool_total -= dynamics.Heat();
			heat_total -= dynamics.Heat();
		} 
		else
		{
			cool_total -= dynamics.Cool();
			heat_total -= dynamics.Cool();
		}
	}

	long ipnt = 0;

	/* heat sources are saved in a 2d square array */
	/* WeakHeatCool set with 'set weakheatcool' command
	 * default is 0.05 */
	for( long i=0; i < LIMELM; i++ )
	{
		for( long j=0; j < LIMELM; j++ )
		{
			/* Skip advective heating */
			if( i == 1 && j == 5 )
				continue;
			if( safe_div( thermal.heating(i,j), heat_total, 0. ) > SMALLFLOAT )
			{
				ipsave[ipnt] = i;
				jpsave[ipnt] = j;
				SaveVal[ipnt] = (realnum)( safe_div( thermal.heating(i,j), heat_total, 0. ) );
				ipnt++;
			}
		}
	}

	/* now check for possible line heating not counted in 1,23
	 * there should not be any significant heating source here
	 * since they would not be counted in derivative correctly */
	for( long i=0; i < thermal.ncltot; i++ )
	{
		if( safe_div( thermal.heatnt[i], heat_total, 0. ) > save.WeakHeatCool )
		{
			realnum awl;
			awl = thermal.collam[i];
			/* force to save wavelength convention as printout */
			if( awl > 100000 )
				awl /= 10000;
			fprintf( io, " Negative coolant was %s %.2f %.2e\n", 
			  thermal.chClntLab[i], awl, safe_div( thermal.heatnt[i], heat_total, 0. ) );
		}
	}

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

	/* this is mainly to keep the compiler from flagging possible paths 
	 * with j not being set */
	long i = INT_MIN;
	long j = INT_MIN;
	/* following loop tries to identify strongest agents and turn to labels */
	for( long k=0; k < ipnt; k++ )
	{
		/* generate labels that make sense in printout 
		 * if not identified with a specific agent, will print indices as [i][j],
		 * heating is thermal.heating(i,j) */
		i = ipsave[k];
		j = jpsave[k];
		/* i >= j indicates agent is one of the first LIMELM elements */
		if( i >= j )
		{
			if( dense.xIonDense[i][j] == 0. && thermal.heating(i,j)>SMALLFLOAT )
				fprintf(ioQQQ,"DISASTER assert about to be thrown - search for hit it\n");
			/*fprintf(ioQQQ,"DEBUG %li %li %.2e %.2e\n", i , j , 
				dense.xIonDense[i][j],
				thermal.heating(i,j));*/
			ASSERT( dense.xIonDense[i][j] > 0. || thermal.heating(i,j)<SMALLFLOAT );
			/* this is case of photoionization of atom or ion */
			chLabel[k] = elementnames.chElementSym[i];
			chLabel[k] += elementnames.chIonStage[j];
		}
		/* notice that in test i and j are swapped from order in heating array */
		else if( i == 0 && j == 1 )
		{
			/* photoionization from all excited states of Hydrogenic species,
			 * heating(0,1) */
			chLabel[k] = "Hn=2";
		}
		else if( i == 0 && j == 3 )
		{
			/* collisional ionization of all iso-seq from all levels, 
			 * heating(0,3) */
			chLabel[k] = "Hion";
		}
		else if( i == 0 && j == 7 )
		{
			/* UTA ionization heating(0,7) */
			chLabel[k] = " UTA";
		}
		else if( i == 0 && j == 8 )
		{
			/* thermal.heating(0,8) is heating due to collisions within 
			 * X of H2, code var hmi.HeatH2Dexc_used */
			chLabel[k] = "H2vH";
		}
		else if( i == 0 && j == 17 )
		{
			/* thermal.heating(0,17) is heating due to photodissociation 
			 * heating of X within H2,
			 * code var hmi.HeatH2Dish_used */
			chLabel[k] = "H2dH";
		}
		else if( i == 0 && j == 9 )
		{
			/* CO dissociation, co.CODissHeat, heating(0,9) */
			chLabel[k] = "COds";
		}
		else if( i == 0 && j == 20 )
		{
			/* extra heat thermal.heating(0,20)*/
			chLabel[k] = "extH";
		}
		else if( i == 0 && j == 21 )
		{
			/* pair heating thermal.heating(0,21)*/
			chLabel[k] = "pair";
		}
		else if( i == 0 && j == 11 )
		{
			/* free free heating, heating(0,11) */
			chLabel[k] = "H FF";
		}
		else if( i == 0 && j == 12 )
		{
			/* heating coolant (not line), physical cooling process, often a bug, heating(0,12) */
			chLabel[k] = "Hcol";
		}
		else if( i == 0 && j == 13 )
		{
			/* grain photoelectric effect, heating(0,13) */
			chLabel[k] = "GrnP";
		}
		else if( i == 0 && j == 14 )
		{
			/* grain collisions, heating(0,14) */
			chLabel[k] = "GrnC";
		}
		else if( i == 0 && j == 15 )
		{
			/* H- heating, heating(0,15) */
			chLabel[k] = "H-  ";
		}
		else if( i == 0 && j == 16 )
		{
			/* H2+ heating, heating(0,16) */
			chLabel[k] = "H2+ ";
		}
		else if( i == 0 && j == 18 )
		{
			/* H2 photoionization heating, heating(0,18) */
			chLabel[k] = "H2ph";
		}
		else if( i == 0 && j == 19 )
		{
			/* Compton heating, heating(0,19) */
			chLabel[k] = "Comp";
		}
		else if( i == 0 && j == 22 )
		{
			/* line heating, heating(0,22) */
			chLabel[k] = "line";
		}
		else if( i == 0 && j == 23 )
		{
			/* iso-sequence line heating - all elements together, 
			 * heating(0,23) */
			chLabel[k] = "Hlin";
		}
		else if( i == 0 && j == 24 )
		{
			/* charge transfer heating, heating(0,24) */
			chLabel[k] = "ChaT";
		}
		else if( i == 0 && j == 27 )
		{
			/* heating due all database species, heating(0,27) */
			chLabel[k] = "DbLn";
		}
		else if( i == 1 && j == 3 )
		{
			/* helium triplet line heating, heating(1,3) */
			chLabel[k] = "He3l";
		}
		else if( false && i == 1 && j == 5 )
		{
			/* advective heating, heating(1,5) */
			chLabel[k] = "adve";
		}
		else if( i == 1 && j == 6 )
		{
			/* cosmic ray heating thermal.heating(1,6)*/
			chLabel[k] = "CR H";
		}
		else if( i == 25 && j == 27 )
		{
			/* Fe 2 line heating, heating(25,27) */
			chLabel[k] = "Fe 2";
		}
		else
		{
			ostringstream oss;
			oss << "[" << i << "][" << j << "]";
			chLabel[k] = oss.str();
		}
	}

	/* now sort by decreasing importance */
	/*spsort netlib routine to sort array returning sorted indices */
	UNUSED int nFail;
	spsort(/* input array to be sorted */
		&SaveVal[0], 
		/* number of values in x */
		ipnt, 
		/* permutation output array */
		&ipOrdered[0], 
		/* flag saying what to do - 1 sorts into increasing order, not changing
		 * the original routine */
		-1, 
		/* error condition, should be 0 */
		&nFail);

	/*>>chng 06 jun 06, change start of save to give same info as cooling 
	 * as per comment by Yumihiko Tsuzuki */
	/* begin the print out with zone number, total heating and cooling */
	fprintf( io, "%.5e\t%.4e\t%.4e\t%.4e\t%.4e", 
		radius.depth_mid_zone, 
		phycon.te, 
		heat_total, 
		cool_total,
		dynamics.Heat() );

	for( long k=0; k < ipnt; k++ )
	{
		int ip = ipOrdered[k];
		i = ipsave[ip];
		j = jpsave[ip];
		ASSERT( i<LIMELM && j<LIMELM );
		if(k > 4 && safe_div( thermal.heating(i,j), heat_total, 0. ) < save.WeakHeatCool )
			break;
		fprintf( io, "\t%s\t%.7f ", 
			 chLabel[ip].c_str(), SaveVal[ip] );
	}
	fprintf( io, " \n" );

	/* a negative pointer in the heating array is probably a problem,
	 * indicating that some line has become a heat source */
	bool lgHeatLine = false;

	/* check if any lines were major heat sources */
	for( i=0; i < ipnt; i++ )
	{
		/* heating(22,0) is line heating - identify line if important */
		if( ipsave[ipOrdered[i]] == 0 && jpsave[ipOrdered[i]] == 22 )
			lgHeatLine = true;
	}

	if( lgHeatLine )
	{
		long level = -1;
		/* a line was a major heat source - identify it */
		TransitionProxy t = FndLineHt(&level);
		if( safe_div( t.Coll().heat(), heat_total, 0. ) > 0.005 )
		{
			ASSERT( t.associated() );
			double TauIn = t.Emis().TauIn();
			double Pump = t.Emis().pump();
			double EscP = t.Emis().Pesc();
			double CS = t.Coll().col_str();
			/* ratio of line to total heating */
			double ColHeat = safe_div( t.Coll().heat(), heat_total, 0. );
			
			fprintf( io, "  LHeat lv%2ld %s TIn%10.2e Pmp%9.1e EscP%9.1e CS%9.1e Hlin/tot%10.2e\n", 
			  level, chLineLbl(t).c_str(), TauIn, Pump, EscP, CS, ColHeat );
		}
	}
	return;
}
