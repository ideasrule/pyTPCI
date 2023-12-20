/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*abund_starburst generate abundance set from Fred Hamann's starburst evolution grid */
#include "cddefines.h"
#include "optimize.h"
#include "input.h"
#include "elementnames.h"
#include "abund.h"
#include "parser.h"

void abund_starburst(Parser &p)
{
	bool lgDebug;
	long int j;
	double sqrzed, 
	  zHigh, 
	  zal, 
	  zar, 
	  zc, 
	  zca, 
	  zcl, 
	  zco, 
	  zed, 
	  zed2, 
	  zed3, 
	  zed4, 
	  zedlog, 
	  zfe, 
	  zh, 
	  zhe, 
	  zmg, 
	  zn, 
	  zna, 
	  zne, 
	  zni, 
	  zo, 
	  zs, 
	  zsi;
	/* this is largest stored metallicity */
	static double zLimit = 35.5;

	DEBUG_ENTRY( "abund_starburst()" );

	// save default abundances since will use as scale factor
	realnum AbunReferSave[LIMELM];
	for( long i=1; i < LIMELM; i++ )
	{
		AbunReferSave[i] = abund.ReferenceAbun[i];
	}

	if( p.nMatch("TRAC") )
	{
		lgDebug = true;
		/* trace keyword did appear
		 * this will not be used, but must trick the sanity test */
		zHigh = zLimit;

		/* if trace option set, print header now, and init ZED */
		fprintf( ioQQQ, " Abundances relative to H, Z\n" );

		/* this is actual header line */
		fprintf( ioQQQ, "     Z   ");

		/* make line of chemical symbol names */
		for(j=0; j < 30; j++)
		{
			fprintf( ioQQQ, "    %2.2s", elementnames.chElementSym[j]);
		}
		fprintf( ioQQQ, "    \n" );

		zed = 1e-3;
	}
	else 
	{
		lgDebug = false;

		/* argument is metallicity  */
		zed = p.getNumberCheckLogLinNegImplLog("metallicity");

		if( zed <= 0. )
		{
			fprintf( ioQQQ, " Z .le.0 not allowed, Z=%10.2e\n", 
						zed );
			cdEXIT(EXIT_FAILURE);
		}

		zHigh = zed;
	}


	/* following if loop only if trace option is set
	 * >>chng 97 jun 17, get rid of go to
	 *999  continue
	 * */
	while( zed <= zHigh )
	{
		if( zed < 1e-3 || zed > zLimit )
		{
			fprintf( ioQQQ, " The metallicity must be between 1E-3 and%10.2e\n", 
			  zLimit );
			cdEXIT(EXIT_FAILURE);
		}
		zed2 = zed*zed;
		zed3 = zed2*zed;
		zed4 = zed3*zed;
		zedlog = log(zed);
		sqrzed = sqrt(zed);
		/* The value of each element as a function of ZED=[Z] */
		zh = 1.081646723 - 0.04600492*zed + 8.6569e-6*zed2 + 1.922e-5*
		  zed3 - 2.0087e-7*zed4;

		/* helium */
		zhe = 0.864675891 + 0.044423807*zed + 7.10886e-5*zed2 - 5.3242e-5*
		  zed3 + 5.70194e-7*zed4;
		abund.ReferenceAbun[1] = (realnum)zhe;

		/* li, b, bo unchanged */
		abund.ReferenceAbun[2] = 1.;
		abund.ReferenceAbun[3] = 1.;
		abund.ReferenceAbun[4] = 1.;

		/* carbon */
		zc = 0.347489799 + 0.016054107*zed - 0.00185855*zed2 + 5.43015e-5*
		  zed3 - 5.3123e-7*zed4;
		abund.ReferenceAbun[5] = (realnum)zc;

		/* nitrogen */
		zn = -0.06549567 + 0.332073984*zed + 0.009146066*zed2 - 0.00054441*
		  zed3 + 6.16363e-6*zed4;
		if( zn < 0.0 )
		{
			zn = 0.05193*zed;
		}
		if( zed < 0.3 )
		{
			zn = -0.00044731103 + 0.00026453554*zed + 0.52354843*zed2 - 
			  0.41156655*zed3 + 0.1290967*zed4;
			if( zn < 0.0 )
			{
				zn = 0.000344828*zed;
			}
		}
		abund.ReferenceAbun[6] = (realnum)zn;

		/* oxygen */
		zo = 1.450212747 - 0.05379041*zed + 0.000498919*zed2 + 1.09646e-5*
		  zed3 - 1.8147e-7*zed4;
		abund.ReferenceAbun[7] = (realnum)zo;

		/* neon */
		zne = 1.110139023 + 0.002551998*zed - 2.09516e-7*zed3 - 0.00798157*
		  POW2(zedlog);
		abund.ReferenceAbun[9] = (realnum)zne;

		/* fluorine, scale from neon */
		abund.ReferenceAbun[8] = abund.ReferenceAbun[9];

		/* sodium */
		zna = 1.072721387 - 0.02391599*POW2(zedlog) + .068644972*
		  zedlog + 0.017030935/sqrzed;
		/* this one is slightly negative at very low Z */
		zna = MAX2(1e-12,zna);
		abund.ReferenceAbun[10] = (realnum)zna;

		/* magnesium */
		zmg = 1.147209646 - 7.9491e-7*POW3(zed) - .00264458*POW2(zedlog) - 
		  0.00635552*zedlog;
		abund.ReferenceAbun[11] = (realnum)zmg;

		/* aluminium */
		zal = 1.068116358 - 0.00520227*sqrzed*zedlog - 0.01403851*
		  POW2(zedlog) + 0.066186787*zedlog;
		/* this one is slightly negative at very low Z */
		zal = MAX2(1e-12,zal);
		abund.ReferenceAbun[12] = (realnum)zal;

		/* silicon */
		zsi = 1.067372578 + 0.011818743*zed - 0.00107725*zed2 + 3.66056e-5*
		  zed3 - 3.556e-7*zed4;
		abund.ReferenceAbun[13] = (realnum)zsi;

		/* phosphorus scaled from silicon */
		abund.ReferenceAbun[14] = abund.ReferenceAbun[13];

		/* sulphur */
		zs = 1.12000;
		abund.ReferenceAbun[15] = (realnum)zs;

		/* chlorine */
		zcl = 1.10000;
		abund.ReferenceAbun[16] = (realnum)zcl;

		/* argon */
		zar = 1.091067724 + 2.51124e-6*zed3 - 0.0039589*sqrzed*zedlog + 
		  0.015686715*zedlog;
		abund.ReferenceAbun[17] = (realnum)zar;

		/* potassium scaled from silicon */
		abund.ReferenceAbun[18] = abund.ReferenceAbun[13];

		/* calcium */
		zca = 1.077553875 - 0.00888806*zed + 0.001479866*zed2 - 6.5689e-5*
		  zed3 + 1.16935e-6*zed4;
		abund.ReferenceAbun[19] = (realnum)zca;

		/* iron */
		zfe = 0.223713045 + 0.001400746*zed + 0.000624734*zed2 - 3.5629e-5*
		  zed3 + 8.13184e-7*zed4;
		abund.ReferenceAbun[25] = (realnum)zfe;

		/* scandium, scaled from iron */
		abund.ReferenceAbun[20] = abund.ReferenceAbun[25];

		/* titanium, scaled from iron */
		abund.ReferenceAbun[21] = abund.ReferenceAbun[25];

		/* vanadium scaled from iron */
		abund.ReferenceAbun[22] = abund.ReferenceAbun[25];

		/* chromium scaled from iron */
		abund.ReferenceAbun[23] = abund.ReferenceAbun[25];

		/* manganese scaled from iron */
		abund.ReferenceAbun[24] = abund.ReferenceAbun[25];

		/* cobalt */
		zco = zfe;
		abund.ReferenceAbun[26] = (realnum)zco;

		/* nickel */
		zni = zfe;
		abund.ReferenceAbun[27] = (realnum)zni;

		/* copper scaled from iron */
		abund.ReferenceAbun[28] = abund.ReferenceAbun[25];

		/* zinc  scaled from iron */
		abund.ReferenceAbun[29] = abund.ReferenceAbun[25];

		/* rescale to true abundances */
		abund.ReferenceAbun[0] = 1.;
		abund.ReferenceAbun[1] = (realnum)(abund.ReferenceAbun[1]*AbunReferSave[1]/
		  zh);

		for( long i=2; i < LIMELM; i++ )
		{
			abund.ReferenceAbun[i] = (realnum)(abund.ReferenceAbun[i]*AbunReferSave[i]*
			  zed/zh);
		}

		if( lgDebug )
		{
			fprintf( ioQQQ, "%10.2e", zed );
			for( long i=0; i < LIMELM; i++ )
			{
				fprintf( ioQQQ, "%6.2f", log10(abund.ReferenceAbun[i]) );
			}
			fprintf( ioQQQ, "\n" );

			if( fp_equal( zed, zLimit ) )
			{
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* this trick makes sure last one is upper limit */
		if( zed < zLimit )
		{
			zed = MIN2(zed*1.5,zLimit);
		}
		else
		{
			zed *= 1.5;
		}
	}

	/* vary option */
	if( optimize.lgVarOn )
	{
		/* this is number of parameters to feed onto the input line */
		optimize.nvarxt[optimize.nparm] = 1;
		strcpy( optimize.chVarFmt[optimize.nparm], "ABUNDANCES STARBURST %f LOG" );
		/* following are upper and lower limits to metallicity range */
		optimize.varang[optimize.nparm][0] = (realnum)log10(1e-3);
		optimize.varang[optimize.nparm][1] = (realnum)log10(zLimit);
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		/* log of metallicity will be argument */
		optimize.vparm[0][optimize.nparm] = (realnum)log10(zed);
		optimize.vincr[optimize.nparm] = 0.2f;
		++optimize.nparm;
	}
	return;
}
