/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtFinal create PrtFinal pages of printout, emission line intensities, etc */
/*StuffComment routine to stuff comments into the stack of comments, def in lines.h */
/*gett2o3 analyze computed [OIII] spectrum to get t^2 */
/*gett2 analyze computed structure to get structural t^2 */
#include "cddefines.h"
#include "iso.h"
#include "cddrive.h"
#include "dynamics.h"
#include "lines.h"
#include "taulines.h"
#include "warnings.h"
#include "phycon.h"
#include "dense.h"
#include "grainvar.h"
#include "h2.h"
#include "hmi.h"
#include "thermal.h"
#include "hydrogenic.h"
#include "rt.h"
#include "atmdat.h"
#include "timesc.h"
#include "opacity.h"
#include "struc.h"
#include "pressure.h"
#include "conv.h"
#include "geometry.h"
#include "called.h"
#include "iterations.h"
#include "version.h"
#include "colden.h"
#include "input.h"
#include "radius.h"
#include "peimbt.h"
#include "ipoint.h"
#include "mean.h"
#include "wind.h"
#include "prt.h"
#include "rfield.h"
#include "freebound.h"
#include "lines_service.h"

// helper routine to center a line in the output
inline void PrintCenterLine(FILE* io,               // file pointer
			    const string& chLine,   // string to print, should not end with '\n'
			    size_t LineLen)         // width of the line
{
	size_t StrLen = chLine.length();
	ASSERT( StrLen < LineLen );
	int pad = (LineLen-StrLen)/2;
	fprintf( io, "%*c%s\n", pad, ' ', chLine.c_str() );
}

/*gett2o3 analyze computed [OIII] spectrum to get t^2 */
STATIC void gett2o3(realnum *tsqr);

/*gett2 analyze computed structure to get structural t^2 */
STATIC void gett2(realnum *tsqr);

/* helper routine for printing averaged quantities */
inline void PrintRatio(double q1, double q2)
{
	double ratio = ( q2 > SMALLFLOAT ) ? q1/q2 : 0.;
	fprintf( ioQQQ, " " );
	fprintf( ioQQQ, PrintEfmt("%9.2e", ratio) );
	return;
}

STATIC void PrintSpectrum ()
{
	vector<realnum> sclsav, scaled;
	vector<long int> ipSortLines;
	vector<double> xLog_line_lumin;
	vector<short int> lgPrt;
	vector<long> Slines;

	long int
	  i,
	  ipEmType,
	  iprnt,
	  nline,
	  j,
	  ipNegIntensity[33], 
	  nNegIntenLines;

	double a,
	  /* N.B. 8 is used for following preset in loop */
	  d[8],
	  snorm;


	DEBUG_ENTRY( "PrintSpectrum()" );


	/*----------------------------------------------------------
	 *
	 * first set scaled lines */

	/* get space for scaled */
	scaled.resize(LineSave.nsum);

	/* get space for xLog_line_lumin */
	xLog_line_lumin.resize(LineSave.nsum);

	/* this is option to not print certain contributions */
	/* gjf 98 jun 10, get space for array lgPrt */
	lgPrt.resize(LineSave.nsum);

	/* get space for sclsav */
	sclsav.resize(LineSave.nsum );

	Slines.resize(LineSave.nsum);

	/* get space for array of indices for lines, for possible sorting */
	ipSortLines.resize(LineSave.nsum );

	ASSERT( LineSave.ipNormWavL >= 0 );

	/* option to also print usual first two sets of line arrays 
	 * but for two sets of cumulative arrays for time-dependent sims too */
	int nEmType = 2;
	if( prt.lgPrintLineCumulative && ! dynamics.isInitialRelaxIteration( iteration ) )
		nEmType = 4;

	for( ipEmType=0; ipEmType<nEmType; ++ipEmType )
	{
		if( ! prt.lgPrintBlockIntrinsic && (ipEmType == 0 || ipEmType == 2) )
			continue;
		if( ! prt.lgPrintBlockEmergent  && (ipEmType == 1 || ipEmType == 3) )
			continue;

		if( ipEmType > 1 && strcmp( rfield.chCumuType, "NONE" ) == 0 )
			continue;

		/* this is the intensity of the line spectrum will be normalized to */
		snorm = LineSave.lines[LineSave.ipNormWavL].SumLine(ipEmType);

		/* check that this line has positive intensity */
		if( ((snorm <= SMALLDOUBLE ) || (LineSave.ipNormWavL < 0)) || (LineSave.ipNormWavL > LineSave.nsum) )
		{
			string wl_str;
			sprt_wl( wl_str, LineSave.lines[LineSave.ipNormWavL].wavelength() );
			fprintf(ioQQQ,
				"\n\n"
				" >>PROBLEM Normalization line (\"%s\" %s) has small or zero intensity, its value was %.2e and its intensity was set to 1.\n"
				" >>Please consider using another normalization line (this is set with the NORMALIZE command).\n",
				LineSave.lines[LineSave.ipNormWavL].chALab(), wl_str.c_str(), snorm);
			fprintf( ioQQQ, " >>The relative intensities will be meaningless, and many lines may appear too faint.\n" );
			snorm = 1.;
		}
		for( i=0; i < LineSave.nsum; i++ )
		{
			/* Do only for emergent lines */
			if( ipEmType == 1 || ipEmType == 3 )
				LineSave.lines[i].checkEmergent( ipEmType );

			/* when normalization line is off-scale small (generally a model
			 * with mis-set parameters) the scaled intensity can be larger than
			 * a realnum - this is not actually a problem since the number will
			 * overflow the format and hence be unreadable */
			double scale = LineSave.lines[i].SumLine(ipEmType)/snorm*LineSave.ScaleNormLine;
			/* this will become a realnum, so limit dynamic range */
			scale = MIN2(BIGFLOAT , scale );
			scale = MAX2( -BIGFLOAT , scale );

			/* find logs of ALL line intensities/luminosities */
			scaled[i] = (realnum)scale;

			// the keyword volatile works around a bug in g++
			// see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=65425
			volatile double line_lumin = LineSave.lines[i].SumLine(ipEmType) * radius.Conv2PrtInten;
			xLog_line_lumin[i] = ( line_lumin > 0. ) ? log10(line_lumin) : 0.;
		}

		/* now find which lines to print, which to ignore because they are the wrong type */
		for( i=0; i < LineSave.nsum; i++ )
		{
			/* never print unit normalization check, at least in main list */
			if( LineSave.lines[i].isUnit() || LineSave.lines[i].isUnitD() )
				lgPrt[i] = false;
			else if( !prt.lgPrnColl && LineSave.lines[i].isCollisional()  )
				lgPrt[i] = false;
			else if( !prt.lgPrnPump && LineSave.lines[i].isPump() )
				lgPrt[i] = false;
			else if( !prt.lgPrnInwd && ( LineSave.lines[i].isInward() ||
							LineSave.lines[i].isInwardContinuum() ||
							LineSave.lines[i].isInwardTotal() )
						   )
				lgPrt[i] = false;
			else if( !prt.lgPrnHeat && LineSave.lines[i].isHeat() )
				lgPrt[i] = false;
			else
				lgPrt[i] = true;
		}

		/* do not print relatively faint lines unless requested */
		nNegIntenLines = 0;

		/* set ipNegIntensity to bomb to make sure set in following */
		for(i=0; i< 32; i++ )
		{
			ipNegIntensity[i] = LONG_MAX;
		}

		for(i=0;i<8;++i)
		{
			d[i] = -DBL_MAX;
		}

		/* create header for blocks of emission line intensities */
		const char chIntensityType[4][100]=
		{"     Intrinsic" , "      Emergent" , "Cumulative intrinsic" , "Cumulative emergent" };
		ASSERT( ipEmType==0 || ipEmType==1 || ipEmType==2 || ipEmType==3 );
		/* if true then printing in 4 columns of lines, this is offset to
		 * center the title */
		fprintf( ioQQQ, "\n" );
		if( prt.lgPrtLineArray )
			fprintf( ioQQQ, "                                              " );
		fprintf( ioQQQ, "%s" , chIntensityType[ipEmType] );
		fprintf( ioQQQ, " line intensities\n" );
		// caution about emergent intensities when outward optical
		// depths are not yet known
		if( ipEmType==1 && iteration==1 )
			fprintf(ioQQQ," Caution: emergent intensities are not reliable on the "
			"first iteration.\n");

		/* option to only print brighter lines */
		if( prt.lgFaintOn )
		{
			iprnt = 0;
			for( i=0; i < LineSave.nsum; i++ )
			{
				/* this insanity can happen when arrays overrun */
				ASSERT( iprnt <= i);
				if( scaled[i] >= prt.TooFaint && lgPrt[i] )
				{
					if( prt.lgPrtLineLog )
					{
						xLog_line_lumin[iprnt] = log10(LineSave.lines[i].SumLine(ipEmType) * radius.Conv2PrtInten);
					}
					else
					{
						xLog_line_lumin[iprnt] = LineSave.lines[i].SumLine(ipEmType) * radius.Conv2PrtInten;
					}
					sclsav[iprnt] = scaled[i];
					/* check that null is correct, string overruns have 
					 * been a problem in the past */
					Slines[iprnt] = i;
					++iprnt;
				}
				else if( -scaled[i] > prt.TooFaint && lgPrt[i] )
				{
					/* negative intensities occur if line absorbs continuum */
					ipNegIntensity[nNegIntenLines] = i;
					nNegIntenLines = MIN2(32,nNegIntenLines+1);
				}
				/* special labels to give id for blocks of lines 
				 * do not add these labels when sorting by wavelength since invalid */
				else if( LineSave.lines[i].isSeparator() &&!prt.lgSortLines )
				{
					xLog_line_lumin[iprnt] = 0.;
					sclsav[iprnt] = 0.;
					Slines[iprnt] = i;
					++iprnt;
				}
			}
		}

		else
		{
			/* print everything */
			iprnt = LineSave.nsum;
			for( i=0; i < LineSave.nsum; i++ )
			{
				if( LineSave.lines[i].isSeparator() )
				{
					xLog_line_lumin[i] = 0.;
					sclsav[i] = 0.;
				}
				else
				{
					sclsav[i] = scaled[i];
				}
				Slines[i] = i;
				if( scaled[i] < 0. )
				{
					ipNegIntensity[nNegIntenLines] = i;
					nNegIntenLines = MIN2(32,nNegIntenLines+1);
				}
			}
		}

		/* reorder lines according to wavelength for comparison with spectrum
		 * including writing out the results */
		if( prt.lgSortLines )
		{
			int lgFlag;
			if( prt.lgSortLineWavelength )
			{
				/* first check if wavelength range specified */
				if( prt.wlSort1 >-0.1 )
				{
					j = 0;
					/* skip over those lines not desired */
					for( i=0; i<iprnt; ++i )
					{
						realnum wavelength=LineSave.lines[Slines[i]].wavelength();
						if( wavelength>= prt.wlSort1 && wavelength<= prt.wlSort2 )
						{
							if( j!=i )
							{
								sclsav[j] = sclsav[i];
								/* >>chng 05 feb 03, add this, had been left out, 
								 * thanks to Marcus Copetti for discovering the bug */
								xLog_line_lumin[j] = xLog_line_lumin[i];
								Slines[j] = Slines[i];
							}
							++j;
						}
					}
					iprnt = j;
				}

				vector<realnum> wavelength(iprnt);
				for( i=0; i<iprnt; ++i )
				{
					wavelength[i] = LineSave.lines[Slines[i]].wavelength();
				}

				/*spsort netlib routine to sort array returning sorted indices */
				if( iprnt > 0 )
				{
					spsort(&wavelength[0], 
							 iprnt, 
							 &ipSortLines[0], 
							 /* flag saying what to do - 1 sorts into increasing order, not changing
							  * the original routine */
							 1, 
							 &lgFlag);
					if( lgFlag ) 
					{
						fprintf(ioQQQ," >>> PROBLEM spsort reports error\n");
						TotalInsanity();
					}
				}
			}
			else if( prt.lgSortLineIntensity )
			{
				if( iprnt > 0 )
				{
					/*spsort netlib routine to sort array returning sorted indices */
					spsort(&sclsav[0], 
							 iprnt, 
						  &ipSortLines[0], 
						  /* flag saying what to do - -1 sorts into decreasing order, not changing
						   * the original routine */
						  -1, 
						  &lgFlag);
					if( lgFlag ) 
						TotalInsanity();
				}
			}
			else
			{
				/* only to keep lint happen, or in case vars clobbered */
				TotalInsanity();
			}
		}
		else
		{
			/* do not sort lines (usual case) simply print in incoming order */
			for( i=0; i<iprnt; ++i )
			{
				ipSortLines[i] = i;
			}
		}

		/* the number of lines to print better be positive */
		if (iprnt == 0)
			fprintf(ioQQQ," >>>> No strong lines found <<<<\n");

		/* print out all lines which made it through the faint filter,
		 * there are iprnt lines to print 
		 * can print in either 4 column (the default ) or one long
		 * column of lines */
		if( prt.lgPrtLineArray )
		{
			/* four or three columns ? - depends on how many sig figs */
			if( LineSave.sig_figs >= 5 )
			{
				nline = (iprnt + 2)/3;
			}
			else
			{
				/* nline will be the number of horizontal lines - 
				* the 4 represents the 4 columns */
				nline = (iprnt + 3)/4;
			}
		}
		else
		{
			/* this option print a single column of emission lines */
			nline = iprnt;
		}

		/* now loop over the spectrum, printing lines */
		for( i=0; i < nline; i++ )
		{

			/* this skips over nline per step, to go to the next column in 
			 * the output */
			for( j=i; j<iprnt; j = j + nline)
			{ 
				/* index with possibly reordered set of lines */
				long ipLin = ipSortLines[j];
				/* this is the actual print statement for the emission line
				 * spectrum*/
				if( LineSave.lines[Slines[ipLin]].isSeparator() )
				{
					/* special labels */
					/*fprintf( ioQQQ, "1111222223333333444444444      " ); */
					/* this string was set in */
					fprintf( ioQQQ, "%s",LineSave.chHoldComments[
							 (int)LineSave.lines[Slines[ipLin]].wavelength()].c_str() ); 
				}
				else
				{
					/* the label for the line */
					LineSave.lines[Slines[ipLin]].prt(ioQQQ);

					/* print the log of the intensity/luminosity of the 
					 * lines, the usual case */
					fprintf( ioQQQ, " " );
					if( prt.lgPrtLineLog )
					{
						fprintf( ioQQQ, "%*.3f", prt_linecol.absint_len, xLog_line_lumin[ipLin] );
					}
					else
					{
						/* print linear quantity instead */
						fprintf( ioQQQ, " %*.2e", prt_linecol.absint_len, xLog_line_lumin[ipLin] );
					}

					/* print scaled intensity with various formats,
					 * depending on order of magnitude.  value is
					 * always the same but the format changes. */
					fprintf( ioQQQ, " " );
					if( sclsav[ipLin] < 9999.5 )
					{
						fprintf( ioQQQ, "%*.4f", prt_linecol.relint_len, sclsav[ipLin] );
					}
					else if( sclsav[ipLin] < 99999.5 )
					{
						fprintf( ioQQQ, "%*.3f", prt_linecol.relint_len, sclsav[ipLin] );
					}
					else if( sclsav[ipLin] < 999999.5 )
					{
						fprintf( ioQQQ, "%*.2f", prt_linecol.relint_len, sclsav[ipLin] );
					}
					else if( sclsav[ipLin] < 9999999.5 )
					{
						fprintf( ioQQQ, "%*.1f", prt_linecol.relint_len, sclsav[ipLin] );
					}
					else
					{
						fprintf( ioQQQ, "%s", prt_linecol.relint_outrange.c_str() );
					}
				}

				/* now end the block with some spaces to set next one off */
				if( j+nline < iprnt )
					fprintf( ioQQQ, "%s", prt_linecol.col_gap.c_str() );
			}
			/* now end the lines */
			fprintf( ioQQQ, "\n" );
		}

		if( nNegIntenLines > 0 )
		{
			fprintf( ioQQQ, " Lines with negative intensities -  Linear intensities relative to normalize line.\n" );
			fprintf( ioQQQ, "          " );

			for( i=0; i < nNegIntenLines; ++i )
			{
				fprintf( ioQQQ, "%ld %s %.1e, ", 
							ipNegIntensity[i], 
							LineSave.lines[ipNegIntensity[i]].label().c_str(), 
							scaled[ipNegIntensity[i]] );
			}
			fprintf( ioQQQ, "\n" );
		}
	}


	/* now find which were the very strongest, more that 5% of cooling */
	j = 0;
	for( i=0; i < LineSave.nsum; i++ )
	{
		a = (double)LineSave.lines[i].SumLine(0)/(double)thermal.totcol;
		if( (a >= 0.05) && LineSave.lines[i].chSumTyp() == 'c' )
		{
			ipNegIntensity[j] = i;
			d[j] = a;
			j = MIN2(j+1,7);
		}
	}

	fprintf( ioQQQ, "\n\n\n %s\n", input.chTitle.c_str() );
	if( j != 0 )
	{
		fprintf( ioQQQ, " Cooling: " );
		for( i=0; i < j; i++ )
		{
			fprintf( ioQQQ, " ");
			LineSave.lines[ipNegIntensity[i]].prt(ioQQQ);
			fprintf( ioQQQ, ":%5.3f", 
				d[i] );
		}
		fprintf( ioQQQ, "  \n" );
	}

	/* now find strongest heating, more that 5% of total */
	j = 0;
	for( i=0; i < LineSave.nsum; i++ )
	{
		a = safe_div((double)LineSave.lines[i].SumLine(0),(double)thermal.power,0.0);
		if( (a >= 0.05) && LineSave.lines[i].chSumTyp() == 'h' )
		{
			ipNegIntensity[j] = i;
			d[j] = a;
			j = MIN2(j+1,7);
		}
	}

	if( j != 0 )
	{
		fprintf( ioQQQ, " Heating: " );
		for( i=0; i < j; i++ )
		{
			fprintf( ioQQQ, " ");
			LineSave.lines[ipNegIntensity[i]].prt(ioQQQ);
			fprintf( ioQQQ, ":%5.3f", 
				d[i] );
		}
		fprintf( ioQQQ, "  \n" );
	}

	return;
}

void PrtFinal(void)
{
	char chCKey[5], 
	  chGeo[7], 
	  chPlaw[21];

	long int ip2500, 
	  ip2kev,
	  j;
	double bacobs, 
	  absint, 
	  bacthn, 
	  hbcab, 
	  hbeta;

	double a,
	  ajmass, 
	  ajmin, 
	  alfox, 
	  bot, 
	  bottom, 
	  bremtk, 
	  coleff, 
	  dmean, 
	  ferode, 
	  he4686, 
	  he5876, 
	  heabun, 
	  hescal, 
	  pion, 
	  plow, 
	  powerl, 
	  qion, 
	  qlow, 
	  rate, 
	  ratio, 
	  reclin, 
	  rjeans, 
	  tmean, 
	  top, 
	  THI,/* HI 21 cm spin temperature */
	  uhel, 
	  uhl, 
	  usp, 
	  wmean, 
	  znit;

	double bac, 
	  tel, 
	  x;

	DEBUG_ENTRY( "PrtFinal()" );

	/* return if not talking */
	if( !called.lgTalk )
	{ 
		return;
	}

	/* print out header, or parts that were saved */

	/* this would be a major logical error */
	ASSERT( LineSave.nsum > 1 );

	/* print name and version number */
	fprintf( ioQQQ, "\f\n");
	fprintf( ioQQQ, "%23c", ' ' );
	int len = t_version::Inst().chVersion.length();
	int repeat = (72-len)/2;
	for( long i=0; i < repeat; ++i )
		fprintf( ioQQQ, "*" );
	fprintf( ioQQQ, "> Cloudy %s <", t_version::Inst().chVersion.c_str() );
	for( long i=0; i < 72-repeat-len; ++i )
		fprintf( ioQQQ, "*" );
	fprintf( ioQQQ, "\n" );

	for( size_t i=0; i < input.crd.size(); i++ )
	{
		if( input.crd[i]->InclLevel > 0 || !input.crd[i]->lgVisible )
			continue;

		/* copy start of command to key, making it into capitals */
		cap4(chCKey,input.crd[i]->chCardSav.c_str());

		/* print if not continue */
		if( strcmp(chCKey,"CONT") != 0 )
			fprintf( ioQQQ, "%23c* %-80s*\n", ' ', input.crd[i]->chCardSav.c_str() );
	}

	/* print log of ionization parameter if greater than zero - U==0 for PDR calcs */
	if( rfield.uh > 0. )
	{
		a = log10(rfield.uh);
	}
	else
	{
		a = -37.;
	}

	fprintf( ioQQQ, 
		"                       *********************************> Log(U):%6.2f <*********************************\n", 
	  a );

	if( warnings.lgWarngs )
	{
		fprintf( ioQQQ, "  \n" );
		fprintf( ioQQQ, "                       >>>>>>>>>> Warning!\n" );
		fprintf( ioQQQ, "                       >>>>>>>>>> Warning!\n" );
		fprintf( ioQQQ, "                       >>>>>>>>>> Warning!  Warnings exist, this calculation has serious problems.\n" );
		fprintf( ioQQQ, "                       >>>>>>>>>> Warning!\n" );
		fprintf( ioQQQ, "                       >>>>>>>>>> Warning!\n" );
		fprintf( ioQQQ, "  \n" );
	}
	else if( warnings.lgCautns )
	{
		fprintf( ioQQQ, 
			"                       >>>>>>>>>> Cautions are present.\n" );
	}

	if( conv.lgBadStop )
	{
		fprintf( ioQQQ, 
			"                       >>>>>>>>>> The calculation stopped for unintended reasons.\n" );
	}

	if( iterations.lgIterAgain )
	{
		fprintf( ioQQQ, 
			"                       >>>>>>>>>> Another iteration is needed.  Use the ITERATE command.\n" );
	}

	/* open or closed geometry?? */
	if( geometry.lgSphere )
	{
		strcpy( chGeo, "Closed" );
	}
	else
	{
		strcpy( chGeo, "  Open" );
	}

	/* now give description of pressure law */
	if( strcmp(dense.chDenseLaw,"CPRE") == 0 )
	{
		strcpy( chPlaw, " Constant  Pressure " );
	}

	else if( strcmp(dense.chDenseLaw,"CDEN") == 0 )
	{
		strcpy( chPlaw, "  Constant  Density " );
	}

	else if( (strcmp(dense.chDenseLaw,"POWD") == 0 || strcmp(dense.chDenseLaw
	  ,"POWR") == 0) || strcmp(dense.chDenseLaw,"POWC") == 0 )
	{
		strcpy( chPlaw, "  Power Law Density " );
	}

	else if( strcmp(dense.chDenseLaw,"SINE") == 0 )
	{
		strcpy( chPlaw, " Rapid Fluctuations " );
	}

	else if( strncmp(dense.chDenseLaw , "DLW" , 3) == 0 )
	{
		strcpy( chPlaw, " Special Density Lw " );
	}

	else if( strcmp(dense.chDenseLaw,"HYDR") == 0 )
	{
		strcpy( chPlaw, " Hydrostatic Equlib " );
	}

	else if( strcmp(dense.chDenseLaw,"WIND") == 0 )
	{
		strcpy( chPlaw, "  Radia Driven Wind " );
	}

	else if( strcmp(dense.chDenseLaw,"DYNA") == 0 )
	{
		strcpy( chPlaw, "  Dynamical Flow    " );
	}

	else if( strcmp(dense.chDenseLaw,"GLOB") == 0 )
	{
		strcpy( chPlaw, " Globule            " );
	}

	else
	{
		strcpy( chPlaw, " UNRECOGNIZED CPRES " );
	}

	fprintf( ioQQQ, 
		"\n                     Emission Line Spectrum. %20.20sModel.  %6.6s geometry.  Iteration %ld of %ld.\n", 
	  chPlaw, chGeo, iteration, iterations.itermx + 1 );

	/* emission lines as logs of total luminosity */
	string chLine, chUnit;
	if(  radius.distance > 0. && radius.lgRadiusKnown && prt.lgPrintFluxEarth )
	{
		if( geometry.iEmissPower == 1 && !geometry.lgSizeSet )
			chUnit = "/arcsec";
		else if( geometry.iEmissPower == 0 && !geometry.lgSizeSet )
			chUnit = "/arcsec^2";
		else
			chUnit = "";

		chLine = "Flux observed at the Earth (erg/s/cm^2" + chUnit + ").";
		PrintCenterLine( ioQQQ, chLine, 130 );
	}
	else if(  prt.lgSurfaceBrightness )
	{
		if( prt.lgSurfaceBrightness_SR )
			chUnit = "sr";
		else
			chUnit = "arcsec^2";

		chLine = "Surface brightness (erg/s/cm^2/" + chUnit + ").";
		PrintCenterLine( ioQQQ, chLine, 130 );
	}
	else if( radius.lgPredLumin )
	{
		if( geometry.iEmissPower == 2 )
			chUnit = "erg/s";
		else if( geometry.iEmissPower == 1 )
			chUnit = "erg/s/cm";
		else if( geometry.iEmissPower == 0 )
			chUnit = "erg/s/cm^2";
		else
			TotalInsanity();

		ostringstream chCoverage;
		if( fp_equal( geometry.covgeo, 1_r ) )
			chCoverage << "with full coverage";
		else
		{
			chCoverage << "with a covering factor of " << fixed;
			chCoverage << setprecision(1) << geometry.covgeo*100. << "%";
		}

		if( radius.lgCylnOn )
			chLine = "Luminosity (" + chUnit + ") emitted by a partial cylinder ";
		else
			chLine = "Luminosity (" + chUnit + ") emitted by a shell ";
		chLine += chCoverage.str() + ".";
		PrintCenterLine( ioQQQ, chLine, 130 );

		if( geometry.iEmissPower != 2 )
		{
			string chAper;
			if( geometry.iEmissPower == 1 )
				chAper = "long slit";
			else if( geometry.iEmissPower == 0 )
				chAper = "pencil beam";
			else
				TotalInsanity();

			ostringstream chCovAper;
			chCovAper << fixed << setprecision(1) << geometry.covaper*100.;

			chLine = "Observed through a " + chAper + " with aperture covering factor ";
			chLine += chCovAper.str() + "%.";
			PrintCenterLine( ioQQQ, chLine, 130 );
		}
	}
	else
	{
		chLine = "Intensity (erg/s/cm^2).";
		PrintCenterLine( ioQQQ, chLine, 130 );
	}

	fprintf( ioQQQ, "\n" );

	/******************************************************************
	 * kill Hummer & Storey case b predictions if outside their table *
	 ******************************************************************/

	/* >>chng 02 aug 29, from lgHCaseBOK to following - caught by Ryan Porter */
	if( !atmdat.lgHCaseBOK[1][ipHYDROGEN] )
	{
		static const int NWLH = 21;
		/* list of all H case b lines */
		realnum wlh[NWLH] = { 6562.80e0f, Hbeta_WavLen, 4340.46e0f, 4101.73e0f, 1.87510e4f, 1.28180e4f,
						  1.09380e4f, 1.00493e4f, 2.62513e4f, 2.16551e4f, 1.94454e4f, 7.45777e4f,
						  4.65247e4f, 3.73951e4f, 4.05113e4f, 7.45777e4f, 3.29607e4f, 1215.67e0f,
						  1025.72e0f, 972.537e0f, 949.743e0f };

		/* table exceeded - kill all H case b predictions*/
		for( long i=0; i < LineSave.nsum; i++ )
		{
			/* >>chng 04 jun 21, kill both case a and b at same time,
			 * actually lgHCaseBOK has separate A and B flags, but
			 * this is simpler */
			if( LineSave.lines[i].isCaseB() || 
				 LineSave.lines[i].isCaseA() )
			{
				realnum errorwave;
				/* this logic must be kept parallel with which H lines are
				 * added as case B in lines_hydro - any automatic hosing of
				 * case b would kill both H I and He II */
				errorwave = WavlenErrorGet( LineSave.lines[i].wavelength(), LineSave.sig_figs );
				for( j=0; j<NWLH; ++j )
				{
					if( fabs(LineSave.lines[i].wavelength()-wlh[j] ) <= errorwave )
					{
						LineSave.lines[i].SumLineZero();
						break;
					}
				}
			}
		}
	}

	if( !atmdat.lgHCaseBOK[1][ipHELIUM] )
	{
		/* table exceeded - kill all He case b predictions*/
		static const int NWLHE = 20;
		realnum wlhe[NWLHE] = {1640.40e0f, 1215.13e0f, 1084.94e0f, 1025.27e0f, 4685.68e0f, 3203.04e0f,
						 2733.24e0f, 2511.15e0f, 1.01233e4f, 6559.91e0f, 5411.37e0f, 4859.18e0f,
						 1.86362e4f, 1.16260e4f, 9344.62, 8236.51e0f, 303.784e0f, 256.317e0f,
						 243.027e0f, 237.331e0f};
		for( long i=0; i < LineSave.nsum; i++ )
		{
			if( LineSave.lines[i].isCaseB() || 
				 LineSave.lines[i].isCaseA() )
			{
				realnum errorwave;
				/* this logic must be kept parallel with which H lines are
				 * added as case B in lines_hydro - any automatic hosing of
				 * case b would kill both H I and He II */
				errorwave = WavlenErrorGet( LineSave.lines[i].wavelength(), LineSave.sig_figs );
				for( j=0; j<NWLHE; ++j )
				{
					if( fabs(LineSave.lines[i].wavelength()-wlhe[j] ) <= errorwave )
					{
						LineSave.lines[i].SumLineZero();
						break;
					}
				}
			}
		}
	}

	/**********************************************************
	 *analyse line arrays for abundances, temp, etc           *
	 **********************************************************/

	/* find apparent helium abundance, will not find these if helium is turned off */

	if( cdLine("H  1",wlAirVac(Hbeta_WavLen),&hbeta,&absint)<=0 )
	{
		if( dense.lgElmtOn[ipHYDROGEN] )
		{
			/* this is a major logical error if hydrogen is turned on */
			fprintf( ioQQQ, " PrtFinal could not find H  1  %.2f with cdLine.\n", Hbeta_WavLen );
			cdEXIT(EXIT_FAILURE);
		}
		else
		{
			hbeta = 0;
			absint = -37.;
		}
	}

	if( dense.lgElmtOn[ipHELIUM] )
	{
		/* get HeI 5876 triplet */
		/* >>chng 06 may 15, changed this so that it works for up to six sig figs. */
		if( cdLine("Blnd",wlAirVac(5875.66),&he5876,&absint)<=0 )
		{
			/* 06 aug 28, from numLevels_max to _local. */
			if( iso_sp[ipHE_LIKE][ipHELIUM].numLevels_local >= 14 )
			{
				/* this is a major logical error if helium is turned on */
				fprintf( ioQQQ, " PrtFinal could not find He 1 5876 with cdLine.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* get HeII 4686 */
		/* >>chng 06 may 15, changed this so that it works for up to six sig figs. */
		if( cdLine("He 2",wlAirVac(4685.68),&he4686,&absint)<=0 )
		{
			/* 06 aug 28, from numLevels_max to _local. */
			if( iso_sp[ipH_LIKE][ipHELIUM].numLevels_local > 5 )
			{
				/* this is a major logical error if helium is turned on 
				 * and the model atom has enough levels */
				fprintf( ioQQQ, " PrtFinal could not find He 2 4686 with cdLine.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	else
	{
		he5876 = 0.;
		absint = 0.;
		he4686 = 0.;
	}

	if( hbeta > 0. )
	{
		heabun = (he4686*0.078 + he5876*0.739)/hbeta;
	}
	else
	{
		heabun = 0.;
	}

	if( dense.lgElmtOn[ipHELIUM] )
	{
		hescal = heabun/(dense.gas_phase[ipHELIUM]/dense.gas_phase[ipHYDROGEN]);
	}
	else
	{
		hescal = 0.;
	}

	/* get temperature from [OIII] spectrum, o may be turned off,
	 * but lines still dumped into stack */
	if( atmdat.lgdBaseSourceExists[ipOXYGEN][2] )
	{
		double o4363, o4959, o5007;

		if( cdLine("O  3",wlAirVac(4958.91),&o4959,&absint)<=0 )
		{
			if( dense.lgElmtOn[ipOXYGEN] )
			{
				/* this is a major logical error if oxygen is turned on */
				fprintf( ioQQQ, " PrtFinal could not find O  3 4959 with cdLine.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		if( cdLine("O  3",wlAirVac(5006.84),&o5007,&absint)<=0 )
		{
			if( dense.lgElmtOn[ipOXYGEN] )
			{
				fprintf( ioQQQ, " PrtFinal could not find O  3 5007 with cdLine.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		if( cdLine("O  3",wlAirVac(4363.21),&o4363,&absint)<=0 )
		{
			if( dense.lgElmtOn[ipOXYGEN] )
			{
				fprintf( ioQQQ, " PrtFinal could not find O  3 4363 with cdLine.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		// Energy ratio of O3 lines -- better to get this from transition structures
		realnum o3enro = (realnum)(4.56e-12/3.98e-12);
		/* first find low density limit OIII zone temp */
		if( o4363 != 0. && o4959 != 0. && o5007 != 0. )
		{
			/* following will assume coll excitation only, so correct
			 * for 4363's that cascade as 5007 */
			bot = o4959+o5007 - o4363/o3enro;

			if( bot > 0. )
			{
				ratio = o4363/bot;
			}
			else
			{
				ratio = 0.178;
			}

			if( ratio > 0.177 )
			{
				/* ratio was above infinite temperature limit */
				peimbt.toiiir = 1e10;
			}
			else
			{
				/* o iii 5007+4959, As 96 NIST */
				/*The cs of the transitions 3P0,1,2 to 1D2 are added together to give oiii_cs3P1D2 */
				/*the cs of the transition 1D2-1S0 is mentioned as oiii_cs1D21S0*/
				/*The cs of the transitions 3P0,1,2 to 1S0 are added together to give oiii_cs3P1S0*/
				realnum o3cs12 = 2.2347f;
				realnum o3cs13 = 0.29811f;
				double a31 = 0.2262;
				double a32 = 1.685;
				realnum o3ex23 = 32947.;
				realnum o3br32 = (realnum)(a32/(a31 + a32));

				/* data for following set in OIII cooling routine
				 * ratio of collision strengths, factor of 4/3 makes 5007+4959
				 * >>chng 96 sep 07, reset cs to values at 10^4K,
				 * had been last temp in model */
				/*>>chng 06 jul 25, update to recent data */
				ratio = ratio/1.3333/(o3cs13/o3cs12);
				/* ratio of energies and branching ratio for 4363 */
				ratio = ratio/o3enro/o3br32;
				peimbt.toiiir = (realnum)(-o3ex23/log(ratio));
			}
		}
		else
		{
			peimbt.toiiir = 0.;
		}
	}

	if( geometry.iEmissPower == 2 )
	{
		/* find temperature from Balmer continuum */
		if( cdLine("Bac",3646.,&bacobs,&absint)<=0 )
		{
			fprintf( ioQQQ, " PrtFinal could not find Bac 3646 with cdLine.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* we pulled hbeta out of stack with cdLine above */
		if( hbeta > 0. )
		{
			bac = bacobs/hbeta;
		}
		else
		{
			bac = 0.;
		}
	}
	else
	{
		bac = 0.;
	}

	if( bac > 0. )
	{
		/*----------------------------------------------------------*
		 ***** TableCurve c:\tcwin2\balcon.for Sep 6, 1994 11:13:38 AM
		 ***** log bal/Hbet
		 ***** X= log temp
		 ***** Y= 
		 ***** Eqn# 6503  y=a+b/x+c/x^2+d/x^3+e/x^4+f/x^5
		 ***** r2=0.9999987190883581
		 ***** r2adj=0.9999910336185065
		 ***** StdErr=0.001705886369042427
		 ***** Fval=312277.1895753243
		 ***** a= -4.82796940090671
		 ***** b= 33.08493347410885
		 ***** c= -56.08886262205931
		 ***** d= 52.44759454803217
		 ***** e= -25.07958990094209
		 ***** f= 4.815046760060006
		 *----------------------------------------------------------*
		 * bac is double precision!!!! */
		x = 1.e0/log10(bac);
		tel = -4.827969400906710 + x*(33.08493347410885 + x*(-56.08886262205931 + 
		  x*(52.44759454803217 + x*(-25.07958990094209 + x*(4.815046760060006)))));

		if( tel > 1. && tel < 5. )
		{
			peimbt.tbac = (realnum)exp10(tel);
		}
		else
		{
			peimbt.tbac = 0.;
		}
	}
	else
	{
		peimbt.tbac = 0.;
	}

	if( geometry.iEmissPower == 2 )
	{
		/* find temperature from optically thin Balmer continuum and case B H-beta */
		if( cdLine("thin",3646.,&bacthn,&absint)<=0 )
		{
			fprintf( ioQQQ, " PrtFinal could not find thin 3646 with cdLine.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* >>chng 06 may 15, changed this so that it works for up to six sig figs. */
		if( cdLine("Ca B",wlAirVac(Hbeta_WavLen),&hbcab,&absint)<=0 )
		{
			fprintf( ioQQQ, " PrtFinal could not find Ca B %.2f with cdLine.\n", Hbeta_WavLen );
			cdEXIT(EXIT_FAILURE);
		}

		if( hbcab > 0. )
		{
			bacthn /= hbcab;
		}
		else
		{
			bacthn = 0.;
		}
	}
	else
	{
		bacthn = 0.;
	}

	if( bacthn > 0. )
	{
		/*----------------------------------------------------------*
		 ***** TableCurve c:\tcwin2\balcon.for Sep 6, 1994 11:13:38 AM
		 ***** log bal/Hbet
		 ***** X= log temp
		 ***** Y= 
		 ***** Eqn# 6503  y=a+b/x+c/x^2+d/x^3+e/x^4+f/x^5
		 ***** r2=0.9999987190883581
		 ***** r2adj=0.9999910336185065
		 ***** StdErr=0.001705886369042427
		 ***** Fval=312277.1895753243
		 ***** a= -4.82796940090671
		 ***** b= 33.08493347410885
		 ***** c= -56.08886262205931
		 ***** d= 52.44759454803217
		 ***** e= -25.07958990094209
		 ***** f= 4.815046760060006
		 *----------------------------------------------------------*
		 * bac is double precision!!!! */
		x = 1.e0/log10(bacthn);
		tel = -4.827969400906710 + x*(33.08493347410885 + x*(-56.08886262205931 + 
		  x*(52.44759454803217 + x*(-25.07958990094209 + x*(4.815046760060006)))));

		if( tel > 1. && tel < 5. )
		{
			peimbt.tbcthn = (realnum)exp10(tel);
		}
		else
		{
			peimbt.tbcthn = 0.;
		}
	}
	else
	{
		peimbt.tbcthn = 0.;
	}

	/* we now have temps from OIII ratio and BAC ratio, now
	 * do Peimbert analysis, getting To and t^2 */

	peimbt.tohyox = (realnum)((8.5*peimbt.toiiir - 7.5*peimbt.tbcthn - 228200. + 
	  sqrt(POW2(8.5*peimbt.toiiir-7.5*peimbt.tbcthn-228200.)+9.128e5*
	  peimbt.tbcthn))/2.);

	if( peimbt.tohyox > 0. )
	{
		peimbt.t2hyox = (realnum)((peimbt.tohyox - peimbt.tbcthn)/(1.7*peimbt.tohyox));
	}
	else
	{
		peimbt.t2hyox = 0.;
	}

	if( prt.lgPrintBlock )
		PrintSpectrum();

	// don't print this text twice...
	if( !prt.lgPrtCitations )
	{
		fprintf( ioQQQ, "\n" );
		DatabasePrintReference();
	}

	/* print out ionization parameters and radiation making it through */
	qlow = 0.;
	plow = 0.;
	for( long i=0; i < (iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon - 1); i++ )
	{
		/* N.B. in following en1ryd prevents overflow */
		plow += (rfield.flux[0][i] + rfield.ConInterOut[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i])*
		  EN1RYD*rfield.anu(i);
		qlow += rfield.flux[0][i] + rfield.ConInterOut[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i];
	}

	qlow *= radius.r1r0sq;
	plow *= radius.r1r0sq;
	if( plow > 0. )
	{
		qlow = log10(qlow * radius.Conv2PrtInten);
		plow = log10(plow * radius.Conv2PrtInten);
	}
	else
	{
		qlow = 0.;
		plow = 0.;
	}

	qion = 0.;
	pion = 0.;
	for( long i=iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1; i < rfield.nflux; i++ )
	{
		/* N.B. in following en1ryd prevents overflow */
		pion += (rfield.flux[0][i] + rfield.ConInterOut[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i])*
		  EN1RYD*rfield.anu(i);
		qion += rfield.flux[0][i] + rfield.ConInterOut[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i];
	}

	qion *= radius.r1r0sq;
	pion *= radius.r1r0sq;

	if( pion > 0. )
	{
		qion = log10(qion * radius.Conv2PrtInten);
		pion = log10(pion * radius.Conv2PrtInten);
	}
	else
	{
		qion = 0.;
		pion = 0.;
	}

	/* derive ionization parameter for spherical geometry */
	if( rfield.qhtot > 0. )
	{
		if( rfield.lgUSphON )
		{
			/* RSTROM is stromgren radius */
			usp = rfield.qhtot/POW2(rfield.rstrom/radius.rinner)/
			  2.998e10/dense.gas_phase[ipHYDROGEN];
			usp = log10(usp);
		}
		else
		{
			/* no stromgren radius found, use outer radius */
			usp = rfield.qhtot/radius.r1r0sq/2.998e10/dense.gas_phase[ipHYDROGEN];
			usp = log10(usp);
		}
	}

	else
	{
		usp = -37.;
	}

	if( rfield.uh > 0. )
	{
		uhl = log10(rfield.uh);
	}
	else
	{
		uhl = -37.;
	}

	if( rfield.uheii > 0. )
	{
		uhel = log10(rfield.uheii);
	}
	else
	{
		uhel = -37.;
	}

	fprintf( ioQQQ, 
		"\n IONIZE PARMET:  U(1-)%8.4f  U(4-):%8.4f  U(sp):%6.2f  "
		"Q(ion):  %7.3f  L(ion)%7.3f    Q(low):%7.3f    L(low)%7.3f\n", 
	  uhl, uhel, usp, qion, pion, qlow, plow );

	a = safe_div(fabs((thermal.power-thermal.totcol)*100.),thermal.power,0.0);
	/* output power and total cooling; can be neg for shocks, collisional ioniz */
	if( thermal.power > 0. )
	{
		powerl = log10(thermal.power * radius.Conv2PrtInten);
	}
	else
	{
		powerl = 0.;
	}

	if( thermal.totcol > 0. )
	{
		thermal.totcol = log10(thermal.totcol * radius.Conv2PrtInten);
	}
	else
	{
		thermal.totcol = 0.;
	}

	if( thermal.FreeFreeTotHeat > 0. )
	{
		thermal.FreeFreeTotHeat = log10(thermal.FreeFreeTotHeat * radius.Conv2PrtInten);
	}
	else
	{
		thermal.FreeFreeTotHeat = 0.;
	}

	/* following is recombination line intensity */
	reclin = totlin('r');
	if( reclin > 0. )
	{
		reclin = log10(reclin * radius.Conv2PrtInten);
	}
	else
	{
		reclin = 0.;
	}

	fprintf( ioQQQ, 
		" ENERGY BUDGET:  Heat:%8.3f  Coolg:%8.3f  Error:%5.1f%%  Rec Lin:%8.3f  F-F  H%7.3f    P(rad/tot)max     ", 
	  powerl, 
	  thermal.totcol, 
	  a, 
	  reclin, 
	  thermal.FreeFreeTotHeat );
	PrintE82( ioQQQ, pressure.RadBetaMax );
	// resolving power in fine continuum for this run
	fprintf(ioQQQ,"    R(F Con):%.3e",
			SPEEDLIGHT / rfield.fine_opac_velocity_width );

	fprintf( ioQQQ, "\n");

	/* effective x-ray column density, from 0.5keV attenuation, no scat
	 * IPXRY is pointer to 73.5 Ryd */
	coleff = opac.TauAbsGeo[0][rt.ipxry-1] - rt.tauxry;
	coleff /= 2.14e-22;

	/* find t^2 from H II structure */
	gett2(&peimbt.t2hstr);

	/* find t^2 from OIII structure */
	gett2o3(&peimbt.t2o3str);

	fprintf( ioQQQ, "\n     Col(Heff):      ");
	PrintE93(ioQQQ, coleff);
	fprintf( ioQQQ, "  snd travl time  ");
	PrintE82(ioQQQ, timesc.sound);
	fprintf( ioQQQ, " sec  Te-low: ");
	PrintE82(ioQQQ, thermal.tlowst);
	fprintf( ioQQQ, "  Te-hi: ");
	PrintE82(ioQQQ, thermal.thist);

	/* this is the intensity of the UV continuum at the illuminated face, relative to the Habing value, as
	 * defined by Tielens & Hollenbach 1985 */
	fprintf( ioQQQ, "  G0TH85: ");
	PrintE82( ioQQQ, hmi.UV_Cont_rel2_Habing_TH85_face );
	/* this is the intensity of the UV continuum at the illuminated face, relative to the Habing value, as
	 * defined by Draine & Bertoldi 1985 */
	fprintf( ioQQQ, "  G0DB96:");
	PrintE82( ioQQQ, hmi.UV_Cont_rel2_Draine_DB96_face );

	fprintf( ioQQQ, "\n");

	fprintf( ioQQQ, "  Emiss Measure    n(e)n(p) dl       ");
	PrintE93(ioQQQ, colden.dlnenp);
	fprintf( ioQQQ, "  n(e)n(He+)dl         ");
	PrintE93(ioQQQ, colden.dlnenHep);
	fprintf( ioQQQ, "  En(e)n(He++) dl         ");
	PrintE93(ioQQQ, colden.dlnenHepp);
	fprintf( ioQQQ, "  ne nC+:");
	PrintE82(ioQQQ, colden.dlnenCp);
	fprintf( ioQQQ, "\n");

	/* following is wl where gas goes thick to bremsstrahlung, in cm */
	if( rfield.EnergyBremsThin > 0. )
	{
		bremtk = RYDLAM*1e-8/rfield.EnergyBremsThin;
	}
	else
	{
		bremtk = 1e30;
	}

	/* apparent helium abundance */
	fprintf( ioQQQ, " He/Ha:");
	PrintE82( ioQQQ, heabun);

	/* he/h relative to correct helium abundance */
	fprintf(ioQQQ, "  =%7.2f*true  Lthin:",hescal);

	/* wavelength were structure is optically thick to bremsstrahlung absorption */
	PrintE82( ioQQQ, bremtk);

	/* this is ratio of conv.nTotalIoniz, the number of times ConvBase 
	 * was called during the model, over the number of zones.
	 * This is a measure of the convergence of the model - 
	 * includes initial search so worse for short calculations.
	 * It is a measure of how hard the model was to converge */
	if( nzone > 0 )
	{
		/* >>chng 07 feb 23, subtract n calls to do initial solution
		 * so this is the number of calls needed to converge the zones.
		 * different is to discount careful approach to molecular solutions
		 * in one zone models */
		znit = (double)(conv.nTotalIoniz-conv.nTotalIoniz_start)/(double)(nzone);
	}
	else
	{
		znit = 0.;
	}
	/* >>chng 02 jan 09, format from 5.3f to 5.2f */
	fprintf(ioQQQ, "  itr/zn:%5.2f",znit);

	/* sort-of same thing for large H2 molecule - number is number of level pop solutions per zone */
	fprintf(ioQQQ, "  H2 itr/zn:%6.2f",h2.H2_itrzn());

	/* say whether we used stored opacities (T) or derived them from scratch (F) */
	fprintf(ioQQQ, "  File Opacity: F" );

	/* log of total mass in grams if spherical, or gm/cm2 if plane parallel */
	/* this is mass per unit inner area */
	double xmass;
	// if spherical change to total mass, if pp leave gm/cm2
	if( radius.pirsq > 0. )
	{
		chUnit = "(gm)";
		xmass = dense.xMassTotal * exp10( (double)radius.pirsq ) ;
	}
	else
	{
		chUnit = "(gm/cm^2)";
		xmass = dense.xMassTotal;
	}
	fprintf(ioQQQ,"  MassTot %.2e  %s", xmass, chUnit.c_str() );
	fprintf(ioQQQ,"\n");

	/* this block is a series of prints dealing with 21 cm quantities
	 * first this is the temperature derived from Lya - 21 cm optical depths
	 * get harmonic mean HI temperature, as needed for 21 cm spin temperature */
	if( cdTemp( "opti",&THI,"volume" ) )
	{
		fprintf(ioQQQ,"\n>>>>>>>>>>>>>>>>> PrtFinal, impossible problem getting 21cm opt.\n");
		TotalInsanity();
	}
	fprintf( ioQQQ, "   Temps(21 cm)   T(21cm/Ly a)  ");
	PrintE82(ioQQQ, THI );

	/* get harmonic mean HI gas kin temperature, as needed for 21 cm spin temperature 
	 * >>chng 06 jul 21, this was over volume but hazy said radius - change to radius,
	 * bug discovered by Eric Pellegrini & Jack Baldwin  */
	/*if( cdTemp( "21cm",0,&THI,"volume" ) )*/
	if( cdTemp( "21cm",&THI,"radius" ) )
	{
		fprintf(ioQQQ,"\n>>>>>>>>>>>>>>>>> PrtFinal, impossible problem getting 21cm temp.\n");
		TotalInsanity();
	}
	fprintf(ioQQQ, "        T(<nH/Tkin>)  ");
	PrintE82(ioQQQ, THI);

	/* get harmonic mean HI 21 21 spin temperature, as needed for 21 cm spin temperature 
	 * for this, always weighted by radius, volume would be ignored were it present */
	if( cdTemp( "spin",&THI,"radius" ) )
	{
		fprintf(ioQQQ,"\n>>>>>>>>>>>>>>>>> PrtFinal, impossible problem getting 21cm spin.\n");
		TotalInsanity();
	}
	fprintf(ioQQQ, "          T(<nH/Tspin>)    ");
	PrintE82(ioQQQ, THI);

	/* now convert this HI spin temperature into a brightness temperature */
	THI *= (1. - sexp( HFLines[0].Emis().TauIn() ) );
	fprintf( ioQQQ, "          TB21cm:");
	PrintE82(ioQQQ, THI);
	fprintf( ioQQQ, "\n");

	fprintf( ioQQQ, "                   N(H0/Tspin)  ");
	PrintE82(ioQQQ, colden.H0_ov_Tspin );
	fprintf( ioQQQ, "        N(OH/Tkin)    ");
	PrintE82(ioQQQ, colden.OH_ov_Tspin );

	/* mean B weighted wrt 21 cm absorption */
	bot = cdB21cm();
	fprintf( ioQQQ, "          B(21cm)          ");
	PrintE82(ioQQQ, bot );

	/* end prints for 21 cm */
	fprintf(ioQQQ, "\n");

	/* timescale (sec here) for photoerosion of Fe-Ni */
	rate = timesc.TimeErode*2e-26;
	if( rate > SMALLFLOAT )
	{
		ferode = 1./rate;
	}
	else
	{
		ferode = 0.;
	}

	/* mean acceleration */
	if( wind.acldr > 0. )
	{
		wind.AccelAver /= wind.acldr;
	}
	else
	{
		wind.AccelAver = 0.;
	}

	/* DMEAN is mean density (gm per cc); mean temp is weighted by mass density */
	/* >>chng 02 aug 21, from radius.depth_x_fillfac to integral of radius times fillfac */
	dmean = colden.TotMassColl/SDIV(radius.depth_x_fillfac);
	tmean = colden.tmas/SDIV(colden.TotMassColl);
	/* mean mass per particle */
	wmean = colden.wmas/SDIV(colden.TotMassColl);

	fprintf( ioQQQ, "   <a>:");
	PrintE82(ioQQQ , wind.AccelAver);
	fprintf( ioQQQ, "  erdeFe");
	PrintE71(ioQQQ , ferode);
	fprintf( ioQQQ, "  Tcompt");
	PrintE82(ioQQQ , timesc.tcmptn);
	fprintf( ioQQQ, "  Tthr");
	PrintE82(ioQQQ , timesc.time_therm_long);
	fprintf( ioQQQ, "  <Tden>: ");
	PrintE82(ioQQQ , tmean);
	fprintf( ioQQQ, "  <dens>:");
	PrintE82(ioQQQ , dmean);
	fprintf( ioQQQ, "  <MolWgt>");
	PrintE82(ioQQQ , wmean);
	fprintf(ioQQQ,"\n");

	/* log of Jeans length and mass - this is mean over model */
	if( tmean > 0. )
	{
		rjeans = (log10(JEANS)+log10(tmean) - log10(dense.wmole) - log10(dense.xMassDensity*
			geometry.FillFac))/2.;
	}
	else
	{
		/* tmean undefined, set rjeans to large value so 0 printed below */
		rjeans = 40.f;
	}

	if( rjeans < 36. )
	{
		rjeans = exp10(rjeans);
		/* AJMASS is Jeans mass in solar units */
		ajmass = 3.*log10(rjeans/2.) + log10(4.*PI/3.*dense.xMassDensity*
		  geometry.FillFac) - log10(SOLAR_MASS);
		if( ajmass < 37 )
		{
			ajmass = exp10(ajmass);
		}
		else
		{
			ajmass = 0.;
		}
	}
	else
	{
		rjeans = 0.;
		ajmass = 0.;
	}

	/* Jeans length and mass - this is smallest over model */
	ajmin = colden.ajmmin - log10(SOLAR_MASS);
	if( ajmin < 37 )
	{
		ajmin = exp10(ajmin);
	}
	else
	{
		ajmin = 0.;
	}

	/* estimate alpha (o-x) */
	if( rfield.anu(rfield.nflux-1) > 150. )
	{
		/* generate pointers to energies where alpha ox will be evaluated */
		ip2kev = ipoint(147.);
		ip2500 = ipoint(0.3648);

		/* now get fluxes there */
		bottom = rfield.flux[0][ip2500-1]*
		  rfield.anu(ip2500-1)/rfield.widflx(ip2500-1);

		top = rfield.flux[0][ip2kev-1]*
		  rfield.anu(ip2kev-1)/rfield.widflx(ip2kev-1);

		/* generate alpha ox if denominator is non-zero */
		if( bottom > 1e-30 && top > 1e-30 )
		{
			ratio = log10(top) - log10(bottom);
			if( ratio < 36. && ratio > -36. )
			{
				ratio = exp10(ratio);
			}
			else
			{
				ratio = 0.;
			}
		}

		else
		{
			ratio = 0.;
		}

		if( ratio > 0. )
		{
			// the separate variable freq_ratio is needed to work around a bug in icc 10.0
			double freq_ratio = rfield.anu(ip2kev-1)/rfield.anu(ip2500-1);
			alfox = log(ratio)/log(freq_ratio);
		}
		else
		{
			alfox = 0.;
		}
	}
	else
	{
		alfox = 0.;
	}

	if( colden.rjnmin < 37 )
	{
		colden.rjnmin = exp10(colden.rjnmin);
	}
	else
	{
		colden.rjnmin = FLT_MAX;
	}

	fprintf( ioQQQ, "     Mean Jeans  l(cm)");
	PrintE82(ioQQQ, rjeans);
	fprintf( ioQQQ, "  M(sun)");
	PrintE82(ioQQQ, ajmass);
	fprintf( ioQQQ, "  smallest:     len(cm):");
	PrintE82(ioQQQ, colden.rjnmin);
	fprintf( ioQQQ, "  M(sun):");
	PrintE82(ioQQQ, ajmin);
	fprintf( ioQQQ, "  a_ox tran:%6.2f\n" , alfox);

	fprintf( ioQQQ, " Rion:");
	double R_ion;
	if( rfield.lgUSphON )
		R_ion = rfield.rstrom;
	else
		R_ion = radius.Radius;
	PrintE93(ioQQQ, R_ion);
	fprintf( ioQQQ, "  Dist:");
	PrintE82(ioQQQ, radius.distance);
	fprintf( ioQQQ, "  Diam:");
	if( radius.distance > 0. )
		PrintE93(ioQQQ, 2.*R_ion*AS1RAD/radius.distance);
	else
		PrintE93(ioQQQ, 0.);
	fprintf( ioQQQ, "\n");

	if( prt.lgPrintTime )
	{
		/* print execution time by default */
		alfox = cdExecTime();
	}
	else
	{
		/* flag set false with no time command, so that different runs can
		 * compare exactly */
		alfox = 0.;
	}

	/* some details about the hydrogen and helium ions */
	fprintf( ioQQQ, 
		" Hatom level%3ld  HatomType:%4.4s  HInducImp %2c"
		"  He tot level:%3ld  He2 level:  %3ld  ExecTime%8.1f\n", 
		/* 06 aug 28, from numLevels_max to _local. */
	  iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_local, 
	  hydro.chHTopType, 
	  TorF(hydro.lgHInducImp), 
	  /* 06 aug 28, from numLevels_max to _local. */
	  iso_sp[ipHE_LIKE][ipHELIUM].numLevels_local,
	  /* 06 aug 28, from numLevels_max to _local. */
	  iso_sp[ipH_LIKE][ipHELIUM].numLevels_local , 
	  alfox );

	/* now give an indication of the convergence error budget */
	fprintf( ioQQQ, 
		" ConvrgError(%%)  <eden>%7.3f  MaxEden%7.3f  <H-C>%7.2f  Max(H-C)%8.2f  <Press>%8.3f  MaxPrs er%7.3f\n", 
		conv.AverEdenError/nzone*100. , 
		conv.BigEdenError*100. ,
		conv.AverHeatCoolError/nzone*100. , 
		conv.BigHeatCoolError*100. ,
		conv.AverPressError/nzone*100. , 
		conv.BigPressError*100. );

	fprintf(ioQQQ,
		"  Continuity(%%)  chng Te%6.1f  elec den%6.1f  n(H2)%7.1f  n(CO)    %7.1f\n",
		phycon.BigJumpTe*100. ,
		phycon.BigJumpne*100. ,
		phycon.BigJumpH2*100. ,
		phycon.BigJumpCO*100. );

	/* print out some average quantities */
	fprintf( ioQQQ, "\n                                                        Averaged Quantities\n" );
	fprintf( ioQQQ, "             Te      Te(Ne)   Te(NeNp)  Te(NeHe+)Te(NeHe2+) Te(NeO+)  Te(NeO2+)"
		 "  Te(H2)     N(H)     Ne(O2+)   Ne(Np)\n" );
	static const char* weight[3] = { "Radius", "Area", "Volume" };
	int dmax = geometry.lgGeoPP ? 1 : 3;
	for( int dim=0; dim < dmax; ++dim )
	{
		fprintf( ioQQQ, " %6s:", weight[dim] );
		// <Te>/<1>
		PrintRatio( mean.TempMean[dim][0], mean.TempMean[dim][1] );
		// <Te*ne>/<ne>
		PrintRatio( mean.TempEdenMean[dim][0], mean.TempEdenMean[dim][1] );
		// <Te*ne*nion>/<ne*nion>
		PrintRatio( mean.TempIonEdenMean[dim][ipHYDROGEN][1][0], mean.TempIonEdenMean[dim][ipHYDROGEN][1][1] );
		PrintRatio( mean.TempIonEdenMean[dim][ipHELIUM][1][0], mean.TempIonEdenMean[dim][ipHELIUM][1][1] );
		PrintRatio( mean.TempIonEdenMean[dim][ipHELIUM][2][0], mean.TempIonEdenMean[dim][ipHELIUM][2][1] );
		PrintRatio( mean.TempIonEdenMean[dim][ipOXYGEN][1][0], mean.TempIonEdenMean[dim][ipOXYGEN][1][1] );
		PrintRatio( mean.TempIonEdenMean[dim][ipOXYGEN][2][0], mean.TempIonEdenMean[dim][ipOXYGEN][2][1] );
		// <Te*nH2>/<nH2>
		PrintRatio( mean.TempIonMean[dim][ipHYDROGEN][2][0], mean.TempIonMean[dim][ipHYDROGEN][2][1] );
		// <nH>/<1>
		PrintRatio( mean.xIonMean[dim][ipHYDROGEN][0][1], mean.TempMean[dim][1] );
		// <ne*nion>/<nion>
		PrintRatio( mean.TempIonEdenMean[dim][ipOXYGEN][2][1], mean.TempIonMean[dim][ipOXYGEN][2][1] );
		PrintRatio( mean.TempIonEdenMean[dim][ipHYDROGEN][1][1], mean.TempIonMean[dim][ipHYDROGEN][1][1] );
		fprintf( ioQQQ, "\n" );
	}

	/* print out Peimbert analysis, tsqden default 1e7, changed
	 * with set tsqden command */
	if( dense.gas_phase[ipHYDROGEN] < peimbt.tsqden )
	{
		fprintf(  ioQQQ, " \n" );

		/* temperature from the [OIII] 5007/4363 ratio */
		fprintf(  ioQQQ, " Peimbert T(OIIIr)");
		PrintE82( ioQQQ , peimbt.toiiir);

		/* temperature from predicted Balmer jump relative to Hbeta */
		fprintf(  ioQQQ, " T(Bac)");
		PrintE82( ioQQQ , peimbt.tbac);

		/* temperature predicted from optically thin Balmer jump rel to Hbeta */
		fprintf(  ioQQQ, " T(Hth)");
		PrintE82( ioQQQ , peimbt.tbcthn);

		/* t^2 predicted from the structure, weighted by H */
		fprintf(  ioQQQ, " t2(Hstrc)");
		fprintf( ioQQQ,PrintEfmt("%9.2e", peimbt.t2hstr));

		/* temperature from both [OIII] and the Balmer jump rel to Hbeta */
		fprintf(  ioQQQ, " T(O3-BAC)");
		PrintE82( ioQQQ , peimbt.tohyox);

		/* t2 from both [OIII] and the Balmer jump rel to Hbeta */
		fprintf(  ioQQQ, " t2(O3-BC)");
		fprintf( ioQQQ,PrintEfmt("%9.2e", peimbt.t2hyox));

		/* structural t2 from the O+2 predicted radial dependence */
		fprintf(  ioQQQ, " t2(O3str)");
		fprintf( ioQQQ,PrintEfmt("%9.2e", peimbt.t2o3str));

		fprintf(  ioQQQ, "\n");

		if( gv.lgDustOn() )
		{
			fprintf( ioQQQ, " Be careful: grains exist.  This spectrum was not corrected for reddening before analysis.\n" );
		}
	}

	/* print average (over radius) grain dust parameters if lgDustOn() is true,
	 * average quantities are incremented in radius_increment, zeroed in IterStart */
	if( gv.lgDustOn() && gv.lgGrainPhysicsOn )
	{
		char chQHeat;
		double AV , AB;
		double total_dust2gas = 0.;

		fprintf( ioQQQ, "\n Average Grain Properties (over radius):\n" );

		for( size_t i0=0; i0 < gv.bin.size(); i0 += 10 ) 
		{
			if( i0 > 0 )
				fprintf( ioQQQ, "\n" );

			/* this is upper limit to how many grain species we will print across line */
			size_t i1 = min(i0+10,gv.bin.size());

			fprintf( ioQQQ, "       " );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				chQHeat = (char)(gv.bin[nd].lgEverQHeat ? '*' : ' ');
				fprintf( ioQQQ, "%-12.12s%c", gv.bin[nd].chDstLab, chQHeat );
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, "    nd:" );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				if( nd != i0 ) fprintf( ioQQQ,"   " );
				fprintf( ioQQQ, "%7ld   ", (unsigned long)nd );
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, " <Tgr>:" );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				if( nd != i0 ) fprintf( ioQQQ,"   " );
				fprintf( ioQQQ,PrintEfmt("%10.3e", gv.bin[nd].avdust/radius.depth_x_fillfac));
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, " <Vel>:" );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				if( nd != i0 ) fprintf( ioQQQ,"   " );
				fprintf( ioQQQ,PrintEfmt("%10.3e", gv.bin[nd].avdft/radius.depth_x_fillfac));
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, " <Pot>:" );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				if( nd != i0 ) fprintf( ioQQQ,"   " );
				fprintf( ioQQQ,PrintEfmt("%10.3e", gv.bin[nd].avdpot/radius.depth_x_fillfac));
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, " <D/G>:" );
			for( size_t nd=i0; nd < i1; nd++ )
			{
				if( nd != i0 ) fprintf( ioQQQ,"   " );
				fprintf( ioQQQ,PrintEfmt("%10.3e", gv.bin[nd].avDGRatio/radius.depth_x_fillfac));
				/* add up total dust to gas mass ratio */
				total_dust2gas += gv.bin[nd].avDGRatio/radius.depth_x_fillfac;
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf(ioQQQ," Dust to gas ratio (by mass):");
		fprintf(ioQQQ,PrintEfmt("%10.3e", total_dust2gas));

		/* total extinction (conv to mags) at V and B per hydrogen, this includes
		 * forward scattering as an extinction process, so is what would be measured
		 * for a star, but not for an extended source where forward scattering
		 * should be discounted */
		AV = rfield.extin_mag_V_point/SDIV(colden.colden[ipCOL_HTOT]);
		AB = rfield.extin_mag_B_point/SDIV(colden.colden[ipCOL_HTOT]);
		/* print A_V/N(Htot) for point and extended sources */
		fprintf(ioQQQ,", A(V)/N(H)(pnt):%.3e, (ext):%.3e", 
			AV,
			rfield.extin_mag_V_extended/SDIV(colden.colden[ipCOL_HTOT]) );

		/* ratio of total to selective extinction */
		fprintf(ioQQQ,", R:");

		/* gray grains have AB - AV == 0 */
		if( fabs(AB-AV)>SMALLFLOAT )
		{
			fprintf(ioQQQ,"%.3e", AV/(AB-AV) );
		}
		else
		{
			fprintf(ioQQQ,"%.3e", 0. );
		}
		fprintf(ioQQQ," AV(ext):%.3e (pnt):%.3e\n",
			rfield.extin_mag_V_extended, 
			rfield.extin_mag_V_point);
	}

	/* option to make short printout */
	if( !prt.lgPrtShort && called.lgTalk )
	{
		/* print log of optical depths, 
		 * calls prtmet if print line optical depths command entered */
		PrtAllTau();

		/* only print mean ionization and emergent continuum on last iteration */
		if( iterations.lgLastIt )
		{
			/* option to print column densities, set with print column density command */
			if( prt.lgPrintColumns )
				PrtColumns(ioQQQ);
			/* print mean ionization fractions for all elements */
			PrtMeanIon('i', false, ioQQQ);
			/* print mean ionization fractions for all elements with density weighting*/
			PrtMeanIon('i', true , ioQQQ);
			/* print mean temperature fractions for all elements */
			PrtMeanIon('t', false, ioQQQ);
			/* print mean temperature fractions for all elements with density weighting */
			PrtMeanIon('t', true , ioQQQ);
		}
	}

	/* print input title for model */
	fprintf( ioQQQ, "%s\n\n", input.chTitle.c_str() );
	fflush(ioQQQ);

	{
		enum { DEBUG_LOC = false };
		if( DEBUG_LOC )
		{
			test_cdTemp_molecules();
		}
	}

	return;
}

/* routine to stuff comments into the stack of comments,
 * return is index to this comment */
long int StuffComment( const char * chComment )
{
	DEBUG_ENTRY( "StuffComment()" );

	/* only do this one time per core load */
	if( LineSave.ipass == 0 )
	{
		if( LineSave.nComment >= NHOLDCOMMENTS )
		{
			fprintf( ioQQQ, " Too many comments have been entered into line array; increase the value of NHOLDCOMMENTS.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		LineSave.chHoldComments[LineSave.nComment] = chComment;
		size_t n = LineSave.chHoldComments[LineSave.nComment].length();
		LineSave.chHoldComments[LineSave.nComment] += string(prt_linecol.column_len-n, '.');
	}

	++LineSave.nComment;
	return( LineSave.nComment-1 );
}

/*gett2 analyze computed structure to get structural t^2 */
STATIC void gett2(realnum *tsqr)
{
	long int i;

	double tmean;
	double a, 
	  as, 
	  b;

	DEBUG_ENTRY( "gett2()" );

	/* get T, t^2 */
	a = 0.;
	b = 0.;

	ASSERT( nzone < struc.nzlim);
	// struc.volstr[] is radius.dVeffAper saved as a function of nzone
	// we need this version of radius.dVeff since we want to compare to
	// emission lines that react to the APERTURE command
	for( i=0; i < nzone; i++ )
	{
		as = (double)(struc.volstr[i])*(double)(struc.hiistr[i])*
		  (double)(struc.ednstr[i]);
		a += (double)(struc.testr[i])*as;
		/* B is used twice below */
		b += as;
	}

	if( b <= 0. )
	{
		*tsqr = 0.;
	}
	else
	{
		/* following is H+ weighted mean temp over vol */
		tmean = a/b;
		a = 0.;

		ASSERT( nzone < struc.nzlim );
		for( i=0; i < nzone; i++ )
		{
			as = (double)(struc.volstr[i])*(double)(struc.hiistr[i])*
			  struc.ednstr[i];
			a += (POW2((double)(struc.testr[i]-tmean)))*as;
		}
		*tsqr = (realnum)(a/(b*(POW2(tmean))));
	}

	return;
}

/*gett2o3 analyze computed [OIII] spectrum to get t^2 */
STATIC void gett2o3(realnum *tsqr)
{
	long int i;
	double tmean;
	double a, 
	  as, 
	  b;

	DEBUG_ENTRY( "gett2o3()" );

	/* get T, t^2 */	a = 0.;
	b = 0.;
	ASSERT(nzone<struc.nzlim);
	// struc.volstr[] is radius.dVeffAper saved as a function of nzone
	// we need this version of radius.dVeff since we want to compare to
	// emission lines that react to the APERTURE command
	for( i=0; i < nzone; i++ )
	{
		as = (double)(struc.volstr[i])*(double)(struc.o3str[i])*
		  (double)(struc.ednstr[i]);
		a += (double)(struc.testr[i])*as;

		/* B is used twice below */
		b += as;
	}

	if( b <= 0. )
	{
		*tsqr = 0.;
	}

	else
	{
		/* following is H+ weighted mean temp over vol */
		tmean = a/b;
		a = 0.;
		ASSERT( nzone < struc.nzlim );
		for( i=0; i < nzone; i++ )
		{
			as = (double)(struc.volstr[i])*(double)(struc.o3str[i])*
			  struc.ednstr[i];
			a += (POW2((double)(struc.testr[i]-tmean)))*as;
		}
		*tsqr = (realnum)(a/(b*(POW2(tmean))));
	}
	return;
}
