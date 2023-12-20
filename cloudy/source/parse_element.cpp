/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseElement parse options on element command */
#include "cddefines.h"
#include "optimize.h"
#include "abund.h"
#include "dense.h"
#include "input.h"
#include "called.h"
#include "iso.h"
#include "mole.h"
#include "elementnames.h"
#include "parser.h"
#include "deuterium.h"


void ParseElement( Parser &p)
{
	/* this will be used to remember how many levels are active in any element we turn off,
	 * so that we can retain state if turned back on  */
	static bool lgFirst = true;
	static long levels[NISO][LIMELM];
	bool lgEnd;

	bool lgForceLog=false, lgForceLinear=false;

	long int nelem;
	double param=0.0;

	DEBUG_ENTRY( "ParseElement()" );

	/* zero out this array if first call */
	if( lgFirst )
	{
		lgFirst = false;
		for( long i=0; i<NISO; ++i )
		{
			for( nelem=0; nelem<LIMELM; ++nelem )
			{
				levels[i][nelem] = 0;
			}
		}
	}
	/* say that abundances have been changed */
	abund.lgAbnReference = false;

	/* find which element - will be used in remainder of routine to 
	 * adjust aspects of this element */
	nelem = p.GetElem( );

	bool lgElementSet = false;
	if( nelem >= ipHYDROGEN )
	{
		lgElementSet = true;
		/* nelem is now the atomic number of the element, and must be correct */
		ASSERT( nelem>=0 && nelem < LIMELM );
	}

	/* look for log or linear keywords */
	if( p.nMatch(" LOG") )
		lgForceLog = true;
	else if( p.nMatch("LINE") )
		lgForceLinear = true;

	if( p.nMatch("SCAL") )
	{
		param = p.FFmtRead();

		/* enter abundance as scale factor, relative to what is in now */
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " There must be a number on this line.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		else
		{
			if( !lgElementSet )
			{
				fprintf( ioQQQ, 
					" ParseElement did not find an element on the following line:\n" );
				p.PrintLine( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}
			/* interpret as log unless forced linear */
			if( (lgForceLog || param <= 0.) && !lgForceLinear )
			{
				/* option for log of scale factor */
				param = exp10(param);
			}
			abund.ScaleElement[nelem] = (realnum)param;
		}
	}

	else if( p.nMatch("ABUN") )
	{
		param = p.FFmtRead();

		/* log of actual abundance */
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " There must be a number on this line.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		else
		{
			if( !lgElementSet )
			{
				fprintf( ioQQQ, 
					" ParseElement did not find an element on the following line:\n" );
				p.PrintLine( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}
			/* interpret as log unless forced linear */
			if( !lgForceLinear )
			{
				/* option for log of scale factor */
				if( param > log10(BIGFLOAT) )
				{
					fprintf(ioQQQ,"log of abundance was entered, but it too large for this cpu.\n");
					cdEXIT(EXIT_FAILURE);
				}
				param = exp10(param);
			}
			abund.ReferenceAbun[nelem] = (realnum)param;

			if( abund.ReferenceAbun[nelem] > 1. && called.lgTalk )
			{
				fprintf( ioQQQ, 
					" Please check the abundance of this element.  It seems high to me.\n" );
			}
		}
	}

	else if( p.nMatch("ISOT") )
	{
		if( p.nMatch(" ALL") )
		{
			for( int i = ipHELIUM; i < LIMELM; i++ )
				mole_global.lgTreatIsotopes[i] = true;
		}
		else
		{
			if( !lgElementSet )
			{
				fprintf( ioQQQ, 
					" ParseElement did not find an element on the following line:\n" );
				p.PrintLine( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}

			if( !dense.lgElmtOn[nelem] )
			{
				fprintf(ioQQQ,"Sorry, you cannot set the isotope fractions of %s since it has been turned off.\n" ,
					elementnames.chElementName[nelem] );
				cdEXIT(EXIT_FAILURE);
			}

			/* option to specify ionization distribution for this element */
			int iso = 0;
			while( true )
			{
				int Aiso = (int)p.FFmtRead();
				if( p.lgEOL() ) break;

				realnum Fiso = (realnum)p.FFmtRead();
				if( p.lgEOL() )
					p.NoNumb( "isotope abundance" );


				int j = abund.IsoAbn[nelem].setAbn( Aiso, Fiso );
				if( j == -1 )
				{
					fprintf(ioQQQ, "Element Isotope: For %s, illegal isotope mass number: %d\n",
						elementnames.chElementName[nelem], Aiso);
					fprintf(ioQQQ, "Choose one of:\t");
					abund.IsoAbn[nelem].prtIso( ioQQQ );
					cdEXIT(EXIT_FAILURE);
				}

				iso++;
			}

			if( iso != abund.IsoAbn[nelem].getNiso() )
			{
				fprintf(ioQQQ, "Element Isotope: %s requires %d isotope pairs to be specified, but found %d\n",
					elementnames.chElementName[nelem], abund.IsoAbn[nelem].getNiso(), iso);
				fprintf(ioQQQ, "Enter fractionations for all isotopes, eg:\n");
				fprintf(ioQQQ, "element %s isotopes", elementnames.chElementName[nelem]);
				abund.IsoAbn[nelem].prtIsoPairs( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}

			if( abund.IsoAbn[nelem].isAnyIllegal() )
			{
				fprintf(ioQQQ, "Element Isotope: Non-positive isotope fractions are illegal.\n");
				fprintf(ioQQQ, "Read: %s\t iso:", elementnames.chElementName[nelem]);
				abund.IsoAbn[nelem].prtIsoPairs( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}

			abund.IsoAbn[nelem].normAbn( );
			//	abund.IsoAbn[nelem].prtIsoPairs( stdout );

			mole_global.lgTreatIsotopes[nelem] = true;
			if( nelem == ipHYDROGEN )
			{
				deut.lgElmtOn = true;
				if( abund.IsoAbn[nelem].getAbn(2) == 0.0 )
					deut.lgElmtOn = false;
			}
		}
	}

	else if( p.nMatch(" OFF") )
	{
		/* keyword LIMIT sets lower limit to abundance of element
		 * that will be included */
		/* option to turn off this element, set abundance to zero */
		if( p.nMatch( "LIMI") )
		{
			param = p.FFmtRead();
			if( p.lgEOL() )
			{
				fprintf( ioQQQ, " There must be an abundances on the ELEMENT OFF LIMIT command.\n" );
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			dense.AbundanceLimit = (realnum)exp10(param);
		}
		else
		{
			if( !lgElementSet )
			{
				fprintf( ioQQQ, 
					" ParseElement did not find an element on the following line:\n" );
				p.PrintLine ( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}
			/* flags saying element off, explicitly set by user */
			dense.lgElmtOn[nelem] = false;
			dense.lgElmtSetOff[nelem] = true;
			/* another flag that element is off */
			dense.IonHigh[nelem] = -1;
			dense.lgSetIoniz[nelem] = false;
			/* set no levels for all elements that are turned off */
			levels[ipH_LIKE][nelem] = iso_sp[ipH_LIKE][nelem].numLevels_max;
			levels[ipHE_LIKE][nelem] = iso_sp[ipHE_LIKE][nelem].numLevels_max;
			iso_sp[ipH_LIKE][nelem].numLevels_max = 0;
			iso_sp[ipHE_LIKE][nelem].numLevels_max = 0;

			if( nelem == ipHYDROGEN )
			{
				fprintf( ioQQQ, " It is not possible to turn hydrogen off.\n" );
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	/* specify an ionization distribution */
	else if( p.nMatch("IONI") )
	{
		bool lgLogSet = false;
		int ion;
		int ihi , low;


		if( !lgElementSet )
		{
			fprintf( ioQQQ, 
				" ParseElement did not find an element on the following line:\n" );
			p.PrintLine( ioQQQ );
			cdEXIT(EXIT_FAILURE);
		}
		if( !dense.lgElmtOn[nelem] )
		{
			fprintf(ioQQQ,"Sorry, you cannot set the ionization of %s since it has been turned off.\n" ,
				elementnames.chElementName[nelem] );
			cdEXIT(EXIT_FAILURE);
		}

		/* option to specify ionization distribution for this element */
		dense.lgSetIoniz[nelem] = true;
		ion = 0;
		while( ion<nelem+2 )
		{
			/* the ionization fractions that are set when above set true,
			 * gas phase abundance is this times total abundance
			 * Ionization fraction for [nelem][ion] */
			dense.SetIoniz[nelem][ion] = (realnum)p.FFmtRead();

			if( p.lgEOL() ) 
				break;

			/* all are log if any are negative */
			if( dense.SetIoniz[nelem][ion] < 0. )
				lgLogSet = true;
			++ion;
		}
		param = dense.SetIoniz[nelem][0];

		/* zero out ones not specified */
		for( long i=ion; i<nelem+2; ++i )
			dense.SetIoniz[nelem][i] = 0.;

		/* convert rest to linear if any were negative */
		if( lgLogSet )
		{
			for( long i=0; i<ion; ++i )
				dense.SetIoniz[nelem][i] = exp10(dense.SetIoniz[nelem][i]);
		}

		/* now check that no zero abundances exist between lowest and highest non-zero values */
		low = 0;
		while( dense.SetIoniz[nelem][low]==0 && low < nelem+1 )
			++low;
		ihi = nelem+1;
		while( dense.SetIoniz[nelem][ihi]==0 && ihi > low )
			--ihi;

		if( ihi==low && dense.SetIoniz[nelem][ihi]==0 )
		{
			fprintf(ioQQQ," element ionization command has all zero ionization fractions.  This is not possible.\n Sorry\n");
			cdEXIT(EXIT_FAILURE);
		}
		realnum AbundSum=0.;
		for(long i=low; i<=ihi; ++i )
		{
			if( dense.SetIoniz[nelem][i]==0 )
			{
				fprintf(ioQQQ," element abundance command has zero abundance between positive values.  This is not possible.\n Sorry\n");
				cdEXIT(EXIT_FAILURE);
			}
			AbundSum += dense.SetIoniz[nelem][i];
		}

		// renormalize so that sum is unity - needed for various sum rules
		for(long i=low; i<=ihi; ++i )
			dense.SetIoniz[nelem][i] /= AbundSum;

	}

	/* option to turn element on - some ini files turn off elements,
	 * this can turn them back on */
	else if( p.nMatch(" ON ") )
	{

		if( !lgElementSet )
		{
			fprintf( ioQQQ, 
				" ParseElement did not find an element on the following line:\n" );
			p.PrintLine( ioQQQ );
			cdEXIT(EXIT_FAILURE);
		}
		/* option to turn off this element, set abundance to zero */
		dense.lgElmtOn[nelem] = true;
		dense.lgElmtSetOff[nelem] = false;
		/* reset levels to default if they were ever turned off with element off command */
		if( levels[ipH_LIKE][nelem] )
		{
			iso_sp[ipH_LIKE][nelem].numLevels_max = levels[ipH_LIKE][nelem];
			iso_sp[ipHE_LIKE][nelem].numLevels_max = levels[ipHE_LIKE][nelem];
		}
	}

	else if( p.nMatch("TABL") )
	{

		if( !lgElementSet )
		{
			fprintf( ioQQQ, 
				" ParseElement did not find an element on the following line:\n" );
			p.PrintLine( ioQQQ );
			cdEXIT(EXIT_FAILURE);
		}
		/* read in table of position-dependent abundances for a particular element. */
		abund.lgAbunTabl[nelem] = true;

		/* general flag saying this option turned on */
		abund.lgAbTaON = true;

		/* does the table give depth or radius?  keyword DEPTH
		 * says it is depth, default is radius */
		if( p.nMatch("DEPT") )
			abund.lgAbTaDepth[nelem] = true;
		else
			abund.lgAbTaDepth[nelem] = false;

		/* make sure not trying to change hydrogen */
		if( nelem == ipHYDROGEN )
		{
			fprintf( ioQQQ, " cannot change abundance of hydrogen.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* read pair giving radius and abundance */
		p.getline();	  
		abund.AbTabRad[0][nelem] = (realnum)p.FFmtRead();
		abund.AbTabFac[0][nelem] = (realnum)p.FFmtRead();

		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " no pairs entered - cannot interpolate\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* number of points in the table */
		abund.nAbunTabl = 2;

		lgEnd = false;
		/* LIMTAB is limit to number of points we can store and is
		 * set to 500 in abundances */
		while( !lgEnd && abund.nAbunTabl < LIMTABD )
		{
			p.getline();
			lgEnd = p.m_lgEOF;
			if( !lgEnd )
			{
				lgEnd = p.hasCommand("END");
			}

			/* lgEnd may have been set within above if, if end line encountered*/
			if( !lgEnd )
			{
				abund.AbTabRad[abund.nAbunTabl-1][nelem] = 
					(realnum)p.FFmtRead();

				abund.AbTabFac[abund.nAbunTabl-1][nelem] = 
					(realnum)p.FFmtRead();
				abund.nAbunTabl += 1;
			}
		}

		abund.nAbunTabl -= 1;

		/* now check that radii are in increasing order */
		for( long i=1; i < abund.nAbunTabl; i++ )
		{
			/* the radius values are assumed to be strictly increasing */
			if( abund.AbTabRad[i][nelem] <= abund.AbTabRad[i-1][nelem] )
			{
				fprintf( ioQQQ, "ParseElement: TABLE ELEMENT TABLE radii "
					"must be in increasing order\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	else
	{
		fprintf( ioQQQ, "ParseElement: ELEMENT command - there must be "
			"a keyword on this line.\n" );
		fprintf( ioQQQ, " The keys I know about are TABLE, SCALE, _OFF, "
			"_ON_, IONIZATION, ABUNDANCE, and ISOTOPES.\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* vary option */
	if( optimize.lgVarOn )
	{
		ASSERT(optimize.nparm < LIMPAR);
		if( p.nMatch("SCAL") )
		{
			/* vary scale factor */
			sprintf( optimize.chVarFmt[optimize.nparm], "ELEMENT %4.4s SCALE %%f LOG", 
			  elementnames.chElementNameShort[nelem] );

		}
		else if( p.nMatch("ABUN") )
		{
			/* vary absolute abundance */
			sprintf( optimize.chVarFmt[optimize.nparm], "ELEMENT %4.4s ABUNDANCE %%f LOG", 
			  elementnames.chElementNameShort[nelem] );
		}
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		if( param > 0. )
			optimize.vparm[0][optimize.nparm] = (realnum)log10(param);
		else
			optimize.vparm[0][optimize.nparm] = (realnum)param;
		optimize.vincr[optimize.nparm] = 0.2f;
		optimize.nvarxt[optimize.nparm] = 1;
		++optimize.nparm;
	}
	return;
}
