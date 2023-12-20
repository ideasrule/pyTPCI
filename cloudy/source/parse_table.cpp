/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseTable parse the table read command */
/*lines_table invoked by table lines command, check if we can find all lines in a given list */
#include "cddefines.h"
#include "cddrive.h"
#include "optimize.h"
#include "rfield.h"
#include "trace.h"
#include "lines.h"
#include "radius.h"
#include "input.h"
#include "stars.h"
#include "prt.h"
#include "parser.h"
#include "save.h"

/*ReadTable called by TABLE READ to read in continuum from PUNCH TRANSMITTED CONTINUUM */
STATIC void ReadTable(const string& fnam);

/* these will become the label and wl for a possible list of lines,
 * obtained when tables lines used */
static string chLINE_LIST;

/* Black's ISM continuum, with He hole filled in */
static const int NISM = 23;
static double tnuism[NISM], 
  fnuism[NISM];

/* z=2 background,
 * >>refer	continuum	background	Haardt, Francesco, & Madau, Piero, 1996, 
 * >>refercon	ApJ, 461, 20 */
static const int NHM96 = 14;
/* log energy in Ryd */
static const double tnuHM96[NHM96]={-8,-1.722735683,-0.351545683,-0.222905683,-0.133385683,
/* changeg these two energies to prevent degeneracy */
-0.127655683,-0.004575683,0.297544317,0.476753,0.476756,0.588704317,
0.661374317,1.500814317,2.245164317};
/*-0.127655683,-0.004575683,0.297544317,0.476754317,0.476754317,0.588704317,*/
/*log J in the units of (erg cm^{-2} s^{-1} Hz^{-1} sr^{-1})*/
static const double fnuHM96[NHM96]={-32.53342863,-19.9789,-20.4204,-20.4443,-20.5756,-20.7546,
-21.2796,-21.6256,-21.8404,-21.4823,-22.2102,-22.9263,-23.32,-24.2865};

/* Mathews and Ferland generic AGN continuum */
static const int NAGN = 8;
static Energy tnuagn[NAGN];
static realnum tslagn[NAGN];

/* table Draine ISM continuum */
static const int NDRAINE = 15;
static double tnudrn[NDRAINE] , tsldrn[NDRAINE];

/* routine that stores values for above vectors */
STATIC void ZeroContin(void);

/* this allows the low energy point of any built in array to be reset to the
 * current low energy point in the code - nothing need be done if this is reset
 * tnu is array of energies, [0] is first, and we want it to be lower
 * fluxlog is flux at tnu, and may or may not be log
 * lgLog says whether it is */
STATIC void resetBltin( Energy *tnu , realnum *fluxlog , bool lgLog )
{
	/* this will multiply low-energy bounds of code and go into element[0]
	 * ensures that energy range is fully covered */
	const double RESETFACTOR = 0.98;
	double power;
	/* this makes sure we are called after emm is defined */
	ASSERT( rfield.emm()  > 0. );

	if( lgLog )
	{
		/* continuum comes in as log of flux */
		/* this is current power-law slope of low-energy continuum */
		power = (fluxlog[1] - fluxlog[0] ) / log10( tnu[1].Ryd()/tnu[0].Ryd() );
		/* this will be new low energy bounds to this continuum */
		tnu[0] = rfield.emm()*RESETFACTOR;
		fluxlog[0] = fluxlog[1] + power * log10( tnu[0].Ryd()/tnu[1].Ryd() );
	}
	else
	{
		/* continuum comes in as linear flux */
		/* this is current power-law slope of low-energy continuum */
		power = log10( fluxlog[1]/fluxlog[0]) / log10( tnu[1].Ryd()/tnu[0].Ryd() );
		/* this will be new low energy bounds to this continuum */
		tnu[0] = rfield.emm()*RESETFACTOR;
		fluxlog[0] = log10(fluxlog[1]) + power * log10( tnu[0].Ryd()/tnu[1].Ryd() );
		/* flux is not really log, we want linear */
		fluxlog[0] = exp10(fluxlog[0]);
	}
	/*fprintf(ioQQQ," power is %f lgLog is %i\n", power, lgLog );*/
	return;
}

void ParseTable(Parser &p)
{
	string chFile;	/*file name for table read */

	IntMode imode = IM_ILLEGAL_MODE;
	bool lgHit, 
	  lgLogSet;
	static bool lgCalled=false;

	long int i, 
	  j, 
	  nstar;

	double alpha, 
	  brakmm, 
	  brakxr, 
	  ConBreak, 
	  fac, 
	  scale, 
	  slopir, 
	  slopxr;

	bool lgNoContinuum = false,
		lgQuoteFound;

	DEBUG_ENTRY( "ParseTable()" );

	/* if first call then set up values for table */
	if( !lgCalled )
	{
		ZeroContin();
		lgCalled = true;
	}

	if( rfield.nShape >= LIMSPC )
	{
		fprintf( ioQQQ, " %ld is too many spectra entered.  Increase LIMSPC\n Sorry.\n", 
		  rfield.nShape );
		cdEXIT(EXIT_FAILURE);
	}

	/* four commands, tables line, read, SED, and star, have quotes on the
	 * lines giving file names.  must get quotes first so that filename
	 * does not confuse parser */
	lgQuoteFound = true;
	if( p.GetQuote( chFile ) )
		lgQuoteFound = false;

	// Backwards compatibility with older versions of build in SEDs,
	// no SED and no file name in quotes, just a keyword 
	bool lgKeyword = false;
	if(!lgQuoteFound)
	{
		if(p.nMatch("AKN1"))
		{
			chFile = "akn120.sed";
			lgKeyword = true;
			lgQuoteFound = true;
		}
		else if(p.nMatch("CRAB"))
		{
			if(p.nMatch("DAVIDSON"))
			{
				chFile = "CrabDavidson.sed";
			}
			else
			{
				chFile = "CrabHester.sed";
			}
			lgKeyword = true;
			lgQuoteFound = true;
		}
		else if(p.nMatch("COOL"))
		{
			chFile = "cool.sed";
			lgKeyword = true;
			lgQuoteFound = true;
		}
		else if(p.nMatch("RUBI"))
		{
			chFile = "Rubin.sed";
			lgKeyword = true;
			lgQuoteFound = true;
		}
		else if(p.nMatch("TRAP"))
		{
			chFile = "Trapezium.sed";
			lgKeyword = true;
			lgQuoteFound = true;
		}
		else if(p.nMatch(" XDR"))
		{
			chFile = "XDR.sed";
			lgKeyword = true;
			lgQuoteFound = true;
		}
	}

	// we reserve space for 1000 table points, no worries if that is too small though...
	ASSERT( rfield.tNu[rfield.nShape].empty() );
	rfield.tNu[rfield.nShape].reserve( 1000 );
	ASSERT( rfield.tslop[rfield.nShape].empty() );
	rfield.tslop[rfield.nShape].reserve( 1000 );
	ASSERT( rfield.tFluxLog[rfield.nShape].empty() );
	rfield.tFluxLog[rfield.nShape].reserve( 1000 );

	/* set flag telling interpolate */
	strcpy( rfield.chSpType[rfield.nShape], "INTER" );

	bool lgHM05 = false, lgHM12 = false;

	/* NB when adding more keys also change the comment at the end */
	if( p.nMatch(" AGN") )
	{
		/* do Mathews and Ferland generic AGN continuum */
		for( i=0; i < NAGN; i++ )
		{
			rfield.tNu[rfield.nShape].emplace_back( tnuagn[i] );
			rfield.tslop[rfield.nShape].emplace_back( tslagn[i] );
		}
		rfield.ncont[rfield.nShape] = NAGN;

		/* optional keyword break, to adjust IR cutoff */
		if( p.nMatch("BREA") )
		{
			ConBreak = p.FFmtRead();

			if( p.lgEOL() )
			{
				/* no break, set to low energy limit of code */
				if( p.nMatch(" NO ") )
				{
					ConBreak = rfield.emm()*1.01f;
				}
				else
				{
					fprintf( ioQQQ, " There must be a number for the break.\n Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}

			if( ConBreak == 0. )
			{
				fprintf( ioQQQ, " The break must be greater than 0.2 Ryd.\n Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			if( p.nMatch("MICR") )
			{
				/*  optional keyword, ``microns'', convert to Rydbergs */
				ConBreak = 0.0912/ConBreak;
			}

			if( ConBreak < 0. )
			{
				/*  option to enter break as LOG10 */
				ConBreak = exp10(ConBreak);
			}

			else if( ConBreak == 0. )
			{
				fprintf( ioQQQ, " An energy of 0 is not allowed.\n Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			if( ConBreak >= rfield.tNu[rfield.nShape][2].Ryd() )
			{
				fprintf( ioQQQ, " The energy of the break cannot be greater than%10.2e Ryd.\n Sorry.\n", 
				  rfield.tNu[rfield.nShape][2].Ryd() );
				cdEXIT(EXIT_FAILURE);
			}

			else if( ConBreak <= rfield.tNu[rfield.nShape][0].Ryd() )
			{
				fprintf( ioQQQ, " The energy of the break cannot be less than%10.2e Ryd.\n Sorry.\n", 
				  rfield.tNu[rfield.nShape][0].Ryd() );
				cdEXIT(EXIT_FAILURE);
			}

			rfield.tNu[rfield.nShape][1].set(ConBreak);

			rfield.tslop[rfield.nShape][1] = 
				(realnum)(rfield.tslop[rfield.nShape][2] + 
			  log10(rfield.tNu[rfield.nShape][2].Ryd()/rfield.tNu[rfield.nShape][1].Ryd()));

			rfield.tslop[rfield.nShape][0] = 
				(realnum)(rfield.tslop[rfield.nShape][1] - 
			  2.5*log10(rfield.tNu[rfield.nShape][1].Ryd()/rfield.tNu[rfield.nShape][0].Ryd()));
		}
	}

	else if( p.nMatchErase("HM96") )
	{
		/* this is the old Haardt & Madau continuum, one set of points
		 * with only the quasars
		 * this command does not include the CMB - do that separately with the CMB command */
		/* set flag telling interpolate */
		strcpy( rfield.chSpType[rfield.nShape], "INTER" );

		/* z=2 background,
		* >>refer	continuum	background	Haardt, Francesco, & Madau, Piero, 1996, ApJ, 461, 20 */
		for( j=0; j < NHM96; j++ )
		{
			/* frequency was stored as log of ryd */
			rfield.tNu[rfield.nShape].emplace_back( Energy(exp10( tnuHM96[j])) );
			rfield.tslop[rfield.nShape].emplace_back( (realnum)fnuHM96[j] );
		}
		rfield.ncont[rfield.nShape] = NHM96;

		/* optional scale factor to change default intensity from their value
		 * assumed to be log if negative, and linear otherwise */
		scale = p.FFmtRead();
		if( scale > 0. )
			scale = log10(scale);

		/* this also sets continuum intensity*/
		if( p.m_nqh >= LIMSPC )
		{
			fprintf( ioQQQ, " %ld is too many continua entered. Increase LIMSPC\n Sorry.\n", 
			p.m_nqh );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that stack of shape and luminosity specifications
		 * is parallel, stop if not - this happens is background comes
		 * BETWEEN another set of shape and luminosity commands */
		if( rfield.nShape != p.m_nqh )
		{
			fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );

		/* this is an isotropic radiation field */
		rfield.lgBeamed[p.m_nqh] = false;
		rfield.Illumination[p.m_nqh] = Illuminate::ISOTROPIC;

		/* this will be flux density at some frequency on the table.  the numbers
		 * are per Hz and sr so must multiply by 4 pi
		 * [2] is not special, could have been any within array*/
		rfield.range[p.m_nqh][0] = exp10( tnuHM96[2] )*1.0001;

		/* convert intensity HM96 give to current units of code */
		rfield.totpow[p.m_nqh] = (fnuHM96[2] + log10(PI4) + scale);

		++p.m_nqh;
	}

	else if( ( lgHM05 = p.nMatchErase("HM05") ) || ( lgHM12 = p.nMatchErase("HM12") ) )
	{
		bool lgQuasar = p.nMatch("QUAS");
		if( lgHM12 && lgQuasar )
		{
			fprintf( ioQQQ, " The QUASAR option is not supported on the TABLE HM12 command.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		// read requested redshift
		double redshift = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("redshift");
		}
		// read optional scale factor, defaults to 1.
		double scale = p.FFmtRead();
		if( scale > 0. )
			scale = log10(scale);
		strcpy( rfield.chSpType[rfield.nShape], "VOLK " );
		UNUSED double zlow, zhigh;
		int version = lgHM05 ? 2005 : 2012;
		// this does the hard work of interpolating the SEDs...
		rfield.ncont[rfield.nShape] = HaardtMadauInterpolate( redshift, version, lgQuasar, &zlow, &zhigh );

		// choose some frequency for normalizing the spectrum, precise value is not important
		long ip = rfield.ipointC(0.95);
		rfield.range[p.m_nqh][0] = rfield.tNu[rfield.nShape][ip].Ryd();

		// the numbers given by Haardt&Madau are per Hz and sr so must multiply by 4 pi
		rfield.totpow[p.m_nqh] = log10(rfield.tslop[rfield.nShape][ip]) + log10(PI4) + scale;

		/* this also sets continuum intensity*/
		if( p.m_nqh >= LIMSPC )
		{
			fprintf( ioQQQ, " %ld is too many continua entered. Increase LIMSPC\n Sorry.\n", p.m_nqh );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that stack of shape and luminosity specifications
		 * is parallel, stop if not - this happens is background comes
		 * BETWEEN another set of shape and luminosity commands */
		if( rfield.nShape != p.m_nqh )
		{
			fprintf( ioQQQ, " This command has come between a previous ordered pair of"
				 " continuum shape and luminosity commands.\n Reorder the commands"
				 " to complete each continuum specification before starting another.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );

		/* this is an isotropic radiation field */
		rfield.lgBeamed[p.m_nqh] = false;
		rfield.Illumination[p.m_nqh] = Illuminate::ISOTROPIC;

		++p.m_nqh;
	}
	else if( p.nMatchErase("KS18") )
	{
		fprintf(ioQQQ, "The TABLE KS18 command has been replaced by the TABLE KS19 command.\n");
		fprintf(ioQQQ, "Please update your script.\n");
		cdEXIT(EXIT_FAILURE);
	}
	else if( p.nMatchErase("KS19") )
	{
		// read requested redshift
		double redshift = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("redshift");
		}
		// read optional scale factor, defaults to 1.
		double scale = p.FFmtRead();
		if( scale > 0. )
			scale = log10(scale);
		// read optional Q value
		int Q = (int)p.FFmtRead();
		if( p.lgEOL() )
			Q = 18;
		if( Q < 14 || Q > 20 )
		{
			fprintf( ioQQQ, " Invalid value Q=%d, should be 14 <= Q <= 20.\n", Q );
			cdEXIT(EXIT_FAILURE);
		}
		strcpy( rfield.chSpType[rfield.nShape], "VOLK " );
		UNUSED double zlow, zhigh;
		// this does the hard work of interpolating the SEDs...
		rfield.ncont[rfield.nShape] = KhaireSrianandInterpolate( redshift, Q, &zlow, &zhigh );

		// choose some frequency for normalizing the spectrum, precise value is not important
		long ip = rfield.ipointC(0.95);
		rfield.range[p.m_nqh][0] = rfield.tNu[rfield.nShape][ip].Ryd();

		// the numbers given by Khaire&Srianand are per Hz and sr so must multiply by 4 pi
		rfield.totpow[p.m_nqh] = log10(rfield.tslop[rfield.nShape][ip]) + log10(PI4) + scale;

		/* this also sets continuum intensity*/
		if( p.m_nqh >= LIMSPC )
		{
			fprintf( ioQQQ, " %ld is too many continua entered. Increase LIMSPC\n Sorry.\n", p.m_nqh );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that stack of shape and luminosity specifications
		 * is parallel, stop if not - this happens is background comes
		 * BETWEEN another set of shape and luminosity commands */
		if( rfield.nShape != p.m_nqh )
		{
			fprintf( ioQQQ, " This command has come between a previous ordered pair of"
				 " continuum shape and luminosity commands.\n Reorder the commands"
				 " to complete each continuum specification before starting another.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );

		/* this is an isotropic radiation field */
		rfield.lgBeamed[p.m_nqh] = false;
		rfield.Illumination[p.m_nqh] = Illuminate::ISOTROPIC;

		++p.m_nqh;
	}
	else if( p.nMatch(" ISM") )
	{
		/* local ISM radiation field from Black 1987, Interstellar Processes */
		/* >>chng 04 mar 16, rm CMB from field so that it can be used at
		 * any redshift */
		rfield.tNu[rfield.nShape].emplace_back( Energy( exp10(6.)/FR1RYD) );
		rfield.tslop[rfield.nShape].emplace_back( (-21.21f - 6.f) );
		for( i=6; i < NISM; i++ )
		{
			/* energies were stored as log Hz and intensity as log nu Fnu, as per John's plot.
			 * convert to Ryd and log Fnu */
			rfield.tNu[rfield.nShape].emplace_back( Energy( exp10(tnuism[i])/FR1RYD ) );
			rfield.tslop[rfield.nShape].emplace_back( (realnum)(fnuism[i] - tnuism[i]) );
		}
		rfield.ncont[rfield.nShape] = NISM -5;

		/* optional scale factor to change default luminosity
		 * from observed value
		 * want final number to be log 
		 * assumed to be log if negative, and linear otherwise unless log option is present */
		scale = p.FFmtRead();
		if( scale > 0. && !p.nMatch(" LOG"))
			scale = log10(scale);

		/* this also sets continuum intensity*/
		if( p.m_nqh >= LIMSPC )
		{
			fprintf( ioQQQ, " %4ld is too many continua entered. Increase LIMSPC\n Sorry.\n", 
			  p.m_nqh );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that stack of shape and luminosity specifications
		 * is parallel, stop if not - this happens is background comes
		 * BETWEEN another set of shape and luminosity commands */
		if( rfield.nShape != p.m_nqh )
		{
			fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );

		/* this is an isotropic radiation field */
		rfield.lgBeamed[p.m_nqh] = false;
		rfield.Illumination[p.m_nqh] = Illuminate::ISOTROPIC;

		/* this will be flux density at 1 Ryd
		 * >>chng 96 dec 18, from 1 Ryd to H mass Rydberg
		 * >>chng 97 jan 10, had HLevNIonRyd but not defined yet */
		rfield.range[p.m_nqh][0] = HIONPOT;

		/* interpolated from Black 1987 */
		rfield.totpow[p.m_nqh] = (-18.517 + scale);

		++p.m_nqh;

		if( optimize.lgVarOn )
		{
			optimize.nvarxt[optimize.nparm] = 1;
			strcpy( optimize.chVarFmt[optimize.nparm], "TABLE ISM LOG %f");
			/*  pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			/*  the scale factor */
			optimize.vparm[0][optimize.nparm] = (realnum)scale;
			optimize.vincr[optimize.nparm] = 0.2f;
			++optimize.nparm;
		}
	}
	else if( p.nMatch("DRAI") )
	{
		rfield.lgMustBlockHIon = true;
		/* local ISM radiation field from equation 23
		 *>>refer	ISM	continuum	Draine & Bertoldi 1996 */
		for( i=0; i < NDRAINE; i++ )
		{
			rfield.tNu[rfield.nShape].emplace_back( Energy(tnudrn[i]) );
			rfield.tslop[rfield.nShape].emplace_back( (realnum)tsldrn[i] );
		}
		rfield.ncont[rfield.nShape] = NDRAINE;

		/* optional scale factor to change default luminosity
		 * from observed value
		 * assumed to be log if negative, and linear otherwise unless log option is present */
		scale = p.FFmtRead();
		if( scale > 0. && !p.nMatch(" LOG") )
			scale = log10(scale);

		/* this also sets continuum intensity*/
		if( p.m_nqh >= LIMSPC )
		{
			fprintf( ioQQQ, " %4ld is too many continua entered. Increase LIMSPC\n Sorry.\n", 
			  p.m_nqh );
			cdEXIT(EXIT_FAILURE);
		}

		/* check that stack of shape and luminosity specifications
		 * is parallel, stop if not - this happens is background comes
		 * BETWEEN another set of shape and luminosity commands */
		if( rfield.nShape != p.m_nqh )
		{
			fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );

		/* this is an isotropic radiation field */
		rfield.lgBeamed[p.m_nqh] = false;
		rfield.Illumination[p.m_nqh] = Illuminate::ISOTROPIC;

		/* continuum normalization given by flux density at first point,
		 * must set energy a bit higher to make sure it is within energy bounds
		 * that results from float arithmetic */
		rfield.range[p.m_nqh][0] = tnudrn[0]*1.01;

		/* this is f_nu at this first point */
		rfield.totpow[p.m_nqh] = tsldrn[0] + scale;

		if( optimize.lgVarOn )
		{
			optimize.nvarxt[optimize.nparm] = 1;
			strcpy( optimize.chVarFmt[optimize.nparm], "TABLE DRAINE LOG %f");
			/*  pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			/*  the scale factor */
			optimize.vparm[0][optimize.nparm] = (realnum)scale;
			optimize.vincr[optimize.nparm] = 0.2f;
			++optimize.nparm;
		}

		++p.m_nqh;
	}

	// match LINES to avoid confusion with LINEAR
	else if( p.nMatch("LINES") )
	{
		/* table lines command - way to check that lines within a data
		 * file are still valid */

		/* say that this is not a continuum command, so don't try to work with unallocated space */
		/* this is not a continuum source - it is to read a table of lines */
		lgNoContinuum = true;

		if( chLINE_LIST.size() > 0 )
		{
			fprintf(ioQQQ," sorry, only one table line per input stream\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* get file name within double quotes, if not present will use default
		 * return value of 1 indicates did not find double quotes on line */
		if( lgQuoteFound && chFile.length() > 0 )
			chLINE_LIST = chFile;
		else
			chLINE_LIST = "LineList_BLR.dat";

		// check if the file exists
		FILE* ioData = open_data( chLINE_LIST, "r", AS_TRY );
		if( ioData == NULL )
		{
			/* did not find file, abort */
			fprintf(ioQQQ,"\n DISASTER PROBLEM ParseTable could not find "
				"line list file %s\n", chLINE_LIST.c_str() );
			fprintf(ioQQQ," Please check the spelling of the file name and that it "
				"is in either the local or data directory.\n\n");
			cdEXIT(EXIT_FAILURE);
		}
		else
		{
			fclose(ioData);
		}
		/* actually reading the data is done in lines_table() */
	}

	else if( p.nMatch("POWE") )
	{
		/* simple power law continuum between 10 micron and 50 keV
		 *  option to read in any slope for the intermediate continuum */
		alpha = p.FFmtRead();

		/* default (no number on line) is f_nu proportional nu^-1 */
		if( p.lgEOL() )
			alpha = -1.;

		/* this is low energy for code */
		rfield.tNu[rfield.nShape].emplace_back( Energy(rfield.emm()) );
		/* and the value of the flux at this point (f_nu units)*/
		rfield.tslop[rfield.nShape].emplace_back( -5.f );

		lgLogSet = false;

		/* option to adjust sub-millimeter break */
		brakmm = p.FFmtRead();

		/* default is 10 microns */
		if( p.lgEOL() )
		{
			lgLogSet = true;
			brakmm = 9.115e-3;
		}

		else if( brakmm == 0. )
		{
			/* if second number on line is zero then set lower limit to
			 * low-energy limit of the code.  Also set linear mode,
			 * so that last number will also be linear. */
			lgLogSet = false;
			brakmm = rfield.tNu[rfield.nShape][0].Ryd()*1.001;
		}

		else if( brakmm < 0. )
		{
			/* if number is negative then this and next are logs */
			lgLogSet = true;
			brakmm = exp10(brakmm);
		}

		/* optional microns keyword - convert to Rydbergs */
		if( p.nMatch("MICR") )
			brakmm = RYDLAM / (1e4*brakmm);

		rfield.tNu[rfield.nShape].emplace_back( Energy(brakmm) );

		/* check whether this is a reasonable mm break */
		if( brakmm > 1. )
			fprintf(ioQQQ,
			" Check the order of parameters on this table power law command - the low-energy break of %f Ryd seems high to me.\n",
			brakmm );

		/* this is spectral slope, in F_nu units, between the low energy limit 
		 * and the break that may have been adjusted above 
		 * this is the slope appropriate for self-absorbed synchrotron, see eq 6.54, p.190
		 *>>refer	continuum	synchrotron	Rybicki, G. B., & Lightman, A.P. 1979, 
		 *>>refercon	Radiative Processes in Astrophysics (New York: Wiley)*/
		slopir = 5./2.;

		/* now extrapolate a flux at this energy using the flux entered for
		 * the first point, and this slope */
		rfield.tslop[rfield.nShape].emplace_back(
		  (realnum)(rfield.tslop[rfield.nShape][0] + 
		  slopir*log10(rfield.tNu[rfield.nShape][1].Ryd()/rfield.tNu[rfield.nShape][0].Ryd()))
			);

		/* option to adjust hard X-ray break */
		brakxr = p.FFmtRead();

		/* default is 50 keV */
		if( p.lgEOL() )
		{
			brakxr = 3676.;
		}

		else if( brakxr == 0. )
		{
			brakxr = rfield.egamry()/1.001;
		}

		else if( lgLogSet )
		{
			/* first number was negative this is a logs */
			brakxr = exp10(brakxr);
		}

		/* note that this second cutoff does not have the micron keyword */
		rfield.tNu[rfield.nShape].emplace_back( Energy(brakxr) );

		/* this is energy of the high-energy limit to code */
		rfield.tNu[rfield.nShape].emplace_back( Energy(rfield.egamry()) );

		/* >>chng 03 jul 19, check that upper energy is greater than lower energy,
		 * quit if this is not the case */
		if( brakmm >= brakxr )
		{
			fprintf( ioQQQ, " HELP!!  The lower energy for the power law, %f, is greater than the upper energy, %f. This is not possible.\n Sorry.\n",
				brakmm , brakxr );
			cdEXIT(EXIT_FAILURE);
		}

		/* alpha was first option on line, is slope of mid-range */
		rfield.tslop[rfield.nShape].emplace_back(
		  (realnum)(rfield.tslop[rfield.nShape][1] + 
		  alpha*log10(rfield.tNu[rfield.nShape][2].Ryd()/rfield.tNu[rfield.nShape][1].Ryd()))
			);

		/* high energy range is nu^-2 */
		slopxr = -2.;

		rfield.tslop[rfield.nShape].emplace_back(
		  (realnum)(rfield.tslop[rfield.nShape][2] + 
		  slopxr*log10(rfield.tNu[rfield.nShape][3].Ryd()/rfield.tNu[rfield.nShape][2].Ryd()))
			);

		/* following is number of portions of continuum */
		rfield.ncont[rfield.nShape] = 4;
	}

	else if( p.nMatch("READ") )
	{
		/* set up eventual read of table of points previously punched by code 
		 * get file name within double quotes, return as null terminated string
		 * also blank out original, chCard  version of name so do not trigger */
		if( !lgQuoteFound )
		{
			fprintf( ioQQQ, " Name of file must appear on TABLE READ.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* set flag saying really just read in continuum exactly as punched */
		strcpy( rfield.chSpType[rfield.nShape], "READ " );

		ReadTable(chFile);

		if ( p.nMatch("SCALE") )
		{
			/* optional scale factor to change default intensity from their value
			 * assumed to be log if negative, and linear otherwise 
			 * increment i to move past the 96 in the keyword */
			scale = p.FFmtRead();
			if( scale > 0. && !p.nMatch(" LOG"))
				scale = log10(scale);
			/* this also sets continuum intensity*/
			if( p.m_nqh >= LIMSPC )
			{
				fprintf( ioQQQ, " %ld is too many continua entered. Increase LIMSPC\n Sorry.\n", 
							p.m_nqh );
				cdEXIT(EXIT_FAILURE);
			}
			/* check that stack of shape and luminosity specifications
			 * is parallel, stop if not - this happens is background comes
			 * BETWEEN another set of shape and luminosity commands */
			if( rfield.nShape != p.m_nqh )
			{
				fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
			strcpy( rfield.chSpNorm[p.m_nqh], "FLUX" );
			// rfield.lgBeamed[p.m_nqh] = false;
			// rfield.Illumination[p.m_nqh] = Illuminate::ISOTROPIC;
		 
			rfield.range[p.m_nqh][0] = rfield.tNu[rfield.nShape][0].Ryd();
			double fmax = -70.;
			for( i=0; i < rfield.ncont[rfield.nShape]; i++ )
			{
				if (rfield.tFluxLog[rfield.nShape][i] > fmax)
				{
					fmax = rfield.tFluxLog[rfield.nShape][i];
					rfield.range[p.m_nqh][0] = rfield.tNu[rfield.nShape][i].Ryd();
				}
			}

			rfield.lgSphericalDilution[rfield.nShape] = true;

			// EN1RYD/FR1RYD == HPLANCK -- EN1RYD scaling already applied when saved
			rfield.totpow[p.m_nqh] = scale+fmax-log10(FR1RYD);
			
			++p.m_nqh;
		}
		/* number of spectra shapes that have been specified */
		++rfield.nShape;		
		
		return;
	}

	else if( p.nMatch("TLUSTY") && !p.nMatch("STAR") )
	{
		/* >>chng 04 nov 30, retired TABLE TLUSTY command, PvH */
		fprintf( ioQQQ, " The TABLE TLUSTY command is no longer supported.\n" );
		fprintf( ioQQQ, " Please use TABLE STAR TLUSTY instead. See Hazy for details.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( p.nMatch("SED") || lgKeyword )
	{
		if( !lgQuoteFound )
		{
			fprintf( ioQQQ, "PROBLEM in TABLE SED: No quotes were found.\n"
					 "The SED table file must be designated in quotes.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		string chPath = "SED" + cpu.i().chDirSeparator() + chFile;
		DataParser d(chPath, ES_STARS_ONLY);

		/* read data and parse optional keywords */
		long entryNum = 0;
		vector<double> tnuInput;
		vector<double> tslopInput;
		const char *chEnergyUnits = nullptr;
		double nuPower = 0.;
		bool lgExtrapolate = false;

		while( d.getline() )
		{
			if( d.lgEODMarker() )
				break;

			double freq, flux;
			d.getToken(freq);
			if( freq <= 0. )
				d.errorAbort("the freq/wavl in the TABLE SED file must be positive");
			tnuInput.emplace_back(freq);
			if( entryNum > 1 && (tnuInput[entryNum-2]-tnuInput[entryNum-1])*
				(tnuInput[entryNum-1]-tnuInput[entryNum]) <= 0. )
				d.errorAbort("freq/wavl values in TABLE SED file must be strictly increasing or decreasing");

			d.getToken(flux);
			if( flux <= 0. )
				d.errorAbort("the flux in the TABLE SED file must be positive");
			tslopInput.emplace_back(log10(flux));
			++entryNum;

			// parse possible keywords
			string key;
			while( d.getKeywordOptional(key) )
			{
				if( key.substr(0,4) == "UNIT" )
				{
					/* change photon units */
					d.getKeyword(key);
					// pad string with spaces to make sure it can be recognized
					key = " " + key + " ";
					chEnergyUnits = StandardEnergyUnitNoAbort(key.c_str());
					if( chEnergyUnits == nullptr )
						d.errorAbort("no energy / wavelength unit was recognized");
				}
				else if( key.substr(0,4) == "NUFN" )
				{
					/* we want F_nu - but entered was nu F_nu */
					nuPower = 1.;
				}
				else if( key.substr(0,4) == "FLAM" )
				{
					/* we want F_nu - but entered was F_lambda */
					nuPower = 2.;
				}
				else if( key.substr(0,4) == "EXTR" )
				{
					/* default is to use exactly the specified continuum, extrapolate option
					 * will extrapolate lowest energy point to low-energy limit of the code */
					lgExtrapolate = true;
				}
				else
				{
					d.warning("keyword not recognized");
				}
			}
		}

		if( entryNum < 2 )
		{
			fprintf(ioQQQ, "PROBLEM in TABLE SED: less than two data pairs found.\n");
			cdEXIT(EXIT_FAILURE);
		}

		// Unit Conversion of frequency, we want linear Rydbergs
		if( chEnergyUnits != nullptr )
		{
			for( long i=0; i < entryNum; i++ )
			{
				Energy unitChange;
				unitChange.set(tnuInput[i], chEnergyUnits);
				tnuInput[i] = unitChange.Ryd();
			}
		}

		/* we want to specify log F_nu - this converts log nuFnu to log F_nu */
		if( nuPower > 0. )
		{
			for( long i=0; i < entryNum; i++ )
				tslopInput[i] -= nuPower*log10(tnuInput[i]);
		}

		if( tnuInput[0] < tnuInput[1] )
		{
			for( long i=0; i < entryNum; i++ )
			{
				rfield.tNu[rfield.nShape].emplace_back( Energy(tnuInput[i]) );
				rfield.tslop[rfield.nShape].emplace_back( (realnum)tslopInput[i] );
			}
		}
		else
		{
			for( long i=0; i < entryNum; i++ )
			{
				rfield.tNu[rfield.nShape].emplace_back( Energy(tnuInput[entryNum-i-1]) );
				rfield.tslop[rfield.nShape].emplace_back( (realnum)tslopInput[entryNum-i-1] );
			}
		}
		rfield.ncont[rfield.nShape] = entryNum;

		/* option to extrapolate lowest energy point specified to low-energy limit of code */
		if( lgExtrapolate )
			resetBltin( rfield.tNu[rfield.nShape].data(), rfield.tslop[rfield.nShape].data(), true );
	}

	/* >>chng 06 jul 10, retired TABLE STARBURST command, PvH */

	else if( p.nMatch("STAR") )
	{
		string chMetalicity, chODFNew, chVaryFlag;
		bool lgHCa = false, lgHNi = false;
		long nval, ndim=0;
		double Tlow = -1., Thigh = -1.;
		double val[MDIM];

		/* >>chng 06 jun 22, add support for 3 and 4-dimensional grids, PvH */
		if( p.nMatchErase("1-DI") )
			ndim = 1;
		else if( p.nMatchErase("2-DI") )
			ndim = 2;
		else if( p.nMatchErase("3-DI") )
			ndim = 3;
		else if( p.nMatchErase("4-DI") )
			ndim = 4;

		if( ndim != 0 )
		{
			/* remember keyword for possible use in optimization command */
			ostringstream oss;
			oss << ndim << "-DIM";
			chVaryFlag = oss.str();
		}

		/* time option to vary only some continua with time */
		rfield.lgTimeVary[p.m_nqh] = p.nMatch( "TIME" );
			
		static const string table[][2] = {
			{"Z+1.0 ", "p10"},
			{"Z+0.75", "p075"},
			{"Z+0.5 ", "p05"},
			{"Z+0.4 ", "p04"},
			{"Z+0.3 ", "p03"},
			{"Z+0.25", "p025"},
			{"Z+0.2 ", "p02"},
			{"Z+0.1 ", "p01"},
			{"Z+0.0 ", "p00"},
			{"Z-0.1 ", "m01"},
			{"Z-0.2 ", "m02"},
			{"Z-0.25", "m025"},
			{"Z-0.3 ", "m03"},
			{"Z-0.4 ", "m04"},
			{"Z-0.5 ", "m05"},
			{"Z-0.7 ", "m07"},
			{"Z-0.75", "m075"},
			{"Z-1.0 ", "m10"},
			{"Z-1.3 ", "m13"},
			{"Z-1.5 ", "m15"},
			{"Z-1.7 ", "m17"},
			{"Z-2.0 ", "m20"},
			{"Z-2.5 ", "m25"},
			{"Z-3.0 ", "m30"},
			{"Z-3.5 ", "m35"},
			{"Z-4.0 ", "m40"},
			{"Z-4.5 ", "m45"},
			{"Z-5.0 ", "m50"},
			{"Z-INF ", "m99"}
		};
	  
		chMetalicity = "p00";	// default
		for (i=0; i < (long)(sizeof(table)/sizeof(*table)); ++i)
		{
			if (p.nMatchErase(table[i][0].c_str()))
			{
				chVaryFlag = table[i][0];
				chMetalicity = table[i][1];
				break;
			}
		}


		/* there are two abundance sets (solar and halo) for CoStar and Rauch H-Ca/H-Ni models.
		 * If halo keyword appears use halo, else use solar unless other metalicity was requested */
		bool lgHalo = p.nMatch("HALO");
		bool lgSolar = ( !lgHalo && chMetalicity == "p00" );

		/* >>chng 05 oct 27, added support for PG1159 Rauch models, PvH */
		bool lgPG1159 = p.nMatchErase("PG1159");

		bool lgList = p.nMatch("LIST");

		if( p.nMatch("AVAI") )
		{
			AtmospheresAvail();
			cdEXIT(EXIT_SUCCESS);
		}

		/* now that all the keywords are out of the way, scan numbers from line image */
		for( nval=0; nval < MDIM; nval++ )
		{
			val[nval] = p.FFmtRead();
			if( p.lgEOL() )
				break;
		}

		if( nval == 0 && !lgList )
		{
			fprintf( ioQQQ, " The stellar temperature MUST be entered.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* match on "LOG " rather than " LOG" to avoid confusion with "LOG(G)" */
		if( p.nMatch("LOG ") )
		{
			if( val[0] < log10(BIGDOUBLE) )
				val[0] = exp10(val[0]);
			else
			{
				fprintf(ioQQQ," DISASTER the log of the parameter was specified, "
						"but the numerical value is too large.\n Sorry.\n\n");
				cdEXIT(EXIT_FAILURE );
			}
		}

		if( lgQuoteFound )
		{
			nstar = GridInterpolate(val,&nval,&ndim,chFile,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("ATLA") )
		{
			/* this sub-branch added by Kevin Volk, to read in large
			 * grid of Kurucz atmospheres */
			/* >>chng 05 nov 19, option TRACE is no longer supported, PvH */

			if( p.nMatch("ODFN") )
				chODFNew = "_odfnew";
			else
				chODFNew = "";

			/* >>chng 05 nov 19, add support for non-solar metalicities, PvH */
			nstar = AtlasInterpolate(val,&nval,&ndim,chMetalicity,chODFNew,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("COST") )
		{
			/* >>chng 99 apr 30, added CoStar stellar atmospheres */
			/* this subbranch reads in CoStar atmospheres, no log(g),
			 * but second parameter is age sequence, a number between 1 and 7,
			 * default is 1 */
			if( p.nMatch("INDE") )
			{
				imode = IM_COSTAR_TEFF_MODID;
				if( nval == 1 )
				{
					val[1] = 1.;
					nval++;
				}

				/* now make sure that second parameter is within allowed range -
				 * we do not have enough information at this time to verify temperature */
				if( val[1] < 1. || val[1] > 7. )
				{
					fprintf( ioQQQ, " The second number must be the id sequence number, 1 to 7.\n" );
					fprintf( ioQQQ, " reset to 1.\n" );
					val[1] = 1.;
				}
			}
			else if( p.nMatch("ZAMS") )
			{
				imode = IM_COSTAR_MZAMS_AGE;
				if( nval == 1 )
				{
					fprintf( ioQQQ, " A second number, the age of the star, must be present.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
			else if( p.nMatch(" AGE") )
			{
				imode = IM_COSTAR_AGE_MZAMS;
				if( nval == 1 )
				{
					fprintf( ioQQQ, " A second number, the ZAMS mass of the star, must be present.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
			else
			{
				if( nval == 1 )
				{
					imode = IM_COSTAR_TEFF_MODID;
					/* default is to use ZAMS models, i.e. use index 1 */
					val[1] = 1.;
					nval++;
				}
				else
				{
					/* Caution: the code in CoStarInterpolate implicitly
					 * assumes that the dependence of log(g) with age is
					 * strictly monotonic. As of June 1999 this is the case. */
					imode = IM_COSTAR_TEFF_LOGG;
				}
			}

			if( ! ( lgSolar || lgHalo ) )
			{
				fprintf( ioQQQ, " Please choose SOLAR or HALO abundances.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			nstar = CoStarInterpolate(val,&nval,&ndim,imode,lgHalo,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("KURU") )
		{
			nstar = Kurucz79Interpolate(val,&nval,&ndim,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("MIHA") )
		{
			nstar = MihalasInterpolate(val,&nval,&ndim,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("RAUC") )
		{
			if( ! ( lgSolar || lgHalo ) )
			{
				fprintf( ioQQQ, " Please choose SOLAR or HALO abundances.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* >>chng 97 aug 23, K Volk added Rauch stellar atmospheres */
			/* >>chng 05 oct 27, added support for PG1159 Rauch models, PvH */
			/* >>chng 06 jun 26, added support for H, He, H+He Rauch models, PvH */
			if( p.nMatch("H-CA") || p.nMatch(" OLD") )
			{
				lgHCa = true;
				nstar = RauchInterpolateHCa(val,&nval,&ndim,lgHalo,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("HYDR") )
			{
				nstar = RauchInterpolateHydr(val,&nval,&ndim,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("HELI") )
			{
				nstar = RauchInterpolateHelium(val,&nval,&ndim,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("H+HE") )
			{
				nstar = RauchInterpolateHpHe(val,&nval,&ndim,lgList,&Tlow,&Thigh);
			}
			else if( lgPG1159 ) /* the keyword PG1159 was already matched and erased above */
			{
				nstar = RauchInterpolatePG1159(val,&nval,&ndim,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("CO W") )
			{
				nstar = RauchInterpolateCOWD(val,&nval,&ndim,lgList,&Tlow,&Thigh);
			}
			else
			{
				lgHNi = true;
				nstar = RauchInterpolateHNi(val,&nval,&ndim,lgHalo,lgList,&Tlow,&Thigh);
			}
		}

		else if( p.nMatch("TLUS") )
		{
			if( p.nMatch("OBST") )
			{
				/* >>chng 09 dec 24, this sub-branch added to read in
				 * merged Tlusty O-star & B-star atmospheres, PvH */
				nstar = TlustyInterpolate(val,&nval,&ndim,TL_OBSTAR,chMetalicity,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("BSTA") )
			{
				/* >>chng 06 oct 19, this sub-branch added to read in
				 * large 2006 grid of Tlusty B-star atmospheres, PvH */
				nstar = TlustyInterpolate(val,&nval,&ndim,TL_BSTAR,chMetalicity,lgList,&Tlow,&Thigh);
			}
			else if( p.nMatch("OSTA") )
			{
				/* >>chng 05 nov 21, this sub-branch added to read in
				 * large 2002 grid of Tlusty O-star atmospheres, PvH */
				nstar = TlustyInterpolate(val,&nval,&ndim,TL_OSTAR,chMetalicity,lgList,&Tlow,&Thigh);
			}
			else
			{
				fprintf( ioQQQ, " There must be a third key on TABLE STAR TLUSTY;" );
				fprintf( ioQQQ, "  the options are OBSTAR, BSTAR, OSTAR.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		else if( p.nMatch("WERN") )
		{
			/* this subbranch added by Kevin Volk, to read in large
			 * grid of hot PN atmospheres computed by Klaus Werner */
			nstar = WernerInterpolate(val,&nval,&ndim,lgList,&Tlow,&Thigh);
		}

		else if( p.nMatch("WMBA") )
		{
			/* >>chng 06 jun 27, this subbranch added to read in
			 * grid of hot atmospheres computed by Pauldrach */
			nstar = WMBASICInterpolate(val,&nval,&ndim,lgList,&Tlow,&Thigh);
		}

		else
		{
			fprintf( ioQQQ, " There must be a second key on TABLE STAR;" );
			fprintf( ioQQQ, "  the options are ATLAS, KURUCZ, MIHALAS, RAUCH, WERNER, and WMBASIC.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* set flag saying really just read in continuum exactly as punched */
		strcpy( rfield.chSpType[rfield.nShape], "VOLK " );

		rfield.ncont[rfield.nShape] = nstar;

		/* vary option */
		if( optimize.lgVarOn )
		{
			optimize.vincr[optimize.nparm] = (realnum)0.1;

			if( lgQuoteFound )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = nval;
				sprintf( optimize.chVarFmt[optimize.nparm], "TABLE STAR \"%s\"", chFile.c_str() );
				strcat( optimize.chVarFmt[optimize.nparm], " %f LOG" );
				for( i=1; i < nval; i++ )
					strcat( optimize.chVarFmt[optimize.nparm], " %f" );
			}

			if( p.nMatch("ATLA") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = ndim;
				strcpy( optimize.chVarFmt[optimize.nparm], "TABLE STAR ATLAS " );
				strcat( optimize.chVarFmt[optimize.nparm], chVaryFlag.c_str() );
				if( p.nMatch("ODFN") )
					strcat( optimize.chVarFmt[optimize.nparm], " ODFNEW" );

				strcat( optimize.chVarFmt[optimize.nparm], " %f LOG %f" );

				if( ndim == 3 )
					strcat( optimize.chVarFmt[optimize.nparm], " %f" );

			}

			else if( p.nMatch("COST") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = 2;

				strcpy( optimize.chVarFmt[optimize.nparm], "TABLE STAR COSTAR " );
				if( lgHalo )
					strcat( optimize.chVarFmt[optimize.nparm], "HALO " );

				if( imode == IM_COSTAR_TEFF_MODID )
				{
					strcat( optimize.chVarFmt[optimize.nparm], "%f LOG , INDEX %f" );
				}
				else if( imode == IM_COSTAR_TEFF_LOGG )
				{
					strcat( optimize.chVarFmt[optimize.nparm], "%f LOG %f" );
				}
				else if( imode == IM_COSTAR_MZAMS_AGE )
				{
					strcat( optimize.chVarFmt[optimize.nparm], "ZAMS %f LOG %f" );
				}
				else if( imode == IM_COSTAR_AGE_MZAMS )
				{
					strcat( optimize.chVarFmt[optimize.nparm], "AGE %f LOG %f" );
					optimize.vincr[optimize.nparm] = (realnum)0.5;
				}
			}

			else if( p.nMatch("KURU") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], 
					"TABLE STAR KURUCZ %f LOG" );
			}

			else if( p.nMatch("MIHA") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( optimize.chVarFmt[optimize.nparm], 
					"TABLE STAR MIHALAS %f LOG" );
			}

			else if( p.nMatch("RAUC") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = ndim;

				strcpy( optimize.chVarFmt[optimize.nparm], "TABLE STAR RAUCH " );

				if( p.nMatch("HYDR") )
					strcat( optimize.chVarFmt[optimize.nparm], "HYDR " );
				else if( p.nMatch("HELI") )
					strcat( optimize.chVarFmt[optimize.nparm], "HELIUM " );
				else if( p.nMatch("H+HE") )
					strcat( optimize.chVarFmt[optimize.nparm], "H+HE " );
				else if( lgPG1159 )
					strcat( optimize.chVarFmt[optimize.nparm], "PG1159 " );
				else if( p.nMatch("CO W") )
					strcat( optimize.chVarFmt[optimize.nparm], "CO WD " );
				else if( lgHCa )
					strcat( optimize.chVarFmt[optimize.nparm], "H-CA " );

				if( ( lgHCa || lgHNi ) && ndim == 2 )
				{
					if( lgHalo )
						strcat( optimize.chVarFmt[optimize.nparm], "HALO " );
					else
						strcat( optimize.chVarFmt[optimize.nparm], "SOLAR " );
				}

				strcat( optimize.chVarFmt[optimize.nparm], "%f LOG %f" );

				if( ndim == 3 )
				{
					if( p.nMatch("H+HE") )
						strcat( optimize.chVarFmt[optimize.nparm], " %f" );
					else
						strcat( optimize.chVarFmt[optimize.nparm], " %f 3-DIM" );
				}
			}

			if( p.nMatch("TLUS") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = ndim;
				strcpy( optimize.chVarFmt[optimize.nparm], "TABLE STAR TLUSTY " );
				if( p.nMatch("OBST") )
					strcat( optimize.chVarFmt[optimize.nparm], "OBSTAR " );
				else if( p.nMatch("BSTA") )
					strcat( optimize.chVarFmt[optimize.nparm], "BSTAR " );
				else if( p.nMatch("OSTA") )
					strcat( optimize.chVarFmt[optimize.nparm], "OSTAR " );
				else
					TotalInsanity();
				strcat( optimize.chVarFmt[optimize.nparm], chVaryFlag.c_str() );
				strcat( optimize.chVarFmt[optimize.nparm], " %f LOG %f" );

				if( ndim == 3 )
					strcat( optimize.chVarFmt[optimize.nparm], " %f" );
			}

			else if( p.nMatch("WERN") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = 2;
				strcpy( optimize.chVarFmt[optimize.nparm], 
					"TABLE STAR WERNER %f LOG %f" );
			}

			else if( p.nMatch("WMBA") )
			{
				/* this is number of parameters to feed onto the input line */
				optimize.nvarxt[optimize.nparm] = 3;
				strcpy( optimize.chVarFmt[optimize.nparm], 
					"TABLE STAR WMBASIC %f LOG %f %f" );
			}

			if( rfield.lgTimeVary[p.m_nqh] )
				strcat( optimize.chVarFmt[optimize.nparm], " TIME" );

			/* pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;

			ASSERT( nval <= LIMEXT );

			/* log of temp will be pointer */
			optimize.vparm[0][optimize.nparm] = (realnum)log10(val[0]);
			for( i=1; i < nval; i++ )
				optimize.vparm[i][optimize.nparm] = (realnum)val[i];

			/* following are upper and lower limits to temperature range */
			optimize.varang[optimize.nparm][0] = (realnum)log10(Tlow);
			optimize.varang[optimize.nparm][1] = (realnum)log10(Thigh);

			/* finally increment this counter */
			++optimize.nparm;
		}
	}

	else
	{
		fprintf( ioQQQ, " There MUST be a keyword on this line. The keys are:"
			"AGN, DRAINE, HM96, HM05, HM12, KS19, _ISM, LINES, POWERlaw, "
			"READ, SED, and STAR.\n  Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* table star and table hm05 / hm12 are special cases
	 * they are not really tables at all,
	 * so return if chSpType is "VOLK " */
	if( strcmp(rfield.chSpType[rfield.nShape],"VOLK ") == 0 )
	{ 
		++rfield.nShape;
		return;
	}

	/* this flag set true if we did not parse a continuum source, thus creating
	 * the arrays that are needed - return at this point */
	if( lgNoContinuum )
		return;

	// must have at least 2 points to interpolate
	ASSERT( rfield.ncont[rfield.nShape] > 1 );

	ASSERT( rfield.tNu[rfield.nShape].size() == (size_t)rfield.ncont[rfield.nShape] );
	ASSERT( rfield.tslop[rfield.nShape].size() == (size_t)rfield.ncont[rfield.nShape] );


	ASSERT( rfield.tNu[rfield.nShape][0].Ryd() > 0. );

	/* tFluxLog will be log(fnu) at that spot, read into tslop */
	for( i=0; i < rfield.ncont[rfield.nShape]; i++ )
	{
		rfield.tFluxLog[rfield.nShape].emplace_back( rfield.tslop[rfield.nShape][i] );
	}

	ASSERT( rfield.tFluxLog[rfield.nShape].size() == (size_t)rfield.ncont[rfield.nShape] );

	/* at this point tslop is the log of the intensity */
	for( i=0; i < rfield.ncont[rfield.nShape]-1; i++ )
	{
		rfield.tslop[rfield.nShape][i] = 
			(realnum)((rfield.tslop[rfield.nShape][i+1] - 
		  rfield.tslop[rfield.nShape][i])/
		  log10(rfield.tNu[rfield.nShape][i+1].Ryd()/
		  rfield.tNu[rfield.nShape][i].Ryd()));
	}
	rfield.tslop[rfield.nShape][rfield.ncont[rfield.nShape]-1] = 0.f;

	if( trace.lgConBug && trace.lgTrace )
	{
		fprintf( ioQQQ, " Table for this continuum; TNU, TFAC, TSLOP\n" );
		for( i=0; i < rfield.ncont[rfield.nShape]; i++ )
		{
			fprintf( ioQQQ, "%12.4e %12.4e %12.4e\n", rfield.tNu[rfield.nShape][i].Ryd(), 
			  rfield.tFluxLog[rfield.nShape][i], rfield.tslop[rfield.nShape][i] );
		}
	}

	/* renormalize continuum so that quantity passes through unity at 1 Ryd
	 * lgHit will be set false when we get a hit in following loop,
	 * then test made to make sure it happened*/
	lgHit = false;
	/*following will be reset when hit occurs*/
	fac = -DBL_MAX;
	/* >>chng 04 mar 16, chng loop so breaks when hit, previously had cycled
	 * until last point reached, so last point always used */
	for( i=0; i < rfield.ncont[rfield.nShape] && !lgHit; i++ )
	{
		if( rfield.tNu[rfield.nShape][i].Ryd() > 1. )
		{
			fac = rfield.tFluxLog[rfield.nShape][i];
			lgHit = true;
		}
	}

	if( lgHit )
	{
		/* do the renormalization */
		for( i=0; i < rfield.ncont[rfield.nShape] ; i++ )
			rfield.tFluxLog[rfield.nShape][i] -= (realnum)fac;
	}

	++rfield.nShape;

	return;
}



/*ZeroContin store sets of built in continua */
STATIC void ZeroContin(void)
{

	DEBUG_ENTRY( "ZeroContin()" );

	/* Draine 1978 ISM continuum */
	/* freq in ryd */
	tnudrn[0] = 0.3676;
	tnudrn[1] = 0.4144;
	tnudrn[2] = 0.4558;
	tnudrn[3] = 0.5064;
	tnudrn[4] = 0.5698;
	tnudrn[5] = 0.6511;
	tnudrn[6] = 0.7012;
	tnudrn[7] = 0.7597;
	tnudrn[8] = 0.8220;
	tnudrn[9] = 0.9116;
	tnudrn[10] = 0.9120;
	tnudrn[11] = 0.9306;
	tnudrn[12] = 0.9600;
	tnudrn[13] = 0.9806;
	/* >>chng 05 aug 03, this energy is so close to 1 ryd that it spills over
	 * into the H-ionizing continuum - move down to 0.99 */
	/* >>chng 05 aug 08, this destabilized pdr_leiden_hack_v4 (!) so put back to
	 * original value and include extinguish command */
	tnudrn[14] = 0.9999;

	/* f_nu from equation 23 of 
	 * >>refer	ism	field	Draine, B.T. & Bertoldi, F. 1996, ApJ, 468, 269 */
	tsldrn[0] = -17.8063;
	tsldrn[1] = -17.7575;
	tsldrn[2] = -17.7268;
	tsldrn[3] = -17.7036;
	tsldrn[4] = -17.6953;
	tsldrn[5] = -17.7182;
	tsldrn[6] = -17.7524;
	tsldrn[7] = -17.8154;
	tsldrn[8] = -17.9176;
	tsldrn[9] = -18.1675;
	tsldrn[10] = -18.1690;
	tsldrn[11] = -18.2477;
	tsldrn[12] = -18.4075;
	tsldrn[13] = -18.5624;
	tsldrn[14] = -18.7722;

	/* generic AGN continuum taken from 
	 * >>refer	AGN	cont	Mathews and Ferland ApJ Dec15 '87
	 * except self absorption at 10 microns, nu**-2.5 below that */
	tnuagn[0] = Energy(1e-5);
	tnuagn[1] = Energy(9.12e-3);
	tnuagn[2] = Energy(.206);
	tnuagn[3] = Energy(1.743);
	tnuagn[4] = Energy(4.13);
	tnuagn[5] = Energy(26.84);
	tnuagn[6] = Energy(7353.);
	tnuagn[7] = Energy(7.4e6);

	tslagn[0] = -3.388_r;
	tslagn[1] = 4.0115_r;
	tslagn[2] = 2.6576_r;
	tslagn[3] = 2.194_r;
	tslagn[4] = 1.819_r;
	tslagn[5] = -.6192_r;
	tslagn[6] = -2.326_r;
	tslagn[7] = -7.34_r;
	resetBltin( tnuagn , tslagn , true );

	/* interstellar radiation field from Black 1987, "Interstellar Processes"
	 * table of log nu, log nu*fnu taken from his figure 2 */
	/* >>chng 99 jun 14 energy range lowered to 1e-8 ryd */
	tnuism[0] = 6.00;
	/*tnuism[0] = 9.00;*/
	tnuism[1] = 10.72;
	tnuism[2] = 11.00;
	tnuism[3] = 11.23;
	tnuism[4] = 11.47;
	tnuism[5] = 11.55;
	tnuism[6] = 11.85;
	tnuism[7] = 12.26;
	tnuism[8] = 12.54;
	tnuism[9] = 12.71;
	tnuism[10] = 13.10;
	tnuism[11] = 13.64;
	tnuism[12] = 14.14;
	tnuism[13] = 14.38;
	tnuism[14] = 14.63;
	tnuism[15] = 14.93;
	tnuism[16] = 15.08;
	tnuism[17] = 15.36;
	tnuism[18] = 15.43;
	tnuism[19] = 16.25;
	tnuism[20] = 17.09;
	tnuism[21] = 18.00;
	tnuism[22] = 23.00;
	/* these are log nu*Fnu */
	fnuism[0] = -16.708;
	/*fnuism[0] = -7.97;*/
	fnuism[1] = -2.96;
	fnuism[2] = -2.47;
	fnuism[3] = -2.09;
	fnuism[4] = -2.11;
	fnuism[5] = -2.34;
	fnuism[6] = -3.66;
	fnuism[7] = -2.72;
	fnuism[8] = -2.45;
	fnuism[9] = -2.57;
	fnuism[10] = -3.85;
	fnuism[11] = -3.34;
	fnuism[12] = -2.30;
	fnuism[13] = -1.79;
	fnuism[14] = -1.79;
	fnuism[15] = -2.34;
	fnuism[16] = -2.72;
	fnuism[17] = -2.55;
	fnuism[18] = -2.62;
	fnuism[19] = -5.68;
	fnuism[20] = -6.45;
	fnuism[21] = -6.30;
	fnuism[22] = -11.3;

	return;
}

/*lines_table invoked by table lines command, check if we can find all lines in a given list */
int lines_table()
{
	DEBUG_ENTRY( "lines_table()" );

	if( chLINE_LIST.empty() )
		return 0;

	vector<LineID> lineids;

	long nLINE_TABLE = cdGetLineList( chLINE_LIST, lineids );

	// the check if the file exists has already been done by the parser
	if( nLINE_TABLE == 0 )
		return 0;

	fprintf( ioQQQ, "lines_table checking lines within data table %s\n", chLINE_LIST.c_str() );
	int miss = 0;

	for( long n=0; n < nLINE_TABLE; ++n )
	{
		if( LineSave.findline( lineids[n] ) <= 0 )
		{
			++miss;
			fprintf( ioQQQ, "lines_table in parse_table.cpp did not find line " );
			prt_line_err(ioQQQ, lineids[n]);
		}
	}
	if( miss )
	{
		/* this is so that we pick up problem automatically */
		fprintf( ioQQQ, "  BOTCHED MONITORS!!!   Botched Monitors!!!"
				  " lines_table could not find a total of %i lines\n\n", miss );
	}
	else
	{
		fprintf( ioQQQ, "lines_table found all lines\n\n" );
	}
	return miss;
}

void clear_lines_table()
{
	/* prevent spurious complaints about more than one table lines command being present in a grid */
	chLINE_LIST.clear();
}

inline double getHexDouble(DataParser& d)
{
	union {
		double x;
		uint64 i;
	} u;

	string s, s2;
	d.getToken(s);
	d.getToken(s2);
	s += s2;
	istringstream iss(s);
	iss >> hex >> u.i;
	return u.x;
}

/*ReadTable called by TABLE READ to read in continuum from SAVE TRANSMITTED CONTINUUM */
STATIC void ReadTable(const string& fnam)
{
	DEBUG_ENTRY( "ReadTable()" );

	DataParser d(fnam, ES_NONE);

	/* line 3: the magic number */
	d.getline();
	d.checkMagic(VERSION_TRNCON);

	/* line 4: the MD5 checksum */
	string md5sum;
	d.getline();
	d.getToken(md5sum);
	
	double mesh_lo, mesh_hi;
	/* line 5: the lower limit of the energy grid */
	d.getline();
	mesh_lo = getHexDouble(d);
	
	/* line 6: the upper limit of the energy grid */
	d.getline();
	mesh_hi = getHexDouble(d);
	
	if( md5sum != rfield.mesh_md5sum() ||
	    !fp_equal_tol( mesh_lo, rfield.emm(), 1.e-11*rfield.emm() ) ||
	    !fp_equal_tol( mesh_hi, rfield.egamry(), 1.e-7*rfield.egamry() ) )
	{
		fprintf( ioQQQ, " the input continuum file has an incompatible energy grid.\n" );
		fprintf( ioQQQ, " please recreate this file using the SAVE TRANSMITTED CONTINUUM command.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* line 7: the energy grid resolution scale factor */
	d.getline();
	rfield.RSFCheck[rfield.nShape] = getHexDouble(d);
	
	/* line 8: the outer radius of the saved model */
	d.getline();
	d.getToken(rfield.TableRadius[rfield.nShape]);

	/* line 9: the number of frequency grid points contained in the file */
	long nflux;
	d.getline();
	d.getToken(nflux);
	
	/* now read in the file of numbers */
	long i = 0;
	/* keep reading until we hit eof */
	while( d.getline() )
	{
		double help[2];
		d.getToken(help[0]);
		d.getToken(help[1]);
		// check that frequency grid matches, frequencies are printed with 6 significant digits
		if( !fp_equal_tol(help[0],rfield.anu(i),1e-5*rfield.anu(i)) )
		{
			fprintf(ioQQQ,"\n\nPROBLEM in TABLE READ continuum frequencies: point %li should "
				"have %e and has %e\n", (unsigned long)i, rfield.anu(i), help[0] );
			fprintf(ioQQQ,"Please recreate this file.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		rfield.tNu[rfield.nShape].emplace_back( Energy(help[0]) );
		// Convert to standard FluxLog units
		if (help[1] > 0.0)
			rfield.tFluxLog[rfield.nShape].emplace_back( (realnum)log10(help[1]/help[0]) );
		else
			rfield.tFluxLog[rfield.nShape].emplace_back( -70_r );
		++i;
	}
	/* put pointer at last good value */
	rfield.ncont[rfield.nShape] = i;

	ASSERT( rfield.tNu[rfield.nShape].size() == (size_t)rfield.ncont[rfield.nShape] );
	ASSERT( rfield.tFluxLog[rfield.nShape].size() == (size_t)rfield.ncont[rfield.nShape] );

	/* check that sane number of values entered */
	if( rfield.ncont[rfield.nShape] != nflux )
	{
		fprintf( ioQQQ, " the input continuum file has the wrong number of points: %ld\n", 
		  rfield.ncont[rfield.nShape] );
		fprintf( ioQQQ, " please recreate this file using the SAVE TRANSMITTED CONTINUUM command.\n" );
		cdEXIT(EXIT_FAILURE);
	}
}
