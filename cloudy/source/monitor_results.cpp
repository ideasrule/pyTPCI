/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*InitMonitorResults, this must be called first, done at startup of ParseCommands*/
/*ParseMonitorResults - parse input stream */
/*lgCheckMonitors, checks asserts, last thing cloudy calls, returns true if all are ok, false if problems */
#include "cddefines.h"
#include "input.h"
#include "conv.h"
#include "optimize.h"
#include "iso.h"
#include "called.h"
#include "atmdat.h"
#include "hcmap.h"
#include "thermal.h"
#include "pressure.h"
#include "struc.h"
#include "wind.h"
#include "h2.h"
#include "colden.h"
#include "dense.h"
#include "lines.h"
#include "secondaries.h"
#include "radius.h"
#include "version.h"
#include "hmi.h"
#include "prt.h"
#include "grainvar.h"
#include "cddrive.h"
#include "elementnames.h"
#include "monitor_results.h"
#include "parser.h"
#include "mole.h"
#include "rfield.h"
#include "lines_service.h"
#include "service.h"
#include "generic_state.h"

bool lgMonitorsOK , lgBigBotch, lgPrtSciNot;
t_monitorresults MonitorResults;

/* flag to remember that InitMonitorResults was called */
static bool lgInitDone=false , 
	/* will be set true when space for asserts is allocated */
	lgSpaceAllocated=false;

/* number of asserts we can handle, used in dim of space */
static const int NASSERTS = 450;

/* default relative error for monitored physical quantities */
static realnum ErrorDefault;

/* default relative error for monitored performance metrics */
static realnum ErrorDefaultPerformance;

static vector<string> chAssertLineLabel;

/* these are the lines read from MONITOR LINE command */
static vector<LineID> lineids;

/* type of line integral to report:
 * 0: intrinsic, 1: emergent, 2: cumul intrinsic, 3: cumul emergent */
static vector<int> iLineType;

static vector<vector<TransitionProxy>*> assertBlends;

/* this will be = for equal, < or > for limit */
static vector<char> chAssertLimit;

/* this will be a two character label identifying which type of monitor */
static vector<string> chAssertType;

/* the values and error in the asserted quantity */
static vector<double> AssertQuantity, AssertQuantity2 ,AssertError;
static multi_arr<double,2> Param;

/* this flag says where we print linear or log quantity */
static vector<bool> lgQuantityLog;
static long nAsserts=0;
static vector<realnum> wavelength;

static vector<string> strAssertSpecies;

void PrtOneMonitor( FILE *ioMONITOR, const string& chAssertType, const string& chAssertLineLabel, 
		const realnum wavelength, const int iLineType, const double PredQuan, const char chAssertLimit, const double AssertQuantity, 
		const double RelError, const double AssertError, const bool lg1OK, const bool lgQuantityLog, const bool lgFound );

inline double ForcePass( char chAssertLimit1 )
{
	// force monitors to pass by returning a safe value
	if( chAssertLimit1 == '=' )
		return 0.;
	else if( chAssertLimit1 == '<' )
		return 1.;
	else if( chAssertLimit1 == '>' )
		return -1.;
	else
		TotalInsanity();
}

/*======================================================================*/
/*InitMonitorResults, this must be called first, done at startup of ParseCommands*/
void InitMonitorResults(void)
{
	/* set flag that init was called, and set number of asserts to zero.
	 * this is done by ParseComments for every model, even when no asserts will
	 * be done, so do not allocate space at this time */
	lgInitDone = true;
	nAsserts = 0;

	// following occur in hazy1 default for monitor commands
	// default error, changed with MONITOR SET ERROR
	ErrorDefault = 0.05f;
	// default performance monitor, changed with MONITOR SET PERFORMANCE ERROR
	ErrorDefaultPerformance = 0.2f;
}

/*======================================================================*/
/*ParseMonitorResults parse the monitor command */
void ParseMonitorResults(Parser &p)
{
	long nelem,
		n2;

	DEBUG_ENTRY( "ParseMonitorResults()" );

	if( !lgInitDone )
	{
		fprintf( ioQQQ, " ParseMonitorResults called before InitAsserResults\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* has space been allocated yet? */
	if( !lgSpaceAllocated )
	{
		/* - no, we must allocate it */
		/* remember that space has been allocated */
		lgSpaceAllocated = true;

		/* create space for the array of labels*/
		chAssertLineLabel.resize(NASSERTS);
		lineids.resize(NASSERTS);

		assertBlends.resize(NASSERTS);

		/* the 2-character string saying what type of monitor */
		chAssertType.resize(NASSERTS);

		/* these are a pair of optional parameters */
		Param.reserve(NASSERTS);

		/* now fill out the 2D arrays */
		for( long i=0; i < NASSERTS; ++i )
		{
			chAssertLineLabel[i] = "unkn" ;

			assertBlends[i] = NULL;

			chAssertType[i] = "un";

			/* these are a pair of optional parameters */
			Param.reserve( i, 5 );
		}

		Param.alloc();

		/* now make space for the asserted quantities  */
		AssertQuantity.resize(NASSERTS);

		AssertQuantity2.resize(NASSERTS);

		/* now the errors */
		AssertError.resize(NASSERTS);

		/* now the line wavelengths */
		wavelength.resize(NASSERTS);

		/* label for a species */
		strAssertSpecies.resize(NASSERTS);

		/* now the flag saying whether should be log */
		lgQuantityLog.resize(NASSERTS);

		/* the flag for upper, lower limit, or equal */
		chAssertLimit.resize(NASSERTS);

		iLineType.resize(NASSERTS, -1);
	}
	/* end space allocation - we are ready to roll */

	/* read asserted values from file if GRID keyword is present */
	vector<double> AssertVector;
	if( p.nMatch( "GRID" ) )
	{
		ASSERT( optimize.nOptimiz >= 0 );

		/* this file should contain all the values that are asserted */
		string chLabel;
		if (p.GetQuote(chLabel))
			p.StringError();

		input_readvector(chLabel, AssertVector);
		if( AssertVector.size() < size_t(optimize.nOptimiz+1) )
		{
			fprintf(ioQQQ,"PROBLEM the file %s does not have enough values. Padding with zeroes.\n", chLabel.c_str() );
			while( AssertVector.size() < size_t(optimize.nOptimiz+1) )
				AssertVector.emplace_back(0.);
		}
	}

	/* false means print linear quantity - will be set true if entered
	 * quantity comes in as a log */
	lgQuantityLog[nAsserts] = false;

	/* all asserts have option for quantity to be a limit, or the quantity itself */
	if( p.nMatch("<" ) )
	{
		chAssertLimit[nAsserts] = '<';
	}
	else if( p.nMatch(">" ) )
	{
		chAssertLimit[nAsserts] = '>';
	}
	else
	{
		chAssertLimit[nAsserts] = '=';
	}

	/* which quantity will we check?, first is */

	if( p.nMatch(" SET" ) )
	{
		/* set an option for the monitor command, not an actual asserted
		 * quantity */

		/* decrement number of asserts since will be incremented below,
		 * this is not an actual asserted quantity */
		if( nAsserts >0 )
			fprintf(ioQQQ," The default monitor error is being changed after"
			" some asserts were entered.  \n This will only affect asserts"
			" that come after this command.\n");
		--nAsserts;

		if( p.nMatch("ERRO" ) )
		{
			if( p.nMatch("PERF" ) )
			{
				// set performance monitor error
				ErrorDefaultPerformance =
					(realnum)p.FFmtRead();
				if( p.lgEOL() )
					p.NoNumb("error");
			}
			/* set default error */
			ErrorDefault = 
				(realnum)p.FFmtRead();
			if( p.lgEOL() )
				p.NoNumb("error");
		}
		else
		{
			/* problem - no recognized quantity */
			fprintf( ioQQQ, 
				" I could not identify an option on this ASSERT SET XXX command.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* monitor mean ionization */
	else if( p.nMatch("IONI" ) )
	{

		/* say that this will be mean ionization fraction */

		/* f will indicate average over radius, F over volume -
		 * check whether keyword radius or volume occurs,
		 * default will be radius */
		if( p.nMatch("VOLU" ) )
		{
			chAssertType[nAsserts] = "F ";
		}
		else
		{
			/* this is default case, Fraction over radius */
			chAssertType[nAsserts] = "f ";
		}

		/* first get element label and make null terminated string*/
		if( (nelem = p.GetElem()) < 0 )
		{
			fprintf( ioQQQ, 
				" I could not identify an element name on this line.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		ASSERT( nelem>= 0);
		ASSERT( nAsserts>= 0);
		/* we now have element name, copy 4-char string (elementnames.chElementNameShort[nelem])
		 * into array that will be used to get ionization after calculation */
		chAssertLineLabel[nAsserts] = elementnames.chElementNameShort[nelem];

		/* now get ionization stage, which will be saved into wavelength */
		wavelength[nAsserts] = (realnum)p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("ionization stage");
		}
		/* ionization stage must be 1 or greater, but not greater than nelem (c scale)+2 */
		if( wavelength[nAsserts] < 1 || wavelength[nAsserts] > nelem+2 )
		{
			fprintf( ioQQQ, 
				"  The ionization stage is inappropriate for this element.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( p.nMatch( "GRID" ) )
		{
			AssertQuantity[nAsserts] = AssertVector[optimize.nOptimiz];
		}
		else
		{
			/* now get ionization fraction, log if number is negative or == 0, 
			 * linear if positive but less than or equal to 1.*/
			AssertQuantity[nAsserts] = p.FFmtRead();
			if( p.lgEOL() )
				p.NoNumb("ionization fraction");
		}

		/* optional error, default available */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
			AssertError[nAsserts] = ErrorDefault;

		/* now make sure we end up with positive linear ionization fraction that
		 * is greater than 0 but less than or equal to 1. */
		if( AssertQuantity[nAsserts] <= 0. )
		{
			/* log since negative or 0 */
			AssertQuantity[nAsserts] = 
				exp10(AssertQuantity[nAsserts] );
			/* entered as a log, so print as a log too */
			lgQuantityLog[nAsserts] = true;
		}
		else if( AssertQuantity[nAsserts] > 1. )
		{
			fprintf( ioQQQ, 
				"  The ionization fraction must be less than one.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* result cannot be zero */
		if( fabs(AssertQuantity[nAsserts]) <= SMALLDOUBLE )
		{
			fprintf( ioQQQ, 
				"  The ionization ionization fraction is too small, or zero.  Check input\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* molecular fraction averaged over model */
	else if( p.nMatch("MOLE" )&&p.nMatch("FRAC" ) )
	{

		/* say that this will be mean molecular fraction */

		/* mf will indicate average over radius, MF over vol -
		 * check whether keyword radius or volume occurs,
		 * default will be radius */
		/** \todo	2	NB this is not used, should do both, and more molecules (H2 only for now) */
		if( p.nMatch("VOLU" ) )
		{
			chAssertType[nAsserts] = "MF";
		}
		else
		{
			/* this is default case, Fraction over radius */
			chAssertType[nAsserts] = "mf";
		}

		if( p.nMatchErase(" H2 " ) )
		{
			chAssertLineLabel[nAsserts] = "H2" ;
			/* increment to get past the label */
		}
		else if( p.nMatch(" CO " ) )
		{
			chAssertLineLabel[nAsserts] = "CO" ;
		}
		else
		{
			fprintf( ioQQQ, 
				" I could not identify CO or H2 on this line.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* not meaningful */
		wavelength[nAsserts] = 0;

		/* now get log of molecular fraction */
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("molecular fraction");
		}
		if( AssertQuantity[nAsserts] <= 0. )
		{
			/* if negative then entered as log, but we will compare with linear */
			AssertQuantity[nAsserts] = 
				exp10(AssertQuantity[nAsserts] );
		}

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
		/* print results as logs */
		lgQuantityLog[nAsserts] = true;
	}

	/* monitor line "LINE"  --  key is ine_ since linear option appears on some commands */
	else if( p.nMatch(" LINE " ) )
	{
		if(  p.nMatch("LUMI") || p.nMatch("INTE"))
		{
			/* say that this is a line luminosity or intensity*/
			chAssertType[nAsserts] = "Ll";

			/* entered as a log, so print as a log too */
			lgQuantityLog[nAsserts] = true;
		}
		else if( p.nMatch( "CASE" ) && p.nMatch( " B ") )
		{
			/* say that this is a case B ratio */
			chAssertType[nAsserts] = "Lb";

			/* entered linear quantity, so print as linear too */
			lgQuantityLog[nAsserts] = false;
		}
		else
		{
			/* say that this is line relative to norm line - this is the default */
			chAssertType[nAsserts] = "Lr";

			/* entered linear quantity, so print as linear too */
			lgQuantityLog[nAsserts] = false;
		}

		iLineType[nAsserts] = 0;

		if( p.nMatch( "EMER" ) )
		{
			iLineType[nAsserts] = 1;
		}

		if( p.nMatch( "CUMU" ) )
		{
			iLineType[nAsserts] += 2;
		}

		/* this will check a line intensity, get the line id */
		LineID line = p.getLineID(false);

		blend_iterator b = blends.find(line.chLabel);
		if( 0 && b != blends.end() )
		{
			assertBlends[nAsserts] = &(b->second);
		}
		else if( line.chLabel.length() > NCHLAB-1 )
		{
			fprintf( ioQQQ, " The label must be no more than %d char long, between double quotes.\n",
				NCHLAB-1 );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* store line into array */
		lineids[nAsserts] = line;
		// need these as well because of PrtOneMonitor()
		chAssertLineLabel[nAsserts] = line.chLabel;
		wavelength[nAsserts] = line.wave;

		/* now get intensity or luminosity - 
		 * rel intensity is linear and intensity or luminosity are log */
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("intensity/luminosity");
		}
		/* luminosity was entered as a log */
		if( lgQuantityLog[nAsserts] )
		{
			if( AssertQuantity[nAsserts] > DBL_MAX_10_EXP || 
				AssertQuantity[nAsserts] < -DBL_MAX_10_EXP )
			{
				fprintf( ioQQQ, 
					" The asserted quantity is a log, but is too large or small, value is %e.\n",
					 AssertQuantity[nAsserts] );
				fprintf( ioQQQ, " I would crash if I tried to use it.\n Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			AssertQuantity[nAsserts] = 
				exp10(AssertQuantity[nAsserts] );
		}
		if( AssertQuantity[nAsserts]< 0. )
		{
			fprintf( ioQQQ, 
				" The relative intensity must be non-negative, and was %e.\n",AssertQuantity[nAsserts] );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* optional error, default available */
		AssertError[nAsserts] =	p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
	}

	/* monitor line predictions relative to case B */
	else if( p.nMatch("CASE" ) )
	{
		/* this is Case B for some element */
		chAssertType[nAsserts] = "CB";
		/* this is relative error */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		AssertQuantity[nAsserts] = 0;
		wavelength[nAsserts] = 0.;

		/* faint option - do not test line if relative intensity is less
		 * than entered value */
		if( p.GetParam("FAINT",&Param[nAsserts][4]) )
		{
			if( p.lgEOL() ) 
			{
				/* did not get 2 numbers */
				fprintf(ioQQQ," The monitor Case B faint option must have a number,"
					" the relative intensity of the fainest line to monitor.\n");
				cdEXIT(EXIT_FAILURE);
			}
			/* number is log if <= 0 */
			if( Param[nAsserts][4]<=0. )
				Param[nAsserts][4] = exp10( Param[nAsserts][4] );
		}
		else
		{
			/* use default - include everything*/
			Param[nAsserts][4] = SMALLFLOAT;
		}

		/* range option - to limit check on a certain wavelength range */
		if( p.GetRange("RANG",&Param[nAsserts][2],&Param[nAsserts][3]) )
		{
			if( p.lgEOL() ) 
			{
				/* did not get 2 numbers */
				fprintf(ioQQQ," The monitor Case B range option must have two numbers,"
					" the lower and upper limit to the wavelengths in Angstroms.\n");
				fprintf(ioQQQ," There must be a total of three numbers on the line,"
					" the relative error followed by the lower and upper limits to the "
					"wavelength in Angstroms.\n");
				cdEXIT(EXIT_FAILURE);
			}
			if( Param[nAsserts][2]>Param[nAsserts][3])
			{
				/* make sure in increasing order */
				double sav = Param[nAsserts][3];
				Param[nAsserts][3] = Param[nAsserts][2];
				Param[nAsserts][2] = sav;
			}
		}
		else
		{
			/* use default - include everything*/
			Param[nAsserts][2] = 0.;
			Param[nAsserts][3] = 1e30;
		}
		/* monitor case b for H - O checking against Hummer & Storey tables */
		if( p.nMatch("H-LI" ) )
		{
			/* H-like - now get an element */
			if( (nelem = p.GetElem()) < 0 )
			{
				/* no name found */
				fprintf(ioQQQ, "monitor H-like case B did not find an element on this line, sorry\n");
				p.PrintLine(ioQQQ);
				cdEXIT(EXIT_FAILURE);
			}
			if( nelem>7 )
			{
				/* beyond reach of tables */
				fprintf(ioQQQ, "monitor H-like cannot do elements heavier than O, sorry\n");
				p.PrintLine(ioQQQ);
				cdEXIT(EXIT_FAILURE);
			}
			Param[nAsserts][0] = ipH_LIKE;
			Param[nAsserts][1] = nelem;
			/* generate string to find simple prediction, as in "Ca B" */
			chAssertLineLabel[nAsserts] = "Ca ";
			if( p.nMatch("CASE A " ) )
				chAssertLineLabel[nAsserts] += "A";
			else
				chAssertLineLabel[nAsserts] += "B";
		}
		else if( p.nMatch("HE-L") )
		{
			/* He-like - only helium itself */
			Param[nAsserts][0] = ipHE_LIKE;
			Param[nAsserts][1] = ipHELIUM;

			chAssertLineLabel[nAsserts] = "Ca B";
		}
		else
		{
			/*no option found */
			fprintf( ioQQQ, 
				" I could not identify an iso-sequence on this Case A/B command.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* monitor departure coefficients */
	else if( p.nMatch("DEPA") )
	{
		string chLabel;
		if ( p.GetQuote(chLabel) == 0 )
		{
			strAssertSpecies[nAsserts] = chLabel;
			// printed label must still be less than NCHLAB in length
			ASSERT(chLabel.length() <= NCHLAB-1); 
			// Pad with spaces
			chAssertLineLabel[nAsserts] = chLabel;
			for (long i=chLabel.length();i<NCHLAB-1;++i)
				chAssertLineLabel[nAsserts] += ' ';
		}

		/* get expected average departure coefficient, almost certainly 1 */
		AssertQuantity[nAsserts] =	p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("average departure coefficient");

		/* this is relative error, max departure from unity of any level or std */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;

		if( !strAssertSpecies[nAsserts].empty() )
		{
			// remember this is departure coefficient for some species 
			chAssertType[nAsserts] = "DC";
		}
		/* H-like key means do one of the hydrogenic ions */
		else if( p.nMatch("H-LI" ) )
		{
			Param[nAsserts][0] = ipH_LIKE;
			chAssertLineLabel[nAsserts] = "dHyd" ;
			/* remember this is departure coefficient for some element */
			chAssertType[nAsserts] = "DI";
			/* now get element number for h ion from element name on card */
			if( (wavelength[nAsserts] = (realnum)p.GetElem()) < 0 )
			{
				fprintf( ioQQQ, 
					" I could not identify an element name on this line.\n");
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			if( p.nMatch("ZEROOK") )
				Param[nAsserts][1] = 1.;
			else
				Param[nAsserts][1] = 0.;
		}

		/* He-like key means do one of the helike ions */
		else if( p.nMatch("HE-L" ) )
		{
			Param[nAsserts][0] = ipHE_LIKE;
			chAssertLineLabel[nAsserts] = "dHel" ;
			/* remember this is departure coefficient for some element */
			chAssertType[nAsserts] = "DI";
			/* now get element number for h ion from element name on card */
			if( (wavelength[nAsserts] = (realnum)p.GetElem()) < 0 )
			{
				fprintf( ioQQQ, 
					" I could not identify an element name on this line.\n");
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			if( p.nMatch("ZEROOK") )
				Param[nAsserts][1] = 1.;
			else
				Param[nAsserts][1] = 0.;
		}

		/* this is H- h minus */
		else if( p.nMatch("HMIN" ) )
		{
			/* label */
			chAssertLineLabel[nAsserts] = "d H-" ;
			/* remember this is departure coefficient for H- */
			chAssertType[nAsserts] = "d-";
			/* the wavelength is meaningless */
			wavelength[nAsserts] = -1;
		}
		else
		{
			fprintf( ioQQQ, 
				" There must be a second key: H-LIke, HMINus, or HE-Like followed by element.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* last check for key "excited" - which means to skip the ground state */
		if( p.nMatch("EXCI" ) )
		{
			/* this is lowest level - do not do 0 */
			AssertQuantity2[nAsserts] = 1.;
		}
		else
		{
			/* do the ground state */
			AssertQuantity2[nAsserts] = 0.;
		}
	}

	/* monitor some results from map */
	else if( p.nMatch(" MAP" ) )
	{

		/* must have heating or cooling, since will check one or the other */
		/* check heating cooling results from map at some temperature */
		if( p.nMatch("HEAT" ) )
		{
			chAssertType[nAsserts] = "mh";
		}
		else if( p.nMatch("COOL" ) )
		{
			chAssertType[nAsserts] = "mc";
		}
		else
		{
			fprintf( ioQQQ, 
				" There must be a second key, HEATing or COOLing.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* now get temperature for AssertQuantity2 array*/
		AssertQuantity2[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("temperature");
		}

		if( AssertQuantity2[nAsserts] <= 10. )
		{
			/* entered as log, but we will compare with linear */
			AssertQuantity2[nAsserts] = 
				exp10(AssertQuantity2[nAsserts] );
		}

		/* print the temperature in the wavelength column */
		wavelength[nAsserts] = (realnum)AssertQuantity2[nAsserts];

		/* heating or cooling, both log, put into error */
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("heating/cooling");
		}

		/* AssertQuantity array will have heating or cooling */
		AssertQuantity[nAsserts] = exp10( AssertQuantity[nAsserts]);

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}

		/* entered as a log, so print as a log too */
		lgQuantityLog[nAsserts] = true;
	}

	/* monitor column density of something */
	else if( p.nMatch("COLU" ) )
	{
		/* this is column density */
		chAssertType[nAsserts] = "cd";

		/* this says to look for molecular column density, also could be ion stage */
		wavelength[nAsserts] = 0;

		string chLabel;
		
		if ( p.GetQuote(chLabel) == 0 )
		{
			// cdColm not yet able to cope with longer species
			ASSERT(chLabel.length() <= NCHLAB-1); 
			trimTrailingWhiteSpace( chLabel );
			chAssertLineLabel[nAsserts] = chLabel;
			for (long i=chLabel.length();i<NCHLAB-1;++i)
				chAssertLineLabel[nAsserts] += ' ';
		}
		else if( p.nMatchErase(" H2 ") )
		{
			chAssertLineLabel[nAsserts] = "H2" ;
		}
		else if( p.nMatchErase("H3+ "))
		{
			chAssertLineLabel[nAsserts] = "H3+" ;
		}
		else if( p.nMatchErase("H2+ "))
		{
			chAssertLineLabel[nAsserts] = "H2+" ;
		}
		else if( p.nMatchErase(" H- "))
		{
			chAssertLineLabel[nAsserts] = "H-" ;
		}
		else if( p.nMatchErase("H2G "))
		{
			chAssertLineLabel[nAsserts] = "H2g" ;
		}
		else if( p.nMatchErase("H2* "))
		{
			chAssertLineLabel[nAsserts] = "H2*" ;
		}
		else if( p.nMatchErase("HEH+"))
		{
			chAssertLineLabel[nAsserts] = "HeH+" ;
		}
		else if( p.nMatchErase(" O2 "))
		{
			chAssertLineLabel[nAsserts] = "O2" ;
		}
		else if( p.nMatchErase("H2O "))
		{
			chAssertLineLabel[nAsserts] = "H2O"  ;
		}
		else if( p.nMatchErase(" C2 "))
		{
			chAssertLineLabel[nAsserts] = "C2" ;
		}
		else if( p.nMatchErase(" C3 "))
		{
			chAssertLineLabel[nAsserts] = "C3" ;
		}
		else if( p.nMatch(" CO "))
		{
			chAssertLineLabel[nAsserts] = "CO" ;
		}
		else if( p.nMatch("SIO "))
		{
			chAssertLineLabel[nAsserts] = "SiO" ;
		}
		else if( p.nMatch(" OH "))
		{
			chAssertLineLabel[nAsserts] = "OH" ;
		}
		else if( p.nMatch(" CN ") )
		{
			chAssertLineLabel[nAsserts] = "CN";
		}
		else if( p.nMatch(" CH ") )
		{
			chAssertLineLabel[nAsserts] = "CH" ;
		}
		else if( p.nMatch(" CH+") )
		{
			chAssertLineLabel[nAsserts] = "CH+" ;
		}
		else
		{
			fprintf( ioQQQ, 
				" I could not identify H2, H3+, H2+, H2g, H2*, H2H+, CO, O2, SiO, OH, C2 or C3 or on this line.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		
		if( p.nMatch( "LEVE" ) )
		{
			if (chAssertLineLabel[nAsserts] != "H2")
			{
				fprintf( ioQQQ, " LEVEL option only implemented for H2.\n");
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			/* this is option for level-specific column density,
			 * next two numbers must be v then J */
			Param[nAsserts][0] = p.FFmtRead();
			if( p.lgEOL() )
				p.NoNumb("level v" );
			Param[nAsserts][1] = p.FFmtRead();
			if( p.lgEOL() )
				p.NoNumb("level J" );
			/* wavelength will be 10. * vib + rot */
			wavelength[nAsserts] = (realnum)(100.*Param[nAsserts][0] + Param[nAsserts][1]);
		}
		else
		{
			/* these are flags saying not to do state specific column densities */
			Param[nAsserts][0] = -1.;
			Param[nAsserts][1] = -1.;
		}
		
		/* i was set above for start of scan */
		/* now get log of column density */
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("column density");
		}
		/* entered as log, but we will compare with linear */
		AssertQuantity[nAsserts] = 
			exp10(AssertQuantity[nAsserts] );

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
		/* the keyword log is special for this case, since H2 and CO column densities can
		 * be so very unstable.  look for work log, in which case the error is log not linear.
		 * main column is always a log */
		if( p.nMatch( " LOG" ) )
		{
			AssertError[nAsserts] = exp10( AssertError[nAsserts] );
		}

		/* entered as a log, so print as a log too although asserted quantity is now linear */
		lgQuantityLog[nAsserts] = true;
	}

	/* monitor rate H2 forms on grain surfaces */
	else if( p.nMatch("GRAI") && p.nMatch(" H2 ") )
	{
		/* this flag will mean h2 form on grains */
		chAssertType[nAsserts] = "g2";
		/* a label */
		chAssertLineLabel[nAsserts] = "R H2" ;
		/* now get the first number on the line, which must be the 2 in H2 */
		nelem = (long int)p.FFmtRead();
		if( nelem!=2 )
		{
			fprintf( ioQQQ, 
				" I did not find a 2, as in H2, as the first number on this line.\n");
			fprintf( ioQQQ, 
				" The rate should be the second number.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		/* now past the 2 in h2, get the real number */
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("grain H2 formation");
		}
		/* if negative (almost certainly the case) then the log of the rate coefficient */
		if( AssertQuantity[nAsserts] <=0. )
			AssertQuantity[nAsserts] = exp10(AssertQuantity[nAsserts]  );
		/* will not use this */
		wavelength[nAsserts] = 0;

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
			/* want to print as a log since so small */
			lgQuantityLog[nAsserts] = true;
		}
	}

	/* monitor grain potential */
	else if( p.nMatch( "GRAI" ) && p.nMatch( "POTE") )
	{
		/* this flag will mean grain potential */
		chAssertType[nAsserts] = "gp";
		/* a label */
		chAssertLineLabel[nAsserts] = "GPot" ;
		/* now get the first number on the line */
		/* grain bin number */
		wavelength[nAsserts] = (realnum)p.FFmtRead();
		/* the potential itself, in volt, always linear */
		AssertQuantity[nAsserts] = p.FFmtRead();

		if( p.lgEOL() )
		{
			p.NoNumb("grain potential");
		}

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
			
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
	}

	/* monitor mean temperature, monitor temperature hydrogen 2 8000 */
	else if( p.nMatch("TEMP") )
	{
		/* say that this will be mean temperature, electron or grain */

		/* t will indicate temperature average over radius, T over volume -
		 * check whether keyword radius or volume occurs,
		 * default will be radius */
		if( p.nMatch("VOLU") )
		{
			chAssertType[nAsserts] = "T ";
		}
		else
		{
			/* this is default case, Fraction over radius */
			chAssertType[nAsserts] = "t ";
		}

		/* first look for keyword Grains, since label silicate may be on it,
		 * and this would trigger the element search */
		if( p.nMatch("GRAI") )
		{
			chAssertLineLabel[nAsserts] = "GTem" ;
		}

		/* face is temperature at illuminated face of cloud */
		else if( p.nMatch("FACE") )
		{
			chAssertLineLabel[nAsserts] = "face" ;
		}

		else if( p.nMatch( "21CM" ) )
		{
			if( p.nMatch( "MEAN" ) )
			{
				chAssertLineLabel[nAsserts] = "21CM" ;
			}
			else if( p.nMatch( "SPIN" ) )
			{
				chAssertLineLabel[nAsserts] = "SPIN" ;
			}
			else if( p.nMatch( "OPTI" ) )
			{
				chAssertLineLabel[nAsserts] = "OPTI" ;
			}
			else
			{
				fprintf( ioQQQ,
					"No option given."
					" One of MEAN, SPIN, and OPTICAL may be used.\n" );
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* get species label -- last option */
		else if( !p.GetQuote( chAssertLineLabel[nAsserts] ) )
		{
			/* nothing to do here */ ;
		}

		else
		{
			fprintf( ioQQQ, 
				" I could not identify a valid option on this line.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( chAssertLineLabel[nAsserts] == "GTem" )
		{
			wavelength[nAsserts] = (realnum)p.FFmtRead();
			if( p.lgEOL() )
			{
				p.NoNumb("grain index");
			}
		
			if( wavelength[nAsserts] < 1 || wavelength[nAsserts] > (realnum)LONG_MAX )
			{
				fprintf( ioQQQ, "  Unacceptable grain index.\n");
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			wavelength[nAsserts] = 0;
		}

		if( p.nMatch( "GRID" ) )
		{
			AssertQuantity[nAsserts] = AssertVector[optimize.nOptimiz];
		}
		else
		{
			/* now get temperature, log if number is <= 10, else linear */
			AssertQuantity[nAsserts] = p.FFmtRead();
			if( p.lgEOL() )
				p.NoNumb("temperature");

			/* 21 is picked up by previous read; repeat */
			if( chAssertLineLabel[nAsserts] == "21CM" )
			{
				AssertQuantity[nAsserts] = p.FFmtRead();
				if( p.lgEOL() )
					p.NoNumb("temperature");
			}
		}

		/* optional error, default available */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
			AssertError[nAsserts] = ErrorDefault;

		/* now make sure we end up with positive linear temperature
		 * number is log if <=10 unless linear keyword appears */
		if( AssertQuantity[nAsserts] <= 10. && !p.nMatch( "LINE" ) )
		{
			/* log since negative or 0 */
			AssertQuantity[nAsserts] = 
				exp10(AssertQuantity[nAsserts] );
			/* entered as a log, so print as a log too */
			lgQuantityLog[nAsserts] = true;
		}
	}

	/* monitor log of helium hydrogen ionization correction factor */
	else if( p.nMatch("HHEI") )
	{
		/* this flag will mean H-He icf */
		chAssertType[nAsserts] = "hh";
		/* say that this is zone numbering */
		chAssertLineLabel[nAsserts] = "HHei" ;

		/* now get the ionization correction factor, it is always the linear
		 * quantity itself, since can be positive or negative*/
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("ionization correction factor");
		}
		/* will not use this */
		wavelength[nAsserts] = 0;

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
	}

	/* this large H2 molecule */
	else if( p.nMatch(" H2 ") )
	{
		/* ortho to para ratio for last computed zone */
		if( p.nMatch("ORTH") )
		{
			/* this flag will mean ortho to para ratio */
			chAssertType[nAsserts] = "or";
			/* say that this is ortho to para density ratio */
			chAssertLineLabel[nAsserts] = "orth" ;
		}
		else
		{
			fprintf( ioQQQ, 
				" I could not identify a second keyword on this line.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* now get the first number, which better be the 2 in H2 */
		n2 = (long)p.FFmtRead();
		if( p.lgEOL() || n2 != 2 )
		{
			p.NoNumb("the 2 in H2 ?!");
		}
		AssertQuantity[nAsserts] = p.FFmtRead();
		/* will not use this */
		wavelength[nAsserts] = 0;

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
	}

	/* monitor we are running in MPI mode */
	else if( p.nMatch(" MPI") )
	{
		/* this flag will mean number of zones */
		chAssertType[nAsserts] = "mp";
		/* say that this is zone numbering */
		chAssertLineLabel[nAsserts] = "mpi" ;

		wavelength[nAsserts] = 0.;
		AssertQuantity[nAsserts] = 1.;
		AssertError[nAsserts] = ErrorDefault;
	}

	/* monitor number of zones */
	else if( p.nMatch("NZON") )
	{
		/* this flag will mean number of zones */
		chAssertType[nAsserts] = "z ";
		/* say that this is zone numbering */
		chAssertLineLabel[nAsserts] = "zone" ;

		/* now get number of zones */
		wavelength[nAsserts] = 0.;
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("zone number");

		/* optional error */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
			AssertError[nAsserts] = ErrorDefault;
	}

	/* monitor (probably upper limit to) error in pressure across model */
	else if( p.nMatch("PRES") && p.nMatch("ERRO") )
	{
		/* this flag indicates ratio of standard deviation to the mean pressure */
		chAssertType[nAsserts] = "pr";
		/* say that this is error in pressure */
		chAssertLineLabel[nAsserts] = "pres" ;

		/* now get the pressure error, which will be saved into wavelength
		 * in nearly all cases this is limit to error */
		wavelength[nAsserts] = 0;
		AssertQuantity[nAsserts] = (double)p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("pressure error");
		}
		else if( AssertQuantity[nAsserts] <= 0.)
		{
			/* number <= 0 is log of error */
			AssertQuantity[nAsserts] = exp10(AssertQuantity[nAsserts]);
		}

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] =	p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
	}

	else if( p.nMatch("PRADMAX") )
	{
		/* monitor pradmax - max ratio of rad to gas pressure */
		/* this flag indicates ratio of rad to gas pressure */
		chAssertType[nAsserts] = "RM";
		/* say that this is error in pressure */
		chAssertLineLabel[nAsserts] = "Prad" ;

		/* now get the pressure error, which will be saved into wavelength
		 * in nearly all cases this is limit to error */
		wavelength[nAsserts] = 0;
		AssertQuantity[nAsserts] = (double)p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("PRADMAX");
		}
		else if( AssertQuantity[nAsserts] <= 0.)
		{
			/* number <= 0 is log of error */
			AssertQuantity[nAsserts] = exp10(AssertQuantity[nAsserts]);
		}

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
	}

	/* monitor secondary ionization rate, csupra */
	else if( p.nMatch("CSUP") )
	{
		/* this flag will mean secondary ionization, entered as log */
		chAssertType[nAsserts] = "sc";
		/* say that this is sec ioniz */
		chAssertLineLabel[nAsserts] = "sion" ;

		/* now get rate, saved into monitor quantity */
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("secondary ionization rate");
		}
		/* entered as log, make linear */
		AssertQuantity[nAsserts] = exp10( AssertQuantity[nAsserts] );

		/* no wavelength */
		wavelength[nAsserts] = 0;

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}

		/* we want to print the log of eden, not linear value */
		lgQuantityLog[nAsserts] = true;
	}

	/* monitor heating rate, erg/cm3/s, htot */
	else if( p.nMatch("HTOT") )
	{
		/* this flag will mean heating, entered as log */
		chAssertType[nAsserts] = "ht";

		/* say that this is heating rate */
		chAssertLineLabel[nAsserts] = "htot" ;

		/* now get rate, saved into monitor quantity */
  		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("heating rate");
		}
		/* entered as log, make linear */
		AssertQuantity[nAsserts] = exp10( AssertQuantity[nAsserts] );

		/* no wavelength */
		wavelength[nAsserts] = 0;

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}

		/* we want to print the log of the heating, not linear value */
		lgQuantityLog[nAsserts] = true;
	}

	/* monitor cooling rate, erg/cm3/s, ctot */
	else if( p.nMatch("CTOT") )
	{
		/* this flag will mean cooling */
		chAssertType[nAsserts] = "ct";

		/* say that this is cooling rate */
		chAssertLineLabel[nAsserts] = "ctot";

		/* Look for GRID command */
		if( p.nMatch( "GRID" ) )
		{
			AssertQuantity[nAsserts] = AssertVector[optimize.nOptimiz];
		}
		else
		{
			AssertQuantity[nAsserts] = p.FFmtRead();
			/* now get rate, saved into monitor quantity */
			if( p.lgEOL() )
			{
				p.NoNumb("cooling rate");
			}
		}
		/* entered as log, make linear */
		AssertQuantity[nAsserts] = exp10( AssertQuantity[nAsserts] );

		/* no wavelength */
		wavelength[nAsserts] = 0;

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}

		/* we want to print the log of the heating, not linear value */
		lgQuantityLog[nAsserts] = true;
	}

	/* monitor number of iterations per zone, a test of convergence */
	else if( p.nMatch("ITRZ") )
	{
		/* this flag will mean number of iterations per zone */
		chAssertType[nAsserts] = "iz";
		/* say that this is iterations per zone  */
		chAssertLineLabel[nAsserts] = "itrz" ;

		/* now get quantity */
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("iterations per zone");
		}
		/* wavelength is meaningless */
		wavelength[nAsserts] = 0;

		/* optional error, default available */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
			AssertError[nAsserts] = ErrorDefaultPerformance;
	}

	/* monitor electron density of the last zone */
	else if( p.nMatch("EDEN") )
	{
		/* this flag will mean electron density of the last zone */
		chAssertType[nAsserts] = "e ";
		/* say that this is electron density */
		chAssertLineLabel[nAsserts] = "eden" ;

		/* now get electron density, which is a log */
		AssertQuantity[nAsserts] = 
			exp10( p.FFmtRead() );
		if( p.lgEOL() )
		{
			p.NoNumb(" electron density of the last zone");
		}

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
		wavelength[nAsserts] = 0;

		/* we want to print the log of eden, not linear value */
		lgQuantityLog[nAsserts] = true;
	}

	/* monitor energy density - Tu - of last zone */
	else if( p.nMatch(" TU ") )
	{
		/* this flag will mean energy density of the last zone */
		chAssertType[nAsserts] = "Tu";
		chAssertLineLabel[nAsserts] = "Tu" ;

		/* get temperature */
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("energy density of last zone");

		/* now make sure we end up with positive linear temperature
		 * number is log if <=10 unless linear keyword appears */
		lgQuantityLog[nAsserts] = false;
		if( AssertQuantity[nAsserts] <= 10. && !p.nMatch( "LINE" ) )
		{
			/* log since negative or 0 */
			AssertQuantity[nAsserts] = 
				exp10(AssertQuantity[nAsserts] );
			/* entered as a log, so print as a log too */
		}

		/* optional error, default available (cannot do before loop since we
		* do not know how many numbers are on line */
		AssertError[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		wavelength[nAsserts] = 0;
	}

	/* monitor thickness or depth of model */
	else if( p.nMatch("THIC") || p.nMatch("DEPT") )
	{
		/* this flag will mean thickness or depth */
		chAssertType[nAsserts] = "th";
		/* say that this is thickness */
		chAssertLineLabel[nAsserts] = "thic" ;

		/* now get thickness, which is a log */
		AssertQuantity[nAsserts] = exp10( p.FFmtRead() );
		if( p.lgEOL() )
		{
			p.NoNumb("thickness or depth of model");
		}

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] =	p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
		wavelength[nAsserts] = 0;

		/* we want to print the log of eden, not linear value */
		lgQuantityLog[nAsserts] = true;
	}

	/* monitor outer radius of model */
	else if( p.nMatch("RADI")  )
	{
		/* this flag will mean radius */
		chAssertType[nAsserts] = "ra";
		/* say that this is radius */
		chAssertLineLabel[nAsserts] = "radi" ;

		/* now get radius, which is a log */
		AssertQuantity[nAsserts] = exp10( p.FFmtRead() );
		if( p.lgEOL() )
		{
			p.NoNumb("outer radius");
		}

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] =	p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
		wavelength[nAsserts] = 0;

		/* we want to print the log of radius, not linear value */
		lgQuantityLog[nAsserts] = true;
	}

	/* monitor number of iterations */
	else if( p.nMatch("NITE") )
	{
		/* this flag will mean number of iterations */
		chAssertType[nAsserts] = "Z ";
		/* say that this is iteration numbering */
		chAssertLineLabel[nAsserts] = "iter" ;

		wavelength[nAsserts] = 0.f;

		/* now get number of iterations, which will be saved into wavelength */
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("number of iterations");
		}

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] =	p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
	}

	/* monitor terminal velocity, at end of calculation
	 * input in km/s and saved that way, even though code uses cm/s */
	else if( p.nMatch("VELO") )
	{
		/* this flag will mean velocity */
		chAssertType[nAsserts] = "Ve";
		/* say that this is velocity */
		chAssertLineLabel[nAsserts] = "vel" ;
		wavelength[nAsserts] = 0;

		/* now get velocity */
		AssertQuantity[nAsserts] = p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("terminal velocity");

		/* optional error, default available (cannot do before loop since we
		 * do not know how many numbers are on line */
		AssertError[nAsserts] =	p.FFmtRead();
		if( p.lgEOL() )
		{
			/* default error was set above */
			AssertError[nAsserts] = ErrorDefault;
		}
	}
	/* monitor nothing - a pacifier */
	else if( p.nMatch("NOTH") )
	{
		chAssertType[nAsserts] = "NO";
		chAssertLineLabel[nAsserts] = "noth" ;
		wavelength[nAsserts] = 0;
		AssertQuantity[nAsserts] = 0.;
		AssertError[nAsserts] = ErrorDefault;
	}
	else
	{
		/* did not recognize a command */
		fprintf( ioQQQ, 
					"  Unrecognized command.  The line image was\n");
		p.PrintLine(ioQQQ);
		fprintf( ioQQQ, 
			"  The options I know about are: ionization, line, departure coefficient, map, column, "
			"temperature, nzone, csupre, htot, itrz, eden, thickness, niter, \n");
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* increment number of asserts and confirm that limit not exceeded */
	++nAsserts;
	if( nAsserts >= NASSERTS )
	{
		fprintf(ioQQQ,
			" ParseMonitorResults: too many asserts, limit is NASSERT=%d\n",
			NASSERTS );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

inline double get_error_ratio ( double pred, double assert )
{
	return 1. - safe_div(pred,assert,1.);
}

/*============================================================================*/
/*lgCheckMonitors, checks asserts, last thing cloudy calls, returns true if all are ok, false if problems */
bool lgCheckMonitors(
	/* this is the file we will write this to, usually standard output, 
	 * but can be save */
	FILE * ioMONITOR )
{
	double PredQuan[NASSERTS] , RelError[NASSERTS];
	/* this will be true if the quantity was found, and false if not.  Used to prevent
	 * big botch flag when quantity not found (as in removed old he atom) */
	bool lgFound[NASSERTS];
	double relint , absint;
	bool lg1OK[NASSERTS];
	long i,j;
	/* This structure is for reporting another close match for asserts of line
	 * intensities only.  The zeroth, first, and second elements for each monitor are,
	 * respectively, the first, second, and third matches the code finds, if any.
	 * A negative number means no match.  A positive number indicates the pointer
	 * in the line stack of that match.  */
	/* use multi_arr here to prevent bogus array bounds violations being reported by pgCC */
	multi_arr<long,2> ipDisambiguate(NASSERTS,3);
	long lgDisambiguate = false;
	char chCaps[NCHLAB], chFind[NCHLAB];
	realnum errorwave;

	DEBUG_ENTRY( "lgCheckMonitors()" );

	/* this is a global variable in monitor_results.h, and can be checked by
	 * other routines to see if asserts are ok - (most runs will not use asserts) */
	lgMonitorsOK = true;

	/* will be used if there were big botched monitors */
	lgBigBotch = false;

	/* the optimize*.in and stars_oppim*.in tests in the test suite include
	 * asserts while optimizing.  We do not want to  test the asserts during
	 * the optimization process, since we will find blown asserts and report
	 * major problems.  We do want to test asserts on the final model however.
	 * Printout will usually be turned off in all except the final model,
	 * so do not to the monitor tests if we are optimizing but not printing */
	if( !called.lgTalk && optimize.lgOptimize )
	{
		/* just return */
		return true;
	}

	/*fprintf(ioQQQ , "DEBUG grid %li\n", optimize.nOptimiz );*/

	/* this will usually just return, but with table lines will check 
	 * existence of some lines */
	if( lines_table() )
	{
		lgBigBotch = true;
		lgMonitorsOK = false;
	}

	/* first check all asserts, there probably are none */
	for(i=0; i<nAsserts; ++i )
	{
		lg1OK[i] = true;
		PredQuan[i] = 0.;
		RelError[i] = 0.;
		for(j=0; j<3; ++j )
			ipDisambiguate[i][j] = -1;

		/* this flag is set false if we don't find the requested quantity */
		lgFound[i] = true;

		/* which type of monitor? */
		/* is it intensity? */
		if( chAssertType[i] == "Lr" )
		{
			/* this is an intensity, get the line, returns false if could not find it */
			ipDisambiguate[i][0] = cdLine( lineids[i], &relint, &absint, iLineType[i] );
			if( ipDisambiguate[i][0] <= 0 )
			{
				fprintf( ioMONITOR, " monitor error: lgCheckMonitors could not find line ");
				prt_line_err( ioMONITOR, lineids[i] );

				fprintf( ioMONITOR, 
					" monitor error: The \"save line labels\" command is a good way to get a list of line labels.\n\n");
				/* go to next line */
				lg1OK[i] = false;
				RelError[i] = 100000.;
				PredQuan[i] = 0;
				lgMonitorsOK = false;
				lgFound[i] = false;
				continue;
			}
			else
			{
				/********* LINE DISAMBIGUATION *************/
				/* Here we look for lines with same label and small wavelength
				 * differences so that we can disambiguate below */

				/* change chLabel to all caps */
				strcpy(chFind, lineids[i].chLabel.c_str());
				caps(chFind);

				/* get the error associated with specified significant figures */
				errorwave = WavlenErrorGet( lineids[i].wave, LineSave.sig_figs );

				/* go through rest of line stack to look for close matches */
				for( j=1; j < LineSave.nsum; j++ )
				{
					/* don't bother with this one, we've already identified it. */
					if( j==ipDisambiguate[i][0] )
						continue;

					/* change chLabel to all caps to be like input chALab */
					cap4(chCaps, LineSave.lines[j].chALab());

					/* look for wavelengths within 3 error bars.
					 * For example, for a line entered in Angstroms with
					 * four significant figures, the error bar is 0.5 Ang.  
					 * So here we will find any lines within 1.5 Angstroms
					 * of the 
					 * asserted wavelength.  check wavelength and chLabel for a match */
					if( fabs(LineSave.lines[j].wavelength()-lineids[i].wave) < 3.f*errorwave )
					{
						/* now see if labels agree */
						if( strcmp(chCaps,chFind) == 0 )
						{
							double relint1, absint1, current_error;

							cdLine_ip( j, &relint1, &absint1, iLineType[i] );
							
							if (AssertError[i] > 0.)
								current_error = fabs(1. - relint1/AssertQuantity[i]);
							else
								current_error = relint1 - AssertQuantity[i];

							if( current_error < 2.*fabs(AssertError[i]) ||
								current_error < 2.*fabs(RelError[i]) )
							{
								lgDisambiguate = true;
								/* if second match (element 1) is already set,
								 * this is third match (element 2).  Set and break out. */
								if( ipDisambiguate[i][1] > 0 )
								{
									ipDisambiguate[i][2] = j;
									break;
								}
								else
								{
									ipDisambiguate[i][1] = j;
								}
							}
						}
					}
				}
			}

			PredQuan[i] = relint;
			if (AssertError[i] > 0.0)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]-AssertQuantity[i] ;
		}

		else if( chAssertType[i] == "Lb" )
		{
			if( cdLine( lineids[i], &relint, &absint, iLineType[i] ) <= 0 )
			{
				fprintf( ioMONITOR, " monitor error: lgCheckMonitors could not find line ");
				prt_line_err( ioMONITOR, lineids[i] );

				fprintf( ioMONITOR, 
					" monitor error: The \"save line labels\" command is a good way to get a list of line labels.\n\n");
				/* go to next line */
				RelError[i] = 10000000.;
				PredQuan[i] = 0;
				lg1OK[i] = false;
				lgFound[i] = false;
				lgMonitorsOK = false;
				continue;
			}

			double relint_cb = 0.,
				absint_cb = 0.;
			if( cdLine( "Ca B", lineids[i].wave, &relint_cb, &absint_cb, iLineType[i] ) <= 0 )
			{
				fprintf( ioMONITOR, " monitor error: lgCheckMonitors could not find line ");
				prt_line_err( ioMONITOR, lineids[i] );

				fprintf( ioMONITOR, 
					" monitor error: The \"save line labels\" command is a good way to get a list of line labels.\n\n");
				/* go to next line */
				RelError[i] = 10000000.;
				PredQuan[i] = 0;
				lg1OK[i] = false;
				lgFound[i] = false;
				lgMonitorsOK = false;
				continue;
			}

			PredQuan[i] = absint / absint_cb;
			if (AssertError[i] > 0.0)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i] - AssertQuantity[i];
			//	printf("i=%ld\t %s\t %g\t pred = %g\t asse = %g\t err = %g\n",
		}

		/*this is line luminosity */
		else if( chAssertType[i] == "Ll" )
		{
			/* this is a luminosity, get the line, returns false if could not find it */
			if( cdLine( lineids[i], &relint, &absint, iLineType[i] ) <= 0 )
			//if (indice <= 0 )
			{
				fprintf( ioMONITOR, " monitor error: lgCheckMonitors could not find line ");
				prt_line_err( ioMONITOR, lineids[i] );

				fprintf( ioMONITOR, 
					" monitor error: The \"save line labels\" command is a good way to get a list of line labels.\n\n");
				/* go to next line */
				RelError[i] = 10000000.;
				PredQuan[i] = 0;
				lg1OK[i] = false;
				lgFound[i] = false;
				lgMonitorsOK = false;
				continue;
			}
			PredQuan[i] = absint;
			if (AssertError[i] > 0.0)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i] - AssertQuantity[i];
		}
		else if( chAssertType[i] == "hh" )
		{
			double hfrac , hefrac;
			/* get H ionization fraction, returns false if could not find it */
			if( cdIonFrac(
				/* four char string, null terminated, giving the element name */
				"hydr", 
				/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number */
				1, 
				/* will be fractional ionization */
				&hfrac, 
				/* how to weight the average, must be "VOLUME" or "RADIUS" */
				"VOLUME" ,
				/* do not want extra factor of density */
				false)  )
			{
				fprintf( ioMONITOR, 
					" monitor error: lgCheckMonitors could not find h ionization fraction \n");
				lg1OK[i] = false;
				RelError[i] = 0;
				PredQuan[i] = 0;
				lgFound[i] = false;
				lgMonitorsOK = false;
				continue;
			}
			if( cdIonFrac(
				/* four char string, null terminated, giving the element name */
				"heli", 
				/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number */
				1, 
				/* will be fractional ionization */
				&hefrac, 
				/* how to weight the average, must be "VOLUME" or "RADIUS" */
				"VOLUME" ,
				/* do not want extra factor of density */
				false)  )
			{
				fprintf( ioMONITOR, 
					" monitor error: lgCheckMonitors could not find h ionization fraction \n");
				lg1OK[i] = false;
				RelError[i] = 0;
				PredQuan[i] = 0;
				lgFound[i] = false;
				lgMonitorsOK = false;
				continue;
			}
			/* the helium hydrogen ionization correction factor */
			PredQuan[i] = hefrac-hfrac;
			/* two icf's in difference, no need to div by mean since already on scale with unity */
			RelError[i] = fabs(AssertQuantity[i] - (hefrac-hfrac) );
		}

		else if( chAssertType[i] == "mp" )
		{
			PredQuan[i] = cpu.i().lgMPI() ? 1. : 0.;
			/* use absolute error */
			RelError[i] = AssertQuantity[i] -  PredQuan[i];
		}

		else if( chAssertType[i] == "z " )
		{
			/* this is the number of zones */
			PredQuan[i] = (double)nzone;
			if( t_version::Inst().lgRelease )
				RelError[i] = ForcePass(chAssertLimit[i]);
			else
			{
				/* two integers in difference */
				if (AssertError[i] > 0.)
					RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
				else
					RelError[i] = PredQuan[i] - AssertQuantity[i];
			}
		}

		else if( chAssertType[i] == "or" )
		{
			/* ortho to para ratio for large H2 molecule in last zone */
			PredQuan[i] = h2.ortho_density / SDIV( h2.para_density );

			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i] - AssertQuantity[i];
		}

		else if( chAssertType[i] == "g2" )
		{
			/* check Jura rate, rate per vol that H2 forms on grain surfaces */
			PredQuan[i] = gv.rate_h2_form_grains_used_total;
			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i] - AssertQuantity[i];
		}

		else if( chAssertType[i] == "RM" )
		{
			/* check Jura rate, rate per vol that H2 forms on grain surfaces */
			PredQuan[i] = pressure.RadBetaMax;
			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i] - AssertQuantity[i];
		}

		else if( chAssertType[i] == "pr" )
		{
			/* standard deviation of the pressure */
			double sumx=0., sumx2=0., average;
			long int n;
			/* do sums to form standard deviation */
			for( n=0; n<nzone; n++ )
			{
				sumx += struc.pressure[n];
				sumx2 += POW2(struc.pressure[n]);
			}
			if( nzone>1 )
			{
				/* this is average */
				average = sumx/nzone;
				/* this is abs std */
				sumx = sqrt( (sumx2-POW2(sumx)/nzone)/(nzone-1) );
				/* save the relative std */
				PredQuan[i] = sumx / average;
			}
			else
			{
				PredQuan[i] = 0.;
			}

			// this is already relative error, do not need 1-ratio
			RelError[i] = PredQuan[i];
		}

		else if( chAssertType[i] == "iz" )
		{
			/* this is number of iterations per zone, a test of convergence properties */
			if( nzone > 0 )
				PredQuan[i] = (double)(conv.nTotalIoniz-conv.nTotalIoniz_start)/(double)(nzone);
			else
				/* something big so monitor will botch. */
				PredQuan[i] = 1e10;

			if( t_version::Inst().lgRelease )
				RelError[i] = ForcePass(chAssertLimit[i]);
			else
			{
				/* this is relative error */
				if (AssertError[i] > 0.)
					RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
				else
					RelError[i] = PredQuan[i]- AssertQuantity[i];
			}
		}

		else if( chAssertType[i] == "e " )
		{
			/* this is electron density of the last zone */
			PredQuan[i] = dense.eden;
			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		else if( chAssertType[i] == "Tu" )
		{
			/* this is radiation energy density of the last zone */
			PredQuan[i] = powpq((rfield.EnergyIncidCont+rfield.EnergyDiffCont)/
					    (4.*STEFAN_BOLTZ),1,4);
			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		else if( chAssertType[i] == "th" )
		{
			/* this is thickness */
			PredQuan[i] = radius.depth;
			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		else if( chAssertType[i] == "ra" )
		{
			/* this is thickness */
			PredQuan[i] = radius.Radius;
			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		else if( chAssertType[i] == "Ve" )
		{
			/* this is final velocity of wind in km/s (code uses cm/s) */
			PredQuan[i] = wind.windv/1e5;
			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		else if( chAssertType[i] == "NO" )
		{
			/* this is nothing */
			PredQuan[i] = 0;
			/* this is relative error */
			RelError[i] = 0.;
		}

		else if( chAssertType[i] == "sc" )
		{
			/* this is secondary ionization rate */
			PredQuan[i] = secondaries.csupra[ipHYDROGEN][0];
			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		else if( chAssertType[i] == "ht" )
		{
			/* this is heating rate */
			PredQuan[i] = thermal.htot;
			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		else if( chAssertType[i] == "ct" )
		{
			/* this is cooling rate */
			PredQuan[i] = thermal.ctot;
			/* this is relative error */
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		else if( chAssertType[i] == "Z " )
		{
			/* this is the number of iterations */
			PredQuan[i] = (double)iteration;
			if( t_version::Inst().lgRelease )
				RelError[i] = ForcePass(chAssertLimit[i]);
			else
			{
				/* two integers in difference */
				if (AssertError[i] > 0.)
					RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
				else
					RelError[i] = PredQuan[i]- AssertQuantity[i];
			}
		}

		else if( chAssertType[i] == "CB" )
		{
			long int nISOCaseB = (long)Param[i][0];
			long int nelemCaseB = (long)Param[i][1];
			string chElemLabelCaseB = chIonLbl( nelemCaseB+1, nelemCaseB+1-nISOCaseB );

			/* sets lowest quantum number index */
			int iCase;
			if( "Ca A" == chAssertLineLabel[i])
				iCase = 0;
			else if( "Ca B" == chAssertLineLabel[i] )
				iCase = 1;
			else
				TotalInsanity();

			iLineType[i] = 0;

			RelError[i] = 0.;
			long nHighestPrinted = iso_sp[nISOCaseB][nelemCaseB].n_HighestResolved_max;
			if( nISOCaseB == ipH_LIKE )
			{
				fprintf(ioMONITOR,"                 Species  nHi nLo  Wl  Computed  Asserted       error\n");
				/* limit of 10 is because that is all we printed and saved in prt_lines_hydro 
				 * wavelengths will come out of atmdat.WaveLengthCaseB - first index is
				 * nelem on C scale, H is 0, second two are configurations of line on
				 * physics scale, so Ha is 3-2 */
				for( long int ipLo=1+iCase; ipLo< MIN2(10,nHighestPrinted-1); ++ipLo )
				{
					for( long int ipHi=ipLo+1; ipHi< MIN2(25,nHighestPrinted); ++ipHi )
					{
						/* monitor the line */
						realnum wl = atmdat.WaveLengthCaseB[nelemCaseB][ipHi][ipLo];
						/* range option to restrict wavelength coverage */
						if( wl < Param[i][2] || wl > Param[i][3] )
							continue;

						double relint, absint, CBrelint, CBabsint;
						/* find the predicted line intensity */
						cdLine( chAssertLineLabel[i], wl, &CBrelint, &CBabsint, iLineType[i] );
						if( CBrelint < Param[i][4]  )
							continue;
						double error;
						/* now find the Case B intensity - may not all be present */
						if( cdLine( chElemLabelCaseB, wl, &relint, &absint, iLineType[i] ) > 0 )
						{
							if (AssertError[i] > 0.)
								error = (CBabsint - absint)/MAX2(CBabsint , absint);
							else
								error = (CBabsint - absint);
							double RelativeError = fabs(error / AssertError[i]);
							/* start of line, flag problems */
							if( RelativeError < 1. )
							{
								if( RelativeError < 0.25 )
								{
									fprintf( ioMONITOR, " ChkMonitor        ");
								}
								else if( RelativeError < 0.50 )
								{
									fprintf( ioMONITOR, " ChkMonitor -      ");
								}
								else if( RelativeError < 0.75 )
								{
									fprintf( ioMONITOR, " ChkMonitor --     ");
								}
								else if( RelativeError < 0.90 )
								{
									fprintf( ioMONITOR, " ChkMonitor ---    ");
								}
								else  if( RelativeError < 0.95 )
								{
									fprintf( ioMONITOR, " ChkMonitor ----   ");
								}
								else  if( RelativeError < 0.98 )
								{
									fprintf( ioMONITOR, " ChkMonitor -----  ");
								}
								else 
								{
									fprintf( ioMONITOR, " ChkMonitor ------ ");
								}

							}
							else
							{
								fprintf( ioMONITOR, " ChkMonitor botch>>");
							}
							fprintf(ioMONITOR," %s %3li %3li ", 
								chElemLabelCaseB.c_str() , ipHi , ipLo );
							prt_wl(ioMONITOR, wl );
							fprintf(ioMONITOR," %.2e %.2e %10.3f", 
								log10(absint) , log10(CBabsint) , error );
						}
						else
							TotalInsanity();
						if( fabs(error) > fabs(AssertError[i])  )
							fprintf(ioMONITOR , " botch \n");
						else
							fprintf(ioMONITOR , "\n");

						PredQuan[i] = 0;
						AssertQuantity[i] = 0.;
						RelError[i] = MAX2( RelError[i] , fabs(error) );

						/* save sum which we will report later */
						MonitorResults.SumErrorCaseMonitor += RelError[i];
						++MonitorResults.nSumErrorCaseMonitor;

					}
				}
				fprintf(ioMONITOR,"\n");
			}
			else if( nISOCaseB == ipHE_LIKE )
			{
				if( !dense.lgElmtOn[ipHELIUM] )
				{
					fprintf(ioQQQ,"PROBLEM monitor case B for a He is requested but He is not "
						"included.\n");
					fprintf(ioQQQ,"Do not turn off He if you want to monitor its spectrum.\n");
					cdEXIT(EXIT_FAILURE);
				}
				
				/* do He I as special case */
				fprintf(ioMONITOR,"                     Wl  Computed  Asserted       error\n");
				for( unsigned int ipLine=0; ipLine< atmdat.CaseBWlHeI.size(); ++ipLine )
				{
					/* monitor the line */
					realnum wl = atmdat.CaseBWlHeI[ipLine];
					/* range option to restrict wavelength coverage */
					if( wl < Param[i][2] || wl > Param[i][3] )
						continue;
					double relint , absint,CBrelint , CBabsint;
					cdLine( chAssertLineLabel[i], wl, &CBrelint, &CBabsint, iLineType[i] );
					if( CBrelint < Param[i][4]  )
						continue;
					double error;
					if( cdLine( chElemLabelCaseB, wl, &relint, &absint, iLineType[i] ) > 0)
					{
						if (AssertError[i] > 0.0)
							error = (CBabsint - absint)/MAX2(CBabsint , absint);
						else
							error = (CBabsint - absint);
						double RelativeError = fabs(error / AssertError[i]);
						/* start of line, flag problems */
						if( RelativeError < 1. )
						{
							if( RelativeError < 0.25 )
							{
								fprintf( ioMONITOR, " ChkMonitor        ");
							}
							else if( RelativeError < 0.50 )
							{
								fprintf( ioMONITOR, " ChkMonitor -      ");
							}
							else if( RelativeError < 0.75 )
							{
								fprintf( ioMONITOR, " ChkMonitor --     ");
							}
							else if( RelativeError < 0.90 )
							{
								fprintf( ioMONITOR, " ChkMonitor ---    ");
							}
							else  if( RelativeError < 0.95 )
							{
								fprintf( ioMONITOR, " ChkMonitor ----   ");
							}
							else  if( RelativeError < 0.98 )
							{
								fprintf( ioMONITOR, " ChkMonitor -----  ");
							}
							else 
							{
								fprintf( ioMONITOR, " ChkMonitor ------ ");
							}

						}
						else
						{
							fprintf( ioMONITOR, " ChkMonitor botch>>");
						}
						prt_wl(ioMONITOR, wl );
						fprintf(ioMONITOR," %.2e %.2e %10.3f", 
							absint , CBabsint , error );
					}
					else
						TotalInsanity();
					if( fabs(error) > fabs(AssertError[i])  )
						fprintf(ioMONITOR , " botch \n");
					else
						fprintf(ioMONITOR , "\n");

					PredQuan[i] = 0;
					AssertQuantity[i] = 0.;
					RelError[i] = MAX2( RelError[i] , fabs(error) );

					/* save sum which we will report later */
					MonitorResults.SumErrorCaseMonitor += RelError[i];
					++MonitorResults.nSumErrorCaseMonitor;
				}
				fprintf(ioMONITOR,"\n");
			}
			else
				TotalInsanity();
		}

		// departure coefficients for a species given in quotes 
		else if( chAssertType[i] == "DC" )
		{
			string this_species = strAssertSpecies[i];
			if( this_species.find( '[' ) == string::npos )
			{
				// keyword EXCITED appeared
				if( AssertQuantity2[i] == 1 )
					this_species += "[2:]";
				else 
					this_species += "[:]";
			}

			vector<long> speciesLevels;
			const molezone *sp = getLevelsGeneric( this_species, false, speciesLevels );
			if( sp == null_molezone )
			{
				fprintf( ioQQQ, "PROBLEM Could not find species between quotes: \"%s\".\n",
						strAssertSpecies[i].c_str() );
				cdEXIT(EXIT_FAILURE);
			}
	
			if( sp->levels == NULL )
			{
				fprintf( ioQQQ, "WARNING Species '%s' has no internal structure."
						" Cannot compute departure coefficient\n",
						strAssertSpecies[i].c_str() );
				// Just report unity?
				PredQuan[i] = 1.;
				RelError[i] = 0.;
			}
			else if( speciesLevels.size() > 0 )
			{
				ASSERT( sp->levels->size() > 0 );
				RelError[i] = 0.;
				PredQuan[i] = 0.;
	
				long numPrintLevels = 0;
				for( vector<long>::const_iterator ilvl = speciesLevels.begin(); ilvl != speciesLevels.end(); ilvl++ )
				{
					const qStateConstProxy& st = (*sp->levels)[ *ilvl ];
					if( st.status() == LEVEL_INACTIVE )
						continue;
					++numPrintLevels;
					PredQuan[i] += st.DepartCoef();
				}
				ASSERT( numPrintLevels > 0 );
				PredQuan[i] /= (double)(numPrintLevels);
				RelError[i] = AssertQuantity[i] - PredQuan[i];

#if 0
				if( 0 &&  fp_equal( Param[i][1], 1. ) && PredQuan[i]==0. )
				{
					// this should only happen if either the present stage or the parent has zero density.
					ASSERT( dense.xIonDense[nelem][nelem-ipISO]==0. || dense.xIonDense[nelem][nelem+1-ipISO]==0. );
					PredQuan[i] = AssertQuantity[i];
					RelError[i] = 0.;
				}
#endif
			}
			else
			{
				fprintf( ioQQQ, "WARNING Requested level(s) in '%s' do not exist\n",
						strAssertSpecies[i].c_str() );
				PredQuan[i] = 0.;
				RelError[i] = AssertQuantity[i]; 
			}
		}

		/* departure coefficients for something in isoelectronic sequences */
		else if( chAssertType[i] == "DI" )
		{
			/* this is departure coefficient for XX-like ion given by wavelength */
			/* stored number was element number on C scale */
			long ipISO = (long)Param[i][0];
			long nelem = (long)wavelength[i];
			if( !dense.lgElmtOn[nelem] )
			{
				fprintf(ioQQQ,"PROBLEM asserted element %ld is not turned on!\n",nelem);
				PredQuan[i] = 0.;
				RelError[i] = 0.;
			}
			else
			{
				RelError[i] = 0.;
				PredQuan[i] = 0.;
				long numPrintLevels = iso_sp[ipISO][nelem].numLevels_local - (long)AssertQuantity2[i];
				for( long n=(long)AssertQuantity2[i]; n<numPrintLevels+(long)AssertQuantity2[i]; ++n )
				{
					PredQuan[i] += iso_sp[ipISO][nelem].st[n].DepartCoef();
				}
				ASSERT( numPrintLevels > 0 );
				PredQuan[i] /= (double)(numPrintLevels);
				RelError[i] = AssertQuantity[i] - PredQuan[i];

				if( fp_equal( Param[i][1], 1. ) && PredQuan[i]==0. )
				{
					// this should only happen if either the present stage or the parent has zero density.
					ASSERT( dense.xIonDense[nelem][nelem-ipISO]==0. || dense.xIonDense[nelem][nelem+1-ipISO]==0. );
					PredQuan[i] = AssertQuantity[i];
					RelError[i] = 0.;
				}
			}
		}

		/* this is H- departure coefficient */
		else if( chAssertType[i] == "d-" )
		{
			PredQuan[i] = hmi.hmidep;
			RelError[i] = AssertQuantity[i] - hmi.hmidep;
		}

		/* this would be ionization fraction */
		else if( chAssertType[i] == "f " || chAssertType[i] == "F " )
		{
			char chWeight[7];
			if( chAssertType[i] == "F " )
			{
				strcpy( chWeight , "VOLUME" );
			}
			else
			{
				/* this is default case, Fraction over radius */
				strcpy( chWeight , "RADIUS" );
			}
			/* get ionization fraction, returns false if could not find it */
			if( cdIonFrac(
				/* four char string, null terminated, giving the element name */
					 chAssertLineLabel[i].c_str(), 
				/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number */
				(long)wavelength[i], 
				/* will be fractional ionization */
				&relint, 
				/* how to weight the average, must be "VOLUME" or "RADIUS" */
				chWeight ,
				/* do not want extra factor of density */
				false)  )
			{
				fprintf( ioMONITOR, 
					" monitor error: lgCheckMonitors could not find a line with label %s %f \n",
							chAssertLineLabel[i].c_str() , wavelength[i] );
				/* go to next line */
				lg1OK[i] = false;
				RelError[i] = 0;
				PredQuan[i] = 0;
				lgFound[i] = false;
				lgMonitorsOK = false;
				continue;
			}
			/* this is ionization fraction */
			PredQuan[i] = relint;
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		/* this would be column density of several molecules */
		else if( chAssertType[i] == "cd" ) 
		{
			/* for H2 column density - total or for a state? */
			if( ( chAssertLineLabel[i] == "H2" ) && (Param[i][0] >= 0.) )
			{
				/* this branch get state specific column density */
				/* get total H2 column density */
				if( (relint = cdH2_colden( (long)Param[i][0] , (long)Param[i][1] ) ) < 0. )
				{
					fprintf(ioQQQ," PROBLEM lgCheckMonitors did not find v=%li, J=%li for H2 column density.\n",
						(long)Param[i][0] , (long)Param[i][1] );
					lg1OK[i] = false;
					RelError[i] = 0;
					PredQuan[i] = 0;
					lgFound[i] = false;
					lgMonitorsOK = false;
					continue;
				}
			}
			else
			{
				/* get ionization fraction, returns 0 if all ok */
				if( cdColm(
					/* four char string, null terminated, giving the element name */
					chAssertLineLabel[i].c_str(),
					/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
					* zero for molecule*/
					(long)wavelength[i],
					/* will be fractional ionization */
					&relint) )
				{
					fprintf( ioMONITOR, 
						" monitor error: lgCheckMonitors could not find a molecule with label %s %f \n",
								chAssertLineLabel[i].c_str() , wavelength[i] );
					/* go to next line */
					lg1OK[i] = false;
					RelError[i] = 0;
					PredQuan[i] = 0;
					lgFound[i] = false;
					lgMonitorsOK = false;
					continue;
				}
			}
			/* this is ionization fraction */
			PredQuan[i] = relint;
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		/* this would be molecular fraction of CO or H2 */
		else if( chAssertType[i] == "MF" || chAssertType[i] == "mf" )
		{
			/* get molecular fraction, returns 0 if all ok */
			relint = 0.;
			if( chAssertLineLabel[i] == "H2" )
			{
				/* get total H2 column density */
				if( cdColm("H2  " , 0, 
					/* will be fractional ionization */
					&relint) )
					TotalInsanity();

				relint = relint / colden.colden[ipCOL_HTOT];
			}
			else
			{
				fprintf( ioMONITOR, 
					" monitor error: lgCheckMonitors could not find a molecule with label %s %f \n",
							chAssertLineLabel[i].c_str() , wavelength[i] );
				/* go to next line */
				lg1OK[i] = false;
				RelError[i] = 0;
				PredQuan[i] = 0;
				lgMonitorsOK = false;
				continue;
			}
			/* this is ionization fraction */
			PredQuan[i] = relint;
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		/* check heating/cooling at some temperature in a thermal map */
		else if( chAssertType[i] == "mh" || chAssertType[i] == "mc" ) 
		{
			/* check heating or cooling (stored in error array) at temperature in monitor results */
			/* check that map was done, and arrays have nmap elements */
			if( hcmap.nMapAlloc == 0 )
			{
				/* this happens if map not done and space for h/c not allocated */
				fprintf( ioMONITOR, 
					" monitor error: lgCheckMonitors cannot check map since map not done.\n");
				/* go to next line */
				lg1OK[i] = false;
				RelError[i] = 0;
				PredQuan[i] = 0;
				lgMonitorsOK = false;
				continue;
			}
			/* now check that requested temperature is within the range of the map we computed */
			if( AssertQuantity2[i]< hcmap.temap[0] || AssertQuantity2[i]> hcmap.temap[hcmap.nmap-1] )
			{
				fprintf( ioMONITOR, 
					" monitor error: lgCheckMonitors cannot check map since temperature not within range.\n");
				/* go to next line */
				lg1OK[i] = false;
				RelError[i] = 0;
				PredQuan[i] = 0;
				lgMonitorsOK = false;
				continue;
			}

			/* we should have valid data - find closest temperature >- requested temperature */
			j = 0;
			while( AssertQuantity2[i]>hcmap.temap[j]*1.001 && j < hcmap.nmap )
			{
				++j; 
			}

			/* j points to correct cell in heating cooling array */
			/* we will not interpolate, just use this value, and clobber te to prove it*/
			if( chAssertType[i] == "mh" )
			{
				/* heating */
				PredQuan[i] = hcmap.hmap[j];
				chAssertLineLabel[i] = "MapH" ;
			}
			else if( chAssertType[i] == "mc" )
			{
				/* cooling */
				PredQuan[i] = hcmap.cmap[j];
				chAssertLineLabel[i] = "MapC" ;
			}
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		/* this will be an average temperature */
		else if( chAssertType[i] == "t " || chAssertType[i] == "T " )
		{
			char chWeight[7];
			if( chAssertType[i] == "T " )
			{
				strcpy( chWeight , "VOLUME" );
			}
			else
			{
				/* this is default case, Fraction over radius */
				strcpy( chWeight , "RADIUS" );
			}

			/* options are average Te for ion, temp at ill face, or temp for grain */
			if( chAssertLineLabel[i] == "GTem" )
			{
				size_t nd;
				/* the minus one is because the grain types are counted from one,
				 * but stuffed into the c array, that counts from zero */
				nd = (size_t)wavelength[i]-1;
				if( nd >= gv.bin.size() )
				{
					fprintf( ioQQQ, "Illegal grain number found: %f\n" , wavelength[i] );
					fprintf( ioQQQ, "Use 1 for first grain that is turned on, " );
					fprintf( ioQQQ, "2 for second, etc....\n" );
					fprintf( ioQQQ, "Old style grain numbers are not valid anymore !!\n" );
					cdEXIT(EXIT_FAILURE);
				}
				relint = gv.bin[nd].avdust/radius.depth_x_fillfac;
			}
			else if( chAssertLineLabel[i] == "face" )
			{
				/* this is the temperature at the illuminated face */
				relint = struc.testr[0];
			}
			else
			{
				/* get temperature, returns false if could not find it */
				if( cdTemp(
					/* string for species & 21cm temperatures */
					chAssertLineLabel[i],
					/* will be mean temperatue */
					&relint,
					/* how to weight the average, must be "VOLUME" or "RADIUS" */
					chWeight ) )
				{
					fprintf( ioMONITOR, 
						" monitor error: lgCheckMonitors could not find an ion with label %s\n",
								chAssertLineLabel[i].c_str() );
					/* go to next line */
					lg1OK[i] = false;
					RelError[i] = 0;
					PredQuan[i] = 0;
					lgFound[i] = false;
					lgMonitorsOK = false;
					continue;
				}
			}
			/* this is the temperature */
			PredQuan[i] = relint;
			if (AssertError[i] > 0.)
				RelError[i] = get_error_ratio( PredQuan[i], AssertQuantity[i] );
			else
				RelError[i] = PredQuan[i]- AssertQuantity[i];
		}

		/* this would be grain potential in volt */
		else if( chAssertType[i] == "gp" )
		{
			/* the minus one is because the grain types are counted from one,
			 * but stuffed into the c array, that counts from zero */
			size_t nd = (size_t)wavelength[i]-1;
			if( nd >= gv.bin.size() )
			{
				fprintf( ioQQQ, "Illegal grain number found: %g\n" , wavelength[i] );
				fprintf( ioQQQ, "Use 1 for first grain that is turned on, " );
				fprintf( ioQQQ, "2 for second, etc....\n" );
				fprintf( ioQQQ, "Old style grain numbers are not valid anymore !!\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* get average grain potential in volt, always averaged over radius */
			PredQuan[i] = gv.bin[nd].avdpot/radius.depth_x_fillfac;
			/* actually absolute error, potential can be zero! */
			RelError[i] = AssertQuantity[i] -  PredQuan[i];
		}

		else
		{
			fprintf( ioMONITOR, 
					 " monitor error: lgCheckMonitors received an insane chAssertType=%s, impossible\n",
					 chAssertType[i].c_str() );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}

		if( chAssertLimit[i] == '=' )
		{
			/* predicted quantity should be within error of expected */
			if( fabs(RelError[i]) > fabs(AssertError[i]) )
			{
				lg1OK[i] = false;
				lgMonitorsOK = false;
			}
		}
		else if( chAssertLimit[i] == '<' )
		{
			/* expected is an upper limit, so PredQuan/AssertQuantity should
			 * be less than one, and so RelError =1-PredQuan[i]/AssertQuantity[i]
			 * should be >= 0 
			 * in case of iterations or zones, iter < iterations should
			 * trigger botched monitor when iter == iterations */
			if( RelError[i] <= 0. )
			{
				lg1OK[i] = false;
				lgMonitorsOK = false;
			}
		}
		else if( chAssertLimit[i] == '>' )
		{
			/* expected is a lower limit, so PredQuan/AssertQuantity should
			 * be greater than one, and so RelError should be negative */
			if( RelError[i] >= 0. )
			{
				lg1OK[i] = false;
				lgMonitorsOK = false;
			}
		}
	}

	/* only print summary if we are talking */
	if( called.lgTalk && nAsserts>0 )
	{ 
		time_t now;

		/* First disambiguate any line identifications */
		if( lgDisambiguate )
		{
			/* change significant figures of WL for this printout */
			long sigfigsav = LineSave.sig_figs;
			double relint1, relint2, absint1;

			LineSave.sig_figs = LineSave.sig_figs_max;

			fprintf( ioMONITOR, "=============Line Disambiguation============================================================\n" );
			fprintf( ioMONITOR, "                  Wavelengths                 ||                  Intensities               \n" );
			fprintf( ioMONITOR, "Label     line     match1   match2   match3   ||   asserted     match1     match2     match3\n" );

			for( i=0; i<nAsserts; ++i )
			{
				if( ipDisambiguate[i][1] > 0 )
				{
					fprintf( ioMONITOR , "%-*s ", NCHLAB-1, chAssertLineLabel[i].c_str() );
					prt_wl( ioMONITOR , wavelength[i] );
					fprintf( ioMONITOR , " " );
					prt_wl( ioMONITOR , LineSave.lines[ipDisambiguate[i][0]].wavelength() );
					fprintf( ioMONITOR , " " );
					prt_wl( ioMONITOR , LineSave.lines[ipDisambiguate[i][1]].wavelength() );
					fprintf( ioMONITOR , " " );
					if( ipDisambiguate[i][2] > 0 )
					{
						prt_wl( ioMONITOR , LineSave.lines[ipDisambiguate[i][2]].wavelength() );
						cdLine_ip( ipDisambiguate[i][2], &relint2, &absint1, iLineType[i] );
					}
					else
					{
						fprintf( ioMONITOR , "--------" );
						relint2 = 0.0;
					}
					fprintf( ioMONITOR , " ||" );

					cdLine_ip( ipDisambiguate[i][1], &relint1, &absint1, iLineType[i] );

					if( lgPrtSciNot )
					{
						fprintf( ioMONITOR , " %10.3e %10.3e %10.3e %10.3e\n",
							AssertQuantity[i],
							PredQuan[i] ,  
							relint1,
							relint2 );
					}
					else
					{
						fprintf( ioMONITOR , " %10.4f %10.4f %10.4f %10.4f\n",
							AssertQuantity[i],
							PredQuan[i] ,  
							relint1,
							relint2 );
					}
				}
			}
			fprintf( ioMONITOR, "\n" );

			/* revert to original significant figures */
			LineSave.sig_figs = sigfigsav;
		}

		/* write start of title and version number of code */
		fprintf( ioMONITOR, "=============Results of monitors: Cloudy %s ",
			 t_version::Inst().chVersion.c_str() );

		/* usually print date and time info - do not if "no times" command entered, 
		 * which set this flag false */
		if( prt.lgPrintTime ) 
		{
			/* now add date of this run */
			now = time(NULL);
			/* now print this time at the end of the string.  the system put cr at the end of the string */
			fprintf(ioMONITOR,"%s", ctime(&now) );
		}
		else 
		{
			fprintf(ioMONITOR,"\n" );
		}

		if( lgMonitorsOK )
		{
			fprintf( ioMONITOR, " No errors were found.  Summary follows.\n");
		}
		else
		{
			fprintf( ioMONITOR, " Errors were found.  Summary follows.\n");
		}

		fprintf( ioMONITOR, 
			"                   %-*s%*s    computed     asserted Rel Err Set err  type         \n",
				NCHLAB, "Label", LineSave.wl_length, "line" );
		/* now print a summary */
		for( i=0; i<nAsserts; ++i )
		{
			PrtOneMonitor( ioMONITOR, chAssertType[i], chAssertLineLabel[i], 
				wavelength[i], iLineType[i], PredQuan[i], chAssertLimit[i], AssertQuantity[i], 
				RelError[i], AssertError[i], lg1OK[i], lgQuantityLog[i], lgFound[i] );
		}
		fprintf( ioMONITOR , " \n");

		/* NB - in following, perl scripts detect these strings - be careful if they
		 * are changed - the scripts may no longer detect major errors */
		if( !lgMonitorsOK && lgBigBotch )
		{
			/* there were big botches */
			fprintf( ioMONITOR, " BIG BOTCHED MONITORS!!!   Big Botched Monitors!!! \n");
		}
		else if( !lgMonitorsOK )
		{
			fprintf( ioMONITOR, " BOTCHED MONITORS!!!   Botched Monitors!!! \n");
		}

		if( MonitorResults.nSumErrorCaseMonitor>0 )
		{
			fprintf(ioMONITOR,"\n The mean of the %li monitor Case A, B relative "
				"residuals is %.2f\n\n" , 
				MonitorResults.nSumErrorCaseMonitor,
				MonitorResults.SumErrorCaseMonitor /MonitorResults.nSumErrorCaseMonitor );
		}

		/* explain how we were compiled, but only if printing time */
		if( prt.lgPrintTime )
		{
			fprintf( ioQQQ, " %s\n\n", t_version::Inst().chInfo.c_str() );
		}
	}
	return lgMonitorsOK;
}

STATIC void prtLineType( FILE *ioMONITOR, const int iLineType )
{
	switch( iLineType )
	{
		case -1:
			/* The monitor is not an atomic transition. */
			fprintf( ioMONITOR, "           " );
			break;
		case 0:
			fprintf( ioMONITOR, "  intr     " );
			break;
		case 1:
			fprintf( ioMONITOR, "  emer     " );
			break;
		case 2:
			fprintf( ioMONITOR, "  intr cumu" );
			break;
		case 3:
			fprintf( ioMONITOR, "  emer cumu" );
			break;
		default:
			fprintf( ioMONITOR, "ERROR: Unrecognized line type: %d\n",
				iLineType );
			cdEXIT( EXIT_FAILURE );
			break;
	}
}

void PrtOneMonitor( FILE *ioMONITOR, const string& chAssertType, const string& chAssertLineLabel, 
		const realnum wavelength, const int iLineType, const double PredQuan, const char chAssertLimit, const double AssertQuantity, 
		const double RelError, const double AssertError, const bool lg1OK, const bool lgQuantityLog, const bool lgFound )
{
	DEBUG_ENTRY( "PrtOneMonitor()" );

	double prtPredQuan, prtAssertQuantity;
	// this is option to print log of quantity rather than linear.  default is
	// linear, and log only for a few such as ionization to see small numbers
	if( lgQuantityLog )
	{
		if( PredQuan > 0. )
			prtPredQuan = log10( MAX2( SMALLDOUBLE , PredQuan ) );
		else
			prtPredQuan = -37.;
		prtAssertQuantity = log10( MAX2( SMALLDOUBLE , AssertQuantity ) );
	}
	else
	{
		prtPredQuan = PredQuan;
		prtAssertQuantity = AssertQuantity;
	}

	/* start of line, flag problems */
	if( lg1OK )
	{
		// the ChkMonitor is a unique label so that we can grep on all output
		// and see what happened, without picking up input stream
		double relative = fabs( RelError / SDIV( fabs(AssertError)));

		if( relative < 0.25 || chAssertLimit != '=' )
		{
			fprintf( ioMONITOR, " ChkMonitor        ");
		}
		else if( relative < 0.50 )
		{
			fprintf( ioMONITOR, " ChkMonitor -      ");
		}
		else if( relative < 0.75 )
		{
			fprintf( ioMONITOR, " ChkMonitor --     ");
		}
		else if( relative < 0.90 )
		{
			fprintf( ioMONITOR, " ChkMonitor ---    ");
		}
		else  if( relative < 0.95 )
		{
			fprintf( ioMONITOR, " ChkMonitor ----   ");
		}
		else  if( relative < 0.98 )
		{
			fprintf( ioMONITOR, " ChkMonitor -----  ");
		}
		else 
		{
			fprintf( ioMONITOR, " ChkMonitor ------ ");
		}
	}
	else
	{
		fprintf( ioMONITOR, " ChkMonitor botch>>");
	}

	fprintf( ioMONITOR , "%-*s ", NCHLAB-1, chAssertLineLabel.c_str() );

	/* special formatting for the emission lines */
	if( chAssertType == "Ll"  || chAssertType == "Lr" || chAssertType == "Lb" )
	{
		prt_wl( ioMONITOR , wavelength );
	}
	else
	{
		fprintf( ioMONITOR , "%*i", LineSave.wl_length, (int)wavelength );
	}

	const char* format = " %10.4f %c %10.4f %7.3f %7.3f ";
	if( lgPrtSciNot )
		format = " %10.3e %c %10.3e %7.3f %7.3f ";

	fprintf( ioMONITOR, format,
		prtPredQuan,  
		chAssertLimit, 
		prtAssertQuantity, 
		RelError, 
		AssertError);

	prtLineType( ioMONITOR, iLineType );

	// if botched and the botch is > 3 sigma, say BIG BOTCH,
	// the lg1OK is needed since some tests (number of zones, etc)
	// are limits, not the quantity, and if limit is large the
	// miss will be big too

	/* >>chng 02 nov 27, added lgFound so don't get big botch when line simply missing */
	if( !lg1OK && (fabs(RelError) > 3.*AssertError) && lgFound )
	{
		/* >>chng 20 aug 03, discern between scope of big botch */
		if( fabs(RelError) > 9.*AssertError )
			fprintf( ioMONITOR , " <<BIG BIG (>9sig) BOTCH!!\n");
		else
			fprintf( ioMONITOR , " <<BIG (3sig) BOTCH!!\n");
		lgBigBotch = true;
	}
	else
	{
		fprintf( ioMONITOR , "\n");
	}

	return;
}
