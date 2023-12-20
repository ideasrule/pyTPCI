/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*cdDrive main routine to call cloudy under all circumstances) */
/*cdReasonGeo write why the model stopped and type of geometry on io file */
/*cdWarnings write all warnings entered into comment stack */
/*cdEms obtain the local emissivity for a line, for the last computed zone */
/*cdColm get the column density for a constituent  */
/*cdLine get the predicted line intensity, also index for line in stack */
/*cdLine_ip get the predicted line intensity, using index for line in stack */
/*cdCautions print out all cautions after calculation, on arbitrary io unit */
/*cdTemp_last routine to query results and return temperature of last zone */
/*cdDepth_depth get depth structure from previous iteration */
/*cdTimescales returns thermal, recombination, and H2 formation timescales */
/*cdSurprises print out all surprises on arbitrary unit number */
/*cdNotes print stack of notes about current calculation */
/*cdPressure_last routine to query results and return pressure of last zone */
/*cdTalk tells the code whether to print results or be silent */
/*cdOutput redirect output to arbitrary file */
/*cdInput redirect input from arbitrary file */
/*cdRead routine to read in command lines when cloudy used as subroutine */
/*cdErrors produce summary of all warnings, cautions, etc, on arbitrary io unit */
/*cdIonFrac get ionization fractions for a constituent */
/*cdTemp get mean electron temperature for any element */
/*cdCooling_last routine to query results and return cooling of last zone */
/*cdHeating_last routine to query results and return heating of last zone */
/*cdEDEN_last return electron density of last zone */
/*cdNoExec call this routine to tell code not to actually execute */
/*cdDate - puts date of code into string */
/*cdVersion produces string that gives version number of the code */
/*cdExecTime any routine can call this, find the time [s] since cdInit was called */
/*cdPrintCommands( FILE *) prints all input commands into file */
/*cdDrive main routine to call cloudy under all circumstances) */
/*cdNwcns get the number of cautions and warnings, to tell if calculation is ok */
/*debugLine provides a debugging hook into the main line array  */
/*cdEms_ip obtain the local emissivity for a line with known index */
/*cdnZone gets number of zones */
/*cdClosePunchFiles closes all the save files that have been used */
/*cdB21cm - returns B as measured by 21 cm */
/*cdPrtWL print line wavelengths in Angstroms in the standard format */

/* CHANGES: (M. Salz 21.05.2013)
 *  - define routine cdEDEN_depth( double )
 *    for obtaining electron density structure 
 *    after the simulation
 *  - cdDenPart_depth
 *    cdDenMass_depth
 *    cdWindVel_depth 
 *    cdCooling_depth
 *    cdHeating_depth 
 *    cdRadAcce_depth*/

#include "cddefines.h"
#include "cddrive.h"

#include "trace.h"
#include "conv.h"
#include "init.h"
#include "lines.h"
#include "pressure.h"
#include "dense.h"
#include "radius.h"
#include "struc.h"
#include "mole.h"
#include "elementnames.h"
#include "mean.h"
#include "phycon.h"
#include "called.h"
#include "parse.h"
#include "input.h"
#include "taulines.h"
#include "version.h"
#include "thermal.h"
#include "grid.h"
#include "timesc.h"
#include "cloudy.h"
#include "warnings.h"
#include "broke.h"
#include "iso.h"
#include "save.h"
#include "rfield.h"
#include "lines_service.h"
#include "service.h"
#include "parser.h"
#include "generic_state.h"
#include "prt.h"

/*************************************************************************
 *
 * cdDrive - main routine to call cloudy - returns 0 if all ok, 1 if problems
 *
 ************************************************************************/

int cdDrive()
{
	bool lgBAD = false;

	DEBUG_ENTRY( "cdDrive()" );
	/*********************************
	 * main routine to call cloudy   *
	 *********************************/

	/* this is set false when code loaded, set true when cdInit called,
	 * this is check that cdInit was called first */
	if( !lgcdInitCalled )
	{
		printf(" cdInit was not called first - this must be the first call.\n");
		cdEXIT(EXIT_FAILURE);
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			"cdDrive: lgOptimr=%1i lgVaryOn=%1i lgNoVary=%1i # of input lines=%li\n",
			 optimize.lgOptimr, optimize.lgVaryOn, optimize.lgNoVary, input.crd.size() );
	}

	if( prt.lgIncludeBlends )
	{
		// include standard blends -- this needs to be done here so that
		// it is also executed when calling cloudy as a subroutine...
		ParseInitFile("blends.ini");
	}

	/* should we call cloudy, or the optimization driver? 
	 * possible to have VARY on line without OPTIMIZE being set 
	 * lgNoVary set with "no optimize" command */
	if( optimize.lgOptimr && optimize.lgVaryOn && !optimize.lgNoVary && !grid.lgInsideGrid )
		optimize.lgVaryOn = true;
	else
		optimize.lgVaryOn = false;

	/* set up the frequency mesh, the SET CONTINUUM RESOLUTION factor was already parsed in cdRead() */
	/* the bool 'false' inidicates that the cell for the unit test is the highest entry in the regular mesh */
	rfield.InitMesh(false);
	rfield.nflux_with_check = rfield.ncells();
	/* there must be a cell above nflux for us to pass unity through the vol integrator */
	/*>>chng 04 oct 10, from nflux = nupper to nupper-1 since vectors must allocate to nupper, but
	 * will address [nflux] for unit continuum test */
	rfield.nflux = rfield.nflux_with_check-1;
	rfield.nPositive = rfield.nflux;

	/* one time initialization of core load - returns if already called 
	 * called here rather than in cdInit since at this point we know if
	 * single sim or grid */
	InitCoreload();

	if( optimize.lgVaryOn )
	{
		/* this branch called if optimizing or grid calculation */
		if( trace.lgTrace )
			fprintf( ioQQQ, "cdDrive: calling grid_do\n");
		/* option to drive optimizer set if OPTIMIZE was in input stream */
		grid_do();
	}
	else
	{
		if( trace.lgTrace )
			fprintf( ioQQQ, "cdDrive: calling cloudy\n");

		/* optimize did not occur, only compute one model, call cloudy */
		lgBAD = cloudy();
	}

	/* reset flag saying that cdInit has not been called */
	lgcdInitCalled = false;

	if( lgBAD )
	{
		if( trace.lgTrace )
			fprintf( ioQQQ, "cdDrive: returning failure during call. \n");
		/* something wrong, so return 1 */
		return 1;
	}
	else
	{
		/* everything is ok, return 0 */
		return 0;
	}
}


/*************************************************************************
 *
 * cdReasonGeo write why the model stopped and type of geometry on io file 
 *
 ************************************************************************/

/*cdReasonGeo write why the model stopped and type of geometry on io file */
void cdReasonGeo(FILE * ioOUT)
{
	DEBUG_ENTRY( "cdReasonGeo()" );

	/*this is the reason the calculation stopped*/
	for( size_t i=0; i < warnings.chRgcln.size(); ++i )
		fprintf( ioOUT, "%s\n", warnings.chRgcln[i].c_str() );
}


/*************************************************************************
 *
 * cdWarnings write all warnings entered into comment stack 
 *
 ************************************************************************/

/*cdWarnings write all warnings entered into comment stack */
void cdWarnings(FILE *ioOUT)
{
	DEBUG_ENTRY( "cdWarnings()" );

	for( size_t i=0; i < warnings.chWarnln.size(); i++ )
		/* these are all warnings that were entered in comment */
		fprintf( ioOUT, "%s\n", warnings.chWarnln[i].c_str() );
}


/*************************************************************************
 *
 * cdCautions print out all cautions after calculation, on arbitrary io unit 
 *
 ************************************************************************/

/*cdCautions print out all cautions after calculation, on arbitrary io unit */
void cdCautions(FILE * ioOUT)
{
	DEBUG_ENTRY( "cdCautions()" );

	for( size_t i=0; i < warnings.chCaunln.size(); i++ )
		fprintf( ioOUT, "%s\n", warnings.chCaunln[i].c_str() );
}

/*************************************************************************
 *
 * cdTimescales returns thermal, recombination, and H2 formation timescales 
 *
 ************************************************************************/

void cdTimescales(
	/* the thermal cooling timescale */
	double *TTherm , 
	/* the hydrogen recombination timescale */
	double *THRecom , 
	/* the H2 formation timescale */
	double *TH2 )
{

	DEBUG_ENTRY( "cdTimescales()" );

	/* these were all evaluated in AgeCheck, which was called by PrtComment */

	/* thermal or cooling timescale */
	*TTherm = timesc.time_therm_long;

	/* the hydrogen recombination timescale */
	*THRecom = timesc.time_Hrecom_long;

	/* longer of the the H2 formation and destruction timescales */
	*TH2 = MAX2( timesc.time_H2_Dest_longest , timesc.time_H2_Form_longest );
	return;
}


/*************************************************************************
 *
 * cdSurprises print out all surprises on arbitrary unit number 
 *
 ************************************************************************/

/*cdSurprises print out all surprises on arbitrary unit number */
void cdSurprises(FILE * ioOUT)
{
	DEBUG_ENTRY( "cdSurprises()" );

	for( size_t i=0; i < warnings.chBangln.size(); i++ )
		fprintf( ioOUT, "%s\n", warnings.chBangln[i].c_str() );

	if ( broke.lgPrintFixits )
	{
		fprintf (ioOUT, "\nActive fixits for this run\n");
		for (vector<string>::iterator it = FixitList::Inst().list.begin();
			  it != FixitList::Inst().list.end(); ++it)
		{
			fprintf(ioOUT,"%s\n",it->c_str());
		}
		fprintf (ioOUT, "\n");
	}
}


/*************************************************************************
 *
 * cdNotes print stack of notes about current calculation 
 *
 ************************************************************************/

/*cdNotes print stack of notes about current calculation */
void cdNotes(FILE * ioOUT)
{
	DEBUG_ENTRY( "cdNotes()" );

	for( size_t i=0; i < warnings.chNoteln.size(); i++ )
		fprintf( ioOUT, "%s\n", warnings.chNoteln[i].c_str() );
}

/*************************************************************************
 *
 * cdCooling_last routine to query results and return cooling of last zone 
 *
 ************************************************************************/

/*cdCooling_last routine to query results and return cooling of last zone */
double cdCooling_last() /* return cooling for last zone */
{
	return thermal.ctot;
}

/*************************************************************************
+ *
+ * cdCooling_depth routine to query results and return cooling structure 
+ *
+ ************************************************************************/

void cdCooling_depth( double Cool_struc[] )
{
  long int nz;
       
  DEBUG_ENTRY( "cdCooling_depth()" );
       
  for( nz = 0; nz<nzone; ++nz )
    {
      Cool_struc[nz] = struc.coolstr[nz];
    }
  return;
}


/*************************************************************************
 *
 * cdVersion - puts version number of code into string 
 * incoming string must have sufficient length and will become null
 * terminated string
 *
 ************************************************************************/

void cdVersion(char chString[])
{
	strcpy( chString , t_version::Inst().chVersion.c_str() );
	return;
}

/*************************************************************************
 *
 * cdHeating_depth routine to query results and return heating structure
 *
 ************************************************************************/

void cdHeating_depth( double Heat_struc[] )
{
       long int nz;

       DEBUG_ENTRY( "cdHeating_depth()" );

       for( nz = 0; nz<nzone; ++nz )
       {
               Heat_struc[nz] = struc.heatstr[nz];
       }
       return;
}

/*************************************************************************
 *
 * cdDenPart_depth get total particle density struc. from previous iteration 
 *
 ************************************************************************/
void cdDenPart_depth( double DenPart[] )
{
        long int nz;

        DEBUG_ENTRY( "cdDenPart_depth()" );

        for( nz = 0; nz<nzone; ++nz )
        {
                DenPart[nz] = struc.DenParticles[nz];
        }
        return;
}

/*************************************************************************
 *
 * cdDenMass_depth get total mass density struc. from previous iteration 
 *
 ************************************************************************/
void cdDenMass_depth( double DenMass[] )
{
        long int nz;

        DEBUG_ENTRY( "cdDenMass_depth()" );

        for( nz = 0; nz<nzone; ++nz )
        {
                DenMass[nz] = struc.DenMass[nz];
        }
        return;
}


/*************************************************************************
 *
 * cdWindVel_depth get total velocity struc. from previous iteration 
 *
 ************************************************************************/
void cdWindVel_depth( double WindVel[] )
{
        long int nz;

        DEBUG_ENTRY( "cdWindVel_depth()" );

        for( nz = 0; nz<nzone; ++nz )
        {
                WindVel[nz] = struc.windv[nz];
        }
        return;
}


/*************************************************************************
 *
 * cdRadAcce_depth get struc. of total radiative acceleration
 * from previous iteration (outward direction)
 *
 ************************************************************************/
void cdRadAcce_depth( double RadAccel[] )
{
        long int nz;

        DEBUG_ENTRY( "cdRadAcce_depth()" );

        for( nz = 0; nz<nzone; ++nz )
        {
                RadAccel[nz] = struc.AccelTotalOutward[nz];
        }
        return;
}

/*************************************************************************
 *
 * cdDate - puts date of code into string 
 * incoming string must have at least 8 char and will become null
 * terminated string
 *
 ************************************************************************/

/* cdDate - puts date of code into string  */
void cdDate(char chString[])
{
	strcpy( chString , t_version::Inst().chDate.c_str() );
	return;
}


/*************************************************************************
 *
 * cdHeating_last routine to query results and return heating of last zone
 *
 ************************************************************************/

/*cdHeating_last routine to query results and return heating of last zone */
double cdHeating_last() /* return heating for last zone */
{
	return thermal.htot;
}


/*************************************************************************
 *
 * cdEDEN_last return electron density of last zone
 *
 ************************************************************************/

/*cdEDEN_last return electron density of last zone */
double cdEDEN_last() /* return electron density for last zone */
{
	return dense.eden;
}

/*************************************************************************
 *
 * cdEDEN_depth get electron density structure from previous iteration 
 *
 ************************************************************************/
void cdEDEN_depth( double cdEDEN[] )
{
        long int nz;

        DEBUG_ENTRY( "cdEDEN_depth()" );

        for( nz = 0; nz<nzone; ++nz )
        {
                cdEDEN[nz] = struc.ednstr[nz];
        }
        return;
}

/*************************************************************************
 *
 * cdNoExec call this routine to tell code not to actually execute
 *
 ************************************************************************/

/*cdNoExec call this routine to tell code not to actually execute */
#include "noexec.h"

void cdNoExec()
{

	DEBUG_ENTRY( "cdNoExec()" );

	/* option to read in all input quantiles and NOT execute the actual model
	 * only check on input parameters - set by calling cdNoExec */
	noexec.lgNoExec = true;

	return;
}


/*************************************************************************
 *
 * cdSetExecTime routine to initialize variables keeping track of time at start of calculation
 *
 ************************************************************************/

/* set to false initially, then to true when cdSetExecTime() is called to
 * initialize the clock */
static bool lgCalled=false;

#if defined(_MSC_VER) 
/* _MSC_VER branches assume getrusage isn't implemented by MS */
struct timeval {
	long    tv_sec;         /* seconds */
	long    tv_usec;        /* microseconds */
};
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

/* will be used to save initial time */
static struct timeval before;

/* cdClock stores time since arbitrary datum in clock_dat           */
STATIC void cdClock(struct timeval *clock_dat)
{
	DEBUG_ENTRY( "cdClock()" );

/* >>chng 06 sep 2 rjrw: use long recurrence, fine grain UNIX clock *
 * -- maintain system dependences in a single routine               */
#if defined(_MSC_VER)
	clock_t clock_val;
	/* >>chng 05 dec 21, from above to below due to negative exec times */
	/*return (double)(clock() - before) / (double)CLOCKS_PER_SEC;*/
	clock_val = clock();
	clock_dat->tv_sec = clock_val/CLOCKS_PER_SEC;
	clock_dat->tv_usec = 1000000*(clock_val-(clock_dat->tv_sec*CLOCKS_PER_SEC))/CLOCKS_PER_SEC;
	/*>>chng 06 oct 05, this produced very large number, time typically 50% too long 
	clock_dat->tv_usec = 0;*/
#else
	struct rusage rusage;
	if(getrusage(RUSAGE_SELF,&rusage) != 0)
	{ 
		fprintf( ioQQQ, "DISASTER cdClock called getrusage with invalid arguments.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	clock_dat->tv_sec = rusage.ru_utime.tv_sec;
	clock_dat->tv_usec = rusage.ru_utime.tv_usec;
#endif

}
/* cdSetExecTime is called by cdInit when everything is initialized,
 * so that every time cdExecTime is called the elapsed time is returned */
void cdSetExecTime()
{
	cdClock(&before);
	lgCalled = true;
}
/* cdExecTime returns the elapsed time cpu time (sec) that has elapsed 
 * since cdInit called cdSetExecTime.*/
double cdExecTime()
{
	DEBUG_ENTRY( "cdExecTime()" );

	struct timeval clock_dat;
	/* check that we were properly initialized */
	if( lgCalled)
	{
		cdClock(&clock_dat);
		/*fprintf(ioQQQ,"\n DEBUG sec %.2e usec %.2e\n",
			(double)(clock_dat.tv_sec-before.tv_sec),
			1e-6*(double)(clock_dat.tv_usec-before.tv_usec));*/
		return (double)(clock_dat.tv_sec-before.tv_sec)+1e-6*(double)(clock_dat.tv_usec-before.tv_usec);
	}
	else
	{
		/* this is a big problem, we were asked for the elapsed time but
		 * the timer was not initialized by calling SetExecTime */
		fprintf( ioQQQ, "DISASTER cdExecTime was called before cdSetExecTime, impossible.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
}

// Maximum memory used, in kB
long cdMemory()
{
#if defined(_MSC_VER)
	return 0L;
#else
	struct rusage usage;
	if(getrusage(RUSAGE_SELF,&usage) != 0)
	{ 
		fprintf( ioQQQ, "DISASTER cdMemory called getrusage with invalid arguments.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return usage.ru_maxrss;
#endif
}

/*************************************************************************
 *
 * cdPrintCommands prints all input commands into file
 *
 ************************************************************************/

/* cdPrintCommands( FILE *)
 * prints all input commands into file */
void cdPrintCommands( FILE * ioOUT )
{
	fprintf( ioOUT, " Input commands follow:\n" );
	fprintf( ioOUT, "c ======================\n" );

	for( size_t i=0; i < input.crd.size(); i++ )
	{
		// exclude lines from init files
		if( input.crd[i]->InclLevel == 0 )
			fprintf( ioOUT, "%s\n", input.crd[i]->chCardSav.c_str() );
	}
	fprintf( ioOUT, "c ======================\n" );
}


/*************************************************************************
 *
 * cdEmis obtain the local emissivity for a line, for the last computed zone
 *
 ************************************************************************/

/** cdEms obtain the local emissivity for a line, for the last computed zone */


void cdEmis(
	const char *chLabel,
	realnum wavelength, 
	/* the vol emissivity of this line in last computed zone */
	double *emiss ,
	// intrinsic or emergent
	bool lgEmergent )
{
	DEBUG_ENTRY( "cdEmis()" );

	return cdEmis(LineID(chLabel, wavelength), emiss, lgEmergent);
}

void cdEmis(
	const LineID& line,
	double *emiss,
	bool lgEmergent )
{
	DEBUG_ENTRY( "cdEmis()" );

	long ipobs = LineSave.findline(line);
	if (ipobs >= 0)
		cdEmis(&LineSave.lines[ipobs],emiss,lgEmergent);
	else
		*emiss = 0.0;
	return;
}

void cdEmis_ip(
	/* index of the line in the stack */
	long int ipLine, 
	/* the vol emissivity of this line in last computed zone */
	double *emiss ,
	// intrinsic or emergent
	bool lgEmergent )
{
	DEBUG_ENTRY( "cdEmis_ip()" );

	/* returns the emissivity in a line - implements save lines structure 
	 * this version uses previously stored line index to save speed */
	ASSERT( ipLine >= 0 && ipLine < LineSave.nsum );
	cdEmis(&LineSave.lines[ipLine],emiss,lgEmergent);
}

/*************************************************************************
 *
 * cdColm get the column density for a constituent - 0 return if ok, 1 if problems 
 *
 ************************************************************************/

int cdColm(
	/* return value is zero if all ok, 1 if errors happened */
	/* 4-char + eol string that is first 
	 * 4 char of element name as spelled by Cloudy, upper or lower case */
	const char *chLabel,	

	/* integer stage of ionization, 1 for atom, 2 for A+, etc, 
	 * 0 is special flag for CO, H2, OH, or excited state */
	long int ion,

	/* the column density derived by the code [cm-2] */
	double *theocl )	
{
	long int nelem;
	/* use following to store local image of character strings */
	char chLABEL_CAPS[NCHLAB];

	DEBUG_ENTRY( "cdColm()" );
	*theocl = 0.;

	if( strlen(chLabel) > NCHLAB-1 )
	{
		fprintf( ioQQQ, " PROBLEM cdColm called with insane chLabel (between quotes) \"%s\", must be no more than %d characters long.\n",
					chLabel, NCHLAB-1 );
		return 1;
	}

	strcpy( chLABEL_CAPS, chLabel );

	/* convert element label to all caps */
	caps(chLABEL_CAPS);
	trimTrailingWhiteSpace( chLABEL_CAPS );

	genericState gs;

	/* zero ionization stage has special meaning.  The quantities recognized are
	 * the molecules, "H2  ", "OH  ", "CO  ", etc 
	 * "CII*" excited state C+ */
	if( ion < 0 )
	{
		fprintf( ioQQQ, " PROBLEM cdColm called with insane ion, =%li\n",
		  ion );
		return 1;
	}

	else if( ion == 0 )
	{
		/* this case molecular column density */
		/* want the *total* molecular hydrogen H2 column density */
		if( strcmp( chLABEL_CAPS , "H2" )==0 )
		{
			*theocl = findspecieslocal("H2")->column + findspecieslocal("H2*")->column;
		}
		/* H2g - ground H2 column density */
		else if( strcmp( chLABEL_CAPS , "H2G" )==0 )
		{
			gs.sp = findspecieslocal("H2");
			*theocl = column(gs);
		}
 		/* H2* - excited H2 - column density */
		else if( strcmp( chLABEL_CAPS , "H2*" )==0 )
		{
			gs.sp = findspecieslocal("H2*");
			*theocl = column(gs);
		}
		/* special option, "H2vJ" */
		else if( strncmp(chLABEL_CAPS , "H2" , 2 ) == 0 &&
					isdigit(chLABEL_CAPS[2]) && isdigit(chLABEL_CAPS[3]) )
		{
			long int iVib = chLABEL_CAPS[2] - '0';
			long int iRot = chLABEL_CAPS[3] - '0';
			if( iVib<0 || iRot < 0 )
			{
				fprintf( ioQQQ, " PROBLEM cdColm called with insane v,J for H2=\"%s\" caps=\"%s\"\n",
				  chLabel , chLABEL_CAPS );
				return 1;
			}
			*theocl = cdH2_colden( iVib , iRot );
		}

		/* ===========================================================*/
		/* end special case molecular column densities, start special cases
		 * excited state column densities */
		else if( strcmp( chLABEL_CAPS , "HE1*" )==0 )
		{
			if (dense.lgElmtOn[ipHELIUM])
			{
				gs.q = iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe2s3S];
			}
			*theocl = column(gs);
		}
		/* general case */
		else 
		{
			vector<genericState> v = matchGeneric(chLabel,false);
			*theocl = 0.;
			for (size_t i=0; i<v.size(); ++i)
			{
				*theocl += column(v[i]);
			}
		}
	}
	else
	{
		/* this case, ionization stage of some element */
		/* find which element this is */
		nelem = 0;
		while( nelem < LIMELM && 
			strncmp(chLABEL_CAPS,elementnames.chElementNameShort[nelem],4) != 0 )
		{
			++nelem;
		}

		/* this is true if we have one of the first 30 elements in the label,
		 * nelem is on C scale */
		if( nelem < LIMELM )
		{

			/* sanity check - does this ionization stage exist?
			 * max2 is to pick up H2 as H 3 */
			if( ion > MAX2(3,nelem + 2) )
			{
				fprintf( ioQQQ, 
				  " cdColm asked to return ionization stage %ld for element \"%s\" but this is not physical.\n", 
				  ion, chLabel );
				return 1;
			}

			/* the column density, ion is on physics scale, but means are on C scale */
			*theocl = mean.xIonMean[0][nelem][ion-1][0];
			/*>>chng 06 jan 23, div by factor of two
			 * special case of H2 when being tricked as H 3 - this stores 2H_2 so that
			 * the fraction of H in H0 and H+ is correct - need to remove this extra
			 * factor of two here */
			if( nelem==ipHYDROGEN && ion==3 )
				*theocl /= 2.;
		}
		else
		{
			fprintf( ioQQQ, 
			  " cdColm did not understand this combination of label \"%s\" and ion %4ld.\n", 
						chLabel, ion );
			return 1;
		}
	}
	return 0;
}


/*************************************************************************
 *
 *cdErrors produce summary of all warnings, cautions, etc, on arbitrary io unit 
 *
 ************************************************************************/

void cdErrors(FILE *ioOUT)
{
	long int nc, 
	  nn, 
	  npe, 
	  ns, 
	  nte, 
	  nw ,
	  nIone,
	  nEdene;

	DEBUG_ENTRY( "cdErrors()" );

	/* first check for number of warnings, cautions, etc */
	cdNwcns(&nw,&nc,&nn,&ns,&nte,&npe, &nIone, &nEdene );

	/* only say something is one of these problems is nonzero */
	if( nw || nc || nte || npe ||	nIone || nEdene )
	{
		/* say the title of the model */
		fprintf( ioOUT, "%75.75s\n", input.chTitle.c_str() );

		/* print warnings on the io unit */
		if( nw != 0 )
		{
			cdWarnings(ioOUT);
		}

		/* print cautions on the io unit */
		if( nc != 0 )
		{
			cdCautions(ioOUT);
		}

		if( nte != 0 )
		{
			fprintf( ioOUT , "Te failures=%4ld\n", nte );
		}

		if( npe != 0 )
		{
			fprintf( ioOUT , "Pressure failures=%4ld\n", npe );
		}

		if( nIone != 0 )
		{
			fprintf( ioOUT , "Ionization failures=%4ld\n", nte );
		}

		if( nEdene != 0 )
		{
			fprintf( ioOUT , "Electron density failures=%4ld\n", npe );
		}
	}
	return;
}

/*************************************************************************
 *
 *cdDepth_depth get depth structure from previous iteration 
 *
 ************************************************************************/
void cdDepth_depth( double cdDepth[] )
{
	long int nz;

	DEBUG_ENTRY( "cdDepth_depth()" );

	for( nz = 0; nz<nzone; ++nz )
	{
		cdDepth[nz] = struc.depth[nz];
	}
	return;
}

/*************************************************************************
 *
 *cdPressure_depth routine to query results and return pressure of last iteration 
 *
 ************************************************************************/

/*
 * cdPressure_depth
 * This returns the pressure and its constituents for the last iteration. 
 * space was allocated in the calling routine for the vectors - 
 * before calling this, cdnZone should have been called to get the number of
 * zones, then space allocated for the arrays */
void cdPressure_depth(
	/* total pressure, all forms*/
	double TotalPressure[],			
	/* gas pressure */
	double GasPressure[],				
	/* line radiation pressure */
	double RadiationPressure[])
{
	long int nz;

	DEBUG_ENTRY( "cdPressure_depth()" );

	for( nz = 0; nz<nzone; ++nz )
	{
		TotalPressure[nz] = struc.pressure[nz];
		GasPressure[nz] = struc.GasPressure[nz];
		RadiationPressure[nz] = struc.pres_radiation_lines_curr[nz];
	}
	return;
}

/*************************************************************************
 *
 *cdPressure_last routine to query results and return pressure of last zone 
 *
 ************************************************************************/

void cdPressure_last(
		double *PresTotal,  /* total pressure, all forms, for the last computed zone*/
		double *PresGas,    /* gas pressure */
		double *PresRad)    /* line radiation pressure */
{

	DEBUG_ENTRY( "cdPressure_last()" );

	*PresGas = pressure.PresGasCurr;
	*PresRad = pressure.pres_radiation_lines_curr;
	*PresTotal = pressure.PresTotlCurr;
	return;
}

/*************************************************************************
 *
 *cdnZone gets number of zones
 *
 ************************************************************************/

/* returns number of zones */
long int cdnZone() 
{
	return nzone;
}

/*************************************************************************
 *
 *cdTemp_last routine to query results and return temperature of last zone 
 *
 ************************************************************************/


double cdTemp_last()
{
	return phycon.te;
}

/*************************************************************************
 *
 *cdIonFrac get ionization fractions for a constituent
 *
 ************************************************************************/

STATIC inline bool set_weighting( const char *chCaller, const char *chWeight, int &dim )
{
	DEBUG_ENTRY( "set_weighting()" );

	string chCARD = chWeight;
	caps(chCARD);

	if( chCARD == "RADIUS" )
	{
		dim = 0;
		return true;
	}
	else if( chCARD == "AREA" )
	{
		dim = 1;
		return true;
	}
	else if( chCARD == "VOLUME" )
	{
		dim = 2;
		return true;
	}
	else
	{
		fprintf( ioQQQ,
			" %s: chWeight=%6.6s makes no sense to me, valid options are RADIUS, AREA, and VOLUME\n",
		  	chCaller, chWeight );
		return false;
	}
}

int cdIonFrac(
	/* four char string, null terminated, giving the element name */
	const char *chLabel, 
	/* IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number,
	 * 0 says special case */
	long int IonStage, 
	/* will be fractional ionization */
	double *fracin, 
	/* how to weight the average, must be "VOLUME", "AREA", or "RADIUS" */
	const char *chWeight ,
	/* if true then weighting also has electron density, if false then only volume or radius */
	bool lgDensity ) 
	/* return value is 0 if element was found, non-zero if failed */
{
	long int ip, 
		ion, /* used as index within aaa vector*/
		nelem;
	realnum aaa[LIMELM + 1];
	/* use following to store local image of character strings */
	string chCARD;

	DEBUG_ENTRY( "cdIonFrac()" );


	int dim = 0;
	if( ! set_weighting( __func__, chWeight, dim ) )
	{
		*fracin = 0.;
		return 1;
	}

	chCARD = chLabel;
	caps(chCARD);
	trimTrailingWhiteSpace(chCARD);
	trimTrailingWhiteSpace(chCARD);

	if( IonStage==0 )
	{
		/* special case */
		if( chCARD == "H2" )
		{
			/* this will be request for H2, treated as third stage of hydrogen */
			nelem = 0;
			IonStage = 3;
		}
		else
		{
			fprintf( ioQQQ, " cdIonFrac: ion stage of zero and element %s makes no sense to me\n", 
				 chCARD.c_str() );
			*fracin = 0.;
			return 1;
		}
	}

	else 
	{
		/* find which element this is */
		nelem = 0;
		while( nelem < LIMELM && 
		       chCARD != elementnames.chElementNameShort[nelem] )
		{
			++nelem;
		}
	}

	/* if element not recognized and above loop falls through, nelem is LIMELM */
	if( nelem >= LIMELM )
	{
		fprintf( ioQQQ, " cdIonFrac called with unknown element chLabel, =%4.4s\n", 
		  chLabel );
		return 1;
	}

	/* sanity check - does this ionization stage exist? 
	 * IonStage is on spectroscopic scale and nelem is on C scale */

	/* ion will be used as pointer within the aaa array that contains average values,
	 * convert to C scale */
	ion = IonStage - 1;

	if( (ion > nelem+1 || ion < 0 ) && !(nelem==ipHYDROGEN&&ion==2))
	{
		fprintf( ioQQQ, " cdIonFrac asked to return ionization stage %ld for element %4.4s but this is not physical.\n", 
		  IonStage, chLabel );
		*fracin = -1.;
		return 1;
	}

	/* get average, aaa is filled in from 0 */
	/* aaa is dim'd LIMELM+1 so largest argument is LIMELM
	 * 'i' means ionization, not temperature */
	/* last argument says whether to include electron density */
	/* MeanIon uses c scale for nelem */
	mean.MeanIon('i',nelem,dim,&ip,aaa,lgDensity);
	*fracin = exp10((double)aaa[ion]);

	/* we succeeded - say so */
	return 0;
}

/*************************************************************************
 *
 * debugLine provides a debugging hook into the main line array 
 *
 ************************************************************************/

 /* debugLine provides a debugging hook into the main line array 
  * loops over whole array and finds every line that matches length,
  * the wavelength, the argument to the function
  * put breakpoint inside if test 
  * return value is number of matches, also prints all matches*/
long debugLine( realnum wavelength )
{
	long j, kount;
	realnum errorwave;

	kount = 0;

	/* get the error associated with specified significant figures */
	errorwave = WavlenErrorGet( wavelength, LineSave.sig_figs );

	for( j=0; j < LineSave.nsum; j++ )
	{
		/* check wavelength and chLabel for a match */
		/* if( fabs(LineSave.lines[j].wavelength- wavelength)/MAX2(DELTA,wavelength) < errorwave ) */
		if( fabs(LineSave.lines[j].wavelength()-wavelength) < errorwave )
		{
			LineSave.lines[j].prt(stdout);
			printf("\n");
			++kount;
		}
	}
	printf(" hits = %li\n", kount );
	return kount;
}

/*************************************************************************
 *
 *cdLine get the predicted line intensity, also index for line in stack 
 *
 ************************************************************************/

// returns array index for line in array stack if we found the line,
// return negative of total number of lines as debugging aid if line not found
long int cdLine(
	const string& chLabel, 
	/* wavelength of line in angstroms, not format printed by code */
	realnum wavelength, 
	/* linear intensity relative to normalization line*/
	double *relint, 
	/* log of luminosity or intensity of line */
	double *absint ,
	// 0 is intrinsic,
	// 1 emergent
	// 2 is intrinsic cumulative,
	// 3 emergent cumulative
	int LineType )
{
	DEBUG_ENTRY( "cdLine()" );

	return cdLine(LineID(chLabel, wavelength), relint, absint, LineType);
}

long int cdLine(
	const LineID& line,
	/* linear intensity relative to normalization line*/
	double *relint, 
	/* log of luminosity or intensity of line */
	double *absint ,
	// 0 is intrinsic,
	// 1 emergent
	// 2 is intrinsic cumulative,
	// 3 emergent cumulative
	int LineType )
{
	DEBUG_ENTRY( "cdLine()" );

	*absint = 0.;
	*relint = 0.;

	long ipobs = LineSave.findline(line);
	if (ipobs >= 0)
		cdLine_ip(ipobs,relint,absint,LineType);
	//	printf("%s\t %g\t ip = %ld\t absint = %g\n", chLabel, wavelength, ipobs, *absint);

	return ipobs;
}

void cdLine_ip(long int ipLine, 
	  /* linear intensity relative to normalization line*/
	  double *relint, 
	  /* log of luminosity or intensity of line */
	  double *absint ,
	  // 0 is intrinsic,
	  // 1 emergent
	  // 2 is intrinsic cumulative,
	  // 3 emergent cumulative
	  int LineType )
{
	DEBUG_ENTRY( "cdLine_ip()" );

	*relint = 0.;
	*absint = 0.;

	if( LineType<0 || LineType>3 )
	{
		fprintf( ioQQQ, " PROBLEM cdLine_ip called with insane nLineType - it must be between 0 and 3.\n");
		return;
	}

	/* this is zero when cdLine called with cdNoExec called too */
	if( LineSave.nsum == 0 )
	{
		return;
	}
	ASSERT( LineSave.ipNormWavL >= 0 );
	ASSERT( LineSave.nsum > 0 );

	/* does the normalization line have a positive intensity */
	if( LineSave.lines[LineSave.ipNormWavL].SumLine(LineType) > 0. )
		*relint = LineSave.lines[ipLine].SumLine(LineType)/
			LineSave.lines[LineSave.ipNormWavL].SumLine(LineType)*
			LineSave.ScaleNormLine;

	/* return log of current line intensity if it is positive */
	if( LineSave.lines[ipLine].SumLine(LineType) > 0. )
		*absint = LineSave.lines[ipLine].SumLine(LineType) *
			radius.Conv2PrtInten;

	return;
}

/*************************************************************************
 *
 *cdNwcns get the number of cautions and warnings, to tell if calculation is ok
 *
 ************************************************************************/

void cdNwcns(
	/* the number of warnings, cautions, notes, and surprises */
	long int *NumberWarnings, 
	long int *NumberCautions, 
	long int *NumberNotes, 
	long int *NumberSurprises, 
	/* the number of temperature convergence failures */
	long int *NumberTempFailures, 
	/* the number of pressure convergence failures */
	long int *NumberPresFailures,
	/* the number of ionization convergence failures */
	long int *NumberIonFailures, 
	/* the number of electron density convergence failures */
	long int *NumberNeFailures )
{

	DEBUG_ENTRY( "cdNwcns()" );

	/* this sub is called after comment, to find the number of various comments */
	*NumberWarnings = warnings.chWarnln.size();
	*NumberCautions = warnings.chCaunln.size();
	*NumberNotes = warnings.chNoteln.size();
	*NumberSurprises = warnings.chBangln.size();

	/* these are counters that were incremented during convergence failures */
	*NumberTempFailures = conv.nTeFail;
	*NumberPresFailures = conv.nPreFail;
	*NumberIonFailures = conv.nIonFail;
	*NumberNeFailures = conv.nNeFail;
	return;
}

void cdOutput( const string& filename, const char *mode )
{
	DEBUG_ENTRY( "cdOutput()" );

	if( ioQQQ != stdout && ioQQQ != NULL )
		fclose(ioQQQ);
	FILE *fp = stdout;
	if( !filename.empty() )
		fp = open_data( filename, mode );
	save.chOutputFile = filename;
	ioQQQ = fp;
}

void cdOutput( const string& filename, FILE* fp )
{
	DEBUG_ENTRY( "cdOutput()" );

	if( ioQQQ != stdout && ioQQQ != NULL )
		fclose(ioQQQ);
	save.chOutputFile = filename;
	ioQQQ = ( fp != NULL ) ? fp : stdout;
}

void cdInput( const string& filename, const char *mode )
{
	DEBUG_ENTRY( "cdInput()" );

	if( ioStdin != stdin && ioStdin != NULL )
		fclose(ioStdin);
	FILE *fp = stdin;
	if( !filename.empty() )
	{
		fp = open_data( filename, mode, AS_LOCAL_ONLY_TRY );
		if( fp == NULL )
		{
			fprintf( ioQQQ, " input file \"%s\" not found\n", filename.c_str() );
			cdEXIT(ES_FAILURE);
		}
	}
	ioStdin = fp;
}

/*************************************************************************
 *
 * cdTalk tells the code whether to print results or be silent
 *
 ************************************************************************/

void cdTalk(bool lgTOn)
{

	DEBUG_ENTRY( "cdTalk()" );

	/* MPI talking rules must be obeyed, otherwise loss of output may result */
	called.lgTalk = lgTOn && cpu.i().lgMPI_talk();
	/* has talk been forced off? */
	called.lgTalkForcedOff = !lgTOn;
	return;
}

/*cdB21cm - returns B as measured by 21 cm */
double cdB21cm()
{
	double ret;

	DEBUG_ENTRY( "cdB21cm()" );

	// return average over radius
	if( mean.TempB_HarMean[0][1] > SMALLFLOAT )
	{
		ret = mean.TempB_HarMean[0][0]/mean.TempB_HarMean[0][1];
	}
	else
	{
		ret = 0.;
	}
	return ret;
}

/*************************************************************************
 *
 * cdTemp get mean electron temperature for any element 
 *
 ************************************************************************/
/* This routine finds the mean electron temperature for any ionization stage 
 * It returns 0 if it could find the species, 1 if it could not find the species.
 * The first argument is a null terminated 4 char string that gives the element
 * name as spelled by Cloudy.  
 * The second argument is ion stage, 1 for atom, 2 for first ion, etc
 * This third argument will be returned as result,
 * Last parameter is either "RADIUS", "AREA", or "VOLUME" to give weighting 
 *
 * if the ion stage is zero then the element label will have a special meaning.
 * The string "21CM" will return the 21 cm temperature.
 * The string "H2  " will return the temperature weighted wrt the H2 density
 * There are several other options as listed below */
/** \todo	2	this should have a last argument like cdIonFrac for whether or not weighting
 * is wrt electron density */

/* return value is o if things are ok and element was found, 
 * non-zero if element not found or there are problems */
int cdTemp(
	/* four char string, null terminated, giving the element name */
	const char *chLabel, 
	/* IonStage is ionization stage, 1 for atom, up to Z+1 where Z is atomic number,
	 * 0 means that chLabel is a special case */
	long int IonStage, 
	/* will be temperature */
	double *TeMean, 
	/* how to weight the average, must be "RADIUS", "AREA, or "VOLUME" */
	const char *chWeight ) 
{
	long int ip, 
		ion, /* used as pointer within aaa vector*/
		nelem;
	realnum aaa[LIMELM + 1];
	/* use following to store local image of character strings */
	string chELEM;

	DEBUG_ENTRY( "cdTemp()" );

	/* make sure strings are all caps */
	chELEM = chLabel;
	caps(chELEM);
	trimTrailingWhiteSpace(chELEM);

	int dim = 0;
	if( ! set_weighting( __func__, chWeight, dim ) )
	{
		*TeMean = 0.;
		return 1;
	}

	if( IonStage == 0 )
	{
		/* return atomic hydrogen weighted harmonic mean of gas kinetic temperature */
		if( chELEM == "21CM" )
		{
			if( mean.TempHarMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempHarMean[dim][0]/mean.TempHarMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* return atomic hydrogen weighted harmonic mean 21 cm spin temperature,
		 * using actual level populations with 1s of H0 */
		else if( chELEM == "SPIN" )
		{
			if( mean.TempH_21cmSpinMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempH_21cmSpinMean[dim][0] / mean.TempH_21cmSpinMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* return temperature deduced from ratio of 21 cm and Lya optical depths */
		else if( chELEM == "OPTI" )
		{
			if( HFLines[0].Emis().TauCon() > SMALLFLOAT )
				*TeMean = 3.84e-7 * iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().TauCon() /
					HFLines[0].Emis().TauCon();
			else
				*TeMean = 0.;
		}
		/* mean temp weighted by H2 density */
		else if( chELEM == "H2" )
		{
			if( mean.TempIonMean[dim][ipHYDROGEN][2][1] > SMALLFLOAT )
				*TeMean = mean.TempIonMean[dim][ipHYDROGEN][2][0] /
					mean.TempIonMean[dim][ipHYDROGEN][2][1];
			else
				*TeMean = 0.;
		}
		/* this is temperature weighted by eden */
		else if( chELEM == "TENE" )
		{
			if( mean.TempEdenMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempEdenMean[dim][0]/mean.TempEdenMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* four spaces mean to return simple mean of Te */
		else if( chELEM == "" )
		{
			if( mean.TempMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempMean[dim][0]/mean.TempMean[dim][1];
			else
				*TeMean = 0.;
		}
		else
		{
			fprintf( ioQQQ, " cdTemp called with ion=0 and unknown quantity, =%4.4s\n", 
			  chLabel );
			return 1;
		}

		/* say things are ok */
		return 0;
	}

	/* find which element this is */
	nelem = 0;
	while( nelem < LIMELM && 
	       chELEM != elementnames.chElementNameShort[nelem] )
	{
		++nelem;
	}

	/* if element not recognized and above loop falls through, nelem is LIMELM */
	if( nelem >= LIMELM )
	{
		fprintf( ioQQQ, " cdTemp called with unknown element chLabel, =%4.4s\n", 
		  chLabel );
		return 1;
	}

	/* sanity check - does this ionization stage exist? 
	 * IonStage on spectroscopic scale, nelem on c */

	/* ion will be used as pointer within the aaa array that contains average values,
	 * done this way to prevent lint from false problem in access to aaa array */
	ion = IonStage - 1;

	if( ion > nelem+1 || ion < 0 )
	{
		fprintf( ioQQQ, " cdTemp asked to return ionization stage %ld for element %4.4s but this is not physical.\n", 
		  IonStage, chLabel );
		return 1;
	}

	/* get average, aaa is filled in from 0 */
	/* aaa is dim'd LIMELM+1 so largest arg is LIMELM */
	/* MeanIon uses C scale for nelem */
	mean.MeanIon('t', nelem,dim,&ip,aaa,false);
	*TeMean = exp10((double)aaa[ion]);
	return 0;
}

int cdTemp(
	/* string giving the species, or a flag */
	const string &chLabel, 
	/* will be temperature */
	double *TeMean, 
	/* how to weight the average, must be "RADIUS", "AREA, or "VOLUME" */
	const char *chWeight ) 
{
	DEBUG_ENTRY( "cdTemp()" );

	/* use following to store local image of character strings */
	string chELEM = chLabel;
	trimTrailingWhiteSpace(chELEM);

	int dim = 0;
	if( ! set_weighting( __func__, chWeight, dim ) )
	{
		*TeMean = 0.;
		return 1;
	}

	if( ! isSpecies( chELEM ) )
	{
		/* make sure strings are all caps */
		caps(chELEM);

		/* return atomic hydrogen weighted harmonic mean of gas kinetic temperature */
		if( chELEM == "21CM" )
		{
			if( mean.TempHarMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempHarMean[dim][0]/mean.TempHarMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* return atomic hydrogen weighted harmonic mean 21 cm spin temperature,
		 * using actual level populations with 1s of H0 */
		else if( chELEM == "SPIN" )
		{
			if( mean.TempH_21cmSpinMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempH_21cmSpinMean[dim][0] / mean.TempH_21cmSpinMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* return temperature deduced from ratio of 21 cm and Lya optical depths */
		else if( chELEM == "OPTI" )
		{
			if( HFLines[0].Emis().TauCon() > SMALLFLOAT )
				*TeMean = 3.84e-7 * iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().TauCon() /
					HFLines[0].Emis().TauCon();
			else
				*TeMean = 0.;
		}
		/* this is temperature weighted by eden */
		else if( chELEM == "TENE" )
		{
			if( mean.TempEdenMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempEdenMean[dim][0]/mean.TempEdenMean[dim][1];
			else
				*TeMean = 0.;
		}
		/* an empty string means to return simple mean of Te */
		else if( chELEM == "" )
		{
			if( mean.TempMean[dim][1] > SMALLFLOAT )
				*TeMean = mean.TempMean[dim][0]/mean.TempMean[dim][1];
			else
				*TeMean = 0.;
		}
		else
		{
			fprintf( ioQQQ, " W- cdTemp called with uknown species '%s'\n", chLabel.c_str() );
			return 1;
		}

		/* say things are ok */
		return 0;
	}

	if( isMolecule( chLabel ) )
	{
		if( mean.MeanMoleculeTemp( chLabel, dim, *TeMean ) )
		{
			return 1;
		}
		return 0;
	}

	long spCharge;

	if( ! parse_chemical( chLabel, chELEM, spCharge ) )
	{
		fprintf( ioQQQ, "cdTemp could not parse species '%s'\n", chLabel.c_str() );
		return 1;
	}

	int nelem = elem_symbol_to_index( chELEM );
	if( nelem < 0 )
	{
		fprintf( ioQQQ, " cdTemp found invalid element (index: %d) in species '%s'\n",
				nelem, chLabel.c_str() );
		return 1;
	}

	/* spCharge will be used as pointer within the aaa array that contains average values,
	 * done this way to prevent lint from false problem in access to aaa array */
	if( ! isAtomicIonValid( nelem, spCharge ) )
	{
		fprintf( ioQQQ, " cdTemp asked to return charge %ld for element '%s' but this is not physical.\n",
				spCharge, chELEM.c_str() );
		return 1;
	}

	/* get average, aaa is filled in from 0 */
	/* aaa is dim'd LIMELM+1 so largest arg is LIMELM */
	/* MeanIon uses C scale for nelem */
	long int ip;
	realnum aaa[LIMELM + 1];

	mean.MeanIon('t', nelem,dim,&ip,aaa,false);
	*TeMean = exp10( double( aaa[ spCharge ] ) );
	return 0;
}

void test_cdTemp_molecules()
{
	DEBUG_ENTRY( "test_cdTemp_molecules()" );

	vector<string> allMolecules;
	getMolecules( allMolecules );

	fprintf( ioQQQ,
		"cdTemp() molecules\n"
		"Molecule\t\tRADIUS\tVOLUME\n"
		"--------\t\t------\t------\n" );

	for( vector<string>::iterator it = allMolecules.begin();
		it != allMolecules.end(); it++ )
	{
		const string this_str = *it;
		double TeMean;
		cdTemp( this_str, &TeMean, "RADIUS" );
		fprintf( ioQQQ, "%10.10s\t%5.2f\t", it->c_str(), TeMean );
		cdTemp( this_str, &TeMean, "VOLUME" );
		fprintf( ioQQQ, "%5.2f\n", TeMean );
	}

	return;
}



/*************************************************************************
 *
 * cdRead routine to read in command lines 
 * called by user when cloudy used as subroutine 
 * called by cdMain and ParseInit when used as standalone program
 *
 ************************************************************************/

/* returns the number of lines available in command stack 
 * this is limit to how many more commands can be read */
int cdRead( const string& chInputLine ) /* the string containing the command */
{
	const bool KEEP_VISIBLE = false;

	DEBUG_ENTRY( "cdRead()" );

	/* this is set false when code loaded, set true when cdInit called,
	 * this is check that cdInit was called first */
	if( !lgcdInitCalled )
	{
		printf(" cdInit was not called first - this must be the first call.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* totally ignore hidden comment lines but want to include visible
	 * comments since they are copied to output */
	if( lgIsCommentSeq( chInputLine, 0, KEEP_VISIBLE ) )
	{
		/* return value is not used, but is kept for backward compatibility
		 * with programs that call cloudy as a subroutine */
		return INT_MAX;
	}

	/* use local version to strip eol characters and hidden comments */
	string chLocal = chInputLine;

	/* erase EOL character and trailing whitespace */
	size_t pp;
	if( (pp = chLocal.find_first_of("\n\r")) != string::npos )
		chLocal.erase(pp);
	trimTrailingWhiteSpace(chLocal);

	/* now kill hidden comments at the end of a command
	 * this stops info user wants ignored from entering after here
	 * also remove possible underscores and square brackets */
	StripComment( chLocal, KEEP_VISIBLE );

	/* save string in master array for later use in output */
	auto ci = new CardInfo;
	ci->chCardSav = chLocal;
	ci->InclLevel = input.curInclLevel;
	input.crd.push_back( ci );

	Parser p;
	p.setline(chLocal.c_str());

	/* check whether this is a trace command, turn on printout if so */
	if( p.hasCommand("TRACE") )
		trace.lgTrace = true;

	/* print input lines if trace specified, need to do this before INIT command is processed */
	if( trace.lgTrace )
		fprintf( ioQQQ,"cdRead [%d] =%s=\n",ci->InclLevel,ci->chCardSav.c_str() );

	if( p.hasCommand("NO") && p.nMatch("VARY") )
		optimize.lgNoVary = true;

	/* prt.lgPrintTime can be used by check_data() before regular parsing starts */
	if( p.hasCommand("NO") && p.nMatch("TIME") )
		prt.lgPrintTime = false;

	if( p.hasCommand("GRID") )
	{
		optimize.lgOptimr = true;
		grid.lgGrid = true;
		++grid.nGridCommands;
	}

	/* NO BUFFERING turn off buffered io for standard output, 
	 * used to get complete output when code crashes */
	if( p.hasCommand("NO") && p.nMatch("BUFF") )
	{
		/* if output has already been redirected (e.g. using cdOutput()) then
		 * ignore this command, a warning will be printed in ParseDont() */
		/* stdout may be a preprocessor macro, so lets be really careful here */
		FILE *test = stdout;
		if( ioQQQ == test )
		{
			fprintf( ioQQQ, " cdRead found NO BUFFERING command, redirecting output to stderr now.\n" );
			/* make sure output is not lost */
			fflush( ioQQQ );
			/* stderr is always unbuffered */
			/* don't use setvbuf() here since we cannot guarantee
			 * that no operations have been performed on stdout */
			ioQQQ = stderr;
			/* will be used to generate comment at end */
			input.lgSetNoBuffering = true;
		}
		else if( !save.chOutputFile.empty() )
		{
			fprintf( ioQQQ, " cdRead found NO BUFFERING command, reopening file %s now.\n",
				 save.chOutputFile.c_str() );
			fclose( ioQQQ );
			// using AS_SILENT_TRY assures that open_data will not write to ioQQQ
			ioQQQ = open_data( save.chOutputFile, "a", AS_SILENT_TRY );
			if( ioQQQ == NULL )
			{
				// ioQQQ is no longer valid, so switch to stderr and abort...
				ioQQQ = stderr;
				fprintf( ioQQQ, " cdRead failed to reopen %s, aborting!\n",
					 save.chOutputFile.c_str() );
				cdEXIT(EXIT_FAILURE);
			}
			if( setvbuf( ioQQQ, NULL, _IONBF, 0 ) != 0 )
				fprintf( ioQQQ, " PROBLEM -- cdRead failed to set NO BUFFERING mode.\n" );
			else
				input.lgSetNoBuffering = true;
		}
	}

	/* option to not include blends */
	if( p.hasCommand("NO") && p.nMatch("BLEN") )
		prt.lgIncludeBlends = false;

	/* optimize command read in */
	if( p.hasCommand("OPTIMIZE") )
		optimize.lgOptimr = true;

	if( p.hasCommand("SET") && p.nMatch("CONT") && p.nMatch("RESO") )
	{
		double factor = p.getNumberDefault("",1.);
		/* negative numbers are always logs */
		if( factor <= 0. )
			factor = exp10(factor);
		rfield.setResolutionScaleFactor(factor);
	}

	/* if the command is an init command, process it immediately
	 * so that the lines in the init file can be inserted in situ */
	if( p.hasCommand("INIT") )
	{
		ParseInit(p);
		input.lgInitPresent = true;
	}

	/* here we check for specific keywords, avoiding the TITLE command since
	 * that contains free text, which could lead to spurious matches... */
	if( ! p.hasCommand("TITLE") )
	{
		/* now check whether VARY is specified */
		if( p.nMatch("VARY") )
			/* a optimize flag was on the line image */
			optimize.lgVaryOn = true;
	}

	return INT_MAX;
}

/* wrapper to close all save files */
void cdClosePunchFiles()
{

	DEBUG_ENTRY( "cdClosePunchFiles()" );

	CloseSaveFiles( true );
	return;
}
