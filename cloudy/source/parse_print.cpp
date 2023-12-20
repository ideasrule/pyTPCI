/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParsePrint parse the print command  */
/*prt_constants print physical constants */
/* PrtMacros - print macros in cddefines.h and their status -  *print macros* */
#include "cddefines.h"
#include "broke.h"
#include "rfield.h"
#include "iso.h"
#include "iterations.h"
#include "called.h"
#include "elementnames.h"
#include "ionbal.h"
#include "prt.h"
#include "lines.h"
#include "h2.h"
#include "parser.h"
#include "version.h"
#include "atmdat.h"
#include "input.h"

/*prt_constants print physical constants */
STATIC void prt_constants(void);


// PrtMacros - print macros in cddefines.h and their status -  *print macros*
STATIC void PrtMacros( void )
{

	DEBUG_ENTRY( "PrtMacros()" );

	fprintf( ioQQQ," PrtMacros:\n FLT_IS_DBL\t");
#	ifdef FLT_IS_DBL
	fprintf( ioQQQ,"defined.\n");
#	endif

	fprintf( ioQQQ , "\n SYS_CONFIG\t");
#	ifdef SYS_CONFIG
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n MPI_GRID_RUN\t");
#	ifdef MPI_GRID_RUN
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n USE_GPROF\t");
#	ifdef USE_GPROF
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n _MSC_VER\t");
#	ifdef _MSC_VER
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n __INTEL_COMPILER\t");
#	ifdef __INTEL_COMPILER
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POWI\t");
#	ifdef HAVE_POWI
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_DOUBLE_INT\t");
#	ifdef HAVE_POW_DOUBLE_INT
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_DOUBLE_LONG\t");
#	ifdef HAVE_POW_DOUBLE_LONG
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_FLOAT_INT\t");
#	ifdef HAVE_POW_FLOAT_INT
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_FLOAT_LONG\t");
#	ifdef HAVE_POW_FLOAT_LONG
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_FLOAT_DOUBLE\t");
#	ifdef HAVE_POW_FLOAT_DOUBLE
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n HAVE_POW_DOUBLE_FLOAT\t");
#	ifdef HAVE_POW_DOUBLE_FLOAT
	fprintf( ioQQQ,"defined");
#	endif

	fprintf( ioQQQ , "\n");

}

STATIC void setIsoNelemFlag(bool t_iso_sp::* cpt, long ipISO, long nelem)
{

	/* enable all iso-sequences for ipISO < 0  */
	long ipISOa=ipISO, ipISOb=ipISO+1;
	if ( ipISO<0 ) {
	ipISOa=ipH_LIKE;
	ipISOb=NISO;
	}
	/* enable all elements for nelem < 0 */
	long nelema=nelem, nelemb=nelem+1;
	if(nelem >=0 && nelem < ipISO )
	{
		fprintf(ioQQQ,"This iso-sequence (%s) and element (%s) are impossible.\n",
				elementnames.chElementName[ipISO],
				elementnames.chElementName[nelem]);
		cdEXIT(EXIT_FAILURE);
	}

	if( nelem<0 )
	{
	// element not specified, do entire sequence
	nelema=ipHYDROGEN;
	nelemb=LIMELM;
	}
	for (long ipISO1=ipISOa; ipISO1 < ipISOb; ++ipISO1)
	{
			for( long nelem1=max(ipISO1, nelema); nelem1 < nelemb; ++nelem1 )
				iso_sp[ipISO1][nelem1].*cpt = true;
	}

}

void ParsePrint(
	/* input line was converted to caps by calling routine */
	Parser &p )
{
	int ipISO;
	long int
	  j, 
	  nelem,
	  num1;
	double a;

	DEBUG_ENTRY( "ParsePrint()" );

	/* >>chng 01 aug 91, had been series of if branches, and could hit more than
	 * one key - dangerous!  changed to else if so only one hit per line possible */
	if( p.nMatch("AGES") )
	{
		/* print all estimates of cloud timescales */
		prt.lgPrnAges = true;
	}

	else if( p.nMatch("ARRA") )
	{
		/* print arrays for ionization balance of heavy elements */
		if( p.nMatch( "IONI" ) && p.nMatch( "LEVE" ) )
		{
			fprintf( ioQQQ, "Please use either IONIZATION or LEVEL keywords\n" );
			cdEXIT( EXIT_FAILURE );
		}
		else if( p.nMatch( "IONI" ) )
		{
			if( p.nMatch( "ONLY"  ) )
			{
				/* returns element number on C scale */
				if( (nelem = p.GetElem())<0 )
				{
					fprintf(ioQQQ,"An element name must appear on this PRINT ARRAYS ONLY xx command.\n");
					cdEXIT(EXIT_FAILURE);
				}
				/* have the element number, turn on its print */
				prt.lgPrtArry[nelem] = true;
			}
			else
			{
				/* this flag, print arrays for all elements */
				for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
				{
					prt.lgPrtArry[nelem] = true;
				}
			}
		}
		else if( p.nMatch( "LEVEL" ) )
		{
			string species;
			p.GetQuote( species );
			if( species == "" )
			{
				fprintf( ioQQQ, "Empty species label: ''\n" );
				cdEXIT( EXIT_FAILURE );
			}

			prt.matrix.setSpecies( species );
		}
		else
		{
			fprintf( ioQQQ, "Please use either the IONIZATION or the LEVEL keywords\n" );
			cdEXIT( EXIT_FAILURE );
		}
	}

	else if( p.nMatch("CITA") )
	{
		prt.lgPrtCitations = true;
	}

	else if( p.nMatch("COLU") && p.nMatch("DENS") )
	{
		if( p.nMatch(" ON ") )
		{
			/* print column densities of elements - this is default */
			prt.lgPrintColumns = true;
		}
		else if( p.nMatch(" OFF") )
		{
			/* print column densities of elements */
			prt.lgPrintColumns = false;
		}
	}

	/* option to print departure coefficients in addition to level pops 
	 * keywords He-like to do He-like sequence element, else do h-like
	 * element name, if not recognized, does hydrogen
	 * so with no options prints hydrogen itself */
	else if( p.nMatch("DEPA") )
	{
		if( p.nMatch("HE-L") )
		{
			ipISO = ipHE_LIKE;
		}
		else if( p.nMatch("H-LI"))
		{
			ipISO = ipH_LIKE;
		}
		else
			ipISO = -1;

		/* now check for element name */
		nelem = p.GetElem( );

		/* set the iso-sequence and element print */
		setIsoNelemFlag(&t_iso_sp::lgPrtDepartCoef, ipISO, nelem);
	}

	else if( p.nMatch("CRIT") )
	{

		if( p.nMatch("HE-L") )
		{
			ipISO = ipHE_LIKE;
		}
		else if( p.nMatch("H-LI"))
		{
			ipISO = ipH_LIKE;
		}
		else
			ipISO = -1;

		/* now check for element name */
		nelem = p.GetElem( );

		/* set the iso-sequence and element print */
		setIsoNelemFlag(&t_iso_sp::lgPrtNCrit, ipISO, nelem);
	}

	else if( p.nMatch("CONS") )
	{
		/* print physical constants, routine is below */
		prt_constants();
	}

	else if( p.nMatch("ERRO") )
	{
		/* print errors to special window */
		lgPrnErr = true;
	}

	else if( p.nMatch("FIXI") )
	{
		broke.lgPrintFixits = true;
	}

	else if( p.nMatch("HEAT") )
	{
		/* print heat arrays */
		prt.lgPrintHeating = true;
	}

	else if( p.nMatch("MODU") )
	{
		/* list loaded modules */
		fprintf(ioQQQ,"Loaded Cloudy modules are:\n");
		fprintf(ioQQQ,"--------------------------\n");
		vector<module*>& mods=module_list::Inst().m_l;
		for (vector<module*>::iterator mi = mods.begin(); mi != mods.end(); ++mi)
		{
			fprintf(ioQQQ,"%s\n",(*mi)->chName());
		}
		fprintf(ioQQQ,"--------------------------\n");
	}
	
	else if( p.nMatch("PATH") )
	{
		/* print the path */
		cpu.i().printDataPath();
	}

	/*else if( p.nMatch("H-LI"))*/
	else if( p.nMatch("POPU"))
	{
		if( p.nMatch("HE-L") )
		{
			ipISO = ipHE_LIKE;
		}
		else if( p.nMatch("H-LI"))
		{
			ipISO = ipH_LIKE;
		}
		else
			ipISO = -1;

		/* now check for element name */
		nelem = p.GetElem( );

		/* set the iso-sequence and element print */
		setIsoNelemFlag(&t_iso_sp::lgPrtLevelPops,ipISO,nelem);

	}

	/* option to only print last iteration */
	else if( p.nMatch("LAST") )
	{
		prt.lgPrtLastIt = true;
	}

	/* the print line command as several options */
	else if( p.nMatch("LINE") )
	{
		if( p.nMatch(" ALL") )
		{
			/* turn on all printed components */
			prt.lgPrnPump = true;
			prt.lgPrnColl = true;
			prt.lgPrnHeat = true;
		}

		else if( p.nMatch("CELL") )
		{
			/* print line cell on physics scale, first cell in continuum is 1
			 * give all lines in this cell */
			prt.lgPrnLineCell = true;
			prt.nPrnLineCell = (long)p.FFmtRead();
			if( p.lgEOL() )
				p.NoNumb("cell for line print" );
			if( prt.nPrnLineCell < 1 )
			{
				/* non-positive cells are not allowed */
				fprintf(ioQQQ , "The cell number on the PRINT LINE CELL command must be positive.\n");
				fprintf(ioQQQ , "The cell number was %li.\n" , prt.nPrnLineCell);
			}
		}

		else if( p.nMatch("COLL") )
		{
			/* either print line collisions or print line iso collapsed */
			/* also print collisional contributions */
			if( p.nMatch(" ISO") )
			{
				if( p.nMatch(" OFF") )
					prt.lgPrnIsoCollapsed = false;
				else
					/* print predictions from collapsed levels of iso sequences */
					prt.lgPrnIsoCollapsed = true;
			}
			else
			{
				/* print line collisions */
				prt.lgPrnColl = true;
			}
		}

		else if( p.nMatch("COLU") )
		{
			/* option to print main line array as a single column */
			prt.lgPrtLineArray = false;
			/* this also has an option - LINEAR - to print linear quantity 
			 * in exponential format */
			if( p.nMatch("LINEAR") )
				prt.lgPrtLineLog = false;
		}

		else if( p.nMatch("FAIN") && !(p.nMatch("OPTI")&&p.nMatch("DEPT")) )
		{
			/* print line faint - above do not trigger on optical depth 
			 * option to adjust intensity of faintest line to print */
			/* >> 01 feb 01, move print faint into print line faint */
			/* faintest line, rel to norm line, to print; either linear of log */
			a = p.FFmtRead();

			/* option for, if no number, keyword=" OFF", to print all lines */
			if( p.lgEOL() )
			{
				if( p.nMatch(" OFF") )
				{
					prt.lgFaintOn = false;
				}
				else
				{
					fprintf( ioQQQ, 
						" There faintest line to print must be on this line, sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}

			prt.lgFntSet = true;
			if( a <= 0. || p.nMatch(" LOG") )
			{
				prt.TooFaint = (realnum)exp10(a);
			}
			else
			{
				prt.TooFaint = (realnum)a;
			}
		}

		else if( p.nMatch("FLUX") && p.nMatch("EART"))
		{
			/* print line flux seen at earth */
			prt.lgPrintFluxEarth = true;
		}

		else if( p.nMatch(" H2") && p.nMatch("ELEC") )
		{
			/* print H2 electronic lines too - -1 since number of electronic
			 * levels is not yet known, will set when H2 actually called */
			h2.nElecLevelOutput = -1;
		}

		else if( p.nMatch("HEAT") )
		{
			/* also print heating contributions */
			prt.lgPrnHeat = true;
		}

		else if( p.nMatch("INWA") )
		{
			/* also print inward contributions */
			prt.lgPrnInwd = true;
		}

		else if( p.nMatch(" OFF") )
		{
			if( p.nMatch( "INTR" ) )
				prt.lgPrintBlockIntrinsic = false;
			else if( p.nMatch( "EMER" ) )
				prt.lgPrintBlockEmergent = false;
			else
				prt.lgPrintBlock = false;
		}

		else if( p.nMatch("OPTI") && p.nMatch("DEPT") )
		{
			/* print line optical depths, with option for smallest to print */
			if( p.nMatch(" OFF") )
			{
				/* turn off or on printing of optical depths - default off */
				prt.lgPrtTau = false;
			}
			else
			{
				prt.lgPrtTau = true;
			}
			if( p.nMatch("FAIN") )
			{
				/* log of faintest optical depth, default is linear value of 0.1 */
				prt.PrtTauFnt = (realnum)exp10(p.FFmtRead());
				if( p.lgEOL() )
				{
					fprintf( ioQQQ, " There must be a number for the FAINT option.  They are HEAD and ZONE.  Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
		}

		else if( p.nMatch("PRECISION"))
		{
			/* print line precision (number of significant figures)
			 * this affects all aspects of reading and writing lines */
			LineSave.sig_figs = (int) p.FFmtRead();
			if (p.lgEOL())
			{
				p.NoNumb("line precision");
			}
			else if (LineSave.sig_figs > LineSave.sig_figs_max)
			{
				fprintf(ioQQQ,
						" print line precision currently only works for up to %ld significant figures.\n",
						LineSave.sig_figs_max);
				cdEXIT(EXIT_FAILURE);
			}
		}

		else if( p.nMatch("PUMP") )
		{
			/* also print pump contributions */
			prt.lgPrnPump = true;
		}

		else if( p.nMatch("SORT") )
		{
			/* >>chng 01 aug 18, print sort command works after all these years,
			 * sort by wavelength or intensity */
			/* turn on sorting with respect to wavelength */
			prt.lgSortLines = true;
			if( p.nMatch("WAVE") )
			{
				/* sort by wavelength */
				/* remember which one to do */
				prt.lgSortLineIntensity = false;
				prt.lgSortLineWavelength = true;

				/* wavelength has range option */
				/* option to only certain print range of lines */
				if( p.nMatch("RANG") )
				{
					prt.wlSort1 = (realnum)p.getWave();

					prt.wlSort2 = (realnum)p.getWave();

					if( p.lgEOL() )
					{
						fprintf( ioQQQ, " There must be two numbers for the RANGE option, the lower and upper wavelength.  Sorry.\n" );
						cdEXIT(EXIT_FAILURE);
					}
					if( prt.wlSort1 <0. || prt.wlSort2 <0. || 
						prt.wlSort1 >= prt.wlSort2 )
					{
						fprintf( ioQQQ, " The lower and upper wavelength must be positive and in the correct order.  Sorry.\n" );
						cdEXIT(EXIT_FAILURE);
					}
				}
				else
				{
					prt.wlSort1 = -1;
					prt.wlSort2 = 1e30f;
				}
			}
			else if( p.nMatch("INTE") )
			{
				/* sort by intensity/luminosity */
				/* remember which one to do */
				prt.lgSortLineIntensity = true;
				prt.lgSortLineWavelength = false;
			}
			else
			{
				fprintf( ioQQQ, "I can sort by wavelength or intensity - one must be specified.\nSorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		else if( p.nMatch(" SUM") )
		{
			/* option to read in set of lines to sum over */
			ParsePrtLineSum( p );
		}

		else if( p.nMatch("SURF") && p.nMatch("BRIG") )
		{
			/* print surface brightness rather than 4pi J */
			prt.lgSurfaceBrightness = true;
			/* default is per sr, arcsec option changes to sq arcsec */
			if( p.nMatch("ARCS" ) )
			{
				/* use sr */
				prt.lgSurfaceBrightness_SR = false;
			}
			else
			{
				/* use sq arcsec */
				prt.lgSurfaceBrightness_SR = true;
			}
		}

		else if( p.nMatch("CUMU") )
		{
			/* print lines cumulative - integral of line emission over time */
			prt.lgPrintLineCumulative = true;
		}

		else if( p.nMatch("VACUUM") )
		{
			/* print lines vacuum - use vacuum wavelengths */
			prt.lgPrintLineAirWavelengths = false;
		}

		else
		{
			fprintf( ioQQQ, "One of the keys should have appeared.  \nPlease consult Hazy.\nSorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* print maser lines when TAV is called */
	else if( p.nMatch("MASE") )
	{
		prt.lgPrtMaser = true;
	}

	else if( p.nMatch("ONLY") )
	{
		if( p.nMatch("ZONE") )
			prt.lgOnlyZone = true;

		else if( p.nMatch("HEAD") )
			prt.lgOnlyHead = true;

		else
		{
			fprintf( ioQQQ, " There must be a keyword for the ONLY option.  They are HEAD and ZONE.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("STAR") )
	{
		/* start printout at specified zone */
		called.lgTalk = false;
		prt.lgPrtStart = true;
		prt.nstart = (long int)p.FFmtRead();
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, 
				" The zone on which the print is to start MUST be entered on this line.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* print continuum command */
	else if( p.nMatch("CONT") )
	{
		if( p.nMatch("BLOC") )
		{
			/* option to print emergent continuum at end of calculation*/
			fprintf(ioQQQ , " PROBLEM The PRINT CONTINUUM BLOCK command has been removed.  Ignored for now.\n");
		}
		else if( p.nMatch("INDI"  ))
		{
			/* option to print lines and continuum that go into each continuum 
			 * index the continuum index is the cell within the continuum 
			 * array - this identifies lines that occur within each 
			 * continuum cell */
			prt.lgPrtContIndices = true;
			/* these are lower and upper limits to the energy range in Rydbergs.
			* they are the first and second number on the command line, lower and
			* upper bounds of the code are used if not specified */
			/* if no number on line then zero is returned, this is fine, since
			 * we want the lower energy bound of the code */
			prt.lgPrtContIndices_lo_E = (realnum)p.FFmtRead();
			prt.lgPrtContIndices_hi_E = (realnum)p.FFmtRead();
			/* if we hit end of line then use high-energy limit of code - that is,
			 * include all energies */
			if( p.lgEOL() )
				prt.lgPrtContIndices_hi_E = (realnum)rfield.egamry();
		}
		else
		{
			/* option to print continuum points within emission lines block */
			fprintf( ioQQQ, " PROBLEM PRINT CONTINUUM command is now the default, and the command has been removed.\n" );
		}
	}

	else if( p.nMatch("COOL") )
	{
		/* print cooling array for a specified one */
		prt.nzdump = (long int)p.FFmtRead();

		/* dump all zones if argument is zero or not present  */
		if( p.lgEOL() )
		{
			prt.nzdump = 0;
		}
	}

	else if( p.nMatch("QUIE") || (p.nMatch(" OFF") && 
		!p.nMatch("FAIN" )) )
	{
		/* in above, there is also a 'print faint off' command
		 * QUIET or OFF means turn off printout */
		called.lgTalk = false;
		input.lgVisibilityStatus = false;
	}

	else if( p.nMatch("MACR") )
	{
		// print status of macros in cddefines.ht */
		PrtMacros();
	}

	else if( p.nMatch(" ON ") )
	{
		/* on means turn on printout, lgTalkIsOK is set false in grid_do.cpp.
		 * this keeps printout quiet during optimize, even when init files are parsed */
		/* called.lgTalkForcedOff was set true with cdTalk(false), if this was
		 * set then do not respect this command.  this is to prevent print on at end
		 * of init file from turning on print in grids when print is turned off */
		if( called.lgTalkIsOK && !called.lgTalkForcedOff )
		{
			called.lgTalk = cpu.i().lgMPI_talk();
		}
		input.lgVisibilityStatus = true;
	}

	else if (p.nMatch("RECOMB"))
	{
		/* option to print recombination rates then exit */
		ionbal.lgRecom_Badnell_print = true;
	}

	else if( p.nMatch("SHOR") )
	{
		/* make short printout, don't print last */
		prt.lgPrtShort = true;
		if( !prt.lgFntSet )
			prt.TooFaint = 0.001f;
	}

	else if( p.nMatch( " UTA" ) && p.nMatch( "REFE" ) )
	{
		atmdat.lgUTAprint = true;
	}

	else if( p.nMatch("VERS") )
	{
		/* print compiler and code version information */
		fprintf( ioQQQ, "\nThis is Cloudy %s\n%s\n\n" ,
			 t_version::Inst().chVersion.c_str(),
			 t_version::Inst().chInfo.c_str() );
	}

	else if( p.nMatch("ZONE") || p.nMatch("EVER") )
	{
		/* print every nth zone - command was originally print every but
		 * is being changed to print zone */
		num1 = (long int)p.FFmtRead();
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " The number of zones to print MUST be entered on this line.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		iterations.IterPrnt[0] = MAX2(num1,1);

		for( j=1; j < iterations.iter_alloc; j++ )
		{
			iterations.IterPrnt[j] = (long int)p.FFmtRead();
			if( p.lgEOL() )
			{
				iterations.IterPrnt[j] = iterations.IterPrnt[j-1];
			}
		}
	}

	/* check if no keywords were recognized. */
	else
	{
		fprintf( ioQQQ, " There MUST be a keyword on the following line.  Sorry.\n" );
		fprintf( ioQQQ, " The PRINT FAINT command is now the PRINT LINE FAINT command.\n" );
		p.PrintLine(ioQQQ);
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

/*prt_constants print physical and machine constants */
STATIC void prt_constants(void)
{
	DEBUG_ENTRY( "prt_constants()" );

	fprintf(ioQQQ,"\n\nPhysical constants used by Cloudy, taken from physconst.h\n");

	prt_phys_constants(ioQQQ);
	fprintf(ioQQQ,"\n");

	fprintf(ioQQQ,"Some other interesting sizes:\n");
	fprintf(ioQQQ,"bool\t%lu\n",(unsigned long)sizeof(bool));
	fprintf(ioQQQ,"char\t%lu\n",(unsigned long)sizeof(char));
	fprintf(ioQQQ,"int\t%lu\n",(unsigned long)sizeof(int));
	fprintf(ioQQQ,"long int\t%lu\n",(unsigned long)sizeof(long int));
	fprintf(ioQQQ,"unsigned int\t%lu\n",(unsigned long)sizeof(unsigned int));
	fprintf(ioQQQ,"float\t%lu\n",(unsigned long)sizeof(sys_float));
	fprintf(ioQQQ,"realnum\t%lu\n",(unsigned long)sizeof(realnum));
	fprintf(ioQQQ,"double\t%lu\n",(unsigned long)sizeof(double));
	fprintf(ioQQQ,"double*\t%lu\n",(unsigned long)sizeof(double*));
	fprintf(ioQQQ,"\n");

	fprintf(ioQQQ,"Some float constants (from float.h and cpu.h):\n");
	/* some constants from float.h */
	fprintf(ioQQQ,"DBL_DIG \t%i\n", DBL_DIG);         /* # of decimal digits of precision */
	fprintf(ioQQQ,"DBL_EPSILON \t%.15g\n",DBL_EPSILON);   /* smallest such that 1.0+DBL_EPSILON != 1.0 */
	fprintf(ioQQQ,"DBL_MANT_DIG\t%i\n",DBL_MANT_DIG); /* # of bits in mantissa */
	fprintf(ioQQQ,"DBL_MAX\t%.15g\n", DBL_MAX);           /* max value */
	fprintf(ioQQQ,"DBL_MAX_10_EXP\t%i\n", DBL_MAX_10_EXP); /* max decimal exponent */
	fprintf(ioQQQ,"DBL_MAX_EXP\t%i\n", DBL_MAX_EXP);  /* max binary exponent */
	fprintf(ioQQQ,"DBL_MIN\t%.15g\n", DBL_MIN);           /* min positive value */

	fprintf(ioQQQ,"FLT_DIG\t%i\n", FLT_DIG);          /* # of decimal digits of precision */
	fprintf(ioQQQ,"FLT_EPSILON\t%.15g\n", FLT_EPSILON);   /* smallest such that 1.0+FLT_EPSILON != 1.0 */
	fprintf(ioQQQ,"FLT_MANT_DIG\t%i\n", FLT_MANT_DIG); /* # of bits in mantissa */
	fprintf(ioQQQ,"FLT_MAX\t%.15g\n", FLT_MAX);            /* max value */
	fprintf(ioQQQ,"FLT_MAX_10_EXP\t%i\n", FLT_MAX_10_EXP);/* max decimal exponent */
	fprintf(ioQQQ,"FLT_MAX_EXP\t%i\n", FLT_MAX_EXP);   /* max binary exponent */
	fprintf(ioQQQ,"FLT_MIN\t%.15g\n", FLT_MIN);            /* min positive value */
	fprintf(ioQQQ,"\n");

	/* some constants from cpu.h */
	fprintf(ioQQQ,"BIGFLOAT\t%.15g\n", BIGFLOAT);
	fprintf(ioQQQ,"SMALLFLOAT\t%.15g\n", SMALLFLOAT);
	fprintf(ioQQQ,"BIGDOUBLE\t%.15g\n", BIGDOUBLE);
	fprintf(ioQQQ,"SMALLDOUBLE\t%.15g\n", SMALLDOUBLE);
	fprintf(ioQQQ,"\n");


	fprintf(ioQQQ,"Some integer constants (from limits.h and cpu.h):\n");
	/* some constants from limits.h */
	fprintf(ioQQQ,"INT_MAX\t%i\n", INT_MAX);
	fprintf(ioQQQ,"INT_MIN\t%i\n", INT_MIN);
	fprintf(ioQQQ,"UINT_MAX\t%u\n", UINT_MAX);
	fprintf(ioQQQ,"LONG_MAX\t%li\n", LONG_MAX);
	fprintf(ioQQQ,"LONG_MIN\t%li\n", LONG_MIN);
	fprintf(ioQQQ,"ULONG_MAX\t%lu\n", ULONG_MAX);
	fprintf(ioQQQ,"LLONG_MAX\t%lli\n", LLONG_MAX);
	fprintf(ioQQQ,"LLONG_MIN\t%lli\n", LLONG_MIN);
	fprintf(ioQQQ,"ULLONG_MAX\t%llu\n", ULLONG_MAX);
	fprintf(ioQQQ,"\n");

	/* some constants from cpu.h */
	fprintf(ioQQQ,"INT16_MAX\t%i\n", INT16_MAX);
	fprintf(ioQQQ,"INT16_MIN\t%i\n", INT16_MIN);
	fprintf(ioQQQ,"UINT16_MAX\t%u\n", UINT16_MAX);
	fprintf(ioQQQ,"INT32_MAX\t%li\n", (long)INT32_MAX);
	fprintf(ioQQQ,"INT32_MIN\t%li\n", (long)INT32_MIN);
	fprintf(ioQQQ,"UINT32_MAX\t%lu\n", (long)UINT32_MAX);
	fprintf(ioQQQ,"INT64_MAX\t%lli\n", (long long) INT64_MAX);
	fprintf(ioQQQ,"INT64_MIN\t%lli\n", (long long) INT64_MIN);
	fprintf(ioQQQ,"UINT64_MAX\t%llu\n", (unsigned long long) UINT64_MAX);
	fprintf(ioQQQ,"\n");


	fprintf(ioQQQ,"\nThis machine has %ld threads.\n", cpu.i().nCPU() );

	return;
}
