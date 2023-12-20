/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseOptimize parse the optimize command lines */
/*GetOptColDen read observed column densities & errors for optimizer */
/*GetOptLineInt parse observed line intensities for optimization routines */
/*GetOptTemp read observed temperatures & errors for optimizer */
#include "cddefines.h"
#include "trace.h"
#include "optimize.h"
#include "prt.h"
#include "flux.h"
#include "predcont.h"
#include "parser.h"
#include "lines_service.h"
#include "lines.h"
#include "elementnames.h"

static const realnum DEFERR = 0.05f;

/*GetOptLineInt parse observed line intensities for optimization routines */
STATIC void GetOptLineInt(Parser &p );

/*GetOptColDen read observed column densities & errors for optimizer */
STATIC void GetOptColDen(Parser &p );

/*GetOptTemp read observed temperatures & errors for optimizer */
STATIC void GetOptTemp(Parser &p );

/* ParseOptimize - called from ParseCommands if OPTIMIZE command found */
void ParseOptimize(
	/* command line, which was changed to all caps in main parsing routine */
	Parser &p)
{
	DEBUG_ENTRY( "ParseOptimize()" );

	/* this must come first so that contents of filename do not trigger wrong command */
	if( p.nMatch("FILE") )
	{
		/* option to send final set of parameters to an input file 
		 * get name within double quotes, 
		 * and set to blanks in chCard and OrgCard */
		/* chCard is all caps at this point.  GetQuote will work with 
		 * original version of command line, which preserves case of
		 * characters.  Also removes string between quotes */
		/* GetQuote will abort if it fails, so we can ignore the return value */
		if (p.GetQuote( chOptimFileName ))
			p.StringError();
	}

	else if( p.nMatch("COLU") )
	{
		/* optimize column density */
		/* read column densities to match */
		GetOptColDen(p);

		optimize.lgOptimize = true;
	}

	else if( p.nMatch("CONT") && p.nMatch("FLUX") )
	{
		/* optimize continuum flux at arbitrary wavelength */
		double energy = p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("energy");
		const char* unit = p.StandardEnergyUnit();
		long ind = t_PredCont::Inst().add( energy, unit );
		Energy E( energy, unit );

		double flux = p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("flux");
		if( flux <= 0. || p.nMatch( " LOG" ) )
			flux = exp10(flux);
		Flux F( E, flux, p.StandardFluxUnit() );

		double relerr = p.FFmtRead();
		if( p.lgEOL() )
			relerr = DEFERR;
		/* value is an upper limit only, use negative error to flag */
		if( p.nMatch("<" ) )
			relerr = -relerr;

		optimize.ContIndex.push_back( ind );
		optimize.ContEner.push_back( E );
		optimize.ContNFnu.push_back( F );
		optimize.ContNFnuErr.push_back( chi2_type(relerr) );
		optimize.lgOptimize = true;
	}

	else if( p.nMatch("CONT") )
	{
		/* set flag saying that optimization should start from continue file */
		optimize.lgOptCont = true;
		optimize.lgOptimize = true;
	}

	else if( p.nMatch("DIAM") )
	{
		/* optimize angular diameter */
		optimize.optDiam = chi2_type( p.FFmtRead() );
		optimize.optDiamErr = chi2_type( p.FFmtRead() );
		if( p.lgEOL() )
			optimize.optDiamErr = chi2_type( DEFERR );
		/* value is an upper limit only, use negative error to flag */
		if( p.nMatch( "<" ) )
			optimize.optDiamErr = -optimize.optDiamErr;

		if( p.nMatch( " LOG" ) )
			optimize.optDiam = exp10(optimize.optDiam);

		optimize.lgOptDiam = true;
		/* angular diameter is default in arcsec, or in cm of keyword present */
		optimize.lgDiamInCM = ( p.nMatch( " CM " ) != 0 );
		optimize.lgOptimize = true;
	}

	else if( p.nMatch("INCR") )
	{
		/* scan off increments for the previously selected parameter */
		if( optimize.nparm > 0 )
		{
			/* also called during optimization process, ignore then */
			optimize.OptIncrm[optimize.nparm-1] = 
				(realnum)p.FFmtRead();
		}
		else
		{
			if( optimize.lgInitialParse && optimize.lgVaryOn )
			{
				fprintf( ioQQQ, " The OPTIMIZE INCREMENT command needs to follow the"
					 " command being varied.\n None was found. Bailing out...\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	else if( p.nMatch("LUMI") || p.nMatch("INTE") )
	{
		/* scan off intensity or luminosity of normalization line */
		optimize.optint = (realnum)p.FFmtRead();
		optimize.optier = (realnum)p.FFmtRead();
		if( p.lgEOL() )
			optimize.optier = DEFERR;

		// default intrinsic intensity, accept emergent
		optimize.nOptLum = p.nMatch("EMER") ? 1 : 0;

		/* set flag to say that intensity or luminosity of line set */
		optimize.lgOptLum = true;
		optimize.lgOptimize = true;
	}

	else if( p.nMatch("ITER") )
	{
		/* scan off number of iterations */
		optimize.nIterOptim = (long)p.FFmtRead();
	}

	else if( p.nMatch("LINE") )
	{
		/* read lines to match */
		GetOptLineInt(p);

		optimize.lgOptimize = true;
	}

	else if( p.nMatch("PHYM") )
	{
		/* use PvH's PHYMIR to optimize parameters */
		strcpy( optimize.chOptRtn, "PHYM" );
#		if defined(__unix) || defined(__APPLE__)
		// turn parallel mode off if explicitly requested or if we are not in default mode
		// the latter usually means that each rank runs its own model in single-CPU mode
		optimize.lgParallel = ! ( p.nMatch("SEQU") || cpu.i().MPIMode() != MS_DEFAULT );
#		else
		optimize.lgParallel = false;
#		endif
		if( optimize.lgParallel ) 
		{
			long dum = (long)p.FFmtRead();
			/* default is the total number of available CPUs */
			optimize.useCPU = p.lgEOL() ? cpu.i().nCPU() : dum;
		}
		else 
		{
			optimize.useCPU = 1;
		}
	}

	else if( p.nMatch("RANG") )
	{
		/* scan off range for the previously selected variable */
		if( optimize.nparm > 0 )
		{
			bool lgFirstOneReal = false;

			optimize.varang[optimize.nparm-1][0] = 
				(realnum)p.FFmtRead();
			if( p.lgEOL() )
			{
				optimize.varang[optimize.nparm-1][0] = -FLT_MAX;
			}
			else
			{
				lgFirstOneReal = true;
			}

			optimize.varang[optimize.nparm-1][1] = 
				(realnum)p.FFmtRead();
			if( p.lgEOL() )
			{
				optimize.varang[optimize.nparm-1][1] = FLT_MAX;
			}
			else if( lgFirstOneReal )
			{
				/* >>chng 06 aug 22, swap if second range parameter is less than the first,
				 * and always increment optimize.nRangeSet */
				++optimize.nRangeSet;
				if( optimize.varang[optimize.nparm-1][1] < optimize.varang[optimize.nparm-1][0] )
				{
					realnum temp = optimize.varang[optimize.nparm-1][0];
					optimize.varang[optimize.nparm-1][0] = optimize.varang[optimize.nparm-1][1];
					optimize.varang[optimize.nparm-1][1] = temp;
				}
			}
		}
	}

	else if( p.nMatch("SUBP") )
	{
		/* use subplex to optimize parameters */
		strcpy( optimize.chOptRtn, "SUBP" );
	}

	/* match a temperature */
	else if( p.nMatch("TEMP") )
	{
		/* read temperatures to match */
		GetOptTemp(p);

		optimize.lgOptimize = true;
	}

	else if( p.nMatch("TOLE") )
	{
		/* scan off tolerance of fit, sum of residuals must be smaller than this
		 * default is 0.10 */
		optimize.OptGlobalErr = (realnum)p.FFmtRead();
	}

	else if( p.nMatch("TRAC") )
	{
		if( p.nMatch("STAR") )
		{
			/* trace start iteration number
			 * turn on trace printout starting on nth call to cloudy */
			optimize.nTrOpt = (long)p.FFmtRead();
			if( p.lgEOL() )
			{
				fprintf( ioQQQ, " optimize trace start command:\n" );
				fprintf( ioQQQ, " The iteration number must appear.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			optimize.lgTrOpt = true;
		}
		else if( p.nMatch("FLOW") )
		{
			/* trace flow
			 * follow logical flow within code */
			optimize.lgOptimFlow = true;
		}
		else
		{
			fprintf( ioQQQ, " optimize trace flow command:\n" );
			fprintf( ioQQQ, " One of the sub keys START or FLOW must appear.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else
	{
		p.PrintLine(ioQQQ);
		fprintf( ioQQQ, " is unrecognized keyword, consult HAZY.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

/*GetOptColDen read observed column densities & errors for optimizer */
STATIC void GetOptColDen(Parser &p )
{

	DEBUG_ENTRY( "GetOptColDen()" );

	/* read observed column densities & errors */

	/* get first line */
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, " Hit EOF while reading column density list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	while( !p.m_lgEOF )
	{
		string chLabel;
		if ( p.GetQuote( chLabel ) == 0 )
		{
			// Read a quote-delimited string, which will be interpreted
			// by cdColm as a generic species label.  It would be better
			// to do the lookup earlier, to validate the input and save
			// repeated effort.
			
			/* order on line is species, column density, err */
			optimize.chColDen_label.push_back(chLabel);
			optimize.ion_ColDen.push_back(0);
		}
		else
		{
			/* order on line is element name (at least 4 characters), ionization stage, column density, err */
			/* copy cap'd version of first 4 char of chCard to chColDen_label */
			// The other branch can handle also ions, with the 'chemical'
			// syntax '"H+"' (in quotes), rather than 'HYDR 2' below
			optimize.chColDen_label.push_back(elementnames.chElementNameShort[p.getElement()]);
			
			/* now get the ion stage, this should be 1 for atom, up to element
			 * number plus one */
			long ion = nint(p.FFmtRead());
			if( p.lgEOL() )
			{
				p.PrintLine( ioQQQ );
				fprintf( ioQQQ, " The ionization stage MUST appear on this line.  Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			
			/* the ion must be 1 or greater unless requesting a special,
			 * like a molecule or excited state population, in which
			 * case ion = 0
			 * can't check on upper limit yet since have not resolved element name */
			if( ion < 0 )
			{
				p.PrintLine( ioQQQ );
				fprintf( ioQQQ, " An ionization stage of %ld does not make sense.  Sorry.\n", ion );
				cdEXIT(EXIT_FAILURE);
			}
			optimize.ion_ColDen.push_back(ion);
		}

		optimize.ColDen_Obs.push_back( realnum(exp10(p.FFmtRead())) );
		if( p.lgEOL() )
		{
			p.PrintLine( ioQQQ );
			fprintf( ioQQQ, " An observed column density MUST be entered.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		realnum err = realnum(p.FFmtRead());
		if( err <= 0.0 )
		{
			/* this is the relative error allowed */
			err = realnum(DEFERR);
		}

		/* check if number is a limit - if '<' appears on the line then it is */
		if( p.nMatch( "<" ) )
		{
			/* value is an upper limit only, use negative error to flag */
			err = -err;
		}
		optimize.ColDen_error.push_back(err);

		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading column density list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( p.hasCommand( "END" ) )
		{
			p.m_lgEOF = true;
		}
	}

	if( trace.lgTrace && optimize.lgTrOpt )
	{
		fprintf( ioQQQ, "%ld columns were entered, they were;\n", 
					(long int) optimize.ColDen_Obs.size() );
		for( long i=0; i < long(optimize.ColDen_Obs.size()); i++ )
		{
			fprintf( ioQQQ, " %s ion=%ld %.2e %.2e\n", 
				 optimize.chColDen_label[i].c_str(), optimize.ion_ColDen[i],
				 optimize.ColDen_Obs[i], optimize.ColDen_error[i] );
		}
	}
	return;
}

/*GetOptLineInt parse observed line intensities for optimization routines */
STATIC void GetOptLineInt(Parser &p)
{
	DEBUG_ENTRY( "GetOptLineInt()" );

	// default intrinsic intensity, accept emergent
	optimize.nEmergent = p.nMatch("EMER") ? 1 : 0;

	/* read observed line fluxes & errors */
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	while( !p.m_lgEOF )
	{
		LineID line = p.getLineID();
		/* order on line is label (col 1-4), wavelength, flux, error */
		optimize.lineids.emplace_back( line );

		/* get the error associated with specified significant figures */
		optimize.errorwave.emplace_back( WavlenErrorGet(line.wave, LineSave.sig_figs ) );

		/* next get the observed intensity */
		realnum xLineInt = realnum(p.FFmtRead());
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " The wavelength and relative intensity MUST be entered on this line.  Sorry.\n" );
			fprintf( ioQQQ, " The command line is the following:\n " );
			p.PrintLine(ioQQQ);
			cdEXIT(EXIT_FAILURE);
		}

		if( xLineInt <= 0. )
		{
			fprintf( ioQQQ, " An observed intensity of %.2e is not allowed.  Sorry.\n", 
				 xLineInt );
			fprintf( ioQQQ, " The command line is the following:\n" );
			p.PrintLine(ioQQQ);
			cdEXIT(EXIT_FAILURE);
		}
		optimize.xLineInt_Obs.push_back(xLineInt);

		/* finally the optional error */
		realnum err = realnum(p.FFmtRead());
		/* most often will use the default */
		if( err <= 0.0 )
			err = realnum(DEFERR);

		/* check if number is a limit - if '<' appears on the line then it is */
		if( p.nMatch( "<" ) )
		{
			/* value is an upper limit only, use negative error to flag */
			err = -err;
		}
		optimize.xLineInt_error.push_back(err);

		if( !p.lgReachedEnd() )
		{
			fprintf( ioQQQ, "GetOptLineInt: found junk at end of input line:\n" );
			p.showLocation();
			cdEXIT(EXIT_FAILURE);
		}

		/* get next line */
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading line list for optimize command; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( p.hasCommand( "END" ) )
			p.m_lgEOF = true;
	}

	if( trace.lgTrace && trace.lgTrOptm )
	{
		fprintf( ioQQQ, "%ld lines were entered, they were:\n", 
			 (long int) optimize.xLineInt_Obs.size() );

		for( long i=0; i < long(optimize.xLineInt_Obs.size()); i++ )
		{
			fprintf( ioQQQ, " %s ", optimize.lineids[i].chLabel.c_str() );
			prt_wl( ioQQQ, optimize.lineids[i].wave );

			fprintf( ioQQQ, " %10.2e%10.2e\n", 
				optimize.xLineInt_Obs[i], 
				optimize.xLineInt_error[i] );
		}
	}
	return;
}

/*GetOptTemp parse observed line intensities for optimization routines */
STATIC void GetOptTemp(Parser &p )
{
	DEBUG_ENTRY( "GetOptTemp()" );

	/* read observed line fluxes & errors - first set total number of observe temps */
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, " Hit EOF while reading line list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* make line caps so we can parse it */
	while( !p.m_lgEOF )
	{
		/* order on line is label (col 1-4), ion, temperature, error */
		optimize.chTempLab.push_back( elementnames.chElementNameShort[p.getElement()] );

		/* next get the ion stage */
		optimize.ionTemp.push_back( nint(p.FFmtRead()) );

		/* next get the observed temperature */
		realnum temp_obs = realnum(p.FFmtRead());
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " The ion stage and temperature MUST be entered on this line.  Sorry.\n" );
			fprintf( ioQQQ, " The command line is the following:\n " );
			p.PrintLine(ioQQQ);
			cdEXIT(EXIT_FAILURE);
		}
		/* temperatures less than or equal to 10 are logs */
		if( temp_obs <= 10. )
			temp_obs = realnum(exp10(  (double)temp_obs ) );
		optimize.temp_obs.push_back(temp_obs);

		/* finally the optional error */
		realnum temp_error = realnum(p.FFmtRead());
		/* most often will use the default */
		if( temp_error <= 0.f )
			temp_error = realnum(DEFERR);
		/* check if number is a limit - if '<' appears on the line then it is */
		if( p.nMatch( "<" ) )
			temp_error = -temp_error;
		optimize.temp_error.push_back(temp_error);

		/* check for radius or volume for how to weight the mean temp 
		 * this will be the default */
		/* >>chng 05 dec 29, from chCard to chCap, unlike much of code in this file,
		 * chCard has been read in by this routine and contains the original form of 
		 * the command line.  It was converted to caps and stored in chCap above.
		 * As it was written only VOLUME was matched, would not match volume.
		 * Bug caught and corrected by Bohdan Melekh */
		if( p.nMatch( "VOLUME" ) )
			optimize.chTempWeight.push_back("volume");
		else
			optimize.chTempWeight.push_back("radius");

		/* get next line */
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, " Hit EOF while reading line list for optimize command; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( p.hasCommand( "END" ) )
			p.m_lgEOF = true;
	}

	if( trace.lgTrace && trace.lgTrOptm )
	{
		fprintf( ioQQQ, "%ld temperatures were entered, they were;\n", 
			 (long int) optimize.temp_obs.size() );

		for( long i=0; i < long(optimize.temp_obs.size()); i++ )
		{
			fprintf( ioQQQ, " %s ", optimize.chTempLab[i].c_str() );
			fprintf( ioQQQ, " %li " , optimize.ionTemp[i] );

			fprintf( ioQQQ, " %.2e %.2e\n", 
				optimize.temp_obs[i], 
				optimize.temp_error[i] );
		}
	}
	return;
}
