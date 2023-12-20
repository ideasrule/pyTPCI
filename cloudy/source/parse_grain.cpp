/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseGrain parse parameters on grains command */
#include "cddefines.h"
#include "grainvar.h"
#include "phycon.h"
#include "input.h"
#include "optimize.h"
#include "parser.h"
#include "grains.h"

void ParseGrain(Parser &p)
{
	bool lgC15 = false,
	  lgC120 = false,
	  lgLinSet, 
	  lgLogLinSet, 
	  lgQuoteFound,
	  lgSizeDistribution;
	GrainPar gp;
	
	/*possible name of input file with opacities */
	string chFile;
	const char *chOption = NULL;

	DEBUG_ENTRY( "ParseGrain()" );

	p.m_lgDSet = true;

	/* >>chng 00 dec 20, initialize chFile to empty string..., PvH */
	chFile = "";
	/* >>chng 01 jan 17, read filename first, it may contain digits that would upset FFmtRead,
	 * as well as other tests that will follow. GetQuote will erase the filename from chCard, PvH */
	lgQuoteFound = p.nMatch("\"");
	if( lgQuoteFound )
	{
		/* this will both scan in whatever label is inside the quotes in OrgCard, 
		 * but also remove the contents there and in parser buffer,
		 * so that following keywords will not trigger off it */
		if (p.GetQuote( chFile ))
			p.StringError();
	}

	if( p.nMatch("GREY") || p.nMatch("GRAY") || nMatch("grey_", chFile.c_str()) )
		gp.lgGreyGrain = true;
	else
		gp.lgGreyGrain = false;

	/* now check for the keyword "function" - this sets the behavior of the grain abundance with depth */
	if( p.nMatch("FUNC") )
	{
		if( p.nMatch("SUBL") )
			gp.nDustFunc = DF_SUBLIMATION;
		else
			gp.nDustFunc = DF_USER_FUNCTION;
	}
	else
		gp.nDustFunc = DF_STANDARD;

	/* check for the keyword "single bin" - 
	 * use resolved grain size distributions if not present 
	 * NB - logic is backwards here */
	if( !p.nMatch("SING") )
		lgSizeDistribution = true;
	else
		lgSizeDistribution = false;

	/* check for the keyword "qheating" - 
	 * quantum heating should be turned on */
	gp.lgForbidQHeating = false;

	if( p.nMatch("QHEA") )
		gp.lgRequestQHeating = true;
	else
		gp.lgRequestQHeating = false;

	/* the keyword "no qheat" always takes precedence */
	if( p.nMatch("NO QH") ) 
	{
		gp.lgForbidQHeating = true;
		gp.lgRequestQHeating = false;
		phycon.lgPhysOK = false;
	}

	/* option to force constant reevaluation of grain physics - 
	 * usually reevaluate grains at all times, but NO REEVALUATE will
	 * save some time but may affect stability */
	gv.lgReevaluate = !p.nMatch(" NO REEV");

	/* option to turn off photoelectric heating by grain, NO HEATING */
	if( p.nMatch("NO HE") )
	{
		phycon.lgPhysOK = false;
		gv.lgDHetOn = false;
	}

	/* option to turn off gas cooling by grain, NO COOLING */
	if( p.nMatch("NO CO") )
	{
		phycon.lgPhysOK = false;
		gv.lgDColOn = false;
	}

	/* these are keywords for PAH's, they need to be read before the depletion factor */
	lgC120 = p.nMatchErase("C120 ");
	if (!lgC120)
	{
		lgC15 = p.nMatchErase("C15 ");
	}

	/* log - linear option for grain abundance */
	lgLogLinSet = false;
	lgLinSet = false;
	if( p.nMatch(" LOG") )
	{
		lgLogLinSet = true;
		lgLinSet = false;
	}
	else if( p.nMatch("LINE") )
	{
		lgLogLinSet = true;
		lgLinSet = true;
	}

	/* get the grain abundance as the first parameter,
	 * returns 0 if no number, ok since interpreted as log, or unity*/
	gp.dep = p.FFmtRead();

	/* was keyword log or linear on the line? */
	if( lgLogLinSet )
	{
		/* log or linear was specified, which was it */
		if( lgLinSet )
		{
			/* linear quantity entered, check it */
			if( gp.dep <= 0. )
			{
				fprintf( ioQQQ, " Impossible value for linear abundance.\n" );
				fprintf( ioQQQ, " Abundance entered was%10.2e\n", gp.dep );
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			gp.dep = exp10(gp.dep);
		}
	}
	else
	{
		/* neither log nor linear specified, check sign
		 * force it to be a log - linear if greater than 0 */
		if( gp.dep <= 0. )
		{
			gp.dep = exp10(gp.dep);
		}
	}

	if( gp.dep < FLT_MIN )
	{
		fprintf( ioQQQ, " Grain abundance entered here (%f) is impossible.\n", gp.dep );
		cdEXIT(EXIT_FAILURE);
	}

	/* it is possible that there is a file name on the command line - 
	 * if so then we want to call correct reoutine, and not look for keywords 
	 * the signature of a keyword is a pair of quotes - is one present? */
	if( lgQuoteFound )
	{
		/* read the file name that was specified */
		chOption = "";
		mie_read_opc(chFile.c_str(),gp);
	}
	else
	{
		if( p.nMatch("ORIO") )
		{
			/* This scales the Orion grain abundance so that the observed 
			 * dust to gas ratio that Cloudy predicts is in agreement with
			 * that observed in the Veil,
			 *>>refer	grain	Abel, N., Brogan, C., Ferland, G., O'Dell, C.R., 
			 *>>refercon	Shaw, G., Troland, T., 2004, ApJ, submitted */
			gp.dep *= 0.85;

			/* optional keyword ORION to use orion curves for large R grains */
			/* only turn on one if graphite or silicate is on line, both if not */
			if( p.nMatch("GRAP") )
			{
				/* only turn on orion graphite */
				chOption = "ORION GRAPHITE ";
				if( lgSizeDistribution )
				{
					mie_read_opc("graphite_orion_10.opc",gp);
				}
				else
				{
					mie_read_opc("graphite_orion_01.opc",gp);
				}
			}
			else if( p.nMatch("SILI") )
			{
				/* only turn on orion silicate */
				chOption = "ORION SILICATE ";
				if( lgSizeDistribution )
				{
					mie_read_opc("silicate_orion_10.opc",gp);
				}
				else
				{
					mie_read_opc("silicate_orion_01.opc",gp);
				}
			}
			else
			{
				/* turn both on */
				chOption = "ORION ";
				if( lgSizeDistribution )
				{
					mie_read_opc("graphite_orion_10.opc",gp);
					mie_read_opc("silicate_orion_10.opc",gp);
				}
				else
				{
					mie_read_opc("graphite_orion_01.opc",gp);
					mie_read_opc("silicate_orion_01.opc",gp);
				}
			}
		}

		else if( p.nMatch(" PAH") )
		{
			/* only turn on the large PAH */
			if( lgC120 )
			{
				chOption = "PAH C120 ";
				mie_read_opc("pah1_c120.opc",gp);
			}
			/* only turn on the small PAH */
			else if( lgC15 )
			{
				chOption = "PAH C15 ";
				mie_read_opc("pah1_c15.opc",gp);
			}
			/* turn on size-distributed PAHs */
			else
			{
				/* the variable abundance for PAHs is hard wired in GrnStdDpth
				 * the function keyword has no effect?? */
				chOption = "PAH ";
				if( lgSizeDistribution )
				{
					mie_read_opc("pah1_ab08_10.opc",gp);
				}
				else
				{
					mie_read_opc("pah1_ab08_01.opc",gp);
				}
			}
		}

		else if( p.nMatch("GREY") || p.nMatch("GRAY") )
		{
			/* grey grains */
			chOption = "GREY ";
			if( lgSizeDistribution )
			{
				mie_read_opc("grey_ism_10.opc",gp);
			}
			else
			{
				mie_read_opc("grey_ism_01.opc",gp);
			}
		}

		else if( p.nMatch(" ISM") )
		{
			if( p.nMatch("GRAP") )
			{
				/* only turn on ism graphite */
				chOption = "ISM GRAPHITE ";
				if( lgSizeDistribution )
				{
					mie_read_opc("graphite_ism_10.opc",gp);
				}
				else
				{
					mie_read_opc("graphite_ism_01.opc",gp);
				}
			}
			else if( p.nMatch("SILI") )
			{
				/* only turn on orion silicate */
				chOption = "ISM SILICATE ";
				if( lgSizeDistribution )
				{
					mie_read_opc("silicate_ism_10.opc",gp);
				}
				else
				{
					mie_read_opc("silicate_ism_01.opc",gp);
				}
			}
			else
			{
				/* turn both ISM graphite and silicate on */
				chOption = "ISM ";
				if( lgSizeDistribution )
				{
					mie_read_opc("graphite_ism_10.opc",gp);
					mie_read_opc("silicate_ism_10.opc",gp);
				}
				else
				{
					mie_read_opc("graphite_ism_01.opc",gp);
					mie_read_opc("silicate_ism_01.opc",gp);
				}
			}
		}

		/* default case */
		else
		{
			/* turn both ISM graphite and silicate on */
			chOption = "";
			if( lgSizeDistribution )
			{
				mie_read_opc("graphite_ism_10.opc",gp);
				mie_read_opc("silicate_ism_10.opc",gp);
			}
			else
			{
				mie_read_opc("graphite_ism_01.opc",gp);
				mie_read_opc("silicate_ism_01.opc",gp);
			}
		}
	}

	/* vary option */
	if( optimize.lgVarOn )
	{
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vparm[0][optimize.nparm] = (realnum)log10(gp.dep);
		optimize.vincr[optimize.nparm] = 1.;

		// now build the input command string with all the necessary options...
		string command( "GRAIN ABUND=%f LOG " );
		if( chFile != "" )
		{
			command += "\"";
			command += chFile;
			command += "\" ";
		}

		// every branch of the if-statement above must set chOption.
		// this sentinel is here in case somebody forgets...
		if( chOption == NULL )
			TotalInsanity();

		command += chOption;

		/* generate keywords to feed to vary form of grains command 
		 * this test makes sure that future versions and values
		 * of nDustFunc are implemented here */
		switch( gp.nDustFunc )
		{
		case DF_STANDARD:
			// nothing to do
			break;
		case DF_USER_FUNCTION:
			command += "FUNCTION ";
			break;
		case DF_SUBLIMATION:
			command += "FUNCTION SUBLIMATION ";
			break;
		default:
			TotalInsanity();
		}

		if( !lgSizeDistribution )
			command += "SINGLE ";

		if( gp.lgForbidQHeating )
			command += "NO QHEAT ";
		else if( gp.lgRequestQHeating )
			command += "QHEAT ";

		if( !gv.lgReevaluate )
			command += "NO REEVALUATE ";
		if( !gv.lgDHetOn )
			command += "NO HEATING ";
		if( !gv.lgDColOn )
			command += "NO COOLING ";

		if( command.length() < static_cast<string::size_type>(FILENAME_PATH_LENGTH_2) )
			strcpy( optimize.chVarFmt[optimize.nparm], command.c_str() );
		else
		{
			fprintf(ioQQQ," grain command string is too long.  This is parse_grain\n");
			TotalInsanity();
		}

		optimize.nvarxt[optimize.nparm] = 1;
		++optimize.nparm;
	}
	return;
}
