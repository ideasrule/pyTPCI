/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseAbundances parse and read in composition as set by abundances command */
#include "cddefines.h"
#include "abund.h"
#include "dense.h"
#include "elementnames.h"
#include "input.h"
#include "parser.h"


STATIC FILE *open_abn_file( const string& chFile )
{
	string chPath = "abundances" + cpu.i().chDirSeparator() + chFile;
	return open_data( chPath, "r" );	// will abort if not found
}

void ParseAbundances(Parser &p)
			/* following set true by grains command,
			 * so this will not set more grains if grains already on. */
{
	bool lgLog;
	long int i;
	double absav[LIMELM];

	DEBUG_ENTRY( "ParseAbundances()" );

	/* abundances no longer set at reference value */
	abund.lgAbnReference = false;

	if( p.nMatch("STAR") )
	{
		/* Fred Hamann's star burst galaxy mixture -- includes a number which isn't an abundance */
		abund_starburst(p);
		return;
	}

	/* GetQuote should be above the other options so other commands don't trip
	 * if there is any overlapping text in the quotes */
	string chFile;	/*file name for table read */
	//records whether or not quotes were found and stores what is between them in chFile
	bool lgMatchFound = true;
	if( p.GetQuote( chFile ) )
		lgMatchFound = false;
	
	bool lgPrint = p.nMatch( "PRINT");
	bool lgIsotp = p.nMatch( "ISOT" );

	if(!lgMatchFound && !lgIsotp)
	{
		lgMatchFound = true;
		if(p.nMatch("CAME"))
			chFile = "Cameron.abn";
		else if (p.nMatch("CRAB"))
			chFile = "Crab.abn";
		else if (p.nMatch("GASS"))
			chFile = "solar_GASS10.abn";
		else if (p.nMatch("H II") || p.nMatch("HII ") || p.nMatch("ORIO"))
			chFile = "HII.abn";
		else if (p.nMatch("ISM"))
			chFile = "ISM.abn";
		else if (p.nMatch("NOVA"))
			chFile = "nova.abn";
		else if (p.nMatch(" AGB") || p.nMatch("AGB ") || p.nMatch("PLAN"))
			chFile = "PN.abn";
		else if (p.nMatch("PRIM"))
			chFile = "primordial.abn";
		else if (p.nMatch("OLD ") && p.nMatch("SOLA"))
			chFile = "solar84.abn";
		else if (p.nMatch("ALLE"))
			chFile = "allen73.abn";
		else
			lgMatchFound = false;
	}
	else if( lgIsotp && chFile.length() == 0 )
	{
		if( p.nMatch( "ASPL" ) )
			chFile = "Asplund09-iso.abn";
		else if( p.nMatch( "LODDERS03" ) )
			chFile = "Lodders03-iso.abn";
		else if( p.nMatch( "LODDERS09" ) )
			chFile = "Lodders09-iso.abn";
		else if( p.nMatch( "ROSM" ) )
			chFile = "Rosman98-iso.abn";
		else
		{
			fprintf(ioQQQ, "Unknown isotope abundances file:\t %s\n", chFile.c_str() );
			cdEXIT(EXIT_FAILURE);
		}
	}

	if(lgMatchFound && !lgIsotp)
	{
		// we have a file name - are other parameters present?
		bool lgGrainsON = true;
		if( p.nMatch("NO GR") != 0 )
			lgGrainsON = false;

		bool lgQHeat = true;
		if( p.nMatch("NO QH") != 0 )
			lgQHeat = false;

		abund.lgAbundancesSet = true;

		// initialization
		for( int nelem=1; nelem < LIMELM; nelem++ )
		{
			/* turn off all elements except Hydrogen,
			 * then turn on each element when found in the file */
			ostringstream chDUMMY;
			chDUMMY << "element " << elementnames.chElementName[nelem] << " off";
			p.setline(chDUMMY.str());
			// need to retain this flag, set by user to explicitly turn off element,
			// overriding our default to turn on elements that appear in abundances list
			bool lgSave = dense.lgElmtSetOff[nelem];
			ParseElement( p );
			dense.lgElmtSetOff[nelem] = lgSave;
		}

		FILE *ioDATA = open_abn_file( chFile );

		/* Sets abundance of Hydrogen to 1.0 in the expectation
		 * that other values are abundances relative to Hydrogen,
		 * although Hydrogen's abundance will still be set to
		 * whatever value is found in the code. This is just in
		 * case it is assumed */
		abund.ReferenceAbun[ipHYDROGEN] = 1.0;//The value found in abund.ReferenceAbun[0] before this point is found to be 1.0 as well

		string chLine;
		while( read_whole_line( chLine, ioDATA ) )
		{
			if( chLine.size() == 0 || chLine[0]=='\n' || chLine[0]=='\r' )
			{
				fprintf(ioQQQ, "PROBLEM in ABUNDANCES: Encountered unexpected empty line.\n");
				cdEXIT(EXIT_FAILURE);
			}

			if( chLine[0]=='*' )
				break;

			/* skip comment */
			if( chLine[0]=='#' )
				continue;

			size_t pp;
			/* erase EOL character */
			if( (pp = chLine.find_first_of("\n\r")) != string::npos )
				chLine.erase(pp);

			// If keyword "grains" is found on a line, calls grain command
			// NO QHEAT may have been on the line in the input deck, or in the abn file
			string chCAPS = chLine;
			caps(chCAPS);
			if( chCAPS.find("GRAINS") != string::npos )
			{
				if( lgPrint )
					fprintf(ioQQQ,"%s\n",chLine.c_str());

				//Makes sure grains have not already been set and are on
				//Either way it skips to the next loop iteration and does
				//not check the line for element
				if( !p.m_lgDSet && lgGrainsON )
				{
					if( !lgQHeat )
						chCAPS += " NO QHEAT";

					p.setline(chCAPS);
					ParseGrain(p);
					continue;
				}
				else
					continue;
			}

			bool lgFound = false;
			for(int nelem=0; nelem<LIMELM; nelem++)
			{
				if( chCAPS.find(elementnames.chElementNameShort[nelem]) != string::npos )
				{
					lgFound = true;
					i = 1;
					bool lgEOL;
					abund.ReferenceAbun[nelem] = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
					if( abund.ReferenceAbun[nelem] <= 0. )
					{
						fprintf(ioQQQ, "PROBLEM in ABUNDANCES: negative abundance not allowed.\n");
						fprintf(ioQQQ, "Non-positive abundance found on this line: %s\n", chLine.c_str());
						cdEXIT(EXIT_FAILURE);
					}

					// do not implicitly turn on an element that has explicitly been turned off
					if( dense.lgElmtSetOff[nelem] )
						break;

					/* turn on any element found while parsing the file */
					ostringstream chDUMMY;
					chDUMMY << "element " << elementnames.chElementName[nelem] << " on";
					p.setline(chDUMMY.str());
					ParseElement( p );

					//we shouldn't need to continue once an element name is found on the line...
					break;
				}
			}
			if( !lgFound )
			{
				fprintf(ioQQQ, "PROBLEM in ABUNDANCES: did not identify element name on this line: %s\n",
					chLine.c_str());
				cdEXIT(EXIT_FAILURE);
			}
		}

		//Normalizes all elements in abund.ReferenceAbun relative to the quantity of Hydrogen
		ASSERT( abund.ReferenceAbun[0]>0. );
		for(int nelem=0; nelem<LIMELM; nelem++)
		{
			abund.ReferenceAbun[nelem] /= abund.ReferenceAbun[ipHYDROGEN];
			if( lgPrint && dense.lgElmtOn[nelem] )
				fprintf(ioQQQ,"%s\t%.3e\t%.3f\n",elementnames.chElementName[nelem],
						abund.ReferenceAbun[nelem] , log10(SDIV(abund.ReferenceAbun[nelem])) );
		}
		fclose( ioDATA );
		return;
	}
	else if( !lgIsotp )
	{
		abund.lgAbundancesSet = true;

		if( p.nMatch(" ALL") )
		{
			/* special option, all abundances will be this number */

			/* Look for a number on the line */
			absav[0] = p.FFmtRead();
			if( p.lgEOL() )
				p.NoNumb("abundance");

			if( absav[0] <= 0. )
				absav[0] = exp10(absav[0]);

			for( i=1; i < LIMELM; i++ )
				abund.ReferenceAbun[i] = (realnum)absav[0];
			return;
		}

		/* regular option, supply He=-1 C=-3 etc... */
		long int ipSolar[LIMELM], npSolar = 0;
		while( true )
		{
			Symbol s = p.getSymbol();
			if( s.toktype == Symbol::EOSTAT )
			{
				// reached EOL, check if there is a continue line
				if( p.peekNextCommand("CONTINUE") )
				{
					p.getline();
					p.hasCommand("CONTINUE"); // skip to first symbol
					s = p.getSymbol();
				}
				else
					break;
			}
			if( s.toktype != Symbol::NAME )
			{
				fprintf( ioQQQ, "Expected an element symbol.\n" );
				p.showLocation();
				cdEXIT(EXIT_FAILURE);
			}
			ipSolar[npSolar] = elem_symbol_to_index(s.value);
			if( ipSolar[npSolar] <= 0 )
			{
				if( ipSolar[npSolar] < 0 )
					fprintf( ioQQQ, "Did not recognize this element symbol.\n" );
				else
					fprintf( ioQQQ, "Setting the hydrogen abundance is not possible.\n" );
				p.showLocation();
				cdEXIT(EXIT_FAILURE);
			}
			s = p.getSymbol();
			if( s.toktype == Symbol::OPERATOR && s.value == "=" )
				s = p.getSymbol();
			if( s.toktype != Symbol::NUMBER )
			{
				fprintf( ioQQQ, "Expected a value for the abundance.\n" );
				p.showLocation();
				cdEXIT(EXIT_FAILURE);
			}
			istringstream iss(s.value);
			iss >> absav[npSolar];
			++npSolar;
		}

		/* are numbers scale factors, or log of abund rel to H?? */
		lgLog = false;
		for( i=0; i < npSolar; i++ )
			if( absav[i] < 0. )
				lgLog = true;

		if( lgLog )
		{
			/* entered as log of number rel to hydrogen */
			for( i=0; i < npSolar; i++ )
				abund.ReferenceAbun[ipSolar[i]] = (realnum)exp10(absav[i]);
		}
		else
		{
			/* scale factors relative to solar */
			for( i=0; i < npSolar; i++ )
				abund.ReferenceAbun[ipSolar[i]] *= (realnum)absav[i];
		}

		/* check whether the abundances are reasonable */
		for( i=1; i < LIMELM; i++ )
		{
			if( abund.ReferenceAbun[i] > 0.2 )
			{
				fprintf( ioQQQ, " Is an abundance of %.3e relative to H reasonable for %2.2s?\n",
					abund.ReferenceAbun[i], elementnames.chElementSym[i] );
			}
		}
		return;
	}
	else if( lgIsotp )
	{
		FILE *ioDATA = open_abn_file( chFile );

		int nRead[LIMELM] = { 0 };
		string chLine;
		while( read_whole_line( chLine, ioDATA ) )
		{
			if( chLine.length() == 0 || chLine[0] == '\n' || chLine[0] == '\r' )
			{
				fprintf(ioQQQ, "PROBLEM in ABUNDANCES ISOTOPES: Encountered unexpected empty line.\n");
				cdEXIT(EXIT_FAILURE);
			}

			if( chLine[0]=='*' )
				break;

			/* skip comment */
			if( chLine[0]=='#' )
				continue;

			long i = 1;
			bool lgEOL;
			int ielem = (int)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL) - 1;
			int Aiso  = (int)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
			realnum Fiso = (realnum)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
			int j = abund.IsoAbn[ielem].setAbn( Aiso, Fiso );
			if( j == -1 )
			{
				fprintf(ioQQQ,
					"PROBLEM in ABUNDANCES ISOTOPES: Could not store isotope fraction (%7.4f) for ^%d %s\n",
					Fiso, Aiso, elementnames.chElementSym[i]);
					cdEXIT( EXIT_FAILURE );
			}
			++nRead[ielem];
		}
		fclose( ioDATA );

		/* Express all isotope fractions in terms of the most abundant isotope. */
		for( int i = 0; i < LIMELM; i++ )
		{
			if( nRead[i] != abund.IsoAbn[i].getNiso() )
			{
				fprintf(ioQQQ, "Abundaces Isotopes: %s requires %d isotope pairs to be specified, but found %d\n",
					elementnames.chElementName[i], abund.IsoAbn[i].getNiso(), nRead[i]);
				cdEXIT(EXIT_FAILURE);
			}
			if( abund.IsoAbn[i].isAnyIllegal() )
			{
				fprintf(ioQQQ, "Abundaces Isotopes: Non-positive isotope fractions are illegal.\n");
				fprintf(ioQQQ, "File: %s\t Read: %s\t iso:", chFile.c_str(), elementnames.chElementName[i]);
				abund.IsoAbn[i].prtIsoPairs( ioQQQ );
				cdEXIT(EXIT_FAILURE);
			}
			abund.IsoAbn[i].normAbn( );
			//	abund.IsoAbn[i].prtIsoPairs( stdout );
		}

		return;
	}
	//return;
	//Should never reach this return, because if & else branch return
}
