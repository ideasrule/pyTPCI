/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*save_average routine to read in list of species to output as averages */
#include "cddefines.h"
#include "cddrive.h"
#include "elementnames.h"
#include "save.h"
#include "parser.h"

/* return value is number of lines, -1 if file could not be opened */
void parse_save_average(Parser &p,
			/* the file we will write to */
			long int ipPun,
			ostringstream& chHeader)
{
	DEBUG_ENTRY( "parse_save_average()" );

	/* ensure the memory from any previous call is freed */
	save.SaveAverageFree(ipPun);
	
	/* use this to count number of species, and will assert equal to
	 * total allocated above */
	size_t nLine = 0;
	/* keep reading until we hit END */
	while( 1 )
	{
		p.getline();
		if( p.m_lgEOF )
		{
			fprintf( ioQQQ, 
						" Save average hit EOF while reading list; use END to end list.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		if (p.hasCommand("END" ))
			break;
		/* count number of species we will save */
		++nLine;
		if( p.nMatch("TEMP" ))
		{
			/* temperature */
			save.chAverageType[ipPun].push_back("TEMP");
		}
		else if( p.nMatch("COLU" ))
		{
			/* column density */
			save.chAverageType[ipPun].push_back("COLU");
		}
		else if( p.nMatch("IONI" ))
		{
			/* ionization fraction */
			save.chAverageType[ipPun].push_back("IONI");
		}
		else
		{
			fprintf(ioQQQ,"PROBLEM one of the jobs TEMPerature, COLUmn density, or IONIzation, must appear.\n");
			cdEXIT(EXIT_FAILURE);
		}
		
		/* get element name, a string we shall pass on to the routine
		 * that computes final quantities */
		int i = p.GetElem( );
		if( i < 0 )
		{
			/* no name found */
			fprintf(ioQQQ, "save average did not see an element on this line, sorry\n");
			p.PrintLine(ioQQQ);
			cdEXIT(EXIT_FAILURE);
		}
		save.chSaveSpecies[ipPun].push_back(elementnames.chElementNameShort[i]);
		
		/* now get ionization stage */
		save.nAverageIonList[ipPun].push_back((int) p.FFmtRead());
		if( p.lgEOL() )
		{
			/* error - needed that ionization stage */
			p.NoNumb("ionization stage" );
		}
		
		/* look for volume keyword, otherwise will be radius 
		 * only used for some options */
		if( p.nMatch( "VOLU" ) )
		{
			/* volume */
			save.nAverage2ndPar[ipPun].push_back(1);
		}
		else
		{
			/* radius */
			save.nAverage2ndPar[ipPun].push_back(0);
		}
	}
	
   ASSERT(nLine == save.chAverageType[ipPun].size());
   ASSERT(nLine == save.chSaveSpecies[ipPun].size());
   ASSERT(nLine == save.nAverageIonList[ipPun].size());
   ASSERT(nLine == save.nAverage2ndPar[ipPun].size());
	save.nAverageList[ipPun] = nLine;
/*#		define PADEBUG*/
#		ifdef PADEBUG
	fprintf(ioQQQ , "DEBUG save_average %li species read in.\n", 
			  save.nAverageList[ipPun] );
#		endif

#		ifdef PADEBUG
	for( i=0; i<nLine ; ++i )
	{
		fprintf(ioQQQ, "PDDEBUG %s %s %i %i\n",
				  save.chAverageType[ipPun][i].c_str(),
				  save.chSaveSpecies[ipPun][i].c_str() ,
				  save.nAverageIonList[ipPun][i] ,
				  save.nAverage2ndPar[ipPun][i] );
	}
#		endif
	
	/* save headers */
	sncatf(chHeader, "#averages");
	for( size_t i=0; i<nLine ; ++i )
	{
		sncatf(chHeader, "\t %s %s %i %i",
				  save.chAverageType[ipPun][i].c_str(),
				  save.chSaveSpecies[ipPun][i].c_str() ,
				  save.nAverageIonList[ipPun][i] ,
				  save.nAverage2ndPar[ipPun][i] );
	}
	sncatf(chHeader, "\n");
}

void save_average( 
	/* the file we will write to */
	long int ipPun)
{
	DEBUG_ENTRY( "save_average()" );

	/* do the output */
	for( long i=0; i<save.nAverageList[ipPun] ; ++i )
	{
		double result;
		char chWeight[7];
		if( save.nAverage2ndPar[ipPun][i] == 0 )
			strcpy( chWeight , "RADIUS");
		else
			strcpy( chWeight , "VOLUME");
		
		if( save.chAverageType[ipPun][i] == "TEMP" )
		{
			/* temperature */
			if( cdTemp( 
					 save.chSaveSpecies[ipPun][i].c_str() ,
					 save.nAverageIonList[ipPun][i] ,
					 &result , 
					 chWeight ) )
			{
				fprintf( ioQQQ, " save average temperature could not identify the species.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else if( save.chAverageType[ipPun][i] == "IONI" )
		{
			/* ionization fraction 
			 * H2 is a special case, HYDRO 0 requests
			 * the H2 fraction, n(H2)/n(H) */
			if( "HYDR" == save.chSaveSpecies[ipPun][i] &&
				 save.nAverageIonList[ipPun][i]== 0 )
				save.chSaveSpecies[ipPun][i] = "H2  ";
			if( cdIonFrac( 
					 save.chSaveSpecies[ipPun][i].c_str() ,
					 save.nAverageIonList[ipPun][i] ,
					 &result , 
					 chWeight ,
					 false 
					 ) )
			{
				fprintf( ioQQQ, " save average ionization fraction could not identify the species.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else if( save.chAverageType[ipPun][i] == "COLU" )
		{
			/* column density */
			if( cdColm( 
					 save.chSaveSpecies[ipPun][i].c_str() ,
					 save.nAverageIonList[ipPun][i] ,
					 &result ) )
			{
				fprintf( ioQQQ, " save average column density fraction could not identify the species.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
			TotalInsanity();
		
		// <=C13 gave log quantities, now do linear but accept log key for backwards compatibility
		fprintf(save.params[ipPun].ipPnunit, "\t %e", PrtLogLin( result ) );
	}
	fprintf(save.params[ipPun].ipPnunit, "\n");
}
