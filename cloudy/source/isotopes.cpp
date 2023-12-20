/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*LoadIsotopes	read in the nuclear isotope data and allocate space */
#include "cddefines.h"
#include "abund.h"
#include "atmdat.h"

/*LoadIsotopes	read in the nuclear isotope data and allocate space */
void LoadIsotopes ( )
{
	DEBUG_ENTRY( "LoadIsotopes()" );

	string chFile = "isotopes.dat";
	FILE* ioDATA = open_data( chFile, "r" );	// will abort if not found

	string chLine;			/*to store file lines*/
	while( read_whole_line( chLine, ioDATA ) )
	{
		if( chLine.length() == 0 || chLine[0]=='\n' || chLine[0]=='\r' )
		{
			fprintf(ioQQQ,
				"PROBLEM in LoadIsotopes: Encountered unexpected empty line in %s.\n",
				chFile.c_str());
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
		ASSERT( ielem >= 0 );

		int Aiso  = (int) FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
		ASSERT( Aiso >  0  );

		FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

		realnum mass = (realnum) FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
		ASSERT( mass >  0. );
	
	
		realnum spin = 0., magm = 0.;
		double tmp = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
		if( tmp )
		{
			spin = (realnum) tmp;
			ASSERT( spin >= 0. );

			magm = (realnum) FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
			ASSERT( magm != 0. );
		}
		//	printf("ielem = %d\t Aiso = %d\t mass = %9.6f\t spin = %3.1f\t magm = %10.6f\n",
		//		ielem, Aiso, mass, spin, magm);

		abund.IsoAbn[ielem].setData( Aiso, mass, spin, magm);
	}

	fclose( ioDATA );
	return;
}
