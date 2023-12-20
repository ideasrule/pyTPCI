/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*cdGetLineList routine to read in master list of emission line wavelengths and ids, for
 * generating loc grids */
#include "cddefines.h"
#include "cddrive.h"
#include "lines.h"
#include "parser.h"

/* return value is number of lines, -1 if file could not be opened */
long int cdGetLineList( 
	/* chFile is optional filename, if void then use BLRLineList,
	 * if not void then use file specified */
	const string& chFile,
	/* array of line identifications */
	vector<LineID>& lineids)
{
	DEBUG_ENTRY( "cdGetLineList()" );

	/* first check that cdInit has been called, since we may have to write
	 * error output */
	if( !lgcdInitCalled )
	{
		fprintf(stderr," cdInit must be called before cdGetLineList.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* use default filename LineList_BLR.dat if void string, else use file specified */
	string chFilename = ( chFile.length() == 0 ) ? "LineList_BLR.dat" : chFile;

	DataParser d( chFilename, ES_STARS_AND_BLANKS, AS_TRY );
	if( !d.isOpen() )
	{
		/* did not find file, return -1 */
		return -1;
	}

	// make sure we have a blank slate
	lineids.clear();

	/* actually read and save the lines */
	while( d.getline() )
	{
		/* this checks for in-file end-of-data markers */
		if( d.lgEODMarker() )
			break;

		LineID line;
		d.getLineID(line);

		/* make sure there is no trailing junk */
		d.checkEOL();

		lineids.emplace_back(line);
	}

	/* return number of lines we found */
	return lineids.size();
}
