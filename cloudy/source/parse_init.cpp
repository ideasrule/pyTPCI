/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseInit bring an initialization file into input stream before parse  */
#include "cddefines.h"
#include "input.h"
#include "trace.h"
#include "cddrive.h"
#include "parser.h"

void ParseInit(Parser &p)
{
	DEBUG_ENTRY( "ParseInit()" );

	/* bring an initialization file into input stream before parsing */

	/* check whether single quote on line, this was used in c90 */
	if( p.nMatch( "\'" ) )
	{
		fprintf( ioQQQ, 
			" ParseInit found a single quote on this line.  This was used"
			 " for file names in C90, but double quotes are used now.\n" );
		fprintf( ioQQQ, " The single quote has been ignored.\n" );
	}

	string chName;
	if( p.nMatch( "\"" ) )
	{
		/* 
		 * if a quote occurs on the line then get the ini file name 
		 * this will also set the name in chCard and OrgCard to spaces
		 * so later keywords do not key off it
		 */
		if( p.GetQuote( chName ) )
			p.StringError();
	}
	else
	{
		/* no quote appeared, so this is the default name, cloudy.ini */
		chName = "cloudy.ini";
	}

	/* at this point we have init file name, now make full name 
	 * this can be a local file, or on the path if the key path appears */
	ParseInitFile(chName);
}

void ParseInitFile(const string& chName)
{
	DEBUG_ENTRY( "ParseInitFile()" );

	/* option to get cloudy.ini from a path */
	FILE * ioInitFile = open_data( chName, "r" );

	++input.curInclLevel;

	/* at this point the init file is open, now bring it into the command stack */
	string chLocal;
	while( read_whole_line(chLocal,ioInitFile) )
	{
		if( lgInputEOF(chLocal) )
			break;

		/* stuff the command line into the internal stack */
		(void)cdRead(chLocal);
	}

	--input.curInclLevel;

	fclose(ioInitFile);
}
