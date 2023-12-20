/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PARSE_H_
#define PARSE_H_

/**ParseCommands main command line parser, called by Cloudy to decode commands, 
 * it then call other routines to parse specific commands */
void ParseCommands(void);

/**initialize save file pointers */
void SaveFilesInit(void);

/**close all open save files 
\param lgFinal - close ALL files, regardless of "no clobber" status when true
*/
void CloseSaveFiles( bool lgFinal );

#endif /* PARSE_H_ */
