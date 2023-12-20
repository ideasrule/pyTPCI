/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* input_readarray read input commands from array where images are stored *
 * returns chCard, which will have <=80 characters before eol    *
 * line image is up and low case                                 */
/*lgInputComment - parse comment - check if argument is comment string */
#include "cddefines.h"
#include "trace.h"
#include "optimize.h"
#include "parser.h"
#include "input.h"

t_input input;

/* lgIsCommentSeq - is the string s[p] pointing to a comment sequence?
 * if lgReportVisible is true, visible comments will be reported, otherwise not */
bool lgIsCommentSeq( const string& s, size_t p, bool lgReportVisible )
{
	DEBUG_ENTRY( "lgIsCommentSeq()" );

	if( s.substr(p,2) == "##" )
		return true;
	else if( lgReportVisible && s[p] == '#' )
		return true;
	else
		return false;
}

/* lgIsExpungedCommentSeq - does the string s start with an old-style comment? */
bool lgIsExpungedCommentSeq( const string& s )
{
	DEBUG_ENTRY( "lgIsExpungedCommentSeq()" );

	string s2 = s.substr(0,2);
	if( s2 == "C " || s2 == "//" || s[0] == '%' || s[0] == '*' )
		return true;
	else
		return false;
}

/* lgInputComment - parse comment - check if argument is comment string, 
 * either upper or lower case -
 * returns true if line is a comment, false if not */
bool lgInputComment( const string& chLine )
{
	DEBUG_ENTRY( "lgInputComment()" );

	return lgIsCommentSeq( chLine, 0, true );
}

/* lgInputEOF - is this line an EOF marker? */
bool lgInputEOF( const string& chLine )
{
	DEBUG_ENTRY( "lgInputEOF()" );

	if( chLine.length() == 0 || chLine[0] == '\n' || chLine[0] == '\r' || chLine[0] == ' ' )
		return true;
	else if( chLine.substr(0,3) == "***" )
		return true;
	else
		return false;
}

/* StripComment- strips comment part off the command line s
 * if lgStripVisible is false, visible comments are retained
 * hidden comments are always stripped
 * this routine also removes underscores and brackets */
void StripComment( string& s, bool lgStripVisible )
{
	DEBUG_ENTRY( "StripComment()" );

	for( size_t p=0; p < s.length(); )
	{
		if( s[p] == '\"' )
		{
			string buf;
			p = GetString( s, p, buf );
		}
		else if( lgIsCommentSeq(s,p,lgStripVisible) )
		{
			s.erase(p);
			break;
		}
		else if( lgIsCommentSeq(s,p,true) )
		{
			break;
		}
		else if( s[p] == '_' )
		{
			s[p++] = ' ';
			input.lgUnderscoreFound = true;
		}
		else if( s[p] == '[' || s[p] == ']' )
		{
			s[p++] = ' ';
			input.lgBracketFound = true;
		}
		else
		{
			++p;
		}
	}
}

/* GetString: retrieve a string between double quotes
 * s : string to be parsed
 * s[p] : place in string where to start parsing, should point to first set of double quotes
 * buf : buffer that will hold the string between double quotes
 * return value : pointer just beyond second set of double quotes, or string::npos in case of failure
 *                (second set of double quotes wasn't found) */
size_t GetString( const string& s, size_t p, string& buf )
{
	DEBUG_ENTRY( "GetString()" );

	ASSERT( s[p] == '\"' );

	buf.clear();
	for( ++p; p < s.length(); )
	{
		if( s[p] == '\\' )
			buf.push_back( GetEscape(s,p) );
		else if( s[p] == '\"' )
			return ++p;
		else
			buf.push_back( s[p++] );
	}
	// no second set of double quotes was found
	buf.clear();
	return string::npos;
}

char GetEscape( const string& s, size_t& p )
{
	DEBUG_ENTRY( "GetEscape()" );

	// This routine is the placeholder for treating character escape
	// sequences. For the moment we treat none. The routine assumes that
	// an esacpe sequence will always generate a single character. This is
	// true for all escape sequences except unicode sequences...
	// on exit, p will point one character beyond the escape sequence.
	return s[p++];
}

/* MakeInputLine: generate input lines for optimizer and grid commands
 * i : index of varied command, should be less than optimize.nvary */
string MakeInputLine(long i)
{
	DEBUG_ENTRY( "MakeInputLine()" );

	ASSERT( i < optimize.nvary );

	int len = 200;
	string cmd;
	while( true )
	{
		int res = 0;
		char* buf = new char[len];
		if( optimize.nvarxt[i] == 1 )
		{
			/* case with 1 parameter */
			res = snprintf( buf, len, optimize.chVarFmt[i], optimize.vparm[0][i] );
		}
		else if( optimize.nvarxt[i] == 2 )
		{
			/* case with 2 parameter */
			res = snprintf( buf, len, optimize.chVarFmt[i], optimize.vparm[0][i], optimize.vparm[1][i]);
		}
		else if( optimize.nvarxt[i] == 3 )
		{
			/* case with 3 parameter */
			res = snprintf( buf, len, optimize.chVarFmt[i], optimize.vparm[0][i], optimize.vparm[1][i],
					optimize.vparm[2][i] );
		}
		else if( optimize.nvarxt[i] == 4 )
		{
			/* case with 4 parameter */
			res = snprintf( buf, len, optimize.chVarFmt[i], optimize.vparm[0][i], optimize.vparm[1][i],
					optimize.vparm[2][i], optimize.vparm[3][i] );
		}
		else if( optimize.nvarxt[i] == 5 )
		{
			/* case with 5 parameter */
			res = snprintf( buf, len, optimize.chVarFmt[i], optimize.vparm[0][i], optimize.vparm[1][i],
					optimize.vparm[2][i], optimize.vparm[3][i], optimize.vparm[4][i]);
		}
		else
		{
			fprintf(ioQQQ,"The number of variable options on this line makes no sense to me\n");
			cdEXIT(EXIT_FAILURE);
		}
		if( res >= len )
		{
			// output in buf was truncated -> try again with bigger buffer
			len = res + 10;
			delete[] buf;
		}
		else
		{
			// output in buf is complete -> return
			cmd = buf;
			delete[] buf;
			break;
		}
	}

	return cmd;
}

void t_input::echo( FILE *ipOUT ) const
{
	for( size_t i=0; i < crd.size(); ++i )
		if( crd[i]->InclLevel == 0 && crd[i]->lgVisible )
			fprintf( ipOUT, "%s\n", crd[i]->chCardSav.c_str() );
}

/** readarray read input commands from array where images are stored *
 * and return this command */
string t_input::readarray(bool *lgEOF)
{
	DEBUG_ENTRY( "t_input::readarray()" );

	/* usual case, reading commands from start of array
	 * nRead points to one plus the array element with the next line, it is
	 * one on the first call, which references line[0] */
	++nRead;

	if( nRead >= long(crd.size()) )
	{
		*lgEOF = true;
		if( trace.lgTrace )
			fprintf( ioQQQ, "t_input::readarray returns EOF\n" );
		return string();
	}
	else
	{
		*lgEOF = false;
		if( trace.lgTrace )
			fprintf( ioQQQ, "t_input::readarray returns=%s=\n", crd[nRead]->chCardSav.c_str() );
		return crd[nRead]->chCardSav;
	}
}

/** peekarray look ahead at the next command, but do not retrieve it yet */
string t_input::peekarray(bool *lgEOF)
{
	DEBUG_ENTRY( "t_input::peekarray()" );

	if( nRead+1 >= long(crd.size()) )
	{
		*lgEOF = true;
		return string();
	}
	else
	{
		*lgEOF = false;
		return crd[nRead+1]->chCardSav;
	}
}

/** input_readvector: read numbers from the file chFile and store them in a vector */
void input_readvector(const string& chFile, /**< file name to read from */
					  vector<double>& vec)  /**< vector - the numbers that were read from the input line(s) */
{
	DEBUG_ENTRY( "input_readvector()" );

	vec.clear();
	DataParser d(chFile, ES_NONE);
	while( d.getline() )
	{
		double tmpValue;
		while( d.getTokenOptional(tmpValue) )
			vec.emplace_back(tmpValue);
	}
}
