/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "date.h"
#include "version.h"
#include "service.h"

static const int CLD_MAJOR = 23;
static const int CLD_MINOR = 1;
static const int CLD_PATCH = 0;

#ifdef REVISION
static const char* revision = REVISION;
#else
static const char* revision = "rev_not_set";
#endif

static const string release_branch = "release";

t_version::t_version()
{
	static const char chMonth[12][4] =
		{ "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

	// REVISION is a dash-delimited string; it may consist of the branch name,
	// followed by a tag, or by a git SHA1 commit string, and a modification
	// identifier

	if( strcmp( revision, "" ) != 0 && strcmp( revision, "rev_not_set" ) != 0 )
	{
		vector<string> Part;
		string rev = revision;
		Split( rev, "-", Part, SPM_RELAX );

		if( Part[ 0 ] == release_branch )
		{
			lgRelease = true;
			Part.erase( Part.begin() );
		}

		string rev_pr = "";
		for( size_t i = 0 ; i < Part.size(); i++ )
		{
			if( i < Part.size()-1 )
				rev_pr += Part[i] + ", ";
			else
				rev_pr += Part[i];
		}
		chVersion = "(" + string( rev_pr ) + ")";
	}
	else
	{
		// create a default version string in case the code is not versioned with git

		lgRelease = true;

		ostringstream oss;
		if( lgRelease )
		{
			oss << setfill('0') << setw(2) << CLD_MAJOR << "." << setw(2) << CLD_MINOR;
			if( CLD_PATCH > 0 )
				oss << " (patch level " << CLD_PATCH << ")";
		}
		else
		{
			oss << setfill('0') << setw(2) << YEAR%100 << "." << setw(2) << MONTH+1;
			oss << "." << setw(2) << DAY;
		}
		chVersion = oss.str();
	}

	ostringstream oss2;
	oss2 << setfill('0') << setw(2) << YEAR%100 << chMonth[MONTH] << setw(2) << DAY;
	chDate = oss2.str();

	char mode[8];
	if( sizeof(int) == 4 && sizeof(long) == 4 && sizeof(long*) == 4 )
		strncpy( mode, "ILP32", sizeof(mode) );
	else if( sizeof(int) == 4 && sizeof(long) == 4 && sizeof(long*) == 8 )
		strncpy( mode, "IL32P64", sizeof(mode) );
	else if( sizeof(int) == 4 && sizeof(long) == 8 && sizeof(long*) == 8 )
		strncpy( mode, "I32LP64", sizeof(mode) );
	else if( sizeof(int) == 8 && sizeof(long) == 8 && sizeof(long*) == 8 )
		strncpy( mode, "ILP64", sizeof(mode) );
	else
		strncpy( mode, "UNKN", sizeof(mode) );

	bool flag[2];
	flag[0] = ( cpu.i().min_float()/2.f > 0.f );
	flag[1] = ( cpu.i().min_double()/2. > 0. );

	/* now generate info on how we were compiled, including compiler version */
	ostringstream oss3;
	oss3 << "Cloudy compiled on " << __DATE__ << " in OS " << __OS << " using the ";
	oss3 << __COMP << " " << __COMP_VER << " compiler. Mode " << mode << ", ";
	oss3 << "denormalized float: " << TorF(flag[0]) << " double: " << TorF(flag[1]) << ".";
	chInfo = oss3.str();
}
