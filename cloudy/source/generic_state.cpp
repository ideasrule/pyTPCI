/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"

#include "generic_state.h"

#include "dynamics.h"
#include "mole.h"
#include "phycon.h"
#include "radius.h"
#include "rfield.h"
#include "species.h"

double column(const genericState&gs)
{
	if ( gs.sp != NULL )
		return gs.sp->column;
	else if (gs.q.associated())
		return gs.q.ColDen();
	else if (gs.val != NULL)
		return *gs.val;
	else
		return 0.;
}

double density(const genericState& gs)
{
	if ( gs.sp != NULL )
		return gs.sp->den;
	else if (gs.q.associated())
		return gs.q.Pop();
	else if (gs.val != NULL)
		return *gs.val;
	else
		return 0.;
}

double depart(const genericState& gs)
{
	if ( gs.sp != NULL )
		return 1.0;
	else if (gs.q.associated())
		return gs.q.DepartCoef();
	else if (gs.val != NULL)
		return *gs.val;
	else
		return 1.0;
}

double energy(const genericState& gs) 
{
	if ( gs.sp != NULL )
		return 0.0;
	else if (gs.q.associated())
		return AnuUnit( gs.q.energy().Ryd() );
	else if (gs.val != NULL)
		return *gs.val;
	else
		return 0.0;
}

double levels(const genericState& gs) 
{
	if ( gs.sp != NULL && gs.sp->levels != NULL )
		return (double) gs.sp->levels->size();
	else if (gs.val != NULL)
		return *gs.val;
	else
		return 0.0;
}

bool genericState::isSpecies() const
{
	if( sp != NULL && sp != null_molezone )
		return true;
	else
		return false;
}

string genericState::label() const
{
	if( isSpecies() )
		return sp->global()->label;
	else if (q.associated())
		return q.chLabel();
	else if (val != NULL)
		return valLabel;
	else
		return "";
}

string genericState::database() const
{
	if( isSpecies() && sp->dbase != NULL )
		return sp->dbase->database;
	else
		return "";
}

bool genericState::associated() const
{
	if( isSpecies() )
		return true;
	else if (q.associated())
		return true;
	else if (val != NULL)
		return true;
	else
		return false;
}

static const long IGNORE_LEVEL = -1;
void extractLevels( const string &chLabel,
			string &species,
			long &nLevelLo,
			long &nLevelHi,
			bool &lgLevels )
{
	DEBUG_ENTRY( "extractLevels()" );

	nLevelLo = -1;
	nLevelHi = -1;
	lgLevels = false;

	size_t lbrac = chLabel.find( "[" );
	size_t rbrac = chLabel.find( "]" );

	if( ( lbrac < string::npos && rbrac == string::npos ) ||
		( lbrac == string::npos && rbrac < string::npos ) )
	{
		fprintf( ioQQQ,
			"PROBLEM: In request for species \"%s\"\n"
			"Excitation levels must be in balanced brackets '[]'\n",
			chLabel.c_str() );
		cdEXIT( EXIT_FAILURE );
	}

	species = chLabel.substr( 0, lbrac );

	// no levels given
	if( lbrac == string::npos )
		return;



	lgLevels = true;

	string levels = chLabel.substr( lbrac+1, rbrac-lbrac-1 );
	size_t colon = levels.find( ":" );

	if( colon == string::npos )
	{
		// case where only one level given, e.g., "C+2[4]"
		nLevelHi = nLevelLo = long( atoi( levels.c_str() ) );
		if( nLevelLo == 0 )
		{
			fprintf( ioQQQ,
				"PROBLEM: In request for species \"%s\"\n"
				"Excitation levels must be at >= 1\n",
				chLabel.c_str() );
			cdEXIT(EXIT_FAILURE);
		}
	}
	else
	{
		// case where at least one level given,
		// e.g., "C+2[1:4]", "C+2[:4]", "C+2[1:]"

		if( levels.substr( 0, colon ).length() == 0 )
			nLevelLo = 1;
		else
		{
			nLevelLo = long( atoi( levels.substr( 0, colon ).c_str() ) );
			if( nLevelLo <= 0 )
			{
				fprintf( ioQQQ,
					"PROBLEM: In request for species \"%s\"\n"
					"Excitation levels must be at >= 1\n",
					chLabel.c_str() );
				cdEXIT(EXIT_FAILURE);
			}
		}
		if( levels.substr( colon+1 ).length() == 0 )
			nLevelHi = 0; // tag for all levels
		else
		{
			nLevelHi = long( atoi( levels.substr( colon+1 ).c_str() ) );
			if( nLevelHi <= 0 )
			{
				fprintf( ioQQQ,
					"PROBLEM: In request for species \"%s\"\n"
					"Excitation levels must be at >= 1\n",
					chLabel.c_str() );
				cdEXIT(EXIT_FAILURE);
			}
			if( nLevelHi < nLevelLo )
			{
				fprintf( ioQQQ,
					"PROBLEM: In request for species \"%s\"\n"
					"Higher level must be >= Lower level\n",
					chLabel.c_str() );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	{
		bool _DEBUG = false;
		if( _DEBUG )
		{
			fprintf( stdout,
				"chLabel = '%s'\t nLevelLo = %ld\t nLevelHi = %ld\t lgLevels = %s\n",
				chLabel.c_str(), nLevelLo, nLevelHi,
				lgLevels ? "true" : "false" );
			cdEXIT( EXIT_FAILURE );
		}
	}
}

const molezone *getLevelsGeneric( const string &chLabel, bool lgValidate, vector<long> &LevelList )
{
	DEBUG_ENTRY( "getLevelsGeneric()" );

	string chLabel_species = "";

	long nLevelLo = -1,
		nLevelHi = -1;
	bool lgLevels = false;

	extractLevels( chLabel,
			chLabel_species,
			nLevelLo,
			nLevelHi,
			lgLevels );


	const molezone *sp = findspecieslocal(chLabel_species.c_str());

	if ( sp == null_molezone)
	{
		if (lgValidate)
		{
			fprintf(ioQQQ,"PROBLEM: Request for unidentified species \"%s\"\n",chLabel_species.c_str());
			cdEXIT(EXIT_FAILURE);
		}
	}
	else if ( lgLevels )
	{
		if ( ! sp->levels )
		{
			if (lgValidate)
			{
				fprintf(ioQQQ,"PROBLEM: Request for level in unresolved species \"%s\"\n",chLabel_species.c_str());
				cdEXIT(EXIT_FAILURE);
			}
			// If no levels, map all population into ground state
			if (nLevelHi == 0)
				nLevelHi = nLevelLo;
			for (int nLevel = nLevelLo; nLevel<=nLevelHi; ++nLevel)
			{
				if ( nLevel == 1 )
					LevelList.push_back( nLevel );
				else
					LevelList.push_back( IGNORE_LEVEL );
			}
		}
		else if ( size_t(nLevelLo) > sp->levels->size() || size_t(nLevelHi) > sp->levels->size() )
		{
			if (lgValidate)
			{
				fprintf(ioQQQ,"PROBLEM: Request for level \"%s\", but species only has %lu levels\n",
					chLabel.c_str() ,(unsigned long)sp->levels->size());
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			// Offset by one to convert from numeric index to C array index
			if (nLevelHi == 0)
				nLevelHi = sp->levels->size();
			for (int nLevel = nLevelLo; nLevel<=nLevelHi; ++nLevel)
			{
				LevelList.push_back( nLevel-1 );
			}
		}
	}

	return sp;
}

vector<genericState> matchGeneric(const string &chLabel, bool lgValidate)
{
	DEBUG_ENTRY("matchGeneric()");
	vector<genericState> v;

	if (chLabel[0] == '*')
	{
		if (chLabel.compare( 1, 5, "depth" ) == 0)
		{
			v.push_back(genericState("*depth",&radius.depth_mid_zone));
			return v;
		}
		else if (chLabel.compare( 1, 2, "AV" ) == 0)
		{
			v.push_back(genericState("*AV",&rfield.extin_mag_V_point));
			return v;
		}
		else if (chLabel.compare( 1, 3, "AVx" ) == 0)
		{
			v.push_back(genericState("*AVx",&rfield.extin_mag_V_extended));
			return v;
		}
		else if (chLabel.compare( 1, 4, "time" ) == 0)
		{
			v.push_back(genericState("*time",&dynamics.time_elapsed));
			return v;
		}
		else if (chLabel.compare( 1, 4, "temp" ) == 0)
		{
			v.push_back(genericState("*temp",&phycon.te));
			return v;
		}
		else if (chLabel.compare( 1, 3, "all" ) == 0)
		{
			for( size_t i=0; i<mole_global.list.size(); ++i )
			{
				v.push_back( genericState(&(mole.species[i])) );
			}
			return v;
		}
	}

	vector<long> LevelList;
	const molezone *sp = getLevelsGeneric( chLabel, lgValidate, LevelList );

	if( sp != null_molezone && LevelList.size() > 0 )
	{
		for( vector<long>::iterator nLevel = LevelList.begin();
			nLevel != LevelList.end(); ++nLevel)
		{
			if( sp->levels )
				v.push_back(genericState((*sp->levels)[ *nLevel ]));
			else if( *nLevel != IGNORE_LEVEL )
				v.push_back(genericState(sp));
			else
				v.push_back(genericState());
		}
	}
	else if( sp != null_molezone )
	{
		v.push_back(genericState(sp));
	}
	return v;
}

genericState getSpeciesGeneric(const string &speciesLabel)
{
	DEBUG_ENTRY( "getSpeciesGeneric()" );

	vector<genericState> v = matchGeneric( speciesLabel, true );
	if( v.size() == 0 || ! v[0].isSpecies() )
	{
		fprintf( ioQQQ, "PROBLEM: Invalid species '%s' requested.\n",
				speciesLabel.c_str() );
		cdEXIT( EXIT_FAILURE );
	}

	return v[0];
}
