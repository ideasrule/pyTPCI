/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef GENERIC_STATE_H_
#define GENERIC_STATE_H_

#include "quantumstate.h"

class molezone;
class genericState
{
public:
	genericState() : sp(NULL), val(NULL) {}
	explicit genericState(const molezone* sp) : sp(sp), val(NULL) {}
	explicit genericState(const qStateConstProxy& q) 
		: sp(NULL), q(q), val(NULL) {}
	explicit genericState(const char *valLabel, const double *val) 
		: sp(NULL), val(val), valLabel(valLabel) {}
	const molezone* sp;
	qStateConstProxy q;
	const double* val;
	string valLabel;
	bool associated() const;
	string label() const;
	string database() const;
	bool isSpecies() const;
};

double column(const genericState&);
double density(const genericState&);
double depart(const genericState&);
double energy(const genericState&);
double levels(const genericState&);

/** extractLevels -- parse string to extract species, low and high levels from, e.g., "Cu+13[1:100]"
 * \param[in]  chLabel   input species & levels string
 * \param[out] species   output species string
 * \param[out] nLevelLo  low level index
 * \param[out] nLevelHi  high level index
 * \param[out] lgLevels  boolean, true if levels given
 */
void extractLevels( const string &chLabel,
                        string &species,
                        long &nLevelLo,
                        long &nLevelHi,
                        bool &lgLevels );

/** getLevelsGeneric -- get the requested species levels
 * \param [in]  chLabel	    species label
 * \param [in]  lgValidate  flag to report errors in the level request
 * \param [out] LevelList   list of species levels
 */
const molezone *getLevelsGeneric(const string &chLabel, bool lgValidate, vector<long> &LevelList);

/** matchGeneric -- get a list of all species that match request
 * \param[in] chLabel	 species label
 * \param[in] lgValidate flag to report errors in the level request
 * \return               list of genericState species
 */
vector<genericState> matchGeneric(const string &chLabel, bool lgValidate);

/** getSpeciesGeneric - acquire the species matching the input string
 *
 * \param [in] speciesLabel  input species string
 *
 * \returns species     reference to requested species
 */
genericState getSpeciesGeneric( const string &speciesLabel );


#endif // GENERIC_STATE_H_
