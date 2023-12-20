/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PARSE_SPECIES_H_
#define PARSE_SPECIES_H_

class Parser;
void ParseSpecies(Parser &p);

class species;
void setProperties(species& sp);

bool speciesOff(const string& label);

void speciesCheck();

#endif // PARSE_SPECIES_H_

