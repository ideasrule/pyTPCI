/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef TAULINES_H_
#define TAULINES_H_

#include "transition.h"
#include "species.h"
#include "container_classes.h"

class CollRateCoeffArray;
class CollSplinesArray;
class StoutCollArray;
class species;

extern bool lgStatesAdded;
extern bool lgLinesAdded;
extern qList AnonStates;

extern char **chSpecies;
extern vector<species> dBaseSpecies;
extern vector<qList> dBaseStates;
extern vector< multi_arr<int,2> > ipdBaseTrans;
extern vector<TransitionList> dBaseTrans;
extern multi_arr<CollRateCoeffArray,2> AtmolCollRateCoeff;
extern vector< multi_arr<CollSplinesArray,3> > AtmolCollSplines;
extern vector< StoutCollArray > StoutCollData;
extern long int nSpecies;
/*************************/
void database_readin( void );
void dBaseTrim();
void dBaseUpdateCollCoeffs(void);
void dBase_solve(void );

/**
 * this is a dummy optical depth array for non-existant lines
 */
extern TransitionProxy::iterator TauDummy;

/** this is the set of extra lines,
 * ExtraLymanLines[ipISO][ipZ][n]*/
extern multi_arr<int,3> ipExtraLymanLines;
extern vector<vector<TransitionList> > ExtraLymanLines;

/** the set of inner shell lines */
extern TransitionList UTALines;

/** this is the number of level 1 lines, and is set in atmdat_readin
 * by counter number of data lines in level1.dat */
extern long int nLevel1;
/**extern transition TauLines[NTAULINES+1];*/

/** these are the public parts of the hyperfine structure line transfer info
 * data gathered from hyperfine.dat using routines in hyperfine.c
 * the structure containing the hfs line information */
/* abundances of these isotopes relative to main species are in hyperfine.h */
extern TransitionList HFLines;

/**
 * main line arrays for hydrogenic ions<br><br>
 *
 * first dimension is atomic number<br>
 * second dim is upper level<br>
 * third dim is lower level<br>
 * nta dim is set of pointers for quantities within line transfer arrays<br>
 * in the forc translation, the upper level was too low by 1, since the<br>
 * fortran was starting at 1.  the lower dim was not changed by translation<br>
 * since it started from ip1s = 0 <br>
 * any place where the third dim has -1 is probably a remnant from forc and is wrong<br>
 */

//extern vector<vector<multi_arr<int,2> > > ipTransitions;
extern vector<vector<TransitionList> > Transitions;

extern vector<TransitionList> AllTransitions;
extern void checkTransitionListOfLists(vector<TransitionList>&); 

extern multi_arr<int,2> ipFe2LevN;
extern TransitionList Fe2LevN;

/** lines forming from doubly excited states */
extern multi_arr<int,3> ipSatelliteLines; /* [ipISO][nelem][level] */
extern vector<vector<TransitionList> > SatelliteLines; /* [ipISO][nelem][level] */

/** this will be set true once space is allocated for the HydroLines array.
 * from then on any HYDROGENIC LEVELS command will be ignored, this is
 * set to false in cddefines.c */
extern bool lgHydroAlloc;

/* all of Dima's level 2 lines */

/**number of level 2 lines, dim for WindLine array */
const int NWINDDIM = 6744;

/** this is set to 0 with no atom_level2 command, normally
 * equal to NWINDDIM, definition is in cddefines.c */
extern long	nWindLine;

/* these are the level two lines themselves */
/** pointers to element and ion, TauLine2[line number][pointer within vector] */
extern TransitionList TauLine2;

extern vector<realnum> cs1_flag_lev2;

/** findTrans_byQuantNumb	Identify transition of a given species by its quantum numbers.
 *
 * \param speciesLabel	(in) label of species, e.g., "H  1"
 * \param n_hi		(in) principal quantum number of upper level in transition
 * \param l_hi		(in) orbital angular momentum of upper level
 * \param S_hi		(in) multiplicity of upper level
 * \param n_lo		(in) principal quantum number of lower level
 * \param l_lo		(in) orbital angular momentum of lower level
 * \param S_lo		(in) multiplicity of lower level
 *
 * \return	pointer to transition in species list of transitions. If no match is
 * 		found, return value is past the end of the entire list.
 */
TransitionList::iterator findTrans_byQuantNumb(
		const string speciesLabel,
		const long n_hi, const long l_hi, const long S_hi,
		const long n_lo, const long l_lo, const long S_lo );

/** findTrans_byWLAng	Identify transition of a given species by wavelength.
 *
 * \param speciesLabel	(in) label of species, e.g., "H  1"
 * \param wavelength	(in) transition wavelength, in Angstrom
 * \param wl_err	(out) error in wavelength match
 *
 * \return	pointer to transition in species list of transitions. If no match is
 * 		found, return value is past the end of the entire list.
 */
TransitionList::iterator findTrans_byWLAng( string speciesLabel, const double wl_Ang,
			double &wl_err );


#endif /* TAULINES_H_ */
