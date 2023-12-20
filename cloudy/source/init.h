/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef INIT_H_
#define INIT_H_

/* initialization routines for Cloudy
 * InitCoreload - one time initialization, called from cdInit
 * InitCoreloadPostparse - one time initialization, called after parser
 * InitDefaultsPreparse initialization at start of simulation, called from cloudy
 * before parser, will be called one time per sim in grid 
 * InitSimPostparse initialization at start of simulation, called from cloudy
 * after parser, will be called one time per sim in grid
 */

/** one time initialization of core load, called from cdDrive, this sets
 * minimum set of values needed for the code to start - called after
 * input lines have been read in and checked for VARY or GRID - so
 * known whether single or multiple sims will be run */
void InitCoreload( void );

/** initialize values that are changed in the parser.  Called just
 * before parser, will be called one time per core load */
void InitDefaultsPreparse( void );

/** initialization after parser, called one time per core load
 * create space needed for code to operate */
void InitCoreloadPostparse( void );

/** initialize values at start of simulation, called after parser, 
 * sets initial or zero values at start of each sim in grid */
void InitSimPostparse( void );

/**zero actively zero out or initialize variables needed for model calculation 
 * this is the old one and should be removed - its vars moved into those above */
void zero(void);

#endif /* INIT_H_ */
