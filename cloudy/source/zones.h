/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ZONES_H_
#define ZONES_H_

/**ZoneStart, ZoneEnd, set variables that change with each zone, like radius, depth 
\param *chMode
*/
void ZoneStart( const char *chMode);
/**ZoneStart, ZoneEnd, set variables that change with each zone, like radius, depth
*/
void ZoneEnd(void);

#endif /* ZONES_H_ */
