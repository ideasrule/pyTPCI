/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef FE_H_
#define FE_H_

/**\file fe.h - vars that are for Fe only */
struct t_fe {

	/** fekhot, fekcld, number of photons in hot and cold iron, per unit vol */
	realnum fekhot, 
	  fekcld;

	/** Fe Ka from iron in grains */
	realnum fegrain;

	/** uv pumping of fe coronal lines */
	long int ipfe10;
	realnum pfe10, 
	  pfe11a, 
	  pfe11b, 
	  pfe14;

	};
extern t_fe fe;


#endif /* FE_H_ */
