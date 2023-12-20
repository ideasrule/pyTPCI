/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef H2_H_
#define H2_H_

#include "h2_priv.h"

extern diatomics h2;
extern diatomics hd;

extern vector<diatomics*> diatoms;
typedef vector<diatomics*>::iterator diatom_iter;

void diatoms_init( void );	

#endif /* H2_H_ */
