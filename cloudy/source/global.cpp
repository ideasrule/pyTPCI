/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "global.h"
#include "dense.h"
#include "collision.h"

//
// The purpose of this file is to make global variable dependencies manifest.
// The order of variable initialization obviously matters.
//
// NB NB NB  Do not move the code below into a different file.
//

t_dense dense;
ColliderList colliders(dense);
