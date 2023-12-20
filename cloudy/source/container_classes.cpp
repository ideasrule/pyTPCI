/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "container_classes.h"

#ifdef _MSC_VER
/* disable "'extern' before template explicit instantiation" */
#	pragma warning( disable : 4231 )
#endif

/* Explicit instantiations for debugging purposes */
INSTANTIATE_MULTI_ARR( bool, lgBOUNDSCHECKVAL )
INSTANTIATE_MULTI_ARR( char, lgBOUNDSCHECKVAL )
INSTANTIATE_MULTI_ARR( int, lgBOUNDSCHECKVAL )
INSTANTIATE_MULTI_ARR( long, lgBOUNDSCHECKVAL )
#ifndef FLT_IS_DBL
INSTANTIATE_MULTI_ARR( realnum,lgBOUNDSCHECKVAL )
#endif
INSTANTIATE_MULTI_ARR( double, lgBOUNDSCHECKVAL )
