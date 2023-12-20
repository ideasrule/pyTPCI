/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "elementnames.h"
t_elementnames elementnames;

/* elem_symbol_to_index - obtain the array index the input element symbol refers to.
 * The first letter is expected to be uppercase, and the second lowercase. */
int elem_symbol_to_index( const string& chSym )
{
	int ielem = -1;

	size_t symlen = chSym.length();
	if( symlen == 0 || symlen >= size_t(CHARS_ELEMENT_SYM) )
		return ielem;

	for( int i = 0; i < LIMELM; ++i )
	{
		if( symlen == 1 )
		{
			if( toupper( chSym[0] ) == elementnames.chElementSym[i][0] &&
				elementnames.chElementSym[i][1] == ' ' )
			{
				ielem = i;
				break;
			}
		}
		else if( toupper( chSym[0] ) == elementnames.chElementSym[i][0] &&
				 tolower( chSym[1] ) == elementnames.chElementSym[i][1] )
		{
			ielem = i;
			break;
		}
	}
	return ielem;
}

/* isElementSym - search for string among element symbols. */
bool isElementSym( const string& chSym )
{
	return ( elem_symbol_to_index( chSym ) != -1 ) ? true : false;
}
