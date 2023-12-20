/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*molcol generate and print molecular column densities */
#include "cddefines.h"
#include "molcol.h"
#include "radius.h"
#include "h2.h"
#include "mole.h"
#include "prt.h"

void molcol(
	const char *chLabel,
	/* file for printout */
	FILE *ioMEAN )
{
	long int i;

	DEBUG_ENTRY( "molcol()" );

	// Call column density routines for individual diatom models
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_Colden( chLabel );

	if( strcmp(chLabel,"ZERO") == 0 )
	{
		/*  zero out the column densities */
		for( i=0; i < mole_global.num_calc; i++ )
		{
			mole.species[i].column = 0.;
		}
	}

	else if( strcmp(chLabel,"ADD ") == 0 )
	{
		/*  add together column densities */
		for( i=0; i < mole_global.num_calc; i++ )
		{
			mole.species[i].column += (realnum)(mole.species[i].den*radius.drad_x_fillfac);
		}
	}

	else if( strcmp(chLabel,"PRIN") == 0 )
	{
		bool lgFirstPrint = true;
		/*  print the molecular column densities */
		int j=0;
		chem_nuclide *heavyAtom = null_nuclide;
		for( i=0; i < mole_global.num_calc; i++ )
		{
			if(mole.species[i].location == NULL && ( !mole_global.list[i]->isMonatomic() || !mole_global.list[i]->lgGas_Phase ) && mole.species[i].global()->n_react > 0 )
			{
				chem_nuclide* lastHeavyAtom = heavyAtom;
				heavyAtom = mole_global.list[i]->heavyAtom();
				int len = 11+CHARS_SPECIES;
				if(j+len > NCOLMAX || heavyAtom != lastHeavyAtom)
				{
					fprintf( ioMEAN, "\n" );
					if (heavyAtom != lastHeavyAtom)
					{
						fprintf ( ioMEAN, "==== %-*.*s compounds ====",
								CHARS_ISOTOPE_SYM, CHARS_ISOTOPE_SYM, heavyAtom->label().c_str() );
						if( lgFirstPrint )
						{
							fprintf ( ioMEAN, "           Log10 column densities [cm^-2]");
							lgFirstPrint = false;
						}
						fprintf(ioMEAN, "\n");
					}
					j = 0;
				}
				fprintf( ioMEAN, "   %-*.*s:", CHARS_SPECIES, CHARS_SPECIES, mole_global.list[i]->label.c_str() );
				fprintf( ioMEAN, "%7.3f",log10(MAX2(SMALLFLOAT,mole.species[i].column )));
				j+=len; 
			}
		}
		if (j != 0)
			fprintf( ioMEAN, "\n" );
	}

	else
	{
		fprintf( ioMEAN, " molcol does not understand the label %4.4s\n", 
		  chLabel );
		cdEXIT(EXIT_FAILURE);
	}
	return;

}
