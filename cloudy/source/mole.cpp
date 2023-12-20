/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "mole.h"
#include "parser.h"

t_mole_global mole_global;
t_mole_local mole;

void t_mole_global::zero(void)
{
	DEBUG_ENTRY( "t_mole_global::zero()" );
	/* flag to turn off molecular network */
	lgNoMole = false;
	lgNoHeavyMole = false;
	/* capture of molecules onto grain surfaces - formation of ices
	 * flag says to include this process - turned off with the
	 * command NO GRAIN MOLECULES */
	lgGrain_mole_deplete = true;
	/* flag saying that H2O water destruction rate went to zero */
	lgH2Ozer = false;

	/* option to turn on the original UMIST rates rather than the corrected network
	 * adopted in the Leiden PDR meeting, by deffault false,
	 * set true with the set UMIST rates command,
	 * these rates are from mole_co_umisthack.dat and also in mole_co_base.dat (UMIST)*/
	lgLeidenHack = false;

	/* option to use diffuse cloud chemical rates from Table 8 of
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is true - changed with set chemistry command
	 * these rates are in both mole_co_federman.dat and mole_co_base.dat,
	 * latter based on UMIST */
	lgFederman = true;

	/* option to use effective temperature as defined in
	 * >> refer Zsargo, J. & Federman, S. R. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	lgNonEquilChem = false;
	/** option to set proton elimination rates to zero
	 * >>refer	CO	chemistry	Huntress, W. T., 1977, ApJS, 33, 495
	 * By default, this is false - changed with set chemistry command */
	lgProtElim = true;
	/** option to not include neutrals in the non-equilibrium scheme
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	lgNeutrals = true;
	/* option to use H2 continuum dissociation cross sections computed by P.C. Stancil
	 * By default, this is true - changed with "set H2 continuum dissociation xxx" command
	 * options are "Stancil" or "AD69" */
	lgStancil = false;
	// all isotopes are currently disabled by default
	lgTreatIsotopes.resize( LIMELM );
	fill( lgTreatIsotopes.begin(), lgTreatIsotopes.end(), false );
	
}
/*=================================================================*/
/*mole_Init called from cdInit to initialize CO routines */
void t_mole_global::init(void)
{
	DEBUG_ENTRY( "Mole::init()" );

	/* prevent memory leaks */
	/* \todo	this is a temporary fix for PR14. We should improve the overall design
	 * of this code to prevent valid pointers being overwritten in a second call to mole_Init */
	static bool lgmole_Init_called = false;
	static long int num_total_allocated = -1;
	if(! lgmole_Init_called )
	{
		/* say that we have been called */
		lgmole_Init_called = true;
		
		make_species();
		mole_make_list();
		mole_make_groups();
		mole.species.resize( num_total );
		num_total_allocated = num_total;
	}

	if( num_total > num_total_allocated )
	{
		/* number of species has increased since last time - this can't happen
		 * tsuite / programs / comp4 has 95 first time, 98 second time */
		fprintf(ioQQQ,"DISASTER - the number of species in the chemistry network has increased.  This is not allowed.\n");
		fprintf(ioQQQ,"This could happen if an element was initially turned off or grains not included, then the element or grains was included.  There are not allowed.\n");
		fprintf(ioQQQ,"Sorry.\n");
		cdEXIT(EXIT_FAILURE);
	}

	for( long i=0; i<num_total; ++i ) 
	{
		mole.species[i].zero();
		mole.species[i].index = i;
	}
	mole.elec = 0.;
}

void t_mole_local::alloc()
{
	source.reserve(LIMELM);
	xMoleChTrRate.reserve(LIMELM);
	for( long nelem=0; nelem < LIMELM; ++nelem )
	{		
		/* chemistry source and sink terms for ionization ladders */
		source.reserve(nelem, nelem+2);
		xMoleChTrRate.reserve(nelem, nelem+2);
		for( long ion=0; ion < nelem+2; ++ion )
		{
			xMoleChTrRate.reserve(nelem, ion, nelem+2);
		}
	}	
	/* these are source and sink terms for heavy element ionization balance from the
	 * chemistry */
	source.alloc();
	sink.alloc(source.clone());
	xMoleChTrRate.alloc();
}

void t_mole_local::zero()
{
	source = 0.;
	sink = 0.;
	xMoleChTrRate = 0.;
	for( long i=0; i < mole_global.num_calc; i++ )
	{
		species[i].column = 0.;
	}
	fill( reaction_rks.begin(), reaction_rks.end(), 0. );
	fill( old_reaction_rks.begin(), old_reaction_rks.end(), 0. );
}

void ParseChemistry(Parser&p)
{
	DEBUG_ENTRY( "ParseChemistry()" );

	Symbol s = p.getSymbol();
	if (s.toktype == Symbol::EOSTAT || 
		 s.toktype == Symbol::ERROR)
	{
		fprintf(ioQQQ,"Error, CHEMISTRY requires option to be specified\n");
		cdEXIT(EXIT_FAILURE);
	}
	
	for(;;)
	{
		if (s.toktype == Symbol::NAME)
		{
			if (!strncmp("REAC",s.value.c_str(),4))
			{
				s = p.getSymbol();
				vector<string> reactions;
				if (s.toktype == Symbol::STRING)
				{
					reactions.push_back(s.value);
				}
				else if (s.toktype == Symbol::OPERATOR && s.value == "(")
				{
					for (;;)
					{
						s = p.getSymbol();
						if (s.toktype != Symbol::STRING)
						{
							fprintf(ioQQQ,"No reactions found for CHEMistry REACtion command\n");
							fprintf(ioQQQ,"Reactions needs to be included in quotes \"\"\n");
							cdEXIT(EXIT_FAILURE);
						}
						reactions.push_back(s.value);
						s = p.getSymbol();
						if (s.toktype == Symbol::OPERATOR)
						{
							if (s.value == ")")
							{
								break;
							}
							else if (s.value != ",")
							{
								fprintf(ioQQQ,"Syntax error for CHEMistry REACtion command\n");
								cdEXIT(EXIT_FAILURE);
							}
						}
						else
						{
							fprintf(ioQQQ,"Syntax error for CHEMistry REACtion command\n");
							cdEXIT(EXIT_FAILURE);
						}
					}
				}
				else
				{
					fprintf(ioQQQ,"No reaction found for CHEMistry REACtion command\n");
					fprintf(ioQQQ,"Reaction needs to be included in quotes \"\"\n");
					cdEXIT(EXIT_FAILURE);
				}
				s = p.getSymbol();
				if (s.toktype == Symbol::NAME)
				{
					if (s.value == "OFF")
					{
						for (vector<string>::iterator it = reactions.begin();
							  it != reactions.end(); ++it)
						{
							mole_global.offReactions[*it] = true;
						}
					}
					else
					{
						fprintf(ioQQQ,"Error, CHEMISTRY REACTION option %s "
								  "not known\n",s.value.c_str());
						cdEXIT(EXIT_FAILURE);
					}
				}
				else
				{
					fprintf(ioQQQ,"Error, CHEMISTRY REACTION needs action\n");
					cdEXIT(EXIT_FAILURE);
				}
			}
		}
		else
		{
			fprintf(ioQQQ,"Error, unknown option to CHEMISTRY: %s\n",
					  s.value.c_str());
			cdEXIT(EXIT_FAILURE);
		}
		s = p.getSymbol();
		if (s.toktype == Symbol::EOSTAT || 
			 s.toktype == Symbol::ERROR)
		{
			break;
		}
		if (s.toktype != Symbol::OPERATOR || s.value != ",")
		{
			fprintf(ioQQQ,"Syntax error, CHEMISTRY option [, option]\n");
			cdEXIT(EXIT_FAILURE);
		}
		s = p.getSymbol();
	}
}

