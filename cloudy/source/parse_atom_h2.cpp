/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseDatabaseH2 parse information from the atom command line */
#include "cddefines.h" 
#include "hmi.h" 
#include "h2.h" 
#include "parser.h" 

/*ParseDatabaseH2 parse information from the atom command line */
void ParseDatabaseH2(Parser &p )
{
	long int j;

	DEBUG_ENTRY( "ParseDatabaseH2()" );

	fixit("this must be generalized!!!");
		// easiest way is to create and populate diatoms before parsing,
		// so we can simply iterator of diatoms and match strings.
		// probably will want to put label on command line in quotes
		// and change command to, for example,
		// diatom "H2"
		// do that work before entering this routine.
	diatomics *diatom = NULL;

	if( p.nMatch(" H2 " ) )
	{
		diatom = &h2;
		/* this command has a 2 in the H2 label - must not parse the two by
		 * accident.  Get the first number off the line image, and confirm that
		 * it is a 2 */
		j = (long int)p.FFmtRead();
		if( j != 2 )
		{
			fprintf( ioQQQ, " Something is wrong with the order of the numbers on this line.\n" );
			fprintf( ioQQQ, " The first number I encounter should be a 2.\n Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	else
	{
		TotalInsanity();
	}

	/* the mere calling of this routine turns the large H2 molecule on */
	diatom->lgEnabled = true;

	if( p.nMatch("LEVE") )
	{
		/* number of electronic levels */

		/* lgREAD_DATA is false at start of calculation, set true when 
		 * space allocated for the H lines.  Once done we must ignore all 
		 * future changes in the number of levels */
		if( !diatom->lgREAD_DATA )
		{
			diatom->n_elec_states = (long int)p.FFmtRead();
			if( p.lgEOL() )
			{
				if( p.nMatch("LARG") )
				{
					/* LARGE is option to use the most number of electronic levels */
					diatom->n_elec_states = N_ELEC;
				}
				else
				{
					p.NoNumb("number of electronic levels");
				}
			}

			/* do not allow fewer than 3 - that includes Lyman & Werner bands */
			if( diatom->n_elec_states < 3 )
			{
				fprintf( ioQQQ, " This would be too few electronic levels - resetting to 3.\n" );
				diatom->n_elec_states = 3;
			}
			/* N_ELEC is the greatest number of elec lev possible */
			else if( diatom->n_elec_states > N_ELEC )
			{
				fprintf( ioQQQ, 
					" This would be too many levels, the limit is %i.\n" , 
					N_ELEC);
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	else if( p.nMatch("LIMI") )
	{
		/* the limit to the H2 / Htot ratio - 
		 * if smaller than this, do not compute large H2 mole */
		diatom->H2_to_H_limit = p.FFmtRead();
		if( p.lgEOL() )
		{
			/* did not find a number, either mistake or key " off" */
			if( p.nMatch( " OFF"  ) )
			{
				/* turn off limit */
				diatom->H2_to_H_limit = -1.;
			}
			else
			{
				p.NoNumb( "limit to the H2 / Htot ratio" );
			}
		}
		else
		{
			/* got a number, check if negative and so a log */
			/* a number <= 0 is the log of the ratio */
			if( diatom->H2_to_H_limit <= 0. )
				diatom->H2_to_H_limit = exp10( diatom->H2_to_H_limit);
		}
	}
	else if( p.nMatch("GBAR" ) )
	{
		/* option to either use, or not use, gbar approximation for low X 
		 * levels with no collision data - by default it is on */
		if( p.nMatch(" OFF" ) )
		{
			diatom->lgColl_gbar = false;
		}
		else if( p.nMatch(" ON " ) )
		{
			diatom->lgColl_gbar = true;
		}
		else
		{
			fprintf( ioQQQ, 
				" The gbar approximation must be off (\" OFF\") or on (\" ON \").\n");
			cdEXIT(EXIT_FAILURE);
		}
	}
	/* option to turn collisional effects off or on */
	else if( p.nMatch("COLL" ) )
	{
		/* option to turn collisional dissociation off or on */
		if( p.nMatch("DISS" ) )
		{
			/* option to turn collisions off */
			if( p.nMatch(" ON " ) )
			{
				/* this is the default, leave collisions off */
				diatom->lgColl_dissoc_coll = true;
			}
			else
			{
				/* default (and only reason for this command) is to turn off collisions */
				diatom->lgColl_dissoc_coll = false;
			}
		}
		/* option to turn collisional dissociation off or on 
		 * >>chng 06 mar 01, had been simply if - so all collisions were turned off
		 * when dissociation collisions turned off - 
		 * due to bucket else at end */
		else if( p.nMatch("ORTH" ) && p.nMatch("PARA" ) )
		{
			/* option to turn ortho - para collisions with particles off */
			if( p.nMatch(" ON " ) )
			{
				/* this is the default, leave collisions off */
				diatom->lgH2_ortho_para_coll_on = true;
			}
			else
			{
				/* default (and only reason for this command) is to turn off 
				 * ortho-para collisions */
				diatom->lgH2_ortho_para_coll_on = false;
			}
		}

		/* option to turn collisional effects off or on */
		else if( p.nMatch("GRAI" ) )
		{
			/* option to turn collisions off */
			if( p.nMatch(" ON" ) )
			{
				/* this is the default, leave collisions off */
				diatom->lgH2_grain_deexcitation = true;
			}
			else
			{
				/* default (and only reason for this command) is to turn off collisions */
				diatom->lgH2_grain_deexcitation = false;
			}
		}
		else if( p.nMatch(" H ") )
		{
			long int nYear = (long) p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("collision data set (year)");
			if (nYear == 1999)
			{
				/* use the coefficients from
				 *>>refer	H2	H collision	Le Bourlot, J., Pineau des Forets,
				 *>>refercon	G., & Flower, D.R. 1999, MNRAS, 305, 802 */
				h2.coll_source[0].filename = "coll_rates_H_99.dat";
			}
			else if (nYear == 2007)
			{
				/* use the coefficients from
				 *>>refer	H2	H collision	Wrathmall, S. A., Gusdorf A.,
				 *>>refercon	& Flower, D.R. 2007, MNRAS, 382, 133 */
				h2.coll_source[0].filename = "coll_rates_H_07.dat";
			}
			else if (nYear == 2015)
			{
				/* use the coefficients from
				 *>>refer	H2	H collision	Lique, F., 2015, MNRAS, 453, 810 */
				h2.coll_source[0].filename = "coll_rates_H_15.dat";
			}
			else 
			{
				/* not an option */
				fprintf(ioQQQ,
					" I did not find the dataset year on this DATABASE H2 H command"
					" - I know the 1999, 2007, and 2015 datasets" );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else if( p.nMatch(" HE " ) )
		{
			if( diatom != &h2 )
			{
				fprintf( ioQQQ, " Sorry. This command only applies to H2.  The keyword H2 must appear on the command.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* atom H2 He collisions ORNL (the default), Le BOURlot, and OFF
			 * which data set for He collisions,
			 * Teck Lee et al. ApJ to be submitted */
			if( p.nMatch(" NEW" ) || p.nMatch("ORNL" ) )
			{
				/* use the new coefficients */
				diatom->lgH2_He_ORNL = true;
				diatom->coll_source[1].filename = "coll_rates_He_ORNL.dat";
			}
			else if( p.nMatch(" OLD" ) || p.nMatch("BOUR" ) )
			{
				/* use the coefficients from
				 *>>refer	H2	collision	Le Bourlot, J., Pineau des Forets, 
				 *>>refercon	G., & Flower, D.R. 1999, MNRAS, 305, 802*/
				diatom->lgH2_He_ORNL = false;
				diatom->coll_source[1].filename = "coll_rates_He_LeBourlot.dat";
			}
			else
			{
				fprintf( ioQQQ, 
					" I did not find a keyword on this DATABASE H2 HE command - I know about the keys ORNL and Le BOURlot\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		/*>>chng 08 feb 27, GS*/
		else if( p.nMatch("ORH2" ) )
		{
			if( diatom != &h2 )
			{
				fprintf( ioQQQ, " Sorry. This command only applies to H2.  The keyword H2 must appear on the command.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* atom H2 H2ortho collisions ORNL (the default), Le BOURlot, and OFF
			 * which data set for H2 collisions,
			 * Teck Lee et al. ApJ to be submitted */
			if( p.nMatch("ORNL" ) )
			{
				/* use the new coefficients */
				diatom->lgH2_ORH2_ORNL = true;
				diatom->coll_source[2].filename = "coll_rates_H2ortho_ORNL.dat";
			}
			else if( p.nMatch("BOUR" ) )
			{
				/* use the coefficients from
				 *>>refer	H2	collision	Le Bourlot, J., Pineau des Forets, 
				 *>>refercon	G., & Flower, D.R. 1999, MNRAS, 305, 802*/
				diatom->lgH2_ORH2_ORNL = false;
				diatom->coll_source[2].filename = "coll_rates_H2ortho_LeBourlot.dat";
			}
			else
			{
				fprintf( ioQQQ, 
					" I did not find a keyword on this SPECIES H2 ohH2 command - I know about the keys ORNL and Le BOURlot\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		else if( p.nMatch("PAH2" ) )
		{
			/* atom H2 H2ortho collisions ORNL (the default), Le BOURlot, and OFF
			 * which data set for He collisions,
			 * Teck Lee et al. ApJ to be submitted */
			if(  p.nMatch("ORNL" ) )
			{
				/* use the new coefficients */
				diatom->lgH2_PAH2_ORNL = true;
				diatom->coll_source[3].filename = "coll_rates_H2para_ORNL.dat";
			}
			else if( p.nMatch("BOUR" ) )
			{
				/* use the coefficients from
				 *>>refer	H2	collision	Le Bourlot, J., Pineau des Forets, 
				 *>>refercon	G., & Flower, D.R. 1999, MNRAS, 305, 802*/
				diatom->lgH2_PAH2_ORNL = false;
				diatom->coll_source[3].filename = "coll_rates_H2para_LeBourlot.dat";
			}
			else
			{
				fprintf( ioQQQ, 
					" I did not find a keyword on this SPECIES H2 paH2 command - I know about the keys ORNL and Le BOURlot\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		else
		{
			/* option to turn all collisions off */
			if( p.nMatch(" ON " ) )
			{
				/* this is the default, leave collisions on */
				diatom->lgColl_deexec_Calc = true;
			}
			else
			{
				/* default (and only reason for this command) is to turn off collisions */
				diatom->lgColl_deexec_Calc = false;
			}
		}
	}

	/* set number of levels in matrix, but not trace matrix option */
	else if( p.nMatch("MATR" ) && !p.nMatch("TRAC" ) )
	{
		/* matrix option sets the number of levels that will
		 * be included in the matrix solution */
		long numLevels = (long)p.FFmtRead();
		if( p.nMatch(" ALL") )
		{
			/* " all" means do all of X, but space has not yet been allocated,
			 * so we do know know how many levels are within X - set special
			 * flag that will be used then this is known */
			numLevels = -1;
		}
		else if( p.lgEOL() && !(p.nMatch(" OFF") || p.nMatch("NONE") ) )
		{
			/* this branch hit eol but OFF or NONE is not on line - this is a mistake */
			fprintf( ioQQQ, 
				" The total number of levels used in the matrix solver must be entered, or keywords ALL or NONE entered.\n Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}

		diatom->set_numLevelsMatrix( numLevels );
		/* cannot check less than total number of levels within x since not yet set 
		 * We do not certify that matrix limits are greater than 1 -
		 * zero or <0 limits just turns if off, as did the off option */
	}
	else if( p.nMatch(" LTE" ) )
	{
		/* LTE option causes code to assume LTE for level populations  */
		diatom->lgLTE = true;
	}

	else if( p.nMatch("TRAC" ) )
	{
		/* these are used to set trace levels of output 
		diatom->n_trace_final = 1;
		diatom->n_trace_iterations = 2;
		diatom->n_trace_full = 3;
		diatom->n_trace_matrix = 4*/

		/* turns on trace printout - there are multiple levels */
		if( p.nMatch("FINA" ) )
		{
			/* FINAL gives only final information when solver exits */
			diatom->nTRACE = diatom->n_trace_final;
		}
		else if( p.nMatch("ITER" ) )
		{
			/* follow iterations within each call */
			diatom->nTRACE = diatom->n_trace_iterations;
		}
		else if( p.nMatch("FULL" ) )
		{
			/* full details of solution - this is also the default*/
			diatom->nTRACE = diatom->n_trace_full;
		}
		else if( p.nMatch("MATR" ) )
		{
			/* print the matrices used for X */
			diatom->nTRACE = diatom->n_trace_matrix;
		}
		else
		{
			/* full details of solution is also the default*/
			diatom->nTRACE = diatom->n_trace_full;
		}
	}
	else if( p.nMatch("NOIS" ) )
	{
		/* check on effects of uncertainties in collision rates */
		diatom->lgH2_NOISE = true;
		diatom->lgH2_NOISECOSMIC = true;

		/* optional mean - default is 0 */
		diatom->xMeanNoise = p.FFmtRead();
		if( p.lgEOL() )
			diatom->xMeanNoise = 0.;

		/* this is the standard deviation for the mole, with default */
		diatom->xSTDNoise = p.FFmtRead();
		if( p.lgEOL() )
			diatom->xSTDNoise = 0.5;
	}

	else if( p.nMatch("THER" ) )
	{
		/* change the treatment of the heating - cooling effects of H2,
		 * options are simple (use TH85 expressions) and full (use large molecule)*/
		if( p.nMatch("SIMP" ) )
		{
			hmi.lgH2_Thermal_BigH2 = false;
		}
		else if( p.nMatch("FULL" ) )
		{
			/* this is the default - use big atom */
			hmi.lgH2_Thermal_BigH2 = true;
		}
	}

	else if( p.nMatch("CHEM" ) )
	{
		/* atom h2 chemistry simple command
		 * change the treatment of the chemistry - formation and destruction,
		 * options are simple (use TH85 expressions) and full (use large molecule)*/
		if( p.nMatch("SIMP" ) )
		{
			hmi.lgH2_Chemistry_BigH2 = false;
		}
		else if( p.nMatch("FULL" ) )
		{
			/* this is the default - use big atom */
			hmi.lgH2_Chemistry_BigH2 = true;
		}
	}

	/* there is no final branch - if we do not find a keyword, simply
	 * turn on the H2 molecule */
	return;
}
