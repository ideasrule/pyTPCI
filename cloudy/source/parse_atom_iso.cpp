/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseDatabaseISO parse information from the atom XX-like command line */
#include "cddefines.h"
#include "elementnames.h"
#include "optimize.h"
#include "hydrogenic.h"
#include "input.h"
#include "iso.h"
#include "parser.h"
#include "phycon.h"
#include "rfield.h"
#include "taulines.h"

/*ParseDatabaseISO parse parameters off the XX-like command */
void ParseDatabaseISO(long ipISO, Parser &p )
{
	long int numLevels;

	DEBUG_ENTRY( "ParseDatabaseISO()" );

	/* look for the name of an element - if we don't find one do the entire
	 * iso sequence - returns negative number if element not found */
	long int nelem = p.GetElem( );

	/* He-like Hydrogen is not possible */
	if( ipISO==ipHE_LIKE && nelem==ipHYDROGEN )
	{
		fprintf(ioQQQ," Sorry, He-like H is unacceptable.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* collisions - don't pick up the levels collapsed command */
	if( p.nMatch("COLL") && !p.nMatch("LEVE"  ) )
	{
		if( p.nMatch("THER") )
		{
			/* This is the "Thermal average" command.  It
			 * causes collisions strengths to be integrated over
			 * a Maxwellian peaked at kT rather than to be
			 * evaluated at kT
			 */
			iso_ctrl.lgCS_therm_ave[ipISO] = true;
		}


		/* option to turn collisions off, all are
		 * set to 1 in zero.  command can accept only one option at a time */
		else if( p.nMatch("EXCI") )
		{
			/* turn off collisional excitation */
			iso_ctrl.lgColl_excite[ipISO] = false;
			phycon.lgPhysOK = false;
		}
		else if( p.nMatch("IONI") )
		{
			/* turn off collisional ionization */
			iso_ctrl.lgColl_ionize[ipISO] = false;
			phycon.lgPhysOK = false;
		}

		else if( p.nMatch("2S2P") || ( p.nMatch("2P2S") && ipISO == ipH_LIKE ) )
		{
			/* >>chng 02 feb 07, change from 2s2p to l-mixing */
			/* this is the atom h-like collision l-mixing command */
			fprintf(ioQQQ,"This command changed to DATABASE H-LIKE COLLISIONS L-MIXING\n");
			fprintf(ioQQQ,"I will parse it for now, but may not in the future.\n");
			/* turn off 2s - 2p collisions */
			iso_ctrl.lgColl_l_mixing[ipISO] = false;
			phycon.lgPhysOK = false;
		}

		else if ( p.nMatch("LEBED"))
		{
			iso_ctrl.lgCS_Lebedev[ipISO] = true;
			iso_ctrl.lgCS_Vriens[ipISO] = false;
			iso_ctrl.lgCS_Fujim[ipISO] = false;
			iso_ctrl.lgCS_vrgm[ipISO] = false;

		}

		else if( p.nMatch("VRIE") )
		{
			/* option to change how collisions for unknown levels affect things */
			iso_ctrl.lgCS_Lebedev[ipISO] =false;
			iso_ctrl.lgCS_Vriens[ipISO] = true;
			iso_ctrl.lgCS_Fujim[ipISO] = false;
			iso_ctrl.lgCS_vrgm[ipISO] = false;

		}
		else if ( p.nMatch("FUJIM"))
		{
			iso_ctrl.lgCS_Fujim[ipISO] = true;
			iso_ctrl.lgCS_Lebedev[ipISO] =false;
			iso_ctrl.lgCS_Vriens[ipISO] = false;
			iso_ctrl.lgCS_vrgm[ipISO] = false;


		}
		else if ( p.nMatch("VAN R"))
		{
			iso_ctrl.lgCS_vrgm[ipISO] = true;
			iso_ctrl.lgCS_Lebedev[ipISO] =false;
			iso_ctrl.lgCS_Vriens[ipISO] = false;
			iso_ctrl.lgCS_Fujim[ipISO] = false;


		}


		else if( p.nMatch("L-MI") )
		{
			/* criteria for degeneration, if criteria are applied
			 * different cross section are obtained for nearly degenarate l and l'.
			 * for B72 refer Brocklehurst MNRAS (1972) 157,211
			 * for PS_deg Refer Penguelly & Seaton MNRAS (1964) 127,165
			 */

			if ((p.nMatch(" ALL") && p.nMatch(" DEG")) || ipISO == ipH_LIKE)
			{
				/* All l's are treated as degenerate
				 * used only for comparisons
				 */
				iso_ctrl.lgCS_Seaton[ipISO] = false;
				iso_ctrl.lgCS_B72[ipISO] = false;
				iso_ctrl.lgCS_PSdeg[ipISO]=false;
				iso_ctrl.lgCS_PSM20[ipISO]=false;
			}
			else if (ipISO == ipHE_LIKE)
			{
				if (p.nMatch(" S62"))
				{
					if (p.nMatch(" OFF"))
					{
						iso_ctrl.lgCS_Seaton[ipISO] = false;
					}
					else
					{
						/* electron collisions dominate for l<3
						 * Proton impact collisions are calculated using formalism from
						 * Seaton 1962 Proc. Phys Soc. 79, 1105
						 * S62 is compatible with other cut offs
						 * DEFAULT
						 */
						iso_ctrl.lgCS_Seaton[ipISO] = true;
					}
				}
				if (p.nMatch(" B72"))
				{
					/* Uses Brocklehurst 1972 MNRAS157, 211
					 * equation 3.12 criterion for degeneracy
					 */
					iso_ctrl.lgCS_B72[ipISO] = true;
					iso_ctrl.lgCS_PSdeg[ipISO]=false;
				}
				else if	(p.nMatch(" NO ") && p.nMatch(" DEG"))
				{
					/* Follows criterion for degeneration
					 * from Pengelly and Seaton (1964) MNRAS 127,165
					 * in eq 23 which has been relaxed to be DE>=hbar/tau
					 */
					iso_ctrl.lgCS_B72[ipISO] = false;
					iso_ctrl.lgCS_PSdeg[ipISO] = true;
				}
				else if(p.nMatch("PSM17"))
					iso_ctrl.lgCS_PSM20[ipISO]=false;
				else if(p.nMatch("PSM20"))
					{
					/*
					 * Default option from Badnell 2021.
					 */
					iso_ctrl.lgCS_PSM20[ipISO]=true;
					}
			}

			/* use l-mix from
			 *>>refer	l-mix	all Vrinceanu, D. & Flannery, M. R. 2001, PhysRevA 63, 032701
			 */
			if( p.nMatch("VRIN") )
			{
				iso_ctrl.lgCS_Vrinceanu[ipISO] = true;
				iso_ctrl.lgCS_PS64[ipISO] = false;
				iso_ctrl.lgCS_VOS12[ipISO] = false;
				iso_ctrl.lgCS_VOS12QM[ipISO]=false;
				if (p.nMatch("THER") || iso_ctrl.lgCS_therm_ave[ipISO])
					iso_ctrl.lgCS_VOS_thermal[ipISO] = true;
			}
			else if ( p.nMatch("VOS12") )
			{
				if (p.nMatch("SEMIC"))
				{
				/* Vrinceanu+ 2012 equation (7) semiclassical */
					iso_ctrl.lgCS_Vrinceanu[ipISO] = false;
					iso_ctrl.lgCS_PS64[ipISO] = false;
					iso_ctrl.lgCS_VOS12[ipISO] = true;
					iso_ctrl.lgCS_VOS12QM[ipISO]=false;
				}
				else //if( p.nMatch("QUANT"))
				{
					/* Vrinceanu+ 2012 equation (2) quantal */
					iso_ctrl.lgCS_Vrinceanu[ipISO] = false;
					iso_ctrl.lgCS_PS64[ipISO] = false;
					iso_ctrl.lgCS_VOS12[ipISO] = false;
					iso_ctrl.lgCS_VOS12QM[ipISO]=true;
					if (p.nMatch("THER") || iso_ctrl.lgCS_therm_ave[ipISO])
						iso_ctrl.lgCS_VOS_thermal[ipISO] = true;
				}
			}
			else if( p.nMatch("PENG") )
			{
				/* Pengelly & Seaton for l-mixing
				 * N. Badnell et al. MNRAS (2021) 507, 2922
				 * THAT'S IS THE DEFAULT
				 */
				iso_ctrl.lgCS_Vrinceanu[ipISO] = false;
				iso_ctrl.lgCS_VOS12[ipISO]=false;
				iso_ctrl.lgCS_VOS12QM[ipISO]=false;
				iso_ctrl.lgCS_PS64[ipISO] = true;
				/*classic Pengelly and Seaton 1964*/
				if (p.nMatch("CLASS"))
					iso_ctrl.lgCS_PSClassic[ipISO] = true;
				/* PSM from Guzman+ 2016, 2017*/
				else if(p.nMatch("PSM17"))
					iso_ctrl.lgCS_PSM20[ipISO]=false;
				else if(p.nMatch("PSM20"))
				{
					/*
					 * Default option from Badnell 2021.
					 */
					iso_ctrl.lgCS_PSM20[ipISO]=true;
				}
			}
			else if( p.nMatch(" OFF"  ) )
			{
				/* this is the atom xx-like collision l-mixing command */
				/* turn off same-n collisions */
				iso_ctrl.lgColl_l_mixing[ipISO] = false;
				phycon.lgPhysOK = false;
				iso_ctrl.lgCS_Vrinceanu[ipISO] = false;
				iso_ctrl.lgCS_VOS12[ipISO]=false;
				iso_ctrl.lgCS_VOS12QM[ipISO]=false;
				iso_ctrl.lgCS_PSM20[ipISO]=false;
			}
			else
			{
				fprintf( ioQQQ, "The database H-like l-mixing command needs a keyword\n"
						" Options are OFF, PENGelly, VRINCeanu, VOS12 (SEMIClassical or Quantal), PSM.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
		else if( p.nMatch(" OFF"  ) )
		{
			/* turn everything off, since no keyword given */
			iso_ctrl.lgColl_excite[ipISO] = false;
			iso_ctrl.lgColl_ionize[ipISO] = false;
			iso_ctrl.lgColl_l_mixing[ipISO] = false;
			phycon.lgPhysOK = false;
		}
		else
		{
			fprintf( ioQQQ, " needs parameter\n"
					" Options are OFF, THERmal, EXCItation off, IONIzation off,2s2p off, l-mixing [options].\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("CONT") && p.nMatch("LOWE") )
	{
		/* disable continuum lowering for this isoelectronic sequence */
		if( p.nMatch("OFF") )
			iso_ctrl.lgContinuumLoweringEnabled[ipISO] = false;
		else
			iso_ctrl.lgContinuumLoweringEnabled[ipISO] = true;
	}

	else if( p.nMatch("DAMP") )
	{
		if( ipISO == ipHE_LIKE )
		{
			fprintf(ioQQQ," Sorry, the DAMPING option is not implemented for the he-like sequence.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* turn off absorption due to Ly alpha damping wings */
		hydro.DampOnFac = 0.;
	}

	else if( p.nMatch("DIEL") )
	{
		if( ipISO == ipH_LIKE )
		{
			fprintf(ioQQQ," Sorry, but dielectronic recombination onto the h-like sequence is not possible.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* This sets which set of data to use for dielectronic recombination.	*/
		if( p.nMatch(" OFF") )
		{
			iso_ctrl.lgDielRecom[ipISO] = false;
		}
		else 
			iso_ctrl.lgDielRecom[ipISO] = true;
	}

	else if( p.nMatch("MASE") && p.nMatch("OFF"))
	{

		/* Do not allow masers in iso-sequences */
		if( nelem<0 )
		{
			// element not specified, do entire sequence
			for( long elms=ipISO; elms<LIMELM; ++elms )
			{
				iso_ctrl.lgNoMaser[ipISO][elms] = true;
			}
		}
		else
			iso_ctrl.lgNoMaser[ipISO][nelem] = true;

	}

	else if( p.nMatch("LEVE") )
	{
		/* the number of levels read in is n, the principal quantum number
		 * only lines with upper levels less than n will be printed */

		/* number of levels for iso-sequence */
		/* there are two options here,
		 * when keyword ELEMENT appears, scan off element name and change levels only for
		 * that one.
		 * when there is no ELEMENT then set all in iso to same number */

		/* lgHydroAlloc is false at start of calculation, set true when space 
		 * allocated for the hydrogen and helium lines.  Once done we must ignore all 
		 * future changes in the number of levels */
		if( p.nMatch("LTE") )
		{
			/* force level ratios to LTE */
			iso_ctrl.lgLTE_levels[ipISO] = true;
		}
		else if( p.nMatch("PRIN") )
		{
			/* only print - do not change levels */
			iso_ctrl.lgPrintNumberOfLevels = true;
		}
		else if( !lgHydroAlloc )
		{
			numLevels = (long int)p.FFmtRead();

			if( !p.lgEOL() )
			{
				if( ipISO == ipH_LIKE && numLevels > NHYDRO_MAX_LEVEL-2 )
				{
					fprintf( ioQQQ, " Not possible to set nhlvl to >NHYDRO_MAX_LEVEL-2= %i\n",
					  NHYDRO_MAX_LEVEL-2 );
					fprintf( ioQQQ, " change NHYDRO_MAX_LEVEL\n");
					cdEXIT(EXIT_FAILURE);
				}

				/* check that alpha transition of highest level is within energy bounds of continuum */
				if( !p.nMatch("COLL") && ipISO == ipH_LIKE &&
					( 2. / POW3((double)numLevels) < rfield.emm() ) )
				{
					fprintf( ioQQQ, " Not possible to set iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_max to such a high value, since "
						"alpha transition not within energy bounds of code\n");

					fprintf( ioQQQ, " lowest energy is %e and corresponding highest level is %li\n" ,
						rfield.emm(), (long)cbrt(2./rfield.emm()) );
					cdEXIT(EXIT_FAILURE);
				}
			}

			if( p.lgEOL() )
			{
				int LevelsResolved=-1 , LevelsCollapsed=10;
				/* no number, so check for either large or small */
				if( p.nMatch("LARG") )
				{
					/* includes all levels with tabulated rec coefficient */
					LevelsResolved = RREC_MAXN;
				}

				/* this is small or compact keyword */
				else if( p.nMatch("SMAL") || p.nMatch("COMP") )
				{
					if( ipISO == ipH_LIKE )
						LevelsResolved = 5;
					else if( ipISO == ipHE_LIKE )
						LevelsResolved = 3;
					else
						TotalInsanity();
				}
				else
					/* punch out if no number */
					p.NoNumb("levels");

				if( nelem<0 )
				{
					// element not specified, do entire sequence
					for( nelem=ipISO; nelem<LIMELM; ++nelem )
					{
						iso_sp[ipISO][nelem].nCollapsed_max =
								MIN2( iso_sp[ipISO][nelem].nCollapsed_max , LevelsCollapsed );
						iso_sp[ipISO][nelem].n_HighestResolved_max =
								MIN2( iso_sp[ipISO][nelem].n_HighestResolved_max , LevelsResolved );
						iso_update_num_levels( ipISO, nelem );
					}
				}
				else
				{
					iso_sp[ipISO][nelem].nCollapsed_max = LevelsCollapsed;
					iso_sp[ipISO][nelem].n_HighestResolved_max = LevelsResolved;
					iso_update_num_levels( ipISO, nelem );
				}
			}

			else if( p.nMatch("COLLAP") )
			{
				// set number of collapsed levels
				if( numLevels < 1 )
				{
					fprintf( ioQQQ, "There must be at least one collapsed level.\n");
					cdEXIT(EXIT_FAILURE);
				}

				if( nelem<0 )
				{
					// element not specified, do entire sequence
					for( nelem=ipISO; nelem<LIMELM; ++nelem )
					{
						iso_sp[ipISO][nelem].nCollapsed_max = numLevels;
						iso_update_num_levels( ipISO, nelem );
					}
				}
				else
				{
					iso_sp[ipISO][nelem].nCollapsed_max = numLevels;
					iso_update_num_levels( ipISO, nelem );
				}
			}
			else if( p.nMatch("RESOLV") )
			{
				// number of resolved levels
				if( ( numLevels < 3 ) && !p.nMatch("COLL") )
				{
					fprintf( ioQQQ, " cannot have fewer than 3 resolved levels, the requested number was %li\n" ,
						numLevels  );
					fprintf( ioQQQ, " Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}

				if( nelem<0 )
				{
					// element not specified, do entire sequence
					for( nelem=ipISO; nelem<LIMELM; ++nelem )
					{
						iso_sp[ipISO][nelem].n_HighestResolved_max = numLevels;
						iso_update_num_levels( ipISO, nelem );
					}
				}
				else
				{
					iso_sp[ipISO][nelem].n_HighestResolved_max = numLevels;
					iso_update_num_levels( ipISO, nelem );
				}
			}
			else
			{
				fprintf(ioQQQ, "I did not recognize a keyword on this atom xx-like levels command."
					"  Should be COLLAPSED or RESOLVED.\n Sorry.\n\n");
					cdEXIT(EXIT_FAILURE);
			}
		}
	}

	else if( p.nMatch("ERRO") && p.nMatch("GENE") )
	{
		/* Rates will be modified by a randomly generated error that falls within
		 * the range specifically set for each rate (or set of rates).	*/
		iso_ctrl.lgRandErrGen[ipISO] = true;

		if( p.nMatch("PESS") )
			iso_ctrl.lgPessimisticErrors = true;
		else
			iso_ctrl.lgPessimisticErrors = false;
	}

	else if( p.nMatch("GBAR") )
	{
		if( ipISO == ipH_LIKE )
		{
			fprintf(ioQQQ," Sorry, the GBAR option is only implemented for the He-like sequence.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* the HEGBAR command - to change cs of higher levels */
		/* first turn all off, one will be turned back on */
		iso_ctrl.lgCS_Vriens[ipISO] = false;
		iso_ctrl.lgCS_None[ipISO] = false;
		iso_ctrl.nCS_new[ipISO] = false;
		iso_ctrl.lgCS_Lebedev[ipISO]=false;
		iso_ctrl.lgCS_Fujim[ipISO] = false;
		iso_ctrl.lgCS_vrgm[ipISO] = false;

		/* now turn one on */
		if( p.nMatch("VRIE") )
		{
			/* option to change how collisions for unknown levels affect things */
			iso_ctrl.lgCS_Vriens[ipISO] = true;
		}
		else if( p.nMatch(" NEW") )
		{
			/* option to change how collisions for unknown levels affect things */
			iso_ctrl.nCS_new[ipISO] = (int)p.FFmtRead();

			/* there are two options, 1 and 2, for which fit - default (no number)
			 * will be 1, the broken power law fit */
			if( p.lgEOL() )
				iso_ctrl.nCS_new[ipISO] = 1;

			ASSERT( iso_ctrl.nCS_new[ipISO] );
		}
		else if( p.nMatch(" OFF") )
		{
			/* option to change how collisions for unknown levels affect things */
			iso_ctrl.lgCS_None[ipISO] = true;
		}
		else
		{
			fprintf( ioQQQ, " needs parameter\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("LYMA") )
	{
		if( ipISO == ipH_LIKE && p.nMatch("PUMP") )
		{
			/* >>chng 05 jul 08, separate out Lyman pump commands */
			if( p.nMatch(" OFF") )
			{
				/* option to turn off all continuum pumping of Lyman lines */
				hydro.lgLymanPumping = false;
			}
			else if( p.nMatch("SCALE") )
			{
				/* multiplicative factor for all continuum pumping of H I Lyman lines,
				 * account for possible emission in the line - only affects H I
				 * not entire H-like iso sequence */
				hydro.xLymanPumpingScaleFactor = 
					(realnum)p.FFmtRead();
				/* scale factor is log if <=0, 
				 * represents line in emission if >1
				 * LOG keyword forces interpretation as a log */
				if( hydro.xLymanPumpingScaleFactor <= 0. ||
					p.nMatch(" LOG") )
				{
					hydro.xLymanPumpingScaleFactor = 
						exp10( hydro.xLymanPumpingScaleFactor );
				}

				/* vary option */
				if( optimize.lgVarOn )
				{
					optimize.nvarxt[optimize.nparm] = 1;
					strcpy( optimize.chVarFmt[optimize.nparm], "DATABASE H-LIKE LYMAN PUMPING SCALE %f LOG" );

					/*  pointer to where to write */
					optimize.nvfpnt[optimize.nparm] = input.nRead;

					/*  current parameters - always log so steps are log  */
					optimize.vincr[optimize.nparm] = 0.1f;
					optimize.vparm[0][optimize.nparm] = (realnum)log10(hydro.xLymanPumpingScaleFactor);
					++optimize.nparm;
				}
			}
			else
			{
				fprintf(ioQQQ," Sorry, I didn\'t recognize an option on this DATABASE H-LIKE LYMAN PUMP command.\n");
				fprintf(ioQQQ," The options are \" OFF\", and \"SCALE\".\n");  
				cdEXIT(EXIT_FAILURE);
			}
		}
		else if( p.nMatch("EXTRA") )
		{
			/* option to set number of "extra" Lyman lines, used for optical depths only */
			iso_ctrl.nLyman_max[ipISO] = (long int)p.FFmtRead();
			iso_ctrl.nLyman[ipISO] = MIN2(iso_ctrl.nLyman[ipISO],
													iso_ctrl.nLyman_max[ipISO]);
			if( p.lgEOL() )
				p.NoNumb("'extra' Lyman lines");
			if( iso_ctrl.nLyman[ipISO] < 2 )
			{
				// Code does not elsewhere protect against values less than 2.
				fprintf(ioQQQ," Sorry, the value on this DATABASE xx-LIKE LYMAN command must be at least 2.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			fprintf(ioQQQ," Sorry, I didn\'t recognize an option on this DATABASE xx-LIKE LYMAN command.\n");
			fprintf(ioQQQ," The options are \"PUMP\", and \"EXTRA\".\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* don't interpolate on look-up tables but compute recombination on the fly instead */
	else if( p.nMatch("RECO") &&
		p.nMatch(" NO ") && p.nMatch("INTE") )
	{
		/* flag set by atom xx-like no recombination interp command, 
		 * says to generate recombination coefficients
		 * on the fly */
		iso_ctrl.lgNoRecombInterp[ipISO] = true;
	}

	else if( p.nMatch("REDI") )
	{
		int ipRedis=0;
		/* there are three functions, PRD_, CRD_, and CRDW,
		 * representing partial redistribution, 
		 * complete redistribution (doppler core only, no wings)
		 * and complete with wings */
		/* partial redistribution */
		if( p.nMatch(" PRD") )
		{
			ipRedis = ipPRD;
		}
		/* complete redistribution no wings */
		else if( p.nMatch(" CRD") )
		{
			ipRedis = ipCRD;
		}
		/* complete redistribution with wings */
		else if( p.nMatch("CRDW") )
		{
			ipRedis = ipCRDW;
		}

		/* if not SHOW option (handled below) then we have a problem */
		else if( !p.nMatch("SHOW") )
		{
			fprintf(ioQQQ," There should have been a second keyword on this command.\n");
			fprintf(ioQQQ," Options are _PRD, _CRD, CRDW (_ is space).  Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* resonance lines - not Lya*/
		if( p.nMatch("ALPH") )
		{
			iso_ctrl.ipLyaRedist[ipISO] = ipRedis;
		}
		/* Lya itself */
		else if( p.nMatch("RESO") )
		{
			iso_ctrl.ipResoRedist[ipISO] = ipRedis;
		}
		/* subordinate lines */
		else if( p.nMatch("SUBO") )
		{
			iso_ctrl.ipSubRedist[ipISO] = ipRedis;
		}
		/* the show option, say what we are assuming */
		else if( p.nMatch("SHOW") )
		{
			fprintf(ioQQQ," Ly a is ");
			if( iso_ctrl.ipLyaRedist[ipISO] ==ipCRDW )
			{
				fprintf(ioQQQ,"complete redistribution with wings\n");
			}
			else if( iso_ctrl.ipLyaRedist[ipISO] ==ipCRD )
			{
				fprintf(ioQQQ,"complete redistribution with core only.\n");
			}
			else if( iso_ctrl.ipLyaRedist[ipISO] ==ipPRD )
			{
				fprintf(ioQQQ,"partial redistribution.\n");
			}
			else if( iso_ctrl.ipLyaRedist[ipISO] ==ipLY_A )
			{
				fprintf(ioQQQ,"special Lya.\n");
			}
			else
			{
				fprintf(ioQQQ," PROBLEM Impossible value for iso_ctrl.ipLyaRedist.\n");
				TotalInsanity();
			}

			fprintf(ioQQQ," Other %s resonance lines are ",
				elementnames.chElementSym[ipISO] );

			if( iso_ctrl.ipResoRedist[ipISO] ==ipCRDW )
			{
				fprintf(ioQQQ,"complete redistribution with wings\n");
			}
			else if( iso_ctrl.ipResoRedist[ipISO] ==ipCRD )
			{
				fprintf(ioQQQ,"complete redistribution with core only.\n");
			}
			else if( iso_ctrl.ipResoRedist[ipISO] ==ipPRD )
			{
				fprintf(ioQQQ,"partial redistribution.\n");
			}
			else
			{
				fprintf(ioQQQ," PROBLEM Impossible value for iso_ctrl.ipResoRedist.\n");
				TotalInsanity();
			}

			fprintf(ioQQQ," %s subordinate lines are ",
				elementnames.chElementSym[ipISO] );

			if( iso_ctrl.ipSubRedist[ipISO] ==ipCRDW )
			{
				fprintf(ioQQQ,"complete redistribution with wings\n");
			}
			else if( iso_ctrl.ipSubRedist[ipISO] ==ipCRD )
			{
				fprintf(ioQQQ,"complete redistribution with core only.\n");
			}
			else if( iso_ctrl.ipSubRedist[ipISO] ==ipPRD )
			{
				fprintf(ioQQQ,"partial redistribution.\n");
			}
			else
			{
				fprintf(ioQQQ," PROBLEM Impossible value for iso_ctrl.ipSubRedist.\n");
				TotalInsanity();
			}
		}
		else
		{
			fprintf(ioQQQ," here should have been another keyword on this command.\n");
			fprintf(ioQQQ," Options are ALPHA, RESONANCE, SUBORDINATE.  Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( ipISO == ipH_LIKE && p.nMatch("KEEP") && p.nMatch("FINE") && p.nMatch("STRU") )
	{
		// option to save fine structure components on the line stack
		// this option is not valid for the He-like sequence
		if( p.nMatch(" OFF") )
			iso_ctrl.lgKeepFS = false;
		else
			iso_ctrl.lgKeepFS = true;
	}

	else if( p.nMatch("TOPO") )
	{
		if( p.nMatch(" OFF") )
			iso_ctrl.lgTopoff[ipISO] = false;
		else
			iso_ctrl.lgTopoff[ipISO] = true;
	}

	else
	{
		fprintf( ioQQQ, " There should have been a keyword on this DATABASE H-LIKE or HE-LIKE command.\n Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
