/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseSet scan parameters off SET command */
#include "cddefines.h"
#include "input.h"
#include "prt.h"
#include "rt.h"
#include "phycon.h"
#include "optimize.h"
#include "hcmap.h"
#include "hmi.h"
#include "iterations.h"
#include "conv.h"
#include "secondaries.h"
#include "rfield.h"
#include "ionbal.h"
#include "numderiv.h"
#include "dynamics.h"
#include "iso.h"
#include "predcont.h"
#include "save.h"
#include "stopcalc.h"
#include "opacity.h"
#include "peimbt.h"
#include "radius.h"
#include "atmdat.h"
#include "continuum.h"
#include "grains.h"
#include "grainvar.h"
#include "parser.h"
#include "lines.h"
#include "monitor_results.h"
#include "mole.h"
#include "dense.h"
#include "hyperfine.h"
#include "service.h"

void ParseSet(Parser &p)
{
	long int ip;
	string chString_quotes_lowercase;
	DEBUG_ENTRY( "ParseSet()" );

	/* first check for any strings within quotes - this will get the string
	 * and blank it out, so that these are not confused with keywords.  if
	 * double quotes not present routine returns unity, zero if found*/
	bool lgQuotesFound = true;
	if (p.GetQuote(chString_quotes_lowercase))
		lgQuotesFound = false;

	/* commands to set certain variables, "SET XXX=" */

	if( p.nMatch("LYA") && p.nMatch("21CM") )
	{
		/* set the shape of the Lya source function at line center.  Important for
		 * 21 cm pumping and the resulting spin temperature
		 */
		if( p.nMatch("EXCI") )
		{
			/* the default, the 2p / 1s excitation temperature */
			hyperfine.LyaSourceFunctionShape_assumed = t_hyperfine::EXCITATION;
		}
		else if( p.nMatch("KINE") )
		{
			/* the gas kinetic temperature */
			hyperfine.LyaSourceFunctionShape_assumed = t_hyperfine::KINETIC;
		}
		else if( p.nMatch("CONS") )
		{
			/* S_nu = constant, as assumed in stellar atmosphere calculations */
			hyperfine.LyaSourceFunctionShape_assumed = t_hyperfine::CONSTANT;
		}
	}

	else if (p.nMatch("MONI") && p.nMatch("SCIE"))
	{
		/* print monitored values using scientific notation.  Useful when an asserted value is 
		 * less than 10^-4 of the normalization line. */
		lgPrtSciNot = true;
	}

	else if( p.nMatch("BLEN") )
	{
		/* option for the user to define a custom blend */
		t_blend blnd;
		if( lgQuotesFound )
		{
			blnd.chLabel = chString_quotes_lowercase.substr(0,NCHLAB-1);
			if( chString_quotes_lowercase.length() > NCHLAB-1 )
				fprintf( ioQQQ, " WARNING blend label is too long, truncated to: \"%s\"\n", blnd.chLabel.c_str() );
		}
		blnd.wave = p.FFmtRead();
		if( p.lgEOL() )
			blnd.wave = -1_r;
		else
			if( blnd.wave < 0_r )
			{
				fprintf( ioQQQ, " PROBLEM invalid wavelength supplied: %g\n", blnd.wave );
				cdEXIT(EXIT_FAILURE);
			}
		/* option to quietly ignore blend if database lines are not loaded for any atomic species in the blend */
		blnd.lgQuiet = p.nMatch("QUIE");
		ParseLineList(p, blnd.component);
		prt.blend.emplace_back(blnd);
	}

	/* check energy conservation on every zone - relatively slow */
	else if (p.nMatch("CHECK") && p.nMatch("ENERGY") && p.nMatch("EVERY")
			&& p.nMatch("ZONE"))
	{
		continuum.lgCheckEnergyEveryZone = true;
	}

	else if (p.nMatch("UPDA") && p.nMatch("COUP") &&
				p.nMatch("EVER") && p.nMatch("ION") )
	{
		conv.lgUpdateCouplings = true;
	}

	/* enable Dere07 collisional ionization rate coeffs */
	else if (p.nMatch(" COLL") && p.nMatch(" IONIZ"))
	{
		if (p.nMatch(" HYBRID"))
		{
			atmdat.CIRCData = t_atmdat::HYBRID;
		}
		else if (p.nMatch(" DIMA"))
		{
			atmdat.CIRCData = t_atmdat::DIMA;
		}
		else
		{
			fprintf(ioQQQ,
					" \nPROBLEM Unrecognized set collisional ionization data option.\n");
			fprintf(ioQQQ, " Valid options are Dima or Hybrid.\n");
			fprintf(ioQQQ, " See Hazy 1 for details.\n");
			cdEXIT( EXIT_FAILURE );
		}
	}

	else if( p.nMatch(" H2 " ) && p.nMatch( "CONT" ) && p.nMatch( "DISS" ))
	{
		if( p.nMatch( "STAN" ) )
		{
			/* This uses on the H2 direct photodissociation cross sections
			 * calculated by Philip Stancil */
			mole_global.lgStancil = true;
		}

		else if( p.nMatch( "AD69" ) )
		{
			 /* Use constant cross section from
			    >>refer	H2	dissoc	Allison, A.C. & Dalgarno, A. 1969,   
			    Atomic Data, 1, 91 */ 

			mole_global.lgStancil = false;
		}
		else
		{
			/* should not have happened ... */
			fprintf( ioQQQ, " There should have been an option on this SET H2 CONTinuum DISSociation command.\n" );
			fprintf( ioQQQ, " consult Hazy to find valid options.\n Sorry.\n" );
			cdEXIT(EXIT_FAILURE);		
		}

	}

	else if( p.nMatch(" CHA") && !p.nMatch( "HO ") && !p.nMatch( "HHE") )
	{
		/* set log of minimum charge transfer rate for high ions and H
		 * default of 1.92e-10 set in zero */
		atmdat.HCTAlex = p.FFmtRead();
		if (p.lgEOL())
		{
			p.NoNumb("minimum charge transfer rate");
		}
		if (atmdat.HCTAlex < 0.)
		{
			atmdat.HCTAlex = exp10(atmdat.HCTAlex);
		}
	}
	else if( p.nMatch("CHEM") && !p.nMatch( "HO ") && !p.nMatch(" HHE") )
	{
		/* turn on Steve Federman's chemistry */
		if (p.nMatch("FEDE"))
		{
			if (p.nMatch(" ON "))
			{
				/* This turns on the diffuse cloud chemical rates of 
				 * >>refer	CO	chemistry	Zsargo, J. & Federman, S. R. 2003, ApJ, 589, 319*/
				mole_global.lgFederman = true;
			}
			else if( p.nMatch( " OFF" ) )
			{
				mole_global.lgFederman = false;
			}
			else
			{
				/* this is the default when command used - true */
				mole_global.lgFederman = true;
			}
		}
		/* >>chng 06 may 30 --NPA.  Turn on non-equilibrium chemistry */
		else if (p.nMatch(" NON") && p.nMatch("EQUI"))
		{

			/* option to use effective temperature as defined in
			 * >>refer	CO	chemistry	Zsargo, J. & Federman, S. R. 2003, ApJ, 589, 319
			 * By default, this is false - changed with set chemistry command */

			mole_global.lgNonEquilChem = true;

			/* >>chng 06 jul 21 -- NPA.  Option to include non-equilibrium
			 * effects for neutral reactions with a temperature dependent rate.
			 * Reasoning is that non-equilibrium chemistry is caused by MHD, and
			 * if magnetic field is only coupled to ions, then neutrals may not be 
			 * affected.  
			 * >>refer	CO	chemistry	Zsargo, J. & Federman, S. R. 2003, ApJ, 589, 319
			 * By default, this is true - changed with set chemistry command */

			if (p.nMatch("NEUT"))
			{
				if (p.nMatch(" ON "))
				{
					/* This turns on the diffuse cloud chemical rates of 
					 * >>refer	CO	chemistry	Zsargo, J. & Federman, S. R. 2003, ApJ, 589, 319*/
					mole_global.lgNeutrals = true;
				}
				else if( p.nMatch( " OFF" ) )
				{
					mole_global.lgNeutrals = false;
				}
				else
				{
					/* this is the default when command used - true */
					mole_global.lgNeutrals = true;
				}
			}
		}

		/* turn off proton elimination rates, which are of the form
		 *
		 *
		 *    A + BH+ -->  AB + H+ or
		 *    AH + B+ -->  AB + H+ 
		 *
		 * the following paper:
		 *
		 * >>refer	CO	chemistry	Huntress, W. T., 1977, ApJS, 33, 495
		 * says reactions of these types are much less likely than
		 * identical reactions which leave the product AB ionized (AB+), 
		 * leaving an H instead of H+ (this is called H atom elimination 
		 * currently we only have one reaction of this type, it is
		 * C+ + OH -> CO + H+ */
		else if (p.nMatch("PROT") && p.nMatch("ELIM"))
		{
			if (p.nMatch(" ON "))
			{
				/* This turns on the diffuse cloud chemical rates of 
				 * >>refer	CO	chemistry	Zsargo, J. & Federman, S. R. 2003, ApJ, 589, 319*/
				mole_global.lgProtElim = true;
			}
			else if( p.nMatch( " OFF" ) )
			{
				mole_global.lgProtElim = false;
			}
			else
			{
				/* this is the default when command used - true */
				mole_global.lgProtElim = true;
			}
		}

		else
		{
			/* should not have happened ... */
			fprintf(ioQQQ,
					" There should have been an option on this SET CHEMISTRY command.\n");
			fprintf(ioQQQ, " consult Hazy to find valid options.\n Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* set collision strength averaging ON / OFF */
	else if( p.nMatch("COLL") && p.nMatch("STRE") && p.nMatch("AVER") )
	{
		if( p.nMatch(" OFF") )
		{
			 iso_ctrl.lgCollStrenThermAver = false;
		}
		else
		{
			/* this is default behavior of this command */
			 iso_ctrl.lgCollStrenThermAver = true;
		}
	}

	else if (p.nMatch("COVE"))
	{
		iterations.lgConverge_set = true;
		/* set coverage - limit number of iterations and zones */
		if (p.nMatch("FAST"))
		{
			iterations.lim_zone = 1;
			iterations.lim_iter = 0;
		}
		else
		{
			iterations.lim_zone = 10;
			iterations.lim_iter = 1;
		}
	}

	else if (p.nMatch("CSUP"))
	{
		/* force H^0 secondary ionization rate, log entered */
		secondaries.SetCsupra = (realnum) p.FFmtRead();
		secondaries.lgCSetOn = true;
		if (p.lgEOL())
		{
			p.NoNumb("secondary ionization rate");
		}
		secondaries.SetCsupra = exp10(secondaries.SetCsupra);
	}

	else if (p.nMatch("CUMU"))
	{
		if (p.nMatch("OFF"))
		{
			strcpy(rfield.chCumuType,"NONE");
		} 
		else if (p.nMatch("MASS"))
		{
			strcpy(rfield.chCumuType,"MASS");
		}
		else if (p.nMatch("FLUX"))
		{
			strcpy(rfield.chCumuType,"FLUX");
		}
		else
		{
			fprintf(ioQQQ,
					  " Did not recognize option on 'set cumulative' command."
					  " It should be FLUX, MASS or OFF.  Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch(" D/H"))
	{
		/* set deuterium abundance, D to H ratio */
		double tmp = p.FFmtRead();
		if (p.lgEOL())
		{
			fprintf(ioQQQ,
				"The command 'set D/H' has been deprecated.\n"
				"Please use 'element hydrogen isotopes' instead.\n");
		}
		else
		{
			if (tmp <= 0. || p.nMatch(" LOG"))
			{
				tmp = exp10(tmp);
			}
			fprintf(ioQQQ,
				"The command 'set D/H' has been deprecated.\n"
				"Please use 'element hydrogen isotopes (1, 1) (2, %g)' instead.\n",
				tmp);
		}
		cdEXIT( EXIT_FAILURE );
	}

	else if( p.nMatch("DENS") && p.nMatch("TOLE") )
	{
		/* set error in total gas-phase density of each element, including molecules */
		conv.GasPhaseAbundErrorAllowed = (realnum)p.FFmtRead();
		if( p.lgEOL() )
		{
			p.NoNumb("density tolerance");
		}
		if( conv.GasPhaseAbundErrorAllowed <= 0. )
		{
			conv.GasPhaseAbundErrorAllowed = exp10(conv.GasPhaseAbundErrorAllowed);
		}
	}

	else if( p.nMatch( " HO " ) && p.nMatch( "CHAR" ) )
	{
		fprintf(ioQQQ, " The SET HO CHAR command is no longer supported.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( p.nMatch( " HHE" ) && p.nMatch( "CHAR" ) )
	{
		fprintf(ioQQQ, " The SET HHE CHAR command is no longer supported.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( p.nMatch("12C1") )
	{
		/* set the 12C to 13C abundance ratio - default is 30 */
		/* first two numbers on the line are 12 and 13 - we don't want them */
		(void) p.FFmtRead();
		(void) p.FFmtRead();

		/* now we can get the ratio */
		realnum c12c13 = (realnum) p.FFmtRead();
		if (p.lgEOL())
		{
			fprintf(ioQQQ, "Deprecated option.  Please use command: \"element carbon isotope\"\n");
		}
		else
		{
			if (c12c13 <= 0. || p.nMatch(" LOG"))
				c12c13 = exp10(c12c13);
			fprintf(ioQQQ, "Deprecated option.  Please use command: \"element carbon isotope (12, %.2f), (13, 1)\"\n",
					c12c13);
		}
		cdEXIT(EXIT_FAILURE);
	}

	/* set dynamics ... */
	else if (p.nMatch("DYNA"))
	{
		/* set dynamics advection length */
		if (p.nMatch("ADVE") && p.nMatch("LENG"))
		{
			/* <0 => relative fraction of length, + val in cm */
			dynamics.AdvecLengthInit = p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("advection length");
			/* if fraction is present, then number was linear fraction, if not present
			 * then physical length in cm, log10 */
			if (p.nMatch("FRAC"))
			{
				/* we scanned in the number, if it is a negative then it is the log of the fraction */
				if (dynamics.AdvecLengthInit <= 0.)
					dynamics.AdvecLengthInit = exp10(dynamics.AdvecLengthInit);

				/* make neg sign as flag for this in dynamics.c */
				dynamics.AdvecLengthInit *= -1.;
			}
			else
			{
				/* fraction did not occur, the number is the log of the length in cm -
				 * convert to linear cm */
				dynamics.AdvecLengthInit = exp10(dynamics.AdvecLengthInit);
			}
		}
		else if (p.nMatch("PRES") && p.nMatch("MODE"))
		{
			dynamics.lgSetPresMode = true;
			if (p.nMatch("SUBS"))
			{
				/* subsonic */
				strcpy(dynamics.chPresMode, "subsonic");
			}
			else if (p.nMatch("SUPE"))
			{
				/* supersonic */
				strcpy(dynamics.chPresMode, "supersonic");
			}
			else if (p.nMatch("STRO"))
			{
				/* strong d */
				strcpy(dynamics.chPresMode, "strongd");
			}
			else if (p.nMatch("ORIG"))
			{
				/* original */
				strcpy(dynamics.chPresMode, "original");
			}
		}
		else if (p.nMatch("ANTI") && p.nMatch("DEPT"))
		{
			dynamics.lgSetPresMode = true;
			strcpy(dynamics.chPresMode, "antishock");
			/* shock depth */
			/* get log of shock depth in cm */
			dynamics.ShockDepth = p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("antishock depth");
			dynamics.ShockDepth = exp10(dynamics.ShockDepth);
		}
		else if (p.nMatch("ANTI") && p.nMatch("MACH"))
		{
			dynamics.lgSetPresMode = true;
			strcpy(dynamics.chPresMode, "antishock-by-mach");
			/* Mach number */
			/* get (isothermal) Mach number where we want antishock to occur */
			dynamics.ShockMach = p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("antishock-by-mach");
		}
		else if (p.nMatch("RELA"))
		{
			/* set how many iterations we will start with, before allowing
			 * changes.  This allows the solution to relax to an equilibrium */
			dynamics.n_initial_relax = (long int) p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("relaxation cycles before start of dynamics");
			else if (dynamics.n_initial_relax < 2)
			{
				fprintf(ioQQQ,
						" First iteration to relax dynamics must be > 1."
							"It was %li. Sorry.\n", dynamics.n_initial_relax);
				cdEXIT(EXIT_FAILURE);
			}
		}
		else if (p.nMatch("SHOC") && p.nMatch("DEPT"))
		{
			dynamics.lgSetPresMode = true;
			strcpy(dynamics.chPresMode, "shock");
			/* shock depth */
			/* get log of shock depth in cm */
			dynamics.ShockDepth = p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("shock depth");
			dynamics.ShockDepth = exp10(dynamics.ShockDepth);
		}
		else if (p.nMatch("POPU") && p.nMatch("EQUI"))
		{
			// set dynamics populations eqauilibrium
			// force equilibrium populations 
			dynamics.lgEquilibrium = true;
		}
		else
		{
			/* should not have happened ... */
			fprintf(ioQQQ,
					" There should have been an option on this SET DYNAMICS command.\n");
			fprintf(ioQQQ, " consult Hazy to find valid options.\n Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch("DIDZ"))
	{
		/* set parameter to do with choice of dr;
		 * par is the largest optical depth to allow in the zone. 
		 * >>chng 96 jan 08 had been two numbers - dtau1 removed */
		radius.drChange = (realnum) p.FFmtRead();
		if (radius.drChange <= 0.)
		{
			radius.drChange = exp10(radius.drChange);
		}
		if (p.lgEOL())
		{
			p.NoNumb("largest optical depth allowed in zone");
		}
	}

	/* something to do with electron density */
	else if (p.nMatch("EDEN"))
	{
		/* set eden convergence sets convergence criterion, keyword parallel with
		 * set temperature convergence,  but also accept (older) error as keyword */
		if (p.nMatch("CONV") || p.nMatch("ERRO"))
		{
			/* keyword is eden convergence  
			 * set tolerance in eden match */
			conv.EdenErrorAllowed = p.FFmtRead();
			if (p.lgEOL())
			{
				p.NoNumb("electron density error allowed");
			}

			if (conv.EdenErrorAllowed < 0.)
			{
				conv.EdenErrorAllowed = exp10(conv.EdenErrorAllowed);
			}
		}
		else if(p.nMatch("FRACTION"))
		{
			/* set eden fraction sets the ratio eden/hden */
			dense.EdenFraction = p.FFmtRead();
			if (p.lgEOL())
			{
				p.NoNumb("electron density fraction");
			}

			if (dense.EdenFraction <= 0. )
			{
				dense.EdenFraction = exp10(dense.EdenFraction);
			}
			/* warn that this model is meaningless */
			phycon.lgPhysOK = false;
		}
		else if (p.nMatch("SOLV"))
		{
			/* options are vWDB (default) and SECAnt */
			if (p.nMatch("VWDB"))
			{
				strcpy( conv.chSolverEden, "vWDB" );
			}
			else if (p.nMatch("SECA"))
			{
				strcpy( conv.chSolverEden, "SECA" );
			}
			else
			{
				fprintf(ioQQQ,
						  "'set eden solver' options are 'vWDB' and 'SECA'.\n");
				cdEXIT( EXIT_FAILURE );
			}
		}
		else
		{
			/* no keyword, assume log of electron density */
			/* set the electron density */
			dense.EdenSet = (realnum) exp10(p.FFmtRead());
			if (p.lgEOL())
			{
				p.NoNumb("electron density");
			}

			/* warn that this model is meaningless */
			phycon.lgPhysOK = false;
		}
	}

	/* Stop using gbar to fill in dBase transitions	 */
	else if( p.nMatch("GBAR"))
	{
		if( p.nMatch(" OFF ") )
		{
			atmdat.lgGbarOn = false;
		}
	}

	/* something to do with electron density */
	else if (p.nMatch("IONI") && p.nMatch("TOLE"))
	{
		/* keyword is ionization tolerance */
		conv.IonizErrorAllowed = p.FFmtRead();
		if (p.lgEOL())
		{
			p.NoNumb("ionization error allowed");
		}
		
		if (conv.IonizErrorAllowed < 0.)
		{
			conv.IonizErrorAllowed = exp10(conv.IonizErrorAllowed);
		}
	}

	else if (p.nMatch("ISOTOPE") || p.nMatch("ISOTOPO") )
	{
		if( p.nMatch(" ALL ") )
		{
			fprintf(ioQQQ,
				"The command 'set isotopes all' has been deprecated.\n"
				"Please use 'element isotopes all' instead.\n");
			cdEXIT( EXIT_FAILURE );
		}
	}


	else if (p.nMatch("FINE") && p.nMatch("CONT"))
	{
		/* set fine continuum resolution - an element name, used to get
		 * thermal width, and how many resolution elements to use to resolve 
		 * a line of this element at 1e4 K 
		 * first get an element name, nelem is atomic number on C scale
		 * default is iron */
		if ((rfield.fine_opac_nelem = p.GetElem()) < 0)
		{
			fprintf(ioQQQ,
					" An element name must appear on this line\n Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
		/* set the number of resolution elements in HWHM at 1e4 K for turbulent 
		 * velocity field, default is one element */
		rfield.fine_opac_nresolv = (long int) p.FFmtRead();
		if (rfield.fine_opac_nresolv < 1)
		{
			fprintf(ioQQQ,
					" The number of resolution elements within FWHM of line must appear\n Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* option to specify the lower and upper limits of the fine continuum
		 * the upper limit is always optional.  */
		if (!p.lgEOL())
		{
			realnum lower_limit = (realnum) p.FFmtRead();
			if (lower_limit < 0)
				lower_limit = exp10(lower_limit);

			if (lower_limit > 0.2f)
				fprintf(
						ioQQQ,
						" The fine continuum lower limit is quite high (%f Ryd). Please check.\n",
						lower_limit);

			/* reset to coarse limits if arguments are beyond them */
			rfield.fine_ener_lo = MAX2(rfield.emm(), lower_limit);

			if (!p.lgEOL())
			{
				realnum upper_limit = (realnum) p.FFmtRead();
				if (upper_limit < 0)
					upper_limit = exp10(upper_limit);

				if (upper_limit < 10.f)
					fprintf(
							ioQQQ,
							" The fine continuum upper limit is quite low (%f Ryd). Please check.\n",
							upper_limit);

				/* reset to coarse limits if arguments are beyond them */
				rfield.fine_ener_hi = MIN2(rfield.egamry(), upper_limit);
			}
		}
	}

	/* set grain command - but not set H2 grain command */
	else if (p.nMatch("GRAI") && !p.nMatch(" H2 "))
	{
		if (p.nMatch("HEAT"))
		{
			/* scale factor to change grain heating as per Allers et al. */
			gv.GrainHeatScaleFactor = (realnum) p.FFmtRead();
			/* warn that this model is not the best we can do */
			phycon.lgPhysOK = false;
			if (p.lgEOL())
			{
				p.NoNumb("grain heating");
			}
		}
		else
		{
			fprintf(ioQQQ,
					" A keyword must appear on the SET GRAIN line - options are HEAT \n Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* these are the "set Leiden hack" commands, used to turn off physics and
	 * sometimes replace with simple approximation */
	else if (p.nMatch("LEID") && p.nMatch("HACK"))
	{
		if (p.nMatch("H2* ") && p.nMatch(" OFF"))
		{
			/* turn off reactions with H2* in the chemistry network */
			hmi.lgLeiden_Keep_ipMH2s = false;
			/* warn that this model is not the best we can do */
			phycon.lgPhysOK = false;
		}
		else if (p.nMatch("CR ") && p.nMatch(" OFF"))
		{
			/* the CR Leiden hack option - turn off CR excitations of H2 */
			hmi.lgLeidenCRHack = false;
		}
		else if( p.nMatch("RATE") &&  p.nMatch("UMIS"))
		{
			/* This command will use the rates given in the UMIST database,
			 * It will set to zero many reactions that are not in UMIST */
			mole_global.lgLeidenHack = true;
		}
	}

	/* set H2 ... */
	else if (p.nMatch(" H2 "))
	{
		ip = (long int) p.FFmtRead();
		if (ip != 2)
		{
			fprintf(ioQQQ,
					" The first number on this line must be the 2 in H2\n Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* SET H2 SMALL MODEL or SOLOMON - which approximation to Solomon process */
		if (p.nMatch("SOLO") || (p.nMatch("SMAL") && p.nMatch("MODE")))
		{
			if (p.nMatch("SOLO"))
			{
				/* warn that Solomon will not be used forever */
				fprintf(ioQQQ,
						"PROBLEM - *set H2 Solomon* has been changed to *set H2 small model*."
							"  This is OK for now but it may not work in a future version.\n");
			}
			if (p.nMatch("TH85"))
			{
				/* rate from is eqn A8 of Tielens and Hollenbach 85a, */
				hmi.chH2_small_model_type = 'T';
			}
			else if (p.nMatch(" BHT"))
			{
				/* the improved H2 formalism given by 
				 *>>refer	H2	dissoc	Burton, M.G., Hollenbach, D.J. & Tielens, A.G.G.M
				 >>refcon	1990, ApJ, 365, 620 */
				hmi.chH2_small_model_type = 'H';
			}
			else if (p.nMatch("BD96"))
			{
				/* this rate is equation 23 of
				 *>>refer	H2	dissoc	Bertoldi, F., & Draine, B.T., 1996, 458, 222 */
				/* this is the default */
				hmi.chH2_small_model_type = 'B';
			}
			else if (p.nMatch("ELWE"))
			{
				/* this rate is equation 23 of
				 *>>refer	H2	dissoc	Elwert et al., in preparation */
				/* this is the default */
				hmi.chH2_small_model_type = 'E';
			}
			else
			{
				fprintf(ioQQQ,
						" One of the keywords TH85, _BHT, BD96 or ELWErt must appear.\n Sorry.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* series of commands that deal with grains */
		/* which approximation to grain formation pumping */
		if (p.nMatch("GRAI") && p.nMatch("FORM") && p.nMatch("PUMP"))
		{
			if (p.nMatch("DB96"))
			{
				/* Draine & Bertoldi 1996 */
				hmi.chGrainFormPump = 'D';
			}
			else if (p.nMatch("TAKA"))
			{
				/* the default from
				 * >>refer	H2	pump	Takahashi, Junko, 2001, ApJ, 561, 254-263 */
				hmi.chGrainFormPump = 'T';
			}
			else if (p.nMatch("THER"))
			{
				/* thermal distribution, upper right column of page 239 of
				 *>>refer	H2	formation	Le Bourlot, J, 1991, A&A, 242, 235 */
				hmi.chGrainFormPump = 't';
			}
			else if (p.nMatch(" OFF"))
			{
				/* disable grain formation pumping */
				hmi.chGrainFormPump = ' ';
			}
			else
			{
				fprintf(ioQQQ,
						" The grain form pump option is wrong.\n Sorry.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* which approximation to Jura rate */
		else if (p.nMatch("JURA"))
		{
			if (p.nMatch("TH85"))
			{
				/* rate from is eqn A8 of Tielens and Hollenbach 85a*/
				hmi.chJura = 'T';
			}
			else if (p.nMatch("CT02"))
			{
				/* this rate is equation Cazeux & Tielens */
				hmi.chJura = 'C';
			}
			else if (p.nMatch("ELRD"))
			{
				/* this rate is equation C.1 from Rollig et al., and includes Eley-Rideal Effect*/
				hmi.chJura = 'E';
			}
			else if (p.nMatch("SN99"))
			{
				/* this rate is from Sternberg & Neufeld 99 */
				hmi.chJura = 'S';
			}
			else if (p.nMatch("RATE"))
			{
				/* set H2 rate - enters log of Jura rate - F for fixed
				 * no dependence on grain properties 
				 * had been C, a bug since triggered Cazeux & Tielens
				 * >>chng 07 jan 10, bug caught by Robin Williams & Fixed by Nick Abel */
				hmi.chJura = 'F';
				hmi.rate_h2_form_grains_set = exp10(p.FFmtRead());
				if (p.lgEOL())
				{
					/* no number on the line so use Jura's value, 3e-17 
					 * >>refer	H2	Jura	Jura, M. 1975, ApJ, 197, 575 */
					hmi.rate_h2_form_grains_set = 3e-17;
				}
			}
			else if (p.nMatch("SCAL"))
			{
				/* this is a scale factor to multiply the Jura rate */
				hmi.ScaleJura = (realnum) p.FFmtRead();
				/* log or negative number means log was entered */
				if (p.nMatch(" LOG") || hmi.ScaleJura < 0.)
				{
					hmi.ScaleJura = exp10(hmi.ScaleJura);
				}
				if (p.lgEOL())
					p.NoNumb("scale for Jura rate");

				/* option to vary scale factor */
				if (optimize.lgVarOn)
				{
					optimize.nvarxt[optimize.nparm] = 1;
					strcpy(optimize.chVarFmt[optimize.nparm],
							"SET H2 JURA SCALE %f LOG");

					/* pointer to where to write */
					optimize.nvfpnt[optimize.nparm] = input.nRead;

					/* log of Jura rate scale factor will be parameter */
					optimize.vparm[0][optimize.nparm] = (realnum) log10(
							hmi.ScaleJura);
					optimize.vincr[optimize.nparm] = 0.3f;

					++optimize.nparm;
				}
			}
			else
			{
				fprintf(ioQQQ, " The Jura rate option is wrong.\n Sorry.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		/* what temperature to use for binding energy, Tad in Le Bourlot, J., 2000, A&A, 360, 656-662  */
		else if (p.nMatch(" TAD"))
		{
			hmi.Tad = (realnum) p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("temperature for binding energy");
			/* log if <= 10. unless linear key appears too */
			if (hmi.Tad <= 10. && !p.nMatch("LINE"))
				hmi.Tad = exp10(hmi.Tad);
		}

		else if (p.nMatch("FRAC"))
		{
			/* this is special option to force H2 abundance to value for testing
			 * this factor will multiply the hydrogen density to become the H2 density
			 * no attempt to conserve particles, or do the rest of the molecular equilibrium
			 * set consistently is made */
			hmi.H2_frac_abund_set = p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("H2 fractional abundance");

			/* a number <= 0 is the log of the ratio */
			if (hmi.H2_frac_abund_set <= 0.)
				hmi.H2_frac_abund_set = exp10(hmi.H2_frac_abund_set);
			/* don't let it exceed 0.5 */
			/* >>chng 03 jul 19, from 0.5 to 0.4999, do not want atomic density exactly zero */
			hmi.H2_frac_abund_set = MIN2(0.49999, hmi.H2_frac_abund_set);
		}
#if 0
		else if( p.nMatch("FORM") && p.nMatch("SCAL") )
		{
			/* this is special option to scale H2 formation rate. 
			 * In the fully molecular or fully ionized limits,
			 * this should supercede the above "FRAC" option because
			 * it allows the same thing without breaking the chemistry
			 * or any conservation checks */
			hmi.H2_formation_scale = p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("H2 formation scale");

			/* a number <= 0 is the log of the ratio */
			if (hmi.H2_formation_scale <= 0.)
				hmi.H2_formation_scale = exp10(hmi.H2_frac_abund_set);
		}
#endif
	}

	/* this is a scale factor that changes the n(H0)*1.7e-4 that is added to the
	 * electron density to account for collisions with atomic H.  it is an order of
	 * magnitude guess, so this command provides ability to see whether it affects results */
	else if (p.nMatch("HCOR"))
	{
		dense.HCorrFac = (realnum) p.FFmtRead();
		if (p.lgEOL())
			p.NoNumb("scale for H0 correction to e- collision rate");
	}

	else if (p.nMatch(" PAH"))
	{
		/* set one of several possible PAH abundance distribution functions */
		if (lgQuotesFound)
		{
			/* abundance depends specified species as fraction of Htot */
			gv.chPAH_abundance = chString_quotes_lowercase;
		}
		else if (p.nMatch("CONS"))
		{
			/* constant PAH abundance */
			gv.chPAH_abundance = "CON";
		}
		else if (p.nMatch("BAKE"))
		{
			/* turn on simple PAH heating from Bakes & Tielens - this is very approximate */
			/*>>>refer	PAH	heating	Bakes, E.L.O., & Tielens, A.G.G.M. 1994, ApJ, 427, 822 */
			gv.lgBakesPAH_heat = true;
			/* warn that this model is not the best we can do: energy conservation is violated */
			phycon.lgPhysOK = false;
		}
		else
		{
			fprintf(ioQQQ,
					" a string, or one of the keywords BAKES, or CONStant must appear.\n Sorry.");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* something to do with pressure */
	else if (p.nMatch("PRES"))
	{
		/* tolerance on pressure convergence */
		if (p.nMatch("CONV"))
		{
			/* keyword is tolerance 
			 * set tolerance or relative error in pressure match */
			conv.PressureErrorAllowed = (realnum) p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("pressure convergence tolerance");

			if (conv.PressureErrorAllowed < 0.)
				conv.PressureErrorAllowed = exp10(conv.PressureErrorAllowed);
		}
		else if( p.nMatch("IONI") )
		{
			/* set limit to number of calls from pressure to ionize solver,
			 * this set limit to conv.nPres2Ioniz */
			conv.limPres2Ioniz = (long) p.FFmtRead();
			if (p.lgEOL())
			{
				p.NoNumb("number of calls from pressure to ion solver");
			}
			else if (conv.limPres2Ioniz <= 0)
			{
				fprintf(ioQQQ, " The limit must be greater than zero.\n Sorry.");
				cdEXIT(EXIT_FAILURE);
			}
		}

		else
		{
			/* >>chng 04 mar 02, printout had been wrong, saying TOLErange
			 * rather than CONVergence.  Nick Abel */
			fprintf(ioQQQ,
					" I didn\'t recognize a key on this SET PRESSURE line.\n");
			fprintf(ioQQQ, " The ones I know about are: CONVergence and IONIze.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}
	else if (p.nMatch("RECOMBIN"))
	{
		/* dielectronic recombination */
		if (p.nMatch("DIELECTR"))
		{
			if (p.nMatch("MEAN"))
			{
				/* change various factors for the dielectronic recombination means */
				if (p.nMatch(" OFF"))
				{
					for (int ion = 0; ion < LIMELM; ++ion)
						ionbal.DR_mean_scale[ion] = 0.;
				}
				/* multiply guess by scale factor */
				else if (p.nMatch("SCALE"))
				{
					// there must be at least one scale factor
					ionbal.DR_mean_scale[0] = (realnum) p.FFmtRead();
					if (p.lgEOL())
					{
						fprintf(
								ioQQQ,
								" There must be at least one scale factor on the SET RECOMBIANTION MEAN command.\n");
						cdEXIT(EXIT_FAILURE);
					}

					for (int ion = 1; ion < LIMELM; ++ion)
					{
						ionbal.DR_mean_scale[ion] = (realnum) p.FFmtRead();
						if (p.lgEOL())
							ionbal.DR_mean_scale[ion]
									= ionbal.DR_mean_scale[ion - 1];
					}
					for (int ion = 0; ion < LIMELM; ++ion)
					{
						if (ionbal.DR_mean_scale[ion] < 0)
						{
							fprintf(ioQQQ,
									" All scale factors on the SET RECOMBIANTION MEAN command must be >=0.\n");
							cdEXIT(EXIT_FAILURE);
						}
					}
				}
				/* option to add log normal noise */
				else if (p.nMatch("NOISE"))
				{
					ionbal.guess_noise = (realnum) p.FFmtRead();
					if (p.lgEOL())
						ionbal.guess_noise = 2.;
				}
				else
				{
					fprintf(ioQQQ, " key OFF or NOISE must appear.\n");
					cdEXIT(EXIT_FAILURE);
				}
			}
			else if( p.nMatch("SUPP") && p.nMatch(" OFF") )
			{
				ionbal.lgDRsup = false;
			}
			else
			{
				fprintf(ioQQQ, " key MEAN or SUPPression must appear.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			fprintf(ioQQQ,
					" key DIELECTRonic must appear on set recombination command.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch(" DR "))
	{
		/* set zone thickness by forcing drmax and drmin */
		/* at this stage linear, but default is log */
		radius.sdrmax = p.FFmtRead();
		if (!p.nMatch("LINE"))
		{
			/* linear was not on command, so default is log */
			radius.sdrmax = exp10(radius.sdrmax);
		}
		if (p.lgEOL())
		{
			p.NoNumb("zone thickness");
		}

		radius.lgSdrmaxRel = p.nMatch("RELA");

		/* NB these being equal are tested in convinittemp to set dr */
		radius.lgFixed = true;	
		radius.sdrmin = radius.sdrmax;
		radius.lgSdrminRel = radius.lgSdrmaxRel;
		if (!radius.lgSdrmaxRel && radius.sdrmax < DEPTH_OFFSET * 1e4)
		{
			fprintf(
					ioQQQ,
					"\n Thicknesses less than about %.0e will NOT give accurate results. If tricking the code\n",
					DEPTH_OFFSET * 1e4);
			fprintf(
					ioQQQ,
					" into computing emissivities instead of intensities, try to instead use a thickness of unity,\n");
			fprintf(
					ioQQQ,
					" and then multiply (divide) the results by the necessary thickness (product of densities).\n\n");
			cdEXIT(EXIT_FAILURE);
		}
		if (radius.lgSdrmaxRel && (radius.sdrmax <= 0. || radius.sdrmax >= 1.))
		{
			fprintf(
					ioQQQ,
					" When using a relative dr, a fraction between 0 and 1 must be entered. Found: %g\n",
					radius.sdrmax);
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch("DRMA"))
	{
		/* set maximum zone thickness" */
		radius.sdrmax = p.FFmtRead();
		if (p.lgEOL())
			p.NoNumb("maximum zone thickness");

		if ((radius.sdrmax < 38. || p.nMatch(" LOG")) && !p.nMatch("LINE"))
			radius.sdrmax = exp10(radius.sdrmax);

		// option to set min relative to current radius
		radius.lgSdrmaxRel = p.nMatch("RELA");
		if (radius.lgSdrmaxRel && (radius.sdrmax <= 0. || radius.sdrmax >= 1.))
		{
			fprintf(
					ioQQQ,
					" When using a relative drmax, a fraction between 0 and 1 must be entered. Found: %g\n",
					radius.sdrmax);
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch("DRMI") && p.nMatch("DEPT"))
	{
		/* option to set minimum zone thickness relative to depth */
		radius.sdrmin_rel_depth = p.FFmtRead();
		if (p.lgEOL())
			p.NoNumb("minimum zone thickness rel to depth");

		if ((radius.sdrmin_rel_depth < 38. || p.nMatch(" LOG")) && !p.nMatch("LINE"))
			radius.sdrmin_rel_depth = exp10(radius.sdrmin_rel_depth);

		if( radius.sdrmin_rel_depth <= 0. || radius.sdrmin_rel_depth >= 1.)
		{
			fprintf( ioQQQ, " When using a relative drmin, a fraction between 0 and 1 must be entered. Found: %g\n",
					radius.sdrmin_rel_depth );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch("DRMI"))
	{
		/* option to set minimum zone thickness */
		radius.sdrmin = p.FFmtRead();
		if (p.lgEOL())
			p.NoNumb("minimum zone thickness");

		if ((radius.sdrmin < 38. || p.nMatch(" LOG")) && !p.nMatch("LINE"))
			radius.sdrmin = exp10(radius.sdrmin);

		// option to set min relative to current radius
		radius.lgSdrminRel = p.nMatch("RELA");
		radius.lgSMinON = true;
		if (radius.lgSdrminRel && (radius.sdrmin <= 0. || radius.sdrmin >= 1.))
		{
			fprintf( ioQQQ, " When using a relative drmin, a fraction between 0 and 1 must be entered. Found: %g\n",
					radius.sdrmin);
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch("FLXF"))
	{
		/* faintest continuum flux to consider */
		rfield.FluxFaint = p.FFmtRead();
		if (p.lgEOL())
		{
			p.NoNumb("faintest continuum flux to consider");
		}
		if (rfield.FluxFaint < 0.)
		{
			rfield.FluxFaint = exp10(rfield.FluxFaint);
		}
	}

	else if (p.nMatch("LINE") && p.nMatch("PREC"))
	{
		fprintf(ioQQQ,"Sorry, this command changed to PRINT LINE PRECISION.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* >>chng 00 dec 08, command added by Peter van Hoof */
	else if (p.nMatch("NFNU"))
	{
		if (p.nMatch(" ADD"))
		{
			/* option to add a continuum point */
			double energy = p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("energy");
			t_PredCont::Inst().add(energy, p.StandardEnergyUnit());
		}
		else
		{
			/* set nFnu [incident_reflected] [incident_transmitted] [diffuse_inward] [diffuse_outward]
			 *   command for specifying what to include in the nFnu entries in LineSave.lines */
			/* >>chng 01 nov 12, also accept form with space */
			/* "incident reflected" keyword */
			prt.lgSourceReflected = p.nMatch("INCIDENT R") || p.nMatch(
					"INCIDENT_R");
			/* "incident transmitted" keyword */
			prt.lgSourceTransmitted = p.nMatch("INCIDENT_T") || p.nMatch(
					"INCIDENT T");
			/* "diffuse inward" keyword */
			prt.lgDiffuseInward = p.nMatch("DIFFUSE_I")
					|| p.nMatch("DIFFUSE I");
			/* "diffuse outward" keyword */
			prt.lgDiffuseOutward = p.nMatch("DIFFUSE_O") || p.nMatch(
					"DIFFUSE O");

			/* total includes all four */
			if( p.nMatch("TOTAL") )
			{
				prt.lgSourceReflected = true;
				prt.lgSourceTransmitted = true;
				prt.lgDiffuseInward =  true;
				prt.lgDiffuseOutward =  true;

			}

			/* at least one of these needs to be set ! */
			if (!(prt.lgSourceReflected || prt.lgSourceTransmitted
					|| prt.lgDiffuseInward || prt.lgDiffuseOutward))
			{
				fprintf(ioQQQ,
						" set nFnu expects one or more of the following keywords:\n");
				fprintf(ioQQQ, " INCIDENT_REFLECTED, INCIDENT_TRANSMITTED, "
					"DIFFUSE_INWARD, DIFFUSE_OUTWARD\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	else if (p.nMatch("IND2"))
	{
		if (p.nMatch(" ON "))
		{
			/* set flag saying to off or on induced two photon processes */
			iso_ctrl.lgInd2nu_On = true;
		}
		else if( p.nMatch(" OFF") )
		{
			iso_ctrl.lgInd2nu_On = false;
		}
		else
		{
			fprintf( ioQQQ, " set ind2 needs either ON or OFF.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch("SPECIES"))
	{
		/* Set the default collision strength for dBase transitions
		 * without collision or radiative data */
		if (p.nMatch(" GBAR"))
		{
			atmdat.collstrDefault = (long)p.FFmtRead();
			if (p.lgEOL())
			{
				p.NoNumb("the default collision strengths when no collision or radiative data are available");
			}
		}

		/* Set the pseudo-continuum range and resolution for given species. */
		else if( p.nMatch("CONT") )
		{
			adjPseudoCont thisCont;

			if( lgQuotesFound )
			{
				string species = chString_quotes_lowercase;
				trimTrailingWhiteSpace( species );
				thisCont.speciesLabel = species;
			}
			else
			{
				fprintf( ioQQQ, "PROBLEM No species specified with the "
						"'set species continuum command'\n" );
				cdEXIT( EXIT_FAILURE );
			}

			/* set lower, upper wavelength range, and number of bins, for
			 * save feii continuum command */
			thisCont.wlLo = (realnum) p.FFmtRead();
			thisCont.wlHi = (realnum) p.FFmtRead();
			thisCont.nBins = (long) p.FFmtRead();

			if (p.lgEOL())
			{
				fprintf(ioQQQ,
						"PROBLEM The 'set species continuum' command must have"
							" three numbers, the lower and upper wavelength range in Angstroms"
							" and the number of bins to divide this into.\n");
				cdEXIT(EXIT_FAILURE);
			}

			if( thisCont.wlLo >= thisCont.wlHi )
			{
				fprintf(ioQQQ,
					"PROBLEM The first two numbers on the 'set "
					"species continuum' command must be the lower and upper "
					"wavelength range in Angstroms and the first must be less "
					"than the second.\n");
				cdEXIT(EXIT_FAILURE);
			}
			if( thisCont.nBins < 2 )
			{
				fprintf(ioQQQ,
					"PROBLEM The third number on the 'set species continuum' "
					"command must be the number of bins to divide the range into and"
					" it must be greater than 1.\n");
				cdEXIT(EXIT_FAILURE);
			}

			save.setPseudoCont.push_back( thisCont );
			{
				enum { DEBUG_SPEC = false };
				if( DEBUG_SPEC )
				{
					fprintf( ioQQQ, "species :\t %s\n", thisCont.speciesLabel.c_str() );
					fprintf( ioQQQ, "conwlLo :\t %g\n", thisCont.wlLo );
					fprintf( ioQQQ, "conwlHi :\t %g\n", thisCont.wlHi );
					fprintf( ioQQQ, "ncon    :\t %ld\n", thisCont.nBins );
				}
			}
		}

		else
		{
			fprintf(ioQQQ, " SET SPECIES takes options GBAR and CONTINUUM.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch("TEMP"))
	{
		/* set something to do with temperature, currently solver, tolerance, floor */
		if (p.nMatch("FLOO"))
		{
			StopCalc.TeFloor = (realnum) p.FFmtRead();
			if (p.lgEOL())
				p.NoNumb("temperature floor");

			/*  linear option */
			if (p.nMatch(" LOG") || (StopCalc.TeFloor <= 10. && !p.nMatch(
					"LINE")))
				StopCalc.TeFloor = exp10(StopCalc.TeFloor);

			if (StopCalc.TeFloor < phycon.TEMP_LIMIT_LOW)
			{
				fprintf(ioQQQ, " TE < %gK, reset to %gK.\n",
						phycon.TEMP_LIMIT_LOW, 1.0001 * phycon.TEMP_LIMIT_LOW);
				StopCalc.TeFloor = realnum(1.0001 * phycon.TEMP_LIMIT_LOW);
			}

			if (StopCalc.TeFloor > phycon.TEMP_LIMIT_HIGH)
			{
				fprintf(ioQQQ,
						" TE > %gK. Cloudy cannot handle this. Bailing out.\n",
						phycon.TEMP_LIMIT_HIGH);
				cdEXIT(EXIT_FAILURE);
			}

			/* option to vary scale factor */
			if (optimize.lgVarOn)
			{
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy(optimize.chVarFmt[optimize.nparm],
						"SET TEMPERATURE FLOOR %f LOG");

				/* pointer to where to write */
				optimize.nvfpnt[optimize.nparm] = input.nRead;

				/* vary log of temperature  */
				optimize.vparm[0][optimize.nparm] = (realnum) log10(
						StopCalc.TeFloor);
				optimize.varang[optimize.nparm][0] = (realnum) log10(
						1.00001 * phycon.TEMP_LIMIT_LOW);
				optimize.varang[optimize.nparm][1] = (realnum) log10(
						0.99999 * phycon.TEMP_LIMIT_HIGH);
				optimize.vincr[optimize.nparm] = 0.3f;

				++optimize.nparm;
			}
		}

		else if (p.nMatch("CONV") || p.nMatch("TOLE"))
		{
			/* error tolerance in heating cooling match, number is error/total */
			conv.HeatCoolRelErrorAllowed = (realnum) p.FFmtRead();
			if (p.lgEOL())
			{
				p.NoNumb("heating cooling tolerance");
			}
			if (conv.HeatCoolRelErrorAllowed <= 0.)
			{
				conv.HeatCoolRelErrorAllowed = exp10(conv.HeatCoolRelErrorAllowed);
			}
		}

		else
		{
			fprintf(ioQQQ,
					"\nI did not recognize a keyword on this SET TEMPERATURE command.\n");
			p.PrintLine(ioQQQ);
			fprintf(ioQQQ, "The keywords are FLOOr and CONVergence.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch("TEST"))
	{
		/* set flag saying to turn on some test - this is in cddefines.h in the global namespace */
		lgTestCodeEnabled = true;
	}

	else if (p.nMatch("TRIM"))
	{
		/* set trim upper or lower, for ionization stage trimming
		 * ion trimming ionization trimming 
		 * in routine TrimStage */
		if (p.nMatch("UPPE"))
		{
			/* set trim upper - either set value or turn off */
			ionbal.trimhi = exp10(p.FFmtRead());
			ionbal.lgTrimhiOn = true;			
			if (p.lgEOL() && p.nMatch(" OFF"))
			{
				/* turn off upward trimming */
				p.setEOL(false);
				ionbal.lgTrimhiOn = false;
				/* reset high limit to proper value */
				ionbal.trimhi = SMALLFLOAT;
			}
		}

		else if (p.nMatch("LOWE"))
		{
			/* set trim lower */
			ionbal.trimlo = exp10(p.FFmtRead());
			ionbal.lgTrimloOn = true;			
			if (p.lgEOL() && p.nMatch(" OFF"))
			{
				/* turn off upward trimming */
				p.setEOL(false);
				ionbal.lgTrimloOn = false;
				/* reset low limit to proper value */
				ionbal.trimlo = SMALLFLOAT;
			}
		}

		/* turn off ionization stage trimming */
		else if (p.nMatch("SMAL") || p.nMatch(" OFF"))
		{
			/* set small limits to both upper and lower limit*/
			ionbal.lgTrimhiOn = false;
			ionbal.lgTrimloOn = false;
			ionbal.trimlo = SMALLFLOAT;
			ionbal.trimhi = SMALLFLOAT;
			p.setEOL(false);
		}

		/* use new trimming approach */
		else if (p.nMatch("NEW"))
		{
			ionbal.lgNewTrim = true;
			p.setEOL(false);
		}

		/* use old trimming approach */
		else if (p.nMatch("NEW"))
		{
			ionbal.lgNewTrim = false;
			p.setEOL(false);
		}

		else
		{
			/* set trim upper */
			ionbal.lgTrimhiOn = true;			
			ionbal.trimhi = exp10(p.FFmtRead());

			/* set trim lower to same number */
			ionbal.trimlo = ionbal.trimhi;
			ionbal.lgTrimloOn = true;			
		}

		if (p.lgEOL())
		{
			p.NoNumb("trimming parameter");
		}

		if (ionbal.trimlo >= 1. || ionbal.trimhi >= 1.)
		{
			fprintf(ioQQQ, " number must be negative since log\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch("SKIP"))
	{
		/* skip save command, for saving every n't point */
		save.ncSaveSkip = (long) p.FFmtRead();
		if (p.lgEOL())
		{
			p.NoNumb("number of points to skip in save");
		}
	}

	else if (p.nMatch(" UTA"))
	{
		if (p.nMatch("KISI"))
		{
			if (p.nMatch(" OFF"))
			{
				/* the default, do not use it */
				atmdat.lgInnerShell_Kisielius = false;
			}
			else if( p.nMatch( " ON " ) )
			{
				atmdat.lgInnerShell_Kisielius = true;
			}
			else
			{
				fprintf( ioQQQ, "Error: SET UTA KISIELIUS expects ON or OFF\n" );
				cdEXIT( EXIT_FAILURE );
			}
		}
		else if( p.nMatch( " OFF" ) )
		{
			/* turn off ALL inner shell absorption ionization */
			atmdat.lgInnerShellLine_on = false;
			phycon.lgPhysOK = false;
		}
	}

	else if (p.nMatch("WEAKH"))
	{
		/* set WeakHeatCool, threshold on save heating and cooling, default 0.05 */
		save.WeakHeatCool = (realnum) p.FFmtRead();

		if (p.lgEOL())
		{
			p.NoNumb("threshold on save heating and cooling");
		}

		if (save.WeakHeatCool < 0.)
		{
			save.WeakHeatCool
					= exp10(save.WeakHeatCool);
		}
	}

	else if (p.nMatch("KSHE"))
	{
		/* upper limit to opacities for k-shell ionization */
		continuum.EnergyKshell = (realnum) p.FFmtRead();
		if (p.lgEOL())
		{
			p.NoNumb("k-shell ionization opacity limit");
		}

		if (continuum.EnergyKshell == 0.)
		{
			/* arg of 0 causes upper limit to energy array */
			continuum.EnergyKshell = rfield.egamry();
		}

		else if (continuum.EnergyKshell < 194.)
		{
			fprintf(ioQQQ, " k-shell energy must be greater than 194 Ryd\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if (p.nMatch("NCHR"))
	{
		/* option to set the number of charge states for grain model */
		double val = p.FFmtRead();

		if (p.lgEOL())
		{
			p.NoNumb("number of charge states");
		}
		else
		{
			long nChrg = nint(val);
			if (nChrg < 2 || nChrg > NCHU)
			{
				fprintf(ioQQQ,
						" illegal value for number of charge states: %ld\n",
						nChrg);
				fprintf(ioQQQ, " choose a value between 2 and %d\n", NCHU);
				fprintf(ioQQQ,
						" or increase NCHS in grainvar.h and recompile\n");
				cdEXIT(EXIT_FAILURE);
			}
			else
			{
				SetNChrgStates(nChrg);
			}
		}
	}

	else if (p.nMatch("NEGO"))
	{
		/* save negative opacities if they occur, set negopac */
		opac.lgNegOpacIO = true;
	}

	else if (p.nMatch("NEND"))
	{
		/* default limit to number of zones to be computed
		 * only do this if nend is NOT currently left at its default
		 * nend is set to nEndDflt in routine zero
		 * this command only has effect if stop zone not entered */
		if (iterations.lgEndDflt)
		{
			/* this is default limit to number of zones, change it to this value */
			iterations.nEndDflt = (long) p.FFmtRead();
			iterations.lgEndDflt = false;
			if (p.lgEOL())
				p.NoNumb("limit to zone number");

			/* now change all limits, for all iterations, to this value */
			for (long i = 0; i < iterations.iter_alloc; i++)
				iterations.nend[i] = iterations.nEndDflt;

			if (iterations.nEndDflt > 2000)
				fprintf(
						ioQQQ,
						"CAUTION - it will take a lot of memory to save"
							" results for %li zones.  Is this many zones really necessary?\n",
						iterations.nEndDflt);
		}
	}

	else if (p.nMatch("TSQD"))
	{
		/* upper limit for highest density considered in the 
		 * Peimbert-style t^2 section of the printout */
		peimbt.tsqden = (realnum) p.FFmtRead();

		if (p.lgEOL())
		{
			p.NoNumb("highest density in t^2 section of printout");
		}
		peimbt.tsqden = exp10(peimbt.tsqden);
	}

	else if (p.nMatch("NMAP"))
	{
		/* how many steps in plot or save of heating-cooling map */
		hcmap.nMapStep = (long) p.FFmtRead();
		if (p.lgEOL())
		{
			p.NoNumb("steps in heating-cooling map");
		}
	}

	else if (p.nMatch("NUME") && p.nMatch("DERI"))
	{
		/* this is an option to use numerical derivatives for heating and cooling */
		NumDeriv.lgNumDeriv = true;
	}

	else if (p.nMatch("PATH"))
	{
		fprintf(ioQQQ, " The SET PATH command is no longer supported.\n");
		fprintf(ioQQQ, " Please set the correct path using the environment variable CLOUDY_DATA_PATH.\n");
		cdEXIT(EXIT_FAILURE);
	}

	else if (p.nMatch("PHFI"))
	{
		/* which version of PHFIT to use, 1995 or 1996 */
		ip = (long) p.FFmtRead();
		if (p.lgEOL())
		{
			p.NoNumb("version of PHFIT");
		}

		if (ip == 1995)
		{
			/* option to go back to old results, pre op */
			t_ADfA::Inst().set_version(PHFIT95);
		}
		else if (ip == 1996)
		{
			/* default is to use newer results, including opacity project */
			t_ADfA::Inst().set_version(PHFIT96);
		}
		else
		{
			fprintf(ioQQQ, " Two possible values are 1995 and 1996.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* set save xxx command */
	else if (p.nMatch("PUNC") || p.nMatch("SAVE"))
	{
		if (p.nMatch("HASH"))
		{
			/* to use the hash option there must be double quotes on line - were there? */
			if (!lgQuotesFound)
			{
				fprintf(ioQQQ,
						" I didn\'t recognize a key on this SET SAVE HASH line.\n");
				cdEXIT(EXIT_FAILURE);
			}
			/* specify the terminator between save output sets - normally a set of hash marks */
			/* 
			 * get any string within double quotes, and return it as
			 * null terminated string
			 * also sets name in OrgCard and chCard to blanks so that
			 * do not trigger off it later */
			if (chString_quotes_lowercase == "return")
			{
				/* special case - return becomes new line */
				save.chHashString = "\n";
			}
			else if (chString_quotes_lowercase == "time")
			{
				/* special case where output time between iterations */
				save.chHashString = "TIME_DEP";
			}
			else
			{
				/* usual case, simply copy what is in quotes */
				save.chHashString = chString_quotes_lowercase;
			}
		}

		else if ((p.nMatch("LINE") && p.nMatch("WIDT")) || p.nMatch("RESO")
				|| p.nMatch("LWID"))
		{
			/* set spectral resolution for contrast in continuum plots */
			if (p.nMatch(" C "))
			{
				fprintf(ioQQQ, " the keyword C is no longer necessary since"
					" energy conservation is now supported by default.\n");
				cdEXIT(EXIT_FAILURE);
			}
			else if (p.nMatch("SUPP"))
				/* suppress emission lines in save output */
				save.Resolution = 0_r;
			else
			{
				double number = p.FFmtRead();
				if (p.lgEOL())
					p.NoNumb("line width or resolution");
				if (number <= 0.)
				{
					fprintf(ioQQQ,
							" line width or resolution must be greater than zero.\n");
					cdEXIT(EXIT_FAILURE);
				}

				if (p.nMatch("RESO"))
					/* resolving power R = lambda/delta(lambda) */
					save.Resolution = realnum(number);
				else
					/* resolution FWHM in km/s */
					save.Resolution = realnum(SPEEDLIGHT / (number * 1.e5));
			}

			// option to do exactly the same thing with absorption lines
			// keyword is absorption 
			if (p.nMatch("ABSO"))
				save.ResolutionAbs = save.Resolution;
		}

		else if (p.nMatch("PREF"))
		{
			if( save.nsave > 0 )
			{
				fprintf(ioQQQ, " The SET SAVE PREFIX command should precede all save commands.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			if( lgQuotesFound )
			{
				// specify a prefix before all save filenames
				save.chFilenamePrefix = chString_quotes_lowercase;
			}
			else
			{
				fprintf( ioQQQ, " PROBLEM No prefix between double quotes was found on this line.\n" );
				cdEXIT( EXIT_FAILURE );
			}
		}

		else if (p.nMatch("FLUS"))
		{
			/* flush the output after every iteration */
			save.lgFLUSH = true;
		}

		else if (p.nMatch("LUMI") && p.nMatch(" OLD") )
		{
			/* use the old style of luminosity per inner cloud area in save continuum */
			save.lgLuminosityOld = true;
		}

		else
		{
			fprintf(ioQQQ,
				" There should have been an option on this command.\n");
			fprintf(ioQQQ,
				" Valid options for SET SAVE are summarized in Hazy 1 Miscellaneous Commands, and are:\n"
				" HASH, LINEWIDTH, RESOLUTION, PREFIX, FLUSH, LUMINOSITY OLD.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* set continuum .... options */
	else if (p.nMatch("CONT"))
	{
		if (p.nMatch("RESO"))
		{
			/* set resolution, get factor that will multiply continuum resolution that
			 * is contained in the file continuum_mesh.ini */
			(void)p.FFmtRead();
			if (p.lgEOL())
			{
				p.NoNumb("continuum resolution scale");
			}
			// no need to do anything further, the resolution scale factor was already set in cdRead()
		}

		else if (p.nMatch("SHIE"))
		{
			/* set continuum shielding function */
			/* these are all possible values of rt.nLineContShield,
			 * first is default, these are set with set continuum shielding */
			/*#define LINE_CONT_SHIELD_PESC	1*/
			/*#define LINE_CONT_SHIELD_FEDERMAN	2*/
			/*#define LINE_CONT_SHIELD_FERLAND	3*/
			if (p.nMatch("PESC"))
			{
				/* this uses an inward looking escape probability */
				rt.nLineContShield = LINE_CONT_SHIELD_PESC;
			}
			else if (p.nMatch("FEDE"))
			{
				/* set continuum shielding Federman,
				 * this is the default, and uses the appendix of
				 * >>refer	RT	continuum shielding	Federman, S.R., Glassgold, A.E., & 
				 * >>refercon	Kwan, J. 1979, ApJ, 227, 466*/
				/* ...with a bug in the wing opacity fixed */
				rt.nLineContShield = LINE_CONT_SHIELD_FEDERMAN;
			}
			else if (p.nMatch("FBUG"))
			{
				/* set continuum shielding Federman,
				 * this is the default, and uses the appendix of
				 * >>refer	RT	continuum shielding	Federman, S.R., Glassgold, A.E., & 
				 * >>refercon	Kwan, J. 1979, ApJ, 227, 466*/
				rt.nLineContShield = LINE_CONT_SHIELD_FEDERMAN_BUG;
			}
			else if (p.nMatch("FERL"))
			{
				rt.nLineContShield = LINE_CONT_SHIELD_FERLAND;
			}
			else if (p.nMatch("RODG"))
			{
				/* set continuum shielding
				 * >>refer Rodgers & Williams	1974JQSRT..14..319R 
				 * as in Draine & Bertoldi */
				rt.nLineContShield = LINE_CONT_SHIELD_RODGERS;
			}
			else if (p.nMatch("INTE"))
			{
				/* set continuum shielding, using best integral
				 * formalism available -- may be slow... */
				rt.nLineContShield = LINE_CONT_SHIELD_INTEGRAL;
			}
			else if (p.nMatch("NONE"))
			{
				/* turn off continuum shielding */
				rt.nLineContShield = 0;
			}
			else
			{
				fprintf(ioQQQ,
						" I didn\'t recognize a key on this SET CONTINUUM SHIELDing line.\n");
				fprintf(ioQQQ,
						" The ones I know about are: PESC, FEDErman, FERLand, RODGers and INTEgral.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		else
		{
			fprintf(ioQQQ,
					" I didn\'t recognize a key on this SET CONTINUUM line.\n");
			fprintf(ioQQQ,
					" The ones I know about are: RESOlution and SHIEld.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}

	else
	{
		fprintf(ioQQQ, " I didn\'t recognize a key on this SET command.\n");
		cdEXIT(EXIT_FAILURE);
	}

	return;
}
