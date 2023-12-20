/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseDont parse the dont command - do not do something */
#include "cddefines.h"
#include "taulines.h"
#include "opacity.h"
#include "phycon.h"
#include "secondaries.h"
#include "pressure.h"
#include "prt.h"
#include "coolheavy.h"
#include "thermal.h"
#include "dynamics.h"
#include "rt.h"
#include "yield.h"
#include "ionbal.h"
#include "atmdat.h"
#include "grainvar.h"
#include "hyperfine.h"
#include "save.h"
#include "parser.h"
#include "iso.h"
#include "mole.h"
#include "rfield.h"
#include "continuum.h"
#include "hydrogenic.h"

void ParseDont(Parser &p )
{
	DEBUG_ENTRY( "ParseDont()" );

	if( (p.nMatch( "21CM"  ) || p.nMatch( "21 CM"  )) && 
		p.nMatch( " LYA"  ) && p.nMatch( "PUMP"  ) )
	{
		/* no Lya 21 cm pump turns off 21 cm pumping of Lya */
		hyperfine.lgLya_pump_21cm = false;
	}
	else if( p.nMatch("ADVE") )
	{
		/*  turn off different aspects of advection */
		if( p.nMatch("H-LI") )
		{
			/* advection for the H-like ion sequence */
			dynamics.lgISO[ipH_LIKE] = 0;
		}
		else if( p.nMatch("HE-L") )
		{
			/* advection for the He-like ion sequence */
			dynamics.lgISO[ipHE_LIKE] = 0;
		}
		else if( p.nMatch("META") )
		{
			/* advection for the everything else - those done in ion_solver */
			dynamics.lgMETALS = 0;
		}
		else if( p.nMatch("COOL") )
		{
			/* turn off cooling - heating due to  advection */
			dynamics.lgCoolHeat = 0;
		}
		else
		{
			/* no sub option, so turn them all off */
			dynamics.lgISO[ipH_LIKE] = 0;
			dynamics.lgISO[ipHE_LIKE] = 0;
			dynamics.lgMETALS = 0;
		}

	}

	else if( p.nMatch("AUGE") )
	{
		/*  turn off auger effect by killing its block data */
		t_yield::Inst().kill_yield();
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("BLEN") )
	{
		// nothing to do here, has already been parsed in cdRead()
		(void)0;
	}

	else if( p.nMatch("BUFF") )
	{
		/* NO BUFFERING turn off buffered io for standard output, 
		 * used to get output when code crashes */

		/* >>chng 06 jun 28, moved handling of NO BUFFERING command to cdRead, PvH */

		/* stderr may be a preprocessor macro, so lets be really careful here */
		FILE *test = stderr;
		if( ioQQQ != test && save.chOutputFile.empty() )
		{
			/* this should really say stdout and not stderr ! */
			fprintf( ioQQQ, " ignored NO BUFFERING command since it could not be done safely.\n" );
		}

	}

	else if( p.nMatch("CHAR") )
	{
		/* turn off all charge transfer interactions */
		atmdat.lgCTOn = false;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("CTHE") )
	{
		/* turn off charge transfer heating */
		atmdat.HCharHeatOn = 0.;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("COMP") )
	{
		/* turn off both recoil ionization and compton heating of free electron */
		rfield.lgComptonOn = false;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("FEII") )
	{
		/* turn off feii ly-alpha pumping - rate evaluated in FeIILyaPump */
		hydro.lgLyaFeIIPumpOn = false;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("FILE") && p.nMatch("OPAC") )
	{
		/* no file opacities, generate them on the fly even if file present */
		opac.lgUseFileOpac = false;
	}

	else if( p.nMatch("FINE") && p.nMatch("OPAC") )
	{
		/* no fine opacities */
		rfield.lgOpacityFine = false;
	}

	else if( p.nMatch("FINE") )
	{
		/* turn off fine structure optical depths */
		rt.lgFstOn = false;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("FREE") )
	{
		/* turn off free free heating and cooling */
		CoolHeavy.lgFreeOn = false;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("GRAI") )
	{
		if( p.nMatch("NEUT") )
		{
			/* turn off ion grain recombination "NO GRAIN NEUTRALIZATION" */
			ionbal.lgGrainIonRecom = false;
			phycon.lgPhysOK = false;
		}
		else if( p.nMatch("GAS ") && p.nMatch("COLL") && p.nMatch("ENER") )
		{
			/* turn off grain - gas collisional energy exchange 
			 * "NO GRAIN GAS COLLISIONAL ENERGY EXCHANGE " */
			gv.lgDColOn = false;
			phycon.lgPhysOK = false;
		}
		else if( p.nMatch("ELEC") )
		{
			/* turn off grain contributions to electrons "NO GRAIN ELECTRONS" */
			gv.lgGrainElectrons = false;
			phycon.lgPhysOK = false;
		}
		else if( p.nMatch("MOLE") )
		{
			/* turn off capture of molecules on grain surfaces "NO GRAIN MOLECULES" */
			mole_global.lgGrain_mole_deplete = false;
			phycon.lgPhysOK = false;
		}
		else if( p.nMatch("QHEA") )
		{
			/* turn off quantum heating of grains "NO GRAIN QHEAT" */
			gv.lgQHeatOn = false;
			phycon.lgPhysOK = false;
		}
		else if( p.nMatch("X-RA") )
		{
			/* revert to WD01 physics "NO GRAIN X-RAY TREATMENT" */
			gv.lgWD01 = true;
		}
		else if( p.nMatch("PHYSICS") )
		{
			/* turn off grain physics "NO GRAIN PHYSICS" */
			gv.lgGrainPhysicsOn = false;
			phycon.lgPhysOK = false;
		}
		else
		{
			fprintf( ioQQQ, " No key recognized on this line.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* no induced processes */
	else if( p.nMatch("INDU") )
	{
		/* turn off induced recombination, stimulated emission,
		 * continuum fluorescent excitation of lines, 
		 * stimulated emission correction to optical depths attenuation */
		rfield.lgInducProcess = false;
	}

	/* no isotropic continua report */
	else if( p.nMatch("ISOT") && p.nMatch("CONT") && p.nMatch("REPO") )
	{
		continuum.lgPrtIsotropicCont = false;
	}

	/* no collisional ionization */
	else if( p.nMatch("COLL") && p.nMatch("IONI") )
	{
		/* turn off collisional ionization */
		atmdat.lgCollIonOn = false;
	}

	else if( p.nMatch("LEVE") )
	{
		/* turn off the set of level 2 lines, safe for lower densities 
		 * this is the upper limit to the counter that is always used,
		 * so no loops will ever occur */
		/* >>chng 06 mar 04 from -1 to 0 so that size_t is zero in state put & get */
		/*nWindLine = -1;*/
		nWindLine = 0;
	}

	/* various no line options */
	else if( p.nMatch("LINE") && !p.nMatch(" OTS") && !p.nMatch("OUTW") )
	{
		if( p.nMatch("DIFF") && p.nMatch("PUMP") )
		{
			/* no diffuse line pumping, 
			* turn off pumping of lines by diffuse continuum*/
			rfield.DiffPumpOn = 0.;
			fprintf( ioQQQ, " This option is disabled.\n Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		else if( p.nMatch("TRAN") )
		{
			/* no line transfer command */
			rfield.lgDoLineTrans = false;
		}
		else if( p.nMatch("ISOT") && p.nMatch("CONT") && p.nMatch("SUBT") )
		{
			/* do NOT subtract continuum from reported line fluxes & emissivities */
			save.lgSubtrCont = false;
		}
		else
		{
			/* missing no line option */
			fprintf( ioQQQ, " There has to be an option on the NO LINE command.\n" );
			fprintf( ioQQQ, " The options are 'DIFFUSE PUMP', 'TRANSFER',"
					" and 'ISOTROPIC CONTINUUM SUBTRACTION'.\n Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}


	else if( p.nMatch("OPAC") && p.nMatch("REEVAL") )
	{
		/* don't constantly reevaluate the opacities */
		rfield.lgOpacityReevaluate = false;
	}

	else if( p.nMatch("IONI") && p.nMatch("REEVAL") )
	{
		/* "no ionization reevaluation" - don't constantly reevaluate the ionization */
		rfield.lgIonizReevaluate = false;
	}

	/* options to kill ots components as debugging aids */
	else if( p.nMatch(" OTS") )
	{
		if( p.nMatch(" LYA") )
		{
			/* turn off Lya ots rates - for debugging oscillations */
			rfield.lgLyaOTS = false;
		}

		else if( p.nMatch("HEII") )
		{
			/* turn off Lya ots rates - for debugging oscillations */
			rfield.lgHeIIOTS = false;
		}

		else if( p.nMatch("LINE") )
		{
			/* turn off line ots rates - for debugging oscillations */
			rfield.lgKillOTSLine = true;
		}
	}

	/* options to kill outward compoents as a debugging aid */
	else if( p.nMatch("OUTW") )
	{
		if( p.nMatch("LINE") )
		{
			/* turn off Lya ots rates - for debugging oscillations */
			rfield.lgKillOutLine = true;
		}
		else if( p.nMatch("CONT") )
		{
			/* turn off Lya ots rates - for debugging oscillations */
			rfield.lgKillOutCont = true;
		}
	}
	else if( p.nMatch("MOLE") )
	{
		/* disable molecule formation, first option is to turn off only high Z part */
		if( p.nMatch("HEAV") )
		{
			/* turn off only Z>=2 molecules */
			mole_global.lgNoHeavyMole = true;
		}
		else
		{
			mole_global.lgNoMole = true;
		}
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("PHOT") )
	{
		/* disable photoionization */
		ionbal.lgPhotoIoniz_On = false;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("RADI") )
	{
		/* don't include line radiation pressure */
		pressure.lgLineRadPresOn = false;
	}

	else if( p.nMatch("RECO") )
	{
		/* disable compton recoil of bound electrons - "no recoil ioniz" */
		ionbal.lgCompRecoil = false;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("SCAT") && p.nMatch("OPAC"))
	{
		/* no scattering opacity, for Compton thick spherical geometry */
		opac.lgScatON = false;
	}

	else if( p.nMatch("SCAT") && p.nMatch("ESCA"))
	{
		/* no electron scattering contribution to line escape probs */
		rt.lgElecScatEscape = false;
	}

	else if( p.nMatch("ABSO") && p.nMatch("ESCA"))
	{
		/* NO ABSORPTION ESCAPE do not count opacity contribution to line escape probs - photon loss, actually */
		rt.lgAbsorLineEscape = false;
	}

	else if( p.nMatch("SCAT") && p.nMatch("INTEN"))
	{
		/* electron scattering is included as a contribution to line escape probs but not the reported line intensity */
		rt.lgElecScatInten = false;
	}

	else if( p.nMatch("SECO") )
	{
		/* turn off secondary electron ionizations */
		secondaries.lgSecOFF = true;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("SPOT") )
	{
		/* no on-the-spot; turn on all ground state rec */
		opac.otsmin = 1.;
	}

	else if( p.nMatch("STAR") )
	{
		/* no stark broadening */
		rt.lgStarkON = false;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("STAT") )
	{
		/* no static opacities - constantly reevaluate them */
		opac.lgOpacStatic = false;
	}

	else if( p.nMatch("TEPR") )
	{
		/* no tepredictor */
		/* turn off prediction of next zone's temperature, as guessed in ZoneStart */
		thermal.lgPredNextTe = false;
	}

	else if( p.nMatch("THRE") )
	{
		/* turn off Cota's three body rec subroutine */
		ionbal.lgNoCota = true;
		phycon.lgPhysOK = false;
	}

	else if( p.nMatch("TIME") )
	{
		/* don't print anything with a time, so that we can expect
		 * perfect agreement between separate runs */
		prt.lgPrintTime = false;
	}

	else if( p.nMatch(" UTA") )
	{
		fprintf( ioQQQ, "Obsolete command.  Please use SET UTA OFF instead.\n" );
		cdEXIT( EXIT_FAILURE );
	}

	/* the no vary command is parsed well before we get to this point,
	 * but we have to do something here or the parser will say that
	 * no command existed on the command line */
	else if( p.nMatch("VARY") )
	{
		/* this is a no-nothing, picked up to stop optimizer */
		((void)0);
	}

	else
	{
		/* end of else if trap */
		fprintf( ioQQQ," I do not recognize a keyword on this NO ... command.\n");
		p.PrintLine(ioQQQ);
		fprintf( ioQQQ, " Sorry.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* this option, if keyword (OK) appears, then do not set warning */
	if( p.nMatch("(OK)") )
	{
		/* say that physical conditions are actually ok */
		phycon.lgPhysOK = true;
	}
	return;
}
