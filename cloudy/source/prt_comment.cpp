/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtComment analyze model, generating comments on its features */
/*chkCaHeps check whether CaII K and H epsilon overlap */
/*prt_smooth_predictions check whether fluctuations in any predicted quantities occurred */
#include "cddefines.h"
#include "prt.h"
#include "cddrive.h"
#include "iso.h"
#include "continuum.h"
#include "stopcalc.h"
#include "hyperfine.h"
#include "grainvar.h"
#include "version.h"
#include "rt.h"
#include "he.h"
#include "taulines.h"
#include "hydrogenic.h"
#include "lines.h"
#include "trace.h"
#include "hcmap.h"
#include "hmi.h"
#include "save.h"
#include "h2.h"
#include "conv.h"
#include "dynamics.h"
#include "opacity.h"
#include "geometry.h"
#include "elementnames.h"
#include "ca.h"
#include "pressure.h"
#include "co.h"
#include "atoms.h"
#include "abund.h"
#include "colden.h"
#include "phycon.h"
#include "timesc.h"
#include "hextra.h"
#include "radius.h"
#include "iterations.h"
#include "fudgec.h"
#include "called.h"
#include "magnetic.h"
#include "wind.h"
#include "secondaries.h"
#include "struc.h"
#include "oxy.h"
#include "input.h"
#include "thermal.h"
#include "atmdat.h"
#include "warnings.h"
#include "mole.h"
#include "rfield.h"
#include "doppvel.h"
#include "freebound.h"
#include "two_photon.h"
#include "dense.h"
#include "lines_service.h"

/*chkCaHeps check whether CaII K and H epsilon overlap */
STATIC void chkCaHeps(double *totwid);

/*prt_smooth_predictions check whether fluctuations in any predicted quantities occurred */
STATIC void prt_smooth_predictions(void);

void PrtComment(void)
{
	char chLine[2*INPUT_LINE_LENGTH];
	
	string chLbl;

	bool lgThick;

	long int i,
	  imas, 
	  ipLo ,
	  ipHi ,
	  ipISO,
	  nelem, 
	  isav, 
	  j, 
	  nline, 
	  nneg;

	double big_ion_jump, 
	  absint ,
	  aj, 
	  alpha, 
	  big, 
	  bigm, 
	  comfrc, 
	  differ, 
	  error, 
	  flur, 
	  freqn, 
	  rate, 
	  ratio, 
	  rel, 
	  small, 
	  tauneg, 
	  ts ,
	  HBeta, 
	  relfl ,
	  relhm,  
	  fedest,
	  GetHeat, 
	  SumNeg, 
	  totwid;

	double VolComputed , VolExpected , ConComputed , ConExpected;

	DEBUG_ENTRY( "PrtComment()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " PrtComment called.\n" );
	}

	/* 
	 * enter all comments cautions warnings and surprises into a large
	 * stack of statements
	 * at end of this routine is master call to printing routines
	 */
	iterations.lgIterAgain = false;

	vector<module*>& mods=module_list::Inst().m_l;
	for (vector<module*>::iterator mi = mods.begin(); mi != mods.end(); ++mi)
	{
		(*mi)->comment(warnings);
	}

	/* say why calculation stopped */
	if( conv.lgBadStop )
	{
		/* this case stop probably was not intended */
		ostringstream oss;
		oss << " W-Calculation stopped because " << StopCalc.chReasonStop << " Iteration ";
		oss << iteration << " of " << iterations.itermx + 1;
		rgcin( oss.str() );
		sprintf( chLine, " W-Calculation stopped because %s", 
		  StopCalc.chReasonStop );
		warnin( chLine );
		warnin( " W-This was not intended." );
	}
	else
	{
		/* for iterate to convergence, print reason why it was not converged on 3rd and higher iterations */
		ostringstream oss;
		if( (conv.lgAutoIt && iteration != iterations.itermx + 1) && 
		  iteration > 2 )
		{
			oss << "   Calculation stopped because " << StopCalc.chReasonStop << " Iteration ";
			oss << iteration << " of " << iterations.itermx + 1 << ", not converged due to ";
			oss << conv.chNotConverged;
			rgcin( oss.str() );
		}
		else
		{
			oss << "   Calculation stopped because " << StopCalc.chReasonStop << " Iteration ";
			oss << iteration << " of " << iterations.itermx + 1;
			rgcin( oss.str() );
		}
	}

	/* check whether stopped because default number of zones hit,
	 * not intended?? */
	if( (!geometry.lgZoneSet) && geometry.lgZoneTrp )
	{
		conv.lgBadStop = true;
		sprintf( chLine, 
			" W-Calculation stopped because default number of zones reached.  Was this intended???" );
		warnin(chLine);
		sprintf( chLine, 
			" W-Default limit can be increased while retaining this check with the SET NEND command." );
		warnin(chLine);
	}

	/* check whether stopped because zones too thin - only happens when glob set
	 * and depth + dr = depth
	 * not intended */
	if( radius.lgDrMinUsed || radius.lgdR2Small )
	{
		conv.lgBadStop = true;
		sprintf( chLine, 
			" W-Calculation stopped zone thickness became too thin.   This was not intended." );
		warnin(chLine);
		sprintf( chLine, 
			" W-The most likely reason was an uncontrolled oscillation." );
		warnin(chLine);
		ShowMe();
	}

	if( radius.lgdR2Small )
	{
		sprintf( chLine, 
			" W-This happened because the globule scale became very small relative to the depth." );
		warnin(chLine);
		sprintf( chLine, 
			" W-This problem is described in Hazy." );
		warnin(chLine);
	}

	/* possible that last zone does not have stored temp - if so
	 * then make it up - this happens for some bad stops */
	ASSERT( nzone < struc.nzlim );

	if( struc.testr[nzone-1] == 0. && nzone > 1)
		struc.testr[nzone-1] = struc.testr[nzone-2];

	if( struc.ednstr[nzone-1] == 0. && nzone > 1)
		struc.ednstr[nzone-1] = struc.ednstr[nzone-2];

	/* give indication of geometry */
	rel = radius.depth/radius.rinner;
	if( rel < 0.1 )
	{
		rgcin( "   The geometry is plane-parallel." );
	}
	else if( rel >= 0.1 && rel < 3. )
	{
		rgcin( "   The geometry is a thick shell." );
	}
	else
	{
		rgcin( "   The geometry is spherical." );
	}
	/* levels of warnings: Warning   (possibly major problems)
	 *                     Caution   (not likely to invalidate the results)
	 *                     [      !] surprise, but not a problem
	 *                     [nothing] interesting note
	 */

	/* incorrect electron density detected */
	if( dense.lgEdenBad )
	{
		if( dense.nzEdenBad == nzone )
		{
			sprintf( chLine, " C-The assumed electron density was incorrect for the last zone." );
			caunin(chLine);
			sprintf( chLine, " C-Did a temperature discontinuity occur??" );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, " W-The assumed electron density was incorrect during the calculation.  This is bad." );
			warnin(chLine);
			ShowMe();
		}
	}

	/* thermal map was done but results were not ok */
	if( hcmap.lgMapDone && !hcmap.lgMapOK )
	{
		sprintf( chLine, "  !The thermal map had changes in slope - check map output." );
		bangin(chLine);
	}

	/* first is greater than zero if fudge factors were entered, second is
	 * true if fudge ever evaluated, even to see if fudge factors are in place */
	if( fudgec.nfudge > 0 || fudgec.lgFudgeUsed )
	{
		sprintf( chLine, "  !Fudge factors were used or were checked.  Why?" );
		bangin(chLine);
	}

	if( cpu.i().foundCSMismatch() )
	{
		sprintf( chLine, "  !Modified data files were used in this simulation."
				 " This is fine if it was done intentionally." );
		bangin(chLine);
	}

	if( dense.gas_phase[ipHYDROGEN] > 1.1e13 )
	{
		if( dense.gas_phase[ipHYDROGEN] > 1e15 )
		{
			sprintf( chLine, " C-Density greater than 10**15, heavy elements are very uncertain." );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, " C-Density greater than 10**13" );
			caunin(chLine);
		}
	}

	/* HBeta is used later in the code to check on line intensities */
	if( cdLine("Pump",Hbeta_WavLen,&relfl,&absint)<=0 )
	{
		fprintf( ioQQQ, " PROBLEM Did not find Pump H-beta, set to unity\n" );
		relfl = 1.;
		absint = 1.;
	}

	/* now find total Hbeta */
	/* >>chng from "totl" Hbeta which was a special entry, to "H  1" Hbeta, which 
	 * is the general case */
	if( cdLine( "H  1",wlAirVac(Hbeta_WavLen),&HBeta,&absint)<=0 )
	{
		fprintf( ioQQQ, " NOTE Did not find H  1 H-beta - set intensity to unity, "
			"will not check on importance of H 1 pumping.\n" );
		HBeta = 1.;
		absint = 1.;
	}
	else 
	{
		/* check on continuum pumping of Balmer lines */
		if( HBeta>SMALLFLOAT )
		{
			flur = relfl/HBeta;
			if( flur > 0.1 )
			{
				sprintf( chLine, "  !Continuum fluorescent production of H-beta was very important." );
				bangin(chLine);
			}
			else if(flur > 0.01 )
			{
				sprintf( chLine, "   Continuum fluorescent production of H-beta was significant." );
				notein(chLine);
			}
		}
	}

	// iterate to convergence - status of this iteration
	if( iteration > 1 && conv.lgAutoIt && iteration<iterations.itermx)
	{
		sprintf( chLine , "   Iteration not converged because %s.",
			conv.chNotConverged );
		notein(chLine);
	}

	/* advice about "iterate to convergence" no converging at end
	 * lgIterAgain -> not converged, not another iteration */
	if( conv.lgAutoIt && (iteration == iterations.itermx + 1) && !iterations.lgOpticalDepthonverged )
	{
		sprintf( chLine, " C-Iterate to convergence did not converge in %li iterations.",
			iteration );
		caunin(chLine);

		if( conv.lgAutoIt && (phycon.te < StopCalc.TempLoStopZone) )
		{
			sprintf( chLine, " C-Iterate to convergence requested but sim stopped due to too-low temperature.");
			caunin(chLine);
			sprintf( chLine, " C-This may have convergence problems due to outer radius changing.");
			caunin(chLine);
			sprintf( chLine, " C-Consider setting a different stop criterion.");
			caunin(chLine);
		}
	}

	/* check if there were problems with initial wind velocity */
	if( wind.lgBallistic() && ((!wind.lgWindOK) || wind.windv < 1e6) )
	{
		sprintf( chLine, " C-Wind velocity below sonic point; solution is not valid." );
		caunin(chLine);
	}

	/* now confirm that mass flux here is correct */
	if( !wind.lgStatic() )
	{
		rel = wind.emdot/(wind.windv*dense.gas_phase[ipHYDROGEN])/radius.r1r0sq;
		if( fabs(1.-rel)> 0.02 )
		{
			sprintf( chLine, " C-Wind mass flux error is %g%%",fabs(1.-rel)*100. );
			caunin(chLine);
			fprintf(ioQQQ,"DEBUG emdot\t%.3e\t%.3e\t%.3e\t%.3e\n",
				wind.emdot , wind.windv*dense.gas_phase[ipHYDROGEN],wind.windv,dense.gas_phase[ipHYDROGEN]);
		}
	}

	/* check that we didn't overrun zone scale */
	if( nzone >= struc.nzlim )
	{
		TotalInsanity();
	}

	// check on energy conservation
	ConserveEnergy();

	/* comments having to do with cosmic rays */
	/* comment if cosmic rays and magnetic field both present */
	if( hextra.cryden*magnetic.lgB > 0. )
	{
		sprintf( chLine, 
			"  !Magnetic field & cosmic rays both present.  Their interactions are not treated." );
		bangin(chLine);
	}

	/* comment if cosmic rays are not included and stop temp has been lowered to go into neutral gas */
	if( hextra.cryden== 0. && StopCalc.TempLoStopZone < phycon.TEMP_STOP_DEFAULT)
	{
		sprintf( chLine, 
			"  !Background cosmic rays are not included - is this physical?  It affects the chemistry." );
		bangin(chLine);
	}

	/* check whether cosmic rays on, but model thick to them */
	if( hextra.cryden > 0. && (findspecieslocal("H")->column/10. + findspecieslocal("H+")->column) > 1e23 )
	{
		sprintf( chLine, 
			" C-Model is thick to cosmic rays, which are on." );
		caunin(chLine);
	}

	/* was ionization rate less than cosmic ray ionization rate in ISM? */
	if( hextra.cryden == 0. && iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].gamnc < 1e-17 )
	{
		sprintf( chLine, 
			"  !Ionization rate fell below background cosmic ray ionization rate.  Should this be added too?" );
		bangin(chLine);
		sprintf( chLine, 
			"  !   Use the COSMIC RAY BACKGROUND command." );
		bangin(chLine);
	}

	/* PrtComment if test code is in place */
	if( lgTestCodeCalled && !t_version::Inst().lgRelease )
	{
		sprintf( chLine, "  !Test code is in place." );
		bangin(chLine);
	}

	/* lgComUndr set to .true. if Compton cooling rate underflows to 0 */
	if( rfield.lgComUndr )
	{
		sprintf( chLine, 
			"  !Compton cooling rate underflows to zero.  Is this important?" );
		bangin(chLine);
	}

	/* make note if input stream contained an underscore, which was converted into a space */
	if( input.lgUnderscoreFound )
	{
		sprintf( chLine, 
			"  !Some input lines contained underscores, these were changed to spaces." );
		bangin(chLine);
	}

	/* make note if input stream contained a left or right bracket, which was converted into a space */
	if( input.lgBracketFound )
	{
		sprintf( chLine, 
			"  !Some input lines contained [ or ], these were changed to spaces." );
		bangin(chLine);
	}

	/* lgHionRad set to .true. if no hydrogen ionizing radiation */
	if( rfield.lgHionRad )
	{
		sprintf( chLine, 
			"  !There is no hydrogen-ionizing radiation.  Was this intended?" );
		bangin(chLine);
	}

	/* check whether certain zones were thermally unstable */
	if( thermal.nUnstable > 0 )
	{
		sprintf( chLine, 
			"   Derivative of net cooling negative and so possibly thermally unstable in%4ld zones.", 
		  thermal.nUnstable );
		notein(chLine);
	}

	/* generate a bang if a large fraction of the zones were unstable */
	if( nzone > 1 && 
		(realnum)(thermal.nUnstable)/(realnum)(nzone) > 0.25 )
	{
		sprintf( chLine, 
			"  !A large fraction of the zones were possibly thermally unstable,%4ld out of%4ld", 
		  thermal.nUnstable, nzone );
		bangin(chLine);
	}

	/* comment if negative coolants were ever significant */
	if( thermal.CoolHeatMax > 0.2 )
	{
		sprintf( chLine, 
			"  !Negative cooling reached %6.1f%% of the local heating, due to %4.4s %.1f", 
		  thermal.CoolHeatMax*100., thermal.chCoolHeatMax, thermal.wlCoolHeatMax );
		bangin(chLine);
	}
	else if( thermal.CoolHeatMax > 0.05 )
	{
		sprintf( chLine, 
			"   Negative cooling reached %6.1f%% of the local heating, due to %4.4s %.2f", 
		  thermal.CoolHeatMax*100., thermal.chCoolHeatMax, thermal.wlCoolHeatMax );
		notein(chLine);
	}

	/* check if advection heating was important */
	if( dynamics.HeatMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Advection heating reached %.2f%% of the local heating.", 
		  dynamics.HeatMax*100. );
		bangin(chLine);
	}
	else if( dynamics.HeatMax > 0.005 )
	{
		sprintf( chLine, 
			"   Advection heating reached %.2f%% of the local heating.", 
		  dynamics.HeatMax*100. );
		notein(chLine);
	}

	/* check if advection cooling was important */
	if( dynamics.CoolMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Advection cooling reached %.2f%% of the local cooling.", 
		  dynamics.CoolMax*100. );
		bangin(chLine);
	}
	else if( dynamics.CoolMax > 0.005 )
	{
		sprintf( chLine, 
			"   Advection cooling reached %.2f%% of the local heating.", 
		  dynamics.CoolMax*100. );
		notein(chLine);
	}

	/* >>chng 06 mar 22, add this comment
	 * check if time dependent ionization front being done with too large a U */
	if( dynamics.lgTimeDependentStatic && dynamics.lgRecom )
	{
		if( rfield.uh > 1. )
		{
			sprintf( chLine, 
				" W-Time dependent ionization front cannot now handle strong-R cases - the ionization parameter is too large." );
			warnin(chLine);
		}
		else if( rfield.uh > 0.1 )
		{
			sprintf( chLine, 
				" C-Time dependent ionization front cannot now handle strong-R cases - the ionization parameter is too large." );
			caunin(chLine);
		}
	}

	/* check if thermal ionization of ground state of hydrogen was important */
	if( hydro.HCollIonMax > 0.10 )
	{
		sprintf( chLine, 
			"  !Thermal collisional ionization of H reached %.2f%% of the local ionization rate.", 
		  hydro.HCollIonMax*100. );
		bangin(chLine);
	}
	else if( hydro.HCollIonMax > 0.005 )
	{
		sprintf( chLine, 
			"   Thermal collisional ionization of H reached %.2f%% of the local ionization rate.", 
		  hydro.HCollIonMax*100. );
		notein(chLine);
	}

	/* check if lookup table for Hummer & Storey case B was exceeded */
	if( !atmdat.lgHCaseBOK[1][ipHYDROGEN]  )
	{
		sprintf( chLine, 
			"   Te-ne bounds of Case B lookup table exceeded, H I Case B line intensities set to zero." );
		notein(chLine);
	}
	if( !atmdat.lgHCaseBOK[1][ipHELIUM]  )
	{
		sprintf( chLine, 
			"   Te-ne bounds of Case B lookup table exceeded, He II Case B line intensities set to zero." );
		notein(chLine);
	}

	if( dense.EdenMax>1e8  )
	{
		sprintf( chLine, 
			"  !The high electron density makes the Nussbaumer/Storey CNO recombination predictions unreliable." );
		bangin(chLine);
	}

	/* check if secondary ionization of hydrogen was important */
	if( secondaries.SecHIonMax > 0.10 )
	{
		sprintf( chLine, 
			"  !Suprathermal collisional ionization of H reached %.2f%% of the local H ionization rate.", 
		  secondaries.SecHIonMax*100. );
		bangin(chLine);
	}
	else if( secondaries.SecHIonMax > 0.005 )
	{
		sprintf( chLine, 
			"   Suprathermal collisional ionization of H reached %.2f%% of the local H ionization rate.", 
		  secondaries.SecHIonMax*100. );
		notein(chLine);
	}

	/* check if H2 vib-deexcitation heating was important */
	if( hmi.HeatH2DexcMax > 0.05 )
	{
		sprintf( chLine, 
			"  !H2 vib deexec heating reached %.2f%% of the local heating.", 
		  hmi.HeatH2DexcMax*100. );
		bangin(chLine);
	}
	else if( hmi.HeatH2DexcMax > 0.005 )
	{
		sprintf( chLine, 
			"   H2 vib deexec heating reached %.2f%% of the local heating.", 
		  hmi.HeatH2DexcMax*100. );
		notein(chLine);
	}

	/* check if H2 vib-deexcitation heating was important */
	if( hmi.CoolH2DexcMax > 0.05 )
	{
		sprintf( chLine, 
			"  !H2 deexec cooling reached %.2f%% of the local heating.", 
		  hmi.CoolH2DexcMax*100. );
		bangin(chLine);
	}
	else if( hmi.CoolH2DexcMax > 0.005 )
	{
		sprintf( chLine, 
			"   H2 deexec cooling reached %.2f%% of the local heating.", 
		  hmi.CoolH2DexcMax*100. );
		notein(chLine);
	}

	/* check if charge transfer ionization of hydrogen was important */
	if( atmdat.HIonFracMax > 0.10 )
	{
		sprintf( chLine, 
			"  !Charge transfer H => H+ reached %.1f%% of the local H ionization rate.", 
		  atmdat.HIonFracMax*100. );
		bangin(chLine);
	}
	else if( atmdat.HIonFracMax > 0.005 )
	{
		sprintf( chLine, 
			"   Charge transfer H => H+ reached %.2f%% of the local H ionization rate.", 
		  atmdat.HIonFracMax*100. );
		notein(chLine);
	}

	/* check if charge transfer heating cooling was important */
	if( atmdat.HCharHeatMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Charge transfer heating reached %.2f%% of the local heating.", 
		  atmdat.HCharHeatMax*100. );
		bangin(chLine);
	}
	else if( atmdat.HCharHeatMax > 0.005 )
	{
		sprintf( chLine, 
			"   Charge transfer heating reached %.2f%% of the local heating.", 
		  atmdat.HCharHeatMax*100. );
		notein(chLine);
	}

	if( atmdat.HCharCoolMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Charge transfer cooling reached %.2f%% of the local heating.", 
		  atmdat.HCharCoolMax*100. );
		bangin(chLine);
	}
	else if( atmdat.HCharCoolMax > 0.005 )
	{
		sprintf( chLine, 
			"   Charge transfer cooling reached %.2f%% of the local heating.", 
		  atmdat.HCharCoolMax*100. );
		notein(chLine);
	}

	/* check whether photo from up level of Mg2 2798 ever important */
	if( atoms.xMg2Max > 0.1 )
	{
		sprintf( chLine, 
			"  !Photoionization of upper level of Mg II 2798 reached %.1f%% of the total Mg+ photo rate.", 
		  atoms.xMg2Max*100. );
		bangin(chLine);
	}
	else if( atoms.xMg2Max > 0.01 )
	{
		sprintf( chLine, 
			"   Photoionization of upper level of Mg II 2798 reached %.1f%% of the total Mg+ photo rate.", 
		  atoms.xMg2Max*100. );
		notein(chLine);
	}

	/* check whether photo from up level of [O I] 6300 ever important */
	if( oxy.poimax > 0.1 )
	{
		sprintf( chLine, 
			"  !Photoionization of upper levels of [O I] reached %.1f%% of the total O destruction rate.", 
		  oxy.poimax*100. );
		bangin(chLine);
	}
	else if( oxy.poimax > 0.01 )
	{
		sprintf( chLine, 
			"   Photoionization of upper levels of [O I] reached %.1f%% of the total O destruction rate.", 
		  oxy.poimax*100. );
		notein(chLine);
	}

	/* check whether photo from up level of [O III] 5007 ever important */
	if( (oxy.poiii2Max + oxy.poiii3Max) > 0.1 )
	{
		sprintf( chLine, 
			"  !Photoionization of upper levels of [O III] reached %.1f%% of the total O++ photo rate.", 
		  (oxy.poiii2Max + oxy.poiii3Max)*100. );
		bangin(chLine);
	}
	else if( (oxy.poiii2Max + oxy.poiii3Max) > 0.01 )
	{
		sprintf( chLine, 
			"   Photoionization of upper levels of [O III] reached %.1f%% of the total O++ photo rate.", 
		  (oxy.poiii2Max + oxy.poiii3Max)*100. );
		notein(chLine);
	}

	/* check whether photoionization of He 2trip S was important */
	if( he.frac_he0dest_23S > 0.1 )
	{
		sprintf( chLine, 
			"  !Destruction of He 2TriS reached %.1f%% of the total He0 dest rate"
			" at zone %li, %.1f%% of that was photoionization.", 
		  he.frac_he0dest_23S*100, 
		  he.nzone, 
		  he.frac_he0dest_23S_photo*100.  );
		bangin(chLine);
	}
	else if( he.frac_he0dest_23S > 0.01 )
	{
		sprintf( chLine, 
			"   Destruction of He 2TriS reached %.1f%% of the total He0 dest rate"
			" at zone %li, %.1f%% of that was photoionization.", 
		  he.frac_he0dest_23S*100, 
		  he.nzone, 
		  he.frac_he0dest_23S_photo*100.  );
		notein(chLine);
	}

	/* check for critical density for l-mixing */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		if( !iso_ctrl.lgCritDensLMix[ipISO] && dense.lgElmtOn[ipISO] )
		{
			sprintf( chLine,
				"   The density is too low to l-mix the lowest %s I collapsed level. "
				" More resolved levels are needed for accurate line ratios.",
				elementnames.chElementSym[ipISO]);
			notein(chLine);
		}
	}

	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* report continuum lowering for xx-like iso-sequence. */
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( iso_sp[ipISO][nelem].lgLevelsLowered && dense.lgElmtOn[nelem] )
			{
				sprintf( chLine, "  !Continuum was lowered into model %s-like %s due to high density.  Highest n is %li",
					elementnames.chElementSym[ipISO],
					elementnames.chElementSym[nelem],
					iso_sp[ipISO][nelem].n_HighestResolved_local+iso_sp[ipISO][nelem].nCollapsed_local);
				bangin(chLine);
			}
			else if( iso_sp[ipISO][nelem].lgLevelsEverLowered && dense.lgElmtOn[nelem] )
			{
				sprintf( chLine, "  !Continuum was lowered into model %s-like %s due to high density at SOME point but NOT at the last zone.",
					elementnames.chElementSym[ipISO],
					elementnames.chElementNameShort[nelem]);
				bangin(chLine);
			}
		}

		/* report pop rescaling for xx-like iso-sequence. */
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( iso_sp[ipISO][nelem].lgPopsRescaled )
			{
				ASSERT( dense.lgSetIoniz[nelem] );
				sprintf( chLine, " C-Populations were rescaled for %s-like %s due to \"element ionization\" command.",
					elementnames.chElementSym[ipISO],
					elementnames.chElementSym[nelem] );
				caunin(chLine);
			}
		}
	}

	/* frequency array may not have been defined for all energies */
	if( !rfield.lgMMok )
	{
		sprintf( chLine, 
			"  !Continuum not defined in extreme infrared - Compton scat, grain heating, not treated properly?" );
		bangin(chLine);
	}

	if( !rfield.lgHPhtOK )
	{
		sprintf( chLine, 
			"  !Continuum not defined at photon energies which ionize excited states of H, important for H- and ff heating." );
		bangin(chLine);
	}

	if( !rfield.lgXRayOK )
	{
		sprintf( chLine, 
			"  !Continuum not defined at X-Ray energies - Compton scattering and Auger ionization wrong?" );
		bangin(chLine);
	}

	if( !rfield.lgGamrOK )
	{
		sprintf( chLine, 
			"  !Continuum not defined at gamma-ray energies - pair production and Compton scattering OK?" );
		bangin(chLine);
	}

	if( continuum.lgCon0 )
	{
		sprintf( chLine, "  !Continuum zero at some energies." );
		bangin(chLine);
	}

	if( continuum.lgCoStarInterpolationCaution )
	{
		sprintf( chLine , " C-CoStarInterpolate interpolated between non-adjoining tracks, this may not be valid." );
		caunin(chLine);
	}

	if( rfield.lgOcc1Hi )
	{
		sprintf( chLine, 
			"  !The continuum occupation number at 1 Ryd is greater than unity." );
		bangin(chLine);
	}

	/* this flag set true it set dr forced first zone to be too big */
	if( radius.lgDR2Big )
	{
		sprintf( chLine, 
			" C-The thickness of the first zone was set larger than optimal by a SET DR command." );
		caunin(chLine);
		/* this is case where did one zone of specified thickness - but it 
		 * was too large */
		if( nzone<=1 )
			sprintf( chLine, 
			" C-Consider using the STOP THICKNESS command instead." );
		caunin(chLine);
	}

	/* check for plasma shielding */
	if( rfield.lgPlasNu )
	{
		sprintf( chLine, 
			"  !The largest plasma frequency was %.2e Ryd = %.2e micron  The continuum is set to 0 below this.", 
		  rfield.plsfrqmax,
		  /* wavelength in microns */
		  RYDLAM/rfield.plsfrqmax/1e4);
		bangin(chLine);
	}

	if( rfield.occmax > 0.1 )
	{
		if( rfield.occmnu > 1e-4 )
		{
			sprintf( chLine, 
				"  !The largest continuum occupation number was %.3e at %.3e Ryd.", 
			  rfield.occmax, rfield.occmnu );
			bangin(chLine);
		}
		else
		{
			/* not surprising if occupation number bigger than 1 around 1e-5 Ryd,
			 * since this is the case for 3K background */
			sprintf( chLine, 
				"   The largest continuum occupation number was %.3e at %.3e Ryd.", 
			  rfield.occmax, rfield.occmnu );
			notein(chLine);
		}
	}

	if( rfield.occmax > 1e4 && rfield.occ1nu > 0. )
	{
		/* occ1nu is energy (ryd) where continuum occupation number falls below 1 */
		if( rfield.occ1nu < 0.0912 )
		{
			sprintf( chLine, 
				"   The continuum occupation number fell below 1 at %.3e microns.", 
			  0.0912/rfield.occ1nu );
			notein(chLine);
		}
		else if( rfield.occ1nu < 1. )
		{
			sprintf( chLine, 
				"   The continuum occupation number fell  below 1 at %.3e Angstroms.", 
			  912./rfield.occ1nu );
			notein(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   The continuum occupation number fell  below 1 at %.3e Ryd.", 
			  rfield.occ1nu );
			notein(chLine);
		}
	}

	if( rfield.tbrmax > 1e3 )
	{
		sprintf( chLine, 
			"  !The largest continuum brightness temperature was %.3eK at %.3e Ryd.", 
		  rfield.tbrmax, rfield.tbrmnu );
		bangin(chLine);
	}

	if( rfield.tbrmax > 1e4 )
	{
		/* tbr4nu is energy (ryd) where continuum bright temp falls < 1e4 */
		if( rfield.tbr4nu < 0.0912 )
		{
			sprintf( chLine, 
				"   The continuum brightness temperature fell below 10000K at %.3e microns.", 
			  0.0912/rfield.tbr4nu );
			notein(chLine);
		}
		else if( rfield.tbr4nu < 1. )
		{
			sprintf( chLine, 
				"   The continuum brightness temperature fell below 10000K at %.3e Angstroms.", 
			  912./rfield.tbr4nu );
			notein(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   The continuum brightness temperature fell below 10000K at %.3e Ryd.", 
			  rfield.tbr4nu );
			notein(chLine);
		}
	}

	/* turbulence AND constant pressure do not make sense */
	if( DoppVel.TurbVel > 0. && strcmp(dense.chDenseLaw,"CPRE") == 0 )
	{
		sprintf( chLine, 
			"  !Both constant pressure and turbulence makes no physical sense?" );
		bangin(chLine);
	}

	/* filling factor AND constant pressure do not make sense */
	if( geometry.FillFac < 1. && strcmp(dense.chDenseLaw,"CPRE") == 0 )
	{
		sprintf( chLine, 
			"  !Both constant pressure and a filling factor makes no physical sense?" );
		bangin(chLine);
	}

	/* grains and solar abundances do not make sense */
	if( gv.lgDustOn() && abund.lgAbnReference )
	{
		sprintf( chLine, 
			"  !Grains are present, but the gas phase abundances were left at the solar default.  This is not physical." );
		bangin(chLine);
	}

	/* check if depletion command set but no grains, another silly thing to do */
	if( abund.lgDepln && !gv.lgDustOn() )
	{
		sprintf( chLine, 
			"  !Grains are not present, but the gas phase abundances were depleted.  This is not physical." );
		bangin(chLine);
	}

	/* check if grains depletion INCREASED the abundance of an element, a feature of Jenkins 2009 */
	for( long int i=0; i < LIMELM; i++ )
	{
		if( abund.DepletionScaleFactor[i] > 1. )
		{
			sprintf( chLine,
				 "  !The grain depletion scale factor was %.2e for %s",
				 abund.DepletionScaleFactor[i] , elementnames.chElementName[i] );
			bangin(chLine);

		}
	}

	if( gv.lgDustOn() )
	{
		long nBin=0L, nFail=0L;
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			if( gv.bin[nd].QHeatFailures > 0L )
			{
				++nBin;
				nFail += gv.bin[nd].QHeatFailures;
			}
		}
		if( nFail > 0 )
		{
			sprintf( chLine,
				 "  !The grain quantum heating treatment failed to converge %ld time(s) in %ld bin(s).", nFail, nBin );
			bangin(chLine);
		}
	}

#if 0
	/* check if PAHs were present in the ionized region */
	/* >>chng 05 jan 01, disabled this code now that PAH's have varying abundances by default, PvH */
	/** \todo	2	this statement needs to be reinstated with a better test for presence in the H II region */
	if( gv.lgDustOn() )
	{
		bool lgPAHsPresent_and_constant = false;
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			lgPAHsPresent_and_constant = lgPAHsPresent_and_constant || 
				/* it is ok to have PAHs in the ionized region if the abundances vary */
				(gv.bin[nd].lgPAHsInIonizedRegion /* && !gv.bin[nd]. lgDustVary */);
		}
		if( lgPAHsPresent_and_constant )
		{
			sprintf( chLine,
				 " C-PAH's were present in the ionized region, this has never been observed in H II regions." );
			caunin(chLine);
		}
	}
#endif

	/* constant temperature greater than continuum energy density temperature */
	if( thermal.lgTemperatureConstant && thermal.ConstTemp*1.0001 < phycon.TEnerDen )
	{
		sprintf( chLine, 
			" C-The continuum energy density temperature (%g K)"
			" is greater than the gas kinetic temperature (%g K).",
			phycon.TEnerDen , thermal.ConstTemp);
		caunin(chLine);
		sprintf( chLine, " C-This is unphysical." );
		caunin(chLine);
	}

	/* remark that grains not present but energy density was low */
	if( !gv.lgDustOn() && phycon.TEnerDen < 800. )
	{
		sprintf( chLine, 
			"   Grains were not present but might survive in this environment (energy density temperature was %.2eK)", 
		  phycon.TEnerDen );
		notein(chLine);
	}

	/* call routine that will check age of cloud */
	AgeCheck();

	/* check on Ca H and H-epsilon overlapping
	 * need to do this since need to refer to lines arrays */
	chkCaHeps(&totwid);
	if( totwid > 121. )
	{
		sprintf( chLine, "   H-eps and Ca H overlap." );
		notein(chLine);
	}

	/* warning that something was turned off */
	if( !phycon.lgPhysOK )
	{
		sprintf( chLine, "  !A physical process has been disabled." );
		bangin(chLine);
	}

	/* check on lifetimes of [O III] against photoionization, only for low den */
	if( dense.gas_phase[ipHYDROGEN] < 1e8 )
	{
		if( oxy.r5007Max > 0.0263f )
		{
			sprintf( chLine, 
				"  !Photoionization of upper level of [O III] 5007 reached %.2e%% of the radiative lifetime.", 
			  oxy.r5007Max*100. );
			bangin(chLine);
		}
		else if( oxy.r5007Max > 0.0263f/10.f )
		{
			sprintf( chLine, 
				"   Photoionization of upper level of [O III] 5007 reached %.2e%% of the radiative lifetime.", 
			  oxy.r5007Max*100. );
			notein(chLine);
		}
		if( oxy.r4363Max > 1.78f )
		{
			sprintf( chLine, 
				"  !Photoionization of upper level of [O III] 4363 reached %.2e%% of the radiative lifetime.", 
			  oxy.r4363Max*100. );
			bangin(chLine);
		}
		else if( oxy.r4363Max > 1.78f/10.f )
		{
			sprintf( chLine, 
				"   Photoionization of upper level of [O III] 4363 reached %.2e%% of the radiative lifetime.", 
			  oxy.r4363Max*100. );
			notein(chLine);
		}
	}

	/* check whether total heating and cooling matched
	 * >>chng 97 mar 28, added GrossHeat, heat in terms normally heat-cool */
	error = fabs(thermal.power-thermal.totcol)/SDIV((thermal.power + thermal.totcol)/2.);
	if( thermal.lgTemperatureConstant )
	{
		if( error > 0.05 )
		{
			sprintf( chLine, 
				"  !Heating - cooling mismatch =%5.1f%%. Caused by constant temperature assumption. ", 
			  error*100. );
			bangin(chLine);
		}
	}

	else
	{
		if( error > 0.05 && error < 0.2 )
		{
			sprintf( chLine, " C-Heating - cooling mismatch =%.1f%%. What\'s wrong?", 
			  error*100. );
			caunin(chLine);
		}
		else if( error >= 0.2 )
		{
			sprintf( chLine, " W-Heating - cooling mismatch =%.2e%%. What\'s wrong????", 
			  error*100. );
			warnin(chLine);
		}
	}

	/* say if Ly-alpha photo of Ca+ excited levels was important */
	if( ca.Ca2RmLya > 0.01 )
	{
		sprintf( chLine, 
			"   Photoionization of Ca+ 2D level by Ly-alpha reached %6.1f%% of the total rate out.", 
		  ca.Ca2RmLya*100. );
		notein(chLine);
	}

	/* check if Lya alpha ever hotter than gas */
	if( hydro.nLyaHot > 0 )
	{
		if( hydro.TLyaMax/hydro.TeLyaMax > 1.05 )
		{
			sprintf( chLine, 
				"  !The excitation temp of Lya exceeded the electron temp, largest value was %.2eK (gas temp there was %.2eK, zone %ld)", 
			  hydro.TLyaMax, hydro.TeLyaMax, hydro.nZTLaMax );
			bangin(chLine);
		}
	}

	/* check if line absorption heating was important */

	/* get all negative lines, check if line absorption significant heat source
	 * this is used in "final" for energy budget print out */
	if( cdLine("Line",0,&SumNeg,&absint)<=0 )
	{
		fprintf( ioQQQ, " did not get sumneg cdLine\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* this is total heating */
	if( cdLine("TotH",0,&GetHeat,&absint)<=0 )
	{
		fprintf( ioQQQ, " did not get GetHeat cdLine\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( GetHeat > 0. )
	{
		SumNeg /= GetHeat;
		if( SumNeg > 0.1 )
		{
			sprintf( chLine, 
				"  !Line absorption heating reached %.2f%% of the global heating.", 
			  SumNeg*100. );
			bangin(chLine);
		}
		else if( SumNeg > 0.01 )
		{
			sprintf( chLine, 
				"   Line absorption heating reached %.2f%% of the global heating.", 
			  SumNeg*100. );
			notein(chLine);
		}
	}

	if( input.lgSetNoBuffering )
	{
		sprintf( chLine, 
			"  !NO BUFFERING command was entered - this increases exec time by LARGE amounts.");
		bangin(chLine);
	}

	/* this is check of extra lines added with g-bar */
	if( thermal.GBarMax > 0.1 )
	{
		ASSERT( thermal.ipMaxExtra > 0 );
		chLbl = chLineLbl(TauLine2[thermal.ipMaxExtra-1]);

		sprintf( chLine, 
			"  !G-bar cooling lines reached %.2f%% of the local cooling.  Line=%.10s", 
					thermal.GBarMax*100., chLbl.c_str() );
		bangin(chLine);
	}

	else if( thermal.GBarMax > 0.01 )
	{
		chLbl = chLineLbl(TauLine2[thermal.ipMaxExtra-1]);

		sprintf( chLine, 
			"   G-bar cooling lines reached %.2f%% of the local cooling.  Line=%.10s", 
					thermal.GBarMax*100., chLbl.c_str() );
		notein(chLine);
	}

	/* this is check of hyperfine structure lines*/
	if( hyperfine.cooling_max > 0.1 )
	{
		sprintf( chLine, 
			"  !Hyperfine structure line cooling reached %.2f%% of the local cooling.", 
		  hyperfine.cooling_max*100.);
		bangin(chLine);
	}

	else if( hyperfine.cooling_max > 0.01 )
	{
		sprintf( chLine, 
			"   Hyperfine structure line cooling reached %.2f%% of the local cooling.", 
		  hyperfine.cooling_max*100. );
		notein(chLine);
	}

	/* line absorption heating reached more than 10% of local heating?
	 * HeatLineMax is largest heating(1,23)/htot */
	if( thermal.HeatLineMax > 0.1 )
	{
		long level = -1;
		TransitionProxy t = FndLineHt(&level);
		chLbl = chLineLbl(t);
		sprintf( chLine, 
			"  !Line absorption heating reached %.2f%% of the local heating - largest by level%2ld line %.10s", 
					thermal.HeatLineMax*100., level, chLbl.c_str() );
		bangin(chLine);
	}

	else if( thermal.HeatLineMax > 0.01 )
	{
		sprintf( chLine, 
			"   Line absorption heating reached %.2f%% of the local heating.", 
		  thermal.HeatLineMax*100. );
		notein(chLine);
	}

	/* check whether any lines in the iso sequences mased */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max. */
				long int nmax = iso_sp[ipISO][nelem].numLevels_local;

				/* minus one here is to exclude highest level */
				for( ipHi=1; ipHi < nmax - 1; ++ipHi )
				{
					for( ipLo=0; ipLo < ipHi; ++ipLo )
					{
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
							continue;

						/* did the line mase */
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauIn() < -0.1 )
						{
							sprintf( chLine, 
								"  !Some iso-structure lines mased: %s-like %s, line %li-%li had optical depth %.2e", 
								elementnames.chElementSym[ipISO],
								elementnames.chElementNameShort[nelem],
								ipHi , ipLo ,
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauIn() );
							bangin(chLine);
						}
					}
				}
			}
		}
	}

	if( dense.gas_phase[ipHYDROGEN] < 1e7 )
	{
		/* check on IR fine structure lines - not necessary if dense since will be in LTE */
		lgThick = false;
		tauneg = 0.;
		alpha = 0.;

		/* now print results, were any fine structure lines optically thick? */
		if( lgThick )
		{
			sprintf( chLine, 
				"  !Some infrared fine structure lines are optically thick:  largest tau was %.2e", 
			  alpha );
			bangin(chLine);
		}
		/* did any fine structure lines mase? */
		if( tauneg < -0.01 )
		{
			sprintf( chLine, 
				"  !Some fine structure lines mased: line %s had optical depth %.2e", 
						chLbl.c_str(), tauneg );
			bangin(chLine);
		}
	}

	/* were any other lines masing? */
	/* this is check that at least a second iteration was done with sphere static,
	 * the test is overridden with the (OK) option on the sphere static command,
	 * which sets geometry.lgStaticNoIt true */
	if( geometry.lgStatic && iterations.lgLastIt && (iteration == 1) && 
		!geometry.lgStaticNoIt)
	{
		sprintf( chLine, " C-I must iterate when SPHERE STATIC is set." );
		caunin(chLine);
		iterations.lgIterAgain = true;
	}

	/* caution if continuum is punched but only one iteration performed */
	if( save.lgPunContinuum && iteration == 1 && iterations.lgLastIt)
	{
		sprintf( chLine, " C-I must iterate when save continuum output is done." );
		caunin(chLine);
		iterations.lgIterAgain = true;
	}

	/** \todo	2	extend to all iso and elem */
	/* how important was induced two photon?? */
	{
		two_photon& tnu = iso_sp[ipH_LIKE][ipHYDROGEN].TwoNu[0];
		if( tnu.induc_dn_max > 1. )
		{
			sprintf( chLine, "  !Rate of induced H 2-photon emission reached %.2e s^-1", 
				tnu.induc_dn_max );
			bangin(chLine);
		}
		else if( tnu.induc_dn_max > 0.01 )
		{
			sprintf( chLine, "   Rate of induced H 2-photon emission reached %.2e s^-1", 
				tnu.induc_dn_max );
			notein(chLine);
		}
	}

	/* how important was induced recombination? */
	if( hydro.FracInd > 0.01 )
	{
		sprintf( chLine, 
			"   Induced recombination was %5.1f%% of the total for H level%3ld", 
		  hydro.FracInd*100., hydro.ndclev );
		notein(chLine);
	}

	if( hydro.fbul > 0.01 )
	{
		sprintf( chLine, 
			"   Stimulated emission was%6.1f%% of the total for H transition%3ld -%3ld", 
		  hydro.fbul*100., hydro.nbul + 1, hydro.nbul );
		notein(chLine);
	}

	/* check whether Fe II destruction of La was important - entry into lines stack 
	 * is in prt_lines_hydro.c */
	if( cdLine("Fe 2",1215.67,&fedest,&absint)<=0 )
	{
		fprintf( ioQQQ, " Did not find Fe II Lya\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* find total Lya for comparison */
	if( cdLine("H  1",1215.67,&relhm,&absint)<=0 )
	{
		fprintf( ioQQQ, " Did not find Lya\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( relhm > 0. )
	{
		ratio = fedest/(fedest + relhm);
		if( ratio > 0.1 )
		{
			sprintf( chLine, "  !Fe II destruction of Ly-a removed %.1f%% of the line.", 
			  ratio *100.);
			bangin(chLine);
		}
		else if( ratio > 0.01 )
		{
			sprintf( chLine, "   Fe II destruction of Ly-a removed %.1f%% of the line.", 
			  ratio );
			notein(chLine);
		}
	}

	if( cdLine("H-CT",6562.80,&relhm,&absint)<=0 )
	{
		fprintf( ioQQQ, " Comment did not find H-CT H-alpha\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( HBeta > 0. )
	{
		if( relhm/HBeta > 0.01 )
		{
			sprintf( chLine, 
				"  !Mutual neutralization production of H-alpha was significant." );
			bangin(chLine);
		}
	}

	/* note about very high population in H n=2 rel to ground, set in hydrogenic */
	if( hydro.lgHiPop2 )
	{
		sprintf( chLine, 
			"   The population of H n=2 reached %.2e relative to the ground state.", 
		  hydro.pop2mx );
		notein(chLine);
	}

	/* check where diffuse emission error */
	for( ipISO=ipH_LIKE; ipISO<=ipHE_LIKE; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( iso_sp[ipISO][nelem].CaseBCheck > 1.5 )
			{
				sprintf( chLine, 
					"   Ratio of computed diffuse emission to case B reached %g for iso %li element %li",
					iso_sp[ipISO][nelem].CaseBCheck , ipISO , nelem+1 );
				notein(chLine);
			}
		}
	}

	/* check whether electrons were relativistic */
	if( thermal.thist > 1e9 )
	{
		/* >>chng 06 feb 19, from 5e9 K for warning to 1e10K.  add test case at 1e10K
		 * and don't want warning in test suite.  nothing is wrong at this temp - eeff
		 * is in correctly for relativistic temps and will eventually dominate cooling */
		if( thermal.thist > 1.0001e10 )
		{
			sprintf( chLine, " W-Electrons were relativistic; High TE=%.2e", 
			  thermal.thist );
			warnin(chLine);
		}
		else
		{
			sprintf( chLine, " C-Electrons were mildly relativistic; High TE=%.2e", 
			  thermal.thist );
			caunin(chLine);
		}
	}

	/* check on timescale for photoerosion of elements */
	rate = timesc.TimeErode*2e-26;
	if( rate > 1e-35 )
	{
		/*  2E-26 is roughly cross section for photoerosion
		 *  see 
		 * >>refer	all	photoerode	Boyd, R., & Ferland, G.J. ApJ, 318, L21. */
		ts = (1./rate)/3e7;
		if( ts < 1e3 )
		{
			sprintf( chLine, "  !Timescale-photoerosion of Fe=%.2e yr", 
			  ts );
			bangin(chLine);
		}
		else if( ts < 1e9 )
		{
			sprintf( chLine, "   Timescale-photoerosion of Fe=%.2e yr", 
			  ts );
			notein(chLine);
		}
	}

	/* check whether Compton heating was significant */
	comfrc = rfield.comtot/SDIV(thermal.power);
	if( comfrc > 0.01 )
	{
		sprintf( chLine, "   Compton heating was %5.1f%% of the total.", 
		  comfrc*100. );
		notein(chLine);
	}

	/* check on relative importance of induced Compton heating */
	if( comfrc > 0.01 && rfield.cinrat > 0.05 )
	{
		sprintf( chLine, 
			"  !Induced Compton heating was %.2e of the total Compton heating.", 
		  rfield.cinrat );
		bangin(chLine);
	}

	/* check whether equilibrium timescales are short rel to Hubble time */
	if( timesc.tcmptn > 5e17 )
	{
		if( comfrc > 0.05 )
		{
			sprintf( chLine, 
				" C-Compton cooling is significant and the equilibrium timescale (%.2e s) is longer than the Hubble time.", 
			  timesc.tcmptn );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   Compton cooling equilibrium timescale (%.2e s) is longer than Hubble time.", 
			  timesc.tcmptn );
			notein(chLine);
		}
	}

	if( timesc.time_therm_long > 5e17 )
	{
		sprintf( chLine, 
			" C-Thermal equilibrium timescale, %.2e s, longer than Hubble time; this cloud is not time-steady.", 
		  timesc.time_therm_long );
		caunin(chLine);
	}

	/* check whether model large relative to Jeans length
	 * DMEAN is mean density (gm per cc)
	 * mean temp is weighted by mass density */
	if( log10(radius.depth) > colden.rjnmin )
	{
		/* AJMIN is minimum Jeans mass, log in grams */
		aj = exp10((double)colden.ajmmin - log10(SOLAR_MASS));
		if( strcmp(dense.chDenseLaw,"CPRE") == 0 )
		{
			sprintf( chLine, 
				" C-Cloud thicker than smallest Jeans length=%8.2ecm; stability problems? (smallest Jeans mass=%8.2eMo)", 
			  exp10(colden.rjnmin), aj );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   Cloud thicker than smallest Jeans length=%8.2ecm; stability problems? (smallest Jeans mass=%8.2eMo)", 
			  exp10(colden.rjnmin), aj );
			notein(chLine);
		}
	}

	/* check whether grains too hot to survive */
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		if( gv.bin[nd].nDustFunc != DF_SUBLIMATION && gv.bin[nd].TeGrainMax > gv.bin[nd].Tsublimat )
		{
			sprintf( chLine, 
			  " W-Maximum temperature of grain %s was %.2eK, above its sublimation temperature, %.2eK.", 
			  gv.bin[nd].chDstLab, gv.bin[nd].TeGrainMax, 
			  gv.bin[nd].Tsublimat );
			warnin(chLine);
		}
		else if( gv.bin[nd].nDustFunc != DF_SUBLIMATION && gv.bin[nd].TeGrainMax > gv.bin[nd].Tsublimat*0.9 )
		{
			sprintf( chLine, 
			  " C-Maximum temperature of grain %s was %.2eK, near its sublimation temperature, %.2eK.", 
			  gv.bin[nd].chDstLab, gv.bin[nd].TeGrainMax, 
			  gv.bin[nd].Tsublimat );
			caunin(chLine);
		}
	}

	if( gv.lgNegGrnDrg )
	{
		sprintf( chLine, "  !Grain drag force <0." );
		bangin(chLine);
	}

	/* largest relative number of electrons donated by grains */
	if( gv.GrnElecDonateMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Grains donated %5.1f%% of the total electrons in some regions.", 
		  gv.GrnElecDonateMax*100. );
		bangin(chLine);
	}
	else if( gv.GrnElecDonateMax > 0.005 )
	{
		sprintf( chLine, 
			"   Grains donated %5.1f%% of the total electrons in some regions.", 
		  gv.GrnElecDonateMax*100. );
		notein(chLine);
	}

	/* largest relative number of electrons on grain surface */
	if( gv.GrnElecHoldMax > 0.05 )
	{
		sprintf( chLine, 
			"  !Grains contained %5.1f%% of the total electrons in some regions.", 
		  gv.GrnElecHoldMax*100. );
		bangin(chLine);
	}
	else if( gv.GrnElecHoldMax > 0.005 )
	{
		sprintf( chLine, 
			"   Grains contained %5.1f%% of the total electrons in some regions.", 
		  gv.GrnElecHoldMax*100. );
		notein(chLine);
	}

	/* is photoelectric heating of gas by photoionization of grains important */
	if( gv.dphmax > 0.5 )
	{
		sprintf( chLine, 
			"  !Local grain-gas photoelectric heating rate reached %5.1f%% of the total.", 
		  gv.dphmax*100. );
		bangin(chLine);
	}
	else if( gv.dphmax > 0.05 )
	{
		sprintf( chLine, 
			"   Local grain-gas photoelectric heating rate reached %5.1f%% of the total.", 
		  gv.dphmax*100. );
		notein(chLine);
	}

	if( gv.TotalDustHeat/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, 
			"   Global grain photoelectric heating of gas was%5.1f%% of the total.", 
		  gv.TotalDustHeat/thermal.power*100. );
		notein(chLine);
		if( gv.TotalDustHeat/thermal.power > 0.25 )
		{
			sprintf( chLine, 
				"  !Grain photoelectric heating is VERY important." );
			bangin(chLine);
		}
	}

	/* grain-gas collisional cooling of gas */
	if( gv.dclmax > 0.05 )
	{
		sprintf( chLine, 
			"   Local grain-gas cooling of gas rate reached %5.1f%% of the total.", 
		  gv.dclmax*100. );
		notein(chLine);
	}

	/* check how H2 chemistry network performed */
	if( h2.renorm_max > 1.05 )
	{
		if( h2.renorm_max > 1.2 )
		{
			sprintf( chLine, 
				"  !The large H2 molecule - main chemistry network renormalization factor reached %.2f.", 
				h2.renorm_max);
			bangin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   The large H2 molecule - main chemistry network renormalization factor reached %.2f.", 
				h2.renorm_max);
			notein(chLine);
		}
	}
	if( h2.renorm_min < 0.95 )
	{
		if( h2.renorm_min < 0.8 )
		{
			sprintf( chLine, 
				"  !The large H2 molecule - main chemistry network renormalization factor reached %.2f.", 
				h2.renorm_min);
			bangin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   The large H2 molecule - main chemistry network renormalization factor reached %.2f.", 
				h2.renorm_min);
			notein(chLine);
		}
	}

	/* check whether photodissociation of H_2^+ molecular ion was important */
	if( hmi.h2pmax > 0.10 )
	{
		sprintf( chLine, 
			"  !The local H2+ photodissociation heating rate reached %5.1f%% of the total heating.", 
		  hmi.h2pmax*100. );
		bangin(chLine);
	}

	else if( hmi.h2pmax > 0.01 )
	{
		sprintf( chLine, 
			"   The local H2+ photodissociation heating rate reached %.1f%% of the total heating.", 
		  hmi.h2pmax*100. );
		notein(chLine);
	}

	/* check whether photodissociation of molecular hydrogen (H2)was important */
	if( hmi.h2dfrc > 0.1 )
	{
		sprintf( chLine, 
			"  !The local H2 photodissociation heating rate reached %.1f%% of the total heating.", 
		  hmi.h2dfrc*100. );
		bangin(chLine);
	}
	else if( hmi.h2dfrc > 0.01 )
	{
		sprintf( chLine, 
			"   The local H2 photodissociation heating rate reached %.1f%% of the total heating.", 
		  hmi.h2dfrc*100. );
		notein(chLine);
	}

	/* check whether cooling by molecular hydrogen (H2) was important */
	if( hmi.h2line_cool_frac > 0.1 )
	{
		sprintf( chLine, 
			"  !The local H2 cooling rate reached %.1f%% of the local cooling.", 
		  hmi.h2line_cool_frac*100. );
		bangin(chLine);
	}
	else if( hmi.h2line_cool_frac > 0.01 )
	{
		sprintf( chLine, 
			"   The local H2 cooling rate reached %.1f%% of the local cooling.", 
		  hmi.h2line_cool_frac*100. );
		notein(chLine);
	}

	if( hmi.h2dtot/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, 
			"   Global H2 photodissociation heating of gas was %.1f%% of the total heating.", 
		  hmi.h2dtot/thermal.power*100. );
		notein(chLine);
		if( hmi.h2dtot/thermal.power > 0.25 )
		{
			sprintf( chLine, "   H2 photodissociation heating is VERY important." );
			notein(chLine);
		}
	}

	/* check whether photodissociation of carbon monoxide (co) was important */
	if( co.codfrc > 0.25 )
	{
		sprintf( chLine, 
			"  !Local CO photodissociation heating rate reached %.1f%% of the total.", 
		  co.codfrc*100. );
		bangin(chLine);
	}
	else if( co.codfrc > 0.05 )
	{
		sprintf( chLine, 
			"   Local CO photodissociation heating rate reached %.1f%% of the total.", 
		  co.codfrc*100. );
		notein(chLine);
	}

	if( co.codtot/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, 
			"   Global CO photodissociation heating of gas was %.1f%% of the total.", 
		  co.codtot/thermal.power*100. );
		notein(chLine);
		if( co.codtot/thermal.power > 0.25 )
		{
			sprintf( chLine, "   CO photodissociation heating is VERY important." );
			notein(chLine);
		}
	}

	if( thermal.lgEdnGTcm )
	{
		sprintf( chLine, 
			"   Energy density of radiation field was greater than the Compton temperature. Is this physical?" );
		notein(chLine);
	}

	/* was cooling due to induced recombination important? */
	if( hydro.cintot/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, "   Induced recombination cooling was %.1f%% of the total.", 
		  hydro.cintot/thermal.power*100. );
		notein(chLine);
	}

	/* check whether free-free heating was significant */
	if( thermal.FreeFreeTotHeat/SDIV(thermal.power) > 0.1 )
	{
		sprintf( chLine, "  !Free-free heating was %.1f%% of the total.", 
		  thermal.FreeFreeTotHeat/thermal.power*100. );
		bangin(chLine);
	}
	else if( thermal.FreeFreeTotHeat/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, "   Free-free heating was %.1f%% of the total.", 
		  thermal.FreeFreeTotHeat/thermal.power*100. );
		notein(chLine);
	}

	/* was heating due to H- absorption important? */
	if( hmi.hmitot/SDIV(thermal.power) > 0.01 )
	{
		sprintf( chLine, "   H- absorption heating was %.1f%% of the total.", 
		  hmi.hmitot/SDIV(thermal.power)*100. );
		notein(chLine);
	}

	/* water destruction rate was zero */
	if( mole_global.lgH2Ozer )
	{
		sprintf( chLine, "   Water destruction rate zero." );
		notein(chLine);
	}

	/* check for negative optical depths,
	 * optical depth in excited state helium lines */
	small = 0.;
	imas = 0;
	isav = 0;
	j = 0;
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* >>chng 06 aug 28, from numLevels_max to _local. */
			for( ipLo=ipH2p; ipLo < (iso_sp[ipH_LIKE][nelem].numLevels_local - 1); ipLo++ )
			{
				/* >>chng 06 aug 28, from numLevels_max to _local. */
				for( ipHi=ipLo + 1; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_local; ipHi++ )
				{
					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauIn() < (realnum)small )
					{
						small = iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().TauIn();
						imas = ipHi;
						j = ipLo;
						isav = nelem;
					}
				}
			}
		}
	}

	if( small < -0.05 )
	{
		sprintf( chLine, 
			"  !Some hydrogenic lines mased, species was %2s%2ld, smallest tau was %.2e, transition %li-%li", 
			elementnames.chElementSym[isav], 
			isav+1,small, imas , j );
		bangin(chLine);
	}

	/* check for negative opacities */
	if( opac.lgOpacNeg )
	{
		sprintf( chLine, "  !Some opacities were negative - the SET NEGOPC command will save which ones." );
		bangin(chLine);
	}

	/* now check continua */
	small = 0.;
	imas = 0;
	isav = 0;
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* >>chng 06 aug 28, from numLevels_max to _local. */
			for( i=0; i < iso_sp[ipH_LIKE][nelem].numLevels_local; i++ )
			{
				if( opac.TauAbsGeo[0][iso_sp[ipH_LIKE][nelem].fb[i].ipIsoLevNIonCon-1] < -0.001 )
				{
					small = MIN2(small,(double)opac.TauAbsGeo[0][iso_sp[ipH_LIKE][nelem].fb[i].ipIsoLevNIonCon-1]);
					imas = i;
					isav = nelem;
				}
			}
		}
	}

	if( small < -0.05 )
	{
		sprintf( chLine, "  !Some hydrogenic (%2s%2ld) continua optical depths were negative; smallest=%.2e level=%3ld", 
			elementnames.chElementSym[isav], 
			isav+1,
		  small, imas );
		bangin(chLine);
	}

	/* check whether any continuum optical depths are negative */
	nneg = 0;
	tauneg = 0.;
	freqn = 0.;
	for( i=0; i < rfield.nflux; i++ )
	{
		if( opac.TauAbsGeo[0][i] < -0.001 )
		{
			nneg += 1;
			/* only remember the smallest freq, and most neg optical depth */
			if( nneg == 1 )
				freqn = rfield.anu(i);
			tauneg = MIN2(tauneg,(double)opac.TauAbsGeo[0][i]);
		}
	}

	if( nneg > 0 )
	{
		sprintf( chLine, "  !Some continuous optical depths <0.  The lowest freq was %.3e Ryd, and a total of%4ld", 
		  freqn, nneg );
		bangin(chLine);
		sprintf( chLine, "  !The smallest optical depth was %.2e", 
		  tauneg );
		bangin(chLine);
	}

	/* say if Balmer continuum optically thick */
	if( opac.TauAbsGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1] > 0.05 )
	{
		sprintf( chLine, "   The Balmer continuum optical depth was %.2e.", 
		  opac.TauAbsGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1] );
		notein(chLine);
	}

	/* was correction for stimulated emission significant? */
	if( opac.stimax[0] > 0.02 && opac.TauAbsGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1] > 0.2 )
	{
		sprintf( chLine, "   The Lyman continuum stimulated emission correction to optical depths reached %.2e.", 
		  opac.stimax[0] );
		notein(chLine);
	}
	else if( opac.stimax[1] > 0.02 && opac.TauAbsGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1] > 0.1 )
	{
		sprintf( chLine, "   The Balmer continuum stimulated emission correction to optical depths reached %.2e.", 
		  opac.stimax[1] );
		notein(chLine);
	}

	/* say if Paschen continuum optically thick */
	if( opac.TauAbsGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[3].ipIsoLevNIonCon-1] > 0.2 )
	{
		sprintf( chLine, 
			"   The Paschen continuum optical depth was %.2e.", 
		  opac.TauAbsGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[3].ipIsoLevNIonCon-1] );
		notein(chLine);
	}

	/* some comments about near IR total optical depth */
	if( opac.TauAbsGeo[0][0] > 1. )
	{
		sprintf( chLine, 
			"   The continuum optical depth at the lowest energy considered (%.3e Ryd) was %.3e.", 
		  rfield.anu(0), opac.TauAbsGeo[0][0] );
		notein(chLine);
	}

	/* comment if optical depth to Rayleigh scattering is big
	 * cs from VAL 76 */
	if( findspecieslocal("H")->column*7e-24 > 0.01 )
	{
		sprintf( chLine, 
			"   The optical depth to Rayleigh scattering at 1300A is %.2e", 
		  findspecieslocal("H")->column*6.71e-24 );
		notein(chLine);
	}

	if( findspecieslocal("H2+")->column*7e-18 > 0.1 )
	{
		sprintf( chLine, 
			"  !The optical depth to the H2+ molecular ion is %.2e", 
		  findspecieslocal("H2+")->column*7e-18 );
		bangin(chLine);
	}
	else if( findspecieslocal("H2+")->column*7e-18 > 0.01 )
	{
		sprintf( chLine, 
			"   The optical depth to the H2+ molecular ion is %.2e", 
		  findspecieslocal("H2+")->column*7e-18 );
		notein(chLine);
	}

	/* warn if optically thick to H- absorption */
	if( opac.thmin > 0.1 )
	{
		sprintf( chLine, 
			"  !Optical depth to negative hydrogen ion is %.2e", 
		  opac.thmin );
		bangin(chLine);
	}
	else if( opac.thmin > 0.01 )
	{
		sprintf( chLine, 
			"   Optical depth to negative hydrogen ion is %.2e", 
		  opac.thmin );
		notein(chLine);
	}

	/* check whether energy density less than background */
	if( phycon.TEnerDen < 2.6 )
	{
		sprintf( chLine, 
			"  !Incident radiation field energy density is less than 2.7K.  Add background with CMB command." );
		bangin(chLine);
	}

	/* check whether CMB set at all */
	if( !rfield.lgCMB_set )
	{
		sprintf( chLine, 
			"  !The CMB was not included.  This is added with the CMB command." );
		bangin(chLine);
	}

	/* incident radiation field is less than background Habing ISM field */
	if( rfield.lgHabing )
	{
		sprintf( chLine, 
			"  !The intensity of the incident radiation field is less than 10 times the Habing diffuse ISM field.  Is this OK?" );
		bangin(chLine);
		sprintf( chLine, 
			"  !   Consider adding diffuse ISM emission with TABLE ISM command." );
		bangin(chLine);
	}

	/* some things dealing with molecules, or molecule formation */

	/* if C/O > 1 then chemistry will be carbon dominated rather than oxygen dominated */
	if( dense.lgElmtOn[ipOXYGEN] && dense.lgElmtOn[ipCARBON] )
	{
		if( dense.gas_phase[ipCARBON]/dense.gas_phase[ipOXYGEN] > 1. )
		{
			sprintf( chLine, "  !The C/O abundance ratio, %.1f, is greater than unity.  The chemistry will be carbon dominated.", 
				dense.gas_phase[ipCARBON]/dense.gas_phase[ipOXYGEN] );
			bangin(chLine);
		}
	}

	bool lgLots_of_moles = false;
	bool lgLotsSolids = false;
	/* largest fraction in any molecule */
	for( i=0; i<mole_global.num_calc; ++i )
	{
		if( mole.species[i].location == NULL && ( mole_global.list[i]->isIsotopicTotalSpecies() || mole_global.list[i]->charge<0 ) )
		{
			if( mole.species[i].xFracLim > 0.1 )
			{
				sprintf( chLine, "  !The fraction of %s in %s reached %.1f%% at some point in the cloud.", 
								 mole.species[i].atomLim->label().c_str(),
								 mole_global.list[i]->label.c_str(),
								 mole.species[i].xFracLim*100. );
				bangin(chLine);
				lgLots_of_moles = true;
				/* check whether molecules are on grains */
				if( !mole_global.list[i]->lgGas_Phase )
					lgLotsSolids = true;
			}
			else if( mole.species[i].xFracLim>0.01 )
			{
				sprintf( chLine, "   The fraction of %s in %s reached %.2f%% at some point in the cloud.", 
								 mole.species[i].atomLim->label().c_str(),
								 mole_global.list[i]->label.c_str(),
								 mole.species[i].xFracLim*100. );
				notein(chLine);
				lgLots_of_moles = true;
				/* check whether molecules are on grains */
				if( !mole_global.list[i]->lgGas_Phase )
					lgLotsSolids = true;
			}
			else if( mole.species[i].xFracLim > 1e-3 )
			{
				sprintf( chLine, "   The fraction of %s in %s reached %.3f%% at some point in the cloud.", 
								 mole.species[i].atomLim->label().c_str(),
								 mole_global.list[i]->label.c_str(),
								 mole.species[i].xFracLim*100. );
				notein(chLine);
				/* check whether molecules are on grains */
				if( !mole_global.list[i]->lgGas_Phase )
					lgLotsSolids = true;
			}
		}
	}

	/* generate comment if molecular fraction was significant but some heavy elements are turned off */
	if( lgLots_of_moles )
	{
		/* find all elements that are turned off */
		for(nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		{
			/* >>chng 05 dec 23, add mole_global.lgElem_in_chemistry */
			if( !dense.lgElmtOn[nelem] )
			{
				/* this triggers if element turned off but it is part of co chem net */
				sprintf( chLine, 
					" C-Molecules are important, but %s, part of the chemistry network, is turned off.", 
					elementnames.chElementName[nelem] );
				caunin(chLine);
			}
#			if 0
				/* this element has been turned off - now check if part of chemistry */
				for( i=NUM_HEAVY_MOLEC+NUM_ELEMENTS; i<NUM_COMOLE_CALC; ++i )
				{
					if( nelem==mole_global.list[i].nelem_den )
					{
						/* this triggers if element turned off but it is part of co chem net */
						sprintf( chLine, 
							" C-Molecules are important, but %s, part of the chemistry network, is turned off.", 
							elementnames.chElementName[nelem] );
						caunin(chLine);
					}
				}
#			endif
		}
	}

	/* say if lots of molecules on grains,
	 * molecules with labels that *GR */
	if( lgLotsSolids ) 
	{
		sprintf( chLine, "  !A significant amount of molecules condensed onto grain surfaces." );
		bangin(chLine);
		sprintf( chLine, "  !These are the molecular species with \"grn\" above." );
		bangin(chLine);
	}

	/* bremsstrahlung optical depth */
	if( rfield.EnergyBremsThin > 0.09 )
	{
		sprintf( chLine, "  !The cloud is optically thick at optical wavelengths, extending to %.3e Ryd =%.3eA", 
		  rfield.EnergyBremsThin, RYDLAM/rfield.EnergyBremsThin );
		bangin(chLine);
	}
	else if( rfield.EnergyBremsThin > 0.009 )
	{
		sprintf( chLine, "   The continuum of the computed structure may be optically thick in the near infrared." );
		notein(chLine);
	}

	/* did model run away to very large radius? */
	if( radius.Radius > 1e23 && radius.Radius/radius.rinner > 10. )
	{
		sprintf( chLine, "   Is an outer radius of %.2e reasonable?", 
		  radius.Radius );
		notein(chLine);
	}

	/* following set true in RT_line_one_tauinc if maser capped at tau = -1 */
	if( rt.lgMaserCapHit )
	{
		sprintf( chLine, "   Laser maser optical depths capped in RT_line_one_tauinc." );
		notein(chLine);
	}

	/* following set true in adius_next if maser cap set dr */
	if( rt.lgMaserSetDR )
	{
		sprintf( chLine, "  !Line maser set zone thickness in some zones." );
		bangin(chLine);
	}

	/* lgPradCap is true if radiation pressure was capped on first iteration
	 * also check that this is a constant total pressure model */
	if( (pressure.lgPradCap && (strcmp(dense.chDenseLaw,"CPRE") == 0)) && 
	  pressure.lgPres_radiation_ON )
	{
		sprintf( chLine, "   Radiation pressure kept below gas pressure on this iteration." );
		notein(chLine);
	}

	if( pressure.RadBetaMax > 0.25 )
	{
		if( pressure.ipPradMax_line == 0 )
		{
			sprintf( chLine, 
				"  !The ratio of radiation to gas pressure reached %.2e at zone %li.  Caused by Lyman alpha.", 
			  pressure.RadBetaMax,
			  pressure.ipPradMax_nzone);
			bangin(chLine);
		}
		else
		{
			sprintf( chLine, 
				"  !The ratio of radiation to gas pressure reached %.2e at zone %li.  "
				"Caused by line number %ld, label %s", 
			  pressure.RadBetaMax, 
			  pressure.ipPradMax_nzone,
			  pressure.ipPradMax_line,
			  pressure.chLineRadPres.c_str() );
			bangin(chLine);
		}
	}

	else if( pressure.RadBetaMax > 0.025 )
	{
		if( pressure.ipPradMax_line == 0 )
		{
			sprintf( chLine, 
				"   The ratio of radiation to gas pressure reached %.2e at zone %li.  Caused by Lyman alpha.", 
			  pressure.RadBetaMax,
			  pressure.ipPradMax_nzone);
			notein(chLine);
		}
		else
		{
			sprintf( chLine, 
				"   The ratio of radiation to gas pressure reached %.2e at zone %li.  "
				"Caused by line number %ld, label %s", 
			  pressure.RadBetaMax, 
			  pressure.ipPradMax_nzone,
			  pressure.ipPradMax_line,
			  pressure.chLineRadPres.c_str() );
			notein(chLine);
		}
	}

	if( opac.telec >= 5. )
	{
		sprintf( chLine, " W-The model is optically thick to electron "
			"scattering; tau=%.2e  Cloudy is NOT intended for this regime.", 
		  opac.telec );
		warnin(chLine);
	}
	else if( opac.telec > 2.0  )
	{
		sprintf( chLine, " C-The model is moderately optically thick to electron scattering; tau=%.1f", 
		  opac.telec );
		caunin(chLine);
	}
	else if( opac.telec > 0.1  )
	{
		sprintf( chLine, "  !The model has modest optical depth to electron scattering; tau=%.2f", 
		  opac.telec );
		bangin(chLine);
	}
	else if( opac.telec > 0.01 )
	{
		sprintf( chLine, "   The optical depth to electron scattering is %.3f", 
		  opac.telec );
		notein(chLine);
	}

	/* optical depth to 21 cm */
	if( HFLines[0].Emis().TauIn() > 0.5 )
	{
		sprintf( chLine, "  !The optical depth in the H I 21 cm line is %.2e",HFLines[0].Emis().TauIn() );
		bangin(chLine);
	}

	/* comment if level2 lines are off - they are used to pump excited states
	 * of ground term by UV light */
	if( nWindLine==0 )
	{
		/* generate comment */
		sprintf( chLine, "  !The level2 lines are disabled.  UV pumping of excited levels within ground terms is not treated." );
		bangin(chLine);
	}

	/* check on optical depth convergence of all hydrogenic lines */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] && !dynamics.lgTimeDependentStatic )
		{
			if( iso_sp[ipH_LIKE][nelem].trans(ipH3p,ipH2s).Emis().TauIn() > 0.2 )
			{
				differ = fabs(1.-iso_sp[ipH_LIKE][nelem].trans(ipH3p,ipH2s).Emis().TauTot()/
					(iso_sp[ipH_LIKE][nelem].trans(ipH3p,ipH2s).Emis().TauIn()*rt.DoubleTau))*100.;

				/* check whether H-alpha optical depth changed by much on last iteration
				 * no tolerance can be finer than autocv, the tolerance on the
				 * iterate to convergence command.  It is 15% */
				if( ((iterations.lgLastIt && iso_sp[ipH_LIKE][nelem].trans(ipH3p,ipH2s).Emis().TauIn() > 0.8) && 
					differ > 20.) && wind.lgStatic() )
				{
					sprintf( chLine, 
						" C-This is the last iteration and %2s%2ld Bal(a) optical depth"
						" changed by %.1f%% (was %.2e). Use the ITERATE command to do more iterations.",
					  elementnames.chElementSym[nelem], 
					  nelem+1, differ, 
					  iso_sp[ipH_LIKE][nelem].trans(ipH3p,ipH2s).Emis().TauTot() );
					caunin(chLine);
					iterations.lgIterAgain = true;
				}

				/* only check on Lya convergence if Balmer lines are thick */
				if( iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn() > 0. )
				{
					differ = fabs(1.-iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot()/
						(iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn()*rt.DoubleTau))*100.;

					/* check whether Lya optical depth changed on last iteration
					 * no tolerance can be finer than autocv, the tolerance on the
					 * iterate to convergence command.  It is 15% */
					if( ((iterations.lgLastIt && iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn() > 0.8) && 
						differ > 25.) && wind.lgStatic() )
					{
						sprintf( chLine, 
							" C-This is the last iteration and %2s%2ld Ly(a) optical depth"
							" changed by %.1f%% (was %.2e). Use the ITERATE command to do more iterations.",
						elementnames.chElementSym[nelem], 
						  nelem+1,differ, iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() );
						caunin(chLine);
						iterations.lgIterAgain = true;
					}
				}
			}
		}
	}

	/* check whether sphere was set if dr/r large */
	if( radius.Radius/radius.rinner > 2. && !geometry.lgSphere )
	{
		sprintf( chLine, " C-R(out)/R(in)=%.2e and SPHERE was not set.", 
		  radius.Radius/radius.rinner );
		caunin(chLine);
	}

	/* check if thin in hydrogen or helium continua, but assumed to be thick */
	if( iterations.lgLastIt && !opac.lgCaseB )
	{

		/* check if thin in Lyman continuum, and assumed thick */
		if( rfield.nflux > iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon )
		{
			if( opac.TauAbsGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1] < 2. && 
				opac.TauAbsGeo[1][iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1] > 2. )
			{
				sprintf( chLine, " C-The H Lyman continuum is thin, and I assumed"
					" that it was thick.  Use the ITERATE command to do more iterations." );
				caunin(chLine);
				iterations.lgIterAgain = true;
			}
		}

		/* check on the He+ ionizing continuum */
		if( dense.lgElmtOn[ipHELIUM] && rfield.nflux > iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].ipIsoLevNIonCon )
		{
			if( (opac.TauAbsGeo[0][iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].ipIsoLevNIonCon-1] < 2. && 
				 opac.TauAbsGeo[1][iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].ipIsoLevNIonCon-1] > 2.) )
			{
				sprintf( chLine, 
					" C-The He II continuum is thin and I assumed that it was thick."
					"  Use the ITERATE command to do more iterations." );
				caunin(chLine);
				iterations.lgIterAgain = true;
			}
		}

		if( dense.lgElmtOn[ipHELIUM] && rfield.nflux > iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon )
		{ 
			if( (opac.TauAbsGeo[0][iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon-1] < 2. && 
				 opac.TauAbsGeo[1][iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon-1] > 2.) )
			{
				sprintf( chLine, 
					" C-The He I continuum is thin and I assumed that it was thick."
					"  Use the ITERATE command to do more iterations." );
				caunin(chLine);
				iterations.lgIterAgain = true;
			}
		}
	}

	/* check whether column density changed by much on this iteration */
	if( iteration > 1 )
	{
		if( colden.colden_old[ipCOL_HTOT] <= 0. )
		{
			fprintf( ioQQQ, " colden_old is insane in PrtComment.\n" );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}

		differ = fabs(1.-colden.colden[ipCOL_HTOT]/
			colden.colden_old[ipCOL_HTOT]);

		if( differ > 0.1 && differ <= 0.3 )
		{
			sprintf( chLine, 
				"   The H column density changed by %.2e%% between this and previous iteration.", 
			  differ*100. );
			notein(chLine);
		}

		else if( differ > 0.3 )
		{
			if( iterations.lgLastIt )
			{
				sprintf( chLine, 
					" C-The H column density changed by %.2e%% and this is the last iteration.  What happened?", 
				  differ*100. );
				caunin(chLine);
			}
			else
			{
				sprintf( chLine, 
					"  !The H column density changed by %.2e%%  What happened?", 
				  differ*100. );
				bangin(chLine);
			}
		}

		/* check on H2 column density, but only if significant fraction of H is molecular */
		if( (findspecieslocal("H2")->column+findspecieslocal("H2*")->column)/SDIV(colden.colden[ipCOL_HTOT]) > 1e-5 )
		{
			differ = fabs(1.-findspecieslocal("H2")->column/
							  SDIV(findspecieslocal("H2")->column_old));

			if( differ > 0.1 && differ <= 0.3 )
			{
				sprintf( chLine, 
					"   The H2 column density changed by %.2e%% between this and previous iteration.", 
				differ*100. );
				notein(chLine);
			}

			else if( differ > 0.3 )
			{
				if( iterations.lgLastIt )
				{
					sprintf( chLine, 
						" C-The H2 column density changed by %.2e%% and this is the last iteration.  What happened?", 
					differ*100. );
					caunin(chLine);
				}
				else
				{
					sprintf( chLine, 
						"  !The H2 column density changed by %.2e%%  What happened?", 
					differ*100. );
					bangin(chLine);
				}
			}
		}
	}

	/* say if rad pressure caused by la and la optical depth changed too much */
	differ = fabs(1.-iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().TauIn()/
		SDIV(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().TauTot()))*100.;

	if( iterations.lgLastIt && (pressure.RadBetaMax > 0.1) && 
		(differ > 50.) && (pressure.ipPradMax_line == 1) && (pressure.lgPres_radiation_ON) && 
		wind.lgStatic() )
	{
		sprintf( chLine, " C-This is the last iteration, radiation pressure was significant, and the L-a optical depth changed by %7.2f%% (was %.2e)", 
			differ, iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().TauTot() );
		caunin(chLine);
	}

	/* caution that 21 cm spin temperature is incorrect when Lya optical depth
	 * scale is overrun */
	if( iterations.lgLastIt &&
		( iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().TauTot() * 1.02 -
		iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().TauIn() ) < 0. )
	{
		sprintf( chLine, " C-The Lya optical depth scale was overrun and this is the last iteration - Tspin(21 cm) is not valid." );
		caunin(chLine);
		sprintf( chLine, " C-Another iteration is needed for Tspin(21 cm) to be valid.  Use the ITERATE command." );
		caunin(chLine);
	}

	/* say if la rad pressure capped by thermalization length */
	if( pressure.lgPradDen )
	{
		sprintf( chLine, "   Line radiation pressure capped by thermalization length." );
		notein(chLine);
	}

	/* print te failures */
	nline = MIN2(conv.nTeFail,10);
	if( conv.nTeFail != 0 )
	{
		if( conv.failmx < 0.1 )
		{
			long _o = sprintf( chLine, "   There were %ld minor temperature failures.  zones:", 
				conv.nTeFail );
			/* don't know how many zones we will save, there are nline,
			 * hence this use of pointer arith */
			for( i=0; i < nline; i++ )
			{
				_o += sprintf( chLine+_o, " %ld", conv.ifailz[i] );
			}
			notein(chLine);
		}
		else
		{
			sprintf( chLine, 
				"  !There were %ld temperature failures, and some were large. The largest was %.1f%%.  What happened?", 
			  conv.nTeFail, conv.failmx*100. );
			bangin(chLine);

			/* don't know how many zones we will save, there are nline,
			 * hence this use of pointer arith */
			long _o = sprintf( chLine , "  !The zones were" );
			for( i=0; i < nline; i++ )
			{
				_o += sprintf( chLine+_o, " %ld", conv.ifailz[i] );
			}
			bangin(chLine);

			if( struc.testr[0] > 8e4 && phycon.te < 6e5 )
			{
				sprintf( chLine, "  !I think they may have been caused by the change from hot to nebular gas phase.  The physics of this is unclear." );
				bangin(chLine);
			}
		}
	}

	/* check for temperature jumps */
	big_ion_jump = 0.;
	j = 0;
	for( i=1; i < nzone; i++ )
	{
		big = fabs(1.-struc.testr[i-1]/struc.testr[i]);
		if( big > big_ion_jump )
		{
			j = i;
			big_ion_jump = big;
		}
	}

	if( big_ion_jump > 0.2 )
	{
		/* this is a sanity check, but only do it if jump detected */
		if( j < 1 )
		{
			fprintf( ioQQQ, " j too small big jump check\n" );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}

		if( big_ion_jump > 0.4 )
		{
			sprintf( chLine, " C-A temperature discontinuity occurred"
				 " from %.2eK (zone %ld) to %.2eK (zone %ld).", 
				 struc.testr[j-1], j, struc.testr[j], j+1 );
			caunin(chLine);
			/* check if the second temperature is between 100 and 1000K */
			/* >>chng 05 nov 07, test second not first temperature since second
			 * will be lower of the two */
			/*if( struc.testr[j-1] < 1000. && struc.testr[j-1]>100. )*/
			if( struc.testr[j]>100. && struc.testr[j] < 1000. )
			{
				sprintf( chLine, " C-This was probably due to a thermal front." );
				caunin(chLine);
			}
		}
		else if( big_ion_jump > 0.2 )
		{
			sprintf( chLine, "  !A temperature discontinuity occurred"
				 " from %.2eK (zone %ld) to %.2eK (zone %ld).", 
				 struc.testr[j-1], j, struc.testr[j], j+1 );
			bangin(chLine);
			/* check if the second temperature is between 100 and 1000K */
			/* >>chng 05 nov 07, test second not first temperature since second
			 * will be lower of the two */
			/*if( struc.testr[j-1] < 1000. && struc.testr[j-1]>100. )*/
			if( struc.testr[j]>100. && struc.testr[j] < 1000. )
			{
				sprintf( chLine, "  !This was probably due to a thermal front." );
				bangin(chLine);
			}
		}
	}

	/* check for largest error in local electron density */
	if( fabs(conv.BigEdenError) > conv.EdenErrorAllowed )
	{
		/* this only produces a warning if not the very last zone */
		if( fabs(conv.BigEdenError) > conv.EdenErrorAllowed*20. && dense.nzEdenBad != 
		  nzone )
		{
			sprintf( chLine, " W-The local error in the electron density reached %.1f%% at zone %ld", 
			  conv.BigEdenError*100, dense.nzEdenBad );
			warnin(chLine);
		}
		else if( fabs(conv.BigEdenError) > conv.EdenErrorAllowed*5. )
		{
			sprintf( chLine, " C-The local error in the electron density reached %.1f%% at zone %ld", 
			  conv.BigEdenError*100, dense.nzEdenBad );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, "   The local error in the electron density reached %.1f%% at zone %ld", 
			  conv.BigEdenError*100, dense.nzEdenBad );
			notein(chLine);
		}
	}

	/* check for temperature oscillations or fluctuations*/
	big_ion_jump = 0.;
	j = 0;
	for( i=1; i < (nzone - 1); i++ )
	{
		big = fabs( (struc.testr[i-1] - struc.testr[i])/struc.testr[i] );
		bigm = fabs( (struc.testr[i] - struc.testr[i+1])/struc.testr[i] );

		/* this is sign of change in temperature, we are looking for change in sign */
		rel = ( (struc.testr[i-1] - struc.testr[i])/struc.testr[i])*
			( (struc.testr[i] - struc.testr[i+1])/struc.testr[i] );

		if( rel < 0. && MIN2( bigm , big ) > big_ion_jump )
		{
			j = i;
			big_ion_jump = MIN2( bigm , big );
		}
	}

	if( big_ion_jump > 0.1 )
	{
		/* only do sanity check if jump detected */
		if( j < 1 )
		{
			fprintf( ioQQQ, " j too small bigjump2 check\n" );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}

		if( big_ion_jump > 0.3 )
		{
			sprintf( chLine, 
				 " C-A temperature oscillation occurred by %.0f%%"
				 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
				 big_ion_jump*100., struc.testr[j-1], j, struc.testr[j], j+1,
				 struc.testr[j+1], j+2 );
			caunin(chLine);
		}
		else if( big_ion_jump > 0.1 )
		{
			sprintf( chLine, 
				"  !A temperature oscillation occurred by %.0f%%"
				 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
				 big_ion_jump*100., struc.testr[j-1], j, struc.testr[j], j+1,
				 struc.testr[j+1], j+2 );
			bangin(chLine);
		}
	}

	/* check for eden oscillations */
	if( strcmp(dense.chDenseLaw,"CDEN") == 0 )
	{
		j = 0;
		big_ion_jump = 0.;
		for( i=1; i < (nzone - 1); i++ )
		{
			big = (struc.ednstr[i-1] - struc.ednstr[i])/struc.ednstr[i];
			if( fabs(big) < conv.EdenErrorAllowed )
				big = 0.;
			bigm = (struc.ednstr[i] - struc.ednstr[i+1])/struc.ednstr[i];
			if( fabs(bigm) < conv.EdenErrorAllowed )
				bigm = 0.;
			if( big*bigm < 0. && 
				fabs(struc.ednstr[i-1]-struc.ednstr[i])/struc.ednstr[i] > big_ion_jump )
			{
				j = i;
				big_ion_jump = fabs(struc.ednstr[i-1]-struc.ednstr[i])/
				  struc.ednstr[i];
			}
		}

		/* only check on j if there was a big jump detected, number must be
		 * smallest jump */
		if( big_ion_jump > conv.EdenErrorAllowed*3. )
		{
			if( j < 1 )
			{
				fprintf( ioQQQ, " j too small bigjump3 check\n" );
				ShowMe();
				cdEXIT(EXIT_FAILURE);
			}

			if( big_ion_jump > conv.EdenErrorAllowed*10. )
			{
				sprintf( chLine, " C-An electron density oscillation occurred by %.0f%%"
					 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
					 big_ion_jump*100., struc.ednstr[j-1], j, struc.ednstr[j], j+1,
					 struc.ednstr[j+1], j+2 );
				caunin(chLine);
			}
			else if( big_ion_jump > conv.EdenErrorAllowed*3. )
			{
				sprintf( chLine, "  !An electron density oscillation occurred by %.0f%%"
					 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
					 big_ion_jump*100., struc.ednstr[j-1], j, struc.ednstr[j], j+1, 
					 struc.ednstr[j+1], j+2 );
				bangin(chLine);
			}
		}
	}

	/*prt_smooth_predictions check whether fluctuations in any predicted quantities occurred */
	/* >>chng 03 dec 05, add this test */
	prt_smooth_predictions();

	/**********************************************************
	 * check that the volume integrates out ok                *
	 **********************************************************/

	/* this was the number 1 fed through the line integrators,
	 * the number 1e-10 is sent to linadd in lineset1 as follows:*/
	/*linadd( 1.e-10 , 1 , "Unit" , 'i' );*/
	i = cdLine( "Unit" , 1 , &rate , &absint );
	ASSERT( i> 0 );

	/* this is now the linear vol, rel to inner radius */
	VolComputed = LineSave.lines[i].SumLine(0) /  1e-10;

	/* spherical or plane parallel case? */
	if( radius.Radius/radius.rinner > 1.0001 )
	{
		/* spherical case, 
		 * geometry.iEmissPower is usually 2,
		 * and can be reset to 1 (long slit) or 0 (beam) with 
		 * slit and beam options on aperture */
		VolExpected = geometry.covaper*geometry.FillFac*radius.rinner/(geometry.iEmissPower+1)*
			( powi( radius.Radius/radius.rinner,geometry.iEmissPower+1 ) - 1. );
	}
	else
	{
		/* plane parallel case */
		/* next can be zero for very thin model, depth is always positive */
		VolExpected = geometry.covaper*geometry.FillFac*(radius.depth-DEPTH_OFFSET);
	}

	/* now get the relative difference between computed and expected volumes */
	error = fabs(VolComputed - VolExpected)/SDIV(VolExpected);

	/* we need to ignore this test if filling factor changes with radius, or
	 * cylinder geometry in place */
	if( radius.lgCylnOn || geometry.filpow!=0. )
	{
		error = 0.;
	}

	/* how large is relative error? */
	if( error > 0.001 )
	{
		sprintf( chLine, 
			" W-PrtComment insanity - Line unit integration did not verify \n");
		warnin(chLine);
		fprintf( ioQQQ,
			" PROBLEM PrtComment insanity - Line unit integration did not verify \n");
		fprintf( ioQQQ,
			" expected, derived vols were %g %g \n",
			VolExpected , VolComputed );
		fprintf( ioQQQ,
			" relative difference is %g, ratio is %g.\n",error,VolComputed/VolExpected);
		TotalInsanity();
	}

	/* next do same thing for fake continuum point propagated in highest energy cell, plus 1 
	 *  = 
	 * the variable rfield.ConEmitLocal[rfield.nflux]
	 * are set to 
	 * the number 1.e-10f * Dilution in RT_diffuse.  this is the outward
	 * local emissivity, per unit vol.  It is then added to the outward beams
	 * by the rest of the code, and then checked here.
	 *
	 * insanity will be detected if diffuse emission is thrown into the outward beam
	 * in MadeDiffuse.  this happens if the range of ionization encompasses the full
	 * continuum array, up to nflux.  */
	ConComputed = rfield.ConInterOut[rfield.nflux]/ 1e-10;
	/* correct for fraction that went out, as set in ZoneStart,
	 * this should now be the volume of the emitting region */
	ConComputed /= ( (1. + geometry.covrt)/2. );

	/* we expect this to add up to the integral of unity over r^-2 */
	if( radius.Radius/radius.rinner < 1.001 )
	{
		/* plane parallel case, no dilution, use thickness */
		ConExpected = (radius.depth-DEPTH_OFFSET)*geometry.FillFac;
	}
	else
	{
		/* spherical case */
		ConExpected = radius.rinner*geometry.FillFac * (1. - radius.rinner/radius.Radius );
	}
	/* this is impossible */
	ASSERT( ConExpected > 0. );

	/* now get the relative difference between computed and expected volumes */
	error = fabs(ConComputed - ConExpected)/ConExpected;

	/* we need to ignore this test if filling factor changes with radius, or
	 * cylinder geometry in place */
	if( radius.lgCylnOn || geometry.filpow!=0. )
	{
		error = 0.;
	}

	/* \todo 2 - These "volumes" seem to be too small by a factor of two.
	 * rfield.ConInterOut[rfield.nflux] (hence ConComputed) and ConExpected  
	 * should be greater by a factor of 2 if comparison is really of "volume"
	 * of 1/cc pencil. */

	/* how large is relative error? */
	if( error > 0.001 )
	{
		sprintf( chLine, 
			" W-PrtComment insanity - Continuum unit integration did not verify \n");
		warnin(chLine);
		fprintf( ioQQQ," PROBLEM PrtComment insanity - Continuum unit integration did not verify \n");
		fprintf( ioQQQ," exact vol= %g, derived vol= %g relative difference is %g \n",
			ConExpected , ConComputed ,error);
		fprintf( ioQQQ," ConInterOut= %g,  \n",
			rfield.ConInterOut[rfield.nflux]);
		TotalInsanity();
	}

	if( called.lgTalk )
	{
		/* print the title of the calculation */
		fprintf( ioQQQ, "   %s\n", input.chTitle.c_str() );
		/* say why the calculation stopped, and indicate the geometry*/
		cdReasonGeo(ioQQQ);
		/* print all warnings */
		cdWarnings(ioQQQ);
		/* all cautions */
		cdCautions(ioQQQ);
		/* surprises, beginning with a ! */
		cdSurprises(ioQQQ);
		/* notes about the calculations */
		cdNotes(ioQQQ);
	}

	/* option to print warnings on special io */
	if( lgPrnErr )
	{
		cdWarnings(ioPrnErr);
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " PrtComment returns.\n" );
	}
	return;
}


/*chkCaHeps check whether CaII K and H epsilon overlap */
STATIC void chkCaHeps(double *totwid)
{
	DEBUG_ENTRY( "chkCaHeps()" );

	*totwid = 0.;

	if( !atmdat.lgdBaseSourceExists[ipCALCIUM][1] )
	{
		return;
	}

	/* pumping of CaH overlapping with H epsilon, 6-2 of H */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_local + 
		iso_sp[ipH_LIKE][ipHYDROGEN].nCollapsed_local >= 6 )
	{
		/* this is 6P */
		long ip6p = iso_sp[ipH_LIKE][ipHYDROGEN].QN2Index(6, 1, 2);

		long id_Ca2 = -1;
		if( (id_Ca2 = atmdat.ipSpecIon[ipCALCIUM][1]) < 0 )
		{
			fprintf(ioQQQ,"PROBLEM: Ca II, the species defined by nelem = %i and ion = %i could not be found.\n",ipCALCIUM,2);
			cdEXIT(EXIT_FAILURE);
		}

		static TransitionList::iterator it;
		static bool lgRunOnce = true;
		if( lgRunOnce )
		{
			for( it = dBaseTrans[id_Ca2].begin(); it != dBaseTrans[id_Ca2].end(); ++it)
			{
				if( it->ipLo()+1 == 1 && it->ipHi()+1 == 5)
				{
					lgRunOnce = false;
					break;
				}
			}
		}

		realnum ca2_3969_TauIn = it->Emis().TauIn();

		if( ca2_3969_TauIn > 0. &&
			iso_sp[ipH_LIKE][ipHYDROGEN].trans(ip6p,ipH2s).Emis().TauIn() >   0. )
		{
			/* casts to double here are to prevent FPE */
			double conca = sqrt(6.1e-5*ca2_3969_TauIn);
			double conalog = log((double)ca2_3969_TauIn);
			conalog = sqrt(MAX2(1., conalog));
			conca = MAX2(conalog,conca);

			conalog = log((double)iso_sp[ipH_LIKE][ipHYDROGEN].trans(ip6p,ipH2s).Emis().TauIn());
			conalog = sqrt(MAX2(1.,conalog));
			double conhe = sqrt(1.7e-6*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ip6p,ipH2s).Emis().TauIn());
			conhe = MAX2(conalog, conhe);

			*totwid = 10.*conhe + 1.6*conca;
		}
	}
	return;
}

/*prt_smooth_predictions check whether fluctuations in any predicted quantities occurred */
STATIC void prt_smooth_predictions(void)
{
	long int i,
		nzone_oscillation,
		nzone_ion_jump,
		nzone_den_jump,
		nelem,
		ion;
	double BigOscillation ,
		big_ion_jump,
		big_jump,
		rel,
		big,
		bigm;

	char chLine[INPUT_LINE_LENGTH];

	DEBUG_ENTRY( "prt_smooth_predictions()" );

	/* check for ionization oscillations or fluctuations and or jumps */
	nzone_oscillation = 0;
	nzone_ion_jump = 0;

	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			for( ion=0; ion<=nelem+1; ++ion) 
			{
				BigOscillation = 0.;
				big_ion_jump = -15.;
				for( i=1; i < nzone-1; i++ )
				{

					/* only do check if all ions are positive */
					if( struc.xIonDense[i-1][nelem][ion]/struc.gas_phase[i-1][nelem]>struc.dr_ionfrac_limit &&
						struc.xIonDense[i  ][nelem][ion]/struc.gas_phase[i  ][nelem]>struc.dr_ionfrac_limit &&
						struc.xIonDense[i+1][nelem][ion]/struc.gas_phase[i+1][nelem]>struc.dr_ionfrac_limit )
					{

						/* this is check for oscillations */
						big = fabs( (struc.xIonDense[i-1][nelem][ion] - struc.xIonDense[i][nelem][ion])/struc.xIonDense[i][nelem][ion] );
						bigm = fabs( (struc.xIonDense[i][nelem][ion]  - struc.xIonDense[i+1][nelem][ion])/struc.xIonDense[i][nelem][ion] );

						/* this is sign of change in ionization, we are looking for change in sign */
						rel = ( (struc.xIonDense[i-1][nelem][ion] - struc.xIonDense[i][nelem][ion]  )/struc.xIonDense[i][nelem][ion])*
							  ( (struc.xIonDense[i][nelem][ion]   - struc.xIonDense[i+1][nelem][ion])/struc.xIonDense[i][nelem][ion] );

						if( rel < 0. && MIN2( bigm , big ) > BigOscillation )
						{
							nzone_oscillation = i;
							BigOscillation = MIN2( bigm , big );
						}

						/* check whether we tripped over an ionization front - a major source
						 * of instability in a complete linearization code like this one */
						/* neg sign picks up only increases in ionization */
						rel = -log10( (struc.xIonDense[i][nelem][ion]/struc.gas_phase[i][nelem]) / 
							(struc.xIonDense[i+1][nelem][ion]/struc.gas_phase[i+1][nelem] ) );
						/* only do significant stages of ionization */
						if( rel > big_ion_jump )
						{
							big_ion_jump = rel;
							nzone_ion_jump = i;
						}
					}
				}
				/* end loop over zones, 
				 * check whether this ion and element underwent fluctuations or jump */

				if( BigOscillation > 0.2 )
				{
					/* only do sanity check if jump detected */
					if( nzone_oscillation < 1 )
					{
						fprintf( ioQQQ, " nzone_oscillation too small bigjump2 check\n" );
						ShowMe();
						cdEXIT(EXIT_FAILURE);
					}

					if( BigOscillation > 3. )
					{
						sprintf( chLine, 
							" C-An ionization oscillation occurred, elem %.2s%2li, by %.0f%%"
							 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
							 elementnames.chElementSym[nelem], ion+1,
							 BigOscillation*100.,
							 struc.xIonDense[nzone_oscillation-1][nelem][ion]/struc.gas_phase[nzone_oscillation-1][nelem],
							 nzone_oscillation,
							 struc.xIonDense[nzone_oscillation][nelem][ion]/struc.gas_phase[nzone_oscillation][nelem],
							 nzone_oscillation+1,
							 struc.xIonDense[nzone_oscillation+1][nelem][ion]/struc.gas_phase[nzone_oscillation+1][nelem],
							 nzone_oscillation+2 ); 
						caunin(chLine);
					}
					else if( BigOscillation > 0.7 )
					{
						sprintf( chLine, 
							"  !An ionization oscillation occurred, elem %.2s%2li, by %.0f%%"
							 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
							 elementnames.chElementSym[nelem], ion+1,
							 BigOscillation*100.,
							 struc.xIonDense[nzone_oscillation-1][nelem][ion]/struc.gas_phase[nzone_oscillation-1][nelem],
							 nzone_oscillation,
							 struc.xIonDense[nzone_oscillation][nelem][ion]/struc.gas_phase[nzone_oscillation][nelem],
							 nzone_oscillation+1,
							 struc.xIonDense[nzone_oscillation+1][nelem][ion]/struc.gas_phase[nzone_oscillation+1][nelem],
							 nzone_oscillation+2 );
						bangin(chLine);
					}
				}

				/* big_ion_jump was a log above, convert to linear quantity */
				/* if no jump occurred then big_ion_jump is small and nzone_ion_jump is 0 */
				big_ion_jump = exp10( big_ion_jump );
				if( big_ion_jump > 1.5 && nzone_ion_jump > 0 )
				{
					if( big_ion_jump > 10. )
					{
						sprintf( chLine, 
							 " C-An ionization jump occurred, elem %.2s%2li, by %.0f%%"
							 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
							 elementnames.chElementSym[nelem], ion+1,
							 big_ion_jump*100.,
							 struc.xIonDense[nzone_ion_jump-1][nelem][ion]/struc.gas_phase[nzone_ion_jump-1][nelem],
							 nzone_ion_jump,
							 struc.xIonDense[nzone_ion_jump][nelem][ion]/struc.gas_phase[nzone_ion_jump][nelem],
							 nzone_ion_jump+1,
							 struc.xIonDense[nzone_ion_jump+1][nelem][ion]/struc.gas_phase[nzone_ion_jump+1][nelem],
							 nzone_ion_jump+2 );
						caunin(chLine);
					}
					else
					{
						sprintf( chLine, 
							 "  !An ionization jump occurred elem %.2s%2li, by %.0f%%"
							 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
							 elementnames.chElementSym[nelem], ion+1,
							 big_ion_jump*100., 
							 struc.xIonDense[nzone_ion_jump-1][nelem][ion]/struc.gas_phase[nzone_ion_jump-1][nelem],
							 nzone_ion_jump,
							 struc.xIonDense[nzone_ion_jump][nelem][ion]/struc.gas_phase[nzone_ion_jump][nelem],
							 nzone_ion_jump+1,
							 struc.xIonDense[nzone_ion_jump+1][nelem][ion]/struc.gas_phase[nzone_ion_jump+1][nelem],
							 nzone_ion_jump+2 );
						bangin(chLine);
					}
				}
			}
		}
	}

	big_jump = -15;
	nzone_den_jump = 0;

	for( i=1; i < nzone-1; i++ )
	{
		/* this first check is on how the total hydrogen density has changed */
		rel = fabs(log10( struc.gas_phase[i][ipHYDROGEN] / 
			struc.gas_phase[i+1][ipHYDROGEN] ) );
		/* only do significant stages of ionization */
		if( rel > big_jump )
		{
			big_jump = rel;
			nzone_den_jump = i;
		}
	}

	/* check how stable density was */
	big_jump = exp10(  big_jump );
	if( big_jump > 1.2 )
	{
		if( big_jump > 3. )
		{
			sprintf( chLine, 
				" C-The H density jumped at by %.0f%%,"
				 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
				 big_jump*100.,
				 struc.gas_phase[nzone_den_jump-1][ipHYDROGEN], 
				 nzone_den_jump,
				 struc.gas_phase[nzone_den_jump][ipHYDROGEN], 
				 nzone_den_jump+1,
				 struc.gas_phase[nzone_den_jump+1][ipHYDROGEN],
				 nzone_den_jump+2 );
			caunin(chLine);
		}
		else
		{
			sprintf( chLine, 
				 "  !An H density jump occurred by %.0f%%"
				 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
				 big_jump*100.,
				 struc.gas_phase[nzone_den_jump-1][ipHYDROGEN],
				 nzone_den_jump,
				 struc.gas_phase[nzone_den_jump][ipHYDROGEN],
				 nzone_den_jump+1,
				 struc.gas_phase[nzone_den_jump+1][ipHYDROGEN],
				 nzone_den_jump+2 );
			bangin(chLine);
		}
	}

	/* now do check on smoothness of radiation pressure */
	big_jump = -15;
	nzone_den_jump = 0;

	/* loop starts on zone 3 since dramatic fall in radiation pressure across first
	 * few zones is normal behavior */
	for( i=3; i < nzone-2; i++ )
	{
		/* this first check is on how the total hydrogen density has changed */
		rel = fabs(log10( SDIV(struc.pres_radiation_lines_curr[i]) / 
			SDIV(0.5*(struc.pres_radiation_lines_curr[i-1]+struc.pres_radiation_lines_curr[i+1])) ) );
		/* only do significant stages of ionization */
		if( rel > big_jump )
		{
			big_jump = rel;
			nzone_den_jump = i;
		}
	}
	/* note that changing log big_jump to linear takes place in next branch */

	/* check how stable radiation pressure was, but only if significant */
	if( pressure.RadBetaMax > 0.01 )
	{
		big_jump = exp10(  big_jump );
		if( big_jump > 1.2 )
		{
			/* only make it a caution is pressure jumped, and we were trying
			* to do a constant pressure model */
			if( big_jump > 3. && strcmp(dense.chDenseLaw,"CPRE") == 0)
			{
				sprintf( chLine, 
					 " C-The radiation pressure jumped by %.0f%%"
					 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
					 big_jump*100.,
					 struc.pres_radiation_lines_curr[nzone_den_jump-1],
					 nzone_den_jump,
					 struc.pres_radiation_lines_curr[nzone_den_jump],
					 nzone_den_jump+1,
					 struc.pres_radiation_lines_curr[nzone_den_jump+1],
					 nzone_den_jump+2 );
				caunin(chLine);
			}
			else
			{
				sprintf( chLine, 
					 "  !The radiation pressure jumped by %.0f%%"
					 " from %.2e (zone %ld) to %.2e (zone %ld) to %.2e (zone %ld)", 
					 big_jump*100.,
					 struc.pres_radiation_lines_curr[nzone_den_jump-1],
					 nzone_den_jump,
					 struc.pres_radiation_lines_curr[nzone_den_jump],
					 nzone_den_jump+1,
					 struc.pres_radiation_lines_curr[nzone_den_jump+1],
					 nzone_den_jump+2 );
				bangin(chLine);
			}
		}
	}

	/* these will be used to check on continuity */
	phycon.BigJumpTe = 0.;
	phycon.BigJumpne = 0.;
	phycon.BigJumpH2 = 0.;
	phycon.BigJumpCO = 0.;

	for( i=1; i < nzone-1; i++ )
	{
		/* check on how much temperature has changed */
		rel = fabs(log10( struc.testr[i] / struc.testr[i+1] ) );
		if( rel > phycon.BigJumpTe )
		{
			phycon.BigJumpTe = (realnum)rel;
		}

		/* check on how much electron density has changed */
		rel = fabs(log10( struc.ednstr[i] / struc.ednstr[i+1] ) );
		if( rel > phycon.BigJumpne )
		{
			phycon.BigJumpne = (realnum)rel;
		}

		/* check on how much H2 density has changed */
		if( (struc.H2_abund[i])>SMALLFLOAT && (struc.H2_abund[i+1]) > SMALLFLOAT 
			/* only do this test if H2 abund is significant */
			&& (struc.H2_abund[i])/struc.gas_phase[i][ipHYDROGEN]>1e-3)
		{
			rel = fabs(log10( (struc.H2_abund[i]) / SDIV(struc.H2_abund[i+1]) ) );
			if( rel > phycon.BigJumpH2 )
			{
				phycon.BigJumpH2 = (realnum)rel;
			}
		}

		int ipCO = findspecies("CO")->index;
		//fprintf(ioQQQ,"PRTCO %d %ld %d\n",ipCO,i,mole_global.num_calc);
		/* check on how much CO density has changed */
		if( ipCO != -1 &&
		    struc.molecules[i][ipCO]>SMALLFLOAT &&
		    struc.molecules[i+1][ipCO]>SMALLFLOAT &&
		    struc.molecules[i][ipCO]/SDIV(struc.gas_phase[i][ipCARBON])>1e-3 )
		{
			rel = fabs(log10( struc.molecules[i][ipCO] / struc.molecules[i+1][ipCO] ) );
			if( rel > phycon.BigJumpCO )
			{
					phycon.BigJumpCO = (realnum)rel;
			}
		}
	}

	/* convert to linear change - subtract 1 to make it the residual difference */
	if( phycon.BigJumpTe > 0. )
		phycon.BigJumpTe = exp10( phycon.BigJumpTe ) - 1.f;

	if( phycon.BigJumpne > 0. )
		phycon.BigJumpne = exp10( phycon.BigJumpne ) - 1.f;

	if( phycon.BigJumpH2 > 0. )
		phycon.BigJumpH2 = exp10( phycon.BigJumpH2 ) - 1.f;

	if( phycon.BigJumpCO > 0. )
		phycon.BigJumpCO = exp10( phycon.BigJumpCO ) - 1.f;
	/*fprintf(ioQQQ,"DEBUG continuity large change %.2e %.2e %.2e %.2e \n",
		phycon.BigJumpTe , phycon.BigJumpne , phycon.BigJumpH2 , phycon.BigJumpCO );*/

	if( phycon.BigJumpTe > 0.3 )
	{
		sprintf( chLine, 
			 " C-The temperature varied by %.1f%% between two zones", 
			 phycon.BigJumpTe*100.);
		caunin(chLine);
	}
	else if( phycon.BigJumpTe > 0.1 )
	{
		sprintf( chLine, 
			 "  !The temperature varied by %.1f%% between two zones", 
			 phycon.BigJumpTe*100.);
		bangin(chLine);
	}

	if( phycon.BigJumpne > 0.3 )
	{
		sprintf( chLine, 
			 " C-The electron density varied by %.1f%% between two zones", 
			 phycon.BigJumpne*100.);
		caunin(chLine);
	}
	else if( phycon.BigJumpne > 0.1 )
	{
		sprintf( chLine, 
			 "  !The electron density varied by %.1f%% between two zones", 
			 phycon.BigJumpne*100.);
		bangin(chLine);
	}

	if( phycon.BigJumpH2 > 0.8 )
	{
		sprintf( chLine, 
			 " C-The H2 density varied by %.1f%% between two zones", 
			 phycon.BigJumpH2*100.);
		caunin(chLine);
	}
	else if( phycon.BigJumpH2 > 0.1 )
	{
		sprintf( chLine, 
			 "  !The H2 density varied by %.1f%% between two zones", 
			 phycon.BigJumpH2*100.);
		bangin(chLine);
	}

	if( phycon.BigJumpCO > 0.8 )
	{
		sprintf( chLine, 
			 " C-The CO density varied by %.1f%% between two zones", 
			 phycon.BigJumpCO*100.);
		caunin(chLine);
	}
	else if( phycon.BigJumpCO > 0.2 )
	{
		sprintf( chLine, 
			 "  !The CO density varied by %.1f%% between two zones", 
			 phycon.BigJumpCO*100.);
		bangin(chLine);
	}

	if( LineSave.lgIsoContSubSignif )
	{
		sprintf( chLine, 
			 "  !Isotropic continuum subtraction significantly"
			 " affects line intensities" );
		bangin(chLine);
	}
	return;
}
