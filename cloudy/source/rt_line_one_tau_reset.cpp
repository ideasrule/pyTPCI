/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_line_one_tau_reset computes average of old and new optical depths for new scale at end of iter,
 * called by update, also FeIILevelPops */
#include "cddefines.h"
#include "opacity.h"
#include "geometry.h"
#include "transition.h"
#include "prt.h"
#include "rt.h"
#include "cosmology.h"

void RT_line_one_tau_reset(const TransitionProxy& t)
{
	DEBUG_ENTRY( "RT_line_one_tau_reset()" );

	if( cosmology.lgDo )
		return;

	if( t.ipCont() <= 0 )
	{
		return;
	}

	/* option to print masing lines, set with print maser */
	if( prt.lgPrtMaser && 
	    ( t.Emis().TauTot() < -0.01 || t.Emis().TauIn() < -0.01 ) )
	{
		fprintf( ioQQQ, " Masing line:%10.10s t(in, out)=%10.2e%10.2e\n", 
		  chLineLbl(t).c_str(), t.Emis().TauIn(), t.Emis().TauTot() );
	}

	realnum TauNext;

	/* lgStatic is false by default, set true with command sphere SPHERE STATIC */
	if( geometry.lgSphere && geometry.lgStatic )
	{
		/* static sphere, both sides interact */
		TauNext = t.Emis().TauIn();
	}
	else
	{
		/* end of iteration, the next estimate of total optical depth 
		 * is now the inward optical depth - set T(ipLnTauTot) to T(1)
		 * DoubleTau normally 1, set to 2 by DoubleTau command in order
		 * to simulate two-sided photoionization */
		TauNext = rt.DoubleTau*t.Emis().TauIn();
	}

	if( 0 && iteration > 1 )
	{
		fprintf(ioQQQ,"DEBUG %s ",chLineLbl(t).c_str());
		fprintf(ioQQQ," %13.6e %13.6e %13.6e\n",t.Emis().TauTot(),TauNext,
				  (t.Emis().TauTot()-TauNext)/(t.Emis().TauTot()+TauNext));
	}

	/* iteration is 1 when this routine is called at end of first iteration
	 * estimate is bad when starting first iteration, so don't track that one */
	/* don't track the 2nd iteration either since that one may also be quite poor */
	if( iteration > tracker::PREV_ITER )
		TauNext = t.Emis().TauTrack().next_val( t.Emis().TauTot(), TauNext );
	
	if( geometry.lgSphere && geometry.lgStatic )
		t.Emis().TauIn() = TauNext/2.f;
	else
		t.Emis().TauIn() = opac.taumin;

	t.Emis().TauInSpecific() = opac.taumin;

	t.Emis().TauTot() = TauNext;

	/* this is escape prob */
	t.Emis().Pesc() = 0.5f*(1.f + 1.f/MAX2(1.f,t.Emis().TauTot()));

	/* this is fraction inward */
	t.Emis().FracInwd() = MIN2(1.f,1.5f-t.Emis().Pesc());

	/*  destruction probability */
	t.Emis().Pdest() = 0.;
	t.Emis().Pelec_esc() = 0.;

	/* optical depth to the continuum source */
	t.Emis().TauCon() = opac.taumin;

	/* >>chng 01 sep 01, zero out some pops and energies */
	(*t.Lo()).Pop() = 0.;
	/* >>chng 97 jul 21, added following zero
	 * population of upper level */
	(*t.Hi()).Pop() = 0.;
	/* population of lower level with correction for stim emission */
	t.Emis().PopOpc() = 0.;
	t.Emis().mult_opac() = 0.;
	/* following two heat exchange excitation, deexcitation */
	t.Coll().cool() = 0.;
	t.Coll().heat() = 0.;
	/* intensity of line */
	t.Emis().xIntensity() = 0.;
	t.Emis().xObsIntensity() = 0.;
	return;
}
