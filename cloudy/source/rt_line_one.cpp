/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_line_escape do line radiative transfer,
 * evaluates escape and destruction probability */
/* NB no wind distinction - that is done where optical depths are incremented
 * with line opacity - rt_line_one_tauinc and in line width for rad pressure */
/*RT_line_fine_opacity do fine opacities for one line */
/*RT_line_pumping pumping by external and locally emitted radiation fields */
#include "cddefines.h"
#include "rfield.h"
#include "opacity.h"
#include "conv.h"
#include "radius.h"
#include "rt_escprob.h"
#include "rt.h"
#include "cosmology.h"
#include "hydrogenic.h"
#include "iso.h"
#include "wind.h"
#include "geometry.h"

/*RT_line_pumping pumping by external and locally emitted radiation fields */
STATIC void RT_line_pumping(
		 /* the em line we will work on  */
		 const TransitionProxy& t ,
		 /* this is option to not include line self shielding across this zone.
		 * this can cause pump to depend on zone thickness, and leads to unstable
		 * feedback in some models with the large H2 molecule, due to Solomon
		 * process depending on zone thickness and level populations. */
		 bool lgShield_this_zone,
		 realnum DopplerWidth)
{
	DEBUG_ENTRY( "RT_line_pumping()" );

	ASSERT( t.ipCont() >= 1 );

	/* pumping by incident and diffuse continua */
	/* option to kill induced processes */
	if( !rfield.lgInducProcess )
	{
		t.Emis().pump() = 0.;
	}
	else if( conv.lgFirstSweepThisZone || t.Emis().iRedisFun() == ipLY_A )
	{
		double EffectiveThickness = radius.drad_x_fillfac;
		if( cosmology.lgDo )
		{
			EffectiveThickness = safe_div( DopplerWidth, wind.dvdr );
			lgShield_this_zone = false;
		}
		double dTau =  (t.Emis().PopOpc() * t.Emis().opacity() / DopplerWidth + 
			opac.opacity_abs[t.ipCont()-1])* EffectiveThickness;

		/* continuum shielding into this function -- now includes
			optional correction for finite zone depth */
		double shield_continuum = RT_continuum_shield_fcn( t, lgShield_this_zone, dTau );

		/* continuum upward pumping rate, A gu/gl abs prob occnum
		* the "no induced" command causes continuum pumping to be set to 0 
		* this includes pumping by diffuse continuum */
		t.Emis().pump() = t.Emis().Aul() * (*t.Hi()).g() / (*t.Lo()).g() * shield_continuum *(
			rfield.OccNumbIncidCont[t.ipCont()-1] + rfield.OccNumbContEmitOut[t.ipCont()-1] );

		if( 0 && t.chLabel() == "H  1 1215.67A" )
			fprintf(ioQQQ,
				"LINE PUMPING:\t label= \"%s\"\t %g\t %g\t %g\n",
				t.chLabel().c_str(),
				shield_continuum,
				rfield.OccNumbIncidCont[t.ipCont()-1] + rfield.OccNumbContEmitOut[t.ipCont()-1],
				t.Emis().pump()
				);
		/* 
		 * This is an option to account for intrinsic absorption or emission line the lyman 
		 * lines.  This is done without changing the incident coarse continuum.  
		 * Check if continuum pumping of H Lyman lines to be multiplied by scale factor
		 * hydro.xLymanPumpingScaleFactor is set with atom h-like Lyman pumping scale command 
		 * Lya pump rate is always updated - this is simplest way to finese intrinsic absorption
		 * or emission 
		 */
		if ( t.ipLo() == 0 && t.systemIs(iso_sp[ipH_LIKE][ipHYDROGEN].tr) )
		{
			t.Emis().pump() *=  hydro.xLymanPumpingScaleFactor;
		}
		/* NB NB line pumping by local diffuse emission is not included. */
	}

	return;
}

/*RT_line_escape do line radiative transfer escape and destruction probabilities 
 * this routine sets */
STATIC void RT_line_escape(
	   /* the em line we will work on  */
	   const TransitionProxy &t,
	   /* Stark escape probability to be added to Pesc */
	   realnum pestrk,
	   realnum DopplerWidth,
	   bool lgGoodTau)
{
	DestType nRedis;

	DEBUG_ENTRY( "RT_line_escape()" );

	/* a runaway maser */
	if( t.Emis().TauIn() < -30. )
	{
		fprintf( ioQQQ, "PROBLEM RT_line_escape called with large negative "
			"optical depth, zone %.2f, aborting.\n",
			fnzone );
		DumpLine(t);
		throw cloudy_abort("large negative optical depth");
	}

	if( cosmology.lgDo )
	{
		/* Sobolev escape */
		if( conv.lgFirstSweepThisZone && lgGoodTau )
		{
			realnum tau_Sobolev =  t.Emis().TauTot();

			if( tau_Sobolev < 1E-5 )
			{
				t.Emis().Pesc() = 1.;
			}
			else
			{
				t.Emis().Pesc() = ( 1.f - exp( -1.f * tau_Sobolev ) )/ tau_Sobolev;
			}

			/* inward escaping fraction */
			t.Emis().FracInwd() = rt.fracin;
		}
		fixit("is this correct?");
		nRedis.t = DestType::ipDEST_K2;
	}
	else if (0) // LVG escape from Castor, Radiation Hydrodynamics
	{
		if( conv.lgFirstSweepThisZone && lgGoodTau )
		{
			long int ipLineCenter = t.Emis().ipFine() + rfield.ipFineConVelShift;

			double OpacityEffective;
			// Don't use fine opacity -- need to make sure a fully updated
			// value is used, or use other means to handle overlap
			if( 0 && 
				 t.Emis().ipFine()>=0 && ipLineCenter>0 && 
				 ipLineCenter<rfield.nfine && rfield.lgOpacityFine )
			{
				OpacityEffective = rfield.fine_opac_zone[ipLineCenter];
			}
			else
			{
				OpacityEffective = t.Emis().VoigtLineCen() * t.Emis().PopOpc() *
							t.Emis().opacity() / DopplerWidth;
			}

			// Castor, Radiation Hydrodynamics, p127
			double exprate, sigma;
			if ( !cosmology.lgDo )
			{
				exprate = radius.depth / wind.windv;
				sigma = exprate * wind.dvdr  - 1.0; // (6.103)
			}
			else
			{
				exprate = 1.0/GetHubbleFactor(cosmology.redshift_current);
				sigma = 0.0;
			}
			double tau = exprate * geometry.FillFac * OpacityEffective * 
				DopplerWidth; // (6.106)			
			t.Emis().Pesc() = (realnum) RT_EscLVG(tau, sigma);
			t.Emis().FracInwd() = 0.5;
		}
		fixit("is this correct?");
		nRedis.t = DestType::ipDEST_K2;
	}
	/* static solution - which type of line will determine
	 * which redistribution function */
	/* iRedisFun() == 1 - alpha resonance line, partial redistribution,
	 * ipPRD == 1 */
	else if( t.Emis().iRedisFun() == ipPRD )
	{
		/* incomplete redistribution with wings */
		if( conv.lgFirstSweepThisZone && lgGoodTau )
		{
			t.Emis().Pesc() = (realnum)esc_PRD( t.Emis().TauIn(), t.Emis().TauTot(), t.Emis().damp() );

			/* >>chng 03 jun 07, do not clobber esp prob when line is masing -
			* this had effect of preventing total escape prob from getting larger than 1 */
			if( pestrk > 0.f && t.Emis().Pesc() < 1.f )
				t.Emis().Pesc() = min( 1.f, t.Emis().Pesc() + pestrk );

			/* inward escaping fraction */
			t.Emis().FracInwd() = rt.fracin;
		}
		nRedis.t = DestType::ipDEST_INCOM;
	}

	/* complete redistribution without wings - t.ipLnRedis is ipCRD == -1 */
	else if( t.Emis().iRedisFun() == ipCRD )
	{
		if( conv.lgFirstSweepThisZone && lgGoodTau )
		{
			/* >>chng 01 mar -6, escsub will call any of several esc prob routines,
			* depending of how core is set.  We always want core-only for this option,
			* so call  esca0k2(tau) directly */
			t.Emis().Pesc() = (realnum)esc_CRDcore( t.Emis().TauIn(), t.Emis().TauTot() );

			if( pestrk > 0.f && t.Emis().Pesc() < 1.f )
				t.Emis().Pesc() = min( 1.f, t.Emis().Pesc() + pestrk );

			/* inward escaping fraction */
			t.Emis().FracInwd() = rt.fracin;
		}
		nRedis.t = DestType::ipDEST_K2;
	}

	/* CRD with wings, = 2 */
	else if( t.Emis().iRedisFun() == ipCRDW )
	{
		/* complete redistribution with damping wings */
		if( conv.lgFirstSweepThisZone && lgGoodTau )
		{
			t.Emis().Pesc() = (realnum)esc_CRDwing( t.Emis().TauIn(), t.Emis().TauTot(), t.Emis().damp() );

			if( pestrk > 0.f && t.Emis().Pesc() < 1.f )
				t.Emis().Pesc() = min( 1.f, t.Emis().Pesc() + pestrk );

			/* inward escaping fraction */
			t.Emis().FracInwd() = rt.fracin;
		}
		nRedis.t = DestType::ipDEST_K2;
	}

	/* Lya is special case */
	else if( t.Emis().iRedisFun() == ipLY_A )
	{
		/* incomplete redistribution with wings, for special case of Lya
		* uses fits to Hummer & Kunasz numerical results 
		* this routine is different because escape and dest probs
		* are evaluated together, so no test of lgDoEsc */
		if( lgGoodTau )
		{
			double dest , esin;

			/* this will always evaluate escape prob, no matter what lgDoEsc is.
			* The destruction prob comes back as dest */
			t.Emis().Pesc() = (realnum)RTesc_lya( &esin, &dest, t.Emis().PopOpc(), t, DopplerWidth );

			if( pestrk > 0.f && t.Emis().Pesc() < 1.f )
				t.Emis().Pesc() = min( 1.f, t.Emis().Pesc() + pestrk );

			/*  this is fraction of line which is inward escaping */
			t.Emis().FracInwd() = rt.fracin;

			nRedis.t = DestType::ipDEST_LYA;
			nRedis.dest = dest;
		}
	}
	else
	{
		fprintf( ioQQQ, " RT_line_escape called with impossible redistribution function %d\n",
				t.Emis().iRedisFun());
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	if( lgGoodTau && t.Emis().opacity() > 0. )
	{
		// Now include electron escape
		RT_DestProb(t,
						/* line width in velocity units */
						DopplerWidth,
						/* redistribution function */
						nRedis);
	}
	return;
}

/*RT_line_fine_opacity do fine opacities for one line */
STATIC void RT_line_fine_opacity(
	/* the em line we will work on  */
	const TransitionProxy& t,
	realnum DopplerWidth)
{
	DEBUG_ENTRY( "RT_line_fine_opacity()" );

	/* this is line center frequency, including bulk motion of gas */
	long int ipLineCenter = t.Emis().ipFine() + rfield.ipFineConVelShift;

	/* define fine opacity fine grid fine mesh */
	/* rfield.lgOpacityFine flag set false with no fine opacities command */
	/* opacities can be negative if masers are allowed */
	ASSERT (conv.lgLastSweepThisZone);
	if( ipLineCenter < 0 || abs(t.Emis().PopOpc()) < SMALLFLOAT ||
		ipLineCenter>rfield.nfine || !rfield.lgOpacityFine )
	{
		return;
	}

	/* number of fine opacity cells corresponding to one doppler width for current
	 * temperature, velocity field, and nuclear mass,
	 * rfield.fine_opac_velocity_width is width per cell, cm/s */
	realnum cells_wide_1x = DopplerWidth/rfield.fine_opac_velocity_width;

	/* line center opacity - type realnum since will add to fine opacity array,
	 * which is realnum */
	realnum opac_line = (realnum)t.Emis().PopOpc() * t.Emis().opacity() / DopplerWidth;

	// this is effective optical depth to this point. Do not do line if 
	// this product is less than SMALLFLOAT
	// negative optical depth due to maser effect are allowed.
	double dTauEffec = opac_line*radius.depth_x_fillfac;
	if( abs(dTauEffec) < SMALLFLOAT )
		return;

	/* core width of optically thick line, do 4x with exponential Doppler core,
	 * must be at least one cell, but profile is symmetric */
	const bool doDamp = dTauEffec*t.Emis().damp()/9. > 0.1;

	long int nCells_damp;
	/* include damping wings if optical depth is large */
	if( doDamp )
	{
		// find number of cells to extend the damping wings - cells_wide_1x is one dop width
		// tests with th85orion, stopping at half -h2 point, showed that 0.01 was
		// needed for outer edge, given the definition of dTauEffec
		realnum x = realnum(sqrt( dTauEffec * t.Emis().damp()*100./SQRTPI ) * cells_wide_1x);
		// test on size of x, which can exceed
		// limits of long in extreme optical depths */
		long max_cells = max(rfield.nfine-ipLineCenter, ipLineCenter+1);
		if( x < realnum(max_cells) )
		{
			nCells_damp = long(x);
		}
		else
		{
			/* x was too big, just set to extreme range, which is
			* number of cells to farthest boundary */
			nCells_damp = max_cells;
		}
	}
	else
	{
		// Doppler core only
		nCells_damp = long(cells_wide_1x*4.f + 1.5f);
	}

	/* want line to be at least one cell wide */
	nCells_damp = max( 1, nCells_damp );

	static vector<realnum> xprofile, profile;
	xprofile.resize(nCells_damp);
	profile.resize(nCells_damp);

	realnum dcell = 1.f/cells_wide_1x;
	for( long int i=0; i < nCells_damp; ++i )
	{
		/* distance from line center in units of doppler width */
		xprofile[i] = i*dcell;
	}
	
	VoigtH(t.Emis().damp(), &xprofile[0], &profile[0], nCells_damp);

	enum { DEBUG_LOC = false };
	if( DEBUG_LOC && t.chLabel() == "Fe17      17.0510A" )
	{
		printf( "line: '%s'\t damp_X_vel= %.4e\t damp= %.4e\t lifetime= %.4e\t phi0= %.4e\n",
			t.chLabel().c_str(),
			t.Emis().dampXvel(),
			t.Emis().damp(),
			(*( t.Hi() )).lifetime(),
			profile[0] );
		exit( 0 );
	}

	// Merging loops with abs() increases dynamics_wind run time by x2 for gcc 4.8.3
	long ilo = max(ipLineCenter-nCells_damp+1,0), ihi = min(ipLineCenter+nCells_damp,rfield.nfine);
	for( long i=ilo; i < ipLineCenter; ++i )
	{
		rfield.fine_opac_zone[i] += profile[ipLineCenter-i]*opac_line;
	}
	for( long i=ipLineCenter; i < ihi; ++i )
	{
		rfield.fine_opac_zone[i] += profile[i-ipLineCenter]*opac_line;
	}

	t.Emis().VoigtLineCen() = profile[0];

	return;
}

/*RT_line_one do rt for emission line structure - calls RT_line_escape or RT_line_wind */
void RT_line_one_escape(
	/* the em line we will work on  */
	const TransitionProxy &t,
	/* this is option to not include line self shielding across this zone.
	 * this can cause pump to depend on zone thickness, and leads to unstable
	 * feedback in some models with the large H2 molecule, due to Solomon
	 * process depending on zone thickness and level populations. */
	bool lgShield_this_zone,
	/* Stark escape probability to be added to Pesc */
	realnum pestrk,
	realnum DopplerWidth )
{
	DEBUG_ENTRY( "RT_line_one_escape()" );

	// do nothing is population and this is not the very first call 
	// skip line transfer if requested with 'no line transfer' command, but never skip Lya
	if( !rfield.lgDoLineTrans && (t.Emis().iRedisFun() != ipLY_A) )
	{
		return;
	}

	/* line damping constant at current temperature  */
	t.Emis().damp() = t.Emis().dampXvel() / DopplerWidth;
	ASSERT( t.Emis().damp() > 0. );

	// do not evaluate if no population 
	if( (*t.Lo()).Pop()<=SMALLFLOAT )
	{
		/* zero population, return after setting everything with side effects */
		t.Emis().Pesc() = 1.f;

		/* inward escaping fraction */
		t.Emis().FracInwd() = 0.5;

		/* pumping rate */
		t.Emis().pump() = 0.;

		/* destruction probability */
		t.Emis().Pdest() = 0.;
		t.Emis().Pelec_esc() = 0.;
	}
	else if( t.EnergyErg() / EN1RYD <= rfield.plsfrq )
	{
		// transition is below plasma frequency - photons not emitted
		t.Emis().Pesc() = SMALLFLOAT;
		t.Emis().Pdest() = SMALLFLOAT;
		t.Emis().Pelec_esc() = SMALLFLOAT;
		t.Emis().pump() = SMALLFLOAT;
	}
	else
	{
		/* this checks if we have overrun the optical depth scale,
		 * in which case the inward optical depth is greater than the
		 * previous iteration's total optical depth.
		 * We do not reevaluate escape probabilities if the optical depth
		 * scale has been overrun due to huge bogus change in solution
		 * that would result */
		bool lgGoodTau = lgTauGood( t );
		RT_line_escape( t, pestrk, DopplerWidth , lgGoodTau);
		RT_line_pumping( t , lgShield_this_zone , DopplerWidth );	
	}

	return;
}

/*RT_line_one do rt for emission line structure - calls RT_line_escape or RT_line_wind */
void RT_line_one_fine(
	/* the em line we will work on  */
	const TransitionProxy &t,
	/* this is option to not include line self shielding across this zone.
	 * this can cause pump to depend on zone thickness, and leads to unstable
	 * feedback in some models with the large H2 molecule, due to Solomon
	 * process depending on zone thickness and level populations. */
	bool /* lgShield_this_zone */,
	/* Stark escape probability to be added to Pesc */
	realnum /* pestrk */ ,
	realnum DopplerWidth )
{
	DEBUG_ENTRY( "RT_line_one_fine()" );

	// do nothing is population and this is not the very first call 
	// skip line transfer if requested with 'no line transfer' command, but never skip Lya
	if( !rfield.lgDoLineTrans && (t.Emis().iRedisFun() != ipLY_A) )
	{
		return;
	}

	/* line damping constant at current temperature  */
	t.Emis().damp() = t.Emis().dampXvel() / DopplerWidth;
	ASSERT( t.Emis().damp() > 0. );

	/* option to keep track of population values during calls,
	 * print out data to make histogram */
	enum {DEBUG_LOC=false};
	if( DEBUG_LOC )
	{
		static long int nTau[100];
		long n;

		if( nzone==0 )
		{
			for(n=0; n<100; ++n )
				nTau[n] = 0;
		}
		if( (*t.Lo()).Pop()<=SMALLFLOAT )
			n = 0;
		else
			n = (long)log10( (*t.Lo()).Pop() )+37;
		n = MIN2( n , 99 );
		n = MAX2( n , 0 );
		++nTau[n];
		if( nzone > 183 )
		{
			for(n=0; n<100; ++n )
				fprintf(ioQQQ,"%li\t%li\n", n , nTau[n] );
			cdEXIT(EXIT_SUCCESS);
		}
	}
	
	if( (*t.Lo()).Pop() > SMALLFLOAT && t.EnergyErg() / EN1RYD > rfield.plsfrq )
	{
		// transition is below plasma frequency - photons not emitted
		// the last sweep through this zone is to do the fine opacities
		// the populations are not updated on the last sweep so the 
		// line transfer details also should not be updated
		RT_line_fine_opacity( t , DopplerWidth );
	}

	return;
}
