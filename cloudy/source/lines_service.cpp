/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*GetGF convert Einstein A into oscillator strength */
/*abscf convert gf into absorption coefficient */
/*RefIndex calculates the index of refraction of air using the line energy in wavenumbers,
 * used to convert vacuum wavelengths to air wavelengths. */
/*eina convert a gf into an Einstein A */
/*WavlenErrorGet - find difference between two wavelengths */
/*linadd enter lines into the line storage array, called once per zone */
/*lindst add local line intensity to line luminosity stack */
/*PntForLine generate pointer for forbidden line */
/*totlin sum total intensity of cooling, recombination, or intensity lines */
/*FndLineHt search through line heat arrays to find the strongest heat source */
/*ConvRate2CS convert down coll rate back into electron cs in case other parts of code need this for reference */
#include "cddefines.h"
#include "lines_service.h"
#include "dense.h"
#include "geometry.h"
#include "ipoint.h"
#include "lines.h"
#include "trace.h"
#include "opacity.h"
#include "radius.h"
#include "rfield.h"
#include "rt.h"
#include "prt.h"
#include "taulines.h"
#include "save.h"

void LineStackCreate()
{
	DEBUG_ENTRY( "LineStackCreate()" );

	// set up emission line intensity stack
	/* there are three types of calls to lines()
	 * ipass = -1, first call, only count number of lines
	 * ipass =  0, second pass, save labels and wavelengths
	 * ipass =  1, integrate intensity*/
	LineSave.ipass = -1;
	lines();
	ASSERT( LineSave.nsum > 0 );

	/* in a grid or MPI run this may not be the first time here,
	 * return old memory and grab new appropriate for this size,
	 * since number of lines to be stored can change */
	LineSave.clear();

	/* this is the large main line array */
	LineSave.resize(LineSave.nsum);

	/* this is only done on first iteration since will integrate over time */
	for( long i=0; i < LineSave.nsum; i++ )
	{
		LineSave.lines[i].SumLineZero();
		LineSave.lines[i].SumLineZeroAccum();
	}

	/* there are three calls to lines()
	 * first call with ipass = -1 was above, only counted number
	 * of lines and allocated space.  This is second call and will do
	 * one-time creation of line labels */
	LineSave.ipass = 0;
	lines();
	/* has to be positive */
	ASSERT( LineSave.nsum > 0);
	/* in the future calls to lines will result in integrations */
	LineSave.ipass = 1;

	if( trace.lgTrace )
		fprintf( ioQQQ, "%7ld lines printed in main line array\n",
		  LineSave.nsum );

	if( false )
	{
		//
		// Go over the map of fine continuum indices just generated
		// and quantify the amount of line overlap present.
		// Currently (Jun 03, 2022) there exist bins with up to 14
		// lines centered on them -- these are mainly due to hyperfine
		// molecular lines.
		//

		// Create an ordered map of the number of lines per bin (key)
		// and store a list (value) of all fine continuum indices that
		// have that many lines overlap.
		//
		map<size_t, vector<long> > overlap;
		for( auto it = rfield.fine_lstack.begin();
			it != rfield.fine_lstack.end(); ++it )
		{
			if( overlap.find( it->second.size() ) == overlap.end() )
			{
				overlap.insert( make_pair( it->second.size(), vector<long>() ) );
			}
			overlap[ it->second.size() ].emplace_back( it->first );
		}

		fprintf( ioQQQ, "Overlap in fine continuum:\n" );
		fprintf( ioQQQ, "==========================\n" );

		const size_t nlines_per_bin = 10;

		// loop over line overlap (number of lines in a fine cont bin)
		//
		for( auto it = overlap.begin(); it != overlap.end(); ++it )
		{
			if( it->first > nlines_per_bin )
			{
				sort( it->second.begin(), it->second.end() );

				// loop over fine continuum indices
				//
				for( auto itf = it->second.begin();
					itf != it->second.end(); ++itf )
				{
					TransitionProxy tr;

					// loop over line stack indices
					//
					for( auto itl = rfield.fine_lstack[ *itf ].begin();
						itl != rfield.fine_lstack[ *itf ].end(); ++itl )
					{
						tr = LineSave.lines[*itl].getTransition();
						fprintf( ioQQQ, "%d\t%s\n",
								*itl, tr.chLabel().c_str() );
					}
					fprintf( ioQQQ, "----------------------------\n" );
				}
			}

			fprintf( ioQQQ, "%d lines per bin: %d bins\n",
					int( it->first ), int( it->second.size() ) );
			fprintf( ioQQQ, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
		}
	}
}

/*eina convert a gf into an Einstein A */
double eina(double gf,
	    double enercm, 
	    double gup)
{
	double eina_v;

	DEBUG_ENTRY( "eina()" );

	/* derive the transition prob, given the following
	 * call to function is gf, energy in cm^-1, g_up
	 * gf is product of g and oscillator strength
	 * eina = ( gf / 1.499e-8 ) / (wl/1e4)**2 / gup  */
	eina_v = (gf/gup)*TRANS_PROB_CONST*POW2(enercm);
	return( eina_v );
}

/*GetGF convert Einstein A into oscillator strength */
double GetGF(double trans_prob, 
	     double enercm, 
	     double gup)
{
	double GetGF_v;

	DEBUG_ENTRY( "GetGF()" );

	ASSERT( enercm > 0. );
	ASSERT( trans_prob > 0. );
	ASSERT( gup > 0.);

	/* derive the transition prob, given the following
	 * call to function is gf, energy in cm^-1, g_up
	 * gf is product of g and oscillator strength
	 * trans_prob = ( GetGF/gup) / 1.499e-8 / ( 1e4/enercm )**2 */
	GetGF_v = trans_prob*gup/TRANS_PROB_CONST/POW2(enercm);
	return( GetGF_v );
}

/* S2Aul convert line strength S into transition probability Aul */
double S2Aul(double S,
	     double waveAng,
	     double gup,
	     const string& transType)
{
	DEBUG_ENTRY( "S2Aul()" );

	static const double BOHR_MAGNETON = ELEM_CHARGE_ESU*H_BAR/2./ELECTRON_MASS/SPEEDLIGHT;
	/** Bohr Magneton, 9.2740096e-21 erg/G */

	static const double EC = pow2(ELEM_CHARGE_ESU)/HPLANCK;
	static const double MC = pow2(BOHR_MAGNETON)/HPLANCK;
	/* common factors */

	if( transType == "E1" )
	{
		/*Convert line strength to Aul for E1 transitions
		 * Aul = 64*Pi^4*e^2*a0^2/3/h*S/gu/WLAng^3  */
		static const double E1Coeff = 64./3.*powi(PI,4)*EC*pow2(BOHR_RADIUS_CM)/pow3(1e-8);
		/* E1Coeff = 2.0261e18 */

		return E1Coeff*S/gup/pow3(waveAng);
	}
	else if( transType == "E2" )
	{
		/*Convert line strength to Aul for E2 transitions
		 * Aul = 64*Pi^6*e^2*a0^4/15/h*S/gu/WLAng^5  */
		static const double E2Coeff = 64./15.*powi(PI,6)*EC*pow4(BOHR_RADIUS_CM)/powi(1e-8,5);
		/* E2Coeff = 1.1199e18 */

		return E2Coeff*S/gup/powi(waveAng,5);
	}
	else if( transType == "E3" )
	{
		/*Convert line strength to Aul for E3 transitions
		 * Aul = 2048*Pi^8*e^2*a0^6/4725/h*S/gu/WLAng^7 */
		static const double E3Coeff = 2048./4725.*powi(PI,8)*EC*powi(BOHR_RADIUS_CM,6)/powi(1e-8,7);
		/* E3Coeff = 3.1444165e17
		 * Atomic Transition Probabilites of Silicon. A Critical Compilation
		 * Kelleher & Podobedova 2006*/

		return E3Coeff*S/gup/powi(waveAng,7);
	}
	else if( transType == "M1" )
	{
		/*Convert line strength to Aul for M1 transitions
		 * Aul = 64*Pi^4*mu_B^2/3/h*S/gu/WLAng^3 */
		static const double M1Coeff = 64./3.*powi(PI,4)*MC/pow3(1e-8);
		/* M1Coeff = 2.697e13 */

		return M1Coeff*S/gup/pow3(waveAng);
	}
	else if( transType == "M2" )
	{
		/*Convert line strength to Aul for M2 transitions
		 * Aul = 64*Pi^6*mu_B^2*a0^2/15/h*S/gu/WLAng^5 */
		static const double M2Coeff = 64./15.*powi(PI,6)*MC*pow2(BOHR_RADIUS_CM)/powi(1e-8,5);
		/* M2Coeff = 1.4909714e13
		 * Tables of Atomic Transition Probabilities for Be and B
		 * Fuhr & Wiese  2010*/

		return M2Coeff*S/gup/powi(waveAng,5);
	}
	else if( transType == "M3" )
	{
		/*Convert line strength to Aul for M3 transitions
		 * Aul = 2048*Pi^8*mu_B^2*a0^4/4725/h*S/gu/WLAng^7 */
		static const double M3Coeff = 2048./4725.*powi(PI,8)*MC*pow4(BOHR_RADIUS_CM)/powi(1e-8,7);
		/* M2Coeff = 4.18610e12
		 * Safronova & Safronova et al 2005*/

		return M3Coeff*S/gup/powi(waveAng,7);
	}
	else
	{
		fprintf( ioQQQ, " Invalid transition type in S2Aul: %s.\n", transType.c_str() );
		cdEXIT(EXIT_FAILURE);
	}
}

/*abscf convert gf into absorption coefficient */
double abscf(double gf, 
	  double enercm, 
	  double gl)
{
	double abscf_v;

	DEBUG_ENTRY( "abscf()" );

	ASSERT(gl > 0. && enercm > 0. && gf >= 0.0 );

	/* derive line absorption coefficient, given the following:
	 * gf, enercm, g_low
	 * gf is product of g and oscillator strength */
	abscf_v = ABSOR_COEFF_CONST * (gf/gl)/enercm;
	return( abscf_v );
}

/* compute wavelength in air or vacuum given hardcoded air wavelengths,
 * option set by parse option PRINT WAVELENGTH VACUUM
 * this allows hardwired air wavelengths (which should not be there
 * in the first place) to be converted air/vacuum automatically
 */
realnum wlAirVac( double wlAir )
{
	DEBUG_ENTRY( "wlAirVac()" );

	// Iterate since EnergyWN depends on wlVac not wlAir but
	// difference should be small
	double RefIndex_v = 1.;
	if( !prt.lgPrintLineAirWavelengths && wlAir > 2000.)
	{
		double wlVacuum = wlAir;
		for( int i=0; i<2; ++i )
		{
			/* WN is wavenumber in microns^-1, WN2 is this squared */
			double WN = 1e4 / wlVacuum;
			double WN2 = WN*WN;

			/* use a formula from
			 *>>refer	air	index refraction	Peck & Reeder 1972, JOSA, 62, 8, 958 */
			RefIndex_v = 1. +
				1e-8 * (8060.51 + 2480990.0 / (132.274 - WN2) + 17455.7 / (39.32957 - WN2));

			wlVacuum = wlAir * RefIndex_v;
		}
	}

	return( (realnum)(wlAir * RefIndex_v) );
}

/*RefIndex calculates the index of refraction of air using the line energy in wavenumbers,
 * by default for STP air, returns index of refraction in vacuum (1) when
 * print line vacuum set
 * used to convert vacuum wavelengths to air wavelengths. */
double RefIndex(double EnergyWN )
{
	DEBUG_ENTRY( "RefIndex()" );

	ASSERT( EnergyWN > 0. );

	double RefIndex_v = 1.0;

	/* only do index of refraction if longward of 2000A */
	if( EnergyWN < 5e4 && prt.lgPrintLineAirWavelengths )
	{
		/* xl is wavenumber in microns^-1, squared */
		double xl = EnergyWN * 1e-4;
		xl *= xl;

		/* use a formula from 
		 *>>refer	air	index refraction	Peck & Reeder 1972, JOSA, 62, 8, 958 */
		RefIndex_v += 1e-8 * (8060.51 + 2480990.0 / (132.274 - xl) + 17455.7 / (39.32957 - xl));
	}

	ASSERT( RefIndex_v >= 1. );
	return( RefIndex_v );
}

/*WavlenErrorGet - given the real wavelength in A for a line
 * routine will find the error expected between the real 
 * wavelength and the wavelength printed in the output, with 6 sig figs,
 * function returns difference between exact and 6 sig fig wl, so 
 * we have found correct line is fabs(d wl) < return */
realnum WavlenErrorGet( realnum wavelength, long sig_figs )
{
	double a;
	realnum errorwave;

	DEBUG_ENTRY( "WavlenErrorGet()" );

	ASSERT( sig_figs <= LineSave.sig_figs_max );

	if( wavelength > 0. )
	{
		/* normal case, positive (non zero) wavelength */
		a = log10(wavelength+FLT_EPSILON);
		a = floor(a);
	}
	else
	{
		/* might be called with wl of zero, this is that case */
		/* errorwave = 5e-6f; */
		a = 0.;
	}

	errorwave = 5_r * (realnum)exp10(a - (double)sig_figs);
	return errorwave;
}

/*linadd enter lines into the line storage array, called once per zone for each line*/
STATIC LinSv* lincom(
  double xEmiss,	/* xEmiss - local emissivity per unit vol, no fill fac */
  double xEmissIsoBkg,	/* xEmissIsoBkg - local emissivity corrected for isotropic backgrounds per unit vol, no fill fac */
  realnum wavelength,	/* realnum wavelength */
  const char *chLab,/* string label for ion */
  // ipnt offset of line in continuum mesh
  long int ipnt, 
  char chInfo,		/* character type of entry for line - given below */
			/* 'c' cooling, 'h' heating, 'i' info only, 'r' recom line, 't' transferred line */
  // *chComment string explaining line 
  const char *chComment,
  // lgAdd says whether we've come in via linadd (true) or lindst (false)
  bool lgAdd,
  const TransitionProxy& tr)
{
	DEBUG_ENTRY( "lincom()" );

	/* main routine to actually enter lines into the line storage array
	 * called at top level within routine lines
	 * called series of times in routine PutLine for lines transferred
	 */

	/* three values, -1 is just counting, 0 if init, 1 for calculation */
	if( LineSave.ipass > 0 )
	{
		/* not first pass, sum lines only
		 * total emission from vol */
		/* LineSave.ipass > 0, integration across simulation, sum lines only 
		 * emissivity, emission per unit vol, for this zone */
		LineSave.lines[LineSave.nsum].SumLineAdd(0,xEmissIsoBkg*radius.dVeffAper);
		/* local emissivity in line */
		/* integrated intensity or luminosity, the emissivity times the volume */
		LineSave.lines[LineSave.nsum].emslinSet(0,xEmiss);

		if (lgAdd)
		{
			if (wavelength > 0 && chInfo == 't' )
			{
				/* no need to increment or set [1] version since this is called with no continuum
				 * index, don't know what to do */
				/* only put informational lines, like "Q(H) 4861", in this stack
				 * many continua have a wavelength of zero and are proper intensities,
				 * it would be wrong to predict their transferred intensity */
				LineSave.lines[LineSave.nsum].emslinThin();
				LineSave.lines[LineSave.nsum].SumLineThin(); 
			}
		}
		else
		{
			if ( ipnt <= rfield.nflux && chInfo == 't' )
			{
				/* emergent_line accounts for destruction by absorption outside
				 * the line-forming region */
				const double saveemis = emergent_line( 
					xEmiss*rt.fracin , xEmiss*(1.-rt.fracin) , ipnt );
				LineSave.lines[LineSave.nsum].emslinSet(1,saveemis);

				const double saveemis_isobkg = emergent_line( 
					xEmissIsoBkg*rt.fracin , xEmissIsoBkg*(1.-rt.fracin) , ipnt );
				LineSave.lines[LineSave.nsum].SumLineAdd(1,saveemis_isobkg*radius.dVeffAper);
			}
		}
	}

	else if( LineSave.ipass == 0 )
	{
		LineSave.init(LineSave.nsum,(char) chInfo,chComment,chLab,lgAdd,wavelength,tr);
		if (!lgAdd)
		{
			// check that line wavelength and continuum index agree to some extent
			// this check cannot be very precise because some lines have 
			// "wavelengths" that are set by common usage rather than the correct
			// wavelength derived from energy and index of refraction of air
			ASSERT( ipnt > 0 );
#		ifndef NDEBUG
			double error = MAX2(0.1*rfield.anu(ipnt-1) , rfield.widflx(ipnt-1) );
			ASSERT( wavelength<=0 ||
					  fabs( rfield.anu(ipnt-1) - RYDLAM / wavelength) < error );
#		endif
		}

		if( tr.associated() )
		{
			auto got = rfield.fine_lstack.find( tr.Emis().ipFine() );
			if( got == rfield.fine_lstack.end() )
			{
				rfield.fine_lstack.insert( make_pair( tr.Emis().ipFine(), vector<int>() ) );
			}
			rfield.fine_lstack[ tr.Emis().ipFine() ].emplace_back( LineSave.nsum );
		}
	}

	/* increment the line counter */
	++LineSave.nsum;

	if (LineSave.ipass == -1)
		return NULL;

	return &LineSave.lines[LineSave.nsum-1];

	/* routine can be called with negative LineSave.ipass, in this case
	 * we are just counting number of lines for current setup */
}

/*linadd enter lines into the line storage array, called once per zone for each line*/
LinSv *linadd(
  double xEmiss,	/* xEmiss - local emissivity per unit vol, no fill fac */
  realnum wavelength,	/* realnum wavelength */
  const char *chLab,/* string label for ion */
  char chInfo,		/* character type of entry for line - given below */
			/* 'c' cooling, 'h' heating, 'i' info only, 'r' recom line, 't' transferred line */
  const char *chComment )
{
	DEBUG_ENTRY( "linadd()" );
	
	// Values added to get common interface with lindst
	const long int ipnt = LONG_MAX;
	
	return lincom( xEmiss, xEmiss, wavelength, chLab, ipnt, chInfo, chComment, true, TransitionProxy() );
}


/*emergent_line find emission from surface of cloud after correcting for
 * extinction due to continuous opacity for inward & outward directed emission */
double emergent_line( 
	/* emiemission in inward direction */
	double emissivity_in , 
	/* emission in outward direction */
	double emissivity_out , 
	/* array index for continuum frequency on fortran scale */
	long int ipCont )
{
	double emergent_in , emergent_out;
	long int i = ipCont-1;

	DEBUG_ENTRY( "emergent_line()" );

	ASSERT( i >= 0 && i < rfield.nflux_with_check-1 );

	/* do nothing if first iteration since we do not know the outward-looking
	 * optical depths.  In version C07.02.00 we assumed an infinite optical
	 * depth in the outward direction, which would be appropriate for a 
	 * HII region on the surface of a molecular cloud.  This converged onto
	 * the correct solution in later iterations, but on the first iteration
	 * this underestimated total emission if the infinite cloud were not
	 * present.  With C07.02.01 we make no assuptions about what is in the
	 * outward direction and simply use the local emission. 
	 * Behavior is unchanged on later iterations */
	if( iteration == 1 )
	{
		/* first iteration - do not know outer optical depths so assume very large 
		 * optical depths */
		emergent_in = emissivity_in*opac.E2TauAbsFace[i];
		emergent_out = emissivity_out;
	}
	else
	{
		if( geometry.lgSphere )
		{
			/* second or later iteration in closed or spherical geometry */
			/* inwardly directed emission must get to central hole then across entire
			 * far side of shell */
			emergent_in = emissivity_in  * opac.E2TauAbsFace[i] *opac.E2TauAbsTotal[i];

			/* E2 is outwardly directed emission to get to outer edge of cloud */
			emergent_out = emissivity_out * opac.E2TauAbsOut[i];
		}
		else
		{
			/* open geometry in second or later iteration, outer optical depths are known 
			 * this is light emitted into the outer direction and backscattered
			 * into the inner */
			double reflected = emissivity_out * opac.albedo[i] * (1.-opac.E2TauAbsOut[i]);
			/* E2 is to get to central hole */
			emergent_in = (emissivity_in + reflected) * opac.E2TauAbsFace[i];
			/* E2 is to get to outer edge */
			emergent_out = (emissivity_out - reflected) * opac.E2TauAbsOut[i];
		}
	}
	/* return the net emission that makes it to the surface */
	return( emergent_in + emergent_out );
}

/* outline_base - calls outline_base_bin after deciding whether to add impulse or resolved line */
void outline_base(double dampXvel, double damp, bool lgTransStackLine, long int ip, double phots, realnum inwd,
						double nonScatteredFraction)
{
	DEBUG_ENTRY( "outline_base()" );

	// Need to consider finite optical depth effects if this option is
	// enabled (see fixit below)
	static const bool DO_PROFILE = false;

	if( !DO_PROFILE || !rfield.lgDoLineTrans || dampXvel == 0.0 )
	{
		outline_base_bin(lgTransStackLine, ip, phots, inwd, nonScatteredFraction);
	}
	else
	{
		ASSERT( damp > 0. );
		double LineWidth = dampXvel/damp;
		LineWidth = MIN2( 0.1 * SPEEDLIGHT, LineWidth );
		double sigma = (LineWidth/SPEEDLIGHT);
		long ip3SigmaRed = ipoint( MAX2( rfield.emm(), rfield.anu(ip) - 3.*sigma*rfield.anu(ip) ) );
		long ip3SigmaBlue = ipoint( MIN2( rfield.egamry(), rfield.anu(ip) + 3.*sigma*rfield.anu(ip) ) );
		ASSERT( ip3SigmaBlue >= ip3SigmaRed );
		long numBins = ip3SigmaBlue - ip3SigmaRed + 1;

		if( numBins < 3 )
		{
			outline_base_bin(lgTransStackLine, ip, phots, inwd, nonScatteredFraction);
		}
		else
		{
			valarray<realnum> x(numBins);
			valarray<realnum> profile(numBins);

			for( long ipBin=ip3SigmaRed; ipBin<=ip3SigmaBlue; ipBin++ )
				x[ipBin] = (rfield.anu(ip) - rfield.anu(ipBin))/rfield.anu(ip)/sigma;
			fixit("Escape from zone will be dominated by line wings at"
					" finite optical depth, so Voigt profile may be incorrect");
			// this profile must have unit normalization to conserve energy -> use U(a,v)
			VoigtU(damp,get_ptr(x),get_ptr(profile),numBins);

			for( long ipBin=ip3SigmaRed; ipBin<=ip3SigmaBlue; ipBin++ )
				outline_base_bin(lgTransStackLine, ipBin, phots*profile[ipBin-ip3SigmaRed], inwd, nonScatteredFraction);
		}
	}
}

/*outline_base_bin - adds line photons to bins of reflin and outlin */
void outline_base_bin(bool lgTransStackLine, long int ip, double phots, realnum inwd,
						double nonScatteredFraction)
{
	DEBUG_ENTRY( "outline_base_bin()" );

	if (lgTransStackLine)
	{
		rfield.DiffuseLineEmission[ip] += 
			(realnum)phots;

		/* the reflected part */
		rfield.reflin[0][ip] +=
			(realnum)(inwd*phots*radius.BeamInIn);
		
		/* inward beam that goes out since sphere set */
		rfield.outlin[0][ip] +=
			(realnum)(inwd*phots*radius.BeamInOut*opac.tmn[ip]*nonScatteredFraction);
		
		/* outward part */
		rfield.outlin[0][ip] +=
			(realnum)((1.-inwd)*phots*radius.BeamOutOut*opac.tmn[ip]*nonScatteredFraction);
	}
	else
	{
		rfield.reflin[0][ip] +=
			(realnum)(phots*radius.dVolReflec);

		rfield.outlin[0][ip] +=
			(realnum)(phots*radius.dVolOutwrd*opac.ExpZone[ip]);
	}
}

/*lindst add line with destruction and outward */
static void lindst1(
  double dampXvel,
  double damp,
  // xEmiss - local emissivity per unit vol
  double xEmiss,
  double xEmissIsoBkg,	/* xEmissIsoBkg - local emissivity corrected for isotropic backgrounds per unit vol, no fill fac */
  // wavelength of line in Angstroms
  realnum wavelength,
  // *chLab string label for ion
  const char *chLab,
  // ipnt offset of line in continuum mesh
  long int ipnt,
  // chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line, 't' transferred line
  char chInfo,
  // lgOutToo should line be included in outward beam?
  bool lgOutToo,
  // *chComment string explaining line
  const char *chComment,
  const TransitionProxy& tr )
{
	DEBUG_ENTRY( "lindst1()" );

	// do not add information lines to outward beam
	ASSERT( !lgOutToo || chInfo!='i' );

	lincom(xEmiss, xEmissIsoBkg, wavelength, chLab, ipnt, chInfo, chComment, false, tr );

	if( LineSave.ipass > 0 )
	{
		/* >>chng 06 feb 08, add test on xEmiss positive, no need to evaluate
		 * for majority of zero */
		if (lgOutToo && xEmiss > 0.)
		{
			/* add line to outward beam
			 * there are lots of lines that are sums of other lines, or
			 * just for info of some sort.  These have flag lgOutToo false.
			 * Note that the EnergyRyd variable only has a rational
			 * value if PntForLine was called just before this routine - in all
			 * cases where this did not happen the flag is false. */
			const bool lgTransStackLine = false;
			const long int ip = ipnt - 1;
			const double phots = xEmiss/(rfield.anu(ipnt-1)*EN1RYD);
			const realnum inwd = (realnum)(1.0-(1.+geometry.covrt)/2.);
			const double nonScatteredFraction = 1.;

			outline_base(dampXvel, damp, lgTransStackLine, ip, phots, inwd, nonScatteredFraction);
		}
	}
}

/*lindst add line with destruction and outward */
void lindst(
  // xEmiss - local emissivity per unit vol
  double xEmiss, 
  // wavelength of line in Angstroms
  realnum wavelength, 
  // *chLab string label for ion
  const char *chLab, 
  // ipnt offset of line in continuum mesh
  long int ipnt, 
  // chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line, 't' transferred line
  char chInfo, 
  // lgOutToo should line be included in outward beam?
  bool lgOutToo,
  // *chComment string explaining line 
  const char *chComment )
{
	DEBUG_ENTRY( "lindst()" );
	lindst(0.0, 0.0, xEmiss, wavelength, chLab, ipnt, chInfo, lgOutToo, chComment);
}

/*lindst add line with destruction and outward */
void lindst(
  double dampXvel,
  double damp,
  // xEmiss - local emissivity per unit vol
  double xEmiss,
  // wavelength of line in Angstroms
  realnum wavelength,
  // *chLab string label for ion
  const char *chLab,
  // ipnt offset of line in continuum mesh
  long int ipnt,
  // chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line, 't' transferred line
  char chInfo,
  // lgOutToo should line be included in outward beam?
  bool lgOutToo,
  // *chComment string explaining line
  const char *chComment )
{
	DEBUG_ENTRY( "lindst()" );
	lindst1(dampXvel,damp,xEmiss,xEmiss,wavelength,chLab,ipnt,chInfo,lgOutToo,chComment,TransitionProxy());
}

/*lindst add line with destruction and outward */
void lindst(
	const TransitionProxy& t,
	const ExtraInten& extra,
	// *chLab string label for ion
	const char *chLab, 
	// chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line, 't' transferred line
	char chInfo, 
	// lgOutToo should line be included in outward beam?
	bool lgOutToo,
	// *chComment string explaining line 
	const char *chComment)
{
	DEBUG_ENTRY( "lindst()" );

	// H2O  212.468m
	if (0 && LineSave.ipass > 0)
		if (strncmp(LineSave.lines[LineSave.nsum].chALab(),"H2O ",4) == 0 &&
			 fabs(LineSave.lines[LineSave.nsum].wavelength()-212.468e4) < 1e4)
			fprintf(ioQQQ,"DEBUG lindst: %ld %4ld %15.8e %15.8e %15.8e %15.8e %15.8e\n",
					  LineSave.nsum,nzone,radius.depth,t.Emis().xObsIntensity(),
					  phots( t ),(*t.Hi()).Pop(),t.Emis().Pesc_total());
	lindst1(t.Emis().dampXvel(),
		t.Emis().damp(),
		t.Emis().xIntensity()+extra.v,
		t.Emis().xObsIntensity()+extra.v,
		t.WLAng(), chLab, t.ipCont(), chInfo, lgOutToo, chComment, t );

}

/*PntForLine generate pointer for forbidden line */
void PntForLine(
  /* wavelength of transition in Angstroms */
  double wavelength, 
  /* label for this line */
  const char *chLabel,
  /* this is array index on the f, not c scale,
   * for the continuum cell holding the line */
  long int *ipnt)
{
	/* 
	 * maximum number of forbidden lines - this is a good bet since
	 * new lines do not go into this group, and lines are slowly 
	 * moving to level 1 
	 */
	const int MAXFORLIN = 1000;
	static long int ipForLin[MAXFORLIN]={0};

	/* number of forbidden lines entered into continuum array */
	static long int nForLin;

	DEBUG_ENTRY( "PntForLine()" );

	/* must be 0 or greater */
	ASSERT( wavelength >= 0. );

	if( wavelength == 0. )
	{
		/* zero is special flag to initialize */
		nForLin = 0;
	}
	else
	{

		if( LineSave.ipass > 0 )
		{
			/* not first pass, sum lines only */
			*ipnt = ipForLin[nForLin];
		}
		else if( LineSave.ipass == 0 )
		{
			/* check if number of lines in arrays exceeded */
			if( nForLin >= MAXFORLIN )
			{
				fprintf( ioQQQ, "PROBLEM %5ld lines is too many for PntForLine.\n", 
				  nForLin );
				fprintf( ioQQQ, " Increase the value of maxForLine everywhere in the code.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* ipLineEnergy will only put in line label if nothing already there */
			const double EnergyRyd = RYDLAM/wavelength;
			ipForLin[nForLin] = ipLineEnergy(EnergyRyd,chLabel , 0);
			*ipnt = ipForLin[nForLin];
		}
		else
		{
			/* this is case where we are only counting lines */
			*ipnt = 0;
		}
		++nForLin;
	}
	return;
}

/*ConvRate2CS convert down coll rate back into electron cs in case other parts of code need this for reference */
double ConvRate2CS( realnum gHi , realnum rate )
{

	double cs;

	DEBUG_ENTRY( "ConvRate2CS()" );

	/* return is collision strength, convert from collision rate from 
	 * upper to lower, this assumes pure electron collisions, but that will
	 * also be assumed by anything that uses cs, for self-consistency */
	cs = rate * gHi / dense.cdsqte;

	/* change assert to non-negative - there can be cases (Iin H2) where cs has
	 * underflowed to 0 on some platforms */
	ASSERT( cs >= 0. );
	return cs;
}

/*ConvCrossSect2CollStr convert collisional deexcitation cross section for into collision strength */
double ConvCrossSect2CollStr( double CrsSectCM2, double gLo, double E_ProjectileRyd, double reduced_mass_grams )
{
	// See e.g. Burgess & Tully 1992 A&A 254, 436, S2.1
	double CollisionStrength;

	DEBUG_ENTRY( "ConvCrossSect2CollStr()" );

	if( ! ( CrsSectCM2 >= 0. && gLo >= 0. && E_ProjectileRyd >= 0. && reduced_mass_grams >= 0. ) )
	{
		fprintf( ioQQQ, "invalid parameter for ConvCrossSect2CollStr\n" );
		cdEXIT(EXIT_FAILURE);
	}

	CollisionStrength = CrsSectCM2 * gLo * E_ProjectileRyd / (PI*BOHR_RADIUS_CM*BOHR_RADIUS_CM);

	// this part is being tested.
#if 1
	CollisionStrength *= reduced_mass_grams / ELECTRON_MASS;
#endif

	ASSERT( CollisionStrength >= 0. );
	return CollisionStrength;
}

/*totlin sum total intensity of cooling, recombination, or intensity lines */
double totlin(
	/* chInfor is 1 char, 
	'i' information, 
	'r' recombination or 
	'c' collision */
	int chInfo)
{
	long int i;
	double totlin_v;

	DEBUG_ENTRY( "totlin()" );

	/* routine goes through set of entered line
	 * intensities and picks out those which have
	 * types agreeing with chInfo.  Valid types are
	 * 'c', 'r', and 'i'
	 *begin sanity check */
	if( (chInfo != 'i' && chInfo != 'r') && chInfo != 'c' )
	{
		fprintf( ioQQQ, " TOTLIN does not understand chInfo=%c\n", 
		  chInfo );
		cdEXIT(EXIT_FAILURE);
	}
	/*end sanity check */

	/* now find sum of lines of type chInfo */
	totlin_v = 0.;
	for( i=0; i < LineSave.nsum; i++ )
	{
		if( LineSave.lines[i].chSumTyp() == chInfo )
		{
			totlin_v += LineSave.lines[i].SumLine(0);
		}
	}
	return( totlin_v );
}


/*FndLineHt search through line heat arrays to find the strongest heat source */
const TransitionProxy FndLineHt(long int *level)
{
	TransitionProxy t;
	DEBUG_ENTRY( "FndLineHt()" );

	double Strong = -1.;
	*level = 0;

	/* now do the level 2 lines */
	for( long i=0; i < nWindLine; i++ )
	{
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
		{
			/* check if a line was the major heat agent */
			if( TauLine2[i].Coll().heat() > Strong )
			{
				*level = 2;
				t = TauLine2[i];
				Strong = TauLine2[i].Coll().heat();
			}
		}
	}

	/* now do the hyperfine structure lines */
	for( size_t i=0; i < HFLines.size(); i++ )
	{
		/* check if a line was the major heat agent */
		if( HFLines[i].Coll().heat() > Strong )
		{
			*level = 3;
			t = HFLines[i];
			Strong = HFLines[i].Coll().heat();
		}
	}

	/* lines from external databases */
	for (int ipSpecies=0; ipSpecies < nSpecies; ++ipSpecies)
	{
		for( EmissionList::iterator em=dBaseTrans[ipSpecies].Emis().begin();
			  em != dBaseTrans[ipSpecies].Emis().end(); ++em)
		{
			/* check if a line was the major heat agent */
			if( (*em).Tran().Coll().heat() > Strong )
			{
				*level = 4;
				t = (*em).Tran();
				Strong = t.Coll().heat();
			}
		}
	}

	fixit("all other line stacks need to be included here.");
	// can we just sweep over line stack?  Is that ready yet?

	ASSERT( t.associated() );
	return t;
}



/*set_xIntensity: compute gross and net number of emitted line photons */
void set_xIntensity( const TransitionProxy& t )
{
	DEBUG_ENTRY( "set_xIntensity()" );

	t.Emis().xIntensity() =
	t.Emis().xObsIntensity() = 0.;

	int j = t.ipCont() - 1;
	if( j < 0 )
		return;

	double Pesc = MIN2(1.0, t.Emis().Pesc_total());
	double nphot = t.Emis().Aul() * Pesc * (*t.Hi()).Pop();
	nphot = MAX2(0., nphot);

	t.Emis().xObsIntensity() =
	t.Emis().xIntensity() = nphot * t.EnergyErg();

	if( ! save.lgSubtrCont )
		return;

	double dnphot = t.Emis().Aul() * Pesc * rfield.flux_isotropic[j]*rfield.convoc[j]
		* ((*t.Hi()).Pop() - (*t.Lo()).Pop() * (*t.Hi()).g() / (*t.Lo()).g());

	/* Test if emitting a comment in the main output is appropriate. */
	if( ! LineSave.lgIsoContSubSignif && fabs( dnphot ) > 0.5 * nphot )
		LineSave.lgIsoContSubSignif = true;

	nphot += dnphot;
	nphot = MAX2(0., nphot);

	t.Emis().xObsIntensity() = nphot * t.EnergyErg();

	if( 0 && t.chLabel() == "Fe24                3068.00m" )
	{
		fprintf(ioQQQ,
			"\"%s\"\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n",
			t.chLabel().c_str(),
			t.Emis().Aul(),
			(*t.Hi()).Pop(),
			(*t.Lo()).Pop(),
			rfield.flux_isotropic[j] * rfield.convoc[j],
			Pesc,
			nphot - dnphot,	// undo the correction
			dnphot,
			t.Emis().xIntensity(),
			t.Emis().xObsIntensity() );
	}

	return;
}
