/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_helium put He-like iso sequence into line intensity stack */
/*TempInterp interpolates on a grid of values to produce predicted value at current Te.*/
#include "cddefines.h"
#include "dense.h"
#include "prt.h"
#include "helike.h"
#include "iso.h"
#include "atmdat.h"
#include "lines.h"
#include "phycon.h"
#include "taulines.h"
#include "thirdparty.h"
#include "trace.h"
#include "freebound.h"
#include "two_photon.h"
#include "lines_service.h"
#include "parser.h"
#include "container_classes.h"

static const int NUMTEMPS = 21;
static const int NUMDENS = 14;

STATIC void GetStandardHeLines(void);
STATIC double TempInterp2( double* TempArray , double* ValueArray, long NumElements, double Te );
STATIC void DoSatelliteLines( long nelem );

static bool lgFirstRun = true;
static double CaBDensities[NUMDENS];
static double CaBTemps[NUMTEMPS];
static long NumLines;
static multi_arr<double,3> CaBIntensity;

inline void setup_multiplet( LinSv *line,
							 vector<long>& multiplet )
{
	if (LineSave.ipass == 0) {
		for( auto tr = multiplet.begin(); tr != multiplet.end(); ++tr )
			line->addComponentID( *tr );
		line->setBlendWavl();
	}
}

inline void randomize_inten( t_iso_sp* sp, long ipLo, long ipHi ) 
{ 
	sp->trans(ipHi,ipLo).Emis().xIntensity() *= sp->ex[ipHi][ipLo].ErrorFactor[IPRAD]; 
	sp->trans(ipHi,ipLo).Emis().xObsIntensity() *= sp->ex[ipHi][ipLo].ErrorFactor[IPRAD]; 
	return; 
} 


void lines_helium()
{
	long ipISO = ipHE_LIKE;
	string chLabel="    ";

	double log10_eden = log10(dense.eden);

	DEBUG_ENTRY( "lines_helium()" );

	if( trace.lgTrace )
		fprintf( ioQQQ, "   prt_lines_helium called\n" );

	// this can be changed with the atom levels command but must be at least 3.
	ASSERT( !dense.lgElmtOn[ipHELIUM] || iso_sp[ipHE_LIKE][ipHELIUM].n_HighestResolved_max >= 3 );

	long i = StuffComment( "He-like iso-sequence" );
	linadd( 0., (realnum)i , "####", 'i',
		" start He-like iso sequence");

	/* read in Case A and B lines from data file	*/
	if( lgFirstRun )
	{
		GetStandardHeLines();
		lgFirstRun = false;
	}

	/* this is the main printout, where line intensities are entered into the stack */
	for( long nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		vector<long> multiplet;

		if( dense.lgElmtOn[nelem] )
		{
			t_iso_sp* sp = &iso_sp[ipHE_LIKE][nelem];

			ASSERT( sp->n_HighestResolved_max >= 3 );

			// add two-photon details here
			if( LineSave.ipass == 0 )
			{
				/* chIonLbl is function that generates a null terminated 4 char string, of form "C  2" 
				 * the result, chLable, is only used when ipass == 0, can be undefined otherwise */
				chLabel = chIonLbl(nelem+1, nelem+1-ipISO);
			}
			for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
			{
				fixit("This was multiplied by Pesc when treated as a line, now what?  Only used for printout?");
				fixit("below should be 'i' instead of 'r' ?");

				string tpc_comment = "";
				if( LineSave.ipass == 0 )
				{
					tpc_comment = " two photon continuum, " +
						iso_comment_tran_levels( ipISO, nelem, tnu->ipLo, tnu->ipHi );
				}

				linadd(	tnu->AulTotal * tnu->E2nu * EN1RYD * (*tnu->Pop), 
					2. * wn2ang( (*sp).trans( tnu->ipHi, tnu->ipLo ).EnergyWN() ),
					chLabel.c_str(), 'r', tpc_comment.c_str() );
			}

			/* here we will create an entry for the three lines 
			 * coming from 2 3P to 1 1S - the entry called TOTL will
			 * appear before the lines of the multiplet */
			for( long i=ipHe2p3P0; i <= ipHe2p3P2; i++ )
			{
				if( sp->trans(i, ipHe1s1S).ipCont() <= 0 )
					continue;
				set_xIntensity( sp->trans(i, ipHe1s1S) );
				// save index into line stack before we call PutLine()
				// this is much more efficient than findline() to get blend components
				multiplet.push_back( LineSave.nsum );
				PutLine( sp->trans(i, ipHe1s1S),
					iso_comment_tran_levels( ipISO, nelem, ipHe1s1S, i ).c_str(),
					chLabel.c_str() );
			}

			if( multiplet.size() > 0 )
			{
				LinSv *lineHe1 =
					linadd( 0.0, 0.0, "Blnd", 'i',
						" total emission in He-like intercombination lines from 2p3P to ground ");
				setup_multiplet( lineHe1, multiplet );
				multiplet.resize( 0 );
			}

			/* set number of levels to print */
			long int nLoop  = iso_Max_Emitting_Level(nelem, ipISO, prt.lgPrnIsoCollapsed);

			/* now do real permitted lines */
			/* NB NB - low and high must be in this order so that all balmer, paschen,
			 * etc series line up correctly in final printout */
			/* >>chng 01 jun 13, bring 23P lines back together */
			for( long ipLo=0; ipLo < ipHe2p3P0; ipLo++ )
			{
				vector<long> EnterTheseLast;
				for( long ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					/* >>chng 01 may 30, do not add fake he-like lines (majority) to line stack */
					/* >>chng 01 dec 11, use variable for smallest A */
					if( sp->trans(ipHi,ipLo).ipCont() < 1 )
						continue;

					if( ipLo==ipHe2s3S && ipHi == ipHe2p3P0 )
					{
						/* here we will create an entry for the three lines 
						 * coming from 2 3P to 2 3S - the entry called TOTL will
						 * appear before the lines of the multiplet 
						 * for He I this is 10830 */

						for( long i=ipHe2p3P0; i <= ipHe2p3P2; i++ )
						{
							if( sp->trans(i, ipLo).ipCont() <= 0 )
								continue;
							set_xIntensity( sp->trans(i, ipLo) );

							multiplet.push_back( LineSave.nsum );
							PutLine( sp->trans(i, ipLo),
								iso_comment_tran_levels( ipISO, nelem, ipLo, i ).c_str(),
								chLabel.c_str() );

							/* >>chng 13-jun-06
							 * correct for isotropic continuum before applying error randomization
							 * to avoid shutting off emission lines (if the correction is applied _after_
							 * the error is computed, it is likely to be higher than the updated intensity)
							 */
							if( iso_ctrl.lgRandErrGen[ipISO] )
							{
								randomize_inten( sp, ipLo, i );
							}
						}

						if( multiplet.size() > 0 )
						{
							LinSv *lineHe1 =
								linadd( 0.0, 0.0, "Blnd", 'i',
									"total emission in He-like lines, use wgt average of three line wavelengths " );
							setup_multiplet( lineHe1, multiplet );
							multiplet.resize( 0 );
						}
					}
					else if( ipHi >= ipHe2p3P0 && ipHi <= ipHe2p3P2 )
					{
						// already done above
						continue;
					}
					else
					{
						set_xIntensity( sp->trans(ipHi,ipLo) );

						if( iso_ctrl.lgRandErrGen[ipISO] )
						{
							randomize_inten( sp, ipLo, ipHi );
						}

						if( abs( L_(ipHi) - L_(ipLo) ) != 1 )
						{
							EnterTheseLast.push_back( ipHi );
							continue;
						}

						/* 
						fprintf(ioQQQ,"1 loop %li %li %.1f\n", ipLo,ipHi, 
							sp->trans(ipHi,ipLo).WLAng() ); */
						string comment_trans = "";
						if( LineSave.ipass == 0 )
						{
							comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi );
						}
						PutLine( sp->trans(ipHi,ipLo), comment_trans.c_str() );

						{
							/* option to print particulars of some line when called
							 * a prettier print statement is near where chSpin is defined below*/
							enum {DEBUG_LOC=false};
							if( DEBUG_LOC )
							{
								if( nelem==1 && ipLo==0 && ipHi==1 )
								{
									fprintf(ioQQQ,"he1 626 %.2e %.2e \n", 
										sp->trans(ipHi,ipLo).Emis().TauIn(),
										sp->trans(ipHi,ipLo).Emis().TauTot()
										);
								}
							}
						}
					}
				}

				// enter these lines last because they are generally weaker quadrupole transitions.
				for( vector<long>::iterator it = EnterTheseLast.begin(); it != EnterTheseLast.end(); it++ )
				{
					string comment_trans = "";
					if( LineSave.ipass == 0 )
					{
						comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, *it );
					}
					PutLine( sp->trans(*it,ipLo), comment_trans.c_str() );
				}
			}

			/* this sum will bring together the three lines going to J levels within 2 3P */
			for( long ipHi=ipHe2p3P2+1; ipHi < nLoop; ipHi++ )
			{
				for( long ipLo=ipHe2p3P0; ipLo <= ipHe2p3P2; ++ipLo )
				{
					if( sp->trans(ipHi,ipLo).ipCont() <= 0 ) 
						continue;

					set_xIntensity( sp->trans(ipHi, ipLo) );

					if( iso_ctrl.lgRandErrGen[ipISO] )
					{
						randomize_inten( sp, ipLo, ipHi );
					}

					multiplet.push_back( LineSave.nsum );
					PutLine( sp->trans(ipHi, ipLo),
						iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi ).c_str(),
						chLabel.c_str() );
				}

				if( multiplet.size() > 0 )
				{
					LinSv *lineHe1 =
						linadd( 0.0, 0.0, "Blnd", 'i',
							"total emission in He-like lines, use wgt average of three line wavelengths " );
					setup_multiplet( lineHe1, multiplet );
					multiplet.resize( 0 );
				}
			}
			for( long ipLo=ipHe2p3P2+1; ipLo < nLoop-1; ipLo++ )
			{
				vector<long> EnterTheseLast;
				for( long ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					/* skip non-radiative lines */
					if( sp->trans(ipHi,ipLo).ipCont() < 1 ) 
						continue;

					set_xIntensity( sp->trans(ipHi,ipLo) );

					if( iso_ctrl.lgRandErrGen[ipISO] )
					{
						randomize_inten( sp, ipLo, ipHi );
					}

					if( N_(ipHi) > sp->n_HighestResolved_max || abs( L_(ipHi) - L_(ipLo) ) != 1 )
					{
						EnterTheseLast.push_back( ipHi );
						continue;
					}

					string comment_trans = "";
					if( LineSave.ipass == 0 )
					{
						comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi );
					}
					PutLine(sp->trans(ipHi,ipLo), comment_trans.c_str() );
				}

				// enter these lines last because they are generally weaker quadrupole transitions.
				for( vector<long>::iterator it = EnterTheseLast.begin(); it != EnterTheseLast.end(); it++ )
				{
					string comment_trans = "";
					if( LineSave.ipass == 0 )
					{
						comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, *it );
					}
					PutLine( sp->trans(*it,ipLo), comment_trans.c_str() );
				}
			}

			/* Now put the satellite lines in */
			if( iso_ctrl.lgDielRecom[ipISO] )
				DoSatelliteLines(nelem);

			if( nelem == ipHELIUM )
			{
				for( long i=0; i< NumLines; i++ )
				{
					double intens_at_Te[NUMDENS];
					for( long ipDens = 0; ipDens < NUMDENS; ++ipDens )
						intens_at_Te[ipDens] = TempInterp2( CaBTemps, &CaBIntensity[i][ipDens][0], NUMTEMPS, phycon.te );
					double intens = linint( CaBDensities, intens_at_Te, NUMDENS, log10_eden );
					intens = exp10(intens) * dense.xIonDense[nelem][nelem+1-ipISO]*dense.eden;
					ASSERT( intens >= 0. );
					linadd( intens, atmdat.CaseBWlHeI[i], "Ca B", 'i', "Case B intensity " );
				}
			}
		}
	}

	if( iso_sp[ipHE_LIKE][ipHELIUM].n_HighestResolved_max >= 4 &&
		iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max+iso_sp[ipH_LIKE][ipHYDROGEN].nCollapsed_max >= 8 )
	{
		t_iso_sp* sp = &iso_sp[ipHE_LIKE][ipHELIUM];
		long ipHe4s3S = iso_sp[ipHE_LIKE][ipHELIUM].QN2Index(4, 0, 3);
		long ipHe4p3P = iso_sp[ipHE_LIKE][ipHELIUM].QN2Index(4, 1, 3);

		/* For info only, add the total photon flux in 3889 and 7065,
		* and 3188, 4713, and 5876. */
		double photons_3889_plus_7065 =
			/* to 2p3P2 */
			phots( sp->trans(ipHe3s3S,ipHe2p3P2) ) +
			phots( sp->trans(ipHe3d3D,ipHe2p3P2) ) +
			phots( sp->trans(ipHe4s3S,ipHe2p3P2) ) +
			/* to 2p3P1 */
			phots( sp->trans(ipHe3s3S,ipHe2p3P1) ) +
			phots( sp->trans(ipHe3d3D,ipHe2p3P1) ) +
			phots( sp->trans(ipHe4s3S,ipHe2p3P1) ) +
			/* to 2p3P0 */
			phots( sp->trans(ipHe3s3S,ipHe2p3P0) ) +
			phots( sp->trans(ipHe3d3D,ipHe2p3P0) ) +
			phots( sp->trans(ipHe4s3S,ipHe2p3P0) ) +
			/* to 2s3S */
			phots( sp->trans(ipHe3p3P,ipHe2s3S) ) +
			phots( sp->trans(ipHe4p3P,ipHe2s3S) ) ;

		long upperIndexofH8 = iso_sp[ipH_LIKE][ipHYDROGEN].QN2Index(8, 1, 2);

		/* Add in photon flux of H8 3889 */
		photons_3889_plus_7065 += 
			phots( iso_sp[ipH_LIKE][ipHYDROGEN].trans(upperIndexofH8,ipH2s) );

		/* now multiply by ergs of normalization line, so that relative flux of
		* this line will be ratio of photon fluxes. */
		if( LineSave.WavLNorm > 0 )
			photons_3889_plus_7065 *= (ERG1CM*1.e8)/LineSave.WavLNorm;
		linadd( photons_3889_plus_7065, 3889., "Pho+", 'i',
			"photon sum given in Porter et al. 2007 (astro-ph/0611579)");
	}

	/* ====================================================
	 * end helium
	 * ====================================================*/

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_helium returns\n" );
	}
	return;
}

STATIC void GetStandardHeLines()
{
	DEBUG_ENTRY( "GetStandardHeLines()" );

	DataParser d("he1_case_b.dat", ES_NONE);

	/* check that magic number is ok */
	d.getline();
	d.checkMagic(CASEBMAGIC);
	d.getToken(NumLines);

	/* get the array of densities */
	d.getline();
	d.getToken(CaBDensities, NUMDENS);
	d.checkEOL();

	/* get the array of temperatures */
	d.getline();
	d.getToken(CaBTemps, NUMTEMPS);
	d.checkEOL();

	/* create space for array of values */
	CaBIntensity.alloc(NumLines, NUMDENS, NUMTEMPS);

	/* now read in the case B line data */
	long line = 0;
	while( d.getline() )
	{
		// the first number is the wavelength, which is not used
		double wl;
		d.getToken(wl);
		long nLo, lLo, sLo, jLo, nHi, lHi, sHi, jHi;
		d.getToken(nLo);
		d.getToken(lLo);
		d.getToken(sLo);
		d.getToken(jLo);
		d.getToken(nHi);
		d.getToken(lHi);
		d.getToken(sHi);
		d.getToken(jHi);
		double Elo = iso_sp[ipHE_LIKE][ipHELIUM].energy(nLo, lLo, sLo, 2*jLo+1);
		double Ehi = iso_sp[ipHE_LIKE][ipHELIUM].energy(nHi, lHi, sHi, 2*jHi+1);
		if( Elo >= 0. && Ehi >= 0. )
		{
			realnum Enerwn = fabs(Ehi - Elo);
			atmdat.CaseBWlHeI.emplace_back( (realnum)wn2ang( (double)Enerwn ) );
		}
		else
			atmdat.CaseBWlHeI.emplace_back( wl );
		d.checkEOL();

		for( long ipDens = 0; ipDens < NUMDENS; ++ipDens )
		{
			d.getline();
			long den;
			d.getToken(den);
			if( den != ipDens + 1 )
				d.errorAbort("invalid density index");
			d.getToken(&CaBIntensity[line][ipDens][0], NUMTEMPS);
			d.checkEOL();
		}
		line++;
	}
	d.checkEOD();

	ASSERT( line == NumLines );
	ASSERT( atmdat.CaseBWlHeI.size() == (unsigned)line );
}

/** \todo	there is a virtually identical routine in helike_recom.cpp -> combine */
STATIC double TempInterp2( double* TempArray , double* ValueArray, long NumElements, double Te )
{
	long int ipTe=-1;
	double rate = 0.;
	long i0;

	DEBUG_ENTRY( "TempInterp2()" );

	/* te totally unknown */
	if( Te <= TempArray[0] )
	{
		return ValueArray[0] + log10( TempArray[0]/Te );
	}
	else if( Te >= TempArray[NumElements-1] )
	{
		return ValueArray[NumElements-1];
	}

	/* now search for temperature */
	ipTe = hunt_bisect( TempArray, NumElements, Te );			

	ASSERT( (ipTe >=0) && (ipTe < NumElements-1)  );

	/* Do a four-point interpolation */
	const int ORDER = 3; /* order of the fitting polynomial */
	i0 = max(min(ipTe-ORDER/2,NumElements-ORDER-1),0);
	rate = lagrange( &TempArray[i0], &ValueArray[i0], ORDER+1, Te );

	return rate;
}

/** \todo	2	say where these come from	*/	
/* For double-ionization discussions, see Lindsay, Rejoub, & Stebbings 2002	*/
/* Also read Itza-Ortiz, Godunov, Wang, and McGuire 2001.	*/
STATIC void DoSatelliteLines( long nelem )
{
	long ipISO = ipHE_LIKE;
	
	DEBUG_ENTRY( "DoSatelliteLines()" );

	ASSERT( iso_ctrl.lgDielRecom[ipISO] && dense.lgElmtOn[nelem] );

	for( long i=0; i < iso_sp[ipISO][nelem].numLevels_max; i++ )
	{
		double dr_rate = iso_sp[ipISO][nelem].fb[i].DielecRecomb;
		TransitionProxy tr = SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][i]];

		tr.Emis().xObsIntensity() = 
		tr.Emis().xIntensity() = 
			dr_rate * dense.eden * dense.xIonDense[nelem][nelem+1-ipISO] * ERG1CM * tr.EnergyWN();
		tr.Emis().pump() = 0.;

		PutLine( tr, "iso satellite line" );
	}

	return;
}
