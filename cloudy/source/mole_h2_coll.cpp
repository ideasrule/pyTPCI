/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*H2_CollidRateRead read collision rates */
/*H2_CollidRateEvalAll - set H2 collision rates */
#include "cddefines.h" 
#include "atmdat.h"
#include "phycon.h"
#include "input.h"

STATIC realnum GbarRateCoeff( long nColl, double ediff );

/*H2_CollidRateEvalAll - set H2 collision rate coefficients */
void diatomics::H2_CollidRateEvalAll( void )
{
	DEBUG_ENTRY( "H2_CollidRateEvalAll()" );

	long int numb_coll_trans = 0;

	if( nTRACE >= n_trace_full ) 
		fprintf(ioQQQ,"%s set collision rates\n", label.c_str());
	/* loop over all possible collisional changes within X 
	* and set collision rates, which only depend on Te
	* will go through array in energy order since coll trans do not
	* correspond to a line 
	* collisional dissociation rate coefficient, units cm3 s-1 */
	H2_coll_dissoc_rate_coef[0][0] = 0.;
	H2_coll_dissoc_rate_coef_H2[0][0] = 0.;
	for( long ipHi=0; ipHi<nLevels_per_elec[0]; ++ipHi )
	{
		/* obtain the proper indices for the upper level */
		long iVibHi = ipVib_H2_energy_sort[ipHi];
		long iRotHi = ipRot_H2_energy_sort[ipHi];

		/* this is a guess of the collisional dissociation rate coefficient -
		 * will be multiplied by the sum of densities of all colliders 
		 * except H2*/
		double energy = H2_DissocEnergies[0] - states[ipHi].energy().WN();
		ASSERT( energy > 0. );
		/* we made this up - Boltzmann factor times rough coefficient */
		H2_coll_dissoc_rate_coef[iVibHi][iRotHi] = 
			1e-14f * (realnum)sexp(energy/phycon.te_wn) * lgColl_dissoc_coll;

		/* collisions with H2 - pre coefficient changed from 1e-8 
		 * (from umist) to 1e-11 as per extensive discussion with Phillip Stancil */
		H2_coll_dissoc_rate_coef_H2[iVibHi][iRotHi] = 
			1e-11f * (realnum)sexp(energy/phycon.te_wn) * lgColl_dissoc_coll;

		/*fprintf(ioQQQ,"DEBUG coll_dissoc_rateee\t%li\t%li\t%.3e\t%.3e\n",
			iVibHi,iRotHi,
			H2_coll_dissoc_rate_coef[iVibHi][iRotHi],
			H2_coll_dissoc_rate_coef_H2[iVibHi][iRotHi]);*/

		for( long ipLo=0; ipLo<ipHi; ++ipLo )
		{
			long iVibLo = ipVib_H2_energy_sort[ipLo];
			long iRotLo = ipRot_H2_energy_sort[ipLo];

			ASSERT( states[ipHi].energy().WN() - states[ipLo].energy().WN() > 0. );

			/* keep track of number of different collision routes */
			++numb_coll_trans;
			for( long nColl=0; nColl<N_X_COLLIDER; ++nColl )
			{
				realnum newrate = H2_CollidRateEvalOne( iVibHi,iRotHi,iVibLo,iRotLo,
					ipHi , ipLo , nColl, phycon.te );
				CollRateCoeff[ipHi][ipLo][nColl] = newrate;
			} 
		}
	}

	fixit("test that this does not matter with new rate table interpolation");
	//Remove if it does not.
#if 0
	/* >>chng 05 nov 30, GS, rates decreases exponentially for low temperature, see Le Bourlot et al. 1999  */
	/* Phillips mail--Apparently, the SD fit is only valid over the range of their  calculations, 100-1000K. 
	* The rate should continue to fall exponentially with decreasing T, something like exp(-3900/T) for 0->1 and  
	* exp[-(3900-170.5)/T] for 1->0. It is definitely, not constant for T  lower than 100 K, as far as we know. 
	* There may actually be a quantum tunneling effect which causes the rate to increase at lower T, but no  
	* one has calculated it (as far as I know) and it might happen below 1K or  so.???*/
	if( phycon.te < 100. )
	{
		/* first term in exp is suggested by Phillip, second temps in paren is to ensure continuity
		* across 100K */
		H2_CollRate[0][1][0][0][0] *= (realnum)(exp(-(3900-170.5)*(1./phycon.te - 1./100.)));
		H2_CollRate[0][3][0][0][0] *= (realnum)(exp(-(3900-1015.1)*(1./phycon.te - 1./100.)));
		H2_CollRate[0][2][0][1][0] *= (realnum)(exp(-(3900-339.3)*(1./phycon.te - 1./100.)));
	}
#endif

	if( nTRACE >= n_trace_full )
		fprintf(ioQQQ,
		" collision rates updated for new temp, number of trans is %li\n",
		numb_coll_trans);

	return;
}

/* compute rate coefficient for a single quenching collision */
realnum diatomics::H2_CollidRateEvalOne( 
	/*returns collision rate coefficient, cm-3 s-1 for quenching collision
	 * from this upper state */
	 long iVibHi, long iRotHi,long iVibLo,
	 /* to this lower state */
	long iRotLo, long ipHi , long ipLo , 
	/* colliders are H, He, H2(ortho), H2(para), and H+ */
	long  nColl,
	double te_K )
{
	DEBUG_ENTRY( "H2_CollidRateEvalOne()" );

	realnum rate = InterpCollRate( RateCoefTable[nColl], ipHi, ipLo, te_K );

	/* this is option to use guess of collision rate coefficient */
	if( (rate==0.f) &&
		/* turn lgColl_gbar on/off with atom h2 gbar on off */
		lgColl_gbar &&
		/* but only if this does not mix ortho and para */
		(H2_lgOrtho[0][iVibHi][iRotHi]==H2_lgOrtho[0][iVibLo][iRotLo] ) )
	{
		double ediff = states[ipHi].energy().WN() - states[ipLo].energy().WN();
		rate = GbarRateCoeff( nColl, ediff );
	}

	rate *= lgColl_deexec_Calc;

	/* >>chng 05 feb 09, add option to kill ortho - para collisions */
	if( !lgH2_ortho_para_coll_on && 
		(H2_lgOrtho[0][iVibHi][iRotHi] != H2_lgOrtho[0][iVibLo][iRotLo]) )
		rate = 0.;

	if( lgH2_NOISE )
		rate *= CollRateErrFac[ipHi][ipLo][nColl];

	return rate;
}

STATIC realnum GbarRateCoeff( long nColl, double ediff )
{
	// these are fits to the existing collision data 
	// used to create g-bar rates 
	double gbarcoll[N_X_COLLIDER][3] = 
	{
		{-9.9265 , -0.1048 , 0.456  },
		{-8.281  , -0.1303 , 0.4931 },
		{-10.0357, -0.0243 , 0.67   },
		{-8.6213 , -0.1004 , 0.5291 },
		{-9.2719 , -0.0001 , 1.0391 }
	};

	// do not let energy difference be smaller than 100 wn, the smallest
	// difference for which we fit the rate coefficients
	ediff = MAX2(100., ediff );

	// the fit is log(K)=y_0+a*((x)^b), where K is the rate coefficient,
	// and x is the energy in wavenumbers
	realnum rate = (realnum)exp10(
		gbarcoll[nColl][0] + gbarcoll[nColl][1] * 
		pow(ediff,gbarcoll[nColl][2]) );

	return rate;
}

void diatomics::H2_CollidRateRead( long int nColl )
{
	DEBUG_ENTRY( "H2_CollidRateRead()" );

	if( coll_source[nColl].filename.size() == 0 && coll_source[nColl].magic == 0 )
		return;

	const char* chFilename = coll_source[nColl].filename.c_str();
	long magic_expect = coll_source[nColl].magic;
	string chPath = path + cpu.i().chDirSeparator() + coll_source[nColl].filename;
	FILE *ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	string chLine;
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " H2_CollidRateRead could not read first line of %s\n", chFilename );
		cdEXIT(EXIT_FAILURE);
	}

	/* magic number */
	long n1 = atoi( chLine.c_str() );
	if( n1 != magic_expect )
	{
		fprintf( ioQQQ, " H2_CollidRateRead: the version of %s is not the current version.\n", chFilename );
		fprintf( ioQQQ, " I expected to find the number %li and got %li instead.\n", magic_expect, n1 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine.c_str() );
		cdEXIT(EXIT_FAILURE);
	}

	FunctPtr func = new FunctDiatoms( *this );
	ReadCollisionRateTable( RateCoefTable[nColl], ioDATA, func, nLevels_per_elec[0] );
	delete func;

	fclose( ioDATA );

	return;
}

void diatomics::GetIndices( long& ipHi, long& ipLo, const char* chLine, long& i ) const
{
	bool lgEOL;
	long iVibHi = (long)FFmtRead( chLine, &i, strlen(chLine), &lgEOL );
	long iRotHi = (long)FFmtRead( chLine, &i, strlen(chLine), &lgEOL );
	long iVibLo = (long)FFmtRead( chLine, &i, strlen(chLine), &lgEOL );
	long iRotLo = (long)FFmtRead( chLine, &i, strlen(chLine), &lgEOL );
	ASSERT( iRotHi >= 0 && iVibHi >= 0 && iRotLo >= 0 && iVibLo >=0 );

	// check that we actually included the levels in the model representation
	// depending on the potential surface, the collision data set
	// may not agree with our adopted model - skip those 
	if(  ( iVibHi > nVib_hi[0] || iVibLo > nVib_hi[0] ) ||
		( iRotHi < Jlowest[0] || iRotLo < Jlowest[0] ) ||
		( iRotHi > nRot_hi[0][iVibHi] || iRotLo > nRot_hi[0][iVibLo] ) ||
		/* some collision rates have the same upper and lower indices - skip them */
		( iVibHi == iVibLo && iRotHi == iRotLo ) )
	{
		ipHi = -1;
		ipLo = -1;
		return;
	}

	ipHi = ipEnergySort[0][iVibHi][iRotHi];
	ipLo = ipEnergySort[0][iVibLo][iRotLo];
	if( ipHi < ipLo )
		swap( ipHi, ipLo );

	return;
}

