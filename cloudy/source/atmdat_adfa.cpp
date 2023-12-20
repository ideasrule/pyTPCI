/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*phfit derive photoionization cross sections for first 30 elements */
#include "cddefines.h"

#include "atmdat_adfa.h"

#include "atmdat.h"
#include "iso.h"
#include "parser.h"

/** constructor: read in all the ADfA data files */
t_ADfA::t_ADfA()
{
	DEBUG_ENTRY( "t_ADfA()" );

	/* this option, use the new atmdat_rad_rec recombination rates */
	version = PHFIT_UNDEF;

	double help[9];
	const long VERSION_MAGIC = 20061204L;

	DataParser d;

	d.open( "phfit.dat", ES_NONE );
	d.getline();
	d.checkMagic(VERSION_MAGIC);

	d.getline();
	for( long i = 0; i < NSHELLS; i++ ) d.getToken( L[i] );

	d.getline();
	for( long i = 0; i < LIMELM; i++ ) d.getToken( NINN[i] );

	d.getline();
	for( long i = 0; i < LIMELM; i++ ) d.getToken( NTOT[i] );

	/* 1995 table
	 * >>refer	all	photo_cs	Verner and Yakovlev, 1995, A&AS, 109, 125
	 */
	while( d.getline() )
	{
		long nshell, nelec, nelem;
		d.getToken(nshell);
		d.getToken(nelec);
		d.getToken(nelem);

		if( nshell < 0 && nelec < 0 && nelem < 0 )
		{
			d.checkEOL();
			break;
		}

		ASSERT( 0 <= nshell && nshell < NSHELLS );
		ASSERT( 0 <= nelem && nelem < LIMELM );
		ASSERT( 0 <= nelec && nelec <= nelem );

		for( long id = 0; id < NFIT_PH1; id++ )
			d.getToken( PH1[nshell][nelec][nelem][id] );
		d.checkEOL();
	}

	/* 1996 table
	 * >>refer	all	photo_cs	Verner, D. A., Ferland, G. J., Korista, K. T., & Yakovlev, D. G. 1996, ApJ, 465, 487.
	 */
	while( d.getline() )
	{
		long nelec, nelem;
		d.getToken(nelec);
		d.getToken(nelem);

		if( nelec < 0 && nelem < 0 )
		{
			d.checkEOL();
			break;
		}

		ASSERT( 0 <= nelem && nelem < LIMELM );
		ASSERT( 0 <= nelec && nelec <= nelem );

		for( long id = 0; id < NFIT_PH2; id++ )
			d.getToken( PH2[nelec][nelem][id] );
		d.checkEOL();
	}

	d.checkEOD();

	//--------------------------------------------

	static const char chFile2[] = "hpfit.dat";

	FILE *io = open_data( chFile2, "r" );

	bool lgErr = false;
	long i=-1, j=-1;

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile2, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	for( i=0; i < NHYDRO_MAX_LEVEL; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le", &help[0], &help[1],
					   &help[2], &help[3], &help[4] ) != 5 );
		for( int l=0; l < 5; ++l )
			PHH[i][l] = (realnum)help[l];
	}

	fclose( io );

	ASSERT( !lgErr );

	static const char chFile3[] = "rec_lines.dat";

	io = open_data( chFile3, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile3, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	for( i=0; i < 110; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le %le %le %le", &help[0], &help[1], &help[2],
					   &help[3], &help[4], &help[5], &help[6], &help[7] ) != 8 );
		for( int l=0; l < 8; ++l )
			P[l][i] = (realnum)help[l];
	}


	for( i=0; i < 405; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le %le %le %le %le", &help[0], &help[1], &help[2],
					   &help[3], &help[4], &help[5], &help[6], &help[7], &help[8] ) != 9 );
		for( int l=0; l < 9; ++l )
			ST[l][i] = (realnum)help[l];
	}

	fclose( io );

	ASSERT( !lgErr );

	static const char chFile4[] = "rad_rec.dat";

	io = open_data( chFile4, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile4, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld %ld", &i, &j ) != 2 );
		if( i == -1 && j == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le", &help[0], &help[1] ) != 2 );
		for( int l=0; l < 2; ++l )
			rrec[i][j][l] = (realnum)help[l];
	}
	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld %ld", &i, &j ) != 2 );
		if( i == -1 && j == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le", &help[0], &help[1],
					   &help[2], &help[3] ) != 4 );
		for( int l=0; l < 4; ++l )
			rnew[i][j][l] = (realnum)help[l];
	}
	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
		if( i == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le %le", &help[0], &help[1], &help[2] ) != 3 );
		for( int l=0; l < 3; ++l )
			fe[i][l] = (realnum)help[l];
	}

	fclose( io );

	ASSERT( !lgErr );

	static const char chFile5[] = "h_rad_rec.dat";

	io = open_data( chFile5, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile5, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	for( i=0; i < NHYDRO_MAX_LEVEL; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le %le %le %le %le", &help[0], &help[1], &help[2],
					   &help[3], &help[4], &help[5], &help[6], &help[7], &help[8] ) != 9 );
		for( int l=0; l < 9; ++l )
			HRF[i][l] = (realnum)help[l];
	}

	fclose( io );

	ASSERT( !lgErr );

	static const char chFile6[] = "h_phot_cs.dat";

	io = open_data( chFile6, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile6, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	for( i=0; i < NHYDRO_MAX_LEVEL; i++ )
	{
		lgErr = lgErr || ( fscanf( io, "%le", &help[0] ) != 1 );
		STH[i] = (realnum)help[0];
	}

	fclose( io );

	ASSERT( !lgErr );

	static const char chFile7[] = "coll_ion.dat";

	io = open_data( chFile7, "r" );

	lgErr = lgErr || ( fscanf( io, "%ld", &i ) != 1 );
	if( lgErr || i != VERSION_MAGIC )
	{
		fprintf( ioQQQ, " File %s has incorrect version: %ld\n", chFile7, i );
		fprintf( ioQQQ, " I expected to find version: %ld\n", VERSION_MAGIC );
		cdEXIT(EXIT_FAILURE);
	}

	while( true )
	{
		lgErr = lgErr || ( fscanf( io, "%ld %ld", &i, &j ) != 2 );
		if( i == -1 && j == -1 )
			break;
		lgErr = lgErr || ( fscanf( io, "%le %le %le %le %le", &CF[i][j][0], &CF[i][j][1],
					   &CF[i][j][2], &CF[i][j][3], &CF[i][j][4] ) != 5 );
	}

	fclose( io );

	ASSERT( !lgErr );

	/*refer	HI	cs	Anderson, H., Ballance, C.P., Badnell, N.R., 
	 *refercon	& Summers, H.P  2000, J Phys B, 33, 1255 */
	d.open( "h_coll_str.dat", ES_NONE );
	d.getline();
	d.checkMagic(VERSION_MAGIC);

	while( d.getline() )
	{
		long nHi, lHi, nLo, lLo;
		d.getToken(nHi);
		if( nHi < 0 )
			break;
		d.getToken(lHi);
		d.getToken(nLo);
		d.getToken(lLo);

		QNPair inPair(nHi, lHi, 2, -1, nLo, lLo, 2, -1);
		auto p = HCS[inPair].data();
		d.getToken(p, NHCSTE);
		d.checkEOL();
	}
}

inline long t_ADfA::set_vshell_index( long nelec, long Z ) const
{
	long nout = NTOT[nelec-1];

	// NB NB
	// 4th row elements start filling 4s (shell 7) before completely
	// filling 3d (shell 6) -- true for all neutrals, and some singly
	// ionized elements, see below
	//
	// The indices below are on C scale, i.e., they start at 0 for H,
	// hence the +1
	//
	if( Z == nelec && Z > ipARGON+1 )
		nout = 7;
	else if( Z == (nelec + 1)
		 && ( Z == ipCALCIUM+1
			|| Z == ipSCANDIUM+1
			|| Z == ipTITANIUM+1
			|| Z == ipMANGANESE+1
			|| Z == ipIRON+1 ) )
		nout = 7;
	return nout;
}

double t_ADfA::getEthresh( long int nshell, long int nelec, long int Z ) const
{
	long nout = set_vshell_index( nelec, Z );
	if( nshell == nout )
		return atmdat.getIonPot(Z-1, Z-nelec) * EVRYD;
	else
		return double(PH1[nshell-1][nelec-1][Z-1][0]);
}

double t_ADfA::phfit(long int nz, 
		     long int ne,
		     long int is, 
		     double e)
{
	long int nint, 
	  nout;
	double a, 
	  b, 
	  crs, 
	  einn, 
	  p1, 
	  q, 
	  x, 
	  y, 
	  z;

	DEBUG_ENTRY( "phfit()" );

	/*** Version 3. October 8, 1996.
	 *** Written by D. A. Verner, verner@pa.uky.edu
	 *** Inner-shell ionization energies of some low-ionized species are slightly
	 *** improved to fit smoothly the experimental inner-shell ionization energies 
	 *** of neutral atoms.
	 ******************************************************************************
	 *** This subroutine calculates partial photoionization cross sections
	 *** for all ionization stages of all atoms from H to Zn (Z=30) by use of
	 *** the following fit parameters:
	 *** Outer shells of the Opacity Project (OP) elements:
	 *** >>refer	all	photo_cs	Verner, D. A., Ferland, G. J., Korista, K. T., & Yakovlev, D. G. 1996, ApJ, 465, 487.
	 *** Inner shells of all elements, and outer shells of the non-OP elements:
	 ***  Verner and Yakovlev, 1995, A&AS, 109, 125
	 *** Input parameters:  nz - atomic number from 1 to 30 (integer) 
	 ***          ne - number of electrons from 1 to iz (integer)
	 ***          is - shell number (integer)
	 ***          e - photon energy, eV 
	 ***          version - enum, PHFIT96 (default): calculates 
	 ***                 new cross sections, PHFIT95: calculates
	 ***                 only old Hartree-Slater cross sections
	 *** Output parameter:  photoionization cross section, Mb
	 *** Shell numbers:
	 *** 1 - 1s, 2 - 2s, 3 - 2p, 4 - 3s, 5 - 3p, 6 - 3d, 7 - 4s. 
	 *** If a species in the ground state has no electrons on the given shell,
	 *** the subroutine returns 0.
	 ****************************************************************************** */

	crs = 0.0;
	if( nz < 1 || nz > 30 )
	{ 
		return crs;
	}

	if( ne < 1 || ne > nz )
	{ 
		return crs;
	}

	nout = set_vshell_index( ne, nz );

	if( is > nout )
	{ 
		return crs;
	}

	if( (is == 6 && (nz == 20 || nz == 19)) && ne >= 19 )
	{ 
		return crs;
	}

	ASSERT( is >= 1 && is <= 7 );

	double ethres = getEthresh( is, ne, nz );

	if( e < ethres )
	{
		return crs;
	}

	nint = NINN[ne-1];
	if( ((nz == 15 || nz == 17) || nz == 19) || (nz > 20 && nz != 26) )
	{
		einn = 0.0;
	}
	else
	{
		if( ne < 3 )
		{
			einn = 1.0e30;
		}
		else
		{
			if( nint == nout )
				einn = ethres;
			else
				einn = double(PH1[nint-1][ne-1][nz-1][0]);
		}
	}

	if( (is <= nint || e >= einn) || version == PHFIT95 )
	{
		p1 = -PH1[is-1][ne-1][nz-1][4];
		y = e/PH1[is-1][ne-1][nz-1][1];
		q = -0.5*p1 - L[is-1] - 5.5;
		a = PH1[is-1][ne-1][nz-1][2]*(POW2(y - 1.0) + 
			POW2(PH1[is-1][ne-1][nz-1][5]));
		b = sqrt(y/PH1[is-1][ne-1][nz-1][3]) + 1.0;
		crs = a*pow(y,q)*pow(b,p1);
	}
	else
	{
		if( (is < nout && is > nint) && e < einn )
		{ 
			return crs;
		}
		p1 = -PH2[ne-1][nz-1][3];
		q = -0.5*p1 - 5.5;
		x = e/PH2[ne-1][nz-1][0] - PH2[ne-1][nz-1][5];
		z = sqrt(x*x+POW2(PH2[ne-1][nz-1][6]));
		a = PH2[ne-1][nz-1][1]*(POW2(x - 1.0) + 
			POW2(PH2[ne-1][nz-1][4]));
		b = 1.0 + sqrt(z/PH2[ne-1][nz-1][2]);
		crs = a*pow(z,q)*pow(b,p1);
	}
	return crs;
}

double t_ADfA::hpfit(long int iz,
		     long int n,
		     double e)
{
	long int l, 
	  m;
	double cs,
	  eth, 
	  ex, 
	  q, 
	  x;

	DEBUG_ENTRY( "hpfit()" );

	/*state specific photoionization cross sections for model hydrogen atom
	 * Version 1, September 23, 1997
	 ******************************************************************************
	 *** This subroutine calculates state-specific photoionization cross sections
	 *** for hydrogen and hydrogen-like ions.
	 *** Input parameters:  iz - atomic number 
	 ***          n  - shell number, from 0 to 400:
	 ***                                    0 - 1s
	 ***                                    1 - 2s
	 ***                                    2 - 2p
	 ***                                    3 - 3 
	 ***                                    ......
	 ***          e  - photon energy, eV
	 *** return value - cross section, cm^(-2)     
	 *******************************************************************************/

	ASSERT( iz > 0 && e>0. );

	if( n >= NHYDRO_MAX_LEVEL )
	{ 
		fprintf( ioQQQ, " hpfit called with too large n, =%li\n" , n );
		cdEXIT(EXIT_FAILURE);
	}

	l = 0;
	if( n == 2 )
	{
		l = 1;
	}
	q = 3.5 + l - 0.5*PHH[n][1];

	if( n == 0 )
	{
		m = 1;
	}
	else
	{
		if( n == 1 )
		{
			m = 2;
		}
		else
		{
			m = n;
		}
	}

	eth = getEthresh(1,1,iz)/POW2((double)m);
	ex = MAX2(1. , e/eth );

	/* Don't just force to be at least one...make sure e/eth is close to one or greater.	*/
	ASSERT( e/eth > 0.95 );

	if( ex < 1.0 )
	{ 
		return 0.;
	}

	x = ex/PHH[n][0];
	cs = (PHH[n][4]*pow(1.0 + ((double)PHH[n][2]/x),(double)PHH[n][3])/
	  pow(x,q)/pow(1.0 + sqrt(x),(double)PHH[n][1])*8.79737e-17/
	  POW2((double)iz));
	return cs;
}

void t_ADfA::rec_lines(double t, 
		       realnum r[][NRECCOEFCNO])
{
	long int i, 
	  j, 
	  ipj1, 
	  ipj2;

	double a, 
	  a1, 
	  dr[4][405], 
	  p1, 
	  p2, 
	  p3, 
	  p4, 
	  p5, 
	  p6, 
	  rr[4][110], 
	  te, 
	  x, 
	  z;

	static long jd[6]={143,145,157,360,376,379};

	static long ip[38]={7,9,12,13,14,16,18,19,20,21,22,44,45,49,50,
	  52,53,54,55,56,57,58,59,60,66,67,78,83,84,87,88,95,96,97,100,
	  101,103,104};

	static long id[38]={7,3,10,27,23,49,51,64,38,47,60,103,101,112,
	  120,114,143,145,157,152,169,183,200,163,225,223,237,232,235,
	  249,247,300,276,278,376,360,379,384};

	DEBUG_ENTRY( "rec_lines()" );

	/*effective recombination coefficients for lines of C, N, O, by D. Verner
	 * Version 2, April 30, 1997
	 ******************************************************************************
	 *** This subroutine calculates effective recombination coefficients
	 *** for 110 permitted recombination lines of C, N, O (Pequignot, Petitjean,
	 *** & Boisson, 1991, A&A, 251, 680) and 405 permitted dielectronic 
	 *** recombination lines (Nussbaumer & Storey, 1984, A&AS, 56, 293)
	 *** Input parameter:   t  - temperature, K
	 *** Output parameters: r(i,j), i=1,471
	 ***          r(i,1) - atomic number
	 ***          r(i,2) - number of electrons
	 ***          r(i,3) - wavelength, angstrom
	 ***          r(i,4) - rate coefficient, cm^3 s^(-1)
	 ****************************************************************************** */

	for( i=0; i < 110; i++ )
	{
		rr[0][i] = P[0][i];
		rr[1][i] = P[1][i];
		rr[2][i] = P[2][i];
		z = P[0][i] - P[1][i] + 1.0;
		te = 1.0e-04*t/z/z;
		p1 = P[3][i];
		p2 = P[4][i];
		p3 = P[5][i];
		p4 = P[6][i];
		if( te < 0.004 )
		{
			a1 = p1*pow(0.004,p2)/(1.0 + p3*pow(0.004,p4));
			a = a1/sqrt(te/0.004);
		}
		else
		{
			if( te > 2.0 )
			{
				a1 = p1*pow(2.0,p2)/(1.0 + p3*pow(2.0,p4));
				a = a1/pow(te/2.0,1.5);
			}
			else
			{
				a = p1*pow(te,p2)/(1.0 + p3*pow(te,p4));
			}
		}
		rr[3][i] = 1.0e-13*z*a*P[7][i];
	}

	for( i=0; i < 405; i++ )
	{
		dr[0][i] = ST[0][i];
		dr[1][i] = ST[1][i];
		dr[2][i] = ST[2][i];
		te = 1.0e-04*t;
		p1 = ST[3][i];
		p2 = ST[4][i];
		p3 = ST[5][i];
		p4 = ST[6][i];
		p5 = ST[7][i];
		p6 = ST[8][i];
		if( te < p6 )
		{
			x = p5*(1.0/te - 1.0/p6);
			if( x > 80.0 )
			{
				a = 0.0;
			}
			else
			{
				a1 = (p1/p6 + p2 + p3*p6 + p4*p6*p6)/powpq(p6,3,2)/exp(p5/p6);
				a = a1/exp(x);
			}
		}
		else
		{
			if( te > 6.0 )
			{
				a1 = (p1/6.0 + p2 + p3*6.0 + p4*36.0)/powpq(6.0,3,2)/exp(p5/6.0);
				a = a1/powpq(te/6.0,3,2);
			}
			else
			{
				a = (p1/te + p2 + p3*te + p4*te*te)/powpq(te,3,2)/exp(p5/te);
			}
		}
		dr[3][i] = 1.0e-12*a;
	}

	for( i=0; i < 6; i++ )
	{
		ipj1 = jd[i];
		ipj2 = ipj1 + 1;
		dr[3][ipj1-1] += dr[3][ipj2-1];
		dr[0][ipj2-1] = 0.0;
	}

	for( i=0; i < 38; i++ )
	{
		ipj1 = ip[i];
		ipj2 = id[i];
		rr[3][ipj1-1] += dr[3][ipj2-1];
		dr[0][ipj2-1] = 0.0;
	}

	for( i=0; i < 110; i++ )
	{
		r[0][i] = (realnum)rr[0][i];
		r[1][i] = (realnum)rr[1][i];
		r[2][i] = (realnum)rr[2][i];
		r[3][i] = (realnum)rr[3][i];
	}

	j = 110;
	for( i=0; i < 405; i++ )
	{
		if( dr[0][i] > 1.0 )
		{
			j += 1;
			r[0][j-1] = (realnum)dr[0][i];
			r[1][j-1] = (realnum)dr[1][i];
			r[2][j-1] = (realnum)dr[2][i];
			r[3][j-1] = (realnum)dr[3][i];
		}
	}
	return;
}

double t_ADfA::rad_rec(long int iz, 
		       long int in, 
		       double t)
{
	/*
	 *** Version 4. June 29, 1999.
	 *** Written by D. A. Verner, verner@pa.uky.edu 
	 ******************************************************************************
	 *** This subroutine calculates rates of radiative recombination for all ions
	 *** of all elements from H through Zn by use of the following fits:
	 *** H-like, He-like, Li-like, Na-like - 
	 *** >>refer	all	rec_coef	Verner, D. A., & Ferland, G. J. 1996, ApJS, 103, 467
	 *** Other ions of C, N, O, Ne - Pequignot et al. 1991, A&A, 251, 680,
	 ***    refitted by Verner & Ferland formula to ensure correct asymptotes
	 *** Fe XVII-XXIII - 
	 *** >>refer	Fe17-23	recom	Arnaud, M., & Raymond, J. 1992, ApJ, 398, 394
	 *** Fe I-XV - refitted by Verner & Ferland formula to ensure correct asymptotes
	 *** Other ions of Mg, Si, S, Ar, Ca, Fe, Ni - 
	 ***                      -
	 *** >>refer	all	recom	Shull, M. J., & Van Steenberg, M. 1982, ApJS, 48, 95
	 *** Other ions of Na, Al - 
	 *** >>refer	Na, Al	recom	Landini, M., & Monsignori-Fossi, B.C. 1990, A&AS, 82, 229
	 *** Other ions of F, P, Cl, K, Ti, Cr, Mn, Co (excluding Ti I-II, Cr I-IV,
	 *** Mn I-V, Co I)        - 
	 *** >>refer	many	recom	Landini, M., & Monsignori-Fossi, B.C. 1991, A&AS, 91, 183
	 *** All other species    - interpolations of the power-law fits
	 *** Input parameters:  iz - atomic number 
	 ***                    in - number of electrons from 1 to iz 
	 ***                    t  - temperature, K
	 *** return result:  - rate coefficient, cm^3 s^(-1)
	 ******************************************************************************
	 */
	double tt;
	double rate;

	DEBUG_ENTRY( "rad_rec()" );

	if( iz < 1 || iz > 30 )
	{
		fprintf( ioQQQ, " rad_rec called with insane atomic number, =%4ld\n", 
		  iz );
		cdEXIT(EXIT_FAILURE);
	}
	if( in < 1 || in > iz )
	{
		fprintf( ioQQQ, " rad_rec called with insane number elec =%4ld\n", 
		  in );
		cdEXIT(EXIT_FAILURE);
	}
	if( (((in <= 3 || in == 11) || (iz > 5 && iz < 9)) || iz == 10) || 
	  (iz == 26 && in > 11) )
	{
		tt = sqrt(t/rnew[in-1][iz-1][2]);
		rate = 
		  rnew[in-1][iz-1][0]/(tt*pow(tt + 1.0,1.0 - rnew[in-1][iz-1][1])*
		  pow(1.0 + sqrt(t/rnew[in-1][iz-1][3]),1.0 + rnew[in-1][iz-1][1]));
	}
	else
	{
		tt = t*1.0e-04;
		if( iz == 26 && in <= 13 )
		{
			rate = fe[in-1][0]/pow(tt,fe[in-1][1] + 
			  fe[in-1][2]*log10(tt));
		}
		else
		{
			rate = rrec[in-1][iz-1][0]/pow(tt,(double)rrec[in-1][iz-1][1]);
		}
	}

	return rate;
}

double t_ADfA::H_rad_rec(long int iz,
			 long int n,
			 double t)
{
	/*
	 * Version 4, October 9, 1997
	 ******************************************************************************
	 *** This subroutine calculates state-specific recombination rates 
	 *** for hydrogen and hydrogen-like ions.
	 *** Input parameters:  iz - atomic number 
	 ***          n  - shell number, from 0 to 400:
	 ***                                    0 - 1s
	 ***                                    1 - 2s
	 ***                                    2 - 2p
	 ***                                    3 - 3 
	 ***                                    ......
	 ***          t  - temperature, K
	 *** Output parameter:  r  - rate coefficient, cm^3 s^(-1)
	 *** If n is negative, the subroutine returns the total recombination 
	 *** rate coefficient
	 ******************************************************************************
	 */
	double rate,
	  TeScaled, 
	  x, 
	  x1, 
	  x2;

	DEBUG_ENTRY( "H_rad_rec()" );

	/* iz is charge, must be 1 or greater */
	ASSERT( iz > 0 );

	/* n is level number, must be less than dim or hydro vectors */
	ASSERT( n < NHYDRO_MAX_LEVEL );

	TeScaled = t/POW2((double)iz);

	if( n < 0 )
	{
		x1 = sqrt(TeScaled/3.148);
		x2 = sqrt(TeScaled/7.036e05);
		rate = 7.982e-11/x1/pow(1.0 + x1,0.252)/pow(1.0 + x2,1.748);
	}
	else
	{
		x = log10(TeScaled);
		rate = (HRF[n][0] + HRF[n][2]*x + HRF[n][4]*
		  x*x + HRF[n][6]*powi(x,3) + HRF[n][8]*powi(x,4))/
		  (1.0 + HRF[n][1]*x + HRF[n][3]*x*x + HRF[n][5]*
		  powi(x,3) + HRF[n][7]*powi(x,4));
		rate = exp10(rate)/TeScaled;
	}
	rate *= iz;

	return rate;
}

/*coll_ion D Verner's routine to compute collisional ionization rate coefficients,
 * returns collisional ionization rate coefficient cm^3 s^-1*/
double t_ADfA::coll_ion(
	/* atomic number, 1 for hydrogen */
	long int iz, 
	/* number of bound electrons before ionization*/
	long int in, 
	/* temperature */
	double t)
{
	double rate, te, u;

	DEBUG_ENTRY( "coll_ion()" );
	/*D Verner's routine to compute collisional ionization rate coefficients
	 * Version 3, April 21, 1997
	 * Cu (Z=29) and Zn (Z=30) are added (fits from Ni, correct thresholds).
	 ******************************************************************************
	 *** This subroutine calculates rates of direct collisional ionization 
	 *** for all ionization stages of all elements from H to Ni (Z=28)
	 *** by use of the fits from
	 *>>refer	all	coll_ion	Voronov, G. S. 1997, At. Data Nucl. Data Tables, 65, 1
	 *** Input parameters:  iz - atomic number on pphysical scale, H is 1
	 ***          in - number of electrons from 1 to iz 
	 ***          t  - temperature, K
	 *** Output parameter:  c  - rate coefficient, cm^3 s^(-1)
	 ****************************************************************************** */

	te = t*EVRYD/TE1RYD;
	u = CF[in-1][iz-1][0]/te;
	if( u > 80.0 )
	{ 
		return 0.;
	}

	rate = (CF[in-1][iz-1][2]*(1.0 + CF[in-1][iz-1][1]*
	  sqrt(u))/(CF[in-1][iz-1][3] + u)*pow(u,CF[in-1][iz-1][4])*
	  exp(-u)) * atmdat.lgCollIonOn;

	return rate;
}

double t_ADfA::coll_ion_wrapper(
	/* (atomic number - 1), 0 for hydrogen */
	long int z,
	/* stage of ionization, 0 for atom */
	long int n,
	/* temperature K */
	double t)
{
	double rate = 0.0;

	static const bool DEBUG_COLL_ION = false;

	DEBUG_ENTRY( "coll_ion_wrapper()" );

	if( z < 0 || z > LIMELM-1 )
	{
		/* return zero rate is atomic number outside range of code */
		return 0.;
	}

	if( n < 0 || n > z )
	{
		/* return zero rate is ion stage is impossible */
		return 0.;
	}

	/*Get the rate coefficients from either Dima or Hybrid.
	 * The enum variable CIRCData is set to a default in init_defaults_preparse
	 * It can be changed parse_set via the input command set collisional ionization XXX */
	if( atmdat.CIRCData == t_atmdat::DIMA)
	{
		/* collisional ionization by thermal electrons
		* >>chng 97 mar 19, to Dima's new routine using
		* >>refer	all	coll_ion	Voronov, G. S. 1997, At. Data Nucl. Data Tables, 65, 1
		*coll_ion(atomic number, number of electrons, temperature)*/
		rate = coll_ion( z+1, z+1-n, t );

		if( DEBUG_COLL_ION )
		{
			fprintf(ioQQQ,"DIMA:%li\t%li\t%e\t%e\n",z,n,t,rate);
		}
	}
	else if( atmdat.CIRCData == t_atmdat::HYBRID)
	{
		// Get the rate coefficient by scaling Dima to Dere.
		rate = coll_ion_hybrid( z, n, t );

		if( DEBUG_COLL_ION )
		{
			fprintf(ioQQQ,"HYBRID:%li\t%li\t%e\t%e\n",z,n,t,rate);
		}
	}
	else
	{
		//CIRCData is an enum so it must be equal to one of the enum values.
		TotalInsanity();
	}

	ASSERT( rate >= 0.0 );
	return rate;
}

/*coll_ion D Verner's routine to compute collisional ionization rate coefficients,
 * returns collisional ionization rate coefficient cm^3 s^-1*/
double t_ADfA::coll_ion_hybrid(
	/* (atomic number - 1)  0 for hydrogen */
	long int nelem,
	/* stage of ionization, 0 for atom */
	long int ion,
	/* temperature */
	double t)
{

	DEBUG_ENTRY( "coll_ion_hybrid()" );

	/** DereRatio ratio of Dere 2007 to Dima data, evaluated where
	  ion is abundant
	  */
	static const double DereRatio[30][30]= 
	    	{{0.9063,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{1.0389,1.0686,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.6466,1.0772,0.963,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.9715,0.7079,0.973,0.9899,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{1.0781,1.3336,1.0032,1.1102,0.9874,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{1.0499,0.913,1.0377,1.299,1.2874,1.0656,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{1.0421,1.1966,1.0842,0.9619,1.0583,1.155,1.0648,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{1.041,1.1181,1.0164,1.0271,0.9789,0.6829,1.1805,1.0672,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.2636,0.658,1.0431,1.0749,0.894,1.059,0.9888,1.1426,1.0633,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.8089,1.1395,1.1041,1.0449,1.1365,1.089,1.0147,0.9135,1.336,1.0429,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.9545,1.1302,1.1105,1.2067,1.2569,1.079,0.8866,0.9924,0.9933,1.0488,1.0602,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.3793,0.9857,1.1152,1.2826,1.3006,1.2416,1.0893,0.93,0.9953,0.9992,1.3479,1.0622,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.7458,0.9975,1.0251,0.9553,1.5023,1.2136,1.2566,1.2417,1.2381,1.3755,1.1113,1.1629,1.0589,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.7328,0.8798,0.4492,0.8221,1.7874,1.5418,1.5777,1.3036,1.124,1.0667,1.0009,1.002,1.1284,1.0673,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.7075,0.9805,0.9451,0.5937,1.2061,1.1772,1.3818,1.4073,1.2643,1.1095,0.9903,1.3716,1.0058,1.0966,1.0517,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{1.3572,0.8925,0.8119,0.8991,0.6633,1.4519,1.2073,1.4058,1.4519,1.2726,1.1428,0.9818,1.3643,1.0074,1.1836,1.1005,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.5412,0.8428,0.9237,0.819,0.9151,0.7909,1.7424,1.2653,1.4513,1.4871,1.2633,1.1373,0.9969,0.9919,1.0091,1.2077,1.105,0,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.9242,0.8644,0.9752,0.8562,0.6998,0.6273,0.8686,2.0326,1.0364,1.4687,1.4812,1.2701,1.1376,0.9932,1.0001,1.0084,1.2036,1.1074,0,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.6381,0.9604,1.3556,0.9364,0.9678,1.0604,1.3067,1.0434,2.0722,0.979,1.5318,1.4968,1.2713,1.1427,1.0036,0.9953,1.03,1.2071,1.1083,0,0,0,0,0,0,0,0,0,0,0},
	    	{0.7652,1.1668,1.0422,0.8705,0.9193,0.9982,1.077,1.3875,1.1544,2.4653,1.3188,1.5375,1.5307,1.2776,1.1427,1.0117,0.9872,1.0118,1.1899,1.1119,0,0,0,0,0,0,0,0,0,0},
	    	{1.0951,0.5891,0.9222,1.2903,1.287,1.0563,1.0444,1.1405,1.4841,1.2583,2.7347,0.9844,1.034,1.0412,1.1062,1.2211,1.6728,1.6347,1.6554,1.8095,1.9561,0,0,0,0,0,0,0,0,0},
	    	{0.9789,0.7389,1.3107,0.9274,0.9921,0.9812,0.8971,0.9936,1.0079,1.0377,1.2073,2.7198,0.9995,1.037,1.0433,1.1093,1.2323,1.6673,1.6231,1.6666,1.8089,1.7302,0,0,0,0,0,0,0,0},
	    	{0.6169,0.4618,1.3566,3.3812,1.1951,1.1746,1.0133,0.9195,1.0548,1.0608,1.1016,1.285,3.1081,0.9986,1.0094,1.0379,1.005,0.9932,1.0651,0.8787,0.952,0.9862,1.0083,0,0,0,0,0,0,0},
	    	{1.0603,0.262,0.9237,2.4323,1.8345,1.2197,1.0924,1.0205,0.9379,1.0946,1.1001,1.1741,1.3608,3.152,0.9692,0.8866,1.0556,1.1127,1.2289,1.6828,1.6298,1.6469,1.8082,1.0605,0,0,0,0,0,0},
	    	{0.9061,0.7287,0.9676,1.9425,1.5785,1.9028,1.3827,1.1136,1.0368,0.9516,1.1389,1.1594,1.1642,1.4161,2.8162,0.9744,0.9246,1.0644,1.1118,1.2332,1.6892,1.6311,1.6347,1.8073,1.0722,0,0,0,0,0},
	    	{0.9904,1.0568,1.824,1.6578,1.5033,0.9704,0.8933,0.8579,1.084,1.0308,0.9637,1.1885,1.2179,1.2979,1.4846,3.0344,0.9879,0.8927,1.0534,1.1233,1.2242,1.6794,1.6255,1.6487,1.8141,1.7629,0,0,0,0},
	    	{0.9667,1.2622,0.959,2.8444,1.304,1.5632,1.709,2.298,1.427,1.0885,1.0285,0.9767,1.2237,1.2632,1.3585,1.5495,3.1385,0.9959,0.9327,1.0904,1.0357,1.0125,1.0027,0.9788,0.8975,1.0539,1.048,0,0,0},
	    	{1.0004,0.8721,1.3073,0.9564,2.7856,1.6197,1.5516,2.229,2.1494,1.4916,1.0871,1.0359,0.9768,1.2269,1.3265,1.4169,1.597,3.4077,0.9979,0.8869,1.0635,1.1063,1.2267,1.6914,1.6401,1.6216,1.0598,1.0363,0,0},
	    	{0.728,1.3939,0.4822,0.626,0.5463,2.2491,1.064,1.1998,1.3083,1.6545,1.1149,0.8826,0.8566,0.8196,1.1173,1.17,1.2445,1.394,2.9832,0.8715,0.7703,0.929,0.968,1.0828,1.4973,1.4513,1.4476,0.9621,0.9353,0},
			{1.0766,0.6459,0.7688,0.2358,0.3709,0.3476,1.6754,0.8288,0.9184,1.0686,1.4546,0.9515,0.7272,0.7048,0.6911,0.9722,1.0308,1.0988,1.1959,2.6404,0.7644,0.6768,0.8181,0.8566,0.9689,1.3344,1.3031,1.3345,0.8676,0.8478}};
	/* Get the hybrid coeffs from the Dima rate coefficient times the average Dere to Dima ratio.
	 */
	ASSERT( nelem>=0 && nelem<LIMELM && ion>=0 && ion<=nelem );
	double rate = coll_ion(nelem+1,nelem+1-ion,t) * DereRatio[nelem][ion];
	ASSERT( rate >=0. );
	return rate;
}
const realnum* t_ADfA::h_coll_str(long nHi, long lHi, long nLo, long lLo) const
{
	DEBUG_ENTRY( "h_coll_str()" );

	ASSERT( nLo < nHi );
	ASSERT( nHi <= 5 );

	QNPair inPair(nHi, lHi, 2, -1, nLo, lLo, 2, -1);
	auto p = HCS.find(inPair);
	ASSERT( p != HCS.end() );

	return p->second.data();
}
