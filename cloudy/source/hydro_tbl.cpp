/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "parser.h"
#include "physconst.h"
#include "dense.h"
#include "thirdparty.h"
#include "hydro_tbl.h"

static const char* tpn_file = "hydro_tpn.dat";
static const char* tpnl_file = "hydro_tpnl.dat";
static const char* tpnn_file = "hydro_tpnn.dat";
static const char* pcs_file = "hydro_pcs.dat";

static const long tpn_magic = 201809260L;
static const long tpnl_magic = 201809271L;
static const long tpnn_magic = 201809281L;
static const long pcs_magic = 201810020L;

// p_stride needs a definition, following the discussion under "Constant static members" here:
//
// https://en.cppreference.com/w/cpp/language/static
//
// p_stride is odr-used when an address is taken in min(), hence a definition at namespace scope
// is required, but it cannot have an initializer here, that should be in the class.
const size_t t_hydro_tbl::p_stride;

void t_hydro_tbl::p_initn(long nu)
{
	DEBUG_ENTRY( "t_hydro_tbl::p_initn()" );

	if( p_nmaxn == 0 )
	{
		DataParser d(tpn_file, ES_NONE);
		d.getline();
		d.checkMagic(tpn_magic);

		// read dimension of the table
		d.getline();
		d.getToken(p_nmaxn);

		p_nmaxn_read = 0;
		p_posn = d.getpos();
	}

	if( nu < long(p_nmaxn_read) )
		return;

	ASSERT( nu <= long(p_nmaxn) );

	while( nu > long(p_nmaxn_read) )
	{
		DataParser d(tpn_file, ES_NONE);
		d.setpos(p_posn);

		// allocate space for data
		auto tpn = new arrn;
		size_t bsize = min(p_stride, p_nmaxn-p_nmaxn_read);
		tpn->reserve(bsize);
		for( size_t i=0; i < bsize; ++i )
		{
			size_t lsize = p_nmaxn_read+i;
			if( lsize > 0 )
				tpn->reserve(i, lsize);
		}
		tpn->alloc();
		p_tpn.push_back(tpn);

		// now read the TP data
		for( size_t k=0; k < bsize; ++k )
		{
			size_t nu = p_nmaxn_read + k + 1;
			for( size_t nl=1; nl < nu; ++nl )
			{
				d.getline();
				size_t n;
				d.getToken(n);
				if( n != nl )
					d.errorAbort("invalid quantum number");
				d.getToken(n);
				if( n != nu )
					d.errorAbort("invalid quantum number");
				d.getToken((*tpn)[k][nl-1]);
				d.checkEOL();
			}
		}

		p_nmaxn_read += bsize;
		p_posn = d.getpos();
	}
}

void t_hydro_tbl::p_initnl(long nu)
{
	DEBUG_ENTRY( "t_hydro_tbl::p_initnl()" );

	if( p_nmaxnl_u == 0 )
	{
		DataParser d(tpnl_file, ES_NONE);
		d.getline();
		d.checkMagic(tpnl_magic);

		// read dimension of the table
		d.getline();
		d.getToken(p_nmaxnl_u);
		d.getline();
		d.getToken(p_nmaxnl_l);

		p_nmaxnl_read = 0;
		p_posnl = d.getpos();
	}

	if( nu < long(p_nmaxn_read) )
		return;

	ASSERT( nu <= long(p_nmaxnl_u) );

	while( nu > long(p_nmaxnl_read) )
	{
		DataParser d(tpnl_file, ES_NONE);
		d.setpos(p_posnl);

		// allocate space for data
		auto tpnl = new arrnl;
		size_t bsize = min(p_stride, p_nmaxnl_u-p_nmaxnl_read);
		tpnl->reserve(bsize);
		for( size_t i=0; i < bsize; ++i )
		{
			size_t im = min(p_nmaxnl_read+i,p_nmaxnl_l);
			if( im > 0 )
				tpnl->reserve(i, im);
			for( size_t j=0; j < im; ++j )
			{
				tpnl->reserve(i, j, im+1);
				for( size_t k=0; k < im+1; ++k )
					tpnl->reserve(i, j, k, 2);
			}
		}
		tpnl->alloc();
		p_tpnl.push_back(tpnl);

		// now read the TP data
		for( size_t k=0; k < bsize; ++k )
		{
			size_t nu = p_nmaxnl_read + k + 1;
			for( size_t nl=1; nl < min(nu,p_nmaxnl_l+1); ++nl )
			{
				for( size_t lu = 0; lu < nu; ++lu )
				{
					for( size_t dl=0; dl < 2; ++dl )
					{
						int ll = int(lu) - 1 + 2*dl;
						if( ll >= 0 && ll < int(nl) )
						{
							d.getline();
							size_t n;
							d.getToken(n);
							if( n != nl )
								d.errorAbort("invalid quantum number");
							d.getToken(n);
							if( n != size_t(ll) )
								d.errorAbort("invalid quantum number");
							d.getToken(n);
							if( n != nu )
								d.errorAbort("invalid quantum number");
							d.getToken(n);
							if( n != lu )
								d.errorAbort("invalid quantum number");
							d.getToken((*tpnl)[k][nl-1][lu][dl]);
							d.checkEOL();
						}
					}
				}
			}
		}

		p_nmaxnl_read += bsize;
		p_posnl = d.getpos();
	}
}

void t_hydro_tbl::p_initnn()
{
	DEBUG_ENTRY( "t_hydro_tbl::p_initnn()" );

	if( p_tpnn.size() > 0 )
		return;

	DataParser d(tpnn_file, ES_NONE);
	d.getline();
	d.checkMagic(tpnn_magic);

	// read dimension of the table
	d.getline();
	d.getToken(p_Zmax);
	d.getline();
	d.getToken(p_nmaxnn);

	// allocate space for data
	p_tpnn.reserve(p_Zmax);
	for( size_t i=0; i < p_Zmax; ++i )
	{
		p_tpnn.reserve(i, p_nmaxnn);
		for( size_t j=0; j < p_nmaxnn; ++j )
		{
			if( j > 0 )
				p_tpnn.reserve(i, j, j);
		}
	}
	p_tpnn.alloc();
	p_wnnn.alloc( p_tpnn.clone() );

	// now read the TP data
	for( long Z=1; Z <= min(p_Zmax,LIMELM); ++Z )
	{
		d.getline();
		long ZZ;
		d.getToken(ZZ);
		if( Z != ZZ )
			d.errorAbort("invalid atomic number");

		for( size_t n=2; n <= 100; ++n )
		{
			for( size_t ll=0; ll < n-1; ++ll )
			{
				d.getline();
				size_t k;
				d.getToken(k);
				if( k != n )
					d.errorAbort("invalid quantum number");
				d.getToken(k);
				if( k != ll )
					d.errorAbort("invalid quantum number");
				d.getToken(k);
				if( k != ll+1 )
					d.errorAbort("invalid quantum number");
				d.getToken(p_tpnn[Z-1][n-1][ll]);
				d.getToken(p_wnnn[Z-1][n-1][ll]);
				d.checkEOL();
			}
		}
	}
}

void t_hydro_tbl::p_initcs()
{
	DEBUG_ENTRY( "t_hydro_tbl::p_initcs()" );

	if( p_en.size() > 0 )
		return;

	DataParser d(pcs_file, ES_NONE);
	d.getline();
	d.checkMagic(pcs_magic);

	// read dimension of the table
	d.getline();
	d.getToken(p_nmaxcs);
	d.getline();
	d.getToken(p_nenrgs);

	// allocate space for data
	p_en.reserve(p_nmaxcs);
	for( size_t i=0; i < p_nmaxcs; ++i )
	{
		p_en.reserve(i, i+1);
		for( size_t j=0; j < i+1; ++j )
		{
			p_en.reserve(i, j, p_nenrgs);
		}
	}
	p_en.alloc();
	p_cs.alloc( p_en.clone() );

	// now read the CS data
	for( size_t n=1; n <= p_nmaxcs; ++n )
	{
		for( size_t l=0; l < n; ++l )
		{
			for( size_t k=0; k < p_nenrgs; ++k )
			{
				d.getline();
				size_t m;
				d.getToken(m);
				if( m != n )
					d.errorAbort("invalid atomic number");
				d.getToken(m);
				if( m != l )
					d.errorAbort("invalid quantum number");
				d.getToken(p_en[n-1][l][k]);
				d.getToken(p_cs[n-1][l][k]);
				d.checkEOL();
			}
		}
	}
}

realnum t_hydro_tbl::p_RM(long Z) const
{
	DEBUG_ENTRY( "t_hydro_tbl::p_RM()" );

	// this is an approximation for the mass of the nucleus...
	double mass_nuc = ( Z == 1 ) ? PROTON_MASS : ATOMIC_MASS_UNIT * dense.AtomicWeight[Z-1];
	return realnum(RYD_INF/(1. + ELECTRON_MASS/mass_nuc));
}

realnum t_hydro_tbl::tp(long nl, long nu, long Z)
{
	DEBUG_ENTRY( "t_hydro_tbl::tp()" );

	p_initn(nu);

	ASSERT( nl > 0 && nl < nu && Z > 0 && Z <= LIMELM );

	size_t k1 = (nu-1)/p_stride;
	size_t k2 = (nu-1)%p_stride;
	return powi(realnum(Z),4)*(*p_tpn[k1])[k2][nl-1]*p_RM(Z);
}

realnum t_hydro_tbl::tp(long nl, long ll, long nu, long lu, long Z)
{
	DEBUG_ENTRY( "t_hydro_tbl::tp()" );

	p_initnl(nu);

	ASSERT( nl > 0 && nl < nu && nl <= long(p_nmaxnl_l) && Z > 0 && Z <= LIMELM );
	ASSERT( ll >= 0 && ll < nl && lu >= 0 && lu < nu && abs(ll-lu) == 1 );

	size_t k1 = (nu-1)/p_stride;
	size_t k2 = (nu-1)%p_stride;
	long dl = ( ll < lu ) ? 0 : 1;
	return powi(realnum(Z),4)*(*p_tpnl[k1])[k2][nl-1][lu][dl]*p_RM(Z);
}

#ifdef NDEBUG
realnum t_hydro_tbl::tp(long n, long ll, long, long Z)
#else
realnum t_hydro_tbl::tp(long n, long ll, long lu, long Z)
#endif
{
	DEBUG_ENTRY( "t_hydro_tbl::tp()" );

	p_initnn();

	ASSERT( n > 0 && n <= long(p_nmaxnn) && Z > 0 && Z <= LIMELM );
	ASSERT( ll >= 0 && ll == lu - 1 && lu < n );

	return p_tpnn[Z-1][n-1][ll];
}

#ifdef NDEBUG
double t_hydro_tbl::wn(long n, long ll, long, long Z)
#else
double t_hydro_tbl::wn(long n, long ll, long lu, long Z)
#endif
{
	DEBUG_ENTRY( "t_hydro_tbl::wn()" );

	p_initnn();

	ASSERT( n > 0 && n <= long(p_nmaxnn) && Z > 0 && Z <= LIMELM );
	ASSERT( ll >= 0 && ll == lu - 1 && lu < n );

	return p_wnnn[Z-1][n-1][ll];
}

double t_hydro_tbl::cs(double e, long n, long l, long Z)
{
	DEBUG_ENTRY( "t_hydro_tbl::cs()" );

	p_initcs();

	ASSERT( n > 0 && n <= long(p_nmaxcs) && l >= 0 && l < n && Z > 0 && Z <= LIMELM );

	if( e < 1. )
		return 0.;
	else
	{
		double elog = max(log(e), 2.*FLT_EPSILON);
		// the tables should extend beyond the highest energy in the Cloudy mesh
		if( elog > p_en.back(n-1,l) )
			TotalInsanity();
		long ilo = min(max(hunt_bisect(&p_en[n-1][l][0], p_nenrgs, elog)-1,0),p_nenrgs-4);
		double cslog = lagrange(&p_en[n-1][l][ilo], &p_cs[n-1][l][ilo], 4, elog);
		return exp(cslog)/double(Z*Z);
	}
}
