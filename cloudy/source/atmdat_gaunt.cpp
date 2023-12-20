/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "phycon.h"
#include "rfield.h"
#include "opacity.h"
#include "iso.h"
#include "dense.h"
#include "mole.h"
#include "vectorize.h"
#include "atmdat_gaunt.h"

// use 3rd-order interpolation
static const size_t N_GFF_INTERP = 4;

static const char* gauntff_file = "gauntff_merged_Zxx.dat";

static const long gaunt_magic = 20140510L;

void t_gaunt::p_read_table(const char* fnam, long nelem)
{
	DEBUG_ENTRY( "t_gaunt::p_read_table()" );

	if( p_gff.size() > 0 && !isnan(p_gff[nelem][0][0][0]) )
		return;

	// generate file name
	string fnam1 = fnam;
	ostringstream oss;
	oss << setw(2) << setfill('0') << nelem+1;
	string::size_type ptr = fnam1.find( "xx" );
	ASSERT( ptr != string::npos );
	fnam1.replace( ptr, 2, oss.str() );

	fstream io;
	open_data( io, fnam1, mode_r );
	string line;
	// read magic number
	getline( io, line );
	istringstream iss( line );
	long magic;
	iss >> magic;
	if( magic != gaunt_magic )
	{
		fprintf( ioQQQ, " t_gaunt::p_read_table() found wrong magic number in file %s.\n", fnam );
		fprintf( ioQQQ, " I expected to find version: %ld\n", gaunt_magic );
		cdEXIT(EXIT_FAILURE);
	}
	// read dimensions of the table
	getline( io, line );
	iss.str(line);
	iss >> p_np_gam2 >> p_np_u;
	ASSERT( p_np_gam2 >= N_GFF_INTERP && p_np_u >= N_GFF_INTERP );
	// read start value for log(gamma^2)
	getline( io, line );
	iss.str(line);
	iss >> p_lg_gam2_min;
	// read start value for log(u)
	getline( io, line );
	iss.str(line);
	iss >> p_lg_u_min;
	// read step size in dex
	getline( io, line );
	iss.str(line);
	iss >> p_step;
	// read atomic number
	getline( io, line );
	iss.str(line);
	int Z;
	iss >> Z;
	ASSERT( Z == nelem+1 );
	
	if( p_gff.size() == 0 )
	{
		p_gff.alloc( LIMELM, N_GFF_INTERP, p_np_u, p_np_gam2 );
		// set to NaN to avoid uninitialized use
		p_gff.invalidate();

		p_lg_gam2.resize( p_np_gam2 );
		for( size_t ipg2=0; ipg2 < p_np_gam2; ++ipg2 )
			p_lg_gam2[ipg2] = p_lg_gam2_min + double(ipg2)*p_step;
		p_lg_gam2_max = p_lg_gam2.back();

		p_lg_u.resize( p_np_u );
		for( size_t ipu=0; ipu < p_np_u; ++ipu )
			p_lg_u[ipu] = p_lg_u_min + double(ipu)*p_step;
		p_lg_u_max = p_lg_u.back();
	}
	else
	{
		// check that table parameters are the same
		ASSERT( p_np_gam2 == p_lg_gam2.size() );
		ASSERT( p_np_u == p_lg_u.size() );
		ASSERT( fp_equal( p_lg_gam2_min, p_lg_gam2[0] ) );
		ASSERT( fp_equal( p_lg_u_min, p_lg_u[0] ) );
		ASSERT( fp_equal( p_step, p_lg_u[1]-p_lg_u[0], 2*int(abs(p_lg_u[0]/p_step)) ) );
	}

	// next two lines are comments
	getline( io, line );
	getline( io, line );

	for( size_t ipu=0; ipu < p_np_u; ++ipu )
	{
		getline( io, line );
		iss.str(line);
		for( size_t ipg2=0; ipg2 < p_np_gam2; ++ipg2 )
		{
			double val;
			iss >> val;
			p_gff[nelem][0][ipu][ipg2] = log(val);
		}
		// reset eof flag, otherwise subsequent reads would fail
		iss.clear( iss.rdstate() & ~iss.eofbit );

		// calculate higher derivatives for Newton interpolation
		for( size_t i=1; i < N_GFF_INTERP; ++i )
			for( size_t ipg2=0; ipg2 < p_np_gam2-i; ++ipg2 )
				p_gff[nelem][i][ipu][ipg2] =
					(p_gff[nelem][i-1][ipu][ipg2+1]-p_gff[nelem][i-1][ipu][ipg2])/
					(double(i)*p_step);
	}

	// the remainder of the file contains the uncertainties of the gff values
	// we will not read those here...

	if( io.fail() || iss.fail() )
	{
		fprintf( ioQQQ, " An error occurred while reading the file %s. Bailing out.\n", fnam );
		cdEXIT(EXIT_FAILURE);
	}
}

double t_gaunt::gauntff(long Z, double Te, double anu)
{
	DEBUG_ENTRY( "t_gaunt::gauntff()" );

	ASSERT( Z > 0 && Z <= LIMELM );

	p_read_table(gauntff_file, Z-1);

	double gam2 = pow2(double(Z))*TE1RYD/Te;
	double u = anu*TE1RYD/Te;

	ASSERT( gam2 > 0. && u > 0. );

	double lg_gam2 = log10(gam2);
	double lg_u = log10(u);

	ASSERT( p_lg_gam2_min <= lg_gam2 && lg_gam2 <= p_lg_gam2_max );
	ASSERT( p_lg_u_min <= lg_u && lg_u <= p_lg_u_max );

	long ipg2 = min(max(size_t((lg_gam2-p_lg_gam2_min)/p_step),1)-1,p_np_gam2-N_GFF_INTERP);
	long ipu = min(max(size_t((lg_u-p_lg_u_min)/p_step),1)-1,p_np_u-N_GFF_INTERP);

	double interp[N_GFF_INTERP];

	for( size_t i=0; i < N_GFF_INTERP; ++i )
		interp[i] = lagrange(get_ptr(&p_lg_gam2[ipg2]), get_ptr(&p_gff[Z-1][0][ipu+i][ipg2]),
				     N_GFF_INTERP, lg_gam2);

	return exp(lagrange(get_ptr(&p_lg_u[ipu]), interp, N_GFF_INTERP, lg_u));
}

void t_gaunt::p_gauntff_vec(long Z, double Te, const double anulog10[], long nflux)
{
	DEBUG_ENTRY( "t_gaunt::p_gauntff_vec()" );

	p_read_table(gauntff_file, Z-1);

	if( fp_equal( Te, p_gff_ion[Z].Te_used )
	    && p_gff_ion[Z].nflux_used >= nflux )
		return;

	if( p_vexp_arg.size() == 0 )
		p_vexp_arg.resize( rfield.nflux_with_check );

	if( p_gff_ion[Z].gff.size() == 0 )
	{
		p_gff_ion[Z].gff.resize( rfield.nflux_with_check );
		set_NaN( get_ptr(p_gff_ion[Z].gff), rfield.nflux_with_check );
	}

	if( !fp_equal( Te, p_gff_ion[Z].Te_used ) )
		p_gauntff_vec_sub( Z, Te, anulog10, 0, nflux );
	else if( p_gff_ion[Z].nflux_used < nflux )
		p_gauntff_vec_sub( Z, Te, anulog10, p_gff_ion[Z].nflux_used, nflux );

	p_gff_ion[Z].Te_used = Te;
	p_gff_ion[Z].nflux_used = nflux;
}

void t_gaunt::p_gauntff_vec_sub(long Z, double Te, const double anulog10[], long nmin, long nmax)
{
	DEBUG_ENTRY( "t_gaunt::p_gauntff_vec_sub()" );

	ASSERT( Z > 0 && Z <= LIMELM );

	vector_avx<realnum>& gff = p_gff_ion[Z].gff;

	double lg_gam2 = log10(pow2(double(Z))*TE1RYD/Te);
	ASSERT( p_lg_gam2_min <= lg_gam2 && lg_gam2 <= p_lg_gam2_max );
	long ipg2 = min(max(size_t((lg_gam2-p_lg_gam2_min)/p_step),1)-1,p_np_gam2-N_GFF_INTERP);

	double lg_u0 = log10(TE1RYD/Te);
	ASSERT( p_lg_u_min <= lg_u0+anulog10[nmin] && lg_u0+anulog10[nmax-1] <= p_lg_u_max );
	size_t ipumin = min(max(size_t((lg_u0+anulog10[nmin]-p_lg_u_min)/p_step),1)-1,p_np_u-N_GFF_INTERP);
	size_t ipumax = min(max(size_t((lg_u0+anulog10[nmax-1]-p_lg_u_min)/p_step),1)-1,p_np_u-N_GFF_INTERP)+N_GFF_INTERP;

	multi_arr<double,2,C_TYPE> interp( N_GFF_INTERP, p_np_u );
	for( size_t ipu=ipumin; ipu < ipumax; ++ipu )
	{
		// Newton interpolation in the gamma^2 direction
		size_t i = 1;
		const double* plu = get_ptr(&p_lg_gam2[ipg2]);
		double val = p_gff[Z-1][0][ipu][ipg2];
		double fac = lg_gam2 - *plu++;
		while( true )
		{
			val += fac*p_gff[Z-1][i][ipu][ipg2];
			if( ++i == N_GFF_INTERP )
				break;
			fac *= lg_gam2 - *plu++;
		}
		interp[0][ipu] = val;
	}

	// fill in higher derivatives for Newton interpolation
	for( size_t i=1; i < N_GFF_INTERP; ++i )
		for( size_t ipu=ipumin; ipu < ipumax-i; ++ipu )
			interp[i][ipu] = (interp[i-1][ipu+1]-interp[i-1][ipu])/(double(i)*p_step);

	size_t ipu = ipumin;

	// this loop burns lots of CPU time, so it should be well optimized
	for( long j=nmin; j < nmax; ++j )
	{
		double lg_u = lg_u0 + anulog10[j];
		while( lg_u >= p_lg_u[ipu+2] )
			ipu++;
		// Newton interpolation in the u direction
		size_t i = 1;
		const double* plu = get_ptr(&p_lg_u[ipu]);
		const double* pin = get_ptr(&interp[0][ipu]);
		double gval = *pin;
		double fac = lg_u - *plu++;
		while( true )
		{
			pin += p_np_u;
			gval += fac*(*pin);
			if( ++i == N_GFF_INTERP )
				break;
			fac *= lg_u - *plu++;
		}
		p_vexp_arg[j] = realnum(gval);
	}
	vexp( get_ptr(p_vexp_arg), get_ptr(gff), nmin, nmax );
}

void t_gaunt::p_setup_brems(long ion, double Te)
{
	DEBUG_ENTRY( "t_gaunt::p_setup_brems()" );

	long limit = min( rfield.ipMaxBolt, rfield.nflux );

	if( ion > 0 && ion < LIMELM+1 )
	{
		if( fp_equal( Te, p_cache_vec[ion].Te_used )
		    && p_cache_vec[ion].nhi >= limit )
			return;

		if( p_cache_vec[ion].brems_vec.size() == 0 )
			p_cache_vec[ion].brems_vec.resize( rfield.nflux_with_check );

		p_gauntff_vec( ion, Te, rfield.anulog10ptr(), limit );

		p_cache_vec[ion].Te_used = Te;
		p_cache_vec[ion].nhi = limit;
		double ion2 = pow2(double(ion));

		for( long i=0; i < p_cache_vec[ion].nhi; ++i )
			p_cache_vec[ion].brems_vec[i] =
				ion2*p_gff_ion[ion].gff[i]*rfield.widflx(i)*rfield.ContBoltz[i];
	}
	else if( ion == -1 )
	{
		// special case for H-
		if( fp_equal( Te, p_hminus_vec.Te_used )
		    && p_hminus_vec.nhi >= limit )
			return;

		if( p_hminus_vec.brems_vec.size() == 0 )
			p_hminus_vec.brems_vec.resize( rfield.nflux_with_check );

		p_gauntff_vec( 1, Te, rfield.anulog10ptr(), limit );

		p_hminus_vec.Te_used = Te;
		p_hminus_vec.nhi = limit;

		for( long i=0; i < p_hminus_vec.nhi; ++i )
			p_hminus_vec.brems_vec[i] =
				p_gff_ion[1].gff[i]*opac.OpacStack[i-1+opac.iphmra]*
				rfield.widflx(i)*rfield.ContBoltz[i];
	}
	else
	{
		TotalInsanity();
	}
}

double t_gaunt::brems_cool(long ion, double Te)
{
	DEBUG_ENTRY( "t_gaunt::brems_cool()" );

	ASSERT( Te > 0. );

	long limit = min( rfield.ipMaxBolt, rfield.nflux );

	if( ion > 0 && ion < LIMELM+1 )
	{
		if( fp_equal( Te, p_cache[ion].Te_used )
		    && p_cache[ion].nlo == rfield.ipEnergyBremsThin
		    && p_cache[ion].nhi == limit )
			return p_cache[ion].brems_sum;

		p_setup_brems( ion, Te );

		p_cache[ion].Te_used = Te;
		p_cache[ion].nlo = rfield.ipEnergyBremsThin;
		p_cache[ion].nhi = limit;

		p_cache[ion].brems_sum = 0.;
		for( long i=p_cache[ion].nlo; i < p_cache[ion].nhi; ++i )
			p_cache[ion].brems_sum += p_cache_vec[ion].brems_vec[i];

		return p_cache[ion].brems_sum;
	}
	else if( ion == -1 )
	{
		// special case for H-
		if( fp_equal( Te, p_hminus.Te_used )
		    && p_hminus.nlo == rfield.ipEnergyBremsThin
		    && p_hminus.nhi == limit )
			return p_hminus.brems_sum;

		p_setup_brems( ion, Te );

		p_hminus.Te_used = Te;
		p_hminus.nlo = rfield.ipEnergyBremsThin;
		p_hminus.nhi = limit;

		p_hminus.brems_sum = 0.;
		for( long i=p_hminus.nlo; i < p_hminus.nhi; ++i )
			// OpacStack contains the ratio of the H- to H brems cross section
			p_hminus.brems_sum += p_hminus_vec.brems_vec[i];

		return p_hminus.brems_sum;
	}
	else
	{
		TotalInsanity();
	}
}

void t_gaunt::brems_rt(long ion, double Te, double abun, vector<double>& arr)
{
	DEBUG_ENTRY( "t_gaunt::brems_rt()" );

	ASSERT( Te > 0. );

	size_t limit = min( rfield.ipMaxBolt, rfield.nflux );

	if( ion > 0 && ion < LIMELM+1 )
	{
		ASSERT( arr.size() >= limit );
		p_setup_brems( ion, Te );
		for( size_t i=0; i < limit; ++i )
			arr[i] += abun*p_cache_vec[ion].brems_vec[i];
	}
	else if( ion == -1 )
	{
		ASSERT( arr.size() >= limit );
		p_setup_brems( ion, Te );
		for( size_t i=0; i < limit; ++i )
			arr[i] += abun*p_hminus_vec.brems_vec[i];
	}
	else
	{
		TotalInsanity();
	}
}

void t_gaunt::brems_opac(long ion, double Te, double abun, vector<double>& arr)
{
	DEBUG_ENTRY( "t_gaunt::brems_opac()" );

	ASSERT( Te > 0. );

	if( ion > 0 && ion < LIMELM+1 )
	{
		p_gauntff_vec( ion, Te, rfield.anulog10ptr(), rfield.nflux );
		for( long i=0; i < rfield.nflux; ++i )
			arr[i] += abun*pow2(double(ion))*p_gff_ion[ion].gff[i];
	}
	else if( ion == -1 )
	{
		p_gauntff_vec( 1, Te, rfield.anulog10ptr(), rfield.nflux );
		for( long i=0; i < rfield.nflux; ++i )
			arr[i] += abun*p_gff_ion[1].gff[i]*opac.OpacStack[i-1+opac.iphmra];
	}
	else
	{
		TotalInsanity();
	}
}

void t_gaunt::brems_sum_ions(t_brems_den& sum) const
{
	DEBUG_ENTRY( "t_gaunt::brems_sum_ions()" );

	sum.den_Hm = iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop();
	sum.den_Hp = dense.xIonDense[ipHYDROGEN][1];
	sum.den_Hp += findspecieslocal("H2+")->den;
	sum.den_Hp += findspecieslocal("H3+")->den;
	sum.den_Hep = dense.xIonDense[ipHELIUM][1];
	sum.den_Hepp = dense.xIonDense[ipHELIUM][2];
	sum.den_ion[0] = 0.;
	for( long ion=1; ion < LIMELM+1; ++ion )
	{
		sum.den_ion[ion] = 0.;
		for( long nelem=ipLITHIUM; nelem < LIMELM; ++nelem )
			if( dense.lgElmtOn[nelem] && ion <= nelem+1 )
				sum.den_ion[ion] += dense.xIonDense[nelem][ion];
	}
	/* add molecular ions */
	for( long ipMol = 0; ipMol < mole_global.num_calc; ipMol++ )
	{
		if( !mole_global.list[ipMol]->isMonatomic() && 
		    mole_global.list[ipMol]->isIsotopicTotalSpecies() &&
		    mole_global.list[ipMol]->charge > 0 &&
		    mole_global.list[ipMol]->label != "H2+" &&
		    mole_global.list[ipMol]->label != "H3+" )
		{	
			ASSERT( mole_global.list[ipMol]->charge < LIMELM+1 );
			sum.den_ion[mole_global.list[ipMol]->charge] += mole.species[ipMol].den;
		}
	}
}
