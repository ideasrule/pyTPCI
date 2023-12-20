/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "trace.h"
#include "iso.h"
#include "freebound.h"
#include "opacity.h"
#include "dense.h"
#include "mesh.h"

void t_mesh::p_ReadResolution()
{
	DEBUG_ENTRY( "p_ReadResolution()" );

	if( trace.lgTrace )
		fprintf( ioQQQ,"p_ReadResolution opening continuum_mesh.ini:");

	FILE* ioDATA = open_data( "continuum_mesh.ini", "r" );

	string chLine;
	/* check that magic number is ok */
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, "p_ReadResolution could not read first line of continuum_mesh.ini.\n");
		cdEXIT(EXIT_FAILURE);
	}

	long i = 1;
	bool lgEOL;
	/* continuum mesh magic number */
	long i1 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long i2 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long i3 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

	bool lgResPower;

	/* the following is the set of numbers that appear at the start of continuum_mesh.ini */
	if( i1 == 1 && i2 == 9 && i3 == 29 )
		// old version of the file (c08 and older), this has pairs: upper limit freq range, resolution
		// this format is still supported to accomodate users with existing continuum_mesh.ini files.
		lgResPower = false;
	else if( i1 == 10 && i2 == 8 && i3 == 8 )
		// new version of the file (c10 and newer), this has pairs: upper limit freq range, resolving power
		// resolving power = 1./resolution
		lgResPower = true;
	else
	{
		fprintf( ioQQQ, 
			"p_ReadResolution: the version of continuum_mesh.ini is not supported.\n" );
		fprintf( ioQQQ, 
			"I found version number %li %li %li.\n" ,
			 i1 , i2 , i3 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine.c_str() );
		cdEXIT(EXIT_FAILURE);
	}

	size_t n = 0;
	while( read_whole_line( chLine, ioDATA ) )
	{
		/* skip comment lines */
		if( chLine[0] != '#')
		{
			i = 1;
			double upper_limit = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
			double val2 = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

			// continuum energy could be 0 to indicate high energy bound of the code
			if( upper_limit == 0. )
				upper_limit = p_egamry;

			// all must be positive
			if( upper_limit <= 0. || val2 <= 0. )
			{
				fprintf(ioQQQ, "DISASTER PROBLEM continuum_mesh.ini has a non-positive number.\n");
				cdEXIT(EXIT_FAILURE);
			}

			// convert resolving power into resolution
			double resolution;
			if( lgResPower )
				// resolving power was entered
				resolution = 1./val2;
			else
				// resolution was entered
				resolution = val2;

			/* this is option to rescale resolution with set resolution command */
			resolution *= p_ResolutionScaleFactor;

			// merge in the major ionization edges
			while( true )
			{
				// if a user-supplied mesh edge agrees within 1/10th of a resolution
				// element with an edge from p_edges[], we will merge and use the latter
				double toler = upper_limit*resolution*0.1;
				if( n < p_edges.size()
				    && fp_equal_tol( p_edges[n].Ryd(), upper_limit, toler ) )
				{
					p_RangeUpperLimit.push_back(p_edges[n].Ryd());
					p_RangeResolution.push_back(resolution);
					++n;
					break;
				}
				else if( n < p_edges.size()
					 && p_edges[n].Ryd() < upper_limit )
				{
					p_RangeUpperLimit.push_back(p_edges[n].Ryd());
					p_RangeResolution.push_back(resolution);
					++n;
				}
				else
				{
					p_RangeUpperLimit.push_back(upper_limit);
					p_RangeResolution.push_back(resolution);
					break;
				}
			}
		}
	}

	ASSERT( n == p_edges.size() );

	fclose( ioDATA );

	/* now verify continuum grid is ok - first are all values but the last positive? */
	for( n=1; n < p_RangeUpperLimit.size(); ++n )
	{
		if( p_RangeUpperLimit[n-1] >= p_RangeUpperLimit[n] )
		{
			fprintf( ioQQQ, 
				"p_ReadResolution: the range upper limits must be in increasing order.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	if( !fp_equal(p_RangeUpperLimit.back(), p_egamry) )
	{
		fprintf( ioQQQ, 
			"p_ReadResolution: the last range upper limit must be zero.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return;
}

void t_mesh::p_SetupMesh(bool lgUnitCell)
{
	DEBUG_ENTRY( "p_SetupMesh()" );

	if( trace.lgTrace )
		fprintf( ioQQQ, "p_SetupMesh() setting up the frequency mesh. ResolutionScaleFactor = %g.\n",
			 p_ResolutionScaleFactor );

	ASSERT( p_RangeUpperLimit.size() > 0 && p_RangeUpperLimit.size() == p_RangeResolution.size() );

	// reserve some space for the mesh, no problem if we overrun this though...
	p_anu_edge.reserve( 8500 );

	// fill in the edges of the frequency cells, these will have a logarithmic spacing
	// in each of the frequency ranges with a resolving power set by the user.
	//
	// the lower bound of the first cell is exactly equal to emm()
	// the upper bound of the last cell is exactly equal to egamry()
	//
	// the anu and widflx arrays will be derived from the cell boundaries
	// widflx[i] will simply be the difference between the upper and lower bound
	// and anu[i] will be exactly in the middle between the two bounds.
	//
	// lgUnitCell == false: use highest cell of regular mesh for unit tests
	//                      the frequency range available for regular RT will NOT reach egamry().
	// lgUnitCell == true: create additional cell above regular mesh for unit tests
	//                     the frequency range available for regular RT will reach egamry().
	//                     this choice currently trips at least one ASSERT.
	double hi_bound = p_emm;
	p_anu_edge.push_back( hi_bound );
	for( size_t n=0; n < p_RangeUpperLimit.size(); ++n )
	{
		double lo_bound = hi_bound;
		hi_bound = p_RangeUpperLimit[n];
		size_t nbin = size_t(ceil(log(hi_bound/lo_bound)/p_RangeResolution[n]));
		double bound = lo_bound;
		double step = exp(log(hi_bound/lo_bound)/nbin);
		for( size_t i=1; i < nbin; ++i )
		{
			bound *= step;
			p_anu_edge.push_back( bound );
		}
		p_anu_edge.push_back( hi_bound );
	}

	if( lgUnitCell )
	{
		// create an extra cell above the upper limit for the unit test
		double edge = pow2(p_anu_edge[p_anu_edge.size()-1])/p_anu_edge[p_anu_edge.size()-2];
		p_anu_edge.push_back( edge );
	}

	// anu[i] should be exactly in the middle of p_anu_edge[i] and p_anu_edge[i+1]
	p_anu.resize( p_anu_edge.size()-1 );
	for( size_t i=0; i < p_anu_edge.size()-1; ++i )
		p_anu[i] = (p_anu_edge[i]+p_anu_edge[i+1])/2.;

	p_widflx.resize( p_anu_edge.size()-1 );
	for( size_t i=0; i < p_anu_edge.size()-1; ++i )
		p_widflx[i] = p_anu_edge[i+1] - p_anu_edge[i];

	p_anu2.resize( p_anu_edge.size()-1 );
	for( size_t i=0; i < p_anu_edge.size()-1; ++i )
		p_anu2[i] = pow2(p_anu[i]);

	p_anu3.resize( p_anu_edge.size()-1 );
	for( size_t i=0; i < p_anu_edge.size()-1; ++i )
		p_anu3[i] = pow3(p_anu[i]);

	p_anusqrt.resize( p_anu_edge.size()-1 );
	for( size_t i=0; i < p_anu_edge.size()-1; ++i )
		p_anusqrt[i] = sqrt(p_anu[i]);

	p_anuln.resize( p_anu_edge.size()-1 );
	for( size_t i=0; i < p_anu_edge.size()-1; ++i )
		p_anuln[i] = log(p_anu[i]);

	p_anulog10.resize( p_anu_edge.size()-1 );
	for( size_t i=0; i < p_anu_edge.size()-1; ++i )
		p_anulog10[i] = p_anuln[i]/LN_TEN;

	return;
}

// Set up special edges that need to be merged into the frequency mesh
void t_mesh::p_SetupEdges()
{
	DEBUG_ENTRY( "SetupEdges()" );

	// a list of major ionization edges that need to be fiddled into the frequency mesh
	// this list will be sorted at the end, so can be entered in any order

	// =================================================================================
	//
	// NB NB - make sure that all values in p_edges[] are slightly lower than the values
	//         used in the code !!! See the comments in ValidateEdges() for details.
	// NB NB - this is especially the case for the H I 1s and He II 2s & 2p edges !!!!
	// NB NB - the values below are (by design) not the most accurate values possible !!!
	// NB NB - to update, use the output from the failed test in ValidateEdges() and
	//         ALWAYS ROUND DOWN to 6 significant digits. It is necessary to round down
	//         to avoid epsilon sensitivity when using ipoint() for the same edge. The 6
	//         digits are needed to avoid problems with FLT_IS_DBL. If the value used by
	//         the code has less than 6 significant digits, you should still make it
	//         slightly smaller (like for the O III edges).
	//
	// =================================================================================

	// H I 1s ionization potential matching the one used in Cloudy, rel. accuracy ~ 6e-6
	p_edges.push_back( Energy(0.999466, "Ryd") );
	// H I 2s ionization potential matching the one used in Cloudy
	p_edges.push_back( Energy(0.249865, "Ryd") );
	// He I 1s2 ionization potential matching the one used in Cloudy
	p_edges.push_back( Energy(1.807139, "Ryd") );
	// He II 1s ionization potential matching the one used in Cloudy, rel. accuracy ~ 8e-5
	p_edges.push_back( Energy(3.9996296, "Ryd") );
	// He II 2s & 2p ionization potential, matching the one used in Cloudy
	// this edge needs to be here to assure that the He II 2s & 2p edges are in a higher
	// cell than the H I 1s edge. The code in cont_createpointers.cpp relies on this!!
	p_edges.push_back( Energy(0.999907, "Ryd") );
	// ionization potential from O III 2p2 1D
	p_edges.push_back( Energy(3.85499, "Ryd") );
	// ionization potential from O III 2p2 1S
	p_edges.push_back( Energy(3.64599, "Ryd") );

	sort( p_edges.begin(), p_edges.end() );
}

// Validate that the p_edges contains values in fair agreement with the rest of the code
void t_mesh::ValidateEdges() const
{
	DEBUG_ENTRY( "ValidateEdges()" );

	vector<Energy> cl_edges;

	// Make sure that the entries here match the ones defined above in SetupEdges()!!
	// These are the edge energies actually used by Cloudy
	cl_edges.push_back( Energy( iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].xIsoLevNIonRyd ) );
	cl_edges.push_back( Energy( iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH2s].xIsoLevNIonRyd ) );
	if( dense.lgElmtOn[ipHELIUM] )
	{
		cl_edges.push_back( Energy( iso_sp[ipHE_LIKE][ipHELIUM].fb[ipHe1s1S].xIsoLevNIonRyd ) );
		cl_edges.push_back( Energy( iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].xIsoLevNIonRyd ) );
		cl_edges.push_back( Energy( iso_sp[ipH_LIKE][ipHELIUM].fb[ipH2s].xIsoLevNIonRyd ) );
	}
	else
	{
		// helium is not enabled, so simply copy numbers from p_SetupEdges()
		cl_edges.push_back( Energy(1.807139*(1.+5.*double(FLT_EPSILON)), "Ryd") );
		cl_edges.push_back( Energy(3.9996296*(1.+5.*double(FLT_EPSILON)), "Ryd") );
		cl_edges.push_back( Energy(0.999907*(1.+5.*double(FLT_EPSILON)), "Ryd") );
	}
	cl_edges.push_back( Energy( opac.o3exc ) );
	cl_edges.push_back( Energy( opac.o3exc3 ) );

	sort( cl_edges.begin(), cl_edges.end() );

	if( cl_edges.size() != p_edges.size() )
	{
		fprintf( ioQQQ, "DISASTER INTERNAL ERROR: the cl_edges and p_edges arrays have different size.\n" );
		fprintf( ioQQQ, "In mesh.cpp, make sure the entries in these arrays match one on one.\n\n" );
		cdEXIT(EXIT_FAILURE);
	}
	bool lgErr = false;
	for( unsigned long i=0; i < p_edges.size(); ++i )
	{
		double rel_diff = 1. - p_edges[i].Ryd()/cl_edges[i].Ryd();
		// the entry in p_edge[i] should be slightly smaller than cl_edge[i]
		// this avoids problems when calculating ipoint(cl_edge[i])...
		// it should be at least a factor 1.5*FLT_EPSILON smaller to avoid problems
		// when using -DFLT_IS_DBL, but is should also not be too large...
		if( rel_diff < 1.5*double(FLT_EPSILON) || rel_diff > 1.e-5 )
		{
			double expected = cl_edges[i].Ryd()*(1. - 5.*double(FLT_EPSILON));
			fprintf( ioQQQ, "DISASTER INTERNAL ERROR: energy mismatch for p_edges[].\n" );
			fprintf( ioQQQ, "I expected %.8g but found %.8g.\n", expected, p_edges[i].Ryd() );
			fprintf( ioQQQ, "Please update p_edges[] to hold %.8g in mesh.cpp.\n\n", expected );
			lgErr = true;
		}
	}

	// make sure the HeII Balmer pointers are above the HI Lyman pointer
	//
	// H and He are special because of their large abundances - the OTS fields they produce are
	// very important in establishing the ionization in the gas. The OTS fields are placed in a
	// single cell at an element's threshold. When that shell's photoionization rate is
	// evaluated the first cell, with the OTS field, is not included - that would be double
	// counting. The self-OTS is taken into account with changes to the recombination
	// coefficient.
	//
	// Note that this code makes implicit assumptions about the abundance pattern, which is a bug
	//
	// The values in p_edges[] should be chosen such that the test below should always pass. If
	// these tests trip, update the energies of the H end He edges in p_SetupEdges(). Normally
	// the tests above should already have caught that.
	if( dense.lgElmtOn[ipHELIUM] )
	{
		if( iso_sp[ipH_LIKE][ipHELIUM].fb[ipH2s].ipIsoLevNIonCon <=
			iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon )
		{
			fprintf( ioQQQ, "DISASTER INTERNAL ERROR: incorrect ionization edge pointers.\n" );
			fprintf( ioQQQ, "The He 2s ionization energy should be in a higher\n" );
			fprintf( ioQQQ, "frequency cell than the H 1s ionization energy.\n" );
			fprintf( ioQQQ, "He H edges were %li %li.\n\n", iso_sp[ipH_LIKE][ipHELIUM].fb[ipH2s].ipIsoLevNIonCon ,
					iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon );
			lgErr = true;
		}
		if( iso_sp[ipH_LIKE][ipHELIUM].fb[ipH2p].ipIsoLevNIonCon <=
			iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon )
		{
			fprintf( ioQQQ, "DISASTER INTERNAL ERROR: incorrect ionization edge pointers.\n" );
			fprintf( ioQQQ, "The He 2p ionization energy should be in a higher\n" );
			fprintf( ioQQQ, "frequency cell than the H 1s ionization energy.\n" );
			fprintf( ioQQQ, "He H edges were %li %li.\n\n", iso_sp[ipH_LIKE][ipHELIUM].fb[ipH2p].ipIsoLevNIonCon ,
					iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon );
			lgErr = true;
		}
	}

	if( lgErr )
	{
		cdEXIT(EXIT_FAILURE);
	}
}

/*CheckMesh perform sanity check confirming that the energy array has been properly filled */
void t_mesh::CheckMesh() const
{
	DEBUG_ENTRY( "CheckMesh()" );

	// is the grid initialized?
	ASSERT( p_anu.size() > 0 );

	// check outer bounds of the complete grid
	ASSERT( fp_equal( anumin(0), emm() ) );
	ASSERT( fp_equal( anumax(ncells()-1), egamry() ) );

	bool lgFail = false;
	double hi_bound = emm();
	for( size_t n=0; n < p_RangeUpperLimit.size(); ++n )
	{
		/* test middle of energy range */
		double lo_bound = hi_bound;
		hi_bound = p_RangeUpperLimit[n];

		double energy = lo_bound*0.5 + hi_bound*0.5;
		size_t ic = ipointC(energy);
		if( energy < anumin(ic) )
		{
			fprintf( ioQQQ, "CheckMesh middle test low fail\n" );
			lgFail = true;
		}
		else if( energy > anumax(ic) )
		{
			fprintf( ioQQQ, "CheckMesh middle test high fail\n" );
			lgFail = true;
		}

		/* test near low bound */
		energy = lo_bound*0.99 + hi_bound*0.01;
		ic = ipointC(energy);
		if( energy < anumin(ic) )
		{
			fprintf( ioQQQ, "CheckMesh low test low fail\n" );
			lgFail = true;
		}
		else if( energy > anumax(ic) )
		{
			fprintf( ioQQQ, "CheckMesh low test high fail\n" );
			lgFail = true;
		}

		/* test near high bound */
		energy = lo_bound*0.01 + hi_bound*0.99;
		ic = ipointC(energy);
		if( energy < anumin(ic) )
		{
			fprintf( ioQQQ, "CheckMesh high test low fail\n" );
			lgFail = true;
		}
		else if( energy > anumax(ic) )
		{
			fprintf( ioQQQ, "CheckMesh high test high fail\n" );
			lgFail = true;
		}
	}

	if( lgFail )
	{
		fprintf( ioQQQ, "CheckMesh: sanity check on frequency mesh failed.\n Bailing out\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return;
}
