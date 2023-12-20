/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MESH_H_
#define MESH_H_

#include "thirdparty.h"
#include "energy.h"

class t_mesh {

	/** ================================================================================= */
	/** the following define the continuum energy scale and its limits */

	/** the energy of the lower limit low-energy limit of the continuum */
	double p_emm;

	/** the energy of the upper limit high-energy limit of the continuum */
	double p_egamry;

	/** factor to reset continuum resolution from continuum_mesh.ini,
	 * default is unity, reset with set continuum resolution command */
	double p_ResolutionScaleFactor;

	/** MD5 sum of the continuum_mesh.ini file */
	string p_mesh_md5sum;

	/** this is information needed to set the energy binning,
	 * full continuum is described by series of ranges where resolution is
	 * constant over that range */
	vector<double> p_RangeUpperLimit;
	vector<double> p_RangeResolution;

	/** a list of major ionization edges that need to be fiddled into the frequency mesh */
	vector<Energy> p_edges;

	/** energy in Ryd of center of cell */
	vector<double> p_anu;

	/** energy in Ryd at the edges of each cell */
	vector<double> p_anu_edge;

	/** width of cells in Rydberg */
	vector<double> p_widflx;

	/** these are log, sqrt, square, and cube of anu array */
	vector<double> p_anulog10;
	vector<double> p_anuln;
	vector<double> p_anusqrt;
	vector<double> p_anu2;
	vector<double> p_anu3;

	/* this routine reads continuum_mesh.ini */
	void p_ReadResolution();

	/* this routine defines the frequency mesh */
	void p_SetupMesh(bool lgUnitCell);

	// Set up special edges that need to be merged into the frequency mesh
	void p_SetupEdges();
public:
	/* set up the frequency mesh */
	void InitMesh(bool lgUnitCell)
	{
		if( lgMeshSetUp() )
			return;

		// these are the low and high energy bounds of the continuum
		Energy Elo( 10., "MHz" );
		p_emm = Elo.Ryd();
		Energy Ehi( 100., "MeV" );
		p_egamry = Ehi.Ryd();

		p_mesh_md5sum = MD5datafile( "continuum_mesh.ini" );

		p_SetupEdges();
		p_ReadResolution();
		p_SetupMesh(lgUnitCell);
	}
	/* perform sanity checks on the frequency mesh */
	void ValidateEdges() const;
	void CheckMesh() const;

	bool lgMeshSetUp() const
	{
		return ( p_anu.size() > 0 );
	}
	// getters and setters
	long ncells() const
	{
		return long(p_anu.size());
	}
	double emm() const
	{
		return p_emm;
	}
	double egamry() const
	{
		return p_egamry;
	}
	double getResolutionScaleFactor() const
	{
		return p_ResolutionScaleFactor;
	}
	void setResolutionScaleFactor(double fac)
	{
		if( !lgMeshSetUp() )
			p_ResolutionScaleFactor = fac;
		else
			ASSERT( fp_equal(fac,p_ResolutionScaleFactor) );
	}
	string mesh_md5sum() const
	{
		return p_mesh_md5sum;
	}
	const double* anuptr() const
	{
		return get_ptr(p_anu);
	}
	double anu(size_t i) const
	{
		return p_anu[i];
	}
	double anu2(size_t i) const
	{
		return p_anu2[i];
	}
	double anu3(size_t i) const
	{
		return p_anu3[i];
	}
	double anuln(size_t i) const
	{
		return p_anuln[i];
	}
	const double* anulog10ptr() const
	{
		return get_ptr(p_anulog10);
	}	
	double anulog10(size_t i) const
	{
		return p_anulog10[i];
	}
	double anusqrt(size_t i) const
	{
		return p_anusqrt[i];
	}
	double anumin(size_t i) const
	{
		return p_anu_edge[i];
	}
	double anumax(size_t i) const
	{
		return p_anu_edge[i+1];
	}
	double widflx(size_t i) const
	{
		return p_widflx[i];
	}
	// return index ic such that p_anu_edge[ic] <= anu < p_anu_edge[ic+1]
	size_t ipointC(double anu) const
	{
		ASSERT( lgMeshSetUp() );
		ASSERT( anu >= p_anu_edge.front() && anu <= p_anu_edge.back() );
		return hunt_bisect( p_anu_edge, anu );
	}
	// same as ipointC, but on fortran scale
	size_t ipointF(double anu) const
	{
		return ipointC(anu) + 1;
	}
	// return smallest index ic such that p_anu[ic] >= threshold
	// this guarantees that rfield.anu(ic) - threshold is always non-negative.
	size_t ithreshC(double threshold) const
	{
		ASSERT( lgMeshSetUp() );
		ASSERT( threshold >= p_anu.front() && threshold <= p_anu.back() );
		size_t ic = hunt_bisect( p_anu, threshold );
		// this if is nearly always taken...
		if( p_anu[ic] < threshold )
			++ic;
		return ic;
	}
	// same as ithreshC, but on fortran scale
	size_t ithreshF(double threshold) const
	{
		return ithreshC(threshold) + 1;
	}
	bool isEnergyBound( Energy en ) const
	{
		return en.Ryd() > emm() && en.Ryd() < egamry();
	}

	// constructor
	t_mesh()
	{
		// these will be initialized in InitMesh()
		p_emm = 0.;
		p_egamry = 0.;
		// this is set with the set continuum resolution command
		p_ResolutionScaleFactor = 1.;
	}
};

#endif /* MESH_H_ */
