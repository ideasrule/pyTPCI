/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ATOMS_H_
#define ATOMS_H_

#include "container_classes.h"
#include "module.h"

class TransitionProxy;

/**atom_level2 do level population and cooling for two level atom 
 \param t
*/
void atom_level2( const TransitionProxy &t );
void atom_level2( const TransitionProxy &t, const bool lgHFS );

/**
atom_levelN - compute populations of arbitrary n-level atom
\param nlev nlev is the number of levels to compute
\param abund ABUND is total abundance of species, used for nth equation
\param g[] G(ndim) is stat weight of levels
\param ex[] EX(ndim) is excitation potential of levels, either wn or deg K 0 for first one, NOT d(ENER), but energy rel to ground 
\param chExUnits this is 'K' for above ex[] as Kelvin deg, is 'w' for wavenumbers
\param pops[] populations of each level as deduced here
\param depart[] departure coefficient derived here
\param AulEscp net transition rate, A * esc prob, s-1
\param AulDest AulDest(ilo,ihi) is destruction rate, from up to low, A * dest prob, [s-1], asserts confirm that ihi,lo is zero
\param AulPump AulPump(lo, hi) is pumping rate, A * occ num, (hi,lo) must be zero, [s-1]  
\param CollRate collision rates, evaluated here and returned for cooling by calling function, unless following flag is true.  
	If true then calling function has already filled in these rates.  CollRate[i][j] is rate from i to j 
\param create this is an additional creation rate, normally zero, units cm-3 s-1
\param destroy[] this is an additional destruction rate to continuum, normally zero, units s-1 
\param lgCollRateDone flag saying whether CollRate already done, or we need to do it here
\param cooltl total cooling, set here but nothing done with it
\param coolder derivative of cooling, set here but nothing done with it
\param chLabel string used to identify calling program in case of error
\param lgNegPop lgNegPop flag indicating what we have done positive if negative populations occurred zero if normal calculation done negative if too cold (for some atoms other routine will be called in this case)
\param lgZeroPop true if populations are zero, either due to zero abundance of very low temperature
\param lgDeBug option to print matrices for debugging
\post atoms.PopLevels[n], atoms.DepLTELevels[n] are set lines added to ots array
\param grnd_excit is the total excitation rate out of ground state
*/ 

class setCollRate
{
	long nLevelAlloc;
	vector<double*> excit;
	vector<double> excit_b;
public:
   setCollRate() : nLevelAlloc(0) {}
   void resize(long int nlev)
	{
		if (nlev > nLevelAlloc)
		{
			excit.resize(nlev);
			excit_b.resize(0.5*nlev*(nlev-1));
			size_t ibase = 0;
			for (long ihi=1; ihi<nlev; ++ihi)
			{
				excit[ihi] = &excit_b[ibase];
				ibase += ihi;
			}
			ASSERT(ibase == excit_b.size());
			nLevelAlloc = nlev;
		}
	}
	void operator()(
		long nlev,
		double TeInverse,
		double **col_str,
		const double ex[],	
		const double g[], 
		double **CollRate
		);
};

class Atom_LevelN
{
	// these are all automatically deallocated when they go out of scope
	valarray<double> bvec;
	multi_arr<double,2,C_TYPE> amat;
public:
	Atom_LevelN() {}
	void operator()(
		long nLevelCalled,
		realnum abund, 
		const vector<double>& g, 
		const vector<double>& ex, 
		char chExUnits,
		vector<double>& pops, 
		vector<double>& depart,
		multi_arr<double,2>& AulEscp, 
		multi_arr<double,2>& AulDest, 
		multi_arr<double,2>& AulPump, 
		const multi_arr<double,2>& CollRate,
		const vector<double>& create,
		const vector<double>& destroy,
		double *cooltl, 
		double *coolder, 
		const char *chLabel, 
		bool lgPrtMatrix,
		bool lgImgMatrix,
		int *nNegPop,
		bool *lgZeroPop,
		bool lgDeBug,
		bool lgLTE=false,
		multi_arr<double,2> *Cool=NULL,
		multi_arr<double,2> *dCooldT=NULL,
		double *grnd_excit = NULL);
};


/** number of levels in PopLevels, DepLTELevels, vectors*/
const long LIMLEVELN = 20L;

struct t_atoms : public module {

	const char *chName() const
	{
		return "atoms";
	}
	void zero();
	void comment(t_warnings&) {}

	/** photo excitation of excited state of NI */
	realnum p2nit, d5200r;

	/** excited state of Mg+ */
	realnum xMg2Max,
		/** its population */
		popMg2,
		/** photoionization rate */
		rateMg2;

	/**
	 * this stores most recently evaluated level populations
	 * resulting from the leven family of routines
	 * PopLevels is population (cm^-3) of the levels
	 * DepLevels is lte departure coef
	 */
	double PopLevels[LIMLEVELN+1], 
	  DepLTELevels[LIMLEVELN+1];

	};
extern t_atoms atoms;

#endif /* ATOMS_H_ */
