/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"

#include "dense.h"

#include "abund.h"
#include "colden.h"
#include "conv.h"
#include "dynamics.h"
#include "elementnames.h"
#include "deuterium.h"
#include "hmi.h"
#include "phycon.h"
#include "radius.h"
#include "struc.h"
#include "thermal.h"
#include "trace.h"
#include "iso.h"
#include "h2.h"
#include "mole.h"
#include "save.h"

void t_dense::updateXMolecules()
{
	total_molecule_elems(m_xMolecules);
}

void t_dense::zero()
{
	for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		dense.SetGasPhaseDensity( nelem, 0. );
		for( long ion=0; ion < LIMELM+1; ion++ )
		{
			dense.xIonDense[nelem][ion] = 0.;
		}
	}
	for (long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem)
		m_xMolecules[nelem] = 0.f;
	// needed for TempChange to work but arrays needed for EdenChange to
	// work are not yet defined
	dense.eden = 1.;

	dense.xMassTotal = 0.;	
	dense.xNucleiTotal = 1.;
	/* WJH */
	dense.xMassDensity0 = -1.0f; 

	/* this is a scale factor that changes the n(H0)*1.7e-4 that is added to the
	* electron density to account for collisions with atomic H.  it is an order of
	* magnitude guess, so this command provides ability to see whether it affects results */
	dense.HCorrFac = 1.f;
}

void ScaleAllDensities(realnum factor)
{
	double edensave = dense.eden;

	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		if (dense.lgElmtOn[nelem])
		{
			ScaleIonDensities( nelem, factor );
			dense.SetGasPhaseDensity( nelem, dense.gas_phase[nelem] * factor);
		}
	}

	EdenChange( dense.eden * factor );
	
	if( trace.lgTrace && trace.lgNeBug )
	{
		fprintf( ioQQQ, " EDEN change PressureChange from to %10.3e %10.3e %10.3e\n", 
					edensave, dense.eden, edensave/dense.eden );
	}

	hmi.H2_total *= factor;
	h2.ortho_density *= factor;
	h2.para_density *= factor;
	for( long mol=0; mol < mole_global.num_total; mol++ )
	{
		mole.species[mol].den *= factor;
	}
	dense.updateXMolecules();
	deut.updateXMolecules();

	ASSERT(lgElemsConserved());
}

void ScaleIonDensities( const long nelem, const realnum factor )
{
	double renorm;
	for( long ion=0; ion < (nelem + 2); ion++ )
	{
		dense.xIonDense[nelem][ion] *= factor;
		if (nelem-ion >= 0 && nelem-ion < NISO)
			iso_renorm(nelem, nelem-ion, renorm);
	}
	
	if( nelem==ipHYDROGEN && deut.lgElmtOn )
		ScaleDensitiesDeuterium( factor );

	return;
}

void t_dense::SetGasPhaseDensity( const long nelem, const realnum density )
{
	gas_phase[nelem] = density;

	if( nelem==ipHYDROGEN && deut.lgElmtOn )
	{
		// fractionation is taken care of inside called routine
		SetGasPhaseDeuterium( density );
	}

	return;
}

bool lgElemsConserved (void)
{
	bool lgOK=true;

	/* this checks all molecules and ions add up to abundance */
	for( ChemNuclideList::iterator nuc = nuclide_list.begin(); nuc != nuclide_list.end(); ++nuc)
	{
		long nelem = (*nuc)->el()->Z-1;
		
		double sum_monatomic = 0.;
		realnum sum_molecules = 0.f;
		realnum sum_gas_phase = 0.f;

		if( (*nuc)->label() == "D" )
		{
			sum_monatomic = deut.xIonDense[0] + deut.xIonDense[1];
			sum_molecules = deut.xMolecules();
			sum_gas_phase = deut.gas_phase;

		}
		/* check that species add up if element is turned on */
		else if( dense.lgElmtOn[nelem] )
		{
			/* this sum of over the molecules did not include the atom and first
			 * ion that is calculated in the molecular solver */
			for( long ion=0; ion<nelem+2; ++ion )
			{
				sum_monatomic += dense.xIonDense[nelem][ion] * (*nuc)->frac;
			}
			sum_molecules = dense.xMolecules(nelem) * (*nuc)->frac;
			sum_gas_phase = dense.gas_phase[nelem] * (*nuc)->frac;
		}
		else
			continue;


		if( sum_monatomic + sum_molecules <= SMALLFLOAT && 
			sum_gas_phase > SMALLFLOAT)
		{
			fprintf(ioQQQ,"PROBLEM non-conservation of nuclei %s\tions %g moles %g error %g of %g\n",
				(*nuc)->label().c_str(),
				sum_monatomic,
				sum_molecules,
				sum_monatomic + sum_molecules-sum_gas_phase,
				sum_gas_phase );
			lgOK=false;
		}

		/* check that sum agrees with total gas phase abundance */
		/** \todo	2	this assert is not passed if error made much smaller.  This
		 * error should be related to a check on convergence of the molecular
		 * networks and their sum rules, with a criteria used here and there */
		if( fabs( sum_monatomic + sum_molecules - sum_gas_phase ) > 
			conv.GasPhaseAbundErrorAllowed * sum_gas_phase )
		{
			fprintf(ioQQQ,"PROBLEM non-conservation of nuclei %s\t nzone %li atoms %.12e moles %.12e sum %.12e tot gas %.12e rel err %.3e\n",
				(*nuc)->label().c_str(),
				nzone,
				sum_monatomic,
				sum_molecules,
				sum_monatomic + sum_molecules,
				sum_gas_phase,
				(sum_monatomic + sum_molecules - sum_gas_phase)/sum_gas_phase );
			fprintf(ioQQQ," Consider increasing value of conv.GasPhaseAbundErrorAllowed with SET DENSITY TOLERANCE as a workaround. "
					"Its current value is %.2e\n" , conv.GasPhaseAbundErrorAllowed );
			fprintf(ioQQQ," Temperature %.2e H density %.2e zone %ld \n", phycon.te, dense.gas_phase[ipHYDROGEN], nzone);
			lgOK=false;
		}
		
		if( 0 )
		{
			// print reactions involving the neutral 
			mole_print_species_reactions( findspecies( elementnames.chElementSym[nelem] ) );
		}
	}

	return lgOK;
}

void lgStatesConserved ( long nelem, long ionStage, qList states, long numStates, realnum err_tol, long loop_ion )
{
	if( 0 && conv.lgConvIoniz() && 
		dense.lgElmtOn[nelem] && 
		dense.xIonDense[nelem][ionStage]/dense.eden > conv.EdenErrorAllowed/10. &&
		dense.xIonDense[nelem][ionStage] > SMALLFLOAT )
	{
		double abund = 0.;
		for (long ipLev=0; ipLev<numStates; ipLev++)
		{
			abund += states[ipLev].Pop();
		}

		double rel_err = fabs(abund-dense.xIonDense[nelem][ionStage])/
			SDIV(dense.xIonDense[nelem][ionStage]);

		//fprintf(ioQQQ,"Iso consistency error %ld %ld %g\n",nelem,ionStage,rel_err);

		if( rel_err > err_tol ) 
		{
			/* we'll set the flag for non-convergence, but don't bother complaining in printout for first sweeps */
			if( 0 && nzone > 0 && loop_ion > 0 )
			{	
				fprintf(ioQQQ,"PROBLEM Inconsistent states/stage pops nzone %3ld loop_ion %2ld nelem %2ld ion %2ld states = %e stage = %e error = %e\n",
					nzone,
					loop_ion,
					nelem,
					ionStage,
					abund,
					dense.xIonDense[nelem][ionStage],
					rel_err );
			}
			ostringstream chConvIoniz;
			chConvIoniz << "States!=stage pops nelem " << nelem << " ion " << ionStage;
			conv.setConvIonizFail( chConvIoniz.str().c_str(), abund,
					       dense.xIonDense[nelem][ionStage] );
		}
	}
}

void SumDensities( void )
{
	realnum den_mole,
		DenseAtomsIons;

	/* find total number of atoms and ions
	 * at this point does not include molecules */
	DenseAtomsIons = 0.;
	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem]  )
		{
			for( long ion=0; ion<=nelem+1; ++ion )
				DenseAtomsIons += dense.xIonDense[nelem][ion];
		}
	}

	/* find number of molecules of the heavy elements in the gas phase.
	 * Atoms and ions were counted above.  Do not include ices, solids
	 * on grains */
	den_mole = total_molecules_gasphase();
	
	/* total number of atoms, ions, and molecules in gas phase per unit vol */
	dense.xNucleiTotal = den_mole + DenseAtomsIons;
	if( dense.xNucleiTotal > BIGFLOAT )
	{
		fprintf(ioQQQ, "PROBLEM DISASTER SumDensities has found "
			"dense.xNucleiTotal with an insane density.\n");
		fprintf(ioQQQ, "The density was %.2e\n",
			dense.xNucleiTotal);
		TotalInsanity();
	}
	ASSERT( dense.xNucleiTotal > 0. );

	/* particle density that enters into the pressure includes electrons
	 * cm-3 */
	dense.pden = (realnum)(dense.eden + dense.xNucleiTotal);
        
	/* dense.wmole is molecular weight, AtomicMassUnit per particle */
	dense.wmole = 0.;
	for( long i=0; i < LIMELM; i++ )
	{
		dense.wmole += dense.gas_phase[i]*(realnum)dense.AtomicWeight[i];
	}
        /* dense.wmole is now molecular weight, AtomicMassUnit per particle */
        dense.wmole /= dense.pden;
        ASSERT( dense.wmole > 0. && dense.pden > 0. );
        
        /* xMassDensity is density in grams per cc */
        /** \todo       2       - should this include mass in grains? */
        /** \todo       2       - should this include mass in grain mantle ice deposits? */
        dense.xMassDensity = (realnum)(dense.wmole*ATOMIC_MASS_UNIT*dense.pden);
        
        /*>>chng 04 may 25, introduce following test on xMassDensity0 
         * WJH 21 may 04, this is the mass density that corresponds to the hden 
         * specified in the init file. It is used to calculate the mass flux in the 
         * dynamics models. It may not necessarily be the same as struc.DenMass[0],
         * since the pressure solver may have jumped to a different density at the
         * illuminated face from that specified.*/   
        if( dense.xMassDensity0 < 0.0 )
                 dense.xMassDensity0 = dense.xMassDensity;

	return;
}

bool AbundChange( )
{
	DEBUG_ENTRY( "AbundChange()" );

	fixit("AbundChange breaks conservation if molecular abundance is finite");

	// Suggested improved procedure:
	// 1. Scale all molecules by the scaling for their most restrictive constituent
	// 2. dense.updateXMolecules();
	// 3. Scale ionic species so as to maintain conservation
	// This will probably lead to material being dumped into the ionic species,
	// but this should slosh back at the next molecular network solve.  This
	// procedure is at least a) conservative; b) consistent with the pure-ion case.

	/* this will be set true if density or abundances change in this zone */
	bool lgChange = false;

	/* >>chng 97 jun 03, added variable abundances for element table command
	 * option to get abundances off a table with element read command */
	if( abund.lgAbTaON )
	{
		lgChange = true;
		for( long nelem=1; nelem < LIMELM; nelem++ )
		{
			if( abund.lgAbunTabl[nelem] )
			{
				double abun = AbundancesTable(radius.Radius,radius.depth,nelem+1)*
				  dense.gas_phase[ipHYDROGEN];

				double hold = abun/dense.gas_phase[nelem];
				dense.SetGasPhaseDensity( nelem, (realnum)abun );

				for( long ion=0; ion < (nelem + 2); ion++ )
				{
					dense.xIonDense[nelem][ion] *= (realnum)hold;
				}
			}
		}
	}

	/* this is set false if fluctuations abundances command entered,
	 * when density variations are in place this is true,
	 * and dense.chDenseLaw is "SINE" */
	double FacAbun = 1.;
	if( !dense.lgDenFlucOn )
	{
		static double FacAbunSav;
		double OldAbun = 0.0;
		/* abundance variations are in place */
		lgChange = true;
		if( nzone > 1 )
		{
			OldAbun = FacAbunSav;
		}

		/* rapid abundances fluctuation */
		if( dense.lgDenFlucRadius )
		{
			/* cycle over radius */
			FacAbunSav = dense.cfirst*cos(radius.depth*dense.flong+
				dense.flcPhase) + dense.csecnd;
		}
		else
		{
			/* cycle over column density */
			FacAbunSav = dense.cfirst*cos( colden.colden[ipCOL_HTOT]*dense.flong+
				dense.flcPhase) + dense.csecnd;
		}

		if( nzone > 1 )
		{
			FacAbun = FacAbunSav/OldAbun;
		}
	}

	if( FacAbun != 1. )
	{
		ASSERT(!dynamics.lgAdvection); // Local abundance change doesn't make sense with advective transport

		/* Li on up is affect by abundance variations */
		for( long nelem=ipLITHIUM; nelem < LIMELM; nelem++ )
		{
			if (dense.lgElmtOn[nelem])
			{
				dense.SetGasPhaseDensity(nelem, dense.gas_phase[nelem] * (realnum) FacAbun);
				ScaleIonDensities( nelem, (realnum)(FacAbun) );
			}
		}

		for( long mol=0; mol < mole_global.num_total; mol++ )
		{
			mole.species[mol].den *= (realnum)(FacAbun);
		}
	}

	if (lgChange)
	{
		/* must call TempChange since ionization has changed, there are some
		 * terms that affect collision rates (H0 term in electron collision) */
		TempChange(phycon.te , false);
	}

	return lgChange;
}

namespace {
	const int SCALE_NEW = 1;
}

realnum scalingDensity(void)
{
	if (SCALE_NEW)
		return dense.xMassDensity/(realnum)ATOMIC_MASS_UNIT;
	else
		return dense.gas_phase[ipHYDROGEN];
}
realnum scalingZoneDensity(long i)
{
	if (SCALE_NEW)
		return struc.DenMass[i]/(realnum)ATOMIC_MASS_UNIT;
	else
		return struc.hden[i];
}
