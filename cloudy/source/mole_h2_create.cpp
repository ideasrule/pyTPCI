/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*H2_Create create variables for the H2 molecule, called by ContCreatePointers after continuum
 * mesh has been set up */
#include "cddefines.h" 
#include "atmdat.h"
#include "taulines.h"
#include "opacity.h" 
#include "hmi.h" 
#include "ipoint.h"
#include "grainvar.h"
#include "h2.h"
#include "mole.h"
#include "lines_service.h"
#include "iso.h"
#include "ran.h"

/* if this is set true then code will print energies and stop */
enum {DEBUG_ENER=false};

/* this is equation 8 of Takahashi 2001, clearer definition is given in
 * equation 5 and following discussion of
 * >>refer	H2	formation	Takahashi, J., & Uehara, H., 2001, ApJ, 561, 843-857
 * 0.27eV, convert into wavenumbers */
static double XVIB[H2_TOP] = { 0.70 , 0.60 , 0.20 };
static double Xdust[H2_TOP] = { 0.04 , 0.10 , 0.40 };

/* this is energy difference between bottom of potential well and 0,0
 * the Takahashi energy scale is from the bottom,
 * 2201.9 wavenumbers  */
static const double energy_off = 0.273*FREQ_1EV/SPEEDLIGHT;

STATIC double EH2_eval(  int ipH2, double DissocEnergy, double energy_wn )
{
	double EH2_here;
	double Evm = DissocEnergy * XVIB[ipH2] + energy_off;

	double Ev = (energy_wn+energy_off);
	/* equation 9 of Takahashi 2001 which is only an approximation
	 * equation 1, 2 of 
	 * Takahashi, Junko, & Uehara, Hideya, 2001, ApJ, 561, 843-857,
	 * this is heat deposited on grain by H2 formation in this state */
	double Edust = DissocEnergy * Xdust[ipH2] *
		( 1. - ( (Ev - Evm) / (DissocEnergy+energy_off-Evm)) *
		( (1.-Xdust[ipH2])/2.) );
	ASSERT( Edust >= 0. );

	/* energy is total binding energy less energy lost on grain surface 
	 * and energy offset */
	EH2_here = DissocEnergy +energy_off - Edust;
	ASSERT( EH2_here >= 0.);

	return EH2_here;
}

/*H2_vib_dist evaluates the vibration distribution for H2 formed on grains */
STATIC double H2_vib_dist( int ipH2 , double EH2, double DissocEnergy, double energy_wn)
{
	double G1[H2_TOP] = { 0.3 , 0.4 , 0.9 };
	double G2[H2_TOP] = { 0.6 , 0.6 , 0.4 };
	double Evm = DissocEnergy * XVIB[ipH2] + energy_off;
	double Fv;
	if( (energy_wn+energy_off) <= Evm )
	{
		/* equation 4 of Takahashi 2001 */
		Fv = sexp( POW2( (energy_wn+energy_off - Evm)/(G1[ipH2]* Evm ) ) );
	}
	else
	{
		/* equation 5 of Takahashi 2001 */
		Fv = sexp( POW2( (energy_wn+energy_off - Evm)/(G2[ipH2]*(EH2 - Evm ) ) ) );
	}
	return Fv;
}

STATIC bool lgRadiative( const TransitionList::iterator &tr )
{
	//quantumState_diatoms *Hi = static_cast<quantumState_diatoms*>( tr.Hi );	
	qList::iterator Lo = (*tr).Lo() ;	
	//if( Lo->n==0 && lgH2_radiative[ ipEnergySort[(*Hi).n][(*Hi).v][(*Hi).J] ][ ipEnergySort[Lo->n][Lo->v][Lo->J] ] )
	if( (*Lo).n()==0 && (*tr).Emis().Aul() > 1.01e-30 )
		return true;
	else
		return false;
}

STATIC bool compareEmis( const TransitionList::iterator& tr1, const TransitionList::iterator &tr2 )
{
	if( lgRadiative( tr1 ) && !lgRadiative( tr2 ) )
		return true;
	else 
		return false; 
}

#if 0
STATIC bool compareEnergies( qStateProxy st1, qStateProxy st2 )
{
	if( st1.energy().Ryd() <= st2.energy().Ryd() )
		return true;
	else 
		return false; 
}
#endif

/* create variables for the H2 molecule, called by
 * ContCreatePointers after continuum mesh has been set up */
void diatomics::init(void)
{
	/* this is flag set above - when true h2 code is not executed - this is way to
	 * avoid this code when it is not working */
	/* only allocate vectors one time per core load */
	if( lgREAD_DATA || !lgEnabled )
		return;

	DEBUG_ENTRY( "H2_Create()" );

	bool lgDebugPrt = false;

	/* print string if H2 debugging is enabled */
	if( lgDebugPrt )
		fprintf(ioQQQ," H2_Create called in DEBUG mode.\n");
		
	// find the species in the chemistry network	
	sp = findspecies( label.c_str() );
	ASSERT( sp != null_mole );

	// find the same species with the excitation marker	
	string strSpStar = shortlabel + "*";
	sp_star = findspecies( strSpStar.c_str() );
	ASSERT( sp != null_mole );
		
	mass_amu = sp->mole_mass / ATOMIC_MASS_UNIT;
	ASSERT( mass_amu > 0. );

	H2_ReadDissocEnergies();

	/* This does more than just read energies... it actually defines the states */
	H2_ReadEnergies();

	// now sort states into energy order
	fixit("this doesn't work!");
	//sort( states.begin(), states.end(), compareEnergies );
	
	ASSERT( n_elec_states > 0 );
	/* create a vector of sorted energies */
	ipEnergySort.reserve(n_elec_states);
	for( long iElecHi=0; iElecHi<n_elec_states; ++iElecHi )
	{
		ipEnergySort.reserve(iElecHi,nVib_hi[iElecHi]+1);
		for( long iVibHi = 0; iVibHi <= nVib_hi[iElecHi]; ++iVibHi )
		{
			ipEnergySort.reserve(iElecHi,iVibHi,nRot_hi[iElecHi][iVibHi]+1);
		}
	}
	ipEnergySort.alloc();

	/* create arrays for energy sorted referencing of e, v, J */
	ipElec_H2_energy_sort.resize( states.size() );
	ipVib_H2_energy_sort.resize( states.size() );
	ipRot_H2_energy_sort.resize( states.size() );
	for( unsigned nEner = 0; nEner < states.size(); ++nEner )
	{
		long iElec = states[nEner].n();
		long iVib = states[nEner].v();
		long iRot = states[nEner].J();
		ipElec_H2_energy_sort[nEner] = iElec;
		ipVib_H2_energy_sort[nEner] = iVib;
		ipRot_H2_energy_sort[nEner] = iRot;
		ipEnergySort[iElec][iVib][iRot] = nEner;
		/* following will print quantum indices and energies */
		/*fprintf(ioQQQ,"%li\t%li\t%li\t%.3e\n", iElec, iVib, iRot, 
			states[nEner].energy().WN() );*/
	}
	
	H2_lgOrtho.alloc( ipEnergySort.clone() );

	/* set all statistical weights - ours is total statistical weight - 
	 * including nuclear spin */
	for( qList::iterator st = states.begin(); st != states.end(); ++st )
	{
		if( this==&h2 )
		{
			/* unlike atoms, for H2 nuclear spin is taken into account - so the
			 * statistical weight of even and odd J states differ by factor of 3 - see page 166, sec par
			 * >>>refer	H2	H2_stat wght	Shull, J.M., & Beckwith, S., 1982, ARAA, 20, 163-188 */
			/* This integer is added to rotation quantum number J for the test of whether
			 * a particular J state is ortho or para - the state is ortho if J+below is odd,
			 * and para if J+below is even */
			const int H2_nRot_add_ortho_para[N_ELEC] = {0, 1, 1, 0, 1, 1, 0};
			if( is_odd( (*st).J() + H2_nRot_add_ortho_para[(*st).n()]) )
			{
				/* ortho */
				H2_lgOrtho[(*st).n()][(*st).v()][(*st).J()] = true;
			}
			else
			{
				/* para */
				H2_lgOrtho[(*st).n()][(*st).v()][(*st).J()] = false;
			}
		}
		else if( this==&hd )
		{
			// No ortho-para distinction, set all of these to false.
			H2_lgOrtho[(*st).n()][(*st).v()][(*st).J()] = false;
		}
		else
			// This will have to change for any additional molecules.
			TotalInsanity();

		realnum rotstat = 2.f*(*st).J()+1.f;
		realnum opstat = 1.f;
		if (H2_lgOrtho[(*st).n()][(*st).v()][(*st).J()])
			opstat = 3.f;

		(*st).g() = opstat*rotstat;
	}

	if( lgDebugPrt )
	{
		for( long iElec=0; iElec<n_elec_states; ++iElec)
		{
			/* print the number of levels within iElec */
			fprintf(ioQQQ,"\t(%li %li)", iElec ,  nLevels_per_elec[iElec] );
		}
		fprintf(ioQQQ,
			" H2_Create: there are %li electronic levels, in each level there are",
			n_elec_states);
		fprintf(ioQQQ,
				  " for a total of %li levels.\n", (long int) states.size() );
	}

	/* now find number of levels in H2g */
	for( long nEner=0; nEner<nLevels_per_elec[0]; ++nEner )
	{
		if( states[nEner].energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
			break;
		nEner_H2_ground = nEner;
	}
	/* need to increment it so that this is the number of levels, not the index
	 * of the highest level */
	++nEner_H2_ground;

	/* this is the number of levels to do with the matrix - set with the
	 * atom h2 matrix command, keyword ALL means to do all of X in the matrix
	 * but number of levels within X was not known when the command was parsed,
	 * so this was set to -1 to defer setting to all until now */
	if( nXLevelsMatrix<0 )
	{
		nXLevelsMatrix = nLevels_per_elec[0];
	}
	else if( nXLevelsMatrix > nLevels_per_elec[0] )
	{
		fprintf( ioQQQ, 
			" The total number of levels used in the matrix solver was set to %li but there are only %li levels in X.\n Sorry.\n",
			nXLevelsMatrix ,
			nLevels_per_elec[0]);
		cdEXIT(EXIT_FAILURE);
	}

	/* at this stage the full electronic, vibration, and rotation energies have been defined,
	 * this is an option to print the energies */
	{
		/* set following true to get printout, false to not print energies */
		if( DEBUG_ENER )
		{
			/* print title for quantum numbers and energies */
			/*fprintf(ioQQQ,"elec\tvib\trot\tenergy\n");*/
			for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
			{
				fprintf(ioQQQ,"%li\t%li\t%li\t%.5e\n", (*st).n(), (*st).v(), (*st).J(), (*st).energy().WN() );
			}
			/* this will exit the program after printing the level energies */
			cdEXIT(EXIT_SUCCESS);
		}
	}

	/* this will prevent data from being read twice */
	lgREAD_DATA = true;

	/* create space for the electronic levels */
	H2_populations_LTE.reserve(n_elec_states);
	pops_per_vib.reserve(n_elec_states);
	H2_dissprob.reserve(n_elec_states);

	for( long iElec = 0; iElec<n_elec_states; ++iElec )
	{

		if( lgDebugPrt )
			fprintf(ioQQQ,"elec %li highest vib= %li\n", iElec , nVib_hi[iElec] );

		ASSERT( nVib_hi[iElec] > 0 );

		/* nVib_hi is now the highest vibrational level before dissociation,
		 * now allocate space to hold the number of rotation levels */
		H2_populations_LTE.reserve(iElec,nVib_hi[iElec]+1);
		pops_per_vib.reserve(iElec,nVib_hi[iElec]+1);

		/* now loop over all vibrational levels, and find out how many rotation levels there are */
		/* ground is special since use tabulated data - there are 14 vib states,
		 * ivib=14 is highest */
		for( long iVib = 0; iVib <= nVib_hi[iElec]; ++iVib )
		{
			/* lastly create the space for the rotation quantum number */
			H2_populations_LTE.reserve(iElec,iVib,nRot_hi[iElec][iVib]+1);
		}
	}

	H2_populations_LTE.alloc();
	H2_old_populations.alloc( H2_populations_LTE.clone() );
	H2_rad_rate_out.alloc( H2_populations_LTE.clone() );
	H2_dissprob.alloc( H2_populations_LTE.clone() );
	H2_disske.alloc( H2_populations_LTE.clone() );
	// this has a different geometry than the above
	pops_per_vib.alloc();

	H2_dissprob.zero();
	H2_disske.zero();

	/* set this one time, will never be set again, but might be printed */
	H2_rad_rate_out.zero();

	/* these do not have electronic levels - all within X */
	H2_ipPhoto.reserve(nVib_hi[0]+1);

	/* space for the vibration levels */
	for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
	{
		/* space for the rotation quantum number */
		H2_ipPhoto.reserve(iVib,nRot_hi[0][iVib]+1);
	}

	H2_ipPhoto.alloc();
	H2_col_rate_in.alloc( H2_ipPhoto.clone() );
	H2_col_rate_out.alloc( H2_ipPhoto.clone() );
	H2_rad_rate_in.alloc( H2_ipPhoto.clone() );
	H2_coll_dissoc_rate_coef.alloc( H2_ipPhoto.clone() );
	H2_coll_dissoc_rate_coef_H2.alloc( H2_ipPhoto.clone() );
	H2_X_colden.alloc( H2_ipPhoto.clone() );
	H2_X_rate_from_elec_excited.alloc( H2_ipPhoto.clone() );
	H2_X_rate_to_elec_excited.alloc( H2_ipPhoto.clone() );
	H2_X_colden_LTE.alloc( H2_ipPhoto.clone() );
	H2_X_formation.alloc( H2_ipPhoto.clone() );
	H2_X_Hmin_back.alloc( H2_ipPhoto.clone() );

	for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
	{
		for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
		{
			/* >>chng 04 jun 14, set these to bad numbers */
			H2_rad_rate_in[iVib][iRot] = -BIGFLOAT;
			H2_coll_dissoc_rate_coef[iVib][iRot] = -BIGFLOAT;
			H2_coll_dissoc_rate_coef_H2[iVib][iRot] = -BIGFLOAT;
		}
	}
	/* zero out the matrices */
	H2_X_colden.zero();
	H2_X_colden_LTE.zero();
	H2_X_formation.zero();
	H2_X_Hmin_back.zero();
	/* rates [cm-3 s-1] from elec excited states into X only vib and rot */
	H2_X_rate_from_elec_excited.zero();
	/* rates [s-1] to elec excited states from X only vib and rot */
	H2_X_rate_to_elec_excited.zero();

	/* distribution function for populations following formation from H minus H- */
	H2_X_hminus_formation_distribution.reserve(nTE_HMINUS);
	for( long i=0; i<nTE_HMINUS; ++i )
	{
		H2_X_hminus_formation_distribution.reserve(i,nVib_hi[0]+1);
		/* space for the vibration levels */
		for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
		{
			H2_X_hminus_formation_distribution.reserve(i,iVib,nRot_hi[0][iVib]+1);
		}
	}
	H2_X_hminus_formation_distribution.alloc();
	H2_X_hminus_formation_distribution.zero();
	H2_Read_hminus_distribution();

	/* >>chng 05 jun 20, do not use this, which is highly processed - use ab initio
	 * rates of excitation to electronic levels instead */
	/* read in cosmic ray distribution information
	H2_Read_Cosmicray_distribution(); */

	/* grain formation matrix */
	H2_X_grain_formation_distribution.reserve(H2_TOP);
	for( long ipH2=0; ipH2<(int)H2_TOP; ++ipH2 )
	{
		H2_X_grain_formation_distribution.reserve(ipH2,nVib_hi[0]+1);

		/* space for the vibration levels */
		for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
		{
			H2_X_grain_formation_distribution.reserve(ipH2,iVib,nRot_hi[0][iVib]+1);
		}
	}
	H2_X_grain_formation_distribution.alloc();
	H2_X_grain_formation_distribution.zero();

	for( long iElec=0; iElec<n_elec_states; ++iElec )
	{
		/* get dissociation probabilities and energies - ground state is stable */
		if( iElec > 0 )
			H2_ReadDissprob(iElec);
	}

	/* >>02 oct 18, add photodissociation, H2 + hnu => 2H + KE */
	/* we now have ro-vib energies, now set up threshold array offsets
	 * for photodissociation */
	for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
	{
		long iElec = (*st).n();
		if( iElec > 0 ) continue;
		long iVib = (*st).v();
		long iRot = (*st).J();
		/* this is energy needed to get up to n=3 electronic continuum 
		 * H2 cannot dissociate following absorption of a continuum photon into the
		 * continuum above X (which would require little energy) 
		 * because that process violates momentum conservation 
		 * these would be the triplet states - permitted are into singlets
		 * the effective full wavelength range of this process is from Lya to 
		 * Lyman limit in shielded regions 
		 * tests show limits are between 850A and 1220A - so Lya is included */
		/** \todo	1	add this as a Lya excitation process */
		/*>>KEYWORD	Allison & Dalgarno; continuum dissociation; */
		double thresh = (H2_DissocEnergies[1] - (*st).energy().WN())*WAVNRYD;
		/*fprintf(ioQQQ,"DEBUG\t%.2f\t%f\n", RYDLAM/thresh , thresh);*/
		/* in theory we should be able to assert that thesh just barely reaches
		 * lya, but actual numbers reach down to 0.749 ryd */
		ASSERT( thresh > 0.74 );
		H2_ipPhoto[iVib][iRot] = ipoint(thresh);
		fixit("this needs to be generalized");
	}

	CollRateCoeff.reserve( nLevels_per_elec[0] );
	for( long j = 1; j < nLevels_per_elec[0]; ++j )
	{
		CollRateCoeff.reserve( j, j );
		for( long k = 0; k < j; ++k )
		{
			CollRateCoeff.reserve( j, k, N_X_COLLIDER );
		}
	}
	CollRateCoeff.alloc();
	CollRateCoeff.zero();

	fixit("Does RateCoefTable only need to be N_X_COLLIDER long?");
	RateCoefTable.resize(ipNCOLLIDER);
	/* now read in the various sets of collision data */
	for( long nColl=0; nColl<N_X_COLLIDER; ++nColl )
	{
		/* ground state has tabulated data */
		H2_CollidRateRead(nColl);
	}

	CollRateErrFac.alloc( CollRateCoeff.clone() );
	CollRateErrFac.zero();

	/* option to add gaussian random mole */
	if( lgH2_NOISE )
	{
		for( long ipHi = 1; ipHi < nLevels_per_elec[0]; ++ipHi )
		{
			for( long ipLo = 0; ipLo < ipHi; ++ipLo )
			{
				for( long nColl=0; nColl<N_X_COLLIDER; ++nColl )
				{
					/* this returns the log of the random noise */
					realnum r = realnum(xMeanNoise + ran.normal()*xSTDNoise);
					CollRateErrFac[ipHi][ipLo][nColl] = exp10f(r);
				}
			}
		}
	}

	/* this will be total collision rate from an upper to a lower level within X */
	H2_X_source.resize( nLevels_per_elec[0] );
	H2_X_sink.resize( nLevels_per_elec[0] );

	H2_X_coll_rate.reserve(nLevels_per_elec[0]);
	/* now expand out to include all lower levels as lower state */
	for( long i=1; i<nLevels_per_elec[0]; ++i )
	{
		H2_X_coll_rate.reserve(i,i);
	}
	H2_X_coll_rate.alloc();

	ipTransitionSort.reserve( states.size() );
	for( unsigned nEner = 1; nEner < states.size(); ++nEner )
		ipTransitionSort.reserve( nEner, nEner );
	ipTransitionSort.alloc();
	ipTransitionSort.zero();
	lgH2_radiative.alloc( ipTransitionSort.clone() );
	lgH2_radiative.zero();
	
	trans.resize( (states.size() * (states.size()-1))/2 );
	AllTransitions.push_back(trans);
	qList initStates(label.c_str(),1);
	TransitionList initlist("H2InitList",&initStates);
	vector<TransitionList::iterator> initptrs;
	initlist.resize(trans.size());
	initlist.states() = &states;
	initptrs.resize(trans.size());

	{
		long lineIndex = 0;
		TransitionList::iterator tr = initlist.begin();
		for( unsigned ipHi=1; ipHi< states.size(); ++ipHi )
		{
			for( unsigned ipLo=0; ipLo<ipHi; ++ipLo )
			{
				(*tr).Junk();
				(*tr).setHi(ipHi);
				(*tr).setLo(ipLo);
				(*tr).Zero();
				initptrs[lineIndex] = tr;
				ipTransitionSort[ipHi][ipLo] = lineIndex;
				lineIndex++;
				++tr;
			}
		}
	}

	/* create the main array of lines */
	H2_SaveLine.reserve(n_elec_states);
	for( long iElecHi=0; iElecHi<n_elec_states; ++iElecHi )
	{
		H2_SaveLine.reserve(iElecHi,nVib_hi[iElecHi]+1);
		for( long iVibHi=0; iVibHi<=nVib_hi[iElecHi]; ++iVibHi )
		{
			H2_SaveLine.reserve(iElecHi,iVibHi,nRot_hi[iElecHi][iVibHi]+1);
			for( long iRotHi=Jlowest[iElecHi]; iRotHi<=nRot_hi[iElecHi][iVibHi]; ++iRotHi )
			{
				/* now the lower levels */
				/* NB - X is the only lower level considered here, since we are only 
				 * concerned with excited electronic levels as a photodissociation process
				 * code exists to relax this assumption - simply change following to iElecHi */
				long int lim_elec_lo = 0;
				H2_SaveLine.reserve(iElecHi,iVibHi,iRotHi,lim_elec_lo+1);
				for( long iElecLo=0; iElecLo<=lim_elec_lo; ++iElecLo )
				{
					H2_SaveLine.reserve(iElecHi,iVibHi,iRotHi,iElecLo,nVib_hi[iElecLo]+1);
					for( long iVibLo=0; iVibLo<=nVib_hi[iElecLo]; ++iVibLo )
					{
						H2_SaveLine.reserve(iElecHi,iVibHi,iRotHi,iElecLo,iVibLo,nRot_hi[iElecLo][iVibLo]+1);
					}
				}
			}
		}
	}

	H2_SaveLine.alloc();
	
	/* zero out array used to save emission line intensities */
	H2_SaveLine.zero();

	/* space for the energy vector is now allocated, must read trans probs from table */
	for( long iElec=0; iElec<n_elec_states; ++iElec )
	{
		/* ground state has tabulated data */
		H2_ReadTransprob(iElec,initlist);
	}

	// sort sys.trans so that radiative lines are at the beginning
	stable_sort( initptrs.begin(), initptrs.end(), compareEmis );
	rad_end = trans.begin();
	// rad_end will be used for the end of the range (non-inclusive) for operations on radiative lines
	{
		TransitionList::iterator tr = trans.begin();
		vector<TransitionList::iterator>::iterator ptr = initptrs.begin();
		for (size_t i=0; i < initptrs.size(); ++i, ++tr, ++ptr)
		{
			(*tr).copy(*(*ptr));
			if (lgRadiative(tr))
			{
				rad_end = tr;
			}
		}
	}
	++rad_end;
	ASSERT( rad_end != trans.end() ); 
	// after above sorting, ipTransitionSort is now invalid.  Fix here
	for( unsigned i = 0; i < trans.size(); ++i )
	{
		qList::iterator Hi = trans[i].Hi();
		qList::iterator Lo = trans[i].Lo();
		ipTransitionSort[ ipEnergySort[(*Hi).n()][(*Hi).v()][(*Hi).J()] ][ ipEnergySort[(*Lo).n()][(*Lo).v()][(*Lo).J()] ] = i;
		trans[i].ipHi() = ipEnergySort[(*Hi).n()][(*Hi).v()][(*Hi).J()];
		trans[i].ipLo() = ipEnergySort[(*Lo).n()][(*Lo).v()][(*Lo).J()];
	} 

	// link level and line stacks to species in the chemistry network	
	mole.species[ sp->index ].levels = &states;
	mole.species[ sp->index ].lines = &trans;

	// after above sorting, ipTransitionSort is now invalid.  Fix here
	for( unsigned i = 0; i < trans.size(); ++i )
	{
		qList::iterator Hi = trans[i].Hi();
		qList::iterator Lo = trans[i].Lo();
		ipTransitionSort[ ipEnergySort[(*Hi).n()][(*Hi).v()][(*Hi).J()] ][ ipEnergySort[(*Lo).n()][(*Lo).v()][(*Lo).J()] ] = i;
                //trans[i].ipHi() = ipEnergySort[(*Hi).n()][(*Hi).v()][(*Hi).J()];
                //trans[i].ipLo() = ipEnergySort[(*Lo).n()][(*Lo).v()][(*Lo).J()];
        }

	// this loop is over all transitions
	for( TransitionList::iterator tr = trans.begin(); tr != trans.end(); ++tr )
	{
		(*tr).EnergyWN() = (realnum)((*(*tr).Hi()).energy().WN() - (*(*tr).Lo()).energy().WN());
		/*wavelength of transition in Angstroms */
		if( (*tr).EnergyWN() > SMALLFLOAT)
			(*tr).WLAng() = (realnum)(1.e8f/(*tr).EnergyWN() / RefIndex( (*tr).EnergyWN() ) );

		(*tr).Coll().col_str() = 0.;
	}

	// this loop is over only radiative transitions.  Notice the end iterator.
	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		/* line redistribution function - will use complete redistribution */
		/* >>chng 04 mar 26, should include damping wings, especially for electronic
		 * transitions, had used doppler core only */
		(*tr).resetEmis();
		(*tr).Emis().iRedisFun() = ipCRDW;
		/* line optical depths in direction towards source of ionizing radiation */
		(*tr).Emis().TauIn() = opac.taumin;
		(*tr).Emis().TauCon() = opac.taumin;
		/* outward optical depth */
		(*tr).Emis().TauTot() = 1e20f;

		(*tr).Emis().dampXvel() = (realnum)( (*tr).Emis().Aul()/(*tr).EnergyWN()/PI4);
		(*tr).Emis().gf() = (realnum)(GetGF( (*tr).Emis().Aul(),(*tr).EnergyWN(), (*(*tr).Hi()).g() ) );

		/* derive the absorption coefficient, call to function is gf, wl (A), g_low */
		(*tr).Emis().opacity() = (realnum)( abscf( (*tr).Emis().gf(), (*tr).EnergyWN(), (*(*tr).Lo()).g()) );
		
		qList::iterator Hi = (*tr).Hi();	
		//qList::iterator Lo = (*tr).Lo();	

		if( (*Hi).n() == 0 )
		{
			/* the ground electronic state, most excitations are not direct pumping 
			 * (rather indirect, which does not count for ColOvTot) */
			(*tr).Emis().ColOvTot() = 1.;
		}
		else
		{
			/* these are excited electronic states, mostly pumped, except for supras */
			/** \todo	2	put supra thermal excitation into excitation of electronic bands */
			(*tr).Emis().ColOvTot() = 0.;
		}

		fixit("why not include this for excitations within X as well?");
		if( (*Hi).n() > 0 )
		{
			/* cosmic ray and non-thermal suprathermal excitation 
			 * to singlet state of H2 (B,C,B',D)
			 * cross section is equation 5 of 
			 *>>refer	H2	cs	Liu, W. & Dalgarno, A. 1994, ApJ, 428, 769
			 * relative to H I Lya cross section 
		 	 * this is used to derive H2 electronic excitations
			 * from the H I Lya rate
			 * the following is dimensionless scale factor for excitation 
			 * relative to H I Lya
			 */
			(*tr).Coll().col_str() = (realnum)(
				( (*tr).Emis().gf()/ (*tr).EnergyWN() ) /
				(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,0).Emis().gf()/
				iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,0).EnergyWN()) );
			ASSERT( (*tr).Coll().col_str()>0.);
		}
	}

	/*  Read continuum photodissociation cross section files */
	if( mole_global.lgStancil )
		Read_Mol_Diss_cross_sections();

	/* define branching ratios for deposition of H2 formed on grain surfaces,
	 * set true to use Takahashi distribution, false to use Draine & Bertoldi */

	/* loop over all types of grain surfaces */
	/* >>chng 02 oct 08, resolved grain types */
	/* number of different grain types H2_TOP is set in grainvar.h,
	 * types are ice, silicate, graphite */
	for( long ipH2=0; ipH2<(int)H2_TOP; ++ipH2 )
	{
		realnum sum = 0., sumj = 0., sumv = 0., sumo = 0., sump = 0.;
		
		/* first is Draine distribution */
		if( hmi.chGrainFormPump == 'D' )
		{
			long iElec = 0;
			/* H2 formation temperature, for equation 19, page 271, of
			* >>refer	H2	formation distribution	Draine, B.T., & Bertoldi, F., 1996, ApJ, 468, 269-289
			*/
			double T_H2_FORM = 50000.;
			for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
			{
				for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
				{
					/* no distinction between grain surface composition */
					H2_X_grain_formation_distribution[ipH2][iVib][iRot] = 
						/* first term is nuclear H2_stat weight */
						(1.f+2.f*H2_lgOrtho[iElec][iVib][iRot]) * (1.f+iVib) *
						(realnum)sexp( states[ ipEnergySort[iElec][iVib][iRot] ].energy().K()/T_H2_FORM );
					sum += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					sumj += iRot * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					sumv += iVib * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					if( H2_lgOrtho[iElec][iVib][iRot] )
					{
						sumo += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					}
					else
					{
						/* >>chng 02 nov 14, [0][iVib][iRot] -> [ipH2][iVib][iRot], PvH */
						sump += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					}
				}
			}
		}
		else if( hmi.chGrainFormPump == 'T' )
		{
			/* Takahashi 2001 distribution */
			double Xrot[H2_TOP] = { 0.14 , 0.15 , 0.15 };
			double Xtrans[H2_TOP] = { 0.12 , 0.15 , 0.25 };
			/* first normalize the vibration distribution function */
			double sumvib = 0.;
			double EH2;
			long iElec = 0;

			for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
			{
				double vibdist;
				EH2 = EH2_eval( ipH2, H2_DissocEnergies[0], states[ ipEnergySort[0][iVib][0] ].energy().WN() );
				vibdist = H2_vib_dist( ipH2 , EH2, H2_DissocEnergies[0], states[ ipEnergySort[0][iVib][0] ].energy().WN() );
				sumvib += vibdist;
			}
			/* this branch, use distribution function from
			* >>refer	grain	physics	Takahashi, Junko, 2001, ApJ, 561, 254-263 */
			for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
			{
				double Ev = states[ ipEnergySort[iElec][iVib][0] ].energy().WN()+energy_off;
				double Fv;
				/* equation 10 of Takahashi 2001, extra term is energy offset between bottom of potential
				 * the 0,0 level */
				double Erot;
				/*fprintf(ioQQQ," Evvvv\t%i\t%li\t%.3e\n", ipH2 ,iVib , Ev*WAVNRYD*EVRYD);*/

				EH2 = EH2_eval( ipH2, H2_DissocEnergies[0], states[ ipEnergySort[0][iVib][0] ].energy().WN() );

				/* equation 3 of Taktahashi & Uehara */
				Erot = (EH2 - Ev) * Xrot[ipH2] / (Xrot[ipH2] + Xtrans[ipH2]);

				/* email exchange with Junko Takahashi - 
				Thank you for your E-mail.
				I did not intend to generate negative Erot.
				I cut off the populations if their energy levels are negative, and made the total
				population be unity by using normalization factors (see, e.g., Eq. 12).

				I hope that my answer is of help to you and your work is going well.
				With best wishes,
				Junko

				>Thanks for the reply.  By cutting off the population, should we set the
				>population to zero when Erot becomes negative, or should we set Erot to
				>a small positive number? 

				I just set the population to zero when Erot becomes negative.
				Our model is still a rough one for the vibration-rotation distribution function
				of H2 newly formed on dust, because we have not yet had any exact
				experimental or theoretical data about it.
				With best wishes,
				Junko

				 */

				if( Erot > 0. )
				{
					/* the vibrational distribution */
					Fv = H2_vib_dist( ipH2 , EH2, H2_DissocEnergies[0], states[ ipEnergySort[0][iVib][0] ].energy().WN() ) / sumvib;
					/*fprintf(ioQQQ," vibbb\t%li\t%.3e\n", iVib , Fv );*/

					for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
					{
						double deltaE = states[ ipEnergySort[iElec][iVib][iRot] ].energy().WN() -
							states[ ipEnergySort[iElec][iVib][0] ].energy().WN();
						/* equation 6 of Takahashi 2001 */
						double gaussian = sexp( POW2( (deltaE - Erot) / (0.5 * Erot) ) );
						/* equation 7 of Takahashi 2001 */
						double thermal_dist = sexp( deltaE / Erot );

						/* take the mean of the two */
						double aver = ( gaussian + thermal_dist ) / 2.;
						/*fprintf(ioQQQ,"rottt\t%i\t%li\t%li\t%.3e\t%.3e\t%.3e\t%.3e\n",
							ipH2,iVib,iRot,
							deltaE*WAVNRYD*EVRYD,
							gaussian, thermal_dist , aver );*/

						/* thermal_dist does become > 1 since Erot can become negative */
						ASSERT( gaussian <= 1. /*&& thermal_dist <= 10.*/ );

						H2_X_grain_formation_distribution[ipH2][iVib][iRot] = (realnum)(
							/* first term is nuclear H2_stat weight */
							(1.f+2.f*H2_lgOrtho[iElec][iVib][iRot]) * Fv * (2.*iRot+1.) * aver );

						sum += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
						sumj += iRot * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
						sumv += iVib * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
						if( H2_lgOrtho[iElec][iVib][iRot] )
						{
							sumo += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
						}
						else
						{
							sump += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
						}

					}
				}
				else
				{
					/* this branch Erot is non-positive, so no distribution */
					for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
					{
						H2_X_grain_formation_distribution[ipH2][iVib][iRot] = 0.;
					}
				}
			}
		}
		else if( hmi.chGrainFormPump == 't' )
		{ 
			/* thermal distribution at 1.5 eV, as suggested by Amiel & Jaques */
			/* thermal distribution, upper right column of page 239 of
			 *>>refer	H2	formation	Le Bourlot, J, 1991, A&A, 242, 235 
			 * set with command
			 * set h2 grain formation pumping thermal */
			double T_H2_FORM = 17329.;
			long iElec = 0;
			for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
			{
				for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
				{
					/* no distinction between grain surface composition */
					H2_X_grain_formation_distribution[ipH2][iVib][iRot] = 
						/* first term is nuclear H2_stat weight */
						states[ ipEnergySort[0][iVib][iRot] ].g() *
						(realnum)sexp( states[ ipEnergySort[0][iVib][iRot] ].energy().K()/T_H2_FORM );
					sum += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					sumj += iRot * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					sumv += iVib * H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					if( H2_lgOrtho[iElec][iVib][iRot] )
					{
						sumo += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					}
					else
					{
						/* >>chng 02 nov 14, [0][iVib][iRot] -> [ipH2][iVib][iRot], PvH */
						sump += H2_X_grain_formation_distribution[ipH2][iVib][iRot];
					}
				}
			}
		}
		// ' ' is no formation pumping
		else if( hmi.chGrainFormPump == ' ' )
		{
			// neglect process
			sumo = sumj = sumv = 0.;
			sump = 1;
			for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
				for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )					/* no distinction between grain surface composition */
					H2_X_grain_formation_distribution[ipH2][iVib][iRot] = 0.;
			// put everything in lowest two levels
			H2_X_grain_formation_distribution[ipH2][0][0] =
				states[ ipEnergySort[0][0][0] ].g();
			H2_X_grain_formation_distribution[ipH2][0][1] =
			   states[ ipEnergySort[0][0][1] ].g();
			sum += H2_X_grain_formation_distribution[ipH2][0][0] +
					H2_X_grain_formation_distribution[ipH2][0][1];
		}
		else
			TotalInsanity();

		if( lgDebugPrt )
			fprintf(ioQQQ, "H2 form grains mean J= %.3f mean v = %.3f ortho/para= %.3f\n", 
				sumj/sum , sumv/sum , sumo/sump );

		//iElec = 0;
		/* now rescale so that integral is unity */
		for( long iVib = 0; iVib <= nVib_hi[0]; ++iVib )
		{
			for( long iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
			{
				H2_X_grain_formation_distribution[ipH2][iVib][iRot] /= sum;
				/* print the distribution function */
				/*if( states[ ipEnergySort[iElec][iVib][iRot] ].energy().WN() < 5200. )
				fprintf(ioQQQ,"disttt\t%i\t%li\t%li\t%li\t%.4e\t%.4e\t%.4e\t%.4e\n",
					ipH2, iVib , iRot, (long)H2_stat[0][iVib][iRot] ,
					states[ ipEnergySort[iElec][iVib][iRot] ].energy().WN(), 
					states[ ipEnergySort[iElec][iVib][iRot] ].energy().K(), 
					H2_X_grain_formation_distribution[ipH2][iVib][iRot],
					H2_X_grain_formation_distribution[ipH2][iVib][iRot]/H2_stat[0][iVib][iRot]
					);*/
			}
		}
	}

	return;
}
