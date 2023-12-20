/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*OpacityAddTotal derive total opacity for this position,
 * called by ConvBase */
#include "cddefines.h"
#include "iso.h"
#include "ipoint.h"
#include "grainvar.h"
#include "ca.h"
#include "rfield.h"
#include "oxy.h"
#include "h2.h"
#include "hmi.h"
#include "atoms.h"
#include "conv.h"
#include "ionbal.h"
#include "trace.h"
#include "phycon.h"
#include "opacity.h"
#include "mole.h"
#include "freebound.h"
#include "dense.h"
#include "atmdat.h"
#include "atmdat_gaunt.h"
#include "cool_eval.h"

void OpacityAddTotal(void)
{
	long int i, 
	  ion,
	  ipHi,
	  ipISO,
	  ipop,
	  limit,
	  low,
	  n; 
	double DepartCoefInv ,
	  fac, 
	  sum;
	realnum SaveOxygen1 , 
		SaveCarbon1;

	DEBUG_ENTRY( "OpacityAddTotal()" );

	/* OpacityZero will zero out scattering and absorption opacities,
	 * and set OldOpacSave to opac to save it */
	OpacityZero();

	/* free electron scattering opacity, Compton recoil energy loss */
	for( i=0; i < rfield.nflux; i++ )
	{
		/* scattering part of total opacity */
		opac.opacity_sct[i] += opac.OpacStack[i-1+opac.iopcom]*
		  dense.eden;
	}

	/* opacity due to Compton bound recoil ionization */
	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{ 
			for( ion=0; ion<nelem+1; ++ion )
			{
				realnum factor = dense.xIonDense[nelem][ion];
				/*>>chng 05 nov 26, add molecular hydrogen assuming same
				 * as two free hydrogen atoms - hence mult density by two
				 *>>KEYWORD	H2 bound Compton ionization */
				if( nelem==ipHYDROGEN )
					factor += hmi.H2_total*2.f;
				if( factor > 0. )
				{
					// loop_min and loop_max are needed to work around a bug in icc 10.0
					long loop_min = ionbal.ipCompRecoil[nelem][ion]-1;
					long loop_max = rfield.nflux;
					/* ionbal.nCompRecoilElec number of electrons in valence shell
					 * that can compton recoil ionize */
					factor *= ionbal.nCompRecoilElec[nelem-ion];
					for( i=loop_min; i < loop_max; i++ )
					{
						/* add in bound hydrogen electron scattering, treated as absorption */
						opac.opacity_abs[i] += opac.OpacStack[i-1+opac.iopcom]*factor;
					}
				}
			}
		}
	}

	/* opacity due to pair production - does not matter what form these
	 * elements are in */
	/** \todo	2	add charged heavy elements */
	sum = dense.gas_phase[ipHYDROGEN] + 4.*dense.gas_phase[ipHELIUM];
	OpacityAdd1Subshell(opac.ioppr,opac.ippr,rfield.nflux,(realnum)sum,'s');

	// only update ee brems spectrum when temperature changes
	// the electron density dependence is explicitly factored in further below
	if( !fp_equal( phycon.te, opac.eeFreeFreeTemp ) )
	{
		eeBremsSpectrum( phycon.te, rfield.eeBremsDif, opac.eeFreeFreeOpacity );
		opac.eeFreeFreeTemp = phycon.te;
	}

	/* free free brems emission for all ions */

	/* Precomputed opacity value in OpacStack is inflated by 30 dex to avoid underflows;
	 * deflation by 1e-30 corrects for that. */
	double bfac = 1e-30 * dense.eden / phycon.sqrte;
	/* >>chng 02 jul 21, use full set of ions and gaunt factor */
	/* gaunt factors depend only on photon energy and ion charge, so do
	 * sum of ions here before entering into loop over photon energy */
	t_brems_den bsum;
	t_gaunt::Inst().brems_sum_ions(bsum);
	// there is no need to keep these separate here...
	bsum.den_ion[1] += bsum.den_Hp + bsum.den_Hep;
	bsum.den_ion[2] += bsum.den_Hepp;

	for( i=0; i < rfield.nflux; i++ )
		opac.FreeFreeOpacity[i] = 0.;

	t_gaunt::Inst().brems_opac( -1, phycon.te, bfac*bsum.den_Hm, opac.FreeFreeOpacity );
	for( ion=1; ion < LIMELM+1; ++ion )
		if( bsum.den_ion[ion] > 0. )
			t_gaunt::Inst().brems_opac( ion, phycon.te, bfac*bsum.den_ion[ion], opac.FreeFreeOpacity );

	for( i=0; i < rfield.nflux; i++ )
	{
		if( rfield.ContBoltz[i] < 0.995 )
			opac.FreeFreeOpacity[i] *= opac.OpacStack[i-1+opac.ipBrems]*(1. - rfield.ContBoltz[i]);
		else
			opac.FreeFreeOpacity[i] *= opac.OpacStack[i-1+opac.ipBrems]*rfield.anu(i)*TE1RYD/phycon.te;
		opac.opacity_abs[i] += opac.FreeFreeOpacity[i] + pow2(dense.eden)*opac.eeFreeFreeOpacity[i];
	}

	/* H minus absorption, with correction for stimulated emission */
	if( hmi.hmidep > SMALLFLOAT )
	{
		DepartCoefInv = 1./hmi.hmidep;
	}
	else
	{
		/* the hmidep departure coef can become vastly small in totally
		 * neutral gas (no electrons) */
		DepartCoefInv = 1.;
	}
	limit = ipoint(atmdat.EIonPot[ipHELIUM][0]);

	const molezone *MHm = findspecieslocal("H-"); 
	for( i=hmi.iphmin-1; i < limit; i++ )
	{
		double factor;
		factor = 1. - rfield.ContBoltz[i]*DepartCoefInv;
		if(factor > 0)
			opac.opacity_abs[i] += opac.OpacStack[i-hmi.iphmin+opac.iphmop]*
				MHm->den*factor;
	}

	/* H2P h2plus bound free opacity */
	{
		const molezone *MH2p = findspecieslocal("H2+");
		limit = opac.ih2pnt[1]; 
		if(limit > rfield.nflux)
			limit = rfield.nflux;
		double frac = MAX2( 0., 1.-hmi.h2plus_exc_frac );
		for( i=opac.ih2pnt[0]-1; i < limit; i++ )
		{
			opac.opacity_abs[i] += MH2p->den*frac*opac.OpacStack[i-opac.ih2pnt[0]+
			  opac.ih2pof];
			ASSERT( !isnan( opac.opacity_abs[i] ) );
		}

		frac = hmi.h2plus_exc_frac;
		limit = opac.ih2pnt_ex[1]; 
		if(limit > rfield.nflux)
			limit = rfield.nflux;
		for( i=opac.ih2pnt_ex[0]-1; i < limit; i++ )
		{
			opac.opacity_abs[i] += MH2p->den*frac*opac.OpacStack[i-opac.ih2pnt_ex[0]+
			  opac.ih2pof_ex];
			ASSERT( !isnan( opac.opacity_abs[i] ) );
		}
	}

	/* H2 continuum dissociation opacity */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
	{
		if( (*diatom)->lgEnabled && mole_global.lgStancil )
		{
			for( vector< diss_tran >::iterator tran = (*diatom)->Diss_Trans.begin(); tran != (*diatom)->Diss_Trans.end(); ++tran )
			{
				long lower_limit = ipoint(tran->energies[0]);
				long upper_limit = ipoint(tran->energies.back());
				upper_limit = MIN2( upper_limit, rfield.nflux-1 );
				for(i = lower_limit; i <= upper_limit; ++i)
				{
					opac.opacity_abs[i] += (*diatom)->MolDissocOpacity( *tran, rfield.anu(i));
				}
			}
		}
	}

	/* get total population of hydrogen ground to do Rayleigh scattering */
	if( dense.xIonDense[ipHYDROGEN][1] <= 0. )
	{
		fac = dense.xIonDense[ipHYDROGEN][0];
	}
	else
	{
		fac = iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop();
	}

	/* Ly a damp wing opac (Rayleigh scattering) */
	limit = iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon;
	if(limit > rfield.nflux)
		limit = rfield.nflux;
	for( i=0; i < limit; i++ )
	{
		opac.opacity_sct[i] += (fac*opac.OpacStack[i-1+opac.ipRayScat]);
	}

	 /* remember largest correction for stimulated emission */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].DepartCoef() > 1e-30 && !conv.lgSearch )
	{
		realnum factor;
		factor = (realnum)(rfield.ContBoltz[iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1]/iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].DepartCoef());
		if(opac.stimax[0] < factor)
			opac.stimax[0] = factor;
	}

	if( iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].DepartCoef() > 1e-30 && !conv.lgSearch )
	{
		realnum factor;
		factor = (realnum)(rfield.ContBoltz[iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH2p].ipIsoLevNIonCon-1]/iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].DepartCoef());
		if(opac.stimax[1] < factor)
			opac.stimax[1] = factor;
	}

#	if 0
	/* check whether hydrogen or Helium singlets mased, if not in search mode */
	if( !conv.lgSearch )
	{
		if( iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).PopOpc() < 0. )
		{
			hydro.lgHLyaMased = true;
		}
	}
#	endif

	/* >>chng 05 nov 25, use Yan et al. H2 photo cs
	 * following reference gives cross section for all energies
	 * >>refer	H2	photo cs	Yan, M., Sadeghpour, H.R., & Dalgarno, A., 1998, ApJ, 496, 1044 
	 * Wilms, J., Allen, A., & McCray, R. 2000, ApJ, 542, 914 */
	/* >>chng 02 jan 16, approximate inclusion of H_2 photoelectric opacity 
	 * include H_2 in total photoelectric opacity as twice H0 cs */
	/* set lower and upper limits to this range */
	/*>>KEYWORD	H2	photoionization opacity */
	low = h2.ip_photo_opac_thresh;
	ipHi = rfield.nflux_with_check;
	ipop = h2.ip_photo_opac_offset;
	/* OpacityAdd1Subshell just returns for static opacities if opac.lgRedoStatic not set*/
	/* >>chng 05 nov 27, change on nov 25 had left 2*density from H0, so
	 * twice the H2 density was used - 	 * also changed to static opacity 
	 * this assumes that all v,J levels contribute the same opacity */
	OpacityAdd1Subshell( ipop , low , ipHi , hmi.H2_total , 's' );

	/*>>KEYWORD	CO photoionization opacity */
	/* include photoionization of CO - assume C and O in CO each have 
	 * same photo cs as atom - this should only be significant in highly
	 * shielded regions where only very hard photons penetrate 
	 * also H2O condensed onto grain surfaces - very important deep in cloud */
	SaveOxygen1 = dense.xIonDense[ipOXYGEN][0];
	SaveCarbon1 = dense.xIonDense[ipCARBON][0];
	/* atomic C and O will include CO during the heating sum loop */
	fixit("CO opacity hack breaks the invariant for total mass");
						/*
						* In any case, is CO the only species which contributes?
						* -- might expect all other molecular species to do so,
						* i.e. the item added should perhaps be
						* dense.xMolecules[nelem]) -- code duplicated in heatsum
						*/
	dense.xIonDense[ipOXYGEN][0] += (realnum) (findspecieslocal("CO")->den + findspecieslocal("H2Ogrn")->den);
	dense.xIonDense[ipCARBON][0] += (realnum) (findspecieslocal("CO")->den);

	/* following loop adds standard opacities for first 30 elements 
	 * most heavy element opacity added here */
	for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		/* this element may be turned off */
		if( dense.lgElmtOn[nelem] )
		{ 
			OpacityAdd1Element(nelem);
		}
	}

	/* now reset the abundances */
	dense.xIonDense[ipOXYGEN][0] = SaveOxygen1;
	dense.xIonDense[ipCARBON][0] = SaveCarbon1;

	/* following are opacities due to specific excited levels */

	/* nitrogen opacity
	 * excited level of N+ */
	OpacityAdd1Subshell(opac.in1[2],opac.in1[0],opac.in1[1],
		dense.xIonDense[ipNITROGEN][0]*atoms.p2nit , 'v' );

	/* oxygen opacity
	 * excited level of Oo */
	OpacityAdd1Subshell(opac.ipo1exc[2],opac.ipo1exc[0],opac.ipo1exc[1],
	  dense.xIonDense[ipOXYGEN][0]*oxy.poiexc,'v');

	/* O2+ excited states */
	OpacityAdd1Subshell(opac.ipo3exc[2],opac.ipo3exc[0],opac.ipo3exc[1],
	  dense.xIonDense[ipOXYGEN][2]*oxy.poiii2,'v');

	OpacityAdd1Subshell(opac.ipo3exc3[2],opac.ipo3exc3[0],opac.ipo3exc3[1],
	  dense.xIonDense[ipOXYGEN][2]*oxy.poiii3,'v');

	/* magnesium opacity
	 * excited level of Mg+ */
	OpacityAdd1Subshell(opac.ipOpMgEx,opac.ipmgex,iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon,
		dense.xIonDense[ipMAGNESIUM][1]* atoms.popMg2,'v');

	/* calcium opacity
	 * photoionization of excited levels of Ca+ */
	OpacityAdd1Subshell(opac.ica2op,opac.ica2ex[0],opac.ica2ex[1],
	  ca.popca2ex,'v');

	/*******************************************************************
	 *
	 * complete evaluation of total opacity by adding in the static part and grains
	 *
	 *******************************************************************/

	/* this loop defines the variable iso_sp[ipH_LIKE][nelem].fb[n].ConOpacRatio,
	 * the ratio of not H to Hydrogen opacity.  for grain free environments
	 * at low densities this is nearly zero.  The correction includes 
	 * stimulated emission correction */
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* this element may be turned off */
			if( dense.lgElmtOn[nelem] )
			{ 
				/* this branch is for startup only */
				if( nzone < 1 )
				{
					/* >>chng 06 aug 17, should go to numLevels_local instead of _max */
					for( n=0; n < iso_sp[ipISO][nelem].numLevels_local; n++ )
					{
						if(iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon < rfield.nflux )
						{
							/*>>chng 04 dec 12, had tested against 1e-30, set to zero if not
							 * greater than this - caused oscillations as opacity fell below
							 * and around this value - change to greater than 0 */
							/*if( opac.opacity_abs[iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon-1] > 1e-30 )*/
							if( opac.opacity_abs[iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon-1] > 0. )
							{
								/* >>chng 02 may 06, use general form of threshold cs */
								/*double t1 = atmdat_sth(n)/(POW2(nelem+1.-ipISO));*/
								long int ip = iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon;
								double t2 = csphot(
									ip ,
									ip ,
									iso_sp[ipISO][nelem].fb[n].ipOpac );

								iso_sp[ipISO][nelem].fb[n].ConOpacRatio = 1.-
									(iso_sp[ipISO][nelem].st[n].Pop()*t2)/
									opac.opacity_abs[ip-1];
							}
							else
							{
								iso_sp[ipISO][nelem].fb[n].ConOpacRatio = 0.;
							}
						}
					}
				}
				/* end branch for startup only, start branch for all zones including startup */
				/* >>chng 06 aug 17, should go to numLevels_local instead of _max */
				for( n=0; n < iso_sp[ipISO][nelem].numLevels_local; n++ )
				{
					/* ratios of other to total opacity for continua of all atoms done with iso model */
					if(iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon < rfield.nflux )
					{
						/*>>chng 04 dec 12, had tested against 1e-30, set to zero if not
						 * greater than this - caused oscillations as opacity fell below
						 * and around this value - change to greater than 0 */
						/*if( opac.opacity_abs[iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon-1] > 1e-30 )*/
						if( opac.opacity_abs[iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon-1] > 0. )
						{
							/* first get departure coef */
							if( iso_sp[ipISO][nelem].st[n].DepartCoef() > 1e-30 && (!conv.lgSearch ) )
							{
								/* this is the usual case, use inverse of departure coef */
								fac = 1./iso_sp[ipISO][nelem].st[n].DepartCoef();
							}
							else if( conv.lgSearch )
							{
								/* do not make correction for stim emission during search
								* for initial temperature solution, since trys are very non-equil */
								fac = 0.;
							}
							else
							{
								fac = 1.;
							}

							/** \todo	1	stupid - why this test on opacity_abs ? - we only get here
							 * if we already passed above test on this very thing */
							/* now get opaicty ratio with correction for stimulated emission */
							/*>>chng 04 dec 12, had tested against 1e-30, set to zero if not
							 * greater than this - caused oscillations as opacity fell below
							 * and around this value - change to greater than 0 */
							/*if( opac.opacity_abs[iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon-1] > 1e-30 )*/
							if( opac.opacity_abs[iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon-1] > 0. )
							{
								/* >>chng 02 may 06, use general form of threshold cs */
								long int ip = iso_sp[ipISO][nelem].fb[n].ipIsoLevNIonCon;

								double t2 = csphot(
									ip ,
									ip ,
									iso_sp[ipISO][nelem].fb[n].ipOpac );

								double opacity_this_species = 
									iso_sp[ipISO][nelem].st[n].Pop()*t2*
									(1. - fac*rfield.ContBoltz[ip-1]);

								double opacity_fraction =  1. - opacity_this_species / opac.opacity_abs[ip-1];
								if(opacity_fraction < 0)
									opacity_fraction = 0.;

								/* use mean of old and new ratios */
								iso_sp[ipISO][nelem].fb[n].ConOpacRatio = 
									iso_sp[ipISO][nelem].fb[n].ConOpacRatio* 0.75 + 0.25*opacity_fraction;

								if(iso_sp[ipISO][nelem].fb[n].ConOpacRatio < 0.)
									iso_sp[ipISO][nelem].fb[n].ConOpacRatio = 0.;
							}
							else
							{
								iso_sp[ipISO][nelem].fb[n].ConOpacRatio = 0.;
							}
						}
						else
						{
							iso_sp[ipISO][nelem].fb[n].ConOpacRatio = 0.;
						}
					}
					else
					{
						iso_sp[ipISO][nelem].fb[n].ConOpacRatio = 0.;
					}
				}
			}
		}
	}

	/* add dust grain opacity if dust present */
	if( gv.lgDustOn() )
	{
		/* generate current grain opacities since may be function of depth */
		/* >>chng 01 may 11, removed code to update grain opacities, already done by GrainChargeTemp */
		for( i=0; i < rfield.nflux; i++ )
		{
			/* units cm-1 */
			opac.opacity_sct[i] += gv.dstsc[i]*dense.gas_phase[ipHYDROGEN];
			opac.opacity_abs[i] += gv.dstab[i]*dense.gas_phase[ipHYDROGEN];
		}
	}

	/* check that opacity is sane */
	for( i=0; i < rfield.nflux; i++ )
	{
		/* OpacStatic was zeroed in OpacityZero, incremented in opacityadd1subshell */
		opac.opacity_abs[i] += opac.OpacStatic[i];
		/* make sure that opacity is positive */
		/*ASSERT( opac.opacity_abs[i] > 0. );*/
	}

	/* compute gas albedo here */
	for( i=0; i < rfield.nflux; i++ )
	{
		opac.albedo[i] = opac.opacity_sct[i]/
			(opac.opacity_sct[i] + opac.opacity_abs[i]);
	}

	/* during search phase set old opacity array to current value */
	if( conv.lgSearch )
		OpacityZeroOld();

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "     OpacityAddTotal returns; grd rec eff (opac) for Hn=1,4%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e", 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ConOpacRatio, 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH2s].ConOpacRatio, 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH2p].ConOpacRatio, 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH3s].ConOpacRatio, 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH3p].ConOpacRatio, 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH3d].ConOpacRatio, 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH4s].ConOpacRatio, 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH4p].ConOpacRatio, 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH4d].ConOpacRatio, 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH4f].ConOpacRatio );
		if( dense.lgElmtOn[ipHELIUM] )
		{
			fprintf( ioQQQ, " HeI,II:%10.2e%10.2e",
					 iso_sp[ipHE_LIKE][ipHELIUM].fb[ipHe1s1S].ConOpacRatio,
					 iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].ConOpacRatio );
		}
		fprintf( ioQQQ, "\n" );
	}

	{
		/* following should be set true to print strongest ots contributors */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && (nzone>=378)/**/ )
		{
			if( nzone > 380 ) 
				cdEXIT( EXIT_FAILURE );
			for( i=0; i<rfield.nflux; ++i )
			{
				fprintf(ioQQQ,"rtotsbugggg\t%li\t%.3e\t%.3e\t%.3e\t%.3e\n",
					conv.nPres2Ioniz,
					rfield.anu(i),
					opac.opacity_abs[i],
					rfield.otscon[i],
					rfield.otslin[i]);
			}
		}
	}
	return;
}
