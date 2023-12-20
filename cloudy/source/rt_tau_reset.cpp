/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_tau_reset after first iteration, updates the optical depths, mirroring this
 * routine but with the previous iteration's variables */
#include "cddefines.h"
#include "taulines.h"
#include "trace.h"
#include "iso.h"
#include "rfield.h"
#include "opacity.h"
#include "h2.h"
#include "geometry.h"
#include "dense.h"
#include "colden.h"
#include "rt.h"
#include "mole.h"

/* ====================================================================== */
/*RT_tau_reset update total optical depth scale, 
 * called by cloudy after iteration is complete */
void RT_tau_reset(void)
{
	long int ipHi, 
	  ipISO,
	  nelem, 
	  ipLo;

	DEBUG_ENTRY( "RT_tau_reset()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " UPDATE estimating new optical depths\n" );
		if( trace.lgHBug && trace.lgIsoTraceFull[ipH_LIKE] )
		{
			fprintf( ioQQQ, " New Hydrogen outward optical depths:\n" );
			for( ipHi=1; ipHi < iso_sp[ipH_LIKE][ trace.ipIsoTrace[ipH_LIKE] ].numLevels_max; ipHi++ )
			{
				fprintf( ioQQQ, "%3ld", ipHi );
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{
					if( iso_sp[ipH_LIKE][trace.ipIsoTrace[ipH_LIKE]].trans(ipHi,ipLo).ipCont() <= 0 )
						fprintf( ioQQQ, "%10.2e", 1e-30 );
					else
						fprintf( ioQQQ, "%10.2e", 
							iso_sp[ipH_LIKE][trace.ipIsoTrace[ipH_LIKE]].trans(ipHi,ipLo).Emis().TauIn() );
				}
				fprintf( ioQQQ, "\n" );
			}
		}
	}

	/* save column densities of H species */
	for( long i=0; i<NCOLD; ++i )
	{
		colden.colden_old[i] = colden.colden[i];
	}
	for( long i=0; i < mole_global.num_calc; i++ )
	{
		mole.species[i].column_old = mole.species[i].column;
	}

	opac.telec = opac.taumin;
	opac.thmin = opac.taumin;

	/* all iso sequences */
	for(ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				if( iso_ctrl.lgDielRecom[ipISO] )
				{
					for( ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
					{
						RT_line_one_tau_reset(SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][ipHi]]);
					}
				}

				for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					/* the bound-bound transitions within the model atoms */
					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						enum {DEBUG_LOC=false};
						if( DEBUG_LOC )
						{
							if( ipISO==1 && nelem==1 && ipHi==iso_ctrl.nLyaLevel[ipISO] && ipLo==0 )
								fprintf(ioQQQ,"DEBUG rt before loop %li %li %li %li tot %.3e in %.3e\n",
								ipISO, nelem, ipHi , ipLo , 
								iso_sp[ipISO][nelem].trans(iso_ctrl.nLyaLevel[ipISO],0).Emis().TauTot() ,
								iso_sp[ipISO][nelem].trans(iso_ctrl.nLyaLevel[ipISO],0).Emis().TauIn());
						}
						/*RT_line_one_tau_reset computes average of old and new optical depths 
						* for new scale at end of iter */
						RT_line_one_tau_reset(iso_sp[ipISO][nelem].trans(ipHi,ipLo));
						if( DEBUG_LOC )
						{
							if( ipISO==1 && nelem==1 && ipHi==iso_ctrl.nLyaLevel[ipISO] && ipLo==0 )
								fprintf(ioQQQ,"DEBUG rt after loop %li %li %li %li tot %.3e in %.3e\n",
								ipISO, nelem, ipHi , ipLo , 
								iso_sp[ipISO][nelem].trans(iso_ctrl.nLyaLevel[ipISO],0).Emis().TauTot() ,
								iso_sp[ipISO][nelem].trans(iso_ctrl.nLyaLevel[ipISO],0).Emis().TauIn());
						}
					}
				}
				/* the extra Lyman lines */
				/* Need all levels, as may have been raised/lowered throughout layer */
				for( ipHi=2; ipHi < iso_ctrl.nLyman_max[ipISO]; ipHi++ )
				{
					/* fully transfer all of the extra lines even though
					 * have not solved for their upper level populations */
					RT_line_one_tau_reset(ExtraLymanLines[ipISO][nelem][ipExtraLymanLines[ipISO][nelem][ipHi]]);
				}
			}
		}
	}

	/* >>>chng 99 nov 11 did not have case b for hydrogenic species on second and
	 * higher iterations */
	/* option to clobber these taus for Lyman lines, if case b is set */
	if( opac.lgCaseB )
	{
		for( nelem=0; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				realnum f;
				/* La may be case B, tlamin set to 1e9 by default with case b command */
				iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn() = opac.tlamin;
				/* >>>chng 99 nov 22, did not reset TauCon */
				iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauCon() = iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn();
				iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauTot() = 
				  2.f*iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn();
				f = opac.tlamin/iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().opacity();

				for( ipHi=3; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
				{
					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).ipCont() <= 0 )
						continue;

					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauIn() = 
						f*iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().opacity();
					/* reset line optical depth to continuum source */
					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauCon() = iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauIn();
					iso_sp[ipH_LIKE][nelem].trans(ipHi,ipH1s).Emis().TauTot() = 
						2.f*iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().TauIn();
				}
			}
		}

		/* now do helium like sequence - different since collapsed levels 
		 * all go to ground */
		for( nelem=1; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				double Aprev;
				realnum ratio;
				/* La may be case B, tlamin set to 1e9 by default with case b command */
				iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().TauIn() = opac.tlamin;

				iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().TauCon() = iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().TauIn();

				iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().TauTot() = 
				  2.f*iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().TauIn();

				ratio = opac.tlamin/iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().opacity();

				/* this will be the trans prob of the previous lyman line, will use this to 
				 * find the next one up in the series */
				Aprev = iso_sp[ipHE_LIKE][nelem].trans(ipHe2p1P,ipHe1s1S).Emis().Aul();
				/* >>chng 02 jan 05, remove explicit list of lyman lines, use As to guess
				 * which are which - this will work for any number of levels */
				for( long i = ipHe2p1P+1; i < iso_sp[ipHE_LIKE][nelem].numLevels_max; i++ )
				{
					if( iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).ipCont() <= 0 )
						continue;

					/* >>chng 02 mar 19 use proper test for resonance collapsed lines */
					if( iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().Aul()> Aprev/10. ||
						iso_sp[ipHE_LIKE][nelem].st[i].S() < 0 )
					{
						Aprev = iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().Aul();
						iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauIn() = 
							ratio*iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().opacity();
						/* reset line optical depth to continuum source */
						iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauCon() = 
							iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauIn();
						iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauTot() = 
							2.f*iso_sp[ipHE_LIKE][nelem].trans(i,ipHe1s1S).Emis().TauIn();
					}
				}
			}
		}
	}

	/* zero out fine opacity fine grid fine mesh array */
	vzero(rfield.fine_opt_depth);
	// Make sure stale zone opacities don't wrap around between
	// iterations.  If this /does/ have an effect, it means that the
	// opacities being used in are one zone stale.  Likely touch points
	// are line overlap & radiation pressure calculations.
	// So zero this isn't ideal, just better than using whatever the
	// value happens to be in the last zone of the previous iteration...
	vzero(rfield.fine_opac_zone);

	/* all level 2 heavy element lines */
	for( long i=0; i < nWindLine; i++ )
	{
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
		{
			RT_line_one_tau_reset(TauLine2[i]);
		}
	}

	/* all hyperfine structure lines */
	for( size_t i=0; i < HFLines.size(); i++ )
	{
		RT_line_one_tau_reset(HFLines[i]);
	}

	/* external database lines */
	for (int ipSpecies=0; ipSpecies < nSpecies; ++ipSpecies)
	{
		for( EmissionList::iterator em=dBaseTrans[ipSpecies].Emis().begin();
			  em != dBaseTrans[ipSpecies].Emis().end(); ++em)
		{
			RT_line_one_tau_reset((*em).Tran());
		}
	}

	/* inner shell lines */
	for( size_t i=0; i < UTALines.size(); i++ )
	{
		/* these are line optical depth arrays
		 * inward optical depth */
		/* heat is special for this array - it is heat per pump */
		double hsave = UTALines[i].Coll().heat();
		RT_line_one_tau_reset(UTALines[i]);
		UTALines[i].Coll().heat() = hsave;
	}

	/* the large H2 molecule */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_RT_tau_reset();

	if( opac.lgCaseB )
	{
		for( long i=0; i < rfield.nflux_with_check; i++ )
		{
			/* DEPABS and SCT are abs and sct optical depth for depth only
			 * we will not change total optical depths, just reset inner to half
			 * TauAbsGeo(i,2) = 2.*TauAbsFace(i) */
			opac.TauAbsGeo[0][i] = opac.TauAbsGeo[1][i]/2.f;
			/* TauScatGeo(i,2) = 2.*TauScatFace(i) */
			opac.TauScatGeo[0][i] = opac.TauScatGeo[1][i]/2.f;
		}
	}
	else if( geometry.lgSphere )
	{
		/* closed or spherical case */
		for( long i=0; i < rfield.nflux_with_check; i++ )
		{
			/* [1] is total optical depth from previous iteration,
			 * [0] is optical depth at current position */
			opac.TauAbsGeo[1][i] = 2.f*opac.TauAbsFace[i];
			opac.TauAbsGeo[0][i] = opac.TauAbsFace[i];
			opac.TauScatGeo[1][i] = 2.f*opac.TauScatFace[i];
			opac.TauScatGeo[0][i] = opac.TauScatFace[i];
			opac.TauTotalGeo[1][i] = opac.TauScatGeo[1][i] + opac.TauAbsGeo[1][i];
			opac.TauTotalGeo[0][i] = opac.TauScatGeo[0][i] + opac.TauAbsGeo[0][i];
		}
	}
	else
	{
		/* open geometry */
		for( long i=0; i < rfield.nflux_with_check; i++ )
		{
			opac.TauTotalGeo[1][i] = opac.TauTotalGeo[0][i];
			opac.TauTotalGeo[0][i] = opac.taumin;
			opac.TauAbsGeo[1][i] = opac.TauAbsGeo[0][i];
			opac.TauAbsGeo[0][i] = opac.taumin;
			opac.TauScatGeo[1][i] = opac.TauScatGeo[0][i];
			opac.TauScatGeo[0][i] = opac.taumin;
		}
	}

	/* same for open and closed geometries */
	for( long i=0; i < rfield.nflux_with_check; i++ )
	{
		/* total optical depth across computed shell */
		opac.TauAbsTotal[i] = opac.TauAbsFace[i];
		/* e2( tau across shell), optical depth from ill face to shielded face of cloud 
		 * not that opac.TauAbsFace is reset to small number just after this */
		opac.E2TauAbsTotal[i] = (realnum)e2( opac.TauAbsTotal[i] );
		/* TauAbsFace and TauScatFace are abs and sct optical depth to ill face */
		opac.TauScatFace[i] = opac.taumin;
		opac.TauAbsFace[i] = opac.taumin;
	}

	/* this is optical depth at x-ray point defining effective optical depth */
	rt.tauxry = opac.TauAbsGeo[0][rt.ipxry-1];
	return;
}
