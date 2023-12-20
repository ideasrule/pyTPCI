/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_line_all do escape and destruction probabilities for all lines in code. */
#include "cddefines.h"
#include "taulines.h"
#include "dense.h"
#include "conv.h"
#include "rfield.h"
#include "wind.h"
#include "iso.h"
#include "h2.h"
#include "opacity.h"
#include "trace.h"
#include "hydrogenic.h"
#include "rt.h"
#include "cosmology.h"
#include "doppvel.h"
#include "atmdat.h"

/*RT_line_all do escape and destruction probabilities for all lines in code. */
void RT_line_all_escape( realnum *error )
{
	DEBUG_ENTRY( "RT_line_all_escape()" );
	/* this flag says whether to update the line escape probabilities */
	bool lgPescUpdate = conv.lgFirstSweepThisZone || conv.lgIonStageTrimed;
	/* find Stark escape probabilities for hydrogen itself */
	if( lgPescUpdate )
		RT_stark();

	bool lgCheck = (error != NULL);
	vector<realnum> oldPesc, oldPdest, oldPelec_esc;
	if (lgCheck)
	{
		for( vector<TransitionList>::iterator it = AllTransitions.begin(); it != AllTransitions.end(); ++it )
		{
			for (TransitionList::iterator tr=it->begin();
				  tr != it->end(); ++tr)
			{
				if ((*tr).ipCont() <= 0)
					continue;
				oldPdest.push_back((*tr).Emis().Pdest());
				oldPesc.push_back((*tr).Emis().Pesc());
				oldPelec_esc.push_back((*tr).Emis().Pelec_esc());
			}
		}
	}
	

	RT_line_all(RT_line_one_escape);

	if( opac.lgCaseB_no_pdest )
	{
		for( long ipISO=ipH_LIKE; ipISO < NISO; ++ipISO )
		{
			/* loop over all iso-electronic sequences */
			for( long nelem=ipISO; nelem < LIMELM; ++nelem )
			{
				/* parent ion stage, for H is 1, for He is 1 for He-like and 
				 * 2 for H-like */
				long ion = nelem+1-ipISO;
				
				/* element turned off */
				if( !dense.lgElmtOn[nelem] )
					continue;
				/* need we consider this ion? */
				if( ion <= dense.IonHigh[nelem] )
				{
					/* loop over all lines */
					/* this is option to not do destruction probabilities
					 * for case b no pdest option */
					long ipLo = 0;
					/* okay to let this go to numLevels_max. */
					for( long ipHi=ipLo+1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ++ipHi )
					{
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() <= 0 )
							continue;
						
						/* hose the previously computed destruction probability */
						iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pdest() = SMALLFLOAT; 
					}
				}
			}
		}
	}

	/* note that pesc and dest are updated no matter what escprob logic we
	 * specify, and that these are not updated if we have overrun the
	 * optical depth scale.  only update here in that case */
	// Logic must be kept consistent with mask for RT_line_escape in rt_line_one.cpp
	if( lgTauGood( iso_sp[ipH_LIKE][ipHYDROGEN].trans(iso_ctrl.nLyaLevel[ipH_LIKE],0)) &&
		  !conv.lgLastSweepThisZone )
	{
		/*fprintf(ioQQQ,"DEBUG fe2 %.2e %.2e\n", hydro.dstfe2lya ,
			hydro.HLineWidth);*/

		/* >>chng 00 jan 06, let dest be large than 1 to desaturate the atom */
		/* >>chng 01 apr 01, add test for tout >= 0., 
		 * must not add to Pdest when it has not been refreshed here */
		iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().Pdest() += hydro.dstfe2lya;
	}

	/* is continuum pumping of H Lyman lines included?  yes, but turned off
	 * with atom h-like Lyman pumping off command */
	if( !hydro.lgLymanPumping )
	{
		long ipISO = ipH_LIKE;
		long nelem = ipHYDROGEN;
		long ipLo = 0;
		for( long ipHi=ipLo+1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ++ipHi )
		{
			iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().pump() = 0.;
		}
	}

	{
		/* following should be set true to print ots contributors for he-like lines*/
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC && nzone>433 /*&& iteration==2*/ )
		{
			/* option to dump a line  */
			DumpLine(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,0) );
		}
	}
#	if 0
	/* this code is a copy of the code within line_one that does cloaking
	 * within this zone.  it can be used to see how a particular line
	 * is being treated. */
	{
#include "doppvel.h"
		double dTau , aa;
		ipISO = 0; nelem = 0;ipLo = 0;
		ipHi = 23;
		dTau =  iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().PopOpc() * 
			iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().opacity() / 
			GetDopplerWidth(dense.AtomicWeight[nelem]) + opac.opacity_abs[iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont()-1];
		dTau *= radius.drad_x_fillfac;
		aa = log(1. + dTau ) / SDIV(dTau);
		fprintf(ioQQQ,"DEBUG dTau\t%.2f\t%.5e\t%.5e\t%.5e\n",fnzone,dTau,
			radius.drad_x_fillfac,
			 aa);
	}
#	endif

	realnum bigerror=0.0;
	TransitionList::iterator trworst;
	long iworst=-1,cworst=-1;
	if (lgCheck)
	{
		long ind=0;
		for( vector<TransitionList>::iterator it = AllTransitions.begin(); it != AllTransitions.end(); ++it )
		{
			for (TransitionList::iterator tr=it->begin();
				  tr != it->end(); ++tr)
			{
				if ((*tr).ipCont() <= 0)
					continue;
				
				const realnum abslim = 1e-3;

				realnum frac = abs(oldPdest[ind]-(*tr).Emis().Pdest())/max((*tr).Emis().Pdest(),abslim);
				if (frac > bigerror)
				{
					bigerror = frac;
					trworst = tr;
					iworst = ind;
					cworst = 1;
				}
				frac = abs(oldPesc[ind]-(*tr).Emis().Pesc())/max((*tr).Emis().Pesc(),abslim);
				if (frac > bigerror)
				{
					bigerror = frac;
					trworst = tr;
					iworst = ind;
					cworst = 2;
				}
				frac = abs(oldPelec_esc[ind]-(*tr).Emis().Pelec_esc())/max((*tr).Emis().Pelec_esc(),abslim);
				if (frac > bigerror)
				{
					bigerror = frac;
					trworst = tr;
					iworst = ind;
					cworst = 3;
				}
				++ind;
			}		 
		}
		if (0 && bigerror > 0.0)
		{
			fprintf(ioQQQ,"RTEscChange nz %4ld %6f ",nzone,bigerror);
			fprintf(ioQQQ,"dst%c %6f=>%6f esc%c %6f=>%6f eesc%c %6f=>%6f -- %s\n",
					  cworst == 1 ? '*':' ',oldPdest[iworst],(*trworst).Emis().Pdest(),
					  cworst == 2 ? '*':' ',oldPesc[iworst],(*trworst).Emis().Pesc(),
					  cworst == 3 ? '*':' ',oldPelec_esc[iworst],(*trworst).Emis().Pelec_esc(),
					  (*trworst).chLabel().c_str());
		}
		*error = bigerror;
	}
}
	
void RT_line_all( linefunc line_one )
{
	long int ion,
		ipISO,
		nelem;
	long ipHi , ipLo;

	DEBUG_ENTRY( "RT_line_all()" );

	if( trace.lgTrace )
		fprintf( ioQQQ, "     RT_line_all called\n" );

	/* rfield.lgDoLineTrans - skip line transfer if requested with no line transfer
	 * conv.nPres2Ioniz is zero during search phase and on first call in this zone */
	if( !rfield.lgDoLineTrans && conv.nPres2Ioniz )
	{
		return;
	}

	/*if( lgUpdateFineOpac )
		fprintf(ioQQQ,"debuggg\tlgUpdateFineOpac in rt_line_all\n");*/
	
	static vector<realnum> DopplerWidth(LIMELM);
	for( nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		DopplerWidth[nelem] = GetDopplerWidth(dense.AtomicWeight[nelem]);	
	}

	for( ipISO=ipH_LIKE; ipISO < NISO; ++ipISO )
	{
		/* loop over all iso-electronic sequences */
		for( nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* parent ion stage, for H is 1, for He is 1 for He-like and 
			 * 2 for H-like */
			ion = nelem+1-ipISO;

			/* element turned off */
			if( !dense.lgElmtOn[nelem] )
				continue;
			/* need we consider this ion? */
			if( ion <= dense.IonHigh[nelem] )
			{
				/* loop over all lines */
				for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_local; ++ipHi )
				{
					for( ipLo=0; ipLo < ipHi; ++ipLo )
					{
						/* negative ipCont means this is not a real line, so do not
						 * transfer it */
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() < 0 ) 
							continue;

						/* generate escape prob, pumping rate, destruction prob, 
						 * inward outward fracs */
						fixit("should this use pestrk_up or pestrk?");
						line_one( iso_sp[ipISO][nelem].trans(ipHi,ipLo),
							     true,(realnum)iso_sp[ipISO][nelem].ex[ipHi][ipLo].pestrk_up,
									 DopplerWidth[nelem]);

						/* set true to print pump rates*/
						enum {DEBUG_LOC=false};
						if( DEBUG_LOC && nelem==1&& ipLo==0  /*&& iteration==2*/ )
						{
							fprintf(ioQQQ,"DEBUG pdest\t%3li\t%.2f\t%.3e\n",
								ipHi ,
								fnzone,
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Pdest());
						}
					}
				}
			}
		}
	}

	if (conv.lgFirstSweepThisZone || conv.lgLastSweepThisZone ) 
	{
		for( ipISO=ipH_LIKE; ipISO < NISO; ++ipISO )
		{
			/* loop over all iso-electronic sequences */
			for( nelem=ipISO; nelem < LIMELM; ++nelem )
			{
				/* parent ion stage, for H is 1, for He is 1 for He-like and 
				 * 2 for H-like */
				ion = nelem+1-ipISO;
				
				/* element turned off */
				if( !dense.lgElmtOn[nelem] )
					continue;
				/* need we consider this ion? */
				if( ion <= dense.IonHigh[nelem] )
				{
					/* loop over all lines */
					ipLo = 0;
					/* these are the extra Lyman lines for the iso sequences */
					/*for( ipHi=2; ipHi < iso_ctrl.nLyman[ipISO]; ipHi++ )*/
					/* only update if significant abundance and need to update fine opac */
					if( dense.xIonDense[nelem][ion] > 1e-30 )
					{
						for( ipHi=iso_sp[ipISO][nelem].st[iso_sp[ipISO][nelem].numLevels_local-1].n()+1; ipHi < iso_ctrl.nLyman[ipISO]; ipHi++ )
						{
							TransitionList::iterator tr = ExtraLymanLines[ipISO][nelem].begin()+ipExtraLymanLines[ipISO][nelem][ipHi];
							/* we just want the population of the ground state */
							(*tr).Emis().PopOpc() = iso_sp[ipISO][nelem].st[0].Pop();
							(*(*tr).Lo()).Pop() =
								iso_sp[ipISO][nelem].st[ipLo].Pop();
							
							/* actually do the work */
							line_one( *tr, true, 0.f, DopplerWidth[nelem]);
						}
					}					
				}/* if nelem if ion <=dense.IonHigh */
			}/* loop over nelem */
		}/* loop over ipISO */

		/* this is a major time sink for this routine - only evaluate on last
		 * sweep when fine opacities are updated since only effect of UTAs is
		 * to pump inner shell lines and add to total opacity */

		for( size_t i=0; i < UTALines.size(); i++ )
		{
			/* these are not defined in cooling routines so must be set here */
			UTALines[i].Emis().PopOpc() = dense.xIonDense[(*UTALines[i].Hi()).nelem()-1][(*UTALines[i].Hi()).IonStg()-1];
			(*UTALines[i].Lo()).Pop() = dense.xIonDense[(*UTALines[i].Hi()).nelem()-1][(*UTALines[i].Hi()).IonStg()-1];
			(*UTALines[i].Hi()).Pop() = 0.;
			line_one( UTALines[i], true,0.f, 
						 DopplerWidth[(*UTALines[i].Hi()).nelem()-1] );
		}
	}


	/* external database lines */
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		if( dBaseSpecies[ipSpecies].lgActive )
		{
			realnum DopplerWidth = GetDopplerWidth( dBaseSpecies[ipSpecies].fmolweight );
			for (TransitionList::iterator tr=dBaseTrans[ipSpecies].begin(); 
				  tr != dBaseTrans[ipSpecies].end(); ++tr)
			{	
				int ipHi = (*tr).ipHi();
				if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
					continue;
				if( (*tr).ipCont() > 0 )
				{
					line_one( *tr, true, 0.f, DopplerWidth );
				}
			}
		}
	}

	/* do satellite lines for iso sequences gt H-like
	 * H-like ions only have one electron, no satellite lines. */
	for( ipISO=ipHE_LIKE; ipISO<NISO; ++ipISO )
	{
		/* loop over all iso-electronic sequences */
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] && iso_ctrl.lgDielRecom[ipISO] )
			{
				for( ipLo=0; ipLo < iso_sp[ipISO][nelem].numLevels_local; ++ipLo )
				{
					line_one( SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][ipLo]], true,0.f, 
						 DopplerWidth[nelem] );
				}
			}
		}
	}

	/* the large H2 molecule */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_RTMake( line_one );

	for( size_t i=0; i < HFLines.size(); i++ )
	{
		line_one( HFLines[i], true,0.f, DopplerWidth[(*HFLines[i].Hi()).nelem()-1] );
	}

	for( long i=0; i < nWindLine; i++ )
	{
		ion = (*TauLine2[i].Hi()).IonStg();
		nelem = (*TauLine2[i].Hi()).nelem();

		if( (dense.lgIonChiantiOn[nelem-1][ion-1] && !atmdat.lgChiantiLvl2Hybrid) ||
				(dense.lgIonStoutOn[nelem-1][ion-1] && !atmdat.lgStoutLvl2Hybrid) )
		{
			/* If a species uses Chianti or Stout and hybrid is off, skip the level 2 lines */
			continue;
		}
		// iso sequence ions are done elsewhere
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO && 
			// dense.maxWN is positive if hybrid chianti is turned on and this element is included
			// in CloudyChianti.ini - zero otherwise
			 TauLine2[i].EnergyWN() >= dense.maxWN[nelem-1][ion-1])
		{
			line_one( TauLine2[i], true,0.f, DopplerWidth[(*TauLine2[i].Hi()).nelem()-1] );
		}
	}

	return;
}

void RT_fine_clear()
{
	DEBUG_ENTRY( "RT_fine_init()" );
	/* zero out fine opacity array */
	/* this array is huge and takes significant time to zero out or update, 
	 * only do so when needed, */
	vzero(rfield.fine_opac_zone);
	
	if( 0 && cosmology.lgDo )
	{
		realnum dVel = (realnum)SPEEDLIGHT * ( 1.f - 1.f/POW2(1.f+cosmology.redshift_start-cosmology.redshift_current) );
		rfield.ipFineConVelShift = -(long int)( dVel/ rfield.fine_opac_velocity_width + 0.5 );
	}
	else
	{
		/* this is offset within fine opacity array caused by Doppler shift
		 * between first zone, the standard of rest, and current position 
		 * in case of acceleration away from star in outflowing wind 
		 * dWind is positive, this means that the rest frame original
		 * velocity is red shifted to lower energy */
		realnum dWind = wind.windv - wind.windv0;
		
		/* will add ipVelShift to index of original energy so redshift has 
		 * to be negative */
		rfield.ipFineConVelShift = -(long int)(dWind / rfield.fine_opac_velocity_width+0.5);
	}
}
