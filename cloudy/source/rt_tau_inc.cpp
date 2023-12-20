/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_tau_inc increment optical depths once per zone, called after radius_increment */
#include "cddefines.h"
#include "elementnames.h"
#include "taulines.h"
#include "iso.h"
#include "rfield.h"
#include "trace.h"
#include "dense.h"
#include "prt.h"
#include "conv.h"
#include "h2.h"
#include "hmi.h"
#include "opacity.h"
#include "cooling.h"
#include "thermal.h"
#include "radius.h"
#include "rt.h"
#include "doppvel.h"
#include "mole.h"


void prt_trans_opc_debug( const char *LineGroup, const TransitionProxy &t )
{
	fprintf( ioQQQ,
		"%s:\t label= '%s'\t nelem= %d ('%s')\t ion= %2d\t dense= %.4e\t TauCon = %.4e\t TauIn= %.4e\t TauTot= %.4e\n",
		LineGroup,
		t.chLabel().c_str(),
		int((*t.Hi()).nelem()),
		elementnames.chElementSym[(*t.Hi()).nelem()-1],
		int((*t.Hi()).IonStg()),
		dense.xIonDense[(*t.Hi()).nelem()-1][(*t.Hi()).IonStg()-1],
		t.Emis().TauCon(),
		t.Emis().TauIn(),
		t.Emis().TauTot() );
}

/*RT_tau_inc increment optical depths once per zone, called after radius_increment */
void RT_tau_inc(void)
{
	DEBUG_ENTRY( "RT_tau_inc()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " RT_tau_inc called.\n" );
	}

	/* call RT_line_all one last time in this zone, to get fine opacities defined */
	ASSERT( !conv.lgFirstSweepThisZone );
	conv.lgLastSweepThisZone = true;
	RT_fine_clear();
	RT_line_all( RT_line_one_fine );

	/* rfield.lgOpacityFine flag set false with no fine opacities command 
	 * tests show that always evaluating this changes fast run of
	 * pn_paris from 26.7 sec to 35.1 sec 
	 * but must always update fine opacities since used for transmission */
	if( rfield.lgOpacityFine )
	{
		/* increment the fine opacity array */
		for( long i=0; i<rfield.nfine; ++i )
		{
			realnum tauzone = rfield.fine_opac_zone[i]*(realnum)radius.drad_x_fillfac;
			rfield.fine_opt_depth[i] += tauzone;
		}
		rfield.trans_coef_total_stale = true;
	}

	/* this may have updated some escape/destruction rates - force update
	 * to all cooling lines */
	CoolEvaluate( &thermal.ctot );

	if( nzone <=1 )
	{
		opac.telec = (realnum)(radius.drad_x_fillfac*dense.eden*6.65e-25);
		opac.thmin = (realnum)(radius.drad_x_fillfac*findspecieslocal("H-")->den*3.9e-17*
			(1. - rfield.ContBoltz[hmi.iphmin-1]/ hmi.hmidep));
	}
	else
	{
		opac.telec += (realnum)(radius.drad_x_fillfac*dense.eden*6.65e-25);
		opac.thmin += (realnum)(radius.drad_x_fillfac*findspecieslocal("H-")->den*3.9e-17*
			(1. - rfield.ContBoltz[hmi.iphmin-1]/ hmi.hmidep));
	}

	/* prevent maser runaway */
	rt.dTauMase = 0;
	rt.mas_species = 0;
	rt.mas_ion = 0;
	rt.mas_hi = 0;
	rt.mas_lo = 0;

	static vector<realnum> DopplerWidth(LIMELM);
	for (long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem)
	{
		DopplerWidth[nelem] = GetDopplerWidth(dense.AtomicWeight[nelem]);
	}

	/* all lines in iso sequences */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( ! dense.lgElmtOn[ nelem ] )
				continue;

			/* this is the parent ion, for HI lines, is 1, 
			 * for element He is 1 for He-like (HeI) and 2 for H-like (HeII) */
			int ion = nelem+1-ipISO;
			/* do not evaluate in case where trivial parent ion */
			if( ion <=dense.IonHigh[nelem] && dense.xIonDense[nelem][ion] > dense.density_low_limit )
			{
				if( iso_ctrl.lgDielRecom[ipISO] )
				{
					// SatelliteLines are indexed by lower level
					for( long ipLo=0; ipLo < iso_sp[ipISO][nelem].numLevels_local; ipLo++ )
					{
						RT_line_one_tauinc(SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][ipLo]], ipISO, nelem, -1, ipLo, 
							 DopplerWidth[nelem] );
					}
				}

				for( long ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
				{
					for( long ipLo=0; ipLo < ipHi; ipLo++ )
					{
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() <= 0 )
							continue;

						/* actually do the work */
						RT_line_one_tauinc(iso_sp[ipISO][nelem].trans(ipHi,ipLo), ipISO, nelem, ipHi, ipLo, 
							DopplerWidth[nelem] );
					}
				}
				/* these are the extra Lyman lines, use all lines so
				 * totals are correct as attribution may change */
				for( long ipHi=2; ipHi < iso_ctrl.nLyman[ipISO]; ipHi++ )
				{
					TransitionList::iterator tr = ExtraLymanLines[ipISO][nelem].begin()+ipExtraLymanLines[ipISO][nelem][ipHi];
					(*tr).Emis().PopOpc() = iso_sp[ipISO][nelem].st[0].Pop();

					/* actually do the work */
					RT_line_one_tauinc(*tr, -1 ,ipISO, nelem, ipHi,
						DopplerWidth[nelem] );
				}
			}
		}
	}

	/* increment optical depths for all heavy element lines
	 * same routine does wind and static */
	/* all lines in cooling with g-bar */
	for( long i=0; i < nWindLine; i++ )
	{
		if( ! dense.lgElmtOn[ (*TauLine2[i].Hi()).nelem()-1 ] )
			continue;

		/* do not include H-like or He-like in the level two lines since
		 * these are already counted in iso sequences */
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
		{
			RT_line_one_tauinc(TauLine2[i], -3, -3, -3, i, DopplerWidth[(*TauLine2[i].Hi()).nelem()-1] );
		}
	}

	/* the block of inner shell lines */
	for( size_t i=0; i < UTALines.size(); i++ )
	{
		if( ! dense.lgElmtOn[ (*UTALines[i].Hi()).nelem()-1 ] )
			continue;

		/* populations have not been set */
		UTALines[i].Emis().PopOpc() = dense.xIonDense[(*UTALines[i].Hi()).nelem()-1][(*UTALines[i].Hi()).IonStg()-1];
		(*UTALines[i].Lo()).Pop() = dense.xIonDense[(*UTALines[i].Hi()).nelem()-1][(*UTALines[i].Hi()).IonStg()-1];
		(*UTALines[i].Hi()).Pop() = 0.;
		RT_line_one_tauinc(UTALines[i], -4 , -4 , -4 , i, DopplerWidth[(*UTALines[i].Hi()).nelem()-1] );

		//	prt_trans_opc_debug( "UTA", UTALines[i] );
	}

	/* all hyper fine structure lines  */
	for( size_t i=0; i < HFLines.size(); i++ )
	{
		if( ! dense.lgElmtOn[ (*HFLines[i].Hi()).nelem()-1 ] )
			continue;

		RT_line_one_tauinc(HFLines[i] , -5 , -5 , -5 , i, DopplerWidth[(*HFLines[i].Hi()).nelem()-1] );
	}

	/* increment optical depth for the H2 molecule */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_RT_tau_inc();

	/* database Lines*/
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		if( dBaseSpecies[ipSpecies].lgActive )
		{
			realnum DopplerWidth = GetDopplerWidth( dBaseSpecies[ipSpecies].fmolweight );
			for (TransitionList::iterator tr=dBaseTrans[ipSpecies].begin(); 
				  tr != dBaseTrans[ipSpecies].end(); ++tr)
			{	
				int ipHi = (*tr).ipHi();
				if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local || (*tr).ipCont() <= 0)
					continue;
				int ipLo = (*tr).ipLo();

				RT_line_one_tauinc( *tr, -10, ipSpecies, ipHi, ipLo, DopplerWidth );
			}
		}
	}
	
	if( trace.lgTrace && trace.lgOptcBug )
	{
		fprintf( ioQQQ, " RT_tau_inc updated optical depths:\n" );
		prtmet();
	}

	if( trace.lgTrace )
		fprintf( ioQQQ, " RT_tau_inc returns.\n" );

	return;
}
