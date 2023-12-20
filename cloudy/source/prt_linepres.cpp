/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtLinePres print line radiation pressures for current conditions */
#include "cddefines.h"
#include "pressure.h"
#include "taulines.h"
#include "iso.h"
#include "dense.h"
#include "h2.h"
#include "prt.h"
#include "doppvel.h"

/* faintest pressure we will bother with */
static const double THRESH = 0.05;

void PrtLinePres( FILE *ioPRESSURE )
{
	long int ip,
	  ipLo,
	  ipHi,
	  nelem;
	int ier;
	double RadPres1;

	// this is limit to number of lines we can store at one time
	const int nLine	= 100;
	/* labels, wavelengths, and fraction of total pressure, for important lines */
	string chLab[nLine];
	realnum wl[nLine] , frac[nLine];
	long int iperm[nLine];

	/* will be used to check on size of opacity, was capped at this value */
	realnum smallfloat=SMALLFLOAT*10.f;

	DEBUG_ENTRY( "PrtLinePres()" );

	/* this routine only called if printout of contributors to line
	 * radiation pressure is desired */

	/* compute radiation pressure due to iso lines */
	ip = 0;

	for(long i=0; i<nLine; ++i)
	{
		frac[i] = FLT_MAX;
		wl[i] = FLT_MAX;
	}

	if( pressure.pres_radiation_lines_curr > 1e-30 )
	{
		/** \todo	1	make this and eval rad pressure same routine, with flag saying to 
		 * print contributors - copy code from other routine - this code has been
		 * left behind */
		for( long ipISO = 0; ipISO<NISO; ipISO++ )
		{
			for( nelem=ipISO; nelem < LIMELM; nelem++ )
			{
				/* does this ion stage exist? */
				if( dense.IonHigh[nelem] >= nelem + 1 - ipISO )
				{
					/* do not include highest levels since maser can occur due to topoff,
					 * and pops are set to small number in this case */
					for( ipHi=1; ipHi <iso_sp[ipISO][nelem].numLevels_local - iso_sp[ipISO][nelem].nCollapsed_local; ipHi++ )
					{
						for( ipLo=0; ipLo < ipHi; ipLo++ ) 
						{
							if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() <= 0 )
								continue;

							ASSERT( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() > iso_ctrl.SmallA );

							/*NB - this code must be kept parallel with code in pressure_total */
							if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().PopOpc() > smallfloat &&
								/* >>chng 01 nov 01, add test that have not overrun optical depth scale */
								( (iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauTot() - iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauIn()) > smallfloat ) )
							{
								RadPres1 =  PressureRadiationLine( iso_sp[ipISO][nelem].trans(ipHi,ipLo), GetDopplerWidth(dense.AtomicWeight[nelem]) );

								if( RadPres1 > pressure.pres_radiation_lines_curr*THRESH )
								{
									wl[ip] = iso_sp[ipISO][nelem].trans(ipHi,ipLo).WLAng();

									/* put null terminated line label into chLab using line array*/
									chLab[ip] = chIonLbl(iso_sp[ipISO][nelem].trans(ipHi,ipLo));
									frac[ip] = (realnum)(RadPres1/pressure.pres_radiation_lines_curr);

									ip = MIN2((long)nLine-1,ip+1);
								}
							}
						}
					}
				}
			}
		}

		for( long i=0; i < nWindLine; i++ )
		{
			if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
			{
				if( (*TauLine2[i].Hi()).Pop() > 1e-30 )
				{
					RadPres1 =  PressureRadiationLine( TauLine2[i], GetDopplerWidth(dense.AtomicWeight[(*TauLine2[i].Hi()).nelem()-1]) );
				}
				else
				{
					RadPres1 = 0.;
				}

				if( RadPres1 > pressure.pres_radiation_lines_curr*THRESH )
				{
					wl[ip] = TauLine2[i].WLAng();

					/* put null terminated line label into chLab using line array*/
					chLab[ip] = chIonLbl(TauLine2[i]);
					frac[ip] = (realnum)(RadPres1/pressure.pres_radiation_lines_curr);

					ip = MIN2((long)nLine-1,ip+1);
				}
			}
		}

		for( size_t i=0; i < HFLines.size(); i++ )
		{
			if( (*HFLines[i].Hi()).Pop() > 1e-30 )
			{
				RadPres1 =  PressureRadiationLine( HFLines[i], GetDopplerWidth(dense.AtomicWeight[(*HFLines[i].Hi()).nelem()-1]) ); 
			}
			else
			{
				RadPres1 = 0.;
			}

			if( RadPres1 > pressure.pres_radiation_lines_curr*THRESH )
			{
				wl[ip] = HFLines[i].WLAng();

				/* put null terminated line label into chLab using line array*/
				chLab[ip] = chIonLbl(HFLines[i]);
				frac[ip] = (realnum)(RadPres1/pressure.pres_radiation_lines_curr);

				ip = MIN2((long)nLine-1,ip+1);
			}
		}

		/* lines from external databases */
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
					if( (*(*tr).Hi()).Pop() > 1e-30 )
					{
						RadPres1 =  PressureRadiationLine( *tr, DopplerWidth );
					}
					else
					{
						RadPres1 = 0.;
					}
					
					if( RadPres1 > pressure.pres_radiation_lines_curr*THRESH )
					{
						wl[ip] = (*tr).WLAng();
						
						/* put null terminated line label into chLab using line array*/
						chLab[ip] = chIonLbl(*tr );
						frac[ip] = (realnum)(RadPres1/pressure.pres_radiation_lines_curr);

						ip = MIN2((long)nLine-1,ip+1);
					}
				}       
			}       
		}       

		/* radiation pressure due to H2 */
		for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		{
			RadPres1 = (*diatom)->H2_RadPress();
			if( RadPres1 > pressure.pres_radiation_lines_curr*THRESH )
			{
				wl[ip] = 0;

				/* put null terminated 4 char line label into chLab using line array*/
				chLab[ip] = "H2  ";
				frac[ip] = (realnum)(RadPres1/pressure.pres_radiation_lines_curr);

				ip = MIN2((long)nLine-1,ip+1);
			}
		}

		/* return if no significant contributors to radiation pressure found */
		if( ip<= 0 )
		{
			fprintf( ioPRESSURE, "\n" );
			return;
		}

		/* sort from strong to weak, then only print strongest */
		spsort(
			/* input array to be sorted */
			frac, 
			/* number of values in x */
			ip, 
			/* permutation output array */
			iperm, 
			/* flag saying what to do - 1 sorts into increasing order, not changing
			* the original routine */
			-1, 
			/* error condition, should be 0 */
			&ier);

		/* now print up to ten strongest contributors to radiation pressure */
		fprintf( ioPRESSURE, " P(Lines):" );
		for( long i=0; i < MIN2(10,ip); i++ )
		{
			int ipline = iperm[i];
			fprintf( ioPRESSURE, "(%4.4s ", chLab[ipline].c_str());
			prt_wl(ioPRESSURE, wl[ipline] );
			fprintf(ioPRESSURE," %.2f) ",frac[ipline]);
		}

		/* finally end the line */
		fprintf( ioPRESSURE, "\n" );
	}
	return;
}
