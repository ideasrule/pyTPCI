/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_hydro put H-like iso sequence into line intensity stack */
#include "cddefines.h"
#include "atmdat.h"
#include "dense.h"
#include "prt.h"
#include "hydrogenic.h"
#include "iso.h"
#include "rfield.h"
#include "geometry.h"
#include "lines.h"
#include "phycon.h"
#include "radius.h"
#include "secondaries.h"
#include "trace.h"
#include "two_photon.h"
#include "lines_service.h"
#include "elementnames.h"
#include "ipoint.h"

static const int nCharL = 21;
static const string chL[nCharL] = {"s","p","d","f","g","h","i","k","l","m","n","o","q","r","t","u","v","w","x","y","z"};

inline string l2str(long l, bool lgShort = false)
{
	if( l < nCharL )
		return chL[l];
	else
	{
		ostringstream oss;
		if( lgShort )
			oss << l;
		else
			oss << "[" << l << "]";
		return oss.str();
	}
}

string iso_comment_tran_levels( long ipISO, long nelem, long ipLo, long ipHi )
{
	string isoSeq = ( ipISO == ipHE_LIKE ) ? "He-like, " : "H-like, ";
	return isoSeq + GenerateTransitionConfiguration( iso_sp[ipISO][nelem].trans(ipHi,ipLo) );
}

void lines_hydro(void)
{
	long ipISO = ipH_LIKE;
	long int i, nelem, ipHi, ipLo;
	string chLabel="    ";

	double hbetab, 
		em , 
		caseb;

	DEBUG_ENTRY( "lines_hydro()" );

	if( trace.lgTrace )
		fprintf( ioQQQ, "   lines_hydro called\n" );

	// this can be changed with the atom levels command but must be at least 3
	ASSERT( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max >= 3 );
	ASSERT( !dense.lgElmtOn[ipHELIUM] || iso_sp[ipH_LIKE][ipHELIUM].n_HighestResolved_max >= 3 );

	i = StuffComment( "H-like iso-sequence" );
	linadd( 0., (realnum)i , "####", 'i',
		" start H -like iso sequence ");

	/*fprintf(ioQQQ," debugg\t%.2e\t%.2e\t%.2e\n", 
		radius.drad,
		iso_sp[ipH_LIKE][ipHYDROGEN].xLineTotCool , 
		iso_sp[ipH_LIKE][ipHYDROGEN].cLya_cool);*/

	/* >>chng 95 jun 25 changed from info to cooling to pick this up in primal.in   */
	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].cLya_cool),1215.67,"Cool",'i',
		"collisionally excited La cooling ");

	linadd(MAX2(0.,-iso_sp[ipH_LIKE][ipHYDROGEN].cLya_cool),1215.67,"Heat",'i',
		"  collisionally de-excited La heating ");

	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].cLyrest_cool),960,"Crst",'i',
		"  cooling due to n>2 Lyman lines ");

	linadd(MAX2(0.,-iso_sp[ipH_LIKE][ipHYDROGEN].cLyrest_cool),960,"Hrst",'i',
		"  heating due to n>2 Lyman lines ");

	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].cBal_cool),Hbeta_WavLen,"Crst",'i',
		"  cooling due to n>3 Balmer lines ");

	linadd(MAX2(0.,-iso_sp[ipH_LIKE][ipHYDROGEN].cBal_cool),Hbeta_WavLen,"Hrst",'i',
		"  heating due to n>3 Balmer lines ");

	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].cRest_cool),0,"Crst",'i',
		"  cooling due to higher Paschen lines ");

	linadd(MAX2(0.,-iso_sp[ipH_LIKE][ipHYDROGEN].cRest_cool),0,"Hrst",'i',
		"  heating due to higher Paschen lines ");

	/* remember largest fractional ionization of H due to secondaries */
	secondaries.SecHIonMax = MAX2( secondaries.SecHIonMax , secondaries.sec2total );

	/* remember fraction of H ionizations due to ct */
	atmdat.HIonFracMax = MAX2( atmdat.HIonFracMax, atmdat.HIonFrac);

	/* remember largest fraction of thermal collisional ionization of H ground state */
	hydro.HCollIonMax = 
		(realnum)MAX2( hydro.HCollIonMax , hydro.H_ion_frac_collis );

	linadd(secondaries.x12tot*iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop()*1.634e-11,1215.67,"LA X" ,'i',
		"Lya contribution from suprathermal secondaries from ground ");

	linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHYDROGEN].coll_ion),0,"CION",'c',
		"collision ionization cooling of hydrogen ");

	linadd(MAX2(-iso_sp[ipH_LIKE][ipHYDROGEN].coll_ion,0.),0,"3bHt",'h',
		"  this is the heating due to 3-body recombination ");

	if( dense.lgElmtOn[ipHELIUM] )
	{
		linadd(MAX2(0.,iso_sp[ipH_LIKE][ipHELIUM].coll_ion),0,"He2C",'c',
			   "collision ionization cooling of He+ ");

		linadd(MAX2(-iso_sp[ipH_LIKE][ipHELIUM].coll_ion,0.),0,"He2H",'h',
			   "  this is the heating due to 3-body recombination onto He+");
	}

	fixit("why is there a zero here?");
	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop()*0.*iso_sp[ipH_LIKE][ipHYDROGEN].ex[ipH2p][ipH1s].pestrk*1.634e-11,1215.67,"Strk",'i',
	  "  Stark broadening contribution to line ");

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH3s].Pop()*iso_sp[ipH_LIKE][ipHYDROGEN].ex[ipH3s][ipH2p].pestrk*3.025e-12,
	  6562.80,"Strk",'i',
	  "  Stark broadening contribution to line ");

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH4s].Pop()*iso_sp[ipH_LIKE][ipHYDROGEN].ex[ipH4s][ipH2p].pestrk*4.084e-12,
	  Hbeta_WavLen,"Strk",'i',
	  "Stark broadening contribution to line ");

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH4p].Pop()*iso_sp[ipH_LIKE][ipHYDROGEN].ex[ipH4p][ipH3s].pestrk*1.059e-12,
	  18751,"Strk",'i',
		   " Stark broadening contribution to line ");

	/* pestrk[5,4] is A[4,5]*pest[4,5] 
	 * Stark broadening contribution to line */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max >= 5 )
	{
		long ip5p = iso_sp[ipH_LIKE][ipHYDROGEN].QN2Index(5, 1, 2);
		linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ip5p].Pop()*iso_sp[ipH_LIKE][ipHYDROGEN].ex[ip5p][ipH4s].pestrk*4.900e-13,40512,"Strk",'i',
			"Stark broadening part of line");
	}
	/* this can fail if RT_line_all never updates the ots rates, a logic error,
	 * but only assert this during actual calculation (ipass>0), */
	ASSERT( LineSave.ipass  <1 ||
		iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().ots()>= 0.);

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().ots()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).EnergyErg(), 1215.67,"Dest",'i',
		"  portion of line lost due to absorp by background opacity ");

	/* portion of line lost due to absorb by background opacity */
	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH3p,ipH2s).Emis().ots()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH3p,ipH2s).EnergyErg(), 6562.80,"Dest",'i',
		"Ha destroyed by background opacity");

	/* portion of line lost due to absorp by background opacity */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max >= 5 )
	{
		long ip5p = iso_sp[ipH_LIKE][ipHYDROGEN].QN2Index(5, 1, 2);
		linadd(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ip5p,ipH4s).Emis().ots()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ip5p,ipH4s).EnergyErg(),40516, "Dest",'i',
			"portion of line lost due to absorb by background opacity");
	}

	/* portion of line lost due to absorb by background opacity */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_max > ipH4p )
		linadd(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH4p,ipH2s).Emis().ots()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH4p,ipH2s).EnergyErg(), Hbeta_WavLen, "Dest",'i',
			"portion of line lost due to absorb by background opacity");

	/* portion of line lost due to absorb by background opacity */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_max > ipH4p )
		linadd(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH4p,ipH3s).Emis().ots()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH4p,ipH3s).EnergyErg() ,18751, "Dest",'i',
			"portion of line lost due to absorb by background opacity");

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop()*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().Aul()*
		hydro.dstfe2lya*iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).EnergyErg() , 1215.67 , "Fe 2" , 'i',
		"Ly-alpha destroyed by overlap with FeII " );

	linadd(iso_sp[ipH_LIKE][ipHYDROGEN].RadRec_caseB*dense.xIonDense[ipHYDROGEN][1]*dense.eden * 1.64e-11,1215.67,"Ca B",'i',
		" simple high-density case b intensity of Ly-alpha, no two photon ");

	/* these entries only work correctly if the APERTURE command is not in effect */
	if( geometry.iEmissPower == 2 )
	{
		/* H-beta computed from Q(H) and specified covering factor */
		if( nzone == 1 )
		{
			/* evaluate the case b emissivity by interpolating on the hummer & storey tables */
			caseb = rfield.qhtot*
				atmdat_HS_caseB( 4 , 2 , 1 , phycon.te , dense.eden, 'b' ) / iso_sp[ipH_LIKE][ipHYDROGEN].RadRec_caseB;
			/* the atmdat_HS_caseB returned -1 if the physical conditions were outside range of validity.  
			 * In this case use simple approximation with no temperature or density dependence */
			if( caseb < 0 )
			{
				caseb = rfield.qhtot*4.75e-13;
			}
			LineSave.lines[LineSave.nsum].SumLineZero();
		}
		else
		{
			caseb = 0.;
		}
		/* H-beta computed from Q(H) and specified covering factor */
		linadd( caseb/radius.dVeffAper*geometry.covgeo , Hbeta_WavLen, "Q(H)" , 'i' ,
			"Case B H-beta computed from Q(H) and specified covering factor");

		if( nzone == 1 )
		{
			// the cast to double prevents an FPE with Solaris Studio 12.4 (limit_compton_hi_t.in)
			caseb = rfield.qhtot*double(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).EnergyErg());
			LineSave.lines[LineSave.nsum].SumLineZero();
		}
		else
		{
			caseb = 0.;
		}
		/* >>chng 02 nov 05, better approximation for Lya for temperature of first zone */
		linadd( caseb/radius.dVeffAper*geometry.covgeo , 1215.67 , "Q(H)" , 'i',
			"Ly-alpha from Q(H), high-dens lim, specified covering factor" );
	}

	/* this is the main printout, where line intensities are entered into the stack */
	for( nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			ASSERT( iso_sp[ipH_LIKE][nelem].n_HighestResolved_max >= 3 );

			for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
			{
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{
					if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					set_xIntensity( iso_sp[ipISO][nelem].trans(ipHi,ipLo) );
				}
			}
		}
	}

	/* option to save the fine structure components of the H-like lines */
	if( iso_ctrl.lgKeepFS )
	{
		for( nelem=0; nelem < LIMELM; nelem++ )
		{
			// next 4 loops make sure that fs components of the same line appear together...
			for( long nLo=1; nLo < iso_sp[ipH_LIKE][nelem].n_HighestResolved_max; nLo++ )
			{
				for( long nHi=nLo+1; nHi <= iso_sp[ipH_LIKE][nelem].n_HighestResolved_max; nHi++ )
				{
					for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
					{
						for( ipLo=0; ipLo < ipHi; ipLo++ )
						{
							// skip the 2s -> 1s M1 line
							if( N_(ipLo) != nLo || N_(ipHi) != nHi || abs(L_(ipHi)-L_(ipLo)) != 1 )
								continue;

							ostringstream oss1, oss2;
							if( LineSave.ipass == 0 )
							{
								string fsLabel = chIonLbl(nelem+1, nelem+1);
								oss1 << fsLabel << " " << l2str(L_(ipHi), true);
								if( L_(ipLo) == L_(ipHi)+1 )
									oss1 << "+";
								else
									oss1 << "-";
								oss2 << fsLabel << " fine structure component " << N_(ipHi);
								oss2 << l2str(L_(ipHi)) << " -> " << N_(ipLo) << l2str(L_(ipLo));
							}
							lindst(iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo), ExtraInten(),
								   oss1.str().c_str(), 't', false, oss2.str().c_str());
						}
					}
				}
			}
		}
	}

	/* create emissivity or intensity for hydrogenic species,
	 * first combine/bring balmer series together */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.IonHigh[nelem] == nelem + 1 )
		{
			/* bring nL - n'L' emission together as n-n' emission. */
			for( ipHi=1; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
			{
				long index_of_nHi_P;

				/* is ipHi is collapsed level, index_of_nHi_P is ipHi */
				if( N_(ipHi) > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max )
					index_of_nHi_P = ipHi;
				else
					index_of_nHi_P = iso_sp[ipH_LIKE][nelem].QN2Index(N_(ipHi), 1, 2);

				/* only need to consider resolved lower level here */
				for( ipLo=0; ipLo < ipHi; ipLo++ )
				{
					/* jump out if ipLo is collapsed 
					 * NB this must be up to n_HighestResolved_local and not n_HighestResolved_max */
					if( N_(ipLo) > iso_sp[ipH_LIKE][nelem].n_HighestResolved_local || N_(ipLo) == N_(ipHi) )
						break;

					long index_of_nLo_S = iso_sp[ipH_LIKE][nelem].QN2Index(N_(ipLo), 0, 2);

					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					/* add everything into nP - n'S, skip if current indices are those levels. */
					if( ipHi == index_of_nHi_P && ipLo == index_of_nLo_S )
						continue;
					else
					{
						/* add resolved line to nP - n'S */
						iso_sp[ipH_LIKE][nelem].trans(index_of_nHi_P,index_of_nLo_S).Emis().xIntensity() +=
							iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().xIntensity();
						iso_sp[ipH_LIKE][nelem].trans(index_of_nHi_P,index_of_nLo_S).Emis().xObsIntensity() +=
							iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().xObsIntensity();
						/* zero out the resolved line */
						iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().xIntensity() = 0;
						iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().xObsIntensity() = 0;
					}
				}
			}
		}
	}

	/* H beta recombination, assuming old case B */
	hbetab = (double)((exp10(-20.89 - 0.10612*POW2(phycon.alogte - 4.4)))/phycon.te);
	/* need to pass this assert if CaBo is to have valid array indices for ipCont */
	/* 06 aug 28, from numLevels_max to _local. */
	/* 06 dec 21, change from numLevels_max to _local was mistake for this entire file.  Undo. */
	ASSERT( iso_sp[ipH_LIKE][ipHYDROGEN].numLevels_max > 4 );
	hbetab *= dense.xIonDense[ipHYDROGEN][1]*dense.eden;

	lindst(hbetab, -Hbeta_WavLen ,"CaBo",
		1 ,'i',false,
		" this is old case b based on Ferland (1980) PASP ");

	if( dense.lgElmtOn[ipHELIUM] )
	{
		/* need to pass this assert if CaBo is to have valid array indices for ipCont */
		/* 06 aug 28, from numLevels_max to _local. */
		/* 06 dec 21, change from numLevels_max to _local was mistake for this entire file.  Undo. */
		ASSERT( iso_sp[ipH_LIKE][ipHELIUM].numLevels_max > 4 );
		/* 1640 1640 1640 */
		em = 2.03e-20/(phycon.te70*phycon.te10*phycon.te03);
		em *= dense.xIonDense[ipHELIUM][2]*dense.eden;

		lindst(em,-1640,"CaBo",
			1,'i',false,
			" old prediction of He II 1640, Case B at low densities");

		/* hydrogenic helium */
		/* old prediction of He II 4686, case B */
		em = 2.52e-20/(pow(phycon.te,1.05881));
		em *= dense.xIonDense[ipHELIUM][2]*dense.eden;

		lindst(em,-4685.64,"CaBo",	1,'i',false,
			   " old prediction of He II 4686, Case B at low densities");
	}

	/* predict case b intensities of hydrogen lines */
	if( LineSave.ipass <= 0 )
	{
		for(nelem=0; nelem<HS_NZ; ++nelem )
		{
			atmdat.lgHCaseBOK[0][nelem] = true;
			atmdat.lgHCaseBOK[1][nelem] = true;
		}
	}
	/* this is the main printout, where line intensities are entered into the stack */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* HS_NZ is limit to charge of elements in HS predictions, now 8 == oxygen */
			/* but don't do the minor elements, Li, Be, B - these were not read in and so should not be
			 * printed - remove equivalent if statement in atmdat_readin.cpp to read them in */
			if( nelem < HS_NZ && (nelem<2 || nelem>4) )
			{
				int iCase;
				for( iCase=0; iCase<2; ++iCase )
				{
					char chAB[2]={'A','B'};
					char chLab[5]="Ca  ";

					/* adding iCase means start from n=1 for case A, n=2 for Case B,
					 * note that principal quantum number is on physics scale, not C */
					/* 06 aug 28, both of these from numLevels_max to _local. */
					/* 06 dec 21, change from numLevels_max to _local was mistake for this entire file.  Undo. */
					/* HS data file extends to ipHi == atmdat.ncut[iCase][nelem]==25, report all possible values for comparisons */
					for( ipLo=1+iCase; ipLo<MIN2(atmdat.ncut[iCase][nelem]-1,iso_sp[ipH_LIKE][nelem].n_HighestResolved_max + iso_sp[ipH_LIKE][nelem].nCollapsed_max); ++ipLo )
					{
						for( ipHi=ipLo+1; ipHi< MIN2(atmdat.ncut[iCase][nelem],iso_sp[ipH_LIKE][nelem].n_HighestResolved_max + iso_sp[ipH_LIKE][nelem].nCollapsed_max+1); ++ipHi )
						{
							/* Put case b predictions into line stack
							 * NB NB NB each Hummer & Storey case b line must be 
							 * explicitly clobbered by hand in routine final if 
							 * atmdat.lgHCaseBOK[iCase][nelem] flag is set false
							 * since this indicates that we exceeded bounds of table,
							 * DO NOT want to print lines in that case */

							/* first do case b emissivity of balmer lines */

							/* get HS predictions */
							double case_b_Intensity = atmdat_HS_caseB(ipHi, ipLo, nelem+1, phycon.te, dense.eden, chAB[iCase]);
							if( case_b_Intensity<=0. )
							{
								atmdat.lgHCaseBOK[iCase][nelem] = false;
								case_b_Intensity = 0.;
							}

							case_b_Intensity *= dense.xIonDense[nelem][nelem+1-ipISO]*dense.eden;

							/* make label either Ca A or Ca B */
							chLab[3] = chAB[iCase];

							/* this is freq in cm^-1 of interpolated case b from HS tables */
							/* NB NB
							 *
							 * 2020-08-16
							 *
							 * This is a temporary hack, until the data are resolved in the input
							 * data file.  See comment with the same date in iso_create.cpp.
							 */
							realnum Enerwn = realnum(hydro_energy(nelem, ipLo, -1, -1, -1) -
										 hydro_energy(nelem, ipHi, -1, -1, -1));
							realnum wl = (realnum)wn2ang( double(Enerwn) );
							atmdat.WaveLengthCaseB[nelem][ipHi][ipLo] = wl;
							long ip = ipoint( Enerwn*WAVNRYD );
							lindst(case_b_Intensity,wl,chLab,ip,'i',false," case a or case b from Hummer & Storey tables" );
						}
					}
				}
			}

			// add two-photon details here
			if( LineSave.ipass == 0 )
			{
				/* chIonLbl is function that generates a null terminated 4 char string, of form "C  2" 
				 * the result, chLable, is only used when ipass == 0, can be undefined otherwise */
				chLabel = chIonLbl(nelem+1, nelem+1-ipISO);
			}
			for( auto tnu = iso_sp[ipH_LIKE][nelem].TwoNu.begin(); tnu != iso_sp[ipH_LIKE][nelem].TwoNu.end(); ++tnu )
			{
				fixit("This was multiplied by Pesc when treated as a line, now what?  Only used for printout?");
				fixit("below should be 'i' instead of 'r' ?");

				string tpc_comment = "";
				if( LineSave.ipass == 0 )
				{
					tpc_comment = " two photon continuum, " +
						iso_comment_tran_levels( ipISO, nelem, (*tnu).ipLo, (*tnu).ipHi );
				}
				linadd(	tnu->AulTotal * tnu->E2nu * EN1RYD * (*tnu->Pop),
					2. * wn2ang( iso_sp[ipH_LIKE][nelem].trans( (*tnu).ipHi, (*tnu).ipLo ).EnergyWN() ),
					chLabel.c_str(), 'r', tpc_comment.c_str() );
			}

			/* NB NB - low and high must be in this order so that all balmer, paschen,
			 * etc series line up correctly in final printout */
			for( ipLo=ipH1s; ipLo < iso_sp[ipH_LIKE][nelem].numLevels_max-1; ipLo++ )
			{
				/* don't bother with decays to 2p since we set them to zero above */
				if( ipLo==ipH2p )
					continue;

				/* set number of levels to print */
				long int nLoop  = iso_Max_Emitting_Level(nelem, ipISO, prt.lgPrnIsoCollapsed);
				long index_of_nLo_S = iso_sp[ipH_LIKE][nelem].QN2Index(N_(ipLo), 0, 2);

				for( ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					// skip non-radiative transitions
					if( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).ipCont() < 1 )
						continue;

					// skip 2s-1s, so that 2p-1s comes first and cdLine finds LyA instead of the M1 transition.	
					if( ipHi==1 && ipLo==0 )
						continue;

					long index_of_nHi_P;
					/* if ipHi is collapsed level, index_of_nHi_P is ipHi */
					if( N_(ipHi) > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max )
						index_of_nHi_P = ipHi;
					else
						index_of_nHi_P = iso_sp[ipH_LIKE][nelem].QN2Index(N_(ipHi), 1, 2);

					bool lgSkip;
					if( N_(ipLo) > iso_sp[ipH_LIKE][nelem].n_HighestResolved_max || N_(ipLo) == N_(ipHi) )
						lgSkip = false;
					else
						lgSkip = !( ipHi == index_of_nHi_P && ipLo == index_of_nLo_S );

					if( lgSkip )
						continue;

					string comment_trans = "";
					if( LineSave.ipass == 0 )
					{
						comment_trans = iso_comment_tran_levels( ipISO, nelem, ipLo, ipHi );
					}
					PutLine(iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo), comment_trans.c_str());
				}
			}
		}
	}
	
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_hydro returns\n" );
	}
	return;
}

