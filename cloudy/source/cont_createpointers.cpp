/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ContCreatePointers set up pointers for lines and continua called by cloudy after input read in 
 * and continuum mesh has been set */
/*ipShells assign continuum energy pointers to shells for all atoms,
 * called by ContCreatePointers */
/*LimitSh sets upper energy limit to subshell integrations */
/*ContBandsCreate - read set of continuum bands to enter total emission into line stack*/
#include "cddefines.h"
#include "iso.h"
#include "secondaries.h"
#include "taulines.h"
#include "elementnames.h"
#include "ionbal.h"
#include "rt.h"
#include "opacity.h"
#include "yield.h"
#include "he.h"
#include "fe.h"
#include "rfield.h"
#include "oxy.h"
#include "trace.h"
#include "hmi.h"
#include "heavy.h"
#include "atmdat_adfa.h"
#include "ipoint.h"
#include "h2.h"
#include "continuum.h"
#include "freebound.h"
#include "two_photon.h"
#include "dense.h"
#include "lines_service.h"

/*LimitSh sets upper energy limit to subshell integrations */
STATIC long LimitSh(long int ion, 
  long int nshell, 
  long int nelem);

STATIC void ipShells(
	/* nelem is the atomic number on the C scale, Li is 2 */
	long int nelem);

/*ContBandsCreate - read set of continuum bands to enter total emission into line*/
STATIC void ContBandsCreate(
	/* chFile is optional filename, if void then use default bands,
	 * if not void then use file specified,
	 * return value is 0 for success, 1 for failure */
	 const char chFile[] );

STATIC inline void print_emline_fine( const char *LineGroup, const TransitionProxy &tr )
{
	fprintf( ioQQQ, "%-10s -> '%s'\t Energy Ang= %.4e\t Ryd= %.4e\t Fine= %ld\t fine_ener= %.4e\n",
			LineGroup,
			tr.chLabel().c_str(),
			tr.EnergyAng(),
			tr.EnergyRyd(),
			tr.Emis().ipFine(),
			rfield.fine_anu[ tr.Emis().ipFine() ] );
}

void ContCreatePointers(void)
{
	/* counter to say whether pointers have ever been evaluated */
	static int nCalled = 0;

	DEBUG_ENTRY( "ContCreatePointers()" );

	/* create the hydrogen atom for this core load, routine creates space then zeros it out
	 * on first call, on second and later calls it only zeros things out */
	iso_create();

	/* create internal static variables needed to do the H2 molecule */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->init();

	/* nCalled is local static variable defined 0 when defined. 
	 * it is incremented below, so that space only allocated one time per coreload. */
	if( nCalled > 0 )
	{
		if( trace.lgTrace )
			fprintf( ioQQQ, " ContCreatePointers called, not evaluating.\n" );
		return;
	}
	else
	{
		if( trace.lgTrace )
			fprintf( ioQQQ, " ContCreatePointers called first time.\n" );
		++nCalled;
	}

	for( long i=0; i < rfield.nflux_with_check; i++ )
	{
		/* this is array of labels for lines and continua, set to blanks at first */
		rfield.chContLabel[i] = "    ";
		rfield.chLineLabel[i] = "    ";
	}

	/* we will generate a set of array indices to ionization edges for
	 * the first thirty elements.  First set all array indices to
	 * totally bogus values so we will crash if misused */
	for( long nelem=0; nelem<LIMELM; ++nelem )
	{
		for( long ion=0; ion<LIMELM; ++ion )
		{
			for( long nshells=0; nshells<7; ++nshells )
			{
				for( long j=0; j<3; ++j )
				{
					opac.ipElement[nelem][ion][nshells][j] = INT_MIN;
				}
			}
		}
	}

	/* pointer to excited state of O+2 */
	opac.o3exc = 3.855;
	opac.ipo3exc[0] = ipContEnergy(opac.o3exc,"O3ex");

	/* main hydrogenic arrays - THIS OCCURS TWICE!! H and He here, then the
	 * remaining hydrogenic species near the bottom.  This is so that H and HeII get
	 * their labels stuffed into the arrays, and the rest of the hydrogenic series 
	 * get whatever is left over after the level 1 lines.
	 * to find second block, search on "ipZ=2" */
	/* NB note that no test for H or He being on exists here - we will always
	 * define the H and He arrays even when He is off, since we need to
	 * know where the he edges are for the bookkeeping that occurs in continuum
	 * binning routines */
	/* this loop is over H, He-like only - DO NOT CHANGE TO NISO */
	for( long ipISO=ipH_LIKE; ipISO<=ipHE_LIKE; ++ipISO )
	{
		/* this will be over HI, HeII, then HeI only (if enabled) */
		for( long nelem=ipISO; nelem < 2; nelem++ )
		{
			if( !dense.lgElmtOn[nelem] )
				continue;

			/* generate label for this ion */
			string chLab = chIonLbl( nelem+1, nelem+1-ipISO );

			/* array index for continuum edges for ground */
			iso_sp[ipISO][nelem].fb[0].ipIsoLevNIonCon = ipContEnergy(iso_sp[ipISO][nelem].fb[0].xIsoLevNIonRyd, chLab.c_str());
			for( long ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
			{
				/* array index for continuum edges for excited levels */
				iso_sp[ipISO][nelem].fb[ipHi].ipIsoLevNIonCon = ipContEnergy(iso_sp[ipISO][nelem].fb[ipHi].xIsoLevNIonRyd, chLab.c_str());

				/* define all line array indices */
				for( long ipLo=0; ipLo < ipHi; ipLo++ )
				{
					if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
						continue;

					/* some lines have negative or zero energy */
					/* >>chng 03 apr 22, this check was if less than or equal to zero,
					 * changed to lowest energy point so that ultra low energy transitions are
					 * not considered.	*/
					if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyRyd() < rfield.emm() )
						continue;

					/* some energies are negative for inverted levels */
					iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() = 
						ipLineEnergy(iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyRyd(), chLab.c_str(),
						iso_sp[ipISO][nelem].fb[ipLo].ipIsoLevNIonCon);
					iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ipFine() = 
						ipFineCont(iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyRyd() );
					/* check that energy scales are the same, to within energy resolution of arrays */
#					ifndef NDEBUG
					if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ipFine() > 0 )
					{
						realnum anuCoarse = rfield.anu(iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont()-1);
						realnum anuFine = rfield.fine_anu[iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ipFine()];
						realnum widCoarse = rfield.widflx(iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont()-1);
						realnum widFine = anuFine - rfield.fine_anu[iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ipFine()-1];
						realnum width = MAX2( widFine , widCoarse );
						/* NB - can only assert more than width of coarse cell 
						 * since close to ionization limit, coarse index is 
						 * kept below next ionization edge 
						 * >>chng 05 mar 02, pre factor below had been 1.5, chng to 2
						 * tripped near H grnd photo threshold */
						ASSERT( fabs(anuCoarse - anuFine) / anuCoarse < 
							2.*width/anuCoarse);
					}
#					endif
				}
			}/* ipHi loop */
		}/* nelem loop */
	}/* ipISO */

	/* we will define an array of either 1 or 0 to show whether photooelectron
	 * energy is large enough to produce secondary ionizations
	 * below 100eV no secondary ionization - this is the threshold*/
	secondaries.ipSecIon = ipoint(7.353);

	/* this is highest energy where k-shell opacities are counted
	 * can be adjusted with "set kshell" command */
	continuum.KshellLimit = ipoint(continuum.EnergyKshell);

	/* pointers for molecules
	 * H2+ dissociation energy 2.647 eV but cs small below 0.638 Ryd */
	opac.ih2pnt[0] = ipContEnergy(0.502,"H2+ ");
	opac.ih2pnt[1] = ipoint(1.03);
	// Excited state H2+
	opac.ih2pnt_ex[0] = ipContEnergy(0.052,"H2+*");
	opac.ih2pnt_ex[1] = ipoint(1.03);

	//pointers for most prominent PAH features
	{
		/* energies given to ipContEnergy are only to put lave in the right place
		 * wavelengths are rough observed values of blends
		 * 7.6 microns */
		(void) ipContEnergy(0.0117, "PAH " );

		/* feature near 6.2 microns */
		(void) ipContEnergy(0.0147, "PAH " );

		/* 3.3 microns */
		(void) ipContEnergy(0.028, "PAH " );

		/* 11.2 microns */
		(void) ipContEnergy(0.0080, "PAH " );

		/* 12.3 microns */
		(void) ipContEnergy(0.0077, "PAH " );

		/* 13.5 microns */
		(void) ipContEnergy(0.0069, "PAH " );
	}

	/* fix pointers for hydrogen and helium */

	/* pointer to Bowen 374A resonance line */
	he.ip660 = ipLineEnergy(1.38,"He 2",0);

	/* pointer to energy defining effect x-ray column */
	rt.ipxry = ipoint(73.5);

	/* pointer to Hminus edge at 0.754eV */
	hmi.iphmin = ipContEnergy(HMINUSIONPOT,"H-  ");

	/* pointer to threshold for H2 photoionization at 15.4 eV */
	fixit("need to generalize ionization energy and label!");
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->ip_photo_opac_thresh = ipContEnergy( 15.4/EVRYD , "H2  ");

	hmi.iheh1 = ipoint(1.6);
	hmi.iheh2 = ipoint(2.3);

	/* pointer to carbon k-shell ionization */
	opac.ipCKshell = ipoint(20.6);

	/* pointer to threshold for pair production */
	opac.ippr = ipoint(7.51155e4) + 1;

	/* pointer to x-ray - gamma ray bound; 100 kev */
	rfield.ipEnerGammaRay = ipoint(rfield.EnerGammaRay);

	/* confirm that labels are in correct location */
	ASSERT( strcmp( rfield.chContLabel[iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1].c_str(), "H  1" ) ==0 );
	ASSERT( strcmp( rfield.chContLabel[iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH2p].ipIsoLevNIonCon-1].c_str(), "H  1" ) ==0 );

	if( dense.lgElmtOn[ipHELIUM] )
		ASSERT( strcmp( rfield.chContLabel[iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].ipIsoLevNIonCon-1].c_str(), "He 2" ) ==0 );

	/* these are indices for centers of B and V filters,
	 * taken from table on page 202 of Allen, AQ, 3rd ed */
	/* the B filter array offset */
	rfield.ipB_filter = ipoint( RYDLAM / WL_B_FILT );
	/* the V filter array offset */
	rfield.ipV_filter = ipoint( RYDLAM / WL_V_FILT );

	/* these are the lower and upper bounds for the G0 radiation field
	 * used by Tielens & Hollenbach in their PDR work */
	rfield.ipG0_TH85_lo =  ipoint(  6.0 / EVRYD );
	rfield.ipG0_TH85_hi =  ipoint( 13.6 / EVRYD );

	/* this is the limits for Draine & Bertoldi Habing field */
	rfield.ipG0_DB96_lo =  ipoint(  RYDLAM / 1110. );
	rfield.ipG0_DB96_hi =  ipoint( RYDLAM / 912. );

	/* this is special form of G0 that could be used in some future version, for now,
	 * use default, TH85 */
	rfield.ipG0_spec_lo = ipoint(  6.0 / EVRYD );
	rfield.ipG0_spec_hi = ipoint( RYDLAM / 912. );

	/* this index is to 1000A to obtain the extinction at 1000A */
	rfield.ip1000A = ipoint( RYDLAM / 1000. );

	/* following order of elements is in roughly decreasing abundance
	 * when ipShells gets a cell for the valence IP it does so through
	 * ipContEnergy, which makes sure that no two ions get the same cell
	 * earliest elements have most precise ip mapping */

	/* set up shell pointers for hydrogen */
	long nelem = ipHYDROGEN;
	long ion = 0;

	/* the number of shells */
	Heavy.nsShells[nelem][0] = 1;

	/*pointer to ionization threshold in energy array*/
	Heavy.ipHeavy[nelem][ion] = iso_sp[ipH_LIKE][nelem].fb[ipH1s].ipIsoLevNIonCon;
	opac.ipElement[nelem][ion][0][0] = iso_sp[ipH_LIKE][nelem].fb[ipH1s].ipIsoLevNIonCon;

	/* upper limit to energy integration */
	opac.ipElement[nelem][ion][0][1] = rfield.nflux_with_check;

	if( dense.lgElmtOn[ipHELIUM] )
	{
		/* set up shell pointers for helium */
		long nelem = ipHELIUM;
		long ion = 0;

		/* the number of shells */
		Heavy.nsShells[nelem][0] = 1;

		/*pointer to ionization threshold in energy array*/
		Heavy.ipHeavy[nelem][ion] = iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon;
		opac.ipElement[nelem][ion][0][0] = iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon;

		/* upper limit to energy integration */
		opac.ipElement[nelem][ion][0][1] = rfield.nflux_with_check;

		/* (hydrogenic) helium ion */
		ion = 1;
		/* the number of shells */
		Heavy.nsShells[nelem][1] = 1;

		/*pointer to ionization threshold in energy array*/
		Heavy.ipHeavy[nelem][ion] = iso_sp[ipH_LIKE][nelem].fb[ipH1s].ipIsoLevNIonCon;
		opac.ipElement[nelem][ion][0][0] = iso_sp[ipH_LIKE][nelem].fb[ipH1s].ipIsoLevNIonCon;

		/* upper limit to energy integration */
		opac.ipElement[nelem][ion][0][1] = rfield.nflux_with_check;
	}

	/* check that ionization potential of neutral carbon valence shell is
	 * positive */
	ASSERT( atmdat.getIonPot(ipCARBON, 0) > 0. );

	/* now fill in all sub-shell ionization array indices for elements heavier than He,
	 * this must be done after previous loop on iso.ipIsoLevNIonCon[ipH_LIKE] since hydrogenic species use
	 * iso.ipIsoLevNIonCon[ipH_LIKE] rather than ipoint in getting array index within continuum array */
	for( long i=NISO; i<LIMELM; ++i )
	{
		/* i is the atomic number on the c scale, 2 for Li */
		if( dense.lgElmtOn[i] )
			ipShells(i);
	}

	/* most of these are set in ipShells, but not h-like or he-like, so do these here */
	Heavy.Valence_IP_Ryd[ipHYDROGEN][0] = atmdat.getIonPot(ipHYDROGEN, 0);
	Heavy.Valence_IP_Ryd[ipHELIUM][0] = atmdat.getIonPot(ipHELIUM, 0);
	Heavy.Valence_IP_Ryd[ipHELIUM][1] = atmdat.getIonPot(ipHELIUM, 1);
	for( long nelem=2; nelem<LIMELM; ++nelem )
	{
		Heavy.Valence_IP_Ryd[nelem][nelem-1] = atmdat.getIonPot(nelem, nelem-1);
		Heavy.Valence_IP_Ryd[nelem][nelem] = atmdat.getIonPot(nelem, nelem);
		if( dense.lgElmtOn[nelem])
		{
			/* now confirm that all are properly set */
			for( long j=0; j<=nelem; ++j )
			{
				ASSERT( Heavy.Valence_IP_Ryd[nelem][j] > 0.05 );
			}
			for( long j=0; j<nelem; ++j )
			{
				ASSERT( Heavy.Valence_IP_Ryd[nelem][j] < Heavy.Valence_IP_Ryd[nelem][j+1]);
			}
		}
	}

	/* array indices to bound Compton electron recoil ionization of all elements */
	for( long nelem=0; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem])
		{
			for( long ion=0; ion<nelem+1; ++ion )
			{
				/* this is the threshold energy to Compton ionize valence shell electrons */
				double energy = sqrt( Heavy.Valence_IP_Ryd[nelem][ion] * EN1RYD * ELECTRON_MASS * SPEEDLIGHT * SPEEDLIGHT ) / EN1RYD;
				/* the array index for this energy */
				ionbal.ipCompRecoil[nelem][ion] = ipoint( energy );
			}
		}
	}

	/* oxygen pointers for excited states
	 * IO3 is pointer to O++ exc state, is set above */
	oxy.i2d = ipoint(1.242);
	oxy.i2p = ipoint(1.367);
	opac.ipo1exc[0] = ipContEnergy(0.856,"O1ex");//2s^2 2p^4, ^1D level, J=2, energy relative to ground level is ~0.1446 Ry
	opac.ipo1exc[1] = ipoint(2.0);// upper limit to range of energies where opacity for 1D absorption is defined

	/* upper limit for excited state photoionization
	 * do not use ipContEnergy since only upper limit */
	opac.ipo3exc[1] = ipoint(5.0);

	/* upper level of 4363 */
	opac.o3exc3 = 3.646;
	opac.ipo3exc3[0] = ipContEnergy(opac.o3exc3,"O3ex");
	opac.ipo3exc3[1] = ipoint(5.0);

	/* >>chng 97 jan 27, move nitrogen after oxygen so that O gets the
	 * most accurate pointers
	 * Nitrogen
	 * in1(1) is thresh for photo from excited state */
	opac.in1[0] = ipContEnergy(0.893,"N1ex");

	/* upper limit */
	opac.in1[1] = ipoint(2.);

	if( (trace.lgTrace && trace.lgConBug) || (trace.lgTrace && trace.lgPointBug) )
	{
		fprintf( ioQQQ, "   ContCreatePointers:%ld energy cells used. N(1R):%4ld", 
				 rfield.nflux_with_check, 
				 iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon );
		if( dense.lgElmtOn[ipHELIUM] )
		{
			fprintf( ioQQQ, " N(1.8):%4ld  N(4Ryd):%4ld", 
					 iso_sp[ipHE_LIKE][ipHELIUM].fb[ipH1s].ipIsoLevNIonCon, 
					 iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].ipIsoLevNIonCon );
		}
		fprintf( ioQQQ, " N(O3)%4ld  N(x-ray):%5ld N(rcoil)%5ld\n", 
				 opac.ipo3exc[0], 
				 opac.ipCKshell, 
				 ionbal.ipCompRecoil[ipHYDROGEN][0] );

		fprintf( ioQQQ, "   ContCreatePointers: ipEnerGammaRay: %5ld IPPRpari produc%5ld\n", 
		  rfield.ipEnerGammaRay, opac.ippr );

		fprintf( ioQQQ, "   ContCreatePointers: H pointers;" );
		for( long i=0; i <= 6; i++ )
		{
			fprintf( ioQQQ, "%4ld%4ld", i, iso_sp[ipH_LIKE][ipHYDROGEN].fb[i].ipIsoLevNIonCon );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, "   ContCreatePointers: Oxy pnters;" );

		for( long i=1; i <= 8; i++ )
		{
			fprintf( ioQQQ, "%4ld%4ld", i, Heavy.ipHeavy[ipOXYGEN][i-1] );
		}
		fprintf( ioQQQ, "\n" );

	}

	/* Magnesium
	 * following is energy for phot of MG+ from exc state producing 2798 */
	opac.ipmgex = ipoint(0.779);

	/* lower, upper edges of Ca+ excited term photoionizaiton */
	opac.ica2ex[0] = ipContEnergy(0.72,"Ca2x");
	opac.ica2ex[1] = ipoint(1.);

	/* set up factors and pointers for Fe continuum fluorescence */
	fe.ipfe10 = ipoint(2.605);

	/* following is WL(CM)**2/(8PI) * conv fac for RYD to NU *A21 */
	fe.pfe10 = (realnum)(2.00e-18/rfield.widflx(fe.ipfe10-1));

	/* this is 353 pump, f=0.032 */
	fe.pfe11a = (realnum)(4.54e-19/rfield.widflx(fe.ipfe10-1));

	/* this is 341.1 f=0.012 */
	fe.pfe11b = (realnum)(2.75e-19/rfield.widflx(fe.ipfe10-1));
	fe.pfe14 = (realnum)(1.15e-18/rfield.widflx(fe.ipfe10-1));

	/* set up energy pointers for line optical depth arrays
	 * this also increments flux, sets other parameters for lines */

	/*Beginning of the dBaseLines*/
	for (int ipSpecies=0; ipSpecies < nSpecies; ++ipSpecies)
	{
		/*Put null terminated line label into chLab*/
		char chLab[NCHLAB];
		strncpy(chLab,dBaseSpecies[ipSpecies].chLabel,NCHLAB-1);
		chLab[NCHLAB-1]='\0';

		for( EmissionList::iterator em=dBaseTrans[ipSpecies].Emis().begin();
			  em != dBaseTrans[ipSpecies].Emis().end(); ++em)
		{
			if( ! rfield.isEnergyBound( Energy( (*em).Tran().EnergyWN(), "cm^-1" ) ) )
			{
				(*em).Tran().ipCont() = -1;
				(*em).ipFine() = -1;
				continue;
			}

			/* upper level lifetime to calculate the damping parameter.*/
			(*em).dampXvel() = (realnum)(1./
					dBaseStates[ipSpecies][em->Tran().ipHi()].lifetime()/em->Tran().EnergyWN()/PI4);
			(*em).damp() = -1000.0;

			// Data intended for ADAS have a large number of transitions with Aul = 1e-30 - these are
			// not real, but are there to keep the ADAS codes happy.  We do not want to transfer them.
			static const double minAul = 1e-29;
			if( (*em).Aul() > minAul )
			{
				(*em).Tran().ipCont() = ipLineEnergy((*em).Tran().EnergyRyd(), chLab ,0);
				(*em).ipFine() = ipFineCont((*em).Tran().EnergyRyd() );
			}
			else
			{
				(*em).Tran().ipCont() = -1;
				(*em).ipFine() = -1;
			}

			// these are edge cases which did pass the minAul above, but still underflowed
			// when converted to (real)gf - conversion depends on energy of transition so
			// Aul - gf mapping is not linear
			if ((*em).gf() <= 0.0)
			{
				(*em).Tran().ipCont() = -1;
				(*em).ipFine() = -1;
			}

			/* derive the abs coefficient, call to function is gf, wl (A), g_low */
			(*em).opacity() = 
				(realnum)(abscf((*em).gf(),(*em).Tran().EnergyWN(),
									 (*(*em).Tran().Lo()).g()));
		}
	}
	/*end of the dBaseLines*/

	/* set the ipCont struc element for the H2 molecule */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_ContPoint();

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* do remaining part of the iso sequences */
		for( long nelem=2; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem])
			{
				string chLab = chIonLbl( nelem+1, nelem+1-ipISO );

				/* array index for continuum edges */
				iso_sp[ipISO][nelem].fb[0].ipIsoLevNIonCon = 
					ipContEnergy(iso_sp[ipISO][nelem].fb[0].xIsoLevNIonRyd, chLab.c_str());

				for( long ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					/* array index for continuum edges */
					iso_sp[ipISO][nelem].fb[ipHi].ipIsoLevNIonCon = ipContEnergy(iso_sp[ipISO][nelem].fb[ipHi].xIsoLevNIonRyd, chLab.c_str());

					/* define all line pointers */
					for( long ipLo=0; ipLo < ipHi; ipLo++ )
					{

						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
							continue;

						/* some lines have negative or zero energy */
						/* >>chng 03 apr 22, this check was if less than or equal to zero,
						 * changed to lowest energy point so that ultra low energy transitions are
						 * not considered.	*/
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyRyd() < rfield.emm() )
							continue;

						iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() = 
							ipLineEnergy(iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyRyd(), chLab.c_str(),
							iso_sp[ipISO][nelem].fb[ipLo].ipIsoLevNIonCon);
						iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().ipFine() = 
							ipFineCont(iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyRyd() );
					}
				}
				iso_sp[ipISO][nelem].fb[0].ipIsoLevNIonCon = ipContEnergy(iso_sp[ipISO][nelem].fb[0].xIsoLevNIonRyd, chLab.c_str());
			}
		}
	}
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* this will be over HI, HeII, then HeI only */
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem])
			{
				/* these are the extra Lyman lines */
				for( long ipHi=2; ipHi < iso_ctrl.nLyman_alloc[ipISO]; ipHi++ )
				{
					long ipLo = 0;
					/* some energies are negative for inverted levels */
					char chLab[NCHLAB];
					strncpy(chLab,"LyEx",NCHLAB-1);
					chLab[NCHLAB-1]='\0';
					TransitionList::iterator tr = ExtraLymanLines[ipISO][nelem].begin()+ipExtraLymanLines[ipISO][nelem][ipHi];
					(*tr).ipCont() = 
						ipLineEnergy((*tr).EnergyRyd() , chLab ,
						iso_sp[ipISO][nelem].fb[ipLo].ipIsoLevNIonCon);

					(*tr).Emis().ipFine() = 
						ipFineCont((*tr).EnergyRyd() );
				}

				if( iso_ctrl.lgDielRecom[ipISO] )
				{
					ASSERT( ipISO>ipH_LIKE );
					for( long ipLo=0; ipLo<iso_sp[ipISO][nelem].numLevels_max; ipLo++ )
					{
						char chLab[NCHLAB];
						strncpy(chLab,"SatL",NCHLAB-1);
						chLab[NCHLAB-1]='\0';
						SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][ipLo]].ipCont() = ipLineEnergy(
							SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][ipLo]].EnergyRyd() , chLab , 
							0);

						SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][ipLo]].Emis().ipFine() =  
							ipFineCont(SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][ipLo]].EnergyRyd() );
					}
				}
			}
		}
	}

	fixit("is this redundant?");
	/* for He-like sequence the majority of the transitions are bogus - A set to special value in this case */
	{
		long ipISO = ipHE_LIKE;
		/* do remaining part of the he iso sequence */
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem])
			{
				for( long ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					for( long ipLo=0; ipLo < ipHi; ipLo++ )
					{
						TransitionProxy tr = iso_sp[ipISO][nelem].trans(ipHi,ipLo);
						if( tr.ipCont() <= 0 )
							continue;
					
						if( fabs(tr.Emis().Aul() - iso_ctrl.SmallA) < SMALLFLOAT )
						{
							/* iso_ctrl.SmallA is value assigned to bogus transitions */
							tr.ipCont() = -1;
							tr.Emis().ipFine() = -1;
						}
					}
				}
			}
		}
	}

	/* inner shell transitions */
	for( size_t i=0; i<UTALines.size(); ++i )
	{
		if( UTALines[i].Emis().Aul() > 0. )
		{

			// dampXvel is derived in atmdat_readin because autoionization rates
			// must be included in the total rate for the damping constant
			ASSERT( UTALines[i].Emis().dampXvel() >0. );

			/* derive the abs coefficient, call to function is gf, wl (A), g_low */
			UTALines[i].Emis().opacity() = 
				(realnum)(abscf( UTALines[i].Emis().gf(), UTALines[i].EnergyWN(), (*UTALines[i].Lo()).g()));

			/* get pointer to energy in continuum mesh */
			UTALines[i].ipCont() = ipLineEnergy(UTALines[i].EnergyRyd(), chIonLbl(UTALines[i]).c_str(),0 );
			UTALines[i].Emis().ipFine() = ipFineCont(UTALines[i].EnergyRyd() );
			{
				enum{ DEBUG_LOC = false };
				if( DEBUG_LOC && UTALines[i].chLabel() == "Ar 7 43.5239A" )
				{
					print_emline_fine( "UTA", UTALines[i] );
				}
			}

			/* find heating per absorption,
			 * first find threshold for this shell in ergs */
			/* ionization threshold in erg */
			double thresh = Heavy.Valence_IP_Ryd[(*UTALines[i].Hi()).nelem()-1][(*UTALines[i].Hi()).IonStg()-1] *EN1RYD;
			UTALines[i].Coll().heat() = (UTALines[i].EnergyErg()-thresh);
			ASSERT( UTALines[i].Coll().heat()> 0. );
		}
	}

	/* level 2 heavy element lines */
	for( long i=0; i < nWindLine; i++ )
	{
		/* derive the A, call to function is gf, wl (A), g_up */
		TauLine2[i].Emis().Aul() = 
			(realnum)(eina(TauLine2[i].Emis().gf(),
		  TauLine2[i].EnergyWN(),(*TauLine2[i].Hi()).g()));

		/* coefficient needed for damping constant - units cm s-1 */
		TauLine2[i].Emis().dampXvel() = 
			(realnum)(TauLine2[i].Emis().Aul()/
		  TauLine2[i].EnergyWN()/PI4);

		/* derive the abs coefficient, call to function is gf, wl (A), g_low */
		TauLine2[i].Emis().opacity() = 
			(realnum)(abscf(TauLine2[i].Emis().gf(),
		  TauLine2[i].EnergyWN(),(*TauLine2[i].Lo()).g()));

		/* get pointer to energy in continuum mesh */
		TauLine2[i].ipCont() = ipLineEnergy(TauLine2[i].EnergyRyd(), chIonLbl(TauLine2[i]).c_str(),0 );
		TauLine2[i].Emis().ipFine() = ipFineCont(TauLine2[i].EnergyRyd() );
		/*if( TauLine2[i].ipCont()==751 )
			fprintf(ioQQQ,"( atom_level2 %s\n", chLab);*/
	}

	/* hyperfine structure lines */
	for( size_t i=0; i < HFLines.size(); i++ )
	{
		ASSERT( HFLines[i].Emis().Aul() > 0. );
		/* coefficient needed for damping constant */
		HFLines[i].Emis().dampXvel() = 
			(realnum)(HFLines[i].Emis().Aul()/
			HFLines[i].EnergyWN()/PI4);
		HFLines[i].Emis().damp() = 1e-20f;

		/* derive the abs coefficient, call to function is gf, wl (A), g_low */
		HFLines[i].Emis().opacity() = 
			(realnum)(abscf(HFLines[i].Emis().gf(),
			HFLines[i].EnergyWN(),
			(*HFLines[i].Lo()).g()));
		/* gf from this and 21 cm do not agree, A for HFS is 10x larger than level1 dat */
		/*fprintf(ioQQQ,"HFLinesss\t%li\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			i,HFLines[i].Emis().opacity() , HFLines[i].Emis().gf() , HFLines[i].Emis().Aul() , HFLines[i].EnergyWN(),(*HFLines[i].Lo()).g());*/

		/* get pointer to energy in continuum mesh */
		HFLines[i].ipCont() = ipLineEnergy(HFLines[i].EnergyRyd() , chIonLbl(HFLines[i]).c_str(),0 );
		HFLines[i].Emis().ipFine() = ipFineCont(HFLines[i].EnergyRyd() );
	}

	/* the group of inner shell fluorescent lines */
	for( long i=0; i < t_yield::Inst().nlines(); ++i )
	{
		string chLab = chIonLbl( t_yield::Inst().nelem(i)+1, t_yield::Inst().ion_emit(i)+1 );
		t_yield::Inst().set_ipoint( i, ipLineEnergy( t_yield::Inst().energy(i) , chLab.c_str() , 0 ) );
	}

	/* ================================================================================== */
	/*        two photon two-photon 2-nu 2 nu 2 photon 2-photon                           */

	/* now loop over the two iso-sequences with two photon continua */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* set up two photon emission */
		for( long nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				// upper level for two-photon emission in H and He iso sequences 
				// The 1+ipISO is not rigorous, it just works for H- and He-like 
				const int TwoS = (1+ipISO);
				double Aul;
				/* 2s two photon */
				if( ipISO==ipH_LIKE )
					Aul =  8.226*powi((double)(nelem+1.),6);
				else 
				{
					ASSERT( ipISO==ipHE_LIKE );
					/* >>refer	Helike	2pho	Derevianko, A., & Johnson, W. R. 1997, Phys. Rev. A, 56, 1288
					 * numbers are not explicitly given in this paper for Z=21-24,26,27,and 29.
					 * So numbers given here are interpolated.	*/
					fixit("where is 51.02 from? Value is 51.3 from the Derevianko & Johnson paper cited above.");
					const double As2nuFrom1S[29] = {51.02,1940.,1.82E+04,9.21E+04,3.30E+05,9.44E+05,2.31E+06,5.03E+06,1.00E+07,
					1.86E+07,3.25E+07,5.42E+07,8.69E+07,1.34E+08,2.02E+08,2.96E+08,4.23E+08,5.93E+08,8.16E+08,
					1.08E+09,1.43E+09,1.88E+09,2.43E+09,3.25E+09,3.95E+09,4.96E+09,6.52E+09,7.62E+09,9.94E+09};
					Aul =  As2nuFrom1S[nelem-1];
				}
	
				TwoPhotonSetup( iso_sp[ipISO][nelem].TwoNu, TwoS, 0, 
					Aul,
					iso_sp[ipISO][nelem].trans(TwoS,0),
					ipISO, nelem );
			}
		}
	}

	// add He-like 2nu 2^3S - 1^1S
	{
		const long ipISO = ipHE_LIKE;
		for( long nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				/* Important clarification, according to Derevianko & Johnson (see ref above), 2^3S can decay
				 * to ground in one of two ways: through a two-photon process, or through a single-photon M1 decay,
				 * but the M1 rates are about 10^4 greater that the two-photon decays throughout the entire
				 * sequence.  Thus these numbers, are much weaker than the effective decay rate, but should probably
				 * be treated in as a two-photon decay at some point	*/
				// >> refer	He	As	Drake, G. W. F., Victor, G. A., & Dalgarno, A. 1969, Physical Review, 180, 25			
				const double As2nuFrom3S[29] = {4.09e-9,1.25E-06,5.53E-05,8.93E-04,8.05E-03,4.95E-02,2.33E-01,8.94E-01,2.95E+00,
					8.59E+00,2.26E+01,5.49E+01,1.24E+02,2.64E+02,5.33E+02,1.03E+03,1.91E+03,3.41E+03,5.91E+03,
					9.20E+03,1.50E+04,2.39E+04,3.72E+04,6.27E+04,8.57E+04,1.27E+05,2.04E+05,2.66E+05,4.17E+05};
		
				TwoPhotonSetup( iso_sp[ipISO][nelem].TwoNu, ipHe2s3S, ipHe1s1S,
					As2nuFrom3S[nelem-1],
					iso_sp[ipISO][nelem].trans(ipHe2s3S,ipHe1s1S),
					ipISO, nelem );
			}
		}
	}

	{
		/* this is an option to print out one of the two photon continua */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{	
			const int nCRS = 21;
			double ener[nCRS]={
			  0.,     0.03738,  0.07506,  0.1124,  0.1498,  0.1875,
			  0.225,  0.263,    0.300,    0.3373,  0.375,   0.4127,
			  0.4500, 0.487,    0.525,    0.5625,  0.6002,  0.6376,
			  0.6749, 0.7126,   0.75};

			long nelem = ipHYDROGEN;
			long ipISO = ipHYDROGEN;
			two_photon& tnu = iso_sp[ipISO][nelem].TwoNu[0];

			double limit = tnu.ipTwoPhoE;

			for( long i=0; i < nCRS; i++ )
			{
				fprintf(ioQQQ,"%.3e\t%.3e\n", ener[i] , 
					atmdat_2phot_shapefunction( ener[i]/0.75, ipISO, nelem ) );
			}

			double xnew = 0.;
			/** \todo	2	what are we trying to print here?	*/
			for( long i=0; i < limit; i++ )
			{
				double fach = tnu.As2nu[i]/2.*rfield.anu2(i)/rfield.widflx(i)*EN1RYD;
				fprintf(ioQQQ,"%.3e\t%.3e\t%.3e\n", 
					rfield.anu(i) , 
					tnu.As2nu[i] / rfield.widflx(i) , 
					fach );
				xnew += tnu.As2nu[i];
			}
			fprintf(ioQQQ," sum is %.3e\n", xnew );
			cdEXIT(EXIT_FAILURE);
		}
	}

	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{	
			for( long i=0; i<11; ++i )
			{
				(*TauDummy).WLAng() = (realnum)(PI * exp10((double)i));
				fprintf(ioQQQ,"%.2f\t%s\n", (*TauDummy).WLAng() , chLineLbl(*TauDummy).c_str());
			}
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* option to print out whole thing with "trace lines" command */
	if( trace.lgTrLine )
	{
		fprintf( ioQQQ, "       WL(Ang)   E(RYD)   IP   gl  gu      gf       A        damp     abs K\n" );

		/*Atomic Or Molecular lines*/
		for (int ipSpecies=0; ipSpecies < nSpecies; ++ipSpecies)
		{
			for( EmissionList::iterator em=dBaseTrans[ipSpecies].Emis().begin();
				  em != dBaseTrans[ipSpecies].Emis().end(); ++em)
			{
				long iWL_Ang = (long)(*em).Tran().WLAng();
				
				if( iWL_Ang > 1000000 )
				{
					iWL_Ang /= 10000;
				}
				else if( iWL_Ang > 10000 )
				{
					iWL_Ang /= 1000;
				}
				fprintf( ioQQQ, " %10.10s%5ld%10.3e %4li%4ld%4ld%10.2e%10.2e%10.2e%10.2e\n", 
							chLineLbl((*em).Tran()).c_str(), iWL_Ang, RYDLAM/(*em).Tran().WLAng(), 
							(*em).Tran().ipCont(), (long)((*(*em).Tran().Lo()).g()), 
							(long)((*(*em).Tran().Hi()).g()),(*em).gf(), 
							(*em).Aul(),(*em).dampXvel(), 
							(*em).opacity());
			}
		}

		for( long i=0; i < nWindLine; i++ )
		{
			long iWL_Ang = (long)TauLine2[i].WLAng();

			if( iWL_Ang > 1000000 )
			{
				iWL_Ang /= 10000;
			}
			else if( iWL_Ang > 10000 )
			{
				iWL_Ang /= 1000;
			}
			fprintf( ioQQQ, " %10.10s%5ld%10.3e %4li%4ld%4ld%10.2e%10.2e%10.2e%10.2e\n", 
			  chLineLbl(TauLine2[i]).c_str(), iWL_Ang, RYDLAM/TauLine2[i].WLAng(), 
			  TauLine2[i].ipCont(), (long)((*TauLine2[i].Lo()).g()), 
			  (long)((*TauLine2[i].Hi()).g()), TauLine2[i].Emis().gf(), 
			  TauLine2[i].Emis().Aul(), TauLine2[i].Emis().dampXvel(), 
			  TauLine2[i].Emis().opacity() );
		}
		for( size_t i=0; i < HFLines.size(); i++ )
		{
			long iWL_Ang = (long)HFLines[i].WLAng();

			if( iWL_Ang > 1000000 )
			{
				iWL_Ang /= 10000;
			}
			else if( iWL_Ang > 10000 )
			{
				iWL_Ang /= 1000;
			}
			fprintf( ioQQQ, " %10.10s%5ld%10.3e %4li%4ld%4ld%10.2e%10.2e%10.2e%10.2e\n", 
			  chLineLbl(HFLines[i]).c_str(), iWL_Ang, RYDLAM/HFLines[i].WLAng(), 
			  HFLines[i].ipCont(), (long)((*HFLines[i].Lo()).g()), 
			  (long)((*HFLines[i].Hi()).g()), HFLines[i].Emis().gf(), 
			  HFLines[i].Emis().Aul(), HFLines[i].Emis().dampXvel(), 
			  HFLines[i].Emis().opacity() );
		}
	}

	/* this is an option to kill fine structure line optical depths */
	if( !rt.lgFstOn )
	{
		/*Atomic or Molecular Lines-Humeshkar Nemala*/
		for (int ipSpecies=0; ipSpecies < nSpecies; ++ipSpecies)
		{
			for( EmissionList::iterator em=dBaseTrans[ipSpecies].Emis().begin();
				  em != dBaseTrans[ipSpecies].Emis().end(); ++em)
			{
				if((*em).Tran().EnergyWN() < 10000. )
				{
					(*em).opacity() = 0.;
				}
			}
		}
	}

	/* read in continuum bands data set */
	ContBandsCreate( "" );

	/* we're done adding lines and states to the stacks.  
	 * This flag is used to make sure we never add them again in this coreload. */
	lgLinesAdded = true;
	lgStatesAdded = true;
	
	checkTransitionListOfLists(AllTransitions);

	return;
}

/*ipShells assign continuum energy pointers to shells for all atoms,
 * called by ContCreatePointers */
STATIC void ipShells(
	/* nelem is the atomic number on the C scale, Li is 2 */
	long int nelem)
{
	long int 
	  imax, 
	  ion, 
	  nelec, 
	  ns, 
	  nshell;
	/* following value cannot be used - will be set to proper threshold */
	double thresh=-DBL_MAX;

	DEBUG_ENTRY( "ipShells()" );

	ASSERT( nelem >= NISO);
	ASSERT( nelem < LIMELM );

	/* fills in pointers to valence shell ionization threshold
	 * PH1(a,b,c,d)
	 * a=1 => thresh, others fitting parameters
	 * b atomic number
	 * c number of electrons
	 * d shell number 7-1 */

	/* threshold in Ryd
	 * ion=0 for atom, up to nelem-1 for helium like, hydrogenic is elsewhere */
	for( ion=0; ion < nelem; ion++ )
	{
		string chLab = chIonLbl( nelem+1, ion+1 );

		/* this is the iso sequence - must not redo sequence if done as iso */
		long int ipISO = nelem-ion;

		/* number of bound electrons */
		nelec = ipISO+1;

		/* nsShells(nelem,ion) is the number of shells for ion with nelec electrons,
		 * physical not c scale */
		imax = Heavy.nsShells[nelem][ion];

		/* loop on all inner shells, valence shell */
		for( nshell=0; nshell < imax; nshell++ )
		{
			/* ionization potential of this shell in rydbergs */
			thresh = t_ADfA::Inst().getEthresh(nshell+1,nelec,nelem+1)/EVRYD;
			if( thresh <= 0.1 )
			{
				/* negative ip shell does not exist, set upper limit
				 * to less than lower limit so this never looped upon
				 * these are used as flags by LimitSh to check whether
				 * this is a real shell - if 1 or 2 is changed - change LimitSh!! */
				opac.ipElement[nelem][ion][nshell][0] = 2;
				opac.ipElement[nelem][ion][nshell][1] = 1;
			}
			else
			{
				/* this is lower bound to energy range for this shell */
				/* >>chng 02 may 27, change to version of ip with label, so that
				 * inner shell edges will appear */
				/*opac.ipElement[nelem][ion][nshell][0] = ipoint(thresh);*/
				opac.ipElement[nelem][ion][nshell][0] = 
					ipContEnergy( thresh , chLab.c_str() );

				/* this is upper bound to energy range for this shell 
				 * LimitSh is an integer function, returns pointer
				 * to threshold of next major shell.  For k-shell it
				 * returns the values KshellLimit, default=7.35e4
				 * >>chng 96 sep 26, had been below, result zero cross sec at 
				 * many energies where opacity project did not produce state specific 
				 * cross section */
				opac.ipElement[nelem][ion][nshell][1] = 
					LimitSh(ion+1,  nshell+1,nelem+1);
				ASSERT( opac.ipElement[nelem][ion][nshell][1] > 0);
			}
		}

		ASSERT( imax > 0 && imax <= 7 );

		/* this will be index pointing to valence edge */
		/* [0] is pointer to threshold in energy array */
		opac.ipElement[nelem][ion][imax-1][0] = 
			ipContEnergy(thresh, chLab.c_str());

		/* pointer to valence electron ionization potential */
		Heavy.ipHeavy[nelem][ion] = opac.ipElement[nelem][ion][imax-1][0];
		ASSERT( Heavy.ipHeavy[nelem][ion]>0 );

		/* ionization potential of valence shell in Ryd 
		 * thresh was evaluated above, now has last value, the valence shell */
		Heavy.Valence_IP_Ryd[nelem][ion] = thresh;

		Heavy.xLyaHeavy[nelem][ion] = 0.;
		if( ipISO >= NISO )
		{
			/* this is set of 3/4 of valence shell IP, this is important
			 * source of ots deep in cloud */
			Heavy.ipLyHeavy[nelem][ion] = 
				ipLineEnergy(thresh*0.75, chLab.c_str(), 0);

			Heavy.ipBalHeavy[nelem][ion] = 
				ipLineEnergy(thresh*0.25, chLab.c_str(), 0);
		}
		else
		{
			/* do not treat this simple way since done exactly with iso 
			 * sequences */
			Heavy.ipLyHeavy[nelem][ion] = -1;
			Heavy.ipBalHeavy[nelem][ion] = -1;
		}
	}

	/* above loop did up to hydrogenic, now do hydrogenic - 
	 * hydrogenic is special since arrays already set up */
	Heavy.nsShells[nelem][nelem] = 1;

	/* this is lower limit to range */
	/* hydrogenic photoionization set to special hydro array 
	 * this is pointer to threshold energy */
	/* this statement is in ContCreatePointers but has not been done when this routine called */
	/*iso_sp[ipH_LIKE][ipZ].fb[ipLo].ipIsoLevNIonCon = ipContEnergy(iso_sp[ipH_LIKE][ipZ].fb[ipLo].xIsoLevNIonRyd,chLab);*/
	/*opac.ipElement[nelem][nelem][0][0] = iso_sp[ipH_LIKE][nelem].fb[ipH1s].ipIsoLevNIonCon;*/
	opac.ipElement[nelem][nelem][0][0] = ipoint( t_ADfA::Inst().getEthresh(1,1,nelem+1)/EVRYD );
	ASSERT( opac.ipElement[nelem][nelem][0][0] > 0 );

	/* this is the high-energy limit */
	opac.ipElement[nelem][nelem][0][1] = continuum.KshellLimit;

	Heavy.ipHeavy[nelem][nelem] = opac.ipElement[nelem][nelem][0][0];

	/* this is for backwards computability with Cambridge code */
	if( trace.lgTrace && trace.lgPointBug )
	{
		for( ion=0; ion < (nelem+1); ion++ )
		{
			fprintf( ioQQQ, "Ion:%3ld%3ld %2.2s%2.2s total shells:%3ld\n", 
			  nelem, ion+1, elementnames.chElementSym[nelem], elementnames.chIonStage[ion]
			  , Heavy.nsShells[nelem][ion] );
			for( ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
			{
				fprintf( ioQQQ, " shell%3ld %2.2s range eV%10.2e-%8.2e\n", 
				  ns+1, Heavy.chShell[ns], rfield.anu(opac.ipElement[nelem][ion][ns][0]-1)*
				  EVRYD, rfield.anu(opac.ipElement[nelem][ion][ns][1]-1)*EVRYD );
			}
		}
	}
	return;
}

/*LimitSh sets upper energy limit to subshell integrations */
STATIC long LimitSh(long int ion, 
  long int nshell, 
  long int nelem)
{
	long int LimitSh_v;

	DEBUG_ENTRY( "LimitSh()" );

	/* this routine returns the high-energy limit to the energy range
	 * for photoionization of a given shell
	 * */
	if( nshell == 1 )
	{
		/* this limit is high-energy limit to code unless changed with set kshell */
		LimitSh_v = continuum.KshellLimit;

	}
	else if( nshell == 2 )
	{
		/* this is 2s shell, upper limit is 1s
		 * >>chng 96 oct 08, up to high-energy limit
		 * LimitSh = ipElement(nelem,ion , 1,1)-1 */
		LimitSh_v = continuum.KshellLimit;

	}
	else if( nshell == 3 )
	{
		/* this is 2p shell, upper limit is 1s
		 * >>chng 96 oct 08, up to high-energy limit
		 * LimitSh = ipElement(nelem,ion , 1,1)-1 */
		LimitSh_v = continuum.KshellLimit;

	}
	else if( nshell == 4 )
	{
		/* this is 3s shell, upper limit is 2p
		 * >>chng 96 oct 08, up to K-shell edge
		 * LimitSh = ipElement(nelem,ion , 3,1)-1 */
		LimitSh_v = opac.ipElement[nelem-1][ion-1][0][0] - 1;

	}
	else if( nshell == 5 )
	{
		/* this is 3p shell, upper limit is 2p
		 * >>chng 96 oct 08, up to K-shell edge
		 * LimitSh = ipElement(nelem,ion , 3,1)-1 */
		LimitSh_v = opac.ipElement[nelem-1][ion-1][0][0] - 1;

	}
	else if( nshell == 6 )
	{
		/* this is 3d shell, upper limit is 2p
		 * >>chng 96 oct 08, up to K-shell edge
		 * LimitSh = ipElement(nelem,ion , 3,1)-1 */
		LimitSh_v = opac.ipElement[nelem-1][ion-1][0][0] - 1;

	}
	else if( nshell == 7 )
	{
		/* this is 4s shell, upper limit is 3d */
		if( opac.ipElement[nelem-1][ion-1][5][0] < 3 )
		{
			/* this is check for empty shell 6, 3d
			 * if so then set to 3p instead */
			LimitSh_v = opac.ipElement[nelem-1][ion-1][4][0] - 
			  1;
		}
		else
		{
			LimitSh_v = opac.ipElement[nelem-1][ion-1][5][0] - 
			  1;
		}
		/* >>chng 96 sep 26, set upper limit down to 2s */
		LimitSh_v = opac.ipElement[nelem-1][ion-1][2][0] - 1;

	}
	else
	{
		fprintf( ioQQQ, " LimitSh cannot handle nshell as large as%4ld\n", 
		  nshell );
		cdEXIT(EXIT_FAILURE);
	}
	return LimitSh_v;
}

/*ContBandsCreate - read set of continuum bands to enter total emission into line*/
STATIC void ContBandsCreate(
	/* chFile is optional filename, if void then use default bands,
	 * if not void then use file specified,
	 * return value is 0 for success, 1 for failure */
	 const char chFile[] )
{
	/* keep track of whether we have been called - want to be
	 * called a total of one time */
	static bool lgCalled=false;

	DEBUG_ENTRY( "ContBandsCreate()" );

	/* do nothing if second or later call*/
	if( lgCalled )
	{
		/* success */
		return;
	}
	lgCalled = true;

	/* use default filename if void string, else use file specified */
	const char* chFilename = ( strlen(chFile) == 0 ) ? "continuum_bands.ini" : chFile;

	/* get continuum band data  */
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " ContBandsCreate opening %s:", chFilename );
	}

	FILE *ioDATA = open_data( chFilename, "r" );

	/* now count how many bands are in the file */
	continuum.nContBand = 0;

	/* first line is a magic number and does not count as a band*/
	string chLine;
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " ContBandsCreate could not read first line of %s.\n", chFilename );
		cdEXIT(EXIT_FAILURE);
	}
	while( read_whole_line( chLine, ioDATA ) )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#' )
			++continuum.nContBand;
	}

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " ContBandsCreate could not rewind %s.\n", chFilename );
		cdEXIT(EXIT_FAILURE);
	}

	continuum.ContBandWavelength.resize(continuum.nContBand);
	continuum.chContBandLabels.resize(continuum.nContBand);
	continuum.ipContBandLow.resize(continuum.nContBand);
	continuum.ipContBandHi.resize(continuum.nContBand);
	continuum.BandEdgeCorrLow.resize(continuum.nContBand);
	continuum.BandEdgeCorrHi.resize(continuum.nContBand);

	/* first line is a versioning magic number - now confirm that it is valid */
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " ContBandsCreate could not read first line of %s.\n", chFilename );
		cdEXIT(EXIT_FAILURE);
	}
	/* bands_continuum magic number here <- this string is in band_continuum.dat
	 * with comment to search for this to find magic number  */

	// the magic number at the start of the data file
	const long myr = 17, mmo = 6, mdy = 27;

	long i = 1;
	bool lgEOL;
	long m1 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long m2 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long m3 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	if( m1 != myr || m2 != mmo || m3 != mdy )
	{
		fprintf( ioQQQ, 
			 " ContBandsCreate: the version of the data file %s I found (%li %li %li)is not the current version (%li %li %li).\n", 
			 chFilename ,
			 m1 , m2 , m3 ,
			 myr , mmo , mdy );
		fprintf( ioQQQ, 
			 " ContBandsCreate: you need to update this file.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* now read in data again, but save it this time */
	long k = 0;
	while( read_whole_line( chLine, ioDATA ) )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#' )
		{
			/* copy 4 char label */
			continuum.chContBandLabels[k] = chLine.substr(0,4);

			/* now get central band wavelength 
			 * >>chng 06 aug 11 from 4 to 6, the first 4 char are labels and
			 * these can contain numbers, next comes a space, then the number */
			i = 6;
			continuum.ContBandWavelength[k] = (realnum)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
			/* >>chng 06 feb 21, multiply by 1e4 to convert micron wavelength into Angstroms,
			 * which is assumed by the code.  before this correction the band centroid 
			 * wavelength was given in the output incorrectly listed as Angstroms.
			 * results were correct just label was wrong */
			continuum.ContBandWavelength[k] *= 1e4f;

			/* these are short and long wave limits, which are high and
			 * low energy limits - these are now wl in microns but are
			 * converted to Angstroms */
			double xHi = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL)*1e4;
			double xLow = FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL)*1e4;
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been 3 numbers on this band line.   Sorry.\n" );
				fprintf( ioQQQ, " string==%s==\n" ,chLine.c_str() );
				cdEXIT(EXIT_FAILURE);
			}

			{
				enum {DEBUG_LOC=false};
				if( DEBUG_LOC )
				{
					fprintf(ioQQQ, "READ:%s\n", chLine.c_str() );
					fprintf(ioQQQ, "GOT: %s %g %g %g\n",continuum.chContBandLabels[k].c_str(),
					continuum.ContBandWavelength[k] , xHi , xLow );
				}
			}

			/* make sure bands bounds are in correct order, shorter - longer wavelength*/
			if( xHi >= xLow )
			{
				fprintf( ioQQQ, " ContBandWavelength band %li "
					"edges are in improper order.\n" ,k);
				fprintf(ioQQQ,"band: %s %.3e %.3e %.3e \n",
						continuum.chContBandLabels[k].c_str(),
						continuum.ContBandWavelength[k],
						xHi , 
						xLow);
				cdEXIT(EXIT_FAILURE);
			}

			// check that central wavelength is indeed between the limits
			// xHi & xLow are hi and low energy limits to band so logic reversed
			if( continuum.ContBandWavelength[k] < xHi ||
				continuum.ContBandWavelength[k] > xLow )
			{
				fprintf( ioQQQ, " ContBandWavelength band number %li, "
					"central wavelength not within band.\n" ,k);
				fprintf(ioQQQ,"band ID:%s WL %.3e microns, band bounds %.3e to %.3e microns\n",
						continuum.chContBandLabels[k].c_str(),
						continuum.ContBandWavelength[k],
						xLow , xHi );
				cdEXIT(EXIT_FAILURE);
			}

			/* get continuum index - RYDLAM is 911.6A = 1 Ryd so 1e4 converts 
			 * micron to Angstrom - xHi is high energy (not wavelength)
			 * edge of the band */
			continuum.ipContBandHi[k] = ipoint( RYDLAM / xHi );
			continuum.ipContBandLow[k] = ipoint( RYDLAM / xLow );

			/* fraction of first and last bin to include */
			continuum.BandEdgeCorrLow[k] = 
				(rfield.anumax(continuum.ipContBandLow[k]-1)-
				(realnum)(RYDLAM/xLow)) /
				rfield.widflx(continuum.ipContBandLow[k]-1);
			ASSERT( continuum.BandEdgeCorrLow[k]>=0. && continuum.BandEdgeCorrLow[k]<=1.);
			continuum.BandEdgeCorrHi[k] = ((realnum)(RYDLAM/xHi) -
				rfield.anumin(continuum.ipContBandHi[k]-1)) /
				rfield.widflx(continuum.ipContBandHi[k]-1);
			ASSERT( continuum.BandEdgeCorrHi[k]>=0. && continuum.BandEdgeCorrHi[k]<=1.);
			/*fprintf(ioQQQ,"DEBUG bands_continuum %s %.3e %li %li \n",
				continuum.chContBandLabels[k].c_str(),
				continuum.ContBandWavelength[k],
				continuum.ipContBandHi[k] , 
				continuum.ipContBandLow[k]);*/

			if( trace.lgTrace && trace.lgConBug )
			{
				if( k==0 )
					fprintf( ioQQQ, "   ContCreatePointer trace bands\n");
				fprintf( ioQQQ, 
						 "     band %ld label %s low wl= %.3e low ipnt= %li "
						 " hi wl= %.3e hi ipnt= %li \n", 
						 k, 
						 continuum.chContBandLabels[k].c_str(),
						 xLow,
						 continuum.ipContBandLow[k],
						 xHi,
						 continuum.ipContBandHi[k] );
			}
#			if 0
			// hazy table giving band properties
#			include "prt.h"
			fprintf(ioQQQ,
					"DEBUG %s & ", 
					continuum.chContBandLabels[k].c_str() );
			prt_wl( ioQQQ , continuum.ContBandWavelength[k] );
			fprintf(ioQQQ," & ");
			prt_wl( ioQQQ , xHi );
			fprintf(ioQQQ," -- ");
			prt_wl( ioQQQ , xLow );
			fprintf(ioQQQ,"\\\\ \n");
#			endif
			++k;
		}
	}
	/* now validate this incoming data */
	for( i=0; i<continuum.nContBand; ++i )
	{
		/* make sure all are positive */
		if( continuum.ContBandWavelength[i] <=0. )
		{
			fprintf( ioQQQ, " ContBandWavelength band %li has non-positive entry.\n",i );
			cdEXIT(EXIT_FAILURE);
		}
	}

	fclose(ioDATA);
	return;
}
