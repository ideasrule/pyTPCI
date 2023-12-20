/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*OpacityCreateAll compute initial set of opacities for all species */
/*OpacityCreate1Element generate ionic subshell opacities by calling t_ADfA::Inst().phfit */
/*Opacity_iso_photo_cs returns photoionization cross section for isoelectronic sequences */
/*hmiopc derive total H- H minus opacity */
/*rayleh compute Rayleigh scattering cross section for Lya */
/*OpacityValenceRescale routine to rescale non-OP valence shell cross sections */
/******************************************************************************
 *NB NB NB  NB NB NB NB NB NB NB  NB NB NB NB
 * everything set here must be written to the opacity store files
 *
 ****************************************************************************** */
#include "cddefines.h"
#include "dense.h"
#include "continuum.h"
#include "iso.h"
#include "hydrogenic.h"
#include "oxy.h"
#include "trace.h"
#include "heavy.h"
#include "rfield.h"
#include "hmi.h"
#include "atmdat_adfa.h"
#include "save.h"
#include "grains.h"
#include "hydro_bauman.h"
#include "opacity.h"
#include "helike_recom.h"
#include "h2.h"
#include "ipoint.h"
#include "mole.h"
#include "freebound.h"
#include "version.h"
#include "prt.h"

/* initial number of opacity cells available in the opacity stack*/
static const long int ndimOpacityStack = 4200000L;

/*OpacityCreate1Element generate opacities for entire element by calling t_ADfA::Inst().phfit */
STATIC void OpacityCreate1Element(long int nelem);

/*hmiopc derive total H- H minus opacity */
STATIC double hmiopc(double freq);

/*rayleh compute Rayleigh scattering cross section for Lya */
STATIC double rayleh(double ener);

/*Opacity_iso_photo_cs returns photoionization cross section for isoelectronic sequences */
STATIC double Opacity_iso_photo_cs( double energy , long ipISO , long nelem , long index );

/*OpacityCreateReilMan generate photoionization cross sections from Reilman and Manson points */
STATIC void OpacityCreateReilMan(long int low, 
  long int ihi, 
  const realnum energ[], 
  const realnum cross[], 
  long int ncr, 
  long int *ipop, 
  const char *chLabl);

/*OpacityCreatePowerLaw generate array of cross sections using a simple power law fit */
STATIC void OpacityCreatePowerLaw(
	/* lower energy limit on continuum mesh */
	long int ilo, 
	/* upper energy limit on continuum mesh */
	long int ihi, 
	/* threshold cross section */
	double cross, 
	/* power law index */
	double s, 
	/* pointer to opacity offset where this starts */
	long int *ip);

/*ofit compute cross sections for all shells of atomic oxygen */
STATIC void ofit(double e, 
	  realnum opart[]);

/*OpacityValenceRescale routine to rescale non-OP valence shell cross sections for atom */
STATIC void OpacityValenceRescale(
	/* element number on C scale */
	long int nelem ,
	/* scale factor, must be >= 0. */
	double scale )
{
	DEBUG_ENTRY( "OpacityValenceRescale()" );

	/* return if element is not turned on 
	 * >>chng 05 oct 19, this had not been done, so low in the opacity offset below was
	 * not set, and opacity index was negative - only problem when K turned off */
	if( !dense.lgElmtOn[nelem] )
	{
		return;
	}

	ASSERT( scale >= 0. );

	long ion = 0;
	/* this is valence shell on C scale */
	long nshell = Heavy.nsShells[nelem][ion] - 1;

	/* set lower and upper limits to this range */
	long low = opac.ipElement[nelem][ion][nshell][0];
	long ihi = opac.ipElement[nelem][ion][nshell][1];
	long ipop = opac.ipElement[nelem][ion][nshell][2];

	/* loop over energy range of this shell */
	for( long ip=low-1; ip < ihi; ip++ )
	{
		opac.OpacStack[ip-low+ipop] *= scale;
	}
	return;
}

void OpacityCreateAll()
{
	long int i, 
	  ipISO ,
	  need ,
	  nelem;

	realnum opart[7];

	double crs, 
	  dx,
	  eps, 
	  thres, 
	  x;

	DEBUG_ENTRY( "OpacityCreateAll()" );

	/* make and print dust opacities
	 * fill in dstab and dstsc, totals, zero if no dust
	 * may be different if different grains are turned on */
	GrainsInit();

	/* flag lgOpacAllocated says whether opacity stack has been generated
	 * only do this one time per core load  */
	if( lgOpacAllocated )
	{
		/* this is not the first time code called */
		if( trace.lgTrace )
		{
			fprintf( ioQQQ, " OpacityCreateAll called but NOT evaluated since already done.\n" );
		}
		return;
	}

	/* create the space for the opacity stack */
	opac.OpacStack.reserve(ndimOpacityStack);
	lgOpacAllocated = true;

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " OpacityCreateAll called, evaluating.\n" );
	}

	/* zero out opac since this array sometimes addressed before OpacityAddTotal called */
	for( i=0; i < rfield.nflux_with_check; i++ )
	{
		opac.opacity_abs[i] = 0.;
	}

	/* nOpacTot is number of opacity cells in OpacStack filled so far by opacity generating routines */
	opac.nOpacTot = 0;

	long ipIonPotHeI = ipoint(atmdat.EIonPot[ipHELIUM][0]);
	long ipIonPotHeII = ipoint(atmdat.EIonPot[ipHELIUM][1]);

	/* photoionization of h, he-like iso electronic sequences */
	for( ipISO=ipH_LIKE; ipISO<=ipHE_LIKE; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			iso_sp[ipISO][nelem].HighestLevelOpacStack.resize(0);
			if( dense.lgElmtOn[nelem] )
			{
				long int nupper;

				/* this is the opacity offset in the general purpose pointer array */
				/* indices are type, shell. ion, element, so this is the inner shell,
				 * NB - this works for H and He, but the second index should be 1 for Li */
				opac.ipElement[nelem][nelem-ipISO][0][2] = opac.nOpacTot + 1;
	
				fixit("opacities really need to be owned by states or species");
				// and stored in STL containers so we don't have to mess
				// with remembering what the upper and lower limits are

				// all iso states go to high-energy limit of code 
				nupper = rfield.nflux_with_check;
				for( long index=0; index < iso_sp[ipISO][nelem].numLevels_max; index++ )
				{
					/* this is array index to the opacity offset */
					iso_sp[ipISO][nelem].fb[index].ipOpac = opac.nOpacTot + 1;

					/* first make sure that first energy point is at least near the limit */
					long ipThresh = iso_sp[ipISO][nelem].fb[index].ipIsoLevNIonCon-1;
					ASSERT( rfield.anumin(ipThresh) <= iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd &&
							rfield.anumax(ipThresh) >= iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd );

					/* number of cells we will need to do this level */
					need = nupper - ipThresh;
					ASSERT( need > 0 );

					for( i=ipThresh; i < nupper; i++ )
					{
						double crs = Opacity_iso_photo_cs( rfield.anu(i) , ipISO , nelem , index );
						opac.OpacStack.emplace_back( crs ); 
						if( index==iso_sp[ipISO][nelem].numLevels_max-1 )
							iso_sp[ipISO][nelem].HighestLevelOpacStack.emplace_back( crs );
					}

					opac.nOpacTot += need;
				}
			}
		}
	}

	/* H2 continuum dissociation opacity */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
	{
		if( (*diatom)->lgEnabled && mole_global.lgStancil )
		{
			for( vector< diss_tran >::iterator tran = (*diatom)->Diss_Trans.begin(); tran != (*diatom)->Diss_Trans.end(); ++tran )
			{
				/* choose to integrate from 0.1 to 4 Ryd, data only extends from 0.7 to ~2 Ryd */
				long lower_limit = ipoint(tran->energies[0]);
 				long upper_limit = ipoint(tran->energies.back());
 				upper_limit = MIN2( upper_limit, rfield.nflux-1 );
 				long num_points = 0;
 
				for(i = lower_limit; i <= upper_limit; ++i)
				{
					opac.OpacStack.emplace_back( MolDissocCrossSection(*tran, rfield.anu(i)) );
					++num_points;
				}
				opac.nOpacTot += num_points;
			}
		}
	}
	
	/* Lyman alpha damping wings - Rayleigh scattering */
	opac.ipRayScat = opac.nOpacTot + 1;
	for( i=0; i < iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon; i++ )
	{
		opac.OpacStack.emplace_back( rayleh(rfield.anu(i)) );
	}
	opac.nOpacTot += iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon;

	/* ==============================================================
	 * this block of code defines the electron scattering cross section
	 * for all energies */

	/* assume Thomson scattering up to ipCKshell, 20.6 Ryd=0.3 keV */
	opac.iopcom = opac.nOpacTot + 1;
	for( i=0; i < opac.ipCKshell; i++ )
	{
		opac.OpacStack.emplace_back( SIGMA_THOMSON );
	}

	/* Klein-Nishina from eqn 7.5, 
	 * >>refer	Klein-Nishina	cs	Rybicki and Lightman */
	for( i=opac.ipCKshell; i < rfield.nflux_with_check; i++ )
	{
		dx = rfield.anu(i)/3.7573e4;

		opac.OpacStack.emplace_back( SIGMA_THOMSON*3.e0/4.e0*((1.e0 + 
		  dx)/POW3(dx)*(2.e0*dx*(1.e0 + dx)/(1.e0 + 2.e0*dx) - log(1.e0+
		  2.e0*dx)) + 1.e0/2.e0/dx*log(1.e0+2.e0*dx) - (1.e0 + 3.e0*
		  dx)/POW3(1.e0 + 2.e0*dx)) );
	}
	opac.nOpacTot += rfield.nflux_with_check - 1 + 1;

	/* ============================================================== */

	/* pair production */
	opac.ioppr = opac.nOpacTot + 1;
	for( i=opac.ippr-1; i < rfield.nflux_with_check; i++ )
	{
		/* pair production heating rate for unscreened H + He
		 * fit to figure 41 of Experimental Nuclear Physics,
		 * Vol1, E.Segre, ed */

		x = rfield.anu(i)/7.512e4*2.;

		opac.OpacStack.emplace_back( 5.793e-28*
		  POW2((-0.46737118 + x*(0.349255416 + x*0.002179893))/(1. + 
		  x*(0.130471301 + x*0.000524906))) );
	}
	opac.nOpacTot += rfield.nflux_with_check - opac.ippr + 1;

	/* brems (free-free) opacity */
	opac.ipBrems = opac.nOpacTot + 1;

	for( i=0; i < rfield.nflux_with_check; i++ )
	{
		/* Inflate precomputed opacity value by 30 dex to avoid underflows */
		/* free free opacity needs g(ff)*(1-exp(hn/kT))/SQRT(T)*1E-30 */
		opac.OpacStack.emplace_back( FREE_FREE_ABS * 1e30 / POW3(rfield.anu(i)) );
	}
	opac.nOpacTot += rfield.nflux_with_check - 1 + 1;

	opac.iphmra = opac.nOpacTot + 1;
	for( i=0; i < rfield.nflux_with_check; i++ )
	{
		/* following is ratio of h minus to neut h bremss opacity */
		opac.OpacStack.emplace_back( 0.1175*rfield.anusqrt(i) );
	}
	opac.nOpacTot += rfield.nflux_with_check - 1 + 1;

	opac.iphmop = opac.nOpacTot + 1;
	for( i=hmi.iphmin-1; i < ipIonPotHeI; i++ )
	{
		/* H- hminus H minus bound-free opacity */
		opac.OpacStack.emplace_back( hmiopc(rfield.anu(i)) );
	}
	opac.nOpacTot += ipIonPotHeI - hmi.iphmin + 1;

	/* ============================================================== */

	/* This check will get us through "H2 photoionization cross section" below.	*/
	/* >>chng 07 oct 10, by Ryan.  Added this check for allotted memory.	*/
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
	{
		(*diatom)->ip_photo_opac_offset = opac.nOpacTot + 1;
		opac.nOpacTot += (*diatom)->OpacityCreate( opac.OpacStack );
	}

	/* H2+ H2P h2plus photoabsorption */
	{
		/* fits to cross section for photo dist of H_2^+ */
		const long nCSH2P = 5;
		static const realnum enh2p[nCSH2P]={6.75f,8.68f,10.54f,12.46f,14.28f};
		static const realnum csh2p[nCSH2P]={0.24f, 2.5f,  7.1f,  6.0f,  2.7f};
		/* >>refer	H2+	photodissoc	Buckingham, R.A., Reid, S., & Spence, R. 1952, MNRAS 112, 382, 0 K temp */
		OpacityCreateReilMan(opac.ih2pnt[0],opac.ih2pnt[1],enh2p,csh2p,nCSH2P,&opac.ih2pof, "H2+ ");

		// Now do an excited superstate of H2+
		const long nCSH2P_ex = 9;
		// These opacities are a rough approximation to figure 4 of the following reference, for v=9 of H2+:
		/* >>refer   H2+     photodissoc     Dunn, G. H. 1968, Phys.Rev. 172, 1 */
		static const realnum enh2p_ex[nCSH2P_ex]={0.69f,0.83f,0.95f,1.03f,1.24f,1.38f,1.77f,2.48f,14.28f};
		static const realnum csh2p_ex[nCSH2P_ex]={1e-5f,1e-4f,0.01f,0.08f, 2.0f,10.0f,20.0f, 8.0f,  1.0f};
		OpacityCreateReilMan(opac.ih2pnt_ex[0],opac.ih2pnt_ex[1],enh2p_ex,csh2p_ex,nCSH2P_ex,&opac.ih2pof_ex, "H2+*");
	}

	if( dense.lgElmtOn[ipHELIUM] )
	{
		/* HeI singlets neutral helium ground */
		opac.iophe1 = opac.nOpacTot + 1;
		opac.ipElement[ipHELIUM][0][0][2] = opac.iophe1;
		for( i=iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon-1; i < rfield.nflux_with_check; i++ )
		{
			crs = t_ADfA::Inst().phfit(2,2,1,rfield.anu(i)*EVRYD);
			opac.OpacStack.emplace_back( crs*1e-18 );
		}
		opac.nOpacTot += rfield.nflux_with_check - iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon + 1;
	}

	/* these are opacity offset points that would be defined in OpacityCreate1Element,
	 * but this routine will not be called for H and He
	 * generate all heavy element opacities, everything heavier than He,
	 * nelem is on the C scale, so Li is 2 */
	/*>>chng 99 jan 27, do not reevaluate hydrogenic opacity below */
	for( nelem=2; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			OpacityCreate1Element(nelem);
		}
	}

	/* option to rescale atoms of some elements that were not done by opacity project
	 * the valence shell - two arguments - element number on C scale, and scale factor */
	/*>>chng 05 sep 26, fudge factor to get atomic K fraction along well defined line of sight
	 * to be observed value - this is ratio of cross sections, actual value is very uncertain since
	 * differences betweeen Verner & opacity project are huge */
	OpacityValenceRescale( ipPOTASSIUM , 5. );

	/* now add on some special cases, where exicted states, etc, come in */
	/* Nitrogen
	 * >>refer	n1	photo	Henry, R., ApJ 161, 1153.
	 * photoionization of excited level of N+ */
	OpacityCreatePowerLaw(opac.in1[0],opac.in1[1],9e-18,1.75,&opac.in1[2]);

	/* atomic Oxygen
	 * only do this if 1996 Verner results are used */
	if( dense.lgElmtOn[ipOXYGEN] && t_ADfA::Inst().get_version() == PHFIT96 )
	{
		/* integrate over energy range of the valence shell of atomic oxygen*/
		for( i=opac.ipElement[ipOXYGEN][0][2][0]-1; i < opac.ipElement[ipOXYGEN][0][2][1]; i++ )
		{
			/* call special routine to evaluate partial cross section for OI shells */
			eps = rfield.anu(i)*EVRYD;
			ofit(eps,opart);

			/* this will be total cs of all processes leaving shell 3 */
			crs = opart[0];
			for( long n=1; n < 6; n++ )
			{
				/* add up table of cross sections */
				crs += opart[n];
			}
			/* convert to cgs and overwrite cross sections set by OpacityCreate1Element */
			crs *= 1e-18;
			/* this should NOT use emplace_back() since it overwrites already allocated elements */
			opac.OpacStack[i-opac.ipElement[ipOXYGEN][0][2][0]+opac.ipElement[ipOXYGEN][0][2][2]] = crs;
		}
	}

	/* Henry nubmers for 1S excit state of OI, OP data very sparse */
	OpacityCreatePowerLaw(opac.ipo1exc[0],opac.ipo1exc[1],4.64e-18,0.,&opac.ipo1exc[2]);

	/* photoionization of excited level of O2+ 1D making 5007
	 * fit to TopBase Opacity Project cs */
	OpacityCreatePowerLaw(opac.ipo3exc[0],opac.ipo3exc[1],3.8e-18,0.,&opac.ipo3exc[2]);

	/* photoionization of excited level of O2+ 1S making 4363 */
	OpacityCreatePowerLaw(opac.ipo3exc3[0],opac.ipo3exc3[1],5.5e-18,0.01,
	  &opac.ipo3exc3[2]);

	/* photoionization to excited states of O+ */
	opac.iopo2d = opac.nOpacTot + 1;
	thres = rfield.anu(oxy.i2d-1);
	for( i=oxy.i2d-1; i < ipIonPotHeII; i++ )
	{
		crs = 3.85e-18*(4.4*powpq(rfield.anu(i)/thres,-3,2) - 3.38*
				powpq(rfield.anu(i)/thres,-5,2));

		opac.OpacStack.emplace_back( crs );
	}
	opac.nOpacTot += ipIonPotHeII - oxy.i2d + 1;

	/* magnesium
	 * photoionization of excited level of Mg+
	 * fit to opacity project data Dima got */
	opac.ipOpMgEx = opac.nOpacTot + 1;
	for( i=opac.ipmgex-1; i < iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon; i++ )
	{
		opac.OpacStack.emplace_back( 
			(0.2602325880970085 + 
			 445.8558249365131*exp(-rfield.anu(i)/0.1009243952792674))*
			1e-18 );
	}
	opac.nOpacTot += iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon - opac.ipmgex + 1;

	/* Calcium
	 * excited states of Ca+ */
	OpacityCreatePowerLaw(opac.ica2ex[0],opac.ica2ex[1],4e-18,1.,&opac.ica2op);

	ASSERT( size_t(opac.nOpacTot) == opac.OpacStack.size() );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			" OpacityCreateAll return OK, number of opacity cells used in OpacStack = %ld\n", 
		  opac.nOpacTot );
	}

	/* option to compile opacities into file for later use 
	 * this is executed if the 'compile opacities' command is entered */
	if( opac.lgCompileOpac )
	{
		fprintf( ioQQQ, "The COMPILE OPACITIES command is currently not supported\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
/*OpacityCreatePowerLaw generate array of cross sections using a simple power law fit */
STATIC void OpacityCreatePowerLaw(
	/* lower energy limit on continuum mesh */
	long int ilo, 
	/* upper energy limit on continuum mesh */
	long int ihi, 
	/* threshold cross section */
	double cross, 
	/* power law index */
	double s, 
	/* pointer to opacity offset where this starts */
	long int *ip)
{
	long int i;
	double thres;

	DEBUG_ENTRY( "OpacityCreatePowerLaw()" );

	/* non-positive cross section is unphysical */
	ASSERT( cross > 0. );

	/* place in the opacity stack where we will stuff cross sections */
	*ip = opac.nOpacTot + 1;
	ASSERT( *ip > 0 );
	ASSERT( ilo > 0 );
	thres = rfield.anu(ilo-1);

	for( i=ilo-1; i < ihi; i++ )
	{
		opac.OpacStack.emplace_back( cross*pow(rfield.anu(i)/thres,-s) );
	}

	opac.nOpacTot += ihi - ilo + 1;
	return;
}

/*OpacityCreateReilMan generate photoionization cross sections from Reilman and Manson points */
STATIC void OpacityCreateReilMan(long int low, 
  long int ihi, 
  const realnum energ[], 
  const realnum cross[], 
  long int ncr, 
  long int *ipop, 
  const char *chLabl)
{
	long int i, 
	  ics, 
	  j; 

	const int NOP = 100;
	realnum cs[NOP], 
	  en[NOP], 
	  slope;

	DEBUG_ENTRY( "OpacityCreateReilMan()" );

	/* this is the opacity entering routine designed for
	 * the Reilman and Manson tables.  It works with incident
	 * photon energy (entered in eV) and cross sections in megabarns
	 * */
	*ipop = opac.nOpacTot + 1;
	ASSERT( *ipop > 0 );

	if( ncr > NOP )
	{
		fprintf( ioQQQ, " Too many opacities were entered into OpacityCreateReilMan.  Increase the value of NOP.\n" );
		fprintf( ioQQQ, " chLabl was %4.4s\n", chLabl );
		cdEXIT(EXIT_FAILURE);
	}
	if( ncr < 2 )
	{
		fprintf( ioQQQ, " Too few opacities were entered into OpacityCreateReilMan.\n" );
		fprintf( ioQQQ, " chLabl was %4.4s\n", chLabl );
		cdEXIT(EXIT_FAILURE);
	}

	/* the array CROSS has ordered pairs of elements.
	 * the first is the energy in eV (not Ryd)
	 * and the second is the cross section in megabarns */
	for( i=0; i < ncr; i++ )
	{
		en[i] = energ[i]/13.6f;
		cs[i] = cross[i]*1e-18f;
	}

	ASSERT( low>0 );
	if( en[0] > rfield.anu(low-1) )
	{
		fprintf( ioQQQ, 
			" OpacityCreateReilMan: The entered opacity energy bandwidth is not large enough (low fail).\n" );
		fprintf( ioQQQ, 
			" The desired energy (Ryd) was%12.5eeV and the lowest entered in the array was%12.5e eV\n", 
		  rfield.anu(low-1)*EVRYD, en[0]*EVRYD );

		fprintf( ioQQQ, " chLabl was %4.4s\n", chLabl );
		fprintf( ioQQQ, " The original energy (eV) and cross section (mb) arrays follow:\n" );
		fprintf( ioQQQ, " " );

		for( i=0; i < ncr; i++ )
		{
			fprintf( ioQQQ, "%11.4e", energ[i] );
			fprintf( ioQQQ, "%11.4e", cross[i] );
		}

		fprintf( ioQQQ, "\n" );
		cdEXIT(EXIT_FAILURE);
	}

	slope = (cs[1] - cs[0])/(en[1] - en[0]);
	ics = 1;

	/* now fill in the opacities using linear interpolation */
	for( i=low-1; i < ihi; i++ )
	{
		if( rfield.anu(i) > en[ics-1] && rfield.anu(i) <= en[ics] )
		{
			opac.OpacStack.emplace_back( cs[ics-1] + slope*(rfield.anu(i) - en[ics-1]) );
		}

		else
		{
			ics += 1;
			if( ics + 1 > ncr )
			{
				fprintf( ioQQQ, " OpacityCreateReilMan: The entered opacity energy bandwidth is not large enough (high fail).\n" );
				fprintf( ioQQQ, " The entered energy was %10.2eeV and the highest in the array was %10.2eeV\n", 
				  rfield.anu(i)*13.6, en[ncr-1]*13.6 );
				fprintf( ioQQQ, " chLabl was %4.4s\n", chLabl
				   );
				fprintf( ioQQQ, " The lowest energy enterd in the array was%10.2e eV\n", 
				  en[0]*13.65 );
				fprintf( ioQQQ, " The highest energy ever needed would be%10.2eeV\n", 
				  rfield.anu(ihi-1)*13.6 );
				fprintf( ioQQQ, " The lowest energy needed was%10.2eeV\n", 
				  rfield.anu(low-1)*13.6 );
				cdEXIT(EXIT_FAILURE);
			}

			slope = (cs[ics] - cs[ics-1])/(en[ics] - en[ics-1]);
			if( rfield.anu(i) > en[ics-1] && rfield.anu(i) <= en[ics] )
			{
				opac.OpacStack.emplace_back( cs[ics-1] + slope*(rfield.anu(i) - en[ics-1]) );
			}
			else
			{
				ASSERT( i > 0);
				fprintf( ioQQQ, " Internal logical error in OpacityCreateReilMan.\n" );
				fprintf( ioQQQ, " The desired energy (%10.2eeV), I=%5ld, is not within the next energy bound%10.2e%10.2e\n", 
				  rfield.anu(i)*13.6, i, en[ics-1], en[ics] );

				fprintf( ioQQQ, " The previous energy (eV) was%10.2e\n", 
				  rfield.anu(i-1)*13.6 );

				fprintf( ioQQQ, " Here comes the energy array.  ICS=%4ld\n", 
				  ics );

				for( j=0; j < ncr; j++ )
				{
					fprintf( ioQQQ, "%10.2e", en[j] );
				}
				fprintf( ioQQQ, "\n" );

				fprintf( ioQQQ, " chLabl was %4.4s\n", chLabl );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	/* >>chng 02 may 09, this was a significant logcal error */
	/* >>chng 02 may 08, by Ryan.  This routine did not update the total slots filled.	*/
	opac.nOpacTot += ihi - low + 1;
	return;
}


/*ofit compute cross sections for all shells of atomic oxygen */
STATIC void ofit(double e, 
	  realnum opart[])
{
	static const double y[7][5] = {
		{8.915,3995.,3.242,10.44,0.0},
		{11.31,1498.,5.27,7.319,0.0},
		{10.5,1.059e05,1.263,13.04,0.0},
		{19.49,48.47,8.806,5.983,0.0},
		{50.,4.244e04,0.1913,7.012,4.454e-02},
		{110.5,0.1588,148.3,-3.38,3.589e-02},
		{177.4,32.37,381.2,1.083,0.0}
	};
	static const double eth[7]={13.62,16.94,18.79,28.48,50.,110.5,538.};
	static const long l[7]={1,1,1,0,1,1,0};

	DEBUG_ENTRY( "ofit()" );
	/*compute cross sections for all shells of atomic oxygen
	 * Photoionization of OI
	 * Input parameter:   e - photon energy, eV
	 * Output parameters: otot - total photoionization cross section, Mb
	 *  opart(1) - 2p-shell photoionization, goes to 4So
	 *  opart(2) - 2p-shell photoionization, goes to 2Do
	 *  opart(3) - 2p-shell photoionization, goes to 2Po
	 *  opart(4) - 2s-shell photoionization
	 *  opart(5) - double photoionization, goes to O++
	 *  opart(6) - triple photoionization, goes to O+++
	 *  opart(7) - 1s-shell photoionization */

	for( int i=0; i < 7; i++ )
	{
		opart[i] = 0.0;
	}

	for( int i=0; i < 7; i++ )
	{
		if( e >= eth[i] )
		{
			// this assert is trivially true, but it helps PGCC
			ASSERT( i < 7 );
			double q = 5.5 - 0.5*y[i][3] + l[i];

			double x = e/y[i][0];

			opart[i] = (realnum)(y[i][1]*(POW2(x - 1.0) + POW2(y[i][4]))/
			  pow(x,q)/pow(1.0 + sqrt(x/y[i][2]),y[i][3]));

		}
	}
	return;
}

/******************************************************************************/

/*OpacityCreate1Element generate ionic subshell opacities by calling t_ADfA::Inst().phfit */
STATIC void OpacityCreate1Element(
		  /* atomic number on the C scale, lowest ever called will be Li=2 */
		  long int nelem)
{
	long int ihi, 
	  ip, 
	  ipop, 
	  low, 
	  nelec, 
	  ion, 
	  nshell;
	double cs; 
	double energy;

	DEBUG_ENTRY( "OpacityCreate1Element()" );

	/* confirm range of validity of atomic number, Li=2 should be the lightest */
	ASSERT( nelem >= 2 );
	ASSERT( nelem < LIMELM );

	/*>>chng 99 jan 27, no longer redo hydrogenic opacity here */
	/* for( ion=0; ion <= nelem; ion++ )*/
	for( ion=0; ion < nelem; ion++ )
	{

		/* will be used for a sanity check on number of hits in a cell*/
		for( ip=0; ip < rfield.nflux_with_check; ip++ )
		{
			opac.opacity_abs[ip] = 0.;
		}

		/* number of bound electrons */
		nelec = nelem+1 - ion;

		/* loop over all shells, from innermost K shell to valence */
		for( nshell=0; nshell < Heavy.nsShells[nelem][ion]; nshell++ )
		{
			/* this is array index for start of this shell within large opacity stack */
			opac.ipElement[nelem][ion][nshell][2] = opac.nOpacTot +  1;

			/* set lower and upper limits to this range */
			low = opac.ipElement[nelem][ion][nshell][0];
			ihi = opac.ipElement[nelem][ion][nshell][1];
			ipop = opac.ipElement[nelem][ion][nshell][2];

			/* make sure indices are within correct bounds,
			 * mainly check on logic for detecting missing shells */
			ASSERT( low <= ihi || low<5 );

			/* loop over energy range of this shell */
			for( ip=low-1; ip < ihi; ip++ )
			{
				/* photo energy MAX so that we never eval below threshold */
				energy = MAX2(rfield.anu(ip)*EVRYD , 
					t_ADfA::Inst().getEthresh(nshell+1,nelec,nelem+1));

				/* the cross section in mega barns */
				cs = t_ADfA::Inst().phfit(nelem+1,nelec,nshell+1,energy);
				/* cannot assert that cs is positive since, at edge of shell,
				 * energy might be slightly below threshold and hence zero,
				 * due to finite size of continuum bins */
				opac.OpacStack.emplace_back( cs*1e-18 );

				/* add this to total opacity, which we will confirm to be greater than zero below */
				opac.opacity_abs[ip] += cs;
			}

			opac.nOpacTot += ihi - low + 1;

			/* save pointers option */
			if( save.lgPunPoint )
			{
				fprintf( save.ipPoint, "%3ld%3ld%3ld%10.2e%10.2e%10.2e%10.2e\n", 
				  nelem, ion, nshell, rfield.anu(low-1), rfield.anu(ihi-1), 
				  opac.OpacStack[ipop-1], opac.OpacStack[ihi-low+ipop-1] );
			}
		}

		ASSERT( Heavy.nsShells[nelem][ion] >= 1 );
		/*confirm that total opacity is greater than zero  */
		for( ip=opac.ipElement[nelem][ion][Heavy.nsShells[nelem][ion]-1][0]-1; 
			ip < continuum.KshellLimit; ip++ )
		{
			ASSERT( opac.opacity_abs[ip] > 0. );
		}

	}
	return;
}

/*Opacity_iso_photo_cs returns photoionization cross section for isoelectronic sequences */
STATIC double Opacity_iso_photo_cs( 
		/* photon energy ryd */
		double EgammaRyd , 
		/* iso sequence */
		long ipISO , 
		/* charge, 0 for H */
		long nelem , 
		/* index */
		long index )
{
	double crs=-DBL_MAX;

	DEBUG_ENTRY( "Opacity_iso_photo_cs()" );

	if( ipISO==ipH_LIKE )
	{
		if( index==0 )
		{
			double EgammaEV = MAX2(EgammaRyd, atmdat.getIonPot(nelem, nelem)) * EVRYD;
			crs = t_ADfA::Inst().phfit(nelem+1,1,1,EgammaEV)* 1e-18;
			/* make sure cross section is reasonable */
			ASSERT( crs > 0. && crs < 1e-10 );
		}
		else if( index < iso_sp[ipISO][nelem].numLevels_max - iso_sp[ipISO][nelem].nCollapsed_max )
		{
			/* photon energy relative to threshold */
			double photon = MAX2( EgammaRyd/iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd, 1. + FLT_EPSILON*2. );

			crs = H_photo_cs( photon , N_(index), L_(index), nelem+1 );
			/* make sure cross section is reasonable */
			ASSERT( crs > 0. && crs < 1e-10 );
		}
		else if( N_(index) <= NHYDRO_MAX_LEVEL )
		{
			/* for first cell, depending on the current resolution of the energy mesh,
			 * the center of the first cell can be below the ionization limit of the
			 * level.  do not let the energy fall below this limit */
			/* This will make sure that we don't call epsilon below threshold,
			 * the factor 1.001 was chosen so that t_ADfA::Inst().hpfit, which works
			 * in terms of Dima's Rydberg constant, is not tripped below threshold */
			EgammaRyd = MAX2( EgammaRyd , iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd*1.001f );

			crs = t_ADfA::Inst().hpfit(nelem+1,N_(index),EgammaRyd*EVRYD);
			/* make sure cross section is reasonable */
			ASSERT( crs > 0. && crs < 1e-10 );
		}
		else
		{
			/* photon energy relative to threshold */
			double photon = MAX2( EgammaRyd/iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd, 1. + FLT_EPSILON*2. );

			/* cross section for collapsed level should be 
			 * roughly equal to cross-section for yrast level,
			 * so third parameter is n - 1. */
			crs = H_photo_cs( photon , N_(index), N_(index)-1, nelem+1 );

			/* make sure cross section is reasonable */
			ASSERT( crs > 0. && crs < 1e-10 );
		}
	}
	else if( ipISO==ipHE_LIKE )
	{
		EgammaRyd = MAX2( EgammaRyd , iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd);
		/* this would be a collapsed level */
		if( index >= iso_sp[ipHE_LIKE][nelem].numLevels_max - iso_sp[ipHE_LIKE][nelem].nCollapsed_max )
		{
			long int nup = iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max + index + 1 -
				(iso_sp[ipHE_LIKE][nelem].numLevels_max - iso_sp[ipHE_LIKE][nelem].nCollapsed_max);

			/* this is a collapsed level - this is hydrogenic routine and
			 * first he-like energy may not agree exactly with threshold for H */
			crs = t_ADfA::Inst().hpfit(nelem,nup ,EgammaRyd*EVRYD);
			/* make sure cross section is reasonable if away from threshold */
			ASSERT( 
				(EgammaRyd < iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd*1.02) ||
				(crs > 0. && crs < 1e-10) );
		}
		else
		{
			long n = N_(index);
			long l = L_(index);
			long S = S_(index);
			/* He_cross_section returns cross section (cm^-2), 
			 * given EgammaRyd, the photon energy in Ryd,
			 * quantum numbers n, l, and S,
			 * nelem is charge, equal to 1 for Helium,
			 * this is a wrapper for cross_section */
			crs = He_cross_section( EgammaRyd, iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd, n, l, S, nelem );

			/* make sure cross section is reasonable */
			ASSERT( crs > 0. && crs < 1e-10 );
		}
	}
	else
		TotalInsanity();
	return(crs);
}

/*hmiopc derive total H- H minus opacity */
static const int NCRS = 33;

STATIC double hmiopc(double freq)
{
	double energy, 
	  hmiopc_v, 
	  x, 
	  y;
	static double y2[NCRS];
	static double crs[NCRS]={0.,0.124,0.398,0.708,1.054,1.437,1.805,
	  2.176,2.518,2.842,3.126,3.377,3.580,3.741,3.851,3.913,3.925,
	  3.887,3.805,3.676,3.511,3.306,3.071,2.810,2.523,2.219,1.898,
	  1.567,1.233,.912,.629,.39,.19};
	static double ener[NCRS]={0.,0.001459,0.003296,0.005256,0.007351,
	  0.009595,0.01201,0.01460,0.01741,0.02044,0.02375,0.02735,0.03129,
	  0.03563,0.04043,0.04576,0.05171,0.05841,0.06601,0.07469,0.08470,
	  0.09638,0.1102,0.1268,0.1470,0.1723,0.2049,0.2483,0.3090,0.4001,
	  0.5520,0.8557,1.7669};
	static bool lgFirst = true;

	DEBUG_ENTRY( "hmiopc()" );

	/* bound free cross section (x10**-17 cm^2) from Doughty et al
	 * 1966, MNRAS 132, 255; good agreement with Wishart MNRAS 187, 59p. */

	/* photoelectron energy, add HMINUSIONPOT to get incoming energy (Ryd) */


	if( lgFirst )
	{
		/* set up coefficients for spline */
		spline(ener,crs,NCRS,2e31,2e31,y2);
		lgFirst = false;
	}

	energy = freq - HMINUSIONPOT;
	if( energy < ener[0] || energy > ener[NCRS-1] )
	{
		hmiopc_v = 0.;
	}
	else
	{
		x = energy;
		splint(ener,crs,y2,NCRS,x,&y);
		hmiopc_v = y*1e-17;
	}
	return( hmiopc_v );
}

/*rayleh compute Rayleigh scattering cross section for Lya */
STATIC double rayleh(double ener)
{
	double rayleh_v;

	DEBUG_ENTRY( "rayleh()" );

	/** \todo	2	update to astro-ph/0308073, Lee, H-W, ApJ in press */
	/* do hydrogen Rayleigh scattering cross sections;
	 * fits to 
	 *>>refer	Ly	scattering	Gavrila, M., 1967, Physical Review 163, 147
	 * and Mihalas radiative damping
	 *
	 * >>chng 96 aug 15, changed logic to do more terms for each part of
	 * rayleigh scattering
	 * if( ener.lt.0.05 ) then
	 *  rayleh = 8.41e-25 * ener**4 * DampOnFac
	 * */
	if( ener < 0.05 )
	{
		rayleh_v = (8.41e-25*powi(ener,4) + 3.37e-24*powi(ener,6))*
		  hydro.DampOnFac;
	}

	else if( ener < 0.646 )
	{
		rayleh_v = (8.41e-25*powi(ener,4) + 3.37e-24*powi(ener,6) + 
		  4.71e-22*powi(ener,14))*hydro.DampOnFac;
	}

	else if( ener >= 0.646 && ener < 1.0 )
	{
		rayleh_v = fabs(0.74959-ener);
		rayleh_v = 1.788e5/POW2(FR1RYD*MAX2(0.001,rayleh_v));
		/*  typical energy between Ly-a and Ly-beta */
		rayleh_v = MAX2(rayleh_v,1e-24)*hydro.DampOnFac;
	}

	else
	{
		rayleh_v = 0.;
	}
	return( rayleh_v );
}
