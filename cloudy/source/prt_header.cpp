/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtHeader print program's header, including luminosities and ionization parameters */
#include "cddefines.h"
#include "phycon.h"
#include "iso.h"
#include "rfield.h"
#include "radius.h"
#include "called.h"
#include "thermal.h"
#include "dense.h"
#include "continuum.h"
#include "ipoint.h"
#include "prt.h"
#include "freebound.h"

void PrtHeader(void)
{
	long int i, 
	  ip2500, 
	  ip2kev;
	double absbol, 
	  absv, 
	  alfox, 
	  avpow, 
	  bolc, 
	  gpowl, 
	  pballog, 
	  pionl, 
	  qballog, 
	  qgaml, 
	  qheiil, 
	  qhel, 
	  ql, 
	  qxl, 
	  radpwl, 
	  ratio, 
	  solar, 
	  tcomp, 
	  tradio;

	DEBUG_ENTRY( "PrtHeader()" );

	if( !called.lgTalk )
	{ 
		return;
	}

	if( prt.lgPrtCitations )
	{
		/* the print citations command */
		CloudyPrintReference();
		fprintf( ioQQQ, "\n" );
		DatabasePrintReference();
		fprintf( ioQQQ, "\n\n" );
	}

	fprintf( ioQQQ, "           %4ldCellPeak",rfield.nflux );
	PrintE82(ioQQQ, rfield.anu(prt.ipeak-1) );
	fprintf( ioQQQ, "   Lo");
	fprintf( ioQQQ,PrintEfmt("%9.2e", rfield.anumin(0) ));
	fprintf( ioQQQ, "=%7.2fm   Hi-Con:", 1e-10*RYDLAM/rfield.anumin(0) );
	PrintE82(ioQQQ,rfield.anumax(rfield.nflux-1));
	fprintf(ioQQQ," Ryd   E(hi):");
	PrintE82(ioQQQ,rfield.egamry());
	fprintf(ioQQQ,"Ryd     E(hi):  %9.2f MeV\n", rfield.egamry()*1e-6*EVRYD );

	if( prt.xpow > 0. )
	{
		prt.xpow = (realnum)(log10(prt.xpow) + radius.pirsq);
		qxl = log10(prt.qx) + radius.pirsq;
	}
	else
	{
		prt.xpow = 0.;
		qxl = 0.;
	}

	if( prt.powion > 0. )
	{
		pionl = log10(prt.powion) + radius.pirsq;
		avpow = prt.powion/rfield.qhtot/EN1RYD;
	}
	else
	{
		pionl = 0.;
		avpow = 0.;
	}

	/* >>chng 97 mar 18, break these two out here, so that returns zero
	 * when no radiation - had been done in the print statement so pirsq was printed */
	if( prt.pbal > 0. )
	{
		pballog = log10(MAX2(1e-30,prt.pbal)) + radius.pirsq;
		qballog = log10(MAX2(1e-30,rfield.qbal)) + radius.pirsq;
	}
	else
	{
		pballog = 0.;
		qballog = 0.;
	}

	if( radius.pirsq > 0. )
	{
		fprintf( ioQQQ, "           L(nu>1ryd):%9.4f   Average nu:",pionl);
		PrintE93(ioQQQ, avpow);
		fprintf( ioQQQ,"   L( X-ray):%9.4f   L(BalC):%9.4f     Q(Balmer C):%9.4f\n", 
		  prt.xpow, pballog, qballog );
		if( pionl > 47. )
		{
			fprintf(ioQQQ,"\n\n WARNING - the continuum has a luminosity %.2e times greater than the sun.\n",
				exp10(  pionl-log10(SOLAR_LUMINOSITY) ) );
			fprintf(ioQQQ," WARNING - Is this correct?  Check the luminosity commands.\n\n\n");
		}
	}
	else
	{
		fprintf( ioQQQ, "           I(nu>1ryd):%9.4f   Average nu:",pionl);
		PrintE93(ioQQQ, avpow);
		fprintf( ioQQQ,"   I( X-ray):%9.4f   I(BalC):%9.4f     Phi(BalmrC):%9.4f\n", 
		  prt.xpow, 
		  log10(MAX2(1e-30,prt.pbal)), 
		  log10(MAX2(1e-30,rfield.qbal)) );
	}

	if( rfield.qhe > 0. )
	{
		qhel = log10(rfield.qhe) + radius.pirsq;
	}
	else
	{
		qhel = 0.;
	}

	if( rfield.qheii > 0. )
	{
		qheiil = log10(rfield.qheii) + radius.pirsq;
	}
	else
	{
		qheiil = 0.;
	}

	if( prt.q > 0. )
	{
		ql = log10(prt.q) + radius.pirsq;
	}
	else
	{
		ql = 0.;
	}

	if( radius.pirsq != 0. )
	{
		fprintf( ioQQQ, 
			"           Q(1.0-1.8):%9.4f   Q(1.8-4.0):%9.4f   Q(4.0-20):"
			"%9.4f   Q(20--):%9.4f     Ion pht flx:", 
		  ql, 
		  qhel, 
		  qheiil, 
		  qxl);
		PrintE93(ioQQQ, rfield.qhtot );
	}
	else
	{
		fprintf( ioQQQ, 
			"           phi(1.0-1.8):%7.4f   phi(1.8-4.0):%7.3f   phi(4.0-20):"
			"%7.3f   phi(20--):%7.3f     Ion pht flx:", 
		  ql, 
		  qhel, 
		  qheiil, 
		  qxl );
		PrintE93(ioQQQ, rfield.qhtot );
	}
	fprintf( ioQQQ,"\n");

	/* estimate alpha (o-x) */
	if( rfield.anu(rfield.nflux-1) > 150. )
	{
		/* the ratio of fluxes is nearly 403.3^alpha ox */
		ip2kev = ipoint(147.);
		ip2500 = ipoint(0.3645);
		if( rfield.flux[0][ip2500-1] > 1e-28 )
		{
			ratio = (rfield.flux[0][ip2kev-1]*rfield.anu(ip2kev-1)/rfield.widflx(ip2kev-1))/
			  (rfield.flux[0][ip2500-1]*rfield.anu(ip2500-1)/rfield.widflx(ip2500-1));
		}
		else
		{
			ratio = 0.;
		}

		if( ratio > 0. )
		{
			alfox = log(ratio)/log(rfield.anu(ip2kev-1)/rfield.anu(ip2500-1));
		}
		else
		{
			alfox = 0.;
		}
	}
	else
	{
		alfox = 0.;
	}

	if( prt.GammaLumin > 0. )
	{
		gpowl = log10(prt.GammaLumin) + radius.pirsq;
		qgaml = log10(prt.qgam) + radius.pirsq;
	}
	else
	{
		gpowl = 0.;
		qgaml = 0.;
	}

	if( prt.pradio > 0. )
	{
		radpwl = log10(prt.pradio) + radius.pirsq;
	}
	else
	{
		radpwl = 0.;
	}

	if( radius.pirsq > 0. )
	{
		fprintf( ioQQQ, 
			"           L(gam ray):%9.4f   Q(gam ray):%9.4f   L(Infred):%9.4f   Alf(ox):%9.4f     Total lumin:%9.4f\n", 
		  gpowl, qgaml, radpwl, alfox, log10(MAX2(SMALLFLOAT,continuum.TotalLumin)) +
		  radius.pirsq );
	}
	else
	{
		fprintf( ioQQQ, 
			"           I(gam ray):%9.4f   phi(gam r):%9.4f   I(Infred):%9.4f   Alf(ox):%9.4f     Total inten:%9.4f\n", 
		  gpowl, qgaml, radpwl, alfox, log10(MAX2(SMALLFLOAT,continuum.TotalLumin)) +
		  radius.pirsq );
	}

	/* magnitudes */
	if( radius.lgPredLumin )
	{
		solar = absbol = -38.;
		if( continuum.TotalLumin> SMALLFLOAT )
		{
			solar = log10(continuum.TotalLumin) + radius.pirsq - log10(SOLAR_LUMINOSITY);
			absbol = 2.5*(log10(MBOL_ZERO_LUMINOSITY/SOLAR_LUMINOSITY) - solar);
		}

		/* absv = 4.79 - 2.5 * (LOG10(MAX(1e-30,continuum.fluxv)) + pirsq - 18.742 -
		 *  1  0.016)
		 * allen 76, page 197 */
		absv = -38.;
		if( continuum.fluxv>SMALLFLOAT )
			absv = -2.5*(log10(MAX2(1e-30,continuum.fluxv)) + radius.pirsq - 20.654202);

		/* >>chng 97 mar 18, following branch so zero returned when no radiation at all */
		if( continuum.fbeta > 0. )
		{
			continuum.fbeta = (realnum)(log10(MAX2(1e-37,continuum.fbeta)) + radius.pirsq);
		}
		else
		{
			continuum.fbeta = 0.;
		}

		bolc = absbol - absv;
		fprintf( ioQQQ, 
			"           log L/Lsun:%9.4f   Abs bol mg:%9.4f   Abs V mag:%9.4f   Bol cor:%9.4f     nuFnu(Bbet):%9.4f\n", 
		  solar, 
		  absbol, 
		  absv, 
		  bolc, 
		  continuum.fbeta );
	}

	rfield.cmcool = 0.;
	rfield.cmheat = 0.;

	for( i=0; i < rfield.nflux; i++ )
	{
		/* CSIGC is Tarter expression times ANU(I)*3.858E-25 */
		rfield.cmcool += (rfield.flux[0][i] + rfield.outlin[0][i] + rfield.outlin_noplot[i] +rfield.ConInterOut[i])*
		  rfield.csigc[i];

		/* Compton heating with correction for induced Compton scattering
		 * CSIGH is Tarter expression times ANU(I)**2 * 3.858E-25 */
		rfield.cmheat += (rfield.flux[0][i] + rfield.outlin[0][i] + rfield.outlin_noplot[i] +rfield.ConInterOut[i])*
		  rfield.csigh[i]*(1. + rfield.OccNumbIncidCont[i]);
	}

	/* all of following need factor of 10^-15 to be true rates */
	rfield.cmheat *= dense.eden*1e-15;
	rfield.cmcool *= dense.eden*4./TE1RYD*1e-15;

	if( rfield.cmcool > 0. )
	{
		rfield.lgComUndr = false;
		tcomp = rfield.cmheat/rfield.cmcool;
	}
	else
	{
		/* underflow if Compt cooling rate is zero */
		rfield.lgComUndr = true;
		tcomp = 0.;
	}

	/* check whether energy density temp is greater than compton temp */
	if( phycon.TEnerDen > (1.05*tcomp) )
	{
		thermal.lgEdnGTcm = true;
	}
	else
	{
		thermal.lgEdnGTcm = false;
	}

	long ipIonPotHeII = ipoint(atmdat.EIonPot[ipHELIUM][1]);

	/* say some ionization parameters */
	fprintf( ioQQQ, "           U(1.0----):");
	PrintE93( ioQQQ, rfield.uh);
	fprintf( ioQQQ, "   U(4.0----):");
	PrintE93( ioQQQ,rfield.uheii );
	fprintf( ioQQQ, "   T(En-Den):");
	PrintE93( ioQQQ,phycon.TEnerDen );
	fprintf( ioQQQ, "   T(Comp):");
	PrintE93( ioQQQ,tcomp );
	fprintf( ioQQQ, "     nuJnu(912A):");
	PrintE93( ioQQQ,prt.fx1ryd );
	fprintf( ioQQQ, "\n");

	/* some occupation numbers */
	fprintf( ioQQQ, "           Occ(FarIR):");
	PrintE93( ioQQQ, rfield.OccNumbIncidCont[0]);
	fprintf( ioQQQ, "   Occ(H n=6):");

	if( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_local + iso_sp[ipH_LIKE][ipHYDROGEN].nCollapsed_local >= 6 )
	{
		long ip6p = iso_sp[ipH_LIKE][ipHYDROGEN].QN2Index(6, 1, 2);
		PrintE93( ioQQQ, rfield.OccNumbIncidCont[iso_sp[ipH_LIKE][ipHYDROGEN].fb[ip6p].ipIsoLevNIonCon-1]);
	}
	else
	{
		PrintE93( ioQQQ, 0.);
	}
	fprintf( ioQQQ, "   Occ(1Ryd):");
	PrintE93( ioQQQ, rfield.OccNumbIncidCont[iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1]);
	fprintf( ioQQQ, "   Occ(4R):");
	PrintE93( ioQQQ, rfield.OccNumbIncidCont[ipIonPotHeII-1]);
	fprintf( ioQQQ, "     Occ (Nu-hi):");
	PrintE93( ioQQQ, rfield.OccNumbIncidCont[rfield.nflux-1] );
	fprintf( ioQQQ, "\n"); 

	/* now print brightness temps */
	tradio = rfield.OccNumbIncidCont[0]*TE1RYD*rfield.anu(0);
	fprintf( ioQQQ, "           Tbr(FarIR):");
	PrintE93( ioQQQ, tradio);
	fprintf( ioQQQ, "   Tbr(H n=6):");
	if( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_local + iso_sp[ipH_LIKE][ipHYDROGEN].nCollapsed_local >= 6 )
	{
		long ip6p = iso_sp[ipH_LIKE][ipHYDROGEN].QN2Index(6, 1, 2);
		PrintE93( ioQQQ, rfield.OccNumbIncidCont[iso_sp[ipH_LIKE][ipHYDROGEN].fb[ip6p].ipIsoLevNIonCon-1]*TE1RYD*rfield.anu(iso_sp[ipH_LIKE][ipHYDROGEN].fb[ip6p].ipIsoLevNIonCon-1));
	}
	else
	{
		PrintE93( ioQQQ, 0.);
	}
	fprintf( ioQQQ, "   Tbr(1Ryd):");
	PrintE93( ioQQQ, rfield.OccNumbIncidCont[iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1]*TE1RYD*rfield.anu(iso_sp[ipH_LIKE][ipHYDROGEN].fb[0].ipIsoLevNIonCon-1));
	fprintf( ioQQQ, "   Tbr(4R):");
	PrintE93( ioQQQ, rfield.OccNumbIncidCont[ipIonPotHeII-1]*TE1RYD*rfield.anu(ipIonPotHeII-1));
	fprintf( ioQQQ, "     Tbr (Nu-hi):");
	PrintE93( ioQQQ, rfield.OccNumbIncidCont[rfield.nflux-1]*TE1RYD*rfield.anu(rfield.nflux-1));
	fprintf( ioQQQ, "\n");

	if( tradio > 1e9 )
	{
		fprintf( ioQQQ, 
			" >>>The radio brightness temperature is very large, %.2eK at %.2ecm.  Is this physical???\n", 
		  tradio, 1e-8*RYDLAM/rfield.anu(0) );
	}

	/* skip a line */
	fprintf( ioQQQ, "  \n" );
	return;
}
