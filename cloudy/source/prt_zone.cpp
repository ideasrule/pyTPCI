/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtZone print out individual zone results, call by iter_end_check at very
 * end of zone calculations */
#include "cddefines.h"
#include "iso.h"
#include "grainvar.h"
#include "pressure.h"
#include "wind.h"
#include "conv.h"
#include "trace.h"
#include "magnetic.h"
#include "called.h"
#include "dynamics.h"
#include "h2.h"
#include "secondaries.h"
#include "opacity.h"
#include "geometry.h"
#include "hmi.h"
#include "thermal.h"
#include "radius.h"
#include "phycon.h"
#include "abund.h"
#include "hydrogenic.h"
#include "ionbal.h"
#include "elementnames.h"
#include "prt.h"
#include "deuterium.h"
#include "mole.h"
#include "rfield.h"
#include "freebound.h"
#include "dense.h"

void PrtZone(void)
{
	char chField7[32];
	char chLet, 
	  chQHMark;
	long int i, 
	  ishift, 
	  nelem ,
	  mol;

	DEBUG_ENTRY( "PrtZone()" );

	if( thermal.lgUnstable )
	{
		chLet = 'u';
	}
	else
	{
		chLet = ' ';
	}

	/* middle of zone for printing
	rmidle = radius.Radius - radius.drad*0.5*radius.dRadSign;
	dmidle = radius.depth - radius.drad*0.5; */

	/* option to print single line when quiet but tracing convergence
	 * with "trace convergence" command */
	if( called.lgTalk || trace.nTrConvg )
	{
		/* print either ####123 or ###1234 */
		if( nzone <= 999 )
		{
			sprintf( chField7, "####%3ld", nzone );
		}
		else
		{
			sprintf( chField7, "###%4ld", nzone );
		}

		fprintf(ioQQQ, " %7.7s %cTe:",chField7, chLet);
		PrintE93(ioQQQ,phycon.te);
		fprintf(ioQQQ," Hden:");
		PrintE93(ioQQQ,dense.gas_phase[ipHYDROGEN]);
		fprintf(ioQQQ," Ne:");
		PrintE93(ioQQQ,dense.eden);
		fprintf(ioQQQ," R:");
		PrintE93(ioQQQ,radius.Radius_mid_zone );
		fprintf(ioQQQ," R-R0:");
		PrintE93(ioQQQ,radius.depth_mid_zone);
		fprintf(ioQQQ," dR:");
		PrintE93(ioQQQ,radius.drad);
		fprintf(ioQQQ," NTR:%3ld Htot:",conv.nPres2Ioniz);
		PrintE93(ioQQQ,thermal.htot);
		fprintf(ioQQQ," T912:");
		fprintf(ioQQQ,PrintEfmt("%9.2e",opac.TauAbsGeo[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1] ));
		fprintf(ioQQQ,"###\n");

		if( trace.nTrConvg )
		{
			fprintf( ioQQQ, " H:%.2e %.2e 2H2/H: %.2e He: %.2e %.2e %.2e\n", 
			  dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN], 
			  dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN], 
			  2.*hmi.H2_total/dense.gas_phase[ipHYDROGEN],
			  dense.xIonDense[ipHELIUM][0]/SDIV(dense.gas_phase[ipHELIUM]), 
			  dense.xIonDense[ipHELIUM][1]/SDIV(dense.gas_phase[ipHELIUM]), 
			  dense.xIonDense[ipHELIUM][2]/SDIV(dense.gas_phase[ipHELIUM])
			  );
		}
	}

	/* now return if not talking */
	if( !called.lgTalk && !trace.nTrConvg )
	{ 
		return;
	}

	/* lgDenFlucOn set to true in zero, only false when variable abundances are on,
	 * lgAbTaON set true when element table used */
	if( !dense.lgDenFlucOn || abund.lgAbTaON )
	{
		fprintf( ioQQQ, " Abun:" );
		for( i=0; i < LIMELM; i++ )
		{
			fprintf( ioQQQ,PrintEfmt("%8.1e", dense.gas_phase[i] ));
		}
		fprintf( ioQQQ, "\n" );
	}

	/*-------------------------------------------------
	 * print wind parameters if windy model */
	if( !wind.lgStatic() )
	{
		double fac;
		/* find denominator for fractional contributions */
		if( wind.AccelTotalOutward == 0. )
			fac = 1.;
		else
			fac = wind.AccelTotalOutward;
		fprintf( ioQQQ, 
			" Dynamics wind V:%.3e km/s a(grav):%.2e a(tot):%.2e Fr(cont):%6.3f "
			"Fr(line):%6.3f \n",
			wind.windv/1e5 ,
			-wind.AccelGravity,
			wind.AccelTotalOutward,
			wind.AccelCont/ fac, 
			wind.AccelLine/fac );

		/* print advection information */
		if( dynamics.lgAdvection || dynamics.lgTimeDependentStatic )
			DynaPrtZone();
	}

	/* print line with radiation pressure if significant */
	if( pressure.pbeta > .05 )
		PrtLinePres(ioQQQ);

	// report fraction of total H in form of H^0 and H^+
	double hatmic = 0.;
	shared_ptr<chem_nuclide> elHydrogen = findnuclide("H");
	for(mol = 0; mol < mole_global.num_calc; mol++) {
		if (mole_global.list[mol]->isIsotopicTotalSpecies() && mole_global.list[mol]->nNuclide.find(elHydrogen) !=
			 mole_global.list[mol]->nNuclide.end())
			hatmic += mole.species[mol].den*mole_global.list[mol]->nNuclide[elHydrogen];
	}
	ASSERT(hatmic > 0.);
	hatmic = (dense.xIonDense[ipHYDROGEN][0] + dense.xIonDense[ipHYDROGEN][1])/hatmic;

	fprintf( ioQQQ, " Hydrogen     ");
	fprintf(ioQQQ,PrintEfmt("%9.2e",dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN]));
	fprintf(ioQQQ,PrintEfmt("%9.2e",dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN]));
	fprintf( ioQQQ, " H+o/Hden");
	fprintf(ioQQQ,PrintEfmt("%9.2e",hatmic ));
	fprintf(ioQQQ,PrintEfmt("%9.2e",findspecieslocal("H-")->den/dense.gas_phase[ipHYDROGEN] ));
	fprintf( ioQQQ, " H-    H2");
	/* this is total H2, the sum of "ground" and excited */
	fprintf(ioQQQ,PrintEfmt("%9.2e",hmi.H2_total/dense.gas_phase[ipHYDROGEN]));
	fprintf(ioQQQ,PrintEfmt("%9.2e",findspecieslocal("H2+")->den/dense.gas_phase[ipHYDROGEN]));
	fprintf( ioQQQ, " H2+ HeH+");
	fprintf(ioQQQ,PrintEfmt("%9.2e",findspecieslocal("HeH+")->den/dense.gas_phase[ipHYDROGEN]));
	fprintf( ioQQQ, " Ho+ ColD");
	fprintf(ioQQQ,PrintEfmt("%9.2e",findspecieslocal("H")->column));
	fprintf(ioQQQ,PrintEfmt("%9.2e",findspecieslocal("H+")->column));
	fprintf( ioQQQ, "\n");

	/* print departure coef if desired */
	if( iso_sp[ipH_LIKE][ipHYDROGEN].lgPrtDepartCoef )
	{
		fprintf( ioQQQ, " Hydrogen     " );
		fprintf(ioQQQ,PrintEfmt("%9.2e",  iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].DepartCoef()));
		fprintf(ioQQQ,PrintEfmt("%9.2e",  1.));
		fprintf( ioQQQ, " H+o/Hden");
		fprintf(ioQQQ,PrintEfmt("%9.2e", (dense.xIonDense[ipHYDROGEN][0] + dense.xIonDense[ipHYDROGEN][1])/dense.gas_phase[ipHYDROGEN]));
		fprintf(ioQQQ,PrintEfmt("%9.2e", hmi.hmidep));
		fprintf( ioQQQ, " H-    H2");
		fprintf(ioQQQ,PrintEfmt("%9.2e", hmi.h2dep));
		fprintf( ioQQQ, "      H2+");
		fprintf(ioQQQ,PrintEfmt("%9.2e", hmi.h2pdep));
		fprintf( ioQQQ, "      H3+");
		fprintf(ioQQQ,PrintEfmt("%9.2e",hmi.h3pdep));
		fprintf( ioQQQ, "\n" );
	}

	if( deut.lgElmtOn )
	{
		fprintf( ioQQQ, " Deuterium    ");
		fprintf(ioQQQ,PrintEfmt("%9.2e",deut.xIonDense[0]/deut.gas_phase));
		fprintf(ioQQQ,PrintEfmt("%9.2e",deut.xIonDense[1]/deut.gas_phase));
		fprintf( ioQQQ, " D+o/Dden");
		fprintf(ioQQQ,PrintEfmt("%9.2e",(deut.xIonDense[0] + deut.xIonDense[1])/deut.gas_phase));
		fprintf(ioQQQ,PrintEfmt("%9.2e",findspecieslocal("D-")->den/deut.gas_phase));
		fprintf( ioQQQ, " D-    HD");
		fprintf(ioQQQ,PrintEfmt("%9.2e",hmi.HD_total/deut.gas_phase));
		fprintf(ioQQQ,PrintEfmt("%9.2e",findspecieslocal("HD+")->den/deut.gas_phase));
		fprintf( ioQQQ, " HD+ HeD+");
		fprintf(ioQQQ,PrintEfmt("%9.2e",findspecieslocal("HeD+")->den/deut.gas_phase));
		fprintf( ioQQQ, " Do+ ColD");
		fprintf(ioQQQ,PrintEfmt("%9.2e",findspecieslocal("D")->column));
		fprintf(ioQQQ,PrintEfmt("%9.2e",findspecieslocal("D+")->column));
		fprintf( ioQQQ, "\n");
	}

	fixit("add a command line option to activate this");
	// print departure coefficients for diatoms species
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
	{
		if( 0 )
			(*diatom)->H2_PrtDepartCoef();
	}

	if( prt.lgPrintHeating )
	{
		fprintf( ioQQQ, "              ");
		fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(0,0)/thermal.htot));
		fprintf( ioQQQ,"         ");
		fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(0,15)/thermal.htot));
		fprintf( ioQQQ,"         ");
		fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(0,16)/thermal.htot));
		fprintf( ioQQQ,"\n");
	}

	/* convert energy flux [erg cm-2 s-1] in attenuated incident continuum, 
	 * diffuse emission, into equivalent temperature */
	double TconTot = powpq((rfield.EnergyIncidCont+rfield.EnergyDiffCont)/(4.*STEFAN_BOLTZ),1,4);
	double Tincid = powpq(rfield.EnergyIncidCont/(4.*STEFAN_BOLTZ),1,4);
	double Tdiff = powpq(rfield.EnergyDiffCont/(4.*STEFAN_BOLTZ),1,4);

	/* convert sum of flux into energy density, then equivalent pressure */
	double Pcon_nT = (Tincid + Tdiff) / SPEEDLIGHT;
	Pcon_nT /= BOLTZMANN;

	if( prt.lgPrintHeating )
	{
		fprintf( ioQQQ, "              ");
		fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(0,1)/thermal.htot ));
		fprintf( ioQQQ, "         ");
		fprintf(ioQQQ,PrintEfmt("%9.2e", 0. ));
		fprintf( ioQQQ, " BoundCom");
		fprintf(ioQQQ,PrintEfmt("%9.2e", ionbal.CompRecoilHeatLocal/ thermal.htot));
		fprintf( ioQQQ, "   Extra:");
		fprintf(ioQQQ,PrintEfmt("%9.2e",thermal.heating(0,20)/thermal.htot));
		fprintf( ioQQQ, "   Pairs:");
		fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(0,21)/ thermal.htot ));
		fprintf( ioQQQ,"  H-lines\n");
	}

	/* Helium */
	if( dense.lgElmtOn[ipHELIUM] )
	{
		fprintf( ioQQQ, " Helium       " );
		for( i=0; i < 3; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", dense.xIonDense[ipHELIUM][i]/dense.gas_phase[ipHELIUM]) );
		}

		fprintf( ioQQQ, " HeI 2s3S");
		fprintf(ioQQQ,PrintEfmt("%9.2e", 
			iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe2s3S].Pop()/dense.gas_phase[ipHELIUM] ));
		fprintf( ioQQQ, " Comp H,C");
		fprintf(ioQQQ,PrintEfmt("%9.2e", rfield.cmheat ));
		fprintf(ioQQQ,PrintEfmt("%9.2e",  rfield.cmcool*phycon.te));
		fprintf( ioQQQ , " Fill Fac");
		fprintf(ioQQQ,PrintEfmt("%9.2e", geometry.FillFac));
		fprintf( ioQQQ , " Gam1/tot");
		fprintf(ioQQQ,PrintEfmt("%9.2e", hydro.H_ion_frac_photo));
		fprintf( ioQQQ, "\n");

		/* option to print departure coef */
		if( iso_sp[ipH_LIKE][ipHELIUM].lgPrtDepartCoef )
		{
			fprintf( ioQQQ, " Helium       " );
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[0].DepartCoef()));
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipH_LIKE][ipHELIUM].st[ipH1s].DepartCoef()));
			fprintf(ioQQQ,PrintEfmt("%9.2e", 1.));

			fprintf( ioQQQ, " Comp H,C");
			fprintf(ioQQQ,PrintEfmt("%9.2e", rfield.cmheat ));
			fprintf(ioQQQ,PrintEfmt("%9.2e", rfield.cmcool*phycon.te ));
			fprintf( ioQQQ , " Fill Fac");
			fprintf(ioQQQ,PrintEfmt("%9.2e", geometry.FillFac ));
			fprintf( ioQQQ , " Gam1/tot");
			fprintf(ioQQQ,PrintEfmt("%9.2e", hydro.H_ion_frac_photo));
			fprintf( ioQQQ, "\n");
		}

		/* print heating from He (and others) if desired
		 * entry "lines" is induced line heating
		 * 1,12 ffheat:  2,3 he triplets, 1,20 compton */
		if( prt.lgPrintHeating )
		{
			/*fprintf( ioQQQ, "            %10.3e%10.3e    Lines:%10.2e%10.2e  Compton:%10.3e FF Heatig%10.3e\n", 
			  thermal.heating(1,0)/thermal.htot, thermal.heating(1,1)/
			  thermal.htot, thermal.heating(0,22)/thermal.htot, thermal.heating(1,2)/
			  thermal.htot, thermal.heating(0,19)/thermal.htot, thermal.heating(0,11)/
			  thermal.htot );*/
			fprintf( ioQQQ, "              ");
			fprintf(ioQQQ,PrintEfmt("%9.2e",thermal.heating(1,0)/thermal.htot));
			fprintf(ioQQQ,PrintEfmt("%9.2e",thermal.heating(1,1)/thermal.htot));
			fprintf( ioQQQ, "   Lines:");
			fprintf(ioQQQ,PrintEfmt("%9.2e",thermal.heating(0,22)/thermal.htot));
			fprintf(ioQQQ,PrintEfmt("%9.2e",thermal.heating(1,2)/thermal.htot));
			fprintf( ioQQQ, " Compton:");
			fprintf(ioQQQ,PrintEfmt("%9.2e",thermal.heating(0,19)/thermal.htot));
			fprintf( ioQQQ, " FFHeatig");
			fprintf(ioQQQ,PrintEfmt("%9.2e",thermal.heating(0,11)/thermal.htot));
			fprintf( ioQQQ, "\n");
		}

		if( dense.lgElmtOn[ipHELIUM] )
		{
			/* helium singlets and triplets relative to total helium gas phase density */
			double fac = 1./dense.gas_phase[ipHELIUM];
			fprintf( ioQQQ, " He singlet n " );
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe1s1S].Pop()*fac ));
			/* singlet n=2 complex */
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe2s1S].Pop()*fac ));
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe2p1P].Pop()*fac ));
			/* singlet n=3 complex */
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe3s1S].Pop()*fac ));
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe3p1P].Pop()*fac ));
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe3d1D].Pop()*fac ));

			fprintf( ioQQQ, " He tripl" );
			/* triplet n=2 complex */
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe2s3S].Pop()*fac ));
			fprintf(ioQQQ,PrintEfmt("%9.2e", 
				iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe2p3P0].Pop()*fac+
				iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe2p3P1].Pop()*fac+
				iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe2p3P2].Pop()*fac ));
			/* triplet n=3 complex */
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe3s3S].Pop()*fac ));
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe3p3P].Pop()*fac ));
			fprintf(ioQQQ,PrintEfmt("%9.2e", iso_sp[ipHE_LIKE][ipHELIUM].st[ipHe3d3D].Pop()*fac ));
			fprintf( ioQQQ, "\n" );
		}
	}

	/* loop over iso sequences to see if any populations
	 * and/or departure coefficients need to be printed */
	for( long ipISO = ipH_LIKE; ipISO < NISO; ipISO++ )
	{
		for( nelem=ipISO; nelem<LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				if( iso_sp[ipISO][nelem].lgPrtLevelPops )
				{
					iso_prt_pops(ipISO, nelem, false);
				}
				if( iso_sp[ipISO][nelem].lgPrtDepartCoef )
				{
					/* true says print departure coefficients
					 * instead of populations. */
					iso_prt_pops(ipISO, nelem, true);
				}
				/* print Critical density if desired */
				if( iso_sp[ipISO][nelem].lgPrtNCrit)
				{

					long int s =2;
					char chSpin[2][9] = {"h-like","he-like"};
					if (ipISO == ipHE_LIKE)
						s = 3;

					/* it can't be critical l-mixing density without l-changing collisions */
					if (!iso_ctrl.lgColl_l_mixing[ipISO])
					{
						fprintf (ioQQQ," It seems that l-mixing collisions are shut down for %s species\n",chSpin[ipISO]);
						fprintf (ioQQQ, "There can't be critical l-mixing density without l-changing collisions\n");
					}
					else
					{
						fprintf( ioQQQ, " %s %s Critical density \n ",elementnames.chElementSym[nelem], chSpin[ipISO]);
						fprintf(ioQQQ,"N \t crit. dens\n" ); // -log (crit dens) \n");
						for (long n=3; n<=iso_sp[ipISO][nelem].n_HighestResolved_max;n++ )
						{
							long ipHi = iso_sp[ipISO][nelem].QN2Index(n, 0, s);
							fprintf(ioQQQ,"%li \t %g\n",n,iso_sp[ipISO][nelem].st[ipHi].NCrit());
								//,-log10(iso_sp[ipISO][nelem].st[ipHi].NCrit()));
						}
					}
				}

			}
		}
	}

	/* >>chng 01 dec 08, move pressure to line before grains, after radiation properties */
	/* gas pressure, pressure due to incident radiation field, rad accel */
	fprintf( ioQQQ, " Pressure      NgasTgas");
	fprintf(ioQQQ,PrintEfmt("%9.2e", pressure.PresGasCurr/BOLTZMANN));
	fprintf( ioQQQ, " P(total)");
	fprintf(ioQQQ,PrintEfmt("%9.2e", pressure.PresTotlCurr));
	fprintf( ioQQQ, " P( gas )");
	fprintf(ioQQQ,PrintEfmt("%9.2e", pressure.PresGasCurr));
	fprintf( ioQQQ, " P(Radtn)");
	fprintf(ioQQQ,PrintEfmt("%9.2e", pressure.pres_radiation_lines_curr));
	fprintf( ioQQQ, " Rad accl");
	fprintf(ioQQQ,PrintEfmt("%9.2e", wind.AccelTotalOutward));
	fprintf( ioQQQ, " ForceMul");
	fprintf(ioQQQ,PrintEfmt("%9.2e", wind.fmul));
	fprintf( ioQQQ, "\n" );

	fprintf( ioQQQ , " Texc(La)     ");
	fprintf(ioQQQ,PrintEfmt("%9.2e",  hydro.TexcLya ));
	/* attenuated incident continuum */
	fprintf( ioQQQ , " T(I con)");
	fprintf(ioQQQ,PrintEfmt("%9.2e",  Tincid ));
	fprintf( ioQQQ , " T(D con)");
	fprintf(ioQQQ,PrintEfmt("%9.2e",  Tdiff ));
	fprintf( ioQQQ , " T(U tot)");
	fprintf(ioQQQ,PrintEfmt("%9.2e",  TconTot ));
	/* print the total radiation density expressed as an equivalent gas pressure */
	fprintf( ioQQQ , " nT (c+d)");
	fprintf(ioQQQ,PrintEfmt("%9.2e", Pcon_nT ));
	/* print the radiation to gas pressure */
	fprintf( ioQQQ , " Prad/Gas");
	fprintf(ioQQQ,PrintEfmt("%9.2e", pressure.pbeta ));
	/* magnetic to gas pressure ratio */
	fprintf( ioQQQ , " Pmag/Gas");
	fprintf(ioQQQ,PrintEfmt("%9.2e",  magnetic.pressure / pressure.PresGasCurr) );
	fprintf( ioQQQ, "\n" );

	if( gv.lgGrainPhysicsOn )
	{
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			/*  Change things so the quantum heated dust species are marked with an
			*  asterisk just after the name (K Volk)
			*  added QHMARK here and in the write statement */
			chQHMark = (char)(( gv.bin[nd].lgQHeat && gv.bin[nd].lgUseQHeat ) ? '*' : ' ');
			fprintf( ioQQQ, "%-12.12s%c  DustTemp",gv.bin[nd].chDstLab, chQHMark);
			fprintf(ioQQQ,PrintEfmt("%9.2e", gv.bin[nd].tedust));
			fprintf( ioQQQ, " Pot Volt");
			fprintf(ioQQQ,PrintEfmt("%9.2e", gv.bin[nd].dstpot*EVRYD));
			fprintf( ioQQQ, " Chrg (e)");
			fprintf(ioQQQ,PrintEfmt("%9.2e", gv.bin[nd].AveDustZ));
			fprintf( ioQQQ, " drf cm/s");
			fprintf(ioQQQ,PrintEfmt("%9.2e", gv.bin[nd].DustDftVel));
			fprintf( ioQQQ, " Heating:");
			fprintf(ioQQQ,PrintEfmt("%9.2e", gv.bin[nd].GasHeatPhotoEl));
			fprintf( ioQQQ, " Frac tot");
			fprintf(ioQQQ,PrintEfmt("%9.2e", gv.bin[nd].GasHeatPhotoEl/thermal.htot));
			fprintf( ioQQQ, "\n" );
		}
	}
	/* >>chng 00 apr 20, moved save-out of quantum heating data to qheat(), by PvH */

	/* heavy element molecules */
	if( findspecieslocal("CO")->den > 0. )
	{
		fprintf( ioQQQ, " Molecules     CH/Ctot:");
		fprintf(ioQQQ,PrintEfmt("%9.2e", findspecieslocal("CH")->den/dense.gas_phase[ipCARBON]));
		fprintf( ioQQQ, " CH+/Ctot");
		fprintf(ioQQQ,PrintEfmt("%9.2e", findspecieslocal("CH+")->den/dense.gas_phase[ipCARBON]));
		fprintf( ioQQQ, " CO/Ctot:");
		fprintf(ioQQQ,PrintEfmt("%9.2e", findspecieslocal("CO")->den/dense.gas_phase[ipCARBON]));
		fprintf( ioQQQ, " CO+/Ctot");
		fprintf(ioQQQ,PrintEfmt("%9.2e", findspecieslocal("CO+")->den/dense.gas_phase[ipCARBON]));
		fprintf( ioQQQ, " H2O/Otot");
		fprintf(ioQQQ,PrintEfmt("%9.2e", findspecieslocal("H2O")->den/dense.gas_phase[ipOXYGEN]));
		fprintf( ioQQQ, " OH/Ototl");
		fprintf(ioQQQ,PrintEfmt("%9.2e", findspecieslocal("OH")->den/dense.gas_phase[ipOXYGEN]));
		fprintf( ioQQQ, "\n");
	}

	/* information about the large H2 molecule - this just returns if not turned on */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_Prt_Zone();

	/* Lithium, Beryllium */
	if( dense.lgElmtOn[ipLITHIUM] || dense.lgElmtOn[ipBERYLLIUM] || 
		(secondaries.csupra[ipHYDROGEN][0]>0.) )
	{
		fprintf( ioQQQ, " Lithium      " );
		for( i=0; i < 4; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", dense.xIonDense[ipLITHIUM][i]/MAX2(1e-35,dense.gas_phase[ipLITHIUM]) ));
		}
		fprintf( ioQQQ, " Berylliu" );
		for( i=0; i < 5; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e",  dense.xIonDense[ipBERYLLIUM][i]/MAX2(1e-35,dense.gas_phase[ipBERYLLIUM])) );
		}

		/* print secondary ionization rate for atomic hydrogen */
		fprintf( ioQQQ, " sec ion:" );
		fprintf(ioQQQ,PrintEfmt("%9.2e", secondaries.csupra[ipHYDROGEN][0]) );
		fprintf( ioQQQ, "\n" );

		/* option to print heating due to these stages*/
		if( prt.lgPrintHeating )
		{
			fprintf( ioQQQ, "              " );
			for( i=0; i < 3; i++ )
			{
				fprintf(ioQQQ,PrintEfmt("%9.2e",  thermal.heating(ipLITHIUM,i)/ thermal.htot) );
			}
			fprintf( ioQQQ, "                    " );

			for( i=0; i < 4; i++ )
			{
				fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(ipBERYLLIUM,i)/thermal.htot ));
			}
			fprintf( ioQQQ, "\n" );
		}
	}

	/* Boron */
	if( dense.lgElmtOn[ipBORON] )
	{
		fprintf( ioQQQ, " Boron        " );
		for( i=0; i < 6; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", dense.xIonDense[ipBORON][i]/MAX2(1e-35,dense.gas_phase[ipBORON]) ));
		}
		fprintf( ioQQQ, "\n" );

		/* option to print heating*/
		if( prt.lgPrintHeating )
		{ 
			fprintf( ioQQQ, "              " );
			for( i=0; i < 5; i++ )
			{
				fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(ipBORON,i)/thermal.htot ));
			}
			fprintf( ioQQQ, "\n" );
		}
	}

	/* Carbon */
	fprintf( ioQQQ, " Carbon       " );
	for( i=0; i < 7; i++ )
	{
		fprintf(ioQQQ,PrintEfmt("%9.2e", dense.xIonDense[ipCARBON][i]/SDIV(dense.gas_phase[ipCARBON])) );
	}
	/* some molecules trail the line */
	fprintf( ioQQQ, " H2O+/O  " );
	fprintf(ioQQQ,PrintEfmt("%9.2e", findspecieslocal("H2O+")->den/MAX2(1e-35,dense.gas_phase[ipOXYGEN]) ));
	fprintf( ioQQQ, " OH+/Otot" );
	fprintf(ioQQQ,PrintEfmt("%9.2e", findspecieslocal("OH+")->den/ MAX2(1e-35,dense.gas_phase[ipOXYGEN]) ));
	/* print extra heating, normally zero */
	fprintf( ioQQQ, " Hex(tot)" );
	fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(0,20) ));
	fprintf( ioQQQ, "\n" );

	/* option to print heating*/
	if( prt.lgPrintHeating )
	{
		fprintf( ioQQQ, "              " );
		for( i=0; i < ipCARBON+1; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(ipCARBON,i)/ thermal.htot) );
		}
		fprintf( ioQQQ, "\n" );
	}

	/* Nitrogen */
	fprintf( ioQQQ, " Nitrogen     " );
	for( i=1; i <= 8; i++ )
	{
		fprintf(ioQQQ,PrintEfmt("%9.2e",dense.xIonDense[ipNITROGEN][i-1]/ SDIV(dense.gas_phase[ipNITROGEN]) ));
	}
	fprintf( ioQQQ, " O2/Ototl" );
	fprintf(ioQQQ,PrintEfmt("%9.2e", findspecieslocal("O2")->den/MAX2(1e-35,dense.gas_phase[ipOXYGEN])));
	fprintf( ioQQQ, " O2+/Otot" );
	fprintf(ioQQQ,PrintEfmt("%9.2e", findspecieslocal("O2+")->den/ MAX2(1e-35,dense.gas_phase[ipOXYGEN]) ));
	fprintf( ioQQQ, "\n" );

	/* option to print heating*/
	if( prt.lgPrintHeating )
	{
		fprintf( ioQQQ, "              " );
		for( i=0; i < ipNITROGEN+1; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(ipNITROGEN,i)/ thermal.htot ));
		}
		fprintf( ioQQQ, "\n" );
	}

#	if 0
	/* Oxygen */
	fprintf( ioQQQ, " Oxygen       " );
	for( i=1; i <= 9; i++ )
	{
		fprintf(ioQQQ,PrintEfmt("%9.2e",dense.xIonDense[ipOXYGEN][i-1]/ SDIV(dense.gas_phase[ipOXYGEN]) ));
	}
	fprintf( ioQQQ, "\n" );

	/* option to print heating*/
	if( prt.lgPrintHeating )
	{
		fprintf( ioQQQ, "              " );
		for( i=0; i < ipOXYGEN+1; i++ )
		{
			fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(ipOXYGEN,i)/ thermal.htot ));
		}
		fprintf( ioQQQ, "\n" );
	}
#	endif

	/* now print rest of elements inside loops */
	/* fluorine through Magnesium */
	for( nelem=ipOXYGEN; nelem < ipALUMINIUM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* print the element name and amount of shift */
			fprintf( ioQQQ, " %10.10s   ", elementnames.chElementName[nelem]);

			for( i=0; i < nelem+2; i++ )
			{
				fprintf(ioQQQ,PrintEfmt("%9.2e", dense.xIonDense[nelem][i]/dense.gas_phase[nelem] ));
			}
			fprintf( ioQQQ, "\n" );

			/* print heating but only if needed */
			if( prt.lgPrintHeating )
			{
				fprintf( ioQQQ, "              " );
				for( i=0; i < nelem+1; i++ )
				{
					fprintf(ioQQQ,PrintEfmt("%9.2e", thermal.heating(nelem,i)/thermal.htot ));
				}
				fprintf( ioQQQ, "\n" );
			}
		}
	}

	/* Aluminium through Zinc */
	for( nelem=ipALUMINIUM; nelem < LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			/* number of ionization stages to print across the page */
			/*@-redef@*/
			enum {LINE=13};
			/*@+redef@*/
			ishift = MAX2(0,dense.IonHigh[nelem]-LINE+1);

			/* print the element name and amount of shift */
			fprintf( ioQQQ, " %10.10s%2ld ", elementnames.chElementName[nelem],ishift );

			for( i=0; i < LINE; i++ )
			{
				fprintf(ioQQQ,PrintEfmt("%9.2e", dense.xIonDense[nelem][i+ishift]/dense.gas_phase[nelem]) );
			}
			fprintf( ioQQQ, "\n" );

			/* print heating but only if needed */
			if( prt.lgPrintHeating )
			{
				fprintf( ioQQQ, "              " );
				for( i=0; i < LINE; i++ )
				{
					fprintf(ioQQQ,
						PrintEfmt("%9.2e", thermal.heating(nelem,i+ishift)/thermal.htot ));
				}
				fprintf( ioQQQ, "\n" );
			}
		}
	}

	return;
}


