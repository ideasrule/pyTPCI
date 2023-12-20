/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ion_photo fill array PhotoRate with photoionization rates for heavy elements */
#include "cddefines.h"
#include "yield.h"
#include "heavy.h"
#include "opacity.h"
#include "thermal.h"
#include "conv.h"
#include "grainvar.h"
#include "elementnames.h"
#include "ionbal.h"
#include "phycon.h"
#include "hmi.h"
#include "rfield.h"
#include "atoms.h"
#include "iso.h"
#include "oxy.h"
#include "atmdat.h"
#include "fe.h"
#include "gammas.h"
#include "mole.h"
#include "freebound.h"
#include "dense.h"
#include "taulines.h"
#include "ipoint.h"

void ion_photo(
	/* nelem is atomic number on C scale, 0 for H */
	long int nelem , 
	/* debugging flag to turn on print */
	bool lgPrintIt )
{
	long int ion, 
	  iphi, 
	  iplow, 
	  ipop, 
	  limit_hi, 
	  limit_lo,
	  ns;

	DEBUG_ENTRY( "ion_photo()" );

	/* IonLow(nelem) and IonHigh(nelem) are bounds for evaluation*/

	/*begin sanity checks */
	ASSERT( nelem < LIMELM );
	ASSERT( dense.IonLow[nelem] >= 0 );
	ASSERT( dense.IonLow[nelem] <= dense.IonHigh[nelem] );
	ASSERT( dense.IonHigh[nelem] <= nelem + 1);
	/*end sanity checks */

	/* NB - in following, nelem is on c scale, so is 29 for Zn */

	/* min since iso-like atom rates produced in iso_photo.
	 * IonHigh and IonLow range from 0 to nelem+1, but for photo rates
	 * we want 0 to nelem since cannot destroy fully stripped ion.  
	 * since iso-seq done elsewhere, want to actually do IonHigh-xxx.*/
	/* >>chng 00 dec 07, logic on limit_hi now precisely identical to ion_solver */
	/* >>chng 02 mar 31, change limit_hi to < in loop (had been <=) and
	 * also coded for arbitrary number of iso sequences */
	/*limit_hi = MIN2(nelem,dense.IonHigh[nelem]);
	limit_hi = MIN2(nelem-2,dense.IonHigh[nelem]-1);*/
	limit_hi = MIN2( dense.IonHigh[nelem] , nelem+1-NISO );

	/* >>chng 03 sep 26, always do atom itself since may be needed for molecules */
	limit_hi = MAX2( 1 , limit_hi );

	/* when grains are present want to use atoms as lower bound to number of stages of ionization,
	 * since atomic rates needed for species within grains */
	if( !conv.nPres2Ioniz && gv.lgDustOn() )
	{
		limit_lo = 0;
	}
	else
	{
		limit_lo = dense.IonLow[nelem];
	}

	/* >>chng 01 dec 11, lower bound now limit_lo */
	/* loop over all ions for this element */
	/* >>chng 02 mar 31, now ion < limit_hi not <= */
	for( ion=limit_lo; ion < limit_hi; ion++ )
	{
		/* loop over all shells for this ion */
		for( ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
		{
			/* always reevaluate the outer shell, and all shells if lgRedoStatic is set */
			if( (ns==(Heavy.nsShells[nelem][ion]-1) || opac.lgRedoStatic) )
			{
				/* option to redo the rates only on occasion */
				iplow = opac.ipElement[nelem][ion][ns][0];
				iphi = opac.ipElement[nelem][ion][ns][1];
				ipop = opac.ipElement[nelem][ion][ns][2];

				t_phoHeat photoHeat;

				/* compute the photoionization rate, ionbal.lgPhotoIoniz_On is 1, set 0
				 * with "no photoionization" command */
				ionbal.PhotoRate_Shell[nelem][ion][ns][0] =
					GammaK(iplow,iphi,
					ipop,t_yield::Inst().elec_eject_frac(nelem,ion,ns,0),
					&photoHeat )*ionbal.lgPhotoIoniz_On;

				/* these three lines must be kept parallel with the lines
				 * in GammaK ion*/

				/* the heating rate */
				ionbal.PhotoRate_Shell[nelem][ion][ns][1] = photoHeat.HeatLowEnr*ionbal.lgPhotoIoniz_On;
				ionbal.PhotoRate_Shell[nelem][ion][ns][2] = photoHeat.HeatHiEnr*ionbal.lgPhotoIoniz_On;
			}
		}

		/* add on compton recoil ionization for atoms to outer shell */
		/* >>chng 02 mar 24, moved here from ion_solver */
		/* this is the outer shell */
		ns = (Heavy.nsShells[nelem][ion]-1);
		/* this must be moved to photoionize and have code parallel to iso_photo code */
		ionbal.PhotoRate_Shell[nelem][ion][ns][0] += ionbal.CompRecoilIonRate[nelem][ion];
		/* add the heat as secondary-ionization capable heating */
		ionbal.PhotoRate_Shell[nelem][ion][ns][2] += ionbal.CompRecoilHeatRate[nelem][ion];
	}

	/* option to print information about these rates for this element */
	if( lgPrintIt )
	{
		/* option to print rates for particular shell */
		ns = 5;
		ion = 1;
		GammaPrt(
		  opac.ipElement[nelem][ion][ns][0], 
		  opac.ipElement[nelem][ion][ns][1], 
		  opac.ipElement[nelem][ion][ns][2], 
		  ioQQQ, /* io unit we will write to */
		  ionbal.PhotoRate_Shell[nelem][ion][ns][0], 
		  0.05);

		/* outer loop is from K to most number of shells present in atom */
		for( ns=0; ns < Heavy.nsShells[nelem][0]; ns++ )
		{
			fprintf( ioQQQ, "\n %s", elementnames.chElementNameShort[nelem] );
			fprintf( ioQQQ, " %s" , Heavy.chShell[ns]);
			/* MB hydrogenic photo rate may not be included in beow */
			for( ion=0; ion < dense.IonHigh[nelem]; ion++ )
			{
				if( Heavy.nsShells[nelem][ion] > ns )
				{
					fprintf( ioQQQ, " %8.1e", ionbal.PhotoRate_Shell[nelem][ion][ns][0] );
				}
				else
				{
					break;
				}
			}
		}
		fprintf(ioQQQ,"\n");
	}

	/* Ly-alpha photoionization of Ca+
	 * valence shell is reevaluated by ion_photo on every call, so this does not double count */
	if( nelem == ipCALCIUM )
	{

		long id_Ca2;

		if( (id_Ca2 = atmdat.ipSpecIon[nelem][1]) < 0 )
		{
			return;
		}

		// H I Lya OTS rates
		double hlgam = rfield.otslin[ iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).ipCont() -1];
		double PhotoRate5 = 1.7e-18*hlgam,
		PhotoRate4 = 8.4e-19*hlgam,
		PhotoRate3 = 7.0e-18*hlgam,
		PhotoRate2 = 4.8e-18*hlgam;

		realnum dstCala = (realnum)(dBaseStates[id_Ca2][4].Pop()*PhotoRate5 +
				dBaseStates[id_Ca2][3].Pop()*PhotoRate4);
		realnum dCaf12 = (realnum)((dBaseStates[id_Ca2][2].Pop()*PhotoRate3 +
				dBaseStates[id_Ca2][1].Pop()*PhotoRate2)*2.31e-12);

		double denseCa2 = 0.;
		for( long i = 0; i < dBaseSpecies[id_Ca2].numLevels_local; ++i)
		{
			denseCa2 += dBaseStates[id_Ca2][i].Pop();
		}

		/* total rate Lalpha destroys CaII,
		 * this is only used in ioncali to increase ionization rate by
		 * adding it to the ct vector */
		if( denseCa2 > 0. )
		{
			dstCala = (realnum)(
				(dstCala + dCaf12/2.31e-12)/denseCa2);
			{
				/*@-redef@*/
				enum {DEBUG_LOC=false};
				/*@+redef@*/
				if( DEBUG_LOC )
				{
					fprintf(ioQQQ," hlgam is %e\n", hlgam);
				}
			}
		}
		else
		{
			dstCala = 0.;
		}

		long ns = 6, ion = 1;
		ionbal.PhotoRate_Shell[nelem][ion][ns][0] += dstCala;
	}

	if( nelem == ipCARBON )
	{
		/* >>chng 05 aug 10, add Leiden hack here to get same C0 photo rate
		 * as UMIST - negates difference in grain opacities */
		if(mole_global.lgLeidenHack)
		{
			int nelem=ipCARBON , ion=0 , ns=2;
			ionbal.PhotoRate_Shell[nelem][ion][ns][0] =
				(hmrate4((1e-10)*3.0,0,0,phycon.te)*(hmi.UV_Cont_rel2_Habing_TH85_face*
				exp(-(3.0*rfield.extin_mag_V_point))/1.66));
			/* heating rates */
			ionbal.PhotoRate_Shell[nelem][ion][ns][1] = 0.;
			ionbal.PhotoRate_Shell[nelem][ion][ns][2] = 0.;
		}
	}

	if( nelem == ipNITROGEN )
	{
		/* photoionization from 2D of NI*/
		t_phoHeat photoHeat;

		long id_n1;
		double pop1 = 0., pop2 = 0.;

		double denseN1 = 0.;

		if( (id_n1 = atmdat.ipSpecIon[nelem][0]) >= 0 )
		{

			for( long i = 0; i < dBaseSpecies[id_n1].numLevels_local; ++i)
			{
				denseN1 += dBaseStates[id_n1][i].Pop();
			}

			if( denseN1 > SMALLFLOAT )
			{
				pop1 = dBaseStates[id_n1][1].Pop();
				pop2 = dBaseStates[id_n1][2].Pop();

				atoms.p2nit = (realnum)(pop1+pop2)/denseN1;
				// photo rate, population atoms.p2nit evaluated in cooling
				atoms.d5200r = (realnum)GammaK(opac.in1[0],opac.in1[1],opac.in1[2],1.,&photoHeat);
				/* valence shell photoionization, followed by heating; [0][6] => atomic nitrogen */
				ionbal.PhotoRate_Shell[ipNITROGEN][0][2][0] = ionbal.PhotoRate_Shell[ipNITROGEN][0][2][0]*
				  (1. - atoms.p2nit) + atoms.p2nit*atoms.d5200r;
				ionbal.PhotoRate_Shell[ipNITROGEN][0][2][1] = ionbal.PhotoRate_Shell[ipNITROGEN][0][2][1]*
				  (1. - atoms.p2nit) + photoHeat.HeatNet*atoms.p2nit;
			}
			else
			{
				atoms.p2nit = 0.;
				atoms.d5200r = 0.;
			}
		}
		else
		{
			atoms.p2nit = 0.;
			atoms.d5200r = 0.;
		}
	}

	if ( nelem == ipMAGNESIUM )
	{
		if( dense.IonLow[ipMAGNESIUM] <= 1 )
		{

			long id_mg2;
			double pop1 = 0., pop2 = 0.;

			double denseMg2 = 0.;

			if( (id_mg2 = atmdat.ipSpecIon[nelem][1]) >= 0 )
			{

				for( long i = 0; i < dBaseSpecies[id_mg2].numLevels_local; ++i)
				{
					denseMg2 += dBaseStates[id_mg2][i].Pop();
				}

				if( denseMg2 > SMALLFLOAT )
				{
					pop1 = dBaseStates[id_mg2][1].Pop();
					pop2 = dBaseStates[id_mg2][2].Pop();

					atoms.popMg2 = (realnum)(pop1+pop2)/denseMg2;

					t_phoHeat dummy;
					/* photoionization from excited upper state of 2798 */
					double rmg2l = GammaK(opac.ipmgex,
					iso_sp[ipH_LIKE][ipHYDROGEN].fb[0].ipIsoLevNIonCon,opac.ipOpMgEx,1., &dummy );
					ionbal.PhotoRate_Shell[ipMAGNESIUM][1][3][0] += rmg2l*atoms.popMg2;

					if( ionbal.PhotoRate_Shell[ipMAGNESIUM][1][3][0] > 1e-30 )
					{
						/* remember max relative photoionization rate for possible comment */
						atoms.xMg2Max = (realnum)(MAX2(atoms.xMg2Max,rmg2l*atoms.popMg2/
								ionbal.PhotoRate_Shell[ipMAGNESIUM][1][3][0]));
					}
				}
				else
				{
					atoms.xMg2Max = 0.;
					atoms.popMg2 = 0.;
				}
			}
			else
			{
				atoms.xMg2Max = 0.;
				atoms.popMg2 = 0.;
			}
		}
		else
		{
			atoms.xMg2Max = 0.;
			atoms.popMg2 = 0.;
		}
	}

	if( nelem == ipOXYGEN )
	{
		t_phoHeat dummy;

		/* O++ */
		oxy.poiii2Max = 0.;
		oxy.poiii3Max = 0.;
		oxy.poiii2 = 0.;
		oxy.poiii3 = 0.;
		oxy.r4363Max = 0.;
		oxy.r5007Max = 0.;
		oxy.d5007r = 0.;

		/* neutral O */
		oxy.poiexc = 0.;
		oxy.poimax = 0.;
		oxy.d6300 = 0.;

		/* other transitions */
		oxy.s3727 = 0.;
		oxy.s7325 = 0.;
		oxy.AugerO3 = 0.;


		/* oxygen, atomic number 8 */
		if( !dense.lgElmtOn[ipOXYGEN] )
		{
			thermal.setHeating(7,9,0.);
			return;
		}

		long id_o3, id_o1;;
		if( (id_o3 = atmdat.ipSpecIon[nelem][2]) >= 0 )
		{
			double	denseO3 = 0.;
			for( long i = 0; i < dBaseSpecies[id_o3].numLevels_local; ++i)
			{
				denseO3 += dBaseStates[id_o3][i].Pop();
			}
			if( denseO3 > SMALLFLOAT )
			{
				oxy.poiii3 = dBaseStates[id_o3][4].Pop()/denseO3;
				oxy.poiii2 = dBaseStates[id_o3][3].Pop()/denseO3;

				double	denseO3_First3 = (dBaseStates[id_o3][0].Pop() +
						dBaseStates[id_o3][1].Pop() +
						dBaseStates[id_o3][2].Pop())/denseO3;

				ASSERT(oxy.poiii3 <= 1. && oxy.poiii3 >= 0.);
				ASSERT(oxy.poiii2 <= 1. && oxy.poiii2 >= 0.);

				/* photoionization from O++ 1D
				 *
				 * estimate gamma function by assuming no frequency dependence
				 * betwen 1D and O++3P edge */
				/* destroy upper level of OIII 5007*/
				oxy.d5007r = (realnum)(GammaK(opac.ipo3exc[0],opac.ipo3exc[1],
				  opac.ipo3exc[2] , 1., &dummy ));

				/* destroy upper level of OIII 4363*/
				realnum oxy_d4363 = (realnum)(GammaK(opac.ipo3exc3[0],opac.ipo3exc3[1],
				  opac.ipo3exc3[2] , 1., &dummy ));

				/* photoionization from excited states */
				/* doubly ionized oxygen destruction */
				ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0] = ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0]*
				  denseO3_First3 + oxy.d5007r*oxy.poiii2 + oxy_d4363*oxy.poiii3;

				if( ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0] > 1e-30 && dense.IonLow[ipOXYGEN] <= 2 )
				{
					if( (oxy.d5007r*oxy.poiii2 + oxy_d4363*oxy.poiii3)/
					  ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0] > (oxy.r4363Max +
					  oxy.r5007Max) )
					{
						oxy.poiii2Max = (realnum)(oxy.d5007r*oxy.poiii2/ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0]);
						oxy.poiii3Max = (realnum)(oxy_d4363*oxy.poiii3/ionbal.PhotoRate_Shell[ipOXYGEN][2][2][0]);
					}
					oxy.r4363Max = (realnum)(MAX2(oxy.r4363Max,oxy_d4363));
					oxy.r5007Max = (realnum)(MAX2(oxy.r5007Max,oxy.d5007r));
				}
			}
		}

		if( (id_o1 = atmdat.ipSpecIon[nelem][0]) >= 0 )
		{
			double denseO1 = 0.;
			for( long i = 0; i < dBaseSpecies[id_o1].numLevels_local; ++i)
			{
				denseO1 += dBaseStates[id_o1][i].Pop();
			}

			if( denseO1 > SMALLFLOAT )
			{
				double	denseO1_First3 = (dBaseStates[id_o1][0].Pop() +
						dBaseStates[id_o1][1].Pop() +
						dBaseStates[id_o1][2].Pop())/denseO1;

				oxy.poiexc = dBaseStates[id_o1][3].Pop()/denseO1;
				ASSERT(oxy.poiexc <= 1. && oxy.poiexc >= 0.);

				/* destroy upper level of OI 6300*/
				oxy.d6300 = (realnum)(GammaK(opac.ipo1exc[0],opac.ipo1exc[1],
				  opac.ipo1exc[2] , 1., &dummy ));

				/* photoionization from excited states */
				/* neutral oxygen destruction */
				ionbal.PhotoRate_Shell[ipOXYGEN][0][2][0] = ionbal.PhotoRate_Shell[ipOXYGEN][0][2][0]*
				  denseO1_First3 + oxy.d6300*oxy.poiexc;

				/* ct into excited states */
				if( dense.IonLow[ipOXYGEN] <= 0 && (ionbal.PhotoRate_Shell[ipOXYGEN][0][2][0] +
				  atmdat.CharExcIonOf[ipHYDROGEN][ipOXYGEN][0]*dense.xIonDense[ipHYDROGEN][1]) > 1e-30 )
				{
					oxy.poimax = (realnum)(MAX2(oxy.poimax,oxy.d6300*oxy.poiexc/
					  (ionbal.PhotoRate_Shell[ipOXYGEN][0][2][0]+
					  atmdat.CharExcIonOf[ipHYDROGEN][ipOXYGEN][0]* dense.xIonDense[ipHYDROGEN][1])));
				}
			}
		}

		/* save atomic oxygen photodestruction rate for 3727 creation */
		if( dense.IonLow[ipOXYGEN] == 0 && oxy.i2d < rfield.nflux )
		{
			oxy.s3727 = (realnum)(GammaK(oxy.i2d,oxy.i2p,opac.iopo2d , 1., &dummy ));

			long iup = ipoint(atmdat.EIonPot[ipHELIUM][1]);
			oxy.s7325 = (realnum)(GammaK(oxy.i2d,iup,opac.iopo2d , 1., &dummy ));

			oxy.s7325 -= oxy.s3727;
			oxy.s3727 = oxy.s3727 + oxy.s7325;

			/* ratio of cross sections */
			oxy.s7325 *= 0.66f;
		}

		oxy.AugerO3 = (realnum)ionbal.PhotoRate_Shell[ipOXYGEN][0][0][0];

		oxy.s3727 *= dense.xIonDense[ipOXYGEN][0];
		oxy.s7325 *= dense.xIonDense[ipOXYGEN][0];
		oxy.AugerO3 *= dense.xIonDense[ipOXYGEN][0];
	}

	if( nelem == ipIRON )
	{
		if( !dense.lgElmtOn[ipIRON] )
		{
			fe.fekcld = 0.;
			fe.fekhot = 0.;
			fe.fegrain = 0.;
		}
		else
		{
			const int NDIM = ipIRON+1;

			static const double fyield[NDIM+1] = {.34,.34,.35,.35,.36,.37,.37,.38,.39,.40,
			  .41,.42,.43,.44,.45,.46,.47,.47,.48,.48,.49,.49,.11,.75,0.,0.,0.};

			long int i, limit, limit2;
			/* ntotal fluorescence production of K-alpha
			 * "cold" iron has M-shell electrons, up to Fe 18 */
			fe.fekcld = 0.;
			limit = MIN2(18,dense.IonHigh[ipIRON]);

			for( i=dense.IonLow[ipIRON]; i < limit; i++ )
			{
				ASSERT( i < NDIM + 1 );
				fe.fekcld +=
					(realnum)(ionbal.PhotoRate_Shell[ipIRON][i][0][0]*dense.xIonDense[ipIRON][i]*
				  fyield[i]);
			}

			/* same sum for hot iron, do not include one and two electron Fe */
			fe.fekhot = 0.;
			limit = MAX2(18,dense.IonLow[ipIRON]);

			limit2 = MIN2(ipIRON-1,dense.IonHigh[ipIRON]);
			ASSERT( limit2 <= LIMELM + 1 );

			for( i=limit; i < limit2; i++ )
			{
				ASSERT( i < NDIM + 1 );
				fe.fekhot +=
					(realnum)(ionbal.PhotoRate_Shell[ipIRON][i][0][0]*dense.xIonDense[ipIRON][i]*
				  fyield[i]);
			}

			/* Fe Ka from grains - Fe in grains assumed to be atomic
			 * gv.elmSumAbund[ipIRON] is number density of iron added over all grain species */
			i = 0;
			/* fyield is 0.34 for atomic fe */
			fe.fegrain = ( gv.lgWD01 ) ? 0.f : (realnum)(ionbal.PhotoRate_Shell[ipIRON][i][0][0]*fyield[i]*
						 gv.elmSumAbund[ipIRON]);
		}
	}

	return;
}
