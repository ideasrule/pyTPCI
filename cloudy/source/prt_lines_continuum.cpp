/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_continuum put energetics, H, and He lines into line intensity stack */
#include "cddefines.h"
#include "taulines.h"
#include "iso.h"
#include "geometry.h"
#include "heavy.h"
#include "dense.h"
#include "prt.h"
#include "opacity.h"
#include "coolheavy.h"
#include "phycon.h"
#include "rfield.h"
#include "predcont.h"
#include "radius.h"
#include "continuum.h"
#include "lines.h"
#include "freebound.h"
#include "lines_service.h"

void lines_continuum(void)
{

	double f1, 
		f2 , 
		bac , 
		flow;
	long i,nBand;

	DEBUG_ENTRY( "lines_continuum()" );

	/* code has all local emissivities zeroed out with cryptic comment about being
	 * situation dependent.  Why?  this is option to turn back on */
	const bool KILL_CONT = false;

	i = StuffComment( "continua" );
	linadd( 0., (realnum)i , "####", 'i',
		" start continua");

	/* these entries only work correctly if the APERTURE command is not in effect */
	if( geometry.iEmissPower == 2 )
	{
		/***********************************************************************
		 * stuff in Bac ratio - continuum above the Balmer Jump 
		 * this is trick, zeroing out saved continuum integrated so far,
		 * and adding the current version, so that the line array gives the 
		 * value in the final continuum 
		 *
		 * reflected continuum is different from others since relative to inner
		 * radius, others for for this radius 
		 *************************************************************************/

		/** \todo	2	this block of lines should have nInu, InwT, InwC like main vector of continuum points */
		/***************************************************************************
		 * "Bac " , 3646,  this is residual continuum at peak of Balmer Jump
		 * flux below - flux above                                   
		 ***************************************************************************/
		/* >>chng 00 dec 02, remove opac.tmn */
		/* >>chng 00 dec 19, remove / radius.GeoDil */
		/* extrapolated continuum above head */
		/* >>chng 01 jul 13, from ConInterOut to ConEmitOut */
		f1 = (rfield.ConEmitOut[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1] + 
			rfield.ConEmitReflec[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1]/radius.r1r0sq )/
			rfield.widflx(iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1);

		/* extrapolated continuum below head */
		/* >>chng 00 dec 19, remove / radius.GeoDil */
		f2 = (rfield.ConEmitOut[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-2]+ 
			rfield.ConEmitReflec[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-2]/radius.r1r0sq )/
			rfield.widflx(iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-2);

		/* convert to nuFnu units */
		f1 = f1*0.250*0.250*EN1RYD*radius.r1r0sq;
		f2 = f2*0.250*0.250*EN1RYD*radius.r1r0sq;
		bac = (f1 - f2);

		/* memory not allocated until ipass >= 0 
		 * clear summed intrinsic and emergent intensity of following
		 * entry - following call to linadd will enter the total and
		 * keep entering the total but is done for each zone hence need to
		 * keep resetting to zero*/
		if( LineSave.ipass > 0 )
		{
			LineSave.lines[LineSave.nsum].SumLineZero();
		}

		linadd(MAX2(0.,bac)/radius.dVeffAper,3646,"Bac ",'i',
		       "residual flux at head of Balmer continuum, nuFnu ");
		/* >>chng 03 feb 06, set to zero */
		/* emslin saves the per unit vol emissivity of a line, which is normally 
		 * what goes into linadd.  We zero this unit emissivity which was set
		 * FOR THE PREVIOUS LINE since it is so situation dependent */
		if( KILL_CONT && LineSave.ipass > 0 )
		{
			LineSave.lines[LineSave.nsum-1].emslinZero();
		}

		/* memory not allocated until ipass >= 0 */
		if( LineSave.ipass > 0 )
		{
			LineSave.lines[LineSave.nsum].SumLineZero();
		}

		linadd(f1/radius.dVeffAper,3645,"nFnu",'i',
		       "total flux above head of Balmer continuum, nuFnu ");
		/* >>chng 03 feb 06, set to zero */
		/* emslin saves the per unit vol emissivity of a line, which is normally 
		 * what goes into linadd.  We zero this unit emissivity which was set
		 * FOR THE PREVIOUS LINE since it is so situation dependent */
		if( KILL_CONT && LineSave.ipass > 0 )
		{
			LineSave.lines[LineSave.nsum-1].emslinZero();
		}

		/* memory not allocated until ipass >= 0 */
		if( LineSave.ipass > 0 )
		{
			LineSave.lines[LineSave.nsum].SumLineZero();
		}

		linadd(f2/radius.dVeffAper,3647,"nFnu",'i',
		       "total flux above head of Balmer continuum, nuFnu ");
		/* >>chng 03 feb 06, set to zero */
		/* emslin saves the per unit vol emissivity of a line, which is normally 
		 * what goes into linadd.  We zero this unit emissivity which was set
		 * FOR THE PREVIOUS LINE since it is so situation dependent */
		if( KILL_CONT && LineSave.ipass > 0 )
		{
			LineSave.lines[LineSave.nsum-1].emslinZero();
		}

		/******************************************************************************
		 * "cout" , 3646,  this is outward residual continuum at peak of Balmer Jump  *
		 * equal to total in spherical geometry, half in opt thin open geometry       *
		 ******************************************************************************/
		/* >>chng 00 dec 02, remove opac.tmn */
		/* >>chng 00 dec 19, remove / radius.GeoDil */
		f1 = rfield.ConEmitOut[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1]/
			rfield.widflx(iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1);

		/* >>chng 00 dec 19, remove / radius.GeoDil */
		f2 = rfield.ConEmitOut[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-2]/
			rfield.widflx(iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-2);

		/* net Balmer jump */
		bac = (f1 - f2)*0.250*0.250*EN1RYD*radius.r1r0sq;

		/* memory not allocated until ipass >= 0 */
		if( LineSave.ipass > 0 )
		{
			LineSave.lines[LineSave.nsum].SumLineZero();
		}

		linadd(MAX2(0.,bac)/radius.dVeffAper,3646,"cout",'i',
		       "residual flux in Balmer continuum, nuFnu ");
		/* >>chng 03 feb 06, set to zero */
		/* emslin saves the per unit vol emissivity of a line, which is normally 
		 * what goes into linadd.  We zero this unit emissivity which was set
		 * FOR THE PREVIOUS LINE since it is so situation dependent */
		if( KILL_CONT && LineSave.ipass > 0 )
		{
			LineSave.lines[LineSave.nsum-1].emslinZero();
		}

		/*********************************************************************
		 * "cref" , 3646,  this is reflected continuum at peak of Balmer Jump*
		 * equal to zero in spherical geometry, half of total in op thin opn *
		 *********************************************************************/
		/* >>chng 00 dec 02, remove opac.tmn */
		/* >>chng 00 dec 19, remove / radius.GeoDil */
		f1 = rfield.ConEmitReflec[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1]/
			rfield.widflx(iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1);

		f2 = rfield.ConEmitReflec[0][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-2]/
			rfield.widflx(iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-2);

		/* net Balmer jump */
		bac = (f1 - f2)*0.250*0.250*EN1RYD;

		/* memory not allocated until ipass >= 0 */
		if( LineSave.ipass > 0 )
		{
			LineSave.lines[LineSave.nsum].SumLineZero();
		}

		linadd(MAX2(0.,bac)/radius.dVeffAper,3646,"cref",'i',
		       "residual flux in Balmer continuum, nuFnu ");
		/* >>chng 03 feb 06, set to zero */
		/* emslin saves the per unit vol emissivity of a line, which is normally 
		 * what goes into linadd.  We zero this unit emissivity which was set
		 * FOR THE PREVIOUS LINE since it is so situation dependent */
		if( KILL_CONT && LineSave.ipass > 0 )
		{
			LineSave.lines[LineSave.nsum-1].emslinZero();
		}

		/*********************************************************************
		 * "thin" , 3646, tot optically thin continuum at peak of Balmer Jump*/
		if( nzone > 0 )
		{
			/* rfield.ConEmitLocal is not defined initially, only evaluate when into model */
			f1 = rfield.ConEmitLocal[nzone][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1]/
				rfield.widflx(iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-1);
			f2 = rfield.ConEmitLocal[nzone][iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-2]/
				rfield.widflx(iso_sp[ipH_LIKE][ipHYDROGEN].fb[2].ipIsoLevNIonCon-2);
		}
		else
		{
			f1 = 0.;
			f2 = 0.;
		}

		bac = (f1 - f2)*0.250*0.250*EN1RYD;

		linadd(MAX2(0.,bac),3646,"thin",'i',
		       "residual flux in Balmer continuum, nuFnu ");

		linadd(continuum.cn4861/radius.dVeffAper,4860,"Inci",'i',
		       "incident continuum nu*f_nu at H-beta, at illuminated face of cloud ");

		linadd(continuum.cn1367/radius.dVeffAper,1367,"Inci",'i',
		       "incident continuum nu*f_nu in FUV 1367A but out of Lya damping wings, at illuminated face of cloud");

		linadd(continuum.cn2066/radius.dVeffAper,2066,"Inci",'i',
		       "incident continuum nu*f_nu in FUV 2066A at illuminated face of cloud");

		linadd(continuum.cn1216/radius.dVeffAper,1215,"Inci",'i',
		       "incident continuum nu*f_nu near Ly-alpha, at illuminated face of cloud");

		if( LineSave.ipass > 0 )
		{
			continuum.cn4861 = 0.;
			continuum.cn1216 = 0.;
			continuum.cn1367 = 0.;
			continuum.cn2066 = 0.;
		}
	}

	flow = (iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH2p].RadRecomb[ipRecRad] + 
		iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH2s].RadRecomb[ipRecRad])*
	  iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH2p].RadRecomb[ipRecEsc]*
	  dense.eden*dense.xIonDense[ipHYDROGEN][1]* 5.45e-12;
	linadd(flow,0,"Ba C",'i',
		"integrated Balmer continuum emission");

	if( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max >= 3 )
	{
		flow = ( iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH3s].RadRecomb[ipRecRad]*
			iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH3s].RadRecomb[ipRecEsc] +
			iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH3p].RadRecomb[ipRecRad]*
			iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH3p].RadRecomb[ipRecEsc] +
			iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH3d].RadRecomb[ipRecRad]*
			iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH3d].RadRecomb[ipRecEsc] ) * 
			dense.eden*dense.xIonDense[ipHYDROGEN][1]*3.53e-12;
	}
	else
	{
		flow = iso_sp[ipH_LIKE][ipHYDROGEN].fb[3].RadRecomb[ipRecRad]*
			iso_sp[ipH_LIKE][ipHYDROGEN].fb[3].RadRecomb[ipRecEsc]*
		  dense.eden*dense.xIonDense[ipHYDROGEN][1]*3.53e-12;
	}
	linadd(flow,0,"PA C",'i',
		"Paschen continuum emission ");

	/* these are a series of continuum bands defined in the file 
	 * continuum_bands.ini - this makes it possible to enter any 
	 * integrated total emission into the emission-line stack */
	/* these entries only work correctly if the APERTURE command is not in effect */
	if( geometry.iEmissPower == 2 )
	{
		for( nBand=0; nBand < continuum.nContBand; ++nBand )
		{
			double EmergentContinuum = 0.;
			double DiffuseEmission = 0.;
			if( LineSave.ipass > 0 )
			{
				/* find total emission over band - units will be erg cm-2 s-1 */
				for( i=continuum.ipContBandLow[nBand]; i<=continuum.ipContBandHi[nBand]; ++i )
				{
					// correction for fraction of low or hi cell
					// that lies within the band
					double EdgeCorrection = 1.;
					if( i==continuum.ipContBandLow[nBand] )
						EdgeCorrection = continuum.BandEdgeCorrLow[nBand];
					else if( i==continuum.ipContBandHi[nBand])
						EdgeCorrection = continuum.BandEdgeCorrHi[nBand];

					double xIntenOut = 
						/* the attenuated incident continuum */
						flux_correct_isotropic( 0, i-1 ) +

						// the outward emitted continuous radiation field 
						(rfield.ConEmitOut[0][i-1] + 

						 /* outward emitted lines */
						 rfield.outlin[0][i-1])*geometry.covgeo;
					xIntenOut *= EdgeCorrection;

					/* div by opac.E2TauAbsFace[i] because ConEmitReflec has already
					 * been corrected for this - emergent_line will introduce a 
					 * further correction so this will cancel the extra factor */
					/* NB: Comparison to 1e-37 suppresses overflows when realnum
					 *     is double (FLT_IS_DBL). */
					double xIntenIn = 0.;
					if( opac.E2TauAbsFace[i-1] > 1e-37 )
						xIntenIn = (double)rfield.ConEmitReflec[0][i-1]/
							(double)opac.E2TauAbsFace[i-1]*geometry.covgeo;
					/* outward emitted lines */
					xIntenIn += rfield.reflin[0][i-1]*geometry.covgeo;
					xIntenIn *= EdgeCorrection;

					/* the fraction of this that gets out */
					EmergentContinuum += rfield.anu(i-1) *
						emergent_line( xIntenIn , xIntenOut , i )
						/ SDIV(opac.tmn[i-1]);

					// diffuse local emission
					DiffuseEmission += (rfield.ConEmitLocal[nzone][i-1] +
							    rfield.DiffuseLineEmission[i-1])*rfield.anu(i-1)*
						EdgeCorrection;
				}
			}
			/* we will call lindst with an energy half way between the two
			 * ends of the band.  This will make an extinction correction that
			 * has already been applied above by multiplying by emergent_line.
			 * Find this factor and remove it before the call */
			double corr = emergent_line( 0.5 , 0.5 , 
						     (continuum.ipContBandLow[nBand]+continuum.ipContBandHi[nBand])/2 );
			/* NB: Comparison to 1e-37 suppresses overflows when realnum
			 *     is double (FLT_IS_DBL). */
			if( corr < 1e-37 )
				EmergentContinuum = 0.;
			else
				EmergentContinuum /= corr;

			/* convert to physical units */
			EmergentContinuum *= EN1RYD*radius.r1r0sq/radius.dVeffAper;
			DiffuseEmission *= EN1RYD;
			/* memory not allocated until ipass >= 0 */
			if( LineSave.ipass > 0 )
			{
				LineSave.lines[LineSave.nsum].SumLineZero();
			}

			lindst( EmergentContinuum,
					// the negative wavelength is a sentinel that the wavelength
					// may be bogus - very often the "center" of the band, as defined
					// by observers, is far to the blue end of the range.  use the
					// observed definition of the wavelength.  This introduces an error
					// since the wavelength is used to determine the transfer of the
					// band against background opacities
					-continuum.ContBandWavelength[nBand],
					continuum.chContBandLabels[nBand].c_str(),
					(continuum.ipContBandLow[nBand]+continuum.ipContBandHi[nBand])/2,
					't', false,
					"continuum bands defined in continuum_bands.ini");

			// emissivity has no meaning for these bands - quantity is net
			// transmitted radiation field
			if( LineSave.ipass > 0 )
			{
				LineSave.lines[LineSave.nsum-1].emslinSet(0,DiffuseEmission);
				LineSave.lines[LineSave.nsum-1].emslinThin();
			}
		}
	}

	linadd(MAX2(0.,CoolHeavy.brems_cool_net),0,"HFFc",'c',
		"net free-free cooling, ALL species, free-free heating subtracted, so nearly cancels with cooling in LTE ");

	linadd(MAX2(0.,-CoolHeavy.brems_cool_net),0,"HFFh",'h',
		"net free-free heating, nearly cancels with cooling in LTE ");

	linadd(CoolHeavy.brems_cool_h,0,"H FF",'i',
		" H brems (free-free) cooling ");

	linadd(CoolHeavy.brems_heat_total,0,"FF H",'i',
		"total free-free heating ");

	linadd(CoolHeavy.brems_cool_he,0,"HeFF",'i',
		"He brems emission ");

	linadd(CoolHeavy.heavfb,0,"MeFB",'c',
		"heavy element recombination cooling ");

	linadd(CoolHeavy.brems_cool_metals,0,"MeFF",'i',
		"heavy elements (metals) brems cooling, heat not subtracted ");

	linadd(CoolHeavy.brems_cool_h+CoolHeavy.brems_cool_he+CoolHeavy.brems_cool_metals,0,"ToFF",'i',
		"total brems emission - total cooling but not minus heating ");

	linadd((CoolHeavy.brems_cool_h+CoolHeavy.brems_cool_he)*sexp(5.8e6/phycon.te),0,"FF X",'i',
		"part of H brems, in x-ray beyond 0.5KeV ");

	linadd(CoolHeavy.eebrm,0,"eeff",'c',
		"electron - electron brems ");

	linadd(CoolHeavy.colmet,0,"Mion",'c',
		" cooling due to collisional ionization of heavy elements" );

	/* predict emitted continuum at series of continuum points */
	/* class is located in predcont.h, 
	 * PredCont - contains vector of pair of Energy and ip where
	 *            we want to predict the continuum,
	 *
	 * the entry nFnu will only be printed if the command
	 * print diffuse continuum
	 * is entered - 
	 *
	 * this code should be kept parallel with that in dopunch, where
	 * save continuum is produced, since two must agree */

	t_PredCont& PredCont = t_PredCont::Inst();
	if( LineSave.ipass == 0 )
		PredCont.set_offset(LineSave.nsum);

	/* these entries only work correctly if the APERTURE command is not in effect */
	if( geometry.iEmissPower == 2 )
	{
		for( i=0; i < long(PredCont.size()); i++ )
		{
			double SourceTransmitted , Cont_nInu;
			double SourceReflected, DiffuseOutward, DiffuseInward;
			double renorm;

			/* put wavelength in Angstroms into dummy structure, so that we can use iWavLen
			 * to get a proper wavelength with units, continuum energies are stored in PredCont */
			(*TauDummy).WLAng() = (realnum)PredCont[i].Angstrom();
			/*lambda = iWavLen(TauDummy , &chUnits , &chShift );*/

			/* >>chng 00 dec 02, there were three occurrences of /opac.tmn which had the
			 * effect of raising the summed continuum by the local opacity correction factor.
			 * in the case of the Lyman continuum this raised the reported value by orders
			 * of magnitude.  There have been commented out in the following for now. */
			/* reflected total continuum (diff+incident emitted inward direction) */

			/* >>chng 00 dec 08, implement the "set nFnu [SOURCE_REFLECTED] ... command, PvH */
			/* >>chng 00 dec 19, remove / radius.GeoDil */
			renorm = rfield.anu2(PredCont[i].ip_C())*EN1RYD/rfield.widflx(PredCont[i].ip_C());

			/* this is the reflected diffuse continuum */
			if( prt.lgDiffuseInward )
			{
				DiffuseInward = rfield.ConEmitReflec[0][PredCont[i].ip_C()]*renorm;
			}
			else
			{
				DiffuseInward = 0.;
			}

			/* the outward diffuse continuum */
			if( prt.lgDiffuseOutward )
			{
				DiffuseOutward = rfield.ConEmitOut[0][PredCont[i].ip_C()]*renorm*radius.r1r0sq;
			}
			else
			{
				DiffuseOutward = 0.;
			}

			/* reflected part of INCIDENT continuum (only incident, not diffuse, which was above) */
			if( prt.lgSourceReflected )
			{
				SourceReflected =  rfield.ConRefIncid[0][PredCont[i].ip_C()]*renorm;
			}
			else
			{
				SourceReflected =  0.;
			}

			/* the attenuated incident continuum */
			if( prt.lgSourceTransmitted )
			{
				SourceTransmitted = rfield.flux[0][PredCont[i].ip_C()]*renorm*radius.r1r0sq;
			}
			else
			{
				SourceTransmitted = 0.;
			}

			/* memory has not been allocated until ipass >= 0, so must not access this element,
			 * this element will be used to save the following quantity */
			if( LineSave.ipass > 0 )
			{
				LineSave.lines[LineSave.nsum].SumLineZero();
			}

			linadd((DiffuseInward+SourceReflected+DiffuseOutward+SourceTransmitted)/radius.dVeffAper,
				(*TauDummy).WLAng(),"nFnu",'i',
				"total continuum at selected energy points " );

			/* emslin saves the per unit vol emissivity of a line, which is normally 
			 * what goes into linadd.  We zero this unit emissivity which was set
			 * FOR THE PREVIOUS LINE since it is so situation dependent */
			if( KILL_CONT && LineSave.ipass > 0 )
			{
				LineSave.lines[LineSave.nsum-1].emslinZero();
			}

			/* this is the normal set to zero to trick the NEXT line into going in properly */
			if( LineSave.ipass > 0 )
			{
				LineSave.lines[LineSave.nsum].SumLineZero();
			}

			/* the nsum-1 -- emslin and nsum -- SumLine is not a bug, look above - they do
			 * different things to different saves */
			Cont_nInu = rfield.flux[0][PredCont[i].ip_C()]*renorm*radius.r1r0sq +
				rfield.ConRefIncid[0][PredCont[i].ip_C()]*renorm;

#			if 0
			/* this code can be used to create assert statements for the continuum shape */
			if( !i )
				fprintf(ioQQQ,"\n");
			string chWL;
			sprt_wl( chWL , (*TauDummy).WLAng() );
			fprintf( ioQQQ,"assert line luminosity \"nInu\" %s  %.3f\n",
				 chWL.c_str(), 
				log10(SDIV(Cont_nInu/radius.dVeffAper) * radius.Conv2PrtInten)  );
#			endif

			linadd( Cont_nInu/radius.dVeffAper,(*TauDummy).WLAng(),"nInu",'i',
				"transmitted and reflected incident continuum at selected energy points " );

			/* emslin saves the per unit volume emissivity of a line, which is normally 
			 * what goes into linadd.  We zero this unit emissivity since it is so situation dependent */
			if( KILL_CONT && LineSave.ipass > 0 )
			{
				LineSave.lines[LineSave.nsum-1].emslinZero();
			}

			/* memory has not been allocated until ipass >= 0 */
			if( LineSave.ipass > 0 )
			{
				LineSave.lines[LineSave.nsum].SumLineZero();
			}

			linadd( (DiffuseInward+SourceReflected)/radius.dVeffAper,(*TauDummy).WLAng(),"InwT",'i',
				"total reflected continuum, total inward emission plus reflected (diffuse) total continuum ");

			if( KILL_CONT && LineSave.ipass > 0 )
			{
				LineSave.lines[LineSave.nsum-1].emslinZero();
			}

			/* memory has not been allocated until ipass >= 0 */
			if( LineSave.ipass > 0 )
			{
				LineSave.lines[LineSave.nsum].SumLineZero();
			}

			linadd(SourceReflected/radius.dVeffAper,(*TauDummy).WLAng(),"InwC",'i',
				"reflected incident continuum (only incident) ");

			if( KILL_CONT && LineSave.ipass > 0 )
			{
				LineSave.lines[LineSave.nsum-1].emslinZero();
			}
		}
	}
	
	i = StuffComment( "RRC" );
	linadd( 0., (realnum)i , "####", 'i',"radiative recombination continua");
	
	// radiative recombination continua, RRC, for iso sequences
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( long n=0; n < iso_sp[ipISO][nelem].numLevels_max; n++ )
				{
					if( LineSave.ipass < 0 )
						// this pass only counting lines
						linadd(0.,0.,"dumy",'i',"radiative recombination continuum");
					else if( LineSave.ipass == 0 )
					{
						// save wavelength and label
						/* chIonLbl generates a null terminated 4 char string, of form "C  2"
						 * the result, chLable, is only used when ipass == 0, can be undefined otherwise */
						realnum wl = (realnum)(RYDLAM / iso_sp[ipISO][nelem].fb[n].xIsoLevNIonRyd);
						wl /= (realnum)RefIndex( 1e8/wl );
						linadd( 0. , wl ,chIonLbl(iso_sp[ipISO][nelem].trans(1,0)).c_str(),'i',
							"radiative recombination continuum");
					}
					else
					{
						// save intensity
						linadd(iso_sp[ipISO][nelem].fb[n].RadRecCon,0,"dumy",'i',
							"radiative recombination continuum");
					}
				}
			}
		}
	}

	// RRC for non iso sequence ions

	/* add recombination continua for elements heavier than those done with iso seq */
	for( long nelem=NISO; nelem < LIMELM; nelem++ )
	{
		/* do not include species with iso-sequence in following */
		/* >>chng 03 sep 09, upper bound was wrong, did not include NISO */
		for( long ion=0; ion < nelem-NISO+1; ion++ )
		{
			if(  dense.lgElmtOn[nelem] )
			{
				if( LineSave.ipass < 0 )
					// this pass only counting lines
					linadd(0.,0.,"dumy",'i',"radiative recombination continuum");
				else if( LineSave.ipass == 0 )
				{
					string chLabel = chIonLbl( nelem+1, ion+1 );
					realnum wl = (realnum)(RYDLAM / Heavy.Valence_IP_Ryd[nelem][ion]);
					wl /= (realnum)RefIndex( 1e8/wl );
					linadd( 0. , wl ,chLabel.c_str(),'i',
						"radiative recombination continuum");
				}
				else
				{
					// save intensity
					linadd(Heavy.RadRecCon[nelem][ion],0,"dumy",'i',
						"radiative recombination continuum");
				}
			}
		}
	}

	return;
}
