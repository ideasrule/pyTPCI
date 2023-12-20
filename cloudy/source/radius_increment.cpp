/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*radius_increment do work associated with geometry increments of this zone, called before RT_tau_inc */
#include "cddefines.h"
#include "iso.h"
#include "hydrogenic.h"
#include "colden.h"
#include "geometry.h"
#include "iterations.h"
#include "opacity.h"
#include "thermal.h"
#include "dense.h"
#include "h2.h"
#include "timesc.h"
#include "hmi.h"
#include "taulines.h"
#include "trace.h"
#include "wind.h"
#include "phycon.h"
#include "pressure.h"
#include "grainvar.h"
#include "molcol.h"
#include "conv.h"
#include "hyperfine.h"
#include "mean.h"
#include "struc.h"
#include "radius.h"
#include "gravity.h"
#include "mole.h"
#include "rfield.h"
#include "doppvel.h"
#include "freebound.h"
#include "dynamics.h"

void radius_increment(void)
{

	long int nelem,
	  ion,
	  nzone_minus_1 , 
	  mol;
	double 
	  avWeight,
	  t;

	double ajmass, 
		Error;

	DEBUG_ENTRY( "radius_increment()" );

	/* when this sub is called radius is the outer edge of zone */

	/* save information about structure of model
	 * first zone is 1 but array starts at 0 - nzone_minus_1 is current zone
	 * max nzone because abort during search phase gets here with nzone = -1 */
	nzone_minus_1 = MAX2( 0, nzone-1 );
	ASSERT(nzone_minus_1>=0 && nzone_minus_1 < struc.nzlim );

	struc.heatstr[nzone_minus_1] = thermal.htot - dynamics.Heat();
	struc.coolstr[nzone_minus_1] = thermal.ctot - dynamics.Cool();
	struc.testr[nzone_minus_1] = (realnum)phycon.te;

	/* number of particles per unit vol */
	struc.DenParticles[nzone_minus_1] = dense.pden;
	/* rjrw: add hden for dilution */
	struc.hden[nzone_minus_1] = (realnum)dense.gas_phase[ipHYDROGEN];
	/* total grams per unit vol */
	struc.DenMass[nzone_minus_1] = dense.xMassDensity;
	struc.volstr[nzone_minus_1] = (realnum)radius.dVeffAper;
	struc.drad[nzone_minus_1] = (realnum)radius.drad;
	struc.drad_x_fillfac[nzone_minus_1] = (realnum)radius.drad_x_fillfac;
	struc.histr[nzone_minus_1] = dense.xIonDense[ipHYDROGEN][0];
	struc.hiistr[nzone_minus_1] = dense.xIonDense[ipHYDROGEN][1];
	struc.ednstr[nzone_minus_1] = (realnum)dense.eden;
	struc.o3str[nzone_minus_1] = dense.xIonDense[ipOXYGEN][2];
	struc.pressure[nzone_minus_1] = (realnum)pressure.PresTotlCurr;
	struc.windv[nzone_minus_1] = (realnum)wind.windv;
	struc.AccelTotalOutward[nzone_minus_1] = wind.AccelTotalOutward;
	struc.AccelGravity[nzone_minus_1] = wind.AccelGravity;
	struc.pres_radiation_lines_curr[nzone_minus_1] = (realnum)pressure.pres_radiation_lines_curr;
	struc.GasPressure[nzone_minus_1] = (realnum)pressure.PresGasCurr;
	struc.depth[nzone_minus_1] = (realnum)radius.depth;
	/* save absorption optical depth from illuminated face to current position */
	struc.xLyman_depth[nzone_minus_1] = opac.TauAbsFace[iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon];
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
	{
		struc.gas_phase[nzone_minus_1][nelem] = dense.gas_phase[nelem];
		for( ion=0; ion<nelem+2; ++ion )
		{
			struc.xIonDense[nzone_minus_1][nelem][ion] = dense.xIonDense[nelem][ion];
		}
	}
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( nelem=ipISO; nelem<LIMELM; ++nelem)
		{
			if( dense.lgElmtOn[nelem] )
			{
				for( long level=0; level < iso_sp[ipISO][nelem].numLevels_max; ++level )
				{
					struc.StatesElem[nzone_minus_1][nelem][nelem-ipISO][level] = (realnum)iso_sp[ipISO][nelem].st[level].Pop();
				}
			}
		}
	}

	/* the hydrogen molecules */
	for(mol=0;mol<mole_global.num_calc;mol++) 
	{
		struc.molecules[nzone_minus_1][mol] = (realnum) mole.species[mol].den;
	}
	struc.H2_abund[nzone_minus_1] = hmi.H2_total;

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			" radius_increment called; radius=%10.3e rinner=%10.3e DRAD=%10.3e drNext=%10.3e ROUTER=%10.3e DEPTH=%10.3e\n", 
		  radius.Radius, radius.rinner, radius.drad, radius.drNext, 
		  iterations.StopThickness[iteration-1], radius.depth );
	}

	/* remember mean and largest errors on electron density */
	Error = fabs(dense.eden - dense.EdenTrue)/SDIV(dense.EdenTrue);
	if( Error > conv.BigEdenError )
	{
		conv.BigEdenError = (realnum)Error;
		dense.nzEdenBad = nzone;
	}
	conv.AverEdenError += (realnum)Error;

	dense.EdenMax = MAX2( dense.eden , dense.EdenMax);
	dense.EdenMin = MIN2( dense.eden , dense.EdenMin);

	/* remember mean and largest errors between heating and cooling */
	Error = fabs(thermal.ctot - thermal.htot) / thermal.ctot;
	conv.BigHeatCoolError = MAX2((realnum)Error , conv.BigHeatCoolError );
	conv.AverHeatCoolError += (realnum)Error;

	/* remember mean and largest pressure errors */
	Error = fabs(pressure.PresTotlError);
	conv.BigPressError = MAX2((realnum)Error , conv.BigPressError );
	conv.AverPressError += (realnum)Error;

	/* integrate total mass over model (relative to 4pi rinner^2) */
	dense.xMassTotal += (realnum)(dense.xMassDensity * radius.dVeffVol);

	/* check cooling time for this zone, remember longest */
	timesc.time_therm_long = MAX2(timesc.time_therm_long,1.5*dense.pden*BOLTZMANN*phycon.te/
	  thermal.ctot);

	/* H 21 cm equilibrium timescale, H21cm returns H (not e) collisional 
	 * deexcitation rate (not cs) */
	t = H21cm_H_atom( phycon.te )* dense.xIonDense[ipHYDROGEN][0] +
		/* >>chng 02 feb 14, add electron term as per discussion in */
		/* >>refer	H1	21cm	Liszt, H., 2001, A&A, 371, 698 */
		H21cm_electron( phycon.te )*dense.eden;

	/* only update time scale if t is significant */
	if( t > SMALLFLOAT )
		timesc.TimeH21cm = MAX2( 1./t, timesc.TimeH21cm );

	/* remember longest CO timescale */
	if( (double)dense.xIonDense[ipCARBON][0]*(double)dense.xIonDense[ipOXYGEN][0] > SMALLFLOAT )
	{
		int ipCO = findspecies("CO")->index;
		/* this is rate CO is destroyed, equal to formation rate in equilibrium */
		if (ipCO != -1)
			timesc.BigCOMoleForm = MAX2( timesc.BigCOMoleForm, 1./SDIV(mole.species[ipCO].snk ));
	}

	/* remember longest H2  destruction timescale timescale */
	timesc.time_H2_Dest_longest = MAX2(timesc.time_H2_Dest_longest, timesc.time_H2_Dest_here );

	/* remember longest H2 formation timescale timescale */
	timesc.time_H2_Form_longest = MAX2( timesc.time_H2_Form_longest , timesc.time_H2_Form_here );

	/* increment counter if this zone possibly thermally unstable
	 * this flag was set in conv_temp_eden_ioniz.cpp, 
	 * derivative of heating and cooling negative */
	if( thermal.lgUnstable )
		thermal.nUnstable += 1;

	/* remember Stromgren radius - where hydrogen ionization falls below half */
	if( !rfield.lgUSphON && dense.xIonDense[ipHYDROGEN][0]/dense.gas_phase[ipHYDROGEN] > 0.49 )
	{
		rfield.rstrom = (realnum)radius.Radius;
		rfield.lgUSphON = true;
	}

	/* remember the largest value */
	wind.AccelMax = (realnum)MAX2(wind.AccelMax,wind.AccelTotalOutward);

	/* keep track of average acceleration */
	wind.AccelAver += wind.AccelTotalOutward*(realnum)radius.drad_x_fillfac;
	wind.acldr += (realnum)radius.drad_x_fillfac;

	/* following is integral of radiative force */
	pressure.pinzon = (realnum)(wind.AccelTotalOutward*dense.xMassDensity*
		radius.drad_x_fillfac*geometry.DirectionalCosin);
	/*fprintf(ioQQQ," debuggg pinzon %.2f %.2e %.2e %.2e\n", 
		fnzone,pressure.pinzon,dense.xMassDensity,wind.AccelTotalOutward);*/
	pressure.PresInteg += pressure.pinzon;

	// the integrated acceleration due to electron scattering, neglecting
	// absorption
	static realnum AccelElecScatZone1;
	if( nzone == 1 )
		AccelElecScatZone1 = wind.AccelElectron;
	pressure.pinzon_PresIntegElecThin = (realnum)(AccelElecScatZone1*dense.xMassDensity/radius.r1r0sq*
		radius.drad_x_fillfac*geometry.DirectionalCosin);
	pressure.PresIntegElecThin += pressure.pinzon_PresIntegElecThin;

	/* integrate gravitational pressure term */
	GravitationalPressure();
	pressure.IntegRhoGravity += pressure.RhoGravity;

	/* sound is sound travel time, sqrt term is sound speed */
	timesc.sound_speed_isothermal = sqrt(pressure.PresGasCurr/dense.xMassDensity);
	/* adiabatic sound speed assuming mono-atomic gas - gamma is 5/3*/
	timesc.sound_speed_adiabatic = sqrt(5.*pressure.PresGasCurr/(3.*dense.xMassDensity) );
	timesc.sound += radius.drad_x_fillfac / timesc.sound_speed_isothermal;

	/* save largest relative change in heating or cooling between this 
	 * iteration and previous iteration at this zone
	 * may be used to set time step in time dependent sims 
	 * nzonePreviousIteration is number of zones in previous iteration,
	 * 1 if only 1 done, while nzone_minus_1 is 1 on first zone */
	if( iteration > 1 && nzone_minus_1 < struc.nzonePreviousIteration )
	{
		/* set largest relative changes in heating/cooling between current
		 * and previous zones */
		realnum TempChange = (realnum)
			fabs( (phycon.te-struc.testr[nzone_minus_1])/phycon.te);

		struc.TempChangeMax = MAX2( struc.TempChangeMax , TempChange );
	}
	else
	{
		/* zero out on first iteration */
		struc.TempChangeMax = 0.;
	}
	/*fprintf(ioQQQ,"DEBUG radius_increment iteration %li Heat %.2e Cool %.2e change max \n",
		iteration , struc.HeatChangeMax , struc.CoolChangeMax);*/

	colden.dlnenp += dense.eden*(double)(dense.xIonDense[ipHYDROGEN][1])*radius.drad_x_fillfac;
	colden.dlnenHep += dense.eden*(double)(dense.xIonDense[ipHELIUM][1])*radius.drad_x_fillfac;
	colden.dlnenHepp += dense.eden*(double)(dense.xIonDense[ipHELIUM][2])*radius.drad_x_fillfac;
	colden.dlnenCp += dense.eden*(double)(dense.xIonDense[ipCARBON][1])*radius.drad_x_fillfac;

	// column densities of all states 
	for( unsigned i = 0; i < mole.species.size(); ++i )
	{
		if( mole.species[i].levels != NULL )
		{
			for( qList::iterator st = mole.species[i].levels->begin(); st != mole.species[i].levels->end(); ++st )
			{
				(*st).ColDen() += radius.drad_x_fillfac * (*st).Pop();
			}
		}
	}

	/* integral of n(H0) / Tspin - related to 21 cm optical depth*/
	if( hyperfine.Tspin21cm > SMALLFLOAT )
		colden.H0_ov_Tspin += (double)(dense.xIonDense[ipHYDROGEN][0]) / hyperfine.Tspin21cm*radius.drad_x_fillfac;

	/* >>chng 05 Mar 07, add integral of n(OH) / Tspin */
	colden.OH_ov_Tspin += (double)(findspecieslocal("OH")->den) / phycon.te*radius.drad_x_fillfac;

	/* this is Lya excitation temperature */
	hydro.TexcLya = (realnum)TexcLine( iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s) );

	/* count number of times Lya excitation temp hotter than gas */
	if( hydro.TexcLya > phycon.te )
	{
		hydro.nLyaHot += 1;
		if( hydro.TexcLya > hydro.TLyaMax )
		{
			hydro.TLyaMax = hydro.TexcLya;
			hydro.TeLyaMax = (realnum)phycon.te;
			hydro.nZTLaMax = nzone;
		}
	}

	/* column densities in various species */
	colden.colden[ipCOL_HTOT] += (realnum)(dense.gas_phase[ipHYDROGEN]*radius.drad_x_fillfac);
	ASSERT( colden.colden[ipCOL_HTOT] >SMALLFLOAT );
	/* >>chng 02 sep 20, from htwo to H2_total */
	/* >>chng 05 mar 14, rather than H2_total, give H2g and H2s */
	/* this is a special form of column density - should be proportional to total shielding */
	colden.coldenH2_ov_vel += hmi.H2_total*(realnum)radius.drad_x_fillfac / GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN]);

	colden.colden[ipCOL_elec] += (realnum)(dense.eden*radius.drad_x_fillfac);

	/* the ortho and para column densities */
	h2.ortho_colden += h2.ortho_density*radius.drad_x_fillfac;
	h2.para_colden += h2.para_density*radius.drad_x_fillfac;
	if( hmi.H2_total > SMALLFLOAT )
		ASSERT( fabs( h2.ortho_density + h2.para_density - hmi.H2_total ) / hmi.H2_total < 1e-4 );
	/*fprintf(ioQQQ,"DEBUG ortho para\t%.3e\t%.3e\ttot\t%.3e\t or pa colden\t%.3e\t%.3e\n", 
		h2.ortho_density, h2.para_density,hmi.H2_total,
		h2.ortho_colden , h2.para_colden);*/

	/*>>chng 27mar, GS, Column density of F=0 and F=1 levels of H0*/
	colden.H0_21cm_upper += ((*HFLines[0].Hi()).Pop()*radius.drad_x_fillfac);
	colden.H0_21cm_lower += ((*HFLines[0].Lo()).Pop()*radius.drad_x_fillfac);
	if (0)
		fprintf(ioQQQ,"DEBUG %g %g %g %g\n",colden.H0_21cm_upper,colden.H0_21cm_lower,
				  (*HFLines[0].Hi()).ColDen(),(*HFLines[0].Lo()).ColDen());
	/*fprintf(ioQQQ,"DEBUG pophi-poplo\t%.3e\t%.3e\radius\t%.3e\t col_hi\t%.3e\t%.3e\n", 
		HFLines[0].PopHi, HFLines[0].PopLo, radius.drad_x_fillfac,
		 HFLines[0].PopHi*radius.drad_x_fillfac,colden.H0_21cm_upper );*/

	/* now add total molecular column densities */
	molcol("ADD ",ioQQQ);

	/* increment forming the mean ionization and temperature */
	mean.MeanInc();

	/*-----------------------------------------------------------------------*/

	/* calculate average atomic weight per hydrogen of the plasma */
	avWeight = 0.;
	for( nelem=0; nelem < LIMELM; nelem++ ) 
	{
		avWeight += dense.gas_phase[nelem]*dense.AtomicWeight[nelem];
	}
	avWeight /= dense.gas_phase[ipHYDROGEN];

	/* compute some average grain properties */
	rfield.opac_mag_B_point = 0.;
	rfield.opac_mag_V_point = 0.;
	rfield.opac_mag_B_extended = 0.;
	rfield.opac_mag_V_extended = 0.;
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		/* this is total extinction in magnitudes at V and B, for a point source 
		 * total absorption and scattering,
		 * does not discount forward scattering to be similar to stellar extinction
		 * measurements made within ism */
		rfield.opac_mag_B_point += (gv.bin[nd].dstab1[rfield.ipB_filter-1] +
			gv.bin[nd].pure_sc1[rfield.ipB_filter-1])*double(gv.bin[nd].dstAbund)*
			double(dense.gas_phase[ipHYDROGEN]) * OPTDEP2EXTIN;

		rfield.opac_mag_V_point += (gv.bin[nd].dstab1[rfield.ipV_filter-1] +
			gv.bin[nd].pure_sc1[rfield.ipV_filter-1])*double(gv.bin[nd].dstAbund)*
			double(dense.gas_phase[ipHYDROGEN]) * OPTDEP2EXTIN;

		/* this is total extinction in magnitudes at V and B, for an extended source 
		 * total absorption and scattering,
		 * DOES discount forward scattering to apply for extended source like Orion */
		rfield.opac_mag_B_extended += (gv.bin[nd].dstab1[rfield.ipB_filter-1] +
			gv.bin[nd].pure_sc1[rfield.ipB_filter-1]*gv.bin[nd].asym[rfield.ipB_filter-1])*
			double(gv.bin[nd].dstAbund)*double(dense.gas_phase[ipHYDROGEN]) * OPTDEP2EXTIN;

		rfield.opac_mag_V_extended += (gv.bin[nd].dstab1[rfield.ipV_filter-1] +
			gv.bin[nd].pure_sc1[rfield.ipV_filter-1]*gv.bin[nd].asym[rfield.ipV_filter-1])*
			double(gv.bin[nd].dstAbund)*double(dense.gas_phase[ipHYDROGEN]) * OPTDEP2EXTIN;

		gv.bin[nd].avdust += gv.bin[nd].tedust*(realnum)radius.drad_x_fillfac;
		gv.bin[nd].avdft += gv.bin[nd].DustDftVel*(realnum)radius.drad_x_fillfac;
		gv.bin[nd].avdpot += (realnum)(gv.bin[nd].dstpot*EVRYD*radius.drad_x_fillfac);
		gv.bin[nd].avDGRatio += (realnum)(gv.bin[nd].dustp[1]*gv.bin[nd].dustp[2]*
			gv.bin[nd].dustp[3]*gv.bin[nd].dustp[4]*gv.bin[nd].dstAbund/avWeight*
			radius.drad_x_fillfac);
	}

	// doing the update outside the loop not only safes a few cycles, but also
	// minimizes the chance of the update getting lost in the numerical precision
	// that could lead to the sim never hitting A_V to go, and getting indefinitely
	// stuck at an epsilon distance before the requested A_V instead...
	rfield.extin_mag_B_point += rfield.opac_mag_B_point*radius.drad_x_fillfac;
	rfield.extin_mag_V_point += rfield.opac_mag_V_point*radius.drad_x_fillfac;
	rfield.extin_mag_B_extended += rfield.opac_mag_B_extended*radius.drad_x_fillfac;
	rfield.extin_mag_V_extended += rfield.opac_mag_V_extended*radius.drad_x_fillfac;

	/* there are some quantities needed to calculation the Jeans mass and radius */
	colden.TotMassColl += dense.xMassDensity*(realnum)radius.drad_x_fillfac;
	colden.tmas += (realnum)phycon.te*dense.xMassDensity*(realnum)radius.drad_x_fillfac;
	colden.wmas += dense.wmole*dense.xMassDensity*(realnum)radius.drad_x_fillfac;

	/* now find minimum Jeans length and mass; length in cm */
	double meanDensity = double(dense.xMassDensity)*geometry.FillFac;
	double rjeans = (log10(JEANS)+phycon.alogte - log10(dense.wmole) - 
						  log10(meanDensity))/2.;

	/* minimum Jeans mass in gm */
	// 0.30103 is log10(2.)
	ajmass = 3.*(rjeans - 0.30103) + log10(4.*PI/3.*meanDensity);

	/* now remember smallest */
	colden.rjnmin = MIN2(colden.rjnmin,(realnum)rjeans);
	colden.ajmmin = MIN2(colden.ajmmin,(realnum)ajmass);

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " radius_increment returns\n" );
	}
	return;
}
