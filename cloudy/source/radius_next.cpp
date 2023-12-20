/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*radius_next use adaptive logic to find next zone thickness */
/*ContRate called by radius_next to find energy of maximum continuum-gas interaction */
/*GrainRateDr called by radius_next to find grain heating rate dr */
/** \todo	2	- this routine is very important since it sets the pace for the calculation,
 * and directly affects the convergence of the code.  Most of the logic is very old and
 * messy.  
 * 1) make sure all test cases have save dr
 * 2) cat all these reasons together into one file and sort on the reason
 * 3) discover what logic is the main pacesetter for the code
 * 4) which are never triggered and so can be removed
 */
#include "cddefines.h"
#include "iso.h"
#include "geometry.h"
#include "h2.h"
#include "hyperfine.h"
#include "opacity.h"
#include "dense.h"
#include "heavy.h"
#include "grainvar.h"
#include "elementnames.h"
#include "dynamics.h"
#include "thermal.h"
#include "hmi.h"
#include "coolheavy.h"
#include "timesc.h"
#include "stopcalc.h"
#include "colden.h"
#include "phycon.h"
#include "rt.h"
#include "trace.h"
#include "wind.h"
#include "save.h"
#include "pressure.h"
#include "iterations.h"
#include "struc.h"
#include "radius.h"
#include "dark_matter.h"
#include "mole.h"
#include "rfield.h"
#include "freebound.h"
#include "lines_service.h"
#include "cosmology.h"

/*ContRate called by radius_next to find energy of maximum continuum-gas interaction */
STATIC void ContRate(double *freqm, 
		     double *opacm);

/*GrainRateDr called by radius_next to find grain heating rate dr */
STATIC void GrainRateDr(double *grfreqm, 
			double *gropacm);

class drChoiceItem
{
	string m_str;
	bool *m_flag;
public:
	drChoiceItem(const string& str, bool *flag) : m_str(str), m_flag(flag) {}
	void select() const
	{
		if (m_flag)
			*m_flag = true;
	}
	const char* c_str() const
	{
		return m_str.c_str();
	}
};

class drList
{
	multimap<double,drChoiceItem> m_map;
public:
	typedef multimap<double,drChoiceItem>::const_iterator const_iterator;
	void insert(double val, const string& str, bool* flag = NULL)
	{
		if (flag != NULL)
			*flag = false;
		m_map.insert(pair<const double,drChoiceItem>(val,drChoiceItem(str,flag)));
	}
	void clear()
		{
		m_map.clear();
	}
	const_iterator begin() const
	{
		return m_map.begin();
	}
	const_iterator end() const
	{
		return m_map.end();
	}
	size_t size() const
	{
		return m_map.size();
	}
};


/*radius_next use adaptive logic to find next zone thickness */
void radius_next()
{
	static double OHIIoHI, 
	  OldHeat = 0., 
	  OldTe = 0.,
	  OlddTdStep = 0.,
	  OldWindVelocity = 0.,
	  Old_H2_heat_cool;

	const double BigRadius = 1e30;

	DEBUG_ENTRY( "radius_next()" );

	/*--------------------------------------------------------------
	 *
	 * this sub determines the thickness of the next zone
	 * if is called one time for each zone
	 *
	 *-------------------------------------------------------------- */

	if( trace.lgTrace )
		fprintf( ioQQQ, "   radius_next called\n" );

	/* >>chng 03 sep 21 - decide whether this is the first call */
	bool lgFirstCall = ( nzone == 0 );

	drList drChoice;

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	 * check on change in hydrogen ionization */
	if( !lgFirstCall && OHIIoHI > 1e-15 &&
	    (dense.xIonDense[ipHYDROGEN][1] > 0. && dense.xIonDense[ipHYDROGEN][0] > 0.) )
	{
		double atomic_frac = dense.xIonDense[ipHYDROGEN][1]/dense.xIonDense[ipHYDROGEN][0];
		double drHion;
		/* ratio of current HII/HI to old value - < 1 when becoming more neutral */
		/* this is now change in atomic fraction */
		double x = fabs(1. - atomic_frac/OHIIoHI);
		if( atomic_frac > 0.05 && atomic_frac < 0.9 )
		{
			/* >>chng 96 feb 23 from 0.7 to 0.3 because punched through H i-front
			 * >>chng 97 may 5, from 0.3 to 0.2 since punched through H i-front */
			/* >>chng 02 jun 05 from 0.2 to 0.05 poorly resolved i-front, also added two-branch logic*/
			drHion = radius.drad*MAX2( 0.2 , 0.05/MAX2(1e-10,x) );
		}
		else
		{
			/* >>chng 96 feb 23 from 0.7 to 0.3 because punched through H i-front
			 * >>chng 97 may 5, from 0.3 to 0.2 since punched through H i-front */
			drHion = radius.drad*MAX2( 0.2 , 0.2/MAX2(1e-10,x) );
		}
		double SaveOHIIoHI = OHIIoHI;
		OHIIoHI = dense.xIonDense[ipHYDROGEN][1]/dense.xIonDense[ipHYDROGEN][0];
		ostringstream oss;
		oss << "change in H ionization old=" << scientific << setprecision(3);
		oss << SaveOHIIoHI << " new=" << OHIIoHI;
		drChoice.insert( drHion, oss.str() );
	}
	else
	{
		if( dense.xIonDense[ipHYDROGEN][1] > 0. && dense.xIonDense[ipHYDROGEN][0] > 0. )
			OHIIoHI = dense.xIonDense[ipHYDROGEN][1]/dense.xIonDense[ipHYDROGEN][0];
		else
			OHIIoHI = 0.;
	}

	if( rt.dTauMase < -0.01 )
	{
		/* maser so powerful, must limit inc in tay
		 * >>chng 96 aug 08, from 0.3 to 0.1 due to large maser crash */
		double drHMase = radius.drad*MAX2(0.1,-0.2/rt.dTauMase);
		ostringstream oss;
		//
		// NB NB -- DON'T ALTER THIS STRING, setting rt.lgMaserSetDR below keys from this!
		// No longer true:
		// >>chng 13 oct 01 rjrw take logical pointer argument for
		// >>optional flag to be set if selected
		//
		oss << "H maser dTauMase=" << scientific << setprecision(2) << rt.dTauMase << " ";
		oss << rt.mas_species << " " << rt.mas_ion << " " << rt.mas_hi << " " << rt.mas_lo;
		drChoice.insert( drHMase, oss.str(), &rt.lgMaserSetDR );
	}

	/* check change in relative ionization - this is to insure a gradual approach to
	 * ionization fronts - were dr allowed to remain constant we would crash through
	 * a front in a single zone */
	double change_heavy_frac_big = -1.;
	double frac_change_heavy_frac_big = -1.;
	const realnum CHANGE_ION_HEAV = 0.2f;
	const realnum CHANGE_ION_HHE = 0.15f;
	if( nzone > 4 )
	{
		long ichange_heavy_nelem = -1L;
		long ichange_heavy_ion = -1L;
		double dr_change_heavy = BigRadius;
		for( int nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				realnum change;
				/* will only track ions with fractional abundances greater than this
				 * He and C are special cases because of special roles in PDRs */
				realnum frac_limit;
				if( nelem<=ipHELIUM || nelem==ipCARBON )
				{
					frac_limit = 1e-4f;
					change = CHANGE_ION_HHE;
				}
				else
				{
					/* struc.dr_ionfrac_limit set in init_coreload, is limit to fractional
					 * abundances of ions that are used in prt_comment to indicate oscillations
					 * or rapid change in ionization */
					frac_limit = struc.dr_ionfrac_limit/2.f;
					change = CHANGE_ION_HEAV;
				}

				for( int ion=0; ion<=nelem+1; ++ion )
				{
					realnum abundnew = dense.xIonDense[nelem][ion]/SDIV( dense.gas_phase[nelem]);
					realnum abundold = struc.xIonDense[nzone-3][nelem][ion]/SDIV( struc.gas_phase[nzone-3][nelem]);
					if( abundnew > frac_limit && abundold > frac_limit )
					{
						/* this test is not done when [nzone-n] <0 */
						realnum abundolder = struc.xIonDense[nzone-4][nelem][ion]/SDIV( struc.gas_phase[nzone-4][nelem]);
						realnum abundoldest = struc.xIonDense[nzone-5][nelem][ion]/SDIV( struc.gas_phase[nzone-5][nelem]);
						/* this is fractional change, use min of old and new abundances, to pick up
						 * rapid changing Ca+ sooner */
						double change_heavy_frac = fabs(abundnew-abundold)/MIN2(abundold,abundnew);
						/* want fractional change to be less than this factor */
						if( (change_heavy_frac > change) && (change_heavy_frac > change_heavy_frac_big) &&
							/* >>chng 03 dec 07, add test that abund is not oscillating */
							/* also test that abundance is increasing - we are headed into a front */
							((abundnew-abundold)>0.)   && 
							((abundold-abundolder)>0.) && 
							((abundolder-abundoldest)>0.) )
						{
							ichange_heavy_nelem = nelem;
							ichange_heavy_ion = ion;
							change_heavy_frac_big = change_heavy_frac;
							frac_change_heavy_frac_big = abundnew;
							/* factor in max has been adjusted to prevent code from running
							 * through various ionization fronts.
							 * >>chng 03 dec 06, from min of 0.5 to min of 0.25, crash into He i-front
							 * in hizqso.in */
							/* >>chng 04 mar 03, min had become 0.1, forced oscillations in nova.in
							 * in silicon, zoning changed greatly, causing change in diffuse lin
							 * pumping.  put back to 0.25 */
							dr_change_heavy = radius.drad * MAX2(0.25, change / change_heavy_frac );
						}
					}
				}
			}
		}
		/* set to -1 before loop over ions, if no significant changes in ionization occurred
		 * then may still be -1
		 */
		ostringstream oss;
		if(ichange_heavy_nelem>=0)
		{
			oss << "change in ionization, element " << elementnames.chElementName[ichange_heavy_nelem];
			oss << " ion (C scale) " << ichange_heavy_ion << " rel change " << scientific << setprecision(2);
			oss << change_heavy_frac_big << " ion frac " << frac_change_heavy_frac_big;
		}
		else
		{
			oss << "none";
			dr_change_heavy = BigRadius;
		}
		drChoice.insert( dr_change_heavy, oss.str() );
	}

	/* if Leiden hacks are on then use increase in dust extinction as control 
	 * >>chng 05 aug 13, add this */
	if( mole_global.lgLeidenHack )
	{
		/* Draine field is only defined over narrow range in FUV - must not let change
		 * in extinction become too large - 
		 * prefactor is change in optical depth */
		double drLeiden_hack = MAX2( 0.02 , 0.05*rfield.extin_mag_V_point ) /
			SDIV( geometry.FillFac * rfield.opac_mag_V_point );
		drChoice.insert( drLeiden_hack, "Leiden hack" );
	}

	// limit to dust optical depth per zone
	static double extin_mag_V_point_old = -1.;
	if( nzone>1 )
	{
		 const double extin_mag_V_limit = 20.+0.01*extin_mag_V_point_old;
		 double dExtin = (rfield.extin_mag_V_point - extin_mag_V_point_old)/extin_mag_V_limit;
		 if( dExtin > SMALLFLOAT )
		 {
			  double drExtintion = radius.drad / dExtin;
				ostringstream oss;
				oss << "change in V mag extinction; current=" << scientific << setprecision(3) <<
					 rfield.extin_mag_V_point;
				oss << fixed << " delta=" << dExtin;
				drChoice.insert( drExtintion, oss.str() );
		 }

	}
	extin_mag_V_point_old = rfield.extin_mag_V_point;

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* check how heating is changing */
	if( nzone > 1 && !thermal.lgTemperatureConstant &&
	    /* if density fluctuations in place then override change in heat for dr set */
	    !( strcmp(dense.chDenseLaw,"SINE") == 0 && dense.lgDenFlucOn ) )
	{
		double dHdStep = safe_div(fabs(thermal.htot-OldHeat),thermal.htot,0.0);
		if( dHdStep > 0. )
		{
			double drdHdStep;
			if( dense.gas_phase[ipHYDROGEN] >= 1e13 )
			{
				drdHdStep = radius.drad*MAX2(0.8,0.075/dHdStep);
			}
			else if( dense.gas_phase[ipHYDROGEN] >= 1e11 )
			{
				drdHdStep = radius.drad*MAX2(0.8,0.100/dHdStep);
			}
			else
			{
				/* changed from .15 to .12 for outer edge of coolhii too steep dT
				 * changed to .10 for symmetry, big change in some rates, 95nov14
				 * changed from .10 to .125 since parispn seemed to waste zones
				 * >>chng 96 may 21, from .125 to .15 since pn's still waste zones */
				drdHdStep = radius.drad*MAX2(0.8,0.15/dHdStep);
			}
			ostringstream oss;
			oss << "change in heating; current=" << scientific << setprecision(3) << thermal.htot;
			oss << fixed << " delta=" << dHdStep;
			drChoice.insert( drdHdStep, oss.str() );
		}
	}
	OldHeat = thermal.htot;

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* pressure due to incident continuum if in eos */
	if( strcmp(dense.chDenseLaw,"CPRE") == 0 && pressure.lgContRadPresOn )
	{
		if( nzone > 2 && pressure.pinzon > 0. )
		{
			/* pinzon is pressrue from acceleration onto previos zone
			 * want this less than some fraction of total pressure */
			/* >>chng 06 feb 01, change from init pres to current total pressure
			 * in const press high U ulirgs current pressure may be quite larger
			 * than init pressure due to continuum absorption */
			double drConPres = 0.05*pressure.PresTotlCurr/
				(wind.AccelTotalOutward*dense.xMassDensity*geometry.FillFac);
			drChoice.insert( drConPres, "change in con accel" );
		}
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* gravitational pressure */
	if( strcmp(dense.chDenseLaw,"CPRE") == 0 && (dark.lgNFW_Set || pressure.gravity_symmetry>=0) )
	{
		if( fabs( pressure.RhoGravity ) > 0. )
		{
			double drGravPres = 0.05 * pressure.PresTotlCurr * radius.drad / fabs( pressure.RhoGravity );
			ASSERT( drGravPres > 0. );
			drChoice.insert( drGravPres, "change in gravitational pressure" );
		}
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* check how temperature is changing */
	double dTdStep = FLT_MAX;
	if( nzone > 1 )
	{
		/* change in temperature; current=	*/
		dTdStep = (phycon.te-OldTe)/phycon.te;
		/* >>chng 02 dec 08, desired change in temperature per zone must not
		 * be smaller than allower error in temperature.  For now use relative error
		 * in heating - - cooling balance.  Better would be to also use c-h deriv wrt t
		 * to get slope */
		/* >>chng 02 dec 11 rjrw change back to 0.03 -- improve dynamics.dRad criterion instead */
		double x = 0.03;
		/* >>chng 02 dec 07, do not do this if there is mild te jitter, 
		 * so that dT is changing sign - this happens
		 * in ism.ini, very stable temperature with slight noise up and down */
		if( dTdStep*OlddTdStep > 0. )
		{
			/* >>chng 96 may 30, variable depending on temp
			 * >>chng 96 may 31, allow 0.7 smaller, was 0.8
			 * >>chng 97 may 05, from 0.7 to 0.5 stop from punching through thermal front
			 * >>chng 04 mar 30, from 0.7 to 0.5 stop from punching through thermal front,
			 * for some reason factor was 0.7, not 0.5 as per previous change */
			double absdt = fabs(dTdStep);
			double drdTdStep = radius.drad*MAX2(0.5,x/absdt);
			ostringstream oss;
			oss << "change in temperature; current=" << scientific << setprecision(3) << phycon.te;
			oss << fixed << " dT/T=" << dTdStep;
			drChoice.insert( drdTdStep, oss.str() );
		}
	}
	OlddTdStep = dTdStep;
	OldTe = phycon.te;

	/* >>chng 02 oct 06, only check on opacity - interaction if not
	 * constant temperature - there were constant temperature models that
	 * extended to infinity but hung with last few photons and this test.
	 * better to ignore this test which is really for thermal models */
	/* >>chng 06 mar 20, do not call if recombination logic in place */
	if( !thermal.lgTemperatureConstant && !dynamics.lgRecom )
	{
		double freqm = 0., opacm = 0.;
		/* find freq where opacity largest and interaction rate is fastest */
		ContRate(&freqm,&opacm);

		/* use optical depth at max interaction energy
		* >>chng 96 jun 06 was drChange=0.15 changed to 0.3 for high Z models
		* taking too many zones
		* drInter = drChange / MAX(1e-30,opacm*FillFac) */

		double drInter = 0.3/MAX2(1e-30,opacm*geometry.FillFac*geometry.DirectionalCosin);
		ostringstream oss;
		oss << "cont inter nu=" << scientific << setprecision(2) << freqm << " opac=" << opacm;
		drChoice.insert( drInter, oss.str() );
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* check whether change in wind velocity constrains DRAD 
	 * WJH 22 May 2004: disable when we are near the sonic point since
	 * the velocity may be jumping all over the place but we just want
	 * to push through it as quickly as we can */
	if( !wind.lgStatic() && !pressure.lgSonicPoint && !pressure.lgStrongDLimbo )
	{
		double v = fabs(wind.windv);
		/* this is relative change in velocity */
		double dVelRelative = fabs(wind.windv-OldWindVelocity)/
			MAX2(v,0.1*timesc.sound_speed_isothermal);

		const double WIND_CHNG_VELOCITY_RELATIVE = 0.01;
		/*fprintf(ioQQQ,"DEBUG rad %.3f vel %.2e dVelRelative %.2e", 
			log10(radius.Radius) , wind.windv ,  dVelRelative );*/

		/* do not use this on first zone since do not have old velocity */
		double winddr;
		if( dVelRelative > WIND_CHNG_VELOCITY_RELATIVE  && nzone>1 )
		{
			/* factor less than one if change larger than WIND_CHNG_VELOCITY_RELATIVE */
			double factor = WIND_CHNG_VELOCITY_RELATIVE / dVelRelative;
			/*fprintf(ioQQQ," DEBUG factor %.2e", factor );*/
			factor = MAX2(0.8 , factor );
			winddr = radius.drad * factor;
		}
		else
		{
			winddr = BigRadius;
		}

		/* set dr from advective term */
		if( dynamics.lgAdvection && ! dynamics.isInitialRelaxIteration( iteration ) )
		{
			winddr = MIN2( winddr , dynamics.dRad );
			/*>>chng 04 oct 06, set dVelRelative to dynamics.dRad since dVelRelative is printed as part
			 * of reason for choosing this criteria, want to reflect valid reason. */
			dVelRelative = dynamics.dRad;
		}
		ostringstream oss;
		oss << "Wind, dVelRelative=" << scientific << setprecision(3);
		oss << dVelRelative << " windv=" << wind.windv;
		drChoice.insert( winddr, oss.str() );
	}
	/* remember old velocity */
	OldWindVelocity = wind.windv;

	const double DNGLOB = 0.10;

	/* inside out globule */
	if( strcmp(dense.chDenseLaw,"GLOB") == 0 )
	{
		if( radius.glbdst < 0. )
		{
			fprintf( ioQQQ, " Globule distance is negative, internal overflow has occured, sorry.\n" );
			fprintf( ioQQQ, " This is routine radius_next, GLBDST is%10.2e\n", 
			  radius.glbdst );
			cdEXIT(EXIT_FAILURE);
		}
		double factor = radius.glbden*pow(radius.glbrad/radius.glbdst,radius.glbpow);
		double fac2 = radius.glbden*pow(radius.glbrad/(radius.glbdst - (realnum)radius.drad),radius.glbpow);
		/* DNGLOB is relative change in density allowed this zone, 0.10 */
		double GlobDr = ( fac2/factor > 1. + DNGLOB ) ? radius.drad*DNGLOB/(fac2/factor - 1.) : BigRadius;
		GlobDr = MIN2(GlobDr,radius.glbdst/20.);
		ostringstream oss;
		oss << "GLOB law, HDEN=" << scientific << setprecision(2) << dense.gas_phase[ipHYDROGEN];
		drChoice.insert( GlobDr, oss.str() );
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	double hdnew = 0.;
	if( strncmp( dense.chDenseLaw , "DLW" , 3) == 0 )
	{
		/* one of the special density laws, first get density at possible next dr */
		if( strcmp(dense.chDenseLaw,"DLW1") == 0 )
		{
			hdnew = dense_fabden(radius.Radius+radius.drad,radius.depth+
			  radius.drad);
		}
		else if( strcmp(dense.chDenseLaw,"DLW2") == 0 )
		{
			hdnew = dense.DLW.tabval(radius.Radius+radius.drad,radius.depth+
			  radius.drad);
		}
		else if( strcmp(dense.chDenseLaw,"DLW3") == 0 )
		{
			hdnew = dense_parametric_wind(radius.Radius+radius.drad);
		}
		else
		{
			fprintf( ioQQQ, " dlw insanity in radius_next\n" );
			cdEXIT(EXIT_FAILURE);
		}
		double drTab = fabs(hdnew-dense.gas_phase[ipHYDROGEN])/MAX2(hdnew,dense.gas_phase[ipHYDROGEN]);
		drTab = radius.drad*MAX2(0.2,0.10/MAX2(0.01,drTab));
		ostringstream oss;
		oss << "spec den law, new old den " << scientific << setprecision(2);
		oss << hdnew << " " << dense.gas_phase[ipHYDROGEN];
		drChoice.insert( drTab, oss.str() );
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* special density law */
	if( strcmp(dense.chDenseLaw,"DLW1") == 0 )
	{
		double dnew = fabs(dense_fabden(radius.Radius+radius.drad,radius.depth+
		  radius.drad)/dense.gas_phase[ipHYDROGEN]-1.);
		/* DNGLOB is relative change in density allowed this zone, 0.10 */
		double SpecDr;
		if( dnew == 0. )
		{
			SpecDr = radius.drad*3.;
		}
		else if( dnew/DNGLOB > 1.0 )
		{
			SpecDr = radius.drad/(dnew/DNGLOB);
		}
		else
		{
			SpecDr = BigRadius;
		}
		ostringstream oss;
		oss << "special law, HDEN=" << scientific << setprecision(2) << dense.gas_phase[ipHYDROGEN];
		drChoice.insert( SpecDr, oss.str() );
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* check grain line heating dominates
	 * this is important in PDR and HII region calculations
	 * >>chng 97 jul 03, added following check */
	if( safe_div(thermal.heating(0,13),thermal.htot,0.0) > 0.2 )
	{
		double grfreqm = 0., gropacm = 0.;
		/* >>chng 01 jan 03, following returns 0 when NO light at all,
		 * had failed with botched monitor */
		GrainRateDr(&grfreqm,&gropacm);
		double DrGrainHeat = 1.0/MAX2(1e-30,gropacm*geometry.FillFac*geometry.DirectionalCosin);
		ostringstream oss;
		oss << "grain heating nu=" << scientific << setprecision(2) << grfreqm << " opac=" << gropacm;
		drChoice.insert( DrGrainHeat, oss.str() );
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* check if line heating dominates
	 * this is important in high metallicity models */
	if( safe_div(thermal.heating(0,22),thermal.htot,0.0) > 0.2 )
	{
		long level = -1;
		TransitionProxy t = FndLineHt(&level);

		if( t.Coll().heat()/thermal.htot > 0.1 )
		{
			ASSERT( t.associated() );

			double TauInwd = t.Emis().TauIn();
			double DopplerWidth = t.Emis().dampXvel()/t.Emis().damp();	
			double TauDTau = t.Emis().VoigtLineCen()*t.Emis().PopOpc()*t.Emis().opacity()/DopplerWidth;

			fixit("all other line stacks need to be included here.");
			// can we just sweep over line stack?  Is that ready yet?

			/* in following logic cannot use usual inverse opacity,
			 * since line heating competes with escape probability,
			 * so is effective at surprising optical depths */
			double drLineHeat = ( TauDTau > 0. ) ? MAX2(1.,TauInwd)*0.4/TauDTau : BigRadius;
			ostringstream oss;
			oss << "level " << level << " line heating, " << chLineLbl(t) << " TauIn " << scientific;
			oss << setprecision(2) << TauInwd << " " << t.Emis().pump() << " " << t.Emis().Pesc();
			drChoice.insert( drLineHeat, oss.str() );
		}
	}

	/* do not let change in cooling/heating due to H2 become too large */
	if( lgFirstCall )
		Old_H2_heat_cool = 0.;
	else if( !thermal.lgTemperatureConstant )
	{
		/* this is case where temperature has not been set */
		/* compare total heating - cooling due to h2 with total due to everything */
		double H2_heat_cool = safe_div((fabs(hmi.HeatH2Dexc_used)+fabs(hmi.HeatH2Dish_used)), thermal.htot,0.0);
		// Leiden hack PDR sims can have H2 heating set by sims equation, does not change
		double dH2_heat_cool = fabs( H2_heat_cool - Old_H2_heat_cool );
		if( H2_heat_cool > 0.1 && Old_H2_heat_cool>0. && dH2_heat_cool>SMALLFLOAT )
		{
			/* >>chng 05 oct 27, had been taking 20% of original radius - this caused zoning
			 * to become very fine and may not have been needed - change from 0.2 to 0.5 */
			double drH2_heat_cool = radius.drad*MAX2(0.5,0.05/dH2_heat_cool);
			ostringstream oss;
			oss << "change in H2 heating/cooling, d(c,h)/H " << scientific << setprecision(2);
			oss << dH2_heat_cool;
			drChoice.insert( drH2_heat_cool, oss.str() );
		}
	}
	Old_H2_heat_cool = safe_div((fabs(hmi.HeatH2Dexc_used)+fabs(hmi.HeatH2Dish_used)) , thermal.htot, 0.0);

	/* >>chng 03 mar 04, add this logic */
	/* do not let change in H2 and heavy element molecular abundances become too large 
	 * "change in heav ele mole abundance, d(mole)/elem" */
	if( nzone >= 4 )
	{
		/* first do H2 abundance */
		/* >>chng 04 dec 15, do not use special logic when large h2 is turned on */
		double Old_H2_abund = struc.H2_abund[nzone-3];
		/* >>chng 03 jun 16, limit from 0.01 to 0.001, some models fell over H2 front due to
		 * large zoning, when large H2 was just this caused oscillations in solomon process */
		/* >>chng 03 nov 22, from > 0.001 to > 3e-4, models that start almost in H2 have
		 * rapid increase in H2 at shallow depths, start sensing this sooner */
		/* >>chng 03 dec 10, from 3e-4 to 1e-4, to get smaller chagned in leiden1.in test */
		/* radius_next keys from change in H2 abundance, d(H2)/H */
		/* >>chng 04 mar 11, start sensing H2 earlier since otherwise step size
		 * needs to become way too small tooo quickly.  limit changed from 1e-4 to 1e-6 */
		/* >>chng 04 jun 29, fromo > 1e-6 to >1e-8, some pdr's had too large a change in H2 */

		if( 2.*hmi.H2_total/dense.gas_phase[ipHYDROGEN] > 1e-8 )
		{
			double fac = 0.2;
			/* this is percentage change in H2 density - "change in H2 abundance" */
			double dH2_abund = 2.*fabs( hmi.H2_total - Old_H2_abund ) / hmi.H2_total;
			/* in testing th85ism the dH2_abund did come out zero */
			/* >>chng 03 jun 16, change d(H2) from 0.05 to 0.1, fine resolution of H2/H exposed
			 * small numerical oscillations in Solomon process */
			/* >>chng 03 nov 22, smallest possible ratio of dr(next)/dr changed from
			 * 0.2 to 0.05, models that started almost in H2 front need to be able to sense it */
			/* >>chng 04 mar 15, with such small possible changes in dr, only 0.05, a thermal front
			 * can easily cause large changes in H2 and T that are not due to zoning, but to the
			 * discontinuity.  make smallest change larger to prevent hald due to dr */
			/* >>chng 05 aug 23, thermal front allowed dr to become much too small
			 * chng from 0.02 to .6 */
			dH2_abund = SDIV(dH2_abund);
			double drH2_abund = radius.drad*MAX2(0.2,fac/dH2_abund );
			ostringstream oss;
			oss << "change in H2 abundance, d(H2)/H " << scientific << setprecision(2) << dH2_abund;
			drChoice.insert( drH2_abund, oss.str() );
		}
	}

	int mole_dr_change = -1;
	if( nzone >= 4 )
	{	
		/* check how molecular fractions of all heavy elements are chaning relative 
		 * to total gas phase abundance */
		double max_change = 0.0;
		/* >>chng 04 jun 02, upper limit had been all species but now limit to real
		 * molecules since do not want this logic to work with the ions */
		for( int i=0; i < mole_global.num_calc; ++i )
		{
			realnum abund, 
				abund1,
				abund_limit;
			if( !mole_global.list[i]->isMonatomic() && mole_global.list[i]->isIsotopicTotalSpecies() )
			{
				/* >>chng 03 sep 21, add CO logic */
				/* >>chng 04 mar 30, generalize to any molecule at all */
				/* >>chng 04 mar 31 lower limit to abund had been 0.01, change
				 * to 0.001 to pick up approach to molecular fronts */
				abund = 0.;
				for( molecule::nNucsMap::iterator atom=mole_global.list[i]->nNuclide.begin();
					  atom != mole_global.list[i]->nNuclide.end(); ++atom)
				{
					long int nelem = atom->first->el()->Z-1;
					if (dense.lgElmtOn[nelem])
					{
						abund1 = (realnum)mole.species[i].den*atom->second/SDIV(dense.gas_phase[nelem]);
						if (abund1 > abund)
						{
							abund = abund1;
						}
					}
				}
				/* is this an ice?  need special logic for ice since density increases
				 * exponentially when grain temp falls below sublimation temperature 
				 * >>chng 05 dec 06 - detect changes for smaller abundances for ices
				 * due to large changes in ice abundances */
				if( mole_global.list[i]->lgGas_Phase )
				{
					abund_limit = 1e-3f;
				}
				else
				{
					/* this is an ice - track its abundance at lower values so that
					 * we resolve the sublimation transition region */
					abund_limit = 1e-5f;
				}

				if( abund > abund_limit )
				{
					/* the relative change in the abundance */
					double relative_change = 
						2.0 * fabs( mole.species[i].den - struc.molecules[nzone-3][i] ) /
						SDIV ( mole.species[i].den + struc.molecules[nzone-3][i] ) ;
					if (relative_change > max_change)
					{
						max_change = relative_change;
						mole_dr_change = i;
					}
				}
			}
		}
		if( mole_dr_change >= 0 )
		{
			double dr_mole_abund = radius.drad * MAX2(0.6, 0.05/SDIV(max_change));
			ostringstream oss;
			oss << "change in abundance species=" << mole_global.list[mole_dr_change]->label << ", ";
			oss << scientific << setprecision(2);
			oss << "old=" << struc.molecules[nzone-3][mole_dr_change] << " new=" << mole.species[mole_dr_change].den;
			drChoice.insert( dr_mole_abund, oss.str() );
		}
	}

	/* some consideration due to big H2 molecule */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		drChoice.insert( (*diatom)->H2_DR(), "change in big H2 Solomon rate line opt depth" );

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* can't make drmax large deep in model since causes feedback
	 * oscillations with changes in heating or destruction rates
	 * >>chng 96 oct 15, change from 5 to 10 */
	/* >>chng 96 nov 11, had been 4 * drad up to 11, change to following
	 * to be similar to c90.01, max of 1.3 between 5 and 11 
	 * >>chng 04 oct 29  geometry.DirectionalCosin was ioncorrect applied to this factor  */
	double drmax = ( nzone < 5 ) ? 4.*radius.drad : 1.3*radius.drad;
	drChoice.insert( drmax, "DRMAX" );

	/* time dependent recombination - conditions become very homogeneous and 
	 * crash into I fronts - use last iteration's  radius as guess of current case */
	double recom_dr_last_iter = BigRadius;
	if( dynamics.lgTimeDependentStatic && dynamics.lgRecom )
	{
		/* first time through nzone == 1 */
		static long int nzone_recom = -1;
		if( nzone <= 1 )
		{
			/* initialization case */
			nzone_recom = 0;
		}
		else if( radius.depth < struc.depth_last[struc.nzonePreviousIteration-1] && 
			nzone_recom < struc.nzonePreviousIteration )
		{
			ASSERT( nzone_recom>=0 && nzone_recom<struc.nzonePreviousIteration );
			/* case where we are within previous computed structure 
			 * first possibly increase nzone_recom so that it points deeper
			 * than this zone */
			while( struc.depth_last[nzone_recom] < radius.depth &&
				nzone_recom < struc.nzonePreviousIteration-1 )
				++nzone_recom;
			recom_dr_last_iter = struc.drad_last[nzone_recom]*3.;
			drChoice.insert( recom_dr_last_iter,
								    "previous iteration recom logic" );
		}
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* change in electron density - radius_next keys from change in elec den,
	 * remember old electron density during first call */
	/* this is low ionization solution */
	if( nzone > 2 )
	{
		/* next is-2 since nzone is on physics not c scale, and we want previous zone */
		double Efrac_old = struc.ednstr[nzone-3]/struc.gas_phase[nzone-3][ipHYDROGEN];
		double Efrac_new = dense.eden/dense.gas_phase[ipHYDROGEN];
		double dEfrac = fabs(Efrac_old-Efrac_new) / Efrac_new;
		double drEfrac;
		if( dEfrac > SMALLFLOAT )
		{
			double fac = 0.04;
			/* >>chng 03 dec 25, use smaller rel change in elec frac when most elec in ipMole or grains */
			/* >>chng 04 sep 14, change to from from metals but comment out */
			/* >>chng 04 sep 17, change to from from metals - uncomment */
			if( dense.eden_from_metals > 0.5 )
			{
				fac = 0.04;
			}
			/* >>chng 04 feb 23, add test on hydrogen being predomintly
			 * recombined due to three-body recom, which is very sensitive
			 * to the electron density - but only do this in partially ionized medium */
			else if( iso_sp[ipH_LIKE][ipHYDROGEN].RecomCollisFrac > 0.8 && 
				dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN] > 0.1 &&
				dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN] < 0.8 )

			{
				fac = 0.02;
			}
			/* >>chng 04 sep 17, change to 0.1 from 0.2 */
			drEfrac = radius.drad*MAX2(0.1,fac/dEfrac);
		}
		else
		{
			drEfrac = BigRadius;
		}
		ostringstream oss;
		oss << "change in elec den, rel chng: " << scientific << setprecision(3) << dEfrac;
		oss << " cur=" << Efrac_new << " old=" << Efrac_old;
		drChoice.insert( drEfrac, oss.str() );
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* do not let thickness get too large */
	if( nzone > 20 )
	{
		/* >>chng 02 nov 05, change from 1/20 to 1/10 wasted zones early on */
		drChoice.insert( radius.depth/10., "relative depth" );
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* case where stopping thickness or edge specified, need to approach slowly */
	double thickness_total = BigRadius;
	double DepthToGo = BigRadius;
	if( StopCalc.HColStop < 5e29 )
	{
		double coleff = SDIV( dense.gas_phase[ipHYDROGEN]*geometry.FillFac);
		DepthToGo = MIN2(DepthToGo ,
			(StopCalc.HColStop-colden.colden[ipCOL_HTOT]) / coleff );
		/* >>chng 03 oct 28, forgot to div col den above by eff density */
		thickness_total = MIN2(thickness_total , StopCalc.HColStop / coleff );
	}

	if( StopCalc.colpls < 5e29 )
	{
		double coleff = (double)SDIV( (dense.xIonDense[ipHYDROGEN][1])*geometry.FillFac);
		DepthToGo = MIN2(DepthToGo ,
			(StopCalc.colpls-findspecieslocal("H+")->column) / coleff );
		thickness_total = MIN2(thickness_total , StopCalc.colpls / coleff );
	}

	if( StopCalc.col_h2 < 5e29 )
	{
		/* >>chng 03 apr 15, add this molecular hydrogen */
		double coleff = (double)SDIV( hmi.H2_total*geometry.FillFac);
		DepthToGo = MIN2(DepthToGo ,
			(StopCalc.col_h2-findspecieslocal("H2")->column-findspecieslocal("H2*")->column) / coleff );
		thickness_total = MIN2(thickness_total , StopCalc.col_h2 / coleff );
	}

	if( StopCalc.col_h2_nut < 5e29 )
	{
		/* neutral column density of H, 2 H_2 + H^0 */
		double coleff = (double)SDIV( (2*hmi.H2_total+dense.xIonDense[ipHYDROGEN][0])*geometry.FillFac);
		DepthToGo = MIN2(DepthToGo ,
			(StopCalc.col_h2_nut-(2*(findspecieslocal("H2")->column+findspecieslocal("H2*")->column)+findspecieslocal("H")->column)) / coleff );
		thickness_total = MIN2(thickness_total , StopCalc.col_h2_nut / coleff );
	}

	if( StopCalc.col_H0_ov_Tspin < 5e29 )
	{
		/* >>chng 05 jan 09, add n(H0)/Tspin */
		double coleff = (double)SDIV( dense.xIonDense[ipHYDROGEN][0] / hyperfine.Tspin21cm*geometry.FillFac );
		DepthToGo = MIN2(DepthToGo ,
			(StopCalc.col_H0_ov_Tspin - colden.H0_ov_Tspin) / coleff  );
		thickness_total = MIN2(thickness_total , StopCalc.col_H0_ov_Tspin / coleff );
	}

	if( StopCalc.col_monoxco < 5e29 )
	{
		/* >>chng 03 apr 15, add this, CO */
		double coleff = (double)SDIV( (findspecieslocal("CO")->den)*geometry.FillFac);
		DepthToGo = MIN2(DepthToGo ,
			(StopCalc.col_monoxco-findspecieslocal("CO")->column) / coleff );
		thickness_total = MIN2(thickness_total , StopCalc.col_monoxco / coleff );
	}

	if( StopCalc.col_species < 5e29 )
	{
		/* arbitrary species */
		static molezone *SpeciesCurrent;
		static bool lgMustFind=true;
		if( lgMustFind )
		{
			lgMustFind = false;
			SpeciesCurrent= findspecieslocal_validate(StopCalc.chSpeciesColumn.c_str());
		}
		double coleff = (double)SDIV( (SpeciesCurrent->den)*geometry.FillFac);
		DepthToGo = MIN2(DepthToGo ,
			(StopCalc.col_species-SpeciesCurrent->column) / coleff );
		thickness_total = MIN2(thickness_total , StopCalc.col_species / coleff );
	}

	if( StopCalc.colnut < 5e29 )
	{
		double coleff = (double)SDIV( (dense.xIonDense[ipHYDROGEN][0])*geometry.FillFac);
		DepthToGo = MIN2(DepthToGo ,
			(StopCalc.colnut - findspecieslocal("H")->column) / coleff );
		thickness_total = MIN2(thickness_total , StopCalc.colnut / coleff );
	}

	/* this is case where outer radius is set */
	if( iterations.StopThickness[iteration-1] < 5e29 )
	{
		thickness_total = MIN2(thickness_total , iterations.StopThickness[iteration-1] );
		DepthToGo = MIN2(DepthToGo ,
			iterations.StopThickness[iteration-1] - radius.depth );
	}

	/* this is case where stopping optical depth was specified */
	if( StopCalc.iptnu != rfield.nflux_with_check )
	{
		/* end optical depth has been specified */
		double dt = SDIV(opac.opacity_abs[StopCalc.iptnu-1]*geometry.FillFac);
		DepthToGo = MIN2(DepthToGo ,
			(StopCalc.tauend-opac.TauAbsGeo[0][StopCalc.iptnu-1])/dt);
	}

	if( DepthToGo <= 0. )
		TotalInsanity();

	/* set dr if one of above tests have triggered */
	if( DepthToGo < BigRadius )
	{
		double depth_min = MIN2( DepthToGo , thickness_total/100. );
		DepthToGo = MAX2( DepthToGo / 10. , depth_min );
		/* want to approach the outer edge slowly,
		 * the need for this logic is most evident in brems.in - 
		 * HI fraction varies across coronal model */
		double drThickness = MIN2( thickness_total/10. , DepthToGo );
		drChoice.insert( drThickness, "depth to go" );
	}

	/* stop AV - usually this is dust, but we consider all opacity sources,
	 * so always include this */
	/* compute some average grain properties */
	double AV_to_go = BigRadius;
	if( rfield.opac_mag_V_extended > SMALLFLOAT && rfield.opac_mag_V_point > SMALLFLOAT )
	{
		double SAFETY = 1. + 10.*DBL_EPSILON;
		/* by default stop av is very large, and opacity can be very small, so ratio
		 * goes to inf - work with logs to see how big the number is */
		/* apply safety margin to avoid updates to the total extinction getting lost
		 * in the numerical precision */
		double ave = log10(SAFETY*StopCalc.AV_extended - rfield.extin_mag_V_extended) - 
			log10(rfield.opac_mag_V_extended);
		double avp = log10(SAFETY*StopCalc.AV_point - rfield.extin_mag_V_point) - 
			log10(rfield.opac_mag_V_point);
		AV_to_go = MIN2( ave , avp );
		if( AV_to_go < 37. )
		{
			AV_to_go = exp10( AV_to_go)/geometry.FillFac;
			/* this is to make sure that we go slightly deeper than AV so that
			 * we trigger this stop */
			AV_to_go *= 1.0001;
		}
		else
			AV_to_go = BigRadius;
		/*fprintf(ioQQQ,"DEBUG next dr %.3e %.3e %.3e\n", AV_to_go , ave, avp );*/
	}
	drChoice.insert( AV_to_go, "A_V to go" );

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* spherical models, do not want delta R/R big */
	double drSphere = radius.Radius*0.04;
	drChoice.insert( drSphere, "sphericity" );

	/* optical depth to electron scattering */
	double dRTaue = radius.drChange/(dense.eden*6.65e-25*geometry.FillFac);
	/* allow larger dr when constant temperature,
	 * to prevent some ct models from taking too many cells */
	if( thermal.lgTemperatureConstant ) 
		dRTaue *= 3.;
	drChoice.insert( dRTaue, "optical depth to electron scattering" );

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	if( dense.flong != 0. )
	{
		double drFluc = PI2/10./2./dense.flong;
		/* >>chng 04 sep 18, caused cautions that ionization jumps occurred.
		 * set to have the value */
		drFluc /= 2.;
		drChoice.insert( drFluc, "density fluctuations" );
	}

	/* keep dr constant in first two zones, if it wants to increase,
	 * but always allow it to decrease.
	 * to guard against large increases in efrac for balmer cont photo dominated models,
	 */
	double drStart = ( nzone <= 1 ) ? radius.drad : BigRadius;
	drChoice.insert( drStart, "capped to old DR in first zone" );

	// lgSdrmaxRel true if sdrmax is relative to current radius, false if limit in cm
	double rfacmax = radius.lgSdrmaxRel ? radius.Radius : 1.;
	drChoice.insert( rfacmax*radius.sdrmax, "sdrmax" );

	/* >>chng 05 aug 05, in case of thermal front, where temp is falling quickly and
	 * conditions change very fast, the zone thickness does not really affect the change
	 * in conditions and can cause zoning to become very very thin, which causes
	 * an abort.  this occurs between 200 and 1000K.  if we are doing temp soln,
	 * temp is between these values, and temp is changing rapidly, do not make sone
	 * thickness much smaller */
	/* do not use thermal front logic in case of recombination time dep cloud */
	if( nzone >= 5 && !dynamics.lgTimeDependentStatic )
	{
		/* >>chng 05 aug 23, upper bound of thermal from from 1000K to 4000K */
		/*if( phycon.te > 200. && phycon.te < 1000. && */
		if( phycon.te > 200. && phycon.te < 3000. && 
			/* >>chng 05 aug 23, from > 10% in zone to to two zones > 5%,
			 * to fix leiden v3 with large H2 */
			(struc.testr[nzone-3] - phycon.te)/phycon.te > 0.02 &&
			(struc.testr[nzone-4] - phycon.te)/phycon.te > 0.02 &&
			(struc.testr[nzone-5] - phycon.te)/phycon.te > 0.02 )
		{
			/* the 0.91 is to make dr unique, so that print statement that
			 * follows will identify this reason */
			double drThermalFront = radius.drad * 0.91;
			drChoice.clear();
			drChoice.insert( drThermalFront, "thermal front logic" );
		}
	}

	ASSERT( drChoice.size() > 0 );

	/* dr = zero is a logical mistake */
	if( drChoice.begin()->first <= 0. )
	{
		drList::const_iterator ptr = drChoice.begin();
		fprintf( ioQQQ, " radius_next finds insane drNext: %.2e\n", ptr->first );
		fprintf( ioQQQ, " all drs follow:\n" );
		while( ptr != drChoice.end() )
		{
			fprintf( ioQQQ, " %.2e %s\n", ptr->first, ptr->second.c_str() );
			++ptr;
		}
		cdEXIT(EXIT_FAILURE);
	}

	/* option to force min drad relative to depth */
	/* see google group thread "[cloudy-dev] sdrmin_rel_depth oscillation" ~2012/6/13 for discussion
	 * of validity of this test */
	if( !radius.lgFixed && drChoice.begin()->first < radius.depth*radius.sdrmin_rel_depth )
	{
		drChoice.clear();
		drChoice.insert( radius.depth*radius.sdrmin_rel_depth, "sdrmin_rel_depth" );
	}

	/* option to force min drad */
	double rfacmin = radius.lgSdrminRel ? radius.Radius : 1.;
	if( drChoice.begin()->first < rfacmin*radius.sdrmin )
	{
		drChoice.clear();
		drChoice.insert( rfacmin*radius.sdrmin, "sdrmin" );
	}

	/* this is general code that prevents zone thickness drNext from
	 * becoming too thin, something that can happen for various bad reasons
	 * HOWEVER we do not want to do this test for certain density laws,
	 * for which very small zone thicknesses are unavoidable
	 * the special cases are:
	 * special density law,
	 * globule density law,
	 * carbon +-0 i front
	 * the fluctuations command
	 * drMinimum was set in radius_first to either sdrmin (set drmin) or
	 * some fraction of the initial radius - it is always set
	 * to something non-trivial.  
	 * sdrmin is set with the "set dr" command - is SMALLFLOAT by default */
	if( strcmp(dense.chDenseLaw,"DLW1") != 0 && 
	    strcmp(dense.chDenseLaw ,"GLOB") != 0 && 
	    dense.flong == 0. &&
	    // stopping on depth to go is not a fault
	    drChoice.begin()->first != DepthToGo &&
		// stopping on Av to go is not a fault
		drChoice.begin()->first != AV_to_go )
	{
		fixit( "drMinimum does not get reevaluated in each iteration,"
			" so time-dependent solutions just use the first drad and Hden."
			" Re-eval here for now to keep this minimum from blowing." );
		radius.drMinimum = (realnum)(radius.drad * dense.gas_phase[ipHYDROGEN]/1e7);

		/* don't let dr get smaller than drMinimum, if this resets drNext
		 * then code stops with warning that zones got too thin
		 * this can happen due to numerical oscillations, which the code
		 * tries to damp out by making the zones thinner.
		 * scale by radius.Radius/radius.rinner to stop very spherical 
		 * simulations from false trigger on dr too small */
		/* drMinimum is drad * hden, to make it proportional to optical depth.
		 * This avoids false trigger across thermal fronts. */
		if( drChoice.begin()->first*radius.Radius/radius.rinner*dense.gas_phase[ipHYDROGEN] < radius.drMinimum )
		{
			/* must decrement nzone, since we will not complete this zone, and will not have
			 * valid structure data for it */
			--nzone;
			fprintf( ioQQQ, 
				"\n DISASTER PROBLEM radius_next finds dr too small and aborts.  "
				"This is zone %ld iteration %ld.  The thickness was %.2e\n", 
				nzone, 
				iteration,
				radius.drNext);
			fprintf( ioQQQ, 
				" If this simulation seems reasonable you can override this limit "
				"with the command SET DRMIN %.2e\n\n", 
				radius.drNext*2);
			throw cloudy_abort("radius_next finds dr too small");
		}
	}

	/* factor to allow for slop in floating numbers */
	const double Z = 1.0001;

	/* following is to make thickness of model exact
	 * n.b., on last zone, drNext can be NEGATIVE!!
	 * DEPTH was incremented at start of zone calc in newrad,
	 * has been outer edge of zone all throughout */
	double drOuterRadius = (iterations.StopThickness[iteration-1]-radius.depth)*Z;
	drChoice.insert( drOuterRadius, "outer radius" );

	if( cosmology.lgDo )
	{
		drChoice.clear();
		double drHz = SPEEDLIGHT / GetHubbleFactor( cosmology.redshift_current );
		drChoice.insert( drHz, "c over H(z)" );
	}

	// choose the smallest dR as the next choice
	radius.drNext = drChoice.begin()->first;

	/* set choice flag if provided -- in practice, just for rt.lgMaserSetDR */
	drChoice.begin()->second.select( );

	/* all this is to only save on last iteration
	 * the save dr command is not really a save command, making this necessary
	 * lgDRon is set true if "save dr" entered */
	/* lgDRPLst was set true if "save" had "last" on it */
	bool lgDoPun = ( save.lgDROn && !( save.lgDRPLst && !iterations.lgLastIt ) );
	if( (trace.lgTrace && trace.lgDrBug) || lgDoPun )
	{
		fprintf( save.ipDRout , "%ld\t%.5e\t%.3e\t%.3e\t%s\n", nzone, radius.depth,
			 radius.drNext, radius.Depth2Go, drChoice.begin()->second.c_str() );
	}

	ASSERT( radius.drNext > 0. );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " radius_next chooses next drad drNext=%.4e; this drad was %.4e\n", 
			 radius.drNext, radius.drad );
	}
}

/*ContRate called by radius_next to find energy of maximum continuum-gas interaction */
STATIC void ContRate(double *freqm, 
		     double *opacm)
{
	long int i, 
	  ipHi,
	  iplow, 
	  limit;
	double FreqH, 
	  Freq_nonIonizing, 
	  Opac_Hion, 
	  Opac_nonIonizing, 
	  Rate_max_Hion, 
	  Rate_max_nonIonizing;

	DEBUG_ENTRY( "ContRate()" );

	/* 
	 * find maximum continuum interaction rate,
	 * these should be reset in following logic without exception,
	 * if they are still zero at the end we have a logical error 
	 */
	Rate_max_nonIonizing = 0.;
	Freq_nonIonizing = 0.;
	Opac_nonIonizing = 0.;

	/* this must be reset to val >= 0 */
	*opacm = -1.;
	*freqm = -1.;

	/* do up to carbon photo edge if carbon is turned on */
	/* >>>chng 00 apr 07, add test for whether element is turned on */
	if( dense.lgElmtOn[ipCARBON] )
	{
		/* carbon is turned on, use carbon 1 edge */
		ipHi = Heavy.ipHeavy[ipCARBON][0] - 1;
	}
	else
	{
		/* carbon turned off, use hydrogen Balmer continuum */
		ipHi = iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH2s].ipIsoLevNIonCon-1;
	}

	/* start at plasma frequency */
	for( i=rfield.ipPlasma; i < ipHi; i++ )
	{
		/* this does not have grain opacity since grains totally passive
		 * at energies smaller than CI edge */
		if( rfield.anu(i)*rfield.flux[0][i]/rfield.widflx(i)*(opac.opacity_abs[i] - 
		  gv.dstab[i]*dense.gas_phase[ipHYDROGEN]) > Rate_max_nonIonizing )
		{
			Rate_max_nonIonizing = rfield.anu(i)*rfield.flux[0][i]/rfield.widflx(i)*
			  (opac.opacity_abs[i] - gv.dstab[i]*dense.gas_phase[ipHYDROGEN]);
			Freq_nonIonizing = rfield.anu(i);
			Opac_nonIonizing = opac.opacity_abs[i] - gv.dstab[i]*dense.gas_phase[ipHYDROGEN];
		}
	}

	/* not every continuum extends beyond C edge-this whole loop can add to zero
	 * use total opacity here
	 * test is to put in fir continuum if free free heating is important */
	if( safe_div(CoolHeavy.brems_heat_total,thermal.htot,0.0) > 0.05 )
	{
		/* this is index for energy where cloud free free optical depth is unity,
		 * is zero if no freq are opt thin */
		iplow = MAX2(1 , rfield.ipEnergyBremsThin );
	}
	else
	{
		/* >>>chng 00 apr 07, from Heavy.ipHeavy[0][5] to ipHi defined above, since
		 * would crash if element not defined */
		iplow = ipHi;
	}
	/* do not go below plasma frequency */
	iplow = MAX2( rfield.ipPlasma , iplow );

	/* this energy range from carbon edge to hydrogen edge */
	limit = MIN2(rfield.nflux,iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1);
	for( i=iplow; i < limit; i++ )
	{
		if( rfield.anu(i)*rfield.flux[0][i]/rfield.widflx(i)*(opac.opacity_abs[i] - 
		  gv.dstab[i]*dense.gas_phase[ipHYDROGEN]) > Rate_max_nonIonizing )
		{
			Rate_max_nonIonizing = rfield.anu(i)*rfield.flux[0][i]/rfield.widflx(i)*
			  (opac.opacity_abs[i] - gv.dstab[i]*dense.gas_phase[ipHYDROGEN]);
			Freq_nonIonizing = rfield.anu(i);
			Opac_nonIonizing = opac.opacity_abs[i] - gv.dstab[i]*dense.gas_phase[ipHYDROGEN];
		}
	}

	/* variables to check continuum interactions over Lyman continuum */
	Rate_max_Hion = 0.;
	Opac_Hion = 0.;
	FreqH = 0.;

	/* not every continuum extends beyond 1 Ryd-this whole loop can add to zero */
	for( i=iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon-1; i < rfield.nflux; i++ )
	{
		if( rfield.anu(i)*rfield.flux[0][i]/rfield.widflx(i)*(opac.opacity_abs[i] - 
		  gv.dstab[i]*dense.gas_phase[ipHYDROGEN]) > Rate_max_Hion )
		{
			/* Rate_max_Hion = anu(i)*flux(i)/widflx(i)*opac(i) */
			Rate_max_Hion = rfield.anu(i)*rfield.flux[0][i]/rfield.widflx(i)*
			  (opac.opacity_abs[i] - gv.dstab[i]*dense.gas_phase[ipHYDROGEN]);
			FreqH = rfield.anu(i);
			Opac_Hion = opac.opacity_abs[i] - gv.dstab[i]*dense.gas_phase[ipHYDROGEN];
		}
	}

	/* use Lyman continuum if its opacity is larger than non-h ion */
	if( Rate_max_nonIonizing < 1e-30  && Opac_Hion > SMALLFLOAT )
	{
		/* this happens for laser source - use Lyman continuum */
		*opacm = Opac_Hion;
		*freqm = FreqH;
	}
	/* >>chng 05 aug 03, add last test on Opac_Hion for case where we go very
	 * deep and very little radiation is left */
	else if( Opac_Hion > Opac_nonIonizing && Rate_max_Hion/SDIV(Rate_max_nonIonizing) > 1e-10 && Opac_Hion > SMALLFLOAT )
	{
		/* use Lyman continuum */
		*opacm = Opac_Hion;
		*freqm = FreqH;
	}
	else
	{
		/* not much rate in the Lyman continuum, stick with low energy */
		*opacm = Opac_nonIonizing;
		*freqm = Freq_nonIonizing;
	}

#	if 0
	/*>>chng 06 aug 12, do not let zones become optically thick to
	* peak of dust emission - one of NA's vastly optically thick dust clouds
	* did not conserve total energy - very opticallly thick to ir dust emission
	* so ir was absorbed and reemitted many times
	* this makes sure the cells remain optically thin to dust emission */
	if( gv.lgDustOn() && gv.lgGrainPhysicsOn )
	{
		double fluxGrainPeak = -1.;
		long int ipGrainPeak = -1;
		long int ipGrainPeak2 = -1;

		i = 0;
		while( rfield.anu(i) < 0.0912  && i < (rfield.nflux-2) )
		{
			/* >>chng 06 aug 23.  Only want to do this test for the IR dust continuum, therefore only 
			 * check on optical depth for wavelengths greater than 1 micron */
			if( gv.GrainEmission[i]*rfield.anu(i)*opac.opacity_abs[i] > fluxGrainPeak )
			{
				ipGrainPeak = i;
				fluxGrainPeak = gv.GrainEmission[i]*rfield.anu(i)*opac.opacity_abs[i];
			}
			++i;
		}
		ASSERT( fluxGrainPeak>=0. && ipGrainPeak >=0 );

		/* >>chng 06 aug 23.  We have just found the wavelength and flux at the peak of the IR emission.
		 * Now we want to find the wavelength, short of the peak, which corresponds to 5% of the 
		 * peak emission.  This wavelength will be where the code checks to make sure the zone has
		 * not become optically thick.  Since dust opacity generally decreases with wavelength,  this assures that
		 * the optical depths are small for all wavelengths where the flux is appreciable.  Tests show that
		 * this allows for flux/luminosity conservation to within ~5% for an AV of 1e4 mag and at least 2 iterations*/
		i = ipGrainPeak;
		/* >>chng 06 aug 23.  Only want to do this test for the IR dust continuum, therefore only 
		 * check on opacities for wavelengths shortward of the peak and greater than 1 micron 
		 * this routine can be called with the dust emission totally zero - it is only evaluated 
		 * late to save time.  don't do the following tests if peak dust emission is zero */
		if( fluxGrainPeak > 0. )
		{
			while( rfield.anu(i) < 0.0912 && i < (rfield.nflux-2) )
			{
				/* find wavelength where flux = 5% of peak flux, shortward of the peak */
				if( gv.GrainEmission[i]*rfield.anu(i)*opac.opacity_abs[i] > 0.05*fluxGrainPeak )
				{
					/* This will be the array number and flux at 10% of the peak */
					ipGrainPeak2 = i;
				}
				++i;
			}
			ASSERT( ipGrainPeak2>=0 );

			/* use this as a limit to the dr if opacity is larger */
			if( opac.opacity_abs[ipGrainPeak2] > *opacm )
			{
				/* scale opacity by factor of 10, which further decreases the zone thickness*/
				*opacm = opac.opacity_abs[ipGrainPeak2]*10.;
				*freqm = rfield.anu(ipGrainPeak2);
				/*fprintf(ioQQQ,"DEBUT opac peak %.2e %.2e \n",
					*opacm , *freqm );*/
			}
		}
	}
#	endif

	{
		/* following should be set true to print contributors */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"conratedebug \t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n", 
			Rate_max_nonIonizing,Freq_nonIonizing,Opac_nonIonizing,
			Rate_max_Hion,FreqH ,Opac_Hion,*freqm,*opacm
			);

		}
	}

	/* these were set to -1 at start, and must have been reset in one of the
	 * two loops.  Logic error if still <0. */
	/* >>chng 05 aug 03, change logic to -1 on entry and check at least zero
	 * here - will be zero if NO radiation field exists, perhaps since totally
	 * attenuated */
	ASSERT( *opacm >= 0. && *freqm >= 0. );
	return;
}

/*GrainRateDr called by radius_next to find grain heating rate dr */
STATIC void GrainRateDr(double *grfreqm, 
			double *gropacm)
{
	long int i, 
	  iplow;
	double xMax;

	DEBUG_ENTRY( "GrainRateDr()" );

	/* in all following changed from anu2 to anu  july 25 95
	 *
	 * find maximum continuum interaction rate */

	/* not every continuum extends beyond C edge-this whole loop can add to zero
	 * use total opacity here
	 * test is to put in fir continuum if free free heating is important */
	if( safe_div(CoolHeavy.brems_heat_total,thermal.htot,0.0) > 0.05 )
	{
		/* this is pointer to energy where cloud free free optical depth is unity,
		 * is zero if no freq are opt thin */
		iplow = MAX2(1 , rfield.ipEnergyBremsThin );
	}
	else
	{
		/* do up to carbon photo edge if carbon is turned on */
		/* >>>chng 00 apr 07, add test for whether element is turned on */
		if( dense.lgElmtOn[ipCARBON] )
		{
			/* carbon is turned on, use carbon 1 edge */
			iplow = Heavy.ipHeavy[ipCARBON][0];
		}
		else
		{
			/* carbon truned off, use hydrogen balmer continuum */
			iplow = iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH2s].ipIsoLevNIonCon;
		}
	}

	xMax = -1.;
	/* integrate up to H edge */
	for( i=iplow-1; i < Heavy.ipHeavy[ipHYDROGEN][0]; i++ )
	{
		if( rfield.anu(i)*rfield.flux[0][i]/rfield.widflx(i)*opac.opacity_abs[i] > xMax )
		{
			xMax = rfield.anu(i)*rfield.flux[0][i]/rfield.widflx(i)*
			  opac.opacity_abs[i];
			*grfreqm = rfield.anu(i);
			*gropacm = opac.opacity_abs[i];
		}
	}
	/* integrate up to heii edge if he is turned on,
	 * this logic will not make sense if grains on but he off, which in itself makes no sense*/
	if( dense.lgElmtOn[ipHELIUM] )
	{
		for( i=Heavy.ipHeavy[ipHYDROGEN][0]-1; i < Heavy.ipHeavy[ipHELIUM][1]; i++ )
		{
			if( rfield.anu(i)*rfield.flux[0][i]/rfield.widflx(i)*opac.opacity_abs[i] > xMax )
			{
				xMax = rfield.anu(i)*rfield.flux[0][i]/rfield.widflx(i)*
				  opac.opacity_abs[i];
				*grfreqm = rfield.anu(i);
				*gropacm = opac.opacity_abs[i];
			}
		}
	}

	/* possible that there is NO ionizing radiation, in extreme cases,
	 * if so then xMax will still be negative */
	if( xMax <= 0. )
	{
		*gropacm = 0.;
		*grfreqm = 0.;
	}
	return;
}
