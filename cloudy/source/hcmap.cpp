/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*map_do produce map of heating-cooling space for specified zone, called as result of
 * map command  */
#include "cddefines.h"
#include "thermal.h"
#include "cooling.h"
#include "called.h"
#include "dense.h"
#include "phycon.h"
#include "trace.h"
#include "pressure.h"
#include "conv.h"
#include "hcmap.h"
#include "ion_trim.h"
#include "hmi.h"

static const double EPS = 0.005;

t_hcmap hcmap;

void map_do(
			FILE *io, 
			const char *chType)
{
	char chLabel[NCOLNT_LAB_LEN+1];
	char units;
	long int i, 
	  ksav, 
	  j, 
	  jsav, 
	  k,
	  nelem;
	realnum wl;
	double cfrac, 
	  ch, 
	  fac, 
	  factor, 
	  hfrac, 
	  oldch, 
	  ratio, 
	  strhet, 
	  strong, 
	  t1, 
	  tlowst, 
	  tmax, 
	  tmin, 
	  torg;

	DEBUG_ENTRY( "map_do()" );

	t1 = phycon.te;
	torg = phycon.te;
	hcmap.lgMapOK = true;
	/* flag indicating that we have computed a map */
	hcmap.lgMapDone = true;

	/* make sure pressure has been evaluated */
	/* this sets values of pressure.PresTotlCurr */
	PresTotCurrent();

	/* print out all coolants if all else fails */
	if( called.lgTalk )
	{
		fprintf( io, "# Cloudy punts, Te=%10.3e HTOT=%10.3e CTOT=%10.3e nzone=%4ld\n", 
		  phycon.te, thermal.htot, thermal.ctot, nzone );
		fprintf( io, "# COOLNG array is\n" );

		if( thermal.ctot > 0. )
		{
			coolpr(io, "ZERO",1,0.,"ZERO");
			for( i=0; i < thermal.ncltot; i++ )
			{
				ratio = thermal.cooling[i]/thermal.ctot;
				if( ratio>EPS )
				{
					coolpr(io, (char*)thermal.chClntLab[i],thermal.collam[i],
					  ratio,"DOIT");
				}
			}

			fprintf( io, "#" );
			coolpr(io, "DONE",1,0.,"DONE");
			fprintf( io, "# Line heating array follows\n" );
			coolpr(io, "ZERO",1,0.,"ZERO");

			for( i=0; i < thermal.ncltot; i++ )
			{
				ratio = thermal.heatnt[i]/thermal.ctot;
				if( ratio>EPS )
				{
					coolpr(io, (char*)thermal.chClntLab[i],thermal.collam[i],
					  ratio,"DOIT");
				}
			}

			fprintf( io, "#" );
			coolpr(io,"DONE",1,0.,"DONE");
		}
	}

	/* map out te-ionization-cooling space before punching out. */
	if( called.lgTalk )
	{
		fprintf( io, "# map of heating, cooling, vs temp, follows.\n");
		fprintf( io, 
			"#    Te\t\t Heat-------------------------------->\tCool----------------------------------------->\t    dH/dT\t    dC/DT\t      Ne\t       NH\t       H2\t HII\t Helium \n" );
	}

	if( strcmp(chType,"punt") == 0 )
	{
		/* this is the original use of punt, we are punting
		 * only map within factor of two of final temperature
		 * fac will the range to either side of punted temperature */
		fac = 1.5;
		tmin = torg/fac;
		tmax = torg*fac;

		/* we want about 20 steps between high and low temperature
		 * default of nMapStep is 20, set with set nmaps command */
		factor = exp10(log10(tmax/tmin)/(double)(hcmap.nMapStep));
		TempChange(tmin , false);
	}

	else if( strcmp(chType," map") == 0 )
	{
		/* create some sort of map of heating-cooling */
		tlowst = MAX2(hcmap.RangeMap[0],phycon.TEMP_LIMIT_LOW);
		tmin = tlowst*0.998;
		tmax = MIN2(hcmap.RangeMap[1],phycon.TEMP_LIMIT_HIGH)*1.002;

		/* we want about nMapStep (=20) steps between high and low temperature */
		factor = exp10(log10(tmax/tmin)/(double)(hcmap.nMapStep));
		double TeNew;
		if( thermal.lgTeHigh )
		{
			/* high te */
			factor = 1./factor;
			/* TeHighest is highest possible temperature, 1E10 */
			TeNew = (MIN2(hcmap.RangeMap[1],phycon.TEMP_LIMIT_HIGH)/factor);
		}

		else
		{
			/* low te */
			TeNew = (tlowst/factor);
		}
		TempChange(TeNew , false);
	}

	else
	{
		/* don't know what to do */
		fprintf( ioQQQ, " PUNT called with insane argument,=%4.4s\n", 
		  chType );
		cdEXIT(EXIT_FAILURE);
	}

	/* now allocate space for te, c, h vectors in map, if not already done */
	if( hcmap.nMapAlloc==0 )
	{
		/* space not allocated, do so now */
		hcmap.nMapAlloc = hcmap.nMapStep+4;

		/* now make the space */
		hcmap.temap.resize(hcmap.nMapStep+4);
		hcmap.cmap.resize(hcmap.nMapStep+4);
		hcmap.hmap.resize(hcmap.nMapStep+4);
	}

	thermal.lgCNegChk = false;
	hcmap.nmap = 0;
	oldch = 0.;
	TempChange(phycon.te *factor , true);
	if( trace.nTrConvg )
		fprintf(ioQQQ, "    MAP called temp range %.4e %.4e in %li stops ===============================================\n",
		tmin,
		tmax,
		hcmap.nmap);

	while( (double)phycon.te < tmax*0.999 && (double)phycon.te > tmin*1.001 )
	{
		/* this sets values of pressure.PresTotlCurr */
		PresTotCurrent();

		/* must reset upper and lower bounds for ionization distributions */
		/* fix number of stages of ionization */
		for( nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				ion_trim_untrim(nelem);
			}
			else
			{
				/* this element is turned off, make stages impossible */
				ion_trim_invalidate(nelem);
			}
		}

		/* this turns on constant reevaluation of everything */
		conv.lgSearch = true;

		if( trace.nTrConvg )
			fprintf(ioQQQ, "    MAP new temp %.4e ===============================================\n",
			phycon.te );

		/* this counts how many times ionize is called in this model after startr,
		 * and is flag used by ionize to understand it is being called the first time*/
		conv.nTotalIoniz = 0;

		/* now get ionization solution for this temperature */
		ConvEdenIoniz();

		/* save results for later prints */
		hcmap.temap[hcmap.nmap] = phycon.te;
		hcmap.cmap[hcmap.nmap] = thermal.ctot;
		hcmap.hmap[hcmap.nmap] = thermal.htot;

		wl = 0.f;
		strong = 0.;

		for( j=0; j < thermal.ncltot; j++ )
		{
			if( thermal.cooling[j] > strong )
			{
				strcpy( chLabel, thermal.chClntLab[j] );
				strong = thermal.cooling[j];
				wl = thermal.collam[j];
			}
		}

		cfrac = strong/thermal.ctot;
		strhet = 0.;
		/* these will be reset in following loop*/
		ksav = -INT_MAX;
		jsav = -INT_MAX;

		for( k=0; k < LIMELM; k++ )
		{
			for( j=0; j < LIMELM; j++ )
			{
				if( thermal.heating(k,j) > strhet )
				{
					strhet = thermal.heating(k,j);
					jsav = j;
					ksav = k;
				}
			}
		}

		ch = thermal.ctot - thermal.htot;
		/* use ratio to check for change of sign since product
		 * can underflow at low densities */
		if( oldch/ch < 0. && called.lgTalk )
		{
			fprintf( io, "# ----------------------------------------------- Probable thermal solution here. --------------------------------------------\n" );
		}

		oldch = ch;
		hfrac = strhet/thermal.htot;
		if( called.lgTalk )
		{
			/* convert to micros if IR transition */
			if( wl < 100000.f )
			{
				units = 'A';
			}

			else
			{
				wl /= 10000.f;
				units = 'm';
			}

			if( trace.lgTrace )
			{
				fprintf( io, "# TRACE: te, htot, ctot%11.3e%11.3e%11.3e\n", 
				  phycon.te, thermal.htot, thermal.ctot );
			}

			/*fprintf( io, 
				"%10.4e%11.4e%4ld%4ld%6.3f%11.4e %4.4s %4ld%c%6.3f%10.2e%11.4e%11.4e%6.2f%6.2f%6.2f%6.2f\n",;*/
			fprintf(io, PrintEfmt("%11.4e\t", phycon.te ) );
			fprintf(io, PrintEfmt("%11.4e\t", thermal.htot ) );
			fprintf(io," [%2ld][%2ld]\t%6.3f\t",
			  ksav, jsav,  
			  hfrac);
			fprintf(io, PrintEfmt("%11.4e\t", thermal.ctot ) );
			fprintf(io," %-10s\t%.1f%c\t%6.3f\t",
			  chLabel , 
			  wl, 
			  units, 
			  cfrac );
			fprintf(io, PrintEfmt("%10.2e\t", thermal.dHeatdT ) );
			fprintf(io, PrintEfmt("%11.2e\t", thermal.dCooldT ) );
			fprintf(io, PrintEfmt("%11.4e\t", dense.eden ) );
			fprintf(io, PrintEfmt("%11.4e\t", dense.gas_phase[ipHYDROGEN] ) );
			fprintf(io, PrintEfmt("%11.4e\t", hmi.H2_total ) );
			if( dense.lgElmtOn[ipHELIUM] )
			{
				fprintf(io,"%6.2f\t%6.2f\t%6.2f\t%6.2f",
				log10(MAX2(1e-9,dense.xIonDense[ipHYDROGEN][1]/dense.gas_phase[ipHYDROGEN])), 
				log10(MAX2(1e-9,dense.xIonDense[ipHELIUM][0]/dense.gas_phase[ipHELIUM])), 
				log10(MAX2(1e-9,dense.xIonDense[ipHELIUM][1]/dense.gas_phase[ipHELIUM])), 
				log10(MAX2(1e-9,dense.xIonDense[ipHELIUM][2]/dense.gas_phase[ipHELIUM])) );
			}
			fprintf(io,"\n");
			fflush(io);
		}

		TempChange(phycon.te*factor , true);
		/* increment nmap but do not exceed nMapAlloc */
		hcmap.nmap = MIN2(hcmap.nMapAlloc,hcmap.nmap+1);

		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
				static int kount = 0;
				factor = 1.;
				TempChange(8674900. , true);
				++kount;
				if( kount >=100 )
				{
					fprintf(ioQQQ," exiting in map_do\n");
					break;
				}
			}
		}
	}

	/* now check whether sharp inflections occurred, and also find the biggest jump
	 * in the heating and cooling */
	hcmap.lgMapOK = true;
	/* >>chng 02 mar 04, lower bound had been 1, so [i-2] below was negative */
	 for( i=2; i< hcmap.nmap-2; ++i )
	 {
		realnum s1,s2,s3;/* the three slopes we will use */
		s1 = hcmap.cmap[i-2] - hcmap.cmap[i-1];
		s2 = hcmap.cmap[i-1] - hcmap.cmap[i];
		s3 = hcmap.cmap[i] - hcmap.cmap[i+1];
		if( s1*s3 > 0. && s2*s3 < 0. )
		{
			 /* of the three points, the outer had the same slope 
			  * (their product was positive) but there was an inflection
			  * between them (the negative product).  The data chain looked like
			  *     2 4
			  *    1 3  or vice versa, either case is wrong,
			  * with the logic in the above test, the problem point will aways be s2 */
			fprintf( io,
				"# cooling curve had double inflection at T=%.2e.  ",
				hcmap.temap[i]);
			fprintf( io,	" Slopes were %.2e %.2e %.2e",	s1, s2, s3);
			if( fabs(s2)/hcmap.cmap[i] > 0.05 )
			{
				fprintf( io,
					" error large, (rel slope of %.2e).\n",
					s2 / hcmap.cmap[i]);
				hcmap.lgMapOK = false;
			}
			else
			{
				fprintf( io,
					" error is small, (rel slope of %.2e).\n",
					s2 / hcmap.cmap[i]);
			}
		}

		s1 = hcmap.hmap[i-2] - hcmap.hmap[i-1];
		s2 = hcmap.hmap[i-1] - hcmap.hmap[i];
		s3 = hcmap.hmap[i]   - hcmap.hmap[i+1];
		if( s1*s3 > 0. && s2*s3 < 0. )
		{
			 /* of the three points, the outer had the same slope 
			  * (their product was positive) but there was an inflection
			  * between them (the negative product).  The data chain looked like
			  *     2 4
			  *    1 3  or vice versa, either case is wrong */
			fprintf( io,
				"# heating curve had double inflection at T=%.2e.\n",
				hcmap.temap[i] );
			hcmap.lgMapOK = false;
		}
	 }

	thermal.lgCNegChk = true;
	TempChange(t1 , false);
	return;
}
