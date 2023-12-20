/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* GammaBn evaluate photoionization rate for single shell with induced recomb */
/* GammaPrt special version of gamma function to print strong contributors */
/* GammaK evaluate photoionization rate for single shell */
/* GammaPrtRate print photo rates for all shells of a ion and element */
/*GammaPrtShells for the element nelem and ion, print total photo rate, subshells,
 * and call GamaPrt for important subshells */
#include "cddefines.h"
#include "gammas.h"
#include "thermal.h"
#include "secondaries.h"
#include "opacity.h"
#include "rfield.h"
#include "ionbal.h"
#include "heavy.h"
#include "phycon.h"
#include "vectorize.h"
#include "dense.h"

/*
 * these are the routines that evaluate the photoionization rates, gammas,
 * throughout cloudy.  a considerable amount of time is spent in these routines,
 * so they should be compiled at the highest possible efficientcy.  
 */

/* removed these two unused routines around r3500 (late October, 2009). */
/* GammaBnPL evaluate photoionization rate for single shell with induced recomb */
/* GammaPL evaluate photoionization rate for power law photo cross section */


/*>>chng 99 apr 16, ConInterOut had not been included in the threshold point for any
 * of these routines.  added it.  also moved thresholds above loop for a few */

double GammaBn(
  /* index for low energy */
  long int ipLoEnr, 
  /* index for high energy */
  long int ipHiEnr, 
  /* offset within the opacity stack */
  long int ipOpac, 
  /* threshold (Ryd) for arbitrary photoionization */
  double thresh, 
  /* induced rec rate */
  double *ainduc, 
  /* rec cooling */
  double *rcool,
  t_phoHeat *photoHeat)
{
	DEBUG_ENTRY( "GammaBn()" );

/**********************************************************************
 *
 * special version of GAMFUN for arbitrary threshold, and induc rec
 * compute ainduc, rate of induced recombinations, rcool, induc rec cool
 *
 **********************************************************************/

	/* special version of GAMFUN for arbitrary threshold, and induc rec
	 * compute ainduc, rate of inducted recombinations, rcool, induc rec cool */

	double bnfun_v = 0.;
	photoHeat->HeatNet = 0.;
	photoHeat->HeatHiEnr = 0.;
	photoHeat->HeatLowEnr = 0.;
	*ainduc = 0.;
	*rcool = 0.;

	if( ipLoEnr >= rfield.nflux || ipLoEnr >= ipHiEnr )
	{
		return( bnfun_v );
	}

	ASSERT( ipLoEnr >= 0 && ipHiEnr >= 0 );

	double GamHi = 0.;
	double g = 0.;
	double RateInducRec = 0.;
	double RateInducRecCool = 0.;

	/* do first point without otscon, which may have diffuse cont,
	 * extra minus one after ipOpac is correct since not present in i = */
	/* ipLoEnr is on fortran scale. test that the pointer is correct */
	long i = ipLoEnr-1;
	ASSERT( rfield.anumin(i) <= thresh && thresh < rfield.anumax(i) );
	/* this is the width of the cell from the threshold upwards */
	double widflx = rfield.anumax(i) - thresh;
	/* fac corrects for the fact that only part of the first frequency cell can photoionize */
	double fac = widflx/rfield.widflx(i);
	/* >>chng 99 apr 16, add ConInterOut */
	/* >>chng 01 jul 04, add *rfield.lgOutOnly logic */
	double phisig = fac*(rfield.flux[0][i] + rfield.otslin[i] + 
		rfield.ConInterOut[i]*rfield.lgOutOnly)*
		opac.OpacStack[i-ipLoEnr+ipOpac];
	double anubar = (thresh+rfield.anumax(i))/2.;

	if (i < secondaries.ipSecIon-1)
	{
		g += phisig;
		/* use average photon energy above threshold in this frequency cell */
		photoHeat->HeatNet += phisig*anubar;
	}
	else
	{
		GamHi += phisig;
		photoHeat->HeatHiEnr += phisig*anubar;
	}

	/* integral part of induced rec rate */
	/* multiply physig with the Boltzmann factor averaged over the ionizing part of the cell */
	double prod = -phisig*phycon.te_ryd*exp(-thresh/phycon.te_ryd)*expm1(-widflx/phycon.te_ryd)/widflx;
	RateInducRec += prod;
	/* induced rec cooling */
	RateInducRecCool += prod*(anubar - thresh);

	long iup = MIN2(ipHiEnr,rfield.nflux);
	iup = MIN2(iup,rfield.nPositive );
	long limit = MIN2(iup,secondaries.ipSecIon-1);

	/* these is no -1 after the ipLoEnr in the following, due to c off by one counting */
	for( i=ipLoEnr; i < limit; i++ )
	{
		/* SummedCon defined in RT_OTS_Update, 
		 * includes flux, otscon, otslin, ConInterOut, outlin */
		phisig = rfield.SummedCon[i]*opac.OpacStack[i-ipLoEnr+ipOpac];

		g += phisig;
		photoHeat->HeatNet += phisig*rfield.anu(i);

		/* phi-sigma times exp(-hnu/kT) */
		prod = phisig*rfield.ContBoltzAvg[i];
		/* induced recombination rate */
		RateInducRec += prod;
		/* induced rec cooling */
		RateInducRecCool += prod*(rfield.anu(i) - thresh);
	}

	/* now do part where secondary ionization is possible */
	/* make sure we don't count twice when ipSecIon = ipLoEnr */
	long ilo = i; // First zone not counted above
	for( i=ilo; i < iup; i++ )
	{
		phisig = rfield.SummedCon[i]*opac.OpacStack[i-ipLoEnr+ipOpac];

		GamHi += phisig;
		photoHeat->HeatHiEnr += phisig*rfield.anu(i);

		/* integral part of induced rec rate */
		prod = phisig*rfield.ContBoltzAvg[i];
		RateInducRec += prod;
		/* induced rec cooling */
		RateInducRecCool += prod*(rfield.anu(i) - thresh);
	}

	/* heating in Rydbergs, no secondary ionization */
	/* >>chng 02 mar 31, added MAX2 due to roundoff error
	 * creating slightly negative number, caught downstream as insanity */
	/* >>chng 16 feb 03, add ASSERT to prevent covering up sins... */
	ASSERT( photoHeat->HeatNet > g*thresh || fp_equal(photoHeat->HeatNet,g*thresh) );
	photoHeat->HeatNet = MAX2(0.,photoHeat->HeatNet - g*thresh);

	/* we will save this heat - the part that does not secondary ionize */
	photoHeat->HeatLowEnr = photoHeat->HeatNet;
	
	/* this is total photo rate, low and high energy parts */
	bnfun_v = g + GamHi;

	/* heating in Rydbergs, uncorrected for secondary ionization */
	photoHeat->HeatHiEnr -= GamHi*thresh;

	/* add corrected high energy heat to total */
	photoHeat->HeatNet += photoHeat->HeatHiEnr*secondaries.HeatEfficPrimary;

	/* final product is heating in erg */
	photoHeat->HeatNet *= EN1RYD;
	photoHeat->HeatHiEnr *= EN1RYD;
	photoHeat->HeatLowEnr *= EN1RYD;

	/* there is an option to turn off induced processes, values are zeroed above */
	if( rfield.lgInducProcess )  		
	{
		*rcool = RateInducRecCool*EN1RYD;
		*ainduc = RateInducRec;
	}

	ASSERT( bnfun_v >= 0. );
	ASSERT( photoHeat->HeatNet>= 0. );
	return bnfun_v;
}

/*GammaPrtShells for the element nelem and ion, print total photo rate, subshells,
 * and call GamaPrt for important subshells */
void GammaPrtShells( long nelem , long ion )
{
	DEBUG_ENTRY( "GammaPrtShells()" );

	if( !dense.lgElmtOn[nelem] )
		return;

	fprintf(ioQQQ," GammaPrtShells nz\t%.2f \t%.2li %.2li ", 
		fnzone ,
		nelem,
		ion
		);
	double sum = 0;
	for( long ns=0; ns < Heavy.nsShells[nelem][ion]; ++ns )
	{
		sum += ionbal.PhotoRate_Shell[nelem][ion][ns][0];
	}
	fprintf(ioQQQ,"\ttot\t%.2e", sum);
	/* loop over all shells for this ion */
	for( long ns=0; ns < Heavy.nsShells[nelem][ion]; ++ns )
	{
		fprintf(ioQQQ,"\t%.2e", ionbal.PhotoRate_Shell[nelem][ion][ns][0] );
#		if 0
		/* always reevaluate the outer shell, and all shells if lgRedoStatic is set */
		if( (ns==(Heavy.nsShells[nelem][ion]-1) || opac.lgRedoStatic) )
		{
			/* option to redo the rates only on occasion */
			iplow = opac.ipElement[nelem][ion][ns][0];
			iphi = opac.ipElement[nelem][ion][ns][1];
			ipop = opac.ipElement[nelem][ion][ns][2];

			/* compute the photoionization rate, ionbal.lgPhotoIoniz_On is 1, set 0
				* with "no photoionization" command */
			ionbal.PhotoRate_Shell[nelem][ion][ns][0] = 
				GammaK(iplow,iphi,
				ipop,t_yield::Inst().elec_eject_frac(nelem,ion,ns,0),
				&phoHeat)*ionbal.lgPhotoIoniz_On;
		}
#		endif
	}
	fprintf(ioQQQ,"\n");
	return;
}

/**********************************************************************
 *
 * GammaPrt special version of gamma function to print strong contributors
 * this only prints info - does not update heating rates properly since
 * this is only a diagnostic
 *
 **********************************************************************/

void GammaPrt(
	  /* first three par are identical to GammaK */
	  long int ipLoEnr, 
	  long int ipHiEnr, 
	  long int ipOpac, 
	  /* io unit we will write to */
	  FILE * ioFILE, 
	  /* total photo rate from previous call to GammaK */
	  double total, 
	  /* we will print contributors that are more than this rate */
	  double threshold)
{
	long int i, 
	  j, 
	  k;
	double flxcor, 
	  phisig;

	DEBUG_ENTRY( "GammaPrt()" );

	/* special special version of GAMFUN to save step-by-step results */
	/* >>chng 99 apr 02, test for lower greater than higher (when shell
	 * does not exist) added.  This caused incorrect photo rates for
	 * non-existant shells */
	if( ipLoEnr >= rfield.nflux || ipLoEnr >= ipHiEnr )
	{ 
		return;
	}

	fprintf( ioFILE, " GammaPrt %.2f from ",fnzone);
	fprintf( ioFILE,PrintEfmt("%9.2e",rfield.anu(ipLoEnr-1)));
	fprintf( ioFILE, " to ");
	fprintf( ioFILE,PrintEfmt("%9.2e",rfield.anu(ipHiEnr-1)));
	fprintf( ioFILE, "R rates >");
	fprintf( ioFILE,PrintEfmt("%9.2e",threshold));
	fprintf( ioFILE, " of total=");
	fprintf( ioFILE,PrintEfmt("%9.2e",total));
	fprintf( ioFILE, " (frac inc, otslin, otscon, ConInterOut, outlin ConOTS_local_OTS_rate ) chL, C\n");

	if( threshold <= 0. || total <= 0. )
	{ 
		return;
	}

	k = ipLoEnr;
	j = MIN2(ipHiEnr,rfield.nflux);

	/* do theshold special, do not pick up otscon */
	i = k-1;

	/* >>chng 01 jul 04, add *rfield.lgOutOnly logic */
	phisig = (rfield.flux[0][i] + rfield.otslin[i]+ rfield.ConInterOut[i]*rfield.lgOutOnly)*
		opac.OpacStack[i-k+ipOpac];
	if( phisig > threshold || phisig < 0.)
	{
		flxcor = rfield.flux[0][i] + rfield.otslin[i] + rfield.ConInterOut[i]*rfield.lgOutOnly;

		/* this really is array index on C scale */
		fprintf( ioFILE, "[%5ld]" , i );
		fprintf( ioFILE, PrintEfmt("%9.2e",rfield.anu(i)));
		fprintf( ioFILE, PrintEfmt("%9.2e",phisig/total));
		fprintf( ioFILE, "%5.2f%5.2f%5.2f%5.2f%5.2f%5.2f %4.4s %4.4s %.2e \n", 
		  rfield.flux[0][i]/SDIV(flxcor), 
		  rfield.otslin[i]/SDIV(flxcor), 
		  /* otscon will appear below, but is not counted here, so do not print it (deceiving) */
		  0./SDIV(flxcor), 
		  rfield.ConInterOut[i]*rfield.lgOutOnly/SDIV(flxcor), 
		  (rfield.outlin[0][i]+rfield.outlin_noplot[i])/SDIV(flxcor), 
		  rfield.ConOTS_local_OTS_rate[i]/SDIV(flxcor), 
		  rfield.chLineLabel[i].c_str(), 
		  rfield.chContLabel[i].c_str(), 
		  opac.OpacStack[i-k+ipOpac]);
	}

	for( i=k; i < j; i++ )
	{
		phisig = rfield.SummedCon[i]*opac.OpacStack[i-k+ipOpac];
		if( phisig > threshold || phisig < 0.)
		{
			/* >>chng 01 jul 04, add *rfield.lgOutOnly logic */
			flxcor = rfield.flux[0][i] + rfield.otslin[i] + rfield.otscon[i] + 
			  rfield.outlin[0][i] + rfield.outlin_noplot[i] +rfield.ConInterOut[i]*rfield.lgOutOnly;

			fprintf( ioFILE, "%5ld", i );
			fprintf(ioFILE,PrintEfmt("%9.2e",rfield.anu(i)));
			fprintf(ioFILE,PrintEfmt("%9.2e",phisig/total));
			fprintf( ioFILE, "%5.2f%5.2f%5.2f%5.2f%5.2f%5.2f %4.4s %4.4s %.2e \n", 
			  rfield.flux[0][i]/SDIV(flxcor), 
			  rfield.otslin[i]/SDIV(flxcor), 
			  rfield.otscon[i]/SDIV(flxcor), 
			  rfield.ConInterOut[i]*rfield.lgOutOnly/SDIV(flxcor), 
			  (rfield.outlin[0][i] + rfield.outlin_noplot[i])/SDIV(flxcor), 
			  rfield.ConOTS_local_OTS_rate[i]/SDIV(flxcor), 
			  rfield.chLineLabel[i].c_str(), 
			  rfield.chContLabel[i].c_str(), 
			  opac.OpacStack[i-k+ipOpac]);
		}
	}
	return;
}

/*GammaK evaluate photoionization rate for single shell */

/* this routine is a major pace setter for the code
 * carefully check anything put into the following do loop */

double GammaK(
	/* ipLoEnr and ipHiEnr are pointers within frequency array to upper and
	 * lower bounds of evaluation */
	long int ipLoEnr, 
	long int ipHiEnr, 
	/* ipOpac is offset to map onto OPSV*/
	long int ipOpac, 
	/* yield1 is fraction of ionizations that emit 1 electron only,
	 * only used for energy balance */
	double yield1,
	t_phoHeat *photoHeat)
{
	DEBUG_ENTRY( "GammaK()" );

	/* evaluate photoioinzation rate and photoelectric heating
	 * returns photoionization rate as GAMK, heating is H	 */

	/* >>chng 99 apr 02, test for lower greater than higher (when shell
	 * does not exist) added.  This caused incorrect photo rates for
	 * non-existant shells */
	if( ipLoEnr >= rfield.nflux || ipLoEnr >= ipHiEnr )
	{
		double gamk_v = 0.;
		photoHeat->HeatNet = 0.;
		photoHeat->HeatHiEnr = 0.;
		photoHeat->HeatLowEnr = 0.;
		return gamk_v;
	}

	long iup = MIN2(ipHiEnr,rfield.nPositive);

	/* anu(i) is threshold, assume each auger electron has this energy
	 * less threshold energy, IP lost in initial photoionizaiton
	 * yield1 is the fraction of photos that emit 1 electron */
	double eauger = rfield.anu(ipLoEnr-1)*yield1;

	/* low energies where secondary ionization cannot occur
	 * will do threshold point, ipLoEnr, later */
	double gamk_v = 0.;
	photoHeat->HeatNet = 0.;

	/* set up total offset for pointer addition not in loop */
	long ipOffset = ipOpac - ipLoEnr;

	/* >>>chng 99 apr 16, this had followed the loop below, moved here to 
	 * be like rest of gamma functions */
	/* first do the threshold point
	 * do not include otscon, which may contain diffuse continuum */
	/* >>chng 01 jul 04, add *rfield.lgOutOnly logic */
	double phisig = (rfield.flux[0][ipLoEnr-1] + rfield.otslin[ipLoEnr-1]+ 
		rfield.ConInterOut[ipLoEnr-1]*rfield.lgOutOnly) * opac.OpacStack[ipOpac-1];
	gamk_v += phisig;
	photoHeat->HeatNet += phisig*rfield.anu(ipLoEnr-1);

	/* this loop only executed for energies than cannot sec ioniz */
	long limit = MIN2(iup,secondaries.ipSecIon-1);
	double help;
	photoHeat->HeatNet +=
		reduce_abc_ab( get_ptr(rfield.SummedCon), &opac.OpacStack[ipOffset], rfield.anuptr(), ipLoEnr, limit, &help );
	gamk_v += help;

	/* correct heating for work function */
	/* >>chng 02 apr 10, from first to sec, due to neg heat in blister.in */
	/* *photoHeat->HeatNet += -gamk_v*eauger;*/
	ASSERT( photoHeat->HeatNet >= 0. );
	photoHeat->HeatNet = MAX2(0.,photoHeat->HeatNet - gamk_v*eauger);

	/* this is the total low-energy heating, in ryd, that cannot secondary ionize */
	photoHeat->HeatLowEnr = photoHeat->HeatNet;

	/* now do part where secondary ionization possible */
	/* make sure we don't count twice when ipSecIon = ipLoEnr */
	long ilo = MAX2(ipLoEnr,secondaries.ipSecIon-1);
	double GamHi;
	photoHeat->HeatHiEnr =
		reduce_abc_ab( get_ptr(rfield.SummedCon), &opac.OpacStack[ipOffset], rfield.anuptr(), ilo, iup, &GamHi );

	/* add on the high energy part */
	gamk_v += GamHi;
	/* correct for work function */
	photoHeat->HeatHiEnr -= GamHi*eauger;

	/* net heating include high energy heat, with correction for
	 * secondary ionization */
	photoHeat->HeatNet += photoHeat->HeatHiEnr*secondaries.HeatEfficPrimary;

	/* final product is heating in erg */
	photoHeat->HeatNet *= EN1RYD;
	photoHeat->HeatLowEnr *= EN1RYD;
	photoHeat->HeatHiEnr *= EN1RYD;

	ASSERT( gamk_v >= 0. );
	ASSERT( photoHeat->HeatNet>= 0. );
	return gamk_v;
}

/* GammaPrtRate will print resulting rates for ion and element */
void GammaPrtRate(
	/* io unit we will write to */
	FILE * ioFILE, 
	/* stage of ionization on C scale, 0 for atom */
	long int ion ,
	/* 0 for H, etc */
	long int nelem,
	/* true - then print photo sources for each shell */
	bool lgPRT )
{
	DEBUG_ENTRY( "GammaPrtRate()" );

	if( !dense.lgElmtOn[nelem] )
		return;

	/* number of shells for this species */
	long nshell = Heavy.nsShells[nelem][ion];

	/* now print subshell photo rate */
	fprintf(ioFILE , "GammaPrtRate: %li %li",ion , nelem );
	for( long ns=nshell-1; ns>=0; --ns )
	{
		fprintf(ioFILE , " %.2e" ,  ionbal.PhotoRate_Shell[nelem][ion][ns][0]);

		/* option to print individual contributors to each shell */
		if( lgPRT )
		{
			fprintf(ioFILE , "\n");
			GammaPrt(
				/* first three par are identical to GammaK */
				opac.ipElement[nelem][ion][ns][0], 
				opac.ipElement[nelem][ion][ns][1], 
				opac.ipElement[nelem][ion][ns][2], 
				/* io unit we will write to */
				ioFILE, 
				/* total photo rate from previous call to GammaK */
				ionbal.PhotoRate_Shell[nelem][ion][ns][0], 
				/* we will print contributors that are more than this rate */
				ionbal.PhotoRate_Shell[nelem][ion][ns][0]*0.05);
		}
	}
	fprintf(ioFILE , "\n");
	return;
}
