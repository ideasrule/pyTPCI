/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*optimize_func actual function called during evaluation of optimization run */
#include "cddefines.h"
#include "init.h"
#include "lines.h"
#include "called.h"
#include "predcont.h"
#include "radius.h"
#include "rfield.h"
#include "input.h"
#include "cloudy.h"
#include "cddrive.h"
#include "grid.h"
#include "flux.h"
/* used below to derive chi2 */
STATIC double chi2_func(double,double,double);

chi2_type optimize_func(const realnum param[],
			int grid_index)
{
	const int MAXCAT = 6;

	static const char name_cat[MAXCAT][13] = 
	{
	  "rel flux    ",
	  "column dens ",
	  "abs flux    ",
	  "mean temp   ",
	  "ang diameter",
	  "photometry  "
	};

	bool lgBAD,
	  lgLimOK;

	long int cat,
	  i, 
	  j, 
	  nfound, 
	  nobs_cat[MAXCAT],
	  np;

	chi2_type chi1, 
	  chi2_cat[MAXCAT],
	  chisq, 
	  help,
	  predin,
	  scld,
	  snorm,
	  theocl,
	  temp_theory;

	DEBUG_ENTRY( "optimize_func()" );

	if( grid_index >= 0 )
		optimize.nOptimiz = grid_index;

	/* This routine is called by optimizer with values of the
	 * variable parameters for CLOUDY in the array p(i). It returns
	 * the value FUNC = SUM (obs-model)**2/sig**2 for the lines
	 * specified in the observational data file, values held in the
	 * common blocks /OBSLIN/ & /OBSINT/
	 * replacement input strings for CLOUDY READR held in /chCardSav/
	 * parameter information for setting chCardSav held in /parmv/
	 * additional variables
	 * Gary's variables
	 */

	if( optimize.lgOptimFlow )
	{
		fprintf( ioQQQ, " trace, optimize_func variables" );
		for( i=0; i < optimize.nvary; i++ )
		{
			fprintf( ioQQQ, "%.2e", param[i] );
		}
		fprintf( ioQQQ, "\n" );
	}

	for( i=0; i < optimize.nvary; i++ )
	{
		optimize.vparm[0][i] = param[i];
	}

	/* call routine to pack variables with appropriate
	 * CLOUDY input lines given the array of variable parameters p(i) */
	vary_input( &lgLimOK, grid_index );

	// nothing more to be done...
	if( strcmp(optimize.chOptRtn,"XSPE") == 0 )
		return 0.;

	/* zero out lots of variables */
	zero();

	for( i=0; i < optimize.nvary; i++ )
	{
		optimize.varmax[i] = max(optimize.varmax[i],min(param[i],optimize.varang[i][1]));
		optimize.varmin[i] = min(optimize.varmin[i],max(param[i],optimize.varang[i][0]));
	}

	if( !lgLimOK )
	{
		/* these parameters are not within limits of parameter search
		 * >>chng 96 apr 26, as per Peter van Hoof comment */
		fprintf( ioQQQ, " Iteration %ld not within range.\n", 
		  optimize.nOptimiz );

		/* always increment nOptimiz, even if parameters are out of bounds,
		 * this prevents optimizer to get stuck in infinite loop */
		++optimize.nOptimiz;

		/* this is error; very bad since not within range of parameters */
		return BIG_CHI2;
	}

	lgBAD = cloudy();
	if( lgBAD )
	{
		fprintf( ioQQQ, " PROBLEM Cloudy returned error condition - what happened?\n" );
	}

	if( grid.lgGrid )
	{
		/* this is the function's return value */
		chisq = 0.;
	}
	else
	{
		/* this branch optimizing, not grid 
		/ * extract line fluxes and compare with observations */
		chisq = 0.0;
		for( i=0; i < MAXCAT; i++ )
		{
			nobs_cat[i] = 0;
			chi2_cat[i] = 0.0;
		}

		if( LineSave.ipNormWavL < 0 )
		{
			fprintf( ioQQQ, 
				" Normalization line array index is bad.  What has gone wrong?\n" );
			cdEXIT(EXIT_FAILURE);
		}

		if( (snorm = LineSave.lines[LineSave.ipNormWavL].SumLine(optimize.nEmergent)) == 0. )
		{
			fprintf( ioQQQ, "\n\n PROBLEM Normalization line has zero intensity.  What has gone wrong?\n" );
			fprintf( ioQQQ, " Is spectrum normalized to a species that does not exist?\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* first print all warnings */
		cdWarnings(ioQQQ);

		/* print header before any of the individual chi2 values are printed */
		int wlab = 4;
		for( i=0; i < long(optimize.chTempLab.size()); i++ )
		{
			wlab = max(wlab,optimize.chTempLab[i].size());
		}
		for( i=0; i < long(optimize.chColDen_label.size()); i++ )
		{
			wlab = max(wlab,optimize.chColDen_label[i].size());
		}
		if( optimize.lgOptimize )
			fprintf( ioQQQ, " %*.*s%7s%12.12s%12.12s%12.12s%12.12s Type\n",
						wlab,wlab,"ID","","Model","Observed","error","chi**2");
		else
			ASSERT( grid.lgGrid );

		/* cycle through the observational values */
		nfound = 0;

		/* first is to optimize relative emission line spectrum */
		if( optimize.xLineInt_Obs.size() > 0 )
		{
			/* set pointers to all optimized lines if first call */
			if( optimize.ipobs.size() == 0 )
			{
				optimize.ipobs.resize( optimize.xLineInt_Obs.size() );
				bool lgHIT = true;
				for( i=0; i < long(optimize.xLineInt_Obs.size()); i++ )
				{
					/* >> chng 06 may 04, use cdLine instead of ad hoc treatment.
					 * no need to complain, cdLine will do it automatically.  */
					/* this is an intensity, get the line, returns false if could not find it */
					/* >> chng 14 aug 23, use findline instead, since intensities not used */
					j = LineSave.findline( optimize.lineids[i] );
					if( j <= 0 )
					{
						fprintf( ioQQQ, "\n" );
						lgHIT = false;
					}
					else
					{
						optimize.ipobs[i] = j;
					}
				}

				/* we did not find the line */
				if( !lgHIT )
				{
					fprintf( ioQQQ, "\n\n Optimizer could not find one or more lines.\n" );
					fprintf( ioQQQ, " Sorry.\n");
					cdEXIT(EXIT_FAILURE);
				}
			}

			for( i=0; i < 10; i++ )
			{
				optimize.SavGenericData[i] = 0.;
			}

			for( i=0; i < long(optimize.xLineInt_Obs.size()); i++ )
			{
				/* and find corresponding model value by straight search */
				nfound += 1;
				scld = (chi2_type)LineSave.lines[optimize.ipobs[i]].SumLine(optimize.nEmergent)/
					(chi2_type)snorm*LineSave.ScaleNormLine;
				chi1 = chi2_func(scld,(chi2_type)optimize.xLineInt_Obs[i],
						 (chi2_type)optimize.xLineInt_error[i]);
				cat = 0;
				nobs_cat[cat]++;
				chi2_cat[cat] += chi1;

				fprintf( ioQQQ, " ");

				LineSave.lines[optimize.ipobs[i]].prt(ioQQQ);

				fprintf( ioQQQ, "%12.5f%12.5f%12.5f%12.2e Relative intensity", 
				  scld, 
				  optimize.xLineInt_Obs[i], 
				  optimize.xLineInt_error[i], 
				  chi1 );

				fprintf( ioQQQ, "\n" );

				if( i<10 )
				{
					optimize.SavGenericData[i] = chi1;
				}
			}
		}

		/* this is to optimize a mean temperature */
		for( i=0; i < long(optimize.temp_obs.size()); i++ )
		{
			if( cdTemp( optimize.chTempLab[i].c_str(), optimize.ionTemp[i],
				    &temp_theory, optimize.chTempWeight[i].c_str()) )
			{
				/* did not find column density */
				fprintf(ioQQQ," optimizer did not find column density %s %li \n",
					optimize.chTempLab[i].c_str(),optimize.ionTemp[i] );
				cdEXIT(EXIT_FAILURE);
			}
			nfound += 1;
			chi1 = chi2_func(temp_theory,(chi2_type)optimize.temp_obs[i],
					 (chi2_type)optimize.temp_error[i]);
			cat = 3;
			nobs_cat[cat]++;
			chi2_cat[cat] += chi1;

			fprintf( ioQQQ, " %*.*s%7ld%12.4e%12.4e%12.5f%12.2e Temperature\n",
						wlab, wlab,
				 optimize.chTempLab[i].c_str(), optimize.ionTemp[i], temp_theory,
				 optimize.temp_obs[i], optimize.temp_error[i], chi1 );
		}

		/* option to optimize column densities */
		for( i=0; i < long(optimize.ColDen_Obs.size()); i++ )
		{
			if( cdColm(optimize.chColDen_label[i].c_str(),optimize.ion_ColDen[i], &theocl) )
			{
				/* did not find column density */
				fprintf(ioQQQ," optimizer did not find column density %s %li \n",
					optimize.chColDen_label[i].c_str(), optimize.ion_ColDen[i] );
				cdEXIT(EXIT_FAILURE);
			}
			nfound++;
			chi1 = chi2_func(theocl,(chi2_type)optimize.ColDen_Obs[i],
					 (chi2_type)optimize.ColDen_error[i]);
			cat = 1;
			nobs_cat[cat]++;
			chi2_cat[cat] += chi1;
			
			fprintf( ioQQQ, " %*.*s%7ld%12.4e%12.4e%12.5f%12.2e Column density\n",
						wlab,wlab,
				 optimize.chColDen_label[i].c_str(), optimize.ion_ColDen[i], theocl, 
				 optimize.ColDen_Obs[i], optimize.ColDen_error[i], chi1 );
		}

		/* option to optimize line flux */
		if( optimize.lgOptLum )
		{
			++nfound;
			if( LineSave.lines[LineSave.ipNormWavL].SumLine(optimize.nOptLum) > 0.f )
			{
				predin = log10(LineSave.lines[LineSave.ipNormWavL].SumLine(optimize.nOptLum) *
						radius.Conv2PrtInten);
				help = exp10(predin-(chi2_type)optimize.optint);
				chi1 = chi2_func(help,1.,(chi2_type)optimize.optier);
			}
			else
			{
				predin = -999.99999;
				chi1 = BIG_CHI2;
			}
			cat = 2;
			nobs_cat[cat]++;
			chi2_cat[cat] += chi1;

			fprintf( ioQQQ, " ");
			LineSave.lines[LineSave.ipNormWavL].prt(ioQQQ);

			fprintf( ioQQQ, "%12.5f%12.5f%12.5f%12.2e Line intensity\n", 
			  predin,
			  optimize.optint,
			  optimize.optier,
			  chi1 );
		}

		/* option to optimize the absolute continuum flux */
		for( size_t k=0; k < optimize.ContIndex.size(); k++ )
		{
			nfound++;
			// there are 4 entries for each wavelength: nFnu, nInu, InwT, InwC
			long ind = t_PredCont::Inst().offset() + 4*optimize.ContIndex[k];
			chi2_type nFnu_model = 0.;
			if( LineSave.lines[ind].SumLine(0) > SMALLFLOAT )
			{
				nFnu_model = chi2_type( LineSave.lines[ind].SumLine(0) * radius.Conv2PrtInten );
			}
			Flux F_model( optimize.ContNFnu[k].E(), nFnu_model );

			chi1 = chi2_func(nFnu_model,optimize.ContNFnu[k].get("erg/s/cm2"),optimize.ContNFnuErr[k]);
			const char* catstr;
			// treat radio continuum flux as absolute flux so that it can be used
			// as a more accurate replacement of the normalization line intensity
			if( optimize.ContEner[k].mm() <= 1. )
			{
				cat = 5;
				catstr = "Photometry";
			}
			else
			{
				cat = 2;
				catstr = "Radio intensity";
			}
			nobs_cat[cat]++;
			chi2_cat[cat] += chi1;

			fprintf( ioQQQ, " ");
			LineSave.lines[ind].prt(ioQQQ);
			string unit = optimize.ContNFnu[k].uu();
			fprintf( ioQQQ, "%12.4g%12.4g%12.5f%12.2e %s [%s]\n", 
				 F_model.get(unit),
				 optimize.ContNFnu[k].get(unit),
				 optimize.ContNFnuErr[k], chi1,
				 catstr, unit.c_str() );
		}

		/* option to optimize angular diamater */
		if( optimize.lgOptDiam )
		{
			nfound++;
			chi2_type diam_model;
			// get diameter in cm
			if( rfield.lgUSphON )
				diam_model = 2.*rfield.rstrom; // ionization bounded -> use Stroemgren radius
			else
				diam_model = 2.*radius.Radius; // density bounded -> use outer radius
			// now convert to arcsec if necessary
			if( !optimize.lgDiamInCM && radius.distance > 0. )
				diam_model *= AS1RAD/radius.distance;

			chi1 = chi2_func(diam_model,optimize.optDiam,optimize.optDiamErr);
			cat = 4;
			nobs_cat[cat]++;
			chi2_cat[cat] += chi1;

			fprintf( ioQQQ, "            %12.4g%12.4g%12.5f%12.2e Angular diameter\n", 
				 diam_model, optimize.optDiam, optimize.optDiamErr, chi1 );
		}

		/* do not have to have line matches if doing grid. */
		if( nfound <= 0 && !grid.lgGrid )
		{
			fprintf( ioQQQ, " WARNING; no line matches found\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* write out chisquared for this iteration */
		fprintf( ioQQQ, "\n" );
		for( i=0; i < MAXCAT; i++ )
		{
			if( nobs_cat[i] > 0 )
			{
				chisq += chi2_cat[i]/nobs_cat[i];
				fprintf( ioQQQ, " Category %s #obs.%3ld  Total Chi**2%11.3e  Average Chi**2%11.3e\n",
				  name_cat[i],nobs_cat[i],chi2_cat[i],chi2_cat[i]/nobs_cat[i] );
			}
		}
		if( nfound )
		{
			fprintf( ioQQQ, "\n Iteration%4ld Chisq=%13.5e\n", optimize.nOptimiz, chisq );
		}
	}

	/* increment nOptimiz, the grid / optimizer counter */
	++optimize.nOptimiz;

	/* only print this if output has been turned on */
	if( called.lgTalk )
	{
		fprintf( ioQQQ, "\n" );
		for( i=0; i < optimize.nvary; i++ )
		{
			np = optimize.nvfpnt[i];

			/* now generate the actual command with parameters */
			input.crd[np]->chCardSav = MakeInputLine(i);

			fprintf( ioQQQ, " Varied command: %s\n", 
				 input.crd[np]->chCardSav.c_str() );
		}
	}

	return min(chisq,BIG_CHI2);
}

/* ============================================================================== */
STATIC chi2_type chi2_func(chi2_type ymodl,
			   chi2_type ymeas,
			   chi2_type yerr)
{
	chi2_type chi2_func_v,
		temp;

	DEBUG_ENTRY( "chi2_func()" );

	/* compute chi**2 by comparing model quantity ymodl with a measured
	 * quantity ymeas with relative error yerr (negative means upper limit)
	 */

	if( ymeas <= 0. )
	{
		fprintf( ioQQQ, "chi2_func: non-positive observed quantity, this should not happen\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( yerr > 0. )
	{
		if( ymodl > 0. )
		{
			temp = pow2((ymodl-ymeas)/(min(ymodl,ymeas)*yerr));
			chi2_func_v = min(temp,BIG_CHI2);
		}
		else
			chi2_func_v = BIG_CHI2;
	}
	else if( yerr < 0. )
	{
		/* value quoted is an upper limit, so add to chisq
		 * only if limit exceeded, otherwise return zero.
		 */
		if( ymodl > ymeas )
		{
			temp = pow2((ymodl-ymeas)/(ymeas*yerr));
			chi2_func_v = min(temp,BIG_CHI2);
		}
		else
			chi2_func_v = 0.;
	}
	else
	{
		fprintf( ioQQQ, "chi2_func: relative error is zero, this should not happen\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return chi2_func_v;
}
