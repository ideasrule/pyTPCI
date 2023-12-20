/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveLineData punches selected line data for all lines transferred in code */
/*Save1LineData save data for one line */
#include "cddefines.h"
#include "taulines.h"
#include "iso.h"
#include "phycon.h"
#include "elementnames.h"
#include "dense.h"
#include "conv.h"
#include "lines.h"
#include "prt.h"
#include "h2.h"
#include "thermal.h"
#include "cooling.h"
#include "save.h"
#include "mole.h"

NORETURN void SaveLineData(FILE * ioPUN)
{
	const long nskip=2; /* number of emission lines per line of output */
	double tot;
	bool lgElemOff=false;

	DEBUG_ENTRY( "SaveLineData()" );

	/* routine punches out (on unit ioPUN) line data
	 * for all recombination lines, and all transitions that are transferred */

	/* say what is happening so we know why we stopped */
	fprintf( ioQQQ, " saving line data, then stopping\n" );

	/* first check that all lines are turned on */
	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		if( !dense.lgElmtOn[nelem] )
		{
			fprintf(ioQQQ," WARNING - I am saving line data but element %s is turned off.\n",
			elementnames.chElementName[nelem]);
			lgElemOff = true;
		}
	}
	if( lgElemOff )
	{
		fprintf(ioQQQ,"Some elements are turned off and save line data requested.\n");
		fprintf(ioQQQ,"Code is now designed to do save line data only with all elements on.\n");
		fprintf(ioQQQ,"Please try again with all elements on.\n");
		fprintf(ioQQQ,"Please try again with all elements on.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* evaluate rec coefficient at constant temperature if this is set, else
	 * use 10,000K */
	double TeNew;
	if( thermal.lgTemperatureConstant )
	{
		TeNew = thermal.ConstTemp;
	}
	else
	{
		TeNew = 1e4;
	}
	TempChange(TeNew , false);

	/* this is set of Dima's recombination lines */
	t_ADfA::Inst().rec_lines(phycon.te,LineSave.RecCoefCNO);
	fprintf( ioPUN, "\n       Recombination lines of C, N, O\n" );
	fprintf( ioPUN, "#Ion\tWL(A)\tCoef\tIon\tWL(A)\tCoef\n" );
	for( long i=0; i<NRECCOEFCNO; i+=nskip)
	{
		/* nskip is set to 2 above */
		long limit = MIN2(NRECCOEFCNO,i+nskip);
		fprintf( ioPUN, "    " );
		for( long j=i; j < limit; j++ )
		{
			string chLabel = chIonLbl( LineSave.RecCoefCNO[0][i], long(LineSave.RecCoefCNO[0][i]-LineSave.RecCoefCNO[1][i]+1.01) );

			fprintf( ioPUN, "%s\t%6ld\t%8.3f\t",
					chLabel.c_str(),
					(long)(LineSave.RecCoefCNO[2][j]+0.5), 
					log10(SDIV(LineSave.RecCoefCNO[3][j]) ) );
		}
		fprintf( ioPUN, "    \n" );
	}
	fprintf( ioPUN, "\n\n" );

	dense.SetGasPhaseDensity( ipHYDROGEN, 1. );
	dense.EdenHCorr = 1.;
	EdenChange( 1. );

	/* want very small neutral fractions so get mostly e- cs */
	dense.xIonDense[ipHYDROGEN][1] = 1.e-5f;
	findspecieslocal("H2")->den = 0.;
	dense.xIonDense[ipHYDROGEN][1] = 1.;

	for( long i=0; i < nWindLine; i++ )
	{
		TauLine2[i].Lo()->Pop() = 1.;
	}

	for( size_t i=0; i < UTALines.size(); i++ )
	{
		(*UTALines[i].Lo()).Pop() = 1.;
	}

	for( long i=0; i < LIMELM; i++ )
	{
		for( long j=0; j < LIMELM+1; j++ )
		{
			dense.xIonDense[i][j] = 1.;
		}
	}

	/* evaluate cooling, this forces evaluation of collision strengths */
	CoolEvaluate(&tot);

	fprintf( ioPUN, "       Level 2 transferred lines\n" );
	PrintLineDataHeader( ioPUN );
	for( long i=0; i < nWindLine; i++ )
	{
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
		{
			Save1LineData( TauLine2[i] , ioPUN , true );
		}
	}

	fprintf( ioPUN, "\n\n\n  end level 2, start inner shell UTA\n" );
	PrintLineDataHeader( ioPUN );
	for( size_t i=0; i < UTALines.size(); i++ )
	{
		Save1LineData( UTALines[i] , ioPUN , true );
	}

	fprintf( ioPUN, "\n\n\n  end inner shell, start H-like iso seq\n" );
	/* h-like iso sequence */
	/* the hydrogen like iso-sequence */
	PrintLineDataHeader( ioPUN );
	for( long nelem=0; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			iso_collide( ipH_LIKE, nelem );
			/* arrays are dim'd iso_sp[ipH_LIKE][nelem].numLevels_max+1 */
			/* keep this limit to iso.numLevels_max, instead of _local.  */
			for( long ipLo=ipH1s; ipLo < iso_sp[ipH_LIKE][nelem].numLevels_max-1; ipLo++ )
			{
				for( long ipHi=ipLo+1; ipHi < iso_sp[ipH_LIKE][nelem].numLevels_max; ipHi++ )
				{
					Save1LineData( iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo) , ioPUN , false );
				}
			}
		}
	}

	fprintf( ioPUN, "\n\n\n  end H-like iso seq, start He-like iso seq\n" );
	PrintLineDataHeader( ioPUN );
	for( long nelem=1; nelem < LIMELM; nelem++ )
	{
		if( nelem < 2 || dense.lgElmtOn[nelem] )
		{
			/* arrays are dim'd iso_sp[ipH_LIKE][nelem].numLevels_max+1  */
			for( long ipLo=ipHe1s1S; ipLo < iso_sp[ipHE_LIKE][nelem].numLevels_max-1; ipLo++ )
			{
				for( long ipHi=ipLo+1; ipHi < iso_sp[ipHE_LIKE][nelem].numLevels_max; ipHi++ )
				{
					Save1LineData( iso_sp[ipHE_LIKE][nelem].trans(ipHi,ipLo) , ioPUN , false );
				}
			}
		}
	}

	fprintf( ioPUN, "\n\n\n   end he-like iso seq, start hyperfine structure lines\n" );
	/* fine structure lines */
	PrintLineDataHeader( ioPUN );
	for( size_t i=0; i < HFLines.size(); i++ )
	{
		Save1LineData( HFLines[i] , ioPUN , true );
	}

	fprintf( ioPUN, "\n\n\n end hyperfine, start database lines\n" );
	/* Databases: Atoms & Molecules*/
	PrintLineDataHeader( ioPUN );
	for (int ipSpecies=0; ipSpecies < nSpecies; ++ipSpecies)
	{
		for( EmissionList::iterator em=dBaseTrans[ipSpecies].Emis().begin();
			  em != dBaseTrans[ipSpecies].Emis().end(); ++em)
		{
			Save1LineData( (*em).Tran() , ioPUN , true );
		}
	}

	fprintf( ioPUN, "\n\n\n  end database, start satellite lines\n" );
	PrintLineDataHeader( ioPUN );
	for( long ipISO = ipHE_LIKE; ipISO < NISO; ipISO++ )
	{
		for( long nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] && iso_ctrl.lgDielRecom[ipISO] ) 
			{
				for( long i=0; i<iso_sp[ipISO][nelem].numLevels_max; i++ )
				{
					Save1LineData( SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][i]],
							ioPUN , true );
				}
			}
		}
	}

	/* want very small ionized fractions so get mostly H2 cs */
	dense.SetGasPhaseDensity( ipHYDROGEN, 1e-6f );
	dense.EdenHCorr = 1e-6f;
	findspecieslocal("H2")->den = 1.;
	findspecieslocal("H2*")->den = 1.;
	dense.xIonDense[ipHYDROGEN][1] = 1e-6;
	EdenChange( 1e-6 );

	/* H2 molecule */
	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end satellite, start H2 lines\n" );

	/* ioPUN unit, and option to print all possible lines - false indicates
	 * save only significant lines */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
	{
		bool lgPopsConverged;
		double old_val, new_val;
		(*diatom)->H2_LevelPops( lgPopsConverged, old_val, new_val );
		(*diatom)->H2_Punch_line_data( ioPUN, false );
	}

	fprintf( ioPUN, "\n\n\n" );
	fprintf( ioPUN, "       end H2\n" );

	/* ChkMonitorstring is searched for by one of the scripts in the nightly run
	 * this run will terminate with no asserts but that is the correct behavior */
	fprintf( ioQQQ , "\n The code is left in a disturbed state after creating the SAVE LINE DATA file.\n"
			" No calculation is actually performed, only the SAVE LINE DATA file is produced.\n"
			" Remove the SAVE LINE DATA command to do the calculation.\n\n ChkMonitorend is ok.\n" );

	/* stop when done, we have done some serious damage to the code */
	cdEXIT(EXIT_SUCCESS);
}

/* Print header for Save1LineData() function */
void PrintLineDataHeader( FILE * ioPUN )
{
	fprintf( ioPUN, "#Ion\tWL\tgl\tgu\tgf\tA\tCS\tn(crt)\tdamp\n" );
}


/*Save1LineData save data for one line */
void Save1LineData(
	const TransitionProxy& t,
	FILE * ioPUN, 
	/* flag saying whether to give collision strength too - in multi level atoms
	 * it will be not valid without a great deal more work */
	bool lgCS_2 )
{
	double CritDen;

	DEBUG_ENTRY( "Save1LineData()" );

	if( t.ipCont() <= 0 )
	{
		// not a real line, just give \n
		return;
	}

	/** \todo	1	define lifetime and collision rate for multi-level species so that the
	 * critical density is derived correctly in this routine.  For now the flag lgCS_2
	 * being true means to save critical den and is only true for two-level systems 
	 * all places where this routine is called with lgCS_2 false need to be fixed */

	/*iWL = iWavLen( t , &chUnits , &chShift );*/
	/* ion label, like C  1 */
	fprintf(ioPUN,"%s\t", chIonLbl( t ).c_str() );

	/* this is the second piece of the line label, pick up string after start */

	/* the wavelength */
	if( strcmp( save.chConSavEnr[save.ipConPun], "labl" )== 0 )
	{
		prt_wl( ioPUN , t.WLAng() );
	}
	else
	{
		/* this converts energy in Rydbergs into any of the other units */
		fprintf( ioPUN , "%.5e", AnuUnit((realnum)(t.EnergyRyd())) );
	}

	fprintf( ioPUN, "\t%3ld\t%3ld",
	  /* lower and upper stat weights */
				(long)((*t.Lo()).g()), 
				(long)((*t.Hi()).g()) );

	/* oscillator strength */
	fprintf( ioPUN,PrintEfmt("\t%9.2e",  t.Emis().gf()));

	/* Einstein A for transition */
	fprintf( ioPUN,PrintEfmt("\t%9.2e",  t.Emis().Aul()));

	/* next collision strengths, use different formats depending on size 
	 * of collision strength */
	if( t.Coll().col_str() > 100. )
	{
		fprintf( ioPUN, "\t%7.1f", t.Coll().col_str() );
	}
	else if( t.Coll().col_str() > 10. )
	{
		fprintf( ioPUN, "\t%7.2f", t.Coll().col_str() );
	}
	else if( t.Coll().col_str() > 1. )
	{
		fprintf( ioPUN, "\t%7.3f", t.Coll().col_str() );
	}
	else if( t.Coll().col_str() > .01 )
	{
		fprintf( ioPUN, "\t%7.4f", t.Coll().col_str() );
	}
	else if( t.Coll().col_str() > 0.0 )
	{
		fprintf( ioPUN, "\t%.3e", t.Coll().col_str() );
	}
	else
	{
		fprintf( ioPUN, "\t%7.4f", 0. );
	}

	/* now print critical density but only if cs is positive 
	 * >>chng 06 mar 24, add flag lgCS_2 - in multi-level systems do not want
	 * to save cs since not computed properly */
	if( lgCS_2 && t.Coll().col_str()> 0. )
	{
		CritDen = t.Emis().Aul() * (*t.Hi()).g()*phycon.sqrte / (t.Coll().col_str()*COLL_CONST);
	}
	else
	{
		CritDen = 0.;
	}
	fprintf( ioPUN, "\t%.3e",CritDen );

	// damping constant for current conditions
	fprintf( ioPUN,PrintEfmt("\t%9.2e",  t.Emis().damp() ));

	fprintf( ioPUN, "\n" );

	return;
}
