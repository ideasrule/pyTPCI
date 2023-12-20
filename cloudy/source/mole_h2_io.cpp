/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*H2_ParseSave parse the save h2 command */
/*H2_PunchDo save some properties of the large H2 molecule */
/*chMolBranch returns a char with the spectroscopic branch of a transition */
/*H2_PunchLineStuff include H2 lines in punched optical depths, etc, called from SaveLineStuff */
/*H2_Punch_line_data save line data for H2 molecule */
/*H2_Read_hminus_distribution read distribution function for H2 population following formation from H minus */
/*H2_ReadDissprob read dissociation probabilities and kinetic energies for all electronic levels */
/*H2_ReadEnergies read energies for all electronic levels */
/*H2_ReadTransprob read transition probabilities */
/*H2_Prt_Zone print H2 info into zone results, called from prtzone for each printed zone */
/*H2_ParseSave parse the save h2 command */
/*H2_Prt_column_density print H2 info into zone results, called from prtzone for each printed zone */
/*H2_LinesAdd add in explicit lines from the large H2 molecule, called by lines_molecules */
 /*cdH2_Line returns 1 if we found the line, 
  * or false==0 if we did not find the line because ohoto-para transition
  * or upper level has lower energy than lower level */
#include "cddefines.h" 
#include "cddrive.h" 
#include "save.h" 
#include "hmi.h"
#include "prt.h"
#include "secondaries.h"
#include "grainvar.h"
#include "input.h"
#include "phycon.h"
#include "rfield.h"
#include "hyperfine.h"
#include "thermal.h"
#include "lines.h"
#include "dense.h"
#include "radius.h"
#include "colden.h"
#include "h2.h"
#include "doppvel.h"
#include "parser.h"
#include "mole.h"
#include "ran.h"

static realnum thresh_punline_h2;

/*H2_LinesAdd add in explicit lines from the large H2 molecule, called by lines_molecules */
void diatomics::H2_LinesAdd(void)
{
	/* H2 not on, so space not allocated */
	if( !lgEnabled )
		return;

	DEBUG_ENTRY( "H2_LinesAdd()" );

	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		qList::iterator Hi = ( (*tr).Hi() );	
		if( (*Hi).n() >= nElecLevelOutput ) continue;
		qList::iterator Lo = ( (*tr).Lo() );	
		/* all ground vib state rotation lines - first is J to J-2 */
		PutLine( *tr, "diatoms lines", label.c_str() );
		if( LineSave.ipass == 0 )
		{
			H2_SaveLine[(*Hi).n()][(*Hi).v()][(*Hi).J()][(*Lo).n()][(*Lo).v()][(*Lo).J()] =  0.;
		}
		else if( LineSave.ipass == 1 )
		{
			H2_SaveLine[(*Hi).n()][(*Hi).v()][(*Hi).J()][(*Lo).n()][(*Lo).v()][(*Lo).J()] +=
				(realnum)( radius.dVeffAper * (*tr).Emis().xObsIntensity() );
		}
	}

	return;
}

/*H2_ParseSave parse the save h2 command */
void diatomics::H2_ParseSave( Parser &p,
			      ostringstream &chHeader)
{
	DEBUG_ENTRY( "H2_ParseSave()" );

	save.whichDiatomToPrint[save.nsave] = &(*this);

	/* this provides info on the large H2 molecule */
	if( p.nMatch("COLU") )
	{
		/* save column density */
		strcpy( save.chSave[save.nsave], "H2cl" );

		/* this is an option to scan off highest vib and rot states 
		 * to save pops - first is limit to vibration, then rotation 
		 * if no number is entered then 0 is set and all levels punched */
		/* now get vib limit */
		save.punarg[save.nsave][0] = (realnum)p.getNumberDefault(
			"H2 vibration state",0.0);

		/* highest rotation */
		save.punarg[save.nsave][1] = (realnum)p.getNumberDefault(
			"H2 rotation state",0.0);
		/* this says whether to save triplets or a matrix for output -
		 * default is triplets, so only check for matrix */
		if( p.nMatch( "MATR"  ) )
		{
			/* matrix */
			save.punarg[save.nsave][2] = 1;
			sncatf( chHeader, "#vib\trot\tcolumn density\n" );
		}
		else
		{
			/* triplets */
			save.punarg[save.nsave][2] = -1;
			sncatf( chHeader,
				"#vib\trot\tEner(K)\tcolden\tcolden/stat wght\tLTE colden\tLTE colden/stat wght\n" );
		}
	}

	else if( p.nMatch("COOL") && p.nMatch("MOLE") )
	{
		/* net cooling rate per particle */
		strcpy( save.chSave[save.nsave], "H2cm" );
		sncatf( chHeader,
			"#Temp\tLTE cooling per molecule\tNet cooling per molecule\n" );
	}

	else if( p.nMatch("COOL") )
	{
		/* heating and cooling rates */
		strcpy( save.chSave[save.nsave], "H2co" );
		sncatf( chHeader,
			"#H2 depth\tTemp\ttot cool\tTH Sol\tBig Sol\tTH pht dis\tpht dis\tTH Xcool\tXcool\tNet cool per H2\n" );
	}

	else if( p.nMatch("CREA") )
	{
		/* H2 creation rates */
		fprintf( ioQQQ, " This command has been superseded by the \"creation\" option of the \"save chemistry rates\" command.\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	else if( p.nMatch("DEST") )
	{
		/* save H2 destruction - output destruction rates */
		fprintf( ioQQQ, " This command has been superseded by the \"destruction\" option of the \"save chemistry rates\" command.\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	else if( p.nMatch("HEAT") )
	{
		/* heating and cooling rates */
		strcpy( save.chSave[save.nsave], "H2he" );
		sncatf( chHeader,
			"#H2 depth\ttot Heat\tHeat(big)\tHeat(TH85)\tDissoc(Big)\tDissoc(TH85) \n" );
	}

	else if( p.nMatch("LEVE") )
	{
		/* save H2 level energies */
		strcpy( save.chSave[save.nsave], "H2le" );
		sncatf( chHeader,
			"#H2 v\tJ\tenergy(wn)\tstat wght\tSum As" );
		for( int nColl=0; nColl<N_X_COLLIDER; ++nColl )
		{
			/* labels for all colliders */
			sncatf(chHeader,
			       "\tq(u,l,%s)",chH2ColliderLabels[nColl]);
		}
		sncatf( chHeader, "\n" );
	}

	else if( p.nMatch("LINE") )
	{
		/* save H2 lines - all in X */
		strcpy( save.chSave[save.nsave], "H2ln" );
		sncatf( chHeader,
			"#H2 line\tEhi\tVhi\tJhi\tElo\tVlo\tJlo\twl(mic)\twl(lab)\tlog L or I\tI/Inorm\tExcit(hi, K)\tg_u h nu * Aul\n" );
		/* first optional number changes the threshold of weakest line to print*/
		/* fe2thresh is intensity relative to normalization line,
		 * normally Hbeta, and is set to zero in zero.c */

		/* threshold for faintest line to save, default is 1e-4 of norm line */
		thresh_punline_h2 = (realnum)p.getNumberDefaultNegImplLog(
			"faintest line to save",1e-4);

		/* lines from how many electronic states?  default is one, just X, and is
		 * obtained with GROUND keyword.  ALL will produce all lines from all levels.
		 * else, if a number is present, will be the number.  if no number, no keyword,
		 * appear then just ground */
		if( p.nMatch( "ELEC"  ) ) 
		{
			if( p.nMatch(" ALL") )
			{
				/* all electronic levels - when done, will set upper limit, the
				 * number of electronic levels actually computed, don't know this yet,
				 * so signify with negative number */
				nElecLevelOutput = -1;
			}
			else if( p.nMatch("GROU") )
			{
				/* just the ground electronic state */
				nElecLevelOutput = 1;
			}
			else
			{
				nElecLevelOutput = (int)p.getNumberDefault(
					"electronic levels for output",1.0);
			}
		}
	}

	else if( p.nMatch(" PDR") )
	{
		/* creation and destruction processes */
		strcpy( save.chSave[save.nsave], "H2pd" );
		sncatf( chHeader, "#depth\tn(o-H2)\tn(p-H2)\tSolomon rate TH85\tSolomon rate BD96\tSolomon rate H2 model\n" );
	}
	else if( p.nMatch("POPU") )
	{
		/* save populations */
		strcpy( save.chSave[save.nsave], "H2po" );

		/* this is an option to scan off highest vib and rot states 
		 * to save pops - first is limit to vibration, then rotation 
		 * if no number is entered then 0 is set and all levels punched */
		/* now get vib lim */
		save.punarg[save.nsave][0] = (realnum)p.getNumberDefault(
			"highest H2 save vibration state",0.0);

		/* this is limit to rotation quantum index */
		save.punarg[save.nsave][1] = (realnum)p.getNumberDefault(
			"highest H2 save rotation state",0.0);

		if( p.nMatch( "ZONE"  ) )
		{
			/* save v=0 pops for each zone, all along one line */
			save.punarg[save.nsave][2] = 0;
			sncatf( chHeader,
				 "#depth\torth\tpar\te=1 rel pop\te=2 rel pop\tv,J rel pops\n" );
		}
		else
		{
			/* will not do zone output, only output at the end of the calculation
			 * now check whether to save triplets or a matrix for output -
			 * default is triplets, so only check for matrix */
			if( p.nMatch( "MATR"  ) )
			{
				/* matrix */
				save.punarg[save.nsave][2] = 1;
				sncatf( chHeader, "#vib\trot\tpops\n" );
			}
			else
			{
				/* triplets */
				save.punarg[save.nsave][2] = -1;
				sncatf( chHeader,
					"#vib\trot\ts\tenergy(wn)\tpops/H2\told/H2\tpops/g/H2\tdep coef\tFin(Col)\tFout(col)\tRCout\tRRout\tRCin\tRRin\n" );
			}
		}
	}

	else if( p.nMatch("RATE") )
	{
		/* save h2 rates - creation and destruction rates */
		strcpy( save.chSave[save.nsave], "H2ra" );
		sncatf( chHeader,
			"#depth\tN(H2)\tN(H2)/u(H2)\tA_V(star)\tn(Eval)"
			"\tH2/Htot\trenorm\tfrm grn\tfrmH-\tdstTH85\tBD96\tELWERT\tBigH2\telec->H2g\telec->H2s"
			"\tG(TH85)\tG(DB96)\tCR\tEleclife\tShield(BD96)\tShield(H2)\tBigh2/G0(spc)\ttot dest"
			"\tHeatH2Dish_TH85\tHeatH2Dexc_TH85\tHeatDish_BigH2\tHeatDexc_BigH2\thtot\n" );
	}
	else if( p.nMatch("SOLO") )
	{
		/* rate of Solomon process then fracs of exits from each v, J level */
		strcpy( save.chSave[save.nsave], "H2so" );
		sncatf( chHeader,
			"#depth\tSol tot\tpump/dissoc\tpump/dissoc BigH2\tavH2g\tavH2s\tH2g chem/big H2\tH2s chem/big H2\tfrac H2g BigH2\tfrac H2s BigH2\teHi\tvHi\tJHi\tvLo\tJLo\tfrac\twl(A)\n" );
	}
	else if( p.nMatch("SPEC") )
	{
		/* special save command*/
		strcpy( save.chSave[save.nsave], "H2sp" );
		sncatf( chHeader,
			"#depth\tspecial\n" );
	}
	else if( p.nMatch("TEMP") )
	{
		/* various temperatures for neutral/molecular gas */
		strcpy( save.chSave[save.nsave], "H2te" );
		sncatf( chHeader,
			"#depth\tH2/H\tn(1/0)\tn(ortho/para)\tT(1/0)\tT(2/0)\tT(3/0)\tT(3/1)\tT(4/0)\tT(kin)\tT(21cm)\tT_sum(1/0)\tT_sum(2/0)\tT_sum(3/0)\tT_sum(3/1)\tT_sum(4/0) \n");
	}
	else if( p.nMatch("THER") )
	{
		/* thermal heating cooling processes involving H2 */
		strcpy( save.chSave[save.nsave], "H2th" );
		sncatf( chHeader,
			"#depth\tH2/H\tn(1/0)\tn(ortho/para)\tT(1/0)\tT(2/0)\tT(3/0)\tT(3/1)\tT(4/0)\tT(kin)\tT(21cm)\tT_sum(1/0)\tT_sum(2/0)\tT_sum(3/0)\tT_sum(3/1)\tT_sum(4/0) \n");
	}
	else
	{
		fprintf( ioQQQ, 
			" There must be a second key; they are  RATE, LINE, COOL, COLUMN, _PDR, SOLOmon, TEMP, and POPUlations\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}


/*H2_Prt_Zone print H2 info into zone results, called from prtzone for each printed zone */
void diatomics::H2_Prt_Zone(void)
{
	/* no print if H2 not turned on, or not computed for these conditions */
	if( !lgEnabled || !nCall_this_zone )
		return;

	DEBUG_ENTRY( "H2_Prt_Zone()" );

	fprintf( ioQQQ, " %s density   ", label.c_str() );
	fprintf(ioQQQ,PrintEfmt("%9.2e", (*dense_total)));

	fprintf( ioQQQ, " orth/par");
	fprintf(ioQQQ,PrintEfmt("%9.2e", ortho_density / SDIV( para_density )));

	fprintf( ioQQQ, " v0 J=0,3");
	fprintf(ioQQQ,PrintEfmt("%9.2e", states[ ipEnergySort[0][0][0] ].Pop() / (*dense_total)));
	fprintf(ioQQQ,PrintEfmt("%9.2e", states[ ipEnergySort[0][0][1] ].Pop() / (*dense_total)));
	fprintf(ioQQQ,PrintEfmt("%9.2e", states[ ipEnergySort[0][0][2] ].Pop() / (*dense_total)));
	fprintf(ioQQQ,PrintEfmt("%9.2e", states[ ipEnergySort[0][0][3] ].Pop() / (*dense_total)));

	fprintf( ioQQQ, " TOTv=0,3");
	fprintf(ioQQQ,PrintEfmt("%9.2e", pops_per_vib[0][0] / (*dense_total)));
	fprintf(ioQQQ,PrintEfmt("%9.2e", pops_per_vib[0][1] / (*dense_total)));
	fprintf(ioQQQ,PrintEfmt("%9.2e", pops_per_vib[0][2] / (*dense_total)));
	fprintf(ioQQQ,PrintEfmt("%9.2e", pops_per_vib[0][3] / (*dense_total)));
	fprintf( ioQQQ, "\n");
	return;
}

void diatomics::H2_PrtDepartCoef(void)
{
	/* no print if H2 not turned on, or not computed for these conditions */
	if( !lgEnabled || !nCall_this_zone )
		return;

	DEBUG_ENTRY( "H2_PrtDepartCoef()" );

	// print departure coefficients
	fprintf( ioQQQ, " %s departure coefficients\n", label.c_str() );
	for( long iElec=0; iElec<n_elec_states; ++iElec )
	{
 		fprintf( ioQQQ, "%li electronic\n", iElec );
		for( long iVib=0; iVib<=nVib_hi[iElec]; ++iVib )
		{
			for( long iRot=0; iRot<Jlowest[iElec]; ++iRot )
				fprintf( ioQQQ, " -----" );
			for( long iRot=Jlowest[iElec]; iRot<=nRot_hi[iElec][iVib]; ++iRot )
			{
				long i = ipEnergySort[iElec][iVib][iRot];
				fprintf( ioQQQ, " %5.3f", depart[i] );
			}
 			fprintf( ioQQQ, "\n" );
		}
 		fprintf( ioQQQ, "\n" );
		if( iElec==0 )
			break;
	}

	return;
}

/*H2_Prt_column_density print H2 info into zone results, called from prtzone for each printed zone */
void diatomics::H2_Prt_column_density(	
	/* this is stream used for io, is stdout when called by final,
	 * is save unit when save output generated */
	 FILE *ioMEAN )

{
	int iVibHi;

	/* no print if H2 not turned on, or not computed for these conditions */
	if( !lgEnabled || !nCall_this_zone )
		return;

	DEBUG_ENTRY( "H2_Prt_column_density()" );

	fprintf( ioMEAN, " H2 total   ");
	fprintf(ioMEAN,"%7.3f", log10(SDIV(ortho_colden + para_colden)));

	fprintf( ioMEAN, " H2 ortho   ");
	fprintf(ioMEAN,"%7.3f", log10(SDIV(ortho_colden)));

	fprintf( ioMEAN, " para");
	fprintf(ioMEAN,"%7.3f", log10(SDIV(para_colden)));

	iVibHi = 0;
	fprintf( ioMEAN, " v0 J=0,3");
	fprintf(ioMEAN,"%7.3f", log10(SDIV(H2_X_colden[iVibHi][0])));
	fprintf(ioMEAN,"%7.3f", log10(SDIV(H2_X_colden[iVibHi][1])));
	fprintf(ioMEAN,"%7.3f", log10(SDIV(H2_X_colden[iVibHi][2])));
	fprintf(ioMEAN,"%7.3f", log10(SDIV(H2_X_colden[iVibHi][3])));

	return;
}


/*H2_ReadTransprob read transition probabilities */
void diatomics::H2_ReadTransprob( long int nelec, TransitionList &trns )
{
	const string cdDATAFILE[N_ELEC] = 
	{
		"transprob_X.dat",
		"transprob_B.dat", 
		"transprob_C_plus.dat",
		"transprob_C_minus.dat", 
		"transprob_B_primed.dat", 
		"transprob_D_plus.dat",
		"transprob_D_minus.dat" 
	};

	DEBUG_ENTRY( "H2_ReadTransprob()" );

	/* now open the data file */
	string chPath = path + cpu.i().chDirSeparator() + cdDATAFILE[nelec];
	FILE *ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	string chLine;
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " H2_ReadTransprob could not read first line of %s\n", cdDATAFILE[nelec].c_str() );
		cdEXIT(EXIT_FAILURE);
	}
	long i = 1;
	bool lgEOL;
	/* magic number */
	long n1 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n2 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n3 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of level1.dat 01 08 10 */
	if( ( n1 != 2 ) || ( n2 != 4 ) || ( n3 != 29 ) )
	{
		fprintf( ioQQQ, 
			 " H2_ReadTransprob: the version of %s is not the current version.\n",
			 cdDATAFILE[nelec].c_str() );
		fprintf( ioQQQ, 
			" I expected to find the number 2 4 29 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		cdEXIT(EXIT_FAILURE);
	}

	while( read_whole_line( chLine, ioDATA ) )
	{
		if( chLine.length() == 0 || chLine[0]=='\n' || chLine[0]=='\r' || chLine[0]==' ' )
			break;

		/* skip comment */
		if( chLine[0]=='#' )
			continue;

		double Aul;
		long int iVibHi , iVibLo , iRotHi , iRotLo , iElecHi , iElecLo;
		if( sscanf(chLine.c_str(),"%li\t%li\t%li\t%li\t%li\t%li\t%le", 
			   &iElecHi , &iVibHi ,&iRotHi , &iElecLo , &iVibLo , &iRotLo , &Aul ) != 7 )
		{
			fprintf( ioQQQ, "failed to read correct number of data values from %s\n", chPath.c_str() );
			cdEXIT(EXIT_FAILURE);
		}
		ASSERT( iElecHi == nelec );
		ASSERT( iElecHi < N_ELEC );
		ASSERT( iElecLo < N_ELEC );

		/* check that we actually included the levels in the model representation */
		if( iVibHi <= nVib_hi[iElecHi] && 
		    iVibLo <= nVib_hi[iElecLo] && 
			iRotHi <= nRot_hi[iElecHi][iVibHi] && 
			iRotLo <= nRot_hi[iElecLo][iVibLo])
		{
			long ipHi = ipEnergySort[iElecHi][iVibHi][iRotHi];
			long ipLo = ipEnergySort[iElecLo][iVibLo][iRotLo];
			double ener = states[ipHi].energy().WN() - states[ipLo].energy().WN();
			long lineIndex = ipTransitionSort[ipHi][ipLo];

			if( lgH2_radiative[ipHi][ipLo] )
			{
				// Could be multiple components (e.g., E2 and M1) for the same transition.
				ASSERT( trns[lineIndex].hasEmis() );
				trns[lineIndex].Emis().Aul() += (realnum)Aul;	
			}
			else
			{
				/* only lines that have real Aul are added to stack.  */
				trns[lineIndex].AddLine2Stack();
				trns[lineIndex].Emis().Aul() = (realnum)Aul;	

				/* say that this line exists */
				lgH2_radiative[ipHi][ipLo] = true;
			}

			/* prints transitions with negative energies  -  should not happen */
			if( ener <= 0. )
			{
				fprintf(ioQQQ,"negative energy H2 transition\t%li\t%li\t%li\t%li\t%.2e\t%.2e\n", 
					iVibHi,iVibLo,iRotHi,iRotLo,Aul,ener);
				ShowMe();
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	fclose( ioDATA );
	return;
}

#if 0
/*H2_Read_Cosmicray_distribution read distribution function for H2 population following cosmic ray collisional excitation */
void H2_Read_Cosmicray_distribution(void)
{
	static const bool CR_PRINT = false;
	static const int CR_X = 1, CR_VIB = 15, CR_J = 10, CR_EXIT = 3; 

	/*>>refer	H2	cr excit	Tine, S., Lepp, S., Gredel, R., & Dalgarno, A. 1997, ApJ, 481, 282 */
	DEBUG_ENTRY( "H2_Read_Cosmicray_distribution()" );

	/* now open the data file */
	string chPath = path + cpu.i().chDirSeparator() + "H2_CosmicRay_collision.dat";
	FILE *ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	string chLine;
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " H2_Read_Cosmicray_distribution could not read first line of %s\n",
			 "H2_Cosmic_collision.dat" );
		cdEXIT(EXIT_FAILURE);
	}

	long i = 1;
	bool lgEOL;	
	/* magic number */
	long n1 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n2 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n3 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of H2_Cosmic_collision.dat 01 21 03 */
	if( ( n1 != 1 ) || ( n2 != 21 ) || ( n3 != 3 ) )
	{
		fprintf( ioQQQ, 
			" H2_Read_Cosmicray_distribution: the version of %s is not the current version.\n",
			 "H2_Cosmic_collision.dat" );
		fprintf( ioQQQ, 
			" I expected to find the number 1 21 3 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		cdEXIT(EXIT_FAILURE);
	}

	long iRot = 1;
	long iVib = 1;
	long neut_frac = 0;
	while( iVib >= 0 )
	{
		/* read until not a comment */
		do
		{
			if( !read_whole_line( chLine, ioDATA ) )
				BadRead();
		}
		while( chLine[0] == '#' );

		long int j_minus_ji;
		double a[10];

		sscanf( chLine.c_str(),"%li\t%li\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", 
			&iVib ,&j_minus_ji , &a[0],&a[1],&a[2],&a[3],&a[4],&a[5],&a[6],&a[7],&a[8],&a[9] );
		/* negative iVib says end of data */
		if( iVib < 0 )
			continue;

		/* cr_rate[CR_X][CR_VIB][CR_J][CR_EXIT];*/
		/* check that we actually included the levels in the model representation */
		ASSERT( iVib < CR_VIB );
		ASSERT( j_minus_ji == -2 || j_minus_ji == +2 || j_minus_ji == 0 );
		ASSERT( neut_frac < CR_X );

		/* now make i_minus_ji an array index */
		j_minus_ji = 1 + j_minus_ji/2;
		ASSERT( j_minus_ji>=0 && j_minus_ji<=2 );

		/* option to add Gaussian random mole */
		for( iRot=0; iRot<CR_J; ++iRot )
		{
			cr_rate[neut_frac][iVib][iRot][j_minus_ji] = (realnum)a[iRot];
		}
		if( lgH2_NOISECOSMIC )
		{
			realnum r = realnum(xMeanNoise + ran.normal()*xSTDNoise);

			for( iRot=0; iRot<CR_J; ++iRot )
			{
				cr_rate[neut_frac][iVib][iRot][j_minus_ji] *= exp10f(r);
			}
		}

		if( CR_PRINT )
		{
			fprintf(ioQQQ,"cr rate\t%li\t%li", iVib , j_minus_ji ); 
			for( iRot=0; iRot<CR_J; ++iRot )
			{ 
				fprintf(ioQQQ,"\t%.3e", cr_rate[neut_frac][iVib][iRot][j_minus_ji] );
			} 
			fprintf(ioQQQ,"\n" );
		}
	}
	fclose( ioDATA );

	return;
}
#endif

class level_tmp
{
public:
	bool operator<( const level_tmp& second ) const
	{
		if( eWN < second.eWN )
			return true;
		else 
			return false; 
	}
	long n, v, J;
	double eWN;
};

/*H2_ReadEnergies read energies for all electronic levels */
void diatomics::H2_ReadEnergies( )
{
	DEBUG_ENTRY( "H2_ReadEnergies()" );

	vector<int> n, v, J;
	vector<double>eWN;
	
	for( long nelec=0; nelec<n_elec_states; ++nelec )
	{
		/* get energies out of files */
		H2_ReadEnergies(nelec,n,v,J,eWN);
	}
	
	vector<level_tmp> levels;
	levels.resize( n.size() );
	ASSERT( levels.size() > 0 );
	for( unsigned i = 0; i < n.size(); ++i )
	{
		levels[i].n = n[i];
		levels[i].v = v[i];
		levels[i].J = J[i];
		levels[i].eWN = eWN[i];
	}

	// now get levels in energy order (comparison done by operator< of level_tmp class) 
	sort( levels.begin(), levels.end() );

	// populate states
	states.init(sp->label.c_str(),0);
	for( vector<level_tmp>::iterator lev = levels.begin(); lev != levels.end(); ++lev )
	{
		states.addone( );
		long i = states.size() - 1;
		states[i].n() = lev->n;
		states[i].v() = lev->v;
		states[i].J() = lev->J;
		states[i].energy().set( lev->eWN, "cm^-1" );
		/* NB this must be kept parallel with nelem and ionstag in transition struc,
		 * since that struc expects to find the abundances here - abund set in hmole.c */
		states[i].nelem() = -1;
		/* this does not mean anything for a molecule */
		states[i].IonStg() = -1;
		ostringstream oss;
		oss << "n=" << lev->n << ',' << "v=" << lev->v << ',' << "J=" << lev->J;
		states[i].chConfig() = oss.str();
	}

	ASSERT( states.size() > 0 );
	ASSERT( states.size() == levels.size() );

	for( long nelec=0; nelec<n_elec_states; ++nelec )
	{
		ASSERT( nLevels_per_elec[nelec] > 0 );
		ASSERT( nVib_hi[nelec] > 0 );
		ASSERT( nVib_hi[nelec] > Jlowest[nelec] );

		nRot_hi[nelec].resize( nVib_hi[nelec]+1 );
		nRot_hi[nelec] = 0;
	}

	for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
	{
		nRot_hi[ (*st).n() ][ (*st).v() ] = MAX2( nRot_hi[ (*st).n() ][ (*st).v() ], (*st).J() );	
	}

	return;
}

void diatomics::H2_ReadEnergies( long int nelec, vector<int>& n, vector<int>& v, vector<int>&J, vector<double>& eWN )
{
	DEBUG_ENTRY( "diatomics::H2_ReadEnergies()" );

	const string cdDATAFILE[N_ELEC] = 
	{
		"energy_X.dat",
		"energy_B.dat", 
		"energy_C_plus.dat",
		"energy_C_minus.dat", 
		"energy_B_primed.dat", 
		"energy_D_plus.dat",
		"energy_D_minus.dat"
	};
	/* now open the data file */
	string chPath= path + cpu.i().chDirSeparator() + cdDATAFILE[nelec];
	FILE *ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	string chLine;
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " H2_ReadEnergies could not read first line of %s\n", cdDATAFILE[nelec].c_str() );
		cdEXIT(EXIT_FAILURE);
	}
	long i = 1;
	bool lgEOL;
	/* magic number */
	long n1 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n2 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n3 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of level1.dat 01 08 10 */
	if( ( n1 != 2 ) || ( n2 != 4 ) || ( n3 != 29 ) )
	{
		fprintf( ioQQQ, 
			 " H2_ReadEnergies: the version of %s is not the current version.\n",
			 cdDATAFILE[nelec].c_str() );
		fprintf( ioQQQ, 
			" I expected to find the number 2 4 29 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		cdEXIT(EXIT_FAILURE);
	}

	/* this will count the number of levels within each electronic state */
	nLevels_per_elec[nelec] = 0;
	nVib_hi[nelec] = 0;
	Jlowest[nelec] = LONG_MAX;
	
	while( read_whole_line( chLine, ioDATA ) )
	{
		if( chLine.length() == 0 || chLine[0] == '\n' || chLine[0] == '\r' || chLine[0] == ' ' )
			break;
		/* skip comment */
		if( chLine[0] == '#' )
			continue;
		long iVib, iRot;
		double energyWN;
		if( sscanf( chLine.c_str(), "%li\t%li\t%le", &iVib, &iRot, &energyWN ) != 3 )
		{
			fprintf( ioQQQ, "failed to read correct number of data values from %s\n", chPath.c_str() );
			cdEXIT(EXIT_FAILURE);
		}
		ASSERT( iVib >= 0 );
		ASSERT( iRot >= 0 );
		ASSERT( energyWN > 0. || (nelec==0 && iVib==0 && iRot==0 ) );

		n.push_back( nelec );
		v.push_back( iVib );
		J.push_back( iRot );
		eWN.push_back( energyWN );

		// update limits		
		nVib_hi[nelec] = MAX2( nVib_hi[nelec], iVib );
		Jlowest[nelec] = MIN2( Jlowest[nelec], iRot );
		/* increment number of levels within this electronic state */
		++nLevels_per_elec[nelec];
	}

	ASSERT( n.size() > 0 );
	ASSERT( nLevels_per_elec[nelec] > 0 );
	ASSERT( nVib_hi[nelec] > 0 );
	ASSERT( nVib_hi[nelec] > Jlowest[nelec] );
	
	fclose( ioDATA );

	return;
}

/*H2_ReadDissocEnergies read energies for all electronic levels */
void diatomics::H2_ReadDissocEnergies( void )
{
	const string cdDATAFILE = "energy_dissoc.dat";

	DEBUG_ENTRY( "H2_ReadDissocEnergies()" );

	/* now open the data file */
	string chPath = path + cpu.i().chDirSeparator() + cdDATAFILE;
	FILE *ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	string chLine;
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " H2_ReadDissocEnergies could not read first line of %s\n", cdDATAFILE.c_str() );
		cdEXIT(EXIT_FAILURE);
	}
	long i = 1;
	bool lgEOL;
	/* magic number */
	long n1 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n2 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n3 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of level1.dat 01 08 10 */
	if( ( n1 != 2 ) || ( n2 != 4 ) || ( n3 != 29 ) )
	{
		fprintf( ioQQQ, 
			" H2_ReadDissocEnergies: the version of %s is not the current version.\n", cdDATAFILE.c_str() );
		fprintf( ioQQQ, 
			" I expected to find the number 2 4 29 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		cdEXIT(EXIT_FAILURE);
	}

	while( read_whole_line( chLine, ioDATA ) )
	{
		if( chLine.length() == 0 || chLine[0] == '\n' || chLine[0] == '\r' || chLine[0] == ' ' )
			break;
		/* skip comment */
		if( chLine[0] == '#' )
			continue;
		long iElec;
		double energyWN;
		if( sscanf(chLine.c_str(),"%li\t%le", &iElec, &energyWN ) != 2 )
		{
			fprintf( ioQQQ, "failed to read correct number of data values from %s\n", chPath.c_str() );
			cdEXIT(EXIT_FAILURE);
		}
		ASSERT( iElec >= 0 );
		ASSERT( iElec < N_ELEC );
		ASSERT( energyWN > 0. );
		H2_DissocEnergies[iElec] = energyWN;
	}
	fclose( ioDATA );

	return;
}

/*H2_ReadDissprob read dissociation probabilities and kinetic energies for all electronic levels */
void diatomics::H2_ReadDissprob( long int nelec )
{
	const string cdDATAFILE[N_ELEC] = 
	{
		"dissprob_X.dat",/* this does not exist and nelec == 0 is not valid */
		"dissprob_B.dat", 
		"dissprob_C_plus.dat",
		"dissprob_C_minus.dat", 
		"dissprob_B_primed.dat", 
		"dissprob_D_plus.dat",
		"dissprob_D_minus.dat"
	};

	DEBUG_ENTRY( "H2_ReadDissprob()" );

	ASSERT( nelec > 0 );

	/* now open the data file */
	string chPath = path + cpu.i().chDirSeparator() + cdDATAFILE[nelec];
	FILE *ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	string chLine;
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " H2_ReadDissprob could not read first line of %s\n", cdDATAFILE[nelec].c_str() );
		cdEXIT(EXIT_FAILURE);
	}
	long i = 1;
	bool lgEOL;
	/* magic number */
	long n1 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n2 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n3 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of level1.dat 01 08 10 */
	if( ( n1 != 3 ) || ( n2 != 2 ) || ( n3 != 11 ) )
	{
		fprintf( ioQQQ, 
			" H2_ReadDissprob: the version of %s is not the current version.\n",
			 cdDATAFILE[nelec].c_str() );
		fprintf( ioQQQ, 
			" I expected to find the number 3 2 11 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		cdEXIT(EXIT_FAILURE);
	}

	while( read_whole_line( chLine, ioDATA ) )
	{
		if( chLine.length() == 0 || chLine[0] == '\n' || chLine[0] == '\r' || chLine[0] == ' ' )
			break;
		/* skip comment */
		if( chLine[0] == '#' )
			continue;

		long iVib, iRot;
		double a, b;
		i = 1;
		sscanf(chLine.c_str(),"%li\t%li\t%le\t%le", 
			&iVib, &iRot, 
			/* dissociation probability */
			&a ,
			/* dissociation kinetic energy - eV not ergs */
			&b);

		/* these have to agree if data file is valid */
		//ASSERT( iVib >= 0 );
		//ASSERT( iVib <= nVib_hi[nelec] );
		//ASSERT( iRot >= Jlowest[nelec] );
		//ASSERT( iRot <= nRot_hi[nelec][iVib] );
		if( ( iVib < 0 ) ||
			( iVib > nVib_hi[nelec] ) ||
			( iRot < Jlowest[nelec] ) ||
			( iRot > nRot_hi[nelec][iVib] ) )
			continue;

		/* dissociation probability */
		H2_dissprob[nelec][iVib][iRot] = (realnum)a;
		/* dissociation kinetic energy - eV not ergs */
		H2_disske[nelec][iVib][iRot] = (realnum)b;
	}
	fclose( ioDATA );
	return;
}


/*H2_Read_hminus_distribution read distribution function for H2 population following formation from H minus */
void diatomics::H2_Read_hminus_distribution(void)
{
	/* set true for lots of printout */
	const bool lgH2HMINUS_PRT = false;

	DEBUG_ENTRY( "H2_Read_hminus_distribution()" );

	/* now open the data file */
	string chPath = path + cpu.i().chDirSeparator() + "hminus_deposit.dat";
	FILE *ioDATA = open_data( chPath, "r" );

	/* read the first line and check that magic number is ok */
	string chLine;
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " H2_Read_hminus_distribution could not read first line of %s\n", chPath.c_str() );
		cdEXIT(EXIT_FAILURE);
	}

	long i = 1;
	bool lgEOL;
	/* magic number */
	long n1 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n2 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
	long n3 = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of H2_hminus_deposit.dat 01 08 10 */
	if( ( n1 != 2 ) || ( n2 != 10 ) || ( n3 != 17 ) )
	{
		fprintf( ioQQQ, 
			" H2_Read_hminus_distribution: the version of %s is not the current version.\n",
			 chPath.c_str() );
		fprintf( ioQQQ, 
			" I expected to find the number 2 10 17 and got %li %li %li instead.\n" ,
			n1 , n2 , n3 );
		cdEXIT(EXIT_FAILURE);
	}

	long iRot = 1;
	long iVib = 1;
	double sumrate[nTE_HMINUS] = {0};
	while( iVib >= 0 )
	{
		/* read until not a comment */
		do
		{
			if( !read_whole_line( chLine, ioDATA ) )
				BadRead();
		}
		while( chLine[0] == '#' );

		double a[nTE_HMINUS], ener;
		sscanf( chLine.c_str(), "%li\t%li\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", 
			&iVib, &iRot, &ener, &a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6] );
		/* negative iVib says end of data */
		if( iVib < 0 )
			continue;

		/* check that we actually included the levels in the model representation */
		ASSERT( iVib <= nVib_hi[0] && 
			iRot <= nRot_hi[0][iVib] );

		if( lgH2HMINUS_PRT )
			fprintf(ioQQQ,"hminusss\t%li\t%li", iVib , iRot );
		for( i=0; i<nTE_HMINUS; ++i )
		{
			H2_X_hminus_formation_distribution[i][iVib][iRot] = (realnum)exp10(-a[i]);
			sumrate[i] += H2_X_hminus_formation_distribution[i][iVib][iRot];
			if( lgH2HMINUS_PRT )
				fprintf(ioQQQ,"\t%.3e", H2_X_hminus_formation_distribution[i][iVib][iRot] );
		}
		if( lgH2HMINUS_PRT )
			fprintf(ioQQQ,"\n" );
	}
	fclose( ioDATA );

	if( lgH2HMINUS_PRT )
	{
		/* print total rate */
		fprintf(ioQQQ," total H- formation rate ");
		/* convert temps to log */
		for(i=0; i<nTE_HMINUS; ++i )
		{
			fprintf(ioQQQ,"\t%.3e" , sumrate[i]);
		}
		fprintf(ioQQQ,"\n" );
	}

	/* convert to dimensionless factors that add to unity */
	for( iVib=0; iVib<=nVib_hi[0]; ++iVib )
	{
		for( iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
		{
			for(i=0; i<nTE_HMINUS; ++i )
			{
				H2_X_hminus_formation_distribution[i][iVib][iRot] /= (realnum)sumrate[i];
			}
		}
	}

	if( lgH2HMINUS_PRT )
	{
		/* print total rate */
		fprintf(ioQQQ,"  H- distribution function ");
		for( iVib=0; iVib<=nVib_hi[0]; ++iVib )
		{
			for( iRot=Jlowest[0]; iRot<=nRot_hi[0][iVib]; ++iRot )
			{
				fprintf(ioQQQ,"%li\t%li", iVib , iRot );
				for(i=0; i<nTE_HMINUS; ++i )
				{
					fprintf(ioQQQ,"\t%.3e", H2_X_hminus_formation_distribution[i][iVib][iRot] );
				}
				fprintf(ioQQQ,"\n" );
			}
		}
	}
	return;
}

/* ===================================================================== */
/*H2_Punch_line_data save line data for H2 molecule */
void diatomics::H2_Punch_line_data(
	/* io unit for save */
	FILE* ioPUN ,
	/* save all levels if true, only subset if false */
	bool lgDoAll )
{
	if( !lgEnabled )
		return;

	DEBUG_ENTRY( "H2_Punch_line_data()" );

	if( lgDoAll )
	{
		fprintf( ioQQQ, 
			" H2_Punch_line_data ALL option not implemented in H2_Punch_line_data yet 1\n" );
		cdEXIT(EXIT_FAILURE);
	}
	else
	{
		fprintf( ioPUN, "#Eu\tVu\tJu\tEl\tVl\tJl\tWL\tgl\tgu\tgf\tA\tCS\tn(crt)\n" );
		/* save line date, looping over all possible lines */
		for( TransitionList::iterator tr = trans.begin(); tr != trans.end(); ++tr )
		{
			if( (*tr).ipCont() <= 0 )
				continue;
			(*tr).Coll().col_str() = 0.;
			qList::iterator Hi = ( (*tr).Hi() );	
			qList::iterator Lo = ( (*tr).Lo() );	
			/* print quantum indices */
			fprintf(ioPUN,"%2li\t%2li\t%2li\t%2li\t%2li\t%2li\t",
				(*Hi).n(), (*Hi).v(), (*Hi).J(),
				(*Lo).n(), (*Lo).v(), (*Lo).J() );
			Save1LineData( *tr, ioPUN, false );
		}

		fprintf( ioPUN , "\n");
	}
	return;
}

/*H2_PunchLineStuff include H2 lines in punched optical depths, etc, called from SaveLineStuff */
void diatomics::H2_PunchLineStuff( FILE * io , realnum xLimit  , long index)
{
	if( !lgEnabled )
		return;

	DEBUG_ENTRY( "H2_PunchLineStuff()" );

	/* loop over all possible lines */
	for( TransitionList::iterator tr = trans.begin(); tr != trans.end(); ++tr )
	{
		if( (*tr).ipCont() <= 0 )
			continue;
		Save1Line( *tr, io, xLimit, index, GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN]));
	}

	return;
}


/*chMolBranch returns a char with the spectroscopic branch of a transition */
STATIC char chMolBranch( long iRotHi , long int iRotLo )
{
	/* these are the spectroscopic branches */
	char chBranch[5] = {'O','P','Q','R','S'};
	/* this is the index within the chBranch array */
	int ip = 2 + (iRotHi - iRotLo);
	if( ip<0 || ip>=5 )
	{
		fprintf(ioQQQ," chMolBranch called with insane iRotHi=%li iRotLo=%li ip=%i\n",
			iRotHi , iRotLo , ip );
		ip = 0;
	}

	return( chBranch[ip] );
}

/*H2_PunchDo save some properties of the large H2 molecule */
void diatomics::H2_PunchDo( FILE* io ,  char chJOB[] , const char chTime[] , long int ipPun )
{
	DEBUG_ENTRY( "H2_PunchDo()" );

	if( !lgEnabled )
	{
		fprintf( io, "%s model is not enabled. This save command is therefore disabled.\n", label.c_str() );
		return;
	}

	/* which job are we supposed to do? This routine is active even when H2 is not turned on
	 * so do not test on lgEnabled initially */

	/* H2 populations computed in last zone - 
	 * give all of molecule in either matrix or triplet format */
	if( (strcmp( chJOB , "H2po" ) == 0) && (strcmp(chTime,"LAST") == 0) &&
		(save.punarg[ipPun][2] != 0) )
	{
		/* >>chng 04 feb 19, do not save if H2 not yet evaluated */
		if( lgEnabled && lgEvaluated )
		{
			long iVibHi= 0;
			long iRotHi = 0;
			long iElecHi=0;
			long LimVib, LimRot;
			/* the limit to the number of vibration levels punched -
			* default is all, but first two numbers on save h2 pops command
			* reset limit */
			/* this is limit to vibration */
			if( save.punarg[ipPun][0] > 0 )
			{
				LimVib = (long)save.punarg[ipPun][0];
			}
			else
			{
				LimVib = nVib_hi[iElecHi];
			}

			/* first save the current ortho, para, and total H2 density */
			fprintf(io,"%i\t%i\t%.3e\tortho\n", 
				103 , 
				103 ,
				ortho_density );
			fprintf(io,"%i\t%i\t%.3e\tpara\n", 
				101 , 
				101 ,
				para_density );
			fprintf(io,"%i\t%i\t%.3e\ttotal\n", 
				0 , 
				0 ,
				(*dense_total) );

			/* now save the actual populations, first part both matrix and triplets */
			for( iVibHi=0; iVibHi<=LimVib; ++iVibHi )
			{
				/* this is limit to rotation quantum index */
				if( save.punarg[ipPun][1] > 0 )
				{
					LimRot = (long)MIN2(
						save.punarg[ipPun][1] , (realnum)nRot_hi[iElecHi][iVibHi]);
				}
				else
				{
					LimRot = nRot_hi[iElecHi][iVibHi];
				}
				if( save.punarg[ipPun][2] > 0 )
				{
					long int i;
					/* this option save matrix */
					if( iVibHi == 0 )
					{
						fprintf(io,"vib\\rot");
						/* this is first vib, so make row of rot numbs */
						for( i=0; i<=LimRot; ++i )
						{
							fprintf(io,"\t%li",i);
						}
						fprintf(io,"\n");
					}
					fprintf(io,"%li",iVibHi );
					for( iRotHi=Jlowest[iElecHi]; iRotHi<=LimRot; ++iRotHi )
					{
						fprintf(io,"\t%.3e", 
							states[ ipEnergySort[iElecHi][iVibHi][iRotHi] ].Pop()/(*dense_total) );
					}
					fprintf(io,"\n" );
				}
				else if( save.punarg[ipPun][2] < 0 )
				{
					/* this option save triplets - the default */
					for( iRotHi=Jlowest[iElecHi]; iRotHi<=LimRot; ++iRotHi )
					{
						/* this will say whether ortho or para,
 						 * H2_lgOrtho is 0 or 1 depending on whether or not ortho, 
						 * so chlgPara[H2_lgOrtho] gives P or O for printing */
						const char chlgPara[2]={'P','O'};
						const long ipHi = ipEnergySort[iElecHi][iVibHi][iRotHi];

						/* intensity, relative to normalization line, for faintest line to save */
						fprintf(io,"%li\t%li\t%c\t%.1f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", 
							/* upper vibration and rotation quantum numbers */
							iVibHi , iRotHi ,
							/* an 'O' or 'P' for ortho or para */
							chlgPara[H2_lgOrtho[iElecHi][iVibHi][iRotHi]],
							/* the level excitation energy in wavenumbers */
							states[ipHi].energy().WN(),
							/* actual population relative to total H2 */
							states[ipHi].Pop()/(*dense_total) ,
							/* old level populations for comparison */
							H2_old_populations[iElecHi][iVibHi][iRotHi]/(*dense_total) ,
							/* populations per h2 and per statistical weight */
						   states[ipHi].Pop()/(*dense_total)/states[ipHi].g() ,
							/* LTE departure coefficient */
							/* >>chng 05 jan 26, missing factor of H2 abundance LTE is norm to unity, not tot abund */
							states[ipHi].Pop()/SDIV(H2_populations_LTE[iElecHi][iVibHi][iRotHi]*(*dense_total) ) ,
							/* fraction of exits that were collisional */
							H2_col_rate_out[iVibHi][iRotHi]/SDIV(H2_col_rate_out[iVibHi][iRotHi]+H2_rad_rate_out[0][iVibHi][iRotHi]) ,
							/* fraction of entries that were collisional */
							H2_col_rate_in[iVibHi][iRotHi]/SDIV(H2_col_rate_in[iVibHi][iRotHi]+H2_rad_rate_in[iVibHi][iRotHi]),
							/* collisions out */
							H2_col_rate_out[iVibHi][iRotHi],
							/* radiation out */
							H2_rad_rate_out[0][iVibHi][iRotHi] ,
							/* radiation out */
							H2_col_rate_in[iVibHi][iRotHi],
							/* radiation in */
							H2_rad_rate_in[iVibHi][iRotHi]
							);
					}
				}
			}
		}
	}
	/* save H2 populations for each zone 
	 * populations of v=0 for each zone */
	else if( (strcmp( chJOB , "H2po" ) == 0) && (strcmp(chTime,"LAST") != 0) &&
		(save.punarg[ipPun][2] == 0) )
	{
		/* >>chng 04 feb 19, do not save if h2 not yet evaluated */
		if( lgEnabled && lgEvaluated )
		{
			fprintf(io,"%.5e\t%.3e\t%.3e", radius.depth_mid_zone , 
				ortho_density , para_density);
			/* rel pops of first two excited electronic states */
			fprintf(io,"\t%.3e\t%.3e", 
				pops_per_elec[1] , pops_per_elec[2]);
			long iElecHi = 0;
			long iVibHi = 0;
			long LimVib, LimRot;
			/* this is limit to vibration quantum index */
			if( save.punarg[ipPun][0] > 0 )
			{
				LimVib = (long)save.punarg[ipPun][1];
			}
			else
			{
				LimVib = nRot_hi[iElecHi][iVibHi];
			}
			LimVib = MIN2( LimVib , nVib_hi[iElecHi] );
			/* this is limit to rotation quantum index */
			if( save.punarg[ipPun][1] > 0 )
			{
				LimRot = (long)save.punarg[ipPun][1];
			}
			else
			{
				LimRot = nRot_hi[iElecHi][iVibHi];
			}
			for( iVibHi = 0; iVibHi<=LimVib; ++iVibHi )
			{
				fprintf(io,"\tv=%li",iVibHi);
				long int LimRotVib = MIN2( LimRot , nRot_hi[iElecHi][iVibHi] );
				for( long iRotHi=Jlowest[iElecHi]; iRotHi<=LimRotVib; ++iRotHi )
				{
					fprintf(io,"\t%.3e", 
						states[ ipEnergySort[iElecHi][iVibHi][iRotHi] ].Pop()/(*dense_total) );
				}
			}
			fprintf(io,"\n");
		}
	}

	/* save column densities */
	else if( (strcmp( chJOB , "H2cl" ) == 0) && (strcmp(chTime,"LAST") == 0) )
	{
		long iVibHi= 0;
		long iRotHi = 0;
		long iElecHi=0;
		long LimVib, LimRot;
		/* the limit to the number of vibration levels punched -
		 * default is all, but first two numbers on save h2 pops command
		 * reset limit */
		/* this is limit to vibration */
		if( save.punarg[ipPun][0] > 0 )
		{
			LimVib = (long)save.punarg[ipPun][0];
		}
		else
		{
			LimVib = nVib_hi[iElecHi];
		}

		/* first save ortho and para populations */
		fprintf(io,"%i\t%i\t%.3e\tortho\n", 
			103 , 
			103 ,
			ortho_colden );
		fprintf(io,"%i\t%i\t%.3e\tpara\n", 
			101 , 
			101 ,
			para_colden );
		/* total H2 column density */
		fprintf(io,"%i\t%i\t%.3e\ttotal\n", 
			0 , 
			0 ,
			mole.species[sp->index].column + (sp_star->index > 0 ? mole.species[sp_star->index].column : 0) );

		/* save level column densities */
		for( iVibHi=0; iVibHi<=LimVib; ++iVibHi )
		{
		if( lgEnabled )
		{
			/* this is limit to rotation quantum index */
			if( save.punarg[ipPun][1] > 0 )
			{
				LimRot = (long)save.punarg[ipPun][1];
			}
			else
			{
				LimRot = nRot_hi[iElecHi][iVibHi];
			}
			if( save.punarg[ipPun][2] > 0 )
			{
				long int i;
				/* save matrix */
				if( iVibHi == 0 )
				{
					fprintf(io,"vib\\rot");
					/* this is first vib, so make row of rot numbs */
					for( i=0; i<=LimRot; ++i )
					{
						fprintf(io,"\t%li",i);
					}
					fprintf(io,"\n");
				}
				fprintf(io,"%li",iVibHi );
				for( iRotHi=Jlowest[iElecHi]; iRotHi<=LimRot; ++iRotHi )
				{
					fprintf(io,"\t%.3e", 
						H2_X_colden[iVibHi][iRotHi]/(*dense_total) );
				}
				fprintf(io,"\n" );
			}
			else
			{
				/* save triplets - the default */
				for( iRotHi=Jlowest[iElecHi]; iRotHi<=LimRot; ++iRotHi )
				{
					fprintf(io,"%li\t%li\t%.1f\t%.3e\t%.3e\t%.3e\t%.3e\n", 
						iVibHi , 
						iRotHi ,
						/* energy relative to 0,0, T1CM converts wavenumber to K */
						states[ ipEnergySort[iElecHi][iVibHi][iRotHi] ].energy().K(),
						/* these are column densities for actual molecule */
						H2_X_colden[iVibHi][iRotHi] ,
						H2_X_colden[iVibHi][iRotHi]/states[ ipEnergySort[iElecHi][iVibHi][iRotHi] ].g() ,
						/* these are same column densities but for LTE populations */
						H2_X_colden_LTE[iVibHi][iRotHi] ,
						H2_X_colden_LTE[iVibHi][iRotHi]/states[ ipEnergySort[iElecHi][iVibHi][iRotHi] ].g());
				}
			}
		}
		}
	}
	else if( (strcmp(chJOB , "H2pd" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* save PDR 
		 * output some PDR information (densities, rates) for each zone */
		fprintf(io,"%.5e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n", 
			/* depth in cm */
			radius.depth_mid_zone ,
			/* the computed ortho and para densities */
			ortho_density , 
			para_density ,
			/* the Lyman Werner band dissociation, Tielens & Hollenbach */
			hmi.H2_Solomon_dissoc_rate_TH85_H2g , 
			/* the Lyman Werner band dissociation, Bertoldi & Draine */
			hmi.H2_Solomon_dissoc_rate_BD96_H2g,
			/* the Lyman Werner band dissociation, big H2 mole */
			Solomon_dissoc_rate_g);
	}
	else if( (strcmp(chJOB , "H2cm" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* save H2 cooling per particle */
		fprintf(io,"%.5e\t%.5e\t%.5e\n", 
			/* temperature in zone */
			phycon.te,
			/* LTE cooling */
			LTE_Cooling_per_H2(),
			/* net H2 cooling per particle */
			-HeatDexc / Abund()
			);
	}

	else if( (strcmp(chJOB , "H2co" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* save H2 cooling - do heating cooling for each zone old new H2 */
		fprintf(io,"%.5e\t%.5e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.5e\n", 
			/* depth in cm */
			radius.depth_mid_zone ,
			/* temperature in zone */
			phycon.te,
			/* total cooling, equal to total heating */
			thermal.ctot ,
			/* H2 destruction by Solomon process, TH85 rate */
			hmi.H2_Solomon_dissoc_rate_TH85_H2g,
			/* H2 destruction by Solomon process, big H2 model rate */
			Solomon_dissoc_rate_g +
				Solomon_dissoc_rate_s,
			/* H2 photodissociation heating, eqn A9 of Tielens & Hollenbach 1985a */
			hmi.HeatH2Dish_TH85,
			/* heating due to dissociation of electronic excited states */
			HeatDiss ,
			/* cooling (usually neg and so heating) due to collisions within X */
			hmi.HeatH2Dexc_TH85,
			HeatDexc,
			-HeatDexc / Abund()
			);
	}

	else if( (strcmp(chJOB , "H2le" ) == 0) && (strcmp(chTime,"LAST") == 0) )
	{
		/* save H2 levels */
		for( long int ipHi=0; ipHi < nLevels_per_elec[0]; ipHi++ )
		{
			long iRotHi = ipRot_H2_energy_sort[ipHi];
			long iVibHi = ipVib_H2_energy_sort[ipHi];
			long int nColl;
			double Asum , Csum[N_X_COLLIDER];
			Asum = 0;
			for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
				Csum[nColl] = 0.;
			for( long int ipLo=0; ipLo<ipHi; ++ipLo )
			{
				/* radiative decays down */
				if( lgH2_radiative[ipHi][ipLo] )
				{
					EmissionProxy em = trans[ ipTransitionSort[ipHi][ipLo] ].Emis();
					Asum += em.Aul() * ( em.Ploss() );
				}
				/* all collisions down */
				mr3ci H2cr = CollRateCoeff.begin(ipHi,ipLo);
				for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
					Csum[nColl] += H2cr[nColl];
			}

			/* save H2 level energies */
			fprintf(io,"%li\t%li\t%.2f\t%li\t%.3e", 
				iVibHi , iRotHi,
				states[ipHi].energy().WN(),
				(long)states[ipHi].g(),
				Asum );
			for( nColl=0; nColl<N_X_COLLIDER; ++nColl )
				/* sum over all lower levels */
				fprintf(io,"\t%.3e",Csum[nColl]);
			fprintf(io,"\n");
		}
	}

	else if( (strcmp(chJOB , "H2ra" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* save h2 rates - some rates and lifetimes */
		double sumpop = 0. , sumlife = 0.;

		/* this block, find lifetime against photo excitation into excited electronic states */
		if( lgEnabled && lgEvaluated )
		{
			/* only do if radiative transition exists */
			for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
			{
				qList::iterator Lo = ( (*tr).Lo() );	
				if( (*Lo).n() > 0 || (*Lo).v() > 0 )
					continue;
				sumlife += (*tr).Emis().pump() * (*(*tr).Lo()).Pop();
				sumpop += (*(*tr).Lo()).Pop(); 
			}
		}

		/* continue output from save h2 rates command */
		/* find photoexcitation rates from v=0 */
		/* PDR information for each zone */
		fprintf(io,
			"%.5e\t%.3e\t%.3e\t%.3e\t%li", 
			/* depth in cm */
			radius.depth_mid_zone ,
			/* the column density (cm^-2) in H2 */
			findspecieslocal("H2")->column+findspecieslocal("H2*")->column,
			/* this is a special form of column density - should be proportional
			 * to total shielding */
			colden.coldenH2_ov_vel ,
			/* visual extinction due to dust alone, of point source (star)*/
			rfield.extin_mag_V_point,
			/* number of large molecule evaluations in this zone */
			nCall_this_zone );
		fprintf(io,
			"\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e",
			/* total H2 fraction */
			(*dense_total)/dense.gas_phase[ipHYDROGEN] ,
			/* chemistry renorm factor */
			H2_renorm_chemistry,
			/* rate H2 forms on grains */
			gv.rate_h2_form_grains_used_total , 
			/* rate H2 forms by H minus route */
			findspecieslocal("H-")->den*1.35e-9,
			/* H2 destruction by Solomon process, TH85 rate */
			hmi.H2_Solomon_dissoc_rate_TH85_H2g + hmi.H2_Solomon_dissoc_rate_TH85_H2s,
			/* H2 destruction by Solomon process, Bertoldi & Draine rate */
			hmi.H2_Solomon_dissoc_rate_BD96_H2g + hmi.H2_Solomon_dissoc_rate_BD96_H2s,
			/* H2 destruction by Solomon process, Elwert et al. in preparation */
			hmi.H2_Solomon_dissoc_rate_ELWERT_H2g + hmi.H2_Solomon_dissoc_rate_ELWERT_H2g,
			/* H2 destruction by Solomon process, big H2 model rate */
			Solomon_dissoc_rate_g + Solomon_dissoc_rate_s,
			/* rate s-1 H2 electronic excit states decay into H2g */
			Solomon_elec_decay_g ,
			/* rate s-1 H2 electronic excit states decay into H2s */
			Solomon_elec_decay_s 
			);
		fprintf(io,
			"\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e", 
			/* The TH85 estimate of the radiation field relative to the Habing value */
			hmi.UV_Cont_rel2_Habing_TH85_depth,
			/* The DB96 estimate of the radiation field relative to the Habing value */
			hmi.UV_Cont_rel2_Draine_DB96_depth,
			/* cosmic ray ionization rate */
			secondaries.csupra[ipHYDROGEN][0]*0.93,
			sumlife/SDIV( sumpop ) ,
			hmi.H2_Solomon_dissoc_rate_BD96_H2g/SDIV(hmi.UV_Cont_rel2_Habing_TH85_depth) ,
			Solomon_dissoc_rate_g/SDIV(hmi.UV_Cont_rel2_Habing_TH85_depth),
			Solomon_dissoc_rate_g/SDIV(hmi.UV_Cont_rel2_Habing_spec_depth),
			hmi.H2_rate_destroy);
		fprintf(io,
			"\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\n",
			hmi.HeatH2Dish_TH85,
			hmi.HeatH2Dexc_TH85,
			HeatDiss, 
			HeatDexc,
			thermal.htot);
	}
	/* save h2 solomon */
	else if( (strcmp(chJOB , "H2so" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* remember as many as nSOL lines contributing to total Solomon process */
		const int nSOL = 100;
		double sum, one;
		long int jlosave[nSOL] , ivlosave[nSOL],
			iehisave[nSOL] ,jhisave[nSOL] , ivhisave[nSOL],
			nsave,
			ipOrdered[nSOL];
		int nFail;
		realnum fsave[nSOL], wlsave[nSOL];
		/* Solomon process, and where it came from */
		fprintf(io,"%.5e\t%.3e", 
			/* depth in cm */
			radius.depth_mid_zone ,
			/* H2 destruction by Solomon process, big H2 model rate */
			Solomon_dissoc_rate_g +
				Solomon_dissoc_rate_s);
		sum = 0.;
		/* find sum of all radiative exits from X into excited electronic states */
		if( lgEnabled && lgEvaluated )
		{
			for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
			{
				qList::iterator Lo = ( (*tr).Lo() );	
				if( (*Lo).n() > 0 )
					continue;
				sum += (*(*tr).Lo()).Pop() * (*tr).Emis().pump();	
			}

			/* make sure it is safe to div by sum */
			sum = SDIV( sum );
			nsave = 0;
			/* now loop back over X and print all those which contribute more than frac of the total */
			const double frac = 0.01;
			/* only do if radiative transition exists */
			for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
			{
				qList::iterator Lo = ( (*tr).Lo() );	
				if( (*Lo).n() > 0 )
					continue;
				one = (*(*tr).Lo()).Pop() * (*tr).Emis().pump();	
				if( one/sum > frac && nsave<nSOL)
				{
					qList::iterator Hi = ( (*tr).Hi() );	
					fsave[nsave] = (realnum)(one/sum);
					jlosave[nsave] = (*Lo).J();
					ivlosave[nsave] = (*Lo).v();
					jhisave[nsave] = (*Hi).J();
					ivhisave[nsave] = (*Hi).v();
					iehisave[nsave] = (*Hi).n();
					wlsave[nsave] = (*tr).WLAng();
					++nsave;
				}
			}
			
			/* now sort by decreasing importance */
			/*spsort netlib routine to sort array returning sorted indices */
			spsort(
				/* input array to be sorted */
				fsave, 
				/* number of values in x */
				nsave, 
				/* permutation output array */
				ipOrdered, 
				/* flag saying what to do - 1 sorts into increasing order, not changing
				* the original routine */
				-1, 
				/* error condition, should be 0 */
				&nFail);

			/* print ratio of pumps to dissociations - this is 9:1 in TH85 */
			/*>>chng 05 jul 20, TE, save average energy in H2s and renormalization factors for H2g and H2s */
			/* >>chng 05 sep 16, TE, chng denominator to do g and s with proper dissoc rates */
			fprintf(io,"\t%.3f\t%.3f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e",
				/* this is sum of photons and CRs */
				(sum + secondaries.csupra[ipHYDROGEN][0]*2.02f)/SDIV((Solomon_dissoc_rate_g * findspecieslocal("H2")->den +
					Solomon_dissoc_rate_s * findspecieslocal("H2*")->den) ), 
				/* this is sum of photons and CRs */
				(sum + secondaries.csupra[ipHYDROGEN][0]*2.02f) /SDIV((Solomon_dissoc_rate_g *  H2_den_g+
					Solomon_dissoc_rate_s * H2_den_s) ),
				average_energy_g, average_energy_s,	
				findspecieslocal("H2")->den/SDIV(H2_den_g), findspecieslocal("H2*")->den/SDIV(H2_den_s),
				H2_den_g/SDIV((*dense_total)), H2_den_s/SDIV((*dense_total))
				);
			for( long i=0; i<nsave; ++i )
			{
				long ip = ipOrdered[i];
				/*lint -e644 not init */
				fprintf(io,"\t%li\t%li\t%li\t%li\t%li\t%.3f\t%.3f", 
					iehisave[ip],ivhisave[ip],jhisave[ip],ivlosave[ip] , jlosave[ip] , fsave[ip] , wlsave[ip] );
				/*lint +e644 not init */
			}
			fprintf(io,"\n"); 
		}
	}

	else if( (strcmp(chJOB , "H2te" ) == 0) && (strcmp(chTime,"LAST") != 0) )
	{
		/* save h2 temperatures */
		double pop_ratio10,pop_ratio20,pop_ratio30,pop_ratio31,pop_ratio40;
		double T10,T20,T30,T31,T40;
		/* subscript"sum" denotes integrated quantities */
		double T10_sum,T20_sum,T30_sum,T31_sum,T40_sum;
		double pop_ratio20_sum,pop_ratio30_sum,pop_ratio31_sum,pop_ratio40_sum;
		if( lgEnabled && nCall_this_zone )
		{
			double pop0 = states[0].Pop();
			double pop1 = states[1].Pop();
			double pop2 = states[2].Pop();
			double pop3 = states[3].Pop();
			double pop4 = states[4].Pop();

			double energyK = states[1].energy().K() - states[0].energy().K();
			/* the ratio of populations of J=1 to 0 */
			pop_ratio10 = pop1/SDIV(pop0);
			double pop_ratio10_sum = H2_X_colden[0][1]/SDIV(H2_X_colden[0][0]);
			/* the corresponding temperature */
			T10 = -170.5/log(SDIV(pop_ratio10) * states[0].g()/states[1].g());
			T10_sum = -170.5/log(SDIV(pop_ratio10_sum) * states[0].g()/states[1].g());

			energyK = states[2].energy().K() - states[0].energy().K();
			pop_ratio20 = pop2/SDIV(pop0);
			T20 = -energyK/log(SDIV(pop_ratio20) * states[0].g()/states[2].g());

			pop_ratio20_sum = H2_X_colden[0][2]/SDIV(H2_X_colden[0][0]);
			T20_sum = -energyK/log(SDIV(pop_ratio20_sum) * states[0].g()/states[2].g());

			energyK = states[3].energy().K() - states[0].energy().K();
			pop_ratio30 = pop3/SDIV(pop0);
			T30 = -energyK/log(SDIV(pop_ratio30) * states[0].g()/states[3].g());

			pop_ratio30_sum = H2_X_colden[0][3]/SDIV(H2_X_colden[0][0]);
			T30_sum = -energyK/log(SDIV(pop_ratio30_sum) * states[0].g()/states[3].g());

			energyK = states[3].energy().K() - states[1].energy().K();
			pop_ratio31 = pop3/SDIV(pop1);
			T31 = -energyK/log(SDIV(pop_ratio31) * states[1].g()/states[3].g());

			pop_ratio31_sum = H2_X_colden[0][3]/SDIV(H2_X_colden[0][1]);
			T31_sum = -energyK/log(SDIV(pop_ratio31_sum) * states[1].g()/states[3].g());

			energyK = states[4].energy().K() - states[0].energy().K();
			pop_ratio40 = pop4/SDIV(pop0);
			T40 = -energyK/log(SDIV(pop_ratio40) * states[0].g()/states[4].g());

			pop_ratio40_sum = H2_X_colden[0][4]/SDIV(H2_X_colden[0][0]);
			T40_sum = -energyK/log(SDIV(pop_ratio40_sum) * states[0].g()/states[4].g());
		}
		else
		{
			pop_ratio10 = 0.;
			T10 = 0.;
			T20 = 0.;
			T30 = 0.;
			T31 = 0.;
			T40 = 0.;
			T10_sum = 0.;
			T20_sum = 0.;
			T30_sum = 0.;
			T31_sum = 0.;
			T40_sum = 0.;
		}

		/* various temperatures for neutral/molecular gas */
		fprintf( io, 
			"%.5e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n" ,
			/* depth in cm */
			radius.depth_mid_zone ,
			/* total H2 fraction */
			(*dense_total)/dense.gas_phase[ipHYDROGEN] ,
			/* ratio of populations of 1 to 0 only */
			pop_ratio10 ,
			/* sum of all ortho and para */
			ortho_density / SDIV(para_density),
			T10,T20,T30,T31,T40,
			phycon.te ,
			hyperfine.Tspin21cm,T10_sum,T20_sum,T30_sum,T31_sum,T40_sum  );
	}
	else if( (strcmp(chJOB , "H2ln" ) == 0) && (strcmp(chTime,"LAST") == 0) )
	{
		/* save H2 lines - output the full emission-line spectrum */
		double thresh;
		double renorm;
		/* first test, is H2 turned on?  Second test, have lines arrays
		 * been set up - nsum is negative if abort occurs before lines
		 * are set up */
		if( lgEnabled && LineSave.nsum > 0)
		{
			ASSERT( LineSave.ipNormWavL >= 0 );
			/* get the normalization line */
			if( LineSave.lines[LineSave.ipNormWavL].SumLine(0) > SMALLFLOAT )
				renorm = LineSave.ScaleNormLine/
				LineSave.lines[LineSave.ipNormWavL].SumLine(0);
			else
				renorm = 1.;

			if( renorm > SMALLFLOAT )
			{
				/* this is threshold for faintest line, normally 0, set with 
				 * number on save H2 command */
				thresh = thresh_punline_h2/(realnum)renorm;
			}
			else
				thresh = 0.f;
						
			/* save H2 line intensities at end of iteration 
			 * nElecLevelOutput is electronic level with 1 for ground, so this loop is < nElecLevelOutput */
			for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
			{
				qList::iterator Hi = ( (*tr).Hi() );	
				qList::iterator Lo = ( (*tr).Lo() );
				long iElecHi = (*Hi).n();
				long iVibHi = (*Hi).v();
				long iRotHi = (*Hi).J();
				long iElecLo = (*Lo).n();
				long iVibLo = (*Lo).v();
				long iRotLo = (*Lo).J();	
				if( iElecHi >= nElecLevelOutput )
					continue;
				if( H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] > thresh )
				{
					/* air wavelength in microns */
					/* WLAng contains correction for index of refraction of air */
					double wl = (*tr).WLAng()/1e4;
					fprintf(io, "%li-%li %c(%li)", iVibHi, iVibLo, chMolBranch( iRotHi, iRotLo ), iRotLo );
					fprintf( io, "\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld", iElecHi , iVibHi , iRotHi , iElecLo , iVibLo , iRotLo);
					/* WLAng contains correction for index of refraction of air */
					fprintf( io, "\t%.7f\t", wl );
					/*prt_wl print floating wavelength in Angstroms, in output format */
					prt_wl( io , (*tr).WLAng() );
					/* the log of the line intensity or luminosity */
					fprintf( io, "\t%.3f\t%.3e", 
						log10( MAX2(1e-37, H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo]*radius.Conv2PrtInten) ), 
						H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo]*renorm );
					/* excitation energy of upper level in K */
					fprintf( io, "\t%.3f", (*Hi).energy().K() );
					/* the product g_hi h nu * Aul */
					fprintf( io, "\t%.3e", (*tr).Emis().Aul() * (*tr).EnergyErg() * (*(*tr).Hi()).g() );
					fprintf( io, "\n");
				}
			}
		}
	}
	else if( (strcmp(chJOB , "H2sp" ) == 0)  )
	{
		fprintf( io, "PUT SOMETHING HERE!\n" );
	}
	return;
}
 /*cdH2_Line determines intensity and luminosity of and H2 line.  The first
  * six arguments give the upper and lower quantum designation of the levels.
  * The function returns 1 if we found the line, 
  * and false==0 if we did not find the line because ohoto-para transition
  * or upper level has lower energy than lower level  */
long int cdH2_Line(
	  /* indices for the upper level */
	  long int iElecHi, 
	  long int iVibHi ,
	  long int iRotHi ,
	  /* indices for lower level */
	  long int iElecLo, 
	  long int iVibLo ,
	  long int iRotLo ,
	  /* linear intensity relative to normalization line*/
	  double *relint, 
	  /* luminosity or intensity of line */
	  double *absint )
{
	DEBUG_ENTRY( "cdH2_Line()" );
	return h2.getLine( iElecHi, iVibHi, iRotHi, iElecLo, iVibLo, iRotLo, relint, absint );
}

long int diatomics::getLine( long iElecHi, long iVibHi, long iRotHi, long iElecLo, long iVibLo, long iRotLo, double *relint, double *absint )
{

	DEBUG_ENTRY( "diatomics::getline()" );

	/* these will be return values if we can't find the line */
	*relint = 0.;
	*absint = 0.;

	/* for now both electronic levels must be zero */
	if( iElecHi!=0 || iElecLo!=0 )
	{
		return 0;
	}
	
	long ipHi = ipEnergySort[iElecHi][iVibHi][iRotHi];
	long ipLo = ipEnergySort[iElecLo][iVibLo][iRotLo];

	/* check that energy of first level is higher than energy of second level */
	if( states[ipHi].energy().WN() < states[ipLo].energy().WN() )
	{
		return 0;
	}

	/* check that ortho-para does not change */
	if( H2_lgOrtho[iElecHi][iVibHi][iRotHi] != H2_lgOrtho[iElecLo][iVibLo][iRotLo] )
	{
		return 0;
	}

	/* exit if lines does not exist */
	if( !lgH2_radiative[ipHi][ipLo] )
	{
		return 0;
	}

	ASSERT( LineSave.ipNormWavL >= 0 );
	/* does the normalization line have a positive intensity*/
	if( LineSave.lines[LineSave.ipNormWavL].SumLine(0) > 0. )
	{
		*relint = H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo]/
			LineSave.lines[LineSave.ipNormWavL].SumLine(0) * LineSave.ScaleNormLine;
	}
	else
	{
		*relint = 0.;
	}

	/* return log of line intensity if it is positive */
	if( H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] > 0. )
	{
		*absint = H2_SaveLine[iElecHi][iVibHi][iRotHi][iElecLo][iVibLo][iRotLo] * 
			radius.Conv2PrtInten;
	}
	else
	{
		*absint = 0.;
	}
	/* this indicates success */
	return 1;
}

void diatomics::set_numLevelsMatrix( long numLevels )
{
	if( !lgREAD_DATA ) 
		nXLevelsMatrix = numLevels;
}

/* Read LTE cooling per molecule */
void diatomics::H2_Read_LTE_cooling_per_H2()
{
	const string cdDATAFILE = "lte_cooling.dat";

	DEBUG_ENTRY( "H2_Read_LTE_cooling_per_H2()" );

	/* now open the data file */
	string chPath = path + cpu.i().chDirSeparator() + cdDATAFILE;
	FILE *ioDATA = open_data( chPath, "r" );

	LTE_Temp.resize( 0 );
	LTE_cool.resize( 0 );

	int nlines = 0;
	string chLine;
	while( read_whole_line( chLine, ioDATA ) )
	{
		if( chLine.length() == 0 || chLine[0] == '\n' || chLine[0] == '\r' || chLine[0] == ' ' )
			break;
		/* skip comment */
		if( chLine[0] == '#' )
			continue;

		double temp, cool;
		if( sscanf(chLine.c_str(), "%lf\t%le", &temp, &cool) != 2 )
			TotalInsanity();

		LTE_Temp.push_back( temp );
		LTE_cool.push_back( cool );

		nlines++;
	}
	fclose( ioDATA );

	{
		enum { DEBUG_LOC = false };
		if( DEBUG_LOC )
		{
			for( int i = 0; i < nlines; i++ )
			{
				fprintf(stdout, "%i\t %f\t %e\n",
					i, log10( LTE_Temp[i] ), LTE_cool[i] );
			}
			cdEXIT( EXIT_SUCCESS );
		}
	}

	return;
}

