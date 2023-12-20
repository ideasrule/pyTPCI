/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "taulines.h"
#include "trace.h"
#include "thirdparty.h"
#include "atmdat.h"
#include "lines_service.h"
#include "parser.h"
#include "parse_species.h"
#include "rfield.h"

typedef vector< pair<double,long> > DoubleLongPairVector;

struct LevelInfo
{
	double nrg;
	long index;
	double stwt;
	string config;
	LevelInfo(double e, long i, double w, string c) : nrg(e), index(i), stwt(w), config(c) {}
	bool operator< ( const LevelInfo& l ) const
	{
		if( nrg < l.nrg )
			return true;
		else if( nrg == l.nrg )
			return ( index < l.index );
		else
			return false;
	}
};

static const bool DEBUGSTATE = false;
// minimum energy difference (wavenumbers) we will accept
const double ENERGY_MIN_WN = 1e-10;

/* convert transition type into integer code */
inline int getCode(const string& transType)
{
	DEBUG_ENTRY( "getCode()" );

	int val = -1;
	if( transType.length() != 2 )
		return -1;
	if( transType == "NS" )
		return 0; // when not set, treat the transition as E1 (this is what NIST does)
	if( transType[0] == 'E' )
		val = 0;
	else if( transType[0] == 'M' )
		val = 3;
	else
		return -1;
	if( transType[1] == '1' )
		(void)0;
	else if( transType[1] == '2' )
		val += 1;
	else if( transType[1] == '3' )
		val += 2;
	else
		return -1;
	return val;
}

inline void processIndices(long ipLoInFile, long ipHiInFile, bool lgIsRegular, const map<long,long>& indexold2new,
			   long& ipLo, long& ipHi)
{
	DEBUG_ENTRY( "processIndices()" );

	if( lgIsRegular )
	{
		ipLo = ipLoInFile-1;
		ipHi = ipHiInFile-1;
	}
	else
	{
		// ignore transitions with an index that is not present in the level file
		// this is not an error as the eof marker in the level file may have been moved upward
		auto plo = indexold2new.find(ipLoInFile);
		if( plo == indexold2new.end() )
		{
			ipLo = ipHi = LONG_MAX;
			return;
		}
		auto phi = indexold2new.find(ipHiInFile);
		if( phi == indexold2new.end() )
		{
			ipLo = ipHi = LONG_MAX;
			return;
		}

		/* Account for reordered energy levels */
		ipLo = plo->second;
		ipHi = phi->second;
	}

	ASSERT( ipLo != ipHi );

	// swap indices if energy levels were not correctly sorted
	if( ipHi < ipLo )
	{
		long swap = ipHi;
		ipHi = ipLo;
		ipLo = swap;
	}
}

void atmdat_STOUT_readin( long intNS, const string& chPrefix )
{
	DEBUG_ENTRY( "atmdat_STOUT_readin()" );

	static const int MAX_NUM_LEVELS = 999;

	// the magic numbers for the stout files
	const long int iyr = 17, imo = 9, idy = 5;

	dBaseSpecies[intNS].lgMolecular = false;
	dBaseSpecies[intNS].lgLTE = false;

	setProperties(dBaseSpecies[intNS]);

	string chUnCaps = chPrefix;
	uncaps(chUnCaps);

	if( dBaseSpecies[intNS].dataset.length() > 0 )
		chUnCaps += "_" + dBaseSpecies[intNS].dataset;

	string chNRGFilename = chUnCaps + ".nrg";
	string chTPFilename = chUnCaps + ".tp";
	string chCOLLFilename = chUnCaps + ".coll";

	DataParser d;

	/******************************************************
	 ***************** Energy Levels File *****************
	 ******************************************************/
	d.open( chNRGFilename, ES_STARS_ONLY );

	long index = 0;
	double nrg = 0.0;
	double istwt, stwt = 0.0;

	/* first line is a version number - now confirm that it is valid */
	d.getline();
	d.checkMagic( iyr, imo, idy );

	/* Create array for holding energies and statistical weights so that
	 * we can put the energies in the correct order before moving them to
	 * dBaseStates */
	vector<LevelInfo> dBaseStatesOrg;

	// check for an end of file sentinel
	bool lgSentinelReached = false;
	while( d.getline() )
	{
		if( d.lgEODMarker() )
		{
			lgSentinelReached = true;
			break;
		}

		d.getToken( index );
		d.getToken( nrg );
		d.getToken( stwt );
		if( stwt <= 0. || modf(stwt, &istwt) != 0. )
			d.errorAbort( "invalid statistical weight" );
		string config;
		d.getQuoteOptional( config );

		d.checkEOL();

		dBaseStatesOrg.emplace_back(nrg, index, stwt, config);
	}

	if( !lgSentinelReached )
	{
		fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n", chNRGFilename.c_str() );
		fprintf( ioQQQ, " Check that stars (***) appear after the last line of data"
			 " and start in the first column of that line.\n");
		cdEXIT(EXIT_FAILURE);
	}

	long nMolLevs = dBaseStatesOrg.size();
	long HighestIndexInFile = nMolLevs;

	dBaseSpecies[intNS].numLevels_max = HighestIndexInFile;

	if( tolower(dBaseSpecies[intNS].chLabel[0]) == 'f' && tolower(dBaseSpecies[intNS].chLabel[1]) == 'e' )
	{
		// Fe is special case with more levels
		nMolLevs = MIN3(nMolLevs, atmdat.nStoutMaxLevelsFe, MAX_NUM_LEVELS );
	}
	else
	{
		nMolLevs = MIN3(nMolLevs, atmdat.nStoutMaxLevels, MAX_NUM_LEVELS );
	}

	//Consider the masterlist specified number of levels as the min. =1 if not specified
	long numMasterlist = MIN2( dBaseSpecies[intNS].numLevels_masterlist, HighestIndexInFile );
	nMolLevs = MAX2(nMolLevs,numMasterlist);

	if (dBaseSpecies[intNS].setLevels != -1)
	{
		if (dBaseSpecies[intNS].setLevels > HighestIndexInFile)
		{
			char chLabelChemical[CHARS_SPECIES] = "";
			spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel );
			if( dBaseSpecies[intNS].setLevels != LONG_MAX )
				fprintf( ioQQQ, "Using STOUT spectrum %s (species: %s) with %li requested,"
					 " only %li energy levels available.\n",
					 dBaseSpecies[intNS].chLabel, chLabelChemical, dBaseSpecies[intNS].setLevels,
					 HighestIndexInFile );
			nMolLevs = HighestIndexInFile;	  
		}
		else
		{
			nMolLevs = dBaseSpecies[intNS].setLevels;
		}
	}

	ASSERT( nMolLevs <= HighestIndexInFile );

	if( atmdat.lgStoutPrint )
	{
		char chLabelChemical[CHARS_SPECIES] = "";
		spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
		fprintf( ioQQQ, "Using STOUT spectrum %s (species: %s) with %li levels of %li available.\n",
				 dBaseSpecies[intNS].chLabel, chLabelChemical, nMolLevs, HighestIndexInFile );
	}

	dBaseSpecies[intNS].numLevels_max = nMolLevs;
	dBaseSpecies[intNS].numLevels_local = dBaseSpecies[intNS].numLevels_max;

	/*Resize the States array*/
	dBaseStates[intNS].init(dBaseSpecies[intNS].chLabel,nMolLevs);
	/*Zero out the maximum wavenumber for each species */
	dBaseSpecies[intNS].maxWN = 0.;

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\nStout Species: %s\n",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Energy Level File: %s\n",chNRGFilename.c_str());
		fprintf(ioQQQ,"Number of Energy Levels in File: %li\n",HighestIndexInFile);
		fprintf(ioQQQ,"Number of Energy Levels Cloudy is Currently Using: %li\n",nMolLevs);
		fprintf(ioQQQ,"Species|File Index|Energy(WN)|StatWT\n");
		for( size_t i=0; i < dBaseStatesOrg.size(); ++i )
		{
			fprintf( ioQQQ, "<%s>\t%li\t%.3f\t%.1f\n", dBaseSpecies[intNS].chLabel,
				 dBaseStatesOrg[i].index, dBaseStatesOrg[i].nrg, dBaseStatesOrg[i].stwt );
		}
	}

	//Sort levels by energy (then by the index in the file if necessary)
	sort(dBaseStatesOrg.begin(),dBaseStatesOrg.end());

	// regular files have the level indices running as 1, 2, 3, ...
	// (without gaps) and the energies are also already sorted
	bool lgIsRegular = true;
	for( size_t i=0; i < dBaseStatesOrg.size(); i++ )
	{
		if( dBaseStatesOrg[i].index != long(i+1) )
		{
			lgIsRegular = false;
			break;
		}
	}

	map<long, long> indexold2new;
	if( !lgIsRegular ) {
		for( size_t i=0; i < dBaseStatesOrg.size(); i++ )
		{
			long old = dBaseStatesOrg[i].index;
			auto p = indexold2new.find(old);
			if( p == indexold2new.end() )
				indexold2new[old] = i;
			else
			{
				fprintf(ioQQQ, "Duplicate index %ld found in file %s\n", old, chNRGFilename.c_str());
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\n\n***Energy levels have been sorted in order of increasing energy.***\n");
		fprintf(ioQQQ,"Species|File Index|Sorted Index|Energy(WN)|StatWT\n");
	}

	/* Store sorted energies in their permanent home */
	for( long i=0; i < nMolLevs; i++ )
	{
		double nrg = dBaseStatesOrg[i].nrg;
		long oldindex = dBaseStatesOrg[i].index;
		double stwt = dBaseStatesOrg[i].stwt;

		if( DEBUGSTATE )
			fprintf( ioQQQ, "<%s>\t%li\t%li\t%.3f\t%.1f\n", dBaseSpecies[intNS].chLabel,
				 oldindex, i+1, nrg, stwt );

		dBaseStates[intNS][i].energy().set(nrg,"cm^-1");
		dBaseStates[intNS][i].g() = stwt;
	}

	/* allocate the Transition array*/
	ipdBaseTrans[intNS].reserve(nMolLevs);
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
		ipdBaseTrans[intNS].reserve(ipHi,ipHi);
	ipdBaseTrans[intNS].alloc();
	dBaseTrans[intNS].resize(ipdBaseTrans[intNS].size());
	dBaseTrans[intNS].states() = &dBaseStates[intNS];
	dBaseTrans[intNS].chLabel() = dBaseSpecies[intNS].chLabel;
	dBaseSpecies[intNS].database = "Stout";

	//This is creating transitions that we don't have collision data for
	//Maybe use gbar or keep all of the Fe 2 even if it was assumed (1e-5)
	int ndBase = 0;
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
	{
		for( long ipLo = 0; ipLo < ipHi; ipLo++)
		{
			ipdBaseTrans[intNS][ipHi][ipLo] = ndBase;
			dBaseTrans[intNS][ndBase].Junk();
			dBaseTrans[intNS][ndBase].setLo(ipLo);
			dBaseTrans[intNS][ndBase].setHi(ipHi);
			++ndBase;
		}
	}

	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	for( auto tr=dBaseTrans[intNS].begin(); tr != dBaseTrans[intNS].end(); ++tr )
	{
		int ipHi = tr->ipHi();
		int ipLo = tr->ipLo();
		ASSERT( ipHi > ipLo &&
			dBaseStates[intNS][ipHi].energy().WN() >= dBaseStates[intNS][ipLo].energy().WN() );
		double fenergyWN = dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN();
		tr->EnergyWN() = fenergyWN;
		if( rfield.isEnergyBound( Energy( fenergyWN, "cm^-1" ) ) )
		{
			tr->WLAng() = (realnum) wn2ang( fenergyWN );
			dBaseSpecies[intNS].maxWN = MAX2(dBaseSpecies[intNS].maxWN,tr->EnergyWN());
		}
		else
			tr->WLAng() = 1e30;
	}

	/******************************************************
	 ************* Transition Probability File ************
	 ******************************************************/
	d.open( chTPFilename, ES_STARS_ONLY );

	/* first line is a version number - now confirm that it is valid */
	d.getline();
	d.checkMagic( iyr, imo, idy );

	/* lgLineStrengthTT functions as a checklist for line strengths
	 * from the various transition types (E1,M2, etc.)
	 * When a line strength is added to the total Aul from a specific
	 * transition type, the corresponding value (based on ipLo, ipHi, and transition type)
	 * is set to true.
	 * lgLineStrengthTT[ipLo][ipHi][k]:
	 * k = 0 => E1
	 * k = 1 => E2
	 * k = 2 => E3
	 * k = 3 => M1
	 * k = 4 => M2
	 * k = 5 => M3
	 */
	static const int intNumCols = 6;
	multi_arr<bool,3,C_TYPE> lgLineStrengthTT(nMolLevs, nMolLevs, intNumCols);
	lgLineStrengthTT = false;

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\nStout Species: %s\n",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Radiative Data File: %s\n",chTPFilename.c_str());
		fprintf(ioQQQ,"Species|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Data Type (A,G,S)|Data\n");
	}

	//Read the remaining lines of the transition probability file
	lgSentinelReached = false;
	while( d.getline() )
	{
		if( d.lgEODMarker() )
		{
			lgSentinelReached = true;
			break;
		}

		string dataType;
		d.getToken( dataType );

		if( dataType != "A" && dataType != "G" && dataType != "S" )
			d.errorAbort( "invalid data type" );

		long ipLoInFile, ipHiInFile, ipLo, ipHi;
		d.getToken( ipLoInFile );
		d.getToken( ipHiInFile );
		processIndices(ipLoInFile, ipHiInFile, lgIsRegular, indexold2new, ipLo, ipHi);

		if( ipHi >= nMolLevs )
			continue;

		double tpData;
		d.getToken( tpData );
		if( tpData < 0. )
			d.errorAbort( "transition probability data must be >= 0" );

		if( tpData < atmdat.aulThreshold && dataType == "A" )
		{
			// skip these lines
			continue;
		}

		string transType;
		if( !d.getTokenOptional( transType ) )
			transType = "NS"; // not set, will be treated as E1
		// sometimes NIST sets the transition type to UT (undefined transition type)
		// this is used for forbidden transitions where they cannot decide if they
		// are M1 or E2, so we treat them here as if they are both...
		if( transType == "UT" )
			transType = "M1+E2";
		if( dataType == "S" )
		{
			if( transType == "NS" )
				d.errorAbort( "line strength specified, but no transition type was found" );
			if( transType.length() != 2 )
				d.errorAbort( "invalid transition type" );
		}

		size_t len = transType.length();
		size_t p = 0;
		// correctly handle transition types like "M1+E2"
		while( p < len )
		{
			string tt = transType.substr(p,2);
			int ttCode = getCode(tt);
			char c = ( p+2 >= len ) ? '+' : transType[p+2];
			if( ttCode < 0 || c != '+' )
				d.errorAbort( "invalid transition type" );
			if( lgLineStrengthTT[ipLo][ipHi][ttCode] )
				d.errorAbort( "this transition already has an Aul value set" );
			lgLineStrengthTT[ipLo][ipHi][ttCode]= true;
			p += 3;
		}

		d.checkEOL();

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"<%s>\t%li:%li\t%li:%li\t%s\t%.2e\n",
				dBaseSpecies[intNS].chLabel,ipLoInFile,ipHiInFile,
				ipLo+1,ipHi+1,dataType.c_str(),tpData);
		}

		TransitionList::iterator tr = dBaseTrans[intNS].begin()+ipdBaseTrans[intNS][ipHi][ipLo];

		/* If we don't already have the transition in the stack, add it and
		 * zero out Aul */
		if( !tr->hasEmis() )
		{
			tr->AddLine2Stack();
			tr->Emis().Aul() = 0.;
			tr->Emis().gf() = 0.;
		}

		//This means last data column has Aul.
		if( dataType == "A" )
		{
			if( tr->EnergyWN() > ENERGY_MIN_WN )
			{
				tr->Emis().Aul() += tpData;
				// use updated total Aul to get gf
				tr->Emis().gf() = (realnum)GetGF(tr->Emis().Aul(), tr->EnergyWN(), tr->Hi()->g());
			}
		}
		//This means last data column has gf.
		else if( dataType == "G" )
		{
			if( tr->EnergyWN() > ENERGY_MIN_WN )
			{
				tr->Emis().gf() += tpData;
				// use updated total gf to get Aul
				tr->Emis().Aul() = (realnum)eina(tr->Emis().gf(), tr->EnergyWN(), tr->Hi()->g());
			}
		}
		//This means last data column has line strengths.
		else if( dataType == "S" )
		{
			if( tr->EnergyWN() > ENERGY_MIN_WN )
			{
				tr->Emis().Aul() += S2Aul(tpData, tr->EnergyAng(), tr->Hi()->g(), transType);
				tr->Emis().gf() = (realnum)GetGF(tr->Emis().Aul(), tr->EnergyWN(), tr->Hi()->g());
			}
		}

		tr->setComment( db_comment_tran_levels( dBaseStatesOrg[ipLo].config, dBaseStatesOrg[ipHi].config ) );
	}

	if( !lgSentinelReached )
	{
		fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n", chTPFilename.c_str() );
		fprintf( ioQQQ, " Check that stars (*****) appear after the last line of data"
			 " and start in the first column of that line.");
		cdEXIT(EXIT_FAILURE);
	}

	/******************************************************
	 ************* Collision Data File ********************
	 ******************************************************/
	d.open( chCOLLFilename, ES_STARS_ONLY );

	/* first line is a version number - now confirm that it is valid */
	d.getline();
	d.checkMagic( iyr, imo, idy );

	/****** Could add ability to count number of temperature changes, electron CS, and proton CS ****/

	/* allocate space for collision strengths */
	StoutCollData[intNS].alloc(nMolLevs,nMolLevs,ipNCOLLIDER);
	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			for( long k=0; k<ipNCOLLIDER; k++ )
			{
				/* initialize all spline variables */
				StoutCollData[intNS].junk(ipHi,ipLo,k);
			}
		}
	}

	int numpoints = 0;
	vector<double> temps;
	long ipCollider = -1;

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\nStout Species: %s\n",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Collision Data File: %s\n",chCOLLFilename.c_str());
		fprintf(ioQQQ,"Species|TEMP|Temperatures (K)\n");
		fprintf(ioQQQ,"Species|Data Type (CS,RATE)|Collider|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Data\n");
	}

	//Read the remaining lines of the collision data file
	lgSentinelReached = false;
	while( d.getline() )
	{
		/* Stop on *** */
		if( d.lgEODMarker() )
		{
			lgSentinelReached = true;
			break;
		}

		string dataType;
		d.getToken( dataType );

		//This is a temperature line
		if( dataType == "TEMP" )
		{
			if( DEBUGSTATE )
				fprintf(ioQQQ,"<%s>\tTEMP\t",dBaseSpecies[intNS].chLabel);

			temps.clear();
			double data;
			while( d.getTokenOptional(data) )
			{
				if( data <= 0. )
					d.errorAbort( "invalid temperature" );
				if( temps.size() > 0 && data <= temps.back() )
					d.errorAbort( "temperatures must be monotonically increasing" );
				temps.emplace_back( data );
				if( DEBUGSTATE )
					fprintf(ioQQQ,"%.2e\t",data);
			}
			if( DEBUGSTATE )
				fprintf(ioQQQ,"\n");
			numpoints = temps.size();
			ASSERT( numpoints > 0 );
		}
		else if( dataType == "CS" || dataType == "RATE" )
		{
			bool isRate = ( dataType == "RATE" );
			string colliderType;
			d.getToken( colliderType );

			if( colliderType == "ELECTRON" )
				ipCollider = ipELECTRON;
			else if( colliderType == "PROTON" || colliderType == "H+" )
				ipCollider = ipPROTON;
			else if( colliderType == "HE+2" )
				ipCollider = ipALPHA;
			else if( colliderType == "HE+" )
				ipCollider = ipHE_PLUS;
			else if( colliderType == "H2ORTHO" )
				ipCollider = ipH2_ORTHO;
			else if( colliderType == "H2PARA" )
				ipCollider = ipH2_PARA;
			else if( colliderType == "H2" )
				ipCollider = ipH2;
			else if( colliderType == "HE" )
				ipCollider = ipATOM_HE;
			else if( colliderType == "H" )
				ipCollider = ipATOM_H;
			else
				d.errorAbort( "invalid type of collider for RATE, I know about"
					"ELECTRON PROTON H+ HE+2 HE+ H2ORTHO H2PARA H2 HE " );

			if( !isRate && ipCollider != ipELECTRON )
				d.errorAbort( "collision strengths are only allowed for electron colliders" );

			if( temps.empty() )
				d.errorAbort( "must specify temperatures before the collision strengths" );

			long ipLoInFile, ipHiInFile, ipLo, ipHi;
			d.getToken( ipLoInFile );
			d.getToken( ipHiInFile );
			processIndices(ipLoInFile, ipHiInFile, lgIsRegular, indexold2new, ipLo, ipHi);

			if( ipHi >= nMolLevs )
				continue;

			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"<%s>\t%s\t%li\t%li:%li\t%li:%li",
					dBaseSpecies[intNS].chLabel,isRate?"RATE":"CS",ipCollider,
					ipLoInFile,ipHiInFile,ipLo+1,ipHi+1);
			}

			/* Set this as a collision strength not a collision rate coefficient*/
			StoutCollData[intNS].lgIsRate(ipHi,ipLo,ipCollider) = isRate;

			ASSERT( numpoints > 0 );
			if( StoutCollData[intNS].ntemps(ipHi,ipLo,ipCollider) > 0 )
				d.errorAbort( "duplicate collisional transition found" );
			StoutCollData[intNS].setpoints(ipHi,ipLo,ipCollider,numpoints);

			for( int j = 0; j < numpoints; j++ )
			{
				StoutCollData[intNS].temps(ipHi,ipLo,ipCollider)[j] = temps[j];
				d.getToken( StoutCollData[intNS].collstrs(ipHi,ipLo,ipCollider)[j] );
				if( StoutCollData[intNS].collstrs(ipHi,ipLo,ipCollider)[j] <= 0. )
					d.errorAbort( "invalid collisional data" );
				if( DEBUGSTATE )
					fprintf(ioQQQ,"\t%.2e",StoutCollData[intNS].collstrs(ipHi,ipLo,ipCollider)[j]);
			}
			if( DEBUGSTATE )
				fprintf(ioQQQ,"\n");
		}
		else
		{
			d.errorAbort( "invalid data type" );
		}

		d.checkEOL();
	}

	if( !lgSentinelReached )
	{
		fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n", chCOLLFilename.c_str() );
		fprintf( ioQQQ, " Check that stars (*****) appear after the last line of data"
			 " and start in the first column.");
		cdEXIT(EXIT_FAILURE);
	}
}

void atmdat_CHIANTI_readin( long intNS, const string& chPrefix )
{
	DEBUG_ENTRY( "atmdat_CHIANTI_readin()" );

	int intsplinepts,intTranType,intxs;
	long int nMolLevs,nMolExpLevs,nElvlcLines,nTheoLevs;// number of experimental and total levels
	FILE *ioElecCollData=NULL, *ioProtCollData=NULL;
	realnum  fstatwt,fenergyWN,fWLAng,fenergy,feinsteina;
	double fScalingParam,fEnergyDiff;
	const char chCommentChianti = '#';

	bool lgProtonData=false;

	// this is the largest number of levels allowed by the chianti format, I3
	static const int MAX_NUM_LEVELS = 999;

	dBaseSpecies[intNS].lgMolecular = false;
	dBaseSpecies[intNS].lgLTE = false;

	string chUnCaps = chPrefix;
	uncaps(chUnCaps);
	
	string chEnFilename = chUnCaps;
	string chTraFilename = chUnCaps;
	string chEleColFilename = chUnCaps;
	string chProColFilename = chUnCaps;

	/*For the CHIANTI DATABASE*/
	/*Open the energy levels file*/
	chEnFilename += ".elvlc";

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chEnFilename.c_str());

	fstream elvlcstream;
	open_data( elvlcstream, chEnFilename, mode_r );

	/*Open the transition probabilities file*/
	chTraFilename += ".wgfa";

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chTraFilename.c_str());

	fstream wgfastream;
	open_data( wgfastream, chTraFilename, mode_r );

	/*Open the electron collision data*/
	chEleColFilename += ".splups";

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chEleColFilename.c_str());

	ioElecCollData = open_data( chEleColFilename, "r" );

	/*Open the proton collision data*/
	chProColFilename += ".psplups";

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chProColFilename.c_str());

	/*We will set a flag here to indicate if the proton collision strengths are available */
	if( ( ioProtCollData = open_data( chProColFilename, "r", AS_TRY ) ) != NULL )
	{
		lgProtonData = true;
	}
	else
	{
		lgProtonData = false;
	}

	/*Loop finds how many theoretical and experimental levels are in the elvlc file */
	//eof_col is used get the first 4 charcters per line to find end of file
	const int eof_col = 5;
	//length (+1) of the nrg in the elvlc file
	const int lvl_nrg_col=16;
	//# of columns skipped from the left to get to nrg start
	const int lvl_skipto_nrg = 40;
	/* # of columns to skip from eof check to nrg start */
	const int lvl_eof_to_nrg = lvl_skipto_nrg - eof_col + 1;
	//# of columns to skip over the rydberg energy, we don't use it
	const int lvl_skip_ryd = 15;
	nElvlcLines = 0;
	nMolExpLevs = 1;
	nTheoLevs = 1;
	if (elvlcstream.is_open())
	{
		int nj = 0;
		char otemp[eof_col];
		char exptemp[lvl_nrg_col],theotemp[lvl_nrg_col];
		double tempexpenergy = 0.,theoenergy = 0.;
		/*This loop counts the number of valid rows within the elvlc file
		  as well as the number of experimental energy levels.*/
		while(nj != -1)
		{
			elvlcstream.get(otemp,eof_col);
			nj = atoi(otemp);
			if( nj == -1)
				break;
			nElvlcLines++;

			elvlcstream.seekg(lvl_eof_to_nrg,ios::cur);
			elvlcstream.get(exptemp,lvl_nrg_col);
			tempexpenergy = (double) atof(exptemp);
			if( tempexpenergy != 0.)
				nMolExpLevs++;

			elvlcstream.seekg(lvl_skip_ryd,ios::cur);
			elvlcstream.get(theotemp,lvl_nrg_col);
			theoenergy = (double) atof(theotemp);
			if( theoenergy != 0. )
				nTheoLevs++;

			elvlcstream.ignore(INT_MAX,'\n');

		}
		elvlcstream.seekg(0,ios::beg);
	}

	//Sometimes the theoretical chianti level data is incomplete.
	//If it is bad use experimental
	bool lgChiaBadTheo = false;
	if( !atmdat.lgChiantiExp && nTheoLevs < nElvlcLines )
	{
		lgChiaBadTheo = true;
		atmdat.lgChiantiExp = true;
		fprintf(ioQQQ,"Warning: The theoretical energy levels for %s are incomplete.",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Switching to the experimental levels for this species.");
	}

	long HighestIndexInFile = -1;

	/* The total number of levels depends on the experimental Chianti switch */
	if( atmdat.lgChiantiExp )
	{
		HighestIndexInFile = nMolExpLevs;
	}
	else
	{
		HighestIndexInFile = nElvlcLines;
	}

	dBaseSpecies[intNS].numLevels_max = HighestIndexInFile;

	setProperties(dBaseSpecies[intNS]);
	
	if( tolower(dBaseSpecies[intNS].chLabel[0]) == 'f' && tolower(dBaseSpecies[intNS].chLabel[1]) == 'e')
	{
		// Fe is special case with more levels
		nMolLevs = MIN3(HighestIndexInFile, atmdat.nChiantiMaxLevelsFe,MAX_NUM_LEVELS );
	}
	else
	{
		nMolLevs = MIN3(HighestIndexInFile, atmdat.nChiantiMaxLevels,MAX_NUM_LEVELS );
	}

	if( nMolLevs <= 0 )
	{
		fprintf( ioQQQ, "The number of energy levels is non-positive in datafile %s.\n", chEnFilename.c_str() );
		fprintf( ioQQQ, "The file must be corrupted.\n" );
		cdEXIT( EXIT_FAILURE );
	}

	//Consider the masterlist specified number of levels as the min. =1 if not specified
	long numMasterlist = MIN2( dBaseSpecies[intNS].numLevels_masterlist , HighestIndexInFile );
	nMolLevs = MAX2(nMolLevs,numMasterlist);

	if (dBaseSpecies[intNS].setLevels != -1)
	{
		if (dBaseSpecies[intNS].setLevels > HighestIndexInFile)
		{
			char chLabelChemical[CHARS_SPECIES] = "";
			spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel );
			if( dBaseSpecies[intNS].setLevels != LONG_MAX )
				fprintf( ioQQQ,"Using CHIANTI spectrum %s (species: %s) with %li requested,"
					 " only %li energy levels available.\n",
					 dBaseSpecies[intNS].chLabel, chLabelChemical, dBaseSpecies[intNS].setLevels,
					 HighestIndexInFile );
			nMolLevs = HighestIndexInFile;
		}
		else
		{
			nMolLevs = dBaseSpecies[intNS].setLevels;
		}
	}

	dBaseSpecies[intNS].numLevels_max = nMolLevs;
	dBaseSpecies[intNS].numLevels_local = dBaseSpecies[intNS].numLevels_max;

	if( atmdat.lgChiantiPrint == true)
	{
		if( atmdat.lgChiantiExp )
		{
			char chLabelChemical[CHARS_SPECIES] = "";
			spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
			fprintf( ioQQQ,"Using CHIANTI spectrum %s (species: %s) with %li experimental energy levels of %li available.\n",
				dBaseSpecies[intNS].chLabel, chLabelChemical, nMolLevs , nMolExpLevs );
		}
		else
		{
			char chLabelChemical[CHARS_SPECIES] = "";
			spectral_to_chemical( chLabelChemical, dBaseSpecies[intNS].chLabel ),
			fprintf( ioQQQ,"Using CHIANTI spectrum %s (species: %s) with %li theoretical energy levels of %li available.\n",
				dBaseSpecies[intNS].chLabel, chLabelChemical, nMolLevs , nElvlcLines );
		}
	}

	/* allocate the States array*/
	dBaseStates[intNS].init(dBaseSpecies[intNS].chLabel,nMolLevs);

	/* allocate the Transition array*/
	ipdBaseTrans[intNS].reserve(nMolLevs);
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
		ipdBaseTrans[intNS].reserve(ipHi,ipHi);
	ipdBaseTrans[intNS].alloc();
	dBaseTrans[intNS].resize(ipdBaseTrans[intNS].size());
	dBaseTrans[intNS].states() = &dBaseStates[intNS];
	dBaseTrans[intNS].chLabel() = dBaseSpecies[intNS].chLabel;
	dBaseSpecies[intNS].database = "Chianti";

	int ndBase = 0;
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
	{
		for( long ipLo = 0; ipLo < ipHi; ipLo++)
		{
			ipdBaseTrans[intNS][ipHi][ipLo] = ndBase;
			dBaseTrans[intNS][ndBase].Junk();
			dBaseTrans[intNS][ndBase].setLo(ipLo);
			dBaseTrans[intNS][ndBase].setHi(ipHi);
			++ndBase;
		}
	}

	/*Keep track of which levels have experimental data and then create a vector
	which relates their indices to the default chianti energy indices.
	 */
	long ncounter = 0;

	//Relate Chianti level indices to a set that only include experimental levels
	vector<long> intExperIndex(nElvlcLines,-1);

	DoubleLongPairVector dBaseStatesEnergy;
	vector<double> dBaseStatesStwt(HighestIndexInFile,-1.0);
	for( long ii = 0; ii < HighestIndexInFile; ii++ )
	{
		dBaseStatesEnergy.push_back(make_pair(-1.0,ii));
	}

	//lvl_skipto_statwt is the # of columns to skip to statwt from left
	const int lvl_skipto_statwt = 37;
	//lvl_statwt_col is the length (+1) of statwt
	const int lvl_statwt_col = 4;
	//Read in stat weight and energy

	//Read in nrg levels to see if they are in order
	for( long ipLev=0; ipLev<nElvlcLines; ipLev++ )
	{
		if(elvlcstream.is_open())
		{
			char gtemp[lvl_statwt_col],thtemp[lvl_nrg_col],obtemp[lvl_nrg_col];
			elvlcstream.seekg(lvl_skipto_statwt,ios::cur);
			elvlcstream.get(gtemp,lvl_statwt_col);
			fstatwt = (realnum)atof(gtemp);
			elvlcstream.get(thtemp,lvl_nrg_col);
			fenergy = (double) atof(thtemp);

			if(fstatwt <= 0.)
			{
				fprintf( ioQQQ, " WARNING: A positive non zero value is expected for the "
						"statistical weight but was not found in %s"
						" level %li\n", chEnFilename.c_str(),ipLev);
				cdEXIT(EXIT_FAILURE);
			}

			if( atmdat.lgChiantiExp )
			{
				/* Go through the entire level list selectively choosing only experimental level energies.
				 * Store them, not zeroes, in order using ncounter to count the index.
				 * Any row on the level list where there is no experimental energy, put a -1 in the relational vector.
				 * If it is a valid experimental energy level store the new ncounter index.
				 */

				if( fenergy != 0. || ipLev == 0 )
				{
					dBaseStatesEnergy.at(ncounter).first = fenergy;
					dBaseStatesEnergy.at(ncounter).second = ncounter;
					dBaseStatesStwt.at(ncounter) = fstatwt;
					intExperIndex.at(ipLev) = ncounter;
					ncounter++;
				}
				else
				{
					intExperIndex.at(ipLev) = -1;
				}
			}
			else
			{
				elvlcstream.seekg(lvl_skip_ryd,ios::cur);
				elvlcstream.get(obtemp,lvl_nrg_col);
				fenergy = (double) atof(obtemp);
				if(fenergy != 0. || ipLev == 0)
				{
					dBaseStatesEnergy.at(ipLev).first = fenergy;
					dBaseStatesEnergy.at(ipLev).second = ipLev;
					dBaseStatesStwt.at(ipLev) = fstatwt;
				}
				else
				{
					dBaseStatesEnergy.at(ipLev).first = -1.;
					dBaseStatesEnergy.at(ipLev).second = ipLev;
					dBaseStatesStwt.at(ipLev) = -1.;
				}
			}

			elvlcstream.ignore(INT_MAX,'\n');
		}
		else
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEnFilename.c_str());
			fclose( ioProtCollData );
			cdEXIT(EXIT_FAILURE);
		}
	}

	elvlcstream.close();

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\nintExperIndex Vector:\n");
		fprintf(ioQQQ,"File Index|Exper Index\n");
		for( vector<long>::iterator i = intExperIndex.begin(); i != intExperIndex.end(); i++ )
		{
			// term on rhs is long in 64 bit, int in 32 bit, print with long format
			long iPrt = (i-intExperIndex.begin())+1;
			fprintf(ioQQQ,"%li\t%li\n",iPrt,(*i)+1);
		}

		for( DoubleLongPairVector::iterator i=dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
		{
			// term on rhs is long in 64 bit, int in 32 bit, print with long format
			long iPrt = (i-dBaseStatesEnergy.begin())+1;
			fprintf(ioQQQ,"PreSort:%li\t%li\t%f\t%f\n",iPrt,
					(i->second)+1,i->first,dBaseStatesStwt.at(i->second));
		}
	}

	//Sort energy levels
	sort(dBaseStatesEnergy.begin(),dBaseStatesEnergy.end());

	std::vector<long> indexold2new(dBaseStatesEnergy.size());
	for( DoubleLongPairVector::const_iterator i = dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
	{
		indexold2new[i->second] = i-dBaseStatesEnergy.begin();
	}

	if( DEBUGSTATE )
	{
		for( DoubleLongPairVector::iterator i=dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
		{
			// term on rhs is long in 64 bit, int in 32 bit, print with long format
			long iPrt = (i-dBaseStatesEnergy.begin())+1;
			if( iPrt > nMolLevs )
				break;
			fprintf(ioQQQ,"PostSort:%li\t%li\t%f\t%f\n",iPrt,
					(i->second)+1,i->first,dBaseStatesStwt.at(i->second));
		}

		fprintf(ioQQQ,"\nChianti Species: %s\n",dBaseSpecies[intNS].chLabel);
		fprintf(ioQQQ,"Energy Level File: %s\n",chEnFilename.c_str());
		if( atmdat.lgChiantiExp )
		{
			fprintf(ioQQQ,"Number of Experimental Energy Levels in File: %li\n",nMolExpLevs);
		}
		else
		{
			fprintf(ioQQQ,"Number of Theoretical Energy Levels in File: %li\n",nElvlcLines);
		}
		fprintf(ioQQQ,"Number of Energy Levels Cloudy is Currently Using: %li\n",nMolLevs);
		fprintf(ioQQQ,"Species|File Index|Cloudy Index|StatWT|Energy(WN)\n");
	}

	vector<long> revIntExperIndex;
	if ( atmdat.lgChiantiExp )
	{
		revIntExperIndex.resize(dBaseStatesEnergy.size());
		for (size_t i = 0; i<dBaseStatesEnergy.size(); ++i)
			revIntExperIndex[i] = -1;
		for ( vector<long>::const_iterator i= intExperIndex.begin();
		      i != intExperIndex.end(); ++i )
		{
			long ipos = intExperIndex[i-intExperIndex.begin()];
			if (ipos >= 0 && ipos < long(dBaseStatesEnergy.size()))
			    revIntExperIndex[ipos] = i-intExperIndex.begin();
		}
	}
	
	for( DoubleLongPairVector::iterator i=dBaseStatesEnergy.begin(); i != dBaseStatesEnergy.end(); i++ )
	{

		long ipLevNew = i - dBaseStatesEnergy.begin();
		long ipLevFile = -1;

		if( ipLevNew >= nMolLevs )
			break;

		if( atmdat.lgChiantiExp )
		{
			ipLevFile = revIntExperIndex[ipLevNew];
		}
		else
		{
			ipLevFile = i->second;
		}

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"<%s>\t%li\t%li\t",dBaseSpecies[intNS].chLabel,ipLevFile+1,ipLevNew+1);
		}

		dBaseStates[intNS][ipLevNew].g() = dBaseStatesStwt.at(i->second);
		dBaseStates[intNS][ipLevNew].energy().set(i->first,"cm^-1");

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"%.1f\t",dBaseStatesStwt.at(i->second));
			fprintf(ioQQQ,"%.3f\n",i->first);
		}
	}

	// highest energy transition in chianti
	dBaseSpecies[intNS].maxWN = 0.;
	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	for(TransitionList::iterator tr=dBaseTrans[intNS].begin();
		 tr!= dBaseTrans[intNS].end(); ++tr)
	{
		int ipHi = tr->ipHi();
		int ipLo = tr->ipLo();
		fenergyWN = (realnum) dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN();

		tr->EnergyWN() = fenergyWN;

		if( rfield.isEnergyBound( Energy( fenergyWN, "cm^-1" ) ) )
		{
			tr->WLAng() = (realnum) wn2ang( fenergyWN );
			dBaseSpecies[intNS].maxWN = MAX2(dBaseSpecies[intNS].maxWN,fenergyWN);
		}
		else
			tr->WLAng() = 1e30;
	}

	/************************************************************************/
	/*** Read in the level numbers, Einstein A and transition wavelength  ***/
	/************************************************************************/

	//Count the number of rows first
	long wgfarows = -1;
	if (wgfastream.is_open())
	{
		int nj = 0;
		char otemp[eof_col];
		while(nj != -1)
		{
			wgfastream.get(otemp,eof_col);
			wgfastream.ignore(INT_MAX,'\n');
			if( otemp[0] == chCommentChianti ) continue;
			nj = atoi(otemp);
			wgfarows++;
		}
		wgfastream.seekg(0,ios::beg);
	}
	else 
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chTraFilename.c_str());


	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\n\nTransition Probability File: %s\n",chTraFilename.c_str());
		fprintf(ioQQQ,"Species|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Wavelength(A)|Ein A\n");
	}


	//line_index_col is the length(+1) of the level indexes in the WGFA file
	const int line_index_col = 6;
	//line_nrg_to_eina is the # of columns to skip from wavelength to eina in WGFA file
	const int line_nrg_to_eina =15;
	//line_eina_col is the length(+1) of einsteinA in WGFA
	const int line_eina_col = 16;
	char lvltemp[line_index_col];
	//Start reading WGFA file
	if (wgfastream.is_open())
	{
		for (long ii = 0;ii<wgfarows;ii++)
		{
			wgfastream.get(lvltemp,line_index_col);
			if( lvltemp[0] == chCommentChianti )
			{
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}

			long ipLoInFile = atoi(lvltemp);
			wgfastream.get(lvltemp,line_index_col);
			long ipHiInFile = atoi(lvltemp);

			// ipLo and ipHi will be manipulated below, want to retain original vals for prints
			long ipLo = ipLoInFile - 1;
			long ipHi = ipHiInFile - 1;

			if( atmdat.lgChiantiExp )
			{
				/* If either upper or lower index is -1 in the relational vector,
				 * skip that line in the wgfa file.
				 * Otherwise translate the level indices.*/
				if( intExperIndex[ipLo] == -1 || intExperIndex[ipHi] == -1 )
				{
					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}
				else
				{
					ipHi = intExperIndex.at(ipHi);
					if (ipHi < long(indexold2new.size()))
					{
						ipHi = indexold2new[ipHi];
					}
					else
					{
						ipHi = -1;
					}
					ipLo = intExperIndex.at(ipLo);
					if (ipLo < long(indexold2new.size()))
					{
						ipLo = indexold2new[ipLo];
					}
					else
					{
						ipLo = -1;
					}
				}
			}
			else
			{
				long testlo = -1, testhi = -1;

				try
				{
					testlo = indexold2new[ipLo];
					testhi = indexold2new[ipHi];
				}
				catch ( out_of_range& /* e */ )
				{
					if( DEBUGSTATE )
					{
						fprintf(ioQQQ,"NOTE: An out of range exception has occurred"
								" reading in data from %s\n",chTraFilename.c_str());
						fprintf(ioQQQ," The line in the file containing the unidentifiable"
								" levels has been ignored.\n");
						fprintf(ioQQQ,"There is no reason for alarm."
								" This message is just for documentation.\n");
					}

					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}

				if(  testlo == -1 || testhi == -1 )
				{
					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}
				else
				{
					ipLo = testlo;
					ipHi = testhi;
				}
			}

			if( ipLo >= nMolLevs || ipHi >= nMolLevs )
			{
				// skip these lines
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}
	
			if( ipHi == ipLo )
			{
				fprintf( ioQQQ," WARNING: Upper level = lower for a radiative transition in %s. Ignoring.\n", chTraFilename.c_str() );
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}

			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"<%s>\t%li:%li\t%li:%li\t",dBaseSpecies[intNS].chLabel,ipLoInFile,ipHiInFile,ipLo+1,ipHi+1);
			}
	
			ASSERT( ipHi != ipLo );
			ASSERT( ipHi >= 0 );
			ASSERT( ipLo >= 0 );

			// sometimes the CHIANTI datafiles list the highest index first as in the middle of these five lines in ne_10.wgfa:
			//    ...
			//    8   10       187.5299      0.000e+00      4.127e+05                 3d 2D1.5 -                  4s 2S0.5           E2
			//    9   10       187.6573      0.000e+00      6.197e+05                 3d 2D2.5 -                  4s 2S0.5           E2
			//   11   10   4842624.0000      1.499e-05      9.423e-06                 4p 2P0.5 -                  4s 2S0.5           E1
			//    1   11         9.7085      1.892e-02      6.695e+11                 1s 2S0.5 -                  4p 2P0.5           E1
			//    2   11        48.5157      6.787e-02      9.618e+10                 2s 2S0.5 -                  4p 2P0.5           E1
			//    ...
			// so, just set ipHi (ipLo) equal to the max (min) of the two indices.
			// NB NB NB it looks like this may depend upon whether one uses observed or theoretical energies.

			//Read in wavelengths
			char trantemp[lvl_nrg_col];
			wgfastream.get(trantemp,lvl_nrg_col);
			fWLAng = (realnum)atof(trantemp);
			if( DEBUGSTATE && atmdat.lgChiantiExp)
			{
				fprintf(ioQQQ,"%.4f\t",fWLAng);
			}

			/* \todo 2 CHIANTI labels the H 1 2-photon transition as z wavelength of zero.
			 * Should we just ignore all of the wavelengths in this file and use the
			 * difference of level energies instead. */

			if( ipHi < ipLo )
			{
				long swap = ipHi;
				ipHi = ipLo;
				ipLo = swap;
			}

			/* If the given wavelength is negative, then theoretical energies are being used.
			 * Take the difference in stored theoretical energies.
			 * It should equal the absolute value of the wavelength in the wgfa file. */
			if( fWLAng <= 0. ) // && !atmdat.lgChiantiExp )
			{
				//if( fWLAng < 0.)
					//fprintf( ioQQQ," WARNING: Negative wavelength for species %6s, indices %3li %3li \n", dBaseSpecies[intNS].chLabel, ipLo, ipHi);
				fWLAng = (realnum)(1e8/abs(dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN()));
			}

			if( DEBUGSTATE && !atmdat.lgChiantiExp)
			{
				fprintf(ioQQQ,"%.4f\t",fWLAng);
			}
			//Skip from end of Wavelength to Einstein A and read in
			wgfastream.seekg(line_nrg_to_eina,ios::cur);
			wgfastream.get(trantemp,line_eina_col);
			feinsteina = (realnum)atof(trantemp);
			if( feinsteina == 0. )
			{
				static bool notPrintedYet = true;
				if( notPrintedYet && atmdat.lgChiantiPrint)
				{
					fprintf( ioQQQ," CAUTION: Radiative rate(s) equal to zero in %s.\n", chTraFilename.c_str() );
					notPrintedYet = false;
				}
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}
			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"%.3e\n",feinsteina);
			}

			fixit("may need to do something with these rates");
			//Read in the rest of the line and look for auto
			string chLine;
			getline( wgfastream, chLine );
			TransitionList::iterator tr = dBaseTrans[intNS].begin()+ipdBaseTrans[intNS][ipHi][ipLo];
			if( chLine.find("auto") != string::npos )
			{
				if( tr->hasEmis() )
				{
					tr->Emis().AutoIonizFrac() =
						feinsteina/(tr->Emis().Aul() + feinsteina);
					ASSERT( tr->Emis().AutoIonizFrac() >= 0. );
					ASSERT( tr->Emis().AutoIonizFrac() <= 1. );
				}
				continue;
			}

			if( tr->hasEmis() )
			{
				fprintf(ioQQQ," PROBLEM duplicate transition found by atmdat_chianti in %s, "
						"wavelength=%f\n", chTraFilename.c_str(),fWLAng);
				wgfastream.close();
				cdEXIT(EXIT_FAILURE);
			}
			tr->AddLine2Stack();
			tr->Emis().Aul() = feinsteina;

			fenergyWN = (realnum)(1e+8/fWLAng);

			// \todo Check the wavelength in the file with the difference in energy levels

			tr->EnergyWN() = fenergyWN;
			if( rfield.isEnergyBound( Energy( fenergyWN, "cm^-1" ) ) )
			{
				tr->WLAng() = (realnum) wn2ang( fenergyWN );
				tr->Emis().gf() = (realnum)GetGF(tr->Emis().Aul(), tr->EnergyWN(), tr->Hi()->g());
			}
			else
				tr->WLAng() = 1e30;

			tr->setComment( db_comment_tran_levels() );
		}
	}
	else fprintf( ioQQQ, " The data file %s is corrupted .\n",chTraFilename.c_str());
	wgfastream.close();

	/* allocate space for splines */
	AtmolCollSplines[intNS].reserve(nMolLevs);
	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		AtmolCollSplines[intNS].reserve(ipHi,nMolLevs);
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			AtmolCollSplines[intNS].reserve(ipHi,ipLo,ipNCOLLIDER);
		}
	}
	AtmolCollSplines[intNS].alloc();

	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			for( long k=0; k<ipNCOLLIDER; k++ )
			{
				/* initialize all spline variables */
				AtmolCollSplines[intNS][ipHi][ipLo][k].nSplinePts = -1; 
				AtmolCollSplines[intNS][ipHi][ipLo][k].intTranType = -1;
				AtmolCollSplines[intNS][ipHi][ipLo][k].EnergyDiff = BIGDOUBLE;
				AtmolCollSplines[intNS][ipHi][ipLo][k].ScalingParam = BIGDOUBLE;
			}
		}
	}

	/************************************/
	/*** Read in the collisional data ***/
	/************************************/

	// ipCollider 0 is electrons, 1 is protons
	for( long ipCollider=0; ipCollider<=1; ipCollider++ )
	{
		string chFilename;

		if( ipCollider == ipELECTRON )
		{
			chFilename = chEleColFilename;
		}
		else if( ipCollider == ipPROTON )
		{
			if( !lgProtonData )
				break;
			chFilename = chProColFilename;
		}
		else
			TotalInsanity();

		if( DEBUGSTATE )
		{
			fprintf(ioQQQ,"\n\nCollision Data File: %s\n",chTraFilename.c_str());
			fprintf(ioQQQ,"Species|File Index (Lo:Hi)|Cloudy Index (Lo:Hi)|Spline Points\n");
		}

		fstream splupsstream;
		open_data( splupsstream, chFilename, mode_r );

		//cs_eof_col is the length(+1) of the first column used for finding the end of file
		const int cs_eof_col = 4;
		//cs_index_col is the length(+1) of the indexes in the CS file
		const int cs_index_col = 4;
		//cs_trantype_col is the length(+1) of the transition type in the CS file
		const int cs_trantype_col = 4;
		//cs_values_col is the length(+1) of the other values in the CS file
		//including: GF, nrg diff, scaling parameter, and spline points
		const int cs_values_col = 11;
		//Determine the number of rows in the CS file
		if (splupsstream.is_open())
		{
			int nj = 0;
			//splupslines is -1 since the loop runs 1 extra time
			long splupslines = -1;
			char otemp[cs_eof_col];
			while(nj != -1)
			{
				splupsstream.get(otemp,cs_eof_col);
				splupsstream.ignore(INT_MAX,'\n');
				nj = atoi(otemp);
				splupslines++;
			}
			splupsstream.seekg(0,ios::beg);
	
			for (int m = 0;m<splupslines;m++)
			{
				if( ipCollider == ipELECTRON )
				{
					splupsstream.seekg(6,ios::cur);
				}

				/* level indices */
				splupsstream.get(otemp,cs_index_col);
				long ipLo = atoi(otemp)-1;
				splupsstream.get(otemp,cs_index_col);
				long ipHi = atoi(otemp)-1;

				long ipLoFile = ipLo;
				long ipHiFile = ipHi;

				/* If either upper or lower index is -1 in the relational vector,
				* skip that line in the splups file.
				* Otherwise translate the level indices.*/
				if( atmdat.lgChiantiExp )
				{
					if( intExperIndex[ipLo] == - 1 || intExperIndex[ipHi] == -1 )
					{
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}
					else
					{
						ipHi = intExperIndex.at(ipHi);
						if (ipHi < long(indexold2new.size()))
						{
							ipHi = indexold2new[ipHi];
						}
						else
						{
							ipHi = -1;
						}
						ipLo = intExperIndex.at(ipLo);
						if (ipLo < long(indexold2new.size()))
						{
							ipLo = indexold2new[ipLo];
						}
						else
						{
							ipLo = -1;
						}
					}
				}
				else
				{
					long testlo = -1, testhi = -1;

					/* With level trimming on it is possible that there can be rows that
					 * have to be skipped when using theoretical
					 * since the levels no longer exist */
					try
					{
						testlo = indexold2new[ipLo];
						testhi = indexold2new[ipHi];
					}
					catch ( out_of_range& /* e */ )
					{
						if( DEBUGSTATE )
						{
							fprintf(ioQQQ,"NOTE: An out of range exception has occurred"
									" reading in data from %s\n",chEleColFilename.c_str());
							fprintf(ioQQQ," The line in the file containing the unidentifiable"
									" levels has been ignored.\n");
							fprintf(ioQQQ,"There is no reason for alarm."
									" This message is for documentation.\n");
						}
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}

					if( testlo == -1 || testhi == -1 )
					{
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}
					else
					{
						ipLo = testlo;
						ipHi = testhi;
					}
				}

				if( ipLo >= nMolLevs || ipHi >= nMolLevs )
				{
					// skip these transitions
					splupsstream.ignore(INT_MAX,'\n');
					continue;
				}

				if( ipHi < ipLo )
				{
					long swap = ipHi;
					ipHi = ipLo;
					ipLo = swap;
				}

				if( DEBUGSTATE )
				{
					fprintf(ioQQQ,"<%s>\t%li:%li\t%li:%li",dBaseSpecies[intNS].chLabel,ipLoFile+1,ipHiFile+1,ipLo+1,ipHi+1);
				}

				/*Transition Type*/
				splupsstream.get(otemp,cs_trantype_col);
				intTranType = atoi(otemp);
				char qtemp[cs_values_col];
				splupsstream.get(qtemp,cs_values_col);
				/*Energy Difference*/
				splupsstream.get(qtemp,cs_values_col);
				fEnergyDiff = atof(qtemp);
				/*Scaling Parameter*/
				splupsstream.get(qtemp,cs_values_col);
				fScalingParam = atof(qtemp);

				ASSERT( ipLo != ipHi );
				ASSERT( ipLo >= 0 && ipLo < nMolLevs );
				ASSERT( ipHi >= 0 && ipHi < nMolLevs );

				while( splupsstream.peek() != '\n' )
				{
					splupsstream.get(qtemp,cs_values_col);
					if( qtemp[0] == ' ' && qtemp[1] == ' ' )
						break;
					double temp = atof(qtemp);
					if( DEBUGSTATE )
					{
						fprintf(ioQQQ,"\t%.3e",temp);
					}
					// intTranType == 6 means log10 of numbers have been fit => allow negative numbers
					// intTranType < 6 means linear numbers have been fit => negative numbers are unphysical
					if( intTranType < 6 )
						temp = max( temp, 0. );
					AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline.push_back( temp );
				}

				if( DEBUGSTATE )
				{
					fprintf(ioQQQ,"\n");
				}

				intsplinepts = int( AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline.size() );
				ASSERT( intsplinepts > 2 );

				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer.resize( intsplinepts );

				/*The zeroth element contains the number of spline points*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].nSplinePts = intsplinepts;
				/*Transition type*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].intTranType = intTranType;
				/*Energy difference between two levels in Rydbergs*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].EnergyDiff = fEnergyDiff;
				/*Scaling parameter C*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].ScalingParam = fScalingParam;

				/*Once the spline points have been filled,fill the second derivatives structure*/
				/*Creating spline points array*/
				vector<double> xs (intsplinepts),
					spl(intsplinepts),
					spl2(intsplinepts);

				double coeff = (double)1/(intsplinepts-1);
				for(intxs=0;intxs<intsplinepts;intxs++)
				{
					xs[intxs] = coeff*intxs;
					spl[intxs] = AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[intxs];
				}

				spline(&xs[0], &spl[0],intsplinepts,2e31,2e31,&spl2[0]);

				/*Filling the second derivative structure*/
				for( long ii=0; ii<intsplinepts; ii++)
				{
					AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer[ii] = spl2[ii];
				}

				splupsstream.ignore(INT_MAX,'\n');
			}
			splupsstream.close();
		}
	}

	// close open file handles
	fclose( ioElecCollData );
	if( lgProtonData )
		fclose( ioProtCollData );

	//Chianti had bad theo level data so we used experimental
	//Changing lgChiantiExp back to false so next speices will use theoretical
	if( lgChiaBadTheo )
	{
		atmdat.lgChiantiExp = false;
	}

	return;
}
