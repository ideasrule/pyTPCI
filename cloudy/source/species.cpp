/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "species.h"
#include "taulines.h"
#include "input.h"
#include "dense.h"
#include "atmdat.h"
#include "elementnames.h"
#include "version.h"
#include "save.h"
#include "mole.h"
#include "service.h"
#include "parser.h"
#include "parse_species.h"
#include "prt.h"
#include "warnings.h"

/*File nemala.cpp was developed by Humeshkar B Nemala as a part of his thesis work during the Summer of 2007*/
/* Initially the code has been developed to read in energy levels,radiative and
 * collisional data from the CHIANTI and LEIDEN databases. The idea is to extend it to more databases.
 * In the case of the Leiden database there is a single .dat file which has the energy levels information,
 * radiative and collisional data, with the data corresponding to each collider coming one after the other.
 * In the case of CHIANTI, the energy levels data, radiative data and collision data are present in seperate files.
 * While LEIDEN gives collisional rate coefficients, CHIANTI gives collisional strengths.
 * In the case of CHIANTI only two colliders are used:electrons and protons. They appear as separate files.
 * The electron collision strengths files are always expected to be there. A flag is set and data processed 
 * if the file on proton collision strengths is available.*/

/* There is an initialization file called species.ini which tells Cloudy what kind of data is to be used */
/* Structures are created separately to hold the transition data,radiative and collisional data */
/* The collisional structures are different for different databases depending upon whether */
/* collisional strengths or collisional rate coefficients are used.Finally a superstructure is constructed to hold */
/* the total collisional rate obtained by considering all the colliders */
/* The colliders considered are electron,proton,Atomic Hydrogen,He,He+,He++,Ortho Molecular Hydrogen,Para Molecular Hydrogen and Molecular Hydrogen */
STATIC void states_popfill(void);
STATIC void states_nelemfill(void);
STATIC void database_prep(int);
STATIC void trim_levels(long);
STATIC void set_fractionation( species *sp );
STATIC void states_propprint(void);

static const bool DEBUGSTATE = false;

void db_basename_to_spectral( const string& chBasename, string &chElement,
				char *chSpecLabel )
{
	DEBUG_ENTRY( "db_basename_to_spectral()" );

	ASSERT( isalpha(chBasename[0]) );

	auto p = chBasename.find("_");
	chElement = chBasename.substr(0,p);
	uncaps(chElement);

	long cursor=0;
	chSpecLabel[0] = chBasename[0];
	if( isalpha(chBasename[1]) )
	{
		chSpecLabel[1] = chBasename[1];
		cursor = 2;
	}
	else
	{
		chSpecLabel[1] = ' ';
		cursor = 1;
	}

	ASSERT( chBasename[cursor]=='_' );
	++cursor;
	ASSERT( isdigit(chBasename[cursor]) );

	if( isdigit(chBasename[cursor+1]) )
	{
		chSpecLabel[2] = chBasename[cursor++];
		chSpecLabel[3] = chBasename[cursor++];
	}
	else
	{
		chSpecLabel[2] = ' ';
		chSpecLabel[3] = chBasename[cursor++];
	}
	chSpecLabel[4] = '\0';
	ASSERT( chBasename[cursor]=='\0' || chBasename[cursor]=='d' );

	// now capitalize the first letter
	chSpecLabel[0] = toupper( chSpecLabel[0] );
}

void database_readin( void )
{
	int i,intNoSp;

	const int MAX_NUM_SPECIES = 1000;

	char chLabels[MAX_NUM_SPECIES][CHARS_SPECIES];
	char chLabelsOrig[MAX_NUM_SPECIES][CHARS_SPECIES];
	string chPaths[MAX_NUM_SPECIES];

	static int nCalled = 0;
	long nSpeciesLAMDA, nSpeciesSTOUT, nSpeciesCHIANTI;

	DEBUG_ENTRY( "database_readin()" );

	/* only do this once. */
	if( nCalled > 0 )
	{
		return;
	}

	/* this is first call, increment the nCalled counter so never do this again */
	++nCalled;

	// read masterlists, count number of species
	nSpecies = 0;

	////////////////////////////////////
	//  
	// Read LAMDA masterlist 
	//
	//////////////////////////////////

	/* count how many lines are in the file, ignoring all lines
	 * starting with '#':This would give the number of molecules */
	nSpeciesLAMDA = 0;

	string sep = cpu.i().chDirSeparator();

	if( atmdat.lgLamdaOn )
	{
		long numModelsNotUsed = 0;
		string chPath = "lamda" + sep + "masterlist" + sep + atmdat.chLamdaFile;
		DataParser d( chPath, ES_NONE );

		while( d.getline() )
		{	
			string chToken;
			d.getToken( chToken );
			if( findspecies( chToken.c_str() ) != null_mole || 
			    ( chToken.length() > 2 && chToken[1]=='-' &&
			      findspecies( chToken.c_str()+2 ) != null_mole ) )
			{
				ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
				ASSERT( nSpeciesLAMDA + 1 <= MAX_NUM_SPECIES );
				ASSERT( chToken.length() < CHARS_SPECIES );
				strcpy( chLabels[nSpecies], chToken.c_str() );
				chLabels[nSpecies][CHARS_SPECIES-1] = '\0';

				// path is, for example, lamda/no.dat
				d.getToken( chToken );
				chPaths[nSpecies] = "lamda" + sep + chToken;
				d.checkEOL();
				++nSpecies;
				++nSpeciesLAMDA;
			}
			else
				++numModelsNotUsed;
		}

		if( !t_version::Inst().lgRelease && prt.lgPrintTime && numModelsNotUsed > 0 )
		{
			ostringstream oss;
			oss << "   " << numModelsNotUsed << " LAMDA models could not be found in chemistry network.";
			notein(oss.str());
		}
	}

	/* Print LAMDA molecule list if save data sources is on*/
	if( save.lgSDSOn && atmdat.lgLamdaOn)
	{
		fprintf(save.ipSDSFile, "##################################################\n");
		fprintf( save.ipSDSFile,"LAMDA (2005, A&A, 432, 369) molecules in this run.\n");
		for( int i=0; i<nSpeciesLAMDA; i++)
		{
			fprintf( save.ipSDSFile,"%s\t\t",chLabels[i]);
			if( (i+1)%5 == 0 )
			{
				fprintf( save.ipSDSFile, "\n");
			}
		}
		fprintf(save.ipSDSFile,"\n\n");
	}

	//////////////////////////////////
	//  
	// Read CDMS/JPL masterlist
	//
	// These data files are in LAMDA format
	//
	//////////////////////////////////

	if( atmdat.lgCalpgmOn )
	{
		long numModelsNotUsed = 0;
		string chPath = "cdms+jpl" + sep + "masterlist";
		DataParser d( chPath, ES_NONE );

		while( d.getline() )
		{	
			string chToken;
			d.getToken( chToken );
			// hacks for alternative dialects...
			if( chToken == "SH" )
				chToken = "HS";
			if( chToken == "SH+" )
				chToken = "HS+";
			if( chToken == "SD" )
				chToken = "DS";
			if( chToken == "CCH" )
				chToken = "C2H";
			if( chToken == "CCD" )
				chToken = "C2D";
			if( chToken == "^17OO" )
				chToken = "O^17O";
			if( chToken == "H^18O" )
				chToken = "^18OH";
			if( chToken == "HCCD" )
				chToken = "C2HD";
			if( chToken == "^13CCCH" )
				chToken = "^13CC2H";
			if( chToken == "CC^13CH" )
				chToken = "C2^13CH";
			if( chToken == "H^13CCCN" )
				chToken = "H^13CC2N";
			if( chToken == "HCC^13CN" )
				chToken = "HC2^13CN";
			if( chToken == "HCCC^15N" )
				chToken = "HC3^15N";
			// this molecule is cyclic, so these two are identical
			if( chToken == "Si^13CC" )
				chToken = "SiC^13C";
			if( findspecies( chToken.c_str() ) != null_mole  || 
			    ( chToken.length() > 2 && chToken[1]=='-' &&
			      findspecies( chToken.c_str()+2 ) != null_mole ) )
			{
				ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
				ASSERT( nSpeciesLAMDA + 1 <= MAX_NUM_SPECIES );
				strcpy( chLabels[nSpecies], chToken.c_str() );
				chLabels[nSpecies][CHARS_SPECIES-1] = '\0';

				d.getToken( chToken );
				chPaths[nSpecies] = "cdms+jpl" + sep + chToken;
				d.checkEOL();
				++nSpecies;
				++nSpeciesLAMDA;
			}
			else
				++numModelsNotUsed;
		}

		if( !t_version::Inst().lgRelease && prt.lgPrintTime && numModelsNotUsed > 0 )
		{
			ostringstream oss;
			oss << "   " << numModelsNotUsed << " CDMS+JPL models could not be found in chemistry network.";
			notein(oss.str());
		}
	}

	////////////////////////////////////
	//
	// Read STOUT masterlist and VERSION
	//
	///////////////////////////////////
	nSpeciesSTOUT = 0;

	//numLevels: index is nSpecies, value is the number of levels
	vector<long> numLevels(MAX_NUM_SPECIES,0L);

	if( atmdat.lgStoutOn )
	{
		// default location of Stout masterlist file
		string chPath = "stout" + sep + "masterlist" + sep + atmdat.chStoutFile;
		DataParser d( chPath, ES_NONE );

		// magic number
		d.getline();
		static const long int nYrST = 11 , nMonST = 10, nDayST = 25;
		d.checkMagic( nYrST, nMonST, nDayST );

		while( d.getline() )
		{
			string chToken;
			// delimiters for tokens are whitespace chars
			// species name can have any number of columns
			// we will split line into two tokens, the name,
			// and the optional number of levels
			d.getToken( chToken );

			ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
			ASSERT( nSpeciesSTOUT + 1 <= MAX_NUM_SPECIES );

			// first token is the species name
			strcpy( chLabels[nSpecies], chToken.c_str() );
			strcpy( chLabelsOrig[nSpecies], chLabels[nSpecies] );

			// second optional token is lower limit to number of levels
			long numLevs;
			if( d.getTokenOptional( numLevs ) )
			{
				if( numLevs > 0 )
					numLevels[nSpecies] = numLevs;
				else
					d.errorAbort( "the limit to the number of levels must be positive" );
			}

			d.checkEOL();

			bool skipSpecies = false;

			//Check for duplicate species within Stout masterlist
			for( int j = nSpeciesLAMDA; j < nSpecies; j++ )
			{
				if( strcmp( chLabelsOrig[j], chLabelsOrig[nSpecies] ) == 0 )
				{
					fprintf(ioQQQ,"%s appears multiple times in %s.\n",
						chLabels[nSpecies],atmdat.chStoutFile);
					skipSpecies = true;
					break;
				}
			}

			if( skipSpecies )
				continue;

			string chElement;
			db_basename_to_spectral( chToken, chElement, chLabels[nSpecies] );

			if( DEBUGSTATE )
				fprintf( ioQQQ, "Stout:\t '%s' -> '%s'\n",
						chToken.c_str(), chLabels[nSpecies] );

			// path is, for example, STOUT/ar/ar_10/ar_10
			// we will append extensions later
			chPaths[nSpecies] = "stout" + sep + chElement + sep + chToken + sep + chToken;

			++nSpecies;
			++nSpeciesSTOUT;
		}
	}


	////////////////////////////////////
	//  
	// Read CHIANTI masterlist and VERSION
	//
	///////////////////////////////////

	nSpeciesCHIANTI = 0;

	if( atmdat.lgChiantiOn )
	{
		string chPath = "chianti" + sep + "VERSION";
		DataParser d;
		d.open( chPath, ES_NONE );
		d.getline();
		// chianti version - string since can contain letters
		d.getToken(atmdat.chVersion);
		d.checkEOL();

		chPath = "chianti" + sep + "masterlist" + sep + atmdat.chCloudyChiantiFile;
		d.open( chPath, ES_NONE );
		d.getline();
		/* magic numbers for this version of Chianti masterlist */
		static const long int nYr = 22, nMon = 7, nDay = 9;
		d.checkMagic( nYr, nMon, nDay );

		while( d.getline() )
		{
			// break line into two chunks, first with species which can have number number of
			// characters, followed by optional chunk with limit to number of levels
			string chToken;
			d.getToken( chToken );

			fixit("insert logic here to exclude some ions");
			// (for example, iso sequences)
			// exclude for now the satellite lines (denoted by a "d" after the label
			if( chToken.back() != 'd' )
			{
				ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
				ASSERT( nSpeciesCHIANTI + 1 <= MAX_NUM_SPECIES );
				strcpy( chLabels[nSpecies], chToken.c_str() );
				strcpy( chLabelsOrig[nSpecies], chLabels[nSpecies]);

				// second optional token is lower limit to number of levels
				// first get full string after arbitrary length species name
				long numLevs;
				// was there a lower bound to the number of levels
				if( d.getTokenOptional( numLevs ) )
				{
					if( numLevs > 0 )
						numLevels[nSpecies] = numLevs;
					else
						d.errorAbort( "the limit to the number of levels must be positive" );
				}

				d.checkEOL();

				bool skipSpecies = false;

				// Check for duplicate species with Stout
				for( int j = nSpeciesLAMDA; j < (nSpecies - nSpeciesCHIANTI); j++ )
				{
					if( strcmp( chLabelsOrig[j], chLabelsOrig[nSpecies] ) == 0 )
					{
						fprintf(ioQQQ,"Skipping the Chianti version of %s, using Stout version\n",
							chLabels[nSpecies]);
						skipSpecies = true;
						break;
					}
				}
				// Check for duplicate species within Chianti masterlist
				for( int j = nSpecies - nSpeciesCHIANTI; j < nSpecies; j++ )
				{
					if( strcmp( chLabelsOrig[j], chLabelsOrig[nSpecies] ) == 0 )
					{
						fprintf(ioQQQ,"%s appears multiple times in %s.\n",
							chLabels[nSpecies],atmdat.chCloudyChiantiFile);
						skipSpecies = true;
						break;
					}
				}
				if( skipSpecies )
					continue;

				string chElement;
				db_basename_to_spectral( chToken, chElement, chLabels[nSpecies] );
				if( DEBUGSTATE )
					fprintf( ioQQQ, "Chianti:\t '%s' -> '%s'\n",
							chToken.c_str(), chLabels[nSpecies] );

				// path is, for example, CHIANTI/ar/ar_10/ar_10
				// we will append extensions later
				chPaths[nSpecies] = "chianti" + sep + chElement + sep + chToken + sep + chToken;

				++nSpecies;
				++nSpeciesCHIANTI;
			}
		}
	}

	/* no species found, nothing to do */
	if( nSpecies==0 )
		return;

	/*Initialization of the dBaseSpecies Structure*/
	dBaseSpecies.resize(nSpecies);

	/*Initialization of the collisional rates array structure*/
	AtmolCollRateCoeff.reserve( nSpecies );
	AtmolCollSplines.resize(nSpecies);
	StoutCollData.resize(nSpecies);

	/*allocating here takes care of the number of colliders*/
	for( i=0; i<nSpecies; i++ )
	{
		AtmolCollRateCoeff.reserve( i, ipNCOLLIDER );
	}
	AtmolCollRateCoeff.alloc();

	// allocate state and transition arrays
	dBaseStates.resize(nSpecies);
	ipdBaseTrans.resize(nSpecies);

	for( i = 0; i < nSpecies; i++ )
	{
		dBaseTrans.push_back(TransitionList("dBaseTrans",&dBaseStates[i]));
		// label should be a minimum of 4 characters long
		size_t los = strlen(chLabels[i]);
		ASSERT( los >= 1 && los <= CHARS_SPECIES );
		dBaseSpecies[i].chLabel = new char[los+1];
		strcpy(dBaseSpecies[i].chLabel,chLabels[i]);
		dBaseSpecies[i].chLabel[los]='\0';
		trimTrailingWhiteSpace( dBaseSpecies[i].chLabel );
		dBaseSpecies[i].lgActive = true;

		// was minimum number of levels specified
		if( numLevels[i] > 0 )
		{
			dBaseSpecies[i].numLevels_masterlist = numLevels[i];
		}

		/* set type and isotopologue fractions */
		set_fractionation( &dBaseSpecies[i] );

		// set_fractionation trims off "p-","o-", etc.  Now have set label.  Check size.
		los = (int)strlen( dBaseSpecies[i].chLabel );
		ASSERT( los < CHARS_SPECIES );

		if( i<nSpeciesLAMDA )
		{
			// Read in LAMDA data files
			atmdat_LAMDA_readin( i, chPaths[i] );
		}
		else if( i < nSpeciesLAMDA + nSpeciesSTOUT )
		{
			// Read in STOUT data files
			atmdat_STOUT_readin( i, chPaths[i] );
		}
		else if( i < nSpeciesLAMDA + nSpeciesSTOUT + nSpeciesCHIANTI )
		{
			// Read in CHIANTI data files
			atmdat_CHIANTI_readin( i, chPaths[i] );
		}
		else
			TotalInsanity();
	}

	speciesCheck();

	states_popfill();
	states_nelemfill();

	for( long i=nSpeciesLAMDA; i<nSpeciesLAMDA+nSpeciesSTOUT; i++ )
	{
		strcpy(atmdat.chdBaseSources[dBaseStates[i][0].nelem()-1][dBaseStates[i][0].IonStg()-1],"Stout");
		atmdat.lgdBaseSourceExists[dBaseStates[i][0].nelem()-1][dBaseStates[i][0].IonStg()-1] = true;
	}
	for( long i=nSpeciesLAMDA+nSpeciesSTOUT; i<nSpeciesLAMDA+nSpeciesSTOUT+nSpeciesCHIANTI; i++ )
	{
		strcpy(atmdat.chdBaseSources[dBaseStates[i][0].nelem()-1][dBaseStates[i][0].IonStg()-1],"Chianti");
		atmdat.lgdBaseSourceExists[dBaseStates[i][0].nelem()-1][dBaseStates[i][0].IonStg()-1] = true;
	}

	if( save.lgSDSOn )
	{
		fprintf(save.ipSDSFile, "##################################################\n");
		fprintf( save.ipSDSFile,"Atomic model for each species used in this run.\n");
		fprintf( save.ipSDSFile,"Chianti (C), Stout(S), Iso-sequences (I), old internal treatment( ).\n\n");

		fprintf( save.ipSDSFile,"Ion");
		for(int i=0; i<LIMELM; i++)
		{
			fprintf( save.ipSDSFile,"%4d",i);
		}
		fprintf( save.ipSDSFile,"\n");

		for( int i=0; i<LIMELM; i++)
		{
			fprintf( save.ipSDSFile,"%s ",elementnames.chElementSym[i]);
			for(int j=0; j<i+1; j++)
			{
				fprintf( save.ipSDSFile,"   %c",atmdat.chdBaseSources[i][j][0]);
			}
			fprintf( save.ipSDSFile,"\n");
		}
	}

	if( DEBUGSTATE )
	{
		fprintf(ioQQQ,"\n\nDEBUG: Below are the contents of chdBaseSources[][]. It should contain a database name for each species.\n");
		fprintf( ioQQQ,"Ion");
		for(int i=0; i<LIMELM; i++)
		{
			fprintf( ioQQQ,"\t%i",i);
		}
		fprintf( ioQQQ,"\n");
		for( int i = 0; i < LIMELM; i++ )
		{
			fprintf( ioQQQ,"%s",elementnames.chElementSym[i]);
			for( int j = 0; j < LIMELM+1; j++ )
			{
				fprintf(ioQQQ,"\t%s",atmdat.chdBaseSources[i][j]);
			}
			fprintf(ioQQQ, "\n");
		}
		fprintf(ioQQQ,"\n\n");
	}

	/*Setting nelem of the states to an arbitrary value*/
	/*Also trim the highest levels if there are no valid Auls */
	/*Loop over species*/
	for( intNoSp=0; intNoSp<nSpecies; intNoSp++ )
	{
		database_prep(intNoSp);
		AllTransitions.push_back(dBaseTrans[intNoSp]);
		trim_levels(intNoSp);
	}

	/*To print the states*/
	if(DEBUGSTATE)
		states_propprint();
	return;
}

/** trim_levels: Trim the highest energy levels until there is a transition out of the highest
 * level with an Aul value that is not the default. The default value of 1e-30 is used when no Aul is
 * given in the database. */
STATIC void trim_levels(long ipSpecies)
{
	DEBUG_ENTRY( "trim_levels()" );

	bool lgLevelsToTrim = true;
	bool lgLevelsTrimmed = false;
	char* spectralLabel = dBaseSpecies[ipSpecies].chLabel;
	string speciesLabel = dBaseStates[ipSpecies].chLabel();

	long totalNumLevels = dBaseSpecies[ipSpecies].numLevels_max;

	while( lgLevelsToTrim )
	{
		long ipHi = dBaseSpecies[ipSpecies].numLevels_max-1;
		lgLevelsToTrim = true;

		if( dBaseSpecies[ipSpecies].numLevels_max == 0)
		{
			fprintf(ioQQQ,"PROBLEM: Spectrum %s (species: %s) has no transition probabilities out of the first %li levels.\n",
				spectralLabel, speciesLabel.c_str(), totalNumLevels);
			fprintf(ioQQQ,"Consider allowing Cloudy to use more levels (see Hazy 1 SPECIES \"NAME\" LEVELS=ALL), add more low-level"
					" transition probabilities, or disable %s in the masterlist.\n\n", spectralLabel);
			cdEXIT(EXIT_FAILURE);
		}

		for( int ipLo = 0; ipLo < ipHi; ipLo++)
		{
			TransitionList::iterator tr = dBaseTrans[ipSpecies].begin()+ipdBaseTrans[ipSpecies][ipHi][ipLo];
			if( DEBUGSTATE )
			{
				fprintf(ioQQQ,"trim_levels():\t%s\t%i\t%li\t%e\n", spectralLabel, ipLo+1, ipHi+1, tr->Emis().Aul());
			}

			if( tr->Emis().Aul() > atmdat.aulThreshold )
			{
				lgLevelsToTrim = false;
				break;
			}

		}
		if( lgLevelsToTrim )
		{
			--dBaseSpecies[ipSpecies].numLevels_max;
			dBaseSpecies[ipSpecies].numLevels_local = dBaseSpecies[ipSpecies].numLevels_max;
			lgLevelsTrimmed = true;
		}
	}

	if( lgLevelsTrimmed && ( atmdat.lgChiantiPrint || atmdat.lgStoutPrint || atmdat.lgLamdaPrint ) )
	{
		fprintf(ioQQQ,"Spectrum %s (species: %s) trimmed to %li levels (original %li) to have positive Aul.\n",
				spectralLabel,
				speciesLabel.c_str(),
				dBaseSpecies[ipSpecies].numLevels_local,
				totalNumLevels);
	}
}

STATIC void set_fractionation( species *sp )
{
	DEBUG_ENTRY( "set_fractionation()" );

	char chToken[3];

	sp->fracIsotopologue = 1.f;
	//types include "p-", "o-", "e-", and "a-"
	strncpy( chToken, sp->chLabel, 2 );
	chToken[2] = '\0';
	if( strcmp( "p-", chToken )==0 )
		sp->fracType = 0.25f;
	else if( strcmp( "o-", chToken )==0 )
		sp->fracType = 0.75f;
	else if( strcmp( "e-", chToken )==0 )
		sp->fracType = 0.5f;
	else if( strcmp( "a-", chToken )==0 ) 
		sp->fracType = 0.5f;
	else
		sp->fracType = 1.0f;

	fixit("what fraction should e-type and a-type Methanol have?  Assume 50/50 for now.");

	// Now scrape the type specifier off the label.
	if( sp->chLabel[1]=='-')
		memmove(sp->chLabel,sp->chLabel+2,strlen(sp->chLabel+2)+1);

	return;
}

/*This function zeros the population of all states */
STATIC void states_popfill( void)
{
	DEBUG_ENTRY( "states_popfill()" );

	for( long i=0; i<nSpecies; i++)
	{
		for( long j=0; j<dBaseSpecies[i].numLevels_max; j++)
		{
			dBaseStates[i][j].Pop() = 0.;
		}
	}
	return;
}

void parsespect(const char* chLabel, long& nelem, long& IonStg)
{
	DEBUG_ENTRY( "parsespect()" );	
	nelem = -1;
	IonStg = -1;
	if ( strlen(chLabel) != 4 || 
		  ! (isalpha(chLabel[0])) ||
		  ! (chLabel[1] == ' ' || isalpha(chLabel[1])) ||
		  ! (chLabel[2] == ' ' || isdigit(chLabel[2])) ||
		  ! (chLabel[3] == ' ' || isdigit(chLabel[3])))
	{
		// Invalid spectrum -- return error state
		return;
	}
	char chToken[3];
	strncpy( chToken, chLabel, 2 );
	chToken[2] = '\0';
	for( long ipElement=0; ipElement<LIMELM; ipElement++ )
	{
		if( strcmp( elementnames.chElementSym[ipElement], chToken )==0 )
		{
			nelem = ipElement;
			break;
		}
	}
	strncpy( chToken, chLabel + 2, 2 );
	IonStg = stol(chToken);
}

bool isAtomicIonValid( const long element_index, const long spCharge )
{
	DEBUG_ENTRY( "isAtomicIonValid()" );

	ASSERT( element_index >= 0 );

	if( spCharge < 0 || spCharge > element_index+1 )
	{
		return false;
	}

	return true;
}

bool isBareNucleus( const long element_index, const long spCharge )
{
	DEBUG_ENTRY( "isBareNucleus()" );

	ASSERT( element_index >= 0 );

	if( spCharge == element_index+1 )
	{
		return true;
	}

	return false;
}

string makeChemical( long nelem, long ion )
{
	DEBUG_ENTRY("makeChemical()");

	string chLabelChemical = elementnames.chElementSym[nelem];
	if( elementnames.chElementSym[nelem][1]==' ' )
		chLabelChemical = elementnames.chElementSym[nelem][0];

	// size has to be > 20 to prevent warnings about writing into an array that may be too small
	char chStage[21] = {'\0'};
	if( ion==1 )
		chStage[0] = '+';
	else if( ion>1 )
		sprintf( chStage, "+%li", ion );

	return chLabelChemical + chStage;
}

void makeChemical(char* chLabelChemical, long nelem, long ion)
{
	DEBUG_ENTRY("makeChemical()");
	string chemLab = makeChemical( nelem, ion );
	strncpy( chLabelChemical, chemLab.c_str(), size_t(CHARS_SPECIES-1) );
	// make sure the string is properly terminated...
	chLabelChemical[CHARS_SPECIES-1] = '\0';
}

STATIC void spectral_to_chemical( char *chLabelChemical, const char* chLabel, long &nelem, long &IonStg )
{
	DEBUG_ENTRY( "spectral_to_chemical()" );

	parsespect( chLabel, nelem, IonStg );
	ASSERT( nelem >= 0 && nelem < LIMELM );
	ASSERT( IonStg >= 1 && IonStg <= nelem+2 );

	//Prevent importing of iso-sequences from Chianti
	if( nelem - (IonStg-1) < NISO )
	{
		fprintf(ioQQQ, " PROBLEM: Cannot use Chianti model for %s%li\n",elementnames.chElementSym[nelem],IonStg);
		fprintf(ioQQQ, "  Iso-sequences are handled by our own model.\n");
		cdEXIT(EXIT_FAILURE);
	}

	makeChemical(chLabelChemical, nelem, IonStg-1);

	return;
}

void spectral_to_chemical( char *chLabelChemical, const char* chLabel )
{
	DEBUG_ENTRY( "spectral_to_chemical()" );

	long nelem, IonStg;
	return spectral_to_chemical( chLabelChemical, chLabel, nelem, IonStg );
}

void spectral_to_chemical( string &chemicalLabel, const char *chLabel )
{
	DEBUG_ENTRY( "spectral_to_chemical()" );

	char speciesLabel[NCHLAB] = "";
	spectral_to_chemical( speciesLabel, chLabel );
	chemicalLabel = speciesLabel;
}

bool parse_chemical( const string &chLabelChem,
				string &chElemMol, long &spCharge )
{
	DEBUG_ENTRY( "parse_chemical()" );

	chElemMol = "";
	spCharge = -1;

	if( chLabelChem == "" )
		return false;

	size_t plus_sign_pos = chLabelChem.find_first_of( '+' );

	{
		enum{ DEBUG_LOC = false };
		if( DEBUG_LOC )
		{
			fprintf( ioQQQ,
				"species= '%s'\t plus_sign_pos= %d\t length= %d\n",
				chLabelChem.c_str(),
	      			int( plus_sign_pos ),
				int( chLabelChem.length() ) );
		}
	}

	if( plus_sign_pos == string::npos )
	{
		chElemMol = chLabelChem;
		spCharge = 0;
	}
	else if( plus_sign_pos+1 == chLabelChem.length() )
	{
		chElemMol = chLabelChem.substr( 0, plus_sign_pos );
		spCharge = 1;
	}
	else
	{
		chElemMol = chLabelChem.substr( 0, plus_sign_pos );
		spCharge = stol( chLabelChem.substr( plus_sign_pos+1 ) );
	}

	return true;
}

void chemical_to_spectral( const string &chLabelChem, string &chLabelSpec )
{
	DEBUG_ENTRY( "chemical_to_spectral()" );

	string chElemMol;
	long spCharge;

	if( parse_chemical( chLabelChem, chElemMol, spCharge ) )
	{
		if( spCharge == 0 )
		{
			/* Both 'H' and 'HCl' go through this branch */
			chLabelSpec = chElemMol;

			if( isElementSym( chElemMol ) )
			{
				if( chElemMol.length() == 1 )
				{
					chLabelSpec += " ";
				}

				chLabelSpec += " 1";
			}
		}
		else
		{
			/* Both 'C+2' and 'LiH+' go through this branch */
			if( isElementSym( chElemMol ) )
			{
				chLabelSpec = chElemMol;

				if( chElemMol.length() == 1 )
				{
					chLabelSpec += " ";
				}

				int nelem = elem_symbol_to_index( chElemMol );
				if( ! isAtomicIonValid( nelem, spCharge ) || isBareNucleus( nelem, spCharge ) )
				{
					chLabelSpec = "";
					return;
				}

				// +1 accounts for the zero-offset of species notation
				if( spCharge+1 < 10 )
					chLabelSpec += " ";
				stringstream tmp;
				tmp << spCharge+1;
				chLabelSpec += tmp.str();
			}
			else
			{
				chLabelSpec = chLabelChem;
			}
		}
	}
	else
	{
		chLabelSpec = "";
	}
}

/* Tell if the given species is enabled */
bool isSpeciesActive( const string &chLabelChem )
{
	DEBUG_ENTRY( "isSpeciesActive()" );

	for( size_t i=0; i<mole_global.list.size(); ++i )
	{
		if( mole.species[i].levels == NULL )
			continue;
		if( mole.species[i].dbase == NULL )
			continue;
		if( mole.species[i].levels->chLabel() == chLabelChem )
			return mole.species[i].dbase->lgActive;
	}

	return false;
}

/*This function fills the nelem and IonStg fields */
STATIC void states_nelemfill(void)
{
	DEBUG_ENTRY( "states_nelemfill()" );

	for( long i=0; i<nSpecies; i++ )
	{
		long nelem = 0, IonStg;
		char chLabelChemical[CHARS_SPECIES] = "";

		if( dBaseSpecies[i].lgMolecular )
		{
			fixit("should never be used if lgMolecular"); 
			/* these should never be used if lgMolecular
			 *set to dangerous values instead of unity. */
			nelem = -1;
			IonStg = -1;
			strcpy( chLabelChemical, dBaseSpecies[i].chLabel );

			dBaseSpecies[i].lgImgMatrix = false;
			if( save.img_matrix.species == dBaseStates[i].chLabel() )
			{
				dBaseSpecies[i].lgImgMatrix = true;
			}
		}
		else
		{
			spectral_to_chemical( chLabelChemical, dBaseSpecies[i].chLabel, nelem, IonStg );
			dBaseStates[i].chLabel_set( chLabelChemical );

			dBaseSpecies[i].fmolweight = dense.AtomicWeight[nelem];

			// do not evaluate our cooling if we are using Chianti for this species
			if( dBaseSpecies[i].database == "Chianti" )
			{
				dense.lgIonChiantiOn[nelem][IonStg-1] = true;
			}
			else if( dBaseSpecies[i].database == "Stout" )
			{
				dense.lgIonStoutOn[nelem][IonStg-1] = true;
			}
			else
			{
				TotalInsanity();
			}

			if( atmdat.lgChiantiLvl2Hybrid || atmdat.lgStoutLvl2Hybrid )
			{
				// used in cool_dima to indicate whether to include line
				// with shorter wl than these databases
				dense.maxWN[nelem][IonStg-1] = dBaseSpecies[i].maxWN;
			}
			else
			{
				dense.maxWN[nelem][IonStg-1] = 0.;
			}

			//Store the value of ipSpecies on C-scale
			//nelem(H) = 0  IonStg(H I) = 0
			atmdat.ipSpecIon[nelem][IonStg-1] = i;

			dBaseSpecies[i].lgPrtMatrix = false;
			if( prt.matrix.species == dBaseStates[i].chLabel() )
			{
				dBaseSpecies[i].lgPrtMatrix = true;
			}

			dBaseSpecies[i].lgImgMatrix = false;
			if( save.img_matrix.species == dBaseStates[i].chLabel() )
			{
				dBaseSpecies[i].lgImgMatrix = true;
			}
		}

		molecule *sp = findspecies(chLabelChemical);
		if( sp == null_mole )
		{
			dBaseSpecies[i].index = INT_MAX;
			if( nelem >= ipHYDROGEN && dense.lgElmtOn[nelem] )
				fprintf(ioQQQ," PROBLEM: could not find species %li - %s\n",i,
					chLabelChemical );
		}
		else
		{
			dBaseSpecies[i].index = sp->index;
			mole.species[ sp->index ].dbase = &dBaseSpecies[i];
			mole.species[ sp->index ].levels = &dBaseStates[i];
			mole.species[ sp->index ].lines = &dBaseTrans[i];			
		}
		
		for( long j=0; j<dBaseSpecies[i].numLevels_max; j++ )
		{
			dBaseStates[i][j].nelem() = nelem+1;
			dBaseStates[i][j].IonStg() = IonStg;
		}
	}
	return;
}

/*This function prints the various properties of states*/
STATIC void states_propprint(void)
{
	DEBUG_ENTRY( "states_propprint()" );

	for( long i=0; i<nSpecies; i++ )
	{
		printf("The species is %s \n",dBaseSpecies[i].chLabel);
		printf("The data output is in the following format \n");
		printf("Label Energy St.wt Pop Lifetime\n");

		for( long j=0; j<dBaseSpecies[i].numLevels_max; j++ )
		{
			printf("This is the %ld state \n",j);
			printf("%s %f %f %f %e \n",dBaseStates[i][j].chLabel().c_str(),
				dBaseStates[i][j].energy().WN(),
				dBaseStates[i][j].g(),
				dBaseStates[i][j].Pop(),
				dBaseStates[i][j].lifetime());
		}
	}
	return;
}

STATIC void database_prep(int intSpIndex)
{
	vector<realnum> fsumAs(dBaseSpecies[intSpIndex].numLevels_max,SMALLFLOAT);

	DEBUG_ENTRY( "database_prep()" );

	/*Get the lifetimes*/
	for( EmissionList::iterator em = dBaseTrans[intSpIndex].Emis().begin();
		  em != dBaseTrans[intSpIndex].Emis().end(); ++em) 
	{
		fsumAs[(*em).Tran().ipHi()] += (*em).Aul();

		// set redistribution functions for all lines
		if( intSpIndex != atmdat.ipSpecIon[ipIRON][1] )
		{
			// default for species is partial redisctribution with wings
			(*em).iRedisFun() = ipPRD;
		}
		else
		{
			// this Fe II is to be the default for cloudy post 2001
			if( em->Tran().ipLo() == 0 )
			{
				// complete redistribution, only core
				em->iRedisFun() = ipCRD;
			}
			else
			{
				/* >>chng 01 feb 27, had been -1, crd with core only,
				 * change to complete redistribution with wings as per discussion with Ivan Hubeny */
				em->iRedisFun() = ipCRDW;
			}
		}
	}
	
	dBaseStates[intSpIndex][0].lifetime()= BIGFLOAT;
	for( int ipHi=1; ipHi < dBaseSpecies[intSpIndex].numLevels_max; ipHi++ )
	{
		dBaseStates[intSpIndex][ipHi].lifetime() = 1./fsumAs[ipHi];
	}
	return;
}
