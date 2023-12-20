/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atmdat_readin read in some data files, but only if this is very first call, 
 * called by Cloudy */
#include "cddefines.h"
#include "taulines.h"
#include "mewecoef.h"
#include "iterations.h"
#include "heavy.h"
#include "yield.h"
#include "trace.h"
#include "lines.h"
#include "struc.h"
#include "dynamics.h"
#include "elementnames.h"
#include "hyperfine.h"
#include "atmdat.h"
#include "iso.h"
#include "save.h"
#include "mole.h"
#include "two_photon.h"
#include "dense.h"
#include "lines_service.h"
#include "parser.h"

/* definition for whether level 2 lines are enabled, will be set to -1 
 * with no level2 command */
/*long nWindLine = NWINDDIM;*/
/*realnum TauLine2[NWINDDIM][NTA];*/
/*realnum **TauLine2;*/

// NB NB - IS_TOP should always be the last entry of this enum!
typedef enum { IS_NONE, IS_K_SHELL, IS_L1_SHELL, IS_L2_SHELL, IS_TOP } exc_type;

struct t_BadnellLevel
{
	string config;
	int irsl;
	int S;
	int L;
	int g; // 2*J+1
	realnum energy; // in cm^-1
	bool lgAutoIonizing;
	exc_type WhichShell;
	t_BadnellLevel() : irsl(0), S(0), L(0), g(0), energy(0.f), lgAutoIonizing(false), WhichShell(IS_NONE) {}
};

STATIC void init_struc();
STATIC void read_level2_lines();
STATIC void read_mewe_gbar();
STATIC void read_Hummer_Storey();

// read Hummer and Storey 98 He1 photoionization cross-section data
STATIC void read_SH98_He1_cross_sections(void);

STATIC void read_Helike_cross_sections(void);

// read autoionization data from Badnell data file
STATIC void ReadBadnellAIData(const string& fnam,      // filename containing the Badnell data
			      long nelem,              // nelem is on C scale
			      long ion,                // ion is on C scale
			      TransitionList& UTA,     // UTA lines will be pushed on this stack
			      bitset<IS_TOP> Skip);    // option to skip transitions from a particular shell

STATIC void validate_magic_number_1arg( const char *chFilename, FILE *ioFile,
				const long magicExp );
STATIC void validate_magic_number_3arg( const char *chFilename, FILE *ioFile,
				const long yearExp, const long monthExp, const long dayExp );
STATIC void read_UTA_lines();
STATIC void read_ionization_potentials();

// simple helper functions for ReadBadnellAIData
inline void InitTransition(const TransitionProxy& t);
inline int irsl2ind(vector<t_BadnellLevel>& level, int irsl);

// UTA lines below this absorption oscillator strength value will be ignored -
// F+13 paper plotted f not gf
const realnum f_cutoff = 1.e-4f;

void atmdat_readin(void)
{
	static bool lgFirstCall = true;

	DEBUG_ENTRY( "atmdat_readin()" );

	/* do nothing if not first call */
	if( !lgFirstCall )
	{
		/* do not do anything, but make sure that number of zones has not increased */
		bool lgTooBig = false;
		for( long j=0; j < iterations.iter_alloc; j++ )
		{
			if( iterations.nend[j]>=struc.nzlim )
				lgTooBig = true;
		}
		if( lgTooBig )
		{
			fprintf(ioQQQ," This is the second or later calculation in a grid.\n");
			fprintf(ioQQQ," The number of zones has been increased beyond what it was on the first calculation.\n");
			fprintf(ioQQQ," This can\'t be done since space has already been allocated.\n");
			fprintf(ioQQQ," Have the first calculation do the largest number of zones so that an increase is not needed.\n");
			fprintf(ioQQQ," Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
		return;
	}

	lgFirstCall = false; /* do not reevaluate again */

	/* make sure that molecules have been initialized - this will fail
	 * if this routine is called before size of molecular network is known */
	if( !mole_global.num_total )
	{
		/* mole_global.num_comole_calc can't be zero */
		TotalInsanity();
	}

	init_struc();

	/* allocate space for some arrays used by dynamics routines, and zero out vars */
	DynaCreateArrays( );

	/*************************************************************
	 *                                                           *
	 * get the level 2 line, opacity project, data set           *
	 *                                                           *
	 *************************************************************/

	/* nWindLine is initialized to the dimension of the vector when it is
	 * initialized in the definition at the start of this file.  
	 * it is set to -1 with the "no level2" command, which
	 * stops us from trying to establish this vector */
	if( nWindLine > 0 )
	{
		read_level2_lines();
	}

	/* the UTA line sets - trace will print summary of various data sources */
	if( atmdat.lgInnerShellLine_on )
		read_UTA_lines();

	/* read in data for the set of hyperfine structure lines, and allocate
	 * space for the transition HFLines[nHFLines] structure */
	HyperfineCreate();

	/* Make sure that if hybrid/chianti is on, then Chianti is also on */
	if( atmdat.lgChiantiLvl2Hybrid && !atmdat.lgChiantiOn)
	{
		TotalInsanity();
	}
	if( atmdat.lgStoutLvl2Hybrid && !atmdat.lgStoutOn )
	{
		TotalInsanity();
	}

	/* read in atomic and molecular models from third-party databases */
	if( atmdat.lgLamdaOn || atmdat.lgChiantiOn || atmdat.lgStoutOn)
		database_readin();
	else
		nSpecies = 0;

	/* initialize the large block of level 1 real lines, and OP level 2 lines */
	lines_setup();

	/* mewe_gbar.dat mewe_gbar.dat mewe_gbar.dat mewe_gbar.dat mewe_gbar.dat mewe_gbar.dat ========*/
	/* read in g-bar data taken from
	 *>>refer	all	gbar	Mewe, R., Gronenschild, E. H. B. M., van den Oord, G. H. J. 1985, A&AS, 62, 197 */
	/* open file with Mewe coefficients */

	read_mewe_gbar();

	 /* This is what remains of the t_yield initialization
	  * this should not be in the constructor of t_yield
	  * since it initializes a different struct */
	for( long nelem=0; nelem < LIMELM; nelem++ )
		for( long ion=0; ion < LIMELM; ion++ )
			Heavy.nsShells[nelem][ion] = LONG_MAX;

	/* now read in all auger yields
	 * will do elements from li on up,
	 * skip any line starting with #
	 * this loop goes from lithium to Zn */
	for( long nelem=2; nelem < LIMELM; nelem++ )
	{
		/* nelem is on the shifted C scale, so 2 is Li */
		for( long ion=0; ion <= nelem; ion++ )
		{
			/* number of bound electrons, = atomic number for neutral */
			long nelec = nelem - ion + 1;
			/* one of dima's routines to determine the number of electrons
			 * for this species, nelem +1 to shift to physical number */
			/* subroutine atmdat_outer_shell(iz,in,imax,ig0,ig1)
			 * iz - atomic number from 1 to 30 (integer) 
			 * in - number of electrons from 1 to iz (integer)
			 * Output: imax - number of the outer shell
			 */
			long imax,
				ig0,
				ig1;
			atmdat_outer_shell(nelem+1,nelec,&imax,&ig0,&ig1);

			ASSERT( imax > 0 && imax <= 10 );

			/* nsShells[nelem][ion] is outer shell number for ion with nelec electrons
			 * on physics scale, with K shell being 1 */
			Heavy.nsShells[nelem][ion] = imax;
		}
	}

	/*************************************************************
	 *                                                           *
	 * get the Auger electron yield data set                     *
	 *                                                           *
	 *************************************************************/

	t_yield::Inst().init_yield();

	/****************************************************************
	 *                                                              *
	 * get the Hummer and Storey model case A, B results, these are *
	 * the two data files e1b.dat and e2b.dat, for H and He         *
	 *                                                              *
	 ****************************************************************/

	read_Hummer_Storey();

	// read cross sections for neutral helium
	read_SH98_He1_cross_sections();
	// read cross sections for some he-like ions
	read_Helike_cross_sections();

	// set up spline coefficients for two-photon continua
	atmdat_2phot_setSplineCoefs();

	// get accurate ionization potentials
	read_ionization_potentials();

	return;
}

STATIC void init_struc()
{
	DEBUG_ENTRY( "init_struc()" );

	/* create space for the structure variables */
	/* nzlim will be limit, and number allocated */
	/* >>chng 01 jul 28, define this var, do all following allocations */
	for( long j=0; j < iterations.iter_alloc; j++ )
	{
		struc.nzlim = MAX2( struc.nzlim , iterations.nend[j] );
	}

	/* sloppy, but add one extra for safety */
	++struc.nzlim;

	struc.coolstr.resize(struc.nzlim);
	struc.heatstr.resize(struc.nzlim);
	struc.testr.resize(struc.nzlim);
	struc.volstr.resize(struc.nzlim);
	struc.drad_x_fillfac.resize(struc.nzlim);
	struc.histr.resize(struc.nzlim);
	struc.hiistr.resize(struc.nzlim);
	struc.ednstr.resize(struc.nzlim);
	struc.o3str.resize(struc.nzlim);
	struc.pressure.resize(struc.nzlim);
	struc.windv.resize(struc.nzlim);
	struc.AccelTotalOutward.resize(struc.nzlim);
	struc.AccelGravity.resize(struc.nzlim);
	struc.pres_radiation_lines_curr.resize(struc.nzlim);
	struc.GasPressure.resize(struc.nzlim);
	struc.hden.resize(struc.nzlim);
	struc.DenParticles.resize(struc.nzlim);
	struc.DenMass.resize(struc.nzlim);
	struc.drad.resize(struc.nzlim);
	struc.depth.resize(struc.nzlim);
	struc.depth_last.resize(struc.nzlim);
	struc.drad_last.resize(struc.nzlim);
	struc.xLyman_depth.resize(struc.nzlim);

	if (mole_global.num_calc != 0)
		struc.molecules.alloc(struc.nzlim, mole_global.num_calc);

	struc.H2_abund.resize(struc.nzlim);

	/* create space for gas phase abundances array, first create space for the elements */
	struc.gas_phase.alloc(struc.nzlim, LIMELM);

	/* create space for struc.xIonDense array, first create space for the zones */
	struc.xIonDense.alloc(struc.nzlim, LIMELM, LIMELM+1);

	struc.StatesElem.reserve(struc.nzlim);
	for( long nz=0; nz<struc.nzlim; ++nz)
	{	
		struc.StatesElem.reserve(nz, LIMELM);
		for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				struc.StatesElem.reserve(nz, nelem, nelem+1);
				for( long ion=0; ion < nelem+1; ion++ )
				{
					long ipISO = nelem-ion;
					if( ipISO < NISO )
					{
						struc.StatesElem.reserve(nz, nelem, ion, iso_sp[ipISO][nelem].numLevels_max);
					}
					else
					{
						fixit("for now, set non-iso ions size to 0");
					}
				}
			}
		}
	}
	struc.StatesElem.alloc();

	/* some structure variables */
	for( long i=0; i < struc.nzlim; i++ )
	{
		struc.testr[i] = 0_r;
		struc.volstr[i] = 0_r;
		struc.drad_x_fillfac[i] = 0_r;
		struc.histr[i] = 0_r;
		struc.hiistr[i] = 0_r;
		struc.ednstr[i] = 0_r;
		struc.o3str[i] = 0_r;
		struc.heatstr[i] = 0.;
		struc.coolstr[i] = 0.;
		struc.pressure[i] = 0_r;
		struc.pres_radiation_lines_curr[i] = 0_r;
		struc.GasPressure[i] = 0_r;
		struc.DenParticles[i] = 0_r;
		struc.depth[i] = 0_r;
		struc.H2_abund[i] = 0_r;
	}
	struc.molecules = 0_r;
}

STATIC void read_level2_lines()
{
	DEBUG_ENTRY( "read_level2_lines()" );

	/* begin level2 level 2 wind block */
	/* open file with level 2 line data */

	/* create the TauLine2 emline array */
	TauLine2.resize(nWindLine);
	AllTransitions.push_back(TauLine2);
	cs1_flag_lev2.resize(nWindLine);

	/* first initialize entire array to dangerously large negative numbers */
	for( long i=0; i< nWindLine; ++i )
	{
		/* >>chng 99 jul 16, from setting all t[] to flt_max, to call for
		 * following, each member of structure set to own type of impossible value */
		TauLine2[i].Junk();

		TauLine2[i].AddHiState();
		TauLine2[i].AddLoState();
		TauLine2[i].AddLine2Stack();
	}

	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_readin reading level2.dat\n");

	FILE *ioDATA = open_data( "level2.dat", "r" );

	validate_magic_number_3arg( "level2.dat", ioDATA, 9, 11, 18 );

	/* now get the actual data */
	string chLine;
	for( long i=0; i < nWindLine; ++i )
	{
		do
		{
			if( !read_whole_line( chLine, ioDATA ) )
			{
				fprintf( ioQQQ, " level2.dat error getting line  %li\n", i );
				cdEXIT(EXIT_FAILURE);
			}
			/* skip any line starting with # */
		} while ( chLine[0] == '#' );

		/* this must be double for sscanf to work below */
		double tt[7] = { 0 };
		/*printf("string is %s\n",chLine.c_str() );*/
		sscanf( chLine.c_str() , "%lf %lf %lf %lf %lf %lf %lf " , 
				  &tt[0] ,
				  &tt[1] ,
				  &tt[2] ,
				  &tt[3] ,
				  &tt[4] ,
				  &tt[5] ,
				  &tt[6] );
		/* these are readjusted into their final form in the structure 
		 * in routine lines_setup*/
		(*TauLine2[i].Hi()).nelem() = (int)tt[0];
		(*TauLine2[i].Hi()).IonStg() = (int)tt[1];
		(*TauLine2[i].Lo()).g() = (realnum)tt[2];
		(*TauLine2[i].Hi()).g() = (realnum)tt[3];
		TauLine2[i].Emis().gf() = (realnum)tt[4];
		TauLine2[i].EnergyWN() = (realnum)tt[5];
		cs1_flag_lev2[i] = (realnum)tt[6];
	}

	/* get magic number off last line */
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " level2.dat error getting last magic number\n" );
		cdEXIT(EXIT_FAILURE);
	}
	long magic;
	sscanf( chLine.c_str() , "%ld" , &magic );
	if( 999 != magic )
	{
		fprintf( ioQQQ, " level2.dat ends will wrong magic number=%ld \n", 
		  magic );
		cdEXIT(EXIT_FAILURE);
	}
	fclose( ioDATA );
	if( trace.lgTrace )
		fprintf( ioQQQ," reading level2.dat OK\n");
}

STATIC void read_mewe_gbar()
{
	DEBUG_ENTRY( "read_mewe_gbar()" );

	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_readin reading mewe_gbar.dat\n");

	FILE *ioDATA = open_data( "mewe_gbar.dat", "r" );

	validate_magic_number_1arg( "mewe_gbar.dat", ioDATA, 9101 );

	string chLine;

	/* now get the actual data, indices are correct for c, in Fort went from 2 to 210 */
	for( long i=1; i < 210; i++ )
	{
		if( !read_whole_line( chLine, ioDATA ) )
		{
			fprintf( ioQQQ, " mewe_gbar.dat error getting line  %li\n", i );
			cdEXIT(EXIT_FAILURE);
		}
		/*printf("%s\n",chLine.c_str());*/
		double help[4];
		sscanf( chLine.c_str(), "%lf %lf %lf %lf ", &help[0], &help[1], &help[2], &help[3] );
		for( int l=0; l < 4; ++l )
			MeweCoef.g[i][l] = (realnum)help[l];
	}

	/* get magic number off last line */
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " mewe_gbar.dat error getting last magic number\n" );
		cdEXIT(EXIT_FAILURE);
	}

	long magic;
	sscanf( chLine.c_str() , "%ld" , &magic );

	if( magic != 9101 )
	{
		fprintf( ioQQQ, " mewe_gbar.dat ends will wrong magic number=%ld \n", 
		  magic );
		cdEXIT(EXIT_FAILURE);
	}

	fclose( ioDATA );

	if( trace.lgTrace )
		fprintf( ioQQQ," reading mewe_gbar.dat OK \n");
}

STATIC void read_Hummer_Storey()
{
	DEBUG_ENTRY( "read_Hummer_Storey()" );

	/* now get Hummer and Storey case b data, loop over H, He */
	/* >>chng 01 aug 08, generalized this to both case A and B, and all elements
	 * up to HS_NZ, now 8 */
	for( long ipZ=0; ipZ<HS_NZ; ++ipZ )
	{
		/* don't do the minor elements, Li-B */
		if( ipZ>1 && ipZ<5 ) continue;

		for( long iCase=0; iCase<2; ++iCase )
		{
			/* open Case A or B data */
			/* >>chng 01 aug 08, add HS_ to start of file names to indicate Hummer Storey origin */
			/* create file name for this charge
			 * first character after e is charge number for this element,
			 * then follows whether this is case A or case B data */
			char chFilename[FILENAME_PATH_LENGTH_2] = { 0 };

			sprintf( chFilename, "HS_e%ld%c.dat", ipZ+1, ( iCase == 0 ) ? 'a' : 'b' );

			if( trace.lgTrace )
				fprintf( ioQQQ," atmdat_readin reading Hummer Storey emission file %s\n",chFilename );

			/* open the file */
			FILE *ioDATA = open_data( chFilename, "r" );

			/* read in the number of temperatures and densities*/
			{
				int nread = fscanf( ioDATA, "%li %li ", 
				  &atmdat.ntemp[iCase][ipZ], &atmdat.nDensity[iCase][ipZ] );

				if (nread != 2)
				{
					fprintf(ioQQQ, "atmdat_readin: bad input file format\n");
					cdEXIT(EXIT_FAILURE);
				}
			}

			/* check that ntemp and nDensity are below NHSDIM, 
			 * set to 15 in atmdat_HS_caseB.h */
			assert (atmdat.ntemp[iCase][ipZ] >0 && atmdat.ntemp[iCase][ipZ] <= NHSDIM );
			assert (atmdat.nDensity[iCase][ipZ] > 0 && atmdat.nDensity[iCase][ipZ] <= NHSDIM);

			/* loop reading in line emissivities for all temperatures*/
			for( long ipTemp=0; ipTemp < atmdat.ntemp[iCase][ipZ]; ipTemp++ )
			{
				for( long ipDens=0; ipDens < atmdat.nDensity[iCase][ipZ]; ipDens++ )
				{
					char cha;
					long int junk, junk2;
					{
					   int nread = fscanf( ioDATA, " %lf %li %lf %c %li %ld ", 
						   &atmdat.Density[iCase][ipZ][ipDens], &junk , 
						   &atmdat.ElecTemp[iCase][ipZ][ipTemp], &cha , &junk2 , 
						   /* highest principal quantum number in table, usually 25 */
						   &atmdat.ncut[iCase][ipZ] );

						if (nread != 6)
						{
							fprintf(ioQQQ, "atmdat_readin: bad input file format\n");
							cdEXIT(EXIT_FAILURE);
						}
					}
						
					long ne = atmdat.ncut[iCase][ipZ]*(atmdat.ncut[iCase][ipZ] - 1)/2;
					ASSERT( ne<=NLINEHS );
					for( long j=0; j < ne; j++ )
					{
					   int nread = fscanf( ioDATA, "%lf ",
							    &atmdat.Emiss[iCase][ipZ][ipTemp][ipDens][j] );
						
						if (nread != 1)
						{
							fprintf(ioQQQ, "atmdat_readin: bad input file format\n");
							cdEXIT(EXIT_FAILURE);
						}
					}
				}
			}

			/*this is end of read-in loop */
			fclose(ioDATA); 
			if( trace.lgTrace )
				fprintf( ioQQQ," reading %s OK\n", chFilename );

#			if 0
			/* print list of densities and temperatures */
			for( ipDens=0; ipDens<atmdat.nDensity[iCase][ipZ]; ipDens++ )
			{
				fprintf(ioQQQ," %e,", atmdat.Density[iCase][ipZ][ipDens]);
			}
			fprintf(ioQQQ,"\n");
			for( ipTemp=0; ipTemp<atmdat.ntemp[iCase][ipZ]; ipTemp++ )
			{
				fprintf(ioQQQ," %e,", atmdat.ElecTemp[iCase][ipZ][ipTemp]);
			}
			fprintf(ioQQQ,"\n");
#			endif
		}
	}
}

STATIC void read_SH98_He1_cross_sections(void)
{
	DEBUG_ENTRY( "read_SH98_He1_cross_sections()" );

	bool lgEOL;

	const int ipNUM_FILES = 10;

	char chFileNames[ipNUM_FILES][10] = {
		"p0202.3se",
		"p0202.3po",
		"p0202.3ge",
		"p0202.3fo",
		"p0202.3de",
		"p0202.1se",
		"p0202.1po",
		"p0202.1ge",
		"p0202.1fo",
		"p0202.1de" };

	HS_He1_Xsectn.reserve(26);
	// use real quantum numbers rather than starting at 0
	for( long in=1; in <= 25; in++ )
	{
		// allocate n values of angular momentum, but not more than 5
		HS_He1_Xsectn.reserve(in, min(5,in));
		for( long il=0; il < MIN2(5,in); il++ )
		{
			// allocate two values of spin
			HS_He1_Xsectn.reserve(in, il, 2);
			HS_He1_Xsectn.reserve(in, il, 0, NUM_HS98_DATA_POINTS);
			HS_He1_Xsectn.reserve(in, il, 1, NUM_HS98_DATA_POINTS);
		}
	}
	HS_He1_Xsectn.alloc();
	HS_He1_Energy.alloc(HS_He1_Xsectn.clone());

	string chDirectory = "sh98_he1" + cpu.i().chDirSeparator() + "pi";

	for( long ipFile=0; ipFile<ipNUM_FILES; ipFile++ )
	{
		long S, L, index, N=0;
		long UNUSED P;

		string chPath = chDirectory + cpu.i().chDirSeparator() + chFileNames[ipFile];
		FILE *ioDATA = open_data( chPath, "r" );

		string chLine;
		while( read_whole_line( chLine, ioDATA ) )
		{
			long i=1, s;
			long i1, i2, i3;
			long numDataPoints;
			
			// first line (read above in while) is not needed except that "0 0 0" marks EOF
			i1 = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			i2 = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			i3 = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			if( i1==0 && i2==0 && i3==0 )
				break;

			// don't need next two lines in each set
			read_whole_line( chLine, ioDATA );
			read_whole_line( chLine, ioDATA );
			
			i=1;
			// 4th line in each set identifies the quantum level
			read_whole_line( chLine, ioDATA );
			S = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			L = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			P = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			index = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );

			//indices start with unity
			ASSERT( index >= 1 );

			// index is energy order in that series and is therefore related
			// to principal quantum number.  For triplet S must add one because
			// series starts at n=2. 
			if( L==0 && S==3 )
				N = L + index + 1;
			else
				N = L + index;

			// data go up to n=25
			ASSERT( N<=25 );

			if( S==1 )
				s=0;
			else if( S==3 )
				s=1;
			else
				TotalInsanity();

			i=1;
			// 5th line in each set contains the number of energies
			read_whole_line( chLine, ioDATA );
			//first value is not needed
			FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			numDataPoints = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			// each set has exactly 811 data points, might as well assert this
			if( numDataPoints != NUM_HS98_DATA_POINTS )
			{
				fprintf( ioQQQ, "wrong number of data points!\n" );
				cdEXIT(EXIT_FAILURE);
			}

			// don't need 6th line either
			read_whole_line( chLine, ioDATA );

			// now begin reading in lines
			// there must be exactly numDataPoints ordered pairs,
			// throw an assert for any deviation
			for( long k=0; k<NUM_HS98_DATA_POINTS; k++ )
			{
				i=1;
				read_whole_line( chLine, ioDATA );
				HS_He1_Energy[N][L][s][k] = (double)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
				HS_He1_Xsectn[N][L][s][k] = (double)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			}
		}

		// we reached end of file, assert last quantum number was 25
		ASSERT( N==25 );

		fclose( ioDATA );
	}

	return;
}

STATIC void read_Helike_cross_sections(void)
{
	DEBUG_ENTRY( "read_Helike_cross_sections()" );

	bool lgEOL;

	const char chFileName[23] = "helike_pcs_topbase.dat";

	// the data only go as high as n=10	
	const int MaxN = 10;
	long last_i1=0;

	// data will be used up to calcium (nelem=19)
	// data exists up to iron but we will not use it because it has strong resonance
	// features, but is nonetheless very nearly hydrogenic
	OP_Helike_Xsectn.reserve(ipCALCIUM+1);
	for( long nelem = ipLITHIUM; nelem <= ipCALCIUM; nelem++ )
	{
		// allocate principal quantum number
		OP_Helike_Xsectn.reserve(nelem, MaxN+1);
		for( long in = 1; in <= MaxN; in++ )
		{
			// allocate angular momentum
			OP_Helike_Xsectn.reserve(nelem, in, in);
			for( long il = 0; il < in; il++ )
			{
				// allocate two values of spin
				OP_Helike_Xsectn.reserve(nelem, in, il, 2);
			}
		}
	}
	OP_Helike_Xsectn.alloc();
	OP_Helike_Energy.alloc(OP_Helike_Xsectn.clone());
	OP_Helike_NumPts.alloc(OP_Helike_Xsectn.clone());

	OP_Helike_NumPts = 0;

	FILE *ioDATA = open_data( chFileName, "r" );

	// Header looks like this:
	// ================================================
	//       I  NZ  NE  ISLP  ILV        E(RYD)      NP
	// ================================================
	//      1   3   2   100    1  -5.53159E+00      55

	// so skip the first three lines.
	string chLine;
	for( long i=0; i<3; i++)
	{
		if( !read_whole_line( chLine, ioDATA ) )
		{
			fprintf( ioQQQ,"PROBLEM corruption in TOPbase Helike pcs datafile.\nSorry\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	while( read_whole_line( chLine, ioDATA ) )
	{
		long i=1;
		long n, l, s;
		long i1, i2, i3, i4, i5, i7;
		double i6;
		
		i1 = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
		i2 = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
		i3 = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
		i4 = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
		i5 = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
		i6 = (double)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
		i7 = (long)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );

		if( lgEOL )
		{
			fprintf( ioQQQ,"PROBLEM corruption in TOPbase Helike pcs datafile.\nSorry\n" );
			cdEXIT(EXIT_FAILURE);
		}

		// this marks end of data.
		if( i1==i2 && i1==i3 && i1==i4 && i1==i5 && i1==i7 && i1==-1 )
			break;

		// first parameter is level index (overall, not just for series or charge)
		// check that it is as expected
		if( ! ( i1 > 0 && i1 == (last_i1 + 1) && i1 <= 795 ) )
		{
			fprintf( ioQQQ, "syntax error found in %s\n", chFileName );
			cdEXIT(EXIT_FAILURE);
		}
		last_i1 = i1;
		// second parameter is nuclear charge
		ASSERT( (i2-1)<=ipCALCIUM && (i2-1)>=ipLITHIUM );
		// third parameter must be 2 for helike.
		ASSERT( i3==2 );
		// fourth parameter is (2S+1)*100+L*10+P 
		ASSERT( i4>=100 && i4<400 );
		if( i4 >= 300 )
			s=1;
		else
			s=0;
		l = (i4 - (2*s+1)*100)/10;
		// data only goes up to l=2
		ASSERT( l<=2 );
		// fifth is index in the series, related to principal quantum number
		ASSERT( i5>=1 && i5<=10 );
		if( s==1 && l==0 )
			n = i5 + 1;
		else
			n = i5 + l;
		ASSERT( l<=MaxN );
		// sixth is threshhold energy, don't need but assert negative
		// \todo 3 save this and renorm cross-section with ratio of actual to recorded Eth?
		if( i6 >= 0. )
		{
			fprintf( ioQQQ, "invalid threshold energy in %s\n", chFileName );
			cdEXIT(EXIT_FAILURE);
		}
		// seventh parameter is number of data points, can be zero, use to allocate otherwise
		OP_Helike_NumPts[i2-1][n][l][s] = i7;
		if( i7 == 0 )
			continue;
		
		OP_Helike_Xsectn[i2-1][n][l][s].resize(i7);
		OP_Helike_Energy[i2-1][n][l][s].resize(i7);

		// now begin reading in lines
		// there must be exactly i7 ordered pairs,
		for( long k=0; k < i7; k++ )
		{
			i=1;
			read_whole_line( chLine, ioDATA );
			OP_Helike_Energy[i2-1][n][l][s][k] = (double)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			OP_Helike_Xsectn[i2-1][n][l][s][k] = (double)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );

			// make sure data is well-behaved
			if( k > 0 )
			{
				ASSERT( OP_Helike_Energy[i2-1][n][l][s][k] > OP_Helike_Energy[i2-1][n][l][s][k-1] );
				ASSERT( OP_Helike_Xsectn[i2-1][n][l][s][k] < OP_Helike_Xsectn[i2-1][n][l][s][k-1] );
			}

			// try to read one more item off the line and verify that lgEOL is true
			FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			ASSERT( lgEOL );
		}
	}

	fclose( ioDATA );

	return;
}

t_yield::t_yield()
{
	 /* preset two arrays that will hold auger data 
	  * set to very large values so that code will blow
	  * if not set properly below */
	for( int nelem=0; nelem < LIMELM; nelem++ )
	{
		for( int ion=0; ion < LIMELM; ion++ )
		{
			for( int ns=0; ns < 7; ns++ )
			{
				n_elec_eject[nelem][ion][ns] = LONG_MAX;
				for( int nelec=0; nelec < 10; nelec++ )
				{
					frac_elec_eject[nelem][ion][ns][nelec] = FLT_MAX;
				}
			}
		}
	}

	lgKillAuger = false;
}

void t_yield::init_yield()
{
	/* following is double for sscanf to work */
	double temp[14];

	DEBUG_ENTRY( "init_yield()" );

	/*************************************************************
	 *                                                           *
	 * get the Auger electron yield data set                     *
	 *                                                           *
	 *************************************************************/

	/* NB NB -- This test of Heavy.nsShells remains here to assure
	 * that it contains meaningful values since needed below, once
	 * t_Heavy is a Singleton, it can be removed !!! */
	ASSERT( Heavy.nsShells[2][0] > 0 );

	/* hydrogen and helium will not be done below, so set yields here*/
	n_elec_eject[ipHYDROGEN][0][0] = 1;
	n_elec_eject[ipHELIUM][0][0] = 1;
	n_elec_eject[ipHELIUM][1][0] = 1;

	frac_elec_eject[ipHYDROGEN][0][0][0] = 1;
	frac_elec_eject[ipHELIUM][0][0][0] = 1;
	frac_elec_eject[ipHELIUM][1][0][0] = 1;

	/* open file auger.dat, yield data file that came from
	 * >>refer	all	auger	Kaastra, J. S., & Mewe, R. 1993, A&AS, 97, 443-482 */
	const char*chFilename = "mewe_nelectron.dat";

	if( trace.lgTrace )
		fprintf( ioQQQ, " init_yield reading %s\n", chFilename );

	/* open the file */
	FILE *ioDATA = open_data( chFilename, "r" );

	/* now read in all auger yields
	 * will do elements from li on up,
	 * skip any line starting with #
	 * this loop goes from lithium to Zn */
	string chLine;
	for( int nelem=2; nelem < LIMELM; nelem++ )
	{
		/* nelem is on the shifted C scale, so 2 is Li */
		for( int ion=0; ion <= nelem; ion++ )
		{
			for( int ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
			{
				char ch1 = '#';
				/* the * is the old comment char, accept it, but really want # */
				while( ch1 == '#' || ch1 == '*' )
				{
					if( !read_whole_line( chLine, ioDATA ) )
					{
						fprintf( ioQQQ, " %s error getting line %i\n", chFilename, ns );
						cdEXIT(EXIT_FAILURE);
					}
					ch1 = chLine[0];
				}
				
				if( sscanf( chLine.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
					    &temp[0], &temp[1], &temp[2], &temp[3], &temp[4],
					    &temp[5], &temp[6], &temp[7], &temp[8], &temp[9],
					    &temp[10],&temp[11],&temp[12],&temp[13] ) != 14 )
				{
					fprintf( ioQQQ, "failed to read correct number of arguments in %s\n",
						 chFilename );
					cdEXIT(EXIT_FAILURE);
				}
				n_elec_eject[nelem][ion][ns] = (long int)temp[3];

				ASSERT( n_elec_eject[nelem][ion][ns] >= 0 && n_elec_eject[nelem][ion][ns] < 11 );

				/* can pop off up to 10 auger electrons, these are the probabilities*/
				for( int j=0; j < 10; j++ )
				{
					frac_elec_eject[nelem][ion][ns][j] = (realnum)temp[j+4];
					ASSERT( frac_elec_eject[nelem][ion][ns][j] >= 0. );
				}
			}
		}
		/* activate this print statement to get yield for k-shell of all atoms */
		/*fprintf(ioQQQ,"yyyield\t%li\t%.4e\n", nelem+1, frac_elec_eject[nelem][0][0][0] );*/
	}

	fclose( ioDATA );

	if( trace.lgTrace )
	{
		/* this is set with no auger command */
		if( lgKillAuger )
			fprintf( ioQQQ, " Auger yields will be killed.\n");
		fprintf( ioQQQ, " reading %s OK\n", chFilename );
	}

	/* open file mewe_fluor.dat, yield data file that came from
	 * >>refer	all	auger	Kaastra, J. S., & Mewe, R. 1993, 97, 443-482 */
	chFilename = "mewe_fluor.dat";

	if( trace.lgTrace )
		fprintf( ioQQQ, " init_yield reading %s\n", chFilename );

	/* open the file */
	ioDATA = open_data( chFilename, "r" );

	/* now read in all auger yields
	 * will do elements from li on up,
	 * skip any line starting with # */
	do 
	{
		if( !read_whole_line( chLine, ioDATA ) )
		{
			fprintf( ioQQQ, " %s error getting line %i\n", chFilename, 0 );
			cdEXIT(EXIT_FAILURE);
		}
	}
	while( chLine[0] == '#' );

	bool lgEOL = false;

	nfl_lines = 0;
	do
	{
		const int NKM = 10;
		int nDima[NKM] = { 0, 1, 2, 2, 3, 4, 4, 5, 5, 6 };
		int nAuger;

		if( nfl_lines >= MEWE_FLUOR )
			TotalInsanity();

		/*printf("string is %s\n",chLine.c_str() );*/
		sscanf( chLine.c_str(), "%lf %lf %lf %lf %lf %lf %lf", 
			&temp[0], &temp[1], &temp[2], &temp[3], &temp[4], 
			&temp[5], &temp[6] );

		/* the atomic number, C is 5 */
		nfl_nelem[nfl_lines] = (int)temp[0]-1;
		ASSERT( nfl_nelem[nfl_lines] >= 0 && nfl_nelem[nfl_lines] < LIMELM );

		/* the ion stage for target, atom is 0 */
		nfl_ion[nfl_lines] = (int)temp[1]-1;
		ASSERT( nfl_ion[nfl_lines] >= 0 && nfl_ion[nfl_lines] <= nfl_nelem[nfl_lines]+1 );

		/* the target's shell */
		nfl_nshell[nfl_lines] = nDima[(long)temp[2]-1];
		ASSERT( nfl_nshell[nfl_lines] >= 0 && 
			/* nsShells is shell number, where K is 1 */
			nfl_nshell[nfl_lines] < Heavy.nsShells[nfl_nelem[nfl_lines]][nfl_ion[nfl_lines]]-1 );

		/* this is the number of Auger electrons ejected */
		nAuger = (int)temp[3];
		/* so this is the spectrum of the photons that are emitted */
		nfl_ion_emit[nfl_lines] = nfl_ion[nfl_lines] + nAuger + 1;
		/* must be gt 0 since at least photoelectron comes off */
		ASSERT( nfl_ion_emit[nfl_lines] > 0 && nfl_ion_emit[nfl_lines] <= nfl_nelem[nfl_lines]+1);

		/* this is the type of line as defined in their paper */
		nfl_nLine[nfl_lines] = (int)temp[4];

		/* energy in Ryd */
		fl_energy[nfl_lines] = (realnum)temp[5] / (realnum)EVRYD;
		ASSERT( fl_energy[nfl_lines] > 0. );

		/* fluor yield */
		fl_yield[nfl_lines] = (realnum)temp[6];
		/* NB cannot assert <=1 since data file has yields around 1.3 - 1.4 */
		ASSERT( fl_yield[nfl_lines] >= 0 );

		++nfl_lines;

		do
		{
			if( !read_whole_line( chLine, ioDATA ) )
				lgEOL = true;
		} 
		while( chLine[0]=='#' && !lgEOL );
	} 
	while( !lgEOL );

	fclose( ioDATA );
	if( trace.lgTrace )
		fprintf( ioQQQ, " reading %s OK\n", chFilename );
}

// read autoionization data from Badnell data file 
STATIC void ReadBadnellAIData(const string& fnam,      // filename containing the Badnell data
			      long nelem,              // nelem is on C scale
			      long ion,                // ion is on C scale
			      TransitionList& UTA,     // UTA lines will be pushed on this stack
			      bitset<IS_TOP> Skip)     // option to skip transitions from a particular shell
{
	DEBUG_ENTRY( "ReadBadnellAIData()" );

	if( trace.lgTrace )
		fprintf( ioQQQ," ReadBadnellAIData reading %s\n", fnam.c_str() );

	fstream ioDATA;
	open_data( ioDATA, fnam, mode_r );

	string line;
	getline( ioDATA, line );
	ASSERT( line.substr(0,4) == "SEQ=" );
	getline( ioDATA, line );
	getline( ioDATA, line );
	// we don't need the parent level data, so we will skip it...
	ASSERT( line.substr(3,21) == "PARENT LEVEL INDEXING" );
	int nParent;
	istringstream iss( line.substr(65,4) );
	iss >> nParent;
	// data lines containing data for all levels of the parent ion will span nMulti lines
	int nMulti = (nParent+5)/6;
	for( int i=0; i < nParent+5; ++i )
		getline( ioDATA, line );

	// here starts the header for the level data we need
	ASSERT( line.substr(3,26) == "IC RESOLVED LEVEL INDEXING" );
	int nLevel;
	istringstream iss2( line.substr(63,6) );
	iss2 >> nLevel;
	// skip rest of header
	for( int i=0; i < 3; ++i )
		getline( ioDATA, line );

	// now get the level data
	vector<t_BadnellLevel> level( nLevel );
	for( int i=0; i < nLevel; ++i )
	{
		getline( ioDATA, line );
		istringstream iss3( line );
		int indx, irsl;
		iss3 >> indx >> irsl;
		level[indx-1].irsl = irsl;
		level[indx-1].config = line.substr(16,20);
		istringstream iss4( line.substr(37,1) );
		iss4 >> level[indx-1].S;
		istringstream iss5( line.substr(39,1) );
		iss5 >> level[indx-1].L;
		istringstream iss6( line.substr(41,4) );
		double J;
		iss6 >> J;
		level[indx-1].g = nint(2.*J + 1.);
		istringstream iss7( line.substr(46,11) );
		iss7 >> level[indx-1].energy;

		// which inner shell has been excited?
		level[indx-1].lgAutoIonizing = ( line[57] == '*' );
		if( level[indx-1].lgAutoIonizing )
		{
			if( level[indx-1].config.find( "1S1" ) != string::npos )
				level[indx-1].WhichShell = IS_K_SHELL;
			else if( level[indx-1].config.find( "2S1" ) != string::npos )
				level[indx-1].WhichShell = IS_L1_SHELL;
			else if( level[indx-1].config.find( "2P5" ) != string::npos )
				level[indx-1].WhichShell = IS_L2_SHELL;
			else
				TotalInsanity();
		}
		else
		{
			level[indx-1].WhichShell = IS_NONE;
		}
	}

	// levels are done, now move on to the lines
	// first search for start of the header
	while( getline( ioDATA, line ) )
	{
		if( line.find( "IRSL  IRSL" ) != string::npos )
			break;
	}
	// skip rest of the header
	for( int i=0; i < nMulti-1; ++i )
		getline( ioDATA, line );

	// blank line that will be pushed on the UTA line stack
	qList BlankStates("BlankStates",1);
	TransitionList BlankList("BlankList",&BlankStates,1);
	TransitionList::iterator BlankLine = BlankList.begin();
	(*BlankLine).Junk();

	// start reading the line data
	while( getline( ioDATA, line ) )
	{
		// have we reached the end of the line data?
		if( line.size() < 10 )
			break;

		// test if there is an autoionization rate on this line
		// if not, it only contains radiative data and we skip it...
		if( line.size() < 50 )
			continue;

		// there may be an asterisk here; wipe it out since we don't need it
		line[19] = ' ';

		int irsl_lo, irsl_hi, dum;
		double edif, Bij, Rji, Aai;
		istringstream iss8( line );
		// irsl_lo: index for lower level of transition
		// irsl_hi: index for upper level of transition
		// edif: energy difference between levels in Ryd
		// Bij: UPWARD Einstein A for transition irsl_lo -> irsl_hi
		// Rji: sum of Aji for all radiative transitions to lower levels
		// Aai: autoionization rate from level irsl_hi
		iss8 >> irsl_lo >> irsl_hi >> dum >> dum >> edif >> Bij >> Rji >> Aai;
		// ind_lo and ind_hi are on C scale
		int ind_lo = irsl2ind( level, irsl_lo );
		int ind_hi = irsl2ind( level, irsl_hi );
		ASSERT( level[ind_hi].lgAutoIonizing );
		// skip rest of partial autoionization rates
		for( int i=0; i < nMulti-1; ++i )
			getline( ioDATA, line );
		// skip this transition if it does not originate from ground level
		// or if the user requested to skip excitations from this inner shell
		if( ind_lo == 0 && !Skip[level[ind_hi].WhichShell] )
		{
			UTA.push_back( *BlankLine );
			InitTransition( UTA.back() );

			// t_emission has nelem and ion on fortran scale...
			(*UTA.back().Hi()).nelem() = nelem+1;
			(*UTA.back().Hi()).IonStg() = ion+1;

			(*UTA.back().Hi()).g() = (realnum)level[ind_hi].g;
			(*UTA.back().Lo()).g() = (realnum)level[ind_lo].g;

			double WavNum = edif*RYD_INF;

			/* wavelength in Angstroms */
			UTA.back().WLAng() = (realnum)(1e8/WavNum);
			UTA.back().EnergyWN() = (realnum)WavNum;

			/* store branching ratio for autoionization */
			double frac_ioniz = Aai/(Rji + Aai);
			ASSERT( frac_ioniz >= 0. &&  frac_ioniz <= 1. );
			UTA.back().Emis().AutoIonizFrac() = (realnum)frac_ioniz;

			/* this is true spontaneous rate for doubly excited state to ground
			 * and is used for pumping, and also relaxing back to inner shell */
			/* Badnell gives UPWARD transition rate Alu, an unusual notation,
			 * convert it here to the normal downward transition rate Aul */
			UTA.back().Emis().Aul() = (realnum)(Bij*(*UTA.back().Lo()).g()/(*UTA.back().Hi()).g());
			UTA.back().Emis().gf() =
				(realnum)GetGF( UTA.back().Emis().Aul(), UTA.back().EnergyWN(), (*UTA.back().Hi()).g() );

			UTA.back().Emis().iRedisFun() = ipPRD;

			UTA.back().Emis().dampXvel() = (realnum)((Rji+Aai)/UTA.back().EnergyWN()/PI4);

			ASSERT( UTA.back().Emis().dampXvel() > 0. );

			// remove this line if it is too weak
			if( UTA.back().Emis().gf() < level[ind_lo].g * f_cutoff )
				UTA.pop_back();
		}
	}

	// perform a sanity check on the tail of the file
	getline( ioDATA, line );
	ASSERT( line.substr(3,7) == "NRSLMX=" );

	ioDATA.close();

	if( trace.lgTrace )
		fprintf( ioQQQ, " reading %s OK\n", fnam.c_str() );
}

inline void InitTransition(const TransitionProxy& t)
{
	t.AddHiState();
	t.AddLoState();
	t.AddLine2Stack();
}

inline int irsl2ind(vector<t_BadnellLevel>& level, int irsl)
{
	for( unsigned int i=0; i < level.size(); ++i )
	{
		if( level[i].irsl == irsl )
			return (int)i;
	}
	// we should never get here...
	TotalInsanity();
}

STATIC void validate_magic_number_1arg( const char *chFilename, FILE *ioFile, const long magicExp )
{
	DEBUG_ENTRY( "validate_magic_number_1arg()" );

	string chLine;
	while( read_whole_line( chLine, ioFile ) )
	{
		// skip any #
		if( chLine[0] != '#' )
			break;
	}

	long magic = -1;
	sscanf( chLine.c_str(), "%ld", &magic );

	if( magic != magicExp )
	{
		fprintf( ioQQQ,
			" atmdat_readin: the version of '%s' is not the current version.\n",
			chFilename );
		fprintf( ioQQQ,
			" I expected to find the number %ld and got %ld instead.\n" ,
			magicExp, magic );
		cdEXIT(EXIT_FAILURE);
	}
}

STATIC void validate_magic_number_3arg( const char *chFilename, FILE *ioFile,
				const long yearExp, const long monthExp, const long dayExp )
{
	DEBUG_ENTRY( "validate_magic_number_3arg()" );

	long year,
		month,
		day;

	string chLine;
	while( read_whole_line( chLine, ioFile ) )
	{
		// skip any #
		if( chLine[0] != '#' )
			break;
	}

	sscanf( chLine.c_str(), "%ld %ld %ld", &year, &month, &day );

	if( year != yearExp || month != monthExp || day != dayExp )
	{
		fprintf( ioQQQ,
			" atmdat_readin: the version of '%s' is not the current version.\n",
			chFilename );
		fprintf( ioQQQ,
			" I expected to find the number %ld %ld %ld and got %ld %ld %ld instead.\n" ,
			yearExp, monthExp, dayExp, year, month, day );
		cdEXIT(EXIT_FAILURE);
	}
}

STATIC void read_UTA_lines()
{
	DEBUG_ENTRY( "read_UTA_lines()" );

	/* reserve space for all data sets, no worries if too small though... */
	UTALines.reserve( 4000 );
	AllTransitions.push_back(UTALines);

	// version of element symbols in lower case and without spaces....
	const char* chElmSymLC[] =
		{ "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na", "mg", "al", "si", "p",
		  "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn" };

	// save cite for UTA sources, insure no double counting
	char chUTA_ref[LIMELM][LIMELM][5];
	for( long nelem=0; nelem < LIMELM; ++nelem )
		for( long ion=0; ion <= nelem; ++ion )
			strcpy( chUTA_ref[nelem][ion] , "" );

	/* first read in the Badnell data */
	for( long ipISO=ipLI_LIKE; ipISO < ipAL_LIKE; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			// ion = 0 for neutral atom
			long ion = nelem - ipISO;
			strcpy( chUTA_ref[nelem][ion] , "B" );

			bitset<IS_TOP> Skip;
			if( ipISO < ipNA_LIKE )
			{
				// construct file name
				ostringstream oss;
				// Badnell calls Li-like series He-like, etc...
				oss << "UTA/nrb00_" << chElmSymLC[ipISO-1] << "_";
				// Badnell uses ion = 1 for neutral atom
				oss << chElmSymLC[nelem] << ion+1 << "ic1-2.dat";
				// now read the data...
				ReadBadnellAIData( oss.str(), nelem, ion, UTALines, Skip );
			}
			else
			{
				// from Na-like onwards both K-shell (ic1-3) and L-shell (ic2-3)
				// excitation are treated in two separate files
				ostringstream oss;
				oss << "UTA/nrb00_" << chElmSymLC[ipISO-1] << "_";
				oss << chElmSymLC[nelem] << ion+1 << "ic1-3.dat";
				ReadBadnellAIData( oss.str(), nelem, ion, UTALines, Skip );

				// Kisielius L1, L2-shell data sets take precedence for Na, Mg, and Al-like iron
				if( atmdat.lgInnerShell_Kisielius && nelem == ipIRON &&
				    ipISO >= ipNA_LIKE && ipISO <= ipAL_LIKE )
				{
					Skip[IS_L1_SHELL] = 1;
					Skip[IS_L2_SHELL] = 1;
				}

				ostringstream oss2;
				oss2 << "UTA/nrb00_" << chElmSymLC[ipISO-1] << "_";
				oss2 << chElmSymLC[nelem] << ion+1 << "ic2-3.dat";
				ReadBadnellAIData( oss2.str(), nelem, ion, UTALines, Skip );
			}
		}
	}

	/* these are the statistical weights for the ground levels of all ions of iron */
	const realnum StatWeightGroundLevelIron[] =
		{ 9.f, 10.f, 9.f, 6.f, 1.f, 4.f, 5.f, 4.f, 1.f, 4.f, 5.f, 4.f, 1.f,
		  2.f, 1.f, 2.f, 1.f, 4.f, 5.f, 4.f, 1.f, 2.f, 1.f, 2.f, 1.f, 2.f };

	// blank line that will be pushed on the UTA line stack
	qList BlankStates("BlankStates",1);
	TransitionList BlankList("BlankList",&BlankStates,1);
	TransitionList::iterator BlankLine = BlankList.begin();
	(*BlankLine).Junk();

	/* next read in the Gu file */
	{
		/* read the Gu et al. (2006) data
		 * >>refer	Fe	UTA	Gu, M. F., Holczer T., Behar E., & Kahn S. M. 2006, ApJ 641, 1227-1232 */
		if( trace.lgTrace )
			fprintf( ioQQQ," atmdat_readin reading UTA_Gu06.dat\n");

		FILE *ioGU06 = open_data( "UTA/UTA_Gu06.dat", "r" );

		validate_magic_number_3arg( "UTA_Gu06.dat", ioGU06, 2007, 1, 23 );

		int nelemGu =-1, ionGu=-1;
		string chLine;

		while( read_whole_line( chLine, ioGU06 ) )
		{
			if( chLine[0] != '#' )
			{
				long ion, i2;
				double EnergyAng, Aul, oscill, Aauto;

				sscanf( chLine.c_str(), "%4li%5li%8lf%13lf%12lf",
					&ion, &i2, &EnergyAng, &Aul, &Aauto );
				sscanf( chLine.c_str()+54, "%13lf", &oscill );

				// avoid duplication of ions: anything upto and including the Mg-like
				// series is covered by the Badnell data set, Al-like is covered by
				// the Kisielius data set if it is turned on.
				int ipISO = ipIRON - ion + 1;
				int ipThres = atmdat.lgInnerShell_Kisielius ? ipAL_LIKE : ipMG_LIKE;
				if( ipISO <= ipThres )
					continue;

				UTALines.push_back( *BlankLine );
				InitTransition( UTALines.back() );

				/* all these are iron, first number was ion stage with 0 the atom */
				(*UTALines.back().Hi()).nelem() = ipIRON+1;

				/* now do stage of ionization */
				ASSERT( ion > 0 && ion <= ipIRON );
				/* the ion stage - 1 is atom - this data file has different
				 * format from other two - others are 0 for atom */
				(*UTALines.back().Hi()).IonStg() = ion;
				if( ipIRON!=nelemGu || ion!=ionGu )
				{
					// one label per ion
					nelemGu = ipIRON;
					ionGu = ion;
					strcpy( chUTA_ref[ipIRON][ion-1] , "G" );
				}

				/* these are the statistical weights 
				 * lower levels are not included in the original data file */
				if( strstr_s( chLine.c_str(), "(J=1/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 2.f;
				else if( strstr_s( chLine.c_str(), "(J=1)" ) != NULL )
					(*UTALines.back().Hi()).g() = 3.f;
				else if( strstr_s( chLine.c_str(), "(J=3/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 4.f;
				else if( strstr_s( chLine.c_str(), "(J=2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 5.f;
				else if( strstr_s( chLine.c_str(), "(J=5/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 6.f;
				else if( strstr_s( chLine.c_str(), "(J=3)" ) != NULL )
					(*UTALines.back().Hi()).g() = 7.f;
				else if( strstr_s( chLine.c_str(), "(J=7/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 8.f;
				else if( strstr_s( chLine.c_str(), "(J=4)" ) != NULL )
					(*UTALines.back().Hi()).g() = 9.f;
				else if( strstr_s( chLine.c_str(), "(J=9/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 10.f;
				else if( strstr_s( chLine.c_str(), "(J=5)" ) != NULL )
					(*UTALines.back().Hi()).g() = 11.f;
				else if( strstr_s( chLine.c_str(), "(J=11/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 12.f;
				else
					TotalInsanity();
				(*UTALines.back().Lo()).g() = StatWeightGroundLevelIron[ion-1];

				/* wavelength in Angstroms */
				double fenergyWN = 1e8/EnergyAng;
				UTALines.back().EnergyWN() = fenergyWN;
				UTALines.back().WLAng() = (realnum) wn2ang( fenergyWN );

				/* store branching ratio for autoionization */
				double frac_ioniz = Aauto/(Aul + Aauto);
				ASSERT( frac_ioniz >= 0. &&  frac_ioniz <= 1. );
				UTALines.back().Emis().AutoIonizFrac() = (realnum)frac_ioniz;

				/* save gf scanned from line */
				UTALines.back().Emis().gf() = (*UTALines.back().Lo()).g() * (realnum)oscill;

				/* this is true spontaneous rate for doubly excited state to inner 
				 * shell UTA, and is used for pumping, and also relaxing back to inner shell */
				UTALines.back().Emis().Aul() =
					(realnum)eina( UTALines.back().Emis().gf(),
					UTALines.back().EnergyWN(), (*UTALines.back().Hi()).g() );

				ASSERT( fp_equal_tol( (realnum)Aul, UTALines.back().Emis().Aul(), 1.e-3f*(realnum)Aul ) );

				UTALines.back().Emis().iRedisFun() = ipPRD;

				UTALines.back().Emis().dampXvel() = (realnum)(
						(UTALines.back().Emis().Aul()+Aauto) /
						UTALines.back().EnergyWN()/PI4);
				ASSERT( UTALines.back().Emis().dampXvel()>0. );

				// ignore line if too weak
				if( UTALines.back().Emis().gf() < StatWeightGroundLevelIron[ion-1] * f_cutoff )
					UTALines.pop_back();
			}
		}

		fclose( ioGU06 );

		if( trace.lgTrace )
			fprintf( ioQQQ, " reading UTA_Gu06.dat OK\n" );
	}

	if( atmdat.lgInnerShell_Kisielius )
	{
		/* last read in the Romas Kisielius data
		 *>>refer	Fe	UTA	Kisielius, R., Hibbert, A.. Ferland, G. J., et al. 2003, MNRAS, 344, 696 */
		if( trace.lgTrace )
			fprintf( ioQQQ," atmdat_readin reading UTA_Kisielius.dat\n");

		FILE *ioROMAS = open_data( "UTA/UTA_Kisielius.dat", "r" );

		validate_magic_number_3arg( "UTA_Kisielius.dat", ioROMAS, 11, 8, 25 );

		long int nRomasUsed = 0 , nRomasTotal = 0;
		int nelemRomas=-1 , ionRomas=-1;
		FILE *ioROMASused=NULL;
		bool lgSaveRomasUsed = false;
		if( lgSaveRomasUsed )
		{
			if( (ioROMASused = open_data("RomasUsed.txt","w")) == NULL )
			{
				fprintf(ioQQQ,"could not open RomasUsed.txt\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		string chLine;

		while( read_whole_line( chLine, ioROMAS ) )
		{
			/* only look at lines without '#' in first col */
			if( chLine[0] != '#' )
			{
				long int i1, i2, i3;
				double f1, f2, oscill;
				double frac_relax;

				++nRomasTotal;
				sscanf( chLine.c_str(), "%li\t%li\t%li\t%lf\t%lf\t%lf\t%lf",
					&i1,&i2,&i3,&f1,&f2,&frac_relax,&oscill );

				// skip line if too weak
				// Fe+13 has 25 000 lines, most are vastly weak and do not add to integrated rate
				// the cutoff of 1e-4 reduces the number of Romas UTA lines from 84224 to 573
				if( oscill < f_cutoff )
					continue;

				if( lgSaveRomasUsed )
					fprintf(ioROMASused , "%s" , chLine.c_str());

				/* For Fe XIV both levels of the 2P* ground term are present in the data.
				 * following true, ignore the transitions from the excited level
				 * false, assume ground term populated by statistical weight if
				 * "fudge 0" also appears in input stream */
				const bool lgAllowSplitFe14 = false;
				if( lgAllowSplitFe14 || i2 == StatWeightGroundLevelIron[i1] )
				{
					++nRomasUsed;
					UTALines.push_back( *BlankLine );
					InitTransition( UTALines.back());

					/* all these are iron, first number was ion stage with 0 the atom */
					(*UTALines.back().Hi()).nelem() = ipIRON+1;

					/* now do stage of ionization */
					ASSERT( i1 >= 0 && i1 < ipIRON );
					/* NB - the plus one is because 1 will be subtracted when used,
					 * in the original data file i1 is 0 for the atom */
					(*UTALines.back().Hi()).IonStg() = i1 + 1;

					realnum facpop = 1.;
					// allow low or high density limit for level populations in ground term of Fe XIV
					if( i1 == 13 )
					{
						// Fe XIV - fudge -1 returns number of fudge parameters, 0 if not specified
						if( lgAllowSplitFe14 )
						{
							if( i2 == StatWeightGroundLevelIron[i1] )
								facpop= 0.333;
							else
								facpop = 0.667;
						}
						else
						{
							if( i2 == StatWeightGroundLevelIron[i1] )
								facpop= 1.;
							else
								facpop = 1e-10;
						}
					}
					if( ipIRON!=nelemRomas || i1!=ionRomas )
					{
						// one label per ion
						nelemRomas = ipIRON;
						ionRomas = i1;
						strcpy( chUTA_ref[ipIRON][i1] , "K" );
					}

					/* these were the statistical weights */
					(*UTALines.back().Hi()).g() = (realnum)i3;
					(*UTALines.back().Lo()).g() = (realnum)i2;

					UTALines.back().WLAng() = (realnum)f1;
					UTALines.back().EnergyWN() = 1.e8f/(realnum)f1;
					// print transitions contributing to 15.5A UTA feature
#					if 0
					if( i1==13 && f1>15.35 && f1<15.55)
					{
						fprintf(ioQQQ,"DEBUG %li\t%.5f\t%.3e\n",i2, f1 , oscill * facpop);
					}
#					endif

					/* this is true spontaneous rate for doubly excited state to inner shell,
					 * and is used for pumping, and also relaxing back to inner shell */
					UTALines.back().Emis().gf() = (*UTALines.back().Lo()).g() * (realnum)oscill*facpop;
					UTALines.back().Emis().Aul() =
						(realnum)eina( UTALines.back().Emis().gf(),
							UTALines.back().EnergyWN(),
							(*UTALines.back().Hi()).g() );

					UTALines.back().Emis().iRedisFun() = ipPRD;

					/* store branching ratio for autoionization */
					ASSERT( frac_relax >= 0.f &&  frac_relax <= 1.f );
					UTALines.back().Emis().AutoIonizFrac() = 1.f-(realnum)frac_relax;

					// Romas data do not have autoionization rates so take typical
					// value of 1e15 s-1, suggested by Badnell OP data
					UTALines.back().Emis().dampXvel() = (realnum)(
							(UTALines.back().Emis().Aul()+1e15) /
							UTALines.back().EnergyWN()/PI4);
					ASSERT( UTALines.back().Emis().dampXvel()>0. );
					//fprintf(ioQQQ,"DEBUGGG %li %.3e\n", nRomasUsed, UTALines.back().Emis().dampXvel() );
				}
			}
		}

		fclose( ioROMAS );
		if( lgSaveRomasUsed )
			fclose( ioROMASused );

		if( trace.lgTrace )
			fprintf( ioQQQ, " reading UTA_Kisielius.dat OK,used %li lines from a total of %li\n" , nRomasUsed , nRomasTotal );
	}

	/* option to dump UTA lines, either save file, or in main output */
	if( trace.lgTrace || save.lgSDSOn || atmdat.lgUTAprint )
	{
		FILE *ioUTA = ioQQQ;
		if( save.lgSDSOn )
			ioUTA = save.ipSDSFile;

		fprintf(ioUTA, "##################################################\n");
		fprintf(ioUTA,"UTA data sources; B=Badnell 05; G==Gu 06, K=Kisielius 03, 13\n");
		fprintf(ioUTA," ion ");
		for( long ion=0; ion<=LIMELM; ++ion )
			fprintf(ioUTA,"%4li",ion);
		fprintf(ioUTA,"\n");
		for( long nelem=0; nelem<LIMELM; ++nelem )
		{
			fprintf(ioUTA,"%4s ", elementnames.chElementSym[nelem] );
			for( long ion=0; ion<=nelem; ++ion )
			{
				fprintf(ioUTA,"%4s",chUTA_ref[nelem][ion] );
			}
			fprintf(ioUTA,"\n");
		}
		fprintf(ioUTA," ion ");
		for( long ion=0; ion<=LIMELM; ++ion )
			fprintf(ioUTA,"%4li",ion);
		fprintf(ioUTA,"\n");
		fprintf(ioUTA,"Badnell 05=2005MNRAS.360..458B; Gu 06=2006ApJ...641.1227G; Kisielius 03, 13= 2003MNRAS.344..696K, 2013ApJ...767..123F\n");
		fprintf(ioUTA, "##################################################\n\n");
	}

	if( false )
	{
		for( size_t iu=0; iu < UTALines.size(); ++iu )
		{
			string chLab = chIonLbl( UTALines[iu] );
			dprintf( ioQQQ, "%5ld %s wavl %7.3f glo %2g gup %2g Aul %.2e gf %.2e ai branch %.3f\n",
				 (long)iu,
				 chLab.c_str(),
				 UTALines[iu].WLAng(),
				 (*UTALines[iu].Lo()).g(),
				 (*UTALines[iu].Hi()).g(),
				 UTALines[iu].Emis().Aul(),
				 UTALines[iu].Emis().gf(),
				 UTALines[iu].Emis().AutoIonizFrac() );
		}
		cdEXIT(EXIT_SUCCESS);
	}
}

static const long IP_MAGIC = 20221115;

STATIC void read_ionization_potentials()
{
	DEBUG_ENTRY( "read_ionization_potentials()" );

	DataParser d( "ionization_potentials.dat", ES_STARS_ONLY );
	d.getline();
	d.checkMagic(IP_MAGIC);
	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		for( long ion=0; ion <= nelem; ++ion )
		{
			d.getline();
			long n1, i0;
			double e, de;
			d.getToken(n1);
			d.getToken(i0);
			d.getToken(e);
			d.getToken(de);
			d.checkEOL();
			ASSERT( n1 == nelem+1 );
			ASSERT( i0 == ion );
			atmdat.EIonPot[nelem][ion] = Energy(e, "cm^-1").Ryd();
		}
		for( long ion=nelem+1; ion < LIMELM; ++ion )
		{
			atmdat.EIonPot[nelem][ion] = -1.;
		}
	}
}
