/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program to call cloudy to generate a large loc data cube */
#include "cddefines.h"
#include "cddrive.h"
#include "lines.h"
#include "prt.h"
#include "version.h"

/* This flag indicates where we are multi-processing.\
 * set with -DMPI on compiler command line
 * if this is true then all output happens at the end,
 * if false then each set of results is printed and flushed when it happens.
 * must be set true on parallel machines, false on pc */
#ifdef MPI 
#	include <mpi.h>
#endif
/*lint -e506  constant value Boolean */


#ifdef _MSC_VER
	/* disable conditional expression is constant */
#	pragma warning( disable : 4127 )
#endif

/* ========================================================================= */
/*
Main output, file ending in .lin 
format, tab separated fields
Cloudy exit condition for that model
number of warnings,
elapsed time, 
hden
flux
par1
par2
*/

#define TITLE "std con Fe2 sol, col 24 " 

/* this is name of resulting output files, FILENAME.lin will be main results */
#define FILENAME "STD_con_sol" 

/* path to directory where we will write, must end in / or \ for current system
 * final location will be PATHNAME + FILENAME */
/*#define PATHNAME "c:\\projects\\cloudy\\tests\\"*/
#define PATHNAME ""

/* this is an additive offset for hden and flux, usually 0 */
#define OFFSET_HDEN (0.)
#define OFFSET_FLUX (0.)

/* these following pairs make 8x8 square grid */
/* lower and upper bounds on hydrogen density - log cm^-3
 * for BLR work usually 7 and 14 */
/*#define HDENINIT (12.5+OFFSET_HDEN)*/
/*#define HDENLIMIT (12.5+OFFSET_HDEN)*/
#define HDENINIT (7.+OFFSET_HDEN)
#define HDENLIMIT (14.+OFFSET_HDEN)

/* low and up bounds on flux of hydrogen ionizing photons - log cm^-2 s-1 
 * for BLR work usually 17 and 24 */
#define FLUXINIT (17.+OFFSET_FLUX) 
#define FLUXLIMIT (24.+OFFSET_FLUX)
/*#define FLUXINIT (20.5+OFFSET_FLUX)
#define FLUXLIMIT (20.5+OFFSET_FLUX)*/

/* the increment on both axis */
/*#define INCREMENT 0.5*/
#define INCREMENT 1./**/
/*#define INCREMENT 2.*/
/*#define INCREMENT 0.05*/

/* this is the number of iterations to do - <0 means to convergence, 1 will do 2 iterations,
 * 0 will do 1  */
#define ITERATIONS  2  
/*#define ITERATIONS  -1  */

/* the column density, set to false if not used */
#define COLDEN 24

/* this says whether or not to include the FeII atom */
#define FEII true

/* quick test - 3 zones, constant temperature, 2 iterations */
#define QUICK_MODEL false

/* if true the only call cdNoExec */
#define VERY_QUICK_MODEL false

/* says whether to make the large cloudy output, if true then output created,
 * if false then not output, only the emission line intensities. */
#define DO_PRINT false

/* print last iteration */
#define PRINT_LAST false

/* set no buffering on output - if true then print as it happens*/
#define NO_BUFFERING false

/* ============================================================================ */

/* this structure will hold results for final printout */
#define NLINTOT 1500 /* the maximum number of lines */

typedef struct
{
	/* save these parameters for each grid point,
	 * density, flux, and full spectrum */
	double hden , flux , pred[NLINTOT] ;
	/* number of cautions and warnings for this calc, all should be zero */
	int nWarnings , nCautions , nTFail;
	/* flag returned by Cloudy, 0 == OK, 1 for failure */
	int exit_status;
	/* the execution time for this model */
	double etime;
	/* two extra parameters */
	double par1 , par2;
} GRIDTAG;

/* this is where we will save results */
#ifdef MPI
void Build_derived_type(GRIDTAG gridtag, MPI_Datatype* message_type_ptr)
{
  int block_lengths[10];
  MPI_Aint displacements[10];
  MPI_Aint addresses[11];
  MPI_Datatype typelist[10];

  /*First specify the types */
  typelist[0] = MPI_DOUBLE;    /*hden*/
  typelist[1] = MPI_DOUBLE;    /*flux*/
  typelist[2] = MPI_DOUBLE;    /*pred*/
  typelist[3] = MPI_INT;       /*nWarnings*/
  typelist[4] = MPI_INT;       /*nCautions*/
  typelist[5] = MPI_INT;       /*nTFail*/
  typelist[6] = MPI_INT;       /*exit_status*/
  typelist[7] = MPI_DOUBLE;    /*etime*/
  typelist[8] = MPI_DOUBLE;    /*par1*/
  typelist[9] = MPI_DOUBLE;    /*par2*/

  /* Specify the number of elements of each type */
  block_lengths[0] = block_lengths[1] = block_lengths[3] = 1;
  block_lengths[4] = block_lengths[5] = block_lengths[6] = 1;
  block_lengths[7] = 1;
  /* this one is special since it is NLINTOT long */
  block_lengths[2] = NLINTOT;
  /* two new parameters are each one element long */
  block_lengths[8] = block_lengths[9] = 1;

  /* Calculate the displacement of the members relative to indata */
  MPI_Address( &gridtag,             &addresses[0]);
  MPI_Address(&(gridtag.hden),       &addresses[1]);
  MPI_Address(&(gridtag.flux),       &addresses[2]);
  MPI_Address(&(gridtag.pred),       &addresses[3]);
  MPI_Address(&(gridtag.nWarnings),  &addresses[4]);
  MPI_Address(&(gridtag.nCautions),  &addresses[5]);
  MPI_Address(&(gridtag.nTFail),     &addresses[6]);
  MPI_Address(&(gridtag.exit_status),&addresses[7]);
  MPI_Address(&(gridtag.etime),      &addresses[8]);
  MPI_Address(&(gridtag.par1),       &addresses[9]);
  MPI_Address(&(gridtag.par2),       &addresses[10]);
 
  /* now far each is from the beginning of the structure */
  displacements[0] = addresses[1] - addresses[0];
  displacements[1] = addresses[2] - addresses[0];
  displacements[2] = addresses[3] - addresses[0];
  displacements[3] = addresses[4] - addresses[0];
  displacements[4] = addresses[5] - addresses[0];
  displacements[5] = addresses[6] - addresses[0];
  displacements[6] = addresses[7] - addresses[0];
  displacements[7] = addresses[8] - addresses[0];
  displacements[8] = addresses[9] - addresses[0];
  displacements[9] = addresses[10]- addresses[0];

  /*Create the derived type */
  MPI_Type_struct(10, block_lengths, displacements, typelist,message_type_ptr);
  
  /*Commit it so that it can be used */
  MPI_Type_commit(message_type_ptr);

} /*Build_derived_type */

# endif

int main( int argc, char *argv[] )
{
	DEBUG_ENTRY( "main()" );

	exit_type exit_status = ES_SUCCESS;

	try {
		int IntegerMod;
		int  myrank ;
		double hden , absolute , relative , flux;
		char chVer[10] , chPath[1000] , chFilename[1000] , chFile[2000] ;

		double *xpar , *ypar, *par1, *par2;

		/* following will become array of line wavelengths and the number of lines 
		 * in this array */
		long int mod, nLines, LimModels, nModels;

		/* these will be passed to cdGetLineList and will become arrays of
		 * labels and wavelengths */
		vector<LineID> lineids;

		long int 
			NumberWarnings1, NumberCautions1, NumberNotes1, NumberSurprises1, NumberTempFailures1,
			NumberPresFailures1, NumberIonFailures1, NumberNeFailures1 ;
		/* number of time Cloudy returned with error condition */
		int nErrorExits;

		FILE *ioDATA ;
		char chLine[100];

		long int n;
		/* number of processors, =1 for single non-MPI */
		int numprocs = 1;

		/* this will hold a single grid point */
		GRIDTAG grid;
		grid.exit_status = ES_SUCCESS;

		/* start MPI if -DMPI on command line */
#		ifdef MPI
		GRIDTAG gridtag;
		GRIDTAG *grids;
		
		MPI_Datatype message_type; /* Arguments to */

		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

		Build_derived_type(gridtag, &message_type);
#		else
		myrank = 0;
#		endif

		/* the following is the path for the directory where the output should appear,
		 * it needs to end with the OS standard directory delimiter, a "/" on unix */
		strcpy( chPath , PATHNAME );

		/* this is the first part of the name of the resulting data file.  The final
		 * file will have this name with a ".lin" after it. */
		strcpy( chFilename , FILENAME );

		/* ====================================================================== */
		/*
		 * make arrays of parameters for this grid */

		/* total number of models */
		LimModels = (long)((FLUXLIMIT - FLUXINIT)/INCREMENT + 1.1);
		LimModels *= (long)((HDENLIMIT - HDENINIT)/INCREMENT+1.1);
		/* add number of processors since this is most we would possibly have to 
		 * add on to make integer number of mpi calls */
		LimModels += numprocs;

		xpar = (double*)malloc(LimModels*sizeof(double) );
		ypar = (double*)malloc(LimModels*sizeof(double) );
		par1 = (double*)malloc(LimModels*sizeof(double) );
		par2 = (double*)malloc(LimModels*sizeof(double) );

		if( xpar==NULL || ypar==NULL || par1==NULL || par2==NULL )
		{
			fprintf(stderr,"malloc failed\n");
#			ifdef MPI 
			MPI_Finalize(); 
#			endif
			cdEXIT(EXIT_FAILURE);
		}
		/* now define set of parameters */
		hden = HDENINIT;
		flux = FLUXINIT;
		nModels = 0;
		nErrorExits = 0;
#		ifdef MPI 
		/* loop along constant U diagonals, starting from
		 * bottom right of loc plane, increasing density, 
		 * and ending on the major diagonal.
		 * this is for load leveling on parallel machines */
		for( hden = HDENLIMIT; hden>= HDENINIT; hden-=INCREMENT )
		{
			for( flux=FLUXINIT; flux<=FLUXLIMIT;  flux+= INCREMENT )
			{
				if( hden + (flux-FLUXINIT) > HDENLIMIT ) break;
				assert( nModels < LimModels );
				xpar[nModels] = hden + (flux-FLUXINIT);
				ypar[nModels] = flux;
				par1[nModels] = COLDEN;
				par2[nModels] = 0.;
				++nModels;
			}
		}
		/* loop increasing flux, starting at lower left plus increment */
		for( flux=FLUXINIT+INCREMENT; flux<=FLUXLIMIT;  flux+= INCREMENT )
		{
			for( hden = HDENINIT; hden<= HDENLIMIT;  hden+=INCREMENT )
			{
				if( flux + (hden-HDENINIT) > FLUXLIMIT ) break;
				assert( nModels < LimModels );
				xpar[nModels] = hden ;
				ypar[nModels] = flux + (hden-HDENINIT);
				par1[nModels] = COLDEN;
				par2[nModels] = 0.;
				++nModels;
			}
		}
#		else
		/* on scalar machines do models in xy order */
		while( hden < 1.00001 * HDENLIMIT )
		{
			while( flux < 1.00001*FLUXLIMIT )
			{
				xpar[nModels] = hden;
				ypar[nModels] = flux;
				par1[nModels] = COLDEN;
				par2[nModels] = 0.;
				++nModels;
				flux += INCREMENT;
			}
			flux = FLUXINIT;
			hden += INCREMENT;
		}
#		endif
		/* increase total number to an integer multiple of number of processors */
		if( nModels%numprocs)
		{
			/* add extra bit */
			IntegerMod = nModels + (numprocs - nModels%numprocs);
		}
		else
		{
			IntegerMod = nModels;
		}

		/* make up data for the "extra" models */
		for( mod=nModels; mod<IntegerMod; ++mod)
		{
			xpar[mod] = xpar[0];
			ypar[mod] = ypar[0];
			par1[mod] = par1[0];
			par2[mod] = par2[0];
		}

		/* ====================================================================== */

#		ifdef MPI
		/* allocate memory for grids */	
		if( (grids = ((GRIDTAG *)malloc( (size_t)(numprocs)*sizeof(GRIDTAG )))) == NULL )
		{
			fprintf(stderr," main could not malloc grid\n");
#			ifdef MPI 
			MPI_Finalize(); 
#			endif
			cdEXIT(EXIT_FAILURE);
		}
#		endif

		/* initialize the code so that cdGetLineList can be called */
		cdInit();
		/* get list of lines from standard data file.  With no arguments this reads
		 * the file BLRLineList, part of the standard data distribution.  There could
		 * have been another file name within the quotes - that file would have been
		 * used instead */
		nLines = cdGetLineList("", lineids);

		/* now copy version into string to print */
		cdVersion(chVer);

		/* open file for large output if this is desired */
		if( DO_PRINT )
		{
			/* now create final filename */
			strcpy( chFile , chPath );
			strcat( chFile , chFilename );
			strcat( chFile , ".out" );
			cdOutput( chFile );
		}
		else
		{
			/* there should be nothing coming here, but something might */
			ioQQQ = stderr;
		}

		/* calculation's results, all the lines in file ending .lin */
		strcpy( chFile , chPath );
		strcat( chFile , chFilename );
		strcat( chFile , ".lin" );
		ioDATA = open_data(chFile,"w");

		/* print the header information before we begin the calculation */
		/* print version number on both headers */
		fprintf(ioDATA,"Cloudy version %s, %s\n",chVer,TITLE);

		for( mod=myrank; mod<IntegerMod; mod+=numprocs )
		{
			hden = xpar[mod];
			flux = ypar[mod];
			/* zero out grid, which will hold all info for a single grid point */
			memset(&grid, 0, sizeof(GRIDTAG));

			/* initialize the code for this run */
			cdInit();

			/* send flag saying whether to create large output */
			cdTalk( DO_PRINT );

			/* title for grid */
			sprintf(chLine,"title %s", TITLE);
			n = cdRead( chLine  );

			/* print only last iteration? */
#			if PRINT_LAST
			n = cdRead("print last iteration ");
#			endif

			/* option for fast model, or only header */
#			if QUICK_MODEL
			n = cdRead( "stop zone 3 "  );/**/
			n = cdRead( "constant temper 4 "  );/**/
#			endif
#			if VERY_QUICK_MODEL
			cdNoExec();
#			endif

			/* iterate if not fast model */
#			if ( (!QUICK_MODEL) && (!VERY_QUICK_MODEL) &&(ITERATIONS!=0) )
#      			if ITERATIONS < 0
			n = cdRead( "iterate convergence "  );/**/
#			else 
			sprintf(chLine,"iterate %i", ITERATIONS);
			n = cdRead( chLine  );
#			endif
#			endif

#			ifdef MPI 
			/* set flag so that exit handler will call MPI_Finalize */
			cpu.set_MPI();
#			endif

			/* option to turn off file buffering if code might crash */
#			ifdef NO_BUFFERING
			n = cdRead("no buffering ");
#			endif

			/* if we have run off end of array of models say not to execute,
			 * still must do something on this processor so that the gather can occur
			 * at the bottom of the loop */
			if( mod >= nModels )
				cdNoExec();

			/* flux and density */
			sprintf(chLine,"phi(h) %f", flux);
			n = cdRead( chLine  );

			sprintf(chLine,"hden %f", hden);
			n = cdRead( chLine  );

			/* turn on large FeII atom, VERY slow */
#			if FEII
			n = cdRead( "atom feii  "  );
#			endif

			/* this is the continuum we used in the large apjs paper */
			n = cdRead( "agn 6.00 -1.40 -0.50 -1.0  "  );/**/

			/* a bare powerlaw with a break at 1 microns */
			/*n = cdRead( "table powerlaw -1.0 1 micron break "  );*/

			/* the much maligned Mathews&Ferland continuum */
			/*n = cdRead( "table agn  "  );*/

			/*n = cdRead( "blackbody 142,000  "  );*/

			/* broken power law used in Fred's paper */
			/*n = cdRead( "interpolate (0.00000001, -14.899), (0.0091126732, 0.000)  "  );
			  n = cdRead( "continue (1.00, -1.2242), (73.54, -4.2106)  "  );
			  n = cdRead( "continue (3675, -5.7395), (7354000, -11.2526)  "  );*/

			/* broken power law used in Mark's paper */
			/*n = cdRead( "interpolate (10.517, -30.850) (13.383, -23.684)   "  );
			  n = cdRead( "continue (16.684, -26.986) (17.383, -27.9996)    "  );
			  n = cdRead( "continue (18.383, -28.8899) (19.124, -29.268)    "  );
			  n = cdRead( "continue (22.383, -35.788)  "  );*/

			/* add cosmic IR background as well */
			n = cdRead( "background z=0  "  );

			sprintf(chLine,"stop column density %f ", par1[mod]);
			n = cdRead( chLine  );/**/

			n = cdRead( "failures 3 "  );

			/* options to change abundances */
#			if 0
			n = cdRead( "metals 10 "  );/**/
			n = cdRead( "element nitrogen scale 10 "  );/**/
			n = cdRead( "element helium scale 1.66 "  );/**/
			n = cdRead( "metals 20 "  );/**/
			n = cdRead( "element nitrogen scale 20 "  );/**/
			n = cdRead( "element helium scale 2.41 "  );/**/
#			endif

			/*sprintf(chLine,"turbulence %f km/sec", vturb);
			  n = cdRead( chLine  );*/

			/* all line intensities will be relative to incident continuum near Lya */
			n = cdRead( "normalize to label \"inci\" 1215 "  );

			/* actually call the code, true indicates error exit 
			 * this is always an abort */
			if( cdDrive() )
			{
				grid.exit_status = ES_FAILURE;
				nErrorExits++;
			}

			/* print header information for very first model only */
			if( mod==0 )
			{
				if( (nLines+nFeIIBands) > NLINTOT )
				{
					fprintf( stderr , "there are too many lines - increase NLINTOT ");
					printf(  "there are too many lines - increase NLINTOT ");
#					ifdef MPI 
					MPI_Finalize(); 
#					endif
					cdEXIT(EXIT_FAILURE);
				}

				/* for very first model (mod==0) send input stream to data file, followed
				 * by keyword showing end and finally set of tables */
				cdPrintCommands( ioDATA );

				/* all the printout happens here */
				/* print header for the data file */
				fprintf(ioDATA,"abort\twarn\tExecT\tdensity\tflux\tcolden\tpar2\t");
				for( n=0; n<nLines; ++n )
				{
					fprintf(ioDATA,"%4s ",lineids[n].chLabel.c_str());
					prt_wl(ioDATA,lineids[n].wave);
					fprintf(ioDATA,"\t");
				}
#				if FEII
				/* the header information for the FeII bands */
				for( n=0; n<nFeIIBands; ++n )
				{
					fprintf(ioDATA,"fe2b %li\t", (long)FeII_Bands[n][0] );
				}
				/* the header information for FeII continuum bins
				 * 1 and 2 are the lower and upper bounds, 0 the center */
				for( n=0; n<nFeIIConBins; ++n )
				{
					fprintf(ioDATA,"fe2c %.1f\t", ( FeII_Cont[n][1] + FeII_Cont[n][2] ) /2. );
				}
#				endif
				fprintf(ioDATA,"\n");
			}

			/* flush this output */
			if( DO_PRINT )
				fflush(ioQQQ);

			/* keep track of all comments about the calculation */
			cdNwcns( 
				&NumberWarnings1 ,   &NumberCautions1,     &NumberNotes1, 
				&NumberSurprises1,   &NumberTempFailures1, &NumberPresFailures1,
				&NumberIonFailures1, &NumberNeFailures1 );
			grid.nWarnings = NumberWarnings1;
			grid.nCautions = NumberCautions1;
			grid.nTFail = NumberTempFailures1;

			/* print the exec time */
			grid.etime = cdExecTime();

			/* for very quick MPI run, remember processor number */
#			if VERY_QUICK_MODEL && MPI
			grid.nWarnings = myrank;
			grid.nCautions = mod;
#			endif

			/* save grid parameters for this model */
			grid.hden = hden;
			grid.flux = flux;
			grid.par1 = par1[mod];
			grid.par2 = par2[mod];

			/* if we have run off end of array of models do not try to pull out
			 * results since they do not exist*/
			if( mod < nModels )
			{
				/* print line intensiies */
#				if !VERY_QUICK_MODEL
				for( n=0; n<nLines; ++n )
				{
					if( (cdLine( lineids[n] , &relative , &absolute ))<=0 )
					{
						fprintf(stderr,"did not find %4s",lineids[n].chLabel.c_str());
						prt_wl(stderr,lineids[n].wave);
						fprintf(stderr,"\n");
						fprintf(ioDATA,"\ndid not find %4s",lineids[n].chLabel.c_str());
						prt_wl(ioDATA,lineids[n].wave);
						fprintf(ioDATA,"\n");
					}
					grid.pred[n] = log10(MAX2(1e-30,relative) );
				}
#				if FEII
				/* first the FeII bands */
				for( n=0; n<nFeIIBands; ++n )
				{
					int lgOK;
					if( (lgOK=cdLine( "Fe2b", FeII_Bands[n][0] , &relative , &absolute ))<=0 )
					{
						fprintf(stderr,"did not find %s %.1f %i\n","Fe2b",FeII_Bands[n][0],lgOK);
						fprintf(ioDATA,"did not find %s %.1f %i\n","Fe2b",FeII_Bands[n][0],lgOK);
					}
					grid.pred[n+nLines] = log10(MAX2(1e-30,relative) );
				}
				/* next the FeII continuum */
				if( cdLine( "inci", 1215 , &relative , &absolute )<=0 )
				{
					fprintf(stderr,"did not find incident continuum\n");
					fprintf(ioDATA,"did not find incident continuum\n");
				}
				for( n=0; n<nFeIIConBins; ++n )
				{
					assert( n+nLines+nFeIIBands < NLINTOT);
					grid.pred[n+nLines+nFeIIBands] = -30.;
					if( absolute > 0. )
						grid.pred[n+nLines+nFeIIBands] = log10(MAX2(1e-30,FeII_Cont[n][0]/ absolute) );
				}
#				endif
#				endif
			}

			/* print the results here if not MPI */
#			ifndef MPI
			/* print exec time and numbers of problems */
			fprintf(ioDATA,"%i\t%i\t%g\t" , grid.exit_status,	grid.nWarnings , grid.etime);

			/* print grid parameters for this model */
			fprintf(ioDATA,"%.3f\t%.3f\t%.3f\t%.3f\t",grid.hden, grid.flux,
				grid.par1 , grid.par2 );

#			if !VERY_QUICK_MODEL
			/* print main set of line intensiies */
			for( n=0; n<nLines; ++n )
			{
				fprintf(ioDATA,"%.3f\t",grid.pred[n] );
			}
#			if FEII
			/* print out the bands and continuum bins */
			for( n=0; n<nFeIIBands+nFeIIConBins; ++n )
			{
				fprintf(ioDATA,"%.3f\t",grid.pred[n+nLines]  );
			}
#			endif
#			endif
			fprintf(ioDATA,"\n");
			/* flush this output */
			fflush(ioDATA );
#			endif

#			ifdef MPI

			MPI_Gather(&grid,1,message_type,grids,1,message_type,0,MPI_COMM_WORLD);

			if(myrank == 0)
			{
				{
					long int i , nmax ;
					/* on last loop can't expect to get numprocs model results */
					nmax = MIN2( numprocs, nModels-mod);
					for( i=0; i<nmax; ++i )
					{
						/* print exec time and numbers of problems */
						fprintf(ioDATA,"%i\t%i\t%g\t",grids[i].exit_status, 
							grids[i].nWarnings,grids[i].etime );
		
						/* print grid parameters for this model */
						fprintf(ioDATA,
							"%7.3f\t%7.3f\t%7.3f\t%7.3f\t",
							grids[i].hden, grids[i].flux ,
							grids[i].par1 , grids[i].par2  );
#						if !VERY_QUICK_MODEL
						/* print line intensiies */
						for( n=0; n<nLines; ++n )
						{
							fprintf(ioDATA,"%.3f\t",grids[i].pred[n] );
						}
#						if FEII
						/* print out the bands and continuum bins */
						for( n=0; n<nFeIIBands+nFeIIConBins; ++n )
						{
							fprintf(ioDATA,"%.3f\t",grids[i].pred[n+nLines]  );
						}
#						endif
#						endif
						fprintf(ioDATA,"\n");
						/* flush this output */
						fflush(ioDATA );
					}
				}
			}
#			endif
		}
	
		free(xpar );
		free(ypar );
		free(par1 );
		free(par2 );

		/* call MPI_Finalize if MPI is set */
#		ifdef MPI 
		free(grids);
		/* only print message if errors occurred */
		if( nErrorExits )
		{
			fprintf( stderr, "processor %i main exit %i aborts\n" ,myrank, nErrorExits );
			if( DO_PRINT )
			{
				fprintf( ioQQQ, "processor %i main exit %i aborts\n" ,myrank, nErrorExits );
			}
		}
		MPI_Finalize(); 
#		else
		fprintf( stderr, "exit in main with %i aborted models\n" , nErrorExits );
		if( DO_PRINT )
		{
			fprintf( ioQQQ, "exit in main with %i aborted models\n" , nErrorExits );
		}
#		endif

		cdEXIT(exit_type(grid.exit_status));
	}
	catch( bad_alloc& )
	{
		fprintf( ioQQQ, " DISASTER - A memory allocation has failed. Most likely your computer "
			 "ran out of memory.\n Try monitoring the memory use of your run. Bailing out...\n" );
		exit_status = ES_BAD_ALLOC;
	}
	catch( out_of_range& e )
	{
		fprintf( ioQQQ, " DISASTER - An out_of_range exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_OUT_OF_RANGE;
	}
	catch( domain_error& e )
	{
		fprintf( ioQQQ, " DISASTER - A vectorized math routine threw a domain_error. Bailing out...\n" );
		fprintf( ioQQQ, " What() = %s", e.what() );
		exit_status = ES_DOMAIN_ERROR;
	}
	catch( bad_assert& e )
	{
		MyAssert( e.file(), e.line() , e.comment() );
		exit_status = ES_BAD_ASSERT;
	}
	catch( bad_signal& e )
	{
		if( e.sig() == SIGILL )
		{
			if( ioQQQ != NULL )
				fprintf( ioQQQ, " DISASTER - An illegal instruction was found. Bailing out...\n" );
			exit_status = ES_ILLEGAL_INSTRUCTION;
		}
		else if( e.sig() == SIGFPE )
		{
			if( ioQQQ != NULL )
				fprintf( ioQQQ, " DISASTER - A floating point exception occurred. Bailing out...\n" );
			exit_status = ES_FP_EXCEPTION;
		}
		else if( e.sig() == SIGSEGV )
		{
			if( ioQQQ != NULL )
				fprintf( ioQQQ, " DISASTER - A segmentation violation occurred. Bailing out...\n" );
			exit_status = ES_SEGFAULT;
		}
#		ifdef SIGBUS
		else if( e.sig() == SIGBUS )
		{
			if( ioQQQ != NULL )
				fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" );
			exit_status = ES_BUS_ERROR;
		}
#		endif
		else
		{
			if( ioQQQ != NULL )
				fprintf( ioQQQ, " DISASTER - A signal %d was caught. Bailing out...\n", e.sig() );
			exit_status = ES_UNKNOWN_SIGNAL;
		}
	}
	catch( cloudy_abort& e )
	{
		fprintf( ioQQQ, " ABORT DISASTER PROBLEM - Cloudy aborted, reason: %s\n", e.comment() );
		exit_status = ES_CLOUDY_ABORT;
	}
	catch( cloudy_exit& e )
	{
		if( ioQQQ != NULL )
		{
			ostringstream oss;
			oss << " [Stop in " << e.routine();
			oss << " at " << e.file() << ":" << e.line();
			if( e.exit_status() == 0 )
				oss << ", Cloudy exited OK]";
			else
				oss << ", something went wrong]";
			fprintf( ioQQQ, "%s\n", oss.str().c_str() );
		}
		exit_status = e.exit_status();
	}
	catch( std::exception& e )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}
	// generic catch-all in case we forget any specific exception above... so this MUST be the last one.
	catch( ... )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught. Bailing out...\n" );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}

	cdPrepareExit(exit_status);

	return exit_status;
}
