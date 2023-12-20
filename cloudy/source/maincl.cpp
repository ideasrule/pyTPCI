/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program that reads input and calls cloudy to compute a single model, or
 * try to optimize an observed model.  Routine returns 0 if model is ok, 
 * and 1 if problems occurred. */
#include "cddefines.h"
#include "cddrive.h"
#include "input.h"
#include "prt.h"
#include "save.h"
#include "called.h"
#include "monitor_results.h"
#include "grid.h"
#include "service.h"
#include "ran.h"

static const uint64 default_seed = 0xc7f8f57fe95956a8ULL;

exit_type cdMain( int argc, const char* argv[] );

inline void print_delimiter(long nOptimiz)
{
	fprintf( ioQQQ, " ************************************************** GRID_DELIMIT" );
	if( nOptimiz >= 0 )
		fprintf( ioQQQ, " -- grid%9.9ld", nOptimiz );
	fprintf( ioQQQ, "\n" );
}

/** main: this is a wrapper around cdMain. It takes care of the MPI stuff
 * for non-MPI runs, this should do nothing more than call cdMain and exit. */
int main( int argc, char *argv[] )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" );

	bool lgMPI = cpu.i().lgMPI();

	cpu.i().set_used_nCPU( 1 );
	if( lgMPI )
	{
		MPI_Init( &argc, &argv );

		int nCPU = 1, nRANK = 0;
		MPI_Comm_size( MPI_COMM_WORLD, &nCPU );
		cpu.i().set_nCPU( nCPU );
		cpu.i().set_used_nCPU( nCPU );
		MPI_Comm_rank( MPI_COMM_WORLD, &nRANK );
		cpu.i().set_nRANK( nRANK );

		// MPI_Init() will have overwritten our signal handlers
		// so we need to set them again....
		cpu.i().set_signal_handlers();
	}

	// save headers are needed for sequential runs
	// saving headers during the initial parsing of grid runs is also OK
	save.SetSaveHeaderNeeded( true );

	cpu.i().set_ExecName( argv[0] );

	uint64 seed = 0;
	bool lgSeedSet = false;
	/* Handle -s argument for setting fixed seed */
	for( int i=1; i < argc; i++ ) 
	{
		string s( argv[i] );
		if( s.substr(0,2) == "-e" )
		{
			// no other flags are allowed after "-e"...
			break;
		}
		else if( s.substr(0,2) == "-s" )
		{
			string t;
			if( s.length() > 2 )
				t = s.substr(2);
			else if( i < argc-1 && argv[i+1][0] != '-' )
				t = argv[++i];
			if( t.empty() )
				seed = default_seed;
			else
			{
				istringstream iss( t );
				iss >> hex >> seed;
				if( iss.fail() )
				{
					fprintf( ioQQQ, "Reading seed failed in \"%s\".\n", t.c_str() );
					exit(EXIT_FAILURE); // this _should_ be exit()...
				}
			}
			lgSeedSet = true;
			break;
		}
	}

	// now intialize the random number generator
	if( cpu.i().lgMaster() )
	{
		if( !lgSeedSet )
			ran.init();
		else
			ran.init(seed, 0);
		seed = ran.get_seed();
	}
	// and copy the seed to other ranks so they can seed as well
	if( lgMPI )
	{
		MPI_Bcast( &seed, 1, MPI_type(seed), 0, MPI_COMM_WORLD );
		// using nRANK guarantees that each rank has a different random number sequence
		ran.init(seed, cpu.i().nRANK());
	}
	
	// this will generate input files for each grid point,
	// or execute the Phymir run, whichever is appropriate
	exit_status = cdMain( argc, (const char**)argv );

	// wait for writing of input files to finish
	if( lgMPI )
		MPI_Barrier( MPI_COMM_WORLD );

	// process the individual grid points
	if( grid.lgGrid && exit_status == ES_SUCCESS )
	{
		// this was set to true after we wrote the last input script
		grid.lgGridDone = false;
		// signal that we are running individual grid models now...
		grid.lgInsideGrid = true;

		// from now on each rank will run its own model
		cpu.i().set_MPIMode( MS_GRID );

		unsigned int nCPU;
		if( lgMPI )
			nCPU = cpu.i().nCPU();
		else
			nCPU = grid.lgParallel ? grid.useCPU : 1;
		cpu.i().set_used_nCPU( nCPU );

		load_balance lb( grid.totNumModels, nCPU );

		// Each MPI rank will get jobs assigned by lb and execute them.
		// If there are no jobs left, lb.next_job() will return -1.
		exit_type retval = ES_SUCCESS;
		while( ( optimize.nOptimiz = lb.next_job() ) >= 0 )
		{
			const char** mpi_argv = new const char*[argc+2];

			string jobName = GridPointPrefix( optimize.nOptimiz );
			for( int i=0; i < argc; ++i )
				mpi_argv[i] = argv[i];
			mpi_argv[argc] = "-g";
			mpi_argv[argc+1] = jobName.c_str();

			// only save header if we are calculating the first grid model
			save.SetSaveHeaderNeeded( optimize.nOptimiz == 0 );

			retval = cdMain( argc+2, mpi_argv );

			exit_status = max( retval, exit_status );
			delete[] mpi_argv;

			++grid.seqNum;
		}

		lb.finalize( exit_status );

		grid.lgGridDone = true;

		if( retval != ES_SUCCESS )
		{
			// Now parse the main input file once again to make sure that
			// process_ouput() has the correct information available. This
			// may be needed if the last grid model of this rank failed
			// during parsing, which could lead to incorrect settings...
			// Parsing is guaranteed to succeed since this was already
			// done once before, and therefore we ignore the return value.
			cpu.i().set_MPIMode( MS_POST_GRID );
			(void)cdMain( argc, (const char**)argv );
		}

		// concatenate the output
		process_output();
	}

	// remove empty output files from slave ranks
	if( lgMPI && cpu.i().lgMaster() )
	{
		for( long n=1; n < cpu.i().nCPU(); ++n )
		{
			ostringstream oss;
			oss << ".err" << setfill('0') << setw(2) << n;
			string slave_output = save.chRedirectPrefix + oss.str();
			FILE *io = open_data( slave_output, "a" );
			bool lgEmpty = ( ftell(io) == 0 );
			fclose( io );
			if( lgEmpty )
				remove( slave_output.c_str() );
		}
	}

	if( lgMPI )
		MPI_Finalize();

	return exit_status;
}

/** cdMain: this is the main entry point for Cloudy */
exit_type cdMain( int argc, const char* argv[] )
{
	/* these will be used to count number of various problems */
	long int NumberWarnings, 
	  NumberCautions, 
	  NumberNotes, 
	  NumberSurprises, 
	  NumberTempFailures, 
	  NumberPresFailures,
	  NumberIonFailures, 
	  NumberNeFailures;

	bool lgEarly_exit=true,
	  lgCommandLineScript=false,
	  lgFileIO;

	int i;
	const char *s, 
	  *prefix = "",
	  *gprefix = "", // grid prefix
	  *pprefix = "", // save prefix
	  *rprefix = ""; // redirect prefix
	string infile("");
	char *outfile = NULL;

	/* indicates that a command line flag to redirect I/O has been used */
	lgFileIO = false;

	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "cdMain()" );

	try {
		/* initialize the code for this run */
		cdInit();

		/* Handle argument input */
		for( i=1; i < argc; i++ ) 
		{
			s = argv[i];
			if( *s != '-' || s[1] == '\0' ) 
			{
				if( infile != "" )
				{
					fprintf( ioQQQ, "%s: only one input file argument allowed\n", argv[0] );
					cdEXIT(ES_FAILURE);
				}
				infile = s;
				if( infile.find( cpu.i().chDirSeparator() ) != string::npos )
				{
					fprintf( ioQQQ, "%s %s: read/write from subdirectories is not supported\n",
								argv[0], infile.c_str() );
					cdEXIT(ES_FAILURE);
				}
				if( infile != "-" )
				{
					FILE *fp = open_data( infile, "r", AS_LOCAL_ONLY_TRY );
					if( fp == NULL )
					{
						fprintf( ioQQQ, " input file \"%s\" not found\n", infile.c_str() );
						cdEXIT(ES_FAILURE);
					}
					fclose(fp);
					size_t suffindex = infile.find_last_of(".");
					if (suffindex != string::npos)
						infile = infile.substr(0,suffindex);
					pprefix = rprefix = infile.c_str();
					lgFileIO = true;
				}
			}
			else
			{
				while( s != NULL && *(++s) )
				{
					exit_type exit = ES_SUCCESS;
					string command;
					switch( *s ) 
					{
					case 'e':
						lgCommandLineScript = true;
						input.clear();
						for( i++; i < argc; i++ )
						{
							string arg = argv[i];
							vector<string> word;
							Split(arg, ";", word, SPM_RELAX);
							if( word.size() > 1 )
							{
								for( size_t k=0; k < word.size(); ++k )
								{
									trimWhiteSpace(word[k]);
									command += word[k];
									if( !command.empty() )
										(void)cdRead(command);
									command.clear();
								}
							}
							else 
							{
								trimWhiteSpace(word[0]);
								command += word[0] + ' ';
							}
						}
						if( !command.empty() )
							(void)cdRead(command);
						break;
					case 'g':
					case 'p':
					case 'r':
						if( s[1] != '\0' )
						{
							prefix = s+1;
						}					
						else
						{
							if( ++i == argc || argv[i][0] == '-' )
							{
								fprintf( ioQQQ, "%s: no argument given for -%c flag\n",
									 argv[0], *s );
								cdEXIT(ES_FAILURE);
							}
							prefix = argv[i];
							if( strchr(prefix, cpu.i().chDirSeparator()[0]) != NULL )
							{
								fprintf( ioQQQ, "%s -%c %s: writing in subdirectories is not supported\n",
									 argv[0], *s, prefix );
								cdEXIT(ES_FAILURE);
							}
						}
						if( *s == 'g' )
							gprefix = prefix;
						else if( *s == 'p' )
						{
							pprefix = prefix;
							rprefix = prefix;
						}
						else if( *s == 'r' )
						{
							// make sure we erase the effects of a possible earlier -p flag
							pprefix = "";
							rprefix = prefix;
						}
						else
							TotalInsanity();
						s = NULL;
						lgFileIO = true;
						break;
					case 's':
						// this was already handled in main() -> skip it here
						if( s[1] == '\0' && i < argc-1 && argv[i+1][0] != '-' )
							++i;
						s = NULL;
						break;
					default:
						fprintf( ioQQQ, "%s: argument %d, `%s': flag -%c not understood\n",
							 argv[0], i, argv[i], *s );
						exit = ES_FAILURE;
						FALLTHROUGH;
					case 'h':
						fprintf( ioQQQ, "\nSupported flags are:\n\n" );
						fprintf( ioQQQ, "-p example\n" );
						fprintf( ioQQQ, "    Cloudy reads the input from the file example.in\n" );
						fprintf( ioQQQ, "    and writes the output to the file example.out.\n" );
						fprintf( ioQQQ, "    Additionally, all file names in SAVE commands are\n" );
						fprintf( ioQQQ, "    prepended with the string \"example\", e.g. the\n" );
						fprintf( ioQQQ, "    output of SAVE DR \".dr\" will be in example.dr.\n" );
						fprintf( ioQQQ, "-r example\n" );
						fprintf( ioQQQ, "    This does the same as the -p switch, except that\n" );
						fprintf( ioQQQ, "    the names used in SAVE commands are not altered.\n" );
						fprintf( ioQQQ, "-g example\n" );
						fprintf( ioQQQ, "    This switch is used internally in MPI grid runs.\n" );
						fprintf( ioQQQ, "    Normal users should not use this switch.\n" );
						fprintf( ioQQQ, "-s [hex number]\n" );
						fprintf( ioQQQ, "    Set the RNG seed to a specified fixed hex value.\n" );
						fprintf( ioQQQ, "    If no number is specified, a default value is used.\n" );
						fprintf( ioQQQ, "-e some command\n" );
						fprintf( ioQQQ, "-e 'some command ; another command'\n" );
						fprintf( ioQQQ, "    This will execute the command \"some command\".\n" );
						fprintf( ioQQQ, "    Multiple commands can be entered separated by ';'.\n" );
						fprintf( ioQQQ, "    No other flags can be entered after this flag.\n" );
						fprintf( ioQQQ, "-h\n" );
						fprintf( ioQQQ, "    Print this message.\n" );
						cdEXIT(exit);
					}
				}
			}
		}

		save.chGridPrefix = gprefix;
		save.chFilenamePrefix = pprefix;
		save.chRedirectPrefix = rprefix;

		/* following should be set true to print to file instead of std output */
		if( lgFileIO )
		{
			string Base = save.chGridPrefix + save.chRedirectPrefix;
			cdInput( Base + ".in" );
			if( cpu.i().MPIMode() == MS_POST_GRID )
			{
				// When here, we parse the input again to make sure that
				// certain variables are correctly set which are needed
				// when postprocessing the output. We are not interested
				// in the output, so we redirect that to a temp file that
				// will be discarded just before we return.
				outfile = new char[Base.length()+8];
				string OutName( Base + ".XXXXXX" );
				strcpy( outfile, OutName.c_str() );
				int fd = mkstemp( outfile );
				cdOutput( outfile, fdopen(fd, "w") );
			}
			else if( cpu.i().lgMPI_talk() )
			{
				cdOutput( Base + ".out" );
			}
			else
			{
				ostringstream oss;
				oss << ".err" << setfill('0') << setw(2) << cpu.i().nRANK();
				cdOutput( Base + oss.str() );
			}
		}

		if( optimize.nOptimiz == 0 && called.lgTalk && cpu.i().MPIMode() == MS_GRID )
			print_delimiter(-1);

		if( !lgCommandLineScript )
		{
			input.clear();
			/* keep reading input lines until end of file */
			string chLine;
			while( read_whole_line(chLine, ioStdin) )
			{
				bool lgReadingOutput = ( chLine.substr(0,25) == "                       * " );

				string chLocal;
				if( lgReadingOutput )
				{
					chLocal = chLine.substr(25);
					/* erase EOL character */
					size_t pp;
					if( (pp = chLocal.find_first_of("\n\r")) != string::npos )
						chLocal.erase(pp);
					trimTrailingWhiteSpace(chLocal);
					// remove final '*' from Cloudy output if present
					if( chLocal.back() == '*' )
						chLocal.pop_back();
				}
				else
				{
					chLocal = chLine;
				}

				if( lgInputEOF(chLocal) )
					break;
			
				/* stuff the command line into the internal stack */
				(void)cdRead(chLocal);
			}
		}

		ASSERT( input.curInclLevel == 0 );

		// optimize.lgVaryOn catches both optimizer and grid runs
		if( ( cpu.i().lgMPI() || optimize.lgVaryOn ) && save.chRedirectPrefix.empty() )
		{
			if( cpu.i().lgMaster() )
			{
				if( cpu.i().lgMPI() )
					fprintf( ioQQQ, " Please use the style \"mpirun -n np /path/to/cloudy.exe -r input\" when doing grid\n"
						" or optimizer runs.  See http://trac.nublado.org/wiki/RunCode for more information.\n" );
				else
					fprintf( ioQQQ, " Please use the style \"/path/to/cloudy.exe -r input\" when doing grid\n"
						" or optimizer runs.  See http://trac.nublado.org/wiki/RunCode for more information.\n" );
			}
			// stop the grid from being executed any further
			grid.lgGrid = false;
			cdEXIT(ES_FAILURE);
		}

		/* actually call the code.  This routine figures out whether the code will do
		 * a single model or be used to optimize on a spectrum, by looking for the
		 * keyword VARY on command lines.  It will call routine cloudy if no vary commands
		 * occur, and optimize_do if VARY does occur.  
		 * cdDrive returns 0 if calculation is ok, 1 if problems happened */
		if( cdDrive() )
			exit_status = ES_FAILURE;

		/* the last line of output will contain some interesting information about the model*/
		cdNwcns(
			/* the number of warnings, cautions, notes, and surprises */
			&NumberWarnings, 
			&NumberCautions, 
			&NumberNotes, 
			&NumberSurprises, 
			/* the number of temperature convergence failures */
			&NumberTempFailures, 
			/* the number of pressure convergence failures */
			&NumberPresFailures,
			/* the number of ionization convergence failures */
			&NumberIonFailures, 
			/* the number of electron density convergence failures */
			&NumberNeFailures );

		ostringstream finalMsg;

		finalMsg << " Cloudy ends: " << nzone << " zone";
		if( nzone > 1 )
			finalMsg << "s";

		finalMsg << ", " << iteration << " iteration";
		if( iteration > 1 )
			finalMsg << "s";

		if( NumberWarnings > 0 )
		{
			finalMsg << ", " << NumberWarnings << " warning";
			if( NumberWarnings > 1 )
				finalMsg << "s";
			/* this indicates error */
			exit_status = ES_FAILURE;
		}

		if( NumberCautions > 0 )
		{
			finalMsg << ", " << NumberCautions << " caution";
			if( NumberCautions > 1 )
				finalMsg << "s";
		}

		/* this flag was set in lgCheckMonitors*/
		if( !lgMonitorsOK )
		{
			finalMsg << ", ";
			/* some botches were three sigma */
			if( lgBigBotch  )
				finalMsg << "BIG ";
			finalMsg << "BOTCHED MONITORS!!!";
			/* this indicates error */
			exit_status = ES_FAILURE;
		}

		if( NumberTempFailures+NumberPresFailures+NumberIonFailures+NumberNeFailures > 0 )
		{
			finalMsg << ". Failures: " << NumberTempFailures << " thermal, ";
			finalMsg << NumberPresFailures << " pressure, ";
			finalMsg << NumberIonFailures << " ionization, ";
			finalMsg << NumberNeFailures << " electron density";
		}

		if( prt.lgPrintTime )
		{
			if( !cpu.i().lgMPI() && cpu.i().used_nCPU() == 1 )
			{
				finalMsg << ". (single thread)";
			}
			else if( !cpu.i().lgMPI() )
			{
				finalMsg << ". (" << cpu.i().used_nCPU() << " forked threads)";
			}
			else if( cpu.i().lgMPI() )
			{
				finalMsg << ". (rank " << cpu.i().nRANK()  << " of " << cpu.i().used_nCPU() << " MPI ranks)";
			}
			if (0)
				finalMsg << " Max memory used " << cdMemory() << "kB.";
			/* NB DO NOT CHANGE ANY ASPECT OF THE FOLLOWING STRINGS - THEY ARE USED TO RECORD
			 * EXEC TIME BY A PERL SCRIPT */
			/* print execution time [s] by default,
			 * need spaces around number so that logging perl script picks up correct number 
			 * ir_extime.pl script will delete through "ExecTime(s)" and remainder of line must be number */
			finalMsg << " ExecTime(s) " << fixed << setprecision(2) << cdExecTime();
		}
		else
		{
			finalMsg << ".";
		}

		if( called.lgTalk )
			fprintf( ioQQQ, "%s\n", finalMsg.str().c_str() );

		lgEarly_exit = false;

		/* cdDrive returned 1 if something bad happened, and 0 if everything is ok.  We will
		 * return 0 if everything is ok, and a non-zero error code if something bad happened.*/
		cdEXIT(exit_status);
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
		if( called.lgTalk )
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

		if ( called.lgTalk && prt.lgPrintHTML )
		{
			fprintf( ioQQQ,"</pre>\n");
			fprintf( ioQQQ,"</body>\n");
			fprintf( ioQQQ,"</html>\n");
		}
		
		exit_status = e.exit_status();
		if( exit_status == ES_FAILURE && !lgEarly_exit )
		{
			// try to make the error code more descriptive
			// when there is no early exit in the code, then these 3
			// should be the only reasons for a non-zero exit code
			if( NumberWarnings > 0 )
				exit_status = ES_WARNINGS;
			if( !lgMonitorsOK )
				exit_status = ES_BOTCHES;
		}
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

	// print backtrace, if there is any...
	cpu.i().PrintBacktrace( "" );

	if( called.lgTalk && cpu.i().MPIMode() == MS_GRID )
		print_delimiter(optimize.nOptimiz);

	cdPrepareExit(exit_status);

	if( outfile != NULL )
	{
		remove( outfile );
		delete[] outfile;
	}

	return exit_status;
}
