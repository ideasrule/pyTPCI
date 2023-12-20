/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*optimize_phymir Peter van Hoof's optimization routine */

#include "cddefines.h"
#include "version.h"
#include "optimize.h"
#include "service.h"
#include "ran.h"
#if defined(__unix) || defined(__APPLE__)
#include <cstddef>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#else
#define pid_t int
#define fork() TotalInsanityAsStub<pid_t>()
#define wait(X) TotalInsanityAsStub<pid_t>()
#endif


/** this is the name of the file that optimize_phymir automatically creates containing
 * information to continue the optimization at a later time.  */
const char* STATEFILE = "continue.pmr";
const char* STATEFILE_BACKUP = "continue.bak";

/* The optimization algorithm Phymir and its subsidiary routines have been
 * written by Peter van Hoof. */


template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_clear1()
{
	DEBUG_ENTRY( "p_clear1()" );

	p_xmax = numeric_limits<X>::max();
	p_ymax = numeric_limits<Y>::max() / Y(2.);
	for( int i=0; i < 2*NP+1; ++i )
	{
		for( int j=0; j < NP; ++j )
			p_xp[i][j] = -p_xmax;
		p_yp[i] = -p_ymax;
	}
	for( int i=0; i < NP; ++i )
	{
		p_absmin[i] = -p_xmax;
		p_absmax[i] = p_xmax;
		p_varmin[i] = p_xmax;
		p_varmax[i] = -p_xmax;
		for( int j=0; j < NP; ++j )
			p_a2[i][j] = -p_xmax;
		p_c1[i] = -p_xmax;
		p_c2[i] = -p_xmax;
		p_xc[i] = -p_xmax;
		p_xcold[i] = -p_xmax;
	}
	p_vers = VRSNEW;
	p_toler = -p_xmax;
	p_dmax = X(0.);
	p_dold = X(0.);
	p_ymin = p_ymax;
	p_dim = int32(NP);
	p_sdim = int32(NSTR);
	p_nvar = int32(0);
	p_noptim = int32(0);
	p_maxiter = int32(0);
	p_jmin = int32(-1);
	p_maxcpu = int32(1);
	p_curcpu = int32(0);
	p_mode = PHYMIR_ILL;
	memset( &p_chState, 0, size_t(NSTR) );
	memset( &p_chStr1, 0, size_t(NSTR) );
	memset( &p_chStr2, 0, size_t(NSTR) );
	memset( &p_chStr3, 0, size_t(NSTR) );
	p_func = NULL;
}

template<class T>
inline bool WR_ITEM(const T& x, FILE* io)
{
	return ( fwrite(&x, sizeof(T), 1, io) != 1 );
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_wr_state(const char *fnam) const
{
	DEBUG_ENTRY( "p_wr_state()" );

	if( cpu.i().lgMaster() && strlen(fnam) > 0 )
	{
		FILE *fdes = open_data( fnam, "wb", AS_TRY );
		bool lgErr = ( fdes == NULL );
		lgErr = lgErr || WR_ITEM(p_xmax, fdes);
		lgErr = lgErr || WR_ITEM(p_ymax, fdes);
		lgErr = lgErr || WR_ITEM(p_xp, fdes);
		lgErr = lgErr || WR_ITEM(p_yp, fdes);
		lgErr = lgErr || WR_ITEM(p_absmin, fdes);
		lgErr = lgErr || WR_ITEM(p_absmax, fdes);
		lgErr = lgErr || WR_ITEM(p_varmin, fdes);
		lgErr = lgErr || WR_ITEM(p_varmax, fdes);
		lgErr = lgErr || WR_ITEM(p_a2, fdes);
		lgErr = lgErr || WR_ITEM(p_c1, fdes);
		lgErr = lgErr || WR_ITEM(p_c2, fdes);
		lgErr = lgErr || WR_ITEM(p_xc, fdes);
		lgErr = lgErr || WR_ITEM(p_xcold, fdes);
		lgErr = lgErr || WR_ITEM(p_vers, fdes);
		lgErr = lgErr || WR_ITEM(p_toler, fdes);
		lgErr = lgErr || WR_ITEM(p_dmax, fdes);
		lgErr = lgErr || WR_ITEM(p_dold, fdes);
		lgErr = lgErr || WR_ITEM(p_ymin, fdes);
		lgErr = lgErr || WR_ITEM(p_dim, fdes);
		lgErr = lgErr || WR_ITEM(p_sdim, fdes);
		lgErr = lgErr || WR_ITEM(p_nvar, fdes);
		lgErr = lgErr || WR_ITEM(p_noptim, fdes);
		lgErr = lgErr || WR_ITEM(p_maxiter, fdes);
		lgErr = lgErr || WR_ITEM(p_jmin, fdes);
		lgErr = lgErr || WR_ITEM(p_maxcpu, fdes);
		lgErr = lgErr || WR_ITEM(p_curcpu, fdes);
		lgErr = lgErr || WR_ITEM(p_mode, fdes);
		lgErr = lgErr || WR_ITEM(p_chState, fdes);
		lgErr = lgErr || WR_ITEM(p_chStr1, fdes);
		lgErr = lgErr || WR_ITEM(p_chStr2, fdes);
		lgErr = lgErr || WR_ITEM(p_chStr3, fdes);
		lgErr = lgErr || ( fclose(fdes) != 0 );
		if( lgErr )
		{
			printf( "p_wr_state: error writing file: %s\n", fnam );
			remove(fnam);
		}
	}
	return;
}

template<class T>
inline bool RD_ITEM(T& x, FILE* io)
{
	return ( fread(&x, sizeof(T), 1, io) != 1 );
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_rd_state(const char *fnam)
{
	DEBUG_ENTRY( "p_rd_state()" );

	FILE *fdes = open_data( fnam, "rb", AS_LOCAL_ONLY );
	bool lgErr = false;
	lgErr = lgErr || RD_ITEM(p_xmax, fdes);
	lgErr = lgErr || RD_ITEM(p_ymax, fdes);
	lgErr = lgErr || RD_ITEM(p_xp, fdes);
	lgErr = lgErr || RD_ITEM(p_yp, fdes);
	lgErr = lgErr || RD_ITEM(p_absmin, fdes);
	lgErr = lgErr || RD_ITEM(p_absmax, fdes);
	lgErr = lgErr || RD_ITEM(p_varmin, fdes);
	lgErr = lgErr || RD_ITEM(p_varmax, fdes);
	lgErr = lgErr || RD_ITEM(p_a2, fdes);
	lgErr = lgErr || RD_ITEM(p_c1, fdes);
	lgErr = lgErr || RD_ITEM(p_c2, fdes);
	lgErr = lgErr || RD_ITEM(p_xc, fdes);
	lgErr = lgErr || RD_ITEM(p_xcold, fdes);
	lgErr = lgErr || RD_ITEM(p_vers, fdes);
	lgErr = lgErr || RD_ITEM(p_toler, fdes);
	lgErr = lgErr || RD_ITEM(p_dmax, fdes);
	lgErr = lgErr || RD_ITEM(p_dold, fdes);
	lgErr = lgErr || RD_ITEM(p_ymin, fdes);
	lgErr = lgErr || RD_ITEM(p_dim, fdes);
	lgErr = lgErr || RD_ITEM(p_sdim, fdes);
	lgErr = lgErr || RD_ITEM(p_nvar, fdes);
	lgErr = lgErr || RD_ITEM(p_noptim, fdes);
	lgErr = lgErr || RD_ITEM(p_maxiter, fdes);
	lgErr = lgErr || RD_ITEM(p_jmin, fdes);
	lgErr = lgErr || RD_ITEM(p_maxcpu, fdes);
	lgErr = lgErr || RD_ITEM(p_curcpu, fdes);
	lgErr = lgErr || RD_ITEM(p_mode, fdes);
	lgErr = lgErr || RD_ITEM(p_chState, fdes);
	lgErr = lgErr || RD_ITEM(p_chStr1, fdes);
	lgErr = lgErr || RD_ITEM(p_chStr2, fdes);
	lgErr = lgErr || RD_ITEM(p_chStr3, fdes);
	lgErr = lgErr || ( fclose(fdes) != 0 );
	if( lgErr )
	{
		printf( "p_rd_state: error reading file: %s\n", fnam );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

template<class X, class Y, int NP, int NSTR>
Y phymir_state<X,Y,NP,NSTR>::p_execute_job(const X x[],
					   int jj,
					   int runNr)
{
	DEBUG_ENTRY( "p_execute_job()" );

	pid_t pid;
	switch( p_mode )
	{
	case PHYMIR_SEQ:
		return p_lgLimitExceeded(x) ? p_ymax : p_func(x,runNr);
	case PHYMIR_FORK:
		p_curcpu++;
		if( p_curcpu > p_maxcpu )
		{
			/* too many processes, wait for one to finish */
			(void)wait(NULL);
			p_curcpu--;
		}
		/* flush all open files */
		fflush(NULL);
		pid = fork();
		if( pid < 0 )
		{
			fprintf( ioQQQ, "creating the child process failed\n" );
			cdEXIT(EXIT_FAILURE);
		}
		else if( pid == 0 )
		{
			/* make sure each child has a unique random number sequence */
			ran.new_rank( p_noptim+1 );
			/* this is child process so execute */
			p_execute_job_parallel( x, jj, runNr );
			/* at this point ioQQQ points to the main output of the parent process,
			 * the child should not close that in cdEXIT(), so wipe the handle here */
			ioQQQ = NULL;
			cdEXIT(EXIT_SUCCESS);
		}
		// the real y-value is not available yet...
		return p_ymax;
	case PHYMIR_MPI:
		/* make sure behavior is the same as for fork runs */
		ran.new_rank( p_noptim+1 );
		if( (jj%cpu.i().nCPU()) == cpu.i().nRANK() )
			p_execute_job_parallel( x, jj, runNr );
		// the real y-value is not available yet...
		return p_ymax;
	default:
		TotalInsanity();
	}
}

// p_execute_job_parallel MUST be const, otherwise changes by child processes may be lost!
template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_execute_job_parallel(const X x[],
						       int jj,
						       int runNr) const
{
	DEBUG_ENTRY( "p_execute_job_parallel()" );

	char fnam1[20], fnam2[20];
	sprintf(fnam1,"yval_%d",jj);
	sprintf(fnam2,"output_%d",jj);
	/* each child should have its own output file */
	FILE *ioQQQ_old = ioQQQ;
	ioQQQ = open_data( fnam2, "w" );
	// fail-safe in case p_func crashes without being caught...
	Y yval = p_ymax;
	wr_block(&yval,sizeof(yval),fnam1);
	if( !p_lgLimitExceeded(x) )
	{
		try
		{
			yval = p_func(x,runNr);
		}
		catch( ... )
		{
			// make sure output is complete since files may not have been closed properly...
			fflush(NULL);
			yval = p_ymax;
		}
		wr_block(&yval,sizeof(yval),fnam1);
	}
	fclose( ioQQQ );
	ioQQQ = ioQQQ_old;
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_barrier(int jlo,
					  int jhi)
{
	DEBUG_ENTRY( "p_barrier()" );

	switch( p_mode )
	{
	case PHYMIR_SEQ:
		// nothing to do...
		break;
	case PHYMIR_FORK:
		/* wait for child processes to terminate */
		while( p_curcpu > 0 )
		{
			(void)wait(NULL);
			p_curcpu--;
		}
		p_process_output( jlo, jhi );
		break;
	case PHYMIR_MPI:
		// wait for all function evaluations to finish so that output is complete
		MPI_Barrier( MPI_COMM_WORLD );
		p_process_output( jlo, jhi );
		// y values are known now, so broadcast to other ranks
		MPI_Bcast( &p_yp[jlo], jhi-jlo+1, MPI_type(p_yp), 0, MPI_COMM_WORLD );
		break;
	default:
		TotalInsanity();
	}
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_process_output(int jlo,
						 int jhi)
{
	DEBUG_ENTRY( "p_process_output()" );

	if( cpu.i().lgMaster() )
	{
		char fnam[20];
		for( int jj=jlo; jj <= jhi; jj++ )
		{
			sprintf(fnam,"yval_%d",jj);
			rd_block(&p_yp[jj],sizeof(p_yp[jj]),fnam);
			remove(fnam);
			sprintf(fnam,"output_%d",jj);
			append_file(ioQQQ,fnam);
			remove(fnam);
		}
	}
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_evaluate_hyperblock()
{
	DEBUG_ENTRY( "p_evaluate_hyperblock()" );

	int jlo = 1, jhi = 0;
	for( int j=0; j < p_nvar; j++ )
	{
		X sgn = X(1.);
		for( int jj=2*j+1; jj <= 2*j+2; jj++ )
		{
			sgn = -sgn;
			for( int i=0; i < p_nvar; i++ )
			{
				p_xp[jj][i] = p_xc[i] + sgn*p_dmax*p_c2[j]*p_a2[j][i];
				p_varmax[i] = max(p_varmax[i],p_xp[jj][i]);
				p_varmin[i] = min(p_varmin[i],p_xp[jj][i]);
			}
			if( !lgMaxIterExceeded() )
			{
				// p_execute_job will not return the correct chi^2 in parallel mode
				// only after p_barrier() has finished will p_yp[jj] contain the correct value
				p_yp[jj] = p_execute_job( p_xp[jj], jj, p_noptim++ );
				jhi = jj;
			}
		}
	}
	/* wait for jobs to terminate */
	p_barrier( jlo, jhi );
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_setup_next_hyperblock()
{
	DEBUG_ENTRY( "p_setup_next_hyperblock()" );

	const Y F0 = Y(1.4142136);
	const X F1 = X(0.7071068);
	const X F2 = X(0.1);

	/* find best model */
	for( int jj=1; jj <= 2*p_nvar; jj++ )
	{
		if( p_yp[jj] < p_ymin )
		{
			p_ymin = p_yp[jj];
			p_jmin = jj;
		}
	}
	bool lgNewCnt = p_jmin > 0;

	/* determine minimum and relative uncertainties */
	bool lgNegd2 = false;
	X xnrm = X(0.);
	X xhlp[NP];
	for( int i=0; i < p_nvar; i++ )
	{
		Y d1 = p_yp[2*i+2] - p_yp[2*i+1];
		Y d2 = Y(0.5)*p_yp[2*i+2] - p_yp[0] + Y(0.5)*p_yp[2*i+1];
		if( d2 <= Y(0.) )
			lgNegd2 = true;
		xhlp[i] = -p_dmax*p_c2[i]*(static_cast<X>(max(min((Y(0.25)*d1)/max(d2,Y(1.e-10)),F0),-F0)) -
					   p_delta(2*i+1,p_jmin) + p_delta(2*i+2,p_jmin));
		xnrm += pow2(xhlp[i]);
	}
	xnrm = static_cast<X>(sqrt(xnrm));
	/* set up new transformation matrix */
	int imax = 0;
	X amax = X(0.);
	for( int j=0; j < p_nvar; j++ )
	{
		for( int i=0; i < p_nvar; i++ )
		{
			if( xnrm > X(0.) )
			{
				if( j == 0 )
				{
					p_a2[0][i] *= xhlp[0];
				}
				else
				{
					p_a2[0][i] += xhlp[j]*p_a2[j][i];
					p_a2[j][i] = p_delta(i,j);
					if( j == p_nvar-1 && abs(p_a2[0][i]) > amax )
					{
						imax = i;
						amax = abs(p_a2[0][i]);
					}
				}
			}
			else
			{
				p_a2[j][i] = p_delta(i,j);
			}
		}
	}
	/* this is to assure maximum linear independence of the base vectors */
	if( imax > 0 )
	{
		p_a2[imax][0] = X(1.);
		p_a2[imax][imax] = X(0.);
	}
	/* apply Gram-Schmidt procedure to get orthogonal base */
	p_phygrm( p_a2, p_nvar );

	/* set up hyperblock dimensions in new base and choose new center */
	for( int i=0; i < p_nvar; i++ )
	{
		p_c2[i] = X(0.);
		for( int j=0; j < p_nvar; j++ )
		{
			p_c2[i] += pow2(p_a2[i][j]/p_c1[j]);
		}
		p_c2[i] = static_cast<X>(1./sqrt(p_c2[i]));
		p_xc[i] = p_xp[p_jmin][i];
		p_xp[0][i] = p_xc[i];
	}
	p_yp[0] = p_yp[p_jmin];
	p_jmin = 0;
	/* choose size of next hyperblock */
	X r1, r2;
	if( lgNegd2 )
	{
		/* this means that the hyperblock either is bigger than the scale
		 * on which the function varies, or is so small that we see noise.
		 * in both cases make the hyperblock smaller. */
		r1 = F1;
		r2 = F1;
	}
	else
	{
		r1 = F2;
		if( lgNewCnt )
		{
			/* when we make progress, p_dmax may become bigger */
			r2 = sqrt(1.f/F1);
		}
		else
		{
			/* when there is no progress force p_dmax smaller, otherwise there
			 * may never be an end */
			r2 = sqrt(F1);
		}
	}
	p_dmax = min(max(xnrm/p_c2[0],p_dmax*r1),p_dmax*r2);
	/* fail-safe against excessive behaviour */
	p_dmax = MIN2(p_dmax,p_dold);
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_reset_hyperblock()
{
	DEBUG_ENTRY( "p_reset_hyperblock()" );

	if( !lgConvergedRestart() )
	{
		/* reset hyperblock so that we can restart the optimization */
		for( int i=0; i < p_nvar; i++ )
		{
			p_xcold[i] = p_xc[i];
			p_c2[i] = p_c1[i];
		}
		p_dmax = p_dold;
		p_reset_transformation_matrix();
	}
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_phygrm(X a[][NP],
					 int n)
{
	DEBUG_ENTRY( "p_phygrm()" );

	/* apply modified Gram-Schmidt procedure to a */
	for( int i=0; i < n; i++ )
	{
		X ip = X(0.);
		for( int k=0; k < n; k++ )
			ip += pow2(a[i][k]);
		ip = sqrt(ip);
		for( int k=0; k < n; k++ )
			a[i][k] /= ip;
		for( int j=i+1; j < n; j++ )
		{
			X ip = X(0.);
			for( int k=0; k < n; k++ )
				ip += a[i][k]*a[j][k];
			for( int k=0; k < n; k++ )
				a[j][k] -= ip*a[i][k];
		}
	}
	return;
}

template<class X, class Y, int NP, int NSTR>
bool phymir_state<X,Y,NP,NSTR>::p_lgLimitExceeded(const X x[]) const
{
	DEBUG_ENTRY( "p_lgLimitExceeded()" );

	for( int i=0; i < p_nvar; i++ )
	{
		if( x[i] < p_absmin[i] )
			return true;
		if( x[i] > p_absmax[i] )
			return true;
	}
	return false;
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::p_reset_transformation_matrix()
{
	DEBUG_ENTRY( "p_reset_transformation_matrix()" );

	/* initialize transformation matrix to unity */
	for( int i=0; i < p_nvar; i++ )
		for( int j=0; j < p_nvar; j++ )
			p_a2[j][i] = p_delta(i,j);
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::init_minmax(const X pmin[], // pmin[nvar]: minimum values for params
					    const X pmax[], // pmax[nvar]: maximum values for params
					    int nvar)       // will not be initialized yet
{
	DEBUG_ENTRY( "init_minmax()" );

	ASSERT( !lgInitialized() );

	for( int i=0; i < nvar; i++ )
	{
		p_absmin[i] = pmin[i];
		p_absmax[i] = pmax[i];
	}
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::init_strings(const string& date,    // date string for inclusion in state file
					     const string& version, // version string for inclusion in state file
					     const char* host_name) // host name for inclusion in state file
{
	DEBUG_ENTRY( "init_strings()" );

	strncpy( p_chStr1, date.c_str(), NSTR );
	p_chStr1[NSTR-1] = '\0';
	strncpy( p_chStr2, version.c_str(), NSTR );
	p_chStr2[NSTR-1] = '\0';
	if( host_name != NULL )
		strncpy( p_chStr3, host_name, NSTR );
	p_chStr3[NSTR-1] = '\0';
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::init_state_file_name(const char* fnam) // name of the state file that will be written
{
	DEBUG_ENTRY( "init_state_file_name()" );

	// use NSTR-1 so that last 0 byte is not overwritten...
	strncpy( p_chState, fnam, NSTR-1 );
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::initial_run(Y (*func)(const X[],int),  // function to be optimized
					    int nvar,                  // number of parameters to be optimized
					    const X start[],           // start[n]: initial estimates for params
					    const X del[],             // del[n]: initial stepsizes for params
					    X toler,                   // absolute tolerance on parameters
					    int maxiter,               // maximum number of iterations
					    phymir_mode mode,          // execution mode: sequential, fork, mpi
					    int maxcpu)                // maximum no. of CPUs to use
{
	DEBUG_ENTRY( "initial_run()" );

	ASSERT( nvar > 0 && nvar <= NP );

	p_func = func;
	p_nvar = nvar;
	p_toler = toler;
	p_maxiter = maxiter;
	p_mode = mode;
	p_maxcpu = maxcpu;
	p_noptim = 0;

	/* initialize hyperblock dimensions and center */
	p_dmax = X(0.);
	for( int i=0; i < p_nvar; i++ )
		p_dmax = max(p_dmax,abs(del[i]));

	p_dold = p_dmax;
	for( int i=0; i < p_nvar; i++ )
	{
		p_xc[i] = start[i];
		p_xcold[i] = p_xc[i] + X(10.)*p_toler;
		p_c1[i] = abs(del[i])/p_dmax;
		p_c2[i] = p_c1[i];
		p_xp[0][i] = p_xc[i];
		p_varmax[i] = max(p_varmax[i],p_xc[i]);
		p_varmin[i] = min(p_varmin[i],p_xc[i]);
	}

	// p_execute_job will not return the correct chi^2 in parallel mode
	// only after p_barrier() has finished will p_yp[0] contain the correct value
	p_yp[0] = p_execute_job( p_xc, 0, p_noptim++ );
	p_barrier( 0, 0 );

	p_ymin = p_yp[0];
	p_jmin = 0;

	p_reset_transformation_matrix();

	p_wr_state( p_chState );
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::continue_from_state(Y (*func)(const X[],int), // function to be optimized	
						    int nvar,                 // number of parameters to be optimized
						    const char* fnam,         // file name of the state file
						    X toler,                  // absolute tolerance on parameters
						    int maxiter,              // maximum number of iterations
						    phymir_mode mode,         // execution mode: sequential, fork, mpi
						    int maxcpu)               // maximum no. of CPUs to use
{
	DEBUG_ENTRY( "continue_from_state()" );

	p_rd_state( fnam );
	// sanity checks
	if( !fp_equal( p_vers, VRSNEW ) )
	{
		printf( "optimize continue - file has incompatible version, sorry\n" );
		cdEXIT(EXIT_FAILURE);
	}
	if( p_dim != NP )
	{
		printf( "optimize continue - arrays have wrong dimension, sorry\n" );
		cdEXIT(EXIT_FAILURE);
	}
	if( p_sdim != NSTR )
	{
		printf( "optimize continue - strings have wrong length, sorry\n" );
		cdEXIT(EXIT_FAILURE);
	}
	if( p_nvar != nvar )
	{
		printf( "optimize continue - wrong number of free parameters, sorry\n" );
		cdEXIT(EXIT_FAILURE);
	}
	// this pointer was not part of the state file, it would be useless...
	p_func = func;
	// Option to override the tolerance set in the state file. This is useful
	// if you want to refine an already finished run to a smaller tolerance.
	p_toler = toler;
	// you ran out of iterations, but wanted to continue...
	p_maxiter = maxiter;
	// it is OK to continue in a different mode
	p_mode = mode;
	p_maxcpu = maxcpu;
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::optimize()
{
	DEBUG_ENTRY( "optimize()" );

	ASSERT( lgInitialized() );

	while( !lgConverged() )
	{
		p_evaluate_hyperblock();
		if( lgMaxIterExceeded() )
			break;
		p_setup_next_hyperblock();
		p_wr_state( p_chState );
	}
}

template<class X, class Y, int NP, int NSTR>
void phymir_state<X,Y,NP,NSTR>::optimize_with_restart()
{
	DEBUG_ENTRY( "optimize_with_restart()" );

	ASSERT( lgInitialized() );

	while( !lgConvergedRestart() )
	{
		optimize();
		if( lgMaxIterExceeded() )
			break;
		p_reset_hyperblock();
	}
}

template<class X, class Y, int NP, int NSTR>
bool phymir_state<X,Y,NP,NSTR>::lgConvergedRestart() const
{
	DEBUG_ENTRY( "lgConvergedRestart()" );

	if( lgConverged() )
	{
		X dist = X(0.);
		for( int i=0; i < p_nvar; i++ )
			dist += pow2(p_xc[i] - p_xcold[i]);
		dist = static_cast<X>(sqrt(dist));
		return ( dist <= p_toler );
	}
	return false;
}

void optimize_phymir(realnum xc[], 
		     const realnum del[], 
		     long int nvarPhymir, 
		     chi2_type *ymin, 
		     realnum toler)
{
	DEBUG_ENTRY( "optimize_phymir()" );

	if( nvarPhymir > LIMPAR )
	{
		fprintf( ioQQQ, "optimize_phymir: too many parameters are varied, increase LIMPAR\n" );
		cdEXIT(EXIT_FAILURE);
	}

	phymir_state<realnum,chi2_type,LIMPAR,STDLEN> phymir;

	// first check if a statefile already exists, back it up in that case
	(void)remove(STATEFILE_BACKUP);
	FILE *tmp = open_data( STATEFILE, "r", AS_LOCAL_ONLY_TRY );
	if( tmp != NULL )
	{
		fclose( tmp );
		// create backup copy of statefile
		FILE *dest = open_data( STATEFILE_BACKUP, "wb", AS_TRY );
		if( dest != NULL )
		{
			append_file( dest, STATEFILE );
			fclose( dest );
		}
	}

	phymir_mode mode = optimize.lgParallel ? ( cpu.i().lgMPI() ? PHYMIR_MPI : PHYMIR_FORK ) : PHYMIR_SEQ;
	long nCPU = optimize.lgParallel ? ( cpu.i().lgMPI() ? cpu.i().nCPU() : optimize.useCPU ) : 1;
	cpu.i().set_used_nCPU( nCPU );
	if( optimize.lgOptCont )
	{
		phymir.continue_from_state( optimize_func, nvarPhymir, STATEFILE, toler,
					    optimize.nIterOptim, mode, nCPU );
	}
	else
	{
		phymir.init_state_file_name( STATEFILE );
		phymir.init_strings( t_version::Inst().chDate, t_version::Inst().chVersion,
				     cpu.i().host_name() );
		phymir.initial_run( optimize_func, nvarPhymir, xc, del, toler,
				    optimize.nIterOptim, mode, nCPU );
	}

	phymir.optimize_with_restart();

	if( phymir.lgMaxIterExceeded() )
	{
		fprintf( ioQQQ, " Optimizer exceeding maximum iterations.\n" );
		fprintf( ioQQQ, " This can be reset with the OPTIMIZE ITERATIONS command.\n" );
	}

	// updates to optimize.nOptimiz and optimize.varmin/max in child processes are lost, so repair here...
	optimize.nOptimiz = phymir.noptim();
	for( int i=0; i < nvarPhymir; i++ )
	{
		xc[i] = phymir.xval(i);
		optimize.varmax[i] = min(phymir.xmax(i),optimize.varang[i][1]);
		optimize.varmin[i] = max(phymir.xmin(i),optimize.varang[i][0]);
	}
	*ymin = phymir.yval();

	// this assures that MPI and fork runs generate the same random numbers in the final model...
	if( optimize.lgParallel && cpu.i().lgMaster() )
		ran.new_rank(optimize.nOptimiz+1);

	return;
}
