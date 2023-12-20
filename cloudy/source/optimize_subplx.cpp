/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "optimize.h"

STATIC void calcc(long,realnum*,long,long,int,realnum[]);
STATIC double cdasum(long,realnum[],long);
STATIC void cdaxpy(long,double,realnum[],long,realnum[],long);
STATIC void cdcopy(long,realnum[],long,realnum[],long);
STATIC void csscal(long,double,realnum[],long);
STATIC double dist(long,realnum[],realnum[]);
STATIC void evalf(long,long[],realnum[],long,realnum[],realnum*,long*);
STATIC void fstats(double,long,int);
STATIC void newpt(long,double,realnum[],realnum[],int,realnum[],int*);
STATIC void order(long,realnum[],long*,long*,long*);
STATIC void partx(long,long[],realnum[],long*,long[]);
STATIC void setstp(long,long,realnum[],realnum[]);
STATIC void simplx(long,realnum[],long,long[],long,int,realnum[],
	realnum*,long*,realnum*,realnum[],long*);
STATIC void sortd(long,realnum[],long[]);
STATIC void start(long,realnum[],realnum[],long,long[],realnum*,int*);
STATIC void subopt(long);

/*lint -e801 goto is deprecated */
/*lint -e64 type mismatch */
/*********************************************************************
 *
 *     NB this is much the original coding and does not use strong typing
 *gjf mod: to avoid conflict with exemplar parallel version of lapack,
 *dasum changed to cdasum
 *daxpy changed to cdaxpy
 *dcopy changed to cdcopy
 *sscal changed to csscal
 *
 *********************************************************************** */
struct t_usubc {
	realnum alpha, 
	  beta, 
	  gamma, 
	  delta, 
	  psi, 
	  omega;
	long int nsmin, 
	  nsmax, 
	  irepl, 
	  ifxsw;
	realnum bonus, 
	  fstop;
	long int nfstop, 
	  nfxe;
	realnum fxstat[4], 
	  ftest;
	int minf, 
	  initx, 
	  newx;
	}	usubc;
struct t_isubc {
	realnum fbonus, 
	  sfstop, 
	  sfbest;
	int IntNew;
	}	isubc;

void optimize_subplex(long int n, 
	  double tol, 
	  long int maxnfe, 
	  long int mode, 
	  realnum scale[], 
	  realnum x[], 
	  realnum *fx, 
	  long int *nfe, 
	  realnum work[], 
	  long int iwork[], 
	  long int *iflag)
{
	static int cmode;
	static long int i, 
	  i_, 
	  ifsptr, 
	  ins, 
	  insfnl, 
	  insptr, 
	  ipptr, 
	  isptr,
	  istep, 
	  istptr, 
	  nnn, 
	  ns, 
	  nsubs;
	static realnum dum,
	  scl[1], 
	  sfx, 
	  xpscl, 
	  xstop,
	  xstop2;

	static const realnum bnsfac[2][3] = { {-1.0f,-2.0f,0.0f}, {1.0f,0.0f,2.0f} };

	DEBUG_ENTRY( "optimize_subplex()" );

	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * Jason Ferguson:
	 * Changed on 8/4/94, in order to incorparate into cloudy
	 * changes made are: double precision to real,
	 * a 2-D function f(n,x) into 1-D f(x),
	 * and the termination check.
	 *
	 * optimize_subplex uses the subplex method to solve unconstrained
	 * optimization problems.  The method is well suited for
	 * optimizing objective functions that are noisy or are
	 * discontinuous at the solution.
	 *
	 * optimize_subplex sets default optimization options by calling the
	 * subroutine subopt.  The user can override these defaults
	 * by calling subopt prior to calling optimize_subplex, changing the
	 * appropriate common variables, and setting the value of
	 * mode as indicated below.
	 *
	 * By default, optimize_subplex performs minimization.
	 *
	 * input
	 *
	 *   ffff   - user supplied function f(n,x) to be optimized,
	 *            declared external in calling routine
	 *   always uses optimize_func - change 97 dec 8
	 *
	 *   n      - problem dimension
	 *
	 *   tol    - relative error tolerance for x (tol .ge. 0.)
	 *
	 *   maxnfe - maximum number of function evaluations
	 *
	 *   mode   - integer mode switch with binary expansion
	 *            (bit 1) (bit 0) :
	 *            bit 0 = 0 : first call to optimize_subplex
	 *                  = 1 : continuation of previous call
	 *            bit 1 = 0 : use default options
	 *                  = 1 : user set options
	 *
	 *   scale  - scale and initial stepsizes for corresponding
	 *            components of x
	 *            (If scale(1) .lt. 0.,
	 *            ABS(scale(1)) is used for all components of x,
	 *            and scale(2),...,scale(n) are not referenced.)
	 *
	 *   x      - starting guess for optimum
	 *
	 *   work   - real work array of dimension .ge.
	 *            2*n + nsmax*(nsmax+4) + 1
	 *            (nsmax is set in subroutine subopt.
	 *            default: nsmax = MIN(5,n))
	 *
	 *   iwork  - integer work array of dimension .ge.
	 *            n + INT(n/nsmin)
	 *            (nsmin is set in subroutine subopt.
	 *            default: nsmin = MIN(2,n))
	 *
	 * output
	 *
	 *   x      - computed optimum
	 *
	 *   fx     - value of f at x
	 *
	 *   nfe    - number of function evaluations
	 *
	 *   iflag  - error flag
	 *            = -2 : invalid input
	 *            = -1 : maxnfe exceeded
	 *            =  0 : tol satisfied
	 *            =  1 : limit of machine precision
	 *            =  2 : fstop reached (fstop usage is determined
	 *                   by values of options minf, nfstop, and
	 *                   irepl. default: f(x) not tested against
	 *                   fstop)
	 *            iflag should not be reset between calls to
	 *            optimize_subplex.
	 *
	 * common
	 * */

	if( (mode%2) == 0 )
	{

		/*       first call, check input
		 * */
		if( n < 1 )
			goto L_120;
		if( tol < 0.0f )
			goto L_120;
		if( maxnfe < 1 )
			goto L_120;
		if( scale[0] >= 0.0f )
		{
			for( i=1; i <= n; i++ )
			{
				i_ = i - 1;
				xpscl = x[i_] + scale[i_];
				if( fp_equal( xpscl, x[i_] ) )
					goto L_120;
			}
		}
		else
		{
			scl[0] = (realnum)fabs(scale[0]);
			for( i=1; i <= n; i++ )
			{
				i_ = i - 1;
				xpscl = x[i_] + scl[0];
				if( fp_equal( xpscl, x[i_] ) )
					goto L_120;
			}
		}
		if( ((mode/2)%2) == 0 )
		{
			subopt(n);
		}
		else if( usubc.alpha <= 0.0f )
		{
			goto L_120;
		}
		else if( usubc.beta <= 0.0f || usubc.beta >= 1.0f )
		{
			goto L_120;
		}
		else if( usubc.gamma <= 1.0f )
		{
			goto L_120;
		}
		else if( usubc.delta <= 0.0f || usubc.delta >= 1.0f )
		{
			goto L_120;
		}
		else if( usubc.psi <= 0.0f || usubc.psi >= 1.0f )
		{
			goto L_120;
		}
		else if( usubc.omega <= 0.0f || usubc.omega >= 1.0f )
		{
			goto L_120;
		}
		else if( (usubc.nsmin < 1 || usubc.nsmax < usubc.nsmin) || 
		  n < usubc.nsmax )
		{
			goto L_120;
		}
		else if( n < ((n - 1)/usubc.nsmax + 1)*usubc.nsmin )
		{
			goto L_120;
		}
		else if( usubc.irepl < 0 || usubc.irepl > 2 )
		{
			goto L_120;
		}
		else if( usubc.ifxsw < 1 || usubc.ifxsw > 3 )
		{
			goto L_120;
		}
		else if( usubc.bonus < 0.0f )
		{
			goto L_120;
		}
		else if( usubc.nfstop < 0 )
		{
			goto L_120;
		}

		/*       initialization
		 * */
		istptr = n + 1;
		isptr = istptr + n;
		ifsptr = isptr + usubc.nsmax*(usubc.nsmax + 3);
		insptr = n + 1;
		if( scale[0] > 0.0f )
		{
			cdcopy(n,scale,1,work,1);
			cdcopy(n,scale,1,&work[istptr-1],1);
		}
		else
		{
			cdcopy(n,scl,0,work,1);
			cdcopy(n,scl,0,&work[istptr-1],1);
		}
		for( i=1; i <= n; i++ )
		{
			i_ = i - 1;
			iwork[i_] = i;
		}
		*nfe = 0;
		usubc.nfxe = 1;
		if( usubc.irepl == 0 )
		{
			isubc.fbonus = 0.0f;
		}
		else if( usubc.minf )
		{
			isubc.fbonus = bnsfac[0][usubc.ifxsw-1]*usubc.bonus;
		}
		else
		{
			isubc.fbonus = bnsfac[1][usubc.ifxsw-1]*usubc.bonus;
		}
		if( usubc.nfstop == 0 )
		{
			isubc.sfstop = 0.0f;
		}
		else if( usubc.minf )
		{
			isubc.sfstop = usubc.fstop;
		}
		else
		{
			isubc.sfstop = -usubc.fstop;
		}
		usubc.ftest = 0.0f;
		cmode = false;
		isubc.IntNew = true;
		usubc.initx = true;
		nnn = 0;
		evalf(nnn,iwork,(realnum*)&dum,n,x,&sfx,nfe);
		usubc.initx = false;

		/*       continuation of previous call
		 * */
	}
	else if( *iflag == 2 )
	{
		if( usubc.minf )
		{
			isubc.sfstop = usubc.fstop;
		}
		else
		{
			isubc.sfstop = -usubc.fstop;
		}
		cmode = true;
		goto L_70;
	}
	else if( *iflag == -1 )
	{
		cmode = true;
		goto L_70;
	}
	else if( *iflag == 0 )
	{
		cmode = false;
		goto L_90;
	}
	else
	{
		return;
	}

	/*     subplex loop
	 * */
L_40:

	for( i=0; i < n; i++ )
	{
		work[i] = (realnum)fabs(work[i]);
	}

	sortd(n,work,iwork);
	partx(n,iwork,work,&nsubs,&iwork[insptr-1]);
	cdcopy(n,x,1,work,1);
	ins = insptr;
	insfnl = insptr + nsubs - 1;
	ipptr = 1;

	/*       simplex loop
	 * */

L_60:
	ns = iwork[ins-1];

L_70:
	simplx(n,&work[istptr-1],ns,&iwork[ipptr-1],maxnfe,cmode,x,&sfx,
	  nfe,&work[isptr-1],&work[ifsptr-1],iflag);

	cmode = false;
	if( *iflag != 0 )
		goto L_121;

	if( ins < insfnl )
	{
		ins += 1;
		ipptr += ns;
		goto L_60;
	}

	/*       end simplex loop
	 * */
	for( i=0; i < n; i++ )
	{
		work[i] = x[i] - work[i];
	}

	/*       check termination
	 * */
L_90:
	/* new stop criterion: calculate distance between
	 * previous and new best model, if distance is
	 * smaller than tol return. */
	istep = istptr-1;
	xstop = 0.f;
	xstop2 = 0.f;
	for( i=0; i < n; i++ )
	{
		realnum ttemp;
		xstop += POW2(work[i]);
		ttemp = (realnum)fabs(work[istep]);
		/* chng to avoid side effects 
		 * xstop2 = MAX2(xstop2,(realnum)fabs(work[istep]));*/
		xstop2 = MAX2( xstop2 , ttemp );
		istep++;
	}

	if( sqrt(xstop) > tol || xstop2*usubc.psi > (realnum)tol )
	{
		setstp(nsubs,n,work,&work[istptr-1]);
		goto L_40;
	}

	/*     end subplex loop
	 * */

	*iflag = 0;
L_121:
	if( usubc.minf )
	{
		*fx = sfx;
	}
	else
	{
		*fx = -sfx;
	}
	return;

	/*     invalid input
	 * */
L_120:

	*iflag = -2;
	return;
}
/********************************************************************* */
STATIC void subopt(long int n)
{

	DEBUG_ENTRY( "subopt()" );


	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * subopt sets options for optimize_subplex.
	 *
	 * input
	 *
	 *   n      - problem dimension
	 *
	 * common
	 * */



	/* subroutines and functions
	 *
	 *   fortran
	 *
	 *-----------------------------------------------------------
	 *
	 ************************************************************
	 * simplex method strategy parameters
	 ************************************************************
	 *
	 * alpha  - reflection coefficient
	 *          alpha .gt. 0
	 * */
	usubc.alpha = 1.0f;

	/* beta   - contraction coefficient
	 *          0 .lt. beta .lt. 1
	 * */
	usubc.beta = .50f;

	/* gamma  - expansion coefficient
	 *          gamma .gt. 1
	 * */
	usubc.gamma = 2.0f;

	/* delta  - shrinkage (massive contraction) coefficient
	 *          0 .lt. delta .lt. 1
	 * */
	usubc.delta = .50f;

	/************************************************************
	 * subplex method strategy parameters
	 ************************************************************
	 *
	 * psi    - simplex reduction coefficient
	 *          0 .lt. psi .lt. 1
	 * */
	usubc.psi = .250f;

	/* omega  - step reduction coefficient
	 *          0 .lt. omega .lt. 1
	 * */
	usubc.omega = .10f;

	/* nsmin and nsmax specify a range of subspace dimensions.
	 * In addition to satisfying  1 .le. nsmin .le. nsmax .le. n,
	 * nsmin and nsmax must be chosen so that n can be expressed
	 * as a sum of positive integers where each of these integers
	 * ns(i) satisfies   nsmin .le. ns(i) .ge. nsmax.
	 * Specifically,
	 *     nsmin*ceil(n/nsmax) .le. n   must be true.
	 *
	 * nsmin  - subspace dimension minimum
	 * */
	usubc.nsmin = MIN2(2,n);

	/* nsmax  - subspace dimension maximum
	 * */
	usubc.nsmax = MIN2(5,n);

	/************************************************************
	 * subplex method special cases
	 ************************************************************
	 * nelder-mead simplex method with periodic restarts
	 *   nsmin = nsmax = n
	 ************************************************************
	 * nelder-mead simplex method
	 *   nsmin = nsmax = n, psi = small positive
	 ************************************************************
	 *
	 * irepl, ifxsw, and bonus deal with measurement replication.
	 * Objective functions subject to large amounts of noise can
	 * cause an optimization method to halt at a false optimum.
	 * An expensive solution to this problem is to evaluate f
	 * several times at each point and return the average (or max
	 * or min) of these trials as the function value.  optimize_subplex
	 * performs measurement replication only at the current best
	 * point. The longer a point is retained as best, the more
	 * accurate its function value becomes.
	 *
	 * The common variable nfxe contains the number of function
	 * evaluations at the current best point. fxstat contains the
	 * mean, max, min, and standard deviation of these trials.
	 *
	 * irepl  - measurement replication switch
	 * irepl  = 0, 1, or 2
	 *        = 0 : no measurement replication
	 *        = 1 : optimize_subplex performs measurement replication
	 *        = 2 : user performs measurement replication
	 *              (This is useful when optimizing on the mean,
	 *              max, or min of trials is insufficient. Common
	 *              variable initx is true for first function
	 *              evaluation. newx is true for first trial at
	 *              this point. The user uses subroutine fstats
	 *              within his objective function to maintain
	 *              fxstat. By monitoring newx, the user can tell
	 *              whether to return the function evaluation
	 *              (newx = .true.) or to use the new function
	 *              evaluation to refine the function evaluation
	 *              of the current best point (newx = .false.).
	 *              The common variable ftest gives the function
	 *              value that a new point must beat to be
	 *              considered the new best point.)
	 * */
	usubc.irepl = 0;

	/* ifxsw  - measurement replication optimization switch
	 * ifxsw  = 1, 2, or 3
	 *        = 1 : retain mean of trials as best function value
	 *        = 2 : retain max
	 *        = 3 : retain min
	 * */
	usubc.ifxsw = 1;

	/* Since the current best point will also be the most
	 * accurately evaluated point whenever irepl .gt. 0, a bonus
	 * should be added to the function value of the best point
	 * so that the best point is not replaced by a new point
	 * that only appears better because of noise.
	 * optimize_subplex uses bonus to determine how many multiples of
	 * fxstat(4) should be added as a bonus to the function
	 * evaluation. (The bonus is adjusted automatically by
	 * optimize_subplex when ifxsw or minf is changed.)
	 *
	 * bonus  - measurement replication bonus coefficient
	 *          bonus .ge. 0 (normally, bonus = 0 or 1)
	 *        = 0 : bonus not used
	 *        = 1 : bonus used
	 * */
	usubc.bonus = 1.0f;

	/* nfstop = 0 : f(x) is not tested against fstop
	 *        = 1 : if f(x) has reached fstop, optimize_subplex returns
	 *              iflag = 2
	 *        = 2 : (only valid when irepl .gt. 0)
	 *              if f(x) has reached fstop and
	 *              nfxe .gt. nfstop, optimize_subplex returns iflag = 2
	 * */
	usubc.nfstop = 0;

	/* fstop  - f target value
	 *          Its usage is determined by the value of nfstop.
	 *
	 * minf   - logical switch
	 *        = .true.  : optimize_subplex performs minimization
	 *        = .false. : optimize_subplex performs maximization
	 * */
	usubc.minf = true;
	return;
}
/**********************************************************************/
/* >>chng 01 jan 03, cleaned up -1 and formatting in this routine */
STATIC void cdcopy(long int n, 
	  realnum dx[], 
	  long int incx, 
	  realnum dy[], 
	  long int incy)
{
	long int i, 
	  ix, 
	  iy, 
	  m;

	DEBUG_ENTRY( "cdcopy()" );

	/* 
	 * copies a vector, x, to a vector, y.
	 * uses unrolled loops for increments equal to one.
	 * Jack Dongarra, lapack, 3/11/78.
	 */

	if( n > 0 )
	{
		if( incx == 1 && incy == 1 )
		{

			/* code for both increments equal to 1 */

			/* first the clean-up loop */
			m = n%7;
			if( m != 0 )
			{
				for( i=0; i < m; i++ )
				{
					dy[i] = dx[i];
				}
				if( n < 7 )
				{ 
					return;
				}
			}

			for( i=m; i < n; i += 7 )
			{
				dy[i] = dx[i];
				dy[i+1] = dx[i+1];
				dy[i+2] = dx[i+2];
				dy[i+3] = dx[i+3];
				dy[i+4] = dx[i+4];
				dy[i+5] = dx[i+5];
				dy[i+6] = dx[i+6];
			}
		}
		else
		{

			/* code for unequal increments or equal increments
			 * not equal to 1 */

			ix = 1;
			iy = 1;
			if( incx < 0 )
				ix = (-n + 1)*incx + 1;
			if( incy < 0 )
				iy = (-n + 1)*incy + 1;
			for( i=0; i < n; i++ )
			{
				dy[iy-1] = dx[ix-1];
				ix += incx;
				iy += incy;
			}
		}
	}
	return;
}
/********************************************************************* */
STATIC void evalf(long int ns, 
	  long int ips[], 
	  realnum xs[], 
	  long int n, 
	  realnum x[], 
	  realnum *sfx, 
	  long int *nfe)
{
	static int newbst;
	static long int i, 
	  i_;
	static realnum fx;
	/* gary delete since in header */
	/*double optimize_func();*/

	DEBUG_ENTRY( "evalf()" );
	/* gary add, not used, so trick compiler notes */
	i_ = n;

	/*>>chng 97 dec 07, implicit nonte, rid of first function argument
	 *
	 *                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * evalf evaluates the function f at a point defined by x
	 * with ns of its components replaced by those in xs.
	 *
	 * input
	 *
	 *   f      - user supplied function f(n,x) to be optimized
	 *   changed to optimize_func - cannot specify arbitrary function now
	 *
	 *   ns     - subspace dimension
	 *
	 *   ips    - permutation vector
	 *
	 *   xs     - real ns-vector to be mapped to x
	 *
	 *   n      - problem dimension
	 *
	 *   x      - real n-vector
	 *
	 *   nfe    - number of function evaluations
	 *
	 * output
	 *
	 *   sfx    - signed value of f evaluated at x
	 *
	 *   nfe    - incremented number of function evaluations
	 *
	 * common
	 * */




	/* local variables
	 * */


	/* subroutines and functions
	 * */

	/*-----------------------------------------------------------
	 * */
	for( i=1; i <= ns; i++ )
	{
		i_ = i - 1;
		x[ips[i_]-1] = xs[i_];
	}
	usubc.newx = isubc.IntNew || usubc.irepl != 2;
	fx = (realnum)optimize_func(x);
	/*      fx = f(n,x) */
	if( usubc.irepl == 0 )
	{
		if( usubc.minf )
		{
			*sfx = fx;
		}
		else
		{
			*sfx = -fx;
		}
	}
	else if( isubc.IntNew )
	{
		if( usubc.minf )
		{
			*sfx = fx;
			newbst = fx < usubc.ftest;
		}
		else
		{
			*sfx = -fx;
			newbst = fx > usubc.ftest;
		}
		if( usubc.initx || newbst )
		{
			if( usubc.irepl == 1 )
				fstats(fx,1,true);
			usubc.ftest = fx;
			isubc.sfbest = *sfx;
		}
	}
	else
	{
		if( usubc.irepl == 1 )
		{
			fstats(fx,1,false);
			fx = usubc.fxstat[usubc.ifxsw-1];
		}
		usubc.ftest = fx + isubc.fbonus*usubc.fxstat[3];
		if( usubc.minf )
		{
			*sfx = usubc.ftest;
			isubc.sfbest = fx;
		}
		else
		{
			*sfx = -usubc.ftest;
			isubc.sfbest = -fx;
		}
	}
	*nfe += 1;
	return;
}

/********************************************************************* */
STATIC void setstp(long int nsubs, 
	  long int n, 
	  realnum deltax[], 
	  realnum step[])
{
	static long int i, 
	  i_;
	static realnum stpfac;

	DEBUG_ENTRY( "setstp()" );


	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * setstp sets the stepsizes for the corresponding components
	 * of the solution vector.
	 *
	 * input
	 *
	 *   nsubs  - number of subspaces
	 *
	 *   n      - number of components (problem dimension)
	 *
	 *   deltax - vector of change in solution vector
	 *
	 *   step   - stepsizes for corresponding components of
	 *            solution vector
	 *
	 * output
	 *
	 *   step   - new stepsizes
	 *
	 * common
	 * */


	/* local variables
	 * */


	/* subroutines and functions
	 *
	 *   blas */
	/*   fortran
	 *
	 *-----------------------------------------------------------
	 *
	 *     set new step
	 * */
	if( nsubs > 1 )
	{	
		double a , b , c;
		a = cdasum(n,deltax,1);
		b = cdasum(n,step,1);
		c = MAX2(a/b ,usubc.omega);
		stpfac = (realnum)MIN2(c , 1.0f/usubc.omega);
	}
	else
	{
		stpfac = usubc.psi;
	}
	csscal(n,stpfac,step,1);

	/*     reorient simplex
	 * */
	for( i=1; i <= n; i++ )
	{
		i_ = i - 1;
		if( deltax[i_] == 0.f )
		{
			step[i_] = -step[i_];
		}
		else
		{
			step[i_] = (realnum)sign(step[i_],deltax[i_]);
		}
	}
	return;
}

/********************************************************************* */
STATIC void sortd(long int n, 
	  realnum xkey[], 
	  long int ix[])
{
	long int i, 
	  i_, 
	  ifirst, 
	  ilast, 
	  iswap, 
	  ixi, 
	  ixip1;

	DEBUG_ENTRY( "sortd()" );


	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * sortd uses the shakersort method to sort an array of keys
	 * in decreasing order. The sort is performed implicitly by
	 * modifying a vector of indices.
	 *
	 * For nearly sorted arrays, sortd requires O(n) comparisons.
	 * for completely unsorted arrays, sortd requires O(n**2)
	 * comparisons and will be inefficient unless n is small.
	 *
	 * input
	 *
	 *   n      - number of components
	 *
	 *   xkey   - real vector of keys
	 *
	 *   ix     - integer vector of indices
	 *
	 * output
	 *
	 *   ix     - indices satisfy xkey(ix(i)) .ge. xkey(ix(i+1))
	 *            for i = 1,...,n-1
	 *
	 * local variables
	 * */

	/*-----------------------------------------------------------
	 * */
	ifirst = 1;
	iswap = 1;
	ilast = n - 1;
	while( ifirst <= ilast )
	{
		for( i=ifirst; i <= ilast; i++ )
		{
			i_ = i - 1;
			ixi = ix[i_];
			ixip1 = ix[i_+1];
			if( xkey[ixi-1] < xkey[ixip1-1] )
			{
				ix[i_] = ixip1;
				ix[i_+1] = ixi;
				iswap = i;
			}
		}
		ilast = iswap - 1;
		for( i=ilast; i >= ifirst; i-- )
		{
			i_ = i - 1;
			ixi = ix[i_];
			ixip1 = ix[i_+1];
			if( xkey[ixi-1] < xkey[ixip1-1] )
			{
				ix[i_] = ixip1;
				ix[i_+1] = ixi;
				iswap = i;
			}
		}
		ifirst = iswap + 1;
	}
	return;
}

/********************************************************************* */
STATIC void partx(long int n, 
	  long int ip[], 
	  realnum absdx[], 
	  long int *nsubs, 
	  long int nsvals[])
{
	static long int i, 
	  limit, 
	  nleft, 
	  ns1, 
	  ns1_, 
	  ns2, 
	  nused;
	static realnum as1, 
	  as1max, 
	  as2, 
	  asleft, 
	  gap, 
	  gapmax;

	DEBUG_ENTRY( "partx()" );


	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * partx partitions the vector x by grouping components of
	 * similar magnitude of change.
	 *
	 * input
	 *
	 *   n      - number of components (problem dimension)
	 *
	 *   ip     - permutation vector
	 *
	 *   absdx  - vector of magnitude of change in x
	 *
	 *   nsvals - integer array dimensioned .ge. INT(n/nsmin)
	 *
	 * output
	 *
	 *   nsubs  - number of subspaces
	 *
	 *   nsvals - integer array of subspace dimensions
	 *
	 * common
	 * */


	/* local variables
	 * */


	/* subroutines and functions
	 *
	 *
	 *-----------------------------------------------------------
	 * */
	*nsubs = 0;
	nused = 0;
	nleft = n;
	asleft = absdx[0];
	for( i=1; i < n; i++ )
	{
		asleft += absdx[i];
	}

	while( nused < n )
	{
		*nsubs += 1;
		as1 = 0.0f;
		for( i=0; i < (usubc.nsmin - 1); i++ )
		{
			as1 += absdx[ip[nused+i]-1];
		}

		gapmax = -1.0f;
		limit = MIN2(usubc.nsmax,nleft);
		for( ns1=usubc.nsmin; ns1 <= limit; ns1++ )
		{
			ns1_ = ns1 - 1;
			as1 += absdx[ip[nused+ns1_]-1];
			ns2 = nleft - ns1;
			if( ns2 > 0 )
			{
				if( ns2 >= ((ns2 - 1)/usubc.nsmax + 1)*usubc.nsmin )
				{
					as2 = asleft - as1;
					gap = as1/ns1 - as2/ns2;
					if( gap > gapmax )
					{
						gapmax = gap;
						nsvals[*nsubs-1] = ns1;
						as1max = as1;
					}
				}
			}
			else if( as1/ns1 > gapmax )
			{
				goto L_21;
			}
		}
		nused += nsvals[*nsubs-1];
		nleft = n - nused;
		asleft -= as1max;
	}
	return;
L_21:
	nsvals[*nsubs-1] = ns1;
	return;
}

/********************************************************************* */

#undef S
#define S(I_,J_)	(*(s+(I_)*(ns)+(J_)))

STATIC void simplx(long int n, 
	  realnum step[], 
	  long int ns, 
	  long int ips[], 
	  long int maxnfe, 
	  int cmode, 
	  realnum x[], 
	  realnum *fx, 
	  long int *nfe, 
	  realnum *s, 
	  realnum fs[], 
	  long int *iflag)
{
	static int small, 
	  updatc;

	static long int i, 
	  icent, 
	  ih, 
	  il, 
	  inew = 0, 
	  is, 
	  itemp, 
	  j, 
	  j_, 
	  npts;

	static realnum dum, 
	  fc, 
	  fe, 
	  fr, 
	  tol;

	DEBUG_ENTRY( "simplx()" );


	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * simplx uses the Nelder-Mead simplex method to minimize the
	 * function f on a subspace.
	 *
	 * input
	 *
	 *   ffff   - function to be minimized, declared external in
	 *            calling routine
	 *   removed - always calls optimize_func
	 *
	 *   n      - problem dimension
	 *
	 *   step   - stepsizes for corresponding components of x
	 *
	 *   ns     - subspace dimension
	 *
	 *   ips    - permutation vector
	 *
	 *   maxnfe - maximum number of function evaluations
	 *
	 *   cmode  - logical switch
	 *            = .true.  : continuation of previous call
	 *            = .false. : first call
	 *
	 *   x      - starting guess for minimum
	 *
	 *   fx     - value of f at x
	 *
	 *   nfe    - number of function evaluations
	 *
	 *   s      - real work array of dimension .ge.
	 *            ns*(ns+3) used to store simplex
	 *
	 *   fs     - real work array of dimension .ge.
	 *            ns+1 used to store function values of simplex
	 *            vertices
	 *
	 * output
	 *
	 *   x      - computed minimum
	 *
	 *   fx     - value of f at x
	 *
	 *   nfe    - incremented number of function evaluations
	 *
	 *   iflag  - error flag
	 *            = -1 : maxnfe exceeded
	 *            =  0 : simplex reduced by factor of psi
	 *            =  1 : limit of machine precision
	 *            =  2 : reached fstop
	 *
	 * common
	 * */




	/* local variables
	 * */


	/* subroutines and functions
	 *
	 *     external optimize_func,calcc,dist,evalf,newpt,order,start */
	/*   blas */

	/*-----------------------------------------------------------
	 * */
	if( cmode )
		goto L_50;
	npts = ns + 1;
	icent = ns + 2;
	itemp = ns + 3;
	updatc = false;
	start(n,x,step,ns,ips,s,&small);
	if( small )
	{
		*iflag = 1;
		return;
	}
	else
	{
		if( usubc.irepl > 0 )
		{
			isubc.IntNew = false;
			evalf(ns,ips,&S(0,0),n,x,&fs[0],nfe);
		}
		else
		{
			fs[0] = *fx;
		}
		isubc.IntNew = true;
		for( j=2; j <= npts; j++ )
		{
			j_ = j - 1;
			evalf(ns,ips,&S(j_,0),n,x,&fs[j_],nfe);
		}
		il = 1;
		order(npts,fs,&il,&is,&ih);
		tol = (realnum)(usubc.psi*dist(ns,&S(ih-1,0),&S(il-1,0)));
	}

	/*     main loop
	 * */
L_20:
	calcc(ns,s,ih,inew,updatc,&S(icent-1,0));
	updatc = true;
	inew = ih;

	/*       reflect
	 * */
	newpt(ns,usubc.alpha,&S(icent-1,0),&S(ih-1,0),true,&S(itemp-1,0),
	  &small);
	if( !small )
	{
		evalf(ns,ips,&S(itemp-1,0),n,x,&fr,nfe);
		if( fr < fs[il-1] )
		{

			/*         expand
			 * */
			newpt(ns,-usubc.gamma,&S(icent-1,0),&S(itemp-1,0),true,
			  &S(ih-1,0),&small);
			if( small )
				goto L_40;
			evalf(ns,ips,&S(ih-1,0),n,x,&fe,nfe);
			if( fe < fr )
			{
				fs[ih-1] = fe;
			}
			else
			{
				cdcopy(ns,&S(itemp-1,0),1,&S(ih-1,0),1);
				fs[ih-1] = fr;
			}
		}
		else if( fr < fs[is-1] )
		{

			/*         accept reflected point
			 * */
			cdcopy(ns,&S(itemp-1,0),1,&S(ih-1,0),1);
			fs[ih-1] = fr;
		}
		else
		{

			/*         contract
			 * */
			if( fr > fs[ih-1] )
			{
				newpt(ns,-usubc.beta,&S(icent-1,0),&S(ih-1,0),true,
				  &S(itemp-1,0),&small);
			}
			else
			{
				newpt(ns,-usubc.beta,&S(icent-1,0),&S(itemp-1,0),false,
				  (realnum*)&dum,&small);
			}
			if( small )
				goto L_40;
			evalf(ns,ips,&S(itemp-1,0),n,x,&fc,nfe);
			if( fc < (realnum)MIN2(fr,fs[ih-1]) )
			{
				cdcopy(ns,&S(itemp-1,0),1,&S(ih-1,0),1);
				fs[ih-1] = fc;
			}
			else
			{

				/*           shrink simplex
				 * */
				for( j=1; j <= npts; j++ )
				{
					j_ = j - 1;
					if( j != il )
					{
						newpt(ns,-usubc.delta,&S(il-1,0),&S(j_,0),
						  false,(realnum*)&dum,&small);
						if( small )
							goto L_40;
						evalf(ns,ips,&S(j_,0),n,x,&fs[j_],nfe);
					}
				}
			}
			updatc = false;
		}
		order(npts,fs,&il,&is,&ih);
	}

	/*       check termination
	 * */

L_40:
	if( usubc.irepl == 0 )
	{
		*fx = fs[il-1];
	}
	else
	{
		*fx = isubc.sfbest;
	}

L_50:
	if( (usubc.nfstop > 0 && *fx <= isubc.sfstop) && usubc.nfxe >= 
	  usubc.nfstop )
		goto L_51;
	if( *nfe >= maxnfe )
		goto L_52;
	if( !(dist(ns,&S(ih-1,0),&S(il-1,0)) <= tol || small) )
		goto L_20;
	*iflag = 0;
	goto L_53;

L_52:
	*iflag = -1;
	goto L_53;

L_51:
	*iflag = 2;

	/*     end main loop, return best point
	 * */

L_53:
	for( i=0; i < ns; i++ )
	{
		x[ips[i]-1] = S(il-1,i);
	}
	return;
}

/********************************************************************* */
STATIC void fstats(double fx, 
	  long int ifxwt, 
	  int reset)
{
	static long int nsv;
	static realnum f1sv, 
	  fscale;

	DEBUG_ENTRY( "fstats()" );


	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * fstats modifies the common /usubc/ variables nfxe,fxstat.
	 *
	 * input
	 *
	 *   fx     - most recent evaluation of f at best x
	 *
	 *   ifxwt  - integer weight for fx
	 *
	 *   reset  - logical switch
	 *            = .true.  : initialize nfxe,fxstat
	 *            = .false. : update nfxe,fxstat
	 *
	 * common
	 * */


	/* local variables
	 * */


	/* subroutines and functions
	 *
	 *
	 *-----------------------------------------------------------
	 * */
	if( reset )
	{
		usubc.nfxe = ifxwt;
		usubc.fxstat[0] = (realnum)fx;
		usubc.fxstat[1] = (realnum)fx;
		usubc.fxstat[2] = (realnum)fx;
		usubc.fxstat[3] = 0.0f;
	}
	else
	{
		nsv = usubc.nfxe;
		f1sv = usubc.fxstat[0];
		usubc.nfxe += ifxwt;
		usubc.fxstat[0] += (realnum)(ifxwt*(fx - usubc.fxstat[0])/usubc.nfxe);
		usubc.fxstat[1] = MAX2(usubc.fxstat[1],(realnum)fx);
		usubc.fxstat[2] = MIN2(usubc.fxstat[2],(realnum)fx);
		fscale = (realnum)MAX3(fabs(usubc.fxstat[1]),fabs(usubc.fxstat[2]),1.);
		usubc.fxstat[3] = (realnum)(fscale*sqrt(((nsv-1)*POW2(usubc.fxstat[3]/
		  fscale)+nsv*POW2((usubc.fxstat[0]-f1sv)/fscale)+ifxwt*
		  POW2((fx-usubc.fxstat[0])/fscale))/(usubc.nfxe-1)));
	}
	return;
}

STATIC double cdasum(long int n, 
	  realnum dx[], 
	  long int incx)
{
	/*
	 *
	 *     takes the sum of the absolute values.
	 *     uses unrolled loops for increment equal to one.
	 *     jack dongarra, lapack, 3/11/78.
	 *     modified to correct problem with negative increment, 8/21/90.
	 *
	 */
	long int i, 
	  ix, 
	  m;
	realnum cdasum_v, 
	  dtemp;

	DEBUG_ENTRY( "cdasum()" );

	cdasum_v = 0.00f;
	dtemp = 0.00f;
	if( n > 0 )
	{
		if( incx == 1 )
		{

			/*        code for increment equal to 1
			 *
			 *
			 *        clean-up loop
			 * */
			m = n%6;
			if( m != 0 )
			{
				for( i=0; i < m; i++ )
				{
					dtemp += (realnum)fabs(dx[i]);
				}
				if( n < 6 )
					goto L_60;
			}

			for( i=m; i < n; i += 6 )
			{
				dtemp += (realnum)(fabs(dx[i]) + fabs(dx[i+1]) + fabs(dx[i+2]) + 
				  fabs(dx[i+3]) + fabs(dx[i+4]) + fabs(dx[i+5]));
			}
L_60:
			cdasum_v = dtemp;
		}
		else
		{

			/*        code for increment not equal to 1
			 * */
			ix = 1;

			if( incx < 0 )
				ix = (-n + 1)*incx + 1;

			for( i=0; i < n; i++ )
			{
				dtemp += (realnum)fabs(dx[ix-1]);
				ix += incx;
			}
			cdasum_v = dtemp;
		}
	}
	return( cdasum_v );
}

STATIC void csscal(long int n, 
	  double da, 
	  realnum dx[], 
	  long int incx)
{
	long int 
	  i, 
	  i_, 
	  m, 
	  mp1, 
	  nincx;

	DEBUG_ENTRY( "csscal()" );

	/*     scales a vector by a constant.
	 *     uses unrolled loops for increment equal to one.
	 *     jack dongarra, lapack, 3/11/78.
	 *     modified 3/93 to return if incx .le. 0.
	 *     modified 12/3/93, array(1) declarations changed to array(*)
	 *     changed to single precisions
	 * */

	if( !(n <= 0 || incx <= 0) )
	{
		if( incx == 1 )
		{

			/*        code for increment equal to 1
			 *
			 *
			 *        clean-up loop
			 * */
			m = n%5;
			if( m != 0 )
			{
				for( i=1; i <= m; i++ )
				{
					i_ = i - 1;
					dx[i_] *= (realnum)da;
				}
				if( n < 5 )
				{ 
					return;
				}
			}
			mp1 = m + 1;
			for( i=mp1; i <= n; i += 5 )
			{
				i_ = i - 1;
				dx[i_] *= (realnum)da;
				dx[i_+1] *= (realnum)da;
				dx[i_+2] *= (realnum)da;
				dx[i_+3] *= (realnum)da;
				dx[i_+4] *= (realnum)da;
			}
		}
		else
		{

			/*        code for increment not equal to 1
			 * */
			nincx = n*incx;
			/*for( i=1, _do0=DOCNT(i,nincx,_do1 = incx); _do0 > 0; i += _do1, _do0-- )*/
			/* gary change forc */
			for( i=0; i<nincx; i=i+incx)
			/*for( i=1, _do0=DOCNT(i,nincx,_do1 = incx); _do0 > 0; i += _do1, _do0-- )*/
			{
				dx[i] *= (realnum)da;
			}
		}
	}
	return;
}

/********************************************************************* */

#undef S
#define S(I_,J_)	(*(s+(I_)*(ns)+(J_)))

STATIC void start(long int, 
	  realnum x[], 
	  realnum step[], 
	  long int ns, 
	  long int ips[], 
	  realnum *s, 
	  int *small)
{
	DEBUG_ENTRY( "start()" );

	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * start creates the initial simplex for simplx minimization.
	 *
	 * input
	 *
	 *   n      - problem dimension (not used)
	 *
	 *   x      - current best point
	 *
	 *   step   - stepsizes for corresponding components of x
	 *
	 *   ns     - subspace dimension
	 *
	 *   ips    - permutation vector
	 *
	 *
	 * output
	 *
	 *   s      - first ns+1 columns contain initial simplex
	 *
	 *   small  - logical flag
	 *            = .true.  : coincident points
	 *            = .false. : otherwise
	 *
	 * local variables
	 * */

	/* subroutines and functions
	 *
	 *   blas */
	/*   fortran */

	/*-----------------------------------------------------------
	 * */
	for( long i=1; i <= ns; i++ )
	{
		long i_ = i - 1;
		S(0,i_) = x[ips[i_]-1];
	}

	for( long j=2; j <= (ns + 1); j++ )
	{
		long j_ = j - 1;
		cdcopy(ns,&S(0,0),1,&S(j_,0),1);
		S(j_,j_-1) = S(0,j_-1) + step[ips[j_-1]-1];
	}

	/* check for coincident points
	 * */
	for( long j=2; j <= (ns + 1); j++ )
	{
		long j_ = j - 1;
		if( (double)(S(j_,j_-1)) == (double)(S(0,j_-1)) )
			goto L_40;
	}
	*small = false;
	return;

	/* coincident points
	 * */
L_40:
	*small = true;
	return;
}

/********************************************************************* */
STATIC void order(long int npts, 
	  realnum fs[], 
	  long int *il, 
	  long int *is, 
	  long int *ih)
{
	long int i, 
	  il0, 
	  j;

	DEBUG_ENTRY( "order()" );


	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * order determines the indices of the vertices with the
	 * lowest, second highest, and highest function values.
	 *
	 * input
	 *
	 *   npts   - number of points in simplex
	 *
	 *   fs     - real vector of function values of
	 *            simplex
	 *
	 *   il     - index to vertex with lowest function value
	 *
	 * output
	 *
	 *   il     - new index to vertex with lowest function value
	 *
	 *   is     - new index to vertex with second highest
	 *            function value
	 *
	 *   ih     - new index to vertex with highest function value
	 *
	 * local variables
	 * */

	/* subroutines and functions
	 *
	 *
	 *-----------------------------------------------------------
	 * */
	il0 = *il;
	j = (il0%npts) + 1;

	if( fs[j-1] >= fs[*il-1] )
	{
		*ih = j;
		*is = il0;
	}
	else
	{
		*ih = il0;
		*is = j;
		*il = j;
	}

	for( i=il0 + 1; i <= (il0 + npts - 2); i++ )
	{
		j = (i%npts) + 1;
		if( fs[j-1] >= fs[*ih-1] )
		{
			*is = *ih;
			*ih = j;
		}
		else if( fs[j-1] > fs[*is-1] )
		{
			*is = j;
		}
		else if( fs[j-1] < fs[*il-1] )
		{
			*il = j;
		}
	}
	return;
}

STATIC double dist(long int n, 
	  realnum x[], 
	  realnum y[])
{
	/*
	 *
	 */
	long int i, 
	  i_;
	realnum absxmy, 
	  dist_v, 
	  scale, 
	  sum;

	DEBUG_ENTRY( "dist()" );

	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * dist calculates the distance between the points x,y.
	 *
	 * input
	 *
	 *   n      - number of components
	 *
	 *   x      - point in n-space
	 *
	 *   y      - point in n-space
	 *
	 * local variables
	 * */

	/* subroutines and functions
	 *
	 *   fortran
	 *
	 *-----------------------------------------------------------
	 * */
	absxmy = (realnum)fabs(x[0]-y[0]);
	if( absxmy <= 1.0f )
	{
		sum = absxmy*absxmy;
		scale = 1.0f;
	}
	else
	{
		sum = 1.0f;
		scale = absxmy;
	}

	for( i=2; i <= n; i++ )
	{
		i_ = i - 1;
		absxmy = (realnum)fabs(x[i_]-y[i_]);
		if( absxmy <= scale )
		{
			sum += POW2(absxmy/scale);
		}
		else
		{
			sum = 1.0f + sum*POW2(scale/absxmy);
			scale = absxmy;
		}
	}
	dist_v = (realnum)(scale*sqrt(sum));
	return( dist_v );
}

/********************************************************************* */

#undef S
#define S(I_,J_)	(*(s+(I_)*(ns)+(J_)))

STATIC void calcc(long int ns, 
	  realnum *s, 
	  long int ih, 
	  long int inew, 
	  int updatc, 
	  realnum c[])
{
	long int i, 
	  j, 
	  j_;
	/* >>chng 99 apr 21, following was not initialized, caught by Peter van Hoof */
	realnum xNothing[1] = { 0.0f };

	DEBUG_ENTRY( "calcc()" );


	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * calcc calculates the centroid of the simplex without the
	 * vertex with highest function value.
	 *
	 * input
	 *
	 *   ns     - subspace dimension
	 *
	 *   s      - real work space of dimension .ge.
	 *            ns*(ns+3) used to store simplex
	 *
	 *   ih     - index to vertex with highest function value
	 *
	 *   inew   - index to new point
	 *
	 *   updatc - logical switch
	 *            = .true.  : update centroid
	 *            = .false. : calculate centroid from scratch
	 *
	 *   c      - centroid of the simplex without vertex with
	 *            highest function value
	 *
	 * output
	 *
	 *   c      - new centroid
	 *
	 * local variables
	 * */
	/*     added to get prototypes to work */

	/* subroutines and functions
	 *
	 *   blas */

	/*-----------------------------------------------------------
	 * */
	if( !updatc )
	{
		/*       call cdcopy (ns,0.0,0,c,1)
		 *       xNothing will not be used since 0 is increment */
		cdcopy(ns,xNothing,0,c,1);
		for( j=1; j <= (ns + 1); j++ )
		{
			j_ = j - 1;
			if( j != ih )
				cdaxpy(ns,1.0f,&S(j_,0),1,c,1);
		}
		csscal(ns,1.0f/ns,c,1);
	}
	else if( ih != inew )
	{
		for( i=0; i < ns; i++ )
		{
			c[i] += (S(inew-1,i) - S(ih-1,i))/ns;
		}
	}
	return;
}

/********************************************************************* */
STATIC void newpt(long int ns, 
	  double coef, 
	  realnum xbase[], 
	  realnum xold[], 
	  int IntNew, 
	  realnum xnew[], 
	  int *small)
{
	int eqbase, 
	  eqold;
	long int i;
	realnum xoldi;

	DEBUG_ENTRY( "newpt()" );


	/*                                         Coded by Tom Rowan
	 *                            Department of Computer Sciences
	 *                              University of Texas at Austin
	 *
	 * newpt performs reflections, expansions, contractions, and
	 * shrinkages (massive contractions) by computing:
	 *
	 * xbase + coef * (xbase - xold)
	 *
	 * The result is stored in xnew if IntNew .eq. .true.,
	 * in xold otherwise.
	 *
	 * use :  coef .gt. 0 to reflect
	 *        coef .lt. 0 to expand, contract, or shrink
	 *
	 * input
	 *
	 *   ns     - number of components (subspace dimension)
	 *
	 *   coef   - one of four simplex method coefficients
	 *
	 *   xbase  - real ns-vector representing base
	 *            point
	 *
	 *   xold   - real ns-vector representing old
	 *            point
	 *
	 *   IntNew    - logical switch
	 *            = .true.  : store result in xnew
	 *            = .false. : store result in xold, xnew is not
	 *                        referenced
	 *
	 * output
	 *
	 *   xold   - unchanged if IntNew .eq. .true., contains new
	 *            point otherwise
	 *
	 *   xnew   - real ns-vector representing new
	 *            point if  new .eq. .true., not referenced
	 *            otherwise
	 *
	 *   small  - logical flag
	 *            = .true.  : coincident points
	 *            = .false. : otherwise
	 *
	 * local variables
	 * */

	/* subroutines and functions
	 *
	 *   fortran */

	/*-----------------------------------------------------------
	 * */
	eqbase = true;
	eqold = true;
	if( IntNew )
	{
		for( i=0; i < ns; i++ )
		{
			xnew[i] = (realnum)(xbase[i] + coef*(xbase[i] - xold[i]));
			eqbase = eqbase && ((double)(xnew[i]) == (double)(xbase[i]));
			eqold = eqold && ((double)(xnew[i]) == (double)(xold[i]));
		}
	}
	else
	{
		for( i=0; i < ns; i++ )
		{
			xoldi = xold[i];
			xold[i] = (realnum)(xbase[i] + coef*(xbase[i] - xold[i]));
			eqbase = eqbase && ((double)(xold[i]) == (double)(xbase[i]));
			eqold = eqold && ((double)(xold[i]) == (double)(xoldi));
		}
	}
	*small = eqbase || eqold;
	return;
}
/********************************************************************* */
STATIC void cdaxpy(long int n, 
	  double da, 
	  realnum dx[], 
	  long int incx, 
	  realnum dy[], 
	  long int incy)
{
	long int i, 
	  i_, 
	  ix, 
	  iy, 
	  m;

	DEBUG_ENTRY( "cdaxpy()" );

	/*     constant times a vector plus a vector.
	 *     uses unrolled loops for increments equal to one.
	 *     jack dongarra, lapack, 3/11/78.
	 * */

	if( n > 0 )
	{
		if( da != 0.00f )
		{
			if( incx == 1 && incy == 1 )
			{

				/*        code for both increments equal to 1
				 *
				 *
				 *        clean-up loop
				 * */
				m = n%4;
				if( m != 0 )
				{
					for( i=1; i <= m; i++ )
					{
						i_ = i - 1;
						dy[i_] += (realnum)(da*dx[i_]);
					}
					if( n < 4 )
					{
						return;
					}
				}

				for( i=m; i < n; i += 4 )
				{
					dy[i] += (realnum)(da*dx[i]);
					dy[i+1] += (realnum)(da*dx[i+1]);
					dy[i+2] += (realnum)(da*dx[i+2]);
					dy[i+3] += (realnum)(da*dx[i+3]);
				}
			}
			else
			{

				/*        code for unequal increments or equal increments
				 *          not equal to 1
				 * */
				ix = 1;
				iy = 1;
				if( incx < 0 )
					ix = (-n + 1)*incx + 1;
				if( incy < 0 )
					iy = (-n + 1)*incy + 1;
				for( i=0; i < n; i++ )
				{
					dy[iy-1] += (realnum)(da*dx[ix-1]);
					ix += incx;
					iy += incy;
				}
			}
		}
	}
	return;
}
/*lint +e801 goto is deprecated */
/*lint +e64 type mismatch */
