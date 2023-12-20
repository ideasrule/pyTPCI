#include "cddefines.h"
#include "thirdparty_quadpack.h"

/* translated by f2c (version 20100827). */

STATIC long xerror_(const char *mesg, const long *nmesg, 
				 const long *nerr, const long *level, long)
{
	DEBUG_ENTRY("xerror()");
	if (*level > 0)
	{
		fprintf(ioQQQ,"Error [%ld]: %*s\n",*nerr,int(*nmesg),mesg);
		cdEXIT(EXIT_FAILURE);
	}
	else
	{
		fprintf(ioQQQ,"Warning [%ld]: %*s\n",*nerr,int(*nmesg),mesg);
	}
	return 0;
}
static sys_float r1mach(long i)
{
	DEBUG_ENTRY("r1mach()");
	//  DOUBLE-PRECISION MACHINE CONSTANTS
   //  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
	//  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
	// D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
	// D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
	// D1MACH( 5) = LOG10(B)
	switch (i)
	{
	case 1:
		return FLT_MIN;
	case 2:
		return FLT_MAX;
	case 3:
		return FLT_EPSILON/FLT_RADIX;
	case 4:
		return FLT_EPSILON;
	case 5:
		return log10(double(FLT_RADIX));
	default:
		fprintf(stderr,"Error in input to r1mach");
		cdEXIT(EXIT_FAILURE);
	}
}
static double d1mach(long i)
{
	DEBUG_ENTRY("d1mach()");
	//  DOUBLE-PRECISION MACHINE CONSTANTS
   //  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
	//  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
	// D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
	// D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
	// D1MACH( 5) = LOG10(B)
	switch (i)
	{
	case 1:
		return DBL_MIN;
	case 2:
		return DBL_MAX;
	case 3:
		return DBL_EPSILON/FLT_RADIX;
	case 4:
		return DBL_EPSILON;
	case 5:
		return log10(double(FLT_RADIX));
	default:
		cdEXIT(EXIT_FAILURE);
	}
}

/* Table of constant values */

const static long c__4 = 4;
const static long c__1 = 1;
const static long c__26 = 26;
const static double c_b20 = 0.;
const static double c_b21 = 1.;
const static long c__2 = 2;
const static long c__0 = 0;
const static double c_b270 = 1.5;
const static double c_b320 = 2.;
const static sys_float c_b390 = 0.f;
const static sys_float c_b391 = 1.f;

static void sgtsl_(const long *, sys_float*, sys_float*, 
						 sys_float*, sys_float *, long *);
static void dgtsl_(const long *, double*, double*, double*, 
						 double *, long *);

void dqage_(const D_fp& f, const double *a, const double *b, const double *epsabs, 
				const double *epsrel, const long *key, const long *limit, 
				double *result, double *abserr, long *neval, long *ier, double *
				alist__, double *blist, double *rlist, double *elist, 
				long *iord, long *last)
{
	/* System generated locals */
	long i__1;
	double d__1, d__2;

	/* Local variables */
	long k;
	double a1, a2, b1, b2, area;
	long keyf;
	double area1, area2, area12, erro12, defab1, defab2;
	long nrmax;
	double uflow;
	long iroff1, iroff2;
	double error1, error2, defabs, epmach, errbnd, resabs, errmax;
	long maxerr;
	double errsum;

	/* ***begin prologue  dqage */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             integrand examinator, globally adaptive, */
	/*             gauss-kronrod */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral   i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            abs(i-reslt).le.max(epsabs,epsrel*abs(i)). */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            key    - long */
	/*                     key for choice of local integration rule */
	/*                     a gauss-kronrod pair is used with */
	/*                          7 - 15 points if key.lt.2, */
	/*                         10 - 21 points if key = 2, */
	/*                         15 - 31 points if key = 3, */
	/*                         20 - 41 points if key = 4, */
	/*                         25 - 51 points if key = 5, */
	/*                         30 - 61 points if key.gt.5. */

	/*            limit  - long */
	/*                     gives an upperbound on the number of subintervals */
	/*                     in the partition of (a,b), limit.ge.1. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed abs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for result and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value */
	/*                             of limit. */
	/*                             however, if this yields no improvement it */
	/*                             is rather advised to analyze the integrand */
	/*                             in order to determine the integration */
	/*                             difficulties. if the position of a local */
	/*                             difficulty can be determined(e.g. */
	/*                             singularity, discontinuity within the */
	/*                             interval) one will probably gain from */
	/*                             splitting up the interval at this point */
	/*                             and calling the integrator on the */
	/*                             subranges. if possible, an appropriate */
	/*                             special-purpose integrator should be used */
	/*                             which is designed for handling the type of */
	/*                             difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                             result, abserr, neval, last, rlist(1) , */
	/*                             elist(1) and iord(1) are set to zero. */
	/*                             alist(1) and blist(1) are set to a and b */
	/*                             respectively. */

	/*            alist   - double precision */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the left */
	/*                      end points of the subintervals in the partition */
	/*                      of the given integration range (a,b) */

	/*            blist   - double precision */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the right */
	/*                      end points of the subintervals in the partition */
	/*                      of the given integration range (a,b) */

	/*            rlist   - double precision */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the */
	/*                      integral approximations on the subintervals */

	/*            elist   - double precision */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the moduli of the */
	/*                      absolute error estimates on the subintervals */

	/*            iord    - int */
	/*                      vector of dimension at least limit, the first k */
	/*                      elements of which are pointers to the */
	/*                      error estimates over the subintervals, */
	/*                      such that elist(iord(1)), ..., */
	/*                      elist(iord(k)) form a decreasing sequence, */
	/*                      with k = last if last.le.(limit/2+2), and */
	/*                      k = limit+1-last otherwise */

	/*            last    - int */
	/*                      number of subintervals actually produced in the */
	/*                      subdivision process */

	/* ***references  (none) */
	/* ***routines called  d1mach,dqk15,dqk21,dqk31, */
	/*                    dqk41,dqk51,dqk61,dqpsrt */
	/* ***end prologue  dqage */




	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                      (alist(i),blist(i)) */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest */
	/*                       error estimate */
	/*           errmax    - elist(maxerr) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       abs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */


	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach  is the largest relative spacing. */
	/*           uflow  is the smallest positive magnitude. */

	/* ***first executable statement  dqage */
	/* Parameter adjustments */
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);

	/*           test on validity of parameters */
	/*           ------------------------------ */

	*ier = 0;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	alist__[1] = *a;
	blist[1] = *b;
	rlist[1] = 0.;
	elist[1] = 0.;
	iord[1] = 0;
	/* Computing MAX */
	d__1 = epmach * 50.;
	if (*epsabs <= 0. && *epsrel < max(d__1,5e-29)) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}

	/*           first approximation to the integral */
	/*           ----------------------------------- */

	keyf = *key;
	if (*key <= 0) {
		keyf = 1;
	}
	if (*key >= 7) {
		keyf = 6;
	}
	*neval = 0;
	if (keyf == 1) {
		dqk15_(f, a, b, result, abserr, &defabs, &resabs);
	}
	if (keyf == 2) {
		dqk21_(f, a, b, result, abserr, &defabs, &resabs);
	}
	if (keyf == 3) {
		dqk31_(f, a, b, result, abserr, &defabs, &resabs);
	}
	if (keyf == 4) {
		dqk41_(f, a, b, result, abserr, &defabs, &resabs);
	}
	if (keyf == 5) {
		dqk51_(f, a, b, result, abserr, &defabs, &resabs);
	}
	if (keyf == 6) {
		dqk61_(f, a, b, result, abserr, &defabs, &resabs);
	}
	*last = 1;
	rlist[1] = *result;
	elist[1] = *abserr;
	iord[1] = 1;

	/*           test on accuracy. */

	/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * fabs(*result);
	errbnd = max(d__1,d__2);
	if (*abserr <= epmach * 50. * defabs && *abserr > errbnd) {
		*ier = 2;
	}
	if (*limit == 1) {
		*ier = 1;
	}
	if (*ier != 0 || ( *abserr <= errbnd && *abserr != resabs ) || *abserr == 0.) 
	{
		goto L60;
	}

	/*           initialization */
	/*           -------------- */


	errmax = *abserr;
	maxerr = 1;
	area = *result;
	errsum = *abserr;
	nrmax = 1;
	iroff1 = 0;
	iroff2 = 0;

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limit;
	for (*last = 2; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with the largest error estimate. */

		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5;
		a2 = b1;
		b2 = blist[maxerr];
		if (keyf == 1) {
			dqk15_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 2) {
			dqk21_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 3) {
			dqk31_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 4) {
			dqk41_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 5) {
			dqk51_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 6) {
			dqk61_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 1) {
			dqk15_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}
		if (keyf == 2) {
			dqk21_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}
		if (keyf == 3) {
			dqk31_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}
		if (keyf == 4) {
			dqk41_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}
		if (keyf == 5) {
			dqk51_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}
		if (keyf == 6) {
			dqk61_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		++(*neval);
		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (defab1 == error1 || defab2 == error2) {
			goto L5;
		}
		if ((d__1 = rlist[maxerr] - area12, fabs(d__1)) <= fabs(area12) * 1e-5 
			 && erro12 >= errmax * .99) {
			++iroff1;
		}
		if (*last > 10 && erro12 > errmax) {
			++iroff2;
		}
	L5:
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(area);
		errbnd = max(d__1,d__2);
		if (errsum <= errbnd) {
			goto L8;
		}

		/*           test for roundoff error and eventually set error flag. */

		if (iroff1 >= 6 || iroff2 >= 20) {
			*ier = 2;
		}

		/*           set error flag in the case that the number of subintervals */
		/*           equals limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at a point of the integration range. */

		/* Computing MAX */
		d__1 = fabs(a1), d__2 = fabs(b2);
		if (max(d__1,d__2) <= (epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) 
		{
			*ier = 3;
		}

		/*           append the newly-created intervals to the list. */

	L8:
		if (error2 > error1) {
			goto L10;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L20;
	L10:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine dqpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the subinterval */
		/*           with the largest error estimate (to be bisected next). */

	L20:
		dqpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (*ier != 0 || errsum <= errbnd) {
			goto L40;
		}
		/* L30: */
	}

	/*           compute final result. */
	/*           --------------------- */

L40:
	*result = 0.;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L50: */
	}
	*abserr = errsum;
L60:
	if (keyf != 1) {
		*neval = (keyf * 10 + 1) * ((*neval << 1) + 1);
	}
	if (keyf == 1) {
		*neval = *neval * 30 + 15;
	}
L999:
	return;
} /* dqage_ */

void dqag_(const D_fp& f, const double *a, const double *b, const double *epsabs, 
			  const double *epsrel, const long *key, double *result, 
			  double *abserr, long *neval, long *ier, long *limit, 
			  const long *lenw, long *last, long *iwork, double *work)
{
	long l1, l2, l3, lvl;

	/* ***begin prologue  dqag */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             integrand examinator, globally adaptive, */
	/*             gauss-kronrod */
	/* ***author  piessens,robert,appl. math. & progr. div - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            abs(i-result)le.max(epsabs,epsrel*abs(i)). */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*            f      - double precision */
	/*                     function subprogam defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration */

	/*            epsabs - double precision */
	/*                     absolute accoracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            key    - int */
	/*                     key for choice of local integration rule */
	/*                     a gauss-kronrod pair is used with */
	/*                       7 - 15 points if key.lt.2, */
	/*                      10 - 21 points if key = 2, */
	/*                      15 - 31 points if key = 3, */
	/*                      20 - 41 points if key = 4, */
	/*                      25 - 51 points if key = 5, */
	/*                      30 - 61 points if key.gt.5. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed abs(i-result) */

	/*            neval  - int */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for result and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*                      error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yield no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulaties. */
	/*                             if the position of a local difficulty can */
	/*                             be determined (i.e.singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or limit.lt.1 or lenw.lt.limit*4. */
	/*                             result, abserr, neval, last are set */
	/*                             to zero. */
	/*                             except when lenw is invalid, iwork(1), */
	/*                             work(limit*2+1) and work(limit*3+1) are */
	/*                             set to zero, work(1) is set to a and */
	/*                             work(limit+1) to b. */

	/*         dimensioning parameters */
	/*            limit - int */
	/*                    dimensioning parameter for iwork */
	/*                    limit determines the maximum number of subintervals */
	/*                    in the partition of the given integration interval */
	/*                    (a,b), limit.ge.1. */
	/*                    if limit.lt.1, the routine will end with ier = 6. */

	/*            lenw  - long */
	/*                    dimensioning parameter for work */
	/*                    lenw must be at least limit*4. */
	/*                    if lenw.lt.limit*4, the routine will end with */
	/*                    ier = 6. */

	/*            last  - int */
	/*                    on return, last equals the number of subintervals */
	/*                    produced in the subdiviosion process, which */
	/*                    determines the number of significant elements */
	/*                    actually in the work arrays. */

	/*         work arrays */
	/*            iwork - int */
	/*                    vector of dimension at least limit, the first k */
	/*                    elements of which contain pointers to the error */
	/*                    estimates over the subintervals, such that */
	/*                    work(limit*3+iwork(1)),... , work(limit*3+iwork(k)) */
	/*                    form a decreasing sequence with k = last if */
	/*                    last.le.(limit/2+2), and k = limit+1-last otherwise */

	/*            work  - double precision */
	/*                    vector of dimension at least lenw */
	/*                    on return */
	/*                    work(1), ..., work(last) contain the left end */
	/*                    points of the subintervals in the partition of */
	/*                     (a,b), */
	/*                    work(limit+1), ..., work(limit+last) contain the */
	/*                     right end points, */
	/*                    work(limit*2+1), ..., work(limit*2+last) contain */
	/*                     the integral approximations over the subintervals, */
	/*                    work(limit*3+1), ..., work(limit*3+last) contain */
	/*                     the error estimates. */

	/* ***references  (none) */
	/* ***routines called  dqage,xerror */
	/* ***end prologue  dqag */



	/*         check validity of lenw. */

	/* ***first executable statement  dqag */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	if (*limit < 1 || *lenw < *limit << 2) {
		goto L10;
	}

	/*         prepare call for dqage. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;

	dqage_(f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, 
			 ier, &work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from dqag ", &c__26, ier, &lvl, 26);
	}
	return;
} /* dqag_ */

void dqagie_(const D_fp& f, const double *bound, const long *inf, 
				 const double *epsabs, const double *epsrel, const long *limit, 
				 double *result, double *abserr, long *neval, long *ier,
				 double *alist__, double *blist, double *rlist, double *elist, 
				 long *iord, long *last)
{
	/* System generated locals */
	long i__1, i__2;
	double d__1, d__2;

	/* Local variables */
	long k;
	double a1, a2, b1, b2;
	long id;
	double area, dres;
	long ksgn;
	double boun;
	long nres;
	double area1, area2, area12;
	double small=0., erro12;
	long ierro;
	double defab1, defab2;
	long ktmin, nrmax;
	double oflow, uflow;
	bool noext;
	long iroff1, iroff2, iroff3;
	double res3la[3], error1, error2, rlist2[52];
	long numrl2;
	double defabs, epmach, erlarg=0., abseps, correc=0., errbnd, resabs;
	long jupbnd;
	double erlast, errmax;
	long maxerr;
	double reseps;
	bool extrap;
	double ertest=0., errsum;

	/* ***begin prologue  dqagie */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a3a1,h2a4a1 */
	/* ***keywords  automatic integrator, infinite intervals, */
	/*             general-purpose, transformation, extrapolation, */
	/*             globally adaptive */
	/* ***author  piessens,robert,appl. math & progr. div - k.u.leuven */
	/*           de doncker,elise,appl. math & progr. div - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            integral   i = integral of f over (bound,+infinity) */
	/*            or i = integral of f over (-infinity,bound) */
	/*            or i = integral of f over (-infinity,+infinity), */
	/*            hopefully satisfying following claim for accuracy */
	/*            abs(i-result).le.max(epsabs,epsrel*abs(i)) */
	/* ***description */

	/* integration over infinite intervals */
	/* standard fortran subroutine */

	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            bound  - double precision */
	/*                     finite bound of integration range */
	/*                     (has no meaning if interval is doubly-infinite) */

	/*            inf    - double precision */
	/*                     indicating the kind of integration range involved */
	/*                     inf = 1 corresponds to  (bound,+infinity), */
	/*                     inf = -1            to  (-infinity,bound), */
	/*                     inf = 2             to (-infinity,+infinity). */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subintervals */
	/*                     in the partition of (a,b), limit.ge.1 */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed abs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                   - ier.gt.0 abnormal termination of the routine. the */
	/*                             estimates for result and error are less */
	/*                             reliable. it is assumed that the requested */
	/*                             accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however,if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. */
	/*                             if the position of a local difficulty can */
	/*                             be determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. */
	/*                             it is assumed that the requested tolerance */
	/*                             cannot be achieved, and that the returned */
	/*                             result is the best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                             result, abserr, neval, last, rlist(1), */
	/*                             elist(1) and iord(1) are set to zero. */
	/*                             alist(1) and blist(1) are set to 0 */
	/*                             and 1 respectively. */

	/*            alist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the left */
	/*                     end points of the subintervals in the partition */
	/*                     of the transformed integration range (0,1). */

	/*            blist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the right */
	/*                     end points of the subintervals in the partition */
	/*                     of the transformed integration range (0,1). */

	/*            rlist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the integral */
	/*                     approximations on the subintervals */

	/*            elist  - double precision */
	/*                     vector of dimension at least limit,  the first */
	/*                     last elements of which are the moduli of the */
	/*                     absolute error estimates on the subintervals */

	/*            iord   - long */
	/*                     vector of dimension limit, the first k */
	/*                     elements of which are pointers to the */
	/*                     error estimates over the subintervals, */
	/*                     such that elist(iord(1)), ..., elist(iord(k)) */
	/*                     form a decreasing sequence, with k = last */
	/*                     if last.le.(limit/2+2), and k = limit+1-last */
	/*                     otherwise */

	/*            last   - long */
	/*                     number of subintervals actually produced */
	/*                     in the subdivision process */

	/* ***references  (none) */
	/* ***routines called  d1mach,dqelg,dqk15i,dqpsrt */
	/* ***end prologue  dqagie */



	/*            the dimension of rlist2 is determined by the value of */
	/*            limexp in subroutine dqelg. */


	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           rlist2    - array of dimension at least (limexp+2), */
	/*                       containing the part of the epsilon table */
	/*                       wich is still needed for further computations */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest error */
	/*                       estimate */
	/*           errmax    - elist(maxerr) */
	/*           erlast    - error on the interval currently subdivided */
	/*                       (before that subdivision has taken place) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       abs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */
	/*           nres      - number of calls to the extrapolation routine */
	/*           numrl2    - number of elements currently in rlist2. if an */
	/*                       appropriate approximation to the compounded */
	/*                       integral has been obtained, it is put in */
	/*                       rlist2(numrl2) after numrl2 has been increased */
	/*                       by one. */
	/*           small     - length of the smallest interval considered up */
	/*                       to now, multiplied by 1.5 */
	/*           erlarg    - sum of the errors over the intervals larger */
	/*                       than the smallest interval considered up to now */
	/*           extrap    - bool variable denoting that the routine */
	/*                       is attempting to perform extrapolation. i.e. */
	/*                       before subdividing the smallest interval we */
	/*                       try to decrease the value of erlarg. */
	/*           noext     - bool variable denoting that extrapolation */
	/*                       is no longer allowed (true-value) */

	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */
	/*           oflow is the largest positive magnitude. */

	/* ***first executable statement  dqagie */
	/* Parameter adjustments */
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = d1mach(c__4);

	/*           test on validity of parameters */
	/*           ----------------------------- */

	*ier = 0;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	alist__[1] = 0.;
	blist[1] = 1.;
	rlist[1] = 0.;
	elist[1] = 0.;
	iord[1] = 0;
	/* Computing MAX */
	d__1 = epmach * 50.;
	if (*epsabs <= 0. && *epsrel < max(d__1,5e-29)) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}


	/*           first approximation to the integral */
	/*           ----------------------------------- */

	/*           determine the interval to be mapped onto (0,1). */
	/*           if inf = 2 the integral is computed as i = i1+i2, where */
	/*           i1 = integral of f over (-infinity,0), */
	/*           i2 = integral of f over (0,+infinity). */

	boun = *bound;
	if (*inf == 2) {
		boun = 0.;
	}
	dqk15i_(f, &boun, inf, &c_b20, &c_b21, result, abserr, &defabs, &
			  resabs);

	/*           test on accuracy */

	*last = 1;
	rlist[1] = *result;
	elist[1] = *abserr;
	iord[1] = 1;
	dres = fabs(*result);
	/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * dres;
	errbnd = max(d__1,d__2);
	if (*abserr <= epmach * 100. * defabs && *abserr > errbnd) {
		*ier = 2;
	}
	if (*limit == 1) {
		*ier = 1;
	}
	if (*ier != 0 || ( *abserr <= errbnd && *abserr != resabs ) || *abserr == 0.) 
	{
		goto L130;
	}

	/*           initialization */
	/*           -------------- */

	uflow = d1mach(c__1);
	oflow = d1mach(c__2);
	rlist2[0] = *result;
	errmax = *abserr;
	maxerr = 1;
	area = *result;
	errsum = *abserr;
	*abserr = oflow;
	nrmax = 1;
	nres = 0;
	ktmin = 0;
	numrl2 = 2;
	extrap = false;
	noext = false;
	ierro = 0;
	iroff1 = 0;
	iroff2 = 0;
	iroff3 = 0;
	ksgn = -1;
	if (dres >= (1. - epmach * 50.) * defabs) {
		ksgn = 1;
	}

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limit;
	for (*last = 2; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with nrmax-th largest error estimate. */

		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5;
		a2 = b1;
		b2 = blist[maxerr];
		erlast = errmax;
		dqk15i_(f, &boun, inf, &a1, &b1, &area1, &error1, &resabs, &
				  defab1);
		dqk15i_(f, &boun, inf, &a2, &b2, &area2, &error2, &resabs, &
				  defab2);

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (defab1 == error1 || defab2 == error2) {
			goto L15;
		}
		if ((d__1 = rlist[maxerr] - area12, fabs(d__1)) > fabs(area12) * 1e-5 ||
			 erro12 < errmax * .99) {
			goto L10;
		}
		if (extrap) {
			++iroff2;
		}
		if (! extrap) {
			++iroff1;
		}
	L10:
		if (*last > 10 && erro12 > errmax) {
			++iroff3;
		}
	L15:
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(area);
		errbnd = max(d__1,d__2);

		/*           test for roundoff error and eventually set error flag. */

		if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
			*ier = 2;
		}
		if (iroff2 >= 5) {
			ierro = 3;
		}

		/*           set error flag in the case that the number of */
		/*           subintervals equals limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at some points of the integration range. */

		/* Computing MAX */
		d__1 = fabs(a1), d__2 = fabs(b2);
		if (max(d__1,d__2) <= (epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) 
		{
			*ier = 4;
		}

		/*           append the newly-created intervals to the list. */

		if (error2 > error1) {
			goto L20;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L30;
	L20:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine dqpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the subinterval */
		/*           with nrmax-th largest error estimate (to be bisected next). */

	L30:
		dqpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		if (errsum <= errbnd) {
			goto L115;
		}
		if (*ier != 0) {
			goto L100;
		}
		if (*last == 2) {
			goto L80;
		}
		if (noext) {
			goto L90;
		}
		erlarg -= erlast;
		if ((d__1 = b1 - a1, fabs(d__1)) > small) {
			erlarg += erro12;
		}
		if (extrap) {
			goto L40;
		}

		/*           test whether the interval to be bisected next is the */
		/*           smallest interval. */

		if ((d__1 = blist[maxerr] - alist__[maxerr], fabs(d__1)) > small) {
			goto L90;
		}
		extrap = true;
		nrmax = 2;
	L40:
		if (ierro == 3 || erlarg <= ertest) {
			goto L60;
		}

		/*           the smallest interval has the largest error. */
		/*           before bisecting decrease the sum of the errors over the */
		/*           larger intervals (erlarg) and perform extrapolation. */

		id = nrmax;
		jupbnd = *last;
		if (*last > *limit / 2 + 2) {
			jupbnd = *limit + 3 - *last;
		}
		i__2 = jupbnd;
		for (k = id; k <= i__2; ++k) {
			maxerr = iord[nrmax];
			errmax = elist[maxerr];
			if ((d__1 = blist[maxerr] - alist__[maxerr], fabs(d__1)) > small) {
				goto L90;
			}
			++nrmax;
			/* L50: */
		}

		/*           perform extrapolation. */

	L60:
		++numrl2;
		rlist2[numrl2 - 1] = area;
		dqelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
		++ktmin;
		if (ktmin > 5 && *abserr < errsum * .001) {
			*ier = 5;
		}
		if (abseps >= *abserr) {
			goto L70;
		}
		ktmin = 0;
		*abserr = abseps;
		*result = reseps;
		correc = erlarg;
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(reseps);
		ertest = max(d__1,d__2);
		if (*abserr <= ertest) {
			goto L100;
		}

		/*            prepare bisection of the smallest interval. */

	L70:
		if (numrl2 == 1) {
			noext = true;
		}
		if (*ier == 5) {
			goto L100;
		}
		maxerr = iord[1];
		errmax = elist[maxerr];
		nrmax = 1;
		extrap = false;
		small *= .5;
		erlarg = errsum;
		goto L90;
	L80:
		small = .375;
		erlarg = errsum;
		ertest = errbnd;
		rlist2[1] = area;
	L90:
		;
	}

	/*           set final result and error estimate. */
	/*           ------------------------------------ */

L100:
	if (*abserr == oflow) {
		goto L115;
	}
	if (*ier + ierro == 0) {
		goto L110;
	}
	if (ierro == 3) {
		*abserr += correc;
	}
	if (*ier == 0) {
		*ier = 3;
	}
	if (*result != 0. && area != 0.) {
		goto L105;
	}
	if (*abserr > errsum) {
		goto L115;
	}
	if (area == 0.) {
		goto L130;
	}
	goto L110;
L105:
	if (*abserr / fabs(*result) > errsum / fabs(area)) {
		goto L115;
	}

	/*           test on divergence */

L110:
	/* Computing MAX */
	d__1 = fabs(*result), d__2 = fabs(area);
	if (ksgn == -1 && max(d__1,d__2) <= defabs * .01) {
		goto L130;
	}
	if (.01 > *result / area || *result / area > 100. || errsum > fabs(area)) {
		*ier = 6;
	}
	goto L130;

	/*           compute global integral sum. */

L115:
	*result = 0.;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L120: */
	}
	*abserr = errsum;
L130:
	*neval = *last * 30 - 15;
	if (*inf == 2) {
		*neval <<= 1;
	}
	if (*ier > 2) {
		--(*ier);
	}
L999:
	return;
} /* dqagie_ */

void dqagi_(const D_fp& f, const double *bound, const long *inf, 
				const double *epsabs, const double *epsrel, double *result, 
				double *abserr, long *neval, long *ier, long *limit, 
				const long *lenw, long *last, long *iwork, double *work)
{
	long l1, l2, l3, lvl;

	/* ***begin prologue  dqagi */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a3a1,h2a4a1 */
	/* ***keywords  automatic integrator, infinite intervals, */
	/*             general-purpose, transformation, extrapolation, */
	/*             globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. -k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            integral   i = integral of f over (bound,+infinity) */
	/*            or i = integral of f over (-infinity,bound) */
	/*            or i = integral of f over (-infinity,+infinity) */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        integration over infinite intervals */
	/*        standard fortran subroutine */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            bound  - double precision */
	/*                     finite bound of integration range */
	/*                     (has no meaning if interval is doubly-infinite) */

	/*            inf    - long */
	/*                     indicating the kind of integration range involved */
	/*                     inf = 1 corresponds to  (bound,+infinity), */
	/*                     inf = -1            to  (-infinity,bound), */
	/*                     inf = 2             to (-infinity,+infinity). */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */


	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                   - ier.gt.0 abnormal termination of the routine. the */
	/*                             estimates for result and error are less */
	/*                             reliable. it is assumed that the requested */
	/*                             accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. if */
	/*                             the position of a local difficulty can be */
	/*                             determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. */
	/*                             it is assumed that the requested tolerance */
	/*                             cannot be achieved, and that the returned */
	/*                             result is the best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                              or limit.lt.1 or leniw.lt.limit*4. */
	/*                             result, abserr, neval, last are set to */
	/*                             zero. exept when limit or leniw is */
	/*                             invalid, iwork(1), work(limit*2+1) and */
	/*                             work(limit*3+1) are set to zero, work(1) */
	/*                             is set to a and work(limit+1) to b. */

	/*         dimensioning parameters */
	/*            limit - long */
	/*                    dimensioning parameter for iwork */
	/*                    limit determines the maximum number of subintervals */
	/*                    in the partition of the given integration interval */
	/*                    (a,b), limit.ge.1. */
	/*                    if limit.lt.1, the routine will end with ier = 6. */

	/*            lenw  - long */
	/*                    dimensioning parameter for work */
	/*                    lenw must be at least limit*4. */
	/*                    if lenw.lt.limit*4, the routine will end */
	/*                    with ier = 6. */

	/*            last  - long */
	/*                    on return, last equals the number of subintervals */
	/*                    produced in the subdivision process, which */
	/*                    determines the number of significant elements */
	/*                    actually in the work arrays. */

	/*         work arrays */
	/*            iwork - long */
	/*                    vector of dimension at least limit, the first */
	/*                    k elements of which contain pointers */
	/*                    to the error estimates over the subintervals, */
	/*                    such that work(limit*3+iwork(1)),... , */
	/*                    work(limit*3+iwork(k)) form a decreasing */
	/*                    sequence, with k = last if last.le.(limit/2+2), and */
	/*                    k = limit+1-last otherwise */

	/*            work  - double precision */
	/*                    vector of dimension at least lenw */
	/*                    on return */
	/*                    work(1), ..., work(last) contain the left */
	/*                     end points of the subintervals in the */
	/*                     partition of (a,b), */
	/*                    work(limit+1), ..., work(limit+last) contain */
	/*                     the right end points, */
	/*                    work(limit*2+1), ...,work(limit*2+last) contain the */
	/*                     integral approximations over the subintervals, */
	/*                    work(limit*3+1), ..., work(limit*3) */
	/*                     contain the error estimates. */
	/* ***references  (none) */
	/* ***routines called  dqagie,xerror */
	/* ***end prologue  dqagi */




	/*         check validity of limit and lenw. */

	/* ***first executable statement  dqagi */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	if (*limit < 1 || *lenw < *limit << 2) {
		goto L10;
	}

	/*         prepare call for dqagie. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;

	dqagie_(f, bound, inf, epsabs, epsrel, limit, result, abserr, neval,
			  ier, &work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from dqagi", &c__26, ier, &lvl, 26);
	}
	return;
} /* dqagi_ */

void dqagpe_(const D_fp& f, const double *a, const double *b, const long *npts2, 
				 const double *points, const double *epsabs, const double *epsrel, 
				 const long *limit, double *result, double *abserr, long *
				 neval, long *ier, double *alist__, double *blist, 
				 double *rlist, double *elist, double *pts, long *iord, 
				 long *level, long *ndin, long *last)
{
	/* System generated locals */
	int i__1, i__2;
	double d__1, d__2;

	/* Local variables */
	long i__, j, k=0;
	double a1, a2, b1, b2;
	long id, ip1, ind1, ind2;
	double area;
	double resa, dres, sign;
	long ksgn;
	double temp;
	long nres, nint, jlow, npts;
	double area1, area2, area12;
	double erro12;
	long ierro;
	double defab1, defab2;
	long ktmin, nrmax;
	double oflow, uflow;
	bool noext;
	long iroff1, iroff2, iroff3;
	double res3la[3];
	long nintp1;
	double error1, error2, rlist2[52];
	long numrl2;
	double defabs, epmach, erlarg, abseps, correc=0., errbnd, resabs;
	long jupbnd;
	double erlast;
	long levmax;
	double errmax;
	long maxerr, levcur;
	double reseps;
	bool extrap;
	double ertest, errsum;

	/* ***begin prologue  dqagpe */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             singularities at user specified points, */
	/*             extrapolation, globally adaptive. */
	/* ***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f over (a,b), hopefully */
	/*            satisfying following claim for accuracy fabs(i-result).le. */
	/*            max(epsabs,epsrel*fabs(i)). break points of the integration */
	/*            interval, where local difficulties of the integrand may */
	/*            occur(e.g. singularities,discontinuities),provided by user. */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration */

	/*            npts2  - long */
	/*                     number equal to two more than the number of */
	/*                     user-supplied break points within the integration */
	/*                     range, npts2.ge.2. */
	/*                     if npts2.lt.2, the routine will end with ier = 6. */

	/*            points - double precision */
	/*                     vector of dimension npts2, the first (npts2-2) */
	/*                     elements of which are the user provided break */
	/*                     points. if these points do not constitute an */
	/*                     ascending sequence there will be an automatic */
	/*                     sorting. */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subintervals */
	/*                     in the partition of (a,b), limit.ge.npts2 */
	/*                     if limit.lt.npts2, the routine will end with */
	/*                     ier = 6. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine. */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. if */
	/*                             the position of a local difficulty can be */
	/*                             determined (i.e. singularity, */
	/*                             discontinuity within the interval), it */
	/*                             should be supplied to the routine as an */
	/*                             element of the vector points. if necessary */
	/*                             an appropriate special-purpose integrator */
	/*                             must be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. it is presumed that */
	/*                             the requested tolerance cannot be */
	/*                             achieved, and that the returned result is */
	/*                             the best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier.gt.0. */
	/*                         = 6 the input is invalid because */
	/*                             npts2.lt.2 or */
	/*                             break points are specified outside */
	/*                             the integration range or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or limit.lt.npts2. */
	/*                             result, abserr, neval, last, rlist(1), */
	/*                             and elist(1) are set to zero. alist(1) and */
	/*                             blist(1) are set to a and b respectively. */

	/*            alist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the left end points */
	/*                     of the subintervals in the partition of the given */
	/*                     integration range (a,b) */

	/*            blist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the right end points */
	/*                     of the subintervals in the partition of the given */
	/*                     integration range (a,b) */

	/*            rlist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the integral */
	/*                     approximations on the subintervals */

	/*            elist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the moduli of the */
	/*                     absolute error estimates on the subintervals */

	/*            pts    - double precision */
	/*                     vector of dimension at least npts2, containing the */
	/*                     integration limits and the break points of the */
	/*                     interval in ascending sequence. */

	/*            level  - long */
	/*                     vector of dimension at least limit, containing the */
	/*                     subdivision levels of the subinterval, i.e. if */
	/*                     (aa,bb) is a subinterval of (p1,p2) where p1 as */
	/*                     well as p2 is a user-provided break point or */
	/*                     integration limit, then (aa,bb) has level l if */
	/*                     fabs(bb-aa) = fabs(p2-p1)*2**(-l). */

	/*            ndin   - long */
	/*                     vector of dimension at least npts2, after first */
	/*                     integration over the intervals (pts(i)),pts(i+1), */
	/*                     i = 0,1, ..., npts2-2, the error estimates over */
	/*                     some of the intervals may have been increased */
	/*                     artificially, in order to put their subdivision */
	/*                     forward. if this happens for the subinterval */
	/*                     numbered k, ndin(k) is put to 1, otherwise */
	/*                     ndin(k) = 0. */

	/*            iord   - long */
	/*                     vector of dimension at least limit, the first k */
	/*                     elements of which are pointers to the */
	/*                     error estimates over the subintervals, */
	/*                     such that elist(iord(1)), ..., elist(iord(k)) */
	/*                     form a decreasing sequence, with k = last */
	/*                     if last.le.(limit/2+2), and k = limit+1-last */
	/*                     otherwise */

	/*            last   - long */
	/*                     number of subintervals actually produced in the */
	/*                     subdivisions process */

	/* ***references  (none) */
	/* ***routines called  d1mach,dqelg,dqk21,dqpsrt */
	/* ***end prologue  dqagpe */




	/*            the dimension of rlist2 is determined by the value of */
	/*            limexp in subroutine epsalg (rlist2 should be of dimension */
	/*            (limexp+2) at least). */


	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           rlist2    - array of dimension at least limexp+2 */
	/*                       containing the part of the epsilon table which */
	/*                       is still needed for further computations */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest error */
	/*                       estimate */
	/*           errmax    - elist(maxerr) */
	/*           erlast    - error on the interval currently subdivided */
	/*                       (before that subdivision has taken place) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */
	/*           nres      - number of calls to the extrapolation routine */
	/*           numrl2    - number of elements in rlist2. if an appropriate */
	/*                       approximation to the compounded integral has */
	/*                       been obtained, it is put in rlist2(numrl2) after */
	/*                       numrl2 has been increased by one. */
	/*           erlarg    - sum of the errors over the intervals larger */
	/*                       than the smallest interval considered up to now */
	/*           extrap    - bool variable denoting that the routine */
	/*                       is attempting to perform extrapolation. i.e. */
	/*                       before subdividing the smallest interval we */
	/*                       try to decrease the value of erlarg. */
	/*           noext     - bool variable denoting that extrapolation is */
	/*                       no longer allowed (true-value) */

	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */
	/*           oflow is the largest positive magnitude. */

	/* ***first executable statement  dqagpe */
	/* Parameter adjustments */
	--ndin;
	--pts;
	--points;
	--level;
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = d1mach(c__4);

	/*            test on validity of parameters */
	/*            ----------------------------- */

	*ier = 0;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	alist__[1] = *a;
	blist[1] = *b;
	rlist[1] = 0.;
	elist[1] = 0.;
	iord[1] = 0;
	level[1] = 0;
	npts = *npts2 - 2;
	/* Computing MAX */
	d__1 = epmach * 50.;
	if (*npts2 < 2 || *limit <= npts || 
		 ( *epsabs <= 0. && *epsrel < max(d__1,5e-29 ))) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}

	/*            if any break points are provided, sort them into an */
	/*            ascending sequence. */

	sign = 1.;
	if (*a > *b) {
		sign = -1.;
	}
	pts[1] = min(*a,*b);
	if (npts == 0) {
		goto L15;
	}
	i__1 = npts;
	for (i__ = 1; i__ <= i__1; ++i__) {
		pts[i__ + 1] = points[i__];
		/* L10: */
	}
L15:
	pts[npts + 2] = max(*a,*b);
	nint = npts + 1;
	a1 = pts[1];
	if (npts == 0) {
		goto L40;
	}
	nintp1 = nint + 1;
	i__1 = nint;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ip1 = i__ + 1;
		i__2 = nintp1;
		for (j = ip1; j <= i__2; ++j) {
			if (pts[i__] <= pts[j]) {
				goto L20;
			}
			temp = pts[i__];
			pts[i__] = pts[j];
			pts[j] = temp;
		L20:
			;
		}
	}
	if (pts[1] != min(*a,*b) || pts[nintp1] != max(*a,*b)) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}

	/*            compute first integral and error approximations. */
	/*            ------------------------------------------------ */

L40:
	resabs = 0.;
	i__2 = nint;
	for (i__ = 1; i__ <= i__2; ++i__) {
		b1 = pts[i__ + 1];
		dqk21_(f, &a1, &b1, &area1, &error1, &defabs, &resa);
		*abserr += error1;
		*result += area1;
		ndin[i__] = 0;
		if (error1 == resa && error1 != 0.) {
			ndin[i__] = 1;
		}
		resabs += defabs;
		level[i__] = 0;
		elist[i__] = error1;
		alist__[i__] = a1;
		blist[i__] = b1;
		rlist[i__] = area1;
		iord[i__] = i__;
		a1 = b1;
		/* L50: */
	}
	errsum = 0.;
	i__2 = nint;
	for (i__ = 1; i__ <= i__2; ++i__) {
		if (ndin[i__] == 1) {
			elist[i__] = *abserr;
		}
		errsum += elist[i__];
		/* L55: */
	}

	/*           test on accuracy. */

	*last = nint;
	*neval = nint * 21;
	dres = fabs(*result);
	/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * dres;
	errbnd = max(d__1,d__2);
	if (*abserr <= epmach * 100. * resabs && *abserr > errbnd) {
		*ier = 2;
	}
	if (nint == 1) {
		goto L80;
	}
	i__2 = npts;
	for (i__ = 1; i__ <= i__2; ++i__) {
		jlow = i__ + 1;
		ind1 = iord[i__];
		i__1 = nint;
		for (j = jlow; j <= i__1; ++j) {
			ind2 = iord[j];
			if (elist[ind1] > elist[ind2]) {
				goto L60;
			}
			ind1 = ind2;
			k = j;
		L60:
			;
		}
		if (ind1 == iord[i__]) {
			goto L70;
		}
		iord[k] = iord[i__];
		iord[i__] = ind1;
	L70:
		;
	}
	if (*limit < *npts2) {
		*ier = 1;
	}
L80:
	if (*ier != 0 || *abserr <= errbnd) {
		goto L210;
	}

	/*           initialization */
	/*           -------------- */

	rlist2[0] = *result;
	maxerr = iord[1];
	errmax = elist[maxerr];
	area = *result;
	nrmax = 1;
	nres = 0;
	numrl2 = 1;
	ktmin = 0;
	extrap = false;
	noext = false;
	erlarg = errsum;
	ertest = errbnd;
	levmax = 1;
	iroff1 = 0;
	iroff2 = 0;
	iroff3 = 0;
	ierro = 0;
	uflow = d1mach(c__1);
	oflow = d1mach(c__2);
	*abserr = oflow;
	ksgn = -1;
	if (dres >= (1. - epmach * 50.) * resabs) {
		ksgn = 1;
	}

	/*           main do-loop */
	/*           ------------ */

	i__2 = *limit;
	for (*last = *npts2; *last <= i__2; ++(*last)) {

		/*           bisect the subinterval with the nrmax-th largest error */
		/*           estimate. */

		levcur = level[maxerr] + 1;
		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5;
		a2 = b1;
		b2 = blist[maxerr];
		erlast = errmax;
		dqk21_(f, &a1, &b1, &area1, &error1, &resa, &defab1);
		dqk21_(f, &a2, &b2, &area2, &error2, &resa, &defab2);

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		*neval += 42;
		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (defab1 == error1 || defab2 == error2) {
			goto L95;
		}
		if ((d__1 = rlist[maxerr] - area12, fabs(d__1)) > fabs(area12) * 1e-5 ||
			 erro12 < errmax * .99) {
			goto L90;
		}
		if (extrap) {
			++iroff2;
		}
		if (! extrap) {
			++iroff1;
		}
	L90:
		if (*last > 10 && erro12 > errmax) {
			++iroff3;
		}
	L95:
		level[maxerr] = levcur;
		level[*last] = levcur;
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(area);
		errbnd = max(d__1,d__2);

		/*           test for roundoff error and eventually set error flag. */

		if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
			*ier = 2;
		}
		if (iroff2 >= 5) {
			ierro = 3;
		}

		/*           set error flag in the case that the number of */
		/*           subintervals equals limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at a point of the integration range */

		/* Computing MAX */
		d__1 = fabs(a1), d__2 = fabs(b2);
		if (max(d__1,d__2) <= (epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) 
		{
			*ier = 4;
		}

		/*           append the newly-created intervals to the list. */

		if (error2 > error1) {
			goto L100;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L110;
	L100:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine dqpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the subinterval */
		/*           with nrmax-th largest error estimate (to be bisected next). */

	L110:
		dqpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (errsum <= errbnd) {
			goto L190;
		}
		/* ***jump out of do-loop */
		if (*ier != 0) {
			goto L170;
		}
		if (noext) {
			goto L160;
		}
		erlarg -= erlast;
		if (levcur + 1 <= levmax) {
			erlarg += erro12;
		}
		if (extrap) {
			goto L120;
		}

		/*           test whether the interval to be bisected next is the */
		/*           smallest interval. */

		if (level[maxerr] + 1 <= levmax) {
			goto L160;
		}
		extrap = true;
		nrmax = 2;
	L120:
		if (ierro == 3 || erlarg <= ertest) {
			goto L140;
		}

		/*           the smallest interval has the largest error. */
		/*           before bisecting decrease the sum of the errors over */
		/*           the larger intervals (erlarg) and perform extrapolation. */

		id = nrmax;
		jupbnd = *last;
		if (*last > *limit / 2 + 2) {
			jupbnd = *limit + 3 - *last;
		}
		i__1 = jupbnd;
		for (k = id; k <= i__1; ++k) {
			maxerr = iord[nrmax];
			errmax = elist[maxerr];
			/* ***jump out of do-loop */
			if (level[maxerr] + 1 <= levmax) {
				goto L160;
			}
			++nrmax;
			/* L130: */
		}

		/*           perform extrapolation. */

	L140:
		++numrl2;
		rlist2[numrl2 - 1] = area;
		if (numrl2 <= 2) {
			goto L155;
		}
		dqelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
		++ktmin;
		if (ktmin > 5 && *abserr < errsum * .001) {
			*ier = 5;
		}
		if (abseps >= *abserr) {
			goto L150;
		}
		ktmin = 0;
		*abserr = abseps;
		*result = reseps;
		correc = erlarg;
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(reseps);
		ertest = max(d__1,d__2);
		/* ***jump out of do-loop */
		if (*abserr < ertest) {
			goto L170;
		}

		/*           prepare bisection of the smallest interval. */

	L150:
		if (numrl2 == 1) {
			noext = true;
		}
		if (*ier >= 5) {
			goto L170;
		}
	L155:
		maxerr = iord[1];
		errmax = elist[maxerr];
		nrmax = 1;
		extrap = false;
		++levmax;
		erlarg = errsum;
	L160:
		;
	}

	/*           set the final result. */
	/*           --------------------- */


L170:
	if (*abserr == oflow) {
		goto L190;
	}
	if (*ier + ierro == 0) {
		goto L180;
	}
	if (ierro == 3) {
		*abserr += correc;
	}
	if (*ier == 0) {
		*ier = 3;
	}
	if (*result != 0. && area != 0.) {
		goto L175;
	}
	if (*abserr > errsum) {
		goto L190;
	}
	if (area == 0.) {
		goto L210;
	}
	goto L180;
L175:
	if (*abserr / fabs(*result) > errsum / fabs(area)) {
		goto L190;
	}

	/*           test on divergence. */

L180:
	/* Computing MAX */
	d__1 = fabs(*result), d__2 = fabs(area);
	if (ksgn == -1 && max(d__1,d__2) <= resabs * .01) {
		goto L210;
	}
	if (.01 > *result / area || *result / area > 100. || errsum > fabs(area)) {
		*ier = 6;
	}
	goto L210;

	/*           compute global integral sum. */

L190:
	*result = 0.;
	i__2 = *last;
	for (k = 1; k <= i__2; ++k) {
		*result += rlist[k];
		/* L200: */
	}
	*abserr = errsum;
L210:
	if (*ier > 2) {
		--(*ier);
	}
	*result *= sign;
L999:
	return;
} /* dqagpe_ */

void dqagp_(const D_fp& f, const double *a, const double *b, const long *npts2, 
				const double *points, const double *epsabs, const double *epsrel, 
				double *result, double *abserr, long *neval, long *ier, 
				const long *leniw, const long *lenw, long *last, long *iwork, 
				double *work)
{
	long l1, l2, l3, l4, lvl, limit;

	/* ***begin prologue  dqagp */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             singularities at user specified points, */
	/*             extrapolation, globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            break points of the integration interval, where local */
	/*            difficulties of the integrand may occur (e.g. */
	/*            singularities, discontinuities), are provided by the user. */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration */

	/*            npts2  - long */
	/*                     number equal to two more than the number of */
	/*                     user-supplied break points within the integration */
	/*                     range, npts.ge.2. */
	/*                     if npts2.lt.2, the routine will end with ier = 6. */

	/*            points - double precision */
	/*                     vector of dimension npts2, the first (npts2-2) */
	/*                     elements of which are the user provided break */
	/*                     points. if these points do not constitute an */
	/*                     ascending sequence there will be an automatic */
	/*                     sorting. */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine. */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. if */
	/*                             the position of a local difficulty can be */
	/*                             determined (i.e. singularity, */
	/*                             discontinuity within the interval), it */
	/*                             should be supplied to the routine as an */
	/*                             element of the vector points. if necessary */
	/*                             an appropriate special-purpose integrator */
	/*                             must be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. */
	/*                             it is presumed that the requested */
	/*                             tolerance cannot be achieved, and that */
	/*                             the returned result is the best which */
	/*                             can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier.gt.0. */
	/*                         = 6 the input is invalid because */
	/*                             npts2.lt.2 or */
	/*                             break points are specified outside */
	/*                             the integration range or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             result, abserr, neval, last are set to */
	/*                             zero. exept when leniw or lenw or npts2 is */
	/*                             invalid, iwork(1), iwork(limit+1), */
	/*                             work(limit*2+1) and work(limit*3+1) */
	/*                             are set to zero. */
	/*                             work(1) is set to a and work(limit+1) */
	/*                             to b (where limit = (leniw-npts2)/2). */

	/*         dimensioning parameters */
	/*            leniw - long */
	/*                    dimensioning parameter for iwork */
	/*                    leniw determines limit = (leniw-npts2)/2, */
	/*                    which is the maximum number of subintervals in the */
	/*                    partition of the given integration interval (a,b), */
	/*                    leniw.ge.(3*npts2-2). */
	/*                    if leniw.lt.(3*npts2-2), the routine will end with */
	/*                    ier = 6. */

	/*            lenw  - long */
	/*                    dimensioning parameter for work */
	/*                    lenw must be at least leniw*2-npts2. */
	/*                    if lenw.lt.leniw*2-npts2, the routine will end */
	/*                    with ier = 6. */

	/*            last  - long */
	/*                    on return, last equals the number of subintervals */
	/*                    produced in the subdivision process, which */
	/*                    determines the number of significant elements */
	/*                    actually in the work arrays. */

	/*         work arrays */
	/*            iwork - long */
	/*                    vector of dimension at least leniw. on return, */
	/*                    the first k elements of which contain */
	/*                    pointers to the error estimates over the */
	/*                    subintervals, such that work(limit*3+iwork(1)),..., */
	/*                    work(limit*3+iwork(k)) form a decreasing */
	/*                    sequence, with k = last if last.le.(limit/2+2), and */
	/*                    k = limit+1-last otherwise */
	/*                    iwork(limit+1), ...,iwork(limit+last) contain the */
	/*                     subdivision levels of the subintervals, i.e. */
	/*                     if (aa,bb) is a subinterval of (p1,p2) */
	/*                     where p1 as well as p2 is a user-provided */
	/*                     break point or integration limit, then (aa,bb) has */
	/*                     level l if fabs(bb-aa) = fabs(p2-p1)*2**(-l), */
	/*                    iwork(limit*2+1), ..., iwork(limit*2+npts2) have */
	/*                     no significance for the user, */
	/*                    note that limit = (leniw-npts2)/2. */

	/*            work  - double precision */
	/*                    vector of dimension at least lenw */
	/*                    on return */
	/*                    work(1), ..., work(last) contain the left */
	/*                     end points of the subintervals in the */
	/*                     partition of (a,b), */
	/*                    work(limit+1), ..., work(limit+last) contain */
	/*                     the right end points, */
	/*                    work(limit*2+1), ..., work(limit*2+last) contain */
	/*                     the integral approximations over the subintervals, */
	/*                    work(limit*3+1), ..., work(limit*3+last) */
	/*                     contain the corresponding error estimates, */
	/*                    work(limit*4+1), ..., work(limit*4+npts2) */
	/*                     contain the integration limits and the */
	/*                     break points sorted in an ascending sequence. */
	/*                    note that limit = (leniw-npts2)/2. */

	/* ***references  (none) */
	/* ***routines called  dqagpe,xerror */
	/* ***end prologue  dqagp */




	/*         check validity of limit and lenw. */

	/* ***first executable statement  dqagp */
	/* Parameter adjustments */
	--points;
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	if (*leniw < *npts2 * 3 - 2 || *lenw < (*leniw << 1) - *npts2 || *npts2 < 
	    2) {
		goto L10;
	}

	/*         prepare call for dqagpe. */

	limit = (*leniw - *npts2) / 2;
	l1 = limit + 1;
	l2 = limit + l1;
	l3 = limit + l2;
	l4 = limit + l3;

	dqagpe_(f, a, b, npts2, &points[1], epsabs, epsrel, &limit, result, 
			  abserr, neval, ier, &work[1], &work[l1], &work[l2], &work[l3], &
			  work[l4], &iwork[1], &iwork[l1], &iwork[l2], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from dqagp", &c__26, ier, &lvl, 26);
	}
	return;
} /* dqagp_ */

void dqagse_(const D_fp& f, const double *a, const double *b, 
				 const double *epsabs, const double *epsrel, const long *limit, 
				 double *result, 
				 double *abserr, long *neval, long *ier, double *alist__,
				 double *blist, double *rlist, double *elist, long *
				 iord, long *last)
{
	/* System generated locals */
	int i__1, i__2;
	double d__1, d__2;

	/* Local variables */
	long k;
	double a1, a2, b1, b2;
	long id;
	double area;
	double dres;
	long ksgn, nres;
	double area1, area2, area12;
	double small=0., erro12;
	long ierro;
	double defab1, defab2;
	long ktmin, nrmax;
	double oflow, uflow;
	bool noext;
	long iroff1, iroff2, iroff3;
	double res3la[3], error1, error2, rlist2[52];
	long numrl2;
	double defabs, epmach, erlarg=0., abseps, correc=0., errbnd, resabs;
	long jupbnd;
	double erlast, errmax;
	long maxerr;
	double reseps;
	bool extrap;
	double ertest=0., errsum;

	/* ***begin prologue  dqagse */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             (end point) singularities, extrapolation, */
	/*             globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upperbound on the number of subintervals */
	/*                     in the partition of (a,b) */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                         = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more sub- */
	/*                             divisions by increasing the value of limit */
	/*                             (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. if */
	/*                             the position of a local difficulty can be */
	/*                             determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is detec- */
	/*                             ted, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour */
	/*                             occurs at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. */
	/*                             it is presumed that the requested */
	/*                             tolerance cannot be achieved, and that the */
	/*                             returned result is the best which can be */
	/*                             obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier. */
	/*                         = 6 the input is invalid, because */
	/*                             epsabs.le.0 and */
	/*                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28). */
	/*                             result, abserr, neval, last, rlist(1), */
	/*                             iord(1) and elist(1) are set to zero. */
	/*                             alist(1) and blist(1) are set to a and b */
	/*                             respectively. */

	/*            alist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the left end points */
	/*                     of the subintervals in the partition of the */
	/*                     given integration range (a,b) */

	/*            blist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the right end points */
	/*                     of the subintervals in the partition of the given */
	/*                     integration range (a,b) */

	/*            rlist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the integral */
	/*                     approximations on the subintervals */

	/*            elist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the moduli of the */
	/*                     absolute error estimates on the subintervals */

	/*            iord   - long */
	/*                     vector of dimension at least limit, the first k */
	/*                     elements of which are pointers to the */
	/*                     error estimates over the subintervals, */
	/*                     such that elist(iord(1)), ..., elist(iord(k)) */
	/*                     form a decreasing sequence, with k = last */
	/*                     if last.le.(limit/2+2), and k = limit+1-last */
	/*                     otherwise */

	/*            last   - long */
	/*                     number of subintervals actually produced in the */
	/*                     subdivision process */

	/* ***references  (none) */
	/* ***routines called  d1mach,dqelg,dqk21,dqpsrt */
	/* ***end prologue  dqagse */




	/*            the dimension of rlist2 is determined by the value of */
	/*            limexp in subroutine dqelg (rlist2 should be of dimension */
	/*            (limexp+2) at least). */

	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           rlist2    - array of dimension at least limexp+2 containing */
	/*                       the part of the epsilon table which is still */
	/*                       needed for further computations */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest error */
	/*                       estimate */
	/*           errmax    - elist(maxerr) */
	/*           erlast    - error on the interval currently subdivided */
	/*                       (before that subdivision has taken place) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left interval */
	/*           *****2    - variable for the right interval */
	/*           last      - index for subdivision */
	/*           nres      - number of calls to the extrapolation routine */
	/*           numrl2    - number of elements currently in rlist2. if an */
	/*                       appropriate approximation to the compounded */
	/*                       integral has been obtained it is put in */
	/*                       rlist2(numrl2) after numrl2 has been increased */
	/*                       by one. */
	/*           small     - length of the smallest interval considered up */
	/*                       to now, multiplied by 1.5 */
	/*           erlarg    - sum of the errors over the intervals larger */
	/*                       than the smallest interval considered up to now */
	/*           extrap    - bool variable denoting that the routine is */
	/*                       attempting to perform extrapolation i.e. before */
	/*                       subdividing the smallest interval we try to */
	/*                       decrease the value of erlarg. */
	/*           noext     - bool variable denoting that extrapolation */
	/*                       is no longer allowed (true value) */

	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */
	/*           oflow is the largest positive magnitude. */

	/* ***first executable statement  dqagse */
	/* Parameter adjustments */
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = d1mach(c__4);

	/*            test on validity of parameters */
	/*            ------------------------------ */
	*ier = 0;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	alist__[1] = *a;
	blist[1] = *b;
	rlist[1] = 0.;
	elist[1] = 0.;
	/* Computing MAX */
	d__1 = epmach * 50.;
	if (*epsabs <= 0. && *epsrel < max(d__1,5e-29)) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}

	/*           first approximation to the integral */
	/*           ----------------------------------- */

	uflow = d1mach(c__1);
	oflow = d1mach(c__2);
	ierro = 0;
	dqk21_(f, a, b, result, abserr, &defabs, &resabs);

	/*           test on accuracy. */

	dres = fabs(*result);
	/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * dres;
	errbnd = max(d__1,d__2);
	*last = 1;
	rlist[1] = *result;
	elist[1] = *abserr;
	iord[1] = 1;
	if (*abserr <= epmach * 100. * defabs && *abserr > errbnd) {
		*ier = 2;
	}
	if (*limit == 1) {
		*ier = 1;
	}
	if (*ier != 0 || ( *abserr <= errbnd && *abserr != resabs ) || *abserr == 0.) 
	{
		goto L140;
	}

	/*           initialization */
	/*           -------------- */

	rlist2[0] = *result;
	errmax = *abserr;
	maxerr = 1;
	area = *result;
	errsum = *abserr;
	*abserr = oflow;
	nrmax = 1;
	nres = 0;
	numrl2 = 2;
	ktmin = 0;
	extrap = false;
	noext = false;
	iroff1 = 0;
	iroff2 = 0;
	iroff3 = 0;
	ksgn = -1;
	if (dres >= (1. - epmach * 50.) * defabs) {
		ksgn = 1;
	}

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limit;
	for (*last = 2; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with the nrmax-th largest error */
		/*           estimate. */

		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5;
		a2 = b1;
		b2 = blist[maxerr];
		erlast = errmax;
		dqk21_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		dqk21_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (defab1 == error1 || defab2 == error2) {
			goto L15;
		}
		if ((d__1 = rlist[maxerr] - area12, fabs(d__1)) > fabs(area12) * 1e-5 ||
			 erro12 < errmax * .99) {
			goto L10;
		}
		if (extrap) {
			++iroff2;
		}
		if (! extrap) {
			++iroff1;
		}
	L10:
		if (*last > 10 && erro12 > errmax) {
			++iroff3;
		}
	L15:
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(area);
		errbnd = max(d__1,d__2);

		/*           test for roundoff error and eventually set error flag. */

		if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
			*ier = 2;
		}
		if (iroff2 >= 5) {
			ierro = 3;
		}

		/*           set error flag in the case that the number of subintervals */
		/*           equals limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at a point of the integration range. */

		/* Computing MAX */
		d__1 = fabs(a1), d__2 = fabs(b2);
		if (max(d__1,d__2) <= (epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) 
		{
			*ier = 4;
		}

		/*           append the newly-created intervals to the list. */

		if (error2 > error1) {
			goto L20;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L30;
	L20:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine dqpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the subinterval */
		/*           with nrmax-th largest error estimate (to be bisected next). */

	L30:
		dqpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (errsum <= errbnd) {
			goto L115;
		}
		/* ***jump out of do-loop */
		if (*ier != 0) {
			goto L100;
		}
		if (*last == 2) {
			goto L80;
		}
		if (noext) {
			goto L90;
		}
		erlarg -= erlast;
		if ((d__1 = b1 - a1, fabs(d__1)) > small) {
			erlarg += erro12;
		}
		if (extrap) {
			goto L40;
		}

		/*           test whether the interval to be bisected next is the */
		/*           smallest interval. */

		if ((d__1 = blist[maxerr] - alist__[maxerr], fabs(d__1)) > small) {
			goto L90;
		}
		extrap = true;
		nrmax = 2;
	L40:
		if (ierro == 3 || erlarg <= ertest) {
			goto L60;
		}

		/*           the smallest interval has the largest error. */
		/*           before bisecting decrease the sum of the errors over the */
		/*           larger intervals (erlarg) and perform extrapolation. */

		id = nrmax;
		jupbnd = *last;
		if (*last > *limit / 2 + 2) {
			jupbnd = *limit + 3 - *last;
		}
		i__2 = jupbnd;
		for (k = id; k <= i__2; ++k) {
			maxerr = iord[nrmax];
			errmax = elist[maxerr];
			/* ***jump out of do-loop */
			if ((d__1 = blist[maxerr] - alist__[maxerr], fabs(d__1)) > small) {
				goto L90;
			}
			++nrmax;
			/* L50: */
		}

		/*           perform extrapolation. */

	L60:
		++numrl2;
		rlist2[numrl2 - 1] = area;
		dqelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
		++ktmin;
		if (ktmin > 5 && *abserr < errsum * .001) {
			*ier = 5;
		}
		if (abseps >= *abserr) {
			goto L70;
		}
		ktmin = 0;
		*abserr = abseps;
		*result = reseps;
		correc = erlarg;
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(reseps);
		ertest = max(d__1,d__2);
		/* ***jump out of do-loop */
		if (*abserr <= ertest) {
			goto L100;
		}

		/*           prepare bisection of the smallest interval. */

	L70:
		if (numrl2 == 1) {
			noext = true;
		}
		if (*ier == 5) {
			goto L100;
		}
		maxerr = iord[1];
		errmax = elist[maxerr];
		nrmax = 1;
		extrap = false;
		small *= .5;
		erlarg = errsum;
		goto L90;
	L80:
		small = (d__1 = *b - *a, fabs(d__1)) * .375;
		erlarg = errsum;
		ertest = errbnd;
		rlist2[1] = area;
	L90:
		;
	}

	/*           set final result and error estimate. */
	/*           ------------------------------------ */

L100:
	if (*abserr == oflow) {
		goto L115;
	}
	if (*ier + ierro == 0) {
		goto L110;
	}
	if (ierro == 3) {
		*abserr += correc;
	}
	if (*ier == 0) {
		*ier = 3;
	}
	if (*result != 0. && area != 0.) {
		goto L105;
	}
	if (*abserr > errsum) {
		goto L115;
	}
	if (area == 0.) {
		goto L130;
	}
	goto L110;
L105:
	if (*abserr / fabs(*result) > errsum / fabs(area)) {
		goto L115;
	}

	/*           test on divergence. */

L110:
	/* Computing MAX */
	d__1 = fabs(*result), d__2 = fabs(area);
	if (ksgn == -1 && max(d__1,d__2) <= defabs * .01) {
		goto L130;
	}
	if (.01 > *result / area || *result / area > 100. || errsum > fabs(area)) {
		*ier = 6;
	}
	goto L130;

	/*           compute global integral sum. */

L115:
	*result = 0.;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L120: */
	}
	*abserr = errsum;
L130:
	if (*ier > 2) {
		--(*ier);
	}
L140:
	*neval = *last * 42 - 21;
L999:
	return;
} /* dqagse_ */

void dqags_(const D_fp& f, const double *a, const double *b, const double *epsabs, 
				const double *epsrel, double *result, double *abserr, 
				long *neval, long *ier, const long *limit, const long *lenw, long *
				last, long *iwork, double *work)
{
	long l1, l2, l3, lvl;

	/* ***begin prologue  dqags */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             (end-point) singularities, extrapolation, */
	/*             globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & prog. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral  i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        double precision version */


	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more sub- */
	/*                             divisions by increasing the value of limit */
	/*                             (and taking the according dimension */
	/*                             adjustments into account. however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. if */
	/*                             the position of a local difficulty can be */
	/*                             determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is detec- */
	/*                             ted, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour */
	/*                             occurs at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. it is presumed that */
	/*                             the requested tolerance cannot be */
	/*                             achieved, and that the returned result is */
	/*                             the best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28) */
	/*                             or limit.lt.1 or lenw.lt.limit*4. */
	/*                             result, abserr, neval, last are set to */
	/*                             zero.except when limit or lenw is invalid, */
	/*                             iwork(1), work(limit*2+1) and */
	/*                             work(limit*3+1) are set to zero, work(1) */
	/*                             is set to a and work(limit+1) to b. */

	/*         dimensioning parameters */
	/*            limit - long */
	/*                    dimensioning parameter for iwork */
	/*                    limit determines the maximum number of subintervals */
	/*                    in the partition of the given integration interval */
	/*                    (a,b), limit.ge.1. */
	/*                    if limit.lt.1, the routine will end with ier = 6. */

	/*            lenw  - long */
	/*                    dimensioning parameter for work */
	/*                    lenw must be at least limit*4. */
	/*                    if lenw.lt.limit*4, the routine will end */
	/*                    with ier = 6. */

	/*            last  - long */
	/*                    on return, last equals the number of subintervals */
	/*                    produced in the subdivision process, detemines the */
	/*                    number of significant elements actually in the work */
	/*                    arrays. */

	/*         work arrays */
	/*            iwork - long */
	/*                    vector of dimension at least limit, the first k */
	/*                    elements of which contain pointers */
	/*                    to the error estimates over the subintervals */
	/*                    such that work(limit*3+iwork(1)),... , */
	/*                    work(limit*3+iwork(k)) form a decreasing */
	/*                    sequence, with k = last if last.le.(limit/2+2), */
	/*                    and k = limit+1-last otherwise */

	/*            work  - double precision */
	/*                    vector of dimension at least lenw */
	/*                    on return */
	/*                    work(1), ..., work(last) contain the left */
	/*                     end-points of the subintervals in the */
	/*                     partition of (a,b), */
	/*                    work(limit+1), ..., work(limit+last) contain */
	/*                     the right end-points, */
	/*                    work(limit*2+1), ..., work(limit*2+last) contain */
	/*                     the integral approximations over the subintervals, */
	/*                    work(limit*3+1), ..., work(limit*3+last) */
	/*                     contain the error estimates. */

	/* ***references  (none) */
	/* ***routines called  dqagse,xerror */
	/* ***end prologue  dqags */





	/*         check validity of limit and lenw. */

	/* ***first executable statement  dqags */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	if (*limit < 1 || *lenw < *limit << 2) {
		goto L10;
	}

	/*         prepare call for dqagse. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;

	dqagse_(f, a, b, epsabs, epsrel, limit, result, abserr, neval, ier, 
			  &work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from dqags", &c__26, ier, &lvl, 26);
	}
	return;
} /* dqags_ */

void dqawce_(const D_fp& f, const double *a, const double *b, const double *c__, 
				 const double *epsabs, const double *epsrel, const long *limit, 
				 double *result, double *abserr, long *neval, long *ier, 
				 double *alist__, double *blist, double *rlist, double 
				 *elist, long *iord, long *last)
{
	/* System generated locals */
	int i__1;
	double d__1, d__2;

	/* Local variables */
	long k;
	double a1, a2, b1, b2, aa, bb;
	long nev;
	double area, area1, area2, area12;
	double erro12;
	long krule, nrmax;
	double uflow;
	long iroff1, iroff2;
	double error1, error2, epmach, errbnd, errmax;
	long maxerr;
	double errsum;

	/* ***begin prologue  dqawce */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1,j4 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             cauchy principal value, clenshaw-curtis method */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***  purpose  the routine calculates an approximation result to a */
	/*              cauchy principal value i = integral of f*w over (a,b) */
	/*              (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying */
	/*              following claim for accuracy */
	/*              fabs(i-result).le.max(epsabs,epsrel*fabs(i)) */
	/* ***description */

	/*        computation of a cauchy principal value */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration */

	/*            c      - double precision */
	/*                     parameter in the weight function, c.ne.a, c.ne.b */
	/*                     if c = a or c = b, the routine will end with */
	/*                     ier = 6. */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subintervals */
	/*                     in the partition of (a,b), limit.ge.1 */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more sub- */
	/*                             divisions by increasing the value of */
	/*                             limit. however, if this yields no */
	/*                             improvement it is advised to analyze the */
	/*                             the integrand, in order to determine the */
	/*                             the integration difficulties. if the */
	/*                             position of a local difficulty can be */
	/*                             determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling */
	/*                             appropriate integrators on the subranges. */
	/*                         = 2 the occurrence of roundoff error is detec- */
	/*                             ted, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour */
	/*                             occurs at some interior points of */
	/*                             the integration interval. */
	/*                         = 6 the input is invalid, because */
	/*                             c = a or c = b or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or limit.lt.1. */
	/*                             result, abserr, neval, rlist(1), elist(1), */
	/*                             iord(1) and last are set to zero. alist(1) */
	/*                             and blist(1) are set to a and b */
	/*                             respectively. */

	/*            alist   - double precision */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the left */
	/*                      end points of the subintervals in the partition */
	/*                      of the given integration range (a,b) */

	/*            blist   - double precision */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the right */
	/*                      end points of the subintervals in the partition */
	/*                      of the given integration range (a,b) */

	/*            rlist   - double precision */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the integral */
	/*                      approximations on the subintervals */

	/*            elist   - double precision */
	/*                      vector of dimension limit, the first  last */
	/*                      elements of which are the moduli of the absolute */
	/*                      error estimates on the subintervals */

	/*            iord    - long */
	/*                      vector of dimension at least limit, the first k */
	/*                      elements of which are pointers to the error */
	/*                      estimates over the subintervals, so that */
	/*                      elist(iord(1)), ..., elist(iord(k)) with k = last */
	/*                      if last.le.(limit/2+2), and k = limit+1-last */
	/*                      otherwise, form a decreasing sequence */

	/*            last    - long */
	/*                      number of subintervals actually produced in */
	/*                      the subdivision process */

	/* ***references  (none) */
	/* ***routines called  d1mach,dqc25c,dqpsrt */
	/* ***end prologue  dqawce */




	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest */
	/*                       error estimate */
	/*           errmax    - elist(maxerr) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */


	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  dqawce */
	/* Parameter adjustments */
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);


	/*           test on validity of parameters */
	/*           ------------------------------ */

	*ier = 6;
	*neval = 0;
	*last = 0;
	alist__[1] = *a;
	blist[1] = *b;
	rlist[1] = 0.;
	elist[1] = 0.;
	iord[1] = 0;
	*result = 0.;
	*abserr = 0.;
	/* Computing MAX */
	d__1 = epmach * 50.;
	if (*c__ == *a || *c__ == *b || ( *epsabs <= 0. && *epsrel < max(d__1,5e-29) )
		) {
		goto L999;
	}

	/*           first approximation to the integral */
	/*           ----------------------------------- */

	aa = *a;
	bb = *b;
	if (*a <= *b) {
		goto L10;
	}
	aa = *b;
	bb = *a;
L10:
	*ier = 0;
	krule = 1;
	dqc25c_(f, &aa, &bb, c__, result, abserr, &krule, neval);
	*last = 1;
	rlist[1] = *result;
	elist[1] = *abserr;
	iord[1] = 1;
	alist__[1] = *a;
	blist[1] = *b;

	/*           test on accuracy */

	/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * fabs(*result);
	errbnd = max(d__1,d__2);
	if (*limit == 1) {
		*ier = 1;
	}
	/* Computing MIN */
	d__1 = fabs(*result) * .01;
	if (*abserr < min(d__1,errbnd) || *ier == 1) {
		goto L70;
	}

	/*           initialization */
	/*           -------------- */

	alist__[1] = aa;
	blist[1] = bb;
	rlist[1] = *result;
	errmax = *abserr;
	maxerr = 1;
	area = *result;
	errsum = *abserr;
	nrmax = 1;
	iroff1 = 0;
	iroff2 = 0;

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limit;
	for (*last = 2; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with nrmax-th largest */
		/*           error estimate. */

		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5;
		b2 = blist[maxerr];
		if (*c__ <= b1 && *c__ > a1) {
			b1 = (*c__ + b2) * .5;
		}
		if (*c__ > b1 && *c__ < b2) {
			b1 = (a1 + *c__) * .5;
		}
		a2 = b1;
		krule = 2;
		dqc25c_(f, &a1, &b1, c__, &area1, &error1, &krule, &nev);
		*neval += nev;
		dqc25c_(f, &a2, &b2, c__, &area2, &error2, &krule, &nev);
		*neval += nev;

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if ((d__1 = rlist[maxerr] - area12, fabs(d__1)) < fabs(area12) * 1e-5 &&
			 erro12 >= errmax * .99 && krule == 0) {
			++iroff1;
		}
		if (*last > 10 && erro12 > errmax && krule == 0) {
			++iroff2;
		}
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(area);
		errbnd = max(d__1,d__2);
		if (errsum <= errbnd) {
			goto L15;
		}

		/*           test for roundoff error and eventually set error flag. */

		if (iroff1 >= 6 && iroff2 > 20) {
			*ier = 2;
		}

		/*           set error flag in the case that number of interval */
		/*           bisections exceeds limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at a point of the integration range. */

		/* Computing MAX */
		d__1 = fabs(a1), d__2 = fabs(b2);
		if (max(d__1,d__2) <= (epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) 
		{
			*ier = 3;
		}

		/*           append the newly-created intervals to the list. */

	L15:
		if (error2 > error1) {
			goto L20;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L30;
	L20:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine dqpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the subinterval */
		/*           with nrmax-th largest error estimate (to be bisected next). */

	L30:
		dqpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (*ier != 0 || errsum <= errbnd) {
			goto L50;
		}
		/* L40: */
	}

	/*           compute final result. */
	/*           --------------------- */

L50:
	*result = 0.;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L60: */
	}
	*abserr = errsum;
L70:
	if (aa == *b) {
		*result = -(*result);
	}
L999:
	return;
} /* dqawce_ */

void dqawc_(const D_fp& f, const double *a, const double *b, const double *c__, 
				const double *epsabs, const double *epsrel, double *result, 
				double *abserr, long *neval, long *ier, long *limit, 
				const long *lenw, long *last, long *iwork, double *work)
{
	long l1, l2, l3, lvl;

	/* ***begin prologue  dqawc */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1,j4 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             cauchy principal value, */
	/*             clenshaw-curtis, globally adaptive */
	/* ***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a */
	/*            cauchy principal value i = integral of f*w over (a,b) */
	/*            (w(x) = 1/((x-c), c.ne.a, c.ne.b), hopefully satisfying */
	/*            following claim for accuracy */
	/*            fabs(i-result).le.max(epsabe,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of a cauchy principal value */
	/*        standard fortran subroutine */
	/*        double precision version */


	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     under limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration */

	/*            c      - parameter in the weight function, c.ne.a, c.ne.b. */
	/*                     if c = a or c = b, the routine will end with */
	/*                     ier = 6 . */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate or the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more sub- */
	/*                             divisions by increasing the value of limit */
	/*                             (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. */
	/*                             if the position of a local difficulty */
	/*                             can be determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling */
	/*                             appropriate integrators on the subranges. */
	/*                         = 2 the occurrence of roundoff error is detec- */
	/*                             ted, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 6 the input is invalid, because */
	/*                             c = a or c = b or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or limit.lt.1 or lenw.lt.limit*4. */
	/*                             result, abserr, neval, last are set to */
	/*                             zero. exept when lenw or limit is invalid, */
	/*                             iwork(1), work(limit*2+1) and */
	/*                             work(limit*3+1) are set to zero, work(1) */
	/*                             is set to a and work(limit+1) to b. */

	/*         dimensioning parameters */
	/*            limit - long */
	/*                    dimensioning parameter for iwork */
	/*                    limit determines the maximum number of subintervals */
	/*                    in the partition of the given integration interval */
	/*                    (a,b), limit.ge.1. */
	/*                    if limit.lt.1, the routine will end with ier = 6. */

	/*           lenw   - long */
	/*                    dimensioning parameter for work */
	/*                    lenw must be at least limit*4. */
	/*                    if lenw.lt.limit*4, the routine will end with */
	/*                    ier = 6. */

	/*            last  - long */
	/*                    on return, last equals the number of subintervals */
	/*                    produced in the subdivision process, which */
	/*                    determines the number of significant elements */
	/*                    actually in the work arrays. */

	/*         work arrays */
	/*            iwork - long */
	/*                    vector of dimension at least limit, the first k */
	/*                    elements of which contain pointers */
	/*                    to the error estimates over the subintervals, */
	/*                    such that work(limit*3+iwork(1)), ... , */
	/*                    work(limit*3+iwork(k)) form a decreasing */
	/*                    sequence, with k = last if last.le.(limit/2+2), */
	/*                    and k = limit+1-last otherwise */

	/*            work  - double precision */
	/*                    vector of dimension at least lenw */
	/*                    on return */
	/*                    work(1), ..., work(last) contain the left */
	/*                     end points of the subintervals in the */
	/*                     partition of (a,b), */
	/*                    work(limit+1), ..., work(limit+last) contain */
	/*                     the right end points, */
	/*                    work(limit*2+1), ..., work(limit*2+last) contain */
	/*                     the integral approximations over the subintervals, */
	/*                    work(limit*3+1), ..., work(limit*3+last) */
	/*                     contain the error estimates. */

	/* ***references  (none) */
	/* ***routines called  dqawce,xerror */
	/* ***end prologue  dqawc */




	/*         check validity of limit and lenw. */

	/* ***first executable statement  dqawc */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	if (*limit < 1 || *lenw < *limit << 2) {
		goto L10;
	}

	/*         prepare call for dqawce. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;
	dqawce_(f, a, b, c__, epsabs, epsrel, limit, result, abserr, neval, 
			  ier, &work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from dqawc", &c__26, ier, &lvl, 26);
	}
	return;
} /* dqawc_ */

void dqawfe_(const D_fp& f, const double *a, const double *omega, 
				 const long *integr, const double *epsabs, const long *limlst, 
				 const long *limit, const long *maxp1, 
				 double *result, double *abserr, long *neval, 
				 long *ier, double *rslst, double *erlst, long *ierlst, 
				 long *lst, double *alist__, double *blist, 
				 double *rlist, double *elist, long *iord, long *nnlog, 
				 double *chebmo)
{
	/* Initialized data */

	double p = .9;
	double pi = 3.1415926535897932384626433832795;

	/* System generated locals */
	int chebmo_dim1, chebmo_offset, i__1;
	double d__1, d__2;

	/* Local variables */
	long l;
	double c1, c2, p1, dl, ep;
	long ll=0;
	double drl=0., eps;
	long nev;
	double fact, epsa;
	long last, nres;
	double psum[52];
	double cycle;
	long ktmin;
	double uflow;
	double res3la[3];
	long numrl2;
	double abseps, correc;
	long momcom;
	double reseps, errsum;

	/* ***begin prologue  dqawfe */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a3a1 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             fourier integrals, */
	/*             integration between zeros with dqawoe, */
	/*             convergence acceleration with dqelg */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           dedoncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a */
	/*            given fourier integal */
	/*            i = integral of f(x)*w(x) over (a,infinity) */
	/*            where w(x)=cos(omega*x) or w(x)=sin(omega*x), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.epsabs. */
	/* ***description */

	/*        computation of fourier integrals */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to */
	/*                     be declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            omega  - double precision */
	/*                     parameter in the weight function */

	/*            integr - long */
	/*                     indicates which weight function is used */
	/*                     integr = 1      w(x) = cos(omega*x) */
	/*                     integr = 2      w(x) = sin(omega*x) */
	/*                     if integr.ne.1.and.integr.ne.2, the routine will */
	/*                     end with ier = 6. */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested, epsabs.gt.0 */
	/*                     if epsabs.le.0, the routine will end with ier = 6. */

	/*            limlst - long */
	/*                     limlst gives an upper bound on the number of */
	/*                     cycles, limlst.ge.1. */
	/*                     if limlst.lt.3, the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subintervals */
	/*                     allowed in the partition of each cycle, limit.ge.1 */
	/*                     each cycle, limit.ge.1. */

	/*            maxp1  - long */
	/*                     gives an upper bound on the number of */
	/*                     chebyshev moments which can be stored, i.e. */
	/*                     for the intervals of lengths fabs(b-a)*2**(-l), */
	/*                     l=0,1, ..., maxp1-2, maxp1.ge.1 */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral x */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - ier = 0 normal and reliable termination of */
	/*                             the routine. it is assumed that the */
	/*                             requested accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine. the */
	/*                             estimates for integral and error are less */
	/*                             reliable. it is assumed that the requested */
	/*                             accuracy has not been achieved. */
	/*            error messages */
	/*                    if omega.ne.0 */
	/*                     ier = 1 maximum number of  cycles  allowed */
	/*                             has been achieved., i.e. of subintervals */
	/*                             (a+(k-1)c,a+kc) where */
	/*                             c = (2*int(fabs(omega))+1)*pi/fabs(omega), */
	/*                             for k = 1, 2, ..., lst. */
	/*                             one can allow more cycles by increasing */
	/*                             the value of limlst (and taking the */
	/*                             according dimension adjustments into */
	/*                             account). */
	/*                             examine the array iwork which contains */
	/*                             the error flags on the cycles, in order to */
	/*                             look for eventual local integration */
	/*                             difficulties. if the position of a local */
	/*                             difficulty can be determined (e.g. */
	/*                             singularity, discontinuity within the */
	/*                             interval) one will probably gain from */
	/*                             splitting up the interval at this polong */
	/*                             and calling appropriate integrators on */
	/*                             the subranges. */
	/*                         = 4 the extrapolation table constructed for */
	/*                             convergence acceleration of the series */
	/*                             formed by the integral contributions over */
	/*                             the cycles, does not converge to within */
	/*                             the requested accuracy. as in the case of */
	/*                             ier = 1, it is advised to examine the */
	/*                             array iwork which contains the error */
	/*                             flags on the cycles. */
	/*                         = 6 the input is invalid because */
	/*                             (integr.ne.1 and integr.ne.2) or */
	/*                              epsabs.le.0 or limlst.lt.3. */
	/*                              result, abserr, neval, lst are set */
	/*                              to zero. */
	/*                         = 7 bad integrand behaviour occurs within one */
	/*                             or more of the cycles. location and type */
	/*                             of the difficulty involved can be */
	/*                             determined from the vector ierlst. here */
	/*                             lst is the number of cycles actually */
	/*                             needed (see below). */
	/*                             ierlst(k) = 1 the maximum number of */
	/*                                           subdivisions (= limit) has */
	/*                                           been achieved on the k th */
	/*                                           cycle. */
	/*                                       = 2 occurrence of roundoff error */
	/*                                           is detected and prevents the */
	/*                                           tolerance imposed on the */
	/*                                           k th cycle, from being */
	/*                                           achieved. */
	/*                                       = 3 extremely bad integrand */
	/*                                           behaviour occurs at some */
	/*                                           points of the k th cycle. */
	/*                                       = 4 the integration procedure */
	/*                                           over the k th cycle does */
	/*                                           not converge (to within the */
	/*                                           required accuracy) due to */
	/*                                           roundoff in the */
	/*                                           extrapolation procedure */
	/*                                           invoked on this cycle. it */
	/*                                           is assumed that the result */
	/*                                           on this interval is the */
	/*                                           best which can be obtained. */
	/*                                       = 5 the integral over the k th */
	/*                                           cycle is probably divergent */
	/*                                           or slowly convergent. it */
	/*                                           must be noted that */
	/*                                           divergence can occur with */
	/*                                           any other value of */
	/*                                           ierlst(k). */
	/*                    if omega = 0 and integr = 1, */
	/*                    the integral is calculated by means of dqagie */
	/*                    and ier = ierlst(1) (with meaning as described */
	/*                    for ierlst(k), k = 1). */

	/*            rslst  - double precision */
	/*                     vector of dimension at least limlst */
	/*                     rslst(k) contains the integral contribution */
	/*                     over the interval (a+(k-1)c,a+kc) where */
	/*                     c = (2*int(fabs(omega))+1)*pi/fabs(omega), */
	/*                     k = 1, 2, ..., lst. */
	/*                     note that, if omega = 0, rslst(1) contains */
	/*                     the value of the integral over (a,infinity). */

	/*            erlst  - double precision */
	/*                     vector of dimension at least limlst */
	/*                     erlst(k) contains the error estimate corresponding */
	/*                     with rslst(k). */

	/*            ierlst - long */
	/*                     vector of dimension at least limlst */
	/*                     ierlst(k) contains the error flag corresponding */
	/*                     with rslst(k). for the meaning of the local error */
	/*                     flags see description of output parameter ier. */

	/*            lst    - long */
	/*                     number of subintervals needed for the integration */
	/*                     if omega = 0 then lst is set to 1. */

	/*            alist, blist, rlist, elist - double precision */
	/*                     vector of dimension at least limit, */

	/*            iord, nnlog - long */
	/*                     vector of dimension at least limit, providing */
	/*                     space for the quantities needed in the subdivision */
	/*                     process of each cycle */

	/*            chebmo - double precision */
	/*                     array of dimension at least (maxp1,25), providing */
	/*                     space for the chebyshev moments needed within the */
	/*                     cycles */

	/* ***references  (none) */
	/* ***routines called  d1mach,dqagie,dqawoe,dqelg */
	/* ***end prologue  dqawfe */





	/*            the dimension of  psum  is determined by the value of */
	/*            limexp in subroutine dqelg (psum must be of dimension */
	/*            (limexp+2) at least). */

	/*           list of major variables */
	/*           ----------------------- */

	/*           c1, c2    - end points of subinterval (of length cycle) */
	/*           cycle     - (2*int(fabs(omega))+1)*pi/fabs(omega) */
	/*           psum      - vector of dimension at least (limexp+2) */
	/*                       (see routine dqelg) */
	/*                       psum contains the part of the epsilon table */
	/*                       which is still needed for further computations. */
	/*                       each element of psum is a partial sum of the */
	/*                       series which should sum to the value of the */
	/*                       integral. */
	/*           errsum    - sum of error estimates over the subintervals, */
	/*                       calculated cumulatively */
	/*           epsa      - absolute tolerance requested over current */
	/*                       subinterval */
	/*           chebmo    - array containing the modified chebyshev */
	/*                       moments (see also routine dqc25f) */

	/* Parameter adjustments */
	--ierlst;
	--erlst;
	--rslst;
	--nnlog;
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;
	chebmo_dim1 = *maxp1;
	chebmo_offset = 1 + chebmo_dim1;
	chebmo -= chebmo_offset;

	/* Function Body */

	/*           test on validity of parameters */
	/*           ------------------------------ */

	/* ***first executable statement  dqawfe */
	*result = 0.;
	*abserr = 0.;
	*neval = 0;
	*lst = 0;
	*ier = 0;
	if ( (*integr != 1 && *integr != 2 ) || *epsabs <= 0. || *limlst < 3) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}
	if (*omega != 0.) {
		goto L10;
	}

	/*           integration by dqagie if omega is zero */
	/*           -------------------------------------- */

	if (*integr == 1) {
		dqagie_(f, &c_b20, &c__1, epsabs, &c_b20, limit, result, abserr,
				  neval, ier, &alist__[1], &blist[1], &rlist[1], &elist[1], &
				  iord[1], &last);
	}
	rslst[1] = *result;
	erlst[1] = *abserr;
	ierlst[1] = *ier;
	*lst = 1;
	goto L999;

	/*           initializations */
	/*           --------------- */

L10:
	l = (int) fabs(*omega);
	dl = (double) ((l << 1) + 1);
	cycle = dl * pi / fabs(*omega);
	*ier = 0;
	ktmin = 0;
	*neval = 0;
	numrl2 = 0;
	nres = 0;
	c1 = *a;
	c2 = cycle + *a;
	p1 = 1. - p;
	uflow = d1mach(c__1);
	eps = *epsabs;
	if (*epsabs > uflow / p1) {
		eps = *epsabs * p1;
	}
	ep = eps;
	fact = 1.;
	correc = 0.;
	*abserr = 0.;
	errsum = 0.;

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limlst;
	for (*lst = 1; *lst <= i__1; ++(*lst)) {

		/*           integrate over current subinterval. */

		epsa = eps * fact;
		dqawoe_(f, &c1, &c2, omega, integr, &epsa, &c_b20, limit, lst, 
				  maxp1, &rslst[*lst], &erlst[*lst], &nev, &ierlst[*lst], &last,
				  &alist__[1], &blist[1], &rlist[1], &elist[1], &iord[1], &
				  nnlog[1], &momcom, &chebmo[chebmo_offset]);
		*neval += nev;
		fact *= p;
		errsum += erlst[*lst];
		drl = (d__1 = rslst[*lst], fabs(d__1)) * 50.;

		/*           test on accuracy with partial sum */

		if (errsum + drl <= *epsabs && *lst >= 6) {
			goto L80;
		}
		/* Computing MAX */
		d__1 = correc, d__2 = erlst[*lst];
		correc = max(d__1,d__2);
		if (ierlst[*lst] != 0) {
			/* Computing MAX */
			d__1 = ep, d__2 = correc * p1;
			eps = max(d__1,d__2);
		}
		if (ierlst[*lst] != 0) {
			*ier = 7;
		}
		if (*ier == 7 && errsum + drl <= correc * 10. && *lst > 5) {
			goto L80;
		}
		++numrl2;
		if (*lst > 1) {
			goto L20;
		}
		psum[0] = rslst[1];
		goto L40;
	L20:
		psum[numrl2 - 1] = psum[ll - 1] + rslst[*lst];
		if (*lst == 2) {
			goto L40;
		}

		/*           test on maximum number of subintervals */

		if (*lst == *limlst) {
			*ier = 1;
		}

		/*           perform new extrapolation */

		dqelg_(&numrl2, psum, &reseps, &abseps, res3la, &nres);

		/*           test whether extrapolated result is influenced by roundoff */

		++ktmin;
		if (ktmin >= 15 && *abserr <= (errsum + drl) * .001) {
			*ier = 4;
		}
		if (abseps > *abserr && *lst != 3) {
			goto L30;
		}
		*abserr = abseps;
		*result = reseps;
		ktmin = 0;

		/*           if ier is not 0, check whether direct result (partial sum) */
		/*           or extrapolated result yields the best integral */
		/*           approximation */

		if (*abserr + correc * 10. <= *epsabs || 
			 ( *abserr <= *epsabs && correc * 10. >= *epsabs ) ) {
			goto L60;
		}
	L30:
		if (*ier != 0 && *ier != 7) {
			goto L60;
		}
	L40:
		ll = numrl2;
		c1 = c2;
		c2 += cycle;
		/* L50: */
	}

	/*         set final result and error estimate */
	/*         ----------------------------------- */

L60:
	*abserr += correc * 10.;
	if (*ier == 0) {
		goto L999;
	}
	if (*result != 0. && psum[numrl2 - 1] != 0.) {
		goto L70;
	}
	if (*abserr > errsum) {
		goto L80;
	}
	if (psum[numrl2 - 1] == 0.) {
		goto L999;
	}
L70:
	if (*abserr / fabs(*result) > (errsum + drl) / (d__1 = psum[numrl2 - 1], 
																	fabs(d__1))) {
		goto L80;
	}
	if (*ier >= 1 && *ier != 7) {
		*abserr += drl;
	}
	goto L999;
L80:
	*result = psum[numrl2 - 1];
	*abserr = errsum + drl;
L999:
	return;
} /* dqawfe_ */

void dqawf_(const D_fp& f, const double *a, const double *omega, const long *integr, 
				const double *epsabs, double *result, double *abserr, 
				long *neval, long *ier, long *limlst, long *lst, const long *
				leniw, const long *maxp1, const long *lenw, long *iwork, double *
				work)
{
	long l1, l2, l3, l4, l5, l6, ll2, lvl, limit;

	/* ***begin prologue  dqawf */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a3a1 */
	/* ***keywords  automatic integrator, special-purpose,fourier */
	/*             integral, integration between zeros with dqawoe, */
	/*             convergence acceleration with dqelg */
	/* ***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            fourier integral i=integral of f(x)*w(x) over (a,infinity) */
	/*            where w(x) = cos(omega*x) or w(x) = sin(omega*x). */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.epsabs. */
	/* ***description */

	/*        computation of fourier integrals */
	/*        standard fortran subroutine */
	/*        double precision version */


	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            omega  - double precision */
	/*                     parameter in the integrand weight function */

	/*            integr - long */
	/*                     indicates which of the weight functions is used */
	/*                     integr = 1      w(x) = cos(omega*x) */
	/*                     integr = 2      w(x) = sin(omega*x) */
	/*                     if integr.ne.1.and.integr.ne.2, the routine */
	/*                     will end with ier = 6. */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested, epsabs.gt.0. */
	/*                     if epsabs.le.0, the routine will end with ier = 6. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine. */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                    if omega.ne.0 */
	/*                     ier = 1 maximum number of cycles allowed */
	/*                             has been achieved, i.e. of subintervals */
	/*                             (a+(k-1)c,a+kc) where */
	/*                             c = (2*int(fabs(omega))+1)*pi/fabs(omega), */
	/*                             for k = 1, 2, ..., lst. */
	/*                             one can allow more cycles by increasing */
	/*                             the value of limlst (and taking the */
	/*                             according dimension adjustments into */
	/*                             account). examine the array iwork which */
	/*                             contains the error flags on the cycles, in */
	/*                             order to look for eventual local */
	/*                             integration difficulties. */
	/*                             if the position of a local difficulty */
	/*                             can be determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling */
	/*                             appropriate integrators on the subranges. */
	/*                         = 4 the extrapolation table constructed for */
	/*                             convergence accelaration of the series */
	/*                             formed by the integral contributions over */
	/*                             the cycles, does not converge to within */
	/*                             the requested accuracy. */
	/*                             as in the case of ier = 1, it is advised */
	/*                             to examine the array iwork which contains */
	/*                             the error flags on the cycles. */
	/*                         = 6 the input is invalid because */
	/*                             (integr.ne.1 and integr.ne.2) or */
	/*                              epsabs.le.0 or limlst.lt.1 or */
	/*                              leniw.lt.(limlst+2) or maxp1.lt.1 or */
	/*                              lenw.lt.(leniw*2+maxp1*25). */
	/*                              result, abserr, neval, lst are set to */
	/*                              zero. */
	/*                         = 7 bad integrand behaviour occurs within */
	/*                             one or more of the cycles. location and */
	/*                             type of the difficulty involved can be */
	/*                             determined from the first lst elements of */
	/*                             vector iwork.  here lst is the number of */
	/*                             cycles actually needed (see below). */
	/*                             iwork(k) = 1 the maximum number of */
	/*                                          subdivisions (=(leniw-limlst) */
	/*                                          /2) has been achieved on the */
	/*                                          k th cycle. */
	/*                                      = 2 occurrence of roundoff error */
	/*                                          is detected and prevents the */
	/*                                          tolerance imposed on the k th */
	/*                                          cycle, from being achieved */
	/*                                          on this cycle. */
	/*                                      = 3 extremely bad integrand */
	/*                                          behaviour occurs at some */
	/*                                          points of the k th cycle. */
	/*                                      = 4 the integration procedure */
	/*                                          over the k th cycle does */
	/*                                          not converge (to within the */
	/*                                          required accuracy) due to */
	/*                                          roundoff in the extrapolation */
	/*                                          procedure invoked on this */
	/*                                          cycle. it is assumed that the */
	/*                                          result on this interval is */
	/*                                          the best which can be */
	/*                                          obtained. */
	/*                                      = 5 the integral over the k th */
	/*                                          cycle is probably divergent */
	/*                                          or slowly convergent. it must */
	/*                                          be noted that divergence can */
	/*                                          occur with any other value of */
	/*                                          iwork(k). */
	/*                    if omega = 0 and integr = 1, */
	/*                    the integral is calculated by means of dqagie, */
	/*                    and ier = iwork(1) (with meaning as described */
	/*                    for iwork(k),k = 1). */

	/*         dimensioning parameters */
	/*            limlst - long */
	/*                     limlst gives an upper bound on the number of */
	/*                     cycles, limlst.ge.3. */
	/*                     if limlst.lt.3, the routine will end with ier = 6. */

	/*            lst    - long */
	/*                     on return, lst indicates the number of cycles */
	/*                     actually needed for the integration. */
	/*                     if omega = 0, then lst is set to 1. */

	/*            leniw  - long */
	/*                     dimensioning parameter for iwork. on entry, */
	/*                     (leniw-limlst)/2 equals the maximum number of */
	/*                     subintervals allowed in the partition of each */
	/*                     cycle, leniw.ge.(limlst+2). */
	/*                     if leniw.lt.(limlst+2), the routine will end with */
	/*                     ier = 6. */

	/*            maxp1  - long */
	/*                     maxp1 gives an upper bound on the number of */
	/*                     chebyshev moments which can be stored, i.e. for */
	/*                     the intervals of lengths fabs(b-a)*2**(-l), */
	/*                     l = 0,1, ..., maxp1-2, maxp1.ge.1. */
	/*                     if maxp1.lt.1, the routine will end with ier = 6. */
	/*            lenw   - long */
	/*                     dimensioning parameter for work */
	/*                     lenw must be at least leniw*2+maxp1*25. */
	/*                     if lenw.lt.(leniw*2+maxp1*25), the routine will */
	/*                     end with ier = 6. */

	/*         work arrays */
	/*            iwork  - long */
	/*                     vector of dimension at least leniw */
	/*                     on return, iwork(k) for k = 1, 2, ..., lst */
	/*                     contain the error flags on the cycles. */

	/*            work   - double precision */
	/*                     vector of dimension at least */
	/*                     on return, */
	/*                     work(1), ..., work(lst) contain the integral */
	/*                      approximations over the cycles, */
	/*                     work(limlst+1), ..., work(limlst+lst) contain */
	/*                      the error extimates over the cycles. */
	/*                     further elements of work have no specific */
	/*                     meaning for the user. */

	/* ***references  (none) */
	/* ***routines called  dqawfe,xerror */
	/* ***end prologue  dqawf */




	/*         check validity of limlst, leniw, maxp1 and lenw. */

	/* ***first executable statement  dqawf */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*result = 0.;
	*abserr = 0.;
	if (*limlst < 3 || *leniw < *limlst + 2 || *maxp1 < 1 || *lenw < (*leniw 
																							<< 1) + *maxp1 * 25) {
		goto L10;
	}

	/*         prepare call for dqawfe */

	limit = (*leniw - *limlst) / 2;
	l1 = *limlst + 1;
	l2 = *limlst + l1;
	l3 = limit + l2;
	l4 = limit + l3;
	l5 = limit + l4;
	l6 = limit + l5;
	ll2 = limit + l1;
	dqawfe_(f, a, omega, integr, epsabs, limlst, &limit, maxp1, result, 
			  abserr, neval, ier, &work[1], &work[l1], &iwork[1], lst, &work[l2]
			  , &work[l3], &work[l4], &work[l5], &iwork[l1], &iwork[ll2], &work[
				  l6]);

	/*         call error handler if necessary */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from dqawf", &c__26, ier, &lvl, 26);
	}
	return;
} /* dqawf_ */

void dqawoe_(const D_fp& f, const double *a, const double *b, const double 
				 *omega, const long *integr, const double *epsabs, 
				 const double *epsrel, 
				 const long *limit, const long *icall, const long *maxp1, 
				 double *result, 
				 double *abserr, long *neval, long *ier, long *last, 
				 double *alist__, double *blist, double *rlist, double 
				 *elist, long *iord, long *nnlog, long *momcom, double *
				 chebmo)
{
	/* System generated locals */
	int chebmo_dim1, chebmo_offset, i__1, i__2;
	double d__1, d__2;

	/* Local variables */
	long k;
	double a1, a2, b1, b2;
	long id, nev;
	double area, dres;
	long ksgn, nres;
	double area1, area2, area12;
	double small, erro12, width, defab1, defab2;
	long ierro, ktmin;
	double oflow;
	long nrmax, nrmom;
	double uflow;
	bool noext;
	long iroff1, iroff2, iroff3;
	double res3la[3], error1, error2, rlist2[52];
	long numrl2;
	double defabs, domega, epmach, erlarg=0., abseps, correc=0., errbnd, 
		resabs;
	long jupbnd;
	bool extall;
	double erlast, errmax;
	long maxerr;
	double reseps;
	bool extrap;
	double ertest=0., errsum;

	/* ***begin prologue  dqawoe */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             integrand with oscillatory cos or sin factor, */
	/*             clenshaw-curtis method, (end point) singularities, */
	/*             extrapolation, globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral */
	/*            i = integral of f(x)*w(x) over (a,b) */
	/*            where w(x) = cos(omega*x) or w(x)=sin(omega*x), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of oscillatory integrals */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration */

	/*            omega  - double precision */
	/*                     parameter in the integrand weight function */

	/*            integr - long */
	/*                     indicates which of the weight functions is to be */
	/*                     used */
	/*                     integr = 1      w(x) = cos(omega*x) */
	/*                     integr = 2      w(x) = sin(omega*x) */
	/*                     if integr.ne.1 and integr.ne.2, the routine */
	/*                     will end with ier = 6. */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subdivisions */
	/*                     in the partition of (a,b), limit.ge.1. */

	/*            icall  - long */
	/*                     if dqawoe is to be used only once, icall must */
	/*                     be set to 1.  assume that during this call, the */
	/*                     chebyshev moments (for clenshaw-curtis integration */
	/*                     of degree 24) have been computed for intervals of */
	/*                     lenghts (fabs(b-a))*2**(-l), l=0,1,2,...momcom-1. */
	/*                     if icall.gt.1 this means that dqawoe has been */
	/*                     called twice or more on intervals of the same */
	/*                     length fabs(b-a). the chebyshev moments already */
	/*                     computed are then re-used in subsequent calls. */
	/*                     if icall.lt.1, the routine will end with ier = 6. */

	/*            maxp1  - long */
	/*                     gives an upper bound on the number of chebyshev */
	/*                     moments which can be stored, i.e. for the */
	/*                     intervals of lenghts fabs(b-a)*2**(-l), */
	/*                     l=0,1, ..., maxp1-2, maxp1.ge.1. */
	/*                     if maxp1.lt.1, the routine will end with ier = 6. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the */
	/*                             requested accuracy has been achieved. */
	/*                   - ier.gt.0 abnormal termination of the routine. */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand, in order to */
	/*                             determine the integration difficulties. */
	/*                             if the position of a local difficulty can */
	/*                             be determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. */
	/*                             it is presumed that the requested */
	/*                             tolerance cannot be achieved due to */
	/*                             roundoff in the extrapolation table, */
	/*                             and that the returned result is the */
	/*                             best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier.gt.0. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or (integr.ne.1 and integr.ne.2) or */
	/*                             icall.lt.1 or maxp1.lt.1. */
	/*                             result, abserr, neval, last, rlist(1), */
	/*                             elist(1), iord(1) and nnlog(1) are set */
	/*                             to zero. alist(1) and blist(1) are set */
	/*                             to a and b respectively. */

	/*            last  -  long */
	/*                     on return, last equals the number of */
	/*                     subintervals produces in the subdivision */
	/*                     process, which determines the number of */
	/*                     significant elements actually in the */
	/*                     work arrays. */
	/*            alist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the left */
	/*                     end points of the subintervals in the partition */
	/*                     of the given integration range (a,b) */

	/*            blist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the right */
	/*                     end points of the subintervals in the partition */
	/*                     of the given integration range (a,b) */

	/*            rlist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the integral */
	/*                     approximations on the subintervals */

	/*            elist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the moduli of the */
	/*                     absolute error estimates on the subintervals */

	/*            iord   - long */
	/*                     vector of dimension at least limit, the first k */
	/*                     elements of which are pointers to the error */
	/*                     estimates over the subintervals, */
	/*                     such that elist(iord(1)), ..., */
	/*                     elist(iord(k)) form a decreasing sequence, with */
	/*                     k = last if last.le.(limit/2+2), and */
	/*                     k = limit+1-last otherwise. */

	/*            nnlog  - long */
	/*                     vector of dimension at least limit, containing the */
	/*                     subdivision levels of the subintervals, i.e. */
	/*                     iwork(i) = l means that the subinterval */
	/*                     numbered i is of length fabs(b-a)*2**(1-l) */

	/*         on entry and return */
	/*            momcom - long */
	/*                     indicating that the chebyshev moments */
	/*                     have been computed for intervals of lengths */
	/*                     (fabs(b-a))*2**(-l), l=0,1,2, ..., momcom-1, */
	/*                     momcom.lt.maxp1 */

	/*            chebmo - double precision */
	/*                     array of dimension (maxp1,25) containing the */
	/*                     chebyshev moments */

	/* ***references  (none) */
	/* ***routines called  d1mach,dqc25f,dqelg,dqpsrt */
	/* ***end prologue  dqawoe */




	/*            the dimension of rlist2 is determined by  the value of */
	/*            limexp in subroutine dqelg (rlist2 should be of */
	/*            dimension (limexp+2) at least). */

	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           rlist2    - array of dimension at least limexp+2 */
	/*                       containing the part of the epsilon table */
	/*                       which is still needed for further computations */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest */
	/*                       error estimate */
	/*           errmax    - elist(maxerr) */
	/*           erlast    - error on the interval currently subdivided */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */
	/*           nres      - number of calls to the extrapolation routine */
	/*           numrl2    - number of elements in rlist2. if an appropriate */
	/*                       approximation to the compounded integral has */
	/*                       been obtained it is put in rlist2(numrl2) after */
	/*                       numrl2 has been increased by one */
	/*           small     - length of the smallest interval considered */
	/*                       up to now, multiplied by 1.5 */
	/*           erlarg    - sum of the errors over the intervals larger */
	/*                       than the smallest interval considered up to now */
	/*           extrap    - bool variable denoting that the routine is */
	/*                       attempting to perform extrapolation, i.e. before */
	/*                       subdividing the smallest interval we try to */
	/*                       decrease the value of erlarg */
	/*           noext     - bool variable denoting that extrapolation */
	/*                       is no longer allowed (true  value) */

	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */
	/*           oflow is the largest positive magnitude. */

	/* ***first executable statement  dqawoe */
	/* Parameter adjustments */
	--nnlog;
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;
	chebmo_dim1 = *maxp1;
	chebmo_offset = 1 + chebmo_dim1;
	chebmo -= chebmo_offset;

	/* Function Body */
	epmach = d1mach(c__4);

	/*         test on validity of parameters */
	/*         ------------------------------ */

	*ier = 0;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	alist__[1] = *a;
	blist[1] = *b;
	rlist[1] = 0.;
	elist[1] = 0.;
	iord[1] = 0;
	nnlog[1] = 0;
	/* Computing MAX */
	d__1 = epmach * 50.;
	if ( ( *integr != 1 && *integr != 2 ) || 
		 ( *epsabs <= 0. && *epsrel < max(d__1, 5e-29) ) || 
		  *icall < 1 || *maxp1 < 1) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}

	/*           first approximation to the integral */
	/*           ----------------------------------- */

	domega = fabs(*omega);
	nrmom = 0;
	if (*icall > 1) {
		goto L5;
	}
	*momcom = 0;
L5:
	dqc25f_(f, a, b, &domega, integr, &nrmom, maxp1, &c__0, result, 
			  abserr, neval, &defabs, &resabs, momcom, &chebmo[chebmo_offset]);

	/*           test on accuracy. */

	dres = fabs(*result);
	/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * dres;
	errbnd = max(d__1,d__2);
	rlist[1] = *result;
	elist[1] = *abserr;
	iord[1] = 1;
	if (*abserr <= epmach * 100. * defabs && *abserr > errbnd) {
		*ier = 2;
	}
	if (*limit == 1) {
		*ier = 1;
	}
	if (*ier != 0 || *abserr <= errbnd) {
		goto L200;
	}

	/*           initializations */
	/*           --------------- */

	uflow = d1mach(c__1);
	oflow = d1mach(c__2);
	errmax = *abserr;
	maxerr = 1;
	area = *result;
	errsum = *abserr;
	*abserr = oflow;
	nrmax = 1;
	extrap = false;
	noext = false;
	ierro = 0;
	iroff1 = 0;
	iroff2 = 0;
	iroff3 = 0;
	ktmin = 0;
	small = (d__1 = *b - *a, fabs(d__1)) * .75;
	nres = 0;
	numrl2 = 0;
	extall = false;
	if ((d__1 = *b - *a, fabs(d__1)) * .5 * domega > 2.) {
		goto L10;
	}
	numrl2 = 1;
	extall = true;
	rlist2[0] = *result;
L10:
	if ((d__1 = *b - *a, fabs(d__1)) * .25 * domega <= 2.) {
		extall = true;
	}
	ksgn = -1;
	if (dres >= (1. - epmach * 50.) * defabs) {
		ksgn = 1;
	}

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limit;
	for (*last = 2; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with the nrmax-th largest */
		/*           error estimate. */

		nrmom = nnlog[maxerr] + 1;
		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5;
		a2 = b1;
		b2 = blist[maxerr];
		erlast = errmax;
		dqc25f_(f, &a1, &b1, &domega, integr, &nrmom, maxp1, &c__0, &
				  area1, &error1, &nev, &resabs, &defab1, momcom, &chebmo[
					  chebmo_offset]);
		*neval += nev;
		dqc25f_(f, &a2, &b2, &domega, integr, &nrmom, maxp1, &c__1, &
				  area2, &error2, &nev, &resabs, &defab2, momcom, &chebmo[
					  chebmo_offset]);
		*neval += nev;

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (defab1 == error1 || defab2 == error2) {
			goto L25;
		}
		if ((d__1 = rlist[maxerr] - area12, fabs(d__1)) > fabs(area12) * 1e-5 ||
			 erro12 < errmax * .99) {
			goto L20;
		}
		if (extrap) {
			++iroff2;
		}
		if (! extrap) {
			++iroff1;
		}
	L20:
		if (*last > 10 && erro12 > errmax) {
			++iroff3;
		}
	L25:
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		nnlog[maxerr] = nrmom;
		nnlog[*last] = nrmom;
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(area);
		errbnd = max(d__1,d__2);

		/*           test for roundoff error and eventually set error flag. */

		if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
			*ier = 2;
		}
		if (iroff2 >= 5) {
			ierro = 3;
		}

		/*           set error flag in the case that the number of */
		/*           subintervals equals limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at a point of the integration range. */

		/* Computing MAX */
		d__1 = fabs(a1), d__2 = fabs(b2);
		if (max(d__1,d__2) <= (epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) 
		{
			*ier = 4;
		}

		/*           append the newly-created intervals to the list. */

		if (error2 > error1) {
			goto L30;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L40;
	L30:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine dqpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the subinterval */
		/*           with nrmax-th largest error estimate (to bisected next). */

	L40:
		dqpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (errsum <= errbnd) {
			goto L170;
		}
		if (*ier != 0) {
			goto L150;
		}
		if (*last == 2 && extall) {
			goto L120;
		}
		if (noext) {
			goto L140;
		}
		if (! extall) {
			goto L50;
		}
		erlarg -= erlast;
		if ((d__1 = b1 - a1, fabs(d__1)) > small) {
			erlarg += erro12;
		}
		if (extrap) {
			goto L70;
		}

		/*           test whether the interval to be bisected next is the */
		/*           smallest interval. */

	L50:
		width = (d__1 = blist[maxerr] - alist__[maxerr], fabs(d__1));
		if (width > small) {
			goto L140;
		}
		if (extall) {
			goto L60;
		}

		/*           test whether we can start with the extrapolation procedure */
		/*           (we do this if we integrate over the next interval with */
		/*           use of a gauss-kronrod rule - see subroutine dqc25f). */

		small *= .5;
		if (width * .25 * domega > 2.) {
			goto L140;
		}
		extall = true;
		goto L130;
	L60:
		extrap = true;
		nrmax = 2;
	L70:
		if (ierro == 3 || erlarg <= ertest) {
			goto L90;
		}

		/*           the smallest interval has the largest error. */
		/*           before bisecting decrease the sum of the errors over */
		/*           the larger intervals (erlarg) and perform extrapolation. */

		jupbnd = *last;
		if (*last > *limit / 2 + 2) {
			jupbnd = *limit + 3 - *last;
		}
		id = nrmax;
		i__2 = jupbnd;
		for (k = id; k <= i__2; ++k) {
			maxerr = iord[nrmax];
			errmax = elist[maxerr];
			if ((d__1 = blist[maxerr] - alist__[maxerr], fabs(d__1)) > small) {
				goto L140;
			}
			++nrmax;
			/* L80: */
		}

		/*           perform extrapolation. */

	L90:
		++numrl2;
		rlist2[numrl2 - 1] = area;
		if (numrl2 < 3) {
			goto L110;
		}
		dqelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
		++ktmin;
		if (ktmin > 5 && *abserr < errsum * .001) {
			*ier = 5;
		}
		if (abseps >= *abserr) {
			goto L100;
		}
		ktmin = 0;
		*abserr = abseps;
		*result = reseps;
		correc = erlarg;
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(reseps);
		ertest = max(d__1,d__2);
		/* ***jump out of do-loop */
		if (*abserr <= ertest) {
			goto L150;
		}

		/*           prepare bisection of the smallest interval. */

	L100:
		if (numrl2 == 1) {
			noext = true;
		}
		if (*ier == 5) {
			goto L150;
		}
	L110:
		maxerr = iord[1];
		errmax = elist[maxerr];
		nrmax = 1;
		extrap = false;
		small *= .5;
		erlarg = errsum;
		goto L140;
	L120:
		small *= .5;
		++numrl2;
		rlist2[numrl2 - 1] = area;
	L130:
		ertest = errbnd;
		erlarg = errsum;
	L140:
		;
	}

	/*           set the final result. */
	/*           --------------------- */

L150:
	if (*abserr == oflow || nres == 0) {
		goto L170;
	}
	if (*ier + ierro == 0) {
		goto L165;
	}
	if (ierro == 3) {
		*abserr += correc;
	}
	if (*ier == 0) {
		*ier = 3;
	}
	if (*result != 0. && area != 0.) {
		goto L160;
	}
	if (*abserr > errsum) {
		goto L170;
	}
	if (area == 0.) {
		goto L190;
	}
	goto L165;
L160:
	if (*abserr / fabs(*result) > errsum / fabs(area)) {
		goto L170;
	}

	/*           test on divergence. */

L165:
	/* Computing MAX */
	d__1 = fabs(*result), d__2 = fabs(area);
	if (ksgn == -1 && max(d__1,d__2) <= defabs * .01) {
		goto L190;
	}
	if (.01 > *result / area || *result / area > 100. || errsum >= fabs(area)) 
	{
		*ier = 6;
	}
	goto L190;

	/*           compute global integral sum. */

L170:
	*result = 0.;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L180: */
	}
	*abserr = errsum;
L190:
	if (*ier > 2) {
		--(*ier);
	}
L200:
	if (*integr == 2 && *omega < 0.) {
		*result = -(*result);
	}
L999:
	return;
} /* dqawoe_ */

void dqawo_(const D_fp& f, const double *a, const double *b, const double *omega, 
				const long *integr, const double *epsabs, const double *epsrel, 
				double *result, double *abserr, long *neval, long *ier, 
				const long *leniw, long *maxp1, const long *lenw, long *last, 
				long *iwork, double *work)
{
	long l1, l2, l3, l4, lvl, limit;
	long momcom;

	/* ***begin prologue  dqawo */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             integrand with oscillatory cos or sin factor, */
	/*             clenshaw-curtis method, (end point) singularities, */
	/*             extrapolation, globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i=integral of f(x)*w(x) over (a,b) */
	/*            where w(x) = cos(omega*x) */
	/*            or w(x) = sin(omega*x), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of oscillatory integrals */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the function */
	/*                     f(x).  the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration */

	/*            omega  - double precision */
	/*                     parameter in the integrand weight function */

	/*            integr - long */
	/*                     indicates which of the weight functions is used */
	/*                     integr = 1      w(x) = cos(omega*x) */
	/*                     integr = 2      w(x) = sin(omega*x) */
	/*                     if integr.ne.1.and.integr.ne.2, the routine will */
	/*                     end with ier = 6. */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if epsabs.le.0 and */
	/*                     epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of  integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                   - ier.gt.0 abnormal termination of the routine. */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             (= leniw/2) has been achieved. one can */
	/*                             allow more subdivisions by increasing the */
	/*                             value of leniw (and taking the according */
	/*                             dimension adjustments into account). */
	/*                             however, if this yields no improvement it */
	/*                             is advised to analyze the integrand in */
	/*                             order to determine the integration */
	/*                             difficulties. if the position of a local */
	/*                             difficulty can be determined (e.g. */
	/*                             singularity, discontinuity within the */
	/*                             interval) one will probably gain from */
	/*                             splitting up the interval at this polong */
	/*                             and calling the integrator on the */
	/*                             subranges. if possible, an appropriate */
	/*                             special-purpose integrator should be used */
	/*                             which is designed for handling the type of */
	/*                             difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some interior points of the */
	/*                             integration interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. it is presumed that */
	/*                             the requested tolerance cannot be achieved */
	/*                             due to roundoff in the extrapolation */
	/*                             table, and that the returned result is */
	/*                             the best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or (integr.ne.1 and integr.ne.2), */
	/*                             or leniw.lt.2 or maxp1.lt.1 or */
	/*                             lenw.lt.leniw*2+maxp1*25. */
	/*                             result, abserr, neval, last are set to */
	/*                             zero. except when leniw, maxp1 or lenw are */
	/*                             invalid, work(limit*2+1), work(limit*3+1), */
	/*                             iwork(1), iwork(limit+1) are set to zero, */
	/*                             work(1) is set to a and work(limit+1) to */
	/*                             b. */

	/*         dimensioning parameters */
	/*            leniw  - long */
	/*                     dimensioning parameter for iwork. */
	/*                     leniw/2 equals the maximum number of subintervals */
	/*                     allowed in the partition of the given integration */
	/*                     interval (a,b), leniw.ge.2. */
	/*                     if leniw.lt.2, the routine will end with ier = 6. */

	/*            maxp1  - long */
	/*                     gives an upper bound on the number of chebyshev */
	/*                     moments which can be stored, i.e. for the */
	/*                     intervals of lengths fabs(b-a)*2**(-l), */
	/*                     l=0,1, ..., maxp1-2, maxp1.ge.1 */
	/*                     if maxp1.lt.1, the routine will end with ier = 6. */

	/*            lenw   - long */
	/*                     dimensioning parameter for work */
	/*                     lenw must be at least leniw*2+maxp1*25. */
	/*                     if lenw.lt.(leniw*2+maxp1*25), the routine will */
	/*                     end with ier = 6. */

	/*            last   - long */
	/*                     on return, last equals the number of subintervals */
	/*                     produced in the subdivision process, which */
	/*                     determines the number of significant elements */
	/*                     actually in the work arrays. */

	/*         work arrays */
	/*            iwork  - long */
	/*                     vector of dimension at least leniw */
	/*                     on return, the first k elements of which contain */
	/*                     pointers to the error estimates over the */
	/*                     subintervals, such that work(limit*3+iwork(1)), .. */
	/*                     work(limit*3+iwork(k)) form a decreasing */
	/*                     sequence, with limit = lenw/2 , and k = last */
	/*                     if last.le.(limit/2+2), and k = limit+1-last */
	/*                     otherwise. */
	/*                     furthermore, iwork(limit+1), ..., iwork(limit+ */
	/*                     last) indicate the subdivision levels of the */
	/*                     subintervals, such that iwork(limit+i) = l means */
	/*                     that the subinterval numbered i is of length */
	/*                     fabs(b-a)*2**(1-l). */

	/*            work   - double precision */
	/*                     vector of dimension at least lenw */
	/*                     on return */
	/*                     work(1), ..., work(last) contain the left */
	/*                      end points of the subintervals in the */
	/*                      partition of (a,b), */
	/*                     work(limit+1), ..., work(limit+last) contain */
	/*                      the right end points, */
	/*                     work(limit*2+1), ..., work(limit*2+last) contain */
	/*                      the integral approximations over the */
	/*                      subintervals, */
	/*                     work(limit*3+1), ..., work(limit*3+last) */
	/*                      contain the error estimates. */
	/*                     work(limit*4+1), ..., work(limit*4+maxp1*25) */
	/*                      provide space for storing the chebyshev moments. */
	/*                     note that limit = lenw/2. */

	/* ***references  (none) */
	/* ***routines called  dqawoe,xerror */
	/* ***end prologue  dqawo */




	/*         check validity of leniw, maxp1 and lenw. */

	/* ***first executable statement  dqawo */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	if (*leniw < 2 || *maxp1 < 1 || *lenw < (*leniw << 1) + *maxp1 * 25) {
		goto L10;
	}

	/*         prepare call for dqawoe */

	limit = *leniw / 2;
	l1 = limit + 1;
	l2 = limit + l1;
	l3 = limit + l2;
	l4 = limit + l3;
	dqawoe_(f, a, b, omega, integr, epsabs, epsrel, &limit, &c__1, 
			  maxp1, result, abserr, neval, ier, last, &work[1], &work[l1], &
			  work[l2], &work[l3], &iwork[1], &iwork[l1], &momcom, &work[l4]);

	/*         call error handler if necessary */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 0;
	}
	if (*ier != 0) {
		xerror_("abnormal return from dqawo", &c__26, ier, &lvl, 26);
	}
	return;
} /* dqawo_ */

void dqawse_(const D_fp& f, const double *a, const double *b, const double *alfa, 
				 const double *beta, const long *integr, const double *epsabs, 
				 const double *epsrel, const long *limit, double *result, 
				 double *abserr, long *neval, long *ier, double *alist__, 
				 double *blist, double *rlist, double *elist, long *iord, 
				 long *last)
{
	/* System generated locals */
	int i__1;
	double d__1, d__2;

	/* Local variables */
	long k;
	double a1, a2, b1, b2, rg[25], rh[25], ri[25], rj[25];
	long nev;
	double area, area1, area2, area12;
	double erro12;
	long nrmax;
	double uflow;
	long iroff1, iroff2;
	double resas1, resas2, error1, error2, epmach, errbnd, centre;
	double errmax;
	long maxerr;
	double errsum;

	/* ***begin prologue  dqawse */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             algebraico-logarithmic end point singularities, */
	/*             clenshaw-curtis method */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f*w over (a,b), */
	/*            (where w shows a singular behaviour at the end points, */
	/*            see parameter integr). */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        integration of functions having algebraico-logarithmic */
	/*        end point singularities */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration, b.gt.a */
	/*                     if b.le.a, the routine will end with ier = 6. */

	/*            alfa   - double precision */
	/*                     parameter in the weight function, alfa.gt.(-1) */
	/*                     if alfa.le.(-1), the routine will end with */
	/*                     ier = 6. */

	/*            beta   - double precision */
	/*                     parameter in the weight function, beta.gt.(-1) */
	/*                     if beta.le.(-1), the routine will end with */
	/*                     ier = 6. */

	/*            integr - long */
	/*                     indicates which weight function is to be used */
	/*                     = 1  (x-a)**alfa*(b-x)**beta */
	/*                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a) */
	/*                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x) */
	/*                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x) */
	/*                     if integr.lt.1 or integr.gt.4, the routine */
	/*                     will end with ier = 6. */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subintervals */
	/*                     in the partition of (a,b), limit.ge.2 */
	/*                     if limit.lt.2, the routine will end with ier = 6. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for the integral and error */
	/*                             are less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                         = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit. however, if this yields no */
	/*                             improvement, it is advised to analyze the */
	/*                             integrand in order to determine the */
	/*                             integration difficulties which prevent the */
	/*                             requested tolerance from being achieved. */
	/*                             in case of a jump discontinuity or a local */
	/*                             singularity of algebraico-logarithmic type */
	/*                             at one or more interior points of the */
	/*                             integration range, one should proceed by */
	/*                             splitting up the interval at these */
	/*                             points and calling the integrator on the */
	/*                             subranges. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 6 the input is invalid, because */
	/*                             b.le.a or alfa.le.(-1) or beta.le.(-1), or */
	/*                             integr.lt.1 or integr.gt.4, or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                             or limit.lt.2. */
	/*                             result, abserr, neval, rlist(1), elist(1), */
	/*                             iord(1) and last are set to zero. alist(1) */
	/*                             and blist(1) are set to a and b */
	/*                             respectively. */

	/*            alist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the left */
	/*                     end points of the subintervals in the partition */
	/*                     of the given integration range (a,b) */

	/*            blist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the right */
	/*                     end points of the subintervals in the partition */
	/*                     of the given integration range (a,b) */

	/*            rlist  - double precision */
	/*                     vector of dimension at least limit,the first */
	/*                      last  elements of which are the integral */
	/*                     approximations on the subintervals */

	/*            elist  - double precision */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the moduli of the */
	/*                     absolute error estimates on the subintervals */

	/*            iord   - long */
	/*                     vector of dimension at least limit, the first k */
	/*                     of which are pointers to the error */
	/*                     estimates over the subintervals, so that */
	/*                     elist(iord(1)), ..., elist(iord(k)) with k = last */
	/*                     if last.le.(limit/2+2), and k = limit+1-last */
	/*                     otherwise form a decreasing sequence */

	/*            last   - long */
	/*                     number of subintervals actually produced in */
	/*                     the subdivision process */

	/* ***references  (none) */
	/* ***routines called  d1mach,dqc25s,dqmomo,dqpsrt */
	/* ***end prologue  dqawse */




	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest */
	/*                       error estimate */
	/*           errmax    - elist(maxerr) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */


	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  dqawse */
	/* Parameter adjustments */
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);

	/*           test on validity of parameters */
	/*           ------------------------------ */

	*ier = 6;
	*neval = 0;
	*last = 0;
	rlist[1] = 0.;
	elist[1] = 0.;
	iord[1] = 0;
	*result = 0.;
	*abserr = 0.;
	/* Computing MAX */
	d__1 = epmach * 50.;
	if (*b <= *a || ( *epsabs == 0. && *epsrel < max(d__1,5e-29) ) || *alfa <= 
	    -1. || *beta <= -1. || *integr < 1 || *integr > 4 || *limit < 2) {
		goto L999;
	}
	*ier = 0;

	/*           compute the modified chebyshev moments. */

	dqmomo_(alfa, beta, ri, rj, rg, rh, integr);

	/*           integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b). */

	centre = (*b + *a) * .5;
	dqc25s_(f, a, b, a, &centre, alfa, beta, ri, rj, rg, rh, &area1, &
			  error1, &resas1, integr, &nev);
	*neval = nev;
	dqc25s_(f, a, b, &centre, b, alfa, beta, ri, rj, rg, rh, &area2, &
			  error2, &resas2, integr, &nev);
	*last = 2;
	*neval += nev;
	*result = area1 + area2;
	*abserr = error1 + error2;

	/*           test on accuracy. */

	/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * fabs(*result);
	errbnd = max(d__1,d__2);

	/*           initialization */
	/*           -------------- */

	if (error2 > error1) {
		goto L10;
	}
	alist__[1] = *a;
	alist__[2] = centre;
	blist[1] = centre;
	blist[2] = *b;
	rlist[1] = area1;
	rlist[2] = area2;
	elist[1] = error1;
	elist[2] = error2;
	goto L20;
L10:
	alist__[1] = centre;
	alist__[2] = *a;
	blist[1] = *b;
	blist[2] = centre;
	rlist[1] = area2;
	rlist[2] = area1;
	elist[1] = error2;
	elist[2] = error1;
L20:
	iord[1] = 1;
	iord[2] = 2;
	if (*limit == 2) {
		*ier = 1;
	}
	if (*abserr <= errbnd || *ier == 1) {
		goto L999;
	}
	errmax = elist[1];
	maxerr = 1;
	nrmax = 1;
	area = *result;
	errsum = *abserr;
	iroff1 = 0;
	iroff2 = 0;

	/*            main do-loop */
	/*            ------------ */

	i__1 = *limit;
	for (*last = 3; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with largest error estimate. */

		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5;
		a2 = b1;
		b2 = blist[maxerr];

		dqc25s_(f, a, b, &a1, &b1, alfa, beta, ri, rj, rg, rh, &area1, &
				  error1, &resas1, integr, &nev);
		*neval += nev;
		dqc25s_(f, a, b, &a2, &b2, alfa, beta, ri, rj, rg, rh, &area2, &
				  error2, &resas2, integr, &nev);
		*neval += nev;

		/*           improve previous approximations integral and error */
		/*           and test for accuracy. */

		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (*a == a1 || *b == b2) {
			goto L30;
		}
		if (resas1 == error1 || resas2 == error2) {
			goto L30;
		}

		/*           test for roundoff error. */

		if ((d__1 = rlist[maxerr] - area12, fabs(d__1)) < fabs(area12) * 1e-5 &&
			 erro12 >= errmax * .99) {
			++iroff1;
		}
		if (*last > 10 && erro12 > errmax) {
			++iroff2;
		}
	L30:
		rlist[maxerr] = area1;
		rlist[*last] = area2;

		/*           test on accuracy. */

		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(area);
		errbnd = max(d__1,d__2);
		if (errsum <= errbnd) {
			goto L35;
		}

		/*           set error flag in the case that the number of interval */
		/*           bisections exceeds limit. */

		if (*last == *limit) {
			*ier = 1;
		}


		/*           set error flag in the case of roundoff error. */

		if (iroff1 >= 6 || iroff2 >= 20) {
			*ier = 2;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at interior points of integration range. */

		/* Computing MAX */
		d__1 = fabs(a1), d__2 = fabs(b2);
		if (max(d__1,d__2) <= (epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) 
		{
			*ier = 3;
		}

		/*           append the newly-created intervals to the list. */

	L35:
		if (error2 > error1) {
			goto L40;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L50;
	L40:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine dqpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the subinterval */
		/*           with largest error estimate (to be bisected next). */

	L50:
		dqpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (*ier != 0 || errsum <= errbnd) {
			goto L70;
		}
		/* L60: */
	}

	/*           compute final result. */
	/*           --------------------- */

L70:
	*result = 0.;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L80: */
	}
	*abserr = errsum;
L999:
	return;
} /* dqawse_ */

void dqaws_(const D_fp& f, const double *a, const double *b, const double *alfa, 
				const double *beta, const long *integr, const double *epsabs, 
				const double *epsrel, double *result, double *abserr, 
				long *neval, long *ier, long *limit, const long *lenw, long *last, 
				long *iwork, double *work)
{
	long l1, l2, l3, lvl;

	/* ***begin prologue  dqaws */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             algebraico-logarithmic end-point singularities, */
	/*             clenshaw-curtis, globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. -k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f*w over (a,b), */
	/*            (where w shows a singular behaviour at the end points */
	/*            see parameter integr). */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        integration of functions having algebraico-logarithmic */
	/*        end point singularities */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*            f      - double precision */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - double precision */
	/*                     lower limit of integration */

	/*            b      - double precision */
	/*                     upper limit of integration, b.gt.a */
	/*                     if b.le.a, the routine will end with ier = 6. */

	/*            alfa   - double precision */
	/*                     parameter in the integrand function, alfa.gt.(-1) */
	/*                     if alfa.le.(-1), the routine will end with */
	/*                     ier = 6. */

	/*            beta   - double precision */
	/*                     parameter in the integrand function, beta.gt.(-1) */
	/*                     if beta.le.(-1), the routine will end with */
	/*                     ier = 6. */

	/*            integr - long */
	/*                     indicates which weight function is to be used */
	/*                     = 1  (x-a)**alfa*(b-x)**beta */
	/*                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a) */
	/*                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x) */
	/*                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x) */
	/*                     if integr.lt.1 or integr.gt.4, the routine */
	/*                     will end with ier = 6. */

	/*            epsabs - double precision */
	/*                     absolute accuracy requested */
	/*            epsrel - double precision */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*         on return */
	/*            result - double precision */
	/*                     approximation to the integral */

	/*            abserr - double precision */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for the integral and error */
	/*                             are less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand, in order to */
	/*                             determine the integration difficulties */
	/*                             which prevent the requested tolerance from */
	/*                             being achieved. in case of a jump */
	/*                             discontinuity or a local singularity */
	/*                             of algebraico-logarithmic type at one or */
	/*                             more interior points of the integration */
	/*                             range, one should proceed by splitting up */
	/*                             the interval at these points and calling */
	/*                             the integrator on the subranges. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 6 the input is invalid, because */
	/*                             b.le.a or alfa.le.(-1) or beta.le.(-1) or */
	/*                             or integr.lt.1 or integr.gt.4 or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or limit.lt.2 or lenw.lt.limit*4. */
	/*                             result, abserr, neval, last are set to */
	/*                             zero. except when lenw or limit is invalid */
	/*                             iwork(1), work(limit*2+1) and */
	/*                             work(limit*3+1) are set to zero, work(1) */
	/*                             is set to a and work(limit+1) to b. */

	/*         dimensioning parameters */
	/*            limit  - long */
	/*                     dimensioning parameter for iwork */
	/*                     limit determines the maximum number of */
	/*                     subintervals in the partition of the given */
	/*                     integration interval (a,b), limit.ge.2. */
	/*                     if limit.lt.2, the routine will end with ier = 6. */

	/*            lenw   - long */
	/*                     dimensioning parameter for work */
	/*                     lenw must be at least limit*4. */
	/*                     if lenw.lt.limit*4, the routine will end */
	/*                     with ier = 6. */

	/*            last   - long */
	/*                     on return, last equals the number of */
	/*                     subintervals produced in the subdivision process, */
	/*                     which determines the significant number of */
	/*                     elements actually in the work arrays. */

	/*         work arrays */
	/*            iwork  - long */
	/*                     vector of dimension limit, the first k */
	/*                     elements of which contain pointers */
	/*                     to the error estimates over the subintervals, */
	/*                     such that work(limit*3+iwork(1)), ..., */
	/*                     work(limit*3+iwork(k)) form a decreasing */
	/*                     sequence with k = last if last.le.(limit/2+2), */
	/*                     and k = limit+1-last otherwise */

	/*            work   - double precision */
	/*                     vector of dimension lenw */
	/*                     on return */
	/*                     work(1), ..., work(last) contain the left */
	/*                      end points of the subintervals in the */
	/*                      partition of (a,b), */
	/*                     work(limit+1), ..., work(limit+last) contain */
	/*                      the right end points, */
	/*                     work(limit*2+1), ..., work(limit*2+last) */
	/*                      contain the integral approximations over */
	/*                      the subintervals, */
	/*                     work(limit*3+1), ..., work(limit*3+last) */
	/*                      contain the error estimates. */

	/* ***references  (none) */
	/* ***routines called  dqawse,xerror */
	/* ***end prologue  dqaws */




	/*         check validity of limit and lenw. */

	/* ***first executable statement  dqaws */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.;
	*abserr = 0.;
	if (*limit < 2 || *lenw < *limit << 2) {
		goto L10;
	}

	/*         prepare call for dqawse. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;

	dqawse_(f, a, b, alfa, beta, integr, epsabs, epsrel, limit, result, 
			  abserr, neval, ier, &work[1], &work[l1], &work[l2], &work[l3], &
			  iwork[1], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from dqaws", &c__26, ier, &lvl, 26);
	}
	return;
} /* dqaws_ */

void dqc25c_(const D_fp& f, const double *a, const double *b, const double *c__, 
				 double *result, double *abserr, long *krul, long *neval)
{
	/* Initialized data */

	double x[11] = { .991444861373810411144557526928563,
									.965925826289068286749743199728897,
									.923879532511286756128183189396788,
									.866025403784438646763723170752936,
									.793353340291235164579776961501299,
									.707106781186547524400844362104849,
									.608761429008720639416097542898164,.5,
									.382683432365089771728459984030399,
									.258819045102520762348898837624048,
									.130526192220051591548406227895489 };

	/* System generated locals */
	double d__1;

	/* Local variables */
	long i__, k;
	double u, p2, p3, p4, cc;
	long kp;
	double ak22, fval[25], res12, res24;
	long isym;
	double amom0, amom1, amom2, cheb12[13], cheb24[25], hlgth, 
		centr;
	double resabs, resasc;

	/* ***begin prologue  dqc25c */
	/* ***date written   810101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a2,j4 */
	/* ***keywords  25-point clenshaw-curtis integration */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f*w over (a,b) with */
	/*            error estimate, where w(x) = 1/(x-c) */
	/* ***description */

	/*        integration rules for the computation of cauchy */
	/*        principal value integrals */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*           f      - double precision */
	/*                    function subprogram defining the integrand function */
	/*                    f(x). the actual name for f needs to be declared */
	/*                    e x t e r n a l  in the driver program. */

	/*           a      - double precision */
	/*                    left end point of the integration interval */

	/*           b      - double precision */
	/*                    right end point of the integration interval, b.gt.a */

	/*           c      - double precision */
	/*                    parameter in the weight function */

	/*           result - double precision */
	/*                    approximation to the integral */
	/*                    result is computed by using a generalized */
	/*                    clenshaw-curtis method if c lies within ten percent */
	/*                    of the integration interval. in the other case the */
	/*                    15-point kronrod rule obtained by optimal addition */
	/*                    of abscissae to the 7-point gauss rule, is applied. */

	/*           abserr - double precision */
	/*                    estimate of the modulus of the absolute error, */
	/*                    which should equal or exceed fabs(i-result) */

	/*           krul   - long */
	/*                    key which is decreased by 1 if the 15-polong */
	/*                    gauss-kronrod scheme has been used */

	/*           neval  - long */
	/*                    number of integrand evaluations */

	/* ....................................................................... */
	/* ***references  (none) */
	/* ***routines called  dqcheb,dqk15w,dqwgtc */
	/* ***end prologue  dqc25c */




	/*           the vector x contains the values cos(k*pi/24), */
	/*           k = 1, ..., 11, to be used for the chebyshev series */
	/*           expansion of f */


	/*           list of major variables */
	/*           ---------------------- */
	/*           fval   - value of the function f at the points */
	/*                    cos(k*pi/24),  k = 0, ..., 24 */
	/*           cheb12 - chebyshev series expansion coefficients, */
	/*                    for the function f, of degree 12 */
	/*           cheb24 - chebyshev series expansion coefficients, */
	/*                    for the function f, of degree 24 */
	/*           res12  - approximation to the integral corresponding */
	/*                    to the use of cheb12 */
	/*           res24  - approximation to the integral corresponding */
	/*                    to the use of cheb24 */
	/*           dqwgtc - external function subprogram defining */
	/*                    the weight function */
	/*           hlgth  - half-length of the interval */
	/*           centr  - mid point of the interval */


	/*           check the position of c. */

	set_NaN( p2 );
	set_NaN( p3 );
	set_NaN( p4 );
	kp = 0;

	/* ***first executable statement  dqc25c */
	cc = (2. * *c__ - *b - *a) / (*b - *a);
	if (fabs(cc) < 1.1) {
		goto L10;
	}

	/*           apply the 15-point gauss-kronrod scheme. */

	--(*krul);
	dqk15w_(f, dqwgtc_, c__, &p2, &p3, &p4, &kp, a, b, result, 
			  abserr, &resabs, &resasc);
	*neval = 15;
	if (resasc == *abserr) {
		++(*krul);
	}
	goto L50;

	/*           use the generalized clenshaw-curtis method. */

L10:
	hlgth = (*b - *a) * .5;
	centr = (*b + *a) * .5;
	*neval = 25;
	d__1 = hlgth + centr;
	fval[0] = f(d__1) * .5;
	fval[12] = f(centr);
	d__1 = centr - hlgth;
	fval[24] = f(d__1) * .5;
	for (i__ = 2; i__ <= 12; ++i__) {
		u = hlgth * x[i__ - 2];
		isym = 26 - i__;
		d__1 = u + centr;
		fval[i__ - 1] = f(d__1);
		d__1 = centr - u;
		fval[isym - 1] = f(d__1);
		/* L20: */
	}

	/*           compute the chebyshev series expansion. */

	dqcheb_(x, fval, cheb12, cheb24);

	/*           the modified chebyshev moments are computed by forward */
	/*           recursion, using amom0 and amom1 as starting values. */

	amom0 = log((d__1 = (1. - cc) / (cc + 1.), fabs(d__1)));
	amom1 = cc * amom0 + 2.;
	res12 = cheb12[0] * amom0 + cheb12[1] * amom1;
	res24 = cheb24[0] * amom0 + cheb24[1] * amom1;
	for (k = 3; k <= 13; ++k) {
		amom2 = cc * 2. * amom1 - amom0;
		ak22 = (double) ((k - 2) * (k - 2));
		if (k / 2 << 1 == k) {
			amom2 -= 4. / (ak22 - 1.);
		}
		res12 += cheb12[k - 1] * amom2;
		res24 += cheb24[k - 1] * amom2;
		amom0 = amom1;
		amom1 = amom2;
		/* L30: */
	}
	for (k = 14; k <= 25; ++k) {
		amom2 = cc * 2. * amom1 - amom0;
		ak22 = (double) ((k - 2) * (k - 2));
		if (k / 2 << 1 == k) {
			amom2 -= 4. / (ak22 - 1.);
		}
		res24 += cheb24[k - 1] * amom2;
		amom0 = amom1;
		amom1 = amom2;
		/* L40: */
	}
	*result = res24;
	*abserr = (d__1 = res24 - res12, fabs(d__1));
L50:
	return;
} /* dqc25c_ */

void dqc25f_(const D_fp& f, const double *a, const double *b, const double 
				 *omega, const long *integr, const long *nrmom, const long *maxp1, 
				 const long *ksave, double *result, double *abserr, long *neval, 
				 double *resabs, double *resasc, long *momcom, double *
				 chebmo)
{
	/* Initialized data */

	double x[11] = { .991444861373810411144557526928563,
									.965925826289068286749743199728897,
									.923879532511286756128183189396788,
									.866025403784438646763723170752936,
									.793353340291235164579776961501299,
									.707106781186547524400844362104849,
									.608761429008720639416097542898164,.5,
									.382683432365089771728459984030399,
									.258819045102520762348898837624048,
									.130526192220051591548406227895489 };

	/* System generated locals */
	int chebmo_dim1, chebmo_offset, i__1;
	double d__1, d__2;

	/* Local variables */
	double d__[25];
	long i__, j, k, m=0;
	double v[28], d1[25], d2[25], p2, p3, p4, ac, an, as, an2, ass,
		par2, conc, asap, par22, fval[25], estc, cons;
	long iers;
	double ests;
	long isym, noeq1;
	double cheb12[13], cheb24[25], resc12, resc24, hlgth, centr;
	double ress12, ress24, oflow;
	long noequ;
	double cospar;
	double parint, sinpar;

	/* ***begin prologue  dqc25f */
	/* ***date written   810101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a2 */
	/* ***keywords  integration rules for functions with cos or sin */
	/*             factor, clenshaw-curtis, gauss-kronrod */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute the integral i=integral of f(x) over (a,b) */
	/*            where w(x) = cos(omega*x) or w(x)=sin(omega*x) and to */
	/*            compute j = integral of fabs(f) over (a,b). for small value */
	/*            of omega or small intervals (a,b) the 15-point gauss-kronro */
	/*            rule is used. otherwise a generalized clenshaw-curtis */
	/*            method is used. */
	/* ***description */

	/*        integration rules for functions with cos or sin factor */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*         on entry */
	/*           f      - double precision */
	/*                    function subprogram defining the integrand */
	/*                    function f(x). the actual name for f needs to */
	/*                    be declared e x t e r n a l in the calling program. */

	/*           a      - double precision */
	/*                    lower limit of integration */

	/*           b      - double precision */
	/*                    upper limit of integration */

	/*           omega  - double precision */
	/*                    parameter in the weight function */

	/*           integr - long */
	/*                    indicates which weight function is to be used */
	/*                       integr = 1   w(x) = cos(omega*x) */
	/*                       integr = 2   w(x) = sin(omega*x) */

	/*           nrmom  - long */
	/*                    the length of interval (a,b) is equal to the length */
	/*                    of the original integration interval divided by */
	/*                    2**nrmom (we suppose that the routine is used in an */
	/*                    adaptive integration process, otherwise set */
	/*                    nrmom = 0). nrmom must be zero at the first call. */

	/*           maxp1  - long */
	/*                    gives an upper bound on the number of chebyshev */
	/*                    moments which can be stored, i.e. for the */
	/*                    intervals of lengths fabs(bb-aa)*2**(-l), */
	/*                    l = 0,1,2, ..., maxp1-2. */

	/*           ksave  - long */
	/*                    key which is one when the moments for the */
	/*                    current interval have been computed */

	/*         on return */
	/*           result - double precision */
	/*                    approximation to the integral i */

	/*           abserr - double precision */
	/*                    estimate of the modulus of the absolute */
	/*                    error, which should equal or exceed fabs(i-result) */

	/*           neval  - long */
	/*                    number of integrand evaluations */

	/*           resabs - double precision */
	/*                    approximation to the integral j */

	/*           resasc - double precision */
	/*                    approximation to the integral of fabs(f-i/(b-a)) */

	/*         on entry and return */
	/*           momcom - long */
	/*                    for each interval length we need to compute the */
	/*                    chebyshev moments. momcom counts the number of */
	/*                    intervals for which these moments have already been */
	/*                    computed. if nrmom.lt.momcom or ksave = 1, the */
	/*                    chebyshev moments for the interval (a,b) have */
	/*                    already been computed and stored, otherwise we */
	/*                    compute them and we increase momcom. */

	/*           chebmo - double precision */
	/*                    array of dimension at least (maxp1,25) containing */
	/*                    the modified chebyshev moments for the first momcom */
	/*                    momcom interval lengths */

	/* ...................................................................... */
	/* ***references  (none) */
	/* ***routines called  d1mach,dgtsl,dqcheb,dqk15w,dqwgtf */
	/* ***end prologue  dqc25f */




	/*           the vector x contains the values cos(k*pi/24) */
	/*           k = 1, ...,11, to be used for the chebyshev expansion of f */

	/* Parameter adjustments */
	chebmo_dim1 = *maxp1;
	chebmo_offset = 1 + chebmo_dim1;
	chebmo -= chebmo_offset;

	set_NaN( p2 );
	set_NaN( p3 );
	set_NaN( p4 );

	/* Function Body */

	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the integration interval */
	/*           hlgth  - half-length of the integration interval */
	/*           fval   - value of the function f at the points */
	/*                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5, k = 0, ..., 24 */
	/*           cheb12 - coefficients of the chebyshev series expansion */
	/*                    of degree 12, for the function f, in the */
	/*                    interval (a,b) */
	/*           cheb24 - coefficients of the chebyshev series expansion */
	/*                    of degree 24, for the function f, in the */
	/*                    interval (a,b) */
	/*           resc12 - approximation to the integral of */
	/*                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a)) */
	/*                    over (-1,+1), using the chebyshev series */
	/*                    expansion of degree 12 */
	/*           resc24 - approximation to the same integral, using the */
	/*                    chebyshev series expansion of degree 24 */
	/*           ress12 - the analogue of resc12 for the sine */
	/*           ress24 - the analogue of resc24 for the sine */


	/*           machine dependent constant */
	/*           -------------------------- */

	/*           oflow is the largest positive magnitude. */

	/* ***first executable statement  dqc25f */
	oflow = d1mach(c__2);

	centr = (*b + *a) * .5;
	hlgth = (*b - *a) * .5;
	parint = *omega * hlgth;

	/*           compute the integral using the 15-point gauss-kronrod */
	/*           formula if the value of the parameter in the integrand */
	/*           is small. */

	if (fabs(parint) > 2.) {
		goto L10;
	}
	dqk15w_(f, dqwgtf_, omega, &p2, &p3, &p4, integr, a, b, 
			  result, abserr, resabs, resasc);
	*neval = 15;
	goto L170;

	/*           compute the integral using the generalized clenshaw- */
	/*           curtis method. */

L10:
	conc = hlgth * cos(centr * *omega);
	cons = hlgth * sin(centr * *omega);
	*resasc = oflow;
	*neval = 25;

	/*           check whether the chebyshev moments for this interval */
	/*           have already been computed. */

	if (*nrmom < *momcom || *ksave == 1) {
		goto L120;
	}

	/*           compute a new set of chebyshev moments. */

	m = *momcom + 1;
	par2 = parint * parint;
	par22 = par2 + 2.;
	sinpar = sin(parint);
	cospar = cos(parint);

	/*           compute the chebyshev moments with respect to cosine. */

	v[0] = sinpar * 2. / parint;
	v[1] = (cospar * 8. + (par2 + par2 - 8.) * sinpar / parint) / par2;
	v[2] = ((par2 - 12.) * 32. * cospar + ((par2 - 80.) * par2 + 192.) * 2. * 
			  sinpar / parint) / (par2 * par2);
	ac = cospar * 8.;
	as = parint * 24. * sinpar;
	if (fabs(parint) > 24.) {
		goto L30;
	}

	/*           compute the chebyshev moments as the solutions of a */
	/*           boundary value problem with 1 initial value (v(3)) and 1 */
	/*           end value (computed using an asymptotic formula). */

	noequ = 25;
	noeq1 = noequ - 1;
	an = 6.;
	i__1 = noeq1;
	for (k = 1; k <= i__1; ++k) {
		an2 = an * an;
		d__[k - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
		d2[k - 1] = (an - 1.) * (an - 2.) * par2;
		d1[k] = (an + 3.) * (an + 4.) * par2;
		v[k + 2] = as - (an2 - 4.) * ac;
		an += 2.;
		/* L20: */
	}
	an2 = an * an;
	d__[noequ - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
	v[noequ + 2] = as - (an2 - 4.) * ac;
	v[3] -= par2 * 56. * v[2];
	ass = parint * sinpar;
	asap = (((((par2 * 210. - 1.) * cospar - (par2 * 105. - 63.) * ass) / an2 
				 - (1. - par2 * 15.) * cospar + ass * 15.) / an2 - cospar + ass * 
				3.) / an2 - cospar) / an2;
	v[noequ + 2] -= asap * 2. * par2 * (an - 1.) * (an - 2.);

	/*           solve the tridiagonal system by means of gaussian */
	/*           elimination with partial pivoting. */

	/* ***        call to dgtsl must be replaced by call to */
	/* ***        double precision version of linpack routine sgtsl */

	dgtsl_(&noequ, d1, d__, d2, &v[3], &iers);
	goto L50;

	/*           compute the chebyshev moments by means of forward */
	/*           recursion. */

L30:
	an = 4.;
	for (i__ = 4; i__ <= 13; ++i__) {
		an2 = an * an;
		v[i__ - 1] = ((an2 - 4.) * ((par22 - an2 - an2) * 2. * v[i__ - 2] - 
											 ac) + as - par2 * (an + 1.) * (an + 2.) * v[i__ - 3]) / (par2 
																														 * (an - 1.) * (an - 2.));
		an += 2.;
		/* L40: */
	}
L50:
	for (j = 1; j <= 13; ++j) {
		chebmo[m + ((j << 1) - 1) * chebmo_dim1] = v[j - 1];
		/* L60: */
	}

	/*           compute the chebyshev moments with respect to sine. */

	v[0] = (sinpar - parint * cospar) * 2. / par2;
	v[1] = (18. - 48. / par2) * sinpar / par2 + (48. / par2 - 2.) * cospar / 
		parint;
	ac = parint * -24. * cospar;
	as = sinpar * -8.;
	if (fabs(parint) > 24.) {
		goto L80;
	}

	/*           compute the chebyshev moments as the solutions of a boundary */
	/*           value problem with 1 initial value (v(2)) and 1 end value */
	/*           (computed using an asymptotic formula). */

	an = 5.;
	i__1 = noeq1;
	for (k = 1; k <= i__1; ++k) {
		an2 = an * an;
		d__[k - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
		d2[k - 1] = (an - 1.) * (an - 2.) * par2;
		d1[k] = (an + 3.) * (an + 4.) * par2;
		v[k + 1] = ac + (an2 - 4.) * as;
		an += 2.;
		/* L70: */
	}
	an2 = an * an;
	d__[noequ - 1] = (an2 - 4.) * -2. * (par22 - an2 - an2);
	v[noequ + 1] = ac + (an2 - 4.) * as;
	v[2] -= par2 * 42. * v[1];
	ass = parint * cospar;
	asap = (((((par2 * 105. - 63.) * ass + (par2 * 210. - 1.) * sinpar) / an2 
				 + (par2 * 15. - 1.) * sinpar - ass * 15.) / an2 - ass * 3. - 
				sinpar) / an2 - sinpar) / an2;
	v[noequ + 1] -= asap * 2. * par2 * (an - 1.) * (an - 2.);

	/*           solve the tridiagonal system by means of gaussian */
	/*           elimination with partial pivoting. */

	/* ***        call to dgtsl must be replaced by call to */
	/* ***        double precision version of linpack routine sgtsl */

	dgtsl_(&noequ, d1, d__, d2, &v[2], &iers);
	goto L100;

	/*           compute the chebyshev moments by means of forward recursion. */

L80:
	an = 3.;
	for (i__ = 3; i__ <= 12; ++i__) {
		an2 = an * an;
		v[i__ - 1] = ((an2 - 4.) * ((par22 - an2 - an2) * 2. * v[i__ - 2] + 
											 as) + ac - par2 * (an + 1.) * (an + 2.) * v[i__ - 3]) / (par2 
																														 * (an - 1.) * (an - 2.));
		an += 2.;
		/* L90: */
	}
L100:
	for (j = 1; j <= 12; ++j) {
		chebmo[m + (j << 1) * chebmo_dim1] = v[j - 1];
		/* L110: */
	}
L120:
	if (*nrmom < *momcom) {
		m = *nrmom + 1;
	}
	if (*momcom < *maxp1 - 1 && *nrmom >= *momcom) {
		++(*momcom);
	}

	/*           compute the coefficients of the chebyshev expansions */
	/*           of degrees 12 and 24 of the function f. */

	d__1 = centr + hlgth;
	fval[0] = f(d__1) * .5;
	fval[12] = f(centr);
	d__1 = centr - hlgth;
	fval[24] = f(d__1) * .5;
	for (i__ = 2; i__ <= 12; ++i__) {
		isym = 26 - i__;
		d__1 = hlgth * x[i__ - 2] + centr;
		fval[i__ - 1] = f(d__1);
		d__1 = centr - hlgth * x[i__ - 2];
		fval[isym - 1] = f(d__1);
		/* L130: */
	}
	dqcheb_(x, fval, cheb12, cheb24);

	/*           compute the integral and error estimates. */

	resc12 = cheb12[12] * chebmo[m + chebmo_dim1 * 13];
	ress12 = 0.;
	k = 11;
	for (j = 1; j <= 6; ++j) {
		resc12 += cheb12[k - 1] * chebmo[m + k * chebmo_dim1];
		ress12 += cheb12[k] * chebmo[m + (k + 1) * chebmo_dim1];
		k += -2;
		/* L140: */
	}
	resc24 = cheb24[24] * chebmo[m + chebmo_dim1 * 25];
	ress24 = 0.;
	*resabs = fabs(cheb24[24]);
	k = 23;
	for (j = 1; j <= 12; ++j) {
		resc24 += cheb24[k - 1] * chebmo[m + k * chebmo_dim1];
		ress24 += cheb24[k] * chebmo[m + (k + 1) * chebmo_dim1];
		*resabs = (d__1 = cheb24[k - 1], fabs(d__1)) + (d__2 = cheb24[k], fabs(
																			d__2));
		k += -2;
		/* L150: */
	}
	estc = (d__1 = resc24 - resc12, fabs(d__1));
	ests = (d__1 = ress24 - ress12, fabs(d__1));
	*resabs *= fabs(hlgth);
	if (*integr == 2) {
		goto L160;
	}
	*result = conc * resc24 - cons * ress24;
	*abserr = (d__1 = conc * estc, fabs(d__1)) + (d__2 = cons * ests, fabs(d__2)
		);
	goto L170;
L160:
	*result = conc * ress24 + cons * resc24;
	*abserr = (d__1 = conc * ests, fabs(d__1)) + (d__2 = cons * estc, fabs(d__2)
		);
L170:
	return;
} /* dqc25f_ */

void dqc25s_(const D_fp& f, const double *a, const double *b, 
				 const double *bl, const double *br, const double *alfa, 
				 const double *beta, const double *ri, const double *rj, 
				 const double *rg, const double *rh, 
				 double *result, double *abserr, double *resasc, 
				 const long *integr, long *nev)
{
	/* Initialized data */

	double x[11] = { .991444861373810411144557526928563,
									.965925826289068286749743199728897,
									.923879532511286756128183189396788,
									.866025403784438646763723170752936,
									.793353340291235164579776961501299,
									.707106781186547524400844362104849,
									.608761429008720639416097542898164,.5,
									.382683432365089771728459984030399,
									.258819045102520762348898837624048,
									.130526192220051591548406227895489 };

	/* System generated locals */
	double d__1, d__2;

	/* Local variables */
	long i__;
	double u, dc, fix, fval[25], res12, res24;
	long isym;
	double cheb12[13], cheb24[25], hlgth, centr;
	double factor, resabs;

	/* ***begin prologue  dqc25s */
	/* ***date written   810101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a2 */
	/* ***keywords  25-point clenshaw-curtis integration */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f*w over (bl,br), with error */
	/*            estimate, where the weight function w has a singular */
	/*            behaviour of algebraico-logarithmic type at the points */
	/*            a and/or b. (bl,br) is a part of (a,b). */
	/* ***description */

	/*        integration rules for integrands having algebraico-logarithmic */
	/*        end point singularities */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*           f      - double precision */
	/*                    function subprogram defining the integrand */
	/*                    f(x). the actual name for f needs to be declared */
	/*                    e x t e r n a l  in the driver program. */

	/*           a      - double precision */
	/*                    left end point of the original interval */

	/*           b      - double precision */
	/*                    right end point of the original interval, b.gt.a */

	/*           bl     - double precision */
	/*                    lower limit of integration, bl.ge.a */

	/*           br     - double precision */
	/*                    upper limit of integration, br.le.b */

	/*           alfa   - double precision */
	/*                    parameter in the weight function */

	/*           beta   - double precision */
	/*                    parameter in the weight function */

	/*           ri,rj,rg,rh - double precision */
	/*                    modified chebyshev moments for the application */
	/*                    of the generalized clenshaw-curtis */
	/*                    method (computed in subroutine dqmomo) */

	/*           result - double precision */
	/*                    approximation to the integral */
	/*                    result is computed by using a generalized */
	/*                    clenshaw-curtis method if b1 = a or br = b. */
	/*                    in all other cases the 15-point kronrod */
	/*                    rule is applied, obtained by optimal addition of */
	/*                    abscissae to the 7-point gauss rule. */

	/*           abserr - double precision */
	/*                    estimate of the modulus of the absolute error, */
	/*                    which should equal or exceed fabs(i-result) */

	/*           resasc - double precision */
	/*                    approximation to the integral of fabs(f*w-i/(b-a)) */

	/*           integr - long */
	/*                    which determines the weight function */
	/*                    = 1   w(x) = (x-a)**alfa*(b-x)**beta */
	/*                    = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a) */
	/*                    = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x) */
	/*                    = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)* */
	/*                                 log(b-x) */

	/*           nev    - long */
	/*                    number of integrand evaluations */
	/* ***references  (none) */
	/* ***routines called  dqcheb,dqk15w */
	/* ***end prologue  dqc25s */




	/*           the vector x contains the values cos(k*pi/24) */
	/*           k = 1, ..., 11, to be used for the computation of the */
	/*           chebyshev series expansion of f. */

	/* Parameter adjustments */
	--rh;
	--rg;
	--rj;
	--ri;

	/* Function Body */

	/*           list of major variables */
	/*           ----------------------- */

	/*           fval   - value of the function f at the points */
	/*                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5 */
	/*                    k = 0, ..., 24 */
	/*           cheb12 - coefficients of the chebyshev series expansion */
	/*                    of degree 12, for the function f, in the */
	/*                    interval (bl,br) */
	/*           cheb24 - coefficients of the chebyshev series expansion */
	/*                    of degree 24, for the function f, in the */
	/*                    interval (bl,br) */
	/*           res12  - approximation to the integral obtained from cheb12 */
	/*           res24  - approximation to the integral obtained from cheb24 */
	/*           dqwgts - external function subprogram defining */
	/*                    the four possible weight functions */
	/*           hlgth  - half-length of the interval (bl,br) */
	/*           centr  - mid point of the interval (bl,br) */

	/* ***first executable statement  dqc25s */
	*nev = 25;
	if (*bl == *a && (*alfa != 0. || *integr == 2 || *integr == 4)) {
		goto L10;
	}
	if (*br == *b && (*beta != 0. || *integr == 3 || *integr == 4)) {
		goto L140;
	}

	/*           if a.gt.bl and b.lt.br, apply the 15-point gauss-kronrod */
	/*           scheme. */


	dqk15w_(f, dqwgts_, a, b, alfa, beta, integr, bl, br, result, 
			  abserr, &resabs, resasc);
	*nev = 15;
	goto L270;

	/*           this part of the program is executed only if a = bl. */
	/*           ---------------------------------------------------- */

	/*           compute the chebyshev series expansion of the */
	/*           following function */
	/*           f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta */
	/*                  *f(0.5*(br-a)*x+0.5*(br+a)) */

L10:
	hlgth = (*br - *bl) * .5;
	centr = (*br + *bl) * .5;
	fix = *b - centr;
	d__1 = hlgth + centr;
	d__2 = fix - hlgth;
	fval[0] = f(d__1) * .5 * pow(d__2, *beta);
	fval[12] = f(centr) * pow(fix, *beta);
	d__1 = centr - hlgth;
	d__2 = fix + hlgth;
	fval[24] = f(d__1) * .5 * pow(d__2, *beta);
	for (i__ = 2; i__ <= 12; ++i__) {
		u = hlgth * x[i__ - 2];
		isym = 26 - i__;
		d__1 = u + centr;
		d__2 = fix - u;
		fval[i__ - 1] = f(d__1) * pow(d__2, *beta);
		d__1 = centr - u;
		d__2 = fix + u;
		fval[isym - 1] = f(d__1) * pow(d__2, *beta);
		/* L20: */
	}
	d__1 = *alfa + 1.;
	factor = pow(hlgth, d__1);
	*result = 0.;
	*abserr = 0.;
	res12 = 0.;
	res24 = 0.;
	if (*integr > 2) {
		goto L70;
	}
	dqcheb_(x, fval, cheb12, cheb24);

	/*           integr = 1  (or 2) */

	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * ri[i__];
		res24 += cheb24[i__ - 1] * ri[i__];
		/* L30: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * ri[i__];
		/* L40: */
	}
	if (*integr == 1) {
		goto L130;
	}

	/*           integr = 2 */

	dc = log(*br - *bl);
	*result = res24 * dc;
	*abserr = (d__1 = (res24 - res12) * dc, fabs(d__1));
	res12 = 0.;
	res24 = 0.;
	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rg[i__];
		res24 = res12 + cheb24[i__ - 1] * rg[i__];
		/* L50: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rg[i__];
		/* L60: */
	}
	goto L130;

	/*           compute the chebyshev series expansion of the */
	/*           following function */
	/*           f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x) */

L70:
	fval[0] *= log(fix - hlgth);
	fval[12] *= log(fix);
	fval[24] *= log(fix + hlgth);
	for (i__ = 2; i__ <= 12; ++i__) {
		u = hlgth * x[i__ - 2];
		isym = 26 - i__;
		fval[i__ - 1] *= log(fix - u);
		fval[isym - 1] *= log(fix + u);
		/* L80: */
	}
	dqcheb_(x, fval, cheb12, cheb24);

	/*           integr = 3  (or 4) */

	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * ri[i__];
		res24 += cheb24[i__ - 1] * ri[i__];
		/* L90: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * ri[i__];
		/* L100: */
	}
	if (*integr == 3) {
		goto L130;
	}

	/*           integr = 4 */

	dc = log(*br - *bl);
	*result = res24 * dc;
	*abserr = (d__1 = (res24 - res12) * dc, fabs(d__1));
	res12 = 0.;
	res24 = 0.;
	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rg[i__];
		res24 += cheb24[i__ - 1] * rg[i__];
		/* L110: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rg[i__];
		/* L120: */
	}
L130:
	*result = (*result + res24) * factor;
	*abserr = (*abserr + (d__1 = res24 - res12, fabs(d__1))) * factor;
	goto L270;

	/*           this part of the program is executed only if b = br. */
	/*           ---------------------------------------------------- */

	/*           compute the chebyshev series expansion of the */
	/*           following function */
	/*           f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa */
	/*                *f(0.5*(b-bl)*x+0.5*(b+bl)) */

L140:
	hlgth = (*br - *bl) * .5;
	centr = (*br + *bl) * .5;
	fix = centr - *a;
	d__1 = hlgth + centr;
	d__2 = fix + hlgth;
	fval[0] = f(d__1) * .5 * pow(d__2, *alfa);
	fval[12] = f(centr) * pow(fix, *alfa);
	d__1 = centr - hlgth;
	d__2 = fix - hlgth;
	fval[24] = f(d__1) * .5 * pow(d__2, *alfa);
	for (i__ = 2; i__ <= 12; ++i__) {
		u = hlgth * x[i__ - 2];
		isym = 26 - i__;
		d__1 = u + centr;
		d__2 = fix + u;
		fval[i__ - 1] = f(d__1) * pow(d__2, *alfa);
		d__1 = centr - u;
		d__2 = fix - u;
		fval[isym - 1] = f(d__1) * pow(d__2, *alfa);
		/* L150: */
	}
	d__1 = *beta + 1.;
	factor = pow(hlgth, d__1);
	*result = 0.;
	*abserr = 0.;
	res12 = 0.;
	res24 = 0.;
	if (*integr == 2 || *integr == 4) {
		goto L200;
	}

	/*           integr = 1  (or 3) */

	dqcheb_(x, fval, cheb12, cheb24);
	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rj[i__];
		res24 += cheb24[i__ - 1] * rj[i__];
		/* L160: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rj[i__];
		/* L170: */
	}
	if (*integr == 1) {
		goto L260;
	}

	/*           integr = 3 */

	dc = log(*br - *bl);
	*result = res24 * dc;
	*abserr = (d__1 = (res24 - res12) * dc, fabs(d__1));
	res12 = 0.;
	res24 = 0.;
	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rh[i__];
		res24 += cheb24[i__ - 1] * rh[i__];
		/* L180: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rh[i__];
		/* L190: */
	}
	goto L260;

	/*           compute the chebyshev series expansion of the */
	/*           following function */
	/*           f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a)) */

L200:
	fval[0] *= log(hlgth + fix);
	fval[12] *= log(fix);
	fval[24] *= log(fix - hlgth);
	for (i__ = 2; i__ <= 12; ++i__) {
		u = hlgth * x[i__ - 2];
		isym = 26 - i__;
		fval[i__ - 1] *= log(u + fix);
		fval[isym - 1] *= log(fix - u);
		/* L210: */
	}
	dqcheb_(x, fval, cheb12, cheb24);

	/*           integr = 2  (or 4) */

	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rj[i__];
		res24 += cheb24[i__ - 1] * rj[i__];
		/* L220: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rj[i__];
		/* L230: */
	}
	if (*integr == 2) {
		goto L260;
	}
	dc = log(*br - *bl);
	*result = res24 * dc;
	*abserr = (d__1 = (res24 - res12) * dc, fabs(d__1));
	res12 = 0.;
	res24 = 0.;

	/*           integr = 4 */

	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rh[i__];
		res24 += cheb24[i__ - 1] * rh[i__];
		/* L240: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rh[i__];
		/* L250: */
	}
L260:
	*result = (*result + res24) * factor;
	*abserr = (*abserr + (d__1 = res24 - res12, fabs(d__1))) * factor;
L270:
	return;
} /* dqc25s_ */

void dqcheb_(const double *x, double *fval, double *cheb12, double *cheb24)
{
	long i__, j;
	double v[12], alam, alam1, alam2, part1, part2, part3;

	/* ***begin prologue  dqcheb */
	/* ***refer to  dqc25c,dqc25f,dqc25s */
	/* ***routines called  (none) */
	/* ***revision date  830518   (yymmdd) */
	/* ***keywords  chebyshev series expansion, fast fourier transform */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  this routine computes the chebyshev series expansion */
	/*            of degrees 12 and 24 of a function using a */
	/*            fast fourier transform method */
	/*            f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)), */
	/*            f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)), */
	/*            where t(k,x) is the chebyshev polynomial of degree k. */
	/* ***description */

	/*        chebyshev series expansion */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*          on entry */
	/*           x      - double precision */
	/*                    vector of dimension 11 containing the */
	/*                    values cos(k*pi/24), k = 1, ..., 11 */

	/*           fval   - double precision */
	/*                    vector of dimension 25 containing the */
	/*                    function values at the points */
	/*                    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24, */
	/*                    where (a,b) is the approximation interval. */
	/*                    fval(1) and fval(25) are divided by two */
	/*                    (these values are destroyed at output). */

	/*          on return */
	/*           cheb12 - double precision */
	/*                    vector of dimension 13 containing the */
	/*                    chebyshev coefficients for degree 12 */

	/*           cheb24 - double precision */
	/*                    vector of dimension 25 containing the */
	/*                    chebyshev coefficients for degree 24 */

	/* ***end prologue  dqcheb */



	/* ***first executable statement  dqcheb */
	/* Parameter adjustments */
	--cheb24;
	--cheb12;
	--fval;
	--x;

	/* Function Body */
	for (i__ = 1; i__ <= 12; ++i__) {
		j = 26 - i__;
		v[i__ - 1] = fval[i__] - fval[j];
		fval[i__] += fval[j];
		/* L10: */
	}
	alam1 = v[0] - v[8];
	alam2 = x[6] * (v[2] - v[6] - v[10]);
	cheb12[4] = alam1 + alam2;
	cheb12[10] = alam1 - alam2;
	alam1 = v[1] - v[7] - v[9];
	alam2 = v[3] - v[5] - v[11];
	alam = x[3] * alam1 + x[9] * alam2;
	cheb24[4] = cheb12[4] + alam;
	cheb24[22] = cheb12[4] - alam;
	alam = x[9] * alam1 - x[3] * alam2;
	cheb24[10] = cheb12[10] + alam;
	cheb24[16] = cheb12[10] - alam;
	part1 = x[4] * v[4];
	part2 = x[8] * v[8];
	part3 = x[6] * v[6];
	alam1 = v[0] + part1 + part2;
	alam2 = x[2] * v[2] + part3 + x[10] * v[10];
	cheb12[2] = alam1 + alam2;
	cheb12[12] = alam1 - alam2;
	alam = x[1] * v[1] + x[3] * v[3] + x[5] * v[5] + x[7] * v[7] + x[9] * v[9]
		+ x[11] * v[11];
	cheb24[2] = cheb12[2] + alam;
	cheb24[24] = cheb12[2] - alam;
	alam = x[11] * v[1] - x[9] * v[3] + x[7] * v[5] - x[5] * v[7] + x[3] * v[
		9] - x[1] * v[11];
	cheb24[12] = cheb12[12] + alam;
	cheb24[14] = cheb12[12] - alam;
	alam1 = v[0] - part1 + part2;
	alam2 = x[10] * v[2] - part3 + x[2] * v[10];
	cheb12[6] = alam1 + alam2;
	cheb12[8] = alam1 - alam2;
	alam = x[5] * v[1] - x[9] * v[3] - x[1] * v[5] - x[11] * v[7] + x[3] * v[
		9] + x[7] * v[11];
	cheb24[6] = cheb12[6] + alam;
	cheb24[20] = cheb12[6] - alam;
	alam = x[7] * v[1] - x[3] * v[3] - x[11] * v[5] + x[1] * v[7] - x[9] * v[
		9] - x[5] * v[11];
	cheb24[8] = cheb12[8] + alam;
	cheb24[18] = cheb12[8] - alam;
	for (i__ = 1; i__ <= 6; ++i__) {
		j = 14 - i__;
		v[i__ - 1] = fval[i__] - fval[j];
		fval[i__] += fval[j];
		/* L20: */
	}
	alam1 = v[0] + x[8] * v[4];
	alam2 = x[4] * v[2];
	cheb12[3] = alam1 + alam2;
	cheb12[11] = alam1 - alam2;
	cheb12[7] = v[0] - v[4];
	alam = x[2] * v[1] + x[6] * v[3] + x[10] * v[5];
	cheb24[3] = cheb12[3] + alam;
	cheb24[23] = cheb12[3] - alam;
	alam = x[6] * (v[1] - v[3] - v[5]);
	cheb24[7] = cheb12[7] + alam;
	cheb24[19] = cheb12[7] - alam;
	alam = x[10] * v[1] - x[6] * v[3] + x[2] * v[5];
	cheb24[11] = cheb12[11] + alam;
	cheb24[15] = cheb12[11] - alam;
	for (i__ = 1; i__ <= 3; ++i__) {
		j = 8 - i__;
		v[i__ - 1] = fval[i__] - fval[j];
		fval[i__] += fval[j];
		/* L30: */
	}
	cheb12[5] = v[0] + x[8] * v[2];
	cheb12[9] = fval[1] - x[8] * fval[3];
	alam = x[4] * v[1];
	cheb24[5] = cheb12[5] + alam;
	cheb24[21] = cheb12[5] - alam;
	alam = x[8] * fval[2] - fval[4];
	cheb24[9] = cheb12[9] + alam;
	cheb24[17] = cheb12[9] - alam;
	cheb12[1] = fval[1] + fval[3];
	alam = fval[2] + fval[4];
	cheb24[1] = cheb12[1] + alam;
	cheb24[25] = cheb12[1] - alam;
	cheb12[13] = v[0] - v[2];
	cheb24[13] = cheb12[13];
	alam = .16666666666666666;
	for (i__ = 2; i__ <= 12; ++i__) {
		cheb12[i__] *= alam;
		/* L40: */
	}
	alam *= .5;
	cheb12[1] *= alam;
	cheb12[13] *= alam;
	for (i__ = 2; i__ <= 24; ++i__) {
		cheb24[i__] *= alam;
		/* L50: */
	}
	cheb24[1] = alam * .5 * cheb24[1];
	cheb24[25] = alam * .5 * cheb24[25];
	return;
} /* dqcheb_ */

void dqelg_(long *n, double *epstab, double *result, 
				double *abserr, double *res3la, long *nres)
{
	/* System generated locals */
	int i__1;
	double d__1, d__2, d__3;

	/* Local variables */
	long i__;
	double e0, e1, e2, e3;
	long k1, k2, k3, ib, ie;
	double ss;
	long ib2;
	double res;
	long num;
	double err1, err2, err3, tol1, tol2, tol3;
	long indx;
	double e1abs, oflow, error;
	double delta1, delta2, delta3, epmach, epsinf;
	long newelm, limexp;

	/* ***begin prologue  dqelg */
	/* ***refer to  dqagie,dqagoe,dqagpe,dqagse */
	/* ***routines called  d1mach */
	/* ***revision date  830518   (yymmdd) */
	/* ***keywords  epsilon algorithm, convergence acceleration, */
	/*             extrapolation */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
	/* ***purpose  the routine determines the limit of a given sequence of */
	/*            approximations, by means of the epsilon algorithm of */
	/*            p.wynn. an estimate of the absolute error is also given. */
	/*            the condensed epsilon table is computed. only those */
	/*            elements needed for the computation of the next diagonal */
	/*            are preserved. */
	/* ***description */

	/*           epsilon algorithm */
	/*           standard fortran subroutine */
	/*           double precision version */

	/*           parameters */
	/*              n      - long */
	/*                       epstab(n) contains the new element in the */
	/*                       first column of the epsilon table. */

	/*              epstab - double precision */
	/*                       vector of dimension 52 containing the elements */
	/*                       of the two lower diagonals of the triangular */
	/*                       epsilon table. the elements are numbered */
	/*                       starting at the right-hand corner of the */
	/*                       triangle. */

	/*              result - double precision */
	/*                       resulting approximation to the integral */

	/*              abserr - double precision */
	/*                       estimate of the absolute error computed from */
	/*                       result and the 3 previous results */

	/*              res3la - double precision */
	/*                       vector of dimension 3 containing the last 3 */
	/*                       results */

	/*              nres   - long */
	/*                       number of calls to the routine */
	/*                       (should be zero at first call) */

	/* ***end prologue  dqelg */


	/*           list of major variables */
	/*           ----------------------- */

	/*           e0     - the 4 elements on which the computation of a new */
	/*           e1       element in the epsilon table is based */
	/*           e2 */
	/*           e3                 e0 */
	/*                        e3    e1    new */
	/*                              e2 */
	/*           newelm - number of elements to be computed in the new */
	/*                    diagonal */
	/*           error  - error = fabs(e1-e0)+fabs(e2-e1)+fabs(new-e2) */
	/*           result - the element in the new diagonal with least value */
	/*                    of error */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           oflow is the largest positive magnitude. */
	/*           limexp is the maximum number of elements the epsilon */
	/*           table can contain. if this number is reached, the upper */
	/*           diagonal of the epsilon table is deleted. */

	/* ***first executable statement  dqelg */
	/* Parameter adjustments */
	--res3la;
	--epstab;

	/* Function Body */
	epmach = d1mach(c__4);
	oflow = d1mach(c__2);
	++(*nres);
	*abserr = oflow;
	*result = epstab[*n];
	if (*n < 3) {
		goto L100;
	}
	limexp = 50;
	epstab[*n + 2] = epstab[*n];
	newelm = (*n - 1) / 2;
	epstab[*n] = oflow;
	num = *n;
	k1 = *n;
	i__1 = newelm;
	for (i__ = 1; i__ <= i__1; ++i__) {
		k2 = k1 - 1;
		k3 = k1 - 2;
		res = epstab[k1 + 2];
		e0 = epstab[k3];
		e1 = epstab[k2];
		e2 = res;
		e1abs = fabs(e1);
		delta2 = e2 - e1;
		err2 = fabs(delta2);
		/* Computing MAX */
		d__1 = fabs(e2);
		tol2 = max(d__1,e1abs) * epmach;
		delta3 = e1 - e0;
		err3 = fabs(delta3);
		/* Computing MAX */
		d__1 = e1abs, d__2 = fabs(e0);
		tol3 = max(d__1,d__2) * epmach;
		if (err2 > tol2 || err3 > tol3) {
			goto L10;
		}

		/*           if e0, e1 and e2 are equal to within machine */
		/*           accuracy, convergence is assumed. */
		/*           result = e2 */
		/*           abserr = fabs(e1-e0)+fabs(e2-e1) */

		*result = res;
		*abserr = err2 + err3;
		/* ***jump out of do-loop */
		goto L100;
	L10:
		e3 = epstab[k1];
		epstab[k1] = e1;
		delta1 = e1 - e3;
		err1 = fabs(delta1);
		/* Computing MAX */
		d__1 = e1abs, d__2 = fabs(e3);
		tol1 = max(d__1,d__2) * epmach;

		/*           if two elements are very close to each other, omit */
		/*           a part of the table by adjusting the value of n */

		if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
			goto L20;
		}
		ss = 1. / delta1 + 1. / delta2 - 1. / delta3;
		epsinf = (d__1 = ss * e1, fabs(d__1));

		/*           test to detect irregular behaviour in the table, and */
		/*           eventually omit a part of the table adjusting the value */
		/*           of n. */

		if (epsinf > 1e-4) {
			goto L30;
		}
	L20:
		*n = i__ + i__ - 1;
		/* ***jump out of do-loop */
		goto L50;

		/*           compute a new element and eventually adjust */
		/*           the value of result. */

	L30:
		res = e1 + 1. / ss;
		epstab[k1] = res;
		k1 += -2;
		error = err2 + (d__1 = res - e2, fabs(d__1)) + err3;
		if (error > *abserr) {
			goto L40;
		}
		*abserr = error;
		*result = res;
	L40:
		;
	}

	/*           shift the table. */

L50:
	if (*n == limexp) {
		*n = (limexp / 2 << 1) - 1;
	}
	ib = 1;
	if (num / 2 << 1 == num) {
		ib = 2;
	}
	ie = newelm + 1;
	i__1 = ie;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ib2 = ib + 2;
		epstab[ib] = epstab[ib2];
		ib = ib2;
		/* L60: */
	}
	if (num == *n) {
		goto L80;
	}
	indx = num - *n + 1;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		epstab[i__] = epstab[indx];
		++indx;
		/* L70: */
	}
L80:
	if (*nres >= 4) {
		goto L90;
	}
	res3la[*nres] = *result;
	*abserr = oflow;
	goto L100;

	/*           compute error estimate */

L90:
	*abserr = (d__1 = *result - res3la[3], fabs(d__1)) + (d__2 = *result - 
																			res3la[2], fabs(d__2)) + (d__3 = *result - res3la[1], fabs(d__3));
	res3la[1] = res3la[2];
	res3la[2] = res3la[3];
	res3la[3] = *result;
L100:
	/* Computing MAX */
	d__1 = *abserr, d__2 = epmach * 5. * fabs(*result);
	*abserr = max(d__1,d__2);
	return;
} /* dqelg_ */

void dqk15_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc)
{
	/* Initialized data */

	double wg[4] = { .129484966168869693270611432679082,
									.27970539148927666790146777142378,
									.381830050505118944950369775488975,
									.417959183673469387755102040816327 };
	double xgk[8] = { .991455371120812639206854697526329,
									 .949107912342758524526189684047851,
									 .864864423359769072789712788640926,
									 .741531185599394439863864773280788,
									 .58608723546769113029414483825873,
									 .405845151377397166906606412076961,
									 .207784955007898467600689403773245,0. };
	double wgk[8] = { .02293532201052922496373200805897,
									 .063092092629978553290700663189204,
									 .104790010322250183839876322541518,
									 .140653259715525918745189590510238,
									 .16900472663926790282658342659855,
									 .190350578064785409913256402421014,
									 .204432940075298892414161999234649,
									 .209482141084727828012999174891714 };

	/* System generated locals */
	double d__1, d__2, d__3;

	/* Local variables */
	long j;
	double fc, fv1[7], fv2[7];
	long jtw;
	double absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	double hlgth, centr, reskh, uflow;
	double epmach, dhlgth;

	/* ***begin prologue  dqk15 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  15-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b), with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           double precision version */

	/*           parameters */
	/*            on entry */
	/*              f      - double precision */
	/*                       function subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the calling program. */

	/*              a      - double precision */
	/*                       lower limit of integration */

	/*              b      - double precision */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - double precision */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 15-polong */
	/*                       kronrod rule (resk) obtained by optimal addition */
	/*                       of abscissae to the7-point gauss rule(resg). */

	/*              abserr - double precision */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should not exceed fabs(i-result) */

	/*              resabs - double precision */
	/*                       approximation to the integral j */

	/*              resasc - double precision */
	/*                       approximation to the integral of fabs(f-i/(b-a)) */
	/*                       over (a,b) */

	/* ***references  (none) */
	/* ***routines called  d1mach */
	/* ***end prologue  dqk15 */



	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 15-point kronrod rule */
	/*                    xgk(2), xgk(4), ...  abscissae of the 7-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 7-point gauss rule */

	/*           wgk    - weights of the 15-point kronrod rule */

	/*           wg     - weights of the 7-point gauss rule */


	/* gauss quadrature weights and kronron quadrature abscissae and weights */
	/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
	/* bell labs, nov. 1981. */





	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 7-point gauss formula */
	/*           resk   - result of the 15-point kronrod formula */
	/*           reskh  - approximation to the mean value of f over (a,b), */
	/*                    i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  dqk15 */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);

	centr = (*a + *b) * .5;
	hlgth = (*b - *a) * .5;
	dhlgth = fabs(hlgth);

	/*           compute the 15-point kronrod approximation to */
	/*           the integral, and estimate the absolute error. */

	fc = f(centr);
	resg = fc * wg[3];
	resk = fc * wgk[7];
	*resabs = fabs(resk);
	for (j = 1; j <= 3; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		d__1 = centr - absc;
		fval1 = f(d__1);
		d__1 = centr + absc;
		fval2 = f(d__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 4; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		d__1 = centr - absc;
		fval1 = f(d__1);
		d__1 = centr + absc;
		fval2 = f(d__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5;
	*resasc = wgk[7] * (d__1 = fc - reskh, fabs(d__1));
	for (j = 1; j <= 7; ++j) {
		*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, fabs(d__1)) + (
											 d__2 = fv2[j - 1] - reskh, fabs(d__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (d__1 = (resk - resg) * hlgth, fabs(d__1));
	if (*resasc != 0. && *abserr != 0.) {
		/* Computing MIN */
		d__3 = *abserr * 200. / *resasc;
		d__1 = 1., d__2 = pow(d__3, c_b270);
		*abserr = *resasc * min(d__1,d__2);
	}
	if (*resabs > uflow / (epmach * 50.)) {
		/* Computing MAX */
		d__1 = epmach * 50. * *resabs;
		*abserr = max(d__1,*abserr);
	}
	return;
} /* dqk15_ */

void dqk15i_(const D_fp& f, const double *boun, const long *inf, 
				 const double *a, const double *b, double *result, double *abserr, 
				 double *resabs, double *resasc)
{
	/* Initialized data */

	double wg[8] = { 0.,.129484966168869693270611432679082,0.,
									.27970539148927666790146777142378,0.,
									.381830050505118944950369775488975,0.,
									.417959183673469387755102040816327 };
	double xgk[8] = { .991455371120812639206854697526329,
									 .949107912342758524526189684047851,
									 .864864423359769072789712788640926,
									 .741531185599394439863864773280788,
									 .58608723546769113029414483825873,
									 .405845151377397166906606412076961,
									 .207784955007898467600689403773245,0. };
	double wgk[8] = { .02293532201052922496373200805897,
									 .063092092629978553290700663189204,
									 .104790010322250183839876322541518,
									 .140653259715525918745189590510238,
									 .16900472663926790282658342659855,
									 .190350578064785409913256402421014,
									 .204432940075298892414161999234649,
									 .209482141084727828012999174891714 };

	/* System generated locals */
	double d__1, d__2, d__3;

	/* Local variables */
	long j;
	double fc, fv1[7], fv2[7], absc, dinf, resg, resk, fsum, absc1,
		absc2, fval1, fval2, hlgth, centr, reskh, uflow;
	double tabsc1, tabsc2, epmach;

	/* ***begin prologue  dqk15i */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a3a2,h2a4a2 */
	/* ***keywords  15-point transformed gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the original (infinite integration range is mapped */
	/*            onto the interval (0,1) and (a,b) is a part of (0,1). */
	/*            it is the purpose to compute */
	/*            i = integral of transformed integrand over (a,b), */
	/*            j = integral of fabs(transformed integrand) over (a,b). */
	/* ***description */

	/*           integration rule */
	/*           standard fortran subroutine */
	/*           double precision version */

	/*           parameters */
	/*            on entry */
	/*              f      - double precision */
	/*                       fuction subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the calling program. */

	/*              boun   - double precision */
	/*                       finite bound of original integration */
	/*                       range (set to zero if inf = +2) */

	/*              inf    - long */
	/*                       if inf = -1, the original interval is */
	/*                                   (-infinity,bound), */
	/*                       if inf = +1, the original interval is */
	/*                                   (bound,+infinity), */
	/*                       if inf = +2, the original interval is */
	/*                                   (-infinity,+infinity) and */
	/*                       the integral is computed as the sum of two */
	/*                       integrals, one over (-infinity,0) and one over */
	/*                       (0,+infinity). */

	/*              a      - double precision */
	/*                       lower limit for integration over subrange */
	/*                       of (0,1) */

	/*              b      - double precision */
	/*                       upper limit for integration over subrange */
	/*                       of (0,1) */

	/*            on return */
	/*              result - double precision */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 15-polong */
	/*                       kronrod rule(resk) obtained by optimal addition */
	/*                       of abscissae to the 7-point gauss rule(resg). */

	/*              abserr - double precision */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should equal or exceed fabs(i-result) */

	/*              resabs - double precision */
	/*                       approximation to the integral j */

	/*              resasc - double precision */
	/*                       approximation to the integral of */
	/*                       fabs((transformed integrand)-i/(b-a)) over (a,b) */

	/* ***references  (none) */
	/* ***routines called  d1mach */
	/* ***end prologue  dqk15i */



	/*           the abscissae and weights are supplied for the interval */
	/*           (-1,1).  because of symmetry only the positive abscissae and */
	/*           their corresponding weights are given. */

	/*           xgk    - abscissae of the 15-point kronrod rule */
	/*                    xgk(2), xgk(4), ... abscissae of the 7-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 7-point gauss rule */

	/*           wgk    - weights of the 15-point kronrod rule */

	/*           wg     - weights of the 7-point gauss rule, corresponding */
	/*                    to the abscissae xgk(2), xgk(4), ... */
	/*                    wg(1), wg(3), ... are set to zero. */





	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc*  - abscissa */
	/*           tabsc* - transformed abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 7-point gauss formula */
	/*           resk   - result of the 15-point kronrod formula */
	/*           reskh  - approximation to the mean value of the transformed */
	/*                    integrand over (a,b), i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  dqk15i */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);
	dinf = (double) min(1,*inf);

	centr = (*a + *b) * .5;
	hlgth = (*b - *a) * .5;
	tabsc1 = *boun + dinf * (1. - centr) / centr;
	fval1 = f(tabsc1);
	if (*inf == 2) {
		d__1 = -tabsc1;
		fval1 += f(d__1);
	}
	fc = fval1 / centr / centr;

	/*           compute the 15-point kronrod approximation to */
	/*           the integral, and estimate the error. */

	resg = wg[7] * fc;
	resk = wgk[7] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 7; ++j) {
		absc = hlgth * xgk[j - 1];
		absc1 = centr - absc;
		absc2 = centr + absc;
		tabsc1 = *boun + dinf * (1. - absc1) / absc1;
		tabsc2 = *boun + dinf * (1. - absc2) / absc2;
		fval1 = f(tabsc1);
		fval2 = f(tabsc2);
		if (*inf == 2) {
			d__1 = -tabsc1;
			fval1 += f(d__1);
		}
		if (*inf == 2) {
			d__1 = -tabsc2;
			fval2 += f(d__1);
		}
		fval1 = fval1 / absc1 / absc1;
		fval2 = fval2 / absc2 / absc2;
		fv1[j - 1] = fval1;
		fv2[j - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[j - 1] * fsum;
		*resabs += wgk[j - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	reskh = resk * .5;
	*resasc = wgk[7] * (d__1 = fc - reskh, fabs(d__1));
	for (j = 1; j <= 7; ++j) {
		*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, fabs(d__1)) + (
											 d__2 = fv2[j - 1] - reskh, fabs(d__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resasc *= hlgth;
	*resabs *= hlgth;
	*abserr = (d__1 = (resk - resg) * hlgth, fabs(d__1));
	if (*resasc != 0. && *abserr != 0.) {
		/* Computing MIN */
		d__3 = *abserr * 200. / *resasc;
		d__1 = 1., d__2 = pow(d__3, c_b270);
		*abserr = *resasc * min(d__1,d__2);
	}
	if (*resabs > uflow / (epmach * 50.)) {
		/* Computing MAX */
		d__1 = epmach * 50. * *resabs;
		*abserr = max(d__1,*abserr);
	}
	return;
} /* dqk15i_ */

void dqk15w_(const D_fp& f, D_fp1 w, const double *p1, const double *p2, 
				 const double *p3, const double *p4, const long *kp, 
				 const double *a, const double *b, 
				 double *result, double *abserr, double *resabs, double *resasc)
{
	/* Initialized data */

	double xgk[8] = { .9914553711208126,.9491079123427585,
									 .8648644233597691,.7415311855993944,.5860872354676911,
									 .4058451513773972,.2077849550078985,0. };
	double wgk[8] = { .02293532201052922,.06309209262997855,
									 .1047900103222502,.1406532597155259,.1690047266392679,
									 .1903505780647854,.2044329400752989,.2094821410847278 };
	double wg[4] = { .1294849661688697,.2797053914892767,
									.3818300505051889,.4179591836734694 };

	/* System generated locals */
	double d__1, d__2, d__3;

	/* Local variables */
	long j;
	double fc, fv1[7], fv2[7];
	long jtw;
	double absc, resg, resk, fsum, absc1, absc2, fval1, fval2;
	long jtwm1;
	double hlgth, centr, reskh, uflow;
	double epmach, dhlgth;

	/* ***begin prologue  dqk15w */
	/* ***date written   810101   (yymmdd) */
	/* ***revision date  830518   (mmddyy) */
	/* ***category no.  h2a2a2 */
	/* ***keywords  15-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f*w over (a,b), with error */
	/*                           estimate */
	/*                       j = integral of fabs(f*w) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           double precision version */

	/*           parameters */
	/*             on entry */
	/*              f      - double precision */
	/*                       function subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the driver program. */

	/*              w      - double precision */
	/*                       function subprogram defining the integrand */
	/*                       weight function w(x). the actual name for w */
	/*                       needs to be declared e x t e r n a l in the */
	/*                       calling program. */

	/*              p1, p2, p3, p4 - double precision */
	/*                       parameters in the weight function */

	/*              kp     - long */
	/*                       key for indicating the type of weight function */

	/*              a      - double precision */
	/*                       lower limit of integration */

	/*              b      - double precision */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - double precision */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 15-polong */
	/*                       kronrod rule (resk) obtained by optimal addition */
	/*                       of abscissae to the 7-point gauss rule (resg). */

	/*              abserr - double precision */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should equal or exceed fabs(i-result) */

	/*              resabs - double precision */
	/*                       approximation to the integral of fabs(f) */

	/*              resasc - double precision */
	/*                       approximation to the integral of fabs(f-i/(b-a)) */


	/* ***references  (none) */
	/* ***routines called  d1mach */
	/* ***end prologue  dqk15w */



	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 15-point gauss-kronrod rule */
	/*                    xgk(2), xgk(4), ... abscissae of the 7-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ... abscissae which are optimally */
	/*                    added to the 7-point gauss rule */

	/*           wgk    - weights of the 15-point gauss-kronrod rule */

	/*           wg     - weights of the 7-point gauss rule */





	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc*  - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 7-point gauss formula */
	/*           resk   - result of the 15-point kronrod formula */
	/*           reskh  - approximation to the mean value of f*w over (a,b), */
	/*                    i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  dqk15w */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);

	centr = (*a + *b) * .5;
	hlgth = (*b - *a) * .5;
	dhlgth = fabs(hlgth);

	/*           compute the 15-point kronrod approximation to the */
	/*           integral, and estimate the error. */

	fc = f(centr) * w(&centr, p1, p2, p3, p4, kp);
	resg = wg[3] * fc;
	resk = wgk[7] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 3; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		absc1 = centr - absc;
		absc2 = centr + absc;
		fval1 = f(absc1) * w(&absc1, p1, p2, p3, p4, kp);
		fval2 = f(absc2) * w(&absc2, p1, p2, p3, p4, kp);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 4; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		absc1 = centr - absc;
		absc2 = centr + absc;
		fval1 = f(absc1) * w(&absc1, p1, p2, p3, p4, kp);
		fval2 = f(absc2) * w(&absc2, p1, p2, p3, p4, kp);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5;
	*resasc = wgk[7] * (d__1 = fc - reskh, fabs(d__1));
	for (j = 1; j <= 7; ++j) {
		*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, fabs(d__1)) + (
											 d__2 = fv2[j - 1] - reskh, fabs(d__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (d__1 = (resk - resg) * hlgth, fabs(d__1));
	if (*resasc != 0. && *abserr != 0.) {
		/* Computing MIN */
		d__3 = *abserr * 200. / *resasc;
		d__1 = 1., d__2 = pow(d__3, c_b270);
		*abserr = *resasc * min(d__1,d__2);
	}
	if (*resabs > uflow / (epmach * 50.)) {
		/* Computing MAX */
		d__1 = epmach * 50. * *resabs;
		*abserr = max(d__1,*abserr);
	}
	return;
} /* dqk15w_ */

void dqk21_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc)
{
	/* Initialized data */

	double wg[5] = { .066671344308688137593568809893332,
									.149451349150580593145776339657697,
									.219086362515982043995534934228163,
									.269266719309996355091226921569469,
									.295524224714752870173892994651338 };
	double xgk[11] = { .995657163025808080735527280689003,
									  .973906528517171720077964012084452,
									  .930157491355708226001207180059508,
									  .865063366688984510732096688423493,
									  .780817726586416897063717578345042,
									  .679409568299024406234327365114874,
									  .562757134668604683339000099272694,
									  .433395394129247190799265943165784,
									  .294392862701460198131126603103866,
									  .14887433898163121088482600112972,0. };
	double wgk[11] = { .011694638867371874278064396062192,
									  .03255816230796472747881897245939,
									  .05475589657435199603138130024458,
									  .07503967481091995276704314091619,
									  .093125454583697605535065465083366,
									  .109387158802297641899210590325805,
									  .123491976262065851077958109831074,
									  .134709217311473325928054001771707,
									  .142775938577060080797094273138717,
									  .147739104901338491374841515972068,
									  .149445554002916905664936468389821 };

	/* System generated locals */
	double d__1, d__2, d__3;

	/* Local variables */
	long j;
	double fc, fv1[10], fv2[10];
	long jtw;
	double absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	double hlgth, centr, reskh, uflow;
	double epmach, dhlgth;

	/* ***begin prologue  dqk21 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  21-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b), with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           double precision version */

	/*           parameters */
	/*            on entry */
	/*              f      - double precision */
	/*                       function subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the driver program. */

	/*              a      - double precision */
	/*                       lower limit of integration */

	/*              b      - double precision */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - double precision */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 21-polong */
	/*                       kronrod rule (resk) obtained by optimal addition */
	/*                       of abscissae to the 10-point gauss rule (resg). */

	/*              abserr - double precision */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should not exceed fabs(i-result) */

	/*              resabs - double precision */
	/*                       approximation to the integral j */

	/*              resasc - double precision */
	/*                       approximation to the integral of fabs(f-i/(b-a)) */
	/*                       over (a,b) */

	/* ***references  (none) */
	/* ***routines called  d1mach */
	/* ***end prologue  dqk21 */



	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 21-point kronrod rule */
	/*                    xgk(2), xgk(4), ...  abscissae of the 10-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 10-point gauss rule */

	/*           wgk    - weights of the 21-point kronrod rule */

	/*           wg     - weights of the 10-point gauss rule */


	/* gauss quadrature weights and kronron quadrature abscissae and weights */
	/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
	/* bell labs, nov. 1981. */





	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 10-point gauss formula */
	/*           resk   - result of the 21-point kronrod formula */
	/*           reskh  - approximation to the mean value of f over (a,b), */
	/*                    i.e. to i/(b-a) */


	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  dqk21 */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);

	centr = (*a + *b) * .5;
	hlgth = (*b - *a) * .5;
	dhlgth = fabs(hlgth);

	/*           compute the 21-point kronrod approximation to */
	/*           the integral, and estimate the absolute error. */

	resg = 0.;
	fc = f(centr);
	resk = wgk[10] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 5; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		d__1 = centr - absc;
		fval1 = f(d__1);
		d__1 = centr + absc;
		fval2 = f(d__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 5; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		d__1 = centr - absc;
		fval1 = f(d__1);
		d__1 = centr + absc;
		fval2 = f(d__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5;
	*resasc = wgk[10] * (d__1 = fc - reskh, fabs(d__1));
	for (j = 1; j <= 10; ++j) {
		*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, fabs(d__1)) + (
											 d__2 = fv2[j - 1] - reskh, fabs(d__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (d__1 = (resk - resg) * hlgth, fabs(d__1));
	if (*resasc != 0. && *abserr != 0.) {
		/* Computing MIN */
		d__3 = *abserr * 200. / *resasc;
		d__1 = 1., d__2 = pow(d__3, c_b270);
		*abserr = *resasc * min(d__1,d__2);
	}
	if (*resabs > uflow / (epmach * 50.)) {
		/* Computing MAX */
		d__1 = epmach * 50. * *resabs;
		*abserr = max(d__1,*abserr);
	}
	return;
} /* dqk21_ */

void dqk31_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc)
{
	/* Initialized data */

	double wg[8] = { .030753241996117268354628393577204,
									.070366047488108124709267416450667,
									.107159220467171935011869546685869,
									.139570677926154314447804794511028,
									.166269205816993933553200860481209,
									.186161000015562211026800561866423,
									.198431485327111576456118326443839,
									.202578241925561272880620199967519 };
	double xgk[16] = { .998002298693397060285172840152271,
									  .987992518020485428489565718586613,
									  .967739075679139134257347978784337,
									  .937273392400705904307758947710209,
									  .897264532344081900882509656454496,
									  .848206583410427216200648320774217,
									  .790418501442465932967649294817947,
									  .724417731360170047416186054613938,
									  .650996741297416970533735895313275,
									  .570972172608538847537226737253911,
									  .485081863640239680693655740232351,
									  .394151347077563369897207370981045,
									  .299180007153168812166780024266389,
									  .201194093997434522300628303394596,
									  .101142066918717499027074231447392,0. };
	double wgk[16] = { .005377479872923348987792051430128,
									  .015007947329316122538374763075807,
									  .025460847326715320186874001019653,
									  .03534636079137584622203794847836,
									  .04458975132476487660822729937328,
									  .05348152469092808726534314723943,
									  .062009567800670640285139230960803,
									  .069854121318728258709520077099147,
									  .076849680757720378894432777482659,
									  .083080502823133021038289247286104,
									  .088564443056211770647275443693774,
									  .093126598170825321225486872747346,
									  .096642726983623678505179907627589,
									  .099173598721791959332393173484603,
									  .10076984552387559504494666261757,
									  .101330007014791549017374792767493 };

	/* System generated locals */
	double d__1, d__2, d__3;

	/* Local variables */
	long j;
	double fc, fv1[15], fv2[15];
	long jtw;
	double absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	double hlgth, centr, reskh, uflow;
	double epmach, dhlgth;

	/* ***begin prologue  dqk31 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  31-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b) with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           double precision version */

	/*           parameters */
	/*            on entry */
	/*              f      - double precision */
	/*                       function subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the calling program. */

	/*              a      - double precision */
	/*                       lower limit of integration */

	/*              b      - double precision */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - double precision */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 31-polong */
	/*                       gauss-kronrod rule (resk), obtained by optimal */
	/*                       addition of abscissae to the 15-point gauss */
	/*                       rule (resg). */

	/*              abserr - double precison */
	/*                       estimate of the modulus of the modulus, */
	/*                       which should not exceed fabs(i-result) */

	/*              resabs - double precision */
	/*                       approximation to the integral j */

	/*              resasc - double precision */
	/*                       approximation to the integral of fabs(f-i/(b-a)) */
	/*                       over (a,b) */

	/* ***references  (none) */
	/* ***routines called  d1mach */
	/* ***end prologue  dqk31 */


	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 31-point kronrod rule */
	/*                    xgk(2), xgk(4), ...  abscissae of the 15-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 15-point gauss rule */

	/*           wgk    - weights of the 31-point kronrod rule */

	/*           wg     - weights of the 15-point gauss rule */


	/* gauss quadrature weights and kronron quadrature abscissae and weights */
	/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
	/* bell labs, nov. 1981. */





	/*           list of major variables */
	/*           ----------------------- */
	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 15-point gauss formula */
	/*           resk   - result of the 31-point kronrod formula */
	/*           reskh  - approximation to the mean value of f over (a,b), */
	/*                    i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */
	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */
	/* ***first executable statement  dqk31 */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);

	centr = (*a + *b) * .5;
	hlgth = (*b - *a) * .5;
	dhlgth = fabs(hlgth);

	/*           compute the 31-point kronrod approximation to */
	/*           the integral, and estimate the absolute error. */

	fc = f(centr);
	resg = wg[7] * fc;
	resk = wgk[15] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 7; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		d__1 = centr - absc;
		fval1 = f(d__1);
		d__1 = centr + absc;
		fval2 = f(d__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 8; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		d__1 = centr - absc;
		fval1 = f(d__1);
		d__1 = centr + absc;
		fval2 = f(d__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5;
	*resasc = wgk[15] * (d__1 = fc - reskh, fabs(d__1));
	for (j = 1; j <= 15; ++j) {
		*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, fabs(d__1)) + (
											 d__2 = fv2[j - 1] - reskh, fabs(d__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (d__1 = (resk - resg) * hlgth, fabs(d__1));
	if (*resasc != 0. && *abserr != 0.) {
		/* Computing MIN */
		d__3 = *abserr * 200. / *resasc;
		d__1 = 1., d__2 = pow(d__3, c_b270);
		*abserr = *resasc * min(d__1,d__2);
	}
	if (*resabs > uflow / (epmach * 50.)) {
		/* Computing MAX */
		d__1 = epmach * 50. * *resabs;
		*abserr = max(d__1,*abserr);
	}
	return;
} /* dqk31_ */

void dqk41_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc)
{
	/* Initialized data */

	double wg[10] = { .017614007139152118311861962351853,
									 .040601429800386941331039952274932,
									 .062672048334109063569506535187042,
									 .083276741576704748724758143222046,
									 .10193011981724043503675013548035,
									 .118194531961518417312377377711382,
									 .131688638449176626898494499748163,
									 .142096109318382051329298325067165,
									 .149172986472603746787828737001969,
									 .152753387130725850698084331955098 };
	double xgk[21] = { .998859031588277663838315576545863,
									  .99312859918509492478612238847132,
									  .981507877450250259193342994720217,
									  .963971927277913791267666131197277,
									  .940822633831754753519982722212443,
									  .912234428251325905867752441203298,
									  .878276811252281976077442995113078,
									  .839116971822218823394529061701521,
									  .795041428837551198350638833272788,
									  .746331906460150792614305070355642,
									  .693237656334751384805490711845932,
									  .636053680726515025452836696226286,
									  .575140446819710315342946036586425,
									  .510867001950827098004364050955251,
									  .44359317523872510319999221349264,
									  .373706088715419560672548177024927,
									  .301627868114913004320555356858592,
									  .227785851141645078080496195368575,
									  .152605465240922675505220241022678,
									  .076526521133497333754640409398838,0. };
	double wgk[21] = { .003073583718520531501218293246031,
									  .008600269855642942198661787950102,
									  .014626169256971252983787960308868,
									  .020388373461266523598010231432755,
									  .025882133604951158834505067096153,
									  .031287306777032798958543119323801,
									  .036600169758200798030557240707211,
									  .041668873327973686263788305936895,
									  .046434821867497674720231880926108,
									  .050944573923728691932707670050345,
									  .055195105348285994744832372419777,
									  .059111400880639572374967220648594,
									  .062653237554781168025870122174255,
									  .065834597133618422111563556969398,
									  .068648672928521619345623411885368,
									  .07105442355344406830579036172321,
									  .073030690332786667495189417658913,
									  .074582875400499188986581418362488,
									  .075704497684556674659542775376617,
									  .076377867672080736705502835038061,
									  .076600711917999656445049901530102 };

	/* System generated locals */
	double d__1, d__2, d__3;

	/* Local variables */
	long j;
	double fc, fv1[20], fv2[20];
	long jtw;
	double absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	double hlgth, centr, reskh, uflow;
	double epmach, dhlgth;

	/* ***begin prologue  dqk41 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  41-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b), with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           double precision version */

	/*           parameters */
	/*            on entry */
	/*              f      - double precision */
	/*                       function subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the calling program. */

	/*              a      - double precision */
	/*                       lower limit of integration */

	/*              b      - double precision */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - double precision */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 41-polong */
	/*                       gauss-kronrod rule (resk) obtained by optimal */
	/*                       addition of abscissae to the 20-point gauss */
	/*                       rule (resg). */

	/*              abserr - double precision */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should not exceed fabs(i-result) */

	/*              resabs - double precision */
	/*                       approximation to the integral j */

	/*              resasc - double precision */
	/*                       approximation to the integal of fabs(f-i/(b-a)) */
	/*                       over (a,b) */

	/* ***references  (none) */
	/* ***routines called  d1mach */
	/* ***end prologue  dqk41 */



	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 41-point gauss-kronrod rule */
	/*                    xgk(2), xgk(4), ...  abscissae of the 20-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 20-point gauss rule */

	/*           wgk    - weights of the 41-point gauss-kronrod rule */

	/*           wg     - weights of the 20-point gauss rule */


	/* gauss quadrature weights and kronron quadrature abscissae and weights */
	/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
	/* bell labs, nov. 1981. */





	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 20-point gauss formula */
	/*           resk   - result of the 41-point kronrod formula */
	/*           reskh  - approximation to mean value of f over (a,b), i.e. */
	/*                    to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  dqk41 */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);

	centr = (*a + *b) * .5;
	hlgth = (*b - *a) * .5;
	dhlgth = fabs(hlgth);

	/*           compute the 41-point gauss-kronrod approximation to */
	/*           the integral, and estimate the absolute error. */

	resg = 0.;
	fc = f(centr);
	resk = wgk[20] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 10; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		d__1 = centr - absc;
		fval1 = f(d__1);
		d__1 = centr + absc;
		fval2 = f(d__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 10; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		d__1 = centr - absc;
		fval1 = f(d__1);
		d__1 = centr + absc;
		fval2 = f(d__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5;
	*resasc = wgk[20] * (d__1 = fc - reskh, fabs(d__1));
	for (j = 1; j <= 20; ++j) {
		*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, fabs(d__1)) + (
											 d__2 = fv2[j - 1] - reskh, fabs(d__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (d__1 = (resk - resg) * hlgth, fabs(d__1));
	if (*resasc != 0. && *abserr != 0.) {
		/* Computing MIN */
		d__3 = *abserr * 200. / *resasc;
		d__1 = 1., d__2 = pow(d__3, c_b270);
		*abserr = *resasc * min(d__1,d__2);
	}
	if (*resabs > uflow / (epmach * 50.)) {
		/* Computing MAX */
		d__1 = epmach * 50. * *resabs;
		*abserr = max(d__1,*abserr);
	}
	return;
} /* dqk41_ */

void dqk51_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc)
{
	/* Initialized data */

	double wg[13] = { .011393798501026287947902964113235,
									 .026354986615032137261901815295299,
									 .040939156701306312655623487711646,
									 .054904695975835191925936891540473,
									 .068038333812356917207187185656708,
									 .080140700335001018013234959669111,
									 .091028261982963649811497220702892,
									 .100535949067050644202206890392686,
									 .108519624474263653116093957050117,
									 .114858259145711648339325545869556,
									 .119455763535784772228178126512901,
									 .122242442990310041688959518945852,
									 .12317605372671545120390287307905 };
	double xgk[26] = { .999262104992609834193457486540341,
									  .995556969790498097908784946893902,
									  .988035794534077247637331014577406,
									  .976663921459517511498315386479594,
									  .961614986425842512418130033660167,
									  .942974571228974339414011169658471,
									  .920747115281701561746346084546331,
									  .894991997878275368851042006782805,
									  .86584706529327559544899696958834,
									  .83344262876083400142102110869357,
									  .797873797998500059410410904994307,
									  .759259263037357630577282865204361,
									  .717766406813084388186654079773298,
									  .673566368473468364485120633247622,
									  .626810099010317412788122681624518,
									  .577662930241222967723689841612654,
									  .52632528433471918259962377815801,
									  .473002731445714960522182115009192,
									  .417885382193037748851814394594572,
									  .361172305809387837735821730127641,
									  .303089538931107830167478909980339,
									  .243866883720988432045190362797452,
									  .183718939421048892015969888759528,
									  .122864692610710396387359818808037,
									  .061544483005685078886546392366797,0. };
	double wgk[26] = { .001987383892330315926507851882843,
									  .005561932135356713758040236901066,
									  .009473973386174151607207710523655,
									  .013236229195571674813656405846976,
									  .016847817709128298231516667536336,
									  .020435371145882835456568292235939,
									  .024009945606953216220092489164881,
									  .027475317587851737802948455517811,
									  .030792300167387488891109020215229,
									  .034002130274329337836748795229551,
									  .03711627148341554356033062536762,
									  .040083825504032382074839284467076,
									  .042872845020170049476895792439495,
									  .04550291304992178890987058475266,
									  .047982537138836713906392255756915,
									  .05027767908071567196332525943344,
									  .052362885806407475864366712137873,
									  .054251129888545490144543370459876,
									  .055950811220412317308240686382747,
									  .057437116361567832853582693939506,
									  .058689680022394207961974175856788,
									  .059720340324174059979099291932562,
									  .060539455376045862945360267517565,
									  .061128509717053048305859030416293,
									  .061471189871425316661544131965264,
									  .061580818067832935078759824240066 };

	/* System generated locals */
	double d__1, d__2, d__3;

	/* Local variables */
	long j;
	double fc, fv1[25], fv2[25];
	long jtw;
	double absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	double hlgth, centr, reskh, uflow;
	double epmach, dhlgth;

	/* ***begin prologue  dqk51 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  51-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b) with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           double precision version */

	/*           parameters */
	/*            on entry */
	/*              f      - double precision */
	/*                       function subroutine defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the calling program. */

	/*              a      - double precision */
	/*                       lower limit of integration */

	/*              b      - double precision */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - double precision */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 51-polong */
	/*                       kronrod rule (resk) obtained by optimal addition */
	/*                       of abscissae to the 25-point gauss rule (resg). */

	/*              abserr - double precision */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should not exceed fabs(i-result) */

	/*              resabs - double precision */
	/*                       approximation to the integral j */

	/*              resasc - double precision */
	/*                       approximation to the integral of fabs(f-i/(b-a)) */
	/*                       over (a,b) */

	/* ***references  (none) */
	/* ***routines called  d1mach */
	/* ***end prologue  dqk51 */



	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 51-point kronrod rule */
	/*                    xgk(2), xgk(4), ...  abscissae of the 25-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 25-point gauss rule */

	/*           wgk    - weights of the 51-point kronrod rule */

	/*           wg     - weights of the 25-point gauss rule */


	/* gauss quadrature weights and kronron quadrature abscissae and weights */
	/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
	/* bell labs, nov. 1981. */



	/*       note: wgk (26) was calculated from the values of wgk(1..25) */


	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 25-point gauss formula */
	/*           resk   - result of the 51-point kronrod formula */
	/*           reskh  - approximation to the mean value of f over (a,b), */
	/*                    i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  dqk51 */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);

	centr = (*a + *b) * .5;
	hlgth = (*b - *a) * .5;
	dhlgth = fabs(hlgth);

	/*           compute the 51-point kronrod approximation to */
	/*           the integral, and estimate the absolute error. */

	fc = f(centr);
	resg = wg[12] * fc;
	resk = wgk[25] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 12; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		d__1 = centr - absc;
		fval1 = f(d__1);
		d__1 = centr + absc;
		fval2 = f(d__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 13; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		d__1 = centr - absc;
		fval1 = f(d__1);
		d__1 = centr + absc;
		fval2 = f(d__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5;
	*resasc = wgk[25] * (d__1 = fc - reskh, fabs(d__1));
	for (j = 1; j <= 25; ++j) {
		*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, fabs(d__1)) + (
											 d__2 = fv2[j - 1] - reskh, fabs(d__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (d__1 = (resk - resg) * hlgth, fabs(d__1));
	if (*resasc != 0. && *abserr != 0.) {
		/* Computing MIN */
		d__3 = *abserr * 200. / *resasc;
		d__1 = 1., d__2 = pow(d__3, c_b270);
		*abserr = *resasc * min(d__1,d__2);
	}
	if (*resabs > uflow / (epmach * 50.)) {
		/* Computing MAX */
		d__1 = epmach * 50. * *resabs;
		*abserr = max(d__1,*abserr);
	}
	return;
} /* dqk51_ */

void dqk61_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc)
{
	/* Initialized data */

	double wg[15] = { .007968192496166605615465883474674,
									 .018466468311090959142302131912047,
									 .028784707883323369349719179611292,
									 .038799192569627049596801936446348,
									 .048402672830594052902938140422808,
									 .057493156217619066481721689402056,
									 .065974229882180495128128515115962,
									 .073755974737705206268243850022191,
									 .08075589522942021535469493846053,
									 .086899787201082979802387530715126,
									 .092122522237786128717632707087619,
									 .09636873717464425963946862635181,
									 .099593420586795267062780282103569,
									 .101762389748405504596428952168554,
									 .102852652893558840341285636705415 };
	double xgk[31] = { .999484410050490637571325895705811,
									  .996893484074649540271630050918695,
									  .991630996870404594858628366109486,
									  .983668123279747209970032581605663,
									  .973116322501126268374693868423707,
									  .960021864968307512216871025581798,
									  .944374444748559979415831324037439,
									  .926200047429274325879324277080474,
									  .905573307699907798546522558925958,
									  .882560535792052681543116462530226,
									  .857205233546061098958658510658944,
									  .829565762382768397442898119732502,
									  .799727835821839083013668942322683,
									  .767777432104826194917977340974503,
									  .733790062453226804726171131369528,
									  .69785049479331579693229238802664,
									  .660061064126626961370053668149271,
									  .620526182989242861140477556431189,
									  .57934523582636169175602493217254,
									  .536624148142019899264169793311073,
									  .492480467861778574993693061207709,
									  .447033769538089176780609900322854,
									  .400401254830394392535476211542661,
									  .352704725530878113471037207089374,
									  .304073202273625077372677107199257,
									  .254636926167889846439805129817805,
									  .204525116682309891438957671002025,
									  .153869913608583546963794672743256,
									  .102806937966737030147096751318001,
									  .051471842555317695833025213166723,0. };
	double wgk[31] = { .00138901369867700762455159122676,
									  .003890461127099884051267201844516,
									  .00663070391593129217331982636975,
									  .009273279659517763428441146892024,
									  .011823015253496341742232898853251,
									  .01436972950704580481245143244358,
									  .016920889189053272627572289420322,
									  .019414141193942381173408951050128,
									  .021828035821609192297167485738339,
									  .024191162078080601365686370725232,
									  .026509954882333101610601709335075,
									  .028754048765041292843978785354334,
									  .030907257562387762472884252943092,
									  .032981447057483726031814191016854,
									  .034979338028060024137499670731468,
									  .036882364651821229223911065617136,
									  .038678945624727592950348651532281,
									  .040374538951535959111995279752468,
									  .04196981021516424614714754128597,
									  .043452539701356069316831728117073,
									  .044814800133162663192355551616723,
									  .046059238271006988116271735559374,
									  .047185546569299153945261478181099,
									  .048185861757087129140779492298305,
									  .049055434555029778887528165367238,
									  .049795683427074206357811569379942,
									  .050405921402782346840893085653585,
									  .050881795898749606492297473049805,
									  .051221547849258772170656282604944,
									  .051426128537459025933862879215781,
									  .051494729429451567558340433647099 };

	/* System generated locals */
	double d__1, d__2, d__3;

	/* Local variables */
	long j;
	double fc, fv1[30], fv2[30];
	long jtw;
	double resg, resk, fsum, fval1, fval2;
	long jtwm1;
	double dabsc, hlgth, centr, reskh, uflow;
	double epmach, dhlgth;

	/* ***begin prologue  dqk61 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  61-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b) with error */
	/*                           estimate */
	/*                       j = integral of dfabs(f) over (a,b) */
	/* ***description */

	/*        integration rule */
	/*        standard fortran subroutine */
	/*        double precision version */


	/*        parameters */
	/*         on entry */
	/*           f      - double precision */
	/*                    function subprogram defining the integrand */
	/*                    function f(x). the actual name for f needs to be */
	/*                    declared e x t e r n a l in the calling program. */

	/*           a      - double precision */
	/*                    lower limit of integration */

	/*           b      - double precision */
	/*                    upper limit of integration */

	/*         on return */
	/*           result - double precision */
	/*                    approximation to the integral i */
	/*                    result is computed by applying the 61-polong */
	/*                    kronrod rule (resk) obtained by optimal addition of */
	/*                    abscissae to the 30-point gauss rule (resg). */

	/*           abserr - double precision */
	/*                    estimate of the modulus of the absolute error, */
	/*                    which should equal or exceed dfabs(i-result) */

	/*           resabs - double precision */
	/*                    approximation to the integral j */

	/*           resasc - double precision */
	/*                    approximation to the integral of dfabs(f-i/(b-a)) */


	/* ***references  (none) */
	/* ***routines called  d1mach */
	/* ***end prologue  dqk61 */



	/*           the abscissae and weights are given for the */
	/*           interval (-1,1). because of symmetry only the positive */
	/*           abscissae and their corresponding weights are given. */

	/*           xgk   - abscissae of the 61-point kronrod rule */
	/*                   xgk(2), xgk(4)  ... abscissae of the 30-polong */
	/*                   gauss rule */
	/*                   xgk(1), xgk(3)  ... optimally added abscissae */
	/*                   to the 30-point gauss rule */

	/*           wgk   - weights of the 61-point kronrod rule */

	/*           wg    - weigths of the 30-point gauss rule */


	/* gauss quadrature weights and kronron quadrature abscissae and weights */
	/* as evaluated with 80 decimal digit arithmetic by l. w. fullerton, */
	/* bell labs, nov. 1981. */




	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           dabsc  - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 30-point gauss rule */
	/*           resk   - result of the 61-point kronrod rule */
	/*           reskh  - approximation to the mean value of f */
	/*                    over (a,b), i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	epmach = d1mach(c__4);
	uflow = d1mach(c__1);

	centr = (*b + *a) * .5;
	hlgth = (*b - *a) * .5;
	dhlgth = fabs(hlgth);

	/*           compute the 61-point kronrod approximation to the */
	/*           integral, and estimate the absolute error. */

	/* ***first executable statement  dqk61 */
	resg = 0.;
	fc = f(centr);
	resk = wgk[30] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 15; ++j) {
		jtw = j << 1;
		dabsc = hlgth * xgk[jtw - 1];
		d__1 = centr - dabsc;
		fval1 = f(d__1);
		d__1 = centr + dabsc;
		fval2 = f(d__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 15; ++j) {
		jtwm1 = (j << 1) - 1;
		dabsc = hlgth * xgk[jtwm1 - 1];
		d__1 = centr - dabsc;
		fval1 = f(d__1);
		d__1 = centr + dabsc;
		fval2 = f(d__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5;
	*resasc = wgk[30] * (d__1 = fc - reskh, fabs(d__1));
	for (j = 1; j <= 30; ++j) {
		*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, fabs(d__1)) + (
											 d__2 = fv2[j - 1] - reskh, fabs(d__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (d__1 = (resk - resg) * hlgth, fabs(d__1));
	if (*resasc != 0. && *abserr != 0.) {
		/* Computing MIN */
		d__3 = *abserr * 200. / *resasc;
		d__1 = 1., d__2 = pow(d__3, c_b270);
		*abserr = *resasc * min(d__1,d__2);
	}
	if (*resabs > uflow / (epmach * 50.)) {
		/* Computing MAX */
		d__1 = epmach * 50. * *resabs;
		*abserr = max(d__1,*abserr);
	}
	return;
} /* dqk61_ */

void dqmomo_(const double *alfa, const double *beta, double *
				 ri, double *rj, double *rg, double *rh, const long *integr)
{
	/* Local variables */
	long i__;
	double an;
	long im1;
	double anm1, ralf, rbet, alfp1, alfp2, betp1, betp2;

	/* ***begin prologue  dqmomo */
	/* ***date written   820101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1,c3a2 */
	/* ***keywords  modified chebyshev moments */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  this routine computes modified chebsyshev moments. the k-th */
	/*            modified chebyshev moment is defined as the integral over */
	/*            (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev */
	/*            polynomial of degree k. */
	/* ***description */

	/*        modified chebyshev moments */
	/*        standard fortran subroutine */
	/*        double precision version */

	/*        parameters */
	/*           alfa   - double precision */
	/*                    parameter in the weight function w(x), alfa.gt.(-1) */

	/*           beta   - double precision */
	/*                    parameter in the weight function w(x), beta.gt.(-1) */

	/*           ri     - double precision */
	/*                    vector of dimension 25 */
	/*                    ri(k) is the integral over (-1,1) of */
	/*                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25. */

	/*           rj     - double precision */
	/*                    vector of dimension 25 */
	/*                    rj(k) is the integral over (-1,1) of */
	/*                    (1-x)**beta*t(k-1,x), k = 1, ..., 25. */

	/*           rg     - double precision */
	/*                    vector of dimension 25 */
	/*                    rg(k) is the integral over (-1,1) of */
	/*                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25. */

	/*           rh     - double precision */
	/*                    vector of dimension 25 */
	/*                    rh(k) is the integral over (-1,1) of */
	/*                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25. */

	/*           integr - long */
	/*                    input parameter indicating the modified */
	/*                    moments to be computed */
	/*                    integr = 1 compute ri, rj */
	/*                           = 2 compute ri, rj, rg */
	/*                           = 3 compute ri, rj, rh */
	/*                           = 4 compute ri, rj, rg, rh */

	/* ***references  (none) */
	/* ***routines called  (none) */
	/* ***end prologue  dqmomo */




	/* ***first executable statement  dqmomo */
	/* Parameter adjustments */
	--rh;
	--rg;
	--rj;
	--ri;

	/* Function Body */
	alfp1 = *alfa + 1.;
	betp1 = *beta + 1.;
	alfp2 = *alfa + 2.;
	betp2 = *beta + 2.;
	ralf = pow(c_b320, alfp1);
	rbet = pow(c_b320, betp1);

	/*           compute ri, rj using a forward recurrence relation. */

	ri[1] = ralf / alfp1;
	rj[1] = rbet / betp1;
	ri[2] = ri[1] * *alfa / alfp2;
	rj[2] = rj[1] * *beta / betp2;
	an = 2.;
	anm1 = 1.;
	for (i__ = 3; i__ <= 25; ++i__) {
		ri[i__] = -(ralf + an * (an - alfp2) * ri[i__ - 1]) / (anm1 * (an + 
																							alfp1));
		rj[i__] = -(rbet + an * (an - betp2) * rj[i__ - 1]) / (anm1 * (an + 
																							betp1));
		anm1 = an;
		an += 1.;
		/* L20: */
	}
	if (*integr == 1) {
		goto L70;
	}
	if (*integr == 3) {
		goto L40;
	}

	/*           compute rg using a forward recurrence relation. */

	rg[1] = -ri[1] / alfp1;
	rg[2] = -(ralf + ralf) / (alfp2 * alfp2) - rg[1];
	an = 2.;
	anm1 = 1.;
	im1 = 2;
	for (i__ = 3; i__ <= 25; ++i__) {
		rg[i__] = -(an * (an - alfp2) * rg[im1] - an * ri[im1] + anm1 * ri[
							i__]) / (anm1 * (an + alfp1));
		anm1 = an;
		an += 1.;
		im1 = i__;
		/* L30: */
	}
	if (*integr == 2) {
		goto L70;
	}

	/*           compute rh using a forward recurrence relation. */

L40:
	rh[1] = -rj[1] / betp1;
	rh[2] = -(rbet + rbet) / (betp2 * betp2) - rh[1];
	an = 2.;
	anm1 = 1.;
	im1 = 2;
	for (i__ = 3; i__ <= 25; ++i__) {
		rh[i__] = -(an * (an - betp2) * rh[im1] - an * rj[im1] + anm1 * rj[
							i__]) / (anm1 * (an + betp1));
		anm1 = an;
		an += 1.;
		im1 = i__;
		/* L50: */
	}
	for (i__ = 2; i__ <= 25; i__ += 2) {
		rh[i__] = -rh[i__];
		/* L60: */
	}
L70:
	for (i__ = 2; i__ <= 25; i__ += 2) {
		rj[i__] = -rj[i__];
		/* L80: */
	}
	/* L90: */
	return;
} /* dqmomo_ */

void dqng_(const D_fp& f, const double *a, const double *b, const double *epsabs, 
			  const double *epsrel, double *result, double *abserr, 
			  long *neval, long *ier)
{
	/* Initialized data */

	double x1[5] = { .973906528517171720077964012084452,
									.865063366688984510732096688423493,
									.679409568299024406234327365114874,
									.433395394129247190799265943165784,
									.14887433898163121088482600112972 };
	double w10[5] = { .066671344308688137593568809893332,
									 .149451349150580593145776339657697,
									 .219086362515982043995534934228163,
									 .269266719309996355091226921569469,
									 .295524224714752870173892994651338 };
	double x2[5] = { .995657163025808080735527280689003,
									.930157491355708226001207180059508,
									.780817726586416897063717578345042,
									.562757134668604683339000099272694,
									.294392862701460198131126603103866 };
	double w21a[5] = { .03255816230796472747881897245939,
									  .07503967481091995276704314091619,
									  .109387158802297641899210590325805,
									  .134709217311473325928054001771707,
									  .147739104901338491374841515972068 };
	double w21b[6] = { .011694638867371874278064396062192,
									  .05475589657435199603138130024458,
									  .093125454583697605535065465083366,
									  .123491976262065851077958109831074,
									  .142775938577060080797094273138717,
									  .149445554002916905664936468389821 };
	double x3[11] = { .999333360901932081394099323919911,
									 .987433402908088869795961478381209,
									 .954807934814266299257919200290473,
									 .900148695748328293625099494069092,
									 .82519831498311415084706673258852,
									 .732148388989304982612354848755461,
									 .622847970537725238641159120344323,
									 .499479574071056499952214885499755,
									 .364901661346580768043989548502644,
									 .222254919776601296498260928066212,
									 .074650617461383322043914435796506 };
	double w43a[10] = { .016296734289666564924281974617663,
										.037522876120869501461613795898115,
										.054694902058255442147212685465005,
										.067355414609478086075553166302174,
										.073870199632393953432140695251367,
										.005768556059769796184184327908655,
										.027371890593248842081276069289151,
										.046560826910428830743339154433824,
										.061744995201442564496240336030883,
										.071387267268693397768559114425516 };
	double w43b[12] = { .001844477640212414100389106552965,
										.010798689585891651740465406741293,
										.021895363867795428102523123075149,
										.032597463975345689443882222526137,
										.042163137935191811847627924327955,
										.050741939600184577780189020092084,
										.058379395542619248375475369330206,
										.064746404951445885544689259517511,
										.069566197912356484528633315038405,
										.072824441471833208150939535192842,
										.074507751014175118273571813842889,
										.074722147517403005594425168280423 };
	double x4[22] = { .999902977262729234490529830591582,
									 .99798989598667874542749632236596,
									 .992175497860687222808523352251425,
									 .981358163572712773571916941623894,
									 .965057623858384619128284110607926,
									 .943167613133670596816416634507426,
									 .91580641468550720959182643072005,
									 .883221657771316501372117548744163,
									 .845710748462415666605902011504855,
									 .803557658035230982788739474980964,
									 .75700573068549555832894279343202,
									 .70627320978732181982409427474084,
									 .651589466501177922534422205016736,
									 .593223374057961088875273770349144,
									 .531493605970831932285268948562671,
									 .46676362304202284487196678165927,
									 .399424847859218804732101665817923,
									 .329874877106188288265053371824597,
									 .258503559202161551802280975429025,
									 .185695396568346652015917141167606,
									 .111842213179907468172398359241362,
									 .037352123394619870814998165437704 };
	double w87a[21] = { .00814837738414917290000287844819,
										.018761438201562822243935059003794,
										.027347451050052286161582829741283,
										.033677707311637930046581056957588,
										.036935099820427907614589586742499,
										.002884872430211530501334156248695,
										.013685946022712701888950035273128,
										.023280413502888311123409291030404,
										.030872497611713358675466394126442,
										.035693633639418770719351355457044,
										9.15283345202241360843392549948e-4,
										.005399280219300471367738743391053,
										.010947679601118931134327826856808,
										.01629873169678733526266570322328,
										.02108156888920383511243306018819,
										.02537096976925382724346799983171,
										.02918969775647575250144615408492,
										.032373202467202789685788194889595,
										.034783098950365142750781997949596,
										.036412220731351787562801163687577,
										.037253875503047708539592001191226 };
	double w87b[23] = { 2.74145563762072350016527092881e-4,
										.001807124155057942948341311753254,
										.00409686928275916486445807068348,
										.006758290051847378699816577897424,
										.009549957672201646536053581325377,
										.01232944765224485369462663996378,
										.015010447346388952376697286041943,
										.0175489679862431910996653529259,
										.019938037786440888202278192730714,
										.022194935961012286796332102959499,
										.024339147126000805470360647041454,
										.026374505414839207241503786552615,
										.02828691078877120065996800298796,
										.030052581128092695322521110347341,
										.031646751371439929404586051078883,
										.033050413419978503290785944862689,
										.034255099704226061787082821046821,
										.035262412660156681033782717998428,
										.036076989622888701185500318003895,
										.036698604498456094498018047441094,
										.037120549269832576114119958413599,
										.037334228751935040321235449094698,
										.037361073762679023410321241766599 };

	/* System generated locals */
	double d__1, d__2, d__3, d__4;

	/* Local variables */
	long k, l;
	double fv1[5], fv2[5], fv3[5], fv4[5];
	long ipx=0.;
	double absc, fval, res10, res21=0., res43=0., res87, fval1, fval2, 
		hlgth, centr, reskh, uflow;
	double epmach, dhlgth, resabs=0., resasc=0., fcentr, savfun[21];

	/* ***begin prologue  dqng */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  810101   (yymmdd) */
	/* ***category no.  h2a1a1 */
	/* ***keywords  automatic integrator, smooth integrand, */
	/*             non-adaptive, gauss-kronrod(patterson) */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl math & progr. div. - k.u.leuven */
	/*           kahaner,david,nbs - modified (2/82) */
	/* ***purpose  the routine calculates an approximation result to a */
	/*            given definite integral i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/* non-adaptive integration */
	/* standard fortran subroutine */
	/* double precision version */

	/*           f      - double precision */
	/*                    function subprogram defining the integrand function */
	/*                    f(x). the actual name for f needs to be declared */
	/*                    e x t e r n a l in the driver program. */

	/*           a      - double precision */
	/*                    lower limit of integration */

	/*           b      - double precision */
	/*                    upper limit of integration */

	/*           epsabs - double precision */
	/*                    absolute accuracy requested */
	/*           epsrel - double precision */
	/*                    relative accuracy requested */
	/*                    if  epsabs.le.0 */
	/*                    and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                    the routine will end with ier = 6. */

	/*         on return */
	/*           result - double precision */
	/*                    approximation to the integral i */
	/*                    result is obtained by applying the 21-polong */
	/*                    gauss-kronrod rule (res21) obtained by optimal */
	/*                    addition of abscissae to the 10-point gauss rule */
	/*                    (res10), or by applying the 43-point rule (res43) */
	/*                    obtained by optimal addition of abscissae to the */
	/*                    21-point gauss-kronrod rule, or by applying the */
	/*                    87-point rule (res87) obtained by optimal addition */
	/*                    of abscissae to the 43-point rule. */

	/*           abserr - double precision */
	/*                    estimate of the modulus of the absolute error, */
	/*                    which should equal or exceed fabs(i-result) */

	/*           neval  - long */
	/*                    number of integrand evaluations */

	/*           ier    - ier = 0 normal and reliable termination of the */
	/*                            routine. it is assumed that the requested */
	/*                            accuracy has been achieved. */
	/*                    ier.gt.0 abnormal termination of the routine. it is */
	/*                            assumed that the requested accuracy has */
	/*                            not been achieved. */
	/*           error messages */
	/*                    ier = 1 the maximum number of steps has been */
	/*                            executed. the integral is probably too */
	/*                            difficult to be calculated by dqng. */
	/*                        = 6 the input is invalid, because */
	/*                            epsabs.le.0 and */
	/*                            epsrel.lt.max(50*rel.mach.acc.,0.5d-28). */
	/*                            result, abserr and neval are set to zero. */

	/* ***references  (none) */
	/* ***routines called  d1mach,xerror */
	/* ***end prologue  dqng */



	/*           the following data statements contain the */
	/*           abscissae and weights of the integration rules used. */

	/*           x1      abscissae common to the 10-, 21-, 43- and 87- */
	/*                   point rule */
	/*           x2      abscissae common to the 21-, 43- and 87-point rule */
	/*           x3      abscissae common to the 43- and 87-point rule */
	/*           x4      abscissae of the 87-point rule */
	/*           w10     weights of the 10-point formula */
	/*           w21a    weights of the 21-point formula for abscissae x1 */
	/*           w21b    weights of the 21-point formula for abscissae x2 */
	/*           w43a    weights of the 43-point formula for abscissae x1, x3 */
	/*           w43b    weights of the 43-point formula for abscissae x3 */
	/*           w87a    weights of the 87-point formula for abscissae x1, */
	/*                   x2, x3 */
	/*           w87b    weights of the 87-point formula for abscissae x4 */


	/* gauss-kronrod-patterson quadrature coefficients for use in */
	/* quadpack routine qng.  these coefficients were calculated with */
	/* 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981. */





	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the integration interval */
	/*           hlgth  - half-length of the integration interval */
	/*           fcentr - function value at mid polong */
	/*           absc   - abscissa */
	/*           fval   - function value */
	/*           savfun - array of function values which have already been */
	/*                    computed */
	/*           res10  - 10-point gauss result */
	/*           res21  - 21-point kronrod result */
	/*           res43  - 43-point result */
	/*           res87  - 87-point result */
	/*           resabs - approximation to the integral of fabs(f) */
	/*           resasc - approximation to the integral of fabs(f-i/(b-a)) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  dqng */
	epmach = d1mach(c__4);
	uflow = d1mach(c__1);

	/*           test on validity of parameters */
	/*           ------------------------------ */

	*result = 0.;
	*abserr = 0.;
	*neval = 0;
	*ier = 6;
	/* Computing MAX */
	d__1 = epmach * 50.;
	if (*epsabs <= 0. && *epsrel < max(d__1,5e-29)) {
		goto L80;
	}
	hlgth = (*b - *a) * .5;
	dhlgth = fabs(hlgth);
	centr = (*b + *a) * .5;
	fcentr = f(centr);
	*neval = 21;
	*ier = 1;

	/*          compute the integral using the 10- and 21-point formula. */

	for (l = 1; l <= 3; ++l) {
		switch (l) {
		case 1:  goto L5;
		case 2:  goto L25;
		case 3:  goto L45;
		}
	L5:
		res10 = 0.;
		res21 = w21b[5] * fcentr;
		resabs = w21b[5] * fabs(fcentr);
		for (k = 1; k <= 5; ++k) {
			absc = hlgth * x1[k - 1];
			d__1 = centr + absc;
			fval1 = f(d__1);
			d__1 = centr - absc;
			fval2 = f(d__1);
			fval = fval1 + fval2;
			res10 += w10[k - 1] * fval;
			res21 += w21a[k - 1] * fval;
			resabs += w21a[k - 1] * (fabs(fval1) + fabs(fval2));
			savfun[k - 1] = fval;
			fv1[k - 1] = fval1;
			fv2[k - 1] = fval2;
			/* L10: */
		}
		ipx = 5;
		for (k = 1; k <= 5; ++k) {
			++ipx;
			absc = hlgth * x2[k - 1];
			d__1 = centr + absc;
			fval1 = f(d__1);
			d__1 = centr - absc;
			fval2 = f(d__1);
			fval = fval1 + fval2;
			res21 += w21b[k - 1] * fval;
			resabs += w21b[k - 1] * (fabs(fval1) + fabs(fval2));
			savfun[ipx - 1] = fval;
			fv3[k - 1] = fval1;
			fv4[k - 1] = fval2;
			/* L15: */
		}

		/*          test for convergence. */

		*result = res21 * hlgth;
		resabs *= dhlgth;
		reskh = res21 * .5;
		resasc = w21b[5] * (d__1 = fcentr - reskh, fabs(d__1));
		for (k = 1; k <= 5; ++k) {
			resasc = resasc + w21a[k - 1] * ((d__1 = fv1[k - 1] - reskh, fabs(
															 d__1)) + (d__2 = fv2[k - 1] - reskh, fabs(d__2))) + w21b[k 
																																		 - 1] * ((d__3 = fv3[k - 1] - reskh, fabs(d__3)) + (d__4 = 
																																																			 fv4[k - 1] - reskh, fabs(d__4)));
			/* L20: */
		}
		*abserr = (d__1 = (res21 - res10) * hlgth, fabs(d__1));
		resasc *= dhlgth;
		goto L65;

		/*          compute the integral using the 43-point formula. */

	L25:
		res43 = w43b[11] * fcentr;
		*neval = 43;
		for (k = 1; k <= 10; ++k) {
			res43 += savfun[k - 1] * w43a[k - 1];
			/* L30: */
		}
		for (k = 1; k <= 11; ++k) {
			++ipx;
			absc = hlgth * x3[k - 1];
			d__1 = absc + centr;
			d__2 = centr - absc;
			fval = f(d__1) + f(d__2);
			res43 += fval * w43b[k - 1];
			savfun[ipx - 1] = fval;
			/* L40: */
		}

		/*          test for convergence. */

		*result = res43 * hlgth;
		*abserr = (d__1 = (res43 - res21) * hlgth, fabs(d__1));
		goto L65;

		/*          compute the integral using the 87-point formula. */

	L45:
		res87 = w87b[22] * fcentr;
		*neval = 87;
		for (k = 1; k <= 21; ++k) {
			res87 += savfun[k - 1] * w87a[k - 1];
			/* L50: */
		}
		for (k = 1; k <= 22; ++k) {
			absc = hlgth * x4[k - 1];
			d__1 = absc + centr;
			d__2 = centr - absc;
			res87 += w87b[k - 1] * (f(d__1) + f(d__2));
			/* L60: */
		}
		*result = res87 * hlgth;
		*abserr = (d__1 = (res87 - res43) * hlgth, fabs(d__1));
	L65:
		if (resasc != 0. && *abserr != 0.) {
			/* Computing MIN */
			d__3 = *abserr * 200. / resasc;
			d__1 = 1., d__2 = pow(d__3, c_b270);
			*abserr = resasc * min(d__1,d__2);
		}
		if (resabs > uflow / (epmach * 50.)) {
			/* Computing MAX */
			d__1 = epmach * 50. * resabs;
			*abserr = max(d__1,*abserr);
		}
		/* Computing MAX */
		d__1 = *epsabs, d__2 = *epsrel * fabs(*result);
		if (*abserr <= max(d__1,d__2)) {
			*ier = 0;
		}
		/* ***jump out of do-loop */
		if (*ier == 0) {
			goto L999;
		}
		/* L70: */
	}
L80:
	xerror_("abnormal return from dqng ", &c__26, ier, &c__0, 26);
L999:
	return;
} /* dqng_ */

void dqpsrt_(const long *limit, long *last, long *maxerr, 
				 double *ermax, double *elist, long *iord, long *nrmax)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	long i__, j, k, ido, ibeg, jbnd, isucc, jupbn;
	double errmin, errmax;

	/* ***begin prologue  dqpsrt */
	/* ***refer to  dqage,dqagie,dqagpe,dqawse */
	/* ***routines called  (none) */
	/* ***revision date  810101   (yymmdd) */
	/* ***keywords  sequential sorting */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  this routine maintains the descending ordering in the */
	/*            list of the local error estimated resulting from the */
	/*            interval subdivision process. at each call two error */
	/*            estimates are inserted using the sequential search */
	/*            method, top-down for the largest error estimate and */
	/*            bottom-up for the smallest error estimate. */
	/* ***description */

	/*           ordering routine */
	/*           standard fortran subroutine */
	/*           double precision version */

	/*           parameters (meaning at output) */
	/*              limit  - long */
	/*                       maximum number of error estimates the list */
	/*                       can contain */

	/*              last   - long */
	/*                       number of error estimates currently in the list */

	/*              maxerr - long */
	/*                       maxerr points to the nrmax-th largest error */
	/*                       estimate currently in the list */

	/*              ermax  - double precision */
	/*                       nrmax-th largest error estimate */
	/*                       ermax = elist(maxerr) */

	/*              elist  - double precision */
	/*                       vector of dimension last containing */
	/*                       the error estimates */

	/*              iord   - long */
	/*                       vector of dimension last, the first k elements */
	/*                       of which contain pointers to the error */
	/*                       estimates, such that */
	/*                       elist(iord(1)),...,  elist(iord(k)) */
	/*                       form a decreasing sequence, with */
	/*                       k = last if last.le.(limit/2+2), and */
	/*                       k = limit+1-last otherwise */

	/*              nrmax  - long */
	/*                       maxerr = iord(nrmax) */

	/* ***end prologue  dqpsrt */


	/*           check whether the list contains more than */
	/*           two error estimates. */

	/* ***first executable statement  dqpsrt */
	/* Parameter adjustments */
	--iord;
	--elist;

	/* Function Body */
	if (*last > 2) {
		goto L10;
	}
	iord[1] = 1;
	iord[2] = 2;
	goto L90;

	/*           this part of the routine is only executed if, due to a */
	/*           difficult integrand, subdivision increased the error */
	/*           estimate. in the normal case the insert procedure should */
	/*           start after the nrmax-th largest error estimate. */

L10:
	errmax = elist[*maxerr];
	if (*nrmax == 1) {
		goto L30;
	}
	ido = *nrmax - 1;
	i__1 = ido;
	for (i__ = 1; i__ <= i__1; ++i__) {
		isucc = iord[*nrmax - 1];
		/* ***jump out of do-loop */
		if (errmax <= elist[isucc]) {
			goto L30;
		}
		iord[*nrmax] = isucc;
		--(*nrmax);
		/* L20: */
	}

	/*           compute the number of elements in the list to be maintained */
	/*           in descending order. this number depends on the number of */
	/*           subdivisions still allowed. */

L30:
	jupbn = *last;
	if (*last > *limit / 2 + 2) {
		jupbn = *limit + 3 - *last;
	}
	errmin = elist[*last];

	/*           insert errmax by traversing the list top-down, */
	/*           starting comparison from the element elist(iord(nrmax+1)). */

	jbnd = jupbn - 1;
	ibeg = *nrmax + 1;
	if (ibeg > jbnd) {
		goto L50;
	}
	i__1 = jbnd;
	for (i__ = ibeg; i__ <= i__1; ++i__) {
		isucc = iord[i__];
		/* ***jump out of do-loop */
		if (errmax >= elist[isucc]) {
			goto L60;
		}
		iord[i__ - 1] = isucc;
		/* L40: */
	}
L50:
	iord[jbnd] = *maxerr;
	iord[jupbn] = *last;
	goto L90;

	/*           insert errmin by traversing the list bottom-up. */

L60:
	iord[i__ - 1] = *maxerr;
	k = jbnd;
	i__1 = jbnd;
	for (j = i__; j <= i__1; ++j) {
		isucc = iord[k];
		/* ***jump out of do-loop */
		if (errmin < elist[isucc]) {
			goto L80;
		}
		iord[k + 1] = isucc;
		--k;
		/* L70: */
	}
	iord[i__] = *last;
	goto L90;
L80:
	iord[k + 1] = *last;

	/*           set maxerr and ermax. */

L90:
	*maxerr = iord[*nrmax];
	*ermax = elist[*maxerr];
	return;
} /* dqpsrt_ */

double dqwgtc_(const double *x, const double *c__, const double *, 
					const double *, const double *, const long *)
{
	/* System generated locals */
	double ret_val;

	/* ***begin prologue  dqwgtc */
	/* ***refer to dqk15w */
	/* ***routines called  (none) */
	/* ***revision date  810101   (yymmdd) */
	/* ***keywords  weight function, cauchy principal value */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  this function subprogram is used together with the */
	/*            routine qawc and defines the weight function. */
	/* ***end prologue  dqwgtc */

	/* ***first executable statement  dqwgtc */
	ret_val = 1. / (*x - *c__);
	return ret_val;
} /* dqwgtc_ */

double dqwgtf_(const double *x, const double *omega, const double *,
					const double *, const double *, const long *integr)
{
	/* System generated locals */
	double ret_val;

	/* Local variables */
	double omx;

	/* ***begin prologue  dqwgtf */
	/* ***refer to   dqk15w */
	/* ***routines called  (none) */
	/* ***revision date 810101   (yymmdd) */
	/* ***keywords  cos or sin in weight function */
	/* ***author  piessens,robert, appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. * progr. div. - k.u.leuven */
	/* ***end prologue  dqwgtf */

	/* ***first executable statement  dqwgtf */
	omx = *omega * *x;
	switch (*integr) {
	case 1:  goto L10;
	case 2:  goto L20;
	}
L10:
	ret_val = cos(omx);
	goto L30;
L20:
	ret_val = sin(omx);
L30:
	return ret_val;
} /* dqwgtf_ */

double dqwgts_(const double *x, const double *a, const double *b, 
					const double *alfa, const double *beta, const long *integr)
{
	/* System generated locals */
	double ret_val;

	/* Local variables */
	double xma, bmx;

	/* ***begin prologue  dqwgts */
	/* ***refer to dqk15w */
	/* ***routines called  (none) */
	/* ***revision date  810101   (yymmdd) */
	/* ***keywords  weight function, algebraico-logarithmic */
	/*             end-point singularities */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  this function subprogram is used together with the */
	/*            routine dqaws and defines the weight function. */
	/* ***end prologue  dqwgts */

	/* ***first executable statement  dqwgts */
	xma = *x - *a;
	bmx = *b - *x;
	ret_val = pow(xma, *alfa) * pow(bmx, *beta);
	switch (*integr) {
	case 1:  goto L40;
	case 2:  goto L10;
	case 3:  goto L20;
	case 4:  goto L30;
	}
L10:
	ret_val *= log(xma);
	goto L40;
L20:
	ret_val *= log(bmx);
	goto L40;
L30:
	ret_val = ret_val * log(xma) * log(bmx);
L40:
	return ret_val;
} /* dqwgts_ */

void qage_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  const sys_float *epsabs, const sys_float *epsrel, const long *key, 
			  const long *limit, sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, sys_float *alist__, sys_float *blist, sys_float *rlist,
			  sys_float *elist, long *iord, long *last)
{
	/* System generated locals */
	int i__1;
	sys_float r__1, r__2;

	/* Local variables */
	long k;
	sys_float a1, a2, b1, b2;
	sys_float area;
	long keyf;
	sys_float area1, area2, area12, erro12, defab1, defab2;
	long nrmax;
	sys_float uflow;
	long iroff1, iroff2;
	sys_float error1, error2, defabs, epmach, errbnd, resabs, errmax;
	long maxerr;
	sys_float errsum;

	/* ***begin prologue  qage */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             integrand examinator, globally adaptive, */
	/*             gauss-kronrod */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral   i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-reslt).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            key    - long */
	/*                     key for choice of local integration rule */
	/*                     a gauss-kronrod pair is used with */
	/*                          7 - 15 points if key.lt.2, */
	/*                         10 - 21 points if key = 2, */
	/*                         15 - 31 points if key = 3, */
	/*                         20 - 41 points if key = 4, */
	/*                         25 - 51 points if key = 5, */
	/*                         30 - 61 points if key.gt.5. */

	/*            limit  - long */
	/*                     gives an upperbound on the number of subintervals */
	/*                     in the partition of (a,b), limit.ge.1. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for result and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value */
	/*                             of limit. */
	/*                             however, if this yields no improvement it */
	/*                             is rather advised to analyze the integrand */
	/*                             in order to determine the integration */
	/*                             difficulties. if the position of a local */
	/*                             difficulty can be determined(e.g. */
	/*                             singularity, discontinuity within the */
	/*                             interval) one will probably gain from */
	/*                             splitting up the interval at this polong */
	/*                             and calling the integrator on the */
	/*                             subranges. if possible, an appropriate */
	/*                             special-purpose integrator should be used */
	/*                             which is designed for handling the type of */
	/*                             difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                             result, abserr, neval, last, rlist(1) , */
	/*                             elist(1) and iord(1) are set to zero. */
	/*                             alist(1) and blist(1) are set to a and b */
	/*                             respectively. */

	/*            alist   - sys_float */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the left */
	/*                      end points of the subintervals in the partition */
	/*                      of the given integration range (a,b) */

	/*            blist   - sys_float */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the right */
	/*                      end points of the subintervals in the partition */
	/*                      of the given integration range (a,b) */

	/*            rlist   - sys_float */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the */
	/*                      integral approximations on the subintervals */

	/*            elist   - sys_float */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the moduli of the */
	/*                      absolute error estimates on the subintervals */

	/*            iord    - long */
	/*                      vector of dimension at least limit, the first k */
	/*                      elements of which are pointers to the */
	/*                      error estimates over the subintervals, */
	/*                      such that elist(iord(1)), ..., */
	/*                      elist(iord(k)) form a decreasing sequence, */
	/*                      with k = last if last.le.(limit/2+2), and */
	/*                      k = limit+1-last otherwise */

	/*            last    - long */
	/*                      number of subintervals actually produced in the */
	/*                      subdivision process */

	/* ***references  (none) */
	/* ***routines called  qk15,qk21,qk31,qk41,qk51,qk61,qpsrt,r1mach */
	/* ***end prologue  qage */




	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                      (alist(i),blist(i)) */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest */
	/*                       error estimate */
	/*           errmax    - elist(maxerr) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */


	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach  is the largest relative spacing. */
	/*           uflow  is the smallest positive magnitude. */

	/* ***first executable statement  qage */
	/* Parameter adjustments */
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);

	/*           test on validity of parameters */
	/*           ------------------------------ */

	*ier = 0;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	alist__[1] = *a;
	blist[1] = *b;
	rlist[1] = 0.f;
	elist[1] = 0.f;
	iord[1] = 0;
	/* Computing MAX */
	r__1 = epmach * 50.f;
	if (*epsabs <= 0.f && *epsrel < max(r__1,5e-15f)) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}

	/*           first approximation to the integral */
	/*           ----------------------------------- */

	keyf = *key;
	if (*key <= 0) {
		keyf = 1;
	}
	if (*key >= 7) {
		keyf = 6;
	}
	*neval = 0;
	if (keyf == 1) {
		qk15_(f, a, b, result, abserr, &defabs, &resabs);
	}
	if (keyf == 2) {
		qk21_(f, a, b, result, abserr, &defabs, &resabs);
	}
	if (keyf == 3) {
		qk31_(f, a, b, result, abserr, &defabs, &resabs);
	}
	if (keyf == 4) {
		qk41_(f, a, b, result, abserr, &defabs, &resabs);
	}
	if (keyf == 5) {
		qk51_(f, a, b, result, abserr, &defabs, &resabs);
	}
	if (keyf == 6) {
		qk61_(f, a, b, result, abserr, &defabs, &resabs);
	}
	*last = 1;
	rlist[1] = *result;
	elist[1] = *abserr;
	iord[1] = 1;

	/*           test on accuracy. */

	/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * fabs(*result);
	errbnd = max(r__1,r__2);
	if (*abserr <= epmach * 50.f * defabs && *abserr > errbnd) {
		*ier = 2;
	}
	if (*limit == 1) {
		*ier = 1;
	}
	if (*ier != 0 || ( *abserr <= errbnd && *abserr != resabs ) || *abserr == 0.f)
	{
		goto L60;
	}

	/*           initialization */
	/*           -------------- */


	errmax = *abserr;
	maxerr = 1;
	area = *result;
	errsum = *abserr;
	nrmax = 1;
	iroff1 = 0;
	iroff2 = 0;

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limit;
	for (*last = 2; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with the largest error estimate. */

		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
		a2 = b1;
		b2 = blist[maxerr];
		if (keyf == 1) {
			qk15_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 2) {
			qk21_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 3) {
			qk31_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 4) {
			qk41_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 5) {
			qk51_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 6) {
			qk61_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		}
		if (keyf == 1) {
			qk15_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}
		if (keyf == 2) {
			qk21_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}
		if (keyf == 3) {
			qk31_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}
		if (keyf == 4) {
			qk41_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}
		if (keyf == 5) {
			qk51_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}
		if (keyf == 6) {
			qk61_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);
		}

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		++(*neval);
		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (defab1 == error1 || defab2 == error2) {
			goto L5;
		}
		if ((r__1 = rlist[maxerr] - area12, fabs(r__1)) <= fabs(area12) * 
			 1e-5f && erro12 >= errmax * .99f) {
			++iroff1;
		}
		if (*last > 10 && erro12 > errmax) {
			++iroff2;
		}
	L5:
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(area);
		errbnd = max(r__1,r__2);
		if (errsum <= errbnd) {
			goto L8;
		}

		/*           test for roundoff error and eventually */
		/*           set error flag. */

		if (iroff1 >= 6 || iroff2 >= 20) {
			*ier = 2;
		}

		/*           set error flag in the case that the number of */
		/*           subintervals equals limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at a point of the integration range. */

		/* Computing MAX */
		r__1 = fabs(a1), r__2 = fabs(b2);
		if (max(r__1,r__2) <= (epmach * 100.f + 1.f) * (fabs(a2) + uflow * 
																		 1e3f)) {
			*ier = 3;
		}

		/*           append the newly-created intervals to the list. */

	L8:
		if (error2 > error1) {
			goto L10;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L20;
	L10:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine qpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the */
		/*           subinterval with the largest error estimate (to be */
		/*           bisected next). */

	L20:
		qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (*ier != 0 || errsum <= errbnd) {
			goto L40;
		}
		/* L30: */
	}

	/*           compute final result. */
	/*           --------------------- */

L40:
	*result = 0.f;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L50: */
	}
	*abserr = errsum;
L60:
	if (keyf != 1) {
		*neval = (keyf * 10 + 1) * ((*neval << 1) + 1);
	}
	if (keyf == 1) {
		*neval = *neval * 30 + 15;
	}
L999:
	return;
} /* qage_ */

void qag_(const E_fp& f, const sys_float *a, const sys_float *b, 
			 const sys_float *epsabs, const sys_float *epsrel, const long *key, 
			 sys_float *result, sys_float *abserr, long *neval, 
			 long *ier, long *limit, const long *lenw, long *last, long *iwork, 
			 sys_float *work)
{
	long l1, l2, l3, lvl;

	/* ***begin prologue  qag */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             integrand examinator, globally adaptive, */
	/*             gauss-kronrod */
	/* ***author  piessens,robert,appl. math. & progr. div - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result)le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*            f      - sys_float */
	/*                     function subprogam defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            key    - long */
	/*                     key for choice of local integration rule */
	/*                     a gauss-kronrod pair is used with */
	/*                       7 - 15 points if key.lt.2, */
	/*                      10 - 21 points if key = 2, */
	/*                      15 - 31 points if key = 3, */
	/*                      20 - 41 points if key = 4, */
	/*                      25 - 51 points if key = 5, */
	/*                      30 - 61 points if key.gt.5. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for result and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*                      error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yield no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulaties. */
	/*                             if the position of a local difficulty can */
	/*                             be determined (i.e.singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or limit.lt.1 or lenw.lt.limit*4. */
	/*                             result, abserr, neval, last are set */
	/*                             to zero. */
	/*                             except when lenw is invalid, iwork(1), */
	/*                             work(limit*2+1) and work(limit*3+1) are */
	/*                             set to zero, work(1) is set to a and */
	/*                             work(limit+1) to b. */

	/*         dimensioning parameters */
	/*            limit - long */
	/*                    dimensioning parameter for iwork */
	/*                    limit determines the maximum number of subintervals */
	/*                    in the partition of the given integration interval */
	/*                    (a,b), limit.ge.1. */
	/*                    if limit.lt.1, the routine will end with ier = 6. */

	/*            lenw  - long */
	/*                    dimensioning parameter for work */
	/*                    lenw must be at least limit*4. */
	/*                    if lenw.lt.limit*4, the routine will end with */
	/*                    ier = 6. */

	/*            last  - long */
	/*                    on return, last equals the number of subintervals */
	/*                    produced in the subdivision process, which */
	/*                    determines the number of significant elements */
	/*                    actually in the work arrays. */

	/*         work arrays */
	/*            iwork - long */
	/*                    vector of dimension at least limit, the first k */
	/*                    elements of which contain pointers to the error */
	/*                    estimates over the subintervals, such that */
	/*                    work(limit*3+iwork(1)),... , work(limit*3+iwork(k)) */
	/*                    form a decreasing sequence with k = last if */
	/*                    last.le.(limit/2+2), and k = limit+1-last otherwise */

	/*            work  - sys_float */
	/*                    vector of dimension at least lenw */
	/*                    on return */
	/*                    work(1), ..., work(last) contain the left end */
	/*                    points of the subintervals in the partition of */
	/*                     (a,b), */
	/*                    work(limit+1), ..., work(limit+last) contain the */
	/*                     right end points, */
	/*                    work(limit*2+1), ..., work(limit*2+last) contain */
	/*                     the integral approximations over the subintervals, */
	/*                    work(limit*3+1), ..., work(limit*3+last) contain */
	/*                     the error estimates. */

	/* ***references  (none) */
	/* ***routines called  qage,xerror */
	/* ***end prologue  qag */




	/*         check validity of lenw. */

	/* ***first executable statement  qag */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	if (*limit < 1 || *lenw < *limit << 2) {
		goto L10;
	}

	/*         prepare call for qage. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;

	qage_(f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, 
			ier, &work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from  qag ", &c__26, ier, &lvl, 26);
	}
	return;
} /* qag_ */

void qagie_(const E_fp& f, const sys_float *bound, const long *inf, 
				const sys_float *epsabs, const sys_float *epsrel, 
				const long *limit, 
				sys_float *result, sys_float *abserr, long *neval, 
				long *ier, sys_float *alist__, sys_float *blist, sys_float *rlist,
				sys_float *elist, long *iord, long *last)
{
	/* System generated locals */
	int i__1, i__2;
	sys_float r__1, r__2;

	/* Local variables */
	long k;
	sys_float a1, a2, b1, b2;
	long id;
	sys_float area;
	sys_float dres;
	long ksgn;
	sys_float boun;
	long nres;
	sys_float area1, area2, area12, small=0.f, erro12;
	long ierro;
	sys_float defab1, defab2;
	long ktmin, nrmax;
	sys_float oflow, uflow;
	bool noext;
	long iroff1, iroff2, iroff3;
	sys_float res3la[3], error1, error2, rlist2[52];
	long numrl2;
	sys_float defabs, epmach, erlarg=0.f, abseps, correc=0.f, errbnd, resabs;
	long jupbnd;
	sys_float erlast, errmax;
	long maxerr;
	sys_float reseps;
	bool extrap;
	sys_float ertest=0.f, errsum;

	/* ***begin prologue  qagie */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a3a1,h2a4a1 */
	/* ***keywords  automatic integrator, infinite intervals, */
	/*             general-purpose, transformation, extrapolation, */
	/*             globally adaptive */
	/* ***author  piessens,robert,appl. math & progr. div - k.u.leuven */
	/*           de doncker,elise,appl. math & progr. div - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            integral   i = integral of f over (bound,+infinity) */
	/*                    or i = integral of f over (-infinity,bound) */
	/*                    or i = integral of f over (-infinity,+infinity), */
	/*                    hopefully satisfying following claim for accuracy */
	/*                    fabs(i-result).le.max(epsabs,epsrel*fabs(i)) */
	/* ***description */

	/* integration over infinite intervals */
	/* standard fortran subroutine */

	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            bound  - sys_float */
	/*                     finite bound of integration range */
	/*                     (has no meaning if interval is doubly-infinite) */

	/*            inf    - sys_float */
	/*                     indicating the kind of integration range involved */
	/*                     inf = 1 corresponds to  (bound,+infinity), */
	/*                     inf = -1            to  (-infinity,bound), */
	/*                     inf = 2             to (-infinity,+infinity). */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subintervals */
	/*                     in the partition of (a,b), limit.ge.1 */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                   - ier.gt.0 abnormal termination of the routine. the */
	/*                             estimates for result and error are less */
	/*                             reliable. it is assumed that the requested */
	/*                             accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however,if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. */
	/*                             if the position of a local difficulty can */
	/*                             be determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. */
	/*                             it is assumed that the requested tolerance */
	/*                             cannot be achieved, and that the returned */
	/*                             result is the best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                             result, abserr, neval, last, rlist(1), */
	/*                             elist(1) and iord(1) are set to zero. */
	/*                             alist(1) and blist(1) are set to 0 */
	/*                             and 1 respectively. */

	/*            alist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the left */
	/*                     end points of the subintervals in the partition */
	/*                     of the transformed integration range (0,1). */

	/*            blist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the right */
	/*                     end points of the subintervals in the partition */
	/*                     of the transformed integration range (0,1). */

	/*            rlist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the integral */
	/*                     approximations on the subintervals */

	/*            elist  - sys_float */
	/*                     vector of dimension at least limit,  the first */
	/*                     last elements of which are the moduli of the */
	/*                     absolute error estimates on the subintervals */

	/*            iord   - long */
	/*                     vector of dimension limit, the first k */
	/*                     elements of which are pointers to the */
	/*                     error estimates over the subintervals, */
	/*                     such that elist(iord(1)), ..., elist(iord(k)) */
	/*                     form a decreasing sequence, with k = last */
	/*                     if last.le.(limit/2+2), and k = limit+1-last */
	/*                     otherwise */

	/*            last   - long */
	/*                     number of subintervals actually produced */
	/*                     in the subdivision process */

	/* ***references  (none) */
	/* ***routines called  qelg,qk15i,qpsrt,r1mach */
	/* ***end prologue  qagie */




	/*            the dimension of rlist2 is determined by the value of */
	/*            limexp in subroutine qelg. */


	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           rlist2    - array of dimension at least (limexp+2), */
	/*                       containing the part of the epsilon table */
	/*                       wich is still needed for further computations */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest error */
	/*                       estimate */
	/*           errmax    - elist(maxerr) */
	/*           erlast    - error on the interval currently subdivided */
	/*                       (before that subdivision has taken place) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */
	/*           nres      - number of calls to the extrapolation routine */
	/*           numrl2    - number of elements currently in rlist2. if an */
	/*                       appropriate approximation to the compounded */
	/*                       integral has been obtained, it is put in */
	/*                       rlist2(numrl2) after numrl2 has been increased */
	/*                       by one. */
	/*           small     - length of the smallest interval considered up */
	/*                       to now, multiplied by 1.5 */
	/*           erlarg    - sum of the errors over the intervals larger */
	/*                       than the smallest interval considered up to now */
	/*           extrap    - bool variable denoting that the routine */
	/*                       is attempting to perform extrapolation. i.e. */
	/*                       before subdividing the smallest interval we */
	/*                       try to decrease the value of erlarg. */
	/*           noext     - bool variable denoting that extrapolation */
	/*                       is no longer allowed (true-value) */

	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */
	/*           oflow is the largest positive magnitude. */

	/* Parameter adjustments */
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = r1mach(c__4);

	/*           test on validity of parameters */
	/*           ----------------------------- */

	/* ***first executable statement  qagie */
	*ier = 0;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	alist__[1] = 0.f;
	blist[1] = 1.f;
	rlist[1] = 0.f;
	elist[1] = 0.f;
	iord[1] = 0;
	/* Computing MAX */
	r__1 = epmach * 50.f;
	if (*epsabs <= 0.f && *epsrel < max(r__1,5e-15f)) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}


	/*           first approximation to the integral */
	/*           ----------------------------------- */

	/*           determine the interval to be mapped onto (0,1). */
	/*           if inf = 2 the integral is computed as i = i1+i2, where */
	/*           i1 = integral of f over (-infinity,0), */
	/*           i2 = integral of f over (0,+infinity). */

	boun = *bound;
	if (*inf == 2) {
		boun = 0.f;
	}
	qk15i_(f, &boun, inf, &c_b390, &c_b391, result, abserr, &defabs, &
			 resabs);

	/*           test on accuracy */

	*last = 1;
	rlist[1] = *result;
	elist[1] = *abserr;
	iord[1] = 1;
	dres = fabs(*result);
	/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dres;
	errbnd = max(r__1,r__2);
	if (*abserr <= epmach * 100.f * defabs && *abserr > errbnd) {
		*ier = 2;
	}
	if (*limit == 1) {
		*ier = 1;
	}
	if (*ier != 0 || ( *abserr <= errbnd && *abserr != resabs ) || *abserr == 0.f)
	{
		goto L130;
	}

	/*           initialization */
	/*           -------------- */

	uflow = r1mach(c__1);
	oflow = r1mach(c__2);
	rlist2[0] = *result;
	errmax = *abserr;
	maxerr = 1;
	area = *result;
	errsum = *abserr;
	*abserr = oflow;
	nrmax = 1;
	nres = 0;
	ktmin = 0;
	numrl2 = 2;
	extrap = false;
	noext = false;
	ierro = 0;
	iroff1 = 0;
	iroff2 = 0;
	iroff3 = 0;
	ksgn = -1;
	if (dres >= (1.f - epmach * 50.f) * defabs) {
		ksgn = 1;
	}

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limit;
	for (*last = 2; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with nrmax-th largest */
		/*           error estimate. */

		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
		a2 = b1;
		b2 = blist[maxerr];
		erlast = errmax;
		qk15i_(f, &boun, inf, &a1, &b1, &area1, &error1, &resabs, &
				 defab1);
		qk15i_(f, &boun, inf, &a2, &b2, &area2, &error2, &resabs, &
				 defab2);

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (defab1 == error1 || defab2 == error2) {
			goto L15;
		}
		if ((r__1 = rlist[maxerr] - area12, fabs(r__1)) > fabs(area12) * 
			 1e-5f || erro12 < errmax * .99f) {
			goto L10;
		}
		if (extrap) {
			++iroff2;
		}
		if (! extrap) {
			++iroff1;
		}
	L10:
		if (*last > 10 && erro12 > errmax) {
			++iroff3;
		}
	L15:
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(area);
		errbnd = max(r__1,r__2);

		/*           test for roundoff error and eventually */
		/*           set error flag. */

		if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
			*ier = 2;
		}
		if (iroff2 >= 5) {
			ierro = 3;
		}

		/*           set error flag in the case that the number of */
		/*           subintervals equals limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at some points of the integration range. */

		/* Computing MAX */
		r__1 = fabs(a1), r__2 = fabs(b2);
		if (max(r__1,r__2) <= (epmach * 100.f + 1.f) * (fabs(a2) + uflow * 
																		 1e3f)) {
			*ier = 4;
		}

		/*           append the newly-created intervals to the list. */

		if (error2 > error1) {
			goto L20;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L30;
	L20:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine qpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the */
		/*           subinterval with nrmax-th largest error estimate (to be */
		/*           bisected next). */

	L30:
		qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		if (errsum <= errbnd) {
			goto L115;
		}
		if (*ier != 0) {
			goto L100;
		}
		if (*last == 2) {
			goto L80;
		}
		if (noext) {
			goto L90;
		}
		erlarg -= erlast;
		if ((r__1 = b1 - a1, fabs(r__1)) > small) {
			erlarg += erro12;
		}
		if (extrap) {
			goto L40;
		}

		/*           test whether the interval to be bisected next is the */
		/*           smallest interval. */

		if ((r__1 = blist[maxerr] - alist__[maxerr], fabs(r__1)) > small) {
			goto L90;
		}
		extrap = true;
		nrmax = 2;
	L40:
		if (ierro == 3 || erlarg <= ertest) {
			goto L60;
		}

		/*           the smallest interval has the largest error. */
		/*           before bisecting decrease the sum of the errors */
		/*           over the larger intervals (erlarg) and perform */
		/*           extrapolation. */

		id = nrmax;
		jupbnd = *last;
		if (*last > *limit / 2 + 2) {
			jupbnd = *limit + 3 - *last;
		}
		i__2 = jupbnd;
		for (k = id; k <= i__2; ++k) {
			maxerr = iord[nrmax];
			errmax = elist[maxerr];
			if ((r__1 = blist[maxerr] - alist__[maxerr], fabs(r__1)) > small) 
			{
				goto L90;
			}
			++nrmax;
			/* L50: */
		}

		/*           perform extrapolation. */

	L60:
		++numrl2;
		rlist2[numrl2 - 1] = area;
		qelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
		++ktmin;
		if (ktmin > 5 && *abserr < errsum * .001f) {
			*ier = 5;
		}
		if (abseps >= *abserr) {
			goto L70;
		}
		ktmin = 0;
		*abserr = abseps;
		*result = reseps;
		correc = erlarg;
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(reseps);
		ertest = max(r__1,r__2);
		if (*abserr <= ertest) {
			goto L100;
		}

		/*            prepare bisection of the smallest interval. */

	L70:
		if (numrl2 == 1) {
			noext = true;
		}
		if (*ier == 5) {
			goto L100;
		}
		maxerr = iord[1];
		errmax = elist[maxerr];
		nrmax = 1;
		extrap = false;
		small *= .5f;
		erlarg = errsum;
		goto L90;
	L80:
		small = .375f;
		erlarg = errsum;
		ertest = errbnd;
		rlist2[1] = area;
	L90:
		;
	}

	/*           set final result and error estimate. */
	/*           ------------------------------------ */

L100:
	if (*abserr == oflow) {
		goto L115;
	}
	if (*ier + ierro == 0) {
		goto L110;
	}
	if (ierro == 3) {
		*abserr += correc;
	}
	if (*ier == 0) {
		*ier = 3;
	}
	if (*result != 0.f && area != 0.f) {
		goto L105;
	}
	if (*abserr > errsum) {
		goto L115;
	}
	if (area == 0.f) {
		goto L130;
	}
	goto L110;
L105:
	if (*abserr / fabs(*result) > errsum / fabs(area)) {
		goto L115;
	}

	/*           test on divergence */

L110:
	/* Computing MAX */
	r__1 = fabs(*result), r__2 = fabs(area);
	if (ksgn == -1 && max(r__1,r__2) <= defabs * .01f) {
		goto L130;
	}
	if (.01f > *result / area || *result / area > 100.f || errsum > fabs(area)
		) {
		*ier = 6;
	}
	goto L130;

	/*           compute global integral sum. */

L115:
	*result = 0.f;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L120: */
	}
	*abserr = errsum;
L130:
	*neval = *last * 30 - 15;
	if (*inf == 2) {
		*neval <<= 1;
	}
	if (*ier > 2) {
		--(*ier);
	}
L999:
	return;
} /* qagie_ */

void qagi_(const E_fp& f, const sys_float *bound, const long *inf, 
			  const sys_float *epsabs, const sys_float *epsrel, 
			  sys_float *result, sys_float *abserr, long *neval, 
			  long *ier, const long *limit, const long *lenw, long *last, 
			  long *iwork, sys_float *work)
{
	long l1, l2, l3, lvl;

	/* ***begin prologue  qagi */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a3a1,h2a4a1 */
	/* ***keywords  automatic integrator, infinite intervals, */
	/*             general-purpose, transformation, extrapolation, */
	/*             globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. -k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            integral   i = integral of f over (bound,+infinity) */
	/*                    or i = integral of f over (-infinity,bound) */
	/*                    or i = integral of f over (-infinity,+infinity) */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        integration over infinite intervals */
	/*        standard fortran subroutine */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            bound  - sys_float */
	/*                     finite bound of integration range */
	/*                     (has no meaning if interval is doubly-infinite) */

	/*            inf    - long */
	/*                     indicating the kind of integration range involved */
	/*                     inf = 1 corresponds to  (bound,+infinity), */
	/*                     inf = -1            to  (-infinity,bound), */
	/*                     inf = 2             to (-infinity,+infinity). */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */


	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                   - ier.gt.0 abnormal termination of the routine. the */
	/*                             estimates for result and error are less */
	/*                             reliable. it is assumed that the requested */
	/*                             accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. if */
	/*                             the position of a local difficulty can be */
	/*                             determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. */
	/*                             it is assumed that the requested tolerance */
	/*                             cannot be achieved, and that the returned */
	/*                             result is the best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                              or limit.lt.1 or leniw.lt.limit*4. */
	/*                             result, abserr, neval, last are set to */
	/*                             zero. exept when limit or leniw is */
	/*                             invalid, iwork(1), work(limit*2+1) and */
	/*                             work(limit*3+1) are set to zero, work(1) */
	/*                             is set to a and work(limit+1) to b. */

	/*         dimensioning parameters */
	/*            limit - long */
	/*                    dimensioning parameter for iwork */
	/*                    limit determines the maximum number of subintervals */
	/*                    in the partition of the given integration interval */
	/*                    (a,b), limit.ge.1. */
	/*                    if limit.lt.1, the routine will end with ier = 6. */

	/*            lenw  - long */
	/*                    dimensioning parameter for work */
	/*                    lenw must be at least limit*4. */
	/*                    if lenw.lt.limit*4, the routine will end */
	/*                    with ier = 6. */

	/*            last  - long */
	/*                    on return, last equals the number of subintervals */
	/*                    produced in the subdivision process, which */
	/*                    determines the number of significant elements */
	/*                    actually in the work arrays. */

	/*         work arrays */
	/*            iwork - long */
	/*                    vector of dimension at least limit, the first */
	/*                    k elements of which contain pointers */
	/*                    to the error estimates over the subintervals, */
	/*                    such that work(limit*3+iwork(1)),... , */
	/*                    work(limit*3+iwork(k)) form a decreasing */
	/*                    sequence, with k = last if last.le.(limit/2+2), and */
	/*                    k = limit+1-last otherwise */

	/*            work  - sys_float */
	/*                    vector of dimension at least lenw */
	/*                    on return */
	/*                    work(1), ..., work(last) contain the left */
	/*                     end points of the subintervals in the */
	/*                     partition of (a,b), */
	/*                    work(limit+1), ..., work(limit+last) contain */
	/*                     the right end points, */
	/*                    work(limit*2+1), ...,work(limit*2+last) contain the */
	/*                     integral approximations over the subintervals, */
	/*                    work(limit*3+1), ..., work(limit*3) */
	/*                     contain the error estimates. */
	/* ***references  (none) */
	/* ***routines called  qagie,xerror */
	/* ***end prologue  qagi */




	/*         check validity of limit and lenw. */

	/* ***first executable statement  qagi */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	if (*limit < 1 || *lenw < *limit << 2) {
		goto L10;
	}

	/*         prepare call for qagie. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;

	qagie_(f, bound, inf, epsabs, epsrel, limit, result, abserr, neval, 
			 ier, &work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from  qagi", &c__26, ier, &lvl, 26);
	}
	return;
} /* qagi_ */

void qagpe_(const E_fp& f, const sys_float *a, const sys_float *b, const long *npts2, 
				const sys_float *points, const sys_float *epsabs, 
				const sys_float *epsrel, const long *limit, sys_float *result, 
				sys_float *abserr, long *neval, long *ier, sys_float *alist__, 
				sys_float *blist, sys_float *rlist, sys_float *elist, 
				sys_float *pts, long *iord, long *level, long *ndin, long *last)
{
	/* System generated locals */
	int i__1, i__2;
	sys_float r__1, r__2;

	/* Local variables */
	long i__, j, k=0;
	sys_float a1, a2, b1, b2;
	long id, ip1;
	long ind1, ind2;
	sys_float area;
	sys_float resa, dres, sign=0.f;
	long ksgn;
	sys_float temp;
	long nres, nint, jlow, npts;
	sys_float area1, area2, area12, erro12;
	long ierro;
	sys_float defab1, defab2;
	long ktmin, nrmax;
	sys_float oflow, uflow;
	bool noext;
	long iroff1, iroff2, iroff3;
	sys_float res3la[3];
	long nintp1;
	sys_float error1, error2, rlist2[52];
	long numrl2;
	sys_float defabs, epmach, erlarg, abseps, correc=0.f, errbnd, resabs;
	long jupbnd;
	sys_float erlast;
	long levmax;
	sys_float errmax;
	long maxerr, levcur;
	sys_float reseps;
	bool extrap;
	sys_float ertest, errsum;

	/* ***begin prologue  qagpe */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             singularities at user specified points, */
	/*             extrapolation, globally adaptive. */
	/* ***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f over (a,b),hopefully */
	/*            satisfying following claim for accuracy fabs(i-result).le. */
	/*            max(epsabs,epsrel*fabs(i)). break points of the integration */
	/*            interval, where local difficulties of the integrand may */
	/*            occur(e.g. singularities,discontinuities),provided by user. */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration */

	/*            npts2  - long */
	/*                     number equal to two more than the number of */
	/*                     user-supplied break points within the integration */
	/*                     range, npts2.ge.2. */
	/*                     if npts2.lt.2, the routine will end with ier = 6. */

	/*            points - sys_float */
	/*                     vector of dimension npts2, the first (npts2-2) */
	/*                     elements of which are the user provided break */
	/*                     points. if these points do not constitute an */
	/*                     ascending sequence there will be an automatic */
	/*                     sorting. */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subintervals */
	/*                     in the partition of (a,b), limit.ge.npts2 */
	/*                     if limit.lt.npts2, the routine will end with */
	/*                     ier = 6. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine. */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. if */
	/*                             the position of a local difficulty can be */
	/*                             determined (i.e. singularity, */
	/*                             discontinuity within the interval), it */
	/*                             should be supplied to the routine as an */
	/*                             element of the vector points. if necessary */
	/*                             an appropriate special-purpose integrator */
	/*                             must be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. it is presumed that */
	/*                             the requested tolerance cannot be */
	/*                             achieved, and that the returned result is */
	/*                             the best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier.gt.0. */
	/*                         = 6 the input is invalid because */
	/*                             npts2.lt.2 or */
	/*                             break points are specified outside */
	/*                             the integration range or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or limit.lt.npts2. */
	/*                             result, abserr, neval, last, rlist(1), */
	/*                             and elist(1) are set to zero. alist(1) and */
	/*                             blist(1) are set to a and b respectively. */

	/*            alist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the left end points */
	/*                     of the subintervals in the partition of the given */
	/*                     integration range (a,b) */

	/*            blist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the right end points */
	/*                     of the subintervals in the partition of the given */
	/*                     integration range (a,b) */

	/*            rlist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the integral */
	/*                     approximations on the subintervals */

	/*            elist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the moduli of the */
	/*                     absolute error estimates on the subintervals */

	/*            pts    - sys_float */
	/*                     vector of dimension at least npts2, containing the */
	/*                     integration limits and the break points of the */
	/*                     interval in ascending sequence. */

	/*            level  - long */
	/*                     vector of dimension at least limit, containing the */
	/*                     subdivision levels of the subinterval, i.e. if */
	/*                     (aa,bb) is a subinterval of (p1,p2) where p1 as */
	/*                     well as p2 is a user-provided break point or */
	/*                     integration limit, then (aa,bb) has level l if */
	/*                     fabs(bb-aa) = fabs(p2-p1)*2**(-l). */

	/*            ndin   - long */
	/*                     vector of dimension at least npts2, after first */
	/*                     integration over the intervals (pts(i)),pts(i+1), */
	/*                     i = 0,1, ..., npts2-2, the error estimates over */
	/*                     some of the intervals may have been increased */
	/*                     artificially, in order to put their subdivision */
	/*                     forward. if this happens for the subinterval */
	/*                     numbered k, ndin(k) is put to 1, otherwise */
	/*                     ndin(k) = 0. */

	/*            iord   - long */
	/*                     vector of dimension at least limit, the first k */
	/*                     elements of which are pointers to the */
	/*                     error estimates over the subintervals, */
	/*                     such that elist(iord(1)), ..., elist(iord(k)) */
	/*                     form a decreasing sequence, with k = last */
	/*                     if last.le.(limit/2+2), and k = limit+1-last */
	/*                     otherwise */

	/*            last   - long */
	/*                     number of subintervals actually produced in the */
	/*                     subdivisions process */

	/* ***references  (none) */
	/* ***routines called  qelg,qk21,qpsrt,r1mach */
	/* ***end prologue  qagpe */




	/*            the dimension of rlist2 is determined by the value of */
	/*            limexp in subroutine epsalg (rlist2 should be of dimension */
	/*            (limexp+2) at least). */


	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           rlist2    - array of dimension at least limexp+2 */
	/*                       containing the part of the epsilon table which */
	/*                       is still needed for further computations */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest error */
	/*                       estimate */
	/*           errmax    - elist(maxerr) */
	/*           erlast    - error on the interval currently subdivided */
	/*                       (before that subdivision has taken place) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */
	/*           nres      - number of calls to the extrapolation routine */
	/*           numrl2    - number of elements in rlist2. if an */
	/*                       appropriate approximation to the compounded */
	/*                       integral has been obtained, it is put in */
	/*                       rlist2(numrl2) after numrl2 has been increased */
	/*                       by one. */
	/*           erlarg    - sum of the errors over the intervals larger */
	/*                       than the smallest interval considered up to now */
	/*           extrap    - bool variable denoting that the routine */
	/*                       is attempting to perform extrapolation. i.e. */
	/*                       before subdividing the smallest interval we */
	/*                       try to decrease the value of erlarg. */
	/*           noext     - bool variable denoting that extrapolation is */
	/*                       no longer allowed (true-value) */

	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */
	/*           oflow is the largest positive magnitude. */

	/* ***first executable statement  qagpe */
	/* Parameter adjustments */
	--ndin;
	--pts;
	--points;
	--level;
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = r1mach(c__4);

	/*            test on validity of parameters */
	/*            ----------------------------- */

	*ier = 0;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	alist__[1] = *a;
	blist[1] = *b;
	rlist[1] = 0.f;
	elist[1] = 0.f;
	iord[1] = 0;
	level[1] = 0;
	npts = *npts2 - 2;
	/* Computing MAX */
	r__1 = epmach * 50.f;
	if (*npts2 < 2 || *limit <= npts || 
		 ( *epsabs <= 0.f && *epsrel < max(r__1, 5e-15f)) ) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L210;
	}

	/*            if any break points are provided, sort them into an */
	/*            ascending sequence. */

	sign = 1.f;
	if (*a > *b) {
		sign = -1.f;
	}
	pts[1] = min(*a,*b);
	if (npts == 0) {
		goto L15;
	}
	i__1 = npts;
	for (i__ = 1; i__ <= i__1; ++i__) {
		pts[i__ + 1] = points[i__];
		/* L10: */
	}
L15:
	pts[npts + 2] = max(*a,*b);
	nint = npts + 1;
	a1 = pts[1];
	if (npts == 0) {
		goto L40;
	}
	nintp1 = nint + 1;
	i__1 = nint;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ip1 = i__ + 1;
		i__2 = nintp1;
		for (j = ip1; j <= i__2; ++j) {
			if (pts[i__] <= pts[j]) {
				goto L20;
			}
			temp = pts[i__];
			pts[i__] = pts[j];
			pts[j] = temp;
		L20:
			;
		}
	}
	if (pts[1] != min(*a,*b) || pts[nintp1] != max(*a,*b)) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}

	/*            compute first integral and error approximations. */
	/*            ------------------------------------------------ */

L40:
	resabs = 0.f;
	i__2 = nint;
	for (i__ = 1; i__ <= i__2; ++i__) {
		b1 = pts[i__ + 1];
		qk21_(f, &a1, &b1, &area1, &error1, &defabs, &resa);
		*abserr += error1;
		*result += area1;
		ndin[i__] = 0;
		if (error1 == resa && error1 != 0.f) {
			ndin[i__] = 1;
		}
		resabs += defabs;
		level[i__] = 0;
		elist[i__] = error1;
		alist__[i__] = a1;
		blist[i__] = b1;
		rlist[i__] = area1;
		iord[i__] = i__;
		a1 = b1;
		/* L50: */
	}
	errsum = 0.f;
	i__2 = nint;
	for (i__ = 1; i__ <= i__2; ++i__) {
		if (ndin[i__] == 1) {
			elist[i__] = *abserr;
		}
		errsum += elist[i__];
		/* L55: */
	}

	/*           test on accuracy. */

	*last = nint;
	*neval = nint * 21;
	dres = fabs(*result);
	/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dres;
	errbnd = max(r__1,r__2);
	if (*abserr <= epmach * 100.f * resabs && *abserr > errbnd) {
		*ier = 2;
	}
	if (nint == 1) {
		goto L80;
	}
	i__2 = npts;
	for (i__ = 1; i__ <= i__2; ++i__) {
		jlow = i__ + 1;
		ind1 = iord[i__];
		i__1 = nint;
		for (j = jlow; j <= i__1; ++j) {
			ind2 = iord[j];
			if (elist[ind1] > elist[ind2]) {
				goto L60;
			}
			ind1 = ind2;
			k = j;
		L60:
			;
		}
		if (ind1 == iord[i__]) {
			goto L70;
		}
		iord[k] = iord[i__];
		iord[i__] = ind1;
	L70:
		;
	}
	if (*limit < *npts2) {
		*ier = 1;
	}
L80:
	if (*ier != 0 || *abserr <= errbnd) {
		goto L999;
	}

	/*           initialization */
	/*           -------------- */

	rlist2[0] = *result;
	maxerr = iord[1];
	errmax = elist[maxerr];
	area = *result;
	nrmax = 1;
	nres = 0;
	numrl2 = 1;
	ktmin = 0;
	extrap = false;
	noext = false;
	erlarg = errsum;
	ertest = errbnd;
	levmax = 1;
	iroff1 = 0;
	iroff2 = 0;
	iroff3 = 0;
	ierro = 0;
	uflow = r1mach(c__1);
	oflow = r1mach(c__2);
	*abserr = oflow;
	ksgn = -1;
	if (dres >= (1.f - epmach * 50.f) * resabs) {
		ksgn = 1;
	}

	/*           main do-loop */
	/*           ------------ */

	i__2 = *limit;
	for (*last = *npts2; *last <= i__2; ++(*last)) {

		/*           bisect the subinterval with the nrmax-th largest */
		/*           error estimate. */

		levcur = level[maxerr] + 1;
		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
		a2 = b1;
		b2 = blist[maxerr];
		erlast = errmax;
		qk21_(f, &a1, &b1, &area1, &error1, &resa, &defab1);
		qk21_(f, &a2, &b2, &area2, &error2, &resa, &defab2);

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		*neval += 42;
		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (defab1 == error1 || defab2 == error2) {
			goto L95;
		}
		if ((r__1 = rlist[maxerr] - area12, fabs(r__1)) > fabs(area12) * 
			 1e-5f || erro12 < errmax * .99f) {
			goto L90;
		}
		if (extrap) {
			++iroff2;
		}
		if (! extrap) {
			++iroff1;
		}
	L90:
		if (*last > 10 && erro12 > errmax) {
			++iroff3;
		}
	L95:
		level[maxerr] = levcur;
		level[*last] = levcur;
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(area);
		errbnd = max(r__1,r__2);

		/*           test for roundoff error and eventually */
		/*           set error flag. */

		if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
			*ier = 2;
		}
		if (iroff2 >= 5) {
			ierro = 3;
		}

		/*           set error flag in the case that the number of */
		/*           subintervals equals limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at a point of the integration range */

		/* Computing MAX */
		r__1 = fabs(a1), r__2 = fabs(b2);
		if (max(r__1,r__2) <= (epmach * 100.f + 1.f) * (fabs(a2) + uflow * 
																		 1e3f)) {
			*ier = 4;
		}

		/*           append the newly-created intervals to the list. */

		if (error2 > error1) {
			goto L100;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L110;
	L100:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine qpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the */
		/*           subinterval with nrmax-th largest error estimate (to be */
		/*           bisected next). */

	L110:
		qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (errsum <= errbnd) {
			goto L190;
		}
		/* ***jump out of do-loop */
		if (*ier != 0) {
			goto L170;
		}
		if (noext) {
			goto L160;
		}
		erlarg -= erlast;
		if (levcur + 1 <= levmax) {
			erlarg += erro12;
		}
		if (extrap) {
			goto L120;
		}

		/*           test whether the interval to be bisected next is the */
		/*           smallest interval. */

		if (level[maxerr] + 1 <= levmax) {
			goto L160;
		}
		extrap = true;
		nrmax = 2;
	L120:
		if (ierro == 3 || erlarg <= ertest) {
			goto L140;
		}

		/*           the smallest interval has the largest error. */
		/*           before bisecting decrease the sum of the errors */
		/*           over the larger intervals (erlarg) and perform */
		/*           extrapolation. */

		id = nrmax;
		jupbnd = *last;
		if (*last > *limit / 2 + 2) {
			jupbnd = *limit + 3 - *last;
		}
		i__1 = jupbnd;
		for (k = id; k <= i__1; ++k) {
			maxerr = iord[nrmax];
			errmax = elist[maxerr];
			/* ***jump out of do-loop */
			if (level[maxerr] + 1 <= levmax) {
				goto L160;
			}
			++nrmax;
			/* L130: */
		}

		/*           perform extrapolation. */

	L140:
		++numrl2;
		rlist2[numrl2 - 1] = area;
		if (numrl2 <= 2) {
			goto L155;
		}
		qelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
		++ktmin;
		if (ktmin > 5 && *abserr < errsum * .001f) {
			*ier = 5;
		}
		if (abseps >= *abserr) {
			goto L150;
		}
		ktmin = 0;
		*abserr = abseps;
		*result = reseps;
		correc = erlarg;
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(reseps);
		ertest = max(r__1,r__2);
		/* ***jump out of do-loop */
		if (*abserr < ertest) {
			goto L170;
		}

		/*           prepare bisection of the smallest interval. */

	L150:
		if (numrl2 == 1) {
			noext = true;
		}
		if (*ier >= 5) {
			goto L170;
		}
	L155:
		maxerr = iord[1];
		errmax = elist[maxerr];
		nrmax = 1;
		extrap = false;
		++levmax;
		erlarg = errsum;
	L160:
		;
	}

	/*           set the final result. */
	/*           --------------------- */


L170:
	if (*abserr == oflow) {
		goto L190;
	}
	if (*ier + ierro == 0) {
		goto L180;
	}
	if (ierro == 3) {
		*abserr += correc;
	}
	if (*ier == 0) {
		*ier = 3;
	}
	if (*result != 0.f && area != 0.f) {
		goto L175;
	}
	if (*abserr > errsum) {
		goto L190;
	}
	if (area == 0.f) {
		goto L210;
	}
	goto L180;
L175:
	if (*abserr / fabs(*result) > errsum / fabs(area)) {
		goto L190;
	}

	/*           test on divergence. */

L180:
	/* Computing MAX */
	r__1 = fabs(*result), r__2 = fabs(area);
	if (ksgn == -1 && max(r__1,r__2) <= resabs * .01f) {
		goto L210;
	}
	if (.01f > *result / area || *result / area > 100.f || errsum > fabs(area)
		) {
		*ier = 6;
	}
	goto L210;

	/*           compute global integral sum. */

L190:
	*result = 0.f;
	i__2 = *last;
	for (k = 1; k <= i__2; ++k) {
		*result += rlist[k];
		/* L200: */
	}
	*abserr = errsum;
L210:
	if (*ier > 2) {
		--(*ier);
	}
	*result *= sign;
L999:
	return;
} /* qagpe_ */

void qagp_(const E_fp& f, const sys_float *a, const sys_float *b, const long *npts2, 
			  const sys_float *points, const sys_float *epsabs, 
			  const sys_float *epsrel, sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, const long *leniw, const long *lenw, 
			  long *last, long *iwork, sys_float *work)
{
	long l1, l2, l3, l4, lvl;
	long limit;

	/* ***begin prologue  qagp */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             singularities at user specified points, */
	/*             extrapolation, globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            break points of the integration interval, where local */
	/*            difficulties of the integrand may occur(e.g. singularities, */
	/*            discontinuities), are provided by the user. */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration */

	/*            npts2  - long */
	/*                     number equal to two more than the number of */
	/*                     user-supplied break points within the integration */
	/*                     range, npts.ge.2. */
	/*                     if npts2.lt.2, the routine will end with ier = 6. */

	/*            points - sys_float */
	/*                     vector of dimension npts2, the first (npts2-2) */
	/*                     elements of which are the user provided break */
	/*                     points. if these points do not constitute an */
	/*                     ascending sequence there will be an automatic */
	/*                     sorting. */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine. */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. if */
	/*                             the position of a local difficulty can be */
	/*                             determined (i.e. singularity, */
	/*                             discontinuity within the interval), it */
	/*                             should be supplied to the routine as an */
	/*                             element of the vector points. if necessary */
	/*                             an appropriate special-purpose integrator */
	/*                             must be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. */
	/*                             it is presumed that the requested */
	/*                             tolerance cannot be achieved, and that */
	/*                             the returned result is the best which */
	/*                             can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier.gt.0. */
	/*                         = 6 the input is invalid because */
	/*                             npts2.lt.2 or */
	/*                             break points are specified outside */
	/*                             the integration range or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             result, abserr, neval, last are set to */
	/*                             zero. exept when leniw or lenw or npts2 is */
	/*                             invalid, iwork(1), iwork(limit+1), */
	/*                             work(limit*2+1) and work(limit*3+1) */
	/*                             are set to zero. */
	/*                             work(1) is set to a and work(limit+1) */
	/*                             to b (where limit = (leniw-npts2)/2). */

	/*         dimensioning parameters */
	/*            leniw - long */
	/*                    dimensioning parameter for iwork */
	/*                    leniw determines limit = (leniw-npts2)/2, */
	/*                    which is the maximum number of subintervals in the */
	/*                    partition of the given integration interval (a,b), */
	/*                    leniw.ge.(3*npts2-2). */
	/*                    if leniw.lt.(3*npts2-2), the routine will end with */
	/*                    ier = 6. */

	/*            lenw  - long */
	/*                    dimensioning parameter for work */
	/*                    lenw must be at least leniw*2-npts2. */
	/*                    if lenw.lt.leniw*2-npts2, the routine will end */
	/*                    with ier = 6. */

	/*            last  - long */
	/*                    on return, last equals the number of subintervals */
	/*                    produced in the subdivision process, which */
	/*                    determines the number of significant elements */
	/*                    actually in the work arrays. */

	/*         work arrays */
	/*            iwork - long */
	/*                    vector of dimension at least leniw. on return, */
	/*                    the first k elements of which contain */
	/*                    pointers to the error estimates over the */
	/*                    subintervals, such that work(limit*3+iwork(1)),..., */
	/*                    work(limit*3+iwork(k)) form a decreasing */
	/*                    sequence, with k = last if last.le.(limit/2+2), and */
	/*                    k = limit+1-last otherwise */
	/*                    iwork(limit+1), ...,iwork(limit+last) contain the */
	/*                     subdivision levels of the subintervals, i.e. */
	/*                     if (aa,bb) is a subinterval of (p1,p2) */
	/*                     where p1 as well as p2 is a user-provided */
	/*                     break point or integration limit, then (aa,bb) has */
	/*                     level l if fabs(bb-aa) = fabs(p2-p1)*2**(-l), */
	/*                    iwork(limit*2+1), ..., iwork(limit*2+npts2) have */
	/*                     no significance for the user, */
	/*                    note that limit = (leniw-npts2)/2. */

	/*            work  - sys_float */
	/*                    vector of dimension at least lenw */
	/*                    on return */
	/*                    work(1), ..., work(last) contain the left */
	/*                     end points of the subintervals in the */
	/*                     partition of (a,b), */
	/*                    work(limit+1), ..., work(limit+last) contain */
	/*                     the right end points, */
	/*                    work(limit*2+1), ..., work(limit*2+last) contain */
	/*                     the integral approximations over the subintervals, */
	/*                    work(limit*3+1), ..., work(limit*3+last) */
	/*                     contain the corresponding error estimates, */
	/*                    work(limit*4+1), ..., work(limit*4+npts2) */
	/*                     contain the integration limits and the */
	/*                     break points sorted in an ascending sequence. */
	/*                    note that limit = (leniw-npts2)/2. */

	/* ***references  (none) */
	/* ***routines called  qagpe,xerror */
	/* ***end prologue  qagp */




	/*         check validity of limit and lenw. */

	/* ***first executable statement  qagp */
	/* Parameter adjustments */
	--points;
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	if (*leniw < *npts2 * 3 - 2 || *lenw < (*leniw << 1) - *npts2 || *npts2 < 
	    2) {
		goto L10;
	}

	/*         prepare call for qagpe. */

	limit = (*leniw - *npts2) / 2;
	l1 = limit + 1;
	l2 = limit + l1;
	l3 = limit + l2;
	l4 = limit + l3;

	qagpe_(f, a, b, npts2, &points[1], epsabs, epsrel, &limit, result, 
			 abserr, neval, ier, &work[1], &work[l1], &work[l2], &work[l3], &
			 work[l4], &iwork[1], &iwork[l1], &iwork[l2], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from  qagp", &c__26, ier, &lvl, 26);
	}
	return;
} /* qagp_ */

void qagse_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *epsabs, const sys_float *epsrel, 
				const long *limit, sys_float *result, sys_float *abserr, 
				long *neval, 
				long *ier, sys_float *alist__, sys_float *blist, 
				sys_float *rlist, sys_float *elist, 
				long *iord, long *last)
{
	/* System generated locals */
	int i__1, i__2;
	sys_float r__1, r__2;

	/* Local variables */
	long k;
	sys_float a1, a2, b1, b2;
	long id;
	sys_float area;
	sys_float dres;
	long ksgn, nres;
	sys_float area1, area2, area12, small=0.f, erro12;
	long ierro;
	sys_float defab1, defab2;
	long ktmin, nrmax;
	sys_float oflow, uflow;
	bool noext;
	long iroff1, iroff2, iroff3;
	sys_float res3la[3], error1, error2, rlist2[52];
	long numrl2;
	sys_float defabs, epmach, erlarg=0.f, abseps, correc=0.f, errbnd, resabs;
	long jupbnd;
	sys_float erlast, errmax;
	long maxerr;
	sys_float reseps;
	bool extrap;
	sys_float ertest=0.f, errsum;

	/* ***begin prologue  qagse */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             (end point) singularities, extrapolation, */
	/*             globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upperbound on the number of subintervals */
	/*                     in the partition of (a,b) */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                         = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more sub- */
	/*                             divisions by increasing the value of limit */
	/*                             (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. if */
	/*                             the position of a local difficulty can be */
	/*                             determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is detec- */
	/*                             ted, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour */
	/*                             occurs at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. */
	/*                             it is presumed that the requested */
	/*                             tolerance cannot be achieved, and that the */
	/*                             returned result is the best which can be */
	/*                             obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier. */
	/*                         = 6 the input is invalid, because */
	/*                             epsabs.le.0 and */
	/*                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28). */
	/*                             result, abserr, neval, last, rlist(1), */
	/*                             iord(1) and elist(1) are set to zero. */
	/*                             alist(1) and blist(1) are set to a and b */
	/*                             respectively. */

	/*            alist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the left end points */
	/*                     of the subintervals in the partition of the */
	/*                     given integration range (a,b) */

	/*            blist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the right end points */
	/*                     of the subintervals in the partition of the given */
	/*                     integration range (a,b) */

	/*            rlist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the integral */
	/*                     approximations on the subintervals */

	/*            elist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the moduli of the */
	/*                     absolute error estimates on the subintervals */

	/*            iord   - long */
	/*                     vector of dimension at least limit, the first k */
	/*                     elements of which are pointers to the */
	/*                     error estimates over the subintervals, */
	/*                     such that elist(iord(1)), ..., elist(iord(k)) */
	/*                     form a decreasing sequence, with k = last */
	/*                     if last.le.(limit/2+2), and k = limit+1-last */
	/*                     otherwise */

	/*            last   - long */
	/*                     number of subintervals actually produced in the */
	/*                     subdivision process */

	/* ***references  (none) */
	/* ***routines called  qelg,qk21,qpsrt,r1mach */
	/* ***end prologue  qagse */




	/*            the dimension of rlist2 is determined by the value of */
	/*            limexp in subroutine qelg (rlist2 should be of dimension */
	/*            (limexp+2) at least). */

	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           rlist2    - array of dimension at least limexp+2 */
	/*                       containing the part of the epsilon table */
	/*                       which is still needed for further computations */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest error */
	/*                       estimate */
	/*           errmax    - elist(maxerr) */
	/*           erlast    - error on the interval currently subdivided */
	/*                       (before that subdivision has taken place) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left interval */
	/*           *****2    - variable for the right interval */
	/*           last      - index for subdivision */
	/*           nres      - number of calls to the extrapolation routine */
	/*           numrl2    - number of elements currently in rlist2. if an */
	/*                       appropriate approximation to the compounded */
	/*                       integral has been obtained it is put in */
	/*                       rlist2(numrl2) after numrl2 has been increased */
	/*                       by one. */
	/*           small     - length of the smallest interval considered */
	/*                       up to now, multiplied by 1.5 */
	/*           erlarg    - sum of the errors over the intervals larger */
	/*                       than the smallest interval considered up to now */
	/*           extrap    - bool variable denoting that the routine */
	/*                       is attempting to perform extrapolation */
	/*                       i.e. before subdividing the smallest interval */
	/*                       we try to decrease the value of erlarg. */
	/*           noext     - bool variable denoting that extrapolation */
	/*                       is no longer allowed (true value) */

	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */
	/*           oflow is the largest positive magnitude. */

	/* ***first executable statement  qagse */
	/* Parameter adjustments */
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = r1mach(c__4);

	/*            test on validity of parameters */
	/*            ------------------------------ */
	*ier = 0;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	alist__[1] = *a;
	blist[1] = *b;
	rlist[1] = 0.f;
	elist[1] = 0.f;
	/* Computing MAX */
	r__1 = epmach * 50.f;
	if (*epsabs <= 0.f && *epsrel < max(r__1,5e-15f)) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}

	/*           first approximation to the integral */
	/*           ----------------------------------- */

	uflow = r1mach(c__1);
	oflow = r1mach(c__2);
	ierro = 0;
	qk21_(f, a, b, result, abserr, &defabs, &resabs);

	/*           test on accuracy. */

	dres = fabs(*result);
	/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dres;
	errbnd = max(r__1,r__2);
	*last = 1;
	rlist[1] = *result;
	elist[1] = *abserr;
	iord[1] = 1;
	if (*abserr <= epmach * 100.f * defabs && *abserr > errbnd) {
		*ier = 2;
	}
	if (*limit == 1) {
		*ier = 1;
	}
	if (*ier != 0 || ( *abserr <= errbnd && *abserr != resabs ) || 
		 *abserr == 0.f)
	{
		goto L140;
	}

	/*           initialization */
	/*           -------------- */

	rlist2[0] = *result;
	errmax = *abserr;
	maxerr = 1;
	area = *result;
	errsum = *abserr;
	*abserr = oflow;
	nrmax = 1;
	nres = 0;
	numrl2 = 2;
	ktmin = 0;
	extrap = false;
	noext = false;
	iroff1 = 0;
	iroff2 = 0;
	iroff3 = 0;
	ksgn = -1;
	if (dres >= (1.f - epmach * 50.f) * defabs) {
		ksgn = 1;
	}

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limit;
	for (*last = 2; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with the nrmax-th largest */
		/*           error estimate. */

		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
		a2 = b1;
		b2 = blist[maxerr];
		erlast = errmax;
		qk21_(f, &a1, &b1, &area1, &error1, &resabs, &defab1);
		qk21_(f, &a2, &b2, &area2, &error2, &resabs, &defab2);

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (defab1 == error1 || defab2 == error2) {
			goto L15;
		}
		if ((r__1 = rlist[maxerr] - area12, fabs(r__1)) > fabs(area12) * 
			 1e-5f || erro12 < errmax * .99f) {
			goto L10;
		}
		if (extrap) {
			++iroff2;
		}
		if (! extrap) {
			++iroff1;
		}
	L10:
		if (*last > 10 && erro12 > errmax) {
			++iroff3;
		}
	L15:
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(area);
		errbnd = max(r__1,r__2);

		/*           test for roundoff error and eventually */
		/*           set error flag. */

		if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
			*ier = 2;
		}
		if (iroff2 >= 5) {
			ierro = 3;
		}

		/*           set error flag in the case that the number of */
		/*           subintervals equals limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at a point of the integration range. */

		/* Computing MAX */
		r__1 = fabs(a1), r__2 = fabs(b2);
		if (max(r__1,r__2) <= (epmach * 100.f + 1.f) * (fabs(a2) + uflow * 
																		 1e3f)) {
			*ier = 4;
		}

		/*           append the newly-created intervals to the list. */

		if (error2 > error1) {
			goto L20;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L30;
	L20:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine qpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the */
		/*           subinterval with nrmax-th largest error estimate (to be */
		/*           bisected next). */

	L30:
		qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (errsum <= errbnd) {
			goto L115;
		}
		/* ***jump out of do-loop */
		if (*ier != 0) {
			goto L100;
		}
		if (*last == 2) {
			goto L80;
		}
		if (noext) {
			goto L90;
		}
		erlarg -= erlast;
		if ((r__1 = b1 - a1, fabs(r__1)) > small) {
			erlarg += erro12;
		}
		if (extrap) {
			goto L40;
		}

		/*           test whether the interval to be bisected next is the */
		/*           smallest interval. */

		if ((r__1 = blist[maxerr] - alist__[maxerr], fabs(r__1)) > small) {
			goto L90;
		}
		extrap = true;
		nrmax = 2;
	L40:
		if (ierro == 3 || erlarg <= ertest) {
			goto L60;
		}

		/*           the smallest interval has the largest error. */
		/*           before bisecting decrease the sum of the errors */
		/*           over the larger intervals (erlarg) and perform */
		/*           extrapolation. */

		id = nrmax;
		jupbnd = *last;
		if (*last > *limit / 2 + 2) {
			jupbnd = *limit + 3 - *last;
		}
		i__2 = jupbnd;
		for (k = id; k <= i__2; ++k) {
			maxerr = iord[nrmax];
			errmax = elist[maxerr];
			/* ***jump out of do-loop */
			if ((r__1 = blist[maxerr] - alist__[maxerr], fabs(r__1)) > small) 
			{
				goto L90;
			}
			++nrmax;
			/* L50: */
		}

		/*           perform extrapolation. */

	L60:
		++numrl2;
		rlist2[numrl2 - 1] = area;
		qelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
		++ktmin;
		if (ktmin > 5 && *abserr < errsum * .001f) {
			*ier = 5;
		}
		if (abseps >= *abserr) {
			goto L70;
		}
		ktmin = 0;
		*abserr = abseps;
		*result = reseps;
		correc = erlarg;
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(reseps);
		ertest = max(r__1,r__2);
		/* ***jump out of do-loop */
		if (*abserr <= ertest) {
			goto L100;
		}

		/*           prepare bisection of the smallest interval. */

	L70:
		if (numrl2 == 1) {
			noext = true;
		}
		if (*ier == 5) {
			goto L100;
		}
		maxerr = iord[1];
		errmax = elist[maxerr];
		nrmax = 1;
		extrap = false;
		small *= .5f;
		erlarg = errsum;
		goto L90;
	L80:
		small = (r__1 = *b - *a, fabs(r__1)) * .375f;
		erlarg = errsum;
		ertest = errbnd;
		rlist2[1] = area;
	L90:
		;
	}

	/*           set final result and error estimate. */
	/*           ------------------------------------ */

L100:
	if (*abserr == oflow) {
		goto L115;
	}
	if (*ier + ierro == 0) {
		goto L110;
	}
	if (ierro == 3) {
		*abserr += correc;
	}
	if (*ier == 0) {
		*ier = 3;
	}
	if (*result != 0.f && area != 0.f) {
		goto L105;
	}
	if (*abserr > errsum) {
		goto L115;
	}
	if (area == 0.f) {
		goto L130;
	}
	goto L110;
L105:
	if (*abserr / fabs(*result) > errsum / fabs(area)) {
		goto L115;
	}

	/*           test on divergence. */

L110:
	/* Computing MAX */
	r__1 = fabs(*result), r__2 = fabs(area);
	if (ksgn == -1 && max(r__1,r__2) <= defabs * .01f) {
		goto L130;
	}
	if (.01f > *result / area || *result / area > 100.f || errsum > fabs(area)
		) {
		*ier = 6;
	}
	goto L130;

	/*           compute global integral sum. */

L115:
	*result = 0.f;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L120: */
	}
	*abserr = errsum;
L130:
	if (*ier > 2) {
		--(*ier);
	}
L140:
	*neval = *last * 42 - 21;
L999:
	return;
} /* qagse_ */

void qags_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  const sys_float *epsabs, const sys_float *epsrel, 
			  sys_float *result, sys_float *abserr, long *neval, long *ier, 
			  const long *limit, const long *lenw, long *last, long *iwork, 
			  sys_float *work)
{
	long l1, l2, l3, lvl;

	/* ***begin prologue  qags */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a1 */
	/* ***keywords  automatic integrator, general-purpose, */
	/*             (end-point) singularities, extrapolation, */
	/*             globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & prog. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral  i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of a definite integral */
	/*        standard fortran subroutine */
	/*        sys_float version */


	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more sub- */
	/*                             divisions by increasing the value of limit */
	/*                             (and taking the according dimension */
	/*                             adjustments into account. however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. if */
	/*                             the position of a local difficulty can be */
	/*                             determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used, which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is detec- */
	/*                             ted, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour */
	/*                             occurs at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. it is presumed that */
	/*                             the requested tolerance cannot be */
	/*                             achieved, and that the returned result is */
	/*                             the best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28) */
	/*                             or limit.lt.1 or lenw.lt.limit*4. */
	/*                             result, abserr, neval, last are set to */
	/*                             zero.except when limit or lenw is invalid, */
	/*                             iwork(1), work(limit*2+1) and */
	/*                             work(limit*3+1) are set to zero, work(1) */
	/*                             is set to a and work(limit+1) to b. */

	/*         dimensioning parameters */
	/*            limit - long */
	/*                    dimensioning parameter for iwork */
	/*                    limit determines the maximum number of subintervals */
	/*                    in the partition of the given integration interval */
	/*                    (a,b), limit.ge.1. */
	/*                    if limit.lt.1, the routine will end with ier = 6. */

	/*            lenw  - long */
	/*                    dimensioning parameter for work */
	/*                    lenw must be at least limit*4. */
	/*                    if lenw.lt.limit*4, the routine will end */
	/*                    with ier = 6. */

	/*            last  - long */
	/*                    on return, last equals the number of subintervals */
	/*                    produced in the subdivision process, detemines the */
	/*                    number of significant elements actually in the work */
	/*                    arrays. */

	/*         work arrays */
	/*            iwork - long */
	/*                    vector of dimension at least limit, the first k */
	/*                    elements of which contain pointers */
	/*                    to the error estimates over the subintervals */
	/*                    such that work(limit*3+iwork(1)),... , */
	/*                    work(limit*3+iwork(k)) form a decreasing */
	/*                    sequence, with k = last if last.le.(limit/2+2), */
	/*                    and k = limit+1-last otherwise */

	/*            work  - sys_float */
	/*                    vector of dimension at least lenw */
	/*                    on return */
	/*                    work(1), ..., work(last) contain the left */
	/*                     end-points of the subintervals in the */
	/*                     partition of (a,b), */
	/*                    work(limit+1), ..., work(limit+last) contain */
	/*                     the right end-points, */
	/*                    work(limit*2+1), ..., work(limit*2+last) contain */
	/*                     the integral approximations over the subintervals, */
	/*                    work(limit*3+1), ..., work(limit*3+last) */
	/*                     contain the error estimates. */



	/* ***references  (none) */
	/* ***routines called  qagse,xerror */
	/* ***end prologue  qags */





	/*         check validity of limit and lenw. */

	/* ***first executable statement  qags */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	if (*limit < 1 || *lenw < *limit << 2) {
		goto L10;
	}

	/*         prepare call for qagse. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;

	qagse_(f, a, b, epsabs, epsrel, limit, result, abserr, neval, ier, &
			 work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from  qags", &c__26, ier, &lvl, 26);
	}
	return;
} /* qags_ */

void qawce_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *c__, const sys_float *epsabs,
				const sys_float *epsrel, const long *limit, 
				sys_float *result, sys_float *abserr, long *
				neval, long *ier, sys_float *alist__, sys_float *blist, 
				sys_float *rlist, sys_float *elist, long *iord, long *last)
{
	/* System generated locals */
	int i__1;
	sys_float r__1, r__2;

	/* Local variables */
	long k;
	sys_float a1, a2, b1, b2, aa, bb;
	long nev;
	sys_float area;
	sys_float area1, area2, area12, erro12;
	long krule, nrmax;
	sys_float uflow;
	long iroff1, iroff2;
	sys_float error1, error2, epmach, errbnd, errmax;
	long maxerr;
	sys_float errsum;

	/* ***begin prologue  qawce */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1,j4 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             cauchy principal value, clenshaw-curtis method */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***  purpose  the routine calculates an approximation result to a */
	/*              cauchy principal value i = integral of f*w over (a,b) */
	/*              (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying */
	/*              following claim for accuracy */
	/*              fabs(i-result).le.max(epsabs,epsrel*fabs(i)) */
	/* ***description */

	/*        computation of a cauchy principal value */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration */

	/*            c      - sys_float */
	/*                     parameter in the weight function, c.ne.a, c.ne.b */
	/*                     if c = a or c = b, the routine will end with */
	/*                     ier = 6. */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subintervals */
	/*                     in the partition of (a,b), limit.ge.1 */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more sub- */
	/*                             divisions by increasing the value of */
	/*                             limit. however, if this yields no */
	/*                             improvement it is advised to analyze the */
	/*                             the integrand, in order to determine the */
	/*                             the integration difficulties. if the */
	/*                             position of a local difficulty can be */
	/*                             determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling */
	/*                             appropriate integrators on the subranges. */
	/*                         = 2 the occurrence of roundoff error is detec- */
	/*                             ted, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour */
	/*                             occurs at some interior points of */
	/*                             the integration interval. */
	/*                         = 6 the input is invalid, because */
	/*                             c = a or c = b or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or limit.lt.1. */
	/*                             result, abserr, neval, rlist(1), elist(1), */
	/*                             iord(1) and last are set to zero. alist(1) */
	/*                             and blist(1) are set to a and b */
	/*                             respectively. */

	/*            alist   - sys_float */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the left */
	/*                      end points of the subintervals in the partition */
	/*                      of the given integration range (a,b) */

	/*            blist   - sys_float */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the right */
	/*                      end points of the subintervals in the partition */
	/*                      of the given integration range (a,b) */

	/*            rlist   - sys_float */
	/*                      vector of dimension at least limit, the first */
	/*                       last  elements of which are the integral */
	/*                      approximations on the subintervals */

	/*            elist   - sys_float */
	/*                      vector of dimension limit, the first  last */
	/*                      elements of which are the moduli of the absolute */
	/*                      error estimates on the subintervals */

	/*            iord    - long */
	/*                      vector of dimension at least limit, the first k */
	/*                      elements of which are pointers to the error */
	/*                      estimates over the subintervals, so that */
	/*                      elist(iord(1)), ..., elist(iord(k)) with k = last */
	/*                      if last.le.(limit/2+2), and k = limit+1-last */
	/*                      otherwise, form a decreasing sequence */

	/*            last    - long */
	/*                      number of subintervals actually produced in */
	/*                      the subdivision process */

	/* ***references  (none) */
	/* ***routines called  qc25c,qpsrt,r1mach */
	/* ***end prologue  qawce */




	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest */
	/*                       error estimate */
	/*           errmax    - elist(maxerr) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */


	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qawce */
	/* Parameter adjustments */
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);


	/*           test on validity of parameters */
	/*           ------------------------------ */

	*ier = 6;
	*neval = 0;
	*last = 0;
	alist__[1] = *a;
	blist[1] = *b;
	rlist[1] = 0.f;
	elist[1] = 0.f;
	iord[1] = 0;
	*result = 0.f;
	*abserr = 0.f;
	/* Computing MAX */
	r__1 = epmach * 50.f;
	if (*c__ == *a || *c__ == *b || 
		 ( *epsabs <= 0.f && *epsrel < max(r__1, 5e-15f))) {
		goto L999;
	}

	/*           first approximation to the integral */
	/*           ----------------------------------- */

	aa = *a;
	bb = *b;
	if (*a <= *b) {
		goto L10;
	}
	aa = *b;
	bb = *a;
L10:
	*ier = 0;
	krule = 1;
	qc25c_(f, &aa, &bb, c__, result, abserr, &krule, neval);
	*last = 1;
	rlist[1] = *result;
	elist[1] = *abserr;
	iord[1] = 1;
	alist__[1] = *a;
	blist[1] = *b;

	/*           test on accuracy */

	/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * fabs(*result);
	errbnd = max(r__1,r__2);
	if (*limit == 1) {
		*ier = 1;
	}
	/* Computing MIN */
	r__1 = fabs(*result) * .01f;
	if (*abserr < min(r__1,errbnd) || *ier == 1) {
		goto L70;
	}

	/*           initialization */
	/*           -------------- */

	alist__[1] = aa;
	blist[1] = bb;
	rlist[1] = *result;
	errmax = *abserr;
	maxerr = 1;
	area = *result;
	errsum = *abserr;
	nrmax = 1;
	iroff1 = 0;
	iroff2 = 0;

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limit;
	for (*last = 2; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with nrmax-th largest */
		/*           error estimate. */

		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
		b2 = blist[maxerr];
		if (*c__ <= b1 && *c__ > a1) {
			b1 = (*c__ + b2) * .5f;
		}
		if (*c__ > b1 && *c__ < b2) {
			b1 = (a1 + *c__) * .5f;
		}
		a2 = b1;
		krule = 2;
		qc25c_(f, &a1, &b1, c__, &area1, &error1, &krule, &nev);
		*neval += nev;
		qc25c_(f, &a2, &b2, c__, &area2, &error2, &krule, &nev);
		*neval += nev;

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if ((r__1 = rlist[maxerr] - area12, fabs(r__1)) < fabs(area12) * 
			 1e-5f && erro12 >= errmax * .99f && krule == 0) {
			++iroff1;
		}
		if (*last > 10 && erro12 > errmax && krule == 0) {
			++iroff2;
		}
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(area);
		errbnd = max(r__1,r__2);
		if (errsum <= errbnd) {
			goto L15;
		}

		/*           test for roundoff error and eventually */
		/*           set error flag. */

		if (iroff1 >= 6 && iroff2 > 20) {
			*ier = 2;
		}

		/*           set error flag in the case that number of interval */
		/*           bisections exceeds limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at a point of the integration range. */

		/* Computing MAX */
		r__1 = fabs(a1), r__2 = fabs(b2);
		if (max(r__1,r__2) <= (epmach * 100.f + 1.f) * (fabs(a2) + uflow * 
																		 1e3f)) {
			*ier = 3;
		}

		/*           append the newly-created intervals to the list. */

	L15:
		if (error2 > error1) {
			goto L20;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L30;
	L20:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine qpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the */
		/*           subinterval with nrmax-th largest error estimate (to be */
		/*           bisected next). */

	L30:
		qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (*ier != 0 || errsum <= errbnd) {
			goto L50;
		}
		/* L40: */
	}

	/*           compute final result. */
	/*           --------------------- */

L50:
	*result = 0.f;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L60: */
	}
	*abserr = errsum;
L70:
	if (aa == *b) {
		*result = -(*result);
	}
L999:
	return;
} /* qawce_ */

void qawc_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  const sys_float *c__, const sys_float *epsabs, 
			  const sys_float *epsrel, sys_float *result, 
			  sys_float *abserr, long *neval, long *ier, 
			  long *limit, const long *lenw, long *last, long *iwork, 
			  sys_float *work)
{
	long l1, l2, l3, lvl;

	/* ***begin prologue  qawc */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1,j4 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             cauchy principal value, */
	/*             clenshaw-curtis, globally adaptive */
	/* ***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a */
	/*            cauchy principal value i = integral of f*w over (a,b) */
	/*            (w(x) = 1/((x-c), c.ne.a, c.ne.b), hopefully satisfying */
	/*            following claim for accuracy */
	/*            fabs(i-result).le.max(epsabe,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of a cauchy principal value */
	/*        standard fortran subroutine */
	/*        sys_float version */


	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     under limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration */

	/*            c      - parameter in the weight function, c.ne.a, c.ne.b. */
	/*                     if c = a or c = b, the routine will end with */
	/*                     ier = 6 . */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate or the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more sub- */
	/*                             divisions by increasing the value of limit */
	/*                             (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand in order to */
	/*                             determine the integration difficulties. */
	/*                             if the position of a local difficulty */
	/*                             can be determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling */
	/*                             appropriate integrators on the subranges. */
	/*                         = 2 the occurrence of roundoff error is detec- */
	/*                             ted, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 6 the input is invalid, because */
	/*                             c = a or c = b or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or limit.lt.1 or lenw.lt.limit*4. */
	/*                             result, abserr, neval, last are set to */
	/*                             zero. exept when lenw or limit is invalid, */
	/*                             iwork(1), work(limit*2+1) and */
	/*                             work(limit*3+1) are set to zero, work(1) */
	/*                             is set to a and work(limit+1) to b. */

	/*         dimensioning parameters */
	/*            limit - long */
	/*                    dimensioning parameter for iwork */
	/*                    limit determines the maximum number of subintervals */
	/*                    in the partition of the given integration interval */
	/*                    (a,b), limit.ge.1. */
	/*                    if limit.lt.1, the routine will end with ier = 6. */

	/*           lenw   - long */
	/*                    dimensioning parameter for work */
	/*                    lenw must be at least limit*4. */
	/*                    if lenw.lt.limit*4, the routine will end with */
	/*                    ier = 6. */

	/*            last  - long */
	/*                    on return, last equals the number of subintervals */
	/*                    produced in the subdivision process, which */
	/*                    determines the number of significant elements */
	/*                    actually in the work arrays. */

	/*         work arrays */
	/*            iwork - long */
	/*                    vector of dimension at least limit, the first k */
	/*                    elements of which contain pointers */
	/*                    to the error estimates over the subintervals, */
	/*                    such that work(limit*3+iwork(1)), ... , */
	/*                    work(limit*3+iwork(k)) form a decreasing */
	/*                    sequence, with k = last if last.le.(limit/2+2), */
	/*                    and k = limit+1-last otherwise */

	/*            work  - sys_float */
	/*                    vector of dimension at least lenw */
	/*                    on return */
	/*                    work(1), ..., work(last) contain the left */
	/*                     end points of the subintervals in the */
	/*                     partition of (a,b), */
	/*                    work(limit+1), ..., work(limit+last) contain */
	/*                     the right end points, */
	/*                    work(limit*2+1), ..., work(limit*2+last) contain */
	/*                     the integral approximations over the subintervals, */
	/*                    work(limit*3+1), ..., work(limit*3+last) */
	/*                     contain the error estimates. */

	/* ***references  (none) */
	/* ***routines called  qawce,xerror */
	/* ***end prologue  qawc */




	/*         check validity of limit and lenw. */

	/* ***first executable statement  qawc */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	if (*limit < 1 || *lenw < *limit << 2) {
		goto L10;
	}

	/*         prepare call for qawce. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;
	qawce_(f, a, b, c__, epsabs, epsrel, limit, result, abserr, neval, 
			 ier, &work[1], &work[l1], &work[l2], &work[l3], &iwork[1], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from  qawc", &c__26, ier, &lvl, 26);
	}
	return;
} /* qawc_ */

void qawfe_(const E_fp& f, const sys_float *a, const sys_float *omega, 
				const long *integr, 
				const sys_float *epsabs, const long *limlst, const long *limit, 
				const long *maxp1, sys_float *result, sys_float *abserr, 
				long *neval, long *ier, sys_float *rslst, sys_float *erlst, 
				long *ierlst, long *lst, sys_float *alist__, sys_float *blist, 
				sys_float *rlist, sys_float *elist, long *iord, long *nnlog, 
				sys_float *chebmo)
{
	/* Initialized data */

	sys_float p = .9f;
	sys_float pi = 3.1415926535897932f;

	/* System generated locals */
	int chebmo_dim1, chebmo_offset, i__1;
	sys_float r__1, r__2;

	/* Local variables */
	long l;
	sys_float c1, c2, p1, dl, ep;
	long ll=0;
	sys_float drl=0.f, eps;
	long nev;
	sys_float fact, epsa;
	long last, nres;
	sys_float psum[52];
	sys_float cycle;
	long ktmin;
	sys_float uflow;
	sys_float res3la[3];
	long numrl2;
	sys_float abseps, correc;
	long momcom;
	sys_float reseps, errsum;

	/* ***begin prologue  qawfe */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a3a1 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             fourier integrals, */
	/*             integration between zeros with dqawoe, */
	/*             convergence acceleration with dqelg */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           dedoncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a */
	/*            given fourier integal */
	/*            i = integral of f(x)*w(x) over (a,infinity) */
	/*             where w(x) = cos(omega*x) or w(x) = sin(omega*x), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.epsabs. */
	/* ***description */

	/*        computation of fourier integrals */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to */
	/*                     be declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            omega  - sys_float */
	/*                     parameter in the weight function */

	/*            integr - long */
	/*                     indicates which weight function is used */
	/*                     integr = 1      w(x) = cos(omega*x) */
	/*                     integr = 2      w(x) = sin(omega*x) */
	/*                     if integr.ne.1.and.integr.ne.2, the routine will */
	/*                     end with ier = 6. */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested, epsabs.gt.0 */
	/*                     if epsabs.le.0, the routine will end with ier = 6. */

	/*            limlst - long */
	/*                     limlst gives an upper bound on the number of */
	/*                     cycles, limlst.ge.1. */
	/*                     if limlst.lt.3, the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subintervals */
	/*                     allowed in the partition of each cycle, limit.ge.1 */
	/*                     each cycle, limit.ge.1. */

	/*            maxp1  - long */
	/*                     gives an upper bound on the number of */
	/*                     chebyshev moments which can be stored, i.e. */
	/*                     for the intervals of lengths fabs(b-a)*2**(-l), */
	/*                     l=0,1, ..., maxp1-2, maxp1.ge.1 */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral x */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - ier = 0 normal and reliable termination of */
	/*                             the routine. it is assumed that the */
	/*                             requested accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine. the */
	/*                             estimates for integral and error are less */
	/*                             reliable. it is assumed that the requested */
	/*                             accuracy has not been achieved. */
	/*            error messages */
	/*                    if omega.ne.0 */
	/*                     ier = 1 maximum number of  cycles  allowed */
	/*                             has been achieved., i.e. of subintervals */
	/*                             (a+(k-1)c,a+kc) where */
	/*                             c = (2*int(fabs(omega))+1)*pi/fabs(omega), */
	/*                             for k = 1, 2, ..., lst. */
	/*                             one can allow more cycles by increasing */
	/*                             the value of limlst (and taking the */
	/*                             according dimension adjustments into */
	/*                             account). */
	/*                             examine the array iwork which contains */
	/*                             the error flags on the cycles, in order to */
	/*                             look for eventual local integration */
	/*                             difficulties. if the position of a local */
	/*                             difficulty can be determined (e.g. */
	/*                             singularity, discontinuity within the */
	/*                             interval) one will probably gain from */
	/*                             splitting up the interval at this polong */
	/*                             and calling appropriate integrators on */
	/*                             the subranges. */
	/*                         = 4 the extrapolation table constructed for */
	/*                             convergence acceleration of the series */
	/*                             formed by the integral contributions over */
	/*                             the cycles, does not converge to within */
	/*                             the requested accuracy. as in the case of */
	/*                             ier = 1, it is advised to examine the */
	/*                             array iwork which contains the error */
	/*                             flags on the cycles. */
	/*                         = 6 the input is invalid because */
	/*                             (integr.ne.1 and integr.ne.2) or */
	/*                              epsabs.le.0 or limlst.lt.3. */
	/*                              result, abserr, neval, lst are set */
	/*                              to zero. */
	/*                         = 7 bad integrand behaviour occurs within one */
	/*                             or more of the cycles. location and type */
	/*                             of the difficulty involved can be */
	/*                             determined from the vector ierlst. here */
	/*                             lst is the number of cycles actually */
	/*                             needed (see below). */
	/*                             ierlst(k) = 1 the maximum number of */
	/*                                           subdivisions (= limit) has */
	/*                                           been achieved on the k th */
	/*                                           cycle. */
	/*                                       = 2 occurrence of roundoff error */
	/*                                           is detected and prevents the */
	/*                                           tolerance imposed on the */
	/*                                           k th cycle, from being */
	/*                                           achieved. */
	/*                                       = 3 extremely bad integrand */
	/*                                           behaviour occurs at some */
	/*                                           points of the k th cycle. */
	/*                                       = 4 the integration procedure */
	/*                                           over the k th cycle does */
	/*                                           not converge (to within the */
	/*                                           required accuracy) due to */
	/*                                           roundoff in the */
	/*                                           extrapolation procedure */
	/*                                           invoked on this cycle. it */
	/*                                           is assumed that the result */
	/*                                           on this interval is the */
	/*                                           best which can be obtained. */
	/*                                       = 5 the integral over the k th */
	/*                                           cycle is probably divergent */
	/*                                           or slowly convergent. it */
	/*                                           must be noted that */
	/*                                           divergence can occur with */
	/*                                           any other value of */
	/*                                           ierlst(k). */
	/*                    if omega = 0 and integr = 1, */
	/*                    the integral is calculated by means of dqagie */
	/*                    and ier = ierlst(1) (with meaning as described */
	/*                    for ierlst(k), k = 1). */

	/*            rslst  - sys_float */
	/*                     vector of dimension at least limlst */
	/*                     rslst(k) contains the integral contribution */
	/*                     over the interval (a+(k-1)c,a+kc) where */
	/*                     c = (2*int(fabs(omega))+1)*pi/fabs(omega), */
	/*                     k = 1, 2, ..., lst. */
	/*                     note that, if omega = 0, rslst(1) contains */
	/*                     the value of the integral over (a,infinity). */

	/*            erlst  - sys_float */
	/*                     vector of dimension at least limlst */
	/*                     erlst(k) contains the error estimate corresponding */
	/*                     with rslst(k). */

	/*            ierlst - long */
	/*                     vector of dimension at least limlst */
	/*                     ierlst(k) contains the error flag corresponding */
	/*                     with rslst(k). for the meaning of the local error */
	/*                     flags see description of output parameter ier. */

	/*            lst    - long */
	/*                     number of subintervals needed for the integration */
	/*                     if omega = 0 then lst is set to 1. */

	/*            alist, blist, rlist, elist - sys_float */
	/*                     vector of dimension at least limit, */

	/*            iord, nnlog - long */
	/*                     vector of dimension at least limit, providing */
	/*                     space for the quantities needed in the subdivision */
	/*                     process of each cycle */

	/*            chebmo - sys_float */
	/*                     array of dimension at least (maxp1,25), providing */
	/*                     space for the chebyshev moments needed within the */
	/*                     cycles */

	/* ***references  (none) */
	/* ***routines called  qagie,qawoe,qelg,r1mach */
	/* ***end prologue  qawfe */





	/*            the dimension of  psum  is determined by the value of */
	/*            limexp in subroutine qelg (psum must be */
	/*            of dimension (limexp+2) at least). */

	/*           list of major variables */
	/*           ----------------------- */

	/*           c1, c2    - end points of subinterval (of length */
	/*                       cycle) */
	/*           cycle     - (2*int(fabs(omega))+1)*pi/fabs(omega) */
	/*           psum      - vector of dimension at least (limexp+2) */
	/*                       (see routine qelg) */
	/*                       psum contains the part of the epsilon */
	/*                       table which is still needed for further */
	/*                       computations. */
	/*                       each element of psum is a partial sum of */
	/*                       the series which should sum to the value of */
	/*                       the integral. */
	/*           errsum    - sum of error estimates over the */
	/*                       subintervals, calculated cumulatively */
	/*           epsa      - absolute tolerance requested over current */
	/*                       subinterval */
	/*           chebmo    - array containing the modified chebyshev */
	/*                       moments (see also routine qc25f) */

	/* Parameter adjustments */
	--ierlst;
	--erlst;
	--rslst;
	--nnlog;
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;
	chebmo_dim1 = *maxp1;
	chebmo_offset = 1 + chebmo_dim1;
	chebmo -= chebmo_offset;

	/* Function Body */

	/*           test on validity of parameters */
	/*           ------------------------------ */

	/* ***first executable statement  qawfe */
	*result = 0.f;
	*abserr = 0.f;
	*neval = 0;
	*lst = 0;
	*ier = 0;
	if ( ( *integr != 1 && *integr != 2 ) || *epsabs <= 0.f || *limlst < 3) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}
	if (*omega != 0.f) {
		goto L10;
	}

	/*           integration by qagie if omega is zero */
	/*           -------------------------------------- */

	if (*integr == 1) {
		qagie_(f, &c_b390, &c__1, epsabs, &c_b390, limit, result, 
				 abserr, neval, ier, &alist__[1], &blist[1], &rlist[1], &elist[
					 1], &iord[1], &last);
	}
	rslst[1] = *result;
	erlst[1] = *abserr;
	ierlst[1] = *ier;
	*lst = 1;
	goto L999;

	/*           initializations */
	/*           --------------- */

L10:
	l = fabs(*omega);
	dl = (sys_float) ((l << 1) + 1);
	cycle = dl * pi / fabs(*omega);
	*ier = 0;
	ktmin = 0;
	*neval = 0;
	numrl2 = 0;
	nres = 0;
	c1 = *a;
	c2 = cycle + *a;
	p1 = 1.f - p;
	eps = *epsabs;
	uflow = r1mach(c__1);
	if (*epsabs > uflow / p1) {
		eps = *epsabs * p1;
	}
	ep = eps;
	fact = 1.f;
	correc = 0.f;
	*abserr = 0.f;
	errsum = 0.f;

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limlst;
	for (*lst = 1; *lst <= i__1; ++(*lst)) {

		/*           integrate over current subinterval. */

		epsa = eps * fact;
		qawoe_(f, &c1, &c2, omega, integr, &epsa, &c_b390, limit, lst, 
				 maxp1, &rslst[*lst], &erlst[*lst], &nev, &ierlst[*lst], &last,
				 &alist__[1], &blist[1], &rlist[1], &elist[1], &iord[1], &
				 nnlog[1], &momcom, &chebmo[chebmo_offset]);
		*neval += nev;
		fact *= p;
		errsum += erlst[*lst];
		drl = (r__1 = rslst[*lst], fabs(r__1)) * 50.f;

		/*           test on accuracy with partial sum */

		if (errsum + drl <= *epsabs && *lst >= 6) {
			goto L80;
		}
		/* Computing MAX */
		r__1 = correc, r__2 = erlst[*lst];
		correc = max(r__1,r__2);
		if (ierlst[*lst] != 0) {
			/* Computing MAX */
			r__1 = ep, r__2 = correc * p1;
			eps = max(r__1,r__2);
		}
		if (ierlst[*lst] != 0) {
			*ier = 7;
		}
		if (*ier == 7 && errsum + drl <= correc * 10.f && *lst > 5) {
			goto L80;
		}
		++numrl2;
		if (*lst > 1) {
			goto L20;
		}
		psum[0] = rslst[1];
		goto L40;
	L20:
		psum[numrl2 - 1] = psum[ll - 1] + rslst[*lst];
		if (*lst == 2) {
			goto L40;
		}

		/*           test on maximum number of subintervals */

		if (*lst == *limlst) {
			*ier = 1;
		}

		/*           perform new extrapolation */

		qelg_(&numrl2, psum, &reseps, &abseps, res3la, &nres);

		/*           test whether extrapolated result is influenced by */
		/*           roundoff */

		++ktmin;
		if (ktmin >= 15 && *abserr <= (errsum + drl) * .001f) {
			*ier = 4;
		}
		if (abseps > *abserr && *lst != 3) {
			goto L30;
		}
		*abserr = abseps;
		*result = reseps;
		ktmin = 0;

		/*           if ier is not 0, check whether direct result (partial */
		/*           sum) or extrapolated result yields the best integral */
		/*           approximation */

		if (*abserr + correc * 10.f <= *epsabs || 
			 ( *abserr <= *epsabs && correc * 10.f >= *epsabs) ) {
			goto L60;
		}
	L30:
		if (*ier != 0 && *ier != 7) {
			goto L60;
		}
	L40:
		ll = numrl2;
		c1 = c2;
		c2 += cycle;
		/* L50: */
	}

	/*         set final result and error estimate */
	/*         ----------------------------------- */

L60:
	*abserr += correc * 10.f;
	if (*ier == 0) {
		goto L999;
	}
	if (*result != 0.f && psum[numrl2 - 1] != 0.f) {
		goto L70;
	}
	if (*abserr > errsum) {
		goto L80;
	}
	if (psum[numrl2 - 1] == 0.f) {
		goto L999;
	}
L70:
	if (*abserr / fabs(*result) > (errsum + drl) / (r__1 = psum[numrl2 - 1], 
																	fabs(r__1))) {
		goto L80;
	}
	if (*ier >= 1 && *ier != 7) {
		*abserr += drl;
	}
	goto L999;
L80:
	*result = psum[numrl2 - 1];
	*abserr = errsum + drl;
L999:
	return;
} /* qawfe_ */

void qawf_(const E_fp& f, const sys_float *a, const sys_float *omega, 
			  const long *integr, 
			  const sys_float *epsabs, sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, const long *limlst, long *lst, 
			  const long *leniw, const long *maxp1, 
			  const long *lenw, long *iwork, sys_float *work)
{
	long l1, l2, l3, l4, l5, l6, ll2, lvl;
	long limit;

	/* ***begin prologue  qawf */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a3a1 */
	/* ***keywords  automatic integrator, special-purpose,fourier */
	/*             integral, integration between zeros with dqawoe, */
	/*             convergence acceleration with dqext */
	/* ***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            fourier integral */
	/*            i = integral of f(x)*w(x) over (a,infinity) */
	/*            where w(x) = cos(omega*x) or w(x) = sin(omega*x). */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.epsabs. */
	/* ***description */

	/*        computation of fourier integrals */
	/*        standard fortran subroutine */
	/*        sys_float version */


	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            omega  - sys_float */
	/*                     parameter in the integrand weight function */

	/*            integr - long */
	/*                     indicates which of the weight functions is used */
	/*                     integr = 1      w(x) = cos(omega*x) */
	/*                     integr = 2      w(x) = sin(omega*x) */
	/*                     if integr.ne.1.and.integr.ne.2, the routine */
	/*                     will end with ier = 6. */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested, epsabs.gt.0. */
	/*                     if epsabs.le.0, the routine will end with ier = 6. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine. */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                    if omega.ne.0 */
	/*                     ier = 1 maximum number of cycles allowed */
	/*                             has been achieved, i.e. of subintervals */
	/*                             (a+(k-1)c,a+kc) where */
	/*                             c = (2*int(fabs(omega))+1)*pi/fabs(omega), */
	/*                             for k = 1, 2, ..., lst. */
	/*                             one can allow more cycles by increasing */
	/*                             the value of limlst (and taking the */
	/*                             according dimension adjustments into */
	/*                             account). examine the array iwork which */
	/*                             contains the error flags on the cycles, in */
	/*                             order to look for eventual local */
	/*                             integration difficulties. */
	/*                             if the position of a local difficulty */
	/*                             can be determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling */
	/*                             appropriate integrators on the subranges. */
	/*                         = 4 the extrapolation table constructed for */
	/*                             convergence accelaration of the series */
	/*                             formed by the integral contributions over */
	/*                             the cycles, does not converge to within */
	/*                             the requested accuracy. */
	/*                             as in the case of ier = 1, it is advised */
	/*                             to examine the array iwork which contains */
	/*                             the error flags on the cycles. */
	/*                         = 6 the input is invalid because */
	/*                             (integr.ne.1 and integr.ne.2) or */
	/*                              epsabs.le.0 or limlst.lt.1 or */
	/*                              leniw.lt.(limlst+2) or maxp1.lt.1 or */
	/*                              lenw.lt.(leniw*2+maxp1*25). */
	/*                              result, abserr, neval, lst are set to */
	/*                              zero. */
	/*                         = 7 bad integrand behaviour occurs within */
	/*                             one or more of the cycles. location and */
	/*                             type of the difficulty involved can be */
	/*                             determined from the first lst elements of */
	/*                             vector iwork.  here lst is the number of */
	/*                             cycles actually needed (see below). */
	/*                             iwork(k) = 1 the maximum number of */
	/*                                          subdivisions (=(leniw-limlst) */
	/*                                          /2) has been achieved on the */
	/*                                          k th cycle. */
	/*                                      = 2 occurrence of roundoff error */
	/*                                          is detected and prevents the */
	/*                                          tolerance imposed on the k th */
	/*                                          cycle, from being achieved */
	/*                                          on this cycle. */
	/*                                      = 3 extremely bad integrand */
	/*                                          behaviour occurs at some */
	/*                                          points of the k th cycle. */
	/*                                      = 4 the integration procedure */
	/*                                          over the k th cycle does */
	/*                                          not converge (to within the */
	/*                                          required accuracy) due to */
	/*                                          roundoff in the extrapolation */
	/*                                          procedure invoked on this */
	/*                                          cycle. it is assumed that the */
	/*                                          result on this interval is */
	/*                                          the best which can be */
	/*                                          obtained. */
	/*                                      = 5 the integral over the k th */
	/*                                          cycle is probably divergent */
	/*                                          or slowly convergent. it must */
	/*                                          be noted that divergence can */
	/*                                          occur with any other value of */
	/*                                          iwork(k). */
	/*                    if omega = 0 and integr = 1, */
	/*                    the integral is calculated by means of dqagie, */
	/*                    and ier = iwork(1) (with meaning as described */
	/*                    for iwork(k),k = 1). */

	/*         dimensioning parameters */
	/*            limlst - long */
	/*                     limlst gives an upper bound on the number of */
	/*                     cycles, limlst.ge.3. */
	/*                     if limlst.lt.3, the routine will end with ier = 6. */

	/*            lst    - long */
	/*                     on return, lst indicates the number of cycles */
	/*                     actually needed for the integration. */
	/*                     if omega = 0, then lst is set to 1. */

	/*            leniw  - long */
	/*                     dimensioning parameter for iwork. on entry, */
	/*                     (leniw-limlst)/2 equals the maximum number of */
	/*                     subintervals allowed in the partition of each */
	/*                     cycle, leniw.ge.(limlst+2). */
	/*                     if leniw.lt.(limlst+2), the routine will end with */
	/*                     ier = 6. */

	/*            maxp1  - long */
	/*                     maxp1 gives an upper bound on the number of */
	/*                     chebyshev moments which can be stored, i.e. for */
	/*                     the intervals of lengths fabs(b-a)*2**(-l), */
	/*                     l = 0,1, ..., maxp1-2, maxp1.ge.1. */
	/*                     if maxp1.lt.1, the routine will end with ier = 6. */
	/*            lenw   - long */
	/*                     dimensioning parameter for work */
	/*                     lenw must be at least leniw*2+maxp1*25. */
	/*                     if lenw.lt.(leniw*2+maxp1*25), the routine will */
	/*                     end with ier = 6. */

	/*         work arrays */
	/*            iwork  - long */
	/*                     vector of dimension at least leniw */
	/*                     on return, iwork(k) for k = 1, 2, ..., lst */
	/*                     contain the error flags on the cycles. */

	/*            work   - sys_float */
	/*                     vector of dimension at least */
	/*                     on return, */
	/*                     work(1), ..., work(lst) contain the integral */
	/*                      approximations over the cycles, */
	/*                     work(limlst+1), ..., work(limlst+lst) contain */
	/*                      the error extimates over the cycles. */
	/*                     further elements of work have no specific */
	/*                     meaning for the user. */

	/* ***references  (none) */
	/* ***routines called  qawfe,xerror */
	/* ***end prologue  qawf */




	/*         check validity of limlst, leniw, maxp1 and lenw. */

	/* ***first executable statement  qawf */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*result = 0.f;
	*abserr = 0.f;
	if (*limlst < 3 || *leniw < *limlst + 2 || *maxp1 < 1 || *lenw < (*leniw 
																							<< 1) + *maxp1 * 25) {
		goto L10;
	}

	/*         prepare call for qawfe */

	limit = (*leniw - *limlst) / 2;
	l1 = *limlst + 1;
	l2 = *limlst + l1;
	l3 = limit + l2;
	l4 = limit + l3;
	l5 = limit + l4;
	l6 = limit + l5;
	ll2 = limit + l1;
	qawfe_(f, a, omega, integr, epsabs, limlst, &limit, maxp1, result, 
			 abserr, neval, ier, &work[1], &work[l1], &iwork[1], lst, &work[l2]
			 , &work[l3], &work[l4], &work[l5], &iwork[l1], &iwork[ll2], &work[
				 l6]);

	/*         call error handler if necessary */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from  qawf", &c__26, ier, &lvl, 26);
	}
	return;
} /* qawf_ */

void qawoe_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *omega, const long *integr, 
				const sys_float *epsabs, const sys_float *epsrel, 
				const long *limit, const long *icall, const long *maxp1, 
				sys_float *result, sys_float *abserr, long *neval, long *
				ier, long *last, sys_float *alist__, sys_float *blist, sys_float *rlist, sys_float *
				elist, long *iord, long *nnlog, long *momcom, sys_float *chebmo)
{
	/* System generated locals */
	int chebmo_dim1, chebmo_offset, i__1, i__2;
	sys_float r__1, r__2;

	/* Local variables */
	long k;
	sys_float a1, a2, b1, b2;
	long id, nev;
	sys_float area;
	sys_float dres;
	long ksgn, nres;
	sys_float area1, area2, area12, small, erro12, defab1, defab2, width;
	long ierro;
	sys_float oflow;
	long ktmin, nrmax, nrmom;
	sys_float uflow;
	bool noext;
	long iroff1, iroff2, iroff3;
	sys_float res3la[3], error1, error2, rlist2[52];
	long numrl2;
	sys_float defabs, domega, epmach, erlarg=0.f, abseps, correc=0.f, errbnd, 
		resabs;
	long jupbnd;
	bool extall;
	sys_float erlast, errmax;
	long maxerr;
	sys_float reseps;
	bool extrap;
	sys_float ertest=0.f, errsum;

	/* ***begin prologue  qawoe */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             integrand with oscillatory cos or sin factor, */
	/*             clenshaw-curtis method, (end point) singularities, */
	/*             extrapolation, globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral */
	/*            i = integral of f(x)*w(x) over (a,b) */
	/*            where w(x) = cos(omega*x) or w(x) = sin(omega*x), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of oscillatory integrals */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration */

	/*            omega  - sys_float */
	/*                     parameter in the integrand weight function */

	/*            integr - long */
	/*                     indicates which of the weight functions is to be */
	/*                     used */
	/*                     integr = 1      w(x) = cos(omega*x) */
	/*                     integr = 2      w(x) = sin(omega*x) */
	/*                     if integr.ne.1 and integr.ne.2, the routine */
	/*                     will end with ier = 6. */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subdivisions */
	/*                     in the partition of (a,b), limit.ge.1. */

	/*            icall  - long */
	/*                     if dqawoe is to be used only once, icall must */
	/*                     be set to 1.  assume that during this call, the */
	/*                     chebyshev moments (for clenshaw-curtis integration */
	/*                     of degree 24) have been computed for intervals of */
	/*                     lenghts (fabs(b-a))*2**(-l), l=0,1,2,...momcom-1. */
	/*                     if icall.gt.1 this means that dqawoe has been */
	/*                     called twice or more on intervals of the same */
	/*                     length fabs(b-a). the chebyshev moments already */
	/*                     computed are then re-used in subsequent calls. */
	/*                     if icall.lt.1, the routine will end with ier = 6. */

	/*            maxp1  - long */
	/*                     gives an upper bound on the number of chebyshev */
	/*                     moments which can be stored, i.e. for the */
	/*                     intervals of lenghts fabs(b-a)*2**(-l), */
	/*                     l=0,1, ..., maxp1-2, maxp1.ge.1. */
	/*                     if maxp1.lt.1, the routine will end with ier = 6. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the */
	/*                             requested accuracy has been achieved. */
	/*                   - ier.gt.0 abnormal termination of the routine. */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand, in order to */
	/*                             determine the integration difficulties. */
	/*                             if the position of a local difficulty can */
	/*                             be determined (e.g. singularity, */
	/*                             discontinuity within the interval) one */
	/*                             will probably gain from splitting up the */
	/*                             interval at this point and calling the */
	/*                             integrator on the subranges. if possible, */
	/*                             an appropriate special-purpose integrator */
	/*                             should be used which is designed for */
	/*                             handling the type of difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. */
	/*                             it is presumed that the requested */
	/*                             tolerance cannot be achieved due to */
	/*                             roundoff in the extrapolation table, */
	/*                             and that the returned result is the */
	/*                             best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier.gt.0. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or (integr.ne.1 and integr.ne.2) or */
	/*                             icall.lt.1 or maxp1.lt.1. */
	/*                             result, abserr, neval, last, rlist(1), */
	/*                             elist(1), iord(1) and nnlog(1) are set */
	/*                             to zero. alist(1) and blist(1) are set */
	/*                             to a and b respectively. */

	/*            last  -  long */
	/*                     on return, last equals the number of */
	/*                     subintervals produces in the subdivision */
	/*                     process, which determines the number of */
	/*                     significant elements actually in the */
	/*                     work arrays. */
	/*            alist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the left */
	/*                     end points of the subintervals in the partition */
	/*                     of the given integration range (a,b) */

	/*            blist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the right */
	/*                     end points of the subintervals in the partition */
	/*                     of the given integration range (a,b) */

	/*            rlist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the integral */
	/*                     approximations on the subintervals */

	/*            elist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the moduli of the */
	/*                     absolute error estimates on the subintervals */

	/*            iord   - long */
	/*                     vector of dimension at least limit, the first k */
	/*                     elements of which are pointers to the error */
	/*                     estimates over the subintervals, */
	/*                     such that elist(iord(1)), ..., */
	/*                     elist(iord(k)) form a decreasing sequence, with */
	/*                     k = last if last.le.(limit/2+2), and */
	/*                     k = limit+1-last otherwise. */

	/*            nnlog  - long */
	/*                     vector of dimension at least limit, containing the */
	/*                     subdivision levels of the subintervals, i.e. */
	/*                     iwork(i) = l means that the subinterval */
	/*                     numbered i is of length fabs(b-a)*2**(1-l) */

	/*         on entry and return */
	/*            momcom - long */
	/*                     indicating that the chebyshev moments */
	/*                     have been computed for intervals of lengths */
	/*                     (fabs(b-a))*2**(-l), l=0,1,2, ..., momcom-1, */
	/*                     momcom.lt.maxp1 */

	/*            chebmo - sys_float */
	/*                     array of dimension (maxp1,25) containing the */
	/*                     chebyshev moments */

	/* ***references  (none) */
	/* ***routines called  qc25f,qelg,qpsrt,r1mach */
	/* ***end prologue  qawoe */




	/*            the dimension of rlist2 is determined by  the value of */
	/*            limexp in subroutine qelg (rlist2 should be of */
	/*            dimension (limexp+2) at least). */

	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           rlist2    - array of dimension at least limexp+2 */
	/*                       containing the part of the epsilon table */
	/*                       which is still needed for further computations */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest */
	/*                       error estimate */
	/*           errmax    - elist(maxerr) */
	/*           erlast    - error on the interval currently subdivided */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */
	/*           nres      - number of calls to the extrapolation routine */
	/*           numrl2    - number of elements in rlist2. if an appropriate */
	/*                       approximation to the compounded integral has */
	/*                       been obtained it is put in rlist2(numrl2) after */
	/*                       numrl2 has been increased by one */
	/*           small     - length of the smallest interval considered */
	/*                       up to now, multiplied by 1.5 */
	/*           erlarg    - sum of the errors over the intervals larger */
	/*                       than the smallest interval considered up to now */
	/*           extrap    - bool variable denoting that the routine is */
	/*                       attempting to perform extrapolation, i.e. before */
	/*                       subdividing the smallest interval we try to */
	/*                       decrease the value of erlarg */
	/*           noext     - bool variable denoting that extrapolation */
	/*                       is no longer allowed (true value) */

	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */
	/*           oflow is the largest positive magnitude. */

	/* ***first executable statement  qawoe */
	/* Parameter adjustments */
	--nnlog;
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;
	chebmo_dim1 = *maxp1;
	chebmo_offset = 1 + chebmo_dim1;
	chebmo -= chebmo_offset;

	/* Function Body */
	epmach = r1mach(c__4);

	/*         test on validity of parameters */
	/*         ------------------------------ */

	*ier = 0;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	alist__[1] = *a;
	blist[1] = *b;
	rlist[1] = 0.f;
	elist[1] = 0.f;
	iord[1] = 0;
	nnlog[1] = 0;
	/* Computing MAX */
	r__1 = epmach * 50.f;
	if ( ( *integr != 1 && *integr != 2) || 
		  (*epsabs <= 0.f && *epsrel < max(r__1,5e-15f)) || 
		  *icall < 1 || *maxp1 < 1) {
		*ier = 6;
	}
	if (*ier == 6) {
		goto L999;
	}

	/*           first approximation to the integral */
	/*           ----------------------------------- */

	domega = fabs(*omega);
	nrmom = 0;
	if (*icall > 1) {
		goto L5;
	}
	*momcom = 0;
L5:
	qc25f_(f, a, b, &domega, integr, &nrmom, maxp1, &c__0, result, 
			 abserr, neval, &defabs, &resabs, momcom, &chebmo[chebmo_offset]);

	/*           test on accuracy. */

	dres = fabs(*result);
	/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dres;
	errbnd = max(r__1,r__2);
	rlist[1] = *result;
	elist[1] = *abserr;
	iord[1] = 1;
	if (*abserr <= epmach * 100.f * defabs && *abserr > errbnd) {
		*ier = 2;
	}
	if (*limit == 1) {
		*ier = 1;
	}
	if (*ier != 0 || *abserr <= errbnd) {
		goto L200;
	}

	/*           initializations */
	/*           --------------- */

	uflow = r1mach(c__1);
	oflow = r1mach(c__2);
	errmax = *abserr;
	maxerr = 1;
	area = *result;
	errsum = *abserr;
	*abserr = oflow;
	nrmax = 1;
	extrap = false;
	noext = false;
	ierro = 0;
	iroff1 = 0;
	iroff2 = 0;
	iroff3 = 0;
	ktmin = 0;
	small = (r__1 = *b - *a, fabs(r__1)) * .75f;
	nres = 0;
	numrl2 = 0;
	extall = false;
	if ((r__1 = *b - *a, fabs(r__1)) * .5f * domega > 2.f) {
		goto L10;
	}
	numrl2 = 1;
	extall = true;
	rlist2[0] = *result;
L10:
	if ((r__1 = *b - *a, fabs(r__1)) * .25f * domega <= 2.f) {
		extall = true;
	}
	ksgn = -1;
	if (dres >= (1.f - epmach * 50.f) * defabs) {
		ksgn = 1;
	}

	/*           main do-loop */
	/*           ------------ */

	i__1 = *limit;
	for (*last = 2; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with the nrmax-th largest */
		/*           error estimate. */

		nrmom = nnlog[maxerr] + 1;
		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
		a2 = b1;
		b2 = blist[maxerr];
		erlast = errmax;
		qc25f_(f, &a1, &b1, &domega, integr, &nrmom, maxp1, &c__0, &
				 area1, &error1, &nev, &resabs, &defab1, momcom, &chebmo[
					 chebmo_offset]);
		*neval += nev;
		qc25f_(f, &a2, &b2, &domega, integr, &nrmom, maxp1, &c__1, &
				 area2, &error2, &nev, &resabs, &defab2, momcom, &chebmo[
					 chebmo_offset]);
		*neval += nev;

		/*           improve previous approximations to integral */
		/*           and error and test for accuracy. */

		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (defab1 == error1 || defab2 == error2) {
			goto L25;
		}
		if ((r__1 = rlist[maxerr] - area12, fabs(r__1)) > fabs(area12) * 
			 1e-5f || erro12 < errmax * .99f) {
			goto L20;
		}
		if (extrap) {
			++iroff2;
		}
		if (! extrap) {
			++iroff1;
		}
	L20:
		if (*last > 10 && erro12 > errmax) {
			++iroff3;
		}
	L25:
		rlist[maxerr] = area1;
		rlist[*last] = area2;
		nnlog[maxerr] = nrmom;
		nnlog[*last] = nrmom;
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(area);
		errbnd = max(r__1,r__2);

		/*           test for roundoff error and eventually */
		/*           set error flag */

		if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
			*ier = 2;
		}
		if (iroff2 >= 5) {
			ierro = 3;
		}

		/*           set error flag in the case that the number of */
		/*           subintervals equals limit. */

		if (*last == *limit) {
			*ier = 1;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at a point of the integration range. */

		/* Computing MAX */
		r__1 = fabs(a1), r__2 = fabs(b2);
		if (max(r__1,r__2) <= (epmach * 100.f + 1.f) * (fabs(a2) + uflow * 
																		 1e3f)) {
			*ier = 4;
		}

		/*           append the newly-created intervals to the list. */

		if (error2 > error1) {
			goto L30;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L40;
	L30:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine qpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the */
		/*           subinterval with nrmax-th largest error estimate (to be */
		/*           bisected next). */

	L40:
		qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (errsum <= errbnd) {
			goto L170;
		}
		if (*ier != 0) {
			goto L150;
		}
		if (*last == 2 && extall) {
			goto L120;
		}
		if (noext) {
			goto L140;
		}
		if (! extall) {
			goto L50;
		}
		erlarg -= erlast;
		if ((r__1 = b1 - a1, fabs(r__1)) > small) {
			erlarg += erro12;
		}
		if (extrap) {
			goto L70;
		}

		/*           test whether the interval to be bisected next is the */
		/*           smallest interval. */

	L50:
		width = (r__1 = blist[maxerr] - alist__[maxerr], fabs(r__1));
		if (width > small) {
			goto L140;
		}
		if (extall) {
			goto L60;
		}

		/*           test whether we can start with the extrapolation */
		/*           procedure (we do this if we integrate over the */
		/*           next interval with use of a gauss-kronrod rule - see */
		/*           subroutine qc25f). */

		small *= .5f;
		if (width * .25f * domega > 2.f) {
			goto L140;
		}
		extall = true;
		goto L130;
	L60:
		extrap = true;
		nrmax = 2;
	L70:
		if (ierro == 3 || erlarg <= ertest) {
			goto L90;
		}

		/*           the smallest interval has the largest error. */
		/*           before bisecting decrease the sum of the errors */
		/*           over the larger intervals (erlarg) and perform */
		/*           extrapolation. */

		jupbnd = *last;
		if (*last > *limit / 2 + 2) {
			jupbnd = *limit + 3 - *last;
		}
		id = nrmax;
		i__2 = jupbnd;
		for (k = id; k <= i__2; ++k) {
			maxerr = iord[nrmax];
			errmax = elist[maxerr];
			if ((r__1 = blist[maxerr] - alist__[maxerr], fabs(r__1)) > small) 
			{
				goto L140;
			}
			++nrmax;
			/* L80: */
		}

		/*           perform extrapolation. */

	L90:
		++numrl2;
		rlist2[numrl2 - 1] = area;
		if (numrl2 < 3) {
			goto L110;
		}
		qelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
		++ktmin;
		if (ktmin > 5 && *abserr < errsum * .001f) {
			*ier = 5;
		}
		if (abseps >= *abserr) {
			goto L100;
		}
		ktmin = 0;
		*abserr = abseps;
		*result = reseps;
		correc = erlarg;
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(reseps);
		ertest = max(r__1,r__2);
		/* ***jump out of do-loop */
		if (*abserr <= ertest) {
			goto L150;
		}

		/*           prepare bisection of the smallest interval. */

	L100:
		if (numrl2 == 1) {
			noext = true;
		}
		if (*ier == 5) {
			goto L150;
		}
	L110:
		maxerr = iord[1];
		errmax = elist[maxerr];
		nrmax = 1;
		extrap = false;
		small *= .5f;
		erlarg = errsum;
		goto L140;
	L120:
		small *= .5f;
		++numrl2;
		rlist2[numrl2 - 1] = area;
	L130:
		ertest = errbnd;
		erlarg = errsum;
	L140:
		;
	}

	/*           set the final result. */
	/*           --------------------- */

L150:
	if (*abserr == oflow || nres == 0) {
		goto L170;
	}
	if (*ier + ierro == 0) {
		goto L165;
	}
	if (ierro == 3) {
		*abserr += correc;
	}
	if (*ier == 0) {
		*ier = 3;
	}
	if (*result != 0.f && area != 0.f) {
		goto L160;
	}
	if (*abserr > errsum) {
		goto L170;
	}
	if (area == 0.f) {
		goto L190;
	}
	goto L165;
L160:
	if (*abserr / fabs(*result) > errsum / fabs(area)) {
		goto L170;
	}

	/*           test on divergence. */

L165:
	/* Computing MAX */
	r__1 = fabs(*result), r__2 = fabs(area);
	if (ksgn == -1 && max(r__1,r__2) <= defabs * .01f) {
		goto L190;
	}
	if (.01f > *result / area || *result / area > 100.f || errsum >= fabs(
			 area)) {
		*ier = 6;
	}
	goto L190;

	/*           compute global integral sum. */

L170:
	*result = 0.f;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L180: */
	}
	*abserr = errsum;
L190:
	if (*ier > 2) {
		--(*ier);
	}
L200:
	if (*integr == 2 && *omega < 0.f) {
		*result = -(*result);
	}
L999:
	return;
} /* qawoe_ */

void qawo_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  const sys_float *omega, const long *integr, 
			  const sys_float *epsabs, const sys_float *epsrel, 
			  sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, const long *leniw, const long *maxp1, 
			  const long *lenw, long *last, long *iwork, sys_float *work)
{
	long l1, l2, l3, l4, lvl;
	long limit, momcom;

	/* ***begin prologue  qawo */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             integrand with oscillatory cos or sin factor, */
	/*             clenshaw-curtis method, (end point) singularities, */
	/*             extrapolation, globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral */
	/*            i = integral of f(x)*w(x) over (a,b) */
	/*            where w(x) = cos(omega*x) or w(x) = sin(omega*x), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        computation of oscillatory integrals */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the function */
	/*                     f(x).  the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration */

	/*            omega  - sys_float */
	/*                     parameter in the integrand weight function */

	/*            integr - long */
	/*                     indicates which of the weight functions is used */
	/*                     integr = 1      w(x) = cos(omega*x) */
	/*                     integr = 2      w(x) = sin(omega*x) */
	/*                     if integr.ne.1.and.integr.ne.2, the routine will */
	/*                     end with ier = 6. */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if epsabs.le.0 and */
	/*                     epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of  integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                   - ier.gt.0 abnormal termination of the routine. */
	/*                             the estimates for integral and error are */
	/*                             less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             (= leniw/2) has been achieved. one can */
	/*                             allow more subdivisions by increasing the */
	/*                             value of leniw (and taking the according */
	/*                             dimension adjustments into account). */
	/*                             however, if this yields no improvement it */
	/*                             is advised to analyze the integrand in */
	/*                             order to determine the integration */
	/*                             difficulties. if the position of a local */
	/*                             difficulty can be determined (e.g. */
	/*                             singularity, discontinuity within the */
	/*                             interval) one will probably gain from */
	/*                             splitting up the interval at this polong */
	/*                             and calling the integrator on the */
	/*                             subranges. if possible, an appropriate */
	/*                             special-purpose integrator should be used */
	/*                             which is designed for handling the type of */
	/*                             difficulty involved. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                             the error may be under-estimated. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some interior points of the */
	/*                             integration interval. */
	/*                         = 4 the algorithm does not converge. */
	/*                             roundoff error is detected in the */
	/*                             extrapolation table. it is presumed that */
	/*                             the requested tolerance cannot be achieved */
	/*                             due to roundoff in the extrapolation */
	/*                             table, and that the returned result is */
	/*                             the best which can be obtained. */
	/*                         = 5 the integral is probably divergent, or */
	/*                             slowly convergent. it must be noted that */
	/*                             divergence can occur with any other value */
	/*                             of ier. */
	/*                         = 6 the input is invalid, because */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or (integr.ne.1 and integr.ne.2), */
	/*                             or leniw.lt.2 or maxp1.lt.1 or */
	/*                             lenw.lt.leniw*2+maxp1*25. */
	/*                             result, abserr, neval, last are set to */
	/*                             zero. except when leniw, maxp1 or lenw are */
	/*                             invalid, work(limit*2+1), work(limit*3+1), */
	/*                             iwork(1), iwork(limit+1) are set to zero, */
	/*                             work(1) is set to a and work(limit+1) to */
	/*                             b. */

	/*         dimensioning parameters */
	/*            leniw  - long */
	/*                     dimensioning parameter for iwork. */
	/*                     leniw/2 equals the maximum number of subintervals */
	/*                     allowed in the partition of the given integration */
	/*                     interval (a,b), leniw.ge.2. */
	/*                     if leniw.lt.2, the routine will end with ier = 6. */

	/*            maxp1  - long */
	/*                     gives an upper bound on the number of chebyshev */
	/*                     moments which can be stored, i.e. for the */
	/*                     intervals of lengths fabs(b-a)*2**(-l), */
	/*                     l=0,1, ..., maxp1-2, maxp1.ge.1 */
	/*                     if maxp1.lt.1, the routine will end with ier = 6. */

	/*            lenw   - long */
	/*                     dimensioning parameter for work */
	/*                     lenw must be at least leniw*2+maxp1*25. */
	/*                     if lenw.lt.(leniw*2+maxp1*25), the routine will */
	/*                     end with ier = 6. */

	/*            last   - long */
	/*                     on return, last equals the number of subintervals */
	/*                     produced in the subdivision process, which */
	/*                     determines the number of significant elements */
	/*                     actually in the work arrays. */

	/*         work arrays */
	/*            iwork  - long */
	/*                     vector of dimension at least leniw */
	/*                     on return, the first k elements of which contain */
	/*                     pointers to the error estimates over the */
	/*                     subintervals, such that work(limit*3+iwork(1)), .. */
	/*                     work(limit*3+iwork(k)) form a decreasing */
	/*                     sequence, with limit = lenw/2 , and k = last */
	/*                     if last.le.(limit/2+2), and k = limit+1-last */
	/*                     otherwise. */
	/*                     furthermore, iwork(limit+1), ..., iwork(limit+ */
	/*                     last) indicate the subdivision levels of the */
	/*                     subintervals, such that iwork(limit+i) = l means */
	/*                     that the subinterval numbered i is of length */
	/*                     fabs(b-a)*2**(1-l). */

	/*            work   - sys_float */
	/*                     vector of dimension at least lenw */
	/*                     on return */
	/*                     work(1), ..., work(last) contain the left */
	/*                      end points of the subintervals in the */
	/*                      partition of (a,b), */
	/*                     work(limit+1), ..., work(limit+last) contain */
	/*                      the right end points, */
	/*                     work(limit*2+1), ..., work(limit*2+last) contain */
	/*                      the integral approximations over the */
	/*                      subintervals, */
	/*                     work(limit*3+1), ..., work(limit*3+last) */
	/*                      contain the error estimates. */
	/*                     work(limit*4+1), ..., work(limit*4+maxp1*25) */
	/*                      provide space for storing the chebyshev moments. */
	/*                     note that limit = lenw/2. */

	/* ***references  (none) */
	/* ***routines called  qawoe,xerror */
	/* ***end prologue  qawo */




	/*         check validity of leniw, maxp1 and lenw. */

	/* ***first executable statement  qawo */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	if (*leniw < 2 || *maxp1 < 1 || *lenw < (*leniw << 1) + *maxp1 * 25) {
		goto L10;
	}

	/*         prepare call for qawoe */

	limit = *leniw / 2;
	l1 = limit + 1;
	l2 = limit + l1;
	l3 = limit + l2;
	l4 = limit + l3;
	qawoe_(f, a, b, omega, integr, epsabs, epsrel, &limit, &c__1, maxp1,
			 result, abserr, neval, ier, last, &work[1], &work[l1], &work[l2],
			 &work[l3], &iwork[1], &iwork[l1], &momcom, &work[l4]);

	/*         call error handler if necessary */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from  qawo", &c__26, ier, &lvl, 26);
	}
	return;
} /* qawo_ */

void qawse_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *alfa, const sys_float *beta, 
				const long *integr, const sys_float *epsabs, 
				const sys_float *epsrel, const long *limit, sys_float *result, 
				sys_float *abserr, long *neval, long *ier, sys_float *alist__, 
				sys_float *blist, sys_float *rlist, sys_float *elist, long *iord, long *last)
{
	/* System generated locals */
	int i__1;
	sys_float r__1, r__2;

	/* Local variables */
	long k;
	sys_float a1, a2, b1, b2, rg[25], rh[25], ri[25], rj[25];
	long nev;
	sys_float area;
	sys_float area1, area2, area12, erro12;
	long nrmax;
	sys_float uflow;
	long iroff1, iroff2;
	sys_float resas1, resas2, error1, error2, epmach, errbnd, centre, 
		errmax;
	long maxerr;
	sys_float errsum;

	/* ***begin prologue  qawse */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             algebraico-logarithmic end point singularities, */
	/*             clenshaw-curtis method */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f*w over (a,b), */
	/*            (where w shows a singular behaviour at the end points, */
	/*            see parameter integr). */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        integration of functions having algebraico-logarithmic */
	/*        end point singularities */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration, b.gt.a */
	/*                     if b.le.a, the routine will end with ier = 6. */

	/*            alfa   - sys_float */
	/*                     parameter in the weight function, alfa.gt.(-1) */
	/*                     if alfa.le.(-1), the routine will end with */
	/*                     ier = 6. */

	/*            beta   - sys_float */
	/*                     parameter in the weight function, beta.gt.(-1) */
	/*                     if beta.le.(-1), the routine will end with */
	/*                     ier = 6. */

	/*            integr - long */
	/*                     indicates which weight function is to be used */
	/*                     = 1  (x-a)**alfa*(b-x)**beta */
	/*                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a) */
	/*                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x) */
	/*                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x) */
	/*                     if integr.lt.1 or integr.gt.4, the routine */
	/*                     will end with ier = 6. */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*            limit  - long */
	/*                     gives an upper bound on the number of subintervals */
	/*                     in the partition of (a,b), limit.ge.2 */
	/*                     if limit.lt.2, the routine will end with ier = 6. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for the integral and error */
	/*                             are less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                         = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit. however, if this yields no */
	/*                             improvement, it is advised to analyze the */
	/*                             integrand in order to determine the */
	/*                             integration difficulties which prevent the */
	/*                             requested tolerance from being achieved. */
	/*                             in case of a jump discontinuity or a local */
	/*                             singularity of algebraico-logarithmic type */
	/*                             at one or more interior points of the */
	/*                             integration range, one should proceed by */
	/*                             splitting up the interval at these */
	/*                             points and calling the integrator on the */
	/*                             subranges. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 6 the input is invalid, because */
	/*                             b.le.a or alfa.le.(-1) or beta.le.(-1), or */
	/*                             integr.lt.1 or integr.gt.4, or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                             or limit.lt.2. */
	/*                             result, abserr, neval, rlist(1), elist(1), */
	/*                             iord(1) and last are set to zero. alist(1) */
	/*                             and blist(1) are set to a and b */
	/*                             respectively. */

	/*            alist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the left */
	/*                     end points of the subintervals in the partition */
	/*                     of the given integration range (a,b) */

	/*            blist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the right */
	/*                     end points of the subintervals in the partition */
	/*                     of the given integration range (a,b) */

	/*            rlist  - sys_float */
	/*                     vector of dimension at least limit,the first */
	/*                      last  elements of which are the integral */
	/*                     approximations on the subintervals */

	/*            elist  - sys_float */
	/*                     vector of dimension at least limit, the first */
	/*                      last  elements of which are the moduli of the */
	/*                     absolute error estimates on the subintervals */

	/*            iord   - long */
	/*                     vector of dimension at least limit, the first k */
	/*                     of which are pointers to the error */
	/*                     estimates over the subintervals, so that */
	/*                     elist(iord(1)), ..., elist(iord(k)) with k = last */
	/*                     if last.le.(limit/2+2), and k = limit+1-last */
	/*                     otherwise form a decreasing sequence */

	/*            last   - long */
	/*                     number of subintervals actually produced in */
	/*                     the subdivision process */

	/* ***references  (none) */
	/* ***routines called  qc25s,qmomo,qpsrt,r1mach */
	/* ***end prologue  qawse */




	/*            list of major variables */
	/*            ----------------------- */

	/*           alist     - list of left end points of all subintervals */
	/*                       considered up to now */
	/*           blist     - list of right end points of all subintervals */
	/*                       considered up to now */
	/*           rlist(i)  - approximation to the integral over */
	/*                       (alist(i),blist(i)) */
	/*           elist(i)  - error estimate applying to rlist(i) */
	/*           maxerr    - pointer to the interval with largest */
	/*                       error estimate */
	/*           errmax    - elist(maxerr) */
	/*           area      - sum of the integrals over the subintervals */
	/*           errsum    - sum of the errors over the subintervals */
	/*           errbnd    - requested accuracy max(epsabs,epsrel* */
	/*                       fabs(result)) */
	/*           *****1    - variable for the left subinterval */
	/*           *****2    - variable for the right subinterval */
	/*           last      - index for subdivision */


	/*            machine dependent constants */
	/*            --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qawse */
	/* Parameter adjustments */
	--iord;
	--elist;
	--rlist;
	--blist;
	--alist__;

	/* Function Body */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);

	/*           test on validity of parameters */
	/*           ------------------------------ */

	*ier = 6;
	*neval = 0;
	*last = 0;
	rlist[1] = 0.f;
	elist[1] = 0.f;
	iord[1] = 0;
	*result = 0.f;
	*abserr = 0.f;
	/* Computing MAX */
	r__1 = epmach * 50.f;
	if (*b <= *a || ( *epsabs == 0.f && *epsrel < max(r__1,5e-15f)) || *alfa <= 
	    -1.f || *beta <= -1.f || *integr < 1 || *integr > 4 || *limit < 2)
	{
		goto L999;
	}
	*ier = 0;

	/*           compute the modified chebyshev moments. */

	qmomo_(alfa, beta, ri, rj, rg, rh, integr);

	/*           integrate over the intervals (a,(a+b)/2) */
	/*           and ((a+b)/2,b). */

	centre = (*b + *a) * .5f;
	qc25s_(f, a, b, a, &centre, alfa, beta, ri, rj, rg, rh, &area1, &
			 error1, &resas1, integr, &nev);
	*neval = nev;
	qc25s_(f, a, b, &centre, b, alfa, beta, ri, rj, rg, rh, &area2, &
			 error2, &resas2, integr, &nev);
	*last = 2;
	*neval += nev;
	*result = area1 + area2;
	*abserr = error1 + error2;

	/*           test on accuracy. */

	/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * fabs(*result);
	errbnd = max(r__1,r__2);

	/*           initialization */
	/*           -------------- */

	if (error2 > error1) {
		goto L10;
	}
	alist__[1] = *a;
	alist__[2] = centre;
	blist[1] = centre;
	blist[2] = *b;
	rlist[1] = area1;
	rlist[2] = area2;
	elist[1] = error1;
	elist[2] = error2;
	goto L20;
L10:
	alist__[1] = centre;
	alist__[2] = *a;
	blist[1] = *b;
	blist[2] = centre;
	rlist[1] = area2;
	rlist[2] = area1;
	elist[1] = error2;
	elist[2] = error1;
L20:
	iord[1] = 1;
	iord[2] = 2;
	if (*limit == 2) {
		*ier = 1;
	}
	if (*abserr <= errbnd || *ier == 1) {
		goto L999;
	}
	errmax = elist[1];
	maxerr = 1;
	nrmax = 1;
	area = *result;
	errsum = *abserr;
	iroff1 = 0;
	iroff2 = 0;

	/*            main do-loop */
	/*            ------------ */

	i__1 = *limit;
	for (*last = 3; *last <= i__1; ++(*last)) {

		/*           bisect the subinterval with largest error estimate. */

		a1 = alist__[maxerr];
		b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
		a2 = b1;
		b2 = blist[maxerr];

		qc25s_(f, a, b, &a1, &b1, alfa, beta, ri, rj, rg, rh, &area1, &
				 error1, &resas1, integr, &nev);
		*neval += nev;
		qc25s_(f, a, b, &a2, &b2, alfa, beta, ri, rj, rg, rh, &area2, &
				 error2, &resas2, integr, &nev);
		*neval += nev;

		/*           improve previous approximations integral and error */
		/*           and test for accuracy. */

		area12 = area1 + area2;
		erro12 = error1 + error2;
		errsum = errsum + erro12 - errmax;
		area = area + area12 - rlist[maxerr];
		if (*a == a1 || *b == b2) {
			goto L30;
		}
		if (resas1 == error1 || resas2 == error2) {
			goto L30;
		}

		/*           test for roundoff error. */

		if ((r__1 = rlist[maxerr] - area12, fabs(r__1)) < fabs(area12) * 
			 1e-5f && erro12 >= errmax * .99f) {
			++iroff1;
		}
		if (*last > 10 && erro12 > errmax) {
			++iroff2;
		}
	L30:
		rlist[maxerr] = area1;
		rlist[*last] = area2;

		/*           test on accuracy. */

		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(area);
		errbnd = max(r__1,r__2);
		if (errsum <= errbnd) {
			goto L35;
		}

		/*           set error flag in the case that the number of interval */
		/*           bisections exceeds limit. */

		if (*last == *limit) {
			*ier = 1;
		}


		/*           set error flag in the case of roundoff error. */

		if (iroff1 >= 6 || iroff2 >= 20) {
			*ier = 2;
		}

		/*           set error flag in the case of bad integrand behaviour */
		/*           at interior points of integration range. */

		/* Computing MAX */
		r__1 = fabs(a1), r__2 = fabs(b2);
		if (max(r__1,r__2) <= (epmach * 100.f + 1.f) * (fabs(a2) + uflow * 
																		 1e3f)) {
			*ier = 3;
		}

		/*           append the newly-created intervals to the list. */

	L35:
		if (error2 > error1) {
			goto L40;
		}
		alist__[*last] = a2;
		blist[maxerr] = b1;
		blist[*last] = b2;
		elist[maxerr] = error1;
		elist[*last] = error2;
		goto L50;
	L40:
		alist__[maxerr] = a2;
		alist__[*last] = a1;
		blist[*last] = b1;
		rlist[maxerr] = area2;
		rlist[*last] = area1;
		elist[maxerr] = error2;
		elist[*last] = error1;

		/*           call subroutine qpsrt to maintain the descending ordering */
		/*           in the list of error estimates and select the */
		/*           subinterval with largest error estimate (to be */
		/*           bisected next). */

	L50:
		qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
		/* ***jump out of do-loop */
		if (*ier != 0 || errsum <= errbnd) {
			goto L70;
		}
		/* L60: */
	}

	/*           compute final result. */
	/*           --------------------- */

L70:
	*result = 0.f;
	i__1 = *last;
	for (k = 1; k <= i__1; ++k) {
		*result += rlist[k];
		/* L80: */
	}
	*abserr = errsum;
L999:
	return;
} /* qawse_ */

void qaws_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  const sys_float *alfa, const sys_float *beta, 
			  const long *integr, const sys_float *epsabs, 
			  const sys_float *epsrel, sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, const long *limit, const long *lenw, 
			  long *last, long *iwork, sys_float *work)
{
	long l1, l2, l3, lvl;

	/* ***begin prologue  qaws */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1 */
	/* ***keywords  automatic integrator, special-purpose, */
	/*             algebraico-logarithmic end-point singularities, */
	/*             clenshaw-curtis, globally adaptive */
	/* ***author  piessens,robert,appl. math. & progr. div. -k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the routine calculates an approximation result to a given */
	/*            definite integral i = integral of f*w over (a,b), */
	/*            (where w shows a singular behaviour at the end points */
	/*            see parameter integr). */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/*        integration of functions having algebraico-logarithmic */
	/*        end point singularities */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*            f      - sys_float */
	/*                     function subprogram defining the integrand */
	/*                     function f(x). the actual name for f needs to be */
	/*                     declared e x t e r n a l in the driver program. */

	/*            a      - sys_float */
	/*                     lower limit of integration */

	/*            b      - sys_float */
	/*                     upper limit of integration, b.gt.a */
	/*                     if b.le.a, the routine will end with ier = 6. */

	/*            alfa   - sys_float */
	/*                     parameter in the integrand function, alfa.gt.(-1) */
	/*                     if alfa.le.(-1), the routine will end with */
	/*                     ier = 6. */

	/*            beta   - sys_float */
	/*                     parameter in the integrand function, beta.gt.(-1) */
	/*                     if beta.le.(-1), the routine will end with */
	/*                     ier = 6. */

	/*            integr - long */
	/*                     indicates which weight function is to be used */
	/*                     = 1  (x-a)**alfa*(b-x)**beta */
	/*                     = 2  (x-a)**alfa*(b-x)**beta*log(x-a) */
	/*                     = 3  (x-a)**alfa*(b-x)**beta*log(b-x) */
	/*                     = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x) */
	/*                     if integr.lt.1 or integr.gt.4, the routine */
	/*                     will end with ier = 6. */

	/*            epsabs - sys_float */
	/*                     absolute accuracy requested */
	/*            epsrel - sys_float */
	/*                     relative accuracy requested */
	/*                     if  epsabs.le.0 */
	/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                     the routine will end with ier = 6. */

	/*         on return */
	/*            result - sys_float */
	/*                     approximation to the integral */

	/*            abserr - sys_float */
	/*                     estimate of the modulus of the absolute error, */
	/*                     which should equal or exceed fabs(i-result) */

	/*            neval  - long */
	/*                     number of integrand evaluations */

	/*            ier    - long */
	/*                     ier = 0 normal and reliable termination of the */
	/*                             routine. it is assumed that the requested */
	/*                             accuracy has been achieved. */
	/*                     ier.gt.0 abnormal termination of the routine */
	/*                             the estimates for the integral and error */
	/*                             are less reliable. it is assumed that the */
	/*                             requested accuracy has not been achieved. */
	/*            error messages */
	/*                     ier = 1 maximum number of subdivisions allowed */
	/*                             has been achieved. one can allow more */
	/*                             subdivisions by increasing the value of */
	/*                             limit (and taking the according dimension */
	/*                             adjustments into account). however, if */
	/*                             this yields no improvement it is advised */
	/*                             to analyze the integrand, in order to */
	/*                             determine the integration difficulties */
	/*                             which prevent the requested tolerance from */
	/*                             being achieved. in case of a jump */
	/*                             discontinuity or a local singularity */
	/*                             of algebraico-logarithmic type at one or */
	/*                             more interior points of the integration */
	/*                             range, one should proceed by splitting up */
	/*                             the interval at these points and calling */
	/*                             the integrator on the subranges. */
	/*                         = 2 the occurrence of roundoff error is */
	/*                             detected, which prevents the requested */
	/*                             tolerance from being achieved. */
	/*                         = 3 extremely bad integrand behaviour occurs */
	/*                             at some points of the integration */
	/*                             interval. */
	/*                         = 6 the input is invalid, because */
	/*                             b.le.a or alfa.le.(-1) or beta.le.(-1) or */
	/*                             or integr.lt.1 or integr.gt.4 or */
	/*                             (epsabs.le.0 and */
	/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
	/*                             or limit.lt.2 or lenw.lt.limit*4. */
	/*                             result, abserr, neval, last are set to */
	/*                             zero. except when lenw or limit is invalid */
	/*                             iwork(1), work(limit*2+1) and */
	/*                             work(limit*3+1) are set to zero, work(1) */
	/*                             is set to a and work(limit+1) to b. */

	/*         dimensioning parameters */
	/*            limit  - long */
	/*                     dimensioning parameter for iwork */
	/*                     limit determines the maximum number of */
	/*                     subintervals in the partition of the given */
	/*                     integration interval (a,b), limit.ge.2. */
	/*                     if limit.lt.2, the routine will end with ier = 6. */

	/*            lenw   - long */
	/*                     dimensioning parameter for work */
	/*                     lenw must be at least limit*4. */
	/*                     if lenw.lt.limit*4, the routine will end */
	/*                     with ier = 6. */

	/*            last   - long */
	/*                     on return, last equals the number of */
	/*                     subintervals produced in the subdivision process, */
	/*                     which determines the significant number of */
	/*                     elements actually in the work arrays. */

	/*         work arrays */
	/*            iwork  - long */
	/*                     vector of dimension limit, the first k */
	/*                     elements of which contain pointers */
	/*                     to the error estimates over the subintervals, */
	/*                     such that work(limit*3+iwork(1)), ..., */
	/*                     work(limit*3+iwork(k)) form a decreasing */
	/*                     sequence with k = last if last.le.(limit/2+2), */
	/*                     and k = limit+1-last otherwise */

	/*            work   - sys_float */
	/*                     vector of dimension lenw */
	/*                     on return */
	/*                     work(1), ..., work(last) contain the left */
	/*                      end points of the subintervals in the */
	/*                      partition of (a,b), */
	/*                     work(limit+1), ..., work(limit+last) contain */
	/*                      the right end points, */
	/*                     work(limit*2+1), ..., work(limit*2+last) */
	/*                      contain the integral approximations over */
	/*                      the subintervals, */
	/*                     work(limit*3+1), ..., work(limit*3+last) */
	/*                      contain the error estimates. */

	/* ***references  (none) */
	/* ***routines called  qawse,xerror */
	/* ***end prologue  qaws */




	/*         check validity of limit and lenw. */

	/* ***first executable statement  qaws */
	/* Parameter adjustments */
	--iwork;
	--work;

	/* Function Body */
	*ier = 6;
	*neval = 0;
	*last = 0;
	*result = 0.f;
	*abserr = 0.f;
	if (*limit < 2 || *lenw < *limit << 2) {
		goto L10;
	}

	/*         prepare call for qawse. */

	l1 = *limit + 1;
	l2 = *limit + l1;
	l3 = *limit + l2;

	qawse_(f, a, b, alfa, beta, integr, epsabs, epsrel, limit, result, 
			 abserr, neval, ier, &work[1], &work[l1], &work[l2], &work[l3], &
			 iwork[1], last);

	/*         call error handler if necessary. */

	lvl = 0;
L10:
	if (*ier == 6) {
		lvl = 1;
	}
	if (*ier != 0) {
		xerror_("abnormal return from  qaws", &c__26, ier, &lvl, 26);
	}
	return;
} /* qaws_ */

void qc25c_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *c__, sys_float *result,
				sys_float *abserr, long *krul, long *neval)
{
	/* Initialized data */

	sys_float x[11] = { .9914448613738104f,.9659258262890683f,
								  .9238795325112868f,.8660254037844386f,.7933533402912352f,
								  .7071067811865475f,.6087614290087206f,.5f,.3826834323650898f,
								  .2588190451025208f,.1305261922200516f };

	/* System generated locals */
	sys_float r__1;

	/* Local variables */
	long i__, k;
	sys_float u, p2, p3, p4, cc;
	long kp;
	sys_float ak22, fval[25], res12, res24;
	long isym;
	sys_float amom0, amom1, amom2, cheb12[13], cheb24[25];
	sys_float hlgth, centr;
	sys_float resabs, resasc;

	/* ***begin prologue  qc25c */
	/* ***date written   810101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a2,j4 */
	/* ***keywords  25-point clenshaw-curtis integration */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f*w over (a,b) with */
	/*            error estimate, where w(x) = 1/(x-c) */
	/* ***description */

	/*        integration rules for the computation of cauchy */
	/*        principal value integrals */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*           f      - sys_float */
	/*                    function subprogram defining the integrand function */
	/*                    f(x). the actual name for f needs to be declared */
	/*                    e x t e r n a l  in the driver program. */

	/*           a      - sys_float */
	/*                    left end point of the integration interval */

	/*           b      - sys_float */
	/*                    right end point of the integration interval, b.gt.a */

	/*           c      - sys_float */
	/*                    parameter in the weight function */

	/*           result - sys_float */
	/*                    approximation to the integral */
	/*                    result is computed by using a generalized */
	/*                    clenshaw-curtis method if c lies within ten percent */
	/*                    of the integration interval. in the other case the */
	/*                    15-point kronrod rule obtained by optimal addition */
	/*                    of abscissae to the 7-point gauss rule, is applied. */

	/*           abserr - sys_float */
	/*                    estimate of the modulus of the absolute error, */
	/*                    which should equal or exceed fabs(i-result) */

	/*           krul   - long */
	/*                    key which is decreased by 1 if the 15-polong */
	/*                    gauss-kronrod scheme has been used */

	/*           neval  - long */
	/*                    number of integrand evaluations */

	/* ***references  (none) */
	/* ***routines called  qcheb,qk15w,qwgtc */
	/* ***end prologue  qc25c */




	/*           the vector x contains the values cos(k*pi/24), */
	/*           k = 1, ..., 11, to be used for the chebyshev series */
	/*           expansion of f */


	/*           list of major variables */
	/*           ---------------------- */
	/*           fval   - value of the function f at the points */
	/*                    cos(k*pi/24),  k = 0, ..., 24 */
	/*           cheb12 - chebyshev series expansion coefficients, */
	/*                    for the function f, of degree 12 */
	/*           cheb24 - chebyshev series expansion coefficients, */
	/*                    for the function f, of degree 24 */
	/*           res12  - approximation to the integral corresponding */
	/*                    to the use of cheb12 */
	/*           res24  - approximation to the integral corresponding */
	/*                    to the use of cheb24 */
	/*           qwgtc - external function subprogram defining */
	/*                    the weight function */
	/*           hlgth  - half-length of the interval */
	/*           centr  - mid point of the interval */


	/*           check the position of c. */

	set_NaN( p2 );
	set_NaN( p3 );
	set_NaN( p4 );
	kp = 0;

	/* ***first executable statement  qc25c */
	cc = (2.f * *c__ - *b - *a) / (*b - *a);
	if (fabs(cc) < 1.1f) {
		goto L10;
	}

	/*           apply the 15-point gauss-kronrod scheme. */

	--(*krul);
	qk15w_(f, qwgtc_, c__, &p2, &p3, &p4, &kp, a, b, result, 
			 abserr, &resabs, &resasc);
	*neval = 15;
	if (resasc == *abserr) {
		++(*krul);
	}
	goto L50;

	/*           use the generalized clenshaw-curtis method. */

L10:
	hlgth = (*b - *a) * .5f;
	centr = (*b + *a) * .5f;
	*neval = 25;
	r__1 = hlgth + centr;
	fval[0] = f(r__1) * .5f;
	fval[12] = f(centr);
	r__1 = centr - hlgth;
	fval[24] = f(r__1) * .5f;
	for (i__ = 2; i__ <= 12; ++i__) {
		u = hlgth * x[i__ - 2];
		isym = 26 - i__;
		r__1 = u + centr;
		fval[i__ - 1] = f(r__1);
		r__1 = centr - u;
		fval[isym - 1] = f(r__1);
		/* L20: */
	}

	/*           compute the chebyshev series expansion. */

	qcheb_(x, fval, cheb12, cheb24);

	/*           the modified chebyshev moments are computed */
	/*           by forward recursion, using amom0 and amom1 */
	/*           as starting values. */

	amom0 = log((r__1 = (1.f - cc) / (cc + 1.f), fabs(r__1)));
	amom1 = cc * amom0 + 2.f;
	res12 = cheb12[0] * amom0 + cheb12[1] * amom1;
	res24 = cheb24[0] * amom0 + cheb24[1] * amom1;
	for (k = 3; k <= 13; ++k) {
		amom2 = cc * 2.f * amom1 - amom0;
		ak22 = (sys_float) ((k - 2) * (k - 2));
		if (k / 2 << 1 == k) {
			amom2 -= 4.f / (ak22 - 1.f);
		}
		res12 += cheb12[k - 1] * amom2;
		res24 += cheb24[k - 1] * amom2;
		amom0 = amom1;
		amom1 = amom2;
		/* L30: */
	}
	for (k = 14; k <= 25; ++k) {
		amom2 = cc * 2.f * amom1 - amom0;
		ak22 = (sys_float) ((k - 2) * (k - 2));
		if (k / 2 << 1 == k) {
			amom2 -= 4.f / (ak22 - 1.f);
		}
		res24 += cheb24[k - 1] * amom2;
		amom0 = amom1;
		amom1 = amom2;
		/* L40: */
	}
	*result = res24;
	*abserr = (r__1 = res24 - res12, fabs(r__1));
L50:
	return;
} /* qc25c_ */

void qc25f_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *omega, const long *integr, 
				const long *nrmom, const long *maxp1, const long *ksave, 
				sys_float *result, 
				sys_float *abserr, long *neval, sys_float *resabs, sys_float *resasc, long *
				momcom, sys_float *chebmo)
{
	/* Initialized data */

	sys_float x[11] = { .9914448613738104f,.9659258262890683f,
								  .9238795325112868f,.8660254037844386f,.7933533402912352f,
								  .7071067811865475f,.6087614290087206f,.5f,.3826834323650898f,
								  .2588190451025208f,.1305261922200516f };

	/* System generated locals */
	int chebmo_dim1, chebmo_offset, i__1;
	sys_float r__1, r__2;

	/* Local variables */
	sys_float d__[25];
	long i__, j, k, m=0;
	sys_float v[28], d1[25], d2[25], p2, p3, p4, ac, an, as, an2, ass, par2,
		conc, asap, par22, fval[25], estc, cons;
	long iers;
	sys_float ests;
	long isym, noeq1;
	sys_float cheb12[13], cheb24[25];
	sys_float resc12, resc24, hlgth, centr, ress12, ress24, oflow;
	long noequ;
	sys_float cospar, sinpar, parint;

	/* ***begin prologue  qc25f */
	/* ***date written   810101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a2 */
	/* ***keywords  integration rules for functions with cos or sin */
	/*             factor, clenshaw-curtis, gauss-kronrod */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute the integral i=integral of f(x) over (a,b) */
	/*            where w(x) = cos(omega*x) or (wx)=sin(omega*x) */
	/*            and to compute j=integral of fabs(f) over (a,b). for small */
	/*            value of omega or small intervals (a,b) 15-point gauss- */
	/*            kronrod rule used. otherwise generalized clenshaw-curtis us */
	/* ***description */

	/*        integration rules for functions with cos or sin factor */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*         on entry */
	/*           f      - sys_float */
	/*                    function subprogram defining the integrand */
	/*                    function f(x). the actual name for f needs to */
	/*                    be declared e x t e r n a l in the calling program. */

	/*           a      - sys_float */
	/*                    lower limit of integration */

	/*           b      - sys_float */
	/*                    upper limit of integration */

	/*           omega  - sys_float */
	/*                    parameter in the weight function */

	/*           integr - long */
	/*                    indicates which weight function is to be used */
	/*                       integr = 1   w(x) = cos(omega*x) */
	/*                       integr = 2   w(x) = sin(omega*x) */

	/*           nrmom  - long */
	/*                    the length of interval (a,b) is equal to the length */
	/*                    of the original integration interval divided by */
	/*                    2**nrmom (we suppose that the routine is used in an */
	/*                    adaptive integration process, otherwise set */
	/*                    nrmom = 0). nrmom must be zero at the first call. */

	/*           maxp1  - long */
	/*                    gives an upper bound on the number of chebyshev */
	/*                    moments which can be stored, i.e. for the */
	/*                    intervals of lengths fabs(bb-aa)*2**(-l), */
	/*                    l = 0,1,2, ..., maxp1-2. */

	/*           ksave  - long */
	/*                    key which is one when the moments for the */
	/*                    current interval have been computed */

	/*         on return */
	/*           result - sys_float */
	/*                    approximation to the integral i */

	/*           abserr - sys_float */
	/*                    estimate of the modulus of the absolute */
	/*                    error, which should equal or exceed fabs(i-result) */

	/*           neval  - long */
	/*                    number of integrand evaluations */

	/*           resabs - sys_float */
	/*                    approximation to the integral j */

	/*           resasc - sys_float */
	/*                    approximation to the integral of fabs(f-i/(b-a)) */

	/*         on entry and return */
	/*           momcom - long */
	/*                    for each interval length we need to compute the */
	/*                    chebyshev moments. momcom counts the number of */
	/*                    intervals for which these moments have already been */
	/*                    computed. if nrmom.lt.momcom or ksave = 1, the */
	/*                    chebyshev moments for the interval (a,b) have */
	/*                    already been computed and stored, otherwise we */
	/*                    compute them and we increase momcom. */

	/*           chebmo - sys_float */
	/*                    array of dimension at least (maxp1,25) containing */
	/*                    the modified chebyshev moments for the first momcom */
	/*                    momcom interval lengths */

	/* ***references  (none) */
	/* ***routines called  qcheb,qk15w,qwgtf,r1mach,sgtsl */
	/* ***end prologue  qc25f */




	/*           the vector x contains the values cos(k*pi/24) */
	/*           k = 1, ...,11, to be used for the chebyshev expansion of f */

	/* Parameter adjustments */
	chebmo_dim1 = *maxp1;
	chebmo_offset = 1 + chebmo_dim1;
	chebmo -= chebmo_offset;

	set_NaN( p2 );
	set_NaN( p3 );
	set_NaN( p4 );

	/* Function Body */

	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the integration interval */
	/*           hlgth  - half-length of the integration interval */
	/*           fval   - value of the function f at the points */
	/*                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5, */
	/*                    k = 0, ..., 24 */
	/*           cheb12 - coefficients of the chebyshev series expansion */
	/*                    of degree 12, for the function f, in the */
	/*                    interval (a,b) */
	/*           cheb24 - coefficients of the chebyshev series expansion */
	/*                    of degree 24, for the function f, in the */
	/*                    interval (a,b) */
	/*           resc12 - approximation to the integral of */
	/*                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a)) */
	/*                    over (-1,+1), using the chebyshev series */
	/*                    expansion of degree 12 */
	/*           resc24 - approximation to the same integral, using the */
	/*                    chebyshev series expansion of degree 24 */
	/*           ress12 - the analogue of resc12 for the sine */
	/*           ress24 - the analogue of resc24 for the sine */


	/*           machine dependent constant */
	/*           -------------------------- */

	/*           oflow is the largest positive magnitude. */

	/* ***first executable statement  qc25f */
	oflow = r1mach(c__2);

	centr = (*b + *a) * .5f;
	hlgth = (*b - *a) * .5f;
	parint = *omega * hlgth;

	/*           compute the integral using the 15-point gauss-kronrod */
	/*           formula if the value of the parameter in the integrand */
	/*           is small. */

	if (fabs(parint) > 2.f) {
		goto L10;
	}
	qk15w_(f, qwgtf_, omega, &p2, &p3, &p4, integr, a, b, result, 
			 abserr, resabs, resasc);
	*neval = 15;
	goto L170;

	/*           compute the integral using the generalized clenshaw- */
	/*           curtis method. */

L10:
	conc = hlgth * cos(centr * *omega);
	cons = hlgth * sin(centr * *omega);
	*resasc = oflow;
	*neval = 25;

	/*           check whether the chebyshev moments for this interval */
	/*           have already been computed. */

	if (*nrmom < *momcom || *ksave == 1) {
		goto L120;
	}

	/*           compute a new set of chebyshev moments. */

	m = *momcom + 1;
	par2 = parint * parint;
	par22 = par2 + 2.f;
	sinpar = sin(parint);
	cospar = cos(parint);

	/*           compute the chebyshev moments with respect to cosine. */

	v[0] = sinpar * 2.f / parint;
	v[1] = (cospar * 8.f + (par2 + par2 - 8.f) * sinpar / parint) / par2;
	v[2] = ((par2 - 12.f) * 32.f * cospar + ((par2 - 80.f) * par2 + 192.f) * 
			  2.f * sinpar / parint) / (par2 * par2);
	ac = cospar * 8.f;
	as = parint * 24.f * sinpar;
	if (fabs(parint) > 24.f) {
		goto L30;
	}

	/*           compute the chebyshev moments as the */
	/*           solutions of a boundary value problem with 1 */
	/*           initial value (v(3)) and 1 end value (computed */
	/*           using an asymptotic formula). */

	noequ = 25;
	noeq1 = noequ - 1;
	an = 6.f;
	i__1 = noeq1;
	for (k = 1; k <= i__1; ++k) {
		an2 = an * an;
		d__[k - 1] = (an2 - 4.f) * -2.f * (par22 - an2 - an2);
		d2[k - 1] = (an - 1.f) * (an - 2.f) * par2;
		d1[k] = (an + 3.f) * (an + 4.f) * par2;
		v[k + 2] = as - (an2 - 4.f) * ac;
		an += 2.f;
		/* L20: */
	}
	an2 = an * an;
	d__[noequ - 1] = (an2 - 4.f) * -2.f * (par22 - an2 - an2);
	v[noequ + 2] = as - (an2 - 4.f) * ac;
	v[3] -= par2 * 56.f * v[2];
	ass = parint * sinpar;
	asap = (((((par2 * 210.f - 1.f) * cospar - (par2 * 105.f - 63.f) * ass) / 
				 an2 - (1.f - par2 * 15.f) * cospar + ass * 15.f) / an2 - cospar + 
				ass * 3.f) / an2 - cospar) / an2;
	v[noequ + 2] -= asap * 2.f * par2 * (an - 1.f) * (an - 2.f);

	/*           solve the tridiagonal system by means of gaussian */
	/*           elimination with partial pivoting. */

	sgtsl_(&noequ, d1, d__, d2, &v[3], &iers);
	goto L50;

	/*           compute the chebyshev moments by means of forward */
	/*           recursion. */

L30:
	an = 4.f;
	for (i__ = 4; i__ <= 13; ++i__) {
		an2 = an * an;
		v[i__ - 1] = ((an2 - 4.f) * ((par22 - an2 - an2) * 2.f * v[i__ - 2] - 
											  ac) + as - par2 * (an + 1.f) * (an + 2.f) * v[i__ - 3]) / (
												  par2 * (an - 1.f) * (an - 2.f));
		an += 2.f;
		/* L40: */
	}
L50:
	for (j = 1; j <= 13; ++j) {
		chebmo[m + ((j << 1) - 1) * chebmo_dim1] = v[j - 1];
		/* L60: */
	}

	/*           compute the chebyshev moments with respect to sine. */

	v[0] = (sinpar - parint * cospar) * 2.f / par2;
	v[1] = (18.f - 48.f / par2) * sinpar / par2 + (48.f / par2 - 2.f) * 
		cospar / parint;
	ac = parint * -24.f * cospar;
	as = sinpar * -8.f;
	if (fabs(parint) > 24.f) {
		goto L80;
	}

	/*           compute the chebyshev moments as the */
	/*           solutions of a boundary value problem with 1 */
	/*           initial value (v(2)) and 1 end value (computed */
	/*           using an asymptotic formula). */

	an = 5.f;
	i__1 = noeq1;
	for (k = 1; k <= i__1; ++k) {
		an2 = an * an;
		d__[k - 1] = (an2 - 4.f) * -2.f * (par22 - an2 - an2);
		d2[k - 1] = (an - 1.f) * (an - 2.f) * par2;
		d1[k] = (an + 3.f) * (an + 4.f) * par2;
		v[k + 1] = ac + (an2 - 4.f) * as;
		an += 2.f;
		/* L70: */
	}
	an2 = an * an;
	d__[noequ - 1] = (an2 - 4.f) * -2.f * (par22 - an2 - an2);
	v[noequ + 1] = ac + (an2 - 4.f) * as;
	v[2] -= par2 * 42.f * v[1];
	ass = parint * cospar;
	asap = (((((par2 * 105.f - 63.f) * ass + (par2 * 210.f - 1.f) * sinpar) / 
				 an2 + (par2 * 15.f - 1.f) * sinpar - ass * 15.f) / an2 - ass * 
				3.f - sinpar) / an2 - sinpar) / an2;
	v[noequ + 1] -= asap * 2.f * par2 * (an - 1.f) * (an - 2.f);

	/*           solve the tridiagonal system by means of gaussian */
	/*           elimination with partial pivoting. */

	sgtsl_(&noequ, d1, d__, d2, &v[2], &iers);
	goto L100;

	/*           compute the chebyshev moments by means of */
	/*           forward recursion. */

L80:
	an = 3.f;
	for (i__ = 3; i__ <= 12; ++i__) {
		an2 = an * an;
		v[i__ - 1] = ((an2 - 4.f) * ((par22 - an2 - an2) * 2.f * v[i__ - 2] + 
											  as) + ac - par2 * (an + 1.f) * (an + 2.f) * v[i__ - 3]) / (
												  par2 * (an - 1.f) * (an - 2.f));
		an += 2.f;
		/* L90: */
	}
L100:
	for (j = 1; j <= 12; ++j) {
		chebmo[m + (j << 1) * chebmo_dim1] = v[j - 1];
		/* L110: */
	}
L120:
	if (*nrmom < *momcom) {
		m = *nrmom + 1;
	}
	if (*momcom < *maxp1 - 1 && *nrmom >= *momcom) {
		++(*momcom);
	}

	/*           compute the coefficients of the chebyshev expansions */
	/*           of degrees 12 and 24 of the function f. */

	r__1 = centr + hlgth;
	fval[0] = f(r__1) * .5f;
	fval[12] = f(centr);
	r__1 = centr - hlgth;
	fval[24] = f(r__1) * .5f;
	for (i__ = 2; i__ <= 12; ++i__) {
		isym = 26 - i__;
		r__1 = hlgth * x[i__ - 2] + centr;
		fval[i__ - 1] = f(r__1);
		r__1 = centr - hlgth * x[i__ - 2];
		fval[isym - 1] = f(r__1);
		/* L130: */
	}
	qcheb_(x, fval, cheb12, cheb24);

	/*           compute the integral and error estimates. */

	resc12 = cheb12[12] * chebmo[m + chebmo_dim1 * 13];
	ress12 = 0.f;
	k = 11;
	for (j = 1; j <= 6; ++j) {
		resc12 += cheb12[k - 1] * chebmo[m + k * chebmo_dim1];
		ress12 += cheb12[k] * chebmo[m + (k + 1) * chebmo_dim1];
		k += -2;
		/* L140: */
	}
	resc24 = cheb24[24] * chebmo[m + chebmo_dim1 * 25];
	ress24 = 0.f;
	*resabs = fabs(cheb24[24]);
	k = 23;
	for (j = 1; j <= 12; ++j) {
		resc24 += cheb24[k - 1] * chebmo[m + k * chebmo_dim1];
		ress24 += cheb24[k] * chebmo[m + (k + 1) * chebmo_dim1];
		*resabs = (r__1 = cheb24[k - 1], fabs(r__1)) + (r__2 = cheb24[k], 
																		fabs(r__2));
		k += -2;
		/* L150: */
	}
	estc = (r__1 = resc24 - resc12, fabs(r__1));
	ests = (r__1 = ress24 - ress12, fabs(r__1));
	*resabs *= fabs(hlgth);
	if (*integr == 2) {
		goto L160;
	}
	*result = conc * resc24 - cons * ress24;
	*abserr = (r__1 = conc * estc, fabs(r__1)) + (r__2 = cons * ests, fabs(
																	 r__2));
	goto L170;
L160:
	*result = conc * ress24 + cons * resc24;
	*abserr = (r__1 = conc * ests, fabs(r__1)) + (r__2 = cons * estc, fabs(
																	 r__2));
L170:
	return;
} /* qc25f_ */

void qc25s_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *bl, const sys_float *br, 
				const sys_float *alfa, const sys_float *beta, 
				sys_float *ri, sys_float *rj, sys_float *rg, sys_float *rh, 
				sys_float *result, sys_float *abserr, sys_float *resasc, 
				const long *integr, long *nev)
{
	/* Initialized data */

	sys_float x[11] = { .9914448613738104f,.9659258262890683f,
								  .9238795325112868f,.8660254037844386f,.7933533402912352f,
								  .7071067811865475f,.6087614290087206f,.5f,.3826834323650898f,
								  .2588190451025208f,.1305261922200516f };

	/* System generated locals */
	sys_float r__1;
	double d__1, d__2;

	/* Local variables */
	long i__;
	sys_float u, dc, fix, fval[25], res12, res24;
	long isym;
	sys_float cheb12[13], cheb24[25];
	sys_float hlgth, centr;
	sys_float factor, resabs;

	/* ***begin prologue  qc25s */
	/* ***date written   810101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a2 */
	/* ***keywords  25-point clenshaw-curtis integration */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f*w over (bl,br), with error */
	/*            estimate, where the weight function w has a singular */
	/*            behaviour of algebraico-logarithmic type at the points */
	/*            a and/or b. (bl,br) is a part of (a,b). */
	/* ***description */

	/*        integration rules for integrands having algebraico-logarithmic */
	/*        end point singularities */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*           f      - sys_float */
	/*                    function subprogram defining the integrand */
	/*                    f(x). the actual name for f needs to be declared */
	/*                    e x t e r n a l  in the driver program. */

	/*           a      - sys_float */
	/*                    left end point of the original interval */

	/*           b      - sys_float */
	/*                    right end point of the original interval, b.gt.a */

	/*           bl     - sys_float */
	/*                    lower limit of integration, bl.ge.a */

	/*           br     - sys_float */
	/*                    upper limit of integration, br.le.b */

	/*           alfa   - sys_float */
	/*                    parameter in the weight function */

	/*           beta   - sys_float */
	/*                    parameter in the weight function */

	/*           ri,rj,rg,rh - sys_float */
	/*                    modified chebyshev moments for the application */
	/*                    of the generalized clenshaw-curtis */
	/*                    method (computed in subroutine dqmomo) */

	/*           result - sys_float */
	/*                    approximation to the integral */
	/*                    result is computed by using a generalized */
	/*                    clenshaw-curtis method if b1 = a or br = b. */
	/*                    in all other cases the 15-point kronrod */
	/*                    rule is applied, obtained by optimal addition of */
	/*                    abscissae to the 7-point gauss rule. */

	/*           abserr - sys_float */
	/*                    estimate of the modulus of the absolute error, */
	/*                    which should equal or exceed fabs(i-result) */

	/*           resasc - sys_float */
	/*                    approximation to the integral of fabs(f*w-i/(b-a)) */

	/*           integr - long */
	/*                    which determines the weight function */
	/*                    = 1   w(x) = (x-a)**alfa*(b-x)**beta */
	/*                    = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a) */
	/*                    = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x) */
	/*                    = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)* */
	/*                                 log(b-x) */

	/*           nev    - long */
	/*                    number of integrand evaluations */

	/* ***references  (none) */
	/* ***routines called  qcheb,qk15w */
	/* ***end prologue  qc25s */




	/*           the vector x contains the values cos(k*pi/24) */
	/*           k = 1, ..., 11, to be used for the computation of the */
	/*           chebyshev series expansion of f. */

	/* Parameter adjustments */
	--rh;
	--rg;
	--rj;
	--ri;

	/* Function Body */

	/*           list of major variables */
	/*           ----------------------- */

	/*           fval   - value of the function f at the points */
	/*                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5 */
	/*                    k = 0, ..., 24 */
	/*           cheb12 - coefficients of the chebyshev series expansion */
	/*                    of degree 12, for the function f, in the */
	/*                    interval (bl,br) */
	/*           cheb24 - coefficients of the chebyshev series expansion */
	/*                    of degree 24, for the function f, in the */
	/*                    interval (bl,br) */
	/*           res12  - approximation to the integral obtained from cheb12 */
	/*           res24  - approximation to the integral obtained from cheb24 */
	/*           qwgts - external function subprogram defining */
	/*                    the four possible weight functions */
	/*           hlgth  - half-length of the interval (bl,br) */
	/*           centr  - mid point of the interval (bl,br) */

	/* ***first executable statement  qc25s */
	*nev = 25;
	if (*bl == *a && (*alfa != 0.f || *integr == 2 || *integr == 4)) {
		goto L10;
	}
	if (*br == *b && (*beta != 0.f || *integr == 3 || *integr == 4)) {
		goto L140;
	}

	/*           if a.gt.bl and b.lt.br, apply the 15-point gauss-kronrod */
	/*           scheme. */


	qk15w_(f, qwgts_, a, b, alfa, beta, integr, bl, br, result, 
			 abserr, &resabs, resasc);
	*nev = 15;
	goto L270;

	/*           this part of the program is executed only if a = bl. */
	/*           ---------------------------------------------------- */

	/*           compute the chebyshev series expansion of the */
	/*           following function */
	/*           f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta */
	/*                  *f(0.5*(br-a)*x+0.5*(br+a)) */

L10:
	hlgth = (*br - *bl) * .5f;
	centr = (*br + *bl) * .5f;
	fix = *b - centr;
	r__1 = hlgth + centr;
	d__1 = (double) (fix - hlgth);
	d__2 = (double) (*beta);
	fval[0] = f(r__1) * .5f * pow(d__1, d__2);
	d__1 = (double) fix;
	d__2 = (double) (*beta);
	fval[12] = f(centr) * pow(d__1, d__2);
	r__1 = centr - hlgth;
	d__1 = (double) (fix + hlgth);
	d__2 = (double) (*beta);
	fval[24] = f(r__1) * .5f * pow(d__1, d__2);
	for (i__ = 2; i__ <= 12; ++i__) {
		u = hlgth * x[i__ - 2];
		isym = 26 - i__;
		r__1 = u + centr;
		d__1 = (double) (fix - u);
		d__2 = (double) (*beta);
		fval[i__ - 1] = f(r__1) * pow(d__1, d__2);
		r__1 = centr - u;
		d__1 = (double) (fix + u);
		d__2 = (double) (*beta);
		fval[isym - 1] = f(r__1) * pow(d__1, d__2);
		/* L20: */
	}
	d__1 = (double) hlgth;
	d__2 = (double) (*alfa + 1.f);
	factor = pow(d__1, d__2);
	*result = 0.f;
	*abserr = 0.f;
	res12 = 0.f;
	res24 = 0.f;
	if (*integr > 2) {
		goto L70;
	}
	qcheb_(x, fval, cheb12, cheb24);

	/*           integr = 1  (or 2) */

	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * ri[i__];
		res24 += cheb24[i__ - 1] * ri[i__];
		/* L30: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * ri[i__];
		/* L40: */
	}
	if (*integr == 1) {
		goto L130;
	}

	/*           integr = 2 */

	dc = log(*br - *bl);
	*result = res24 * dc;
	*abserr = (r__1 = (res24 - res12) * dc, fabs(r__1));
	res12 = 0.f;
	res24 = 0.f;
	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rg[i__];
		res24 = res12 + cheb24[i__ - 1] * rg[i__];
		/* L50: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rg[i__];
		/* L60: */
	}
	goto L130;

	/*           compute the chebyshev series expansion of the */
	/*           following function */
	/*           f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x) */

L70:
	fval[0] *= log(fix - hlgth);
	fval[12] *= log(fix);
	fval[24] *= log(fix + hlgth);
	for (i__ = 2; i__ <= 12; ++i__) {
		u = hlgth * x[i__ - 2];
		isym = 26 - i__;
		fval[i__ - 1] *= log(fix - u);
		fval[isym - 1] *= log(fix + u);
		/* L80: */
	}
	qcheb_(x, fval, cheb12, cheb24);

	/*           integr = 3  (or 4) */

	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * ri[i__];
		res24 += cheb24[i__ - 1] * ri[i__];
		/* L90: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * ri[i__];
		/* L100: */
	}
	if (*integr == 3) {
		goto L130;
	}

	/*           integr = 4 */

	dc = log(*br - *bl);
	*result = res24 * dc;
	*abserr = (r__1 = (res24 - res12) * dc, fabs(r__1));
	res12 = 0.f;
	res24 = 0.f;
	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rg[i__];
		res24 += cheb24[i__ - 1] * rg[i__];
		/* L110: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rg[i__];
		/* L120: */
	}
L130:
	*result = (*result + res24) * factor;
	*abserr = (*abserr + (r__1 = res24 - res12, fabs(r__1))) * factor;
	goto L270;

	/*           this part of the program is executed only if b = br. */
	/*           ---------------------------------------------------- */

	/*           compute the chebyshev series expansion of the */
	/*           following function */
	/*           f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa */
	/*                *f(0.5*(b-bl)*x+0.5*(b+bl)) */

L140:
	hlgth = (*br - *bl) * .5f;
	centr = (*br + *bl) * .5f;
	fix = centr - *a;
	r__1 = hlgth + centr;
	d__1 = (double) (fix + hlgth);
	d__2 = (double) (*alfa);
	fval[0] = f(r__1) * .5f * pow(d__1, d__2);
	d__1 = (double) fix;
	d__2 = (double) (*alfa);
	fval[12] = f(centr) * pow(d__1, d__2);
	r__1 = centr - hlgth;
	d__1 = (double) (fix - hlgth);
	d__2 = (double) (*alfa);
	fval[24] = f(r__1) * .5f * pow(d__1, d__2);
	for (i__ = 2; i__ <= 12; ++i__) {
		u = hlgth * x[i__ - 2];
		isym = 26 - i__;
		r__1 = u + centr;
		d__1 = (double) (fix + u);
		d__2 = (double) (*alfa);
		fval[i__ - 1] = f(r__1) * pow(d__1, d__2);
		r__1 = centr - u;
		d__1 = (double) (fix - u);
		d__2 = (double) (*alfa);
		fval[isym - 1] = f(r__1) * pow(d__1, d__2);
		/* L150: */
	}
	d__1 = (double) hlgth;
	d__2 = (double) (*beta + 1.f);
	factor = pow(d__1, d__2);
	*result = 0.f;
	*abserr = 0.f;
	res12 = 0.f;
	res24 = 0.f;
	if (*integr == 2 || *integr == 4) {
		goto L200;
	}

	/*           integr = 1  (or 3) */

	qcheb_(x, fval, cheb12, cheb24);
	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rj[i__];
		res24 += cheb24[i__ - 1] * rj[i__];
		/* L160: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rj[i__];
		/* L170: */
	}
	if (*integr == 1) {
		goto L260;
	}

	/*           integr = 3 */

	dc = log(*br - *bl);
	*result = res24 * dc;
	*abserr = (r__1 = (res24 - res12) * dc, fabs(r__1));
	res12 = 0.f;
	res24 = 0.f;
	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rh[i__];
		res24 += cheb24[i__ - 1] * rh[i__];
		/* L180: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rh[i__];
		/* L190: */
	}
	goto L260;

	/*           compute the chebyshev series expansion of the */
	/*           following function */
	/*           f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a)) */

L200:
	fval[0] *= log(hlgth + fix);
	fval[12] *= log(fix);
	fval[24] *= log(fix - hlgth);
	for (i__ = 2; i__ <= 12; ++i__) {
		u = hlgth * x[i__ - 2];
		isym = 26 - i__;
		fval[i__ - 1] *= log(u + fix);
		fval[isym - 1] *= log(fix - u);
		/* L210: */
	}
	qcheb_(x, fval, cheb12, cheb24);

	/*           integr = 2  (or 4) */

	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rj[i__];
		res24 += cheb24[i__ - 1] * rj[i__];
		/* L220: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rj[i__];
		/* L230: */
	}
	if (*integr == 2) {
		goto L260;
	}
	dc = log(*br - *bl);
	*result = res24 * dc;
	*abserr = (r__1 = (res24 - res12) * dc, fabs(r__1));
	res12 = 0.f;
	res24 = 0.f;

	/*           integr = 4 */

	for (i__ = 1; i__ <= 13; ++i__) {
		res12 += cheb12[i__ - 1] * rh[i__];
		res24 += cheb24[i__ - 1] * rh[i__];
		/* L240: */
	}
	for (i__ = 14; i__ <= 25; ++i__) {
		res24 += cheb24[i__ - 1] * rh[i__];
		/* L250: */
	}
L260:
	*result = (*result + res24) * factor;
	*abserr = (*abserr + (r__1 = res24 - res12, fabs(r__1))) * factor;
L270:
	return;
} /* qc25s_ */

void qcheb_(sys_float *x, sys_float *fval, sys_float *cheb12, sys_float *cheb24)
{
	long i__, j;
	sys_float v[12], alam, alam1, alam2, part1, part2, part3;

	/* ***begin prologue  qcheb */
	/* ***refer to  qc25c,qc25f,qc25s */
	/* ***routines called  (none) */
	/* ***revision date  830518   (yymmdd) */
	/* ***keywords  chebyshev series expansion, fast fourier transform */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  this routine computes the chebyshev series expansion */
	/*            of degrees 12 and 24 of a function using a */
	/*            fast fourier transform method */
	/*            f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)), */
	/*            f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)), */
	/*            where t(k,x) is the chebyshev polynomial of degree k. */
	/* ***description */

	/*        chebyshev series expansion */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*          on entry */
	/*           x      - sys_float */
	/*                    vector of dimension 11 containing the */
	/*                    values cos(k*pi/24), k = 1, ..., 11 */

	/*           fval   - sys_float */
	/*                    vector of dimension 25 containing the */
	/*                    function values at the points */
	/*                    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24, */
	/*                    where (a,b) is the approximation interval. */
	/*                    fval(1) and fval(25) are divided by two */
	/*                    (these values are destroyed at output). */

	/*          on return */
	/*           cheb12 - sys_float */
	/*                    vector of dimension 13 containing the */
	/*                    chebyshev coefficients for degree 12 */

	/*           cheb24 - sys_float */
	/*                    vector of dimension 25 containing the */
	/*                    chebyshev coefficients for degree 24 */

	/* ***end prologue  qcheb */



	/* ***first executable statement  qcheb */
	/* Parameter adjustments */
	--cheb24;
	--cheb12;
	--fval;
	--x;

	/* Function Body */
	for (i__ = 1; i__ <= 12; ++i__) {
		j = 26 - i__;
		v[i__ - 1] = fval[i__] - fval[j];
		fval[i__] += fval[j];
		/* L10: */
	}
	alam1 = v[0] - v[8];
	alam2 = x[6] * (v[2] - v[6] - v[10]);
	cheb12[4] = alam1 + alam2;
	cheb12[10] = alam1 - alam2;
	alam1 = v[1] - v[7] - v[9];
	alam2 = v[3] - v[5] - v[11];
	alam = x[3] * alam1 + x[9] * alam2;
	cheb24[4] = cheb12[4] + alam;
	cheb24[22] = cheb12[4] - alam;
	alam = x[9] * alam1 - x[3] * alam2;
	cheb24[10] = cheb12[10] + alam;
	cheb24[16] = cheb12[10] - alam;
	part1 = x[4] * v[4];
	part2 = x[8] * v[8];
	part3 = x[6] * v[6];
	alam1 = v[0] + part1 + part2;
	alam2 = x[2] * v[2] + part3 + x[10] * v[10];
	cheb12[2] = alam1 + alam2;
	cheb12[12] = alam1 - alam2;
	alam = x[1] * v[1] + x[3] * v[3] + x[5] * v[5] + x[7] * v[7] + x[9] * v[9]
		+ x[11] * v[11];
	cheb24[2] = cheb12[2] + alam;
	cheb24[24] = cheb12[2] - alam;
	alam = x[11] * v[1] - x[9] * v[3] + x[7] * v[5] - x[5] * v[7] + x[3] * v[
		9] - x[1] * v[11];
	cheb24[12] = cheb12[12] + alam;
	cheb24[14] = cheb12[12] - alam;
	alam1 = v[0] - part1 + part2;
	alam2 = x[10] * v[2] - part3 + x[2] * v[10];
	cheb12[6] = alam1 + alam2;
	cheb12[8] = alam1 - alam2;
	alam = x[5] * v[1] - x[9] * v[3] - x[1] * v[5] - x[11] * v[7] + x[3] * v[
		9] + x[7] * v[11];
	cheb24[6] = cheb12[6] + alam;
	cheb24[20] = cheb12[6] - alam;
	alam = x[7] * v[1] - x[3] * v[3] - x[11] * v[5] + x[1] * v[7] - x[9] * v[
		9] - x[5] * v[11];
	cheb24[8] = cheb12[8] + alam;
	cheb24[18] = cheb12[8] - alam;
	for (i__ = 1; i__ <= 6; ++i__) {
		j = 14 - i__;
		v[i__ - 1] = fval[i__] - fval[j];
		fval[i__] += fval[j];
		/* L20: */
	}
	alam1 = v[0] + x[8] * v[4];
	alam2 = x[4] * v[2];
	cheb12[3] = alam1 + alam2;
	cheb12[11] = alam1 - alam2;
	cheb12[7] = v[0] - v[4];
	alam = x[2] * v[1] + x[6] * v[3] + x[10] * v[5];
	cheb24[3] = cheb12[3] + alam;
	cheb24[23] = cheb12[3] - alam;
	alam = x[6] * (v[1] - v[3] - v[5]);
	cheb24[7] = cheb12[7] + alam;
	cheb24[19] = cheb12[7] - alam;
	alam = x[10] * v[1] - x[6] * v[3] + x[2] * v[5];
	cheb24[11] = cheb12[11] + alam;
	cheb24[15] = cheb12[11] - alam;
	for (i__ = 1; i__ <= 3; ++i__) {
		j = 8 - i__;
		v[i__ - 1] = fval[i__] - fval[j];
		fval[i__] += fval[j];
		/* L30: */
	}
	cheb12[5] = v[0] + x[8] * v[2];
	cheb12[9] = fval[1] - x[8] * fval[3];
	alam = x[4] * v[1];
	cheb24[5] = cheb12[5] + alam;
	cheb24[21] = cheb12[5] - alam;
	alam = x[8] * fval[2] - fval[4];
	cheb24[9] = cheb12[9] + alam;
	cheb24[17] = cheb12[9] - alam;
	cheb12[1] = fval[1] + fval[3];
	alam = fval[2] + fval[4];
	cheb24[1] = cheb12[1] + alam;
	cheb24[25] = cheb12[1] - alam;
	cheb12[13] = v[0] - v[2];
	cheb24[13] = cheb12[13];
	alam = .16666666666666666f;
	for (i__ = 2; i__ <= 12; ++i__) {
		cheb12[i__] *= alam;
		/* L40: */
	}
	alam *= .5f;
	cheb12[1] *= alam;
	cheb12[13] *= alam;
	for (i__ = 2; i__ <= 24; ++i__) {
		cheb24[i__] *= alam;
		/* L50: */
	}
	cheb24[1] = alam * .5f * cheb24[1];
	cheb24[25] = alam * .5f * cheb24[25];
	return;
} /* qcheb_ */

void qelg_(long *n, sys_float *epstab, sys_float *result, sys_float *
			  abserr, sys_float *res3la, long *nres)
{
	/* System generated locals */
	int i__1;
	sys_float r__1, r__2, r__3;

	/* Local variables */
	long i__;
	sys_float e0, e1, e2, e3;
	long k1, k2, k3, ib, ie;
	sys_float ss;
	long ib2;
	sys_float res;
	long num;
	sys_float err1, err2, err3, tol1, tol2, tol3;
	long indx;
	sys_float e1abs, oflow, error, delta1, delta2, delta3;
	sys_float epmach, epsinf;
	long newelm, limexp;

	/* ***begin prologue  qelg */
	/* ***refer to  qagie,qagoe,qagpe,qagse */
	/* ***routines called  r1mach */
	/* ***revision date  830518   (yymmdd) */
	/* ***keywords  epsilon algorithm, convergence acceleration, */
	/*             extrapolation */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
	/* ***purpose  the routine determines the limit of a given sequence of */
	/*            approximations, by means of the epsilon algorithm of */
	/*            p. wynn. an estimate of the absolute error is also given. */
	/*            the condensed epsilon table is computed. only those */
	/*            elements needed for the computation of the next diagonal */
	/*            are preserved. */
	/* ***description */

	/*           epsilon algorithm */
	/*           standard fortran subroutine */
	/*           sys_float version */

	/*           parameters */
	/*              n      - long */
	/*                       epstab(n) contains the new element in the */
	/*                       first column of the epsilon table. */

	/*              epstab - sys_float */
	/*                       vector of dimension 52 containing the elements */
	/*                       of the two lower diagonals of the triangular */
	/*                       epsilon table. the elements are numbered */
	/*                       starting at the right-hand corner of the */
	/*                       triangle. */

	/*              result - sys_float */
	/*                       resulting approximation to the integral */

	/*              abserr - sys_float */
	/*                       estimate of the absolute error computed from */
	/*                       result and the 3 previous results */

	/*              res3la - sys_float */
	/*                       vector of dimension 3 containing the last 3 */
	/*                       results */

	/*              nres   - long */
	/*                       number of calls to the routine */
	/*                       (should be zero at first call) */

	/* ***end prologue  qelg */


	/*           list of major variables */
	/*           ----------------------- */

	/*           e0     - the 4 elements on which the */
	/*           e1       computation of a new element in */
	/*           e2       the epsilon table is based */
	/*           e3                 e0 */
	/*                        e3    e1    new */
	/*                              e2 */
	/*           newelm - number of elements to be computed in the new */
	/*                    diagonal */
	/*           error  - error = fabs(e1-e0)+fabs(e2-e1)+fabs(new-e2) */
	/*           result - the element in the new diagonal with least value */
	/*                    of error */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           oflow is the largest positive magnitude. */
	/*           limexp is the maximum number of elements the epsilon */
	/*           table can contain. if this number is reached, the upper */
	/*           diagonal of the epsilon table is deleted. */

	/* ***first executable statement  qelg */
	/* Parameter adjustments */
	--res3la;
	--epstab;

	/* Function Body */
	epmach = r1mach(c__4);
	oflow = r1mach(c__2);
	++(*nres);
	*abserr = oflow;
	*result = epstab[*n];
	if (*n < 3) {
		goto L100;
	}
	limexp = 50;
	epstab[*n + 2] = epstab[*n];
	newelm = (*n - 1) / 2;
	epstab[*n] = oflow;
	num = *n;
	k1 = *n;
	i__1 = newelm;
	for (i__ = 1; i__ <= i__1; ++i__) {
		k2 = k1 - 1;
		k3 = k1 - 2;
		res = epstab[k1 + 2];
		e0 = epstab[k3];
		e1 = epstab[k2];
		e2 = res;
		e1abs = fabs(e1);
		delta2 = e2 - e1;
		err2 = fabs(delta2);
		/* Computing MAX */
		r__1 = fabs(e2);
		tol2 = max(r__1,e1abs) * epmach;
		delta3 = e1 - e0;
		err3 = fabs(delta3);
		/* Computing MAX */
		r__1 = e1abs, r__2 = fabs(e0);
		tol3 = max(r__1,r__2) * epmach;
		if (err2 > tol2 || err3 > tol3) {
			goto L10;
		}

		/*           if e0, e1 and e2 are equal to within machine */
		/*           accuracy, convergence is assumed. */
		/*           result = e2 */
		/*           abserr = fabs(e1-e0)+fabs(e2-e1) */

		*result = res;
		*abserr = err2 + err3;
		/* ***jump out of do-loop */
		goto L100;
	L10:
		e3 = epstab[k1];
		epstab[k1] = e1;
		delta1 = e1 - e3;
		err1 = fabs(delta1);
		/* Computing MAX */
		r__1 = e1abs, r__2 = fabs(e3);
		tol1 = max(r__1,r__2) * epmach;

		/*           if two elements are very close to each other, omit */
		/*           a part of the table by adjusting the value of n */

		if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
			goto L20;
		}
		ss = 1.f / delta1 + 1.f / delta2 - 1.f / delta3;
		epsinf = (r__1 = ss * e1, fabs(r__1));

		/*           test to detect irregular behaviour in the table, and */
		/*           eventually omit a part of the table adjusting the value */
		/*           of n. */

		if (epsinf > 1e-4f) {
			goto L30;
		}
	L20:
		*n = i__ + i__ - 1;
		/* ***jump out of do-loop */
		goto L50;

		/*           compute a new element and eventually adjust */
		/*           the value of result. */

	L30:
		res = e1 + 1.f / ss;
		epstab[k1] = res;
		k1 += -2;
		error = err2 + (r__1 = res - e2, fabs(r__1)) + err3;
		if (error > *abserr) {
			goto L40;
		}
		*abserr = error;
		*result = res;
	L40:
		;
	}

	/*           shift the table. */

L50:
	if (*n == limexp) {
		*n = (limexp / 2 << 1) - 1;
	}
	ib = 1;
	if (num / 2 << 1 == num) {
		ib = 2;
	}
	ie = newelm + 1;
	i__1 = ie;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ib2 = ib + 2;
		epstab[ib] = epstab[ib2];
		ib = ib2;
		/* L60: */
	}
	if (num == *n) {
		goto L80;
	}
	indx = num - *n + 1;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		epstab[i__] = epstab[indx];
		++indx;
		/* L70: */
	}
L80:
	if (*nres >= 4) {
		goto L90;
	}
	res3la[*nres] = *result;
	*abserr = oflow;
	goto L100;

	/*           compute error estimate */

L90:
	*abserr = (r__1 = *result - res3la[3], fabs(r__1)) + (r__2 = *result - 
																			res3la[2], fabs(r__2)) + (r__3 = *result - res3la[1], fabs(r__3));
	res3la[1] = res3la[2];
	res3la[2] = res3la[3];
	res3la[3] = *result;
L100:
	/* Computing MAX */
	r__1 = *abserr, r__2 = epmach * 5.f * fabs(*result);
	*abserr = max(r__1,r__2);
	return;
} /* qelg_ */

void qk15_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  sys_float *result, sys_float *abserr, 
			  sys_float *resabs, sys_float *resasc)
{
	/* Initialized data */

	sys_float xgk[8] = { .9914553711208126f,.9491079123427585f,
									.8648644233597691f,.7415311855993944f,.5860872354676911f,
									.4058451513773972f,.2077849550078985f,0.f };
	sys_float wgk[8] = { .02293532201052922f,.06309209262997855f,
									.1047900103222502f,.1406532597155259f,.1690047266392679f,
									.1903505780647854f,.2044329400752989f,.2094821410847278f };
	sys_float wg[4] = { .1294849661688697f,.2797053914892767f,
								  .3818300505051189f,.4179591836734694f };

	/* System generated locals */
	sys_float r__1, r__2;
	double d__1;

	/* Local variables */
	long j;
	sys_float fc, fv1[7], fv2[7];
	long jtw;
	sys_float absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	sys_float hlgth, centr, reskh, uflow;
	sys_float epmach, dhlgth;

	/* ***begin prologue  qk15 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  15-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b), with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           sys_float version */

	/*           parameters */
	/*            on entry */
	/*              f      - sys_float */
	/*                       function subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the calling program. */

	/*              a      - sys_float */
	/*                       lower limit of integration */

	/*              b      - sys_float */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - sys_float */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 15-polong */
	/*                       kronrod rule (resk) obtained by optimal addition */
	/*                       of abscissae to the7-point gauss rule(resg). */

	/*              abserr - sys_float */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should not exceed fabs(i-result) */

	/*              resabs - sys_float */
	/*                       approximation to the integral j */

	/*              resasc - sys_float */
	/*                       approximation to the integral of fabs(f-i/(b-a)) */
	/*                       over (a,b) */

	/* ***references  (none) */
	/* ***routines called  r1mach */
	/* ***end prologue  qk15 */



	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 15-point kronrod rule */
	/*                    xgk(2), xgk(4), ...  abscissae of the 7-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 7-point gauss rule */

	/*           wgk    - weights of the 15-point kronrod rule */

	/*           wg     - weights of the 7-point gauss rule */



	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 7-point gauss formula */
	/*           resk   - result of the 15-point kronrod formula */
	/*           reskh  - approximation to the mean value of f over (a,b), */
	/*                    i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qk15 */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);

	centr = (*a + *b) * .5f;
	hlgth = (*b - *a) * .5f;
	dhlgth = fabs(hlgth);

	/*           compute the 15-point kronrod approximation to */
	/*           the integral, and estimate the absolute error. */

	fc = f(centr);
	resg = fc * wg[3];
	resk = fc * wgk[7];
	*resabs = fabs(resk);
	for (j = 1; j <= 3; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 4; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5f;
	*resasc = wgk[7] * (r__1 = fc - reskh, fabs(r__1));
	for (j = 1; j <= 7; ++j) {
		*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, fabs(r__1)) + (
											 r__2 = fv2[j - 1] - reskh, fabs(r__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (r__1 = (resk - resg) * hlgth, fabs(r__1));
	if (*resasc != 0.f && *abserr != 0.f) {
		/* Computing MIN */
		d__1 = (double) (*abserr * 200.f / *resasc);
		r__1 = 1.f, r__2 = pow(d__1, c_b270);
		*abserr = *resasc * min(r__1,r__2);
	}
	if (*resabs > uflow / (epmach * 50.f)) {
		/* Computing MAX */
		r__1 = epmach * 50.f * *resabs;
		*abserr = max(r__1,*abserr);
	}
	return;
} /* qk15_ */

void qk15i_(const E_fp& f, const sys_float *boun, const long *inf, const sys_float *a,
				const sys_float *b, sys_float *result, sys_float *abserr, 
				sys_float *resabs, sys_float *resasc)
{
	/* Initialized data */

	sys_float xgk[8] = { .9914553711208126f,.9491079123427585f,
									.8648644233597691f,.7415311855993944f,.5860872354676911f,
									.4058451513773972f,.2077849550078985f,0.f };
	sys_float wgk[8] = { .02293532201052922f,.06309209262997855f,
									.1047900103222502f,.1406532597155259f,.1690047266392679f,
									.1903505780647854f,.2044329400752989f,.2094821410847278f };
	sys_float wg[8] = { 0.f,.1294849661688697f,0.f,.2797053914892767f,0.f,
								  .3818300505051189f,0.f,.4179591836734694f };

	/* System generated locals */
	sys_float r__1, r__2;
	double d__1;

	/* Local variables */
	long j;
	sys_float fc, fv1[7], fv2[7], absc, dinf, resg, resk, fsum, absc1, 
		absc2, fval1, fval2, hlgth, centr, reskh, uflow;
	sys_float tabsc1, tabsc2, epmach;

	/* ***begin prologue  qk15i */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a3a2,h2a4a2 */
	/* ***keywords  15-point transformed gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  the original (infinite integration range is mapped */
	/*            onto the interval (0,1) and (a,b) is a part of (0,1). */
	/*            it is the purpose to compute */
	/*            i = integral of transformed integrand over (a,b), */
	/*            j = integral of fabs(transformed integrand) over (a,b). */
	/* ***description */

	/*           integration rule */
	/*           standard fortran subroutine */
	/*           sys_float version */

	/*           parameters */
	/*            on entry */
	/*              f      - sys_float */
	/*                       fuction subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the calling program. */

	/*              boun   - sys_float */
	/*                       finite bound of original integration */
	/*                       range (set to zero if inf = +2) */

	/*              inf    - long */
	/*                       if inf = -1, the original interval is */
	/*                                   (-infinity,bound), */
	/*                       if inf = +1, the original interval is */
	/*                                   (bound,+infinity), */
	/*                       if inf = +2, the original interval is */
	/*                                   (-infinity,+infinity) and */
	/*                       the integral is computed as the sum of two */
	/*                       integrals, one over (-infinity,0) and one over */
	/*                       (0,+infinity). */

	/*              a      - sys_float */
	/*                       lower limit for integration over subrange */
	/*                       of (0,1) */

	/*              b      - sys_float */
	/*                       upper limit for integration over subrange */
	/*                       of (0,1) */

	/*            on return */
	/*              result - sys_float */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 15-polong */
	/*                       kronrod rule(resk) obtained by optimal addition */
	/*                       of abscissae to the 7-point gauss rule(resg). */

	/*              abserr - sys_float */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should equal or exceed fabs(i-result) */

	/*              resabs - sys_float */
	/*                       approximation to the integral j */

	/*              resasc - sys_float */
	/*                       approximation to the integral of */
	/*                       fabs((transformed integrand)-i/(b-a)) over (a,b) */

	/* ***references  (none) */
	/* ***routines called  r1mach */
	/* ***end prologue  qk15i */



	/*           the abscissae and weights are supplied for the interval */
	/*           (-1,1).  because of symmetry only the positive abscissae and */
	/*           their corresponding weights are given. */

	/*           xgk    - abscissae of the 15-point kronrod rule */
	/*                    xgk(2), xgk(4), ... abscissae of the 7-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 7-point gauss rule */

	/*           wgk    - weights of the 15-point kronrod rule */

	/*           wg     - weights of the 7-point gauss rule, corresponding */
	/*                    to the abscissae xgk(2), xgk(4), ... */
	/*                    wg(1), wg(3), ... are set to zero. */





	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc*  - abscissa */
	/*           tabsc* - transformed abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 7-point gauss formula */
	/*           resk   - result of the 15-point kronrod formula */
	/*           reskh  - approximation to the mean value of the transformed */
	/*                    integrand over (a,b), i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qk15i */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);
	dinf = (sys_float) min(1,*inf);

	centr = (*a + *b) * .5f;
	hlgth = (*b - *a) * .5f;
	tabsc1 = *boun + dinf * (1.f - centr) / centr;
	fval1 = f(tabsc1);
	if (*inf == 2) {
		r__1 = -tabsc1;
		fval1 += f(r__1);
	}
	fc = fval1 / centr / centr;

	/*           compute the 15-point kronrod approximation to */
	/*           the integral, and estimate the error. */

	resg = wg[7] * fc;
	resk = wgk[7] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 7; ++j) {
		absc = hlgth * xgk[j - 1];
		absc1 = centr - absc;
		absc2 = centr + absc;
		tabsc1 = *boun + dinf * (1.f - absc1) / absc1;
		tabsc2 = *boun + dinf * (1.f - absc2) / absc2;
		fval1 = f(tabsc1);
		fval2 = f(tabsc2);
		if (*inf == 2) {
			r__1 = -tabsc1;
			fval1 += f(r__1);
		}
		if (*inf == 2) {
			r__1 = -tabsc2;
			fval2 += f(r__1);
		}
		fval1 = fval1 / absc1 / absc1;
		fval2 = fval2 / absc2 / absc2;
		fv1[j - 1] = fval1;
		fv2[j - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[j - 1] * fsum;
		*resabs += wgk[j - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	reskh = resk * .5f;
	*resasc = wgk[7] * (r__1 = fc - reskh, fabs(r__1));
	for (j = 1; j <= 7; ++j) {
		*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, fabs(r__1)) + (
											 r__2 = fv2[j - 1] - reskh, fabs(r__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resasc *= hlgth;
	*resabs *= hlgth;
	*abserr = (r__1 = (resk - resg) * hlgth, fabs(r__1));
	if (*resasc != 0.f && *abserr != 0.f) {
		/* Computing MIN */
		d__1 = (double) (*abserr * 200.f / *resasc);
		r__1 = 1.f, r__2 = pow(d__1, c_b270);
		*abserr = *resasc * min(r__1,r__2);
	}
	if (*resabs > uflow / (epmach * 50.f)) {
		/* Computing MAX */
		r__1 = epmach * 50.f * *resabs;
		*abserr = max(r__1,*abserr);
	}
	return;
} /* qk15i_ */

void qk15w_(const E_fp& f, E_fp1 w, const sys_float *p1, const sys_float *p2, 
				const sys_float *p3, const sys_float *p4, const long *kp,
				const sys_float *a, const sys_float *b, sys_float *result, 
				sys_float *abserr, sys_float *resabs, sys_float *resasc)
{
	/* Initialized data */

	sys_float xgk[8] = { .9914553711208126f,.9491079123427585f,
									.8648644233597691f,.7415311855993944f,.5860872354676911f,
									.4058451513773972f,.2077849550078985f,0.f };
	sys_float wgk[8] = { .02293532201052922f,.06309209262997855f,
									.1047900103222502f,.1406532597155259f,.1690047266392679f,
									.1903505780647854f,.2044329400752989f,.2094821410847278f };
	sys_float wg[4] = { .1294849661688697f,.2797053914892767f,
								  .3818300505051889f,.4179591836734694f };

	/* System generated locals */
	sys_float r__1, r__2;
	double d__1;

	/* Local variables */
	long j;
	sys_float fc, fv1[7], fv2[7];
	long jtw;
	sys_float absc, resg, resk, fsum, absc1, absc2, fval1, fval2;
	long jtwm1;
	sys_float hlgth, centr, reskh, uflow;
	sys_float epmach, dhlgth;

	/* ***begin prologue  qk15w */
	/* ***date written   810101   (yymmdd) */
	/* ***revision date  830518   (mmddyy) */
	/* ***category no.  h2a2a2 */
	/* ***keywords  15-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f*w over (a,b), with error */
	/*                           estimate */
	/*                       j = integral of fabs(f*w) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           sys_float version */

	/*           parameters */
	/*             on entry */
	/*              f      - sys_float */
	/*                       function subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the driver program. */

	/*              w      - sys_float */
	/*                       function subprogram defining the integrand */
	/*                       weight function w(x). the actual name for w */
	/*                       needs to be declared e x t e r n a l in the */
	/*                       calling program. */

	/*              p1, p2, p3, p4 - sys_float */
	/*                       parameters in the weight function */

	/*              kp     - long */
	/*                       key for indicating the type of weight function */

	/*              a      - sys_float */
	/*                       lower limit of integration */

	/*              b      - sys_float */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - sys_float */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 15-polong */
	/*                       kronrod rule (resk) obtained by optimal addition */
	/*                       of abscissae to the 7-point gauss rule (resg). */

	/*              abserr - sys_float */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should equal or exceed fabs(i-result) */

	/*              resabs - sys_float */
	/*                       approximation to the integral of fabs(f) */

	/*              resasc - sys_float */
	/*                       approximation to the integral of fabs(f-i/(b-a)) */

	/* ***references  (none) */
	/* ***routines called  r1mach */
	/* ***end prologue  qk15w */



	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 15-point gauss-kronrod rule */
	/*                    xgk(2), xgk(4), ... abscissae of the 7-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ... abscissae which are optimally */
	/*                    added to the 7-point gauss rule */

	/*           wgk    - weights of the 15-point gauss-kronrod rule */

	/*           wg     - weights of the 7-point gauss rule */





	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc*  - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 7-point gauss formula */
	/*           resk   - result of the 15-point kronrod formula */
	/*           reskh  - approximation to the mean value of f*w over (a,b), */
	/*                    i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qk15w */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);

	centr = (*a + *b) * .5f;
	hlgth = (*b - *a) * .5f;
	dhlgth = fabs(hlgth);

	/*           compute the 15-point kronrod approximation to the */
	/*           integral, and estimate the error. */

	fc = f(centr) * w(&centr, p1, p2, p3, p4, kp);
	resg = wg[3] * fc;
	resk = wgk[7] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 3; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		absc1 = centr - absc;
		absc2 = centr + absc;
		fval1 = f(absc1) * w(&absc1, p1, p2, p3, p4, kp);
		fval2 = f(absc2) * w(&absc2, p1, p2, p3, p4, kp);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 4; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		absc1 = centr - absc;
		absc2 = centr + absc;
		fval1 = f(absc1) * w(&absc1, p1, p2, p3, p4, kp);
		fval2 = f(absc2) * w(&absc2, p1, p2, p3, p4, kp);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5f;
	*resasc = wgk[7] * (r__1 = fc - reskh, fabs(r__1));
	for (j = 1; j <= 7; ++j) {
		*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, fabs(r__1)) + (
											 r__2 = fv2[j - 1] - reskh, fabs(r__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (r__1 = (resk - resg) * hlgth, fabs(r__1));
	if (*resasc != 0.f && *abserr != 0.f) {
		/* Computing MIN */
		d__1 = (double) (*abserr * 200.f / *resasc);
		r__1 = 1.f, r__2 = pow(d__1, c_b270);
		*abserr = *resasc * min(r__1,r__2);
	}
	if (*resabs > uflow / (epmach * 50.f)) {
		/* Computing MAX */
		r__1 = epmach * 50.f * *resabs;
		*abserr = max(r__1,*abserr);
	}
	return;
} /* qk15w_ */

void qk21_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result, 
			  sys_float *abserr, sys_float *resabs, sys_float *resasc)
{
	/* Initialized data */

	sys_float xgk[11] = { .9956571630258081f,.9739065285171717f,
									 .9301574913557082f,.8650633666889845f,.7808177265864169f,
									 .6794095682990244f,.5627571346686047f,.4333953941292472f,
									 .2943928627014602f,.1488743389816312f,0.f };
	sys_float wgk[11] = { .01169463886737187f,.03255816230796473f,
									 .054755896574352f,.07503967481091995f,.09312545458369761f,
									 .1093871588022976f,.1234919762620659f,.1347092173114733f,
									 .1427759385770601f,.1477391049013385f,.1494455540029169f };
	sys_float wg[5] = { .06667134430868814f,.1494513491505806f,
								  .219086362515982f,.2692667193099964f,.2955242247147529f };

	/* System generated locals */
	sys_float r__1, r__2;
	double d__1;

	/* Local variables */
	long j;
	sys_float fc, fv1[10], fv2[10];
	long jtw;
	sys_float absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	sys_float hlgth, centr, reskh, uflow;
	sys_float epmach, dhlgth;

	/* ***begin prologue  qk21 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  21-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b), with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           sys_float version */

	/*           parameters */
	/*            on entry */
	/*              f      - sys_float */
	/*                       function subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the driver program. */

	/*              a      - sys_float */
	/*                       lower limit of integration */

	/*              b      - sys_float */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - sys_float */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 21-polong */
	/*                       kronrod rule (resk) obtained by optimal addition */
	/*                       of abscissae to the 10-point gauss rule (resg). */

	/*              abserr - sys_float */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should not exceed fabs(i-result) */

	/*              resabs - sys_float */
	/*                       approximation to the integral j */

	/*              resasc - sys_float */
	/*                       approximation to the integral of fabs(f-i/(b-a)) */
	/*                       over (a,b) */

	/* ***references  (none) */
	/* ***routines called  r1mach */
	/* ***end prologue  qk21 */



	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 21-point kronrod rule */
	/*                    xgk(2), xgk(4), ...  abscissae of the 10-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 10-point gauss rule */

	/*           wgk    - weights of the 21-point kronrod rule */

	/*           wg     - weights of the 10-point gauss rule */





	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 10-point gauss formula */
	/*           resk   - result of the 21-point kronrod formula */
	/*           reskh  - approximation to the mean value of f over (a,b), */
	/*                    i.e. to i/(b-a) */


	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qk21 */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);

	centr = (*a + *b) * .5f;
	hlgth = (*b - *a) * .5f;
	dhlgth = fabs(hlgth);

	/*           compute the 21-point kronrod approximation to */
	/*           the integral, and estimate the absolute error. */

	resg = 0.f;
	fc = f(centr);
	resk = wgk[10] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 5; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 5; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5f;
	*resasc = wgk[10] * (r__1 = fc - reskh, fabs(r__1));
	for (j = 1; j <= 10; ++j) {
		*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, fabs(r__1)) + (
											 r__2 = fv2[j - 1] - reskh, fabs(r__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (r__1 = (resk - resg) * hlgth, fabs(r__1));
	if (*resasc != 0.f && *abserr != 0.f) {
		/* Computing MIN */
		d__1 = (double) (*abserr * 200.f / *resasc);
		r__1 = 1.f, r__2 = pow(d__1, c_b270);
		*abserr = *resasc * min(r__1,r__2);
	}
	if (*resabs > uflow / (epmach * 50.f)) {
		/* Computing MAX */
		r__1 = epmach * 50.f * *resabs;
		*abserr = max(r__1,*abserr);
	}
	return;
} /* qk21_ */

void qk31_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result,
			  sys_float *abserr, sys_float *resabs, sys_float *resasc)
{
	/* Initialized data */

	sys_float xgk[16] = { .9980022986933971f,.9879925180204854f,
									 .9677390756791391f,.9372733924007059f,.8972645323440819f,
									 .8482065834104272f,.7904185014424659f,.72441773136017f,
									 .650996741297417f,.5709721726085388f,.4850818636402397f,
									 .3941513470775634f,.2991800071531688f,.2011940939974345f,
									 .1011420669187175f,0.f };
	sys_float wgk[16] = { .005377479872923349f,.01500794732931612f,
									 .02546084732671532f,.03534636079137585f,.04458975132476488f,
									 .05348152469092809f,.06200956780067064f,.06985412131872826f,
									 .07684968075772038f,.08308050282313302f,.08856444305621177f,
									 .09312659817082532f,.09664272698362368f,.09917359872179196f,
									 .1007698455238756f,.1013300070147915f };
	sys_float wg[8] = { .03075324199611727f,.07036604748810812f,
								  .1071592204671719f,.1395706779261543f,.1662692058169939f,
								  .1861610000155622f,.1984314853271116f,.2025782419255613f };

	/* System generated locals */
	sys_float r__1, r__2;
	double d__1;

	/* Local variables */
	long j;
	sys_float fc, fv1[15], fv2[15];
	long jtw;
	sys_float absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	sys_float hlgth, centr, reskh, uflow;
	sys_float epmach, dhlgth;

	/* ***begin prologue  qk31 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  31-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b) with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           sys_float version */

	/*           parameters */
	/*            on entry */
	/*              f      - sys_float */
	/*                       function subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the calling program. */

	/*              a      - sys_float */
	/*                       lower limit of integration */

	/*              b      - sys_float */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - sys_float */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 31-polong */
	/*                       gauss-kronrod rule (resk), obtained by optimal */
	/*                       addition of abscissae to the 15-point gauss */
	/*                       rule (resg). */

	/*              abserr - sys_float */
	/*                       estimate of the modulus of the modulus, */
	/*                       which should not exceed fabs(i-result) */

	/*              resabs - sys_float */
	/*                       approximation to the integral j */

	/*              resasc - sys_float */
	/*                       approximation to the integral of fabs(f-i/(b-a)) */
	/*                       over (a,b) */

	/* ***references  (none) */
	/* ***routines called  r1mach */
	/* ***end prologue  qk31 */


	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 31-point kronrod rule */
	/*                    xgk(2), xgk(4), ...  abscissae of the 15-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 15-point gauss rule */

	/*           wgk    - weights of the 31-point kronrod rule */

	/*           wg     - weights of the 15-point gauss rule */



	/*           list of major variables */
	/*           ----------------------- */
	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 15-point gauss formula */
	/*           resk   - result of the 31-point kronrod formula */
	/*           reskh  - approximation to the mean value of f over (a,b), */
	/*                    i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */
	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qk31 */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);

	centr = (*a + *b) * .5f;
	hlgth = (*b - *a) * .5f;
	dhlgth = fabs(hlgth);

	/*           compute the 31-point kronrod approximation to */
	/*           the integral, and estimate the absolute error. */

	fc = f(centr);
	resg = wg[7] * fc;
	resk = wgk[15] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 7; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 8; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5f;
	*resasc = wgk[15] * (r__1 = fc - reskh, fabs(r__1));
	for (j = 1; j <= 15; ++j) {
		*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, fabs(r__1)) + (
											 r__2 = fv2[j - 1] - reskh, fabs(r__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (r__1 = (resk - resg) * hlgth, fabs(r__1));
	if (*resasc != 0.f && *abserr != 0.f) {
		/* Computing MIN */
		d__1 = (double) (*abserr * 200.f / *resasc);
		r__1 = 1.f, r__2 = pow(d__1, c_b270);
		*abserr = *resasc * min(r__1,r__2);
	}
	if (*resabs > uflow / (epmach * 50.f)) {
		/* Computing MAX */
		r__1 = epmach * 50.f * *resabs;
		*abserr = max(r__1,*abserr);
	}
	return;
} /* qk31_ */

void qk41_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result,
			  sys_float *abserr, sys_float *resabs, sys_float *resasc)
{
	/* Initialized data */

	sys_float xgk[21] = { .9988590315882777f,.9931285991850949f,
									 .9815078774502503f,.9639719272779138f,.9408226338317548f,
									 .9122344282513259f,.878276811252282f,.8391169718222188f,
									 .7950414288375512f,.7463319064601508f,.6932376563347514f,
									 .636053680726515f,.5751404468197103f,.5108670019508271f,
									 .4435931752387251f,.3737060887154196f,.301627868114913f,
									 .2277858511416451f,.1526054652409227f,.07652652113349733f,0.f };
	sys_float wgk[21] = { .003073583718520532f,.008600269855642942f,
									 .01462616925697125f,.02038837346126652f,.02588213360495116f,
									 .0312873067770328f,.0366001697582008f,.04166887332797369f,
									 .04643482186749767f,.05094457392372869f,.05519510534828599f,
									 .05911140088063957f,.06265323755478117f,.06583459713361842f,
									 .06864867292852162f,.07105442355344407f,.07303069033278667f,
									 .07458287540049919f,.07570449768455667f,.07637786767208074f,
									 .07660071191799966f };
	sys_float wg[10] = { .01761400713915212f,.04060142980038694f,
									.06267204833410906f,.08327674157670475f,.1019301198172404f,
									.1181945319615184f,.1316886384491766f,.1420961093183821f,
									.1491729864726037f,.1527533871307259f };

	/* System generated locals */
	sys_float r__1, r__2;
	double d__1;

	/* Local variables */
	long j;
	sys_float fc, fv1[20], fv2[20];
	long jtw;
	sys_float absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	sys_float hlgth, centr, reskh, uflow;
	sys_float epmach, dhlgth;

	/* ***begin prologue  qk41 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  41-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b), with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           sys_float version */

	/*           parameters */
	/*            on entry */
	/*              f      - sys_float */
	/*                       function subprogram defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the calling program. */

	/*              a      - sys_float */
	/*                       lower limit of integration */

	/*              b      - sys_float */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - sys_float */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 41-polong */
	/*                       gauss-kronrod rule (resk) obtained by optimal */
	/*                       addition of abscissae to the 20-point gauss */
	/*                       rule (resg). */

	/*              abserr - sys_float */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should not exceed fabs(i-result) */

	/*              resabs - sys_float */
	/*                       approximation to the integral j */

	/*              resasc - sys_float */
	/*                       approximation to the integal of fabs(f-i/(b-a)) */
	/*                       over (a,b) */

	/* ***references  (none) */
	/* ***routines called  r1mach */
	/* ***end prologue  qk41 */



	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 41-point gauss-kronrod rule */
	/*                    xgk(2), xgk(4), ...  abscissae of the 20-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 20-point gauss rule */

	/*           wgk    - weights of the 41-point gauss-kronrod rule */

	/*           wg     - weights of the 20-point gauss rule */



	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 20-point gauss formula */
	/*           resk   - result of the 41-point kronrod formula */
	/*           reskh  - approximation to mean value of f over (a,b), i.e. */
	/*                    to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qk41 */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);

	centr = (*a + *b) * .5f;
	hlgth = (*b - *a) * .5f;
	dhlgth = fabs(hlgth);

	/*           compute the 41-point gauss-kronrod approximation to */
	/*           the integral, and estimate the absolute error. */

	resg = 0.f;
	fc = f(centr);
	resk = wgk[20] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 10; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 10; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5f;
	*resasc = wgk[20] * (r__1 = fc - reskh, fabs(r__1));
	for (j = 1; j <= 20; ++j) {
		*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, fabs(r__1)) + (
											 r__2 = fv2[j - 1] - reskh, fabs(r__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (r__1 = (resk - resg) * hlgth, fabs(r__1));
	if (*resasc != 0.f && *abserr != 0.f) {
		/* Computing MIN */
		d__1 = (double) (*abserr * 200.f / *resasc);
		r__1 = 1.f, r__2 = pow(d__1, c_b270);
		*abserr = *resasc * min(r__1,r__2);
	}
	if (*resabs > uflow / (epmach * 50.f)) {
		/* Computing MAX */
		r__1 = epmach * 50.f * *resabs;
		*abserr = max(r__1,*abserr);
	}
	return;
} /* qk41_ */

void qk51_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result, 
			  sys_float *abserr, sys_float *resabs, sys_float *resasc)
{
	/* Initialized data */

	sys_float xgk[26] = { .9992621049926098f,.9955569697904981f,
									 .9880357945340772f,.9766639214595175f,.9616149864258425f,
									 .9429745712289743f,.9207471152817016f,.8949919978782754f,
									 .8658470652932756f,.833442628760834f,.7978737979985001f,
									 .7592592630373576f,.7177664068130844f,.6735663684734684f,
									 .6268100990103174f,.577662930241223f,.5263252843347192f,
									 .473002731445715f,.4178853821930377f,.3611723058093878f,
									 .3030895389311078f,.2438668837209884f,.1837189394210489f,
									 .1228646926107104f,.06154448300568508f,0.f };
	sys_float wgk[26] = { .001987383892330316f,.005561932135356714f,
									 .009473973386174152f,.01323622919557167f,.0168478177091283f,
									 .02043537114588284f,.02400994560695322f,.02747531758785174f,
									 .03079230016738749f,.03400213027432934f,.03711627148341554f,
									 .04008382550403238f,.04287284502017005f,.04550291304992179f,
									 .04798253713883671f,.05027767908071567f,.05236288580640748f,
									 .05425112988854549f,.05595081122041232f,.05743711636156783f,
									 .05868968002239421f,.05972034032417406f,.06053945537604586f,
									 .06112850971705305f,.06147118987142532f,.06158081806783294f };
	sys_float wg[13] = { .01139379850102629f,.02635498661503214f,
									.04093915670130631f,.05490469597583519f,.06803833381235692f,
									.08014070033500102f,.09102826198296365f,.1005359490670506f,
									.1085196244742637f,.1148582591457116f,.1194557635357848f,
									.12224244299031f,.1231760537267155f };

	/* System generated locals */
	sys_float r__1, r__2;
	double d__1;

	/* Local variables */
	long j;
	sys_float fc, fv1[25], fv2[25];
	long jtw;
	sys_float absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	sys_float hlgth, centr, reskh, uflow;
	sys_float epmach, dhlgth;

	/* ***begin prologue  qk51 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  51-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b) with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*           integration rules */
	/*           standard fortran subroutine */
	/*           sys_float version */

	/*           parameters */
	/*            on entry */
	/*              f      - sys_float */
	/*                       function subroutine defining the integrand */
	/*                       function f(x). the actual name for f needs to be */
	/*                       declared e x t e r n a l in the calling program. */

	/*              a      - sys_float */
	/*                       lower limit of integration */

	/*              b      - sys_float */
	/*                       upper limit of integration */

	/*            on return */
	/*              result - sys_float */
	/*                       approximation to the integral i */
	/*                       result is computed by applying the 51-polong */
	/*                       kronrod rule (resk) obtained by optimal addition */
	/*                       of abscissae to the 25-point gauss rule (resg). */

	/*              abserr - sys_float */
	/*                       estimate of the modulus of the absolute error, */
	/*                       which should not exceed fabs(i-result) */

	/*              resabs - sys_float */
	/*                       approximation to the integral j */

	/*              resasc - sys_float */
	/*                       approximation to the integral of fabs(f-i/(b-a)) */
	/*                       over (a,b) */

	/* ***references  (none) */
	/* ***routines called  r1mach */
	/* ***end prologue  qk51 */



	/*           the abscissae and weights are given for the interval (-1,1). */
	/*           because of symmetry only the positive abscissae and their */
	/*           corresponding weights are given. */

	/*           xgk    - abscissae of the 51-point kronrod rule */
	/*                    xgk(2), xgk(4), ...  abscissae of the 25-polong */
	/*                    gauss rule */
	/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
	/*                    added to the 25-point gauss rule */

	/*           wgk    - weights of the 51-point kronrod rule */

	/*           wg     - weights of the 25-point gauss rule */



	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 25-point gauss formula */
	/*           resk   - result of the 51-point kronrod formula */
	/*           reskh  - approximation to the mean value of f over (a,b), */
	/*                    i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qk51 */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);

	centr = (*a + *b) * .5f;
	hlgth = (*b - *a) * .5f;
	dhlgth = fabs(hlgth);

	/*           compute the 51-point kronrod approximation to */
	/*           the integral, and estimate the absolute error. */

	fc = f(centr);
	resg = wg[12] * fc;
	resk = wgk[25] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 12; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 13; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5f;
	*resasc = wgk[25] * (r__1 = fc - reskh, fabs(r__1));
	for (j = 1; j <= 25; ++j) {
		*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, fabs(r__1)) + (
											 r__2 = fv2[j - 1] - reskh, fabs(r__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (r__1 = (resk - resg) * hlgth, fabs(r__1));
	if (*resasc != 0.f && *abserr != 0.f) {
		/* Computing MIN */
		d__1 = (double) (*abserr * 200.f / *resasc);
		r__1 = 1.f, r__2 = pow(d__1, c_b270);
		*abserr = *resasc * min(r__1,r__2);
	}
	if (*resabs > uflow / (epmach * 50.f)) {
		/* Computing MAX */
		r__1 = epmach * 50.f * *resabs;
		*abserr = max(r__1,*abserr);
	}
	return;
} /* qk51_ */

void qk61_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result, 
			  sys_float *abserr, sys_float *resabs, sys_float *resasc)
{
	/* Initialized data */

	sys_float xgk[31] = { .9994844100504906f,.9968934840746495f,
									 .9916309968704046f,.9836681232797472f,.9731163225011263f,
									 .9600218649683075f,.94437444474856f,.9262000474292743f,
									 .9055733076999078f,.8825605357920527f,.8572052335460611f,
									 .8295657623827684f,.7997278358218391f,.7677774321048262f,
									 .7337900624532268f,.6978504947933158f,.660061064126627f,
									 .6205261829892429f,.5793452358263617f,.5366241481420199f,
									 .4924804678617786f,.4470337695380892f,.4004012548303944f,
									 .3527047255308781f,.3040732022736251f,.2546369261678898f,
									 .2045251166823099f,.1538699136085835f,.102806937966737f,
									 .0514718425553177f,0.f };
	sys_float wgk[31] = { .001389013698677008f,.003890461127099884f,
									 .006630703915931292f,.009273279659517763f,.01182301525349634f,
									 .0143697295070458f,.01692088918905327f,.01941414119394238f,
									 .02182803582160919f,.0241911620780806f,.0265099548823331f,
									 .02875404876504129f,.03090725756238776f,.03298144705748373f,
									 .03497933802806002f,.03688236465182123f,.03867894562472759f,
									 .04037453895153596f,.04196981021516425f,.04345253970135607f,
									 .04481480013316266f,.04605923827100699f,.04718554656929915f,
									 .04818586175708713f,.04905543455502978f,.04979568342707421f,
									 .05040592140278235f,.05088179589874961f,.05122154784925877f,
									 .05142612853745903f,.05149472942945157f };
	sys_float wg[15] = { .007968192496166606f,.01846646831109096f,
									.02878470788332337f,.03879919256962705f,.04840267283059405f,
									.05749315621761907f,.0659742298821805f,.07375597473770521f,
									.08075589522942022f,.08689978720108298f,.09212252223778613f,
									.09636873717464426f,.09959342058679527f,.1017623897484055f,
									.1028526528935588f };

	/* System generated locals */
	sys_float r__1, r__2;
	double d__1;

	/* Local variables */
	long j;
	sys_float fc, fv1[30], fv2[30];
	long jtw;
	sys_float absc, resg, resk, fsum, fval1, fval2;
	long jtwm1;
	sys_float hlgth, centr, reskh, uflow;
	sys_float epmach, dhlgth;

	/* ***begin prologue  qk61 */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a2 */
	/* ***keywords  61-point gauss-kronrod rules */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  to compute i = integral of f over (a,b) with error */
	/*                           estimate */
	/*                       j = integral of fabs(f) over (a,b) */
	/* ***description */

	/*        integration rule */
	/*        standard fortran subroutine */
	/*        sys_float version */


	/*        parameters */
	/*         on entry */
	/*           f      - sys_float */
	/*                    function subprogram defining the integrand */
	/*                    function f(x). the actual name for f needs to be */
	/*                    declared e x t e r n a l in the calling program. */

	/*           a      - sys_float */
	/*                    lower limit of integration */

	/*           b      - sys_float */
	/*                    upper limit of integration */

	/*         on return */
	/*           result - sys_float */
	/*                    approximation to the integral i */
	/*                    result is computed by applying the 61-polong */
	/*                    kronrod rule (resk) obtained by optimal addition of */
	/*                    abscissae to the 30-point gauss rule (resg). */

	/*           abserr - sys_float */
	/*                    estimate of the modulus of the absolute error, */
	/*                    which should equal or exceed fabs(i-result) */

	/*           resabs - sys_float */
	/*                    approximation to the integral j */

	/*           resasc - sys_float */
	/*                    approximation to the integral of fabs(f-i/(b-a)) */


	/* ***references  (none) */
	/* ***routines called  r1mach */
	/* ***end prologue  qk61 */



	/*           the abscissae and weights are given for the */
	/*           interval (-1,1). because of symmetry only the positive */
	/*           abscissae and their corresponding weights are given. */

	/*           xgk   - abscissae of the 61-point kronrod rule */
	/*                   xgk(2), xgk(4)  ... abscissae of the 30-polong */
	/*                   gauss rule */
	/*                   xgk(1), xgk(3)  ... optimally added abscissae */
	/*                   to the 30-point gauss rule */

	/*           wgk   - weights of the 61-point kronrod rule */

	/*           wg    - weigths of the 30-point gauss rule */


	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the interval */
	/*           hlgth  - half-length of the interval */
	/*           absc   - abscissa */
	/*           fval*  - function value */
	/*           resg   - result of the 30-point gauss rule */
	/*           resk   - result of the 61-point kronrod rule */
	/*           reskh  - approximation to the mean value of f */
	/*                    over (a,b), i.e. to i/(b-a) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qk61 */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);

	centr = (*b + *a) * .5f;
	hlgth = (*b - *a) * .5f;
	dhlgth = fabs(hlgth);

	/*           compute the 61-point kronrod approximation to the */
	/*           integral, and estimate the absolute error. */

	resg = 0.f;
	fc = f(centr);
	resk = wgk[30] * fc;
	*resabs = fabs(resk);
	for (j = 1; j <= 15; ++j) {
		jtw = j << 1;
		absc = hlgth * xgk[jtw - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtw - 1] = fval1;
		fv2[jtw - 1] = fval2;
		fsum = fval1 + fval2;
		resg += wg[j - 1] * fsum;
		resk += wgk[jtw - 1] * fsum;
		*resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
		/* L10: */
	}
	for (j = 1; j <= 15; ++j) {
		jtwm1 = (j << 1) - 1;
		absc = hlgth * xgk[jtwm1 - 1];
		r__1 = centr - absc;
		fval1 = f(r__1);
		r__1 = centr + absc;
		fval2 = f(r__1);
		fv1[jtwm1 - 1] = fval1;
		fv2[jtwm1 - 1] = fval2;
		fsum = fval1 + fval2;
		resk += wgk[jtwm1 - 1] * fsum;
		*resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
		/* L15: */
	}
	reskh = resk * .5f;
	*resasc = wgk[30] * (r__1 = fc - reskh, fabs(r__1));
	for (j = 1; j <= 30; ++j) {
		*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, fabs(r__1)) + (
											 r__2 = fv2[j - 1] - reskh, fabs(r__2)));
		/* L20: */
	}
	*result = resk * hlgth;
	*resabs *= dhlgth;
	*resasc *= dhlgth;
	*abserr = (r__1 = (resk - resg) * hlgth, fabs(r__1));
	if (*resasc != 0.f && *abserr != 0.f) {
		/* Computing MIN */
		d__1 = (double) (*abserr * 200.f / *resasc);
		r__1 = 1.f, r__2 = pow(d__1, c_b270);
		*abserr = *resasc * min(r__1,r__2);
	}
	if (*resabs > uflow / (epmach * 50.f)) {
		/* Computing MAX */
		r__1 = epmach * 50.f * *resabs;
		*abserr = max(r__1,*abserr);
	}
	return;
} /* qk61_ */

void qmomo_(const sys_float *alfa, const sys_float *beta, sys_float *ri, 
				sys_float *rj, sys_float *rg, sys_float *rh, const long *integr)
{
	/* System generated locals */
	double d__1;

	/* Local variables */
	long i__;
	sys_float an;
	long im1;
	sys_float anm1, ralf, rbet, alfp1, alfp2, betp1, betp2;

	/* ***begin prologue  qmomo */
	/* ***date written   810101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a2a1,c3a2 */
	/* ***keywords  modified chebyshev moments */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  this routine computes modified chebsyshev moments. the k-th */
	/*            modified chebyshev moment is defined as the integral over */
	/*            (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev */
	/*            polynomial of degree k. */
	/* ***description */

	/*        modified chebyshev moments */
	/*        standard fortran subroutine */
	/*        sys_float version */

	/*        parameters */
	/*           alfa   - sys_float */
	/*                    parameter in the weight function w(x), alfa.gt.(-1) */

	/*           beta   - sys_float */
	/*                    parameter in the weight function w(x), beta.gt.(-1) */

	/*           ri     - sys_float */
	/*                    vector of dimension 25 */
	/*                    ri(k) is the integral over (-1,1) of */
	/*                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25. */

	/*           rj     - sys_float */
	/*                    vector of dimension 25 */
	/*                    rj(k) is the integral over (-1,1) of */
	/*                    (1-x)**beta*t(k-1,x), k = 1, ..., 25. */

	/*           rg     - sys_float */
	/*                    vector of dimension 25 */
	/*                    rg(k) is the integral over (-1,1) of */
	/*                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25. */

	/*           rh     - sys_float */
	/*                    vector of dimension 25 */
	/*                    rh(k) is the integral over (-1,1) of */
	/*                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25. */

	/*           integr - long */
	/*                    input parameter indicating the modified */
	/*                    moments to be computed */
	/*                    integr = 1 compute ri, rj */
	/*                           = 2 compute ri, rj, rg */
	/*                           = 3 compute ri, rj, rh */
	/*                           = 4 compute ri, rj, rg, rh */
	/* ***references  (none) */
	/* ***routines called  (none) */
	/* ***end prologue  qmomo */




	/* ***first executable statement  qmomo */
	/* Parameter adjustments */
	--rh;
	--rg;
	--rj;
	--ri;

	/* Function Body */
	alfp1 = *alfa + 1.f;
	betp1 = *beta + 1.f;
	alfp2 = *alfa + 2.f;
	betp2 = *beta + 2.f;
	d__1 = (double) alfp1;
	ralf = pow(c_b320, d__1);
	d__1 = (double) betp1;
	rbet = pow(c_b320, d__1);

	/*           compute ri, rj using a forward recurrence relation. */

	ri[1] = ralf / alfp1;
	rj[1] = rbet / betp1;
	ri[2] = ri[1] * *alfa / alfp2;
	rj[2] = rj[1] * *beta / betp2;
	an = 2.f;
	anm1 = 1.f;
	for (i__ = 3; i__ <= 25; ++i__) {
		ri[i__] = -(ralf + an * (an - alfp2) * ri[i__ - 1]) / (anm1 * (an + 
																							alfp1));
		rj[i__] = -(rbet + an * (an - betp2) * rj[i__ - 1]) / (anm1 * (an + 
																							betp1));
		anm1 = an;
		an += 1.f;
		/* L20: */
	}
	if (*integr == 1) {
		goto L70;
	}
	if (*integr == 3) {
		goto L40;
	}

	/*           compute rg using a forward recurrence relation. */

	rg[1] = -ri[1] / alfp1;
	rg[2] = -(ralf + ralf) / (alfp2 * alfp2) - rg[1];
	an = 2.f;
	anm1 = 1.f;
	im1 = 2;
	for (i__ = 3; i__ <= 25; ++i__) {
		rg[i__] = -(an * (an - alfp2) * rg[im1] - an * ri[im1] + anm1 * ri[
							i__]) / (anm1 * (an + alfp1));
		anm1 = an;
		an += 1.f;
		im1 = i__;
		/* L30: */
	}
	if (*integr == 2) {
		goto L70;
	}

	/*           compute rh using a forward recurrence relation. */

L40:
	rh[1] = -rj[1] / betp1;
	rh[2] = -(rbet + rbet) / (betp2 * betp2) - rh[1];
	an = 2.f;
	anm1 = 1.f;
	im1 = 2;
	for (i__ = 3; i__ <= 25; ++i__) {
		rh[i__] = -(an * (an - betp2) * rh[im1] - an * rj[im1] + anm1 * rj[
							i__]) / (anm1 * (an + betp1));
		anm1 = an;
		an += 1.f;
		im1 = i__;
		/* L50: */
	}
	for (i__ = 2; i__ <= 25; i__ += 2) {
		rh[i__] = -rh[i__];
		/* L60: */
	}
L70:
	for (i__ = 2; i__ <= 25; i__ += 2) {
		rj[i__] = -rj[i__];
		/* L80: */
	}
	/* L90: */
	return;
} /* qmomo_ */

void qng_(const E_fp& f, sys_float *a, sys_float *b, sys_float *epsabs, sys_float *
			 epsrel, sys_float *result, sys_float *abserr, long *neval, long *ier)
{
	/* Initialized data */

	sys_float x1[5] = { .9739065285171717f,.8650633666889845f,
								  .6794095682990244f,.4333953941292472f,.1488743389816312f };
	sys_float x2[5] = { .9956571630258081f,.9301574913557082f,
								  .7808177265864169f,.5627571346686047f,.2943928627014602f };
	sys_float x3[11] = { .9993333609019321f,.9874334029080889f,
									.9548079348142663f,.9001486957483283f,.8251983149831142f,
									.732148388989305f,.6228479705377252f,.4994795740710565f,
									.3649016613465808f,.2222549197766013f,.07465061746138332f };
	sys_float x4[22] = { .9999029772627292f,.9979898959866787f,
									.9921754978606872f,.9813581635727128f,.9650576238583846f,
									.9431676131336706f,.9158064146855072f,.8832216577713165f,
									.8457107484624157f,.803557658035231f,.7570057306854956f,
									.7062732097873218f,.6515894665011779f,.5932233740579611f,
									.5314936059708319f,.4667636230420228f,.3994248478592188f,
									.3298748771061883f,.2585035592021616f,.1856953965683467f,
									.1118422131799075f,.03735212339461987f };
	sys_float w10[5] = { .06667134430868814f,.1494513491505806f,
									.219086362515982f,.2692667193099964f,.2955242247147529f };
	sys_float w21a[5] = { .03255816230796473f,.07503967481091995f,
									 .1093871588022976f,.1347092173114733f,.1477391049013385f };
	sys_float w21b[6] = { .01169463886737187f,.054755896574352f,
									 .09312545458369761f,.1234919762620659f,.1427759385770601f,
									 .1494455540029169f };
	sys_float w43a[10] = { .01629673428966656f,.0375228761208695f,
									  .05469490205825544f,.06735541460947809f,.07387019963239395f,
									  .005768556059769796f,.02737189059324884f,.04656082691042883f,
									  .06174499520144256f,.0713872672686934f };
	sys_float w43b[12] = { .001844477640212414f,.01079868958589165f,
									  .02189536386779543f,.03259746397534569f,.04216313793519181f,
									  .05074193960018458f,.05837939554261925f,.06474640495144589f,
									  .06956619791235648f,.07282444147183321f,.07450775101417512f,
									  .07472214751740301f };
	sys_float w87a[21] = { .008148377384149173f,.01876143820156282f,
									  .02734745105005229f,.03367770731163793f,.03693509982042791f,
									  .002884872430211531f,.0136859460227127f,.02328041350288831f,
									  .03087249761171336f,.03569363363941877f,9.152833452022414e-4f,
									  .005399280219300471f,.01094767960111893f,.01629873169678734f,
									  .02108156888920384f,.02537096976925383f,.02918969775647575f,
									  .03237320246720279f,.03478309895036514f,.03641222073135179f,
									  .03725387550304771f };
	sys_float w87b[23] = { 2.741455637620724e-4f,.001807124155057943f,
									  .004096869282759165f,.006758290051847379f,.009549957672201647f,
									  .01232944765224485f,.01501044734638895f,.01754896798624319f,
									  .01993803778644089f,.02219493596101229f,.02433914712600081f,
									  .02637450541483921f,.0282869107887712f,.0300525811280927f,
									  .03164675137143993f,.0330504134199785f,.03425509970422606f,
									  .03526241266015668f,.0360769896228887f,.03669860449845609f,
									  .03712054926983258f,.03733422875193504f,.03736107376267902f };

	/* System generated locals */
	sys_float r__1, r__2, r__3, r__4;
	double d__1;

	/* Local variables */
	long k, l;
	sys_float fv1[5], fv2[5], fv3[5], fv4[5];
	long ipx=0;
	sys_float absc, fval, res10, res21=0., res43=0., res87, fval1, fval2, hlgth, 
		centr, reskh, uflow;
	sys_float epmach, dhlgth, resabs=0., resasc=0., fcentr, savfun[21];

	/* ***begin prologue  qng */
	/* ***date written   800101   (yymmdd) */
	/* ***revision date  830518   (yymmdd) */
	/* ***category no.  h2a1a1 */
	/* ***keywords  automatic integrator, smooth integrand, */
	/*             non-adaptive, gauss-kronrod(patterson) */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl math & progr. div. - k.u.leuven */
	/*           kahaner,david,nbs - modified (2/82) */
	/* ***purpose  the routine calculates an approximation result to a */
	/*            given definite integral i = integral of f over (a,b), */
	/*            hopefully satisfying following claim for accuracy */
	/*            fabs(i-result).le.max(epsabs,epsrel*fabs(i)). */
	/* ***description */

	/* non-adaptive integration */
	/* standard fortran subroutine */
	/* sys_float version */

	/*           f      - sys_float version */
	/*                    function subprogram defining the integrand function */
	/*                    f(x). the actual name for f needs to be declared */
	/*                    e x t e r n a l in the driver program. */

	/*           a      - sys_float version */
	/*                    lower limit of integration */

	/*           b      - sys_float version */
	/*                    upper limit of integration */

	/*           epsabs - sys_float */
	/*                    absolute accuracy requested */
	/*           epsrel - sys_float */
	/*                    relative accuracy requested */
	/*                    if  epsabs.le.0 */
	/*                    and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
	/*                    the routine will end with ier = 6. */

	/*         on return */
	/*           result - sys_float */
	/*                    approximation to the integral i */
	/*                    result is obtained by applying the 21-polong */
	/*                    gauss-kronrod rule (res21) obtained by optimal */
	/*                    addition of abscissae to the 10-point gauss rule */
	/*                    (res10), or by applying the 43-point rule (res43) */
	/*                    obtained by optimal addition of abscissae to the */
	/*                    21-point gauss-kronrod rule, or by applying the */
	/*                    87-point rule (res87) obtained by optimal addition */
	/*                    of abscissae to the 43-point rule. */

	/*           abserr - sys_float */
	/*                    estimate of the modulus of the absolute error, */
	/*                    which should equal or exceed fabs(i-result) */

	/*           neval  - long */
	/*                    number of integrand evaluations */

	/*           ier    - ier = 0 normal and reliable termination of the */
	/*                            routine. it is assumed that the requested */
	/*                            accuracy has been achieved. */
	/*                    ier.gt.0 abnormal termination of the routine. it is */
	/*                            assumed that the requested accuracy has */
	/*                            not been achieved. */
	/*           error messages */
	/*                    ier = 1 the maximum number of steps has been */
	/*                            executed. the integral is probably too */
	/*                            difficult to be calculated by dqng. */
	/*                        = 6 the input is invalid, because */
	/*                            epsabs.le.0 and */
	/*                            epsrel.lt.max(50*rel.mach.acc.,0.5d-28). */
	/*                            result, abserr and neval are set to zero. */

	/* ***references  (none) */
	/* ***routines called  r1mach,xerror */
	/* ***end prologue  qng */



	/*           the following data statements contain the */
	/*           abscissae and weights of the integration rules used. */

	/*           x1      abscissae common to the 10-, 21-, 43- */
	/*                   and 87-point rule */
	/*           x2      abscissae common to the 21-, 43- and */
	/*                   87-point rule */
	/*           x3      abscissae common to the 43- and 87-polong */
	/*                   rule */
	/*           x4      abscissae of the 87-point rule */
	/*           w10     weights of the 10-point formula */
	/*           w21a    weights of the 21-point formula for */
	/*                   abscissae x1 */
	/*           w21b    weights of the 21-point formula for */
	/*                   abscissae x2 */
	/*           w43a    weights of the 43-point formula for */
	/*                   abscissae x1, x3 */
	/*           w43b    weights of the 43-point formula for */
	/*                   abscissae x3 */
	/*           w87a    weights of the 87-point formula for */
	/*                   abscissae x1, x2, x3 */
	/*           w87b    weights of the 87-point formula for */
	/*                   abscissae x4 */


	/*           list of major variables */
	/*           ----------------------- */

	/*           centr  - mid point of the integration interval */
	/*           hlgth  - half-length of the integration interval */
	/*           fcentr - function value at mid polong */
	/*           absc   - abscissa */
	/*           fval   - function value */
	/*           savfun - array of function values which */
	/*                    have already been computed */
	/*           res10  - 10-point gauss result */
	/*           res21  - 21-point kronrod result */
	/*           res43  - 43-point result */
	/*           res87  - 87-point result */
	/*           resabs - approximation to the integral of fabs(f) */
	/*           resasc - approximation to the integral of fabs(f-i/(b-a)) */

	/*           machine dependent constants */
	/*           --------------------------- */

	/*           epmach is the largest relative spacing. */
	/*           uflow is the smallest positive magnitude. */

	/* ***first executable statement  qng */
	epmach = r1mach(c__4);
	uflow = r1mach(c__1);

	/*           test on validity of parameters */
	/*           ------------------------------ */

	*result = 0.f;
	*abserr = 0.f;
	*neval = 0;
	*ier = 6;
	/* Computing MAX */
	r__1 = 5e-15f, r__2 = epmach * 50.f;
	if (*epsabs <= 0.f && *epsrel < max(r__1,r__2)) {
		goto L80;
	}
	hlgth = (*b - *a) * .5f;
	dhlgth = fabs(hlgth);
	centr = (*b + *a) * .5f;
	fcentr = f(centr);
	*neval = 21;
	*ier = 1;

	/*          compute the integral using the 10- and 21-point formula. */

	for (l = 1; l <= 3; ++l) {
		switch (l) {
		case 1:  goto L5;
		case 2:  goto L25;
		case 3:  goto L45;
		}
	L5:
		res10 = 0.f;
		res21 = w21b[5] * fcentr;
		resabs = w21b[5] * fabs(fcentr);
		for (k = 1; k <= 5; ++k) {
			absc = hlgth * x1[k - 1];
			r__1 = centr + absc;
			fval1 = f(r__1);
			r__1 = centr - absc;
			fval2 = f(r__1);
			fval = fval1 + fval2;
			res10 += w10[k - 1] * fval;
			res21 += w21a[k - 1] * fval;
			resabs += w21a[k - 1] * (fabs(fval1) + fabs(fval2));
			savfun[k - 1] = fval;
			fv1[k - 1] = fval1;
			fv2[k - 1] = fval2;
			/* L10: */
		}
		ipx = 5;
		for (k = 1; k <= 5; ++k) {
			++ipx;
			absc = hlgth * x2[k - 1];
			r__1 = centr + absc;
			fval1 = f(r__1);
			r__1 = centr - absc;
			fval2 = f(r__1);
			fval = fval1 + fval2;
			res21 += w21b[k - 1] * fval;
			resabs += w21b[k - 1] * (fabs(fval1) + fabs(fval2));
			savfun[ipx - 1] = fval;
			fv3[k - 1] = fval1;
			fv4[k - 1] = fval2;
			/* L15: */
		}

		/*          test for convergence. */

		*result = res21 * hlgth;
		resabs *= dhlgth;
		reskh = res21 * .5f;
		resasc = w21b[5] * (r__1 = fcentr - reskh, fabs(r__1));
		for (k = 1; k <= 5; ++k) {
			resasc = resasc + w21a[k - 1] * ((r__1 = fv1[k - 1] - reskh, fabs(
															 r__1)) + (r__2 = fv2[k - 1] - reskh, fabs(r__2))) + w21b[
																 k - 1] * ((r__3 = fv3[k - 1] - reskh, fabs(r__3)) + (r__4 
																																		= fv4[k - 1] - reskh, fabs(r__4)));
			/* L20: */
		}
		*abserr = (r__1 = (res21 - res10) * hlgth, fabs(r__1));
		resasc *= dhlgth;
		goto L65;

		/*          compute the integral using the 43-point formula. */

	L25:
		res43 = w43b[11] * fcentr;
		*neval = 43;
		for (k = 1; k <= 10; ++k) {
			res43 += savfun[k - 1] * w43a[k - 1];
			/* L30: */
		}
		for (k = 1; k <= 11; ++k) {
			++ipx;
			absc = hlgth * x3[k - 1];
			r__1 = absc + centr;
			r__2 = centr - absc;
			fval = f(r__1) + f(r__2);
			res43 += fval * w43b[k - 1];
			savfun[ipx - 1] = fval;
			/* L40: */
		}

		/*          test for convergence. */

		*result = res43 * hlgth;
		*abserr = (r__1 = (res43 - res21) * hlgth, fabs(r__1));
		goto L65;

		/*          compute the integral using the 87-point formula. */

	L45:
		res87 = w87b[22] * fcentr;
		*neval = 87;
		for (k = 1; k <= 21; ++k) {
			res87 += savfun[k - 1] * w87a[k - 1];
			/* L50: */
		}
		for (k = 1; k <= 22; ++k) {
			absc = hlgth * x4[k - 1];
			r__1 = absc + centr;
			r__2 = centr - absc;
			res87 += w87b[k - 1] * (f(r__1) + f(r__2));
			/* L60: */
		}
		*result = res87 * hlgth;
		*abserr = (r__1 = (res87 - res43) * hlgth, fabs(r__1));
	L65:
		if (resasc != 0.f && *abserr != 0.f) {
			/* Computing MIN */
			d__1 = (double) (*abserr * 200.f / resasc);
			r__1 = 1.f, r__2 = pow(d__1, c_b270);
			*abserr = resasc * min(r__1,r__2);
		}
		if (resabs > uflow / (epmach * 50.f)) {
			/* Computing MAX */
			r__1 = epmach * 50.f * resabs;
			*abserr = max(r__1,*abserr);
		}
		/* Computing MAX */
		r__1 = *epsabs, r__2 = *epsrel * fabs(*result);
		if (*abserr <= max(r__1,r__2)) {
			*ier = 0;
		}
		/* ***jump out of do-loop */
		if (*ier == 0) {
			goto L999;
		}
		/* L70: */
	}
L80:
	xerror_("abnormal return from  qng ", &c__26, ier, &c__0, 26);
L999:
	return;
} /* qng_ */

void qpsrt_(const long *limit, long *last, long *maxerr, 
				sys_float *ermax, sys_float *elist, long *iord, long *nrmax)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	long i__, j, k, ido, ibeg, jbnd, isucc, jupbn;
	sys_float errmin, errmax;

	/* ***begin prologue  qpsrt */
	/* ***refer to  qage,qagie,qagpe,qagse,qawce,qawse,qawoe */
	/* ***routines called  (none) */
	/* ***keywords  sequential sorting */
	/* ***description */

	/* 1.        qpsrt */
	/*           ordering routine */
	/*              standard fortran subroutine */
	/*              sys_float version */

	/* 2.        purpose */
	/*              this routine maintains the descending ordering */
	/*              in the list of the local error estimates resulting from */
	/*              the interval subdivision process. at each call two error */
	/*              estimates are inserted using the sequential search */
	/*              method, top-down for the largest error estimate */
	/*              and bottom-up for the smallest error estimate. */

	/* 3.        calling sequence */
	/*              call qpsrt(limit,last,maxerr,ermax,elist,iord,nrmax) */

	/*           parameters (meaning at output) */
	/*              limit  - long */
	/*                       maximum number of error estimates the list */
	/*                       can contain */

	/*              last   - long */
	/*                       number of error estimates currently */
	/*                       in the list */

	/*              maxerr - long */
	/*                       maxerr points to the nrmax-th largest error */
	/*                       estimate currently in the list */

	/*              ermax  - sys_float */
	/*                       nrmax-th largest error estimate */
	/*                       ermax = elist(maxerr) */

	/*              elist  - sys_float */
	/*                       vector of dimension last containing */
	/*                       the error estimates */

	/*              iord   - long */
	/*                       vector of dimension last, the first k */
	/*                       elements of which contain pointers */
	/*                       to the error estimates, such that */
	/*                       elist(iord(1)),... , elist(iord(k)) */
	/*                       form a decreasing sequence, with */
	/*                       k = last if last.le.(limit/2+2), and */
	/*                       k = limit+1-last otherwise */

	/*              nrmax  - long */
	/*                       maxerr = iord(nrmax) */

	/* 4.        no subroutines or functions needed */
	/* ***end prologue  qpsrt */


	/*           check whether the list contains more than */
	/*           two error estimates. */

	/* ***first executable statement  qpsrt */
	/* Parameter adjustments */
	--iord;
	--elist;

	/* Function Body */
	if (*last > 2) {
		goto L10;
	}
	iord[1] = 1;
	iord[2] = 2;
	goto L90;

	/*           this part of the routine is only executed */
	/*           if, due to a difficult integrand, subdivision */
	/*           increased the error estimate. in the normal case */
	/*           the insert procedure should start after the */
	/*           nrmax-th largest error estimate. */

L10:
	errmax = elist[*maxerr];
	if (*nrmax == 1) {
		goto L30;
	}
	ido = *nrmax - 1;
	i__1 = ido;
	for (i__ = 1; i__ <= i__1; ++i__) {
		isucc = iord[*nrmax - 1];
		/* ***jump out of do-loop */
		if (errmax <= elist[isucc]) {
			goto L30;
		}
		iord[*nrmax] = isucc;
		--(*nrmax);
		/* L20: */
	}

	/*           compute the number of elements in the list to */
	/*           be maintained in descending order. this number */
	/*           depends on the number of subdivisions still */
	/*           allowed. */

L30:
	jupbn = *last;
	if (*last > *limit / 2 + 2) {
		jupbn = *limit + 3 - *last;
	}
	errmin = elist[*last];

	/*           insert errmax by traversing the list top-down, */
	/*           starting comparison from the element elist(iord(nrmax+1)). */

	jbnd = jupbn - 1;
	ibeg = *nrmax + 1;
	if (ibeg > jbnd) {
		goto L50;
	}
	i__1 = jbnd;
	for (i__ = ibeg; i__ <= i__1; ++i__) {
		isucc = iord[i__];
		/* ***jump out of do-loop */
		if (errmax >= elist[isucc]) {
			goto L60;
		}
		iord[i__ - 1] = isucc;
		/* L40: */
	}
L50:
	iord[jbnd] = *maxerr;
	iord[jupbn] = *last;
	goto L90;

	/*           insert errmin by traversing the list bottom-up. */

L60:
	iord[i__ - 1] = *maxerr;
	k = jbnd;
	i__1 = jbnd;
	for (j = i__; j <= i__1; ++j) {
		isucc = iord[k];
		/* ***jump out of do-loop */
		if (errmin < elist[isucc]) {
			goto L80;
		}
		iord[k + 1] = isucc;
		--k;
		/* L70: */
	}
	iord[i__] = *last;
	goto L90;
L80:
	iord[k + 1] = *last;

	/*           set maxerr and ermax. */

L90:
	*maxerr = iord[*nrmax];
	*ermax = elist[*maxerr];
	return;
} /* qpsrt_ */

sys_float qwgtc_(const sys_float *x, const sys_float *c__, const sys_float *, 
					  const sys_float *, const sys_float *, const long *)
{
	/* System generated locals */
	sys_float ret_val;

	/* ***begin prologue  qwgtc */
	/* ***refer to qk15w */
	/* ***routines called  (none) */
	/* ***revision date  810101   (yymmdd) */
	/* ***keywords  weight function, cauchy principal value */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  this function subprogram is used together with the */
	/*            routine qawc and defines the weight function. */
	/* ***end prologue  qwgtc */

	/* ***first executable statement */
	ret_val = 1.f / (*x - *c__);
	return ret_val;
} /* qwgtc_ */

sys_float qwgtf_(const sys_float *x, const sys_float *omega, const sys_float *,
					  const sys_float *, const sys_float *, const long *integr)
{
	/* System generated locals */
	sys_float ret_val;

	/* Local variables */
	sys_float omx;

	/* ***begin prologue  qwgtf */
	/* ***refer to   qk15w */
	/* ***routines called  (none) */
	/* ***revision date 810101   (yymmdd) */
	/* ***keywords  cos or sin in weight function */
	/* ***author  piessens,robert, appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. * progr. div. - k.u.leuven */
	/* ***end prologue  qwgtf */

	/* ***first executable statement */
	omx = *omega * *x;
	switch (*integr) {
	case 1:  goto L10;
	case 2:  goto L20;
	}
L10:
	ret_val = cos(omx);
	goto L30;
L20:
	ret_val = sin(omx);
L30:
	return ret_val;
} /* qwgtf_ */

sys_float qwgts_(const sys_float *x, const sys_float *a, const sys_float *b,
					  const sys_float *alfa, const sys_float *beta, 
					  const long *integr)
{
	/* System generated locals */
	sys_float ret_val;
	double d__1, d__2, d__3, d__4;

	/* Local variables */
	sys_float xma, bmx;

	/* ***begin prologue  qwgts */
	/* ***refer to qk15w */
	/* ***routines called  (none) */
	/* ***revision date  810101   (yymmdd) */
	/* ***keywords  weight function, algebraico-logarithmic */
	/*             end-point singularities */
	/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
	/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
	/* ***purpose  this function subprogram is used together with the */
	/*            routine qaws and defines the weight function. */
	/* ***end prologue  qwgts */

	/* ***first executable statement */
	xma = *x - *a;
	bmx = *b - *x;
	d__1 = (double) xma;
	d__2 = (double) (*alfa);
	d__3 = (double) bmx;
	d__4 = (double) (*beta);
	ret_val = pow(d__1, d__2) * pow(d__3, d__4);
	switch (*integr) {
	case 1:  goto L40;
	case 2:  goto L10;
	case 3:  goto L20;
	case 4:  goto L30;
	}
L10:
	ret_val *= log(xma);
	goto L40;
L20:
	ret_val *= log(bmx);
	goto L40;
L30:
	ret_val = ret_val * log(xma) * log(bmx);
L40:
	return ret_val;
} /* qwgts_ */

void sgtsl_(const long *n, sys_float *c__, sys_float *d__, 
				sys_float *e, sys_float *b, long *info)
{
	/* System generated locals */
	int i__1;
	sys_float r__1, r__2;

	/* Local variables */
	long k;
	sys_float t;
	long kb, kp1, nm1, nm2;


	/*     sgtsl given a general tridiagonal matrix and a right hand */
	/*     side will find the solution. */

	/*     on entry */

	/*        n       long */
	/*                is the order of the tridiagonal matrix. */

	/*        c       sys_float(n) */
	/*                is the subdiagonal of the tridiagonal matrix. */
	/*                c(2) through c(n) should contain the subdiagonal. */
	/*                on output c is destroyed. */

	/*        d       sys_float(n) */
	/*                is the diagonal of the tridiagonal matrix. */
	/*                on output d is destroyed. */

	/*        e       sys_float(n) */
	/*                is the superdiagonal of the tridiagonal matrix. */
	/*                e(1) through e(n-1) should contain the superdiagonal. */
	/*                on output e is destroyed. */

	/*        b       sys_float(n) */
	/*                is the right hand side vector. */

	/*     on return */

	/*        b       is the solution vector. */

	/*        info    long */
	/*                = 0 normal value. */
	/*                = k if the k-th element of the diagonal becomes */
	/*                    exactly zero.  the subroutine returns when */
	/*                    this is detected. */

	/*     linpack. this version dated 08/14/78 . */
	/*     jack dongarra, argonne national laboratory. */

	/*     no externals */
	/*     fortran abs */

	/*     internal variables */

	/*     begin block permitting ...exits to 100 */

	/* Parameter adjustments */
	--b;
	--e;
	--d__;
	--c__;

	/* Function Body */
	*info = 0;
	c__[1] = d__[1];
	nm1 = *n - 1;
	if (nm1 < 1) {
		goto L40;
	}
	d__[1] = e[1];
	e[1] = 0.f;
	e[*n] = 0.f;

	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;

		/*              find the largest of the two rows */

		if ((r__1 = c__[kp1], fabs(r__1)) < (r__2 = c__[k], fabs(r__2))) {
			goto L10;
		}

		/*                 interchange row */

		t = c__[kp1];
		c__[kp1] = c__[k];
		c__[k] = t;
		t = d__[kp1];
		d__[kp1] = d__[k];
		d__[k] = t;
		t = e[kp1];
		e[kp1] = e[k];
		e[k] = t;
		t = b[kp1];
		b[kp1] = b[k];
		b[k] = t;
	L10:

		/*              zero elements */

		if (c__[k] != 0.f) {
			goto L20;
		}
		*info = k;
		/*     ............exit */
		goto L100;
	L20:
		t = -c__[kp1] / c__[k];
		c__[kp1] = d__[kp1] + t * d__[k];
		d__[kp1] = e[kp1] + t * e[k];
		e[kp1] = 0.f;
		b[kp1] += t * b[k];
		/* L30: */
	}
L40:
	if (c__[*n] != 0.f) {
		goto L50;
	}
	*info = *n;
	goto L90;
L50:

	/*           back solve */

	nm2 = *n - 2;
	b[*n] /= c__[*n];
	if (*n == 1) {
		goto L80;
	}
	b[nm1] = (b[nm1] - d__[nm1] * b[*n]) / c__[nm1];
	if (nm2 < 1) {
		goto L70;
	}
	i__1 = nm2;
	for (kb = 1; kb <= i__1; ++kb) {
		k = nm2 - kb + 1;
		b[k] = (b[k] - d__[k] * b[k + 1] - e[k] * b[k + 2]) / c__[k];
		/* L60: */
	}
L70:
L80:
L90:
L100:

	return;
} /* sgtsl_ */

void dgtsl_(const long *n, double *c__, double *d__, 
				double *e, double *b, long *info)
{
	/* System generated locals */
	int i__1;
	double d__1, d__2;

	/* Local variables */
	long k;
	double t;
	long kb, kp1, nm1, nm2;


	/*     dgtsl given a general tridiagonal matrix and a right hand */
	/*     side will find the solution. */

	/*     on entry */

	/*        n       long */
	/*                is the order of the tridiagonal matrix. */

	/*        c       double precision(n) */
	/*                is the subdiagonal of the tridiagonal matrix. */
	/*                c(2) through c(n) should contain the subdiagonal. */
	/*                on output c is destroyed. */

	/*        d       double precision(n) */
	/*                is the diagonal of the tridiagonal matrix. */
	/*                on output d is destroyed. */

	/*        e       double precision(n) */
	/*                is the superdiagonal of the tridiagonal matrix. */
	/*                e(1) through e(n-1) should contain the superdiagonal. */
	/*                on output e is destroyed. */

	/*        b       double precision(n) */
	/*                is the right hand side vector. */

	/*     on return */

	/*        b       is the solution vector. */

	/*        info    long */
	/*                = 0 normal value. */
	/*                = k if the k-th element of the diagonal becomes */
	/*                    exactly zero.  the subroutine returns when */
	/*                    this is detected. */

	/*     linpack. this version dated 08/14/78 . */
	/*     jack dongarra, argonne national laboratory. */

	/*     no externals */
	/*     fortran dabs */

	/*     internal variables */

	/*     begin block permitting ...exits to 100 */

	/* Parameter adjustments */
	--b;
	--e;
	--d__;
	--c__;

	/* Function Body */
	*info = 0;
	c__[1] = d__[1];
	nm1 = *n - 1;
	if (nm1 < 1) {
		goto L40;
	}
	d__[1] = e[1];
	e[1] = 0.;
	e[*n] = 0.;

	i__1 = nm1;
	for (k = 1; k <= i__1; ++k) {
		kp1 = k + 1;

		/*              find the largest of the two rows */

		if ((d__1 = c__[kp1], fabs(d__1)) < (d__2 = c__[k], fabs(d__2))) {
			goto L10;
		}

		/*                 interchange row */

		t = c__[kp1];
		c__[kp1] = c__[k];
		c__[k] = t;
		t = d__[kp1];
		d__[kp1] = d__[k];
		d__[k] = t;
		t = e[kp1];
		e[kp1] = e[k];
		e[k] = t;
		t = b[kp1];
		b[kp1] = b[k];
		b[k] = t;
	L10:

		/*              zero elements */

		if (c__[k] != 0.) {
			goto L20;
		}
		*info = k;
		/*     ............exit */
		goto L100;
	L20:
		t = -c__[kp1] / c__[k];
		c__[kp1] = d__[kp1] + t * d__[k];
		d__[kp1] = e[kp1] + t * e[k];
		e[kp1] = 0.;
		b[kp1] += t * b[k];
		/* L30: */
	}
L40:
	if (c__[*n] != 0.) {
		goto L50;
	}
	*info = *n;
	goto L90;
L50:

	/*           back solve */

	nm2 = *n - 2;
	b[*n] /= c__[*n];
	if (*n == 1) {
		goto L80;
	}
	b[nm1] = (b[nm1] - d__[nm1] * b[*n]) / c__[nm1];
	if (nm2 < 1) {
		goto L70;
	}
	i__1 = nm2;
	for (kb = 1; kb <= i__1; ++kb) {
		k = nm2 - kb + 1;
		b[k] = (b[k] - d__[k] * b[k + 1] - e[k] * b[k + 2]) / c__[k];
		/* L60: */
	}
L70:
L80:
L90:
L100:

	return;
} /* dgtsl_ */

