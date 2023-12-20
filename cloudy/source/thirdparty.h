/* This file contains routines (perhaps in modified form) by third parties.
 * Use and distribution of these works are determined by their respective copyrights. */

#ifndef THIRDPARTY_H_
#define THIRDPARTY_H_

#include "physconst.h"

/*============================================================================*/

/* these are the routines in thirdparty.cpp */

bool linfit(
	long n,
	const double x[], /* x[n] */
	const double y[], /* y[n] */
	double &a,
	double &siga,
	double &b,
	double &sigb
);

/** number of predefined factorials, so largest factorial
 * that can be represented as double is (NPRE_FACTORIAL-1)! */
static const int NPRE_FACTORIAL = 171;

/* largest value of fct function corresponding to factorial in six j sjs calculation */
static const int MXDSF = 340;

/** factorial: compute n! by lookup in table of predefined factorials */
double factorial(long n);

/** lfactorial: compute log10(n!), this sroutine cahes its results for efficiency */
double lfactorial(long n);

complex<double> cdgamma(complex<double> x);

double bessel_j0(double x);
double bessel_y0(double x);
double bessel_j1(double x);
double bessel_y1(double x);
double bessel_jn(int n, double x);
double bessel_yn(int n, double x);

double ellpk(double x);

/** expn, returns exponential integral,
 \param n is order, 1 for first integral integral
 \param x is argument, must be positive
 */
double expn(int n, double x);

/** erfce(a) returns erfc(a)*exp(a^2) */
double erfce(double);

double igam (double a, double x);
double igamc (double a, double x);
double igamc_scaled(double a, double x);

double bessel_k0(double x);
double bessel_k0_scaled(double x);
double bessel_k1(double x);
double bessel_k1_scaled(double x);
void bessel_k0_k1(double x, double* k0val, double* k1val);
void bessel_k0_k1_scaled(double x, double* k0val, double* k1val);

double bessel_i0(double x);
double bessel_i0_scaled(double x);
double bessel_i1(double x);
double bessel_i1_scaled(double x);
void bessel_i0_i1(double x, double* k0val, double* k1val);
void bessel_i0_i1_scaled(double x, double* k0val, double* k1val);

#ifndef HAVE_SINCOS
// this is a GNU extension
inline void sincos(double x, double* s, double* c)
{
	*s = sin(x);
	*c = cos(x);
}
#endif

/**e1 first exponential integral
 \param x optical depth argument */
double e1(double x);

/**e1_scaled is exp(x)*e1(x) */
double e1_scaled(double x);

/**e2 second exponential integral 
 \param t optical depth argument */
double e2(double t);

/* random number generator written by David Blackman and Sebastiano Vigna (vigna@acm.org) */
void xoroshiro128plus(uint64* pool, size_t size, uint64 state[], size_t ns);
/* this call is equivalent to skipping ahead 2^64 random variates in the sequence */
void xoroshiro128plus_jump(uint64& state0, uint64& state1);
/* another random number generator written by David Blackman and Sebastiano Vigna (vigna@acm.org) */
void xoshiro256starstar(uint64* pool, size_t size, uint64 state[], size_t ns);
/* this call is equivalent to skipping ahead 2^128 random variates in the sequence */
void xoshiro256starstar_jump(uint64& state0, uint64& state1, uint64& state2, uint64& state3);

/* simple random number generator written by Sebastiano Vigna (vigna@acm.org) */
uint64 splitmix64(uint64& state);

/*============================================================================*/

/* these are the routines in thirdparty_interpolate.cpp */

void spline_cubic_set( long n, const double t[], const double y[], double ypp[],
		       int ibcbeg, double ybcbeg, int ibcend, double ybcend );
void spline_cubic_val( long n, const double t[], double tval, const double y[], const double ypp[],
		       double *yval, double *ypval, double *yppval );

/** spline: set of routines to do spline interpolation
 * call spline first to set coefficients up,
 * then call splint to interpolate
 * spldrv gives dy/dx(x) rather than y(x)
 * splint_safe and spldrv_safe check whether x is within bounds
\param x[] - array of x values given
\param y[] - array of y values given
\param n - number of points in these arrays
\param yp1 some sort of boundary condition, set to 2e31
\param yp2 some sort of boundary condition, set to 2e31
\param y2a[] - array n long to save coefficients
*/
inline void spline(const double x[], 
		   const double y[], 
		   long int n, 
		   double yp1, 
		   double ypn, 
		   double y2a[])
{
	int ibcbeg = ( yp1 > 0.99999e30 ) ? 2 : 1;
	double ybcbeg = ( yp1 > 0.99999e30 ) ? 0. : yp1;
	int ibcend = ( ypn > 0.99999e30 ) ? 2 : 1;
	double ybcend = ( ypn > 0.99999e30 ) ? 0. : ypn;
	spline_cubic_set( n, x, y, y2a, ibcbeg, ybcbeg, ibcend, ybcend );
	return;
}

/** splint: do spline interpolation
\param x[] - array of x values given
\param y[] - array of y valus
\param y2a[] - array of save values found above
\param n - number of points in these arrays
\param x value where x is desired
\param y2a[] - array n long to save coefficients
*/
inline void splint(const double xa[], 
		   const double ya[], 
		   const double y2a[], 
		   long int n, 
		   double x, 
		   double *y)
{
	spline_cubic_val( n, xa, x, ya, y2a, y, NULL, NULL );
	return;
}

inline void spldrv(const double xa[], 
		   const double ya[], 
		   const double y2a[], 
		   long int n, 
		   double x, 
		   double *y)
{
	spline_cubic_val( n, xa, x, ya, y2a, NULL, y, NULL );
	return;
}

/* wrapper routine for splint that checks whether x-value is within bounds
 * if the x-value is out of bounds, a flag will be raised and the function
 * will be evaluated at the nearest boundary */
/* >>chng 03 jan 15, added splint_safe, PvH */
inline void splint_safe(const double xa[], 
			const double ya[], 
			const double y2a[], 
			long int n, 
			double x, 
			double *y,
			bool *lgOutOfBounds)
{
	double xsafe;

	const double lo_bound = MIN2(xa[0],xa[n-1]);
	const double hi_bound = MAX2(xa[0],xa[n-1]);
	const double SAFETY = MAX2(hi_bound-lo_bound,1.)*10.*DBL_EPSILON;

	DEBUG_ENTRY( "splint_safe()" );

	if( x < lo_bound-SAFETY )
	{
		xsafe = lo_bound;
		*lgOutOfBounds = true;
	}
	else if( x > hi_bound+SAFETY )
	{
		xsafe = hi_bound;
		*lgOutOfBounds = true;
	}
	else
	{
		xsafe = x;
		*lgOutOfBounds = false;
	}

	splint(xa,ya,y2a,n,xsafe,y);
	return;
}

/* wrapper routine for spldrv that checks whether x-value is within bounds
 * if the x-value is out of bounds, a flag will be raised and the function
 * will be evaluated at the nearest boundary */
/* >>chng 03 jan 15, added spldrv_safe, PvH */
inline void spldrv_safe(const double xa[],
			const double ya[],
			const double y2a[],
			long int n,
			double x,
			double *y,
			bool *lgOutOfBounds)
{
	double xsafe;

	const double lo_bound = MIN2(xa[0],xa[n-1]);
	const double hi_bound = MAX2(xa[0],xa[n-1]);
	const double SAFETY = MAX2(fabs(lo_bound),fabs(hi_bound))*10.*DBL_EPSILON;

	DEBUG_ENTRY( "spldrv_safe()" );

	if( x < lo_bound-SAFETY )
	{
		xsafe = lo_bound;
		*lgOutOfBounds = true;
	}
	else if( x > hi_bound+SAFETY )
	{
		xsafe = hi_bound;
		*lgOutOfBounds = true;
	}
	else
	{
		xsafe = x;
		*lgOutOfBounds = false;
	}

	spldrv(xa,ya,y2a,n,xsafe,y);
	return;
}

/** lagrange: do lagrange interpolation of order n on x[], y[]
 * use with caution, especialy for high order n!
 * using spline interpolation above is preferred
\param x[]
\param y[]
\param n
\param xval
*/
double lagrange(const double x[], /* x[n] */
		const double y[], /* y[n] */
		long n,
		double xval);

/** find index ilo such that x[ilo] <= xval < x[ilo+1] using bisection
 * this version implicitly assumes that x is monotonically increasing */
template<class T>
inline long hunt_bisect(const T x[], /* x[n] */
			long n,
			T xval)
{
	/* do bisection hunt */
	long ilo = 0, ihi = n-1;
	while( ihi-ilo > 1 )
	{
		long imid = (ilo+ihi)/2;
		if( xval < x[imid] )
			ihi = imid;
		else
			ilo = imid;
	}
	return ilo;
}

/** vector<> versions of the above */
template<class T>
inline long hunt_bisect(const vector<T>& x,
			T xval)
{
	return hunt_bisect(get_ptr(x), x.size(), xval);
}

template<class T>
inline long hunt_bisect(const vector<T>& x,
			size_t n,
			T xval)
{
	ASSERT( n <= x.size() );
	return hunt_bisect(get_ptr(x), n, xval);
}

/** find index ilo such that x[ilo] >= xval > x[ilo+1] using bisection
 * this version implicitly assumes that x is monotonically decreasing */
template<class T>
inline long hunt_bisect_reverse(const T x[], /* x[n] */
				long n,
				T xval)
{
	/* do bisection hunt */
	long ilo = 0, ihi = n-1;
	while( ihi-ilo > 1 )
	{
		long imid = (ilo+ihi)/2;
		if( xval <= x[imid] )
			ilo = imid;
		else
			ihi = imid;
	}
	return ilo;
}

/** vector<> versions of the above */
template<class T>
inline long hunt_bisect_reverse(const vector<T>& x,
				T xval)
{
	return hunt_bisect_reverse(get_ptr(x), x.size(), xval);
}

template<class T>
inline long hunt_bisect_reverse(const vector<T>& x,
				size_t n,
				T xval)
{
	ASSERT( n <= x.size() );
	return hunt_bisect_reverse(get_ptr(x), n, xval);
}

/** do linear interpolation on x[], y[];
 * it is assumed that x[] is strictly monotonically increasing
\param x[]
\param y[]
\param n
\param xval
*/
template<class T>
T linint(const T x[], /* x[n] */
	 const T y[], /* y[n] */
	 long n,
	 T xval)
{
	T yval;

	ASSERT( n >= 2 );

	if( xval <= x[0] )
		yval = y[0];
	else if( xval >= x[n-1] )
		yval = y[n-1];
	else
	{
		long ilo = hunt_bisect( x, n, xval );
		T deriv = (y[ilo+1]-y[ilo])/(x[ilo+1]-x[ilo]);
		yval = y[ilo] + deriv*(xval-x[ilo]);
	}
	return yval;
}

/*============================================================================*/

/* these are the routines in thirdparty_lapack.cpp */

/* there are wrappers for lapack linear algebra routines.
 * there are two versions of the lapack routines - a fortran
 * version that I converted to C with forc to use if nothing else is available
 * (included in the Cloudy distribution),
 * and an option to link into an external lapack library that may be optimized
 * for your machine.  By default the tralated C routines will be used.
 * To use your machine's lapack library instead, define the macro
 * LAPACK and link into your library.  This is usually done with a command line
 * switch "-DLAPACK" on the compiler command, and the linker option "-llapack"
 */

/**getrf_wrapper return value is zero for success, non-zero is error condition 
\param M
\param N
\param *A
\param lda
\param *ipiv
\param *info
*/
void getrf_wrapper(long M, long N, double *A, long lda, int32 *ipiv, int32 *info);

/**getrs_wrapper return value is zero for success, non-zero is error condition 
\param trans
\param N
\param nrhs
\param *A
\param lda
\param *ipiv
\param *B
\param ldb
\param *info
*/
void getrs_wrapper(char trans, long N, long nrhs, double *A, long lda, int32 *ipiv, double *B, long ldb, int32 *info);

double dawson(double x);

void humlik(int n, const realnum x[], realnum y, realnum k[]);

realnum FastVoigtH(realnum a, realnum v);
void FastVoigtH(realnum a, const realnum v[], realnum y[], size_t n);

// calculates y[i] = H(a,v[i]) as defined in Eq. 9-44 of Mihalas
inline void VoigtH(realnum a, const realnum v[], realnum y[], int n)
{
	if( a <= 0.1f )
	{
		FastVoigtH( a, v, y, n );
	}
	else
	{
		humlik( n, v, a, y );
	}
}

// calculates y[i] = U(a,v[i]) as defined in Eq. 9-45 of Mihalas
inline void VoigtU(realnum a, const realnum v[], realnum y[], int n)
{
	VoigtH( a, v, y, n );
	for( int i=0; i < n; ++i )
		y[i] /= realnum(SQRTPI);
}

// VoigtH0(a) returns the value H(a,0) following Eq. 9-44 of Mihalas
inline double VoigtH0(double a)
{
	return erfce(a);
}

// VoigtU0(a) returns the value U(a,0) following Eq. 9-45 of Mihalas
inline double VoigtU0(double a)
{
	return erfce(a)/SQRTPI;
}

/** the size of an MD5 checksum in characters */
static const unsigned int NMD5 = 32;

/** calculate the MD5 sum of a file */
string MD5file(const char* fnam, access_scheme scheme=AS_DEFAULT);
/** non-standard MD5 algorithm that skips eol characters and comments lines */
string MD5datafile(const char* fnam, access_scheme scheme=AS_DEFAULT);
/** same as MD5datafile(), but operates on an already open fstream */
string MD5datastream(fstream& ioFile);
/** calculate the MD5 sum of a string */
string MD5string(const string& str);
void MD5string(const string& str, uint64 md5sum[2]);

void test_expn();
void chbfit(double a, double b, vector<double>& c, double (*func)(double));

void svd(const int nm, const int m, const int n, double *a, 
		  double *w, bool matu, double *u, bool matv, 
		  double *v, int *ierr, double *rv1);

// Evaluate the Gegenbauer (aka ultraspherical) polynomial C_n^(alpha)(x)
double gegenbauer(long n, double al, double x);

// Diagonal generating function for ultraspherical polynomials,
// following Badnell analysis.  Integer argument to constructor is
// sum of raised and lowered index to the ultraspherical function.
class UltraGen
{
	long   m_a, m_n;
	double m_x, m_c1, m_c2;
public:
UltraGen(long a, double x) : m_a(a), m_n(0), m_x(x), m_c1(0.), m_c2(1.) 
	{}
	void step()
	{
		double c1 = m_c1;
		--m_a;
		m_c1 = 2.*m_a*(m_c2-m_x*c1)/(m_n+2*m_a);
		m_c2 = 2.*m_a*(m_x*m_c2-c1)/(m_n+1);
		++m_n;
	}
	double val() const
	{
		return m_c2;
	}
};

// this is the triangular inequality for angular momenta in 2*J format
inline bool Triangle2(long J1, long J2, long J3)
{
	return ( J1 >= 0 && J2 >= 0 && J3 >= 0 &&
		 J1 >= abs(J2-J3) && J1 <= J2+J3 &&
		 J1%2 == (J2+J3)%2 );
}

// 6j Wigner evaluation, original routine Fortran Nigel Badnell 6j for Autostructure
double sjs(long j1, long j2, long j3, long l1, long l2, long l3);
// version of the above with less range in j, is exported for testing purposes
double SixJFull(long j1, long j2, long j3, long j4, long j5, long j6);
// Implementation utlizing recursion relation for vector of sixj values,
// should be robust to high j
void rec6j(double *sixcof, double l2, double l3, 
			  double l4, double l5, double l6, double *l1min,
			  double *l1max, double *lmatch, long ndim, long *ier);

// the following code was taken (and altered) from:
// http://www.geeksforgeeks.org/find-all-shortest-unique-prefixes-to-represent-each-word-in-a-given-list/
static const size_t TRIESZ = 128;

struct trieNode
{
	trieNode* child[TRIESZ];
	size_t freq;  // To store frequency
	trieNode()
	{
		memset(child, 0, TRIESZ*sizeof(trieNode*));
		freq = 0;
	}
	trieNode(const trieNode&) = delete;
	trieNode& operator= (const trieNode&) = delete;
	~trieNode()
	{
		for( size_t i=0; i < TRIESZ; i++ )
			delete child[i];
	}
};
 
void insertToken(trieNode* root, const string& token);
size_t findUniqueLen(trieNode* root, const string& token);

// the following code was taken (and altered) from:
// https://en.wikipedia.org/wiki/Levenshtein_distance

size_t LevenshteinDistance(const string& s, const string& t);

#endif /* THIRDPARTY_H_ */
