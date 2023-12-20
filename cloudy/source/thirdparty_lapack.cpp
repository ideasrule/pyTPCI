/* These are wrappers for lapack linear algebra routines.
 * There are two versions of the lapack routines - a fortran
 * version that I converted to C with forc to use if nothing else is available
 * (included in the Cloudy distribution),
 * and an option to link into an external lapack library that may be optimized
 * for your machine.  By default the tralated C routines will be used.
 * To use your machine's lapack library instead, define the macro
 * LAPACK and link into your library.  This is usually done with a command line
 * switch "-DLAPACK" on the compiler command, and the linker option "-llapack"
 */
#include "cddefines.h"
#include "thirdparty.h"
/*lint -e725 expected pos indentation */
/*lint -e801 goto is deprecated */

#ifdef LAPACK
/*********************The functions for FORTRAN version of the LAPACK calls *******************/
/* dgetrf, dgetrs: lapack FORTRAN general full matrix solution using LU decomposition */

extern "C" void dgetrf_(int32 *M, int32 *N, double *A, int32 *LDA, int32 *IPIV, int32 *INFO);
extern "C" void dgetrs_(char *TRANS, int32 *N, int32 *NRHS, double *A, int32 *LDA, int32 *iPiv, double *B,
	     int32 *LDB, int32 *INFO, int32 translen);
extern "C" void dgtsv_(int32 *n, int32 *nrhs, double *dl, double *d, double *du, double *b, int32 *ldb, int32 *info);

#else

/*********************The functions for C version of the LAPACK calls *******************/
/*
 * these are the routines that are, part of lapack, Some had their names slightly
 * changed so as to not conflict with the special lapack that exists on our exemplar.
 */

/* DGETRF lapack, perform LU decomposition */
STATIC void DGETRF(int32,int32,double*,int32,int32[],int32*);

/* DGETRS lapack, solve linear system */
STATIC void DGETRS(int32 TRANS,int32 N,int32 NRHS,double *A,int32 LDA,int32 IPIV[],double *B,int32 LDB,int32 *INFO);

/* DGTSV lapack, solve tridiagonal system */
/*static int32 DGTSV(int32 *n,int32 *nrhs,double *dl,double *d__,double *du,double *b,int32 *ldb,int32 *info);*/

#endif


/**********************************************************************************************/
/* returns zero if successful termination */
void getrf_wrapper(long M, long N, double *A, long lda, int32 *ipiv, int32 *info)
{
	if( *info == 0 )
	{
		int32 M_loc, N_loc, lda_loc;

		ASSERT( M < INT32_MAX && N < INT32_MAX && lda < INT32_MAX );

		M_loc = (int32)M;
		N_loc = (int32)N;
		lda_loc = (int32)lda;

#		ifdef  LAPACK
		/* Calling the special version in library */
		dgetrf_(&M_loc, &N_loc, A , &lda_loc, ipiv, info);
#		else
		/* Calling the old slower one, included with cloudy */
		DGETRF(M_loc, N_loc, A, lda_loc, ipiv, info);
#		endif
	}
}

void getrs_wrapper(char trans, long N, long nrhs, double *A, long lda, int32 *ipiv, 
		   double *B, long ldb, int32 *info)
{
	if( *info == 0 )
	{
		int32 N_loc, nrhs_loc, lda_loc, ldb_loc;

		ASSERT( N < INT32_MAX && nrhs < INT32_MAX && lda < INT32_MAX && ldb < INT32_MAX );

		N_loc = (int32)N;
		nrhs_loc = (int32)nrhs;
		lda_loc = (int32)lda;
		ldb_loc = (int32)ldb;

#		ifdef LAPACK
		/* Calling the special version in library */
		dgetrs_(&trans, &N_loc, &nrhs_loc, A, &lda_loc, ipiv, B, &ldb_loc, info, sizeof(char));
#		else
		/* Calling the old slower one, included with cloudy */
		DGETRS(trans, N_loc, nrhs_loc, A, lda_loc, ipiv, B, ldb_loc, info);
#		endif
	}
}

#if 0
void dgtsv_wrapper(long N, long nrhs, double *dl, double *d__, double *du, double *b, long ldb, int32 *info)
{
	printf("Inside dgtsv\n");
	cdEXIT( EXIT_FAILURE );
	if( *info == 0 )
	{
		int32 N_loc, nrhs_loc, ldb_loc;

		ASSERT( N < INT32_MAX && nrhs < INT32_MAX && ldb < INT32_MAX );

		N_loc = (int32)N;
		nrhs_loc = (int32)nrhs;
		ldb_loc = (int32)ldb;

#		ifdef LAPACK
		/* Calling the special version in library */
		dgtsv_(&N_loc, &nrhs_loc, dl, d__, du, b, &ldb_loc, info);
#		else
		/* Calling the old slower one, included with cloudy */
		/* DGTSV always returns zero, so it is safe to ignore the return value */
		(void)DGTSV(&N_loc, &nrhs_loc, dl, d__, du, b, &ldb_loc, info);
#		endif
	}
}
#endif


#ifndef LAPACK

#define	ONE	1.0e0
#define	ZERO	0.0e0

#ifdef AA
#	undef AA
#endif
#ifdef BB
#	undef BB
#endif
#ifdef CC
#	undef CC
#endif

#define AA(I_,J_)	(*(A+(I_)*(LDA)+(J_)))
#define BB(I_,J_)	(*(B+(I_)*(LDB)+(J_)))
#define CC(I_,J_)	(*(C+(I_)*(LDC)+(J_)))

/* 
 * these are the routines that are, part of lapack, Some had their names slightly
 * changed so as to not conflict with the special lapack that exits on our exemplar.
 */

/* dgtsv, dgetrf, dgetrs: lapack general tridiagonal solution */
/*int32 DGTSV(int32 *n, int32 *nrhs, double *dl, 
	double *d__, double *du, double *b, int32 *ldb, int32 *info);*/

/* DGETRF lapack service routine */
/*void DGETRF(int32,int32,double*,int32,int32[],int32*);*/

/*DGETRS lapack matrix inversion routine */
/*void DGETRS(int32 TRANS, 
	  int32 N, 
	  int32 NRHS, 
	  double *A, 
	  int32 LDA, 
	  int32 IPIV[], 
	  double *B, 
	  int32 LDB, 
	  int32 *INFO);
*/
/* DGEMM matrix inversion helper routine*/
STATIC void DGEMM(int32 TRANSA, 
	  int32 TRANSB, 
	  int32 M, 
	  int32 N, 
	  int32 K, 
	  double ALPHA, 
	  double *A, 
	  int32 LDA, 
	  double *B, 
	  int32 LDB, 
	  double BETA, 
	  double *C, 
	  int32 LDC);

/*LSAME LAPACK auxiliary routine  */
STATIC int32 LSAME(int32 CA, 
	  int32 CB);

/*IDAMAX lapack service routine */
STATIC int32 IDAMAX(int32 n, 
	  double dx[], 
	  int32 incx);

/*DTRSM lapack service routine */
STATIC void DTRSM(int32 SIDE, 
	  int32 UPLO, 
	  int32 TRANSA, 
	  int32 DIAG, 
	  int32 M, 
	  int32 N, 
	  double ALPHA, 
	  double *A, 
	  int32 LDA, 
	  double *B, 
	  int32 LDB);

/* ILAENV lapack helper routine */
STATIC int32 ILAENV(int32 ISPEC, 
	  const char *NAME, 
	  /*char *OPTS, */
	  int32 N1, 
	  int32 N2, 
	  /*int32 N3, */
	  int32 N4);

/*DSWAP lapack routine */
STATIC void DSWAP(int32 n, 
	  double dx[], 
	  int32 incx, 
	  double dy[], 
	  int32 incy);

/*DSCAL lapack routine */
STATIC void DSCAL(int32 n, 
	  double da, 
	  double dx[], 
	  int32 incx);

/*DLASWP  -- LAPACK auxiliary routine (version 2.0) --*/
STATIC void DLASWP(int32 N, 
	  double *A, 
	  int32 LDA, 
	  int32 K1, 
	  int32 K2, 
	  int32 IPIV[], 
	  int32 INCX);

/*DGETF2 lapack service routine */
STATIC void DGETF2(int32 M, 
	  int32 N, 
	  double *A, 
	  int32 LDA, 
	  int32 IPIV[], 
	  int32 *INFO);

/*DGER service routine for matrix inversion */
STATIC void DGER(int32 M, 
	  int32 N, 
	  double ALPHA, 
	  double X[], 
	  int32 INCX, 
	  double Y[], 
	  int32 INCY, 
	  double *A, 
	  int32 LDA);

/*XERBLA  -- LAPACK auxiliary routine (version 2.0) -- */
STATIC void XERBLA(const char *SRNAME, 
		   int32 INFO)
{

	DEBUG_ENTRY( "XERBLA()" );

	/*  -- LAPACK auxiliary routine (version 2.0) --
	 * Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
	 * Courant Institute, Argonne National Lab, and Rice University
	 * September 30, 1994
	 *
	 * .. Scalar Arguments .. */
	/* ..
	 *
	 *  Purpose
	 *  =======
	 *
	 *  XERBLA  is an error handler for the LAPACK routines.
	 *  It is called by an LAPACK routine if an input parameter has an
	 *  invalid value.  A message is printed and execution stops.
	 *
	 *  Installers may consider modifying the STOP statement in order to
	 *  call system-specific exception-handling facilities.
	 *
	 *  Arguments
	 *  =========
	 *
	 *  SRNAME  (input) CHARACTER*6
	 *  The name of the routine which called XERBLA.
	 *
	 *  INFO    (input) INTEGER
	 *  The position of the invalid parameter in the parameter list
	 *  of the calling routine.
	 *
	 * =====================================================================
	 *
	 * .. Executable Statements ..
	 * */
	fprintf( ioQQQ, " ** On entry to %6.6s parameter number %2ld had an illegal value\n",
		 SRNAME, (long)INFO );

	cdEXIT(EXIT_FAILURE);
}


STATIC void DGETRF(
		/* number of rows of the matrix */
	  int32 M, 
		/* number of columns of the matrix
		 * M=N for square matrix */
	  int32 N, 
	  /* double precision matrix */
	  double *A, 
	  /* LDA is right dimension of matrix */
	  int32 LDA, 
	  /* following must dimensions the smaller of M or N */
	  int32 IPIV[], 
	  /* following is zero for successful exit */
	  int32 *INFO)
{

	char chL1, 
	  chL2, 
	  chL3, 
	  chL4;
	int32 I, 
	  IINFO, 
	  I_, 
	  J, 
	  JB, 
	  J_, 
	  NB, 

	  limit, 
	  limit2;
	/*double _d_l;*/

	DEBUG_ENTRY( "DGETRF()" );

	/*  -- LAPACK routine (version 2.0) --
	 * Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
	 * Courant Institute, Argonne National Lab, and Rice University
	 * March 31, 1993
	 *
	 *  Purpose
	 *  =======
	 *
	 *  DGETRF computes an LU factorization of a general M-by-N matrix A
	 *  using partial pivoting with row interchanges.
	 *
	 *  The factorization has the form
	 * A = P * L * U
	 *  where P is a permutation matrix, L is lower triangular with unit
	 *  diagonal elements (lower trapezoidal if m > n), and U is upper
	 *  triangular (upper trapezoidal if m < n).
	 *
	 *  This is the right-looking Level 3 BLAS version of the algorithm.
	 *
	 *  Arguments
	 *  =========
	 *
	 *  M       (input) INTEGER
	 *  The number of rows of the matrix A.  M >= 0.
	 *
	 *  N       (input) INTEGER
	 *  The number of columns of the matrix A.  N >= 0.
	 *
	 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
	 *  On entry, the M-by-N matrix to be factored.
	 *  On exit, the factors L and U from the factorization
	 *  A = P*L*U; the unit diagonal elements of L are not stored.
	 *
	 *  LDA     (input) INTEGER
	 *  The leading dimension of the array A.  LDA >= MAX(1,M).
	 *
	 *  IPIV    (output) INTEGER array, dimension (MIN(M,N))
	 *  The pivot indices; for 1 <= i <= MIN(M,N), row i of the
	 *  matrix was interchanged with row IPIV(i).
	 *
	 *  INFO    (output) INTEGER
	 *  = 0:  successful exit
	 *  < 0:  if INFO = -i, the i-th argument had an illegal value
	 *  > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
	 *      has been completed, but the factor U is exactly
	 *      singular, and division by zero will occur if it is used
	 *      to solve a system of equations.
	 *
	 *  =====================================================================
	 *
	 * .. Parameters .. */
	/* ..
	 * .. Local Scalars .. */
	/* ..
	 * .. External Subroutines .. */
	/* ..
	 * .. External Functions .. */
	/* ..
	 * .. Intrinsic Functions .. */
	/* ..
	 * .. Executable Statements ..
	 *
	 * Test the input parameters.
	 * */
	*INFO = 0;
	if( M < 0 )
	{
		*INFO = -1;
	}
	else if( N < 0 )
	{
		*INFO = -2;
	}
	else if( LDA < MAX2(1,M) )
	{
		*INFO = -4;
	}
	if( *INFO != 0 )
	{
		XERBLA("DGETRF",-*INFO);
		/* XERBLA does not return */
	}

	/* Quick return if possible
	 * */
	if( M == 0 || N == 0 )
	{ 
		return;
	}

	/* Determine the block size for this environment.
	 * */
		/* >>chng 01 oct 22, remove two parameters since not used */
	NB = ILAENV(1,"DGETRF",/*" ",*/M,N,/*-1,*/-1);
	if( NB <= 1 || NB >= MIN2(M,N) )
	{
		/*  Use unblocked code.
		 * */
		DGETF2(M,N,A,LDA,IPIV,INFO);
	}
	else
	{

		/*  Use blocked code.
		 * */
		limit = MIN2(M,N);
		/*for( J=1, _do0=DOCNT(J,limit,_do1 = NB); _do0 > 0; J += _do1, _do0-- )*/
		/*do J = 1, limit , NB */
		for( J=1; J<=limit; J += NB )
		{
			J_ = J - 1;
			JB = MIN2(MIN2(M,N)-J+1,NB);

			/* Factor diagonal and subdiagonal blocks and test for exact
			 * singularity.
			 * */
			DGETF2(M-J+1,JB,&AA(J_,J_),LDA,&IPIV[J_],&IINFO);

			/* Adjust INFO and the pivot indices.
			 * */
			if( *INFO == 0 && IINFO > 0 )
				*INFO = IINFO + J - 1;
			limit2 = MIN2(M,J+JB-1);
			for( I=J; I <= limit2; I++ )
			{
				I_ = I - 1;
				IPIV[I_] += J - 1;
			}

			/* Apply interchanges to columns 1:J-1.
			 * */
			DLASWP(J-1,A,LDA,J,J+JB-1,IPIV,1);

			if( J + JB <= N )
			{

				/*    Apply interchanges to columns J+JB:N.
				 * */
				DLASWP(N-J-JB+1,&AA(J_+JB,0),LDA,J,J+JB-1,IPIV,1);

				/*    Compute block row of U.
				 * */
				chL1 = 'L';
				chL2 = 'L';
				chL3 = 'N';
				chL4 = 'U';
				/*    CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, */
				DTRSM(chL1,chL2,chL3,chL4,JB,N-J-JB+1,ONE,&AA(J_,J_),
				  LDA,&AA(J_+JB,J_),LDA);
				if( J + JB <= M )
				{

					/*       Update trailing submatrix.
					 * */
					chL1 = 'N';
					chL2 = 'N';
					/*       CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, */
					DGEMM(chL1,chL2,M-J-JB+1,N-J-JB+1,JB,-ONE,&AA(J_,J_+JB),
					  LDA,&AA(J_+JB,J_),LDA,ONE,&AA(J_+JB,J_+JB),LDA);
				}
			}
		}
	}
	return;

	/* End of DGETRF
	 * */
#undef	A
}

/*DGETRS lapack matrix inversion routine */
/*****************************************************************
 *****************************************************************
 *
 * matrix inversion routine
 *
 * solves Ax=B  A is an nXn matrix, C and B are nX1 matrices
 * dim A(n,n), B(n,1)  C overwrites B.
 * integer ipiv(n)
 * integer info  see below: 
 *  INFO    (output) INTEGER
 *  = 0:  successful exit
 *  < 0:  if INFO = -i, the i-th argument had an illegal value
 *  > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
 *      has been completed, but the factor U is exactly
 *      singular, and division by zero will occur if it is used
 *      to solve a system of equations.
 *
 *
 *
 *must call the routines in the following order:
 * call dgetrf(n,n,A,n,ipiv,info)
 * call dgetrs('N',n,1,A,n,ipiv,B,n,info)
 *
 *
 *************************************************************************** */

STATIC void DGETRS(
	  /* 1 ch var saying what to do */
	  int32 TRANS, 
	  /* order of the matrix */
	  int32 N, 
	  /* number of right hand sides */
	  int32 NRHS, 
	  /* double [N][LDA] */
	  double *A, 
	  /* second dim of array */
	  int32 LDA, 
	  /* helper vector, dimensioned N*/
	  int32 IPIV[], 
	  /* on input the ri=hs vector, on output, the result */
	  double *B, 
	  /* dimension of B, 1 if one vector */
	  int32 LDB, 
	  /* = 0 if ok */
	  int32 *INFO)
{
/*#define A(I_,J_)	(*(A+(I_)*(LDA)+(J_)))*/
/*#define B(I_,J_)	(*(B+(I_)*(LDB)+(J_)))*/
	int32 NOTRAN;
	char chL1, 
	  chL2, 
	  chL3, 
	  chL4;

	DEBUG_ENTRY( "DGETRS()" );

	/*  -- LAPACK routine (version 2.0) --
	 * Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
	 * Courant Institute, Argonne National Lab, and Rice University
	 * March 31, 1993
	 *
	 *
	 *  Purpose
	 *  =======
	 *
	 *  DGETRS solves a system of linear equations
	 * A * X = B  or  A' * X = B
	 *  with a general N-by-N matrix A using the LU factorization computed
	 *  by DGETRF.
	 *
	 *  Arguments
	 *  =========
	 *
	 *  TRANS   (input) CHARACTER*1
	 *  Specifies the form of the system of equations:
	 *  = 'N':  A * X = B  (No transpose)
	 *  = 'T':  A'* X = B  (Transpose)
	 *  = 'C':  A'* X = B  (Conjugate transpose = Transpose)
	 *
	 *  N       (input) INTEGER
	 *  The order of the matrix A.  N >= 0.
	 *
	 *  NRHS    (input) INTEGER
	 *  The number of right hand sides, i.e., the number of columns
	 *  of the matrix B.  NRHS >= 0.
	 *
	 *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
	 *  The factors L and U from the factorization A = P*L*U
	 *  as computed by DGETRF.
	 *
	 *  LDA     (input) INTEGER
	 *  The leading dimension of the array A.  LDA >= MAX(1,N).
	 *
	 *  IPIV    (input) INTEGER array, dimension (N)
	 *  The pivot indices from DGETRF; for 1<=i<=N, row i of the
	 *  matrix was interchanged with row IPIV(i).
	 *
	 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
	 *  On entry, the right hand side matrix B.
	 *  On exit, the solution matrix X.
	 *
	 *  LDB     (input) INTEGER
	 *  The leading dimension of the array B.  LDB >= MAX(1,N).
	 *
	 *  INFO    (output) INTEGER
	 *  = 0:  successful exit
	 *  < 0:  if INFO = -i, the i-th argument had an illegal value
	 *
	 *  =====================================================================
	 *
	 * .. Parameters .. */
	/* ..
	 * .. Local Scalars .. */
	/* ..
	 * .. External Functions .. */
	/* ..
	 * .. External Subroutines .. */
	/* ..
	 * .. Intrinsic Functions .. */
	/* ..
	 * .. Executable Statements ..
	 *
	 * Test the input parameters.
	 * */
	*INFO = 0;
	NOTRAN = LSAME(TRANS,'N');
	if( (!NOTRAN && !LSAME(TRANS,'T')) && !LSAME(TRANS,'C') )
	{
		*INFO = -1;
	}
	else if( N < 0 )
	{
		*INFO = -2;
	}
	else if( NRHS < 0 )
	{
		*INFO = -3;
	}
	else if( LDA < MAX2(1,N) )
	{
		*INFO = -5;
	}
	else if( LDB < MAX2(1,N) )
	{
		*INFO = -8;
	}
	if( *INFO != 0 )
	{
		XERBLA("DGETRS",-*INFO);
		/* XERBLA does not return */
	}

	/* Quick return if possible
	 * */
	if( N == 0 || NRHS == 0 )
	{ 
		return;
	}

	if( NOTRAN )
	{

		/*  Solve A * X = B.
		 *
		 *  Apply row interchanges to the right hand sides.
		 * */
		DLASWP(NRHS,B,LDB,1,N,IPIV,1);

		/*  Solve L*X = B, overwriting B with X.
		 * */
		chL1 = 'L';
		chL2 = 'L';
		chL3 = 'N';
		chL4 = 'U';
		/*  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, */
		DTRSM(chL1,chL2,chL3,chL4,N,NRHS,ONE,A,LDA,B,LDB);

		/*  Solve U*X = B, overwriting B with X.
		 * */
		chL1 = 'L';
		chL2 = 'U';
		chL3 = 'N';
		chL4 = 'N';
		/*  CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, */
		DTRSM(chL1,chL2,chL3,chL4,N,NRHS,ONE,A,LDA,B,LDB);
	}
	else
	{

		/*  Solve A' * X = B.
		 *
		 *  Solve U'*X = B, overwriting B with X.
		 * */
		chL1 = 'L';
		chL2 = 'U';
		chL3 = 'T';
		chL4 = 'N';
		/*  CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, */
		DTRSM(chL1,chL2,chL3,chL4,N,NRHS,ONE,A,LDA,B,LDB);

		/*  Solve L'*X = B, overwriting B with X.
		 * */
		chL1 = 'L';
		chL2 = 'L';
		chL3 = 'T';
		chL4 = 'U';
		/*  CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, */
		DTRSM(chL1,chL2,chL3,chL4,N,NRHS,ONE,A,LDA,B,LDB);

		/*  Apply row interchanges to the solution vectors.
		 * */
		DLASWP(NRHS,B,LDB,1,N,IPIV,-1);
	}

	return;

	/* End of DGETRS
	 * */
#undef	B
#undef	A
}

/*LSAME LAPACK auxiliary routine  */

STATIC int32 LSAME(int32 CA, 
	  int32 CB)
{
	/*
	 *
	 *  -- LAPACK auxiliary routine (version 2.0) --
	 * Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
	 * Courant Institute, Argonne National Lab, and Rice University
	 * September 30, 1994
	 *
	 * .. Scalar Arguments ..
	 */
	int32 LSAME_v;
	int32 INTA, 
	  INTB, 
	  ZCODE;

	DEBUG_ENTRY( "LSAME()" );
	/* ..
	 *
	 *  Purpose
	 *  =======
	 *
	 *  LSAME returns .true. if CA is the same letter as CB regardless of
	 *  case.
	 *
	 *  Arguments
	 *  =========
	 *
	 *  CA      (input) CHARACTER*1
	 *  CB      (input) CHARACTER*1
	 *  CA and CB specify the single characters to be compared.
	 *
	 * =====================================================================
	 *
	 * .. Intrinsic Functions .. */
	/* ..
	 * .. Local Scalars .. */
	/* ..
	 * .. Executable Statements ..
	 *
	 * Test if the characters are equal
	 * */
	LSAME_v = CA == CB;
	if( LSAME_v )
	{ 
		return LSAME_v;
	}

	/* Now test for equivalence if both characters are alphabetic.
	 * */
	ZCODE = 'Z';

	/* Use 'Z' rather than 'A' so that ASCII can be detected on Prime
	 * machines, on which ICHAR returns a value with bit 8 set.
	 * ICHAR('A') on Prime machines returns 193 which is the same as
	 * ICHAR('A') on an EBCDIC machine.
	 * */
	INTA = (CA);
	INTB = (CB);

	if( ZCODE == 90 || ZCODE == 122 )
	{

		/*  ASCII is assumed - ZCODE is the ASCII code of either lower or
		 *  upper case 'Z'.
		 * */
		if( INTA >= 97 && INTA <= 122 )
			INTA -= 32;
		if( INTB >= 97 && INTB <= 122 )
			INTB -= 32;

	}
	else if( ZCODE == 233 || ZCODE == 169 )
	{

		/*  EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
		 *  upper case 'Z'.
		 * */
		if( ((INTA >= 129 && INTA <= 137) || (INTA >= 145 && INTA <= 
		  153)) || (INTA >= 162 && INTA <= 169) )
			INTA += 64;
		if( ((INTB >= 129 && INTB <= 137) || (INTB >= 145 && INTB <= 
		  153)) || (INTB >= 162 && INTB <= 169) )
			INTB += 64;

	}
	else if( ZCODE == 218 || ZCODE == 250 )
	{

		/*  ASCII is assumed, on Prime machines - ZCODE is the ASCII code
		 *  plus 128 of either lower or upper case 'Z'.
		 * */
		if( INTA >= 225 && INTA <= 250 )
			INTA -= 32;
		if( INTB >= 225 && INTB <= 250 )
			INTB -= 32;
	}
	LSAME_v = INTA == INTB;

	/* RETURN
	 *
	 * End of LSAME
	 * */
	return LSAME_v;
}

/*IDAMAX lapack service routine */

STATIC int32 IDAMAX(int32 n, 
	  double dx[], 
	  int32 incx)
{
	/*
	 * finds the index of element having max. absolute value.
	 * jack dongarra, lapack, 3/11/78.
	 * modified 3/93 to return if incx .le. 0.
	 * modified 12/3/93, array(1) declarations changed to array(*)
	 *
	 */
	int32 IDAMAX_v, 
	  i, 
	  ix;
	double dmax;

	DEBUG_ENTRY( "IDAMAX()" );

	IDAMAX_v = 0;

	if( n < 1 || incx <= 0 )
	{ 
		return IDAMAX_v;
	}

	IDAMAX_v = 1;

	if( n == 1 )
	{ 
		return IDAMAX_v;
	}

	if( incx == 1 )
		goto L_20;

	/*  code for increment not equal to 1
	 * */
	ix = 1;
	dmax = fabs(dx[0]);
	ix += incx;
	for( i=2; i <= n; i++ )
	{
		/*  if(ABS(dx(ix)).le.dmax) go to 5 */
		if( fabs(dx[ix-1]) > dmax )
		{
			IDAMAX_v = i;
			dmax = fabs(dx[ix-1]);
		}
		ix += incx;
	}
	return IDAMAX_v;

	/*  code for increment equal to 1
	 * */
L_20:
	dmax = fabs(dx[0]);
	for( i=1; i < n; i++ )
	{
		/*  if(ABS(dx(i)).le.dmax) go to 30 */
		if( fabs(dx[i]) > dmax )
		{
			IDAMAX_v = i+1;
			dmax = fabs(dx[i]);
		}
	}
	return IDAMAX_v;
}

/*DTRSM lapack service routine */
STATIC void DTRSM(int32 SIDE, 
	  int32 UPLO, 
	  int32 TRANSA, 
	  int32 DIAG, 
	  int32 M, 
	  int32 N, 
	  double ALPHA, 
	  double *A, 
	  int32 LDA, 
	  double *B, 
	  int32 LDB)
{
	int32 LSIDE, 
	  NOUNIT, 
	  UPPER;
	int32 I, 
	  INFO, 
	  I_, 
	  J, 
	  J_, 
	  K, 
	  K_, 
	  NROWA;
	double TEMP;

	DEBUG_ENTRY( "DTRSM()" );
	/* .. Scalar Arguments .. */
	/* .. Array Arguments .. */
	/* ..
	 *
	 *  Purpose
	 *  =======
	 *
	 *  DTRSM  solves one of the matrix equations
	 *
	 * op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
	 *
	 *  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
	 *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
	 *
	 * op( A ) = A   or   op( A ) = A'.
	 *
	 *  The matrix X is overwritten on B.
	 *
	 *  Parameters
	 *  ==========
	 *
	 *  SIDE   - CHARACTER*1.
	 * On entry, SIDE specifies whether op( A ) appears on the left
	 * or right of X as follows:
	 *
	 *    SIDE = 'L' or 'l'   op( A )*X = alpha*B.
	 *
	 *    SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
	 *
	 * Unchanged on exit.
	 *
	 *  UPLO   - CHARACTER*1.
	 * On entry, UPLO specifies whether the matrix A is an upper or
	 * lower triangular matrix as follows:
	 *
	 *    UPLO = 'U' or 'u'   A is an upper triangular matrix.
	 *
	 *    UPLO = 'L' or 'l'   A is a lower triangular matrix.
	 *
	 * Unchanged on exit.
	 *
	 *  TRANSA - CHARACTER*1.
	 * On entry, TRANSA specifies the form of op( A ) to be used in
	 * the matrix multiplication as follows:
	 *
	 *    TRANSA = 'N' or 'n'   op( A ) = A.
	 *
	 *    TRANSA = 'T' or 't'   op( A ) = A'.
	 *
	 *    TRANSA = 'C' or 'c'   op( A ) = A'.
	 *
	 * Unchanged on exit.
	 *
	 *  DIAG   - CHARACTER*1.
	 * On entry, DIAG specifies whether or not A is unit triangular
	 * as follows:
	 *
	 *    DIAG = 'U' or 'u'   A is assumed to be unit triangular.
	 *
	 *    DIAG = 'N' or 'n'   A is not assumed to be unit
	 *                        triangular.
	 *
	 * Unchanged on exit.
	 *
	 *  M      - INTEGER.
	 * On entry, M specifies the number of rows of B. M must be at
	 * least zero.
	 * Unchanged on exit.
	 *
	 *  N      - INTEGER.
	 * On entry, N specifies the number of columns of B.  N must be
	 * at least zero.
	 * Unchanged on exit.
	 *
	 *  ALPHA  - DOUBLE PRECISION.
	 * On entry,  ALPHA specifies the scalar  alpha. When  alpha is
	 * zero then  A is not referenced and  B need not be set before
	 * entry.
	 * Unchanged on exit.
	 *
	 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
	 * when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
	 * Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
	 * upper triangular part of the array  A must contain the upper
	 * triangular matrix  and the strictly lower triangular part of
	 * A is not referenced.
	 * Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
	 * lower triangular part of the array  A must contain the lower
	 * triangular matrix  and the strictly upper triangular part of
	 * A is not referenced.
	 * Note that when  DIAG = 'U' or 'u',  the diagonal elements of
	 * A  are not referenced either,  but are assumed to be  unity.
	 * Unchanged on exit.
	 *
	 *  LDA    - INTEGER.
	 * On entry, LDA specifies the first dimension of A as declared
	 * in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
	 * LDA  must be at least  MAX( 1, m ),  when  SIDE = 'R' or 'r'
	 * then LDA must be at least MAX( 1, n ).
	 * Unchanged on exit.
	 *
	 *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
	 * Before entry,  the leading  m by n part of the array  B must
	 * contain  the  right-hand  side  matrix  B,  and  on exit  is
	 * overwritten by the solution matrix  X.
	 *
	 *  LDB    - INTEGER.
	 * On entry, LDB specifies the first dimension of B as declared
	 * in  the  calling  (sub)  program.   LDB  must  be  at  least
	 * MAX( 1, m ).
	 * Unchanged on exit.
	 *
	 *
	 *  Level 3 Blas routine.
	 *
	 *
	 *  -- Written on 8-February-1989.
	 * Jack Dongarra, Argonne National Laboratory.
	 * Iain Duff, AERE Harwell.
	 * Jeremy Du Croz, Numerical Algorithms Group Ltd.
	 * Sven Hammarling, Numerical Algorithms Group Ltd.
	 *
	 *
	 * .. External Functions .. */
	/* .. External Subroutines .. */
	/* .. Intrinsic Functions .. */
	/* .. Local Scalars .. */
	/* .. Parameters .. */
	/* ..
	 * .. Executable Statements ..
	 *
	 * Test the input parameters.
	 * */
	LSIDE = LSAME(SIDE,'L');
	if( LSIDE )
	{
		NROWA = M;
	}
	else
	{
		NROWA = N;
	}
	NOUNIT = LSAME(DIAG,'N');
	UPPER = LSAME(UPLO,'U');

	INFO = 0;
	if( (!LSIDE) && (!LSAME(SIDE,'R')) )
	{
		INFO = 1;
	}
	else if( (!UPPER) && (!LSAME(UPLO,'L')) )
	{
		INFO = 2;
	}
	else if( ((!LSAME(TRANSA,'N')) && (!LSAME(TRANSA,'T'))) && (!LSAME(TRANSA,
	  'C')) )
	{
		INFO = 3;
	}
	else if( (!LSAME(DIAG,'U')) && (!LSAME(DIAG,'N')) )
	{
		INFO = 4;
	}
	else if( M < 0 )
	{
		INFO = 5;
	}
	else if( N < 0 )
	{
		INFO = 6;
	}
	else if( LDA < MAX2(1,NROWA) )
	{
		INFO = 9;
	}
	else if( LDB < MAX2(1,M) )
	{
		INFO = 11;
	}
	if( INFO != 0 )
	{
		XERBLA("DTRSM ",INFO);
		/* XERBLA does not return */
	}

	/* Quick return if possible.
	 * */
	if( N == 0 )
		{ 
		return;}

	/* And when  alpha.eq.zero.
	 * */
	if( ALPHA == ZERO )
	{
		for( J=1; J <= N; J++ )
		{
			J_ = J - 1;
			for( I=1; I <= M; I++ )
			{
				I_ = I - 1;
				BB(J_,I_) = ZERO;
			}
		}
		return;
	}

	/* Start the operations.
	 * */
	if( LSIDE )
	{
		if( LSAME(TRANSA,'N') )
		{

			/* Form  B := alpha*inv( A )*B.
			 * */
			if( UPPER )
			{
				for( J=1; J <= N; J++ )
				{
					J_ = J - 1;
					if( ALPHA != ONE )
					{
						for( I=1; I <= M; I++ )
						{
							I_ = I - 1;
							BB(J_,I_) *= ALPHA;
						}
					}
					for( K=M; K >= 1; K-- )
					{
						K_ = K - 1;
						if( BB(J_,K_) != ZERO )
						{
							if( NOUNIT )
								BB(J_,K_) /= AA(K_,K_);
							TEMP = -BB(J_,K_);
							for( I=0; I < K-1; I++ )
							{
								BB(J_,I) += TEMP*AA(K_,I);
							}
						}
					}
				}
			}
			else
			{
				for( J=1; J <= N; J++ )
				{
					J_ = J - 1;
					if( ALPHA != ONE )
					{
						for( I=1; I <= M; I++ )
						{
							I_ = I - 1;
							BB(J_,I_) *= ALPHA;
						}
					}
					for( K=1; K <= M; K++ )
					{
						K_ = K - 1;
						if( BB(J_,K_) != ZERO )
						{
							if( NOUNIT )
								BB(J_,K_) /= AA(K_,K_);
							TEMP = -BB(J_,K_);
							for( I=K; I < M; I++ )
							{
								BB(J_,I) += TEMP*AA(K_,I);
							}
						}
					}
				}
			}
		}
		else
		{

			/* Form  B := alpha*inv( A' )*B.
			 * */
			if( UPPER )
			{
				for( J=1; J <= N; J++ )
				{
					J_ = J - 1;
					for( I=1; I <= M; I++ )
					{
						I_ = I - 1;
						TEMP = ALPHA*BB(J_,I_);
						for( K=1; K <= (I - 1); K++ )
						{
							K_ = K - 1;
							TEMP += -AA(I_,K_)*BB(J_,K_);
						}
						if( NOUNIT )
							TEMP /= AA(I_,I_);
						BB(J_,I_) = TEMP;
					}
				}
			}
			else
			{
				for( J=1; J <= N; J++ )
				{
					J_ = J - 1;
					for( I=M; I >= 1; I-- )
					{
						I_ = I - 1;
						TEMP = ALPHA*BB(J_,I_);
						for( K=I + 1; K <= M; K++ )
						{
							K_ = K - 1;
							TEMP += -AA(I_,K_)*BB(J_,K_);
						}
						if( NOUNIT )
							TEMP /= AA(I_,I_);
						BB(J_,I_) = TEMP;
					}
				}
			}
		}
	}
	else
	{
		if( LSAME(TRANSA,'N') )
		{

			/* Form  B := alpha*B*inv( A ).
			 * */
			if( UPPER )
			{
				for( J=1; J <= N; J++ )
				{
					J_ = J - 1;
					if( ALPHA != ONE )
					{
						for( I=1; I <= M; I++ )
						{
							I_ = I - 1;
							BB(J_,I_) *= ALPHA;
						}
					}
					for( K=1; K <= (J - 1); K++ )
					{
						K_ = K - 1;
						if( AA(J_,K_) != ZERO )
						{
							for( I=1; I <= M; I++ )
							{
								I_ = I - 1;
								BB(J_,I_) += -AA(J_,K_)*BB(K_,I_);
							}
						}
					}
					if( NOUNIT )
					{
						TEMP = ONE/AA(J_,J_);
						for( I=1; I <= M; I++ )
						{
							I_ = I - 1;
							BB(J_,I_) *= TEMP;
						}
					}
				}
			}
			else
			{
				for( J=N; J >= 1; J-- )
				{
					J_ = J - 1;
					if( ALPHA != ONE )
					{
						for( I=1; I <= M; I++ )
						{
							I_ = I - 1;
							BB(J_,I_) *= ALPHA;
						}
					}
					for( K=J + 1; K <= N; K++ )
					{
						K_ = K - 1;
						if( AA(J_,K_) != ZERO )
						{
							for( I=1; I <= M; I++ )
							{
								I_ = I - 1;
								BB(J_,I_) += -AA(J_,K_)*BB(K_,I_);
							}
						}
					}
					if( NOUNIT )
					{
						TEMP = ONE/AA(J_,J_);
						for( I=1; I <= M; I++ )
						{
							I_ = I - 1;
							BB(J_,I_) *= TEMP;
						}
					}
				}
			}
		}
		else
		{

			/* Form  B := alpha*B*inv( A' ).
			 * */
			if( UPPER )
			{
				for( K=N; K >= 1; K-- )
				{
					K_ = K - 1;
					if( NOUNIT )
					{
						TEMP = ONE/AA(K_,K_);
						for( I=1; I <= M; I++ )
						{
							I_ = I - 1;
							BB(K_,I_) *= TEMP;
						}
					}
					for( J=1; J <= (K - 1); J++ )
					{
						J_ = J - 1;
						if( AA(K_,J_) != ZERO )
						{
							TEMP = AA(K_,J_);
							for( I=1; I <= M; I++ )
							{
								I_ = I - 1;
								BB(J_,I_) += -TEMP*BB(K_,I_);
							}
						}
					}
					if( ALPHA != ONE )
					{
						for( I=1; I <= M; I++ )
						{
							I_ = I - 1;
							BB(K_,I_) *= ALPHA;
						}
					}
				}
			}
			else
			{
				for( K=1; K <= N; K++ )
				{
					K_ = K - 1;
					if( NOUNIT )
					{
						TEMP = ONE/AA(K_,K_);
						for( I=1; I <= M; I++ )
						{
							I_ = I - 1;
							BB(K_,I_) *= TEMP;
						}
					}
					for( J=K + 1; J <= N; J++ )
					{
						J_ = J - 1;
						if( AA(K_,J_) != ZERO )
						{
							TEMP = AA(K_,J_);
							for( I=1; I <= M; I++ )
							{
								I_ = I - 1;
								BB(J_,I_) += -TEMP*BB(K_,I_);
							}
						}
					}
					if( ALPHA != ONE )
					{
						for( I=1; I <= M; I++ )
						{
							I_ = I - 1;
							BB(K_,I_) *= ALPHA;
						}
					}
				}
			}
		}
	}

	return;

	/* End of DTRSM .
	 * */
#undef	B
#undef	A
}

/* ILAENV lapack helper routine */

/* >>chng 01 oct 22, remove two parameters since not used */
STATIC int32 ILAENV(int32 ISPEC, 
	  const char *NAME, 
	  /*char *OPTS, */
	  int32 N1, 
	  int32 N2, 
	  /*int32 N3, */
	  int32 N4)
{
	/*
	 *
	 *  -- LAPACK auxiliary routine (version 2.0) --
	 * Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
	 * Courant Institute, Argonne National Lab, and Rice University
	 * September 30, 1994
	 *
	 * .. Scalar Arguments ..
	 */
	char C2[3], 
	  C3[4], 
	  C4[3], 
	  SUBNAM[7];
	int32 CNAME, 
	  SNAME;
	char C1;
	int32 I, 
	  IC, 
	  ILAENV_v, 
	  IZ, 
	  NB, 
	  NBMIN, 
	  NX;

	DEBUG_ENTRY( "ILAENV()" );
	/* ..
	 *
	 *  Purpose
	 *  =======
	 *
	 *  ILAENV is called from the LAPACK routines to choose problem-dependent
	 *  parameters for the local environment.  See ISPEC for a description of
	 *  the parameters.
	 *
	 *  This version provides a set of parameters which should give good,
	 *  but not optimal, performance on many of the currently available
	 *  computers.  Users are encouraged to modify this subroutine to set
	 *  the tuning parameters for their particular machine using the option
	 *  and problem size information in the arguments.
	 *
	 *  This routine will not function correctly if it is converted to all
	 *  lower case.  Converting it to all upper case is allowed.
	 *
	 *  Arguments
	 *  =========
	 *
	 *  ISPEC   (input) INTEGER
	 *  Specifies the parameter to be returned as the value of
	 *  ILAENV.
	 *  = 1: the optimal blocksize; if this value is 1, an unblocked
	 *     algorithm will give the best performance.
	 *  = 2: the minimum block size for which the block routine
	 *     should be used; if the usable block size is less than
	 *     this value, an unblocked routine should be used.
	 *  = 3: the crossover point (in a block routine, for N less
	 *     than this value, an unblocked routine should be used)
	 *  = 4: the number of shifts, used in the nonsymmetric
	 *     eigenvalue routines
	 *  = 5: the minimum column dimension for blocking to be used;
	 *     rectangular blocks must have dimension at least k by m,
	 *     where k is given by ILAENV(2,...) and m by ILAENV(5,...)
	 *  = 6: the crossover point for the SVD (when reducing an m by n
	 *     matrix to bidiagonal form, if MAX(m,n)/MIN(m,n) exceeds
	 *     this value, a QR factorization is used first to reduce
	 *     the matrix to a triangular form.)
	 *  = 7: the number of processors
	 *  = 8: the crossover point for the multishift QR and QZ methods
	 *     for nonsymmetric eigenvalue problems.
	 *
	 *  NAME    (input) CHARACTER*(*)
	 *  The name of the calling subroutine, in either upper case or
	 *  lower case.
	 *
	 *  OPTS    (input) CHARACTER*(*)
	 *  The character options to the subroutine NAME, concatenated
	 *  into a single character string.  For example, UPLO = 'U',
	 *  TRANS = 'T', and DIAG = 'N' for a triangular routine would
	 *  be specified as OPTS = 'UTN'.
	 *
	 *  N1      (input) INTEGER
	 *  N2      (input) INTEGER
	 *  N3      (input) INTEGER
	 *  N4      (input) INTEGER
	 *  Problem dimensions for the subroutine NAME; these may not all
	 *  be required.
	 *
	 * (ILAENV) (output) INTEGER
	 *  >= 0: the value of the parameter specified by ISPEC
	 *  < 0:  if ILAENV = -k, the k-th argument had an illegal value.
	 *
	 *  Further Details
	 *  ===============
	 *
	 *  The following conventions have been used when calling ILAENV from the
	 *  LAPACK routines:
	 *  1)  OPTS is a concatenation of all of the character options to
	 *  subroutine NAME, in the same order that they appear in the
	 *  argument list for NAME, even if they are not used in determining
	 *  the value of the parameter specified by ISPEC.
	 *  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
	 *  that they appear in the argument list for NAME.  N1 is used
	 *  first, N2 second, and so on, and unused problem dimensions are
	 *  passed a value of -1.
	 *  3)  The parameter value returned by ILAENV is checked for validity in
	 *  the calling subroutine.  For example, ILAENV is used to retrieve
	 *  the optimal blocksize for STRTRI as follows:
	 *
	 *  NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
	 *  IF( NB.LE.1 ) NB = MAX( 1, N )
	 *
	 *  =====================================================================
	 *
	 * .. Local Scalars .. */
	/* ..
	 * .. Intrinsic Functions .. */
	/* ..
	 * .. Executable Statements ..
	 * */
	switch( ISPEC )
	{
		case 1: 
			{
				goto L_100;
			}
		case 2: 
			{
				goto L_100;
			}
		case 3: 
			{
				goto L_100;
			}
		case 4: 
			{
				goto L_400;
			}
		case 5: 
			{
				goto L_500;
			}
		case 6: 
			{
				goto L_600;
			}
		case 7: 
			{
				goto L_700;
			}
		case 8: 
			{
				goto L_800;
			}
		/* this is impossible, added by gjf to stop lint from complaining */
		default: 
		{
			/* Invalid value for ISPEC
			 * */
			ILAENV_v = -1;
			return ILAENV_v;
		}
	}

L_100:

	/* Convert NAME to upper case if the first character is lower case.
	 * */
	ILAENV_v = 1;
	strncpy( SUBNAM, NAME, 6 );
	SUBNAM[6] = '\0';
	IC = (SUBNAM[0]);
	IZ = 'Z';
	if( IZ == 90 || IZ == 122 )
	{

		/*  ASCII character set
		 * */
		if( IC >= 97 && IC <= 122 )
		{
			SUBNAM[0] = (char)(IC - 32);
			for( I=2; I <= 6; I++ )
			{
				IC = (SUBNAM[I-1]);
				if( IC >= 97 && IC <= 122 )
					SUBNAM[I - 1] = (char)(IC - 32);
			}
		}

	}
	else if( IZ == 233 || IZ == 169 )
	{

		/*  EBCDIC character set
		 * */
		if( ((IC >= 129 && IC <= 137) || (IC >= 145 && IC <= 153)) || 
		  (IC >= 162 && IC <= 169) )
		{
			SUBNAM[0] = (char)(IC + 64);
			for( I=2; I <= 6; I++ )
			{
				IC = (SUBNAM[I-1]);
				if( ((IC >= 129 && IC <= 137) || (IC >= 145 && IC <= 
				  153)) || (IC >= 162 && IC <= 169) )
					SUBNAM[I - 1] = (char)(IC + 64);
			}
		}

	}
	else if( IZ == 218 || IZ == 250 )
	{

		/*  Prime machines:  ASCII+128
		 * */
		if( IC >= 225 && IC <= 250 )
		{
			SUBNAM[0] = (char)(IC - 32);
			for( I=2; I <= 6; I++ )
			{
				IC = (SUBNAM[I-1]);
				if( IC >= 225 && IC <= 250 )
					SUBNAM[I - 1] = (char)(IC - 32);
			}
		}
	}

	C1 = SUBNAM[0];
	SNAME = C1 == 'S' || C1 == 'D';
	CNAME = C1 == 'C' || C1 == 'Z';
	if( !(CNAME || SNAME) )
	{ 
		return ILAENV_v;
	}

#	if 0
	strncpy( C2, SUBNAM+1, 2 );
	strncpy( C3, SUBNAM+3, 3 );
	strncpy( C4, C3+1, 2 );
#	endif

	/* >>chng 00 nov 08, from above to below, insure had run
	 * off end of string*/
	strncpy( C2, SUBNAM+1, 2 );
	C2[2] = '\0'; 
	strncpy( C3, SUBNAM+3, 3 );
	C3[3] = '\0'; 
	strncpy( C4, C3+1, 2 );
	C4[2] = '\0'; 

	switch( ISPEC )
	{
		case 1: goto L_110;
		case 2: goto L_200;
		case 3: goto L_300;
	}

L_110:

	/* ISPEC = 1:  block size
	 *
	 * In these examples, separate code is provided for setting NB for
	 * real and complex.  We assume that NB will take the same value in
	 * single or double precision.
	 * */
	NB = 1;

	if( strcmp(C2,"GE") == 0 )
	{
		if( strcmp(C3,"TRF") == 0 )
		{
			if( SNAME )
			{
				NB = 64;
			}
			else
			{
				NB = 64;
			}
		}
		else if( ((strcmp(C3,"QRF") == 0 || strcmp(C3,"RQF") == 0) || 
		  strcmp(C3,"LQF") == 0) || strcmp(C3,"QLF") == 0 )
		{
			if( SNAME )
			{
				NB = 32;
			}
			else
			{
				NB = 32;
			}
		}
		else if( strcmp(C3,"HRD") == 0 )
		{
			if( SNAME )
			{
				NB = 32;
			}
			else
			{
				NB = 32;
			}
		}
		else if( strcmp(C3,"BRD") == 0 )
		{
			if( SNAME )
			{
				NB = 32;
			}
			else
			{
				NB = 32;
			}
		}
		else if( strcmp(C3,"TRI") == 0 )
		{
			if( SNAME )
			{
				NB = 64;
			}
			else
			{
				NB = 64;
			}
		}
	}
	else if( strcmp(C2,"PO") == 0 )
	{
		if( strcmp(C3,"TRF") == 0 )
		{
			if( SNAME )
			{
				NB = 64;
			}
			else
			{
				NB = 64;
			}
		}
	}
	else if( strcmp(C2,"SY") == 0 )
	{
		if( strcmp(C3,"TRF") == 0 )
		{
			if( SNAME )
			{
				NB = 64;
			}
			else
			{
				NB = 64;
			}
		}
		else if( SNAME && strcmp(C3,"TRD") == 0 )
		{
			NB = 1;
		}
		else if( SNAME && strcmp(C3,"GST") == 0 )
		{
			NB = 64;
		}
	}
	else if( CNAME && strcmp(C2,"HE") == 0 )
	{
		if( strcmp(C3,"TRF") == 0 )
		{
			NB = 64;
		}
		else if( strcmp(C3,"TRD") == 0 )
		{
			NB = 1;
		}
		else if( strcmp(C3,"GST") == 0 )
		{
			NB = 64;
		}
	}
	else if( SNAME && strcmp(C2,"OR") == 0 )
	{
		if( C3[0] == 'G' )
		{
			if( (((((strcmp(C4,"QR") == 0 || strcmp(C4,"RQ") == 0) || 
			  strcmp(C4,"LQ") == 0) || strcmp(C4,"QL") == 0) || strcmp(C4
			  ,"HR") == 0) || strcmp(C4,"TR") == 0) || strcmp(C4,"BR") == 
			  0 )
			{
				NB = 32;
			}
		}
		else if( C3[0] == 'M' )
		{
			if( (((((strcmp(C4,"QR") == 0 || strcmp(C4,"RQ") == 0) || 
			  strcmp(C4,"LQ") == 0) || strcmp(C4,"QL") == 0) || strcmp(C4
			  ,"HR") == 0) || strcmp(C4,"TR") == 0) || strcmp(C4,"BR") == 
			  0 )
			{
				NB = 32;
			}
		}
	}
	else if( CNAME && strcmp(C2,"UN") == 0 )
	{
		if( C3[0] == 'G' )
		{
			if( (((((strcmp(C4,"QR") == 0 || strcmp(C4,"RQ") == 0) || 
			  strcmp(C4,"LQ") == 0) || strcmp(C4,"QL") == 0) || strcmp(C4
			  ,"HR") == 0) || strcmp(C4,"TR") == 0) || strcmp(C4,"BR") == 
			  0 )
			{
				NB = 32;
			}
		}
		else if( C3[0] == 'M' )
		{
			if( (((((strcmp(C4,"QR") == 0 || strcmp(C4,"RQ") == 0) || 
			  strcmp(C4,"LQ") == 0) || strcmp(C4,"QL") == 0) || strcmp(C4
			  ,"HR") == 0) || strcmp(C4,"TR") == 0) || strcmp(C4,"BR") == 
			  0 )
			{
				NB = 32;
			}
		}
	}
	else if( strcmp(C2,"GB") == 0 )
	{
		if( strcmp(C3,"TRF") == 0 )
		{
			if( SNAME )
			{
				if( N4 <= 64 )
				{
					NB = 1;
				}
				else
				{
					NB = 32;
				}
			}
			else
			{
				if( N4 <= 64 )
				{
					NB = 1;
				}
				else
				{
					NB = 32;
				}
			}
		}
	}
	else if( strcmp(C2,"PB") == 0 )
	{
		if( strcmp(C3,"TRF") == 0 )
		{
			if( SNAME )
			{
				if( N2 <= 64 )
				{
					NB = 1;
				}
				else
				{
					NB = 32;
				}
			}
			else
			{
				if( N2 <= 64 )
				{
					NB = 1;
				}
				else
				{
					NB = 32;
				}
			}
		}
	}
	else if( strcmp(C2,"TR") == 0 )
	{
		if( strcmp(C3,"TRI") == 0 )
		{
			if( SNAME )
			{
				NB = 64;
			}
			else
			{
				NB = 64;
			}
		}
	}
	else if( strcmp(C2,"LA") == 0 )
	{
		if( strcmp(C3,"UUM") == 0 )
		{
			if( SNAME )
			{
				NB = 64;
			}
			else
			{
				NB = 64;
			}
		}
	}
	else if( SNAME && strcmp(C2,"ST") == 0 )
	{
		if( strcmp(C3,"EBZ") == 0 )
		{
			NB = 1;
		}
	}
	ILAENV_v = NB;
	return ILAENV_v;

L_200:

	/* ISPEC = 2:  minimum block size
	 * */
	NBMIN = 2;
	if( strcmp(C2,"GE") == 0 )
	{
		if( ((strcmp(C3,"QRF") == 0 || strcmp(C3,"RQF") == 0) || strcmp(C3
		  ,"LQF") == 0) || strcmp(C3,"QLF") == 0 )
		{
			if( SNAME )
			{
				NBMIN = 2;
			}
			else
			{
				NBMIN = 2;
			}
		}
		else if( strcmp(C3,"HRD") == 0 )
		{
			if( SNAME )
			{
				NBMIN = 2;
			}
			else
			{
				NBMIN = 2;
			}
		}
		else if( strcmp(C3,"BRD") == 0 )
		{
			if( SNAME )
			{
				NBMIN = 2;
			}
			else
			{
				NBMIN = 2;
			}
		}
		else if( strcmp(C3,"TRI") == 0 )
		{
			if( SNAME )
			{
				NBMIN = 2;
			}
			else
			{
				NBMIN = 2;
			}
		}
	}
	else if( strcmp(C2,"SY") == 0 )
	{
		if( strcmp(C3,"TRF") == 0 )
		{
			if( SNAME )
			{
				NBMIN = 8;
			}
			else
			{
				NBMIN = 8;
			}
		}
		else if( SNAME && strcmp(C3,"TRD") == 0 )
		{
			NBMIN = 2;
		}
	}
	else if( CNAME && strcmp(C2,"HE") == 0 )
	{
		if( strcmp(C3,"TRD") == 0 )
		{
			NBMIN = 2;
		}
	}
	else if( SNAME && strcmp(C2,"OR") == 0 )
	{
		if( C3[0] == 'G' )
		{
			if( (((((strcmp(C4,"QR") == 0 || strcmp(C4,"RQ") == 0) || 
			  strcmp(C4,"LQ") == 0) || strcmp(C4,"QL") == 0) || strcmp(C4
			  ,"HR") == 0) || strcmp(C4,"TR") == 0) || strcmp(C4,"BR") == 
			  0 )
			{
				NBMIN = 2;
			}
		}
		else if( C3[0] == 'M' )
		{
			if( (((((strcmp(C4,"QR") == 0 || strcmp(C4,"RQ") == 0) || 
			  strcmp(C4,"LQ") == 0) || strcmp(C4,"QL") == 0) || strcmp(C4
			  ,"HR") == 0) || strcmp(C4,"TR") == 0) || strcmp(C4,"BR") == 
			  0 )
			{
				NBMIN = 2;
			}
		}
	}
	else if( CNAME && strcmp(C2,"UN") == 0 )
	{
		if( C3[0] == 'G' )
		{
			if( (((((strcmp(C4,"QR") == 0 || strcmp(C4,"RQ") == 0) || 
			  strcmp(C4,"LQ") == 0) || strcmp(C4,"QL") == 0) || strcmp(C4
			  ,"HR") == 0) || strcmp(C4,"TR") == 0) || strcmp(C4,"BR") == 
			  0 )
			{
				NBMIN = 2;
			}
		}
		else if( C3[0] == 'M' )
		{
			if( (((((strcmp(C4,"QR") == 0 || strcmp(C4,"RQ") == 0) || 
			  strcmp(C4,"LQ") == 0) || strcmp(C4,"QL") == 0) || strcmp(C4
			  ,"HR") == 0) || strcmp(C4,"TR") == 0) || strcmp(C4,"BR") == 
			  0 )
			{
				NBMIN = 2;
			}
		}
	}
	ILAENV_v = NBMIN;
	return ILAENV_v;

L_300:

	/* ISPEC = 3:  crossover point
	 * */
	NX = 0;
	if( strcmp(C2,"GE") == 0 )
	{
		if( ((strcmp(C3,"QRF") == 0 || strcmp(C3,"RQF") == 0) || strcmp(C3
		  ,"LQF") == 0) || strcmp(C3,"QLF") == 0 )
		{
			if( SNAME )
			{
				NX = 128;
			}
			else
			{
				NX = 128;
			}
		}
		else if( strcmp(C3,"HRD") == 0 )
		{
			if( SNAME )
			{
				NX = 128;
			}
			else
			{
				NX = 128;
			}
		}
		else if( strcmp(C3,"BRD") == 0 )
		{
			if( SNAME )
			{
				NX = 128;
			}
			else
			{
				NX = 128;
			}
		}
	}
	else if( strcmp(C2,"SY") == 0 )
	{
		if( SNAME && strcmp(C3,"TRD") == 0 )
		{
			NX = 1;
		}
	}
	else if( CNAME && strcmp(C2,"HE") == 0 )
	{
		if( strcmp(C3,"TRD") == 0 )
		{
			NX = 1;
		}
	}
	else if( SNAME && strcmp(C2,"OR") == 0 )
	{
		if( C3[0] == 'G' )
		{
			if( (((((strcmp(C4,"QR") == 0 || strcmp(C4,"RQ") == 0) || 
			  strcmp(C4,"LQ") == 0) || strcmp(C4,"QL") == 0) || strcmp(C4
			  ,"HR") == 0) || strcmp(C4,"TR") == 0) || strcmp(C4,"BR") == 
			  0 )
			{
				NX = 128;
			}
		}
	}
	else if( CNAME && strcmp(C2,"UN") == 0 )
	{
		if( C3[0] == 'G' )
		{
			if( (((((strcmp(C4,"QR") == 0 || strcmp(C4,"RQ") == 0) || 
			  strcmp(C4,"LQ") == 0) || strcmp(C4,"QL") == 0) || strcmp(C4
			  ,"HR") == 0) || strcmp(C4,"TR") == 0) || strcmp(C4,"BR") == 
			  0 )
			{
				NX = 128;
			}
		}
	}
	ILAENV_v = NX;
	return ILAENV_v;

L_400:

	/* ISPEC = 4:  number of shifts (used by xHSEQR)
	 * */
	ILAENV_v = 6;
	return ILAENV_v;

L_500:

	/* ISPEC = 5:  minimum column dimension (not used)
	 * */
	ILAENV_v = 2;
	return ILAENV_v;

L_600:

	/* ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
	 * */
	ILAENV_v = (int32)((realnum)(MIN2(N1,N2))*1.6e0);
	return ILAENV_v;

L_700:

	/* ISPEC = 7:  number of processors (not used)
	 * */
	ILAENV_v = 1;
	return ILAENV_v;

L_800:

	/* ISPEC = 8:  crossover point for multishift (used by xHSEQR)
	 * */
	ILAENV_v = 50;
	return ILAENV_v;

	/* End of ILAENV
	 * */
}

/*DSWAP lapack routine */

STATIC void DSWAP(int32 n, 
	  double dx[], 
	  int32 incx, 
	  double dy[], 
	  int32 incy)
{
	int32 i, 
	  ix, 
	  iy, 
	  m;
	double dtemp;

	DEBUG_ENTRY( "DSWAP()" );

	/* interchanges two vectors.
	 * uses unrolled loops for increments equal one.
	 * jack dongarra, lapack, 3/11/78.
	 * modified 12/3/93, array(1) declarations changed to array(*)
	 * */

	if( n <= 0 )
		{ 
		return;}
	if( incx == 1 && incy == 1 )
		goto L_20;

	/* code for unequal increments or equal increments not equal
	 * to 1
	 * */
	ix = 1;
	iy = 1;

	if( incx < 0 )
		ix = (-n + 1)*incx + 1;

	if( incy < 0 )
		iy = (-n + 1)*incy + 1;

	for( i=0; i < n; i++ )
	{
		dtemp = dx[ix-1];
		dx[ix-1] = dy[iy-1];
		dy[iy-1] = dtemp;
		ix += incx;
		iy += incy;
	}
	return;

	/* code for both increments equal to 1
	 *
	 *
	 * clean-up loop
	 * */
L_20:
	m = n%3;
	if( m == 0 )
		goto L_40;

	for( i=0; i < m; i++ )
	{
		dtemp = dx[i];
		dx[i] = dy[i];
		dy[i] = dtemp;
	}

	if( n < 3 )
	{ 
		return;
	}

L_40:
	for( i=m; i < n; i += 3 )
	{
		dtemp = dx[i];
		dx[i] = dy[i];
		dy[i] = dtemp;
		dtemp = dx[i+1];
		dx[i+1] = dy[i+1];
		dy[i+1] = dtemp;
		dtemp = dx[i+2];
		dx[i+2] = dy[i+2];
		dy[i+2] = dtemp;
	}
	return;
}

/*DSCAL lapack routine */

STATIC void DSCAL(int32 n, 
	  double da, 
	  double dx[], 
	  int32 incx)
{
	int32 i, 
	  nincx;

	DEBUG_ENTRY( "DSCAL()" );

	/* scales a vector by a constant.
	 * uses unrolled loops for increment equal to one.
	 * jack dongarra, lapack, 3/11/78.
	 * modified 3/93 to return if incx .le. 0.
	 * modified 12/3/93, array(1) declarations changed to array(*)
	 * */

	if( n <= 0 || incx <= 0 )
		{ 
		return;}
	if( incx == 1 )
		goto L_20;

	/*  code for increment not equal to 1
	 * */
	nincx = n*incx;
	/*for( i=1, _do0=DOCNT(i,nincx,_do1 = incx); _do0 > 0; i += _do1, _do0-- )*/
	for( i=0; i<nincx; i = i + incx)
	{
		dx[i] *= da;
	}
	return;

	/*  code for increment equal to 1
	 *
	 *
	 *  clean-up loop
	 * */
L_20:
	for( i=0; i < n; i++ )
	{
		dx[i] *= da;
	}
	return;
}

/*DLASWP  -- LAPACK auxiliary routine (version 2.0) --*/

STATIC void DLASWP(int32 N, 
	  double *A, 
	  int32 LDA, 
	  int32 K1, 
	  int32 K2, 
	  int32 IPIV[], 
	  int32 INCX)
{
	int32 I, 
	  IP, 
	  IX, 
	  I_;

	DEBUG_ENTRY( "DLASWP()" );

	/*  -- LAPACK auxiliary routine (version 2.0) --
	 * Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
	 * Courant Institute, Argonne National Lab, and Rice University
	 * October 31, 1992
	 *
	 * .. Scalar Arguments .. */
	/* ..
	 * .. Array Arguments .. */
	/* ..
	 *
	 *  Purpose
	 *  =======
	 *
	 *  DLASWP performs a series of row interchanges on the matrix A.
	 *  One row interchange is initiated for each of rows K1 through K2 of A.
	 *
	 *  Arguments
	 *  =========
	 *
	 *  N       (input) INTEGER
	 *  The number of columns of the matrix A.
	 *
	 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
	 *  On entry, the matrix of column dimension N to which the row
	 *  interchanges will be applied.
	 *  On exit, the permuted matrix.
	 *
	 *  LDA     (input) INTEGER
	 *  The leading dimension of the array A.
	 *
	 *  K1      (input) INTEGER
	 *  The first element of IPIV for which a row interchange will
	 *  be done.
	 *
	 *  K2      (input) INTEGER
	 *  The last element of IPIV for which a row interchange will
	 *  be done.
	 *
	 *  IPIV    (input) INTEGER array, dimension (M*ABS(INCX))
	 *  The vector of pivot indices.  Only the elements in positions
	 *  K1 through K2 of IPIV are accessed.
	 *  IPIV(K) = L implies rows K and L are to be interchanged.
	 *
	 *  INCX    (input) INTEGER
	 *  The increment between successive values of IPIV.  If IPIV
	 *  is negative, the pivots are applied in reverse order.
	 *
	 * =====================================================================
	 *
	 * .. Local Scalars .. */
	/* ..
	 * .. External Subroutines .. */
	/* ..
	 * .. Executable Statements ..
	 *
	 * Interchange row I with row IPIV(I) for each of rows K1 through K2.
	 * */
	if( INCX == 0 )
		{ 
		return;}
	if( INCX > 0 )
	{
		IX = K1;
	}
	else
	{
		IX = 1 + (1 - K2)*INCX;
	}
	if( INCX == 1 )
	{
		for( I=K1; I <= K2; I++ )
		{
			I_ = I - 1;
			IP = IPIV[I_];
			if( IP != I )
				DSWAP(N,&AA(0,I_),LDA,&AA(0,IP-1),LDA);
		}
	}
	else if( INCX > 1 )
	{
		for( I=K1; I <= K2; I++ )
		{
			I_ = I - 1;
			IP = IPIV[IX-1];
			if( IP != I )
				DSWAP(N,&AA(0,I_),LDA,&AA(0,IP-1),LDA);
			IX += INCX;
		}
	}
	else if( INCX < 0 )
	{
		for( I=K2; I >= K1; I-- )
		{
			I_ = I - 1;
			IP = IPIV[IX-1];
			if( IP != I )
				DSWAP(N,&AA(0,I_),LDA,&AA(0,IP-1),LDA);
			IX += INCX;
		}
	}

	return;

	/* End of DLASWP
	 * */
#undef	A
}

/*DGETF2 lapack service routine */

STATIC void DGETF2(int32 M, 
	  int32 N, 
	  double *A, 
	  int32 LDA, 
	  int32 IPIV[], 
	  int32 *INFO)
{
	int32 J, 
	  JP, 
	  J_, 
	  limit;

	DEBUG_ENTRY( "DGETF2()" );

	/*  -- LAPACK routine (version 2.0) --
	 * Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
	 * Courant Institute, Argonne National Lab, and Rice University
	 * June 30, 1992
	 *
	 * .. Scalar Arguments .. */
	/* ..
	 * .. Array Arguments .. */
	/* ..
	 *
	 *  Purpose
	 *  =======
	 *
	 *  DGETF2 computes an LU factorization of a general m-by-n matrix A
	 *  using partial pivoting with row interchanges.
	 *
	 *  The factorization has the form
	 * A = P * L * U
	 *  where P is a permutation matrix, L is lower triangular with unit
	 *  diagonal elements (lower trapezoidal if m > n), and U is upper
	 *  triangular (upper trapezoidal if m < n).
	 *
	 *  This is the right-looking Level 2 BLAS version of the algorithm.
	 *
	 *  Arguments
	 *  =========
	 *
	 *  M       (input) INTEGER
	 *  The number of rows of the matrix A.  M >= 0.
	 *
	 *  N       (input) INTEGER
	 *  The number of columns of the matrix A.  N >= 0.
	 *
	 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
	 *  On entry, the m by n matrix to be factored.
	 *  On exit, the factors L and U from the factorization
	 *  A = P*L*U; the unit diagonal elements of L are not stored.
	 *
	 *  LDA     (input) INTEGER
	 *  The leading dimension of the array A.  LDA >= MAX(1,M).
	 *
	 *  IPIV    (output) INTEGER array, dimension (MIN(M,N))
	 *  The pivot indices; for 1 <= i <= MIN(M,N), row i of the
	 *  matrix was interchanged with row IPIV(i).
	 *
	 *  INFO    (output) INTEGER
	 *  = 0: successful exit
	 *  < 0: if INFO = -k, the k-th argument had an illegal value
	 *  > 0: if INFO = k, U(k,k) is exactly zero. The factorization
	 *     has been completed, but the factor U is exactly
	 *     singular, and division by zero will occur if it is used
	 *     to solve a system of equations.
	 *
	 *  =====================================================================
	 *
	 * .. Parameters .. */
	/* ..
	 * .. Local Scalars .. */
	/* ..
	 * .. External Functions .. */
	/* ..
	 * .. External Subroutines .. */
	/* ..
	 * .. Intrinsic Functions .. */
	/* ..
	 * .. Executable Statements ..
	 *
	 * Test the input parameters.
	 * */
	*INFO = 0;
	if( M < 0 )
	{
		*INFO = -1;
	}
	else if( N < 0 )
	{
		*INFO = -2;
	}
	else if( LDA < MAX2(1,M) )
	{
		*INFO = -4;
	}
	if( *INFO != 0 )
	{
		XERBLA("DGETF2",-*INFO);
		/* XERBLA does not return */
	}

	/* Quick return if possible
	 * */
	if( M == 0 || N == 0 )
		{ 
		return;}

	limit = MIN2(M,N);
	for( J=1; J <= limit; J++ )
	{
		J_ = J - 1;

		/*  Find pivot and test for singularity.
		 * */
		JP = J - 1 + IDAMAX(M-J+1,&AA(J_,J_),1);
		IPIV[J_] = JP;
		if( AA(J_,JP-1) != ZERO )
		{
			/* Apply the interchange to columns 1:N.
			 * */
			if( JP != J )
				DSWAP(N,&AA(0,J_),LDA,&AA(0,JP-1),LDA);

			/* Compute elements J+1:M of J-th column.
			 * */
			if( J < M )
				DSCAL(M-J,ONE/AA(J_,J_),&AA(J_,J_+1),1);
		}
		else if( *INFO == 0 )
		{
			*INFO = J;
		}

		if( J < MIN2(M,N) )
		{
			/* Update trailing submatrix.
			 * */
			DGER(M-J,N-J,-ONE,&AA(J_,J_+1),1,&AA(J_+1,J_),LDA,&AA(J_+1,J_+1),
			  LDA);
		}
	}
	return;

	/* End of DGETF2
	 * */
#undef	A
}

/*DGER service routine for matrix inversion */

STATIC void DGER(int32 M, 
	  int32 N, 
	  double ALPHA, 
	  double X[], 
	  int32 INCX, 
	  double Y[], 
	  int32 INCY, 
	  double *A, 
	  int32 LDA)
{
	int32 I, 
	  INFO, 
	  IX, 
	  I_, 
	  J, 
	  JY, 
	  J_, 
	  KX;
	double TEMP;

	DEBUG_ENTRY( "DGER()" );
	/* .. Scalar Arguments .. */
	/* .. Array Arguments .. */
	/* ..
	 *
	 *  Purpose
	 *  =======
	 *
	 *  DGER   performs the rank 1 operation
	 *
	 * A := alpha*x*y' + A,
	 *
	 *  where alpha is a scalar, x is an m element vector, y is an n element
	 *  vector and A is an m by n matrix.
	 *
	 *  Parameters
	 *  ==========
	 *
	 *  M      - INTEGER.
	 * On entry, M specifies the number of rows of the matrix A.
	 * M must be at least zero.
	 * Unchanged on exit.
	 *
	 *  N      - INTEGER.
	 * On entry, N specifies the number of columns of the matrix A.
	 * N must be at least zero.
	 * Unchanged on exit.
	 *
	 *  ALPHA  - DOUBLE PRECISION.
	 * On entry, ALPHA specifies the scalar alpha.
	 * Unchanged on exit.
	 *
	 *  X      - DOUBLE PRECISION array of dimension at least
	 * ( 1 + ( m - 1 )*ABS( INCX ) ).
	 * Before entry, the incremented array X must contain the m
	 * element vector x.
	 * Unchanged on exit.
	 *
	 *  INCX   - INTEGER.
	 * On entry, INCX specifies the increment for the elements of
	 * X. INCX must not be zero.
	 * Unchanged on exit.
	 *
	 *  Y      - DOUBLE PRECISION array of dimension at least
	 * ( 1 + ( n - 1 )*ABS( INCY ) ).
	 * Before entry, the incremented array Y must contain the n
	 * element vector y.
	 * Unchanged on exit.
	 *
	 *  INCY   - INTEGER.
	 * On entry, INCY specifies the increment for the elements of
	 * Y. INCY must not be zero.
	 * Unchanged on exit.
	 *
	 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
	 * Before entry, the leading m by n part of the array A must
	 * contain the matrix of coefficients. On exit, A is
	 * overwritten by the updated matrix.
	 *
	 *  LDA    - INTEGER.
	 * On entry, LDA specifies the first dimension of A as declared
	 * in the calling (sub) program. LDA must be at least
	 * MAX( 1, m ).
	 * Unchanged on exit.
	 *
	 *
	 *  Level 2 Blas routine.
	 *
	 *  -- Written on 22-October-1986.
	 * Jack Dongarra, Argonne National Lab.
	 * Jeremy Du Croz, Nag Central Office.
	 * Sven Hammarling, Nag Central Office.
	 * Richard Hanson, Sandia National Labs.
	 *
	 *
	 * .. Parameters .. */
	/* .. Local Scalars .. */
	/* .. External Subroutines .. */
	/* .. Intrinsic Functions .. */
	/* ..
	 * .. Executable Statements ..
	 *
	 * Test the input parameters.
	 * */
	INFO = 0;
	if( M < 0 )
	{
		INFO = 1;
	}
	else if( N < 0 )
	{
		INFO = 2;
	}
	else if( INCX == 0 )
	{
		INFO = 5;
	}
	else if( INCY == 0 )
	{
		INFO = 7;
	}
	else if( LDA < MAX2(1,M) )
	{
		INFO = 9;
	}
	if( INFO != 0 )
	{
		XERBLA("DGER  ",INFO);
		/* XERBLA does not return */
	}

	/* Quick return if possible.
	 * */
	if( ((M == 0) || (N == 0)) || (ALPHA == ZERO) )
		{ 
		return;}

	/* Start the operations. In this version the elements of A are
	 * accessed sequentially with one pass through A.
	 * */
	if( INCY > 0 )
	{
		JY = 1;
	}
	else
	{
		JY = 1 - (N - 1)*INCY;
	}
	if( INCX == 1 )
	{
		for( J=1; J <= N; J++ )
		{
			J_ = J - 1;
			if( Y[JY-1] != ZERO )
			{
				TEMP = ALPHA*Y[JY-1];
				for( I=0; I < M; I++ )
				{
					AA(J_,I) += X[I]*TEMP;
				}
			}
			JY += INCY;
		}
	}
	else
	{
		if( INCX > 0 )
		{
			KX = 1;
		}
		else
		{
			KX = 1 - (M - 1)*INCX;
		}
		for( J=1; J <= N; J++ )
		{
			J_ = J - 1;
			if( Y[JY-1] != ZERO )
			{
				TEMP = ALPHA*Y[JY-1];
				IX = KX;
				for( I=1; I <= M; I++ )
				{
					I_ = I - 1;
					AA(J_,I_) += X[IX-1]*TEMP;
					IX += INCX;
				}
			}
			JY += INCY;
		}
	}

	return;

	/* End of DGER  .
	 * */
#undef	A
}

/* DGEMM matrix inversion helper routine*/

STATIC void DGEMM(int32 TRANSA, 
	  int32 TRANSB, 
	  int32 M, 
	  int32 N, 
	  int32 K, 
	  double ALPHA, 
	  double *A, 
	  int32 LDA, 
	  double *B, 
	  int32 LDB, 
	  double BETA, 
	  double *C, 
	  int32 LDC)
{
	int32 NOTA, 
	  NOTB;
	int32 I, 
	  INFO, 
	  J, 
	  L, 
	  NROWA, 
	  NROWB;
	double TEMP;

	DEBUG_ENTRY( "DGEMM()" );
	/* .. Scalar Arguments .. */
	/* .. Array Arguments .. */
	/* ..
	 *
	 *  Purpose
	 *  =======
	 *
	 *  DGEMM  performs one of the matrix-matrix operations
	 *
	 * C := alpha*op( A )*op( B ) + beta*C,
	 *
	 *  where  op( X ) is one of
	 *
	 * op( X ) = X   or   op( X ) = X',
	 *
	 *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
	 *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
	 *
	 *  Parameters
	 *  ==========
	 *
	 *  TRANSA - CHARACTER*1.
	 * On entry, TRANSA specifies the form of op( A ) to be used in
	 * the matrix multiplication as follows:
	 *
	 *    TRANSA = 'N' or 'n',  op( A ) = A.
	 *
	 *    TRANSA = 'T' or 't',  op( A ) = A'.
	 *
	 *    TRANSA = 'C' or 'c',  op( A ) = A'.
	 *
	 * Unchanged on exit.
	 *
	 *  TRANSB - CHARACTER*1.
	 * On entry, TRANSB specifies the form of op( B ) to be used in
	 * the matrix multiplication as follows:
	 *
	 *    TRANSB = 'N' or 'n',  op( B ) = B.
	 *
	 *    TRANSB = 'T' or 't',  op( B ) = B'.
	 *
	 *    TRANSB = 'C' or 'c',  op( B ) = B'.
	 *
	 * Unchanged on exit.
	 *
	 *  M      - INTEGER.
	 * On entry,  M  specifies  the number  of rows  of the  matrix
	 * op( A )  and of the  matrix  C.  M  must  be at least  zero.
	 * Unchanged on exit.
	 *
	 *  N      - INTEGER.
	 * On entry,  N  specifies the number  of columns of the matrix
	 * op( B ) and the number of columns of the matrix C. N must be
	 * at least zero.
	 * Unchanged on exit.
	 *
	 *  K      - INTEGER.
	 * On entry,  K  specifies  the number of columns of the matrix
	 * op( A ) and the number of rows of the matrix op( B ). K must
	 * be at least  zero.
	 * Unchanged on exit.
	 *
	 *  ALPHA  - DOUBLE PRECISION.
	 * On entry, ALPHA specifies the scalar alpha.
	 * Unchanged on exit.
	 *
	 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
	 * k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
	 * Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
	 * part of the array  A  must contain the matrix  A,  otherwise
	 * the leading  k by m  part of the array  A  must contain  the
	 * matrix A.
	 * Unchanged on exit.
	 *
	 *  LDA    - INTEGER.
	 * On entry, LDA specifies the first dimension of A as declared
	 * in the calling (sub) program. When  TRANSA = 'N' or 'n' then
	 * LDA must be at least  MAX( 1, m ), otherwise  LDA must be at
	 * least  MAX( 1, k ).
	 * Unchanged on exit.
	 *
	 *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
	 * n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
	 * Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
	 * part of the array  B  must contain the matrix  B,  otherwise
	 * the leading  n by k  part of the array  B  must contain  the
	 * matrix B.
	 * Unchanged on exit.
	 *
	 *  LDB    - INTEGER.
	 * On entry, LDB specifies the first dimension of B as declared
	 * in the calling (sub) program. When  TRANSB = 'N' or 'n' then
	 * LDB must be at least  MAX( 1, k ), otherwise  LDB must be at
	 * least  MAX( 1, n ).
	 * Unchanged on exit.
	 *
	 *  BETA   - DOUBLE PRECISION.
	 * On entry,  BETA  specifies the scalar  beta.  When  BETA  is
	 * supplied as zero then C need not be set on input.
	 * Unchanged on exit.
	 *
	 *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
	 * Before entry, the leading  m by n  part of the array  C must
	 * contain the matrix  C,  except when  beta  is zero, in which
	 * case C need not be set on entry.
	 * On exit, the array  C  is overwritten by the  m by n  matrix
	 * ( alpha*op( A )*op( B ) + beta*C ).
	 *
	 *  LDC    - INTEGER.
	 * On entry, LDC specifies the first dimension of C as declared
	 * in  the  calling  (sub)  program.   LDC  must  be  at  least
	 * MAX( 1, m ).
	 * Unchanged on exit.
	 *
	 *
	 *  Level 3 Blas routine.
	 *
	 *  -- Written on 8-February-1989.
	 * Jack Dongarra, Argonne National Laboratory.
	 * Iain Duff, AERE Harwell.
	 * Jeremy Du Croz, Numerical Algorithms Group Ltd.
	 * Sven Hammarling, Numerical Algorithms Group Ltd.
	 *
	 *
	 * .. External Functions .. */
	/* .. External Subroutines .. */
	/* .. Intrinsic Functions .. */
	/* .. Local Scalars .. */
	/* .. Parameters .. */
	/* ..
	 * .. Executable Statements ..
	 *
	 * Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
	 * transposed and set  NROWA, NROWB  as the number of rows
	 * and  columns of  A  and the  number of  rows  of  B  respectively.
	 * */
	NOTA = LSAME(TRANSA,'N');
	NOTB = LSAME(TRANSB,'N');

	if( NOTA )
	{
		NROWA = M;
	}
	else
	{
		NROWA = K;
	}

	if( NOTB )
	{
		NROWB = K;
	}
	else
	{
		NROWB = N;
	}

	/* Test the input parameters.
	 * */
	INFO = 0;
	if( ((!NOTA) && (!LSAME(TRANSA,'C'))) && (!LSAME(TRANSA,'T')) )
	{
		INFO = 1;
	}
	else if( 
		((!NOTB) && (!LSAME(TRANSB,'C'))) && (!LSAME(TRANSB,'T')) )
	{
		INFO = 2;
	}

	else if( M < 0 )
	{
		INFO = 3;
	}

	else if( N < 0 )
	{
		INFO = 4;
	}

	else if( K < 0 )
	{
		INFO = 5;
	}

	else if( LDA < MAX2(1,NROWA) )
	{
		INFO = 8;
	}

	else if( LDB < MAX2(1,NROWB) )
	{
		INFO = 10;
	}

	else if( LDC < MAX2(1,M) )
	{
		INFO = 13;
	}

	if( INFO != 0 )
	{
		XERBLA("DGEMM ",INFO);
		/* XERBLA does not return */
	}

	/* Quick return if possible.
	 * */
	if( ((M == 0) || (N == 0)) || (((ALPHA == ZERO) || (K == 0)) && 
	  (BETA == ONE)) )
	{ 
		return;
	}

	/* And if  alpha.eq.zero. */
	if( ALPHA == ZERO )
	{
		if( BETA == ZERO )
		{
			for( J=0; J < N; J++ )
			{
				for( I=0; I < M; I++ )
				{
					CC(J,I) = ZERO;
				}
			}
		}

		else
		{
			for( J=0; J < N; J++ )
			{
				for( I=0; I < M; I++ )
				{
					CC(J,I) *= BETA;
				}
			}
		}
		return;
	}

	/* Start the operations.*/
	if( NOTB )
	{

		if( NOTA )
		{

			/* Form  C := alpha*A*B + beta*C.
			 * */
			for( J=0; J < N; J++ )
			{
				if( BETA == ZERO )
				{
					for( I=0; I < M; I++ )
					{
						CC(J,I) = ZERO;
					}
				}

				else if( BETA != ONE )
				{
					for( I=0; I < M; I++ )
					{
						CC(J,I) *= BETA;
					}
				}

				for( L=0; L < K; L++ )
				{
					if( BB(J,L) != ZERO )
					{
						TEMP = ALPHA*BB(J,L);
						for( I=0; I < M; I++ )
						{
							CC(J,I) += TEMP*AA(L,I);
						}
					}
				}
			}
		}
		else
		{

			/* Form  C := alpha*A'*B + beta*C */
			for( J=0; J < N; J++ )
			{
				for( I=0; I < M; I++ )
				{
					TEMP = ZERO;
					for( L=0; L < K; L++ )
					{
						TEMP += AA(I,L)*BB(J,L);
					}

					if( BETA == ZERO )
					{
						CC(J,I) = ALPHA*TEMP;
					}
					else
					{
						CC(J,I) = ALPHA*TEMP + BETA*CC(J,I);
					}
				}
			}
		}
	}
	else
	{
		if( NOTA )
		{

			/* Form  C := alpha*A*B' + beta*C
			 * */
			for( J=0; J < N; J++ )
			{
				if( BETA == ZERO )
				{
					for( I=0; I < M; I++ )
					{
						CC(J,I) = ZERO;
					}
				}

				else if( BETA != ONE )
				{
					for( I=0; I < M; I++ )
					{
						CC(J,I) *= BETA;
					}
				}

				for( L=0; L < K; L++ )
				{
					if( BB(L,J) != ZERO )
					{
						TEMP = ALPHA*BB(L,J);
						for( I=0; I < M; I++ )
						{
							CC(J,I) += TEMP*AA(L,I);
						}
					}
				}
			}
		}

		else
		{

			/* Form  C := alpha*A'*B' + beta*C */
			for( J=0; J < N; J++ )
			{

				for( I=0; I < M; I++ )
				{
					TEMP = ZERO;

					for( L=0; L < K; L++ )
					{
						TEMP += AA(I,L)*BB(L,J);
					}

					if( BETA == ZERO )
					{
						CC(J,I) = ALPHA*TEMP;
					}
					else
					{
						CC(J,I) = ALPHA*TEMP + BETA*CC(J,I);
					}

				}
			}
		}
	}

	return;

	/* End of DGEMM .*/
#undef	C
#undef	B
#undef	A
}
#undef AA
#undef BB
#undef CC

/* Subroutine */ 
#if 0
STATIC int32 DGTSV(int32 *n, int32 *nrhs, double *dl, 
	double *d__, double *du, double *b, int32 *ldb, int32 
	*info)
{
    /* System generated locals */
    int32 b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    double fact, temp;
    int32 i__, j;


#define b_ref(a_1,a_2) b[(a_2)*(b_dim1) + (a_1)]


/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1999   


    Purpose   
    =======   

    DGTSV  solves the equation   

       A*X = B,   

    where A is an n by n tridiagonal matrix, by Gaussian elimination with   
    partial pivoting.   

    Note that the equation  A'*X = B  may be solved by interchanging the   
    order of the arguments DU and DL.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    DL      (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, DL must contain the (n-1) sub-diagonal elements of   
            A.   

            On exit, DL is overwritten by the (n-2) elements of the   
            second super-diagonal of the upper triangular matrix U from   
            the LU factorization of A, in DL(1), ..., DL(n-2).   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, D must contain the diagonal elements of A.   

            On exit, D is overwritten by the n diagonal elements of U.   

    DU      (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, DU must contain the (n-1) super-diagonal elements   
            of A.   

            On exit, DU is overwritten by the (n-1) elements of the first   
            super-diagonal of U.   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N by NRHS matrix of right hand side matrix B.   
            On exit, if INFO = 0, the N by NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, U(i,i) is exactly zero, and the solution   
                 has not been computed.  The factorization has not been   
                 completed unless i = N.   

    =====================================================================   


       Parameter adjustments */
    --dl;
    --d__;
    --du;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    if(*n < 0) {
	*info = -1;
    } else if(*nrhs < 0) {
	*info = -2;
    } else if(*ldb < *n && *ldb < 1) {
	*info = -7;
    }
    if(*info != 0) {
	i__1 = -(*info);
	XERBLA("DGTSV ", i__1);
	/* XERBLA does not return */
    }

    if(*n == 0) {
	return 0;
    }

    if(*nrhs == 1) {
	i__1 = *n - 2;
	for(i__ = 1; i__ <= i__1; ++i__) {
	    if(fabs(d__[i__]) >= fabs(dl[i__])) {

/*              No row interchange required */

		if(d__[i__] != 0.) {
		    fact = dl[i__] / d__[i__];
		    d__[i__ + 1] -= fact * du[i__];
		    b_ref(i__ + 1, 1) = b_ref(i__ + 1, 1) - fact * b_ref(i__, 
			    1);
		} else {
		    *info = i__;
		    return 0;
		}
		dl[i__] = 0.;
	    } else {

/*              Interchange rows I and I+1 */

		fact = d__[i__] / dl[i__];
		d__[i__] = dl[i__];
		temp = d__[i__ + 1];
		d__[i__ + 1] = du[i__] - fact * temp;
		dl[i__] = du[i__ + 1];
		du[i__ + 1] = -fact * dl[i__];
		du[i__] = temp;
		temp = b_ref(i__, 1);
		b_ref(i__, 1) = b_ref(i__ + 1, 1);
		b_ref(i__ + 1, 1) = temp - fact * b_ref(i__ + 1, 1);
	    }
/* L10: */
	}
	if(*n > 1) {
	    i__ = *n - 1;
	    if(fabs(d__[i__]) >= fabs(dl[i__])) {
		if(d__[i__] != 0.) {
		    fact = dl[i__] / d__[i__];
		    d__[i__ + 1] -= fact * du[i__];
		    b_ref(i__ + 1, 1) = b_ref(i__ + 1, 1) - fact * b_ref(i__, 
			    1);
		} else {
		    *info = i__;
		    return 0;
		}
	    } else {
		fact = d__[i__] / dl[i__];
		d__[i__] = dl[i__];
		temp = d__[i__ + 1];
		d__[i__ + 1] = du[i__] - fact * temp;
		du[i__] = temp;
		temp = b_ref(i__, 1);
		b_ref(i__, 1) = b_ref(i__ + 1, 1);
		b_ref(i__ + 1, 1) = temp - fact * b_ref(i__ + 1, 1);
	    }
	}
	if(d__[*n] == 0.) {
	    *info = *n;
	    return 0;
	}
    } else {
	i__1 = *n - 2;
	for(i__ = 1; i__ <= i__1; ++i__) {
	    if(fabs(d__[i__]) >= fabs(dl[i__])) {

/*              No row interchange required */

		if(d__[i__] != 0.) {
		    fact = dl[i__] / d__[i__];
		    d__[i__ + 1] -= fact * du[i__];
		    i__2 = *nrhs;
		    for(j = 1; j <= i__2; ++j) {
			b_ref(i__ + 1, j) = b_ref(i__ + 1, j) - fact * b_ref(
				i__, j);
/* L20: */
		    }
		} else {
		    *info = i__;
		    return 0;
		}
		dl[i__] = 0.;
	    } else {

/*              Interchange rows I and I+1 */

		fact = d__[i__] / dl[i__];
		d__[i__] = dl[i__];
		temp = d__[i__ + 1];
		d__[i__ + 1] = du[i__] - fact * temp;
		dl[i__] = du[i__ + 1];
		du[i__ + 1] = -fact * dl[i__];
		du[i__] = temp;
		i__2 = *nrhs;
		for(j = 1; j <= i__2; ++j) {
		    temp = b_ref(i__, j);
		    b_ref(i__, j) = b_ref(i__ + 1, j);
		    b_ref(i__ + 1, j) = temp - fact * b_ref(i__ + 1, j);
/* L30: */
		}
	    }
/* L40: */
	}
	if(*n > 1) {
	    i__ = *n - 1;
	    if( fabs(d__[i__]) >= fabs(dl[i__])) 
		{
			if(d__[i__] != 0.) 
			{
				fact = dl[i__] / d__[i__];
				d__[i__ + 1] -= fact * du[i__];
				i__1 = *nrhs;
				for(j = 1; j <= i__1; ++j) {
				b_ref(i__ + 1, j) = b_ref(i__ + 1, j) - fact * b_ref(
					i__, j);
	/* L50: */
				}
			} 
			else 
			{
				*info = i__;
				return 0;
			}
	    } else {
		fact = d__[i__] / dl[i__];
		d__[i__] = dl[i__];
		temp = d__[i__ + 1];
		d__[i__ + 1] = du[i__] - fact * temp;
		du[i__] = temp;
		i__1 = *nrhs;
		for(j = 1; j <= i__1; ++j) {
		    temp = b_ref(i__, j);
		    b_ref(i__, j) = b_ref(i__ + 1, j);
		    b_ref(i__ + 1, j) = temp - fact * b_ref(i__ + 1, j);
/* L60: */
		}
	    }
	}
	if(d__[*n] == 0.) {
	    *info = *n;
	    return 0;
	}
    }

/*     Back solve with the matrix U from the factorization. */

    if(*nrhs <= 2) {
	j = 1;
L70:
	b_ref(*n, j) = b_ref(*n, j) / d__[*n];
	if(*n > 1) {
	    b_ref(*n - 1, j) = (b_ref(*n - 1, j) - du[*n - 1] * b_ref(*n, j)) 
		    / d__[*n - 1];
	}
	for(i__ = *n - 2; i__ >= 1; --i__) {
	    b_ref(i__, j) = (b_ref(i__, j) - du[i__] * b_ref(i__ + 1, j) - dl[
		    i__] * b_ref(i__ + 2, j)) / d__[i__];
/* L80: */
	}
	if(j < *nrhs) {
	    ++j;
	    goto L70;
	}
    } else {
	i__1 = *nrhs;
	for(j = 1; j <= i__1; ++j) {
	    b_ref(*n, j) = b_ref(*n, j) / d__[*n];
	    if(*n > 1) {
		b_ref(*n - 1, j) = (b_ref(*n - 1, j) - du[*n - 1] * b_ref(*n, 
			j)) / d__[*n - 1];
	    }
	    for(i__ = *n - 2; i__ >= 1; --i__) {
		b_ref(i__, j) = (b_ref(i__, j) - du[i__] * b_ref(i__ + 1, j) 
			- dl[i__] * b_ref(i__ + 2, j)) / d__[i__];
/* L90: */
	    }
/* L100: */
	}
    }

    return 0;

/*     End of DGTSV */

} /* dgtsv_ */
#endif
#undef b_ref
#endif
/*lint +e725 expected pos indentation */
/*lint +e801 goto is deprecated */
