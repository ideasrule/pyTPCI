/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef THIRDPARTY_QUADPACK_H_
#define THIRDPARTY_QUADPACK_H_

class E_fp
{
public:
	virtual sys_float operator()(sys_float) const = 0;
	E_fp() {}
	E_fp(const E_fp&) = default;
	E_fp& operator= (const E_fp&) = default;
protected:
	~E_fp() {}
};

class E_fp_fp : public E_fp
{
	sys_float (*m_f)(sys_float);
public:
	E_fp_fp(sys_float (*f)(sys_float)) : m_f(f) {}
	sys_float operator()(sys_float x) const
	{
		return m_f(x);
	}
};

class D_fp
{
public:
	virtual double operator()(double) const = 0;
	D_fp() {}
	D_fp(const D_fp&) = default;
	D_fp& operator= (const D_fp&) = default;
protected:
	~D_fp() {}
};

class D_fp_fp : public D_fp
{
	double (*m_f)(const double);
public:
	D_fp_fp(double (*f)(double)) : m_f(f) {}
	double operator()(double x) const
	{
		return m_f(x);
	}
};

typedef double (*D_fp1)(
	const double *, const double *, const double *, const double *, 
	const double *, const long *);
typedef sys_float (*E_fp1)(
	const sys_float *, const sys_float *, const sys_float *, const sys_float *,
	const sys_float *, const long *);

void dqage_(const D_fp& f, const double *a, const double *b, const double *epsabs, 
				const double *epsrel, const long *key, const long *limit, 
				double *result, double *abserr, long *neval, long *ier, double *
				alist__, double *blist, double *rlist, double *elist, 
				long *iord, long *last);
void dqag_(const D_fp& f, const double *a, const double *b, 
			  const double *epsabs, const double *epsrel, const long *key, 
			  double *result, 
			  double *abserr, long *neval, long *ier, long *limit, 
			  const long *lenw, long *last, long *iwork, double *work);
void dqagie_(const D_fp& f, const double *bound, const long *inf, 
				 const double *epsabs, const double *epsrel, const long *limit, 
				 double *result, double *abserr, long *neval, long *ier, double *
				 alist__, double *blist, double *rlist, double *elist, 
				 long *iord, long *last);
void dqagi_(const D_fp& f, const double *bound, const long *inf, 
				const double *epsabs, const double *epsrel, double *result, 
				double *abserr, long *neval, long *ier, long *limit, 
				const long *lenw, long *last, long *iwork, double *work);
void dqagpe_(const D_fp& f, const double *a, const double *b, 
				 const long *npts2, const double *polongs, const double *epsabs, 
				 const double *epsrel, 
				 const long *limit, double *result, double *abserr, long *
				 neval, long *ier, double *alist__, double *blist, 
				 double *rlist, double *elist, double *pts, long *iord, 
				 long *level, long *ndin, long *last);
void dqagp_(const D_fp& f, const double *a, const double *b, const long *npts2, 
				const double *polongs, const double *epsabs, const double *epsrel, 
				double *result, double *abserr, long *neval, long *ier, 
				const long *leniw, const long *lenw, long *last, long *iwork, 
				double *work);
void dqagse_(const D_fp& f, const double *a, const double *b, const double *epsabs, 
				 const double *epsrel, const long *limit, double *result, 
				 double *abserr, long *neval, long *ier, double *alist__,
				 double *blist, double *rlist, double *elist, long *
				 iord, long *last);
void dqags_(const D_fp& f, const double *a, const double *b, const double *epsabs, 
				const double *epsrel, double *result, double *abserr, 
				long *neval, long *ier, const long *limit, const long *lenw, long *last, 
				long *iwork, double *work);
void dqawce_(const D_fp& f, const double *a, const double *b, const double *c__, 
				 const double *epsabs, const double *epsrel, const long *limit, 
				 double *result, double *abserr, long *neval, long *ier, 
				 double *alist__, double *blist, double *rlist, double 
				 *elist, long *iord, long *last);
void dqawc_(const D_fp& f, const double *a, const double *b, const double *c__, 
				const double *epsabs, const double *epsrel, double *result, 
				double *abserr, long *neval, long *ier, long *limit, 
				const long *lenw, long *last, long *iwork, double *work);
void dqawfe_(const D_fp& f, const double *a, const double *omega, 
				 const long *integr, const double *epsabs, const long *limlst, 
				 const long *limit, const long *maxp1, 
				 double *result, double *abserr, long *
				 neval, long *ier, double *rslst, double *erlst, long *
				 ierlst, long *lst, double *alist__, double *blist, 
				 double *rlist, double *elist, long *iord, long *nnlog, 
				 double *chebmo);
void dqawf_(const D_fp& f, const double *a, const double *omega, const long *integr, 
				const double *epsabs, double *result, double *abserr, 
				long *neval, long *ier, long *limlst, long *lst, 
				const long *leniw, const long *maxp1, const long *lenw, 
				long *iwork, double *work);
void dqawoe_(const D_fp& f, const double *a, const double *b, const double *omega, 
				 const long *integr, const double *epsabs, const double *epsrel, 
				 const long *limit, const long *icall, const long *maxp1, 
				 double *result, 
				 double *abserr, long *neval, long *ier, long *last, 
				 double *alist__, double *blist, double *rlist, double 
				 *elist, long *iord, long *nnlog, long *momcom, double *
				 chebmo);
void dqawo_(const D_fp& f, const double *a, const double *b, const double *omega, 
				const long *integr, const double *epsabs, const double *epsrel, 
				double *result, double *abserr, long *neval, long *ier, 
				const long *leniw, long *maxp1, const long *lenw, long *last, long 
				*iwork, double *work);
void dqawse_(const D_fp& f, const double *a, const double *b, const double *alfa, 
				 const double *beta, const long *integr, const double *epsabs, 
				 const double *epsrel, const long *limit, double *result, 
				 double *abserr, long *neval, long *ier, double *alist__, 
				 double *blist, double *rlist, double *elist, long *iord, 
				 long *last);
void dqaws_(const D_fp& f, const double *a, const double *b, const double *alfa, 
				const double *beta, const long *integr, const double *epsabs, 
				const double *epsrel, double *result, double *abserr, 
				long *neval, long *ier, long *limit, const long *lenw, long *last, 
				long *iwork, double *work);
void dqc25c_(const D_fp& f, const double *a, const double *b, const double *c__, 
				 double *result, double *abserr, long *krul, long *neval);
void dqc25f_(const D_fp& f, const double *a, const double *b, const double *omega, 
				 const long *integr, const long *nrmom, const long *maxp1, 
				 const long *ksave, double *result, double *abserr, long *neval, 
				 double *resabs, double *resasc, long *momcom, double *
				 chebmo);
void dqc25s_(const D_fp& f, const double *a, const double *b, const double *bl, 
				 const double *br, const double *alfa, const double *beta, 
				 const double *ri, const double *rj, const double *rg, 
				 const double *rh, 
				 double *result, double *abserr, double *resasc, 
				 const long *integr, long *nev);
void dqcheb_(const double *x, double *fval, double *cheb12, double *cheb24);
void dqelg_(long *n, double *epstab, double *result, 
				double *abserr, double *res3la, long *nres);
void dqk15_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc);
void dqk15i_(const D_fp& f, const double *boun, const long *inf, 
				 const double *a, const double *b, double *result, double *abserr, 
				 double *resabs, double *resasc);
void dqk15w_(const D_fp& f, D_fp1 w, const double *p1, const double *p2, 
				 const double *p3, const double *p4, const long *kp, 
				 const double *a, 
				 const double *b, double *result, double *abserr, double *
				 resabs, double *resasc);
void dqk21_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc);
void dqk31_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc);
void dqk41_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc);
void dqk51_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc);
void dqk61_(const D_fp& f, const double *a, const double *b, double *
				result, double *abserr, double *resabs, double *resasc);
void dqmomo_(const double *alfa, const double *beta, double *
				 ri, double *rj, double *rg, double *rh, const long *integr);
void dqng_(const D_fp& f, const double *a, const double *b, const double *epsabs, 
			  const double *epsrel, double *result, double *abserr, 
			  long *neval, long *ier);
void dqpsrt_(const long *limit, long *last, long *maxerr, 
				 double *ermax, double *elist, long *iord, long *nrmax);
void qage_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  const sys_float *epsabs, const sys_float *epsrel, 
			  const long *key, const long *limit, 
			  sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, sys_float *alist__, sys_float *blist, sys_float *rlist,
			  sys_float *elist, long *iord, long *last);
void qag_(const E_fp& f, const sys_float *a, const sys_float *b, 
			 const sys_float *epsabs, const sys_float *epsrel, const long *key, 
			 sys_float *result, sys_float *abserr, long *neval, 
			 long *ier, long *limit, const long *lenw, long *last, long *
			 iwork, sys_float *work);
void qagie_(const E_fp& f, const sys_float *bound, const long *inf, 
				const sys_float *epsabs, const sys_float *epsrel, 
				const long *limit, sys_float *result, sys_float *abserr, 
				long *neval, long *ier, sys_float *alist__, sys_float *blist, 
				sys_float *rlist, sys_float *elist, long *iord, long *last);
void qagi_(const E_fp& f, const sys_float *bound, const long *inf, 
			  const sys_float *epsabs, const sys_float *epsrel, 
			  sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, const long *limit, const long *lenw, 
			  long *last, long *iwork, sys_float *work);
void qagpe_(const E_fp& f, const sys_float *a, const sys_float *b, const long *npts2, 
				const sys_float *polongs, const sys_float *epsabs, 
				const sys_float *epsrel, const long *limit, sys_float *result, 
				sys_float *abserr, long *neval, long *ier, sys_float *alist__, 
				sys_float *blist, sys_float *rlist, sys_float *elist, 
				sys_float *pts, long *iord, long *level, long *ndin, long *last);
void qagp_(const E_fp& f, const sys_float *a, const sys_float *b, const long *npts2,
			  const sys_float *polongs, const sys_float *epsabs, 
			  const sys_float *epsrel, sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, const long *leniw, const long *lenw, 
			  long *last, long *iwork, sys_float *work);
void qagse_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *epsabs, const sys_float *epsrel, 
				const long *limit, sys_float *result, sys_float *abserr, 
				long *neval, 
				long *ier, sys_float *alist__, sys_float *blist, sys_float *rlist, sys_float *elist, 
				long *iord, long *last);
void qags_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  const sys_float *epsabs, const sys_float *epsrel, 
			  sys_float *result, sys_float *abserr, long *neval, long *ier, 
			  const long *limit, const long *lenw, long *last, long *iwork,
			  sys_float *work);
void qawce_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *c__, const sys_float *epsabs,
				const sys_float *epsrel, const long *limit, 
				sys_float *result, sys_float *abserr, long *
				neval, long *ier, sys_float *alist__, sys_float *blist, 
				sys_float *rlist, sys_float *elist, long *iord, long *last);
void qawc_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  const sys_float *c__, const sys_float *epsabs, 
			  const sys_float *epsrel, sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, long *limit, const long *lenw, 
			  long *last, long *iwork, 
			  sys_float *work);
void qawfe_(const E_fp& f, const sys_float *a, const sys_float *omega, 
				const long *integr, const sys_float *epsabs, const long *limlst, 
				const long *limit, const long *maxp1, 
				sys_float *result, sys_float *abserr, long *neval, 
				long *ier, sys_float *rslst, sys_float *erlst, 
				long *ierlst, long *lst, sys_float *alist__, sys_float *blist, 
				sys_float *rlist, sys_float *elist, long *iord, 
				long *nnlog, sys_float *chebmo);
void qawf_(const E_fp& f, const sys_float *a, const sys_float *omega, 
			  const long *integr, 
			  const sys_float *epsabs, sys_float *result, 
			  sys_float *abserr, long *neval, long *ier, 
			  const long *limlst, long *lst, const long *leniw, const long *maxp1, 
			  const long *lenw, long *iwork, sys_float *work);
void qawoe_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *omega, const long *integr, 
				const sys_float *epsabs, const sys_float *epsrel, 
				const long *limit, const long *icall, const long *maxp1, 
				sys_float *result, sys_float *abserr, long *neval, long *
				ier, long *last, sys_float *alist__, sys_float *blist, sys_float *rlist, sys_float *
				elist, long *iord, long *nnlog, long *momcom, sys_float *chebmo);
void qawo_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  const sys_float *omega, const long *integr, 
			  const sys_float *epsabs, const sys_float *epsrel, 
			  sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, const long *leniw, const long *maxp1, 
			  const long *lenw, long *last, long *iwork, sys_float *work);
void qawse_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *alfa, const sys_float *beta, 
				const long *integr, const sys_float *epsabs, 
				const sys_float *epsrel, const long *limit, 
				sys_float *result, sys_float *abserr, long *neval, 
				long *ier, sys_float *alist__, 
				sys_float *blist, sys_float *rlist, sys_float *elist, long *iord, long *last);
void qaws_(const E_fp& f, const sys_float *a, const sys_float *b, 
			  const sys_float *alfa, const sys_float *beta, 
			  const long *integr, const sys_float *epsabs, 
			  const sys_float *epsrel, sys_float *result, sys_float *abserr, 
			  long *neval, long *ier, const long *limit, const long *lenw, 
			  long *last, long *iwork, sys_float *work);
void qc25c_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *c__, sys_float *result,
				sys_float *abserr, long *krul, long *neval);
void qc25f_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *omega, const long *integr, const long *nrmom, 
				const long *maxp1, const long *ksave, sys_float *result, 
				sys_float *abserr, long *neval, sys_float *resabs, sys_float *resasc, long *
				momcom, sys_float *chebmo);
void qc25s_(const E_fp& f, const sys_float *a, const sys_float *b, 
				const sys_float *bl, const sys_float *br, 
				const sys_float *alfa, const sys_float *beta, 
				sys_float *ri, sys_float *rj, sys_float *rg, sys_float *rh, 
				sys_float *result, sys_float *abserr, sys_float *resasc, 
				const long *integr, long *nev);
void qcheb_(sys_float *x, sys_float *fval, sys_float *cheb12, sys_float *cheb24);
void qelg_(long *n, sys_float *epstab, sys_float *result, sys_float *
			  abserr, sys_float *res3la, long *nres);
void qk15_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result, 
			  sys_float *abserr, sys_float *resabs, sys_float *resasc);
void qk15i_(const E_fp& f, const sys_float *boun, const long *inf, const sys_float *a,
				const sys_float *b, sys_float *result, sys_float *abserr, 
				sys_float *resabs, sys_float *resasc);
void qk15w_(const E_fp& f, E_fp1 w, const sys_float *p1, const sys_float *p2, 
				const sys_float *p3, const sys_float *p4, const long *kp, 
				const sys_float *a, const sys_float *b, sys_float *result, 
				sys_float *abserr, sys_float *resabs, sys_float *resasc);
void qk21_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result, 
			  sys_float *abserr, sys_float *resabs, sys_float *resasc);
void qk31_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result, 
			  sys_float *abserr, sys_float *resabs, sys_float *resasc);
void qk41_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result, 
			  sys_float *abserr, sys_float *resabs, sys_float *resasc);
void qk51_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result, 
			  sys_float *abserr, sys_float *resabs, sys_float *resasc);
void qk61_(const E_fp& f, const sys_float *a, const sys_float *b, sys_float *result, 
			  sys_float *abserr, sys_float *resabs, sys_float *resasc);
void qmomo_(const sys_float *alfa, const sys_float *beta, sys_float *ri, 
				sys_float *rj, sys_float *rg, sys_float *rh, const long *integr);
void qng_(const E_fp& f, sys_float *a, sys_float *b, sys_float *epsabs, sys_float *
			 epsrel, sys_float *result, sys_float *abserr, long *neval, long *ier);
void qpsrt_(const long *limit, long *last, long *maxerr, 
				sys_float *ermax, sys_float *elist, long *iord, long *nrmax);
double dqwgtc_(const double *x, const double *c__, const double *p2, 
					const double *p3, const double *p4, const long *kp);
double dqwgtf_(const double *x, const double *omega, const double *p2, 
					const double *p3, const double *p4, const long *integr);
double dqwgts_(const double *x, const double *a, const double *b, 
					const double *alfa, const double *beta, const long *integr);
sys_float qwgtc_(const sys_float *x, const sys_float *c__, const sys_float *p2,
					  const sys_float *p3, const sys_float *p4, const long *kp);
sys_float qwgtf_(const sys_float *x, const sys_float *omega, 
					  const sys_float *p2, const sys_float *p3, const sys_float *p4,
					  const long *integr);
sys_float qwgts_(const sys_float *x, const sys_float *a, const sys_float *b, 
					  const sys_float *alfa, const sys_float *beta, 
					  const long *integr);

#endif /* THIRDPARTY_QUADPACK_H_ */
