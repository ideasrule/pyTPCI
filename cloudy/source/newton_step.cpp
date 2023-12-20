/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*mole_newton_step line-search improvement in chemical network */
/* >> chng 02 nov 7 rjrw, Mole Moreliano:
 *   changes to linearized iterative form */
/*lint -e778 const express evaluates to 0 */
/*lint -e725 expect positive indentation */
#include "cddefines.h"
#include "newton_step.h"
#include "conv.h"
#include "thirdparty.h"
#include "mole.h"
#include "mole_priv.h"
#include "trace.h"
#include "save.h"

#	ifdef MAT
#		undef MAT
#	endif
#	define MAT(a,I_,J_)	((a)[(I_)*(n)+(J_)])

STATIC void mole_system_error(long n, long merror, 
											 const valarray<double> &a, 
											 const valarray<double> &b);

enum {PRINTSOL = false};

/* mole_newton_step -- improve balance in chemical network along
 * descent direction, step limited to ensure improvement */
bool newton_step(GroupMap &MoleMap, const valarray<double> &b0vec, valarray<double> &b2vec, 
					  realnum *eqerror, realnum *error, const long n, double *rlimit, double *rmax, 
					  valarray<double> &escale,
					  void (*jacobn)(GroupMap &MoleMap, 
							const valarray<double> &b2vec,
							double * const ervals, double * const amat,
							const bool lgJac, bool *lgConserve))
{
	bool lgOK=true;
	long int i, 
		loop;

	double 
		f1,
		f2,
		error2,
		error1,
		error0,
		erroreq,
		grad,
		pred;

	valarray<double> 
		amat(n*n),
		b1vec(n),
		ervals(n),
		ervals1(n),
		ervals0(n);
	int32 merror;

	DEBUG_ENTRY( "newton_step()" );

	for( i=0; i < n; i++ )
	{
		b2vec[i] = b0vec[i];
		ervals0[i] = 0.0;
		ervals1[i] = 0.0;
	}

	pred = error0 = error2 = erroreq = grad = f1 = f2 = 0.;

	static ConvergenceCounter cctr=conv.register_("NEWTON");
	++cctr;
	const int LOOPMAX = 40;
	for (loop=0;loop<LOOPMAX;loop++) 
	{
		bool lgConserve;
		jacobn(MoleMap, b2vec,get_ptr(ervals),get_ptr(amat),loop==0,&lgConserve);

		if (loop == 0) 
		{
			for( i=0; i < n; i++ )
			{
				escale[i] = MAT(amat,i,i);
			}
			// First time through this case, set rlimit by guesswork based on
			// maintaining matrix condition
			*rmax = 0.0;
			for( i=0; i<n; ++i )
			{
				if (-escale[i]>*rmax)
				{
					*rmax = -escale[i];
				}
			}
			if (*rlimit < 0.0) 
			{	
				// large generally more stable, small more quickly convergent
				*rlimit = 1e-19 * (*rmax);
			}
			else if (*rlimit > *rmax) 
			{	
				// large generally more stable, small more quickly convergent
				*rlimit = *rmax;
			}
			for( i=0; i < n; i++ )
			{
				if (! lgConserve || ( (! groupspecies[i]->isMonatomic()) || groupspecies[i]->charge < 0 ) ) 
				{
					// Only apply *rlimit to rows which haven't been
					// overwritten by conservation constraints
					MAT(amat,i,i) -= *rlimit;
				}
			}
		}

		// External error limit based on difference to state relaxed to
		// equilibrium
		erroreq = 0.;
		for( i=0; i < n; i++ )
		{
			double etmp = ervals[i]/(SMALLABUND*(*rmax)+fabs(b0vec[i]*escale[i]));
			erroreq += etmp*etmp;
		}

		if (trace.lgTrace)
		{
			fprintf(ioQQQ,"Newton step loop %ld error %15.8g\n",loop,erroreq);
		}

		// Internal (step) limit based on finite-(pseudo-)timestep solution
		for (i=0; i < n; i++)
		{
			ervals1[i] = ervals[i];
			if (! lgConserve || ( (! groupspecies[i]->isMonatomic()) || groupspecies[i]->charge < 0 ) ) 
				ervals1[i] -= (*rlimit)*(b2vec[i]-b0vec[i]);
		}

		error1 = 0.;
		long maxi = -1;
		double emax0 = 0.0, emax1 = 0.0;
		for( i=0; i < n; i++ )
		{
			if (! lgConserve || ( (! groupspecies[i]->isMonatomic()) || groupspecies[i]->charge < 0 ) ) 
			{
				/* Scale the errors weighting to ensure trace species are accurate */
				double etmp = ervals1[i]/(SMALLABUND*(*rmax)+fabs(b0vec[i]*escale[i]));
				double etmp0 = ervals0[i]/(SMALLABUND*(*rmax)+fabs(b0vec[i]*escale[i]));
				etmp *= etmp;
				etmp0 *= etmp0;
				if (fabs(etmp-etmp0) > fabs(emax1-emax0) )
				{
					maxi = i;
					emax1 = etmp;
					emax0 = etmp0;
				}
				error1 += etmp;
			}
		}

		if (loop == 0) 
		{
			for( i=0; i < n; i++ )
			{
				ervals0[i] = ervals1[i];
				b1vec[i] = ervals1[i];
			}
			merror = solve_system(amat,b1vec,n,mole_system_error);
			
			if (merror != 0) {
			  *eqerror = *error = 1e30f;
			  lgOK = false;
			  return lgOK;
			}
			error0 = error1;
			grad = -2*error0; 
			f1 = 1.0;

		} else {
			//fprintf(ioQQQ,"Newt1 %ld stp %11.4g err %11.4g erreq %11.4g rlimit %11.4g grd %11.4g fdg %11.4g e0 %11.4g de %11.4g\n",
			//		  loop,f1,error1,erroreq,*rlimit,grad,(error1-error0)/f1,error0,error1-error0);
			if (error1 < (1-2e-4*f1)*error0 || error1 < 1e-20)
			{
				break;
			} 
			// Backtrack using quadratic or cubic fit, ref Dennis & Schnabel 1996
			if (loop == 1) 
			{
				f2 = f1;
				f1 *= -0.5*f1*grad/(error1-error0-f1*grad);
				pred = error0+0.5*grad*f1;
			}
			else
			{
				if (0)
				{
					fprintf(ioQQQ,"Newt %ld stp %11.4g err %11.4g erreq %11.4g rlimit %11.4g grd %11.4g fdg %11.4g e0 %11.4g de %11.4g pred %11.4g\n",
							  loop,f1,error1,erroreq,*rlimit,grad,(error1-error0)/f1,error0,error1-error0,pred			  
						);
					if (maxi != -1)
						fprintf(ioQQQ,"Maxi %ld %s emax from %11.4g of %11.4g -> %11.4g of  %11.4g, diff %11.4g\n",
								  maxi,groupspecies[maxi]->label.c_str(),emax0,error0, emax1,error1,emax1-emax0);
				}
				// NB we must be careful how a is calculated because terms can be nearly equal and/or
 				// vastly different from other terms.
				double a = (error1-error0)/(f1*f1) - (error2-error0)/(f2*f2);
				a += grad * (1./f2 - 1./f1);
				// similarly calculate b
				double b = -f2*(error1-error0)/(f1*f1) + f1*(error2-error0)/(f2*f2);
				b += grad * (f2/f1 - f1/f2);
				double ft = f1;
				//fprintf(ioQQQ,"CHK    %15.8g %15.8g %15.8g %15.8g\n",error0,grad,a,b);
				//fprintf(ioQQQ,"Pred 1 %15.8g %15.8g %15.8g\n",error1,f1,error0+ft*(grad+ft/(ft-f2)*(b+a*ft)));
				//fprintf(ioQQQ,"Pred 2 %15.8g %15.8g %15.8g\n",error2,f2,error0+f2*(grad+f2/(ft-f2)*(b+a*f2)));				
				//f1 = (-b+sqrt(b*b-3.*a*grad*(ft-f2)))/(3.*a);
				//fprintf(ioQQQ,"Pred x %g %g\n",f1,error0+f1*(grad+f1/(ft-f2)*(b+a*f1)));
				if ( a != 0.0 )
				{
					f1 = 1.-3.*(a/b)*(grad/b)*(ft-f2);
					if (f1 > 0.) 
						f1 = b/(3.*a)*(sqrt(f1)-1.);
					else
						f1 = -b/(3.*a);
				}
				else
				{
					f1 = -grad/(2.*b);
				}
				//fprintf(ioQQQ,"Pred n %g %g\n",f1,error0+f1*(grad+f1/(ft-f2)*(b+a*f1)));				
				pred = error0+f1*(grad+f1/(ft-f2)*(b+a*f1));
				f2 = ft;
			}
			error2 = error1;
			if (f1 > 0.5*f2 || f1 < 0.)
				f1 = 0.5*f2;
			else if (f1 < 0.03*f2)
				f1 = 0.03*f2;			
		}
		
		/*
		 * b1vec[] is descent direction
		 * new densities in b2vec[]
		 */
		
		// Count number of attempts required to find new position
		static ConvergenceCounter cctrl=conv.register_("NEWTON_LOOP");
		++cctrl;
		if (f1 > 1e-6)
		{
			for( i=0; i < n; i++ )
			{
				b2vec[i] = b0vec[i]-b1vec[i]*f1;				
			}
		}
		else
		{
			// Try again with larger rlimit
			lgOK = false;
			for( i=0; i < n; i++ )
			{
				b2vec[i] = b0vec[i];
			}
			break;
		}
	}
	if (0 && LOOPMAX == loop)
	{
		double rvmax = 0., rval;
		int imax=0;
		for( i=0; i < n; ++i )
		{
			if (b0vec[i] != 0.)
				rval = b1vec[i]/b0vec[i];
			else
				rval = 0.;
			fprintf(ioQQQ,"%7s %11.4g %11.4g\n",
							groupspecies[i]->label.c_str(),rval,b0vec[i]);
			if (fabs(rval) > rvmax) 
			{
				rvmax = fabs(rval);
				imax = i;
			}
		}
		fprintf(ioQQQ,"Biggest is %s\n",groupspecies[imax]->label.c_str());
		if (0)
		{ // Verify Jacobian
			long j;
			for( j=0; j < n; j++ )
			{
				b1vec[j] = b0vec[j];
			}
			bool lgConserve;
			jacobn(MoleMap, b1vec,get_ptr(ervals1),get_ptr(amat),false,&lgConserve);
			for( i=0; i < n; i++ )
			{
				for( j=0; j < n; j++ )
				{
					b1vec[j] = b0vec[j];
				}
				double db = 1e-3*fabs(b0vec[i])+1e-9;
				b1vec[i] += db;
				db = b1vec[i]-b0vec[i];
				jacobn(MoleMap, b1vec,get_ptr(escale),get_ptr(amat),false,&lgConserve);
				for( j=0; j < n; j++ )
				{
					double e1 = MAT(amat,i,j);
					double e2 = (escale[j]-ervals1[j])/db;
					if (fabs(e1-e2) > 0.01*fabs(e1+e2))
						fprintf(ioQQQ,"%7s %7s %11.4g %11.4g %11.4g\n",
										groupspecies[i]->label.c_str(),groupspecies[j]->label.c_str(),e1,e2,
										ervals1[j]/db);
				}
			}
		}
		exit(-1);
	}

	*error = (realnum) MIN2(error1,1e30);
	*eqerror = (realnum) MIN2(erroreq,1e30);

	return lgOK;
}

STATIC void mole_system_error(long n, long merror, 
											 const valarray<double> &a, const valarray<double> &b)
{
	fprintf( ioQQQ, " CO_solve getrf_wrapper error %ld",(long int) merror );
	if( merror > 0 && merror <= n )
	{
		fprintf( ioQQQ," - problem with species %s\n\n", groupspecies[merror-1]->label.c_str() );
		fprintf( ioQQQ,"index \t Row A(i,%li)\t Col A(%li,j) \t B \t Species\n", merror, merror );
		for( long index=1; index<=n; index++ )
		{
			fprintf( ioQQQ,"%li\t%+.4e\t%+.4e\t%+.4e\t%s\n", index,
					a[(merror-1)*n + index - 1],
					a[(index -1)*n + merror- 1],
					b[index-1],
					groupspecies[index-1]->label.c_str() );
		}

		mole_print_species_reactions( groupspecies[merror-1] );
	}

	fprintf( ioQQQ,"\n" );
}

int32 solve_system(const valarray<double> &a, valarray<double> &b, 
									 long int n, error_print_t error_print)
{
	int32 merror;
	valarray<int32> ipiv(n);
	
	const int nrefine=3;
	int i, j, k;
	valarray<double> lufac(n*n),oldb(n),err(n);

	ASSERT(a.size() == size_t(n*n));
	ASSERT(b.size() == size_t(n));

	DEBUG_ENTRY( "solve_system()" );

	lufac = a;

	if (nrefine>0) 
	{
		oldb = b;
	}

	merror = 0;
	getrf_wrapper(n,n,get_ptr(lufac),n,get_ptr(ipiv),&merror);
	if ( merror != 0 )
	{
		if (error_print != NULL)
			error_print(n,merror,a,b);
		else
			fprintf(ioQQQ," PROBLEM singular matrix in solve_system\n");
		return merror; //cdEXIT(EXIT_FAILURE);
	}

	getrs_wrapper('N',n,1,get_ptr(lufac),n,get_ptr(ipiv),get_ptr(b),n,&merror);

	if ( merror != 0 )
	{
		fprintf( ioQQQ, " solve_system: dgetrs finds singular or ill-conditioned matrix\n" );
		return merror; //cdEXIT(EXIT_FAILURE);
	}

	for (k=0;k<nrefine;k++) 
	{
		for (i=0;i<n;i++)
		{
			err[i] = oldb[i];
		}
		for (j=0;j<n;j++) 
		{
			for (i=0;i<n;i++)
			{
				err[i] -= a[i+j*n]*b[j];
			}
		}
		getrs_wrapper('N',n,1,get_ptr(lufac),n,get_ptr(ipiv),get_ptr(err),n,&merror);
		if (0)
		{
			// Quick-and-dirty condition number estimate
			// see Golub & Van Loan, 3rd Edn, p128
			double maxb=0., maxe=0.;
			for (i=0;i<n;i++)
			{
				if (fabs(b[i])>maxb)
					maxb = fabs(b[i]);
				if (fabs(err[i])>maxe)
					maxe = fabs(err[i]);
			}
			if (k == 0)
			{
				fprintf(ioQQQ,"Max error %g, Condition %g\n",maxe,1e16*maxe/maxb);
			}
		}
		for (i=0;i<n;i++)
		{
			b[i] += err[i];
		}
	}

	return merror;
}
