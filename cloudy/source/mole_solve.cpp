/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CO_step fills in matrix for heavy elements molecular routines */
#include "cddefines.h"
#include "deuterium.h"
#include "ionbal.h"
#include "phycon.h"
#include "hmi.h"
#include "dynamics.h"
#include "conv.h"
#include "trace.h"
#include "grainvar.h"
#include "newton_step.h"
#include "h2.h"
#include "mole_priv.h"
#include "mole.h"
#include "dense.h"
/* Nick Abel between July and October of 2003 assisted Dr. Ferland in
 * improving the heavy element molecular network in Cloudy. Before
 * this routine would predict negative abundances if the fraction of
 * carbon in the form of molecules came close to 100%. A reorganizing
 * of the reaction network detected several bugs.  Treatment of
 * "coupled reactions", in which both densities in the reaction rate
 * were being predicted by Cloudy, were also added.  Due to these
 * improvements, Cloudy can now perform calculations where 100% of the
 * carbon is in the form of CO without predicting negative abundances
 *
 * Additional changes were made in November of 2003 so that our reaction 
 * network would include all reactions from the TH85 paper.  This involved 
 * adding silicon to the chemical network.  Also the reaction rates were
 * labeled to make identification with the reaction easier and the matrix 
 * elements of atomic C, O, and Si are now done in a loop, which makes 
 * the addition of future chemical species (like N or S) easy.
 * */
void check_co_ion_converge(void);

STATIC void funjac(GroupMap &MoleMap, const valarray<double> &b2vec,
						 double * const ervals, double * const amat, const bool lgJac, bool *lgConserve);

STATIC void mole_h_fixup(void);

STATIC void grouped_elems(const double bvec[], double mole_elems[]);

#define SMALLABUND 1e-24

double mole_solve()
{
	int nBad, nPrevBad, i;
	realnum
		eqerror, error;
	valarray<double> oldmols(mole_global.num_compacted),
		newmols(mole_global.num_compacted);
	GroupMap MoleMap( nuclide_list.size() );

	DEBUG_ENTRY( "mole_solve()" );

	if (hmi.H2_frac_abund_set>0.)
	{
		mole_h_fixup();
		fixit("Need to treat hmi.H2_frac_abund_set in revised network"); 
		fprintf(stderr,"Need to treat hmi.H2_frac_abund_set in revised network\n");
		fprintf(stderr,"%g\n",hmi.H2_frac_abund_set);
		exit(-1);
	}

	ASSERT(lgElemsConserved());

	/* will test against this error for quitting */
	const double BIGERROR = 1e-8;
	/* >>chng 04 jul 20, upper limit had been 6, why?  change to 20
	 * to get closer to soln */
	const int LIM_LOOP = 100;
	/* loop until run to limit, or no fix ups needed and error is small */
	/* will be used to keep track of neg solns */
	nBad = nPrevBad = 0;
	eqerror = -1.;

	valarray<double> escale(mole_global.num_compacted);
	double rlimit=-1., rmax=0.0;
		
	// Run enough times through to converge nonlinear system.
	for(i=0; i < LIM_LOOP;i++) 
	{
		if (trace.lgTrace)
		{
			fprintf(ioQQQ," Mole solve loop %d\n",i);
		}

		{
			/* option to print deduced abundances */
			/*@-redef@*/
			enum {DEBUG_LOC=false};
			/*@+redef@*/
			if( DEBUG_LOC && (nzone>140) )
			{
				fprintf(ioQQQ,"DEBUG hmole in \t%.2f\t%.5e\t%.5e",
					fnzone,
					phycon.te,
					dense.eden);
				for( long k=0; k<mole_global.num_calc; k++ )
					fprintf(ioQQQ,"\t%.2e", mole.species[k].den );
				fprintf(ioQQQ,"\n" );
			}
		}
		/* nBad returns number of times negative abundances were fixed
		 * for latest step, should be zero at return for valid soln */

		nPrevBad = nBad;

		MoleMap.setup(get_ptr(oldmols));

		long j;

		// Need to allow enough loops for j that rlimit increases until
		// no species get to -ve densities.  Often, though, this will
		// just be the first time through.
		static ConvergenceCounter cctr=conv.register_("MOLE_SOLVE");
		++cctr;
		for (j=0; j<30; j++) 
		{
			static ConvergenceCounter cctrs=conv.register_("MOLE_SOLVE_STEPS");
			++cctrs;
			bool lgOK = newton_step(MoleMap, oldmols, newmols, &eqerror, &error, mole_global.num_compacted, &rlimit, &rmax, escale, funjac); 
		
			/* check for negative populations */
			if (lgOK)
			{
				nBad = 0; 
				double oldMolsTmp = 0.;
				double newMolsTmp = 0.;
				long iworst = -1;
				for( long mol=0; mol < mole_global.num_compacted; mol++ )
				{
					if( newmols[mol] < 0. ) 
					{
						if( -newmols[mol]*oldMolsTmp >= oldmols[mol]*newMolsTmp )
						{
							
							oldMolsTmp = abs(oldmols[mol]);
							newMolsTmp = abs(newmols[mol]);
							iworst = mol;
						}
						++nBad;
						newmols[mol] = 0.;
					}
				}
				if ( 0 != nBad)
				{
					lgOK = false;
					if (0)
					{
						fprintf(ioQQQ,"-ve density in inner chemistry loop %ld, worst is %s rlimit %g\n",j,groupspecies[iworst]->label.c_str(),rlimit);
					}
				}
			}
			if ( lgOK )
			{
				//fprintf(ioQQQ,"OK at z %ld l %d j %ld t %g e %g\n",nzone,i,j,1./rlimit,eqerror);
				break;
			}

			ASSERT( rlimit < BIGFLOAT );
			ASSERT( rlimit > 0.0 );
			rlimit = 10.*rlimit;
		}
		
		//fprintf(stdout,"Mole zone %ld -- %7.2f error %15.8g eqerror %15.8g rlimit %15.8g nBad %d ninner %ld loop %d\n",nzone,fnzone,error,eqerror,rlimit,nBad,j+1,i);

		// If the accuracy of the solution obtained was strongly limited
		// by the pseudo-timestep limit, extend it
		if ( error < 0.01*eqerror )
			rlimit = 0.1*rlimit;
		
		MoleMap.updateMolecules( newmols );
			
		/* Stop early if good enough */
		if (eqerror < BIGERROR && nBad == 0 && nPrevBad == 0) 
			break;
	}

	//fprintf(stdout,"Mole: zn %ld -- %7.2f err %13.6g rlimit %13.6g nBad %d loop %d\n",nzone,fnzone,eqerror,rlimit,nBad,i);

	if( (i == LIM_LOOP && eqerror > BIGERROR)  || nBad != 0 )
	{
		conv.setConvIonizFail("Chem conv", eqerror, nBad);
	}
	
	// Effect_delta determines whether changes due to the molecular
	// network rescale values enough to invalidate other solutions
	realnum effect_delta = mole_return_cached_species(MoleMap);

	realnum delta_threshold = 0.2*conv.IonizErrorAllowed;
	if (effect_delta > delta_threshold)
	{
		conv.setConvIonizFail("chem eft chg", effect_delta, 0.0);
	}		

	if (trace.lgTrace)
	{
		fprintf(ioQQQ,"Mole zn %3ld -- %7.2f\n",nzone,fnzone); 
		fprintf(ioQQQ,"Internal error %15.8g nBad %4d loops %3d\n",
				  eqerror,nBad,i);
		fprintf(ioQQQ,"Scaling delta on ions %15.8g; threshold %15.8g\n", 
				  effect_delta, delta_threshold);
	}

	check_co_ion_converge();

	{
		/* option to print deduced abundances */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC /*&& (nzone>68)*/ )
		{
			fprintf(ioQQQ,"DEBUG mole out\t%i\t%.2f\t%.5e\t%.5e",
							i,
							fnzone,
							phycon.te,
							dense.eden);
			fprintf(ioQQQ,
							"\terror\t%.4e\tH0\t%.4e\tH+\t%.4e\tsink\t%.4e\t%.4e\tsour\t%.4e\t%.4e\tion\t%.4e\trec\t%.4e", 
							eqerror,
							dense.xIonDense[ipHYDROGEN][0],
							dense.xIonDense[ipHYDROGEN][1],
							mole.sink[ipHYDROGEN][0],
							mole.sink[ipHYDROGEN][1],
							mole.source[ipHYDROGEN][0] , 
							mole.source[ipHYDROGEN][1] ,
							ionbal.RateIonizTot(ipHYDROGEN,0),
							ionbal.RateRecomTot[ipHYDROGEN][0]);
			
			/*for( j=0; j<mole_global.num_calc; j++ )
				fprintf(ioQQQ,"\t%.4e", HMOLEC(j) );*/
			fprintf(ioQQQ,"\n" );
		}
	}
	
	return eqerror;

/*lint +e550 */
/*lint +e778 constant expression evaluates to 0 in operation '-' */
}


void check_co_ion_converge(void)
{
	DEBUG_ENTRY( "check_co_ion_converge()" );
	/* check whether ion and chem solvers agree yet */
	if( dense.lgElmtOn[ipCARBON] &&
		fabs(dense.xIonDense[ipCARBON][0]- findspecieslocal("C")->den)/SDIV(dense.gas_phase[ipCARBON]) >0.001 )
	{
		conv.setConvIonizFail("CO C0 con",
									 dense.xIonDense[ipCARBON][0],
									 findspecieslocal("C")->den);
	}
	else if( dense.lgElmtOn[ipCARBON] &&
		fabs(dense.xIonDense[ipCARBON][1]- findspecieslocal("C+")->den)/SDIV(dense.gas_phase[ipCARBON]) >0.001 )
	{
		conv.setConvIonizFail("CO C1 con",
									 dense.xIonDense[ipCARBON][1],
									 findspecieslocal("C+")->den);
	}
	else if( dense.lgElmtOn[ipOXYGEN] &&
		fabs(dense.xIonDense[ipOXYGEN][0]- findspecieslocal("O")->den)/SDIV(dense.gas_phase[ipOXYGEN]) >0.001 )
	{
		conv.setConvIonizFail(
			"CO O0 con",dense.xIonDense[ipOXYGEN][0],
			findspecieslocal("O")->den);
	}
	else if( dense.lgElmtOn[ipOXYGEN] &&
		fabs(dense.xIonDense[ipOXYGEN][1]- findspecieslocal("O+")->den)/SDIV(dense.gas_phase[ipOXYGEN]) >0.001 )
	{
		conv.setConvIonizFail(
			"CO O1 con", dense.xIonDense[ipOXYGEN][1],
			findspecieslocal("O+")->den);
	}
}

STATIC void mole_h_fixup(void)
{
	long int mol;
	
	/* there are two "no molecules" options, the no co, which turns off everything,
	 * and the no n2, which only turns off the h2.  in order to not kill the co
	 * part we still need to compute the hydrogen network here, and then set h2 to
	 * small values */
	/* >> chng 03 jan 15 rjrw -- suddenly switching off molecules confuses the solvers... */
	DEBUG_ENTRY( "mole_h_fixup()" );

	/* >>chng 02 jun 19, add option to force H2 abundance, for testing h2 molecules,
	 * hmi.H2_frac_abund_set is fraction in molecules that is set by set h2 fraction command */
	if( hmi.H2_frac_abund_set>0.)
	{
		for(mol=0;mol<mole_global.num_calc;mol++) 
		{
			mole.species[mol].den = 0.;
		}
		/* >>chng 03 jul 19, from 0 to SMALLFLOAT, to pass asserts in ConvBase,
		 * problem is that ion range has not been reset for hydrogen */
		dense.xIonDense[ipHYDROGEN][0] = dense.xIonDense[ipHYDROGEN][1] = 
			2.f*SMALLFLOAT*dense.gas_phase[ipHYDROGEN];
		/* put it all in the ground state */
		findspecieslocal("H2")->den = (realnum)(dense.gas_phase[ipHYDROGEN] * hmi.H2_frac_abund_set);
		findspecieslocal("H2*")->den = 0.;

		hmi.H2_total = findspecieslocal("H2")->den + findspecieslocal("H2*")->den;
		/* first guess at ortho and para densities */
		h2.ortho_density = 0.75*hmi.H2_total;
		h2.para_density = 0.25*hmi.H2_total;
		{
			hmi.H2_total_f = (realnum)hmi.H2_total;
			h2.ortho_density_f = (realnum)h2.ortho_density;
			h2.para_density_f = (realnum)h2.para_density;
		}

		hmi.hmihet = 0.;
		hmi.h2plus_exc_frac = 0.;
		hmi.h2plus_heatcoef = 0.;
		hmi.h2plus_heat = 0.;
		hmi.H2Opacity = 0.;
		hmi.hmicol = 0.;
		hmi.HeatH2Dish_TH85 = 0.;
		hmi.HeatH2Dexc_TH85 = 0.;
		hmi.deriv_HeatH2Dexc_TH85 = 0.;
		hmi.hmidep = 1.;

		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			gv.bin[nd].rate_h2_form_grains_used = 0.;
		}

		return;
	}
}

#define ABSLIM  1e-12
#define ERRLIM  1e-12
#	ifdef MAT
#		undef MAT
#	endif
#	define MAT(a,I_,J_)	((a)[(I_)*(mole_global.num_compacted)+(J_)])


STATIC void mole_eval_dynamic_balance(long int num_total, double *b, bool lgJac, multi_arr<double,2> &c);
enum {PRINTSOL = false};

STATIC void funjac(GroupMap &MoleMap, const valarray<double> &b2vec, double * const ervals,
						 double * const amat, const bool lgJac, bool *lgConserve)
{
	static multi_arr<double,2> c(mole_global.num_total,mole_global.num_total);
	bool printsol = PRINTSOL;
	valarray<double> b(mole_global.num_total);

	DEBUG_ENTRY( "funjac()" );

	MoleMap.updateMolecules( b2vec );		
		
	if( dynamics.doNonEquilibriumSolve( iteration ) )
	{
		ASSERT(dynamics.Rate > 0.0);
		*lgConserve = false; 
	}
	else
	{
		*lgConserve = true; 
	}		

	/* Generate chemical balance vector (mole.b[]) and Jacobian array 
		 (c[][], first iteration) from reaction list */
	
	mole_eval_dynamic_balance(mole_global.num_total,get_ptr(b),lgJac,c);	
	
	/*------------------------------------------------------------------ */
	if(printsol || (trace.lgTrace && trace.lgTraceMole ))
	{
		/* print the full matrix */
		fprintf( ioQQQ, "                ");
		for( long i=0; i < mole_global.num_calc; i++ )
		{
			fprintf( ioQQQ, "      %-4.4s", mole_global.list[i]->label.c_str() );
		}
		fprintf( ioQQQ, " \n" );
		
		fprintf(ioQQQ,"       MOLE old abundances\t%.2f",fnzone);
		for( long i=0; i<mole_global.num_calc; i++ )
			fprintf(ioQQQ,"\t%.2e", mole.species[i].den );
		fprintf(ioQQQ,"\n" );
		
		for( long i=0; i < mole_global.num_calc; i++ )
		{
			fprintf( ioQQQ, "       MOLE%2ld %-4.4s", i ,mole_global.list[i]->label.c_str());
			for( long j=0; j < mole_global.num_calc; j++ )
			{
				fprintf( ioQQQ, "%10.2e", c[j][i] );
			}
			fprintf( ioQQQ, "%10.2e", b[i] );
			fprintf( ioQQQ, "\n" );
		}
	}

	/* add positive ions and neutral atoms: ratios are set by ion_solver,
	 * we determine abundance of the group as a whole here */
	
	if (lgJac) {
		for(long i=0;i<mole_global.num_calc;i++) 
		{
			ASSERT(mole_global.list[i]->index == i);
		}
		for (unsigned long j=0; j<nuclide_list.size(); ++j )
		{
			vector<int> &jlist = nuclide_list[j]->ipMl;
			if (jlist[0] != -1)
			{
				for(long i=0;i<mole_global.num_calc;i++) 
				{
					c[jlist[0]][i] *= MoleMap.fion[j][0];
				}
				for (unsigned long ion=1;ion<jlist.size();ion++) 
				{
					double fion = MoleMap.fion[j][ion];
					if (jlist[ion] != -1 && fion != 0.0)
					{
						for(long i=0;i<mole_global.num_calc;i++) 
						{
							c[jlist[0]][i] += fion*c[jlist[ion]][i];
						}
					}
				}
			}
		}
		for (unsigned long j=0; j<nuclide_list.size(); ++j )
		{
			vector<int> &jlist = nuclide_list[j]->ipMl;
			if (jlist[0] != -1)
			{
				for(long i=0;i<mole_global.num_calc;i++) 
				{
					double sum = 0.0;

					for (unsigned long ion=0;ion<jlist.size();ion++)
					{
						if (jlist[ion] != -1)
						{
							sum += c[i][jlist[ion]];
						}
					}
					c[i][jlist[0]] = sum; 
				}					
			}
		}
	}
	
	for (unsigned long j=0; j<nuclide_list.size(); ++j )
	{
		//if (j == 8)
		//	fprintf(ioQQQ,"Nsum1 %s %ld %d %d\n",
		//			  nuclide_list[j]->label().c_str(), j, groupspecies[29]->index,nuclide_list[j]->ipMl[0]);
		vector<int> &jlist = nuclide_list[j]->ipMl;
		if (jlist[0] != -1)
		{
			double sum = 0.0;
			for (unsigned long ion=0;ion<jlist.size();ion++)
			{	
				if (jlist[ion] != -1)
				{
					//			if (j == 8)
					//	fprintf(ioQQQ,"Nsum %d %15.8g\n",nuclide_list[j]->ipMl[ion],
					//			  b[nuclide_list[j]->ipMl[ion]]);
					sum += b[jlist[ion]];
					b[jlist[ion]] = 0.0;
				}
			}
			b[jlist[0]] = sum;
		}
	}
	
	/* Species are now grouped -- only mole_global.num_compacted elements now active */

	if (*lgConserve)
	{
		// Replace atom rows with conservation constraints
		if (lgJac)
		{
			for ( unsigned long j = 0; j < nuclide_list.size(); ++j )
			{
				long ncons = nuclide_list[j]->ipMl[0];
				if (ncons != -1)
				{
					ASSERT( mole_global.list[ncons]->isMonatomic() );
					double scale=1.0/MoleMap.molElems[j]; //fabs(c[ncons][ncons]);// 
					for( long i=0;i<mole_global.num_compacted;i++)
					{
						if( groupspecies[i]->nNuclide.find(nuclide_list[j]) != groupspecies[i]->nNuclide.end() )
							c[groupspecies[i]->index][ncons] = groupspecies[i]->nNuclide[nuclide_list[j]]*scale;
						else
							c[groupspecies[i]->index][ncons] = 0.;
						
					}
				}
			}
		}

		valarray<double> molnow(nuclide_list.size());
		grouped_elems(get_ptr(b2vec), get_ptr(molnow));
		for ( unsigned long j = 0; j < nuclide_list.size(); ++j )
		{
			long ncons = nuclide_list[j]->ipMl[0];
			if (ncons != -1)
			{
				ASSERT( mole_global.list[ncons]->isMonatomic() );
				double scale = c[ncons][ncons];
				b[ncons] = (molnow[j]-MoleMap.molElems[j])*scale;
				if (false)
					if (b[ncons] != 0.0)
						fprintf(ioQQQ,"Cons %s err %g rel %g\n",
								  nuclide_list[j]->label().c_str(),b[ncons],
								  (molnow[j]-MoleMap.molElems[j])/SDIV(MoleMap.molElems[j]));
			}
		}
	}
	
	for( long i=0; i < mole_global.num_compacted; i++ )
	{
		ervals[i] = b[groupspecies[i]->index];
	}
	
	if (lgJac)
	{

		for( long i=0; i < mole_global.num_compacted; i++ )
		{
			for( long j=0; j < mole_global.num_compacted; j++ )
			{
				MAT(amat,i,j) = c[groupspecies[i]->index][groupspecies[j]->index];
			}			 
		}

		// Replace rows for species with no sources and sinks with
		// an identity
		for(long i=0;i<mole_global.num_compacted;i++)
		{
			double sum1 = 0.;
			for (long j=0;j<mole_global.num_compacted;j++)
			{
				sum1 += fabs(MAT(amat,j,i));
			}
			if (sum1 == 0.0)
			{
				ASSERT(ervals[i] == 0.0);
				MAT(amat,i,i) = 1.0;
				// fprintf(ioQQQ,"Fixing %ld %s\n",i,groupspecies[i]->label.c_str());
			}
		}
	}
}

STATIC void grouped_elems(const double bvec[], double mole_elems[])
{
	map<chem_nuclide*, long> nuclide_to_index;
	for (unsigned long j=0; j<nuclide_list.size(); ++j )
	{
		mole_elems[j] = 0.;
		nuclide_to_index[nuclide_list[j].get()] = j;
	}
	for( long i=0; i < mole_global.num_compacted; i++ )
	{
		if( groupspecies[i]->isIsotopicTotalSpecies() == false )
			continue;
		
		for (molecule::nNucsMap::const_iterator el = groupspecies[i]->nNuclide.begin(); 
			  el != groupspecies[i]->nNuclide.end(); ++el)
		{
			mole_elems[nuclide_to_index[el->first.get()]] += bvec[i]*el->second;
		}
	}
}

void GroupMap::setup(double *b0vec)
{
	valarray<double> calcv(mole_global.num_total);
	bool lgSet;

	for( long i=0;i<mole_global.num_total;i++)
	{
		calcv[i] = mole.species[i].den;
	}

	for (unsigned long j=0; j<nuclide_list.size(); ++j )
	{
		if (nuclide_list[j]->ipMl[0] != -1)
		{
			double sum = 0.;
			for (unsigned long ion=0; ion<nuclide_list[j]->ipMl.size(); ion++)
			{
				if (nuclide_list[j]->ipMl[ion] != -1)
					sum += calcv[nuclide_list[j]->ipMl[ion]];
			}
			if (sum > SMALLFLOAT)
			{
				double factor = 1./sum;
				for (unsigned long ion=0; ion<nuclide_list[j]->ipMl.size(); ion++)
				{
					if (nuclide_list[j]->ipMl[ion] != -1)
						fion[j][ion] = calcv[nuclide_list[j]->ipMl[ion]]*factor;
					else
						fion[j][ion] = 0.;
				}
			}
			else
			{
				lgSet = false;
				for (unsigned long ion=0; ion<nuclide_list[j]->ipMl.size(); ion++)
				{
					if (nuclide_list[j]->ipMl[ion] != -1 && !lgSet)
					{
						fion[j][ion] = 1.0;
						lgSet = true;
					}
					else
					{
						fion[j][ion] = 0.;
					}
				}
			}
			
			lgSet = false;
			for (unsigned long ion=0; ion<nuclide_list[j]->ipMl.size(); ion++)
			{
				if (nuclide_list[j]->ipMl[ion] != -1)
				{
					if (!lgSet) 
						calcv[nuclide_list[j]->ipMl[ion]] = sum;
					else
						calcv[nuclide_list[j]->ipMl[ion]] = 0.;
					lgSet = true;
				}
			}
		}
	}	

	for( long i=0; i < mole_global.num_compacted; i++ )
	{
		b0vec[i] = calcv[groupspecies[i]->index];
	}
	
	grouped_elems(get_ptr(b0vec),get_ptr(molElems));

	for (unsigned long i = 0; i<nuclide_list.size(); ++i)
	{
		double densAtom = 0.;
		// deuterium is special case
		if( nuclide_list[i]->el()->Z==1 && nuclide_list[i]->A==2 )
		{
			ASSERT( deut.lgElmtOn );
			densAtom = deut.gas_phase;
		}
		// skip other isotopes
		else if( !nuclide_list[i]->lgMeanAbundance() )
			continue;
		else
		{
			int nelem = nuclide_list[i]->el()->Z-1;
			densAtom = dense.gas_phase[nelem];
		}
		bool lgTest = 
			( densAtom < SMALLABUND && molElems[i] < SMALLABUND ) ||
			( fabs(molElems[i]- densAtom) <= conv.GasPhaseAbundErrorAllowed*densAtom );
		if( !(lgTest || conv.lgSearch) )
		{
			fprintf( ioQQQ, "PROBLEM molElems BAD  %li\t%s\t%.12e\t%.12e\t%.12e\n", 
						i, nuclide_list[i]->label().c_str(), nuclide_list[i]->frac, densAtom, molElems[i] );
			fprintf( ioQQQ, "molElems[i] %.12e densAtom %.12e diff %.12e error %.12e lgSearch2 %i\n",
					molElems[i] , densAtom,
					fabs(molElems[i]- densAtom) , conv.GasPhaseAbundErrorAllowed*densAtom , conv.lgSearch );
		}
		
		ASSERT( lgTest || conv.lgSearch );
		molElems[i] = densAtom;
	}

}

void GroupMap::updateMolecules(const valarray<double> & b2 )
{
	DEBUG_ENTRY( "updateMolecules()" );

	for (long mol=0;mol<mole_global.num_calc;mol++)
	{
		mole.species[mol].den = 0.;
	}
	for (long mol=0;mol<mole_global.num_compacted;mol++)
	{
		mole.species[ groupspecies[mol]->index ].den = b2[mol];	/* put derived abundances back into appropriate molecular species */
	}

	// calculate isotopologue densities
	for (long mol=0;mol<mole_global.num_calc;mol++)
	{
		if( mole_global.list[mol]->parentIndex >= 0 )
		{
			ASSERT( !mole_global.list[mol]->isIsotopicTotalSpecies() );
			long parentIndex = mole_global.list[mol]->parentIndex;
			mole.species[mol].den = mole.species[parentIndex].den;
			for( nNucs_i it = mole_global.list[mol]->nNuclide.begin(); it != mole_global.list[mol]->nNuclide.end(); ++it )
			{
				if( !it->first->lgMeanAbundance() )
				{
					mole.species[mol].den *= pow( it->first->frac, it->second );
				}
			} 
		}
	}
	
	// calculate densities for monatomic species that had been collapsed into a group
	for (unsigned long j=0; j<nuclide_list.size(); ++j )
	{
		if (nuclide_list[j]->ipMl[0] != -1)
		{
			double grouptot = mole.species[nuclide_list[j]->ipMl[0]].den;
			double sum = 0.0;
			for (unsigned long ion=0;ion<nuclide_list[j]->ipMl.size();ion++)
			{
				if (nuclide_list[j]->ipMl[ion] != -1)
				{
					mole.species[nuclide_list[j]->ipMl[ion]].den = grouptot * fion[j][ion];
					sum += mole.species[nuclide_list[j]->ipMl[ion]].den;
				}
			}
			ASSERT(fabs(sum-grouptot) <= 1e-10 * fabs(grouptot));
		}
	}
	
	mole.set_isotope_abundances();

	return;
}

STATIC void mole_eval_dynamic_balance(long int num_total, double *b, bool lgJac, multi_arr<double,2> &c)
{
	double source;
	
	DEBUG_ENTRY( "mole_eval_dynamic_balance()" );

	mole_eval_balance(num_total, b, lgJac, c);

	/* >>chng 06 mar 17, comment out test on old full depth - keep old solution if overrun scale */
	if( dynamics.doNonEquilibriumSolve( iteration ) )
	{
		/* Don't use conservation form in matrix solution -- dynamics rate terms make c[][] non-singular */		

		for( long i=0;i<mole_global.num_calc;i++)
		{
			if (lgJac)
				c[i][i] -= dynamics.Rate;

			if( !mole_global.list[i]->isIsotopicTotalSpecies() )
				continue;

			b[i] -= mole.species[i].den*dynamics.Rate;

			if (!mole_global.list[i]->isMonatomic() || mole_global.list[i]->charge < 0 || ! mole_global.list[i]->lgGas_Phase ||
				( mole_global.list[i]->isMonatomic() && !mole_global.list[i]->nNuclide.begin()->first->lgMeanAbundance() ) )
			{
				b[i] += dynamics.molecules[i];
			}
			else if (mole_global.list[i]->charge == 0)
			{
				ASSERT( mole_global.list[i]->isMonatomic() );
				ASSERT( (int)mole_global.list[i]->nNuclide.size() == 1 );
				const shared_ptr<chem_nuclide>& atom = mole_global.list[i]->nNuclide.begin()->first;
				long nelem = atom->el()->Z-1;
				if( nelem >= LIMELM )						
					continue;
				source = 0.0;
				for ( long ion=dense.IonLow[nelem]; ion<=dense.IonHigh[nelem]; ++ion )
				{
					source += dynamics.Source[nelem][ion] * atom->frac;
				}
				b[i] += source;
			}
		}
	}
}
