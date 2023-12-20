/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CO_step fills in matrix for heavy elements molecular routines */
#include "cddefines.h"
#include "dense.h"
#include "trace.h"
#include "atmdat.h"
#include "mole_priv.h"
#include "mole.h"
/* Nick Abel between July and October of 2003 assisted Dr. Ferland in improving the heavy element 
 * molecular network in Cloudy. Before this routine would predict negative abundances if 
 * the fraction of carbon in the form of molecules came close to 100%. A reorganizing of 
 * the reaction network detected several bugs.  Treatment of "coupled reactions",
 * in which both densities in the reaction rate were being predicted by Cloudy, were also 
 * added.  Due to these improvements, Cloudy can now perform calculations
 * where 100% of the carbon is in the form of CO without predicting negative abundances
 *
 * Additional changes were made in November of 2003 so that our reaction 
 * network would include all reactions from the TH85 paper.  This involved 
 * adding silicon to the chemical network.  Also the reaction rates were
 * labeled to make identification with the reaction easier and the matrix 
 * elements of atomic C, O, and Si are now done in a loop, which makes 
 * the addition of future chemical species (like N or S) easy.
 * */
/* Robin Williams in August 2006 onwards reorganized the coding to cut down repetitions.  
 * This isolated several further bugs, and allows a sigificant number of lines of
 * code to be eliminated.  The balance of S2/S2+ amd ClO/ClO+ seems highly sensitive
 * (with small log scale results varying significantly if the order of arithmetic
 * operations is changed) -- I suspect this may imply a bug somewhere.
 * */
/*lint -e778 constant expression evaluatess to 0 in operation '-' */
/*=================================================================*/

STATIC bool lgNucleiConserved(const multi_arr<double,2> &c);

void mole_eval_balance(long int num_total, double *b, bool lgJac, multi_arr<double,2> &c)
{
	long int i, j;
	mole_reaction *rate;
	double rate_tot, rate_deriv[MAXREACTANTS], rk;
	molecule *sp;

	DEBUG_ENTRY( "mole_eval_balance()" );
	/* zero out array used for formation rates */
	for( i=0; i < num_total; i++ )
	{
		b[i] = 0.;
	}
	if (lgJac)
		c.zero();

	if (trace.lgTrace && trace.lgTraceMole)
	{
		for(long mol=0;mol<mole_global.num_calc;mol++) 
		{
			fprintf(ioQQQ," TrChemSp %20.13g %s\n",mole.species[mol].den,
					  mole_global.list[mol]->label.c_str());
		}
		for(mole_reaction_i p=mole_priv::reactab.begin(); 
			 p != mole_priv::reactab.end(); ++p) 
		{
			rate = &(*p->second);
			rk = mole.reaction_rks[ rate->index ];
			fprintf(ioQQQ," TrChem %20.13g %s\n",rk,rate->label.c_str());
		}
	}


	for(mole_reaction_i p=mole_priv::reactab.begin(); 
			p != mole_priv::reactab.end(); ++p) 
	{
		rate = &(*p->second);
		rk = mole.reaction_rks[ rate->index ];

		rate_tot = rk;
		for(i=0;i<rate->nreactants;i++)
		{
			rate_tot *= mole.species[ rate->reactants[i]->index ].den;
		}		
		
		if (trace.lgTrace && trace.lgTraceMole)
		{
			fprintf(ioQQQ," TrChemTot %20.13g %20.13g %s\n",rate_tot,rk,rate->label.c_str());			
			for(i=0;i<rate->nreactants;i++)
			{
				fprintf(ioQQQ, " %20.13g", mole.species[ rate->reactants[i]->index ].den);
			}
			fprintf(ioQQQ,"\n");
		}
		for(i=0;i<rate->nreactants;i++)
		{	
			sp = rate->reactants[i];
			if (rate->rvector[i] == NULL)
			{
				b[sp->index] -= rate_tot;
			}
		}
		
		for(i=0;i<rate->nproducts;i++)
		{
			sp = rate->products[i];
			if (rate->pvector[i] == NULL)
			{
				b[sp->index] += rate_tot;
			}
		}
		
		if (lgJac)
		{
			for(i=0;i<rate->nreactants;i++)
			{
				rate_deriv[i] = rk;
				for(j=0;j<rate->nreactants;j++)
				{
					if(i!=j)
					{
						rate_deriv[i] *= mole.species[ rate->reactants[j]->index ].den;
					}
				}
			}
			for(j=0;j<rate->nreactants;j++)
			{
				sp = rate->reactants[j];
				const double rated = rate_deriv[j];
				for(i=0;i<rate->nreactants;i++)
				{
					if (rate->rvector[i] == NULL)
						c[sp->index][rate->reactants[i]->index] -= rated;
				}
				for(i=0;i<rate->nproducts;i++)
				{
					if (rate->pvector[i] == NULL)
						c[sp->index][rate->products[i]->index] += rated;
				}
			}
		}
	}

	if (lgJac)
	{
		ASSERT( lgNucleiConserved(c) );
		// Make sure there appears to be a call-point
		// even when NDEBUG is off
		if (0) { void(lgNucleiConserved(c)); };
	}

	//mole_dominant_rates(findspecies("H+"),ioQQQ);
	return;
}

void mole_eval_sources(long int num_total)
{
	long int i, j, ion, ion2;
	mole_reaction *rate;
	double rate_tot, rate_deriv[MAXREACTANTS], rk;
	molecule *sp;

	DEBUG_ENTRY( "mole_eval_sources()" );
	/* zero out array used for formation rates */
	for( i=0; i < num_total; i++ )
	{
		mole.species[i].src = mole.species[i].snk = 0.;
	}

	for( long nelem=0; nelem< LIMELM; ++nelem )
	{
		/* these have one more ion than above */
		for( ion=0; ion<nelem+2; ++ion )
		{
			/* zero out the transfer array */
			for( ion2=0; ion2<nelem+2; ++ion2 )
			{
				mole.xMoleChTrRate[nelem][ion][ion2] = 0.;
			}
		}
	}

	for(mole_reaction_i p=mole_priv::reactab.begin(); 
			p != mole_priv::reactab.end(); ++p) 
	{
		rate = &(*p->second);
		rk = mole.reaction_rks[ rate->index ];
		
		for(i=0;i<rate->nreactants;i++)
		{
			rate_deriv[i] = rk;
			for(j=0;j<rate->nreactants;j++)
			{
				if(i!=j)
				{
					rate_deriv[i] *= mole.species[ rate->reactants[j]->index ].den;
				}
			}
		}
		
		rate_tot = rate_deriv[0] * mole.species[ rate->reactants[0]->index ].den;
		
		for(i=0;i<rate->nreactants;i++)
		{	
			sp = rate->reactants[i];
			if (rate->rvector[i] == NULL)
			{
				mole.species[sp->index].snk += rate_deriv[i];
			}
			else
			{
				if ( atmdat.lgCTOn )
				{
					for( molecule::nNucsMap::iterator nuc_i = sp->nNuclide.begin(); 
						  nuc_i != sp->nNuclide.end(); ++nuc_i)
					{
						if (! nuc_i->first->lgHasLinkedIon() )
							continue;
						ASSERT(nuc_i->second != 0);
						if (rate->rvector[i]->charge != sp->charge)
						{
							long nelem = nuc_i->first->el()->Z-1;
							mole.xMoleChTrRate[nelem][sp->charge][rate->rvector[i]->charge] +=
								(realnum) rate_deriv[i];
							break;
						}
					}
				}
			}
		}
		
		for(i=0;i<rate->nproducts;i++)
		{
			sp = rate->products[i];
			if (rate->pvector[i] == NULL)
			{
				mole.species[sp->index].src += rate_tot;
			}
		}
	}

	for (ChemNuclideList::iterator atom = nuclide_list.begin();
		  atom != nuclide_list.end(); ++atom)
	{
		if (!(*atom)->lgHasLinkedIon())
			continue;
		const long int nelem=(*atom)->el()->Z-1;
		if( !dense.lgElmtOn[nelem] )
			continue;

		for (long int ion=0;ion<nelem+2;ion++) 
		{
			if ((*atom)->ipMl[ion] != -1)
			{
				mole.source[nelem][ion] = mole.species[(*atom)->ipMl[ion]].src;
				mole.sink[nelem][ion] = mole.species[(*atom)->ipMl[ion]].snk;
			}
			else
			{
				mole.source[nelem][ion] = 0.0;
				mole.sink[nelem][ion] = 0.0;
			}
		}
	}

	//mole_dominant_rates(findspecies("H+"),ioQQQ);
	return;
}

STATIC bool lgNucleiConserved(const multi_arr<double,2> &c)
{
	DEBUG_ENTRY ("lgNucleiConserved()" );
	bool checkAllOK = true;

	for (unsigned long j=0; j<nuclide_list.size(); ++j )
	{
		nuclide_list[j]->index = j;
	}

	size_t size1=0;
	for (long j=0;j<mole_global.num_calc;j++) 
	{
		size1 += mole_global.list[j]->nNuclide.size();
	}
	vector<long> natoms(size1), nNucs(size1), jend(mole_global.num_calc);
	for (long j=0,i=0;j<mole_global.num_calc;j++) 
	{
		for (molecule::nNucsMap::const_iterator el = mole_global.list[j]->nNuclide.begin(); 
			  el != mole_global.list[j]->nNuclide.end(); ++el)
		{
			natoms[i] = el->first->index;
			nNucs[i] = el->second;
			++i;
		}
		jend[j] = long(i);
	}

	vector<double> test(nuclide_list.size()),
		tot(nuclide_list.size());
	for (long i=0;i<mole_global.num_calc;i++) 
	{
		for( unsigned long natom=0; natom < nuclide_list.size(); ++natom)
		{
			test[natom] = tot[natom] = 0.0;
		}
		for (long j=0,jstart=0;j<mole_global.num_calc;j++) 
		{
			long lim = jend[j];
			if (c[i][j] != 0.0)
			{
				for (long pos=jstart; pos < lim; ++pos)
				{
					const long natom = natoms[pos];
					const int nNuclidej = nNucs[pos];
					const double term = c[i][j] * nNuclidej;
					test[natom] += term;
					tot[natom] += fabs(term);
				}
			}
			jstart = lim;
		}

		for( unsigned long natom=0; natom < nuclide_list.size(); ++natom)
		{
			/*
			 * The following is Robin's response to a request to adjust the limit below
			 * (to 2e-9, from 1e-9) in order for a failing sim to complete successfully
			 * (Cloudy bailed with a conservation error in the Li H chemical network).
			 * The initial e-mail was sent on June 1, 2015, under the subject 'molecular
			 * balance tolerance', and Robin's response was posted the following day.
			 * The thread may be found on the cloudy-dev Google group.
			 *
			 * "Where this happens, there is generally a very rapid equilibrium reaction
			 *  between two species, so where all the destruction rates are added on the
			 *  diagonal of the Jacobian, you lose significance in the value of the slower
			 *  secular rate.
			 *
			 * "When the elements are added to make the total consistency check, generally
			 *  in a different order, this suffers from the problem that
			 *  ((A+B)-A)-B != (A+B)-(A+B) != 0. in finite precision mathematics.
			 *
			 * "As the size of the differences in rates can be very large, I used to worry
			 *  that this test was going to gradually loosen into meaninglessness.  On a
			 *  previous occasion, a version of the test was disabled because it was more
			 *  convenient to do this than deal with the real errors it was demonstrating,
			 *  and it took quite a while to pull things back.  But this has been stable
			 *  for long enough now that I don't think that loosening by a factor of 2
			 *  would be a major problem.
			 *
			 * "The only numerically robust way around this would be to ensure you always
			 *  add the rates in order from the smallest to the largest amplitude -- but
			 *  I'm not sure how to do this without effectively pushing the test back to
			 *  being one on the incoming rates, and generally making it too complex to
			 *  be a robust test.  [Cf. the testing of the HST mirror: this was figured
			 *  to a test which was very accurate but complex -- and in the end broken;
			 *  a simpler, more robust test which would have demonstrated the problem was
			 *  (reportedly) not used because they couldn't see the point in running a less
			 *  accurate test.]"
			 *
			 */
			/*>>chng 16 dec 16, from 2e-9 to 3e-9.  igm_primal fails with this print:
			 * PROBLEM Network conservation error Li H 7.10467e-41 2.46016e-09 1.56344e-05 3.65606e-07
			 * this occurred after changing iso n-changing to P&R.  Network is very sensitive
			 * to small changes for this sim, probably due to the small Li abundance?
			 * This first fail was in default gcc compile on radegund / gcc 6.2.  Fault is
			 * currently on gcc 6.2, all double, with default arith OK, and on icc
			 */
			const bool checkOK = 
				( fabs(test[natom]) <= MAX2(3e-9*tot[natom], 1e10*DBL_MIN) );
			if ( UNLIKELY(!checkOK) )
			{
				chem_nuclide *atom = nuclide_list[natom].get();
				fprintf( ioQQQ, " PROBLEM Network conservation error %s %s %g %g %g %g\n",
						  atom->label().c_str(),
						  mole_global.list[i]->label.c_str(),
						  test[natom],
						  test[natom]/tot[natom],
						  mole.species[atom->ipMl[0]].den,
						  mole.species[atom->ipMl[1]].den);
				//fprintf(stdout,"Problem at %s\n",rate->label);
				checkAllOK = false;
			}
		}
	}
	return checkAllOK;
}

