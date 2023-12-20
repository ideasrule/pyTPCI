/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolDima compute cooling due to level 2 lines */
/*ColStrGBar generate g-bar collision strengths for level 2 line2 */
#include "cddefines.h"
#include "taulines.h"
#include "dense.h"
#include "rt.h"
#include "doppvel.h"
#include "phycon.h"
#include "mewecoef.h"
#include "atoms.h"
#include "atmdat.h"
#include "thermal.h"
#include "cooling.h"

/*ColStrGBar generate g-bar collision strengths for level 2 line2 */
STATIC double ColStrGBar(const TransitionProxy::iterator& t , realnum cs1 );

void CoolDima(void)
{
	long int i, 
	  ion,
	  nelem;
	double cs;

	DEBUG_ENTRY( "CoolDima()" );

	/* no level2 command sets nWindLine to -1 */
	if( nWindLine<0 )
		return;

	static vector<realnum> DopplerWidth(LIMELM);
	for(nelem = ipHYDROGEN; nelem < LIMELM; ++nelem)
		DopplerWidth[nelem] = GetDopplerWidth(dense.AtomicWeight[nelem]);

	for( i=0; i < nWindLine; i++ )
	{
		ion = (*TauLine2[i].Hi()).IonStg();
		nelem = (*TauLine2[i].Hi()).nelem();

		if( (dense.lgIonChiantiOn[nelem-1][ion-1] && !atmdat.lgChiantiLvl2Hybrid) ||
				(dense.lgIonStoutOn[nelem-1][ion-1] && !atmdat.lgStoutLvl2Hybrid) )
		{
			/* If a species uses Chianti or Stout and hybrid is off, skip the level 2 lines */
			continue;
		}
		// iso sequence ions are done elsewhere
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO && 
			// dense.maxWN is positive if hybrid chianti is turned on and this element is included
			// in CloudyChianti.ini - zero otherwise
			 TauLine2[i].EnergyWN() > dense.maxWN[nelem-1][ion-1])
		{
			/* only evaluate cs if positive abundance */
			if( dense.xIonDense[nelem-1][ion-1] > 0. )
			{
				/* now generate the collision strength */
				cs = ColStrGBar(TauLine2.begin()+i , cs1_flag_lev2[i] );
			}
			else
			{
				cs = 1.;
			}
			/* now put the cs into the line array */
			PutCS(cs,TauLine2[i] );
			RT_line_one_escape( TauLine2[i], true,0.f, DopplerWidth[(*TauLine2[i].Hi()).nelem()-1] );
			atom_level2(TauLine2[i] );
			thermal.dima += TauLine2[i].Coll().cool();
		}
	}

	return;
}

/*ColStrGBar generate g-bar collision strengths for level 2 line2 */
STATIC double ColStrGBar(const TransitionList::iterator &t , realnum cs1 )
{
	long int i, 
	  j;
	double ColStrGBar_v, 
	  a, 
	  b, 
	  c, 
	  d, 
	  e1, 
	  gb, 
	  x, 
	  y;
	double xx, 
	  yy;

	DEBUG_ENTRY( "ColStrGBar()" );

	/* Calculation of the collision strengths of multiplets.
	 * Neutrals are recalculated from 
	 * >>refer cs	gbar	Fisher et al. (1996)
	 * >>refer cs	gbar	Gaetz & Salpeter (1983, ApJS 52, 155) and 
	 * >>refer cs	gbar	Mewe (1972, A&A 20, 215) 
	 * fits for ions. */

	/* routine to implement g-bar data taken from
	 *>>refer	cs	gbar	Mewe, R.; Gronenschild, E. H. B. M.; van den Oord, G. H. J., 1985,
	 *>>refercon	A&AS, 62, 197 */

	/* zero hydrogenic lines since they are done by iso-sequence */
	if( (*(*t).Hi()).nelem() == (*(*t).Hi()).IonStg() )
	{
		ColStrGBar_v = 0.0;
		return( ColStrGBar_v );
	}

	/*was the block data linked in? */
	ASSERT( MeweCoef.g[1][0] != 0.);

	/* which type of transition is this? cs1 is the flag */

	/* >>chng 01 may 30 - cs1 < 0 means a forced collision strength */
	if( cs1 < 0. )
	{
		ColStrGBar_v = -cs1;
		return( ColStrGBar_v );
	}

	/* >>chng 99 feb 27, change to assert */
	ASSERT( cs1 >= 0.05 );

	/* excitation energy over kT */
	y = (*t).EnergyK()/phycon.te;
	if( cs1 < 1.5 )
	{
		xx = -log10(y);

		if( cs1 < 0.5 )
		{
			yy = (1.398813573838321 + xx*(0.02943050869177121 + xx*
			  (-0.4439783893114510 + xx*(0.2316073358577902 + xx*(0.001870493481643103 + 
			  xx*(-0.008227246351067403))))))/(1.0 + xx*(-0.6064792600526370 + 
			  xx*(0.1958559534507252 + xx*(-0.02110452007196644 + 
			  xx*(0.01348743933722316 + xx*(-0.0001944731334371711))))));
		}

		else
		{
			yy = (1.359675968512206 + xx*(0.04636500015069853 + xx*
			  (-0.4491620298246676 + xx*(0.2498199231048967 + xx*(0.005053803073345794 + 
			  xx*(-0.01015647880244268))))))/(1.0 + xx*(-0.5904799485819767 + 
			  xx*(0.1877833737815317 + xx*(-0.01536634911179847 + 
			  xx*(0.01530712091180953 + xx*(-0.0001909176790831023))))));
		}

		ColStrGBar_v = exp10(yy)*(*t).Emis().gf()/((*t).EnergyRyd() * 13.6);
	}
	else
	{
		i = (long int)cs1;

		if( i < 26 )
		{
			e1 = log(1.0+1.0/y) - 0.4/POW2(y + 1.0);
			a = MeweCoef.g[i-1][0];
			b = MeweCoef.g[i-1][1];
			c = MeweCoef.g[i-1][2];
			d = MeweCoef.g[i-1][3];
			x = (double)(*(*t).Hi()).nelem() - 3.0;

			if( i == 14 )
			{
				a *= 1.0 - 0.5/x;
				b = 1.0 - 0.8/x;
				c *= 1.0 - 1.0/x;
			}

			else if( i == 16 )
			{
				a *= 1.0 - 0.9/x;
				b *= 1.0 - 1.7/x;
				c *= 1.0 - 2.1/x;
			}

			else if( i == 18 )
			{
				a *= 1.0 + 2.0/x;
				b *= 1.0 - 0.7/x;
			}

			gb = a + (b*y - c*y*y + d)*e1 + c*y;

			/*  ipLnRyd is exciation energy in Rydbergs */
			ColStrGBar_v = 14.510395*(*t).Emis().gf()*gb/((*t).EnergyRyd() );
			/* following i>=26 */
		}

		else
		{
			/* 210 is the dimem of g, so [209] is largest val */
			if( i < 210 )
			{
				j = (long)(MeweCoef.g[i-1][3]);
				if( j == 1 )
				{
					ColStrGBar_v = (*(*t).Lo()).g()*MeweCoef.g[i-1][0]*
					  pow(phycon.te/exp10((double)MeweCoef.g[i-1][2]),(double)MeweCoef.g[i-1][1]);
				}
				else
				{
					ColStrGBar_v = (*(*t).Lo()).g()*MeweCoef.g[i-1][0]*
					  sexp(MeweCoef.g[i-1][1]*(exp10((double)MeweCoef.g[i-1][2])/
					  phycon.te));
				}
			}

			else
			{
				/* This is for AlII 1670 line only!
				 *    ColStrGBar=0.0125*te**0.603 */
				/* 98 dec 27, this is still in use */
				ColStrGBar_v = 0.0125*phycon.sqrte*phycon.te10*
				  phycon.te003;
			}
		}
	}

	/* following to make sure that negative values not returned */
	ColStrGBar_v = MAX2(ColStrGBar_v,1e-10);
	return( ColStrGBar_v );
}
