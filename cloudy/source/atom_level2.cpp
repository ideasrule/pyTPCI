/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atom_level2 do level population and cooling for two level atom,
 * side effects:
 * set elements of transition struc
 * cooling via 	CoolAdd( chLab, (long)t->WLAng , t->cool());
 * cooling derivative */
#include "cddefines.h"
#include "phycon.h"
#include "transition.h"
#include "rfield.h"
#include "thermal.h"
#include "cooling.h"
#include "atoms.h"
#include "ionbal.h"
#include "dense.h"
#include "lines_service.h"
#include "ionbal.h"


/* There are two special uses of this function:
 * - With hyperfine transitions:
 * 	Then, IndirectRate are out of the ground state.  The de-excitations to
 * 	the ground state are distributed between the two levels according to
 * 	their statistical weights.
 * - With Opacity Project lines:
 *	IndirectRate is given as 0, causing the previous total to be incremented
 *	by the upward rate of OP transitions.
 * By default, the function is invoked with the indirect rate set to 0. */
void atom_level2( const TransitionProxy &t, const bool lgHFS )
{
	long int ion, 
	  ip, 
	  nelem;
	double AbunxIon, 
	  a21, 
	  boltz, 
	  col12, 
	  col21, 
	  coolng, 
	  g1, 
	  g2, 
	  omega, 
	  pfs1,
	  pfs2,
	  r, 
	  rate12, 
	  ri21;


	DEBUG_ENTRY( "atom_level2()" );

	/* result is density (cm-3) of excited state times A21
	 * result normalized to N1+N2=ABUND
	 * routine increments dCooldT, call CoolAdd
	 * CDSQTE is EDEN / SQRTE * 8.629E-6
	 */

	ion = (*t.Hi()).IonStg();
	nelem = (*t.Hi()).nelem();

	/* dense.xIonDense[nelem][i] is density of ith ionization stage (cm^-3) */
	AbunxIon = dense.xIonDense[nelem-1][ion-1];


	/* continuum pointer */
	ip = t.ipCont();

	/* approximate Boltzmann factor to see if results zero */
	boltz = rfield.ContBoltz[ip-1];


	/* statistical weights of upper and lower levels */
	g1 = (*t.Lo()).g();
	g2 = (*t.Hi()).g();


	/* indirect excitations from lower level that populate upper level */
	double IndirLU = 0.0;

	if( lgHFS )
		IndirLU = ionbal.ExcitationGround[nelem-1][ion-1] * g2 / (g1 + g2);

	/* collision strength for this transition, omega is zero for hydrogenic
	 * species which are done in special hydro routines */
	omega = t.Coll().col_str();

	/* ROUGH check whether upper level has significant population,*/
	r = (boltz*dense.cdsqte + t.Emis().pump() + IndirLU)/(dense.cdsqte + t.Emis().Aul());

	/* following first test needed for 32 bit cpu on search phase
	 * >>chng 96 jul 02, was 1e-30 crashed on dec, change to 1e-25
	 * if( AbunxIon.lt.1e-25 .or. boltz.gt.30. ) then
	 * >>chng 96 jul 11, to below, since can be strong pumping when
	 * Boltzmann factor essentially zero */
	/* omega in following is zero for hydrogenic species, since done
	 * in hydro routines, so this should cause us to quit on this test */
	/* >>chng 99 nov 29, from 1e-25 to 1e-30 to keep same result for
	 * very low density models, where AbunxIon is very small but still significant*/
	/*if( omega*AbunxIon < 1e-25 || r < 1e-25 )*/
	if( omega*AbunxIon < 1e-30 || r < 1e-25 )
	{
		/* put in pop since possible just too cool */
		(*t.Lo()).Pop() = AbunxIon;
		t.Emis().PopOpc() = AbunxIon;
		(*t.Hi()).Pop() = 0.;
		t.Emis().xIntensity() = 0.;
		t.Emis().xObsIntensity() = 0.;
		t.Coll().cool() = 0.;
		t.Emis().ots() = 0.;
		t.Emis().ColOvTot() = 0.;
		t.Coll().heat() = 0.;
		/* level populations */
		atoms.PopLevels[0] = AbunxIon;
		atoms.PopLevels[1] = 0.;
		atoms.DepLTELevels[0] = 1.;
		atoms.DepLTELevels[1] = 0.;
		return;
	}

	/* net rate down A21*(escape + destruction) */
	a21 = t.Emis().Aul()*(t.Emis().Ploss());

	/* now get real Boltzmann factor */
	boltz = t.EnergyK()/phycon.te;

	ASSERT( boltz > 0. );
	boltz = sexp(boltz);

	ASSERT( g1 > 0. && g2 > 0. );

	/* this lacks the upper statistical weight */
	col21 = dense.cdsqte*omega;
	/* upward coll rate s-1 */
	col12 = col21/g1*boltz;
	/* downward coll rate s-1 */
	col21 /= g2;


	/* rate 1 to 2 is both collisions and pumping */
	/* the total excitation rate from lower to upper, collisional and pumping */
	rate12 = col12 + t.Emis().pump() + IndirLU;

	/* include Opacity Project excitations */
	if( ! lgHFS )
		ionbal.ExcitationGround[nelem-1][ion-1] += rate12;

	/* induced emissions down */
	ri21 = (t.Emis().pump()+IndirLU)*g1/g2;

	/* this is the ratio of lower to upper level populations */
	r = (a21 + col21 + ri21)/rate12;

	/* upper level pop */
	pfs2 = AbunxIon/(r + 1.);
	atoms.PopLevels[1] = pfs2;
	(*t.Hi()).Pop() = pfs2;

	/* pop of ground */
	pfs1 = pfs2*r;
	atoms.PopLevels[0] = pfs1;

	/* compute ratio Aul/(Aul+Cul) */
	/* level population with correction for stimulated emission */
	(*t.Lo()).Pop() = atoms.PopLevels[0];


	t.Emis().PopOpc() = (atoms.PopLevels[0] - atoms.PopLevels[1]*g1/g2 );

	/* departure coef of excited state rel to ground */
	atoms.DepLTELevels[0] = 1.;
	if( boltz > 1e-20 && atoms.PopLevels[1] > 1e-20 )
	{
		/* this line could have zero boltz factor but radiatively excited
		 * dec alpha does not obey () in fast mode?? */
		atoms.DepLTELevels[1] = (atoms.PopLevels[1]/atoms.PopLevels[0])/
		  (boltz*g2/g1);
	}
	else
	{
		atoms.DepLTELevels[1] = 0.;
	}


	/* number of escaping line photons, used elsewhere for outward beam
	 * and line intensity */
	set_xIntensity( t );


	/* ratio of collisional to total (collisional + pumped) excitation */
	t.Emis().ColOvTot() = (col12 + IndirLU) /rate12;


	/* two cases - collisionally excited (usual case) or 
	 * radiatively excited - in which case line can be a heat source
	 * following are correct heat exchange, they will mix to get correct deriv 
	 * the sum of heat-cool will add up to EnrUL - EnrLU - this is a trick to
	 * keep stable solution by effectively dividing up heating and cooling,
	 * so that negative cooling does not occur */

	//double Enr12 = plower*col12*t.EnergyErg;
	//double Enr21 = pfs2*col21*t.EnergyErg;

	/* energy exchange due to this level
	 * net cooling due to excit minus part of de-excit -
	 * note that ColOvTot cancels out in the sum heat - cool */
	//coolng = Enr12 - Enr21*t.Emis().ColOvTot();

	/* this form of coolng is guaranteed to be non-negative */
	//coolng = t.EnergyErg()*AbunxIon*col12*(a21 + ri21)/(a21 + col21 + ri21 + rate12);
	coolng = t.EnergyErg()*(pfs1*col12-pfs2*col21);
	//ASSERT( coolng >= 0. );

	t.Coll().cool() = MAX2(0.,coolng);

	/* net heating is remainder */
	//t.Coll().heat() = t.EnergyErg()*AbunxIon*col21*(t.Emis().pump())/(a21 + col21 + ri21 + rate12);
	t.Coll().heat() = MAX2(0.,-coolng);

	/* expression pre jul 3 95, changed for case where line heating dominates
	 * coolng = (plower*col12 - pfs2*col21)*t.t(ipLnEnrErg)
	 * t.t(ipLnCool) = cooling */

	/* add to cooling - heating part was taken out above,
	 * and is not added in here - it will be added to thermal.heating(0,22)
	 * in CoolSum */
	CoolAdd( chIonLbl(t).c_str(), t.WLAng() , t.Coll().cool());
	thermal.elementcool[nelem-1] += MAX2( 0., t.Coll().cool() );

	/* derivative of cooling function */
	thermal.dCooldT += coolng * (t.EnergyK() * thermal.tsq1 - thermal.halfte );
	return;
}



void atom_level2( const TransitionProxy &t )
{
	atom_level2( t, false );

	return;
}
