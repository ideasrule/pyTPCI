/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ChargTranEval fill in the CharExcIonOf[ipHYDROGEN] and Rec arrays with Kingdon's fitted CT with H */
/*ChargTranSumHeat sum net heating due to charge transfer, called in HeatSum */
/*HCTIon H charge transfer ionization*/ 
/*HCTRecom H charge transfer recombination */
/*ChargTranPun, save charge transfer coefficient */
/*MakeHCTData holds data for charge transfer fits */
#include "cddefines.h"
#include "phycon.h"
#include "abund.h"
#include "dense.h"
#include "iso.h"
#include "thermal.h"
#include "elementnames.h"
#include "heavy.h"
#include "trace.h"
#include "conv.h"
#include "atmdat.h"
#include "ion_trim.h"
#include "mole.h"
#include "parser.h"

/*HCTion H charge transfer ionization, H+ + A => H + A+ */
STATIC double HCTIon(long int ion, 
  long int nelem);

/*HCTRecom H charge transfer recombination, H + A+ => H+ + A */
STATIC double HCTRecom(long int ion, 
  long int nelem);

/* the translated block data */
STATIC void MakeHCTData(void);

/* the structures for storing the charge transfer data, these are filled in
 * at the end of this file, in what used to be a block data  */
static double CTIonData[LIMELM][4][8];
static double CTRecombData[LIMELM][4][7];
/* this will be flag for whether or not charge transfer data
 * have been initialized */
static bool lgCTDataDefined = false;

/*ChargTranEval update charge trans rate coefficients if temperature has changed */
void ChargTranEval( void )
{
	long int ion, 
	  nelem;
	double a, b, c, a_op, b_op, c_op, d_op, e_op, f_op, a_o, 
			b_o, c_o, d_o, e_o, f_o, g_o; 

	static double TeUsed = -1.;

	DEBUG_ENTRY( "ChargTranEval()" );

	/* first is to force reevaluation on very first call */
	if( !conv.nTotalIoniz || !fp_equal(phycon.te,TeUsed) )
	{
		/* refresh hydrogen charge transfer arrays */
		/* >>chng 01 apr 25, lower limit had been 2, lithium, changed to 1, helium */
		for( nelem=ipHELIUM; nelem < LIMELM; nelem++ )
		{
			for( ion=0; ion <= nelem; ion++ )
			{
				/* >>chng 01 apr 28, add factors  to turn off ct,
				 * had previously been handled downstream */

				/* CharExcIonOf[ipHYDROGEN][nelem][ion]*hii  is ion => ion+1 for nelem */
				/* charge transfer ionization O^0 + H+ -> O+ + H0 
				 * is CharExcIonOf[ipHYDROGEN][ipOXYGEN][0]*dense.xIonDense[ipHYDROGEN][1]
				 * charge transfer recombination of atomic hydrogen is
				 * CharExcIonOf[ipHYDROGEN][ipOXYGEN][0]*dense.xIonDense[ipOXYGEN][0] */
				atmdat.CharExcIonOf[ipHYDROGEN][nelem][ion] = HCTIon(ion,nelem);

				/* CharExcRecTo[ipHYDROGEN][nelem][ion]*hi is ion+1 => ion of nelem */
				/* charge transfer recombination O+ + H0 -> O^0 + H+ is
				 * CharExcRecTo[ipHYDROGEN][ipOXYGEN][0]*dense.xIonDense[ipHYDROGEN][0]
				 * charge transfer ionization of atomic hydrogen is
				 * CharExcRecTo[ipHYDROGEN][ipOXYGEN][0]*dense.xIonDense[ipOXYGEN][1] */
				atmdat.CharExcRecTo[ipHYDROGEN][nelem][ion] = HCTRecom(ion,nelem);
			}

			/* zero out helium charge transfer arrays */
			for( ion=0; ion < LIMELM; ion++ )
			{
				atmdat.CharExcIonOf[ipHELIUM][nelem][ion] = 0.;
				atmdat.CharExcRecTo[ipHELIUM][nelem][ion] = 0.;
			}
		}

		/* The above included only the radiative charge transfer from
		 * Stancil et al 1998.  must explicitly add on the ct fitted by Kingdon & Ferland,
		 * The process H0 + He+ -> He0 + H+ */
		atmdat.CharExcRecTo[ipHYDROGEN][ipHELIUM][0] += 7.47e-15*pow(phycon.te/1e4,2.06)*
			(1.+9.93*sexp(3.89e-4*phycon.te) );

		/* >>chng 04 jun 21 -- NPA.  Put in the charge transfer rate for:
			He+ + H => He + H+ as defined in the UMIST database.  This is only
			used if the "set UMIST rates" command is used */
		if(mole_global.lgLeidenHack)
		{
			atmdat.CharExcRecTo[ipHYDROGEN][ipHELIUM][0] = 4.85e-15*pow(phycon.te/300, 0.18);
		}
		/* >>chng 04 jun 21 -- NPA.  update the charge transfer rates between hydrogen
		   and oxygen to:  
		   >>refer	O	CT	Stancil, et al. 1999, A&AS, 140, 225-234 
		   Instead of using the UMIST rate, the program TableCurve was used
		   to generate a fit to the data listed in Tables 2, 3, and 4 of the 
		   aforementioned reference. The following fitting equations agree 
		   very well with the published data. */

		/* At or below 10K, just use the value the formula's below give
		   at 10K.*/
		/* do both O+ -> O and O -> O+ for low T limit */
		if(phycon.te <= 10. )
		{
			atmdat.CharExcRecTo[ipHYDROGEN][ipOXYGEN][0] = 3.744e-10;
			atmdat.CharExcIonOf[ipHYDROGEN][ipOXYGEN][0] = 4.749e-20;
		}

		/* this does O+ -> O for all higher temperatures */
		if( phycon.te > 10.)
		{
			a_op = 2.3344302e-10;
			b_op = 2.3651505e-10;
			c_op = -1.3146803e-10;
			d_op = 2.9979994e-11;
			e_op = -2.8577012e-12;
			f_op = 1.1963502e-13;

			double lnte = phycon.alnte;

			/* equation rank 53 of TableCurve */
			atmdat.CharExcRecTo[ipHYDROGEN][ipOXYGEN][0] =
				((((f_op*lnte + e_op)*lnte + d_op)*lnte + c_op)*lnte + b_op)*lnte + a_op;
		}

		/* now do O -> O+
		 * The next two equations were chosen to match up at 200K, so that there
		 *are no discontinuities */
		if((phycon.te > 10.) && (phycon.te <= 190.))
		{
			a = -21.134531;
			b = -242.06831;
			c = 84.761441;

			/* equation rank 2 of TableCurve */
			atmdat.CharExcIonOf[ipHYDROGEN][ipOXYGEN][0] = exp((a + b/SDIV(phycon.te) + c/SDIV(phycon.tesqrd)));
		}

		else if((phycon.te > 190.) && (phycon.te <= 200.))
		{

			/* We are using two fitting formula's for this rate, in order to assure no
			sudden "jumps" in the rate, the rate between 190-200K is made to 
			increase linearly.  The formula below gets the same answer as the equation
			above at 190K, and gets the same answer as the the formula below this one 
			at 200K*/
			atmdat.CharExcIonOf[ipHYDROGEN][ipOXYGEN][0] = 2.18733e-12*(phycon.te-190) + 1.85823e-10;
		}

		else if(phycon.te > 200.)
		{

			a_o = -7.6767404e-14;
			b_o = -3.7282001e-13;
			c_o = -1.488594e-12;
			d_o = -3.6606214e-12; 
			e_o = 2.0699463e-12;
			f_o = -2.6139493e-13;
			g_o = 1.1580844e-14;

			double lnte = phycon.alnte;

			/* equation rank 120 of TableCurve */
			atmdat.CharExcIonOf[ipHYDROGEN][ipOXYGEN][0] = 
				(((((g_o*lnte + f_o)*lnte + e_o)*lnte + d_o)*lnte + c_o)*lnte + b_o)*lnte + a_o;
		}

		/* use UMIST rates for charge transfer if UMIST command is used - disagree
		 * by 20% at 5000K and diverge at high temperature */
		if(mole_global.lgLeidenHack)
		{
			atmdat.CharExcIonOf[ipHYDROGEN][ipOXYGEN][0] = hmrate4(7.31e-10,0.23,225.9,phycon.te);
			atmdat.CharExcRecTo[ipHYDROGEN][ipOXYGEN][0] = hmrate4(5.66e-10,0.36,-8.60,phycon.te);
		}

		/* >>chng 01 may 07, following had all been +=, ok if above was zero.
		 * changed to = and added HCTOn */
		/* >>chng 01 jan 30, add following block of CT reactions */
		/* ionization, added as per Phillip Stancil OK, number comes from 
		 * >>refer	Fe	CT	Tielens, A. G. G. M., & Hollenbach, D. 1985a, ApJ, 294, 722-746 */
		/* >>refer	Fe	CT	Prasad, S. S., & Huntress, W. T. 1980, ApJS, 43, 1-35 */
		/* the actual rate comes from the following paper: */
		/* >>refer	Fe	CT	Pequignot, D., & Aldrovandi, S. M. V. 1986, A&A, 161, 169-176 */
		/* Fe0 + H+ => Fe+ + H0 */
		/*>>chng 05 sep 15, GS, old rate had problem in predicting observed Fe I column density along HD185418.
		 *>> refer Private communication with Stancil, data taken from ORNL web site,
		 * "There is a well known problem with the Fe charge transfer rate  coefficients: i.e., there are no accurate calculations nor or there
			any experiments. For Fe + H+ -> Fe+ + H, I estimated for Gary a few  years ago the value of 5.4e-9. So mid way between the two values
			you are using. I have some notes on it in my office, but not with me.  See: http://cfadc.phy.ornl.gov/astro/ps/data/home.html
			 value changed from 3e-9 to 5.4e-9 */

		atmdat.CharExcIonOf[ipHYDROGEN][ipIRON][0] = 5.4e-9;
		/*>>chng 06 sep 20 - following sets removes Fe ionization ct to be similar to Mg */
		/*atmdat.CharExcIonOf[ipHYDROGEN][ipIRON][0] = 0.;broken();rm this line */

		/* all remaining entries are from Pequignot & Aldrovandi*/
		/* >>refer	Al	CT	Pequignot, D., & Aldrovandi, S. M. V. 1986, A&A, 161, 169-176 */
		/* Al0 + H+ => Al+ + H0 */
		atmdat.CharExcIonOf[ipHYDROGEN][ipALUMINIUM][0] = 3e-9;

		/* >>refer	P	CT	Pequignot, D., & Aldrovandi, S. M. V. 1986, A&A, 161, 169-176 */
		/* P0 + H+ => P+ + H0 */
		atmdat.CharExcIonOf[ipHYDROGEN][ipPHOSPHORUS][0] = 1e-9;

		/* >>refer	Cl	CT	Pequignot, D., & Aldrovandi, S. M. V. 1986, A&A, 161, 169-176 */
		/* Cl0 + H+ => Cl+ + H0 */
		atmdat.CharExcIonOf[ipHYDROGEN][ipCHLORINE][0] = 1e-9;

		/* >>refer	Ti	CT	Pequignot, D., & Aldrovandi, S. M. V. 1986, A&A, 161, 169-176 */
		/* Ti0 + H+ => Ti+ + H0 */
		atmdat.CharExcIonOf[ipHYDROGEN][ipTITANIUM][0] = 3e-9;

		/* >>refer	Mn	CT	Pequignot, D., & Aldrovandi, S. M. V. 1986, A&A, 161, 169-176 */
		/* Mn0 + H+ => Mn+ + H0 */
		atmdat.CharExcIonOf[ipHYDROGEN][ipMANGANESE][0] = 3e-9;

		/* >>refer	Ni	CT	Pequignot, D., & Aldrovandi, S. M. V. 1986, A&A, 161, 169-176 */
		/* Ni0 + H+ => Ni+ + H0 */
		atmdat.CharExcIonOf[ipHYDROGEN][ipNICKEL][0] = 3e-9;

		/* >>chng 01 feb 02, add following CT reaction from */
		/* >>refer	Na0	CT	Dutta, C. M., Nordlander, P., Kimura, M., & Dalgarno, A. 2001, Phys. Rev. A, 63, 022709 */
		/* this is roughly their number around 500K - they do not give explicit values, rather
		 * a small figure.  Previous calculations were 5 orders of mag smaller at this temp.  
		 * ND this deposits electron into n=2 */
		/* Na0 + H+ => Na+ + H0(n=2) */
		atmdat.CharExcIonOf[ipHYDROGEN][ipSODIUM][0] = 7e-12;

		/* >>chng 05 sep 15,GS, add following CT reaction from */
		/* >>refer	Na0	CT	Watanabe, A., Dutta, C. M., Nordlander, P., et al. 2002, Phys. Rev. A, 66, 044701 */
		/* this is roughly their number around 50K - they do not give explicit values, rather
		 * a small figure. this deposits electron into n=1
		 * Na0 + H+ => Na+ + H0(n=1) 
		 * add to previous rate which was for population of n=2 */
		atmdat.CharExcIonOf[ipHYDROGEN][ipSODIUM][0] += 0.7e-12;

		/* >>chng 05 sep 15, GS, add following CT reaction from 
		 * >>refer	K0	CT	Watanabe, A., Dutta, C. M., Nordlander, P., et al. 2002, Phys. Rev. A, 66, 044701 
		 * this is roughly their number around 50K - they do not give explicit values, rather
		 * a small figure. 
		 * K0 + H+ => K+ + H0(n=1) */
		atmdat.CharExcIonOf[ipHYDROGEN][ipPOTASSIUM][0] = 1.25e-12;

		/* >>chng 05 sep 15, GS, add following CT reaction from 
		 * >>refer	S0	CT	ORNL data base for charge transfer	
		 * This rate is valid for 1e3 to 1e4. Due to the small value, I did not put any limit on temp. 
		 * Earlier, other reactions also assume constant value
		 * S0 + H+ => H + S+ */
		atmdat.CharExcIonOf[ipHYDROGEN][ipSULPHUR][0] = 1.e-14;

		if( phycon.te < 1e5 )
		{

			/* >>chng 05 sep 15, GS, add following CT reaction from 
			 * >>refer	Mg0	CT	ORNL data base for charge transfer,
			 * this rate is valid for temp 5e3 to 3e4, The rate goes down very fast in low temp. So I did not put a lower cut of for temp	
			 * Mg0 + H+ => H + Mg+ */
			atmdat.CharExcIonOf[ipHYDROGEN][ipMAGNESIUM][0] = 9.76e-12*pow((phycon.te/1e4),3.14)*(1. + 55.54*sexp(1.12*phycon.te/1e4));
			/*>>chng 06 jul 20, do not allow this to fall below UMIST rate - above fit not intended for 
			 * very low temperatures */
			/*>>chng 06 aug 01, UMIST is bogus - email exchange with Phillip Stancil, late July 2006 */
			/*atmdat.CharExcIonOf[ipHYDROGEN][ipMAGNESIUM][0] = MAX2( 1.1e-9 , atmdat.CharExcIonOf[ipHYDROGEN][ipMAGNESIUM][0]);*/
			/*>>chng 06 sep 20 - following sets Mg ionization ct to Fe */
			/*atmdat.CharExcIonOf[ipHYDROGEN][ipMAGNESIUM][0] = 5.4e-9;broken(); rm this line */

			/* >>chng 05 sep 15, GS, add following CT reaction from 
			 * >>refer	Si0	CT	ORNL data base for charge transfer
			 * this rate is valid for temp 1e3 to 2e5, The rate goes down very fast in low temp. So I did not put a lower cut of for temp
			 * Si0 + H+ => H + Si+ */
			atmdat.CharExcIonOf[ipHYDROGEN][ipSILICON][0] = 0.92e-12*pow((phycon.te/1e4),1.15)*(1. + 0.80*sexp(0.24*phycon.te/1e4));
			/*>>chng 06 jul 20, do not allow this to fall below UMIST rate - above fit not intended for 
			 * very low temperatures */
			/** \todo	1	update ct to Kimura et al. (1996) */
			/*>>chng 06 aug 01, UMIST rate is bogus as per Phillip Stancil emails of late July 2006 */
			/*atmdat.CharExcIonOf[ipHYDROGEN][ipSILICON][0] = MAX2( 9.9e-10 , atmdat.CharExcIonOf[ipHYDROGEN][ipSILICON][0]);*/

			/* >>chng 05 sep 15, GS, add following CT reaction from 
			 * >>refer	Li0	CT	ORNL data base for charge transfer
			 * this rate is valid for temp 1e2 to 1e4, The rate goes down very fast in low temp. So I did not put a lower cut of for temp
			 * Li0 + H+ => H + Li+ */
			atmdat.CharExcIonOf[ipHYDROGEN][ipLITHIUM][0] = 2.84e-12*pow((phycon.te/1e4),1.99)*(1. + 375.54*sexp(54.07*phycon.te/1e4));
			/** \todo	1	above rate not intended for very low temperatures - find ref for low-T rate,
			 * probably is 1e-9 like above */
		}
		else
		{
			/** \todo	0	these should be values at 1e5 K */
			atmdat.CharExcIonOf[ipHYDROGEN][ipMAGNESIUM][0] = 0.;
			atmdat.CharExcIonOf[ipHYDROGEN][ipSILICON][0] = 0.;
			atmdat.CharExcIonOf[ipHYDROGEN][ipLITHIUM][0] = 0.;
		}

		{
			/*>>chng 06 jul 07, Terry Yun add these charge transfer reactions */
			/*>>refer	N0	CT	Lin, C. Y., Stancil, P. C., Gu, J. P., et al. 2005, Phys. Rev. A, 71, 062708 
			 * and combined with data from 
			 *>>refer	N0	CT	Butler, S. E., & Dalgarno, A. 1979, ApJ, 234, 765 */

			/* natural log of te */
			double tefac = phycon.te * phycon.alnte;

			/* N(4S) + H+ -> N+(3P) + H */
			/* >>chng 06 jul 10, add exp for endoergic reaction */
			double ct_from_n0grhp_to_npgrh0 = (1.64e-16*phycon.te - 8.76e-17*tefac + 2.41e-20*phycon.tesqrd + 9.83e-13*phycon.alnte )*
				sexp( 10985./phycon.te ) * atmdat.lgCTOn;

			/* N(2D) + H+ -> N+(3P) + H */
			/** \todo	2	not currently used - include as deexcitation process */
			/*double ct_from_n0exhp_to_npgrh0 = 1.51e-15*phycon.te -1.61e-16*tefac + 7.74e-21*phycon.tesqrd + 1.34e-16*phycon.alnte;*/

			/* N+(3P) + H -> N(4S) + H+ endoergic */
			double ct_from_npgrh0_to_n0grhp = (1.56e-15*phycon.te - 1.79e-16*tefac + 1.15e-20*phycon.tesqrd + 1.08e-13*phycon.alnte) * atmdat.lgCTOn;

			/* N+(3P) + H0 -> N(2D) + H+ */
			/* >>chng 06 jul 10, add exp for endoergic reaction */
			atmdat.HCharExcRecTo_N0_2D = (6.83e-16*phycon.te - 7.40e-17*tefac + 3.73e-21*phycon.tesqrd + 1.75e-15*phycon.alnte)*
				 sexp( 16680./phycon.te ) * atmdat.lgCTOn;

			/* these rates are from the ground state into all possible states of the
			 * species that is produced */
			atmdat.CharExcIonOf[ipHYDROGEN][ipNITROGEN][0] = ct_from_n0grhp_to_npgrh0;
			atmdat.CharExcRecTo[ipHYDROGEN][ipNITROGEN][0] = ct_from_npgrh0_to_n0grhp + atmdat.HCharExcRecTo_N0_2D;
		}

		/*>>chng 06 aug 01, update O++ and N++ -- H0 CT recombination 
		 *>>refer	O3	CT	Barragan, P., Errea, L. F., Mendez, L., et al. 2006, ApJ, 636, 544 */
		/* O+2 + H -> O+ + H+ */
		if( phycon.te <= 1500. )
		{
			atmdat.CharExcRecTo[ipHYDROGEN][ipOXYGEN][1] = 0.5337e-9*pow( (phycon.te/100.) ,-0.076);
		}
		else
		{
			atmdat.CharExcRecTo[ipHYDROGEN][ipOXYGEN][1] = 0.4344e-9 +
				0.6340e-9*pow( log10(phycon.te/1500.) ,2.06 );
		}

		/* N+2 + H -> N+ + H+ */
		if( phycon.te <= 1500. )
		{
			atmdat.CharExcRecTo[ipHYDROGEN][ipNITROGEN][1] = 0.8692e-9*pow( (phycon.te/1500.) ,0.17);
		}
		else if( phycon.te <= 20000. )
		{
			atmdat.CharExcRecTo[ipHYDROGEN][ipNITROGEN][1] = 0.9703e-9*pow( (phycon.te/10000.) ,0.058);
		}
		else
		{
			atmdat.CharExcRecTo[ipHYDROGEN][ipNITROGEN][1] = 1.0101e-9 +
				1.4589e-9*pow( log10(phycon.te/20000.) ,2.06 );
		}

		/* ===================== helium charge transfer ====================*/

		/* atmdat.CharExcIonOf[ipHELIUM] is ionization, */
		/* [0] is Atom^0 + He+ => Atom+1 + He0
		 * [n] is Atom^+n + He+ => Atom^+n-1 + He0 */

		/* atmdat.CharExcRecTo[ipHELIUM] is recombination */
		/* [0] is Atom^+1 + He0 => Atom^0 + He^+
		 * [n] is Atom^+n+1 + He0 => Atom^+n + He^+ */

		/* Carbon */
		/* recombination */
		/* C+3 + He => C+2 + He+ */
		/* >>refer	C3	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipCARBON][2] = 4.6e-19*phycon.tesqrd;

		/* C+4 + He => C+3 + He+ */
		/* >>refer	C4	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipCARBON][3] = 1e-14;

		/* ionization */
		/* C0 + He+ => C+ + He0 */
		/* unknown reference, from older version of the code */
		/*atmdat.CharExcIonOf[ipHELIUM][ipCARBON][0] = 4.17e-17*(phycon.te/phycon.te03);*/

		/* >>chng 04 jun 21 -- update this rate to that given in the UMIST database - NPA */
		atmdat.CharExcIonOf[ipHELIUM][ipCARBON][0] = 6.3e-15*pow((phycon.te/300),0.75);

		/* C+1 + He+ => C+2 + He */
		/* >>refer	C1	CT	Butler, S. E., Heil, T. G., & Dalgarno, A. 1980, ApJ, 241, 442*/
		atmdat.CharExcIonOf[ipHELIUM][ipCARBON][1] = 
			5e-20*phycon.tesqrd*sexp(0.07e-4*phycon.te)*sexp(6.29/phycon.te_eV);

		/* nitrogen */
		/* recombination */
		/* N+2 => N+ Butler and Dalgarno 1980B
		 * ct with update
		 * >>refer	N2	CT	Sun, Y., Sadeghpour, H.R., Kirby, K., et al. CfA preprint 4208
		 * this agrees with exp 
		 * >>refer	N2	CT	Fang, Z., & Kwong, V. H. S. 1997, ApJ, 474, 529 */
		atmdat.CharExcRecTo[ipHELIUM][ipNITROGEN][1] = 0.8e-10;

		/* N+3 => N+2 */
		/* >>refer	N3	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipNITROGEN][2] = 1.5e-10;

		/* ce rate from quantal calc of ce,
		 * >>refer	N4	CT	Feickert, C. A., Blint, R. J., Surratt, G. T., & Watson, W.D. 1984, ApJ, 286, 371,
		 * >>refer	N4	CT	Rittby, M., Elander, N., Brandas, E., & Barany, A. 1984, J. Phys. B, 17, L677.
		 * CR = 1.E-9 + 8E-12 * TE10 * SQRTE */
		atmdat.CharExcRecTo[ipHELIUM][ipNITROGEN][3] = 2e-9;

		/* ionization */
		/* N+1 + He+ => N+2 + He */
		/* >>refer	N1	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcIonOf[ipHELIUM][ipNITROGEN][1] = 
			3.7e-20*phycon.tesqrd*sexp(0.063e-4*phycon.te)*sexp(1.44/phycon.te_eV);

		/* oxygen */
		/* recombination */
		/* O+2 + He  => O+1 + He+ */
		/* >>refer	O2	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipOXYGEN][1] = 3.2e-14*phycon.te/phycon.te05;
		/* O+3 + He => O+2 + He+ */
		/* >>refer	O3	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipOXYGEN][2] = 1e-9;
		/* O+4 + He  => O+3 + He+ */
		/* >>refer	O4	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipOXYGEN][3] = 6e-10;

		/* ionization */
		/* O0 + He+ => O+ + He0 */
		/* >>refer	O0	CT	Zhao, L. B., Stancil, P. C., Gu, J. P., et al. 2004, ApJ, 615, 1063 */
		atmdat.CharExcIonOf[ipHELIUM][ipOXYGEN][0] = 
			4.991E-15 * pow( phycon.te / 1e4, 0.3794 )* sexp( phycon.te/1.121E6 ) +
			2.780E-15 * pow( phycon.te / 1e4, -0.2163 )* exp( -1. * MIN2(1e7, phycon.te)/(-8.158E5) );

		/* neon */
		/* recombination */
		/* Ne+2 + He  => Ne+1 + He+ */
		/* >>refer	Ne2	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipNEON][1] = 1e-14;
		/* Ne+3 + He => Ne+2 + He+ */
		/* >>refer	Ne3	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipNEON][2] = 1e-16*phycon.sqrte;
		/* Ne+4 + He  => Ne+3 + He+ */
		/* >>refer	Ne4	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipNEON][3] = 1.7e-11*phycon.sqrte;

		/* magnesium */
		/* recombination */
		/* Mg+3 + Heo => Mg+2 */
		/* >>refer	Mg3	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipMAGNESIUM][2] = 7.5e-10;
		/* Mg+4 + Heo => Mg+3 */
		/* >>refer	Mg4	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipMAGNESIUM][3] = 1.4e-10*phycon.te30;


		/* silicon */
		/* recombination */
		/* Si+3 +He => Si+2 + He+ */
		/* >>refer	Si3	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 
		 * scale by 1.3 to bring into agreement with
		 * >>refer	Si3	CT	Fang, Z., & Kwong, V. H. S. 1997, ApJ, 483, 527 */
		atmdat.CharExcRecTo[ipHELIUM][ipSILICON][2] += phycon.sqrte*phycon.te10*phycon.te10*
		  1.3*1.5e-12;

		/* Si+4 + Heo => Si+3
		 * >>refer	Si4	CT	Opradolce, L., McCarrol, R., & Valiron, P. 1985, A&A, 148, 229 */
		atmdat.CharExcRecTo[ipHELIUM][ipSILICON][3] = 2.54e-11*phycon.sqrte/phycon.te03/
		  phycon.te01/phycon.te01;

		/* ionization */
		/* Si0 + He+ => Si+ + He0 */
		/* >>refer	Si0	CT	Prasad, S. S., & Huntress, W. T. 1980, ApJS, 43, 1-35 */
		atmdat.CharExcIonOf[ipHELIUM][ipSILICON][0] = 3.3e-9;

		/* Si+1 + He+ => Si+2 + He */
		/* >>refer	Si1	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcIonOf[ipHELIUM][ipSILICON][1] = 
			1.5e-11*phycon.te20*phycon.te05*sexp(6.91/phycon.te_eV);

		/* Si+2 + He+ => Si+3 + He */
		/* >>refer	Si2	CT	Gargaud, M., McCarroll, R., & Valiron, P. 1982, A&AS, 45, 603 */
		atmdat.CharExcIonOf[ipHELIUM][ipSILICON][2] = 
			1.15e-11*phycon.sqrte*sexp(8.88/phycon.te_eV);

		/* sulphur */
		/* recombination */
		/* S+3 + Heo => S+2 */
		/* >>refer	S3	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipSULPHUR][2] = phycon.sqrte*1.1e-11;

		/* S+4 + Heo => S+3 */
		/* >>refer	S4	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		/* >>chng 04 jul 01, from [ipSULPHUR][2] to [3] - bug */
		atmdat.CharExcRecTo[ipHELIUM][ipSULPHUR][3] = 4.8e-14*phycon.te30;

		/* ionization */
		/* S+1 + He+ => S+2 + He */
		/* >>refer	S1	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcIonOf[ipHELIUM][ipSULPHUR][1] = 
			4.4e-16*phycon.te*phycon.te20*sexp(0.036e-4*phycon.te)*sexp(9.2/phycon.te_eV);

		/* S+2 + He+ => S+3 + He */
		/* >>refer	S2	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcIonOf[ipHELIUM][ipSULPHUR][2] = 
			5.5e-18*phycon.te*phycon.sqrte*phycon.te10*sexp(0.046e-4*phycon.te)*sexp(10.5/phycon.te_eV);

		/* Argon */
		/* recombination */
		/* >>refer	Ar2	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipARGON][1] = 1.3e-10;

		/* >>refer	Ar3	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipARGON][2] = 1.e-14;

		/* >>refer	Ar4	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcRecTo[ipHELIUM][ipARGON][3] = 1.6e-8/phycon.te30;

		/* ionization */
		/* Ar+1 + He+ => Ar+2 + He0 */
		/* >>refer	Ar1	CT	Butler, S. E., & Dalgarno, A. 1980b, ApJ, 241, 838 */
		atmdat.CharExcIonOf[ipHELIUM][ipARGON][1] = 1.1e-10;

		TeUsed = phycon.te;

		if(mole_global.lgLeidenHack)
		{
			/* Set all charge transfer rates equal to zero that do not appear
			   in the UMIST database.  This if statement is only performed
			   if the "set UMIST rates" command is used */

			atmdat.CharExcIonOf[ipHYDROGEN][ipHELIUM][0] = 0;
			atmdat.CharExcIonOf[ipHYDROGEN][ipCARBON][0] = 0;
			atmdat.CharExcRecTo[ipHYDROGEN][ipCARBON][0] = 0;

			atmdat.CharExcRecTo[ipHELIUM][ipCARBON][0] = 0;
			atmdat.CharExcIonOf[ipHELIUM][ipOXYGEN][0] = 0;
			atmdat.CharExcRecTo[ipHELIUM][ipOXYGEN][0] = 0;
		}


		/* this is set false with the no charge transfer command */
		if( !atmdat.lgCTOn )
		{
			for( nelem=0; nelem< LIMELM; ++nelem )
			{
				for( ion=0; ion<LIMELM; ++ion )
				{
					atmdat.CharExcIonOf[ipHYDROGEN][nelem][ion] = 0.;
					atmdat.CharExcRecTo[ipHYDROGEN][nelem][ion] = 0.;
					atmdat.CharExcIonOf[ipHELIUM][nelem][ion] = 0.;
					atmdat.CharExcRecTo[ipHELIUM][nelem][ion] = 0.;
				}
			}
		}
	}

	return;
}

/*================================================================================*
 *================================================================================*/
double ChargTranSumHeat(void)
{
	long int ion, 
	  nelem;
	double SumCTHeat_v;

	DEBUG_ENTRY( "ChargTranSumHeat()" );

	/* second dimension is ionization stage,
	 * 1=+1 for parent, etc
	 * third dimension is atomic weight of atom */

	/* make sure data are initialized */
	ASSERT( lgCTDataDefined );

	SumCTHeat_v = 0.;
	/* >>chng 01 apr 25, lower limit had been 0 should have been 1 (helium) */
	for( nelem=ipHELIUM; nelem < LIMELM; nelem++ )
	{
		/* >>chng >>01 apr 25, loops had been to LIMELM, which may have done no harm
		 * since extra array elements were set to zero, but is incorrect since the physical
		 * limit is the number of stages of ionization */
		int limit = MIN2(4, nelem);
		/* this first group of lower stages of ionization have exact rate coefficients */
		for( ion=0; ion < limit; ion++ )
		{
			/* CTIonData[nelem][ion][7] and CTRecombData[nelem][ion][6] are the energy deficits in eV,
			 * atmdat.CharExcIonOf[ipHYDROGEN][nelem][ion] and atmdat.CharExcIonOf[ipHYDROGEN][nelem][ion] 
			 * save the rate coefficients 
			 * this is sum of heat exchange in eV s^-1 cm^-3 */
			SumCTHeat_v += 

				/* heating due to ionization of heavy element, recombination of hydrogen */
				atmdat.CharExcIonOf[ipHYDROGEN][nelem][ion]*CTIonData[nelem][ion][7]*
				dense.xIonDense[ipHYDROGEN][1]*
				dense.xIonDense[nelem][ion] + 

				/* heating due to recombination of heavy element, ionization of hydrogen */
				atmdat.CharExcRecTo[ipHYDROGEN][nelem][ion]*CTRecombData[nelem][ion][6]*
				//iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop()*
				dense.xIonDense[ipHYDROGEN][0]*
				dense.xIonDense[nelem][ion+1];
		}

		/* >>chng >>01 apr 25, following loop had been to LIMELM, change to nelem */
		/* following do not have exact energies, so use 2.86*(Z-1) */
		for( ion=4; ion < nelem; ion++ )
		{
			SumCTHeat_v += 
				atmdat.CharExcRecTo[ipHYDROGEN][nelem][ion]* 2.86*(double)ion *
				dense.xIonDense[ipHYDROGEN][0]*
				dense.xIonDense[nelem][ion+1];
		}
	}

#if 0
	/* charge exchange with helium */
	for( nelem=ipLITHIUM; nelem < LIMELM; nelem++ )
	{
		/* >>chng >>01 apr 25, loops had been to LIMELM, which may have done no harm
		 * since extra array elements were set to zero, but is incorrect since the physical
		 * limit is the number of stages of ionization */
		int limit = MIN2(4, nelem);
		/* this first group of lower stages of ionization have exact rate coefficients */
		for( ion=0; ion < limit; ion++ )
		{
			fixit(); // fill in these barriers
			double barrier_rec_eV = CTIonData[nelem][ion][7];
			double barrier_ion_eV = CTRecombData[nelem][ion][6];
			
			/* atmdat.CharExcIonOf[ipHELIUM][nelem][ion] and atmdat.CharExcIonOf[ipHELIUM][nelem][ion] 
			 * save the rate coefficients 
			 * this is sum of heat exchange in eV s^-1 cm^-3 */
			SumCTHeat_v += 

				/* heating due to ionization of heavy element, recombination of helium */
				atmdat.CharExcIonOf[ipHELIUM][nelem][ion]*barrier_rec_eV*
				dense.xIonDense[ipHELIUM][1]*
				dense.xIonDense[nelem][ion] + 

				/* heating due to recombination of heavy element, ionization of helium */
				atmdat.CharExcRecTo[ipHELIUM][nelem][ion]*barrier_ion_eV*
				//iso_sp[ipHE_LIKE][ipHELIUM].st[ipH1s].Pop()*
				dense.xIonDense[ipHELIUM][0]*
				dense.xIonDense[nelem][ion+1];
		}

		/* >>chng >>01 apr 25, following loop had been to LIMELM, change to nelem */
		/* following do not have exact energies, so use 2.86*(Z-1) */
		for( ion=4; ion < nelem; ion++ )
		{
			SumCTHeat_v += 
				atmdat.CharExcRecTo[ipHELIUM][nelem][ion]* 2.86*(double)ion *
				dense.xIonDense[ipHELIUM][0]*
				dense.xIonDense[nelem][ion+1];
		}
	}
#endif

	/* convert from eV to ergs, HCharHeatOn usually 1, set to 0 with no CTHeat,  
	 * EN1EV is ergs in 1 eV, 1.602176e-012*/
	SumCTHeat_v *= EN1EV * atmdat.HCharHeatOn;

	if( thermal.htot > 1e-35 )
	{
		/* remember largest fractions of heating and cooling for comment */
		atmdat.HCharHeatMax = MAX2(atmdat.HCharHeatMax,
			SumCTHeat_v/thermal.htot );

		atmdat.HCharCoolMax = MAX2(atmdat.HCharCoolMax,
			-SumCTHeat_v/thermal.htot);
	}

	/* debug code to print out the contributors to total CT heating */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC)
		{
#			define FRAC	0.1
			for( nelem=ipHELIUM; nelem < LIMELM; nelem++ )
			{
				/* >>chng >>01 apr 25, loops had been to LIMELM, which may have done no harm
				 * since extra array elements were set to zero, but is incorrect since the physical
				 * limit is the number of stages of ionization */
				int limit = MIN2(4, nelem);
				/* this first group of lower stages of ionization have exact rate coefficients */
				for( ion=0; ion < limit; ion++ )
				{
					if(
						/* heating due to ionization of heavy element, recombination of hydrogen */
						(atmdat.CharExcIonOf[ipHYDROGEN][nelem][ion]*CTIonData[nelem][ion][7]*
						(double)dense.xIonDense[ipHYDROGEN][1]*(double)dense.xIonDense[nelem][ion] + 

						/* heating due to recombination of heavy element, ionization of hydrogen */
						atmdat.CharExcRecTo[ipHYDROGEN][nelem][ion]*CTRecombData[nelem][ion][6]*
						iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop()*
						(double)dense.xIonDense[nelem][ion+1])/SumCTHeat_v> FRAC )

						fprintf(ioQQQ,"DEBUG ct %li %li %.3f\n", nelem, ion, 
							(atmdat.CharExcIonOf[ipHYDROGEN][nelem][ion]*CTIonData[nelem][ion][7]*
							(double)dense.xIonDense[ipHYDROGEN][1]*(double)dense.xIonDense[nelem][ion] + 

							/* heating due to recombination of heavy element, ionization of hydrogen */
							atmdat.CharExcRecTo[ipHYDROGEN][nelem][ion]*CTRecombData[nelem][ion][6]*
							iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop()*
							(double)dense.xIonDense[nelem][ion+1])  );
				}

				for( ion=4; ion < nelem; ion++ )
				{
					if(
						(atmdat.CharExcRecTo[ipHYDROGEN][nelem][ion]* 2.86*(double)ion *
						(double)dense.xIonDense[ipHYDROGEN][0]*(double)dense.xIonDense[nelem][ion+1])/SumCTHeat_v> FRAC )
						fprintf(ioQQQ,"DEBUG ct %li %li %.3f\n", nelem, ion, 
						(atmdat.CharExcRecTo[ipHYDROGEN][nelem][ion]* 2.86*(double)ion *
						(double)dense.xIonDense[ipHYDROGEN][0]*(double)dense.xIonDense[nelem][ion+1]) );
				}
			}
#			undef FRAC
			fprintf(ioQQQ,"DEBUG ct tot %.3e\n", SumCTHeat_v / thermal.htot );
		}
	}
	return( SumCTHeat_v );
}

/*================================================================================*
 *================================================================================*/
STATIC double HCTIon(
	/* ion is stage of ionization on C scale, 0 for atom */
	long int ion,		
	/* nelem is atomic number of element on C scale, 1 to 29 */
	/* HCTIon(1,5) is C+ + H+ => C++ + H */
	long int nelem)	
{
	double HCTIon_v, 
	  tused;

	DEBUG_ENTRY( "HCTIon()" );
	/* H charge transfer ionization, using Jim Kingdon's ctdata.for */

	/* set up the rate coefficients if this is first call */
	if( !lgCTDataDefined )
	{
		if( trace.lgTrace )
		{
			fprintf( ioQQQ,"       HCTIon doing 1-time init of charge transfer data\n");
		}
		lgCTDataDefined = true;
		MakeHCTData();
	}

	/* check that data have been linked in,
	 * error here would mean that data below had not been loaded */
	ASSERT( CTIonData[2][0][0] > 0. );

	/* return zero for highly ionized species */
	if( ion >= 3 )
	{
		HCTIon_v = 0.;
		return( HCTIon_v );
	}

	/*begin sanity checks */
	/* check that ionization stage is ok for this element*/
	ASSERT( ion >= 0);
	ASSERT( ion <= nelem );

	/* now check the element is valid, this is routine HCTIon */
	ASSERT( nelem > 0 );
	ASSERT( nelem < LIMELM );

	/* may be no entry, first coefficient is zero in this case */
	if( CTIonData[nelem][ion][0] <= 0. )
	{
		HCTIon_v = 0.;

	}

	else
	{
		/* Make sure te is between temp. boundaries; set constant outside of range */
		tused = MAX2((double)phycon.te,CTIonData[nelem][ion][4]);
		tused = MIN2(tused,CTIonData[nelem][ion][5]);
		tused *= 1e-4;

		// make sure that the Boltzmann factor always uses the real temperature, even
		// outside the range of validity, so that the CT rate properly goes to zero
		// when the the gas temperature goes to zero.
		HCTIon_v = CTIonData[nelem][ion][0]*1e-9*(pow(tused,CTIonData[nelem][ion][1]))*
		  (1. + CTIonData[nelem][ion][2]*exp(CTIonData[nelem][ion][3]*tused))*
		  exp(-CTIonData[nelem][ion][6]*1.e4/phycon.te);
	}
	return( HCTIon_v );
}

/*================================================================================*
 *================================================================================*/
STATIC double HCTRecom(
	/* ion is stage of ionization on C scale, 0 for rec to atom */
	long int ion,  	
	/* nelem is atomic number of element on C scale, 1 = he up to LIMELM */
	/* HCTRecom(1,5) would be C+2 + H => C+ + H+ */
	long int nelem)	
{
	double HCTRecom_v, 
	  tused;

	DEBUG_ENTRY( "HCTRecom()" );
	/* 
	 * H charge transfer recombination using Jim Kingdon's block ctdata.for
	 */

	/* set up the rate coefficients if this is first call */
	if( !lgCTDataDefined )
	{
		if( trace.lgTrace )
		{
			fprintf( ioQQQ,"       HCTIon doing 1-time init of charge transfer data\n");
		}
		lgCTDataDefined = true;
		MakeHCTData();
	}

	/* this is check that data have been set up properly, will
	 * fail if arrays are not initialized properly */
	ASSERT( CTRecombData[1][0][0] > 0. );

	/* use Dalgarno estimate for highly ionized species, number reset with
	 * set charge transfer command */
	if( ion > 3 )
	{
		/* >>chng 96 nov 25, added this option, default is 1.92e-9
		 * Dalgarno's charge transfer */
		HCTRecom_v = atmdat.HCTAlex*((double)ion+1.);
		return( HCTRecom_v );
	}

	/* check that ion stage within bound for this atom */
	ASSERT( ion >= 0 && ion <= nelem );

	/* now check the element is valid, this is routine HCTIon */
	ASSERT( nelem > 0 && nelem < LIMELM );

	tused = MAX2((double)phycon.te,CTRecombData[nelem][ion][4]);
	tused = MIN2(tused,CTRecombData[nelem][ion][5]);
	tused *= 1e-4;

	if( tused == 0. )
	{
		HCTRecom_v = 0.;
		return( HCTRecom_v );
	}

	/* the interpolation equation */
	HCTRecom_v = CTRecombData[nelem][ion][0]*1e-9*(pow(tused,CTRecombData[nelem][ion][1]))*
	  (1. + CTRecombData[nelem][ion][2]*sexp(-CTRecombData[nelem][ion][3]*tused));

	/* in sexp negative sign not typo - there are negative signs already
	 * in coefficient, and sexp has implicit negative sign */
	return( HCTRecom_v );
}

/*================================================================================*
 *================================================================================*/
/*block data with Jim Kingdon's charge transfer data */
/* >>refer	H	CT	Kingdon, J. B., & Ferland, G. J. 1996, ApJS, 106, 205 */
/* 
 * first dimension is atomic number of atom, 0 for H 
 * second dimension is ionization stage,
 * 1=+0 for parent, etc
 * third dimension is atomic number of atom 
 * second dimension is ionization stage,
 * 1=+1 for parent, etc
 */

/* digital form of the fits to the charge transfer
 * ionization rate coefficients 
 *
 * Note: First parameter is in units of 1e-9!
 * Note: Seventh parameter is in units of 1e4 K */

/* digital form of the fits to the charge transfer
 * recombination rate coefficients (total)
 *
 * Note: First parameter is in units of 1e-9!
 * recombination 
 */

/* holds data for charge transfer fits */
STATIC void MakeHCTData()
{
	DEBUG_ENTRY( "MakeHCTData()" );

	/* 
	 * following are coefficients for charge transfer ionization,
	 * H+ + A => H + A+
	 */
	DataParser d("ctiondata.dat", ES_NONE);
	d.getline();
	d.checkMagic(201903041);
	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		for( long i=0; i < 4; ++i )
		{
			d.getline();
			d.getToken(&CTIonData[nelem][i][0], 8);
			d.checkEOL();
		}
	}
	d.checkEOD();

	/* CT recombination, A+n + H => A+n-1 + H+ */
	/* >>chng 01 may 03, first coefficient multiplied by 0.25, as per comment in
	 * >>refer	Li	CT	Stancil, P. C., & Zygelman, B. 1996, ApJ, 472, 102
	 * which corrected the error in 
	 * >>refer	He	CT	Zygelman, B., Dalgarno, A., Kimura, M., & Lane, N. F. 1989, Phys. Rev. A, 40, 2340
	 * this was used in the original Kingdon & Ferland paper so no correction required
	 * >>chng 04 apr 27, He was in error above as well, factor of 4, noted in 
	 * >>refer	He	CT	Stancil, P. C., Lepp, S., & Dalgarno, A. 1998, ApJ, 509, 1
	 */
	d.open( "ctrecombdata.dat", ES_NONE );
	d.getline();
	d.checkMagic(201903042);
	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
	{
		for( long i=0; i < 4; ++i )
		{
			d.getline();
			d.getToken(&CTRecombData[nelem][i][0], 7);
			d.checkEOL();
		}
	}
	d.checkEOD();
}

/* ChargTranPun, save charge transfer coefficient */
void ChargTranPun( FILE* ipPnunit , char* chSave )
{
	long int j, jj;
	/* restore temp when done with this routine	*/
	double TempSave = phycon.te;

	DEBUG_ENTRY( "ChargTranPun()" );

	/* this is the usual charge transfer option */
	if( strcmp( chSave,"CHAR") == 0 )
	{
		/* charge exchange rate coefficients, entered with
		 * save charge transfer command.  Queries Jim Kingdon's
		 * CT fits and routines to get H - He and higher CT rates
		 * rates are evaluated at the current temperature, which can
		 * be specified with the constant temperature command */
		/* first group is charge transfer recombination,
		 * the process A+x + H => A+x-1 + H+ */
		fprintf( ipPnunit, "#element\tion\n");
		for( j=1; j < LIMELM; j++ )
		{
			/* first number is atomic number, so add 1 to get onto physical scale */
			/* >>chng 00 may 09, caught by Jon Slavin */
			fprintf( ipPnunit, "%s\t", elementnames.chElementSym[j] );
			/*fprintf( ipPnunit, "%3ld\t", j+1 );*/

			for( jj=0; jj < j; jj++ )
			{
				fprintf( ipPnunit, "%.2e\t", 
				  HCTRecom(jj,j) );
			}
			fprintf( ipPnunit, "\n" );
		}

		/* second group is charge transfer ionization,
		 * the process A+x + H+ => A+x+1 + H */
		fprintf( ipPnunit, "\n#ionization rates, atomic number\n");
		for( j=1; j < LIMELM; j++ )
		{
			fprintf( ipPnunit, "%s\t", elementnames.chElementSym[j] );
			for( jj=0; jj < j; jj++ )
			{
				fprintf( ipPnunit, "%.2e\t", 
				  HCTIon(jj,j) );
			}
			fprintf( ipPnunit, "\n" );
		}
	}

	/* this is the charge transfer to produce output for AGN3,
	 * invoked with the save charge transfer AGN command */
	else if( strcmp( chSave,"CHAG") == 0 )
	{
		/* this is boundary between two tables */
		double BreakEnergy = 100./13.0;
		realnum teinit = 5000.;
		realnum tefinal = 20000.;
		realnum te1;
		long int nelem, ion;

		te1 = teinit;
		fprintf(ipPnunit,"H ioniz\n X+i\\Te");
		while( te1 <= tefinal )
		{
			fprintf(ipPnunit,"\t%.0f K",te1);
			te1 *= 2.;
		}
		fprintf(ipPnunit,"\n");

		/* make sure rates already evaluated at least one time */
		ChargTranEval();

		/* loop over all elements, H charge transfer ionization */
		for( nelem=ipHELIUM; nelem<LIMELM; ++nelem )
		{
			/* this list of elements included in the AGN tables is defined in zeroabun.c */
			if( abund.lgAGN[nelem] )
			{
				ion_trim_untrim(nelem);

				for( ion=0; ion<=nelem; ++ion )
				{
					/* skip high ionization CT */
					if( Heavy.Valence_IP_Ryd[nelem][ion] > BreakEnergy )
						break;
					/* most of these are actually zero */
					if( atmdat.CharExcIonOf[ipHYDROGEN][nelem][ion] == 0 )
						continue;

					/* print chemical symbol */
					fprintf(ipPnunit,"%s", 
						elementnames.chElementSym[nelem]);
					/* now ionization stage */
					if( ion==0 )
					{
						fprintf(ipPnunit,"0 ");
					}
					else if( ion==1 )
					{
						fprintf(ipPnunit,"+ ");
					}
					else
					{
						fprintf(ipPnunit,"+%li",ion);
					}

					/* fully define the new temperature */
					TempChange(teinit , false);

					while( phycon.te <= tefinal )
					{
						ASSERT (dense.IonLow[nelem] == 0 &&	dense.IonHigh[nelem] == nelem+1);
						ChargTranEval();

						fprintf(ipPnunit,"\t%.2e",atmdat.CharExcIonOf[ipHYDROGEN][nelem][ion]);
						TempChange(phycon.te *2.f , false);
					}
					fprintf(ipPnunit,"\n");
				}
				fprintf(ipPnunit,"\n");
			}
		}

		te1 = teinit;
		fprintf(ipPnunit,"H recom\n X+i\\Te");
		while( te1 <= tefinal )
		{
			fprintf(ipPnunit,"\t%.0f K",te1);
			te1 *= 2.;
		}
		fprintf(ipPnunit,"\n");

		/* loop over all elements, H charge transfer recombination */
		for( nelem=ipHELIUM; nelem<LIMELM; ++nelem )
		{
			/* this list of elements included in the AGN tables is defined in zeroabun.c */
			if( abund.lgAGN[nelem] )
			{
				ion_trim_untrim(nelem);
				for( ion=0; ion<=nelem; ++ion )
				{
					/* skip high ionization CT */
					if( Heavy.Valence_IP_Ryd[nelem][ion] > BreakEnergy )
						break;
					/* most of these are actually zero */
					if( atmdat.CharExcRecTo[ipHYDROGEN][nelem][ion] == 0 )
						continue;

					/* print chemical symbol */
					fprintf(ipPnunit,"%s", 
						elementnames.chElementSym[nelem]);
					/* now ionization stage */
					if( ion==0 )
					{
						fprintf(ipPnunit,"0 ");
					}
					else if( ion==1 )
					{
						fprintf(ipPnunit,"+ ");
					}
					else
					{
						fprintf(ipPnunit,"+%li",ion);
					}

					/* fully define the new temperature */
					TempChange(teinit , false);
					while( phycon.te <= tefinal )
					{
						ASSERT (dense.IonLow[nelem] == 0 &&	dense.IonHigh[nelem] == nelem+1);
						ChargTranEval();

						fprintf(ipPnunit,"\t%.2e",atmdat.CharExcRecTo[ipHYDROGEN][nelem][ion]);
						TempChange(phycon.te *2.f , false);
					}
					fprintf(ipPnunit,"\n");
				}
				fprintf(ipPnunit,"\n");
			}
		}

#		if 0
		te1 = teinit;
		fprintf(ipPnunit,"He recom\n Elem\\Te");
		while( te1 <= tefinal )
		{
			fprintf(ipPnunit,"\t%.0f",te1);
			te1 *= 2.;
		}
		fprintf(ipPnunit,"\n");

		/* loop over all elements, H charge transfer recombination */
		for( nelem=ipHELIUM; nelem<LIMELM; ++nelem )
		{
			/* this list of elements included in the AGN tables is defined in zeroabun.c */
			if( abund.lgAGN[nelem] )
			{
				ion_trim_untrim(nelem);
				for( ion=0; ion<=nelem; ++ion )
				{
					/* most of these are actually zero */
					if( atmdat.CharExcRecTo[ipHELIUM][nelem][ion] == 0 )
						continue;
					fprintf(ipPnunit,"%.2s%.2s", 
						elementnames.chElementSym[nelem],
						elementnames.chIonStage[ion]);

					/* fully define the new temperature */
					TempChange(teinit , false);
					while( phycon.te <= tefinal )
					{
						ASSERT (dense.IonLow[nelem] == 0 &&	dense.IonHigh[nelem] == nelem+1);
						ChargTranEval();

						fprintf(ipPnunit,"\t%.2e",atmdat.CharExcRecTo[ipHELIUM][nelem][ion]);
						TempChange(phycon.te *2.f , false);
					}
					fprintf(ipPnunit,"\n");
				}
				fprintf(ipPnunit,"\n");
			}
		}


		te1 = teinit;
		fprintf(ipPnunit,"He ioniz\n Elem\\Te");
		while( te1 <= tefinal )
		{
			fprintf(ipPnunit,"\t%.0f",te1);
			te1 *= 2.;
		}
		fprintf(ipPnunit,"\n");

		/* loop over all elements, H charge transfer recombination */
		for( nelem=ipHELIUM; nelem<LIMELM; ++nelem )
		{
			/* this list of elements included in the AGN tables is defined in zeroabun.c */
			if( abund.lgAGN[nelem] )
			{
				ion_trim_untrim(nelem);
				for( ion=0; ion<=nelem; ++ion )
				{
					/* most of these are actually zero */
					if( atmdat.CharExcIonOf[ipHELIUM][nelem][ion] == 0 )
						continue;
					fprintf(ipPnunit,"%.2s%.2s", 
						elementnames.chElementSym[nelem],
						elementnames.chIonStage[ion]);

					/* fully define the new temperature */
					TempChange(teinit , false);
					while( phycon.te <= tefinal )
					{
						ASSERT (dense.IonLow[nelem] == 0 &&	dense.IonHigh[nelem] == nelem+1);
						ChargTranEval();

						fprintf(ipPnunit,"\t%.2e",atmdat.CharExcIonOf[ipHELIUM][nelem][ion]);
						TempChange(phycon.te*2.f , true);
					}
					fprintf(ipPnunit,"\n");
				}
				fprintf(ipPnunit,"\n");
			}
		}
#		endif
	}
	else
	{
		fprintf( ioQQQ, " save charge keyword insane\n" );
		cdEXIT(EXIT_FAILURE);
	}

	TempChange(TempSave , false);
	return;
}
