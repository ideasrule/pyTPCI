/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* HyperfineCreat establish space for hf arrays, reads atomic data from hyperfine.dat */
/* HyperfineCS - returns collision strengths for hyperfine struc transitions */
/*H21cm computes rate for H 21 cm from upper to lower excitation by atomic hydrogen */ 
/*h21_t_ge_20 compute rate for H 21 cm from upper to lower excitation by atomic hydrogen */ 
/*h21_t_lt_20 compute rate for H 21 cm from upper to lower excitation by atomic hydrogen */ 
/*H21cm_electron compute H 21 cm rate from upper to lower excitation by electrons - call by CoolEvaluate */
/*H21cm_H_atom - evaluate H atom spin changing collision rate, called by CoolEvaluate */
/*H21cm_proton - evaluate proton spin changing H atom collision rate, */
#include "cddefines.h"
#include "abund.h"
#include "conv.h"
#include "phycon.h"
#include "dense.h"
#include "rfield.h"
#include "taulines.h"
#include "iso.h"
#include "trace.h"
#include "hyperfine.h"
#include "lines_service.h"
#include "service.h"

/* H21_cm_pops - fine level populations for 21 cm with Lya pumping included 
 * called in CoolEvaluate */
void H21_cm_pops( void )
{
	/*atom_level2( HFLines[0] );*/
	/*return;*/
	/*
	things we know on entry to this routine:
	total population of 2p: iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop
	total population of 1s: iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop
	continuum pumping rate (lo-up) inside 21 cm line: HFLines[0].pump()
	upper to lower collision rate inside 21 cm line: HFLines[0].cs*dense.cdsqte
	occupation number inside Lya: OccupationNumberLine( &iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s) )

	level populations (cm-3) must be computed:
	population of upper level of 21cm: HFLines[0].Hi->Pop
	population of lower level of 21cm: (*HFLines[0].Lo()).Pop
	stimulated emission corrected population of lower level: HFLines[0].Emis->PopOpc()
	*/

	double	PopTot = iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop();

	/* population can be zero in certain tests where H is turned off,
	 * also if initial solver does not see any obvious source of ionization
	 * also possible to set H0 density to zero with element ionization command,
	 * as is done in func_set_ion test case */
	if( PopTot <0 )
		TotalInsanity();
	else if( PopTot == 0 )
	{
		/*return after zeroing local variables */
		(*HFLines[0].Hi()).Pop() = 0.;
		(*HFLines[0].Lo()).Pop() = 0.;
		HFLines[0].Emis().PopOpc() = 0.;
		HFLines[0].Emis().xIntensity() = 0.;
		HFLines[0].Emis().xObsIntensity() = 0.;
		HFLines[0].Emis().ColOvTot() = 0.;
		hyperfine.Tspin21cm = 0.;
		return;
	}

	double e1 = 0.;
	double e2 = HFLines[ 0 ].EnergyWN();
	/* The 2p fine stucture energies are current with NIST as of May 14, 2019. */
	double e2p12 = 82258.9191133;
	double e2p32 = 82259.2850014;
	/* The hyperfine splittings of the 2p fine structure levels are from
	 * >>refer	HI	Bethe & Salpeter (1977) Section 22, page 110. */
	double e2p12_splitting = e2 / 24.;
	double e2p32_splitting = e2 / 60.;

	/* The hyperfine states have statistical weights 2F+1, so they differ from
	 * the unsplit level energy by:
	 *   El = E - dE * gu / (gu+gl)
	 *   Eu = E + dE * gl / (gu+gl)
	 * where E is the unsplit energy, dE the hyperfine splitting, and gu, gl the
	 * statistical weights of the hyperfine states.
	 * For 2p1/2: g(F=1) = 3, g(F=0) = 1
	 * For 2p3/2: g(F=2) = 5, g(F=1) = 3
	 * The levels of interest here are the 2p1/2(F=1) and 2p3/2(F=1), the top
	 * and bottom levels of the hyperfine states, resp.
	 * >>refer	HI	Deguchi & Watson 1985 ApJ, 290, 578
	 *	refcon		see their Fig. 1
	 */
	double e3 = e2p12 + 0.25 * e2p12_splitting;
	double e4 = e2p32 - 0.625 * e2p32_splitting;

	double de31 = e3 - e1;
	double de32 = e3 - e2;
	double de41 = e4 - e1;
	double de42 = e4 - e2;

	if( false )
	{
		fprintf( ioQQQ, "-------\n" );
		fprintf( ioQQQ, "de32 = %.9e\n", de32 );
		fprintf( ioQQQ, "de31 = %.9e\n", de31 );
		fprintf( ioQQQ, "de42 = %.9e\n", de42 );
		fprintf( ioQQQ, "de41 = %.9e\n", de41 );
		fprintf( ioQQQ, "-------\n" );
	}

	double a31 = 2.08e8;   /* Einstein co-efficient for transition 1p1/2 to 0s1/2 */
	double a32 = 4.16e8;   /* Einstein co-efficient for transition 1p1/2 to 1s1/2 */
	double a41 = 4.16e8;   /* Einstein co-efficient for transition 1p3/2 to 0s1/2 */
	double a42 = 2.08e8;   /* Einstein co-efficient for transition 1p3/2 to 1s1/2 */
	/* These A values are determined from eqn. 17.64 of "The theory of Atomic structure
	 * and Spectra" by R. D. Cowan 
	 * A hyperfine level has degeneracy Gf=(2F + 1)
	 * a2p1s = 6.24e8;  Einstein co-efficient for transition 2p to 1s */
	double a21 = HFLines[0].Emis().Aul();	/* Einstein co-efficient for transition 1s1/2 to 0s1/2 */

	/* above is spontaneous rate - the net rate is this times escape and destruction
	 * probabilities */
	a21 *= HFLines[0].Emis().Ploss();
	ASSERT( a21 >= 0. );

	/* hyperfine.lgLya_pump_21cm is option to turn off Lya pump
	 * of 21 cm, with 'no 21cm lya pump' command */
	double occnu_lya = OccupationNumberLine( iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s) ) *
		hyperfine.lgLya_pump_21cm;

	if( !conv.lgSearch && occnu_lya < 0. )
	{
		occnu_lya = 0.;
		fixit( "PopOpc <0 but Pesc > 0: We may need to review when Pesc is computed to get non-negative occupation numbers" );
	}

	/* Lya occupation number for the hyperfine levels 0S1/2 and 1S1/2 can be different
	 * this is related to the "Wouthuysen-Field coupling",
	 * https://en.wikipedia.org/wiki/Wouthuysen%E2%80%93Field_coupling
	 * which assumes that the variation of the Lya source function is the gas kinetic temperature.
	 * Following Adams 1971 we assume variation is line excitation temperature.
	 * Third possibility is that given in stellar atmosphere texts, S_nu = constant
	 */
	double occnu_lya_23 = occnu_lya,
	       occnu_lya_13 = 0.,
	       occnu_lya_24 = 0.,
	       occnu_lya_14 = 0.;

	/* selected with SET LYA 21CM command */
	if( hyperfine.LyaSourceFunctionShape_assumed == t_hyperfine::EXCITATION ||
	    hyperfine.LyaSourceFunctionShape_assumed == t_hyperfine::KINETIC )
	{
		double Temp = phycon.te;

		if( hyperfine.LyaSourceFunctionShape_assumed == t_hyperfine::EXCITATION )
			Temp = TexcLine( iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s) );

		if( occnu_lya * Temp > 0. )
		{
			/* If the continuum is described by a Planck function, then the continuum
			 * within Lya seen by the two levels is not exactly of the same brightness.
			 * They differ by the exp when Lya is on the Wien tail of black body, which
			 * must be true if 21 cm is important. */

			double texc1 = sexp( HFLines[0].EnergyK() / Temp );
			double texc2 = sexp( ((e2p32-e2p12)*T1CM) / Temp );

			occnu_lya_23 = occnu_lya;
			occnu_lya_13 = occnu_lya * texc1;
			occnu_lya_24 = occnu_lya * texc2;
			occnu_lya_14 = occnu_lya_13 * texc2;
		}

		enum { DEBUG_SPEC = false };
		if( DEBUG_SPEC )
		{
			fprintf(ioQQQ,"DEBUG texc %12.3e excitation %12.3e kinetic %12.3e\n",
				Temp, TexcLine( iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s) ) , phycon.te );
		}
	}
	else if( hyperfine.LyaSourceFunctionShape_assumed == t_hyperfine::CONSTANT )
	{
		occnu_lya_23 = occnu_lya;
		occnu_lya_13 = powi(de32/de31, 3) * occnu_lya;
		occnu_lya_24 = powi(de32/de42, 3) * occnu_lya;
		occnu_lya_14 = powi(de32/de41, 3) * occnu_lya;
	}
	else
		TotalInsanity();

	if( false )
	{
		fprintf( ioQQQ, "=======\n" );
		fprintf( ioQQQ, "oc32 = %.9e\n", occnu_lya_23 );
		fprintf( ioQQQ, "oc31 = %.9e\n", occnu_lya_13 );
		fprintf( ioQQQ, "oc42 = %.9e\n", occnu_lya_24 );
		fprintf( ioQQQ, "oc41 = %.9e\n", occnu_lya_14 );
		fprintf( ioQQQ, "=======\n" );
	}

	/* this is the 21 cm upward continuum pumping rate [s-1] for the attenuated incident and
	 * local continuum and including line optical depths */
	double pump12 = HFLines[0].Emis().pump();
	double pump21 = pump12 * (*HFLines[0].Lo()).g() / (*HFLines[0].Hi()).g();

	/* collision rates s-1 within 1s,
	 * were multiplied by collider density when evaluated in CoolEvaluate */
	/* ContBoltz is Boltzmann factor for wavelength of line */
	ASSERT( HFLines[0].Coll().col_str()>0. );
	double coll12 = HFLines[0].Coll().col_str()*dense.cdsqte/(*HFLines[0].Lo()).g()*rfield.ContBoltz[HFLines[0].ipCont()-1];
	double coll21 = HFLines[0].Coll().col_str()*dense.cdsqte/(*HFLines[0].Hi()).g();

	/* set up rate (s-1) equations
	 * all process out of 1 that eventually go to 2 */
	double rate12 = 
		/* collision rate (s-1) from 1 to 2 */
		coll12 + 
		/* direct external continuum pumping (s-1) in 21 cm line - usually dominated by CMB */
		pump12 +
		/* pump rate (s-1) up to 3, times fraction that decay to 2, hence net 1-2 */
		3.*a31*occnu_lya_13 *a32/(a31+a32)+
		/* pump rate (s-1) up to 4, times fraction that decay to 2, hence net 1-2 */
		/* >>chng 05 apr 04, GS, degeneracy corrected from 6 to 3 */
		3.*a41*occnu_lya_14 *a42/(a41+a42);

	/* set up rate (s-1) equations
	 * all process out of 2 that eventually go to 1 */
	/* spontaneous + induced 2 -> 1 by external continuum inside 21 cm line */
	/* >>chng 04 dec 03, do not include spontaneous decay, for numerical stability */
	double rate21 = 
		/* collisional deexcitation */
		coll21 +
		/* net spontaneous decay plus external continuum pumping in 21 cm line */
		pump21 +
		/* rate from 2 to 3 time fraction that go back to 1, hence net 2 - 1 */
		/* >>chng 05 apr 04,GS, degeneracy corrected from 2 to unity */
		occnu_lya_23*a32 * a31/(a31+a32)+
		occnu_lya_24*a42*a41/(a41+a42);

	/* x = (*HFLines[0].Hi()).Pop/(*HFLines[0].Lo()).Pop */
	double x = rate12 / SDIV(a21 + rate21);
	ASSERT( x > 0. );

	/* the Transitions term is the total population of 1s */
	(*HFLines[0].Hi()).Pop() = (x/(1.+x))* PopTot;
	(*HFLines[0].Lo()).Pop() = (1./(1.+x))* PopTot;

	/* the population with correction for stimulated emission */
	HFLines[0].Emis().PopOpc() = (*HFLines[0].Lo()).Pop()*((3*rate21- rate12) + 3*a21)/SDIV(3*(a21+ rate21));

	/* ratio of collisional to total (collisional + pumped) excitation */
	HFLines[0].Emis().ColOvTot() = 0.;
	if( rate12 > 0. )
		HFLines[0].Emis().ColOvTot() = coll12 / rate12;

	/* set number of escaping line photons, used elsewhere for outward beam
	 * and line intensity
	 * NB: continuum subtraction is performed within PutLine() */
	set_xIntensity(HFLines[0]);

	/* finally save the spin temperature */
	hyperfine.Tspin21cm = HyperfineTspin( HFLines[0] );

	return;
}

/*H21cm_electron computes rate for H 21 cm from upper to lower excitation by electrons - call by CoolEvaluate
 * >>refer	H1	CS	Smith, F. J. 1966, Planet. Space Sci., 14, 929 */
double H21cm_electron( double temp )
{
	temp = MIN2(1e4 , temp );

	/* following fit is from */
	/* >>refer	H1	21cm	Liszt, H. 2001, A&A, 371, 698 */
	return	exp10( -9.607 + log10( sqrt(temp)) * sexp( powpq(log10(temp),9,2) / 1800. ));
}

/* computes rate for H 21 cm from upper to lower excitation by atomic hydrogen 
 * from 
 * >>refer	H1	CS	Allison, A. C., & Dalgarno A. 1969, ApJ 158, 423 */
/* the following is the best current survey of 21 cm excitation */
/* >>refer	H1	21cm	Liszt, H. 2001, A&A, 371, 698 */
#if 0
STATIC double h21_t_ge_20( double temp )
{
	double y;
	double x1,
		teorginal = temp;
	/* data go up to 1,000K must not go above this */
	temp = MIN2( 1000.,temp );
	x1 =1.0/sqrt(temp);
	y =-21.70880995483007-13.76259674006133*x1;
	y = exp(y);

	/* >>chng 02 feb 14, extrapolate above 1e3 K as per Liszt 2001 recommendation 
	 * page 699 of */
	/* >>refer	H1	21cm	Liszt, H. 2001, A&A, 371, 698 */
	if( teorginal > 1e3 )
	{
		y *= pow(teorginal/1e3 , 0.33 );
	}

	return( y );
}

/* this branch for T < 20K, data go down to 1 K */
STATIC double h21_t_lt_20( double temp )
{
	double y;
	double x1;

	/* must not go below 1K */
	temp = MAX2( 1., temp );
	x1 =temp*log(temp);
	y =9.720710314268267E-08+6.325515312006680E-08*x1;
	return(y*y);
}
#endif

/* >> chng 04 dec 15, GS. The fitted rate co-efficients (cm3s-1) in the temperature range 1K to 300K is from
 * >>refer	H1	CS	Zygelman, B. 2005, ApJ, 622, 1356 
 * The rate is 4/3 times the Dalgarno (1969) rate for the 
 temperature range 300K to 1000K. Above 1000K, the rate is extrapolated according to Liszt 2001.*/
STATIC double h21_t_ge_10( double temp )
{
	double teorginal = temp;

	/* data go up to 300K  */
	temp = MIN2( 300., temp );

	double y = 1.4341127e-9
		 + 9.4161077e-15 * temp
		 - 9.2998995e-9  / log(temp)
		 + 6.9539411e-9  / sqrt(temp)
		 + 1.7742293e-8  * log(temp)/pow2(temp);
	if( teorginal > 300. )
	{
		/* data go up to 1000*/
		temp = MIN2( 1000., teorginal );
		y = -21.70880995483007 - 13.76259674006133 / sqrt(temp);
		y = 1.236686*exp(y);
	}
	if( teorginal > 1e3 )
	{
		/*data go above 1000*/
		y *= pow( teorginal/1e3 , 0.33 );
	}
	return( y );
}
/* this branch for T < 10K, data go down to 1 K */
STATIC double h21_t_lt_10( double temp )
{
	/* must not go below 1K */
	temp = MAX2(1., temp );
	return	  8.5622857e-10
		+ 2.331358e-11  * temp
		+ 9.5640586e-11 * pow2(log(temp))
		- 4.6220869e-10 * sqrt(temp)
		- 4.1719545e-10 / sqrt(temp);
}

/*H21cm_H_atom - evaluate H atom spin changing H atom collision rate, 
 * called by CoolEvaluate 
 * >>refer	H1	CS	Allison, A. C. & Dalgarno, A. 1969, ApJ 158, 423 
 */
double H21cm_H_atom( double temp )
{
	double hold;
	if( temp >= 10. )
	{
		hold = h21_t_ge_10( temp );
	}
	else
	{
		hold = h21_t_lt_10( temp );
	}

	return hold;
}

/*H21cm_proton - evaluate proton spin changing H atom collision rate, 
* called by CoolEvaluate */
double H21cm_proton( double temp )
{
	/*>>refer	21cm	p coll	Furlanetto, S. R. & Furlanetto, M. R. 2007, MNRAS, 379, 130
	 * previously had used proton rate, which is 3.2 times H0 rate according to
	 *>>refer	21cm	CS	Liszt, H. 2001, A&A, 371, 698 */
	/* fit to table 1 of first paper */
	/*--------------------------------------------------------------*
	TableCurve Function: c:\storage\litera~1\21cm\atomic~1\p21cm.c Jun 20, 2007 3:37:50 PM
	proton coll deex
	X= temperature (K)
	Y= rate coefficient (1e-9 cm3 s-1)
	Eqn# 4419  y=a+bx+cx^2+dx^(0.5)+elnx/x
	r2=0.9999445384690351
	r2adj=0.9999168077035526
	StdErr=5.559328579039901E-12
	Fstat=49581.16793656295
	a= 9.588389834316704E-11
	b= -5.158891920816405E-14
	c= 5.895348443553458E-19
	d= 2.05304960232429E-11
	e= 9.122617940315725E-10
	*--------------------------------------------------------------*/

	/* only fit this range, did not include T = 1K point which 
	 * causes an inflection */
	temp = MAX2( 2. , temp );
	temp = MIN2( 2e4 , temp );

	/* within range of fitted rate coefficients */
	return	  9.588389834316704E-11
		- 5.158891920816405E-14 * temp
		+ 5.895348443553458E-19 * temp * temp
		+ 2.053049602324290E-11 * sqrt(temp)
		+ 9.122617940315725E-10 * log(temp) / temp;
}

/* 
 * HyperfineCreate, HyperfineCS written July 2001
 * William Goddard for Gary Ferland
 * This code calculates line intensities for known
 * hyperfine transitions.
 */

/* two products, the transition structure HFLines, which contains all information for the lines,
 * and nHFLines, the number of these lines.  
 *
 * these are in taulines.h
 *
 * info to create them contained in hyperfine.dat
 *
 * abundances of nuclei are also in hyperfine.dat, stored in 
 */

/* Ion contains varying temperatures, specified above, used for	*/
/* calculating collision strengths.				*/	
static int Ntemp = -1;
static vector<double> csTemp;

typedef struct 
{
	vector<double> cs;
	vector<double> cs2d;
} t_ColStr;

static vector<t_ColStr> colstr;

const double ENERGY_MIN_WN = 1e-10;


/* HyperfineCreate establish space for HFS arrays, reads atomic data from hyperfine.dat */
void HyperfineCreate(void)
{
	vector<string>	data;

	DEBUG_ENTRY( "HyperfineCreate()" );

	/* list of ion collision strengths for the temperatures listed in table */
	/* HFLines containing all the data in Hyperfine.dat, and transition is	*/
	/* defined in cddefines.h						*/

	/*transition *HFLines;*/

	/* get the line data for the hyperfine lines */
	if( trace.lgTrace )
		fprintf( ioQQQ," Hyperfine opening hyperfine.dat:");

	FILE *ioDATA = open_data( "hyperfine.dat", "r" );

	/* first line is a version number and does not count */
	string chLine;
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " Hyperfine could not read first line of hyperfine.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* count lines in the file, ignoring lines starting with '#',
	 * and get temperature array for HFS collision strengths */
	size_t nHFLines = 0;
	while( read_whole_line( chLine, ioDATA ) )
	{
		if( chLine[0] == '#' )
		{
			continue;
		}
		else if( chLine.find("TDATA:") == string::npos )
		{
			Split(chLine, "\t", data, SPM_STRICT);
			int Aiso  = atoi(data[0].c_str());
			int nelem = atoi(data[1].c_str());
			double  wavelength = atof(data[3].c_str());

			if( abund.IsoAbn[nelem-1].getAbn( Aiso ) > 0 &&
				WAVNRYD/wavelength > rfield.emm() )
				++nHFLines;
		}
		else
		{
			Split(chLine, " ", data, SPM_STRICT);

			if (data.size() <= 1)
			{
				fprintf(ioQQQ, "HyperfineCreate: Error: Invalid number of temperatures in 'TDATA:': %d\n",
						(int) data.size());
				cdEXIT(EXIT_FAILURE);
			}

			Ntemp = data.size() - 1;
			csTemp.resize(Ntemp);

			int i = 0;
			for (std::vector<string>::const_iterator it = data.begin(); it != data.end() && i <= Ntemp; it++, i++)
			{
				if(i == 0)
					continue;
				csTemp[i-1] = atof((*it).c_str());
			}
		}

		data.resize(0);
	}

	ASSERT(nHFLines > 0 && Ntemp > 0);
	for(int i = 0; i < Ntemp; i++)
	{
		ASSERT( csTemp[i] > phycon.TEMP_LIMIT_LOW &&
			csTemp[i] < phycon.TEMP_LIMIT_HIGH );
		if( i > 0 )
			ASSERT(csTemp[i] > csTemp[i-1]);
		//	printf("i=%d\t t = %g\n", i, csTemp[i]);
	}

	/* allocate the transition HFLines array */
	HFLines.resize(nHFLines);
	AllTransitions.push_back(HFLines);

	/* initialize array to impossible values to make sure eventually done right */
	for( size_t i=0; i< HFLines.size(); ++i )
	{
		HFLines[i].Junk();
		HFLines[i].AddHiState();
		HFLines[i].AddLoState();
		HFLines[i].AddLine2Stack();
	}

	colstr.resize(HFLines.size());
	for (size_t j = 0; j < HFLines.size(); j++)
	{
		colstr[j].cs.resize(Ntemp);
		colstr[j].cs2d.resize(Ntemp);
	}

	hyperfine.HFLabundance.resize(HFLines.size());

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " Hyperfine could not rewind hyperfine.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* check that magic number is ok, read the line */
	if( !read_whole_line( chLine, ioDATA ) )
	{
		fprintf( ioQQQ, " Hyperfine could not read first line of hyperfine.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* check that magic number is ok, scan it in */
	{
		long j = 1;
		bool lgEOL;
		int year  = (int) FFmtRead(chLine.c_str(),&j,chLine.length(),&lgEOL);
		int month = (int) FFmtRead(chLine.c_str(),&j,chLine.length(),&lgEOL);
		int day   = (int) FFmtRead(chLine.c_str(),&j,chLine.length(),&lgEOL);

		/* the following is the set of numbers that appear at the start of hyperfine.dat 13 02 09 */
		const int iYR=13, iMN=10, iDY=18;
		if( ( year != iYR ) || ( month != iMN ) || ( day != iDY ) )
		{
			fprintf( ioQQQ, 
				" Hyperfine: the version of hyperfine.dat in the data directory is not the current version.\n" );
			fprintf( ioQQQ, 
				" I expected to find the number %i %i %i and got %i %i %i instead.\n" ,
				iYR, iMN , iDY ,
				year , month , day );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* 
	 * scan the string taken from Hyperfine.dat, parsing into
	 * needed variables.
	 * nelem is the atomic number.
	 * IonStg is the ionization stage.  Atom = 1, Z+ = 2, Z++ = 3, etc.
	 * Aul is used to find the oscillator strength in the function GetGF.
	 * most of the variables are floats.
	 */

	size_t j = 0;
	while( j < HFLines.size() && read_whole_line( chLine, ioDATA ) )
	{
		/* skip lines starting with '#' or containing the temperature array */
		if( chLine[0] == '#' || chLine.find("TDATA:") != string::npos )
			continue;

		Split(chLine, "\t", data, SPM_STRICT);
		int Aiso  = atoi(data[0].c_str());
		int nelem = atoi(data[1].c_str());
		double wavelength = atof(data[3].c_str());

		/* Ignore lines that fall beyond the lowest energy. */
		if( ! ( abund.IsoAbn[nelem-1].getAbn( Aiso ) > 0 &&
			WAVNRYD/wavelength > rfield.emm() ) )
		{
			data.resize(0);
			continue;
		}

		(*HFLines[j].Hi()).nelem() = nelem;
		ASSERT((*HFLines[j].Hi()).nelem()  > 0);

		(*HFLines[j].Hi()).IonStg() = atoi(data[2].c_str());
		ASSERT((*HFLines[j].Hi()).IonStg() > 0);

		hyperfine.HFLabundance[j] = abund.IsoAbn[nelem-1].getAbn( Aiso );
		ASSERT(hyperfine.HFLabundance[j] >= 0.0 && hyperfine.HFLabundance[j] <= 1.0);

		HFLines[j].Emis().Aul()	= (realnum) atof(data[4].c_str());

		(*HFLines[j].Hi()).g() = (realnum) (2*(abund.IsoAbn[nelem-1].getSpin( Aiso ) + .5) + 1);
		(*HFLines[j].Lo()).g() = (realnum) (2*(abund.IsoAbn[nelem-1].getSpin( Aiso ) - .5) + 1);

		/* account for inverted levels */
		if( abund.IsoAbn[nelem-1].getMagMom( Aiso ) < 0 )
		{
			double	tmp = (*HFLines[j].Hi()).g();
			(*HFLines[j].Hi()).g() = (*HFLines[j].Lo()).g();
			(*HFLines[j].Lo()).g() = tmp;
		}

		double fenergyWN = MAX2(ENERGY_MIN_WN, 1.0 / wavelength);
		HFLines[j].WLAng() = (realnum)(wavelength * 1e8f);
		HFLines[j].EnergyWN() = (realnum) fenergyWN;

		HFLines[j].Emis().dampXvel() = (realnum)( HFLines[j].Emis().Aul()
							/ HFLines[j].EnergyWN() / PI4 );
		HFLines[j].Emis().damp() = 1e-20f;

		HFLines[j].Emis().gf() = (realnum)(GetGF(HFLines[j].Emis().Aul(), fenergyWN, (*HFLines[j].Hi()).g()));
		ASSERT(HFLines[j].Emis().gf() > 0.0);

		(*HFLines[j].Lo()).nelem() = (*HFLines[j].Hi()).nelem();
		(*HFLines[j].Lo()).IonStg() = (*HFLines[j].Hi()).IonStg();

		//	printf("line %3ld:\t A= %2d\t Z= %2d\t Spin= %3.1f\t MagMom= %8.5f\t IonStg= %2d\t Frac= %6.4f\t"
		//		" Wl= %7.4f\t Aul= %.4e\t glo= %1.0f\t ghi= %1.0f\n",
		//		j, Aiso, nelem,
		//		abund.IsoAbn[nelem-1].getSpin( Aiso ),
		//		abund.IsoAbn[nelem-1].getMagMom( Aiso ),
		//		(*HFLines[j].Hi()).IonStg(),
		//		hyperfine.HFLabundance[j], wavelength, HFLines[j].Emis().Aul(),
		//		(*HFLines[j].Lo()).g(), (*HFLines[j].Hi()).g());

 
		if( data.size() > 6 )
		{
			//	printf("data for line %ld\t %d\t %d:\t", j, nelem, (*HFLines[j].Hi()).IonStg());
			for (int ij = 6, ii = 0; ij < (int) data.size() && ii < Ntemp; ij++, ii++)
			{
				colstr[j].cs[ii] = atof(data[ij].c_str());
				ASSERT(colstr[j].cs[ii] >= 0.0);
				//	printf("%g\t", colstr[j].cs[ii]);
			}
			//	printf("\n");
			spline(csTemp.data(), colstr[j].cs.data(), Ntemp, 2e31, 2e31, colstr[j].cs2d.data());
		}
		else
		{
			MakeCS( HFLines[j] );
			colstr[j].cs.clear();
			colstr[j].cs2d.clear();
		}

		data.resize(0);

		j++;
	}
	fclose(ioDATA);

	ASSERT( j == HFLines.size() );

	/* Discard no-longer needed nuclear data */
	for( long nelem = 0; nelem < LIMELM; nelem++ )
		abund.IsoAbn[nelem].rm_nuc_data();


#	if 0
	/* for debugging and developing only */
	/* calculating the luminosity for each isotope */
	for(int i = 0; i < HFLines.size(); i++)
	{
		N = dense.xIonDense[(*HFLines[i].Hi()).nelem()-1][(*HFLines[i].Hi()).IonStg()-1];
		Ne = dense.eden;

		h = 6.626076e-27;			/* erg * sec */
		c = 3e10;					/* cm / sec	 */
		k = 1.380658e-16;			/* erg / K   */

		upsilon = HyperfineCS(i);
				/*statistical weights must still be identified */
		q21 = COLL_CONST * upsilon / (phycon.sqrte * (*HFLines[i].Hi()).g());

		q12 = (*HFLines[i].Hi()).g()/ (*HFLines[i].Lo()).g() * q21 * exp(-1 * h * c * HFLines[i].EnergyWN / (k * phycon.te)); 

		x = Ne * q12 / (HFLines[i].Emis().Aul() * (1 + Ne * q21 / HFLines[i].Aul()));
		HFLines[i].xIntensity() = N * HFLines[i].Emis().Aul() * x / (1.0 + x) * h * c / (HFLines[i].EnergyAng() / 1e8);

	}
#	endif

	return;
}


/*HyperfineCS returns interpolated collision strength for transition index i */
double HyperfineCS( size_t i )
{
	double 	upsilon;

	DEBUG_ENTRY( "HyperfineCS()" );

	ASSERT( i >= 0. && i < HFLines.size() );

	if( colstr[i].cs.size() == 0 )
		return	HFLines[i].Coll().col_str();

	if( phycon.te <= csTemp[0] )
	{
		/* constant CS, if temperature below bounds of table */
		upsilon = colstr[i].cs[0];
	}
	else if( phycon.te >= csTemp[Ntemp-1] )
	{
		/* extrapolate, if temperature above bounds of table */
		int j = Ntemp - 1;
		double slope = log10(colstr[i].cs[j-1]/colstr[i].cs[j]) / log10(csTemp[j-1]/csTemp[j]);
		upsilon = log10(phycon.te/csTemp[j])*slope + log10(colstr[i].cs[j]);
		upsilon = exp10( upsilon);
	}
	else
	{
		splint( csTemp.data(), colstr[i].cs.data(), colstr[i].cs2d.data(), Ntemp, phycon.te, &upsilon );
	}

	return upsilon;
}

double HyperfineTspin( const TransitionProxy &t )
{
	DEBUG_ENTRY( "HyperfineTspin()" );

	double Tspin = phycon.te;

	if( (*t.Hi()).Pop() > SMALLFLOAT )
	{
		Tspin = TexcLine( t );

		/* this line MUST be non-zero
		 *
		 * H I 21cm does strongly mase in limit_compton_hi_t;
		 * in that sim pop ratio goes to unity for a float and
		 * TexcLine() returns zero */
		if( Tspin == 0. )
			Tspin = phycon.te;
	}

	ASSERT( Tspin > 0. or Tspin < 0. );

	return Tspin;
}
