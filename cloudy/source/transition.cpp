/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "version.h"
#include "dense.h"
#include "elementnames.h"
#include "lines.h"
#include "opacity.h"
#include "phycon.h"
#include "radius.h"
#include "rfield.h"
#include "rt.h"
#include "taulines.h"
#include "conv.h"
#include "lines_service.h"
#include "prt.h"

/*outline - adds line photons to reflin and outlin */
/*PutLine enter local line intensity into the intensity stack for eventual printout */
/*PutExtra enter and 'extra' intensity source for some line */
/*DumpLine print various information about an emission line vector, 
 * used in debugging, print to std out, ioQQQ */
/*TexcLine derive excitation temperature of line from contents of line array */
/*transition::Zero zeros out transition */
/*LineConvRate2CS convert down coll rate back into electron cs in case other parts of code need this for reference */
/*OccupationNumberLine - derive the photon occupation number at line center for any line */
/*MakeCS compute collision strength by g-bar approximations */
/*gbar1 compute g-bar collision strength using Mewe approximations */
/*gbar0 compute g-bar gaunt factor for neutrals */
/*emit_frac returns fraction of populations the produce emission */
/*chIonLbl use information in line array to generate a null terminated ion label in "Fe 2" */
/*chLineLbl use information in line transfer arrays to generate a line label */
/*PutCS enter a collision strength into an individual line vector */

map<std::string,std::vector<TransitionProxy> > blends;

/*outline - adds line photons to reflin and outlin */
void TransitionProxy::outline_resonance( ) const
{
	bool lgDoChecks = true;
	outline(Emis().ColOvTot(), lgDoChecks);
}

/*outline - adds line photons to reflin and outlin */
void TransitionProxy::outline( double nonScatteredFraction, 
										 bool lgDoChecks  ) const
{
	long int ip = ipCont()-1;

	DEBUG_ENTRY( "TransitionProxy::outline()" );

	if ( 0 && lgDoChecks)
	{
		// nothing to do if very small photon flux or below plasma frequency
		ASSERT( Emis().xIntensity() >= 0. );
		if( Emis().xIntensity() < SMALLFLOAT * EnergyErg() || EnergyErg() / EN1RYD <= rfield.plsfrq )
			return;

		ASSERT( Emis().FracInwd() >= 0. );
		ASSERT( radius.BeamInIn >= 0. );
		ASSERT( radius.BeamInOut >= 0. );
#ifndef NDEBUG
		double Ptot = Emis().Pesc_total();
		double PhotEmit = Emis().Aul()*Ptot*(*Hi()).Pop();
		// do not assert accuracy if close to bounds of fp precision
		double error = MAX2( SMALLFLOAT*1e3 , 3e-1*PhotEmit );
		// do not assert if do not have valid solution
		bool lgGoodSolution = conv.lgConvEden && conv.lgConvIoniz() &&
			conv.lgConvPops && conv.lgConvPres && conv.lgConvTemp;
#endif
		// see ticket #135 rt inconsistent results
		// this assert trips on a regular basis.  The fix is to rewrite
		// the RT and level population solvers so that the level population
		// and OTS rates are done simultaneously rather than sequentially
		// For now do not throw the assert on a release version - we know
		// about this problem
		ASSERT( t_version::Inst().lgRelease || !lgGoodSolution || 
			Ptot < 1e-3 ||  fp_equal_tol( Emis().xIntensity()/EnergyErg(), PhotEmit , error ));
	}

	bool lgTransStackLine = true;
	outline_base(Emis().dampXvel(), Emis().damp(), lgTransStackLine, ip, Emis().xIntensity()/EnergyErg(), Emis().FracInwd(),
					 nonScatteredFraction);
}

/*emit_frac returns fraction of populations the produce emission */
double emit_frac(const TransitionProxy& t)
{
	DEBUG_ENTRY( "emit_frac()" );

	if (! t.associated())
		return 0.0;

	ASSERT( t.ipCont() > 0 );

	/* collisional deexcitation and destruction by background opacities
	 * are loss of photons without net emission */
	double deexcit_loss = t.Coll().col_str() * dense.cdsqte + t.Emis().Aul()*t.Emis().Pdest();
	/* this is what is observed */
	double rad_deexcit = t.Emis().Aul()*(t.Emis().Pesc_total());
	return rad_deexcit/(deexcit_loss + rad_deexcit);
}

/*GetLineRec return rec coef*hnu*eden*n_ion for C, N, or O recombination lines from Dima's list,
 * also zero's line in master stack so not entered second time in later dump of all rec lines */
double GetLineRec(
	/* this is the number of the emission line in the stack of lines, on the C scale */
	long int ip,
	/* the multiplet wavelength */
  long int lWl)
{
	double GetLineRec_v;

	DEBUG_ENTRY( "GetLineRec()" );

	if( (long)(LineSave.RecCoefCNO[2][ip]+0.5) != lWl )
	{
		fprintf( ioQQQ, " GetLineRec called with incorrect wavelength.\n" );
		fprintf( ioQQQ, " index, call and get wl are %5ld%5ld%5ld\n",
		  ip, lWl, (long)(LineSave.RecCoefCNO[2][ip]+0.5) );
		cdEXIT(EXIT_FAILURE);
	}

	/* final product is vol emissivity in line */
	GetLineRec_v = LineSave.RecCoefCNO[3][ip]*dense.eden*
		dense.xIonDense[(long)(LineSave.RecCoefCNO[0][ip])-1][(long)(LineSave.RecCoefCNO[0][ip]-LineSave.RecCoefCNO[1][ip]+2)-1]*
		HC_ERG_ANG/LineSave.RecCoefCNO[2][ip];

	/* zero out rec coefficient so that not used again in master dump
	 * this routine cannot be called twice on same line */
	LineSave.RecCoefCNO[3][ip] = 0.;
	return( GetLineRec_v );
}

/*DumpLine print various information about an emission line vector, 
 * used in debugging, print to std out, ioQQQ */
void DumpLine(const TransitionProxy& t)
{
	DEBUG_ENTRY( "DumpLine()" );

	ASSERT( t.ipCont() > 0 );

	/* routine to print contents of line arrays */
	string chLbl = "DEBUG "+chLineLbl(t);

	fprintf( ioQQQ, 
		"%10.10s Te%.2e eden%.1e CS%.2e Aul%.1e Tex%.2e cool%.1e het%.1e conopc%.1e albdo%.2e\n", 
	  chLbl.c_str(), 
	  phycon.te, 
	  dense.eden, 
	  t.Coll().col_str(), 
	  t.Emis().Aul(), 
	  TexcLine(t), 
	  t.Coll().cool(), 
	  t.Coll().heat() ,
	  opac.opacity_abs[t.ipCont()-1],
	  opac.albedo[t.ipCont()-1]);

	fprintf( ioQQQ, 
		"Tin%.1e Tout%.1e Esc%.1e eEsc%.1e DesP%.1e Pump%.1e OTS%.1e PopL,U %.1e %.1e PopOpc%.1e\n", 
	  t.Emis().TauIn(), 
	  t.Emis().TauTot(), 
	  t.Emis().Pesc(), 
	  t.Emis().Pelec_esc(), 
	  t.Emis().Pdest(), 
	  t.Emis().pump(), 
	  t.Emis().ots(), 
	  (*t.Lo()).Pop(), 
	  (*t.Hi()).Pop() ,
	  t.Emis().PopOpc() );
	return;
}


/*OccupationNumberLine - derive the photon occupation number at line center for any line */
double OccupationNumberLine(const TransitionProxy& t)
{
	double OccupationNumberLine_v;

	DEBUG_ENTRY( "OccupationNumberLine()" );

	ASSERT( t.ipCont() > 0 );

	/* routine to evaluate line photon occupation number - 
	 * return negative number if line is maser */
	if( fabs(t.Emis().PopOpc()) > SMALLFLOAT )
	{
		/* the lower population with correction for stimulated emission */
		/* If the line mases PopOpc() will be negative, but Pesc() will > 1 as well,
		 * so the occupation number will be positive.
		 * However, this may be not the case if either value is stale. */
		OccupationNumberLine_v = ( (*t.Hi()).Pop() / (*t.Hi()).g() ) /
			( t.Emis().PopOpc() / (*t.Lo()).g() )  *
			(1. - t.Emis().Pesc());
	}
	else
	{
		OccupationNumberLine_v = 0.;
	}
	/* return value is not guaranteed to be positive - negative if
	 * line mases */
	return( OccupationNumberLine_v );
}

/*TexcLine derive excitation temperature of line from contents of line array */
double TexcLine(const TransitionProxy& t)
{
	double TexcLine_v;

	DEBUG_ENTRY( "TexcLine()" );

	/* routine to evaluate line excitation temp using contents of line array
	 * */
	if( (*t.Hi()).Pop() * (*t.Lo()).Pop() > 0. )
	{
		TexcLine_v = ( (*t.Hi()).Pop() / (*t.Hi()).g() )/( (*t.Lo()).Pop() / (*t.Lo()).g() );
		TexcLine_v = log(TexcLine_v);
		/* protect against infinite temp limit */
		if( fabs(TexcLine_v) > SMALLFLOAT )
		{
			TexcLine_v = - t.EnergyK() / TexcLine_v;
		}
	}
	else
	{
		TexcLine_v = 0.;
	}
	return( TexcLine_v );
}

/*chIonLbl use information in line array to generate a null terminated ion label in "Fe 2" */
string chIonLbl(const TransitionProxy& t)
{
	DEBUG_ENTRY( "chIonLbl()" );

	/* function to use information within the line array
	 * to generate an ion label, giving element and
	 * ionization stage
	 * */
	string chIonLbl_v;
	if( (*t.Hi()).nelem() <= 0 )
	{
		/* this line is to be ignored */
		chIonLbl_v = t.list()->states->chLabel();
		if( chIonLbl_v[0]=='\0' )
			chIonLbl_v = "Dumy";
	}
	else
	{
		chIonLbl_v = chIonLbl( (*t.Hi()).nelem(), (*t.Hi()).IonStg() );
	}
	/* chIonLbl is four char null terminated string */
	return chIonLbl_v;
}

string chIonLbl(const long& nelem, const long& IonStg)
{
	DEBUG_ENTRY( "chIonLbl()" );

	// NB - nelem passed here is expected to be on the physical scale, not C (so hydrogen = 1).
	ASSERT( nelem >= 1 && nelem <= LIMELM );
	ASSERT( IonStg >= 1 && IonStg <= nelem + 1 );
	/* ElmntSym.chElementSym is null terminated, 2 ch + null, string giving very
	 * short form of element name */
	string chIonLbl_v = elementnames.chElementSym[nelem-1];
	/* chIonStage is two char null terminated string, starting with "_1" */
	chIonLbl_v += elementnames.chIonStage[IonStg-1];
	return chIonLbl_v;
}

/*chLineLbl use information in line transfer arrays to generate a line label */
/* ContCreatePointers has test this with full range of wavelengths */
string TransitionProxy::chLabel() const
{
	DEBUG_ENTRY( "chLabel()" );

	string chSpecies;
	if( (*Hi()).nelem() < 1 && (*Hi()).IonStg() < 1 )
	{
		chSpecies = (*list()).chLabel;
	}
	else
	{
		chSpecies = chIonLbl( (*Hi()).nelem(), (*Hi()).IonStg() ); 
	}

	chSpecies.resize( NCHLAB-1, ' ' );

	/* NB this function is profoundly slow due to sprintf statement
	 * also - it cannot be evaluated within a write statement itself*/
	string chWavLen;
	sprt_wl(chWavLen, WLAng());
	return chSpecies + " " + chWavLen;
}

/*PutCS enter a collision strength into an individual line vector */
void PutCS(double cs, 
  const TransitionProxy& t)
{
	DEBUG_ENTRY( "PutCS()" );

	/* collision strength must be non-negative */
	ASSERT( cs > 0. );

	t.Coll().col_str() = (realnum)cs;

	return;
}

string GenerateTransitionConfiguration( const TransitionProxy &t )
{
	return	(*t.Lo()).chConfig() + " - " + (*t.Hi()).chConfig();
}

/*PutLine enter local line intensity into the intensity stack for eventual printout */
void PutLine(const TransitionProxy& t, const char *chComment, const char *chLabelTemp, const ExtraInten& extra)
{
	DEBUG_ENTRY( "PutLine()" );

	string chLabel;
	double xIntensity,
		other,
		xIntensity_in;
		
	/* routine to use line array data to generate input
	 * for emission line array */
	ASSERT( t.ipCont() > 0 );
		
	if (chLabelTemp)
	{
		chLabel = chLabelTemp;
		ASSERT (chLabel.length() <= NCHLAB-1);
	}

	/* if ipass=0 then we must generate label info since first pass
	 * gt.0 then only need line intensity data */
	if( LineSave.ipass == 0 )
	{
		if (!chLabelTemp)
		{
			/* these variables not used by linadd if ipass ne 0 */
			/* chIonLbl is function that generates a null terminated 4 char string, of form "C  2" */
			chLabel = chIonLbl(t);
		}
		xIntensity = 0.;
	}
	else
	{
		/* both the counting and integrating modes comes here */
		/* not actually used so set to safe value */
		chLabel = "";

		/* total line intensity or luminosity 
		 * these may not be defined in initial calls so define here */
		xIntensity = t.Emis().xIntensity() + extra.v;
	}

	/* initial counting case, where ipass == -1, just ignored above, call linadd below */
	
	/* ExtraInten is option that allows extra intensity (i.e., recomb)
	 * to be added to this line  with Call PutExtra( exta )
	 * in main lines this extra
	 * contribution must be identified explicitly */
	/*linadd(xIntensity,wl,chLabel,'i');*/
	/*lindst add line with destruction and outward */
	rt.fracin = t.Emis().FracInwd();
	lindst(t, extra,
			 chLabel.c_str(), 
			 /* this is information only - has been counted in cooling already */
			 't', 
			 /* do not add to outward beam, also done separately */
			 false,
			 chComment);
	rt.fracin = 0.5;

	/* inward part of line - do not move this away from previous lines
	 * since xIntensity is used here */
	xIntensity_in = xIntensity*t.Emis().FracInwd();
	ASSERT( xIntensity_in>=0. );
	linadd(xIntensity_in,t.WLAng(),"Inwd",'i',chComment);
	
	/* cooling part of line */
	other = t.Coll().cool();
	linadd(other,t.WLAng(),"Coll",'i',chComment);
	
	/* fluorescent excited part of line */
	double radiative_branching;
	enum { lgNEW = true };
	if (lgNEW)
	{
		// Improved two-level version of radiative branching ratio
		const double AulEscp = t.Emis().Aul()*(t.Emis().Pesc_total());
		// Would be better to include all outward transition processes from the
		// line, to cater for the general non-two-level case
		const double sinkrate = t.Emis().Aul()*t.Emis().Ploss() + t.Coll().ColUL( colliders );
		if (sinkrate > 0.0) 
		{
			radiative_branching = AulEscp/sinkrate;
		}
		else
		{
			radiative_branching = 0.;
		}
	}
	else
	{
		// This is the excitation ratio not the de-excitation ratio according
		// to its specification
		radiative_branching = (1.-t.Emis().ColOvTot());
	}

	other = (*t.Lo()).Pop() * t.Emis().pump() * radiative_branching * t.EnergyErg();
	linadd(other,t.WLAng(),"Pump",'i',chComment);
		

	/* heating part of line */
	other = t.Coll().heat();
	linadd(other,t.WLAng(),"Heat",'i',chComment);

	return;
}

/*PutLine enter local line intensity into the intensity stack for eventual printout */
void PutLine(const TransitionProxy& t, const char *chComment, const char *chLabelTemp)
{	
	DEBUG_ENTRY( "PutLine()" );
	ExtraInten extra_s;
	PutLine(t, chComment, chLabelTemp, extra_s);
}

/*PutLine enter local line intensity into the intensity stack for eventual printout */
void PutLine(const TransitionProxy& t,
	     const char *chComment)
{
	const char *chLabelTemp = NULL;
	DEBUG_ENTRY( "PutLine()" );
	ExtraInten extra_s;
	PutLine(t, chComment, chLabelTemp, extra_s);
}

void TransitionProxy::Junk() const
{

	DEBUG_ENTRY( "TransitionProxy::Junk()" );

		/* wavelength, usually in A, used for printout */
	WLAng() = -FLT_MAX;

	/* transition energy in wavenumbers */
	EnergyWN() = -FLT_MAX;

	/* array offset for radiative transition within continuum array 
	 * is negative if transition is non-radiative. */
	ipCont() = -10000;

	CollisionJunk( Coll() );

	/* set these equal to NULL first. That will cause the code to crash if
	 * the variables are ever used before being deliberately set. */
	ipEmis() = -1;
	
	setLo(-1);
	setHi(-1);
	return;
}

/*TransitionZero zeros out transition array at start of calculation, sets
 * optical depths to initial values */
void TransitionProxy::Zero() const
{

	DEBUG_ENTRY( "TransitionProxy::Zero()" );

	CollisionZero( Coll() );

	::Zero( *Lo() );
	::Zero( *Hi() );
	EmLineZero( Emis() );
	TauZero( Emis() );

	return;
}

/*LineConvRate2CS convert down coll rate back into electron cs in case other parts of code need this for reference */
void LineConvRate2CS( const TransitionProxy& t , realnum rate )
{

	DEBUG_ENTRY( "LineConvRate2CS()" );

	/* return is collision strength, convert from collision rate from 
	 * upper to lower, this assumes pure electron collisions, but that will
	 * also be assumed by anything that uses cs, for self-consistency */
	t.Coll().col_str() = rate * (*t.Hi()).g() / (realnum)dense.cdsqte;

	/* change assert to non-negative - there can be cases (Iin H2) where cs has
	 * underflowed to 0 on some platforms */
	ASSERT( t.Coll().col_str() >= 0. );
	return;
}

/*gbar0 compute g-bar gaunt factor for neutrals */
STATIC void gbar0(double ex, 
	  realnum *g)
{
	double a, 
	  b, 
	  c, 
	  d, 
	  y;

	DEBUG_ENTRY( "gbar0()" );

	/* written by Dima Verner
	 *
	 * Calculation of the effective Gaunt-factor by use of 
	 * >>refer	gbar	cs	Van Regemorter, H., 1962, ApJ 136, 906
	 * fits for neutrals
	 *  Input parameters: 
	 * ex - energy ryd - now K
	 * t  - temperature in K
	 *  Output parameter:
	 * g  - effective Gaunt factor
	 * */

	/* y = ex*157813.7/te */
	y = ex/phycon.te;
	if( y < 0.01 )
	{
		*g = (realnum)(0.29*(log(1.0+1.0/y) - 0.4/POW2(y + 1.0))/exp(y));
	}
	else
	{
		if( y > 10.0 )
		{
			*g = (realnum)(0.066/sqrt(y));
		}
		else
		{
			a = 1.5819068e-02;
			b = 1.3018207e00;
			c = 2.6896230e-03;
			d = 2.5486007e00;
			d = log(y/c)/d;
			*g = (realnum)(a + b*exp(-0.5*POW2(d)));
		}
	}
	return;
}

/*gbar1 compute g-bar collision strength using Mewe approximations */
STATIC void gbar1(double ex, 
	  realnum *g)
{
	double y;

	DEBUG_ENTRY( "gbar1()" );

	/*	*written by Dima Verner
	 *
	 * Calculation of the effective Gaunt-factor by use of 
	 * >>refer	gbar	cs	Mewe,R., 1972, A&A 20, 215
	 * fits for permitted transitions in ions MgII, CaII, FeII (delta n = 0)
	 * Input parameters: 
	 * ex - excitation energy in Ryd - now K
	 * te  - temperature in K
	 * Output parameter:
	 * g  - effective Gaunt factor
	 */

	/* y = ex*157813.7/te */
	y = ex/phycon.te;
	*g = (realnum)(0.6 + 0.28*(log(1.0+1.0/y) - 0.4/POW2(y + 1.0)));
	return;
}

/*MakeCS compute collision strength by g-bar approximations */
void MakeCS(const TransitionProxy& t)
{
	long int ion;
	double Abun, 
	  cs;
	realnum
	  gbar;

	DEBUG_ENTRY( "MakeCS()" );

	/* routine to get cs from various approximations */

	/* check if abundance greater than 0 */
	ion = (*t.Hi()).IonStg();

	//This is the oscillator strength limit where larger values are assumed to be for allowed transitions.
	const double gfLimit = 1e-8;

	Abun = dense.xIonDense[ (*t.Hi()).nelem() -1 ][ ion-1 ];
	if( Abun <= 0. )
	{
		gbar = 1.;
	}
	else if( t.Emis().gf() >= gfLimit )
	{
		/* check if neutral or ion */
		if( ion == 1 )
		{
			/* neutral - compute gbar for eventual cs */
			gbar0(t.EnergyK(),&gbar);
		}
		else
		{
			/* ion - compute gbar for eventual cs */
			gbar1(t.EnergyK(),&gbar);
		}
	}
	else
	{
		//Mewe72 provides a gbar estimate for forbidden transitions.
		gbar = 0.15;
	}

	/* above was g-bar, convert to cs */
	cs = gbar*(14.5104/WAVNRYD)*t.Emis().gf()/t.EnergyWN();

	/* stuff the cs in place */
	t.Coll().col_str() = (realnum)cs;
	// t.Coll().is_gbar() = 1;
	return;
}

void TransitionProxy::AddLine2Stack() const
{
	DEBUG_ENTRY( "AddLine2Stack()" );

	ASSERT( lgLinesAdded == false );

	size_t newsize = m_list->Emis.size()+1; 
	m_list->Emis.resize(newsize);
	ipEmis() = newsize-1;
	this->resetEmis();
}

void TransitionProxy::AddLoState() const
{
	DEBUG_ENTRY( "AddLoState()" );

	ASSERT( !lgStatesAdded );

	m_list->states->addone();

	setLo(m_list->states->size()-1);
}

void TransitionProxy::AddHiState() const
{
	DEBUG_ENTRY( "AddHiState()" );

	ASSERT( !lgStatesAdded );

	m_list->states->addone();

	setHi(m_list->states->size()-1);
}
