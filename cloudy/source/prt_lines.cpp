/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines main routine to put emission line intensities into line stack,
 * calls lineset1, 2, 3, 4 */
/*FindStrongestLineLabels find strongest lines contributing to point in continuum energy mesh, output in some save commands */
#include "cddefines.h"
#include "taulines.h"
#include "thermal.h"
#include "yield.h"
#include "ipoint.h"
#include "ionbal.h"
#include "cddrive.h"
#include "trace.h"
#include "prt.h"
#include "rt.h"
#include "rfield.h"
#include "phycon.h"
#include "iso.h"
#include "hyperfine.h"
#include "hydrogenic.h"
#include "atmdat.h"
#include "lines.h"
#include "radius.h"
#include "dense.h"
#include "lines_service.h"
#include "mole.h"
#include "oxy.h"
#include "continuum.h"
#include "fe.h"
#include "species.h"
#include "generic_state.h"

STATIC void lines_iron_Ka();

STATIC void getTransition(const LineID& line, TransitionProxy& tr)
{
	if( LineSave.ipass == 0 )
	{
		long id = LineSave.findline(line);
		if( id <= 0 )
		{
			fprintf( ioQQQ, "getTransition: the following line was not found: \"%s\" %.3f\n",
					 line.chLabel.c_str(), line.wave );
			cdEXIT(EXIT_FAILURE);
		}
		tr = LineSave.lines[id].getTransition();
		if( !tr.associated() )
		{
			fprintf( ioQQQ, "getTransition: the following line is not associated with a transition: \"%s\" %.3f\n",
					 line.chLabel.c_str(), line.wave );
			cdEXIT(EXIT_FAILURE);
		}
	}
}

void lines()
{
	long int i, 
	  ipnt,
	  nelem;
	double f2, sum; 

	DEBUG_ENTRY( "lines()" );

	/* LineSave.ipass
	 * -1 - space not yet allocated - just count number of lines entered into stack
	 *  0 - first call with space allocated - must create labels and add in wavelengths
	 * +1 - later calls - just add intensity 
	 */

	/* major routines used here:
	 *
	 * PutLine( tarray )
	 * this uses information in tarray to give various
	 * contributions to lines, and their intensities
	 *
	 * PutExtra( extra )
	 * extra is some extra intensity to add to the line
	 * it will go into the totl contribution put out by PutLine,
	 * and this contribution should be indicated by independent
	 * call to linadd
	 * */

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " lines called\n" );
	}

	/* total luminosity radiated by this model - will be compared with energy in incident
	 * continuum when calculation is complete */
	thermal.power += thermal.htot*radius.dVeffAper;

	/* remember the total free-free heating */
	fixit("need to get rid of brems_heat_total"); 
	// get rid of brems_heat_total entirely (only net heating is added into stack, so use net here to avoid nonsense ratios elsewhere)  
	//thermal.FreeFreeTotHeat += CoolHeavy.brems_heat_total*radius.dVeffAper;
	thermal.FreeFreeTotHeat += thermal.heating(0,11)*radius.dVeffAper;

	/* total Compton heating - cooling */
	rfield.comtot += rfield.cmheat*radius.dVeffAper;
	thermal.totcol += thermal.ctot*radius.dVeffAper;

	/* up up induced recombination cooling */
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		hydro.cintot += iso_sp[ipH_LIKE][nelem].RecomInducCool_Rate*radius.dVeffAper;
	}

	/* nsum is line pointer within large stack of line intensities */
	LineSave.nsum = 0;
	LineSave.nComment = 0;

	/* this is used by lindst to proportion inward and outward.  should be 50% for
	 * optically thin line.  putline sets this to actual value for particular line
	 * and calls lindst then rests to 50% */
	rt.fracin = 0.5;

	/* last arg in call to lindst and linadd is info on line
	 * info is char variable indicating type of line this is
	 * 'c' cooling
	 * 'h' heating
	 * 'i' information only
	 * 'r' recombination line
	 *
	 * all components of lines are entered into the main line stack here
	 * when printing, filters exist to not print Inwd component */

	/* initialize routine that can generate pointers for forbidden lines,
	 * these are lines that are not transferred otherwise,
	 * in following routines there will be pairs of calls, first to
	 * PntForLine to get pointer, then lindst to add to stack */
	PntForLine(0.,"FILL",&i);

	/* evaluate rec coefficient for rec lines of C, N, O
	 * some will be used in LineSet2 and then zeroed out,
	 * others left alone and used below */
	t_ADfA::Inst().rec_lines(phycon.te,LineSave.RecCoefCNO);

	/* put in something impossible in element 0 of line stack */
	linadd(0.f,0,"zero",'i' , "null placeholder");

	/* this is compared with true volume in final.  The number can't
	 * actually be unity since this would overflow on a 32 bit machine */
	/* integrate the volume as a sanity check */
	linadd( 1.e-10 , 1 , "Unit" , 'i' , "unit integration placeholder");
	static long int ipOneAng=-1;
	if( LineSave.ipass<0 )
		ipOneAng = ipoint( RYDLAM );
	lindst( 1.e-10 , 1. , "UntD" , ipOneAng , 'i' , false,"unit integration placeholder");

	/* initial set of general properties */
	lines_general();

	/* do all continua */
	lines_continuum();

	/* information about grains */
	lines_grains();

	/* update all satellite lines */
	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
		iso_satellite_update(nelem);

	/* do all hydrogenic ions */
	lines_hydro();

	/* enter He-iso sequence lines */
	lines_helium();

	/* next come the extra Lyman lines */
	i = StuffComment( "extra Lyman" );
	linadd( 0., (realnum)i , "####", 'i' ,
			  "extra Lyman lines");
	
	for( long ipISO=ipH_LIKE; ipISO < NISO; ++ipISO )
	{
		/* loop over all iso-electronic sequences */
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( ! dense.lgElmtOn[nelem] )
				continue;
			for( long ipHi=iso_sp[ipISO][nelem].numLevels_max; ipHi < iso_ctrl.nLyman_max[ipISO]; ipHi++ )
			{
				if (ExtraLymanLines[ipISO][nelem][ipExtraLymanLines[ipISO][nelem][ipHi]].ipCont() > 0)
					PutLine(ExtraLymanLines[ipISO][nelem][ipExtraLymanLines[ipISO][nelem][ipHi]],
							  "extra Lyman line");
			}
		}
	}

#if	0
	/* This is Ryan's code for dumping lots of Helium lines according to
	 * quantum number rather than wavelength, principally for comparison with Rob
	 * Bauman's code. */
	if( iteration > 1 )
	{
		fprintf( ioQQQ,"ipHi\tipLo\tnu\tlu\tsu\tnl\tll\tsl\tWL\tintens\n" );
		for( long ipHi=5; ipHi<= iso_sp[ipHE_LIKE][ipHELIUM].numLevels_local - iso_sp[ipHE_LIKE][ipHELIUM].nCollapsed_local; ipHi++ )
		{
			for( long ipLo=0; ipLo<ipHi; ipLo++ )
			{
				if( iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHi,ipLo).ipCont() > 0 )
				{
					double relint, absint;

					if( cdLine("He 1", 
						iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHi,ipLo).WLAng(),
						&relint, &absint ) )
					{
						//iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHi,ipLo).Hi()->chLabel

						//if( iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHi,ipLo).WLAng() < 1.1E4 &&
						//	iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHi,ipLo).WLAng() > 3.59E3 &&
						//	ipLo!=3 && ipLo!=4 && relint >= 0.0009 )
						long n=iso_sp[ipHE_LIKE][ipHELIUM].st[ipHi].n();
						long np=iso_sp[ipHE_LIKE][ipHELIUM].st[ipLo].n();
						long l=iso_sp[ipHE_LIKE][ipHELIUM].st[ipHi].l();
						long lp=iso_sp[ipHE_LIKE][ipHELIUM].st[ipLo].l();
						long s=iso_sp[ipHE_LIKE][ipHELIUM].st[ipHi].S();
						long sp=iso_sp[ipHE_LIKE][ipHELIUM].st[ipLo].S();
						if ( ((l== lp+1 || l == lp-1) && l==n-3 && s== sp))//|| ((l == n-1) && s==sp ))
						{
							fprintf( ioQQQ,"lines %li\t%li\t%li\t%li\t%li\t%li\t%li\t%li\t%e\t%e\n",
								ipHi,
								ipLo,
								iso_sp[ipHE_LIKE][ipHELIUM].st[ipHi].n(),
								iso_sp[ipHE_LIKE][ipHELIUM].st[ipHi].l(),
								iso_sp[ipHE_LIKE][ipHELIUM].st[ipHi].S(),
								iso_sp[ipHE_LIKE][ipHELIUM].st[ipLo].n(),
								iso_sp[ipHE_LIKE][ipHELIUM].st[ipLo].l(),
								iso_sp[ipHE_LIKE][ipHELIUM].st[ipLo].S(),
								iso_sp[ipHE_LIKE][ipHELIUM].trans(ipHi,ipLo).WLAng(),
								relint );
						}
					}
				}
			}
		}
	}
#endif

	// these must come before the old level 1 or 2 lines.  we now have the option
	// to use the external database by default, or turn it off (set chianti off)
	// and fall back to the old line treatment.  When database is off those lines
	// do not exist.  But when database is on the level 1 lines are still evaluated
	// but with the ion abundance set to zero.  If the level 1 lines were entered
	// into the stack then searches for the line with cdLine would turn up the level 1
	// line, which would be zero.
	// The long term goal is to have all lines be external databases & rm level 1&2 liens

	/* external database lines */
	i = StuffComment( "database lines" );
	linadd( 0., (realnum)i , "####", 'i' ,
		"database lines");
	for (int ipSpecies=0; ipSpecies < nSpecies; ++ipSpecies)
	{
		for( EmissionList::iterator em=dBaseTrans[ipSpecies].Emis().begin();
			  em != dBaseTrans[ipSpecies].Emis().end(); ++em)
		{
			if( (*em).Tran().ipCont() > 0)
			{
				/* Set the comment to the database name followed by the
				 * lower and upper indices of the transition energy levels
				 * as they appear in the database input file. */
				string chComment = "";
				if( LineSave.ipass == 0 )
				{
					chComment  = dBaseSpecies[ipSpecies].database
						+ ", " + (*em).Tran().getComment();
				}
				PutLine((*em).Tran(), chComment.c_str(), dBaseSpecies[ipSpecies].chLabel);
			}
		}
	}

	/* level 1 lines */
	ipnt = StuffComment( "level 1 lines" );
	linadd( 0., (realnum)ipnt , "####", 'i',
		"start level 1 lines" );

	/* do iron Ka */
	lines_iron_Ka();

	/* next come some recombination lines */
	i = StuffComment( "recombination" );
	linadd( 0., (realnum)i , "####", 'i' ,
		"recombination lines");

	/***********************************************************************
	 * large number of C, N, and O recombination lines                     *
	 *************************************************************************/

	for( i=0; i < NRECCOEFCNO; i++ )
	{
		string chLabel;
		/* generate label for the line if ipass is -1 or 0 - saved in arrays
		 * so no need to do this during production */
		if( LineSave.ipass <= 0 )
		{
			chLabel = chIonLbl( LineSave.RecCoefCNO[0][i], long(LineSave.RecCoefCNO[0][i]-LineSave.RecCoefCNO[1][i]+1.01) );
		}
		else
			chLabel = "    ";

		/* number of rec per unit vol
		 * do not predict these simple reccombination intensities at high densities
		 * since lines neglect collisional deexciation and line optical depth effects.
		 * They were not intended for high densities or column densities.
		 * As a result they become unphysically bright at high densities and
		 * violate the black body limit.  There would be major
		 * energy conservation problems if they were added in the outward beam in
		 * dense simulations.
		 * */
		if( dense.eden < 1e8  )
		{
			nelem = (long)LineSave.RecCoefCNO[0][i]-1;
			long int ion = (long)(LineSave.RecCoefCNO[0][i]-LineSave.RecCoefCNO[1][i]+2)-1;
			f2 = LineSave.RecCoefCNO[3][i]*dense.eden*
				dense.xIonDense[nelem][ion];

			/* convert to intensity */
			f2 = f2*HC_ERG_ANG/LineSave.RecCoefCNO[2][i];
		}
		else
		{
			f2 = 0.;
		}
		/* stuff it into the stack */
		PntForLine(LineSave.RecCoefCNO[2][i], chLabel.c_str(), &ipnt);
		lindst(f2,wlAirVac(LineSave.RecCoefCNO[2][i]), chLabel.c_str(), ipnt, 'r', true,
			"recombination line");
	}

	/* next come the atom_level2 lines */
	i = StuffComment( "level2 lines" );
	linadd( 0., (realnum)i , "####", 'i' ,
		"level2 lines");

	/* add in all the other level 2 wind lines
	 * Dima's 6k lines */
	double ExtraCool = 0.;
	double BigstExtra = 0.;
	for( i=0; i < nWindLine; i++ )
	{
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
		{
			PutLine(TauLine2[i],"level 2 line");
			if( TauLine2[i].Coll().cool() > BigstExtra )
			{
				BigstExtra = TauLine2[i].Coll().cool();
				thermal.ipMaxExtra = i+1;
			}
			ExtraCool += TauLine2[i].Coll().cool();
		}
	}
	/* keep track of how important this is */
	thermal.GBarMax = MAX2(thermal.GBarMax,(realnum)(ExtraCool/thermal.ctot));

	/* next come the hyperfine structure lines */
	i = StuffComment( "hyperfine structure" );
	linadd( 0., (realnum)i , "####", 'i' ,
		"hyperfine structure lines");

	/* this is total cooling due to all HF lines */
	linadd( hyperfine.cooling_total, 0., "hfin", 'i' ,
		"total cooling all hyperfine structure lines");

	/* remember largest local cooling for possible printout in comments */
	hyperfine.cooling_max = (realnum)MAX2(hyperfine.cooling_max,hyperfine.cooling_total/thermal.ctot);

	/* the hyperfine lines */
	for( size_t ih=0; ih < HFLines.size(); ih++ )
	{
		PutLine(HFLines[ih],
			"hyperfine structure line");
	}

	/* next come the inner shell fluorescent lines */
	i = StuffComment( "inner shell" );
	linadd( 0., (realnum)i , "####", 'i' ,
		"inner shell lines");

	// Opacity Project satellite lines from data / UTA
	// the heat component of these lines is the heat per pump, a constant
	for( size_t i=0; i < UTALines.size(); i++ )
	{
		PutLine(UTALines[i], "OP satellite");
	}

	/* the group of inner shell fluorescent lines
	 * >>refer all	auger	Kaastra, J. S. & Mewe, R. 1993, A&AS, 97, 443-482
	 * http://adsabs.harvard.edu/abs/1993A%26AS...97..443K */
	for( i=0; i < t_yield::Inst().nlines(); ++i )
	{
		double xInten = 
			/* density of parent ion, cm-3 */
			dense.xIonDense[t_yield::Inst().nelem(i)][t_yield::Inst().ion(i)] *
			/* photo rate per atom per second, s-1 */
			ionbal.PhotoRate_Shell[t_yield::Inst().nelem(i)][t_yield::Inst().ion(i)][t_yield::Inst().nshell(i)][0]*
			/* fluor yield - dimensionless */
			t_yield::Inst().yield(i) *
			/* photon energy - ryd, converted into ergs */
			t_yield::Inst().energy(i) * EN1RYD;

		/* create label if initializing line stack */
		string chLabel;
		if( LineSave.ipass == 0 )
		{
			chLabel = chIonLbl( t_yield::Inst().nelem(i)+1, t_yield::Inst().ion_emit(i)+1 );
#			if 0
			/* only print yields for atoms */
			if( t_yield::Inst().ion(i) == 0 && t_yield::Inst().nelem()(i) == ipIRON )
			fprintf(ioQQQ,"DEBUGyeild\t%s\t%.3f\t%.3e\n",
				/* line designation, energy in eV, yield */
				chLabel.c_str() , t_yield::Inst().energy()(i)*EVRYD, t_yield::Inst().yield(i) );
#			endif
		}

		/* the group of inner shell fluorescent lines */
		lindst(
			/* intensity of line */
			xInten,
			/* wavelength of line in Angstroms */
			(realnum)RYDLAM / t_yield::Inst().energy(i),
			/* label */
			chLabel.c_str() ,
			/* continuum array offset for line as set in ipoint */
			t_yield::Inst().ipoint(i), 
			/* type of line - count as a recombination line */
			'r',
			/* include line in continuum? */
			true ,
			"inner shell line");
	}

	/* now do all molecules - do last since so many H2 lines */
	lines_molecules();

	SpeciesPseudoContAccum();
	SpeciesBandsAccum();

	/* blends of other lines */
	i = StuffComment( "miscellaneous" );

	linadd( 0., (realnum)i , "####", 'i' ,	"miscellaneous");

	/* add blends */
	for( size_t i=0; i < prt.blend.size(); ++i )
	{
		if( prt.blend[i].lgQuiet && LineSave.ipass < 0 )
		{
			for( size_t j=0; j < prt.blend[i].component.size(); ++j )
			{
				string speciesLabel;
				long nelem, ion;
				parsespect(prt.blend[i].component[j].chLabel.c_str(), nelem, ion);
				if( nelem < 0 || ion <= 0 )
					// assume molecular species
					speciesLabel = prt.blend[i].component[j].chLabel;
				else
				{
					char chLabelChemical[10];
					// NB parsespect() returns nelem on C scale, but ion on fortran scale!
					makeChemical(chLabelChemical, nelem, ion-1);
					speciesLabel = chLabelChemical;
				}
				vector<genericState> v = matchGeneric( speciesLabel, false );
				if( v.size() == 1 )
				{
					if( v[0].sp->lines == nullptr )
					{
						prt.blend[i].lgIgnore = true;
						break;
					}
				}
				else
				{
					prt.blend[i].lgIgnore = true;
					break;
				}
			}
		}
		if( !prt.blend[i].lgIgnore )
		{
			LinSv *UserBlnd = linadd(0.0,prt.blend[i].wave,prt.blend[i].chLabel.c_str(),'i',"Blend" );
			if (UserBlnd)
			{
				for( size_t j=0; j < prt.blend[i].component.size(); ++j )
					UserBlnd->addComponent(prt.blend[i].component[j]);
				if( prt.blend[i].wave < 0_r )
					UserBlnd->setBlendWavl();
			}
		}
	}

	/*************Calcium *******************/
	if( atmdat.lgdBaseSourceExists[ipCALCIUM][1] )
	{
		double eff = dense.eden*dense.xIonDense[ipCALCIUM][2]*5.4e-21/(phycon.te/
		  phycon.te10/phycon.te10);
		PntForLine( 3933., "Ca2R", &ipnt);
		lindst(eff, 3933., "Ca2R", ipnt, 't', false, "recombination emission" );
	}

	/*** CARBON ***/
	/* Recombination and Pump lines */
	double pump = 0, fac = 0;

	/* C 1 1656 */
	if( atmdat.lgdBaseSourceExists[ipCARBON][0] )
	{
		static TransitionProxy lineCIa;
		getTransition(LineID("C  1",1657.01), lineCIa);
		/* >>chng 97 may 02, added better rec coefficient
		 * C I 1656 recombination, all agents */
		double rec = LineSave.ipass <= 0 ? 0.0 :
			GetLineRec(3,1657)*emit_frac(lineCIa);
		PntForLine( 1656., "C 1R", &ipnt);
		lindst(rec, 1656., "C 1R", ipnt, 't', false, "recombination" );
	}

	/* C 1 9850 */
	if( atmdat.lgdBaseSourceExists[ipCARBON][0] )
	{
		/* >>chng 97 may 02, added better rec coefficient
		 * C 1 9850, recombination contribution rec coefficient from
		 * >>refer	C1	rec	Escalante, Vladimir, & Victor, G.A., 1990, ApJS 73, 513.
		 * r9850 is correction for collisional deexcitation as in carb cool
		 * >>chng 97 aug 2, had factor of rec, changed to r9850, this
		 * was a big mistake
		 * >>chng 12 nov 11 Now using Escalante 1990 data */
		double A,B,C,t4;
		t4 = 1e-4*phycon.te;
		A = 3.10e-17;
		B = 0.25;
		C = -0.41;
		double c19850WL = 9850.26;
		double recCoeff = A*pow(t4,-1*B*(1+C*log10(t4)));

		double volEmis = recCoeff*dense.eden*dense.xIonDense[ipCARBON][1]*HC_ERG_ANG/c19850WL;

		static TransitionProxy lineCIb;
		getTransition(LineID("C  1",9850.26), lineCIb);
		double rec = LineSave.ipass <= 0 ? 0.0 :
			volEmis*emit_frac(lineCIb); // 9850.26
		PntForLine( wlAirVac(9850.), "C 1R", &ipnt);
		lindst(rec, wlAirVac(9850.), "C 1R", ipnt, 't', false, "recombination" );
	}

	/* C 2 2326 */
	if( atmdat.lgdBaseSourceExists[ipCARBON][1] )
	{
		linadd(ionbal.PhotoRate_Shell[ipCARBON][0][1][0]*dense.xIonDense[ipCARBON][0]*0.1*8.6e-12,wlAirVac(2326.),"C 2H",'i' ,
			"photoproduction, Hofmann and Trefftz");
		// see 1980A&A....82..256H, 1983A&A...126..415H
	}

	/* C 2 1335 */
	if( atmdat.lgdBaseSourceExists[ipCARBON][1] )
	{
		static TransitionProxy lineCIIa;
		getTransition(LineID("C  2",1335.71), lineCIIa);

		/* >>chng 97 may 02, better rec coef */
		/* >>chng 02 jul 01, add function to return emission probability */
		double rec = LineSave.ipass <= 0 ? 0.0 :
			GetLineRec(11,1335)*emit_frac(lineCIIa); // 1335.71
		PntForLine( 1335., "C 2R", &ipnt);
		lindst(rec, 1335., "C 2R", ipnt, 't', false, "recombination" );
	}

	/* C 2 3920 */
	if( atmdat.lgdBaseSourceExists[ipCARBON][1] )
	{
		/* the CII 3918.98/3920.68 and 6578.05/6582.88 multiplets,
		 * contributions by both continuum pumping through XUV line
		 * and recombination */
		/* this is the driving line, pump is photons cm^-3 s^-1 */
		if( nWindLine > 0 )
		{
			pump = TauLine2[186].Emis().pump()*TauLine2[186].Emis().PopOpc();
		}
		else
		{
			pump = 0.;
		}

		// pumped 3920, count as recomb since does remove energy
		PntForLine(3920.,"C 2P",&ipnt);
		lindst(pump*0.387 * 5.08e-12/(1.+dense.eden/1e12) ,3920,"C 2P",ipnt,'r',true ,
			"this is only pumped, no recombination part");

		/* convert UV pump rate to intensity with branching ratio and hnu */
		pump *= 0.305 * 0.387 * 3.02e-12;

		linadd(pump/(1.+dense.eden/1e12),wlAirVac(6580.),"C 2P",'i',
			"excitation by pumping" );

		static TransitionProxy lineCIIc;
		getTransition(LineID("C  2",6578.05), lineCIIc);
		/* C 2 6580 */
		double rec = LineSave.ipass <= 0 ? 0.0 :
			GetLineRec(8, 6580 )*emit_frac(lineCIIc);
		PntForLine( wlAirVac(6580.), "C 2R", &ipnt);
		lindst(rec/(1.+dense.eden/1e12), wlAirVac(6580.), "C 2R", ipnt, 't', false, "recombination" );
	}

	/* C 3 977
	 * recombination contribution from nussbaumer and story 84 */
	if( atmdat.lgdBaseSourceExists[ipCARBON][2] )
	{
		static TransitionProxy lineCIIIc;
		getTransition(LineID("C  3",977.020), lineCIIIc);
		/* >>chng 02 jul 01, add function to compute emission fraction */
		double rec = LineSave.ipass <= 0 ? 0.0 :
			GetLineRec(179,977)*emit_frac(lineCIIIc);
		PntForLine( 977., "C 3R", &ipnt);
		lindst(rec, 977., "C 3R", ipnt, 't', false, "dielectronic recombination" );
	}

	/* C 3 1909 */
	if( atmdat.lgdBaseSourceExists[ipCARBON][2] )
	{
		static TransitionProxy lineCIII;
		getTransition(LineID("C  3",1908.73), lineCIII);
		/* >>chng 02 jul 01, add function to compute emission fraction */
		double corr =  LineSave.ipass <= 0 ? 0.0 :
			emit_frac(lineCIII); // 1908.73
		fac = dense.eden*dense.xIonDense[ipCARBON][3]/(phycon.te/phycon.te10);

		PntForLine( 1909., "C 3R", &ipnt );
		lindst( 3.1e-19*fac*corr, 1909., "C 3R", ipnt, 't', false,
			"recombination from Storey" );

		lindst(ionbal.PhotoRate_Shell[ipCARBON][1][1][0]*dense.xIonDense[ipCARBON][1]*0.62*corr*1.05e-11,1909,"C 3H",ipnt,'t',false,
			"production following relax following inner shell photoionization" );
	}

	/* C 3 1175 */
	if( atmdat.lgdBaseSourceExists[ipCARBON][2] )
	{
		static TransitionProxy lineCIIIb;
		getTransition(LineID("C  3",1175.71), lineCIIIb);
		/* >>chng 97 may 02, better rec ocef */
		/* >>chng 02 jul 01, add function to compute emission fraction */
		double rec = LineSave.ipass <= 0 ? 0.0 :
			GetLineRec(178,1176)*emit_frac(lineCIIIb); // 1175.71
		PntForLine( 1175., "C 3R", &ipnt );
		lindst(rec, 1175., "C 3R", ipnt, 't', false, "dielectronic recombination" );
	}

	/* C 4 1549 */
	if( atmdat.lgdBaseSourceExists[ipCARBON][3] )
	{
		static TransitionProxy lineCIV;
		getTransition(LineID("C  4",1548.19), lineCIV);
		/* recombination C 4 1549 from C 5
		 * >>chng 97 may 02, better rec coef */
		/* >>chng 02 jul 01, add function to compute emission fraction */
		double rec = LineSave.ipass <= 0 ? 0.0 :
			GetLineRec(25,1549)*emit_frac(lineCIV); // 1548.19
		PntForLine( 1549., "C 4R", &ipnt );
		lindst(rec, 1549., "C 4R", ipnt, 't', false, "recombination" );
	}

	/*** End of Carbon ***/

	/*** NITROGEN ***/
	/* Recombination and Pump lines */

	/***********************/
	/****** N 1 5199 ******/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipNITROGEN][0] )
	{
		/**** Terry's addition, recombination from N+ **************/
		/* rate coefficient (cm3 s-1) from Table 3 of
		 *>>refer	NI	rec	Pequignot, D., Petijean, P. & Boisson, C. 1991, A&A, 251, 680 */
		double eff_recrate_2D = 1.108e-13 * pow(phycon.te*1e-4, -0.6085) /
			(1. - 0.0041 * pow(phycon.te*1e-4, -0.3975) );
		double eff_recrate_2P = 0.659e-13 * pow(phycon.te*1e-4, -0.6158);
		double fac_n1;
		if( dense.xIonDense[ipNITROGEN][0] > 0. )
			fac_n1 = dense.eden * dense.xIonDense[ipNITROGEN][1] / dense.xIonDense[ipNITROGEN][0];
		else
			fac_n1 = 0.;

		// assume levels are populated according to statistical weight
		double rec14 = eff_recrate_2P * fac_n1 * 2./6.;
		double rec15 = eff_recrate_2P * fac_n1 * 4./6.;
		double rec13 = eff_recrate_2D * fac_n1 * 4./10.;
		double rec12 = eff_recrate_2D * fac_n1 * 6./10.;

		static TransitionProxy lineNI5199;
		getTransition(LineID("N  1",5197.90), lineNI5199);

		double emit_frac_5197 = LineSave.ipass <= 0 ? 0.0 : 
			emit_frac(lineNI5199); // 5197.90

		// estimate of recombination contribution to intensity of [NI] 5199
		double rec = (rec12+rec13+0.9710*rec14+0.9318*rec15) * dense.xIonDense[ipNITROGEN][0] *
			3.82e-12*emit_frac_5197;

		/**** Terry's addition **************/
		PntForLine( wlAirVac(5199.), "N  1", &ipnt);
		lindst(rec, wlAirVac(5199.), "N 1R", ipnt, 't', true,
			"estimate of production by recombination" );

		/* this is upper limit to production of 5200 by chemistry - assume every photo dissociation
		 * populates upper level
		 * co.nitro_dissoc_rate is the total N photo dissociation rate, cm-3 s-1 */
		// count as recombination since removes energy but is not coolng
		double chem = mole.dissoc_rate("N") * 3.82e-12 * emit_frac_5197;
		lindst( chem, wlAirVac(5199.), "N 1C", ipnt, 'r', true,
			"upper limit to production by chemistry" );

		/* this is upper limit to production of 5200 by charge transfer -
		 * atmdat.HCharExcRecTo_N0_2D is the rate coefficient (cm3 s-1) for N+(3P) + H0 -> H+ + N0(2D) */
		double ctRate = atmdat.HCharExcRecTo_N0_2D*dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipNITROGEN][1] *
			3.82e-12 * emit_frac_5197;
		lindst( ctRate, wlAirVac(5199.), "N 1T", ipnt, 'r', true,
			"upper limit to production by charge transfer" );

		//	// estimate of pumping contribution to 5200 doublet
		//	lindst( nitro.pump5199, 5199, "Pump", ipnt, 'r', true,
		//		"estimate of production by FUV pumping" );
	}

	/***********************/
	/****** N 2 6584 ******/
	/*********************/

	/* May be a combination of lines 1D - 3P  */
	double efficn2 = 4e-3/(4e-3 + 5.18e-6*dense.eden/phycon.sqrte);
	double rec = 8e-22/(phycon.te70/phycon.te03/phycon.te03)*efficn2;
	PntForLine( 6584., "N 2R", &ipnt );
	lindst( rec*dense.xIonDense[ipNITROGEN][2]*dense.eden, 6584., "N 2R", ipnt, 't', false,
		"N 2 6584 alone, recombination" );

	/***********************/
	/****** N 2 5755 ******/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipNITROGEN][1] )
	{
		/* helium charge transfer from
		>>refer	n2	CT	Sun Sadeghpour, Kirby Dalgarno and Lafyatis, CfA preprint 4208 */
		double ctRate = 1.8e-11*dense.xIonDense[ipHELIUM][0]*dense.xIonDense[ipNITROGEN][2]*1.146/(1.146 +
		  0.87*dense.cdsqte)*3.46e-12;

		PntForLine(5755.,"N  2",&ipnt);

		lindst(ctRate,wlAirVac(5755.),"N 2T",ipnt,'r',true,
			"charge transfer" );

		static TransitionProxy lineNII5755;
		getTransition(LineID("N  2",5754.59), lineNII5755);

		/* >>chng 01 jul 09, add recombination contribution to 5755 */
		/* >>refer	n2	rec	Liu, X.W., Storey, P.J., Barlow, M.J., Danziger, I.J., Cohen, M.,
		 * >>refercon	& Bryce, M., 2000, MNRAS, 312, 585 */
		/* they give intensity in terms of hbeta intensity as their equation 1 */
		if( dense.xIonDense[ipHYDROGEN][1] > SMALLFLOAT )
		{
			/* this test on >0 is necessary because for sims with no H-ionizing radiation
			 * the H+ density is initially zero */
			/* H beta recombination, assuming old case B, needed since HS tables have
			 * only a narrow temperature range - at this point units are ergs cm^3 s-1 */
			double HBeta = (exp10(-20.89 - 0.10612*POW2(phycon.alogte - 4.4)))/phycon.te;

			/* now convert to ergs cm-3 s-1
			 * >>chng 05 mar 17, this step was missing, so recombination intensity off by density squared,
			 * bug reported by Marcelo Castellanos */
			HBeta *= dense.eden * dense.xIonDense[ipHYDROGEN][1];

			/* CoolHeavy.xN2_A3_tot is fraction of excitations that produce a photon
			 * and represents the correction for collisional deexcitation */
			/*>>chng 05 dec 16, Liu et al. (2000) eqn 1 uses t = Te/10^4 K, not Te so phycon.te30
			 * is too large: (10^4)^0.3 = 16 - div by 15.8489 - bug caught by Kevin Blagrave */
			rec = LineSave.ipass <= 0 ? 0.0 :
				emit_frac(lineNII5755)* // 5754.61 *
				HBeta * 3.19 * phycon.te30 / 15.84893 * 
				dense.xIonDense[ipNITROGEN][2]/dense.xIonDense[ipHYDROGEN][1];
		}
		else
		{
			rec = 0.;
		}

		lindst( rec ,wlAirVac(5755.),"N 2R",ipnt,'t',true, "recombination" );
	}

	/***********************/
	/****** N 2 1085 ******/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipNITROGEN][1] )
	{
		static TransitionProxy lineNII1085;
		getTransition(LineID("N  2",1085.70), lineNII1085);

		double rec = LineSave.ipass <= 0 ? 0.0 :
			GetLineRec(201,1085)*emit_frac(lineNII1085); // 1085.70
		/* Collisional suppression from emit_frac_db of 1085A may not be accurate.
		 * It is based on the strongest line in the blend.*/
		PntForLine( 1085., "N 2R", &ipnt );
		lindst( MAX2(0.,rec), 1085., "N 2R", ipnt, 't', false, "dielectronic recombination" );
	}

	/***********************/
	/****** N 3 990 *******/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipNITROGEN][2] )
	{
		static TransitionProxy lineNIIIb;
		getTransition(LineID("N  3",989.799), lineNIIIb);

		double rec = LineSave.ipass <= 0 ? 0.0 :
			GetLineRec(216,991)*emit_frac(lineNIIIb); // 989.799

		PntForLine( 990., "N 3R", &ipnt );
		lindst(rec, 990., "N 3R", ipnt, 't', false, "recombination" );
	}

	/***********************/
	/****** N 4 765* ******/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipNITROGEN][3] )
	{
		static TransitionProxy lineNIV765;
		getTransition(LineID("N  4",765.147), lineNIV765);

		/* >>chng 97 may 02, better expression for dielectronic recombination */
		/* >>chng 02 jul 01, add function to get emission fraction */
		double rec = LineSave.ipass <= 0 ? 0.0 :
			GetLineRec(287,765)*emit_frac(lineNIV765); // 765.147

		/* dielectronic recombination contribution from Nussbaumer and Storey 1984 */
		PntForLine( 765., "N 4R", &ipnt );
		lindst( MAX2(0.,rec), 765., "N 4R", ipnt, 't', false,
			"recombination" );
	}
	/*** End of Nitrogen ***/

	/*** OXYGEN***/
	/* Recombination and Pump lines */

	/***********************/
	/****** O 2 Setup *****/
	/*********************/

	double rec7323 , rec7332, rec3730 , rec3726 , rec2471, reco23tot , reco22tot;

	static const bool debug2471 = false;
	static const bool debug7323 = false;
	static const bool debug7332 = false;

	if( atmdat.lgdBaseSourceExists[ipOXYGEN][1] )
	{
		/***********************/
		/****** O 2 2471  *****/
		/*********************/
		static TransitionProxy lineOII2471a;
		getTransition(LineID("O  2",2470.22), lineOII2471a);
		static TransitionProxy lineOII2471b;
		getTransition(LineID("O  2",2470.34), lineOII2471b);
		/**********************/

		/***********************/
		/****** O 2 3726  *****/
		/*********************/
		static TransitionProxy lineOII3726;
		getTransition(LineID("O  2",3726.03), lineOII3726);
		/*********************/

		/***********************/
		/****** O 2 3728  *****/
		/*********************/
		static TransitionProxy lineOII3728;
		getTransition(LineID("O  2",3728.81), lineOII3728);
		/********************/

		/***********************/
		/****** O 2 7323  *****/
		/*********************/
		static TransitionProxy lineOII7323a;
		getTransition(LineID("O  2",7318.92), lineOII7323a);
		static TransitionProxy lineOII7323b;
		getTransition(LineID("O  2",7319.99), lineOII7323b);
		/*********************/

		/***********************/
		/****** O 2 7332  *****/
		/*********************/
		static TransitionProxy lineOII7332a;
		getTransition(LineID("O  2",7329.67), lineOII7332a);
		static TransitionProxy lineOII7332b;
		getTransition(LineID("O  2",7330.73), lineOII7332b);
		/*********************/

		/* total recombination to 2P^o, the highest two levels of the 5-level atom,
		 * which produces the 7325 multiplet, last factor accounts for coll deexcitation
		 * this implements equation 2 of
		 * refer	o2	rec	Liu, X-W., Storey, P.J., Barlow, M.J., Danziger, I.J.,
		 * refercon	Cohen, M., & Bryce, M., 2000, MNRAS, 312, 585 */
		/* >>chng 05 dec 29, from first eqn, or unknown origin, to second, from indicated
		 * reference.  They agreed within 20% */

		/* Get lines in ergs/s to replace CoolHeavy.Oxxxx */
		double chO2471 = 0.0;
		double chO7323 = 0.0; 
		double chO7332 = 0.0; 
		double chO3726 = 0.0;
		double chO3730 = 0.0; 
		double O2_Lev45_rad = 0.0;
		double O2_Lev45_coll = 0.0;
		double O2_Lev23_rad = 0.0;
		double O2_Lev23_coll = 0.0;
		if (LineSave.ipass > 0 )
		{
			TransitionProxy tr;
			tr = lineOII2471a;
			chO2471 += tr.Emis().xObsIntensity();

			tr = lineOII2471b;
			chO2471 += tr.Emis().xObsIntensity();

			tr = lineOII7323a;
			chO7323 += tr.Emis().xObsIntensity();

			tr = lineOII7323b;
			chO7323 += tr.Emis().xObsIntensity();

			tr = lineOII7332a;
			chO7332 += tr.Emis().xObsIntensity();

			tr = lineOII7332b;
			chO7332 += tr.Emis().xObsIntensity();

			tr = lineOII3726;
			chO3726 += tr.Emis().xObsIntensity();

			tr = lineOII3728;
			chO3730 += tr.Emis().xObsIntensity();

			//Radiative decays per second of 7319.99 and 7329.67
			//Collisional decays per second of 7319.99 and 7329.67
			tr = lineOII7323b;
			O2_Lev45_rad += tr.Hi()->Pop()*tr.Emis().Aul();
			O2_Lev45_coll += tr.Hi()->Pop()*tr.Coll().ColUL(colliders);

			tr = lineOII7332a;
			O2_Lev45_rad += tr.Hi()->Pop()*tr.Emis().Aul();
			O2_Lev45_coll += tr.Hi()->Pop()*tr.Coll().ColUL(colliders);

			//Radiative decays per second of 3726.03 and 3728.81
			//Collisional decays per second of 3726.03 and 3728.81
			tr = lineOII3726;
			O2_Lev23_rad += tr.Hi()->Pop()*tr.Emis().Aul();
			O2_Lev23_coll += tr.Hi()->Pop()*tr.Coll().ColUL(colliders);

			tr = lineOII3728;
			O2_Lev23_rad += tr.Hi()->Pop()*tr.Emis().Aul();
			O2_Lev23_coll += tr.Hi()->Pop()*tr.Coll().ColUL(colliders);
		}

		//This replaces CoolHeavy.O2_A3_tot
		double O2_Lev45_radTot = O2_Lev45_rad/SDIV(O2_Lev45_rad + O2_Lev45_coll);

		//This replaces CoolHeavy.O2_A2_tot
		//FIX THIS: both rad and coll should have 3 lines
		double O2_Lev23_radTot = O2_Lev23_rad/SDIV(O2_Lev23_rad + O2_Lev23_coll);

		/* this test is necessary because for sims with no H-ionizing radiation
		 * the H+ density is initially zero */
		if( dense.xIonDense[ipHYDROGEN][1] > SMALLFLOAT )
		{
			double HBeta = (exp10(-20.89 - 0.10612*POW2(phycon.alogte - 4.4)))/phycon.te;
			HBeta *= dense.eden * dense.xIonDense[ipHYDROGEN][1];
			reco23tot = O2_Lev45_radTot * HBeta *
				9.36 * phycon.te40*phycon.te04 / 57.544 * dense.xIonDense[ipOXYGEN][2]/dense.xIonDense[ipHYDROGEN][1];
		}
		else
		{
			reco23tot = 0.;
		}

		if( debug2471 || debug7323 || debug7332 )
		{
			fprintf(ioQQQ,"reco23tot\t%e\n",reco23tot);
			fprintf(ioQQQ,"O2_Lev45_radTot\t%e\t%e\t%e\n",O2_Lev45_radTot,O2_Lev45_rad,O2_Lev45_coll);
		}

		sum = chO2471*2471./7325. + chO7323 + chO7332;
		if( sum > SMALLFLOAT )
		{
			/* assume effective branching ratio according to predicted intensities from 5-lev atom*/
			reco23tot /= sum;
		}
		else
		{
			reco23tot = 0.;
		}
		/* these are now ergs per sec unit vol for each transition */

		if( debug2471 || debug7323 || debug7332 )
		{
			fprintf(ioQQQ,"reco23tot\t%e\tsum\t%e\n",reco23tot,sum);
		}

		rec7332 = reco23tot * chO7332;
		linadd(rec7332,wlAirVac(7332.),"O 2R",'i',"recombination, P1/2-D3/2 and P3/2-D3/2 together" );

		rec7323 = reco23tot * chO7323;
		linadd(rec7323,wlAirVac(7323.),"O 2R",'i',"recombination, P1/2-D5/2 and P3/2-D5/2 together" );

		/* total recombination to 2D^o, the middle two levels of the 5-level atom,
		 * which produces the 3727 multiplet, last factor accounts for coll deexcit */
		reco22tot = 1.660e-10 / ( phycon.sqrte * phycon.te03 * phycon.te005 ) *
			dense.eden * dense.xIonDense[ipOXYGEN][2] * O2_Lev23_radTot;
		/* assume effective branching ratio according to predicted intensities from 5-lev atom*/
		sum = chO3726 + chO3730;
		// It would be better to have a physically-based branching ratio as fallback...
		realnum fracO3726 = (realnum) safe_div(chO3726, sum, 0.5),
			fracO3730 = (realnum) safe_div(chO3730, sum, 0.5);
		/* these are now ergs per sec unit vol for each transition */
		rec3726 = reco22tot * fracO3726 * 5.34e-12;
		rec3730 = reco22tot * fracO3730 * 5.34e-12;

		/***********************/
		/****** O 2 3727  *****/
		/*********************/

		/* O II 3727 produced by photoionization of O0 */
		oxy.s3727 = (realnum)((oxy.s3727 + oxy.s7325*0.5)*5.34e-12*
		  9.7e-5/(9.7e-5 + dense.eden*1.15e-6/phycon.sqrte));

		linadd(oxy.s3727,3727,"O 2H",'i',"line produced by photoionization of O0" );

		linadd( rec3726 ,3726,"O 2R",'i',"recombination, D3/2 - S3/2 transition" );

		PntForLine( 3729., "O 2R", &ipnt );
		// refer	o2	rec	Liu, X-W., Storey, P.J., Barlow, M.J., Danziger, I.J.,refercon	Cohen, M., & Bryce, M., 2000, MNRAS, 312, 585
		lindst( rec3730 ,3729., "O 2R", ipnt, 't', false,
				"recombination, Liu et al. 2000, MNRAS, 312, 585, D5/2 - S3/2 transition" );

		/***********************/
		/****** O 2 2471 Cont**/
		/*********************/

		rec2471 = reco23tot * chO2471*2471./7325. * 8.05e-12/2.72e-12;
		linadd(rec2471,wlAirVac(2471.),"O 2R",'i', "recombination, both 2P 1/2 and 3/2 to ground" );

		/***********************/
		/****** O 2 7325  *****/
		/*********************/

		oxy.s7325 = (realnum)(oxy.s7325*2.72e-12*0.34/(0.34 + dense.eden*
			  6.04e-6/phycon.sqrte));

		linadd(oxy.s7325,7325,"O 2H",'i',"line produced by photoionization of O0" );

		/***********************/
		/****** O 2 P and R ***/
		/*********************/

		/* the OII multiplets,
		 * contributions by both continuum pumping through XUV line
		 * and recombination */
		/* this is the driving line, pump is photons cm^-3 s^-1 */
		if( nWindLine > 0 )
		{
			pump = TauLine2[387].Emis().pump()*TauLine2[387].Emis().PopOpc();
		}
		else
		{
			pump = 0.;
		}

		PntForLine(3120.,"O  2",&ipnt);
		lindst(pump*0.336 * 6.37e-12/(1.+dense.eden/1e12) ,3120,"O 2P",ipnt,'r',true,
			"OII 3113.62 - 3139.68 (8 lines) are only pumped, no recombination part" );

		PntForLine(3300.,"O  2",&ipnt);
		lindst(pump*0.147 * 6.03e-12/(1.+dense.eden/1e12) ,3300,"O 2P",ipnt,'r',true,
			"OII 3277.56 - 3306.45 (6 lines) are only pumped, no recombination part" );

		PntForLine(3762.,"O  2",&ipnt);
		lindst(pump*0.087 * 5.29e-12/(1.+dense.eden/1e12) ,3762,"O 2P",ipnt,'r',true,
			"OII 3739.76/3762.47/3777.42 (3 lines) are only pumped, no recombination part" );

		/* recombination and specific pump for OII 4638.86-4696.35 (8 lines) */
		rec = GetLineRec(82, 4651 );
		PntForLine( 4651.,"O  2",&ipnt);
		lindst(rec, 4651.,"O 2R",ipnt,'t',true,
			"total recombination, 4638.86-4696.35 (8 lines)" );

		/* convert UV pump rate to intensity with branching ratio and hnu recombination
		 * part of O II 4651 line */
		linadd(pump* 0.336 * 0.933 * 4.27e-12/(1.+dense.eden/1e12),4651,"O 2P",'i',
			"excitation by pumping" );

		/***********************/
		/****** O 2 4341  *****/
		/*********************/

		/* recombination and specific pump for OII 4317.14-4366.89 (6 lines) */
		rec = GetLineRec(83, 4341 );

		PntForLine( wlAirVac(4341.), "O 2R", &ipnt );
		lindst(rec/(1.+dense.eden/1e12), wlAirVac(4341.), "O 2R", ipnt, 't', false,
			"recombination" );

		linadd(pump* 0.147 * 0.661 * 4.58e-12/(1.+dense.eden/1e12),wlAirVac(4341.),"O 2P",'i',
			"excitation by pumping" );

		/***********************/
		/****** O 2 3736  *****/
		/*********************/
		/* recombination and specific pump for OII 3712.74/3727.32/3749.48 (3 lines) */
		double rec = GetLineRec(84, 3736 );
		/* convert UV pump rate to intensity with branching ratio and hnu */

		PntForLine( wlAirVac(3736.), "O 2R", &ipnt );
		lindst(rec/(1.+dense.eden/1e12), wlAirVac(3736.), "O 2R", ipnt, 't', false,
			"recombination" );
		linadd(pump* 0.087 * 0.763 * 5.33e-12/(1.+dense.eden/1e12),wlAirVac(3736.),"O 2P",'i',
			   "excitation by pumping" );
	}

	/***********************/
	/****** O 3 1666  *****/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipOXYGEN][2] )
	{
		static TransitionProxy lineOIII1666a;
		getTransition(LineID("O  3",1666.15), lineOIII1666a);
		static TransitionProxy lineOIII1666b;
		getTransition(LineID("O  3",1660.81), lineOIII1666b);

		double efac = 0.0;
		if (LineSave.ipass > 0)
		{
			efac = (emit_frac(lineOIII1666a) + emit_frac(lineOIII1666b))*0.5;
		}
			
		linadd(ionbal.PhotoRate_Shell[ipOXYGEN][3][1][0]*dense.xIonDense[ipOXYGEN][1]*0.3*1.20e-11*efac,1665,"O 3H",'i',
			"due to inner shell (2s^2) ionization" );

		linadd(oxy.AugerO3*1.20e-11*efac*0.27,1665,"O 3A",'i',
				"due to K-shell ionization" );
	}

	/***********************/
	/****** O 3 5007  *****/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipOXYGEN][2] )
	{
		static TransitionProxy lineOIII5007;
		getTransition(LineID("O  3",5006.84), lineOIII5007);

		/* o iii 5007+4959, As 96 NIST */
		/*The cs of the transitions 3P0,1,2 to 1D2 are added together to give oiii_cs3P1D2 */
		/*the cs of the transition 1D2-1S0 is mentioned as oiii_cs1D21S0*/
		/*The cs of the transitions 3P0,1,2 to 1S0 are added together to give oiii_cs3P1S0*/
		double aeff = 0.027242 + oxy.d5007r;
		//double a21 = aeff - oxy.d5007r;

		realnum d5007t = 0.0;
		if (LineSave.ipass > 0)
		{
			TransitionProxy tr = lineOIII5007;
			d5007t += (realnum)(tr.Emis().xObsIntensity()*oxy.d5007r/aeff);
		}

		linadd(d5007t/1.25,5007,"LOST",'i',"O III 5007 lost through excited state photo");
	}

	/***********************/
	/****** O 3 4363  *****/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipOXYGEN][2] )
	{
		/* collisional quenching ratio */
		double effec = 1.6/(1.6 + 0.9*dense.cdsqte);

		/* O III 4363 recombination, coefficient from Burgess and Seaton */
		/*
		double r43630 = 6.3e-21/(phycon.te70*phycon.te10)*dense.eden*dense.xIonDense[ipOXYGEN][3]*effec;

		 * Hi Gary, Daniel, Francisco and Veronica,
			I was wondering why the Liu's fit to the recombination of OIII4363 was so different from what I obtained, from PPB91 or BS60. The point is that the dielectronic recombination is also included in Liu's formula. It is taken from Nussbaumer and Storey 1984. If I add this contribution to PPB91, I obtain the same result than Liu, for Te> 5,000K. At lower Te, the Liu's fit is not correct as only the recombination is acting and the slope changes.

			A full correct recombination computation of 4363 must then take these 2 effects into account, using something like:
			3.0e-22 * T**-0.6 + 6.37e-25 * exp(-8000/T)
			i.e. in C++:
			update:  published in  2020MNRAS.498L..82G
			*/
		double r4363 = (3.0e-22/(phycon.te40*phycon.te20) + 6.37e-25*exp(-8000/phycon.te) ) *
										   dense.eden*dense.xIonDense[ipOXYGEN][3]*effec;

			/*This is fitting PPB91+NS84 within 3% from 200 to 30,000K.
			I'm not sure that this is what you want to include in Cloudy?

			Best regards,
			Christophe Morisset
		 */


		/* charge exchange,
		 * >>refer	O3	CT	Dalgarno+Sternberg ApJ Let 257, L87.
		 * scaled to agree with
		 * >>refer	O3	CT	Gargaud et al AA 208, 251, (1989) */
		double ct4363 = phycon.sqrte*1.3e-12*4.561e-12*dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipOXYGEN][3]*
		  effec;

		PntForLine(4363.,"O  3",&ipnt);

		lindst(r4363,wlAirVac(4363.),"O 3R",ipnt,'t',true,
			"recombination, coefficient from Burgess and Seaton" );

		linadd(ct4363,wlAirVac(4363.),"O 3C",'i' ,
			"charge exchange, Dalgarno+Sternberg ApJ Let 257, L87");
	}

	/***********************/
	/****** O 3 5592  *****/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipOXYGEN][2] )
	{
		linadd(dense.xIonDense[ipHYDROGEN][0]*dense.xIonDense[ipOXYGEN][3]*0.225*3.56e-12*1.34e-11*phycon.sqrte,
				wlAirVac(5592.),"O 3C",'i',"charge exchange rate, D+S" );
	}

	/***********************/
	/****** O 3 835  ******/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipOXYGEN][2] )
	{
		static TransitionProxy lineOIII835;
		getTransition(LineID("O  3",835.059), lineOIII835);

		/* >>chng 97 may 02, better rec contribution */
		double rec = LineSave.ipass <= 0 ? 0.0 : GetLineRec(331,835)*emit_frac(lineOIII835);

		PntForLine( 835., "O 3R", &ipnt );
		lindst(MAX2(0.,rec), 835., "O 3R", ipnt, 't', false, "dielectronic recombination only" );
	}

	/***********************/
	/****** O 4 1402 ******/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipOXYGEN][3] )
	{
		linadd(ionbal.PhotoRate_Shell[ipOXYGEN][2][1][0]*dense.xIonDense[ipOXYGEN][2]*0.43*1.42e-11,1401,"O 4H",'i',
				"inner shell photoionization, relaxation" );
	}

	/***********************/
	/****** O 4 789  ******/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipOXYGEN][3] )
	{
		static TransitionProxy lineOIV789;
		getTransition(LineID("O  4",790.201), lineOIV789);

		/* >>chng 97 may 02, better rec contribution */
		double rec = LineSave.ipass <= 0 ? 0.0 : GetLineRec(378,789)*emit_frac(lineOIV789);
		PntForLine( 789., "O 4R", &ipnt );
		lindst(MAX2(0.,rec), 789., "O 4R", ipnt, 't', false, "dielectronic recombination only" );
	}

	/***********************/
	/****** O 5 630  ******/
	/*********************/
	if( atmdat.lgdBaseSourceExists[ipOXYGEN][4] )
	{
		static TransitionProxy lineOV630;
		getTransition(LineID("O  5",629.732), lineOV630);

		/* >>chng 97 may 02, better rec contribution */
		double rec = LineSave.ipass <= 0 ? 0.0 : GetLineRec(466,630)*emit_frac(lineOV630);

		PntForLine( 630., "O 5R", &ipnt );
		lindst(MAX2(0.,rec), 630., "O 5R", ipnt, 't', false,
			"dielectronic recombination only" );
	}

	/* add up line intensities for certain set of lines,
	 * must come after all lines defined */
	sum = PrtLineSum();
	/* zero out the location that will receive this information,
	 * remember that memory not allocated until ipass >= 0 */
	if( LineSave.ipass > 0 )
		LineSave.lines[LineSave.nsum].SumLineZero();
	/* optional sum of certain emission lines, set with "print sum" */
	linadd(sum/radius.dVeffAper,0,"Stoy",'i' ,
		"Stoy method energy sum");

	/* >>chng 06 jan 03, confirm that number of lines never changes once we have
	 * created the labels */
	{
		static long nLineSave=-1 , ndLineSave=-1;
		if( LineSave.ipass == 0 )
		{
			nLineSave = LineSave.nsum;
			ndLineSave = LineSave.nsum;
			LineSave.setSortWL();
		}
		else if( LineSave.ipass > 0 )
		{
			/* this can't happen */
			if( nLineSave<= 0 || ndLineSave < 0 )
				TotalInsanity();

			/* now make sure that we have the same number of lines as we had previously
			 * created labels.  This would not pass if we did not add in exactly the same
			 * number of lines on each pass */
			if( nLineSave != LineSave.nsum )
			{
				fprintf( ioQQQ, "DISASTER number of lines in LineSave.nsum changed between pass 0 and 1 - this is impossible\n" );
				fprintf( ioQQQ, "DISASTER LineSave.nsum is %li and nLineSave is %li\n",
					LineSave.nsum , 
					nLineSave);
				ShowMe();
				cdEXIT(EXIT_FAILURE);
			}
			if( ndLineSave != LineSave.nsum )
			{
				fprintf( ioQQQ, "DISASTER number of lines in LineSave.nsum changed between pass 0 and 1 - this is impossible\n" );
				fprintf( ioQQQ, "DISASTER LineSave.nsum is %li and ndLineSave is %li\n",
					LineSave.nsum , 
					ndLineSave);
				ShowMe();
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " lines returns\n" );
	}
	return;
}

STATIC void lines_iron_Ka()
{
	DEBUG_ENTRY( "cool_iron_Ka()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   cool_iron_Ka called\n" );
	}

	double FeKaHLike = 0. , FeKaHeLike=0.;
	/* one and two electron Ka */
	if( dense.lgElmtOn[ipIRON] )
	{
		/* H-like one-electron Ka */
		FeKaHLike = iso_sp[ipH_LIKE][ipIRON].trans(ipH2p,ipH1s).Emis().xIntensity();

		/* He-like two-electron Ka */
		FeKaHeLike =
				iso_sp[ipHE_LIKE][ipIRON].trans(ipHe2p1P,ipHe1s1S).Emis().xIntensity()+
				iso_sp[ipHE_LIKE][ipIRON].trans(ipHe2p3P0,ipHe1s1S).Emis().xIntensity()+
				iso_sp[ipHE_LIKE][ipIRON].trans(ipHe2p3P1,ipHe1s1S).Emis().xIntensity()+
				iso_sp[ipHE_LIKE][ipIRON].trans(ipHe2p3P2,ipHe1s1S).Emis().xIntensity();

		/* total intensity of K-alpha line, cold, grain, hot, 1 and two electron
		 * 19 sep 21 segfault if evaluated with dense.lgElmtOn[ipIRON] false
		 * since continuum index not defined */
		lindst((fe.fekcld+fe.fegrain)*1.03e-8+fe.fekhot*1.11e-8+FeKaHLike+FeKaHeLike,1.78f,"FeKa",
			iso_sp[ipH_LIKE][ipIRON].trans(ipH2p,ipH1s).ipCont(),'i',false,
			   "total intensity of Fe K-alpha line, grain, cold, hot, 1 and 2 electron" );
	}

	linadd(FeKaHLike,1.78177,"FeK1",'i' ,
		"H-like one-electron Ka");

	linadd(FeKaHeLike,1.85,"FeK2",'i' ,
		"He-like two-electron Ka");

	linadd( fe.fekhot*1.11e-8 ,1.8,"FeKH",'i' ,
		"fluorescent hot iron, Fe 18 - 23 times ionized");

	linadd(fe.fekcld*1.03e-8,1.75,"FeKC",'i',
		"fluorescent cold iron, less than or 17 times ionized" );

	linadd(fe.fegrain*1.03e-8,1.75,"FeKG",'i' ,
		"grain production of cold iron");

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   cool_iron_Ka returns\n" );
	}
	return;
}
