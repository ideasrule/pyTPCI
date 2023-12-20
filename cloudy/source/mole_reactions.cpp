/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include "cddefines.h"
#include "prt.h"
#include "version.h"
#include "hmi.h"
#include "conv.h"
#include "thermal.h"
#include "grainvar.h"
#include "ionbal.h"
#include "opacity.h"
#include "radius.h"
#include "atmdat.h"
#include "trace.h"
#include "deuterium.h"
#include "grains.h"
#include "mole_priv.h"
#include "gammas.h"
#include "mole.h"
#include "freebound.h"
#include "dense.h"
#include "ipoint.h"

/*
 * HOWTO:- add a reaction to the new CO network, as at 2006 December 18.
 *
 * add a line of the form
 *   O,C2=>CO,C:hmrate:SCALE:B:C # Data source
 * to the relevant data file for the new reaction.  The first component is 
 * the chemical reaction, the second is a string naming the function which 
 * is used to evaluate the rate coefficients.
 *
 * SCALE is an overall constant by which the reaction rate is scaled,
 * the remaining arguments constants used by this function.
 *
 * If all the species have previously been defined, and the parser in
 * mole_co_etc.c can understand the reaction string, then that's it
 * (unless the rate function is of a different form to those already
 * defined).
 *
 * The sources and sinks to other networks will need to be defined, 
 * as this hasn't been made generic yet.
 *
 */

/** set UDFA = 1 to make UDFA comparison and stop
 * http://udfa.ajmarkwick.net/index.php */
enum {UDFA=0};

STATIC void newreact(const char label[], const char fun[], double a, double b, double c);
STATIC long parse_reaction( shared_ptr<mole_reaction>& rate, const char label[] );
STATIC string canonicalize_reaction_label( const char label[] );
STATIC void canonicalize_reaction( shared_ptr<mole_reaction>& rate );
STATIC bool lgReactionTrivial( shared_ptr<mole_reaction>& rate );
STATIC void register_reaction_vectors( shared_ptr<mole_reaction> rate );
/* run through all reactions and report which ones do not have a reverse reaction */
STATIC void mole_check_reverse_reactions(void);
STATIC double mole_get_equilibrium_condition( const char buf[] );
STATIC double mole_get_equilibrium_condition( const mole_reaction* const rate );
STATIC double mole_partition_function( const molecule* const sp);
STATIC void mole_generate_isotopologue_reactions( string atom_old, string atom_new );

STATIC double sticking_probability_H_func( double T_gas, double T_grain );
// Calculate sticking probability of H on grains according to Hollenback and McKee 1979.
STATIC double sticking_probability_H_HM79( double T_gas, double T_grain );

/* Save PBM of network coupling matrix fill */
STATIC void	plot_sparsity(void);

#ifndef NDEBUG
/* Check invariants for defined reaction (species and charge totals) */
STATIC bool lgReactBalance(const shared_ptr<mole_reaction> &rate);
#endif

/* Functions to read UDFA database and compare results with Cloudy's internal reaction set */
STATIC void read_data(const char file[], void (*parse)(char *s));
STATIC void parse_base(char *s);
STATIC void parse_udfa(char *s);
STATIC void compare_udfa(const shared_ptr<mole_reaction> &rate);

/* Nick Abel, 06 Nov 27 states: 
	 "In all the crnu reactions, the factor of two is the albedo factor 1/(1-albedo)."
	 Can this be obtained from the grain physics? */
static realnum albedo = 0.5;


namespace {
/* Define ratefunc type to make specification of newfunc and findfunc clearer */
	typedef shared_ptr<mole_reaction> ratefunc;
	template<class T> void newfunc()
	{
		ratefunc fun = shared_ptr<mole_reaction>(new T);
		ASSERT( mole_priv::functab.find(fun->name()) == mole_priv::functab.end() );
		mole_priv::functab[fun->name()] = fun;
	}
	ratefunc findfunc(const char name[]);
	
/* The functions which define the reaction rates -- must all have same template */
	class mole_reaction_hmrate;
	double hmrate(const mole_reaction *rate);
	class mole_reaction_th85rate;
	double th85rate(const mole_reaction *rate);
	class mole_reaction_crnurate_noalbedo;
	double crnurate_noalbedo(const mole_reaction *rate);
/*  >>chng 07 Dec 11, add photodesorption of molecules frozen on grain surfaces */
/*  >>chng 07 Dec 11, add grain surface reactions to chemistry.  Now two molecules adsorbed on a grain can interact to form a new molecule */
	class mole_reaction_grn_react;
	double grn_react(const mole_reaction *rate);
	class mole_reaction_h2_collh_excit;
	double h2_collh_excit(const mole_reaction *rate);
	class mole_reaction_h2_collh2_excit;
	double h2_collh2_excit(const mole_reaction *rate);
	class mole_reaction_h2_collh_deexc;
	double h2_collh_deexc(const mole_reaction *rate);
	class mole_reaction_h2_collh2_deexc;
	double h2_collh2_deexc(const mole_reaction *rate);
	class mole_reaction_rh2g_dis_h;
	double rh2g_dis_h(const mole_reaction *rate);
	class mole_reaction_rh2s_dis_h;
	double rh2s_dis_h(const mole_reaction *rate);
	class mole_reaction_rh2g_dis_h2;
	double rh2g_dis_h2(const mole_reaction *rate);
	class mole_reaction_rh2s_dis_h2;
	double rh2s_dis_h2(const mole_reaction *rate);
	class mole_reaction_rh2s_dis_h2_nodeexcit;
	double rh2s_dis_h2_nodeexcit(const mole_reaction *rate);
	class mole_reaction_bh2g_dis_h;
	double bh2g_dis_h(const mole_reaction *rate);
	class mole_reaction_bh2s_dis_h;
	double bh2s_dis_h(const mole_reaction *rate);
	class mole_reaction_bh2g_dis_h2;
	double bh2g_dis_h2(const mole_reaction *rate);
	class mole_reaction_hneut;
	double hneut(const mole_reaction *rate);
	class mole_reaction_cionhm;
	double cionhm(const mole_reaction *rate);
}

/*
 * Functions to specify chemical rates -- note that the rate->a overall scale
 * factor is applied in mole_update_rks
 *
 */

#include "phycon.h"
#include "doppvel.h"

namespace
{
	class mole_reaction_null : public mole_reaction
	{
		typedef mole_reaction_null T;  
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "null";}
		double rk() const
			{
				ASSERT( false );
				return 0.0;
			}
	};
}


/* Add in non-equilibrium chemistry.  This is done by assuming
 * that turbulence reduces the temperature barrier, thereby
 * enhancing reaction rates for molecules such as CH+.  The 
 * "effective temperature is defined according to 
 * >>refer Federman et al. 1996, MNRAS. L41-46  */

namespace {
/* The effective temperature is defined as:
 * T(effective) = T + (1/3)*reduced_mass*turbulent_velocity^2/BOLTZMANN_CONSTANT
 */
	STATIC double noneq_offset(const mole_reaction *rate)
	{	
		/* This logic could be cached by using distinct rate functions in newreact */
		int nreact, n;
		bool lgFact;
		
		DEBUG_ENTRY( "noneq_offset()" );
		
		lgFact = false;
		if(mole_global.lgNonEquilChem)
		{ 
			if(mole_global.lgNeutrals) 
			{
				lgFact = true;
			}
			else
			{
				nreact = rate->nreactants;
				for(n=0;n<nreact;n++) 
				{
					if (rate->reactants[n]->charge != 0)
					{
						lgFact = true;
						break;
					}
				}
			}
		}
		
		if( lgFact ) 
			return 0.333f*POW2(DoppVel.TurbVel)/BOLTZMANN*rate->reduced_mass;
	else
		return 0.;
	}
	/*hmrate - evalurate UMIST expression for rate coefficient */
	double hmrate(const mole_reaction *rate) 
	{
		double te;
		
		DEBUG_ENTRY( "hmrate()" );
		/* the UMIST equation is
		 * rate = alpha \times (T/300)^beta \times exp( -gamma/T 
		 alpha==rate->a; beta== rate->b; gamma==rate->c */  
		
		te = phycon.te+noneq_offset(rate);
		/* UMIST rates are simple temperature power laws that
	 	 * can become large at the high temperatures Cloudy
	 	 * may encounter. Do not extrapolate to above T>2.5e3K */ 
		/* rate-b is the power beta in (T/300)^beta, positive beta
		 * can diverge at high temperatures */
		/* THIS CODE MUST BE KEPT PARALLEL WITH HMRATE4 IN MOLE.H */ 
		if( rate->b > 0.)	
			te = min(te, 2500.);
		/* rate->c is gamma in expontntial */
		if( rate->c < 0. )
			ASSERT( -rate->c/te < 10. );

		double r = 1.;
		if( rate->b != 0. )
			r *= pow(te/300.,rate->b);
		if( rate->c != 0. )
			r *= exp(-rate->c/te);
		return r;
	}
	
	class mole_reaction_hmrate_exo : public mole_reaction
	{
		typedef mole_reaction_hmrate_exo T;
	public:
		virtual T* Create() const {return new T;}
		
		virtual const char* name() {return "hmrate_exo";}
		
		double rk() const
			{
				double te;
				
				DEBUG_ENTRY( "hmrate_exo()" );
				
				te = phycon.te+noneq_offset(this);
				
				if( this->c < 0. )
					te = MIN2( te, -10. * this->c );

				return pow(te/300.,this->b)*exp(-te/this->c);
			}
	};
	
	class mole_reaction_hmrate : public mole_reaction
	{
		typedef mole_reaction_hmrate T;  
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "hmrate";}
		double rk() const
			{
				return hmrate(this);
			}
	};
	
	class mole_reaction_constrate : public mole_reaction
	{
		typedef mole_reaction_constrate T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "constrate";}
		double rk() const
			{
				return 1.;
			}
	};
	
	/* hmi.UV_Cont_rel2_Habing_TH85_depth is field relative to Habing background, dimensionless */
	/* >>chng 04 apr 01, move from TH85 to DB96, also correct leading coef to use
	 * UMIST database value */
	/* CO_photo_dissoc_rate = 1.67e-10f*hmi.UV_Cont_rel2_Habing_TH85_depth;*/
	
	/* TRY MOVING PHOTORATES OUT OF LOOP */

	/* >>chng 02 jul 04 -- The following are dissociation rates for various molecular species
	For right now we will calculate this rate by the standard form:
		(alpha)*Go*exp(-Beta*AV)
	when the command "set Leiden hack UMIST rates" is used.  Otherwise we
	will just let cloudy calculate the value of the UV radiation field */
}
#include "rfield.h"
namespace {
	double th85rate(const mole_reaction *rate) 
	{
		double rk;

		DEBUG_ENTRY( "th85rate()" );
		
		if (!mole_global.lgLeidenHack || rate->c == 0.0)
		{
			rk = hmi.UV_Cont_rel2_Habing_TH85_depth/1.66;
		}
		else 
		{
			rk = hmi.UV_Cont_rel2_Habing_TH85_face/1.66*exp(-(rate->c*rfield.extin_mag_V_point));
		}
		
		return rk;
	}
	class mole_reaction_th85rate : public mole_reaction
	{
		typedef mole_reaction_th85rate T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "th85rate";}
		double rk() const
			{
				return th85rate(this);
			}
	};
}	

#include "secondaries.h"
namespace {
	/* >> chng aug 24, 05 NPA This is the cosmic ray ionization rate used in the molecular network.  
	 * TH85 and the e-mail from Amiel Sternberg has each cosmic ray ionization rate as a 
	 * leading coefficient multiplied by a scale factor.
	 * The leading coefficient is defined as the cosmic ray rate for H + CRPHOT = > H+
	 * + e- .  For molecules in the heavy element molecular network, this scale
	 * factor is derived by taking the rate for:

	X + CRPHOT => Y + Z 
	and dividing it by the rate:
	H + CRPHOTP => H+ + e-

	This scale factor is 2.17 for all cosmic ray reactions in this network 
		crnu_rate = secondaries.csupra[ipHYDROGEN][0];*/
	double crnurate_noalbedo(const mole_reaction *) 
	{
		return 2.17*secondaries.csupra[ipHYDROGEN][0];
	}
	class mole_reaction_crnurate_noalbedo : public mole_reaction
	{
		typedef mole_reaction_crnurate_noalbedo T;
	public:
		virtual T* Create() const {return new T;}
		
		virtual const char* name() {return "crnurate_noalbedo";}
		double rk() const
			{
				return crnurate_noalbedo(this);
			}
	};
	
	class mole_reaction_crnurate : public mole_reaction
	{
		typedef mole_reaction_crnurate T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "crnurate";}
		double rk() const
			{
				return crnurate_noalbedo(this)/(1.0-albedo);
			}
	};
}
#include "hextra.h"
namespace {
/* Can this be found directly from radiation field?? */
	class mole_reaction_cryden_ov_bg : public mole_reaction
	{
		typedef mole_reaction_cryden_ov_bg T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "cryden_ov_bg";}
		double rk() const
			{
				return hextra.cryden_ov_background;
			}
	};
	
	class mole_reaction_co_lnu_c_o_lnu : public mole_reaction
	{
		typedef mole_reaction_co_lnu_c_o_lnu T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "co_lnu_c_o_lnu";}
		double rk() const
			{
				double val = 0;
				int ns, ion;
				/* inner shell photoionization of CO, assume rates are same as K 1s and 2s
				 * shell of C and O */
				/* >>chng 04 may 26, upper limit should be ns<2, had been ns<2 so picked up
			 * valence shell, which was incorrect */
				
				DEBUG_ENTRY( "co_lnu_c_o_lnu()" );
				
				for( ns=0; ns<2; ++ns )
				{
					ion = 0;
					val += ionbal.PhotoRate_Shell[ipCARBON][ion][ns][0];
					val += ionbal.PhotoRate_Shell[ipOXYGEN][ion][ns][0];
				}
			
				return val;
			}
	};
	
	/******************************** Gas-Grain Chemistry**********************************/

	/*  The Gas-Grain Chemistry rates are taken from: 
	>>refer	Hasegawa, T. I. & Herbst, E. 1993, MNRAS, 261, 83

	So far only CO depletion onto grains is considered, however, this code can be generalized 
	if desired to other molecules, using the data in the above paper.  There are three important reactions 
	to determine the abundance of a molecule on grain surfaces deep in molecular clouds.  The
	rate of accretion of molecule X onto grains is

	R(accretion) = PI*[grain_radius]^2*[thermal velocity]*[density_of_grains]

	Two processes remove molecule X from the grain surface, and that is thermal evaporation, due
	to the heat of the grain, and cosmic ray deabsorption.  The first of these rates come from the 
	above paper, and depends primarily on the dust temperature.  The cosmic ray rate is a constant,
	calculated in Hasegawa and Herbst.  

	For each molecule desired, I have created a new species which is the density of that molecule
	on the grain surface */


	/* evaporation rate of molecule on grain is:

	k(evap) = [vibrational absorption frequency]*exp[-binding_energy/dust_temperature]

	The binding energies come from Hasegawa and Herbst, Table 4.  The vibrational frequency comes from
	equation 3 of >>refer Hasegawa, T. I., Herbst, E., & Leung, C. M. 1992, ApJSS, 82, 167

	[vibrational absorption frequency] = 
	SQRT[[2*number_of_sites_for_grain_absorption*binding_energy]/[PI^2*mass_of_molecule]]

	**********************************************************************************************/

	class mole_reaction_vib_evap : public mole_reaction
	{
		typedef mole_reaction_vib_evap T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "vib_evap";}
		double rk() const
			{
				double binding_energy,  exponent, vib_freq;
			
				DEBUG_ENTRY( "vib_evap()" );
				
				exponent = 0.0;
				
				binding_energy = this->b;
				double bin_total=0.0;
				for( size_t nd=0; nd < gv.bin.size() ; nd++ )
				{
					double bin_area = gv.bin[nd].IntArea*gv.bin[nd].cnv_H_pCM3;
					exponent += exp(-binding_energy/gv.bin[nd].tedust)*bin_area;
					bin_total += bin_area;
				}
				exponent /= bin_total;
				const double surface_density_of_sites = 1.5e15;
				
				/* >> chng 07 dec 08 rjrw -- add 0.3*BOLTZMANN factor from Hasegawa & Herbst */
				vib_freq = sqrt(2*surface_density_of_sites*(0.3*BOLTZMANN)*binding_energy/(PI*PI*this->reactants[0]->mole_mass));
				//fprintf(ioQQQ,"Vib freq %g for mass %g\n",vib_freq,rate->reactants[0]->mole_mass);
				
				/*>>chng 06 jan 11, NPA - In some H+ regions the grain temperatures are so low
				  that molecular freeze out occurs.  This should not happen, because the ices
				  should sublimate in such a harsh environment.  Therefore, we introduce an
				  artificial sublimation rate to destroy ices.  THIS IS NOT A PHYSICAL RATE!!!!
				  only a rate that gets the desired, realistic result */
				/*>>chng 06 sep 03 rjrw -- include this in standard evaporation rate coeff (the artificial part is the sexp term) */
				/** \todo	0	find physical theory for this process */
				/* Rate comes from Table curve and assumes that rate is high (~1) in H+
				 * region and negligible ( < 1e-30) in molecular cloud - designed to stop
				 * freeze-out above 100 K */
				
				return vib_freq*exponent +sexp( 555.89/phycon.sqrte - 5.55 );
		}
	};
	
	class mole_reaction_grn_abs : public mole_reaction
	{
		typedef mole_reaction_grn_abs T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "grn_abs";}
		double rk() const
			{
				double mass;
				
				DEBUG_ENTRY( "grn_abs()" );
				
				/* Can't rely on order, as it is standardized using mole_cmp */
				if (this->reactants[0]->n_nuclei() != 0)
					mass = this->reactants[0]->mole_mass;
				else
					mass = this->reactants[1]->mole_mass;
				
				return sqrt(8.*BOLTZMANN*phycon.te/(PI*mass));
			}
	};
	
	class mole_reaction_grn_react : public mole_reaction
	{
		typedef mole_reaction_grn_react T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "grn_react";}
		double rk() const
			{
				return grn_react(this);
			}
	};
	double grn_react(const mole_reaction *rate)
	{
		/*  This function computes the formation of a new molecule on the surface of 
		 *  of a grain due to the interaction between two species on the grain surface.
		 *  We already do this type of reaction for H + H --> H2, but now we are doing 
 	 *  it for other species as well. The rate depends on:
	 
	 1) The ability of each species to move on the surface
	 2) The ability of each species to interact and form a new molecule
	 3) The ability of each species to not evaporate away before interacting
	 
	 * The variable "number_of_sites" tells us how many sites there are which 
	 * a molecule can attach itself to a grain.  quant_barrier tells us the size
 	 * of the quantum mechanical barrier between the two atoms.  E_i and E_j 
	 * relates to how well a species can stay on a grain of temperature T(dust).  Finally,
	 * activ_barrier provides information on how well two atoms on the surface of a grain 
	 * can actually interact to form a new molecule. The formalism below is taken from
	 * the Hasegawa, Herbst, and Leung paper from 1992 (referenced above)*/
		
		DEBUG_ENTRY( "grn_react()" );
		
		double quant_barrier = 1e-8; /* 1 Angstrom rectangular barrier */
		double surface_density_of_sites = 1.5e15; /* number of sites available for adsorption */
		fixit("rate->a is _always_ overall scaling factor");
		ASSERT( rate->nreactants == 2 ); // this routine seems coded to only work with 2 reactants
		double E_i = rate->reactants[0]->form_enthalpy; /* adsorption energy barrier (w/ grain) for first reactant */
		double E_j = rate->reactants[1]->form_enthalpy; /* adsorption energy barrier (w/ grain)for second reactant */
		double activ_barrier = rate->c; /* energy barrier between atoms on the grain */
		
		/* Each species has a vibrational frequency, dependent on its adsorption energy barrier */
		fixit("Ordering of reactants may have switched");
		double vib_freq_i = sqrt(2*surface_density_of_sites*(0.3*BOLTZMANN)*E_i/(PI*PI*rate->reactants[0]->mole_mass));
		double vib_freq_j = sqrt(2*surface_density_of_sites*(0.3*BOLTZMANN)*E_j/(PI*PI*rate->reactants[1]->mole_mass));
		
		double Exp_i = 0.0;
		double Exp_j = 0.0;
		double dust_density = 0.0;
		
		/* Compute thermal evaporation terms along with total dust number density */
		fixit("should there be an area weighting to exponential terms?");
		
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			double bin_density = gv.bin[nd].IntArea*gv.bin[nd].cnv_H_pCM3;
			Exp_i += exp(-E_i/gv.bin[nd].tedust)*bin_density;
			Exp_j += exp(-E_j/gv.bin[nd].tedust)*bin_density;
			dust_density += bin_density/(4*1e-10);
		}
		
		ASSERT(fp_equal((realnum)dust_density,(realnum)(mole.grain_area/1e-10)));
		
		double total_sites = 4.*mole.grain_area*surface_density_of_sites;
		
		/* rates of migration on surface for each species */
		double diff_i = vib_freq_i*Exp_i/total_sites;
		double diff_j = vib_freq_j*Exp_j/total_sites;
		
		/* Exponential term models the likelyhood that two species, once they meet on the surface, will interact.  Larger the barrier, the smaller
		 * the chance of a new molecule */
		double scale = exp(-2*(quant_barrier/(HBAReV*EN1EV))*sqrt(2.0*rate->reduced_mass*0.3*BOLTZMANN*activ_barrier));
		
		/* Total Rate */
		return scale*(diff_i + diff_j)/SDIV(dust_density);
	}
	
	class mole_reaction_grn_photo : public mole_reaction
	{
		typedef mole_reaction_grn_photo T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "grn_photo";}
		double rk() const
			{
				/* This function models the rate at which UV photons remove
				 * molecules adsorbed on the surface of grains, using the treatment given in:
				 * >>refer Bergin, E. A., Langer, W. D., & Goldsmith, P. F. 1995, ApJ, 441, 222 
				 * The following treatment uses their equation 3 */
				
				DEBUG_ENTRY( "grn_photo()" );
			
				/* This is the number of molecules removed from the grain per
				 * incident photon use same value of 1e-4 as used in Bergin, include
				 * conversion to radiation field units to get the dimensions right
				 * (cf the d'Hendercourt reference in Bergin et al) */
				
				fixit("grainn photoionization"); 
				// 2e-15 is mean adsorbed species cross section, from
				// >>refer Greenberg, L. T., 1973, IAUS, 52, 413
				// 1.232e7f * 1.71f is from DB96 referred to below 
				// will be multiplied by yield factor (~1e-4) to get overall photodesorption rate
				return 2e-15 * hmi.UV_Cont_rel2_Draine_DB96_depth *(1.232e7f * 1.71f);
			}
	};
}	
#include "rt_escprob.h"
namespace {
	class mole_reaction_th85rate_co : public mole_reaction
	{
	typedef mole_reaction_th85rate_co T;
	public:
		virtual T* Create() const {return new T;}
		
		virtual const char* name() {return "th85rate_co";}
		
		double rk() const
			{
				double esc_co, column;
				/******************************************************************************************
				 *	   First define the rate, which is of the form:
				 *	   
				 *		R = (Ro)*(Go*exp(-3.2Av))*Beta(tau(CO))
				 *
				 *	   where:
				 *
				 *	   Ro = 1.67*e-10
				 *	   (Go*exp(-3.2Av)) = hmi.UV_Cont_rel2_Habing_TH85_depth
				 *	   tauCO = 4.4e-15 * findspecies("CO")->column / (DopplerWidth/1e5) /
				 *		(1. + phycon.sqrte*0.6019); 
				 *       tauC = 1.6*e17*N(C)
				 *	   Beta(tau(CO)) = esca0k2(esc_co) 
				 ********************************************************************************************/   
				/* eqn 12 of 
				 * >>refer	CO	dissoc	Hollenbach, D.J., Takahashi, T., & Tielens, A. 1991, ApJ, 377, 192
				 * based on
				 * >>refer	CO	dissoc	Black, J.H., & van Dishoeck, E.F. 1988, ApJ, 334, 771 */
				
				if (this->reactants[0]->n_nuclei() != 0)
					column = mole.species[ this->reactants[0]->index ].column;
				else
					column = mole.species[ this->reactants[1]->index ].column;
				
				esc_co = 4.4e-15 * column / 
					/* the line width in km/s */
					(GetDopplerWidth(dense.AtomicWeight[5]+dense.AtomicWeight[7])/1e5) /
					/* this term accounts for populations within ground elec state */
					(1. + phycon.sqrte*0.6019);
				return esca0k2(esc_co)*th85rate(this);
			}
	};
	
	class mole_reaction_oh_c2h2_co_ch3 : public mole_reaction
	{
		typedef mole_reaction_oh_c2h2_co_ch3 T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "oh_c2h2_co_ch3";}
		double rk() const
			{
				/* This rate will blow up if the temperature gets too low, therefore 
				 * set this rate for T < 500 equal to the rate at 500 K */
				if(phycon.te > 500)
				{		
					return hmrate(this);
				}
				else
				{
					return 6.3E-18;
				}
			}
	};
	
	
	class mole_reaction_h_hnc_hcn_h : public mole_reaction
	{
		typedef mole_reaction_h_hnc_hcn_h T;
	public:
		virtual T* Create() const {return new T;}
		
		virtual const char* name() {return "h_hnc_hcn_h";}
		double rk() const
			{
				if(phycon.te > 100)
				{		
					return hmrate(this);
				}
				else
				{
					return 1e-15;
				}
			}
	};
}
#include "iso.h"
#include "h2.h"
double frac_H2star_hminus(void)
{
	/* >>chng 03 sep 09, get ratio of excited to ground state H2 */
	if( h2.lgEnabled  && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		return hmi.H2star_forms_hminus / 
			SDIV(hmi.H2star_forms_hminus+hmi.H2_forms_hminus);
		
		/* option print statement for above */
		/* printf( " H2s frac h- %.3e f(H2g) %.3e\n",frac_H2star_hminus ,
			hmi.H2_forms_hminus/SDIV(hmi.H2star_forms_hminus+hmi.H2_forms_hminus));*/
	}
	else
	{
		/* the large H2 molecule was not evaluated, so we can't use exact
		 * branching ratios.  These are the distribution fractions for around 500K */
		/*These depend on temperature and distribution function and the definition of ENERGY_H2_STAR.
		  So reset the values properly*/
		/* >>chng 05 jul 13, TE, with the new definition of H2s these are both 1 */
		/* >>chng 05 jul 31, activate above print, rest for current 0.5 ev defn */
		return 1. - 4.938e-6;
	}
}

namespace {
	double frac_H2star_grains(void)
	{
		/* >>chng 03 sep 09, get ratio of excited to ground state H2 */
		if( h2.lgEnabled  && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			return hmi.H2star_forms_grains / 
				SDIV(hmi.H2star_forms_grains+hmi.H2_forms_grains);
			
			/* option print statement for above */
			/*printf( "DEBUG H2s frac grain %.3e f(H2g) %.3e ",frac_H2star_grains ,
			  hmi.H2_forms_grains/SDIV(hmi.H2star_forms_grains+hmi.H2_forms_grains) ); */
		}
		else
		{
			/* the large H2 molecule was not evaluated, so we can't use exact
			 * branching ratios.  These are the distribution fractions for around 500K */
			/*These depend on temperature and distribution function and the definition of ENERGY_H2_STAR.
			  So reset the values properly*/
			/* >>chng 05 jul 13, TE, with the new definition of H2s these are both 1 */
			/* >>chng 05 jul 31, activate above print, rest for current 0.5 ev defn */
			return 0.9416;
		}
	}
	
	class mole_reaction_h2ph3p : public mole_reaction
	{
		typedef mole_reaction_h2ph3p T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2ph3p";}
		double rk() const
			{
				/* H2 + H2+ => H + H3+ */
				/** \todo	1	equivalent reaction for H2* is not included in chemistry, 
				 * Big h2 does not include this reaction, what to do? GS */
				
				if (!mole_global.lgLeidenHack)
					return 1.40e-9*(1. - sexp(9940./phycon.te));
				else
					return 2.08e-9;
		}
	};
}

	
namespace {
	class mole_reaction_hopexch : public mole_reaction
	{
		typedef mole_reaction_hopexch T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "hopexch";}
		double rk() const
			{
				return atmdat.CharExcRecTo[ipHYDROGEN][ipOXYGEN][0];
		}
	};
	
	class mole_reaction_hpoexch : public mole_reaction
	{
		typedef mole_reaction_hpoexch T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "hpoexch";}
		double rk() const
			{
				return atmdat.CharExcIonOf[ipHYDROGEN][ipOXYGEN][0];
		}
	};

	class mole_reaction_hmattach : public mole_reaction
	{
		typedef mole_reaction_hmattach T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "hmattach";}
		double rk() const
			{
				return hmirat(phycon.te) + (hmi.HMinus_induc_rec_rate)/SDIV(dense.eden);
			}
	};

	class mole_reaction_hmihph2p : public mole_reaction
	{
		typedef mole_reaction_hmihph2p T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "hmihph2p";}
		double rk() const
			{
				/* H- + H+ => H2+ + e
				 * equation (H6) from 
				 * >>refer	H2	chemistry	Galli,D., & Palla, F. 1998, A&A, 335, 403-420 
				 * hmihph2p = 6.9e-9f*(Tg)^(-0.35) for Tg<=8000
				 * hmihph2p = 6.9e-9f*(Tg)^(-0.9) for Tg>=8000  */
				if(phycon.te <= 7891.)
				{
					/*hmihph2p = 6.9e-9*pow(phycon.te , -0.35);*/
					return 6.9e-9 / (phycon.te30 * phycon.te05);
				}
				else 
				{
					/* >>chng 02 nov 18, had typo for leading coefficient here */
					/*hmihph2p = 9.6e-7*pow(phycon.te , -0.9);*/
					return 9.6e-7 / phycon.te90;
				}
			}
	};

	class mole_reaction_hmphoto : public mole_reaction
	{
		typedef mole_reaction_hmphoto T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "hmphoto";}
		double rk() const
			{
				return hmi.HMinus_photo_rate;
			}
	};

	double cionhm(const mole_reaction *)
	{
		/* collisional ionization of H-, rate from Janev, Langer et al.
		 * added factor hmi.exphmi for Boltzmann factor at low T
		 * Janev+ figure ends at 1e3 K */
		if( phycon.te < 3074. )
		{
			return 1.46e-32*(powi(phycon.te,6))*phycon.sqrte*hmi.exphmi;
		}
		else if( phycon.te >= 3074. && phycon.te < 30000. )
		{

			/* >>chng 03 mar 07, from above to below */
			/*cionhm = 5.9e-19*phycon.te*phycon.te*phycon.sqrte*phycon.te03*
			  phycon.te01*phycon.te01;*/
			return 5.9e-19*phycon.tesqrd*phycon.sqrte*phycon.te05;
		}
		else
		{
			return 1.54e-7;
		}

	}
	class mole_reaction_cionhm : public mole_reaction
	{
		typedef mole_reaction_cionhm T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "cionhm";}
		double rk() const
			{
				return cionhm(this);
			}
	};

	double assoc_detach(void)
	{
		/* form molecular hydrogen from H minus,
		 * associative detachment:  H- + H => H2 + E */
		/* make H2 from H- 
		 * associative detachment; H- + H => H2: 
		 * >>referold	H2	rates	Browne & Dalgarno J PHys B 2, 885 */
		/* rate coefficient from 
		 * >>refer	H2	form	Launay, J.R., Le Dourneuf, M., & Zeippen, C.J., 
		 * >>refercon	1991, A&A, 252, 842-852*/
		/* >>chng 02 oct 17, temp dependent fit to rate, updated reference,
		 * about 40% larger than before */
		double y , x;
		x = MAX2(10., phycon.te );
		x = MIN2(1e4, x );
		y=545969508.1323510+x*71239.23653059864;
		return 1./y;
	}

	class mole_reaction_c3bod : public mole_reaction
	{
		typedef mole_reaction_c3bod T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "c3bod";}
		double rk() const
			{
				return cionhm(this)*hmi.rel_pop_LTE_Hmin;
			}
	};

	class mole_reaction_asdfg : public mole_reaction
	{
		typedef mole_reaction_asdfg T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "asdfg";}
		double rk() const
			{
				return assoc_detach()*(1-frac_H2star_hminus());
			}
	};

	class mole_reaction_asdfs : public mole_reaction
	{
		typedef mole_reaction_asdfs T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "asdfs";}
		double rk() const
			{
				return assoc_detach()*frac_H2star_hminus();
			}
	};

	class mole_reaction_asdbg : public mole_reaction
	{
		typedef mole_reaction_asdbg T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "asdbg";}
		double rk() const
			{
				double ratio = mole_get_equilibrium_condition( this->label.c_str() );
				return assoc_detach() * ratio * (1.-frac_H2star_hminus());
			}
	};

	class mole_reaction_asdbs : public mole_reaction
	{
		typedef mole_reaction_asdbs T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "asdbs";}
		double rk() const
			{
				double ratio = mole_get_equilibrium_condition(this);
				return assoc_detach() * ratio * frac_H2star_hminus();
			}
	};

	class mole_reaction_bhneut : public mole_reaction
	{
		typedef  mole_reaction_bhneut T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "bhneut";}
		double rk() const
			{
				if( phycon.te > 1000. && dense.xIonDense[ipHYDROGEN][0] > 0.0 )
				{
					double ratio = mole_get_equilibrium_condition(this);
					/* HBN(3,1) is defined; when <HydTempLimit then set to 1 */
					/* mutual neut, mostly into n=3; rates from Janev et al
					 * H + H(n=3) => H- + H+ */
					/** \todo	2	process is net ionization term for H(n=3) states */
					/* this is the back reaction, forming H- from Ho */
					return (hneut(this)*ratio*
							  (iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH3s].Pop()+
								iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH3p].Pop()+
								iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH3d].Pop()) /
							  SDIV(dense.xIonDense[ipHYDROGEN][0]));
				}
				else
				{
					return 0.;
				}
			}
	};

	double hneut(const mole_reaction *)
	{
		if( phycon.te < 14125. )
		{
			/* the fit in Lepp et al. explodes at high temperature,
			 * Te = 14,125 is the temp where the rates reaches its lowest value */
			return 1.4e-7*pow(phycon.te/300,-0.487)*exp(phycon.te/29300);
		}
		else
		{
			return 3.4738192887404660e-008;
		}
	}
	class mole_reaction_hneut : public mole_reaction
	{
		typedef mole_reaction_hneut T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "hneut";}

		double rk() const
			{
				return hneut(this);
			}
	};

	class mole_reaction_h2_spon_diss : public mole_reaction
	{
		typedef mole_reaction_h2_spon_diss T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2_spon_diss";}
		double rk() const
			{
				return h2.spon_diss_tot;
			}
	};

	class mole_reaction_h2_ind_rad_diss : public mole_reaction
	{
		typedef mole_reaction_h2_ind_rad_diss T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2_ind_rad_diss";}
		double rk() const
			{
				return 0.;
			}
	};

	double grnh2tot(const mole_reaction *)
	{
		DEBUG_ENTRY( "grnh2tot()" );
		fixit("Remove factor of dense.gas_phase[ipHYDROGEN] factor from"
				"derivation of rate_h2_form_grains_used to avoid"
				"division here"); 
		if( mole.grain_area*dense.xIonDense[ipHYDROGEN][0]>0 )
			return gv.rate_h2_form_grains_used_total/(mole.grain_area*dense.xIonDense[ipHYDROGEN][0]);
		else
			return 0.;
	}
	class mole_reaction_grnh2tot : public mole_reaction
	{
	public:
		virtual const char* name() {return "grnh2tot";}
		double rk() const
			{
				return grnh2tot(this);
			}
	};

	class mole_reaction_grnh2 : public mole_reaction
	{
		typedef mole_reaction_grnh2 T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "grnh2";}
		double rk() const
			{
				return grnh2tot(this)*(1.-frac_H2star_grains());
			}
	};

	class mole_reaction_grnh2s : public mole_reaction
	{
		typedef mole_reaction_grnh2s T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "grnh2s";}
		double rk() const
			{
				return grnh2tot(this)*frac_H2star_grains();
			}
	};

	class mole_reaction_radasc : public mole_reaction
	{
		typedef mole_reaction_radasc T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "radasc";}
		double rk() const
			{
				// excited atom radiative association,
				// H(n=2) + H(n=1) => H2 + hnu
			
				if( dense.xIonDense[ipHYDROGEN][0] > 0. ) 
				{
					return hmrate(this) * 
						(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop() / dense.xIonDense[ipHYDROGEN][0] ) *
						(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop() + iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2s].Pop())/
						dense.xIonDense[ipHYDROGEN][0]; 
				}
				else 
					return 0.; 
			}
	};

	class mole_reaction_assoc_ion : public mole_reaction
	{
		typedef mole_reaction_assoc_ion T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "assoc_ion";}
		double rk() const
			{
				// excited atom associative ionization,
				// H(n=2) + H(n=1) => H2+ + e-
				if( dense.xIonDense[ipHYDROGEN][0] > 0. ) 
				{
					return hmrate(this) * 
						(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].Pop() / dense.xIonDense[ipHYDROGEN][0] ) *
						(iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop() + iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2s].Pop())/
						dense.xIonDense[ipHYDROGEN][0]; 
				}
				else 
					return 0.; 
			}
	};


	double rh2g_dis_h(const mole_reaction *)
	{
		if( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			return h2.Average_collH_dissoc_g;
		}
		else
		{
			double corr = MIN2(6.,14.44-phycon.alogte*3.08);

			if(corr > 0.)
				corr = exp10(corr*findspecieslocal("H")->den/(findspecieslocal("H")->den+1.6e4));
			else
				corr = 1.;
			/* must kill H- when UMIST is in place since they do not consider it */
			return 1.55e-8/phycon.sqrte*sexp(65107./phycon.te)* corr;
		}
	}
	class mole_reaction_rh2g_dis_h : public mole_reaction
	{
		typedef mole_reaction_rh2g_dis_h T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "rh2g_dis_h";}
		double rk() const
			{
				return rh2g_dis_h(this);
			}
	};

	double rh2s_dis_h(const mole_reaction *rate)
	{
		DEBUG_ENTRY( "rh2s_dis_h()" );
		if( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			return h2.Average_collH_dissoc_s;
		}
		else
		{
			if( ! fp_equal( rate->a, 1. ) )
			{
				fprintf( ioQQQ, "invalid parameter for rh2s_dis_h\n" );
				cdEXIT(EXIT_FAILURE);
			}
			return hmrate4(4.67e-7,-1.,5.5e4,phycon.te);
		}
	}
	class mole_reaction_rh2s_dis_h : public mole_reaction
	{
		typedef mole_reaction_rh2s_dis_h T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "rh2s_dis_h";}
		double rk() const
			{
				return rh2s_dis_h(this);
			}
	};

	double rh2g_dis_h2(const mole_reaction *rate)
	{
		DEBUG_ENTRY( "rh2g_dis_h2()" );
		if( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			return h2.Average_collH2_dissoc_g;
		}
		else
		{
			/* >>refer	H2	chemistry Palla, F., Salpeter, E.E., & Stahler, S.W., 1983, ApJ,271, 632-641 + detailed balance relation */
			if( ! fp_equal( rate->a, 1. ) )
			{
				fprintf( ioQQQ, "invalid parameter for rh2g_dis_h2\n" );
				cdEXIT(EXIT_FAILURE);
			}
			return hmrate4(5.5e-29*0.5/(SAHA*3.634e-5)*sqrt(300.),0.5,5.195e4,phycon.te); 
		}
	}
	class mole_reaction_rh2g_dis_h2 : public mole_reaction
	{
		typedef mole_reaction_rh2g_dis_h2 T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "rh2g_dis_h2";}
		double rk() const
			{
				return rh2g_dis_h2(this);
			}
	};

	double rh2s_dis_h2(const mole_reaction *rate)
	{
		DEBUG_ENTRY( "rh2s_dis_h2()" );
		if( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			return h2.Average_collH2_dissoc_s;
		}
		else
		{
			if( ! fp_equal( rate->a, 1. ) )
			{
				fprintf( ioQQQ, "invalid parameter for rh2s_dis_h2\n" );
				cdEXIT(EXIT_FAILURE);
			}
			return hmrate4(1e-11,0.,0.,phycon.te);
		}
	}
	class mole_reaction_rh2s_dis_h2 : public mole_reaction
	{
		typedef mole_reaction_rh2s_dis_h2 T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "rh2s_dis_h2";}
		double rk() const
			{
				return rh2s_dis_h2(this);
			}
	};

	double rh2s_dis_h2_nodeexcit(const mole_reaction *rate)
	{
		DEBUG_ENTRY( "rh2s_dis_h2_nodeexcit()" );
		if( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			return h2.Average_collH2_dissoc_s;
		}
		else
		{
			if( ! fp_equal( rate->a, 1. ) )
			{
				fprintf( ioQQQ, "invalid parameter for rh2s_dis_h2_nodeexcit\n" );
				cdEXIT(EXIT_FAILURE);
			}
			return hmrate4(1e-11,0.,2.18e4,phycon.te);
		}
	}
	class mole_reaction_rh2s_dis_h2_nodeexcit : public mole_reaction
	{
		typedef mole_reaction_rh2s_dis_h2_nodeexcit T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "rh2s_dis_h2_nodeexcit";}
		double rk() const
			{
				return rh2s_dis_h2_nodeexcit(this);
			}
	};

	double bh2g_dis_h(const mole_reaction *rate)
	{
		return rh2g_dis_h(rate)*hmi.rel_pop_LTE_H2g;
	}
	class mole_reaction_bh2g_dis_h : public mole_reaction
	{
		typedef mole_reaction_bh2g_dis_h T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "bh2g_dis_h";}
		double rk() const
			{
				return bh2g_dis_h(this);
			}
	};

	double bh2s_dis_h(const mole_reaction *rate)
	{
		return rh2s_dis_h(rate)*hmi.rel_pop_LTE_H2s;
	}
	class mole_reaction_bh2s_dis_h : public mole_reaction
	{
		typedef mole_reaction_bh2s_dis_h T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "bh2s_dis_h";}
		double rk() const
			{
				return bh2s_dis_h(this);
			}
	};

	double bh2g_dis_h2(const mole_reaction *rate)
	{
		return rh2g_dis_h2(rate)*hmi.rel_pop_LTE_H2g;
	}
	class mole_reaction_bh2g_dis_h2 : public mole_reaction
	{
		typedef mole_reaction_bh2g_dis_h2 T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "bh2g_dis_h2";}
		double rk() const
			{
				return bh2g_dis_h2(this);
			}
	};

	class mole_reaction_bh2s_dis_h2 : public mole_reaction
	{
		typedef mole_reaction_bh2s_dis_h2 T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "bh2s_dis_h2";}
		double rk() const
			{
				return rh2s_dis_h2(this)*hmi.rel_pop_LTE_H2s;
			}
	};

	class mole_reaction_h2photon : public mole_reaction
	{
		typedef mole_reaction_h2photon T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2photon";}
		double rk() const
			{
				return h2.photoionize_rate;
			}
	};

	class mole_reaction_h2crphot : public mole_reaction
	{


		typedef mole_reaction_h2crphot T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2crphot";}
		double rk() const
			{
				return secondaries.csupra[ipHYDROGEN][0]*2.02;
			}
	};

	class mole_reaction_h2crphh : public mole_reaction
	{
		typedef mole_reaction_h2crphh T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2crphh";}
		double rk() const
			{
				if( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
				{
					/* cosmic ray & secondary electron excitation to triplets
					 * big molecule scale factor changed 3-> 10
					 *>>chng 07 apr 08, from 3 to 10 to better capture results
					 * of Dalgarno et al 99 */
					//return secondaries.x12tot*10.;
					/*>>chng 20 jul 04, from 10 to 1.4 as discussed in
					 * https://ui.adsabs.harvard.edu/abs/2020RNAAS...4...78S/abstract
					 */
					return secondaries.x12tot*1.4;
				}
				else
				{
					/* the original TH85 approximation */
					return secondaries.x12tot*3.;
				}
			}
	};

	class mole_reaction_h2scrphh : public mole_reaction
	{
		typedef mole_reaction_h2scrphh T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2scrphh";}
		double rk() const
			{
				if( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
				{
					//return secondaries.x12tot*10.;
					//return secondaries.x12tot*10.;
					/*>>chng 20 jul 04, from 10 to 1.4 as discussed in
					 * https://ui.adsabs.harvard.edu/abs/2020RNAAS...4...78S/abstract
					 */
					return secondaries.x12tot*1.4;
				}
				else
				{
					return secondaries.x12tot*3.;
				}
			}
	};

	class mole_reaction_radath : public mole_reaction
	{
		typedef mole_reaction_radath T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "radath";}
		double rk() const
			{
				return MAX2(0.,2.325*MIN2(5000.,phycon.te)-1375.)*1e-20;
			}
	};

	class mole_reaction_gamtwo : public mole_reaction
	{
		typedef mole_reaction_gamtwo T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "gamtwo";}
		double rk() const
			{
				t_phoHeat dummyGrd, dummyExc;
				hmi.h2plus_exc_frac = dsexp( 10000./phycon.te );
				double fracGrd = MAX2( 0., 1. - hmi.h2plus_exc_frac );
				double coefGrd = GammaK(opac.ih2pnt[0],opac.ih2pnt[1],opac.ih2pof,1.,&dummyGrd);
				double coefExc = GammaK(opac.ih2pnt_ex[0],opac.ih2pnt_ex[1],opac.ih2pof_ex,1.,&dummyExc);
				hmi.h2plus_heatcoef = dummyGrd.HeatNet * fracGrd + dummyExc.HeatNet * hmi.h2plus_exc_frac;
				double coef = coefGrd * fracGrd + coefExc * hmi.h2plus_exc_frac;
				return coef;
			}
	};

	class mole_reaction_hlike_phot : public mole_reaction
	{
		typedef mole_reaction_hlike_phot T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "hlike_phot";}
		double rk() const
			{
				/* photoionization by hard photons, crossection =2*HI (wild guess)
				 * -- rjrw: where do they go??? 
				 * -- H3+ + hv => H2+ + H+ + e, best guess (P. Stancil, priv comm) */
			
				// on first pass this has not been set, do it here, but only do it once.
				if( !conv.nTotalIoniz )
					iso_photo( ipH_LIKE, ipHYDROGEN );
				return iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].gamnc;
			}
	};

	class mole_reaction_h2s_sp_decay : public mole_reaction
	{
		typedef mole_reaction_h2s_sp_decay T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2s_sp_decay";}
		double rk() const
			{
				/* >>chng 05 jul 11, TE, rename to use in punch file*/
				/* >>chng 05 jul 9, GS, use average A calculated from Big H2 */
				if( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
				{
					return h2.Average_A;
				}
				else
				{
					return 2e-7;
				}
			}
	};

	double h2_collh2_deexc(const mole_reaction *)
	{
		if( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			return h2.Average_collH2_deexcit;
		}
		else
		{
			return ((1.4e-12*phycon.sqrte * sexp( 18100./(phycon.te + 1200.) ))/6.); 
		}
	}
	class mole_reaction_h2_collh2_deexc : public mole_reaction
	{
		typedef mole_reaction_h2_collh2_deexc T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2_collh2_deexc";}
		double rk() const
			{
				return h2_collh2_deexc(this);
			}
	};

	double h2_collh_deexc(const mole_reaction *)
	{
		if ( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			return h2.Average_collH_deexcit;
		}
		else
		{
			return ((1e-12*phycon.sqrte * sexp(1000./phycon.te ))/6.); 
		}
	}
	class mole_reaction_h2_collh_deexc : public mole_reaction
	{
		typedef mole_reaction_h2_collh_deexc T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2_collh_deexc";}
		double rk() const
			{
				return h2_collh_deexc(this);
			}
	};

	double h2_collh2_excit(const mole_reaction *rate)
	{
		/* >>chng 05 jul 10, GS, use average collisional rate calculated from Big H2 */
		if ( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			return h2.Average_collH2_excit;
		}
		else
		{
			return h2_collh2_deexc(rate)*sexp( 30172./phycon.te); /* deexc_htwo*Boltz_fac_H2_H2star */
		}
	}
	class mole_reaction_h2_collh2_excit : public mole_reaction
	{
		typedef mole_reaction_h2_collh2_excit T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2_collh2_excit";}
		double rk() const
			{
				return h2_collh2_excit(this);
			}
	};

	double h2_collh_excit(const mole_reaction *rate)
	{
		/* >>chng 05 jul 10, GS, use average collisional rate calculated from Big H2 */
		if ( h2.lgEnabled && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
		{
			return h2.Average_collH_excit;
		}
		else
		{
			return h2_collh_deexc(rate)*sexp( 30172./phycon.te); /* deexc_hneut*Boltz_fac_H2_H2star */
		}
	}
	class mole_reaction_h2_collh_excit : public mole_reaction
	{
		typedef mole_reaction_h2_collh_excit T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2_collh_excit";}
		double rk() const
			{
				return h2_collh_excit(this);
			}
	};

	class mole_reaction_h2gexcit : public mole_reaction
	{
		typedef mole_reaction_h2gexcit T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2gexcit";}
		double rk() const
			{
				return hmi.H2_H2g_to_H2s_rate_used;
			}
	};

	class mole_reaction_h2sdissoc : public mole_reaction
	{
		typedef mole_reaction_h2sdissoc T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2sdissoc";};
		double rk() const
			{
				/* >>chng 03 sep 11, add this process */
				/* photo-destroy H2* by Solomon process at same rate as H2ground dissociation,
					see above eqn A12 in TH85 */
				/* >>chng 00 nov 25 factor of 0.1, assume pump is total, and 10% destroy H2 */
				/* >>chng 03 mar 07, had factor of 0.1 for branching ratio from H2** to H+H, 
				 * but branching is now already included */
				return hmi.H2_photodissoc_used_H2s +hmi.H2_Solomon_dissoc_rate_used_H2s;
			}
	};

	class mole_reaction_h2gdissoc : public mole_reaction
	{
		typedef mole_reaction_h2gdissoc T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "h2gdissoc";}
		double rk() const
			{
				/* >>chng 03 sep 11, add this process */
				/* photo-destroy H2* by Solomon process at same rate as H2ground dissociation,
					see above eqn A12 in TH85 */
				/* >>chng 00 nov 25 factor of 0.1, assume pump is total, and 10% destroy H2 */
				/* >>chng 03 mar 07, had factor of 0.1 for branching ratio from H2** to H+H, 
				 * but branching is now already included */
				return hmi.H2_photodissoc_used_H2g+hmi.H2_Solomon_dissoc_rate_used_H2g;
			}
	};

	class mole_reaction_hd_photodissoc : public mole_reaction
	{
		typedef mole_reaction_hd_photodissoc T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "hd_photodissoc";}
		double rk() const
			{
				return hd.photodissoc_BigH2_H2g + hd.Solomon_dissoc_rate_g;
			}
	};
}

/*hmirat compute radiative association rate for H- */
double hmirat(double te)
{
	double hmirat_v;
	
	DEBUG_ENTRY( "hmirat()" );
	
	// Compare published literature:
	// >>refer H- form deJong 1972A+A....20..263D
	// >>refer H- form Dalgarno 1973ApJ...181...95D
	// >>refer H- form Rawlings 1988MNRAS.230..695R

	/* fits to radiative association rate coefficients */
	if( te < 31.62 )
	{
		hmirat_v = 8.934e-18*phycon.sqrte*phycon.te003*phycon.te001*
			phycon.te001;
	}
	else if( te < 90. )
	{
		hmirat_v = 5.159e-18*phycon.sqrte*phycon.te10*phycon.te03*
				phycon.te03*phycon.te003*phycon.te001;
	}
	else if( te < 1200. )
	{
		hmirat_v = 2.042e-18*te/phycon.te10/phycon.te03;
	}
	else if( te < 3800. )
	{
		hmirat_v = 8.861e-18*phycon.te70/phycon.te03/phycon.te01*
			phycon.te003;
	}
	else if( te <= 1.4e4 )
	{
		/* following really only optimized up to 1e4 */
		hmirat_v = 8.204e-17*phycon.sqrte/phycon.te10/phycon.te01*
			phycon.te003;
	}
	else
	{
		/* >>chng 00 sep 28, add this branch */
		/* >>chng 14 mar 26, rate slightly discontinuous */
		hmirat_v = 5.70e-16*phycon.te20/phycon.te01;
	}
	
	return( hmirat_v );
}
STATIC void mole_h2_grain_form(void);

/**hmole_reactions - evaluates hydrogen chemistry reactions */
STATIC void mole_h_reactions(void);


namespace {
	class mole_reaction_gamheh : public mole_reaction
	{
		typedef mole_reaction_gamheh T;
	public:
		virtual T* Create() const {return new T;}
		virtual const char* name() {return "gamheh";}
		double rk() const
			{
				double retval = 0.;
				long int limit,
					i;
				/* photodissociation through 1.6->2.3 continuum */
			
				/* why is this in a look instead of GammaK?
				 * to fix must set opacities into stack */
			
				limit = MIN2(hmi.iheh2-1 , rfield.nflux );
				for( i=hmi.iheh1-1; i < limit; i++ )
				{
					retval += rfield.flux[0][i] + rfield.ConInterOut[i]+ rfield.outlin[0][i] + rfield.outlin_noplot[i];
				}
				retval *= 4e-18;
			
				/* hard radiation */
				retval += 3.*iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].gamnc;
			
				return retval;
			}
	};

	enum {exclude, base, umisthack, federman, lithium, deuterium, ti, misc, in_code, generated};
	static int source;


	static bool lgReactInitialized = false;
}

void mole_create_react( void )
{
	/* Should adaptively select reactions rather than list them explicitly here */
	
	/* prevent memory leaks */
	/* \todo	this is a temporary fix for PR14. We should improve the overall design
	 * of this code to prevent valid pointers being overwritten in a second call to mole_create_react */
	if( lgReactInitialized )
		return;
	
	lgReactInitialized = true;
	
	DEBUG_ENTRY( "mole_create_react()" );
	
	/* UDFA set 0 (false) above, UDFA = 1 (true) to make UDFA comparison and stop */
	if( UDFA )
		read_data("rate05_feb17.csv",parse_udfa);
	
	/* Set up registry of rate functions -- could do something intelligent about
	 * caching rate values by extending this API... */
	newfunc<mole_reaction_null>();
	{
		map<string,bool> canonical;
		for (map<string,bool>::iterator it=mole_global.offReactions.begin();
			  it != mole_global.offReactions.end(); ++it)
		{
			canonical[canonicalize_reaction_label(it->first.c_str())] = true;
		}
		mole_global.offReactions = canonical;
	}

	newfunc<mole_reaction_hmrate>();
	newfunc<mole_reaction_hmrate_exo>();
	newfunc<mole_reaction_constrate>();
	newfunc<mole_reaction_th85rate>();
	newfunc<mole_reaction_crnurate>();
	newfunc<mole_reaction_crnurate_noalbedo>();
	newfunc<mole_reaction_co_lnu_c_o_lnu>();
	newfunc<mole_reaction_vib_evap>();
	newfunc<mole_reaction_th85rate_co>();
	newfunc<mole_reaction_grn_abs>();
	/*  >chng 07 Dec 11, add grain surface chemical reactions */
	newfunc<mole_reaction_grn_react>();
	/* >>chng 07 Dec 10, add photodesorption of grain molecules */
	newfunc<mole_reaction_grn_photo>();
	newfunc<mole_reaction_oh_c2h2_co_ch3>();
	newfunc<mole_reaction_h_hnc_hcn_h>();
	
	newfunc<mole_reaction_gamheh>();
	newfunc<mole_reaction_hd_photodissoc>();
	newfunc<mole_reaction_h2gdissoc>();
	newfunc<mole_reaction_h2sdissoc>();
	newfunc<mole_reaction_h2gexcit>();
	newfunc<mole_reaction_h2_collh_excit>();
	newfunc<mole_reaction_h2_collh2_excit>();
	newfunc<mole_reaction_h2_collh_deexc>();
	newfunc<mole_reaction_h2_collh2_deexc>();	
	newfunc<mole_reaction_h2s_sp_decay>();
	newfunc<mole_reaction_hlike_phot>();
	newfunc<mole_reaction_gamtwo>();
	newfunc<mole_reaction_h2ph3p>();
	newfunc<mole_reaction_radath>();
	newfunc<mole_reaction_cryden_ov_bg>();
	newfunc<mole_reaction_h2scrphh>();
	newfunc<mole_reaction_h2crphh>();
	newfunc<mole_reaction_h2photon>();
	newfunc<mole_reaction_h2crphot>();
	newfunc<mole_reaction_rh2g_dis_h>();
	newfunc<mole_reaction_rh2s_dis_h>();
	newfunc<mole_reaction_rh2g_dis_h2>();
	newfunc<mole_reaction_rh2s_dis_h2>();
	newfunc<mole_reaction_rh2s_dis_h2_nodeexcit>();
	newfunc<mole_reaction_bh2g_dis_h>();
	newfunc<mole_reaction_bh2s_dis_h>();
	newfunc<mole_reaction_bh2g_dis_h2>();
	newfunc<mole_reaction_bh2s_dis_h2>();
	newfunc<mole_reaction_radasc>();
	newfunc<mole_reaction_assoc_ion>();
	newfunc<mole_reaction_grnh2>();
	newfunc<mole_reaction_grnh2s>();
	newfunc<mole_reaction_h2_spon_diss>();
	newfunc<mole_reaction_h2_ind_rad_diss>();
	newfunc<mole_reaction_bhneut>();
	newfunc<mole_reaction_hneut>();
	newfunc<mole_reaction_asdbs>();
	newfunc<mole_reaction_asdbg>();
	newfunc<mole_reaction_asdfs>();
	newfunc<mole_reaction_asdfg>();
	newfunc<mole_reaction_c3bod>();
	newfunc<mole_reaction_cionhm>();
	newfunc<mole_reaction_hmphoto>();
	newfunc<mole_reaction_hmihph2p>();
	newfunc<mole_reaction_hmattach>();
	//newfunc<mole_reaction_hpoexch>();
	//newfunc<mole_reaction_hopexch>();
	
	source = base;
	read_data("mole_co_base.dat",parse_base);
	/* default is true to use Federman rates instead of UMIST */
	if (mole_global.lgFederman)
	{
		source = federman;
		read_data("mole_co_federman.dat",parse_base);
	}
	/* default is false so we do read in umist hack, the Leiden workshop corrected UMIST data */
	if (!mole_global.lgLeidenHack) 
	{
		source = umisthack;
		read_data("mole_co_umisthack.dat",parse_base);
	}
	
	source = lithium;
	read_data("mole_lithium.dat",parse_base);
	
	source = deuterium;
	read_data("mole_deuterium.dat",parse_base);
	
#if 0
	source = ti;
	read_data("mole_ti.dat",parse_base);
#endif
	
	source = misc;
	read_data("mole_misc.dat",parse_base);

	/* Load null reaction to delete real rate from database */
	if (!mole_global.lgProtElim) 
	{
		source = exclude;
		newreact("C+,OH=>CO,H+","hmrate",0.,0.,0.); 
	}

	source = in_code;
	
	newreact("H,H,grn=>H2,grn","grnh2",1.,0.,0.);
	newreact("H,H,grn=>H2*,grn","grnh2s",1.,0.,0.);
	newreact("H-,PHOTON=>H,e-","hmphoto",1.,0.,0.);
	
	/* mutual neutralization with heavies, rate from Dalgarno and McCray
	 * all charged ions contribute equally,
		 * H- + A+ => H + A */
	fixit("this should be atom_list instead of unresolved_element_list, but we have not defined ionized species of all isotopes yet!!!");
	//for(vector<shared_ptr<chem_atom> >::iterator atom=atom_list.begin(); atom != atom_list.end(); ++atom)
	for(ChemNuclideList::iterator atom=nuclide_list.begin(); atom != nuclide_list.end(); ++atom)
	{
		if( !(*atom)->lgHasLinkedIon() )
			continue;
		long nelem = (*atom)->el()->Z-1;
		if( nelem >= ipHELIUM && dense.lgElmtOn[nelem] )
		{
			char react[32];
			sprintf(react,"H-,%s+=>H,%s", (*atom)->label().c_str(), (*atom)->label().c_str() );
			newreact(react,"hmrate",4e-6/sqrt(300.),-0.5,0.);
		}
	}
	
	newreact("H,e-=>H-,PHOTON","hmattach",1.,0.,0.);
	newreact("H-,H+=>H2+,e-","hmihph2p",1.,0.,0.);
	newreact("H-,e-=>H,e-,e-","cionhm",1.,0.,0.);
	newreact("H,e-,e-=>H-,e-","c3bod",1.,0.,0.);
	newreact("H,H-=>H2,e-","asdfg",1.,0.,0.);
	newreact("H,H-=>H2*,e-","asdfs",1.,0.,0.);
	newreact("H2,e-=>H,H-","asdbg",1.,0.,0.);
	newreact("H2*,e-=>H,H-","asdbs",1.,0.,0.);
	newreact("H-,H+=>H,H","hneut",1.,0.,0.);
	newreact("H,H=>H-,H+","bhneut",1.,0.,0.);

	if( hmi.lgLeiden_Keep_ipMH2s )
	{
		//newreact("H2*=>H,H,PHOTON","h2_spon_diss",1.,0.,0.);
		newreact("H2*,PHOTON=>H,H,PHOTON","h2_ind_rad_diss",1.,0.,0.);
		newreact("H,H2+=>H+,H2","hmrate",0.,0.,0.); // This process yields H2* plus, so if that species is tracked, the rate into H2 is zero.
	}
	else
	{
		//newreact("H2=>H,H,PHOTON","h2_spon_diss",1.,0.,0.);
		newreact("H2,PHOTON=>H,H,PHOTON","h2_ind_rad_diss",1.,0.,0.);
		newreact("H,H2+=>H+,H2","hmrate",6.4e-10f,0.,0.); /* this process populates v=4,no J information assume into J=0 -> H2s not H2g */
	}

	if (!mole_global.lgLeidenHack)
	{	
		newreact("H,H=>H2,PHOTON","radasc",3e-14,0.,0.);
		newreact("H,H=>H2+,e-","assoc_ion",3.27e-10,-0.35,17829.); /* >>refer	H2	rates	Rawlings J.M.C, Drew J.E, Barlow M. J., 1993, MNRAS 265, 968 */
		newreact("H2,H=>H,H,H","rh2g_dis_h",1.,0.,0.);
		newreact("H2,H2=>H,H,H2","rh2g_dis_h2",1.,0.,0.);
		newreact("H,H,H=>H2,H","bh2g_dis_h",1.,0.,0.); /* back rate, three body recombination, 2H + S => H_2 + S */
		newreact("H,H,H2=>H2,H2","bh2g_dis_h2",1.,0.,0.);
		//newreact("H,H,H2=>H2,H2","hmrate",5.5e-29/(8*300.),-1.,0.); /* >>refer	H2	chemistry Palla, F., Salpeter, E.E., & Stahler, S.W., 1983, ApJ,271, 632-641 */
		newreact("H2,PHOTON=>H2+,e-","h2photon",1.,0.,0.);
		newreact("H2,CRPHOT=>H2+,e-","h2crphot",1.,0.,0.);
		newreact("H2*,PHOTON=>H2+,e-","h2photon",1.,0.,0.);
		newreact("H2*,CRPHOT=>H2+,e-","h2crphot",1.,0.,0.);
		newreact("H2,CRPHOT=>H,H","h2crphh",1.,0.,0.); /* >>refer	H2	k	Shaw+20  2020RNAAS...4...78S */
		newreact("H2,CRPHOT=>H+,H-","cryden_ov_bg",3.9e-21,0.,0.); /* >>refer	H2	k	Millar, T.J., et.al, 1997,A&AS, 121, 139 */
		newreact("H2*,CRPHOT=>H+,H-","cryden_ov_bg",3.9e-21,0.,0.); /* >>refer	H2	k	Millar, T.J., et.al, 1997,A&AS, 121, 139 */
		newreact("H2,CRPHOT=>H+,H,e-","cryden_ov_bg",2.2e-19,0.,0.); /* >>refer	H2	k	Millar, T.J., et.al, 1997,A&AS, 121, 139 */
		newreact("H2*,CRPHOT=>H+,H,e-","cryden_ov_bg",2.2e-19,0.,0.); /* >>refer	H2	k	Millar, T.J., et.al, 1997,A&AS, 121, 139 */
		newreact("H+,H=>H2+,PHOTON","radath",1.,0.,0.);
		newreact("H2+,PHOTON=>H,H+","gamtwo",1.,0.,0.);
		newreact("H2+,CRPHOT=>H,H+","hlike_phot",1.,0.,0.);
		/* >> chng 05 aug 05, NPA comment.  This reaction is not in UMIST, for the case 
			of hard photons.  Turn if off for the comparison. */
		newreact("H3+,CRPHOT=>H2+,H+,e-","hlike_phot",2.,0.,0.);
		newreact("H,H+,e-=>H2+,e-","hmrate",2e-7 * SAHA * (4./(2.*1.)) * 3.634e-5 * pow(300.,-1.50),-1.50,0.);
	}
	else
	{
		newreact("H2,CRPHOT=>H2+,e-","hmrate",4.4e-17,0.,0.);
		newreact("H2,CRPHOT=>H,H+,e-","hmrate",1e-19,0.,0.);
		newreact("H2*,CRPHOT=>H,H+,e-","hmrate",1e-19,0.,0.);
		newreact("H2*,CRPHOT=>H2+,e-","hmrate",4.4e-17,0.,0.);
		newreact("H2,CRPHOT=>H,H","hmrate",5e-18,0.,0.);
		newreact("e-,H3+=>H2,H","hmrate",2.5e-8,-0.3,0.);	
		newreact("e-,H3+=>H,H,H","hmrate",7.5e-8,-0.3,0.);
		newreact("H+,H=>H2+,PHOTON","hmrate",5.3e-19,1.85,0);
		/* H2+ + e=> H + H; 
		 * equation (6,table 4) from
		 * >>refer	H2	l	Maloney et.al, 1996,ApJ, 466, 561 */
		newreact("H2+,e-=>H,H","hmrate",1.6e-8,-0.43,0.);
		newreact("H2+,PHOTON=>H,H+","th85rate",5.7e-10,0.,1.9);
	}
	
	newreact("H2*,CRPHOT=>H,H","h2scrphh",1.,0.,0.); /* >>refer	H2	k	Millar, T.J., et.al, 1997,A&AS, 121, 139 */
	newreact("H2,H2+=>H,H3+","h2ph3p",1.,0.,0.);
	newreact("H2*=>H2,PHOTON","h2s_sp_decay",1.,0.,0.);
	newreact("H2*,H2=>H2,H2","h2_collh2_deexc",1.,0.,0.);
	newreact("H2*,H=>H2,H","h2_collh_deexc",1.,0.,0.);
	newreact("H2,H2=>H2*,H2","h2_collh2_excit",1.,0.,0.);
	newreact("H2,H=>H2*,H","h2_collh_excit",1.,0.,0.);
	
	// collisional dissociation of H2*
	newreact("H2*,H=>H,H,H","rh2s_dis_h",1.,0.,0.);
	newreact("H2*,H2=>H2,H,H","rh2s_dis_h2",1.,0.,0.);
	newreact("H2*,H2*=>H2,H,H","rh2s_dis_h2",1.,0.,0.);
	newreact("H2*,H2*=>H2*,H,H","rh2s_dis_h2_nodeexcit",1.,0.,0.);
	
#if 0
	// more collision dissociation of H2
	newreact("H2,He=>He,H,H","rh2g_dis_h",1.,0.,0.);
	newreact("H2,H+=>H+,H,H","rh2g_dis_h",1.,0.,0.);
	newreact("H2,H3+=>H3+,H,H","rh2g_dis_h",1.,0.,0.);
	newreact("H2,e-=>e-,H,H","rh2g_dis_h",1.,0.,0.);
	newreact("H2*,He=>He,H,H","rh2s_dis_h",1.,0.,0.);
	newreact("H2*,H+=>H+,H,H","rh2s_dis_h",1.,0.,0.);
	newreact("H2*,H3+=>H3+,H,H","rh2s_dis_h",1.,0.,0.);
	newreact("H2*,e-=>e-,H,H","rh2s_dis_h",1.,0.,0.);
	
	// back reactions
	newreact("He,H,H=>H2,He","bh2g_dis_h",1.,0.,0.); 
	newreact("H+,H,H=>H2,H+","bh2g_dis_h",1.,0.,0.); 
	newreact("H3+,H,H=>H2,H3+","bh2g_dis_h",1.,0.,0.); 
	newreact("e-,H,H=>H2,e-","bh2g_dis_h",1.,0.,0.); 
	newreact("H,H,H=>H2*,H","bh2s_dis_h",1.,0.,0.);
	newreact("H2,H,H=>H2*,H2","bh2s_dis_h2",1.,0.,0.);
	newreact("H2,H,H=>H2*,H2*","bh2s_dis_h2",1.,0.,0.);
	newreact("H2*,H,H=>H2*,H2*","bh2s_dis_h2",1.,0.,0.);
	newreact("He,H,H=>H2*,He","bh2s_dis_h",1.,0.,0.); 
	newreact("H+,H,H=>H2*,H+","bh2s_dis_h",1.,0.,0.); 
	newreact("H3+,H,H=>H2*,H3+","bh2s_dis_h",1.,0.,0.); 
	newreact("e-,H,H=>H2*,e-","bh2s_dis_h",1.,0.,0.); 
#endif
	
	newreact("H2,PHOTON=>H2*","h2gexcit",1.,0.,0.);
	newreact("H2*,PHOTON=>H,H","h2sdissoc",1.,0.,0.);
	newreact("H2,PHOTON=>H,H","h2gdissoc",1.,0.,0.);
	newreact("HeH+,PHOTON=>H,He+","gamheh",1.,0.,0.);
	
	if(0)
	{
			// grain surface reactions
		newreact("OHgrn,Hgrn=>H2Ogrn","grn_react",1.,0.,0.);
	}
	
	if(UDFA)
	{
		fprintf(stderr,"Finished testing against UDFA database\n");
		cdEXIT(EXIT_FAILURE);
	}
	
	if (0) 
		plot_sparsity();
	
	source = generated;

	mole_check_reverse_reactions();
	
	if( deut.lgElmtOn )
		mole_generate_isotopologue_reactions( "H", "D" );
	
	long index = 0;
	for( mole_reaction_i it = mole_priv::reactab.begin(); it != mole_priv::reactab.end(); ++it, ++index )
		it->second->index = index;
	
	mole.reaction_rks.resize( index );
	mole.old_zone = -1;
	if( index > 0 )
		memset( &mole.reaction_rks[0], 0, (unsigned long)(index)*sizeof(double) );
	
	//  label catalytic, intra-group, and excitition/deexcitation vectors for all reactions
	for( mole_reaction_i it = mole_priv::reactab.begin(); it != mole_priv::reactab.end(); ++it )
		register_reaction_vectors( it->second );
}

STATIC void mole_generate_isotopologue_reactions( string atom_old, string atom_new )
{
	DEBUG_ENTRY( "mole_generate_isotopologue_reactions()" );
	
	bool lgDebug = false;

	shared_ptr<chem_nuclide> atomOld = findnuclide(atom_old.c_str());
	// store new reaction strings (and pointer to reaction derived from) and add all to reactab map after following iterator
	vector<string> newReactionStrings;
	vector<mole_reaction*> oldReactions;
	vector<realnum> branchingRatios;
	
	// iterate over all reactions and generate new strings by replacing atom_old with atom_new
	for( mole_reaction_ci it = mole_priv::reactab.begin(); it != mole_priv::reactab.end(); ++it )
	{
		bool lgFound = false;

		// Find number of atom sites in products
		int numSites = 0;
		for( long j=0; j<it->second->nproducts; ++j )
		{
			// Note that we must search before accessing the map because attempting to access an undeclared key will write to the map.
			if( it->second->products[j]->nNuclide.find( atomOld ) != it->second->products[j]->nNuclide.end() )
				numSites += it->second->products[j]->nNuclide[ atomOld ];
		}

		for( long i=0; i<it->second->nreactants; ++i )
		{
			// search for atom_old among reactants
			for( nNucs_i atom=it->second->reactants[i]->nNuclide.begin(); atom != it->second->reactants[i]->nNuclide.end(); ++atom )
			{
				if( atom->first->label() == atom_old )
				{
					lgFound = true;
					continue;
				}
			}
			if( lgFound )
				continue;
		}
		
		if( !lgFound )
			continue;
		
		if( lgDebug )
			fprintf( ioQQQ, "DEBUGGG mole_generate_isotopologue_reactions %s ..........\n", it->first.c_str() );
		
		for( long i=0; i<it->second->nreactants; ++i )
		{
			// ignore reactants with no nucleons (e-, PHOTON, CRPHOT, ... )
			if( it->second->reactants[i]->mole_mass < 0.99 * ATOMIC_MASS_UNIT )
				continue;
			
			vector<string> react_iso_labels;
			ChemNuclideList atomsLeftToRight;
			vector< int > numAtoms;
			string embellishments;
			// parse reactant label
			bool lgParseOK = parse_species_label( it->second->reactants[i]->label.c_str(), atomsLeftToRight, numAtoms, embellishments );
			if( !lgParseOK )
				TotalInsanity();
			// generate isotopologue labels
			create_isotopologues(
				atomsLeftToRight,
				numAtoms,
				atom_old,
				atom_new,
				embellishments,
				react_iso_labels );	
			for( unsigned j=0; j<react_iso_labels.size(); ++j )
			{	
				int numAtomsTot = 0;	
				for( long k=0; k<it->second->nproducts; ++k )
				{
					// ignore products with no nucleons (e-, PHOTON, CRPHOT, ... )
					if( it->second->products[k]->mole_mass < 0.99 * ATOMIC_MASS_UNIT )
						continue;
				
					atomsLeftToRight.clear();
					numAtoms.clear();
					embellishments.clear();
					// parse product label
					lgParseOK = parse_species_label( it->second->products[k]->label.c_str(), atomsLeftToRight, numAtoms, embellishments );
					ASSERT( lgParseOK == true );
					// Loop over all positions in this product
					for( unsigned position = 0; position < atomsLeftToRight.size(); ++position )
					{
						string prod_iso_label;
						// generate isotopologue labels
						create_isotopologues_one_position(
							position,
							atomsLeftToRight,
							numAtoms,
							atom_old,
							atom_new,
							embellishments,
							prod_iso_label );	

						if( prod_iso_label.empty() )
							continue;

						// Generate new reaction string
						string react_string;
						// first write reactants
						for( long i1=0; i1<i; ++i1 )
						{
							react_string += it->second->reactants[i1]->label;
							if( i1 != it->second->nreactants-1 )
								react_string += ",";
						}
						react_string += react_iso_labels[j];
						if( i != it->second->nreactants-1 )
							react_string += ",";
						for( long i2=i+1; i2<it->second->nreactants; ++i2 )
						{
							react_string += it->second->reactants[i2]->label;
							if( i2 != it->second->nreactants-1 )
								react_string += ",";
						}
						
						react_string += "=>";
						// now write products
						for( long k1=0; k1<k; ++k1 )
						{
							react_string += it->second->products[k1]->label;
							if( k1 != it->second->nproducts-1 )
								react_string += ",";
						}
						react_string += prod_iso_label;
						if( k != it->second->nproducts-1 )
							react_string += ",";
						for( long k2=k+1; k2<it->second->nproducts; ++k2 )
						{
							react_string += it->second->products[k2]->label;
							if( k2 != it->second->nproducts-1 )
								react_string += ",";
						}
						
						string canon_react_string = canonicalize_reaction_label( react_string.c_str() );
						// store new reaction string and pointer to old
						newReactionStrings.push_back( canon_react_string );
						oldReactions.push_back( it->second.get() );
						// This is the number of unique product lists given a particular reactant list.  Will divide rate below.
						branchingRatios.push_back( numAtoms[position]/numSites );
				
						if( lgDebug )
							fprintf( ioQQQ, "DEBUGGG mole_generate_isotopologue_reactions .................... %s\t\t(%2i/%2i).\n",
								canon_react_string.c_str(), numAtoms[position], numSites );

						numAtomsTot += numAtoms[position];
					}
				}
				ASSERT( numAtomsTot == numSites );
			}
		}
	}
	
	ASSERT( oldReactions.size() == newReactionStrings.size() );
	
	// now declare new reactions
	vector<mole_reaction*>::const_iterator it2 = oldReactions.begin();
	vector<realnum>::const_iterator it3 = branchingRatios.begin();
	for( vector<string>::const_iterator it1 = newReactionStrings.begin(); it1 != newReactionStrings.end(); ++it1, ++it2, ++it3 )
	{
		fixit("make adjustments to a for mass?");
		
		// don't overwrite existing reaction with these new auto-generated details
		if( mole_priv::reactab.find( it1->c_str() ) == mole_priv::reactab.end() )
		{
			ASSERT( *it3 <= 1. + FLT_EPSILON );
			fixit("multiply by *it3 here.  ASSERT at mole_reactions.cpp:1190 will blow currently.");
			newreact( it1->c_str(), (*it2)->name(), (*it2)->a /* * (*it3) */, (*it2)->b, (*it2)->c );
		}
	}
	
	return;
}

STATIC void mole_check_reverse_reactions(void)
{
	DEBUG_ENTRY( "mole_check_reverse_reactions()" );
	
	char chLabel[50], chLabelSave[50];
	int exists;
	
	for(mole_reaction_i p=mole_priv::reactab.begin(); 
			 p != mole_priv::reactab.end(); ++p) 
	{
		mole_reaction_i q = p;
		strcpy( chLabel, p->second->label.c_str() );
		strcpy( chLabelSave, p->second->label.c_str() );
		char *str = chLabel;
		const char *delim = "=>";
		char *chNewProducts = strtok( str, delim );
		char *chNewReactants = strtok( NULL, delim );
		char chNewReaction[50];
		
		strcpy( chNewReaction, chNewReactants );
		strcat( chNewReaction, "=>" );
		strcat( chNewReaction, chNewProducts );
		
		
		q = mole_priv::reactab.find(chNewReaction);
		
		exists = (q != mole_priv::reactab.end());
		if ( !exists )
		{
			if( trace.lgTraceMole )
			{
				fprintf(ioQQQ,"Warning! No reverse reaction for %30s.  ", p->second->label.c_str() );
				fprintf( ioQQQ,"\n" );
			}
			
			fixit("NB reverse reactions should be generated here");
		}
	}
	
	return;
}

STATIC double mole_get_equilibrium_condition( const char buf[] )
{
	DEBUG_ENTRY( "mole_get_equilibrium_condition()" );
	
	mole_reaction *rate = mole_findrate_s(buf);
	double result =	mole_get_equilibrium_condition( rate );
	return result;
}

STATIC double mole_get_equilibrium_condition( const mole_reaction* const rate )
{
	DEBUG_ENTRY( "mole_get_equilibrium_condition()" );

	if( !rate )
		return 0.;
 
	// multiply by reactant partition functions, then divide by products'
	double ln_result = 0.;
	for( long i=0; i<rate->nproducts; ++i)
	{
		double fac = mole_partition_function(rate->products[i]);
		if( fac==0. )
			return 0.;
		ln_result += log(fac);
	}
	for( long i=0; i<rate->nreactants; ++i)
	{
		double fac = mole_partition_function(rate->reactants[i]);
		if( fac==0. )
			return (double) BIGFLOAT;
		ln_result -= log(fac);
	}
	// Prevent overflow 
	double result = exp( MIN2( SEXP_LIMIT, ln_result ) );

	//fprintf( ioQQQ, "DEBUGGG equilibrium %20s\t%e\n", rate->label.c_str(), result );
	return result;
}

STATIC double mole_partition_function( const molecule* const sp)
{
	DEBUG_ENTRY( "mole_partition_function()" );
	
	if( sp->label == "PHOTON" || sp->label == "CRPHOT" )
	{
		fixit("How can we adapt existing structures to have a photon energy or range available here?");
		fixit("include 2hnu^3/c^2.  Understand HNU3C2 macro!");
		return 1.; //sexp( energy_ryd/phycon.te_ryd );
	}
	else if( sp->label == "CRP" || sp->label == "grn" )
		return 1.;
	
	fixit("need to figure out stat weight for any given particle");
	double q = 1.;
	// last factors convert kJ/mol to Kelvin/particle
	double deltaH = sp->form_enthalpy * (1e10/AVOGADRO/BOLTZMANN);
	ASSERT( sp->mole_mass > 0. );
	double part_fun = powpq(phycon.te*sp->mole_mass/(HION_LTE_POP*ELECTRON_MASS),3,2) * q * dsexp(deltaH/phycon.te);		
	ASSERT( part_fun < BIGFLOAT );	
	ASSERT( part_fun >= 0. );
	
	return part_fun;
}

void mole_cmp_num_in_out_reactions()
{
	DEBUG_ENTRY( "mole_cmp_num_in_out_reactions()" );
	
	vector<long> numForm, numDest;
	numForm.resize( mole_global.num_total );
	numDest.resize( mole_global.num_total );
	
	for(mole_reaction_i p=mole_priv::reactab.begin(); p != mole_priv::reactab.end(); p++) 
	{
		shared_ptr<mole_reaction> rate = p->second;
		for( long i=0; i<rate->nreactants; ++i)
		{
			++numDest[ rate->reactants[i]->index ];
		}
		
		for( long i=0; i<rate->nproducts; ++i)
		{
			++numForm[ rate->products[i]->index ];
		}
	}
	
	for( unsigned i=0; i<numForm.size(); ++i )
	{
		if( numForm[i]==0 && numDest[i]==0 )
			continue;
		if( numForm[i]>1 && numDest[i]>1 )
			continue;
		if( mole_global.list[i]->isMonatomic() )
			continue;
		fprintf( ioQQQ, "DEBUGGG mole_cmp_num_in_out_reactions %*s: in %4li out %4li\n", CHARS_SPECIES, mole_global.list[i]->label.c_str(), numForm[i], numDest[i] );
	}
	
	return;
}

STATIC char *getcsvfield(char **s,char c);
STATIC void parse_base(char *s)
{
	char *label, *reactstr, *f;
	double a, b, c;
	label = getcsvfield(&s,':');
	reactstr = getcsvfield(&s,':');
	f = getcsvfield(&s,':');
	a = atof(f);
	f = getcsvfield(&s,':');
	b = atof(f);
	f = getcsvfield(&s,':');
	c = atof(f);
	
	newreact(label,reactstr,a,b,c);
	
}

STATIC void newreact(const char label[], const char fun[], double a, double b, double c)
{
	DEBUG_ENTRY( "newreact()" );

	ratefunc rate = findfunc(fun);
	if (rate.get() == NULL) 
	{
		fprintf(stderr,"Rate function %s not known for reaction %s.  Aborting.  Sorry.\n",fun,label);
		cdEXIT( EXIT_FAILURE );
	}

	if( a < 0. )
	{
		fprintf(ioQQQ,"\n PROBLEM Reaction %s has negative pre-coefficient.  Aborting.  Sorry.\n", label);
		cdEXIT( EXIT_FAILURE );
	}

	rate->label = label;
	rate->a = a;
	rate->b = b;
	rate->c = c;
	rate->source = source;
	
	rate->photon = 0;
	
	if( parse_reaction( rate, label ) == 0 )
		return;
	
	canonicalize_reaction( rate );
	
	if( lgReactionTrivial( rate ) )
		return;

	if (mole_global.offReactions.find(rate->label)
		 != mole_global.offReactions.end())
	{
		fprintf(ioQQQ," W-reaction %s disabled\n",rate->label.c_str());
		phycon.lgPhysOK = false;
		return;
	}

	const char *rateLabelPtr = rate->label.c_str();
	
	ASSERT(lgReactBalance(rate)); /* Verify rate conserves particles and charge */
	
	rate->udfastate = ABSENT;
	if( UDFA)
	{
		compare_udfa(rate);
		if (rate->udfastate == ABSENT) 
		{
			fprintf(stderr,"Reaction %s not in UDFA\n",rateLabelPtr);
		}
	}
	
	/* >> chng 06 Oct 10 rjrw: use 1/(1/m1+1/m2) for reduced mass to prevent underflow */
	if (rate->nreactants == 2 && rate->reactants[0]->mole_mass!=0.0 && rate->reactants[1]->mole_mass!=0.0 )
	{
		rate->reduced_mass = 1./(1./rate->reactants[0]->mole_mass+1./rate->reactants[1]->mole_mass);
	}
	else
	{
		rate->reduced_mass = 0.;
	}
	
	/* If everything is OK, can register the reaction */
	mole_reaction_i p = mole_priv::reactab.find(rateLabelPtr);
	int exists = (p != mole_priv::reactab.end());
	// do not comment in release or beta version, or if NO TIMES entered
	if( exists && !t_version::Inst().lgRelease && prt.lgPrintTime )
	{
		/* Replace old rate */
		fprintf(ioQQQ,"Attention: duplicate reaction %s -- using new version\n",rateLabelPtr);
	}
	mole_priv::reactab[rateLabelPtr] = rate;
}

STATIC long parse_reaction( shared_ptr<mole_reaction>& rate, const char label[] )
{
	DEBUG_ENTRY( "parse_reaction()" );
	
	for (int i=0;i<MAXREACTANTS;i++)
	{
		rate->reactants[i] = NULL;
	}
	rate->nreactants = 0;
		
	for (int i=0;i<MAXPRODUCTS;i++)
	{
		rate->products[i] = NULL;
	}
	rate->nproducts = 0;
	
	bool lgProd = false;
	string buf = "";
	for(int i=0;!i || label[i-1]!='\0';i++) 
	{
		if(label[i] == ',' || label[i] == '=' || label[i] == '\0') 
		{
			molecule *sp = findspecies(buf.c_str());
			if( sp == null_mole || !sp->isEnabled ) 
			{
				if( trace.lgTraceMole )
					fprintf(ioQQQ,"Mole_reactions: ignoring reaction %s (species %s not active)\n",label,buf.c_str());
				return 0;
			}
			buf = "";
			if(! lgProd) 
			{
				if (rate->nreactants >= MAXREACTANTS) 
				{
					fprintf(stderr,"Mole_reactions: Too many reactants in %s, only %d allowed\n",label,MAXREACTANTS);
					cdEXIT(EXIT_FAILURE);
				}
				rate->reactants[rate->nreactants] = sp;
				rate->nreactants++;
			} 
			else 
			{
				if (rate->nproducts >= MAXPRODUCTS) 
				{
					fprintf(stderr,"Mole_reactions: Too many products in %s, only %d allowed\n",label,MAXPRODUCTS);
					cdEXIT(EXIT_FAILURE);
				}
				rate->products[rate->nproducts] = sp;
				rate->nproducts++;
			}
			if(label[i] == '=') 
			{
				i++; /* skip '>' as well */
				if (label[i] != '>')
				{
					fprintf(ioQQQ,"Format error in reaction %s\n",label);
					cdEXIT(EXIT_FAILURE);
				}
				lgProd = true;
			}
		} 
		else 
		{
			buf += label[i];
		}
	}

	for (int i=0;i<MAXREACTANTS;i++)
		if( rate->reactants[i] != NULL )
			++rate->reactants[i]->n_react;
	for (int i=0;i<MAXPRODUCTS;i++)
		if( rate->products[i] != NULL )
			++rate->products[i]->n_react;

	ASSERT( rate->nreactants );
	ASSERT( rate->nproducts );
	
	return 1;
}

STATIC string canonicalize_reaction_label( const char label[] )
{
	DEBUG_ENTRY( "canonicalize_reaction_label()" );
	
	// set up a dummy reaction
	ratefunc rate = findfunc("null");
	rate->label = label;
	parse_reaction( rate, label );
	canonicalize_reaction( rate );
	
	//if( !conv.nTotalIoniz && strcmp( label, rate->label.c_str() ) != 0 )
	//	fprintf( ioQQQ, "DEBUGGG reaction label %20s canonicalized to %20s\n", label, rate->label.c_str() );
	
	return rate->label;
}

STATIC void canonicalize_reaction( shared_ptr<mole_reaction>& rate )
{
	DEBUG_ENTRY( "canonicalize_reaction()" );
	
	/* Put species in canonical order to make sure reactions are unique.
		Can cause problems when a consistent ordering of species is
		required (look for references to "reactants[0]" for examples) */
	t_mole_global::sort(rate->reactants,rate->reactants+rate->nreactants);
	t_mole_global::sort(rate->products,rate->products+rate->nproducts);
	
	// now reorder label in same (new) order
	string newLabel;
	for( long i=0; i<rate->nreactants; ++i )
	{
		newLabel += rate->reactants[i]->label;
		if( i != rate->nreactants-1 )
			newLabel += ",";
	}
	newLabel += "=>";
	for( long i=0; i<rate->nproducts; ++i )
	{
		newLabel += rate->products[i]->label;
		if( i != rate->nproducts-1 )
			newLabel += ",";
	}
	
	// overwrite original label with canonical
	rate->label = newLabel;
	
	return;
}

// Returns true if reactant and product vectors are identical.
STATIC bool lgReactionTrivial( shared_ptr<mole_reaction>& rate )
{
	DEBUG_ENTRY( "lgReactionTrivial()" );

	bool lgTrivial = false;
	if( rate->nreactants == rate->nproducts )
	{
		lgTrivial = true;
		for( int k=0; k<rate->nreactants; ++k )
		{
			if( rate->reactants[k] != rate->products[k] )
			{
				lgTrivial = false;
				break;
			}
		}
	}

	return lgTrivial;
}

STATIC void register_reaction_vectors( shared_ptr<mole_reaction> rate )
{
	DEBUG_ENTRY( "register_reaction_vectors()" );
	
	for (long k=0;k<rate->nreactants;k++) 
	{
		rate->rvector[k] = NULL;
		rate->rvector_excit[k] = NULL;
	}
	
	for (long k=0;k<rate->nproducts;k++) 
	{
		rate->pvector[k] = NULL;
		rate->pvector_excit[k] = NULL;
	}
	
	/* Label catalytic species */
	for (long i=0;i<rate->nproducts;i++)
	{
		if (rate->pvector[i] == NULL)
		{
			for (long k=0;k<rate->nreactants;k++) 
			{
				if (rate->rvector[k] == NULL) 
				{
					if (rate->products[i] == rate->reactants[k])
					{
						rate->rvector[k] = rate->products[i];
						rate->pvector[i] = rate->reactants[k];
						break;
						}
				}
			}
		}
	}
	
	/* Label other intra-group transfers */
	for (long i=0;i<rate->nproducts;i++)
	{
		if (rate->pvector[i] == NULL)
		{
			for (long k=0;k<rate->nreactants;k++) 
			{
				if (rate->rvector[k] == NULL) 
				{
					if (rate->products[i]->groupnum != -1 &&
						 rate->products[i]->groupnum == 
						 rate->reactants[k]->groupnum)
					{
						rate->rvector[k] = rate->products[i];
						rate->pvector[i] = rate->reactants[k];
						break;
					}
				}
			}
		}
	}
	
	/* Label excited/deexcited pairs */
	for (long i=0;i<rate->nproducts;i++)
	{
		if (rate->pvector[i] == NULL && rate->pvector_excit[i] == NULL)
		{
			for (long k=0;k<rate->nreactants;k++) 
			{
				if (rate->rvector[k] == NULL && rate->rvector_excit[k] == NULL ) 
				{
					if ( lgDifferByExcitation( *rate->products[i], *rate->reactants[k] ) )
					{
						rate->rvector_excit[k] = rate->products[i];
						rate->pvector_excit[i] = rate->reactants[k];
						break;
					}
				}
			}
		}
	}

	return;
}


STATIC void plot_sparsity(void)
{ 
	FILE *sparsefp;
	int i, j, nb, ch;
	long int ratej;
	
	multi_arr<double,2> c(mole_global.num_total, mole_global.num_total);
	
	for(mole_reaction_i p=mole_priv::reactab.begin(); 
		 p != mole_priv::reactab.end(); ++p) 
	{
		mole_reaction &rate = *p->second;
		
		for (j=0;j<rate.nreactants;j++)
		{
			ratej = rate.reactants[j]->index;
			for (i=0;i<rate.nreactants;i++)
			{
				if (rate.rvector[i] == NULL)
					c[ratej][rate.reactants[i]->index] = 1.0;
			}
			for (i=0;i<rate.nproducts;i++)
			{
				if (rate.pvector[i] == NULL)
					c[ratej][rate.products[i]->index] = 1.0;
			}
		}
	}
	
	sparsefp = open_data("sparse.pbm","w");
	fprintf(sparsefp,"P4\n%d %d\n",
			  mole_global.num_total,mole_global.num_total);
	
	for ( j=0; j < mole_global.num_total; j++ )
	{
		nb = ch = 0;
		for ( i=0; i < mole_global.num_total; i++ )
		{
			ch = (ch << 1) | (c[i][j] != 0.0);
			nb++;
			if (nb == 8)
			{
				fputc(ch,sparsefp);
				nb = ch = 0;
			}
		}
		if (nb != 0) 
		{
			ch <<= 8-nb;
			fputc(ch,sparsefp);
		}
	}
	fclose(sparsefp);
}

#ifndef NDEBUG
STATIC bool lgReactBalance(const shared_ptr<mole_reaction> &rate)
{
	molecule::nNucsMap nel;
	int dcharge, n, sign;
	bool lgOK = true;
	
	dcharge = 0;
	for (n=0;n<rate->nreactants;n++) 
	{
		for( nNucs_i it = rate->reactants[n]->nNuclide.begin(); it != rate->reactants[n]->nNuclide.end(); ++it )
			nel[it->first] += it->second;
		dcharge +=  rate->reactants[n]->charge;
	}
	for (n=0;n<rate->nproducts;n++) 
	{
		for( nNucs_i it = rate->products[n]->nNuclide.begin(); it != rate->products[n]->nNuclide.end(); ++it )
			nel[it->first] -= it->second;
		dcharge -= rate->products[n]->charge;
	}
	if (dcharge != 0) 
	{
		fprintf(stderr,"Reaction %s charge out of balance by %d\n",
				  rate->label.c_str(),dcharge);
		lgOK = false;
	}
	
	for( nNucs_i it = nel.begin(); it != nel.end(); ++it )
	{
		if(it->second != 0)
		{
			if(it->second > 0)
				sign = 1;
			else
				sign = -1;
			fprintf(stderr,"Error: reaction %s %s %d of element %s\n",
					  rate->label.c_str(),sign==1?"destroys":"creates",
					  sign*it->second,
					  it->first->label().c_str() );
			lgOK = false;
		}
	}
	return lgOK;
}
#endif

enum {BUFSIZE=256};

namespace
{
	class formula_species {
	public:
		molecule *reactants[MAXREACTANTS], *products[MAXPRODUCTS];
	};
	
	bool operator< (const formula_species &a, const formula_species &b)
	{
		int i;
		for (i=0;i<MAXREACTANTS;i++)
		{
			if (a.reactants[i]<b.reactants[i])
				return true;
			else if (a.reactants[i]>b.reactants[i])
				return false;
		}
		for (i=0;i<MAXPRODUCTS;i++)
		{
			if (a.products[i]<b.products[i])
				return true;
			else if (a.products[i]>b.products[i])
				return false;
		}
		return false;
	}

	class udfa_reaction {
	public:
		int index;
		formula_species l;
		char source;            /* Calculated, Estimated, Literature compilation, Measured */
		double a, b, c, trange[2]; /* Overall scale and valid temperature range */
	};
}

static map <formula_species,shared_ptr<udfa_reaction> > udfatab;

STATIC void read_data(const char file[], void (*parse)(char *s))
{
	DEBUG_ENTRY( "read_data()" );
	char buf[BUFSIZE];

	FILE *fp = open_data(file,"r");
	if (!fp)
	{
		fprintf(stderr,"Error, could not read %s\n",file);
		cdEXIT(EXIT_FAILURE);
	}
	
	fixit("this seg-faults if file ends in blank line!");
	while(fgets(buf,BUFSIZE,fp)) 
	{
		if( buf[0] == '#' )
			continue;
		parse(buf);
	}
	fclose(fp);
}
#define FLTEQ(a,b) (fabs((a)-(b)) <= 1e-6*fabs((a)+(b)))
STATIC void parse_udfa(char *s)
{
	char *f;
	unsigned int havespecies = 1, i, n;
	/* lgCRPHOT is true for CRPHOT reactions as we change the data
	 * format for them */
	int lgCRPHOT=0;

	shared_ptr<udfa_reaction>r(new udfa_reaction);
	f = getcsvfield(&s,',');
	r->index = atoi(f);
	
	/* Load reactants */
	for (n=0;n<MAXREACTANTS;n++) 
	{
		r->l.reactants[n] = NULL; 
	}
	i = 0;
	for (n=0;n<MAXREACTANTS;n++) 
	{
		f = getcsvfield(&s,',');
		if (f[0] != '\0') 
		{
			i++;
			r->l.reactants[n] = findspecies(f);			
			if (r->l.reactants[n] == null_mole)
				havespecies = 0;
			if (!strncmp(f,"CRPHOT",6))
				lgCRPHOT = 1;
		}
	}
	t_mole_global::sort(r->l.reactants,r->l.reactants+i);

	/* Load products */
	for (n=0;n<MAXPRODUCTS;n++) 
	{  
		r->l.products[n] = NULL; /* Sentinel */
	}
	i = 0;
	for (n=0;n<MAXPRODUCTS;n++) 
	{  
		f = getcsvfield(&s,',');
		if (f[0] != '\0') 
		{
			i++;
			r->l.products[n] = findspecies(f);
			if (r->l.products[n] == null_mole)
				havespecies = 0;
		}
	}
	
	t_mole_global::sort(r->l.products,r->l.products+i);

	/* Load rate parameters */
	f = getcsvfield(&s,',');
	r->a = atof(f);
	f = getcsvfield(&s,',');
	r->b = atof(f);
	f = getcsvfield(&s,',');
	r->c = atof(f);

	if( lgCRPHOT) 
	{
		/* UDFA has a uniform value for the cosmic ray field independent of 
			circumstances which we correct -- verify they haven't changed 
			anything and move the multiplicative constant into our usual place 
			for it. */
		ASSERT (FLTEQ(r->a,1.3e-17));
		r->a = r->c;
		r->c = 0.;
	}

	/* Load data confidence and range of temperature validity */
	f = getcsvfield(&s,',');
	r->source = f[0]?f[0]:'?';
	for (n=0;n<2;n++) {
		f = getcsvfield(&s,',');
		r->trange[n] = atof(f);
	}

	if (havespecies)
	{
		if (udfatab.find(r->l) != udfatab.end())
		{
			fprintf(stderr,"Duplicate reaction\n");
		}
		udfatab[r->l] = r;
	}
}
STATIC char *getcsvfield(char **s, char c)
{
	char *sv, *f;
  
	sv = strchr(*s,c);
	if (sv) {
		*sv++ = '\0';
	}
	f = *s;
	*s = sv;
	return f;
}
STATIC void compare_udfa(const shared_ptr<mole_reaction> &rate)
{
	formula_species s;

	for (int n=0;n<MAXREACTANTS;n++) 
	{
		s.reactants[n] = rate->reactants[n]; 
	}
	for (int n=0;n<MAXPRODUCTS;n++) 
	{
		s.products[n] = rate->products[n]; 
	}
	auto p = udfatab.find(s);
	if (p == udfatab.end() )
	{
		/* fprintf(stdout,"Did not find reaction %s\n",rate->label); */
		return;
	}
	else
	{
		shared_ptr<udfa_reaction>& u = p->second;
		if (FLTEQ(rate->a,u->a) && FLTEQ(rate->b,u->b) && FLTEQ(rate->c,u->c))
		{
			rate->udfastate = CORRECT;
			/* fprintf(stdout,"Reaction %s agrees\n",rate->label); */
		} 
		else
		{
			rate->udfastate = CONFLICT;
			/* fprintf(stdout,"Reaction %18.18s clashes: a %9.2e %9.2e|b %9.2e %9.2e|c %9.2e %9.2e\n",
				rate->label,rate->a,u->a,rate->b,u->b,rate->c,u->c); */
		}
	}
}

namespace
{
	ratefunc findfunc (const char name[])
	{
		return shared_ptr<mole_reaction>(mole_priv::functab[name]->Create());
	}
}

void mole_update_rks(void)
{
	enum { DEBUG_MOLE = false };

	DEBUG_ENTRY( "mole_update_rks()" );

	mole_h2_grain_form();
	
	mole_h_reactions();
	
	for (mole_reaction_i p
			  =mole_priv::reactab.begin(); p != mole_priv::reactab.end(); ++p) 
	{
		mole_reaction &rate = *p->second;
		long index = rate.index;
		realnum newrk = rate.a*rate.rk();
		if (DEBUG_MOLE)
		{
			realnum oldrk = (realnum)mole.reaction_rks[index];
			if (fabs(newrk-oldrk) > 0.1*newrk)
				fprintf(ioQQQ,"%s: %15.8g => %15.8g\n",
						  rate.label.c_str(),oldrk,newrk);
		}
		mole.reaction_rks[index] = newrk;
	}	
}
void mole_rk_bigchange(void)
{
	DEBUG_ENTRY( "mole_rk_bigchange()" );
	enum { DEBUG_MOLE = false };

	if ( mole.old_reaction_rks.size() == 0 )
	{
		mole.old_zone = -1;
		mole.old_reaction_rks.resize(mole.reaction_rks.size());
	}

	if (nzone > 1)
	{
		ASSERT(mole.old_zone == nzone - 1);
		double bigchange = 0.;
		unsigned long bigindex = ULONG_MAX;
		for (unsigned long index=0; index<mole.reaction_rks.size(); ++index)
		{
			double oldrk = mole.old_reaction_rks[index],
				newrk = mole.reaction_rks[index],
				sum = oldrk+newrk, diff = newrk-oldrk;
			if (sum > 0.)
			{
				double change = fabs(diff)/sum;
				if (change > bigchange)
				{
					bigchange = change;
					bigindex  = index;
				}
			}
		}

		for (mole_reaction_i p
				  =mole_priv::reactab.begin(); p != mole_priv::reactab.end(); ++p) 
		{
			mole_reaction &rate = *p->second;
			if (bigindex == (unsigned) rate.index)
			{
				double oldrk = mole.old_reaction_rks[bigindex],
					newrk = mole.reaction_rks[bigindex];
				double pct = 0.;
				if (oldrk > 0.)
					pct = 100.*(newrk-oldrk)/oldrk;
				fprintf(ioQQQ, "Zone %ld, big chemistry rate change %s:"
						  " %15.8g => %15.8g (%6.2g%%)\n",
						  nzone,rate.label.c_str(),oldrk,newrk,pct);
				break;
			}
		}
	}
	
	mole.old_zone = nzone;
	for (unsigned long index=0; index<mole.reaction_rks.size(); ++index)
	{
		mole.old_reaction_rks[index] = mole.reaction_rks[index];
	}
}

STATIC void mole_h2_grain_form(void)
{
	DEBUG_ENTRY( "mole_h2_grain_form()" );

	double T_ortho_para_crit,xi_ELRD, beta_alpha_ELRD, recombination_efficiency_ELRD, Td;
	realnum AveVelH = GetAveVelocity( dense.AtomicWeight[ipHYDROGEN] );
	realnum AveVelH2 = GetAveVelocity( 2.f * dense.AtomicWeight[ipHYDROGEN] );

	/* H2 formation on grains;
	 * rate from 
	 * >>refer	H2	grain formation	Hollenbach, D., & McKee, C.F., 1979, ApJS, 41, 555 eq 3.4 3.8 */
	if( gv.lgDustOn() )
	{

		if (ENABLE_QUANTUM_HEATING)
		{
			fixit("Is this still necessary?");
			/* hmole is called before grains, so assure that all the grain stuff is properly initialized */
			GrainDrive();
		}

		/* these are rates (s-1) H2 will be deactivated by collisions with grains 
		 * will be incremented below 
		 * H2 ortho - para conversion on grain surface */
		h2.rate_grain_op_conserve = 0.;
		/* rate (s-1) v=0, J=1 level goes to 0 */
		h2.rate_grain_J1_to_J0 = 0.;

		/* loop over all grain species */
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			/* >>chng 02 feb 15, removed check tedust > 1.01, change in GrainsInit
			 * guarantees that all relevant parameters are initialized, PvH */

			double sticking_probability_H = sticking_probability_H_func( phycon.te, gv.bin[nd].tedust );

			bool lgUseQHeat = ENABLE_QUANTUM_HEATING && 
				gv.lgGrainPhysicsOn && gv.bin[nd].lgQHeat;
			long k, qnbin=0;
			vector<double> qtemp, qprob;
			if ( lgUseQHeat )
			{
				/* >>chng 04 feb 21, included quantum heating in calculation of formation rate, PvH */
				qtemp.resize(NQGRID);
				qprob.resize(NQGRID);
				
				qheat(qtemp,qprob,&qnbin,nd);
				
				if( gv.bin[nd].lgUseQHeat )
				{
					ASSERT( qnbin > 0 );
				}
				else
				{
					qnbin = 1;
					qprob[0] = 1.;
					qtemp[0] = gv.bin[nd].tedust;
				}
				
				gv.bin[nd].rate_h2_form_grains_HM79 = 0.;
				
				for( k=0; k < qnbin; k++ )
				{
					/* fraction of impacts that produce H2 before evaporation from grain surface.
					 * this is equation 3.4 of
					 * >>refer	grain	phys	Hollenbach, D.J., & McKee, C.F., 1979, ApJS, 41, 555
					 * 1e4 is ratio of total absorption sites to appropriate sites 
					 * 920 is D_H and chosen to get f_a = 0.5 at 100 K.
					 * factor of 0.6252 needed to obtain std ism rate to be 3e-17 at 100 K,
					 * the value deduced by
					 * >>refer	H2	grain physics	Jura, M., 1974, ApJ, 197, 581 */
					double conversion_efficiency_HM79 = 1/(1. + 1e4*sexp(920./qtemp[k]));
					sticking_probability_H = sticking_probability_H_func( phycon.te, qtemp[k] );
					
					gv.bin[nd].rate_h2_form_grains_HM79 += qprob[k] * sticking_probability_H *
						conversion_efficiency_HM79;
				}
				
				/* NB IntArea is total, not projected, area, must div by 4 */
				/* gv.bin[nd].rate_h2_form_grains_HM79 has units s^-1 since gv.bin[nd].cnv_H_pCM3 has units cm-3 */
				/* cnv_H_pCM3 converts <unit>/H (default depletion) -> <unit>/cm^3 (actual depletion), units are cm-3 */
				gv.bin[nd].rate_h2_form_grains_HM79 *= 0.5 * AveVelH *
					gv.bin[nd].IntArea/4. * gv.bin[nd].cnv_H_pCM3;
				
				ASSERT( gv.bin[nd].rate_h2_form_grains_HM79 > 0. );
			}
			else
			{
				/* fraction of impacts that produce H2 before evaporation from grain surface.
				 * this is equation 3.4 of
				 * >>refer	grain	phys	Hollenbach, D.J., & McKee, C.F., 1979, ApJS, 41, 555
				 * 1e4 is ratio of total absorption sites to appropriate sites 
				 * 920 is D_H and chosen to get f_a = 0.5 at 100 K.
				 * factor of 0.6252 needed to obtain std ism rate to be 3e-17 at 100 K,
				 * the value deduced by
				 * >>refer	H2	grain physics	Jura, M., 1974, ApJ, 197, 581 */
				double conversion_efficiency_HM79 = 1/(1. + 1e4*sexp(920./gv.bin[nd].tedust));

				/* NB IntArea is total area per H for default abundances, not projected area, must div by 4 
				 * units s^-1 since gv.bin[nd].cnv_H_pCM3 has units H cm-3 
				 * final units are cm s-1*/
				gv.bin[nd].rate_h2_form_grains_HM79 = 0.5 * AveVelH * gv.bin[nd].IntArea/4. * 
					/* cnv_H_pCM3 converts <unit>/H (default depletion) -> <unit>/cm^3 (actual depletion), units are cm-3 */
					gv.bin[nd].cnv_H_pCM3 * sticking_probability_H * conversion_efficiency_HM79;
				ASSERT( gv.bin[nd].rate_h2_form_grains_HM79 > 0. );
			}

			if( lgUseQHeat )
			{
				/* H2 formation on grains from 
				 * >>refer	H2	form	Cazaux, S., & Tielens, A.G.G.M., 2002, ApJ, 575, L29 */
				/* number of monolayers per second - only affects efficiency at very low or high temperatures */
				double f = 1e-10;
				/* equation 17 
					double sqrt_term = POW2( 1. + sqrt( (10000.-200.)/(600.-200.) ) );*/
				double sqrt_term = 35.399494936611667;

				gv.bin[nd].rate_h2_form_grains_CT02 = 0.;

				for( k=0; k < qnbin; k++ )
				{
					double beta_alpha = 0.25 * sqrt_term *sexp(200./qtemp[k] );
					/* equation 16 */
					double xi =  1./ (1. + 1.3e13*sexp(1.5*1e4/qtemp[k])*sqrt_term/(2.*f) );
					/* expression for beta comes from just after equation 5 */
					double beta = 3e12 * sexp( 320. / qtemp[k] );
					/* recombination efficiency given by their equation 15, they call
					 * this epsilon_H2 */
					double recombination_efficiency_CT02 = xi / (1. + 0.005*f/2./SDIV(beta) + beta_alpha );
					sticking_probability_H = sticking_probability_H_func( phycon.te, qtemp[k] );

					/* printf( " k %ld Td %.6e re*sp %.6e\n", k, qtemp[k], recombination_efficiency_CT02* */
					/* sticking_probability_H ); */

					gv.bin[nd].rate_h2_form_grains_CT02 += qprob[k] * sticking_probability_H *
						recombination_efficiency_CT02;
				}

				/* gv.bin[nd].IntArea integrated grain surface area Int(4pi*a^2), normalized per H, in cm^2/H,
				 * so x/4 is projected area of circle */
				/* gv.bin[nd].cnv_H_pCM3 is H density [cm-3] times grain depletion factor */
				/* gv.bin[nd].rate_h2_form_grains_CT02 units s-1 */
				gv.bin[nd].rate_h2_form_grains_CT02 *= 0.5 * AveVelH * 
					gv.bin[nd].IntArea/4. * gv.bin[nd].cnv_H_pCM3;

				ASSERT( gv.bin[nd].rate_h2_form_grains_CT02 > 0. );
			}
			else
			{
				/* H2 formation on grains from 
				 * >>refer	H2	form	Cazaux, S., & Tielens, A.G.G.M., 2002, ApJ, 575, L29 */
				/* number of monolayers per second - only affects efficiency at very low or high temperatures */
				double f = 1e-10;
				/* equation 17 
					double sqrt_term = POW2( 1. + sqrt( (10000.-200.)/(600.-200.) ) );*/
				double sqrt_term = 35.399494936611667;
				double beta_alpha = 0.25 * sqrt_term * exp(-200./gv.bin[nd].tedust);
				/* equation 16 */
				double xi =  1./ (1. + 1.3e13*exp(-1.5e4/gv.bin[nd].tedust)*sqrt_term/(2.*f) );
				/* expression for beta comes from just after equation 5 */
				double beta = 3e12 * exp(-320./gv.bin[nd].tedust);
				/* recombination efficiency given by their equation 15, they call
				 * this epsilon_H2 */
				double recombination_efficiency_CT02 = beta*xi / (beta + 0.0025*f + beta*beta_alpha);

				/* gv.bin[nd].IntArea integrated grain surface area Int(4pi*a^2), normalized per H, in cm^2/H,
				 * so x/4 is projected area of circle */
				/* gv.bin[nd].cnv_H_pCM3 is H density [cm-3] times grain depletion factor */
				/* units s-1 */
				gv.bin[nd].rate_h2_form_grains_CT02 = 0.5 * AveVelH * gv.bin[nd].IntArea/4. * 
					gv.bin[nd].cnv_H_pCM3 * sticking_probability_H * recombination_efficiency_CT02;
				ASSERT( gv.bin[nd].rate_h2_form_grains_CT02 >= 0. );
			}

			if( lgUseQHeat )
			{
				/* H2 formation on grains from 
 				 ** >>refer	H2	form	Rollig, M. et al., 2013, A&A, 549, A85 
 				 ** >>refer	H2	form	Cazaux, S.; Tielens, A. G. G. M. 2010, 2010ApJ...715..698C
				 ** improved treatment modifying CT rate above to include Eley-Rideal effect
 				 ** data and formalism comes from Appendix C of above reference
 				 ** this branch does quantum heating assuming rates simply scale as grain temperature */
				/* new rate variable ending stands for "Eley Rideal" */

				gv.bin[nd].rate_h2_form_grains_ELRD= 0.;

				if( gv.bin[nd].matType == MAT_CAR || gv.bin[nd].matType == MAT_CAR2 ||
				    gv.bin[nd].matType == MAT_SIC || gv.bin[nd].matType == MAT_PAH ||
				    gv.bin[nd].matType == MAT_PAH2 )
				{
					for( k=0; k < qnbin; k++ )
					{
						Td = qtemp[k];

						beta_alpha_ELRD = exp(-800./Td) / (0.5389970511202651 * exp(-540./Td) + 5.6333909478365e-14*sqrt(Td) ) ;
						xi_ELRD= 1./(1. + 4.61e24*sexp(45000./Td));
						recombination_efficiency_ELRD= (1. / (1. + beta_alpha_ELRD))*xi_ELRD;
								sticking_probability_H = 1./(1. + 0.04*sqrt(Td+phycon.te) +
								0.002*phycon.te + 8e-6*phycon.te*phycon.te);

						gv.bin[nd].rate_h2_form_grains_ELRD+= qprob[k] * sticking_probability_H *
								recombination_efficiency_ELRD;
					}

					/* gv.bin[nd].IntArea integrated grain surface area Int(4pi*a^2), normalized per H, in cm^2/H,
 					** so x/4 is projected area of circle */
					/* gv.bin[nd].cnv_H_pCM3 is H density [cm-3] times grain depletion factor */
					/* gv.bin[nd].rate_h2_form_grains_ELRD units s-1 */
	
					gv.bin[nd].rate_h2_form_grains_ELRD*= 0.5 * AveVelH * 
						gv.bin[nd].IntArea/4. * gv.bin[nd].cnv_H_pCM3;
	
					ASSERT( gv.bin[nd].rate_h2_form_grains_ELRD> 0. );
				}

				else if( gv.bin[nd].matType == MAT_SIL || gv.bin[nd].matType == MAT_SIL2 )
				{
 					for( k=0; k < qnbin; k++ )
					{
						Td = qtemp[k];

						beta_alpha_ELRD = exp(-450./Td) / (0.4266153643741715*exp(-340./Td) + 2.5335919594255e-14*sqrt(Td) ) ;
						xi_ELRD= 1./(1. + 7.00e24*sexp(45000./Td));
						recombination_efficiency_ELRD= (1. / (1. + beta_alpha_ELRD))*xi_ELRD;
								sticking_probability_H = 1./(1. + 0.04*sqrt(Td+phycon.te) +
								0.002*phycon.te + 8e-6*phycon.te*phycon.te);

						gv.bin[nd].rate_h2_form_grains_ELRD+= qprob[k] * sticking_probability_H *
						recombination_efficiency_ELRD;
					}

					/* gv.bin[nd].IntArea integrated grain surface area Int(4pi*a^2), normalized per H, in cm^2/H,
 					** so x/4 is projected area of circle */
					/* gv.bin[nd].cnv_H_pCM3 is H density [cm-3] times grain depletion factor */
					/* gv.bin[nd].rate_h2_form_grains_ELRD units s-1 */
					gv.bin[nd].rate_h2_form_grains_ELRD*= 0.5 * AveVelH * 
						gv.bin[nd].IntArea/4. * gv.bin[nd].cnv_H_pCM3;
	
					ASSERT( gv.bin[nd].rate_h2_form_grains_ELRD> 0. );
				}
			}
			else
			{
				/* H2 formation on grains from 
 				** >>refer	H2	form	Rollig, M. et al., 2013, A&A, 549, A85 
 				** improved treatment modifying CT rate above to include Eley-Rideal effect
 				** data and formalism comes from Appendix C of above reference */
				/* new rate variable ending stands for "Eley Rideal" */
				gv.bin[nd].rate_h2_form_grains_ELRD= 0.;

				if( gv.bin[nd].matType == MAT_CAR || gv.bin[nd].matType == MAT_CAR2 ||
				    gv.bin[nd].matType == MAT_SIC || gv.bin[nd].matType == MAT_PAH ||
				    gv.bin[nd].matType == MAT_PAH2 )
				{
					Td = gv.bin[nd].tedust;

					beta_alpha_ELRD = exp(-800./Td) / (0.5389970511202651 * exp(-540./Td) + 5.6333909478365e-14*sqrt(Td) ) ;
					xi_ELRD= 1./(1. + 4.61e24*sexp(45000./Td));
					recombination_efficiency_ELRD = (1. / (1. + beta_alpha_ELRD))*xi_ELRD;

					/* gv.bin[nd].IntArea integrated grain surface area Int(4pi*a^2), normalized per H, in cm^2/H,
 					** so x/4 is projected area of circle */
					/* gv.bin[nd].cnv_H_pCM3 is H density [cm-3] times grain depletion factor */
					/* units s-1 */
					gv.bin[nd].rate_h2_form_grains_ELRD = 0.5 * AveVelH * gv.bin[nd].IntArea/4. * 
					gv.bin[nd].cnv_H_pCM3 * sticking_probability_H * recombination_efficiency_ELRD;

					ASSERT( gv.bin[nd].rate_h2_form_grains_ELRD > 0. );
				}

				else if( gv.bin[nd].matType == MAT_SIL || gv.bin[nd].matType == MAT_SIL2 )
				{
					Td = gv.bin[nd].tedust;

					beta_alpha_ELRD = exp(-450./Td) / (0.4266153643741715*exp(-340./Td) + 2.5335919594255e-14*sqrt(Td) ) ;
					xi_ELRD= 1./(1. + 7.00e24*sexp(45000./Td));
					recombination_efficiency_ELRD = (1. / (1. + beta_alpha_ELRD))*xi_ELRD;

					/* gv.bin[nd].IntArea integrated grain surface area Int(4pi*a^2), normalized per H, in cm^2/H,
 					** so x/4 is projected area of circle */
					/* gv.bin[nd].cnv_H_pCM3 is H density [cm-3] times grain depletion factor */
					/* units s-1 */
					gv.bin[nd].rate_h2_form_grains_ELRD = 0.5 * AveVelH * gv.bin[nd].IntArea/4. * 
					gv.bin[nd].cnv_H_pCM3 * sticking_probability_H * recombination_efficiency_ELRD;

					ASSERT( gv.bin[nd].rate_h2_form_grains_ELRD > 0. );
				}
			}

			if (ENABLE_QUANTUM_HEATING)
			{
				/* reset sticking probability for code below */
				sticking_probability_H = sticking_probability_H_func( phycon.te, gv.bin[nd].tedust );
			}

			/* rate (s-1) all H2 v,J levels go to 0 or 1, preserving nuclear spin */
			/* ortho to para on grain surfaces, taken from 
			 *>refer	H2	sticking	Le Bourlot, J., 2000, A&A, 360, 656-662 
			 * >chng 05 apr 30, GS, hmi.H2_total/dense.gas_phase[ipHYDROGEN] is removed
			 * This is used in h2.c.
			 * NB IntArea is total are per H, not projected area, must div by 4 
			 * gv.bin[nd].cnv_H_pCM3 has units H cm-3 to product with above
			 * is cm2/H H/cm3 or cm-1 or an opacity
			 * multiply by velocity of H2, cm s-1, so product 
			 * h2.rate_grain_op_conserve has units s^-1  */
			h2.rate_grain_op_conserve += AveVelH2 * gv.bin[nd].IntArea/4. *
				gv.bin[nd].cnv_H_pCM3 * sticking_probability_H;

			/* ortho to para on grain surfaces, taken from 
			 *>refer	H2	sticking	Le Bourlot, J., 2000, A&A, 360, 656-662 
			 * For all grain temperatures, this process corresponds to high J going to
			 * either 0 or 1 preserving nuclear spin.  All ortho go to 1 and para go to 0.
			 * When the dust temperature is below Tcrit all 1 go to 0 and so all J go to 0.

			 * this temperature depends on grain composition, discussion left column of page 657,
			 * this is for a bare grain */
			/** \todo	2	- put in actual composition dependent Tad - this is only valid 
			 * for bare surfaces - not ice - for ice Tad is 555K 
			 * hmi.Tad is binding energy expressed as a temperature 
			 * note that hmi.Tad is set to 800. in zero 
			 * tau_nu the first equation in section 2.5
			 * equation one paragraph before equation 2 
			 * at low grain temperatures all end in para, J=0 */

			/* AveVelH2 is average speed of H2 molecules 
			 * for now assume that sticking probability for H2 on the grain is equal to
			 * that for H 
			 * efficiency factor efficiency_opr is vary fast function of t dust - 
			 * large at low Td and small at Td > T_ortho_para_crit
			 * start evaluating just above the critical temperature 
			 * T_ortho_para_crit this is roughly 24.345 K,GS */
			T_ortho_para_crit = 2. * hmi.Tad / log( POW2(60. *1.1e11)*hmi.Tad); 
			if( gv.bin[nd].tedust < T_ortho_para_crit )
			{
				double efficiency_opr = sexp(60.*1.1e11*sqrt(hmi.Tad)*sexp(hmi.Tad/gv.bin[nd].tedust));
				/* rate (s-1) all v,J levels go to 0, regardless of nuclear spin
				 * see above discussion for how units work out */
				h2.rate_grain_J1_to_J0 += AveVelH2 * gv.bin[nd].IntArea/4. * 
					gv.bin[nd].cnv_H_pCM3 * sticking_probability_H * efficiency_opr;
			}
		}
		/*fprintf(ioQQQ," H2 grain form rate HM79 %.2e  %.2e CT02 %.2e  %.2e O-P grn %.2e %.2e\n", 
		  gv.bin[nd].rate_h2_form_grains_HM79 , 
		  gv.bin[nd].rate_h2_form_grains_HM79 ,
		  gv.bin[nd].rate_h2_form_grains_CT02 , 
		  gv.bin[nd].rate_h2_form_grains_CT02 , 
		  h2.rate_grain_J1_to_J0,
		  hmi.rate_h2_allX_2_J1_grains
		  );*/
		/* options to turn off grain collision with atom h2 collisions grains off command */
		h2.rate_grain_op_conserve *= h2.lgH2_grain_deexcitation;
		h2.rate_grain_J1_to_J0 *= h2.lgH2_grain_deexcitation;

	}
	else
	{
		/* grains are not enabled, set these to zero */
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			gv.bin[nd].rate_h2_form_grains_CT02 = 0.;
			gv.bin[nd].rate_h2_form_grains_HM79 = 0.;
		}
		/* rate all H2 goes to either 0 or 1 depending on ortho/para */
		h2.rate_grain_op_conserve = 0.;
		/* at low temp, rate all H2 goes to J=0 */
		h2.rate_grain_J1_to_J0 = 0.;
	}

	/* the H2 catalysis rate on grains that is actually used in calculations
	 * hmi.ScaleJura is scale factor set with set Jura scale command 
	 * units are s-1 
	 * default is 'C' Cazaux & Tielens */
	gv.rate_h2_form_grains_used_total = 0.;
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		if( hmi.chJura == 'C' )
		{
			/* use the new rate by 
			 * >>refer	H2	form	Cazaux, S., & Tielens, A.G.G.M., 2002, ApJ, 575, L29 
			 * units are s-1*/
			gv.bin[nd].rate_h2_form_grains_used = 
				gv.bin[nd].rate_h2_form_grains_CT02*hmi.ScaleJura;
			gv.rate_h2_form_grains_used_total += gv.bin[nd].rate_h2_form_grains_used;
		}

		if( hmi.chJura == 'E' )
		{
			/* use the new revised CT02 rate from
 			** >>refer	H2	form	Rollig, M. et al., 2013, A&A, 549, A85 
 			** units are s-1 */
			gv.bin[nd].rate_h2_form_grains_used = 
				gv.bin[nd].rate_h2_form_grains_ELRD*hmi.ScaleJura;
			gv.rate_h2_form_grains_used_total += gv.bin[nd].rate_h2_form_grains_used;
		}

		else if( hmi.chJura == 'T' )
		{
			/* rate from Hollenbach & McKee 1979  */
			gv.bin[nd].rate_h2_form_grains_used = 
				gv.bin[nd].rate_h2_form_grains_HM79*hmi.ScaleJura;
			gv.rate_h2_form_grains_used_total += gv.bin[nd].rate_h2_form_grains_used;
		}
		else if( hmi.chJura == 'S' )
		{
			/* H2 formation rate from 
			 * >>refer	H2	form	Sternberg, A. & Neufeld, D.A. 1999, ApJ, 516, 371 */
			gv.bin[nd].rate_h2_form_grains_used = 
				3e-18 * phycon.sqrte / gv.bin.size() * dense.gas_phase[ipHYDROGEN]*hmi.ScaleJura;
			/* this is simple rate from Sternberg & Neufeld 99 */
			gv.rate_h2_form_grains_used_total += gv.bin[nd].rate_h2_form_grains_used;
		}
		/*>>chng 07 jan 10, this had been C for constant, and so could never have been triggered.
		 * caught by robin Williams, fixed by nick Abel, error was in sense that any set jura rate
		 * would use Cazaux & Tielens */
		else if( hmi.chJura == 'F' )
		{
			/* command "set H2 rate" - enters log of Jura rate - C for constant,
			 * no dependence on grain properties */
			gv.bin[nd].rate_h2_form_grains_used = hmi.rate_h2_form_grains_set*dense.gas_phase[ipHYDROGEN] / gv.bin.size();
			gv.rate_h2_form_grains_used_total += gv.bin[nd].rate_h2_form_grains_used;
		}
	}
	ASSERT( gv.rate_h2_form_grains_used_total >= 0. );

	if (ENABLE_QUANTUM_HEATING)
	{
		fprintf(ioQQQ, " fnzone %.2f H2 rate %.4e\n", fnzone, gv.rate_h2_form_grains_used_total );
	}

	/* print rate coefficient */
	/*fprintf(ioQQQ," total grain h2 form rate %.3e\n",gv.rate_h2_form_grains_used_total);*/

}
/*mole_h_reactions update mole reactions for H */
STATIC void mole_h_reactions( void )
{
	static double teused=-1;
	double exph2,
		exph2s,
		exphp,
		ex3hp;
	long i;
	double h2esc, 
		th2,
		cr_H2s ,
		cr_H2dis,
		cr_H2dis_ELWERT_H2g,
		cr_H2dis_ELWERT_H2s;

	DEBUG_ENTRY( "mole_h_reactions()" );

	/* everything here depends only on temperature - don't do anything if we don't
	 * need to */

	bool need_update = ! fp_equal( phycon.te, teused );

	if( need_update )
	{
		teused = phycon.te;

		/* get LTE populations */
		/* related to population of H- in LTE
		 * IP is 0.754 eV */
		hmi.exphmi = sexp(8.745e3/phycon.te);
		if( hmi.exphmi > 0. )
		{
			/* these are ratio n*(H-)/[  n*(ne) n*(Ho)  ] */
			hmi.rel_pop_LTE_Hmin = SAHA/(phycon.te32*hmi.exphmi)*(1./(2.*2.));
		}
		else
		{
			hmi.rel_pop_LTE_Hmin = 0.;
		}

		/* population of H2+ in LTE, hmi.rel_pop_LTE_H2p is H_2+/H / H+
		 * dissociation energy is 2.647 */
		exphp = sexp(3.072e4/phycon.te);
		if( exphp > 0. )
		{
			/* stat weight of H2+ is 4
			 * last factor was put in ver 85.23, missing before */
			hmi.rel_pop_LTE_H2p = SAHA/(phycon.te32*exphp)*(4./(2.*1.))*3.634e-5;
		}
		else
		{
			hmi.rel_pop_LTE_H2p = 0.;
		}

		/* related to population of H3+ in LTE, hmi.rel_pop_LTE_H3p is H_3+/( H2+ H+ )
		 * dissociation energy is 2.647 */
		ex3hp = sexp(1.882e4/phycon.te);
		if( ex3hp > 0. )
		{
			/* stat weight of H2+ is 4
			 * last factor was put in ver 85.23, missing before */
			hmi.rel_pop_LTE_H3p = SAHA/(phycon.te32*ex3hp)*(4./(2.*1.))*3.634e-5;
		}
		else
		{
			hmi.rel_pop_LTE_H3p = 0.;
		}
	}
	/* end constant temperature - */

	// Big H2 rates are dependent on population as well as temperature
	/* population of H2 in LTE
	 * dissociation energy of H2g is 4.477eV, for TH85 model */
	/* >>chng 05 oct 17, GS, averaged in big H2 in order to consider correct statistical weight*/
	if( h2.lgEnabled  && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		/* the terms on the right were computed in the large molecule */
		hmi.rel_pop_LTE_H2g = h2.rel_pop_LTE_g;
		hmi.rel_pop_LTE_H2s = h2.rel_pop_LTE_s;
	}
	else
	{
		if (need_update)
		{
			/* H2 ground */
			exph2 = sexp((5.195e4)/phycon.te);
			/* >>chng 05 oct 17, GS, note that statical weight of H2g is assumed to be 1 if big H2 is not turned on*/
			
			if( exph2 > 0. ) 
			{
				/* >>chng 05 oct 17, GS, note that statical weight of H2g is assumed to be 1 if big H2 is not turned on*/
				hmi.rel_pop_LTE_H2g = SAHA/(phycon.te32*exph2)*(1./(2.*2.))*3.634e-5;
			}
			else
			{
				hmi.rel_pop_LTE_H2g = 0.;
			}
			
			/* H2 star */
			/* population of H2s in LTE
			 * dissociation energy is 1.877eV, if h2s = 2.6eV, assumed for TH85 model */
			/* >>chng 05 oct 17, GS, averaged in big H2 in order to consider correct statistical weight*/
			exph2s = sexp(2.178e4/phycon.te);
			
			if( exph2s > 0. ) 
			{
				/* >>chng 05 oct 17, GS, note that statical weight of H2s is assumed to be 1 if big H2 is not turned on*/
				hmi.rel_pop_LTE_H2s = SAHA/(phycon.te32*exph2s)*(1./(2.*2.))*3.634e-5;
			}
			else
			{
				hmi.rel_pop_LTE_H2s = 0.;
			}
		}
	}
	{
		/*@-redef@*/
		/* often the H- route is the most efficient formation mechanism for H2,
		 * will be through rate called Hmolec_old[ipMH]*hmi.assoc_detach
		 * this debug print statement is to trace h2 oscillations */
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && nzone>187&& iteration > 1/**/)
		{
			/* rapid increase in H2 density caused by rapid increase in hmi.rel_pop_LTE_H2g */
			fprintf(ioQQQ,"ph2lteee\t%.2e\t%.1e\t%.1e\n", 
					  hmi.rel_pop_LTE_H2g, 
					  sexp(2.178e4/phycon.te),
					  phycon.te);
		}
	}
	
	/* cooling due to radiative attachment */
	hmi.hmicol = hmirat(phycon.te)*EN1RYD*phycon.te*1.15e-5;

	fixit("Wasted cycles if we don't use Stancil's rates below?  Why not put this down there/if used?");
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
	{
		if( (*diatom)->lgEnabled && mole_global.lgStancil )
			(*diatom)->Mol_Photo_Diss_Rates();
	}

	/*fprintf(ioQQQ,"%.2e %.2e %.2e %.2e\n", phycon.te, hmi.hminus_rad_attach , hmi.hmicol,
	  hmi.hmicol/(hmi.hminus_rad_attach*EN1RYD*phycon.te*1.15e-5) );*/

	/* get per unit vol */
	hmi.hmicol *= dense.eden*findspecieslocal("H")->den;

	/* ================================================================= */
	/* evaluate H- photodissociation rate, induced rec and rec cooling rates */
	/* >>chng 00 dec 24, add test so that photo rate only reevaluated two times per zone.
	 * in grain-free models this was sometimes dominated by Lya and so oscillated.  
	 * especially bad in primal.in - change 2 to 4 and primal.in will stop due to Lya oscil */
	/** \todo	2	following always true, why?  either remove test or use it - 
	 * it is here to save time - this step routine is called very often */
	/* >>chng 02 feb 16, add damper on H- photo rate, wild oscillations in Lya photo rate in 
	 * grain free models */

	t_phoHeat photoHeat;

	long ipIonPotHeI = ipoint(atmdat.EIonPot[ipHELIUM][0]);
	hmi.HMinus_photo_rate = GammaBn( hmi.iphmin, ipIonPotHeI, opac.iphmop, HMINUSIONPOT,
									 &hmi.HMinus_induc_rec_rate, &hmi.HMinus_induc_rec_cooling, &photoHeat );

	/* save H- photodissociation heating */
	hmi.HMinus_photo_heat = photoHeat.HeatNet;

	{
		/* following should be set true to print populations */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC)
		{
			fprintf(ioQQQ,"hminphoto\t%li\t%li\t%.2e\n", nzone, conv.nPres2Ioniz , hmi.HMinus_photo_rate );
		}
	}

	/* induced recombination */
	hmi.HMinus_induc_rec_rate *= hmi.rel_pop_LTE_Hmin*dense.eden;

	/* induced recombination cooling per unit volume */
	/** \todo	2	this should be done with new populations after converged soln */
	hmi.HMinus_induc_rec_cooling *= hmi.rel_pop_LTE_Hmin*dense.eden*findspecieslocal("H")->den;

	{
		/* following should be set true to debug H- photoionization rates */
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && nzone>400/*&& iteration > 1*/)
		{
			fprintf(ioQQQ,"hmoledebugg %.2f ",fnzone);
			GammaPrt(
				hmi.iphmin-1 , ipIonPotHeI , opac.iphmop ,
				/* io unit we will write to */
				ioQQQ, 
				/* total photo rate from previous call to GammaK */
				hmi.HMinus_photo_rate, 
				/* we will print contributors that are more than this rate */
				hmi.HMinus_photo_rate*0.05);
		}
	}
	/* add on high energy ionization, assume hydrogen cross section
	 * n.b.; HGAMNC with secondaries */
	/* >>chng 00 dec 24, above goes to HeI edge, no need for this, and was not important*/
	/*hmi.HMinus_photo_rate += iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].gamnc;*/

	/* ================================================================= */
	/* photodissociation by Lyman band absorption: esc prob treatment,
	 * treatment based on 
	 * >>refer	HI	abs	Tielens & Hollenbach 1985 ApJ 291, 722. */
	/* do up to carbon photo edge if carbon is turned on */
	/* >>>chng 00 apr 07, add test for whether element is turned on */
	hmi.UV_Cont_rel2_Habing_TH85_depth = 0.;
	/* >>chng 00 apr 07 from explicit ipHeavy to ipLo */
	/* find total intensity over carbon-ionizing continuum */
	/* >>chng 03 jun 09, use exact bounds rather than CI photo threshold for lower bound */
	/*for( i=ipLo; i < iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH1s].ipIsoLevNIonCon; i++ )*/
	/* the integral is from 6eV to 13.6, or 2060 - 912 Ang */
	for( i=rfield.ipG0_TH85_lo; i < rfield.ipG0_TH85_hi; ++i )
	{
		hmi.UV_Cont_rel2_Habing_TH85_depth += ((rfield.flux[0][i-1]) + (rfield.ConInterOut[i-1])+ 
															(rfield.outlin[0][i-1])+ (rfield.outlin_noplot[i-1]))*rfield.anu(i-1);
	}

	/* now convert to Habing ISM units
	 * UV_Cont_rel2_Habing_TH85_face is FUV continuum relative to Habing value 
	 * 1.6e-3 ergs em-2 s-1 is the Habing 1968 field, quoted on page 723, end of first
	 * full paragraph on left */
	hmi.UV_Cont_rel2_Habing_TH85_depth = (realnum)(hmi.UV_Cont_rel2_Habing_TH85_depth*EN1RYD/1.6e-3);
	/* if start of calculation remember G0 at illuminated face */
	if( nzone<=1 )
	{
		hmi.UV_Cont_rel2_Habing_TH85_face = hmi.UV_Cont_rel2_Habing_TH85_depth;
	}


	/* >>chng 05 jan 09, add special ver of G0 */
	hmi.UV_Cont_rel2_Habing_spec_depth = 0.; 
	for( i=rfield.ipG0_spec_lo; i < rfield.ipG0_spec_hi; ++i )
	{
		hmi.UV_Cont_rel2_Habing_spec_depth += ( rfield.flux[0][i-1] +  rfield.ConInterOut[i-1]+ 
															 rfield.outlin[0][i-1]+  rfield.outlin_noplot[i-1])*rfield.anu(i-1);
	}
	hmi.UV_Cont_rel2_Habing_spec_depth = (realnum)(hmi.UV_Cont_rel2_Habing_spec_depth*EN1RYD/1.6e-3);

	/* the Draine & Bertoldi version of the same thing, defined over their energy range */
	/* >>chng 04 feb 07, only evaluate at the illuminated face */
	if( !conv.nTotalIoniz )
	{
		hmi.UV_Cont_rel2_Draine_DB96_face = 0.f;
		/* this is sum of photon number between 912A and 1110, as per BD96 */
		for( i=rfield.ipG0_DB96_lo; i < rfield.ipG0_DB96_hi; ++i )
		{
			hmi.UV_Cont_rel2_Draine_DB96_face += ( rfield.flux[0][i-1] +  rfield.ConInterOut[i-1]+ 
																rfield.outlin[0][i-1]+  rfield.outlin_noplot[i-1]);
		}
		/* Convert into scaled ISM background field, total number of photons over value for 1 ISM field,
		 * the coefficient 1.232e7 is the number of photons over this wavelength range for 1H and is
		 * given in BD96, page 225, 4th line from start of section 4, also page 272, table 1, 2nd line
		 * from bottom */
		/* >>chng 04 mar 16, introduce the 1.71 */
		/* equation 20 of DB96 gives chi as flux over 1.21e7, to produce one Habing field.
		 * to get the Draine field we need to further divide by 1.71 as stated on the first
		 * line after equation 23 */
		hmi.UV_Cont_rel2_Draine_DB96_face = hmi.UV_Cont_rel2_Draine_DB96_face/(1.232e7f * 1.71f);
	}

	/* escape prob takes into account line shielding, 
	 * next is opacity then optical depth in H2 UV lines, using eqn A7 of TH85 */
	hmi.H2Opacity = (realnum)(1.2e-14*(1e5/GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN])));
	/* the typical Lyman -Werner H2 line optical depth eq A7 of TH85a */
	th2 = (findspecieslocal("H2")->column+ findspecieslocal("H2*")->column)*hmi.H2Opacity;
	/* the escape probability - chance that continuum photon will penetrate to
	 * this depth to pump the Lyman Werner bands */
	h2esc = esc_PRD_1side(th2,1e-4);

	/* cross section is eqn A8 of 
	 * >>refer	H2	dissoc	Tielens, A.G.G.M., & Hollenbach, D., 1985, ApJ, 291, 722
	 * branching ratio of 10% is already included, so 10x smaller than their number
	 * 10% lead to dissociation through H_2 + h nu => 2H */
	/* >>chng 05 mar 10, by TE, break into 2g and 2s 
	 * note use of same shielding column in below - can do far better */
	hmi.H2_Solomon_dissoc_rate_TH85_H2g = 3.4e-11 * hmi.UV_Cont_rel2_Habing_TH85_depth * h2esc;
	hmi.H2_Solomon_dissoc_rate_TH85_H2s = 3.4e-11 * hmi.UV_Cont_rel2_Habing_TH85_depth * h2esc;
	hmi.H2_H2g_to_H2s_rate_TH85 = hmi.H2_Solomon_dissoc_rate_TH85_H2g*9.;

	/* these are Burton et al. 1990 rates */
	hmi.H2_Solomon_dissoc_rate_BHT90_H2g = 3.4e-11 * hmi.UV_Cont_rel2_Habing_TH85_depth * h2esc;
	hmi.H2_Solomon_dissoc_rate_BHT90_H2s = 3.4e-11 * hmi.UV_Cont_rel2_Habing_TH85_depth * h2esc;
	hmi.H2_H2g_to_H2s_rate_BHT90 = hmi.H2_Solomon_dissoc_rate_BHT90_H2g*9.;

	{ 
		/* the following implements Drain & Bertoldi 1996, equation 37 from
		 * >>refer	H2	destruction	Draine, B.T., & Bertoldi, F., 1996, ApJ, 468, 269-289
		 * but the constant 4.6e-11 comes from Bertoldi & Draine equation 23,
		 * this is the normalized H2 column density */
		double x = (findspecieslocal("H2")->column+findspecieslocal("H2*")->column) / 5e14;
		double sqrtx = sqrt(1.+x);
		/* Doppler with of H2 in km/s */
		double b5 = GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN])/1e5;
		/* the molecular hydrogen line self-shielding factor */
		double fshield = 0.965/POW2(1.+x/b5) + 0.035/sqrtx *
			sexp(8.5e-4*sqrtx);

		/*double fshield = pow( MAX2(1.,colden.colden[ipCOLH2]/1e14) , -0.75 );*/
		/* this is the Bertoldi & Draine version of the Habing field,
		 * with dilution of radiation and extinction due to grains */
		/* >>chng 04 apr 18, moved fshield, the line shielding factor, from this defn to
		 * the following defn of dissociation rate, since following should measure continuum */
		hmi.UV_Cont_rel2_Draine_DB96_depth = hmi.UV_Cont_rel2_Draine_DB96_face * 
			(realnum)(sexp( opac.TauAbsFace[rfield.ip1000A-1] )/radius.r1r0sq);

		/* the following comes from Bertoldi & Draine 1996, equation 23,
		 * hmi.UV_Cont_rel2_Draine_DB96_depth already includes a factor of 1.71 to correct back from TH85 */
		/* >>chng 05 mar 10, TE, break into 2s and 2s */
		if( !mole_global.lgLeidenHack )
		{
			/* this is default, when set Leiden hack UMIST rates not entered */
			hmi.H2_Solomon_dissoc_rate_BD96_H2g = 4.6e-11 * hmi.UV_Cont_rel2_Draine_DB96_depth * fshield;
			hmi.H2_Solomon_dissoc_rate_BD96_H2s = 4.6e-11 * hmi.UV_Cont_rel2_Draine_DB96_depth * fshield;
		}
		else
		{
			/* done when set Leiden hack UMIST rates command entered */
			hmi.H2_Solomon_dissoc_rate_BD96_H2g = 5.18e-11* (hmi.UV_Cont_rel2_Habing_TH85_face/1.66f)
				*sexp(3.02*rfield.extin_mag_V_point)* fshield;
			hmi.H2_Solomon_dissoc_rate_BD96_H2s = 5.18e-11* (hmi.UV_Cont_rel2_Habing_TH85_face/1.66f)
				*sexp(3.02*rfield.extin_mag_V_point)* fshield;
		}

		/* BD do not give an excitation rate, so used 9 times the dissociation
		 * rate by analogy with 90% going back into X, as per TH85 */
		/*>>chng 04 feb 07, had used 90% relax into X from TH85,
		 * BD say that 15% dissociate, so 85/15 = 5.67 is factor */
		hmi.H2_H2g_to_H2s_rate_BD96 = 5.67* hmi.H2_Solomon_dissoc_rate_BD96_H2g;
	}

	/* do Elwert approximation to the dissociation rate */
	if( hmi.UV_Cont_rel2_Draine_DB96_face > SMALLFLOAT )
	{
		/* this evaluates the new H2 dissociation rate by Torsten Elwert */
		/* chng 05 jun 23, TE
		 * >>chng 05 sep 13, update master source with now approximation */

		/* Doppler with of H2 in km/s */
		double b5 = GetDopplerWidth(2.f*dense.AtomicWeight[ipHYDROGEN])/1e5;

		/* split the Solomon rates in H2g and H2s */
		/* >>chng 05 oct 19, 
		 * >>chng 05 dec 05, TE, define new approximation for the heating due to the destruction of H2
		 *	use this approximation for the specified cloud parameters, otherwise
		 * use limiting cases for 1 <= G0, G0 >= 1e7, n >= 1e7, n <= 1 */

		double x_H2g, x_H2s,
			fshield_H2g, fshield_H2s,
			f_H2s;
		static double a_H2g, a_H2s,
			e1_H2g, e1_H2s,
			e2_H2g,
			b_H2g,
			sl_H2g, sl_H2s,
			k_f_H2s,
			k_H2g_to_H2s, 
			log_G0_face = -1;

		/* define parameter range for the new approximation
		 * test for G0 
		 *BUGFIX - this tested on lg_G0_face < 0 for initialization needed so did not work
		 * in grids - change to evaluate in zone 0 */
		/* >>chng 07 feb 24, BUGFIX crash when G0=0 at start and radiation
		 * field builds up due to diffuse fields - soln is to always reevaluate */
		/*if( !nzone )*/
		{
			if(hmi.UV_Cont_rel2_Draine_DB96_face <= 1.) 
			{ 
				log_G0_face = 0.;
			}
			else if(hmi.UV_Cont_rel2_Draine_DB96_face >= 1e7) 
			{ 
				log_G0_face = 7.;
			}
			else 
			{ 
				log_G0_face = log10(hmi.UV_Cont_rel2_Draine_DB96_face); 
			}

			/* terms only dependent on G0_face */

			/* coefficients and exponents */
			a_H2g = 0.06 * log_G0_face + 1.32;
			a_H2s = 0.375 * log_G0_face + 2.125;

			e1_H2g = -0.05 * log_G0_face + 2.25;
			e1_H2s = -0.125 * log_G0_face + 2.625;

			e2_H2g = -0.005 * log_G0_face + 0.625;

			b_H2g = -4.0e-3  * log_G0_face + 3.2e-2;

			/* scalelength for H2g and H2s */
			sl_H2g = 4.0e14;
			sl_H2s = 9.0e15;

			/* coefficient for 2nd term of Solomon H2s */
			k_f_H2s = MAX2(0.1,2.375 * log_G0_face - 1.875 );

			/* coefficient for branching ratio */
			k_H2g_to_H2s =  MAX2(1.,-1.75 * log_G0_face + 11.25);

			/*fprintf( ioQQQ, "e1_H2g%.2e, e1_H2s%.2e, e2_H2g%.2e, b_H2g%.2e, a_H2g%.2e, a_H2s%.2e,sl_H2g: %.2e,sl_H2s: %.2e\n",
			  e1_H2g, e1_H2s, e2_H2g, b_H2g, a_H2g, a_H2s, sl_H2g, sl_H2s);
			*/
		}

		/* Solomon H2s ~G0^0.2 at large depth*/
		f_H2s = k_f_H2s * pow((double)hmi.UV_Cont_rel2_Draine_DB96_depth, 0.2 );

		/* scale length for absorption of UV lines */
		x_H2g = (findspecieslocal("H2")->column) / sl_H2g;
		x_H2s = (findspecieslocal("H2*")->column) / sl_H2s;

		/* the molecular hydrogen line self-shielding factor */
		fshield_H2g = 0.965/pow(1.+x_H2g/b5,e1_H2g) + b_H2g/pow(1.+x_H2g/b5,e2_H2g);
		fshield_H2s = 0.965/pow(1.+x_H2s/b5,e1_H2s);

		/* the Elwert Solomon rates for H2g and H2s  hmi.chH2_small_model_type == 'E' */
		hmi.H2_Solomon_dissoc_rate_ELWERT_H2g = a_H2g * 4.6e-11 * fshield_H2g * hmi.UV_Cont_rel2_Draine_DB96_depth;
		hmi.H2_Solomon_dissoc_rate_ELWERT_H2s = a_H2s * 4.6e-11 * fshield_H2s * (hmi.UV_Cont_rel2_Draine_DB96_depth + f_H2s);

		/* assume branching ratio dependent on G0*/
		hmi.H2_H2g_to_H2s_rate_ELWERT = k_H2g_to_H2s * hmi.H2_Solomon_dissoc_rate_ELWERT_H2g;

		/* use G0_BD96 as this definition declines faster with depth which is physical as
		 * the longer wavelengths in the definition of G0_TH85 cannot dissociate
		 * H2s directly */
		hmi.H2_photodissoc_ELWERT_H2s = hmi.UV_Cont_rel2_Draine_DB96_depth*1e-11;
		hmi.H2_photodissoc_ELWERT_H2g = hmi.H2_photodissoc_ELWERT_H2s * 1.0e-10;
	}
	else
	{
		hmi.H2_Solomon_dissoc_rate_ELWERT_H2g = 0.;
		hmi.H2_Solomon_dissoc_rate_ELWERT_H2s = 0.;
		hmi.H2_photodissoc_ELWERT_H2s = 0.;
		hmi.H2_photodissoc_ELWERT_H2g = 0.;
	}

	/* this is rate of photodissociation of H2*, A12 of TH85 */
	hmi.H2_photodissoc_TH85 = hmi.UV_Cont_rel2_Habing_TH85_depth*1e-11;

	/* dissociation rate from Burton et al. 1990 */
	hmi.H2_photodissoc_BHT90 = hmi.UV_Cont_rel2_Habing_TH85_depth*1e-11;

	/* rates for cosmic ray excitation of singlet bound electronic bound excited states 
	 * only add this to small molecule since automatically included in large 
	 *>>refer	H2	cr excit	Dalgarno, A., Yan, Min, & Liu, Weihong 1999, ApJS, 125, 237
	 * this is excitation of H2* */
	/* >>chng 05 sep 13, do not include this process when Leiden hacks are in place */
	cr_H2s = secondaries.x12tot*0.9 / 2. * hmi.lgLeidenCRHack;
	/* this is the fraction that dissociate */
	/* >>chng 05 sep 13, do not include this process when Leiden hacks are in place */
	cr_H2dis = secondaries.x12tot*0.1 / 2. * hmi.lgLeidenCRHack;

	/* >>chng 05 sep 13, TE, synchronize treatment of CR */
	/* cosmic ray rates for dissociation of ground and H2s 
	 * two factors done to agree with large H2 deep in the cloud where
	 * cosmic rays are important */
	cr_H2dis_ELWERT_H2g = secondaries.x12tot*5e-8 * hmi.lgLeidenCRHack; 
	cr_H2dis_ELWERT_H2s = secondaries.x12tot*4e-2 * hmi.lgLeidenCRHack;

	/* at this point there are two or three independent estimates of the H2 dissociation rate.
	 * if the large H2 molecule is on, then H2 Solomon rates has been defined in the last
	 * call to the large molecule.  Just above we have defined hmi.H2_Solomon_dissoc_rate_TH85,
	 * the dissociation rate from Tielens & Hollenbach 1985, and hmi.H2_Solomon_dissoc_rate_BD96,
	 * the rate from Bertoldi & Draine 1996.  We can use any defined rate.  If the big H2
	 * molecule is on, use its rate.  If not, for now use the TH85 rate, since that is the
	 * rate we always have used in the past.
	 * The actual rate we will use is given by hmi.H2_Solomon_dissoc_rate_used
	 */
	/* this is the Solomon process dissociation rate actually used */
	if( h2.lgEnabled  && h2.lgEvaluated && hmi.lgH2_Chemistry_BigH2 )
	{
		/* only update after big H2 molecule has been evaluated,
		 * when very little H2 and big molecule not used, leave at previous (probably TH85) value,
		 * since that value is always known */

		/* Solomon process rate from X into the X continuum with units s-1
		 * rates are total rate, and rates from H2g and H2s */ 
		hmi.H2_Solomon_dissoc_rate_used_H2g = h2.Solomon_dissoc_rate_g;
		ASSERT( hmi.H2_Solomon_dissoc_rate_used_H2g >= 0. );

		hmi.H2_Solomon_dissoc_rate_used_H2s = h2.Solomon_dissoc_rate_s;
		ASSERT( hmi.H2_Solomon_dissoc_rate_used_H2s >= 0. );

		/* photoexcitation from H2g to H2s */
		hmi.H2_H2g_to_H2s_rate_used = h2.gs_rate();
		ASSERT( hmi.H2_H2g_to_H2s_rate_used >= 0. );

		/* add up H2s + hnu (continuum) => 2H + KE, continuum photodissociation,
		 * this is not the Solomon process, true continuum, units s-1 */
		/* only include rates from H2s since this is only open channel, this process is well
		 * shielded against Lyman continuum destruction by atomic hydrogen */
		hmi.H2_photodissoc_used_H2s = h2.photodissoc_BigH2_H2s;
		/* NPA - 07/24/09 - logic to use Stancil photodissociation rate for H2s */
		if( mole_global.lgStancil && h2.lgEnabled )
			hmi.H2_photodissoc_used_H2s = h2.Cont_Dissoc_Rate_H2s;
		ASSERT( hmi.H2_photodissoc_used_H2s >= 0. );

		/* >>chng 05 mar 24, TE, continuum photodissociation rate of H2g, small correction factor accounts
		 * for unfavorable wavelength range of G0*/
		hmi.H2_photodissoc_used_H2g = h2.photodissoc_BigH2_H2g;
		/* NPA - 07/24/09 - logic to use Stancil photodissociation rate for H2g */
		if( mole_global.lgStancil && h2.lgEnabled )
			hmi.H2_photodissoc_used_H2g = h2.Cont_Dissoc_Rate_H2g;
		ASSERT( hmi.H2_photodissoc_used_H2g >= 0. );
	}
	else if( hmi.chH2_small_model_type == 'T' )
	{
		/* the TH85 rate  */
		/*>>chng 05 jun 23, add cosmic rays */
		hmi.H2_Solomon_dissoc_rate_used_H2g = hmi.H2_Solomon_dissoc_rate_TH85_H2g + cr_H2dis;
		/* >>chng 05 sep 13, cr_H2dis was not included */
		hmi.H2_Solomon_dissoc_rate_used_H2s = hmi.H2_Solomon_dissoc_rate_TH85_H2s + cr_H2dis;
		hmi.H2_H2g_to_H2s_rate_used = hmi.H2_H2g_to_H2s_rate_TH85 + cr_H2s;

		/* continuum photodissociation H2s + hnu => 2H, ,
		 * this is not the Solomon process, true continuum, units s-1 */
		hmi.H2_photodissoc_used_H2s = hmi.H2_photodissoc_TH85;
		/* >>chng 05 mar 24, TE, continuum photodissociation rate of H2g, small correction factor accounts
		 * for unfavorable wavelength range of G0*/
		hmi.H2_photodissoc_used_H2g = hmi.H2_photodissoc_TH85*1.0e-10f;
	}

	else if( hmi.chH2_small_model_type == 'H' )
	{
		/* the improved H2 formalism given by 
		 *>>refer	H2	dissoc	Burton, M.G., Hollenbach, D.J. & Tielens, A.G.G.M 
		 >>refcon	1990, ApJ, 365, 620 */
		hmi.H2_Solomon_dissoc_rate_used_H2g = hmi.H2_Solomon_dissoc_rate_BHT90_H2g + cr_H2dis;
		/* >>chng 05 sep 13, cr_H2dis was not included */
		hmi.H2_Solomon_dissoc_rate_used_H2s = hmi.H2_Solomon_dissoc_rate_BHT90_H2s + cr_H2dis;
		hmi.H2_H2g_to_H2s_rate_used = hmi.H2_H2g_to_H2s_rate_BHT90 + cr_H2s;

		/* continuum photodissociation H2s + hnu => 2H, ,
		 * this is not the Solomon process, true continuum, units s-1 */
		hmi.H2_photodissoc_used_H2s = hmi.H2_photodissoc_BHT90;
		/* >>chng 05 mar 24, TE, continuum photodissociation rate of H2g, small correction factor accounts
		 * for unfavorable wavelength range of G0*/
		hmi.H2_photodissoc_used_H2g = hmi.H2_photodissoc_BHT90*1.0e-10f;
	}

	else if( hmi.chH2_small_model_type == 'B' )
	{
		/* the Bertoldi & Draine rate - this is the default */
		/*>>chng 05 jun 23, add cosmic rays */
		hmi.H2_Solomon_dissoc_rate_used_H2g = hmi.H2_Solomon_dissoc_rate_BD96_H2g + cr_H2dis;
		/* >>chng 05 sep 13, cr_H2dis was not included */
		hmi.H2_Solomon_dissoc_rate_used_H2s = hmi.H2_Solomon_dissoc_rate_BD96_H2s + cr_H2dis;
		/* they did not do the excitation or dissoc rate, so use TH85 */
		hmi.H2_H2g_to_H2s_rate_used = hmi.H2_H2g_to_H2s_rate_BD96 + cr_H2s;


		/* continuum photodissociation H2s + hnu => 2H, ,
		 * this is not the Solomon process, true continuum, units s-1 */
		hmi.H2_photodissoc_used_H2s = hmi.H2_photodissoc_TH85;
		/* >>chng 05 mar 24, TE, continuum photodissociation rate of H2g, small correction factor accounts
		 * for unfavorable wavelength range of G0*/
		hmi.H2_photodissoc_used_H2g = hmi.H2_photodissoc_TH85*1.0e-10f;
	}
	else if( hmi.chH2_small_model_type == 'E' )
	{
		/* the Elwert et al. rate 
		 *>>chng 05 jun 23, add cosmic rays */
		hmi.H2_Solomon_dissoc_rate_used_H2g = hmi.H2_Solomon_dissoc_rate_ELWERT_H2g + cr_H2dis_ELWERT_H2g;
		hmi.H2_Solomon_dissoc_rate_used_H2s = hmi.H2_Solomon_dissoc_rate_ELWERT_H2s + cr_H2dis_ELWERT_H2s;
		hmi.H2_H2g_to_H2s_rate_used = hmi.H2_H2g_to_H2s_rate_ELWERT + cr_H2s;


		/* continuum photodissociation H2s + hnu => 2H, ,
		 * this is not the Solomon process, true continuum, units s-1 */
		hmi.H2_photodissoc_used_H2s = hmi.H2_photodissoc_ELWERT_H2s;
		hmi.H2_photodissoc_used_H2g = hmi.H2_photodissoc_ELWERT_H2g;
	}
	else
		TotalInsanity();

	{
		/*@-redef@*/
		enum {DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && h2.lgEnabled )
		{
			fprintf(ioQQQ," Solomon H2 dest rates: TH85 %.2e BD96 %.2e Big %.2e excit rates: TH85 %.2e Big %.2e\n",
					  hmi.H2_Solomon_dissoc_rate_TH85_H2g,
					  hmi.H2_Solomon_dissoc_rate_BD96_H2g,
					  h2.Solomon_dissoc_rate_g,
					  hmi.H2_H2g_to_H2s_rate_TH85 ,
					  h2.gs_rate() );
		}
	}

	if( !hd.lgEnabled )
		hd.photodissoc_BigH2_H2g = hmi.H2_photodissoc_used_H2g;

	return;
}
mole_reaction *mole_findrate_s(const char buf[])
{
	DEBUG_ENTRY( "mole_findrate_s()" );

	string newbuf = canonicalize_reaction_label(buf);

	mole_reaction_i p = mole_priv::reactab.find(newbuf);
		
	if(p != mole_priv::reactab.end())
		return &(*p->second);
	else
		return NULL;
}

double t_mole_local::findrk(const char buf[]) const
{
	DEBUG_ENTRY( "t_mole_local::findrk()" );

	mole_reaction *rate = mole_findrate_s(buf);

	if(!rate)
		return 0.;

	/* check for NaN */
	ASSERT( !isnan( reaction_rks[ rate->index ] ) );

	return reaction_rks[ rate->index ];
}
double t_mole_local::findrate(const char buf[]) const
{
	double ret;
	int i;

	DEBUG_ENTRY( "t_mole_local::findrate()" );

	mole_reaction *rate = mole_findrate_s(buf);
	if(!rate)
	{
		return 0.;
	}

	ret = reaction_rks[ rate->index ];
	for(i=0;i<rate->nreactants;i++)
		ret *= species[ rate->reactants[i]->index ].den;

	return ret;
}
/* Calculate rate at which molecular network abstracts species */

/* Need to check reactants vs reactant behaviour */
double t_mole_local::sink_rate_tot(const char chSpecies[]) const
{
	DEBUG_ENTRY( "t_mole_local::sink_rate_tot()" );

	const molecule* const sp = findspecies(chSpecies);
	double ratev = sink_rate_tot(sp);

	return ratev;
}
double t_mole_local::sink_rate_tot(const molecule* const sp) const
{
	DEBUG_ENTRY( "t_mole_local::sink_rate_tot()" );
	double ratev = 0;

	for(mole_reaction_i p=mole_priv::reactab.begin(); 
			p != mole_priv::reactab.end(); ++p) 
	{
		mole_reaction &rate = *p->second;
		ratev += sink_rate( sp, rate );
	}	

	return ratev;
}

double t_mole_local::sink_rate(const molecule* const sp, const char buf[]) const
{
	const mole_reaction* const rate = mole_findrate_s(buf);
	return sink_rate( sp, *rate );
}

double t_mole_local::sink_rate(const molecule* const sp, const mole_reaction& rate) const
{
	DEBUG_ENTRY( "t_mole_local::sink_rate()" );

	int ipthis = -1;
	for(int i=0;i<rate.nreactants && ipthis == -1;i++)
	{
		if(rate.reactants[i] == sp && rate.rvector[i]==NULL && rate.rvector_excit[i]==NULL ) 
		{
			ipthis = i;
		}
	}
	if(ipthis != -1) 
	{
		double ratevi = rate.a * rate.rk();
		for(int i=0;i<rate.nreactants;i++)
		{
			if(i!=ipthis)
			{
				ratevi *= species[ rate.reactants[i]->index ].den;
			}
		}
		return ratevi;
	}
	else
		return 0.;
}

/** returns the photodissociation rate per unit volume [cm^-3 s^-1] of species chSpecies */
double t_mole_local::dissoc_rate(const char chSpecies[]) const
{
	DEBUG_ENTRY( "t_mole_local::dissoc_rate()" );

	molecule *sp = findspecies(chSpecies);
	if (sp == null_mole)
		return 0.0;
	ASSERT(sp->isMonatomic());
	const chem_nuclide *tgt = sp->nNuclide.begin()->first.get();
	molecule *ph = findspecies("PHOTON");
	double ratev = 0.0;

	for (mole_reaction_i p
				 =mole_priv::reactab.begin(); p != mole_priv::reactab.end(); ++p) 
	{
		mole_reaction &rate = *p->second;
		
		// Must have a photon in to be a dissociation rate
		int ipph = 0;
		for (int i=0;i<rate.nreactants;i++)
		{
			if (rate.reactants[i] == ph) 
				ipph++;
		}
		if (!ipph)
			continue;

		// ipsp is number of *specific* species of interest, 
		// ipfree is number in same ionization ladder, including X-
		int ipspin = 0, ipfreein = 0;
		for (int i=0;i<rate.nreactants;i++)
		{
			if (rate.reactants[i] == sp)
				++ipspin;
			if (rate.reactants[i]->isMonatomic() && tgt == sp->nNuclide.begin()->first.get()) 
				++ipfreein;
		}
		int ipspout = 0, ipfreeout = 0;
		for (int i=0;i<rate.nproducts;i++)
		{
			if (rate.products[i] == sp) 
				++ipspout;
			if (rate.products[i]->isMonatomic() && tgt == sp->nNuclide.begin()->first.get()) 
				++ipfreeout;
		}

		// Must produce the species requested
		int newsp = ipspout-ipspin;
		if (newsp <= 0)
			continue;

		// And must do so by breaking bonds
		int nbondsbroken = ipfreeout-ipfreein;
		if (nbondsbroken <= 0)
			continue;
		// Fraction of the generated monatomic species which were *originally* bound
		double fracbroken = nbondsbroken/((double)ipfreeout);
		ASSERT( fracbroken <= 1.0 );
		
		double ratevi = reaction_rks[ rate.index ];
		for (int i=0;i<rate.nreactants;i++)
		{
			ratevi *= species[ rate.reactants[i]->index ].den;
		}

		// Photoproduction rate is rate of production of the species
		// which has not come from an initially monatomic source

		double ratesp = ratevi*newsp; // This is the total production
												// rate of the specific species
		ratesp *= fracbroken; // Scale back for any initially unbound
									 // monatoms

		ratev += ratesp;
	}	
	return ratev;
}
double t_mole_local::source_rate_tot(const char chSpecies[]) const
{
	DEBUG_ENTRY( "t_mole_local::source_rate_tot()" );

	molecule *sp = findspecies(chSpecies);
	double ratev = source_rate_tot(sp);

	return ratev;
}
double t_mole_local::source_rate_tot(const molecule* const sp) const
{
	DEBUG_ENTRY( "t_mole_local::source_rate_tot()" );
	double ratev = 0;

	for (mole_reaction_i p =mole_priv::reactab.begin(); p != mole_priv::reactab.end(); ++p) 
	{
		mole_reaction &rate = *p->second;
		int ipthis = 0;
		for(int i=0;i<rate.nproducts;i++)
		{
			if( rate.products[i] == sp && rate.pvector[i]==NULL && rate.pvector_excit[i]==NULL ) 
			{
				ipthis++;
			}
		}
		if(ipthis) 
		{
			double ratevi = rate.a * rate.rk();
			for(int i=0;i<rate.nreactants;i++)
			{
				ratevi *= species[ rate.reactants[i]->index ].den;
			}
			ratev += ipthis*ratevi;
		}
	}	

	return ratev;
}

double t_mole_local::chem_heat(void) const
{
	/* >>chng 07, Feb 11 NPA.  Calculate the chemical heating rate.  This is defined as the net energy of the
	 * reaction, which is:
	 *
	 * Energy = SUM[formation energies of reactants] - SUM[formation energies of products]
	 * 
	 * Now take the energy, and multiply by the densities of the reactants and the rate constant, finally
	 * you have to multiply by 1.66e-14, which is the conversion factor to go from kJ/mol to erg/atom 
	 * this gives the units in the form of erg/atom*cm3/s*cm-3*cm-3 = erg/cm-3/s/atom, which is  
	 * a heating rate 
	 */

	DEBUG_ENTRY( "t_mole_local::chem_heat()" );
	
	double heating = 0.;
	map<double,string> heatMap;
	molecule *ph = findspecies("PHOTON");
	molecule *crph = findspecies("CRPHOT");
	molecule *grn = findspecies("grn");

	/* loop over all reactions */
	for (mole_reaction_i p
				 =mole_priv::reactab.begin(); p != mole_priv::reactab.end(); ++p) 
	{
		mole_reaction &rate = *p->second;
	
		// If PHOTON appears, assume it accounts for energy difference
		bool lgCanSkip = false;
		for (int i=0;i<rate.nproducts;i++)
		{
			if( rate.products[i] == ph || rate.products[i] == crph ) 
				lgCanSkip = true;
		}
		for (int i=0;i<rate.nreactants;i++)
		{
			if( rate.reactants[i] == ph || rate.reactants[i] == crph ) 
				lgCanSkip = true;
		}
		// grain catalyst reactions are handled in grain physics, don't double count here
		for (int i=0;i<rate.nreactants;i++)
		{
			if( rate.reactants[i] == grn && rate.rvector[i] != NULL ) 
				lgCanSkip = true;
		}

		if( lgCanSkip )
			continue;

		/* This loop calculates the product of the rate constant times the densities*/
		double rate_tot = reaction_rks[ rate.index ];
		for( long i=0; i < rate.nreactants; ++i )
		{
			rate_tot *= species[ rate.reactants[i]->index ].den;
		}

		realnum reaction_enthalpy = 0.;

		/* Calculate the sum of the formation energies for the reactants */
		for( long i=0; i < rate.nreactants; ++i )
		{
			reaction_enthalpy += rate.reactants[i]->form_enthalpy;
		}

		/* Subtract from that the sum of the formation energies of the products */
		for( long i=0; i < rate.nproducts; ++i ) 
		{
			reaction_enthalpy -= rate.products[i]->form_enthalpy;
		}

		/* this is the chemical heating rate. */
	       	/** \todo  Once the H chem is merged with the C chem, then
		 * we will have the chemical heating rate for all reactions.
		 * This is only a subset and, thusfar, not actually used in
		 * getting the total heating.  Tests with pdr_leiden_hack_f1.in
		 * show that this heating rate can be up to 10% of the total
		 * heating */

		double heat = reaction_enthalpy*rate_tot*(1e10/AVOGADRO); /* 1.66e-14f; */
		heatMap[heat] = rate.label;
		heating += heat;
	}

	// use reverse iterator to print out biggest contributors
	long index = 0;
	// this should be a const_reverse_iterator, but pgCC 12.2-0 64-bit cannot handle this.
	// it appears there is no const version of heatMap.rend(); as a result the compiler
	// cannot find a suitable version of operator != in "it != heatMap.rend()"
	// The Solaris Studio compiler version 12.3 has the same problem
	for( map<double,string>::reverse_iterator it = heatMap.rbegin(); it != heatMap.rend(); ++it, ++index )
	{
		fprintf( ioQQQ, "DEBUGGG heat %li\t%li\t%.6e\t%s\n", index, nzone, it->first, it->second.c_str() );
		if( index==2 )
			break;
	}
	index = 0;
	for( map<double,string>::iterator it = heatMap.begin(); it != heatMap.end(); ++it, ++index )
	{
		if( it->first >= 0. )
			break;
		fprintf( ioQQQ, "DEBUGGG cool %li\t%li\t%.6e\t%s\n", index, nzone, it->first, it->second.c_str() );
		if( index==2 )
			break;
	}

	return heating;
}

STATIC double sticking_probability_H_func( double T_gas, double T_grain )
{
	DEBUG_ENTRY( "sticking_probability_H_func()" );
	double S = sticking_probability_H_HM79( T_gas, T_grain );
	return S;
}

STATIC double sticking_probability_H_HM79( double T_gas, double T_grain )
{
	DEBUG_ENTRY( "sticking_probability_H_HM79()" );
		
	/* sticking probability, 2H + grain equation 3.7 of
	 * >>refer	grain	phys	Hollenbach, D.J., & McKee, C.F., 1979, ApJS, 41, 555,
	 * fraction of H impacts on grain surface that stick */
	/* this sticking probability is used for both HM79 and CT02 */
	double T2 = T_gas / 100.;
	double Tg2 = T_grain / 100.;
	double S = 1./(1. + 0.4*sqrt(Tg2 + T2) + 0.2*T2 + 0.08*T2*T2);
	return S;
}

