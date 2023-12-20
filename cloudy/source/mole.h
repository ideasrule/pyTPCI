/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MOLE_H_
#define MOLE_H_

/* mole.h */

#include "module.h"
#include "container_classes.h"

class TransitionList;
class qList;
class species;
class molezone;

static const double SMALLABUND = 1e-24;

enum mole_state {MOLE_NULL, MOLE_PASSIVE, MOLE_ACTIVE};

class chem_nuclide;

class chem_element {
	explicit chem_element(); // Do not implement
	chem_element &operator=(const chem_element&); // Do not implement
public:
	explicit chem_element(int Z, const char*label) : Z(Z), label(label)
	{}
	~chem_element() throw()
	{}
	const int Z;
	const string label;
	map<int, shared_ptr<chem_nuclide> > isotopes;
	//(first -> Atomic A; second -> chem_nuclide )
	//(first -> -1 for bogus isotope, i.e. where no
	//	isotopes have been explicitly defined)
};

typedef map<int, shared_ptr<chem_nuclide> >::iterator isotopes_i;

class chem_nuclide {
public:
	// Link back to basic element for convenience -- not a shared_ptr
	// as this would lead to a reference cycle (and hence the
	// destructors never being called).  Many-to-one relation suggests
	// that the weak link should be this way around.
	chem_element* wel;
	int A;               /* mass number */
	unsigned long index; /* index in list for temporary use */
	vector<int> ipMl;    /* Atom and ion species in molecule arrays */
	realnum mass_amu;    /* mass of isotope in AMU */
	double frac;         /* fraction of element in this isotope */

	chem_element* el() const
	{
		return wel;
	}
	bool lgMeanAbundance( void ) const
	{
		return ( A < 0 );
	}

	// May have independent logic in the future...
	bool lgHasLinkedIon( void ) const
	{
		return lgMeanAbundance();
	}

	/* Chemical symbols for elements */
	string label( void ) const
	{
		if( lgMeanAbundance() )
			return el()->label;
		else if( el()->Z==1 && A==2 )
		{
			// Deuterium is a special case
			return "D";
		}
		else
		{
			char str[20];
			sprintf(str,"^%d",A);
			return ( str + el()->label );
		}
	}
	
	int compare ( const chem_nuclide &b ) const
	{
		// sort by proton number first
		if ( el()->Z < b.el()->Z )
			return -1;
		if ( el()->Z > b.el()->Z )
			return 1;

		if ( A < b.A )
			return -1;
		if ( A > b.A )
			return 1;

		return 0;
	}
};
inline bool operator< (const chem_nuclide &a, const chem_nuclide &b)
{
	return a.compare(b) < 0;
}
inline bool operator> (const chem_nuclide &a, const chem_nuclide &b)
{
	return a.compare(b) > 0;
}
inline bool operator<= (const chem_nuclide &a, const chem_nuclide &b)
{
	return a.compare(b) <= 0;
}
inline bool operator>= (const chem_nuclide &a, const chem_nuclide &b)
{
	return a.compare(b) >= 0;
}
inline bool operator== (const chem_nuclide &a, const chem_nuclide &b)
{
	return a.compare(b) == 0;
}
inline bool operator!= (const chem_nuclide &a, const chem_nuclide &b)
{
	return !(a == b);
}

typedef vector< shared_ptr<chem_nuclide> > ChemNuclideList;
extern ChemNuclideList nuclide_list;
typedef vector< shared_ptr<chem_element> > ChemElementList;
extern ChemElementList element_list;
extern chem_nuclide *null_nuclide;

class element_pointer_value_less
{
public:
	bool operator()(const shared_ptr<chem_nuclide>& a,
						const shared_ptr<chem_nuclide>& b) const
	{
		return *a < *b;
	}
};

/* Structure containing molecule data, initially only CO */
class molecule {
public:
	typedef map<const shared_ptr<chem_nuclide>, int,
		element_pointer_value_less> nNucsMap;

	string parentLabel;
	bool isIsotopicTotalSpecies() const
	{
		return parentLabel.empty();
	}
	int parentIndex;
	bool isEnabled;   /* Is it enabled? */

	/* Species physical data */
	string    label;        /** name */
	nNucsMap nNuclide;  	/** number of each nuclide in molecule */
	int     charge;         /** Charge on species/number of e- liberated by formation */
	bool    lgExcit;        /** Is species excited (e.g. H2*) */
	bool    lgGas_Phase;    /** Solid or gas phase? */
	int     n_nuclei(void) const       /** total number of nuclei */
	{
		int num = 0;
		for (nNucsMap::const_iterator el = nNuclide.begin(); 
			  el != nNuclide.end(); ++el)
		{
			num += el->second;
		}
		return num;
	}
	int nElement(int nelem)
	{
		int num = 0;
		for (nNucsMap::const_iterator el = nNuclide.begin(); 
			  el != nNuclide.end(); ++el)
		{
			if (el->first->el()->Z-1 == nelem)
				num += el->second;
		}
		return num;
	}

	/**isMonatomic -- Tell if a species is composed of only one atom
	 *
	 * \return 		true, if one nucleus in species; false, otherwise
	 */
	bool isMonatomic(void) const       /** total number of nuclei */
	{
		if (nNuclide.size() == 1 && nNuclide.begin()->second == 1)
			return true;
		return false;
	}

	/**isMolecule -- Tell if a species is a molecule, that is, if it is composed
	 * 		of more than one nuclei.
	 *
	 * \return 		true, if one nucleus in species; false, otherwise
	 */
	bool isMolecule() const
	{
		if( nNuclide.size() > 1 || ( nNuclide.size() == 1 && nNuclide.begin()->second > 1 ) )
			return true;
		return false;
	}

	realnum form_enthalpy;  /** formation enthalpy for the molecule (at 0K), in units of KJ/mol */
	realnum mole_mass;      /** Mass of molecule */

	/* Parameters as computational object */
	enum mole_state state;
	int index, groupnum;
	int n_react;		/** number of reactions that involve this molecule */

	class molezone *local(void) const;

	chem_nuclide *heavyAtom(void) //const
	{
		for( nNucsMap::reverse_iterator it=nNuclide.rbegin(); it!=nNuclide.rend(); ++it )
		{
			if (0 != it->second )
			{
				return it->first.get();
			}
		}
		return null_nuclide;
	}

	int compare(const molecule &mol2) const
	{
		nNucsMap::const_reverse_iterator it1, it2;

		for( it1 = nNuclide.rbegin(), it2 = mol2.nNuclide.rbegin();
			it1 != nNuclide.rend() && it2 != mol2.nNuclide.rend(); ++it1, ++it2 )
		{
			if( *(it1->first) > *(it2->first) )
				return 1;
			else if( *(it1->first) < *(it2->first) )
				return -1;
			else if( it1->second > it2->second)
				return 1;
			else if( it1->second < it2->second)
				return -1;
		}

		if( it1 != nNuclide.rend() && it2 == mol2.nNuclide.rend() )
			return 1;
		else if( it1 == nNuclide.rend() && it2 != mol2.nNuclide.rend() )
			return -1;
		else
			ASSERT( it1 == nNuclide.rend() && it2 == mol2.nNuclide.rend() );

		// sort by label if falls through to here       
		return ( label.compare(mol2.label) );
			
	}
}; 

/* iterators on nNuclide */	
typedef molecule::nNucsMap::iterator nNucs_i;
typedef molecule::nNucsMap::reverse_iterator nNucs_ri;
typedef molecule::nNucsMap::const_reverse_iterator nNucs_cri;

/**mole_drive main driver for chemical equilibrium routines */
extern void mole_drive(void);

/**mole_create_react build reaction structures */
extern void mole_create_react(void);

class mole_reaction;

mole_reaction *mole_findrate_s(const char buf[]);

extern molecule *null_mole;

inline bool exists (const molecule* m)
{
	return m != null_mole;
}

/** version for species we know are valid, but may not be present in calculation */
extern molecule *findspecies(const char buf[]);
/** version for user-supplied species which may not be valid */
extern molecule *findspecies_validate(const char buf[]);

/** version for species we know are valid, but may not be present in calculation */
extern molezone *findspecieslocal(const char buf[]);
/** version for user-supplied species which may not be valid */
extern molezone *findspecieslocal_validate(const char buf[]);

extern shared_ptr<chem_nuclide> findnuclide(const char buf[]);

/** \verbatim >>chng 03 feb 09, rm ipH3P_hev, since not used, and decrement NUM_HEAVY_MOLEC to 17 
 >>chng 03 aug 04, rm ipCTWO and ipC2P from den since not included in balance,
 and always finds zero column density, so NUM_HEAVY_MOLEC from 17 to 15 
 >>chng 03 aug 05, rm ch2 and ch3, so n from 15 to 13 
 >>chng 03 nov 14  add Si chemistry & CH3+, so that now every
     reaction that is in the TH85 chemical network is also included
     in Cloudy.  Additionally, there are also reactions taken from other
     papers (mostly Hollenbach and McKee...see co.c).  In all 20 molecular
     species are calculated, along with the atomic and first ionization 
	 stages of C, O, and Si
 >>chng 04 May 13, Nick Abel.  Add CH3, CH4, CH4+, and CH5+ to network in order 
	to get the same chemical abundances vs. depth as other PDR codes in the Leiden
	meeting.  With changes we now can predict molecular abundances for 24 C, O, 
	and Si bearing molecules. 

 >>chng 04 jul 13, Nick Abel.  Add nitrogen and sulphur bearing molecules
    to the chemical network.  First added to generate a chemical model for
    eta carinae, but is applicable to all molecular clouds 

 >>chng 05 mar 11, Nick Abel.  Add C2 and C2+ to chemistry, reactions 
    involving these species affects the abundance of C

 >>chng 05 mar 23, Nick Abel.  Add Chlorine to chemistry 
 \endverbatim */

/** this includes the atomic and first ionized species
   of each element that can combine to form molecules.  
   This is the number of molecules, ions, and atoms that the co network uses
   This is used in comole
   to improve the calculation, as deep in molecular regions reactions with molecules
   can be important to the ionization balance */

class Properties;


class t_mole_global : public module {

public:
	const char *chName() const
	{
		return "mole_global";
	}
	void zero();
	void comment(t_warnings&) {}
	/**mole_zero allocate + initialize workspace */
	void init(void);

	void make_species(void);
		
	/** flag to turn off all molecules, set with no molecules command */
	bool lgNoMole;

	/** flag to turn off heavy molecules, set with no heavy molecules command */
	bool lgNoHeavyMole;

	/** flag set true if H2O destruction rate went to zero */
	bool lgH2Ozer;

	/** option to turn on the original UMIST rates rather than the corrected network
	 * adopted in the Leiden PDR meeting, by deffault false,
	 * set true with the set UMIST rates command,
	 * these rates are from mole_co_umisthack.dat and also in mole_co_base.dat (UMIST)*/
	bool lgLeidenHack;

	/** option to use diffuse cloud chemical rates from Table 8 of
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is true - changed with set chemistry command
	 * these rates are in both mole_co_federman.dat and mole_co_base.dat,
	 * latter based on UMIST */
	bool lgFederman;

	bool lgStancil;

	/** option to use effective temperature as defined in
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	bool lgNonEquilChem;

	/** option to set proton elimination rates to zero
	 * >>refer	CO	chemistry	Huntress, W. T., 1977, ApJS, 33, 495
	 * By default, this is false - changed with set chemistry command */
	bool lgProtElim;

	/** option to not include neutrals in the non-equilibrium scheme
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	bool lgNeutrals;

	 /** do we include capture of molecules onto grain surfaces?  default is true,
	  * turned off with NO GRAIN MOLECULES command */
	bool lgGrain_mole_deplete;

	// flag saying whether to model isotopes (and isotopologues) of a given element
	vector<bool> lgTreatIsotopes;

	/** flag saying whether an element is in the chemistry network */
	int num_total, num_calc, num_compacted;

	typedef vector<shared_ptr<molecule> > MoleculeList;
	MoleculeList list;

	map<string,bool> offReactions;
	map<string,Properties > speciesProperties;

	static void sort(MoleculeList::iterator start,
				 MoleculeList::iterator end);
	static void sort(molecule **start, molecule **end);

};

extern t_mole_global mole_global;

class molezone {
public:
	molezone()
	{
		init();
	}
	void init (void)
	{
		location = NULL;
		levels = NULL;
		lines = NULL;
		dbase = NULL;
		zero();
	}
	void zero (void)
	{
		src = 0.;
		snk = 0.;
		den = 0.;
		column = 0.;
		atomLim = null_nuclide;
		xFracLim = 0.;
		column_old = 0.;
		index = -1;
	}
	int index;
	const molecule* global() const
	{
		if (index == -1)
			return null_mole;
		return mole_global.list[index].get();
	}
	
	double *location;      /** Location of density in non-molecule code, NULL if none exists */

	/** rate s-1 for molecular charge transfer, nelem from to */
	double src, snk;

	species* dbase;
	qList* levels;
	TransitionList* lines;

	/* Current zone data */
	double  den;       /** density (cm-3) */
	realnum column;    /** total column density in this iteration */
	chem_nuclide* atomLim; /* atom which is closest to limiting molecule density */
	realnum xFracLim;  /** The fraction of that element in this species */

	/* Historical solution data */
	realnum column_old;   /** total column density in previous iteration */
};

extern molezone *null_molezone;

class t_mole_local : public module
{
public:
	void alloc();
	void zero();
	void comment(t_warnings&) {}
	const char *chName() const
	{
		return "mole";
	}
	void set_ion_locations( );
	void set_isotope_abundances( void );
	double sink_rate_tot(const char chSpecies[]) const;
	double sink_rate_tot(const molecule* const sp) const;
	double sink_rate(const molecule* const sp, const mole_reaction& rate) const;
	double sink_rate(const molecule* const sp, const char buf[]) const;
	double source_rate_tot(const char chSpecies[]) const;
	double source_rate_tot(const molecule* const sp) const;
	/** returns the photodissociation rate per unit volume [cm^-3 s^-1] producing monatomic species chSpecies.
	 * *Excludes* photoionizations of other monatomic species, e.g. N-,PHOTON=>N,e- */
	double dissoc_rate(const char chSpecies[]) const;
	double chem_heat(void) const;
	double findrk(const char buf[]) const;
	double findrate(const char buf[]) const;

	double grain_area, grain_density, grain_saturation;

	/** total charge in molecules */
	double elec;

	/** these are source and sink terms for the ionization ladder, for chemical
	 * processes that remove and add species */
	multi_arr<double,2> source, sink;
	
	multi_arr<realnum,3> xMoleChTrRate;

	valarray<class molezone> species;

	vector<double> reaction_rks;
	vector<double> old_reaction_rks;
	long old_zone;
};

extern t_mole_local mole;

inline bool exists (const molezone* m)
{
	return m != null_molezone;
}

inline molezone *molecule::local(void) const
{
	if ( index != -1)
		return &mole.species[index];
	else
		return null_molezone;
}

extern void total_molecule_elems(realnum total[LIMELM]);
extern void total_molecule_deut(realnum &total);

extern realnum total_molecules(void);

extern realnum total_molecules_gasphase(void);

/**species_gasphase_density - Report gas phase density of requested species
 * \param chSpecies [in]	species to process
 *
 * \return density of species in the current zone
 */
realnum species_gasphase_density( const string &chSpecies );

extern void mole_make_list(void);
extern void mole_make_groups(void);

void mole_cmp_num_in_out_reactions(void);

bool lgDifferByExcitation( const molecule &mol1, const molecule &mol2 );

extern void mole_update_species_cache(void);

void mole_update_sources(void);

void mole_rk_bigchange(void);

void create_isotopologues(
	ChemNuclideList& atoms,
	vector< int >& numAtoms,
	string atom_old,
	string atom_new,
	string embellishments,
	vector<string>& newLabels );

void create_isotopologues_one_position(
	unsigned position,
	ChemNuclideList& atoms,
	vector< int >& numAtoms,
	string atom_old,
	string atom_new,
	string embellishments,
	string& newLabel );

bool parse_species_label( const char label[], ChemNuclideList &atomsLeftToRight, vector<int> &numAtoms, string &embellishments );
bool parse_species_label( const char mylab[], ChemNuclideList &atomsLeftToRight, vector<int> &numAtoms, string &embellishments,
	bool &lgExcit, int &charge, bool &lgGas_Phase );

/*HMRATE compile molecular rates using Hollenbach and McKee fits */
/* #define HMRATE(a,b,c) ( ((b) == 0 && (c) == 0) ? (a) : \
 *	( ((c) == 0) ? (a)*pow(phycon.te/300.,(b)) : \
 *	( ((c)/phycon.te > 50.) ? 0. : ( ((b) == 0) ?  (a)*exp(-(c)/phycon.te) : \
 *					 (a)*pow(phycon.te/300.,(b))*exp(-(c)/phycon.te) ) ) ) ) */

/* hmrate4 is clone of hmrate but takes a, b, c and te as explicit arguments */
inline double hmrate4( double a, double b, double c, double te )
{
	/* UMIST rates are simple temperaturr power laws that
	 * can become large at the high temperatures Cloudy
	 * may encounter. Do not extrapolate rates above T>2.5e3K 
	 * THIS CODE MUST BE KEPT PARALLEL WITH HMRATE IN MOLE_REACTIONS.CPP */ 
	if (b >0.)	
		te = min(te, 2500.);	
	if( b == 0. && c == 0. )
		return a;
	else if( c == 0. )
		return a*pow(te/300.,b);
	else if( b == 0. )
		return a*sexp(c/te);
	else
		return a*pow(te/300.,b)*sexp(c/te);
}

class Parser;
void ParseChemistry(Parser &p);

/**isSpecies - Tell if the input string corresponds to a valid & active species.
 *
 * \param chSpecies [in]	Species string
 *
 * \return		bool value: true, if a valid species; false, otherwise
 */
bool isSpecies( const string &chSpecies );

/**isMolecule - Tell if the input string corresponds to a valid & active molecule.
 *
 * \param chSpecies [in]	Species string
 *
 * \return		bool value: true, if a valid molecule; false, otherwise
 */
bool isMolecule( const string &chSpecies );

/**test_isMolecule() - Test function isMolecule()
 */
void test_isMolecule();

/**getMolecules - Get a list of all active molecules.
 *
 * \param allMolecules [in,out]		List of molecule strings.
 */
void getMolecules( vector<string> &allMolecules );

#endif /* MOLE_H_ */
