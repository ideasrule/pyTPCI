/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CO_Init called from cdInit to initialize co routines */
/*CO_update_chem_rates update rate coefficients, only temp part - in mole_co_etc.c */
#include "cdstd.h"
#include "cddefines.h"
#include "deuterium.h"
#include "elementnames.h"
#include "h2.h"
#include "hmi.h"
#include "dense.h"
#include "grainvar.h"
#include "trace.h"
#include "abund.h"
#include "mole_priv.h"
#include "mole.h"
#include "parse_species.h"
/*lint -e778 constant expression evaluates to 0 in operation '-' */

/* CO_update_chem_rates update rate coefficients, only temp part - in
 * mole_co_etc.c called in conv_base before any chemistry or
 * ionization is done */

enum spectype {MOLECULE, OTHER};

STATIC void read_species_file( string filename, bool lgCreateIsotopologues = true );
STATIC void newelement(const char label[], int ipion);
STATIC void newisotope( const shared_ptr<chem_element> &el, int massNumberA,
	realnum mass_amu, double frac );
STATIC realnum MeanMassOfElement( const shared_ptr<chem_element> &el );
// newspecies is overloaded.  The first one just calls the second with the additional lgCreateIsotopologues set true
STATIC molecule *newspecies(const char label[], enum spectype type,
	enum mole_state state, realnum form_enthalpy);
STATIC molecule *newspecies(const char label[], enum spectype type,
	enum mole_state state, realnum form_enthalpy, bool lgCreateIsotopologues);
STATIC bool isactive(const molecule &mol);
STATIC bool ispassive(const molecule &mol);
STATIC void SetIsotopeFractions( const vector<bool>& lgResolveNelem );

namespace mole_priv {
	map <string,shared_ptr<molecule> > spectab;
	map <string,shared_ptr<mole_reaction> > reactab;
	map <string,shared_ptr<chem_element> > elemtab;
	map <string,shared_ptr<mole_reaction> > functab;
}

map <string,size_t > nuclidetab;
typedef map<string,size_t >::iterator chem_nuclide_i;

molecule *null_mole;
molezone *null_molezone;
molezone null_molezone_impl;
chem_nuclide *null_nuclide;
chem_nuclide null_nuclide_impl;
ChemElementList element_list;
ChemNuclideList nuclide_list;
vector<molecule *> groupspecies;


namespace
{
	class MoleCmp : public binary_function<const shared_ptr<molecule>,
										   const shared_ptr<molecule>,bool>
	{
	public:
		bool operator()(const shared_ptr<molecule> &mol1, 
						const shared_ptr<molecule> &mol2) const
		{
			return mol1->compare(*mol2) < 0;
		}
		bool operator()(const molecule *mol1, const molecule *mol2) const
		{
			return mol1->compare(*mol2) < 0;
		}
	};
}

void t_mole_global::sort(t_mole_global::MoleculeList::iterator start, 
					 t_mole_global::MoleculeList::iterator end)
{
	std::sort(start,end,MoleCmp());	
}
void t_mole_global::sort(molecule **start, molecule **end)
{
	std::sort(start,end,MoleCmp());	
}

STATIC void SetIsotopeFractions( const vector<bool>& lgResolveNelem )
{
	DEBUG_ENTRY( "SetIsotopeFractions()" );

	for( int ielm = 0; ielm < LIMELM; ielm++ )
	{
		long Z = ielm;
		for( int j = 0; j < abund.IsoAbn[ielm].getNiso(); j++ )
		{
			long A = abund.IsoAbn[ielm].getAiso( j );
			realnum mass = abund.IsoAbn[ielm].getMass( A );
			realnum frac = abund.IsoAbn[ielm].getAbn( A );

			if( (unsigned)Z <= lgResolveNelem.size() && lgResolveNelem[Z] )
				newisotope( element_list[Z], A, mass, frac );
			// always do this to continue history of predicting 13CO
			else if( Z==ipCARBON )
				newisotope( element_list[Z], A, mass, frac );
		}
	}	

	return;
}

void t_mole_global::make_species(void)
{
	DEBUG_ENTRY( "mole_global::make_species()" );

	long int i;
	molecule *sp;

	null_nuclide = &null_nuclide_impl;
	null_molezone = &null_molezone_impl;
	null_molezone->den = 0.;

	/* set up concordance of elemental species to external Cloudy indices */
	for (i=0;i<LIMELM;i++) 
	{
		newelement(elementnames.chElementSym[i], i);
	}

	// flip this to treat deuterium
	if( deut.lgElmtOn )
		lgTreatIsotopes[ipHYDROGEN] = true;

	// read and define isotopes, set default fractionations
	SetIsotopeFractions( lgTreatIsotopes );

	if( lgTreatIsotopes[ipHYDROGEN] )
	{
		SetDeuteriumFractionation( element_list[ipHYDROGEN]->isotopes[2]->frac );
		SetGasPhaseDeuterium( dense.gas_phase[ipHYDROGEN] );
		InitDeuteriumIonization();
	}

	for( long nelem=0; nelem<LIMELM; ++nelem )
	{
		realnum mass_amu = MeanMassOfElement( element_list[nelem] );
		//define generic mean-abundance isotopes
		newisotope( element_list[nelem], -1, mass_amu, 1.0 );
	}
	
	/* set up properties of molecular species -- chemical formulae,
		 array indices, elementary components (parsed from formula), 
		 status within CO network, location of stored value external 
		 to CO network (must be floating point). */

	/* Sizes of different parts of network are calculated by increments in newspecies */
	mole_global.num_total = mole_global.num_calc = 0;
	/* Enthalpies of formation taken from
	 * >>refer    mol     Le Teuff, Y. E., Millar, T. J., & Markwick, A. J.,2000, A&AS, 146, 157
	 */
	
	/* Zero density pseudo-species to return when molecule is switched off */
	null_mole = newspecies("      ",OTHER,MOLE_NULL, 0.);
	null_mole->index = -1;

	read_species_file( "chem_species.dat" );

	if(gv.lgDustOn() && mole_global.lgGrain_mole_deplete )
	{		
		/* What are formation enthalpies of COgrn, H2Ogrn, OHgrn?  For 
			 present, take grn as standard state, and neglect adsorption enthalpy.

			 -- check e.g. http://www.arxiv.org/abs/astro-ph/0702322 for CO adsorption energy.
		*/
		if (0)
		{
			read_species_file( "chem_species_grn.dat" );
		}
		else
		{
			(void) newspecies("COgrn ",MOLECULE,MOLE_ACTIVE, -113.8f);
			(void) newspecies("H2Ogrn",MOLECULE,MOLE_ACTIVE, -238.9f);
			(void) newspecies("OHgrn ",MOLECULE,MOLE_ACTIVE, 38.4f);
			//sp = newspecies("Hgrn ",MOLECULE,MOLE_ACTIVE, 216.f);
		}
	}

	/* Add passive species to complete network */
	sp = newspecies("e-    ",OTHER,MOLE_PASSIVE, 0.0f);
	sp->charge = -1;	sp->mole_mass = (realnum)ELECTRON_MASS; /* Augment properties for this non-molecular species */
	(void) newspecies("grn   ",OTHER,MOLE_PASSIVE, 0.0f);
	(void) newspecies("PHOTON",OTHER,MOLE_PASSIVE, 0.0f);
	(void) newspecies("CRPHOT",OTHER,MOLE_PASSIVE, 0.0f);
	(void) newspecies("CRP   ",OTHER,MOLE_PASSIVE, 0.0f);

	if (!mole_global.lgLeidenHack)
		(void) newspecies("H-    ",MOLECULE,MOLE_ACTIVE, 143.2f);
	if (hmi.lgLeiden_Keep_ipMH2s) 
	{
		(void) newspecies("H2*   ",MOLECULE,MOLE_ACTIVE,
							 h2.ENERGY_H2_STAR * KJMOL1CM);
	}

	if( deut.lgElmtOn )
	{
		read_species_file( "chem_species_deuterium.dat", false );
	}

	/* Add species for all other elements and their first ions -- couple to network at least via H- */
	for( ChemNuclideList::iterator atom = nuclide_list.begin();
		  atom != nuclide_list.end(); ++atom) 
	{
		if ( !(*atom)->lgHasLinkedIon() )
			continue;
		unsigned long int nelem = (*atom)->el()->Z-1;

		for( unsigned long ion=0; ion<=nelem+1; ion++ )
		{
			char str[CHARS_ISOTOPE_SYM+CHARS_ION_STAGE+1];
			char temp[CHARS_ION_STAGE+1];
			if( ion==0 )
				temp[0] = '\0';
			else if( ion==1 )
				sprintf( temp, "+" );
			else if( ion < 100 )
				sprintf( temp, "+%ld", ion );
			else
				TotalInsanity();
			sprintf( str, "%.*s%s", CHARS_ISOTOPE_SYM, (*atom)->label().c_str(), temp );
			if (findspecies(str) == null_mole)
			{
				sp = newspecies(str,MOLECULE,MOLE_ACTIVE, 0.f);
				fixit("populate these in a local update");
#if 0
				if( sp != NULL )
				{
					sp->levels = NULL;
					sp->numLevels = 0;
				}
#endif
			}
		}
	}

	return;
}

void mole_make_list()
{
	DEBUG_ENTRY( "mole_make_list()" );

	/* Create linear list of species and populate it... */
	mole_global.list.resize(mole_global.num_total);

	/* ...first active species */
	long int i = 0;
	for (molecule_i p = mole_priv::spectab.begin(); p != mole_priv::spectab.end(); ++p) 
	{
		if (isactive(*(p->second)))
			mole_global.list[i++] = p->second;
	}
	ASSERT (i == mole_global.num_calc); 

	/* Sort list into a standard ordering */
	t_mole_global::sort(mole_global.list.begin(),mole_global.list.begin()+mole_global.num_calc);

	for (molecule_i p = mole_priv::spectab.begin(); p != mole_priv::spectab.end(); ++p) 
	{
		if (ispassive(*(p->second)))
			mole_global.list[i++] = p->second;
	}
	ASSERT (i == mole_global.num_total); 

	/* Set molecule indices to order of list just created */
	for(i=0;i<mole_global.num_total;i++) 
	{
		mole_global.list[i]->index = i;
	}

	for(i=0;i<mole_global.num_total;i++)
	{
		if( !mole_global.list[i]->parentLabel.empty() )
		{
			long parentIndex = findspecies( mole_global.list[i]->parentLabel.c_str() )->index;
			mole_global.list[i]->parentIndex = parentIndex;
		}
		else
			mole_global.list[i]->parentIndex = -1;
	}

	/* Register the atomic ladders */
	for(i=0;i<mole_global.num_total;i++) 
	{
		molecule *sp = &(*mole_global.list[i]);
		if (sp->isMonatomic())
		{
			ASSERT( (int)sp->nNuclide.size() == 1 );
			shared_ptr<chem_nuclide> atom = sp->nNuclide.begin()->first;
			ASSERT(sp->charge <= atom->el()->Z);
			if(sp->charge >= 0 && sp->lgGas_Phase) 
			{
				atom->ipMl[sp->charge] = i;
			}
		}
	}

	return;

}

STATIC void read_species_file( string filename, bool lgCreateIsotopologues )
{
	DEBUG_ENTRY( "read_species_file()" );
	
	fstream ioDATA;
	open_data( ioDATA, filename, mode_r );
	string line;

	while( getline( ioDATA,line ) )
	{
		if( line.empty() )
			break;
		if( line[0] == '#' )
			continue;
		istringstream iss( line );
		string species;
		double formation_enthalpy;
		iss >> species;
		iss >> formation_enthalpy;
		ASSERT( iss.eof() );
		newspecies( species.c_str(), MOLECULE,MOLE_ACTIVE, formation_enthalpy, lgCreateIsotopologues );
		//fprintf( ioQQQ, "DEBUGGG read_species_file %s\t%f\n", species.c_str(), formation_enthalpy );
	}

	return;
}
/*lint +e778 constant expression evaluates to 0 in operation '-' */

void create_isotopologues(
	ChemNuclideList& atoms,
	vector< int >& numNuclides,
	string atom_old,
	string atom_new,
	string embellishments,
	vector<string>& newLabels )
{
	DEBUG_ENTRY( "create_isotopologues()" );

	fixit("make sure atom_new and atom_old are isotopes");
	fixit("make sure atom_new is not already present");

	//for( ChemNuclideList::iterator it = atoms.begin(); it != atoms.end(); ++it )
	for( unsigned position = 0; position < atoms.size(); ++position )
	{
		string newLabel;
		create_isotopologues_one_position( position, atoms, numNuclides, atom_old, atom_new, embellishments, newLabel );
		if( !newLabel.empty() )
			newLabels.push_back( newLabel );
	}

	return;
}

void create_isotopologues_one_position(
	unsigned position,
	ChemNuclideList& atoms,
	vector< int >& numNuclides,
	string atom_old,
	string atom_new,
	string embellishments,
	string& newLabel )
{
	DEBUG_ENTRY( "create_isotopologues_one_position()" );

	stringstream str;
	if( atoms[position]->label() == atom_old )
	{
		for( unsigned i=0; i<position; ++i )
		{
			str << atoms[i]->label();
			if( numNuclides[i]>1 )
				str << numNuclides[i];
		}
			
		if( numNuclides[position] > 1 )
		{
			str << atom_old;
			if( numNuclides[position] > 2 )
				str << numNuclides[position]-1;
		}
			
		if( position+1 == atoms.size() )
			str << atom_new;

		for( unsigned i=position+1; i<atoms.size(); ++i )
		{
			if( i==position+1 )
			{
				// remove adjacent duplicates
				if( atom_new == atoms[i]->label() )
				{
					str << atoms[i]->label();
					ASSERT( numNuclides[i] + 1 > 1 );
					str << numNuclides[i] + 1;
				}
				else
				{
					str << atom_new;
					str << atoms[i]->label();
					if( numNuclides[i] > 1 )
						str << numNuclides[i];
				}
			}
			else
			{
				str << atoms[i]->label();
				if( numNuclides[i] > 1 )
					str << numNuclides[i];
			}
		}

		// add on charge, grn, and excitation embellishments
		str << embellishments;

		newLabel = str.str();
		//fprintf( ioQQQ, "DEBUGGG create_isotopologues_one_position %s\n", newLabel.c_str() );
	}

	return;
}

/* Fill element linking structure */
STATIC void newelement(const char label[], int ipion)
{
	DEBUG_ENTRY( "newelement()" );

	/* Create private workspace for label; copy and strip trailing whitespace */
	string mylab(label);
	auto p = mylab.find(' ');
	if( p != string::npos )
		mylab.erase(p);

	int exists = (mole_priv::elemtab.find(mylab) != mole_priv::elemtab.end());
	if (!exists)
	{
		shared_ptr<chem_element> element(new chem_element(ipion+1,mylab.c_str()));
		mole_priv::elemtab[element->label] = element;
		element_list.push_back(element);
	}
	return;
}

/* Fill isotope lists */
STATIC void newisotope( const shared_ptr<chem_element>& el, int massNumberA, 
	realnum mass_amu, double frac )
{

	DEBUG_ENTRY( "newisotope()" );

	ASSERT( massNumberA < 3 * LIMELM && ( massNumberA > 0 || massNumberA == -1 ) );
	ASSERT( mass_amu < 3. * LIMELM && mass_amu > 0. );
	ASSERT( frac <= 1. + FLT_EPSILON && frac >= 0. );

	shared_ptr<chem_nuclide> isotope(new chem_nuclide);
	isotope->A = massNumberA;
	isotope->mass_amu = mass_amu;
	isotope->frac = frac;
	isotope->ipMl.resize(el->Z+1);
	isotope->wel = el.get();
	for (long int ion = 0; ion < el->Z+1; ion++)
		isotope->ipMl[ion] = -1; 	/* Chemical network species indices not yet defined */

	//int exists = (mole_priv::nuclidetab.find( isotope->label() ) != mole_priv::nuclidetab.end());
	nuclidetab[ isotope->label() ] = nuclide_list.size();
	nuclide_list.push_back(isotope);
	// register with 'parent' element
	el->isotopes[massNumberA] = isotope;
}

STATIC realnum MeanMassOfElement( const shared_ptr<chem_element>& el )
{
	DEBUG_ENTRY( "MeanMassOfElement()" );

	// if no isotopes have been defined, just use the mean mass stored elsewhere
	if( el->isotopes.size()==0 )
		return dense.AtomicWeight[el->Z-1];

	realnum averageMass = 0., fracsum = 0.;
	for( isotopes_i it = el->isotopes.begin(); it != el->isotopes.end(); ++it )
	{
		fracsum += it->second->frac;
		averageMass += it->second->mass_amu * it->second->frac;
	}
	ASSERT( fp_equal( fracsum, 1_r ) );

	return averageMass;
}

STATIC molecule *newspecies(
	const char label[],
	enum spectype type,
	enum mole_state state,
	realnum form_enthalpy)
{
	return newspecies( label, type, state, form_enthalpy, true);
}

/* Parse species string to find constituent atoms, charge etc. */
STATIC molecule *newspecies(
	const char label[], 
	enum spectype type, 
	enum mole_state state, 
	realnum form_enthalpy,
	bool lgCreateIsotopologues )/* formation enthalpy at 0K */
{
	DEBUG_ENTRY( "newspecies()" );

	int exists;
	ChemNuclideList nuclidesLeftToRight;
	vector< int > numNuclides;
	string embellishments;
	shared_ptr<molecule> mol(new molecule);
	
	mol->parentLabel.clear();
	mol->isEnabled = true;
	mol->charge = 0;
	mol->lgExcit = false;
	mol->mole_mass = 0.0;
	mol->state = state;
	mol->lgGas_Phase = true;
	mol->form_enthalpy = form_enthalpy;
	mol->groupnum = -1;
	mol->n_react = 0;

	/* Create private workspace for label; copy and strip trailing whitespace */
	string mylab(label);
	auto p = mylab.find(' ');
	if( p != string::npos )
		mylab.erase(p);
	mol->label = mylab;

	if(type == MOLECULE)
	{
		if( parse_species_label( mylab.c_str(), nuclidesLeftToRight, numNuclides, embellishments, mol->lgExcit, mol->charge, mol->lgGas_Phase ) == false )
			return NULL;		
		
		for( unsigned i = 0; i < nuclidesLeftToRight.size(); ++i )
		{
			mol->nNuclide[ nuclidesLeftToRight[i] ] += numNuclides[i];
			mol->mole_mass += numNuclides[i] * nuclidesLeftToRight[i]->mass_amu * (realnum)ATOMIC_MASS_UNIT;
		}
	}
	
	// we also kill H- if molecules are disabled.  This is less than ideal,
	// but physically motivated by the fact that one of the strongest H- sinks
	// involves formation of H2 (disabled by "no molecules"), while one the 
	// fastest sources is e- recombination (which would still be allowed).
	if ( (mol->n_nuclei() > 1 || (mol->isMonatomic() && mol->charge==-1) ) && mole_global.lgNoMole) 
	{
		if( trace.lgTraceMole )
			fprintf(ioQQQ,"No species %s as molecules off\n",label);
		return NULL;
	}

	if (mol->n_nuclei() > 1 && mole_global.lgNoHeavyMole)
	{
		for( nNucs_ri it=mol->nNuclide.rbegin(); it != mol->nNuclide.rend(); ++it )
		{
			if( it->first->el()->Z-1 != ipHYDROGEN )
			{
				ASSERT( it->second > 0 );
				if( trace.lgTraceMole )
					fprintf(ioQQQ,"No species %s as heavy molecules off\n",label);
				return NULL;
			}
		}
	}

	if (speciesOff(mol->label))
	{
		fprintf( ioQQQ,"Species %s disabled\n",mol->label.c_str() );
		return NULL;
	}

	/* Insert species into hash table */
	exists = (mole_priv::spectab.find(mol->label) != mole_priv::spectab.end());
	if( exists )
	{
		fprintf( ioQQQ,"Warning: duplicate species %s - using first one\n", 
					mol->label.c_str() );
		return NULL;
	}
	else
		mole_priv::spectab[mol->label] = mol;

	// all map entries should have strictly positive number of nuclei
	for( nNucs_i it=mol->nNuclide.begin(); it != mol->nNuclide.end(); ++it )
		ASSERT( it->second > 0 );

	if (state != MOLE_NULL)
		mole_global.num_total++;
	if (state == MOLE_ACTIVE)
		mole_global.num_calc++;

	// this is a special case to always treat 13CO (as has long been done) 
	if( lgCreateIsotopologues && type==MOLECULE && mol->label.compare("CO")==0 ) 
	{
		molecule *sp = newspecies( "^13CO", MOLECULE, mol->state, mol->form_enthalpy, false );
		sp->parentLabel = mol->label;
	}

	// create singly-substituted isotopologues	
	if( lgCreateIsotopologues && type==MOLECULE && !mol->isMonatomic() )
	{
		for( nNucs_i it1 = mol->nNuclide.begin(); it1 != mol->nNuclide.end(); ++it1 )
		{
			for( map<int, shared_ptr<chem_nuclide> >::iterator it2 = it1->first->el()->isotopes.begin(); 
				 it2 != it1->first->el()->isotopes.end(); ++it2 )
			{
				// we don't want to create ^1H isotopologues (only substitute D for H)
				if( it1->first->el()->Z-1 == ipHYDROGEN && it2->second->A != 2 )
					continue;
				if( !mole_global.lgTreatIsotopes[it1->first->el()->Z-1] )
					continue;
				if( it2->second->lgMeanAbundance() )
					continue;
				vector<string> isoLabs;
				create_isotopologues( nuclidesLeftToRight, numNuclides, it1->first->label(), it2->second->label(), embellishments, isoLabs );
				//fprintf( ioQQQ, " DEBUGGG %10s isotopologues of %10s:", it1->first->label().c_str(), mol->label.c_str() );
				//for( vector<string>::iterator lab = isoLabs.begin(); lab != isoLabs.end(); ++ lab )
				//	fprintf( ioQQQ, " %10s", lab->c_str() );
				//fprintf( ioQQQ, "\n" );
				for( vector<string>::iterator newLabel = isoLabs.begin(); newLabel != isoLabs.end(); ++newLabel )
				{
					molecule *sp = newspecies( newLabel->c_str(), MOLECULE, mol->state, mol->form_enthalpy, false );
					// D is special -- don't set parentLabel
					if( sp!=NULL && it2->second->A != 2 )
						sp->parentLabel = mol->label;
				}
			}
		}
	}
		
	return &(*mol);
}
bool parse_species_label( const char label[], ChemNuclideList &nuclidesLeftToRight, vector<int> &numNuclides, string &embellishments ) 
{
	bool lgExcit, lgGas_Phase; 
	int charge;
	bool lgOK = parse_species_label( label, nuclidesLeftToRight, numNuclides, embellishments, lgExcit, charge, lgGas_Phase );
	return lgOK; 
}
bool parse_species_label( const char label[], ChemNuclideList &nuclidesLeftToRight, vector<int> &numNuclides, string &embellishments, 
	bool &lgExcit, int &charge, bool &lgGas_Phase )
{
	long int i, n, ipAtom;
	char thisAtom[CHARS_ISOTOPE_SYM];
	shared_ptr<chem_nuclide> atom;
	const char *s;
	int iend = strlen(label);
	
	/* Excitation... */
	s = strpbrk(label,"*");
	if(s)
	{
		lgExcit = true;
		embellishments = s;
		iend = s-label;
	} 

	/* ...Charge */
	s = strpbrk(label,"+-");
	if(s)
	{
		if(isdigit(*(s+1))) 
			n = atoi(s+1);
		else
			n = 1;
		if(*s == '+')
			charge = n;
		else
			charge = -n;
		embellishments = s + embellishments;
		iend = s-label;
	}
	/* ...Grain */
	s = strstr(label,"grn");
	if(s) 
	{
		lgGas_Phase = false;
		embellishments = s + embellishments;
		iend = s-label;
	} 
	else 
	{
		lgGas_Phase = true;
	}
	//fprintf( ioQQQ, "DEBUGGG parse_species_label %s\t%d\t%s\t%d\t%s\t\n", 
	//			label, (int)strlen(label), label, iend, embellishments.c_str() );

	/* Now analyse chemical formula */
	i = 0;
	while (i < iend && label[i] != ' ' && label[i] != '*') 
	{
		ipAtom = 0;
		/* Select next nuclide in species, matches regexp (^\d+)?[A-Z][a-z]? */
		/* look for isotope prefix */
		if(label[i]=='^')
		{
			thisAtom[ipAtom++] = label[i++];
			ASSERT( isdigit(label[i]) );
			thisAtom[ipAtom++] = label[i++];
			if(isdigit(label[i]))
			{
				thisAtom[ipAtom++] = label[i++];
			}
		}
		// should be first character of an element symbol
		thisAtom[ipAtom++] = label[i++];
		if( i<iend && islower(label[i])) 
		{
			thisAtom[ipAtom++] = label[i++];
		}
		thisAtom[ipAtom] = '\0';
		ASSERT(ipAtom <= CHARS_ISOTOPE_SYM);

		atom = findnuclide(thisAtom);
		if(atom.get() == NULL) 
		{
			fprintf(stderr,"Did not recognize atom at %s in \"%s \"[%ld]\n",thisAtom,label,i);
			exit(-1);
		}
		if(!dense.lgElmtOn[atom->el()->Z-1])
		{
			if( trace.lgTraceMole )
				fprintf(ioQQQ,"No species %s as element %s off\n",label,atom->el()->label.c_str() );
			return false;
		}
			
		if(i < iend && isdigit(label[i])) /* If there is >1 of this atom */
		{
			n = 0;
			do {
				n = 10*n+(long int)(label[i]-'0');
				i++;
			} while (i < iend && isdigit(label[i]));
		}
		else
		{
			n = 1;
		}
		nuclidesLeftToRight.push_back( atom );
		numNuclides.push_back( n );
	}

	return true;
}
STATIC bool isactive(const molecule &mol)
{
	DEBUG_ENTRY( "isactive()" );
	return mol.state == MOLE_ACTIVE;
}
STATIC bool ispassive(const molecule &mol)
{

	DEBUG_ENTRY( "ispassive()" );
	return mol.state == MOLE_PASSIVE;
}

bool lgDifferByExcitation( const molecule &mol1, const molecule &mol2 )
{
	if( mol1.label == mol2.label + "*" )
		return true;
	else if( mol2.label == mol1.label + "*" )
		return true;
	else 
		return false;
}

molecule *findspecies(const char buf[])
{
	DEBUG_ENTRY( "findspecies()" );

	// strip string of the first space and anything after it
	string s;
	for (const char *pb = buf; *pb && *pb != ' '; ++pb)
	{
		s += *pb;
	}

	const molecule_i p = mole_priv::spectab.find(s);

	if(p != mole_priv::spectab.end()) 
		return &(*p->second);
	else 
		return null_mole;
}

molecule *findspecies_validate(const char buf[])
{
	DEBUG_ENTRY( "findspecies_validate()" );

	// strip string of the first space and anything after it
	string s;
	for (const char *pb = buf; *pb && *pb != ' '; ++pb)
	{
		s += *pb;
	}

	const molecule_i p = mole_priv::spectab.find(s);

	if(p != mole_priv::spectab.end())
		return &(*p->second);
	else
	{
		fprintf(ioQQQ," PROBLEM specified species, '%s', is not valid or disabled.\n", s.c_str() );
		fprintf(ioQQQ," The SAVE SPECIES LABELS command will generate a list of species.  Species names are case sensitive.\n");
		cdEXIT( EXIT_FAILURE );
	}
}

molezone *findspecieslocal(const char buf[])
{
	DEBUG_ENTRY( "findspecieslocal()" );

	molecule *sp = findspecies(buf);

	if(sp != null_mole) 
		return &mole.species[ sp->index ];
	else 
		return null_molezone;
}

molezone *findspecieslocal_validate(const char buf[])
{
	DEBUG_ENTRY( "findspecieslocal_validate()" );

	molecule *sp = findspecies_validate(buf);

	return &mole.species[ sp->index ];
}

shared_ptr<chem_nuclide> findnuclide(const char buf[])
{
	chem_nuclide_i p;

	DEBUG_ENTRY( "findnuclide()" );

	p = nuclidetab.find(buf);

	if(p != nuclidetab.end())  
		return nuclide_list[p->second];
	else 
		return shared_ptr<chem_nuclide>(NULL);
}

void mole_update_species_cache(void)
{
	int i;
	double den_times_area, den_grains, adsorbed_density;

	DEBUG_ENTRY( "mole_update_species_cache()" );

	enum { DEBUG_MOLE = false };
 
	/* For rates that are dependent on grain physics.  This includes grain density, 
	cross sectional area, and dust temperature of each constituent.  Note that 
	
	gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3
	
	is the integrated projected grain surface area per cm^3 of gas for each grain size bin */

	/* >>chng 06 feb 28, turn off this rate when no grain molecules */
	/* >>chng 06 dec 05 rjrw: do this in newreact rather than rate */
	if( gv.lgDustOn() )
	{
		den_times_area = den_grains = 0.0;
		for( size_t nd=0; nd < gv.bin.size(); nd++ )
		{
			/* >>chng 06 mar 04, update expression for projected grain surface area, PvH */
			den_times_area += gv.bin[nd].IntArea/4.*gv.bin[nd].cnv_H_pCM3;
			den_grains += gv.bin[nd].cnv_GR_pCM3;
		}
		
		adsorbed_density = 0.0;
		for (i=0;i<mole_global.num_total;i++) 
		{
			if( !mole_global.list[i]->lgGas_Phase && mole_global.list[i]->isIsotopicTotalSpecies() ) 
				adsorbed_density += mole.species[i].den;
		}

		mole.grain_area = den_times_area;
		mole.grain_density = den_grains;

		double mole_cs = 1e-15;
		if (4*den_times_area <= mole_cs*adsorbed_density)
			mole.grain_saturation = 1.0;
		else
			mole.grain_saturation = ((mole_cs*adsorbed_density)/(4.*den_times_area));
	}
	else
	{
		mole.grain_area = 0.0;
		mole.grain_density = 0.0;
		mole.grain_saturation = 1.0;
	}

	if (DEBUG_MOLE)
		fprintf(ioQQQ,"Dust: %f %f %f\n",
						mole.grain_area,mole.grain_density,mole.grain_saturation);

	for (i=0;i<mole_global.num_total;i++) 
	{
		if(mole.species[i].location != NULL) 
		{
			ASSERT( mole_global.list[i]->isIsotopicTotalSpecies() );
			mole.species[i].den = *(mole.species[i].location);
			if (DEBUG_MOLE)
				fprintf(ioQQQ,"Updating %s: %15.8g\n",mole_global.list[i]->label.c_str(),mole.species[i].den);
		}
	}
	
	mole.set_isotope_abundances();
}

realnum mole_return_cached_species(const GroupMap &) 
// Finding the total atom density from MoleMap.molElems w
{
	int i;

	/* These two updates should together maintain the abundance invariant */

	// Assert invariant
	ASSERT(lgElemsConserved());

	// Update total of non-ladder species
	dense.updateXMolecules();
	deut.updateXMolecules();

	/* charge on molecules */
	mole.elec = 0.;
	for(i=0;i<mole_global.num_calc;i++)
	{
		if (mole.species[i].location == NULL && mole_global.list[i]->isIsotopicTotalSpecies())
			mole.elec += mole.species[i].den*mole_global.list[i]->charge;
	}

	// Update ionization ladder species
	realnum delta = 0.0;
	long ncpt = 0;

	for (i=0;i<mole_global.num_total;i++) 
	{
		if(mole.species[i].location && mole_global.list[i]->state == MOLE_ACTIVE) 
		{
			realnum new_pop = mole.species[i].den;

			if( mole_global.list[i]->isMonatomic() )
			{
				realnum old_pop = *(mole.species[i].location);
				long nelem = mole_global.list[i]->nNuclide.begin()->first->el()->Z-1;
				realnum frac = (new_pop-old_pop)/SDIV(new_pop+old_pop+1e-8*
																  dense.gas_phase[nelem]);
				delta += frac*frac;
				++ncpt;
			}

			//if( iteration >= 3 && nzone >= 100  )
			//	fprintf( ioQQQ, "DEBUGGG mole_return_ %i\t%s\t%.12e\t%.12e\t%.12e\t%.12e\t%li\n", 
			//	i, mole_global.list[i]->label.c_str(), new_pop, old_pop, frac, delta, ncpt );
			*(mole.species[i].location) = new_pop;
		}
	}

	// Assert invariant
	ASSERT(lgElemsConserved());
	return ncpt>0 ? sqrt(delta/ncpt) : 0.f;
}

void t_mole_local::set_isotope_abundances( void )
{
	DEBUG_ENTRY( "t_mole_local::set_isotope_abundances()" );

	// loop over unresolved elements
	for(ChemNuclideList::iterator atom = nuclide_list.begin(); atom != nuclide_list.end(); ++atom)
	{
		if ( !(*atom)->lgMeanAbundance() )
			continue;

		// loop over all isotopes of each element
		for( isotopes_i it = (*atom)->el()->isotopes.begin(); it != (*atom)->el()->isotopes.end(); ++it )
		{
			// skip mean-abundance "isotopes"
			if( it->second->lgMeanAbundance() )
				continue;
			
			// loop over all ions of the isotope
			for( unsigned long ion=0; ion<it->second->ipMl.size(); ++ion )
			{
				if( it->second->ipMl[ion] != -1 && 
					(species[ it->second->ipMl[ion] ].location == NULL ) && (*atom)->ipMl[ion] != -1 )
				{
					species[ it->second->ipMl[ion] ].den = species[ (*atom)->ipMl[ion] ].den * it->second->frac;
				}
			}
		}
	}
	
	return;
}

void t_mole_local::set_ion_locations()
{
	DEBUG_ENTRY( "t_mole_local::set_ion_locations()" );

	for( ChemNuclideList::iterator atom = nuclide_list.begin();
		  atom != nuclide_list.end(); ++atom) 
	{
		if ( !(*atom)->lgHasLinkedIon() )
			continue;
		long nelem = (*atom)->el()->Z-1;
		for( long ion=0; ion<nelem+2; ++ion )
		{
			long mole_index = (*atom)->ipMl[ion];
			// element not enabled if index is -1
			if( mole_index != -1 )
			{
				ASSERT( mole_index < mole_global.num_total );
				species[mole_index].location = &(dense.xIonDense[nelem][ion]);
				dense.ionMole[nelem][ion] = &species[mole_index];
			}
		}
	}

	if( deut.lgElmtOn )
	{
		findspecieslocal("D")->location = &(deut.xIonDense[0]);
		findspecieslocal("D+")->location = &(deut.xIonDense[1]);
	}	
}

void total_molecule_deut( realnum &total_f )
{
	DEBUG_ENTRY( "total_molecule_deut()" );
	
	double total = 0.;

	if( !deut.lgElmtOn )
		return;
	
	for (long int i=0;i<mole_global.num_calc;++i) 
	{
		if (mole.species[i].location == NULL && mole_global.list[i]->isIsotopicTotalSpecies() )
		{
			for( molecule::nNucsMap::iterator atom=mole_global.list[i]->nNuclide.begin();
				  atom != mole_global.list[i]->nNuclide.end(); ++atom)
			{
				long int nelem = atom->first->el()->Z-1;
				if( nelem==0 && atom->first->A==2 )
				{
					total += mole.species[i].den*atom->second;
				}
			}
		}
	}

	total_f = (realnum)total;

	return;
}
void total_molecule_elems(realnum total[LIMELM])
{
	DEBUG_ENTRY( "total_molecule_elems()" );

	/* now set total density of each element locked in molecular species */
	for ( long int nelem=ipHYDROGEN;nelem<LIMELM; ++nelem )
	{
		total[nelem] = 0.;
	}
	for (long int i=0;i<mole_global.num_calc;++i) 
	{
		if (mole.species[i].location == NULL && mole_global.list[i]->isIsotopicTotalSpecies() )
		{
			for( molecule::nNucsMap::iterator atom=mole_global.list[i]->nNuclide.begin();
				  atom != mole_global.list[i]->nNuclide.end(); ++atom)
			{
				ASSERT( atom->second > 0 );
				long int nelem = atom->first->el()->Z-1;
				if( atom->first->lgMeanAbundance() )
					total[nelem] += (realnum) mole.species[i].den*atom->second;
			}
		}
	}
}
realnum total_molecules(void)
{
	long int i;
	realnum total;

	DEBUG_ENTRY( "total_molecules()" );

	total = 0.;
	for (i=0;i<mole_global.num_calc;++i) 
	{
		if (mole.species[i].location == NULL && mole_global.list[i]->isIsotopicTotalSpecies())
			total += (realnum) mole.species[i].den;
	}
	return total;
}
realnum total_molecules_gasphase(void)
{
	long int i;
	realnum total;

	DEBUG_ENTRY( "total_molecules_gasphase()" );

	total = 0.;
	for (i=0;i<mole_global.num_calc;++i) 
	{
		if (mole_global.list[i]->lgGas_Phase && mole.species[i].location == NULL && mole_global.list[i]->isIsotopicTotalSpecies())
			total += (realnum) mole.species[i].den;
	}
	return total;
}
realnum species_gasphase_density( const string &chSpecies )
{
	DEBUG_ENTRY( "species_gasphase_density()" );

	const molecule *sp = findspecies( chSpecies.c_str() );

	if( sp != null_mole )
		return mole.species[ sp->index ].den;
	else
		return 0.;
}
/*lint +e778 const express eval to 0 */
/*lint +e725 expect positive indentation */

void mole_make_groups(void)
{
	long int i, j;
	/* Neutrals and positive ions will be treated as single species inside 
		 molecular equilibrium solver, to facilitate coupling with ionization
		 solvers */
	DEBUG_ENTRY ("mole_make_groups()" );
	if (mole_global.num_calc == 0)
	{
		groupspecies.resize(0);
		mole_global.num_compacted = 0;
		return;
	}
	groupspecies.resize(mole_global.num_calc);
	for (i=0,j=0;i<mole_global.num_calc;i++) 
	{
		if( mole_global.list[i]->isIsotopicTotalSpecies() && ( !mole_global.list[i]->isMonatomic() || mole_global.list[i]->charge <= 0 || ! mole_global.list[i]->lgGas_Phase ) )
		{
			/* Compound molecules and negative ions are represented individually */
			mole_global.list[i]->groupnum = j;
			groupspecies[j++] = &(*mole_global.list[i]);
		}
		else
		{
			/* All positive ions are collapsed into single macrostate (i.e. H+ -> H) */
			/* Need to increase constant if higher ions are included */
			ASSERT (mole_global.list[i]->charge < LIMELM+1);
			ASSERT (mole_global.list[i]->groupnum == -1);
		}
	}
	mole_global.num_compacted = j;
	groupspecies.resize(mole_global.num_compacted);

	for (i=0;i<mole_global.num_calc;i++) 
	{
		if (mole_global.list[i]->groupnum == -1)
		{
			if( mole_global.list[i]->isMonatomic() && mole_global.list[i]->isIsotopicTotalSpecies() )
			{
				for( nNucs_i it = mole_global.list[i]->nNuclide.begin(); it != mole_global.list[i]->nNuclide.end(); ++it )
				{
					ASSERT( it->second > 0 );
					mole_global.list[i]->groupnum = mole_global.list[it->first->ipMl[0]]->groupnum;
					break;
				}
			}
			else
			{
				ASSERT( !mole_global.list[i]->isIsotopicTotalSpecies() );
				mole_global.list[i]->groupnum = mole_global.list[ mole_global.list[i]->parentIndex ]->groupnum;
			}
		}
	
		ASSERT( mole_global.list[i]->groupnum != -1);
	}

	return;
}

bool isSpecies( const string &chSpecies )
{
	DEBUG_ENTRY ("isSpecies()" );

	const molecule *this_mole = findspecies( chSpecies.c_str() );

	{
		enum { DEBUG_LOC = false };
		if( DEBUG_LOC )
		{
			printf( "chSpecies= '%s'\t", chSpecies.c_str() );
			if( this_mole != null_mole )
				printf( "label= '%s'\n", this_mole->label.c_str() );
			else
				printf( "NULL\n" );
		}
	}

	if( this_mole == null_mole )
		return false;
	if( this_mole->isMonatomic() || this_mole->isMolecule() )
		return true;
	return false;
}

bool isMolecule( const string &chSpecies )
{
	DEBUG_ENTRY ("isMolecule()" );

	const molecule *this_mole = findspecies( chSpecies.c_str() );

	{
		enum { DEBUG_LOC = false };
		if( DEBUG_LOC )
		{
			printf( "chSpecies= '%s'\t", chSpecies.c_str() );
			if( this_mole != null_mole )
				printf( "label= '%s'\n", this_mole->label.c_str() );
			else
				printf( "NULL\n" );
		}
	}

	if( this_mole == null_mole )
		return false;
	if( this_mole->isMolecule() )
		return true;
	return false;
}

void test_isMolecule()
{
	DEBUG_ENTRY( "test_isMolecule()" );

	for( size_t iSpec = 0; iSpec < mole_global.list.size(); iSpec++ )
	{
		fprintf( ioQQQ,
			"'%s'\t%s\n",
			mole_global.list[ iSpec ]->label.c_str(),
			isMolecule( mole_global.list[ iSpec ]->label ) ? "TRUE" : "FALSE" );
	}

	cdEXIT( EXIT_SUCCESS );
}

void getMolecules( vector<string> &allMolecules )
{
	DEBUG_ENTRY( "getMolecules()" );

	for( size_t iSpec = 0; iSpec < mole_global.list.size(); iSpec++ )
	{
		if( isMolecule( mole_global.list[ iSpec ]->label ) )
			allMolecules.push_back( mole_global.list[ iSpec ]->label );
	}
}
