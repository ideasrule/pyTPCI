/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MOLE_PRIV_H_
#define MOLE_PRIV_H_

#include "container_classes.h"

class molecule;
class mole_reaction;
class chem_element;
class chem_nuclide;

namespace mole_priv {
	extern map <string,shared_ptr<molecule> > spectab;
	extern map <string,shared_ptr<mole_reaction> > reactab;
	extern map <string,shared_ptr<chem_element> > elemtab;
	extern map <string,shared_ptr<mole_reaction> > functab;
}

class GroupMap {
public:
	multi_arr<double,2> fion;
	valarray<double> molElems;
	void updateMolecules(const valarray<double> &b2);
	void setup(double *b0vec);
	GroupMap (size_t size)
	{
		molElems.resize( size );
		fion.reserve(size);
		for( unsigned long i = 0; i < size; ++i )
			fion.reserve(i,size+1);
		fion.alloc();
	}
};

typedef map<string,shared_ptr<mole_reaction> >::iterator mole_reaction_i;
typedef map<string,shared_ptr<mole_reaction> >::const_iterator mole_reaction_ci;
typedef map<string,shared_ptr<molecule> >::iterator molecule_i;
typedef map<string,shared_ptr<chem_element> >::iterator chem_element_i;

extern vector<molecule *> groupspecies;

#define MAXREACTANTS 3
#define MAXPRODUCTS  4

/* Structure containing reaction data */
class mole_reaction {
 public:
	string label;
	int nreactants, nproducts, photon;
	molecule *reactants[MAXREACTANTS];
	molecule *rvector[MAXREACTANTS];
	molecule *rvector_excit[MAXREACTANTS];
	molecule *products[MAXPRODUCTS];
	molecule *pvector[MAXPRODUCTS];
	molecule *pvector_excit[MAXPRODUCTS];
	double reduced_mass, a, b, c;
	int udfastate;
	int source;
	long index;
	virtual double rk() const = 0;
	virtual mole_reaction* Create() const = 0;
	virtual const char *name() = 0;
	virtual ~mole_reaction() {};
};

enum udfastate {ABSENT, CORRECT, CONFLICT};

/** mole_eval_balance
*/
extern void mole_eval_balance(long int n, double *b, bool lgJac, multi_arr<double,2> &c);
/**mole_solve fills in matrix for heavy elements molecular routines 
*/
extern double mole_solve( void );

extern void mole_eval_sources(long int num_total);

extern realnum mole_return_cached_species(const GroupMap &MoleMap);

extern double frac_H2star_hminus();

/**mole_update_rks update rate coefficients, only temp part */
extern void mole_update_rks( void );

#endif /* MOLE_PRIV_H_ */

