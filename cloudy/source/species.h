/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef SPECIES_H_
#define SPECIES_H_

/*The species structure is used to hold information about a particular atom,ion or molecule
mentioned in the species.ini file.The name of the atom/ion/molecule is used to obtain the density
of molecules in the case of the Leiden Database and along with atomic number and ion stage the 
density of atoms/ions in the case of the CHIANTI database */
class species
{
public:
	species()
	{
		chLabel = NULL;
		database = "";
		dataset = "";
		index = INT_MAX;
		numLevels_max = -INT_MAX;
		setLevels = -1;
		numLevels_local = -INT_MAX;
		numLevels_masterlist = -1;
		fmolweight = realnum(0.);
		lgMolecular = false;
		fracType = realnum(0.);
		fracIsotopologue = realnum(0.);
		CoolTotal = 0.;
		lgActive = false;
		maxWN = 0.;
		lgLTE = false;
		lgPrtMatrix = false;
		lgImgMatrix = false;
	}
	~species()
	{
		delete[] chLabel;
	}
	/*Name of the atom/ion/ molecule*/
	char *chLabel;
	/* Database origin */
	string database;
	/* nickname of alternate dataset, empty string for default set */
	string dataset;
	// index in chemistry
	long index;
	/*Actual Number of energy levels in the data file*/
	long numLevels_max;
	/*Number of levels defined in "species" command */
	long setLevels;
	/*Number of energy levels used locally*/
	long numLevels_local;
	/** Masterlist defined minimum number of levels, needed for
	 *  species like NI which require significant number of
	 *  levels to get pumped contributions to important lines */
	long numLevels_masterlist;
	/*Molecular weight*/
	realnum fmolweight;
	/* is molecular? */
	bool lgMolecular;
	/* fraction in this "type" (e.g. para, ortho) */
	realnum fracType;
	/* chemical fractionation */
	realnum fracIsotopologue;
	/** total cooling due to this species */
	double CoolTotal;
	/** are populations currently defined? */
	bool lgActive;
	/* maximum wavenumber in chianti */
	double maxWN;
	/** should level pops be done in LTE? */
	bool lgLTE;
	/** print Matrix input to solver */
	bool lgPrtMatrix;
	/** create image for Matrix input to solver */
	bool lgImgMatrix;
};

struct t_pseudo_cont
{
	realnum wlLo;
	realnum wlHi;
	long nBins;
};

extern t_pseudo_cont pseudoContDef;

/**
 * db_basename_to_spectral - Convert the species basename in Stout or Chianti
 * 				to a spectral label.
 * \param chBasename[in]	File basename, e.g., "c_2"
 * \param chElement[out]	Extracted element name, e.g., "C"
 * \param chSpecLabel[out]	Spectral label, e.g., "C  2"
 */
void db_basename_to_spectral( const string& chBasename, string &chElement,
                                char *chSpecLabel );

void parsespect(const char* chLabel, long& nelem, long& IonStg);
string makeChemical(long nelem, long ion);
void makeChemical(char* chLabelChemical, long nelem, long ion);

/**
 * spectral_to_chemical - Convert a spectral label to chemical.
 *
 * \param chLabelSpec[in]       Spectral label, e.g., "C  2"
 * \param chLabelChem[out]      Chemical label, e.g., "C+"
 */
void spectral_to_chemical( char *chLabelSpec, const char* chLabelChem );

/**
 * spectral_to_chemical - Convert a spectral label to chemical.
 *
 * \param chLabelSpec[in]       Spectral label, e.g., "C  2"
 * \param chLabelChem[out]      Chemical label, e.g., "C+"
 */
void spectral_to_chemical( string &chLabelSpec, const char* chLabelChem );

/**
 * parse_chemical - Parse the element and charge of a chemical label.
 *
 * \param chLabelChem[in]       Chemical label, e.g., "C+2", or "HCl"
 * \param chElemMol[out]        Element or molecule, e.g., "C", or "HCl"
 * \param spCharge[out]         Species charge, e.g., 2, or 0
 *
 * \return bool			True, if successful; false, otherwise.
 */
bool parse_chemical( const string &chLabelChem,
			string &chElemMol, long &spCharge );

/**
 * chemical_to_spectral - Convert a chemical label to spectral.
 *
 * \param chLabelChem[in]       Chemical label, e.g., "C+"
 * \param chLabelSpec[out]      Spectral label, e.g., "C  2"
 */
void chemical_to_spectral( const string &chLabelChem, string &chLabelSpec );

/**
 * isSpeciesActive - Tell if the species is enabled
 *
 * An atomic/ionic species may be inactive if the relevant element is disabled.
 * Likewise, a molecular species may be inactive if one of its constituents
 * elements is disabled.  This function reports on the given species.
 * If not found, the result is false.
 *
 * \param chLabelChem [in]	Chemical label, e.g., "C+2", or "HCl"
 *
 * \return bool			True, if species active; false, otherwise.
 */
bool isSpeciesActive( const string &chSpecLabel );

/** isAtomicIonValid - Tell if an atomic ion is has a charge consistent
 * 			with its atomic number.
 *
 * \param element_index[in]	Zero-offset index for element; 0 for H, 29 for Zn
 * \param spCharge[in]		Species charge
 *
 * \return bool			True, if ion is valid; false, otherwise.
 */
bool isAtomicIonValid( const long element_index, const long spCharge );

/** isBareNucleus - Tell if the given ion is a bare nucleus.
 *
 * \param element_index[in]	Zero-offset index for element; 0 for H, 29 for Zn
 * \param spCharge[in]		Species charge
 *
 * \return bool			True, if ion is bare nucleus; false, otherwise.
 */
bool isBareNucleus( const long element_index, const long spCharge );

#endif /* SPECIES_ */
