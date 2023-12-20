/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ELEMENTNAMES_H_
#define ELEMENTNAMES_H_

enum {CHARS_ELEMENT_NAME=11, CHARS_ELEMENT_NAME_SHORT=5, CHARS_ELEMENT_SYM=3,
		CHARS_ION_STAGE=3, CHARS_ION_ROMAN=7};

/**\file elementnames.h */
/** set of names of the chemical elements, long and short forms */
struct t_elementnames {

	/** following used for prints in each zone, full name.
	 * the LIMELEM element is 12CO, 
	 * +1 is 13CO, +2 is H2 */
	char chElementName[LIMELM][CHARS_ELEMENT_NAME];

	/** labels for match on element name
	 * this must be caps for present logic in matches */
	char chElementNameShort[LIMELM][CHARS_ELEMENT_NAME_SHORT];

	/** two letter very short form of element name, used to make
	 * emission line labels */
	char chElementSym[LIMELM][CHARS_ELEMENT_SYM];

	/** this is series of two char numbers, beginning with " 1" and
	 * ending with "31" */
	char chIonStage[LIMELM+1][CHARS_ION_STAGE];

	/** string giving ionization stage as roman numerals */
	char chIonRoman[LIMELM+1][CHARS_ION_ROMAN];

	/* here lies the first C++ code, written by PvH, 2006 Nov 20, after
	 * the conversion from C to C++.  It is a constructor. */
	t_elementnames() {
		strncpy( chElementName[0],  "Hydrogen  ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[1],  "Helium    ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[2],  "Lithium   ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[3],  "Beryllium ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[4],  "Boron     ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[5],  "Carbon    ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[6],  "Nitrogen  ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[7],  "Oxygen    ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[8],  "Fluorine  ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[9],  "Neon      ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[10], "Sodium    ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[11], "Magnesium ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[12], "Aluminium ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[13], "Silicon   ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[14], "Phosphorus", CHARS_ELEMENT_NAME );
		strncpy( chElementName[15], "Sulphur   ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[16], "Chlorine  ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[17], "Argon     ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[18], "Potassium ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[19], "Calcium   ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[20], "Scandium  ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[21], "Titanium  ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[22], "Vanadium  ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[23], "Chromium  ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[24], "Manganese ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[25], "Iron      ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[26], "Cobalt    ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[27], "Nickel    ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[28], "Copper    ", CHARS_ELEMENT_NAME );
		strncpy( chElementName[29], "Zinc      ", CHARS_ELEMENT_NAME );

		strncpy( chElementNameShort[0],  "HYDR", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[1],  "HELI", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[2],  "LITH", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[3],  "BERY", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[4],  "BORO", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[5],  "CARB", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[6],  "NITR", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[7],  "OXYG", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[8],  "FLUO", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[9],  "NEON", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[10], "SODI", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[11], "MAGN", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[12], "ALUM", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[13], "SILI", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[14], "PHOS", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[15], "SULP", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[16], "CHLO", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[17], "ARGO", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[18], "POTA", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[19], "CALC", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[20], "SCAN", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[21], "TITA", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[22], "VANA", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[23], "CHRO", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[24], "MANG", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[25], "IRON", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[26], "COBA", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[27], "NICK", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[28], "COPP", CHARS_ELEMENT_NAME_SHORT );
		strncpy( chElementNameShort[29], "ZINC", CHARS_ELEMENT_NAME_SHORT );

		strncpy( chElementSym[0],  "H ", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[1],  "He", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[2],  "Li", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[3],  "Be", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[4],  "B ", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[5],  "C ", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[6],  "N ", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[7],  "O ", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[8],  "F ", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[9],  "Ne", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[10], "Na", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[11], "Mg", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[12], "Al", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[13], "Si", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[14], "P ", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[15], "S ", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[16], "Cl", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[17], "Ar", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[18], "K ", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[19], "Ca", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[20], "Sc", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[21], "Ti", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[22], "V ", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[23], "Cr", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[24], "Mn", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[25], "Fe", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[26], "Co", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[27], "Ni", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[28], "Cu", CHARS_ELEMENT_SYM );
		strncpy( chElementSym[29], "Zn", CHARS_ELEMENT_SYM );

		strncpy( chIonStage[0],  " 1", CHARS_ION_STAGE );
		strncpy( chIonStage[1],  " 2", CHARS_ION_STAGE );
		strncpy( chIonStage[2],  " 3", CHARS_ION_STAGE );
		strncpy( chIonStage[3],  " 4", CHARS_ION_STAGE );
		strncpy( chIonStage[4],  " 5", CHARS_ION_STAGE );
		strncpy( chIonStage[5],  " 6", CHARS_ION_STAGE );
		strncpy( chIonStage[6],  " 7", CHARS_ION_STAGE );
		strncpy( chIonStage[7],  " 8", CHARS_ION_STAGE );
		strncpy( chIonStage[8],  " 9", CHARS_ION_STAGE );
		strncpy( chIonStage[9],  "10", CHARS_ION_STAGE );
		strncpy( chIonStage[10], "11", CHARS_ION_STAGE );
		strncpy( chIonStage[11], "12", CHARS_ION_STAGE );
		strncpy( chIonStage[12], "13", CHARS_ION_STAGE );
		strncpy( chIonStage[13], "14", CHARS_ION_STAGE );
		strncpy( chIonStage[14], "15", CHARS_ION_STAGE );
		strncpy( chIonStage[15], "16", CHARS_ION_STAGE );
		strncpy( chIonStage[16], "17", CHARS_ION_STAGE );
		strncpy( chIonStage[17], "18", CHARS_ION_STAGE );
		strncpy( chIonStage[18], "19", CHARS_ION_STAGE );
		strncpy( chIonStage[19], "20", CHARS_ION_STAGE );
		strncpy( chIonStage[20], "21", CHARS_ION_STAGE );
		strncpy( chIonStage[21], "22", CHARS_ION_STAGE );
		strncpy( chIonStage[22], "23", CHARS_ION_STAGE );
		strncpy( chIonStage[23], "24", CHARS_ION_STAGE );
		strncpy( chIonStage[24], "25", CHARS_ION_STAGE );
		strncpy( chIonStage[25], "26", CHARS_ION_STAGE );
		strncpy( chIonStage[26], "27", CHARS_ION_STAGE );
		strncpy( chIonStage[27], "28", CHARS_ION_STAGE );
		strncpy( chIonStage[28], "29", CHARS_ION_STAGE );
		strncpy( chIonStage[29], "30", CHARS_ION_STAGE );
		strncpy( chIonStage[30], "31", CHARS_ION_STAGE );

		strncpy( chIonRoman[0],  "I", CHARS_ION_ROMAN );
		strncpy( chIonRoman[1],  "II", CHARS_ION_ROMAN );
		strncpy( chIonRoman[2],  "III", CHARS_ION_ROMAN );
		strncpy( chIonRoman[3],  "IV", CHARS_ION_ROMAN );
		strncpy( chIonRoman[4],  "V", CHARS_ION_ROMAN );
		strncpy( chIonRoman[5],  "VI", CHARS_ION_ROMAN );
		strncpy( chIonRoman[6],  "VII", CHARS_ION_ROMAN );
		strncpy( chIonRoman[7],  "VIII", CHARS_ION_ROMAN );
		strncpy( chIonRoman[8],  "IX", CHARS_ION_ROMAN );
		strncpy( chIonRoman[9],  "X", CHARS_ION_ROMAN );
		strncpy( chIonRoman[10], "XI", CHARS_ION_ROMAN );
		strncpy( chIonRoman[11], "XII", CHARS_ION_ROMAN );
		strncpy( chIonRoman[12], "XIII", CHARS_ION_ROMAN );
		strncpy( chIonRoman[13], "XIV", CHARS_ION_ROMAN );
		strncpy( chIonRoman[14], "XV", CHARS_ION_ROMAN );
		strncpy( chIonRoman[15], "XVI", CHARS_ION_ROMAN );
		strncpy( chIonRoman[16], "XVII", CHARS_ION_ROMAN );
		strncpy( chIonRoman[17], "XVIII", CHARS_ION_ROMAN );
		strncpy( chIonRoman[18], "XIX", CHARS_ION_ROMAN );
		strncpy( chIonRoman[19], "XX", CHARS_ION_ROMAN );
		strncpy( chIonRoman[20], "XXI", CHARS_ION_ROMAN );
		strncpy( chIonRoman[21], "XXII", CHARS_ION_ROMAN );
		strncpy( chIonRoman[22], "XXIII", CHARS_ION_ROMAN );
		strncpy( chIonRoman[23], "XXIV", CHARS_ION_ROMAN );
		strncpy( chIonRoman[24], "XXV", CHARS_ION_ROMAN );
		strncpy( chIonRoman[25], "XXVI", CHARS_ION_ROMAN );
		strncpy( chIonRoman[26], "XXVII", CHARS_ION_ROMAN );
		strncpy( chIonRoman[27], "XXVIII", CHARS_ION_ROMAN );
		strncpy( chIonRoman[28], "XXIX", CHARS_ION_ROMAN );
		strncpy( chIonRoman[29], "XXX", CHARS_ION_ROMAN );
		strncpy( chIonRoman[30], "XXXI", CHARS_ION_ROMAN );
	};

};
extern t_elementnames elementnames;

/** isElementSym - search for string among element symbols.
 *
 * \param chSym	input string
 * \return	TRUE is an element, FALSE otherwise 
 */
bool isElementSym( const string& chSym );

/* elem_symbol_to_index - obtain the array index the input element symbol refers to.
 * String length must match element symbol string length.
 * The first letter is expected to be uppercase, and the second lowercase.
 *
 * \param chSym	input string
 * \return	the element index, or -1 if not found
 */
int elem_symbol_to_index( const string& chSym );

#endif /* ELEMENTNAMES_H_ */
