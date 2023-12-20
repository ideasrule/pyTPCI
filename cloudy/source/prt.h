/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PRT_H_
#define PRT_H_

#include "module.h"
#include "container_classes.h"
#include "lines.h"

class TransitionProxy;

//* Maximum number of columns in output print
const long NCOLMAX=132;

/**PrtZone print out individual zone results */
void PrtZone(void);

/**PrtComment analyze model, generating comments on its features */
void PrtComment(void);

/**PrtFinal create final pages of printout, emission line intensities, etc */
void PrtFinal(void);


/**SetPrintLineCol	set main line block & wl printing formats */
void SetPrintLineCol ();

/**prt_wl write wavelength to io 
\param *io
\param wavelength
*/
void prt_wl( 
	FILE *io , 
	realnum wavelength );

/**sprt_wl write wavelength to string - must be kept parallel with prt_wl 
\param *chString
\param wl
*/
void sprt_wl( 
	string& chString,
	realnum wl );

/** prt_line_err produce an error message containing the line label and wavelength,
 *              followed, if given, by the wavelength of the closest line of the same label 
\param *ioOUT		output file handle
\param *label		line label
\param wvlng		line wavelength
 */
void prt_line_err( FILE *ioOUT, const string& label, realnum wvlng );
void prt_line_err( FILE *ioOUT, const LineID& lineid );

/* prt_line_inlist print line suitable for output list, label not enclosed in quotation marks
\param *ioOUT		output file handle
\param *label		line label
\param wvlng		line wavelength
 */
void prt_line_inlist ( FILE *ioOUT, const char *label, realnum wvlng );

/**PrtHeader print large block of incident continuum numbers at start, 
 just after echoing input commands */
void PrtHeader(void);

/**prt_LineLabels save all labels and wavelengths for emission line array 
 *
 * \param io            file handle to write output
 * \param lgPrintAll    print all if true, if false then do not print parts of 
 *                      transferred lines
 * \param lgPrintIndex  print the index of the line in the line stack
*/
void prt_LineLabels(
	FILE * io,
	bool lgPrintAll,
	bool lgPrintIndex
	);

/**prtmet print all line optical depths at end of iteration */
void prtmet(void);

/**PrtMeanIon print mean ionization fractions for all elements,
 * output will go to stream pointed to by argument  
 * chTyp is either 'i' or 't' for mean ionization or temperature 
 \param chType
 \param lgDensity true include density, false do not
 */
void PrtMeanIon( char chType , 
			bool lgDensity,
			FILE *);

/**PrtLineSum parse print line sum command to enter set of lines into sum  
\param chDo the job to do, either " SUM" or "READ"
*/
double PrtLineSum(void);

/**PrtLinePres print line radiation pressures for current conditions 
 * output goes top openned file handle */
void PrtLinePres(FILE *ioPRESSURE);

/**PrtColumns print column densities of all elements in standard output
\param ioMEAN this is stream used for io, is stdout when called by final,
       is save unit when save output generated
*/
void PrtColumns(
	 FILE *ioMEAN  );

/** CloudyPrintReference print preferred citation to Cloudy */
void CloudyPrintReference();

/** DatabasePrintReference print some database references */
void DatabasePrintReference();

/**PrtAllTau master routine controlling printout of optical depths at
 end of calculation */
void PrtAllTau(void);

class t_prt_matrix : public module {
public:
	/** species element and ionization stage set with print array command to print
	  * matrixes input to solvers */
	string species;
	string speciesLevels;
	vector<long> speciesLevelList;
	bool lgLevelsResolved = false;

	void zero();
	void comment(t_warnings&) {}

	const char *chName() const
	{
		return "prt_matrix";
	}

	void setSpecies( const string &sspec );
	void resolveLevels();
	void prtRates( const long numLevels,
			const multi_arr<double,2,C_TYPE> &matrix,
			valarray<double> &b );
};

/** struct for holding user-defined blend */
struct t_blend {
	string chLabel;
	realnum wave;
	bool lgQuiet;
	bool lgIgnore;
	vector<LineID> component;
	t_blend() : chLabel("Blnd"), wave(0_r), lgQuiet(false), lgIgnore(false) {}
};

struct t_prt {

	/** lgPrintBlock, option to turn off printing of the main line blocks */
	bool lgPrintBlock;

	/** lgPrintBlockIntrinsic, option to turn off printing of the intrinsic
	  * line blocks */
	bool lgPrintBlockIntrinsic;

	/** lgPrintBlockEmergent, option to turn off printing of the emergent
	  * line blocks */
	bool lgPrintBlockEmergent;

	/** lgSortLines, option to sort lines by wavelength- print sort
	 command */
	bool lgSortLines;

	/** if above is set, then one of the following must also be set,
	 * say whether to sort by wavelength or intensity */
	bool lgSortLineWavelength , lgSortLineIntensity;

	/** lower and upper wavelength bounds for printed spectrum,
	 * range option on print sort command */
	realnum wlSort1 , wlSort2;

	/** print hydrogenic level populations, 
	 * set with print hydrogenic command
	bool lgPrintHLevPops; */

	/** print column densities */
	bool lgPrintColumns;

	/** should we print execution time?  normally true, but set false
	 * with no times command so that different runs can compare exactly */
	bool lgPrintTime;

	/** print ages command tells code to print various timescales */
	bool lgPrnAges;

	/** option to print maser lines (true) normally false
	 * print maser turns on */
	bool lgPrtMaser;

	/** lgPrtTau tells whether to print line optical depths */
	bool lgPrtTau;

	/** lgPrintFluxEarth says to print flux of lines at Earth, 
	 * if luminosity can be predicted */
	bool lgPrintFluxEarth;

	/** print line surface brightness command, units either sr or sq arcsec,
	 * default is SR, set to arcsec with arcsec option */
	bool lgSurfaceBrightness , lgSurfaceBrightness_SR;

	/** PrtTauFnt is smallest line optical depth to print */
	realnum PrtTauFnt;

	/** these are various contributors to the line output,
	 * and are changed with the print line or print continuum commands
	 * in prtfinal code uses these to make a final filter over what lines
	 * will be printed */
	bool lgPrnPump, 
	  lgPrnHeat, 
	  lgPrnColl, 
	  lgPrnInwd;

	/** print predictions from collapsed levels of iso sequences,
	 * print line iso collapsed */
	bool lgPrnIsoCollapsed;

	/* flag set with print continuum index command, to identify all lines
	 * that lie within a continuum cell */
	bool lgPrtContIndices;
	/* these are lower and upper limits to the energy range in Rydbergs.
	 * they are the first and second number on the command line, lower and
	 * upper bounds of the code are used if not specified */
	realnum lgPrtContIndices_lo_E , 
		lgPrtContIndices_hi_E;

	/** flags for determining what is included in nFnu */
	bool lgSourceReflected;
	bool lgSourceTransmitted;
	bool lgDiffuseInward;
	bool lgDiffuseOutward;

	/** flag set with print departure coefficients */
	bool lgPrtBN;

	/** if true then print only last iteration */
	bool lgPrtLastIt;

	/** flag set with print short command */
	bool lgPrtShort;

	/** lgOnlyZone set with print only zones */
	bool lgOnlyZone;
	/** lgOnlyHead set with print only header */
	bool lgOnlyHead;

	/** lgPrtStart is option to start printout at certain zone */
	bool lgPrtStart;

	/**nstart is zone number, set with print start command */
	long int nstart;

	/** flag to turn on printout of heat sources */
	bool lgPrintHeating;

	/** flag set with print array command to print ionization recombination arrays */
	bool lgPrtArry[LIMELM];

	t_prt_matrix matrix;

	/** logical lgFaintOn normally true, says to not print very faint lines
	  set false with print faint off command
	 lines fainter than TooFaint will not be printed.  This is set in 
	 zerologic and reset with print line faint command  */
	realnum TooFaint;
	bool lgFaintOn;

	/** flag set true if print faint command entered,
	 * only used to not override it with print short */
	bool lgFntSet;

	/** these implement the print line cell commmand, 
	 * flag saying to do this */
	bool lgPrnLineCell;
	/** the cell number, on the physics scale, counts from 1, 
	 * will print labels of all lines that lie within that cell */
	long int nPrnLineCell;

	/** option to print main block of lines as a single column
	 * instead of the normal array.  if true then usual array */
	bool lgPrtLineArray;

	/** printing as a column also has an option to print linear quantity 
	 * in exponential format */
	bool lgPrtLineLog;

	/** flag set by print line cumulative command, also print large set of
	 * emission line integrated intensities over time depend model */
	bool lgPrintLineCumulative;

	/** use air wavelengths for wl > 2000A, as per atomic physics tradition
	 * dating back to 19th century.  SDSS does not follow this tradition and
	 * uses vacuum wavelengths - PRINT LINE VACUUM will do this.
	 */
	bool lgPrintLineAirWavelengths;

	/* Generate output in HTML format */
	bool lgPrintHTML;

	/** quantities to do with radiation field and printed in header */
	realnum qx, 
	  powion, 
	  xpow, 
	  pbal, 
	  q, 
	  qgam, 
	  pradio, 
	  fx1ryd;
	long int ipeak;
	realnum GammaLumin;

	long int nzdump;

	/** print citations command flag */
	bool lgPrtCitations;

	/** array of blends from data/blends.ini */
	vector<t_blend> blend;

	/** should blends.ini be included? */
	bool lgIncludeBlends;

	t_prt()
	{
		// make sure this has the correct value before main() starts
		// this is needed by check_data() and possibly others
		lgPrintTime = true;
		// this needs to be set before the code starts reading the input script
		lgIncludeBlends = true;
	}
};
extern t_prt prt;



struct t_line_col : public module
{
	const char* chName() const
	{
		return "prt_linecol";
	}
	void zero();
	void comment(t_warnings&) {}

	/** width of absolute line intensity string in output column */
	int	absint_len;

	/** width of relative line intensity string in output column */
	int	relint_len;

	/** width of line column, compiled from the number of widths
	 * of the label, the wavelength, and the two intensities above */
	int	column_len;

	/** width of space between columns in line output, used unless
	 * single column requested */
	int	col_gap_len;

	/** string of stars signifies out-of-range values in relint */
	string	relint_outrange;

	/** gap bw columns, based on column_gap_len above */
	string	col_gap;
};
extern struct t_line_col prt_linecol;

#endif /* PRT_H_ */
