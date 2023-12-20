/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef SAVE_H_
#define SAVE_H_

#include "energy.h"
#include "lines.h"
#include "prt.h"

class diatomics;
class TransitionProxy;

/* save.h */
static const long LIMPUN = 100L;

/** magic version number for the transmitted continuum output file */
static const long VERSION_TRNCON = 20190210L;

/**SaveDo produce save output during calculation 
\param chTime chTime is null terminated 4 char string, either "MIDL" or "LAST" 
*/
void SaveDo(
	const char *chTime); 

/** SAVE command has option LOG to print log quantities as in <= C13 */
double PrtLogLin( double value );

/** save one line, called by SaveLineStuff 
\param t
\param io
\param xLimit
\param index
\param DopplerWidth
*/
void Save1Line( 
	const TransitionProxy & t , 
	FILE * io , 
	realnum xLimit , 
	long index,
	realnum DopplerWidth);

/**SaveLineData punches selected line data for all lines transferred in code 
\param io
*/
NORETURN void SaveLineData(
	FILE * io);

/**save_opacity save total opacity in any element, save opacity command 
\param io
\param np
*/
void save_opacity(
	FILE * io, 
	long int np);

/**SaveSpecial generate output for the save special command 
\param io
\param chTime
*/
void SaveSpecial(
	FILE* io , 
  const char *chTime);


/**SaveAllSpeciesLabelsLevels -- generate output of all species labels & levels
 * \param ioPUN		io unit for save
*/
void SaveAllSpeciesLabelsLevels( FILE *ioPUN );


/**SaveSpecies generate output for the save species command 
\param ioPUN io unit for save
\param ipPun arrayin save array
*/
void SaveSpecies(
		FILE* ioPUN,
		long int ipPun );

/** SaveSpeciesPseudoCont - save species pseudo-continuum
 *
 * \param ipPun		save file index
 * \param speciesLabel	species emission to be stored
 */
void SaveSpeciesPseudoCont( const long ipPun, const string &speciesLabel );

/** SaveSpeciesBands - save species bands emission
 *
 * \param ipPun		[in]	save file index
 * \param speciesLabel	[in]	species emission to be stored
 * \param fileBands	[in]	file that holds the bands' definitions
 * \param lgEmergent	[in]	report emergent emission; intrinsic by default
 */
void SaveSpeciesBands( const long ipPun, const string &speciesLabel,
			const string &fileBands, const bool lgEmergent );

/**SaveSpeciesOptDep -- generate output for 'save species optical depth'
 *
 * \param ipPun		array index in save array
 * \param speciesLabel	label to acquire species
 */
void SaveSpeciesOptDep( const long int ipPun, const string &speciesLabel );

class molecule;

extern void mole_save(FILE *punit, const char speciesname[], const char args[], bool lgHeader, bool lgData, bool lgCoef, double depth);

// debug_list -- list of species to print -- add multiple species to get net
//    rate of change for all together, useful for understanding net sources
//    sinks to networks
// ioOut -- stream
// lgPrintReagents -- print abundances of reagents
// NPRINT -- max number of reactions
// fprint -- minimum fraction of total rate which is worth printing
extern void mole_dominant_rates( const vector<const molecule *>& debug_list, 
											FILE *ioOut,
											bool lgPrintReagents, size_t NPRINT, double fprint );

extern void mole_print_species_reactions( molecule *speciesToPrint );


/* Print header for Save1LineData() function
 * \param ioPUN		file pointer to write to
 */
void PrintLineDataHeader( FILE * ioPUN );


/**Save1LineData data for save one line 
\param t
\param io
\param lgCS_2 this flag says whether collision strength should be saveed - should be false
	for multi level atoms since sums are not done properly 
*/
void Save1LineData( 
	const TransitionProxy & t , 
	FILE * io,
	bool lgCS_2 );

/**save_line do the save output 
\param ip the file we will write to 
\param chDo
\param lgEmergent intrinsic or emergent emission
\param ipPun
*/
void save_line(
	FILE * ip,
	const char *chDo, 
	bool lgEmergent,
	long ipPun );

/**save_average parse save average command, or actually do the save output 
\param ipPun - array index for file for final save output 
*/
void save_average( 
	long int ipPun);

/**Save_Line_RT parse the save line rt command - read in a set of lines 
\param ip the file we will write to
*/
void Save_Line_RT(
	FILE * ip);

/** Save spectra to a FITS compatible file. 
\param io the file we will write to
\param option choose type of spectrum to save
\param Elo if set, lower bound of spectrum to save
\param Ehi if set, upper bound of spectrum to save
\param Enorm if set, normalize spectrum to 1 at this frequency
*/
void saveFITSfile( 
	FILE* io, 
	int option,
	realnum Elo = 0.f,
	realnum Ehi = 0.f,
	realnum Enorm = 0.f);

/**Save rate matrix to a FITS file.  Mainly a debugging aid.
 * \param io       [in]	the file we will write to
 * \param extName  [in]	extension name
 * \param units    [in]	vector units
 * \param nPixels  [in]	dimension of (square) image
 * \param image    [in]	image to write in FITS extension
 */
void saveFITSimg( FILE *io, const string &extName, const string &units,
		const long nPixels,
		const multi_arr<double,2,C_TYPE> &image );

/**Save rate matrix to a FITS file.  Mainly a debugging aid.
 * \param io       [in]	the file we will write to
 * \param extName  [in]	extension name
 * \param units    [in]	vector units
 * \param nCols    [in]	length of vector
 * \param vector   [in]	vector to write in FITS extension
 */
void saveFITSimg( FILE *io, const string &extName, const string &units,
		const long nCols, const valarray<double> &vec );

/**SaveHeat save contributors to local heating, with save heat command, called by save_do 
\param io the file we will write to
*/
void SaveHeat(FILE* io);

/**CoolSave save coolants, called by save_do 
\param io the file we will write to
\param *chJob
*/
void CoolSave(FILE * io, const char chJob[]);

/**SaveGrid implements the SAVE GRID command
\param pnunit : file descriptor to write to
\param status : exit status of the grid point
*/
void SaveGrid(FILE* pnunit, exit_type status);

class SaveParams
{
public:
	/** this is the file where we will direct the output */
	FILE *ipPnunit;
};

class adjPseudoCont
{
public:
	string speciesLabel;
	realnum wlLo,
		wlHi;
	long nBins;
};

class save_img_matrix : public t_prt_matrix
{
public:
	bool lgImgRates;
	bool lgFITS;
	long iteration;
	long zone;

	void zero();
	void comment(t_warnings&) {}

	const char *chName() const
	{
		return "img_matrix";
	}

	inline bool matchIteration( const long this_iteration ) const
	{
		return ( ( iteration > 0 && iteration == this_iteration ) ||
				! iteration );
	}

	inline bool matchZone( const long this_zone ) const
	{
		return ( ( zone > 0 && zone == this_zone ) || ! zone );
	}

	/** Create image (FITS or PPM) using the given species.
	 *  Useful for generating FITS images when negative populations are
	 *  computed.
	 */
	void createImage( const string &species,
			const long iteration,
			const long nzone,
			const long numLevels,
			const multi_arr<double,2,C_TYPE> &matrix,
			const valarray<double> &creation,
			bool haveNegPop = false );

	/** Create image (FITS or PPM) using the species given in the
	 *  'save arrays levels' command.
	 */
	void createImage(
			const long iteration,
			const long nzone,
			const long numLevels,
			const multi_arr<double,2,C_TYPE> &matrix,
			const valarray<double> &creation,
			bool haveNegPop = false );

	/** Add FITS image to an existing file using the given species.
	 *  Useful for generating FITS images when negative populations are
	 *  computed.
	 */
	void addImagePop_FITS( const string &species,
			const long iteration,
			const long nzone,
			const long numLevels,
			const valarray<double> &pop,
			bool haveNegPop = false );

	/** Add FITS image to an existing file using the species given in the
	 *  'save arrays levels' command.
	 */
	void addImagePop_FITS(
			const long iteration,
			const long nzone,
			const long numLevels,
			const valarray<double> &pop,
			bool haveNegPop = false );

private:
	string set_basename( const string &species,
			const long iteration,
			const long nzone,
			bool haveNegPop = false );

	void createImage_PPM( const string &basename,
			const long numLevels,
			const multi_arr<double,2,C_TYPE> &matrix );

	void createImage_FITS( const string &basename,
			const long numLevels,
			const multi_arr<double,2,C_TYPE> &matrix,
			const valarray<double> &creation );
};

struct t_save {

	t_save()
	{
		for( long i=0; i < LIMPUN; ++i )
		{
			nLineList[i] = -1;
			lgFITS[i] = false;
			lgXSPEC[i] = false;
			FITStype[i] = -1;
			nAverageList[i] = -1;
			p_lgSaveTitleDone[i] = false;
			p_lgSaveHeaderDone[i] = false;
		}
		p_lgSaveHeaderNeeded = true;
		ipSaveGrid = -1;
	}

	~t_save()
	{
		for( long i=0; i < LIMPUN; ++i )
		{
			SaveLineListFree(i);
			SaveAverageFree(i);
		}
	}

	void SaveLineListFree(long i)
	{
		LineList[i].clear();
	}
	void SaveAverageFree(long i)
	{
		nAverageList[i] = 0;
		chAverageType[i].clear();
		chSaveSpecies[i].clear();
		/* nAverageIonList is set of ions for averages */
		nAverageIonList[i].clear();
		/* nAverage2ndPar is set of second parameters for averages */
		nAverage2ndPar[i].clear();
	}

	string chGridDelimeter(long n)
	{
		ostringstream oss;
		// NB NB -- do not change this format, tests in
		// check_grid_file() and fix_grid_file() depend on this
		oss << chHashString << " GRID_DELIMIT -- grid";
		oss << setfill('0') << setw(9) << n;
		return oss.str();
	}

	/** following are for save LineList option
	 * nLineList is number of em lines, -1 if not defined */
	long nLineList[LIMPUN];
	/** LineID's for line list */
	vector<LineID> LineList[LIMPUN];
	/** flag saying whether to take ratio (true) of pairs */
	bool lgLineListRatio[LIMPUN];

	/** following are for save averages option
	* nAverageList is number of averages, -1 if not defined */
	long nAverageList[LIMPUN];
	/** chAverageType is label for type of average */
	vector<string> chAverageType[LIMPUN];
	/** nAverageIonList is set of ions for averages */
	vector<int> nAverageIonList[LIMPUN];
	/** nAverage2ndPar is set of second parameters for averages */
	vector<int> nAverage2ndPar[LIMPUN];

	/** this is the file where we will direct the output */
	SaveParams params[LIMPUN]; 

	/** the filename is needed during the file gathering phase in grid runs */
	vector<string> chFileName;

	/** pointer into chFileName for the name of the SAVE GRID file */
	long ipSaveGrid;

	/** option to append instead of overwrite */
	bool lgNoClobber[LIMPUN];

	/** print old style logs */
	bool lgPrtOldStyleLogs[LIMPUN];

	/** these are for some special cases, same purpose as previous no clobber */
	bool lgPunConv_noclobber;
	bool lgDROn_noclobber;
	bool lgDTOn_noclobber;
	bool lgPunPoint_noclobber;
	bool lgioRecom_noclobber;
	bool lgQHSaveFile_noclobber;
	bool lgTraceConvergeBase_noclobber;
	bool lgSaveGrid_noclobber;

	/** is this a real save command, or one of the similar options like
	 * save dr, which is not done in save files */
	bool lgRealSave[LIMPUN]; 

	/** emergent or intrinsic line intensities */
	bool lgEmergent[LIMPUN];

	/** cumulative line intensities for time dependent sim */
	bool lgCumulative[LIMPUN];

	/** user-defined option to subtract continuum in reported line flux */
	bool lgSubtrCont;

	/** option to remove isotropic continua from the output of 'save continuum' */
	bool lgPrtIsotropicCont[LIMPUN];

	/** number of save commands entered */
	long int nsave;

	/**chSave - what is it we want to save? set in GetPunch, used in DoPunch */
	char chSave[LIMPUN][5];
	diatomics *whichDiatomToPrint[LIMPUN];

	/** which opacity to save out  */
	char chOpcTyp[LIMPUN][5];

	/** this flag tells us whether to save results of a grid to separate files
	 * for each grid point or all to the same file.  Different for different 
	 * save commands */
	bool lgSaveToSeparateFiles[LIMPUN];

	/** option to not insert end-of-iteration separator - used for save files
	 * that create one line per iteration */
	bool lg_separate_iterations[LIMPUN];

private:
	/** flag saying whether the title has been saved already.
	 * Used in grid runs so that the title only gets saved once. */
	bool p_lgSaveTitleDone[LIMPUN];

	/** flag saying whether the header has been saved already.
	 * Used in grid runs so that the header only gets saved once. */
	bool p_lgSaveHeaderDone[LIMPUN];

	/** flag saying whether a header is needed at all
	 * Used in grid runs so that the header only gets saved once. */
	bool p_lgSaveHeaderNeeded;

public:
	void SetSaveHeaderNeeded(bool lgVal) { p_lgSaveHeaderNeeded = lgVal; }

	/** is it OK to save title now? */
	bool lgSaveTitle(int ipPun) const
	{
		return !p_lgSaveTitleDone[ipPun] && p_lgSaveHeaderNeeded;
	}
	void SaveTitleDone(int ipPun) { p_lgSaveTitleDone[ipPun] = true; }

	/** is it OK to save header now? */
	bool lgSaveHeader(int ipPun) const
	{
		return !p_lgSaveHeaderDone[ipPun] && p_lgSaveHeaderNeeded;
	}
	void SaveHeaderDone(int ipPun) { p_lgSaveHeaderDone[ipPun] = true; }

	/** flag saying whether any save continuum commands were entered
	 * set true in parsecontinuum when save continuum entered,
	 * used in PrtComment to warn if continuum punched with no iterations */
	bool lgPunContinuum;

	/** punarg is set of optional arguments for the save command */
	realnum punarg[LIMPUN][3];

	/** punarg is set of optional arguments for the save command */
	string optname[LIMPUN];

	/** implement save every option - lgSaveEveryZone true if want to save
	 *every zone, nSaveEveryZone is number of zones to save */
	bool lgSaveEveryZone[LIMPUN];
	long int nSaveEveryZone[LIMPUN];

	/** set of optional arguments for save command, but as a string */
	char chSaveArgs[LIMPUN][5];

	vector<string> chSaveSpecies[LIMPUN];

	/** lg flag lgPunLstIter for this save option,
	 do we only want to save on last iteration? */
	bool lgPunLstIter[LIMPUN];

	/** flag saying that this save file is in FITS format */
	bool lgFITS[LIMPUN];

	/** flag saying that this save file is in XSPEC format */
	bool lgXSPEC[LIMPUN];

	/** which FITS type is in this file */
	int FITStype[LIMPUN];

	/** option to say whether any FITS output should be punched, 
	 * initialized to false, but turned to true after last grid exec.  */
	bool lgPunchFits;

	/**chConPunEnr - units of continuum in save output */
	const char *chConSavEnr[LIMPUN];

	/** this global variable is index of save command loop in dopunch */
	long int ipConPun;

	/** should hash marks be printed after every iteration?
	 * default is yes, set no with no hash option on save command */
	bool lgHashEndIter[LIMPUN];

	/** this is the hash string, normally a set of hash marks, can be reset with
	 * set save hash command */
	string chHashString;

	/** flush file after every iteration */
	bool lgFLUSH;

	/** in save continuum output, in the luminosity case, use old style of
	 * luminosity per inner cloud area rather than true luminosity */
	bool lgLuminosityOld;

	/** this is a prefix that will be used at the start of all file names
	 * when doing an MPI grid run, normally an empty string, set with -g flag */
	string chGridPrefix;

	/** this is a prefix that will be set at the start of all save file names
	 * normally an empty string, set with PUNCH PREFIX command or -p flag */
	string chFilenamePrefix;

	/** this is the prefix that will be set at the start of the input and output file
	 * normally an empty string, set with -p or -r flag */
	string chRedirectPrefix;

	/** this is the name of the main output file we are writing to
	 * it is set by calling cdOutput() */
	string chOutputFile;

	/** set with save line intensities and save results commands,
	 * says whether results arrays produced by routine PunResults1Line should
	 * be column or array */
	char chPunRltType[7];

	/** option to save out pointers with save pointers command
	 * ipPoint is save file handle, lgPunPoint says whether we will do it
	 */
	FILE* ipPoint;
	bool lgPunPoint;

	/** unit number, and flag, for saving reason for continued iterations */
	bool lgPunConv;
	FILE* ipPunConv;


	/** Identify whether Save Data Sources is on */
	bool lgSDSOn;
	FILE* ipSDSFile;

	/** these control saving choice of dr - this is not really a save command
	 * ipDRout is io unit, lgDROn says saving dr logic has been set,
	 * and lgDRPLst says to save the last iteration */
	FILE * ipDRout;
	bool lgDROn, 
	  lgDRPLst,
	  lgDRHash;

        /** Likewise for saving choice of dt */
	FILE * ipDTout;
	bool lgDTOn,
		lgDTPLst,
		lgDTHash;

	/* set true save convergence base */
	bool lgTraceConvergeBase,
		lgTraceConvergeBaseHash;
	FILE * ipTraceConvergeBase;

	/** option to save recombination coefficients to external file */
	FILE* ioRecom;
	bool lgioRecom;

	/** option to save line intensities every for every zone
	 * logical variable says whether LinEvery was set */
	long int LinEvery;
	bool lgLinEvery;

	/** set skip sets this variable, which says how many cells to skip in save */
	long int ncSaveSkip;

	/** threshold for faintest cooling or heating to be punched
	 * default is set to 0.05 in scalar, can be reset with 'set weakheatcool' */
	realnum WeakHeatCool;

	/** contrast factor for lines to continuum in save output
	* default is 1 (gives correct line intensities) and changed with
	* set width command (enters width in km/sec) */
	realnum Resolution;
	// same thing but for absorption lines - default is unity, set to
	// SaveLWidth if ABSORPTION keyword occurs on save line width command
	realnum ResolutionAbs;

	string chSpeciesDominantRates[LIMPUN];

	// the frequency at which the continuum volume emissivity should be saved
	Energy emisfreq[LIMPUN];
	long ipEmisFreq[LIMPUN];

	// Switches for save Database command
	// Wavenumbers vs WLAng, Gf vs A, Collision Rates or not
	bool lgSaveDataWn, lgSaveDataGf, lgSaveDataRates;

	/** Parameters for saving the pseudocontinuum */
	vector<string> contSaveSpeciesLabel;
	vector<adjPseudoCont> setPseudoCont;

	/** Parameters for species bands */
	string SpeciesBandFile[LIMPUN];

	save_img_matrix img_matrix;
};

extern t_save save;

#endif /* SAVE_H_ */
