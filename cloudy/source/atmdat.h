/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ATMDAT_H_
#define ATMDAT_H_

#include "container_classes.h"
#include "module.h"

/*This structure is specifically used to hold the collision data in the format given in the LEIDEN Database
The data is available as collision rate coefficients(cm3 s-1) over different temperatures*/
class CollRateCoeffArray
{
public:
	/*Array of temps*/
	vector<double> temps;
	/*Matrix of collision rates(temp,up,lo)*/
	multi_arr<double,3> collrates;
	
}  ;

/*This structure is specifically used to hold the collision data in the format given in the CHIANTI Database
The data is available as spline fits to the Maxwellian averaged collision strengths */
class CollSplinesArray
{
public:
	/*Matrix of spline fits(hi,lo,spline index)*
	 *The first five columns gives the no of spline pts,transition type,gf value,delta E
	 *& Scaling parameter ,in the specified order*/
	/*The transition type basically tells how the temperature and collision
	strengths have been scaled*/
	vector<double> collspline;
	vector<double> SplineSecDer;

	long nSplinePts = 0L;
	long intTranType = 0L;
	double EnergyDiff;
	double ScalingParam;
	CollSplinesArray() : EnergyDiff(0.), ScalingParam(0.) {}
};

/*This structure is specifically used to hold the collision data in the format given in the STOUT Database
The data are available as collision strengths and rates over different temperatures*/
struct StoutColls
{
private:
	long m_offset, m_ntemps;
public:
	/*Number of temps*/
	long ntemps()
	{
		return m_ntemps;
	}
	/*Reset values to safe defaults*/
	void junk()
	{
		m_offset = -1;
		m_ntemps = -1;
		lgIsRate = false;
	}
	// Set range of tabulated data corresponding to this level & collider
	void setslice(long offset, long ntemps)
	{
		m_offset = offset;
		m_ntemps = ntemps;
	}
	long offset() const
	{
		return m_offset;
	}
	/*Is this a deexcitation rate or collision strength*/
	bool lgIsRate;

};

class StoutCollArray
{
	multi_arr<StoutColls,2> m_a;
	vector<double> m_temps,m_collstrs;
	static long ilev(long ipHi, long ipLo)
	{
		// Can use triangular arrangement, as ipHi > ipLo
		return (ipHi*(ipHi-1))/2+ipLo;
	}
public:
	// Allocates background data structures
	void alloc(long nHi, long , long nCollider)
	{
		// Size of array is first level which *isn't* required
		m_a.alloc(ilev(nHi+1,0),nCollider);
		m_temps.resize(0);
		m_collstrs.resize(0);
	}
	// Initializes per-collision data to safe values
	void junk(long ipHi, long ipLo, long ipCollider)
	{
		m_a[ilev(ipHi,ipLo)][ipCollider].junk();
	}
	bool& lgIsRate(long ipHi, long ipLo, long ipCollider)
	{
		return m_a[ilev(ipHi,ipLo)][ipCollider].lgIsRate;
	}
	// Reserves space for the number of table points for this level and
	// collider
	void setpoints(long ipHi, long ipLo, long ipCollider,long npoints)
	{
		long offset = m_temps.size();
		long length = offset+npoints;
		// testing shows that 1200 is a reasonable initial guess for the size
		m_temps.resize(max(length,1200));
		m_collstrs.resize(max(length,1200));
		m_a[ilev(ipHi,ipLo)][ipCollider].setslice(offset,npoints);
	}
	// Returns number of table points
	long ntemps(long ipHi, long ipLo, long ipCollider)
	{
		return m_a[ilev(ipHi,ipLo)][ipCollider].ntemps();
	}
	// Returns pointer to temperature ordinates
	double* temps(long ipHi, long ipLo, long ipCollider)
	{
		return &m_temps[m_a[ilev(ipHi,ipLo)][ipCollider].offset()];
	}
	// Returns pointer to collision data
	double* collstrs(long ipHi, long ipLo, long ipCollider)
	{
		return &m_collstrs[m_a[ilev(ipHi,ipLo)][ipCollider].offset()];
	}
};

/**
 * LoadIsotopes	read in the nuclear isotope data and allocate space
 */
void LoadIsotopes ( );


 /**
  atmdat_2phot_shapefunction two photon emission function for all atomic and ionic species 
  \param  EbyE2nu 
  \param  ipISO
  \param  nelem 
 */ 
double atmdat_2phot_shapefunction( double EbyE2nu, long ipISO, long nelem );

 /**
  atmdat_readin read in some data files, but only if this is very first call 
 */ 
void atmdat_readin(void);

 /**
 atmdat_STOUT_readin read in data from STOUT database files
 \param intNS
 \param chFileName
 */
void atmdat_STOUT_readin( long intNS, const string& chFileName );

 /**
  atmdat_CHIANTI_readin read in data from CHIANTI database files
  \param intNS
  \param chFileName
 */ 
void atmdat_CHIANTI_readin( long intNS, const string& chFileName );

 /**
  atmdat_LAMDA_readin read in data from LAMDA database files
  \param intNS
  \param chFileName
 */ 
void atmdat_LAMDA_readin( long intNS, const string& chFileName );


 /**
  atmdat_outer_shell determine outer shell, and statistical weights of that and higher ion, for any ion
  written by Dima Verner
  \param [in] iz  atomic number from 1 to 30
  \param [in] in  number of electrons from 1 to iz
  \param [out] *imax  number of the outer shell
  \param [out] *ig0   statistical weight of (iz,in) ground state
  \param [out] *ig1 statistical weight of (iz,in-1) ground state
  \author Dima Verner
 */ 
void atmdat_outer_shell(
  long int iz, 
  long int in,
  long int *imax,
  long int *ig0, 
  long int *ig1);

 /**
  atmdat fill in the CharExcIonOf[ipHYDROGEN] and Rec arrays with Kingdon's fitted CT with H, 
 */ 
void ChargTranEval( void );

/**
 sum up the charge transfer heating
 \return 
*/ 
double ChargTranSumHeat(void);

/*ChargTranPun save charge transfer rate coefficients */
 /** 
  save charge transfer rate coefficients
  \param ipPnunit 
  \param chSave 
 */ 
void ChargTranPun( FILE* ipPnunit , char* chSave );

/** CHIANTI_Upsilon converts Chianti collision splines to collision strengths */
double CHIANTI_Upsilon(long, long, long, long,double);

 /**
  atmdat_dielrec_fe Dielectronic recombination rates for Fe from Arnaud & Raymond 1992
  \param  ion
  \param  t
 */ 
double atmdat_dielrec_fe(long int ion, double t);

/** this initializes the arrays containing the fitting coefficients,
 * called by OpacityCreateAll, done once per coreload */
void atmdat_H_phot_cs(void);

/**atmdat_3body derive three-body recombination coefficients */
void atmdat_3body(void);

 /** 
 general utility to read in line emissivities from the Storey & Hummer tables
 of case B emissivities.  
\param iHi the principal quantum numbers, .	 
\param iLo upper and lower levels in any order
\param iZ charge of ion, only 1 and 2 for now	 
\param TempIn temperature, must lie within the range of the table, which depends on the ion charge, and is 500 - 30,000K for hydrogen
\param DenIn the density and must lie within the range of the table
\param chCase case - 'a' or 'b'
 */ 
double atmdat_HS_caseB( 
	long int iHi, long int iLo, long int iZ, double TempIn,
	double DenIn,
	char chCase 
	);

// arrays for Hummer & Storey 98 He1 cross sections and energies
extern multi_arr<double,4> HS_He1_Xsectn;
extern multi_arr<double,4> HS_He1_Energy;

// arrays for TOPbase Helike cross sections and energies
extern multi_arr<vector<double>,4> OP_Helike_Xsectn;
extern multi_arr<vector<double>,4> OP_Helike_Energy;
extern multi_arr<long,4> OP_Helike_NumPts;

/* these are the vectors that store the original Hummer and Storey case B
 * line data for H and He - the declaration for the interpolator follows */
#define NHSDIM 15 /**< used for following vectors*/
#define NLINEHS 300  /**< dimension of array with lines*/
#define HS_NZ 8 /**< number of elements that can be read in */
#define NHCSTE	8 /**< number of temperature points in h_coll_str arrays */
#define NUM_HS98_DATA_POINTS 811

struct t_atmdat : public module {
	const char *chName() const
	{
		return "atmdat";
	}
	void zero();
	void comment(t_warnings&) {}
	/**
	 * ion, nelem
	 * these arrays save the charge transfer ionization and recombination
	 * rates for the heavy elements onto hydrogen.  ionization is
	 * of the heavy element, and so is a recombination for hydrogen
	 * 
	 * CharExcIonOf[ipHYDROGEN]( ion , nelem ), CharExcRecTo[ipHYDROGEN]( ion , nelem )
	 * charge transfer ionization of atomic oxygen = CharExcIonOf[ipHYDROGEN][ipOXYGEN][0]*hii
	 * charge transfer recombination of ionized oxygen = CharExcRecTo[ipHYDROGEN][ipOXYGEN][0]*hi
	 * HCharHeatMax, HCharCoolMax are largest fractions of local heating
	 * or cooling due to ct
	 * HCharHeatOn usually 1, set to 0 with no CTHeat command
	 */

	enum {NCX=2};

	/** accurate ionization potentials in Ryd */
	double EIonPot[LIMELM][LIMELM];

	/** getIonPot -- get ionization potential, in Ryd
	 *
	 * \param nelem [in]   element index (0 for H, 29 for Zn)
	 * \param ion   [in]   ion charge (0 for neutral, nelem for H-like)
	 * \return double      ionization potential, in Ry
	 */
	inline double getIonPot( long int nelem, long int ion )
	{
		return EIonPot[nelem][ion];
	}

	/** CharExcIon is ionization, */
	/** [0] is Atom^0 + H+ => Atom+1 + H0
  	  * [n] is Atom^+n + H+ => Atom^+n-1 + H0 */
	/** CharExcRec is recombination */
	/** [0] is Atom^+1 + H0 => Atom^0 + H^+
	  * [n] is Atom^+n+1 + H0 => Atom^+n + H^+ */
	double CharExcIonOf[NCX][LIMELM][LIMELM+1], //(cm3 s-1)
	  CharExcRecTo[NCX][LIMELM][LIMELM+1],		//(cm3 s-1)
	  HCharHeatMax, 
	  HCharCoolMax, 
	  HCharHeatOn;

	/* rate coefficient (cm3 s-1) for N+(3P) + H+ -> N(2D) + H+ charge transfer*/
	double HCharExcRecTo_N0_2D;

	/** this is total rate (s-1) for ct ionization and recombination of nelem */
	double CharExcIonTotal[NCX],
		CharExcRecTotal[NCX];

	/** this is the current ratio of ct ionization of H, relative to total dest rate*/
	double HIonFrac;

	/** this is the largest ratio of ct ionization of H, relative to total dest rate*/
	double HIonFracMax;

	/** Dalgarno H charge transfer rate coefficient for high stages of ionization
	 * default is 1.92e-9 in zero, reset with 'set charge transfer' command */
	double HCTAlex;

	/** variable to turn on or off ct ionization-recombination of
	 * all elements - set off with no charge transfer command */
	bool lgCTOn;

	/** these are the density and temperature mesh points on the
	 * original Hummer & Storey data, for H[0] and He[1], */
	double Density[2][HS_NZ][NHSDIM], 
		ElecTemp[2][HS_NZ][NHSDIM],
		/**emiss[ipTemp][ipDens][ipLevel]*/
		Emiss[2][HS_NZ][NHSDIM][NHSDIM][NLINEHS];

	/** the number of density temperature mesh points for H&S Case A, B */
	long int nDensity[2][HS_NZ] , ntemp[2][HS_NZ];
	/** highest principal quantum number they report, usually 25, their paper says 50 */
	long int ncut[2][HS_NZ];

	/** following will be set false if we ever stop over bounds of HS table
	 * for any element.  first index is case A [0] or case B [1] -
	 * second is element number */
	bool lgHCaseBOK[2][HS_NZ];

	/** related to highest stage of ionization needed for Cota recom */
	long int nsbig;

	/** by default, include collisional ionization, option to not include it,
	 * with "no collisional ionization" command */
	bool lgCollIonOn;

	/** wavelengths of Hummer & Storey case B lines for H - O 
	 * first dimension is atomic number of C scale, H is 0
	 * next two are upper and lower configurations on physics 
	 * scale - Lya is 2-1, Lyb is 3-1, Ha is 3-2, etc */
	realnum WaveLengthCaseB[8][25][24];

	/** wavelengths for HeI case b */
	vector<realnum> CaseBWlHeI;

	const long nDefaultPhotoLevelsFe;
	/** Default number of non-iron levels when not using the coronal command */
	const long nDefaultPhotoLevels;
	/** Default number of iron levels for collisional ionization cases using the coronal command */
	const long nDefaultCollLevelsFe;
	/** Default number of non-iron levels for collisional ionization cases using the coronal command */
	const long nDefaultCollLevels;
	/** Default number of molecular levels */
	const long nDefaultMolLevels;
	// Set the constant member variables to the default values.

	/** true if CHIANTI database is enabled **/
	bool lgChiantiOn;
	/** true if CHIANTI database is supplemented by opacity project lines */
	bool lgChiantiLvl2Hybrid;
	/** true if Cloudy will print which Chianti species are being used as well as number of levels */
	bool lgChiantiPrint;
	/** true if Cloudy will use no theoretical energy levels from Chianti, only experimental. False means that only theoretical energy levels are used */
	bool lgChiantiExp;
	/**CloudyChianti filename variable **/
	char chCloudyChiantiFile[FILENAME_PATH_LENGTH];
	/** The maximum number of chianti energy levels used for Fe */
	long nChiantiMaxLevelsFe;
	/** The maximum number of chianti energy levels used for all other species */
	long nChiantiMaxLevels;
	/** Flag to determine whether nChiantiMaxLevelsFe has been set by the user */
	bool lgChiantiLevelsSet;

	/** true if LAMDA database is enabled **/
	bool lgLamdaOn;
	/** true if Cloudy will print which Lamda species are being used as well as number of levels */
	bool lgLamdaPrint;
	/**LAMDA filename variable **/
	char chLamdaFile[FILENAME_PATH_LENGTH];
	/** maximum number of lamda energy levels */
	long nLamdaMaxLevels;
	/** Flag to determine whether nLamdaMaxLevels has been set by the user */
	bool lgLamdaLevelsSet;

	/** true if Stout database is enabled **/
	bool lgStoutOn;
	/** true if Stout database supplements opacity project (Level 2) lines */
	bool lgStoutLvl2Hybrid;
	/** true if Cloudy will print which Stout species are being used as well as number of levels */
	bool lgStoutPrint;
	/**Stout filename variable **/
	char chStoutFile[FILENAME_PATH_LENGTH];
	/** The default maximum number of stout energy levels for Fe */
	long nStoutMaxLevelsFe;
	/** default maximum number of stout energy levels used for other
	 *  species */
	long nStoutMaxLevels;
	/** Flag to determine whether nStoutMaxLevelsFe has been set by the user */
	bool lgStoutLevelsSet;
	/** Default number of iron levels when not using the coronal command */

	/** true if CDMS/JPL database is enabled **/
	bool lgCalpgmOn;
	/** true if dBase transitions have collision strengths supplemented by gbar **/
	bool lgGbarOn;

	/** The default collision strength used for dBaseTrans
	 * no collision or radiative data are available so 
	 * conventional g-bar is unavailable **/ 
	double collstrDefault;

	/** this says whether to include inner shell absorption lines */
	bool lgInnerShellLine_on;
	/** says whether to include the new Romas data set */
	bool lgInnerShell_Kisielius;
	/** flag to print out the UTA data */
	bool lgUTAprint;

	/** The threshold for Aul. For species where the highest level has no
	 * transitions with Aul > aulThreshold, the highest levels will be trimmed until a real Aul is found.*/
	static const double aulThreshold;

	/* type and enum for determining what collisional ionization rate coefficient data to use.
	* Dima: Cloudy original data from Voronov97
	* Hybrid: Dima version scaled by the ratio of Dere07 to Dima */
	typedef enum { DIMA, HYBRID } CollIonRC;
	CollIonRC CIRCData;

	/** Chianti version read from /data/chianti/VERSION **/
	string chVersion;

	/** Array of Data Sources **/
	char chdBaseSources[LIMELM][LIMELM+1][10];
	/** this indicates that a model atom has been set up for this ion **/
	bool lgdBaseSourceExists[LIMELM][LIMELM+1];

	/** Array to store species labels and ipSpecies values **/
	multi_arr<long,2> ipSpecIon;

	t_atmdat() : nDefaultPhotoLevelsFe(25), nDefaultPhotoLevels(15), nDefaultCollLevelsFe(100),
			nDefaultCollLevels(50), nDefaultMolLevels(70)
	{
		//set all sources to blank and false.
		//Sources will be added from masterlist files in species.cpp
		for( int nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
		{
			for( int ion=0; ion < nelem+2; ++ion )
			{
				strcpy(chdBaseSources[nelem][ion]," ");
				lgdBaseSourceExists[nelem][ion] = false;
			}
		}

		/** set atmdat.ipSpecIon to uninit flag **/
		ipSpecIon.alloc(LIMELM,LIMELM);
		ipSpecIon = -1;
	}
};

extern t_atmdat atmdat;

class Funct 
{
public:
	virtual void operator()( long&, long&, const char*, long&) = 0;
	virtual ~Funct() = 0;
};
 
inline Funct::~Funct() {}

typedef Funct* FunctPtr;

class FunctLAMDA : public Funct 
{
public:
	explicit FunctLAMDA(void) { }
	virtual void operator()( long& ipHi, long& ipLo, const char* chLine, long& i) 
	{ 
		DEBUG_ENTRY( "FunctLAMDA()" );
		bool lgEOL;
		long index = (long)FFmtRead( chLine, &i, strlen(chLine), &lgEOL );
		if( index <= 0 )
		{
			fprintf( ioQQQ, " PROBLEM with LAMDA readin.\n");
			cdEXIT( EXIT_FAILURE );
		}
		ipHi = (long)FFmtRead( chLine, &i, strlen(chLine), &lgEOL ) - 1;
		ipLo = (long)FFmtRead( chLine, &i, strlen(chLine), &lgEOL ) - 1;
		return;
	}
private:
};

#include "h2_priv.h"

class FunctDiatoms : public Funct 
{
public:
	explicit FunctDiatoms(const diatomics& diatom) : diatom_(diatom) { }
	virtual void operator()( long& ipHi, long& ipLo, const char* chLine, long& i) 
	{ 
		diatom_.GetIndices( ipHi, ipLo, chLine, i );
	}
private:
	const diatomics& diatom_;
};

void ReadCollisionRateTable( CollRateCoeffArray& coll_rate_table, FILE* io, FunctPtr GetIndices, long nMolLevs, long nTemps = -1, long nTrans = -1 );

double InterpCollRate( const CollRateCoeffArray& rate_table, const long& ipHi, const long& ipLo, const double& ftemp);


/* Form a comment for a transition obtained from an external database.
 * The comment is the configuration of the lower and upper level */
inline string db_comment_tran_levels(const string& clo = "", const string& chi = "")
{
	stringstream Comment;
	if( clo.length() > 0 && chi.length() > 0 )
		Comment << clo << " -- " << chi;
	return Comment.str();
}

#endif /* ATMDAT_H_ */
