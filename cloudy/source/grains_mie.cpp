/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "elementnames.h"
#include "dense.h"
#include "called.h"
#include "version.h"
#include "grainvar.h"
#include "rfield.h"
#include "atmdat_adfa.h"
#include "grains.h"

/*=======================================================*
 *
 * Mie code for spherical grains.
 *
 * Calculates <pi*a^2*Q_abs>, <pi*a^2*Q_sct>, and (1-<g>)
 * for arbitrary grain species and size distributions.
 *
 * This code is derived from the program cmieuvx.f
 *
 * Written by: P.G. Martin (CITA), based on the code described in
 * >>refer	grain	physics	Hansen, J. E., Travis, L. D. 1974, Space Sci. Rev., 16, 527
 *
 * Adapted for Cloudy by Peter A.M. van Hoof (University of Kentucky,
 *   Canadian Institute for Theoretical Astrophysics,
 *   Queen's University of Belfast,
 *   Royal Observatory of Belgium)
 *
 *=======================================================*/


/* these are the magic numbers for the .rfi, .szd, .opc, and .mix files
 * the first digit is file type, the rest is date (YYMMDD) */
static const long MAGIC_RFI = 1030103L;
static const long MAGIC_SZD = 2010403L;
static const long MAGIC_OPC = 3100827L;
static const long MAGIC_MIX = 4030103L;

/* >>chng 02 may 28, by Ryan, moved struct complex to cddefines.h to make it available to entire code. */

/* these are the absolute smallest and largest grain sizes we will
 * consider (in micron). the lower limit gives a grain with on the
 * order of one atom in it, so it is physically motivated. the upper
 * limit comes from the series expansions used in the mie theory,
 * they will have increasingly more problems converging for larger
 * grains, so this limit is numerically motivated */
static const double SMALLEST_GRAIN = 0.0001*(1.-10.*DBL_EPSILON);
static const double LARGEST_GRAIN = 10.*(1.+10.*DBL_EPSILON);

/* maximum no. of parameters for grain size distribution */
static const int NSD = 7;

/* these are the indices into the parameter array a[NSD],
 * NB NB -- the numbers defined below should range from 0 to NSD-1 */
static const int ipSize  = 0; /**< single size */
static const int ipBLo   = 0; /**< lower bound */
static const int ipBHi   = 1; /**< upper bound */
static const int ipExp   = 2; /**< exponent for powerlaw */
static const int ipBeta  = 3; /**< beta parameter for powerlaw */
static const int ipSLo   = 4; /**< scale size for lower exp. cutoff */
static const int ipSHi   = 5; /**< scale size for upper exp. cutoff */
static const int ipAlpha = 6; /**< alpha parameter for exp. cutoff */
static const int ipGCen  = 2; /**< center of gaussian distribution */
static const int ipGSig  = 3; /**< 1-sigma width of gaussian distribution */

/* these are the types of refractive index files we recognize */
typedef enum {
	RFI_TABLE, OPC_TABLE, OPC_GREY, OPC_PAH1, OPC_PAH2N, OPC_PAH2C, OPC_PAH3N, OPC_PAH3C
} rfi_type;

/* these are the types of EMT's we recognize */
typedef enum {
	FARAFONOV00, STOGNIENKO95, BRUGGEMAN35, MAXWELL_GARNETT04
} emt_type;

/* these are all the size distribution cases we support */
typedef enum {
	SD_ILLEGAL, SD_SINGLE_SIZE, SD_POWERLAW, SD_EXP_CUTOFF1, SD_EXP_CUTOFF2,
	SD_EXP_CUTOFF3, SD_LOG_NORMAL, SD_LIN_NORMAL, SD_TABLE, SD_NR_CARBON
} sd_type;

class sd_data {
	void p_clear1()
	{
		xx.clear();
		aa.clear();
		rr.clear();
		ww.clear();
		ln_a.clear();
		ln_a4dNda.clear();
	}
public:
	double a[NSD];            /**< parameters for size distribution */
	double lim[2];            /**< holds lower and upper size limit for entire distribution */
	double clim[2];           /**< holds lower and upper size limit for current bin */
	vector<double> xx;        /**< xx[nn]: abcissas for Gauss quadrature on [-1,1] */
	vector<double> aa;        /**< aa[nn]: weights for Gauss quadrature on [-1,1] */
	vector<double> rr;        /**< rr[nn]: abcissas for Gauss quadrature */
	vector<double> ww;        /**< ww[nn]: weights for Gauss quadrature */
	double unity;             /**< normalization for integrals over size distribution */
	double unity_bin;         /**< normalization for integrals over size distribution */
	double cSize;             /**< the grain size currently being used for the calculations */
	double radius;            /**< average grain radius for current bin <a> */
	double area;              /**< average grain surface area for current bin <4pi*a^2> */
	double vol;               /**< average grain volume for current bin <4pi/3*a^3> */
	vector<double> ln_a;      /**< ln(a)[npts]: log of grain radii for user-supplied size distr */
	vector<double> ln_a4dNda; /**< ln(a^4*dN/da)[npts]: log of user-supplied size distr */
	sd_type sdCase;           /**< SD_SINGLE_SIZE, SD_POWERLAW, ... */
	long int nCarbon;         /**< number of carbon atoms requested in SD_NR_CARBON file */
	long int magic;           /**< magic number */
	long int cPart;           /**< current partition no. for size distribution */
	long int nPart;           /**< total no. of partitions for size distribution */
	long int nmul;            /**< multiplier for obtaining no. of abscissas in gaussian quadrature */
	long int nn;              /**< no. of abscissas used in gaussian quadrature (per partition) */
	long int npts;            /**< no. of points in user-supplied size distr */
	bool lgLogScale;          /**< use logarithmic mesh for integration over size ? */
	void clear()
	{
		p_clear1();
	}
};

/* maximum no. of principal axes for crystalline grains */
static const int NAX = 3;
static const int NDAT = 4;

class grain_data {
	void p_clear0()
	{
		nAxes = 0;
		nOpcCols = 0;
	}
	void p_clear1()
	{
		for( int j=0; j < NAX; j++ ) 
		{
			wavlen[j].clear();
			n[j].clear();
			nr1[j].clear();
		}
		opcAnu.clear();
		for( int j=0; j < NDAT; j++ )
			opcData[j].clear();
	}
public:
	vector<double> wavlen[NAX];      /**< wavelength grid for rfi for all axes (micron) */
	vector< complex<double> > n[NAX];/**< refractive index n for all axes */
	vector<double> nr1[NAX];         /**< re(n)-1 for all axes */
	vector<double> opcAnu;           /**< energies for data points in OPC_TABLE file */
	vector<double> opcData[NDAT];    /**< data values from OPC_TABLE file */
	double wt[NAX];                  /**< relative weight of each axis */
	double abun;                     /**< abundance of grain molecule rel. to hydrogen for max depletion */
	double depl;                     /**< depletion efficiency */
	double elmAbun[LIMELM];          /**< abundances of constituent elements rel. to hydrogen */
	double mol_weight;               /**< molecular weight of grain molecule (amu) */
	double atom_weight;              /**< molecular weight of grain molecule per atom (amu) */
	double rho;                      /**< specific weight (g/cm^3) */
	double norm;                     /**< number of protons in plasma per average grain */
	double work;                     /**< work function (Ryd) */
	double bandgap;                  /**< gap between valence and conduction band (Ryd) */
	double therm_eff;                /**< efficiency of thermionic emission, between 0 and 1 */
	double subl_temp;                /**< sublimation temperature (K) */
	long int magic;                  /**< magic number */
	long int cAxis;                  /**< number of axis currently being used */
	long int nAxes;                  /**< no. of principal axes for this grain */
	long int ndata[NAX];             /**< no. of wavelength points for each axis */
	long int nOpcCols;               /**< no. of data columns in OPC_TABLE file */
	long int nOpcData;               /**< no. of data points in OPC_TABLE file */
	long int charge;                 /**< grain charge, needed for charge dependent cross sections (e) */
	mat_type matType;                /**< material type, determines enthalpy function, etc. */
	rfi_type rfiType;                /**< type of data in rfi file: rfi table, grey grain, pah, etc. */
	void clear()
	{
		p_clear1();
		p_clear0();
	}
	grain_data()
	{
		p_clear0();
	}
};

/* maximum size for grain type labels */
static const int LABELSUB1 = 3;
static const int LABELSUB2 = 5;
static const int LABELSIZE = LABELSUB1 + LABELSUB2 + 4;

/* this is the number of data points used to set up table of optical constants for a mixed grain */
static const long MIX_TABLE_SIZE = 2000L;

STATIC void mie_auxiliary(/*@partial@*/sd_data*,/*@in@*/const grain_data*,/*@in@*/const char*);
STATIC bool mie_auxiliary2(/*@partial@*/vector<int>&,/*@partial@*/multi_arr<double,2>&,
						   /*@partial@*/multi_arr<double,2>&,/*@partial@*/multi_arr<double,2>&,long,long);
STATIC void mie_integrate(/*@partial@*/sd_data*,double,double,/*@out@*/double*);
STATIC void mie_cs_size_distr(double,/*@partial@*/sd_data*,/*@in@*/const grain_data*,
							  void(*)(double,/*@in@*/const sd_data*,/*@in@*/const grain_data*,
									  /*@out@*/double*,/*@out@*/double*,/*@out@*/double*,/*@out@*/int*),
							  /*@out@*/double*,/*@out@*/double*,/*@out@*/double*,/*@out@*/int*);
STATIC void mie_step(double,/*@partial@*/sd_data*,/*@in@*/const grain_data*,
					 void(*)(double,/*@in@*/const sd_data*,/*@in@*/const grain_data*,
							 /*@out@*/double*,/*@out@*/double*,/*@out@*/double*,/*@out@*/int*),
					 /*@partial@*/double*,/*@partial@*/double*,/*@in@*/const double[],/*@out@*/double*,
					 /*@out@*/double*,/*@out@*/double*,/*@out@*/int*);
STATIC void mie_cs(double,/*@in@*/const sd_data*,/*@in@*/const grain_data*,/*@out@*/double*,/*@out@*/double*,
				   /*@out@*/double*,/*@out@*/int*);
STATIC void ld01_fun(/*@in@*/void(*)(double,const sd_data*,const grain_data[],double*,double*,double*,int*),
					 /*@in@*/double,/*@in@*/double,double,/*@in@*/const sd_data*,/*@in@*/const grain_data[],
					 /*@out@*/double*,/*@out@*/double*,/*@out@*/double*,/*@out@*/int*);
inline void car1_fun(double,/*@in@*/const sd_data*,/*@in@*/const grain_data[],/*@out@*/double*,/*@out@*/double*,
					 /*@out@*/double*,/*@out@*/int*);
STATIC void pah1_fun(double,/*@in@*/const sd_data*,/*@in@*/const grain_data*,/*@out@*/double*,/*@out@*/double*,
					 /*@out@*/double*,/*@out@*/int*);
inline void car2_fun(double,/*@in@*/const sd_data*,/*@in@*/const grain_data[],/*@out@*/double*,/*@out@*/double*,
					 /*@out@*/double*,/*@out@*/int*);
STATIC void pah2_fun(double,/*@in@*/const sd_data*,/*@in@*/const grain_data*,/*@out@*/double*,/*@out@*/double*,
					 /*@out@*/double*,/*@out@*/int*);
inline void car3_fun(double,/*@in@*/const sd_data*,/*@in@*/const grain_data[],/*@out@*/double*,/*@out@*/double*,
					 /*@out@*/double*,/*@out@*/int*);
STATIC void pah3_fun(double,/*@in@*/const sd_data*,/*@in@*/const grain_data*,/*@out@*/double*,/*@out@*/double*,
					 /*@out@*/double*,/*@out@*/int*);
inline double Drude(double,double,double,double);
STATIC void tbl_fun(double,/*@in@*/const sd_data*,/*@in@*/const grain_data*,/*@out@*/double*,/*@out@*/double*,
					/*@out@*/double*,/*@out@*/int*);
STATIC double size_distr(double,/*@in@*/const sd_data*);
STATIC double search_limit(double,double,double,sd_data);
STATIC void mie_calc_ial(/*@in@*/const grain_data*,long,/*@out@*/vector<double>&,/*@in@*/const string&,/*@in@*/bool*);
STATIC void mie_repair(/*@in@*/const string&,long,int,int,/*@in@*/const double[],double[],/*@in@*/vector<int>&,
					   bool,/*@in@*/bool*);
STATIC double mie_find_slope(/*@in@*/const double[],/*@in@*/const double[],/*@in@*/const vector<int>&,
							 long,long,int,bool,/*@in@*/bool*);
STATIC void mie_read_ocn(/*@in@*/const string&,/*@out@*/grain_data*);
STATIC void mie_read_rfi(/*@in@*/const string&,/*@out@*/grain_data*);
STATIC void mie_read_mix(/*@in@*/const string&,/*@out@*/grain_data*);
STATIC void init_eps(double,long,/*@in@*/const vector<grain_data>&,/*@out@*/vector< complex<double> >&);
STATIC complex<double> cnewton(
	void(*)(complex<double>,const vector<double>&,const vector< complex<double> >&,
			long,complex<double>*,double*,double*),
	const vector<double>&,const vector< complex<double> >&,long,complex<double>,double);
STATIC void Stognienko(complex<double>,const vector<double>&,const vector< complex<double> >&,
					   long,complex<double>*,double*,double*);
STATIC void Bruggeman(complex<double>,const vector<double>&,const vector< complex<double> >&,
					  long,complex<double>*,double*,double*);
STATIC void mie_read_szd(/*@in@*/const string&,/*@out@*/sd_data*);
STATIC void mie_read_long(/*@in@*/const string&,/*@in@*/const string&,/*@out@*/long int*,bool,long int);
STATIC void mie_read_realnum(/*@in@*/const string&,/*@in@*/const string&,/*@out@*/realnum*,bool,long int);
STATIC void mie_read_double(/*@in@*/const string&,/*@in@*/const string&,/*@out@*/double*,bool,long int);
STATIC void mie_read_form(/*@in@*/const string&,/*@out@*/double[],/*@out@*/double*,/*@out@*/double*);
STATIC void mie_write_form(/*@in@*/const double[],/*@out@*/string&);
STATIC void mie_read_word(/*@in@*/const string&,/*@out@*/string&,bool);
STATIC void mie_next_data(/*@in@*/const string&,/*@in@*/FILE*,/*@out@*/string&,/*@in@*/long int*);
STATIC void mie_next_line(/*@in@*/const string&,/*@in@*/FILE*,/*@out@*/string&,/*@in@*/long int*);

/*=======================================================*/
/* the following five routines are the core of the Mie code supplied by Peter Martin */

STATIC void sinpar(double,double,double,/*@out@*/double*,/*@out@*/double*,/*@out@*/double*,
		   /*@out@*/double*,/*@out@*/double*,/*@out@*/long*);
STATIC void anomal(double,/*@out@*/double*,/*@out@*/double*,/*@out@*/double*,/*@out@*/double*,double,double);
STATIC void bigk(complex<double>,/*@out@*/complex<double>*);
STATIC void ritodf(double,double,/*@out@*/double*,/*@out@*/double*);
STATIC void dftori(/*@out@*/double*,/*@out@*/double*,double,double);

void mie_write_opc(/*@in@*/ const char *rfi_file,
				   /*@in@*/ const char *szd_file,
				   long int nbin)
{
	string chString;

	/* no. of logarithmic intervals in table printout of size distribution function */
	const long NPTS_TABLE = 100L;

	DEBUG_ENTRY( "mie_write_opc()" );

	sd_data sd;

	mie_read_szd( szd_file , &sd );

	sd.nPart = ( sd.sdCase == SD_SINGLE_SIZE || sd.sdCase == SD_NR_CARBON ) ? 1 : nbin;
	if( sd.nPart <= 0 || sd.nPart >= 100 ) 
	{
		fprintf( ioQQQ, " Illegal number of size distribution bins: %ld\n", sd.nPart );
		fprintf( ioQQQ, " The number should be between 1 and 99.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	sd.lgLogScale = true;

	vector<grain_data> gdArr(2);
	grain_data& gd = gdArr[0];
	grain_data& gd2 = gdArr[1];

	mie_read_ocn( rfi_file, &gd );

	if( gd.rfiType == OPC_TABLE && sd.nPart > 1 )
	{
		fprintf( ioQQQ, " Illegal number of size distribution bins: %ld\n", sd.nPart );
		fprintf( ioQQQ, " The number should always be 1 for OPC_TABLE files.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	if( gd.rho <= 0. )
	{
		fprintf( ioQQQ, " Illegal value for the specific density: %.4e\n", gd.rho );
		cdEXIT(EXIT_FAILURE);
	}
	if( gd.mol_weight <= 0. )
	{
		fprintf( ioQQQ, " Illegal value for the molecular weight: %.4e\n", gd.mol_weight );
		cdEXIT(EXIT_FAILURE);
	}

	bool lgWarning = false;

	/* generate output file name from input file names */
	string chFile = rfi_file;
	auto pp = chFile.find( "." );
	chFile.erase(pp);

	string chFile2 = szd_file;
	pp = chFile2.find( "." );
	chFile2.erase(pp);

	if( sd.sdCase != SD_SINGLE_SIZE && sd.sdCase != SD_NR_CARBON ) 
	{
		ostringstream oss;
		oss << setfill('0') << setw(2) << nbin;
		chFile += "_" + chFile2 + "_" + oss.str() + ".opc";
	}
	else 
	{
		chFile += "_" + chFile2 + ".opc";
	}

	mie_auxiliary(&sd,&gd,"init");

	/* number of protons in plasma per average grain volume */
	gd.norm = sd.vol*gd.rho/(ATOMIC_MASS_UNIT*gd.mol_weight*gd.abun*gd.depl);
	double volnorm = sd.vol;
	double volfrac = 1.;

	multi_arr<double,2> acs_abs( sd.nPart, rfield.nflux_with_check );
	multi_arr<double,2> acs_sct( acs_abs.clone() );
	multi_arr<double,2> a1g( acs_abs.clone() );
	vector<double> inv_att_len( rfield.nflux_with_check );

	fprintf( ioQQQ, "\n Starting mie_write_opc, output will be written to %s\n\n", chFile.c_str() );

	FILE *fdes = open_data( chFile, "w" );
	bool lgErr = false;

	time_t timer;
	(void)time(&timer);
	lgErr = lgErr || ( fprintf(fdes,"# this file was created by Cloudy %s (%s) on %s",
				   t_version::Inst().chVersion.c_str(),t_version::Inst().chDate.c_str(),
				   ctime(&timer)) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"# ===========================================\n#\n") < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%12ld # magic number opacity file\n",MAGIC_OPC) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%12ld # magic number rfi/mix file\n",gd.magic) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%12ld # magic number szd file\n",sd.magic) < 0 );

	/* generate grain label for Cloudy output
	 * adjust LABELSIZE in mie.h when the format defined below is changed ! */
	string hlp1 = chFile.substr(0,LABELSUB1+1);
	pp = hlp1.find("-");
	if( pp != string::npos )
		hlp1.erase(pp);

	string hlp2 = chFile2.substr(0,LABELSUB2+1);
	pp = hlp2.find("-");
	if( pp != string::npos )
		hlp2.erase(pp);

	string chGrainLabel = " ";
	if( sd.nPart > 1 ) 
	{
		chGrainLabel += hlp1.substr(0,LABELSUB1) + "-" + hlp2.substr(0,LABELSUB2) + "xx";
		lgErr = lgErr || ( fprintf(fdes,"%-12.12s # grain type label, xx will be replaced by bin no.\n",
					   chGrainLabel.c_str()) < 0 );
	}
	else 
	{
		chGrainLabel += hlp1 + "-" + hlp2;
		lgErr = lgErr || ( fprintf(fdes,"%-12.12s # grain type label\n", chGrainLabel.c_str()) < 0 );
	}

	lgErr = lgErr || ( fprintf(fdes,"%.6e # specific weight (g/cm^3)\n",gd.rho) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # molecular weight of grain molecule (amu)\n",gd.mol_weight) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # average molecular weight per atom (amu)\n", gd.atom_weight) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # abundance of grain molecule at max depletion\n",gd.abun) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # depletion efficiency\n",gd.depl) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # average grain radius <a^3>/<a^2>, full size distr (cm)\n",
			   3.*sd.vol/sd.area) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # average grain surface area <4pi*a^2>, full size distr (cm^2)\n",
			   sd.area) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # average grain volume <4/3pi*a^3>, full size distr (cm^3)\n",
			   sd.vol) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # total grain radius Int(a) per H, full size distr (cm/H)\n",
			   sd.radius/gd.norm) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # total grain area Int(4pi*a^2) per H, full size distr (cm^2/H)\n",
			   sd.area/gd.norm) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # total grain vol Int(4/3pi*a^3) per H, full size distr (cm^3/H)\n",
			   sd.vol/gd.norm) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # work function (Ryd)\n",gd.work) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # gap between valence and conduction band (Ryd)\n",gd.bandgap) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # efficiency of thermionic emission\n",gd.therm_eff) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%.6e # sublimation temperature (K)\n",gd.subl_temp) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%12d # material type, 1=carbonaceous, 2=silicate\n",gd.matType) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"#\n# abundances of constituent elements rel. to hydrogen\n#\n") < 0 );

	for( long nelem=0; nelem < LIMELM; nelem++ ) 
	{
		lgErr = lgErr || ( fprintf(fdes,"%.6e # %s\n",gd.elmAbun[nelem],
					   elementnames.chElementSym[nelem]) < 0 );
	}

	if( sd.sdCase != SD_SINGLE_SIZE && sd.sdCase != SD_NR_CARBON )
	{
		lgErr = lgErr || ( fprintf(fdes,"#\n# entire size distribution, amin=%.5f amax=%.5f micron\n",
					   sd.lim[ipBLo],sd.lim[ipBHi]) < 0 );
		lgErr = lgErr || ( fprintf(fdes,"#\n%.6e # ratio a_max/a_min in each size bin\n",
					   pow(sd.lim[ipBHi]/sd.lim[ipBLo],1./(double)sd.nPart) ) < 0 );
		lgErr = lgErr || ( fprintf(fdes,"#\n# size distribution function\n#\n") < 0 );
		lgErr = lgErr || ( fprintf(fdes,"%12ld # number of table entries\n#\n",NPTS_TABLE+1) < 0 );
		lgErr = lgErr || ( fprintf(fdes,"# ============================\n") < 0 );
		lgErr = lgErr || ( fprintf(fdes,"# size (micr) a^4*dN/da (cm^3/H)\n#\n") < 0 );
		for( long i=0; i <= NPTS_TABLE; i++ )
		{
			double radius, a4dNda;
			radius = sd.lim[ipBLo]*exp((double)i/(double)NPTS_TABLE*log(sd.lim[ipBHi]/sd.lim[ipBLo]));
			radius = max(min(radius,sd.lim[ipBHi]),sd.lim[ipBLo]);
			a4dNda = POW4(radius)*size_distr(radius,&sd)/gd.norm*1.e-12/sd.unity;
			lgErr = lgErr || ( fprintf(fdes,"%.6e %.6e\n",radius,a4dNda) < 0 );
		}
	}
	else
	{
		lgErr = lgErr || ( fprintf(fdes,"#\n") < 0 );
		lgErr = lgErr || ( fprintf(fdes,"%.6e # a_max/a_min = 1 for single sized grain\n", 1. ) < 0 );
		lgErr = lgErr || ( fprintf(fdes,"%12ld # no size distribution table\n",0L) < 0 );
	}

	union {
		double x;
		uint32 i[2];
	} u;

	lgErr = lgErr || ( fprintf(fdes,"#\n") < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%s # check 1\n",rfield.mesh_md5sum().c_str()) < 0 );
	u.x = rfield.emm();
	if( cpu.i().big_endian() )
		lgErr = lgErr || ( fprintf(fdes,"%23.8x %8.8x # check 2\n",u.i[0],u.i[1]) < 0 );
	else
		lgErr = lgErr || ( fprintf(fdes,"%23.8x %8.8x # check 2\n",u.i[1],u.i[0]) < 0 );
	u.x = rfield.egamry();
	if( cpu.i().big_endian() )
		lgErr = lgErr || ( fprintf(fdes,"%23.8x %8.8x # check 3\n",u.i[0],u.i[1]) < 0 );
	else
		lgErr = lgErr || ( fprintf(fdes,"%23.8x %8.8x # check 3\n",u.i[1],u.i[0]) < 0 );
	u.x = rfield.getResolutionScaleFactor();
	if( cpu.i().big_endian() )
		lgErr = lgErr || ( fprintf(fdes,"%23.8x %8.8x # check 4\n",u.i[0],u.i[1]) < 0 );
	else
		lgErr = lgErr || ( fprintf(fdes,"%23.8x %8.8x # check 4\n",u.i[1],u.i[0]) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%32ld # rfield.nflux_with_check\n",rfield.nflux_with_check) < 0 );
	lgErr = lgErr || ( fprintf(fdes,"%32ld # number of size distr. bins\n#\n",sd.nPart) < 0 );

	if( gd.rfiType == OPC_PAH1 )
	{
		gd2.clear();
		mie_read_rfi("graphite.rfi",&gd2);
	}
	else if( gd.rfiType == OPC_PAH2N || gd.rfiType == OPC_PAH2C ||
		 gd.rfiType == OPC_PAH3N || gd.rfiType == OPC_PAH3C )
	{
		gd2.clear();
		mie_read_rfi("gdraine.rfi",&gd2);
	}

	vector<int> ErrorIndex( rfield.nflux_with_check );

	for( long p=0; p < sd.nPart; p++ ) 
	{
		sd.cPart = p;

		mie_auxiliary(&sd,&gd,"step");

		if( sd.nPart > 1 ) 
		{
			/* >>chng 01 mar 20, creating mie_integrate introduced a change in the normalization
			 * of sd.radius, sd.area, and sd.vol; they now give average quantities for this bin.
			 * gd.norm converts average quantities to integrated quantities per H assuming the
			 * number of grains for the entire size distribution, hence multiplication by frac is
			 * needed to convert to the number of grains in this particular size bin, PvH */
			double frac = sd.unity_bin/sd.unity;
			volfrac = sd.vol*frac/volnorm;
			fprintf( ioQQQ, " Starting size bin %ld, amin=%.5f amax=%.5f micron\n",
				 p+1,sd.clim[ipBLo],sd.clim[ipBHi] );
			lgErr = lgErr || ( fprintf(fdes,"# size bin %ld, amin=%.5f amax=%.5f micron\n",
						   p+1,sd.clim[ipBLo],sd.clim[ipBHi]) < 0 );
			lgErr = lgErr || ( fprintf(fdes,"%.6e # average grain ",3.*sd.vol/sd.area) < 0 );
			lgErr = lgErr || ( fprintf(fdes,"radius <a^3>/<a^2>, this bin (cm)\n") < 0 );
			lgErr = lgErr || ( fprintf(fdes,"%.6e # average ",sd.area) < 0 );
			lgErr = lgErr || ( fprintf(fdes,"grain area <4pi*a^2>, this bin (cm^2)\n") < 0 );
			lgErr = lgErr || ( fprintf(fdes,"%.6e # average ",sd.vol) < 0 );
			lgErr = lgErr || ( fprintf(fdes,"grain volume <4/3pi*a^3>, this bin (cm^3)\n") < 0 );
			lgErr = lgErr || ( fprintf(fdes,"%.6e # total grain ",sd.radius*frac/gd.norm) < 0 );
			lgErr = lgErr || ( fprintf(fdes,"radius Int(a) per H, this bin (cm/H)\n") < 0 );
			lgErr = lgErr || ( fprintf(fdes,"%.6e # total grain area ",sd.area*frac/gd.norm) < 0 );
			lgErr = lgErr || ( fprintf(fdes,"Int(4pi*a^2) per H, this bin (cm^2/H)\n") < 0 );
			lgErr = lgErr || ( fprintf(fdes,"%.6e # total grain volume ",sd.vol*frac/gd.norm) < 0 );
			lgErr = lgErr || ( fprintf(fdes,"Int(4/3pi*a^3) per H, this bin (cm^3/H)\n#\n") < 0 );
		}

		bool lgErrorOccurred = false;

		/* calculate the opacity data */
		for( long i=0; i < rfield.nflux_with_check; i++ ) 
		{
			double wavlen = WAVNRYD/rfield.anu(i)*1.e4;

			ErrorIndex[i] = 0;
			acs_abs[p][i] = 0.;
			acs_sct[p][i] = 0.;
			a1g[p][i] = 0.;

			int Error = 0;
			double cosb, cs_abs, cs_sct;

			switch( gd.rfiType )
			{
			case RFI_TABLE:
				for( gd.cAxis=0; gd.cAxis < gd.nAxes; gd.cAxis++ ) 
				{
					mie_cs_size_distr(wavlen,&sd,&gd,mie_cs,&cs_abs,&cs_sct,&cosb,&Error);
					ErrorIndex[i] = max(ErrorIndex[i],Error);
					acs_abs[p][i] += cs_abs*gd.wt[gd.cAxis];
					acs_sct[p][i] += cs_sct*gd.wt[gd.cAxis];
					a1g[p][i] += cs_sct*(1.-cosb)*gd.wt[gd.cAxis];
				}
				lgErrorOccurred = mie_auxiliary2(ErrorIndex,acs_abs,acs_sct,a1g,p,i);
				break;
			case OPC_TABLE:
				gd.cAxis = 0;
				mie_cs_size_distr(wavlen,&sd,&gd,tbl_fun,&cs_abs,&cs_sct,&cosb,&Error);
				ErrorIndex[i] = min(Error,2);
				lgErrorOccurred = lgErrorOccurred || ( Error > 0 );
				acs_abs[p][i] = cs_abs*gd.norm;
				acs_sct[p][i] = cs_sct*gd.norm;
				a1g[p][i] = 1.-cosb;
				break;
			case OPC_GREY:
				ErrorIndex[i] = 0;
				acs_abs[p][i] = 1.3121e-23*volfrac*gd.norm;
				acs_sct[p][i] = 2.6242e-23*volfrac*gd.norm;
				a1g[p][i] = 1.;
				break;
			case OPC_PAH1:
				gd.cAxis = 0;
				for( gd2.cAxis=0; gd2.cAxis < gd2.nAxes; gd2.cAxis++ ) 
				{
					mie_cs_size_distr(wavlen,&sd,&gd,car1_fun,&cs_abs,&cs_sct,&cosb,&Error);
					ErrorIndex[i] = max(ErrorIndex[i],Error);
					acs_abs[p][i] += cs_abs*gd2.wt[gd2.cAxis];
					acs_sct[p][i] += 0.1*cs_abs*gd2.wt[gd2.cAxis];
					a1g[p][i] += 0.1*cs_abs*1.*gd2.wt[gd2.cAxis];
				}
				lgErrorOccurred = mie_auxiliary2(ErrorIndex,acs_abs,acs_sct,a1g,p,i);
				break;
			case OPC_PAH2N:
			case OPC_PAH2C:
				gd.cAxis = 0;
				// any non-zero charge will do in the OPC_PAH2C case
				gd.charge = ( gd.rfiType == OPC_PAH2N ) ? 0 : 1;
				for( gd2.cAxis=0; gd2.cAxis < gd2.nAxes; gd2.cAxis++ ) 
				{
					mie_cs_size_distr(wavlen,&sd,&gd,car2_fun,&cs_abs,&cs_sct,&cosb,&Error);
					ErrorIndex[i] = max(ErrorIndex[i],Error);
					acs_abs[p][i] += cs_abs*gd2.wt[gd2.cAxis];
					acs_sct[p][i] += 0.1*cs_abs*gd2.wt[gd2.cAxis];
					a1g[p][i] += 0.1*cs_abs*1.*gd2.wt[gd2.cAxis];
				}
				lgErrorOccurred = mie_auxiliary2(ErrorIndex,acs_abs,acs_sct,a1g,p,i);
				break;
			case OPC_PAH3N:
			case OPC_PAH3C:
				gd.cAxis = 0;
				// any non-zero charge will do in the OPC_PAH3C case
				gd.charge = ( gd.rfiType == OPC_PAH3N ) ? 0 : 1;
				for( gd2.cAxis=0; gd2.cAxis < gd2.nAxes; gd2.cAxis++ ) 
				{
					mie_cs_size_distr(wavlen,&sd,&gd,car3_fun,&cs_abs,&cs_sct,&cosb,&Error);
					ErrorIndex[i] = max(ErrorIndex[i],Error);
					acs_abs[p][i] += cs_abs*gd2.wt[gd2.cAxis];
					acs_sct[p][i] += 0.1*cs_abs*gd2.wt[gd2.cAxis];
					a1g[p][i] += 0.1*cs_abs*1.*gd2.wt[gd2.cAxis];
				}
				lgErrorOccurred = mie_auxiliary2(ErrorIndex,acs_abs,acs_sct,a1g,p,i);
				break;
			default:
				TotalInsanity();
			}
		}

		/* extrapolate/interpolate for missing data */
		if( lgErrorOccurred ) 
		{
			chString = "absorption cs";
			mie_repair(chString,rfield.nflux_with_check,2,0,rfield.anuptr(),&acs_abs[p][0],ErrorIndex,false,&lgWarning);
			chString = "scattering cs";
			mie_repair(chString,rfield.nflux_with_check,2,1,rfield.anuptr(),&acs_sct[p][0],ErrorIndex,false,&lgWarning);
			chString = "asymmetry parameter";
			mie_repair(chString,rfield.nflux_with_check,1,1,rfield.anuptr(),&a1g[p][0],ErrorIndex,true,&lgWarning);
		}

		for( long i=0; i < rfield.nflux_with_check; i++ ) 
		{
			acs_abs[p][i] /= gd.norm;
			/* >>chng 02 dec 30, do not multiply with (1-g) and write this factor out
			 * separately; this is useful for calculating extinction properties of grains, PvH */
			acs_sct[p][i] /= gd.norm;
		}
	}

	lgErr = lgErr || ( fprintf(fdes,"#\n") < 0 );
	lgErr = lgErr || ( fprintf(fdes,"# ===========================================\n") < 0 );
	lgErr = lgErr || ( fprintf(fdes,"# anu (Ryd) abs_cs_01 (cm^2/H) abs_cs_02.....\n#\n") < 0 );

	for( long i=0; i < rfield.nflux_with_check; i++ ) 
	{
		lgErr = lgErr || ( fprintf(fdes,"%.6e ",rfield.anu(i)) < 0 );
		for( long p=0; p < sd.nPart; p++ ) 
		{
			lgErr = lgErr || ( fprintf(fdes,"%.6e ",acs_abs[p][i]) < 0 );
		}
		lgErr = lgErr || ( fprintf(fdes,"\n") < 0 );
	}

	lgErr = lgErr || ( fprintf(fdes,"#\n") < 0 );
	lgErr = lgErr || ( fprintf(fdes,"# ===========================================\n") < 0 );
	lgErr = lgErr || ( fprintf(fdes,"# anu (Ryd) sct_cs_01 (cm^2/H) sct_cs_02.....\n#\n") < 0 );

	for( long i=0; i < rfield.nflux_with_check; i++ ) 
	{
		lgErr = lgErr || ( fprintf(fdes,"%.6e ",rfield.anu(i)) < 0 );
		for( long p=0; p < sd.nPart; p++ ) 
		{
			lgErr = lgErr || ( fprintf(fdes,"%.6e ",acs_sct[p][i]) < 0 );
		}
		lgErr = lgErr || ( fprintf(fdes,"\n") < 0 );
	}

	lgErr = lgErr || ( fprintf(fdes,"#\n") < 0 );
	lgErr = lgErr || ( fprintf(fdes,"# ===========================================\n") < 0 );
	lgErr = lgErr || ( fprintf(fdes,"# anu (Ryd) (1-g)_bin_01 (1-g)_bin_02.....\n#\n") < 0 );

	for( long i=0; i < rfield.nflux_with_check; i++ ) 
	{
		lgErr = lgErr || ( fprintf(fdes,"%.6e ",rfield.anu(i)) < 0 );
		for( long p=0; p < sd.nPart; p++ ) 
		{
			// cap of 1-g is needed when g is negative...
			lgErr = lgErr || ( fprintf(fdes,"%.6e ",min(a1g[p][i],1.)) < 0 );
		}
		lgErr = lgErr || ( fprintf(fdes,"\n") < 0 );
	}

	fprintf( ioQQQ, " Starting calculation of inverse attenuation length\n" );
	chString = "inverse attenuation length";
	if( gd.rfiType != RFI_TABLE ) 
	{
		/* >>chng 02 sep 18, added graphite case for special files like PAH's, PvH */
		gd2.clear();
		ial_type icase = gv.which_ial[gd.matType];
		switch( icase )
		{
		case IAL_CAR:
			mie_read_rfi("graphite.rfi",&gd2);
			mie_calc_ial(&gd2,rfield.nflux_with_check,inv_att_len,chString,&lgWarning);
			break;
		case IAL_SIL:
			mie_read_rfi("silicate.rfi",&gd2);
			mie_calc_ial(&gd2,rfield.nflux_with_check,inv_att_len,chString,&lgWarning);
			break;
		default:
			fprintf( ioQQQ, " mie_write_opc detected unknown ial type: %d\n" , icase );
			cdEXIT(EXIT_FAILURE);
		}
	}
	else 
	{
		mie_calc_ial(&gd,rfield.nflux_with_check,inv_att_len,chString,&lgWarning);
	}

	lgErr = lgErr || ( fprintf(fdes,"#\n") < 0 );
	lgErr = lgErr || ( fprintf(fdes,"# ===========================================\n") < 0 );
	lgErr = lgErr || ( fprintf(fdes,"# anu (Ryd) inverse attenuation length (cm^-1)\n#\n") < 0 );

	for( long i=0; i < rfield.nflux_with_check; i++ ) 
	{
		lgErr = lgErr || ( fprintf(fdes,"%.6e %.6e\n",rfield.anu(i),inv_att_len[i]) < 0 );
	}

	fclose(fdes);

	if( lgErr ) 
	{
		fprintf( ioQQQ, "\n Error writing file: %s\n", chFile.c_str() );
		if( remove(chFile.c_str()) == 0 )
		{
			fprintf( ioQQQ, " The file has been removed\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	else 
	{
		fprintf( ioQQQ, "\n Opacity file %s written succesfully\n\n", chFile.c_str() );
		if( lgWarning )
		{
			fprintf( ioQQQ, "\n !!! Warnings were detected !!!\n\n" );
		}
	}
	return;
}

STATIC void mie_auxiliary(/*@partial@*/ sd_data *sd,
						  /*@in@*/ const grain_data *gd,
						  /*@in@*/ const char *auxCase)
{
	double amin,
	  amax,
	  delta,
	  oldvol,
	  step;

	/* desired relative accuracy of integration over size distribution */
	const double TOLER = 1.e-3;

	DEBUG_ENTRY( "mie_auxiliary()" );
	if( strcmp(auxCase,"init") == 0 )
	{
		double mass, radius, CpMolecule;
		/* this is the initial estimate for the multiplier needed to get the
		 * number of abscissas in the gaussian quadrature, the correct value
		 * will be iterated below */
		sd->nmul = 1;

		/* calculate average grain surface area and volume over size distribution */
		switch( sd->sdCase ) 
		{
		case SD_SINGLE_SIZE:
			sd->radius = sd->a[ipSize]*1.e-4;
			sd->area = 4.*PI*pow2(sd->a[ipSize])*1.e-8;
			sd->vol = 4./3.*PI*pow3(sd->a[ipSize])*1.e-12;
			break;
		case SD_NR_CARBON:
			if( gd->elmAbun[ipCARBON] == 0. )
			{
				fprintf( ioQQQ, "\n This size distribution can only be combined with"
					 " carbonaceous material, bailing out...\n" );
				cdEXIT(EXIT_FAILURE);
			}
			// calculate number of C atoms per grain molecule
			CpMolecule = gd->elmAbun[ipCARBON]/(gd->abun*gd->depl);
			// now calculate the mass of the whole grain in gram
			mass = (double)sd->nCarbon/CpMolecule*gd->mol_weight*ATOMIC_MASS_UNIT;
			radius = cbrt(3.*mass/(PI4*gd->rho));
			sd->a[ipSize] = radius*1.e4;
			sd->radius = radius;
			sd->area = 4.*PI*pow2(radius);
			sd->vol = 4./3.*PI*pow3(radius);
			break;
		case SD_POWERLAW:
		case SD_EXP_CUTOFF1:
		case SD_EXP_CUTOFF2:
		case SD_EXP_CUTOFF3:
		case SD_LOG_NORMAL:
		case SD_LIN_NORMAL:
		case SD_TABLE:
			/* set up Gaussian quadrature for entire size range,
			 * first estimate no. of abscissas needed */
			amin = sd->lgLogScale ? log(sd->lim[ipBLo]) : sd->lim[ipBLo];
			amax = sd->lgLogScale ? log(sd->lim[ipBHi]) : sd->lim[ipBHi];

			sd->clim[ipBLo] = sd->lim[ipBLo];
			sd->clim[ipBHi] = sd->lim[ipBHi];

			/* iterate nmul until the integrated volume has converged sufficiently */
			oldvol= 0.;
			do 
			{
				sd->nmul *= 2;
				mie_integrate(sd,amin,amax,&sd->unity);
				delta = fabs(sd->vol-oldvol)/sd->vol;
				oldvol = sd->vol;
			} while( sd->nmul <= 1024 && delta > TOLER );

			if( delta > TOLER )
			{
				fprintf( ioQQQ, " could not converge integration of size distribution\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* we can safely reduce nmul by a factor of 2 and
			 * still reach a relative accuracy of TOLER */
			sd->nmul /= 2;
			mie_integrate(sd,amin,amax,&sd->unity);
			break;
		case SD_ILLEGAL:
		default:
			TotalInsanity();
		}
	}
	else if( strcmp(auxCase,"step") == 0 ) 
	{
		/* calculate average grain surface area and volume over size bin */
		switch( sd->sdCase ) 
		{
		case SD_SINGLE_SIZE:
		case SD_NR_CARBON:
			break;
		case SD_POWERLAW:
		case SD_EXP_CUTOFF1:
		case SD_EXP_CUTOFF2:
		case SD_EXP_CUTOFF3:
		case SD_LOG_NORMAL:
		case SD_LIN_NORMAL:
		case SD_TABLE:
			amin = sd->lgLogScale ? log(sd->lim[ipBLo]) : sd->lim[ipBLo];
			amax = sd->lgLogScale ? log(sd->lim[ipBHi]) : sd->lim[ipBHi];
			step = (amax - amin)/(double)sd->nPart;
			amin = amin + (double)sd->cPart*step;
			amax = min(amax,amin + step);

			sd->clim[ipBLo] = sd->lgLogScale ? exp(amin) : amin;
			sd->clim[ipBHi] = sd->lgLogScale ? exp(amax) : amax;

			sd->clim[ipBLo] = max(sd->clim[ipBLo],sd->lim[ipBLo]);
			sd->clim[ipBHi] = min(sd->clim[ipBHi],sd->lim[ipBHi]);

			mie_integrate(sd,amin,amax,&sd->unity_bin);

			break;
		case SD_ILLEGAL:	
		default:
			TotalInsanity();
		}
	}
	else 
		TotalInsanity();
	return;
}

STATIC bool mie_auxiliary2(vector<int>& ErrorIndex,
						   multi_arr<double,2>& acs_abs,
						   multi_arr<double,2>& acs_sct,
						   multi_arr<double,2>& a1g,
						   long p,
						   long i)
{
	DEBUG_ENTRY( "mie_auxiliary2()" );

	bool lgErrorOccurred = false;
	if( ErrorIndex[i] > 0 ) 
	{
		ErrorIndex[i] = min(ErrorIndex[i],2);
		lgErrorOccurred = true;
	}

	switch( ErrorIndex[i] ) 
	{
		/*lint -e616 */
	case 2:
		acs_abs[p][i] = 0.;
		acs_sct[p][i] = 0.;
		/*lint -fallthrough */
		/* controls is supposed to flow to the next case */
	case 1:
		a1g[p][i] = 0.;
		break;
		/*lint +e616 */
	case 0:
		a1g[p][i] /= acs_sct[p][i];
		break;
	default:
		TotalInsanity();
	}

	/* sanity checks */
	if( ErrorIndex[i] < 2 )
		ASSERT( acs_abs[p][i] > 0. && acs_sct[p][i] > 0. );
	if( ErrorIndex[i] < 1 )
		ASSERT( a1g[p][i] > 0. );

	return lgErrorOccurred;
}

STATIC void mie_integrate(/*@partial@*/ sd_data *sd,
						  double amin,
						  double amax,
						  /*@out@*/ double *normalization)
{
	long int j;
	double unity;

	DEBUG_ENTRY( "mie_integrate()" );

	/* set up Gaussian quadrature for size range,
	 * first estimate no. of abscissas needed */
	sd->nn = sd->nmul*((long)(2.*log(sd->clim[ipBHi]/sd->clim[ipBLo])) + 1);
	sd->nn = min(max(sd->nn,2*sd->nmul),4096);
	sd->xx.resize(sd->nn);
	sd->aa.resize(sd->nn);
	sd->rr.resize(sd->nn);
	sd->ww.resize(sd->nn);
	gauss_legendre(sd->nn,sd->xx,sd->aa);
	gauss_init(sd->nn,amin,amax,sd->xx,sd->aa,sd->rr,sd->ww);

	/* now integrate surface area and volume */
	unity = 0.;
	sd->radius = 0.;
	sd->area = 0.;
	sd->vol = 0.;

	for( j=0; j < sd->nn; j++ ) 
	{
		double weight;

		/* use extra factor of size in weights when we use logarithmic mesh */
		if( sd->lgLogScale ) 
		{
			sd->rr[j] = exp(sd->rr[j]);
			sd->ww[j] *= sd->rr[j];
		}
		weight = sd->ww[j]*size_distr(sd->rr[j],sd);
		unity += weight;
		sd->radius += weight*sd->rr[j];
		sd->area += weight*pow2(sd->rr[j]);
		sd->vol += weight*pow3(sd->rr[j]);
	}
	*normalization = unity;
	sd->radius *= 1.e-4/unity;
	sd->area *= 4.*PI*1.e-8/unity;
	sd->vol *= 4./3.*PI*1.e-12/unity;
	return;
}

/* read in the *.opc file with opacities and other relevant information */
void mie_read_opc(/*@in@*/const char *chFile,
				  /*@in@*/const GrainPar& gp)
{
	int res,
	  lgDefaultQHeat;
	long int dl,
	  help,
	  i,
	  nelem,
	  j,
	  magic,
	  nbin,
	  nup;
	size_t nd,
	  nd2;
	realnum RefAbund[LIMELM],
	  VolTotal;
	double anu;
	double RadiusRatio;
	string chLine;
	string md5sum;

	/* if a_max/a_min in a single size bin is less than
	 * RATIO_MAX quantum heating will be turned on by default */
	const double RATIO_MAX = cbrt(100.);

	DEBUG_ENTRY( "mie_read_opc()" );

	FILE *io2 = open_data( chFile, "r" );

	/* include the name of the file we are reading in the Cloudy output */
	chLine.assign( 40, ' ' );
	if( strlen(chFile) <= 40 )
	{
		chLine.replace( 0, strlen(chFile), chFile );
	}
	else
	{
		chLine.replace( 0, 37, chFile );
		chLine.replace( 37, 3, "..." );
	}
	if( called.lgTalk )
		fprintf( ioQQQ, "                       * #<<< mie_read_opc reading file -- %40s >>>> *\n",
			 chLine.c_str() );

	/* >>chng 02 jan 30, check if file has already been read before, PvH */
	for( size_t p=0; p < gv.ReadRecord.size(); ++p )
	{
		if( gv.ReadRecord[p] == chFile )
		{
			fprintf( ioQQQ, " File %s has already been read before, was this intended ?\n", chFile );
			break;
		}
	}
	gv.ReadRecord.push_back( chFile );

	/* allocate memory for first bin */
	gv.bin.emplace_back( GrainBin() );
	nd = gv.bin.size()-1;

	dl = 0; /* line counter for input file */

	/* first read magic numbers */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&magic,true,dl);
	if( magic != MAGIC_OPC ) 
	{
		fprintf( ioQQQ, " Opacity file %s has obsolete magic number\n",chFile );
		fprintf( ioQQQ, " I found magic number %ld, but expected %ld on line #%ld\n",
			 magic,MAGIC_OPC,dl );
		fprintf( ioQQQ, " Please recompile this file\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* the following two magic numbers are for information only */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&magic,true,dl);

	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&magic,true,dl);

	/* the grain scale factor is set equal to the abundances scale factor
	 * that might have appeared on the grains command.  Later, in grains.c, 
	 * it will be further multiplied by gv.GrainMetal, the scale factor that
	 * appears on the metals & grains command.  That command may, or may not,
	 * have been parsed yet, so can't do it at this stage. */
	gv.bin[nd].dstfactor = (realnum)gp.dep;

	/* grain type label */
	mie_next_data(chFile,io2,chLine,&dl);
	strncpy(gv.bin[nd].chDstLab,chLine.c_str(),(size_t)LABELSIZE);
	gv.bin[nd].chDstLab[LABELSIZE] = '\0';

	/* specific weight (g/cm^3) */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].dustp[0],true,dl);
	/* constant needed in the evaluation of the electron escape length */
	gv.bin[nd].eec = pow((double)gv.bin[nd].dustp[0],-0.85);

	/* molecular weight of grain molecule (amu) */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].dustp[1],true,dl);

	/* average molecular weight per atom (amu) */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].atomWeight,true,dl);

	/* abundance of grain molecule for max depletion */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].dustp[2],true,dl);

	/* depletion efficiency */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].dustp[3],true,dl);

	/* fraction of the integrated volume contained in this bin */
	gv.bin[nd].dustp[4] = 1.;

	/* average grain radius <a^3>/<a^2> for entire size distr (cm) */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].AvRadius,true,dl);
	gv.bin[nd].eyc = 1./gv.bin[nd].AvRadius + 1.e7;

	/* average grain area <4pi*a^2> for entire size distr (cm^2) */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].AvArea,true,dl);

	/* average grain volume <4/3pi*a^3> for entire size distr (cm^3) */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].AvVol,true,dl);

	/* total grain radius Int(a) per H for entire size distr (cm/H) */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].IntRadius,true,dl);

	/* total grain area Int(4pi*a^2) per H for entire size distr (cm^2/H) */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].IntArea,true,dl);

	/* total grain vol Int(4/3pi*a^3) per H for entire size distr (cm^3/H) */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].IntVol,true,dl);

	/* work function, in Rydberg */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].DustWorkFcn,true,dl);

	/* bandgap, in Rydberg */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].BandGap,false,dl);

	/* efficiency of thermionic emissions, between 0 and 1 */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].ThermEff,true,dl);

	/* sublimation temperature in K */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_realnum(chFile,chLine,&gv.bin[nd].Tsublimat,true,dl);

	/* material type, determines enthalpy function, etc... */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&help,true,dl);
	gv.bin[nd].matType = (mat_type)help;

	for( nelem=0; nelem < LIMELM; nelem++ ) 
	{
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_realnum(chFile,chLine,&RefAbund[nelem],false,dl);

		gv.bin[nd].elmAbund[nelem] = RefAbund[nelem];

		/* this coefficient is defined at the end of appendix A.10 of BFM */
		gv.bin[nd].AccomCoef[nelem] = 2.*gv.bin[nd].atomWeight*dense.AtomicWeight[nelem]/
			pow2(gv.bin[nd].atomWeight+dense.AtomicWeight[nelem]);
	}

	/* ratio a_max/a_min for grains in a single size bin */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_double(chFile,chLine,&RadiusRatio,true,dl);

	gv.bin[nd].nDustFunc = gp.nDustFunc;
	lgDefaultQHeat = ( RadiusRatio < RATIO_MAX && !gp.lgGreyGrain );
	gv.bin[nd].lgQHeat = ( gp.lgForbidQHeating ) ? false : ( gp.lgRequestQHeating || lgDefaultQHeat );
	gv.bin[nd].cnv_H_pGR = gv.bin[nd].AvVol/gv.bin[nd].IntVol;
	gv.bin[nd].cnv_GR_pH = 1./gv.bin[nd].cnv_H_pGR;

	/* this is capacity per grain, in Farad per grain */
	gv.bin[nd].Capacity = PI4*ELECTRIC_CONST*gv.bin[nd].IntRadius/100.*gv.bin[nd].cnv_H_pGR;

	/* skip the table of the size distribution function (if present) */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&nup,false,dl);
	for( i=0; i < nup; i++ )
		mie_next_data(chFile,io2,chLine,&dl);

	/* read checksum of continuum_mesh.ini */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_word(chLine,md5sum,false);

	union {
		double x;
		uint32 i[2];
	} u;
	double mesh_lo, mesh_hi;

	/* read lower limit of frequency mesh in hex form */
	mie_next_data(chFile,io2,chLine,&dl);
	if( cpu.i().big_endian() )
		sscanf( chLine.c_str(), "%x %x", &u.i[0], &u.i[1] );
	else
		sscanf( chLine.c_str(), "%x %x", &u.i[1], &u.i[0] );
	mesh_lo = u.x;

	/* read upper limit of frequency mesh in hex form */
	mie_next_data(chFile,io2,chLine,&dl);
	if( cpu.i().big_endian() )
		sscanf( chLine.c_str(), "%x %x", &u.i[0], &u.i[1] );
	else
		sscanf( chLine.c_str(), "%x %x", &u.i[1], &u.i[0] );
	mesh_hi = u.x;

	if( md5sum != rfield.mesh_md5sum() ||
	    !fp_equal_tol( mesh_lo, rfield.emm(), 1.e-11*rfield.emm() ) ||
	    !fp_equal_tol( mesh_hi, rfield.egamry(), 1.e-7*rfield.egamry() ) )
	{
		fprintf( ioQQQ, " Opacity file %s has an incompatible energy grid.\n", chFile );
		fprintf( ioQQQ, " Please recompile this file using the COMPILE GRAINS command.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* read mesh resolution scale factor in hex form */
	mie_next_data(chFile,io2,chLine,&dl);
	if( cpu.i().big_endian() )
		sscanf( chLine.c_str(), "%x %x", &u.i[0], &u.i[1] );
	else
		sscanf( chLine.c_str(), "%x %x", &u.i[1], &u.i[0] );
	/* this number is checked later since it may not have been set yet by the input script */
	gv.bin[nd].RSFCheck = u.x;

	/* nup is number of frequency bins stored in file, this should match rfield.nflux_with_check */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&nup,true,dl);

	/* no. of size distribution bins */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&nbin,true,dl);

	/* now update the fields for a resolved size distribution */
	if( nbin > 1 )
	{
		/* remember this number since it will be overwritten below */
		VolTotal = gv.bin[nd].IntVol;

		for( i=0; i < nbin; i++ )
		{
			if( i == 0 )
				nd2 = nd;
			else
			{
				/* allocate memory for remaining bins */
				/* and do a straight copy of all the fields ... */
				gv.bin.emplace_back( GrainBin(gv.bin[nd]) );
				nd2 = gv.bin.size()-1;
			}

			/* ... then update anything that needs updating */

			/* average grain radius <a^3>/<a^2> for this bin (cm) */
			mie_next_data(chFile,io2,chLine,&dl);
			mie_read_realnum(chFile,chLine,&gv.bin[nd2].AvRadius,true,dl);

			/* average grain area in this bin (cm^2) */
			mie_next_data(chFile,io2,chLine,&dl);
			mie_read_realnum(chFile,chLine,&gv.bin[nd2].AvArea,true,dl);

			/* average grain volume in this bin (cm^3) */
			mie_next_data(chFile,io2,chLine,&dl);
			mie_read_realnum(chFile,chLine,&gv.bin[nd2].AvVol,true,dl);

			/* total grain radius Int(a) per H for this bin (cm/H) */
			mie_next_data(chFile,io2,chLine,&dl);
			mie_read_realnum(chFile,chLine,&gv.bin[nd2].IntRadius,true,dl);

			/* total grain area Int(4pi*a^2) per H for this bin (cm^2/H) */
			mie_next_data(chFile,io2,chLine,&dl);
			mie_read_realnum(chFile,chLine,&gv.bin[nd2].IntArea,true,dl);

			/* total grain vol Int(4/3pi*a^3) per H for this bin (cm^3/H) */
			mie_next_data(chFile,io2,chLine,&dl);
			mie_read_realnum(chFile,chLine,&gv.bin[nd2].IntVol,true,dl);

			gv.bin[nd2].cnv_H_pGR = gv.bin[nd2].AvVol/gv.bin[nd2].IntVol;
			gv.bin[nd2].cnv_GR_pH = 1./gv.bin[nd2].cnv_H_pGR;

			/* this is capacity per grain, in Farad per grain */
			gv.bin[nd2].Capacity =
				PI4*ELECTRIC_CONST*gv.bin[nd2].IntRadius/100.*gv.bin[nd2].cnv_H_pGR;

			/* dustp[4] gives the fraction of the grain abundance that is
			 * contained in a particular bin. for unresolved distributions it is
			 * by definition 1, for resolved distributions it is smaller than 1. */
			gv.bin[nd2].dustp[4] = gv.bin[nd2].IntVol/VolTotal;
			for( nelem=0; nelem < LIMELM; nelem++ )
				gv.bin[nd2].elmAbund[nelem] = RefAbund[nelem]*gv.bin[nd2].dustp[4];
		}

		/* this must be in a separate loop! */
		for( i=0; i < nbin; i++ )
		{
			nd2 = nd + i;
			/* modify grain labels */
			char *str = strstr_s(gv.bin[nd2].chDstLab,"xx");
			if( str != NULL )
				sprintf(str,"%02hu",(unsigned short)(i+1));
		}
	}

	/* allocate memory for arrays */
	for( i=0; i < nbin; i++ )
	{
		nd2 = nd + i;
		gv.bin[nd2].dstab1.resize(nup);
		gv.bin[nd2].pure_sc1.resize(nup);
		gv.bin[nd2].asym.resize(nup);
		gv.bin[nd2].dstab1_x_anu.resize(nup);
		gv.bin[nd2].inv_att_len.resize(nup);
	}

	/* skip the next 5 lines */
	for( i=0; i < 5; i++ )
		mie_next_line(chFile,io2,chLine,&dl);

	/* now read absorption opacities */
	for( i=0; i < nup; i++ ) 
	{
		/* read in energy scale and then opacities */
		if( (res = fscanf(io2,"%le",&anu)) != 1 ) 
		{
			fprintf( ioQQQ, " Read failed on %s\n",chFile );
			if( res == EOF )
				fprintf( ioQQQ, " EOF reached prematurely\n" );
			cdEXIT(EXIT_FAILURE);
		}
		// check that frequency grid matches, frequencies are printed with 7 significant digits
		if( !fp_equal_tol(anu,rfield.anu(i),1e-6*rfield.anu(i)) )
		{
			fprintf(ioQQQ,"\n\n PROBLEM while reading frequencies: point %li should "
				"have value %e, but actually has %e\n", (unsigned long)i, rfield.anu(i), anu );
			fprintf(ioQQQ," Please recompile the grain opacity file %s.\n", chFile );
			cdEXIT(EXIT_FAILURE);
		}
		for( j=0; j < nbin; j++ ) 
		{
			nd2 = nd + j;
			if( (res = fscanf(io2,"%le",&gv.bin[nd2].dstab1[i])) != 1 ) 
			{
				fprintf( ioQQQ, " Read failed on %s\n",chFile );
				if( res == EOF )
					fprintf( ioQQQ, " EOF reached prematurely\n" );
				cdEXIT(EXIT_FAILURE);
			}
			ASSERT( gv.bin[nd2].dstab1[i] > 0. );
			gv.bin[nd2].dstab1_x_anu[i] = gv.bin[nd2].dstab1[i]*rfield.anu(i);
		}
	}

	/* skip to end-of-line and then skip next 4 lines */
	for( i=0; i < 5; i++ )
		mie_next_line(chFile,io2,chLine,&dl);

	/* now read scattering opacities */
	for( i=0; i < nup; i++ ) 
	{
		if( (res = fscanf(io2,"%le",&anu)) != 1 ) 
		{
			fprintf( ioQQQ, " Read failed on %s\n",chFile );
			if( res == EOF )
				fprintf( ioQQQ, " EOF reached prematurely\n" );
			cdEXIT(EXIT_FAILURE);
		}
		for( j=0; j < nbin; j++ ) 
		{
			nd2 = nd + j;
			if( (res = fscanf(io2,"%le",&gv.bin[nd2].pure_sc1[i])) != 1 ) 
			{
				fprintf( ioQQQ, " Read failed on %s\n",chFile );
				if( res == EOF )
					fprintf( ioQQQ, " EOF reached prematurely\n" );
				cdEXIT(EXIT_FAILURE);
			}
			ASSERT( gv.bin[nd2].pure_sc1[i] > 0. );
		}
	}

	/* skip to end-of-line and then skip next 4 lines */
	for( i=0; i < 5; i++ )
		mie_next_line(chFile,io2,chLine,&dl);

	/* now read asymmetry factor */
	for( i=0; i < nup; i++ ) 
	{
		if( (res = fscanf(io2,"%le",&anu)) != 1 ) 
		{
			fprintf( ioQQQ, " Read failed on %s\n",chFile );
			if( res == EOF )
				fprintf( ioQQQ, " EOF reached prematurely\n" );
			cdEXIT(EXIT_FAILURE);
		}
		for( j=0; j < nbin; j++ ) 
		{
			nd2 = nd + j;
			if( (res = fscanf(io2,"%le",&gv.bin[nd2].asym[i])) != 1 ) 
			{
				fprintf( ioQQQ, " Read failed on %s\n",chFile );
				if( res == EOF )
					fprintf( ioQQQ, " EOF reached prematurely\n" );
				cdEXIT(EXIT_FAILURE);
			}
			ASSERT( gv.bin[nd2].asym[i] > 0. );
			// just in case we read an old opacity file...
			gv.bin[nd2].asym[i] = min(gv.bin[nd2].asym[i],1.);
		}
	}

	/* skip to end-of-line and then skip next 4 lines */
	for( i=0; i < 5; i++ )
		mie_next_line(chFile,io2,chLine,&dl);

	/* now read inverse attenuation length */
	for( i=0; i < nup; i++ ) 
	{
		double help;
		if( (res = fscanf(io2,"%le %le",&anu,&help)) != 2 ) 
		{
			fprintf( ioQQQ, " Read failed on %s\n",chFile );
			if( res == EOF )
				fprintf( ioQQQ, " EOF reached prematurely\n" );
			cdEXIT(EXIT_FAILURE);
		}
		gv.bin[nd].inv_att_len[i] = (realnum)help;
		ASSERT( gv.bin[nd].inv_att_len[i] > 0. );

		for( j=1; j < nbin; j++ ) 
		{
			nd2 = nd + j;
			gv.bin[nd2].inv_att_len[i] = gv.bin[nd].inv_att_len[i];
		}
	}

	fclose(io2);
	return;
}

/* calculate average absorption, scattering cross section (i.e. pi a^2 Q) and
 * average asymmetry parameter g for an arbitrary grain size distribution */
STATIC void mie_cs_size_distr(double wavlen, /* micron */
							  /*@partial@*/ sd_data *sd,
							  /*@in@*/ const grain_data *gd,
							  void(*cs_fun)(double,/*@in@*/const sd_data*,
											/*@in@*/const grain_data*,
											/*@out@*/double*,/*@out@*/double*,
											/*@out@*/double*,/*@out@*/int*),
							  /*@out@*/ double *cs_abs, /* cm^2, average */
							  /*@out@*/ double *cs_sct, /* cm^2, average */
							  /*@out@*/ double *cosb,
							  /*@out@*/ int *error)
{
	DEBUG_ENTRY( "mie_cs_size_distr()" );

	/* sanity checks */
	ASSERT( wavlen > 0. );
	ASSERT( gd->cAxis >= 0 && gd->cAxis < gd->nAxes && gd->cAxis < NAX );

	bool lgADLused;
	const double TOLER = 1.e-3;
	double rr, h, absval, sctval, mycosb, toler[3];

	switch( sd->sdCase ) 
	{
	case SD_SINGLE_SIZE:
	case SD_NR_CARBON:
		/* do single sized grain */
		ASSERT( sd->a[ipSize] > 0. );
		sd->cSize = sd->a[ipSize];
		(*cs_fun)(wavlen,sd,gd,cs_abs,cs_sct,cosb,error);
		break;
	case SD_POWERLAW:
		/* simple powerlaw distribution */
	case SD_EXP_CUTOFF1:
	case SD_EXP_CUTOFF2:
	case SD_EXP_CUTOFF3:
		/* powerlaw distribution with exponential cutoff */
	case SD_LOG_NORMAL:
		/* gaussian distribution in ln(a) */
	case SD_LIN_NORMAL:
		/* gaussian distribution in a */
	case SD_TABLE:
		/* user supplied table of a^4*dN/da */
		ASSERT( sd->lim[ipBLo] > 0. && sd->lim[ipBHi] > 0. && sd->lim[ipBHi] > sd->lim[ipBLo] );
		lgADLused = false;
		*cs_abs = 0.;
		*cs_sct = 0.;
		*cosb = 0.;
		*error = 0;
		// integrate from the upper end downwards since the upper end contributes most
		rr = sd->clim[ipBHi];
		h = min(0.03*rr*sqrt(TOLER),sd->clim[ipBHi]-sd->clim[ipBLo]);
		toler[0] = -TOLER;
		toler[1] = -TOLER;
		toler[2] = -TOLER;
		do
		{
			mie_step(wavlen,sd,gd,cs_fun,&rr,&h,toler,&absval,&sctval,&mycosb,error);
			if( *error >= 2 ) 
			{
				/* mie_cs or mie_step failed to converge -> integration is invalid */
				*cs_abs = -1.;
				*cs_sct = -1.;
				*cosb = -2.;
				return;
			}
			else if( *error == 1 ) 
			{
				/* anomalous diffraction limit used -> g is not valid */
				lgADLused = true;
			}
			*cs_abs += absval;
			*cs_sct += sctval;
			*cosb += mycosb;
			toler[0] = (*cs_abs)*TOLER;
			toler[1] = (*cs_sct)*TOLER;
			toler[2] = abs(*cosb)*TOLER;
			if( rr-h < sd->clim[ipBLo] )
				h = rr-sd->clim[ipBLo];
		}
		while( !fp_equal( rr, sd->clim[ipBLo] ) );
		if( lgADLused )
		{
			*error = 1;
			*cosb = -2.;
		}
		else 
		{
			*error = 0;
			*cosb /= *cs_sct;
		}
		*cs_abs /= sd->unity;
		*cs_sct /= sd->unity;
		break;
	case SD_ILLEGAL:
	default:
		TotalInsanity();
	}
	/* sanity checks */
	if( *error < 2 )
		ASSERT( *cs_abs > 0. && *cs_sct > 0. );
	if( *error < 1 )
		ASSERT( fabs(*cosb) <= 1.+10.*DBL_EPSILON );
	return;
}

STATIC void mie_step(double wavlen, /* micron */
					 /*@partial@*/ sd_data* sd,
					 /*@in@*/ const grain_data* gd,
					 void(*cs_fun)(double,/*@in@*/const sd_data*,
								   /*@in@*/const grain_data*,
								   /*@out@*/double*,/*@out@*/double*,
								   /*@out@*/double*,/*@out@*/int*),
					 /*@partial@*/ double* x,
					 /*@partial@*/ double* h,
					 /*@in@*/ const double toler[], /* toler[3] */
					 /*@out@*/ double* absval,
					 /*@out@*/ double* sctval,
					 /*@out@*/ double* cosb,
					 /*@out@*/ int* error)
{
	DEBUG_ENTRY( "mie_step()" );

	const int MAX_ITER = 8;
	const int MAX_DIVS = 3;
	const int MAX_ABS = (1<<MAX_DIVS) + 1;

	double y1[MAX_DIVS+1], y2[MAX_DIVS+1], y3[MAX_DIVS+1];
	double y1a[MAX_ABS], y2a[MAX_ABS], y3a[MAX_ABS];

	for( int k=0; k < MAX_ITER; ++k )
	{
		for( int i=0; i <= MAX_DIVS; ++i )
		{
			int nabs = 1<<i;
			int stride = 1<<(MAX_DIVS-i);
			double h1 = (*h)/double(nabs);
			// for k > 0 the outer points at 0 and MAX_DIVS are already set up...
			if( k == 0 || i > 0 )
			{
				for( int j=0; j <= nabs; ++j )
				{
					int index = j*stride;
					// y1a[index], etc, may already have been initialized on a previous
					// iteration of the loop over i. If so, do not repeat the calculation
					if( i == 0 || j%2 == 1 )
					{
						sd->cSize = max(*x - h1*double(j),sd->clim[ipBLo]);
						int myerror;
						double myabsval, mysctval, myg;
						(*cs_fun)(wavlen,sd,gd,&myabsval,&mysctval,&myg,&myerror);
						double weight = size_distr(sd->cSize,sd);
						y1a[index] = weight*myabsval;
						y2a[index] = weight*mysctval;
						y3a[index] = weight*mysctval*myg;
						*error = max(*error,myerror);
						if( *error >= 2 )
							return;
					}
				}
			}
			// use trapezium rule
			y1[i] = (y1a[0]+y1a[nabs*stride])/2.;
			y2[i] = (y2a[0]+y2a[nabs*stride])/2.;
			y3[i] = (y3a[0]+y3a[nabs*stride])/2.;
			for( int j=1; j < nabs; ++j )
			{
				y1[i] += y1a[j*stride];
				y2[i] += y2a[j*stride];
				y3[i] += y3a[j*stride];
			}
			y1[i] *= h1;
			y2[i] *= h1;
			y3[i] *= h1;
			// check for convergence by comparing the integrals with stepsizes h1 and 2*h1
			if( i > 0 )
			{
				double err;
				// the asymmetry parameter requires special treatment since it can be
				// quite noisy due to cancellation error when it is very close to zero
				const double TINY_G = 1.e-15;
				if( toler[0] < 0. )
				{
					// -toler is relative precision
					double e1 = safe_div(abs((y1[i]-y1[i-1])/toler[0]),y1[i],0.);
					double e2 = safe_div(abs((y2[i]-y2[i-1])/toler[1]),y2[i],0.);
					err = max(e1,e2);
					if( *error == 0 )
					{
						double e3 = abs(y3[i]-y3[i-1])/max(abs(toler[2]*y3[i]),TINY_G);
						err = max(err,e3);
					}
				}
				else
				{
					// toler is absolute precision
					double e1 = abs(y1[i]-y1[i-1])/toler[0];
					double e2 = abs(y2[i]-y2[i-1])/toler[1];
					err = max(e1,e2);
					if( *error == 0 )
					{
						double e3 = abs(y3[i]-y3[i-1])/max(toler[2],TINY_G);
						err = max(err,e3);
					}
				}
				if( err <= 1. )
				{
					*absval = y1[i];
					*sctval = y2[i];
					*cosb = y3[i];
					double fac = 0.9*sqrt(1./max(err,0.1));
					*x -= *h;
					*h = 2.*h1*fac;
					return;
				}
			}
		}
		// convergence failed (e.g. due to the presence of a resonance in
		// the region [x,x+h]). Here we reduce the stepsize and try again...
		y1a[MAX_ABS-1] = y1a[1];
		y2a[MAX_ABS-1] = y2a[1];
		y3a[MAX_ABS-1] = y3a[1];
		*h /= double(MAX_ABS-1);
	}
	if( false )
	{
		fprintf(ioQQQ, "=================================================\n");
		for( int j=0; j <= MAX_DIVS; ++j )
			fprintf(ioQQQ, "%d %e %e %e\n",j,y1[j],y2[j],y3[j]);
		fprintf(ioQQQ, "\n");
		for( int j=0; j <= (1<<MAX_DIVS); ++j )
			fprintf(ioQQQ, "%d %e %e %e\n",j,y1a[j],y2a[j],y3a[j]);
		fprintf(ioQQQ, "=================================================\n");
	}
	// no convergence -> declare failure
	*error = 2;
	return;
}

/* calculate absorption, scattering cross section (i.e. pi a^2 Q) and
 * asymmetry parameter g (=cosb) for a single sized grain defined by gd */
STATIC void mie_cs(double wavlen,  /* micron */
				   /*@in@*/ const sd_data *sd,
				   /*@in@*/ const grain_data *gd,
				   /*@out@*/ double *cs_abs, /* cm^2 */
				   /*@out@*/ double *cs_sct, /* cm^2 */
				   /*@out@*/ double *cosb,
				   /*@out@*/ int *error)
{
	bool lgOutOfBounds;
	long int iflag,
	  ind;
	double area,
	  aqabs,
	  aqext,
	  aqphas,
	  beta,
	  ctbrqs = -DBL_MAX,
	  delta,
	  frac,
	  nim,
	  nre,
	  nr1,
	  qback,
	  qext = -DBL_MAX,
	  qphase,
	  qscatt = -DBL_MAX,
	  x,
	  xistar;

	DEBUG_ENTRY( "mie_cs()" );

	/* sanity checks, should already have been checked further upstream */
	ASSERT( wavlen > 0. );
	ASSERT( sd->cSize > 0. );
	ASSERT( gd->ndata[gd->cAxis] > 1 && (long)gd->wavlen[gd->cAxis].size() == gd->ndata[gd->cAxis] );

	/* first interpolate optical constants */
	/* >>chng 02 oct 22, moved calculation of optical constants from mie_cs_size_distr to mie_cs,
	 * this increases overhead, but makes the code in mie_cs_size_distr more transparent, PvH */
	find_arr(wavlen,gd->wavlen[gd->cAxis],gd->ndata[gd->cAxis],&ind,&lgOutOfBounds);

	if( lgOutOfBounds ) 
	{
		*error = 3;
		*cs_abs = -1.;
		*cs_sct = -1.;
		*cosb = -2.;
		return;
	}

	frac = (wavlen-gd->wavlen[gd->cAxis][ind])/(gd->wavlen[gd->cAxis][ind+1]-gd->wavlen[gd->cAxis][ind]);
	ASSERT( frac > 0.-10.*DBL_EPSILON && frac < 1.+10.*DBL_EPSILON );
	nre = (1.-frac)*gd->n[gd->cAxis][ind].real() + frac*gd->n[gd->cAxis][ind+1].real();
	ASSERT( nre > 0. );
	nim = (1.-frac)*gd->n[gd->cAxis][ind].imag() + frac*gd->n[gd->cAxis][ind+1].imag();
	ASSERT( nim > 0. );
	nr1 = (1.-frac)*gd->nr1[gd->cAxis][ind] + frac*gd->nr1[gd->cAxis][ind+1];
	ASSERT( fabs(nre-1.-nr1) < 10.*max(nre,1.)*DBL_EPSILON );

	/* size in micron, area in cm^2 */
	area = PI*pow2(sd->cSize)*1.e-8;

	x = sd->cSize/wavlen*2.*PI;

	/* note that in the following, only nre, nim are used in sinpar
	 * and also nr1 in anomalous diffraction limit */

	sinpar(nre,nim,x,&qext,&qphase,&qscatt,&ctbrqs,&qback,&iflag);

	/* iflag=0 normal exit, 1 failure to converge, 2 not even tried
	 * for exit 1,2, see whether anomalous diffraction is available */

	if( iflag == 0 ) 
	{
		*error = 0;
		*cs_abs = area*(qext - qscatt);
		*cs_sct = area*qscatt;
		*cosb = ctbrqs/qscatt;
	}
	else 
	{
		/* anomalous diffraction -- x >> 1 and |m-1| << 1 but any phase shift */
		if( x >= 100. && sqrt(nr1*nr1+nim*nim) <= 0.001 ) 
		{
			delta = -nr1;
			beta = nim;

			anomal(x,&aqext,&aqabs,&aqphas,&xistar,delta,beta);

			/* cosb is invalid */
			*error = 1;
			*cs_abs = area*aqabs;
			*cs_sct = area*(aqext - aqabs);
			*cosb = -2.;
		}
		/* nothing works */
		else 
		{
			*error = 2;
			*cs_abs = -1.;
			*cs_sct = -1.;
			*cosb = -2.;
		}
	}
	if( *error < 2 ) 
	{
		if( *cs_abs <= 0. || *cs_sct <= 0. ) 
		{
			fprintf( ioQQQ, " illegal opacity found: wavl=%.4e micron," , wavlen );
			fprintf( ioQQQ, " abs_cs=%.2e, sct_cs=%.2e\n" , *cs_abs , *cs_sct );
			fprintf( ioQQQ, " please check refractive index file...\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	if( *error < 1 ) 
	{
		if( fabs(*cosb) > 1.+10.*DBL_EPSILON ) 
		{
			fprintf( ioQQQ, " illegal asymmetry parameter found: wavl=%.4e micron," , wavlen );
			fprintf( ioQQQ, " cosb=%.2e\n" , *cosb );
			fprintf( ioQQQ, " please check refractive index file...\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	return;
}

/* this routine calculates the absorption cross sections of carbonaceous grains, it is based on Eq. 2 of:
 * >>refer	grain	physics	Li, A., & Draine, B.T. 2001 ApJ, 554, 778 */
STATIC void ld01_fun(/*@in@*/ void(*pah_fun)(double,const sd_data*,const grain_data[],double*,double*,double*,int*),
					 /*@in@*/ double q_gra,   /* defined in LD01 */
					 /*@in@*/ double wmin,    /* below wmin use pure graphite */
					 /*@in@*/ double wavl,    /* in micron */
					 /*@in@*/ const sd_data *sd,
					 /*@in@*/ const grain_data gdArr[], /* gdArr[2] */
					 /*@out@*/ double *cs_abs,
					 /*@out@*/ double *cs_sct,
					 /*@out@*/ double *cosb,
					 /*@out@*/ int *error)
{
	DEBUG_ENTRY( "ld01_fun()" );

	// this implements Eqs. 2 & 3 of LD01; it creates a gradual change from PAH-like behavior at
	// small radii (less than 50A) to graphite-like at large radii. The factor q_gra is non-zero
	// to include a small amount of "continuum" absorption even for very small grains.

	const double a_xi = 50.e-4;
	double xi_PAH, cs_abs1, cs_abs2;
	if( wavl >= wmin )
	{
		(*pah_fun)(wavl,sd,&gdArr[0],&cs_abs1,cs_sct,cosb,error);
		xi_PAH = (1.-q_gra)*min(1.,pow3(a_xi/sd->cSize));
	}
	else
	{
		cs_abs1 = 0.;
		xi_PAH = 0.;
	}
	// ignore cs_sct, cosb, and error from pah2_fun and return the graphite ones instead.
	// pah2_fun never returns errors and the other two values are ignored by the upstream code
	mie_cs(wavl,sd,&gdArr[1],&cs_abs2,cs_sct,cosb,error);
	*cs_abs = xi_PAH*cs_abs1 + (1.-xi_PAH)*cs_abs2;
	return;
}

/* this routine calculates the absorption cross sections of carbonaceous grains, it creates a gradual
 * change from pah1_fun defined below at small grain radii to graphite-like behavior at large radii */
inline void car1_fun(double wavl,    /* in micron */
					 /*@in@*/ const sd_data *sd,
					 /*@in@*/ const grain_data gdArr[], /* gdArr[2] */
					 /*@out@*/ double *cs_abs,
					 /*@out@*/ double *cs_sct,
					 /*@out@*/ double *cosb,
					 /*@out@*/ int *error)
{
	ld01_fun(pah1_fun,0.,0.,wavl,sd,gdArr,cs_abs,cs_sct,cosb,error);
}

/* this routine calculates the absorption cross sections of PAH molecules, it is based on:
 * >>refer	grain	physics	Desert, F.-X., Boulanger, F., Puget, J. L. 1990, A&A, 237, 215
 *
 * the original version of this routine was written by Kevin Volk (University Of Calgary) */

static const double pah1_strength[7] = { 1.4e-21,1.8e-21,1.2e-20,6.0e-21,4.0e-20,1.9e-20,1.9e-20 };
static const double pah1_wlBand[7] = { 3.3, 6.18, 7.7, 8.6, 11.3, 12.0, 13.3 };
static const double pah1_width[7] = { 0.024, 0.102, 0.24, 0.168, 0.086, 0.174, 0.174 };

STATIC void pah1_fun(double wavl,    /* in micron */
					 /*@in@*/ const sd_data *sd,
					 /*@in@*/ const grain_data *gd,
					 /*@out@*/ double *cs_abs,
					 /*@out@*/ double *cs_sct,
					 /*@out@*/ double *cosb,
					 /*@out@*/ int *error)
{
	long int j;
	double cval1,
	  pah1_fun_v,
	  term,
	  term1,
	  term2,
	  term3,
	  x;

	const double p1 = 4.0e-18;
	const double p2 = 1.1e-18;
	const double p3 = 3.3e-21;
	const double p4 = 6.0e-21;
	const double p5 = 2.4e-21;
	const double wl1a = 5.0;
	const double wl1b = 7.0;
	const double wl1c = 9.0;
	const double wl1d = 9.5;
	const double wl2a = 11.0;
	const double delwl2 = 0.3;
	/* this is the rise interval for the second plateau */
	const double wl2b = wl2a + delwl2;
	const double wl2c = 15.0;
	const double wl3a = 3.2;
	const double wl3b = 3.57;
	const double wl3m = (wl3a+wl3b)/2.;
	const double wl3sig = 0.1476;
	const double x1 = 4.0;
	const double x2 = 5.9;
	const double x2a = RYD_INF/1.e4;
	const double x3 = 0.1;

	/* grain volume */
	double vol = 4.*PI/3.*pow3(sd->cSize)*1.e-12;
	/* number of carbon atoms in PAH molecule */
	double xnc = floor(vol*gd->rho/(ATOMIC_MASS_UNIT*dense.AtomicWeight[ipCARBON]));
	/* number of hydrogen atoms in PAH molecule */
	/* >>chng 02 oct 18, use integral number of hydrogen atoms instead of fractional number */
	double xnh = floor(sqrt(6.*xnc));
	/* this is the hydrogen over carbon ratio in the PAH molecule */
	double xnhoc = xnh/xnc;
	/* ftoc3p3 is the feature to continuum ratio at 3.3 micron */
	double ftoc3p3 = 100.;

	double csVal1 = 0.;
	double csVal2 = 0.;

	DEBUG_ENTRY( "pah1_fun()" );

	x = 1./wavl;

	if( x >= x2a )
	{
		/* >>chng 02 oct 18, use atomic cross sections for energies larger than 1 Ryd */
		double anu_ev = x/x2a*EVRYD;

		/* use Hartree-Slater cross sections */
		t_ADfA::Inst().set_version( PHFIT95 );

		term1 = t_ADfA::Inst().phfit(ipHYDROGEN+1,1,1,anu_ev);
		term2 = 0.;
		for( j=1; j <= 3; ++j )
			term2 += t_ADfA::Inst().phfit(ipCARBON+1,6,j,anu_ev);

		csVal1 = (xnh*term1 + xnc*term2)*1.e-18;
	}

	if( x <= 2.*x2a )
	{
		cval1 = log(sqrt(xnc)*ftoc3p3/1.2328)/12.2;

		term = pow2(min(x,x1))*(3.*x1 - 2.*min(x,x1))/pow3(x1);

		term1 = (0.1*x + 0.41)*pow2(max(x-x2,0.));

		/* The following is an exponential cut-off in the continuum for 
		 * wavelengths larger than 2500 Angstroms; it is exponential in 
		 * wavelength so it varies as 1/x here.  This replaces the 
		 * sharp cut-off at 8000 Angstroms in the original paper.
		 *
		 * This choice of continuum shape is also arbitrary.  The continuum
		 * is never observed at these wavelengths.  For the "standard" ratio 
		 * value at 3.3 microns the continuum level in the optical is not that 
		 * much smaller than it was in the original paper.  If one wants to 
		 * exactly reproduce the original optical opacity, one can change 
		 * the x1 value to a value of 1.125.  Then there will be a discontinuity
		 * in the cross-section at 8000 Angstroms.
		 *
		 * My judgement was that the flat cross-section in the optical used by 
		 * Desert, Boulander, and Puget (1990) is just a rough value that is not 
		 * based upon much in the way of direct observations, and so I could 
		 * change the cross-section at wavelengths above 2500 Angstroms.  It is
		 * likely that one should really build in the ERE somehow, judging from 
		 * the spectrum of the Red Rectangle, and there is no trace of this in
		 * the original paper.  The main concern in adding this exponential 
		 * drop-off in the continuum was to have a finite infrared continuum 
		 * between the features. */
		term2 = exp(cval1*(1. - (x1/min(x,x1))));

		term3 = p3*exp(-pow2(x/x3))*min(x,x3)/x3;

		csVal2 = xnc*((p1*term + p2*term1)*term2 + term3);
	}

	if( x2a <= x && x <= 2.*x2a )
	{
		/* create gradual change from Desert et al to atomic cross sections */
		double frac = pow2(2.-x/x2a);
		pah1_fun_v = exp((1.-frac)*log(csVal1) + frac*log(csVal2));
	}
	else
	{
		/* one of these will be zero */
		pah1_fun_v = csVal1 + csVal2;
	}

	/* now add in the three plateau features. the first two are based upon
	 * >>refer	grain	physics	Schutte, Tielens, and Allamandola (1993) ApJ, 415, 397. */
	if( wl1a <= wavl && wl1d >= wavl )
	{
		if( wavl < wl1b )
			term = p4*(wavl - wl1a)/(wl1b - wl1a);
		else
			term = ( wavl > wl1c ) ? p4*(wl1d - wavl)/(wl1d - wl1c) : p4;
		pah1_fun_v += term*xnc;
	}
	if( wl2a <= wavl && wl2c >= wavl )
	{
		term = ( wavl < wl2b ) ? p5*((wavl - wl2a)/delwl2) : p5*sqrt(1.-pow2((wavl-wl2a)/(wl2c-wl2a)));
		pah1_fun_v += term*xnc;
	}
	if( wl3m-10.*wl3sig <= wavl && wavl <= wl3m+10.*wl3sig )
	{
		/* >>chng 02 nov 08, replace top hat distribution with gaussian, PvH */
		term = 1.1*pah1_strength[0]*exp(-0.5*pow2((wavl-wl3m)/wl3sig));
		pah1_fun_v += term*xnh;
	}

	/* add in the various discrete features in the infrared: 3.3, 6.2, 7.6, etc. */
	for( j=0; j < 7; j++ )
	{
		term1 = (wavl - pah1_wlBand[j])/pah1_width[j];
		term = 0.;
		if( j == 2 )
		{
			/* This assumes linear interpolation between the points, which are
			 * located at -1, -0.5, +1.5, and +3 times the width, or a fine spacing that
			 * well samples the profile. Otherwise there will be an error in the total
			 * feature strength of order 50%  */
			if( term1 >= -1. && term1 < -0.5 )
			{
				term = pah1_strength[j]/(3.*pah1_width[j]);
				term *= 2. + 2.*term1;
			}
			if( term1 >= -0.5 && term1 <= 1.5 )
				term = pah1_strength[j]/(3.*pah1_width[j]);
			if( term1 > 1.5 && term1 <= 3. )
			{
				term = pah1_strength[j]/(3.*pah1_width[j]);
				term *= 2. - term1*2./3.;
			}
		}
		else
		{
			/* This assumes linear interpolation between the points, which are
			 * located at -2, -1, +1, and +2 times the width, or a fine spacing that
			 * well samples the profile. Otherwise there will be an error in the total
			 * feature strength of order 50%  */
			if( term1 >= -2. && term1 < -1. )
			{
				term = pah1_strength[j]/(3.*pah1_width[j]);
				term *= 2. + term1;
			}
			if( term1 >= -1. && term1 <= 1. )
				term = pah1_strength[j]/(3.*pah1_width[j]);
			if( term1 > 1. && term1 <= 2. )
			{
				term = pah1_strength[j]/(3.*pah1_width[j]);
				term *= 2. - term1;
			}
		}
		if( j == 0 || j > 2 )
			term *= xnhoc;
		pah1_fun_v += term*xnc;
	}

	*cs_abs = pah1_fun_v;
	/* the next two numbers are completely arbitrary, but the code requires them... */
	/* >>chng 02 oct 18, cs_sct was 1.e-30, but this is very high for X-ray photons, PvH */
	*cs_sct = 0.1*pah1_fun_v;
	*cosb = 0.;
	*error = 0;

	return;
}

/* this routine calculates the absorption cross sections of carbonaceous grains, it is based on Eq. 2 of:
 * >>refer	grain	physics	Li, A., & Draine, B.T. 2001 ApJ, 554, 778 */
inline void car2_fun(double wavl,    /* in micron */
					 /*@in@*/ const sd_data *sd,
					 /*@in@*/ const grain_data gdArr[], /* gdArr[2] */
					 /*@out@*/ double *cs_abs,
					 /*@out@*/ double *cs_sct,
					 /*@out@*/ double *cosb,
					 /*@out@*/ int *error)
{
	ld01_fun(pah2_fun,0.01,1./17.25,wavl,sd,gdArr,cs_abs,cs_sct,cosb,error);
}

// these values are taken from Table 1 of LD01
static const double pah2_wavl[14] = { 0.0722, 0.2175, 3.3, 6.2, 7.7, 8.6, 11.3,
				      11.9, 12.7, 16.4, 18.3, 21.2, 23.1, 26.0 };
static const double pah2_width[14] = { 0.195, 0.217, 0.012, 0.032, 0.091, 0.047, 0.018,
				       0.025, 0.024, 0.010, 0.036, 0.038, 0.046, 0.69 };
static const double pah2n_strength[14] = { 7.97e-13, 1.23e-13, 1.97e-18, 1.96e-19, 6.09e-19, 3.47e-19, 4.27e-18,
					   7.27e-19, 1.67e-18, 5.52e-20, 6.04e-20, 1.08e-19, 2.78e-20, 1.52e-19 };
static const double pah2c_strength[14] = { 7.97e-13, 1.23e-13, 4.47e-19, 1.57e-18, 5.48e-18, 2.42e-18, 4.00e-18,
					   6.14e-19, 1.49e-18, 5.52e-20, 6.04e-20, 1.08e-19, 2.78e-20, 1.52e-19 };

/* this routine calculates the absorption cross sections of PAH molecules, it is based on Eqs. 6-11 of:
 * >>refer	grain	physics	Li, A., & Draine, B.T. 2001 ApJ, 554, 778 */
STATIC void pah2_fun(double wavl,    /* in micron */
					 /*@in@*/ const sd_data *sd,
					 /*@in@*/ const grain_data *gd,
					 /*@out@*/ double *cs_abs,
					 /*@out@*/ double *cs_sct,
					 /*@out@*/ double *cosb,
					 /*@out@*/ int *error)
{
	DEBUG_ENTRY( "pah2_fun()" );

	// these choices give the best fit to the observed spectrum of the diffuse ISM
	// setting all these to 1 reproduces the laboratory measured band strengths
	const double E62 = 3.;
	const double E77 = 2.;
	const double E86 = 2.;

	// grain volume
	double vol = 4.*PI/3.*pow3(sd->cSize)*1.e-12;
	// number of carbon atoms in PAH molecule
	double xnc = vol*gd->rho/(ATOMIC_MASS_UNIT*dense.AtomicWeight[ipCARBON]);
	// this is the hydrogen over carbon ratio in the PAH molecule
	double xnhoc;
	// this is Eq. 4 of LD01
	if( xnc <= 25. )
		xnhoc = 0.5;
	else if( xnc <= 100. )
		xnhoc = 2.5/sqrt(xnc);
	else
		xnhoc = 0.25;

	double x = 1./wavl;

	double pah2_fun_v = 0.;
	if( x < 3.3 )
	{
		double M = ( xnc <= 40. ) ? 0.3*xnc : 0.4*xnc;
		// calculate cutoff wavelength in micron, this is Eq. A3 of LD01
		double cutoff;
		if( gd->charge == 0 )
			cutoff = 1./(3.804/sqrt(M) + 1.052);
		else
			cutoff = 1./(2.282/sqrt(M) + 0.889);
		double y = cutoff/wavl;
		// this is Eq. A2 of LD01
		double cutoff_fun = atan(1.e3*pow3(y-1.)/y)/PI + 0.5;
		// this is Eq. 11 of LD01
		pah2_fun_v = 34.58*exp10(  -18.-3.431/x )*cutoff_fun;
		for( int j=2; j < 14; ++j )
		{
			double strength = ( gd->charge == 0 ) ? pah2n_strength[j] : pah2c_strength[j];
			if( j == 2 )
				strength *= xnhoc;
			else if( j == 3 )
				strength *= E62;
			else if( j == 4 )
				strength *= E77;
			else if( j == 5 )
				strength *= E86*xnhoc;
			else if( j == 6 || j == 7 || j == 8 )
				strength *= xnhoc/3.;
			pah2_fun_v += Drude( wavl, pah2_wavl[j], pah2_width[j], strength );
		}
	}
	else if( x < 5.9 )
	{
		// this is Eq. 10 of LD01, strength is identical for charged and neutral PAHs
		pah2_fun_v = Drude( wavl, pah2_wavl[1], pah2_width[1], pah2n_strength[1] );
		pah2_fun_v += (1.8687 + 0.1905*x)*1.e-18;
	}
	else if( x < 7.7 )
	{
		// this is Eq. 9 of LD01, strength is identical for charged and neutral PAHs
		pah2_fun_v = Drude( wavl, pah2_wavl[1], pah2_width[1], pah2n_strength[1] );
		double y = x - 5.9;
		pah2_fun_v += (1.8687 + 0.1905*x + pow2(y)*(0.4175 + 0.04370*y))*1.e-18;
	}
	else if( x < 10. )
	{
		// this is Eq. 8 of LD01
		pah2_fun_v = (((-0.1057*x + 2.950)*x - 24.367)*x + 66.302)*1.e-18;
	}
	else if( x < 15. )
	{
		// this is Eq. 7 of LD01, strength is identical for charged and neutral PAHs
		pah2_fun_v = Drude( wavl, pah2_wavl[0], pah2_width[0], pah2n_strength[0] );
		pah2_fun_v += (-3. + 1.35*x)*1.e-18;
	}
	else if( x < 17.26 )
	{
		// this is Eq. 6 of LD01
		pah2_fun_v = (126.0 - 6.4943*x)*1.e-18;
	}
	else
		// this routine should never be called for wavelengths shorter than 1/17.25 micron
		// graphite opacities should be used in that case; this is handled in car2_fun()
		TotalInsanity();

	// normalize cross section per PAH molecule
	pah2_fun_v *= xnc;

	*cs_abs = pah2_fun_v;
	// the next two numbers are completely arbitrary
	*cs_sct = 0.1*pah2_fun_v;
	*cosb = 0.;
	*error = 0;

	return;
}

/* this routine calculates the absorption cross sections of carbonaceous grains, it is based on Eqs. 5-7 of:
 * >>refer	grain	physics	Draine, B.T., & Li, A., 2007 ApJ, 657, 810 */
inline void car3_fun(double wavl,    /* in micron */
					 /*@in@*/ const sd_data *sd,
					 /*@in@*/ const grain_data gdArr[], /* gdArr[2] */
					 /*@out@*/ double *cs_abs,
					 /*@out@*/ double *cs_sct,
					 /*@out@*/ double *cosb,
					 /*@out@*/ int *error)
{
	ld01_fun(pah3_fun,0.01,1./17.25,wavl,sd,gdArr,cs_abs,cs_sct,cosb,error);
}

// these values are taken from Table 1 of DL07
static const double pah3_wavl[30] = { 0.0722, 0.2175, 1.05, 1.26, 1.905, 3.3, 5.27, 5.7, 6.22, 6.69, 7.417, 7.598,
				      7.85, 8.33, 8.61, 10.68, 11.23, 11.33, 11.99, 12.62, 12.69, 13.48, 14.19,
				      15.9, 16.45, 17.04, 17.375, 17.87, 18.92, 15. };
static const double pah3_width[30] = { 0.195, 0.217, 0.055, 0.11, 0.09, 0.012, 0.034, 0.035, 0.03, 0.07, 0.126,
				       0.044, 0.053, 0.052, 0.039, 0.02, 0.012, 0.032, 0.045, 0.042, 0.013, 0.04,
				       0.025, 0.02, 0.014, 0.065, 0.012, 0.016, 0.1, 0.8 };
static const double pah3n_strength[30] = { 7.97e-13, 1.23e-13, 0., 0., 0., 3.94e-18, 2.5e-20, 4.e-20, 2.94e-19,
					   7.35e-20, 2.08e-19, 1.81e-19, 2.19e-19, 6.94e-20, 2.78e-19, 3.e-21,
					   1.89e-19, 5.2e-19, 2.42e-19, 3.5e-19, 1.3e-20, 8.e-20, 4.5e-21,
					   4.e-22, 5.e-21, 2.22e-20, 1.1e-21, 6.7e-22, 1.e-21, 5.e-19 };
static const double pah3c_strength[30] = { 7.97e-13, 1.23e-13, 2.e-16, 7.8e-17, -1.465e-18, 8.94e-19, 2.e-19,
					   3.2e-19, 2.35e-18, 5.9e-19, 1.81e-18, 1.63e-18, 1.97e-18, 4.8e-19,
					   1.94e-18, 3.e-21, 1.77e-19, 4.9e-19, 2.05e-19, 3.1e-19, 1.3e-20,
					   8.e-20, 4.5e-21, 4.e-22, 5.e-21, 2.22e-20, 1.1e-21, 6.7e-22, 1.7e-21,
					   5.e-19 };
// should we multiply the strength with H/C ?
static const bool pah3_hoc[30] = { false, false, false, false, false, true, false, false, false, false, false,
				   false, false, true, true, true, true, true, true, true, true, true, false,
				   false, false, false, false, false, false, false };

/* this routine calculates the absorption cross sections of PAH molecules, it is based on
 * >>refer	grain	physics	Draine, B.T., & Li, A., 2007 ApJ, 657, 810 */
STATIC void pah3_fun(double wavl,    /* in micron */
					 /*@in@*/ const sd_data *sd,
					 /*@in@*/ const grain_data *gd,
					 /*@out@*/ double *cs_abs,
					 /*@out@*/ double *cs_sct,
					 /*@out@*/ double *cosb,
					 /*@out@*/ int *error)
{
	DEBUG_ENTRY( "pah3_fun()" );

	// grain volume
	double vol = 4.*PI/3.*pow3(sd->cSize)*1.e-12;
	// number of carbon atoms in PAH molecule
	double xnc = vol*gd->rho/(ATOMIC_MASS_UNIT*dense.AtomicWeight[ipCARBON]);
	// this is the hydrogen over carbon ratio in the PAH molecule
	double xnhoc;
	// this is Eq. 4 of DL07
	if( xnc <= 25. )
		xnhoc = 0.5;
	else if( xnc <= 100. )
		xnhoc = 2.5/sqrt(xnc);
	else
		xnhoc = 0.25;

	double x = 1./wavl;

	// this is Eq. 2 of DL07
	double pah3_fun_v = ( gd->charge == 0 ) ? 0. : 3.5*exp10(-19.-1.45/x)*exp(-0.1*pow2(x));

	if( x < 3.3 )
	{
		double M = ( xnc <= 40. ) ? 0.3*xnc : 0.4*xnc;
		// calculate cutoff wavelength in micron, this is Eq. A3 of LD01
		double cutoff;
		if( gd->charge == 0 )
			cutoff = 1./(3.804/sqrt(M) + 1.052);
		else
			cutoff = 1./(2.282/sqrt(M) + 0.889);
		double y = cutoff/wavl;
		// this is Eq. A2 of LD01
		double cutoff_fun = atan(1.e3*pow3(y-1.)/y)/PI + 0.5;
		// this is Eq. 11 of LD01
		pah3_fun_v += 34.58*exp10(  -18.-3.431/x )*cutoff_fun;
		for( int j=2; j < 30; ++j )
		{
			double strength = ( gd->charge == 0 ) ? pah3n_strength[j] : pah3c_strength[j];
			if( pah3_hoc[j] )
				strength *= xnhoc;
			pah3_fun_v += Drude( wavl, pah3_wavl[j], pah3_width[j], strength );
		}
	}
	else if( x < 5.9 )
	{
		// this is Eq. 10 of LD01, strength is identical for charged and neutral PAHs
		pah3_fun_v += Drude( wavl, pah3_wavl[1], pah3_width[1], pah3n_strength[1] );
		pah3_fun_v += (1.8687 + 0.1905*x)*1.e-18;
	}
	else if( x < 7.7 )
	{
		// this is Eq. 9 of LD01, strength is identical for charged and neutral PAHs
		pah3_fun_v += Drude( wavl, pah3_wavl[1], pah3_width[1], pah3n_strength[1] );
		double y = x - 5.9;
		pah3_fun_v += (1.8687 + 0.1905*x + pow2(y)*(0.4175 + 0.04370*y))*1.e-18;
	}
	else if( x < 10. )
	{
		// this is Eq. 8 of LD01
		pah3_fun_v += (((-0.1057*x + 2.950)*x - 24.367)*x + 66.302)*1.e-18;
	}
	else if( x < 15. )
	{
		// this is Eq. 7 of LD01, strength is identical for charged and neutral PAHs
		pah3_fun_v += Drude( wavl, pah3_wavl[0], pah3_width[0], pah3n_strength[0] );
		pah3_fun_v += (-3. + 1.35*x)*1.e-18;
	}
	else if( x < 17.26 )
	{
		// this is Eq. 6 of LD01
		pah3_fun_v += (126.0 - 6.4943*x)*1.e-18;
	}
	else
		// this routine should never be called for wavelengths shorter than 1/17.25 micron
		// graphite opacities should be used in that case; this is handled in car3_fun()
		TotalInsanity();

	// normalize cross section per PAH molecule
	pah3_fun_v *= xnc;

	*cs_abs = pah3_fun_v;
	// the next two numbers are completely arbitrary
	*cs_sct = 0.1*pah3_fun_v;
	*cosb = 0.;
	*error = 0;

	return;
}

// Drude: helper function to calculate Drude profile, return value is cross section in cm^2/C
inline double Drude(double lambda,  // wavelength (in micron)
					double lambdac, // central wavelength of feature (in micron)
					double gamma,   // width of the feature (dimensionless)
					double sigma)   // strength of the feature (in cm/C)
{
	double x = lambda/lambdac;
	// this is Eq. 12 of LD01
	return 2.e-4/PI*gamma*lambdac*sigma/(pow2(x - 1./x) + pow2(gamma));
}

STATIC void tbl_fun(double wavl,    /* in micron */
					/*@in@*/ const sd_data *sd, /* NOT USED -- MUST BE KEPT FOR COMPATIBILITY */
					/*@in@*/ const grain_data *gd,
					/*@out@*/ double *cs_abs,
					/*@out@*/ double *cs_sct,
					/*@out@*/ double *cosb,
					/*@out@*/ int *error)
{
	bool lgOutOfBounds;
	long int ind;
	double anu = WAVNRYD/wavl*1.e4;

	DEBUG_ENTRY( "tbl_fun()" );

	/* >>chng 02 nov 17, add this test to prevent warning that this var not used */
	if( sd == NULL )
		TotalInsanity();

	/** \todo	2	include code for interpolating inv_att_len somewhere!! */

	find_arr(anu,gd->opcAnu,gd->nOpcData,&ind,&lgOutOfBounds);
	if( !lgOutOfBounds ) 
	{
		double a1g;
		double frac = log(anu/gd->opcAnu[ind])/log(gd->opcAnu[ind+1]/gd->opcAnu[ind]);

		*cs_abs = exp((1.-frac)*log(gd->opcData[0][ind])+frac*log(gd->opcData[0][ind+1]));
		ASSERT( *cs_abs > 0. );
		if( gd->nOpcCols > 1 )
			*cs_sct = exp((1.-frac)*log(gd->opcData[1][ind])+frac*log(gd->opcData[1][ind+1]));
		else
			*cs_sct = 0.1*(*cs_abs);
		ASSERT( *cs_sct > 0. );
		if( gd->nOpcCols > 2 )
			a1g = exp((1.-frac)*log(gd->opcData[2][ind])+frac*log(gd->opcData[2][ind+1]));
		else
			a1g = 1.;
		ASSERT( a1g > 0. );
		*cosb = 1. - a1g;
		*error = 0;
	}
	else
	{
		*cs_abs = -1.;
		*cs_sct = -1.;
		*cosb = -2.;
		*error = 3;
	}
	return;
}

STATIC double size_distr(double size,
						 /*@in@*/ const sd_data *sd)
{
	bool lgOutOfBounds;
	long ind;
	double frac,
	  res,
	  x;

	DEBUG_ENTRY( "size_distr()" );

	if( size >= sd->lim[ipBLo] && size <= sd->lim[ipBHi] )
		switch( sd->sdCase ) 
		{
		case SD_SINGLE_SIZE:
		case SD_NR_CARBON:
			res = 1.; /* should really not be used in this case */
			break;
		case SD_POWERLAW:
			/* simple powerlaw */
		case SD_EXP_CUTOFF1:
		case SD_EXP_CUTOFF2:
		case SD_EXP_CUTOFF3:
			/* powerlaw with exponential cutoff, inspired by Greenberg (1978)
			 * Cosmic Dust, ed. J.A.M. McDonnell, Wiley, p. 187 */
			res = pow(size,sd->a[ipExp]);
			if( sd->a[ipBeta] < 0. )
				res /= (1. - sd->a[ipBeta]*size);
			else if( sd->a[ipBeta] > 0. )
				res *= (1. + sd->a[ipBeta]*size);
			if( size < sd->a[ipBLo] && sd->a[ipSLo] > 0. )
				res *= exp(-powi((sd->a[ipBLo]-size)/sd->a[ipSLo],nint(sd->a[ipAlpha])));
			if( size > sd->a[ipBHi] && sd->a[ipSHi] > 0. )
				res *= exp(-powi((size-sd->a[ipBHi])/sd->a[ipSHi],nint(sd->a[ipAlpha])));
			break;
		case SD_LOG_NORMAL:
			x = log(size/sd->a[ipGCen])/sd->a[ipGSig];
			res = exp(-0.5*pow2(x))/size;
			break;
		case SD_LIN_NORMAL:
			x = (size-sd->a[ipGCen])/sd->a[ipGSig];
			res = exp(-0.5*pow2(x))/size;
			break;
		case SD_TABLE:
			find_arr(log(size),sd->ln_a,sd->npts,&ind,&lgOutOfBounds);
			if( lgOutOfBounds )
			{
				fprintf( ioQQQ, " size distribution table has insufficient range\n" );
				fprintf( ioQQQ, " requested size: %.5f table range %.5f - %.5f\n",
					 size, exp(sd->ln_a[0]), exp(sd->ln_a[sd->npts-1]) );
				cdEXIT(EXIT_FAILURE);
			}				
			frac = (log(size)-sd->ln_a[ind])/(sd->ln_a[ind+1]-sd->ln_a[ind]);
			ASSERT( frac > 0.-10.*DBL_EPSILON && frac < 1.+10.*DBL_EPSILON );
			res = (1.-frac)*sd->ln_a4dNda[ind] + frac*sd->ln_a4dNda[ind+1];
			/* convert from a^4*dN/da to dN/da */
			res = exp(res)/POW4(size);
			break;
		case SD_ILLEGAL:
		default:
			TotalInsanity();
		}
	else
		res = 0.;
	return res;
}

/* search for upper/lower limit lim of size distribution such that
 * lim^4 * dN/da(lim) < rel_cutoff * ref^4 * dN/da(ref)
 * the initial estimate of lim = ref + step
 * step may be both positive (upper limit) or negative (lower limit) */ 
STATIC double search_limit(double ref,
						   double step,
						   double rel_cutoff,
						   // do NOT use sd_data* since that would trash sd.lim[ipBLo]
						   sd_data sd)
{
	long i;
	double f1,
	  f2,
	  fmid,
	  renorm,
	  x1,
	  x2 = DBL_MAX,
	  xmid = DBL_MAX;

	/* TOLER is the relative accuracy with which lim is determined */
	const double TOLER = 1.e-6;

	DEBUG_ENTRY( "search_limit()" );

	/* sanity check */
	ASSERT( rel_cutoff > 0. && rel_cutoff < 1. );

	if( step == 0. )
	{
		return ref;
	}

	/* these need to be set in order for size_distr to work...
	 * NB - since this is a local copy of sd, it will not
	 * upset anything in the calling routine */
	sd.lim[ipBLo] = 0.;
	sd.lim[ipBHi] = DBL_MAX;

	x1 = ref;
	/* previous assert guarantees that f1 > 0. */
	f1 = -log(rel_cutoff);
	renorm = f1 - log(POW4(x1)*size_distr(x1,&sd));

	/* bracket solution */
	f2 = 1.;
	for( i=0; i < 20 && f2 > 0.; ++i )
	{
		x2 = max(ref + step,SMALLEST_GRAIN);
		f2 = log(POW4(x2)*size_distr(x2,&sd)) + renorm;
		if( f2 >= 0. )
		{
			x1 = x2;
			f1 = f2;
		}
		step *= 2.;
	}
	if( f2 > 0. )
	{
		fprintf( ioQQQ, " Could not bracket solution\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* do bisection search */
	while( 2.*fabs(x1-x2)/(x1+x2) > TOLER )
	{
		xmid = (x1+x2)/2.;
		fmid = log(POW4(xmid)*size_distr(xmid,&sd)) + renorm;

		if( fmid == 0. )
			break;

		if( f1*fmid > 0. )
		{
			x1 = xmid;
			f1 = fmid;
		}
		else
		{
			x2 = xmid;
//			f2 = fmid;
		}
	}
	return (x1+x2)/2.;
}

/* calculate the inverse attenuation length for given refractive index data */
STATIC void mie_calc_ial(/*@in@*/ const grain_data *gd,
						 long int n,
						 /*@out@*/ vector<double>& invlen,  /* invlen[n] */
						 /*@in@*/ const string& chString,
						 /*@in@*/ bool *lgWarning)
{
	bool lgErrorOccurred=true,
	  lgOutOfBounds;
	long int i,
	  ind,
	  j;
	double frac,
	  InvDep,
	  nim,
	  wavlen;

	DEBUG_ENTRY( "mie_calc_ial()" );

	/* sanity check */
	ASSERT( gd->rfiType == RFI_TABLE );

	vector<int> ErrorIndex( rfield.nflux_with_check );

	for( i=0; i < n; i++ ) 
	{
		wavlen = WAVNRYD/rfield.anu(i)*1.e4;

		ErrorIndex[i] = 0;
		lgErrorOccurred = false;
		invlen[i] = 0.;

		for( j=0; j < gd->nAxes; j++ ) 
		{
			/* first interpolate optical constants */
			find_arr(wavlen,gd->wavlen[j],gd->ndata[j],&ind,&lgOutOfBounds);
			if( lgOutOfBounds ) 
			{
				ErrorIndex[i] = 3;
				lgErrorOccurred = true;
				invlen[i] = 0.;
				break;
			}
			frac = (wavlen-gd->wavlen[j][ind])/(gd->wavlen[j][ind+1]-gd->wavlen[j][ind]);
			nim = (1.-frac)*gd->n[j][ind].imag() + frac*gd->n[j][ind+1].imag();
			/* this is the inverse of the photon attenuation depth,
			 * >>refer	grain	physics	Weingartner & Draine, 2000, ApJ, ... */
			InvDep = PI4*nim/wavlen*1.e4;
			ASSERT( InvDep > 0. );

			invlen[i] += InvDep*gd->wt[j];
		}
	}

	if( lgErrorOccurred ) 
	{
		mie_repair(chString,n,3,3,rfield.anuptr(),&invlen[0],ErrorIndex,false,lgWarning);
	}

	return;
}

/* this is the number of x-values we use for extrapolating functions */
const int NPTS_DERIV = 8;

/* extrapolate/interpolate mie data to fill in the blanks */
STATIC void mie_repair(/*@in@*/ const string& chString,
					   long int n,
					   int val,
					   int del,
					   /*@in@*/ const double anu[],     /* anu[n] */
					   double data[],                    /* data[n] */
					   /*@in@*/ vector<int>& ErrorIndex, /* ErrorIndex[n] */
					   bool lgRound,
					   /*@in@*/ bool *lgWarning)
{
	bool lgExtrapolate,
	  lgVerbose;
	long int i1,
	  i2,
	  ind1,
	  ind2,
	  j;
	double dx,
	  sgn,
	  slp1,
	  xlg1,
	  xlg2,
	  y1lg1,
	  y1lg2;

	/* interpolating over more that this number of points results in a warning */
	const long BIG_INTERPOLATION = 10;

	DEBUG_ENTRY( "mie_repair()" );

	lgVerbose = ( chString.length() != 0 );

	for( ind1=0; ind1 < n; ) 
	{
		if( ErrorIndex[ind1] == val ) 
		{
			/* search for region with identical error index */
			ind2 = ind1;
			bool lgValid = true;
			do
			{
				while( ind2 < n && ErrorIndex[ind2] == val )
					ind2++;
				// check if the following points can be used to determine slope
				// this check is only needed for the low energy extrapolation,
				// hence the test for ind1 == 0...
				for( j=ind2; ind1 == 0 && j < min(ind2+NPTS_DERIV,n); j++ )
				{
					if( ErrorIndex[j] == val )
					{
						lgValid= false;
						break;
					}
				}
				for( j=ind2; !lgValid && j < min(ind2+NPTS_DERIV,n); j++ )
				{
					if( ErrorIndex[j] == val )
					{
						ind2 = j+1;
						break;
					}
					else
					{
						ErrorIndex[j] = val;
					}
				}
			}
			while( !lgValid && ind2 < n );

			if( lgVerbose )
				fprintf( ioQQQ, "    %s", chString.c_str() );

			if( ind1 == 0 ) 
			{
				/* low energy extrapolation */
				i1 = ind2;
				i2 = ind2+NPTS_DERIV-1;
				lgExtrapolate = true;
				sgn = +1.;
				if( lgVerbose ) 
				{
					fprintf( ioQQQ, " extrapolated below %.4e Ryd\n",anu[i1] );
				}
			}
			else if( ind2 == n ) 
			{
				/* high energy extrapolation */
				i1 = ind1-NPTS_DERIV;
				i2 = ind1-1;
				lgExtrapolate = true;
				sgn = -1.;
				if( lgVerbose ) 
				{
					fprintf( ioQQQ, " extrapolated above %.4e Ryd\n",anu[i2] );
				}
			}
			else 
			{
				/* interpolation */
				i1 = ind1-1;
				i2 = ind2;
				lgExtrapolate = false;
				sgn = 0.;
				if( lgVerbose ) 
				{
					fprintf( ioQQQ, " interpolated between %.4e and %.4e Ryd\n",
						 anu[i1],anu[i2] );
				}
				if( i2-i1-1 > BIG_INTERPOLATION )
				{
					if( lgVerbose )
					{
						fprintf( ioQQQ, " ***Warning: extensive interpolation used\n" );
					}
					*lgWarning = true;
				}
			}

			if( i1 < 0 || i2 >= n ) 
			{
				fprintf( ioQQQ, " Insufficient data for extrapolation\n" );
				cdEXIT(EXIT_FAILURE);
			}

			xlg1 = log(anu[i1]);
			y1lg1 = log(data[i1]);
			/* >>chng 01 jul 30, replace simple-minded extrapolation with more robust routine, PvH */
			if( lgExtrapolate )
				slp1 = mie_find_slope(anu,data,ErrorIndex,i1,i2,val,lgVerbose,lgWarning);
			else
			{
				xlg2 = log(anu[i2]);
				y1lg2 = log(data[i2]);
				slp1 = (y1lg2-y1lg1)/(xlg2-xlg1);
			}
			if( lgRound && lgExtrapolate && sgn > 0. ) 
			{
				/* in low energy extrapolation, 1-g is very close to 1 and almost constant
				 * hence slp1 is very close to 0 and can even be slightly negative
				 * to prevent 1-g becoming greater than 1, the following is necessary */
				slp1 = max(slp1,0.);
			}
			/* >>chng 02 oct 22, changed from sgn*slp1 <= 0. to accomodate grey grains */
			else if( lgExtrapolate && sgn*slp1 < 0. ) 
			{
				fprintf( ioQQQ, " Unphysical value for slope in extrapolation %.6e\n", slp1 );
				fprintf( ioQQQ, " The most likely cause is that your refractive index or "
					 "opacity data do not extend to low or high enough frequencies. "
					 "See Hazy 1 for more details.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			for( j=ind1; j < ind2; j++ ) 
			{
				dx = log(anu[j]) - xlg1;
				data[j] = exp(y1lg1 + dx*slp1);
				ErrorIndex[j] -= del;
			}

			ind1 = ind2;
		}
		else 
		{
			ind1++;
		}
	}
	/* sanity check */
	for( j=0; j < n; j++ )
	{
		if( ErrorIndex[j] > val-del ) 
			TotalInsanity();
	}
	return;
}

STATIC double mie_find_slope(/*@in@*/ const double anu[],
							 /*@in@*/ const double data[],
							 /*@in@*/ const vector<int>& ErrorIndex,
							 long i1,
							 long i2,
							 int val,
							 bool lgVerbose,
							 /*@in@*/ bool *lgWarning)
{
	/* threshold for standard deviation in the logarithmic derivative to generate warning,
	 * this corresponds to an uncertainty of a factor 10 for a typical extrapolation */
	const double LARGE_STDEV = 0.2;

	DEBUG_ENTRY( "mie_find_slope()" );

	/* sanity check */
	for( long i=i1; i <= i2; i++ )
	{
		if( ! ( ErrorIndex[i] < val && anu[i] > 0. && data[i] > 0. ) )
		{
			fprintf( ioQQQ, "invalid parameter for mie_find_slope\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* calculate the logarithmic derivative for every possible combination of data points */
	vector<double> slp1;
	for( long i=i1; i < i2; i++ )
		for( long j=i+1; j <= i2; j++ )
			slp1.push_back( log(data[j]/data[i])/log(anu[j]/anu[i]) );
	/* sort the values; we want the median */
	size_t n = slp1.size()/2;
	partial_sort( slp1.begin(), slp1.begin()+n+1, slp1.end() );
	/* now calculate the median value */
	double slope = ( slp1.size()%2 == 1 ) ? slp1[n] : (slp1[n-1]+slp1[n])/2.;

	/* and finally calculate the standard deviation of all slopes */
	double s1 = 0.;
	double s2 = 0.;
	for( size_t i=0; i < slp1.size(); i++ )
	{
		s1 += slp1[i];
		s2 += pow2(slp1[i]);
	}
	/* >>chng 06 jul 12, protect against roundoff error, PvH */
	double stdev = sqrt(max(s2/(double)slp1.size() - pow2(s1/(double)slp1.size()),0.));

	/* print warning if standard deviation is large */
	if( stdev > LARGE_STDEV )
	{
		if( lgVerbose )
			fprintf( ioQQQ, " ***Warning: slope for extrapolation may be unreliable\n" );
		*lgWarning = true;
	}
	return slope;
}

/* read optical constants using either a .rfi or .mix file */
STATIC void mie_read_ocn(/*@in@*/  const string& chFile,
						 /*@out@*/ grain_data *gd)
{
	if( chFile.find( ".rfi" ) != string::npos )
	{
		mie_read_rfi( chFile, gd );
	}
	else if( chFile.find( ".mix" ) != string::npos )
	{
		mie_read_mix( chFile, gd );
	}
	else
	{
		fprintf( ioQQQ, " Refractive index file name %s has wrong extention\n", chFile.c_str() );
		fprintf( ioQQQ, " It should have extention .rfi or .mix.\n" );
		cdEXIT(EXIT_FAILURE);
	}
}

/* read in the file with optical constants and other relevant information */
STATIC void mie_read_rfi(/*@in@*/  const string& chFile,
						 /*@out@*/ grain_data *gd)
{
	bool lgLogData = false;
	long int dl, help, i, nelem, j, nridf, sgn = 0;
	double eps1, eps2, LargestLog, molw, nAtoms, nr, ni, tmp1, tmp2, total = 0.;
	string chLine, chWord;

	DEBUG_ENTRY( "mie_read_rfi()" );

	FILE *io2 = open_data( chFile, "r" );

	dl = 0; /* line counter for input file */

	/* first read magic number */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&gd->magic,true,dl);
	if( gd->magic != MAGIC_RFI ) 
	{
		fprintf( ioQQQ, " Refractive index file %s has obsolete magic number\n",chFile.c_str() );
		fprintf( ioQQQ, " I found %ld, but expected %ld on line #%ld\n",gd->magic,MAGIC_RFI,dl );
		fprintf( ioQQQ, " Please replace this file with an up to date version\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* get chemical formula of the grain, e.g., Mg0.4Fe0.6SiO3 */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_word(chLine,chWord,false);
	mie_read_form(chWord,gd->elmAbun,&nAtoms,&molw);

	/* molecular weight, in atomic units */
	gd->mol_weight = molw;
	gd->atom_weight = gd->mol_weight/nAtoms;

	/* determine abundance of grain molecule assuming max depletion */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_double(chFile,chLine,&gd->abun,true,dl);

	/* default depletion of grain molecule */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_double(chFile,chLine,&gd->depl,true,dl);
	if( gd->depl > 1. ) 
	{
		fprintf( ioQQQ, " Illegal value for default depletion in %s\n",chFile.c_str() );
		fprintf( ioQQQ, " Line #%ld, depl=%14.6e\n",dl,gd->depl);
		cdEXIT(EXIT_FAILURE);
	}

	for( nelem=0; nelem < LIMELM; nelem++ )
		gd->elmAbun[nelem] *= gd->abun*gd->depl;

	/* material density, to get cross section per unit mass: rho in cgs */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_double(chFile,chLine,&gd->rho,false,dl);

	/* material type, determines enthalpy function: 1 -- carbonaceous, 2 -- silicate */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&help,true,dl);
	gd->matType = (mat_type)help;
	if( gd->matType >= MAT_TOP ) 
	{
		fprintf( ioQQQ, " Illegal value for material type in %s\n",chFile.c_str() );
		fprintf( ioQQQ, " Line #%ld, type=%d\n",dl,gd->matType);
		cdEXIT(EXIT_FAILURE);
	}

	/* work function, in Rydberg */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_double(chFile,chLine,&gd->work,true,dl);

	/* bandgap, in Rydberg */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_double(chFile,chLine,&gd->bandgap,false,dl);
	if( gd->bandgap >= gd->work ) 
	{
		fprintf( ioQQQ, " Illegal value for bandgap in %s\n",chFile.c_str() );
		fprintf( ioQQQ, " Line #%ld, bandgap=%.4e, work function=%.4e\n",dl,gd->bandgap,gd->work);
		fprintf( ioQQQ, " Bandgap should always be less than work function\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* efficiency of thermionic emission, between 0 and 1 */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_double(chFile,chLine,&gd->therm_eff,true,dl);
	if( gd->therm_eff > 1.f ) 
	{
		fprintf( ioQQQ, " Illegal value for thermionic efficiency in %s\n",chFile.c_str() );
		fprintf( ioQQQ, " Line #%ld, value=%.4e\n",dl,gd->therm_eff);
		fprintf( ioQQQ, " Allowed values are 0. < efficiency <= 1.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* sublimation temperature in K */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_double(chFile,chLine,&gd->subl_temp,true,dl);

	/* >>chng 02 sep 18, add keyword for special files (grey grains, PAH's, etc.), PvH */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_word(chLine,chWord,true);

	if( chWord.find( "RFI_" ) != string::npos )
		gd->rfiType = RFI_TABLE;
	else if( chWord.find( "OPC_" ) != string::npos )
		gd->rfiType = OPC_TABLE;
	else if( chWord.find( "GREY" ) != string::npos )
		gd->rfiType = OPC_GREY;
	else if( chWord.find( "PAH1" ) != string::npos )
		gd->rfiType = OPC_PAH1;
	else if( chWord.find( "PH2N" ) != string::npos )
		gd->rfiType = OPC_PAH2N;
	else if( chWord.find( "PH2C" ) != string::npos )
		gd->rfiType = OPC_PAH2C;
	else if( chWord.find( "PH3N" ) != string::npos )
		gd->rfiType = OPC_PAH3N;
	else if( chWord.find( "PH3C" ) != string::npos )
		gd->rfiType = OPC_PAH3C;
	else
	{
		fprintf( ioQQQ, " Illegal keyword in %s\n",chFile.c_str() );
		fprintf( ioQQQ, " Line #%ld, value=%s\n",dl,chWord.c_str());
		fprintf( ioQQQ, " Allowed values are: RFI_TBL, OPC_TBL, GREY, PAH1, PH2N, PH2C, PH3N, PH3C\n");
		cdEXIT(EXIT_FAILURE);
	}

	switch( gd->rfiType )
	{
	case RFI_TABLE:
		/* nridf is for choosing ref index or diel function input
		 * case 2 allows greater accuracy reading in, when nr is close to 1. */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_long(chFile,chLine,&nridf,true,dl);
		if( nridf > 3 ) 
		{
			fprintf( ioQQQ, " Illegal data code in %s\n",chFile.c_str() );
			fprintf( ioQQQ, " Line #%ld, data code=%ld\n",dl,nridf);
			cdEXIT(EXIT_FAILURE);
		}

		/* no. of principal axes, always 1 for amorphous grains,
		 * maybe larger for crystalline grains */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_long(chFile,chLine,&gd->nAxes,true,dl);
		if( gd->nAxes > NAX ) 
		{
			fprintf( ioQQQ, " Illegal no. of axes in %s\n",chFile.c_str() );
			fprintf( ioQQQ, " Line #%ld, number=%ld\n",dl,gd->nAxes);
			cdEXIT(EXIT_FAILURE);
		}

		/* now get relative weights of axes */
		mie_next_data(chFile,io2,chLine,&dl);
		switch( gd->nAxes ) 
		{
		case 1:
			mie_read_double(chFile,chLine,&gd->wt[0],true,dl);
			total = gd->wt[0];
			break;
		case 2:
			if( sscanf( chLine.c_str(), "%lf %lf", &gd->wt[0], &gd->wt[1] ) != 2 ) 
			{
				fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
				fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
				cdEXIT(EXIT_FAILURE);
			}
			if( gd->wt[0] <= 0. || gd->wt[1] <= 0. ) 
			{
				fprintf( ioQQQ, " Illegal data in %s\n",chFile.c_str());
				fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
				cdEXIT(EXIT_FAILURE);
			}
			total = gd->wt[0] + gd->wt[1];
			break;
		case 3:
			if( sscanf( chLine.c_str(), "%lf %lf %lf", &gd->wt[0], &gd->wt[1], &gd->wt[2] ) != 3 ) 
			{
				fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
				fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
				cdEXIT(EXIT_FAILURE);
			}
			if( gd->wt[0] <= 0. || gd->wt[1] <= 0. || gd->wt[2] <= 0. ) 
			{
				fprintf( ioQQQ, " Illegal data in %s\n",chFile.c_str());
				fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
				cdEXIT(EXIT_FAILURE);
			}
			total = gd->wt[0] + gd->wt[1] + gd->wt[2];
			break;
		default:
			TotalInsanity();
		}
		for( j=0; j < gd->nAxes; j++ ) 
		{
			gd->wt[j] /= total;

			/* read in optical constants for each principal axis. */
			mie_next_data(chFile,io2,chLine,&dl);
			mie_read_long(chFile,chLine,&gd->ndata[j],false,dl);
			if( gd->ndata[j] < 2 ) 
			{
				fprintf( ioQQQ, " Illegal number of data points in %s\n",chFile.c_str() );
				fprintf( ioQQQ, " Line #%ld, number=%ld\n",dl,gd->ndata[j]);
				cdEXIT(EXIT_FAILURE);
			}

			/* allocate space for wavelength and optical constants arrays */
			gd->wavlen[j].resize( gd->ndata[j] );
			gd->n[j].resize( gd->ndata[j] );
			gd->nr1[j].resize( gd->ndata[j] );

			for( i=0; i < gd->ndata[j]; i++ ) 
			{
				/* read in the wavelength in microns
				 * and the complex refractive index or dielectric function of material */
				mie_next_data(chFile,io2,chLine,&dl);
				if( sscanf( chLine.c_str(), "%lf %lf %lf", &gd->wavlen[j][i], &nr, &ni ) != 3 ) 
				{
					fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
					fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
					cdEXIT(EXIT_FAILURE);
				}
				if( gd->wavlen[j][i] <= 0. ) 
				{
					fprintf( ioQQQ, " Illegal value for wavelength in %s\n",chFile.c_str());
					fprintf( ioQQQ, " Line #%ld, wavl=%14.6e\n",dl,gd->wavlen[j][i]);
					cdEXIT(EXIT_FAILURE);
				}
				/* the data in the input file should be sorted on wavelength, either
				 * strictly monotonically increasing or decreasing, check this here... */
				if( i == 1 )
				{
					sgn = sign3(gd->wavlen[j][1]-gd->wavlen[j][0]);
					if( sgn == 0 ) 
					{
						fprintf( ioQQQ, " Illegal value for wavelength in %s\n",chFile.c_str());
						fprintf( ioQQQ, " Line #%ld, wavl=%14.6e\n",dl,gd->wavlen[j][i]);
						cdEXIT(EXIT_FAILURE);
					}
				}
				else if( i > 1 ) 
				{
					if( sign3(gd->wavlen[j][i]-gd->wavlen[j][i-1]) != sgn ) 
					{
						fprintf( ioQQQ, " Illegal value for wavelength in %s\n",chFile.c_str());
						fprintf( ioQQQ, " Line #%ld, wavl=%14.6e\n",dl,gd->wavlen[j][i]);
						cdEXIT(EXIT_FAILURE);
					}
				}
				/* this version reads in real and imaginary parts of the refractive
				 * index, with imaginary part positive (nridf = 3) or nr-1 (nridf = 2) or
				 * real and imaginary parts of the dielectric function (nridf = 1) */
				switch( nridf ) 
				{
				case 1:
					eps1 = nr;
					eps2 = ni;
					dftori(&nr,&ni,eps1,eps2);
					gd->nr1[j][i] = nr - 1.;
					break;
				case 2:
					gd->nr1[j][i] = nr;
					nr += 1.;
					break;
				case 3:
					gd->nr1[j][i] = nr - 1.;
					break;
				default:
					TotalInsanity();
				}
				gd->n[j][i] = complex<double>(nr,ni);

				/* sanity checks */
				if( nr <= 0. || ni < 0. ) 
				{
					fprintf( ioQQQ, " Illegal value for refractive index in %s\n",chFile.c_str());
					fprintf( ioQQQ, " Line #%ld, (nr,ni)=(%14.6e,%14.6e)\n",dl,nr,ni);
					cdEXIT(EXIT_FAILURE);
				}
				ASSERT( fabs(nr-1.-gd->nr1[j][i]) < 10.*nr*DBL_EPSILON );
			}
		}
		break;
	case OPC_TABLE:
		/* no. of data columns in OPC_TABLE file:
		 * 1: absorption cross sections only
		 * 2: absorption + scattering cross sections
		 * 3: absorption + pure scattering cross sections + asymmetry factor
		 * 4: absorption + pure scattering cross sections + asymmetry factor +
		 *      inverse attenuation length */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_long(chFile,chLine,&gd->nOpcCols,true,dl);
		if( gd->nOpcCols > NDAT ) 
		{
			fprintf( ioQQQ, " Illegal no. of data columns in %s\n",chFile.c_str() );
			fprintf( ioQQQ, " Line #%ld, number=%ld\n",dl,gd->nOpcCols);
			cdEXIT(EXIT_FAILURE);
		}

		/* keyword indicating whether the table contains linear or logarithmic data */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_word(chLine,chWord,true);

		if( chWord.find( "LINE" ) != string::npos )
			lgLogData = false;
		else if( chWord.find( "LOG" ) != string::npos )
			lgLogData = true;
		else
		{
			fprintf( ioQQQ, " Keyword not recognized in %s\n",chFile.c_str() );
			fprintf( ioQQQ, " Line #%ld, keyword=%s\n",dl,chWord.c_str());
			cdEXIT(EXIT_FAILURE);
		}


		/* read in number of data points supplied. */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_long(chFile,chLine,&gd->nOpcData,false,dl);
		if( gd->nOpcData < 2 ) 
		{
			fprintf( ioQQQ, " Illegal number of data points in %s\n",chFile.c_str() );
			fprintf( ioQQQ, " Line #%ld, number=%ld\n",dl,gd->nOpcData);
			cdEXIT(EXIT_FAILURE);
		}

		/* allocate space for frequency and data arrays */
		gd->opcAnu.resize(gd->nOpcData);
		for( j=0; j < gd->nOpcCols; j++ ) 
			gd->opcData[j].resize(gd->nOpcData);

		tmp1 = -log10(1.01*DBL_MIN);
		tmp2 = log10(0.99*DBL_MAX);
		LargestLog = min(tmp1,tmp2);

		/* now read the data... each line should contain:
		 *
		 * if gd->nOpcCols == 1, anu, abs_cs
		 * if gd->nOpcCols == 2, anu, abs_cs, sct_cs
		 * if gd->nOpcCols == 3, anu, abs_cs, sct_cs, inv_att_len
		 * 
		 * the frequencies in the table should be either monotonically increasing or decreasing.
		 * frequencies should be in Ryd, cross sections in cm^2/H, and the inverse attenuation length
		 * in cm^-1. If lgLogData is true, each number should be the log10 of the data value */
		for( i=0; i < gd->nOpcData; i++ ) 
		{
			mie_next_data(chFile,io2,chLine,&dl);
			switch( gd->nOpcCols )
			{
			case 1:
				if( sscanf( chLine.c_str(), "%lf %lf", &gd->opcAnu[i], &gd->opcData[0][i] ) != 2 ) 
				{
					fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
					fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
					cdEXIT(EXIT_FAILURE);
				}
				break;
			case 2:
				if( sscanf( chLine.c_str(), "%lf %lf %lf", &gd->opcAnu[i], &gd->opcData[0][i],
					    &gd->opcData[1][i] ) != 3 ) 
				{
					fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
					fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
					cdEXIT(EXIT_FAILURE);
				}
				break;
			case 3:
				if( sscanf( chLine.c_str(), "%lf %lf %lf %lf", &gd->opcAnu[i], &gd->opcData[0][i],
					    &gd->opcData[1][i], &gd->opcData[2][i] ) != 4 ) 
				{
					fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
					fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
					cdEXIT(EXIT_FAILURE);
				}
				break;
			case 4:
				if( sscanf( chLine.c_str(), "%lf %lf %lf %lf %lf", &gd->opcAnu[i], &gd->opcData[0][i],
					    &gd->opcData[1][i], &gd->opcData[2][i], &gd->opcData[3][i] ) != 5 ) 
				{
					fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
					fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
					cdEXIT(EXIT_FAILURE);
				}
				break;
			default:
				TotalInsanity();
			}
			if( lgLogData )
			{
				/* this test will guarantee there will be neither under- nor overflows */
				if( fabs(gd->opcAnu[i]) > LargestLog )
				{
					fprintf( ioQQQ, " Illegal value for frequency in %s\n",chFile.c_str() );
					fprintf( ioQQQ, " Line #%ld, freq=%14.6e\n",dl,gd->opcAnu[i] );
					cdEXIT(EXIT_FAILURE);
				}
				gd->opcAnu[i] = exp10(gd->opcAnu[i]);
				for( j=0; j < gd->nOpcCols; j++ )
				{
					if( fabs(gd->opcData[j][i]) > LargestLog )
					{
						fprintf( ioQQQ, " Illegal data value in %s\n",chFile.c_str() );
						fprintf( ioQQQ, " Line #%ld, value=%14.6e\n",dl,gd->opcData[j][i] );
						cdEXIT(EXIT_FAILURE);
					}
					gd->opcData[j][i] = exp10(gd->opcData[j][i]);
				}
			}
			if( gd->opcAnu[i] <= 0. )
			{
				fprintf( ioQQQ, " Illegal value for frequency in %s\n",chFile.c_str() );
				fprintf( ioQQQ, " Line #%ld, freq=%14.6e\n",dl,gd->opcAnu[i] );
				cdEXIT(EXIT_FAILURE);
			}
			for( j=0; j < gd->nOpcCols; j++ )
			{
				if( gd->opcData[j][i] <= 0. )
				{
					fprintf( ioQQQ, " Illegal data value in %s\n",chFile.c_str() );
					fprintf( ioQQQ, " Line #%ld, value=%14.6e\n",dl,gd->opcData[j][i] );
					cdEXIT(EXIT_FAILURE);
				}
			}
			/* the data in the input file should be sorted on frequency, either
			 * strictly monotonically increasing or decreasing, check this here... */
			if( i == 1 )
			{
				sgn = sign3(gd->opcAnu[1]-gd->opcAnu[0]);
				if( sgn == 0 ) 
				{
					double dataVal = lgLogData ? log10(gd->opcAnu[1]) : gd->opcAnu[1];
					fprintf( ioQQQ, " Illegal value for frequency in %s\n",chFile.c_str() );
					fprintf( ioQQQ, " Line #%ld, freq=%14.6e\n",dl,dataVal );
					cdEXIT(EXIT_FAILURE);
				}
			}
			else if( i > 1 ) 
			{
				if( sign3(gd->opcAnu[i]-gd->opcAnu[i-1]) != sgn ) 
				{
					double dataVal = lgLogData ? log10(gd->opcAnu[i]) : gd->opcAnu[i];
					fprintf( ioQQQ, " Illegal value for frequency in %s\n",chFile.c_str());
					fprintf( ioQQQ, " Line #%ld, freq=%14.6e\n",dl,dataVal);
					cdEXIT(EXIT_FAILURE);
				}
			}
		}
		gd->nAxes = 1;
		break;
	case OPC_GREY:
	case OPC_PAH1:
	case OPC_PAH2N:
	case OPC_PAH2C:
	case OPC_PAH3N:
	case OPC_PAH3C:
		/* nothing much to be done here, the opacities
		 * will be calculated without any further data */
		gd->nAxes = 1;
		break;
	default:
		TotalInsanity();
	}

	fclose(io2);
	return;
}

/* construct optical constants for mixed grain using a specific EMT */
STATIC void mie_read_mix(/*@in@*/  const string& chFile,
						 /*@out@*/ grain_data *gd)
{
	emt_type EMTtype;
	long int dl, i, j, k, l, nelem, nMaterial, sumAxes;
	double maxIndex = DBL_MAX, minIndex = DBL_MAX, nAtoms, sum, sum2, wavHi, wavLo, wavStep;
	complex<double> eps_eff(-DBL_MAX,-DBL_MAX);
	string chLine, chWord;

	DEBUG_ENTRY( "mie_read_mix()" );

	FILE *io2 = open_data( chFile, "r" );

	dl = 0; /* line counter for input file */

	/* first read magic number */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&gd->magic,true,dl);
	if( gd->magic != MAGIC_MIX ) 
	{
		fprintf( ioQQQ, " Mixed grain file %s has obsolete magic number\n",chFile.c_str() );
		fprintf( ioQQQ, " I found %ld, but expected %ld on line #%ld\n",gd->magic,MAGIC_MIX,dl );
		fprintf( ioQQQ, " Please replace this file with an up to date version\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* default depletion of grain molecule */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_double(chFile,chLine,&gd->depl,true,dl);
	if( gd->depl > 1. ) 
	{
		fprintf( ioQQQ, " Illegal value for default depletion in %s\n",chFile.c_str() );
		fprintf( ioQQQ, " Line #%ld, depl=%14.6e\n",dl,gd->depl);
		cdEXIT(EXIT_FAILURE);
	}

	/* read number of different materials contained in this grain */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&nMaterial,true,dl);
	if( nMaterial < 2 ) 
	{
		fprintf( ioQQQ, " Illegal number of materials in mixed grain file %s\n",chFile.c_str() );
		fprintf( ioQQQ, " I found %ld on line #%ld\n",nMaterial,dl );
		fprintf( ioQQQ, " This number should be at least 2\n" );
		cdEXIT(EXIT_FAILURE);
	}

	vector<double> frac(nMaterial);
	vector<double> frac2(nMaterial);
	vector<grain_data> gdArr(nMaterial);

	sum = 0.;
	sum2 = 0.;
	sumAxes = 0;
	for( i=0; i < nMaterial; i++ )
	{
		string chFile2;

		/* each line contains relative fraction of volume occupied by each material,
		 * followed by the name of the refractive index file between double quotes */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&frac[i],true,dl);
		if( frac[i] <= 0. )
		{
			fprintf( ioQQQ, " Invalid volume fraction was found on line #%ld of file %s\n",dl,chFile.c_str() );
			fprintf( ioQQQ, " Please supply a positive value\n" );
			cdEXIT(EXIT_FAILURE);
		}

		sum += frac[i];

		auto pp = chLine.find( '\"' );
		if( pp != string::npos )
		{
			chFile2 = chLine.substr(++pp);
			pp = chFile2.find( '\"' );
			if( pp != string::npos )
				chFile2.erase(pp);
		}
		if( pp == string::npos )
		{
			fprintf( ioQQQ, " No pair of double quotes was found on line #%ld of file %s\n",dl,chFile.c_str() );
			fprintf( ioQQQ, " Please supply the refractive index file name between double quotes\n" );
			cdEXIT(EXIT_FAILURE);
		}

		mie_read_ocn( chFile2, &gdArr[i] );
		if( gdArr[i].rfiType != RFI_TABLE )
		{
			fprintf( ioQQQ, " Input error on line #%ld of file %s\n",dl,chFile.c_str() );
			fprintf( ioQQQ, " File %s is not of type RFI_TBL, this is illegal\n",chFile2.c_str() );
			cdEXIT(EXIT_FAILURE);
		}

		frac2[i] = ( gdArr[i].mol_weight > 0. ) ? frac[i] : 0.;
		sum2 += frac2[i];

		sumAxes += gdArr[i].nAxes;
	}

	/* read keyword to determine which EMT to use */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_word(chLine,chWord,true);

	if( chWord.find( "FA00" ) != string::npos )
		EMTtype = FARAFONOV00;
	else if( chWord.find( "ST95" ) != string::npos )
		EMTtype = STOGNIENKO95;
	else if( chWord.find( "BR35" ) != string::npos )
		EMTtype = BRUGGEMAN35;
	else if( chWord.find( "MG04" ) != string::npos )
		EMTtype = MAXWELL_GARNETT04;
	else
	{
		fprintf( ioQQQ, " Keyword not recognized in %s\n",chFile.c_str() );
		fprintf( ioQQQ, " Line #%ld, keyword=%s\n",dl,chWord.c_str());
		cdEXIT(EXIT_FAILURE);
	}

	if( EMTtype == FARAFONOV00 )
	{
		for( i=0; i < nMaterial; i++ )
		{
			if( gdArr[i].nAxes > 1 )
			{
				fprintf( ioQQQ, " ERROR: layers cannot be crystalline for the Farafonov EMT.\n" );
				cdEXIT(EXIT_FAILURE);
			}			
		}
	}
	if( EMTtype == FARAFONOV00 || EMTtype == STOGNIENKO95 || EMTtype == BRUGGEMAN35 )
	{
		if( gdArr[nMaterial-1].mol_weight == 0. )
		{
			fprintf( ioQQQ, " ERROR: last entry cannot be vacuum for the Farafonov, Stognienko, or Bruggeman EMT.\n" );
			fprintf( ioQQQ, " Please move this entry to an earlier position in the file\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	if( EMTtype == MAXWELL_GARNETT04 )
	{
		if( gdArr[0].mol_weight == 0. )
		{
			fprintf( ioQQQ, " ERROR: the matrix material cannot be vacuum for the Maxwell Garnett EMT.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		if( gdArr[0].nAxes > 1 )
		{
			fprintf( ioQQQ, " ERROR: the matrix material cannot be crystalline for the Maxwell Garnett EMT.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* normalize sum of fractional volumes to 1 */
	for( i=0; i < nMaterial; i++ )
	{
		frac[i] /= sum;
		frac2[i] /= sum2;
		/* renormalize elmAbun to chemical formula */
		for( nelem=0; nelem < LIMELM; nelem++ )
		{
			gdArr[i].elmAbun[nelem] /= gdArr[i].abun*gdArr[i].depl;
		}
	}

	wavLo = 0.;
	wavHi = DBL_MAX;
	gd->abun = DBL_MAX;
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		gd->elmAbun[nelem] = 0.;
	}
	gd->mol_weight = 0.;
	gd->rho = 0.;
	gd->work = DBL_MAX;
	gd->bandgap = DBL_MAX;
	gd->therm_eff = 0.;
	gd->subl_temp = DBL_MAX;
	gd->nAxes = 1;
	gd->wt[0] = 1.;
	gd->ndata[0] = MIX_TABLE_SIZE;
	gd->rfiType = RFI_TABLE;

	for( i=0; i < nMaterial; i++ )
	{
		for( k=0; k < gdArr[i].nAxes; k++ )
		{
			double wavMin = min(gdArr[i].wavlen[k][0],gdArr[i].wavlen[k][gdArr[i].ndata[k]-1]);
			double wavMax = max(gdArr[i].wavlen[k][0],gdArr[i].wavlen[k][gdArr[i].ndata[k]-1]);
			wavLo = max(wavLo,wavMin);
			wavHi = min(wavHi,wavMax);
		}
		minIndex = DBL_MAX;
		maxIndex = 0.;
		for( nelem=0; nelem < LIMELM; nelem++ )
		{
			gd->elmAbun[nelem] += frac[i]*gdArr[i].elmAbun[nelem];
			if( gd->elmAbun[nelem] > 0. )
			{
				minIndex = min(minIndex,gd->elmAbun[nelem]);
			}
			maxIndex = max(maxIndex,gd->elmAbun[nelem]);
		}
		gd->mol_weight += frac[i]*gdArr[i].mol_weight;
		gd->rho += frac[i]*gdArr[i].rho;
		/* ignore parameters for vacuum */
		if( gdArr[i].mol_weight > 0. )
		{
			gd->abun = min(gd->abun,gdArr[i].abun/frac[i]);
			switch( EMTtype )
			{
			case FARAFONOV00:
				/* this is appropriate for a layered grain */
				gd->work = gdArr[i].work;
				gd->bandgap = gdArr[i].bandgap;
				gd->therm_eff = gdArr[i].therm_eff;
				gd->matType = gdArr[i].matType;
				gd->subl_temp = min(gd->subl_temp,gdArr[i].subl_temp);
				break;
			case STOGNIENKO95:
			case BRUGGEMAN35:
				/* this is appropriate for a randomly mixed grain */
				gd->work = min(gd->work,gdArr[i].work);
				gd->bandgap = min(gd->bandgap,gdArr[i].bandgap);
				gd->therm_eff += frac2[i]*gdArr[i].therm_eff;
				gd->matType = gdArr[i].matType;
				gd->subl_temp = min(gd->subl_temp,gdArr[i].subl_temp);
				break;
			case MAXWELL_GARNETT04:
				/* this is appropriate for a matrix with some small inclusions, index 0 is the matrix */
				gd->work = gdArr[0].work;
				gd->bandgap = gdArr[0].bandgap;
				gd->therm_eff = gdArr[0].therm_eff;
				gd->matType = gdArr[0].matType;
				gd->subl_temp = gdArr[0].subl_temp;
				break;
			default:
				TotalInsanity();
			}
		}
	}

	if( gd->rho <= 0. )
	{
		fprintf( ioQQQ, " Illegal value for the density: %.3e\n", gd->rho );
		cdEXIT(EXIT_FAILURE);
	}
	if( gd->mol_weight <= 0. )
	{
		fprintf( ioQQQ, " Illegal value for the molecular weight: %.3e\n", gd->mol_weight );
		cdEXIT(EXIT_FAILURE);
	}
	if( maxIndex <= 0. )
	{
		fprintf( ioQQQ, " No atoms were found in the grain molecule\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* further sanity checks */
	ASSERT( wavLo > 0. && wavHi < DBL_MAX && wavLo < wavHi );
	ASSERT( gd->abun > 0. && gd->abun < DBL_MAX );
	ASSERT( gd->work > 0. && gd->work < DBL_MAX );
	ASSERT( gd->bandgap >= 0. && gd->bandgap < gd->work );
	ASSERT( gd->therm_eff > 0. && gd->therm_eff <= 1. );
	ASSERT( gd->subl_temp > 0. && gd->subl_temp < DBL_MAX );
	ASSERT( minIndex > 0. && minIndex < DBL_MAX );

	/* apply safety margin */
	wavLo *= 1. + 10.*DBL_EPSILON;
	wavHi *= 1. - 10.*DBL_EPSILON;

	/* renormalize the chemical formula such that the lowest index is 1 */
	nAtoms = 0.;
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		gd->elmAbun[nelem] /= minIndex;
		nAtoms += gd->elmAbun[nelem];
	}
	ASSERT( nAtoms > 0. );
	gd->abun *= minIndex;
	gd->mol_weight /= minIndex;
	/* calculate average weight per atom */
	gd->atom_weight = gd->mol_weight/nAtoms;

	mie_write_form(gd->elmAbun,chWord);
	fprintf( ioQQQ, "\n The chemical formula of the new grain molecule is: %s\n", chWord.c_str() );
	fprintf( ioQQQ, " The abundance wrt H at maximum depletion of this molecule is: %.3e\n",
		 gd->abun );
	fprintf( ioQQQ, " The abundance wrt H at standard depletion of this molecule is: %.3e\n",
		 gd->abun*gd->depl );

	/* finally renormalize elmAbun back to abundance relative to hydrogen */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		gd->elmAbun[nelem] *= gd->abun*gd->depl;
	}

	vector<double> delta(sumAxes);
	vector<double> frdelta(sumAxes);
	vector< complex<double> > eps(sumAxes);

	l = 0;
	for( i=0; i < nMaterial; i++ )
	{
		for( k=0; k < gdArr[i].nAxes; k++ )
		{
			frdelta[l] = gdArr[i].wt[k]*frac[i];
			delta[l] = ( l == 0 ) ? frdelta[l] : delta[l-1] + frdelta[l];
			++l;
		}
	}
	ASSERT( l == sumAxes && fabs(delta[l-1]-1.) < 10.*DBL_EPSILON );

	/* allocate space for wavelength and optical constants arrays */
	gd->wavlen[0].resize( gd->ndata[0] );
	gd->n[0].resize( gd->ndata[0] );
	gd->nr1[0].resize( gd->ndata[0] );

	wavStep = log(wavHi/wavLo)/(double)(gd->ndata[0]-1);

	switch( EMTtype )
	{
	case FARAFONOV00:
		/* this implements the EMT described in
		 * >>refer	grain	physics	Voshchinnikov N.V., Mathis J.S., 1999, ApJ, 526, 257
		 * based on the theory described in
		 * >>refer	grain	physics	Farafonov V.G., 2000, Optics & Spectroscopy, 88, 441
		 *
		 * NB - note that Eq. 3 in Voshchinnikov & Mathis is incorrect! */
		for( j=0; j < gd->ndata[0]; j++ )
		{
			double nre,nim;
			complex<double> a1,a2,a1c,a2c,a11,a12,a21,a22,ratio;

			gd->wavlen[0][j] = wavLo*exp((double)j*wavStep);

			init_eps(gd->wavlen[0][j],nMaterial,gdArr,eps);

			ratio = eps[0]/eps[1];

			a1 = (ratio+2.)/3.;
			a2 = (1.-ratio)*delta[0];

			for( l=1; l < sumAxes-1; l++ )
			{
				ratio = eps[l]/eps[l+1];

				a1c = a1;
				a2c = a2;
				a11 = (ratio+2.)/3.;
				a12 = (2.-2.*ratio)/(9.*delta[l]);
				a21 = (1.-ratio)*delta[l];
				a22 = (2.*ratio+1.)/3.;

				a1 = a11*a1c + a12*a2c;
				a2 = a21*a1c + a22*a2c;
			}

			a1c = a1;
			a2c = a2;
			a11 = 1.;
			a12 = 1./3.;
			a21 = eps[sumAxes-1];
			a22 = -2./3.*eps[sumAxes-1];

			a1 = a11*a1c + a12*a2c;
			a2 = a21*a1c + a22*a2c;

			ratio = a2/a1;
			dftori(&nre,&nim,ratio.real(),ratio.imag());

			gd->n[0][j] = complex<double>(nre,nim);
			gd->nr1[0][j] = nre-1.;
		}
		break;
	case STOGNIENKO95:
	case BRUGGEMAN35:
		for( j=0; j < gd->ndata[0]; j++ )
		{
			const double EPS_TOLER = 10.*DBL_EPSILON;
			double nre,nim;
			complex<double> eps0;

			gd->wavlen[0][j] = wavLo*exp((double)j*wavStep);

			init_eps(gd->wavlen[0][j],nMaterial,gdArr,eps);

			/* get initial estimate for effective dielectric function */
			if( j == 0 )
			{
				/* use simple average as first estimate */
				eps0 = 0.;
				for( l=0; l < sumAxes; l++ )
					eps0 += frdelta[l]*eps[l];
			}
			else
			{
				/* use solution from previous wavelength as first estimate */
				eps0 = eps_eff;
			}

			if( EMTtype == STOGNIENKO95 )
				/* this implements the EMT described in
				 * >>refer	grain	physics	Stognienko R., Henning Th., Ossenkopf V., 1995, A&A, 296, 797 */
				eps_eff = cnewton( Stognienko, frdelta, eps, sumAxes, eps0, EPS_TOLER );
			else if( EMTtype == BRUGGEMAN35 )
				/* this implements the classical Bruggeman rule
				 * >>refer	grain	physics	Bruggeman D.A.G., 1935, Ann. Phys. (5th series), 24, 636 */
				eps_eff = cnewton( Bruggeman, frdelta, eps, sumAxes, eps0, EPS_TOLER );
			else
				TotalInsanity();

			dftori(&nre,&nim,eps_eff.real(),eps_eff.imag());

			gd->n[0][j] = complex<double>(nre,nim);
			gd->nr1[0][j] = nre-1.;
		}
		break;
	case MAXWELL_GARNETT04:
		/* this implements the EMT described in
		 * >>refer	grain	physics	Maxwell Garnett, J.C., 1904, Philos. Trans. R. Soc. London, Ser. A 203, 385
		 *
		 * here we implement the generalization to spherical inclusions of multiple materials given in Eq. 26 of
		 * >>refer	grain	physics	Markel, V.A., 2016, JOSA A 33, 1244 */
		for( j=0; j < gd->ndata[0]; j++ )
		{
			gd->wavlen[0][j] = wavLo*exp((double)j*wavStep);

			init_eps(gd->wavlen[0][j],nMaterial,gdArr,eps);

			complex<double> sum = 0.;
			for( l=1; l < sumAxes; l++ )
				sum += frdelta[l]*(eps[l]-eps[0])/(eps[l]+2.*eps[0]);
			complex<double> epsMG = (1. + 2.*sum)/(1. - sum)*eps[0];

			double nre,nim;
			dftori(&nre,&nim,epsMG.real(),epsMG.imag());

			gd->n[0][j] = complex<double>(nre,nim);
			gd->nr1[0][j] = nre-1.;
		}
		break;
	default:
		TotalInsanity();
	}

	fclose(io2);
	return;
}

/* helper routine for mie_read_mix, initializes the array of dielectric functions */
STATIC void init_eps(double wavlen,
					 long nMaterial,
					 /*@in@*/ const vector<grain_data>& gdArr, /* gdArr[nMaterial] */
					 /*@out@*/ vector< complex<double> >& eps) /* eps[sumAxes] */
{
	DEBUG_ENTRY( "init_eps()" );

	long l = 0;
	for( long i=0; i < nMaterial; i++ )
	{
		for( long k=0; k < gdArr[i].nAxes; k++ )
		{
			bool lgErr;
			long ind;
			double eps1,eps2,frc,nim,nre;

			find_arr(wavlen,gdArr[i].wavlen[k],gdArr[i].ndata[k],&ind,&lgErr);
			ASSERT( !lgErr );
			frc = (wavlen-gdArr[i].wavlen[k][ind])/(gdArr[i].wavlen[k][ind+1]-gdArr[i].wavlen[k][ind]);
			ASSERT( frc > 0.-10.*DBL_EPSILON && frc < 1.+10.*DBL_EPSILON );
			nre = (1.-frc)*gdArr[i].n[k][ind].real() + frc*gdArr[i].n[k][ind+1].real();
			ASSERT( nre > 0. );
			nim = (1.-frc)*gdArr[i].n[k][ind].imag() + frc*gdArr[i].n[k][ind+1].imag();
			ASSERT( nim >= 0. );
			ritodf(nre,nim,&eps1,&eps2);
			eps[l++] = complex<double>(eps1,eps2);
		}
	}
	return;
}

/*******************************************************
 *  This routine is derived from the routine Znewton   *
 * --------------------------------------------------- *
 *  Reference; BASIC Scientific Subroutines, Vol. II   *
 *  by F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981.      *
 *                                                     *
 *              C++ version by J-P Moreau, Paris.      *
 ******************************************************/
/* find complex root of fun using the Newton-Raphson algorithm */
STATIC complex<double> cnewton(
	void(*fun)(complex<double>,const vector<double>&,const vector< complex<double> >&,
			   long,complex<double>*,double*,double*),
	/*@in@*/ const vector<double>& frdelta,        /* frdelta[sumAxes] */
	/*@in@*/ const vector< complex<double> >& eps, /* eps[sumAxes] */
	long sumAxes,
	complex<double> x0,
	double tol)
{
	const int LOOP_MAX = 100;
	const double TINY = 1.e-12;

	DEBUG_ENTRY( "cnewton()" );
	for( int i=0; i < LOOP_MAX; i++ )
	{
		complex<double> x1,y;
		double dudx,dudy,norm2;

		(*fun)(x0,frdelta,eps,sumAxes,&y,&dudx,&dudy);

		norm2 = pow2(dudx) + pow2(dudy);
		/* guard against norm2 == 0 */
		if( norm2 < TINY*norm(y) )
		{
			fprintf( ioQQQ, " cnewton - zero divide error\n" );
			cdEXIT(EXIT_FAILURE);
		}
		x1 = x0 - complex<double>( y.real()*dudx-y.imag()*dudy, y.imag()*dudx+y.real()*dudy )/norm2;

		/* check for convergence */
		if( fabs(x0.real()/x1.real()-1.) + fabs(x0.imag()/x1.imag()-1.) < tol )
		{
			return x1;
		}

		x0 = x1;
	}

	fprintf( ioQQQ, " cnewton did not converge\n" );
	cdEXIT(EXIT_FAILURE);
}

/* this evaluates the function defined in Eq. 3 of
 * >>refer	grain	physics	Stognienko R., Henning Th., Ossenkopf V., 1995, A&A, 296, 797
 * and its derivatives */
STATIC void Stognienko(complex<double> x,
					   /*@in@*/ const vector<double>& frdelta,        /* frdelta[sumAxes] */
					   /*@in@*/ const vector< complex<double> >& eps, /* eps[sumAxes] */
					   long sumAxes,
					   /*@out@*/ complex<double> *f,
					   /*@out@*/ double *dudx,
					   /*@out@*/ double *dudy)
{
	static const double L[4] = { 0., 1./2., 1., 1./3. };
	static const double fl[4] = { 5./9., 2./9., 2./9., 1. };

	DEBUG_ENTRY( "Stognienko()" );
	*f = complex<double>(0.,0.);
	*dudx = 0.;
	*dudy = 0.;
	for( long l=0; l < sumAxes; l++ )
	{
		complex<double> hlp = eps[l] - x;
		double h1 = eps[l].imag()*x.real() - eps[l].real()*x.imag();

		for( long i=0; i < 4; i++ )
		{
			double f1 = fl[i]*frdelta[l];
			double xx = ( i < 3 ) ? sin(PI*frdelta[l]) : cos(PI*frdelta[l]);
			complex<double> f2 = f1*xx*xx;
			complex<double> one = x + hlp*L[i];
			complex<double> two = f2*hlp/one;
			double h2 = norm(one);
			*f += two;
			*dudx -= f2.real()*(eps[l].real()*h2 + h1*2.*one.imag()*(1.-L[i]))/pow2(h2);
			*dudy -= f2.real()*(eps[l].imag()*h2 - h1*2.*one.real()*(1.-L[i]))/pow2(h2);
		}
	}
	return;
}

/* this evaluates the classical Bruggeman rule and its derivatives
 * >>refer	grain	physics	Bruggeman D.A.G., 1935, Ann. Phys. (5th series), 24, 636 */
STATIC void Bruggeman(complex<double> x,
					  /*@in@*/ const vector<double>& frdelta,        /* frdelta[sumAxes] */
					  /*@in@*/ const vector< complex<double> >& eps, /* eps[sumAxes] */
					  long sumAxes,
					  /*@out@*/ complex<double> *f,
					  /*@out@*/ double *dudx,
					  /*@out@*/ double *dudy)
{
	static const double L = 1./3.;

	DEBUG_ENTRY( "Bruggeman()" );
	*f = complex<double>(0.,0.);
	*dudx = 0.;
	*dudy = 0.;
	for( long l=0; l < sumAxes; l++ )
	{
		complex<double> hlp = eps[l] - x;
		double h1 = eps[l].imag()*x.real() - eps[l].real()*x.imag();
		complex<double> f2 = frdelta[l];
		complex<double> one = x + hlp*L;
		complex<double> two = f2*hlp/one;
		double h2 = norm(one);
		*f += two;
		*dudx -= f2.real()*(eps[l].real()*h2 + h1*2.*one.imag()*(1.-L))/pow2(h2);
		*dudy -= f2.real()*(eps[l].imag()*h2 - h1*2.*one.real()*(1.-L))/pow2(h2);
	}
	return;
}

/* read in the file with optical constants and other relevant information */
STATIC void mie_read_szd(/*@in@*/  const string& chFile,
						 /*@out@*/ sd_data *sd)
{
	bool lgTryOverride = false;
	long int dl,
	  i;
	double mul = 0.,
	  ref_neg = DBL_MAX,
	  ref_pos = DBL_MAX,
	  step_neg = DBL_MAX,
	  step_pos = DBL_MAX;
	string chLine,
	  chWord;
	FILE *io2;

	/* these constants are used to get initial estimates for the cutoffs (lim)
	 * in the SD_EXP_CUTOFFx and SD_xxx_NORMAL cases, they are iterated by
	 * search_limit such that
	 * lim^4 * dN/da(lim) == FRAC_CUTOFF * ref^4 * dN/da(ref)
	 * where ref as an appropriate reference point for each of the cases */
	const double FRAC_CUTOFF = 1.e-4;
	const double MUL_CO1 = -log(FRAC_CUTOFF);
	const double MUL_CO2 = sqrt(MUL_CO1);
	const double MUL_CO3 = cbrt(MUL_CO1);
	const double MUL_LND = sqrt(-2.*log(FRAC_CUTOFF));
	const double MUL_NRM = MUL_LND;

	DEBUG_ENTRY( "mie_read_szd()" );

	io2 = open_data( chFile, "r" );

	dl = 0; /* line counter for input file */

	/* first read magic number */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_long(chFile,chLine,&sd->magic,true,dl);
	if( sd->magic != MAGIC_SZD ) 
	{
		fprintf( ioQQQ, " Size distribution file %s has obsolete magic number\n",chFile.c_str() );
		fprintf( ioQQQ, " I found %ld, but expected %ld on line #%ld\n",sd->magic,MAGIC_SZD,dl );
		fprintf( ioQQQ, " Please replace this file with an up to date version\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* size distribution case */
	mie_next_data(chFile,io2,chLine,&dl);
	mie_read_word(chLine,chWord,true);

	if( chWord.find( "SSIZ" ) != string::npos )
	{
		sd->sdCase = SD_SINGLE_SIZE;
	}
	else if( chWord.find( "NCAR" ) != string::npos )
	{
		sd->sdCase = SD_NR_CARBON;
	}
	else if( chWord.find( "POWE" ) != string::npos )
	{
		sd->sdCase = SD_POWERLAW;
	}
	else if( chWord.find( "EXP1" ) != string::npos )
	{
		sd->sdCase = SD_EXP_CUTOFF1;
		sd->a[ipAlpha] = 1.;
		mul = MUL_CO1;
	}
	else if( chWord.find( "EXP2" ) != string::npos )
	{
		sd->sdCase = SD_EXP_CUTOFF2;
		sd->a[ipAlpha] = 2.;
		mul = MUL_CO2;
	}
	else if( chWord.find( "EXP3" ) != string::npos )
	{
		sd->sdCase = SD_EXP_CUTOFF3;
		sd->a[ipAlpha] = 3.;
		mul = MUL_CO3;
	}
	else if( chWord.find( "LOGN" ) != string::npos )
	{
		sd->sdCase = SD_LOG_NORMAL;
		mul = MUL_LND;
	}
	/* this one must come after LOGNORMAL */
	else if( chWord.find( "NORM" ) != string::npos )
	{
		sd->sdCase = SD_LIN_NORMAL;
		mul = MUL_NRM;
	}
	else if( chWord.find( "TABL" ) != string::npos )
	{
		sd->sdCase = SD_TABLE;
	}
	else
	{
		sd->sdCase = SD_ILLEGAL;
	}

	switch( sd->sdCase ) 
	{
	case SD_SINGLE_SIZE:
		/* single sized grain */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipSize],true,dl);
		if( sd->a[ipSize] < SMALLEST_GRAIN || sd->a[ipSize] > LARGEST_GRAIN ) 
		{
			fprintf( ioQQQ, " Illegal value for grain size\n" );
			fprintf( ioQQQ, " Grain sizes should be between %.5f and %.0f micron\n",
				 SMALLEST_GRAIN, LARGEST_GRAIN );
			fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
			cdEXIT(EXIT_FAILURE);
		}
		break;
	case SD_NR_CARBON:
		/* single sized PAH with fixed number of carbon atoms */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_long(chFile,chLine,&sd->nCarbon,true,dl);
		break;
	case SD_POWERLAW:
		/* simple power law distribution, first get lower limit */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipBLo],true,dl);
		if( sd->a[ipBLo] < SMALLEST_GRAIN || sd->a[ipBLo] > LARGEST_GRAIN ) 
		{
			fprintf( ioQQQ, " Illegal value for grain size (lower limit)\n" );
			fprintf( ioQQQ, " Grain sizes should be between %.5f and %.0f micron\n",
				 SMALLEST_GRAIN, LARGEST_GRAIN );
			fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
			cdEXIT(EXIT_FAILURE);
		}

		/* upper limit */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipBHi],true,dl);
		if( sd->a[ipBHi] < SMALLEST_GRAIN || sd->a[ipBHi] > LARGEST_GRAIN ||
		    sd->a[ipBHi] <= sd->a[ipBLo] ) 
		{
			fprintf( ioQQQ, " Illegal value for grain size (upper limit)\n" );
			fprintf( ioQQQ, " Grain sizes should be between %.5f and %.0f micron\n",
				 SMALLEST_GRAIN, LARGEST_GRAIN );
			fprintf( ioQQQ, " and upper limit should be greater than lower limit\n" );
			fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
			cdEXIT(EXIT_FAILURE);
		}

		/* slope */
		mie_next_data(chFile,io2,chLine,&dl);
		if( sscanf( chLine.c_str(), "%lf", &sd->a[ipExp] ) != 1 ) 
		{
			fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
			fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
			cdEXIT(EXIT_FAILURE);
		}

		sd->a[ipBeta] = 0.;
		sd->a[ipSLo] = 0.;
		sd->a[ipSHi] = 0.;

		sd->lim[ipBLo] = sd->a[ipBLo];
		sd->lim[ipBHi] = sd->a[ipBHi];
		break;
	case SD_EXP_CUTOFF1:
	case SD_EXP_CUTOFF2:
	case SD_EXP_CUTOFF3:
		/* powerlaw with first/second/third order exponential cutoff, inspired by
		 * Greenberg (1978), Cosmic Dust, ed. J.A.M. McDonnell, Wiley, p. 187 */
		/* "lower limit", below this the exponential cutoff sets in */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipBLo],false,dl);

		/* "upper" limit, above this the exponential cutoff sets in */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipBHi],false,dl);

		/* exponent for power law */
		mie_next_data(chFile,io2,chLine,&dl);
		if( sscanf( chLine.c_str(), "%lf", &sd->a[ipExp] ) != 1 ) 
		{
			fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
			fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
			cdEXIT(EXIT_FAILURE);
		}

		/* beta parameter, for extra curvature in the powerlaw region */
		mie_next_data(chFile,io2,chLine,&dl);
		if( sscanf( chLine.c_str(), "%lf", &sd->a[ipBeta] ) != 1 ) 
		{
			fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
			fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
			cdEXIT(EXIT_FAILURE);
		}

		/* scale size for lower exponential cutoff, zero indicates normal cutoff */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipSLo],false,dl);

		/* scale size for upper exponential cutoff, zero indicates normal cutoff */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipSHi],false,dl);

		ref_neg = sd->a[ipBLo];
		step_neg = -mul*sd->a[ipSLo];
		ref_pos = sd->a[ipBHi];
		step_pos = mul*sd->a[ipSHi];
		lgTryOverride = true;
		break;
	case SD_LOG_NORMAL:
		/* log-normal distribution, first get center of gaussian */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipGCen],true,dl);

		/* 1-sigma width */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipGSig],true,dl);

		/* ref_pos, ref_neg is the grain radius at which a^4*dN/da peaks */
		ref_neg = ref_pos = sd->a[ipGCen]*exp(3.*pow2(sd->a[ipGSig]));
		step_neg = sd->a[ipGCen]*(exp(-mul*sd->a[ipGSig]) - 1.);
		step_pos = sd->a[ipGCen]*(exp(mul*sd->a[ipGSig]) - 1.);
		lgTryOverride = true;
		break;
	case SD_LIN_NORMAL:
		/* normal gaussian distribution, first get center of gaussian */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipGCen],true,dl);

		/* 1-sigma width */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipGSig],true,dl);

		/* ref_pos, ref_neg is the grain radius at which a^4*dN/da peaks */
		ref_neg = ref_pos = (sd->a[ipGCen] + sqrt(pow2(sd->a[ipGCen]) + 12.*pow2(sd->a[ipGSig])))/2.;
		step_neg = -mul*sd->a[ipGSig];
		step_pos = mul*sd->a[ipGSig];
		lgTryOverride = true;
		break;
	case SD_TABLE:
		/* user-supplied table of a^4*dN/da vs. a, first get lower limit on a */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipBLo],true,dl);
		if( sd->a[ipBLo] < SMALLEST_GRAIN || sd->a[ipBLo] > LARGEST_GRAIN ) 
		{
			fprintf( ioQQQ, " Illegal value for grain size (lower limit)\n" );
			fprintf( ioQQQ, " Grain sizes should be between %.5f and %.0f micron\n",
				 SMALLEST_GRAIN, LARGEST_GRAIN );
			fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
			cdEXIT(EXIT_FAILURE);
		}

		/* upper limit */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&sd->a[ipBHi],true,dl);
		if( sd->a[ipBHi] < SMALLEST_GRAIN || sd->a[ipBHi] > LARGEST_GRAIN ||
		    sd->a[ipBHi] <= sd->a[ipBLo] ) 
		{
			fprintf( ioQQQ, " Illegal value for grain size (upper limit)\n" );
			fprintf( ioQQQ, " Grain sizes should be between %.5f and %.0f micron\n",
				 SMALLEST_GRAIN, LARGEST_GRAIN );
			fprintf( ioQQQ, " and upper limit should be greater than lower limit\n" );
			fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
			cdEXIT(EXIT_FAILURE);
		}

		/* number of user supplied points */
		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_long(chFile,chLine,&sd->npts,true,dl);
		if( sd->npts < 2 )
		{
			fprintf( ioQQQ, " Illegal value for no. of points in table\n" );
			fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
			cdEXIT(EXIT_FAILURE);
		}

		/* allocate space for the table */
		sd->ln_a.resize(sd->npts);
		sd->ln_a4dNda.resize(sd->npts);

		/* and read the table */
		for( i=0; i < sd->npts; ++i )
		{
			double help1, help2;

			mie_next_data(chFile,io2,chLine,&dl);
			/* read data pair: a (micron), a^4*dN/da (arbitrary units) */
			if( sscanf( chLine.c_str(), "%le %le", &help1, &help2 ) != 2 ) 
			{
				fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
				fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
				cdEXIT(EXIT_FAILURE);
			}

			if( help1 <= 0. || help2 <= 0. )
			{
				fprintf( ioQQQ, " Reading table failed on line #%ld of %s\n",dl,chFile.c_str() );
				fprintf( ioQQQ, " Illegal data value %.6e or %.6e\n", help1, help2 );
				cdEXIT(EXIT_FAILURE);
			}

			sd->ln_a[i] = log(help1);
			sd->ln_a4dNda[i] = log(help2);

			if( i > 0 && sd->ln_a[i] <= sd->ln_a[i-1] )
			{
				fprintf( ioQQQ, " Reading table failed on line #%ld of %s\n",dl,chFile.c_str() );
				fprintf( ioQQQ, " Grain radii should be monotonically increasing\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}

		sd->lim[ipBLo] = sd->a[ipBLo];
		sd->lim[ipBHi] = sd->a[ipBHi];
		break;
	case SD_ILLEGAL:
	default:
		fprintf( ioQQQ, " unimplemented case for grain size distribution in file %s\n", chFile.c_str() );
		fprintf( ioQQQ, " Line #%ld: value %s\n",dl,chWord.c_str());
		cdEXIT(EXIT_FAILURE);
	}

	/* >>chng 01 feb 12, use a^4*dN/da instead of dN/da to determine limits,
	 * this assures that upper limit gives negligible mass fraction, PvH */
	/* in all cases where search_limit is used to determine lim[ipBLo] and lim[ipBHi],
	 * the user can override these values in the last two lines of the size distribution
	 * file. these inputs are mandatory, and should be given in the sequence lower
	 * limit, upper limit. a value <= 0 indicates that search_limit should be used. */
	if( lgTryOverride )
	{
		double help;

		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&help,false,dl);
		sd->lim[ipBLo] = ( help <= 0. ) ? search_limit(ref_neg,step_neg,FRAC_CUTOFF,*sd) : help;

		mie_next_data(chFile,io2,chLine,&dl);
		mie_read_double(chFile,chLine,&help,false,dl);
		sd->lim[ipBHi] = ( help <= 0. ) ? search_limit(ref_pos,step_pos,FRAC_CUTOFF,*sd) : help;

		if( sd->lim[ipBLo] < SMALLEST_GRAIN || sd->lim[ipBHi] > LARGEST_GRAIN ||
		    sd->lim[ipBHi] <= sd->lim[ipBLo] ) 
		{
			fprintf( ioQQQ, " Illegal size limits: lower %.5f and/or upper %.5f\n",
				 sd->lim[ipBLo], sd->lim[ipBHi] );
			fprintf( ioQQQ, " Grain sizes should be between %.5f and %.0f micron\n",
				 SMALLEST_GRAIN, LARGEST_GRAIN );
			fprintf( ioQQQ, " and upper limit should be greater than lower limit\n" );
			fprintf( ioQQQ, " Please alter the size distribution file\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	fclose(io2);
	return;
}

STATIC void mie_read_long(/*@in@*/  const string& chFile,
						  /*@in@*/  const string& chLine,
						  /*@out@*/ long int *data,
						  bool lgZeroIllegal,
						  long int dl)
{
	DEBUG_ENTRY( "mie_read_long()" );
	if( sscanf( chLine.c_str(), "%ld", data ) != 1 )
	{
		fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
		fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
		cdEXIT(EXIT_FAILURE);
	}
	if( *data < 0 || (*data == 0 && lgZeroIllegal) )
	{
		fprintf( ioQQQ, " Illegal data value in %s\n",chFile.c_str());
		fprintf( ioQQQ, " Line #%ld: %ld\n",dl,*data);
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

STATIC void mie_read_realnum(/*@in@*/  const string& chFile,
							 /*@in@*/  const string& chLine,
							 /*@out@*/ realnum *data,
							 bool lgZeroIllegal,
							 long int dl)
{
	DEBUG_ENTRY( "mie_read_realnum()" );
	double help;
	if( sscanf( chLine.c_str(), "%lf", &help ) != 1 )
	{
		fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
		fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
		cdEXIT(EXIT_FAILURE);
	}
	*data = (realnum)help;
	if( *data < 0. || (*data == 0. && lgZeroIllegal) )
	{
		fprintf( ioQQQ, " Illegal data value in %s\n",chFile.c_str());
		fprintf( ioQQQ, " Line #%ld: %14.6e\n",dl,*data);
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

STATIC void mie_read_double(/*@in@*/  const string& chFile,
							/*@in@*/  const string& chLine,
							/*@out@*/ double *data,
							bool lgZeroIllegal,
							long int dl)
{
	DEBUG_ENTRY( "mie_read_double()" );
	if( sscanf( chLine.c_str(), "%lf", data ) != 1 )
	{
		fprintf( ioQQQ, " Syntax error in %s\n",chFile.c_str());
		fprintf( ioQQQ, " Line #%ld: %s\n",dl,chLine.c_str());
		cdEXIT(EXIT_FAILURE);
	}
	if( *data < 0. || (*data == 0. && lgZeroIllegal) )
	{
		fprintf( ioQQQ, " Illegal data value in %s\n",chFile.c_str());
		fprintf( ioQQQ, " Line #%ld: %14.6e\n",dl,*data);
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

STATIC void mie_read_form(/*@in@*/  const string& chWord,
						  /*@out@*/ double elmAbun[],    /* elmAbun[LIMELM] */
						  /*@out@*/ double *no_atoms,
						  /*@out@*/ double *mol_weight)
{
	DEBUG_ENTRY( "mie_read_form()" );

	*no_atoms = 0.;
	*mol_weight = 0.;
	for( long nelem=0; nelem < LIMELM; nelem++ ) 
	{
		double frac = 0.;
		string chElmName(elementnames.chElementSym[nelem]);
		if( chElmName[1] == ' ' )
			chElmName.pop_back();
		string::size_type ptr = chWord.find( chElmName );
		if( ptr != string::npos )
		{
			long len = chElmName.length();
			/* prevent spurious match, e.g. F matches Fe */
			if( !islower((unsigned char)chWord[ptr+len]) ) 
			{
				if( isdigit((unsigned char)chWord[ptr+len]) ) 
				{
					sscanf(chWord.c_str()+ptr+len,"%lf",&frac);
				}
				else 
				{
					frac = 1.;
				}
			}
		}
		elmAbun[nelem] = frac;
		/* >>chng 02 apr 22, don't count hydrogen in PAH's, PvH */
		if( nelem != ipHYDROGEN )
			*no_atoms += frac;
		*mol_weight += frac*dense.AtomicWeight[nelem];
	}
	/* prevent division by zero when no chemical formula was supplied */
	if( *no_atoms == 0. ) 
		*no_atoms = 1.;
	return;
}

STATIC void mie_write_form(/*@in@*/  const double elmAbun[], /* elmAbun[LIMELM] */
						   /*@out@*/ string& chWord)
{
	DEBUG_ENTRY( "mie_write_form()" );

	ostringstream oss;
	for( long nelem=0; nelem < LIMELM; nelem++ ) 
	{
		if( elmAbun[nelem] > 0. )
		{
			string chElmName(elementnames.chElementSym[nelem]);
			if( chElmName[1] == ' ' )
				chElmName.pop_back();

			long index100 = nint(100.*elmAbun[nelem]);
			if( index100 == 100 )
				oss << chElmName;
			else if( index100%100 == 0 )
				oss << chElmName << index100/100;
			else
			{
				double xIndex = (double)index100/100.;
				oss << chElmName << fixed << setprecision(2) << xIndex;
			}
		}
	}
	chWord = oss.str();
}

STATIC void mie_read_word(/*@in@*/  const string& chLine,
						  /*@out@*/ string& chWord,
						  bool lgToUpper)
{
	DEBUG_ENTRY( "mie_read_word()" );

	chWord.clear();
	/* skip leading spaces or double quotes */
	size_t ip = 0;
	while( ip < chLine.length() && ( chLine[ip] == ' ' || chLine[ip] == '"' ) )
		ip++;
	/* now copy string until we hit next space or double quote */
	while( ip < chLine.length() && chLine[ip] != ' ' && chLine[ip] != '"' )
		if( lgToUpper )
			chWord.push_back( toupper(chLine[ip++]) );
		else
			chWord.push_back( chLine[ip++] );
}

/*=====================================================================*/
STATIC void mie_next_data(/*@in@*/  const string& chFile,
						  /*@in@*/  FILE *io,
						  /*@out@*/ string& chLine,
						  /*@in@*/  long int *dl)
{
	DEBUG_ENTRY( "mie_next_data()" );

	/* lines starting with a pound sign are considered comments and are skipped,
	 * lines not starting with a pound sign are considered to contain useful data.
	 * however, comments may still be appended to the line and will be erased. */

	do
	{
		mie_next_line(chFile,io,chLine,dl);
	}
	while( chLine[0] == '#' );
	/* erase comment part of the line */
	size_t pp;
	if( (pp = chLine.find("#")) != string::npos )
		chLine.erase(pp);
}

/*=====================================================================*/
STATIC void mie_next_line(/*@in@*/  const string& chFile,
						  /*@in@*/  FILE *io,
						  /*@out@*/ string& chLine,
						  /*@in@*/  long int *dl)
{
	DEBUG_ENTRY( "mie_next_line()" );

	if( !read_whole_line( chLine, io ) ) 
	{
		fprintf( ioQQQ, " Could not read from %s\n",chFile.c_str());
		if( feof(io) )
			fprintf( ioQQQ, " EOF reached\n");
		fprintf( ioQQQ, " This grain data file does not have the expected format.\n");
		cdEXIT(EXIT_FAILURE);
	}
	/* erase EOL character */
	size_t pp;
	if( (pp = chLine.find_first_of("\n\r")) != string::npos )
		chLine.erase(pp);
	(*dl)++;
}

/*=====================================================================*
 *
 * The routines gauss_init and gauss_legendre were derived from the
 * program cmieuvx.f.
 *
 * Written by: P.G. Martin (CITA), based on the code described in
 * >>refer	grain	physics	Hansen, J. E., Travis, L. D. 1974, Space Sci. Rev., 16, 527
 *
 * The algorithm in gauss_legendre was modified by Peter van Hoof to
 * avoid FP overflow for large values of nn.
 *
 *=====================================================================*/
/* set up Gaussian quadrature for arbitrary interval */
void gauss_init(long int nn,
				double xbot,
				double xtop,
				const vector<double>& x, /* x[nn]  */
				const vector<double>& a, /* a[nn]  */
				vector<double>& rr,      /* rr[nn] */
				vector<double>& ww)      /* ww[nn] */
{
	long int i;
	double bma,
	  bpa;

	DEBUG_ENTRY( "gauss_init()" );

	bpa = (xtop+xbot)/2.;
	bma = (xtop-xbot)/2.;

	for( i=0; i < nn; i++ ) 
	{
		rr[i] = bpa + bma*x[nn-1-i];
		ww[i] = bma*a[i];
	}
	return;
}

/*=====================================================================*/
/* set up abscissas and weights for Gauss-Legendre intergration of arbitrary even order */
void gauss_legendre(long int nn,
					vector<double>& x, /* x[nn] */
					vector<double>& a) /* a[nn] */
{
	long int i,
	  iter,
	  j;
	double cc,
	  csa,
	  d,
	  dp1,
	  dpn = 0.,
	  dq,
	  fj,
	  fn,
	  pn,
	  pn1 = 0.,
	  q,
	  xt = 0.;

	const double SAFETY = 5.;


	DEBUG_ENTRY( "gauss_legendre()" );

	if( nn%2 == 1 ) 
	{
		fprintf( ioQQQ, " Illegal number of abcissas\n" );
		cdEXIT(EXIT_FAILURE);
	}

	vector<double> c(nn);

	fn = (double)nn;
	csa = 0.;
	cc = 2.;
	for( j=1; j < nn; j++ ) 
	{
		fj = (double)j;
		/* >>chng 01 apr 10, prevent underflows in cc, pn, pn1, dpn and dp1 for large nn
		 * renormalize c[j] -> 4*c[j],  cc -> 4^(nn-1)*cc,  hence cc = O(1), etc...
		 * Old code: c[j] = pow2(fj)/(4.*(fj-0.5)*(fj+0.5)); */
		c[j] = pow2(fj)/((fj-0.5)*(fj+0.5));
		cc *= c[j];
	}

	for( i=0; i < nn/2; i++ ) 
	{
		switch( i ) 
		{
		case 0:
			xt = 1. - 2.78/(4. + pow2(fn));
			break;
		case 1:
			xt = xt - 4.1*(1. + 0.06*(1. - 8./fn))*(1. - xt);
			break;
		case 2:
			xt = xt - 1.67*(1. + 0.22*(1. - 8./fn))*(x[0] - xt);
			break;
		default:
			xt = 3.*(x[i-1] - x[i-2]) + x[i-3];
		}
		d = 1.;
		for( iter=1; (iter < 20) && (fabs(d) > DBL_EPSILON); iter++ ) 
		{
			/* >>chng 01 apr 10, renormalize pn -> 2^(nn-1)*pn, dpn -> 2^(nn-1)*dpn
			 * pn1 -> 2^(nn-2)*pn1, dp1 -> 2^(nn-2)*dp1
			 * Old code: pn1 = 1.; */
			pn1 = 0.5;
			pn = xt;
			dp1 = 0.;
			dpn = 1.;
			for( j=1; j < nn; j++ )
			{
				/* >>chng 01 apr 10, renormalize pn -> 2^(nn-1)*pn, dpn -> 2^(nn-1)*dpn
				 * Old code: q = xt*pn - c[j]*pn1;  dq = xt*dpn - c[j]*dp1 + pn; */
				q = 2.*xt*pn - c[j]*pn1;
				dq = 2.*xt*dpn - c[j]*dp1 + 2.*pn;
				pn1 = pn;
				pn = q;
				dp1 = dpn;
				dpn = dq;
			}
			d = pn/dpn;
			xt -= d;
		}
		x[i] = xt;
		x[nn-1-i] = -xt;
		/* >>chng 01 apr 10, renormalize dpn -> 2^(nn-1)*dpn, pn1 -> 2^(nn-2)*pn1
		 * Old code: a[i] = cc/(dpn*pn1); */
		a[i] = cc/(dpn*2.*pn1);
		a[nn-1-i] = a[i];
		csa += a[i];
	}

	/* this routine has been tested for every even nn between 2 and 4096
	 * it passed the test for each of those cases with SAFETY < 3.11 */
	if( fabs(1.-csa) > SAFETY*fn*DBL_EPSILON ) 
	{
		fprintf( ioQQQ, " gauss_legendre failed to converge: delta = %.4e\n", fabs(1.-csa) );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}

/* find index ind such that min(xa[ind],xa[ind+1]) <= x <= max(xa[ind],xa[ind+1]).
 * xa is assumed to be strictly monotically increasing or decreasing.
 * if x is outside the range spanned by xa, lgOutOfBounds is raised and ind is set to -1
 * n is the number of elements in xa. */
void find_arr(double x,
			  const vector<double>& xa,
			  long int n,
			  /*@out@*/ long int *ind,
			  /*@out@*/ bool *lgOutOfBounds)
{
	long int i1,
		i2,
		i3,
		sgn,
		sgn2;

	DEBUG_ENTRY( "find_arr()" );
	/* this routine works for strictly monotically increasing
	 * and decreasing arrays, sgn indicates which case it is */
	if( n < 2 ) 
	{
		fprintf( ioQQQ, " Invalid array\n");
		cdEXIT(EXIT_FAILURE);
	}

	i1 = 0;
	i3 = n-1;
	sgn = sign3(xa[i3]-xa[i1]);
	if( sgn == 0 ) 
	{
		fprintf( ioQQQ, " Ill-ordered array\n");
		cdEXIT(EXIT_FAILURE);
	}

	*lgOutOfBounds = x < min(xa[0],xa[n-1]) || x > max(xa[0],xa[n-1]);
	if( *lgOutOfBounds ) 
	{
		*ind = -1;
		return;
	}

	i2 = (n-1)/2;
	while( (i3-i1) > 1 ) 
	{
		sgn2 = sign3(x-xa[i2]);
		if( sgn2 != 0 )
		{
			if( sgn == sgn2 ) 
			{
				i1 = i2;
			}
			else 
			{
				i3 = i2;
			}
			i2 = (i1+i3)/2;
		}
		else 
		{
			*ind = i2;
			return;
		}
	}
	*ind = i1;
	return;
}

/*=====================================================================*
 *
 * The routines sinpar, anomal, bigk, ritodf, and dftori were derived
 * from the program cmieuvx.f.
 *
 * Written by: P.G. Martin (CITA), based on the code described in
 * >>refer	grain	physics	Hansen, J. E., Travis, L. D. 1974, Space Sci. Rev., 16, 527
 *
 *=====================================================================*/

/* Oct 1988 for UV - X-ray extinction, including anomalous diffraction check
 *     this version reads in real and imaginary parts of the refractive
 *     index, with imaginary part positive (nridf = 3) or nr-1 (nridf = 2) or
 *     real and imaginary parts of the dielectric function (nridf = 1)
 * Dec 1988: added qback; approximation for small x;
 * qphase, better convergence checking
 *
 * in anomalous diffraction: qext and qabs calculated - qscatt by subtraction
 * in rayleigh-gans:         qscatt and qabs calculated
 * in mie:                   qext and qscatt calculated
 * */

/* sinpar.f
 * consistency checks updated july 1999
 * t1 updated mildly 19 oct 1992
 * utility for mieuvx.f and mieuvxsd.f */
static const int NMXLIM = 80000;

STATIC void sinpar(double nre,
				   double nim,
				   double x,
				   /*@out@*/ double *qext,
				   /*@out@*/ double *qphase,
				   /*@out@*/ double *qscat,
				   /*@out@*/ double *ctbrqs,
				   /*@out@*/ double *qback,
				   /*@out@*/ long int *iflag)
{
	long int n,
	  nmx1,
	  nmx2,
	  nn,
	  nsqbk;
	double ectb,
	  eqext,
	  eqpha,
	  eqscat,
	  error=0.,
	  error1=0.,
	  rx,
	  t1,
	  t2,
	  t3,
	  t4,
	  t5,
	  tx,
	  x3,
	  x5=0.,
	  xcut,
	  xrd;
	complex<double> cdum1,
	  cdum2,
	  ci,
	  eqb,
	  nc,
	  nc2,
	  nc212,
	  qbck,
	  rrf,
	  rrfx,
	  sman,
	  sman1,
	  smbn,
	  smbn1,
	  tc1,
	  tc2,
	  wn,
	  wn1,
	  wn2;

	DEBUG_ENTRY( "sinpar()" );

	*iflag = 0;
	ci = complex<double>(0.,1.);
	nc = complex<double>(nre,-nim);
	nc2 = nc*nc;
	rrf = 1./nc;
	rx = 1./x;
	rrfx = rrf*rx;

	/*  t1 is the number of terms nmx2 that will be needed to obtain convergence
	 *  try to minimize this, because the a(n) downwards recursion has to
	 *  start at nmx1 larger than this
	 *
	 * major loop series is summed to nmx2, or less when converged
	 * nmx1 is used for a(n) only, n up to nmx2.
	 * must start evaluation sufficiently above nmx2 that a(nmx1)=(0.,0.)
	 * is a good approximation
	 *
	 *
	 *orig with slight modification for extreme UV and X-ray, n near 1., large x
	 *orig      t1=x*dmax1( 1.1d0,dsqrt(nr*nr+ni*ni) )*
	 *orig     1(1.d0+0.02d0*dmax1(dexp(-x/100.d0)*x/10.d0,dlog10(x)))
	 *
	 * rules like those of wiscombe 1980 are slightly more efficient */
	xrd = cbrt(x);
	/* the final number in t1 was 1., 2. for large x, and 3. is needed sometimes
	 * see also idnint use below */
	t1 = x + 4.*xrd + 3.;
	/*      t1=t1+0.05d0*xrd
	 * was 0., then 1., then 2., now 3. for intermediate x
	 * 19 oct 1992 */
	if( !(x <= 8. || x >= 4200.) )
		t1 += 0.05*xrd + 3.;
	t1 *= 1.01;

	/* the original rule of dave for starting the downwards recursion was
	 * to start at 1.1*|mx| + 1, i.e. with the original version of t1
	 *orig      nmx1=1.10d0*t1
	 *
	 * try a simpler, less costly one, as in bohren and huffman, p 478
	 * this is the form for use with wiscombe rules for t1
	 * tests: it produces the same results as the more costly version
	 * */
	t4 = x*sqrt(nre*nre+nim*nim);
	nmx1 = nint(max(t1,t4)) + 15;

	if( nmx1 < NMXLIM ) 
	{
		nmx2 = nint(t1);
		/*orig      if( nmx1  .gt. 150 ) go to 22
		 *orig      nmx1 = 150
		 *orig      nmx2 = 135
		 *
		 * try a more efficient scheme */
		if( nmx2 <= 4 ) 
		{
			nmx2 = 4;
			nmx1 = nint(max(4.,t4)) + 15;
		}

		vector< complex<double> > a(nmx1+1);

		/* downwards recursion for logarithmic derivative */
		a[nmx1] = 0.;

		/* note that with the method of lentz 1976 (appl opt 15, 668), it would be
		 * possible to find a(nmx2) directly, and start the downwards recursion there
		 * however, there is not much in it with above form for nmx1 which uses just */
		for( n=0; n < nmx1; n++ )
		{
			nn = nmx1 - n;
			a[nn-1] = (double)(nn+1)*rrfx - 1./((double)(nn+1)*rrfx+a[nn]);
		}

		sincos(x,&t2,&t1);
		wn2 = complex<double>(t1,-t2);
		wn1 = complex<double>(t2,t1);
		wn = rx*wn1 - wn2;
		tc1 = a[0]*rrf + rx;
		tc2 = a[0]*nc + rx;
		sman = (tc1*wn.real() - wn1.real())/(tc1*wn - wn1);
		smbn = (tc2*wn.real() - wn1.real())/(tc2*wn - wn1);

		/* small x; above calculations subject to rounding errors
		 * see bohren and huffman p 131
		 * wiscombe 1980 appl opt 19, 1505 gives alternative formulation */
		xcut = 3.e-04;
		if( x < xcut ) 
		{
			nc212 = (nc2-1.)/(nc2+2.);
			x3 = pow3(x);
			x5 = x3*pow2(x);
			/* note change sign convention for m = n - ik here */
			sman = ci*2.*x3*nc212*(1./3.+x*x*0.2*(nc2-2.)/(nc2+2.)) + 4.*x5*x*nc212*nc212/9.;
			smbn = ci*x5*(nc2-1.)/45.;
		}

		sman1 = sman;
		smbn1 = smbn;
		t1 = 1.5;
		sman *= t1;
		smbn *= t1;
		/* case n=1; note previous multiplication of sman and smbn by t1=1.5 */
		*qext = 2.*(sman.real() + smbn.real());
		*qphase = 2.*(sman.imag() + smbn.imag());
		nsqbk = -1;
		qbck = -2.*(sman - smbn);
		*qscat = (norm(sman) + norm(smbn))/.75;

		*ctbrqs = 0.0;
		n = 2;

		/************************* Major loop begins here ************************/
		while( true ) 
		{
			t1 = 2.*(double)n - 1.;
			t3 = 2.*(double)n + 1.;
			wn2 = wn1;
			wn1 = wn;
			wn = t1*rx*wn1 - wn2;
			cdum1 = a[n-1];
			cdum2 = n*rx;
			tc1 = cdum1*rrf + cdum2;
			tc2 = cdum1*nc + cdum2;
			sman = (tc1*wn.real() - wn1.real())/(tc1*wn - wn1);
			smbn = (tc2*wn.real() - wn1.real())/(tc2*wn - wn1);

			/* small x, n=2
			 * see bohren and huffman p 131 */
			if( x < xcut && n == 2 ) 
			{
				/* note change sign convention for m = n - ik here */
				sman = ci*x5*(nc2-1.)/(15.*(2.*nc2+3.));
				smbn = 0.;
			}

			eqext = t3*(sman.real() + smbn.real());
			*qext += eqext;
			eqpha = t3*(sman.imag() + smbn.imag());
			*qphase += eqpha;
			nsqbk = -nsqbk;
			eqb = t3*(sman - smbn)*(double)nsqbk;
			qbck += eqb;
			tx = norm(sman) + norm(smbn);
			eqscat = t3*tx;
			*qscat += eqscat;
			t2 = (double)(n - 1);
			t5 = (double)n;
			t4 = t1/(t5*t2);
			t2 = (t2*(t5 + 1.))/t5;
			ectb = t2*(sman1.real()*sman.real()+sman1.imag()*sman.imag() + smbn1.real()*smbn.real() +
				   smbn1.imag()*smbn.imag()) +
				t4*(sman1.real()*smbn1.real()+sman1.imag()*smbn1.imag());
			*ctbrqs += ectb;

			/* check convergence
			 * could decrease for large x and small m-1 in UV - X-ray; probably negligible */
			if( tx < 1.e-14 ) 
			{
				/* looks good but check relative convergence */
				eqext = fabs(eqext/ *qext);
				eqpha = fabs(eqpha/ *qphase);
				eqscat = fabs(eqscat/ *qscat);
				ectb = ( n == 2 ) ? 0. : fabs(ectb/ *ctbrqs);
				eqb = complex<double>( fabs(eqb.real()/qbck.real()), fabs(eqb.imag()/qbck.imag()) );
				/* leave out eqb.re/im, which are sometimes least well converged */
				error = MAX4(eqext,eqpha,eqscat,ectb);
				/* put a milder constraint on eqb.re/im */
				error1 = max(eqb.real(),eqb.imag());
				if( error < 1.e-07 && error1 < 1.e-04 )
					break;

				/* not sufficiently converged
				 *
				 * cut out after n=2 for small x, since approximation is being used */
				if( x < xcut )
					break;
			}

			smbn1 = smbn;
			sman1 = sman;
			n++;
			if( n > nmx2 ) 
			{
				*iflag = 1;
				break;
			}
		}
		/* renormalize */
		t1 = 2.*pow2(rx);
		*qext *= t1;
		*qphase *= t1;
		*qback = norm(qbck)*pow2(rx);
		*qscat *= t1;
		*ctbrqs *= 2.*t1;
	}
	else 
	{
		*iflag = 2;
	}
	return;
}

STATIC void anomal(double x,
				   /*@out@*/ double *qext,
				   /*@out@*/ double *qabs,
				   /*@out@*/ double *qphase,
				   /*@out@*/ double *xistar,
				   double delta,
				   double beta)
{
	/*
	 *
	 * in anomalous diffraction: qext and qabs calculated - qscatt by subtraction
	 * in rayleigh-gans:         qscatt and qabs calculated
	 * in mie:                   qext and qscatt calculated
	 *
	 */
	double xi,
	  xii;
	complex<double> cbigk,
	  ci,
	  cw;

	DEBUG_ENTRY( "anomal()" );
	/* anomalous diffraction: x>>1 and |m-1|<<1, any xi,xii
	 * original approach see Martin 1970. MN 149, 221 */
	xi = 2.*x*delta;
	xii = 2.*x*beta;
	/* xistar small is the basis for rayleigh-gans, any x, m-1 */
	*xistar = sqrt(pow2(xi)+pow2(xii));
	/* alternative approach see martin 1978 p 23 */
	ci = complex<double>(0.,1.);
	cw = -complex<double>(xi,xii)*ci;
	bigk(cw,&cbigk);
	*qext = 4.*cbigk.real();
	*qphase = 4.*cbigk.imag();
	cw = 2.*xii;
	bigk(cw,&cbigk);
	*qabs = 2.*cbigk.real();
	/* ?? put g in here - analytic version not known */
	return;
}

STATIC void bigk(complex<double> cw,
				 /*@out@*/ complex<double> *cbigk)
{
	/*
	 * see martin 1978 p 23
	 */

	DEBUG_ENTRY( "bigk()" );
	/* non-vax; use generic function */
	if( abs(cw) < 1.e-2 ) 
	{
		/* avoid severe loss of precision for small cw; expand exponential
		 * coefficients are 1/n! - 1/(n+1)! = 1/(n+1)(n-1)!;n=2,3,4,5,6,7
		 * accurate to  (1+ order cw**6) */
		*cbigk = cw*((1./3.)-cw*((1./8.)-cw*((1./30.)-cw*((1./144.)-cw*((1./840.)-cw*(1./5760.))))));
	}
	else 
	{
		*cbigk = 0.5 + (exp(-cw)*(1.+cw)-1.)/(cw*cw);
	}
	return;
}

/* utility for use with mieuvx/sd */
STATIC void ritodf(double nr,
				   double ni,
				   /*@out@*/ double *eps1,
				   /*@out@*/ double *eps2)
{
	DEBUG_ENTRY( "ritodf()" );
	/* refractive index to dielectric function */
	*eps1 = nr*nr - ni*ni;
	*eps2 = 2.*nr*ni;
	return;
}

/* utility for use with mieuvx/sd */
STATIC void dftori(/*@out@*/ double *nr,
				   /*@out@*/ double *ni,
				   double eps1,
				   double eps2)
{
	double eps;

	DEBUG_ENTRY( "dftori()" );
	/* dielectric function to refractive index  */
	eps = sqrt(eps2*eps2+eps1*eps1);
	*nr = sqrt((eps+eps1)/2.);
	ASSERT( *nr > 0. );
	/* >>chng 03 jan 02, old expression for ni suffered
	 * from cancellation error in the X-ray regime, PvH */
	/* *ni = sqrt((eps-eps1)/2.); */
	*ni = eps2/(2.*(*nr));
	return;
}
