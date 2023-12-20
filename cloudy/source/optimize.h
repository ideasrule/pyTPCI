/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef OPTIMIZE_H_
#define OPTIMIZE_H_

#include "energy.h"
#include "flux.h"
#include "lines.h"

/** called instead of cloudy to do optimize routine */
void optimize_do();

/**vary_input sets input lines to feed into cloudy in optimization runs 
\param *lgLimOK
*/
void vary_input(bool *lgLimOK, int grid_index);

/** optimize_subplex is the main driver, and only exposed, routine for the
 * cowan downhill simplex routine 
 \param n
 \param tol
 \param maxnfe
 \param mode
 \param scale[]
 \param x[]
 \param *fx
 \param *nfe
 \param work[]
 \param iwork[]
 \param *iflag
 */
void optimize_subplex(long int n, 
	  double tol, 
	  long int maxnfe, 
	  long int mode, 
	  realnum scale[], 
	  realnum x[], 
	  realnum *fx, 
	  long int *nfe, 
	  realnum work[], 
	  long int iwork[], 
	  long int *iflag);

const realnum VRSNEW = realnum(5.00);

typedef double chi2_type;

const chi2_type BIG_CHI2 = chi2_type(FLT_MAX);

void optimize_phymir(realnum[],const realnum[],long,chi2_type*,realnum);

/**optimize_func actual function called during evaluation of optimization run */
chi2_type optimize_func(const realnum param[], int grid_index = -1);

/* varypar.h */
/** the limit to the number of numbers on the command line */
const long LIMEXT = 5L;
const long LIMPAR = 20L;

extern string chOptimFileName; /* name of output file from optimization run */

typedef enum { PHYMIR_ILL, PHYMIR_SEQ, PHYMIR_FORK, PHYMIR_MPI } phymir_mode;

// this struct is written straight to a binary file, so it should be
// self-contained and NOT contain any pointers to other data!
// this also implies that only POD types are allowed here
template<class X, class Y = X, int NP = 32, int NSTR = 32>
class phymir_state
{
	X p_xmax;
	Y p_ymax;
	X p_xp[2*NP+1][NP];
	Y p_yp[2*NP+1];
	X p_absmin[NP];
	X p_absmax[NP];
	X p_varmin[NP];
	X p_varmax[NP];
	X p_a2[NP][NP];
	X p_c1[NP];
	X p_c2[NP];
	X p_xc[NP];
	X p_xcold[NP];
	X p_vers;
	X p_toler;
	X p_dmax;
	X p_dold;
	Y p_ymin;
	int32 p_dim;
	int32 p_sdim;
	int32 p_nvar;
	int32 p_noptim;
	int32 p_maxiter;
	int32 p_jmin;
	int32 p_maxcpu;
	int32 p_curcpu;
	phymir_mode p_mode;
	char p_chState[NSTR];
	char p_chStr1[NSTR];
	char p_chStr2[NSTR];
	char p_chStr3[NSTR];

	// ONLY data items BEFORE this item will be written to the state file
	Y (*p_func)(const X[],int);

	// private methods
	void p_clear1();
	void p_wr_state(const char*) const;
	void p_rd_state(const char*);
	Y p_execute_job( const X[], int, int );
	// p_execute_job_parallel MUST be const, otherwise changes by child processes may be lost!
	void p_execute_job_parallel( const X[], int, int ) const;
	void p_barrier( int, int );
	void p_process_output(int, int );
	void p_evaluate_hyperblock();
	void p_setup_next_hyperblock();
	void p_reset_hyperblock();
	void p_phygrm( X[][NP], int );
	bool p_lgLimitExceeded(const X[]) const;
	X p_delta( int i, int j ) const { return ( i == j ) ? X(1.) : X(0.); }
	void p_reset_transformation_matrix();

public:
	// public methods, there should be no data here...
	phymir_state() { p_clear1(); }
	void clear() { p_clear1(); }
	void init_minmax( const X[], const X[], int );
	void init_state_file_name( const char* );
	void init_strings( const string&, const string&, const char* );
	void initial_run( Y (*)(const X[],int), int, const X[], const X[], X, int, phymir_mode, int );
	void continue_from_state( Y (*)(const X[],int), int, const char*, X, int, phymir_mode, int );
	void optimize();
	void optimize_with_restart();
	bool lgMaxIterExceeded() const { return ( p_noptim >= p_maxiter ); }
	bool lgInitialized() const { return ( p_nvar > 0 ); }
	bool lgConverged() const { return ( p_dmax <= p_toler ); }
	bool lgConvergedRestart() const;
	X xval( int i ) const { return p_xc[i]; }
	X xmin( int i ) const { return max(p_varmin[i],p_absmin[i]); }
	X xmax( int i ) const { return min(p_varmax[i],p_absmax[i]); }
	Y yval() const { return p_ymin; }
	int32 noptim() const { return p_noptim; }
};


/**
 *
 *logical variable says whether current line image has vary option
 *
 *optimize increments, deltas for changing optimized variables
 *range for optimize command
 *io unit for final best parameters from optimizer
 *limit to number of iterations for optimizer, set with 
 *optimize iterations command
 *set with optimize tolerance command, used for global match to fit
 * default set in scalar to 0.10
 *current counter for the number of calls to the optimizer<BR>
 *lgTrOpt flag set with optimization trace command<BR>
 *nTrOpt is which call to cloudy to turn on trace<BR>
 *flags set if we are to optimize lines, luminosity, or colums<BR>
 *
 *
 *     labels for column densities on vary command<BR>
 *     this specifies the optimization routine<BR>
 *     'amoe', 'powe', 'bubr'
 */
struct t_optimize {

	/** flag set true if vary option used */
	bool lgVaryOn;
	/** flag set true if no vary command entered */
	bool lgNoVary;

	/** flag set true if optimize command entered */
	bool lgOptimr;

	/** flag to indicate that we are doing the initial parsing in grid_do()
	 *  rather than parsing an individual model during the optimization */
	bool lgInitialParse;

	bool lgOptimizeAsLinear[LIMPAR];

	/** min and max of range of variation */
	realnum varmax[LIMPAR];
	realnum varmin[LIMPAR];

	/** vparm is the value of the parameters on the line ,
	 * there can be as many as LIMEXT parameters - only first is varied */
	realnum vparm[LIMEXT][LIMPAR];

	/** the increment */
	realnum vincr[LIMPAR];

	/** the number of numbers on the command line */
	long int nvarxt[LIMPAR], 
	  nvfpnt[LIMPAR];

	realnum OptIncrm[LIMPAR],
	  varang[LIMPAR][2];

	long int nvary, 
	  nparm,
	  nRangeSet;
	bool lgVarOn;

	long int nIterOptim;

	/** parameters for the OPTIMIZE COLUMN DENSITY command */
	vector<string> chColDen_label;
	vector<long> ion_ColDen; 
	vector<realnum> ColDen_Obs;         
	vector<realnum> ColDen_error;

	/** parameters for the OPTIMIZE LINES command */
	int nEmergent;
	vector<LineID> lineids;
	/** error on the wavelength */
	vector<realnum> errorwave;
	vector<long> ipobs;
	vector<realnum> xLineInt_Obs; 
	vector<realnum> xLineInt_error;

	/** parameters for the OPTIMIZE TEMPERATURE command */
	vector<string> chTempLab;
	vector<long> ionTemp;
	vector<realnum> temp_obs;
	vector<realnum> temp_error;
	vector<string> chTempWeight;

	/** parameters for the OPTIMIZE DIAMETER command */
	bool lgOptDiam;
	bool lgDiamInCM;
	chi2_type optDiam;
	chi2_type optDiamErr;

	/** parameters for the OPTIMIZE CONTINUUM FLUX command */
	vector<long> ContIndex;
	vector<Energy> ContEner;
	vector<Flux> ContNFnu;
	vector<chi2_type> ContNFnuErr;

	realnum OptGlobalErr;

	/** counter for number of models in a grid - zero for first simulation */
	long int nOptimiz;

	bool lgOptimFlow;
	realnum optint;
	realnum optier;
	long int nTrOpt;
	bool lgTrOpt;
	bool lgOptimize;
	/** this flag says we are optimizing on luminosity */
	bool lgOptLum;
	int nOptLum;

	/** the following is needed by PHYMIR */
	bool lgParallel;
	bool lgOptCont;
	long useCPU;

	char chVarFmt[LIMPAR][FILENAME_PATH_LENGTH_2];
	char chOptRtn[5];

	double SavGenericData[10];

	t_optimize()
	{
		for( long i=0; i < LIMPAR; ++i )
			lgOptimizeAsLinear[i] = false;
	}
};
extern t_optimize optimize;

#endif /* OPTIMIZE_H_ */
