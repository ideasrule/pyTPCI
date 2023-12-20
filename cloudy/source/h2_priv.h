/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef H2_PRIV_H_
#define H2_PRIV_H_

#include "transition.h"
#include "container_classes.h"

typedef void (*linefunc)(const TransitionProxy& t, 
								 bool lgShield_this_zone,
								 realnum pestrk,
								 realnum DopplerWidth);
class Parser;
class molecule;
class CollRateCoeffArray;

/** the number of different types of colliders
*  */
const int N_X_COLLIDER = 5;
/** labels for the colliders */
const int chN_X_COLLIDER = 10;
// colliders are (see h2.cpp diatoms_init)
// 0 H
// 1 He
// 2 H2 ortho
// 3 H2 para
// 4 H+

/** the number of temperature points in the data file */
const int nTE_HMINUS = 7;
		
/** this is the number of electronic levels */
const int N_ELEC = 7;

/** log10 of temperatures where H- distribution are set */
const realnum H2_logte_hminus[nTE_HMINUS] = {1.,1.47712,2.,2.47712,3.,3.47712,4.};

struct diss_level 
{
	long n, v, j;
};

class diss_tran 
{
public:
	explicit diss_tran( diss_level a, diss_level b )
	{
		initial = a;
		final = b;
		energies.clear();
		xsections.clear();
		rate_coeff = 0.;
	};
	diss_level initial, final;
	vector<double> energies;	
	vector<double> xsections;
	double rate_coeff;	
};

struct t_coll_source
{
	t_coll_source()
	{
		magic = 0;
		filename = "";
	};
	long magic;
	string filename;
};

class diatomics
{
public:
	double Abund() const
	{
		return *dense_total;
	}
	void GetIndices( long& ipHi, long& ipLo, const char* chLine, long& i ) const;
	void CalcPhotoionizationRate(void);
	double (*photoion_opacity_fun)( double energy );
	long OpacityCreate( vector<double>& stack );
	double GetHeatRate( const diss_tran& tran );
	double GetDissociationRate( const diss_tran& tran );
	double MolDissocOpacity( const diss_tran& tran, const double& Mol_Ene );
	double Cont_Diss_Heat_Rate( void );
	void Mol_Photo_Diss_Rates( void );
	void Read_Mol_Diss_cross_sections(void);
	void SolveExcitedElectronicLevels(void);
	void SolveSomeGroundElectronicLevels(void);
	double GetExcitedElecDensity(void);
	realnum GetXColden( long iVib, long iRot );

	long int getLine( long iElecHi, long iVibHi, long iRotHi, long iElecLo, long iVibLo, long iRotLo, double *relint, double *absint );

	/* compute rate coefficient for a single quenching collision */
	realnum H2_CollidRateEvalOne( long iVibHi, long iRotHi, long iVibLo, long iRotLo, long ipHi, long ipLo, long nColl, double temp_K );

	void H2_Calc_Average_Rates( void );

	void H2_X_sink_and_source( void );
	/*H2_X_coll_rate_evaluate find collisional rates within X - 
	 * this is one time upon entry into H2_LevelPops */
	void H2_X_coll_rate_evaluate( void );

	/*H2_Level_low_matrix evaluate lower populations within X */
	/* total abundance within matrix */
	void H2_Level_low_matrix(realnum abundance );

	/*H2_Read_Cosmicray_distribution read distribution function for H2 population following cosmic ray collisional excitation
	void H2_Read_Cosmicray_distribution(void); */

	// read energies for all electronic levels
	void H2_ReadEnergies();
	void H2_ReadEnergies( long int nelec, vector<int>& n, vector<int>& v, vector<int>&J, vector<double>& eWN );

	/** read dissociation probabilities and kinetic energies for all electronic levels
	\param nelec
	*/
	void H2_ReadDissprob( long int nelec );

	/** H2_CollidRateEvalAll - set H2 collision rates */
	void H2_CollidRateEvalAll( void );

	/** read collision rates
	\param nColl
	*/
	void H2_CollidRateRead( long int nColl );

	/** read transition probabilities
	\param nelec
	*/
	void H2_ReadTransprob( long int nelec, TransitionList &trans );

	/**H2_Read_hminus_distribution read distribution function for H2 population following formation from H minus */
	void H2_Read_hminus_distribution(void);

	/**mole_H2_form find state specific rates grains and H- form H2 */
	void mole_H2_form( void );

	/**mole_H2_LTE sets Boltzmann factors and LTE unit population of large H2 molecular */
	void mole_H2_LTE( void );

	/**H2_Solomon_rate find rates between H2s and H2g and other levels,
	 * for use in the chemistry */
	void H2_Solomon_rate( void );

	/** gs_rate evaluate rate between ground and star states */
	double gs_rate( void );

	/** H2_zero_pops_too_low - zero out some H2 variables if we decide not to compute
	 * the full sim, called by H2_LevelPops*/
	void H2_zero_pops_too_low( void );

	/** create H2 molecules, called by ContCreatePointers */
	void init(void);

	/** set the ipCont struc element for the H2 molecule, called by ContCreatePointers */
	void H2_ContPoint( void );

	/**H2_DR choose next zone thickness based on H2 big molecule */
	double H2_DR(void);

	/** radiative acceleration due to H2 called in rt_line_driving */
	double H2_Accel(void);

	/**H2_RT_OTS - add H2 ots fields */
	void H2_RT_OTS( void );

	/** rad pre due to h2 lines called in PresTotCurrent*/
	double H2_RadPress(void);

	/** add in explicit lines from the large H2 molecule, called by lines_molecules */
	void H2_LinesAdd(void);

	/**H2_Reset called to reset variables that are needed after an iteration */
	void H2_Reset( void );

	/** internal energy of H2 called in PresTotCurrent */
	double H2_InterEnergy(void);

	/**H2_Colden maintain H2 column densities within X
	\param *chLabel
	*/
	void H2_Colden( const char *chLabel );

	/**H2_cooling evaluate cooling and heating due to H2 molecule */
	void H2_Cooling(void);

	/** LTE_Cooling_per_H2 compute exact cooling at LTE, currently for transitions
	 *  within the X band only. */
	double LTE_Cooling_per_H2();

	/** save H2 line data
	\param ioPUN io unit for save
	\param lgDoAll save all levels if true, only subset if false
	*/
	void H2_Punch_line_data(
		FILE* ioPUN ,
		bool lgDoAll );

	/** include H2 lines in punched optical depths, etc, called from SaveLineStuff
	\param io
	\param xLimit
	\param index
	*/
	void H2_PunchLineStuff( FILE * io , realnum xLimit  , long index);

	/** do emission from H2 - called from RT_diffuse */
	void H2_RT_diffuse(void);

	/** do RT for H2 lines
	*/
	void H2_RTMake( linefunc line_one );

	/** increment optical depth for the H2 molecule, called from RT_tau_inc */
	void H2_RT_tau_inc(void);

	/**H2_Prt_Zone print H2 info into zone results, called from prtzone for each printed zone */
	void H2_Prt_Zone(void);

	// print departure coefficients for all X levels
	void H2_PrtDepartCoef(void);

	/** initialize optical depths in H2, called from RT_tau_init */
	void H2_LineZero( void );

	/** the large H2 molecule, called from RT_tau_reset */
	void H2_RT_tau_reset( void );

	/** do level populations for H2, called by iso_solve */
	void H2_LevelPops( bool &lgPopsConverged, double &old_value, double &new_value );

	/** save some properties of the large H2 molecule
	\param io
	\param chJOB[]
	\param chTime[]
	\param ipPun
	*/
	void H2_PunchDo( FILE* io , char chJOB[] , const char chTime[] , long int ipPun );

	/**H2_ParseSave parse the save h2 command */
	void H2_ParseSave( Parser &p, ostringstream& chHeader );

	/**H2_itrzn - average number of H2 pop evaluations per zone */
	double H2_itrzn( void );

	/**H2_Prt_column_density print H2 info into zone results, called from prtzone for each printed zone
	\param *ioMEAN this is stream used for io, is stdout when called by final,
		is save unit when save output generated
	*/
	void H2_Prt_column_density( FILE *ioMEAN );

	void set_numLevelsMatrix( long numLevels );

	void H2_ReadDissocEnergies( void );

	/** flag saying whether molecular data have been read in yet */
	bool lgREAD_DATA;

	double photoionize_rate;
	double photo_heat_soft;
	double photo_heat_hard;
	double photodissoc_BigH2_H2s;
	double photodissoc_BigH2_H2g;

	/* total spontaneous dissociation rate [s-1],
	 * summed over excited electronic states, weighted by pops */ 
	double spon_diss_tot;

	double Solomon_dissoc_rate_g;
	double Solomon_dissoc_rate_s;
	
	/** these are decay rates from electronic levels into g and s */
	double Solomon_elec_decay_g;
	double Solomon_elec_decay_s;

	/** rate H2 goes from all X into either J=1 (ortho) 
	 * or (J=0) para on grain surfaces - units s-1*/
	double rate_grain_op_conserve;
	double rate_grain_J1_to_J0;

	/** H2 continuum photodissociation rate coefficient (not scaled by density) from P.C. Stancil data */
	double Cont_Dissoc_Rate_H2s;
	double Cont_Dissoc_Rate_H2g;
	multi_arr<double,3> Cont_Dissoc_Rate;
	
	/** LTE pops of g and s used for H- back reactions */
	double rel_pop_LTE_g;
	double rel_pop_LTE_s;

	/** average energy level of H2g and H2s */
	double average_energy_g;
	double average_energy_s;

	double HeatDiss;
	double HeatDexc;
	double HeatDexc_old;
	double HeatDexc_deriv;
	double HeatChangeOld, HeatChange;

	/** Average Einstein A for H2s to H2g transition*/
	double Average_A;
	/** Average noreactive collisional rate for H2s to H2g transition*/
	double Average_collH2_deexcit;
	double Average_collH_deexcit;
	double Average_collH2_excit;
	double Average_collH_excit;
	/**Average collisional dissociation of H2g and H2s by H and H2 */
	double Average_collH_dissoc_g;
	double Average_collH_dissoc_s;
	double Average_collH2_dissoc_g;
	double Average_collH2_dissoc_s;

	/** says whether model has ever been evaluated in this run - if it has
	 * not been then use TH85 physics for mole balance and cooling */
	bool lgEvaluated;

	/** index for threshold for photoionization */
	long ip_photo_opac_thresh;
	long ip_photo_opac_offset;

	t_coll_source coll_source[N_X_COLLIDER];

	/** the density (cm-3) of ortho H2 */
	double ortho_density,
	/** the density (cm-3) of para H2 */
		para_density;

	// single precision versions of the above	
	realnum ortho_density_f,
		para_density_f;

	/** column density in ortho and para H2 */
	double ortho_colden ,
		para_colden;

	/* old and older ortho - para ratios, used to determine whether soln is converged */
	double ortho_para_old, ortho_para_older, ortho_para_current;

	// these remember the largest and smallest factors needed to
	// renormalize the H2 chemistry 
	double renorm_max ,
		renorm_min;

	// this will say how many times the large H2 molecule has been called in this zone -
	// if not called (due to low H2 abundance) then not need to update its line arrays 
	long int nCall_this_zone;

	// flag saying whether to bother with the large molecule at all,
	// default is false, set true with atom h2 on command 
	bool lgEnabled;

	// this is the number of electronic levels to include in the output - default is 1,
	// only X.  changed with PRINT LINES H2 ELECTRONIC and option on PUNCH H2 LINES commands 
	int nElecLevelOutput;

	/** this is option to use estimates of the collision rates from g-bar approximations */
	/** turn mole.lgColl_gbar on/off with atom h2 gbar on off */
	bool lgColl_gbar;

	/** this is option to turn off the calculated collision rates */
	bool lgColl_deexec_Calc;

	/** this is option to turn off guesses of collisional dissociation rates */
	bool lgColl_dissoc_coll;

	// include collision rates that come from real calculations,
	// off with atom h2 collisions off command 
	bool lgH2_grain_deexcitation;

	/** flag to force LTE level populations, atom H2 LTE */
	bool lgLTE;

	/** option to turn off ortho-para collisions, command SPECIES H2 COLLISIONS ORTHO PARA OFF */
	bool lgH2_ortho_para_coll_on;

	// which set of He - H2 collisions to use? default is ORNL, other
	// is Le BOURlet 
	bool lgH2_He_ORNL;

	// flag saying whether (true) or not to use ORNL H2 - H2 collisions
	bool lgH2_ORH2_ORNL;
	bool lgH2_PAH2_ORNL;

	/** put noise into collision rates */
	bool lgH2_NOISE;
	/** noise for the CR collisions */
	bool lgH2_NOISECOSMIC;
			
	long int loop_h2_oscil;
	long int nzoneEval;
	
	/** std and mean for the noise, log normal distribution */
	double xMeanNoise , xSTDNoise;

	/** limit to the ratio H2/Htot - if ratio is below this, large atom is not called */
	double H2_to_H_limit;

	realnum mass_amu;

	/** turn on trace information */
	int nTRACE;
	
	/** this sets how fine a trace we want for atom  h2 trace */
	int n_trace_final , 
		n_trace_iterations , 
		n_trace_full,
		n_trace_matrix;

	// the number of electronic quantum states to include.
	// To do both Lyman and Werner bands want nelec = 3 
	long int n_elec_states;
		
	/* this is fraction of population that is within levels done with matrix */
	double frac_matrix;

	/* used to recall the temperature used for last set of Boltzmann factors */
	double TeUsedBoltz;
	double TeUsedColl;
	
	explicit diatomics( const string& a, const double& e_star, const double* const abund, double (*fun)(double) ) ;

	molecule *sp;	
	molecule *sp_star;
	qList states;	
	TransitionList trans;	
	TransitionList::iterator rad_end;
	vector< diss_tran > Diss_Trans;

private:
	string label;
	string shortlabel;
	string path;

	/** this is the energy (in cm-1), above which levels are considered to be H2*, 
	 * and below which they are H2 */
	/* >> chng 05 jul 15, TE, H2g = sum (v=0, J=0,1) */
	/* >>chng 05 jul 29, to 0.5 eV, this goes up to J=8 for v=0 */
	/* >>chng 05 aug 03, slight upward change in energy to include the J=8 level,
	 * also give energy in waveumbers for simplicity (save h2 levels give energy in ryd) */
	/*#define	ENERGY_H2_STAR	(0.5/EVRYD/WAVNRYD)*/
	/* energy of v=0, J=8 is 4051.73, J=9 is 5001.97
	 * v=1, J=0 is 4161.14 */
public:
	const double ENERGY_H2_STAR;

private:
	// pointer to the density of the species	
	const double* const dense_total;	
	
	char chH2ColliderLabels[N_X_COLLIDER][chN_X_COLLIDER];

	/* these vars are private for H2 but uses same style as all other header files -
	 * the extern is extern in all except cddefines */

	/** number of levels in H2g */
	long int nEner_H2_ground;

	/** total population in each vib state */
	multi_arr<double,2> pops_per_vib;

	/** the renorm factor for this H2 to the chemistry - should be unity */
	double H2_renorm_chemistry;

	/** rate [s-1] for collisions from ihi to ilo */
	multi_arr<realnum,2> H2_X_coll_rate;

	//int H2_nRot_add_ortho_para[N_ELEC];
	double H2_DissocEnergies[N_ELEC];
	/** number of vib states within electronic states */
	long int nVib_hi[N_ELEC];
	/** number of rotation levels within each elec - vib */
	valarray<long> nRot_hi[N_ELEC];
	/** this gives the first rotational state for each electronic state - J=0 does
	 * not exist when Lambda = 1 */
	long int Jlowest[N_ELEC];
	/** the number of ro-vib levels in each elec state */
	long int nLevels_per_elec[N_ELEC];
	/** the total population in each elec state */
	double pops_per_elec[N_ELEC];
	multi_arr<realnum,3> CollRateCoeff;
	multi_arr<realnum,3> CollRateErrFac; 
	vector<CollRateCoeffArray> RateCoefTable;

	// quantities dealing with chemistry 
#if 1
#endif

	// these are quantities for each state
#if 1
	/** these will mostly become xxx[elec][vib][rot] */
	multi_arr<realnum,3> H2_dissprob;
	multi_arr<realnum,3> H2_disske;
	multi_arr<double,3> H2_rad_rate_out;

	/** these will mostly become xxx[elec][vib][rot] */
	multi_arr<double,3> H2_old_populations;
	multi_arr<double,3> H2_populations_LTE;
	/** this is true if state is para, false if ortho */
	multi_arr<bool,3> H2_lgOrtho;
#endif

	long int nzoneAsEval , iterationAsEval;

	// these are quantities for states with X
#if 1
	multi_arr<int,2> H2_ipPhoto;
	multi_arr<double,2> H2_col_rate_in;
	multi_arr<double,2> H2_col_rate_out;
	multi_arr<double,2> H2_rad_rate_in;
	/** formation into specific states within X only vib and rot,
	 * includes both H- and H2 routes */
	multi_arr<realnum,2> H2_X_formation;
	/** backwards destruction of v,J levels due to the H- route */
	multi_arr<realnum,2> H2_X_Hmin_back;
	/** column density within X only vib and rot */
	multi_arr<realnum,2> H2_X_colden;
	/** LTE column density within X only vib and rot */
	multi_arr<realnum,2> H2_X_colden_LTE;
	/** rates [cm-3 s-1] from elec excited states into X only vib and rot */
	multi_arr<double,2> H2_X_rate_from_elec_excited;
	/** rates [s-1] to elec excited states from X only vib and rot */
	multi_arr<double,2> H2_X_rate_to_elec_excited;
	/** save rate coef (cm3 s-1) for collisional dissociation */
	multi_arr<realnum,2> H2_coll_dissoc_rate_coef;

	/** save rate coef (cm3 s-1) for collisional dissociation with H2g and H2s*/
	multi_arr<realnum,2> H2_coll_dissoc_rate_coef_H2;
#endif

	valarray<realnum> H2_X_source;
	valarray<realnum> H2_X_sink;

	/** distribution function for formation on grain surfaces,
	 * vib, rot, last dim is grain type */
	multi_arr<realnum,3> H2_X_grain_formation_distribution;

	/** density of H2s and H2g during current iteration */
	double H2_den_s , H2_den_g;

	/** vib, rot, last dim is temperature */
	multi_arr<realnum,3> H2_X_hminus_formation_distribution;

	valarray<long> ipVib_H2_energy_sort;
	valarray<long> ipElec_H2_energy_sort;
	valarray<long> ipRot_H2_energy_sort;
	multi_arr<long int,3> ipEnergySort;
	multi_arr<long int,2> ipTransitionSort;

	/** number of levels within X which are done with matrix solver,
	 * set with atom h2 matrix command */
	long int nXLevelsMatrix;
	long int ndim_allocated;
	multi_arr<double,2> AulEscp,
		AulDest, 
		AulPump,
		CollRate_levn;
	vector<double> pops, create, destroy, depart, stat_levn, excit;
	
	long int levelAsEval;
	bool lgFirst;
	long int nzone_eval;
	long int iteration_evaluated;

	/** this is array of accumulated line intensities, used for save he lines command */
	multi_arr<realnum,6> H2_SaveLine;

	/** fully defined array saying whether (true) or not (false) a radiative decay
	 * is defined by the standard emission line structure */
	multi_arr<bool,2> lgH2_radiative;

	/** counters used by H2_itrzn to find number of calls of h2 per zone */
	long int nH2_pops;
	long int nH2_zone;

	/** this is used to establish zone number for evaluation of number of levels in matrix */
	long int nzone_nlevel_set;

	/** the number of times the H2 molecules has been called in this iteration.  For the
	 * very first call we will use lte for the level populations, for later calls
	 * use the last solution */
	long int nCall_this_iteration;

	/* LTE cooling per molecule */
	vector<double> LTE_Temp, LTE_cool;

	/** print Matrix input to solver */
	bool lgPrtMatrix;
	/** create image for Matrix input to solver */
	bool lgImgMatrix;

public:
	/* Read LTE cooling per molecule */
	void H2_Read_LTE_cooling_per_H2();

	/* interpolate over LTE cooling array */
	double interpolate_LTE_Cooling( double Temp );
};

/* compute H2 continuum dissociation cross sections */
double MolDissocCrossSection( const diss_tran& tran, const double& Mol_Ene );

double Yan_H2_CS( double energy_ryd /* photon energy in ryd */);

/* compute H2 continuum dissoication opacities */
//double MolDissocOpacity( const diss_tran& tran, const double& Mol_Ene );

#endif /* H2_PRIV_H_ */

