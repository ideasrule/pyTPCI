/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef GRAINVAR_H_
#define GRAINVAR_H_

/* grainvar.h */

#include "container_classes.h"

/** flag that determines if quantum heating is to be taken into account in H2
 * grain surface formation rate, set this true to enable quantum heating treatment, PvH */
const bool ENABLE_QUANTUM_HEATING=false;

/** number of grid points for which grain emissivity is evaluated in InitEmissivities */
const int NDEMS = 200;

/** these are used for setting up grain emissivities in InitEmissivities() */
const double GRAIN_TMIN = 1.e-3;
const double GRAIN_TMID = 5.e3;
const double GRAIN_TMAX = 1.2e9;

/** maximum number of discrete grain charge states that will be kept in memory for possible reuse */
const int NCHS = 50;
/** this is the largest number of charge bins that may be in use at any one time; the
 * remaining NCHS-NCHU charge bins are used for backing up data for possible later use */
const int NCHU = NCHS/5;
/** default number of charge states to be used */
const int NCHRG_DEFAULT= 2;

/** maximum size of temperature grid for quantum heating routine <BR>
 >>chng 02 aug 01, changed 10000 -> 20000, because of laqheat2.in test<BR> 
 >>chng 04 nov 09, changed back to 10000 to save CPU time, laqheat2.in still OK */
const int NQGRID = 10000;

/** all grains should conserve charge and energy to within this precision */
const double CONSERV_TOL = 1.e-3;

/** the grain depletion functions - possible values of gv.nDustFunc */
typedef enum {
	DF_STANDARD,      /**< standard, PAHs vary, but classical grains do not */
	DF_USER_FUNCTION, /**< user supplied function, can be anything... */
	DF_SUBLIMATION    /**< grain abundance is suppressed above sublimation temperature */
} df_type;

/** the following constants are used to define the enthalpy function */
typedef enum {
	ENTH_CAR,   /**< graphite, Guhathakurtha & Draine 1989, ApJ, 345, 230 */
	ENTH_CAR2,  /**< graphite, Draine & Li 2001, ApJ, 551, 807 */
	ENTH_SIL,   /**< silicate, Guhathakurtha & Draine 1989, ApJ, 345, 230 */
	ENTH_SIL2,  /**< silicate, Draine & Li 2001, ApJ, 551, 807 */
	ENTH_PAH,   /**< PAH molecules, Dwek et al. 1997, ApJ, 345, 230 */
	ENTH_PAH2,  /**< PAH molecules, Draine & Li 2001, ApJ, 551, 807 */
	ENTH_SIC    /**< SiC, Chekhovskoy V. Ya., 1971, J. Chem. Thermodynamics, 3, 289 */
} enth_type;

/** the following constants are used to define the expression for Z_min */
typedef enum {
	ZMIN_CAR,   /**< graphitic material, Weingartner & Draine */
	ZMIN_SIL    /**< silicate grains, grey grains, Weingartner & Draine */
} zmin_type;

/** the following constants are used to define the expression for the ionization potential */
typedef enum {
	POT_CAR,    /**< graphitic material, Weingartner & Draine */
	POT_SIL     /**< silicate grains, grey grains, Weingartner & Draine */
} pot_type;

/** the following constants are used to define the expression for the inverse attenuation length */
typedef enum {
	IAL_CAR,    /**< graphitic material, Weingartner & Draine */
	IAL_SIL     /**< silicate grains, grey grains, Weingartner & Draine */
} ial_type;

/** the following constants are used to define the expression for the photo-electric effect */
typedef enum {
	PE_CAR,     /**< graphitic material, Weingartner & Draine */
	PE_SIL      /**< silicate grains, grey grains, Weingartner & Draine */
} pe_type;

/** the following constants are used to define where the emitted spectrum is stored */
typedef enum {
	STRG_CAR,   /**< graphitic spectrum */
	STRG_SIL    /**< silicate spectrum */
} strg_type;

/** the following constants are used to define which H2 ro-vib distribution should be used at formation, see h2.c */
typedef enum {
	H2_ICE,     /**< ice mantles */
	H2_SIL,     /**< silicate grains */
	H2_CAR,     /**< graphitic material */

	/** H2_TOP is the number of different grain surfaces that are treated -
	 * used for H_2 processes that occur on grain surfaces */

	H2_TOP      /**< NB NB, this MUST always be the last definition !! */
} H2_type;

/** the following constants are used to define material types<BR>
 * NB NB - whenever you define a new material type, fill in the physical properties
 *         in GrainZero (i.e., make initializations for gv.which_enth[MAT_xxx], etc...) !! */
typedef enum {
	MAT_USR=0,  /**< reserved */
	MAT_CAR,    /**< graphitic material, Guhathakurtha & Draine enthalpy */
	MAT_SIL,    /**< silicate grains, grey grains, Guhathakurtha & Draine enthalpy */
	MAT_PAH,    /**< PAH molecules, Dwek et al. enthalpy */
	MAT_CAR2,   /**< graphitic material, Draine & Li enthalpy */
	MAT_SIL2,   /**< silicate grains, grey grains, Draine & Li enthalpy */
	MAT_PAH2,   /**< PAH molecules, Draine & Li enthalpy */
	MAT_SIC,    /**< SiC molecules, Chekhovskoy enthalpy */

	/** MAT_TOP is the number of different material types that are treated */

	MAT_TOP     /**< NB NB, this MUST always be the last definition !! */
} mat_type;

struct GrainPar {
	double dep;
	df_type nDustFunc;
	bool lgForbidQHeating,
	  lgGreyGrain,
	  lgRequestQHeating;
};

class ShellData
{
	void p_clear0();
	void p_clear1();

public:
	ShellData()
	{
		p_clear1();
	}
	void clear()
	{
		p_clear0();
		p_clear1();
	}

	long nelem;             /**< element number, -1 = band, 0 = hydrogen, etc. */
	long ns;                /**< shell number, -1 = band, 0 = K shell, etc. */
	double ionPot;          /**< ionization potential for this shell, work function for band */
	long ipLo;              /**< offset into p, y0, y1, and y01 arrays */
	flex_arr<realnum> p;    /**< probability that photon absorption occurs in this shell */
	flex_arr<realnum> y01;  /**< product y0*y1 for primary electron */
	long nData;             /**< number of Auger electrons with different energies */
	vector<double> AvNr;    /**< AvNr[nData[ns]]: no. of electrons per primary ionization */
	vector<double> Ener;    /**< Ener[nData[ns]]: energy of electron in Ryd */
	vector< flex_arr<realnum> > y01A;/**< y0*y1 for Auger electrons */
};

/** this struct stores the data for the energy spectrum
 * of emitted Auger electrons for a single element
 * The data are currently taken from Table 4.1 and 4.2 of
 * >>refer	grain	Dwek, E., & Smith, R.K., 1996, ApJ, 459, 686 */
class AEInfo
{
	void p_clear0();
	void p_clear1();

public:
	AEInfo()
	{
		p_clear1();
	}
	void clear()
	{
		p_clear0();
		p_clear1();
	}

	unsigned int nSubShell;  /**< number of subshells of this element for which data is present */
	vector<unsigned int> nData;/**< nData[nSubShell]: number of data points for each subshell */
	vector<double> IonThres; /**< IonThres[nSubShell]: ionization threshold for electron in this subshell in Ryd */
	vector< vector<double> > AvNumber;
	                         /**< AvNumber[nSubShell][nData[ns]]: no. of electrons per primary ionization */
	vector< vector<double> > Energy;/**< Energy[nSubShell][nData[ns]]: energy of electron in Ryd */
};

/** NB NB NB NB NB NB
 *
 * this is the data structure for all grain data that depends on the charge state
 * (i.e. all data that used to have an [NCHS] dependance in days of old),
 *
 * each data item will be referenced as:   gv.bin[nd]->chrg[nz]->data_item
 *
 * this structure is allocated for each charge state at run time.
 *
 * Data items that do not depend on charge state, or are summed over all
 * charge states, should go in gv or gv.bin[nd] below !!
 *
 * NB NB !! every field that is in ChargeBin should be initialized or invalidated in UpdatePot1/2 !! NB NB */

class ChargeBin
{
	void p_clear0();
	void p_clear1();

public:
	ChargeBin()
	{
		p_clear1();
	}
	void clear()
	{
		p_clear0();
		p_clear1();
	}

	/** grain charging */
	long DustZ,             /**< grain charge, in e */
	  ipThresInf,           /**< pointer to ThresInf in anu array */
	  ipThresInfVal,        /**< pointer to ThresInfVal in anu array */
	  nfill;                /**< remember how far the flex_arr's were filled in */
	double FracPop,         /**< fractional population for lower and upper charge state */
	  Emin,                 /**< negative threshold to be overcome by outgoing electron, in Ryd */
	  EminInc,              /**< same as Emin, but for incoming electrons (using Zg instead of Zg+1), in Ryd */
	  PotSurf,              /**< potential difference between grain surface and infinity (phi_g), in Ryd */
	  PotSurfInc,           /**< same as PotSurf, but for incoming electrons (using Zg instead of Zg+1), in Ryd */
	  ThresInf,             /**< threshold for removing electron from grain well to infinity (phi), in Ryd */
	  ThresInfInc,          /**< same as ThresInf, but for incoming electrons (using Zg instead of Zg+1), in Ryd */
	  ThresInfVal,          /**< threshold for removing electron from valence band to infinity (phi^v), in Ryd */
	  ThresSurf,            /**< threshold for removing electron from grain well to surface (phi_s), in Ryd */
	  ThresSurfInc,         /**< same as ThresSurf, but for incoming electrons (using Zg instead of Zg+1), in Ryd */
	  ThresSurfVal,         /**< threshold for removing electron from valence band to surface (phi_s^v), in Ryd */
	  ThermRate;            /**< thermionic rate, in e/cm^2/s */
	flex_arr<realnum> yhat; /**< electron yield per absorbed photon (incl. Auger and secondary electrons) */
	flex_arr<realnum> yhat_primary;/**< electron yield per absorbed photon (only primary electrons) */
	flex_arr<realnum> ehat; /**< average energy per photo-electron (incl. Auger and secondary electrons), Ryd */
	flex_arr<double> cs_pdt;/**< photo-detachment cross section, for default depl */

	long RecomZ0[LIMELM][LIMELM+1];/**< ionization stage the atom/ion recombines to upon impact */

	double eta[LIMELM+2],   /**< cache for GrainScreen results */
	  xi[LIMELM+2];         /**< cache for GrainScreen results */

	double RSum1;           /**< cache for electron gain from colliding electrons */
	double RSum2;           /**< cache for electron gain from colliding ions */
	double ESum1a;          /**< cache for PE rate from band and inner shells (incl. Auger elec. & secondaries) */
	double ESum1b;          /**< cache for PE rate from conduction band */
	double ESum2;           /**< cache for electron loss due to recombination with colliding ions */

	/** grain heating */
	realnum tedust;         /**< equilibrium temperature for this zone */
	double hcon1;           /**< reminder of heating integral over incident flux, Ryd/H/s at default depl */
	double hots1;           /**< reminder of heating integral over diffuse fields, Ryd/H/s at default depl */
	double bolflux1;        /**< reminder of heating integral over all fields, Ryd/H/s at default depl */
	double pe1;             /**< reminder of photoelectric heating integral, Ryd/H/s at default depl */
	flex_arr<double> fac1;  /**< auxiliary data for GrainTemperature */
	flex_arr<double> fac2;  /**< auxiliary data for GrainTemperature */

	realnum RecomEn[LIMELM][LIMELM+1];/**< chemical energy released into grain upon impact, in Ryd */
	realnum ChemEn[LIMELM][LIMELM+1];/**< net contribution of ion recomb to grain heating, in Ryd */

	/** heating/cooling balance, all entries are valid for current zone, actual depl, and are in erg/cm^3/s */
	double BolFlux,         /**< total photon flux absorbed, used for energy conservation test */
	  GrainHeat,            /**< total heating of current grain type */
	  GrainHeatColl,        /**< collisional heating of current grain type */
	  GasHeatPhotoEl,       /**< photoelectric heating of the gas, added in GrGH 0 */
	  GasHeatTherm,         /**< heating due to thermionic emission */
	  GrainCoolTherm,       /**< grain cooling due to thermionic emissions, summed over charge states */
	  ChemEnIon,            /**< net amount of energy donated by recombining ions */
	  ChemEnH2;             /**< net amount of energy donated by H2 formation on grain surface */

	/** quantum heating */
	double HeatingRate2;    /**< quantum heating by electron recomb - thermionic cooling, erg/H/s, default depl */
};

/** NB NB NB NB NB NB
 *
 * this is the data structure for all grain data that depends on grain type
 * (i.e. all data that can differ from one grain bin to the next),
 *
 * each data item will be referenced as:   gv.bin[nd]->data_item
 *
 * this structure is allocated for each grain bin at run time.
 *
 * Data items that are generic for all grain types, or are summed over all
 * grain types, should go in gv below !!
 *
 * NB NB NB NB NB NB */

class GrainBin
{
	void p_clear0();
	void p_clear1();

	vector<ChargeBin> pool; /**< pool for storing charge dependent data, includes currently unused bins */
public:
	GrainBin()
	{
		p_clear1();
	}
	void clear()
	{
		p_clear0();
		p_clear1();
	}

	/** grain logic */
	df_type nDustFunc;      /**< has user requested custom grain abundance function ? */
	bool lgPAHsInIonizedRegion;/**< were PAHs present in the ionized region ? */
	bool lgIterStart;       /**< are we converging the grain physics for the first time this iteration? */

	/** general information on the grains */
	char chDstLab[13];      /**< label for the species */
	double eec;             /**< pow(dustp[0],-0.85), needed for electron esacpe length */
	double eyc;             /**< 1./AvRadius + 1.e7, needed for electron yield */
	realnum dustp[5],       /**< 0 = specific weight (g/cm^3), 1 = mol. weight (amu), 2 = default abundance,
	                         *   3 = default depletion, 4 = fraction of the mass in this grain bin */
	  AvRadius,             /**< average grain radius, <a^3>/<a^2>, in cm */
	  AvArea,               /**< average grain surface area, <4pi*a^2>, in cm^2, CURRENTLY NOT USED */
	  AvVol,                /**< average grain volume, <4/3pi*a^3>, in cm^3 */
	  IntRadius,            /**< integrated grain radius Int(a), normalized per H, in cm/H */
	  IntArea,              /**< integrated grain surface area Int(4pi*a^2), normalized per H, in cm^2/H */
	  IntVol,               /**< integrated grain volume Int(4/3pi*a^3), normalized per H, in cm^3/H */
	  elmAbund[LIMELM],     /**< chemical composition, abundance at default depl, see comment below */ 
	  atomWeight,           /**< molecular weight per atom, in amu */
	  Tsublimat,            /**< sublimation temperature */
	  DustWorkFcn,          /**< work function, in Ryd */
	  BandGap,              /**< gap between valence and conduction band, in Ryd */
	  ThermEff,             /**< efficiency of thermionic emission, between 0 and 1 */
	  avDGRatio;            /**< Integral(D/G*dReff) for average dust to gas ratio, OUTPUT ONLY */
	mat_type matType;       /**< material type, 1 = carbonaceous, 2 = silicate, 3 = PAH, etc... */

	/** chemical abundances of grains in elmAbund are defined as follows:
	 * e.g. for MgSiFeO4 the oxygen abundance would be given by
	 *    elmAbund[nd][7] = 4.*dustp[nd][2]*dustp[nd][3]*dustp[nd][4] */
	/** chemical abundances of grains are normalized such that in any given zone
	 * the total fractional abundance of an element Z locked up in grains is given by
	 *    SumFracAbund(nelem) = Sum_over_nd ( gv.bin[nd]->elmAbund[nelem]*gv.dstAbund[nd] ) */

	/** grain depletion, normalization conversion factors */
	realnum dstfactor,      /**< grain depletion factor, dep from GRAINS command */
	  dstAbund,             /**< grain abundance in zone, dstfactor*GrainMetal*GrnVryDpth(radius) */
	  GrnDpth;              /**< grain abundance scale factor in zone, GrnStdDpth(radius),
				 *   used by set PAH constant / H0 commands */
	double cnv_H_pGR,       /**< grain unit conversion, \<unit\>/H (default depl) -> \<unit\>/grain */
	  cnv_H_pCM3,           /**< grain unit conversion, \<unit\>/H (default depl) -> \<unit\>/cm^3 (actual depl) */
	  cnv_CM3_pGR,          /**< grain unit conversion, \<unit\>/cm^3 (actual depl) -> \<unit\>/grain */
	  cnv_CM3_pH,           /**< grain unit conversion, \<unit\>/cm^3 (actual depl) -> \<unit\>/H (default depl) */
	  cnv_GR_pH,            /**< grain unit conversion, \<unit\>/grain -> \<unit\>/H (default depl) */
	  cnv_GR_pCM3;          /**< grain unit conversion, \<unit\>/grain -> \<unit\>/cm^3 (actual depl) */

	/** grain opacities */
	double RSFCheck;        /**< save resolution scale factor for later check */

	/** >>chng 02 dec 30, separated scattering cross section and asymmetry factor (1-g),
	 * NB NB NB -- note that pure_sc1 DOES NOT contain the asymmetry factor, while gv.dstsc DOES !!! */
	vector<double> dstab1;  /**< absorption cross section per grain species, for default depl */
	vector<double> pure_sc1;/**< scattering cross section per grain species, for default depl */
	vector<double> asym;    /**< asymmetry factor (1-g) */
	vector<double> dstab1_x_anu; /**< helper array, dstab1[i]*anu[i] for default depl */

	/** equilibrium temperature */
	double dstems[NDEMS],   /**< grain emissivity at dsttmp[], default depl, normalized per H */
	  dstslp[NDEMS],        /**< auxiliary array for spline interpolation */
	  dstslp2[NDEMS];       /**< auxiliary array for inverse spline interpolation */

	bool lgTdustConverged;  /**< is dust temperature converged ? */
	realnum tedust,         /**< equilibrium temperature for this zone */
	  TeGrainMax,           /**< highest equilibrium temperature as a function of radius */
	  avdust;               /**< Integral(Tdust*dReff) for average equilibrium temperature, OUTPUT ONLY */

	/** grain charging, photoelectric effect, thermionic emissions
	 *
	 * all charge and energy rates will be calculated by resolving the charge distribution into nChrg
	 * integral charge states. To implement this, certain parameters have been moved into the ChargeBin
	 * structure, currently limiting the maximum number of charge states the code can handle. For details see:
	 * >>refer	grain	physics	van Hoof et al., 2001, ASP Conf. Series 247, p. 353 (astroph/0107183) */
	long LowestZg;          /**< lowest charge a grain can ever have, in e */
	long nfill;             /**< remember how far the flex_arr's in the ShellData were filled in */
	vector<ShellData> sd;   /**< specific data for each atomic shell in this grain material */
	vector<realnum> y0b06;  /**< bulk yield for band according to Eq. 9 of WDB06 */
	double AveDustZ;        /**< average charge per grain, in electrons */
	long ZloSave;           /**< grain charge at the start of the iteration, in e */
	double Capacity;        /**< grain capacity, in Farad/grain */
	double dstpot,          /**< grain potential in Ryd */
	  RateUp,               /**< total charging rate up, used for balance check, in e/cm^2/s */
	  RateDn,               /**< total charging rate down, used for balance check, in e/cm^2/s */
	  StickElecNeg,         /**< sticking efficiency for electrons on negative or neutral grains */
	  StickElecPos;         /**< sticking efficiency for electrons on positive grains */
	realnum avdpot,         /**< Integral(Vg*dReff) for average grain potential, OUTPUT ONLY */
	  le_thres;             /**< threshold for using X-ray prescription for l_e, in Ryd */
	vector<realnum> inv_att_len;/**< inverse attenuation length (in cm) */
	double AccomCoef[LIMELM];/**< accommodation coefficient, needed for collisional heating of grain */

	/** heating/cooling balance, all entries are valid for current zone, actual depl, and are in erg/cm^3/s */
	double BolFlux,         /**< total photon flux absorbed, used for energy conservation test */
	  GrainCoolTherm,       /**< grain cooling due to thermionic emissions, summed over charge states */
	  GasHeatPhotoEl,       /**< photoelectric heating of the gas, added in GrGH 0 */
	  GrainHeat,            /**< total heating of current grain type */
	  GrainHeatColl,        /**< collisional heating of current grain type */
	  GrainGasCool,         /**< gas cooling due to collisions with grains */
	  ChemEn,               /**< net amount of energy donated by recombining ions */
	  ChemEnH2,             /**< net amount of energy donated by H2 formation on grain surface */
	  thermionic;           /**< heating due to thermionic emission */

	/** quantum heating physics */
	bool lgQHeat,           /**< is quantum heating turned on ? */
	  lgUseQHeat,           /**< should quantum heating be used for this zone ? */
	  lgEverQHeat,          /**< was quantum heating used in any zone ? */
	  lgQHTooWide;          /**< is probability distribution too wide to fit in NQGRID array elements ? */
	long QHeatFailures,     /**< counter for number of times qheat algorithm failed */
	  qnflux,               /**< like rfield.nflux, but may point to higher energy, for phiTilde and Phi */
	  qnflux2;              /**< like rfield.nflux, only for max electron energy, for phiTilde and Phi */
	double qtmin;           /**< lowest grain temperature used in calculations, set per zone */
	double qtmin_zone1;     /**< lowest grain temperature used in calculations, initial zone */
	double HeatingRate1;    /**< quantum heating due to molecule/ion collisions, erg/H/s, default depl */
	double DustEnth[NDEMS], /**< grain enthalpy at dsttmp[], in Ryd/grain */
	  EnthSlp[NDEMS],       /**< auxiliary array for spline interpolation */
	  EnthSlp2[NDEMS];      /**< auxiliary array for inverse spline interpolation */

	/** H2 physics - each has units s^-1 */
	double rate_h2_form_grains_HM79;/**< H2 formation rate, Hollenbach & McKee 79, units s^-1, actual depl */
	double rate_h2_form_grains_CT02;/**< H2 formation rate, Cazaux & Tielens 02, units s^-1, actual depl */
	double rate_h2_form_grains_ELRD;/**< H2 formation rate, Rollig et al. 2013 with Eley-Rideal effect,
				 * units s^-1, actual depl */
	double rate_h2_form_grains_used;/**< H2 rate actually used, evaluated in hmole.c, units s^-1, actual depl
                                 * when multiplied with hden, this is formation rate in H2-molecules/cm^3/s */

	/** grain drift */
	realnum DustDftVel,     /**< grain drift velocity for this zone */
	  avdft;                /**< Integral(vdrift*dReff) for average drift velocity, OUTPUT ONLY */

	long nChrgOrg;          /**< number of charge states at the start of the iteration */
	long nChrg;             /**< number of charge states used for the current zone */
	long ichrg[NCHS];       /**< list of indices into pool[], the first nChrg values are in use */
	ChargeBin& chrg(long nz) { return pool[ichrg[nz]]; }
	const ChargeBin& chrg(long nz) const { return pool[ichrg[nz]]; }
};

/** NB NB NB NB NB NB
 *
 * this is the master data structure for grain physics, it is statically allocated
 *
 * all entries that depend on grain type should go in the GrainBin structure above
 *
 * NB NB NB NB NB NB */

class GrainVar
{
	void p_clear0();
	void p_clear1();

public:
	GrainVar()
	{
		p_clear1();
	}
	void clear()
	{
		p_clear0();
		p_clear1();
	}

	/** grain logic */
	bool lgDustOn() const          /**< have any grains been switched on ? */
	{
		return ( bin.size() > 0 );
	}
	bool lgWD01,                   /**< use physics as described in Weingartner & Draine (2001) */
	  lgReevaluate,                /**< reevaluate keyword on grains cmd => always reevalutae grain 
	                                *   conditions during calculation */
	  lgGrainPhysicsOn,            /**< this is option to turn off all grain physics while leaving
	                                *   the opacity in, set false with no grain physics command */
	  lgAnyDustVary,               /**< do any of the grain abundances vary with depth */
	  lgBakesPAH_heat;             /**< turn on simple formula for PAH heating of gas */
	bool lgNegGrnDrg;              /**< flag set if negative grin drag force encountered */

	/** test logic */
	bool lgDHetOn,                 /**< default true, turned off with GRAIN NO HEATING */
	  lgQHeatOn;                   /**< default true, turned off with GRAIN NO QHEAT */

	/** turn off grain gas collisional energy exchange, set false with 
	 * NO GRAIN COLLISIONAL ENERGY EXCHANGE */
	bool lgDColOn;                 /**< default true, turned off with GRAIN NO COOLING */

	/** should electrons from/to grains be included in the total electron sum? 
	 * del true, set false with no grain electrons command */
	bool lgGrainElectrons;

	long nCalledGrainDrive;        /**< count how many times GrainDrive has been called */

	vector<string> ReadRecord;     /**< record of all the files read by mie_read_opc */

	string chPAH_abundance;        /**< which functions describes the PAH abundances, changed with SET PAH */

	/** grain depletion */

	/** chemical abundances of grains are normalized such that in any given zone
	 * the total number density of an element nelem locked up in grains is given by
	 *    elmSumAbund(nelem) = Sum_over_nd ( elmAbund[nd][nelem]*dstAbund[nd]*dense.hden ) */

	realnum GrainMetal,            /**< grain depletion factor, from METALS xxx GRAINS command, usually 1 */
	  elmSumAbund[LIMELM];         /**< number density of elements summed over all grain bins, at actual depl */ 

	/** these arrays define the physical properties of material MAT_PAH, MAT_CAR, MAT_SIL, etc... */
	/** they are initialized in the constructor */
	enth_type which_enth[MAT_TOP]; /**< defines expression for the enthalpy function */
	zmin_type which_zmin[MAT_TOP]; /**< defines expression for Z_min */
	pot_type which_pot[MAT_TOP];   /**< defines expression for the ionization potential */
	ial_type which_ial[MAT_TOP];   /**< defines expression for the inverse attenuation length */
	pe_type which_pe[MAT_TOP];     /**< defines expression for the photo-electric effect */
	strg_type which_strg[MAT_TOP]; /**< defines where the emitted spectrum is stored */
	H2_type which_H2distr[MAT_TOP];/**< defines expression for H2 ro-vib distribution at formation */

	/** grain opacities */
	long nzone;                    /**< remember in what zone grain quantities were last updated */
	vector<double> dstab;          /**< total absorption cross section, current depl is factored in */
	vector<double> dstsc;          /**< total scattering cross section, current depl and asymmetry factored in */

	/** grain charging */
	double TotalEden;              /**< contribution to eden from all grain species, a positive number means
	                                *   that the grains contribute to the free electron pool, in cm^-3 */
	realnum GrnElecDonateMax;      /**< largest local fraction of electrons donated by grains */
	realnum GrnElecHoldMax;        /**< largest local fraction of electrons contained by grains */
	double GrnRecomTe;             /**< the electron temperature that was used to calculate the recom rates */
	long nChrgRequested;           /**< number of charge states requested by user, default value is 2 */

	AEInfo AugerData[LIMELM];      /**< store data about energy spectrum of Auger electrons */

	/** grain surface recombination rate for element nelem, ionization stage "ion-from" to "ion-to"
	 * already multiplied by grain area, actual depletion, GrainChTrRate[nelem][ion-from][ion-to], units s^-1,
	 * this can be both ionization and neutralization process, so ion-to can be > or < ion-from */
	realnum GrainChTrRate[LIMELM][LIMELM+1][LIMELM+1];
	
	/** heating/cooling balance, all entries are valid for current zone, actual depl, and are in erg/cm^3/s */
	double GasCoolColl,            /**< cooling of the gas by collisions with all grains, summed in GrGC 0 */
	  GasHeatPhotoEl,              /**< heating of the gas by photoelectric effect of all grains */
	  GasHeatTherm,                /**< heating of the gas by thermionic emissions of all grains */
	  GasHeatNet,                  /**< net heating of the gas: PhotoEl + Therm - CoolColl */
	  GrainHeatSum,                /**< total heating of all grain types, added in GraT 0 */
	  GrainHeatLya,                /**< heating of all grains by Lya, added in GraL 1216 */
	  GrainHeatDif,                /**< heating of all grains by all diffuse fields (incl Lya), added in GraD 0 */
	  GrainHeatInc,                /**< heating of all grains by incident continuum, added in GraI 0 */
	  GrainHeatCollSum,            /**< collisional heating of all grains, added in GraC 0 */
	  GrainHeatChem;               /**< net amount of energy donated by recombining ions to all grains */

	double dHeatdT;                /**< gas heating derivative for all grains, in erg/cm^3/s/K */
	
	realnum GrainHeatScaleFactor;  /**< scale factor for PE heating as per Allers et al. 2005 (SET GRAINS HEAT) */

	realnum TotalDustHeat,         /**< total PE heating integrated over model, erg/s or erg/cm^2/s, actual depl */
	  dphmax,                      /**< largest fraction of local heating due to grain PE heating */
	  dclmax;                      /**< largest local cooling of gas by collisions with grains */

	/** equilibrium temperature */
	double dsttmp[NDEMS];          /**< grain temperature grid for dstems (in GrainBin) */

	/** quantum heating physics */
	realnum dstAbundThresholdNear; /**< min grain abundance for enabling quantum heating near I-front */
	realnum dstAbundThresholdFar;  /**< min grain abundance for enabling quantum heating in molecular regions */
	bool lgQHeatAll;               /**< use quantum heating for all grains, WILL BECOME DEFAULT */
	bool lgQHPunLast;              /**< only save quantum heating information on last iteration */
	FILE *QHSaveFile;              /**< file pointer for PUNCH QHEAT command */

	/** H2 physics */
	double rate_h2_form_grains_used_total; /**< rate H2 forms on grains, summed over bins, units s^-1, actual depl
	                                * when multiplied with hden, this is formation rate in H2-molecules/cm^3/s */

	/** grain emission */
	vector<realnum> GrainEmission; /**< total emission from this zone, per unit vol, for add outward in metfic */
	vector<realnum> GraphiteEmission;/**< graphite emission from this zone only */
	vector<realnum> SilicateEmission;/**< silicate emission from this zone only */

	/** per bin grain data */
	vector<GrainBin> bin;          /**< pointers to memory allocated for bins */
};

extern GrainVar gv;

#endif /* GRAINVAR_H_ */
