/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*InitDefaultsPreparse initialization at start of simulation, called from cloudy
* before parser, sets initial values of quantities changed by parser 
* called for each point in a grid but one time in multi-iteration sims */
#include "cddefines.h" 
#include "init.h"
#include "phycon.h"
#include "radius.h"
#include "trace.h"
#include "dynamics.h"
#include "geometry.h"
#include "noexec.h"
#include "opacity.h"
#include "stopcalc.h"
#include "rt.h"
#include "fudgec.h"
#include "abund.h"
#include "ionbal.h"
#include "hextra.h"
#include "wind.h"
#include "atmdat.h"
#include "pressure.h"
#include "thermal.h"
#include "continuum.h"
#include "save.h"
#include "hcmap.h"
#include "prt.h"
#include "grid.h"
#include "iso.h"
#include "rfield.h"
#include "flux.h"
#include "dense.h"
#include "species.h"
#include "grainvar.h"

/*InitDefaultsPreparse initialization at start of simulation, called from cloudy
 * before parser, will now be called one time per sim in grid but long term
 * goal is to call parser only one time in a grid.  This would be called 
 * first to set defaults that may be changed by the parser.  This routine
 * should only set defaults for variables that change in the parser.
 * It does not initialize variables for the start of a calculation.  That
 * is done by InitSimPostparse 
 * called one time in multi iteration sim */
void InitDefaultsPreparse( void )
{
	long int i, 
		ipISO,
		nelem;

	DEBUG_ENTRY( "InitDefaultsPreparse()" );

	/* init vars before parsing commands for each sim */

	/* option to turn off collisional ionization with "no collisional ionization" cmmnd */
	atmdat.lgCollIonOn = true;

	// default, use Chianti, use Chianti-level2 hybrid, do not print
	atmdat.lgChiantiOn = true;
	atmdat.lgChiantiLvl2Hybrid = true;
	atmdat.lgChiantiPrint = false;
	//Use gbar to fill in dBase transitions if they lack collision strengths
	atmdat.lgGbarOn = true;
	//Tells Cloudy to exclusively use experimental energies in Chianti.
	atmdat.lgChiantiExp = true;
	// Set the default number of Chianti energy levels to use for Fe for photoionization case
	atmdat.nChiantiMaxLevelsFe = atmdat.nDefaultPhotoLevelsFe;
	// Set the default number of Chianti energy levels to use for all other elements
	atmdat.nChiantiMaxLevels = atmdat.nDefaultPhotoLevels;
	// ChiantiLevelsSet is false until the user specifies the number of chianti levels to use.
	atmdat.lgChiantiLevelsSet = false;
	/* nChiantiMaxLevelsFe and nChiantiMaxLevels are defaulted to 100 and 50 respectively if
	 * the coronal command is used. See parse_coronal.cpp */

	// Set the default number of Stout energy levels to use for Fe for photoionization case
	atmdat.nStoutMaxLevelsFe = atmdat.nDefaultPhotoLevelsFe;
	// Set the default number of Stout energy levels to use for all other elements
	atmdat.nStoutMaxLevels = atmdat.nDefaultPhotoLevels;
	// StoutLevelsSet is false until the user specifies the number of stout levels to use.
	atmdat.lgStoutLevelsSet = false;
	/* nStoutMaxLevelsFe and nStoutMaxLevels are defaulted to 100 and 50 respectively if
	 * the coronal command is used. See parse_coronal.cpp */

	/* Stout on by default, use level2 OP hybrid,  optional printout off */
	atmdat.lgStoutOn = true;
	atmdat.lgStoutPrint = false;
	atmdat.lgStoutLvl2Hybrid = true;

	/* Lamda on by default, optional printout off */
	atmdat.lgLamdaOn = true;
	atmdat.lgLamdaPrint = false;
	atmdat.nLamdaMaxLevels = atmdat.nDefaultMolLevels;
	// LamdaLevelsSet is false until the user specifies the number of Lamda levels to use.
	atmdat.lgLamdaLevelsSet = false;

	// By default use Dima (Voronov97) data for collisional ionization rate coefficients
	// possible options are DIMA or HYBRID
	atmdat.CIRCData = t_atmdat::HYBRID;

	/** Set the default collision strength for dBase transitions when there is no radiative data */
	atmdat.collstrDefault = 1e-10;

#	ifdef USE_CDMS
	atmdat.lgCalpgmOn = true;
#	else
	atmdat.lgCalpgmOn = false;
#	endif
	strcpy(atmdat.chCloudyChiantiFile, "CloudyChianti.ini");
	strcpy(atmdat.chStoutFile, "Stout.ini");
	strcpy(atmdat.chLamdaFile, "Lamda.ini");

	/* drChange was reset to get orion flux in h-beta correct
	 * drChange is really tau of current zone */
	radius.drChange = 0.15f;

	radius.glbdst = 0.;
	radius.glbrad = 0.;

	// energy conservation check is slow, turn on with set check energy every zone
	continuum.lgCheckEnergyEveryZone = false;

	/* option to read in all input quantities and NOT execute the actual model
	* only check on input parameters - set by calling cdNoExec */
	noexec.lgNoExec = false;

	/* constant for the extinguish command */
	rfield.ExtinguishColumnDensity = 0.;
	rfield.ExtinguishLeakage = 0.;
	rfield.ExtinguishLowEnergyLimit = 1.;

	/* parameters having to do with thermal map */
	hcmap.RangeMap[0] = 10.f;
	hcmap.RangeMap[1] = .99e10f;
	/* zone where map is to be done */
	hcmap.MapZone = -1;
	hcmap.nMapStep = 20;

	thermal.ConstGrainTemp = 0.;
	thermal.lgTemperatureConstant = false;
	thermal.lgTemperatureConstantCommandParsed = false;
	thermal.ConstTemp = 0.;
	thermal.lgTeHigh = false;
	thermal.lgTeBD96 = false;
	thermal.lgTeTLaw = false;
	thermal.lgTLaw = false;
	thermal.lgTeSN99 = false;
	thermal.tlaw.clear();
	/* try toe predict next zone's temperature in constant density models,
	* as done in ZoneStart.  Turned off with no tepred command */
	thermal.lgPredNextTe = true;

	/* turbulent heating - set with hextra command */
	hextra.TurbHeat = 0.;
	/* set true with TIME option on hextra command for time dependent sims */
	hextra.lgTurbHeatVaryTime = false;
	/* options for for extra heating to depend on scale radius */
	hextra.lgHextraDepth = false;
	/* options for for extra heating to depend on density */
	hextra.lgHextraDensity = false;
	/* options for alpha disk model heating */
	hextra.lgHextraSS = false;

	/* options set by cosmic ray command */
	hextra.cryden = 0.;
	hextra.cryden_ov_background = 0.;
	hextra.lg_CR_B_equipartition = false;
	hextra.crtemp = 0.;
	hextra.crpowr = 0.;
	hextra.cr_energydensity = 0;

	/* set with coronal equilibrium init time command */
	dynamics.lg_coronal_time_init = false;
	dynamics.lgTracePrint = false;

	/* parameters to do with wind */
	wind.lgWindOK = true;
	wind.DiskRadius = 0;
	wind.lgDisk = false;
	wind.windv0 = 0.;
	wind.setStatic();
	wind.comass = 0.;
	wind.windv = 0.;
	wind.dvdr = 0.;
	wind.emdot = 0.;
	wind.AccelAver = 0.;
	wind.acldr = 0.;
	wind.AccelGravity = 0.;
	wind.AccelTotalOutward = 0.;
	wind.AccelCont = 0.;
	wind.AccelElectron = 0.;
	wind.AccelLine = 0.;
	wind.AccelMax = 0.;
	wind.fmul = 0.;
	wind.lgVelPos = true;

	/* argument on ELEMENT LIMIT OFF XXX command, lowest abundance 
	 * element to include in the calculation */
	dense.AbundanceLimit = 0.;

	/* controls density fluctuations, when true variations are due to density changing,
	* when false are due to abundances changing if changes are on */
	dense.lgDenFlucOn = true;
	dense.lgDenFlucRadius = true;
	dense.flong = 0.;
	dense.flcPhase = 0.;

	/* this says keep initial density constant, 
	* so pressure from iter to iter not really const */
	dense.lgDenseInitConstant = true;

	// pressure does not vary with time by default
	dense.lgPressureVaryTime = false;
	//  required number is timescale for time variation
	dense.PressureVaryTimeTimescale = -1.;
	//  optional number is index for time variation
	dense.PressureVaryTimeIndex = 0.;

	/* extra electron density, set with eden command */
	dense.EdenExtra = 0.;

	/* forced electron density, set with set eden command */
	dense.EdenSet = 0.;

	/* option to set electron fraction, n_e/n_H */
	dense.EdenFraction = 0.;

	/* individual terms for the pressure equation of state */
	/* >>chng 05 jan 01, all three are set true at start *
	* code default is constant density, so all this is ignored
	* is turned to constant pressure then these flags must be adjusted so
	* that we get the pressure we expect.  with all true, all three
	* contributors to pressure will be counted - with constant gas pressure
	* these are turned off */
	pressure.lgPres_radiation_ON = true;
	pressure.lgPres_magnetic_ON = true;
	pressure.lgPres_ram_ON = true;
	/* flag for constant pressure, include continuum too */
	pressure.lgContRadPresOn = true;
	/* constant density is the default */
	strcpy( dense.chDenseLaw, "CDEN" );
	dense.DLW.clear();
	/* number on line is log of nT - option to specify initial pressure */
	pressure.lgPressureInitialSpecified = false;
	/* this is log of nT product - if not present then set zero */
	pressure.PressureInitialSpecified = 0;

	/* select certain atomic data, changed with set atomic data command 
	* this says to use Zeippen 1982 [OII] transition probabilities */
	dense.lgAsChoose[ipOXYGEN][1] = false;

	abund.lgAbnReference = false;
	abund.lgAbundancesSet = false;

	/* option to turn off an element */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		/* set of scale factors for changing abundances with elements command */
		abund.ScaleElement[nelem] = 1.;

		// default scale factors for SET DIELECTRONIC RECOMBINATION KLUDGE SCALE
		ionbal.DR_mean_scale[nelem] = 1.;
	}

	abund.lgAbTaON = false;

	/* option to turn off an element */
	for( nelem=0; nelem < LIMELM; nelem++ )
	{
		/* option to have abundances from table */
		abund.lgAbunTabl[nelem] = false;
	}

	/* threshold for faintest heating cooling to save with save heating or 
	 * save cooling commands, reset with PUNCH WEAKHEATCOOL command */
	save.WeakHeatCool = 0.05f;

	/* set of variables used to control save command */
	save.nsave = 0;
	save.lgPunContinuum = false;
	save.lgDRPLst = false;
	save.lgDRHash = true;
	save.lgDTPLst = false;
	save.lgDTHash = true;
	save.lgTraceConvergeBaseHash = true;
	/* this is the string that will appear after each model in the save output,
	 * reset with the "set save hash" command */
	save.chHashString = "###########################";
	/* save every one continuum point, set skip save command */
	save.ncSaveSkip = 1;
	/* flush file after every iteration - set save flush */
	save.lgFLUSH = false;
	/* do not use old style, luminosity per unit inner cloud area */
	save.lgLuminosityOld = false;

	for( i=0; i<LIMPUN; ++i )
	{
		save.lgHashEndIter[i] = true;
		/* set false for time dependent calculations*/
		save.lg_separate_iterations[i] = true;
		save.lgSaveEveryZone[i] = false;
		save.nSaveEveryZone[i] = -1;
		save.emisfreq[i].set( -1. );
		save.ipEmisFreq[i] = -1;
	}

	/* default is to conserve energy, reset with
	 * set save line width / resolution command */
	save.Resolution = realnum(-1.);
	save.ResolutionAbs = realnum(-1.);

	/* punch dominant rates variables */
	for( long i=0; i<LIMPUN; i++ )
		save.chSpeciesDominantRates[i] = "" ;

	/* subtract continuum from line fluxes, etc */
	save.lgSubtrCont = true;

	// Switches for save Database command
	// Wavenumbers vs WLAng, Gf vs A, Collision Rates or not
	save.lgSaveDataWn = false;
	save.lgSaveDataGf = false;
	save.lgSaveDataRates = false;

	/* default no printing of optical depths, TooFaint is .1 */
	prt.lgPrtTau = false;
	prt.PrtTauFnt = 0.1f;
	prt.lgPrtShort = false;
	prt.TooFaint = 1e-3f;
	prt.lgFaintOn = true;
	prt.lgFntSet = false;
	prt.lgPrnLineCell = false;
	prt.nPrnLineCell = -1000;
	prt.lgPrtCitations = false;

	/* change angle of illumination 
	* this is angle away from the normal, so 0 is a normal ray, the default*/
	geometry.DirectionalCosin = 1.;
	geometry.size = 1.f;
	geometry.lgSizeSet = false;
	geometry.covaper = -1.f;

	/* if true then print main block of lines as array,
	 * set false with print lines column, will then
	 * do a single column of lines */
	prt.lgPrtLineArray = true;

	/* when printing a column this is option to print linear rather than log */
	prt.lgPrtLineLog = true;

	/* print ages */
	prt.lgPrnAges = false;

	/* print column densities */
	prt.lgPrintColumns = true;

	/* option to turn off printing of the main line block */
	prt.lgPrintBlock = true;
	prt.lgPrintBlockEmergent = true;
	prt.lgPrintBlockIntrinsic = true;

	/* option to sort lines by wavelength, print sort command */
	prt.lgSortLines = false;

	prt.lgPrtMaser = false;

	prt.lgPrtContIndices = false;
	prt.lgPrnPump = false;
	prt.lgPrnInwd = false;
	prt.lgPrnColl = false;
	prt.lgPrnHeat = false;

	/* >>chng 00 dec 08, these determine the standard items included in "nFnu", PvH */
	prt.lgSourceReflected = true;
	prt.lgSourceTransmitted = false;
	prt.lgDiffuseInward = true;
	prt.lgDiffuseOutward = true;

	prt.lgPrtLastIt = false;
	prt.lgOnlyZone = false;
	prt.lgOnlyHead = false;
	prt.lgPrtStart = false;
	prt.nstart = 0;
	/* print predictions from collapsed levels of iso sequences */
	prt.lgPrnIsoCollapsed = true;

	/* turn off printing of heating agents */
	prt.lgPrintHeating = false;

	/* flag saying to print all ionization balance elements
	 * set with PRINT ARRAYS IONIZATION command */
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		prt.lgPrtArry[nelem] = false;
	}

	/* print line flux at earth */
	prt.lgPrintFluxEarth = false;

	/* print line surface brightness, def sr, option arcsec */
	prt.lgSurfaceBrightness = false;
	prt.lgSurfaceBrightness_SR = true;

	/* print line cumulative sets true, print integrated line intensity over
	 * time in temp dependent simulation */
	prt.lgPrintLineCumulative = false;

	prt.lgPrintLineAirWavelengths = true;

	prt.lgPrintHTML = false;

	prt.nzdump = -100;

	prt.blend.clear();

	trace.lgSecIon = false;
	trace.lgTrOvrd = true;
	trace.lgOpacBug = false;
	trace.nTrConvg = 0;
	trace.lgTr8446 = false;
	trace.lgTrLevN = false;
	trace.lgOptcBug = false;
	trace.lgTrace3Bod = false;
	trace.lgOTSBug = false;
	trace.lgESOURCE = false;
	trace.lgTraceMole = false;
	trace.lgHeatBug = false;
	trace.lgHeavyBug = false;
	trace.lgBug2nu = false;
	trace.lgDrBug = false;
	trace.lgWind = false;
	trace.lgDustBug = false;
	trace.lgComBug = false;
	trace.lgHeBug = false;
	trace.lgCarBug = false;
	trace.lgCalBug = false;
	trace.lgConBug = false;
	trace.lgNeBug = false;
	trace.lgFeBug = false;
	trace.lgHBug = false;
	trace.lgTrLine = false;
	trace.nznbug = 10000;
	trace.npsbug = 10000;
	trace.lgTrace = false;
	trace.lgPointBug = false;
	trace.lgNeonBug = false;
	trace.lgCoolTr = false;
	trace.lgTrDiff = false;
	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		trace.lgIsoTraceFull[ipISO] = false;

	/* variables used in stop ... command */

	/* various criteria for stopping model */
	/* >>chng 04 dec 21, remove from here and init to 1e30 in zero */
	/*StopCalc.tauend = 0.;*/
	StopCalc.tauend = 1e30f;

	/* >>chng 05 nov 22 - NPA.  Stop calculation when fraction of oxygen frozen
	* out on grains ices gets too high - formation of ices */
	/*StopCalc.StopDepleteFrac = 0.99f;*/
	/* >>chng 05 dec 16, with revised ion solver logic, code should be able to
	* converge away from situation where ices have disturbed the chemistry and
	* net negative atomic abundances result.  now we say solution not converged and
	* soldier on 
	* this test should not be necessary */
	StopCalc.StopDepleteFrac = 1.02f;

	StopCalc.xMass = 0.;
	StopCalc.taunu = 0.;
	StopCalc.iptnu = -1;
	/* stopping AV */
	StopCalc.AV_extended = 1e30f;
	StopCalc.AV_point = 1e30f;
	/* highest allowed temperature */
	StopCalc.TempHiStopZone = (realnum)phycon.TEMP_LIMIT_HIGH;
	StopCalc.TempHiStopIteration = (realnum)phycon.TEMP_LIMIT_HIGH;

	/* the floor sets a limit to the temperature in the calculation -
	* if te falls below this, we do a constant temperature cloud at
	* this temperature */
	StopCalc.TeFloor = 0.;

	/* stop zone calculations when Te falls below this,
	 * TEMP_STOP_DEFAULT in phycon.h and is 4000 */
	StopCalc.TempLoStopZone = (realnum)phycon.TEMP_STOP_DEFAULT;
	/* stop iterations, used to stop time dependent command */
	StopCalc.TempLoStopIteration = -1.;

	/* ending column densities */
	StopCalc.HColStop = COLUMN_INIT;
	StopCalc.colpls = COLUMN_INIT;
	StopCalc.colnut = COLUMN_INIT;
	StopCalc.col_h2 = COLUMN_INIT;
	StopCalc.col_h2_nut = COLUMN_INIT;
	StopCalc.col_H0_ov_Tspin = COLUMN_INIT;
	StopCalc.col_monoxco = COLUMN_INIT;

	StopCalc.lgStopSpeciesColumn = false;
	StopCalc.chSpeciesColumn[0] = '\0';
	StopCalc.col_species = COLUMN_INIT;

	/* stopping electron density */
	StopCalc.StopElecDensity = -COLUMN_INIT;

	/* stopping electron and molecular fractions */
	StopCalc.StopElecFrac = -FLT_MAX;
	StopCalc.StopHPlusFrac = -FLT_MAX;
	/* stopping molecular fraction has opposite sign - want to stop when 2H_2/NH gt this */
	StopCalc.StopH2MoleFrac = FLT_MAX;
	/* this flag says that 21cm line optical depth is the stop quantity */
	StopCalc.lgStop21cm = false;
	/* stop when absolute value of velocity falls below this */
	StopCalc.StopVelocity = 0.;
	/* number of stop line commands entered */
	StopCalc.nstpl = 0;

	/* initialize some variables for the optimizer */
	optimize.nIterOptim = 400;
	optimize.OptGlobalErr = 0.10f;
	optimize.nEmergent = 0;
	optimize.lineids.clear();
	optimize.errorwave.clear();
	optimize.ipobs.clear();
	optimize.xLineInt_Obs.clear(); 
	optimize.xLineInt_error.clear();	
	optimize.chTempLab.clear();
	optimize.ionTemp.clear();
	optimize.temp_obs.clear(); 
	optimize.temp_error.clear();
	optimize.chTempWeight.clear();
	optimize.chColDen_label.clear();
	optimize.ion_ColDen.clear();
	optimize.ColDen_Obs.clear();
	optimize.ColDen_error.clear();
	optimize.ContIndex.clear();
	optimize.ContEner.clear();
	optimize.ContNFnu.clear();
	optimize.ContNFnuErr.clear();
	optimize.nRangeSet = 0;
	strcpy( optimize.chOptRtn, "PHYM" );

	/* flags says what is to be matched */
	optimize.lgOptLum = false;
	optimize.lgOptDiam = false;
	optimize.lgOptimize = false;
	optimize.lgInitialParse = false;

	/* trace flag for optimization process */
	optimize.lgTrOpt = false;

	optimize.lgOptimFlow = false;
	optimize.optint = 0.;
	optimize.optier = 0.;
#	if defined(__unix) || defined(__APPLE__)
	optimize.lgParallel = ( cpu.i().MPIMode() == MS_DEFAULT );
	grid.lgParallel = true;
#	else
	optimize.lgParallel = false;
	grid.lgParallel = false;
#	endif
	if( optimize.lgParallel )
		optimize.useCPU = cpu.i().nCPU();
	else
		optimize.useCPU = 1;
	optimize.lgOptCont = false;

	grid.lgNegativeIncrements = false;
	grid.lgSaveXspec = false;
	grid.lgKeepMainOutputSeparate = false;
	if( grid.lgParallel )
		grid.useCPU = cpu.i().nCPU();
	else
		grid.useCPU = 1;
	grid.nCycle = 1;
	grid.lgCrash = false;
	grid.lgCrashEval = false;

	/* the fudge factors command */
	fudgec.nfudge = 0;
	fudgec.lgFudgeUsed = false;
	for( i=0; i < NFUDGC; i++ )
		fudgec.fudgea[i] = 0.;

	EmLineZero( DummyEmis );
	TauZero( DummyEmis );
	DummyEmis.iRedisFun() = 0;
	DummyEmis.ipFine() = -1;
	DummyEmis.gf() = 0.;
	DummyEmis.damp() = 0.;
	DummyEmis.dampXvel() = 0.;
	DummyEmis.opacity() = 0.;
	DummyEmis.Aul() = 1e-30_r;

	/* following were block data logic */
	rt.lgStarkON = true;

	/* by default use Federman form of shielding function */
	rt.nLineContShield = LINE_CONT_SHIELD_FEDERMAN;

	/* parameters set with Case A and Case B commands */
	/* this is flag for turning on case b */
	opac.lgCaseB = false;

	/* this is separate flag for turning off collisions from n=2 */
	opac.lgCaseB_HummerStorey = false;

	/* this is separate flag for turning off excited state photoionization */
	opac.lgCaseB_no_photo = false;
	/* another case b option, turn off background opacities, no Pdest */
	opac.lgCaseB_no_pdest = false;

	/* smallest allowed line and Lya optical depths, reset with
	* Case B command */
	opac.tlamin = 0.f;

	/* taumin command minimum optical depths for lines default 1e-20 */
	opac.taumin = 0.f;

	opac.eeFreeFreeTemp = -1.;

	/* set false with no induced processes */
	rfield.lgInducProcess = true;

	/* this flag says that CMB has been set */
	rfield.lgCMB_set = false;
	
	rfield.lgComptonOn = true;

	for( i=0; i < LIMSPC; i++ )
	{
		/* this is set true if particular continuum source can vary with time 
		* set true if TIME appears on intensity / luminosity command line */
		rfield.lgTimeVary[i] = false;
		/* most continua enter as a beam rather than isotropic */
		rfield.lgBeamed[i] = true;
		// default is radiation from the "illuminated" face
		rfield.Illumination[i] = Illuminate::FORWARD;
		// optical depth = normal optical depth * this scale factor,
		// is 1 / cos theta
		rfield.OpticalDepthScaleFactor[i] = 1.;
		/* default energy range is H-ionizing radiation */
		rfield.range[i][0] = HIONPOT;
		rfield.range[i][1] = rfield.egamry();
		rfield.RSFCheck[i] = 0.;
		rfield.tNu[i].clear();
		rfield.tslop[i].clear();
		rfield.tFluxLog[i].clear();
		rfield.TableRadius[i] = -1.;
		rfield.lgSphericalDilution[i] = false;
		rfield.ncont[i] = 0;
	}

	/* line overlap opacity, turn off with no fine opacity command */
	rfield.lgOpacityFine = true;

	pseudoContDef.wlLo = 1000.;
	pseudoContDef.wlHi = 7000.;
	pseudoContDef.nBins = 1000;

	/* this is the faintest the high-energy tail of the continuum be */
	rfield.FluxFaint = 0.; // 1e-10;

	/* >>chng 01 jul 26, moved next statement from below loop to avoid bug in gcc 2.95.3, PvH */
	/* default diffuse fields is outward only */
	strcpy( rfield.chDffTrns, "OU2" );
	rfield.lgOutOnly = true;

	/* flags for whether continuum is defined over all energies */
	rfield.lgMMok = true;
	rfield.lgHPhtOK = true;
	rfield.lgXRayOK = true;
	rfield.lgGamrOK = true;

	/* set logical flags saying whether to include element in AGN tables */
	/* first set all false, since most not included */
	for( long int i=0; i < LIMELM; i++ )
	{
		abund.lgAGN[i] = false;
	}
	abund.lgAGN[ipHYDROGEN] = true;
	abund.lgAGN[ipHELIUM] = true;
	abund.lgAGN[ipCARBON] = true;
	abund.lgAGN[ipNITROGEN] = true;
	abund.lgAGN[ipOXYGEN] = true;
	abund.lgAGN[ipNEON] = true;
	abund.lgAGN[ipMAGNESIUM] = true;
	abund.lgAGN[ipSILICON] = true;
	abund.lgAGN[ipSULPHUR] = true;
	abund.lgAGN[ipARGON] = true;
	abund.lgAGN[ipIRON] = true;

	for( long int i=0; i < LIMELM; i++ )
	{
		abund.IsoAbn[i].init();
	}

	gv.clear();

	clear_lines_table();

	return;
}
