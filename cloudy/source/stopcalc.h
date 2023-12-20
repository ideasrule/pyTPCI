/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef STOPCALC_H_
#define STOPCALC_H_

class Flux;

/** maximum number of stop line commands */
const int MXSTPL = 10;

/** this is the stopping column density that is "default" - if still at
 * this value then stop column is not used */
const realnum COLUMN_INIT = 1e30f;

const int nCHREASONSTOP = 100;

/** stopcalc.h */
struct t_StopCalc {
	/** this has the various ending criteria for stopping a model */

	/** optical depth to stop calculation */
	realnum tauend;
	/** this is the energy within the continuum mesh, we will stop when the
	 * optical depth at this energy reaches tauend */
	realnum taunu;

	/** this points to taunu within the continuum mesh */
	long int iptnu;

	/** this provides a "floor" for the temperature - when the temperature
	 * falls to this limit, go to constant temperature solution */
	double TeFloor;

	/** highest allowed zone temperature, set with stop temperature exceeds command */
	realnum TempHiStopZone;
	/** highest allowed iteration temperature, set with stop time temperature exceeds command */
	realnum TempHiStopIteration;

	/** TempLoStopZone is lowest temperature to allow in radial zone integrations, 
	 * set with stop temperature command */
	realnum TempLoStopZone;
	/** TempLoStopIteration is lowest temperature to allow in iterations, 
	 * set with stop time temperature command, used to stop time dependent sims */
	realnum TempLoStopIteration;

	/** STOP EFRAC sets this limiting ratio of electron to H densities */
	realnum StopElecFrac;

	/** stop at a hydrogen molecular fraction, relative to total hydrogen,
	 * this is 2H_2 / H_total*/
	realnum StopH2MoleFrac;

	/** stop at a ionized hydrogen fraction, this was put in to simulate 
	 * the contribution of the H+ region to the overall intensity of PDR lines */
	realnum StopHPlusFrac;

	/** stop when a certain fraction of O frozen out on grains is reached, 
	 * this was put in to stop Cloudy from going into regimes where time 
	 * dependence must be considered */
	realnum StopDepleteFrac;

	/** stop when absolute value of velocity falls below this value.  
	 * Entered in km/s but converted to cm/s */
	realnum StopVelocity;

	/** limit column densities set with stop column command
	 * HColStop is total hydrogen column, others are ionized and neutral */
	realnum HColStop, 
	  colpls, 
	  colnut;

	/** column density in molecular hydrogen */
	realnum col_h2;

	/** stop at some mass */
	realnum xMass;

	/** column density in molecular + neutral hydrogen */
	realnum col_h2_nut;

	/** integrated n()H0) / Tspin */
	realnum col_H0_ov_Tspin;

	/** column density in carbon monoxide */
	realnum col_monoxco;

	/** stop at AV, point or extended */
	realnum AV_point , AV_extended;

	/** stopping electron density, set with stop eden command */
	realnum StopElecDensity;

	/** MXSTPL (above, equal to 10) is maximum number of stop line commands 
	 *
	 * parameters for stop line command
	 * LineStopWl - line wavelength for stop line command
	 * ipStopLin1 is pointer to first line of pair in stack of lines 
	 * ipStopLin2 is pointer to second line of pair in stack of lines, 
	 * usually Hbeta.  initially holds wavelength, converted to pointer in 
	 * iter_startend.cpp
	 */
	realnum stpint[MXSTPL];

	/** line indices for first and second lines in stop line command */
	long int ipStopLin1[MXSTPL], 
	  ipStopLin2[MXSTPL];
	/** number of stop line commands entered, 0 if none */
	long int nstpl;

	realnum StopLineWl1[MXSTPL],
		StopLineWl2[MXSTPL];

	/** line labels */
	char chStopLabel1[MXSTPL][NCHLAB] ,
		chStopLabel2[MXSTPL][NCHLAB];
	int nEmergent[MXSTPL];

	/** flag saying to stop at 21cm line optical depth */
	bool lgStop21cm;

	/** parameters for STOP CONTINUUM FLUX */
	vector<long> ContIndex;
	vector<Flux> ContNFnu;

	/** the reason the calculation stopped */
	char chReasonStop[nCHREASONSTOP];

	bool lgStopSpeciesColumn;
	string chSpeciesColumn;
	realnum col_species;

	};

extern t_StopCalc StopCalc;


#endif /* STOPCALC_H_ */
