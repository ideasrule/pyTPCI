/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef TRACE_H_
#define TRACE_H_

/* trace.h */

struct t_trace {

	/** flag saying that trace has been set */
	bool lgTrace;

	/**nznbug zone to turn on trace, set by trace command, npsbug is iteratoin number */
	long int nznbug;

	/** which iteration to turn on trace */
	long int npsbug;

	/**flag set to true by trace heavy command */
	bool lgHeavyBug;

	/**flag set to true by trace esource convergence command, will identify sources of electrons */
	bool lgESOURCE;

	/**trace convergence, number indicates level of trace needed */
	int nTrConvg;

	/** flag to indicate level of detail in trace, now only in trace convergence */
	/** flag set with trace leveln command, for n level atom */
	bool lgTrLevN;

	/** flag set with trace pointers command */
	bool lgPointBug;

	/** trace compton flag */
	bool lgComBug;

	/** flag set with trace neon command */
	bool lgNeonBug;

	/** flag set with trace line command */
	bool lgTrLine; 

	/** flag set with trace iron bug */
	bool lgFeBug;

	/** flag set if negative opacities every occured */
	bool lgOptcBug;

	/** trace 3 body recombination routines, trace three body */
	bool lgTrace3Bod;

	/** flag set with trace molecules command */
	bool lgTraceMole;

	/** flag set by trace heating command */
	bool lgHeatBug;

	/** flag set with trace dr command */
	bool lgDrBug;

	/** flag set with trace optimizer command */
	bool lgTrOptm;

	/** flag set with trace difuse fields command */
	bool lgTrDiff;

	/** flag set with trace beta command */
	bool lgTr8446;

	/** flag set with trace opacity */
	bool lgOpacBug;

	/** flag set with trace grains command */
	bool lgDustBug;

	/** flags set with trace helium (lgHeBug) */
	bool lgHeBug;

	/** lgHBug set with trace hydrogen command */
	bool lgHBug;

	/** set full trace with trace h-like or he-like full command */
	bool lgIsoTraceFull[NISO];

	/** ipIsoTrace is atomic number for iso-electronic species with full trace */
	long int ipIsoTrace[NISO];

	/** trace calcium atom flag */
	bool lgCalBug;

	/** trace carbon flag */
	bool lgCarBug;

	/** flag set with trace continuum command */
	bool lgConBug;

	/** flag set with trace ots command */
	bool lgOTSBug;

	/** flag set with trace two photon command */
	bool lgBug2nu;

	/** flag set with trace wind command */
	bool lgWind;

	/** flag set with trace cooling */
	bool lgCoolTr;

	/** set true if trace eden is entered */
	bool lgNeBug;

	/** debug level for use with dbg_printf command (in servicce.c)*/
	int debug_level;

	/** flag affecting which iteration to turn on trace */
	bool lgTrOvrd;

	/** trace secondary ionizaiton */
	bool lgSecIon;

	};

extern t_trace trace;

#endif /* TRACE_H_ */
