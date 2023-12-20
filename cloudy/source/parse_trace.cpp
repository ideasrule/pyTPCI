/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseTrace read in options off the trace command line */
#include "cddefines.h"
#include "iterations.h"
#include "geometry.h"
#include "trace.h"
#include "iso.h"
#include "parser.h"

void ParseTrace(Parser &p )
{
	DEBUG_ENTRY( "ParseTrace()" );

	/* turn on trace at a certain zone; .LE.0 or none for starting trace
	 * optional second number is iteration to start debugger */

	/* generate string that says debug turned on - this is caught by perl script 
	 * if "no print" occurs on line do not print it - there is one sim in the
	 * test suite which tests debug print and do not want to trigger comment
	 * that debug prints are accidentally turned on */
	if( !p.nMatch("NO PR") )
		fprintf(ioQQQ,"DEBUG trace output turned on.\n");

	/* set initially false, in case we do not turn on trace until
	 * later iteration or zone */
	trace.lgTrace = false;
	/* this is which zone to turn on */
	trace.nznbug = (long)p.FFmtRead();
	if( p.lgEOL() )
		trace.lgTrace = true;

	/* this is which iteration to turn on */
	trace.npsbug = (long)p.FFmtRead();
	if( p.lgEOL() )
		trace.npsbug = 1;

	/* turn trace on now if no numbers on line */
	if( trace.nznbug == 0 && trace.npsbug <= 1 )
	{
		trace.lgTrace = true;
		geometry.nprint = 1;
		iterations.IterPrnt[0] = 1;
	}
	
	/* trace convergence is a special command, 
	 * only convergence loops, not full trace */
	if( p.nMatch("CONV") )
	{
		/* check for keyword, if not present
		 * then set to very high level of debugging - initially set to negative number, a sign
		 * that trace is not on yet, but to turn on trace convergence when we hit the right zone */
		/* 1 ConvPresTempEdenIoniz */
		if( p.nMatch("PRES") )
			trace.nTrConvg = -1;
		/* 2 ConvTempEdenIoniz*/
		else if( p.nMatch("TEMP") )
			trace.nTrConvg = -2;
		/* 3 ConvEdenIoniz*/
		else if( p.nMatch("EDEN") )
			trace.nTrConvg = -3;
		/* 4 ConvIoniz*/
		else if( p.nMatch("IONI") )
			trace.nTrConvg = -4;
		/* 5 below ConvBase*/
		/* > 5 all levels*/
		else
			trace.nTrConvg = -100;

		/* above set trace level to negative number - this will not trigger
		 * trace output - turn trace on now if no zone or iteration on line */
		if( trace.nznbug == 0 && trace.npsbug <= 1 )
			trace.nTrConvg *= -1;

		/* turn off normal trace parameters, this is a special case */
		trace.lgTrace = false;
		/*trace.nznbug = 10000;*/
		geometry.nprint = 10000;
		iterations.IterPrnt[0] = 10000;

		/* this is an option to also turn on ots rate debug prints */
		if( p.nMatch(" OTS") )
			trace.lgOTSBug = true;

		/* this is an option to also turn on electron density source debug prints,
		 * key is ESOURCE */
		if( p.nMatch("ESOU") )
			trace.lgESOURCE = true;
	}

	/* trace he-like and h-like must come early since they may have name of element */
	/* the trace h-like hydrogenic species command, with lots of options */
	if( p.nMatch("H-LI") )
	{
		/* turn on short trace for h-like species */
		trace.lgHBug = true;

		/* option to turn on full printout */
		if( p.nMatch("FULL") )
		{
			trace.lgIsoTraceFull[ipH_LIKE] = true;
		}
		else
		{
			trace.lgIsoTraceFull[ipH_LIKE] = false;
		}

		/* look for one of the element names on the line*/
		trace.ipIsoTrace[ipH_LIKE] = p.GetElem();

		/* if no element appears on the line GetElem fcn returns -1,
		 * in this case we want to do hydrogen */
		trace.ipIsoTrace[ipH_LIKE] = MAX2(0, trace.ipIsoTrace[ipH_LIKE] );
	}

	/* the trace h-like hydrogenic species command, with lots of options */
	if( p.nMatch("HE-L") )
	{
		/* turn on short trace for helium - like species */
		trace.lgHeBug = true;

		/* option to turn on full printout */
		if( p.nMatch("FULL") )
			trace.lgIsoTraceFull[ipHE_LIKE] = true;
		else
			trace.lgIsoTraceFull[ipHE_LIKE] = false;

		/* look for one of the element names on the line*/
		trace.ipIsoTrace[ipHE_LIKE] = p.GetElem();

		/* if no element appears on the line fcn returns -1,
		 * in this case we want to do helium */
		trace.ipIsoTrace[ipHE_LIKE] = MAX2(1, trace.ipIsoTrace[ipHE_LIKE] );
	}

	/* were there any keywords on the line? */
	if( p.nMatch("BETA") )
		trace.lgTr8446 = true;

	if( p.nMatch("CARB") )
		trace.lgCarBug = true;

	if( p.nMatch("CALC") )
		trace.lgCalBug = true;

	if( p.nMatch("COMP") )
		trace.lgComBug = true;

	if( p.nMatch("CONT") )
		trace.lgConBug = true;

	if( p.nMatch("COOL") )
		trace.lgCoolTr = true;

	if( p.nMatch("DIFF") )
		trace.lgTrDiff = true;

	if( p.nMatch(" DR ") )
		trace.lgDrBug = true;

	if( p.nMatch("EDEN") || p.nMatch("ELECTRON") )
		trace.lgNeBug = true;

	if( p.nMatch("GRAI") )
		trace.lgDustBug = true;

	if( p.nMatch("HEAV") )
		trace.lgHeavyBug = true;

	if( p.nMatch("HEAT") )
		trace.lgHeatBug = true;

	/* trace helium, but not h-like or he-like */
	if( p.nMatch("HELI") && !p.nMatch("H-LI")  && !p.nMatch("HE-L") )
		trace.lgHeBug = true;

	/* the simple trace hydrogen command */
	if( p.nMatch("HYDR") && !p.nMatch("H-LI"))
	{
		trace.lgHBug = true;
		trace.lgIsoTraceFull[ipH_LIKE] = false;
		/* this says which element, on the C scale (H=0), to trace */
		trace.ipIsoTrace[ipH_LIKE] = 0;
	}

	if( p.nMatch("IRON") )
		trace.lgFeBug = true;

	if( p.nMatch("LEVELN") )
		trace.lgTrLevN = true;

	if( p.nMatch("LINE") )
		trace.lgTrLine = true;

	if( p.nMatch("NEON") )
		trace.lgNeonBug = true;

	
	if( p.nMatch("HMOL") || p.nMatch("CMOL") )
	{
		fprintf( ioQQQ," This command is deprecated. Please use TRACE MOLE.  Sorry.\n" );
		cdEXIT( EXIT_FAILURE );
	}
	/* turn on molecular trace */
	if( p.nMatch(" MOLE") )
	{
		trace.lgTraceMole = true;
	}

	/* trace pointers */
	if( p.nMatch("POIN") )
		trace.lgPointBug = true;

	/* following two are optical, optimize */
	if( p.nMatch("OPTIC") )
		trace.lgOptcBug = true;

	if( p.nMatch("OPTIM") )
		trace.lgTrOptm = true;

	if( p.nMatch(" OTS") )
		trace.lgOTSBug = true;

	if( p.nMatch("SECO") && p.nMatch("IONI") )
		/* secondary ionization */
		trace.lgSecIon = true;

	if( p.nMatch("THRE") )
		trace.lgTrace3Bod = true;

	/* two photon emission, spontaneous and induced */
	if( p.nMatch(" TWO") )
		trace.lgBug2nu = true;

	/* wind geometry */
	if( p.nMatch("WIND") )
		trace.lgWind = true;

	/* falling through is fine - just turn on minimal trace */
	return;
}
