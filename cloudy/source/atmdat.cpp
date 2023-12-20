/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "atmdat.h"
#include "thirdparty.h"

t_atmdat atmdat;

void t_atmdat::zero()
{
	DEBUG_ENTRY( "t_atmdat::zero()" );
	nsbig = 0;

	/***************************************************
	 * charge transfer ionization and recombination 
	 ***************************************************/
	/* HCharHeatMax, HCharCoolMax are largest fractions of heating in cur zone
	 * or cooling due to ct */
	HCharHeatMax = 0.;
	HCharCoolMax = 0.;

	for ( long nelem=0; nelem < t_atmdat::NCX; ++nelem)
	{
		CharExcIonTotal[nelem] = 0.;
		CharExcRecTotal[nelem] = 0.;
	}
	HIonFrac = 0.;
	HIonFracMax = 0.;
	/* option to turn off all charge transfer, turned off with no charge transfer command */
	lgCTOn = true;

	/* flag saying that charge transfer heating should be included,
	 * turned off with no CTHeat commmand */
	HCharHeatOn = 1.;
	for( long nelem1=0; nelem1 < t_atmdat::NCX; ++nelem1)
	{
		for( long nelem=0; nelem< LIMELM; ++nelem )
		{
			for( long ion=0; ion<LIMELM; ++ion )
			{
				CharExcIonOf[nelem1][nelem][ion] = 0.;
				CharExcRecTo[nelem1][nelem][ion] = 0.;
			}
		}
	}

	/* >>chng 97 jan 6, from 0 to 8.5e-10*q as per Alex Dalgarno e-mail
	 * >>chng 97 feb 6, from 8.5e-10*q 1.92e-9x as per Alex Dalgarno e-mail */
	HCTAlex = 1.92e-9;

	lgInnerShellLine_on = true;
	lgInnerShell_Kisielius = true;
	lgUTAprint = false;
}

/** The threshold for Aul. For species where the highest level has no
	 transitions with Aul > aulThreshold, the highest levels will be
	 trimmed until a real Aul is found.*/
const double t_atmdat::aulThreshold = 1e-29;

multi_arr<double,4> HS_He1_Xsectn;
multi_arr<double,4> HS_He1_Energy;
multi_arr<vector<double>,4> OP_Helike_Xsectn;
multi_arr<vector<double>,4> OP_Helike_Energy;
multi_arr<long,4> OP_Helike_NumPts;

void ReadCollisionRateTable( CollRateCoeffArray& coll_rate_table, FILE* io, FunctPtr GetIndices, long nMolLevs, long nTemps, long nTrans )
{
	DEBUG_ENTRY( "ReadCollisionRateTable()" );

	// negative nTrans and/or nTemps will be signal that the numbers are not already known
	// They should never be zero.
	ASSERT( nTemps != 0 && nTrans != 0 );

	// Get the row of temperatures 
	string chLine;
	while( read_whole_line( chLine, io ) )
	{
		if( chLine[0] == '!' || chLine[0] == '#' || chLine[0] == '\n' || chLine[0] == '\r' )
			continue;
		else
			break;
	}
	ASSERT( chLine.length() > 0 );

	// fill the temperature array
	istringstream iss( chLine );
	while( true )
	{
		double help;
		iss >> help;
		if( !iss.good() )
			break;
		coll_rate_table.temps.push_back( help );
	}

	if( nTemps < 0 )
		nTemps = coll_rate_table.temps.size();
	else
		ASSERT( (unsigned)nTemps == coll_rate_table.temps.size() );

	ASSERT( coll_rate_table.collrates.size() == 0 );
	coll_rate_table.collrates.reserve( nMolLevs );
	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		coll_rate_table.collrates.reserve( ipHi, nMolLevs );
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			coll_rate_table.collrates.reserve( ipHi, ipLo, nTemps );
		}
	}
	coll_rate_table.collrates.alloc();
	// initialize to zero
	coll_rate_table.collrates.zero();

	long ipTrans = 0;
	while( read_whole_line( chLine, io ) )
	{
		if( chLine[0] == '!' || chLine[0] == '#' || chLine[0] == '\n' || chLine[0] == '\r' )
			continue;

		long i = 1;
		long ipHi = -1, ipLo = -1;
		(*GetIndices)( ipHi, ipLo, chLine.c_str(), i );
		ipTrans++;

		// sentinel that transition will not be used for whatever reason
		if( ipHi == -1 && ipLo == -1 )
			continue;

		// skip these lines
		if( ipLo >= nMolLevs || ipHi >= nMolLevs )
		{
			if( nTrans > 0 && ipTrans == nTrans )
				break;
			else
				continue;
		}

		/* Indices between the very highest levels seem to be reversed */
		if( ipHi < ipLo )
		{
			ASSERT( ipLo == nMolLevs - 1);
			long temp = ipHi;
			ipHi = ipLo;
			ipLo = temp;
		}

		ASSERT( ipHi >= 0 );
		ASSERT( ipLo >= 0 );

		bool lgEOL = false;
		for( long j = 0; j < nTemps; ++j )
		{
			coll_rate_table.collrates[ipHi][ipLo][j] = 
				(double)FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
			ASSERT( !lgEOL );
		}
	
		// try to read one more and assert that it gets lgEOL	
		FFmtRead( chLine.c_str(), &i, chLine.length(), &lgEOL );
		ASSERT( lgEOL );

		{
			enum {DEBUG_LOC=false};
			if( DEBUG_LOC )
			{
				printf("The values of up and lo are %ld & %ld \n",ipHi,ipLo);
				printf("The collision rates are");
				for( long i = 0; i < nTemps; ++i )
				{
					printf( "\n %e", coll_rate_table.collrates[ipHi][ipLo][i]);
				}
				printf("\n");
			}
		}

		if( nTrans > 0 && ipTrans == nTrans )
			break;
	}
	ASSERT( ipTrans > 0 );
	if( nTrans > 0 )
		ASSERT( ipTrans == nTrans );

	return;
}

double InterpCollRate( const CollRateCoeffArray& rate_table, const long& ipHi, const long& ipLo, const double& ftemp)
{
	DEBUG_ENTRY( "InterpCollRate()" );
	double ret_collrate = 0.;

	if( rate_table.temps.size() == 0 )
	{
		return 0.;
	}

	if( ftemp < rate_table.temps[0] )
	{
		ret_collrate = rate_table.collrates[ipHi][ipLo][0];
	}
	else if( ftemp > rate_table.temps.back() )
	{
		ret_collrate = rate_table.collrates[ipHi][ipLo][ rate_table.temps.size()-1 ];
	}
	else if( rate_table.temps.size() == 1 )
	{
		// lamda data files can have only one temperature (see cn.dat as of Feb. 10, 2009)
		ret_collrate = rate_table.collrates[ipHi][ipLo][0];
	}
	else
	{
		ret_collrate = linint(&rate_table.temps[0],
			&rate_table.collrates[ipHi][ipLo][0],
			rate_table.temps.size(),
			ftemp);
	}

	ASSERT( !isnan( ret_collrate ) );
	return(ret_collrate);
}

