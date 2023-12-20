/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "save.h"
#include "cddrive.h"
#include "grid.h"
#include "rfield.h"
#include "prt.h"
#include "input.h"
#include "version.h"
#include "service.h"
#include "container_classes.h"

static const int RECORDSIZE = 2880;
static const int LINESIZE = 80;

#if defined(_BIG_ENDIAN) 
	/* the value of A will not be manipulated */
#	define HtoNL(A) (A)	
/*
#	define HtoNS(A) (A)
#	define NtoHS(A) (A)
#	define NtoHL(A) (A)
*/
#else
/* defined(_LITTLE_ENDIAN) */
/* the value of A will be byte swapped */
#	define HtoNL(A) ((((A) & 0xff000000) >> 24) | \
		(((A) & 0x00ff0000) >> 8) | \
		(((A) & 0x0000ff00) << 8) | \
		(((A) & 0x000000ff) << 24))
/*
#	define HtoNS(A) ((((A) & 0xff00) >> 8) | (((A) & 0x00ff) << 8))
#	define NtoHS HtoNS
#	define NtoHL HtoNL
*/
/*#else
error One of BIG_ENDIAN or LITTLE_ENDIAN must be #defined.*/
#endif

#define ByteSwap5(x) ByteSwap((unsigned char *) &x,sizeof(x))

#if !defined(_BIG_ENDIAN) 
STATIC void ByteSwap(unsigned char * b, int n)
{
	int i = 0;
	int j = n-1;
	while (i<j)
	{
		char temp = b[i];
		b[i] = b[j];
		b[j] = temp;
		/* std::swap(b[i], b[j]); */
		i++, j--;
	}
	return;
}
#endif

static FILE *ioFITS_OUTPUT;
static long bytesAdded = 0;
static long bitpix = 8;
static long pcount = 0;
static long gcount = 1;
static long maxParamValues = 0;
const char ModelUnits[2][17] = {"'dimensionless '", "'photons/cm^2/s'" };

static const long bitpix_img = -64;

STATIC void punchFITS_PrimaryHeader( bool lgAddModel, bool lgNormalize );
STATIC void punchFITS_ParamHeader( /* long *numParamValues, */ long nintparm, long naddparm );
STATIC void punchFITS_ParamData( const vector<string>& paramNames, vector<long>& paramMethods,
								 const multi_arr<realnum,2>& paramRange, const multi_arr<realnum,2>& paramData,
								 long nintparm, long naddparm, long *numParamValues );
STATIC void punchFITS_EnergyHeader( long numEnergies );
STATIC void punchFITS_EnergyData( long ipLo, long ipHi );
STATIC void punchFITS_SpectraHeader( bool lgAdditiveModel, bool lgNormalize, long nintparm, long naddparm,
				     long totNumModels, long numEnergies );
STATIC void punchFITS_SpectraData( const multi_arr<realnum,2>& interpParameters, multi_arr<realnum,3>& theSpectrum,
					 int option, long totNumModels, long ipLo, long ipHi, long ipNorm, long nintparm, long naddparm );
STATIC void punchFITS_GenericHeader();
STATIC void punchFITS_GenericData();
STATIC void writeCloudyDetails( void );
STATIC long addComment( const string& CommentToAdd );
STATIC long addKeyword_txt( const char *theKeyword, const void *theValue, const char *theComment, long Str_Or_Log );
STATIC long addKeyword_num( const char *theKeyword, long theValue, const char *theComment);
STATIC long addKeyword_dbl( const char *theKeyword, const double theValue, const char *theComment);
inline string int2string(int val);

void saveFITSfile( FILE* ioPUN, int option, realnum Elo, realnum Ehi, realnum Enorm )
{
	DEBUG_ENTRY( "saveFITSfile()" );

	ASSERT( option >= 0 && option <= NUM_OUTPUT_TYPES );

	if( !grid.lgGridDone && option != NUM_OUTPUT_TYPES )
	{
		// save spectrum into intermediate binary file, these results will
		// be gathered at the end of the grid run into a proper FITS file
		GridRetrieveXSPECData(option);
		wr_block(&grid.Spectra[option][optimize.nOptimiz][0],
			 size_t(rfield.nflux)*sizeof(decltype(grid.Spectra[0][0][0])),
			 ioPUN);
		return;
	}

	ioFITS_OUTPUT = ioPUN;

	if( false )
	{
		FILE* asciiDump = open_data( "gridspectra.con", "w" );
		for( long i=0; i < rfield.nflux; i++ )
		{
			fprintf( asciiDump, "%.7e\t", rfield.anu(i) );
			for( long j=0; j < grid.totNumModels; j++ )
			{
				fprintf( asciiDump, "%.7e\t", grid.Spectra[4][j][i] );
			}
			fprintf( asciiDump, "\n" );
		}
		fclose( asciiDump );
	}

	/* This is generic FITS option */
	if( option == NUM_OUTPUT_TYPES )
	{
		punchFITS_PrimaryHeader( false, false );
		punchFITS_GenericHeader();
		punchFITS_GenericData();
	}
	/* These are specially designed XSPEC outputs. */
	/* the code below will only be executed during the gather phase of the grid */
	else if( option < NUM_OUTPUT_TYPES )
	{
		/* option 10 is exp(-tau). */
		/* false says not an additive model */
		bool lgAdditiveModel = ( option != 10 );
		bool lgNormalize = ( Enorm > 0.f );

		long ipLo = ( Elo > 0.f ) ? rfield.ipointC(Elo) : 0;
		long ipHi = ( Ehi > 0.f ) ? rfield.ipointC(Ehi) : rfield.nflux-1;
		long ipNorm = ( Enorm > 0.f ) ? rfield.ipointC(Enorm) : -1;
		long numEnergies = ipHi-ipLo+1;

		punchFITS_PrimaryHeader( lgAdditiveModel, lgNormalize );

		for( long i=0; i<grid.nintparm+grid.naddparm; i++ )
		{
			maxParamValues = MAX2( maxParamValues, grid.numParamValues[i] );
		}

		ASSERT( maxParamValues >= 2 );

		punchFITS_ParamHeader( /* grid.numParamValues, */ grid.nintparm, grid.naddparm );
		punchFITS_ParamData( grid.paramNames, grid.paramMethods, grid.paramRange, grid.paramData,
			grid.nintparm, grid.naddparm, grid.numParamValues );
		punchFITS_EnergyHeader( numEnergies );
		punchFITS_EnergyData( ipLo, ipHi );
		punchFITS_SpectraHeader( lgAdditiveModel, lgNormalize, grid.nintparm, grid.naddparm,
			grid.totNumModels, numEnergies );
		punchFITS_SpectraData( grid.interpParameters, grid.Spectra, option, grid.totNumModels,
				       ipLo, ipHi, ipNorm, grid.nintparm, grid.naddparm );
	}
}

STATIC void punchFITS_PrimaryHeader( bool lgAddModel, bool lgNormalize )
{
	const char *ModelName = "'CLOUDY'";

	DEBUG_ENTRY( "punchFITS_PrimaryHeader()" );

	int iunit = ( lgAddModel && !lgNormalize ) ? 1 : 0;

	bytesAdded = 0;

	fixit("bitpix is wrong when realnum is double?");

	bytesAdded += addKeyword_txt( "SIMPLE"	, "T",					"file does conform to FITS standard", 1 );
	bytesAdded += addKeyword_num( "BITPIX"	, bitpix,				"number of bits per data pixel" );
	bytesAdded += addKeyword_num( "NAXIS"	, 0,					"number of data axes" );
	bytesAdded += addKeyword_txt( "EXTEND"	, "T",					"FITS dataset may contain extensions", 1 );
	bytesAdded += addKeyword_txt( "CONTENT" , "'MODEL   '",			"spectrum file contains time intervals and event", 0 );
	bytesAdded += addKeyword_txt( "MODLNAME", ModelName,			"Model name", 0 );
	bytesAdded += addKeyword_txt( "MODLUNIT", ModelUnits[iunit],		"Model units", 0 );
	bytesAdded += addKeyword_txt( "REDSHIFT", "T",				"If true then redshift will be included as a par", 1 );
	if( lgAddModel == true )
	{
		bytesAdded += addKeyword_txt( "ADDMODEL", "T",				"If true then this is an additive table model", 1 );
	}
	else
	{
		bytesAdded += addKeyword_txt( "ADDMODEL", "F",				"If true then this is an additive table model", 1 );
	}

	/* bytes are added here as well */
	writeCloudyDetails();

	bytesAdded += addKeyword_txt( "HDUCLASS", "'OGIP    '",			"Format conforms to OGIP/GSFC conventions", 0 );
	bytesAdded += addKeyword_txt( "HDUCLAS1", "'XSPEC TABLE MODEL'","Extension contains an image", 0 );
	bytesAdded += addKeyword_txt( "HDUVERS"	, "'1.0.0   '",			"Version of format (OGIP memo OGIP-92-001)", 0 );
	/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	/* Now add blanks */
	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITS_ParamHeader( /* long *numParamValues, */ long nintparm, long naddparm )
{
	long numFields = 10;
	long naxis, naxis1, naxis2;

	DEBUG_ENTRY( "punchFITS_ParamHeader()" );

	ASSERT( nintparm+naddparm <= LIMPAR );

	/* Make sure the previous blocks are the right size */
	ASSERT( bytesAdded%RECORDSIZE == 0 );

	naxis = 2;
	/* >>chng 06 aug 23, change to maximum number of parameter values. */
	naxis1 = 44+4*maxParamValues;
	naxis2 = nintparm+naddparm;

	bytesAdded += addKeyword_txt( "XTENSION", "'BINTABLE'",			"binary table extension", 0  );
	bytesAdded += addKeyword_num( "BITPIX"	, bitpix,					"8-bit bytes" );
	bytesAdded += addKeyword_num( "NAXIS"	, naxis,					"2-dimensional binary table" );
	bytesAdded += addKeyword_num( "NAXIS1"	, naxis1,				"width of table in bytes" );
	bytesAdded += addKeyword_num( "NAXIS2"	, naxis2,					"number of rows in table" );
	bytesAdded += addKeyword_num( "PCOUNT"	, pcount,					"size of special data area" );
	bytesAdded += addKeyword_num( "GCOUNT"	, gcount,					"one data group (required keyword)" );
	bytesAdded += addKeyword_num( "TFIELDS"	, numFields,			"number of fields in each row" );
	bytesAdded += addKeyword_txt( "TTYPE1"	, "'NAME    '",			"label for field   1", 0  );
	bytesAdded += addKeyword_txt( "TFORM1"	, "'12A     '",			"data format of the field: ASCII Character", 0  );
	bytesAdded += addKeyword_txt( "TTYPE2"	, "'METHOD  '",			"label for field   2", 0  );
	bytesAdded += addKeyword_txt( "TFORM2"	, "'J       '",			"data format of the field: 4-byte INTEGER", 0  );
	bytesAdded += addKeyword_txt( "TTYPE3"	, "'INITIAL '",			"label for field   3", 0  );
	bytesAdded += addKeyword_txt( "TFORM3"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE4"	, "'DELTA   '",			"label for field   4", 0  );
	bytesAdded += addKeyword_txt( "TFORM4"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE5"	, "'MINIMUM '",			"label for field   5", 0  );
	bytesAdded += addKeyword_txt( "TFORM5"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE6"	, "'BOTTOM  '",			"label for field   6", 0  );
	bytesAdded += addKeyword_txt( "TFORM6"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE7"	, "'TOP     '",			"label for field   7", 0  );
	bytesAdded += addKeyword_txt( "TFORM7"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE8"	, "'MAXIMUM '",			"label for field   8", 0  );
	bytesAdded += addKeyword_txt( "TFORM8"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE9"	, "'NUMBVALS'",			"label for field   9", 0  );
	bytesAdded += addKeyword_txt( "TFORM9"	, "'J       '",			"data format of the field: 4-byte INTEGER", 0  );
	bytesAdded += addKeyword_txt( "TTYPE10"	, "'VALUE   '",			"label for field  10", 0  );

	/* >>chng 06 aug 23, use maxParamValues instead of numParamValues */
	/* The size of this array is dynamic, set to size of the maximum of the numParamValues vector */
	string theValue = int2string(maxParamValues);
	bytesAdded += addKeyword_txt( "TFORM10"	, theValue.c_str(),		"data format of the field: 4-byte REAL", 0  );

	bytesAdded += addKeyword_txt( "EXTNAME"	, "'PARAMETERS'",		"name of this binary table extension", 0  );
	bytesAdded += addKeyword_txt( "HDUCLASS", "'OGIP    '",			"Format conforms to OGIP/GSFC conventions", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS1", "'XSPEC TABLE MODEL'","model spectra for XSPEC", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS2", "'PARAMETERS'",		"Extension containing paramter info", 0  );
	bytesAdded += addKeyword_txt( "HDUVERS"	, "'1.0.0   '",			"Version of format (OGIP memo OGIP-92-001)", 0  );
	bytesAdded += addKeyword_num( "NINTPARM", nintparm,				"Number of interpolation parameters" );
	bytesAdded += addKeyword_num( "NADDPARM", naddparm,				"Number of additional parameters" );
	/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	/* Now add blanks */
	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITS_ParamData(const vector<string>& paramNames,
								vector<long>& paramMethods,
								const multi_arr<realnum,2>& paramRange,
								const multi_arr<realnum,2>& paramData,
								long nintparm,
								long naddparm,
								long *numParamValues )
{
	long i, j;

	DEBUG_ENTRY( "punchFITS_ParamData()" );

	ASSERT( nintparm+naddparm <= LIMPAR );

	/* Now add the parameters data */
	for( i=0; i<nintparm+naddparm; i++ )
	{
		int32 numTemp;

#define	LOG2LINEAR 0

		paramMethods[i] = HtoNL(paramMethods[i]);
		/* >>chng 06 aug 23, numParamValues is now an array.  */
		numTemp = HtoNL(numParamValues[i]);

#if LOG2LINEAR
		/* change to linear */
		paramRange[i][0] = (realnum)exp10(  (double)paramRange[i][0] );
		paramRange[i][1] = (realnum)exp10(  (double)paramRange[i][1] );
		paramRange[i][2] = (realnum)exp10(  (double)paramRange[i][2] );
		paramRange[i][3] = (realnum)exp10(  (double)paramRange[i][3] );
		paramRange[i][4] = (realnum)exp10(  (double)paramRange[i][4] );
		paramRange[i][5] = (realnum)exp10(  (double)paramRange[i][5] );
#endif

#if !defined(_BIG_ENDIAN) 
		ByteSwap5( paramRange[i][0] );
		ByteSwap5( paramRange[i][1] );
		ByteSwap5( paramRange[i][2] );
		ByteSwap5( paramRange[i][3] );
		ByteSwap5( paramRange[i][4] );
		ByteSwap5( paramRange[i][5] );
#endif

		/* >>chng 06 aug 23, numParamValues is now an array.  */
		for( j=0; j<numParamValues[i]; j++ )
		{

#if LOG2LINEAR
			paramData[i][j] = (realnum)exp10(  (double)paramData[i][j] );
#endif

#if !defined(_BIG_ENDIAN) 
			ByteSwap5( paramData[i][j] );
#endif
		}

		bytesAdded += fprintf(ioFITS_OUTPUT, "%-12s", paramNames[i].substr(0,12).c_str() );
		bytesAdded += (long)fwrite( &paramMethods[i],	1,				  sizeof(int32),   ioFITS_OUTPUT );
		bytesAdded += (long)fwrite( &paramRange[i][0],	1,				6*sizeof(realnum),   ioFITS_OUTPUT );
		bytesAdded += (long)fwrite( &numTemp,			1,				  sizeof(int32),   ioFITS_OUTPUT );
		/* >>chng 06 aug 23, numParamValues is now an array.  */
		bytesAdded += (long)fwrite( &paramData[i][0],	1, (unsigned)numParamValues[i]*sizeof(realnum), ioFITS_OUTPUT );

		for( j=numParamValues[i]+1; j<=maxParamValues; j++ )
		{
			realnum filler = -10.f;
			bytesAdded += (long)fwrite( &filler,		1, sizeof(realnum),   ioFITS_OUTPUT );
		}
	}

	/* Switch the endianness again */
	for( i=0; i<nintparm+naddparm; i++ )
	{
		paramMethods[i] = HtoNL(paramMethods[i]);

#if !defined(_BIG_ENDIAN) 
		ByteSwap5( paramRange[i][0] );
		ByteSwap5( paramRange[i][1] );
		ByteSwap5( paramRange[i][2] );
		ByteSwap5( paramRange[i][3] );
		ByteSwap5( paramRange[i][4] );
		ByteSwap5( paramRange[i][5] );
#endif

		/* >>chng 06 aug 23, numParamValues is now an array.  */
		for( j=0; j<numParamValues[i]; j++ )
		{
#if !defined(_BIG_ENDIAN) 
			ByteSwap5( paramData[i][j] );
#endif
		}
	}

	while( bytesAdded%RECORDSIZE > 0 )
	{
		int	tempInt = 0;
		bytesAdded += (long)fwrite( &tempInt, 1, 1,   ioFITS_OUTPUT );
	}
	return;
}

STATIC void punchFITS_EnergyHeader( long numEnergies )
{
	long numFields = 2;
	long naxis, naxis1, naxis2;

	DEBUG_ENTRY( "punchFITS_EnergyHeader()" );

	/* Make sure the previous blocks are the right size */
	ASSERT( bytesAdded%RECORDSIZE == 0 );

	naxis = 2;
	naxis1 = 2*sizeof(realnum);
	naxis2 = numEnergies;

	bytesAdded += addKeyword_txt( "XTENSION", "'BINTABLE'",			"binary table extension", 0 );
	bytesAdded += addKeyword_num( "BITPIX"	, bitpix,				"8-bit bytes" );
	bytesAdded += addKeyword_num( "NAXIS"	, naxis,				"2-dimensional binary table" );
	bytesAdded += addKeyword_num( "NAXIS1"	, naxis1,				"width of table in bytes" );
	bytesAdded += addKeyword_num( "NAXIS2"	, naxis2,				"number of rows in table" );
	bytesAdded += addKeyword_num( "PCOUNT"	, pcount,				"size of special data area" );
	bytesAdded += addKeyword_num( "GCOUNT"	, gcount,				"one data group (required keyword)" );
	bytesAdded += addKeyword_num( "TFIELDS"	, numFields,			"number of fields in each row" );
	bytesAdded += addKeyword_txt( "TTYPE1"	, "'ENERG_LO'",			"label for field   1", 0  );
	bytesAdded += addKeyword_txt( "TFORM1"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE2"	, "'ENERG_HI'",			"label for field   2", 0  );
	bytesAdded += addKeyword_txt( "TFORM2"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "EXTNAME"	, "'ENERGIES'",			"name of this binary table extension", 0  );
	bytesAdded += addKeyword_txt( "HDUCLASS", "'OGIP    '",			"Format conforms to OGIP/GSFC conventions", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS1", "'XSPEC TABLE MODEL'","model spectra for XSPEC", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS2", "'ENERGIES'",			"Extension containing energy bin info", 0  );
	bytesAdded += addKeyword_txt( "HDUVERS"	, "'1.0.0   '",			"Version of format (OGIP memo OGIP-92-001)", 0  );
	/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITS_EnergyData( long ipLo, long ipHi )
{
	DEBUG_ENTRY( "punchFITS_EnergyData()" );

	/* Now add the energies data */
	for( long i=ipLo; i <= ipHi; i++ )
	{
		/* Convert to kev */
		realnum EnergyLow = realnum(EVRYD*rfield.anumin(i)/1000.);
		realnum EnergyHi = realnum(EVRYD*rfield.anumax(i)/1000.);

#if !defined(_BIG_ENDIAN) 
		ByteSwap5(EnergyLow);
		ByteSwap5(EnergyHi);
#endif

		bytesAdded += (long)fwrite( &EnergyLow, 1, sizeof(realnum), ioFITS_OUTPUT );
		bytesAdded += (long)fwrite( &EnergyHi, 1, sizeof(realnum), ioFITS_OUTPUT );
	}

	int tempInt = 0;
	while( bytesAdded%RECORDSIZE > 0 )
		bytesAdded += (long)fwrite( &tempInt, 1, 1, ioFITS_OUTPUT );
}

STATIC void punchFITS_SpectraHeader( bool lgAddModel, bool lgNormalize, long nintparm, long naddparm,
				     long totNumModels, long numEnergies )
{
	long i, numFields = 2+naddparm;
	long naxis, naxis1, naxis2;
	char theKeyword1[30];
	char theKeyword2[30];
	char theKeyword3[30];
	char theComment1[47];

	DEBUG_ENTRY( "punchFITS_SpectraHeader()" );

	ASSERT( nintparm + naddparm <= LIMPAR );

	/* Make sure the previous blocks are the right size */
	ASSERT( bytesAdded%RECORDSIZE == 0 );

	naxis = 2;
	naxis1 = ( numEnergies*(naddparm+1) + nintparm ) * (long)sizeof(realnum);
	naxis2 = totNumModels; 
	int iunit = ( lgAddModel && !lgNormalize ) ? 1 : 0;

	bytesAdded += addKeyword_txt( "XTENSION", "'BINTABLE'",			"binary table extension", 0  );
	bytesAdded += addKeyword_num( "BITPIX"	, bitpix,				"8-bit bytes" );
	bytesAdded += addKeyword_num( "NAXIS"	, naxis,				"2-dimensional binary table" );
	bytesAdded += addKeyword_num( "NAXIS1"	, naxis1,				"width of table in bytes" );
	bytesAdded += addKeyword_num( "NAXIS2"	, naxis2,				"number of rows in table" );
	bytesAdded += addKeyword_num( "PCOUNT"	, pcount,				"size of special data area" );
	bytesAdded += addKeyword_num( "GCOUNT"	, gcount,				"one data group (required keyword)" );
	bytesAdded += addKeyword_num( "TFIELDS"	, numFields,			"number of fields in each row" );

	/******************************************/
	/* These are the interpolation parameters */
	/******************************************/
	bytesAdded += addKeyword_txt( "TTYPE1"	, "'PARAMVAL'",			"label for field   1", 0 );
	/* The size of this array is dynamic, set to size of nintparm */
	string theValue2 = int2string(nintparm);
	bytesAdded += addKeyword_txt( "TFORM1"	, theValue2.c_str(),		"data format of the field: 4-byte REAL", 0  );

	/******************************************/
	/* This is the interpolated spectrum      */	
	/******************************************/
	bytesAdded += addKeyword_txt( "TTYPE2"	, "'INTPSPEC'",	"label for field 2", 0  );
	/* The size of this array is dynamic, set to size of numEnergies */
	string theValue = int2string(numEnergies);
	bytesAdded += addKeyword_txt( "TFORM2"	, theValue.c_str(),		"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TUNIT2"	, ModelUnits[iunit],		"physical unit of field", 0  );

	/******************************************/
	/* These are the additional parameters    */
	/******************************************/
	for( i=1; i<=naddparm; i++ )
	{
		sprintf( theKeyword1,	"%s%ld", "TTYPE", i+2 );
		sprintf( theKeyword2,	"%s%ld", "TFORM", i+2 );
		sprintf( theKeyword3,	"%s%ld", "TUNIT", i+2 );

		ostringstream theValue1;
		theValue1 << "'ADDSP" << setw(2) << setfill('0') << i << "'";

		sprintf( theComment1,	"%s%ld", "label for field ", i+2 );

		bytesAdded += addKeyword_txt( theKeyword1	, theValue1.str().c_str(),	theComment1, 0  );
		bytesAdded += addKeyword_txt( theKeyword2	, theValue.c_str(),		"data format of the field: 4-byte REAL", 0  );
		bytesAdded += addKeyword_txt( theKeyword3	, ModelUnits[iunit],	"physical unit of field", 0  );
	}

	bytesAdded += addKeyword_txt( "EXTNAME"	, "'SPECTRA '",			"name of this binary table extension", 0  );
	bytesAdded += addKeyword_txt( "HDUCLASS", "'OGIP    '",			"Format conforms to OGIP/GSFC conventions", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS1", "'XSPEC TABLE MODEL'","model spectra for XSPEC", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS2", "'MODEL SPECTRA'",	"Extension containing model spectra", 0  );
	bytesAdded += addKeyword_txt( "HDUVERS"	, "'1.0.0   '",			"Version of format (OGIP memo OGIP-92-001)", 0  );
	/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITS_SpectraData( const multi_arr<realnum,2>& interpParameters, multi_arr<realnum,3>& theSpectrum,
				   int option, long totNumModels, long ipLo, long ipHi, long ipNorm, long nintparm, long naddparm )
{
	long i;
	long naxis2 = totNumModels;

	DEBUG_ENTRY( "punchFITS_SpectraData()" );

	ASSERT( nintparm + naddparm <= LIMPAR );

	/* Now add the spectra data */
	for( i=0; i < naxis2; i++ )
	{
		realnum fluxNorm = 0.f;
		if( ipNorm >= 0 )
		{
			// normalize spectrum to 1 photons/cm^2/s/keV
			realnum binwidth_keV = realnum(Energy(rfield.widflx(ipNorm)).keV());
			fluxNorm = theSpectrum[option][i][ipNorm]/binwidth_keV;
		}
		flex_arr<realnum> flux(ipLo, ipHi+1);
		for( long j=ipLo; j <= ipHi; j++ )
		{
			flux[j] = theSpectrum[option][i][j];
			if( fluxNorm > 0.f )
				flux[j] /= fluxNorm;
		}

#if !defined(_BIG_ENDIAN)
		for( long j=ipLo; j <= ipHi; j++ )
		{
			ByteSwap5( flux[j] );
		}

		for( long j = 0; j < nintparm; j++ )
		{
			ByteSwap5( interpParameters[i][j] );
		}
#endif

		/* The interpolation parameters vector */
		bytesAdded += (long)fwrite( &interpParameters[i][0], 1, (unsigned)nintparm*sizeof(realnum), ioFITS_OUTPUT );
		/* The interpolated spectrum */
		bytesAdded += (long)fwrite( &flux[ipLo], 1, (unsigned)(ipHi-ipLo+1)*sizeof(realnum), ioFITS_OUTPUT );

#if !defined(_BIG_ENDIAN)
		/* Switch the endianness back to native. */
		for( long j=0; j < nintparm; j++ )
		{
			ByteSwap5( interpParameters[i][j] );
		}
#endif

		/* >>chng 06 aug 23, disable additional parameters for now */
		if( naddparm > 0 )
		{
			/* The additional parameters */
			/* \todo 2	must create another array if we are to save additional parameter information. */
			fprintf( ioQQQ, " Additional parameters not currently supported.\n" );
			cdEXIT( EXIT_FAILURE );
		}
	}

	int tempInt = 0;
	while( bytesAdded%RECORDSIZE > 0 )
		bytesAdded += (long)fwrite( &tempInt, 1, 1, ioFITS_OUTPUT );
}

STATIC void punchFITS_GenericHeader()
{
	long numFields = 2;
	long naxis, naxis1, naxis2;

	DEBUG_ENTRY( "punchFITS_GenericHeader()" );

	/* Make sure the previous blocks are the right size */
	ASSERT( bytesAdded%RECORDSIZE == 0 );

	naxis = 2;
	naxis1 = numFields*(long)sizeof(realnum);
	naxis2 = rfield.nflux;

	bytesAdded += addKeyword_txt( "XTENSION", "'BINTABLE'",			"binary table extension", 0 );
	bytesAdded += addKeyword_num( "BITPIX"	, bitpix,				"8-bit bytes" );
	bytesAdded += addKeyword_num( "NAXIS"	, naxis,				"2-dimensional binary table" );
	bytesAdded += addKeyword_num( "NAXIS1"	, naxis1,				"width of table in bytes" );
	bytesAdded += addKeyword_num( "NAXIS2"	, naxis2,				"number of rows in table" );
	bytesAdded += addKeyword_num( "PCOUNT"	, pcount,				"size of special data area" );
	bytesAdded += addKeyword_num( "GCOUNT"	, gcount,				"one data group (required keyword)" );
	bytesAdded += addKeyword_num( "TFIELDS"	, numFields,			"number of fields in each row" );
	bytesAdded += addKeyword_txt( "TTYPE1"	, "'ENERGY  '",			"label for field   1", 0  );
	bytesAdded += addKeyword_txt( "TFORM1"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "TTYPE2"	, "'TRN_SPEC'",			"label for field   2", 0  );
	bytesAdded += addKeyword_txt( "TFORM2"	, "'E       '",			"data format of the field: 4-byte REAL", 0  );
	bytesAdded += addKeyword_txt( "EXTNAME"	, "'SPECTRA '",			"name of this binary table extension", 0  );
	bytesAdded += addKeyword_txt( "HDUCLASS", "'OGIP    '",			"Format conforms to OGIP/GSFC conventions", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS1", "'XSPEC TABLE MODEL'","model spectra for XSPEC", 0  );
	bytesAdded += addKeyword_txt( "HDUCLAS2", "'ENERGIES'",			"Extension containing energy bin info", 0  );
	bytesAdded += addKeyword_txt( "HDUVERS"	, "'1.0.0   '",			"Version of format (OGIP memo OGIP-92-001)", 0  );
	/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITS_GenericData()
{
	DEBUG_ENTRY( "punchFITS_GenericData()" );

	vector<realnum> TransmittedSpectrum(rfield.nflux);

	cdSPEC2( 8, get_ptr(TransmittedSpectrum) );

	/* Now add the energies data */
	for( long j=0; j < rfield.nflux; j++ )
	{
		realnum Energy = rfield.anu(j);

#if !defined(_BIG_ENDIAN) 
		ByteSwap5(Energy);
		ByteSwap5(TransmittedSpectrum[j]);
#endif

		bytesAdded += (long)fwrite( &Energy, 1, sizeof(realnum), ioFITS_OUTPUT );
		bytesAdded += (long)fwrite( &TransmittedSpectrum[j], 1, sizeof(realnum), ioFITS_OUTPUT );
	}

	int tempInt = 0;
	while( bytesAdded%RECORDSIZE > 0 )
		bytesAdded += (long)fwrite( &tempInt, 1, 1, ioFITS_OUTPUT );
}

STATIC void writeCloudyDetails( void )
{
	char timeString[30]="";
	char tempString[70];
	time_t now;

	/* usually print date and time info - do not if "no times" command entered, 
	 * which set this flag false */
	now = time(NULL);
	if( prt.lgPrintTime ) 
	{
		/* now add date of this run */
		/* now print this time at the end of the string.  the system put cr at the end of the string */
		strcpy( timeString , ctime(&now) );
	}
	/* ctime puts a carriage return at the end, but we can't have that in a fits file.
	 * remove the carriage return here. */
	for( long i=0; i<30; i++ )
	{
		if( timeString[i] == '\n' )
		{
			timeString[i] = ' ';
		}
	}

	bytesAdded += addComment( "Generated by Cloudy " + t_version::Inst().chVersion );
	bytesAdded += addComment( t_version::Inst().chInfo );
	strcpy( tempString, "--- " );
	strcat( tempString, timeString );
	bytesAdded += addComment( tempString );
	bytesAdded += addComment( "Input string was as follows: " );
	for( size_t i=0; i < input.crd.size(); i++ )
	{
		char firstLine[70], extraLine[64];

		// exclude lines from init files
		if( input.crd[i]->InclLevel > 0 )
			continue;

		size_t j = input.crd[i]->chCardSav.length();
		const char* chCardSav = input.crd[i]->chCardSav.c_str();

		strncpy(firstLine, chCardSav, sizeof(firstLine));
		size_t k = sizeof(firstLine)-1;
		firstLine[k] = '\0';
		bytesAdded += addComment( firstLine );
		while( j > k )
		{
			strncpy(extraLine, chCardSav+k, sizeof(extraLine));
			size_t l = sizeof(extraLine)-1;
			extraLine[l] = '\0';
			strcpy( tempString, "cont> " );
			strcat( tempString, extraLine );
			bytesAdded += addComment( tempString );
			k += l;
		}
	}
}

STATIC long addKeyword_txt( const char *theKeyword, const void *theValue, const char *theComment, long Str_Or_Log )
{
	long numberOfBytesWritten = 0;

	DEBUG_ENTRY( "addKeyword_txt()" );

	/* False means string, true means logical */
	if( Str_Or_Log == 0 )
	{
		numberOfBytesWritten = fprintf(ioFITS_OUTPUT, "%-8s%-2s%-20s%3s%-47s",
			theKeyword,
			"= ",
			(char *)theValue,
			" / ",
			theComment );
	}
	else
	{
		ASSERT( Str_Or_Log == 1 );
		numberOfBytesWritten = fprintf(ioFITS_OUTPUT, "%-8s%-2s%20s%3s%-47s",
			theKeyword,
			"= ",
			(char *)theValue,
			" / ",
			theComment );
	}

	ASSERT( numberOfBytesWritten%LINESIZE == 0 );
	return numberOfBytesWritten;
}

STATIC long addKeyword_num( const char *theKeyword, long theValue, const char *theComment)
{
	long numberOfBytesWritten = 0;

	DEBUG_ENTRY( "addKeyword_num()" );

	numberOfBytesWritten = fprintf(ioFITS_OUTPUT, "%-8s%-2s%20ld%3s%-47s",
		theKeyword,
		"= ",
		theValue,
		" / ",
		theComment );

	ASSERT( numberOfBytesWritten%LINESIZE == 0 );
	return numberOfBytesWritten;
}

STATIC long addKeyword_dbl( const char *theKeyword, const double theValue, const char *theComment)
{
	long numberOfBytesWritten = 0;

	DEBUG_ENTRY( "addKeyword_dbl()" );

	numberOfBytesWritten = fprintf(ioFITS_OUTPUT, "%-8s%-2s%20.4e%3s%-47s",
		theKeyword,
		"= ",
		theValue,
		" / ",
		theComment );

	ASSERT( numberOfBytesWritten%LINESIZE == 0 );
	return numberOfBytesWritten;
}

long addComment( const string& CommentToAdd )
{
	DEBUG_ENTRY( "addComment()" );

	// pad with spaces to make sure that tempString.length() >= 80
	string tempString = "COMMENT   " + CommentToAdd.substr(0, 69) + string(70, ' ');

	/* tabs violate FITS standard, replace them with spaces. */
	for( size_t i=10; i < LINESIZE; i++ )
	{
		if( tempString[i] == '\t' )
		{
			tempString[i] = ' ';
		}
	}

	long numberOfBytesWritten = fprintf(ioFITS_OUTPUT, "%s", tempString.substr(0,LINESIZE).c_str() );

	ASSERT( numberOfBytesWritten%LINESIZE == 0 );
	return numberOfBytesWritten;
}

inline string int2string(int val)
{
	DEBUG_ENTRY( "int2string()" );

	ostringstream oss1;
	oss1 << val << "E";
	ostringstream oss2;
	oss2 << "'" << left << setw(8) << oss1.str() << "'";
	return oss2.str();
}

//------------------------------------------------------------------------------
STATIC void get_data_minmax( const long nPixels,
				const multi_arr<double,2,C_TYPE> &a,
				double &datamin, double &datamax )
{
	DEBUG_ENTRY( "get_data_minmax()" );

	datamax = -1e300;
	datamin = 1e300;

	for( long irow = 0; irow < nPixels; irow++ )
	{
		for( long icol = 0; icol < nPixels; icol++ )
		{
			if( a[irow][icol] > datamax )
				datamax = a[irow][icol];
			if( a[irow][icol] < datamin )
				datamin = a[irow][icol];
		}
	}

	return;
}

STATIC void get_data_minmax( const long nCols, const valarray<double> &a,
				double &datamin, double &datamax )
{
	DEBUG_ENTRY( "get_data_minmax()" );

	datamax = -1e300;
	datamin = 1e300;

	for( long icol = 0; icol < nCols; icol++ )
	{
		if( a[icol] > datamax )
			datamax = a[icol];
		if( a[icol] < datamin )
			datamin = a[icol];
	}

	return;
}

STATIC void punchFITSimg_PrimaryHeader()
{
	DEBUG_ENTRY( "punchFITSimg_PrimaryHeader()" );

	bytesAdded = 0;

	bytesAdded += addKeyword_txt( "SIMPLE", "T",	"file does conform to FITS standard", 1 );
	bytesAdded += addKeyword_num( "BITPIX", bitpix,	"number of bits per data pixel" );
	bytesAdded += addKeyword_num( "NAXIS" , 0,	"number of data axes" );
	bytesAdded += addKeyword_txt( "EXTEND", "T",	"FITS dataset may contain extensions", 1 );

	/* bytes are added here as well */
	writeCloudyDetails();

	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded % LINESIZE == 0 );

	/* Now add blanks */
	while( bytesAdded % RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITSimg_ExtensionHeader( const string &extName,
					const string &units,
					const long nCols, const long nRows,
					const double datamin, const double datamax )
{
	DEBUG_ENTRY( "punchFITSimg_ExtensionHeader()" );

	string name = "'" + extName + "'";
	string bunits = "'" + units + "'";

	bytesAdded += addKeyword_txt( "XTENSION", "'IMAGE   '",	"IMAGE extension", 0 );
	bytesAdded += addKeyword_num( "BITPIX"  , bitpix_img,	"number of bits per data pixel" );
	bytesAdded += addKeyword_num( "NAXIS"   , 2,		"number of data axes" );
	bytesAdded += addKeyword_num( "NAXIS1"  , nCols,	"length of data axis 1" );
	bytesAdded += addKeyword_num( "NAXIS2"  , nRows,	"length of data axis 2" );
	bytesAdded += addKeyword_num( "PCOUNT"  , 0,		"required keyword; must = 0" );
	bytesAdded += addKeyword_num( "GCOUNT"  , 1,		"required keyword; must = 1" );
	bytesAdded += addKeyword_txt( "EXTNAME" , name.c_str(),	"Extension Name", 0 );
	bytesAdded += addKeyword_num( "EXTVER"  , 1,		"Extension Version" );
	bytesAdded += addKeyword_dbl( "DATAMIN" , datamin,	"Lowest non-zero array value" );
	bytesAdded += addKeyword_dbl( "DATAMAX" , datamax,	"Highest array value" );
	bytesAdded += addKeyword_txt( "BUNIT"   , bunits.c_str(),"Image Units", 0 );

	/* After everything else */
	bytesAdded += fprintf(ioFITS_OUTPUT, "%-80s", "END" );

	ASSERT( bytesAdded%LINESIZE == 0 );

	/* Now add blanks */
	while( bytesAdded%RECORDSIZE > 0 )
	{
		bytesAdded += fprintf(ioFITS_OUTPUT, "%-1s", " " );
	}
	return;
}

STATIC void punchFITSimg_ExtensionData_matrix( const long nPixels,
					const multi_arr<double,2,C_TYPE> &image )
{
	DEBUG_ENTRY( "punchFITSimg_ExtensionData_matrix()" );

	for( long irow = 0; irow < nPixels; irow++ )
	{
		for( long icol = 0; icol < nPixels; icol++ )
		{
			double v = image[irow][icol];
#if !defined(_BIG_ENDIAN) 
			ByteSwap5(v);
#endif
			bytesAdded += (long)fwrite( &v, 1, sizeof(double), ioFITS_OUTPUT );
		}
	}

	int tempInt = 0;
	while( bytesAdded % RECORDSIZE > 0 )
		bytesAdded += (long)fwrite( &tempInt, 1, 1, ioFITS_OUTPUT );

	return;
}

STATIC void punchFITSimg_ExtensionData_vector( const long nCols,
					const valarray<double> &vec )
{
	DEBUG_ENTRY( "punchFITSimg_ExtensionData_vec()" );

	for( long icol = 0; icol < nCols; icol++ )
	{
		double v = vec[icol];
#if !defined(_BIG_ENDIAN) 
		ByteSwap5(v);
#endif
		bytesAdded += (long)fwrite( &v, 1, sizeof(double), ioFITS_OUTPUT );
	}

	int tempInt = 0;
	while( bytesAdded % RECORDSIZE > 0 )
		bytesAdded += (long)fwrite( &tempInt, 1, 1, ioFITS_OUTPUT );

	return;
}

void saveFITSimg( FILE *ioPUN, const string &extName, const string &units,
		const long nPixels, const multi_arr<double,2,C_TYPE> &image )
{
	DEBUG_ENTRY( "saveFITSimg()" );

	double datamin,
	       datamax;
	get_data_minmax( nPixels, image, datamin, datamax );

	ioFITS_OUTPUT = ioPUN;
	punchFITSimg_PrimaryHeader();
	punchFITSimg_ExtensionHeader( extName, units, nPixels, nPixels,
					datamin, datamax );
	punchFITSimg_ExtensionData_matrix( nPixels, image );
}

void saveFITSimg( FILE *ioPUN, const string &extName, const string &units,
			const long nCols, const valarray<double> &vec )
{
	DEBUG_ENTRY( "saveFITSimg()" );

	double datamin,
	       datamax;
	get_data_minmax( nCols, vec, datamin, datamax );

	ioFITS_OUTPUT = ioPUN;
	punchFITSimg_ExtensionHeader( extName, units, nCols, 1, datamin, datamax );
	punchFITSimg_ExtensionData_vector( nCols, vec );
}
