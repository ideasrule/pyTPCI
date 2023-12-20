/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "save.h"
t_save save;


void save_img_matrix::zero()
{
	t_prt_matrix::zero();
	lgImgRates = false;
	lgFITS = false;
	iteration = 0;
	zone = 0;
}

string save_img_matrix::set_basename( const string &this_species,
					const long this_iteration,
					const long this_zone,
	       				bool haveNegPop )
{	
	stringstream ss;
	ss << this_species << "_it" << this_iteration << "_nz" << this_zone;
	string basename = ss.str();

	if( haveNegPop )
	{
		basename += "_NEG_POP";
	}

	return basename;
}

void save_img_matrix::createImage(
				const long this_iteration,
				const long this_zone,
				const long numLevels,
				const multi_arr<double,2,C_TYPE> &matrix,
	       			const valarray<double> &creation,
				bool haveNegPop )
{
	createImage( species, this_iteration, this_zone, numLevels, matrix, creation, haveNegPop );
}

void save_img_matrix::createImage( const string &this_species,
				const long this_iteration,
				const long this_zone,
				const long numLevels,
				const multi_arr<double,2,C_TYPE> &matrix,
	       			const valarray<double> &creation,
				bool haveNegPop )
{
	DEBUG_ENTRY( "save_img_matrix::createImage()" );

	if( not lgLevelsResolved )
	{
		resolveLevels();
	}

	if( speciesLevelList.size() == 0 )
		return;

	long nlev = 0;
	multi_arr<double,2,C_TYPE> submatrix;
	valarray<double> creation_subv;

	if( speciesLevelList.front() == 1 &&
		speciesLevelList.back() == numLevels )
	{
		nlev = numLevels;
		submatrix = matrix;
		creation_subv = creation;
	}
	else
	{
		long ifront = speciesLevelList.front();
		if( ifront >= numLevels )
			ifront = numLevels-1;

		long ilast = speciesLevelList.back();
		if( ilast >= numLevels )
			ilast = numLevels - 1;

		nlev = ilast - ifront + 1;

		submatrix.alloc( nlev, nlev );
		creation_subv.resize( nlev );

		for( vector<long>::iterator ipLo = speciesLevelList.begin();
			 ipLo != speciesLevelList.end(); ++ipLo )
		{
			if( *ipLo >= numLevels )
				continue;

			long iy = *ipLo - ifront;
			if( iy >= nlev )
				continue;

			creation_subv[ iy ] = creation[ *ipLo ];

			for( vector<long>::iterator ipHi = speciesLevelList.begin();
				ipHi != speciesLevelList.end(); ++ipHi )
			{
				if( *ipHi >= numLevels )
					continue;

				long ix = *ipHi - ifront;
				if( ix >= nlev )
					continue;

				submatrix[ iy ][ ix ] = matrix[ *ipLo ][ *ipHi ];
			}
		}
	}

	string basename = set_basename( this_species, this_iteration, this_zone, haveNegPop );

	if( lgFITS )
	{
		createImage_FITS( basename, nlev, submatrix, creation_subv );
	}
	else
	{
		createImage_PPM( basename, nlev, submatrix );
	}
}


void save_img_matrix::createImage_FITS( const string &basename,
				const long numLevels,
				const multi_arr<double,2,C_TYPE> &matrix,
	       			const valarray<double> &creation )
{
	DEBUG_ENTRY( "save_img_matrix::createImage_PPM()" );

	string filename = basename + ".fits";

	FILE *fp = open_data( filename, "w" );
	saveFITSimg( fp, "MATRIX", "s^{-1}", numLevels, matrix );
	saveFITSimg( fp, "CREATION", "cm^{-3} s^{-1}", numLevels, creation );
	fclose( fp );
}

void save_img_matrix::addImagePop_FITS(
				const long this_iteration,
				const long this_zone,
				const long numLevels,
				const valarray<double> &pop,
				bool haveNegPop )
{
	addImagePop_FITS( species, this_iteration, this_zone, numLevels, pop, haveNegPop );
}

void save_img_matrix::addImagePop_FITS( const string &this_species,
				const long this_iteration,
				const long this_zone,
				const long numLevels,
				const valarray<double> &pop,
				bool haveNegPop )
{
	string basename = set_basename( this_species, this_iteration, this_zone, haveNegPop );
	string filename = basename + ".fits";

	FILE *fp = open_data( filename, "a" );
	saveFITSimg( fp, "POPULATIONS", "cm^{-3}", numLevels, pop );
	fclose( fp );
}

void save_img_matrix::createImage_PPM( const string &basename,
				const long numLevels,
				const multi_arr<double,2,C_TYPE> &matrix )
{
	DEBUG_ENTRY( "save_img_matrix::createImage_PPM()" );

	double amax = 0.0;

	for( long iy = 0; iy < numLevels; iy++ )
	{
		for( long ix = 0; ix < numLevels; ix++ )
		{
			amax = MAX2(amax, fabs( matrix[ iy ][ ix ] ) );
		}
	}

	string filename = basename + ".ppm";

	FILE *fp = open_data( filename, "w" );
	const int ipix=4;
	fprintf(fp, "P6\n%ld %ld\n255\n", ipix*numLevels, ipix*numLevels ); // PGM

	char buf[3*ipix+1];
	buf[3*ipix] = '\0';

	for( long iy = 0; iy < numLevels; iy++ )
	{
		for( int i1 = 0; i1 < ipix; ++i1 )
		{
			for( long ix = 0; ix < numLevels; ix++ )
			{
				unsigned char ch[3];
				if ( matrix[ iy ][ ix ] == 0.0 )
				{
					ch[0] = ch[1] = ch[2] = 0;
				}
				else
				{
					double val = 1.+log10(fabs( matrix[ iy ][ ix ] )/amax)/16.;
					ch[0] = (unsigned int)(MAX2(1.,MIN2(256.*val,255.)));
					ch[1] = ch[0];
					ch[2] = (255+ch[0])/2;
				}
				for( int i2=0; i2 < ipix; ++i2 )
				{
					buf[3*i2] = ch[0];
					buf[3*i2+1] = ch[1];
					buf[3*i2+2] = ch[2];
				}
				fwrite(buf,3*ipix,sizeof(char),fp);
			}
		}
	}
	fclose(fp);
}
