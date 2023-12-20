/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SetPrintLineCol set some parameters for the main line block & wl format */
/*prt_LineLabels save all labels and wavelengths for emission line array */
/*sprt_wl write wavelength to string - must be kept parallel with prt_wl */
/*prt_wl - print floating wavelength in Angstroms, in output format */
/* prt_line print the line label, followed by the line wl, and the wavelength of the closest match, if given */
/* prt_line_inlist print line suitable for output list, label not enclosed in quotation marks */
#include "cddefines.h"
#include "lines.h"
#include "prt.h"
#include "generic_state.h"
#include "save.h"

t_prt prt;
t_line_col prt_linecol;

void t_line_col::zero()
{
	DEBUG_ENTRY( "t_line_col::zero()" );
	/* set the format of the main output line list */
	absint_len = 8;
	relint_len = 9;
	col_gap_len = 6;
}


/*SetPrintLineCol set some parameters for the main line block & wl format */
void SetPrintLineCol ()
{
	LineSave.wl_length = (int) LineSave.sig_figs + 2;

	/* set the format of the main output line list */
	prt_linecol.column_len = (NCHLAB-1) + LineSave.wl_length + prt_linecol.absint_len + prt_linecol.relint_len + 3;

	prt_linecol.col_gap.assign( prt_linecol.col_gap_len, ' ' );
	prt_linecol.relint_outrange.assign( prt_linecol.relint_len , '*' );

	return;
}



/*prt_wl print floating wavelength in Angstroms, in output format */
void prt_wl( FILE *ioOUT , realnum wl )
{
	DEBUG_ENTRY( "prt_wl()" );

	string chString;
	sprt_wl( chString , wl );

	fprintf(ioOUT, "%.*s", LineSave.wl_length, chString.c_str() );
	return;
}

/* write wavelength to string */
void sprt_wl( string& chString , realnum wl )
{
	string chUnits;

	DEBUG_ENTRY( "sprt_wl()" );

	/* print in A unless > 1e4, then use microns */
	if( wl > 1e8 )
	{
		/* centimeters */
		chUnits = "c";
		wl /= 1e8;
	}
	else if( wl > 1e4 )
	{
		/* microns */
		chUnits = "m";
		wl /= 1e4;
	}
	else if( wl == 0. )
	{
		chUnits = " ";
	}
	else
	{
		/* Angstroms units */
		chUnits = "A";
	}

	ostringstream oss;
	/* want total of LineSave.sig_figs sig figs */
	if( wl==0. )
	{
		oss << setw(LineSave.wl_length-1) << 0;
	}
	else
	{
		int n = LineSave.sig_figs - 1 - (int)log10(wl);
		if( n > 0 )
		{
			oss << fixed << setw(LineSave.sig_figs+1) << setprecision(n) << wl;
		}
		else if (wl < (realnum)INT_MAX)
		{
			oss << setw(LineSave.sig_figs+1) << (int)wl;
		}
		else
		{
			oss << setw(LineSave.sig_figs+1) << "*";
		}
	}
	chString = oss.str() + chUnits;
	return;
}

/*prt_LineLabels save all labels and wavelengths for emission line array */
void prt_LineLabels(
	/* io file handle */
	FILE * ioOUT ,
	/* print all if true, if false then do not print parts of 
	 * transferred lines */
	bool lgPrintAll,
       	/* print index of line in line stack */
	bool lgPrintIndex )
{
	long int i;

	DEBUG_ENTRY( "prt_LineLabels()" );

	for( i=0; i < LineSave.nsum; i++ )
	{
		if( LineSave.lines[i].isSeparator() )
		{
			fprintf(ioOUT,"####\t%s",LineSave.chHoldComments[(int)LineSave.lines[i].wavelength()].c_str()); 
		}
		else
		{
			if( !lgPrintAll &&
				 ( LineSave.lines[i].isInward() ||
				   LineSave.lines[i].isCollisional() ||
				   LineSave.lines[i].isPump() ||
				   LineSave.lines[i].isHeat() )
				)
				/* option to do not print lots of redundant labels 
				 * lgPrintAll is false by default set true with LONG option
				 * on save line labels command */
				continue;

			/* this format chosen to be identical to that used by final */
			if( lgPrintIndex )
				fprintf( ioOUT, "%li\t", i );

			fprintf( ioOUT, "%s\t", LineSave.lines[i].label().c_str() );

			/* skip over leading spaces - a formatting problem */
			long int j = 0;
			string comment = LineSave.lines[i].chComment();
			while( comment[j] != '\0' && comment[j] == ' ' )
				++j;
			/* comment entered when line intensity generated and other useful information */
			fprintf(ioOUT, "\t# type: %c, ", LineSave.lines[i].chSumTyp() );
			auto tr = LineSave.lines[i].getTransition();
			if( tr.associated() )
			{
				fprintf(ioOUT, "index=%d, %d ", tr.ipLo()+1, tr.ipHi()+1);
				fprintf(ioOUT, "Elow=%.7g   ", tr.Lo()->energy().WN());
			}
			fprintf(ioOUT, "%s" , comment.substr(j).c_str());
		}
		fprintf( ioOUT, "\n" );
	}
	return;
}

/* prt_line_err produce an error message containing the line label and wavelength,
 * 		followed, if given, by the wavelength of the closest line of the same label */
void prt_line_err( FILE *ioOUT, const LineID& line )
{
	fprintf( ioOUT, "with label (between quotes) \"%s\" and wavelength ", line.chLabel.c_str() );
	prt_wl ( ioOUT, line.wave );
	fprintf( ioOUT, ".\n" );
	return;
}
void prt_line_err( FILE *ioOUT, const string& label, realnum wvlng )
{
	prt_line_err( ioOUT, LineID(label, wvlng) );
}


/* prt_line_inlist print line suitable for output list, label not enclosed in quotation marks */
void prt_line_inlist ( FILE *ioOUT, const char *label, realnum wvlng )
{
	fprintf( ioOUT, "%-*s\t", NCHLAB-1, label );
	prt_wl ( ioOUT, wvlng );

	return;
}


void t_prt_matrix::zero()
{
	DEBUG_ENTRY( "t_prt_matrix::zero()" );

	species = "";
	speciesLevels = "";
}

void t_prt_matrix::setSpecies( const string &sspec )
{
	DEBUG_ENTRY( "t_prt_matrix::setSpecies()" );

	size_t lbrac = sspec.find( "[" );
	size_t rbrac = sspec.find( "]" );

	if( ( lbrac < string::npos && rbrac == string::npos ) ||
		( lbrac == string::npos && rbrac < string::npos ) )
	{
		fprintf( ioQQQ, "PROBLEM: Unbalanced brackets '[]' in '%s'\n",
				sspec.c_str() );
		cdEXIT( EXIT_FAILURE );
	}

	speciesLevels = sspec;
	species = sspec.substr( 0, lbrac );

	if( lbrac == string::npos )
		speciesLevels += "[:]";
}

void t_prt_matrix::resolveLevels()
{
	DEBUG_ENTRY( "t_prt_matrix::`resolveLevels()" );

	if( speciesLevels.length() == 0 )
		return;

	getLevelsGeneric( speciesLevels, true, speciesLevelList );
	lgLevelsResolved = true;
}

void t_prt_matrix::prtRates( const long numLevels,
				const multi_arr<double,2,C_TYPE> &matrix,
				valarray<double> &b )
{
	DEBUG_ENTRY( "t_prt_matrix::prtRates()" );

	if( not lgLevelsResolved )
	{
		resolveLevels();
	}

	if( speciesLevelList.size() == 0 )
		return;

	fprintf( ioQQQ, "'%s' lvl / creation /=>rates", species.c_str() );
	for( vector<long>::iterator ipLo = speciesLevelList.begin();
		 ipLo != speciesLevelList.end(); ++ipLo )
	{
		if( *ipLo >= numLevels )
			continue;
		if( ipLo == speciesLevelList.begin() )
			fprintf( ioQQQ, "\t%3ld", *ipLo+1 );
		else
			fprintf( ioQQQ, "\t%11ld", *ipLo+1 );
	}
	fprintf( ioQQQ, "\n" );

	for( vector<long>::iterator ipLo = speciesLevelList.begin();
		 ipLo != speciesLevelList.end(); ++ipLo )
	{
		if( *ipLo >= numLevels )
			continue;
		fprintf( ioQQQ, "%3ld\t %.4e", *ipLo+1, b[ *ipLo ] );
		for( vector<long>::iterator ipHi = speciesLevelList.begin();
			 ipHi != speciesLevelList.end(); ++ipHi )
		{
			if( *ipHi >= numLevels )
				continue;
			fprintf( ioQQQ, "\t%11.4e", matrix[ *ipLo ][ *ipHi ] );
		}
		fprintf( ioQQQ, "\n" );
	}
}
