/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*save_line parse save lines command, or actually do the save lines output */
/*Save_Line_RT parse the save line rt command - read in a set of lines */
#include "cddefines.h"
#include "cddrive.h"
#include "radius.h"
#include "opacity.h"
#include "phycon.h"
#include "dense.h"
#include "lines.h"
#include "h2.h"
#include "prt.h"
#include "iso.h"
#include "parser.h"
#include "save.h"

/* implement the save line xxx command.  cumulative, structure, and
 * emissivity all use same code base and variables, so only one can be used
 * at present */

class SaveLineList
{
public:
	vector<LineID> lineids;
	vector<long> ipLine;
	bool lgRelativeIntensity, lgMustGetLines, lgMustPrintFirstTime;
	SaveLineList() : lgRelativeIntensity(false), lgMustGetLines(true), lgMustPrintFirstTime(true) {}
};

static shared_ptr<SaveLineList> linelist[LIMPUN];

void parse_save_line(Parser &p, 
					 /* true, return rel intensity, false, log of luminosity or intensity I */
					 bool lgLog3,
					 ostringstream& chHeader,
					 long ipPun)
{
	DEBUG_ENTRY( "parse_save_line()" );

	/* very first time this routine is called, chDo is "READ" and we read
	 * in lines from the input stream.  The line labels and wavelengths
	 * are store locally, and output in later calls to this routine
	 * following is flag saying whether to do relative intensity or
	 * absolute emissivity */
	linelist[ipPun] = shared_ptr<SaveLineList>(new SaveLineList);
	linelist[ipPun]->lgRelativeIntensity = lgLog3;

	ParseLineList(p, linelist[ipPun]->lineids);
	linelist[ipPun]->ipLine.resize(linelist[ipPun]->lineids.size());

	chHeader << "#depth";
	for( size_t i=0; i < linelist[ipPun]->lineids.size(); i++ )
	{
		string chTemp;
		sprt_wl( chTemp, linelist[ipPun]->lineids[i].wave );
		chHeader << "\t" << linelist[ipPun]->lineids[i].chLabel << " " << chTemp;
	}
	chHeader << endl;
}

void save_line(FILE * ioPUN, /* the file we will write to */
			   const char *chDo,
			   // intrinsic or emergent line emission?
			   bool lgEmergent,
			   long ipPun)
{
	long int i;

	DEBUG_ENTRY( "save_line()" );

	long nLinesNow = linelist[ipPun]->lineids.size();
	vector<double> a(nLinesNow);

	bool lgBadLine = false;
	if( nzone <= 1 && linelist[ipPun]->lgMustGetLines )
	{
		for( i=0; i < nLinesNow; i++ )
		{
			linelist[ipPun]->ipLine[i] = LineSave.findline(linelist[ipPun]->lineids[i]);
			if( linelist[ipPun]->ipLine[i] <= 0 )
			{
				// missed line - ignore if H2 line since large model may be off
				if( !h2.lgEnabled && linelist[ipPun]->lineids[i].chLabel == "H2" )
				{
					if( linelist[ipPun]->lgMustPrintFirstTime )
					{
						/* it's an H2 line and H2 is not being done - ignore it */
						fprintf( ioQQQ,"\nPROBLEM Did not find an H2 line, the large model is not "
									"included, so I will ignore it.  Log intensity set to -30.\n" );
						fprintf( ioQQQ,"I will totally ignore any future missed H2 lines\n\n");
						linelist[ipPun]->lgMustPrintFirstTime = false;
					}
					/* flag saying to ignore this line */
					linelist[ipPun]->ipLine[i] = -2;
				}
				else
				{
					fprintf( ioQQQ, " save_line could not find line ");
					prt_line_err( ioQQQ, linelist[ipPun]->lineids[i] );
					linelist[ipPun]->ipLine[i] = -1;
					lgBadLine = true;
				}
			}
		}
		linelist[ipPun]->lgMustGetLines = false;
		if( lgBadLine )
		{
			cdEXIT(EXIT_FAILURE);
		}
	}

	if( strcmp(chDo,"PUNS") == 0 )
	{
		/* save lines emissivity command */
		/* save lines structure command */
		for( i=0; i < nLinesNow; i++ )
		{
			/* this version of cdEmis uses index, does not search, do not call if line could not be found */
			if( linelist[ipPun]->ipLine[i] > 0 )
				cdEmis_ip(linelist[ipPun]->ipLine[i],&a[i],lgEmergent);
			else
				a[i] = 1e-30f;
		}
	}

	else if( strcmp(chDo,"PUNC") == 0 )
	{
		/* save lines cumulative command */		
		for( i=0; i < nLinesNow; i++ )
		{
			if( linelist[ipPun]->ipLine[i]<=0 )
			{
				a[i] = 0.;
			}
			else
			{
				double absint, relint;
				cdLine_ip(linelist[ipPun]->ipLine[i],&relint,&absint,lgEmergent);
				if( linelist[ipPun]->lgRelativeIntensity )
					/* relative intensity case */
					a[i] = relint;
				else
					/* emissivity or luminosity case */
					a[i] = absint;
			}
		}		
	}
	else if( strcmp(chDo,"PUNO") == 0 )
	{
		/* save lines optical depth some command */		
		for( i=0; i < nLinesNow; i++ )
		{
			if( linelist[ipPun]->ipLine[i]<=0 )
			{
				a[i] = 0.;
			}
			else
			{
				TransitionProxy tr = LineSave.lines[linelist[ipPun]->ipLine[i]].getTransition();
				if (tr.associated())
					a[i] = tr.Emis().TauIn()*SQRTPI;
				else
					a[i] = 0.;
			}
		}		
	}
	else
	{
		fprintf( ioQQQ, 
			" unrecognized key for save_line=%4.4s\n", 
		  chDo );
		cdEXIT(EXIT_FAILURE);
	}

	fprintf( ioPUN, "%.5e", radius.depth_mid_zone );
	
	for( i=0; i < nLinesNow; i++ )
	{
		fprintf( ioPUN, "\t%.4e", a[i] );
	}
	fprintf( ioPUN, "\n" );

	return;
}

static const int LIMLINE = 10;
static long int line_RT_type[LIMLINE] = 
	{LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, 
	 LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN },
	line_RT_ipISO[LIMLINE] =  
	{LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, 
	 LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN },
	line_RT_nelem[LIMLINE] =  
	{LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, 
	 LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN },
	line_RT_ipHi[LIMLINE] =  
	{LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, 
	 LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN },
	line_RT_ipLo[LIMLINE] = 
	{LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, 
	 LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN };
static bool lgMustPrintHeader = true;
static long int nLine = -1;

/*Save_Line_RT parse the save line rt command - read in a set of lines */
void Parse_Save_Line_RT(Parser &p)
{
	/* save line rt parameters */
	DEBUG_ENTRY( "Parse_Save_Line_RT()" );

	/* very first time this routine is called, chDo is "READ" and we read
	 * in lines from the input stream.  The line labels and wavelengths
	 * are store locally, and output in later calls to this routine */
	
	/* say that we must print the header */
	lgMustPrintHeader = true;
	
	/* get the next line, and check for eof */
	nLine = 0;
	p.getline();
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, 
					" Hit EOF while reading line list; use END to end list.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	do
	{
		if( nLine >= LIMLINE )
		{
			fprintf(ioQQQ," PUNCH RT has too many lines - increase LIMLINE in save_line.cpp\n");
			cdEXIT(EXIT_FAILURE);
		}
		
		/* right now it just does lines in the iso sequences */
		line_RT_type[nLine] = (long int)p.FFmtRead();
		line_RT_ipISO[nLine] = (long int)p.FFmtRead();
		line_RT_nelem[nLine] = (long int)p.FFmtRead();
		line_RT_ipHi[nLine] = (long int)p.FFmtRead();
		line_RT_ipLo[nLine] = (long int)p.FFmtRead();
		
		if( p.lgEOL() )
		{
			fprintf( ioQQQ, " there must be five numbers on this line\n");
			p.PrintLine(ioQQQ);
			cdEXIT(EXIT_FAILURE);
		}
		
		/* now increment number of lines */
		++nLine;
		
		/* now get next line until we hit eof or the word END */
		p.getline();
	} while( !p.m_lgEOF && !p.nMatch( "END") );
	if( p.m_lgEOF )
	{
		fprintf( ioQQQ, " Save_Line_RT hit end of file looking for END of RT lines\n");
		p.PrintLine(ioQQQ);
		cdEXIT(EXIT_FAILURE);
	}
}

void Save_Line_RT( 
	FILE * ioPUN )
{
	/* save line rt parameters */

	DEBUG_ENTRY( "Save_Line_RT()" );


	static char chLabel[LIMLINE][30];
	long int n;
	if( lgMustPrintHeader )
	{
		fprintf( ioPUN, "Line\tP(con,inc)\tAul\tgl\tgu\n");
		for( n=0; n<nLine; ++n )
		{
			TransitionProxy tr = iso_sp[line_RT_ipISO[n]][line_RT_nelem[n]].trans(line_RT_ipHi[n],line_RT_ipLo[n]);
			/* print info for header of file, line id and pump rate */
			sprintf( chLabel[n], "%s ", 
					chLineLbl(tr).c_str() );
			fprintf( ioPUN, "%s ", chLabel[n] );
			fprintf( ioPUN, "%.4e ",
					tr.Emis().pump());
			fprintf( ioPUN, "%.4e ",
					tr.Emis().Aul());
			fprintf( ioPUN, "%.0f ",
					(*tr.Lo()).g());
			fprintf( ioPUN, "%.0f ",
					(*tr.Hi()).g());
			fprintf( ioPUN, "\n");
			
			if( line_RT_type[n]!=0. )
			{
				/* for now code only exists for H He like iso - assert this */
				fprintf( ioQQQ, 
							" Save_Line_RT only H, He like allowed for now\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
		fprintf( ioPUN, "Line\tTauIn\tPopLo\tPopHi\tCul\tk(line)\tk(con,abs)\tk(con,scat)\n");
		lgMustPrintHeader = false;
	}
	
	fprintf(ioPUN, "RADIUS\t%e\tDEPTH\t%e\tTe\t%e\tNe\t%e\n",
			  radius.Radius_mid_zone, 
			  radius.depth_mid_zone, 
			  phycon.te, 
			  dense.eden );
	for( n=0; n<nLine; ++n )
	{
		TransitionProxy tr = iso_sp[line_RT_ipISO[n]][line_RT_nelem[n]].trans(line_RT_ipHi[n],line_RT_ipLo[n]);

		/* index for line within continuum array */
		long int ipCont = tr.ipCont();
		fprintf( ioPUN, "%s ", chLabel[n] );
		fprintf( ioPUN, "\t%e\t%e\t%e",
					tr.Emis().TauIn(), 
					(*tr.Lo()).Pop(),
					(*tr.Hi()).Pop()
			);
		fprintf( ioPUN, "\t%e",
					tr.Coll().ColUL( colliders ) / dense.EdenHCorr
			);
		
		fprintf( ioPUN, "\t%e\t%e\t%e\n",
					tr.Emis().PopOpc(),
					opac.opacity_abs[ipCont-1], 
					opac.opacity_sct[ipCont-1]
			);
	}
}
