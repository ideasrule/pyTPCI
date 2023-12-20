/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*prtmet print all line optical depths at end of iteration */
#include "cddefines.h"
#include "taulines.h"
#include "h2.h"
#include "iso.h"
#include "dense.h"
#include "prt.h"
#include "trace.h"

STATIC void prme(
  const bool lgReset,
  const TransitionProxy &t);

STATIC void prt_header_cols();

/*prtmet print all line optical depths at end of iteration */
void prtmet(void)
{
	DEBUG_ENTRY( "prtmet()" );

	/* default is to not print optical depths, turn on with
	 * print optical depths on command */
	if( prt.lgPrtTau  || (trace.lgTrace && trace.lgOptcBug) )
	{
		fprintf( ioQQQ, "\n\n                                                 Mean Line Optical Depths\n");

		prt_header_cols();

		/* iso sequences */
		for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( long nelem=ipISO; nelem < LIMELM; nelem++ )
			{
				if( dense.lgElmtOn[nelem] )
				{
					if( (*iso_sp[ipISO][nelem].trans(1,0).Lo()).ColDen() <= 0. )
						continue;

					/* print Lyman, Balmer, Paschen, etc sequence optical depths */	
					for( long ipLo=0; ipLo < iso_sp[ipISO][nelem].numLevels_local-1; ipLo++ )
					{
						for( long ipHi=ipLo+1; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
						{
							prme(false,iso_sp[ipISO][nelem].trans(ipHi,ipLo));
						}
					}
				}
			}
		}

		/* print main lines optical depths */
		for( long i=0; i < nWindLine; i++ )
		{
			if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
			{
				prme(false,TauLine2[i]);
			}
		}

		for( size_t i=0; i < UTALines.size(); i++ )
		{
			prme(false,UTALines[i]);
		}

		/* print H2 line optical depths */
		for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		{
			for( TransitionList::iterator tr = (*diatom)->trans.begin();
				tr != (*diatom)->trans.end(); ++tr )
			{
				if( (*tr).ipCont() <= 0 )
					continue;
				prme( false, *tr );
			}
		}

		for( size_t i=0; i < HFLines.size(); i++ )
		{
			prme(false,HFLines[i]);
		}

		/* data base lines */
		for (int ipSpecies=0; ipSpecies < nSpecies; ++ipSpecies)
		{
			for( EmissionList::iterator em=dBaseTrans[ipSpecies].Emis().begin();
				  em != dBaseTrans[ipSpecies].Emis().end(); ++em)
			{
				prme(false,(*em).Tran());
			}
		}

		fprintf( ioQQQ, "\n");
	}
	return;
}

STATIC void prt_header_cols()
{
	string columns = "Species Wavelength  Total     Single  ";
	size_t nchars = columns.length();
	size_t header_len = 0;

	while( header_len < NCOLMAX )
	{
		fprintf( ioQQQ, "%s", columns.c_str() );
		header_len += nchars;
		if( header_len > NCOLMAX )
			continue;
		fprintf( ioQQQ, "%s", prt_linecol.col_gap.c_str() );
		header_len += prt_linecol.col_gap_len;
	}

	fprintf( ioQQQ, "\n");
}


/* prme - print line optical depth */
STATIC void prme(
  const bool lgReset,
  const TransitionProxy &t)
{
	static long int n ;

	DEBUG_ENTRY( "prme()" );

	if( lgReset )
		n = 0;

	if( t.ipCont() <= 0 )
	{
		/* line is not transferred */
		return;
	}

	if( (*t.Lo()).ColDen() <= 0. )
		return;

	/* print optical depth if greater than lower limit, or significantly negative
	 * PrtTauFnt is threshold for printing it
	 * */
	if( t.Emis().TauIn()*SQRTPI > prt.PrtTauFnt || t.Emis().TauIn()*SQRTPI < -1e-5 )
	{
		// throw CR after printing NCOLMAX characters
		const int NOPCMAX = 10;
		string label = t.chLabel();
		long len = 2*NOPCMAX + label.size();
		
		// Reset n if would overflow
		if(n+len > NCOLMAX)
		{
			n = 0;
			fprintf( ioQQQ, "\n");
		}
		fprintf( ioQQQ, "%s", label.c_str() );
		/*>> chng 12 jul 25, print mean optical depths, rather than line center */
		fprintf( ioQQQ, "%*.2e", NOPCMAX, t.Emis().TauIn()*SQRTPI );

		fprintf( ioQQQ, "%*.2e", NOPCMAX, t.Emis().TauInSpecific()*SQRTPI );

		// Add wide gap between columns 
		if( n + 2*len + prt_linecol.col_gap_len < NCOLMAX )
		{
			fprintf( ioQQQ, "%s", prt_linecol.col_gap.c_str() );
			len += prt_linecol.col_gap_len;
		}
		else
		{
			n = NCOLMAX;
		}
		n += len;
	}

	return;
}
