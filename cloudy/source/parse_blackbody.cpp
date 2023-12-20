/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseBlackbody parse parameters off black body command */
#include "cddefines.h"
#include "optimize.h"
#include "input.h" 
#include "rfield.h"
#include "radius.h"
#include "parser.h"

void ParseBlackbody(
  /* input command line, already changed to caps */
  Parser &p)
{
	bool 
		lgIntensitySet=false;
	double a, 
	  dil, 
	  rlogl;

	string chParamType;
	int nParam = 0;

	DEBUG_ENTRY( "ParseBlackbody()" );

	set_NaN( rlogl );

	/* type is blackbody */
	strcpy( rfield.chSpType[rfield.nShape], "BLACK" );
	strcpy( rfield.chSpNorm[p.m_nqh], "LUMI" );

	/* these two are not used for this continuum shape */
	rfield.cutoff[rfield.nShape][0] = 0.;
	rfield.cutoff[rfield.nShape][1] = 0.;

	/* get the blackbody temperature */
	rfield.slope[rfield.nShape] = p.FFmtRead();
	if( p.lgEOL() )
		p.NoNumb("blackbody temperature");

	/* this is the temperature - make sure its linear in the end
	 * there are two keys, LINEAR and LOG, that could be here,
	 * else choose which is here by which side of 10 */
	if( (rfield.slope[rfield.nShape] <= 10. && !p.nMatch("LINE")) || 
		p.nMatch(" LOG") )
	{
		/* log option */
		if( rfield.slope[rfield.nShape]>log10(BIGFLOAT) )
		{
			fprintf(ioQQQ,"PROBLEM The specified log of the temperature, %.3e, is too large.\nSorry.\n",
					rfield.slope[rfield.nShape]);
			cdEXIT(EXIT_FAILURE);
		}
		rfield.slope[rfield.nShape] = exp10(rfield.slope[rfield.nShape]);
	}

	/* check that temp is not too low - could happen if log misused */
	if( rfield.slope[rfield.nShape] < 1e4 )
	{
		fprintf( ioQQQ, " Is T(star)=%10.2e correct???\n", 
		  rfield.slope[rfield.nShape] );
	}

	/* now check that temp not too low - would peak below low
	 * energy limit of the code
	 * factor is temperature of 1 Ryd, egamry is high-energy limit of code */
	if( rfield.slope[rfield.nShape]/TE1RYD < rfield.emm() )
	{
		fprintf( ioQQQ, " This temperature is very low - the blackbody will have significant flux low the low energy limit of the code, presently %10.2e Ryd.\n", 
		  rfield.emm() );
		fprintf( ioQQQ, " Was this intended?\n" );
	}

	/* now check that temp not too high - would extend beyond high
	 * energy limit of the code
	 * factor is temperature of 1 Ryd, egamry is high-energy limit of code */
	if( rfield.slope[rfield.nShape]/TE1RYD*2. > rfield.egamry() )
	{
		fprintf( ioQQQ, " This temperature is very high - the blackbody will have significant flux above the high-energy limit of the code,%10.2e Ryd.\n", 
		  rfield.egamry() );
		fprintf( ioQQQ, " Was this intended?\n" );
	}

	/* also possible to input log(total luminosity)=real log(l) */
	a = p.FFmtRead();

	/* there was not a second number on the line; check if LTE or STE */
	if( p.nMatch(" LTE") || p.nMatch("LTE ") ||
	    p.nMatch(" STE") || p.nMatch("STE ") )
	{
		/* set energy density to the STE - strict thermodynamic equilibrium - value */
		chParamType = "STE";
		nParam = 1;

		if( !p.lgEOL() )
		{
			fprintf(ioQQQ,"PROBLEM the luminosity was specified on "
				"the BLACKBODY K STE command.\n");
			fprintf(ioQQQ,"Do not specify the luminosity since STE does this.\n");
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* use blackbody relations to get intensity from temperature */
		rlogl = log10(4.*STEFAN_BOLTZ) + 4.*log10(rfield.slope[rfield.nShape]);

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		lgIntensitySet = true;

		if (p.nMatch(" STE") || p.nMatch("STE "))
			rfield.Illumination[p.m_nqh] = Illuminate::ISOTROPIC;
	}

	/* a second number was entered, what was it? */
	else if( p.nMatch("LUMI") )
	{
		chParamType = "LUMINOSITY";
		nParam = 2;
		rlogl = a;
		strcpy( rfield.chRSpec[p.m_nqh], "4 PI" );
		if( p.lgEOL() )
			p.NoNumb("luminosity" );
		lgIntensitySet = true;
	}

	else if( p.nMatch("RADI") )
	{
		chParamType = "RADIUS";
		nParam = 2;
		/* radius was entered, convert to total luminosity */
		rlogl = -3.147238 + 2.*a + 4.*log10(rfield.slope[rfield.nShape]);
		strcpy( rfield.chRSpec[p.m_nqh], "4 PI" );
		if( p.lgEOL() )
			p.NoNumb("radius" );
		lgIntensitySet = true;
	}

	else if( p.nMatch("DENS") )
	{
		chParamType = "ENERGY DENSITY";
		nParam = 2;
		/* number was temperature to deduce energy density
		 * number is linear if greater than 10, or if LINEAR appears on line
		 * want number to be log of temperature at end of this */
		if( !p.nMatch(" LOG") && (p.nMatch("LINE") || a > 10.) )
		{
			a = log10(a);
		}
		rlogl = log10(4.*STEFAN_BOLTZ) + 4.*a;
		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		if( p.lgEOL() )
			p.NoNumb("energy density");
		lgIntensitySet = true;
	}

	else if( p.nMatch("DILU") )
	{
		chParamType = "DILUTION FACTOR";
		nParam = 2;
		/* number is dilution factor, if negative then its log */
		if( a <= 0. )
			dil = a;
		else
			dil = log10(a);

		if( dil > 0. )
			fprintf( ioQQQ, "PROBLEM Is the dilution factor > 1 on this "
			"blackbody command physical?\n" );

		/* intensity from black body relations and temperature */
		rlogl = log10(4.*STEFAN_BOLTZ) + 4.*log10(rfield.slope[rfield.nShape]);

		/* add on dilution factor */
		rlogl += dil;

		strcpy( rfield.chRSpec[p.m_nqh], "SQCM" );
		if( p.lgEOL() )
			p.NoNumb("dilution factor" );
		lgIntensitySet = true;
	}

	else if( p.nMatch("DISK") )
	{
		if( p.lgEOL() )
			p.NoNumb("disk Te" );

		chParamType = "DISK";
		nParam = 2;
		
		rfield.cutoff[rfield.nShape][0] = a;
		/* this is the temperature - make sure its linear in the end
		 * there are two keys, LINEAR and LOG, that could be here,
		 * else choose which is here by which side of 10 */
		if( (rfield.cutoff[rfield.nShape][0] <= 10. && !p.nMatch("LINE")) ||
			p.nMatch(" LOG") )
		{
			/* log option */
			rfield.cutoff[rfield.nShape][0] = exp10(rfield.cutoff[rfield.nShape][0]);
		}
		a = log10( rfield.cutoff[rfield.nShape][0] );

		strcpy( rfield.chSpType[rfield.nShape], "DISKB" );
		lgIntensitySet = false;
	}

	if( lgIntensitySet )
	{
		/* a luminosity option was specified
		 * check that stack of shape and luminosity specifications
		 * is parallel, stop if not - this happens is background comes
		 * BETWEEN another set of shape and luminosity commands */
		if( rfield.nShape != p.m_nqh )
		{
			fprintf( ioQQQ, " This command has come between a previous ordered pair of continuum shape and luminosity commands.\n Reorder the commands to complete each continuum specification before starting another.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		rfield.range[p.m_nqh][0] = rfield.emm();
		rfield.range[p.m_nqh][1] = rfield.egamry();
		rfield.totpow[p.m_nqh] = rlogl;
		++p.m_nqh;
	}
	/* vary option */
	if( optimize.lgVarOn )
	{
		/* this test no option on blackbody command */
		if( chParamType.length() == 0 )
		{
			/* no luminosity options on vary */
			optimize.nvarxt[optimize.nparm] = 1;
			strcpy( optimize.chVarFmt[optimize.nparm], "BLACKbody= %f LOG" );
		}
		else
		{
			char chHold[100];
			/* there was an option - honor it */
			if( nParam==1 )
			{
				optimize.nvarxt[optimize.nparm] = 1;
				strcpy( chHold , "BLACKbody= %f LOG ");
				strcat( chHold , chParamType.c_str() );
			}
			else if( nParam==2 )
			{
				optimize.nvarxt[optimize.nparm] = 2;
				optimize.vparm[1][optimize.nparm] = (realnum)a;
				strcpy( chHold , "BLACKbody= %f LOG %f ");
				strcat( chHold , chParamType.c_str() );
			}
			else
				TotalInsanity();
			strcpy( optimize.chVarFmt[optimize.nparm], chHold );
		}

		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		/* log of temp stored here  */
		optimize.vparm[0][optimize.nparm] = (realnum)log10(rfield.slope[rfield.nShape]);
		/* the increment in the first steps away from the original value */
		optimize.vincr[optimize.nparm] = 0.5f;
		++optimize.nparm;
	}

	/* increment SED indices */
	++rfield.nShape;
	if( rfield.nShape >= LIMSPC )
	{
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return;
}
