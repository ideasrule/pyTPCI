/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*save_opacity save total opacity in any element, save opacity command */
#include "cddefines.h"
#include "iso.h"
#include "rfield.h"
#include "ipoint.h"
#include "continuum.h"
#include "opacity.h"
#include "dense.h"
#include "yield.h"
#include "atmdat_adfa.h"
#include "abund.h"
#include "heavy.h"
#include "elementnames.h"
#include "save.h"
#include "freebound.h"
#include "species.h"

/* print final information about where opacity files are */
STATIC void prtPunOpacSummary(void);

void save_opacity(FILE * ioPUN, 
  long int ipPun)
{
	/* this will be the file name for opacity output */
	char chFileName[FILENAME_PATH_LENGTH_2];	

	/* this is its pointer */
	FILE *ioFILENAME;

	char chNumbers[31][3] = {"0","1","2","3","4","5","6","7","8","9",
		"10","11","12","13","14","15","16","17","18","19",
		"20","21","22","23","24","25","26","27","28","29","30"};

	long int i, 
	  ilow, 
	  ion, 
	  ipS, 
	  j, 
	  nelem;

	double ener, 
	  ener3;

	DEBUG_ENTRY( "save_opacity()" );

	/* this flag says to redo the static opacities */
	opac.lgRedoStatic = true;

	/* save total opacity in any element, save opacity command
	 * ioPUN is ioPUN unit number, ipPun is pointer within stack of punches */

	if( strcmp(save.chOpcTyp[ipPun],"TOTL") == 0 )
	/* total opacity */
	{
		for( j=0; j < rfield.nflux; j++ )
		{
			fprintf( ioPUN, "%12.4e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%4.4s\n", 
			  AnuUnit(rfield.anu(j)),
			  opac.opacity_abs[j]+opac.OpacStatic[j] + opac.opacity_sct[j], 
			  opac.opacity_abs[j]+opac.OpacStatic[j], 
			  opac.opacity_sct[j], 
			  opac.opacity_sct[j]/MAX2(1e-37,opac.opacity_sct[j]+opac.opacity_abs[j]+opac.OpacStatic[j]), 
			  rfield.chContLabel[j].c_str() );
		}

		fprintf( ioPUN, "\n" );
	}

	else if( strcmp(save.chOpcTyp[ipPun],"BREM") == 0 )
	/* bremsstrahlung opacity */
	{
		for( j=0; j < rfield.nflux; j++ )
		{
			fprintf( ioPUN, "%.5e\t%.3e\t%.3e\n", 
				 AnuUnit(rfield.anu(j)),
				 opac.FreeFreeOpacity[j],
				 pow2(dense.eden)*opac.eeFreeFreeOpacity[j] );
		}
	}

	/* subshell photo cross sections */
	else if( strcmp(save.chOpcTyp[ipPun],"SHEL") == 0 )
	{
		nelem = (long)save.punarg[ipPun][0];
		ion = (long)save.punarg[ipPun][1];
		ipS = (long)save.punarg[ipPun][2];

		if( !dense.lgElmtOn[nelem] )
		{
			fprintf( ioQQQ, "SAVE OPACITY SHELL: requested element is not active. Bailing out.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		for( i=opac.ipElement[nelem-1][ion-1][ipS-1][0]-1; 
		  i < opac.ipElement[nelem-1][ion-1][ipS-1][1]; i++ )
		{
			fprintf( ioPUN, 
				"%11.3e %11.3e\n", rfield.anu(i), 
			  opac.OpacStack[i-opac.ipElement[nelem-1][ion-1][ipS-1][0]+
			  opac.ipElement[nelem-1][ion-1][ipS-1][2]] );
		}
	}

	else if( strcmp(save.chOpcTyp[ipPun],"FINE") == 0 )
	{
		/* the fine opacity array */
		realnum sum;
		long int k, nu_hi , nskip;
		if( save.punarg[ipPun][0] > 0. )
		{
			/* possible lower bounds to energy range - will be zero if not set */
			j = ipFineCont( save.punarg[ipPun][0] );
		}
		else
		{
			j = 0;
		}

		/* upper limit */
		if( save.punarg[ipPun][1]> 0. )
		{
			nu_hi = ipFineCont( save.punarg[ipPun][1]);
		}
		else
		{
			nu_hi = rfield.nfine;
		}

		nskip = (long)abs(save.punarg[ipPun][2]);
		nskip = MAX2( 1, nskip );

		for( i=j; i<nu_hi; i+=nskip )
		{
			sum = 0;
			for( k=0; k<nskip; ++k )
			{
				sum += rfield.fine_opac_zone[i+k];
			}
			sum /= nskip;
			// ALL key -> report all, even zero
			if( sum > 0.  ||  save.punarg[ipPun][2]<0)
				fprintf(ioPUN,"%.7e\t%.7e\n",
					AnuUnit( rfield.fine_anu[i+k/2] ),	sum );
		}
	}

	/* figure for hazy */
	else if( strcmp(save.chOpcTyp[ipPun],"FIGU") == 0 )
	{
		nelem = 0;
		for( i=iso_sp[ipH_LIKE][ipHYDROGEN].fb[0].ipIsoLevNIonCon-1; i < (iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon - 1); i++ )
		{
			ener = rfield.anu(i)*0.01356;
			ener3 = 1e24*POW3(ener);
			fprintf( ioPUN, 
				"%12.4e%12.4e%12.4e%12.4e%12.4e\n", 
			  rfield.anu(i), ener, 
			  opac.OpacStack[i-iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon+ iso_sp[ipH_LIKE][nelem].fb[ipH1s].ipOpac]*ener3, 
			  0., 
			  (opac.opacity_abs[i]+opac.OpacStatic[i])/ dense.gas_phase[ipHYDROGEN]*ener3 );
		}

		if( dense.lgElmtOn[ipHELIUM] )
		{
			for( i=iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon-1; i < rfield.nflux_with_check; i++ )
			{
				ener = rfield.anu(i)*0.01356;
				ener3 = 1e24*POW3(ener);
				fprintf( ioPUN, 
						 "%12.4e%12.4e%12.4e%12.4e%12.4e\n", 
						 rfield.anu(i), 
						 ener, 
						 opac.OpacStack[i-iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon+ iso_sp[ipH_LIKE][nelem].fb[ipH1s].ipOpac]*ener3, 
						 opac.OpacStack[i-iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon+ opac.iophe1]*dense.gas_phase[ipHELIUM]/dense.gas_phase[ipHYDROGEN]*ener3, 
						 (opac.opacity_abs[i]+opac.OpacStatic[i])/dense.gas_phase[ipHYDROGEN]*ener3 );
			}
		}
	}

	/* photoionization data table for AGN */
	else if( strcmp(save.chOpcTyp[ipPun]," AGN") == 0 )
	{
		/* now do loop over temp, but add elements */
		for( nelem=ipLITHIUM; nelem<LIMELM; ++nelem )
		{
			/* this list of elements included in the AGN tables is defined in zeroabun.c */
			if( dense.lgElmtOn[nelem] && abund.lgAGN[nelem] )
			{
				for( ion=0; ion<=nelem; ++ion )
				{
					/* number of bound electrons */
					long nelec = nelem+1 - ion;
					fprintf( ioPUN, "%s ", makeChemical( nelem, nelec ).c_str() );
					/*fprintf(ioPUN,"\t%.2f\n", Heavy.Valence_IP_Ryd[nelem][ion] );*/

					/* now loop over all shells */
					for( long nshell=0; nshell < Heavy.nsShells[nelem][ion]; nshell++ )
					{
						/* shell designation */
						fprintf(ioPUN,"\t%s",Heavy.chShell[nshell] );

						/* ionization potential of shell */
						fprintf(ioPUN,"\t%.2f" ,
							t_ADfA::Inst().getEthresh(nshell+1,nelec,nelem+1)/EVRYD);

						/* set lower and upper limits to this range */
						long ipop = opac.ipElement[nelem][ion][nshell][2];
						fprintf(ioPUN,"\t%.2f",opac.OpacStack[ipop-1]/1e-18);
						for( i=0; i < t_yield::Inst().nelec_eject(nelem,ion,nshell); ++i )
						{
							fprintf(ioPUN,"\t%.2f",
								t_yield::Inst().elec_eject_frac(nelem,ion,nshell,i));
						}
						fprintf(ioPUN,"\n");
					}

				}
			}
		}
	}

	/* hydrogen */
	else if( strcmp(save.chOpcTyp[ipPun],"HYDR") == 0 )
	{
		nelem = ipHYDROGEN;
		/* zero out the opacity arrays */
		OpacityZero();

		OpacityAdd1SubshellInduc(iso_sp[ipH_LIKE][nelem].fb[ipH1s].ipOpac,iso_sp[ipH_LIKE][nelem].fb[ipH1s].ipIsoLevNIonCon,
		  rfield.nflux_with_check,1.,0.,'v');
		ilow = Heavy.ipHeavy[nelem][0];

		/* start filename for output */
		strcpy( chFileName , elementnames.chElementNameShort[0] );

		/* continue filename with ionization stage */
		strcat( chFileName , chNumbers[1] );

		/* end it with string ".opc", name will be of form HYDR.opc */
		strcat( chFileName , ".opc" );

		/* now try to open it for writing */
		ioFILENAME = open_data( chFileName, "w" );
		for( j=ilow; j <= rfield.nflux_with_check; j++ )
		{
			/* photon energy in rydbergs */
			PrintE93( ioFILENAME , rfield.anu(j-1)*EVRYD );
			fprintf( ioFILENAME , "\t");
			/* cross section in megabarns */
			PrintE93( ioFILENAME, (opac.opacity_abs[j-1]+opac.OpacStatic[j-1])/1e-18 );
			fprintf( ioFILENAME , "\n");
		}

		fclose( ioFILENAME );
		prtPunOpacSummary();
		cdEXIT(EXIT_SUCCESS);
	}

	/* helium */
	else if( strcmp(save.chOpcTyp[ipPun],"HELI") == 0 )
	{
		if( !dense.lgElmtOn[ipHELIUM] )
		{
			fprintf( ioQQQ, "SAVE ELEMENT OPACITY: requested element is not active. Bailing out.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* atomic helium first, HELI1.opc */
		nelem = ipHELIUM;
		OpacityZero();
		OpacityAdd1SubshellInduc(opac.iophe1,iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon,rfield.nflux,1.,0.,'v');
		ilow = Heavy.ipHeavy[nelem][0];

		/* start filename for output */
		strcpy( chFileName , elementnames.chElementNameShort[1] );

		/* continue filename with ionization stage */
		strcat( chFileName , chNumbers[1] );

		/* end it wil .opc, name will be of form HYDR.opc */
		strcat( chFileName , ".opc" );

		/* now try to open it for writing */
		ioFILENAME = open_data( chFileName, "w" );
		for( j=ilow; j <= rfield.nflux_with_check; j++ )
		{
			/* photon energy in rydbergs */
			PrintE93( ioFILENAME , rfield.anu(j-1)*EVRYD );
			fprintf( ioFILENAME , "\t");
			/* cross section in megabarns */
			PrintE93( ioFILENAME, (opac.opacity_abs[j-1]+opac.OpacStatic[j-1])/1e-18 );
			fprintf( ioFILENAME , "\n");
		}
		fclose( ioFILENAME );

		/* now do helium ion, HELI2.opc */
		OpacityZero();
		OpacityAdd1SubshellInduc(iso_sp[ipH_LIKE][1].fb[ipH1s].ipOpac,iso_sp[ipH_LIKE][1].fb[ipH1s].ipIsoLevNIonCon,rfield.nflux_with_check,1.,0.,'v');
		ilow = Heavy.ipHeavy[nelem][1];

		/* start filename for output */
		strcpy( chFileName , elementnames.chElementNameShort[1] );

		/* continue filename with ionization stage */
		strcat( chFileName , chNumbers[2] );

		/* end it wil .opc, name will be of form HYDR.opc */
		strcat( chFileName , ".opc" );

		/* now try to open it for writing */
		ioFILENAME = open_data( chFileName, "w" );
		for( j=ilow; j <= rfield.nflux_with_check; j++ )
		{
			/* photon energy in rydbergs */
			PrintE93( ioFILENAME , rfield.anu(j-1)*EVRYD );
			fprintf( ioFILENAME , "\t");
			/* cross section in megabarns */
			PrintE93( ioFILENAME, (opac.opacity_abs[j-1]+opac.OpacStatic[j-1])/1e-18 );
			fprintf( ioFILENAME , "\n");
		}

		prtPunOpacSummary();
		fclose( ioFILENAME );
		cdEXIT(EXIT_SUCCESS);
	}

	else
	{
		/* check for hydrogen through zinc, nelem is atomic number on the c scale */
		nelem = -1;
		i = 0;
		while( i < LIMELM )
		{
			if( strcmp(save.chOpcTyp[ipPun],elementnames.chElementNameShort[i]) == 0 )
			{
				nelem = i;
				break;
			}
			++i;
		}

		/* nelem is still negative if above loop fell through */
		if( nelem < 0 )
		{
			fprintf( ioQQQ, " Unidentified opacity key=%4.4s\n", 
			  save.chOpcTyp[ipPun] );
			cdEXIT(EXIT_FAILURE);
		}

		/* pc lint did not pick up the above logice an warned possible negative array index */
		ASSERT( nelem>=0 );

		if( !dense.lgElmtOn[nelem] )
		{
			fprintf( ioQQQ, "SAVE ELEMENT OPACITY: requested element is not active. Bailing out.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* generic driving of OpacityAdd1Element */
		iso_sp[ipH_LIKE][nelem].st[ipH1s].Pop() = 1.;
		iso_sp[ipH_LIKE][nelem].st[ipH1s].DepartCoef() = 0.;
		if( nelem > ipHYDROGEN )
		{
			iso_sp[ipHE_LIKE][nelem].st[ipH1s].Pop() = 1.;
			iso_sp[ipHE_LIKE][nelem].st[ipH1s].DepartCoef() = 0.;
		}

		for( ion=0; ion <= nelem; ion++ )
		{
			for( j=0; j < (nelem + 2); j++ )
			{
				dense.xIonDense[nelem][j] = 0.;
			}

			dense.xIonDense[nelem][ion] = 1.;

			OpacityZero();

			/* generate opacity with standard routine - this is the one
			 * called in OpacityAddTotal to make opacities in usual calculations */
			OpacityAdd1Element(nelem);

			/* start filename for output */
			strcpy( chFileName , elementnames.chElementNameShort[nelem] );

			/* continue filename with ionization stage */
			strcat( chFileName , chNumbers[ion+1] );

			/* end it wil .opc, name will be of form HYDR.opc */
			strcat( chFileName , ".opc" );

			/* now try to open it for writing */
			ioFILENAME = open_data( chFileName, "w" );

			ilow = Heavy.ipHeavy[nelem][ion];

			for( j=ilow-1; j < MIN2(rfield.nflux,continuum.KshellLimit); j++ )
			{
				/* photon energy in rydbergs */
				PrintE93( ioFILENAME , rfield.anu(j)*EVRYD );
				fprintf( ioFILENAME , "\t");

				/* cross section in megabarns */
				PrintE93( ioFILENAME, (opac.opacity_abs[j]+opac.OpacStatic[j])/1e-18 );
				fprintf( ioFILENAME , "\n");
			}
			/* close this ionization stage */
			fclose( ioFILENAME );
		}

		prtPunOpacSummary();
		cdEXIT(EXIT_SUCCESS);
	}

	return;
}

/* print final information about where opacity files are */
STATIC void prtPunOpacSummary(void)
{
	fprintf(ioQQQ,"\n\nThe opacity files have been successfully created.\n");
	fprintf(ioQQQ,"The files have names that start with the first 4 characters of the element name.\n");
	fprintf(ioQQQ,"There is one file per ion and the number after the element name indicates the ion.\n");
	fprintf(ioQQQ,"The energies are in eV and the cross sections in megabarns.\n");
	fprintf(ioQQQ,"All end in \".opc\"\n");
	fprintf(ioQQQ,"The data only extend to the highest energy in this continuum source.\n");
}
