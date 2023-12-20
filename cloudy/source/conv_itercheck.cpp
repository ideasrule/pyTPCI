/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ConvIterCheck check whether model has converged or whether more iterations
 * are needed - implements the iter to converg comnd */
#include "cddefines.h"
#include "taulines.h"
#include "iso.h"
#include "phycon.h"
#include "cddrive.h"
#include "mole.h"
#include "elementnames.h"
#include "dynamics.h"
#include "stopcalc.h"
#include "dense.h"
#include "iterations.h"
#include "colden.h"
#include "save.h"
#include "rt.h"
#include "conv.h"
#include "lines.h"

/*ConvIterCheck check whether model has converged or whether more iterations
 * are needed - implements the iterate to convergence command */
void ConvIterCheck( void )
{
	long int nelem, 
		i,
		ipISO,
		ipHi, ipLo;

	DEBUG_ENTRY( "ConvIterCheck()" );

	/* =======================================================================*/
	/* this is an option to keep iterating until it converges
	 * iterate to convergence option
	 * autocv is percentage difference in optical depths allowed,
	 * =0.20 in block data
	 * checking on Ly and Balmer lines */
	/*>>chng 04 oct 19, promote loop to do all iso-electronic series */
	iterations.lgOpticalDepthonverged = true;
	strcpy( conv.chNotConverged, "Converged!" );

	// set up intensities used to converge outward intensity of Hb
	static double HbFracOutOld=-1. , HbFracOutNew=-1.;
	HbFracOutOld = HbFracOutNew;

	double a, total, BeamedIn;
	long int ipTotal = cdLine( "H  1" , Hbeta_WavLen, &a , &total );
	long int ipInwd  = cdLine( "Inwd" , Hbeta_WavLen, &a , &BeamedIn );

	/* 2014 aug 23, mchatzikos
	 * when cdLine returned -37 for log of zero intensity, the ratio below evaluated to 0;
	 * preserve this behaviour now that cdLine returns linear intensities */
	HbFracOutNew = 0.;
	if( total > 0. )
		HbFracOutNew = 1. - BeamedIn / total;

	ASSERT( iteration == 1 || (HbFracOutNew>=0 && HbFracOutNew<=1.) );
	// this disables the test on the outward Hb
	ipInwd = -1;

	if( save.lgPunConv )
	{
		fprintf( save.ipPunConv, " iteration %li of %li\n" ,
				iteration, iterations.itermx);
	}

	bool lgReasonGiven = false;
	if( iteration > 1 && conv.lgAutoIt )
	{
		if( nzone>3 && ipInwd>=0 && ipTotal>=0 )
		{
			// check whether outward intensity of Hb has converged
			if( fabs(HbFracOutNew-HbFracOutOld)/HbFracOutNew> conv.autocv )
			{
				iterations.lgOpticalDepthonverged = false;
				sprintf( conv.chNotConverged, "change in outward Hb");
				if( save.lgPunConv )
				{
					lgReasonGiven = true;
					fprintf( save.ipPunConv, " Change in outward Hb, "
						"old=%.3e new=%.3e \n",
						HbFracOutOld , HbFracOutNew);
				}
			}
		}
		for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem < LIMELM; nelem++ )
			{
				if( dense.lgElmtOn[nelem] )
				{
					/* now check if major subordinate line is converged - for H-like this will
					 * be Ha, and for He-like, the 23P - 23S transition - this will not work for
					 * NISO > 2 so must check against this */
					if(ipISO==ipH_LIKE )
					{
						ipHi = ipH3p;
						ipLo = ipH2s;
					}
					else if( ipISO==ipHE_LIKE )
					{
						ipHi = ipHe2p3P2;
						ipLo = ipHe2s3S;
					}
					else
						/* fails when NISO increased, add more sequences */
						TotalInsanity();

					/* check both H-alpha and Ly-alpha for all species - 
					 * only if Balmer lines thick 
					 * so check if Ha optical depth significant */
					if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauIn() > 0.5 )
					{
						/* test if Lya converged, nLyaLevel is upper level of Lya for iso seq */
						if( fabs(iso_sp[ipISO][nelem].trans(iso_ctrl.nLyaLevel[ipISO],0).Emis().TauTot() -
							 iso_sp[ipISO][nelem].trans(iso_ctrl.nLyaLevel[ipISO],0).Emis().TauIn()*rt.DoubleTau) > 
						    conv.autocv*fabs(iso_sp[ipISO][nelem].trans(iso_ctrl.nLyaLevel[ipISO],0).Emis().TauIn()*rt.DoubleTau) )
						{
							/* not converged to within AUTOCV, normally 15 percent */
							iterations.lgOpticalDepthonverged = false;

							/* for iterate to convergence, print reason why it was not converged 
							* on 3rd and higher iterations */
							sprintf( conv.chNotConverged, "%s-like Lya",elementnames.chElementSym[ipISO] );

							if( save.lgPunConv )
							{
								lgReasonGiven = true;
								fprintf( save.ipPunConv, " %s-like Lya thick, "
									"nelem= %s iteration %li old %.3e new %.3e \n",
									elementnames.chElementSym[ipISO] ,
									elementnames.chElementSym[nelem], 
									iteration,
									iso_sp[ipISO][nelem].trans(iso_ctrl.nLyaLevel[ipISO],0).Emis().TauTot() ,
									iso_sp[ipISO][nelem].trans(iso_ctrl.nLyaLevel[ipISO],0).Emis().TauIn());
							}
						}

						if( fabs(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauTot() -
							 iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauIn()*rt.DoubleTau) >
						    conv.autocv*fabs(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauIn()*rt.DoubleTau) )
						{
							/* not converged to within AUTOCV, normally 15 percent */
							iterations.lgOpticalDepthonverged = false;

							/* for iterate to convergence, print reason why it was not converged 
							* on 3rd and higher iterations */
							sprintf( conv.chNotConverged, "%s-like subord",elementnames.chElementSym[ipISO] );

							if( save.lgPunConv )
							{
								lgReasonGiven = true;
								fprintf( save.ipPunConv, " %s-like subord, nelem= %s iteration %li old %.3e new %.3e \n" ,
									elementnames.chElementSym[ipISO],
									elementnames.chElementSym[nelem], 
									iteration,
									iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauTot() ,
									iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().TauIn()	);
							}
						}
					}
				}
			}
		}

		if(0)
		{
			// all database lines
			for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
			{
				if( dBaseSpecies[ipSpecies].lgActive )
				{
					for (TransitionList::iterator tr=dBaseTrans[ipSpecies].begin();
						  tr != dBaseTrans[ipSpecies].end(); ++tr)
					{
						if( (*tr).Emis().TauIn() > 1. &&
							 fabs((*tr).Emis().TauTot() - (*tr).Emis().TauIn()*rt.DoubleTau) >
							 conv.autocv*fabs((*tr).Emis().TauIn()*rt.DoubleTau) )
						{
							/* not converged to within AUTOCV, normally 15 percent */
							iterations.lgOpticalDepthonverged = false;
							
							/* for iterate to convergence, print reason why it was not converged
							 * on 3rd and higher iterations */
							sprintf( conv.chNotConverged, "database %s %.1f",
										dBaseSpecies[ipSpecies].chLabel ,
										(*tr).WLAng() );
							
							if( save.lgPunConv )
							{
								lgReasonGiven = true;
								fprintf( save.ipPunConv, " database  %s %.1f iteration %li old %.3e new %.3e \n",
											dBaseSpecies[ipSpecies].chLabel ,
											(*tr).WLAng(),
											iteration,
											(*tr).Emis().TauTot() ,
											(*tr).Emis().TauIn());
							}
						}
					}
				}
			}
		}

		if( conv.lgAllTransitions )
		{
			realnum baddiff=0., badscale=0.;
			bool lgTauConv=true;
			TransitionList::iterator badtr;

			// all lines
			for( vector<TransitionList>::iterator it = AllTransitions.begin(); it != AllTransitions.end(); ++it )
			{
				for (TransitionList::iterator tr=it->begin();
					  tr != it->end(); ++tr)
				{
					if ( tr->ipCont() <= 0 )
						continue;
					realnum diff = fabs((*tr).Emis().TauTot() - (*tr).Emis().TauIn()*rt.DoubleTau);
					realnum scale = fabs((*tr).Emis().TauIn()*rt.DoubleTau);
					if( (*tr).Emis().TauIn() > 1. && diff > conv.autocv*scale )
					{
						/* not converged to within AUTOCV, normally 15 percent */
						iterations.lgOpticalDepthonverged = false;
						if ( lgTauConv || diff*badscale > scale*baddiff )
						{
							lgTauConv = false;
							badscale = scale;
							baddiff = diff;
							badtr = tr;
						}
						if( save.lgPunConv )
						{
							lgReasonGiven = true;
							fprintf( save.ipPunConv, " database %s line %s iteration %li old %.3e new %.3e \n",
										(*it).chLabel().c_str(),
										chLineLbl(*tr).c_str() ,
										iteration,
										(*tr).Emis().TauTot() ,
										(*tr).Emis().TauIn()*rt.DoubleTau);
						}
					}
				}
			}
			if (!lgTauConv)
			{
				/* for iterate to convergence, print reason why it was not converged
				 * on 3rd and higher iterations */
				sprintf( conv.chNotConverged, "%s line '%s' %.3e=>%.3e",
							(*badtr).system().chLabel.c_str(),
							chLineLbl(*badtr).c_str(), 
							(*badtr).Emis().TauTot(), (*badtr).Emis().TauIn()*rt.DoubleTau);
			}
		}

		/* >>chng 03 sep 07, add this test */
		/* check on changes in major column densities */
		for( i=0; i<NCOLD; ++i )
		{
			/* was the species column density significant relative to
			 * the total H column density, and was its abundance changing? */
			if( colden.colden[i]/colden.colden[ipCOL_HTOT] > 1e-5 &&
			    fabs(colden.colden_old[i]-colden.colden[i]) > conv.autocv*colden.colden[i] )
			{
				/* not converged to within conv.autocv, normally 20 percent */
				iterations.lgOpticalDepthonverged = false;

				/* for iterate to convergence, print reason why it was not converged 
				 * on 3rd and higher iterations */
				strcpy( conv.chNotConverged, "H mole col" );

				if( save.lgPunConv )
				{
					lgReasonGiven = true;
					fprintf( save.ipPunConv, " H mole col species %li iteration %li old %.2e new %.2e H col den %.2e\n",
						i,iteration,
						colden.colden_old[i],
						colden.colden[i],
						colden.colden[ipCOL_HTOT] );
				}
			}
		}

		double biggestDiffer = 0.;
		/* >>chng 03 sep 07, add this test */
		/* check on changes in major column densities */
		for( i=0; i<mole_global.num_calc; ++i )
		{
 			if(mole_global.list[i]->isMonatomic())
 				continue;

 			/* was the species abundance and changing? */
 			double differ = (double)fabs(mole.species[i].column_old-mole.species[i].column) ;
 			if( (mole.species[i].column/colden.colden[ipCOL_HTOT] > 1e-5) &&
				(differ > conv.autocv*mole.species[i].column) )
			{
				/* not converged to within conv.autocv, normally 20 percent */
				iterations.lgOpticalDepthonverged = false;

				/* for iterate to convergence, print reason why it was not converged 
				 * on 3rd and higher iterations */
				if( differ > biggestDiffer )
				{
					strcpy( conv.chNotConverged, mole_global.list[i]->label.c_str() );
					strcat( conv.chNotConverged, " column" );
					/*fprintf(ioQQQ,"debugggreset\t CO mole %li %li %.2e %.2e\n",
						i,iteration,mole.species[i].column_old,mole.species[i].column);*/
					biggestDiffer = differ;
				}

				if( save.lgPunConv )
				{
					lgReasonGiven = true;
					fprintf( save.ipPunConv, "%s, old:%.3e new:%.3e\n" ,
						mole_global.list[i]->label.c_str(),
						mole.species[i].column_old ,
						mole.species[i].column );
				}
			}
		}

		/* check on dynamical convergence in wind model with negative velocity */
		if( dynamics.lgAdvection )
		{
			/* >>chng 02 nov 29, as per Will Henney email */
			/* NB NB
			 *
			 * Use of +1 here makes sure that one iteration is done after the
			 * initial relaxation iterations, using dynamical sources in the
			 * solvers and without declaring convergence prematurely.
			 */
			if( iteration <= dynamics.n_initial_relax+1 ||
			    dynamics.convergence_error > conv.autocv*dynamics.error_scale2*dynamics.convergence_tolerance ||
			    dynamics.discretization_error > conv.autocv*dynamics.error_scale2 )
			{
				iterations.lgOpticalDepthonverged = false;
				/* for iterate to convergence, print reason why it was not converged 
				 * on 3rd and higher iterations */
				strcpy( conv.chNotConverged, "Dynamics  " );
				if( save.lgPunConv )
				{
					lgReasonGiven = true;
					fprintf( save.ipPunConv, " Dynamics\n" );
				}
			}
		}

		if( save.lgPunConv )
		{
			if( iterations.lgOpticalDepthonverged )
				fprintf( save.ipPunConv, " conv_itercheck exits converged\n" );
			else
				fprintf( save.ipPunConv, " conv_itercheck exits NOT converged\n" );
		}

		/* lower limit to number of iterations if converged */
		if( iterations.lgOpticalDepthonverged )
			iterations.itermx = MIN2(iterations.itermx,iteration);

		/* test for stopping on first zone due to too-low temperature */
		if( phycon.te < StopCalc.TempLoStopZone && nzone == 1 )
		{
			iterations.lgOpticalDepthonverged = true;
			strcpy( conv.chNotConverged, "          " );
			iterations.itermx = MIN2(iterations.itermx,iteration);
		}

		/* Fails if we have not fully implemented save convergence reason -
		 * should generate string stating reason, and set this flat,
		 * if lgOpticalDepthonverged is set true.
		 * These tests are only done when save output is requested
		 */
		if( save.lgPunConv )
			if( !iterations.lgOpticalDepthonverged && !lgReasonGiven )
				TotalInsanity();
	}
	return;
}
