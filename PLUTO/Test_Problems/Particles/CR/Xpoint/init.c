/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief X-point test particle acceleration.

  Test-particle acceleration near an X-type magnetic
  reconnection region (see sect. 4.6 of Mignone et al. 2018).

  Configurations #1, #2 correspond, respectively, to the cases
  without and with guide field.
  They follow the normalization indicated in the
  paper (vA = 1 = C/100, t = 1/Omega_L).
 
  Configurations #3, #4 replicate #1 and #2 with the resistive
  RMHD module and are re-normalized so that C = 1, t = 1/Omega_L.

  \author A. Mignone (mignone@to.infn.it),
          G. Mattia

  \b References
    - "MAGNETOHYDRODYNAMIC-PARTICLE-IN-CELL METHOD FOR COUPLING COSMIC RAYS 
       WITH A THERMAL PLASMA: APPLICATION TO NON-RELATIVISTIC SHOCKS"\n
       Bai et al., ApJ (2015) 809, 55

  \date   June 13, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x, double y, double z)
/*! 
 *
 *
 *********************************************************************** */
{
  double alpha = 4.0/(g_domEnd[JDIR] - g_domBeg[JDIR]);
  double beta  = 4.0/(g_domEnd[IDIR] - g_domBeg[IDIR]);
  double Bz = g_inputParam[BMAG_Z];
  double Ez = g_inputParam[EMAG_Z];
  double B0 = PARTICLES_CR_C/100.0;

  v[RHO] = 1.0;
  v[PRS] = 1.0;

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = alpha*y*y/2.0 - beta*x*x/2.0;

/* -- Mori et al 1998 -- */

  v[BX1] = B0*alpha*y;
  v[BX2] = B0*beta*x;
  v[BX3] = B0*Bz;
  
  v[EX1] = 0.0;   
  v[EX2] = 0.0;
  v[EX3] = B0*B0*Ez; /* If c != 1, the electric field is actually cE */

/* -- Zharkova 2011 (page 386-387, exepct for a minus sign) -- */  
/*
  v[BX1] = tanh(-alpha*y);
  v[BX2] = beta*x;
  v[BX3] = Bz;

  v[EX1] = 0.0;
  v[EX2] = 0.0;
  v[EX3] = Ez;
*/

  g_smallPressure = 1.e-5;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 *
 *
 *********************************************************************** */
{
  
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *
 *********************************************************************** */
{
  
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif

