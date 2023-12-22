/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Circularly polarized Alfven waves for RMHD

  Setup the initial conditions for the large-amplitude circularly
  polarized Alfven wave test, as in section 4.1 of Del Zanna (2007).

  Note: in 2D the solution is rotated around the $z$-axis by
  assuming one  wavelength in both directions.
  This implies that \f$ k_x = 2\pi/L_x\f$ and \f$ k_y = 2\pi/L_y\f$.
  The transformation from the 1D (unrotated) system with primed
  coordinates to the actual computational system is given by
  \f[
      \vec{V}  = R\vec{V}' \qquad{\rm where}\quad
      R = \left(\begin{array}{ccc}
              \cos\alpha  & -\sin\alpha   &  0 \\ \noalign{\medskip}
              \sin\alpha  &  \cos\alpha   &  0 \\ \noalign{\medskip}
                     0  &           0   &  1 \\ \noalign{\medskip}
          \end{array}\right) 
  \f]
  while the invers transormation is
  \f[
      \vec{V}'  = R^{-1}\vec{V} \qquad{\rm where}\quad
      R^{-1} = \left(\begin{array}{ccc}
                \cos\alpha  &  \sin\alpha   &  0 \\ \noalign{\medskip}
               -\sin\alpha  &  \cos\alpha   &  0 \\ \noalign{\medskip}
                       0  &           0   &  1 \\ \noalign{\medskip}
          \end{array}\right) 
  \f]
  Note that the wave phase
  \f[
      \phi = \vec{k}\cdot\vec{x} = \vec{k}'\cdot\vec{x}'
  \f]
  is invariant under rotations.

  \author A. Mignone (mignone@to.infn.it)
  \date   Feb 08, 2021

  \b References: 
     - Del Zanna et al, A&A (2007) 473, 11-30
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x, double y, double z)
/*
 *
 *
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double Bx, By, Bz, B0;
  double vx, vy, vz, vA;
  double eta, w, scrh;
  double c,s;
  double Lx = g_domEnd[IDIR] - g_domBeg[IDIR];
  double Ly = g_domEnd[JDIR] - g_domBeg[JDIR];
  double kx = 2.0*CONST_PI/Lx;  /* Choose wavevector so that we have one  */
  double ky = 2.0*CONST_PI/Ly;  /* wavelenght in each direction.          */
  double t  = ky/kx;
  double phi, yp;

  g_gamma = 4.0/3.0;

  B0 = eta = 1.0;
  us[RHO] = 1.0;
  us[PRS] = 1.0;
  
  w  = us[RHO] + g_gamma*us[PRS]/(g_gamma - 1.0);
  vA  = B0*B0/(w + B0*B0*(1.0 + eta*eta));
  vA /= 0.5*(1.0 + sqrt(1.0 - 4.0*eta*eta*vA*vA));
  vA  = sqrt(vA);

  #if DIMENSIONS == 1
  c   = 1.0;
  s   = 0.0;
  phi = kx*x;
  #elif DIMENSIONS == 2
  c   = 1.0/sqrt(1.0 + t*t);
  s   = t*c;
  phi = kx*x + ky*y;
  #endif

/* -------------------------------------------------
   1. Define solution in the unrotated (1D) frame
   ------------------------------------------------- */

  Bx = B0;
  By = eta*B0*cos(phi);
  Bz = eta*B0*sin(phi);
  
  vx = 0.0;
  vy = -vA*By/B0;
  vz = -vA*Bz/B0;

  yp = -x*s + y*c;  /* y' coordinate */

/* -------------------------------------------------
   2. Rotate vector
   ------------------------------------------------- */

  us[VX1] = vx*c - vy*s;
  us[VX2] = vx*s + vy*c;
  us[VX3] = vz;
 
  us[BX1] = Bx*c - By*s;
  us[BX2] = Bx*s + By*c;
  us[BX3] = Bz;

  double kmod = sqrt(kx*kx + ky*ky);
  us[AX1] = 0.0;
  us[AX2] = 0.0;
  us[AX3] = B0*( yp - eta*sin(phi)/kmod);

  if (first_call == 1){
    print ("o vA    = %18.12e\n",vA);
    print ("o omega = %18.12e\n",kmod*vA);
    print ("o T     = %18.12e\n",c/vA);
    first_call = 0;
  }
/*
restart;
B[t] := eta*B[0]*cos(k*(c*x+s*y));
B[x] := B[0]*c - B[t]*s;
B[y] := B[0]*s + B[t]*c;

A[z]:= B[0]*((-x*s+y*c) - eta*sin(k*(c*x+s*y))/k);
simplify( diff(A[z],y) - B[x]);
simplify(-diff(A[z],x) - B[y]);
*/


}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ************************************************************** */
void BACKGROUND_FIELD (real x1, real x2, real x3, real *B0)
/* 
 *
 * PURPOSE
 *
 *   Define the component of a static, curl-free background 
 *   magnetic field.
 *
 *
 * ARGUMENTS
 *
 *   x1, x2, x3  (IN)    coordinates
 *
 *   B0         (OUT)    vector component of the background field.
 *
 *
 **************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in/out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies on which side boundary conditions need 
 *                    to be assigned. side can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{ }

