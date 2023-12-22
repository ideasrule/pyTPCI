/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Numerical benchmarks for Guiding Center method (fluid background).

  Set plasma (background) initial conditions for standard reference tests
  for the guiding center approximation (GCA).

  The macro \c SETUP (set in \c definitions.h) is used to select the
  configuration:

  - \c SETUP = 1: simple gyration setup
  - \c SETUP = 2: ExB drift
  - \c SETUP = 3: Gradient drift
  - \c SETUP = 4: Curvature drift
  - \c SETUP = 5: E // B test
 
  \author A. Mignone   (mignone@to.infn.it)
          H. Haudemand

  \date   March 4, 2021
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 *
 *********************************************************************** */
{
  #if SETUP == 1     /* Simple gyration */ 
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 1.e6;
  #endif

  #if SETUP == 2     /* ExB Drift */
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 1.0;
  
  v[EX1] = 1. - 5.e-5;
  v[EX2] = 0.0;
  v[EX3] = 0.0;
  #endif
  
  #if SETUP == 3     /* Gradient drift */
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 1.*(1.0 + x1/1.e-3);
  //v[BX3] = 1.*(1.0 + x1/1.e40);
  #endif
 
  #if SETUP == 4    /* Magnetic null */
  double L = 1., B_0 = 1.;
  double Bz = g_inputParam[BZ_GUIDE];
  
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
 
  v[BX1] = B_0*x2/L;
  v[BX2] = B_0*x1/L;
  v[BX3] = Bz;
  #endif
  
  #if SETUP == 5   /* E // B test (ResRMHD is mandatory ) */
  #if PHYSICS != ResRMHD
  print("! Init(): PHYSICS must be set to  ResRMHD!\n");
  QUIT_PLUTO(1);
  #endif
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
 
  v[BX1] = 1.;
  v[BX2] = 0.0;
  v[BX3] = 0.0;
  
  v[EX1] = 1.;
  v[EX2] = 0.;
  v[EX3] = 0.;
  #endif

  v[RHO] = 1.0;
  #if HAVE_ENERGY
  v[PRS] = 1.0;
  #endif
  v[TRC] = 0.0;
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
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
  int dir;
  static int first_call = 1, counter = 0;
  int np_glob = RuntimeGet()->Nparticles_glob;
  double t, x, y, z, vx, vy, vz, v_mod;
  double c2 = PARTICLES_CR_C*PARTICLES_CR_C;
  char file_name[32];
  double gamma, gamma_ex, err_gamma;
  double mu, mu_ex, err_mu;
  double E, B;
  double x_ex, err_x, v_ex, err_v;
  Particle *p;
  particleNode *curr;
  FILE *fp;
  static double ***w;

  if (w == NULL){
    w = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);
  }
  
/* ----------------------------------------------
   1. Loop over particles
   --------------------------------------------- */

  PARTICLES_LOOP(curr, d->PHead){
    p = &(curr->p);

  /* -- Get particle coordinates (x and v) in Lab frame -- */
    
    x  = p->coord[IDIR];  y = p->coord[JDIR];  z = p->coord[KDIR];
    vx = p->speed[IDIR]; vy = p->speed[JDIR]; vz = p->speed[KDIR];
    t  = g_time;

  /* -- Setup #1 (const B): check error on gamma -- */

    #if SETUP == 1 && PARTICLES_CR_GC == YES
    gamma_ex = 1.e6;
    B = d->Vc[BX3][KBEG][JBEG][IBEG];
    mu_ex     = 0.5*(1. - 5.e-13)*(1. - 5.e-13)*1.e12/B;
    gamma     = sqrt(2.*vz*B + 1.);
    mu        = vz;
    err_gamma = fabs(gamma - gamma_ex)/gamma_ex;
    err_mu    = fabs(mu - mu_ex)/mu_ex;
    #elif SETUP ==  1 && PARTICLES_CR_GC == NO
    gamma_ex = 1.e6;
    B = d->Vc[BX3][KBEG][JBEG][IBEG];
    mu_ex     = 0.5*(1. - 5.e-13)*(1. - 5.e-13)*1.e12/B;
    v_mod     = sqrt(vx*vx + vy*vy + vz*vz);
    gamma     = sqrt(1. + v_mod*v_mod);
    err_gamma = fabs(gamma - gamma_ex)/gamma_ex;
    mu        = 0.5*(v_mod*v_mod)/B;
    err_mu    = fabs(mu - mu_ex)/mu_ex;
    #endif

  /* -- Setup #2 (ExB setup): checking for mu conservation -- */

    #if SETUP == 2 && PARTICLES_CR_GC == YES
    E = d->Vc[EX1][KBEG][JBEG][IBEG];
    B = d->Vc[BX3][KBEG][JBEG][IBEG];

    gamma_ex    = 1./sqrt(1. - E*E/(B*B));  /* boost gamma to get mustar */    
    gamma       = vy;                       /* gamma  */
    err_gamma   = fabs(gamma - gamma_ex)/gamma_ex; 
    err_x       = fabs(fabs( y ) - fabs(t*E))/(t*E+1.e-12);  /* y position error */
    #elif SETUP == 2 && PARTICLES_CR_GC == NO
    err_mu = err_gamma = err_x = 0.0; /* No comparison needed with Boris */
    #endif
	
  /* -- Setup #3 (grad(B)): comparing non-rel velocity to analytic -- */

    #if SETUP == 3
    err_mu = p->color;  /* We store color into err_mu for convenience */
    #endif

    #if SETUP == 4
    continue;
    #endif

  /* -- Setup #5 (E // B): check error on gamma and position -- */

    #if SETUP == 5 && PARTICLES_CR_GC == YES
    E = d->Vc[EX1][KBEG][JBEG][IBEG];
    gamma     = sqrt(vx*vx + 1.);
    gamma_ex  = sqrt(1. + t*t*E*E);
    x_ex      = (gamma_ex - 1.)/E;
    v_ex      = E*t;
    err_gamma = fabs(gamma - gamma_ex)/gamma_ex;
    err_x     = fabs(x - x_ex)/x_ex;
    err_v     = fabs(v_ex - vx)/v_ex;
    #elif SETUP == 5 && PARTICLES_CR_GC == NO
    E = d->Vc[EX1][KBEG][JBEG][IBEG];
    gamma = sqrt(1. + vx*vx);//(sqrt(1. - vx)*sqrt(1. + vx));
    v_ex   = E*t;
    gamma_ex = sqrt(1. + t*t*E*E);
    x_ex = (gamma_ex - 1.)/E;
    err_v = fabs(v_ex - vx);
    err_x = fabs(x_ex - x)/x_ex;
    err_gamma = fabs(gamma - gamma_ex);
    #endif
    

    #if PARTICLES_CR_GC == YES
      #if PARTICLES_CR_GC_INTEGRATOR == RK4
      sprintf (file_name,"setup%02d_id%02d_gcRK4.dat", SETUP, p->id);
      #else
      sprintf (file_name,"setup%02d_id%02d_gcRK2.dat", SETUP, p->id);
      #endif
    #else
      sprintf (file_name,"setup%02d_id%02d_boris.dat", SETUP, p->id);
    #endif
    if (first_call) { 
      fp = fopen (file_name,"w");
      fprintf (fp,"#  %9s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
                   "t","x","y","z","vx","vy","vz","err_mu","err_gamma", "err_x");
      fprintf (fp,"# ---------------------------------------");
      fprintf (fp,"-----------------------------------------");
      fprintf (fp,"-----------------------------------------\n");  
    } else fp = fopen (file_name,"a");

    fprintf (fp,"%12.6e  %12.6e  %12.6e  %12.6e  ",t, x, y, z);
    fprintf (fp,"%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n",  
              vx, vy, vz, err_mu, err_gamma, err_x);
    fclose (fp);
  }
  
  first_call = 0;
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
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
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
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
