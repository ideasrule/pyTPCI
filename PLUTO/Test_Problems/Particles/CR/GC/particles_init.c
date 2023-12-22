/* ///////////////////////////////////////////////////////////////////// */
/*!
 \file
 \brief Numerical benchmarks for Guiding Center method (particles).
 
  Set particle initial conditions for standard reference tests
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
void Particles_Init(Data *d, Grid *grid)
/*!
 *  Sets initial conditions on particles.
 *
 *  \param [in]    d       Pointer to the PLUTO data structure.
 *  \param [in]    grid    Pointer to the PLUTO grid structure.
 *
 *********************************************************************** */
{
  int i,j,k, np, dir, nc;
  int np_glob = RuntimeGet()->Nparticles_glob;
  int np_cell = RuntimeGet()->Nparticles_cell;
  static int first_call = 1;
  double xbeg[3], xend[3];
  Particle p;
  double theta, gamma;

  if (first_call){
    RandomSeed(time(NULL),0);
    first_call = 0;
  }

/* --------------------------------------------------------------
   1. Global initialization
   -------------------------------------------------------------- */

  if (np_glob > 0){

    for (dir = 0; dir < 3; dir++){
      xbeg[dir] = grid->xbeg_glob[dir];
      xend[dir] = grid->xend_glob[dir];
    }

    for (np = 0; np < np_glob; np++){
      Particles_LoadUniform(np, np_glob, xbeg, xend, p.coord);

      #if SETUP == 1 && PARTICLES_CR_GC_4VEL == NO  /* Simple gyration */ 
      double v_perp = 1. - 5.e-13;
      p.coord[IDIR] = 1.0;
      p.coord[JDIR] = 0.0;
      p.coord[KDIR] = 0.0;
      p.speed[IDIR] = 0.0;
      p.speed[JDIR] = -v_perp;
      p.speed[KDIR] = 0.0;
      #elif SETUP == 1 && PARTICLES_CR_GC_4VEL == YES
      double v_perp = 1.e6*(1. - 5.e-13);
      p.coord[IDIR] = 1.0;
      p.coord[JDIR] = 0.0;
      p.coord[KDIR] = 0.0;
      p.speed[IDIR] = 0.0;
      p.speed[JDIR] = -v_perp;
      p.speed[KDIR] = 0.0;
      #endif
          
      #if SETUP == 2        /* ExB Drift */
      p.coord[IDIR] = 0.0;
      p.coord[JDIR] = 0.0;
      p.coord[KDIR] = 0.0;
      p.speed[IDIR] = 0.0;
      p.speed[JDIR] = 0.0;
      p.speed[KDIR] = 0.0;
      #endif
          
      #if SETUP == 3        /* Gradient drift */
      double max_vel = 1.e-4;       //v is not 4vel, but error is 10^-9
      double v_0 = (max_vel*(np + 1.))/(np_glob);
      p.coord[IDIR] = 0.0;
      p.coord[JDIR] = 0.0;
      p.coord[KDIR] = 0.0;
      p.speed[IDIR] = -v_0;
      p.speed[JDIR] = 0.0;
      p.speed[KDIR] = 0.0;
      p.color       = v_0;
      #endif

      #if SETUP == 4      /* Curvature drift */
      double theta = 2*CONST_PI*np/np_glob;
      double v_0 = 0.1;
      v_0 /= sqrt(1. - v_0*v_0/(PARTICLES_CR_C*PARTICLES_CR_C)); //convert to 4vel
      p.coord[IDIR] = cos(theta);
      p.coord[JDIR] = sin(theta);
      p.coord[KDIR] = 0.0;
      p.speed[IDIR] = -v_0*cos(theta);
      p.speed[JDIR] = -v_0*sin(theta);
      p.speed[KDIR] = 0.0;
      p.color       = theta; /* Use color to store initial angle */ 
      #endif
          
      #if SETUP == 5  /* E // B test  */
      p.coord[IDIR] = 0.0;
    //p.coord[IDIR] = ( sqrt(1. + 5.e2*5.e2) - 1.);         //Modified initial position to test Boris scheme
      p.coord[JDIR] = 0.0;
      p.coord[KDIR] = 0.0;
      p.speed[IDIR] = 0.0;
      p.speed[JDIR] = 0.0;
      p.speed[KDIR] = 0.0;
      #endif

      Particles_Insert (&p, d, PARTICLES_CREATE, grid);
    }
  }

/* ------------------------------------------------------------------
   2. Cell by cell initialization.
      Note: You may use Particles_LoadRandom() to initialize
            velocity components but not spatial coordinates.
   ------------------------------------------------------------------ */

  if (np_cell > 0){

    DOM_LOOP(k,j,i){
      xbeg[IDIR] = grid->xl[IDIR][i]; xend[IDIR] = grid->xr[IDIR][i];
      xbeg[JDIR] = grid->xl[JDIR][j]; xend[JDIR] = grid->xr[JDIR][j];
      xbeg[KDIR] = grid->xl[KDIR][k]; xend[KDIR] = grid->xr[KDIR][k];


      for (np = 0; np < np_cell; np++){

      /* -- Spatial distribution -- */

        Particles_LoadUniform(np, np_cell, xbeg, xend, p.coord);

     /* -- Velocity distribution -- */

        #if (PARTICLES == PARTICLES_CR)
        p.speed[IDIR] = 0.1;
        p.speed[JDIR] = 0.0001; 
        p.speed[KDIR] = 0.0;
        #endif
        Particles_Insert (&p, d, PARTICLES_CREATE, grid);
      }
    }
  }  
  Particles_SetID(d->PHead);
}

/* ********************************************************************* */
void Particles_Inject(Data *data, Grid *grid)
/*!
 *  Inject particles as you wish.
 *
 *  \param [in]  data    Pointer to the PLUTO data structure.
 *  \param [in]  grid    Pointer to the PLUTO grid structure.
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Particles_UserDefBoundary(Data *d, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
  int    dir;
  double xbeg[3], xend[3];
  particleNode *curr = d->PHead, *next;
  Particle *p;
  
  for (dir = 0; dir < 3; dir++) {
    xbeg[dir] = grid->xbeg_glob[dir];
    xend[dir] = grid->xend_glob[dir];
  }

  if (side == X1_BEG){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X1_END){
    dir = IDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X2_BEG){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X2_END){
    dir = JDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X3_BEG){
    dir = KDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] < xbeg[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }

  if (side == X3_END){
    dir = KDIR;
    while (curr != NULL){
      p    = &(curr->p);  
      next = curr->next; 
      if (p->coord[dir] > xend[dir]) Particles_Destroy (curr,d);
      curr = next;
    }
  }
}
