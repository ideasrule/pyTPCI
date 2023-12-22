/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  #if SETUP == 1 
  /* Simple gyration */ 
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 1.e6;
  #endif

  #if SETUP == 2
  /* ExB Drift */
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
  
  #if SETUP == 3
  /* Gradient drift */
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 1.*(1.0 + x1/1.e-3);
  //v[BX3] = 1.*(1.0 + x1/1.e40);
  #endif

  #if SETUP == 4
  /* Curvature drift */
  double L = 1., B_0 = 1.;
  double Bz = 0.;
  
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
 
  v[BX1] = B_0*x2/L;
  v[BX2] = B_0*x1/L;
  v[BX3] = Bz;
  #endif
  
  #if SETUP == 5
  /* Resistive E drift */
  /* ResRMHD is mandatory */
  #if PHYSICS != ResRMHD
	print("\n! Error in Init.c: PHYSICS is not ResRMHD!\n");
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
	double buf[7];	//needed to string buffer quantities to raw binary file
	double gamma, analiticgamma, gammaerror;
	double mu, analiticmu, muerror;
	double E, B;
	double analiticx, xerror, analiticv, verror;
	Particle *p;
	particleNode *curr;
	FILE *fp;
	static double ***w;
	if (w == NULL){
    w = ARRAY_BOX (-1, 1, -1, 1, -1, 1, double);
  }
  
	/*Loop over particles*/
  PARTICLES_LOOP(curr, d->PHead){
    p = &(curr->p);

  /* -- Get particle coordinates (x and v) in Lab frame -- */
    
    x  = p->coord[IDIR];  y = p->coord[JDIR];  z = p->coord[KDIR];
    vx = p->speed[IDIR]; vy = p->speed[JDIR]; vz = p->speed[KDIR];
    t  = g_time;
	
	/*Compute quantities*/
	
	//-----------------E setup, checking for gamma and position error
	#if SETUP == 5 && PARTICLES_CR_GC == YES
	E = d->Vc[EX1][KBEG][JBEG][IBEG];
	gamma = sqrt(vx*vx + 1.);
	analiticgamma = sqrt(1. + t*t*E*E);
	analiticx = (analiticgamma - 1.)/E;
	analiticv = E*t;
	gammaerror = fabs(gamma - analiticgamma)/analiticgamma;
	xerror = fabs(x - analiticx)/analiticx;
	verror = fabs(analiticv - vx)/analiticv;
	#elif SETUP == 5 && PARTICLES_CR_GC == NO
	E = d->Vc[EX1][KBEG][JBEG][IBEG];
	gamma = sqrt(1. + vx*vx);//(sqrt(1. - vx)*sqrt(1. + vx));
	analiticv = E*t;
	analiticgamma = sqrt(1. + t*t*E*E);
	analiticx = (analiticgamma - 1.)/E;
	verror = fabs(analiticv - vx);
	xerror = fabs(analiticx - x)/analiticx;
	gammaerror = fabs(gamma - analiticgamma);
	//print("t\t%12.6e\n\n",E);
	#endif
	
	//-----------------B setup, checking for gamma error
	#if SETUP == 1 && PARTICLES_CR_GC == YES
	analiticgamma = 1.e6;
	B = d->Vc[BX3][KBEG][JBEG][IBEG];
	analiticmu = 0.5*(1. - 5.e-13)*(1. - 5.e-13)*1.e12/B;
	gamma = sqrt(2.*vz*B + 1.);
	mu = vz;
	gammaerror = fabs(gamma - analiticgamma)/analiticgamma;
	muerror = fabs(mu - analiticmu)/analiticmu;
	#elif SETUP ==  1 && PARTICLES_CR_GC == NO
	analiticgamma = 1.e6;
	B = d->Vc[BX3][KBEG][JBEG][IBEG];
	analiticmu = 0.5*(1. - 5.e-13)*(1. - 5.e-13)*1.e12/B;
	v_mod = sqrt(vx*vx + vy*vy + vz*vz);
	gamma = sqrt(1. + v_mod*v_mod);
	gammaerror = fabs(gamma - analiticgamma)/analiticgamma;
	mu = 0.5*(v_mod*v_mod)/B;
	muerror = fabs(mu - analiticmu)/analiticmu;
	#endif
	
	//-----------------ExB setup, checking for mu conservation
	#if SETUP == 2 && PARTICLES_CR_GC == YES
	double analiticgammaE = 100.;
	E = d->Vc[EX1][KBEG][JBEG][IBEG];
	B = d->Vc[BX3][KBEG][JBEG][IBEG];
	//boost gamma to get mustar
	analiticgamma = 1./sqrt(1. - E*E);
	analiticmu = 0.5*E*E/(B*analiticgamma);	// v/c = E/B -> v = E
	mu = vy;	//it is really mustar
	muerror = fabs(mu - analiticmu)/analiticmu;
	xerror = fabs(y + t*E)/t*E;	//y position error
	#elif SETUP == 2 && PARTICLES_CR_GC == NO
	double k = 100.;
	B = d->Vc[BX3][KBEG][JBEG][IBEG];
	E = d->Vc[EX1][KBEG][JBEG][IBEG];
	double Bprime = B/100.;
	v_mod = sqrt(vx*vx + vy*vy + vz*vz);
	gamma = sqrt(v_mod*v_mod + 1.);	//vmod is fourvelocity
	//!!!!!!!!v is now velocity
	vy/=gamma;
	vx/=gamma;
	double gammaprime = k*gamma*(1. - vy*E*B);	//p91 Ripperda phd (2018)
	//compute vx and vy in drift frame
	double vxprime = vx*(1./(1. - vy*E*B))/k;
	double vyprime = (vy/k - E*B + k*vy*E*B*E*B/(k + 1.))/(1. - E*B*vy);
	double R_c = gammaprime*sqrt(vxprime*vxprime + vyprime*vyprime)/Bprime;
	#endif
	
	//-----------------B gradient setup, comparing non-rel velocity to analytic
	double max_vel = 1.e-5;
	double n_particles = np_glob;	//20
	#if SETUP == 3 && PARTICLES_CR_GC == NO
		Particles_GetWeights(p, p->cell, w, grid); 	
		B = Particles_Interpolate(d->Vc[BX3], w, p->cell);
		analiticv = (n_particles - p->id)*max_vel/n_particles;	//make sure particle number and max v are correct, depends on initial setup
		analiticv *= sqrt(vy*vy + vx*vx)/(2*B*1.e-3);
		analiticx = analiticv*t;
		xerror = fabs(analiticx - y)/analiticx;
		analiticv = (n_particles - p->id)*max_vel/(n_particles*B);
	#elif SETUP == 3 && PARTICLES_CR_GC == YES
		Particles_GetWeights(p, p->cell, w, grid); 
		B = Particles_Interpolate(d->Vc[BX3], w, p->cell);
		analiticv = (n_particles - p->id)*max_vel/n_particles; //is Rc =v0
		analiticv *= analiticv/(2*B*1.e-3);
		analiticx = analiticv*t;	//x is actually y
		xerror = fabs(analiticx - y)/analiticx;
		analiticv = (n_particles - p->id)*max_vel/(n_particles*B);
	#endif
	buf[0] = t;
	buf[1] = x;
	buf[2] = y;
	buf[3] = z;
	buf[4] = vx;
	buf[5] = vy;
	buf[6] = vz;

	//-----------------B null setup, comparing difference in position
	#if SETUP == 4
	double theta = 2*CONST_PI*(np_glob - p->id)/np_glob;
	buf[4] = theta;
	#endif


	/*!-------------------CONVERGENCE RK4 TEST------------------*/
	#if CONVERGENCE_TEST == YES	
	int overwrite = 0;
		#if PARTICLES_CR_GC_RK4 == YES
	sprintf (file_name,"EconvergenceRK4.dat");
		#else
	sprintf (file_name,"Econvergence.dat");
		#endif
    if (overwrite) { 
      fp = fopen(file_name,"w");
      fprintf (fp,"#  %9s  %12s\n",
                   "t","gammaerror");
      fprintf (fp,"# ---------------------------------------");
      fprintf (fp,"-----------------------------------------");
      fprintf (fp,"-----------------------------------------\n");  
    }
	if(g_stepNumber > 1){	//prints only last result
		fp = fopen(file_name,"a");
		fprintf (fp,"%15.8e  %15.8e\n",g_dt, gammaerror);
	}
    fclose(fp);	
	
	/*!-------------------RAW BINARY OUTPUT------------------*/
	#elif  BINARY_OUT == YES
	sprintf (file_name,"bin.particle.%02d.dat", p->id);
	if(first_call) {
		fp = fopen(file_name,"w");
		fclose(fp);
	}
	else if (g_stepNumber > 1){
		FileWriteArray(buf, 0, 7, sizeof(double), file_name);
	}
	
	/*!-------TEST 4, BINARY PARTICLE DATA-------*/
	#elif SETUP == 4
		#if PARTICLES_CR_GC == NO
		sprintf (file_name, "bin.particle.test4.dat");
		#else 
		sprintf (file_name, "bin.particle.test4.GCA.dat");
		#endif
	if (first_call) {
		fp = fopen(file_name,"w");
		fclose (fp);
	} else if (g_stepNumber > 1){
		FileWriteArray (buf, 0, 7, sizeof(double), file_name);
	}
	fclose (fp);
	
	/*-------------------PRINT USEFUL RESULTS------------------*/
	#else
		#if PARTICLES_CR_GC == YES
	sprintf (file_name,"GCA.particle.%02d.dat", p->id);
		#else 
	sprintf (file_name,"particle.%02d.dat", p->id);
		#endif
	if (first_call) { 
      fp = fopen (file_name,"w");
      fprintf (fp,"#  %9s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
                   "t","x","y","z","theta","analiticx","xerror","gamma","analiticgamma","gammaerror","muerror");
      fprintf (fp,"# ---------------------------------------");
      fprintf (fp,"-----------------------------------------");
      fprintf (fp,"-----------------------------------------\n");  
    }else if (g_stepNumber > 1)    fp = fopen (file_name,"a");
 
    fprintf (fp,"%12.6e  %12.6e  %12.6e  %12.6e  ",t, x, y, z);
    fprintf (fp,"%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n"
	, vx, analiticx, xerror, gamma, analiticgamma, muerror, muerror);
    fclose (fp);
	#endif
	
	
	/*!----------------------TEST ON OLD AND NEW BORIS------------*/
	#if BORIS_TEST == YES
	sprintf (file_name,"boris.dat", p->id);
	if (first_call) { 
      fp = fopen(file_name,"w");
      fprintf (fp,"#  %9s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
                   "t","x","vx", "analiticv","deltav","gamma","analiticgamma","gammaerror","xerror");
      fprintf (fp,"# ---------------------------------------");
      fprintf (fp,"-----------------------------------------");
      fprintf (fp,"-----------------------------------------\n");  
    }else     fp = fopen(file_name,"a");
 
    fprintf (fp,"%12.6e  %12.6e  ",t, x);
    fprintf (fp,"%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n"
	, vx, analiticv, verror, gamma, analiticgamma, gammaerror, xerror);
    fclose(fp);
	#endif
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
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){}
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
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
