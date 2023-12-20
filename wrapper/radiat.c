/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute right hand side for Tabulated cooling

  \authors A. Mignone (mignone@ph.unito.it)\n
           M. Sormani\n

 \b References

  \date   Apr 09, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

extern double gCooling_x1, gCooling_x2, gCooling_x3;

/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */
{
  static int ntab;
  static double *radii, *heating_rates;
  static double heating_unit = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
  
/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */
  //printf("r=%lf\n", gCooling_x1);
  if (radii == NULL){
    FILE *fcool;
    printLog (" > Reading table from disk...\n");
    fcool = fopen("curr_heat.dat","r");
    if (fcool == NULL){
      printLog ("! Radiat: cooltable.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    radii = ARRAY_1D(20000, double);
    heating_rates = ARRAY_1D(20000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", radii + ntab, 
                                       heating_rates + ntab)!=EOF) {
      ntab++;
    }
    fclose(fcool);
  }

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

  int    klo=0, khi=ntab-1, kmid;
  double r = gCooling_x1, r_mid, dr;
  if (r > radii[khi] || r < radii[klo]){
    rhs[RHOE] = 0.0;
    //printf("WARNING: out of bounds %lf\n", r);
    return;
    QUIT_PLUTO(1);
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    r_mid = radii[kmid];
    if (r <= r_mid){
      khi = kmid;
    }else if (r > r_mid){
      klo = kmid;
    }
  }

/* -----------------------------------------------
    Compute r.h.s
   ----------------------------------------------- */

  dr        = radii[khi] - radii[klo];
  rhs[RHOE] = heating_rates[klo] * (radii[khi] - r)/dr + heating_rates[khi]*(r - radii[klo])/dr;
  rhs[RHOE] = rhs[RHOE] * (v[RHO] * UNIT_DENSITY) * heating_unit;
  //printf("Setting heating: r=%lf, rate=%e, unit=%e\n", gCooling_x1, rhs[RHOE], heating_unit);
}

