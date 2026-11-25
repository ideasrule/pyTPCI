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
  /* THERMO: curr_heat.dat rows are
       r  H/rho  C/rho  rho_tab  theta_tab  bhT  bhr  bcT  bcr
     (H/rho, C/rho per gram in cgs; rho_tab, theta_tab = P/rho in code units; b = local exponents).
     The rate is rescaled to the cell's current state so the thermostat responds between Cloudy calls. */
  #define NCOL 9
  static int ntab;
  static double *tab[NCOL];
  static double heating_unit = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
  int n;

  if (tab[0] == NULL){
    FILE *fcool;
    printLog (" > Reading table from disk...\n");
    fcool = fopen("curr_heat.dat","r");
    if (fcool == NULL){
      printLog ("! Radiat: curr_heat.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    for (n = 0; n < NCOL; n++) tab[n] = ARRAY_1D(20000, double);
    ntab = 0;
    while (fscanf(fcool, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", tab[0] + ntab, tab[1] + ntab,
                  tab[2] + ntab, tab[3] + ntab, tab[4] + ntab, tab[5] + ntab, tab[6] + ntab,
                  tab[7] + ntab, tab[8] + ntab) == NCOL) {
      ntab++;
    }
    fclose(fcool);
    if (ntab < 2){
      printLog ("! Radiat: curr_heat.dat has %d rows of %d columns.\n", ntab, NCOL);
      QUIT_PLUTO(1);
    }
  }

  int    klo=0, khi=ntab-1, kmid;
  double r = gCooling_x1, w, q[NCOL];
  double *radii = tab[0];
  if (r > radii[khi] || r < radii[klo]){
    rhs[RHOE] = 0.0;
    return;
  }
  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    if (r <= radii[kmid]) khi = kmid; else klo = kmid;
  }
  w = (r - radii[klo])/(radii[khi] - radii[klo]);
  for (n = 1; n < NCOL; n++) q[n] = tab[n][klo]*(1.0 - w) + tab[n][khi]*w;

  /* current state relative to the table state; clamped to the range the trigger allows (~10%) */
  double theta = (g_gamma - 1.0)*v[RHOE]/v[RHO];
  double x = theta/q[4], y = v[RHO]/q[3];
  x = MAX(0.85, MIN(1.18, x));
  y = MAX(0.85, MIN(1.18, y));

  rhs[RHOE] = q[1]*pow(x, q[5])*pow(y, q[6]) - q[2]*pow(x, q[7])*pow(y, q[8]);
  rhs[RHOE] = rhs[RHOE] * (v[RHO] * UNIT_DENSITY) * heating_unit;
}
