/* /////////////////////////////////////////////////////////////////// */
/*! \file  
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* ************************************************************************** */
void Visc_nu(double *v, double x1, double x2, double x3,
                        double *nu1, double *nu2)
/*! 
 *
 *  \param [in]      v  pointer to data array containing cell-centered quantities
 *  \param [in]      x1 real, coordinate value 
 *  \param [in]      x2 real, coordinate value 
 *  \param [in]      x3 real, coordinate value 
 *  \param [in, out] nu1  pointer to first viscous coefficient
 *  \param [in, out] nu2  pointer to second viscous coefficient
 *
 *  \return This function has no return value.
 * ************************************************************************** */
{
  double Re = 1;
  double U = 1e6 / UNIT_VELOCITY; //typical velocity at low radii
  double Delta = 0.0002; //grid spacing at low radii
  //double unit_viscosity = UNIT_DENSITY * UNIT_VELOCITY * UNIT_LENGTH;
  *nu1 = v[RHO] * (U * Delta / Re); 
  *nu2 = 0.0;
}
