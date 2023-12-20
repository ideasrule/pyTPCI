/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef GRAINS_H_
#define GRAINS_H_

/** GrainDrive main routine to converge grains thermal solution */
void GrainDrive();

/** GrainDrift computes grains drift velocity */
void GrainDrift();

/** this routine is called by IterStart() */
void GrainStartIter();

/** this routine is called by IterRestart() */
void GrainRestartIter();

/** this routine is called by ParseSet() */
void SetNChrgStates(long);

/** startup routine for grains, called before first calculations, but after parsecommands */
void GrainsInit();

/** main routine for generating the grain diffuse emission */
void GrainMakeDiffuse();

/** main routine for quantum heating */
void qheat(/*@out@*/vector<double>&,/*@out@*/vector<double>&,/*@out@*/long*,size_t);

/** initialize interpolation arrays for grain enthalpy */
void InitEnthalpy();

struct GrainPar;

/** mie_write_opc
 \param [in] *rfi_file 
 \param [in] *szd_file 
*/
void mie_write_opc(/*@in@*/const char*,/*@in@*/const char*,long int);
/**read in the *.opc file with opacities and other relevant information
\param *chFile
\param gp
*/
void mie_read_opc(/*@in@*/const char*,const GrainPar&);
/**set up Gaussian quadrature for arbitrary interval
\param nn
\param xbot
\param xtop
\param x[]
\param a[]
\param rr[]
\param ww[]
*/
void gauss_init(long int,double,double,const vector<double>&,const vector<double>&,vector<double>&,vector<double>&);
/** set up abscissas and weights for Gauss-Legendre intergration of arbitrary even order 
\param nn
\param x[]
\param a[]
*/
void gauss_legendre(long int,vector<double>&,vector<double>&);
/** find index ind such that min(xa[ind],xa[ind+1]) <= x <= max(xa[ind],xa[ind+1]).
 * xa is assumed to be strictly monotically increasing or decreasing.
 * if x is outside the range spanned by xa, lgOutOfBounds is raised and ind is set to -1
 * n is the number of elements in xa. 
 \param x
 \param xa[]
 \param n
 \param [out] *ind
 \param [out] *lgOutOfBounds
*/
void find_arr(double,const vector<double>&,long int,/*@out@*/long int*,/*@out@*/bool*);

#endif /* GRAINS_H_ */
