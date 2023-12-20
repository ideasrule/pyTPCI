/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ATMDAT_ADFA_H_
#define ATMDAT_ADFA_H_

#include "iso.h"
#include "atmdat.h"

const int NRECCOEFCNO=471;
typedef enum { PHFIT_UNDEF, PHFIT95, PHFIT96 } phfit_version;

class t_ADfA : public Singleton<t_ADfA>
{
	friend class Singleton<t_ADfA>;
protected: 
	t_ADfA();
private:
	phfit_version version;
	/* phfit.dat */
	static const int NSHELLS = 7;
	long int L[NSHELLS];
	long int NINN[LIMELM];
	long int NTOT[LIMELM];
	static const int NFIT_PH1 = 6,
			NFIT_PH2 = 7;
	realnum PH1[NSHELLS][LIMELM][LIMELM][NFIT_PH1];
	realnum PH2[LIMELM][LIMELM][NFIT_PH2];

	/** set_vshell_index -- Set valence shell index for given ion
	 * The index is used to access the photoionization cross-section
	 * data of Verner et al. (1995) or (1996).
	 *
	 * \param [in] ne   number of electrons
	 * \param [in] Z    atomic number of element
	 *
	 * \returns long    index of valence shell
	 */
	inline long set_vshell_index( long ne, long Z ) const;

	/* hpfit.dat */
	realnum PHH[NHYDRO_MAX_LEVEL][5];
	/* rec_lines.dat */
	realnum P[8][110];
	realnum ST[9][405];
	/* rad_rec.dat */
	realnum rrec[LIMELM][LIMELM][2];
	realnum rnew[LIMELM][LIMELM][4];
	realnum fe[13][3];
	/* h_rad_rec */
	realnum HRF[NHYDRO_MAX_LEVEL][9];
	/* h_phot_cs.dat */
	/** array of cross sections for photoionization of hydrogen at threshold,
	 * 0 is 1s, 1 is 2s, 2 is 2p, up to 400 */
	realnum STH[NHYDRO_MAX_LEVEL];
	/* coll_ion.dat */
	double CF[LIMELM][LIMELM][5];
	/* h_coll_str.dat */
	/** array of EIE cross sections for hydrogen atom.  
	 * For all E1 transitions nl - n'l', with n' < n <= 5 */
	/* >>refer	H1	cs	Anderson, H., Ballance, C.P., Badnell, N.R., & Summers, H.P., 
	 * >>refercon	2000, J Phys B, 33, 1255; erratum, 2002 */
	map<QNPair, array<realnum,NHCSTE>> HCS;
public:
	/** set_version set version of phfit data to be used
	    \param val
	*/
	void set_version(phfit_version val) { version = val; }

	/** get_version which version of phfit data should be used? */
	phfit_version get_version() const { return version; }

	/** ph1 access elements of PH1 data block with parameters for photoionization cross section fits
	    \param i
	    \param j
	    \param k
	    \param l
	*/
	realnum ph1(int i, int j, int k, int l) const { return PH1[i][j][k][l]; }

	/** getEthresh -- get ionization threshold for shell
	 * For inner shells use data of Verner et al. (1995) or (1996);
	 * for valence shells use NIST ionization potentials.
	 *
	 * \param nshell [in]  shell index (1 - 1s, 2 - 2s, etc, <= 7)
	 * \param nel    [in]  number of electrons in ion
	 * \param Z      [in]  atomic number of element
	 * \return ionization threshold, in eV
	 */
	double getEthresh( long nshell, long nel, long Z ) const;

	/** sth array of cross sections for photoionization of hydrogen at threshold,
	    0 is 1s, 1 is 2s, 2 is 2p, up to 400
	    \param i
	*/
	realnum sth(int i) const { return STH[i]; }

	/** phfit this subroutine calculates partial photoionization cross sections
	    for all ionization stages of all atoms from H to Zn (Z=30)
	    \param nz
	    \param ne
	    \param is
	    \param e
	    \author Dima Verner
	*/
	double phfit(long int nz, long int ne, long int is, double e);

	/** hpfit state specific photoionization cross sections for model hydrogen atom
	    \param iz 
	    \param n 
	    \param e 
	    \author Dima Verner
	*/ 
	double hpfit(long int iz, long int n, double e);

	/** rec_lines effective recombination coefficients for lines of C, N, O, by D. Verner
	    \param  t 
	    \param  r
	    \author Dima Verner
	*/ 
	void rec_lines(double t, realnum r[][NRECCOEFCNO]);

	/** rad_rec calculates rates of radiative recombination for all ions
	    \param iz nuclear number on physics scale
	    \param in number of recombined electrons
	    \param t temperature K
	    \author Dima Verner
	*/ 
	double rad_rec(long int iz, long int in, double t);

	/** H_rad_rec calculates state-specific recombination rates for H and H-like ions
	    \param iz 
	    \param n 
	    \param t
	    \author Dima Verner
	*/ 
	double H_rad_rec(long int iz, long int n, double t);

	/** coll_ion D Verner's routine to compute collisional ionization rate coefficients 
	    \param iz 
	    \param in 
	    \param t
	    \author Dima Verner
	*/
	double coll_ion(long int iz, long int in, double t);

	double coll_ion_wrapper(long int z, long int n, double t);

	/* coll_ion_hybrid computes hybrid collisional ionization rates */
	double coll_ion_hybrid(	long int z, long int n,	double t);

	/** h_coll_str routine to grab H cross sections from Anderson et al. 2002.
	    \param ipLo 
	    \param ipHi 
	    \param ipTe
	*/
	const realnum* h_coll_str(long nHi, long lHi, long nLo, long lLo) const;
};

#endif
