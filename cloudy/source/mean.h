/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MEAN_H_
#define MEAN_H_

#include "container_classes.h"
#include "module.h"

/** used to store information for mean ionization and electron temperature */

struct t_mean : public module
{
	/** xIonMean[dim][nelem][ion][n]
	 * dim = 0 is radius integration, 1=area integration, 2=vol integration
	 * n=0 for Sum(quant*norm) and n=1 for Sum(norm) */
	multi_arr<double,4> xIonMean;
	/** following includes electron density */
	multi_arr<double,4> xIonEdenMean;

	/** TempIonMean[dim][nelem][stage][n]
	 * dim = 0 is radius integration, 1=area integration, 2=vol integration
	 * n=0 for Sum(quant*norm) and n=1 for Sum(norm) */
	multi_arr<double,4> TempIonMean;
	/** following includes electron density */
	multi_arr<double,4> TempIonEdenMean;

	/** mean B weighted by n(H0) / T as in 21 cm measures */
	multi_arr<double,2> TempB_HarMean;
	/** this is used to get the harmonic mean temperature with respect to
	 * atomic hydrogen density, for comparison with 21cm data */
	multi_arr<double,2> TempHarMean;
	/** this is mean of computed 21 cm spin temperature */
	multi_arr<double,2> TempH_21cmSpinMean;
	/** used to get mean temp */
	multi_arr<double,2> TempMean;
	/** used to get mean temp weighted by eden */
	multi_arr<double,2> TempEdenMean;

	map<string, multi_arr<double,2> > molecules;

	t_mean();

	/**MeanZero zero mean of ionization fractions array */
	void zero();
	void comment(t_warnings&) {}
	const char *chName() const
	{
		return "mean";
	}

	/**setup_molecules - Allocate memory for molecular species */
	void setup_molecules();

	/**MeanInc update the running averages */
	void MeanInc();

	/**MeanIon do mean of ionization fractions of any element 
	   \param chType either 'i' or 't' for ionization or temperature 
	   \param nelem atomic number on physical, no c, scale
	   \param dim type of average: 0=radius, 1=area, 2=volume
	   \param *n this will say how many of arlog have non-zero values
	   \param arlog[] array of values, log both cases 
	   \param lgDensity true, include electron density, false do not
	*/
	void MeanIon( char chType, long nelem, long dim, long *n, realnum arlog[], bool lgDensity ) const;

	/**MeanMoleculeTemp - compute mean temperature of molecular species
	 * \param chSpecies [in]	species label
	 * \param dim [in]		type of average: 0=radius, 1=area, 2=volume
	 * \param TeMean [out]		mean temperature of molecule
	 *
	 * \return bool			true if success, false otherwise
	 */
	bool MeanMoleculeTemp( const string &chSpecies, long dim, double &TeMean );
};

extern t_mean mean;


#endif /* MEAN_H_ */
