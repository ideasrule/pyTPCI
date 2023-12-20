/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */


#ifndef ABUND_H_
#define ABUND_H_

#include "module.h"

 /**
  AbundancesSet sets initial abundances after parameters are entered by reading input
 */ 
void AbundancesSet(void);

/**
  AbundancesPrt print all abundances, both gas phase and grains
 */ 
void AbundancesPrt( void );

/**
  generate abundance set from Fred Hamann's starburst evolution grid 
 \param chCard
*/
class Parser;
void abund_starburst(Parser &p);

 /**
  AbundancesTable interpolate on table of points to do 'element table' command, 
  \param  r0 
  \param  depth 
  \param  iel
 */ 
double AbundancesTable(double r0, 
  double depth, 
  long int iel);



typedef struct isotope
{
private:
	int	Niso;		/* Number of isotopes for element */
	vector<int>	Aiso;
	vector<realnum>	abnd, mass, spin, magmom;

	/* Functions to access abundances */
	int getIsoInd (int Anum)
	{
		int j = -1;
		for( int i = 0; i < Niso; i++ )
		{
			if( Aiso[i] == Anum )
			{
				j = i;
				break;
			}
		}
		return	j;
	}
public:
	void init ( )
	{
		Niso = 0;
		return;
	}
	int getNiso ( )
	{
		return	Niso;
	}
	void setData (int Anum, realnum nmass, realnum nspin, realnum nmagmom)
	{
		Aiso.push_back( Anum );
		abnd.push_back( realnum( 0. ) );
		mass.push_back( nmass );
		spin.push_back( nspin );
		magmom.push_back( nmagmom );
		Niso++;
	}
	int getAiso (int iso)
	{
		int Aval = -1;
		if( iso >= 0 && iso < Niso )
			Aval = Aiso[iso];
		return	Aval;
	}
	realnum getMass (int Anum)
	{
		realnum	nmass = -1.0;
		int j = getIsoInd( Anum );
		if( j != -1 )
			nmass = mass[j];
		return	nmass;
	}
	realnum getSpin (int Anum)
	{
		realnum	nspin = -1.0;
		if( spin.size() == 0 )
			return	nspin;
		int j = getIsoInd( Anum );
		if( j != -1 )
			nspin = spin[j];
		return	nspin;
	}
	realnum getMagMom (int Anum)
	{
		realnum	nmagmom = -1.0;
		if( magmom.size() == 0 )
			return	nmagmom;
		int j = getIsoInd( Anum );
		if( j != -1 )
			nmagmom = magmom[j];
		return	nmagmom;
	}
	int setAbn (int Anum, realnum abn)
	{
		int j = getIsoInd( Anum );
		if( j != -1 )
			abnd[j] = abn;
		return	j;
	}
	realnum getAbn (int Anum)
	{
		realnum	abn = -1.0;
		int j = getIsoInd( Anum );
		if( j != -1 )
			abn = abnd[j];
		return	abn;
	}
	/* Express all isotope fractions in terms of the most abundant isotope. */
	void normAbn ( )
	{
		realnum fracsum = 0.;
		for( int j = 0; j < Niso; j++ )
		{
			fracsum += abnd[j];
		}

		fracsum = 1.0 / fracsum;
		for( int j = 0; j < Niso; j++ )
		{
			abnd[j] *= fracsum;
		}
		return;
	}
	bool isAnyIllegal (void)
	{
		bool isZero = false;
		for( int j = 0; j < Niso; j++ )
		{
			if( abnd[j] < 0.0 )
			{
				isZero = true;
				break;
			}
		}
		return	isZero;
	}
	void prtIso (FILE *fp)
	{
		for( int j = 0; j < Niso; j++ )
		{
			fprintf(fp, "%2d", Aiso[j]);
			if( j != Niso-1 )
				fprintf(fp, ", ");
		}
		fprintf(fp, "\n");
		return;
	}
	void prtIsoPairs (FILE *fp)
	{
		for( int j = 0; j < Niso; j++ )
		{
			fprintf(fp, "  (%2d, %.4e)", Aiso[j], abnd[j]);
		}
		fprintf(fp, "\n");
		return;
	}
	void rm_nuc_data ( )
	{
		spin.resize( 0 );
		magmom.resize( 0 );
		return;
	}
} t_isotope;



/** abund.h */
struct t_abund : public module {
	const char *chName() const
	{
		return "abund";
	}

	void zero();
	void comment(t_warnings&) {}
	/** logical flag saying whether to include this element in save output for AGN tables */
	bool lgAGN[LIMELM];

	/** indicates whether abundances still at reference value */
	bool lgAbnReference;

	/** set true if abundances set with command, if false then
	 * must set default solar abundances
	 */
	bool lgAbundancesSet;

	/** initial elemental abundances for the current calculation */
	realnum ReferenceAbun[LIMELM];

	/** final gas phase abundances for the current calculation */
	realnum GasPhase[LIMELM];

	/** number of atoms depleted into grains */
	realnum AtomsDepleted;

	/* isotope fractions */
	t_isotope IsoAbn[LIMELM];

	/**lgAbunTabl says whether this element is to have its abundance
	 *determined from a table (true) or stored constant (false)
	 *set true with element table command */
	bool lgAbunTabl[LIMELM], 

	  /** lgAbTaDepth says whether depth or radius, true is depth */
	  lgAbTaDepth[LIMELM], 

	  /** general flag saying this option turned on */
	  lgAbTaON;

#	define	LIMTABD	500

	/**AbTabFac abundances for element table*/
	realnum AbTabFac[LIMTABD][LIMELM], 

	/**AbTabRad depth scale 
	 *parameters for dlaw table command*/
	  AbTabRad[LIMTABD][LIMELM];

	long int nAbunTabl;

	/** scale factors to alter abundances of elements, set with element scale */
	realnum ScaleElement[LIMELM];

	/** grain depletion scale factor is unity unless METALS DEPLETE is used */
	realnum DepletionScaleFactor[LIMELM];

	/** lgDepln is true if depln used */
	bool lgDepln;

	/** scale factor for metals, set with metals command	 */
	realnum ScaleMetals;

	/** sum up depletion across all elements */
	double SumDepletedAtoms();
	};
extern t_abund abund;



#endif /* ABUND_H_ */
