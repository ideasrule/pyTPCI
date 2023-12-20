/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef YIELD_H_
#define YIELD_H_

/* yield.h */

/** there are fewer lines than this in the mewe fluores file */
#define	MEWE_FLUOR	12000

class t_yield : public Singleton<t_yield>
{
	friend class Singleton<t_yield>;
protected:
	t_yield();
private:
	/** block data of yields from Mewe paper
	 * frac_elec_eject( nelem, nstage, nshell, nelec )
	 * nelem is element number, 6 for carbon
	 * nstage is stage of ionization, 1 for atom
	 * nshell is shell number in dima notation
	 * nelec is fraction that yield that number of electrons
	 * n_elec_eject is number of electrons at most */
	realnum frac_elec_eject[LIMELM][LIMELM][7][10];
	long int n_elec_eject[LIMELM][LIMELM][7];

	/** there are fewer than MEWE_FLUOR different lines in the fluores file,
	 * remember the atomic number, ionization stage */
	int nfl_nelem[MEWE_FLUOR];
	int nfl_ion[MEWE_FLUOR];
	int nfl_nshell[MEWE_FLUOR];
	int nfl_ion_emit[MEWE_FLUOR];
	int nfl_nLine[MEWE_FLUOR];
	realnum fl_energy[MEWE_FLUOR];
	/** fluorescense yield */
	realnum fl_yield[MEWE_FLUOR];
	long int nfl_ipoint[MEWE_FLUOR];

	/** this is the total number of fluorescent lines */
	long int nfl_lines;

	/** this is set true with the "no auger" command, normally false */
	bool lgKillAuger;
public:
	void init_yield();

	realnum elec_eject_frac( long n, long i, long ns, long ne ) const
	{
		if( lgKillAuger )
			return ( ne == 0 ) ? 1.f : 0.f;
		else
			return frac_elec_eject[n][i][ns][ne];
	}
	long nelec_eject( long n, long i, long ns ) const
	{
		return lgKillAuger ? 1 : n_elec_eject[n][i][ns];
	}
	int nelem( long n ) const { return nfl_nelem[n]; }
	int ion( long n ) const { return nfl_ion[n]; }
	int nshell( long n ) const { return nfl_nshell[n]; }
	int ion_emit( long n ) const { return nfl_ion_emit[n]; }
	realnum energy( long n ) const { return fl_energy[n]; }

	realnum yield( long n ) const
	{
		if( lgKillAuger )
			return 0.;
		else
			return fl_yield[n];
	}

	void set_ipoint( long n, long val ) { nfl_ipoint[n] = val; }
	int ipoint( long n ) const { return nfl_ipoint[n]; }

	int nlines() const { return nfl_lines; }

	void kill_yield() { lgKillAuger = true; }
	void reset_yield() { lgKillAuger = false; }
};

#endif /* YIELD_H_ */
