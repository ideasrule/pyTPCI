/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef LINES_H_
#define LINES_H_

#include <cstdio>
#include "module.h"
#include "transition.h"
#include "atmdat_adfa.h" // For NRECCOEFCNO

class cdstream;

class LineID
{
public:
	string chLabel;
	realnum wave;
	// the remaining parameters are optional
	int indLo;
	int indHi;
	realnum ELo;
	LineID() : wave(-1_r), indLo(-1), indHi(-1), ELo(-1_r) {}
	LineID(string lbl, realnum wv) : chLabel(lbl), wave(wv), indLo(-1), indHi(-1), ELo(-1_r) {}
	LineID(string lbl, realnum wv, realnum e) : chLabel(lbl), wave(wv), indLo(-1), indHi(-1), ELo(e) {}
	LineID(string lbl, realnum wv, int ilo, int ihi) : chLabel(lbl), wave(wv), indLo(ilo), indHi(ihi), ELo(-1_r) {}
};

/**lines main routine to put emission line intensities into line stack */
void lines(void);

/** general information at start of lines */
void lines_general(void);

/** the hydrogenic iso-sequence */
void lines_hydro(void);

/** create vectors to save line intensities */
void LineStackCreate(void);

/** information about grains */
void lines_grains(void);

/**lines_setup convert level 1 and level 2 line parameters and pointers 
 * into internal form used by code */
void lines_setup(void);

/** enter all continua */
void lines_continuum(void);

/** enter all molecules into emission line stack */
void lines_molecules(void);

/** enter all helium iso seq into emission line stack */
void lines_helium(void);

/**lines_lv1_li_ne place lines of elements lithium through neon into lines storage stack */
void lines_lv1_li_ne(void);

/**lines_lv1_k_zn place lines of elements potassium and heavier into lines storage stack */
void lines_lv1_k_zn(void);

/** routine to stuff comments into the stack of comments,
 * return is index to this comment */
long int StuffComment( const char * chComment );

/** lines_table invoked by table lines command, check if we can find all lines in a given list
 * returns 0 if ok, n is n lines not found */
int lines_table();

/** clear the name of the table read by the table lines command. */
void clear_lines_table();

class LinSv;

void cdEmis(
	const LinSv* line,
	/* the vol emissivity of this line in last computed zone */
	double *emiss ,
	// intrinsic or emergent
	bool lgEmergent );

static const int NHOLDCOMMENTS = 100;

extern const realnum Hbeta_WavLen;

/** this struc is different from following since they are only pointer here, will be allocated 
 * to form a large array after number of lines is counted, but this is the final form */
struct t_LineSave : public module {
	const char *chName() const
	{
		return "LineSave";
	}
	void zero();
	void comment(t_warnings&) {}
	/** number of emission lines in main stack */

	/** nsum is current index, to last line entered in the array */
	/** linadd increments nsum before doing anything else */
	long int nsum;

	/** index to number of comments printed within the block of lines */
	long int nComment;

	/** there are three types of calls to lines() 
	 * ipass = -1, first call, only count number off lines
	 * ipass =  0, second pass, save labels and wavelengths
	 * ipass =  1, integrate intensity*/
	long int ipass;

	/** holds comment strings associated with various blocks of output lines */
	string chHoldComments[NHOLDCOMMENTS];

	/** NormWL is array index for emission line on normalize command */
	long int ipNormWavL;

	/** WavLNorm is wavelength of emission line on normalize command */
	realnum WavLNorm;

	/** number of significant figures for lines
	 * this affects all aspects of reading and writing lines */
	long int sig_figs;
	static const long sig_figs_max = 6;

	/** length of string wl not including units
	 *  typically, sig_figs+2 (dot & unit) */
	int	wl_length;

	/** ScaleNormLine is the scale factor for its appearance */
	double ScaleNormLine;

	/** chNormLab is optional label */
	char chNormLab[NCHLAB];

	/** flag saying whether norm has been set */
	bool lgNormSet;

	/** save rec coefficient data for recombination lines of C, N, O */
	realnum RecCoefCNO[4][NRECCOEFCNO];

	long findline(const LineID& line);

	/** number of lines allocated in emission line stack
	 * must not change between iterations or grid points */
	vector<LinSv> lines;
	vector<realnum> m_wavelength;
	vector<size_t> SortWL;
	void clear()
	{
		lines.clear();
		m_wavelength.clear();
	}
	void resize(long nlines);
	
	void setSortWL();
	void init(long index, char chSumTyp, const char *chComment, const char *label,
				 bool lgAdd, realnum wavelength, const TransitionProxy& tr);
	realnum wavelength(long index)
	{
		return m_wavelength[index];
	}

	void resetWavelength( long index, realnum wl )
	{
		m_wavelength[index] = wl;
	}

	bool lgIsoContSubSignif;
};
extern t_LineSave LineSave;

/** this struc is different from above since only pointer here, will be allocated 
 * to form a large array after number of lines is counted. 
 * these are the main line save arrays */ 
class LinSv {
	long m_index;
	char m_chSumTyp;
	char m_chALab[NCHLAB];
	char m_chCLab[NCHLAB];
	double m_SumLine[4];	
	/** the emissivity, per unit vol, for current conditions, */
	double m_emslin[2];
	string m_chComment;
	vector<long> m_component;
	TransitionProxy m_tr;
	enum { DEFAULT, SEPARATOR, UNIT, UNITD, INWARD, INWARDTOTAL, INWARDCONTINUUM,
			 COLLISIONAL, PUMP, HEAT, CASEA, CASEB, NINU, NFNU, PHOPLUS, PCON, QH }
		m_type;
	// private accessors
	void chALabSet(const char *that);
public:

	/** Return the TransitionProxy so that transition data can be extracted (col_str, Aul, etc.) */
	TransitionProxy getTransition()
	{
		return m_tr;
	}

	/** one char saying whether heat 'h', cooling 'c', information, 'i' */
	char chSumTyp() const
	{
		return m_chSumTyp;
	}

	/** the four char string label for the line */
	const char *chALab() const
	{
		return m_chALab;
	}
	/** the four char string label for the line, all caps */
	const char *chCLab() const
	{
		return m_chCLab;
	}
	/** integrated intensity of the line, 
	 * [0] is intrinsic, 
	 * [1] emergent 
	 * [2] is intrinsic, 
	 * [3] emergent 
	 */
	bool isBlend() const
	{
		return !m_component.empty();
	}
	char LineType() const
	{
		return	m_chSumTyp;
	}
	void addComponent(const string& species,const double wavelength);
	void addComponent(const LineID& line);
	void addComponentID(long id);
	void makeBlend(const char* species, const double wavelength, 
				   const double width);
	void setBlendWavl();

	const TransitionProxy getComponent(long ind)
	{
		long id = m_component[ind];
		return LineSave.lines[id].getTransition();
	}
	double SumLine(int i) const
	{
		if (!isBlend())
		{
			return m_SumLine[i];
		}
		else
		{
			double sum = 0.;
			for (size_t j=0; j<m_component.size(); ++j)
			{
				long id = m_component[j];
				sum += LineSave.lines[id].SumLine(i);
			}
			return sum;
		}
	}
	void SumLineAdd(int i, double val)
	{
		if (!isBlend())
			m_SumLine[i] += val;
	}
	void SumLineAccum(double cumulative_factor)
	{
		if (!isBlend())
		{
			for( long nEmType=0; nEmType<2; ++nEmType )
			{
				m_SumLine[nEmType+2] += cumulative_factor*m_SumLine[nEmType];
			}
		}
	}
	void SumLineZero()
	{
		if (!isBlend())
			m_SumLine[0] = m_SumLine[1] = 0.0;
	}
	void SumLineZeroAccum()
	{
		if (!isBlend())
			m_SumLine[2] = m_SumLine[3] = 0.0;
	}
	void SumLineThin()
	{
		if (!isBlend() && m_chSumTyp == 't')
			m_SumLine[1] = m_SumLine[0];
	}

	/** the emissivity, per unit vol, for current conditions, */
	double emslin(int i) const
	{
		if (!isBlend())
		{
			return m_emslin[i];
		}
		else
		{
			double sum = 0.;
			for (size_t j=0; j<m_component.size(); ++j) 
			{
				long id = m_component[j];
				sum += LineSave.lines[id].emslin(i);
			}
			return sum;
		}
	}
	void emslinZero()
	{
		if (!isBlend())
		{
			m_emslin[0] = 0.0;
			m_emslin[1] = 0.0;
		}
	}
	void emslinSet(int i, double v)
	{
		if (!isBlend())
			m_emslin[i] = v;
	}
	void emslinThin()
	{
		if (!isBlend() && m_chSumTyp == 't')
			m_emslin[1] = m_emslin[0];
	}

	/** the wavelength of the line */
	realnum wavelength() const
	{
		return LineSave.wavelength(m_index);
	}

	/** comment describing the line */
	string chComment() const;
	void init(long index, char chSumTyp, const char *chComment, const char *label,
				 const TransitionProxy& tr);

	void prt(FILE *fp) const;
	string label() const;
	string biglabel() const;
	bool isCat(const char *s) const;
	bool isSeparator() const
	{
		return m_type == SEPARATOR;
	}
	bool isUnit() const
	{
		return m_type == UNIT;
	}
	bool isUnitD() const
	{
		return m_type == UNITD;
	}
	bool isInward() const
	{
		return m_type == INWARD;
	}
	bool isInwardTotal() const
	{
		return m_type == INWARDTOTAL;
	}
	bool isInwardContinuum() const
	{
		return m_type == INWARDCONTINUUM;
	}
	bool isCollisional() const
	{
		return m_type == COLLISIONAL;
	}
	bool isPump() const
	{
		return m_type == PUMP;
	}
	bool isHeat() const
	{
		return m_type == HEAT;
	}
	bool isCaseA() const
	{
		return m_type == CASEA;
	}
	bool isCaseB() const
	{
		return m_type == CASEB;
	}
	bool isNInu() const
	{
		return m_type == NINU;
	}
	bool isNFnu() const
	{
		return m_type == NFNU;
	}
	bool isPhoPlus() const
	{
		return m_type == PHOPLUS;
	}
	bool isPcon() const
	{
		return m_type == PCON;
	}
	bool isQH() const
	{
		return m_type == QH;
	}

#ifndef NDEBUG
	void checkEmergent( const long ipEmType ) const
	{
		/* Pick instantaneous or cumulative */
		int ipIntr, ipEmer;
		if( ipEmType == 1 )
		{
			ipIntr = 0;
			ipEmer = 1;
		}
		else if( ipEmType == 3 )
		{
			ipIntr = 2;
			ipEmer = 3;
		}
		else
			return;

		if( !isBlend() )
		{
			if( m_chSumTyp != 't' || LineSave.wavelength(m_index) <= realnum(0.) )
				ASSERT( m_SumLine[ipEmer] == 0. );
			else
				ASSERT( m_SumLine[ipIntr] > m_SumLine[ipEmer] ||
					fp_equal( m_SumLine[ipIntr], m_SumLine[ipEmer], 10 ) );
		}
	}
#else
	void checkEmergent( const long ) const
	{
		(void)0;
	}
#endif
};

inline void t_LineSave::init(long index, char chSumTyp, const char *chComment, const char *label,
									  bool lgAdd, realnum wavelength, const TransitionProxy& tr)
{
	if (!lgAdd)
	{
		// number of lines OK, set parameters for first pass
		// negative wavelength means it is just label, possibly not correct
		wavelength = fabs(wavelength);
	}
	
	m_wavelength[index] = wavelength;
	lines[index].init(index,chSumTyp,chComment,label,tr);
}
inline void t_LineSave::resize(long nlines)
{
	lines.resize(nlines);
	m_wavelength.resize(nlines);
}

#endif /* LINES_H_ */
