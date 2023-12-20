/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef TRANSITION_H_
#define TRANSITION_H_

#include "emission.h"
#include "collision.h"
#include "quantumstate.h"

/** \todo 2 bring these two together. */
/**AddLine2Stack add generic emission line to GenericLines and return pointer to that state. */

/* create a dummy emission structure.  Non-radiative transitions will point to this */
extern EmissionProxy DummyEmis;

// Proxy class provides 'object' access to struct-of-vectors
class TransitionList;
class TransitionListImpl;
class TransitionConstProxy;
class ExtraInten;

class TransitionProxy
{
public:
	typedef TransitionListImpl list_type;
	typedef ProxyIterator<TransitionProxy,TransitionConstProxy> iterator;
	typedef ProxyIterator<TransitionConstProxy,TransitionConstProxy> const_iterator;  
private:
	friend class  ProxyIterator<TransitionProxy,TransitionConstProxy>;
	list_type *m_list;
	int m_index;
public:
	explicit TransitionProxy(): m_list(0), m_index(-1) {}
	explicit TransitionProxy(TransitionListImpl *list, int index) : 
	m_list(list), m_index(index) {}
	void copy(const TransitionProxy& other) const;
	qList::iterator Lo() const;
	qList::iterator Hi() const;
	void setLo(int ipLo) const;
	void setHi(int ipHi) const;
	void AddLine2Stack() const;
	EmissionList::reference Emis() const;	
	int &ipEmis() const;
	string chLabel() const;
	void setComment( const string &comment ) const;
	string &getComment() const;
	bool associated() const
	{
		return m_list != NULL && m_index >= 0;
	}
	bool isSameAs(const TransitionProxy& other) const
	{
		return m_list == other.m_list && m_index == other.m_index;
	}
	bool hasEmis() const
	{
		return ipEmis() != -1;
	}
	void resetEmis() const
	{
		Emis().ipTran() = m_index;
	}
	bool systemIs(const TransitionList *query) const;
	const TransitionListImpl& system() const
	{
		return *m_list;
	}
	void check() const
	{
		ASSERT(!hasEmis() || Emis().ipTran() == m_index);
	}
	CollisionProxy Coll() const;

	/** wavelength, usually in Angstroms, used for printout, can be any units */
	realnum &WLAng() const;

	/** transition energy in degrees kelvin*/
	realnum EnergyK() const
	{
		return (realnum)T1CM*EnergyWN();
	}
	/** transition energy in ergs */
	realnum EnergyErg() const
	{
		return (realnum)ERG1CM*EnergyWN();
	}
	/** transition energy in ergs */
	double EnergyRyd() const
	{
		return WAVNRYD*EnergyWN();
	}
	/** vacuum wavelength in Angstroms */
	realnum EnergyAng() const
	{
		return 1e8f/EnergyWN();
	}


	/** transition energy in wavenumbers */
	realnum &EnergyWN() const;

	/** index for line within continuum array,
	 * this is on the f, not c, scale,
	 * negative ipCont means this is not a radiative transition, 
	 * and is used as a sentnecl */
	long &ipCont() const;

   
	/**set all elements of transition struc to dangerous values  
		\param *t
	*/
	void Junk() const;

	/**TransitionZero set all elements of transition struc to zero 
		\param *t
	*/
	void Zero() const;

	/**outline - adds line photons to reflin and outlin */
	void outline( double nonScatteredFraction, bool lgDoChecks ) const ;

	   
	/** outline_resonance - adds line photons to reflin and outlin,
	  * setting nonScatteredFraction as default for resonance lines */
	void outline_resonance(  ) const ;
	int &ipLo() const;
	int &ipHi() const;
	/**AddState2Stack add generic quantum state to GenericStates and return pointer to that state. */
	void AddHiState() const;
	void AddLoState() const;
	realnum width() const;
	list_type *list() const
	{
		return m_list;
	}
};

class TransitionConstProxy
{
public:
	typedef const TransitionListImpl list_type;
	typedef ProxyIterator<TransitionConstProxy,TransitionConstProxy> iterator;
	typedef ProxyIterator<TransitionConstProxy,TransitionConstProxy> const_iterator;  
private:
	friend class  ProxyIterator<TransitionConstProxy,TransitionConstProxy>;
	const list_type *m_list;
	int m_index;
public:
	explicit TransitionConstProxy(): m_list(0), m_index(-1) {}
	explicit TransitionConstProxy(const TransitionListImpl *list, int index) : 
	m_list(list), m_index(index) {}
	void copy(const TransitionConstProxy& other) const;
	qList::iterator Lo() const;
	qList::iterator Hi() const;
	void AddLine2Stack() const;
	EmissionList::reference Emis() const;	
	int ipEmis() const;
	string getComment() const;
	bool associated() const
	{
		return m_list != NULL && m_index >= 0;
	}
	bool hasEmis() const
	{
		return ipEmis() != -1;
	}
	void check() const
	{
		ASSERT(!hasEmis() || Emis().ipTran() == m_index);
	}
	CollisionProxy Coll() const;

	/** wavelength, usually in Angstroms, used for printout, can be any units */
	realnum WLAng() const;

	/** transition energy in degrees kelvin*/
	realnum EnergyK() const
	{
		return (realnum)T1CM*EnergyWN();
	}
	/** transition energy in ergs */
	realnum EnergyErg() const
	{
		return (realnum)ERG1CM*EnergyWN();
	}
	/** transition energy in Ryd */
	double EnergyRyd() const
	{
		return WAVNRYD*EnergyWN();
	}
	/** vacuum wavelength in Angstroms */
	realnum EnergyAng() const
	{
		return 1e8f/EnergyWN();
	}


	/** transition energy in wavenumbers */
	realnum EnergyWN() const;

	/** index for line within continuum array,
	 * this is on the f, not c, scale,
	 * negative ipCont means this is not a radiative transition, 
	 * and is used as a sentnecl */
	long ipCont() const;

	/**outline - adds line photons to reflin and outlin */
	void outline( double nonScatteredFraction, bool lgDoChecks ) const ;

	/**outline_resonance - adds line photons to reflin and outlin,
	   setting nonScatteredFraction as default for resonance lines */
	void outline_resonance(  ) const ;
	int ipLo() const;
	int ipHi() const;
};


// Structure-of-vectors for transition data
class TransitionListImpl
{
	vector<int> ipHi, ipLo;
	vector<long> ipCont;
	CollisionList Coll;
	vector<realnum> EnergyWN, WLAng;
	vector<string> chComment;
	// DO NOT IMPLEMENT
	TransitionListImpl(const TransitionListImpl&);
	TransitionListImpl& operator=(const TransitionListImpl&);
public:
	friend class TransitionProxy;
	friend class TransitionConstProxy;
	string chLabel;
	qList *states; // List of individual states
	EmissionList Emis;
	vector<int> ipEmis;
	explicit TransitionListImpl(
		const string &chLabel,
		qList *states) : chLabel(chLabel), states(states), Emis(this)
	{}
	explicit TransitionListImpl(
		const string &chLabel,
		qList *states,
		size_t size) : chLabel(chLabel), states(states), Emis(this)
	{
		resize(size);
	}
	void resize(size_t newsize);
	void reserve(size_t newsize);
	typedef TransitionProxy::iterator iterator;
	typedef TransitionConstProxy::iterator const_iterator;
	typedef TransitionProxy reference;
	reference operator[](size_t i)
	{
		return TransitionProxy(this,i);
	}
	size_t size(void) const
	{
		return ipCont.size();
	}
	void pop_back(void)
	{
		resize(size()-1);
	}
	iterator begin(void)
	{
		return iterator(this,0);
	}
	const_iterator begin(void) const
	{
		return const_iterator(this,0);
	}
	iterator end(void)
	{
		return iterator(this,size());
	}
	const_iterator end(void) const
	{
		return const_iterator(this,size());
	}
	void push_back(const TransitionProxy &tr)
	{
		int newsize=size()+1;
		resize(newsize);
		(*this)[newsize-1].copy(tr);
	}
	const TransitionProxy back(void)
	{
		return *(end()-1);
	}
	realnum width() const
	{
		return states->width();
	}
};

class TransitionList
{
	// Internal vectors all need to be sized consistently (see three
	// functions below)
	shared_ptr<TransitionListImpl> p_impl;
public:
	typedef TransitionProxy::iterator iterator;
	typedef TransitionConstProxy::iterator const_iterator;
	explicit TransitionList(const string &chLabel, qList *states, size_t size=0)
		: p_impl(new TransitionListImpl(chLabel, states, size))
	{}
	void resize(size_t newsize)
	{
		p_impl->resize(newsize);
	}
	void reserve(size_t newsize)
	{
		p_impl->reserve(newsize);
	}
	TransitionProxy operator[](size_t i)
	{
		return (*p_impl)[i];
	}
	size_t size(void) const
	{
		return p_impl->size();
	}
	void pop_back(void)
	{
		p_impl->pop_back();
	}
	iterator begin(void)
	{
		return p_impl->begin();
	}
	iterator end(void)
	{
		return p_impl->end();
	}
	void push_back(const TransitionProxy &tr)
	{
		p_impl->push_back(tr);
	}
	const TransitionProxy back(void)
	{
		return p_impl->back();
	}
	string &chLabel()
	{
		return p_impl->chLabel;
	}
	qList *&states()
	{
		return p_impl->states;
	}
	EmissionList &Emis()
	{
		return p_impl->Emis;
	}
	vector<int> &ipEmis()
	{
		return p_impl->ipEmis;
	}
	bool isSame (const TransitionListImpl *other) const
	{
		return p_impl.get() == other;
	}
	realnum width() const
	{
		return p_impl->width();
	}
};

inline bool TransitionProxy::systemIs(const TransitionList *query) const
{
	return query->isSame(m_list);
}

// Must include all internal vector elements in these three functions 
inline void TransitionListImpl::resize(size_t newsize)
{
	ipLo.resize(newsize);
	ipHi.resize(newsize);
	ipCont.resize(newsize);
	Coll.resize(newsize);
	EnergyWN.resize(newsize);
	WLAng.resize(newsize);
	ipEmis.resize(newsize,-1);
	chComment.resize(newsize);
}
inline void TransitionListImpl::reserve(size_t newsize)
{
	ipLo.reserve(newsize);
	ipHi.reserve(newsize);
	ipCont.reserve(newsize);
	Coll.reserve(newsize);
	EnergyWN.reserve(newsize);
	WLAng.reserve(newsize);
	ipEmis.reserve(newsize);
}
inline void TransitionProxy::copy(const TransitionProxy& other) const
{
	m_list->ipLo[m_index] = other.m_list->ipLo[other.m_index];
	m_list->ipHi[m_index] = other.m_list->ipHi[other.m_index];
	m_list->ipCont[m_index] = other.m_list->ipCont[other.m_index];
	m_list->Coll[m_index].copy(other.m_list->Coll[other.m_index]);
	m_list->EnergyWN[m_index] = other.m_list->EnergyWN[other.m_index];
	m_list->WLAng[m_index] = other.m_list->WLAng[other.m_index];
	if (other.m_list->ipEmis[other.m_index] == -1)
	{
		m_list->ipEmis[m_index] = -1;
	}
	else
	{
		ASSERT (m_list->ipEmis[m_index] == -1);
		AddLine2Stack();
		m_list->Emis[m_list->ipEmis[m_index]].copy( 
			other.m_list->Emis[other.m_list->ipEmis[other.m_index]]);
	}
}
// End of region needing consistency with TransitionListImpl class

// Handle accessors need to see the structure of the TransitionList
inline qList::iterator TransitionProxy::Lo() const
{
	return m_list->states->begin()+m_list->ipLo[m_index];
}
inline qList::iterator TransitionProxy::Hi() const
{
	return m_list->states->begin()+m_list->ipHi[m_index];
}
inline void TransitionProxy::setLo(int ipLo) const
{
	m_list->ipLo[m_index] = ipLo;
}
inline void TransitionProxy::setHi(int ipHi) const
{
	m_list->ipHi[m_index] = ipHi;
}
inline EmissionList::reference TransitionProxy::Emis() const
{
	int ipEmis = m_list->ipEmis[m_index];
	if (ipEmis == -1)
		return DummyEmis;
	else
		return m_list->Emis[ipEmis];
}
inline int& TransitionProxy::ipEmis() const
{
	return m_list->ipEmis[m_index];
}
inline int TransitionConstProxy::ipEmis() const
{
	return m_list->ipEmis[m_index];
}
inline CollisionProxy TransitionProxy::Coll() const
{
	return m_list->Coll[m_index];
}
/** wavelength, usually in Angstroms, used for printout, can be any units */
inline  realnum &TransitionProxy::WLAng() const
{
	return m_list->WLAng[m_index];
}
inline  realnum TransitionConstProxy::WLAng() const
{
	return m_list->WLAng[m_index];
}
/** transition energy in wavenumbers */
inline realnum &TransitionProxy::EnergyWN() const
{
	return m_list->EnergyWN[m_index];
}
inline realnum TransitionConstProxy::EnergyWN() const
{
	return m_list->EnergyWN[m_index];
}
/** index for line within continuum array,
 * this is on the f, not c, scale,
 * negative ipCont means this is not a radiative transition, 
 * and is used as a sentnecl */
inline long &TransitionProxy::ipCont() const
{
	return m_list->ipCont[m_index];
}
inline long TransitionConstProxy::ipCont() const
{
	return m_list->ipCont[m_index];
}
inline int &TransitionProxy::ipLo() const
{
	return m_list->ipLo[m_index];
}
inline int TransitionConstProxy::ipLo() const
{
	return m_list->ipLo[m_index];
}
inline int &TransitionProxy::ipHi() const
{
	return m_list->ipHi[m_index];
}
inline int TransitionConstProxy::ipHi() const
{
	return m_list->ipHi[m_index];
}

inline TransitionProxy EmissionProxy::Tran() const
{
	TransitionProxy t = TransitionProxy(m_list->m_tlist,ipTran());
	t.check();
	return t;
}
inline TransitionConstProxy EmissionConstProxy::Tran() const
{
	TransitionConstProxy t = TransitionConstProxy(m_list->m_tlist,ipTran());
	t.check();
	return t;
}

inline void TransitionProxy::setComment( const string &comment ) const
{
	m_list->chComment[m_index] = comment;
}
inline string &TransitionProxy::getComment() const
{
	return m_list->chComment[m_index];
}
inline string TransitionConstProxy::getComment() const
{
	return m_list->chComment[m_index];
}

inline realnum TransitionProxy::width() const
{
	return m_list->width();
}

/** enter lines into the line storage array, called once per zone for each line
\param xInten xInten - local emissivity per unit vol, no fill fac
\param wavelength lam integer wavelength
\param *chLab string label for ion
\param chInfo character type of entry for line - 'c' cooling, 'h' heating, 'i' info only, 'r' recom line
\param *chComment string explaining line 
*/

/**PutLine enter local line intensity into the intensity stack for eventual printout 
\param *t transition structure for line
\param *chComment a description of the line
*/
void PutLine(const TransitionProxy &t, const char *chComment);

/**PutLine enter local line intensity into the intensity stack for eventual printout 
\param *t transition structure for line
\param *chComment a description of the line
\param *chLabel the line label
*/
void PutLine(const TransitionProxy &t, const char *chComment, const char *chLabel);

void PutLine(const TransitionProxy& t, const char *chComment, const char *chLabel, const ExtraInten& extra);

/**TexcLine derive excitation temperature of line from contents of line array 
\param *t
*/
double TexcLine(const TransitionProxy &t);

/**DumpLine print various information about an emission line vector, used in debugging 
\param *t
*/
void DumpLine(const TransitionProxy &t);

/** returns fraction of populations that produce emission
\param *t
*/
double emit_frac(const TransitionProxy &t);

double GetLineRec(
	/* this is the number of the emission line in the stack of lines, on the C scale */
	long int ip,
	/* the multiplet wavelength */
	long int lWl);

/** generate null terminated line label from contents of line trans array 
\param *t
*/
string chIonLbl(const TransitionProxy &t);
string chIonLbl(const long& nelem, const long& IonStg);

/**chLineLbl use information in line transfer arrays to generate a line label<BR>
 this label is null terminated 
 \param *t
 */
inline string chLineLbl(const TransitionProxy &t)
{
	return t.chLabel();
}

/**PutCS enter a collision strength into an individual line struc 
\param cs
\param *t  the line struc 
*/
void PutCS(double cs, 
	   const TransitionProxy & t);

/**GenerateTransitionConfiguration - given const TransitionList::iterator &t, writes a label
 * t->Lo->chConfig() - t->Hi->chConfig() (i.e., 2^3S - 2^3P)
 * \param t	transition
 * \return	comment string
 */
string GenerateTransitionConfiguration( const TransitionProxy &t );

/**OccupationNumberLine - derive the photon occupation number at line center for any line 
\param *t
*/
double OccupationNumberLine(const TransitionProxy &t);

/**PutExtra enter and 'extra' intensity source for some line 
\param Extra
*/
class ExtraInten
{
public:
	double v;
	ExtraInten() : v(0.0) {}
	explicit ExtraInten( double extra ) : v(extra) {}
};

/** convert down coll rate back into electron cs in case other parts of code need this for reference 
\param *t - line struct collision strength is stored in t->cs 
\param rate - deexcitation rate, units s-1 
*/
void LineConvRate2CS( const TransitionProxy & t , realnum rate );

/**lgTauGood returns true is we have good (positive) outward optical depths
 * not true if we have overrun optical depth scale from previous iteration
\param *t
*/
inline bool lgTauGood( const TransitionProxy& t )
{
	bool lgOverrunOK = true;
	if ( lgOverrunOK )
	{
		return true;
	}
	else
	{
		// first iteration only use inward optical depths so scale good
		return ( iteration == 1 || 
					// maser - optical depths also ok (but bizarre) 
					t.Emis().TauIn() <= 0. || 
					// TauIn < TauTot means outward optical depth is positive, so OK
					t.Emis().TauIn() < t.Emis().TauTot() );
	}
}

/**MakeCS compute collision strength by g-bar approximations 
\param *t
*/
void MakeCS(const TransitionProxy & t );

extern map<std::string,std::vector<TransitionProxy> > blends;
typedef map<std::string,std::vector<TransitionProxy> >::iterator blend_iterator;

inline double phots( const TransitionProxy &t )
{
	return t.Emis().xIntensity() / t.EnergyErg();
}

#endif // _TRANSITION_H_
