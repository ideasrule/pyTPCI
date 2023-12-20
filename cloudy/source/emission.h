/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef EMISSION_H_
#define EMISSION_H_

#include "proxy_iterator.h"
#include "iter_track.h"

class TransitionListImpl;
class TransitionProxy;
class TransitionConstProxy;
class EmissionList;
class EmissionConstProxy;

//typedef iter_track_basic<realnum> tracker;
typedef secant_track tracker;

class EmissionProxy
{
public:
	typedef EmissionList list_type;
	typedef ProxyIterator<EmissionProxy,EmissionConstProxy> iterator;
	typedef ProxyIterator<EmissionConstProxy,EmissionConstProxy> const_iterator;
private:
	friend class ProxyIterator<EmissionProxy,EmissionConstProxy>;
	EmissionList *m_list;
	int m_index;
public:
	explicit EmissionProxy() : m_list(NULL), m_index(-1) {}
	explicit EmissionProxy(EmissionList *list, int index) : 
	m_list(list), m_index(index) {}
	void copy(const EmissionProxy& other);

	void check() const
	{
		ASSERT(ipTran() >= 0);
	}

	/**< all of these are set to dangerous values by EmLineJunk */

	/** type of redistribution function, 
	-1 complete redis Doppler core only, 
	 0 insanity
	 1 incomplete redistribution with damping wings, 
	 2 complete redistribution with damping wings,
	*/
	int &iRedisFun() const;

	/** index for line within fine continuum array, for line overlap */
	long int &ipFine() const;

	/** optical depths - the escape probability is derived only from TauIn
	 * on first iteration.  on later iterations it is the average of esc prob
	 * in inward (TauIn) and outward (TauTot - TauIn) directions */

	/** TauIn - the total inward line optical depth [Napier], in the direction of the 
	 * continuum source.
	 * This is measured from the illuminated face of 
	 * the cloud to the current position on the first iteration, and on second
	 * and later iterations in an open geometry or in an expanding closed geometry.
	 * For a static spherical geometry TauIn is the sum of the optical depth from the
	 * current position to the illuminated face plus the optical depth on the "other
	 * side".
	 * At the illuminated face in a static spherical geometry TauIn is TauTot / 2.
	 * It includes the effects of line overlap. */
	realnum &TauIn() const;

	/** TauInSpecific - the inward optical line depth [Napier], in the direction of
	 *  the continuum source.
	 *  This is much like TauIn(), except that it does not contain contributions
	 *  from line overlap. */
	realnum &TauInSpecific() const;

	/** TauTot - total line optical depth [Napier] through the cloud. 
	 * TauTot is not used on the first iteration since it is not known.  
	 * On second and later iterations in an open or expanding closed geometry 
	 * this is the total optical depth through the computed structure.  
	 * For a static closed geometry this is twice the computed structure since the 
	 * geometry is assumed to be symmetric.
	 *
	 * when the double command is entered TauTot is set to twice the normal value
	 * to simulate the presence of material beyond the computed structure.
	 * */
	realnum &TauTot() const;

	/** TauTrack - track convergence of TauIn / TauTot
	 * This class will detect oscillations and adjust the next estimate for
	 * TauIn / TauTot if this occurs
	 * */
	tracker& TauTrack() const;

	/** TauCon - line optical depth [Napier] to the continuum source from the
	 * illuminated face to the current position.
	 * For an open or expanding closed geometry TauCon is equal to TauIn.
	 * For a static closed geometry TauCon is optical depth from the illuminated face
	 * to the current depth */
	realnum &TauCon() const;

	/** inward fraction [dimensionless] of total line emission*/
	realnum &FracInwd() const;

	/** continuum pumping rate [s-1] from lower to upper level, A*occ num * g_up/g_lo,
	 * this is evaluated in RTMakeStat and RTMakeWind, which are called by HydroPEsc, RT_line_all */
	double &pump() const;

	/** line intensity per unit time and vol [erg s-1 cm-3] */
	double &xIntensity() const;

	/** observed line intensity per unit time and vol [erg s-1 cm-3]; 
	 *  relative to xIntensity(), it includes correction for isotropic background radiation */ 
	double &xObsIntensity() const;

	/** gf value [dimensionless] */
	realnum &gf() const;

	/** escape prob [dimensionless] */
	realnum &Pesc() const;

	/** electron scattering escape prob [dimensionless] */
	realnum &Pelec_esc() const;

	/** destruction probs [dimensionless] */
	realnum &Pdest() const;

	/** total escape prob, from line and continuum scattering [dimensionless] */
	realnum Pesc_total() const { return Pesc() + Pelec_esc(); }

	/** total loss from trapped line -- escape & destruction [dimensionless] */
	realnum Ploss() const { return Pesc_total()+Pdest(); }

	/** damping constant is dampXvel divided by line width 
	 * units are velocity, since becomes dimensionless when div by line width in cm/2
	 * [cm s-1] */
	realnum &dampXvel() const;

	/** [dimensionless] damping constant */
	realnum &damp() const;

	/** [dimensionless] ratio of collisional to radiative excitation, C_lu / ( C_lu + pump )*/
	double &ColOvTot() const;

	/** [dimensionless] branching ratio to auto-ionization, Sum(Aai) / ( Sum(Aul) + Sum(Aai) ) */
	realnum &AutoIonizFrac() const;

	/** atomic constant part of line opacity per atom, 
	 * units: cm^3 / s
	 * divide by line width in cm/s,
	 * to get line center opacity per atom, or absorption cross section, with units cm^2
	 * multiply by PopOpc to get PopOpc/dopper width, the true opacity (cm-1),
	 * then by length to get optical depth */
	realnum &opacity() const;
	// This is the summed opacity in a multiplet.
	double &mult_opac() const;

	/** Population that enters net opacity after correction for stimulated emission [cm-3] */
	double &PopOpc() const;

	/** This variable is the Voigt profile value at line center.
	 *  For lines with damp < 1, the normalized line center profile value is ~(1-damp)
	 *  >>refer	RT	Rutten 2003 (online book) */
	double &VoigtLineCen() const;

	/** transition prob, Einstein A upper to lower [s-1] */
	realnum &Aul() const;

	/** ots rate [cm-3 s-1] */
	double &ots() const;

	int &ipTran() const;

	TransitionProxy Tran() const;
};
class EmissionConstProxy
{
public:
	typedef const EmissionList list_type;
	typedef ProxyIterator<EmissionConstProxy,EmissionConstProxy> iterator;
	typedef ProxyIterator<EmissionConstProxy,EmissionConstProxy> const_iterator;
private:
	friend class ProxyIterator<EmissionConstProxy,EmissionConstProxy>;
	const EmissionList *m_list;
	int m_index;
public:
	explicit EmissionConstProxy() : m_list(NULL), m_index(-1) {}
	explicit EmissionConstProxy(const EmissionList *list, int index) : 
	m_list(list), m_index(index) {}
	void copy(const EmissionConstProxy& other);

	void check() const
	{
		ASSERT(ipTran() >= 0);
	}

	/**< all of these are set to dangerous values by EmLineJunk */

	/** type of redistribution function, 
	-1 complete redis Doppler core only, 
	 0 insanity
	 1 incomplete redistribution with damping wings, 
	 2 complete redistribution with damping wings,
	*/
	int iRedisFun() const;

	/** index for line within fine continuum array, for line overlap */
	long int ipFine() const;

	/** optical depths - the escape probability is derived only from TauIn
	 * on first iteration.  on later iterations it is the average ofesc prob
	 * in inward (TauIn) and outward (TauTot - TauIn) directions */

	/** TauIn - the total inward line optical depth [Napier], in the direction of the 
	 * continuum source.  
	 * This is measured from the illuminated face of 
	 * the cloud to the current position on the first iteration, and on second
	 * and later iterations in an open geometry or in an expanding closed geometry.
	 * For a static spherical geometry TauIn is the sum of the optical depth from the
	 * current position to the illuminated face plus the optical depth on the "other
	 * side".
	 * At the illuminated face in a static spherical geometry TauIn is TauTot / 2.
	 * It includes the effects of line overlap. */
	realnum TauIn() const;

	/** TauInSpecific - the inward optical line depth [Napier], in the direction of
	 *  the continuum source.
	 *  This is much like TauIn(), except that it does not contain contributions
	 *  from line overlap. */
	realnum TauInSpecific() const;

	/** TauTot - total line optical depth [Napier] through the cloud. 
	 * TauTot is not used on the first iteration since it is not known.  
	 * On second and later iterations in an open or expanding closed geometry 
	 * this is the total optical depth through the computed structure.  
	 * For a static closed geometry this is twice the computed structure since the 
	 * geometry is assumed to be symmetric.
	 *
	 * when the double command is entered TauTot is set to twice the normal value
	 * to simulate the presence of material beyond the computed structure.
	 * */
	realnum TauTot() const;

	/** TauTrack - track convergence of TauIn / TauTot
	 * This class will detect oscillations and adjust the next estimate for
	 * TauIn / TauTot if this occurs
	 * */
	const tracker& TauTrack() const;

	/** TauCon - line optical depth [Napier] to the continuum source from the
	 * illuminated face to the current position.
	 * For an open or expanding closed geometry TauCon is equal to TauIn.
	 * For a static closed geometry TauCon is optical depth from the illuminated face
	 * to the current depth */
	realnum TauCon() const;

	/** inward fraction [dimensionless] of total line emission*/
	realnum FracInwd() const;

	/** continuum pumping rate [s-1] from lower to upper level, A*occ num * g_up/g_lo,
	 * this is evaluated in RTMakeStat and RTMakeWind, which are called by HydroPEsc, RT_line_all */
	double pump() const;

	/** line intensity per unit time and vol [erg s-1 cm-3] */
	double xIntensity() const;

	/** observed line intensity per unit time and vol [erg s-1 cm-3]; 
	 *  relative to xIntensity(), it includes correction for isotropic background radiation */ 
	double xObsIntensity() const;

	/** gf value [dimensionless] */
	realnum gf() const;

	/** escape prob [dimensionless] */
	realnum Pesc() const;

	/** electron scattering escape prob [dimensionless] */
	realnum Pelec_esc() const;

	/** destruction probs [dimensionless] */
	realnum Pdest() const;

	/** total escape prob, from line and continuum scattering [dimensionless] */
	realnum Pesc_total() const { return Pesc() + Pelec_esc(); }

	/** total loss from trapped line -- escape & destruction [dimensionless] */
	realnum Ploss() const { return Pesc_total()+Pdest(); }

	/** damping constant is dampXvel divided by line width 
	 * units are velocity, since becomes dimensionless when div by line width in cm/2
	 * [cm s-1] */
	realnum dampXvel() const;

	/** [dimensionless] damping constant */
	realnum damp() const;

	/** [dimensionless] ratio of collisional to radiative excitation, C_lu / ( C_lu + pump )*/
	double ColOvTot() const;

	/** [dimensionless] branching ratio to auto-ionization, Sum(Aai) / ( Sum(Aul) + Sum(Aai) ) */
	realnum AutoIonizFrac() const;

	/** atomic constant part of line opacity per atom, 
	 * units: cm^3 / s
	 * divide by line width in cm/s,
	 * to get line center opacity per atom, or absorption cross section, with units cm^2
	 * multiply by PopOpc to get PopOpc/dopper width, the true opacity (cm-1),
	 * then by length to get optical depth */
	realnum opacity() const;
	// This is the summed opacity in a multiplet.
	double mult_opac() const;

	/** Population that enters net opacity after correction for stimulated emission [cm-3] */
	double PopOpc() const;

	/** This variable is the Voigt profile value at line center.
	 *  For lines with damp < 1, the normalized line center profile value is ~(1-damp)
	 *  >>refer	RT	Rutten 2003 (online book) */
	double VoigtLineCen() const;

	/** transition prob, Einstein A upper to lower [s-1] */
	realnum Aul() const;

	/** ots rate [cm-3 s-1] */
	double ots() const;

	int ipTran() const;

	TransitionConstProxy Tran() const;
};

class EmissionList
{
	TransitionListImpl *m_tlist;
	vector<realnum> m_Aul;
	vector<realnum> m_AutoIonizFrac;
	vector<double> m_ColOvTot;
	vector<realnum> m_damp;
	vector<realnum> m_dampXvel;
	vector<realnum> m_FracInwd;
	vector<realnum> m_gf;
	vector<int> m_iRedisFun;
	vector<long> m_ipFine;
	vector<realnum> m_opacity;
	vector<double> m_mult_opac;
	vector<double> m_ots;
	vector<realnum> m_Pdest;
	vector<realnum> m_Pesc;
	vector<realnum> m_Pelec_esc;
	vector<double> m_PopOpc;
	vector<double> m_VoigtLineCen;
	vector<double> m_pump;
	vector<realnum> m_TauCon;
	vector<realnum> m_TauIn;
	vector<realnum> m_TauInSpecific;
	vector<realnum> m_TauTot;
	vector<tracker > m_TauTrack;
	vector<double> m_xIntensity;
	vector<double> m_xObsIntensity;
	vector<int> m_ipTran;
	friend class EmissionProxy;
	friend class EmissionConstProxy;
public:
	typedef EmissionProxy reference;
	typedef EmissionProxy::iterator iterator;
	typedef EmissionConstProxy::iterator const_iterator;
	explicit EmissionList(TransitionListImpl *tlist, size_t i) : m_tlist(tlist)
	{
		resize(i);
	}
	explicit EmissionList(TransitionListImpl *tlist) : m_tlist(tlist) {}
	reference operator[](size_t i)
	{
		return EmissionProxy(this,i);
	}
	size_t size(void) const
	{
		return m_Aul.size();
	}
	void resize(size_t i);
	iterator begin()
	{
		return iterator(this,0);
	}
	const_iterator begin() const
	{
		return const_iterator(this,0);
	}
	iterator end()
	{
		return iterator(this,size());
	}
	const_iterator end() const
	{
		return const_iterator(this,size());
	}
};

/**EmLineJunk set all elements of emission struc to dangerous values 
\param *t
*/
void EmLineJunk( EmissionList::reference t );
/**EmLineZero set all elements of emission struc to zero 
\param *t
*/
void EmLineZero( EmissionList::reference t );
/**TauZero set initial values of inward, outward, and local-to-continuum-source
 * optical depths
 \param *t
 * */
void TauZero( EmissionList::reference t );

inline void EmissionList::resize(size_t i)
{
	size_t oldsize = m_Aul.size();
	m_Aul.resize(i);
	m_AutoIonizFrac.resize(i);
	m_ColOvTot.resize(i);
	m_damp.resize(i);
	m_dampXvel.resize(i);
	m_gf.resize(i);
	m_FracInwd.resize(i);
	m_ipFine.resize(i);
	m_iRedisFun.resize(i);
	m_ots.resize(i);
	m_opacity.resize(i);
	m_mult_opac.resize(i);
	m_Pdest.resize(i);
	m_Pelec_esc.resize(i);
	m_Pesc.resize(i);
	m_PopOpc.resize(i);
	m_VoigtLineCen.resize(i);
	m_TauCon.resize(i);
	m_TauIn.resize(i);
	m_TauInSpecific.resize(i);
	m_TauTot.resize(i);
	m_TauTrack.resize(i);
	m_pump.resize(i);
	m_xIntensity.resize(i);
	m_xObsIntensity.resize(i);
	m_ipTran.resize(i,-1);
	for (size_t newelem=oldsize; newelem < size(); ++newelem)
	{
		EmLineJunk((*this)[newelem]);
	/* 
		\todo 2 Does doing EmLineZero here defeat the purpose of EmLineJunk? 
		* maybe we should pass full set of Emis components, fill everything in 
		* here, and THEN use EmLineZero?  */
		EmLineZero((*this)[newelem]);
		TauZero((*this)[newelem]);
	}
}

inline int &EmissionProxy::iRedisFun() const
{
	return m_list->m_iRedisFun[m_index];
}

inline int EmissionConstProxy::iRedisFun() const
{
	return m_list->m_iRedisFun[m_index];
}

inline long int &EmissionProxy::ipFine() const
{
	return m_list->m_ipFine[m_index];
}

inline long int EmissionConstProxy::ipFine() const
{
	return m_list->m_ipFine[m_index];
}

inline realnum &EmissionProxy::TauIn() const
{
	return m_list->m_TauIn[m_index];
}

inline realnum EmissionConstProxy::TauIn() const
{
	return m_list->m_TauIn[m_index];
}

inline realnum &EmissionProxy::TauInSpecific() const
{
	return m_list->m_TauInSpecific[m_index];
}

inline realnum EmissionConstProxy::TauInSpecific() const
{
	return m_list->m_TauInSpecific[m_index];
}

inline realnum &EmissionProxy::TauTot() const
{
	return m_list->m_TauTot[m_index];
}

inline realnum EmissionConstProxy::TauTot() const
{
	return m_list->m_TauTot[m_index];
}

inline tracker&EmissionProxy:: TauTrack() const
{
	return m_list->m_TauTrack[m_index];
}

inline const tracker&EmissionConstProxy:: TauTrack() const
{
	return m_list->m_TauTrack[m_index];
}

inline realnum &EmissionProxy::TauCon() const
{
	return m_list->m_TauCon[m_index];
}

inline realnum EmissionConstProxy::TauCon() const
{
	return m_list->m_TauCon[m_index];
}

inline realnum &EmissionProxy::FracInwd() const
{
	return m_list->m_FracInwd[m_index];
}

inline realnum EmissionConstProxy::FracInwd() const
{
	return m_list->m_FracInwd[m_index];
}

inline double &EmissionProxy::pump() const
{
	return m_list->m_pump[m_index];
}

inline double EmissionConstProxy::pump() const
{
	return m_list->m_pump[m_index];
}

inline double &EmissionProxy::xIntensity() const
{
	return m_list->m_xIntensity[m_index];
}

inline double &EmissionProxy::xObsIntensity() const
{
	return m_list->m_xObsIntensity[m_index];
}

inline double EmissionConstProxy::xIntensity() const
{
	return m_list->m_xIntensity[m_index];
}

inline double EmissionConstProxy::xObsIntensity() const
{
	return m_list->m_xObsIntensity[m_index];
}

inline int &EmissionProxy::ipTran() const
{
	return m_list->m_ipTran[m_index];
}

inline int EmissionConstProxy::ipTran() const
{
	return m_list->m_ipTran[m_index];
}

inline realnum &EmissionProxy::gf() const
{
	return m_list->m_gf[m_index];
}

inline realnum EmissionConstProxy::gf() const
{
	return m_list->m_gf[m_index];
}

inline realnum &EmissionProxy::Pesc() const
{
	return m_list->m_Pesc[m_index];
}

inline realnum EmissionConstProxy::Pesc() const
{
	return m_list->m_Pesc[m_index];
}

inline realnum &EmissionProxy::Pelec_esc() const
{
	return m_list->m_Pelec_esc[m_index];
}

inline realnum EmissionConstProxy::Pelec_esc() const
{
	return m_list->m_Pelec_esc[m_index];
}

inline realnum &EmissionProxy::Pdest() const
{
	return m_list->m_Pdest[m_index];
}

inline realnum EmissionConstProxy::Pdest() const
{
	return m_list->m_Pdest[m_index];
}

inline realnum &EmissionProxy::dampXvel() const
{
	return m_list->m_dampXvel[m_index];
}

inline realnum EmissionConstProxy::dampXvel() const
{
	return m_list->m_dampXvel[m_index];
}

inline realnum &EmissionProxy::damp() const
{
	return m_list->m_damp[m_index];
}

inline realnum EmissionConstProxy::damp() const
{
	return m_list->m_damp[m_index];
}

inline double &EmissionProxy::ColOvTot() const
{
	return m_list->m_ColOvTot[m_index];
}

inline double EmissionConstProxy::ColOvTot() const
{
	return m_list->m_ColOvTot[m_index];
}

inline realnum &EmissionProxy::AutoIonizFrac() const
{
	return m_list->m_AutoIonizFrac[m_index];
}

inline realnum EmissionConstProxy::AutoIonizFrac() const
{
	return m_list->m_AutoIonizFrac[m_index];
}

inline realnum &EmissionProxy::opacity() const
{
	return m_list->m_opacity[m_index];
}

inline realnum EmissionConstProxy::opacity() const
{
	return m_list->m_opacity[m_index];
}

inline double &EmissionProxy::mult_opac() const
{
	return m_list->m_mult_opac[m_index];
}

inline double EmissionConstProxy::mult_opac() const
{
	return m_list->m_mult_opac[m_index];
}

inline double &EmissionProxy::PopOpc() const
{
	return m_list->m_PopOpc[m_index];
}

inline double EmissionConstProxy::PopOpc() const
{
	return m_list->m_PopOpc[m_index];
}

inline double &EmissionProxy::VoigtLineCen() const
{
	return m_list->m_VoigtLineCen[m_index];
}

inline double EmissionConstProxy::VoigtLineCen() const
{
	return m_list->m_VoigtLineCen[m_index];
}

inline realnum &EmissionProxy::Aul() const
{
	return m_list->m_Aul[m_index];
}

inline realnum EmissionConstProxy::Aul() const
{
	return m_list->m_Aul[m_index];
}

inline double &EmissionProxy::ots() const
{
	return m_list->m_ots[m_index];
}

inline double EmissionConstProxy::ots() const
{
	return m_list->m_ots[m_index];
}

inline void EmissionProxy::copy(const EmissionProxy& other)
{
	iRedisFun() = other.iRedisFun();
	ipFine() = other.ipFine();
	TauIn() = other.TauIn();
	TauInSpecific() = other.TauInSpecific();
	TauTot() = other.TauTot();
	TauCon() = other.TauCon();
	FracInwd() = other.FracInwd();
	gf() = other.gf();
	Pesc() = other.Pesc();
	Pelec_esc() = other.Pelec_esc();
	Pdest() = other.Pdest();
	dampXvel() = other.dampXvel();
	damp() = other.damp();
	AutoIonizFrac() = other.AutoIonizFrac();
	opacity() = other.opacity();
	mult_opac() = other.mult_opac();
	Aul() = other.Aul();
	TauTrack() = other.TauTrack();
	pump() = other.pump();
	xIntensity() = other.xIntensity();
	xObsIntensity() = other.xObsIntensity();
	ColOvTot() = other.ColOvTot();
	PopOpc() = other.PopOpc();
	VoigtLineCen() = other.VoigtLineCen();
	ots() = other.ots();
	ipTran() = other.ipTran();
}

#endif // EMISSION_H_
