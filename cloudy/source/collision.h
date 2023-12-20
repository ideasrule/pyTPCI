/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef COLLISION_H_
#define COLLISION_H_

#include "global.h"

 /* these are flags for various colliders that are used across the code */
enum collider {
	ipELECTRON,
	ipPROTON,
	ipHE_PLUS,
	ipALPHA,
	ipATOM_H,
	ipATOM_HE,
	ipH2_ORTHO,
	ipH2_PARA,
	ipH2,
	ipNCOLLIDER
};

class t_collider 
{
public:
	long charge;
	double *density;
	realnum mass_amu;

	t_collider()
	{
		charge = LONG_MAX;
		density = NULL;
		mass_amu = FLT_MAX;
	}
};

class ColliderList
{
public:
	vector<t_collider> list;
	ColliderList(const t_dense &d);
	void init();
};
class ColliderDensities
{
	double m_dens[ipNCOLLIDER];
public:
	ColliderDensities(const ColliderList& colls)
	{
		ASSERT( colls.list.size() == ipNCOLLIDER );
		for( unsigned i = 0; i < ipNCOLLIDER; ++i )
		{
			m_dens[i] = (*colls.list[i].density);
		}
	}
	double density(long ipCollider) const
	{
		return m_dens[ipCollider];
	}
};

class collision_rates
{
	double m_rate_coef_ul[ipNCOLLIDER];
public:
	collision_rates()
	{
		for( long i=0; i<ipNCOLLIDER; i++ )
			m_rate_coef_ul[i] = 0.;
	}
	/** collisional de-excitation rate coefficients for individual colliders [cm3 s-1] */
	double *rate_coef_ul_set()
	{
		return m_rate_coef_ul;
	}
	const double *rate_coef_ul() const
	{
		return m_rate_coef_ul;
	}
};

class CollisionList;
class CollisionProxy
{
	CollisionList *m_list;
	int m_index;
public:
	explicit CollisionProxy(CollisionList *list, int index) 
		: m_list(list), m_index(index) {}

	/** [dimensionless] collision strength of rates for transition */
	realnum &col_str() const;
	/** is the collision strength created from gbar */
	int &is_gbar() const;
	/** collisional de-excitation rate coefficients for individual colliders [cm3 s-1] */
	double *rate_coef_ul_set() const;
	const double *rate_coef_ul() const;
	realnum &rate_lu_nontherm_set() const;
	realnum rate_lu_nontherm() const;
	/** cooling and heating due to collisional excitation [erg s-1 cm-3] */
	double &cool() const;
	double &heat() const;

	/** collisional de-excitation rate, [s-1] */
	double ColUL( const ColliderList& colls ) const
	{
		double rate = 0.;
		ASSERT( colls.list.size() == ipNCOLLIDER );
		for( unsigned i = 0; i < colls.list.size(); ++i )
		{
			ASSERT( rate_coef_ul()[i] >= 0.0 );
			rate += rate_coef_ul()[i] * (*colls.list[i].density);
		}
		ASSERT( rate >= 0. );
		return rate;
	}
	double ColUL( const ColliderDensities& colld ) const
	{
		double rate = 0.;
		for( unsigned i = 0; i < ipNCOLLIDER; ++i )
		{
			ASSERT( rate_coef_ul()[i] >= 0.0 );
			rate += rate_coef_ul()[i] * colld.density(i);
		}
		ASSERT( rate >= 0. );
		return rate;
	}
	
	void copy(CollisionProxy other)
	{
		col_str() = other.col_str();
		cool() = other.cool();
		heat() = other.heat();
		for (int i=0; i<ipNCOLLIDER; ++i)
		{
			rate_coef_ul_set()[i] = other.rate_coef_ul()[i];
		}
		rate_lu_nontherm_set() = other.rate_lu_nontherm();
	}		
};

class  CollisionList
{
	vector<collision_rates> m_rates;
	vector<realnum> m_col_str;
	vector<realnum> m_rate_lu_nontherm;
	vector<int> m_is_gbar;
	vector<double> m_cool;
	vector<double> m_heat;
	// DO NOT IMPEMENT
	CollisionList &operator=(const CollisionList&);
	CollisionList(const CollisionList&);
public:
	friend class CollisionProxy;
	typedef CollisionProxy reference;
	explicit CollisionList(size_t i)
	{
		resize(i);
	}
	explicit CollisionList() {}
	CollisionProxy operator[](size_t i)
	{
		return CollisionProxy(this,i);
	}
	size_t size(void) const
	{
		return m_rates.size();
	}
	void resize(size_t i)
	{
		m_rates.resize(i);
		m_rate_lu_nontherm.resize(i);
		m_col_str.resize(i);
		m_is_gbar.resize(i);
		m_cool.resize(i);
		m_heat.resize(i);
	}
	void reserve(size_t i)
	{
		m_rates.reserve(i);
		m_rate_lu_nontherm.reserve(i);
		m_col_str.reserve(i);
		m_is_gbar.reserve(i);
		m_cool.reserve(i);
		m_heat.reserve(i);
	}
};

/** [dimensionless] collision strength of rates for transition */
inline realnum &CollisionProxy::col_str() const
{
	return m_list->m_col_str[m_index];
}
/** is the collision strength created from gbar
 *  -1 = undefined, 0 = Not gbar, 1 = gbar */
inline int &CollisionProxy::is_gbar() const
{
	return m_list->m_is_gbar[m_index];
}
/** collisional de-excitation rate coefficients for individual colliders [cm3 s-1] */
inline double *CollisionProxy::rate_coef_ul_set() const
{
	return m_list->m_rates[m_index].rate_coef_ul_set();
}
inline const double *CollisionProxy::rate_coef_ul() const
{
	return m_list->m_rates[m_index].rate_coef_ul();
}
// nonthermal collisional excitation rate [s-1]
inline realnum &CollisionProxy::rate_lu_nontherm_set() const
{
	return m_list->m_rate_lu_nontherm[m_index];
}
inline realnum CollisionProxy::rate_lu_nontherm() const
{
	return m_list->m_rate_lu_nontherm[m_index];
}
/** cooling and heating due to collisional excitation [erg s-1 cm-3] */
inline double &CollisionProxy::cool() const
{
	return m_list->m_cool[m_index];
}
inline double &CollisionProxy::heat() const
{
	return m_list->m_heat[m_index];
}

/**CollisionJunk set all elements of emission struc to dangerous values  
\param *t
*/
void CollisionJunk( const CollisionProxy & t );

/**CollisionZero set all elements of collision struc to zero 
\param *t
*/
void CollisionZero( const CollisionProxy & t );

#endif // COLLISION_H_
