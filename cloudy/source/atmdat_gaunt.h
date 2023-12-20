/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ATMDAT_GAUNT_H
#define ATMDAT_GAUNT_H

#include "container_classes.h"
#include "vectorize.h"

struct t_gff
{
	double Te_used;
	long nflux_used;
	vector_avx<realnum> gff;
	t_gff() : Te_used(-1.), nflux_used(0) {}
};

struct t_brems_sum
{
	double Te_used;
	long nlo;
	long nhi;
	double brems_sum;
	t_brems_sum() : Te_used(-1.), nlo(0), nhi(0), brems_sum(0.) {}
};

struct t_brems_vec
{
	double Te_used;
	long nhi;
	vector<double> brems_vec;
	t_brems_vec() : Te_used(-1.), nhi(0) {}
};

struct t_brems_den
{
	double den_Hm;
	double den_Hp;
	double den_Hep;
	double den_Hepp;
	// this will exclude hydrogen and helium
	double den_ion[LIMELM+1];
};

class t_gaunt : public Singleton<t_gaunt>
{
	friend class Singleton<t_gaunt>;

	size_t p_np_gam2;
	size_t p_np_u;

	double p_lg_gam2_min;
	double p_lg_gam2_max;
	double p_lg_u_min;
	double p_lg_u_max;
	double p_step;

	vector<double> p_lg_gam2;
	vector<double> p_lg_u;
	multi_arr<double,4> p_gff;
	vector_avx<realnum> p_vexp_arg;

	t_gff p_gff_ion[LIMELM+1];
	t_brems_sum p_cache[LIMELM+1];
	t_brems_sum p_hminus;
	t_brems_vec p_cache_vec[LIMELM+1];
	t_brems_vec p_hminus_vec;

	void p_read_table(const char* fnam, long nelem);
	void p_gauntff_vec_sub(long Z, double Te, const double anulog10[], long nmin, long nmax);
	void p_gauntff_vec(long Z, double Te, const double anulog10[], long nflux);
	void p_setup_brems(long ion, double Te);
protected:
	t_gaunt() {}
public:
	double gauntff(long Z, double Te, double anu);
	double brems_cool(long ion, double Te);
	void brems_rt(long ion, double Te, double abun, vector<double>& arr);
	void brems_opac(long ion, double Te, double abun, vector<double>& arr);
	void brems_sum_ions(t_brems_den& sum) const;
};

#endif
