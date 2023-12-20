/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef ATMDAT_HYDR_TP_H
#define ATMDAT_HYDR_TP_H

#include "container_classes.h"

class t_hydro_tbl : public Singleton<t_hydro_tbl>
{
	friend class Singleton<t_hydro_tbl>;

	// The n -> n' and n,l -> n',l' files are much too big to be read in
	// its entirety for each sim. The vast majority of the data would
	// never be used. To reduce I/O time and memory use we initially read
	// data with 0 <= n_up < p_stride. If data for higher levels are
	// needed, we read the next chunk: p_stride <= n_up < 2*n_stride, etc.
	// Each chunk is kept in a separate data structure of type arrn or
	// arrnl. The pointers to these structures are kept in the vectors
	// p_tpn and p_tpnl. Pointers into the ascii files are kept in p_posn
	// and p_posnl so that for each chunk we can immediately start reading
	// in the correct position in the file.

	static const size_t p_stride = 64;

	size_t p_nmaxn;
	size_t p_nmaxn_read;
	long p_posn;
	typedef multi_arr<realnum,2> arrn;
	vector<arrn*> p_tpn;

	size_t p_nmaxnl_l;
	size_t p_nmaxnl_u;
	size_t p_nmaxnl_read;
	long p_posnl;
	typedef multi_arr<realnum,4> arrnl;
	vector<arrnl*> p_tpnl;

	size_t p_Zmax;
	size_t p_nmaxnn;
	multi_arr<realnum,3> p_tpnn;
	multi_arr<double,3> p_wnnn;

	size_t p_nmaxcs;
	size_t p_nenrgs;
	multi_arr<double,3> p_en;
	multi_arr<double,3> p_cs;

	void p_initn(long nu);
	void p_initnl(long nu);
	void p_initnn();
	void p_initcs();

	realnum p_RM(long Z) const;
protected:
	t_hydro_tbl()
	{
		p_nmaxn = 0;
		p_nmaxnl_l = 0;
		p_nmaxnl_u = 0;
	}
	~t_hydro_tbl()
	{
		for( size_t i=0; i < p_tpn.size(); ++i )
			delete p_tpn[i];
		for( size_t i=0; i < p_tpnl.size(); ++i )
			delete p_tpnl[i];
	}
public:
	size_t nmaxn() { p_initn(-1); return p_nmaxn; }
	size_t nmaxnl_l() { p_initnl(-1); return p_nmaxnl_l; }
	size_t nmaxnl_u() { p_initnl(-1); return p_nmaxnl_u; }
	size_t nmaxnn() { p_initnn(); return p_nmaxnn; }
	size_t nmaxcs() { p_initcs(); return p_nmaxcs; }

	realnum tp(long nl, long nu, long Z);
	realnum tp(long nl, long ll, long nu, long lu, long Z);
	realnum tp(long n, long ll, long lu, long Z);
	double wn(long n, long ll, long lu, long Z);
	double cs(double e, long n, long l, long Z);
};

#endif
