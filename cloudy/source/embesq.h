/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef EMBESQ_H_
#define EMBESQ_H_

/** embesq.h */
struct t_embesq {
	realnum em13C1910,
	  em1218, 
	  em1486,
	  em2669, 
	  em1895;
	};
extern t_embesq embesq;


#endif /* EMBESQ_H_ */
