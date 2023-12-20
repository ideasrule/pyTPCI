/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef NOEXEC_H_
#define NOEXEC_H_

/* noexec.h */

/**
 *flag saying no to compute model, just generate continuum and stop
 *set with cdNoex call
 */
struct t_noexec {
	bool lgNoExec;
	};
extern t_noexec noexec;



#endif /* NOEXEC_H_ */
