/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef FUDGEC_H_
#define FUDGEC_H_

/* fudgec.h */

#define	NFUDGC	10
struct t_fudgec {
	/**
	 *parameters set with fudge command
	 *used to pass parameters to temporary parts of code
	 */
	realnum fudgea[NFUDGC];
	/** number of parameters */
	long int nfudge;
	/** set true if fudge ever used */
	bool lgFudgeUsed;
	};
extern t_fudgec fudgec;


#endif /* FUDGEC_H_ */
