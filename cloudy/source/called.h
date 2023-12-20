/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CALLED_H_
#define CALLED_H_

/** called.h */
struct t_called {

	/** should we print info?  usually true, set false with
	 * "print off" or "print quiet" commands */
	bool lgTalk;	

	/** saves inital value of lgTalk in case it needs to be reset */
	bool lgTalkSave;

	/** this is set true when cdTalk is called with false,
	 * means do not ever turn printout on again */
	bool lgTalkForcedOff;

	/** set true in cdInit, optimize_do sets false when 
	 * we want to make sure that print on does not turn on print */
	bool lgTalkIsOK;

};

extern t_called called;

#endif /* CALLED_H_ */
