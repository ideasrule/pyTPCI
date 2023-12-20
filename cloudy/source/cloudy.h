/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CLOUDY_H_
#define CLOUDY_H_


/**cloudy the main routine, this IS Cloudy,
  \return 0 for no-crash, 1 for bad start */
bool cloudy(void);

/** SanityCheck confirm that various parts of the code still work
 * \param *chJob is either "begin" or "final"<BR>
 * begin is before code starts up<BR>
 * final is after model is complete */
void SanityCheck(const char *chJob);


#endif /* CLOUDY_H_ */
