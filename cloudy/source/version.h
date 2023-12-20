/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VERSION_H_
#define VERSION_H_

/** version.h */

class t_version : public Singleton<t_version> 
{
	friend class Singleton<t_version>;
protected:
	t_version();
public:
	/** date of this version of the code as a string */
	string chDate;

	/** version string of this version of the code */
	string chVersion;

	/** is this a release version?  if so do not print some internal comments */
	bool lgRelease;

	/** information about when and how the code was compiled, including 
	 * compiler version */
	string chInfo;
};

#endif /* VERSION_H_ */
