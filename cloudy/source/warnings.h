/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef WARNINGS_H_
#define WARNINGS_H_

/* warnings.h */

#include "module.h"

class t_warnings : public module {
public:
	const char *chName() const
	{
		return "warnings";
	}

	/**wcnint initialize stack or warnings, cautions, notes */
	void zero(void);
	void comment(t_warnings&) {}

	/** a comment about the geometry after model stops */
	vector<string> chRgcln;

	/** these are the strings that contain the warnings, cautions,
	 * and notes about the calculation */
	vector<string> chWarnln; 
	vector<string> chCaunln;
	vector<string> chBangln;
	vector<string> chNoteln;

	/** flags set if warnings or cautions present */
	bool lgWarngs;
	bool lgCautns;

	/**rgcin enter comment about the geometry after model stops
	   \param *chLine comment to be printed
	*/
	void rgcin(const string& chLine) { chRgcln.emplace_back( chLine ); }

	/**warnin enter warnings at the end of the calculations into large stack 
	   \param *chLine warning to be printed
	*/
	void warnin(const string& chLine) { lgWarngs = true; chWarnln.emplace_back( chLine ); }
	
	/**caunin called by comment to enter caution into comment stack 
	   \param *chLine caution to be printed
	*/
	void caunin(const string& chLine) { lgCautns = true; chCaunln.emplace_back( chLine ); }
	
	/**bangin called by routine comment to enter surprise into comment stack 
	   \param *chLine surprise to be printed
	*/
	void bangin(const string& chLine) { chBangln.emplace_back( chLine ); }

	/**notein enter a note about calculation into comment array 
	   \param *chLine note to be printed
	*/
	void notein(const string& chLine) { chNoteln.emplace_back( chLine ); }
};
extern t_warnings warnings;


/**rgcin enter comment about the geometry after model stops
   \param *chLine comment to be printed
*/
inline void rgcin(const string& chLine)
{
	warnings.rgcin(chLine);
}

/**warnin enter warnings at the end of the calculations into large stack 
   \param *chLine warning to be printed
*/
inline void warnin(const string& chLine)
{
	warnings.warnin(chLine);
}

/**caunin called by comment to enter caution into comment stack 
   \param *chLine caution to be printed
*/
inline void caunin(const string& chLine)
{
	warnings.caunin(chLine);
}

/**bangin called by routine comment to enter surprise into comment stack 
   \param *chLine surprise to be printed
*/
inline void bangin(const string& chLine)
{
	warnings.bangin(chLine);
}

/**notein enter a note about calculation into comment array 
   \param *chLine note to be printed
*/
inline void notein(const string& chLine)
{
	warnings.notein(chLine);
}

#endif /* WARNINGS_H_ */
