/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef INPUT_H_
#define INPUT_H_

/* input.h */

#include "module.h"

/** lgIsCommentSeq - is the string s[p] pointing to a comment sequence?
 * if lgReportVisible is true, visible comments will be reported, otherwise not */
bool lgIsCommentSeq( const string& s, size_t p, bool lgReportVisible );

/** lgIsExpungedCommentSeq - does the string s start with an old-style comment? */
bool lgIsExpungedCommentSeq( const string& s );

/** lgInputComment - parse comment - check if argument is comment string, 
 * either upper or lower case -
 * returns true if line is a comment, false if not 
 \param *chLine the input line string
 */
bool lgInputComment( const string& chLine );

/** lgInputEOF - is this line an EOF marker? */
bool lgInputEOF( const string& chLine );

/** StripComment- strips comment part off the command line s
 * if lgStripVisible is false, visible comments are retained
 * hidden comments are always stripped */
void StripComment( string& s, bool lgStripVisible );

/** GetString: retrieve a string between double quotes
 * s : string to be parsed
 * s[p] : place in string where to start parsing, should point to first set of double quotes
 * buf : buffer that will hold the string between double quotes
 * return value : pointer just beyond second set of double quotes, or string::npos in case of failure
 *                (second set of double quotes wasn't found) */
size_t GetString( const string& s, size_t p, string& buf );

/** GetEscape:  This routine is the placeholder for treating character escape sequences.
 * For the moment we treat none. On exit, *p will point one character beyond the escape sequence. */
char GetEscape( const string& s, size_t& p );

/** MakeInputLine: generate input lines for optimizer and grid commands
 * i : index of varied command, should be less than optimize.nvary */
string MakeInputLine(long i);

/** input_readvector: read numbers from the file chFile and store them in a vector */
void input_readvector(const string& chFile, /**< file name to read from */
					  vector<double> &vec); /**< vector - the numbers that were read from the input line(s) */

struct CardInfo {
	/** we will save the original (not caped) image of the line here */
	string chCardSav;
	/** keep track of the include level of each input line
	 * level 0 - main input file
	 * level 1 - init file included from main input file
	 * level 2 - init file included from level 1 init file
	 *   etc... */
	int InclLevel;
	/** is this line visible, or should it be hidden due to the
	 * HIDE keyword and/or the PRINT ON / OFF commands... */
	bool lgVisible;
	CardInfo() : InclLevel(-1), lgVisible(true) {}
};

struct t_input : public module {

	const char *chName() const
	{
		return "input";
	}

	void comment(t_warnings&) {}
	void zero() {}

	/** pointers to structs holding the input lines and associated info */
	vector<CardInfo*> crd;

	/** title entered with the title command */
	string chTitle;

	/** curInclLevel is the current value of the include level stored in CardInfo */
	int curInclLevel;

	/** current visibility status due to the PRINT ON / OFF commands */
	bool lgVisibilityStatus;

	/** this points to the command we are now parsing, within the stack of commands */
	long int nRead;

	/** this is set true if an init file was used in the input deck */
	bool lgInitPresent;

	/** this is set true if underscore present in input stream, which was
	 * set to space */
	bool lgUnderscoreFound;

	/** this is set true if left or right bracket, [ or ], present in input stream, which was
	 * set to space */
	bool lgBracketFound;

	/** set true with no buffering command, used to print comment at end */
	bool lgSetNoBuffering;

	void clear()
	{
		for( size_t i=0; i < crd.size(); ++i )
			delete crd[i];
		crd.clear();
		curInclLevel = 0;
		lgInitPresent = false;
		lgUnderscoreFound = false;
		lgBracketFound = false;
	}

	t_input()
	{
		/* will be set true if no buffering command entered */
		lgSetNoBuffering = false;
	}

	~t_input()
	{
		clear();
	}

private:
	// friend CodeSmell
	friend class Parser;

	/** return the next input command off the command stack 
	 * if this exists, return input line and set lgEOF false,
	 * if none left, return empty string and set lgEOF true.
	 \param *lgEOF true if end of file hit */
	string readarray(bool *lgEOF);
	/** similar to readarray, except that it reads the next
	 * command without actually retrieving it, this is useful
	 * to determine if the next line is a continuation of the
	 * current command */
	string peekarray(bool *lgEOF);

public:
	void echo(FILE *ipOUT) const;

	/** called when 'init' command hit, to reset counters for
	 * placing line images within the storage array */
	void init()
	{
		nRead = -1;
	}
};
extern t_input input;

#endif /* INPUT_H_ */
