/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PARSER_H_
#define PARSER_H_

 /**nWord determine whether match to a keyword occurs on command line,
   return value is 0 if no match, and position of match within string if hit 
	  \param *chKey
	  \param *chCard
 */ 

#include "service.h"
#include "thirdparty.h"

const char * nWord(const char *chKey, 
	    const char *chCard);

class Parser;
class DepthTable;
class LineID;

typedef void (*OptionParser)(Parser &);

struct CloudyCommand {
	const char *name;
	OptionParser action;
};

bool isBoundaryChar(char c);

class Symbol {
public:
	enum tokens { ERROR, EOSTAT, NUMBER, NAME, STRING, OPERATOR, NTOKS };
	enum tokens toktype;
	string value;
	Symbol(enum tokens t, const string& v) : toktype(t), value(v) {}
};

/** Parser class holds pointer to string currently being analysed */
class Parser
{
	string m_card;         // all-caps version of m_card_raw
	string m_card_raw;     // raw input line with all comments stripped
	string m_card_comment; // input line with visible comments included, used for output
	long int m_len;
	long int m_off;
	bool m_lgEOL;
	const CloudyCommand * const m_Commands;
	std::map<string,double> m_symtab;
	std::map<string,size_t> m_uniqueLen;
public:
	long int m_nqh;
	bool m_lgDSet, m_lgEOF;

	explicit Parser(void) :	m_Commands(NULL)
	{
		init();
	}
	explicit Parser(const CloudyCommand *commands) : m_Commands(commands) 
	{
		init();
	}
private:
	void init()
	{
		m_nqh = 0;
		m_lgDSet = m_lgEOF = false;
		setline("");
		if( m_Commands != NULL )
		{
			trieNode* root = new trieNode;
			for( size_t i=0; m_Commands[i].name != NULL; ++i )
				insertToken(root, m_Commands[i].name);
			for( size_t i=0; m_Commands[i].name != NULL; ++i )
			{
				size_t ul = findUniqueLen(root, m_Commands[i].name);
				if( ul < 4 )
					ul = min(strlen(m_Commands[i].name),4);
				m_uniqueLen[m_Commands[i].name] = ul;
			}
			delete root;
		}
	}
	void newlineProcess();
	bool at_end() const;
	char current( void ) const
	{
		return m_card[m_off];
	}
	char current_raw() const;
	void skip_whitespace();
	string m_getCommandToken() const;
	size_t m_getUniqueLen(const string& s2) const
	{
		size_t uniqueLen;
		auto p = m_uniqueLen.find(s2);
		if( p == m_uniqueLen.end() )
			uniqueLen = s2.length();
		else
			uniqueLen = p->second;
		return uniqueLen;
	}
public:
	bool getline();
	void setline(const string& card);
	void set_point(long int ipnt)
	{
		m_off = ipnt;
	}
	const char * nWord(const char *chKey) const;
	bool lgReachedEnd();
	void showLocation(FILE *io = ioQQQ) const;
public:
	long int GetElem( void ) const;
	double FFmtRead( void );
	double getNumberPlain( const char *chDesc );
	double getNumberCheck( const char *chDesc );
	double getNumberDefault( const char *chDesc, double fdef );
	double getNumberCheckLogLinNegImplLog( const char *chDesc );
	double getNumberCheckAlwaysLog( const char *chDesc );
	double getNumberCheckAlwaysLogLim( const char *chDesc, double flim );
	double getNumberDefaultAlwaysLog( const char *chDesc, double fdef );
	double getNumberDefaultNegImplLog( const char *chDesc, double fdef );
	bool lgEOL(void) const
	{
		return m_lgEOL;
	}
	void setEOL(bool val)
	{
		m_lgEOL = val;
	}
	NORETURN void NoNumb(const char *chDesc) const;
private:
	int nMatch1(const char *chKey) const
	{
		const char *p=chKey;

		while (isspace(*p))
			++p;

		for (const char *q=p; *q; ++q)
			ASSERT(!islower(*q));

		if ( !isBoundaryChar(*p))
		{
			const char *c = m_card.c_str();
			const char *q = ::nWord(p, c);
			if (NULL == q)
				return 0;
			else
				return q-c+1;
		}
		else
		{
			// If the keyword starts with a member of the boundary character
			// set, can't require it to be preceded by one so revert to explicit
			// matching
			return ::nMatch(chKey, m_card.c_str());
		}
	}
public:
	bool nMatch(const char *chKey) const
	{
		return nMatch1(chKey) != 0;
	}
	bool GetParam(const char *chKey, double *val)
	{
		int i = nMatch1(chKey);
		if (i > 0) {
			m_off = i-1;
			*val = FFmtRead();
		}
		return i>0;
	}
	bool GetRange(const char *chKey, double *val1, double *val2)
	{
		int i = nMatch1(chKey);
		if (i > 0) {
			m_off = i-1;
			*val1 = FFmtRead();
			*val2 = FFmtRead();
		}
		return i>0;
	}
	bool nMatchErase(const char *chKey)
	{
		const char *p=chKey;
		while (isspace(*p))
			++p;
		int i = nMatch1(p);
		bool found = (i != 0);
		if(found) 
		{
			const long len = strlen(p);
			/* erase this keyword, it upsets FFmtRead */
			for (long j=0; j<len; ++j)
			{
				m_card[i+j-1] = ' ';
			}
		}		
		return found;
	}
	bool hasCommand(const string& s2);
	bool peekNextCommand(const string& s2);
	bool Command(const char *name, OptionParser doOpts)
	{
		bool lgFound = hasCommand(name);
		if ( lgFound )
			(*doOpts)(*this);
		return lgFound;
	}
	string ClosestMatch(const string& token) const;
	bool isComment(void) const;
	bool isVar(void) const;
	std::string getVarName(void);
	void doSetVar(void);
	void echo(void) const;
	bool last(void) const;
	int PrintLine(FILE *fp) const
	{
		return fprintf( fp, " ==%-.80s==\n", m_card_comment.c_str());
	}
	NORETURN void CommandError( void ) const;
	NORETURN void Error( const char *msg ) const;
	NORETURN void StringError( ) const;
	int GetQuote( string& chLabel );
	const char *StandardEnergyUnit(void) const;
	string StandardFluxUnit(void) const;
	string getFirstChunk(long i);
	string getFirstChunkRaw(long i);
	string getRawTail()
	{
		return m_card_raw.substr(m_off);
	}
	void help(FILE *fp) const;
	double getWave();
	double getWaveOpt();
	LineID getLineID(bool lgAtStart=true);
	Symbol getSymbol();
	int getElement();
	void getPairs(vector<double>& a, vector<double> & b);
	void readList(vector<string>& list, const char *chName);
	void readLaw(DepthTable& table);
};

// helper functions for the DataParser class, do not call these directly
template<typename T>
inline void getTokenOptionalImpl(istringstream& iss, const string&, T& var)
{
	iss >> var;
	if( iss.fail() )
		var = T();
}

// optimized specializations for the most common types
template<>
inline void getTokenOptionalImpl(istringstream& iss, const string& s, double& var)
{
	FPRead(iss, s, var);
}

template<>
inline void getTokenOptionalImpl(istringstream& iss, const string& s, sys_float& var)
{
	FPRead(iss, s, var);
}

template<>
inline void getTokenOptionalImpl(istringstream& iss, const string& s, long long& var)
{
	IntRead(iss, s, var);
}

template<>
inline void getTokenOptionalImpl(istringstream& iss, const string& s, unsigned long long& var)
{
	IntRead(iss, s, var);
}

template<>
inline void getTokenOptionalImpl(istringstream& iss, const string& s, long& var)
{
	IntRead(iss, s, var);
}

template<>
inline void getTokenOptionalImpl(istringstream& iss, const string& s, unsigned long& var)
{
	IntRead(iss, s, var);
}

template<>
inline void getTokenOptionalImpl(istringstream& iss, const string& s, int& var)
{
	IntRead(iss, s, var);
}

template<>
inline void getTokenOptionalImpl(istringstream& iss, const string& s, unsigned int& var)
{
	IntRead(iss, s, var);
}

//! ES_NONE means that neither blank lines nor a field of stars are end-of-data (EOD)
//!         markers a blank line is considered a comment and a field of stars is ignored
//! ES_STARS_ONLY means that a field of stars is an EOD marker
//!               a blank line is considered a comment
//! ES_STARS_AND_BLANKS means that both blank lines and a field of stars are EOD markers

enum eod_style { ES_INVALID, ES_NONE, ES_STARS_ONLY, ES_STARS_AND_BLANKS };

class DataParser {
	string p_filename;   //! the name of the data file
	fstream p_io;        //! stream for reading data file
	eod_style p_es;      //! what are the allowed EOD markers?
	string p_line;       //! the current line being read
	size_t p_nr;         //! number of the line we are parsing
	istringstream p_ls;  //! stream for reading current line
	bool p_lgEOF;        //! have we passed beyond the EOF?

	void p_open(const string& name, eod_style es, access_scheme as);
	void p_close();
	// prepare a recently obtained line for parsing
	void p_newlineProcess();
	// returns true if the current line is concidered comment only
	bool p_isComment() const
	{
		if( p_es == ES_NONE || p_es == ES_STARS_ONLY )
			return ( p_blankLine() || p_line[0] == '#' );
		else if( p_es == ES_STARS_AND_BLANKS )
			return ( p_line[0] == '#' );
		else
			TotalInsanity();
	}
	// returns true if the line contains only whitespace, possibly followed by a comment
	// NB NB -- this routine returns false if the comment character is in column 0...
	// this is needed because comments starting in column 0 may have a different meaning
	// than blank lines, depending on the value of p_es...
	bool p_blankLine() const
	{
		for( size_t i=0; i < p_line.length(); ++i )
		{
			if( i > 0 && p_line[i] == '#' )
				return true;
			else if( !isspace(p_line[i]) )
				return false;
		}
		return true;
	}
	// get the position on the current line
	size_t p_pos()
	{
		long p = p_ls.tellg();
		if( p < 0 )
			return p_line.length();
		else
			return size_t(p);
	}
	// set the position on the current line (i.e. skip to position p)
	void p_pos(size_t p)
	{
		p_ls.seekg(p);
		if( p >= p_line.length() )
			p_ls.setstate(ios_base::eofbit);
	}
	// skip whitespace on the current line
	void p_skipWS()
	{
		while( p_ls.good() )
		{
			char c = p_ls.peek();
			if( !isspace(c) )
				break;
			(void)p_ls.get();
		}
	}
	bool p_isSeparator(char c)
	{
		return ( c == '=' || c == ',' );
	}
	// replace separators with spaces
	void p_replaceSep()
	{
		bool lgModified = false;
		for( size_t i=0; i < p_line.length(); ++i )
			if( p_isSeparator(p_line[i]) )
			{
				p_line[i] = ' ';
				lgModified = true;
			}
		if( lgModified )
			p_ls.str(p_line);
	}
	// helper routine showing where parsing stopped
	void p_showLocation(size_t p, FILE *io);
	// this implements reading a quoted string
	void p_getQuoteOptional(string& str);
public:
	// default constructor
	DataParser() : p_es(ES_INVALID), p_nr(0), p_lgEOF(false) {}	
	// constructor
	// name: the name of the data file to be read
	// es: which inline EOD markers are allowed? (see also above)
	// as: access scheme for searching the data file (see cpu.h)
	DataParser(const string& name, eod_style es, access_scheme as = AS_DEFAULT)
	{
		p_open(name, es, as);
	}
	// open data file for parsing, normally done via ctor
	void open(const string& name, eod_style es, access_scheme as = AS_DEFAULT)
	{
		p_close();
		p_open(name, es, as);
	}
	// close data file, normally done via dtor
	void close()
	{
		p_close();
	}
	// check whether opening the file succeeded, only relevant for optional access schemes
	bool isOpen() const { return p_io.is_open(); }
	// check magic number, aborts if mismatch is found
	void checkMagic(long i0);
	void checkMagic(long i0, long i1);
	void checkMagic(long i0, long i1, long i2);
	void checkMagic(long i0, long i1, long i2, long i3);
	// read next non-comment line, returns false on EOF or failure
	bool getline();
	// feed a data line to the parser
	void setline(const string& line)
	{
		p_line = line;
		p_newlineProcess();
	}
	// try to read optional token, returns false if operation failed
	template<typename T>
	bool getTokenOptional(T& var)
	{
		getTokenOptionalImpl(p_ls, p_line, var);
		if( p_ls.fail() )
			return false;
		else if( !isspace(p_ls.peek()) && !p_ls.eof() )
			errorAbort("found trailing junk after token");
		else
			return true;
	}
	// try to read vector of optional tokens, returns number of tokens read successfully
	template<typename T>
	size_t getTokenOptional(T var[], size_t n)
	{
		static_assert( !SameType<T,char>::value &&
					   !SameType<T,signed char>:: value &&
					   !SameType<T,unsigned char>:: value,
					   "Reading C-style strings is not supported, use C++ strings instead." );
		for( size_t i=0; i < n; ++i )
		{
			if( !getTokenOptional(var[i]) )
				return i;
		}
		return n;
	}
	// read token, will abort on failure
	template<typename T>
	void getToken(T& var)
	{
		if( !getTokenOptional(var) )
			errorAbort("failed to read valid data");
	}
	// read vector of tokens, will abort on failure
	template<typename T>
	void getToken(T var[], size_t n)
	{
		static_assert( !SameType<T,char>::value &&
					   !SameType<T,signed char>:: value &&
					   !SameType<T,unsigned char>:: value,
					   "Reading C-style strings is not supported, use C++ strings instead." );
		for( size_t i=0; i < n; ++i )
			getToken(var[i]);
	}
	// read an optional quoted string
	bool getQuoteOptional(string& str)
	{
		p_getQuoteOptional(str);
		if( p_ls.fail() )
			return false;
		else if( !isspace(p_ls.peek()) && !p_ls.eof() )
			errorAbort("found trailing junk after token");
		else
			return true;
	}
	// read a quoted string
	void getQuote(string& str)
	{
		if( !getQuoteOptional(str) )
			errorAbort("failed to read a quoted string");
	}
	// read an optional keyword
	bool getKeywordOptional(string& str)
	{
		bool lgSuccess = getTokenOptional(str);
		caps(str);
		return lgSuccess;
	}
	// read a keyword
	void getKeyword(string& str)
	{
		if( !getKeywordOptional(str) )
			errorAbort("failed to read a keyword");
	}
	// read line label plus wavelength
	void getLineID(LineID& line);
	// skip to a specified position on the line
	void skipTo(size_t p)
	{
		if( p < p_pos() )
			errorAbort("skipping to requested position failed");
		p_pos(p);
	}
	// skip to the first position just after the next instance of a specified string
	void skipAfter(const string& s)
	{
		auto cp = p_pos();
		auto p = p_line.substr(cp).find(s);
		if( p != string::npos )
			skipTo(cp + p + s.length());
		else
		{
			ostringstream oss;
			oss << "skipAfter could not find string =" << s << "=";
			errorAbort(oss.str());
		}
	}
	// returns true if the end of the line is reached
	bool lgEOL()
	{
		// skipping whitespace allows the trailing junk to be shown
		// correctly by showLocation() if the test for EOL fails...
		p_skipWS();
		return p_ls.eof();
	}
	// check if EOL is reached and abort if this is not the case
	void checkEOL()
	{
		if( !lgEOL() )
			errorAbort("found trailing junk at the end of this line");
	}
	// returns true if the end of the file is reached
	bool lgEOF() const { return p_lgEOF; }
	// returns true if the current line is an inline end-of-data marker
	bool lgEODMarker() const;
	// check if EOD is reached and abort if this is not the case
	void checkEOD()
	{
		if( getline() && !lgEODMarker() )
			errorAbort("found surplus input at the end of this file");
	}
	// get the current position in the stream buffer
	long getpos() { return p_io.tellg(); }
	// jump to specific position in the stream buffer (previously returned by getpos())
	void setpos(long pos) { p_io.seekg(pos, p_io.beg); }
	// CODE SMELL -- rewind file to the beginning
	void rewind()
	{
		p_io.clear();
		p_io.seekg(0);
		p_line.clear();
		p_nr = 0;
		p_ls.str("");
		p_ls.clear();
		p_lgEOF = false;
	}
	// abort with specific error message
	NORETURN void errorAbort(const string& msg, FILE *io = ioQQQ);
	// non-fatal warning message
	void warning(const string& msg, FILE *io = ioQQQ);
};

/** Links text string to an action on a specified argument */
template <typename V>
class KeyAction {
	const char * const m_keyword;
	V m_action;
public:
	KeyAction(const char *keyword, const V &action) :
		m_keyword(keyword), m_action(action) {}
		
	const char *key(void) const
	{
		return m_keyword;
	}
	void operator()(realnum *v) const
	{
		m_action(v);
	}
};

/** Helper template to make it easier to generate KeyActions */
template <typename V>
inline KeyAction<V> MakeKeyAction(const char *keyword, const V &action)
{
	return KeyAction<V>(keyword, action);
}

/** Generator for functors which convert the unit of their argument */
class UnitConverter
{
	const realnum m_unit;
public:
	UnitConverter ( double unit ) : m_unit((realnum)unit) {}
		
	void operator()( realnum *t ) const
	{
		*t *= m_unit;
	}
};

/** Interate through a list of KeyActions: apply the first which
	 matches and then quit */
template <typename T, typename V>
bool parserProcess(Parser &p, T *list, unsigned long nlist, V *value)
{
	bool lgFound = false;
	for (unsigned long option=0; option < nlist; ++option)
	{
		if( p.nWord( list[option].key() ) )
		{
			list[option]( value );
			lgFound = true;
			break;
		}
	}
	return lgFound;
}

/**ParseCosmicRays parse the cosmic rays command 
\param *chCard
*/
void ParseCosmicRays( Parser &p );

/**ParseCosmology parse the cosmology command 
\param *chCard
*/
void ParseCosmology( Parser &p );

/**ParseAbundances parse and read in composition as set by abundances command 
\param *chCard
\param lgDSet
*/
void ParseAbundancesNonSolar(Parser &p);

void ParseAbundances(Parser &p);

/**ParseDont parse the dont command */
void ParseDont(Parser &p);

/**ParseSave parse the save command 
\param *chCard
*/
void ParseSave(Parser &p);

void parse_save_line(Parser &p, 
		     /* true, return rel intensity, false, log of luminosity or intensity I */
		     bool lgLog3,
			  ostringstream& chHeader,
							long int ipPun
);

void parse_save_average(Parser &p,
			/* the file we will write to */
			long int ipPun, 
			ostringstream& chHeader);

void parse_save_colden(Parser &p,
		       /* the header for the file, a list of identifications */
		       ostringstream& chHeader);

void Parse_Save_Line_RT(Parser &p);

/**ParseAge - parse the age command */
void ParseAge(Parser &p);

/**ParseAgn parse parameters for the AGN continuum shape command 
\param *chCard
*/
void ParseAgn(Parser &p);

/** parse the blackbody command 
\param *chCard input command line, already changed to caps
\param *nqh counter for which continuum source this is
\param *ar1 optional area that might be set here
*/
void ParseBlackbody(Parser &p);					

/**ParseCompile compile werner or kurucz model atmospheres into cloudy format, by K Volk 
\param *chCard
*/
void ParseCompile(Parser &p );

/**ParseConstant parse the constant ... command */
void ParseConstant(Parser &p);

/**ParseDLaw parse parameters on the dlaw command so set some density vs depth law
\param *chCard
*/
void ParseDLaw(Parser &p );

/**ParseTLaw parse parameters on the tlaw command to set some temperature vs depth 
\param *chCard
*/
void ParseTLaw(Parser &p);

/**ParseGrain parse parameters on Peter's version of the grains command 
\param *chCard
\param *lgDset
*/
void ParseGrain(Parser &p);

/**ParseFluc parse the fluctuations command */
void ParseFluc(Parser &p);

/**ParseHDEN parse the HDEN command */
void ParseHDEN(Parser &p);

/**ParseDatabaseISO parse the atom XX-like command, to set options for iso sequences 
\param ipISO
\param *chCard
*/
void ParseDatabaseISO(long ipISO, Parser &p);

/**ParseDatabaseH2 parse information from the rotor command line 
\param *chCard
*/
void ParseDatabaseH2(Parser &p );

/**ParseGrid parse the grid command 
\param *chCard
*/
void ParseGrid(Parser &p);

/**ParseInit parse the init command */
void ParseInit(Parser &p);

/**ParseInitFile helper routine for ParseInit */
void ParseInitFile(const string& chName);

/**ParseInterp parse parameters on interpolate command 
\param *chCard
\param *lgEOF
*/
void ParseInterp(Parser &p);

/**ParseIonParI parse the ionization parameter command (IONI variant)
\param *nqh
\param *chCard
\param *chType
*/
void ParseIonParI(Parser &p);

/**ParseIonParX parse the ionization parameter command (XI variant)
\param *nqh
\param *chCard
\param *chType
*/

void ParseIonParX(Parser &p);
/**ParseIonPar parse the ionization parameter command 
\param *nqh
\param *chCard
\param *chType
*/
void ParseIonPar(Parser &p,
					  char chType);

/**ParseNorm parse parameters on the normalize command 
\param *chCard
*/
void ParseNorm(Parser &p);

/**ParseOptimize parse the optimize command 
\param *chCard
*/
void ParseOptimize(Parser &p);

/**ParsePrint parse the print command  
\param *chCard
*/
void ParsePrint(Parser &p );

/**ParseRadius parse the radius command */
void ParseRadius(Parser &p);

/**ParseSet parse the set command */
void ParseSet(Parser &p);

/**ParseTable parse the table read command 
\param *nqh
\param *chCard
\param *ar1
*/
void ParseTable(Parser &p);

/**ParseTrace parse the trace command */
void ParseTrace(Parser &p);

/*ParseExtinguish parse the extinguish command */
void ParseExtinguish( Parser &p );

/*ParseIlluminate parse the illuminate command */
void ParseIlluminate( Parser &p );

/*ParseCaseB - parse the Case B command */
void ParseCaseB(Parser &p );

/**ParseTest parse the test command */
void ParseTest(Parser &p);

/**ParseAbsMag parse the absolute magnitude command */
void ParseAbsMag(Parser &p);

/**ParseBackgrd parse the background continuum command */
void ParseBackgrd(Parser &p);

/**ParseCoronal parse the cronal equilibrum command */
void ParseCoronal(Parser &p);

/**ParseElement parse options on element command */
void ParseElement(Parser &p);

/**ParseCMB parse parameters from fireball command 
\param z
\param *nqh
\param *ar1
*/
void ParseCMB(double z, 
  long int *nqh);

/**ParseF_nu parse intensity command parameters 
\param *chCard
\param *nqh
\param *ar1
\param *chType
\param lgNU2
*/
void ParseF_nu(
  Parser &p, 
  const char *chType, 
  bool lgNU2);

/**ParseGlobule parse parameters off the globule command 
\param *chCard
*/
void ParseGlobule(Parser &p);

/**ParseRangeOption parse the range option on the luminosity command */
void ParseRangeOption(Parser &p);

/**ParseMap parse map command to produce map of heating and cooling */
void ParseMap(Parser &p);

/**ParseMetal parse parameters on metal command */
void ParseMetal(Parser &p);

void ParseLineList(Parser &p, vector<LineID>& lines);

void ParsePrtLineSum(Parser &p);

/**ParsePowerlawContinuum parse the power law continuum command */
void ParsePowerlawContinuum(Parser &p);

/**ParseRatio parse the ratio command */
void ParseRatio(Parser &p);

/**ParseSphere parse the sphere command */
void ParseSphere(Parser &p);

/**ParseStop parse the stop command */
void ParseStop(Parser &p);

/**ParseCrashDo any of several tests to check that the code can crash 
\param *chCard
*/
void ParseCrashDo(Parser &p);
 
class Option {
public:
	enum Opts { BOOL, LONG, REAL, OPTION, STRING, INVALID } opttype;
	enum Quoted { NOTQUOTED, QUOTED};
	union {
		bool   l;
		long   i;
		double r;
	};
	string s;
	Option(bool val)
	{
		opttype=BOOL;
		l = val;
	}
	Option(long val)
	{
		opttype=LONG;
		i = val;
	}
	Option(double val)
	{
		opttype=REAL;
		r = val;
	}
	Option(const string& val, enum Quoted q)
	{
		if (q == QUOTED)
			opttype=STRING;
		else
			opttype=OPTION;
		s = val;
	}
};

class Properties
{
	bool m_lgDone;
public:
	vector< pair<string, shared_ptr<Option> > > p;
   Properties() : m_lgDone(false) {}
	void setDone() 
	{
		m_lgDone = true;
	}
	bool isDone() const
	{
		return m_lgDone;
	}
};

#endif // _PARSER_H_
