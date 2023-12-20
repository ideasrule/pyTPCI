/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "lines.h"
#include "prt.h"
#include "lines_service.h"
#include "version.h"
#include "service.h"

/* >>chng 21 april 7, H I wl from 4861.33A to 4861.34A, replacing average energy with np energies*/
/* >>chng 21 august 18, H I wl from 4861.34A to 4861.32A, revert to H energies computed by Peter van Hoof*/
const realnum Hbeta_WavLen = 4861.32f;

t_LineSave LineSave;
/* these are the definitions of the line save arrays in lines.h */

void t_LineSave::zero()
{
	DEBUG_ENTRY( "t_LineSave::zero()" );

	/* index within the line in the line stack 
	 * default is Hbeta total - the third line in the stack
	 * 0th is a zero for sanity, 1st is unit, 2nd is a comment */
	/* >>chng 02 apr 22 from 2 to 3 since added unit at 1 */
	/* >>chng 06 mar 11, from 3 to -1 will now set to "H  1" 4861 */
	ipNormWavL = -1;
	WavLNorm = Hbeta_WavLen;
	lgNormSet = false;
	sig_figs = sig_figs_max;

	/* the label for the normalization line */
	strcpy( chNormLab, "    " );

	/* The scale factor for the normalization line */
	ScaleNormLine = 1.;

}

void LinSv::prt(FILE* ioPUN) const
{
	DEBUG_ENTRY( "LinSv::prt()" );
	fprintf(ioPUN,"%s",label().c_str());
}

string LinSv::label() const
{
	DEBUG_ENTRY( "LinSv::label()" );
	string val = chALab();
	val.resize( NCHLAB-1, ' ' );
	val += " ";
	string buf;
	sprt_wl(buf, wavelength());
	val += buf;
	return val;
}

string LinSv::biglabel() const
{
	DEBUG_ENTRY( "LinSv::biglabel()" );
	string val = label();
	if (m_tr.associated())
	{
		val += " ";
		val += "[";
		val += (*m_tr.Hi()).chConfig();
		val += "->";
		val += (*m_tr.Lo()).chConfig();
		val += "]";
	}
	return val;
}

bool LinSv::isCat(const char *s) const
{
	DEBUG_ENTRY( "LinSv::isCat()" );

	const char* lbl = chALab();
	while( *s != '\0' )
		if( *lbl++ != *s++ )
			return false;
	return true;
}

void LinSv::chALabSet(const char *that)
{
	/* check that null is correct, string overruns have 
	 * been a problem in the past */
	ASSERT( (int) strlen( that )< NCHLAB );
	strncpy(m_chALab,that,NCHLAB-1);
	m_chALab[NCHLAB-1] = '\0';
	trimTrailingWhiteSpace( m_chALab );

	strncpy(m_chCLab,that,NCHLAB-1);
	m_chCLab[NCHLAB-1] = '\0';
	trimTrailingWhiteSpace( m_chCLab );
	caps(m_chCLab);
}

void LinSv::init(long index, char chSumTyp, const char *chComment, const char *label,
					  const TransitionProxy& tr)
{
	DEBUG_ENTRY( "LinSv::init()" );
	/* first call to stuff lines in array, confirm that label is one of
	 * the four correct ones */
	ASSERT( (chSumTyp == 'c') || (chSumTyp == 'h') || (chSumTyp == 'i') || (chSumTyp == 'r')  || (chSumTyp == 't') );
	ASSERT(!m_tr.associated());
	ASSERT(index >= 0);
	m_index = index;
	/* then save it into array */
	m_chSumTyp = chSumTyp;
	emslinZero();
	m_chComment = chComment;
	chALabSet( label );
	if (isCat("####"))
	{
		m_type = SEPARATOR;
	}
	else if (isCat("Unit"))
	{
		m_type = UNIT;
	}
	else if (isCat("UntD"))
	{
		m_type = UNITD;
	}
	else if (isCat("Inwd"))
	{
		m_type = INWARD;
	}
	else if (isCat("InwC"))
	{
		m_type = INWARDCONTINUUM;
	}
	else if (isCat("InwT"))
	{
		m_type = INWARDTOTAL;
	}
	else if (isCat("Coll"))
	{
		m_type = COLLISIONAL;
	}
	else if (isCat("Pump"))
	{
		m_type = PUMP;
	}
	else if (isCat("Heat"))
	{
		m_type = HEAT;
	}
	else if (isCat("Ca A"))
	{
		m_type = CASEA;
	}
	else if (isCat("Ca B"))
	{
		m_type = CASEB;
	}
	else if (isCat("nInu"))
	{
		m_type = NINU;
	}
	else if (isCat("nFnu"))
	{
		m_type = NFNU;
	}
	else if (isCat("Pho+"))
	{
		m_type = PHOPLUS;
	}
	else if (isCat("Pcon"))
	{
		m_type = PCON;
	}
	else if (isCat("Q(H)"))
	{
		m_type = QH;
	}
	else
	{
		m_type = DEFAULT;
	}

	SumLineZero();
	m_tr = tr;
}

void LinSv::addComponentID(long id)
{
	if ( LineSave.ipass == 0 )
	{
		if (id <= 0)
		{
			fprintf( ioQQQ, "ERROR: line blend component has invalid id: %ld.\n", id );
			cdEXIT( EXIT_FAILURE );
		}
		m_component.push_back(id);
		if (m_component.size() == 1)
			m_chComment += ": ";
		else
			m_chComment += "+";
		m_chComment += "\""+LineSave.lines[id].label()+"\"";
	}
}

void LinSv::addComponent(const LineID& line)
{
	if ( LineSave.ipass == 0 )
	{
		long id = LineSave.findline(line);
		if (id <= 0)
		{
			fprintf( ioQQQ, "ERROR: line blend component \"%s\" %.3f was not found.\n",
					 line.chLabel.c_str(), line.wave );
			cdEXIT( EXIT_FAILURE );
		}
		if (LineSave.lines[id].LineType() != 't')
		{
			fprintf( ioQQQ, "ERROR: line blend component \"%s\" %.3f does not have type 't'.\n",
					 line.chLabel.c_str(), line.wave );
			cdEXIT( EXIT_FAILURE );
		}
		auto tr = LineSave.lines[id].getTransition();
		if (!tr.associated())
		{
			fprintf( ioQQQ, "ERROR: line blend component \"%s\" %.3f is not associated with a transition.\n",
					 line.chLabel.c_str(), line.wave );
			cdEXIT( EXIT_FAILURE );
		}
		addComponentID(id);
	}
}
void LinSv::addComponent(const string& species,const double wavelength1)
{
	addComponent(LineID(species, wlAirVac(wavelength1)));
}

// Automatically generate blend for specified species, at wavelength +/- width
void LinSv::makeBlend(const char* chLabel,const double wavelength1, const double width)
{
	DEBUG_ENTRY("LinSv::makeBlend()");
	if ( LineSave.ipass == 0 )
	{
		realnum wlo = wlAirVac(wavelength1-width),
			whi = wlAirVac(wavelength1+width);
			
		/* check that chLabel[4] is null - supposed to be 4 char + end */
		if( strlen(chLabel) > NCHLAB-1 )
		{
			fprintf( ioQQQ, " makeBlend called with insane species \"%s\", must be %d or less characters long.\n",
						chLabel, NCHLAB-1 );
			cdEXIT(EXIT_FAILURE);
		}
		
		string chCARD(chLabel);
		
		/* make sure chLabel is all caps */
		caps(chCARD);/* convert to caps */
		
		/* this replaces tabs with spaces. */
		/* \todo	2 keep this in, do it elsewhere, just warn and bail? */
		for( size_t i=0; i < chCARD.length(); ++i )
		{
			if( chCARD[i] == '\t' )
			{
				chCARD[i] = ' ';
			}
		}

		for( long j=1; j < LineSave.nsum; j++ )
		{
			/* use pre-capitalized version of label to be like input chLineLabel */
			const char *chCaps = LineSave.lines[j].chCLab();

			if (LineSave.wavelength(j) >= wlo && 
				 LineSave.wavelength(j) <= whi &&
			    strcmp(chCaps,chCARD.c_str()) == 0)
			{
				addComponentID(j);
				fprintf( ioQQQ,"Adding \"%s\" to blend\n", 
							LineSave.lines[j].label().c_str() );
			}
		}
	}
}

void LinSv::setBlendWavl()
{
	if( LineSave.ipass == 0 && isBlend() )
	{
		realnum	sum_wn_num = 0.,
			sum_wn_den = 0.;
		for (size_t j=0; j<m_component.size(); ++j)
		{
			long id = m_component[j];
			TransitionProxy tr = LineSave.lines[id].getTransition();
			realnum weight = tr.Hi()->g() * tr.Emis().Aul();
			sum_wn_den += weight;
			sum_wn_num += weight * tr.EnergyWN();
		}
		double wl = 0.;
		if( sum_wn_den > 0. )
			sum_wn_num /= sum_wn_den;
		if( sum_wn_num > 0. )
			wl = wn2ang( sum_wn_num );

		LineSave.resetWavelength( m_index, wl );
	}
}

string LinSv::chComment() const
{
	return m_chComment;
}

static bool wavelength_compare (long a, long b)
{
	LinSv* a1 = &LineSave.lines[a];
	LinSv* b1 = &LineSave.lines[b];
	/* comparison is b-a so we get inverse wavelength order (increasing energy order) */
	if( b1->wavelength() < a1->wavelength() )
		return true;
	else 
		return false;
}

static bool wavelength_compare_realnum (size_t a, realnum wavelength)
{
	/* comparison is b-a so we get inverse wavelength order (increasing energy order) */
	if( wavelength < LineSave.wavelength(a) )
		return true;
	else 
		return false;
}

void t_LineSave::setSortWL()
{
	DEBUG_ENTRY( "t_LineSave::setSortWL()" );
	SortWL.resize((unsigned)nsum);
	for (long nlin=0; nlin < nsum; ++nlin)
	{
		SortWL[nlin] = nlin;
	}
	stable_sort(SortWL.begin(), SortWL.end(), wavelength_compare);
}

long t_LineSave::findline(const LineID& line)
{
	DEBUG_ENTRY( "t_LineSave::findline()" );

	ASSERT( nsum >= 0 );
	if (nsum == 0)
		return -1;

	bool lgDEBUG = false;
	
	if( line.chLabel.length() > NCHLAB-1 )
	{
		fprintf( ioQQQ, " findline called with insane chLabel (between quotes) \"%s\","
				 " must be no more than %d characters long.\n",
				 line.chLabel.c_str(), NCHLAB-1 );
		return -2;
	}
 
	// make static to avoid constant allocation and deallocation of memory for this var
	static string chCARD;
	chCARD = line.chLabel;

	/* make sure chLabel is all caps */
	caps(chCARD);/* convert to caps */
	trimTrailingWhiteSpace( chCARD );

	/* this replaces tabs with spaces. */
	/* \todo	2 keep this in, do it elsewhere, just warn and bail? */
	for( size_t i=0; i < chCARD.length(); ++i )
	{
		if( chCARD[i] == '\t' )
		{
			chCARD[i] = ' ';
		}
	}

	/* get the error associated with specified significant figures; add
	 * FLT_EPSILON*wavelength to broaden bounds enough to allow for
	 * cancellation error
	 */
	realnum errorwave = WavlenErrorGet( line.wave, LineSave.sig_figs ) + FLT_EPSILON*line.wave, 
		smallest_error=BIGFLOAT,
		smallest_error_w_correct_label=BIGFLOAT;

	// Possibly more 'user friendly' line search mode -- not active
	// NB falls through to previous implementation if ipass < 1, as
	// lines will not yet be sorted
	if (ipass == 1)
	{
		string wlbuf;
		// Need to search for lines within allowed band, and plausible confusions

		// Find position in list of lines
		vector<size_t>::iterator first =
			lower_bound(SortWL.begin(),SortWL.end(),
						line.wave+errorwave, wavelength_compare_realnum);

		// first is now first line below upper limit
		if (first == SortWL.end())
		{
			sprt_wl(wlbuf,line.wave);
			fprintf(ioQQQ,"Didn't find anything at %s\n",wlbuf.c_str());
			cdEXIT(EXIT_FAILURE);
		}

		// look for first line below lower limit -- may be the same as first if there are no matches
		vector<size_t>::iterator second;
		for(second=first; second != SortWL.end(); ++second)
		{
			if (wavelength(*second) < line.wave-errorwave)
				break;
		}

		vector<size_t>::iterator found = SortWL.end();
		int nmatch=0;
		realnum dbest = 0.;
		for (vector<size_t>::iterator pos = first; pos != second; ++pos)
		{
			bool lgMatch = ( strcmp(lines[*pos].chCLab(),chCARD.c_str()) == 0 );
			// do line disambiguation, if requested...
			if( lgMatch && line.indLo > 0 && line.indHi > 0 )
			{
				auto tr = lines[*pos].getTransition();
				if( tr.associated() )
					lgMatch = ( line.indLo == tr.ipLo()+1 && line.indHi == tr.ipHi()+1 );
				else
					lgMatch = false;
			}
			if( lgMatch && line.ELo >= 0_r )
			{
				auto tr = lines[*pos].getTransition();
				if( tr.associated() )
					lgMatch = ( fabs(tr.Lo()->energy().WN()-line.ELo) <= 1.e-6_r*line.ELo );
				else
					lgMatch = false;
			}
			if ( lgMatch )
			{
				++nmatch;
				realnum dwl = wavelength(*pos)-line.wave;
				if ( nmatch >= 2 )
				{
					if ( nmatch == 2 )
					{
						sprt_wl(wlbuf,line.wave);
						fprintf(ioQQQ,"WARNING: multiple matching lines found in search for \"%s\" %s\n",
								line.chLabel.c_str(),wlbuf.c_str());
						fprintf(ioQQQ,"WARNING: match 1 is \"%s\" (dwl=%gA)\n",
								lines[*found].biglabel().c_str(),wavelength(*found)-line.wave);
					}
					fprintf(ioQQQ,"WARNING: match %d is \"%s\" (dwl=%gA)\n",
							nmatch, lines[*pos].biglabel().c_str(),dwl);
				}
				if ( found == SortWL.end() )
				{
					found = pos;
					dbest = fabs(wavelength(*pos)-line.wave);
				}
				else if ( fabs(dwl) < dbest )
				{
					found = pos;
					dbest = fabs(dwl);
				}
			}
		}
		if ( found != SortWL.end() )
		{
			if (0)
				fprintf(ioQQQ,"Found %s\n", lines[*found].label().c_str());
			return *found;
		}

		// do this case first as output below would be confusing when line disambiguation is used
		if( ( line.indLo > 0 && line.indHi > 0 ) || line.ELo >= 0_r )
		{
			fprintf(ioQQQ, "PROBLEM: Line disambiguation was requested, but line was not found: \"%s\" ",
					line.chLabel.c_str());
			prt_wl(ioQQQ, line.wave);
			if( line.indLo > 0 && line.indHi > 0 )
				fprintf(ioQQQ, " index=%d, %d", line.indLo, line.indHi);
			if( line.ELo >= 0_r )
				fprintf(ioQQQ, " Elow=%g", line.ELo);
			fprintf(ioQQQ, "\nCheck the SAVE LINE LABELS output to find the correct match.\n");
			return -4;
		}

		sprt_wl(wlbuf,line.wave);
		fprintf(ioQQQ,"WARNING: no exact matching lines found for \"%s\" %s\n",line.chLabel.c_str(),wlbuf.c_str());
		for (vector<size_t>::iterator pos = first; pos != second; ++pos)
		{
			fprintf(ioQQQ,"WARNING: Line with incorrect label found close \"%s\"\n",
					  lines[*pos].label().c_str());
		}
		// Haven't found a match with correct label
		vector<size_t>::iterator best = SortWL.end();
		realnum besterror = 0.;
		for (;;)
		{
			realnum errordown = wavelength(*(first-1))-line.wave;
			realnum errorup = line.wave - (second == SortWL.end() ? 0.0 : wavelength(*second)) ;
			realnum error = 0.;
			vector<size_t>::iterator next;
			if ( errordown < errorup || second == SortWL.end())
			{
				error = errordown;
				--first;
				next = first;
			}
			else
			{
				error = errorup;
				next = second;
				++second;
			}
			if ( strcmp(lines[*next].chCLab(),chCARD.c_str()) == 0 )
			{
				if (best == SortWL.end())
				{
					best = next;
					besterror = error;
					fprintf(ioQQQ,"Taking best match as \"%s\"\n",lines[*next].label().c_str());
				}
				else
				{
					fprintf(ioQQQ,"Possible ambiguity with \"%s\"\n",lines[*next].label().c_str());
				}
			}
			// Assume this is clearly unambiguous
			if (best != SortWL.end() && error > 100.*besterror)
				break;
			// Assume this is clearly unmatched
			if (error > 0.01*line.wave)
				break;
		}
		if (best != SortWL.end())
			return *best;

		sprt_wl(wlbuf,line.wave);
		fprintf(ioQQQ,"PROBLEM: no matching line found in search for \"%s\" %s\n",line.chLabel.c_str(),wlbuf.c_str());
		cdEXIT(EXIT_FAILURE);
	}

	// Falls through to previous implementation if ipass < 1.  Can
	// be more restrictive with match handling, as this will be from
	// a request internal to the code infrastructure, rather than
	// user data

	long int j, index_of_closest=LONG_MIN,
		index_of_closest_w_correct_label=-1;
	
	for( j=1; j < nsum; j++ )
	{
		realnum current_error = (realnum)fabs(wavelength(j)-line.wave);
		/* use pre-capitalized version of label to be like input chLineLabel */
		const char *chCaps = lines[j].chCLab();

		// do line disambiguation, if requested...
		if( line.indLo > 0 && line.indHi > 0 )
		{
			auto tr = lines[j].getTransition();
			if( !tr.associated() || line.indLo != tr.ipLo()+1 || line.indHi != tr.ipHi()+1 )
				continue;
		}
		if( line.ELo >= 0_r )
		{
			auto tr = lines[j].getTransition();
			if( !tr.associated() || fabs(tr.Lo()->energy().WN()-line.ELo) > 1.e-6_r*line.ELo )
				continue;
		}

		if( current_error < smallest_error )
		{
			index_of_closest = j;
			smallest_error = current_error;
		}
			
		if( current_error < smallest_error_w_correct_label &&
		    (strcmp(chCaps,chCARD.c_str()) == 0) )
		{
			index_of_closest_w_correct_label = j;
			smallest_error_w_correct_label = current_error;
		}

		/* check wavelength and chLabel for a match */
		if( lgDEBUG && (current_error <= errorwave || 
				fp_equal( line.wave + errorwave, wavelength(j) ) ||
				fp_equal( line.wave - errorwave, wavelength(j) ))  
		    && strcmp(chCaps,chCARD.c_str()) == 0 )
		{
			/* match, so set emiss to emissivity in line */
			/* and announce success by returning line index within stack */
			printf("Matched %s %15.8g %ld %18.11g %s\n",
					 chCaps,line.wave,j,wavelength(j),
					 lines[j].biglabel().c_str());
		}
	}

	// Didn't find line, handle error
	if( index_of_closest_w_correct_label == -1 ||
		smallest_error_w_correct_label > errorwave )
	{
		// do this case first as output below would be confusing when line disambiguation is used
		if( ( line.indLo > 0 && line.indHi > 0 ) || line.ELo >= 0_r )
		{
			fprintf(ioQQQ, "PROBLEM: Line disambiguation was requested, but line was not found: \"%s\" ",
					line.chLabel.c_str());
			prt_wl(ioQQQ, line.wave);
			if( line.indLo > 0 && line.indHi > 0 )
				fprintf(ioQQQ, " index=%d, %d", line.indLo, line.indHi);
			if( line.ELo >= 0_r )
				fprintf(ioQQQ, " Elow=%g", line.ELo);
			fprintf(ioQQQ, "\nCheck the SAVE LINE LABELS output to find the correct match.\n");
			return -5;
		}

		/* >>chng 05 dec 21, report closest line if we did not find exact match, note that
		 * exact match returns above, where we will return negative number of lines in stack */
		fprintf( ioQQQ," PROBLEM findline did not find line " );
		prt_line_err( ioQQQ, chCARD.c_str(), line.wave );
		if( index_of_closest >= 0 )
		{
			fprintf( ioQQQ,"  The closest line (any label) was   \"%s\"\n", 
						lines[index_of_closest].label().c_str() );
			if( index_of_closest_w_correct_label >= 0 )
			{
				fprintf( ioQQQ,"  The closest with correct label was \"%s\"\n", 
							lines[index_of_closest_w_correct_label].label().c_str() );
				fprintf( ioQQQ,"  Error was %15.8g vs. tolerance %15.8g\n", 
							smallest_error_w_correct_label, errorwave );
			}
			else
				fprintf( ioQQQ,"\n  No line found with label \"%s\".\n", chCARD.c_str() );
			fprintf( ioQQQ,"\n" );
		}
		else
		{
			fprintf( ioQQQ,".\n PROBLEM No close line was found\n" );
			TotalInsanity();
		}
		return -3;
	}

	if (lgDEBUG)
		fprintf(ioQQQ,"Identified %ld\n",index_of_closest_w_correct_label);

	return index_of_closest_w_correct_label;
}

/*************************************************************************
 *
 * cdEmis obtain the local emissivity for a line with known index
 *
 ************************************************************************/
void cdEmis(
	const LinSv* line,
	/* the vol emissivity of this line in last computed zone */
	double *emiss ,
	// intrinsic or emergent
	bool lgEmergent )
{
	DEBUG_ENTRY( "cdEmis()" );

	if (line)
		*emiss = line->emslin(lgEmergent);
	else
		*emiss = 0.;
}
