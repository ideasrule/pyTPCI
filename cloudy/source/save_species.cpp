/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveSpecies generate output for the save species command */
#include "cddefines.h"
#include "taulines.h"
#include "radius.h"
#include "save.h"
#include "mole.h"
#include "mole_priv.h"
#include "generic_state.h"
#include "iterations.h"
#include "dense.h"
#include "prt.h"



STATIC void SaveSpeciesLines( FILE *ioPUN, const vector<genericState> &speciesList );

/* save results for one particular species */
STATIC void SaveSpeciesOne( 
	const vector<genericState>& SpeciesList,
	double (*job)(const genericState&),
	bool lgZonal,
	const char* chFmt,
	FILE *ioPUN );

STATIC void SaveSpeciesHeader( 
	const vector<genericState>& SpeciesList,
	const char *chJob,
	bool lgZonal,
	FILE *ioPUN,
	size_t maxLevels );

/*SaveSpecies generate output for the save species command */
void SaveSpecies(
		FILE* ioPUN,
		long int ipPun)
{
	DEBUG_ENTRY( "SaveSpecies()" );

	/* first branch; save results for all species if list is empty */
	vector<genericState> speciesList;
	if( save.chSaveSpecies[ipPun].size() == 0 )
	{
		// loop over species
		for( size_t i=0; i<mole_global.list.size(); ++i )
		{
			speciesList.push_back( genericState(&(mole.species[i])) );
		}
	}
	else
	{
		for (vector<string>::iterator name=save.chSaveSpecies[ipPun].begin();
			  name != save.chSaveSpecies[ipPun].end(); ++name)
		{
			vector<genericState> v = matchGeneric( *name, false );
			if( v.size() == 0 )
			{
				fprintf( ioQQQ,"Could not match species '%s', "
							"use SAVE SPECIES LABELS ALL to get a list of all species."
							"\nSorry.\n", name->c_str() );
				cdEXIT(EXIT_FAILURE);
			}
			speciesList.insert(speciesList.end(),v.begin(),v.end());
		}
	}

	if( strcmp( save.chSaveArgs[ipPun], "LABE" )==0 )
	{
		if( save.lgSaveHeader(ipPun) )
		{
			/* save list of species labels */
			fprintf( ioPUN, "#Species label\tDatabase\n" );
			for( size_t i=0; i<speciesList.size(); ++i )
			{
				fprintf( ioPUN, "%s\t%s\n",
					speciesList[i].label().c_str(),
					speciesList[i].database().c_str() );
			}
			save.SaveHeaderDone(ipPun);
		}
		return;
	}

	if( strcmp( save.chSaveArgs[ipPun], "DATA" ) == 0 )
	{
		SaveSpeciesLines( ioPUN, speciesList );
		return;
	}

	/* remaining options are column densities, densities, and energies */
	// max number of levels, for header print
	size_t mostLevels = 0;
	for (vector<genericState>::iterator sp=speciesList.begin();
		  sp != speciesList.end(); ++sp)
	{
		const molezone *saveSpecies = sp->sp;
		if( saveSpecies != NULL && saveSpecies != null_molezone &&
			 saveSpecies->levels != NULL )		
			mostLevels = MAX2(mostLevels, saveSpecies->levels->size() );
	}
	
	ASSERT( mostLevels < 10000 );

	bool lgZonal = true;
	const char* chJob=NULL, *chFmt=NULL;
	double (*job)(const genericState&)=NULL;
	if( strcmp( save.chSaveArgs[ipPun], "ENER" )==0 )
	{
		lgZonal = false;
		chJob = "energies";
		job = energy;
		chFmt = "%.5e";
	}
	else if( strcmp( save.chSaveArgs[ipPun], "DEPA" )==0 )
	{
		lgZonal = true;
		chJob = "departure coefficients";
		job = depart;
		chFmt = "%.5e";
	}
	else if( strcmp( save.chSaveArgs[ipPun], "DENS" )==0 )
	{
		lgZonal = true;
		chJob = "densities";
		chFmt = "%.5e";
		job = density;
	}
	else if( strcmp( save.chSaveArgs[ipPun], "LEVL" )==0 )
	{
		lgZonal = true;
		chJob = "levels";
		chFmt = "%4.0f";
		job = levels;
	}
	else if( strcmp( save.chSaveArgs[ipPun], "COLU" )==0 )
	{
		lgZonal = false;
		chJob = "column density";
		chFmt = "%.5e";
		job = column;
	}
	else
	{
		fprintf(ioQQQ,"PROBLEM, save species job type \"%s\" not known\n",save.chSaveArgs[ipPun]);
		TotalInsanity();
	}
	
	if( save.lgSaveHeader(ipPun) )
	{
		SaveSpeciesHeader( speciesList, chJob, lgZonal, ioPUN, mostLevels );
		save.SaveHeaderDone(ipPun);
	}
	SaveSpeciesOne( speciesList, job, lgZonal, chFmt, ioPUN );

	return;
}

/* print 0.000e+00 as simply 0 */
STATIC void PrintShortZero( FILE *ioPUN , const char* chFmt, double arg )
{
	DEBUG_ENTRY( "PrintShortZero()" );
	if( arg==0. )
		fprintf(ioPUN,"0");
	else
		fprintf(ioPUN,chFmt, arg);

}

// Switch from initial format to single row per zone
const bool lgRowPerZone = true;
/* save results for one particular species */
STATIC void SaveSpeciesHeader( 
	const vector<genericState>& speciesList,
	const char *chJob,
	bool lgZonal,
	FILE *ioPUN,
	size_t maxLevels )
{
	DEBUG_ENTRY( "SaveSpeciesHeader()" );

	if (!lgRowPerZone)
	{
		fprintf( ioPUN, "#%sspecies %s", lgZonal ? "depth\t":"", chJob );
		for( size_t i = 0; i < maxLevels; ++i )
		{
			fprintf( ioPUN, "\t%lu", (unsigned long)i );
		}
		fprintf( ioPUN, "\n");
	}
	else
	{
		int col=0;
		if (lgZonal)
		{
			fprintf( ioPUN, "#depth %s", chJob );
			++col;
		}
		else
		{
			fprintf( ioPUN, "#%s ", chJob );
		}
		for (vector<genericState>::const_iterator gs=speciesList.begin();
			  gs != speciesList.end(); ++gs)
		{
			if (col == 0)
				fprintf( ioPUN, "%s",gs->label().c_str() );
			else
				fprintf( ioPUN, "\t%s",gs->label().c_str() );
			++col;
		}
		fprintf( ioPUN,"\n" );
	}
}

STATIC void SaveSpeciesOne( 
	const vector<genericState>& speciesList,
	double (*job)(const genericState&),
	bool lgZonal,
	const char* chFmt,
	FILE *ioPUN)
{
	DEBUG_ENTRY( "SaveSpeciesOne()" );

	int col = 0;
	if (lgRowPerZone && lgZonal)
	{
		fprintf( ioPUN, "%.5e", radius.depth_mid_zone );
		++col;
	}
	for (vector<genericState>::const_iterator gs=speciesList.begin();
		  gs != speciesList.end(); ++gs)
	{
		if (!lgRowPerZone)
		{
			if (lgZonal)
			{
				fprintf( ioPUN, "%.5e", radius.depth_mid_zone );
				++col;
			}
			fprintf( ioPUN, "\t%s", gs->label().c_str() );
			++col;
		}
		if (col > 0)
			fprintf(ioPUN,"\t");
		PrintShortZero( ioPUN, chFmt, job(*gs) );
		++col;
		if (!lgRowPerZone)
		{
			fprintf(ioPUN,"\n");
			col = 0;
		}
	}
	
	if (lgRowPerZone)
		fprintf(ioPUN,"\n");
	return;
}

/**SaveAllSpeciesLabelsLevels -- generate output of all species labels & levels */
STATIC void SaveAllSpeciesLabelsLevels( FILE *ioPUN, const vector<genericState> &speciesList )
{
	//Print out number of levels used
	fprintf( ioPUN, "# The number of levels used in each species.\n" );
	fprintf( ioPUN, "# Species\tSpectrum\tUsed\tMax.\tDatabase\n" );
	for( size_t ipSpecies=0; ipSpecies < speciesList.size(); ++ipSpecies )
	{
		const molezone *this_mole = speciesList[ ipSpecies ].sp;
		if( this_mole == NULL || this_mole == null_molezone )
			continue; 

		species *this_species = (*this_mole).dbase;
		if( this_species == NULL )
			continue;

		fprintf( ioPUN, "%-8s", speciesList[ipSpecies].label().c_str() );

		string spectralLabel;
		chemical_to_spectral( speciesList[ ipSpecies ].label(), spectralLabel );
		fprintf( ioPUN, "\t%s", spectralLabel.c_str() );

		fprintf( ioPUN, "\t%li", (*this_species).numLevels_local );
		fprintf( ioPUN, "\t%li", (*this_species).numLevels_max );
		fprintf( ioPUN, "\t%s", speciesList[ipSpecies].database().c_str() );
		fprintf( ioPUN, "\n" );
	}
}

STATIC void SaveSpeciesLines( FILE *ioPUN, const vector<genericState> &speciesList )
{
	static vector<string> saved_species;

	vector<genericState> speciesList_local;

	for( auto spcit = speciesList.begin(); spcit != speciesList.end(); ++spcit )
	{
		bool found = false;

		for( auto sit = saved_species.begin(); sit != saved_species.end();
			++sit )
		{
			if( *sit == spcit->label() )
			{
				found = true;
				break;
			}
		}

		if( not found )
		{
			speciesList_local.insert( speciesList_local.end(), *spcit );
			saved_species.push_back( spcit->label() );
		}
	}

	if( speciesList_local.size() == 0 )
		return;

	SaveAllSpeciesLabelsLevels( ioPUN, speciesList_local );

	fprintf( ioPUN,"\n\n");

	//Species  ipLo ipHi gLo gHi wavelen EinA CS Rates
	fprintf( ioPUN,"Spectrum");

	if( save.lgSaveDataWn )
	{
		fprintf( ioPUN,"\tWavenumbers");
	}
	else
	{
		fprintf( ioPUN,"\tWavelength");
	}

	fprintf( ioPUN,"\tLo");
	fprintf( ioPUN,"\tHi");
	fprintf( ioPUN,"\tgLo");
	fprintf( ioPUN,"\tgHi");

	if( save.lgSaveDataGf )
	{
		fprintf( ioPUN,"\t   gf   ");
	}
	else
	{
		fprintf( ioPUN,"\tEinstein A");
	}
	fprintf( ioPUN,"\tColl_Str");
	fprintf( ioPUN,"\tgbar");

	if( save.lgSaveDataRates )
	{
		fprintf( ioPUN,"\tRate electron");
		fprintf( ioPUN,"\tRate proton");
		fprintf( ioPUN,"\tRate He+   ");
		fprintf( ioPUN,"\tRate Alpha");
		fprintf( ioPUN,"\tRate Atom H");
		fprintf( ioPUN,"\tRate Atom He");
		fprintf( ioPUN,"\tRate Ortho");
		fprintf( ioPUN,"\tRate Para");
		fprintf( ioPUN,"\tRate H2");
	}

	fprintf( ioPUN,"\n");


	for( size_t ipSpecies=0; ipSpecies < speciesList_local.size(); ++ipSpecies )
	{
		const molezone *this_mole = speciesList_local[ ipSpecies ].sp;
		if( this_mole == NULL || this_mole == null_molezone )
			continue;
		if( (*this_mole).lines == NULL )
			continue;

		species *this_species = (*this_mole).dbase;
		if( this_species == NULL )
			continue;

		for (TransitionList::iterator tr = (*this_mole).lines->begin();
			tr != (*this_mole).lines->end(); ++tr )
		{
			long ipLo = tr->ipLo() +1;
			long ipHi = tr->ipHi() +1;
			int nelem = tr->Hi()->nelem();

			if( nelem != -1 && !dense.lgElmtOn[nelem-1] )
			{
				continue;
			}

			if( ipHi >= (*this_species).numLevels_local )
			{
				continue;
			}

			string spectralLabel;
			chemical_to_spectral( speciesList_local[ ipSpecies ].label(), spectralLabel );
			fprintf( ioPUN,"%-8s", spectralLabel.c_str() );

			if( save.lgSaveDataWn )
			{
				fprintf( ioPUN,"\t%.3e", tr->EnergyWN() );
			}
			else
			{
				fprintf( ioPUN, "\t" );
				prt_wl( ioPUN, realnum( tr->WLAng() ) );
			}

			fprintf( ioPUN,"\t%li", ipLo);
			fprintf( ioPUN,"\t%li", ipHi);
			fprintf( ioPUN,"\t%i", int( tr->Lo()->g() ) );
			fprintf( ioPUN,"\t%i", int( tr->Hi()->g() ) );

			if( save.lgSaveDataGf )
			{
				fprintf( ioPUN,"\t%.3e", tr->Emis().gf() );
			}
			else
			{
				fprintf( ioPUN,"\t%.3e", tr->Emis().Aul() );
			}

			fprintf( ioPUN,"\t%.3e", tr->Coll().col_str());
			fprintf( ioPUN,"\t%i", tr->Coll().is_gbar() );

			if( save.lgSaveDataRates )
			{
				for( long intCollNo=0; intCollNo<ipNCOLLIDER; intCollNo++)
				{
					fprintf( ioPUN,"\t%.3e",tr->Coll().rate_coef_ul()[intCollNo]);
				}
			}

			fprintf( ioPUN,"\n");
		}
	}
}

void SaveSpeciesOptDep( const long int ipPun, const string &speciesLabel )
{
	DEBUG_ENTRY( "SaveSpeciesOptDep()" );

	if( ! iterations.lgLastIt )
		return;

	if( save.lgSaveHeader(ipPun) )
	{
		fprintf( save.params[ipPun].ipPnunit,
			"#hi\tlo\tWavelength(A)\ttau\n" );
		save.SaveHeaderDone( ipPun );
	}

	genericState species = getSpeciesGeneric( speciesLabel );

	if( species.sp->lines == NULL )
	{
		fprintf( ioQQQ,
			"WARNING: Species '%s' does not have any data for 'save species optical depth'.\n",
			species.label().c_str() );
		return;
	}

	for( TransitionList::iterator tr = (*species.sp->lines).begin();
		tr != (*species.sp->lines).end(); ++tr)
	{
		if( (*tr).ipCont() <= 0 )
			continue;

		fprintf( save.params[ipPun].ipPnunit,
			"%i\t%i\t%.5e\t%.5e\n",
			(*tr).ipHi()+1,
			(*tr).ipLo()+1,
			(*tr).WLAng(),
			(*tr).Emis().TauIn() * SQRTPI );
	}
}

class Field
{
public:
	const char* label;
	const char* fmt;
private:
	double m_value;
public:
	Field(const char* label, const char* fmt, double value) :
		label(label), fmt(fmt), m_value(value)
		{}
	double value() const
		{
			return m_value;
		}
};

STATIC void doHeader(FILE *punit, const Field& f)
{
	DEBUG_ENTRY( "doHeader()" );
	fprintf(punit,"%s",f.label);
}
STATIC void doData(FILE *punit, const Field& f)
{
	DEBUG_ENTRY( "doData()" );
	fprintf(punit,f.fmt,f.value());
}

STATIC inline bool isCatalystReactant(const mole_reaction& rate, int i)
{
	return rate.rvector[i]!=NULL;
}
STATIC inline bool isCatalystProduct(const mole_reaction& rate, int i)
{
	return rate.pvector[i]!=NULL;
}
STATIC inline bool isDestroyed(const mole_reaction& rate, int i)
{
	return !isCatalystReactant(rate,i) && rate.rvector_excit[i]==NULL;
}
STATIC inline bool isCreated(const mole_reaction& rate, int i)
{
	return !isCatalystProduct(rate,i) && rate.pvector_excit[i]==NULL;
}

/* Save all rates involving specified species */
void mole_save(FILE *punit, const char speciesname[], const char args[], bool lgHeader, bool lgData, bool lgCoef, double depth)
{
	DEBUG_ENTRY( "mole_save()" );

	molecule *sp = findspecies(speciesname);
	
	if( sp == null_mole )
	{
		fprintf( punit, " Species %s could not be found. Note that labels are case-sensitive in this context.\n", speciesname );
		return;
	}
		
	void (*doTask)(FILE *punit, const Field& f);

	if (lgHeader)
		doTask=doHeader;
	else
		doTask=doData;
	
	doTask(punit,Field("#Depth","%.5e",depth));

	for (mole_reaction_i p = mole_priv::reactab.begin(); p != mole_priv::reactab.end(); ++p)
	{
		const mole_reaction &rate = *p->second;
		bool ipthis = false;

		for (int i=0;i<rate.nreactants && !ipthis;i++)
		{
			if( rate.reactants[i] == sp )
			{
				if( ( strcmp( args, "DEST" )==0 && isDestroyed(rate,i) ) ||
					 ( strcmp( args, "DFLT" )==0 && isDestroyed(rate,i) ) ||
					 ( strcmp( args, "CATA" )==0 && isCatalystReactant(rate,i) ) ||
					strcmp( args, "ALL " )==0 )
					ipthis = true;
			}
		}
		
		for(int i=0;i<rate.nproducts && !ipthis;i++)
		{
			if( rate.products[i] == sp )
			{
				if( ( strcmp( args, "CREA" )==0 && isCreated(rate,i) ) ||
					( strcmp( args, "DFLT" )==0 && isCreated(rate,i) ) ||
					( strcmp( args, "CATA" )==0 && isCatalystProduct(rate,i) ) ||
					strcmp( args, "ALL " )==0 )
					ipthis = true;
			}
		}

		if(ipthis) 
		{
			double ratevi=0.0;
			if (lgData)
			{
				ratevi = mole.reaction_rks[ rate.index ];
				if( !lgCoef )
				{
					for(int i=0;i<rate.nreactants;i++)
					{
						ratevi *= mole.species[ rate.reactants[i]->index ].den;
					}
				}
			}
			fprintf(punit,"\t");

			doTask(punit,Field(rate.label.c_str(),"%.3e",ratevi));
		}
	}
	fprintf(punit,"\n");
}	

// Helper class for sorting rates
class RateCmp
{
	const vector<double>& m_stack;
public:
	RateCmp(const vector<double>& stack) : m_stack(stack) {}
	bool operator()(size_t a, size_t b)
		{
			return m_stack[a] > m_stack[b];
		}
};

void mole_dominant_rates( const vector<const molecule*>& debug_list, 
								  FILE *ioOut,
								  bool lgPrintReagents, size_t NPRINT, double fprint )
{
	vector<double> snkx, srcx;
	vector<mole_reaction *> ratesnk, ratesrc;

	double src_total = 0.0, snk_total = 0.0;

	for(mole_reaction_i p=mole_priv::reactab.begin();
			p != mole_priv::reactab.end(); ++p)
	{
		mole_reaction *rate = &(*p->second);
		double rk = mole.reaction_rks[ rate->index ];

		double rate_tot = rk;
		for(long i=0;i<rate->nreactants;i++)
		{
			rate_tot *= mole.species[ rate->reactants[i]->index ].den;
		}
		
		int nrate=0;
		for(size_t s=0; s<debug_list.size(); ++s)
		{
			for(long i=0;i<rate->nproducts;++i)
			{
				molecule *sp = rate->products[i];
				if (sp == debug_list[s] && rate->pvector[i] == NULL)
				{
					++nrate;
				}
			}
			for(long i=0;i<rate->nreactants;++i)
			{
				molecule *sp = rate->reactants[i];
				if (sp == debug_list[s] && rate->rvector[i] == NULL)
				{
					--nrate;
				}
			}
		}
		if (nrate > 0)
		{
			srcx.push_back(nrate*rate_tot);
			ratesrc.push_back(rate);
			src_total += nrate*rate_tot;
		}
		else if (nrate < 0)
		{
			snkx.push_back(-nrate*rate_tot);
			ratesnk.push_back(rate);
			snk_total -= nrate*rate_tot;
		}
	}

	if (!ratesrc.empty())
	{
		vector<size_t> isrc(ratesrc.size());
		for (size_t i=0; i<ratesrc.size(); ++i)
			isrc[i] = i;
		sort(isrc.begin(),isrc.end(),RateCmp(srcx));

		fprintf( ioOut, "Src %13.7g: ",src_total);
		for (size_t i=0; i<ratesrc.size(); ++i)
		{
			if (i == NPRINT || srcx[isrc[i]] < fprint*src_total)
			{
				break;
			}
			fprintf( ioOut, "%20.20s %13.7g",
						ratesrc[isrc[i]]->label.c_str(),srcx[isrc[i]]);
			if (lgPrintReagents)
			{
				fprintf( ioOut, " [");
				for (long j=0;j<ratesrc[isrc[i]]->nreactants;j++)
				{
					if (j)
					{
						fprintf( ioOut, "," );
					}
					fprintf( ioOut, "%-6.6s %13.7g",
								ratesrc[isrc[i]]->reactants[j]->label.c_str(),
								mole.species[ ratesrc[isrc[i]]->reactants[j]->index ].den);
				}
				fprintf( ioOut, "]" );
			}
			else
			{
				fprintf( ioOut, ";" );
			}
		}
	}
	if (!ratesnk.empty())
	{		
		vector<size_t> isnk(ratesnk.size());
		for (size_t i=0; i<ratesnk.size(); ++i)
			isnk[i] = i;
		sort(isnk.begin(),isnk.end(),RateCmp(snkx));

		fprintf( ioOut, " Snk %13.7g: ", snk_total);
		for (size_t i=0; i<ratesnk.size(); ++i)
		{
			if (i == NPRINT || snkx[isnk[i]] < fprint*snk_total)
			{
				break;
			}
			fprintf( ioOut, "%20.20s %13.7g",
						ratesnk[isnk[i]]->label.c_str(), snkx[isnk[i]] );
			if (lgPrintReagents)
			{
				fprintf( ioOut, " [" );
				for (long j=0;j<ratesnk[isnk[i]]->nreactants;j++)
				{
					if (j)
					{
						fprintf( ioOut, "," );
					}
					fprintf( ioOut, "%-6.6s %13.7g",
								ratesnk[isnk[i]]->reactants[j]->label.c_str(),
								mole.species[ ratesnk[isnk[i]]->reactants[j]->index ].den);
				}
				fprintf( ioOut, "]" );
			}
			else
			{
				fprintf( ioOut, ";" );
			}
		}
	}
	fprintf( ioOut, "\n" );

	return;
}

void mole_print_species_reactions( molecule *speciesToPrint )
{
	if( speciesToPrint==NULL )
	{	
		fprintf( ioQQQ,"\n NULL species found in mole_print_species_reactions.\n" );
		return;
	}

	long numReacts = 0;

	fprintf( ioQQQ,"\n Reactions involving species %s:\n", speciesToPrint->label.c_str() );
	for( mole_reaction_i p=mole_priv::reactab.begin(); p!=mole_priv::reactab.end(); ++p )
	{
		mole_reaction &rate = *p->second;
		for( long i=0; i<rate.nreactants; i++ )
		{
			molecule *sp = rate.reactants[i];
			if(rate.rvector[i]==NULL && sp==speciesToPrint )
			{
				double drate = mole.reaction_rks[ rate.index ];
				for (long j=0; j<rate.nreactants; ++j)
				{
					drate *= mole.species[ rate.reactants[j]->index ].den;
				}
				fprintf( ioQQQ,"%s rate = %g\n", rate.label.c_str(), drate );	
				numReacts++;
			}
		}
		for( long i=0; i<rate.nproducts; i++ )
		{
			molecule *sp = rate.products[i];
			if(rate.pvector[i]==NULL && sp==speciesToPrint )
			{
				double drate = mole.reaction_rks[ rate.index ];
				for (long j=0; j<rate.nreactants; ++j)
				{
					drate *= mole.species[ rate.reactants[j]->index ].den;
				}
				fprintf( ioQQQ,"%s rate = %g\n", rate.label.c_str(), drate );	
				numReacts++;
			}
		}
	}
	fprintf( ioQQQ," End of reactions involving species %s.  There were %li.\n", speciesToPrint->label.c_str(), numReacts );
	
	return;
}
