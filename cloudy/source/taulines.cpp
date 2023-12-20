/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "taulines.h"
#include "atmdat.h"

INSTANTIATE_MULTI_ARR( qList, lgBOUNDSCHECKVAL )

vector<TransitionList> AllTransitions;

bool lgStatesAdded = false;
bool lgLinesAdded = false;
multi_arr<qList,2> StatesElemNEW;
char **chSpecies;
vector<species> dBaseSpecies;
vector<qList> dBaseStates;
vector< multi_arr<int,2> > ipdBaseTrans;
vector<TransitionList> dBaseTrans;
multi_arr<CollRateCoeffArray,2> AtmolCollRateCoeff;
vector< multi_arr<CollSplinesArray,3> > AtmolCollSplines;
vector< StoutCollArray > StoutCollData;
long int nSpecies;
qList AnonStates("AnonStates",1);
multi_arr<int,3> ipExtraLymanLines;
vector<vector<TransitionList> > ExtraLymanLines;
TransitionList UTALines("UTALines", &AnonStates);
long int nLevel1;
TransitionList HFLines("HFLines", &AnonStates);
//vector<vector<multi_arr<int,2> > > ipTransitions;
vector<vector<TransitionList> > Transitions;
multi_arr<int,2> ipFe2LevN;
static qList Fe2LevNStates;
TransitionList Fe2LevN("Fe2LevN", &Fe2LevNStates);
multi_arr<int,3> ipSatelliteLines; /* [ipISO][nelem][level] */
vector<vector<TransitionList> > SatelliteLines; /* [ipISO][nelem][level] */
TransitionList TauLine2("TauLine2", &AnonStates);
vector<realnum> cs1_flag_lev2;

extern void checkTransitionListOfLists(vector<TransitionList>&list)
{ 
	for (vector<TransitionList>::iterator it=list.begin(); 
		  it != list.end(); ++it)
	{
		for (TransitionList::iterator tr = it->begin();
			  tr != it->end(); ++tr)
		{
			(*tr).check();
		}
		for (EmissionList::iterator em = it->Emis().begin(); 
			  em != it->Emis().end(); ++em)
		{
			(*em).check();
		}
	}
}

TransitionList::iterator findTrans_byQuantNumb( const string speciesLabel,
		const long n_hi, const long l_hi, const long S_hi,
		const long n_lo, const long l_lo, const long S_lo )
{
	TransitionList::iterator matchedTrans = AllTransitions.back().end();
	bool lgFound = false;

	for (vector<TransitionList>::iterator it = AllTransitions.begin(); 
		  it != AllTransitions.end(); ++it)
	{
		if( (*(it->begin())).chLabel().find( speciesLabel ) == string::npos )
			continue;

		for (TransitionList::iterator tr = it->begin();
			  tr != it->end(); ++tr)
		{
			if( (*tr).Lo()->n() == n_lo && (*tr).Lo()->l() == l_lo &&
			    (*tr).Lo()->S() == S_lo &&
			    (*tr).Hi()->n() == n_hi && (*tr).Hi()->l() == l_hi &&
			    (*tr).Hi()->S() == S_hi )
			{
				matchedTrans = tr;
				lgFound = true;
				break;
			}
		}

		if( lgFound )
			break;
	}

	return matchedTrans;
}

TransitionList::iterator findTrans_byWLAng( string speciesLabel, const double wl_Ang,
			double &wl_err )
{
	TransitionList::iterator matchedTrans = AllTransitions.back().end();
	double dwl = 1e30;

	if( wl_Ang < 0. )
		return matchedTrans;

	for (vector<TransitionList>::iterator it = AllTransitions.begin(); 
		  it != AllTransitions.end(); ++it)
	{
		if( (*(it->begin())).chLabel().find( speciesLabel ) == string::npos )
			continue;

		for (TransitionList::iterator tr = it->begin();
			  tr != it->end(); ++tr)
		{
			if( fabs( (*tr).WLAng() - wl_Ang ) < dwl )
			{
				wl_err = (*tr).WLAng() - wl_Ang;
				dwl = fabs( wl_err );
				matchedTrans = tr;
			}
		}
	}

	return matchedTrans;
}

TransitionProxy::iterator TauDummy;
EmissionProxy DummyEmis;

namespace 
{
	class Init
	{
		EmissionList DummyEmisList;
		TransitionListImpl TauDummyTrans;
	public:
		Init(qList*states) : 
			DummyEmisList(&TauDummyTrans, 1), TauDummyTrans("TauDummy",states, 1)
		{
			DummyEmis = DummyEmisList[0];
			TauDummy=TauDummyTrans.begin();
		};
	};
	qList TauDummyStates("TauDummyStates",1);
	Init TauDummyInit(&TauDummyStates);
}
