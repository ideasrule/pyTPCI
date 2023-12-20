/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef MONITOR_RESULTS_H_
#define MONITOR_RESULTS_H_


 /**
  ParseMonitorResults - parse input stream
 */ 
class Parser;
void ParseMonitorResults(Parser &);

 /**
  must be called before rest, initializes assert variables
 */ 
void InitMonitorResults(void);


 /**
  lgCheckMonitors		 
  \param *ioMONITORS this is unit we will write output to
 */ 
bool lgCheckMonitors(
	FILE *ioMONITORS );

/** these flags are set in lgCheckMonitors */
extern bool lgMonitorsOK , lgBigBotch, lgPrtSciNot;

struct t_monitorresults {
	double SumErrorCaseMonitor;
	long int nSumErrorCaseMonitor;

};
extern t_monitorresults MonitorResults;

#endif /* MONITOR_RESULTS_H_ */
