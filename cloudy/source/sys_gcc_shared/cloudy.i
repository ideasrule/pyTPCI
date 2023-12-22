%module cloudy
%{
#include "cddefines.h"
#include "cddrive.h"
%}

FILE *fopen(char *, char *);
int fclose(FILE *);

extern void cdInit(void);
extern void cdTalk(bool);
extern void cdOutput(const char *, const char *);
void cdInput(char *, char * );
void cdDepth_depth( double cdDepth[] );
long int cdnZone(void );
double cdB21cm( void );
extern int cdRead(const char *);
void cdPrtWL( FILE *io , realnum wl );
long debugLine( realnum wavelength );
void cdNoExec(void);
int cdDrive(void);
void cdErrors(FILE* );void cdNwcns(
  bool *lgAbort ,
  long int *NumberWarnings, 
  long int *NumberCautions, 
  long int *NumberNotes, 
  long int *NumberSurprises, 
  long int *NumberTempFailures, 
  long int *NumberPresFailures,
  long int *NumberIonFailures, 
  long int *NumberNeFailures );
void cdReasonGeo(FILE*);
void cdWarnings(FILE*);
void cdCautions(FILE*);
void cdSurprises(FILE*);
void cdNotes(FILE*);
long int cdLine(
	const char *chLabel, 
	realnum wavelength, 
	double *relint, 
	double *absint);
void cdLine_ip(long int ipLine, 
	  double *relint, 
	  double *absint );
long int cdDLine(char *chLabel, 
	  realnum wavelength, 
	  double *relint, 
	  double *absint );
int cdColm(const char*, long, double* );
double cdH2_colden( long iVib , long iRot );
double cdCO_colden( long isotope , long iRot );
long int cdEmis(
	char *chLabel,
	realnum wavelength, 
	double *emiss );
void cdEmis_ip(
	long int ipLine, 
	double *emiss );
double cdCooling_last(void);
double cdHeating_last(void);
double cdEDEN_last(void);
void cdPressure_last(
	double *TotalPressure,
	double *GasPressure,
	double *RadiationPressure);
void cdPressure_depth(
	double TotalPressure[],
	double GasPressure[],
	double RadiationPressure[]);
double cdTemp_last(void);
int cdIonFrac(
	const char *chLabel, 
	long int IonStage, 
	double *fracin, 
	const char *chWeight ,
	bool lgDensity );
void cdVersion(char chString[] );
void cdDate(char chString[] );
void cdSetExecTime(void);
double cdExecTime(void);
long int cdGetLineList(
	 char chFile[] ,
	char ***chLabels ,
	realnum **wl );
void cdTimescales(
	double *TTherm , 
	double *THRecom , 
	double *TH2 );
extern long int nFeIIBands;
extern long int nFeIIConBins;
extern realnum **FeII_Bands; 
extern realnum **FeII_Cont; 
void cdSPEC( 
	int Option ,
    double EnergyLow[] , 
    long int nEnergy ,
    double ReturnedSpectrum[] );
void cdSPEC2( 
	int Option ,
    long int nEnergy ,
	long ipLoEnergy,
	long ipHiEnergy,	
    realnum ReturnedSpectrum[] );
int cdTemp(
	const char *chLabel, 
	long int IonStage, 
	double *TeMean, 
	const char *chWeight );
void cdPrintCommands( FILE * );
void cdClosePunchFiles( void );
long int cdH2_Line(
	  /* indices for the upper level */
	  long int iElecHi, 
	  long int iVibHi ,
	  long int iRotHi ,
	  /* indices for lower level */
	  long int iElecLo, 
	  long int iVibLo ,
	  long int iRotLo ,
	  /* linear intensity relative to normalization line*/
	  double *relint, 
	  /* log of luminosity or intensity of line */
	  double *absint );
extern bool lgcdInitCalled;
