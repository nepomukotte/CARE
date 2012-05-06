
#ifndef Trace_Generator
#define Trace_Generator

#include <TROOT.h>
#include <TRandom3.h>
#include <TH1F.h>		
#include <vector>
#include <iostream>
#include <sstream>
#include "TelescopeData.h"
#include "ReadConfig.h"
#include "Display.h"
#include "GOrderedGrid.h"
#include "GOrderedGridSearch.h"


using namespace std;

class TraceGenerator {

 public:
 
  TraceGenerator(ReadConfig *readConfig, int telType,  TRandom3 *generator, Bool_t debug = kFALSE,  Display *display = NULL);
  void             LoadCherenkovPhotons( std::vector< float > *v_f_X,std::vector< float > *v_f_Y,std::vector< float > *v_f_time,std::vector< float > *v_f_lambda, Float_t delay, Double_t dEfficiencyFactor);      //Loads Photons into traces
  void             GenerateNSB();
  Float_t          GetAreaToPeakConversion(){return fAreaToPeakConversion;};

  vector<Float_t>  GetLowGainTrace(Int_t PixelID);


  void     ShowAverageSinglePEPulse();
  //Need to generate vector with times for each sample


  void     SetTelData(TelescopeData *telData);

  void     PrintHowOftenTheTraceWasTooShort();

 protected:


  void     SetParametersFromConfigFile(ReadConfig *readConfig );
  
  void     AddPEToTrace(Int_t PixelID, Float_t time, Int_t NumPE = 1);
  
  Int_t    SimulateSiPM(Int_t PixelID, Int_t NumPE);

  void     ResetSiPM(); 

                                                            
  void     SetGaussianPulse(Float_t fwhm);
  void     SetSigmaofSinglePEPulseHeightDistribution(Float_t sigma);
  void     SetSinglePESamlingWidth(Float_t width);
  void     SetTraceSampleWidthAndLength(Float_t t,Float_t length=100,Float_t=30);
 
  
            //Sets the NSB rate per pixel, units kHz  
  void     SetNSBRatePerPixel(Float_t rate){ fNSBRatePerPixel = rate; };

           //Use NSB 
  void     SetNSBUsage(Bool_t use){ bUseNSB = use;};

           //Defines if afterpulsing is used in the simulations, Fit function exp(a+b*x), where a=constant b=slope
  void     SetAfterPulsing(Bool_t afterpulsing, Float_t constant = -1, Float_t slope = -1);

           //Sets the Factor by which the Cherenkov Photons need to be scaled down
  void     SetWinstonConeEfficiency(Float_t CherenkovPhotonScaling);
  
           //Set the single pe pulse shape from an external file
  void     SetSinglePEShapeFromFile(TString sfilename);

           //Set the low gain single pe pulse shape from an external file
  void     SetLowGainSinglePEShapeFromFile(TString sfilename);

            //returns a TH1F histogram with the average single PE pulse shape used in the simulation
  TH1F*    GetAverageSinglePEPulse();
  TH1F*    GetAverageSinglePELowGainPulse();
  //Variables

  TelescopeData *telData;

  Bool_t bDebug;
  Display *debugDisplay;
  TRandom3 *rand;                                  //Our random number generator

  GOrderedGridSearch *gridsearch;

  Int_t  iTelType; 
 
  Int_t NumTraceIsTooShortLowEnd;
  Int_t NumTraceIsTooShortHighEnd;

 
  //PMT signal related variables
  Float_t  fSigmaSinglePEPulseHeightDistribution;  //sigma of the single PE pulse height distribution
  
  vector<Float_t> fAverageSinglePEPulse;           //Holds the average digitized single PE pulse
                                                   //Max amplitude is normalized to one 20 nsec pulses are simulated
  Float_t fAreaToPeakConversion;                   //the conversion factor from area to peak value for a single pe pulse.
  Int_t iNumSamplesAverageSinglePEPulse;           //the number of samples in the average pulse  
  Float_t fStartTimeAveragePulse;                  //How long before the maximum the pulse is being sampled in nanoseconds
  Float_t fStopTimeAveragePulse;                   //How long after the maximum the pulse is being sampled in nanoseconds
  Float_t fFWHMofSinglePEPulse;                    //if set to null the pulse shape is read from an external file
  

  //SiPM related Variable
  Bool_t bSiPM;                                    //flag to figure out if we use SiPMs
  vector<Int_t> vNumCellsPerSIPM;                  //the number of cells in one SiPM
  vector<vector<Bool_t> > vSiPMCellFired;          //boolean array to simulate the dynamic behavior of the SiPM
  vector<Float_t > vSiPMOpticalCrosstalk;          //the optical crosstalk of the SiPM

  //The low gain pulse shape 
  vector<Float_t> fAverageLowGainSinglePEPulse;    //The average pulse shape for the low gain channel
  Float_t fLowGainAreaToPeakConversion;            
  Float_t fLowGainStartTimeAveragePulse;
  Float_t fLowGainStopTimeAveragePulse;
  Int_t iNumSamplesLowGainAverageSinglePEPulse;

  Float_t         fSamplingTime;                           //The sampling rate or resolution of the simulated trace
  Float_t         fSamplingTimeAveragePulse;               //The sampling time of the average PE pulse shape
  Float_t         fTraceLength;                            //the length of the simulated trace per group
  Float_t         fStartSamplingBeforeAverageTime;         //The offset from the average photon arrival time, when the trace gets sampled

  //Afterpulsing
  Float_t fAPconstant;                             // values for afterpulsing from rate vs. threshold curve
                                                   //Fit function exp(a+b*x), where a=constant b=slope
  Float_t fAPslope;
  Bool_t  bAfterPulsing;                           //Do we simulate Afterpulsing
 
  //Crosstalk between camera pixel
  Bool_t  bCrosstalk;                    //Use crosstalk between pixel
  Float_t fCrosstalk;                    //Crosstalk value


  //Camera layout
  Int_t iNumPixels;                    //Number of pixels in the camera 
  vector<vector<int> > vNeighbors;     //neighbors of each pixel 

  vector<Float_t> qe;                  //holds the QE values. the index number of each entry gives the wavlength, 
                                       //i.e. qe[400] is the qe at 400 nm


  //NSB
  Float_t fNSBRatePerPixel;                        //the NSB rate kHz per mm squared in the focal plane;
  Bool_t  bUseNSB;

  //Shower photons
  Float_t fWinstonConeEfficiency;                 //The efficiency of the Winstoncone

  Float_t fCherenkovPhotonScaling;                //Fraction of photons that have been recorded into the PE file


  vector<double>  fXTubeMM;
  vector<double>  fYTubeMM;
  vector<double>  fSizeTubeMM;
  vector<double>  fRotAngle;
  vector<int>     iTubeSides;



};

#endif

