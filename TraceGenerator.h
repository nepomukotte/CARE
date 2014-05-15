
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

  vector<Float_t>  GetLowGainTrace(Int_t PixelID);

  vector<Float_t>  GetTrace(Int_t PixelID,Bool_t bLowGain=kFALSE);

  Float_t GetHighGainAreaToPeak(Int_t pulse){ return fHighGainAreaToPeakConversion[pulse]; };

  void BuildAllHighGainTraces();

  void BuildTrace(Int_t PixelID,Bool_t bLowGain=kFALSE);

  void BuildLowGainTrace(Int_t PixelID);

  void     ShowPulseShape(UInt_t n);
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

           //Set the pulse shapes from an external file
  void     SetPulseShapesFromFile(TString sfilename,Bool_t bLowGain);

            //returns a TH1F histogram with the average single PE pulse shape used in the simulation
  TH1F*    GetHighGainPulse(UInt_t n);
  TH1F*    GetLowGainPulse(UInt_t n);
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
  
  vector< vector<Float_t> > fHighGainPulse;        //Holds the average digitized single PE pulse
                                                   //Max amplitude is normalized to one 20 nsec pulses are simulated
  vector<Float_t> fHighGainStartTime;         //How long before the maximum the pulse is being sampled in nanoseconds
  vector<Float_t> fHighGainStopTime;          //How long after the maximum the pulse is being sampled in nanoseconds
  vector<Float_t> fHighGainAreaToPeakConversion;            
  vector<Float_t> fHighGainLinearAmplitude;
  vector<Float_t> fHighGainNonLinearityFactor;

  Float_t fFWHMofSinglePEPulse;                    //if set to null the pulse shape is read from an external file
  
  Float_t fLinearGainmVPerPE;                     //The mv/pe conversion factor at the input of the FADC in the linear regime

  Float_t fPileUpWindow;                    //The width of the sliding window used to find the number of photons that pile up around any given photon 

  //SiPM related Variable
  Bool_t bSiPM;                                    //flag to figure out if we use SiPMs
  vector<Int_t> vNumCellsPerSIPM;                  //the number of cells in one SiPM
  vector<vector<Bool_t> > vSiPMCellFired;          //boolean array to simulate the dynamic behavior of the SiPM
  vector<Float_t > vSiPMOpticalCrosstalk;          //the optical crosstalk of the SiPM

  //The low gain pulse shape 
  vector< vector<Float_t> > fLowGainPulse;    //The average pulse shape for the low gain channel
  vector<Float_t> fLowGainAreaToPeakConversion;            
  vector<Float_t> fLowGainStartTime;
  vector<Float_t> fLowGainStopTime;
  vector<Float_t> fLowGainLinearAmplitude;
  vector<Float_t> fLowGainNonLinearityFactor;

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

