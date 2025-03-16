#ifndef TriggerTelescopeBase_H_
#define TriggerTelescopeBase_H_


#include <TROOT.h>
#include <vector>
#include <TH1F.h>
#include <TRandom3.h>
#include <TVectorT.h>

#include "TelescopeData.h"
#include "TraceGenerator.h"
#include "ReadConfig.h"
#include "Display.h"

using namespace std;

class TriggerTelescopeBase {

 public:

  TriggerTelescopeBase(ReadConfig *readConfig, int telType, TRandom3 *generator, Bool_t debug = kFALSE, Display *display = NULL);


  virtual void     LoadEvent(TelescopeData *TelData);
  virtual Bool_t   RunTrigger();
  void     RunBiasCurve(UInt_t Trials,Float_t LowerBoundary,Float_t UpperBoundary,Float_t StepWidth,TraceGenerator *tracegenerator,TelescopeData *TelData);

  //configuration of a leading edge and constant fraction discriminator
  void     SetDiscriminatorThresholdAndWidth(Float_t threshold, Float_t width);
  void     SetDiscriminatorDelayAndAttenuation(Float_t delay, Float_t attenuation);


                   //returns the bias curve and its error
  TVectorD GetBiasCurve(){ return fBiasCurve; };
  TVectorD GetBiasCurveError(){ return fBiasCurveErr; };
  TVectorD GetBiasCurveScanPoints(){ return fBiasCurveScanPoints; };

  TVectorD GetTrigPixRateVsThreshold(){ return fTrigPixRateVsThreshold; };
  TVectorD GetTrigPixRateVsThresholdError(){ return fTrigPixRateVsThresholdErr; };

  void ShowTrace(int TrigPixID);
  TH1F GetTraceHistogramCFD(int TrigPixID);
  TH1F GetTraceHistogramThreshold(int TrigPixID);

 protected:

  void  SetParametersFromConfigFile(ReadConfig *readConfig );
  
  Bool_t GroupInPatch(Int_t TrigPixID,Int_t PatchID);

  void CreateTraces();

  Int_t CalcCluster(Int_t TrigPixID, Int_t ClusterID, Int_t PatchID = -1);

  Bool_t RunDiscriminator(Int_t TrigPixID);

  Bool_t RunL2WithPatches();

  Bool_t  RunL2Patch(Int_t PatchNumber,Float_t *fPatchTriggerTimes);


  TRandom3 *rand;                   //Our random number generator

  Bool_t bDebug;

  Display *debugDisplay;

  TelescopeData *telData;

  Int_t iTelType;


  Int_t iNumTriggerPixels;              //Holds the number of trigger pixels
  vector< vector<int> > iTrigPixNeighbors;          //will hold the neighbors of each trigger pixel;
  vector< vector<int> > iTrigPixMembers;            //will hold the members of each trigger pixel;
  vector< vector<int> > iPixelTriggeredInPatch;
  //Conversion factors
  Float_t fFADCSamplingWidth;          //the sampling time of the FADC in ns
  Float_t fFADCconversion;             //the mv to DC conversion factor of the FADC
  Float_t fPEtomVConversion;           //the conversion factor from pe to mV [mV/pe]



  //Pattern trigger
  vector< vector<int> > vPatch;        //!< pattern trigger patches
  Bool_t bUsePatches;                  //do we use patches in L2 yes or no
  
  //Discriminator Related Variables
  Float_t fDiscThreshold;              //Discriminator threshold of pixel
  Float_t fWidthDiscriminator;         //Width of Discriminator output
  Bool_t  bDiscCFDUsage;
  Float_t fDiscDelay;                  //internal delay of the inverted signal in the disc. in ns
  Float_t fDiscConstantFractionAttenuation;       //constant fraction of the discriminator

  Float_t fClippingLevel;              //The level in mV at which the signals are clipped 
  Bool_t  bDoClipping;                 //Do we clip the signals before summing

  
  Int_t *iClusterID;                   //holds the ClusterID of each sumgroup; -1 if the pixel 
                                       //is not assigned to a cluster
  Int_t iMultiplicity;                 //How many groups need to be in a cluster for a trigger  
 
  //Biascurve variables
  TVectorD fBiasCurve;                 //holds the output of the bias curve 
  TVectorD fBiasCurveErr;              //holds the error of the bias curve 
  TVectorD fBiasCurveScanPoints;
  TVectorD fTrigPixRateVsThreshold;      //holds the rate vs. Threshold for one triggerpixel
  TVectorD fTrigPixRateVsThresholdErr;   //holds the Error of the rate vs. Threshold for one group

  vector<Float_t> *fTracesInTriggerPixels;         //Stores the signal traces for each sumgroup
  vector<Float_t> *fTracesInTriggerPixelsConstantFraction;         //Stores the signal traces for each sumgroup
  vector<Float_t> *fTracesInTriggerPixelsNSBOnly;  //holds the NSB traces for each sumgroup
  Float_t fTraceLength;                   //the length of the simulated trace per group
  Float_t fStartSamplingBeforeAverageTime; //the offset from the average photon arrival time when the trace starts to be sampled
  Float_t fSamplingTime;                //The sampling rate or resolution of the simulated trace
  Float_t fSamplingTimeAveragePulse;    //The sampling time of the average PE pulse shape

  vector<vector<int> > *vTrigPixsInCluster;                //pixel that are in one cluster of triggered pixel 
  vector<vector<float> > *vTriggerTimesInCluster;        //the trigger times of all pixels in the cluster


};
#endif // TriggerTelescopeBase_H_

