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

class TriggerTelescopeNextNeighbor {

 public:

  TriggerTelescopeNextNeighbor(ReadConfig *readConfig, int telType, TRandom3 *generator, Bool_t debug = kFALSE, Display *display = NULL);


  void     LoadEvent(TelescopeData *TelData);
  Bool_t   RunTrigger();
  void     RunBiasCurve(UInt_t Trials,Float_t LowerBoundary,Float_t UpperBoundary,Float_t StepWidth,TraceGenerator *tracegenerator,TelescopeData *TelData);

  void     SetDiscriminatorThresholdAndWidth(Float_t threshold, Float_t width);
  void     SetDiscriminatorDelayAndAttenuation(Float_t delay, Float_t attenuation);
  void     SetDiscriminatorRFBConstant(Float_t rfb);
  void     SetDiscriminatorRFBDynamic(Float_t rfb);
  void     SetDiscriminatorRFBUsage(Bool_t rfbuse);


           //Defines how many next neighbors are required to trigger the telescope
  void     SetMultiplicity(Int_t multiplicity){ iMultiplicity = multiplicity; };

  Float_t  GetDiscZeroCrossingRate(){ return lZeroCrossings /(lNumEvents* fTraceLength*1e-3*iNumSumPixGroups); };

  Float_t  GetDiscRFBDynamicValue(){ return fDiscRFBDynamic; };


                   //returns the bias curve and its error
  TVectorD GetBiasCurve(){ return fBiasCurve; };
  TVectorD GetBiasCurveError(){ return fBiasCurveErr; };
  TVectorD GetBiasCurveScanPoints(){ return fBiasCurveScanPoints; };

  TVectorD GetGroupRateVsThreshold(){ return fGroupRateVsThreshold; };
  TVectorD GetGroupRateVsThresholdError(){ return fGroupRateVsThresholdErr; };

  void ShowTrace(int GroupID, bool NSBOnly=kFALSE);
  TH1F GetTraceHistogramCFD(int GroupID);
  TH1F GetTraceHistogramThreshold(int GroupID, bool NSBOnly=kFALSE);

 protected:

  void  SetParametersFromConfigFile(ReadConfig *readConfig );
  
  void CreateTraces();

  Int_t CalcCluster(Int_t GroupID, Int_t ClusterID, Int_t PatchID = -1);

  Bool_t RunDiscriminator(Int_t GroupID);

  Long_t GetNumZeroCrossings();

  Bool_t RunL2WithPatches();

  Bool_t  RunL2Patch(Int_t PatchNumber,Float_t *fPatchTriggerTimes);

  Bool_t GroupInPatch(Int_t GroupID,Int_t PatchID);

  TRandom3 *rand;                   //Our random number generator

  Bool_t bDebug;

  Display *debugDisplay;

  TelescopeData *telData;

  Int_t iTelType;

  Int_t iNumSumPixGroups;              //Holds the number of SumPixGroups
  vector< vector<int> > iSumGroupNeighbors;          //will hold the neighbors of each sumgroup;
  vector< vector<int> > iSumGroupMembers;            //will hold the members of each sumgroup;

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
  Float_t fDiscRFBConstant;            //RFB of the discriminator in pe/MHz
  Float_t fDiscRFBDynamic;             //The dynamic value of the RFB feedback in pe
  Bool_t  bDiscCFDUsage;
  Bool_t  bDiscRFBUsage;
  Float_t fDiscDelay;                  //internal delay of the inverted signal in the disc. in ns
  Float_t fDiscConstantFractionAttenuation;       //constant fraction of the discriminator

  Float_t fClippingLevel;              //The level in mV at which the signals are clipped 
  Bool_t  bDoClipping;                 //Do we clip the signals before summing

  

  Long_t lZeroCrossings;               //holds the number of zerocrossings for one event counted over all pixels in the camera

  Int_t *iClusterID;                   //holds the ClusterID of each sumgroup; -1 if the pixel 
                                       //is not assigned to a cluster
  Int_t iMultiplicity;                 //How many groups need to be in a cluster for a trigger  
 
  //Biascurve variables
  Long_t lNumEvents;                   //holds the number of events filled into traces
  TVectorD fBiasCurve;                 //holds the output of the bias curve 
  TVectorD fBiasCurveErr;              //holds the error of the bias curve 
  TVectorD fBiasCurveScanPoints;
  TVectorD fGroupRateVsThreshold;      //holds the rate vs. Threshold for one group
  TVectorD fGroupRateVsThresholdErr;   //holds the Error of the rate vs. Threshold for one group

  vector<Float_t> *fTracesInSumGroups;         //Stores the signal traces for each sumgroup
  vector<Float_t> *fTracesInSumGroupsConstantFraction;         //Stores the signal traces for each sumgroup
  vector<Float_t> *fTracesInSumGroupsNSBOnly;  //holds the NSB traces for each sumgroup
  Float_t fTraceLength;                   //the length of the simulated trace per group
  Float_t fSamplingTime;                //The sampling rate or resolution of the simulated trace
  Float_t fSamplingTimeAveragePulse;    //The sampling time of the average PE pulse shape

  vector<vector<int> > *vGroupsInCluster;                //pixel that are in one cluster of triggered pixel 
  vector<vector<float> > *vTriggerTimesInCluster;        //the trigger times of all pixels in the cluster


};
