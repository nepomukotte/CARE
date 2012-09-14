
#ifndef Telescope_Data
#define Telescope_Data

#include <TROOT.h>
#include <TRandom3.h>
#include <TH1F.h>		
#include <vector>
#include <iostream>
#include <sstream>

#include "ReadConfig.h"

using namespace std;

class TelescopeData {

 public:
 

  TelescopeData( ReadConfig *readConfig, Int_t telTID = 0, Bool_t debug = kFALSE);

  //Trace related Functions
  Float_t          GetAverageArrivalTime(){return fAveragePhotonArrivalTime;};


  //Trigger related functions
  Bool_t           CherenkovPhotonsInCamera(){return bCherenkovPhotonsInCamera;};

  //Display related functions
  vector<Int_t>*   GetDisplayTraces(){ return iDisplayTraceInPixel; };
  vector<Int_t>    GetDisplayTrace(Int_t pixelID){return iDisplayTraceInPixel[pixelID]; };
  Bool_t           GetPixelLowGainSwitch(Int_t pixelID){ return bInLoGain[pixelID]; };

  vector<Int_t>    GetFADCTrace(Int_t pixelID){ return iFADCTraceInPixel[pixelID]; };
  vector<Int_t>    GetQDCValues(){return iQDCInPixel; };  

  vector<Int_t>    GetPEInPixels(){return iPEInPixel; };  

  //Telescope related stuff
  Int_t            GetTelescopeType(){return iTelType;};
  Int_t            GetTelescopeID(){return iTelID; };
  Int_t            GetTelescopeIDinSuperArray(){return iTelIDinSuperArray; };

  //trigger stuff
                   //returns vector with the trigger bits of each group
  vector< Bool_t>  GetTriggeredGroups(){ return bTriggeredGroups; };
  Int_t            GetNumTriggeredGroups(){ return iNumTriggeredGroups; };
                   //time when telescope has triggered
  Float_t          GetTelescopeTriggerTime(){ return fTelescopeTriggerTime; };
                   //Returns true if telescope has triggered, false otherwise
  Bool_t           GetTelescopeTrigger(){ return bTelescopeHasTriggered; };
  Bool_t           GetGroupTrigger(Int_t iGroupID){return  bTriggeredGroups[iGroupID]; };
                   //Returns a vector with a list of Groups that are in the Trigger cluster
  vector<int>      GetTriggerCluster(){ return vTriggerCluster; };

  //general functions to maintain object
  void ResetTraces(); //Sets all vectors and numbers to initial values.



  //////////////////////////////////////////////////////////////////////////////////////////////////////
 
  // Variables



  //Trace related parameters  analog signal
  vector<Float_t> *fTraceInPixel;
  vector<Float_t> *fTraceInPixelNSBOnly;
  vector<Float_t> *fTimesInPixel;
  vector<Float_t> *fAmplitudesInPixel;
  Float_t         fAveragePhotonArrivalTime;               //Holds the average photon arrival time of all photons in one event
  Double_t        mean;                                     //trace mean
  Bool_t          bCherenkovPhotonsInCamera;

  Int_t           iNumPixels;
  Int_t           iNumSamplesPerTrace;                     //the number of samples in the analog trace for each sum group

  //FADC
  vector<Int_t>   *iFADCTraceInPixel;
  Int_t           iNumFADCSamples;

  //Display trace
  vector<Bool_t>  bInLoGain;
  vector<Int_t>   *iDisplayTraceInPixel;
  Int_t           iNumDisplaySamples; 

  //QDC
  vector<Int_t>   iQDCInPixel;

  vector<Int_t>   iPEInPixel;                          //the number of Cherenkov photoelectrons in each pixel

  //Trigger related numbers
  Float_t         fDiscriminatorThreshold;              
  Float_t         fTriggerTime;    
  vector<Bool_t> bTriggeredGroups;                      //holds the information whether a group has triggered or not;
  vector<Float_t> fDiscriminatorTime;                   //holds the time when the group has triggered;
  Float_t fTelescopeTriggerTime;                        //The time when the telescope triggered
  Int_t iNumTriggeredGroups;                            //holds the number of triggered groups

  Bool_t bTelescopeHasTriggered;                        //If telescope has triggered
  vector<int> vTriggerCluster;                           //the IDs of the groups that are in the triggered cluster

  //Pixel related variables
  vector< Float_t > fRelQE;
  vector< Float_t > fRelQEwWC;                           //with Winston cone efficiency
  vector< Float_t > fRelGain;

  //Blur optical PSF and optical efficieny
  Float_t fWinstonConeEfficiency;                       //The efficiency of the Winstoncone
  Bool_t  bBlurPSF;                      
  Float_t fBlurSigma;                                   //sigma in mm by which the optical PSF 
                                                        //is blured furthery of the Winstoncone
  //Other Telescope Data
  Float_t TelXpos;                                      //x coordinate of the telescope position, comes from the photon input file  
  Float_t TelYpos;                                      //y coordinate of the telescope position, comes from the photon input file
  Float_t TelZpos;                                      //y coordinate of the telescope position, comes from the photon input file
  Float_t OpticsTransitTime;                            //Transit time for the photons through the optics
private:

  void     SetParametersFromConfigFile(ReadConfig *readConfig );

  void SetupArrays();   //initializes all arrays

  Bool_t bDebug;

  Int_t iTelID;
  Int_t iTelIDinSuperArray;
  Int_t iTelType;   //What telescope type we have

};

#endif

