#include <TROOT.h>
#include <vector>

#include "ReadConfig.h"

using namespace std;

class ArrayTrigger {

 public:
  ~ArrayTrigger();
  ArrayTrigger(ReadConfig *readConfig,Bool_t debug);
  void SetTelescopeTriggerBitsAndTimes(vector<Bool_t> TelTriggerBit, vector<Float_t> TelTriggerTimes);
  void SetInterTelTransitTimes(vector<float> vTimes);
  bool RunTrigger();
 
  Float_t GetL3DeltaT(){ return fDeltaTL3; };
  Bool_t GetArrayTriggerBitForTelescope(int TelID){ return bTelArrayTriggerBits[TelID]; };
  Float_t GetArrayTriggerTime(int TelID){ return vTelTriggerTimesAfterArrayTrigger[TelID]; };
  Float_t GetEarliestArrayTriggerTime(){ return fEarliestTriggerTime; };

 protected:

  Int_t CalcCluster(Int_t TelID, Int_t ClusterID);
  bool RunTriggerWithNextNeighborRequirement();
  void SetParametersFromConfigFile( ReadConfig *readConfig );

  Bool_t bDebug;

  vector<Bool_t>       bTelArrayTriggerBits;              //stores the information whether a Telescope has been triggered by the array or not
  vector<Bool_t>       bTelTriggerBit;                    //stores the information whether a Telescope has triggered or not
  vector<Float_t>      fTelTriggerTimes;                  //stores the trigger times of each telescope
  Int_t                *iClusterID;                       //Helper variable that assigns a cluster ID to each triggered telescope

  Float_t              *fDeltaL3Cluster;                  //Stores the smallest time difference between L2 signals contributing to an L3 trigger
 
  Float_t              fEarliestTriggerTime;              //The earliest time the array trigger triggered for this event.
 
  Int_t                iNumTel;                           //The number of telescopes
  Int_t                iMultiplicity;                     //How many telescopes need to trigger
  Bool_t               bNextNeighborReq;                  //If the multiplicity requires next neighbors
  Float_t              fCoincidence;                      //the coincidence window between telescope triggers

  Float_t              fDeltaTL3;                         //The time difference between two telescopes L2 signals that triggered L3     

  vector<float>         vTelTransitTimes;                 //the Inter telescope transit times
   vector<float>       vTelTriggerTimesAfterArrayTrigger; //The trigger times of the telescopes after the arraytrigger
  vector< vector < int > > neighbors;            //Has the neigbors of each telescope to be considered in a next neighbor trigger
};
