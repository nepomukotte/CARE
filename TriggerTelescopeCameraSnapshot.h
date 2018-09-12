#include "TriggerTelescopeNextNeighbor.h"

using namespace std;

class TriggerTelescopeCameraSnapshot : public TriggerTelescopeNextNeighbor {

 public:

  TriggerTelescopeCameraSnapshot(ReadConfig *readConfig, int telType, TRandom3 *generator, Bool_t debug = kFALSE, Display *display = NULL);
  void     SetDiscriminatorThresholdAndWidth(Float_t threshold, Float_t width);
  Bool_t   RunTrigger();
  void     RunBiasCurve(UInt_t Trials,Float_t LowerBoundary,Float_t UpperBoundary,Float_t StepWidth, TraceGenerator *tracegenerator,TelescopeData *TelData);

 protected:

  // Set the internal parameters using *readConfig (overwrite TriggerTelescopeNextNeighbor::SetParametersFromConfigFile)
  void  SetParametersFromConfigFile(ReadConfig *readConfig );
  // Set the internal parameters through *readConfig (overwrite TriggerTelescopeNextNeighbor::SetFromConfigFile)
  void  SetTriggerLogicSettings( ReadConfig *readConfig );
  // Digitized amplitudes at iPositionInAnalogTrace in the traces, find clusters and discriminate them.
  //void  DiscriminateCameraSnapshot( Int_t iPositionInAnalogTrace );
  vector<Int_t> DiscriminateCameraSnapshot( Int_t iPositionInAnalogTrace );
  //Check if the telescope is ready to fire. Return false if not enough snapshots (yet) checked or if the trigger has already fired!
  Bool_t  ReadyToFire();
  // Check if the trigger is firing, using the selected 'iComboMode' mode
  Bool_t  TriggerIsFiring();
  //Check if the 'telData' parameters so far set are consistent before using a mode. This method is meant to be used insise a mode (following) as sanity checks
  Bool_t  IsConsistent(); //Is it necessary?
  //iComboMode==0: if consecutive 'iSamplingsWindow's have at least a matching dicriminated clusters, the trigger fires!
  Bool_t  BlindMode();
  //iComboMode==1: if edges of 2 readoout windows 'iSamplingsWindow' long with discriminated clusters are contained inside each other, the trigger fires!
  Bool_t  EdgeMode();
  //iComboMode==2: ...not sure...
  Bool_t  LevelMode();

  //Special snapshot parameters
  Int_t  iResolution;                  //number of bit resolution for the camera snapshot
  Int_t  iScalingDivisor;              //the scaling divisor to reduce the ADC counts of the group sum
  Int_t  iOffset;                      //offset t obe added to the digitized signals (in ADC counts)
  Int_t  iCircles;                     //number of circles to form the pattern for the camera snapshot (SST-1M)
  UInt_t uNeighbors;                   //number of neighboirs to let a group to form a pattern for the camera snapshot (SST-1M)
  //Int_t  iNslices;                     //how many slices are expected in the trace length
  Int_t  iComboMode;                   //combination type: 0=sequential samplings, 1=TBD, 2=TBD
  Int_t  iSamplingsWindow;             //Samplings window
  //Single PE pulse shape parameters
  Float_t fPEtoDCconversion;                       //the gain calibration constant DC per PE; defined for the amplitude of a single pe pulse

  //FADC parameters
  Float_t fFADCSamplingWidth;                      //the sampling time of the FADC in ns
  Float_t fFADCconversion;                         //the mv to DC conversion factor of the FADC
  Int_t   iResolutionRange;                        //the counts dynamic range of the digitizer: Power(2,iResolution)-1
  Float_t fHiLoGainThreshold;                      //dc counts at which the HiLoGain switch is activated
  Float_t fLowHiGainRatio;                         //Gain ratio logain/higain

  //vector < vector < vector<int> > > vDiscriminatedClustersInSnapshots; //Container of found clusters (of groups) for snapshots: vClustersInSnapshots[iSnapshot][iCluster][iGroup]
  
};
