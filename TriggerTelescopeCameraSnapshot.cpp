/* \file TriggerTelescopeCameraSnapshot.cpp
   Implementation of a camera snapshot trigger for SST-1M
   1. Pulses sampled every FADCsample ns (a snapshot)
   2. In a snapshot, signals of X pixels (a patch) are summed, digitized and clipped with a dynamic range of 2^B-1 (B in bits)
   2. In a snapshot, the sum of patches in a cluster (built through surrounding neighbors from a central one,
      cicle by cicler for C circles) is discriminated
   3. Combinig snapshots together:
     3a. Blind mode: If found the same discirminated cluster for N consecutive snapshots, the trigger fires
     3b. ...
*/

#include "TriggerTelescopeCameraSnapshot.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>

#include <TMath.h>
#include <TTimer.h>
#include <Getline.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLine.h>

using namespace std;

//---------------------------------------------------------------------------------------
//Constructor
TriggerTelescopeCameraSnapshot::TriggerTelescopeCameraSnapshot(ReadConfig *readConfig, int telType, TRandom3 *generator,Bool_t debug,Display *display)
  : TriggerTelescopeNextNeighbor(nullptr, telType, generator, debug, display)
{
  //trigger logic
  iCircles = 0;
  uNeighbors = 0;
  iComboMode = -1;
  iSamplingsWindow = 0;

  //signal digitization
  fPEtoDCconversion=-1;
  fFADCSamplingWidth=-1;
  fFADCconversion=-1;
  fHiLoGainThreshold=-1;
  fLowHiGainRatio=-1;
  iResolution = -1;
  iResolutionRange = 0;
  iScalingDivisor = -1;
  iOffset = 0;
  
  cout<<"You just bought an excellent camera snapshot trigger"<<endl;
  SetParametersFromConfigFile(readConfig);
}


//-----------------------------------------------------------------------------------------------------------------------
//This does the actual triggering
//It returns true if the telescope has triggered and false if not
Bool_t  TriggerTelescopeCameraSnapshot::RunTrigger()
{
  
  //Get all the zerocrossings for the RFB
  lZeroCrossings += GetNumZeroCrossings();
  lNumEvents++;
  
  //Set the updated RFB feedback. Does it need to be adapted?
  if(bDiscRFBUsage)
    {
      if(bDebug)
	    cout<<lNumEvents<<" RFB dynamic value  "<<fDiscRFBConstant * lZeroCrossings /(lNumEvents* (fTraceLength-fDiscDelay)*1e-3*iNumSumPixGroups)<<endl;
      SetDiscriminatorRFBDynamic(fDiscRFBConstant * lZeroCrossings /(lNumEvents* (fTraceLength-fDiscDelay)*1e-3*iNumSumPixGroups));
    }
  
  //Reset trigger parameters
  telData->iNumTriggeredGroups = 0;
  telData->bTelescopeHasTriggered = kFALSE;
  telData->vTriggerCluster.clear();
  telData->vSnapshotsDiscriminatedClusters.clear();
  Int_t  iPositionInAnalogTrace;

  for( Int_t i=0; i<telData->iSnapshots; i++ )
  {
    iPositionInAnalogTrace = (Int_t) ( i * fFADCSamplingWidth / fSamplingTime );
    if ( iPositionInAnalogTrace >= telData->iNumSamplesPerTrace)
	{
	  break;
	}
      telData->iSnapshotsDiscriminatedGroups[i].clear();
      telData->iSnapshotsDiscriminatedGroups[i] = DiscriminateCameraSnapshot( iPositionInAnalogTrace );
      if ( TriggerIsFiring() )
	  {
	    telData->iNumTriggeredGroups = telData->iSnapshotsDiscriminatedGroups[i].size();
	    telData->fTelescopeTriggerTime = (i-iSamplingsWindow+1) * fFADCSamplingWidth + telData->fAveragePhotonArrivalTime - fStartSamplingBeforeAverageTime;
	  }
  }

  //debug
  if(telData->bTelescopeHasTriggered && bDebug)
    {
      cout<<"Event triggered telescope"<<endl; 
      cout<<"at time "<< telData->fTelescopeTriggerTime<<endl;
    }
  
  return telData->bTelescopeHasTriggered;
}

//To doscriminate the clusters in a snapshot
vector<Int_t> TriggerTelescopeCameraSnapshot::DiscriminateCameraSnapshot( Int_t iPositionInAnalogTrace )
{
  Int_t iClusterDigitalSum;
  Float_t fGain = 1.;  //For future introduction of the low gain option
  vector < vector<int> > vDiscriminatedClusters; //container of the discriminated clusters found for this snapshot 
  //Find all possible clusters above the discriminator threshold (in dc counts!)
  for(Int_t group=0; group<iNumSumPixGroups; group++)
    {
      //DigiCam FIRST digitize the group, THEN sum: if not "fFADCconversion"ed each group signal independently, the total sum might have different approximation...
      //Since the conversion is operated on a negative pulse shape, pedestal-signal*conversion is the correct procedure.
      //Then must be devided by the number of pixels in the group, because the scaling factor (pedestal hsould be already scaled!) is per group...
      Int_t gADC = (Int_t)(iOffset-fGain*fTracesInSumGroups[group][iPositionInAnalogTrace]*fFADCconversion);
      gADC = gADC<iResolutionRange ? gADC : iResolutionRange;
      iClusterDigitalSum = gADC>0 ? gADC : 0;
      vector<int> vGroupsInCluster; //groups forming the cluster built around 'group'
      vGroupsInCluster.push_back(group); //...at least 'group' is inside the cluster!
      //Look for surrounding groups of the central one to form a cluster or iCircles circles
      Int_t circles = iCircles>0 ? iCircles : 0; //how many circles must be in the cluster? Note: at least 0, iCircles<0 should NEVER happen!
      Bool_t bWrongConditions = kFALSE; //initial assumption: everything will be fine!
      Bool_t bClusterFound = kTRUE; //initial assumption: a cluster will be found
      while ( circles-- ) //build the pattern circle by circle. If iCircle is 0, this loop is skipped!
	{
	  Int_t iCurrentClusterSize = vGroupsInCluster.size(); //how many groups must be cheched to add this new external circle?
	  for ( int i=0; i<iCurrentClusterSize; i++ ) //check neighbors of each group already part of the cluster to enlarge it of one circle
	    {
	      Int_t g = vGroupsInCluster[i]; //group ID currently to be checked
	      if ( iSumGroupNeighbors[g].size() != uNeighbors ) //minimum condition: at least SNAPSHOTSURROUNDINGNEIGHBORS neighbors of the group 'g'
		{
		  bWrongConditions = kTRUE; //no cluster can be formed starting from 'group' as central group becasue the minimum condition!
		  break;
		}
	      for ( UInt_t n=0; n<uNeighbors; n++ ) //looping over the 'g' neighbors
		{
		  Int_t neighbor = iSumGroupNeighbors[g][n];
		  if ( find(vGroupsInCluster.begin(), vGroupsInCluster.end(), neighbor) == vGroupsInCluster.end() ) //if 'neighbor' is NOT in already in 'vGroupsInCluster'...
		    {
		      //...sum up the signal and include the 'neighbor' in the cluster
		      Int_t gNeighADC = (Int_t)(iOffset-fGain*fTracesInSumGroups[neighbor][iPositionInAnalogTrace]*fFADCconversion);
		      gNeighADC = gNeighADC<iResolutionRange ? gNeighADC : iResolutionRange;
		      iClusterDigitalSum += gNeighADC>=0 ? gNeighADC : 0; //empty signals might have a bad memory allocation: give positive large numbers => iGroupADC<0
		      vGroupsInCluster.push_back(neighbor); //the found neighbor is part of the vGroupsInCluster vector, now!
		    }
		}
	    }
	  if ( bWrongConditions ) //Strating from the current central group is not possible to form the wanted cluster => break the loop!
	    {
	      bClusterFound = kFALSE; //no cluster can be found!
	      break;
	    }
	}
      if ( bClusterFound ) //if false, means that bWrongConditions was true for this group => nothing to do, go to the next!
	{
	  iClusterDigitalSum = (Int_t)(iClusterDigitalSum/vGroupsInCluster.size()); //Rescaling resolution of the digitized groups
	  iClusterDigitalSum = iClusterDigitalSum>iResolutionRange ? iResolutionRange : iClusterDigitalSum;
	  if ( iClusterDigitalSum >= (Int_t)fDiscThreshold ) //The cluster is discriminated...
	    {
	      sort( vGroupsInCluster.begin(), vGroupsInCluster.end() ); //...sort the group IDs in the cluster for a better comparison of clusters between snapshosts
	      vDiscriminatedClusters.push_back(vGroupsInCluster); //...include this cluster int the clusters container for this snapshot
	    }
	}
    }
  telData->vSnapshotsDiscriminatedClusters.push_back(vDiscriminatedClusters); //no clusters found or none discriminated => empty array for this snapshot
  vector<int> vDiscriminatedGroups; //array of all triggered groups in this snapshot
  for ( vector< vector<Int_t> >::iterator itCluster=vDiscriminatedClusters.begin();
	itCluster!=vDiscriminatedClusters.end(); ++itCluster ) //skipped if no clusters found or none discriminated
    {
      for ( vector<Int_t>::iterator itGroup=itCluster->begin();
	    itGroup!=itCluster->end(); ++itGroup ) //At least a group is inside...
	{
	  if ( find(vDiscriminatedGroups.begin(), vDiscriminatedGroups.end(), *itGroup) == vDiscriminatedGroups.end() ) //'vTriggeredGroups' must contain unique values...
	    {
	      vDiscriminatedGroups.push_back(*itGroup);//...if 'groupID is not already there, push back!
	    }
	}
    }
  sort( vDiscriminatedGroups.begin(), vDiscriminatedGroups.end() ); //sort 'vTriggeredGroups' for a better show in a root file
  return vDiscriminatedGroups;
}

//Combining the snapshot so far accumulated according to the selected mode 
Bool_t   TriggerTelescopeCameraSnapshot::TriggerIsFiring()
{
  if ( ReadyToFire() ) //if the telescope is already triggered, no need to check again
    {
      Bool_t bTriggered;
      switch ( iComboMode )
	{
	case (0) :
	  bTriggered = BlindMode();
	  break;
	case (1) :
	  bTriggered = EdgeMode();
	  break;
	case (2) :
	  bTriggered = LevelMode();
	  break;
	default: //This must NEVER happen!
	  cout << "Unknown mode is requested! This forced exit should never happen because and unknown mode should be set by default to 0 in the ReadConfig class!" << endl;
	  exit(1);
	}
      if ( bTriggered )
	{
	  telData->bTelescopeHasTriggered = kTRUE;
	}
      return bTriggered;
    }
  return kFALSE; //if snapshots not yet enough or telescope is already triggered (triggered paramenter already set => forbit their overwriting!)
}

//Blind mode logic to combine the snapshots together
Bool_t   TriggerTelescopeCameraSnapshot::BlindMode()
{
  if ( !IsConsistent() )
    {
      return kFALSE;
    }
  Int_t iLatest = telData->vSnapshotsDiscriminatedClusters.size()-1;
  vector< vector<Int_t> > vLatestSnapshot = telData->vSnapshotsDiscriminatedClusters[iLatest];
  Bool_t bMatchingCluster = kTRUE; //if vLatestSnapshot is empty, doesn't care, if iSamplingsWindow==1, must be true for each cluster in vLatestSnapshot
  Bool_t bMatchingGroup; //If not filled at least once, there a problem: good if an exception occurs!
  for ( vector< vector<Int_t> >::iterator itCluster=vLatestSnapshot.begin();
	itCluster!=vLatestSnapshot.end(); ++itCluster ) //Loop over the latest dicriminated clusters
    {
      for ( Int_t iBack=1; iBack<iSamplingsWindow; iBack++ ) //From 1 snapshot to iSamplingsWindow-1 back
	{
	  vector< vector<int> > vBackSnapshot = telData->vSnapshotsDiscriminatedClusters[iLatest-iBack];
	  bMatchingCluster = kFALSE; //Not yet found in this previous snapshot
	  for ( vector< vector<Int_t> >::iterator itBackCluster = vBackSnapshot.begin();
		itBackCluster!=vBackSnapshot.end(); ++itBackCluster ) // loop over the dicriminated clusters in the previous snapshot
	    {
	      for ( vector<Int_t>::iterator itGroup = itCluster->begin();
		    itGroup!=itCluster->end(); ++itGroup) // loop over the groups in the dicriminated cluster of the latest snapshot. AT LEAST 1 is there!
		{
		  bMatchingGroup = find( itBackCluster->begin(), itBackCluster->end(), *itGroup ) != itBackCluster->end();
		  if ( !bMatchingGroup ) //If a group is not in itBackCluster, the cluster is not matched: try next cluster in the previous snapshot!
		    {
		      break;
		    }
		}
	      if ( bMatchingGroup ) //If all groups are found, the *itCluster is found. It cannot be found more then once!
		{
		  bMatchingCluster = kTRUE;
		  break;
		}
	    }
	  if ( !bMatchingCluster ) //If no matching cluster is found (bMatchingCluster==kFALSE) in this previous snapshot, no trigger! Pointelss to look further back...
	    {
	      break;
	    }
	}
      if ( bMatchingCluster )
	{
	  for ( vector<Int_t>::iterator itGroup = itCluster->begin();
		itGroup!=itCluster->end(); ++itGroup ) //Loop over the contained groups to fill the 'vTriggerCluster' array
	    {
	      if ( find(telData->vTriggerCluster.begin(), telData->vTriggerCluster.end(), *itGroup) == telData->vTriggerCluster.end() ) //Group already stored?
		{
		  telData->vTriggerCluster.push_back(*itGroup); //Array filled with unique group IDs (only when trigger is fired for the first time!)
		}
	    }
	}
    }
  sort( telData->vTriggerCluster.begin(), telData->vTriggerCluster.end() ); //The groups IDs are shown sorted in the root file!
  return telData->vTriggerCluster.size()>0 ? kTRUE : kFALSE;
}

//Edge mode logic to combine the snapshots together
Bool_t   TriggerTelescopeCameraSnapshot::EdgeMode()
{
  if ( !IsConsistent() )
    {
      return kFALSE;
    }
  cout << "'TriggerTelescopeCameraSnapshot::EdgeMode()' is not yet available"<<endl;
  exit(1);
}

//Level mode logic to combine the snapshots together
Bool_t   TriggerTelescopeCameraSnapshot::LevelMode()
{
  if ( !IsConsistent() )
    {
      return kFALSE;
    }
  cout << "'TriggerTelescopeCameraSnapshot::LevelMode()' is not yet available"<<endl;
  exit(1);
}

//Check if the data stored so far are consistent with the structure of the logical operation. If everything is used consistently, it returns kTRUE always
Bool_t   TriggerTelescopeCameraSnapshot::ReadyToFire()
{
  if ( telData->vSnapshotsDiscriminatedClusters.size()<(UInt_t)iSamplingsWindow || //if there are not at least iSamplingsWindow snapshots stored, not ready to fire! 
       telData->bTelescopeHasTriggered ) //trigger already fired: cannot fire anymore!
    {
      return kFALSE;
    }
  if ( telData->vTriggerCluster.size()>0 ) //Consistency condition: if the method is used correctly, this check is not necessary!
    {
      cout<<"The trigger cluster ('vTriggerCluster': groups that contributed to fire the trigger) is not empty!"<<endl;
      cout<<" ...Forgot to reset it? Reset (wrongly) somewhere 'bTelescopHasTerigger' bit?"<<endl;
      cout<<"    -> Forced reset of the vTriggerCluster vector, but it is better to check what is going on!!!"<<endl;
      telData->vTriggerCluster.clear();
    }
  return kTRUE;
}

//Check if the data stored so far are consistent with the structure of the logical operation. If everything is used consistently, it returns kTRUE always
Bool_t   TriggerTelescopeCameraSnapshot::IsConsistent()
{
  if ( telData->vSnapshotsDiscriminatedClusters.size()<(UInt_t)iSamplingsWindow ) // Sanity check! At this point, the number of sotred snapshots MUST BE AT LEAST iSamplingsWindow...
    {
      cout << "To check a logic, the snapshots stored MUST BE AT LEAST "<<iSamplingsWindow
	   <<", but they are "<<telData->vSnapshotsDiscriminatedClusters.size()
	   <<"! THIS SHOULD NEVER HAPPEN!" << endl;
      exit(1);
    }
  if ( telData->bTelescopeHasTriggered ) //Consistency condition: if the method is used correctly, this check is not necessary!
    {
      cout<<"Trigger already fired for this event! No check necessary... NOTE: Are you using this method wrong???"<<endl;
      return kFALSE;
    }
  if ( telData->vTriggerCluster.size()>0 ) //Consistency condition: if the method is used correctly, this check is not necessary!
    {
      cout<<"The trigger cluster ('vTriggerCluster': groups that contributed to fire the trigger) is not empty!"<<endl;
      cout<<" ...Forgot to reset it? Reset (wrongly) somewhere 'bTelescopHasTerigger' bit?"<<endl;
      cout<<"    -> Forced reset of the vTriggerCluster vector, but it is better to check what is going on!!!"<<endl;
      telData->vTriggerCluster.clear();
    }
  return kTRUE;
}

//Reads in  the config file and sets all variables
void   TriggerTelescopeCameraSnapshot::SetParametersFromConfigFile(ReadConfig *readConfig )
{
  cout <<endl<< "TriggerTelescopeCameraSnapshot::SetParametersFromConfigFile " << endl;
  
  //Set configuration commond to any trigger logic
  SetCommonSettings(readConfig);
  
  //Set specific configuration for the trigger logic: new trigger logic (= child of this class) might require it's own SetTriggerSettings
  SetTriggerLogicSettings(readConfig);

  //Create the traces for the summed pixels
  CreateTraces();
}

void   TriggerTelescopeCameraSnapshot::SetTriggerLogicSettings( ReadConfig *readConfig )
{

  //Forced variables
  if ( bDiscCFDUsage )
    {
      bDiscCFDUsage = kFALSE;
      cout<<"The use of CFD is forced to be false for the camera snapshot logic: signal is digitized at trigger level!"<<endl;
    }
  //The number of circles
  iCircles = readConfig->GetSnapshotCircle(iTelType);   
  cout<<"The number of circles to find the pattern: "<<iCircles<<endl;

  //Number of neigbors
  uNeighbors = readConfig->GetSnapshotNeighbors(iTelType);   
  cout<<"The number of neighbors to let a group to form a pattern: "<<uNeighbors<<endl;

  //Combination type
  iComboMode = readConfig->GetSnapshotComboMode(iTelType);   
  cout<<"The combo mode selected (ID number): "<<iComboMode<<endl;

  //Samplings window
  iSamplingsWindow = readConfig->GetSnapshotSamplingWindow(iTelType);   
  cout<<"The samplings window: "<<iSamplingsWindow<<endl;

  //the sampling width of the FADC
  fFADCSamplingWidth = readConfig->GetFADCSampleTime(iTelType);
  cout<<"Digital sampling width in ns: "<<fFADCSamplingWidth<<endl;
  
  //The resolution in bits
  iResolution = readConfig->GetSnapshotBitsResolution(iTelType);
  iResolution = iResolution>0 ? iResolution : 1;
  cout<<"The snapshoot resolution in bits (reminder: <=0 means forced to 1!): "<<iResolution<<endl;
  //the dynamic range of the FADC
  iResolutionRange = (Int_t)TMath::Power(2,iResolution)-1;
  cout<<"The dynamic range of the trigger resolution in ADC: "<<iResolutionRange<<endl;

  //the scaling divisor
  iScalingDivisor = readConfig->GetSnapshotScalingDivisor(iTelType);
  iScalingDivisor = iScalingDivisor>0 ? iScalingDivisor : 1;
  cout<<"The rescaling divisor (reminder: <=0 means forced to 1!): "<<iScalingDivisor<<endl;

  //the offset (pedestal) to be addes to the digitized signal of a group in DC
  iOffset = readConfig->GetSnapshotFADCOffset(iTelType);
  iOffset = iOffset>0 ? iOffset : 0;
  cout<<"The offset (pedestal) to be addes to the digitized signal of a group in DC (reminder: <0 means forced to 0!): "<<iOffset<<endl;

  //the conversions between PE, mV and ADC
  fPEtomVConversion = readConfig->GetDiscriminatorConversionFactormVperPE(iTelType);
  cout<<"The conversion factor from pe to mV is (Amplitude) [mV/pe]: "<<fPEtomVConversion<<endl;

  //Gain drop option
  bUseGainDrop = readConfig->GetGainDropUsage();
  if (bUseGainDrop) {
    cout << "Using gain drop" << endl;
    fPEtoDCconversion = readConfig->GetFADCDCtoPEconversionFactor(iTelType)/(1.+ 85.e-11 *readConfig->GetNSBRate(iTelType)*1.e3);
    cout<<"The conversion factor from pe to ADC is (Amplitude) [DC/pe]: "<<fPEtoDCconversion<<endl;
    fFADCconversion = fPEtoDCconversion/fPEtomVConversion/iScalingDivisor;
    cout<<"The conversion factor from mV to ADC per group (scaled with the 'scaling divisor') is "<<fFADCconversion<<endl;
  }

  else {
  fPEtoDCconversion = readConfig->GetFADCDCtoPEconversionFactor(iTelType);
  cout<<"The conversion factor from pe to ADC is (Amplitude) [DC/pe]: "<<fPEtoDCconversion<<endl;
  fFADCconversion = fPEtoDCconversion/fPEtomVConversion/iScalingDivisor;
  cout<<"The conversion factor from mV to ADC per group (scaled with the 'scaling divisor') is "<<fFADCconversion<<endl;
  }



  //the threshold when the low gain Channel is activated
  fHiLoGainThreshold = readConfig->GetFADCHiLoGainThreshold(iTelType);
  cout<<"High/Low gain threshold: "<<fHiLoGainThreshold<<endl;
  
  //the gain ratio of the low gain channel / high gain channel 
  fLowHiGainRatio = readConfig->GetFADCLowHiGainRatio(iTelType);
  cout<<"Low/High gain ratio: "<<fLowHiGainRatio<<endl;

  //The threshold of the discriminator
  fDiscThreshold = readConfig->GetDiscriminatorThreshold(iTelType);   
  cout<<"The Threshold of the discriminator (per cluster!) in ADC counts "<<fDiscThreshold<<endl;
  
  //The width of the output signal of the discriminator
  fWidthDiscriminator = readConfig->GetDiscriminatorOutputWidth(iTelType);
  cout<<"The width of the discriminator output signal in ns "<<fWidthDiscriminator<<endl;
  
  //The delay of the inverted signal in the CFD
  fDiscDelay = readConfig->GetDiscriminatorDelay(iTelType);
  cout<<"The delay of the inverted signal in the CFD in ns "<<fDiscDelay<<endl;
  
  //The attenuation of the non-inverted signal in the CFD
  fDiscConstantFractionAttenuation = readConfig->GetDiscriminatorAttenuation(iTelType);
  cout<<"The attenuation of the non-inverted signal in the CFD "<<fDiscConstantFractionAttenuation<<endl;
  
  //Do we use the RFB circuit? FOR NOW FORCED TO BE FALSE!!! Is it true?
  bDiscRFBUsage = readConfig->GetRFBUsage(iTelType);
  if ( bDiscRFBUsage )
    {
      bDiscRFBUsage = kFALSE;
      cout<<"The use of RFB is forced to be false for the camera snapshot logic: not (yet?) implemented!"<<endl;
    }
  else
    {
      cout<<"Do we use the RFB circuit: "<<bDiscRFBUsage<<endl;
    }
  
  //CDF part should not be used in this trigger logic, since the signal is digityzed snapshot by snapshot...
  bDiscCFDUsage = kFALSE;
  
  //The Constant Value in the RFB in the discriminator
  fDiscRFBConstant = readConfig->GetDiscriminatorRFBConstant(iTelType);
  cout<<"If the RFB circuit is used this is the magnitude of the rate dependend feedback  in mV/MHz "<<fDiscRFBConstant<<endl;
  
  //The dyamic value in the RFB in the discriminator
  fDiscRFBDynamic = readConfig->GetDiscriminatorRFBDynamic(iTelType);
  cout<<"The initial dynamic value in the RFB in mV (will be mulitplied with 0.18 in the sims). If RFB circuit is not used this is not used "<<fDiscRFBDynamic<<endl;
  
  SetDiscriminatorThresholdAndWidth(fDiscThreshold,fWidthDiscriminator);
  SetDiscriminatorDelayAndAttenuation(fDiscDelay, fDiscConstantFractionAttenuation);
  SetDiscriminatorRFBConstant(fDiscRFBConstant);
  SetDiscriminatorRFBDynamic(fDiscRFBDynamic);
  
}

void TriggerTelescopeCameraSnapshot::SetDiscriminatorThresholdAndWidth(Float_t threshold, Float_t width)
{ 

  if(threshold<=0 || width<=0)
    {
      cout<<"SetDiscriminatorThresholdAndWidth: the threshold and the width have to be > 0"<<endl;
      cout<<"You tried to set them to "<<threshold<<" and "<<width<<endl;
      exit(1);
    }

  fDiscThreshold = threshold; 
  fWidthDiscriminator =  width;

}

//-----------------------------------------------------------------------------------------
//Calculate a bias curve for the current trigger settings
//needs direct acces to tracegenerator to produce NSB by itself

void TriggerTelescopeCameraSnapshot::RunBiasCurve(UInt_t Trials,Float_t LowerBoundary,Float_t UpperBoundary,Float_t StepWidth,
						  TraceGenerator *tracegenerator,TelescopeData *TelData)
{
  cout<<"Doing a bias curve"<<endl;
  if ( fSamplingTime <= 0 )
    {
      cout<<"You need to set the sampling rate with SetSampleWidth(Float_t t)"<<endl;
      exit(0);
    }
  if ( StepWidth<0 )
    {
      cout<<"RunBiasCurve:  step size must be positive!. Something you should fix first"<<endl;
      exit(1);
    }
  if(UpperBoundary<0 || LowerBoundary<0 || UpperBoundary<LowerBoundary)
    {
      cout<<"RunBiasCurve:  Boundaries not compatible: they must be both positive and start smaller than stop! Something you should fix first"<<endl;
      exit(1);
    }

  //Setting the class pointer to the TelescopeData container
  telData = TelData;
  
  //Get the timing right. We do not work with showers here
  telData->fAveragePhotonArrivalTime = 0;
  
  Int_t NumScanPoints = Int_t((UpperBoundary-LowerBoundary)/StepWidth+1);
  Float_t TriggerWindow = fFADCSamplingWidth*1e-9;//fFADCSamplingWidth*iSamplingsWindow*1e-9;
  fBiasCurve.ResizeTo(NumScanPoints);
  fBiasCurve.Zero();
  fBiasCurveErr.ResizeTo(NumScanPoints);
  fBiasCurveErr.Zero();
  fBiasCurveScanPoints.ResizeTo(NumScanPoints);
  fBiasCurveScanPoints.Zero();
  fGroupRateVsThreshold.ResizeTo(NumScanPoints);
  fGroupRateVsThreshold.Zero();
  fGroupRateVsThresholdErr.ResizeTo(NumScanPoints);
  fGroupRateVsThresholdErr.Zero();
  
  cout<<"Going in loop"<<endl;
  
  for(UInt_t i=1;i<=Trials;i++)
    {
      telData->ResetTraces();
      //Load the NSB into the Traces and adjust the RFB Feedback
      //generate traces with trace generator
      tracegenerator->GenerateNSB();
      tracegenerator->BuildAllHighGainTraces();
      //Load event with trace from trace generator
      LoadEvent(telData);
      
      for(Int_t t = 0; t<NumScanPoints; t++)
	  {
	    if(bDebug)
	      cout<<"Lower boundary: "<<LowerBoundary<<" upper boundary: "<<UpperBoundary<<" Step: "<<StepWidth<<" being at "<<LowerBoundary+StepWidth*t<<endl;

	    SetDiscriminatorThresholdAndWidth(LowerBoundary+StepWidth*t,fWidthDiscriminator);
	    RunTrigger();

  	    if( telData->GetNumTriggeredGroups()!=0 )
	      {
	      fGroupRateVsThreshold[t]+=telData->GetNumTriggeredGroups();
	      if(telData->bTelescopeHasTriggered==kTRUE)
		{
		  fBiasCurve[t]++;
		}
	    }
	  else //will not find anything at higher thresholds, saves a lot of time
	    break;
	}
      
      //Output the state of the art every 1000 events
      if(i%1000 == 0 && i > 0 )
	{
	  cout<<"Events :"<<i<<"  simulated time "<<i*TriggerWindow<<" s"<<endl;
	  //cout<<"Dynamic Value of the RFB: "<<fDiscRFBDynamic<<" mV. The rate of zero crossings in MHz: "
	  //    << lZeroCrossings /(i* (fTraceLength-fDiscDelay)*1e-3*iNumSumPixGroups)<<endl;
	  for(Int_t t = 0; t<NumScanPoints; t++)
	    {
	      cout<<"NSB Trigger rate at "<<LowerBoundary+StepWidth*t<<" ADC threshold: triggers "
		  <<fBiasCurve[t]<<" rate:  "<<fBiasCurve[t]/(i*TriggerWindow)<<" Hz"<<endl;
	    }
	  cout<<endl;
	}
    }
  
  //Get the units right
  Float_t Norm = Trials*TriggerWindow;
  for(Int_t t = 0; t<NumScanPoints; t++)
    {
      fBiasCurveScanPoints[t]=LowerBoundary+StepWidth*t;
      fBiasCurveErr[t]=sqrt(1.0*fBiasCurve[t])/Norm;
      fBiasCurve[t]=fBiasCurve[t]/Norm;
      
      fGroupRateVsThresholdErr[t]=sqrt(fGroupRateVsThreshold[t])/iNumSumPixGroups/Norm;
      fGroupRateVsThreshold[t]=fGroupRateVsThreshold[t]/iNumSumPixGroups/Norm;
      cout<<"NSB Telescope Trigger rate at "<<LowerBoundary+StepWidth*t<<" ADC threshold "<<fBiasCurve[t]<<"+-"<<fBiasCurveErr[t]<<" Hz"<<endl;
      cout<<"Group rate :"<<fGroupRateVsThreshold[t]<<"+-"<<fGroupRateVsThresholdErr[t]<<" Hz"<<endl;
    }  
}
