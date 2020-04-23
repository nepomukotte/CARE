/* \file TriggerTelescopeNextNeighbor.cpp
   Implementation of a next neighbor trigger for AGIS
   1. Signal of X pixels are summed
   2. Discrimination of summed signal
   3. next neighbor logic (multiplicity of Y) 
*/

#include "TriggerTelescopeNextNeighbor.h"
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
TriggerTelescopeNextNeighbor::TriggerTelescopeNextNeighbor(ReadConfig *readConfig, int telType, TRandom3 *generator,Bool_t debug,Display *display)
{
  cout<<"You just bought an excellent next neighbor trigger"<<endl;
  bDebug = debug;
  debugDisplay = display;
  iTelType = telType;
  fTracesInSumGroups = NULL;
  fTracesInSumGroupsConstantFraction = NULL;	
  iClusterID = NULL;
  rand = generator;
  bUsePatches = kFALSE;
  iNumSumPixGroups= -1;
  fDiscThreshold = -1;
  fDiscRFBConstant = 0;  
  fDiscRFBDynamic = 0;
  bDiscRFBUsage= 0;
  fDiscDelay = -1; 
  fDiscConstantFractionAttenuation = -1;
  iMultiplicity = -1;
  fSamplingTime = -1;
  fSamplingTimeAveragePulse = -1; //sampling time of the average PE pulse
  fTraceLength = 100; //Length of simulated trace in ns
  fStartSamplingBeforeAverageTime = 30;

  lZeroCrossings = 0;               //holds the number of zerocrossings for one event counted over all pixels in the camera
  lNumEvents = 0;   

  SetParametersFromConfigFile(readConfig);	
}



void TriggerTelescopeNextNeighbor::CreateTraces()
{

  ///////////////////////////////////////////////////
  //
  //  Creating the arrays for the traces


  //delete trace arrays and create new ones
  if(fTracesInSumGroups)
    {
      //for(Int_t i=0;i<iNumSumPixGroups;i++)
	delete []  fTracesInSumGroups;
    }
  //create the array with all the vectors, one vector for each sumgroup
  fTracesInSumGroups = new vector<Float_t>[iNumSumPixGroups];
   

  if(fTracesInSumGroupsConstantFraction)
    {
      //for(Int_t i=0;i<iNumSumPixGroups;i++)
      delete []  fTracesInSumGroupsConstantFraction;
    }
  //create the array with all the vectors, one vector for each sumgroup
  fTracesInSumGroupsConstantFraction = new vector<Float_t>[iNumSumPixGroups];
  
  cout<<"Number of groups/pixel "<<iNumSumPixGroups<<endl;

}





TH1F TriggerTelescopeNextNeighbor::GetTraceHistogramThreshold(int GroupID ){

 
  TString title;
  title.Form("Trace in pixel %i",GroupID);
  TH1F hTrace ("hTrace",title,
	       telData->iNumSamplesPerTrace+1,telData->fAveragePhotonArrivalTime-fStartSamplingBeforeAverageTime-fSamplingTime*0.5,
	       telData->fAveragePhotonArrivalTime+fTraceLength-fStartSamplingBeforeAverageTime+fSamplingTime*0.5);
  hTrace.GetXaxis()->SetTitle("Time [ns]");
  hTrace.GetYaxis()->SetTitle("Amplitude [mV]");
	     
  for(Int_t i = 0; i<telData->iNumSamplesPerTrace; i++)
    {
      Float_t signal = fTracesInSumGroups[GroupID][i];
      hTrace.SetBinContent(i+1,signal);
    }

  return hTrace;

}


TH1F TriggerTelescopeNextNeighbor::GetTraceHistogramCFD(int GroupID){

 
  TString title;
  title.Form("Trace in pixel %i",GroupID);
  TH1F hTraceCFD("hTraceCFD",title,
		    telData->iNumSamplesPerTrace+1,telData->fAveragePhotonArrivalTime-fStartSamplingBeforeAverageTime-fSamplingTime*0.5,
					      telData->fAveragePhotonArrivalTime+fTraceLength-fStartSamplingBeforeAverageTime+fSamplingTime*0.5);
  hTraceCFD.GetXaxis()->SetTitle("Time [ns]");
  hTraceCFD.GetYaxis()->SetTitle("Amplitude [mV]");
	     
  hTraceCFD.SetLineColor(kRed);
  for(Int_t i = 0; i<telData->iNumSamplesPerTrace; i++)
    {
      hTraceCFD.SetBinContent(i+1,fTracesInSumGroupsConstantFraction[GroupID][i]);
    }

 return hTraceCFD;

}


void TriggerTelescopeNextNeighbor::ShowTrace(int GroupID){

 	TCanvas cTraceTriggeredGroup("cTraceTriggeredGroup","Trace in triggered group",700,500);
        cTraceTriggeredGroup.Divide(1,2);
	cTraceTriggeredGroup.cd(1);
	gPad->SetGrid();
      
	TH1F hTrace = GetTraceHistogramThreshold(GroupID);
	TH1F hTraceCFD = GetTraceHistogramCFD(GroupID);

	
	hTrace.SetMaximum(hTraceCFD.GetMaximum()>hTrace.GetMaximum() ? 
                                 hTraceCFD.GetMaximum()+10:hTrace.GetMaximum()+10 );
	hTrace.SetMinimum(hTraceCFD.GetMinimum()<hTrace.GetMinimum() ? 
                                 hTraceCFD.GetMinimum()-10:hTrace.GetMinimum()-10 );

	hTrace.Draw();
	hTraceCFD.Draw("same");
    
    TLine lDisc(telData->fAveragePhotonArrivalTime-fStartSamplingBeforeAverageTime-fSamplingTime*0.5,fDiscThreshold,
	         telData->fAveragePhotonArrivalTime+fTraceLength-fStartSamplingBeforeAverageTime+fSamplingTime*0.5, fDiscThreshold);
	lDisc.Draw();

	TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
	
	//
	// While reading the input process gui events asynchronously
	//
	timer.TurnOn();
	TString input = Getline("Type <return> to go on: ");
	timer.TurnOff();  

}



//-----------------------------------------------------------------------------------------------------------------------
//Loads the traces from the pixels (sums them up)
void  TriggerTelescopeNextNeighbor::LoadEvent(TelescopeData *TelData)

//vector<Float_t>* vNSBOnly, vector<Float_t>* vNSBandCherenkov, Float_t AverageArrivalTime)
{

 telData = TelData;

 telData->bTriggeredGroups.assign(iNumSumPixGroups,kFALSE);

 telData->fDiscriminatorTime.assign(iNumSumPixGroups,-1e6);

 //can be removed is only needed to check the Time over threshold of all pixels in a trigger.
 telData->fTimeOverThreshold.assign(iNumSumPixGroups,0);

 //The start sample for the delayed trace;
 Int_t iStartSample = (int)(fDiscDelay/fSamplingTime)+1;

 if(bDebug)
     {
        cout<<"Start sample for delayed pulse "<<iStartSample<<endl;
     }

 Float_t fOffsetDueToRFBinCFD = 0.0;
   if(bDiscRFBUsage) fOffsetDueToRFBinCFD = -0.18*fDiscRFBDynamic; 

 for(Int_t g=0;g<iNumSumPixGroups;g++)
    {
      fTracesInSumGroups[g].assign(telData->iNumSamplesPerTrace,0.0);
      fTracesInSumGroupsConstantFraction[g].assign(telData->iNumSamplesPerTrace,fOffsetDueToRFBinCFD);
      //loop over all group members and add their trace to the trace of the sumgroup
      for(UInt_t n = 0; n<iSumGroupMembers[g].size();n++)
         {      
           Int_t memberID = iSumGroupMembers[g][n];
           for(Int_t i = 0 ; i<telData->iNumSamplesPerTrace; i++)
	          {
	             Float_t fsignal = telData->fTraceInPixel[memberID][i]*fPEtomVConversion ;
                 //cout<<fsignal<<"  "<<telData->fTraceInPixel[memberID][i]<<"  "<<fPEtomVConversion<<endl;
                 fsignal = fsignal < fClippingLevel && bDoClipping == kTRUE   ? fClippingLevel : fsignal;	      
	             Float_t fsigDelayed = i<iStartSample ? telData->fTraceInPixel[memberID][0]*fPEtomVConversion : telData->fTraceInPixel[memberID][i-iStartSample]*fPEtomVConversion;
                 fsigDelayed = fsigDelayed > fClippingLevel && bDoClipping == kTRUE  ? fClippingLevel : fsigDelayed;	      
	      
	             Float_t CFDsignal = fsignal*fDiscConstantFractionAttenuation-fsigDelayed;
	             fTracesInSumGroups[g][i]+=fsignal;
	             fTracesInSumGroupsConstantFraction[g][i]+=CFDsignal;
	         }
            if(bDebug)
                {
                     cout<<"Summed group "<<g<<" added pixel "<<memberID<<endl; 
	             //    ShowTrace(g, false);
               }
           
	     }
    }
  

  //cout<<"Finished loading the event"<<endl;

}


//-----------------------------------------------------------------------------------------------------------------------
//This does the actual triggering
//It returns true if the telescope has triggered and false if not
Bool_t  TriggerTelescopeNextNeighbor::RunTrigger()
{

    //Get all the zerocrossings for the RFB
    lZeroCrossings += GetNumZeroCrossings();
    lNumEvents++;

    //    cout<<"Have done the zerocrossings "<<lZeroCrossings<<endl;

    //Set the updated RFB feedback
    if(bDiscRFBUsage)
      {

	   if(bDebug)
	     cout<<lNumEvents<<" RFB dynamic value  "<<fDiscRFBConstant * lZeroCrossings /(lNumEvents* (fTraceLength-fDiscDelay)*1e-3*iNumSumPixGroups)<<endl;

       SetDiscriminatorRFBDynamic(fDiscRFBConstant * lZeroCrossings /(lNumEvents* (fTraceLength-fDiscDelay)*1e-3*iNumSumPixGroups));
      }

   //Run the CFD simulation
   telData->iNumTriggeredGroups=0;
   for(Int_t i=0;i<iNumSumPixGroups;i++)
    {
      if(RunDiscriminator(i))
	     telData->iNumTriggeredGroups++;
    }

  //  cout<<"Done with the CFD"<<endl<<endl;

  telData->vTriggerCluster.clear();
  //Below is the L2 simulation.
   if(bUsePatches==kFALSE)
     {
      if(bDebug)
       cout<<"Not using trigger patches in the simulation"<<endl;

      iPixelTriggeredInPatch.resize(1);
      telData->bTelescopeHasTriggered =  RunL2Patch(0,&(telData->fTelescopeTriggerTime));//important, first argument has to be 0
      telData->vTriggerCluster = iPixelTriggeredInPatch[0];
     }   
   else //If we use patches
     {
      if(bDebug)
       cout<<"Using trigger patches in the simulation"<<endl;
	   
      iPixelTriggeredInPatch.resize(vPatch.size());
      telData->bTelescopeHasTriggered = RunL2WithPatches();
     }

   
   if(telData->bTelescopeHasTriggered && bDebug)
     {
       cout<<"Event triggered telescope"<<endl; 
       cout<<"at time "<< telData->fTelescopeTriggerTime<<endl;
     }
  
  return telData->bTelescopeHasTriggered;
}


//------------
//
// An L2 trigger that does care about patches
Bool_t  TriggerTelescopeNextNeighbor::RunL2WithPatches()
{
  // cout<<"Patches"<<endl;
  //Go over all patches and find those which would fullfil the trigger conditions
  vector< Bool_t > bPatchTrigger(vPatch.size(),kFALSE);
  Float_t *fPatchTriggertimes = new Float_t[vPatch.size()];

  for(UInt_t p=0;p<vPatch.size();p++)
    {
      bPatchTrigger[p] = RunL2Patch(p,fPatchTriggertimes);
    }
  
  //Find earliest time a patch triggered
  Int_t* indext = new Int_t[ vPatch.size()  ];
  TMath::Sort((Int_t)vPatch.size(),fPatchTriggertimes,indext,kFALSE);
  telData->fTelescopeTriggerTime = fPatchTriggertimes[ indext[0] ];

  delete [] fPatchTriggertimes;
                         
  //cout<<"Done with patches"<<endl;

  Bool_t triggered = bPatchTrigger[ indext[0] ];

  if(triggered == kTRUE)
   {
      telData->vTriggerCluster = iPixelTriggeredInPatch[ indext[0] ];
      if(bDebug)
       {
         cout<<"we triggered: "<<triggered<<"; Patch that triggered: "<<indext[0]<<"  out of a total of "<<vPatch.size()<<" Patches"<<endl;
	     cout<<"Trigger time: "<<telData->fTelescopeTriggerTime<<"; Patch that triggered: "<<indext[0]<<"  out of a total of "<<vPatch.size()<<" Patches"<<endl;
       }
   }

  delete [] indext;

  return triggered;

}

//----------------------------------------------------------------------------------------
//
// An L2 trigger that does care about patches
Bool_t  TriggerTelescopeNextNeighbor::RunL2Patch(Int_t PatchNumber,Float_t *fPatchTriggerTimes)
{
  if(bDebug)
   {
     cout<<"Running L2 Patch "<<PatchNumber<<endl;
   }

  if(iClusterID)
    delete [] iClusterID;
  iClusterID = new Int_t[iNumSumPixGroups];

  Int_t *iNumGroupsInCluster = new Int_t[iNumSumPixGroups]; //How many groups are in each cluster of this patch

  for(Int_t i=0;i<iNumSumPixGroups;i++)
    {
      iClusterID[i] = -1;
      iNumGroupsInCluster[i]=0; 
    }

  vector< int > i_pix;
  vGroupsInCluster = new vector< vector<int> >;
  vGroupsInCluster->assign(iNumSumPixGroups , i_pix );

  vector< float > f_pix;
  vTriggerTimesInCluster = new vector< vector<float> >;
  vTriggerTimesInCluster->assign(iNumSumPixGroups , f_pix );

  //cout<<"allocated arrays and vectors"<<endl;
  //cout<<"Number of summed groups "<<iNumSumPixGroups<<endl;

  Int_t NGroups = bUsePatches==kTRUE ? vPatch[PatchNumber].size() :  iNumSumPixGroups;
  if(bDebug)
  cout<<"Number of groups in this patch: "<<NGroups<<endl;
  for(Int_t i=0;i<NGroups;i++)
    {
      Int_t GroupID = bUsePatches==kTRUE ? vPatch[PatchNumber][i] : i;
    
      if (telData->bTriggeredGroups[GroupID]==kTRUE) //(summed) pixel has triggered
	{
	  iNumGroupsInCluster[GroupID] = CalcCluster(GroupID,GroupID,PatchNumber);
	  // cout<<iNumGroupsInCluster[i]<<endl;       
	}
	  
    }

  //cout<<"Finished finding all clusters"<<endl<<endl;
   
  //Sort all cluster descending
  Int_t* index = new Int_t[iNumSumPixGroups];
  TMath::Sort(iNumSumPixGroups,iNumGroupsInCluster,index);
 
/*  for(int n=0;n<iNumSumPixGroups;n++)
	 cout<<iNumGroupsInCluster[index[n]]<<"  ";
	 
	 cout<<endl;
  */
  //Find the cluster that triggered first
  Float_t fPatchTriggerTime = 1e6;
  Int_t t = 0;
  Int_t tmin = 0;
  while(iNumGroupsInCluster[index[t]]>=iMultiplicity )
    {
      Int_t igroups = iNumGroupsInCluster[index[t]];

      if(bDebug)
      cout<<"cluster "<<index[t]<<" has "<<igroups<<" groups"<<endl;

      //Find the L2 trigger time of this cluster
      Float_t *times = new Float_t[ igroups ];
      for(Int_t i=0;i<igroups;i++)
	  {
           times[i] = (vTriggerTimesInCluster->at(index[t]))[i];
           if(bDebug)
           cout<<"triggered group in cluster "<<vGroupsInCluster->at(index[t])[i]<<" t: "<<times[i]<<endl;
	  }
      // cout<<endl;
      //sort the times in increasing order
      Int_t* indext = new Int_t[ igroups  ];
      TMath::Sort(igroups,times,indext,kFALSE);
      //Trigger time is determined by the /pixel summed group that triggered last (iMultiplicity-1) 
      Float_t triggertime = times[ indext[iMultiplicity-1] ];
      if(triggertime<fPatchTriggerTime)
	    {
	       fPatchTriggerTime = triggertime;
	       tmin=t;
	    }
       t++;

       delete [] times;
       delete [] indext;

       if(t==iNumSumPixGroups)
          break;

    }
  fPatchTriggerTime+= telData->fAveragePhotonArrivalTime;
  
  fPatchTriggerTimes[PatchNumber] = fPatchTriggerTime;

//move the vector with the groups of pixels that are in the trigger cluster to a separate vector
  //which can be retrieved with GetTriggerCluster()
  //cout<< vGroupsInCluster->at(index[tmin]).size()<<endl;
  //cout<<(telData->vTriggerCluster).size()<<endl;
  iPixelTriggeredInPatch[PatchNumber]=vGroupsInCluster->at(index[tmin]);
  //cout<<"Have transferred the TriggerCluster"<<telData->vTriggerCluster.size()<<endl;
  //Now trigger if we have enough groups in the largest cluster
 

  Bool_t bPatchTrigger = kFALSE;

  if(iNumGroupsInCluster[index[0]]>=iMultiplicity)
    {
      //cout<<iNumGroupsInCluster[index[0]]<<endl;
      bPatchTrigger= kTRUE;
    }

  delete [] iNumGroupsInCluster;
  delete vGroupsInCluster;
  delete vTriggerTimesInCluster;
  delete [] index;

  return bPatchTrigger;
}



//---------------------------------------------------------------------------------------
//Get the numbers of zerocrossings for the loaded event from negative to positive
Long_t TriggerTelescopeNextNeighbor::GetNumZeroCrossings()
{

  Long_t lCrossings=0;

  Int_t iStartSample = (int)(fDiscDelay/fSamplingTime)+1;

   for(Int_t g=0;g<iNumSumPixGroups;g++)
    {
      
      Bool_t bAboveZero = kFALSE;
      for(Int_t i = iStartSample; i<telData->iNumSamplesPerTrace; i++)
	{
	  //cout<<g<<"  "<<i<<endl;
	  Bool_t bCFDOut = fTracesInSumGroupsConstantFraction[g][i]>=0 ? kTRUE : kFALSE;

	  if(bCFDOut && !bAboveZero )
	    {
	      lCrossings++;
	      bAboveZero=kTRUE;
	    }
	  else if(!bCFDOut && bAboveZero)
	    bAboveZero=kFALSE;
	}

    }
   

   return lCrossings;

}

//-----------------------------------------------------------------------------------------------------------------------
// Does the constant fraction discrimination for a sum group using a VERITAS Discriminator Design
//This is a two level discriminator
//The first one is a threshold discriminator that stays up as long as the signal is above threshold
//The second one gives an output if the sum of the attenuated (0.4) signal plus offset and the delayed and inverted copy of the
//input signal cross zero 
//The whole thing gives an output if both discriminators give logic one
Bool_t TriggerTelescopeNextNeighbor::RunDiscriminator(Int_t GroupID)
{

 
  Int_t iStartSample = (int)(fDiscDelay/fSamplingTime)+1;

  //1. Threshold discriminator: output as long as signal is above threshold
  //2. Output as long if  CFD trace attenuated copy plus inverted copy plus constant offset output is negativ
  //3. Trigger if 1. and 2. are high
  telData->bTriggeredGroups[GroupID] = kFALSE;
  telData->fDiscriminatorTime[GroupID] = -1e6;

 

  for(Int_t i = iStartSample; i<telData->iNumSamplesPerTrace; i++)
    {
	  Float_t fsignal = fTracesInSumGroups[GroupID][i];

	  Bool_t bCFDOut = kTRUE; //always true if we do not use the CFD part
      if(bDiscCFDUsage) bCFDOut=fTracesInSumGroupsConstantFraction[GroupID][i]>=0 ? kTRUE : kFALSE;
	
	  //In case both discriminators are firing we have a trigger
	  if( bCFDOut && fsignal<=fDiscThreshold &&  telData->bTriggeredGroups[GroupID] == kFALSE ) 
        {
	      telData->fDiscriminatorTime[GroupID] = i*fSamplingTime-fStartSamplingBeforeAverageTime;
	      telData->bTriggeredGroups[GroupID] = kTRUE;

		  //loop over the rest of the trace to find out how long the discriminator is above threshold
		  //this can be deactivated because it is only used in simulation studies done for the 1m sst.
		  Int_t iToT=0;
		  while(fsignal<=fDiscThreshold && i<telData->iNumSamplesPerTrace-1)
			{
			  iToT++;

			  i++;
			  fsignal = fTracesInSumGroups[GroupID][i];
			}
          telData->fTimeOverThreshold[GroupID] = iToT*fSamplingTime;
		  break; //exit loop because we triggered the trace and we are only interested in the first trigger, saves some compute time
	    }
    }


     if(bDebug && telData->bTriggeredGroups[GroupID] == kTRUE)
      {
        cout<<"Trace on which the discriminator was run and triggered; groupid "<<GroupID<<endl;
        debugDisplay->AddDiscriminatorTraces(telData->GetTelescopeID(),GroupID,fDiscThreshold,GetTraceHistogramThreshold(GroupID),GetTraceHistogramCFD(GroupID));

      }

    return telData->bTriggeredGroups[GroupID];
}


// --------------------------------------------------------------------------------------------
//
//Helper function to loop over all neighbours to find triggered pixels belonging to the cluster
//returns the number of pixel in that cluster

Int_t TriggerTelescopeNextNeighbor::CalcCluster(Int_t GroupID, Int_t ClusterID, Int_t PatchID)
{

  if(bDebug)
	  cout<<"visiting group "<<GroupID<<endl;

  // If we have visited this group in this round ... do nothing.
  if (iClusterID[GroupID]==ClusterID)
    return 0;
  // Assign the new cluster ID this Group
  iClusterID[GroupID]=ClusterID;
  //Put the group in the list of groups that belong to this cluster centered around the group with ID ClusterID
  vGroupsInCluster->at(ClusterID).push_back(GroupID);
  
  
  //Add the trigger time of the group to the vector of this cluster
  vTriggerTimesInCluster->at(ClusterID).push_back(telData->fDiscriminatorTime[GroupID]);
  if(bDebug)
  {
	 cout<<"Pixel number: "<<vTriggerTimesInCluster->at(ClusterID).size()<<endl;
	 cout<<"added group "<<GroupID<<", which triggered at time "<<telData->fDiscriminatorTime[GroupID]<<endl;
  }

  // Need this to store the number of groups in the cluster
  Int_t NumGroupsInCluster = 1;

  // Now do the same with all its neighbors
  for(UInt_t n = 0; n<iSumGroupNeighbors[GroupID].size();n++)
    {
      Int_t GroupIDNeighbor = iSumGroupNeighbors[GroupID][n];
      if (telData->bTriggeredGroups[GroupIDNeighbor] 
           && fabs(telData->fDiscriminatorTime[ClusterID]-telData->fDiscriminatorTime[GroupIDNeighbor])<fWidthDiscriminator 
           && GroupInPatch(GroupIDNeighbor,PatchID))
	    NumGroupsInCluster += CalcCluster(GroupIDNeighbor,ClusterID,PatchID);
    }
 // return the number of groups in this cluster
  return NumGroupsInCluster;

}

//--------------------------------------------------------------------------------------------
//
// returns true if group is in patch
// if patches are not used it always returns true
Bool_t TriggerTelescopeNextNeighbor::GroupInPatch(Int_t GroupID,Int_t PatchID)
{

  if(!bUsePatches)
    return kTRUE;
 
  for(UInt_t i=0;i<vPatch[PatchID].size();i++)
    {
      if(GroupID==vPatch[PatchID][i])
	 return kTRUE;
    }

  return kFALSE;
}


//-----------------------------------------------------------------------------------------
//Calculate a bias curve for the current trigger settings
//needs direct acces to tracegenerator to produce NSB by itself
void TriggerTelescopeNextNeighbor::RunBiasCurve(UInt_t Trials,Float_t LowerBoundary,Float_t UpperBoundary,Float_t StepWidth, TraceGenerator *tracegenerator,TelescopeData *TelData)
{

  cout<<"Doing a bias curve"<<endl;

  if(fSamplingTime <= 0)
    {
      cout<<"You need to set the sampling rate with SetSampleWidth(Float_t t)"<<endl;
      exit(0);
    }

  if(UpperBoundary<LowerBoundary && StepWidth>0)
    {
      cout<<"RunBiasCurve:  Upper boundary is less than the Lower Boundary and step size is larger 0. Something you should fix first"<<endl;
      exit(1);
    }

  if(UpperBoundary>LowerBoundary && StepWidth<0)
    {
      cout<<"RunBiasCurve:  Upper boundary is larger than the Lower Boundary and step size is < 0. Something you should fix first"<<endl;
      exit(1);
    }

  //Setting the class pointer to the TelescopeData container
  telData = TelData;

  //Get the timing right. We do not work with showers here
  telData->fAveragePhotonArrivalTime = 0;

  Int_t NumScanPoints = Int_t((UpperBoundary-LowerBoundary)/StepWidth+1);
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
      //cout<<"Number of groups that triggered "<<telData->GetNumTriggeredGroups()<<endl;
  	  if( telData->GetNumTriggeredGroups()!=0 )
	    {
	      fGroupRateVsThreshold[t]+=telData->GetNumTriggeredGroups();
	      if(telData->bTelescopeHasTriggered==kTRUE)
	          fBiasCurve[t]++;
	    }
	  else //will not find anything at higher thresholds, saves a lot of time
	    break;


	}

      //Output the state of the art every 1000 events
      if(i%1000 == 0 && i > 0 )
	{
	  cout<<"Events :"<<i<<"  simulated time "<<i*fTraceLength*1e-9<<" s"<<endl;
	  cout<<"Dynamic Value of the RFB: "<<fDiscRFBDynamic<<" mV. The rate of zero crossings in MHz: "<< lZeroCrossings /(i* (fTraceLength-fDiscDelay)*1e-3*iNumSumPixGroups)<<endl;
	  for(Int_t t = 0; t<NumScanPoints; t++)
	    {
	      cout<<"NSB Trigger rate at "<<LowerBoundary+StepWidth*t<<" mV threshold: triggers "
                   <<fBiasCurve[t]<<" rate:  "<<fBiasCurve[t]/(i*fTraceLength*1e-9)<<" Hz"<<endl;
	    }
	  cout<<endl;
	}

    }

  //Get the units right
  for(Int_t t = 0; t<NumScanPoints; t++)
    {
      fBiasCurveScanPoints[t]=LowerBoundary+StepWidth*t;
      fBiasCurveErr[t]=sqrt(1.0*fBiasCurve[t])/(Trials*fTraceLength*1e-9);
      fBiasCurve[t]=fBiasCurve[t]/(Trials*fTraceLength*1e-9);

      fGroupRateVsThresholdErr[t]=sqrt(fGroupRateVsThreshold[t])/iNumSumPixGroups/(Trials*fTraceLength*1e-9);
      fGroupRateVsThreshold[t]=fGroupRateVsThreshold[t]/iNumSumPixGroups/(Trials*fTraceLength*1e-9);
      cout<<"NSB Telescope Trigger rate at "<<LowerBoundary+StepWidth*t<<" mV threshold "<<fBiasCurve[t]<<"+-"<<fBiasCurveErr[t]<<" Hz"<<endl;
      cout<<"Group rate :"<<fGroupRateVsThreshold[t]<<"+-"<<fGroupRateVsThresholdErr[t]<<" Hz"<<endl;
    }
   


}

//----------------------------------------------------------------------------------------------------------
//
// Set the threshold mV of the discriminator and the width of the output pulse (ns) 
void TriggerTelescopeNextNeighbor::SetDiscriminatorThresholdAndWidth(Float_t threshold, Float_t width)
{ 

  if(threshold>=0 || width<=0)
    {
      cout<<"SetDiscriminatorThresholdAndWidth: the threshold has to be < 0 and the width have to be > 0"<<endl;
      cout<<"You tried to set them to "<<threshold<<" and "<<width<<endl;
      exit(1);
    }

  fDiscThreshold = threshold; 
  fWidthDiscriminator =  width;

}

//------------------------------------------------------------------------------------------------------
//
// Set the rate feedback of the discriminator
// Putting it to 0 will turn it off 
void TriggerTelescopeNextNeighbor::SetDiscriminatorRFBConstant(Float_t rfb)
{ 

  if(rfb<0)
    {
      cout<<"SetDiscriminatorRFBConstant: The RFB is set to  a value < 0 put it to either 0 or a lager value"<<endl;
      cout<<"You tried to set them to "<<rfb<<endl;
      exit(1);
    }

  fDiscRFBConstant = rfb; 

}

//------------------------------------------------------------------------------------------------------
//
// Set the dynamic value of the rate feedback of the discriminator
// Putting it to 0 will turn it off 
void TriggerTelescopeNextNeighbor::SetDiscriminatorRFBDynamic(Float_t rfb)
{ 

  if(rfb<0)
    {
      cout<<"SetDiscriminatorRFBDynamic: The dynamic value for the RFB is set to  a value < 0 put it to either 0 or a lager value"<<endl;
      cout<<"You tried to set them to "<<rfb<<endl;
      exit(1);
    }

  fDiscRFBDynamic = rfb; 

}

//------------------------------------------------------------------------------------------------------
//
// Set the dynamic value of the rate feedback of the discriminator
// Putting it to 0 will turn it off 
void TriggerTelescopeNextNeighbor::SetDiscriminatorRFBUsage(Bool_t rfbuse)
{ 

  bDiscRFBUsage = rfbuse; 

}


//------------------------------------------------------------------------------------------------------
//
// Set the delay of the inverted signal in the discriminator and the attenuation of the non-inverted signal
// 
void TriggerTelescopeNextNeighbor::SetDiscriminatorDelayAndAttenuation(Float_t delay, Float_t attenuation)
{ 

  if(delay<0 || attenuation <0)
    {
      cout<<"SetDiscriminatorDelayAndConstantFraction: The Delay or the Attenuation of the discriminator have been set to a value smaller 0"<<endl;
      cout<<"You tried to set them to "<<delay<<"  "<<attenuation<<endl;
      exit(1);
    }

  fDiscDelay = delay;
  fDiscConstantFractionAttenuation = attenuation;

}


//Reads in  the config file and sets all variables
void   TriggerTelescopeNextNeighbor::SetParametersFromConfigFile(ReadConfig *readConfig ){

   cout <<endl<< "TriggerTelescopeNextNeighbor::SetParametersFromConfigFile " << endl;

   //Trace length in ns
   fTraceLength = readConfig->GetTraceLength(iTelType);
   cout<<"The analog trace is "<<fTraceLength<<" ns long"<<endl;

   //Offset from average photon arrival time
   fStartSamplingBeforeAverageTime = readConfig->GetStartSamplingTimeOffsetFromAveragePhotonTime(iTelType);
   cout<<"The analog trace starts to be sampled "<<fStartSamplingBeforeAverageTime<<" ns before the average photon arrival time"<<endl; 

   //Width of one sample  of the trace
   fSamplingTime = readConfig->GetSamplingTime(iTelType);
   cout<<"The sampling time of the analog trace is "<<fSamplingTime<<" ns"<<endl;
 
   //The threshold of the discriminator
   fDiscThreshold = readConfig->GetDiscriminatorThreshold(iTelType);   
   cout<<"The Threshold of the discriminator in mV "<<fDiscThreshold<<endl;
   
   //The width of the output signal of the discriminator
   fWidthDiscriminator = readConfig->GetDiscriminatorOutputWidth(iTelType);
   cout<<"The width of the discriminator output signal in ns "<<fWidthDiscriminator<<endl;

   //The delay of the inverted signal in the CFD
   fDiscDelay = readConfig->GetDiscriminatorDelay(iTelType);
   cout<<"The delay of the inverted signal in the CFD in ns "<<fDiscDelay<<endl;

   //The attenuation of the non-inverted signal in the CFD
   fDiscConstantFractionAttenuation = readConfig->GetDiscriminatorAttenuation(iTelType);
   cout<<"The attenuation of the non-inverted signal in the CFD "<<fDiscConstantFractionAttenuation<<endl;
    
   //Do we use the RFB circuit?
   bDiscRFBUsage = readConfig->GetRFBUsage(iTelType);
   cout<<"Do we use the RFB circuit: "<<bDiscRFBUsage<<endl;

   //Do we use the CFD part of the discriminator?
   bDiscCFDUsage = readConfig->GetCFDUsage(iTelType);
   cout<<"Do we use the CFD circuit: "<<bDiscCFDUsage<<endl;

   //The Constant Value in the RFB in the discriminator
   fDiscRFBConstant = readConfig->GetDiscriminatorRFBConstant(iTelType);
   cout<<"If the RFB circuit is used this is the magnitude of the rate dependend feedback  in mV/MHz "<<fDiscRFBConstant<<endl;
       
   //The dyamic value in the RFB in the discriminator
   fDiscRFBDynamic = readConfig->GetDiscriminatorRFBDynamic(iTelType);
   cout<<"The initial dynamic value in the RFB in mV (will be mulitplied with 0.18 in the sims). If RFB circuit is not used this is not used "<<fDiscRFBDynamic<<endl;
     
   //How many goups need to be in a cluster for a telescope trigger
   iMultiplicity = readConfig->GetGroupMultiplicity(iTelType);
   cout<<"The multiplicity requirement for a telescope trigger is "<<iMultiplicity<<endl;

   //The conversion factor at the input of the discriminator in units of  mV / pe peak ampltude
   fPEtomVConversion = readConfig->GetDiscriminatorConversionFactormVperPE(iTelType);
   cout<<"The conversion factor from pe to mV is (Amplitude) [mV/pe]: "<<fPEtomVConversion<<endl;

   //Do we want to use the Patches 
   bUsePatches = readConfig->GetTriggerPatchUsage(iTelType);
   cout<<"Will use patches : "<<bUsePatches<<endl;
 
   //Patches
   vPatch = readConfig->GetTriggerPatches(iTelType);

   //Do we clip?
   bDoClipping = readConfig->GetClippingUsage(iTelType);
   cout<<"Will use clipping: "<<bDoClipping<<endl;

   //the clipping level in mV
   fClippingLevel = readConfig->GetClippingLevel(iTelType);
   cout<<"the clipping level in mV is: "<<fClippingLevel<<endl;

   if(readConfig->GetUseSumTrigger(iTelType))
     {
         cout<<"Will simulate a sumtrigger"<<endl;
        iNumSumPixGroups = readConfig->GetNumberGroups(iTelType);
        iSumGroupNeighbors = readConfig->GetNeighborsOfGroup(iTelType);
        iSumGroupMembers = readConfig->GetMembersOfGroups(iTelType);
        if(iNumSumPixGroups <=0)
          {
             cout<<"The sumtrigger simulation was turned on but the number of summed groups has been put to "<<iNumSumPixGroups<<endl;
             cout<<"do something about it!!"<<endl;
             exit(1);
          }
     }
   else //If no sumtrigger is used
     {
        cout<<"This simulation does not use a sumtrigger"<<endl;
		iNumSumPixGroups= readConfig->GetNumberPixels(iTelType);
		iSumGroupNeighbors = readConfig->GetNeighbors(iTelType);
		vector<int> vmembers;
        for(int i = 0; i<iNumSumPixGroups; i++)
         {
			iSumGroupMembers.push_back(vmembers);
            iSumGroupMembers[i].assign(1,i);
         }
      }
      
   SetDiscriminatorThresholdAndWidth(fDiscThreshold,fWidthDiscriminator);
   SetDiscriminatorDelayAndAttenuation(fDiscDelay, fDiscConstantFractionAttenuation);
   SetDiscriminatorRFBConstant(fDiscRFBConstant);
   SetDiscriminatorRFBDynamic(fDiscRFBDynamic);

  //Create the traces for the summed pixels
  CreateTraces();
}

