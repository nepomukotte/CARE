/* \file TriggerTelescopeBase.cpp
   Implementation of a Telescope Trigger Baseclass
   The base class provides 
    1. an implementation of a discriminator (leading edge and CFD)
    2. Summation of the signals of all pixels in a "group" to form one trigger pixel
    3. Evaluation if X neighboring trigger pixels are above the trigger threshold within a coincidence window
    4. Definition of "patches", Groups can belong to one or more patches. The trigger condition can only be met if all pixels belong to the same patch.
    5. Execution of bias curves

    Derived classes overload all the functions that require adapting to meet a specific trigger logic/topology
     
*/

#include "TriggerTelescopeBase.h"
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
TriggerTelescopeBase::TriggerTelescopeBase(ReadConfig *readConfig, int telType, TRandom3 *generator,Bool_t debug,Display *display)
{
  cout<<"Initiating Base Class of the Telescope Trigger"<<endl;
  bDebug = debug;
  debugDisplay = display;
  iTelType = telType;
  fTracesInTriggerPixels = NULL;
  fTracesInTriggerPixelsConstantFraction = NULL;	
  iClusterID = NULL;
  rand = generator;
  bUsePatches = kFALSE;
  iNumTriggerPixels= -1;
  fDiscThreshold = -1;
  fDiscDelay = -1; 
  fDiscConstantFractionAttenuation = -1;
  iMultiplicity = -1;
  fSamplingTime = -1;
  fSamplingTimeAveragePulse = -1; //sampling time of the average PE pulse
  fTraceLength = 100; //Length of simulated trace in ns
  fStartSamplingBeforeAverageTime = 30;


  SetParametersFromConfigFile(readConfig);	
}



void TriggerTelescopeBase::CreateTraces()
{

  ///////////////////////////////////////////////////
  //
  //  Creating the arrays for the traces


  //delete trace arrays and create new ones
  if(fTracesInTriggerPixels)
    {
      //for(Int_t i=0;i<iNumTriggerPixels;i++)
	delete []  fTracesInTriggerPixels;
    }
  //create the array with all the vectors, one vector for each sumgroup
  fTracesInTriggerPixels = new vector<Float_t>[iNumTriggerPixels];
   

  if(fTracesInTriggerPixelsConstantFraction)
    {
      //for(Int_t i=0;i<iNumTriggerPixels;i++)
      delete []  fTracesInTriggerPixelsConstantFraction;
    }
  //create the array with all the vectors, one vector for each sumgroup
  fTracesInTriggerPixelsConstantFraction = new vector<Float_t>[iNumTriggerPixels];
  
  cout<<"Number of trigger pixels "<<iNumTriggerPixels<<endl;

}





TH1F TriggerTelescopeBase::GetTraceHistogramThreshold(int TrigPixID ){

 
  TString title;
  title.Form("Trace in pixel %i",TrigPixID);
  TH1F hTrace ("hTrace",title,
	       telData->iNumSamplesPerTrace+1,telData->fAveragePhotonArrivalTime-fStartSamplingBeforeAverageTime-fSamplingTime*0.5,
	       telData->fAveragePhotonArrivalTime+fTraceLength-fStartSamplingBeforeAverageTime+fSamplingTime*0.5);
  hTrace.GetXaxis()->SetTitle("Time [ns]");
  hTrace.GetYaxis()->SetTitle("Amplitude [mV]");
	     
  for(Int_t i = 0; i<telData->iNumSamplesPerTrace; i++)
    {
      Float_t signal = fTracesInTriggerPixels[TrigPixID][i];
      hTrace.SetBinContent(i+1,signal);
    }

  return hTrace;

}


TH1F TriggerTelescopeBase::GetTraceHistogramCFD(int TrigPixID){

 
  TString title;
  title.Form("Trace in pixel %i",TrigPixID);
  TH1F hTraceCFD("hTraceCFD",title,
		    telData->iNumSamplesPerTrace+1,telData->fAveragePhotonArrivalTime-fStartSamplingBeforeAverageTime-fSamplingTime*0.5,
					      telData->fAveragePhotonArrivalTime+fTraceLength-fStartSamplingBeforeAverageTime+fSamplingTime*0.5);
  hTraceCFD.GetXaxis()->SetTitle("Time [ns]");
  hTraceCFD.GetYaxis()->SetTitle("Amplitude [mV]");
	     
  hTraceCFD.SetLineColor(kRed);
  for(Int_t i = 0; i<telData->iNumSamplesPerTrace; i++)
    {
      hTraceCFD.SetBinContent(i+1,fTracesInTriggerPixelsConstantFraction[TrigPixID][i]);
    }

 return hTraceCFD;

}


void TriggerTelescopeBase::ShowTrace(int TrigPixID){

 	TCanvas cTraceTriggeredGroup("cTraceTriggeredGroup","Trace in triggered triggerpixel",700,500);
        cTraceTriggeredGroup.Divide(1,2);
	cTraceTriggeredGroup.cd(1);
	gPad->SetGrid();
      
	TH1F hTrace = GetTraceHistogramThreshold(TrigPixID);
	TH1F hTraceCFD = GetTraceHistogramCFD(TrigPixID);

	
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
void  TriggerTelescopeBase::LoadEvent(TelescopeData *TelData)

{

 telData = TelData;

 telData->bTriggeredTriggerPixels.assign(iNumTriggerPixels,kFALSE);

 telData->fDiscriminatorTime.assign(iNumTriggerPixels,-1e6);

 //can be removed is only needed to check the Time over threshold of all pixels in a trigger.
 telData->fTimeOverThreshold.assign(iNumTriggerPixels,0);

 //The start sample for the delayed trace;
 Int_t iStartSample = (int)(fDiscDelay/fSamplingTime)+1;

 if(bDebug)
     {
        cout<<"Start sample for delayed pulse "<<iStartSample<<endl;
     }


 for(Int_t g=0;g<iNumTriggerPixels;g++)
    {
      fTracesInTriggerPixels[g].assign(telData->iNumSamplesPerTrace,0.0);
      fTracesInTriggerPixelsConstantFraction[g].assign(telData->iNumSamplesPerTrace,0.0);
      //loop over all group members and add their trace to the trace of the sumgroup
      for(UInt_t n = 0; n<iTrigPixMembers[g].size();n++)
         {      
           Int_t memberID = iTrigPixMembers[g][n];
           for(Int_t i = 0 ; i<telData->iNumSamplesPerTrace; i++)
	          {
	             Float_t fsignal = telData->fTraceInPixel[memberID][i]*fPEtomVConversion ;
                 //cout<<fsignal<<"  "<<telData->fTraceInPixel[memberID][i]<<"  "<<fPEtomVConversion<<endl;
                 fsignal = fsignal < fClippingLevel && bDoClipping == kTRUE   ? fClippingLevel : fsignal;	      
	             Float_t fsigDelayed = i<iStartSample ? telData->fTraceInPixel[memberID][0]*fPEtomVConversion : telData->fTraceInPixel[memberID][i-iStartSample]*fPEtomVConversion;
                 fsigDelayed = fsigDelayed > fClippingLevel && bDoClipping == kTRUE  ? fClippingLevel : fsigDelayed;	      
	      
	             Float_t CFDsignal = fsignal*fDiscConstantFractionAttenuation-fsigDelayed;
	             fTracesInTriggerPixels[g][i]+=fsignal;
	             fTracesInTriggerPixelsConstantFraction[g][i]+=CFDsignal;
	         }
            if(bDebug)
                {
                     cout<<"Summed group "<<g<<" added pixel "<<memberID<<endl; 
	                // ShowTrace(g);
               }
           
	     }
    }
  

  //cout<<"Finished loading the event"<<endl;

}


//-----------------------------------------------------------------------------------------------------------------------
//This does the actual triggering
//It returns true if the telescope has triggered and false if not
Bool_t  TriggerTelescopeBase::RunTrigger()
{
   //cout<<"TriggerTelescopeBase::RunTrigger"<<endl;
   //Run the CFD simulation
   telData->iNumTriggeredTriggerPixels=0;
   for(Int_t i=0;i<iNumTriggerPixels;i++)
    {
      if(RunDiscriminator(i))
	     telData->iNumTriggeredTriggerPixels++;
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



//-----------------------------------------------------------------------------------------------------------------------
// Does the constant fraction discrimination for a sum group 
//This is a two level discriminator
//The first one is a threshold discriminator that stays up as long as the signal is above threshold
//The second one gives an output if the sum of the attenuated (0.4) signal plus offset and the delayed and inverted copy of the
//input signal cross zero 
//The whole thing gives an output if both discriminators give logic one
Bool_t TriggerTelescopeBase::RunDiscriminator(Int_t TrigPixID)
{

  //cout<<"TriggerTelescopeBase::RunDiscriminator"<<endl;

  Int_t iStartSample = (int)(fDiscDelay/fSamplingTime)+1;

  //1. Threshold discriminator: output as long as signal is above threshold
  //2. Output as long if  CFD trace attenuated copy plus inverted copy plus constant offset output is negativ
  //3. Trigger if 1. and 2. are high
  telData->bTriggeredTriggerPixels[TrigPixID] = kFALSE;
  telData->fDiscriminatorTime[TrigPixID] = -1e6;

  for(Int_t i = iStartSample; i<telData->iNumSamplesPerTrace; i++)
    {
	  Float_t fsignal = fTracesInTriggerPixels[TrigPixID][i];

	  Bool_t bCFDOut = kTRUE; //always true if we do not use the CFD part
      if(bDiscCFDUsage) bCFDOut=fTracesInTriggerPixelsConstantFraction[TrigPixID][i]>=0 ? kTRUE : kFALSE;
          //cout<<i<<"  "<<fsignal<<" threshold "<<fDiscThreshold<<" "<<bCFDOut<< endl;
	
	  //In case both discriminators are firing we have a trigger
	  if( bCFDOut && fsignal<=fDiscThreshold &&  telData->bTriggeredTriggerPixels[TrigPixID] == kFALSE ) 
            {
	      telData->fDiscriminatorTime[TrigPixID] = i*fSamplingTime-fStartSamplingBeforeAverageTime;
	      telData->bTriggeredTriggerPixels[TrigPixID] = kTRUE;

	      //loop over the rest of the trace to find out how long the discriminator is above threshold
	      //this can be deactivated because it is only used in simulation studies done for the 1m sst.
	      Int_t iToT=0;
	      while(fsignal<=fDiscThreshold && i<telData->iNumSamplesPerTrace-1)
		{
		  iToT++;
		  i++;
		  fsignal = fTracesInTriggerPixels[TrigPixID][i];
		}
              telData->fTimeOverThreshold[TrigPixID] = iToT*fSamplingTime;
       	      break; //exit loop because we triggered the trace and we are only interested in the first trigger, saves some compute time
	    }
    }

     if(bDebug && telData->bTriggeredTriggerPixels[TrigPixID] == kTRUE)
     if(bDebug)
      {
        ShowTrace(TrigPixID);
        cout<<"Trace on which the discriminator was run and triggered; groupid "<<TrigPixID<<endl;
        debugDisplay->AddDiscriminatorTraces(telData->GetTelescopeID(),TrigPixID,fDiscThreshold,GetTraceHistogramThreshold(TrigPixID),GetTraceHistogramCFD(TrigPixID));
                       debugDisplay->ShowSelectedDiscriminatorPixels();

      }

    return telData->bTriggeredTriggerPixels[TrigPixID];
}


//--------------------------------------------------------------------------------------------
//
// returns true if group is in patch
// if patches are not used it always returns true
Bool_t TriggerTelescopeBase::GroupInPatch(Int_t TrigPixID,Int_t PatchID)
{

  if(!bUsePatches)
    return kTRUE;
 
  for(UInt_t i=0;i<vPatch[PatchID].size();i++)
    {
      if(TrigPixID==vPatch[PatchID][i])
	 return kTRUE;
    }

  return kFALSE;
}


//-----------------------------------------------------------------------------------------
//Calculate a bias curve for the current trigger settings
//needs direct acces to tracegenerator to produce NSB by itself
void TriggerTelescopeBase::RunBiasCurve(UInt_t Trials,Float_t LowerBoundary,Float_t UpperBoundary,Float_t StepWidth, TraceGenerator *tracegenerator,TelescopeData *TelData)
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
  fTrigPixRateVsThreshold.ResizeTo(NumScanPoints);
  fTrigPixRateVsThreshold.Zero();
  fTrigPixRateVsThresholdErr.ResizeTo(NumScanPoints);
  fTrigPixRateVsThresholdErr.Zero();

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
      //cout<<"Number of groups that triggered "<<telData->GetNumTriggeredTriggerPixels()<<endl;
  	  if( telData->GetNumTriggeredTriggerPixels()!=0 )
	    {
	      fTrigPixRateVsThreshold[t]+=telData->GetNumTriggeredTriggerPixels();
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
//	  cout<<"Dynamic Value of the RFB: "<<fDiscRFBDynamic<<" mV. The rate of zero crossings in MHz: "<< lZeroCrossings /(i* (fTraceLength-fDiscDelay)*1e-3*iNumTriggerPixels)<<endl;
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

      fTrigPixRateVsThresholdErr[t]=sqrt(fTrigPixRateVsThreshold[t])/iNumTriggerPixels/(Trials*fTraceLength*1e-9);
      fTrigPixRateVsThreshold[t]=fTrigPixRateVsThreshold[t]/iNumTriggerPixels/(Trials*fTraceLength*1e-9);
      cout<<"NSB Telescope Trigger rate at "<<LowerBoundary+StepWidth*t<<" mV threshold "<<fBiasCurve[t]<<"+-"<<fBiasCurveErr[t]<<" Hz"<<endl;
      cout<<"Group rate :"<<fTrigPixRateVsThreshold[t]<<"+-"<<fTrigPixRateVsThresholdErr[t]<<" Hz"<<endl;
    }
   


}

//----------------------------------------------------------------------------------------------------------
//
// Set the threshold mV of the discriminator and the width of the output pulse (ns) 
void TriggerTelescopeBase::SetDiscriminatorThresholdAndWidth(Float_t threshold, Float_t width)
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
// Set the delay of the inverted signal in the discriminator and the attenuation of the non-inverted signal
// 
void TriggerTelescopeBase::SetDiscriminatorDelayAndAttenuation(Float_t delay, Float_t attenuation)
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

//----------------------------------------------------------------------------------------
//
// An L2 trigger that does care about patches
Bool_t  TriggerTelescopeBase::RunL2Patch(Int_t PatchNumber,Float_t *fPatchTriggerTimes)
{
  if(bDebug)
   {
     cout<<"Running L2 Patch "<<PatchNumber<<endl;
   }

  if(iClusterID)
    delete [] iClusterID;
  iClusterID = new Int_t[iNumTriggerPixels];

  Int_t *iNumTriggerPixelsInCluster = new Int_t[iNumTriggerPixels]; //How many groups are in each cluster of this patch

  for(Int_t i=0;i<iNumTriggerPixels;i++)
    {
      iClusterID[i] = -1;
      iNumTriggerPixelsInCluster[i]=0; 
    }

  vector< int > i_pix;
  vTrigPixsInCluster = new vector< vector<int> >;
  vTrigPixsInCluster->assign(iNumTriggerPixels , i_pix );

  vector< float > f_pix;
  vTriggerTimesInCluster = new vector< vector<float> >;
  vTriggerTimesInCluster->assign(iNumTriggerPixels , f_pix );

  //cout<<"allocated arrays and vectors"<<endl;
  //cout<<"Number of summed groups "<<iNumTriggerPixels<<endl;

  Int_t NTrigPixs = bUsePatches==kTRUE ? vPatch[PatchNumber].size() :  iNumTriggerPixels;
  if(bDebug)
  cout<<"Number of trigger pixels in this patch: "<<NTrigPixs<<endl;
  for(Int_t i=0;i<NTrigPixs;i++)
    {
      Int_t TrigPixID = bUsePatches==kTRUE ? vPatch[PatchNumber][i] : i;
    
      if (telData->bTriggeredTriggerPixels[TrigPixID]==kTRUE) //(summed) pixel has triggered
	{
	  iNumTriggerPixelsInCluster[TrigPixID] = CalcCluster(TrigPixID,TrigPixID,PatchNumber);
	  // cout<<iNumTriggerPixelsInCluster[i]<<endl;       
	}
	  
    }

  //cout<<"Finished finding all clusters"<<endl<<endl;
   
  //Sort all cluster descending
  Int_t* index = new Int_t[iNumTriggerPixels];
  TMath::Sort(iNumTriggerPixels,iNumTriggerPixelsInCluster,index);
 
/*  for(int n=0;n<iNumTriggerPixels;n++)
	 cout<<iNumTriggerPixelsInCluster[index[n]]<<"  ";
	 
	 cout<<endl;
  */
  //Find the cluster that triggered first
  Float_t fPatchTriggerTime = 1e6;
  Int_t t = 0;
  Int_t tmin = 0;
  while(iNumTriggerPixelsInCluster[index[t]]>=iMultiplicity )
    {
      Int_t igroups = iNumTriggerPixelsInCluster[index[t]];

      if(bDebug)
      cout<<"cluster "<<index[t]<<" has "<<igroups<<" groups"<<endl;

      //Find the L2 trigger time of this cluster
      Float_t *times = new Float_t[ igroups ];
      for(Int_t i=0;i<igroups;i++)
	  {
           times[i] = (vTriggerTimesInCluster->at(index[t]))[i];
           if(bDebug)
           cout<<"triggered group in cluster "<<vTrigPixsInCluster->at(index[t])[i]<<" t: "<<times[i]<<endl;
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

       if(t==iNumTriggerPixels)
          break;

    }
  fPatchTriggerTime+= telData->fAveragePhotonArrivalTime;
  
  fPatchTriggerTimes[PatchNumber] = fPatchTriggerTime;

//move the vector with the groups of pixels that are in the trigger cluster to a separate vector
  //which can be retrieved with GetTriggerCluster()
  //cout<< vGroupsInCluster->at(index[tmin]).size()<<endl;
  //cout<<(telData->vTriggerCluster).size()<<endl;
  iPixelTriggeredInPatch[PatchNumber]=vTrigPixsInCluster->at(index[tmin]);
  //cout<<"Have transferred the TriggerCluster"<<telData->vTriggerCluster.size()<<endl;
  //Now trigger if we have enough groups in the largest cluster
 

  Bool_t bPatchTrigger = kFALSE;

  if(iNumTriggerPixelsInCluster[index[0]]>=iMultiplicity)
    {
      //cout<<iNumTriggerPixelsInCluster[index[0]]<<endl;
      bPatchTrigger= kTRUE;
    }

  delete [] iNumTriggerPixelsInCluster;
  delete vTrigPixsInCluster;
  delete vTriggerTimesInCluster;
  delete [] index;

  return bPatchTrigger;
}

//------------
//
// An L2 trigger that does care about patches
Bool_t  TriggerTelescopeBase::RunL2WithPatches()
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

// --------------------------------------------------------------------------------------------
//
//Helper function to loop over all neighbours to find triggered pixels belonging to the cluster
//returns the number of pixel in that cluster

Int_t TriggerTelescopeBase::CalcCluster(Int_t TrigPixID, Int_t ClusterID, Int_t PatchID)
{

  if(bDebug)
	  cout<<"visiting group "<<TrigPixID<<endl;

  // If we have visited this group in this round ... do nothing.
  if (iClusterID[TrigPixID]==ClusterID)
    return 0;
  // Assign the new cluster ID this Group
  iClusterID[TrigPixID]=ClusterID;
  //Put the group in the list of groups that belong to this cluster centered around the group with ID ClusterID
  vTrigPixsInCluster->at(ClusterID).push_back(TrigPixID);
  
  
  //Add the trigger time of the group to the vector of this cluster
  vTriggerTimesInCluster->at(ClusterID).push_back(telData->fDiscriminatorTime[TrigPixID]);
  if(bDebug)
  {
	 cout<<"Pixel number: "<<vTriggerTimesInCluster->at(ClusterID).size()<<endl;
	 cout<<"added group "<<TrigPixID<<", which triggered at time "<<telData->fDiscriminatorTime[TrigPixID]<<endl;
  }

  // Need this to store the number of groups in the cluster
  Int_t NumTrigPixsInCluster = 1;

  // Now do the same with all its neighbors
  for(UInt_t n = 0; n<iTrigPixNeighbors[TrigPixID].size();n++)
    {
      Int_t TrigPixIDNeighbor = iTrigPixNeighbors[TrigPixID][n];
      if (telData->bTriggeredTriggerPixels[TrigPixIDNeighbor] 
           && fabs(telData->fDiscriminatorTime[ClusterID]-telData->fDiscriminatorTime[TrigPixIDNeighbor])<fWidthDiscriminator 
           && GroupInPatch(TrigPixIDNeighbor,PatchID))
	    NumTrigPixsInCluster += CalcCluster(TrigPixIDNeighbor,ClusterID,PatchID);
    }
 // return the number of groups in this cluster
  return NumTrigPixsInCluster;

}


//Reads in  the config file and sets all variables
void   TriggerTelescopeBase::SetParametersFromConfigFile(ReadConfig *readConfig ){

   cout <<endl<< "TriggerTelescopeBase::SetParametersFromConfigFile " << endl;

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
    
   //Do we use the CFD part of the discriminator?
   bDiscCFDUsage = readConfig->GetCFDUsage(iTelType);
   cout<<"Do we use the CFD circuit: "<<bDiscCFDUsage<<endl;

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
         cout<<"Will confine pixels to groups"<<endl;
        iNumTriggerPixels = readConfig->GetNumberGroups(iTelType);
        iTrigPixNeighbors = readConfig->GetNeighborsOfGroup(iTelType);
        iTrigPixMembers = readConfig->GetMembersOfGroups(iTelType);
        if(iNumTriggerPixels <=0)
          {
             cout<<"The sumtrigger simulation was turned on but the number of summed groups has been put to "<<iNumTriggerPixels<<endl;
             cout<<"do something about it!!"<<endl;
             exit(1);
          }
     }
   else //If no sumtrigger or groups are used
     {
        cout<<"This simulation does not use a sumtrigger"<<endl;
		iNumTriggerPixels= readConfig->GetNumberPixels(iTelType);
		iTrigPixNeighbors = readConfig->GetNeighbors(iTelType);
		vector<int> vmembers;
        for(int i = 0; i<iNumTriggerPixels; i++)
         {
			iTrigPixMembers.push_back(vmembers);
            iTrigPixMembers[i].assign(1,i);
         }
      }
      
   SetDiscriminatorThresholdAndWidth(fDiscThreshold,fWidthDiscriminator);
   SetDiscriminatorDelayAndAttenuation(fDiscDelay, fDiscConstantFractionAttenuation);
  //Create the traces for the summed pixels
  CreateTraces();
}

