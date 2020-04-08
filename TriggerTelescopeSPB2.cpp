/* \file TriggerTelescopeVERITAS.cpp
   Implementation of a next neighbor trigger for VERITAS
   1. Signal of X pixels are summed
   2. Discrimination of summed signal
   3. next neighbor logic (multiplicity of Y) 
   4. Implementaton of the rate feedback in the discriminator path.
*/

#include "TriggerTelescopeVERITAS.h"
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
TriggerTelescopeVERITAS::TriggerTelescopeVERITAS(ReadConfig *readConfig, int telType, TRandom3 *generator,Bool_t debug,Display *display) : TriggerTelescopeBase(readConfig,telType,generator,debug,display)
{
  cout<<"Initializing the Telescope Trigger for VERITAS"<<endl;
  fDiscRFBConstant = 0;  
  fDiscRFBDynamic = 0;
  bDiscRFBUsage= 0;

  lZeroCrossings = 0;               //holds the number of zerocrossings for one event counted over all pixels in the camera
  lNumEvents = 0;   

  SetParametersFromConfigFile(readConfig);	
}

//-----------------------------------------------------------------------------------------------------------------------
//Loads the traces from the pixels (sums them up)
void  TriggerTelescopeVERITAS::LoadEvent(TelescopeData *TelData)

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


//Needs to be moved to VERITAS Trigger class
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


//This does the actual triggering
//It returns true if the telescope has triggered and false if not
Bool_t  TriggerTelescopeVERITAS::RunTrigger()
{
    //cout<<"TriggerTelescopeVERITAS::RunTrigger"<<endl;
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




//---------------------------------------------------------------------------------------
//Get the numbers of zerocrossings for the loaded event from negative to positive
Long_t TriggerTelescopeVERITAS::GetNumZeroCrossings()
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


//------------------------------------------------------------------------------------------------------
//
// Set the rate feedback of the discriminator
// Putting it to 0 will turn it off 
void TriggerTelescopeVERITAS::SetDiscriminatorRFBConstant(Float_t rfb)
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
void TriggerTelescopeVERITAS::SetDiscriminatorRFBDynamic(Float_t rfb)
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
void TriggerTelescopeVERITAS::SetDiscriminatorRFBUsage(Bool_t rfbuse)
{ 

  bDiscRFBUsage = rfbuse; 

}


//Reads in  the config file and sets all variables
void   TriggerTelescopeVERITAS::SetParametersFromConfigFile(ReadConfig *readConfig ){

   cout <<endl<< "TriggerTelescopeVERITAS::SetParametersFromConfigFile " << endl;

   //Do we use the RFB circuit?
   bDiscRFBUsage = readConfig->GetRFBUsage(iTelType);
   cout<<"Do we use the RFB circuit: "<<bDiscRFBUsage<<endl;

   //The Constant Value in the RFB in the discriminator
   fDiscRFBConstant = readConfig->GetDiscriminatorRFBConstant(iTelType);
   cout<<"If the RFB circuit is used this is the magnitude of the rate dependend feedback  in mV/MHz "<<fDiscRFBConstant<<endl;
       
   //The dyamic value in the RFB in the discriminator
   fDiscRFBDynamic = readConfig->GetDiscriminatorRFBDynamic(iTelType);
   cout<<"The initial dynamic value in the RFB in mV (will be mulitplied with 0.18 in the sims). If RFB circuit is not used this is not used "<<fDiscRFBDynamic<<endl;
     
   SetDiscriminatorRFBConstant(fDiscRFBConstant);
   SetDiscriminatorRFBDynamic(fDiscRFBDynamic);
}

