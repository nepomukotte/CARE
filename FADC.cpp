/* \file FADC.cpp
   Implementation of the FADC for VERITAS
*/

#include "FADC.h"
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
FADC::FADC( ReadConfig *readConfig, TraceGenerator *traceGenerator, int telType, TRandom3 *generator,Bool_t debug,Display *display)
{
  cout<<"Initializing the FADC"<<endl;

  bDebug = debug;

  debugDisplay = display;
 
  rand = generator;
  tracegenerator = traceGenerator;

  iTelType = telType;

  iNumTimesOutOfAnalogueTraceUpperEnd = 0;
  iNumTimesOutOfAnalogueTraceLowerEnd = 0;

  //Read the config file
  SetParametersFromConfigFile( readConfig );
  
}


////////////////////////////////////////////////////////////////////////////////////
// Does the FADC simulations
void FADC::RunFADC(TelescopeData *TelData)
{

  telData = TelData;

  //Check if we have all we need, e.g. that the analog trace is simulated sufficiently long and can be fully digitized

  //Take care of the case that the array triggers but this telescope has no Cherenkov photons in the focal plane
  if( telData->bTelescopeHasTriggered )
    fTimeStartFirstSample = telData->fTriggerTime - fOffset + ( rand->Uniform() - 0.5 ) *  fFADCSamplingWidth;
  else if( telData->bCherenkovPhotonsInCamera == kTRUE )
    fTimeStartFirstSample = telData->fAveragePhotonArrivalTime - fOffset + ( rand->Uniform() - 0.5 ) *  fFADCSamplingWidth;
  else
    fTimeStartFirstSample = - fOffset + ( rand->Uniform() - 0.5 ) *  fFADCSamplingWidth;

  if(bDebug)
	 cout<<"The time of the first sample to be readout is: "<<fTimeStartFirstSample<<endl;

  telData->bInLoGain.assign(telData->iNumPixels,kFALSE);

  //loop over pixels and digitize signals
   for(Int_t g=0;g<telData->iNumPixels;g++)
    {

    

        DigitizePixel( g );
        if(bDebug)
		  {
			  cout<<"done pixel "<<g<<endl; 
		      //debugDisplay->Show(telData->GetTelescopeID(),g);
	          //ShowTraceOfPixel(g);
		  }

    }
  

}

//////////////////////////////////////////////////////////////////////////////////
// Digitize the analog trace in one pixel
void FADC::DigitizePixel( Int_t PixelID )
{
  

    //clear trace
    telData->iFADCTraceInPixel[PixelID].assign(iFADCSamples,0);

    //Check if the Lo Gain is active
    Float_t fGain = 1;
    Bool_t bLowGain = kFALSE;
    Float_t LowGainthresholdInPE =  (fHiLoGainThreshold - fPedestal) / fDCtoPEconversion ;

    //At the same time do the QDC
    Float_t fQDC = 0;

    if(bDebug)
      {
      cout<<endl<<"Pixel: "<<PixelID<<", analog samples in trace  "<<iFADCSamples*fFADCSamplingWidth/fTraceSamplingTime
          <<", time of first sample  "<<fTimeStartFirstSample<<" ns , average photon arrival time [ns] "<<telData->fAveragePhotonArrivalTime<<endl;
      cout<<"Sample; Pos in analog trace; amplitude -pe; ampl in DC; digitized value after cutting to dynamic range"<<endl;
      }

    for(int i =0;i<iFADCSamples*fFADCSamplingWidth/fTraceSamplingTime;i++)
      {
  	     Int_t  iPositionInAnalogTrace = (Int_t)
	     (( fTimeStartFirstSample +  i * fTraceSamplingTime - (telData->fAveragePhotonArrivalTime - fStartSamplingBeforeAverageTime) ) 
	        / fTraceSamplingTime);

	     if(iPositionInAnalogTrace>=(Int_t)(telData->fTraceInPixel[PixelID].size()) )
               {                         
                   iPositionInAnalogTrace = iPositionInAnalogTrace  % (Int_t)telData->fTraceInPixel[PixelID].size();
               }
             else if(iPositionInAnalogTrace<0)
               {
                  iPositionInAnalogTrace = abs((Int_t)telData->fTraceInPixel[PixelID].size()-iPositionInAnalogTrace) % (Int_t)telData->fTraceInPixel[PixelID].size() ;                  
               }

	     if( -1 * telData->fTraceInPixel[PixelID][iPositionInAnalogTrace] > LowGainthresholdInPE ) 
	       {   
	          if(bDebug)
	              cout<<"Do low gain for Pixel "<<PixelID<<endl;
	          fGain = fLowHiGainRatio;
	          bLowGain = kTRUE;
	          telData->bInLoGain[PixelID] = kTRUE;
	      }

         fQDC += telData->fTraceInPixel[PixelID][iPositionInAnalogTrace];
      }

    fQDC*=-1*fTraceSamplingTime*fDCtoPEconversion*tracegenerator->GetAreaToPeakConversion();

    telData->iQDCInPixel[PixelID] = (Int_t)(fQDC+fPedestal);

    if(bDebug)
	  cout<<"QDC: "<<fQDC<<" before pedestal additio,n "<<telData->iQDCInPixel[PixelID]<<" after pedestal addition"<<endl;

    //Get the right trace
    vector<Float_t> trace;
    if(bLowGain)
      {
        tracegenerator->SetTelData(telData);
        trace = tracegenerator->GetLowGainTrace(PixelID);
      }
    else
      trace = telData->fTraceInPixel[PixelID];

    //Write the FADC trace
    Float_t fConversionFactor = fGain * fDCtoPEconversion  ;
    bool outofbound = false;
    for(int i =0;i<iFADCSamples;i++)
      {
	
	    Int_t  iPositionInAnalogTrace = (Int_t)
	      ( ( fTimeStartFirstSample +  i * fFADCSamplingWidth -(telData->fAveragePhotonArrivalTime-fStartSamplingBeforeAverageTime) ) 
	      / fTraceSamplingTime);

        //Float_t fAnalogValue = iPositionInAnalogTrace >= (Int_t)telData->fTraceInPixel[PixelID].size() ? 0 : trace[iPositionInAnalogTrace];   
	   if(iPositionInAnalogTrace >= (Int_t)telData->fTraceInPixel[PixelID].size())
		  {
                       if(! outofbound && PixelID==1)
                          {
                            iNumTimesOutOfAnalogueTraceUpperEnd++;
                            outofbound = true;
	                    cout<<"Simulate a longer analog trace at least a length of "<<fTimeStartFirstSample + iFADCSamples*fFADCSamplingWidth  -(telData->fAveragePhotonArrivalTime-fStartSamplingBeforeAverageTime) <<endl;
			   cout<<endl<<"Pixel: "<<PixelID<<", samples in trace  "<<iFADCSamples*fFADCSamplingWidth/fTraceSamplingTime
			        <<", time of first sample  "<<fTimeStartFirstSample<<", average photon arrival time "<<telData->fAveragePhotonArrivalTime<<endl;
			    cout<<"E: "<<fenergy<<" Tel "<<ftelid<<" Zenith "<<fzenith<<" Az  "<<fazimuth<<endl<<endl;
			    cout<<"I am going back to the very beginning of the trace and record from there"<<endl;
                            cout<<"Lets hope that there are no Cherenkov photons"<<endl;
                          }
                          iPositionInAnalogTrace = iPositionInAnalogTrace  % (Int_t)telData->fTraceInPixel[PixelID].size();
		  }
             else if(iPositionInAnalogTrace<0)
               {
                  if( !outofbound && PixelID==1 )
                     {
                       iNumTimesOutOfAnalogueTraceLowerEnd++;
                       outofbound = true;
                       cout<<endl<<"The readout start before sample 0. Wanted position in analog trace is "<<iPositionInAnalogTrace<<endl;
                       cout<<"I go to the end of the trace and load the part of the trace there, hopefully there is only noise"<<endl;
		       cout<<"E: "<<fenergy<<" Tel "<<ftelid<<" Zenith "<<fzenith<<" Az  "<<fazimuth<<endl<<endl;
                     }
                  iPositionInAnalogTrace = abs((Int_t)telData->fTraceInPixel[PixelID].size()-iPositionInAnalogTrace) % (Int_t)telData->fTraceInPixel[PixelID].size() ;                  
               }
	 
	   Float_t fDigitizedValue =  -1 * trace[iPositionInAnalogTrace] * fConversionFactor + fPedestal;

	    if(bDebug)
	       cout<<i<<"  "<<iPositionInAnalogTrace<<"  "<<trace[iPositionInAnalogTrace]<<"  "<<fDigitizedValue;

	    //convert to integer and check if out of FADC dynamic range
	    fDigitizedValue = fDigitizedValue > iDynamicRange ? iDynamicRange : (Int_t)fDigitizedValue;
	    fDigitizedValue = fDigitizedValue < 0 ? 0 : (Int_t)fDigitizedValue;

	    if(bDebug)
	      cout<<"  "<<fDigitizedValue<<endl;

	    telData->iFADCTraceInPixel[PixelID][i] = (Int_t)fDigitizedValue;

      }
}

//////////////////////////////////////////////////////
//
// Show the digitized trace of one pixel
TH1F FADC::GetTraceHistogram(Int_t PixNum){

  TString title;
  title.Form("Trace in pixel %i",PixNum);
  TH1F hTrace ("hTrace",title,
	       iFADCSamples,-0.5,iFADCSamples-0.5);
  hTrace.GetXaxis()->SetTitle("Sample Number");
  hTrace.GetYaxis()->SetTitle("Digitized Amplitude [dc]");
  
  for(Int_t i = 0; i<iFADCSamples; i++)
    {
   
      hTrace.SetBinContent(i+1,telData->iFADCTraceInPixel[PixNum][i]);
      cout<<telData->iFADCTraceInPixel[PixNum][i]<<endl;
    }
  
 return  hTrace;
 
}



//////////////////////////////////////////////////////
//
// Show the digitized trace of one pixel
void FADC::ShowTraceOfPixel(Int_t PixNum){

  cout<<"Outputing FADC trace for pixel "<<PixNum<<endl;

  TCanvas cTraceTriggeredGroup("cTraceTriggeredGroup","Digitized Trace in pixel",700,500);
  cTraceTriggeredGroup.cd();
  gPad->SetGrid();
 
  TH1F hTrace = GetTraceHistogram(PixNum);
  
  hTrace.Draw();
  
  TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
  
  //
  // While reading the input process gui events asynchronously
  //
  timer.TurnOn();
  TString input = Getline("Type <return> to go on: ");
  timer.TurnOff();  

}

//----------------------------------------------------------------------------------------
//  Prints how often the trace was too short
void   FADC::PrintHowOftenTheTraceWasTooShort(){

  cout<<"The trace was too short at the high end: "<<iNumTimesOutOfAnalogueTraceUpperEnd<<" times for the FADC"<<endl;
  cout<<"The trace was too short at the lower end: "<<iNumTimesOutOfAnalogueTraceLowerEnd<<" times for the FADC"<<endl;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Reads in  the config file and sets all variables
void  FADC::SetParametersFromConfigFile( ReadConfig *readConfig ){

   cout <<endl<< "FADC::SetParametersFromConfigFile " << endl;

   fTraceLength=-1;
   fStartSamplingBeforeAverageTime=-1;
   fTraceSamplingTime=-1;
   fFADCSamplingWidth=-1;
   iFADCSamples=-1;
   iDynamicRange=-1;
   fOffset=-1;
   fDCtoPEconversion=-1;
   fHiLoGainThreshold=-1;
   fLowHiGainRatio=-1;
   fPedestal=-1;      

   //Length of the trace
   fTraceLength = readConfig->GetTraceLength(iTelType);

   //Offset from average photon arrival time
    fStartSamplingBeforeAverageTime = readConfig->GetStartSamplingTimeOffsetFromAveragePhotonTime(iTelType);

   //Width of one sample  of the trace
   fTraceSamplingTime = readConfig->GetSamplingTime(iTelType);

   //the sampling width of the FADC
   fFADCSamplingWidth = readConfig->GetFADCSampleTime(iTelType);

   //the number of sample in a FADC trace
   iFADCSamples = readConfig->GetFADCSamples(iTelType);

   //the dynamic range of the FADC
   iDynamicRange = readConfig->GetFADCDynamicRange(iTelType);

   //the Offset of the beginning of the readout window from the trigger decision
   fOffset = readConfig->GetFADCOffsetFromTrigger(iTelType);

   //the conversion from DC to PE
   fDCtoPEconversion = readConfig->GetFADCDCtoPEconversionFactor(iTelType);

   //the threshold when the low gain Channel is activated
   fHiLoGainThreshold = readConfig->GetFADCHiLoGainThreshold(iTelType);

   //the gain ratio of the low gain channel / high gain channel 
   fLowHiGainRatio = readConfig->GetFADCLowHiGainRatio(iTelType);

   //the FADC pedestal 
   fPedestal = readConfig->GetFADCPedestal(iTelType);

    cout<<endl<<"FADC settings"<<endl;
    cout<<"-------------"<<endl;
	cout<<"Trace length in ns set to "<<fTraceLength<<endl;
        cout<<"The analog trace starts to be sampled "<<fStartSamplingBeforeAverageTime<<" ns before the average photon arrival time"<<endl;
	cout<<"Width of one sample in ns of the Trace set to "<<fTraceSamplingTime<<endl;
	cout<<"The sampling time of the FADC in ns is: "<<fFADCSamplingWidth<<endl;
	cout<<"The number of FADC samples is: "<<iFADCSamples<<endl;
	cout<<"The dynamic range of the FADC is: "<<iDynamicRange<<endl;
	cout<<"The offset between the trigger and the beginning of the readout [ns]: "<<fOffset<<endl;
	cout<<"The DC to PE conversion factor is: "<<fDCtoPEconversion<<endl;
	cout<<"The Low gain is activated at [dc]: "<<fHiLoGainThreshold<<endl;
	cout<<"The gain ratio of the low gain channel to the high gain channel is: "<<fLowHiGainRatio<<endl;
	cout<<"The FADC pedestal is [dc]: "<<fPedestal<<endl<<endl;;      
}

void FADC::SetDebugInfo(Float_t energy, Int_t telid, Float_t zenith, Float_t azimuth){

fenergy = energy;
ftelid = telid;
fzenith = zenith;
fazimuth = azimuth;

}
