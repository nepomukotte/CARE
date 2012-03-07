/* \file TelescopeData.cpp
   Implementation of a telescope data class
   This class stores the Traces of all telescope that contain NSB and Cherenkov photons
   It also stores all the trigger bits and other information that is needed to write a raw
   data file etc.
*/

#include "TelescopeData.h"
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
TelescopeData::TelescopeData( ReadConfig *readConfig, Int_t telID, Bool_t debug)
{

  bDebug = debug;

  iTelID = telID;

  iNumPixels = 0;
  fTraceInPixel = NULL;
  fTraceInPixelNSBOnly = NULL;
  fTimesInPixel  = NULL;
  fAmplitudesInPixel = NULL;
  iNumSamplesPerTrace = -1; 

  iFADCTraceInPixel = NULL;
  iNumFADCSamples = -1;

  TelXpos = 0;
  TelYpos = 0;
  TelZpos = 0;
  OpticsTransitTime = 0;

  vTriggerCluster.resize(0);

  bTelescopeHasTriggered = kFALSE;                        //If telescope has triggered

  SetParametersFromConfigFile( readConfig );

  SetupArrays();

}


///////////////////////////////////////////////////
//
//  Creating the arrays for the traces

void   TelescopeData::SetupArrays(){



  //delete traces from previous event
  if(fTraceInPixel)
    {
	delete []  fTraceInPixel;
    }
  //create the array with all the vectors, one vector for each pixel
  fTraceInPixel = new vector<Float_t>[iNumPixels];

  if(fTraceInPixelNSBOnly )
    {
	delete []  fTraceInPixelNSBOnly ;
    }
   //create the array with all the vectors, one vector for each pixel
  fTraceInPixelNSBOnly  = new vector<Float_t>[iNumPixels];


  if(fTimesInPixel )
    {
      delete []  fTimesInPixel ;
    }
  //create the array with all the vectors, one vector for each pixel
  fTimesInPixel  = new vector<Float_t>[iNumPixels];
   
  if(fAmplitudesInPixel )
    {
      delete []  fAmplitudesInPixel ;
    }
  //create the array with all the vectors, one vector for each pixel
  fAmplitudesInPixel  = new vector<Float_t>[iNumPixels];

  if(iFADCTraceInPixel)
    {
    delete []  iFADCTraceInPixel;
    }
  //create the array with all the vectors, one vector for each pixel
      iFADCTraceInPixel = new vector<Int_t>[iNumPixels];

 
   ResetTraces();

}


///////////////////////////////////////////////////////////////////////////////
// 
// Reset all the traces in the container to zero, to be ready for a new event
//
//////////////////////////////////////////////////////////////////////////////

void TelescopeData::ResetTraces()
{

  fAveragePhotonArrivalTime = 0.0;
  mean = 0.0;
  bCherenkovPhotonsInCamera = kFALSE;

  iQDCInPixel.assign(iNumPixels,0);
  bInLoGain.assign(iNumPixels,kFALSE);
  bTelescopeHasTriggered = kFALSE;

  //Clear the trace arrays
  for(Int_t g=0;g<iNumPixels;g++)
    {
      fTraceInPixel[g].assign(iNumSamplesPerTrace,0.0);
      fTraceInPixelNSBOnly[g].assign(iNumSamplesPerTrace,0.0);
      fTimesInPixel[g].clear();
      fAmplitudesInPixel[g].clear();
      iFADCTraceInPixel[g].assign(iNumFADCSamples,0);
    }

}



//----------------------------------------------------------------------------------------
//Reads in  the config file and sets all variables
void   TelescopeData::SetParametersFromConfigFile( ReadConfig *readConfig ){

   cout << "TelescopeData::SetParametersFromConfigFile " << endl;
   cout<<  "------------------------------"<<endl;
   cout<< "Telescope ID: "<<iTelID<<endl;


   iTelType = readConfig->GetTelescopeType(iTelID);
   cout<< "Telescope type: "<<iTelType<<endl;


   iTelIDinSuperArray = readConfig->GetTelescopeIDinSuperArray(iTelID);
   cout<<"Telescope ID in superarray: "<<iTelIDinSuperArray<<endl;

   iNumPixels = readConfig->GetNumberPixels(iTelType);


   //The sampling width used in the sample single pe pulse
   Float_t fSamplingTimeAveragePulse = readConfig->GetSampleWidthAveragePulse(iTelType);
   cout<<"The sampling width used in the  single pe pulse in ns is "<<fSamplingTimeAveragePulse<<endl;

   //Length of the trace
   Float_t fTraceLength = readConfig->GetTraceLength(iTelType);
   cout<<"Trace length in ns set to "<<fTraceLength<<endl;

   //Width of one sample  of the trace
   Float_t fSamplingTime = readConfig->GetSamplingTime(iTelType);
   if(fSamplingTime<fSamplingTimeAveragePulse)
       {
   	      cout<<"Sampling time has to be larger than the sampling time of the single pe pulse "<<fSamplingTimeAveragePulse<<endl;
           exit(0);
       }
   cout<<"Width of one sample in ns of the Trace set to "<<fSamplingTime<<endl;

   //Setting the number of samples per trace (This is not the FADC trace)
   iNumSamplesPerTrace =  int(fTraceLength/fSamplingTime);
   cout<<"The number of samples in one trace: "<<iNumSamplesPerTrace<<endl; 

   //the number of sample in a FADC trace
   iNumFADCSamples = readConfig->GetFADCSamples(iTelType);
   
   fRelQE = readConfig->GetRelQE(iTelType);
   fRelQEwWC = fRelQE; 
   fRelGain = readConfig->GetRelGain(iTelType,fRelQE);

   fWinstonConeEfficiency = readConfig->GetWinstonConeEfficiency(iTelID) ;   //The efficiency of the Winstoncone
   cout<<"Winstoncone efficiency: "<<fWinstonConeEfficiency<<endl;
   //add the WinstonConeQE to the rel QE
   for(UInt_t i = 0; i<fRelQEwWC.size();i++)
    fRelQEwWC[i]*=fWinstonConeEfficiency;

   bBlurPSF =   readConfig->GetOpticalPSFBlurBit(iTelID);        
   fBlurSigma = readConfig->GetOpticalPSFBlurSigma(iTelID);     //sigma in mm by which the optical PSF 
                                                               //is blured furthery of th
   cout<<"PSF will be blurred: "<<bBlurPSF<<" with a sigma of [mm] "<<fBlurSigma<<endl;
}
