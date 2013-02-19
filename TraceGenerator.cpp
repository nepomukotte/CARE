/* \file TraceGenerator.cpp
   Implementation of a trace generator for VERITAS
   This class generates Traces that contain NSB and Cherenkov photons which can be used in a trigger simulation
   or output it in VERITAS FADC format
*/

#include "TraceGenerator.h"
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
TraceGenerator::TraceGenerator( ReadConfig *readConfig , int telType, TRandom3 *generator,Bool_t debug,Display *display)
{
  cout<<endl<<endl<<"Initializing the Trace Generator"<<endl;
  cout<<"================================"<<endl<<endl;

  bDebug = debug;
  debugDisplay = display;

  telData = NULL;

  rand = NULL;

  iTelType = telType;

  bSiPM = kFALSE;

  bCrosstalk = kFALSE;                    //Use crosstalk between pixel
  fCrosstalk = 0.0;                    //Crosstalk value

  NumTraceIsTooShortLowEnd=0;
  NumTraceIsTooShortHighEnd=0;
  
  fSamplingTime = -1;
  fSamplingTimeAveragePulse = -1; //sampling time of the average PE pulse
  fNSBRatePerPixel = -1; 
  fStopTimeAveragePulse = 35; //time in ns 
  fStartTimeAveragePulse = 5; //time in ns 
  fTraceLength = 100; //Length of simulated trace in ns
  fStartSamplingBeforeAverageTime=30;
  iNumSamplesAverageSinglePEPulse=-1; 
  bAfterPulsing = kFALSE;
  fFWHMofSinglePEPulse = -1;
  fAreaToPeakConversion = 0;
  //Read the config file
  SetParametersFromConfigFile( readConfig );
  
  rand = generator;
string s = "";                                                                             

gridsearch = new GOrderedGridSearch(fXTubeMM,fYTubeMM,fSizeTubeMM,iTubeSides,fRotAngle,19,19,1,s);
  if(bDebug) cout<<"Done initializing the TraceGenerator"<<endl;

  if(bDebug)
    ShowAverageSinglePEPulse();
}


///////////////////////////////////////////////////
//
//  Setting the TelescopeData container

void   TraceGenerator::SetTelData(TelescopeData *TelData){

telData = TelData;

}


//-----------------------------------------------------------------------------------------------------------------------
//Loads the PEs into the pixels. If UseNSB is set true NSB is added to the trace
//which requires that SetNSBRatePerPixel(Float_t rate) has been set before
void  TraceGenerator::LoadCherenkovPhotons(std::vector< float > *v_f_X,std::vector< float > *v_f_Y,std::vector< float > *v_f_time,std::vector< float > *v_f_lambda, Float_t delay, Double_t dEfficiencyFactor)
{
   

  if(iNumPixels == 0)
    {
      cout<<"You need to set the pixel coordinates first, maybe they are empty"<<endl;
      exit(1);
    }

  if(fSamplingTime <= 0)
    {
      cout<<"You need to set the sampling rate with SetSampleWidth(Float_t t)"<<endl;
      exit(0);
    }

  if(bDebug)
    cout<<"Camera has "<<iNumPixels<<" pixel"<<endl;


    //Load the NSB into the Traces, will be skipped in function if no NSB generation is wanted
    GenerateNSB();
 
  if(bDebug)
   cout<<"NSB done"<<endl;

  //Find the average time of all PEs
  telData->fAveragePhotonArrivalTime = 0;
  Float_t fMinPhotonArrivalTime = 1e6;
  Float_t fMaxPhotonArrivalTime = -1e6;
  for(UInt_t p=0; p<v_f_time->size(); p++)
    {
      telData->fAveragePhotonArrivalTime+=v_f_time->at(p);
      fMinPhotonArrivalTime = fMinPhotonArrivalTime > v_f_time->at(p) ? v_f_time->at(p) : fMinPhotonArrivalTime;
      fMaxPhotonArrivalTime = fMaxPhotonArrivalTime < v_f_time->at(p) ? v_f_time->at(p) : fMaxPhotonArrivalTime;
    }

  if(v_f_time->size()>0)
   telData->fAveragePhotonArrivalTime/=v_f_time->size();

  if(bDebug)
    cout<<"Average Photon arrival time "<<telData->fAveragePhotonArrivalTime<<" with "<<v_f_time->size()<<" Photons "<<endl;

  //Set if we have Cherenkov photons
  telData->bCherenkovPhotonsInCamera = v_f_time->size()>0 ? kTRUE: kFALSE;

  if(bDebug)
  cout<<"Min, Average, Max photon arrival time "<<fMinPhotonArrivalTime<<"  "<<telData->fAveragePhotonArrivalTime<<"  "<<fMaxPhotonArrivalTime<<endl;

   if(telData->fAveragePhotonArrivalTime-fMinPhotonArrivalTime > fStartSamplingBeforeAverageTime)
      NumTraceIsTooShortLowEnd++;
   if(fMaxPhotonArrivalTime - telData->fAveragePhotonArrivalTime > telData->iNumSamplesPerTrace*fSamplingTime-fStartSamplingBeforeAverageTime )
     NumTraceIsTooShortHighEnd++;

    /*{
      cout<<endl<<"Ups the trace is not long enough to save all Cherenkov Photons"<<endl;
      cout<<telData->fAveragePhotonArrivalTime-fMinPhotonArrivalTime<<"  "<<fMaxPhotonArrivalTime - telData->fAveragePhotonArrivalTime <<endl;
      cout<<"We have "<<telData->iNumSamplesPerTrace*fSamplingTime<<" ns, but we need "<<fMaxPhotonArrivalTime-fMinPhotonArrivalTime<<endl; 
    }
*/

  //Load Cherenkov photons into the traces of each summed group
  for(UInt_t p=0; p<v_f_time->size(); p++)
    {

      float fx = v_f_X->at(p);
      float fy = v_f_Y->at(p);


      if(bDebug)
        cout<<endl<<"Photon "<<p<<" position: "<<fx<<"  "<<fy<<endl;

      //add smearing of photon position in focal plane
      if(telData->bBlurPSF==kTRUE)
        {
          fx +=  rand->Gaus(0.0,telData->fBlurSigma);
          fy +=  rand->Gaus(0.0,telData->fBlurSigma);
        }

      if(bDebug)
		   cout<<"position after smearing "<<fx<<"   "<<fy<<endl;

      Int_t pixID = gridsearch->getElemNumber(fx,fy);


      if(pixID<0)   //Photon does not hit a pixel
	   {
		  if(bDebug) cout<<"no pixel hit "<<fx<<"  "<<fy<<endl;
		  continue;
       }

      
      if(pixID>telData->iNumPixels)
         cout<<"Trying to access something that is out of range of the available pixel"<<endl;

      if(pixID>=0)
	    {                                     
           UInt_t lambda = (UInt_t)(v_f_lambda->at(p));

           //Lets see if the Photon will be detected after the Winston cone and the PMT
           //QE includes the winston cone and the cherenkov scaling factor, see function that
           //writes qe[]
           Float_t eff = qe[lambda]*telData->fRelQEwWC[pixID]/dEfficiencyFactor;
           if(bDebug)
 	   cout<<"eff: "<<eff<<" lambda "<<lambda<<"  qe[lambda] "<<qe[lambda]<<" telData->fRelQEwWC[pixID] "<<telData->fRelQEwWC[pixID]<<" efficiency factor: "<<dEfficiencyFactor<<endl;
           if(rand->Uniform()<eff)
                 {
	          AddPEToTrace(pixID, v_f_time->at(p)-(telData->fAveragePhotonArrivalTime-fStartSamplingBeforeAverageTime)); //Start filling fStartSamplingBeforeAverageTime ns before the average time
                  telData->iPEInPixel[pixID]++;
                 }
	    }
    }

  if(v_f_time->size()==0)
    delay = 0;

  if(bDebug)
    cout<<"Delay corrected arrival time:"<< telData->fAveragePhotonArrivalTime<<"+"<<delay<<"="<<telData->fAveragePhotonArrivalTime+delay<<endl;
  
  telData->fAveragePhotonArrivalTime+=delay;

  //cout<<"Finished loading the event"<<endl;

}


//-------------------------------------------------------------------------------------------
//
// Get simulated average single pe pulse. Returns a root histogram
TH1F* TraceGenerator::GetAverageSinglePEPulse()
{

  if(fAverageSinglePEPulse.empty())
    {
      cout<<"Pulse has not been set. Set a pulse shape first"<<endl;
      return NULL;
    }

  TH1F *h = new TH1F("h","Average Single PE Pulse",iNumSamplesAverageSinglePEPulse,-1*fStartTimeAveragePulse- fSamplingTimeAveragePulse*0.5,fStopTimeAveragePulse+ fSamplingTimeAveragePulse*0.5);
  h->GetXaxis()->SetTitle("Time [ns]");
  h->GetYaxis()->SetTitle("Amplitude normalized to peak");
  for(Int_t i = 0; i<iNumSamplesAverageSinglePEPulse; i++)
    {
      cout<<i<<"  "<<fAverageSinglePEPulse[i]<<endl;
      h->SetBinContent(i+1,fAverageSinglePEPulse[i]);
    }


  return h;
}

//-------------------------------------------------------------------------------------------
//
// Get simulated average single pe pulse. Returns a root histogram
TH1F* TraceGenerator::GetAverageSinglePELowGainPulse()
{

  if(fAverageLowGainSinglePEPulse.empty())
    {
      cout<<"Low Gain pulse has not been set. Set a pulse shape first"<<endl;
      return NULL;
    }
  
  TH1F *h = new TH1F("h","Average LowGain Single PE Pulse",iNumSamplesLowGainAverageSinglePEPulse,-1*fLowGainStartTimeAveragePulse- fSamplingTimeAveragePulse*0.5,fLowGainStopTimeAveragePulse+ fSamplingTimeAveragePulse*0.5);
  h->GetXaxis()->SetTitle("Time [ns]");
  h->GetYaxis()->SetTitle("Amplitude normalized to peak");
  for(Int_t i = 0; i<iNumSamplesLowGainAverageSinglePEPulse; i++)
    {
      h->SetBinContent(i+1,fAverageLowGainSinglePEPulse[i]);
    }


  return h;
}



//-------------------------------------------------------------------------------------------
//
// shows the average single pe pulse shape
void  TraceGenerator::ShowAverageSinglePEPulse()
{

  TH1F *h = GetAverageSinglePEPulse();
  h->GetXaxis()->SetTitle("Time [ns]");
  h->GetYaxis()->SetTitle("Amplitude [photoelectrons]");

  TH1F *hL = GetAverageSinglePELowGainPulse();
  hL->GetXaxis()->SetTitle("Time [ns]");
  hL->GetYaxis()->SetTitle("Amplitude [photoelectrons]");

  TCanvas *cSinglePEShape = new TCanvas("cSinglePEShape","The single pe pulse Shape use in the simulation",1000,700);
  
  cSinglePEShape->Divide(2,1);
  cSinglePEShape->Draw();
  cSinglePEShape->cd(1);
  gPad->SetGrid();
	   
  h->Draw();  
  cSinglePEShape->cd(2);
  gPad->SetGrid();
	   
  hL->Draw();  
	      
  TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
	      
  //
  // While reading the input process gui events asynchronously
  //
  timer.TurnOn();
  TString input = Getline("Type <return> to go on: ");
  timer.TurnOff(); 
	      
}

//-------------------------------------------------------------------------------------
//
// Function to reset the SiPM
void TraceGenerator::ResetSiPM() 
{

   if(bDebug)
	        cout<<"Resetting the SiPM"<<endl;

   vSiPMCellFired.clear();
   vector< Bool_t > bCells;
   for(int p = 0; p <iNumPixels; p++)
   {
	  if(bDebug)
		  cout<<"This SiPM "<<p<<" has "<<vNumCellsPerSIPM[p]<<" cells"<<endl;
      bCells.assign(vNumCellsPerSIPM[p],kFALSE);
      vSiPMCellFired.push_back(bCells);
   } 

}

//---------------------------------------------------------------------------------
//
//Function to simulate the proper dynamic behavior of a SiPM
//This includes optical crosstalk and the number of cells 
//available in one SiPM
Int_t TraceGenerator::SimulateSiPM(Int_t PixelID, Int_t NumPE)
{
    if(bDebug)
	{
	 cout<<"We are in camera pixel: "<<PixelID<<endl;   
     cout<<"simulating optical crosstalk of a SiPM"<<endl;
	 cout<<" Probability "<<vSiPMOpticalCrosstalk[PixelID]<<endl;
     }
     if(vNumCellsPerSIPM[PixelID]<=0)
       {
         cout<<"SiPM of Pixel "<<PixelID<<" has no cells, will skip the SiPM simulation....not sure if you want that!!! You might want to check your config file"<<endl;
         return NumPE;
       } 

	 int MeasuredPE = 0;
 
     //loop over all primary pe's
     //figure out for each pe if it fires a cell
     //and how many cells are fired in addition due
     //to optical crosstalk
	 for(int p = 0; p<NumPE; p++ )
	 {
       //in which cell of the SiPM is the pe created
	   int iCellID = (int)(rand->Uniform()*vNumCellsPerSIPM[PixelID]);

	   //if cell not fired, fire it and do optical crosstalk
	   if(!vSiPMCellFired[PixelID][iCellID])
		 {
           MeasuredPE++;
		   vSiPMCellFired[PixelID][iCellID]=kTRUE;
           if(bDebug)
			   cout<<"One SiPM cell fired"<<endl;
		   //Do optical crosstalk
		   while(1)
		   {

              if(rand->Uniform()>vSiPMOpticalCrosstalk[PixelID]) 
			   break;   //no additional crosstalk photon fires another cell
              
              //ok, another cell is supposed to fire, if that is possible
              //this is a simplification. In principle the distance between the 
              //cell from which the optical crosstalk photon originates and
              //the cell in which it converts is correlated. Here we assume
              //that each cell in the SiPM can be fired with the same probability
              //But this should be almost exactly right.   
              iCellID = (int)(rand->Uniform()*vNumCellsPerSIPM[PixelID]);

              if(bDebug)
				  cout<<"fired due to O-Xtalk "<<vSiPMCellFired[PixelID][iCellID]<<endl;

              if(vSiPMCellFired[PixelID][iCellID])
                break; //The cell was fired before -> no more crosstalk

              //ok a new cell of the SiPM is fired
              MeasuredPE++;
              vSiPMCellFired[PixelID][iCellID]=kTRUE;

              //start all over with the optical crosstalk
            } //end optical crosstalk loop
	      } //end treating this one primary pe

     } //End looping over all primary pe's


    if(bDebug)
	      cout<<"in total in this SiPM "<<MeasuredPE<<" pes fire"<<endl;


     return MeasuredPE; //update the number of measured pe's
}

//--------------------------------------------------------------------------------------------
//
//Adds one PE signal to the trace of pixel PixelID at time 
// time is in units of nanoseconds
void TraceGenerator::AddPEToTrace(Int_t PixelID, Float_t time, Int_t NumPE)
{

  if(bDebug)
   cout<<"Adding a "<<NumPE<<" pe signal to Pixel: "<<PixelID<<" at time "<<time<<endl;

  //Get the right amplitude for the simulation of an SiPM 
  if(bSiPM)
  {
    NumPE = SimulateSiPM(PixelID, NumPE);
    if(NumPE==0)
      return;
  }

  //Add some fluctuations to the amplitude
  Float_t newAmpl=0.;
  while(newAmpl<=0)
    newAmpl = rand->Gaus((Float_t)NumPE,sqrt((Float_t)NumPE)*fSigmaSinglePEPulseHeightDistribution);

  if(bDebug)
    cout<<"Amplitude of signal after PMT fluctuations "<<newAmpl<<endl;

  //Multiply with the relative gain of this pixel
  newAmpl *= telData->fRelGain[PixelID]; 

  if(bDebug)
    cout<<"Amplitude of signal times absolute pixel gain "<<newAmpl<<endl;
   
  //save the time and the amplitude for a later generation of the low gain trace
  telData->fTimesInPixel[PixelID].push_back(time);
  telData->fAmplitudesInPixel[PixelID].push_back(newAmpl);

 
  //In Reference of Trace
  Float_t StartTime = time-fStartTimeAveragePulse;
  Int_t StartSample = Int_t(StartTime/fSamplingTime+1);
  Int_t StopSample = Int_t((StartTime+fStartTimeAveragePulse+fStopTimeAveragePulse)/fSamplingTime+1);

  //The time between the start of the Trace and the first sampled value
  Float_t TimeAveragePulse= StartSample*fSamplingTime-StartTime; 

  //Catch the case when we get out of the trace window
  //Trace fully out of window
  if(StartTime<-fStartTimeAveragePulse-fStopTimeAveragePulse)
    return;

  if(StartSample<0)
    {
      TimeAveragePulse = -1.0*StartTime;
      StartTime = 0;
      StartSample = 0;
    }
  if(StopSample>telData->iNumSamplesPerTrace)
    {
      StopSample = telData->iNumSamplesPerTrace;
    }

  //  cout<<"Fill pe from: "<<StartSample<<" to "<<StopSample<<" Start time is "<<StartTime<<endl;

  //Fill the PE into the trace
  TimeAveragePulse = TimeAveragePulse / fSamplingTimeAveragePulse;
  Float_t step = fSamplingTime/fSamplingTimeAveragePulse;
  

  for(Int_t i=StartSample;i<StopSample;i++)
    {
      Int_t s = (Int_t)(TimeAveragePulse);
      telData->fTraceInPixel[PixelID][i]+=fAverageSinglePEPulse[s]*newAmpl;

      //loop over neighboring pixel if crosstalk
      if(bCrosstalk)
        {
          if(bDebug)
            cout<<"Adding Crosstalk "<<fCrosstalk<<" to neighboring pixel"<<endl;

          for(UInt_t n = 0;n<vNeighbors[PixelID].size(); n++)
            {
			if(bDebug)
			   cout<<"visiting pixel "<<vNeighbors[PixelID][n]<<endl;
               telData->fTraceInPixel[ vNeighbors[PixelID][n] ][i]+=fAverageSinglePEPulse[s]*newAmpl*fCrosstalk;
            }
        } 

      TimeAveragePulse+=step;
    }


  if(bDebug)
   {
     //debugDisplay->Show(telData->GetTelescopeID(),PixelID);
   }

}

//---------------------------------------------------------------------
//Add NSB to all traces
void TraceGenerator::GenerateNSB()
{
   if(bDebug)
    cout<<"Generating NSB"<<endl;

   if(rand == 0)
    {
      cout<<"FillPixelsWithNSB: You need to set the random number Generator first"<<endl;
      exit(1);
    }

   if(fNSBRatePerPixel<=0 && bUseNSB)
     {
      cout<<"FillPixelsNSB: You need to set the NSB rate in kHz per pixel with SetNSBRatePerPixel(Float_t rate)"<<endl;
      exit(1);
    }


   if(bSiPM)  //if we have SiPMs we need to reset them
     {
       ResetSiPM();
     }


  //cout<<iNumPixInSumGroup[0]<<endl;
  //Loop over all pixel
  if(bUseNSB)
  {
   for(Int_t i=0;i<iNumPixels;i++)
    {
      if(bDebug)
      cout<<"Pixel "<<i<<endl;
      //convert to counts per ns remember the rate is in units kHz 
      Float_t fNSBRate = fNSBRatePerPixel*1e-6 * telData->fRelQE[i]; ; 

      Float_t t = -1.0*fStartTimeAveragePulse;
      
      while(1)
	{
	  t+= -1*log(rand->Uniform())/fNSBRate;
	  
	 if(t>fTraceLength+fStartTimeAveragePulse)
	   break;

      int l = 1;
            
	  if(bAfterPulsing==kTRUE)
	  {
	    //mix in some afterpulsing
	    double m = rand->Uniform();

            //convert this into the proper afterpulsing amplitude if we 
            //are above 1.5 photoelectrons
            if( m < exp(fAPconstant + fAPslope * 1.5) )
              l = int( ( log(m) - fAPconstant ) / fAPslope +1 ) ;
        if(bDebug)
          cout<<"Afterpulsing added: "<<l<<endl;
	  }
                 
	  AddPEToTrace(i,t,l);
	}
    }
   } //end bUseNSB


   //Shift all traces down by the global (camera) mean, i.e. AC coupling
   telData->mean = 0.0;
    for(Int_t i=0;i<iNumPixels;i++)
    {
      Double_t SumTrace = 0.0;
      for(Int_t t=0;t<telData->iNumSamplesPerTrace;t++)
	{
	  SumTrace += telData->fTraceInPixel[i][t];
	}
      telData->mean+= SumTrace / (1.0*telData->iNumSamplesPerTrace * iNumPixels);
    }

    if(bDebug)
	cout<<"Pedestal "<<telData->mean<<endl;

    //shift trace up such that the mean is zero (AC coupling) 
    for(Int_t i=0;i<iNumPixels;i++)
      {
	for(Int_t t=0;t<telData->iNumSamplesPerTrace;t++)
	  {
	    telData->fTraceInPixel[i][t]=telData->fTraceInPixel[i][t]-telData->mean;
	  }
      }

  
   //copy all traces into the array of NSB traces only 

  //Copy vectors
  for(Int_t g=0;g<iNumPixels;g++)
    {
      telData->fTraceInPixelNSBOnly[g] = telData->fTraceInPixel[g]; 
    }

}



//--------------------------------------------------------------------------------------------------
//
void TraceGenerator::SetAfterPulsing(Bool_t afterpulsing, Float_t constant, Float_t slope)
{

 if(constant > 0 && afterpulsing==1)
    {
      cout<<"SetAFterPulsing: fAPconstant set to a value >0, bad!"<<endl;
      exit(1);
    }

  if(slope > 0 && afterpulsing==1)
    {
      cout<<"SetAFterPulsing: fAPslope set to a value >0, bad!"<<endl;
      exit(1);
    }

  fAPconstant=constant; 
  fAPslope=slope;
  bAfterPulsing=afterpulsing;   

}

//---------------------------------------------------------------------------------------------
//
// Set a Gaussian pulse for the single pe pulse shape. the argument is the FWHM auf the Gauss
// 2.35 times sigma in ns
void TraceGenerator::SetGaussianPulse(Float_t fwhm)
{ 


 if(fSamplingTimeAveragePulse <= 0)
    {
      cout<<"SetGaussianPulse: You need to set the sampling rate with SetSampleWidth(Float_t t)"<<endl;
      exit(0);
    }


 if(fwhm <= 0)
    {
      cout<<"SetGaussionPulse You need to set a FWHM > 0"<<endl;
      exit(0);
    }

 fStartTimeAveragePulse = 3*fwhm;
 fStopTimeAveragePulse = 3*fwhm;

 iNumSamplesAverageSinglePEPulse=Int_t((fStartTimeAveragePulse+fStopTimeAveragePulse)/fSamplingTimeAveragePulse)+1; 

  fAverageSinglePEPulse.assign(iNumSamplesAverageSinglePEPulse,0.0);

  //convert fwhm (in ns) into sampling units
  fwhm/=fSamplingTimeAveragePulse;
  Float_t sigma=fwhm/(2*sqrt(-2*log(0.5)));
  cout<<"Simulate a Gaussian with sigma "<<sigma*fSamplingTimeAveragePulse<<" ns"<<endl;

  for(Int_t i = 0; i<iNumSamplesAverageSinglePEPulse; i++)
    {
      Float_t t = -1*Int_t(fStartTimeAveragePulse/fSamplingTimeAveragePulse)+i;
      fAverageSinglePEPulse[i] = -1.0*exp(-1*t*t/(2*sigma*sigma));
      // cout<<i<<" "<<fAverageSinglePEPulse[i]<<endl;
    }

fAreaToPeakConversion = 1./(sqrt(2*TMath::Pi())*sigma) ;
cout<<"The area to peak conversion is: "<<fAreaToPeakConversion<<endl;

//copy everything to the low gain pulse shape
 fAverageLowGainSinglePEPulse=fAverageSinglePEPulse;
 fLowGainAreaToPeakConversion = fAreaToPeakConversion ;
 
 fLowGainStartTimeAveragePulse= fStartTimeAveragePulse;
 fLowGainStopTimeAveragePulse= fStopTimeAveragePulse;

 iNumSamplesLowGainAverageSinglePEPulse= iNumSamplesAverageSinglePEPulse;

}



//////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the shape of a single pe pulse from an external file
// the format of the file is plain text with two columns
// the first column holds the time in nanoseconds
// the second column holds the amplitude of the file in arbitrary units
// this function converts the pulse shape such that the peak amplitude is -1
void   TraceGenerator::SetSinglePEShapeFromFile(TString sfilename)
{ 
  
  if(fSamplingTimeAveragePulse <= 0)
    {
      cout<<"SetSinglePEShapeFromFile: You need to set the sampling rate with SetSampleWidth(Float_t t)"<<endl;
      exit(0);
    }

  std::ifstream PulseShapefile(sfilename.Data());
 
  if(!PulseShapefile)
    {
      cout<<"SetSinglePEShapeFromFile:could not open file with the sample Pulse Shape"<<endl;
      exit(1);
    }

  cout<<"Setting the high gain pulse Shape from file"<<endl;
  
  vector<Double_t> ts;
  vector<Double_t> y;
  



  while(PulseShapefile.good()) {
    Double_t time = 0;
    Double_t amplitude = 0;

    PulseShapefile >> time;
    PulseShapefile >> amplitude;

    // cout<<time<<"  "<<amplitude<<endl;

    ts.push_back(time);
    y.push_back(amplitude);
  }


 
  //Find the time where the pulse shape is minimal and get the conversion factor
  Int_t imin = 0;
  Double_t dmin = 1e9;
  Double_t integral = 0;
  for(UInt_t i = 1;i<ts.size()-1;i++)
    {
      integral+=y[i]*0.5*(ts[i+1]-ts[i-1]);;
      if(y[i]<dmin)
	{
	  imin=i;
	  dmin=y[i];
	}
    }

  fAreaToPeakConversion = dmin/(integral) ;
  cout<<"The area to peak conversion factor is: "<<fAreaToPeakConversion<<endl;

  dmin=fabs(dmin);

  cout<<"Found that the amplitude is minimal at time: "<<imin<<" =  "<<ts[imin]<<". The amplitude is: "<<y[imin]<<endl;
 
  cout<<"Number of samples in the sample pulse: "<<ts.size()<<endl;
  fStartTimeAveragePulse=fabs(ts[imin]-ts[0]);
  fStopTimeAveragePulse=ts[ts.size()-2]-ts[imin];

  cout<<"Start of sample pulse shape: "<<fStartTimeAveragePulse<<endl;
  cout<<"Stop  of sample pulse shape: "<<fStopTimeAveragePulse<<endl;

  //now fill the vector with the normalized sample pulse shape
  iNumSamplesAverageSinglePEPulse=Int_t((fStartTimeAveragePulse+fStopTimeAveragePulse)/fSamplingTimeAveragePulse); 
  
  fAverageSinglePEPulse.assign(iNumSamplesAverageSinglePEPulse,0.0);

//cout<<iNumSamplesAverageSinglePEPulse<<endl; 
  for(Int_t i = 0; i<iNumSamplesAverageSinglePEPulse; i++)
    {
      Float_t t = -1*Int_t(fStartTimeAveragePulse/fSamplingTimeAveragePulse)+i;
      t*=fSamplingTimeAveragePulse;

      unsigned d=0;
      while(ts[d]-ts[imin]<t+1e-6 && d<ts.size()-1)
	d++;

      //average between two samples
      Double_t amplitude = y[d-1] + (t- (ts[d-1]-ts[imin]) ) * (y[d]-y[d-1])/(ts[d]-ts[d-1]);
      //normalizing so the peak is one
      amplitude = amplitude / dmin ; 

      fAverageSinglePEPulse[i] = amplitude;
      //cout<<i<<" "<<t<<"  "<<d<<"  "<<fAverageSinglePEPulse[i]<<endl;
    }

 
  cout<<"Filled in the sample pulse shape"<<endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the shape of a low gain single pe pulse from an external file
// the format of the file is plain text with two columns
// the first column holds the time in nanoseconds
// the second column holds the amplitude of the file in arbitrary units
// this function converts the pulse shape such that the peak amplitude is -1
void   TraceGenerator::SetLowGainSinglePEShapeFromFile(TString sfilename)
{ 
  
  cout<<"Reading in the low gain single pe pulse shape"<<endl;

  if(fSamplingTimeAveragePulse <= 0)
    {
      cout<<"SetLowGainSinglePEShapeFromFile: You need to set the sampling rate with SetSampleWidth(Float_t t)"<<endl;
      exit(0);
    }

  std::ifstream PulseShapefile(sfilename.Data());
 
  if(!PulseShapefile)
    {
      cout<<"SetSinglePEShapeFromFile:could not open file with the low gain sample Pulse Shape: "<<sfilename.Data()<<endl;
      exit(1);
    }

  cout<<"Setting the high gain pulse Shape from file"<<endl;

  
  vector<Double_t> ts;
  vector<Double_t> y;
  
  while(PulseShapefile.good()) {
    Double_t time = 0;
    Double_t amplitude = 0;

    PulseShapefile >> time;
    PulseShapefile >> amplitude;

    // cout<<time<<"  "<<amplitude<<endl;

    ts.push_back(time);
    y.push_back(amplitude);
  }


  //Find the time where the pulse shape is minimal and get the conversion factor
  Int_t imin = 0;
  Double_t dmin = 1e9;
  Double_t integral = 0;
  for(UInt_t i = 1;i<ts.size()-1;i++)
    {
      integral+=y[i]*0.5*(ts[i+1]-ts[i-1]);
      if(y[i]<dmin)
	{
	  imin=i;
	  dmin=y[i];
	}
    }

  fLowGainAreaToPeakConversion = dmin/(integral) ;
  cout<<"The area to peak conversion factor is: "<<fLowGainAreaToPeakConversion<<endl;

  dmin=fabs(dmin);

  cout<<"Found that the amplitude is minimal at time: "<<imin<<" =  "<<ts[imin]<<". The amplitude is: "<<y[imin]<<endl;
 
  cout<<ts[0]<<endl;
  fLowGainStartTimeAveragePulse=fabs(ts[imin]-ts[0]);
  fLowGainStopTimeAveragePulse=ts[ts.size()-2]-ts[imin];

  cout<<"Start of sample pulse shape: "<<fLowGainStartTimeAveragePulse<<endl;
  cout<<"Stop  of sample pulse shape: "<<fLowGainStopTimeAveragePulse<<endl;

  //now fill the vector with the normalized sample pulse shape
  iNumSamplesLowGainAverageSinglePEPulse=Int_t((fLowGainStartTimeAveragePulse+fLowGainStopTimeAveragePulse)/fSamplingTimeAveragePulse); 
  
  fAverageLowGainSinglePEPulse.assign(iNumSamplesLowGainAverageSinglePEPulse,0.0);
  cout<<"done assigning"<<endl;
 
  for(Int_t i = 0; i<iNumSamplesLowGainAverageSinglePEPulse; i++)
    {
      Float_t t = -1*Int_t(fLowGainStartTimeAveragePulse/fSamplingTimeAveragePulse)+i;
      t*=fSamplingTimeAveragePulse;

      unsigned d=0;
      while(ts[d]-ts[imin]<t+1e-6 && d<ts.size()-1)
	d++;

      //average between two samples
      Double_t amplitude = y[d-1] + (t- (ts[d-1]-ts[imin]) ) * (y[d]-y[d-1])/(ts[d]-ts[d-1]);

      //normalizing so the peak is one
      amplitude = amplitude / dmin ; 


      fAverageLowGainSinglePEPulse[i] = amplitude;
      // cout<<i<<" "<<t<<"  "<<fAverageLowGainSinglePEPulse[i]<<endl;
    }

 
  cout<<"Filled in the low gain sample pulse shape"<<endl;
}




//----------------------------------------------------------------------------------------
//
//and the sampling stepwidth
void TraceGenerator::SetTraceSampleWidthAndLength(Float_t t,Float_t length,Float_t offset)
{

  if(t<fSamplingTimeAveragePulse)
    {
      cout<<"Sampling time has to be larger than the sampling time of the single pe pulse "<<fSamplingTimeAveragePulse<<endl;
      exit(0);
    }

  if(length<0)
    {
      cout<<"You try to set the length of the simulated trace below zero "<<length<<endl;
      exit(0);
    }

  if(length<offset)
    {
      cout<<"You want to simulate a trace that is less than "<<offset<<"ns long. This is less then the time the trace gots sampled before the average photon arrival time. I suggest you do at least 100ns. Right now you want to simulated "<<length<<" ns"<<endl;
      exit(0);
    }

 fSamplingTime = t;

 fTraceLength=length;

 fStartSamplingBeforeAverageTime=offset;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////


vector<Float_t> TraceGenerator::GetLowGainTrace(Int_t PixelID){

  vector<Float_t> trace(telData->iNumSamplesPerTrace,-1*telData->mean);


  for(UInt_t g=0;g<telData->fTimesInPixel[PixelID].size();g++)
    {
      Float_t amplitude = telData->fAmplitudesInPixel[PixelID][g];
      Float_t time = telData->fTimesInPixel[PixelID][g];
      //In Reference of Trace
      Float_t StartTime = time-fLowGainStartTimeAveragePulse;
      Int_t StartSample = Int_t(StartTime/fSamplingTime+1);
      Int_t StopSample = Int_t((StartTime+fLowGainStartTimeAveragePulse+fLowGainStopTimeAveragePulse)/fSamplingTime+1);

      //The time between the start of the Trace and the first sampled value
      Float_t TimeAveragePulse= StartSample*fSamplingTime-StartTime;

      //Catch the case when we get out of the trace window
      //Trace fully out of window
      if(StartTime<-fLowGainStartTimeAveragePulse-fLowGainStopTimeAveragePulse)
    continue;

      if(StartSample<0)
    {
      TimeAveragePulse = -1.0*StartTime;
      StartTime = 0;
      StartSample = 0;
    }
      if(StopSample>telData->iNumSamplesPerTrace)
    {
      StopSample = telData->iNumSamplesPerTrace;
    }

      // cout<<"Fill pe from: "<<StartSample<<" to "<<StopSample<<" Start time is "<<StartTime<<endl;


      //Fill the PE into the trace
      TimeAveragePulse = TimeAveragePulse / fSamplingTimeAveragePulse;
      Float_t step = fSamplingTime/fSamplingTimeAveragePulse;
      for(Int_t i=StartSample;i<StopSample;i++)
    {
          Int_t s = (Int_t)(TimeAveragePulse);
      trace[i]+= fAverageLowGainSinglePEPulse[s]*amplitude;
      TimeAveragePulse+=step;
    }

    }

  return trace;

}


//--------------------------------------------------------------------------------------------------------
// Set the sampling width of the single
void   TraceGenerator::SetSinglePESamlingWidth(Float_t width){


  if(width<0)
    {
      cout<<"You try to set the sampling time of the single PE pulse shape below zero "<<width<<endl;
      exit(0);
    }

 fSamplingTimeAveragePulse = width;

}



//--------------------------------------------------------------------------------------------
//
// Set the sigma of the pulse height distribution of the single PEs
//
void TraceGenerator::SetSigmaofSinglePEPulseHeightDistribution(Float_t sigma)
{ 
 if(sigma <= 0)
    {
      cout<<"SetSigmaofSinglePEPulseHeightDistribution: You need to set the sigma of the singple pe's with SetSigmaofSinglePEPulseHeightDistribution(Float_t sigma)"<<endl;
      exit(0);
    }


 fSigmaSinglePEPulseHeightDistribution=sigma;

}

//----------------------------------------------------------------------------------------
//  Prints how often the trace was too short
void   TraceGenerator::PrintHowOftenTheTraceWasTooShort(){

  cout<<"The trace was too short on the low end: "<<NumTraceIsTooShortLowEnd<<" times"<<endl;
  cout<<"The trace was too short on the high end: "<<NumTraceIsTooShortHighEnd<<" times"<<endl;

}



//----------------------------------------------------------------------------------------
//Reads in  the config file and sets all variables
void   TraceGenerator::SetParametersFromConfigFile(ReadConfig *readConfig){

   cout << "TraceGenerator::SetParametersFromConfigFile " << endl;
   cout<<  "--------------------------------------------"<<endl;


  SetAfterPulsing(readConfig->GetAfterPulsingUsage(iTelType), 
               readConfig->GetAfterPulsingConstant(iTelType), readConfig->GetAfterPulsingSlope(iTelType));


  SetSinglePESamlingWidth(readConfig->GetSampleWidthAveragePulse(iTelType));
  SetSigmaofSinglePEPulseHeightDistribution(readConfig->GetSigmaofSinglePEPulseHeightDistribution(iTelType));
 

  SetTraceSampleWidthAndLength(readConfig->GetSamplingTime(iTelType),readConfig->GetTraceLength(iTelType),readConfig->GetStartSamplingTimeOffsetFromAveragePhotonTime(iTelType));
                                            

  Float_t fFWHMofSinglePEPulse = readConfig->GetFWHMofSinglePEPulse(iTelType);

  if(fFWHMofSinglePEPulse<=0)
	{
	   SetSinglePEShapeFromFile(readConfig->GetNameofSinglePEPulseFile(iTelType));
	   cout<<"The single pe pulse shape was read from: "<<readConfig->GetNameofSinglePEPulseFile(iTelType)<<endl;

	   SetLowGainSinglePEShapeFromFile(readConfig->GetNameofLowGainSinglePEPulseFile(iTelType));
	   cout<<"The low gain single pe pulse shape was read from: "<<readConfig->GetNameofLowGainSinglePEPulseFile(iTelType)<<endl;
	}
	else
	{
          SetGaussianPulse(fFWHMofSinglePEPulse);
	}
   if(fAverageSinglePEPulse.size()==0)
	  {
	     cout<<"Something went wrong, you did neither provide an external file with the pulse shape nor did you set a fwhm that can be used to generate a sample Gaus pulse"<<endl;
	     exit(1);
	  }

   //Set the Quantum efficienicy
   vector<Float_t> wlOrig = readConfig->GetWavelengthsOfQEValues(iTelType);
   vector<Float_t> qeOrig = readConfig->GetQEValues(iTelType);
   UInt_t maxwl = (UInt_t)wlOrig[wlOrig.size()-1];
   //Set the qe values in the array
   qe.assign(maxwl,0);
   cout<<"Setting the QE values"<<endl;
   for(UInt_t w = wlOrig[0]; w<maxwl; w++)
     {
        Int_t i= 0;
        while(wlOrig[i]<w)  i++;
           i--;

        Float_t Qeff = qeOrig[i]+ (qeOrig[i+1]-qeOrig[i])/(wlOrig[i+1]-wlOrig[i])*(w-wlOrig[i]);

        qe[w] = Qeff;

        //cout<<w<<"  "<<qe[w]<<endl;
     }


  //NSB
  fNSBRatePerPixel = readConfig->GetNSBRate(iTelType);      //the NSB rate kHz per mm squared in the focal plane;
  bUseNSB = readConfig->GetNSBUsage();

  iNumPixels = readConfig->GetNumberPixels(iTelType); //Has to be filled with telescope type
  vNeighbors = readConfig->GetNeighbors(iTelType);


  //SiPM
  bSiPM = readConfig->GetSiPMUsage(iTelType);
  vNumCellsPerSIPM = readConfig->GetNumCellsPerSiPM(iTelType);                  //the number of cells in one SiPM
  vSiPMOpticalCrosstalk = readConfig->GetSiPMOpticalCrosstalk(iTelType);          // 


  //Crosstalk between pixel
  bCrosstalk = readConfig->GetCrosstalkUsage(iTelType);                    //Use crosstalk between pixel
  fCrosstalk = readConfig->GetCrosstalkValue(iTelType);                    //Crosstalk value


  //for the Gridsearch algorithm
  fXTubeMM =   readConfig->GetXMM(iTelType);
  fYTubeMM =   readConfig->GetYMM(iTelType);
  fSizeTubeMM =readConfig->GetTubeSizeMM(iTelType);
  iTubeSides = readConfig->GetTubeSides(iTelType);
  fRotAngle =  readConfig->GetTubeRotAngle(iTelType);

}
