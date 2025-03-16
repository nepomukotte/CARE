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
  fTraceLength = 100; //Length of simulated trace in ns
  fStartSamplingBeforeAverageTime=30;
  bAfterPulsing = kFALSE;
  fFWHMofSinglePEPulse = -1;
  fLinearGainmVPerPE = -1;
  fPileUpWindow = -1;
  //Read the config file
  SetParametersFromConfigFile( readConfig );
  
  rand = generator;
  string s = "";                                                                             

  gridsearch = new GOrderedGridSearch(fXTubeMM,fYTubeMM,fSizeTubeMM,iTubeSides,fRotAngle,19,19,1,s);
  if(bDebug) cout<<"Done initializing the TraceGenerator"<<endl;

  if(bDebug)
   { 
     for(UInt_t i=0;i<max(fLowGainPulse.size(),fHighGainPulse.size());i++){
      ShowPulseShape(i);
     }
   }
}


///////////////////////////////////////////////////
//
//  Setting the TelescopeData container

void   TraceGenerator::SetTelData(TelescopeData *TelData){

telData = TelData;


}


//-----------------------------------------------------------------------------------------------------------------------
//Loads the PEs into the pixels. If UseNSB is set true NSB is added to the trace
void  TraceGenerator::LoadCherenkovPhotons(std::vector< float > *v_f_X,std::vector< float > *v_f_Y,std::vector< float > *v_f_time,std::vector< float > *v_f_lambda, Float_t delay, Double_t dEfficiencyFactor)
{
  

  if(iNumPixels == 0)
    {
      cout<<"You need to set the pixel coordinates first, maybe they are empty"<<endl;
      exit(1);
    }

  if(fSamplingTime <= 0)
    {
      cout<<"You need to set the sampling rate with SetTraceSampleWidthAndLength(Float_t t)"<<endl;
      exit(0);
    }

  if(bDebug)
   {
    cout<<"Camera has "<<iNumPixels<<" pixel"<<endl;
   }

    //Load the NSB into the Traces, will be skipped in function if no NSB generation is wanted
    GenerateNSB();
 
  if(bDebug)
   cout<<"NSB done"<<endl;

  //Find the average time of all PEs
  telData->fAveragePhotonArrivalTime = 0;
  Float_t fMinPhotonArrivalTime = 1e8;
  Float_t fMaxPhotonArrivalTime = -1e8;
  for(UInt_t p=0; p<v_f_time->size(); p++)
    {
	  //add transit time spread
	  v_f_time->at(p)+=rand->Gaus(0.0,telData->fTransitTimeSpread);
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
           Float_t eff = qe[iTubeType[pixID]][lambda]*telData->fRelQEwWC[pixID]/dEfficiencyFactor;
           if(bDebug)
 	   cout<<"eff: "<<eff<<" lambda "<<lambda<<" pixel "<<pixID<<"  qe[iTubeType[pixID][lambda] "<<qe[iTubeType[pixID]][lambda]<<" telData->fRelQEwWC[pixID] "<<telData->fRelQEwWC[pixID]<<" efficiency factor: "<<dEfficiencyFactor<<endl;
           if(rand->Uniform()<eff)
                 {
	          AddPEToTrace(pixID, v_f_time->at(p)-(telData->fAveragePhotonArrivalTime-fStartSamplingBeforeAverageTime)); //Start filling fStartSamplingBeforeAverageTime ns before the average time
                  telData->iPEInPixel[pixID]++;
                  telData->fSumTimeInPixel[pixID]+=v_f_time->at(p)-telData->fAveragePhotonArrivalTime;
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

  //loop over neighboring pixel if crosstalk
  if(bCrosstalk)
    {
     if(bDebug==1)
      {
        cout<<"Adding Crosstalk "<<fCrosstalk<<" to neighboring pixel"<<endl;
      }
        for(UInt_t n = 0;n<vNeighbors[PixelID].size(); n++)
            {
		if(bDebug)
		   cout<<"visiting pixel "<<vNeighbors[PixelID][n]<<endl;

                telData->fTimesInPixel[ vNeighbors[PixelID][n] ].push_back(time);
                telData->fAmplitudesInPixel[ vNeighbors[PixelID][n] ].push_back(newAmpl*fCrosstalk);
            }
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

   if(vNSBRatePerPixel[0]<=0 && bUseNSB) //this only checks correctly if there is one PMT/SiPM type in the telescope
     {
      cout<<"FillPixelsNSB: You need to set the NSB rate in kHz per pixel"<<endl;
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
      //if(bDebug)
	//   cout<<"Pixel "<<i<<" NSB "<<vNSBRatePerPixel[i]<<endl;
      //convert to counts per ns remember the rate is in units kHz 
      
      Float_t fNSBRate = vNSBRatePerPixel[i]*1e-6 * telData->fRelQE[i];
      if(bDebug)
           cout<<"Pixel "<<i<<" PMT type: "<<iTubeType[i]<<" NSB Rate "<< vNSBRatePerPixel[i]*1e-6 * telData->fRelQE[i] <<endl;
      Float_t t = -1.0*fHighGainStopTime[0];
      
      while(1)
	{
	  t+= -1*log(rand->Uniform())/fNSBRate;
	  
	 if(t>fTraceLength+fHighGainStartTime[0])
	   break;

      int l = 1;
            
	  if(bAfterPulsing==kTRUE && vAPslope[iTubeType[i]]!=0)
	  {
	    //mix in some afterpulsing
	    double m = rand->Uniform();

            //convert this into the proper afterpulsing amplitude if we 
            //are above 1.5 photoelectrons
            if( m < exp(vAPconstant[iTubeType[i]] + vAPslope[iTubeType[i]] * 1.5) )
              l = int( ( log(m) - vAPconstant[iTubeType[i]] ) / vAPslope[iTubeType[i]] +1 ) ;
        if(bDebug)
          cout<<"Afterpulsing added: "<<l<<endl;
	  }
                 
	  AddPEToTrace(i,t,l);
	}
    }
   } //end bUseNSB

}



//--------------------------------------------------------------------------------------------------
//
void TraceGenerator::SetAfterPulsing(Bool_t afterpulsing, vector<Float_t> constant, vector<Float_t> slope)
{

 for(uint i = 0; i<constant.size();i++)
  {
    if(constant[i] > 0 && afterpulsing==1)
      {
        cout<<"SetAFterPulsing: vAPconstant "<<i<<" set to a value >0, bad!"<<endl;
        exit(1);
      }

    if(slope[i] > 0 && afterpulsing==1)
      {
        cout<<"SetAFterPulsing: vAPslope "<<i<<" set to a value >0, bad!"<<endl;
        exit(1);
      }
  }
  vAPconstant=constant; 
  vAPslope=slope;
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
      cout<<"SetGaussianPulse: You need to set the sampling rate with SetTraceSampleWidthAndLength"<<endl;
      exit(0);
    }


 if(fwhm <= 0)
    {
      cout<<"SetGaussionPulse You need to set a FWHM > 0"<<endl;
      exit(0);
    }
         
  Float_t fStartTime = 3*fwhm;
  Float_t fStopTime = 3*fwhm;
  

  Int_t iNumSamples=Int_t((fStartTime+fStopTime)/fSamplingTimeAveragePulse)+1; 

  vector< Float_t > vPulse;
  vPulse.assign(iNumSamples,0.0);

  //convert fwhm (in ns) into sampling units
  fwhm/=fSamplingTimeAveragePulse;
  Float_t sigma=fwhm/(2*sqrt(-2*log(0.5)));
  cout<<"Simulate a Gaussian pulse shape with sigma "<<sigma*fSamplingTimeAveragePulse<<" ns"<<endl;

  for(Int_t i = 0; i<iNumSamples; i++)
    {
      Float_t t = -1*Int_t(fStartTime/fSamplingTimeAveragePulse)+i;
      vPulse[i] = -1.0*exp(-1*t*t/(2*sigma*sigma));
    }

  Float_t fAreaToPeakConversion = 1./(sqrt(2*TMath::Pi())*sigma) ;
  cout<<"The area to peak conversion is: "<<fAreaToPeakConversion<<endl;

  //copy everything where it belongs
  fHighGainAreaToPeakConversion.push_back(fAreaToPeakConversion);
  fHighGainPulse.push_back(vPulse);
  fHighGainLinearAmplitude.push_back(1.0);
  fHighGainNonLinearityFactor.push_back(1.0);
  fHighGainStartTime.push_back(fStartTime);
  fHighGainStopTime.push_back(fStopTime);

  fLowGainAreaToPeakConversion.push_back(fAreaToPeakConversion);
  fLowGainPulse.push_back(vPulse);
  fLowGainLinearAmplitude.push_back(1.0);
  fLowGainNonLinearityFactor.push_back(1.0);
  fLowGainStartTime.push_back(fStartTime);
  fLowGainStopTime.push_back(fStopTime);

}

//////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the pulse shapes and the non linearities with an external file
// the format of the file is plain text with two columns
// the very first line gives the amplitude in mV for which the pulse shape is valid and the 
// relative deviation from linearity measured/expected
// the first column holds the time in nanoseconds
// the second column holds the amplitude of the file in arbitrary units
// this function converts the pulse shape such that the peak amplitude is -1
void   TraceGenerator::SetPulseShapesFromFile(TString sfilename,Bool_t bLowGain=kFALSE)
{ 
  
  cout<<"Reading in the pulse shapes; LowGain="<<bLowGain<<endl;

  if(fSamplingTimeAveragePulse <= 0)
    {
      cout<<"SetPulseShapesFromFile: You need to set the sampling rate with SetTraceSampleWidthAndLength"<<endl;
      exit(0);
    }

  std::ifstream PulseShapefile(sfilename.Data());
 
  if(!PulseShapefile)
    {
      cout<<"SetPulseShapesFromFile: I could not open file with the Pulse Shapes: "<<sfilename.Data()<<endl;
      exit(1);
    }


  //looping over the file getting one pulse shape after the other
  Int_t iPulseNumber=0;
  string iline;
  getline( PulseShapefile, iline ) ;
  while(PulseShapefile.good()) {

    vector<Double_t> ts;
    vector<Double_t> y;

    Double_t time = 0;
    Double_t amplitude = 0;

  
    //the first number in the first line gives the amplitude in mV 
    //for the pulse shape if the system would be linear.
    //The second number gives the factor by which the true amplitude deviates from linearity
    istringstream i_stream( iline );
    Float_t fLinearAmplitude;
    Float_t fNonLinearityFactor;
    i_stream>>fLinearAmplitude;
    i_stream>>fNonLinearityFactor; 

    cout<<"Pulse Number "<<iPulseNumber<<", amplitude if linear: "<<fLinearAmplitude<<" [mV], "<<" non-linearity: "<<fNonLinearityFactor<<endl;

    //loop over rest of lines to read in the pulse. The end of the pulse shape is marked with a *
    getline( PulseShapefile, iline );
    while( (iline.substr( 0, 1 ) != "*") && PulseShapefile.good() ){

      istringstream(iline)>>time>>amplitude;
      cout<<time<<"  "<<amplitude<<endl;

      ts.push_back(time);
      y.push_back(amplitude);
      getline( PulseShapefile, iline );
    }

    //process the pulse shape by normalizing it to the minimum amplitude (assuming negative pulse shape) and filling it into the telescope data container
    Int_t imin = 0;
    Int_t itzero = 0;
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
        if(ts[i+1]>0 && ts[i-1]<0)
	 {
	  itzero=i;
	}

    }

  Float_t fAreaToPeakConversion = dmin/(integral) ;
  cout<<"The area to peak conversion factor is: "<<fAreaToPeakConversion<<endl;

  dmin=fabs(dmin);

  cout<<"Found that the amplitude is minimal at time: "<<imin<<" =  "<<ts[imin]<<". The amplitude is: "<<y[imin]<<endl;
 
  cout<<"Number of samples in the sample pulse: "<<ts.size()<<endl;
  Float_t fStartTime=fabs(ts[itzero]-ts[0]);
  Float_t fStopTime=ts[ts.size()-2]-ts[itzero];

  cout<<"Start of sample pulse shape: "<<fStartTime<<endl;
  cout<<"Stop  of sample pulse shape: "<<fStopTime<<endl;

  //now fill the vector with the normalized sample pulse shape
  Int_t iNumSamples=Int_t((fStartTime+fStopTime)/fSamplingTimeAveragePulse);
  vector<Float_t> fPulse;  
  fPulse.assign(iNumSamples,0.0);

  for(Int_t i = 0; i<iNumSamples; i++)
    {
      Float_t t = -1*Int_t(fStartTime/fSamplingTimeAveragePulse)+i;
      t*=fSamplingTimeAveragePulse;

      //find the sample of the raw pulse shape that is the next one after time t
      unsigned d=0;
      while(ts[d]-ts[itzero]<t+1e-6 && d<ts.size()-1)
	d++;

      //average between two samples
      Double_t amplitude = y[d-1] + (t- (ts[d-1]-ts[itzero]) ) * (y[d]-y[d-1])/(ts[d]-ts[d-1]);
      //normalizing so the peak is one
      amplitude = amplitude / dmin ; 

      fPulse[i] = amplitude;
      //cout<<i<<" "<<t<<"  "<<d<<"  "<<fAverageSinglePEPulse[i]<<endl;
    }


    if(bLowGain)
      {
         fLowGainPulse.push_back(fPulse);
         fLowGainAreaToPeakConversion.push_back(fAreaToPeakConversion);
         fLowGainStartTime.push_back(fStartTime);
         fLowGainStopTime.push_back(fStopTime);
         fLowGainLinearAmplitude.push_back(fLinearAmplitude);
         fLowGainNonLinearityFactor.push_back(fNonLinearityFactor);

      }
    else
      {

         //if we are the first pulse lets get the gain conversion factor
         fHighGainPulse.push_back(fPulse);
         fHighGainAreaToPeakConversion.push_back(fAreaToPeakConversion);
         fHighGainStartTime.push_back(fStartTime);
         fHighGainStopTime.push_back(fStopTime);
         fHighGainLinearAmplitude.push_back(fLinearAmplitude);
         fHighGainNonLinearityFactor.push_back(fNonLinearityFactor);
      }

  if(bDebug)
   {
     //ShowPulseShape(iPulseNumber);
   }


    //increment the number of pulses by one
    iPulseNumber++;
  getline( PulseShapefile, iline ) ;
  }//end with the while loop, going to the next pulse shape in the file

 //make consistency check that all pulse shapes are in increasing order.
 vector< Float_t> *la = &fHighGainLinearAmplitude;
 if(bLowGain)
   la = &fLowGainLinearAmplitude;
 if(la->size()>1)
  {
    for(UInt_t i=1;i<la->size();i++)
      {
        if(la->at(i)<la->at(i-1))
           {
             cout<<"Pulse shape "<<i<<" is for a smaller amplitude ("<<la->at(i)<<") than the previous one ("<<la->at(i-1)<<")"<<endl;     
             exit(0);
           } 
      }
  }
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
//
// Assemble the trace for a given pixel in low or high gain
//

void TraceGenerator::BuildTrace(Int_t PixelID,Bool_t bLowGain){

  //First Check if the fPileUpAmplitudeForPhoton array is zero. If it is we have to generate it first
  if(telData->fPileUpAmplitudeForPhoton[PixelID].size()==0 && telData->fTimesInPixel[PixelID].size()!=0)
    {
       telData->fPileUpAmplitudeForPhoton[PixelID].assign(telData->fTimesInPixel[PixelID].size(),0);
       //loop over all the photons in that pixel and determine for each photon how many photons are next to it
       //within a time window of +- fPileUpWindow
        for(UInt_t g=0;g<telData->fTimesInPixel[PixelID].size();g++)
         { 
           //start with the photons own amplitude
           telData->fPileUpAmplitudeForPhoton[PixelID][g]+=telData->fAmplitudesInPixel[PixelID][g];

           for(UInt_t l=g+1;l<telData->fTimesInPixel[PixelID].size();l++)          
              { 
                 Float_t fDeltaT = telData->fTimesInPixel[PixelID][g]-telData->fTimesInPixel[PixelID][l];

                 if(fDeltaT<0)
                    fDeltaT*=-1.0;
 
                 if(fDeltaT<fPileUpWindow)
                     {
                        telData->fPileUpAmplitudeForPhoton[PixelID][g]+=telData->fAmplitudesInPixel[PixelID][l];
                        telData->fPileUpAmplitudeForPhoton[PixelID][l]+=telData->fAmplitudesInPixel[PixelID][g];
                     }
              }
         }
    }//end building the fPileUpAmplitudeForPhoton array

   //add electronic noise to the high gain trace
   if(telData->fSigmaElectronicNoise>0 && !bLowGain)
     {
       for(Int_t i=0; i<telData->iNumSamplesPerTrace;i++)
         {
               telData->fTraceInPixel[PixelID][i]=rand->Gaus(0.0,telData->fSigmaElectronicNoise);
         }
     }
   else
     telData->fTraceInPixel[PixelID].assign(telData->iNumSamplesPerTrace,0.0);


  //Now prepare to fill in the photons in the trace

  //get vectors to low gain or high gain pulses and information respectively
  vector< vector<Float_t> > *vPulse; 
  //vector<Float_t> *vAreaToPeakConversion;
  vector<Float_t> *vStartTime;
  vector<Float_t> *vStopTime;
  vector<Float_t> *vLinearAmplitude;
  vector<Float_t> *vNonLinearityFactor;

  if(bLowGain)
   {
    vPulse = &fLowGainPulse; 
    //vAreaToPeakConversion = &fLowGainAreaToPeakConversion;
    vStartTime = &fLowGainStartTime;
    vStopTime = &fLowGainStopTime;
    vLinearAmplitude = &fLowGainLinearAmplitude;
    vNonLinearityFactor = &fLowGainNonLinearityFactor;
   }
  else
   {
    vPulse = &fHighGainPulse;    
    //vAreaToPeakConversion = &fHighGainAreaToPeakConversion;
    vStartTime = &fHighGainStartTime;
    vStopTime = &fHighGainStopTime;
    vLinearAmplitude = &fHighGainLinearAmplitude;
    vNonLinearityFactor = &fHighGainNonLinearityFactor;
   }

  //loop over all the photons in the trace
  for(UInt_t g=0;g<telData->fTimesInPixel[PixelID].size();g++)
    {
      Float_t time = telData->fTimesInPixel[PixelID][g];

      Float_t fNonLinearity = 1;
      
      //Figure out which pulse shape we ought to use
      UInt_t uPulseShape = 0;
      
      //if we only have one pulse shape spare the pain of going through the following
      if(vPulse->size()>1)
       {
        //1. Get the total amplitude in pe's that eventually make up the entire pulse shape 
        Float_t fPileUpAmplitude = telData->fPileUpAmplitudeForPhoton[PixelID][g];

        //2. Get the expected amplitude in mV at the FADC if everything is perfectly linear
        Float_t fLinAmplitude = fPileUpAmplitude*fLinearGainmVPerPE;

        if(bDebug)
        cout<<"time (samples): "<<time/fSamplingTime<<" signal [pe]: "<<telData->fAmplitudesInPixel[PixelID][g]<<" fPileUpAmplitude "<<fPileUpAmplitude<<" fLinAmplitude "<<fLinAmplitude<<endl;

        //3. loop over pulse shapes until we get the one that is closest to the expected amplitude
        while(fLinAmplitude>vLinearAmplitude->at(uPulseShape) && uPulseShape<vLinearAmplitude->size()-1)
           uPulseShape++;

        if(uPulseShape>0 && uPulseShape!=vPulse->size()-1 )
         {
           //get the nonlinearity factor
           //fNonLinearity = vNonLinearityFactor->at(uPulseShape-1)+
           //                      (fLinAmplitude-vLinearAmplitude->at(uPulseShape-1))/
           //                      (vLinearAmplitude->at(uPulseShape)-vLinearAmplitude->at(uPulseShape-1))
           //                      *(vNonLinearityFactor->at(uPulseShape)-vNonLinearityFactor->at(uPulseShape-1));
         
           //find out which is the better fit, this one or the previous pulse shape
           if((fLinAmplitude - vLinearAmplitude->at(uPulseShape-1))/
              (vLinearAmplitude->at(uPulseShape)-vLinearAmplitude->at(uPulseShape-1)) < 0.5)
           uPulseShape--;
           //Go back to pulse number one if you get the linear pulse shape in the above lines. This is because the linear amplitude given for the first pulse shape gives
           //the upper boundary of the linear regime. For all the other pulse shapes the linear amplitude gives the value for which the pulse shape was recorded.    
           if(uPulseShape==0)
               uPulseShape++;
           fNonLinearity = vNonLinearityFactor->at(uPulseShape);
         }
        //if(uPulseShape==vPulse->size()-1)
        // {
        //   fNonLinearity = vNonLinearityFactor->at(uPulseShape); 
        // }
       }

      if(bDebug)
      cout<<"Use pulse number: "<<uPulseShape<<" non linearity "<<fNonLinearity<<endl;       

      //figure out where we start filling the signal into the Trace
      Float_t StartTimeInTrace = time-vStartTime->at(uPulseShape);
      Int_t StartSample = Int_t(StartTimeInTrace/fSamplingTime+1);
      Int_t StopSample = Int_t((StartTimeInTrace+vStartTime->at(uPulseShape)+vStopTime->at(uPulseShape))/fSamplingTime+1);

      //The time between the start of the Trace and the first sampled value
      Float_t TimeAveragePulse= StartSample*fSamplingTime-StartTimeInTrace;

      //Catch the case when we get out of the trace window
      //Trace fully out of window
      if(StartTimeInTrace<-vStartTime->at(uPulseShape)-vStopTime->at(uPulseShape))
        continue;

      if(StartSample<0)
     {
      TimeAveragePulse = -1.0*StartTimeInTrace;
      StartTimeInTrace = 0;
      StartSample = 0;
     }
      if(StopSample>telData->iNumSamplesPerTrace)
     {
      StopSample = telData->iNumSamplesPerTrace;
     }

      // cout<<"Fill pe from: "<<StartSample<<" to "<<StopSample<<" Start time is "<<StartTime<<endl;

      //Finally fill the PE into the trace
      TimeAveragePulse = TimeAveragePulse / fSamplingTimeAveragePulse;
      Float_t step = fSamplingTime/fSamplingTimeAveragePulse;
      Float_t amplitude = telData->fAmplitudesInPixel[PixelID][g]*fNonLinearity;
      for(Int_t i=StartSample;i<StopSample;i++)
      {
       Int_t s = (Int_t)(TimeAveragePulse);
       telData->fTraceInPixel[PixelID][i]+= vPulse->at(uPulseShape)[s]*amplitude;
       TimeAveragePulse+=step;
      }

    }//go to the next photon.

  if(bDebug && telData->fTimesInPixel[PixelID].size()!=0)
   {
     cout<<"Pixel ID "<<PixelID<<endl;
     debugDisplay->Show(telData->GetTelescopeID(),PixelID);
   }

}

/////////////////////////////////////////////////////////////////
//
//  Assemble the high gain traces for all pixels
//  
void TraceGenerator::BuildAllHighGainTraces(){

 if(bDebug)
   cout<<"Generating all high gain traces"<<endl;
 //go over all pixels and generate the traces
 for(Int_t i=0;i<iNumPixels;i++)
    {
       BuildTrace(i,kFALSE);
    }
    
 //shift mean to zero, only needed if NSB is simulated
 telData->mean = 0.0;
 if(bUseNSB)
   {
     Float_t pix = 0;
     for(Int_t i=0;i<iNumPixels;i++)
       {
         //use trace if it has 1 or less Cherenkov photons in it. 
         if(telData->iPEInPixel[i]<2)
            {
               pix++;
               Double_t SumTrace = 0.0;
               for(Int_t t=0;t<telData->iNumSamplesPerTrace;t++)
	         {
	            SumTrace += telData->fTraceInPixel[i][t];
	         }
               telData->mean+= SumTrace/ (1.0*telData->iNumSamplesPerTrace);
            }
        }

    if(pix!=0)
       telData->mean=telData->mean/pix;
    else
      {
        cout<<"TraceGenerator: Warning, Could not find a pixel with less then 1 Cherenkov Photon to calculate the mean of the traces"<<endl; 
        telData->mean = 0.0;
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
   }//end shifting the mean to zero if we simulate NSB

}



/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Building the low gain trace for one pixel
//

void TraceGenerator::BuildLowGainTrace(Int_t PixelID){

  BuildTrace(PixelID,kTRUE);

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
               readConfig->GetAfterPulsingConstant(), readConfig->GetAfterPulsingSlope());


  SetSinglePESamlingWidth(readConfig->GetSampleWidthAveragePulse(iTelType));
  SetSigmaofSinglePEPulseHeightDistribution(readConfig->GetSigmaofSinglePEPulseHeightDistribution(iTelType));
 

  SetTraceSampleWidthAndLength(readConfig->GetSamplingTime(iTelType),readConfig->GetTraceLength(iTelType),readConfig->GetStartSamplingTimeOffsetFromAveragePhotonTime(iTelType));
                                            
  
  fLinearGainmVPerPE =  readConfig->GetFADCConversionFactormVperPE(iTelType);

  fPileUpWindow = readConfig->GetPileUpWindow(iTelType);  

  Float_t fFWHMofSinglePEPulse = readConfig->GetFWHMofSinglePEPulse(iTelType);

  if(fFWHMofSinglePEPulse<=0)
	{
	   SetPulseShapesFromFile(readConfig->GetNameofHighGainPulseFile(iTelType),kFALSE);
	   cout<<"Completed reading in the high gain pulse shapes from: "<<readConfig->GetNameofHighGainPulseFile(iTelType)<<endl;

	   SetPulseShapesFromFile(readConfig->GetNameofLowGainPulseFile(iTelType),kTRUE);
	   cout<<"Completed reading in the low gain pulse shapes from: "<<readConfig->GetNameofLowGainPulseFile(iTelType)<<endl;
	}
	else
	{
          SetGaussianPulse(fFWHMofSinglePEPulse);
	}
   if( fHighGainPulse.size()==0 || fLowGainPulse.size()==0 )
	  {
	     cout<<"Something went wrong, you did neither provide an external file with the pulse shape nor did you set a fwhm that can be used to generate a sample Gaus pulse"<<endl;
	     exit(1);
	  }

   //Set the Quantum efficienicy
   //loop over different types of PMTs/SiPMs and fill the array
   vector< float > f_tel;
   cout<<"Setting the QE values"<<endl;
   for(Int_t iPMTType = 0; iPMTType<readConfig->GetNumberOfPMTTypes(); iPMTType++)
    {
      qe.push_back( f_tel );
      vector<Float_t> wlOrig = readConfig->GetWavelengthsOfQEValues(iPMTType);
      vector<Float_t> qeOrig = readConfig->GetQEValues(iPMTType);
      UInt_t maxwl = (UInt_t)wlOrig[wlOrig.size()-1];
      //Set the qe values in the array
      qe[iPMTType].assign(maxwl,0);
      cout<<"for PMT/SiPM type "<<iPMTType<<endl;
      for(UInt_t w = wlOrig[0]; w<maxwl; w++)
        {
          Int_t i= 0;
          while(wlOrig[i]<w)  i++;
             i--;

          Float_t Qeff = qeOrig[i]+ (qeOrig[i+1]-qeOrig[i])/(wlOrig[i+1]-wlOrig[i])*(w-wlOrig[i]);

          qe[iPMTType][w] = Qeff;

          //cout<<w<<"  "<<qe[iPMTType][w]<<endl;
        }
    }//end looping over all different PMT types and loading the QEs 

  //NSB
  vNSBRatePerPixel = readConfig->GetNSBRate(iTelType);      //the NSB rate in kHz per camera pixel;
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

  iTubeType = readConfig->GetTubeType(iTelType);

}

//-------------------------------------------------------------------------------------------
//
// Get simulated average single pe pulse. Returns a root histogram
TH1F* TraceGenerator::GetHighGainPulse(UInt_t n)
{

  if(fHighGainPulse.empty())
    {
      cout<<"High gain pulses have not been set. Set a pulse shape first"<<endl;
      return NULL;
    }

 if(fHighGainPulse.size()<=n)
    {
      cout<<"You are asking for a pulse shape that does not exist"<<endl;
      return NULL;
    }


  TString title;
  title.Form("High gain pulse %i, start: %2.1f, stop: %2.1f, AreaToPeak: %1.2f, amplitude: %2.2f, non-linearity: %2.2f"
     ,n,fHighGainStartTime[n],fHighGainStopTime[n],fHighGainAreaToPeakConversion[n],fHighGainLinearAmplitude[n],fHighGainNonLinearityFactor[n]);

  TString name;
  name.Form("h%i",n);

  TH1F *h = new TH1F(name.Data(),title.Data(),fHighGainPulse[n].size(),-1*fHighGainStartTime[n]-fSamplingTimeAveragePulse*0.5,fHighGainStopTime[n]+ fSamplingTimeAveragePulse*0.5);
  h->GetXaxis()->SetTitle("Time [ns]");
  h->GetYaxis()->SetTitle("Amplitude normalized to peak");
  for(UInt_t i = 0; i<fHighGainPulse[n].size(); i++)
    {
      //cout<<i<<"  "<<fHighGainPulse[n][i]<<endl;
      h->SetBinContent(i+1,fHighGainPulse[n][i]);
    }


  return h;
}

//-------------------------------------------------------------------------------------------
//
// Get low gain pulse. Returns a root histogram
TH1F* TraceGenerator::GetLowGainPulse(UInt_t n)
{

 if(fLowGainPulse.empty())
    {
      cout<<"Low gain pulses have not been set. Set a pulse shape first"<<endl;
      return NULL;
    }

 if(fLowGainPulse.size()<=n)
    {
      cout<<"You are asking for a pulse shape that does not exist"<<endl;
      return NULL;
    }

  TString title;
  title.Form("Low gain pulse %i, start: %2.1f, stop: %2.1f, AreaToPeak: %1.2f, amplitude: %2.2f, non-linearity: %2.2f"
     ,n,fLowGainStartTime[n],fLowGainStopTime[n],fLowGainAreaToPeakConversion[n],fLowGainLinearAmplitude[n],fLowGainNonLinearityFactor[n]);

  TString name;
  name.Form("h%i",n);

  TH1F *h = new TH1F(name.Data(),title.Data(),fLowGainPulse[n].size(),-1*fLowGainStartTime[n]-fSamplingTimeAveragePulse*0.5,fLowGainStopTime[n]+ fSamplingTimeAveragePulse*0.5);
  h->GetXaxis()->SetTitle("Time [ns]");
  h->GetYaxis()->SetTitle("Amplitude normalized to peak");
  for(UInt_t i = 0; i<fLowGainPulse[n].size(); i++)
    {
      //cout<<i<<"  "<<fLowGainPulse[n][i]<<endl;
      h->SetBinContent(i+1,fLowGainPulse[n][i]);
    }


  return h;
}



//-------------------------------------------------------------------------------------------
//
// shows the average single pe pulse shape
void  TraceGenerator::ShowPulseShape(UInt_t n)
{



  TCanvas *cSinglePEShape = new TCanvas("cSinglePEShape","The pulse shape used in the simulation",1000,700);
  
  cSinglePEShape->Divide(2,1);
  cSinglePEShape->Draw();
  cSinglePEShape->cd(1);
  gPad->SetGrid();

  TH1F *h = NULL;
  if(fHighGainPulse.size()>n)
  {	   
      h = GetHighGainPulse(n);
      h->GetXaxis()->SetTitle("Time [ns]");
      h->GetYaxis()->SetTitle("Amplitude [photoelectrons]");
      h->Draw();  
  }

  cSinglePEShape->cd(2);
  gPad->SetGrid();

  TH1F *hL = NULL;
  if(fLowGainPulse.size()>n)
  {
     hL = GetLowGainPulse(n);
     hL->GetXaxis()->SetTitle("Time [ns]");
     hL->GetYaxis()->SetTitle("Amplitude [photoelectrons]");
     hL->Draw();  
  }
	      
  TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
	      
  //
  // While reading the input process gui events asynchronously
  //
  timer.TurnOn();
  TString input = Getline("Type <return> to go on: ");
  timer.TurnOff(); 
	      
}


