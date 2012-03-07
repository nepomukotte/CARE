/* \file Display.cpp
   Implementation of a class to view events, mainly for debugging
*/
#include <iostream>
#include <string>
#include <vector>
#include <Getline.h>

#include "Display.h"
#include <TCanvas.h>
#include <TH1F.h>
#include <TLine.h>
#include <TROOT.h>
#include <TApplication.h>

using namespace std;

//---------------------------------------------------------------------------------------
//Constructor
Display::Display( int argc, char **argv,ReadConfig *readConfig )
{
  cout<<"Initializing the Display"<<endl;


   //needed to display stuff while CARE executes
   TROOT root("DisplayEvts","Display Results");
   theApp = new TApplication("App",&argc,argv);

  configuration = readConfig;  

  cTraceNSB = NULL;
  hTraceNSB = NULL; 


}

Bool_t Display::HandleInput()
{
  TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
  while (1)
    {
      //
      // While reading the input process gui events asynchronously
      //
      timer.TurnOn();
      TString input = Getline("Type 'q' to exit, <return> to go on: ");
      timer.TurnOff();

      if (input=="q\n")
        return kFALSE;

      if (input=="\n")
        return kTRUE;
    };

  return kFALSE;
}

void Display::Show(int TelID = 0,int PixID = 0 )
{
    cout<<"Showing TelID "<<TelID<<"  Pixel: "<<PixID<<endl;


    //if we have not yet created the plots showing the camera parameters
    if(vTelParCanvases.size()==0)
        MakeCameraParameterPlots();

	
	int iNumAnalogSamples = allTelData[TelID]->iNumSamplesPerTrace;
	{
       //NSB Trace
	   if(!cTraceNSB)
		 cTraceNSB = new TCanvas("cTraceNSB","NSB Trace",500,500);
	   
	   cTraceNSB->Draw();
	   cTraceNSB->cd();
	   if(hTraceNSB)
		  hTraceNSB->Delete();
       hTraceNSB = new TH1F("hTraceNSB","NSB Trace",iNumAnalogSamples,0,iNumAnalogSamples);
	   vector<Float_t> fTraceInPixel = allTelData[TelID]->fTraceInPixel[PixID];
	   for(int i = 0; i<iNumAnalogSamples;i++)
		  hTraceNSB->SetBinContent(i,fTraceInPixel[i]);
	   hTraceNSB->Draw();
	   cTraceNSB->Modified();
	   cTraceNSB->Update();

       for(UInt_t i=0;i<vCanvasesToShow.size();i++)
         {
           vCanvasesToShow[i]->Draw();
         }

       for(UInt_t i=0;i<vTelParCanvases.size();i++)
         {
           vTelParCanvases[i]->Draw();
         }

       //hold the code;
	   TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
	   timer.TurnOn();
	   TString input = Getline("Type <return> to go on: ");
	   timer.TurnOff();

    }


}

void Display::AddCanvas(TCanvas *c)
{
  vCanvasesToShow.push_back(c);
}

void Display::MakeCameraParameterPlots()
{

  cout<<"Plotting the camera distribution plots"<<endl;
  
  for(int t=0;t<configuration->GetNumberOfTelescopes();t++)
  {
     cout<<endl<<endl;
     TString title;
     title.Form("cTelDistr%i",t);
     TString name;
     name.Form("Distributions of Telescope %i",t);
     TCanvas *cTelDistr = new TCanvas(title.Data(),name.Data(),1000,1000);
     vTelParCanvases.push_back(cTelDistr);

     cTelDistr->Divide(2,2);

     //The relative QE distribution
     cTelDistr->cd(1);
     title.Form("hRelQE%i",t);
     TH1F *hRelQE = new TH1F(title.Data(),"relative QE distribution",100,-0.2,3);
     for(int p=0;p<allTelData[t]->iNumPixels;p++)
       hRelQE->Fill(allTelData[t]->fRelQE[p]);

     cout<<"Fitting rel. QE distribution of telescope 0"<<endl;
     hRelQE->Fit("gaus");
     hRelQE->Draw();

     //The relative gain distribution
     cTelDistr->cd(2);
     title.Form("hRelGain%i",t);
     TH1F *hRelGain = new TH1F(title.Data(),"relative absolute Gain distribution",100,-0.2,3);
     for(int p=0;p<allTelData[t]->iNumPixels;p++)
       hRelGain->Fill(allTelData[t]->fRelGain[p]);

     cout<<"Fitting rel. Gain distribution of telescope 0"<<endl;
     hRelGain->Fit("gaus");
     hRelGain->Draw();

     //The relative gain times QE distribution
     cTelDistr->cd(3);
     title.Form("hRelQEtimesGain%i",t);
     TH1F *hRelQEtimesGain = new TH1F(title.Data(),"relative absolute Gain times rel. QE distribution",100,0,2);
     for(int p=0;p<allTelData[t]->iNumPixels;p++)
        hRelQEtimesGain->Fill(allTelData[t]->fRelGain[p]*allTelData[t]->fRelQE[p]);

     cout<<"Fitting rel. Gain times rel QE distribution of telescope 0"<<endl;
     hRelQEtimesGain->Fit("gaus");
     hRelQEtimesGain->Draw();

     //The relative QE distribution  times the WC efficiency
     cTelDistr->cd(4);
     title.Form("hRelQEwWC%i",t);
     TH1F *hRelQEwWC = new TH1F(title.Data(),"relative QE times the Winston cone efficiency",100,-0.2,3);
     for(int p=0;p<allTelData[t]->iNumPixels;p++)
        hRelQEwWC->Fill(allTelData[t]->fRelQEwWC[p]);

     cout<<"Fitting rel. QE distribution  times the winston cone efficiency of telescope 0"<<endl;
     hRelQEwWC->Fit("gaus");
     hRelQEwWC->Draw();

     cout<<endl<<endl;
  }
 

}

