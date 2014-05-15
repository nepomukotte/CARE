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

  cTrace = NULL;
  cFADCTrace = NULL;
  hFADCTrace = NULL; 
  hTrace = NULL; 


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

void Display::ResetTriggerTraces()
{
 vDiscThreshold.clear();
 vDiscTelID.clear();
 vDiscPixel.clear();
 vHCFDTrace.clear();
 vHThreshTrace.clear();
}


//Add one discriminator pixel that is supposed to be displayed
void Display::AddDiscriminatorTraces(int TelID, int triggerpixel,float threshold, TH1F hThresholdTrace,TH1F hCFDTrace)
{
   //Add the two histograms to a vector of histograms
   vDiscThreshold.push_back(threshold);
   vDiscTelID.push_back(TelID);
   vDiscPixel.push_back(triggerpixel);
   vHCFDTrace.push_back(hCFDTrace);
   vHThreshTrace.push_back(hThresholdTrace);
}

void Display::ShowSelectedDiscriminatorPixels()
{
   vector<TCanvas*> vCanvases;
   vector<TLine> vThreshLvls;
   vThreshLvls.resize(vDiscPixel.size());
   for(int t=0;t<configuration->GetNumberOfTelescopes();t++)
    {
     TString title;
     title.Form("cDisc%i",t);
     vCanvases.push_back(new TCanvas(title.Data(),"Traces selected from Discriminator",1000,500));
     title.Form("Traces selected from Discriminator for Telescope %i",t);
     vCanvases[t]->SetTitle(title.Data());
     cout<<"Pixel that triggered Telescope "<<t<<endl;
     vector<int> trigclust = allTelData[t]->GetTriggerCluster();
     vCanvases[t]->Divide(5,(trigclust.size()-1)/5+1);
     int iPix=0;
     for(unsigned i=0;i<trigclust.size();i++)
      {  
        cout<<trigclust[i]<<", ";

        for(unsigned p=0;p<vDiscPixel.size();p++)
          {
            if(vDiscPixel[p]==trigclust[i] && vDiscTelID[p]==t)
              {
                cout<<iPix<<" plotting "<<vDiscPixel[p]<<" telescope "<<vDiscTelID[p]<<endl;
                vCanvases[t]->cd(iPix+1);
                gPad->SetGrid();
                vHThreshTrace[p].SetMaximum(vHCFDTrace[p].GetMaximum()>vHThreshTrace[p].GetMaximum() ? 
                                  vHCFDTrace[p].GetMaximum()+10:vHThreshTrace[p].GetMaximum()+10 );
                vHThreshTrace[p].SetMinimum(vHCFDTrace[p].GetMinimum()<vHThreshTrace[p].GetMinimum() ? 
                                  vHCFDTrace[p].GetMinimum()-10:vHThreshTrace[p].GetMinimum()-10 );

                vHThreshTrace[p].Draw();
                vHCFDTrace[p].Draw("same");
                //cout<<vHThreshTrace[i].GetBinLowEdge(1)<<"  "<<vHThreshTrace[i].GetBinLowEdge(vHThreshTrace[i].GetNbinsX())<<" "<<vDiscThreshold[i]<<endl;
                vThreshLvls[p]=TLine(vHThreshTrace[i].GetBinLowEdge(1),vDiscThreshold[i],
	        vHThreshTrace[p].GetBinLowEdge(vHThreshTrace[i].GetNbinsX()), vDiscThreshold[i]);
                vThreshLvls[p].Draw();
                iPix++;
               }
    
      
    }
      }
     cout<<endl<<endl;
  }
  //hold the code;
  TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
  timer.TurnOn();
  TString input = Getline("Type <return> to go on: ");
  timer.TurnOff();
  vCanvases.clear();
}

void Display::ShowFADC(int TelID = 0,int PixID = 0 )
{

    //if we have not yet created the plots showing the camera parameters
    if(vTelParCanvases.size()==0)
        MakeCameraParameterPlots();
	
      int iNumSamples = allTelData[TelID]->iNumFADCSamples;
	{
	   if(!cFADCTrace)
		 cFADCTrace = new TCanvas("cFADCTrace","Trace",500,500);
	   
	   cFADCTrace->Draw();
	   cFADCTrace->cd();
	   if(hFADCTrace)
		  hFADCTrace->Delete();
           hFADCTrace = new TH1F("hFADCTrace","FADC Trace",iNumSamples,0,iNumSamples);
	   vector<Int_t> iFADCTrace = allTelData[TelID]->iFADCTraceInPixel[PixID];
	   for(int i = 0; i<iNumSamples;i++)
		  hFADCTrace->SetBinContent(i+1,iFADCTrace[i]);
	   hFADCTrace->Draw();
	   cFADCTrace->Modified();
	   cFADCTrace->Update();

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

void Display::Show(int TelID = 0,int PixID = 0 )
{
    cout<<"Showing TelID "<<TelID<<"  Pixel: "<<PixID<<endl;


    //if we have not yet created the plots showing the camera parameters
    if(vTelParCanvases.size()==0)
        MakeCameraParameterPlots();

	
	int iNumAnalogSamples = allTelData[TelID]->iNumSamplesPerTrace;
	{
	   if(!cTrace)
		 cTrace = new TCanvas("cTrace","Trace",500,500);
	   
	   cTrace->Draw();
	   cTrace->cd();
	   if(hTrace)
		  hTrace->Delete();
           hTrace = new TH1F("hTrace","Trace",iNumAnalogSamples,0,iNumAnalogSamples);
	   vector<Float_t> fTraceInPixel = allTelData[TelID]->fTraceInPixel[PixID];
	   for(int i = 0; i<iNumAnalogSamples;i++)
		  hTrace->SetBinContent(i+1,fTraceInPixel[i]);
	   hTrace->Draw();
	   cTrace->Modified();
	   cTrace->Update();

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

