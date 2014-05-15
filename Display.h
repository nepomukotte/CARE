#ifndef Display_Class
#define Display_Class


#include <TH1F.h>		
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TGraph.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TTimer.h"
#include "TMath.h"

#include "TelescopeData.h"

using namespace std;

class Display {

 public:
 
  Display( int argc, char **argv,ReadConfig *readConfig);
  void ShowSelectedDiscriminatorPixels();
  void Show(int TelID, int pixel  );
  void ShowFADC(int TelID, int pixel  );
  void SetTelescopeData(TelescopeData **AllTelescopes){ allTelData = AllTelescopes;};
  void AddCanvas(TCanvas *c);
  void AddDiscriminatorTraces(int TelID, int triggerpixel,float threshold,TH1F hThresholdTrace,TH1F hCFDTrace);
  void ResetTriggerTraces();
 protected:

 void MakeCameraParameterPlots();

 TApplication *theApp;
 
 ReadConfig *configuration;

 Bool_t HandleInput();
 TCanvas *cTrace;
 TCanvas *cFADCTrace;
 TH1F *hTrace;   //Histogram to show the trace of a pixel   
 TH1F *hFADCTrace;   //Histogram to show the FADC trace of a pixel   

 TelescopeData **allTelData;
 vector<TCanvas*> vCanvasesToShow;
 vector<TCanvas*> vTelParCanvases;

 vector<float>  vDiscThreshold;
 vector<int>  vDiscTelID;
 vector<int>  vDiscPixel;
 vector<TH1F> vHCFDTrace;
 vector<TH1F> vHThreshTrace;

};

#endif

