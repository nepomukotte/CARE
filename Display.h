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
  void Show(int TelID, int pixel  );
  void SetTelescopeData(TelescopeData **AllTelescopes){ allTelData = AllTelescopes;};
  void AddCanvas(TCanvas *c);

 protected:

 void MakeCameraParameterPlots();

 TApplication *theApp;
 
 ReadConfig *configuration;

 Bool_t HandleInput();
 TCanvas *cTraceNSB;
 TH1F *hTraceNSB;   //Histogram to show the NSB trace of a pixel   

 TelescopeData **allTelData;
 vector<TCanvas*> vCanvasesToShow;
 vector<TCanvas*> vTelParCanvases;

};

#endif

