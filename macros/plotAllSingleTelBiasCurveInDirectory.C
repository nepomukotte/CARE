#include <TVectorD>
//////////////////////////////////////


void readDir(string dirname,vector<string> *names) {

    string file;
  void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(dirname.c_str()));
    const char* strptr;
    if (dir) {
        while ((strptr = gSystem->GetDirEntry(dir))) {
          file = strptr;
          //if (file=="." || file==".." || file=="log") continue;
          //names->push_back(file);
          size_t found;
          found=file.find(".root");
          if (found==string::npos) continue;
          names->push_back(file);
        }
        gSystem->FreeDirectory(dir);
    }
     // for(size_t i(0); i<names->size(); i++)
     //   Printf("%s",(*names)[i].c_str());
}

/////////////////////////////////////////

void plotAllSingleTelBiasCurveInDirectory(string DataDir,const Int_t iMaxTelescopeMultiplicity = 4, Float_t fCoincidence = 50, Int_t iNumTelNextNeighborGroup =4)
{

  //const Int_t iMaxTelescopeMultiplicity = 4;
  //Float_t fCoincidence = 50; //coincidence window in ns
  

  //Get all the filenames
  vector<string> * filenames = new vector<string>; 
  readDir(DataDir,filenames);

 
  stringstream fnamestream;
  fnamestream<<DataDir<<"/"<<(*filenames)[0].c_str();

  string fstring(fnamestream.str());
  cout<<"opening file: "<<fstring<<endl;

  gStyle->SetCanvasColor(0);


  TVectorD BiasCurveScanPoints;


  TFile *f = new TFile(fstring.c_str(),"READ");
  
  //  f->cd("BiasCurve");
 
  Int_t NumPoints = ((TVectorD*)f->Get("BiasCurve/TVecBiasCurve"))->GetNoElements();
  cout<<"Number of points per curve "<<NumPoints<<endl;

   BiasCurveScanPoints.ResizeTo(NumPoints);
   BiasCurveScanPoints.Zero();

  TVectorF *TVecNumTelescopes = (TVectorF*)f->Get("BiasCurve/TVecNumTelescopes");
  Int_t NumTelescopes = (Int_t)((*TVecNumTelescopes)[0]+0.5);


  TVectorD *BiasCurveScanPointsThis = (TVectorD*)f->Get("BiasCurve/TVecBiasCurveScanPoints");
  BiasCurveScanPoints.SetElements(BiasCurveScanPointsThis->GetMatrixArray());  

  /*
      TVectorF TVecXpos = ((TVectorF*)f->Get("BiasCurve/TVecXpos"))->Clone();
      TVectorF TVecYpos = ((TVectorF*)f->Get("BiasCurve/TVecYpos"))->Clone();
  */

  f->Close();

  TVectorD BiasCurveValues(NumPoints);
  BiasCurveValues.Zero();
 
 
  TVectorD GroupRateVsThreshold(NumPoints);
  GroupRateVsThreshold.Zero();
  TVectorD GroupRateVsThresholdError(NumPoints);
  GroupRateVsThresholdError.Zero();


  TVectorD  FakeXErrors(NumPoints);
  FakeXErrors.Zero();
  
  //Loop over all files

  const Int_t iNumFiles = 50;

  //filenames->size();

  TGraphErrors *grSingleTelBiasCurve[iNumFiles];

  Int_t NumFiles = 0;
  for(size_t i(0); i<filenames->size(); i++)
    {

       stringstream fnamestream;
      fnamestream<<DataDir<<"/"<<(*filenames)[i].c_str();
      string fstring(fnamestream.str());
      cout<<"opening file: "<<fstring<<endl;
      f = new TFile(fstring.c_str(),"READ");
      if( f->IsZombie() )
         {
           cout << "error opening root input file: " <<fnamestream.str().c_str()  << endl;
           cout << "...skipping" << endl;
           continue; //exit( -1 );
         }

      TVectorD *BiasCurveValuesThis = (TVectorD*)f->Get("BiasCurve/TVecBiasCurve");
      TVectorD *BiasCurveErrorThis  = (TVectorD*)f->Get("BiasCurve/TVecBiasCurveError");
   
      TVectorD *GroupRateVsThresholdThis = (TVectorD*)f->Get("BiasCurve/TVecGroupRateVsThreshold");
      TVectorD *GroupRateVsThresholdErrorThis  = (TVectorD*)f->Get("BiasCurve/TVecGroupRateVsThresholdError");
     
      grSingleTelBiasCurve[i] = new TGraphErrors(BiasCurveScanPoints,*BiasCurveValuesThis,FakeXErrors,*BiasCurveErrorThis);
      grSingleTelBiasCurve[i]->SetMarkerStyle(29);
      grSingleTelBiasCurve[i]->SetMarkerSize(1.6);
      grSingleTelBiasCurve[i]->SetLineWidth(2);    

    }

  Double_t *BiasCurveError = grSingleTelBiasCurve[0]->GetEY();

  Int_t lowestBin=0;
  while(BiasCurveError[lowestBin]>0 && lowestBin<NumPoints-1)
    lowestBin++;
  
  lowestBin--;

 

  TCanvas *c = new TCanvas("c","",1000,800);
  c->SetLogy();
  c->SetGrid();
  gPad->SetFrameBorderMode(0);
  c->Draw();
  

  TH1D *h = new TH1D("h","Single Tel Bias Curves 120 MHZ NSB, 9ns Coinc, 6ns FWHM signals, XP2970 with AP",1000,1*0.95,6*1.1);
  h->SetMaximum(1e7*5);
  h->SetMinimum(1);
  h->GetXaxis()->SetTitle("CFD Threshold [pe]");
  h->GetYaxis()->SetTitle("Trigger Rate [Hz]");
  h->SetStats(0);
  h->Draw();
  TLine *lMaxAccidentals = new TLine( 1*0.95 , 10 , 6*1.1, 10);
  lMaxAccidentals->SetLineStyle(7);
  lMaxAccidentals->SetLineWidth(2);
  lMaxAccidentals->Draw();

   TLegend *legend = new TLegend(0.55,0.55,0.88,0.88);
 

  
 for(size_t i(0); i<filenames->size(); i++)
    {
      
      cout<<i<<endl;
      grSingleTelBiasCurve[i]->SetMarkerStyle(20+i);
      grSingleTelBiasCurve[i]->SetMarkerColor(2+i);
      grSingleTelBiasCurve[i]->SetMarkerSize(1.4);
      grSingleTelBiasCurve[i]->SetLineWidth(2);
      grSingleTelBiasCurve[i]->SetLineColor(2+i);
      grSingleTelBiasCurve[i]->Draw("PL");

      
      TString label;
      label.Form("RFB ");
      label+=(*filenames)[i].c_str());
        legend->AddEntry(grSingleTelBiasCurve[i],label,"pl");
    
      
    }

  
 
TGraphErrors *grMeasured =new TGraphErrors();
grMeasured->SetPoint(0,2.5,8e5);
grMeasured->SetPointError(0,2.5*0.1,8e5*0.6);
grMeasured->SetPoint(1,3.3,2.5e4);
grMeasured->SetPointError(1,3.3*0.1,2.5e4*0.6);
grMeasured->SetPoint(2,4.2,1.5e3);
grMeasured->SetPointError(2,4.2*0.1,1.5e3*0.6);

grMeasured->SetMarkerStyle(29);
//grMeasured->SetMarkerColor(kRed-7);
grMeasured->SetMarkerSize(2.9);
grMeasured->SetLineWidth(2);
//grMeasured->SetLineColor(kRed-7);
grMeasured->Draw("P");

 legend->AddEntry(grMeasured,"Measured Bias Curve Run 50639","p");
legend->Draw();
}
