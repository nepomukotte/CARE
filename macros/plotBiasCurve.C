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

void plotBiasCurve(string DataDir,const Int_t iMaxTelescopeMultiplicity = 2, Float_t fCoincidence = 50, Int_t iNumTelNextNeighborGroup =8)
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
  TVectorD BiasCurveError(NumPoints);
  BiasCurveError.Zero();
 
  TVectorD GroupRateVsThreshold(NumPoints);
  GroupRateVsThreshold.Zero();
  TVectorD GroupRateVsThresholdError(NumPoints);
  GroupRateVsThresholdError.Zero();


  TVectorD  FakeXErrors(NumPoints);
  FakeXErrors.Zero();
  
  //Loop over all files
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
     

      if(BiasCurveValuesThis!=0 && BiasCurveErrorThis!=0 && BiasCurveScanPointsThis!=0 
         && GroupRateVsThresholdThis!=0 && GroupRateVsThresholdErrorThis!=0 )
	{
	  for(Int_t p = 0; p<NumPoints; p++)
	    {
	      BiasCurveValues[p]+= (*BiasCurveValuesThis)[p];
	      BiasCurveError[p]+= (*BiasCurveErrorThis)[p] * (*BiasCurveErrorThis)[p];

	      GroupRateVsThreshold[p]+= (*GroupRateVsThresholdThis)[p];
	      GroupRateVsThresholdError[p]+= (*GroupRateVsThresholdErrorThis)[p] * (*GroupRateVsThresholdErrorThis)[p];

	    } 
	  NumFiles++;
	}

      f->Close();

    }

  for(Int_t p = 0; p<NumPoints; p++)
    {
      BiasCurveValues[p] = BiasCurveValues[p]/NumFiles;
      if( BiasCurveError[p]>0 )
	BiasCurveError[p] =  sqrt(BiasCurveError[p])/NumFiles;

      GroupRateVsThreshold[p] = GroupRateVsThreshold[p]/NumFiles;
      if( GroupRateVsThresholdError[p]>0 )
	GroupRateVsThresholdError[p] =  sqrt(GroupRateVsThresholdError[p])/NumFiles;

       BiasCurveScanPoints[p]=-1*BiasCurveScanPoints[p];

      cout<<p<<" rate:  "<<BiasCurveValues[p]<<"+-"<<BiasCurveError[p]<<"  "<<BiasCurveScanPoints[p]<<endl;
      cout<<"Group rate "<<GroupRateVsThreshold[p]<<"+-"<<GroupRateVsThresholdError[p]<<"  "<<endl;
  } 

  Int_t lowestBin=0;
  while(BiasCurveError[lowestBin]>0 && lowestBin<NumPoints-1)
    lowestBin++;
  
  lowestBin--;

  TCanvas *c = new TCanvas("c","",1000,800);
  c->SetLogy();
  c->SetGrid();
  gPad->SetFrameBorderMode(0);
  c->Draw();
  

  TH1D *h = new TH1D("h","Bias Curve",1000,BiasCurveScanPoints.Min()*0.95,BiasCurveScanPoints[lowestBin]*1.1);
  h->SetMaximum(BiasCurveValues.Max()*5);
  h->SetMinimum(1);
  h->GetXaxis()->SetTitle("Discriminator Threshold [mV]");
  h->GetYaxis()->SetTitle("Trigger Rate [Hz]");
  h->SetStats(0);
  h->Draw();
  TLine *lMaxAccidentals = new TLine( BiasCurveScanPoints.Min()*0.95 , 10 , BiasCurveScanPoints[lowestBin]*1.1, 10);
  lMaxAccidentals->SetLineStyle(7);
  lMaxAccidentals->SetLineWidth(2);
  lMaxAccidentals->Draw();

  TGraphErrors *gr = new TGraphErrors(BiasCurveScanPoints,BiasCurveValues,FakeXErrors,BiasCurveError);
  gr->SetMarkerStyle(29);
  gr->SetMarkerSize(1.6);
  gr->SetLineWidth(2);

  TLegend *legend = new TLegend(0.55,0.55,0.88,0.88);
  legend->AddEntry(gr,"Single Telescope Rate 3 fold multiplicity","pl");

  TGraphErrors *grNoNextNeighbor[iMaxTelescopeMultiplicity-1];

  TGraphErrors *grCompactNextNeighbor[iMaxTelescopeMultiplicity-1];

  TGraphErrors *grLooseNextNeighbor[iMaxTelescopeMultiplicity-1];
  Double_t dAverageTelescopes[] = { 8, 11, 13.3 };

  for(Int_t m =2;m<=iMaxTelescopeMultiplicity;m++)
    {

      //Array Trigger Bias Curves
  
      //no next neighbor
      grNoNextNeighbor[m-2] = new TGraphErrors();
      grNoNextNeighbor[m-2]->SetMarkerStyle(20+(m-2));
      grNoNextNeighbor[m-2]->SetMarkerColor(2+(m-2));
      grNoNextNeighbor[m-2]->SetMarkerSize(1.4);
      grNoNextNeighbor[m-2]->SetLineWidth(2);
      grNoNextNeighbor[m-2]->SetLineColor(2+(m-2));

      grCompactNextNeighbor[m-2] = (TGraphErrors*)grNoNextNeighbor[m-2]->Clone();
      grCompactNextNeighbor[m-2]->SetLineStyle(9);
      grCompactNextNeighbor[m-2]->SetMarkerStyle(24+(m-2));
      grCompactNextNeighbor[m-2]->SetMarkerSize(2.3);

      grLooseNextNeighbor[m-2] = (TGraphErrors*)grNoNextNeighbor[m-2]->Clone();
      grLooseNextNeighbor[m-2]->SetLineStyle(2);
      grLooseNextNeighbor[m-2]->SetMarkerStyle(2+(m-2));

      TString label;
      label.Form("Any telescope Combination; Telescope multiplicity %i",m);
      //  legend->AddEntry(grNoNextNeighbor[m-2],label,"pl");
      label.Form("Loose group;       multiplicity %i",m);
      //  legend->AddEntry(grLooseNextNeighbor[m-2],label,"pl");
      label.Form("Compact group;  multiplicity %i",m);
      legend->AddEntry(grCompactNextNeighbor[m-2],label,"pl");
     

      //cout<<NumTelescopes<<"  "<<iNumTelNextNeighborGroup<<endl;

      Double_t lowestrate = 0;
      Double_t highestrate = 0;
      for(Int_t p = 0; p<NumPoints;p++)
	{
	  Double_t rate = BiasCurveValues[p];
	  Double_t rateCompactNextNeighbor = BiasCurveValues[p];
	  Double_t rateLooseNextNeighbor = BiasCurveValues[p];
	  for(Int_t t=1;t<m;t++)
	    {
	      cout<<m<<endl;
	      Double_t prob = 1.0-exp(-BiasCurveValues[p]*fCoincidence*1e-9*(NumTelescopes-t));
	      rate*=prob;

	      //compact next neighbor
	      prob = 1.0-exp(-BiasCurveValues[p]*fCoincidence*1e-9*(iNumTelNextNeighborGroup-t));
	      rateCompactNextNeighbor*=prob;

	      //loose next neighbor
	      prob = 1.0-exp(-BiasCurveValues[p]*fCoincidence*1e-9*dAverageTelescopes[t-1]);
	      rateLooseNextNeighbor*=prob;

	    }
      
	  if(rate!=0)
	    lowestrate=rate;
      
	  //cout<<rate<<"  "<<rateCompactNextNeighbor<<endl;
	  grNoNextNeighbor[m-2]->SetPoint(p,BiasCurveScanPoints[p],rate);
	  if(BiasCurveValues[p]>0)
	    grNoNextNeighbor[m-2]->SetPointError(p,0,rate*BiasCurveError[p]/BiasCurveValues[p]);

	  //compact next neighbor
	  grCompactNextNeighbor[m-2]->SetPoint(p,BiasCurveScanPoints[p],rateCompactNextNeighbor);
	  if(BiasCurveValues[p]>0)
	    grCompactNextNeighbor[m-2]->SetPointError(p,0,rateCompactNextNeighbor*BiasCurveError[p]/BiasCurveValues[p]);

	  //loose next neighbor
	  grLooseNextNeighbor[m-2]->SetPoint(p,BiasCurveScanPoints[p],rateLooseNextNeighbor);
	  if(BiasCurveValues[p]>0)
	    grLooseNextNeighbor[m-2]->SetPointError(p,0,rateLooseNextNeighbor*BiasCurveError[p]/BiasCurveValues[p]);

      
	}
       //grNoNextNeighbor[m-2]->Draw("PL");
       grCompactNextNeighbor[m-2]->Draw("PL");
      // grLooseNextNeighbor[m-2]->Draw("PL");
    }


  


  gr->Draw("PL");
 
  legend->Draw();
  
  //Single Group Rate s. Threshold

  lowestBin=0;
  while(GroupRateVsThresholdError[lowestBin]>0 && lowestBin<NumPoints-1)
    lowestBin++;
  
  lowestBin--;

  TCanvas *cg = new TCanvas("cg","",1000,800);
  cg->SetLogy();
  cg->SetGrid();
  gPad->SetFrameBorderMode(0);
  cg->Draw();
  

  TH1D *hg = new TH1D("hg","Bias Curve of a single group",1000,BiasCurveScanPoints.Min()*0.95,BiasCurveScanPoints[lowestBin]*1.1);
  hg->SetMaximum(GroupRateVsThreshold.Max()*5);
  hg->SetMinimum(1);
  hg->GetXaxis()->SetTitle("Discriminator Threshold [mV]");
  hg->GetYaxis()->SetTitle("Trigger Rate [Hz]");
  hg->SetStats(0);
  hg->Draw();


  TGraphErrors *grg = new TGraphErrors(BiasCurveScanPoints,GroupRateVsThreshold,FakeXErrors,GroupRateVsThresholdError);
  grg->SetMarkerStyle(29);
  grg->SetMarkerSize(1.6);
  grg->SetLineWidth(2);
  grg->Draw("PL");

  TLegend *legendg = new TLegend(0.55,0.55,0.88,0.88);
  legendg->AddEntry(grg,"Single Group Rate","pl");
  legendg->Draw();

      //compact next neighbor
      
      //next neighbor

}
