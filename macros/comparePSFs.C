void comparePSFs()
{

//string sFileData = "T2.root";
string sFileData = "/nv/hp11/aotte6/data/analysis_projects/comissioning/opticalPSF/results/T1.root";
//string sFileSims ="PSFgraph.root";
string sFileSims ="/nv/hp11/aotte6/code/GrOptics/PSFgraph.root";

TFile *fData = new TFile(sFileData.c_str());
fData->ls();

//TCanvas *cData = (TCanvas*)fData->Get("T4FB_120913_051907_fits");
//TCanvas *cData = (TCanvas*)fData->Get("T3FB_120913_052548_fits");
//TCanvas *cData = (TCanvas*)fData->Get("T2FB_120913_052453_fits");
TCanvas *cData = (TCanvas*)fData->Get("T1FB_120913_071517_fits");
cData->Draw();

cData->ls();
TH1F *hr1data = cData->FindObject("hr1");
TH1F *hr1sims = (TH1F*)hr1data->Clone("hr1sims");

TH1F *hr2data = cData->FindObject("hr2");
TH1F *hr2sims = (TH1F*)hr1data->Clone("hr2sims");

TFile *fSims = new TFile(sFileSims.c_str());
fSims->ls();

TH2D *hSims = (TH2D*)fSims->Get("hist0");

TCanvas *cSims = new TCanvas("cSims","Simulated PSF",1000,1000);
cSims->Divide(2,2);
cSims->Draw();
cSims->cd(1);
gStyle->SetPalette(1);
hSims->SetStats(0);
hSims->Draw("colz");

//Get the radial distribution of the PSF
hr1sims->Reset();
for(unsigned int x = 0; x< hSims->GetNbinsX(); x++)
{
  for(unsigned int y = 0; y< hSims->GetNbinsY(); y++)
    {
      Double_t BinContent = hSims->GetBinContent(x,y);
      Double_t xCoord = hSims->GetXaxis()->GetBinCenter(x);
      xCoord=xCoord/60.; //convert from arcmin to degreees
      Double_t yCoord = hSims->GetYaxis()->GetBinCenter(y);
      yCoord=yCoord/60.; //convert from arcmin to degreees
      Double_t dDist=sqrt(xCoord*xCoord+yCoord*yCoord);
      Int_t nBin=hr1sims->FindBin(dDist);
     // cout<<nBin<<"  "<<dDist<<"  "<<BinContent<<endl;
      //Fill radial intensity histogram
      hr1sims->SetBinContent(nBin,hr1sims->GetBinContent(nBin)+BinContent);
    }
}

//normalize
Float_t normalize = 1/hr1sims->Integral(0,hr1sims->FindBin(0.1)-1);
hr1sims->Scale(normalize);
hr1sims->SetLineWidth(2);
normalize = 1/hr1data->Integral(0,hr1data->FindBin(0.1)-1);
hr1data->Scale(normalize);
hr1data->SetLineColor(kRed);
hr1data->SetLineWidth(2);
hr1data->SetFillStyle(0);
cSims->cd(2);
hr1sims->Draw();
hr1data->Draw("same");

//Get the radial distribution
hr2sims->Reset();

TH1F *hr2dataOverlay = (TH1F*)hr2data->Clone("hr2dataOverlay");
hr2dataOverlay->SetFillStyle(0);
hr2dataOverlay->Reset();
for(unsigned int x = 0; x< hr2sims->GetNbinsX(); x++)
{

  //sims
  Double_t dIntegral = hr1sims->Integral(0,x);
  hr2sims->SetBinContent(x,dIntegral);
  dIntegral = hr1data->Integral(0,x);
  hr2dataOverlay->SetBinContent(x,dIntegral);
}
cSims->cd(3);
hr2sims->SetLineWidth(2);
hr2sims->Draw();
hr2dataOverlay->SetLineWidth(2);
hr2dataOverlay->SetLineColor(kRed);
hr2dataOverlay->Draw("same");

TLegend *legend = new TLegend(0.5,0.15,0.8,0.35);
legend->AddEntry(hr2sims,"simulation","l");
legend->AddEntry(hr2dataOverlay,"data","l");
legend->Draw();

}
