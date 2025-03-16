#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>

int plotCompareBiasCurves(std::string fileNameA, std::string fileNameB){
	TFile *f = new TFile(fileNameA.c_str(),"READ");
	TFile *fb = new TFile(fileNameB.c_str(),"READ");

	TGraph *bfBias = (TGraph*)f->Get("biasTel");
	TGraph *pxBias = (TGraph*)f->Get("biasPx");
	TGraph *sfBias = (TGraph*)fb->Get("biasTel");
	TLine *l = new TLine(1,10,20,10);
	TCanvas *c = new TCanvas ("c","",800,500);
	TH1 *h = new TH1F("h","",100,1,20);
	h->Draw();
	h->SetMinimum(1);
	h->SetMaximum(5e6);
	//bfBias->GetHistogram()->Draw();
	h->GetXaxis()->SetTitle("Threshold [pe]");
	h->GetXaxis()->SetRangeUser(1,20);
	h->GetYaxis()->SetTitle("Trigger Rate [Hz]");
	bfBias->GetHistogram()->SetMinimum(1);
	c->SetLogy();
	bfBias->SetMarkerStyle(29);
	bfBias->SetMarkerSize(1.6);
	bfBias->SetLineColor(kBlack);
	bfBias->SetMarkerColor(kBlack);
	bfBias->SetLineWidth(3);
	bfBias->Draw("LP");

        pxBias->SetMarkerStyle(29);
        pxBias->SetMarkerSize(1.6);
        pxBias->SetLineColor(kRed);
        pxBias->SetMarkerColor(kRed);
        pxBias->SetLineWidth(3);
        pxBias->Draw("LP");
	
	sfBias->GetHistogram()->GetXaxis()->SetRangeUser(1,20);
        sfBias->SetMarkerStyle(29);
        sfBias->SetMarkerSize(1.6);
        sfBias->SetLineColor(kBlue);
        sfBias->SetMarkerColor(kBlue);
        sfBias->SetLineWidth(3);
        sfBias->Draw("LP");
	
	auto legend = new TLegend(0.7,0.7,0.9,0.9);
	legend->AddEntry(bfBias,"BiFocal Telescope Rate","lp");
	legend->AddEntry(pxBias,"Single Pixel Rate","lp");
	legend->AddEntry(sfBias,"Single Focus Telescope Rate","lp");
	legend->Draw();
	
	l->SetLineStyle(2);
	l->SetLineWidth(2);
	l->Draw();
	return 0;
}
