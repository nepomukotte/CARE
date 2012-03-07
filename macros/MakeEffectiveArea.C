void MakeEffectiveArea(string directory,double radius = 750,double zd = 20.)
{

  gStyle->SetOptStat(0);
  gStyle->SetFillColor(0);


 void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(directory.c_str()));
  const char* strptr;

  int fNumFiles = 0;

  if (!dir) {

     cout<<"Could not open directory: "<<directory.c_str()<<endl;

  }

 vector<string> vFiles;



  while ((strptr = gSystem->GetDirEntry(dir))) {
    //cout<<"next file"<<endl;
    string file = strptr;
    size_t found;
    found=file.find("root");
    if (found==string::npos) continue;
    vFiles.push_back(file);
    fNumFiles++;
  }


 cout<<"Have found "<<fNumFiles<<" files"<<endl;


 TH1D *hEvents = new TH1D("hEvents","",30,-1.5,1);
 hEvents->GetXaxis()->SetTitle("Energy [ log10(TeV) ]");
 hEvents->GetYaxis()->SetTitle("Effective Area [ m^{2} ]");
 hEvents->Sumw2();
 hEvents->SetMarkerStyle(23);
 TH1D *hTriggered = hEvents->Clone("hTriggered");
 TH1D *hEffArea = hEvents->Clone("hEffArea");

 TH1D *hRadialEventDistr = new TH1D("hRadialEventDistr","Radial Event distribution",200,0,1000);

 for(int i =0;i<vFiles.size();i++)
   {
     stringstream fnamestream;
     fnamestream<<directory<<"/"<<vFiles[i].c_str();

     cout<<"opening file: "<<fnamestream.str()<<endl;
 
     TFile *fCurrentFile = new TFile(fnamestream.str().c_str(),"READ");
     if( fCurrentFile->IsZombie() )
         {
           cout << "error opening root input file: " <<fnamestream.str().c_str()  << endl;
           cout << "...skipping" << endl;
           continue; //exit( -1 );
         }
  
     TTree *t = (TTree*)fCurrentFile->Get("Events/tSimulatedEvents" );

     if( !t )
       {
	 cout << "error: tree tSimulatedEvents not found in " << fnamestream.str().c_str()<< endl;
	 cout << "...skipping" << endl;
	 continue; //exit( -1 );
       }

     Float_t energy;
     t->SetBranchAddress("energy", &energy );
     Char_t arrayTriggerBit;
     t->SetBranchAddress("arrayTriggerBit", &arrayTriggerBit);
     Float_t xcore;
     t->SetBranchAddress("xcore", &xcore);
     Float_t ycore;
     t->SetBranchAddress("ycore", &ycore);

     //Looping over events in File
     for( int n = 0; n < t->GetEntries() ; n++ )
       {
	 
	 t->GetEntry(n);

	 hEvents->Fill(log10(energy));
	 
	 if( arrayTriggerBit == 1)
	   {
	     hTriggered->Fill(log10(energy));
	     Float_t r = sqrt(xcore*xcore+ycore*ycore);
	     //cout<<radius<<endl;
	     if(r!=0);
	     hRadialEventDistr->Fill(ycore);//,1/(TMath::Pi()*radius*radius))
	   }

       }//end looping over events in file

    fCurrentFile->Close();
    fCurrentFile->Delete();

   }//end looping over all files


 for(int i=1;i<=hEvents->GetNbinsX();i++)
   {

     Double_t dTriggered = hTriggered->GetBinContent(i);
     Double_t dNEvents = hEvents->GetBinContent(i);

     if(dNEvents!=0)
       {
	 hEffArea->SetBinContent(i,dTriggered/dNEvents*radius*radius*TMath::Pi());
	 hEffArea->SetBinError(i,sqrt(dTriggered)/dNEvents*radius*radius*TMath::Pi());
       }
     else
        hEffArea->SetBinContent(i,0);
   }


 TCanvas *c = new TCanvas("c","",800,500);
 c->SetLogy();
 c->SetGrid();
 c->Draw();
 hEffArea->Draw();

 TCanvas *cRadialDistribution = new TCanvas("cRadialDistribution","",800,500);
 cRadialDistribution->Draw();
 hRadialEventDistr->Draw();
}
