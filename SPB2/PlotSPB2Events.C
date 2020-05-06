//You need to make sure that this is the actual number of pixels in the camera
const int iNumPixels = 512; 
vector< vector<Int_t> *>   iFADCTraceInPixel;
vector<Int_t>   *iPEInPixel;
int iLastPix = -1;
TLatex *text = 0;
TCanvas *cDisplay = 0;
TH1F *hPixelTrace =0;

Bool_t HandleInput()
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

void ShowPixeltrace(int iPix)
{
   if(hPixelTrace==0)
     {
      hPixelTrace = new TH1F("hPixelTrace","Pixel Trace",512,-0.5,511.5);
      hPixelTrace->SetStats(0);
      hPixelTrace->GetXaxis()->SetTitle("ADC sample");
      hPixelTrace->GetYaxis()->SetTitle("ADC counts");
      
     }
   hPixelTrace->Clear();
   TString title;
   title.Form("Trace of Pixel %i",iPix);
   hPixelTrace->SetTitle(title);
   for(int s=0;s<iFADCTraceInPixel[iPix]->size();s++)       
     {
        hPixelTrace->SetBinContent(s+1,(iFADCTraceInPixel[iPix])->at(s));
     }
   cDisplay->cd(2);
   hPixelTrace->Draw();
   gPad->Modified();
   gPad->Update();
  
}

int FindPixel(int x, int y)
{
 //cout<<x<<"  "<<y<<endl;
 //column    
 int MUSIC_column = x/2;

 //cout<<"MUSIC column "<<MUSIC_column<<endl;
 //row
 int MUSIC_row = y/4;
 //cout<<"MUSIC row "<<MUSIC_row<<endl;

 //MUSIC ID
 int MUSIC_ID = MUSIC_column+MUSIC_row*16;
 //cout<<"MUSIC_ID "<<MUSIC_ID<<endl;

 //Channel in MUSIC
 int MUSIC_Channel = y%4+4*(x%2); 
 //cout<<"MUSIC_Channel "<<MUSIC_Channel<<endl;

 return MUSIC_row*8*16+MUSIC_column*8+MUSIC_Channel;
}

void ShowInfoAtCursor(int x, int y)
{
 //cout<<x<<"  "<<y<<endl;
 //column    
 int MUSIC_column = x/2;

 //cout<<"MUSIC column "<<MUSIC_column<<endl;
 //row
 int MUSIC_row = y/4;
 //cout<<"MUSIC row "<<MUSIC_row<<endl;

 //MUSIC ID
 int MUSIC_ID = MUSIC_column+MUSIC_row*16;
 //cout<<"MUSIC_ID "<<MUSIC_ID<<endl;

 //Channel in MUSIC
 int MUSIC_Channel = y%4+4*(x%2); 
 //cout<<"MUSIC_Channel "<<MUSIC_Channel<<endl;

 int PixID = MUSIC_row*8*16+MUSIC_column*8+MUSIC_Channel;

 //cout<<"MUSIC_ID: "<<MUSIC_ID<<" MUSIC_Channel: "<<MUSIC_Channel<<" Pixel ID: "<<PixID<<endl;
 TString statusline;
 statusline.Form("MUSIC_ID: %i MUSIC_Channel: %i Pixel ID: %i  PEs in Pixel: %i",MUSIC_ID,MUSIC_Channel,PixID,iPEInPixel->at(PixID));
 if(text!=0)
   text->Delete();
 TLatex T1;
 text = T1.DrawLatexNDC(0.1,0.95,statusline.Data());
 gPad->Modified();
 gPad->Update();
}

void FindBin(int iPix,int *nx, int *ny)
{
 //MUSIC ID
 int MUSIC_ID = iPix/8;

 //Channel in MUSIC
 int MUSIC_Channel = iPix%8; 

 //row
 int MUSIC_row = MUSIC_ID/16;
 
 //column    
 int MUSIC_column = MUSIC_ID%16;

 *ny = 4*MUSIC_row+MUSIC_Channel%4+1; 
 *nx = 2*MUSIC_column+MUSIC_Channel/4+1; 
 //cout<<"Pixel "<<iPix<<" x: "<<*nx<<" y: "<<*ny<<endl;

}

void PixelClicked() {
   //this action function is called whenever you move the mouse
   //it just prints the id of the picked triangle
   //you can add graphics actions instead
   int event = gPad->GetEvent();
   TObject *o = gPad->GetSelected();
   if (!o) return;
   if (!(o->InheritsFrom("TH2")))
       return;
   TH2F *h = (TH2F*)o;
   int px = gPad->GetEventX();
   int py = gPad->GetEventY();
   //cout<<px<<"  "<<py<<endl;
   Float_t xx = gPad->AbsPixeltoX(px);
   Float_t yy = gPad->AbsPixeltoY(py);
   Float_t x = 0.5+gPad->PadtoX(xx);
   Float_t y = 0.5+gPad->PadtoY(yy);
   //cout<<x<<"  "<<y<<endl;
   int pix = FindPixel((int)x,(int)y);
   if(pix!=iLastPix)
    {
      iLastPix=pix;
      ShowInfoAtCursor((int) x, (int) y);
    }

   if (event == 11){
      //cout<<pix<<endl;  
      ShowPixeltrace(pix);
   }
}

void DrawMUSICBoundaries()
{
  TBox *b = new TBox(-.5,-0.5,1.5,3.5);
  b->SetFillStyle(0);
  b->SetLineColor(kRed);
  b->Draw(); 
  TPoint p;
  for(int i=0;i<iNumPixels/8;i++)
   {
     TBox *bn = (TBox*)b->Clone();
     bn->SetX1((i%16)*2-0.5);
     bn->SetX2((i%16)*2+1.5);
     bn->SetY1((i/16)*4-0.5);
     bn->SetY2((i/16)*4+3.5);
     bn->Draw();
    
   } 
}

void PlotSPB2Events(string fInputFileName = "test.root")
{
       //Open the CARE file
       TFile *fO = new TFile( fInputFileName.c_str(), "READ" );
       if( fO->IsZombie() )
	 {
	   cout << "error opening root input file: " << fInputFileName << endl;
	   cout << "...exiting" << endl;
	   exit( -1 );
	 }

       cout<<"Have opened the file with the simulated events: "<<fInputFileName.c_str()<<endl;
      
       UInt_t uNumTelescopes;
       Float_t energy;
       UInt_t eventNumber ; 
       Float_t xcore ;
       Float_t ycore ;
       Float_t azPrim ;
       Float_t znPrim ;
       Bool_t arrayTriggerBit ;
       std::vector< Bool_t > *vTelescopeTriggerBits =0 ;
       Float_t DeltaTL3 ;

       // read tree with all the event by event information
       cout<<"Looking for Tree tSimulatedEvents"<<endl;
       TTree *tSimulatedEvents = (TTree*)fO->Get( "Events/tSimulatedEvents" );
	   if( !tSimulatedEvents )
	     {
	       cout << "error: tree tSimulatedEvents not found in " << fInputFileName << endl;
	       cout << "...exiting" << endl;
	       exit( -1 );
	     }

       tSimulatedEvents->SetBranchAddress("energy",&energy);
       tSimulatedEvents->SetBranchAddress("ZnPrim",&znPrim);
       tSimulatedEvents->SetBranchAddress("AzPrim",&azPrim);
       tSimulatedEvents->SetBranchAddress("xcore",&xcore);
       tSimulatedEvents->SetBranchAddress("ycore",&ycore);
       tSimulatedEvents->SetBranchAddress("arrayTriggerBit",&arrayTriggerBit);
       tSimulatedEvents->SetBranchAddress("uNumTelescopes",&uNumTelescopes);
       tSimulatedEvents->SetBranchAddress("eventNumber",&eventNumber);
       tSimulatedEvents->SetBranchAddress("DeltaTL3",&DeltaTL3);
       tSimulatedEvents->SetBranchAddress("vTelescopeTriggerBits",&vTelescopeTriggerBits);

       cout<<"Looking for tree with camera output"<<endl;
       TTree *T0 = (TTree*)fO->Get( "Events/T0" );
	   if( !T0 )
	     {
	       cout << "error: tree T0 not found in " << fInputFileName << endl;
	       cout << "...exiting" << endl;
	       exit( -1 );
	     }

       vector<int> *vTriggerCluster;
       vector<Float_t> *fTimeOverThreshold;
       vector<Float_t>   *fSumTimeInPixel; 
       vector<Int_t>   *iQDCInPixel;
       Int_t           iNumPhotonsInFocalPlane;
       Float_t         fAzTel ;
       Float_t         fZnTel ;
       vector<Bool_t>  *bInLoGain;
       iFADCTraceInPixel.assign(iNumPixels,0);

      T0->SetBranchAddress("vGroupsInTriggerCluster",&vTriggerCluster);
      T0->SetBranchAddress("vTimeOverThreshold",&fTimeOverThreshold);
      T0->SetBranchAddress("vSumTimeInPixel", &fSumTimeInPixel);
      T0->SetBranchAddress("vPEInPixel", &iPEInPixel);
      T0->SetBranchAddress("vQDCValue", &iQDCInPixel);
      T0->SetBranchAddress("iPhotonsInFocalPlane", &iNumPhotonsInFocalPlane);
      T0->SetBranchAddress("fAzTel", &fAzTel);
      T0->SetBranchAddress("fZnTel", &fZnTel);
      T0->SetBranchAddress("vHiLoGainBit", &bInLoGain);
      TString name;
      for(int g=0;g<iNumPixels;g++)       
         {                                                                                                     
            name.Form("vFADCTraces%i",g);         
            T0->SetBranchAddress(name,&(iFADCTraceInPixel[g]));
         }

      //Stuff for visualizing
      cDisplay = new TCanvas("cDisplay","Display",2000,500);
      cDisplay->Divide(2,1);
      cDisplay->cd(1);
      gPad->AddExec("ex","PixelClicked()");
      TH2F *hDisplay = new TH2F("hDisplay","Charge",32,-0.5,31.5,16,-0.5,15.5);
      hDisplay->SetStats(0);
      hDisplay->Draw("colz");
      DrawMUSICBoundaries(); 
       //looping over events doing something
       for(int n=0;n<tSimulatedEvents->GetEntries();n++)
         {
           tSimulatedEvents->GetEntry( n );
           if(arrayTriggerBit)
             {
                hDisplay->Clear();
                cout<<"Event "<<n<<" is triggered"<<endl; 
                T0->GetEntry(n);
                for(int g=0;g<iNumPixels;g++)       
                   {                                             
                
                    //cout<<"Pixel "<<g<<endl;
                    //cout<<"PEs: "<<iPEInPixel->at(g)<<endl;
                    //cout<<"QDC: "<<iQDCInPixel->at(g)<<endl;
                    int nx, ny;
                    FindBin(g,&nx,&ny);
                    hDisplay->SetBinContent(nx,ny,iPEInPixel->at(g));
                    //cout<<endl<<endl;
                   }
               cDisplay->Modified();
               cDisplay->Update();
               HandleInput();
              }
           }
 
}
