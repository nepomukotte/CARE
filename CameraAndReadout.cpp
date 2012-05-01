/*! \file CameraAndReadout.cpp
    \Main program of CARE a camera simulation for Cherenkov telescopes

    \author Nepomuk Otte
*/

#include <iostream>
#include <string>
#include <vector>

#include "TriggerTelescopeNextNeighbor.h"
#include "ArrayTrigger.h"
#include "ReadConfig.h"
#include "TraceGenerator.h"
#include "FADC.h"
#include "Display.h"
#include "TelescopeData.h"

enum SimulationPackageCodes   {KNOWNNOT,LEEDS,GRISU,KASCADE,CORSIKA,UCLA};
#include "VG_writeVBF.h"



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

using namespace std;


static const Bool_t DEBUG_TRACE = kFALSE;
static const Bool_t DEBUG_FADC = kFALSE;
static const Bool_t DEBUG_TELTRIGGER = kFALSE;
static const Bool_t DEBUG_ARRAYTRIGGER = kFALSE;
static const Bool_t DEBUG_MAIN = kFALSE;
static const Bool_t DEBUG_TELESCOPEDATA = kFALSE;


void getRotatedShowerDirection( double ze, double az, double y, double x, double &rze, double &raz )
{

    double degrad = 57.29577;

// get all directions in [rad]
    x /= degrad;
    y /= (-1*degrad);

// assume all telescopes point in same directions
    double el = (90.-ze)/degrad;
    az = az/degrad;
// these are the resulting directions

    double r = sqrt( 1. + x*x + y*y );
    double cx = x/r;
    double cy = 1./r;
    double cz = y/r;

// rotate telescope around elevation axis
    double ex = cx;
    double ey = cy * cos( el ) - cz * sin( el );
    double ez = cy * sin( el ) + cz * cos( el );
// rotate around azimuth
    double rx, ry, rz;
    rx =     ex * cos( az ) + ey * sin( az );
    ry = -1.*ex * sin( az ) + ey * cos( az );
    rz = ez;
// calculate new azimuth, zenith
    r = sqrt( rx*rx + ry*ry );
// small value check
    if( fabs(r) < 1.e-10 ) r = 0.;
    if( fabs(rx) < 1.e-10 ) rx = 0.;
    if( fabs(ry) < 1.e-10 ) ry = 0.;
    if( fabs(rz) < 1.e-10 ) rz = 0.;

    if( r == 0. ) raz = az * degrad;
    else
    {
        raz = (TMath::Pi()/2.-atan2(ry,rx))*degrad;
        if( raz > 180. )  raz = -1.*(360.-raz);
        if( raz < -180. ) raz *= -1.;
    }
    if( rz == 0. ) rze = 90. - el*degrad;
    else
    {
        rze = 90.-atan2(rz,r)*degrad;
    }
}


void help()
{
     cout << "CameraAndReadout <rng seed> <Care cfg file> <basename output files> <pe input file> <VBF run number>" << endl;
     cout << endl;
     cout << "command line options: " << endl;
     cout << "\t <rng seed> \t Seed of the random number generator" << endl;
     cout << "\t <Care cfg file> \t configuration file of the Camera and electronics" << endl;
     cout << "\t <basename output files> \t the base name for the vbf and root output files" << endl;
     cout << "\t <pe input file> \t the input pe file" << endl;
     cout << "\t <VBF run number> \t the runnumber written into the vbf file" << endl;
     cout << "\t <Pedestal flag> \t if 1 the pedestal events will be written into the output file" << endl;
     cout << endl;

    exit( 0 );
}

int main( int argc, char **argv )
{
   if( argc < 3 ) help();
   
   //configuration File
   
   UInt_t uSeed = (UInt_t)(atoi( argv[1] ) );
   string sConfigFileName = argv[2]; 
   string sOutputFileName = argv[3];


   string fInputFileName;
   if(argc >= 5)
     fInputFileName = argv[4];
  
   long lVBFRunNum;
   if(argc >= 6)
     lVBFRunNum  = (long)(atoi( argv[5] ) );
   else
     lVBFRunNum = 99999;
   
   Int_t iPedestalWriteFlag = 0;
   if(argc >=7)
	  iPedestalWriteFlag = (int)(atoi( argv[6] ) );

   Int_t iDebugLevel=0;
    if(argc == 8)
      iDebugLevel  = atoi( argv[7] );

   gROOT->ProcessLine("#include <vector>");
   gStyle->SetPalette(1);
   gStyle->SetOptStat(0);
   
   ///////////////////////////////////////////
   //
   //Initialize Random Number Generator
   //
   //////////////////////////////////////////

   TRandom3 *rand = new TRandom3(uSeed);


   ///////////////////////////////////////////////////////////////////////////////////////////////
   // 
   // Getting Array and Camera information
   //
   ///////////////////////////////////////////////////////////////////////////////////////////////


   ReadConfig *readConfig = new ReadConfig(rand);
   readConfig->ReadConfigFile(sConfigFileName);

   ///////////////////////////////////////////////////////
   //
   // Debug Display
   //
   //////////////////////////////////////////////////////


   Display *display = NULL;
   if(DEBUG_TRACE || DEBUG_FADC || DEBUG_TELTRIGGER || DEBUG_ARRAYTRIGGER || DEBUG_MAIN || DEBUG_TELESCOPEDATA)
      display = new Display(argc,argv,readConfig);




   //////////////////////
   //  Array Information
   /////////////////////

   //Get the information about the arrays
   //coordinates etc...
   //Get the subbarray that we want to study

   //The number of telescopes in the Subarray
   UInt_t uNumTelescopes = readConfig->GetNumberOfTelescopes();

   cout<<"This array has "<<uNumTelescopes<<" telescopes"<<endl;

   UInt_t uNumTelescopeTypes = readConfig->GetNumberOfTelescopeTypes(); //number of telescope types for CTA array

   cout<<"This array has "<<uNumTelescopeTypes<<" telescope types"<<endl;

   ////////////////////////////////////////////////////
   //
   //Initialize Trace Generator
   //
   ///////////////////////////////////////////////////


   //needs to be changed for telescope types
   TraceGenerator *traceGenerator[uNumTelescopeTypes];

   for(UInt_t t = 0 ; t<uNumTelescopeTypes; t++)
     {

       cout<<endl<<"Starting TraceGenerator for telescope"<<t<<endl;
       traceGenerator[t] = new TraceGenerator(readConfig,t,rand,DEBUG_TRACE,display);
  
       //traceGenerator->ShowAverageSinglePEPulse();

       cout<<"Telescope tracegenerator is initialized"<<endl;



     }

   ////////////////////////////////////////////////////
   //
   //Initialize Telescope Trigger
   //
   ///////////////////////////////////////////////////


   //needs to be changed for telescope types
   TriggerTelescopeNextNeighbor *Teltrigger[uNumTelescopeTypes]; 
   for(UInt_t t = 0 ; t<uNumTelescopeTypes; t++)
     {
	 
       cout<<endl<<"Starting Trigger for telescopetype"<<t<<endl;
       Teltrigger[t] = new TriggerTelescopeNextNeighbor(readConfig, t, rand, DEBUG_TELTRIGGER, display ); 
       
       cout<<"Telescope trigger is initialized"<<endl;
     }

   ////////////////////////////////////////////////////
   //
   //Initialize Array Trigger
   //
   ///////////////////////////////////////////////////

   ArrayTrigger *arraytrigger = new ArrayTrigger(readConfig,DEBUG_ARRAYTRIGGER);

   ///////////////////////////////////////////////////////
   //
   // Initialize FADCs
   //
   //////////////////////////////////////////////////////
   
   FADC *fadc[uNumTelescopeTypes];

   for(UInt_t t = 0 ; t<uNumTelescopeTypes; t++)
     {
       cout<<endl<<"starting FADC for telescope "<<t<<endl;

       fadc[t] = new FADC(readConfig,traceGenerator[t], t, rand,DEBUG_FADC, display);

       cout<<"FADC is initialized"<<endl;
     }

   ////////////////////////////////////////////////////////
   //
   // Create the Telescope Storage containers
   //
   ///////////////////////////////////////////////////////


   TelescopeData *telData[uNumTelescopes];

   for(UInt_t t = 0; t<uNumTelescopes; t++)
    {
      cout<<endl<<"creating Telescope Data container for telescope id: "<<t<<endl;
      //Note that argument contains telescope type and number pixels. Needs to be done cleaner with config reader
      telData[t] = new TelescopeData(readConfig,t,DEBUG_TELESCOPEDATA);
    }
	//Setting the telescope data for the debug display
   
   if(display)
    { 
      display->SetTelescopeData(telData);
    }

   ////////////////////////////////////////////////////////
   //
   // Open root file for output
   //
   ////////////////////////////////////////////////////////

   TString outname =  sOutputFileName.c_str();
     outname+=".root";
    
	TFile *fOut = new TFile( outname.Data(),"RECREATE");


   if( !fOut->IsOpen() )
	 {
	   cout << "error opening root output file: " << outname << endl;
	   cout << "...exiting" << endl;
	   exit( -1 );
	 }


       cout<<"Have opened the root outputfile: "<<outname<<endl;


       //Here comes the Trigger and array configuration used in the simulation
       fOut->cd();
       fOut->mkdir("TriggerArrayConfiguration");
       fOut->cd("TriggerArrayConfiguration");

       //change back to the root root directory
       gROOT->cd();


   ////////////////////////////////////////////////////////
   //
   //   Bias curve 
   //
   /////////////////////////////////////////////////////////


   if(readConfig->GetBiasCurveBit())
     {
       cout<<"Run a bias curve"<<endl;

       //move to Trigger class

       UInt_t trials = readConfig->GetNumberOfTrialsForBiasCurve();

       if(trials<=0)
	 {
	   cout<<"The number of trials for the bias curve is set wrong "<<trials<<endl;
	   exit(1);
	 }

       Float_t start = readConfig->GetBiasCurveStartScanRange();

       Float_t stop = readConfig->GetBiasCurveStopScanRange();

       Float_t step = readConfig->GetBiasCurveStep();
       if(step==0)
	 {
	   cout<<"The step in the scan of the discriminator values is set wrong "<<step<<endl;
	   exit(1);
	 }

      
      UInt_t TelID = readConfig->GetBiasCurveTelescopeID(); 

       if(TelID<0)
     {
       cout<<"The telescope ID has to be >= 0 but is "<<TelID<<endl;
       exit(1);
     }

      UInt_t TelType = readConfig->GetTelescopeType(TelID); //need to get this out of the configuration file if there are more than one telescope type

	  traceGenerator[TelType]->SetTelData(telData[TelID]);
//	  Teltrigger[TelType]->LoadEvent(telData[TelID]); //necessary in order to pass the telData container
      Teltrigger[TelType]->RunBiasCurve(trials,start,stop,step,traceGenerator[TelType],telData[TelID]); 


      //reset discriminator values to the ones in the config file for going over the events
      Teltrigger[TelType]->SetDiscriminatorThresholdAndWidth(readConfig->GetDiscriminatorThreshold(TelType),
						     readConfig->GetDiscriminatorOutputWidth(TelType));

       fOut->cd();
       fOut->mkdir("BiasCurve");
       fOut->cd("BiasCurve");

       Teltrigger[TelType]->GetBiasCurve().Write("TVecBiasCurve");
       Teltrigger[TelType]->GetBiasCurveError().Write("TVecBiasCurveError");
       Teltrigger[TelType]->GetBiasCurveScanPoints().Write("TVecBiasCurveScanPoints");
       Teltrigger[TelType]->GetGroupRateVsThresholdError().Write("TVecGroupRateVsThresholdError");
       Teltrigger[TelType]->GetGroupRateVsThreshold().Write("TVecGroupRateVsThreshold");
       TVectorF NumTelescopes(1);
       NumTelescopes[0]= uNumTelescopes;
       NumTelescopes.Write("TVecNumTelescopes");

       //change back to the root root directory
       gROOT->cd();

     }



   /////////////////////////////////////////////////////////////////////////////
   // 
   //    Going into the simulated showers
   //
   //////////////////////////////////////////////////////////////////////////////


   if(readConfig->GetLoopOverEventsBit())
     {

       
       VG_writeVBF *VBFwrite = NULL;

       if(readConfig->GetVBFwriteBit())
	 {
	   //this is probably only needed to get something into the constructor. 
	   //the number of pixel are set again when the traces are being filled.
       vector<unsigned> numPixTel(uNumTelescopes,readConfig->GetNumberPixels(0)); //this will be problematic when simulating hybrid arrays.
       
       //Initiate the vbf output file
       string outname =  sOutputFileName.c_str();
       outname+=".vbf";

	   VBFwrite = new VG_writeVBF(uNumTelescopes,
				      numPixTel,outname,
				      lVBFRunNum,
				      iDebugLevel);

	 }
   
       //initiate the root output file
       fOut->cd();
       fOut->mkdir("Events");
       fOut->cd("Events");

       TTree tSimulatedEvents("tSimulatedEvents","Tree that holds the output after trigger");
       Float_t energy;
       UInt_t eventNumber ; 
       Float_t xcore ;
       Float_t ycore ;
       Bool_t arrayTriggerBit ;
       std::vector< Bool_t > vTelescopeTriggerBits ;
       vTelescopeTriggerBits.assign(uNumTelescopes , 0 ) ;

       std::vector< int > vGroupsInTriggerCluster[uNumTelescopes];
       
       Float_t DeltaTL3 ;

       tSimulatedEvents.Branch("energy",&energy,"energy/F");
       tSimulatedEvents.Branch("xcore",&xcore,"xcore/F");
       tSimulatedEvents.Branch("ycore",&ycore,"ycore/F");
       tSimulatedEvents.Branch("arrayTriggerBit",&arrayTriggerBit,"arrayTriggerBit/B");
       tSimulatedEvents.Branch("uNumTelescopes",&uNumTelescopes,"uNumTelescopes/l");
       tSimulatedEvents.Branch("eventNumber",&eventNumber,"eventNumber/l");
       tSimulatedEvents.Branch("vTelescopeTriggerBits",&vTelescopeTriggerBits);
       tSimulatedEvents.Branch("DeltaTL3",&DeltaTL3,"DeltaTL3/F");

       for(UInt_t i = 0; i<uNumTelescopes; i++)
	 {
	   TString namevar;
	   namevar.Form("vGroupsInTriggerClusterTel%i",i);
	   tSimulatedEvents.Branch(namevar.Data(),&(vGroupsInTriggerCluster[i]));
	 }

       //Open the photon input file
       TFile *fO = new TFile( fInputFileName.c_str(), "READ" );
       if( fO->IsZombie() )
	 {
	   cout << "error opening root input file: " << fInputFileName << endl;
	   cout << "...exiting" << endl;
	   exit( -1 );
	 }

       cout<<"Have opened the file with the simulated events: "<<fInputFileName.c_str()<<endl;
       
       // read tree with general infos
	   cout<<"Looking for Tree allT"<<endl;
	   TTree *tGeneralInfo = (TTree*)fO->Get( "allT" );
	   if( !tGeneralInfo )
	     {
	       cout << "error: tree allT not found in " << fInputFileName << endl;
	       cout << "...exiting" << endl;
	       exit( -1 );
	     }
        Double_t dObsHeight;
        Double_t dGlobalPhotonEffic;
	    string *fileheader = 0;
        std::vector<int>  *telIDVector = 0;
		std::vector<float> *telLocXGCVector = 0;
		std::vector<float> *telLocYGCVector = 0;
		std::vector<float> *telLocZGCVector = 0;
		std::vector<float> *transitTimeVector = 0;
    
        TBranch *b_telIDVector;
        TBranch *b_telLocXGCVector;
        TBranch *b_telLocYGCVector;
        TBranch *b_telLocZGCVector;
        TBranch *b_transitTimeVector;
        TBranch *b_fileheader;

        tGeneralInfo->SetBranchAddress("telIDVector",&telIDVector,&b_telIDVector);
        tGeneralInfo->SetBranchAddress("telLocXVector",&telLocXGCVector,&b_telLocXGCVector);
		tGeneralInfo->SetBranchAddress("telLocYVector",&telLocYGCVector,&b_telLocYGCVector);
		tGeneralInfo->SetBranchAddress("telLocZVector",&telLocZGCVector,&b_telLocZGCVector);
		tGeneralInfo->SetBranchAddress("transitTimeVector",&transitTimeVector,&b_transitTimeVector);
		tGeneralInfo->SetBranchAddress("obsHgt", &dObsHeight );
        tGeneralInfo->SetBranchAddress("globalEffic", &dGlobalPhotonEffic );
        tGeneralInfo->SetBranchAddress("fileHeader", &fileheader,&b_fileheader );
        tGeneralInfo->GetEntry( 0 );
        cout<<"Observatory Height "<<dObsHeight<<endl;
        cout<<"Global Photon Efficiency "<<dGlobalPhotonEffic<<endl;
        cout<<"File header "<<fileheader->c_str()<<endl; 

        if((unsigned)readConfig->GetNumberOfTelescopes() != telIDVector->size())
		{
			cout<<"The number of telescopes in the input file does not match the number of telescopes in the configuration file"<<endl;
			exit(1);
		}

         vector<float> vTelTransitTimes;
        //loop over telescopes in CARE
		for(int c =0 ; c<readConfig->GetNumberOfTelescopes() ; c++)
		   {
               //loop over telescopes in GrOptics file to find the right match
               for(unsigned i = 0;i<telIDVector->size();i++)
            	   {
		              int iGrOpticsTelID = telIDVector->at(i);
                      if(readConfig->GetTelescopeIDinSuperArray(c) == iGrOpticsTelID)
			              {
				             telData[c]->TelXpos = telLocXGCVector->at(i);
				             telData[c]->TelYpos = telLocYGCVector->at(i);
				             telData[c]->TelZpos = telLocZGCVector->at(i);
							 cout<<"Telescope "<<c<<" location x: "<< telData[c]->TelXpos<<" y: "<<telData[c]->TelYpos<<" z: "<<telData[c]->TelZpos<<endl; 
							 telData[c]->OpticsTransitTime = transitTimeVector->at(i);
                             cout<<"transit time throught the telescope optics "<<telData[c]->OpticsTransitTime<<endl;
                             vTelTransitTimes.push_back(transitTimeVector->at(i)); 
                             break;
						  }
				    }     
			   }
        arraytrigger->SetInterTelTransitTimes(vTelTransitTimes);
       // read trees from file one for each telescope
       TTree **t = new TTree*[uNumTelescopes];
       for(UInt_t i = 0; i<uNumTelescopes; i++)
	 {
	   UInt_t uTrueTelID = readConfig->GetTelescopeIDinSuperArray(i); 
	   char hname[400];
	   sprintf( hname, "T%d", uTrueTelID );
	   cout<<"Looking for Tree "<<hname<<endl;
	   t[i] = (TTree*)fO->Get( hname );
	   if( !t[i] )
	     {
	       cout << "error: tree " << hname << " not found in " << fInputFileName << endl;
	       cout << "...exiting" << endl;
	       exit( -1 );
	     }
	 }

       UInt_t fEventNumber = 0;
       float fPrimaryEnergy = 0.;
       UInt_t iPrimaryType = 0;
       float fXcore = 0.;
       float fYcore = 0.;
       float fXcos = 0.;
       float fYcos = 0.;
       float fXsource = 0.;
       float fYsource = 0.;
       float fDelay = 0.;
       std::vector< float > *v_f_x = 0;
       std::vector< float > *v_f_y = 0;
       std::vector< float > *v_f_time = 0;
       std::vector< float > *v_f_lambda = 0;
       TBranch *b_v_f_x;
       TBranch *b_v_f_y;
       TBranch *b_v_f_time = 0;
       TBranch *b_v_f_lambda = 0;
       
       cout << "total number of entries: " << t[0]->GetEntries() << endl;
     

       //Looping over the events
       vector<Float_t>fTelTriggerTimes;
       fTelTriggerTimes.assign(uNumTelescopes,0);
       
       Int_t NumTriggeredEvents = 0 ;
       Int_t NumSkippeddEvents = 0 ;

       vector<Bool_t> *GroupTriggerBits = new vector<Bool_t>[uNumTelescopes];


       //Write the Simulation Header and do the pedestal events

       //we need this information for the VBF file simulation header and pedestals
       t[0]->SetBranchAddress("primaryType", &iPrimaryType );
       t[0]->SetBranchAddress("Xcos", &fXcos );
       t[0]->SetBranchAddress("Ycos", &fYcos );
       t[0]->GetEntry( 0 );

       if(readConfig->GetVBFwriteBit())
	 {
	   cout<<endl<<"Write the simulation header "<<endl;
	   
	   // set particle IDtype for both writers
	   VBFwrite->setParticleIDType(iPrimaryType);
		       
	   // create simulation headers
	   
	   //telescope and pixel locations
	   for (UInt_t tel = 0; tel < uNumTelescopes; tel++) {
	     float telSouth = -telData[tel]->TelYpos;
	     float telEast  = telData[tel]->TelXpos;
	     float telUp    = telData[tel]->TelZpos;
	     VBFwrite->setTelescopeLocations(tel,
					     telSouth,
					     telEast,
					     telUp);
	 
	     //needs to be fixed for different telescope types 
	     for (UInt_t pix = 0; pix < readConfig->GetNumberPixels(0);pix++) {
	       float pixEast = readConfig->GetXUnrotated(0)[pix];  
	       float pixUp = readConfig->GetYUnrotated(0)[pix];
	       float pixRad = readConfig->GetTubeSize(0)[pix];
	       VBFwrite->setPixelLocations(tel,pix,
					   pixEast,pixUp,pixRad);
	     } 
	   }
 
	 
	   u_int32_t simulationPackage = GRISU;
		   
	   std::string simulator = readConfig->GetSimulatorName().Data();

       string sSimInfoForHeader = sConfigFileName;
	   sSimInfoForHeader.append(*fileheader); 

       cout<<"Dumping the header written into the VBF file"<<endl;
       cout<<sSimInfoForHeader.c_str()<<endl;

	   VBFwrite->makeSimulationHeader( sSimInfoForHeader,
					   simulator,
					   dObsHeight,
					   simulationPackage,
					   readConfig->GetAtmosphericModel());


	   // set pedestal parameters for data writeVBF object
	   bool makePedsFlag = false;
	   if( iPedestalWriteFlag == 1)
	     makePedsFlag = true;
	   int pedTimeInterval = 1; //sec
	   std::string pedFileName("");

	   VBFwrite->setPedestalParameters(makePedsFlag,
					   pedTimeInterval,
					   pedFileName);

	   //std::cout << " setFirstEventTimeString for data" << std::endl;
		   
	   std::string fFirstValidEventTimeStr;
	   fFirstValidEventTimeStr = readConfig->GetDayOfSimulatedEvents().Data();
	   fFirstValidEventTimeStr += " 02:00:00.000000000 GPS";
	   cout<<"First event time stamp: "<<fFirstValidEventTimeStr<<endl;
	   VBFwrite->setFirstEventTimeString(fFirstValidEventTimeStr);
		   
	   double EventRatePerSec = 100;
	   VBFwrite->setShowerFrequencyHz(EventRatePerSec);
		   
	   // makeEventTimeClasses for data VBF writer
	   VBFwrite->makeEventTimeClasses();
	 } //That is it with the simulation header   

     
       cout<<"Doing the pedestals"<<endl;
       //Running the pedestal events
	   Int_t NumPedestalsToSimulate = readConfig->GetNumberOfPedestalEvents() 
	                                   + readConfig->GetNumberOfPedestalEventsToStabilize();
       
	   for(Int_t p = 0 ;p<  NumPedestalsToSimulate;p++)
	 {
	   if(p%100==0)
	     cout<<endl<<"Done "<<p<<" pedestal events"<<endl;

	   for (UInt_t tel=0;tel<uNumTelescopes;tel++){
		 
	     //generate traces with trace generator
         Int_t telType = telData[tel]->GetTelescopeType();

	   if(p%100==0)
	     cout<<"RFB:"<<Teltrigger[telType]->GetDiscRFBDynamicValue()<<endl;

         telData[tel]->ResetTraces();
         traceGenerator[telType]->SetTelData(telData[tel]);
	     traceGenerator[telType]->GenerateNSB();
		 
	     //Load event with trace from trace generator.
	     Teltrigger[telType]->LoadEvent(telData[tel]);
	     Teltrigger[telType]->RunTrigger();
         
	     //Running the FADC
	     telData[tel]->fTriggerTime = telData[tel]->fAveragePhotonArrivalTime ;
         fadc[telType]->RunFADC(telData[tel]);
	   }

	   //save to VBF file if we want to do that
	   if(readConfig->GetVBFwriteBit() 
		         && p>=readConfig->GetNumberOfPedestalEventsToStabilize() 
				 && iPedestalWriteFlag==1)
	     {
	       VBFwrite->setPedestalEvent();  // tell the writer this is peds packet
	       for (UInt_t tel1=0;tel1<uNumTelescopes;tel1++) {

		 Float_t Az = atan2( fXcos, fYcos );
		 Float_t Ze = 1. - (fXcos*fXcos+fYcos*fYcos);
		 if( Ze < 0. ) Ze = 0.;  // should never happen
		 Ze = acos( sqrt( Ze ) );

		 VBFwrite->setAzimElevTelDeg(tel1,Az*57.29577,90.0-Ze*57.29577);//az, elev in deg
        
	       }
	       VBFwrite->makePacket(); // make packet,
	       
	       VBFwrite->makeArrayEvent();  // make array event
	       
	       //pedVBF->setDebugLevel(1);
	       VBFwrite->makeArrayTrigger(); // number of triggering telescopes = all
	       // telescopes, 
	       
	       for (UInt_t tel=0;tel<uNumTelescopes;tel++){
        
		 VBFwrite->setTelescope(tel);  // set telescope number, starting from zero
		 
         // these two methods will set the numFadcsamples and the MaxNumberChannels which
         //    I also set to the MaxNumber
        VBFwrite->setMaxNumberChannels(telData[tel]->iNumPixels); // set channelNum for this telescope
        VBFwrite->setNumFadcSamples( telData[tel]->iNumFADCSamples);  // set numFadcSamples for this telescope

		 VBFwrite->makeEvent();  // make event for this telescope
       
		 for(UInt_t pix=0;pix<readConfig->GetNumberPixels( telData[tel]->GetTelescopeType() );pix++){
		   
		   VBFwrite->setPixel(pix);  // set pixel number

		   // set charge and pedestal to zero here
		   VBFwrite->setChargePedHigain();
		   /*A starting point is randomly selected along the 
		     simulated noise records*/
			   
          
		   /*Write the FADC trace out*/
		   vector<Int_t> trace = telData[tel]->GetFADCTrace(pix);
		   for(Int_t i=0;i<(Int_t)trace.size();i++){
		     VBFwrite->storeSample(i,trace[i]);  // store the sample
		   }
		 } // end pixel loop
		 
		 VBFwrite->storeEvent();
        
	       }   // end telescope loop
  
	       VBFwrite->storeArrayTrigger();
		       
	       VBFwrite->makeVSimulationData();  // take the parameter defaults 
	       VBFwrite->storeArrayEvent();      // store the arrayevent
		       
	       VBFwrite->storePacket();          // store the packet
	     }//end writing the pedestal event into the vbf file
	 }//end of loop generating all pedestal events

       cout<<"Done generating all the pedestal events"<<endl;
       
       ///////////////////////////////////////////////////////////////////////////////////////////

       //Going into the events
       for( int i = 0; i < t[0]->GetEntries() ; i++ )
	 {
	   if(DEBUG_MAIN)
	   cout<<endl<<endl<<endl<<"Event "<<i<<endl;
	   else
	     if(i%1000==0) cout<<endl<<"Event "<<i<<endl;

	   bool isArrayTriggered=false;

	   //Check if the minimum number of photons arrive in the focal plane if not skip the event
	   Bool_t bSkipEvent = kTRUE;
           for(UInt_t n = 0; n<uNumTelescopes; n++)
	     {   
	       if(DEBUG_MAIN)
	       cout<<"Telescope "<<n<<endl;
	       t[n]->SetBranchAddress("time", &v_f_time, &b_v_f_time );
	       t[n]->SetBranchAddress("eventNumber", &fEventNumber );
	       t[n]->SetBranchAddress("primaryEnergy", &fPrimaryEnergy ); 
	       t[n]->SetBranchAddress("Xcore", &fXcore );
	       t[n]->SetBranchAddress("Ycore", &fYcore );
	       t[n]->SetBranchAddress("Xcos", &fXcos );
	       t[n]->SetBranchAddress("Ycos", &fYcos );
	       t[n]->SetBranchAddress("Xsource", &fXsource );
	       t[n]->SetBranchAddress("Ysource", &fYsource );
       	       t[n]->GetEntry( i );

	       //General things we want to have in the root output file characterizing the event
	       energy = fPrimaryEnergy;
	       eventNumber = fEventNumber;
	       xcore = fXcore;
	       ycore = fYcore;

           Int_t telType = telData[n]->GetTelescopeType();

	       if(v_f_time->size()>readConfig->GetRequestedMinNumberOfPhotonsInCamera(telType))
		 bSkipEvent = kFALSE;

	       if(n == 0 && DEBUG_MAIN )
		 {
		   cout << "entry " << i << "\t, event number: " << fEventNumber << ", primary energy [TeV]: " << fPrimaryEnergy << endl;
		   cout << "\t core position [m]: " << fXcore << ", " << fYcore << endl;
		   cout << "\t source direction and offsets: " << fXcos << "\t" << fYcos << "\t" << fXsource << "\t" << fYsource << endl;
		 }
	     }//end checking if the event is skipped
	 
	   if(bSkipEvent)
	     {
	       NumSkippeddEvents++;
	       arrayTriggerBit = kFALSE;

	       if(DEBUG_MAIN)
		 {
		   cout<<"no photons in focal plane "<<NumSkippeddEvents<<"  "<<i<<endl;
		 }	
	     }
	   else
	     {
	       //Loop over the telescopes and see if they have triggered
	       for(UInt_t n = 0; n<uNumTelescopes; n++)
		 {   
	       
		   if(DEBUG_MAIN)
		     cout<<"Telescope "<<n<<endl;
		   t[n]->SetBranchAddress("eventNumber", &fEventNumber );
		   t[n]->SetBranchAddress("primaryEnergy", &fPrimaryEnergy );
		   t[n]->SetBranchAddress("primaryType", &iPrimaryType );
		   t[n]->SetBranchAddress("delay", &fDelay );
		   t[n]->SetBranchAddress("Xcore", &fXcore );
		   t[n]->SetBranchAddress("Ycore", &fYcore );
		   t[n]->SetBranchAddress("Xcos", &fXcos );
		   t[n]->SetBranchAddress("Ycos", &fYcos );
		   t[n]->SetBranchAddress("Xsource", &fXsource );
		   t[n]->SetBranchAddress("Ysource", &fYsource );
		   t[n]->SetBranchAddress("photonX", &v_f_x, &b_v_f_x ); 
		   t[n]->SetBranchAddress("photonY", &v_f_y, &b_v_f_y ); 
		   t[n]->SetBranchAddress("time", &v_f_time, &b_v_f_time );
		   t[n]->SetBranchAddress("wavelength", &v_f_lambda, &b_v_f_lambda );

		   t[n]->GetEntry( i );
		
		if( v_f_time->size() != v_f_x->size() )
		     {
		       
		       cout<<"Vectors do not have the same size, should never happen, isn't it?"<<endl;
		     }
		   

           //generate traces with trace generator
           Int_t telType = telData[n]->GetTelescopeType();

           telData[n]->ResetTraces();
           traceGenerator[telType]->SetTelData(telData[n]);

		   traceGenerator[telType]->LoadCherenkovPhotons( v_f_x ,  v_f_y, v_f_time, v_f_lambda, fDelay,dGlobalPhotonEffic);	
		   
		   //   TelData->ShowTrace(0,kTRUE); 
           		   

		   //run the telescope trigger
	       if(DEBUG_MAIN)
		   cout<<"Load the traces into trigger"<<endl;
           Teltrigger[telType]->LoadEvent(telData[n]);
 
	       if(DEBUG_MAIN)
		   cout<<"run trigger"<<endl;
           Teltrigger[telType]->RunTrigger();

		   // Teltrigger->ShowTrace(0,0);

	       if(DEBUG_MAIN)
		   cout<<"done trigger, RFB:"<<Teltrigger[telType]->GetDiscRFBDynamicValue()<<endl;

           fTelTriggerTimes[n] = telData[n]->GetTelescopeTriggerTime(); 
		   vGroupsInTriggerCluster[n] = telData[n]->GetTriggerCluster();
		   GroupTriggerBits[n] =  telData[n]->GetTriggeredGroups();
		   vTelescopeTriggerBits[n]=telData[n]->GetTelescopeTrigger();
		 } //end looping over telescopes doing trigger
       
     
	       //Go into the array trigger                       
	       if(DEBUG_MAIN)
		 {
		   cout<<endl<<"Event trigger status"<<endl;
		   cout<<"Telescope trigger bits: "<<vTelescopeTriggerBits[0]<<vTelescopeTriggerBits[1]<<vTelescopeTriggerBits[2]<<vTelescopeTriggerBits[3]<<endl;
		   cout<<"Telescope trigger times: "<<fTelTriggerTimes[0]<<"  "<<fTelTriggerTimes[1]<<"  "<<fTelTriggerTimes[2]<<"  "<<fTelTriggerTimes[3]<<"  "<<endl;
		 }
	   
	       arraytrigger->SetTelescopeTriggerBitsAndTimes(vTelescopeTriggerBits,fTelTriggerTimes);
	       isArrayTriggered = arraytrigger->RunTrigger(); //think of changing this if no array trigger is used
	       arrayTriggerBit = 0;
	       DeltaTL3=1e6;

		   //Run FADC
	       if(isArrayTriggered)
		 {
		   NumTriggeredEvents++;
		   arrayTriggerBit = 1;
		   DeltaTL3 = arraytrigger->GetL3DeltaT();
		   
		   if(DEBUG_MAIN)
		     cout<<"Array triggered"<<endl<<endl;
		   
		   //Readout Telescopes
		   for(UInt_t l =0 ; l<uNumTelescopes; l++)
		     {
		       
		       if(DEBUG_MAIN)
			 {
			   cout<<"Telescope "<<l<<endl;
			   cout<<"Average arrival time of photons: "<<telData[l]->GetAverageArrivalTime()<<endl;
			   cout<<"Telescope trigger time: "<<fTelTriggerTimes[l]<<endl;
			   cout<<"Array Trigger Time for this telescope: "<<arraytrigger->GetArrayTriggerTime(l)<<endl;
			 }

			  //This should only happen if something is wrong in real life with the array trigger and this telescope gets dropped
                if(  arraytrigger->GetArrayTriggerTime(l)>5e5 )
                     telData[l]->bTelescopeHasTriggered = kFALSE;
                else
                     telData[l]->bTelescopeHasTriggered = kTRUE;
					  

               //Do this
               Int_t telType = telData[l]->GetTelescopeType();
		       fadc[telType]->SetDebugInfo(fPrimaryEnergy,l,0,0);
		       //needs to change to work with option no Array trigger
		       telData[l]->fTriggerTime = arraytrigger->GetArrayTriggerTime(l) ;
		       fadc[telType]->RunFADC(telData[l]);
		    

			   if(DEBUG_MAIN)
				 cout<<"FADC completed"<<endl;
		     } 
		 }
	     }
         
		 if(DEBUG_MAIN)
			cout<<"writing event into file"<<endl; 
	
	     tSimulatedEvents.Fill(); 

	     //Write the event into the VBF File
	     if(readConfig->GetVBFwriteBit())
	       {
		 // check numbers to avoid atan2(0.,-0.) = 3.14159
		 if( fabs(fYcos) < 1.e-10 ) fYcos = 0.;

		 Float_t Az = atan2( fXcos, fYcos );
		 Float_t Ze = 1. - (fXcos*fXcos+fYcos*fYcos);
		 if( Ze < 0. ) Ze = 0.;  // should never happen
		 Ze = acos( sqrt( Ze ) );
		 // set telescope azimuth and elevation vectors in vbf writer
		 for (UInt_t tel1=0;tel1 < uNumTelescopes;tel1++){     
		   VBFwrite->setAzimElevTelDeg(tel1,Az*57.29577,90.0-Ze*57.29577);
		 }              
    
		 VBFwrite->setDataEvent();  // tell vw that the event is a data event and
		 //      not a pedestal event
		 
		 VBFwrite->makePacket();   // create an empty packet
	       
		 double rze;
		 double raz;
		 getRotatedShowerDirection( Ze*57.29577,Az*57.29577,fYsource,fXsource,rze,raz );
	       
		 // make simulation data header and add to packet
		 VBFwrite->makeVSimulationData( fPrimaryEnergy * 1000.0,
					      rze, //primary zenith angle fix
					      raz, //primary az angle
					      fXcore, //x:primary core east
					      -fYcore,//-y : primary core south
					      dObsHeight, //observatory altitude
					      iPrimaryType);

		 // if no array trigger, that is all we need to write to the VBF file, 
		 // needs to change if there is no array trigger used
		 if(isArrayTriggered==false) {
	     
		   // the storePacket method stores the packet and increments the eventtime.
		   VBFwrite->storePacket();           // do a writePacket
	     
		   continue;
		 
		 }//Finished writing a non-triggered event to the VBF file
	       
		 //now the case the event has triggered
	    
		 //we read out all telescope in the event of a trigger
		 VBFwrite->setTriggeredReadout(true);
	                                            
		 VBFwrite->makeArrayEvent();  // make an arrayevent, will add to it and store later

		 // setup trigger vector of triggered telescopes
		 std::vector<int> arrayTriggerV(uNumTelescopes,0); 
		 for(UInt_t tel1=0;tel1< uNumTelescopes;tel1++){
		   if ( telData[tel1]->bTelescopeHasTriggered ) {
		     arrayTriggerV[tel1] = 1;
		   }
    
		 }	    
		 // make an array trigger
		 VBFwrite->makeArrayTrigger(arrayTriggerV);
	   
		 /*Create the tracking information record for each telecope*/
		 /*We then write FADC output, trigger bit and hi lo gain switch bit*/
		 for (UInt_t tel=0;tel< uNumTelescopes;tel++){
	    
        
		   VBFwrite->setTelescope(tel);  // set current telescope in VBF writer
		   //vw->setDebugLevel(1);
            // these two methods will set the numFadcsamples and the MaxNumberChannels which
             //    I also set to the MaxNumber
            VBFwrite->setMaxNumberChannels(telData[tel]->iNumPixels); // set channelNum for this telescope
            VBFwrite->setNumFadcSamples( telData[tel]->iNumFADCSamples);  // set numFadcSamples for this telescope

		   VBFwrite->makeEvent();        // make the individual telescope event
		   //vw->setDebugLevel(0);
		 
		   for(UInt_t pix=0;pix<readConfig->GetNumberPixels( telData[tel]->GetTelescopeType() );pix++){
		     VBFwrite->setPixel(pix);         // set the pixel number in the VBF writer
		     VBFwrite->setChargePedHigain();  // initialize settings
		     
		     
		     //Set the Hi Lo gain bin true = low gain; false = high gain
		     VBFwrite->setHiLoGain( telData[tel]->GetPixelLowGainSwitch(pix) );
		     
		     VBFwrite->setTriggerBit( telData[tel]->GetGroupTrigger(pix) );
	       
		     /*Write the FADC trace out*/
		     vector<Int_t> trace = telData[tel]->GetFADCTrace(pix);
		     for(Int_t i=0;i<(Int_t)trace.size();i++){
		       VBFwrite->storeSample(i,trace[i]);  // store the sample
		     }	  
		  
		   } // pixel loop
		   if(DEBUG_MAIN)
			   cout<<"writing event into vbf format"<<endl;
		   VBFwrite->storeEvent();  // store the event for this telescope
		 } // end loop over telescopes all fadc traces and trigger bits passed to the vbfwriter
		 
       
		 VBFwrite->storeArrayTrigger();     // store in arrayevent
		 VBFwrite->storeArrayEvent();       // store in packet	    
		 VBFwrite->storePacket();           // do a writePacket, increment event time.

	       }                                                   
	     //Finished writing the event into the VBF file
	 
      if(DEBUG_MAIN)
		    cout<<"Event finished going into next event"<<endl;
	 }//Loop over to the next event



       
       fOut->cd("Events");
       tSimulatedEvents.Write();
       //change back to the root root directory
       gROOT->cd();

       cout<<"Have "<<NumTriggeredEvents<<" triggered events!"<<endl;
       cout<<"Have "<<NumSkippeddEvents<<" events that are skipped because no telescope had the min required number of Cherenkov photons in the focal plane"<<endl;
       fO->Close();



       //Finish up and close the vbf file
       if(readConfig->GetVBFwriteBit())
	 {
	   cout<<" number of vbf packets in data file: "<<VBFwrite->getEventNumber()<<endl;
	   cout<<" number of vbf data packets: "<<VBFwrite->getEventNumber() -  VBFwrite->getPedEventCount() - 1<<endl;
	   cout<<" number of vbf pedestal packets: "<<VBFwrite->getPedEventCount()<<endl;

	   VBFwrite->finishVBF();
	   delete  VBFwrite;
	 }
       
     }// Done looping over the events
   
   cout<<"Close output file"<<endl;

   fOut->Close();


   cout<<"Done everything... Have a pleasant day"<<endl;

   //Delete objects
/*
   delete fDeg;
   delete readConfig;
   delete arraytrigger;
   for(UInt_t t = 0 ; t<uNumTelescopes; t++)
   {
	   delete traceGenerator[t];
	   delete Teltrigger[t];
	   delete fadc[t];
   }
  */
}
