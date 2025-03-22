#include <iomanip>

//Macro that writes the Trinity Demonstrator camera configuration for CARE
//Writes a configuration for rectangular pixel

int iYrows = 4; //siab rows in the camera
int iSIABID[] = {17, 16, 1, 0, 25, 24, 9, 8, 19, 18, 3, 2, 27, 26, 11, 10, 21, 20, 5, 4, 29, 28, 13, 12, 23, 22, 7, 6, 31, 30, 15, 14};  //number of siabs in each row
//from bottom to top of camera the NSB in each row
//float fNSBPerRow[] = {70, 52, 38, 28, 21, 16, 12, 8.6, 6.4, 4.7, 3.5, 2.6, 1.9, 1.4, 1.0, 0.8};
float fNSBPerRow[] = {21, 21, 19, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
int iModulesInRow[] = {4, 4, 4, 4};  //number of siabs in each row
float fModulePitch = 2.5; //cm //pitch of siabs
float fPixelPitch = 0.625;   //[cm] pitch of pixels within module
int iNPixelsInRow = 4;   //number of pixels in one row in a module (number of pixel is this number squared

int iNSummedPixel = 8; //number of pixel that get summed to one SIAB

int iTelID = 0;


ofstream *fConfigFile = new ofstream("TrinityDemonstrator_Camera.cfg");     //! ascii file


void WriteCameraGeometryTrinityDemonstrator()
{
 cout<<"this macro writes out the Geometry of a camera in a format that can be put into CARE"<<endl;


   //Trigger patch list records (PATCH)
   //   -Telescope ID
   //   -Patch identification number
   //   -A flag indicating if the patch is enabled (1) or not (0).
   //   -List of id numbers for groups in patch (0 indicates no pixel connected)
   //    The order of the 19 pixel list is important.
   //    Note that this is not pixels but summed pixels, i.e. groups
   //
   //   * PATCH 1 1 1  12  13  14  15  16 11    4   5   6  17  10   3   1   7  18   9   2  19   8
   

   gSystem->Load("libGeom");
   new TGeoManager("world", "the simplest geometry");

   TGeoMaterial *matAl = new TGeoMaterial("Al",0,0,0);
   TGeoMedium *medAl = new TGeoMedium("Al",1,matAl);
   TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
   TGeoMedium *med = new TGeoMedium("Vacuum",1,mat);
   TGeoVolume *top = gGeoManager->MakeBox("Top",med,37.,37.,50.);
   top->SetLineColor(kMagenta);
   TGeoVolume *PMTModule = gGeoManager->MakeBox("SIAB",med,fModulePitch/2.,fModulePitch/2,10.);
   PMTModule->SetLineColor(kRed);
   PMTModule->SetLineWidth(2);
   TGeoVolume *PMTPixel = gGeoManager->MakeBox("PMTPixel",medAl,fPixelPitch/2,fPixelPitch/2,1.);
   PMTPixel->SetLineColor(kGreen);
   top->Raytrace();
   PMTPixel->Raytrace();



   //Create SIAB         
   //Coordinates of center of pixels with respect to center of SIAB in mm
   vector<float> vXCoordPixelInModule;
   vector<float> vYCoordPixelInModule;
   
   //vector that stores the sumgroup ID of the pixel within the SIAB 
   vector<int> vSumGroupIDinModule;

   for(int x = 0; x<iNPixelsInRow; x++)
   {
	 cout<<x<<endl;
	 float xPixelCoord = (x-iNPixelsInRow/2)*fPixelPitch + fPixelPitch/2;  
     for(int y = 0; y<iNPixelsInRow; y++)
	 {
      float yPixelCoord = (y-iNPixelsInRow/2)*fPixelPitch + fPixelPitch/2;
      cout<<"x: "<<xPixelCoord<<" y: "<<yPixelCoord<<endl;
TGeoTranslation *tr = new TGeoTranslation(xPixelCoord,yPixelCoord, 0.);
	  PMTModule->AddNode(PMTPixel,y*10+x,tr);
      vXCoordPixelInModule.push_back(xPixelCoord*10);
      vYCoordPixelInModule.push_back(yPixelCoord*10);

      int sumgroup = x/2;
	  vSumGroupIDinModule.push_back(sumgroup);
	  cout<<sumgroup<<endl;
	 }
    }
    TGeoTranslation *tr = new TGeoTranslation(-35,-35, 0.);
    top->AddNode(PMTModule,100000,tr);

   
   //Create modules
   cout<<endl<<"Creating Modules"<<endl;
   vector<float> vXCoordModule;
   vector<float> vYCoordModule;
   for(int y = 0; y<iYrows; y++)
   {
	  cout<<y<<endl;
	  float yModuleCoord = (y-iYrows/2)*fModulePitch+fModulePitch/2.;
	  for(int x = 0; x<iModulesInRow[y];x++)
	  {
		 float xModuleCoord = (x-iModulesInRow[y]/2)*fModulePitch+fModulePitch/2.; 
         cout<<"y: "<<yModuleCoord<<" x: "<<xModuleCoord<<endl;
         vXCoordModule.push_back(xModuleCoord*10);
         vYCoordModule.push_back(yModuleCoord*10);
         TGeoTranslation *tr = new TGeoTranslation(xModuleCoord,yModuleCoord, 0.);
         top->AddNode(PMTModule,y*15+x,tr);
	  }
   }

   gGeoManager->SetTopVolume(top);
   gGeoManager->CloseGeometry();
   gGeoManager->SetVisLevel(2);

///   gGeoManager->SetTopVisible();
   top->Draw();
   //Write the SiPM pixel parameters
   UInt_t uPixID = 0;
   for(UInt_t m=0; m<vYCoordModule.size();m++)
   {
     for(UInt_t p=0; p<vYCoordPixelInModule.size();p++)
	 {
		 *fConfigFile <<"* SIPMPARS "<<iTelID<<"  "<<uPixID<<" 14000 0.02"<<endl;
		uPixID++;
	 }

   }

   *fConfigFile<<endl<<endl;


   //Write the list with the pixel coordinates
   //   * PMPIX 1  1  0.000000  0.000000 1 14.8 1 1 1 
   uPixID = 0;
   vector<float> vXCoordPixel;
   vector<float> vYCoordPixel;

   UInt_t uSumGroupID = 0;
   vector<int> vSumGroupIDs; //The sumgroupID of each pixel

   for(UInt_t m=0; m<vYCoordModule.size();m++)
   {
     for(UInt_t p=0; p<vYCoordPixelInModule.size();p++)
	 {
		float xc = vXCoordModule[m] + vXCoordPixelInModule[p];
		float yc = vYCoordModule[m] + vYCoordPixelInModule[p];
                int row = 7+int((yc+5*fPixelPitch)/(10*fPixelPitch));
                cout<<"row "<<row<<" NSB "<<fNSBPerRow[row]<<endl;
                *fConfigFile <<"* PMPIX "<<"0 "<<iTelID<<" 1  "<<uPixID<<" "<<setw(5)<< setprecision(5)<<fNSBPerRow[row]<<" "<<setw(10)<< setprecision(5)<<xc<<"   "<<yc<<" "<<fPixelPitch/sqrt(2.0)/2*10<<" 0.0"<<endl;
		vXCoordPixel.push_back(xc);
		vYCoordPixel.push_back(yc);

		//The sum groups
		int sumgroupid = vSumGroupIDinModule[p]+m*2;
               // cout<<"pix: "<<uPixID<<" GroupIDinModule: "<<vSumGroupIDinModule[p]<<" m: "<<m<<" sumgroupid "<<sumgroupid<<"  "<<uPixID/16<<endl;
	    vSumGroupIDs.push_back(iSIABID[sumgroupid]);
	  
		uPixID++;
        
	 }

   }

   *fConfigFile<<endl<<endl;


   //Write neighbors of each pixel
   //Neighbors list records (NGHBR)
   //   -Telescope identification number
   //   -Pixel identification number
   //   -Number of neighbors
   //   -List of neighbors identification numbers
   //
   //    * NGHBR 1     1     6      2     3     4     5     6     7
   //               
   for(UInt_t p=0;p<vXCoordPixel.size();p++)
   {
	   float xp = vXCoordPixel[p];
	   float yp = vYCoordPixel[p];
	   vector<int> nghbrIDs;
       for(UInt_t np=0;np<vXCoordPixel.size();np++)
	   {
		   if(fabs(vXCoordPixel[np]-xp)<=1.5*fPixelPitch*10 && fabs(vYCoordPixel[np]-yp)<=1.5*fPixelPitch*10 && np!=p )
			   nghbrIDs.push_back(np);
	   }

	   *fConfigFile<<"* PIXNGHBR "<<iTelID<<"  "<<p<<"  "<<nghbrIDs.size();
	   
       for(UInt_t np=0;np<nghbrIDs.size();np++)
              *fConfigFile<<"  "<<nghbrIDs[np];
	   *fConfigFile<<endl;
   }
   *fConfigFile<<endl<<endl<<endl;



   //Summed pixels groups
   // - Telescope ID
   // -Group identification number
   // -number of pixel in group
   // -Pixel ids
   // * GROUP 1 1 2 3 4 ...
   vector<float> vXCoordGroup;
   vector<float> vYCoordGroup;

   UInt_t numbergroups = vYCoordPixelInModule.size()/iNSummedPixel*vYCoordModule.size();
   for(UInt_t g = 0; g<numbergroups; g++)
   {
	   vector<int> pixelingroup;
	   float xc = 0.0;
	   float yc = 0.0;
	   for(UInt_t p=0;p<vSumGroupIDs.size();p++)
	   {
		   if(vSumGroupIDs[p]==g)
	        {
				pixelingroup.push_back(p);
				xc+=vXCoordPixel[p];
				yc+=vYCoordPixel[p];
            }
	   }

       vXCoordGroup.push_back(xc/iNSummedPixel);
       vYCoordGroup.push_back(yc/iNSummedPixel);
	   *fConfigFile<<"* GROUP "<<iTelID<<" "<<pixelingroup.size()<<" "<<g<<" ";
	   for(UInt_t p=0;p<pixelingroup.size();p++)
		       *fConfigFile<<pixelingroup[p]<<"  ";
	   *fConfigFile<<endl;
   }
	   *fConfigFile<<endl<<endl<<endl;



   //Write neighbors of groups
   // - Telesocpe ID
   // - Group ID
   // - Number of neighbors
   // - IDs of neighboring groups
   //   * GRPNGHBR 1 1 3 4 5 7 
   for(UInt_t g=0;g<vXCoordGroup.size();g++)
   {
	   float xg = vXCoordGroup[g];
	   float yg = vYCoordGroup[g];
cout<<xg<<" "<<yg<<endl;
	   vector<int> nghbrIDs;
       for(UInt_t ng=0;ng<vXCoordGroup.size();ng++)
	   {
		   if(fabs(vXCoordGroup[ng]-xg)<=(iNPixelsInRow/2)*fPixelPitch*10 && fabs(vYCoordGroup[ng]-yg)<=(iNPixelsInRow)*fPixelPitch*10 && ng!=g )
                    {
			   nghbrIDs.push_back(ng);
                           cout<<g<<"  "<<ng<<"  "<<vXCoordGroup[ng]<<"  "<<vYCoordGroup[ng]<<endl;
                     }
	   }
	   *fConfigFile<<"* GRPNGHBR "<<iTelID<<"  "<<g<<"  "<<nghbrIDs.size();
       cout<<g<<endl;	   
       for(UInt_t ng=0;ng<nghbrIDs.size();ng++)
       {
cout<<ng<<" "<<iSIABID[ng]<<"; ";
              *fConfigFile<<"  "<<nghbrIDs[ng];
       }
cout<<(iNPixelsInRow+1)*fPixelPitch*10<<endl;
	   *fConfigFile<<endl;
       cout<<endl;
   }
   *fConfigFile<<endl<<endl<<endl;


  
   fConfigFile->close();
}
