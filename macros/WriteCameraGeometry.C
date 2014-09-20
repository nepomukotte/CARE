#include <iomanip>

//Macro that writes a camera configuration for CARE
//Writes a configuration for rectangular pixel

int iYrows = 15; //module rows in the camera
int iModulesInRow[] = {5, 9, 11, 13, 13, 15, 15, 15, 15, 15, 13, 13, 11, 9, 5};  //number of modules in each row
float fModulePitch = 5.4; //cm //pitch of modules
float fPixelPitch = 0.6506;   //pitch of pixels within module
int iNPixelsInRow = 8;   //number of pixels in one row in a module (number of pixel is this number squared

int iNSummedPixel = 4; //number of pixel that get summed to one group

int iTelID = 0;


ofstream *fConfigFile = new ofstream("SC_Camera.cfg");     //! ascii file


void WriteCameraGeometry()
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
   TGeoVolume *PMTModule = gGeoManager->MakeBox("PMTModule",med,fModulePitch/2.,fModulePitch/2,10.);
   PMTModule->SetLineColor(kRed);
   PMTModule->SetLineWidth(2);
   TGeoVolume *PMTPixel = gGeoManager->MakeBox("PMTPixel",medAl,fPixelPitch/2,fPixelPitch/2,1.);
   PMTPixel->SetLineColor(kGreen);
   top->Raytrace();
   PMTPixel->Raytrace();



   //Create Module             
   //Coordinates of center of pixels with respect to center of module in mm
   vector<float> vXCoordPixelInModule;
   vector<float> vYCoordPixelInModule;
   
   //vector that stores the sumgroup ID of the pixel within the module
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

      int sumgroup = y/((int)sqrt(iNSummedPixel))+x/((int)sqrt(iNSummedPixel))*iNSummedPixel;
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
	  float yModuleCoord = (y-iYrows/2)*fModulePitch;
	  for(int x = 0; x<iModulesInRow[y];x++)
	  {
		 float xModuleCoord = (x-iModulesInRow[y]/2)*fModulePitch; 
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
		 *fConfigFile <<"* SIPMPARS "<<iTelID<<"  "<<uPixID<<" 3464 0.4"<<endl;
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
		 *fConfigFile <<"* PMPIX "<<iTelID<<" 1  "<<uPixID<<" "<<setw(10)<< setprecision(5)<<xc<<"   "<<yc<<" "<<fPixelPitch*sqrt(2.0)*10<<" 0.0"<<endl;
		vXCoordPixel.push_back(xc);
		vYCoordPixel.push_back(yc);
		uPixID++;

		//The sum groups
		int sumgroupid = vSumGroupIDinModule[p]+vYCoordPixelInModule.size()/iNSummedPixel*m;
	    vSumGroupIDs.push_back(sumgroupid);
	  
        
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
	   vector<int> nghbrIDs;
       for(UInt_t ng=0;ng<vXCoordGroup.size();ng++)
	   {
		   if(fabs(vXCoordGroup[ng]-xg)<=(sqrt(iNSummedPixel)+1)*fPixelPitch*10 && fabs(vYCoordGroup[ng]-yg)<=(sqrt(iNSummedPixel)+1)*fPixelPitch*10 && ng!=g )
			   nghbrIDs.push_back(ng);
	   }
	   *fConfigFile<<"* GRPNGHBR "<<iTelID<<"  "<<g<<"  "<<nghbrIDs.size();
	   
       for(UInt_t ng=0;ng<nghbrIDs.size();ng++)
              *fConfigFile<<"  "<<nghbrIDs[ng];
	   *fConfigFile<<endl;
   }
   *fConfigFile<<endl<<endl<<endl;


  
   fConfigFile->close();
}
