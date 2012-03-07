/* \file TriggerTelescopeNextNeighbor.cpp
   Implementation of a next neighbor trigger for AGIS
   1. Signal of X pixels are summed
   2. Discrimination of summed signal
   3. next neighbor logic (multiplicity of Y) 
*/

#include "ArrayTrigger.h"
#include <iostream>
#include <math.h>
#include <TMath.h>
#include <sstream>
#include <fstream>

using namespace std;


//Constructor
ArrayTrigger::ArrayTrigger(ReadConfig *readConfig,Bool_t debug)
{
  cout<<"Welcome to the Array Trigger"<<endl;
  bDebug = debug;
  iNumTel = -1;
  iMultiplicity = -1;
  bNextNeighborReq = false;
  iClusterID = NULL;
  fDeltaL3Cluster = NULL;        
  fCoincidence = -1;
  fDeltaTL3 = -10000;

  SetParametersFromConfigFile( readConfig );
  vTelTransitTimes.assign(iNumTel,0.0);
}

//Destructor
ArrayTrigger::~ArrayTrigger()
{
}

//-----------------------------------------------------------------------------------
// Set the inter telescope transit times
void ArrayTrigger::SetInterTelTransitTimes(vector<float> vTimes)
{
  if(vTimes.size() != (unsigned)iNumTel)
   {
     cout<<"The number of telescopes does not match the number of transit times we have in the vector"<<endl;
     exit(1);
   }
   vTelTransitTimes = vTimes;
}


//-----------------------------------------------------------------------------------
// Set the Trigger Bit and Time vectors which store the information whether a telescope
// has triggered and what the time of trigger is
void ArrayTrigger::SetTelescopeTriggerBitsAndTimes(vector<Bool_t> TelTriggerBit, vector<Float_t> TelTriggerTimes)
{ 

  if( TelTriggerBit.size() != (UInt_t)iNumTel || TelTriggerTimes.size() != (UInt_t)iNumTel )
    {
      cout<<"Vector sizes with telescope trigger bits and times  do not match with the number of telescopes in the array!!"<<endl;
      cout<<"Telescopes: "<<iNumTel<<" TriggerBits: "<<TelTriggerBit.size()<<" Times: "<<TelTriggerTimes.size()<<endl;
      exit(1);
    }

  if( TelTriggerBit.empty() || TelTriggerTimes.empty() )
    {
      cout<<"Vector with telescope trigger bits or times  are empty!!"<<endl;
      cout<<"TriggerBits: "<<TelTriggerBit.size()<<" Times: "<<TelTriggerTimes.size()<<endl;
      exit(1);
    }

  bTelTriggerBit = TelTriggerBit; 
  fTelTriggerTimes = TelTriggerTimes;

  //Subtracting the inter telescope transit time from the Trigger times
  for(unsigned i=0;i<fTelTriggerTimes.size();i++)
   {
     fTelTriggerTimes[i] = fTelTriggerTimes[i] - vTelTransitTimes[i];
   }

}


//----------------------------------------------------------------------------
// After the trigger bits for each telescope have been set with SetTelescopeTriggerBits,
// this function is called. It takes care that the trigger is evaluated either with
// the next neighbor requirement or not. It returns a boolean. 1 The array has triggered
// 0 the array has not triggered

bool ArrayTrigger::RunTrigger()
{

  if(bDebug)
	  cout<<"In the Array trigger"<<endl;

  if(iNumTel<=0)
    {
      cout<<"Have no telescopes to work with. Do SetTelescopeCoordinates(vector<float> x, vector<float> y) first"<<endl;
      exit(1);
    }

  if(fCoincidence<0)
    {
      cout<<"array coincidence window is <0, BAD!!"<<endl;
      exit(1);
    }

  bool ArrayHasTriggered = false;

  Int_t iRefTelescope = 0;
  while(!bTelTriggerBit[iRefTelescope] &&  iRefTelescope<iNumTel-1)
    iRefTelescope++;


  //No telescope triggered, I guess we are done
  if(!bTelTriggerBit[iRefTelescope])
    return  false;

  if(bDebug)
  cout<<"Correcting or not correcting for the propagation delay"<<endl;

  if(bDebug)
  cout<<"Going into the array trigger"<<endl;  


  fDeltaTL3 = 1e6;


  vTelTriggerTimesAfterArrayTrigger.assign(iNumTel,1e6);
  //The case if we do not have an array trigger
  if(iMultiplicity==1)
  {
      for( Int_t i = 0; i<iNumTel; i++)        
	   {
	      if(bTelTriggerBit[iRefTelescope])
		  {
			if(bDebug)
				cout<<"Telescope "<<i<<" triggered at time "<<fTelTriggerTimes[i]<<endl;
	        vTelTriggerTimesAfterArrayTrigger[i]= fTelTriggerTimes[i];
			ArrayHasTriggered = true;
		  }
       }

  }
  //Loop over Telescopes and find clusters if we require next neighbors
  else if(bNextNeighborReq)
    {
      ArrayHasTriggered = RunTriggerWithNextNeighborRequirement();
    }
  else //run array trigger with multiplicity requirement but could be any combination of telescopes in the array
    {   

      if(bDebug)
		  cout<<"No next neighbor coincidence and multipl > 1"<<endl;
		  
      Float_t fArrayTriggerTime = 1e6;
      
      for(Int_t i = 0; i < iNumTel; i++)
	{
	  if(bTelTriggerBit[i])
	    {

	      Float_t *vTelTriggerTimes = new Float_t[iNumTel];
	      vTelTriggerTimes[0] = fTelTriggerTimes[i];
              
	      Int_t NumTelFullfilTrigger = 1;
	      for(Int_t a = 0; a < iNumTel; a++)
		{
		  Float_t fDeltaT = fabs(fTelTriggerTimes[i]-fTelTriggerTimes[a]);
		  if( bTelTriggerBit[a] && i!=a && fDeltaT<fCoincidence)
		    {
		      fDeltaTL3 = fDeltaT < fDeltaTL3 ?  fDeltaT : fDeltaTL3;
		      if(bDebug)
			  {
			  cout<<fDeltaT<<endl;
		      cout<<i<<" i "<<fTelTriggerTimes[i]<<endl;
		      cout<<a<<" a "<<fTelTriggerTimes[a]<<endl;
               }
		      vTelTriggerTimes[NumTelFullfilTrigger]=fTelTriggerTimes[a];

		      NumTelFullfilTrigger++;
		    }
		}
	      if(NumTelFullfilTrigger>=iMultiplicity)
		{  
		  if(bDebug)
		     cout<<NumTelFullfilTrigger<<"  "<<iMultiplicity<<endl;
		  Int_t* index = new Int_t[NumTelFullfilTrigger];
		  TMath::Sort(NumTelFullfilTrigger,vTelTriggerTimes,index,kFALSE);
                  fArrayTriggerTime = vTelTriggerTimes[index[iMultiplicity-1]] < fArrayTriggerTime ? 
                                                                 vTelTriggerTimes[index[iMultiplicity-1]] : fArrayTriggerTime ;
		  ArrayHasTriggered = true;
          vTelTriggerTimesAfterArrayTrigger.assign(iNumTel,fArrayTriggerTime);
		  if(bDebug)
		    cout<<"Trigger time "<<fArrayTriggerTime<<endl;
		  delete [] index;
		}
	      delete [] vTelTriggerTimes;
	    }

	}
    }


  //adding back the inter telescope delay
  for(unsigned i=0;i<vTelTriggerTimesAfterArrayTrigger.size();i++)
   {
     if(bDebug)
       cout<<i<<" Adding back inter telescope transit time "<<vTelTransitTimes[i]<<endl;
     vTelTriggerTimesAfterArrayTrigger[i] += vTelTransitTimes[i];
   }

  

  if(bDebug)
  cout<<"Final answer "<<ArrayHasTriggered<<endl;
  
  return ArrayHasTriggered;

}


// ----------------------------------------------------------------------------
//
//Helper function to loop over all neighbours of a telescope to find telescopes belonging to the cluster
//returns the number of telescope in that cluster

Int_t ArrayTrigger::CalcCluster(Int_t TelID, Int_t ClusterID)
{
 
  // If a cluster number was already assigned to this telescope... do nothing.
  if (iClusterID[TelID]==ClusterID)
      return 0;

  // Assign the new cluster ID this telescope
  iClusterID[TelID]=ClusterID;

  // Need this to store the number of groups in the cluster
  Int_t NumTelInCluster = 1;
 
  // Now do the same with all its neighbors
  for(UInt_t n = 0; n<neighbors[TelID].size();n++)
    {
      Int_t GroupIDNeighbor = neighbors[TelID][n];
      Float_t DeltaT = fabs( fTelTriggerTimes[ClusterID]- fTelTriggerTimes[GroupIDNeighbor]);
      if (bTelTriggerBit[GroupIDNeighbor] &&  DeltaT<fCoincidence)
	{
	  fDeltaL3Cluster[ClusterID] = DeltaT< fDeltaL3Cluster[ClusterID] ? DeltaT : fDeltaL3Cluster[ClusterID];
          NumTelInCluster += CalcCluster(GroupIDNeighbor,ClusterID);

	}
    }


  // return the number of groups in this cluster
  return NumTelInCluster;
}


//----------------------------------------------------------------------
//This function is called if the trigger requirements include that the
//telescopes that can positively contribute to the trigger decision
//have to be next neighbors

bool ArrayTrigger::RunTriggerWithNextNeighborRequirement()
{

  bool bArrayHasTriggered = false;

  //loop over all telescopes. If that telescopes receives enough triggers
  //including from itself within the coincidence window then it will readout. 
  //the time for the readout is determined by the time the multiplicity requirement is reached
  for(Int_t i=0;i<iNumTel;i++)
    {
      if(bDebug)
		  cout<<"Array Trigger With next neighbor: Checking Telescope "<<i<<endl;

      //initiate the array for the triggertimes 
	  int iNumTriggeredTelescopes=0;
      Float_t *fTriggerTimesInGroup = new Float_t[iNumTel];

	  //add trigger time if telescope itself triggered
      if (bTelTriggerBit[i])
	  {
		  iNumTriggeredTelescopes++;
		  fTriggerTimesInGroup[0]=fTelTriggerTimes[i];
		  if(bDebug)
			  cout<<" telescope itself triggered at time "<<fTriggerTimesInGroup[0]<<endl;
	  }

      //loop over all neighbors of this telescope
      for(UInt_t n = 0; n<neighbors[i].size();n++)
        {
           Int_t GroupIDNeighbor = neighbors[i][n];
		   if(bDebug)
		      cout<<"checking telescope "<<GroupIDNeighbor<<endl; 
           if (bTelTriggerBit[GroupIDNeighbor])
		     {
		         fTriggerTimesInGroup[iNumTriggeredTelescopes]=fTelTriggerTimes[GroupIDNeighbor];
                 if(bDebug)
					 cout<<"Telescope triggered at time "<<fTriggerTimesInGroup[iNumTriggeredTelescopes]<<endl;
		         iNumTriggeredTelescopes++;
             }
	    }
	
      //Sort all times ascending
      Int_t* index = new Int_t[iNumTel];
      TMath::Sort(iNumTriggeredTelescopes,fTriggerTimesInGroup,index,false);

      //Now go through the times and find when we have reached the necessary multiplicity earliest
      int iFirstToTest = 0;
      while(iFirstToTest+iMultiplicity-1<iNumTriggeredTelescopes)
	  {
       Float_t fDeltaT = fTriggerTimesInGroup[index[iFirstToTest+iMultiplicity-1]] 
	                   - fTriggerTimesInGroup[ index[iFirstToTest] ];
       
       if(bDebug)
       cout<<"checking between "<<iFirstToTest<<" and "<<iFirstToTest+iMultiplicity-1<<" with times "<<fTriggerTimesInGroup[ index[iFirstToTest] ]<<" and "<<fTriggerTimesInGroup[index[iFirstToTest+iMultiplicity-1]]<<endl;

       if(fDeltaT<2*fCoincidence)
		   {
		     fDeltaTL3 = fDeltaT;

			 vTelTriggerTimesAfterArrayTrigger[i] = fTriggerTimesInGroup[index[iFirstToTest+iMultiplicity-1]];
             bArrayHasTriggered = true;
			 if(bDebug)
				 cout<<"Ok, array triggered and this telescope is readout at time "<<vTelTriggerTimesAfterArrayTrigger[i]<<endl;
             break;                  
		   }
	   iFirstToTest++;
	  }

      delete [] index;
      delete [] fTriggerTimesInGroup;
	 
	}
                          
    return bArrayHasTriggered;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Reads in  the config file and sets all variables
void  ArrayTrigger::SetParametersFromConfigFile( ReadConfig *readConfig ){

   cout <<endl<< "ArrayTrigger::SetParametersFromConfigFile " <<  endl;

   iMultiplicity = -1;
   bNextNeighborReq = false;
   fCoincidence = -1;
      
   //How many Telescopes need to trigger to get an array trigger
   iMultiplicity = readConfig->GetTelescopeMultiplicity();
   cout<<"The multiplicity requirement for an array trigger is "<<  iMultiplicity<<endl;

   //Defines wether the array trigger requires next neighbor requirement 1 or not 0
   bNextNeighborReq  = readConfig->GetArrayTriggerNextNeighborRequirement();
   cout<<"The array next neighbor requirement is "<<bNextNeighborReq<<endl;

   //The coincidence window on array level in ns
   fCoincidence = readConfig->GetArrayCoincidenceWindow();
   cout<<"The coincidence window in ns on array level is "<<fCoincidence<<endl;

   //Set the telescope neighbors
   neighbors = readConfig->GetTelescopeNeighbors();

   //Set the number of telescopes in the simulated array
   iNumTel = readConfig->GetNumberOfTelescopes();


}
