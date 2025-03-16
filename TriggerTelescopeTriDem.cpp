/* \file TriggerTelescopeTriDem.cpp
   Implementation of a next neighbor trigger for TriDem
   1. Signal of X pixels are summed
   2. Discrimination of summed signal
   3. next neighbor logic (multiplicity of Y) 
   4. Implementaton of the rate feedback in the discriminator path.
*/

#include "TriggerTelescopeTriDem.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>

#include <TMath.h>
#include <TTimer.h>
#include <Getline.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLine.h>

using namespace std;

//---------------------------------------------------------------------------------------
//Constructor
TriggerTelescopeTriDem::TriggerTelescopeTriDem(ReadConfig *readConfig, int telType, TRandom3 *generator,Bool_t debug,Display *display) : TriggerTelescopeBase(readConfig,telType,generator,debug,display)
{
  cout<<"Initializing the Telescope Trigger for TriDem"<<endl;

  if(readConfig->GetUseSumTrigger(telType))
   {
     cout<<"Disable the sumtrigger in the configuration file!!"<<endl;
     exit(-1);
   }

  iNumPixels = readConfig->GetNumberPixels(telType);
 
  cout<<"Setting up the groups"<<endl;
  iNumGroups = readConfig->GetNumberGroups(telType);
  iGroupNeighborIDs = readConfig->GetNeighborsOfGroup(telType);
  iGroupMembers = readConfig->GetMembersOfGroups(telType);
  if(iNumGroups <=0)
      {
        cout<<"You did not seem to have defined a proper number of groups: "<<iNumGroups<<endl;
        cout<<"do something about it!!"<<endl;
        exit(1);
      }
 cout<<"The camera is divided into "<<iNumGroups<<" groups."<<endl;


  //Creating the vectors that hold the IDs of triggered pixels and their times
  vector< int > i_pix;
  vTriggeredPixelInGroup = new vector< vector<int> >;
  vTriggeredPixelInGroup->assign(iNumGroups , i_pix );

  vector< float > f_pix;
  vTriggerTimesInGroup = new vector< vector<float> >;
  vTriggerTimesInGroup->assign(iNumGroups , f_pix );

 
}

TriggerTelescopeTriDem::~TriggerTelescopeTriDem()
{
  delete vTriggeredPixelInGroup;
  delete vTriggerTimesInGroup;
}

//---------------------------------------------------------------------------------------------
// Identifies the pixels in each group/MUSIC that triggered and saves their IDs and
// trigger times
void TriggerTelescopeTriDem::FindTriggeredGroups()
{
 
  //Loop over all groups
  for(int g=0;g<iNumGroups;g++)
    {
       vTriggeredPixelInGroup->at(g).clear();
       vTriggerTimesInGroup->at(g).clear();
       //cout<<"MUSIC "<<g<<endl;
         //Loop over all pixels in a group
         for(unsigned p=0;p<iGroupMembers[g].size();p++)
           {
             //cout<<"pixel "<<iGroupMembers[g][p]<<endl;
             if(telData->bTriggeredTriggerPixels[iGroupMembers[g][p]])
               {
                 //save the pixel ID and trigger time if the pixel is triggered
                // cout<<g<<"  "<<iGroupMembers[g][p]<<"  "<<telData->fDiscriminatorTime[iGroupMembers[g][p]]<<endl;
                 vTriggeredPixelInGroup->at(g).push_back(iGroupMembers[g][p]);
                 vTriggerTimesInGroup->at(g).push_back(telData->fDiscriminatorTime[iGroupMembers[g][p]]);
                 //cout<<iGroupMembers[g][p]<<"  "<<telData->fDiscriminatorTime[iGroupMembers[g][p]]<<endl;
               }
           }
    }
}


//----------------------------------------------------------------------------------------------
// Identify the trigger time if a trigger accured
//
bool TriggerTelescopeTriDem::FindTrigger()
{

  telData->bTelescopeHasTriggered = false;
 
 Float_t *vTriggerTimes = new Float_t[iNumGroups]; //How many groups are in each cluster of this patch
  for(int i=0;i<iNumGroups;i++)
   vTriggerTimes[i]=1e10;

  iPixelTriggeredInPatch.resize(1);
 
  //Loop over all groups and identify if a pixel in a neighboring group fired in coincidence
  for(int g=0;g<iNumGroups;g++)
    {
       //cout<<"Group "<<g<<" has "<<vTriggeredPixelInGroup->at(g).size()<<" triggered pixel"<<endl;
       //loop over all pixels in this group that triggered and find the pixel
       //that triggered first.
       for(unsigned p=0;p<vTriggeredPixelInGroup->at(g).size();p++)
          {
             //get triggered time
             vTriggerTimes[g]  = vTriggerTimesInGroup->at(g)[p]< vTriggerTimes[g] ? vTriggerTimesInGroup->at(g)[p] : vTriggerTimes[g] ;
             //we have a trigger, save the trigger time and
             telData->bTelescopeHasTriggered = true;
          }

    } 


   //Find the first trigger time
   telData->fTelescopeTriggerTime = 1e10;
   if(telData->bTelescopeHasTriggered)
     {
       //Sort all trigger times descending
       Int_t* index = new Int_t[iNumGroups];
       TMath::Sort(iNumGroups,vTriggerTimes,index,false);
        //cout<<vCoincidenceTriggerTimes[index[t]]<<"  "<<vTrigPixsInCluster->at(index[t])[0]<<"  "<<vTrigPixsInCluster->at(index[t])[1]<<endl;
        telData->fTelescopeTriggerTime = vTriggerTimes[index[0]];
        telData->vTriggerCluster.clear();
        telData->vTriggerCluster.push_back(index[0]);
        //cout<<"Telescope trigger at "<<telData->fTelescopeTriggerTime<<" trigger bit "<<telData->bTelescopeHasTriggered<<endl;
      }
   delete vTrigPixsInCluster;
   delete [] vTriggerTimes;

   return telData->bTelescopeHasTriggered;
}

//This does the actual triggering
//It returns true if the telescope has triggered and false if not
Bool_t  TriggerTelescopeTriDem::RunTrigger()
{
  // cout<<"TriggerTelescopeTriDem::RunTrigger"<<endl;


   //Run the CFD simulation on all pixels 
   telData->iNumTriggeredTriggerPixels=0;
   for(Int_t i=0;i<iNumTriggerPixels;i++)
    {
      if(RunDiscriminator(i))
	     telData->iNumTriggeredTriggerPixels++;
    }

   //cout<<"Done with the CFD"<<endl<<endl;
   //cout<<telData->iNumTriggeredTriggerPixels<<" Trigger pixel triggered"<<endl;

  //Implementation of the TriDem Trigger logic
  //First figure out which MUSIC chips have triggered
  FindTriggeredGroups();

  telData->bTelescopeHasTriggered = FindTrigger();

  

  if(telData->bTelescopeHasTriggered && bDebug)
     {
       cout<<"Event triggered telescope"<<endl; 
       cout<<"at time "<< telData->fTelescopeTriggerTime<<endl;
     }

  return telData->bTelescopeHasTriggered;
}



