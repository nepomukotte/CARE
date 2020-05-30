/* \file TriggerTelescopeSPB2.cpp
   Implementation of a next neighbor trigger for SPB2
   1. Signal of X pixels are summed
   2. Discrimination of summed signal
   3. next neighbor logic (multiplicity of Y) 
   4. Implementaton of the rate feedback in the discriminator path.
*/

#include "TriggerTelescopeSPB2.h"
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
TriggerTelescopeSPB2::TriggerTelescopeSPB2(ReadConfig *readConfig, int telType, TRandom3 *generator,Bool_t debug,Display *display) : TriggerTelescopeBase(readConfig,telType,generator,debug,display)
{
  cout<<"Initializing the Telescope Trigger for SPB2"<<endl;

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

TriggerTelescopeSPB2::~TriggerTelescopeSPB2()
{
  delete vTriggeredPixelInGroup;
  delete vTriggerTimesInGroup;
}

//---------------------------------------------------------------------------------------------
// Identifies the pixels in each group/MUSIC that triggered and saves their IDs and
// trigger times
void TriggerTelescopeSPB2::FindTriggeredGroups()
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
// Execute the Coincidence logic. Identify the trigger time if a trigger accured
//
bool TriggerTelescopeSPB2::RunCoincidenceLogic()
{

  telData->bTelescopeHasTriggered = false;
 
  Float_t *vCoincidenceTriggerTimes = new Float_t[iNumGroups]; //How many groups are in each cluster of this patch
  for(int i=0;i<iNumGroups;i++)
   vCoincidenceTriggerTimes[i]=0;

  iPixelTriggeredInPatch.resize(1);
 
  vector< int > i_pix;
  vTrigPixsInCluster = new vector< vector<int> >;
  vTrigPixsInCluster->assign(iNumGroups , i_pix );

  //Loop over all groups and identify if a pixel in a neighboring group fired in coincidence
  for(int g=0;g<iNumGroups;g++)
    {
       //cout<<"Group "<<g<<" has "<<vTriggeredPixelInGroup->at(g).size()<<" triggered pixel"<<endl;
       //loop over all pixels in this group that triggered
       for(unsigned p=0;p<vTriggeredPixelInGroup->at(g).size();p++)
          {
             //get triggered time
             float fTriggerTime = vTriggerTimesInGroup->at(g)[p];
             //cout<<"Trigger time this pixel "<<fTriggerTime<<endl;
             //loop over the triggered pixels in the neighboring groups 
             for(unsigned ng=0;ng<iGroupNeighborIDs[g].size();ng++)
               {
                 int GroupNeighborID = iGroupNeighborIDs[g][ng];
                 //cout<<"Checking out neighbor group "<<GroupNeighborID<<endl;
                 if(GroupNeighborID>g) //only do this if we have not already visited the group
                   {
                     for(unsigned np=0;np<vTriggerTimesInGroup->at(GroupNeighborID).size();np++)
                        {
                          //cout<<vTriggerTimesInGroup->at(GroupNeighborID)[np]<<endl;
                          //for each of the triggered pixel in a neighbor group find out if it is in coincidence
                          //we only do this if we dont already have a trigger,
                          //this is not entirely correct we would have to check
                          //which pixels triggered first to catch the case where
                          //several pixels trigger at different times
                          if(fabs(fTriggerTime-vTriggerTimesInGroup->at(GroupNeighborID)[np])<fWidthDiscriminator && vTrigPixsInCluster->at(g).size()==0)
                            {
                               // cout<<"trigger"<<endl;
                               //we have a trigger, save the trigger time and
                               telData->bTelescopeHasTriggered = true;
                               vCoincidenceTriggerTimes[g] = fTriggerTime>vTriggerTimesInGroup->at(GroupNeighborID)[np] ? fTriggerTime : vTriggerTimesInGroup->at(GroupNeighborID)[np];
                               vTrigPixsInCluster->at(g).push_back(g);
                               vTrigPixsInCluster->at(g).push_back(GroupNeighborID);
                               //cout<<"Trigger time "<<vCoincidenceTriggerTimes[g]<<"  "<<g<<"  "<<GroupNeighborID<<endl;
                            }
                        }
                   }
               }
          }   
    } 


   //Find the first trigger time
   telData->fTelescopeTriggerTime = 1e10;
   if(telData->bTelescopeHasTriggered)
     {
       //for(int i=0;i<iNumGroups;i++)
       //  cout<<i<<"  "<<vCoincidenceTriggerTimes[i]<<"  "<<vTrigPixsInCluster->at(i).size()<<endl;
       //Sort all trigger times descending
       Int_t* index = new Int_t[iNumGroups];
       TMath::Sort(iNumGroups,vCoincidenceTriggerTimes,index,false);
       int t=0;
       while(vTrigPixsInCluster->at(index[t]).size()==0)
         t++;
        //cout<<vCoincidenceTriggerTimes[index[t]]<<"  "<<vTrigPixsInCluster->at(index[t])[0]<<"  "<<vTrigPixsInCluster->at(index[t])[1]<<endl;
        telData->fTelescopeTriggerTime = vCoincidenceTriggerTimes[index[t]];
        telData->vTriggerCluster = vTrigPixsInCluster->at(index[t]);
        //cout<<"Telescope trigger at "<<telData->fTelescopeTriggerTime<<" trigger bit "<<telData->bTelescopeHasTriggered<<endl;
      }
   delete vTrigPixsInCluster;
   delete [] vCoincidenceTriggerTimes;

   return telData->bTelescopeHasTriggered;
}

//This does the actual triggering
//It returns true if the telescope has triggered and false if not
Bool_t  TriggerTelescopeSPB2::RunTrigger()
{
  // cout<<"TriggerTelescopeSPB2::RunTrigger"<<endl;


   //Run the CFD simulation on all pixels 
   telData->iNumTriggeredTriggerPixels=0;
   for(Int_t i=0;i<iNumTriggerPixels;i++)
    {
      if(RunDiscriminator(i))
	     telData->iNumTriggeredTriggerPixels++;
    }

   //cout<<"Done with the CFD"<<endl<<endl;
   //cout<<telData->iNumTriggeredTriggerPixels<<" Trigger pixel triggered"<<endl;

  //Implementation of the SPB2 Trigger logic
  //First figure out which MUSIC chips have triggered
  FindTriggeredGroups();

  RunCoincidenceLogic();

  

  if(telData->bTelescopeHasTriggered && bDebug)
     {
       cout<<"Event triggered telescope"<<endl; 
       cout<<"at time "<< telData->fTelescopeTriggerTime<<endl;
     }

  return telData->bTelescopeHasTriggered;
}



