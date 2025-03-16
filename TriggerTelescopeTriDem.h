#include "TriggerTelescopeBase.h"

#include <TROOT.h>
#include <vector>
#include <TH1F.h>
#include <TRandom3.h>
#include <TVectorT.h>


#include "TelescopeData.h"
#include "TraceGenerator.h"
#include "ReadConfig.h"
#include "Display.h"

using namespace std;

class TriggerTelescopeTriDem : public TriggerTelescopeBase
{

 public:

  TriggerTelescopeTriDem(ReadConfig *readConfig, int telType, TRandom3 *generator, Bool_t debug = kFALSE, Display *display = NULL);
  ~TriggerTelescopeTriDem();

  Bool_t   RunTrigger();

 private:
   
  void FindTriggeredGroups(); 
  bool FindTrigger(); 
 
  Int_t iNumGroups;                  //all the MUSIC chips in the camera group = one music
  vector< vector<int> > iGroupNeighborIDs;          //holds the neighbors of each group;
  vector< vector<int> > iGroupMembers;            //holds the members of each group;
  Int_t iNumPixels;                  //all the pixels in the camera

  vector<vector<int> > *vTriggeredPixelInGroup;        //pixel in each group that triggered 
  vector<vector<float> > *vTriggerTimesInGroup;        //times when discriminator fired for each pixel in a group
};
