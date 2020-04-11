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

class TriggerTelescopeVERITAS : public TriggerTelescopeBase
{

 public:

  TriggerTelescopeVERITAS(ReadConfig *readConfig, int telType, TRandom3 *generator, Bool_t debug = kFALSE, Display *display = NULL);

  void     LoadEvent(TelescopeData *TelData);

  Bool_t   RunTrigger();

  void     SetDiscriminatorRFBConstant(Float_t rfb);
  void     SetDiscriminatorRFBDynamic(Float_t rfb);
  void     SetDiscriminatorRFBUsage(Bool_t rfbuse);
  Float_t  GetDiscRFBDynamicValue(){ return fDiscRFBDynamic; };
  Float_t  GetDiscZeroCrossingRate(){ return lZeroCrossings /(lNumEvents* fTraceLength*1e-3*iNumTriggerPixels); };

           //Defines how many next neighbors are required to trigger the telescope
  void     SetMultiplicity(Int_t multiplicity){ iMultiplicity = multiplicity; };

 protected:

  //is also in base class need to adapt functions in both classes. Which version
  //is used when?
  void  SetParametersFromConfigFile(ReadConfig *readConfig );

  Long_t GetNumZeroCrossings();
  
  Float_t fDiscRFBConstant;            //RFB of the discriminator in pe/MHz
  Float_t fDiscRFBDynamic;             //The dynamic value of the RFB feedback in pe
  Bool_t  bDiscRFBUsage;

  Long_t lNumEvents;                   //holds the number of events filled into traces
  Long_t lZeroCrossings;               //holds the number of zerocrossings for one event counted over all pixels in the camera

  Int_t iMultiplicity;                 //How many groups need to be in a cluster for a trigger  
 
};
