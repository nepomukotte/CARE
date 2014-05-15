#include <TROOT.h>
#include <TRandom3.h>
#include <TH1F.h>		
#include <vector>

#include "TelescopeData.h"
#include "TraceGenerator.h"
#include "ReadConfig.h"
#include "Display.h"

using namespace std;

class FADC {

 public:
 
  FADC(ReadConfig *readConfig, TraceGenerator *traceGenerator, int telType, TRandom3 *generator,Bool_t debug,Display *display);

  //Event handling
  void RunFADC(TelescopeData *telData);

  void SetDebugInfo(Float_t energy, Int_t telid, Float_t zenith, Float_t azimuth);

  void     PrintHowOftenTheTraceWasTooShort();



 protected:

  void DigitizePixel( Int_t PixelID );

  void     SetParametersFromConfigFile(ReadConfig *readConfig );

  //Variables

  Bool_t bDebug;

  Display *debugDisplay;

  TRandom3 *rand;                                  //Our random number generator

  TelescopeData  *telData;
  TraceGenerator *tracegenerator;

  Int_t iTelType;

  Int_t iNumPixels;                                //Number of pixels in the camera 

  Int_t iNumTimesOutOfAnalogueTraceUpperEnd;
  Int_t iNumTimesOutOfAnalogueTraceLowerEnd;

  //Single PE pulse shape parameters
  Float_t fDCtoPEconversion;                       //the gain calibration constant DC per PE; defined for the amplitude of a single pe pulse

  //FADC parameters
  Float_t fFADCSamplingWidth;                      //Sampling time
  Int_t   iFADCSamples;                              //Number of samples per trace
  Float_t fOffset;                                 //time offset of the trigger from the beginning of the readout window
  Int_t   iDynamicRange;                           //the count range of the digitizer
  Float_t fFADCconversion;                         //the conversion of mV to DC
  Float_t fHiLoGainThreshold;                      //dc counts at which the HiLoGain switch is activated
  Float_t fLowHiGainRatio;                         //Gain ratio logain/higain
  Float_t fHighGainPedestal;                       //High Gain Pedestal offset in the FADC
  Float_t fLowGainPedestal;                        //Low Gain Pedestal offset in the FADC

  //Event dependent parameters
  Float_t fTimeStartFirstSample;
   
  Float_t fTraceSamplingTime;                     //The sampling time steps for the analog trace
  Float_t fTraceLength;                           //The total length of an analog trace
  Float_t fStartSamplingBeforeAverageTime;        //the time before the average photon arrival time when the trace gets started to be sampled

  //Make Function that allows to scan the position of the Cherenkov pulse in the Trace
  Float_t fenergy;
  Int_t ftelid;
  Float_t fzenith;
  Float_t fazimuth;

};
