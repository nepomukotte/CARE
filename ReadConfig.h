//! TriggerRead reading a configuration file for the trigger simulation
#ifndef READCONFIG_H
#define READCONFIG_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <TROOT.h>
#include <TRandom3.h>

using namespace std;

class ReadConfig
{

 public:
         
  ReadConfig(TRandom3 *random);
  ~ReadConfig() {};
  Bool_t  ReadConfigFile( string iFile );
  void ReadCommandLine( int argc, char **argv);
  Float_t GetStartSamplingTimeOffsetFromAveragePhotonTime(UInt_t telType){ return fStartSamplingBeforeAverageTime[telType]; };
  Float_t GetTraceLength(UInt_t telType){ return fTraceLength[telType]; };
  Float_t GetSamplingTime(UInt_t telType){ return fSamplingTime[telType]; };

  Bool_t  GetUseSumTrigger(UInt_t telType){ return bUseSumTrigger[telType]; };
  Float_t GetDiscriminatorThreshold(UInt_t telType){ return fDiscThreshold[telType]; };
  Float_t GetDiscriminatorOutputWidth(UInt_t telType){ return fDiscWidth[telType]; };
  Float_t GetDiscriminatorDelay(UInt_t telType){ return fDiscDelay[telType]; };
  Float_t GetDiscriminatorAttenuation(UInt_t telType){ return fDiscAttenuation[telType]; };
  Float_t GetDiscriminatorRFBConstant(UInt_t telType){ return fDiscRFBConstant[telType]; };
  Float_t GetDiscriminatorRFBDynamic(UInt_t telType){ return fDiscRFBDynamic[telType]; };
  Float_t GetDiscriminatorConversionFactormVperPE(UInt_t telType){ return fDiscPEtomVConversion[telType]; };
  Bool_t  GetRFBUsage(UInt_t telType){ return  bDiscUseRFBCircuit[telType]; };
  Bool_t  GetCFDUsage(UInt_t telType){ return  bDiscUseCFD[telType]; };
  Bool_t  GetClippingUsage(UInt_t telType){ return bDoClipping[telType]; };
  Float_t GetClippingLevel(UInt_t telType){ return fClippingLevel[telType]; };
  Bool_t  GetTriggerPatchUsage(UInt_t telType){ return bUsePatches[telType]; };
  vector< vector<int> > GetTriggerPatches(UInt_t telType){ return vPatch[telType]; };        // pattern trigger patches
  
  //Camera configuration
  vector< vector<int> > GetNeighbors(UInt_t telType){ return fNeighbour[telType]; }; 
  UInt_t GetNumberPixels(UInt_t telType){ return fCNChannels[telType]; }; 
  UInt_t GetNumberGroups(UInt_t telType){ return iNumberGroups[telType]; };
  vector< vector<int> > GetNeighborsOfGroup(UInt_t telType){ return fNeighbourGroups[telType]; }; 
  vector< vector<int> > GetMembersOfGroups(UInt_t telType){ return fPixelInGroup[telType]; }; 
  vector<float>         GetXUnrotated(UInt_t telType) const { return fXTube[telType];}   //!< get x-position of tube for current telescop
  vector<float>         GetYUnrotated(UInt_t telType) const { return fYTube[telType];}    //!< get y-position of tube f
  const vector<double>         GetXMM(UInt_t telType) const { return fXTubeMM[telType];}   //!< get x-position of tube for current telescop
  const vector<double>         GetYMM(UInt_t telType) const { return fYTubeMM[telType];}    //!< get y-position of tube f
  vector<float>         GetTubeSize(UInt_t telType) const { return fSizeTube[telType];}                   //!< get tube size in mm
  const vector<double>         GetTubeSizeMM(UInt_t telType) const { return fSizeTubeMM[telType];}                   //!< get tube size in mm
  const vector<double>         GetTubeRotAngle(UInt_t telType) const { return fRotAngle[telType];}         //!< get the rotation angle of the tube
  const vector<int>           GetTubeSides(UInt_t telType) const { return iTubeSides[telType];}  //!< get the number of sides per tube

  Int_t   GetGroupMultiplicity(UInt_t telType){ return iGroupMultiplicity[telType]; };
  Float_t GetFWHMofSinglePEPulse(UInt_t telType){ return fFWHMofSinglePEPulse[telType]; };
  TString  GetNameofHighGainPulseFile(UInt_t telType){ return sHighGainPulseShapeFile[telType]; };
  TString  GetNameofLowGainPulseFile(UInt_t telType){ return sLowGainPulseShapeFile[telType]; };
  Float_t GetSigmaofSinglePEPulseHeightDistribution(UInt_t telType){ return fSigmaSinglePEPulseHeightDistribution[telType]; };
  Float_t GetSampleWidthAveragePulse(UInt_t telType){ return fSampleWidthAveragePulse[telType]; };
  Float_t GetNSBRate(UInt_t telType){ return fNSBRatePerPixel[telType]; };
  Bool_t  GetNSBUsage(){ return bUseNSB; };
  Bool_t  GetAfterPulsingUsage(UInt_t telType){ return  bUseAfterPulsing[telType]; };
  Float_t GetAfterPulsingConstant(UInt_t telType){ return  fAfterPulsingConstant[telType]; };
  Float_t GetAfterPulsingSlope(UInt_t telType){ return  fAfterPulsingSlope[telType]; };


  vector< Float_t > GetRelQE(UInt_t telType);
  vector< Float_t > GetRelGain(UInt_t telType,UInt_t telID, const vector< Float_t > vRelQE);


  Int_t   GetTelescopeMultiplicity(){ return iTelescopeMultiplicity; };
  Bool_t  GetArrayTriggerNextNeighborRequirement(){ return  bArrayTriggerRequiresNextNeighbor; };
  Float_t GetArrayCoincidenceWindow(){ return fArrayCoincidence; };
  vector< vector < int > > GetTelescopeNeighbors(){ return vTelescopeNeighbors; };

  UInt_t  GetRequestedMinNumberOfPhotonsInCamera(UInt_t telType){ return uMinNumPhotonsRequired[telType]; };

  Bool_t  GetBiasCurveBit(){ return bMakeBiasCurve; };
  UInt_t  GetNumberOfTrialsForBiasCurve(){ return uBiasCurveTrials; };
  Float_t GetBiasCurveStartScanRange(){ return fBiasCurveStart; };
  Float_t GetBiasCurveStopScanRange(){ return fBiasCurveStop; };
  Float_t GetBiasCurveStep(){ return fBiasCurveStep; };
  Int_t   GetBiasCurveTelescopeID(){ return iBiasCurveTelescopeID; };  

  Bool_t  GetLoopOverEventsBit(){ return bLoopOverEvents; };
  TString GetSimulatorName(){return sSimulatorName; };
  Int_t   GetAtmosphericModel(){ return iAtmosphericModel; };
  TString GetDayOfSimulatedEvents(){return sDayOfSimulatedEvents; };
  Int_t   GetNumberOfPedestalEvents(){return iNumberPedestalEvents; };
  Int_t   GetNumberOfPedestalEventsToStabilize(){return iNumberPedestalEventsToStabilize; };
  Bool_t  GetVBFwriteBit(){return bWriteVFB; };          


  Int_t   GetFADCSamples(UInt_t telType){return iFADCSamples[telType]; };
  Int_t   GetFADCDynamicRange(UInt_t telType){return iFADCDynamicRange[telType]; };
  Float_t GetFADCSampleTime(UInt_t telType){return fFADCSamplingWidth[telType]; };
  Float_t GetFADCOffsetFromTrigger(UInt_t telType){return fFADCTimeOffsetFromTrigger[telType]; };
  Float_t GetFADCDCtoPEconversionFactor(UInt_t telType){return fFADCDCtoPEconversion[telType]; };
  Float_t GetFADCHiLoGainThreshold(UInt_t telType){return fFADCHiLoGainThreshold[telType]; };
  Float_t GetFADCLowHiGainRatio(UInt_t telType){return fFADCLowHiGainRatio[telType]; };
  Float_t GetFADCHighGainPedestal(UInt_t telType){return fFADCHighGainPedestal[telType]; };
  Float_t GetFADCLowGainPedestal(UInt_t telType){return fFADCLowGainPedestal[telType]; };
  Float_t GetFADCConversionFactormVperPE(UInt_t telType){ return fFADCdctomVConversion[telType]*fFADCDCtoPEconversion[telType]; };

  Float_t GetPileUpWindow(UInt_t telType){ return fPileUpWindow[telType]; };



  vector<Float_t> GetWavelengthsOfQEValues(UInt_t telType){return wl[telType]; };
  vector<Float_t> GetQEValues(UInt_t telType){return qe[telType]; };

  Bool_t            GetSiPMUsage(UInt_t telType){return bSiPM[telType]; };                              //do we use SiPM or not
  vector<Int_t>     GetNumCellsPerSiPM(UInt_t telType){return vNumCellsPerSIPM[telType]; };             //return the number of cells in one SiPM
  vector<Float_t >  GetSiPMOpticalCrosstalk(UInt_t telType){return vSiPMOpticalCrosstalk[telType]; };   //return the optical crosstalk of the SiPM


  Bool_t            GetCrosstalkUsage(UInt_t telType){return bCrosstalk[telType]; };                     //Use crosstalk between pixel
  Float_t           GetCrosstalkValue(UInt_t telType){return fCrosstalkValue[telType]; };  

  Int_t   GetNumberOfTelescopes(){return iNumberOfTelescopes; };
  Int_t   GetNumberOfTelescopeTypes(){return iNumberOfTelescopeTypes; };
  Int_t   GetTelescopeType(UInt_t telID){return iTelType[telID]; };
  Int_t   GetTelescopeIDinSuperArray(UInt_t telID){return iTelIDInSuperArray[telID]; };
  Float_t GetWinstonConeEfficiency(UInt_t telID){ return fWinstonConeEfficiency[telID]; };
  Float_t GetRelativeTelescopeGain(UInt_t telID){ return fRelativeTelescopeGain[telID]; };
  Float_t GetSigmaElectronicNoise(UInt_t telID){ return fSigmaElectronicNoise[telID]; };
  Bool_t  GetOpticalPSFBlurBit(UInt_t telID){return bBlurPSF[telID];};
  Float_t GetOpticalPSFBlurSigma(UInt_t telID){return fBlurSigma[telID];};
  Float_t GetTransitTimeSpread(UInt_t telID){return fTransitTimeSpread[telID];};



 protected:

  void resetCamVectors();
  void resetNeighbourGroupLists();
  void resetTelVectors();
  void resetTelTypeVectors();
 
  void convertMMtoDeg();

  void ReadLine(string iline, std::ifstream *inFileStream);

  Bool_t fDebug;

  TRandom3 *rand;

  //Telescope trigger configuration
  vector<Bool_t>  bUseSumTrigger;              //Sum pixels before discriminator
  vector<Float_t> fClippingLevel;              //The level in mV at which the signals are clipped 
  vector<Bool_t>  bDoClipping;                 //Do we clip the signals before summing
  vector<Float_t> fDiscThreshold;              //Discriminator threshold of pixel
  vector<Float_t> fDiscWidth;                  //Width of Discriminator output
  vector<Float_t> fDiscDelay;                  //Delay of the inverted signal in the CFD
  vector<Float_t> fDiscAttenuation;            //Attenuation of the non-inverted signal in the CFD
  vector<Float_t> fDiscRFBConstant;            //The constant in the RFB feedback units pe/MHz
  vector<Float_t> fDiscPEtomVConversion;       //The conversion factor at the input of the Discriminator mV per pe.
  vector<Float_t> fFADCdctomVConversion;       //The conversion factor at the input of the FADC mV per pe.
  vector<Float_t> fDiscRFBDynamic;             //The value in the RFB feedback in units pe applied as offset.
                                               //If the RFB circuit is used this is just a start value 
  vector<Bool_t>  bDiscUseCFD;                 //Do we use the CFD part of the discriminator
  vector<Bool_t>  bDiscUseRFBCircuit;          //Do we use the RFB circuit
  vector<Int_t>   iGroupMultiplicity;          //How many groups need to be in a cluster for a trigger  

  vector<Bool_t>  bUsePatches;                 //Is the trigger topology divided into patches
  vector<vector< vector<int> > > vPatch;        // pattern trigger patches


  //Trace related variables
  vector<Float_t> fFWHMofSinglePEPulse;        //The full width at hald maximum of the single pe pulse
  vector<TString> sHighGainPulseShapeFile;         //name of the file that stores the single pe pulse shape
  vector<TString> sLowGainPulseShapeFile;  //name of the file that stores the single pe pulse shape
  vector<Float_t> fSigmaSinglePEPulseHeightDistribution; //The sigma of the single PE pulse height distribution
  vector<Float_t> fSampleWidthAveragePulse;    //The sample width used in the average single pe pulse
  vector<Float_t> fNSBRatePerPixel;            //the NSB rate per pixel in the focal plane;
  Bool_t  bUseNSB;                     //set to true if we want to use NSB in the simulation
  vector<Bool_t>  bUseAfterPulsing;            //Do we simulate Afterpulsing: true yes false else
  vector<Float_t> fAfterPulsingConstant;       //Constant of a fit to the rate vs. threshold curve of a single pe
                                       //Fit function exp(a+b*x), where a=constant b=slope
  vector<Float_t> fAfterPulsingSlope;          //Slope of a fit to the rate vs. threshold curve of a single pe
                                       //Fit function exp(a+b*x), where a=constant b=slope
  vector<Float_t>         fSamplingTime;       //The sampling rate or resolution of the simulated trace
  vector<Float_t>         fTraceLength;        //the length of the simulated trace per group
  vector<Float_t>         fStartSamplingBeforeAverageTime;   //start of sampling the trace before the average photon arrival time
  vector<Float_t>         fPileUpWindow;   //Half width of the window used to determine the right pulse shape

  //Array trigger configuration
  Int_t   iTelescopeMultiplicity;      //How many telescopes need to be in a cluster for a trigger
  Bool_t  bArrayTriggerRequiresNextNeighbor;  //Does the array trigger require next neighbor coincidence 
  Float_t fArrayCoincidence;           //Coincidence window on array trigger level
  vector< vector< int >  >   vTelescopeNeighbors; //Has the neighbors of each telescope which are considered in an array trigger

  //Biascurve parameters
  Bool_t  bMakeBiasCurve;              //Set to 1 if we want to make a bias curve
  UInt_t  uBiasCurveTrials;            //Trials that go into the production of the bias curve
  Float_t fBiasCurveStart;             //Start of the discriminator scan range for the bias curve
  Float_t fBiasCurveStop;              //Stop of the discriminator scan range for the bias curve
  Float_t fBiasCurveStep;              //Step in the scan of the discriminator values
  Int_t   iBiasCurveTelescopeID;

  //Looping over events
  Bool_t  bLoopOverEvents;
  Int_t   iNumberPedestalEvents;              //The number of pedestal events that will be simulated
  Int_t   iNumberPedestalEventsToStabilize;   //The number of pedestal events that will be simulated
                                              //to stabilize the discriminator

  //Cherenkovphoton throughput
  vector<Float_t> fWinstonConeEfficiency;                 //The efficiency of the Winstoncone

  //Relative gain of telescopes
  vector<Float_t> fRelativeTelescopeGain;                 //The relative gain of each telescope

  vector<Float_t> fSigmaElectronicNoise;                 //The electronic noise in the readout chain

  //Quantum/PDE of the Photon detectors
  vector<vector<Float_t> > wl;
  vector<vector<Float_t> > qe;
  vector< Bool_t >  bFlatfieldCamera;                 //Flatfield the camera response
  vector< Float_t >  fGainSigma;                       //sigma of the gain distribution
  vector< Float_t >  fQESigma;                         //sigma of the QE distribution

  vector< UInt_t >  uMinNumPhotonsRequired;


  //Optical PSF bluring 
  vector< Bool_t >  bBlurPSF;
  vector<Float_t> fBlurSigma;                    //sigma in mm by which the optical PSF is blured furthery of the Winstoncone

  //PMT/SiPM Transit time spread
  vector<Float_t> fTransitTimeSpread;           //time spread of the photelectrons making it through the sensor. RMS in ns.
 
  
  vector< Bool_t >  bCrosstalk;                     //Use crosstalk between pixel
  vector< Float_t > fCrosstalkValue;  


  //SiPM related variables
  vector< Bool_t > bSiPM;                                    //flag to figure out if we use SiPMs
  vector< vector<Int_t> > vNumCellsPerSIPM;                  //the number of cells in one SiPM
  vector< vector<Float_t > > vSiPMOpticalCrosstalk;          //the optical crosstalk of the SiPM
 

  //VBF related variables and configs
  Bool_t  bWriteVFB;                   //Is a VBF file written
  TString sSimulatorName;              //Name of the person executing the simulation
  Int_t   iAtmosphericModel;           //The atmospheric model used in the simulation
  TString sDayOfSimulatedEvents;       //The datum that will be attached to each simulated event

  //FADC parameters
  vector<Int_t>   iFADCSamples;                //The number of FADC samples
  vector<Float_t> fFADCSamplingWidth;          //The sampling width of the FADC
  vector<Float_t> fFADCDCtoPEconversion;           //The DC to PE conversion
  vector<Int_t>   iFADCDynamicRange;           //The dynamic range of the FADC
  vector<Float_t> fFADCTimeOffsetFromTrigger;  //The offset between the trigger time and the start of the readout window
  vector<Float_t> fFADCHiLoGainThreshold;      //The Threshold in FADC counts at which switching to the low gain is activated
  vector<Float_t> fFADCLowHiGainRatio;         //Gain ratio logain/higain
  vector<Float_t> fFADCHighGainPedestal;        //The FADC high gain pedestal in dc counts
  vector<Float_t> fFADCLowGainPedestal;         //The FADC low gain pedestal in dc counts


  Int_t  iNumberOfTelescopes;          //The number of Telescopes in the array
  Int_t  iNumberOfTelescopeTypes;      //The number of different Telescope types in the array

  //Camera configuration 
  vector< unsigned int > fCNChannels;               //!< number of channels
  vector< unsigned int > iNumberGroups;             //the number of groups of summed pixel
  vector<float> fMirFocalLength;                    //the focal length of the mirror in m
  
  vector< vector<float> > fXTube;                   //!< x-position of tube in [deg] (in camera coordinates)
  vector< vector<float> > fYTube;                   //!< y-position of tube in [deg] (in camera coordiantes)
  vector< vector<float> > fSizeTube;                //!< tube radius in [deg]
 
  vector< vector<double> > fXTubeMM;                 //!< x-position of tube in [mm]
  vector< vector<double> > fYTubeMM;                 //!< y-position of tube in [mm]
  vector< vector<double> > fSizeTubeMM;              //!< tube radius in [mm]
  vector< vector<double> > fRotAngle;                //!< tube rotation angle in rad (converted from deg when reading in)
  vector< vector<int> > iTubeSides;                 //!< the number of sides the pixel has
  vector< vector< vector<int> > > fNeighbour;       //!< neighbour identifier

  vector< vector< vector<int> > > fPixelInGroup;    //!< The pixel that are in one group
  vector< vector< vector<int> > > fNeighbourGroups; //!< neighbour identifier

  vector< int > iTelIDInSuperArray;
  vector< int > iTelType;

};
    
#endif
