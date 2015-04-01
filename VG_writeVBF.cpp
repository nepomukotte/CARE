#include "VBF/VBankFileWriter.h"
#include "VBF/VPacket.h"
#include "VBF/VArrayEvent.h"
#include "VBF/VDatum.h"
#include "VBF/VEventType.h"
#include "VBF/VSimulationHeader.h"
#include "VBF/VSimulationData.h"
#include "VBF/VCorsikaSimulationData.h"

#include "VBF/VConfigMaskUtil.h"
#include "VATime.h"

using namespace VConfigMaskUtil;

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <time.h>
#include <vector>

#include "TROOT.h"
#include "TRandom3.h"
#include "VG_writeVBF.h"

#define DEBUG(x) std::cerr << #x << " = " << x << std::endl

VG_writeVBF::VG_writeVBF(const int &numTel,
                         const std::vector<unsigned> &numPixTel,
                         const std::string &vbfFileName,
                         const long &vbfRunNumber,
                         const int &debugLevel):
  fnbr_tel(numTel),
  fnumPixV(numPixTel),
  fVbfFileName(vbfFileName),
  fRunNumber(vbfRunNumber),
  fDebugLevel(debugLevel)
{

  if (fnbr_tel ==0) {
    std::cerr << "vbf constructor: number of telescope = 0" << std::endl;
    std:: cerr << "   shutting down" << std::endl;
    exit(0);
  }

  if (fnumPixV.size() != fnbr_tel) {
    std::cerr << "vbf constructor: number of pixel totals: " 
              << fnumPixV.size() << " !=  number of telescopes"
              << std::endl;
    std:: cerr << "   shutting down" << std::endl;
    exit(0);
  }
  
  initialize();

 // make the configuration mask (brute force)
  std::ostringstream os;
  for (unsigned i=0;i<fnbr_tel;i++) {
    os << i << ",";    // may have an extra comma? at the end, works ok though
  }
  
  fConfigMaskAll = os.str();  // for all telescopes
  fCMaskAll = (uint8_t)toDAQMask(parseConfigMask(fConfigMaskAll.c_str()));

  if (fDebugLevel > 0) {
    std::cerr << "************* in constructor ***********" << std::endl;
    std::cerr << "  fnbr_tel: " << fnbr_tel << std::endl;
    std::cerr << "  number of pixels in each telescope" << std::endl;
    for (unsigned tel=0;tel<fnbr_tel;tel++) {
      std::cerr << "    " << fnumPixV[tel] << "  ";

    }
    std::cerr << std::endl;

    std::cerr << " parameters for VBankFileWriter  " << std::endl;
    std::cerr <<"    fVbfFileName:    " << fVbfFileName << std::endl;
    std::cerr <<"    fRunNumber:     " << fRunNumber << std::endl;
    std::cerr <<"    fConfigMaskAll: " << fConfigMaskAll << std::endl;
    std::cerr <<"  writing VBankFileWriter " << std::endl;
  }

  // open vbf file for writing
  pfWriter = new VBankFileWriter(fVbfFileName,fRunNumber,
                                 parseConfigMask(fConfigMaskAll.c_str()));

  if (pfWriter==NULL){
    printf(" vbf file failed to open: %s\n",fVbfFileName.c_str());
    exit(EXIT_FAILURE);
  }
  
}

//********************** initialize() ********************
void VG_writeVBF::initialize() {

  frand = new TRandom3();
  setSeed();   // set to default seed

  fCorsikaParticleIDType = 0; // initialize to kascade ID type

  fFirstArgon = true;
  fEventNumber = 0;
  fPedEventCount = 0;

  fPedFlag = true;
  fPedFileName = "";
  fPedTimeInterval = 1;
  fFirstValidEventTimeStr = "2006-11-28 02:00:00.000000000 GPS";
  fEventRatePerSec = 100.0;
  fMeanTimeBetweenEventsSec = 1/fEventRatePerSec;

  fEventTime = 0;
  fPedEventTime = 0;
  packet = 0;
  ae  = 0;   
  event = 0;
  at = 0;
  simu_data = 0;
  cors_data = 0;
  
  fEvType = new VEventType();
  fEvType->setNewStyleCode(1);

  fnumFADCSamples = 24;
  fnumPmts = 0; 
  ftriggered_readout = false;
  fnum_vevent = 0;

  // initializations for VERITAS, have to change for different cameras
  fmaxNumChannels = 500;  // for veritas

  // create vector containing all pixel locations
  fTelPixV = new std::vector<  std::vector<VPixelLocation> >(fnbr_tel);

  // reserve sufficient space for ease in addressing vector later
  for (unsigned tel=0;tel<fnbr_tel;tel++) {
    (*fTelPixV)[tel].resize(fnumPixV[tel]);
  }

  // reserve sufficient space in telescope vectors
  fTelSouth.resize(fnbr_tel);
  fTelEast.resize(fnbr_tel);
  fTelUp.resize(fnbr_tel);
  
  aztel.resize(fnbr_tel);
  eltel.resize(fnbr_tel);

}

//************************** setTelescopeLocations *************************************
void VG_writeVBF::setTelescopeLocations(const int & telId, 
                                           const float & telSouthLoc, 
                             const float & telEastLoc,
                             const float & telUpLoc) {
  fTelSouth[telId] = telSouthLoc;
  fTelEast[telId]  = telEastLoc;
  fTelUp[telId]    = telUpLoc;
  
}

//**************************** setPixelLocations ***********************************
void VG_writeVBF::setPixelLocations(const int &telId, const int &pixId,
                                       const float &pixEastLocDeg,
                                       const float &pixUpLocDeg,
                                       const float &pixRadDeg) {
  (*fTelPixV)[telId][pixId].fPixLocEastAtStowDeg = pixEastLocDeg;
  (*fTelPixV)[telId][pixId].fPixLocUpAtStowDeg = pixUpLocDeg;
  (*fTelPixV)[telId][pixId].fPixRadiusDeg      = pixRadDeg;

}

//******************* makeSimulationHeader ********************************************
bool VG_writeVBF::makeSimulationHeader(std::string &configFilePlus,
                                          const std::string &simulator,
                                          const float &obserAltitudeM,
                                          const uint32_t &simulationPackage,
                                          const uint32_t &atmosphericModel) {

  fSimConfigFile = configFilePlus;
  std::string fSimulator =  simulator;
  uword32 fSimulationPackage = simulationPackage;
  uword32 fAtmosphericModel = atmosphericModel;
  fObsAlt = obserAltitudeM;
 // create and store simulation header in event 0.
 
  // use new packet variable, don't have to call make_packet 
  VPacket *packet=new VPacket();

  // Get present time
  uint32_t y,m,d;
  time_t t1 = time(NULL);
  tm *now = localtime(&t1);
  y = now->tm_year + 1900;
  m = now->tm_mon + 1;
  d = now->tm_mday;
  uword32 fDateOfSimsUTC = y*10000ULL + m*100ULL + d;
  uword32 fDateOfArrayForSims = fDateOfSimsUTC;

 // create VArrayConfiguration vector
  std::vector<VArrayConfiguration> fArray;
  for (unsigned tel=0;tel<fnbr_tel;tel++) {
    
    VArrayConfiguration fArTel;
    fArTel.fRelTelLocSouthM = fTelSouth[tel];
    fArTel.fRelTelLocEastM  = fTelEast[tel];
    fArTel.fRelTelLocUpM    = fTelUp[tel];

    fArTel.fCamera = (*fTelPixV)[tel];
    fArray.push_back(fArTel);
  }
  

  if (fDebugLevel>0) {
    std::cerr << "************* in makeSimulationHeader *******" 
              << std::endl;
    std::cerr << "   parameters for VSimulationHeader Constructor" 
              << std::endl;
    std::cerr << "     fDateOfSimsUTC:      " << fDateOfSimsUTC 
              << std::endl;
    std::cerr << "     fSimulationPackage:  " << fSimulationPackage 
              << std::endl;
    std::cerr << "     fSimulator:          " << fSimulator << std::endl;
    std::cerr << "     fDateOfArrayForSims: " << fDateOfArrayForSims 
              << std::endl;
    std::cerr << "     fObsAlt:             " << fObsAlt << std::endl;
    std::cerr << "     fSimConfigFile:      " << fSimConfigFile 
              << std::endl;
    std::cerr << "print fArray with fDebugLevel >= 2" << std::endl;
  } 
  
  if (fDebugLevel > 1) {
    std::cerr << "**************in makeSimulationHeader *******" 
         << std::endl; 
    std::cerr << "      printing fArray, this may take some space" 
              << std::endl;
    
    std::cerr << "       telescope locations south east up" << std::endl 
              << std::endl;
    for (unsigned tel = 0; tel < fnbr_tel;tel++) {
      
      std::cerr << tel << " " << fArray[tel].fRelTelLocSouthM << " " 
                << fArray[tel].fRelTelLocEastM << " " 
                << fArray[tel].fRelTelLocUpM << std::endl << std::endl;
    }
    std::cerr << std::endl;
    
    std::cerr << "     pixel locations, telescope, East,Up Radius"
              << std::endl;
  
    for (unsigned tel = 0; tel < fnbr_tel;tel++) {
      for (unsigned pix = 0;pix < fnumPixV[tel];pix++) {
        float east = fArray[tel].fCamera[pix].fPixLocEastAtStowDeg;
        float up = fArray[tel].fCamera[pix].fPixLocEastAtStowDeg;
        float rad = fArray[tel].fCamera[pix].fPixLocEastAtStowDeg;
        std::cerr << tel << " " << pix << " " 
                  << east << " " << up << " " << rad << std::endl;
      }
      std::cerr << std::endl;
    }
  }
   


  
  VSimulationHeader * pfSimHead =
    new VSimulationHeader(fDateOfSimsUTC,fSimulationPackage,
                          fSimulator,fDateOfArrayForSims,
                          fAtmosphericModel,fObsAlt,fArray,
                          fSimConfigFile);

  // Put the simulation header data into the packet
  packet->put(VGetSimulationHeaderBankName(), pfSimHead);
  if (!packet->has(VGetSimulationHeaderBankName())  )  {
    showXErrorVbfWriter("VG_writeVBF: No SimulationHeader bank in packet \
      but we just put one in");
 }

  // finally, write the packet into the file
  
  pfWriter->writePacket(fEventNumber, packet);

  // dispose of the packet, so that we don't leak memory
  delete packet;
  
  fEventNumber++;  // increment event number

  return true;
}

/************************ setPedestalParameters **********************************/
void VG_writeVBF::setPedestalParameters(const bool &makePedsFlag,
                            const int &pedTimeIntervalS,
                                            const std::string &pedFileName) {
   fPedFlag = makePedsFlag;
   fPedTimeInterval = pedTimeIntervalS;
   fPedFileName = pedFileName;

   if (fDebugLevel > 0) {
     std::cerr << "******** in setPedestalParameters ********" << std::endl;
     DEBUG(fPedFlag); DEBUG(fPedTimeInterval); DEBUG(fPedFileName); 
   }
 }

/************************** setFirstEventTimeString ********************************/
void VG_writeVBF::setFirstEventTimeString(const std::string &firstEventTime) {

  fFirstValidEventTimeStr = firstEventTime;
  if (fDebugLevel > 0) {
    std::cerr << "********** in setFirstEventTimeString ******"
              << std::endl;
    DEBUG(fFirstValidEventTimeStr);
    std::cerr << std::endl;
  }
}

/************************* setShowerFrequencyHz *********************************/
void VG_writeVBF::setShowerFrequencyHz(const double &showerfreq) {

  fEventRatePerSec = showerfreq;
  fMeanTimeBetweenEventsSec = 1.0/fEventRatePerSec;

  if (fDebugLevel > 0) {
    std::cerr << "******** in setShowerFrequencyHz ***** " << std::endl;
    DEBUG(fEventRatePerSec);
  }
}

/************************* makeEventTimeClasses *********************************/
bool VG_writeVBF::makeEventTimeClasses() {

    // set initial event time
  fEventTime = new VATime();
  fEventTime->setFromString(fFirstValidEventTimeStr);
  
  if (fDebugLevel > 0) {
    std::cerr << "********** in makeEventTimeClasses ***********" << std::endl;
    std::cerr << "    fPedFlag: " << fPedFlag << std::endl;

    std::cerr << "       calling test_vatime for fEventTime" << std::endl;
    testVATime(fEventTime);
  }

  if (fPedFlag) {
    fPedEventTime = new VATime();
    fPedEventTime->setFromString(fFirstValidEventTimeStr);
   // make sure ns = 0
    uint32_t fYear,fMonth,fDay,H,M,S,NS;
    
    fPedEventTime->getCalendarDate(fYear,fMonth,fDay);
    fPedEventTime->getTime(H,M,S,NS);
    NS=0;   //On the tick!
    fPedEventTime->setFromCalendarDateAndTime(fYear,fMonth,fDay,
                                              H,M,S,NS);
    if (fDebugLevel > 0) {
      std::cerr << " pedestal event time " << std::endl;
      testVATime(fPedEventTime);
    }

  }
    
  // this ensures that the first event will be a pedestal event
  incrementEventTime();

    return true;
}

/************************* setAzimElevTelDeg *********************************/
void VG_writeVBF::setAzimElevTelDeg(const unsigned &telId, 
                                    const float &azimDeg,
                                      const float &elevDeg) {

 aztel[telId] = azimDeg;
 eltel[telId] = elevDeg;
}

/************************ makePacket **********************************/
bool VG_writeVBF::makePacket() {

  if (fDebugLevel > 0) {
    std::cerr << "*********in makePacket ******** " << std::endl;
    std::cerr << "     fpedestalevent: " << fpedestalevent 
              << std::endl;
  }


  if (packet !=0) {
   showXErrorVbfWriter("attempting to make packet with packet != 0");
  }

  packet = new VPacket();
  
  // initialize parameters
  fCMaskTrig = 0;

  // set time for packet, depends on pedestal or data event
  // all events and array trigger will have this time.

  // set event type and GPS variables for later use.
  if (fpedestalevent) {
    fPedEventTime->getForVBF(fGPSYear,5,fGPSWords);
    fEvType->trigger = VEventType::PED_TRIGGER;
  }
  else {
    fEventTime->getForVBF(fGPSYear,5,fGPSWords);
    fEvType->trigger = VEventType::L2_TRIGGER;
  }
  // time will be incremented in store_packet.

  return true;

}

 /************************* storePacket *********************************/
bool VG_writeVBF::storePacket() {

  if (packet==0) {
    showXErrorVbfWriter(" can't write packet, packet = 0"); 
  }
  
  pfWriter->writePacket(fEventNumber,packet);
  
  delete packet;
  packet = 0; // set pointers to zero
  ae = 0;
  at = 0;
  
  fEventNumber++;
  if (fpedestalevent) {
    incrementPedTime();
    fPedEventCount++;
  }
  else {
    incrementEventTime();
  }
  
  if (events.size()) events.clear();
  
  return true;
}

/********************** makeEvent ************************************/
bool VG_writeVBF::makeEvent() {

  if (fDebugLevel > 0) {
    std::cerr << "******** in make_event *******" << std::endl;
    DEBUG(fEventNumber);
    DEBUG(fcurrent_tel);
    DEBUG(fnumFADCSamples);
    DEBUG(fmaxNumChannels);
  }

  events.push_back(new VEvent()); 
  event = events.back();

  uint16_t rct = 1;
  event->resizeClockTrigData(rct);
  event->setEventNumber(fEventNumber);
  event->setCompressedBit(true);
  event->setNodeNumber((uint8_t)fcurrent_tel);

  
  if (fDebugLevel > 0) {
    std::cerr << "TRIGGERMASKEVENT: fCMaskTrig " 
              << (unsigned)fCMaskTrig << std:: endl;
  }
  event->setTriggerMask(fCMaskTrig);
  
  // initialize clock trigger boards
  for (unsigned k=0;k<event->getNumClockTrigBoards();k++) {
    for (unsigned l=0;l<7;l++) {
      event->getClockTrigData(k)[l]=0;
    }
  }
  
  event->resizeChannelData((uint16_t)fnumFADCSamples,
                           (uint16_t)fmaxNumChannels);

  event->resizeChannelBits((uint16_t)fmaxNumChannels);


  event->setFlags(1);  // enable compression

  // store event type and store GPS times
  event->setEventType(*fEvType);     // set in make_packet

  event->getGPSTime()[0] = fGPSWords[0];
  event->getGPSTime()[1] = fGPSWords[1];
  event->getGPSTime()[2] = fGPSWords[2];
  event->getGPSTime()[3] = fGPSWords[3];
  event->getGPSTime()[4] = fGPSWords[4];
  event->setGPSYear(fGPSYear);

  // initialize event pixel variables, same for peds and data
  for (unsigned k=0;k<fmaxNumChannels;k++) {
    event->setTriggerBit(k,false);  // no channels triggered
    event->setHitBit(k,true);       // all passed zero suppression
    setPixel( (int)k);
    setChargePedHigain();

    // zero fadc samples also
    for (unsigned fs=0;fs<fnumFADCSamples;fs++) {
      int samnum = (int)fs;
      int initzero = 0;
      storeSample(samnum,initzero);
    }
  }
  
  return true;
}

/*********************** storeEvent ***********************************/
bool VG_writeVBF::storeEvent() {

  if (fDebugLevel > 0) {
    std::cerr << "******** in store_event *******" << std::endl;
    
  }

  if (ae!=0) {
    ae->addEvent(event);
  }
  else {
    showXErrorVbfWriter("no arrayevent in store_event");
  }
  // increment count of stored events
  fnum_vevent++;
  
  return true;
}

/************************* makeArrayEvent *********************************/
bool VG_writeVBF::makeArrayEvent() {
  
  if(fDebugLevel > 0) {
    std::cerr << "************ in make_arrayevent() *****" << std::endl;
    DEBUG(fEventNumber);
  }

  if (events.size()) {
    std::cerr << "   events.size non-zero, in make_arrayevent " <<
      fEventNumber << std::endl;
    events.clear();
  }
  fnum_vevent = 0;

  if (ae == 0) {
    ae = new VArrayEvent(fRunNumber);
  }
  else {
    showXErrorVbfWriter("trying to create arrayevent with ae !=0");
  }
  return true;
}

/**************************** storeArrayEvent ******************************/
bool VG_writeVBF::storeArrayEvent() {

  packet->putArrayEvent(ae);

  if (!packet->hasArrayEvent()) {
    showXErrorVbfWriter("packet has no arrayevent even thou one was written");
  }
  
  return true;
}

/*************************** storeSample *******************************/
void VG_writeVBF::storeSample(int &ti,int &iDC) {

  // store time slice in event.
  if (fDebugLevel > 3) {
    std::cerr << "********* store_sample tel pix ti iDC:  " 
              << fcurrent_tel << " " << fcurrent_pix 
              << " " << ti << " " 
              << iDC << std::endl;
  }

  unsigned uti = (unsigned) ti;
  uint8_t iDC_unsign = (uint8_t) iDC;
  event->setSample(fcurrent_pix, uti,iDC_unsign);
 
}

/*************************** setChargePedHigain *******************************/
void VG_writeVBF::setChargePedHigain() {

  uword16 chped = 0;
  event->setCharge(fcurrent_pix,chped);
  event->setPedestal(fcurrent_pix,chped);
  event->setHiLo(fcurrent_pix,false);  // initialize HiLo gain here also
  
}

/**************************** makeVSimulationData1 ******************************/
bool VG_writeVBF::makeVSimulationData1() {

  float fEnergyGeV;             // parameters for VSimulationData
  float fObservationZenithDeg;
  float fObservationAzimuthDeg;
  float fPrimaryZenithDeg;
  float fPrimaryAzimuthDeg;
  float fCoreEastM;
  float fCoreSouthM;
  float fCoreElevationMASL;
  uword32 fCORSIKAParticleID;
  
  // fdim computes the positive difference between the two arguments
  fObservationZenithDeg  = fdimf(90.0,eltel[0]);
  fObservationAzimuthDeg = aztel[0];
  
  // set these to zero for now for pedestal events
  
  fEnergyGeV=0.0;             // parameters for VSimulationData
  fPrimaryZenithDeg=0.0;
  fPrimaryAzimuthDeg=0.0;
  fCoreEastM=0.0;
  fCoreSouthM=0.0;
  fCoreElevationMASL=0.0;
  fCORSIKAParticleID=1;

  // kludges
  float fRefZenithDeg=fPrimaryZenithDeg; //assumes point source
  float fRefAzimuthDeg=fPrimaryAzimuthDeg; //assumes point sourcex
  float fRefPositionAngleDeg=0;
  
  
  simu_data = new VSimulationData(fCORSIKAParticleID,fEnergyGeV,
				  fObservationZenithDeg, fObservationAzimuthDeg,
				  fPrimaryZenithDeg, fPrimaryAzimuthDeg, 
                                  fRefZenithDeg, 
				  fRefAzimuthDeg, fRefPositionAngleDeg,
				  fCoreEastM,fCoreSouthM,fCoreElevationMASL);
  packet->putSimulationData(simu_data);
  
  return true;
}

/**************************** makeVSimulationData ******************************/
 bool VG_writeVBF::makeVSimulationData(const float &energyGeV,
                          const float &primaryZenithAngleDeg,
                          const float &primaryAzimuthDeg,
                          const float &primaryCoreEastM,
                          const float &primaryCoreSouthM,
                          const float &primaryCoreElevMASL,
                          const float &firstInteractionHeight, 
                          const float &firstInteractionDepth, 
                          const uint32_t &corsikaParticleID,
                          const u_int32_t &corsikaRunNumber,
 			  const u_int32_t &fShowerID) {

  float fEnergyGeV;             // parameters for VSimulationData
  float fObservationZenithDeg;
  float fObservationAzimuthDeg;
  float fPrimaryZenithDeg;
  float fPrimaryAzimuthDeg;
  float fCoreEastM;
  float fCoreSouthM;
  float fCoreElevationMASL;
  float fFirstInteractionHeight;
  float fFirstInteractionDepth;
  uword32 fCORSIKAParticleID;
  u_int32_t uShowerID;
  u_int32_t uCorsikaRunNumber; 


  fEnergyGeV              = energyGeV;              // 0.0 default
  fPrimaryZenithDeg       = primaryZenithAngleDeg;  // 0.0 default
  fPrimaryAzimuthDeg      = primaryAzimuthDeg;      // 0.0 default
  fCoreEastM              = primaryCoreEastM;       // 0.0 default
  fCoreSouthM             = primaryCoreSouthM;      // 0.0 default
  fCoreElevationMASL      = primaryCoreElevMASL;    // 0.0 default
  fFirstInteractionHeight = firstInteractionHeight; // 0.0 default
  fFirstInteractionDepth  = firstInteractionDepth;  // 0.0 default
  fCORSIKAParticleID      = corsikaParticleID;      // 1   default
  uShowerID               = fShowerID;              // 0   default
  uCorsikaRunNumber       = corsikaRunNumber;       // 0   default



  // if fCoreElevationMASL < 0.0, set to observatory height
  if (fCoreElevationMASL < 0.0) {
    fCoreElevationMASL = fObsAlt; 
  }
    
  fObservationZenithDeg  = 90.0 - eltel[0];
  fObservationAzimuthDeg = aztel[0];
  
  if (fCorsikaParticleIDType == 0) {
    fCORSIKAParticleID = KascadeType2CorsikaType(fCORSIKAParticleID);
  }

  float fRefZenithDeg=fPrimaryZenithDeg; //assumes point source
  float fRefAzimuthDeg=fPrimaryAzimuthDeg; //assumes point sourcex
  float fRefPositionAngleDeg=0;
  
  if (fDebugLevel > 0) {
    std::cerr << "************** in makeVSimulationData ********" << std::endl;
    std::cerr << "    VSimulationData Parameters " << std::endl;
    DEBUG(fCORSIKAParticleID); DEBUG(fEnergyGeV);
    DEBUG(fObservationZenithDeg); DEBUG(fObservationAzimuthDeg);
    DEBUG(fPrimaryZenithDeg); DEBUG(fPrimaryAzimuthDeg);
    DEBUG(fRefZenithDeg);  DEBUG(fRefAzimuthDeg);
    DEBUG(fRefPositionAngleDeg);  DEBUG(fCoreEastM);
    DEBUG(fCoreSouthM);  DEBUG(fCoreElevationMASL);
    DEBUG(fFirstInteractionHeight); DEBUG(fFirstInteractionDepth);
    DEBUG(uShowerID); DEBUG(uCorsikaRunNumber);
  }
 
  simu_data = new VSimulationData(fCORSIKAParticleID,fEnergyGeV,
				  fObservationZenithDeg, fObservationAzimuthDeg,
				  fPrimaryZenithDeg, fPrimaryAzimuthDeg, 
                                  fRefZenithDeg, 
				  fRefAzimuthDeg, fRefPositionAngleDeg,
				  fCoreEastM,fCoreSouthM,fCoreElevationMASL);

  packet->putSimulationData(simu_data);


  cors_data = new VCorsikaSimulationData( fFirstInteractionHeight, fFirstInteractionDepth, uCorsikaRunNumber, uShowerID );

  packet->putCorsikaSimulationData( cors_data );


  if (!packet->has(VGetSimulationDataBankName()) ) {
    showXErrorVbfWriter(" packet has no simulation_data_bank as it should ");
  }
  
  return true;
 }

/*************************** makeArrayTrigger *******************************/
bool VG_writeVBF::makeArrayTrigger(const std::vector<int> &localTriggerV) {

  fCMaskTrig = 0;

  if (fDebugLevel > 0) {
    std::cerr << "********* in makeArrayTrigger(localTriggerV) ******" 
              << std::endl;
    std::cerr << "    making VArrayTrigger for data " << std::endl;

    std::cerr << "         trigger vector size: " 
              << localTriggerV.size() << std::endl;
    std::cerr << "         local trigger vector " << std::endl;
    for (unsigned i = 0;i<fnbr_tel;i++) {
      std::cerr << "       " << i << " " << localTriggerV[i] << std::endl;
    }
    std::cerr << std::endl;
  }

  if (localTriggerV.size() != fnbr_tel ) {
    showXErrorVbfWriter("size of localTriggerV != fnbr_tel");
  }

  // make an array trigger when need indiv. telescope triggering
  // information
  if (at == 0) {
    at = new VArrayTrigger();
  }
  else {
    showXErrorVbfWriter("can't make VArrayTrigger since at != 0 ");
  }
    
  // make event types, for non-forced and forced
  VEventType evtype_trig(VEventType::L2_TRIGGER,false,false,
                          VEventType::NOT_CALIBRATION,false);
  ubyte rawcode_trig = evtype_trig.getBestNewStyleCode();

  VEventType evtype_force(VEventType::L2_TRIGGER,true,false,
                          VEventType::NOT_CALIBRATION,false);
  
  ubyte rawcode_force = evtype_force.getBestNewStyleCode();
    
  // set up trig_tel and data_tel vectors
  std::vector<unsigned> trig_tel; 
  std::vector<unsigned> data_tel;
  
  // trig_bool[tel] = true if telescope,tel, triggered
  std::vector<bool> trig_bool(fnbr_tel,false);

  for (unsigned tel = 0;tel<fnbr_tel;tel++) {
    if (localTriggerV[tel] > 0.5) {
      trig_tel.push_back(tel);
      data_tel.push_back(tel);  

      trig_bool[tel] = true;
    }
    else { 
      if (ftriggered_readout) {
        data_tel.push_back(tel);  
      }
    }
  }

  fnumTelsWithData = data_tel.size(); 
  fnumTrigTels     = trig_tel.size();

  if (fDebugLevel > 0) {
    std::cerr << "  array trigger set parameters " << std::endl;
    std::cerr << "    fnumTelsWithData: " << " " 
              << fnumTelsWithData << std::endl;
    std::cerr << "    fnumTrigTels: " 
              << fnumTrigTels << std::endl;
    std::cerr << " fEventNumber: " << fEventNumber 
              << std::endl;
    std::cerr << " fRunNumber: " <<fRunNumber 
              << std::endl;
    std::cerr << "      eventtype: " << (unsigned int)rawcode_force 
              << std::endl;
  }
 
  at->resizeSubarrayTelescopes((uint8_t)fnumTelsWithData);
  at->resizeTriggerTelescopes((uint8_t)fnumTrigTels);
  at->setEventNumber((uint32_t)fEventNumber);
  at->setNodeNumber(255); 
  at->setRunNumber((uint32_t)fRunNumber);
  at->setFlags(0);  
  at->setATFlags(0);
  
  // telescope pointing done in optics.h prior to this call
  
  // set gps time, if pedestal, gps arrays set in make_packet
  at->getGPSTime()[0]=fGPSWords[0];
  at->getGPSTime()[1]=fGPSWords[1];
  at->getGPSTime()[2]=fGPSWords[2];
  at->getGPSTime()[3]=fGPSWords[3];
  at->getGPSTime()[4]=fGPSWords[4];
  
  at->setGPSYear(fGPSYear);
  
  // set rawevent code, always the same for non-ped array trigger
  at->setRawEventTypeCode(rawcode_trig);

  for (unsigned i=0;i < fnumTelsWithData;i++) {
    unsigned tel = data_tel.at(i); // triggered telescope number
    at->setSubarrayTelescopeId(i,(uint32_t)tel);
    
    if (trig_bool[tel]) {
      at->setSpecificRawEventTypeCode(i,rawcode_trig);

      at->setShowerDelay(i,0);
      at->setCompDelay(i,0);
      at->setAltitude(i,eltel[tel]);
      at->setAzimuth(i,aztel[tel]);
      at->setTDCTime(i,0);
    }
    else {
      if (ftriggered_readout) {
        at->setSpecificRawEventTypeCode(i,rawcode_force);
        at->setShowerDelay(i,0);
        at->setCompDelay(i,0);
        at->setAltitude(i,eltel[tel]);
        at->setAzimuth(i,aztel[tel]);
        at->setTDCTime(i,0);
      }
    }
  }
  
  // set trigger telescope ids and trigger mask
  std::ostringstream os;
  for (unsigned i=0;i<fnumTrigTels;i++) {
    unsigned tel = trig_tel.at(i); // triggered telescope number
    at->setTriggerTelescopeId(i,(uint32_t)tel);
    os << tel << ",";
  }
  
  std::string trigger_config_mask (os.str());
  //std::cout << "trigger_config_mask: " << trigger_config_mask << std::endl;
  
  fCMaskTrig =
    toDAQMask(parseConfigMask(trigger_config_mask.c_str()));

  at->setTriggerMask(fCMaskTrig);

  // set configuration mask, always use full array of telescopes
  //uint8_t fCMask=(uint8_t)toDAQMask(parseConfigMask(fConfigMaskAll.c_str()));
  at->setConfigMask(fCMaskAll);

  if (fDebugLevel > 0) {
    std::cerr << "SETTRIGMASK DATA: " << (unsigned)at->getTriggerMask() 
              << std::endl;
    std::cerr << "      trigger_config_mask: " << trigger_config_mask 
              << std::endl;
    std::cerr << "       fCMaskTrig: " << (unsigned)fCMaskTrig 
              << std::endl; 
    std::cerr << "      fCMaskAll:   " << (unsigned)fCMaskAll << std::endl;
    std::cerr << "************* end of makeArrayTrigger()" << std::endl;
  }


  return true;
}

/*************************** makeArrayTrigger *******************************/
bool VG_writeVBF::makeArrayTrigger() {

  // make an array trigger where all telescopes trigger
  // appropriate for making pedestal events

  // set trigger mask to all telescope triggering
  // then VEvent triggermask will be correct.
  fCMaskTrig = fCMaskAll;

  if (fDebugLevel > 0) {
    std::cerr << "******** in makeArrayTrigger() ********** " << std::endl;
    std::cerr << "    making VArrayTrigger for peds " << std::endl;
  }

  if (at == 0) {
    at = new VArrayTrigger();
  }
  else {
    showXErrorVbfWriter("can't make new VArrayTrigger, at !=0");
  }

  std::vector<unsigned> ftrig_tel;   //!< list of triggered telescopes
  std::vector<unsigned> fdata_tel;

  // all telescopes trigger for ped evt
  for (unsigned tel = 0;tel<fnbr_tel;tel++) {  
      ftrig_tel.push_back(tel);
      fdata_tel.push_back(tel);
  }
    fnumTelsWithData = fnbr_tel;  // ok for peds

  if (fDebugLevel > 0) {
    std::cerr << "  array trigger set parameters " << std::endl;
    std::cerr << "    fnumTelsWithData: " << fnumTelsWithData 
              << std::endl;
    std::cerr << " fEventNumber: " << fEventNumber 
              << std::endl;
    std::cerr << " fRunNumber: " <<fRunNumber 
              << std::endl;
  }

    at->resizeSubarrayTelescopes((uint8_t)fnumTelsWithData);
    at->resizeTriggerTelescopes((uint8_t)fnumTelsWithData);
    at->setEventNumber((uint32_t)fEventNumber);
    at->setNodeNumber(255); 
    at->setRunNumber((uint32_t)fRunNumber);
    at->setFlags(0);  // on set ATflags?
    at->setATFlags(0);
   
    // use same config as for simulation header
    //uint8_t fCMask= (uint8_t)toDAQMask(parseConfigMask(fConfigMaskAll.c_str()));
    at->setConfigMask(fCMaskAll);

    // set gps time, if pedestal, gps arrays set in make_packet
    at->getGPSTime()[0]=fGPSWords[0];
    at->getGPSTime()[1]=fGPSWords[1];
    at->getGPSTime()[2]=fGPSWords[2];
    at->getGPSTime()[3]=fGPSWords[3];
    at->getGPSTime()[4]=fGPSWords[4];

    at->setGPSYear(fGPSYear);

    VEventType evtype_ped(VEventType::PED_TRIGGER,true,false,
                          VEventType::NOT_CALIBRATION,false);

    ubyte ped_evtype = evtype_ped.getBestNewStyleCode();
    if (fDebugLevel > 0) {
      std::cerr << "      ped_eventtype (10): " << (unsigned int)ped_evtype 
                << std::endl;

    }

    at->setRawEventTypeCode(ped_evtype);
 
    uint32_t iz = 0;
    for (unsigned tel=0;tel < fnumTelsWithData;tel++) {
      at->setSubarrayTelescopeId(tel,(uint32_t)tel);
      at->setTriggerTelescopeId(tel,(uint32_t)tel);
      at->setSpecificRawEventTypeCode(tel,ped_evtype);
      at->setShowerDelay(tel,iz);
      at->setCompDelay(tel,iz);
      at->setAltitude(tel,eltel[tel]);
      at->setAzimuth(tel,aztel[tel]);
      at->setTDCTime(tel,iz);
    }
    
    at->setTriggerMask(fCMaskAll);
    if (fDebugLevel > 0) {
      std::cerr << "fCMaskAll: " << (unsigned)fCMaskAll << std::endl;
      std::cerr << "PED at->getTriggerMask: " << (unsigned)at->getTriggerMask() 
                << std::endl;
      std::cerr << "************* end of makeArrayTrigger()" << std::endl;
    }
    
  return true;
}

/*************************** storeArrayTrigger *******************************/
  bool VG_writeVBF::storeArrayTrigger() {

    if (fDebugLevel > 0) {
      std::cerr << "*********** in store_arraytrigger() *****" << std::endl;
    }
 
   // put array trigger into arrayevent
    if (ae != 0) {
      ae->setTrigger(at);
    }
    else {
      showXErrorVbfWriter(" can't store arraytrigger, no arrayevent");
    }
   
  return true;
}

/************************ getStats **********************************/
void VG_writeVBF::getStats() {

  std::cerr << "Number of normal events written to vbf file: "
            << fEventNumber - fPedEventCount << std::endl;
  std::cerr << "Number of pedestal events written to vbf file: "
            << fPedEventCount << std::endl;
  
}

/************************ finishVBF **********************************/
void VG_writeVBF::finishVBF() {
  
  // finish up, this creates the index and writes the checksum
  pfWriter->finish();
  
}

/**************************** KascadeType2CorsikaType ******************************/
uint32_t VG_writeVBF::KascadeType2CorsikaType(uint32_t fKType) 

// ************************************************************************
// Convert Kascade primary particle type to Corika particle type
// ************************************************************************
//From CORSIKA manual Table 4 (pg 80 in current manual)
//
//gamma        1
//e+           2
//e-           3
//muon+        5
//muon-        6
//neutron      13
//proton       14
//anti-proton  15
//ION(Z,A)     A x 100 + Z
//He(2,4)      402
//Fe(26,56)    5626
// **************************************************************************
// Kascade Particle species codes.
//         1:Gamma
//         2:positron
//         3:electron
//         4:muon (+)
//         5:muon (-)
//         6:pion (0)
//         7:pion (+)
//         8:pion (-)
//         9:kaon (+)
//        10:kaon (-)
//        11:kaon (0long)
//        12:kaon (0short)
//        13:proton
//        14:neutron
//        15:neutrino(electron)
//        16:anti-neutrino(electron)
//        17:neutrino(muon)
//        18:anti-neutrino(muon)
//        ION(Z,A)=A+20
//        He4=24;
//        Fe56=60;
// ****************************************************************************
{
  switch(fKType)
    {
    case(1):return 1;//gamma
    case(2):return 2;//electron
    case(3):return 3;//positron
    case(4):return 5;//mu+
    case(5):return 6;//mu-
    case(13):return 14;//proton
    case(14):return 13;//neutron
    default: 
      {
	if(fKType<13 || (fKType>14 &&fKType<22))
	  {
	    std::cerr<<"ksAomega: KSEvent:Primary Particle type: "<<fKType
		     <<" not in table."<<std::endl;
	    exit(1);
	  }
	else
	  {
	    int fZNuclei;
	    double fXMass;
	    int fKascadeHeavyType=fKType-20;
	    MassNumber2ChargeAndMass(fKascadeHeavyType,fZNuclei,fXMass);
	    return 100*fKascadeHeavyType+fZNuclei;
	  }
      }
    }
}

// **************************** MassNumber2ChargeAndMass **********************************************
void VG_writeVBF::MassNumber2ChargeAndMass(
                  int fIAtomicNum,int& fCharge, double& fMass)
// **************************************************************************
//   Determine charge Z of most stable nuclei of nuclear number A
//   Determine mass from Z and A using semi-empirical mass formula
// ***************************************************************************
//   Uses fromulae from "Physics of Nuclei and Particles, Vol 1", Marmier and 
//   Sheldon,1969, Acedemic Press, pg 36-38, formula: 2-35,2-38,2-39,2-40
// ***************************************************************************
//   This is from the liquid drop model of stable nuclei.

//  Written by:
//  Glenn Sembroski
//  Physics Dept.
//  Purdue Univ.
//  W. Lafayette, IN USA 47907

// Modified:
//	SRC 3/15/06 Converted from FORTRAN to C++
{
  double fAtomicNum  = (double)fIAtomicNum;
  double fAtomicNumTwoThirds = pow(fAtomicNum,(2./3.));

  const double  fMassProton  = 938.256e-6;  //Mass of proton (TeV)
  const double  fMassNeutron = 939.550e-6;  //Mass of neutron(TeV)
  const double  fPaircoef    = 33.5e-6;     //Pairing coef.(TeV)
  const double  fVolcoef     = 14.1e-6;     //Volume coef (TeV)
  const double  fSurfcoef    = 13.0e-6;     //Surface coef(TeV)
  const double  fCoulcoef    = 0.595e-6;    //Coulmb coef.(TeV)
  const double  fAsymcoef    = 19e-6;       //Assymetry coef(TeV)


  // ***********************************************************************
  //Correct our formula for elements up to a=56(Fe) which is as high as 
  //NUC_LIB goes.
  // ***********************************************************************
  switch(fIAtomicNum) 
    {  
    case 18:fCharge= 8;break;	  //Force Oxygen isotope
    case 24:fCharge=12;break;     //Force Magnesium
    case 28:fCharge=14;break;     //Force silicon
    case 32:fCharge=16;break;     //Force Sulpher
    case 33:fCharge=16;break;     //Force Sulpher
    case 35:fCharge=17;break;     //Force Chlorine
    case 39:fCharge=19;break;     //Force Potassium
    case 40:fCharge=18;           //Force Argon //Could he been calcium 40.
            if(fFirstArgon)
	      { 
		std::cerr<<"Warning--Forcing Argon for all atomic masses of 40"
			 <<std::endl;
	  
		fFirstArgon=false;
	      }
	    break;
    case 56:fCharge=26;break;         //Force Iron.
    default:
      {
	double fChrg=(fAtomicNum/(1.98+(0.0155*fAtomicNumTwoThirds))); 
	fCharge= NearestInt(fChrg);
	//Use nearest integer 
      }
    }
  // *******************************************************************
  // Determine Mass from liquid drop model
  // First determine pairing mass term
  // *******************************************************************
  double fPairingMass;
  if((fCharge % 2)==0)   //even
    {
      if((fIAtomicNum % 2)==0)
	{
	  fPairingMass=-fPaircoef*(pow(fAtomicNum,-.75));   //even-even nuclei
	}
      else
	{
	  fPairingMass=0.;              //even-odd nuclei
	}
    }
  else
    { //charge odd
      if((fIAtomicNum % 2)==0)
	{
	  fPairingMass=0.;	  //Odd-even nuclei
	}
      else
	{
	  fPairingMass=fPaircoef*(pow(fAtomicNum,-.75)); //Odd-odd  nuclei
	}             
    }
  
  fMass = (float)( fCharge*fMassProton 
		   +(fAtomicNum-fCharge)*fMassNeutron 
		   -fVolcoef*fAtomicNum 
		   +fSurfcoef*fAtomicNumTwoThirds  
		   +fCoulcoef*(fCharge*fCharge/pow(fAtomicNum,(1./3.))) 
		   +fAsymcoef*pow((fAtomicNum-2*fCharge),2)/fAtomicNum  
		   +fPairingMass );
  return;
}

// ************************ NearestInt ***************************************************
int VG_writeVBF::NearestInt(double fX)
{
  // Return nearest integer
  int fN=(int)floor((double)fX);
  if((fX-fN)<(fN+1-fX))
    {
      return fN;
    }
  else
    {
      return fN+1;
    }
}

// **************************** getPedEventTime ************************************************
VATime *VG_writeVBF::getPedEventTime() {
  return fPedEventTime;
}

// **************************** setPedEventTime ************************************************
void VG_writeVBF::setPedEventTime(VATime *eventtime) {
  // set pedestal event time = eventtime
  // don't use pointer arithmetic, use VATime equality operator
  
 *fPedEventTime = *eventtime;
  
}

// ****************************** getEventTime **********************************************
VATime *VG_writeVBF::getEventTime() {
  return fEventTime;
}

// **************************** setEventTime ************************************************
void VG_writeVBF::setEventTime(VATime *eventtime) {
  // set event time = eventtime
  // don't use pointer arithmetic, use VATime equality operator
  
  *fEventTime = *eventtime;
}

// ********************************* incrementEventTime *******************************************
bool VG_writeVBF::incrementEventTime() {

  //std::cout << "setting expdev to 0.5 " << std::endl;
  float tmp = 0.8;

  double fEventTimeMJD=fEventTime->getMJDDbl();
  //double increm1 = fMeanTimeBetweenEventsSec*expdev();
  double increm1 = fMeanTimeBetweenEventsSec*tmp;
  //double increm1 = fMeanTimeBetweenEventsSec*frand->Exp(1.0); 
  double increm;

  increm = increm1/(60.*60.*24.); 

  if (fDebugLevel > 0) {
    std::cerr << "************ in incrementEventTime() *******" << std::endl;
    testVATime(fEventTime);
    std::cerr << "   increm: " << increm << std::endl;
  }

  fEventTimeMJD+= increm;
  fEventTime->setFromMJDDbl(fEventTimeMJD);

  if (fDebugLevel > 0) {
    std::cerr << "   test_vatime after increm" << std::endl;
    testVATime(fEventTime);
  }
  
  return true;
}

// ***************************** makePedestalNow ***********************************************
bool VG_writeVBF::makePedestalNow()
{
  // determine if it's time to write a pedestal but only if
  // we need to produce peds

  fmake_pedestal_now_flag = false;
  
  if (fPedFlag) {
  
    double fEventTimeMJD=fEventTime->getMJDDbl();
    double fPedEventTimeMJD=fPedEventTime->getMJDDbl();
    
    if (fEventTimeMJD > fPedEventTimeMJD ) {
      fmake_pedestal_now_flag = true;
    } 

  }
  return fmake_pedestal_now_flag;
}

// ******************************** incrementPedTime ********************************************
bool VG_writeVBF::incrementPedTime() {
  // add the pedestal time increment, with NS = 0
  
  uint32_t fYear,fMonth,fDay,H,M,S,NS;

  // make sure pedtime exceeds current eventtime
  while ((*fPedEventTime) < (*fEventTime)) {
  
    fPedEventTime->getCalendarDate(fYear,fMonth,fDay);
    fPedEventTime->getTime(H,M,S,NS);
    NS = 0;
    S += fPedTimeInterval;
    
    fPedEventTime->setFromCalendarDateAndTime(fYear,fMonth,fDay,
                                              H,M,S,NS);
  }
    
  return true;
}

// **************************** testVATime ************************************************
void VG_writeVBF::testVATime(VATime *vat) {

  uint8_t fGPSYear;
  uint16_t fGPSWords[5];
  uint32_t fYear,fMonth,fDay,H,M,S,NS;
  
  vat->getCalendarDate(fYear,fMonth,fDay);
  vat->getTime(H,M,S,NS);
  vat->getForVBF(fGPSYear,5,fGPSWords);

  std::cerr << std::endl;
  std::cerr << " --------  entering testVATime " << std::endl;
  std::cerr << "    date and time: ";
  std::cerr << "    " << fYear << " " << fMonth << " " << fDay << "  ";
  std::cerr << H << " " << M << " " << S << " " << NS << std::endl;

  std::cerr << "    GPS Year " << (int)fGPSYear << " ";
  for (int i=0;i<5;i++) {
    std::cerr << (int)fGPSWords[i] << " ";
  }
  std::cerr << std::endl;
  std::cerr << " --------  leaving testVATime " << std::endl;  

}
// **********************************************
void VG_writeVBF::setHiLoGain(const bool  &bhilo) {
  // set hilo gain switchf or current pixel and telescope
  event->setHiLo(fcurrent_pix,bhilo);
}

//************************* setTriggerBit *******************
void VG_writeVBF::setTriggerBit(bool  btrigger) {

  // set trigger bit if cvf fired
  event->setTriggerBit(fcurrent_pix,btrigger);
}

//************************ setTriggeredReadout ********************
void VG_writeVBF::setTriggeredReadout(const bool &triggered_readout) {

  ftriggered_readout = triggered_readout;
  
}

//********************* setSeed ***********************
void VG_writeVBF::setSeed(const unsigned &iseed) {

  UInt_t seed = (UInt_t)iseed;
  std::cout << "setting TRandom3 seed in vbf writer: " << seed << std::endl;
  frand->SetSeed(seed);

}

//************************ setDebugLevel ********************
void VG_writeVBF::setDebugLevel(const int &debugLevel) {

  fDebugLevel=debugLevel;
  //std::cout << "fDebugLevel set to: " << fDebugLevel << std::endl;
}

//********************************************
 int VG_writeVBF::showErrorVbfWriter(const char *msg) {
  fprintf(stderr,"%s\n",msg);
  return -1;
}

//*********************** showXErrorVbfWriter *********************
 int VG_writeVBF::showXErrorVbfWriter(const char *msg) {
  fprintf(stderr,"%s\n",msg);
  exit(-1);
}
