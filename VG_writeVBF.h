/*! \brief VG_writeVBF class for writing vbf simulation files
  Version 1.0
  Charlie Duke
  Grinnell College
  May 5, 2010

  This code mirrors the older GrISU_writeVBF class used with the grisudet 
  detector code in GrISU. However, this code has a generic interface
  independent of the grisudet code and can be used with any simulation 
  detector code.
 */

#ifndef VG_WRITEVBF
#define VG_WRITEVBF

class VBankFileWriter; //!< forward declarations of VBF classes and structure
struct VArrayConfiguration;
class VPixelLocation;
class VATime;
struct VEventType;
class  VPacket;
class  VArrayEvent;
class  VEvent;
class  VArrayTrigger;
struct VSimulationData;
class  VCorsikaSimulationData;

class TRandom3;  //!< forward declaration of random number generator

class VG_writeVBF {
 private:

  bool fPedFlag;  //!< if true, pedestals are generated
  std::string fPedFileName; //!< filename for separate peds file, if "" no file
 
  std::vector<  std::vector<VPixelLocation> > *fTelPixV; //!< vector containing pixel location vector for each telescope 

  std::vector<float> fTelSouth; //<! vector containing telescope locations, South direction
  std::vector<float> fTelEast; //<! vector containing telescope locations, East direction
  std::vector<float> fTelUp; //<! vector containing telescope locations, Up direction

  TRandom3 *frand; //<! pointer to random number generator from ROOT

  unsigned fnbr_tel;    //!< number of telescopes
  float fObsAlt;        //!< observatory altitude
  std::vector<unsigned> fnumPixV;  //!< vector of number of pixels in each camera

  std::string fVbfFileName; //!< name of vbf file for writing
  long fRunNumber;      //!< runnumber (from pilot file or default)

  bool fpedestalevent;  //!< true for pedestal event, false for data event
  int fDebugLevel; //!< debug print level; if 0, no debug printing
  
  u_int32_t fEventNumber;    //!< current event number
  u_int32_t fPedEventCount;  //!< current pedestal count 
  
  VATime *fEventTime;           //!< event time for next event
  VATime *fPedEventTime;        //!< event time for next pedestal event
  u_int8_t fGPSYear;             //!< GPS variables used by all VATime classes
  u_int16_t fGPSWords[5];        //!< GPS words used by VATime classes
  
  unsigned fnumFADCSamples;     //!< number of fadc time bins
  unsigned fnumPmts;            //!< number of pixels

  unsigned fnumTelsWithData;    //!< number of telescopes with data
  unsigned fnumTrigTels;        //!< number of triggered telescopes
  bool ftriggered_readout;         //!< if true, readout all telescopes
  
  VEventType *fEvType;  //!< pointer to eventtype created in make_packet
  
  std::vector<float> aztel;    //!< telescope azimuth in degree units
  std::vector<float> eltel;    //!< telescope elevation in degree units

  std::string fConfigMaskPed;  //!< config.mask string for ped events
  std::string fConfigMaskAll;  //!< config.mask string for all telescopes

  unsigned short fCMaskPed;  //<! vbf config.make for pedestal events
  unsigned short fCMaskTrig; //<! vbf config.mask for data, local triggers
  unsigned short fCMaskAll;  //<! vbf config.mask for all telescopes

  unsigned fmaxNumChannels;  //!< set by default to 500 (why not 499?)
  
  unsigned fcurrent_tel;     //!< current telescope id
  unsigned fcurrent_pix;     //!< current pixel id
  unsigned fnum_vevent;      //!< number of VEvent's in shower (for debugging)
  
  bool fmake_pedestal_now_flag; //!< set when event time > ped time
  
  std::string fSimConfigFile;          //!< string written to VSimulationHeader

  std::string fFirstValidEventTimeStr; //!< string to initialize event times
  double fEventRatePerSec;          //!< data rate in Hertz
  double fMeanTimeBetweenEventsSec; //!< mean time between events in seconds
  int fPedTimeInterval;     //!< time interval between pedestal events, in seconds (integer)

  u_int32_t fCorsikaParticleIDType;  //<! Corika style particle type number
  bool  fFirstArgon;       //!< needed for Kascade2Corsika particle type conversion
    
  VBankFileWriter *pfWriter; //!< pointer to vbf bank class
  VPacket *packet;    //!< pointer to vbf packet class
  VArrayEvent *ae;    //!< pointer to vbf arrayevent class
  VEvent *event;      //!< pointer to vbf event class
  std::vector<VEvent*> events;  //!< vector of pointers to VEvent class
  VArrayTrigger *at;  //!< pointer to vbf arraytrigger class
  VSimulationData *simu_data; //!< pointer to vbf simul.data class
  VCorsikaSimulationData * cors_data;

  /*! Used internally in converting kascade id to Corsika id
   */
  void MassNumber2ChargeAndMass(
                   int fIAtomicNum,int& fCharge, double& fMass);

  /*! Used internally to convert double to nearest integer
   */
  int NearestInt(double fX);

  /*! Used internally, for initialization 
   */
  void initialize();

 public:
  
  /*! Constructor
    \param numTel  number of telescopes
    \param numPixTel vector of number of pixels per telescope
    \param vbfFileName  vbf output filename, default "photon.vbf"
    \param vbfRunNumber  vbf run number, default 99999
    \param debugLevel  debug level for testing, default 0 (no printing)
   */
  VG_writeVBF(const int &numTel,
                 const std::vector<unsigned> &numPixTel,
                 const std::string &vbfFileName="photon.vbf",
                 const long &vbfRunNumber=99999,
                 const int &debugLevel=0);

  /*! Destructor
   */
  ~VG_writeVBF() {};

  /*!  set telescope location
    \param telId  telescope identification number, numbering from zero
    \param telSouthLoc  relative location south of array origin
    \param telEastLoc   relative location east of array origin
    \param telUpLoc     relative location up from array origin
   */
  void setTelescopeLocations(const int & telId, const float & telSouthLoc, 
                             const float & telEastLoc,
                             const float & telUpLoc);

  /*!  set pixel locations
    \param telId  telescope identification, numbering from zero
    \param pixId  pixel number for that telescope, numbering from zero
    \param pixEastLocDeg  pixel East location in Degrees, telescope in stow
    \param pixUpLocDeg    pixel Up location in Degrees, telescope in stow
    \param pixRadDeg      pixel radius in Degrees
   */
  void setPixelLocations(const int &telId, const int &pixId,
                         const float &pixEastLocDeg,
                         const float &pixUpLocDeg,
                         const float &pixRadDeg);

  /*! Make vbf VSimulationHeader class and put into vbf packet
    \param configFilePlus  string with user supplied information
    \param simulator  name of simulator
    \param obserAltitudeM observatory altitude in meters above mean sea level
    \param simulationPackage see enum SimulationPackageCodes for possible values
    \param atmosphericModel  use 0, only current value
    \return true
   */
bool makeSimulationHeader(std::string &configFilePlus,
                          const std::string &simulator,
                          const float &obserAltitudeM,
                          const u_int32_t &simulationPackage,
                          const u_int32_t &atmosphericModel=0);

 
/*! set pedestal parameters
  \param makePedsFlag make pedestals if true, default true
  \param pedTimeIntervalS time interval (seconds between peds, default 1 sec
  \param pedFileName name of separate peds file, default "" (no file)
*/
 void setPedestalParameters(const bool &makePedsFlag=true,
                            const int &pedTimeIntervalS= 1,
                            const std::string &pedFileName="");
 
/*! set first event time string
  \param firstEventTime, first event time: example format "2006-11-28 02:00:00.000000000 GPS"
  \return true
*/
 void setFirstEventTimeString(const std::string &firstEventTime);

/*! set shower frequency
  \param showerfreq shower frequency in hz, default 100.0
*/
 void setShowerFrequencyHz(const double &showerfreq= 100.0);

 /*! create VATime classes for showers and peds
   \return true
 */
 bool makeEventTimeClasses();

 /*! set the azimuth(degrees) and elevation(degrees) for each telescope
   \param telId  telescope number, starting from zero
   \param azimDeg  azimuth in degrees
   \param elevDeg  elevation in degrees
   \return true
 */
 void setAzimElevTelDeg(const unsigned &telId, const float &azimDeg,
                        const float &elevDeg);
 
 /*! create new vbf VPacket class
   \return true
 */
 bool makePacket();
 
 /*! store packet in vbf file
   \return true    
 */  
 bool storePacket();
 
 /*! make new vbf VEvent class
   \return true    
 */
 bool makeEvent();
 
 /*! store VEvent class in vbf arrayevent 
   \return true    
 */
 bool storeEvent();
 
 /*! make new vbf VArrayEvent class
   \return true    
 */
 bool makeArrayEvent();
 
 /*! store VarrayEvent class in vbf packet
   \return true    
 */
 bool storeArrayEvent();
 
 /*! store time bin using event->SetSample method
   \param ti time bin index
   \param iDC time bin data in digital counts
   \param iPC time bin pedestal in digital counts    
 */
 void storeSample(int &ti,int &iDC);
 
 /*! make a VSimulationData class using default values
   and put into the vbf packet. Normally used with pedestal packets
   \return true    
 */
 // NEED TO FIX THIS, USE DEFAULTS BELOW IF POSSIBLE.
 bool makeVSimulationData1();
 
 /*! make a VSimulationData class 
     use no arguments in call for simulation data for pedestals
   the vbf packet. Normally used with data packets
   \param energyGeV primary energy in GeV
   \param primaryZenithAngleDeg primary zenith angle in degrees
   \param primaryAzimuthDe primary azimuth in degrees
   \param primaryCoreEastM core location East
   \param primaryCoreSouthM core location South
   \param primaryCoreElevMASL, negative value replaces observatory height
                                zero value sets Core Elevation to zero
   \param corsikaParticleID

   \return true    
 */
 bool makeVSimulationData(const float &energyGeV=0.0,
                          const float &primaryZenithAngleDeg=0.0,
                          const float &primaryAzimuthDeg=0.0,
                          const float &primaryCoreEastM=0.0,
                          const float &primaryCoreSouthM=0.0,
                          const float &primaryCoreElevMASL=0.0,
                          const float &firstInteractionHeight=0.0, 
                          const float &firstInteractionDepth=0.0, 
                          const u_int32_t &corsikaParticleID=1,
                          const u_int32_t &corsikaRunNumber = 0,
                          const u_int32_t &fShowerID = 0
);
  
 /*! make a VArrayTrigger class, get triggering information 
   \param localTriggerV  vector of local trigger 0/1 for each telescope
   \return true    
 */
 bool makeArrayTrigger(const std::vector<int> &localTriggerV);
 
 /*! make a VArrayTrigger, used for pedestals only 
   \return true    
 */
 bool makeArrayTrigger();
 
 /*! store the VArrayTrigger class in VArrayEvent
   \return true    
 */
 bool storeArrayTrigger();
 
 /*! create vbf file index and write checksum    
 */
 void finishVBF();
 
 /*! labels event as a pedestal event
   \return true    
 */
 void setPedestalEvent() {fpedestalevent=true;};
 
 /*! labels event as a data event
  */
 void setDataEvent() {fpedestalevent=false;};
 
 /*! sets the current telescope number (starting at 0)
   \param tel current telescope number
 */
 void setTelescope(const int &tel) {fcurrent_tel = tel;};
 
 /*! sets the current pixel number (starting at 0)
   \param pix current pixel number
 */
 void setPixel(const int &pix) {fcurrent_pix = pix;};
 
 /*! gets the current vbf event number
   \return eventnumber
 */
 u_int32_t getEventNumber() {return fEventNumber;};
 
 /*!  gets the current vbf pedestal event count
   \return current pedestal event count
 */
 u_int32_t getPedEventCount() {return fPedEventCount; };
 
 /*! prints to cout: event number and pedestal event number
  */
 void getStats();
 
 /*! converts kascade particle id to Corsika particle id
   \param fKType kascade particle id
   \return Corsika particle id
 */
 u_int32_t KascadeType2CorsikaType(u_int32_t fKType);
 
 // time management methods (can connect two vbf writers, e.g. for peds)
 
 /*! gets current pedestal event time as pointer to VATime class
   \return VATime class pointer for current pedestal event time
 */
 VATime *getPedEventTime();
 
 /*! sets pedestal event time
   \param eventtime pointer to VATime class
   \return true    
 */
 void setPedEventTime(VATime *eventtime);
 
 /*! sets eventtime
   \param eventtime pointer to VAtime class   
 */
 void setEventTime(VATime *eventtime);
 
 /*! gets current event time as pointer to VATime class
   \return true    
 */
 VATime *getEventTime();
 
 /*! increment event time 
  */
 bool incrementEventTime();
 
 /*! increment pedestal time
   \return true    
 */
 bool incrementPedTime();
 
 /*! flag that, if true, indicates it's time to make a pedestal event
   \return bool    
 */
 bool makePedestalNow();

 /*! print year, month, day, hour, min, sec, nanosec for testing
   \param vat VATime class pointer
 */
 void testVATime(VATime *vat);
 
 /*! set change and pedestal for the current pixel to zero
  */
 void setChargePedHigain();
 
 /*! set high/low gain
   \param bhilo if true, set to high gain; false, low gain
  */ 
 void setHiLoGain(const bool &bhilo);
 
 /*! set pixel trigger bit, need to setPixel(int pix) first
   /param btrigger set to true if pixel fired
  */ 
 void setTriggerBit(bool  btrigger);
 
 /*!
  */ 
 void setTriggeredReadout(const bool &triggered_readout);
 
 /*! set seed for TRandom3, unsigned integer
   /param iseed, seed: default = 111
  */ 
 void setSeed(const unsigned &iseed=111);
 
 /*! set particleIDType
   \param particleIDType, 0 for kascade, 1 for corsika
 */
 void setParticleIDType(const u_int32_t &particleIDType) {
   fCorsikaParticleIDType = particleIDType;
 };

 /*! set debug print level, 0: no printout, 1: print all parameters
   \param debugLevel  
 */
 void setDebugLevel(const int &debugLevel);

 /*! setMaxNumberChannels: internally set to 500, set only if !=500. same for all cameras
   \param maxNumberChannels
 */
 void setMaxNumberChannels(const unsigned &maxNumberChannels) {
   fmaxNumChannels = maxNumberChannels;
 };

 /*! setNumFadcSamples: internally set to 24, set only if !=24, all fadcs have same number of channels
   \param numFadcSamples
 */
 void setNumFadcSamples(const unsigned &numFadcSamples) {
   fnumFADCSamples = numFadcSamples;
 };

 /*!  print error message to stderr, do not exit
  */
int showErrorVbfWriter(const char *msg);

 /*! print error message to stderr, exit code
  */
int showXErrorVbfWriter(const char *msg);

};

#endif
