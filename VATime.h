//-*-mode:c++; mode:font-lock;-*-

/**
 * \file VATime.h
 *
 * Time conversion/manipulation class for VEGAS
 * 
 * Original Author: Stephen Fegan
 * $Author: nepomuk $
 * $Date: 2010/11/13 01:32:22 $
 * $Revision: 1.2 $
 * $Tag$
 *
 **/

#ifndef VATIME_H
#define VATIME_H

#ifdef __APPLE__
#ifdef __CINT__
#undef __GNUC__
#endif
#endif

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstdio>
#include <sys/time.h>
#include <stdint.h>

#ifndef NOROOT
#include <TObject.h>
#endif

#ifdef __APPLE__
#if  __GNUC__ < 4
//#warning Fix for POSIX non-complient Mac/BSD should be removed when no longer needed
typedef long suseconds_t;
#endif
#endif

class VATime
#ifndef NOROOT
  : public TObject
#endif
{
public:
  enum TimeSystem { TS_TAI, TS_UTC, TS_GPS };
  
  VATime(TimeSystem ts=TS_UTC): 
    fDecodedOK(false), fGPSStatus(), fTimeSystem(ts), fMJD(), fDayNS()
  { }

  VATime(const std::string& time):
    fDecodedOK(false), fGPSStatus(), fTimeSystem(TS_UTC), fMJD(), fDayNS()
  { setFromString(time) ; }

  inline VATime(uint8_t gps_year, 
		unsigned num_gps_words, const uint16_t* gps_words);
  virtual ~VATime();
  
  // -------
  // SETTERS
  // -------
  inline bool setFromVBF(uint8_t gps_year, 
			 unsigned num_gps_words, const uint16_t* gps_words);
  inline bool setFromMJDIntAndNS(uint32_t mjd, uint64_t day_ns);
  inline bool setFromPOSIXTimeT(const time_t* timep, uint32_t ns=0);
  inline bool setFromPOSIXTimeVal(const struct timeval* tv);
  inline bool setFromSystemTime();
  inline bool setFromMJDDbl(double mjd);
  inline bool setFromCalendarDateAndTime(uint32_t Y, uint32_t M, uint32_t D,
					 uint32_t HRS, uint32_t MIN, 
					 uint32_t SEC, uint32_t NS);
  inline bool setFromCalendarDOYAndTime(uint32_t Y, uint32_t DOY,
					uint32_t HRS, uint32_t MIN, 
					uint32_t SEC, uint32_t NS);
  inline bool setFromMSTimeStamp(uint64_t timestamp);
  inline bool setFromDBTimeStamp(uint64_t timestamp);
  inline bool setFromString(const std::string& str);
  
  inline void setTimeSystemWithoutConvertingTime(TimeSystem s);
  
  // ---------------
  // TIME CONVERSION
  // ---------------
  inline bool convertTimeToSystem(TimeSystem s);

  // -------
  // GETTERS
  // -------
  inline bool getForVBF(uint8_t& gps_year, 
			unsigned num_gps_words, uint16_t* gps_words) const;

  //! Return the integer component of the modified Julian date
  inline uint32_t getMJDInt() const { return fMJD; }

  //! Return the number of nanoseconds since midnight as an unsigned
  //  64-bit integer
  inline uint64_t getDayNS() const { return fDayNS; }
  
  inline double getMJDDbl() const;
  
  inline time_t getPOSIXTimeT(time_t* t=0) const;
  inline void getPOSIXTimeVal(struct timeval* tv) const;
  
  inline void getCalendarDate(uint32_t& y, uint32_t& m, uint32_t& d) const;
  inline void getTime(uint32_t& h, uint32_t& m, uint32_t& s, 
		      uint32_t& ns) const;
  
  inline uint32_t getYear() const;
  inline uint32_t getMonth() const;
  inline uint32_t getDay() const;
  inline uint32_t getDayOfYear() const;
  inline uint32_t getDayOfWeek() const;
  inline uint32_t getHour() const;
  inline uint32_t getMin() const;
  inline uint32_t getSec() const;
  inline uint32_t getNanoSec() const;

  inline uint64_t getMSTimeStamp() const;
  inline uint64_t getDBTimeStamp() const;
  
  inline void getString(std::string& str) const;
  //! Convenience function which calls getString(std::string&)
  inline std::string getString() const 
  { std::string s; getString(s); return s;}
  //! Convenience function which calls getString(std::string&)
  inline std::string toString() const { return getString();}

  inline TimeSystem getTimeSystem() const;
  
  //! Return true if the GPS clock indicated it was fly wheeling (unlocked
  // from GPS time source).
  bool isFlyWheeling() const { return fGPSStatus&0x0001; }

  //! Return true if the GPS clock indicated there was likely a large
  // offset in the time it reports with respect to UTC.
  bool isLargeTimeOffset() const { return fGPSStatus&0x0002; }

  //! See Symmetricon GPS specifications
  bool isLargeFreqOffset() const { return fGPSStatus&0x0004; }

  //! See Symmetricon GPS specifications
  bool isBatteryFailed() const { return fGPSStatus&0x0008; }

  // -----------------
  // LOGICAL OPERATORS
  // -----------------
  inline bool isGood() const;
  inline operator bool() const;
  
  inline bool operator <= (const VATime& o) const;
  inline bool operator == (const VATime& o) const;
  inline bool operator >= (const VATime& o) const;
  inline bool operator < (const VATime& o) const;
  inline bool operator > (const VATime& o) const;
  inline bool operator != (const VATime& o) const;
  
  // ----------------------
  // MATHEMATICAL OPERATORS
  // ----------------------
  inline VATime& operator += (const int64_t& ns);
  inline VATime& operator -= (const int64_t& ns);
  inline int64_t operator - (const VATime& o) const;
    
  // ----------------
  // STATIC FUNCTIONS
  // ----------------
  static inline VATime now();
  static inline bool isLeapYear(uint32_t Y);
  static inline uint32_t calendarDateToDOY(uint32_t Y, uint32_t M, uint32_t D);

  static bool getIgnoreLeapSeconds() { return sIgnoreLeapSeconds; }
  static void setIgnoreLeapSeconds(bool i=true) { sIgnoreLeapSeconds=i; }

  // ---------------------
  // TIME SYSTEM FUNCTIONS
  // ---------------------
  static inline uint64_t getDayLength(uint32_t mjd, TimeSystem ts=TS_UTC);
  static inline void setZeroOffsetToTAI(int offset, const std::string& name,
					  TimeSystem ts=TS_UTC);
  static inline void setOffsetToTAI(uint32_t mjd, int offset,
				    TimeSystem ts=TS_UTC);
  static inline void clearAllOffsets(TimeSystem ts=TS_UTC);
  static inline void dumpTAIOffsets(std::ostream& stream,
				    TimeSystem ts=TS_UTC);
  static inline int getOffsetToTAI(uint32_t mjd, TimeSystem ts=TS_UTC);
  
  static void configure(/* VSOptions* options, ConfigInfo* info */);
  
private:
  bool       fDecodedOK;
  uint32_t   fGPSStatus;
  TimeSystem fTimeSystem;
  uint32_t   fMJD;
  uint64_t   fDayNS;
  
  struct TimeSystemDefinition
  {
    TimeSystem              fTimeSystem;
    std::string             fTimeSystemName;
    int                     fZeroOffsetToTAI;
    std::map<uint32_t, int> fOffsetToTAI;
  };

  static bool                              sIgnoreLeapSeconds; //!
  static std::vector<TimeSystemDefinition> sTimeSystems; //!

  static const uint64_t JULIAN_NS_U = 86400000000000ULL;
  static const int64_t  JULIAN_NS_I = 86400000000000LL;
  static const uint64_t HOUR_NS_U   = 3600000000000ULL;
  static const int64_t  HOUR_NS_I   = 3600000000000LL;
  static const uint64_t MIN_NS_U    = 60000000000ULL;
  static const int64_t  NIN_NS_I    = 60000000000LL;
  static const uint64_t SECOND_NS_U = 1000000000ULL;
  static const int64_t  SECOND_NS_I = 1000000000LL;
  
  static const uint32_t JULIAN_S_U  = 86400;
  static const int32_t  JULIAN_S_I  = 86400;    

#ifndef NOROOT  
  ClassDef(VATime,1); // VATime
#endif
};

// ============================================================================
// Inline static functions used below
// ============================================================================

/**
 * Static function which returns "true" if the year is a leap year.
 **/ 
inline bool VATime::isLeapYear(uint32_t Y)
{
  return (Y%4==0)&&((Y%100!=0)||(Y%400==0));
}

/**
 * Static function to calculate the day of year given the calendar 
 * year, month and date.
 **/ 
inline uint32_t VATime::calendarDateToDOY(uint32_t Y, uint32_t M, uint32_t D)
{
  // Julian Date to Day Of Year Conversion Algorithm
  // Reference: John Meeus, Astronomical Algorithms, Second Edition, 
  // Willmann-Bell, 1998, Chapter 7

  uint32_t K   = isLeapYear(Y)?1:2;
  uint32_t DOY = 
    uint32_t(floor(double(275*M)/9.0)) - 
    K * uint32_t(floor(double(M+9)/12.0)) + D - 30;
  return DOY;
}

/**
 * Static function which returns the length of the day (in nanoseconds)
 * given the modified Julian date and the time system.
 * year, month and date.
 **/ 
inline uint64_t VATime::getDayLength(uint32_t mjd, TimeSystem ts)
{
  uint32_t sec_len = JULIAN_S_U;
  if((!sIgnoreLeapSeconds)&&(int(ts)<(int)sTimeSystems.size()))
    {
      std::map<uint32_t, int>::const_iterator i = 
	sTimeSystems[int(ts)].fOffsetToTAI.find(mjd+1);
      if(i != sTimeSystems[int(ts)].fOffsetToTAI.end())
	{
	  if(i==sTimeSystems[int(ts)].fOffsetToTAI.begin())
	    sec_len += i->second-sTimeSystems[int(ts)].fZeroOffsetToTAI;
	  else
	    sec_len += i->second-(--i)->second;
	}
    }
  return uint64_t(sec_len)*SECOND_NS_U;
}

/**
 * Static function which is not yet tested or documented
 **/ 
inline void VATime::setZeroOffsetToTAI(int offset, 
                                       const std::string& name, 
                                       TimeSystem ts)
{
  if(int(ts) >= (int)sTimeSystems.size())sTimeSystems.resize(int(ts)+1);
  sTimeSystems[int(ts)].fTimeSystem = ts;
  sTimeSystems[int(ts)].fTimeSystemName = name;
  sTimeSystems[int(ts)].fZeroOffsetToTAI = offset;
}

/**
 * Static function which is not yet tested or documented
 **/ 
inline void VATime::setOffsetToTAI(uint32_t mjd, 
                                   int offset, 
                                   TimeSystem ts)
{
  if(int(ts) >= (int)sTimeSystems.size())sTimeSystems.resize(int(ts)+1);
  sTimeSystems[int(ts)].fOffsetToTAI[mjd] = offset;  
}

/**
 * Static function which is not yet tested or documented
 **/ 
inline void VATime::clearAllOffsets(TimeSystem ts)
{
  if(int(ts) < (int)sTimeSystems.size())
    {
      sTimeSystems[int(ts)].fZeroOffsetToTAI=0;
      sTimeSystems[int(ts)].fOffsetToTAI.clear();
    }
}

/**
 * Static function which is not yet tested or documented
 **/ 
inline void VATime::dumpTAIOffsets(std::ostream& stream, TimeSystem ts)
{
  if(int(ts) < (int)sTimeSystems.size())sTimeSystems.resize(int(ts)+1);
    {
      stream << "zero  " << sTimeSystems[int(ts)].fZeroOffsetToTAI 
	     << std::endl;
      std::map<uint32_t, int>::const_iterator i = 
	sTimeSystems[int(ts)].fOffsetToTAI.begin();
      while(i != sTimeSystems[int(ts)].fOffsetToTAI.end())
	{
	  stream << i->first << ' ' << i->second << std::endl;
	  i++;
	}
    }
}

/**
 * Static function which is not yet tested or documented
 **/ 
inline int VATime::getOffsetToTAI(uint32_t mjd, TimeSystem ts)
{
  int offset=0;
  if(int(ts) < (int)sTimeSystems.size())
    {
      offset=sTimeSystems[int(ts)].fZeroOffsetToTAI;
      std::map<uint32_t, int>::const_iterator i = 
	sTimeSystems[int(ts)].fOffsetToTAI.upper_bound(mjd);
      if(i!=sTimeSystems[int(ts)].fOffsetToTAI.begin())
	offset = (--i)->second;
    }
  return offset;
}

// ============================================================================
// Setters
// ============================================================================

/**
 * Set the VATime given the calendar date (Y, M, D) and time (HRS,
 * MIN, SEC, NS)
 **/ 
inline bool VATime::setFromCalendarDateAndTime(uint32_t Y, 
                                               uint32_t M, 
                                               uint32_t D,
                                               uint32_t HRS, 
                                               uint32_t MIN, 
                                               uint32_t SEC, 
                                               uint32_t NS)
{
  fDecodedOK = true;
  fGPSStatus = 0;
  
  // Calendar Date to Julian Date Conversion Algorithms
  // Reference: John Meeus, Astronomical Algorithms, Second Edition, 
  // Willmann-Bell, 1998, Chapter 7

  if(M<=2) { Y--; M+=12; }
  uint32_t A   = Y/100;       // integer division
  uint32_t B   = 2 - A + A/4; // B is negative -- seems to work as uint32_t :-)

//warning CLEAN ME UP
#if 0
  fMJD        = 
    uint32_t(floor(365.25*double(Y+4716))) + 
    uint32_t(floor(30.6001*double(M+1))) + D + B -
    uint32_t(floor(1524.5+2400000.5));
  std::cerr << "MJD 1: " << fMJD << std::endl;

  fMJD          = 
    uint32_t(floor(365.25*double(Y+4716))) + 
    uint32_t(floor(30.6001*double(M+1))) + D + B - 2401525;
  std::cerr << "MJD 2: " << fMJD << std::endl;
#endif

  // From Meeus:
  //
  // JD = INT(365.25(Y+4716))+INT(30.6001(M+1))+D+B-1524.5
  //
  // Therefore:
  //
  // MJD = INT(365.25(Y+4716))+INT(30.6001(M+1))+D+B-1524.5-2400000.5
  //     = INT(365.25(Y+4716))+INT(30.6001(M+1))+D+B-2401525
  //     = INT(365.25(Y-1856))+INT(30.6001(M+1))+D+B-1102

  fMJD         = 
    uint32_t(floor(365.25*double(Y-1856))) + 
    uint32_t(floor(30.6001*double(M+1))) + D + B - 1102;

#if 0
  std::cerr << "MJD 3: " << fMJD << std::endl;
#endif

  fDayNS       = HRS*HOUR_NS_U + MIN*MIN_NS_U + SEC*SECOND_NS_U + NS;

  return fDecodedOK;
}

/**
 * Set the VATime given the calendar year (Y), day of year (DOY) and
 * time (HRS, MIN, SEC, NS)
 **/ 
inline bool VATime::setFromCalendarDOYAndTime(uint32_t Y, 
                                              uint32_t DOY,
                                              uint32_t HRS, 
                                              uint32_t MIN, 
                                              uint32_t SEC, 
                                              uint32_t NS)
{
  // Day Of Year to Calendar Date Algorithm
  // Reference: John Meeus, Astronomical Algorithms, Second Edition, 
  // Willmann-Bell, 1998, Chapter 7

  uint32_t K   = isLeapYear(Y)?1:2;
  uint32_t M   = (DOY<32)?1:uint32_t(floor(double(9*(K+DOY))/275.0 + 0.98));
  uint32_t D   = (DOY - uint32_t(floor(double(275*M)/9.0)) 
		  + K*uint32_t(double(M+9)/12.0) + 30);
  
#if 0
  std::cerr << "DOY:" << DOY 
	    << " Y:" << Y 
	    << " M:" << M 
	    << " D:" << D 
	    << " HRS:" << HRS 
	    << " MIN:" << MIN
	    << " SEC:" << SEC
	    << " NS:" << NS
	    << std::endl;
#endif

  return setFromCalendarDateAndTime(Y,M,D,HRS,MIN,SEC,NS);
}

/**
 * Set the VATime given the number of years since 2000 (gps_year) and
 * the Symmetricom time words (gps_words), which presumably comes from
 * VBF. The "num_gps_words" parameter describes the number of words in
 * the array and MUST be set to a value of 4.
 **/ 
inline bool VATime::setFromVBF(uint8_t gps_year, 
                               unsigned num_gps_words, 
                               const uint16_t* gps_words)
{
  if(num_gps_words!=5){ fDecodedOK=false; return false; }
  else fDecodedOK=true;

  // ============== Adapted from GPSDecoder (in turn from Scott) ==============
  uint32_t STAT = (gps_words[0]>>4 & 0x000F);

  uint32_t Y    = 2000+gps_year;
  //if(STAT)Y-=100; // Flag to make user notice GPS time is wrong
  
  uint32_t DOY  = 
    (gps_words[0]>>0  & 0x000F) * 100 + 
    (gps_words[1]>>12 & 0x000F) * 10 +
    (gps_words[1]>>8  & 0x000F);

  uint64_t HRS  = (gps_words[1]>>4  & 0x000F) *10 + (gps_words[1]>>0 & 0x000F);
  uint64_t MIN  = (gps_words[2]>>12 & 0x000F) *10 + (gps_words[2]>>8 & 0x000F);
  uint64_t SEC  = (gps_words[2]>>4  & 0x000F) *10 + (gps_words[2]>>0 & 0x000F);
  
  uint64_t NS   = 
    (gps_words[3]>>12 & 0x000F) * 100000000ULL +
    (gps_words[3]>>8  & 0x000F) * 10000000ULL +
    (gps_words[3]>>4  & 0x000F) * 1000000ULL +
    (gps_words[3]>>0  & 0x000F) * 100000ULL +
    (gps_words[4]>>12 & 0x000F) * 10000ULL +
    (gps_words[4]>>8  & 0x000F) * 1000ULL +
    (gps_words[4]>>4  & 0x000F) * 100ULL;
  // ==========================================================================

  bool retval  = setFromCalendarDOYAndTime(Y, DOY, HRS, MIN, SEC, NS);
  fGPSStatus   = STAT;
  return retval;
}

/**
 * Set the VATime given an integer modified Julian date (mjd) and number
 * of nanoseconds since midnight (day_ns).
 **/ 
inline bool VATime::setFromMJDIntAndNS(uint32_t mjd, uint64_t day_ns)
{
  fDecodedOK = true;
  fGPSStatus = 0;
  fMJD       = mjd;
  fDayNS     = day_ns;
  return fDecodedOK;
}

/**
 * Set the VATime given a POSIX (UNIX) time pointer (timep) and the
 * number of nanoseconds since the start of the second (ns).
 **/ 
inline bool VATime::setFromPOSIXTimeT(const time_t* timep, uint32_t ns)
{
  // See note in getPOSIXTimeT() below
  if(!timep)
    {
      fDecodedOK = false;
    }
  else
    {
      fDecodedOK = true;
      fGPSStatus = 0;
      fMJD       = (*timep/JULIAN_S_I)+40587;
      fDayNS     = uint64_t(*timep%JULIAN_S_I)*SECOND_NS_U+uint64_t(ns);
    }
  return fDecodedOK;
}

/**
 * Set the VATime given a pointer to a POSIX (UNIX) time structure
 * (tv).
 **/ 
inline bool VATime::setFromPOSIXTimeVal(const struct timeval* tv)
{
  if(!tv)
    {
      fDecodedOK = false;
      return false;
    }

#ifdef __APPLE__
  //#warning Fix for POSIX non-complient Mac/BSD should be removed when no longer needed
  time_t tv_sec = tv->tv_sec;
  return setFromPOSIXTimeT(&tv_sec, uint32_t(tv->tv_usec)*1000);
#else
  return setFromPOSIXTimeT(&tv->tv_sec, uint32_t(tv->tv_usec)*1000);
#endif
}

/**
 * Set the VATime directly from the computers system clock.
 **/ 
inline bool VATime::setFromSystemTime()
{
  struct timeval tv;
  if(gettimeofday(&tv,0)>=0)setFromPOSIXTimeVal(&tv);
  else fDecodedOK=false;
  return fDecodedOK;
}

/**
 * Set the VATime from a double precision modified Julian date.
 **/ 
inline bool VATime::setFromMJDDbl(double mjd)
{
  uint64_t day_length = getDayLength(fMJD);
  fDecodedOK = true;
  fGPSStatus = 0;
  fMJD       = uint32_t(floor(mjd));
  fDayNS     = uint64_t(floor(fmod(mjd,1.0)*double(day_length)));
  return fDecodedOK;
}

/**
 * Set the VATime from an integer time stamp with millisecond
 * accuracy, in the form YYYYMMDDhhmmssfff (timestamp), where the
 * final three digits give the fraction of seconds.
 **/ 
inline bool VATime::setFromMSTimeStamp(uint64_t timestamp)
{
  if((timestamp < 10000000000000000ULL)||(timestamp >= 100000000000000000ULL))
    return false;
  uint32_t Y   = timestamp/10000000000000ULL;
  uint32_t M   = timestamp/  100000000000ULL % 100ULL;
  uint32_t D   = timestamp/    1000000000ULL % 100ULL;
  uint32_t HRS = timestamp/      10000000ULL % 100ULL;
  uint32_t MIN = timestamp/        100000ULL % 100ULL;
  uint32_t SEC = timestamp/          1000ULL % 100ULL;
  uint32_t NS  = timestamp                   % 1000ULL * 1000000;
  return setFromCalendarDateAndTime(Y,M,D,HRS,MIN,SEC,NS);
}

/**
 * Set the VATime from an integer time stamp with second accuracy, in
 * the form YYYYMMDDhhmmss (timestamp). This form of timestamp is used
 * by the MySQL database system to represent times.
 **/ 
inline bool VATime::setFromDBTimeStamp(uint64_t timestamp)
{
  if((timestamp < 10000000000000ULL)||(timestamp >= 100000000000000ULL))
    return false;
  uint32_t Y   = timestamp/10000000000ULL;
  uint32_t M   = timestamp/  100000000ULL % 100ULL;
  uint32_t D   = timestamp/    1000000ULL % 100ULL;
  uint32_t HRS = timestamp/      10000ULL % 100ULL;
  uint32_t MIN = timestamp/        100ULL % 100ULL;
  uint32_t SEC = timestamp                % 100ULL;
  return setFromCalendarDateAndTime(Y,M,D,HRS,MIN,SEC,0);
}

/**
 * Set the VATime from a stringified time stamp in on the the
 * following forms:
 * 
 * "YYYY-MM-DD"
 * "YYYY-MM-DD hh:mm:ss"
 * "YYYY-MM-DD hh:mm:ss.fffffffff"
 **/ 
inline bool VATime::setFromString(const std::string& str)
{
  uint32_t Y(0),M(0),D(0),HRS(0),MIN(0),SEC(0),NS(0);
  int num_scan = sscanf(str.c_str(),"%4u-%2u-%2u %2u:%2u:%2u.%9u",
			&Y,&M,&D,&HRS,&MIN,&SEC,&NS);
  if((num_scan != 3)&&(num_scan != 6)&&(num_scan != 7))return false;

  if(num_scan < 7)NS=0;
  if(num_scan < 6)HRS=MIN=SEC=0;

  for(std::vector<TimeSystemDefinition>::const_iterator i=sTimeSystems.begin();
      i!=sTimeSystems.end(); i++)
    {
      if(str.substr(str.length()-i->fTimeSystemName.length())==
	 i->fTimeSystemName)fTimeSystem=i->fTimeSystem;
    }
  return setFromCalendarDateAndTime(Y,M,D,HRS,MIN,SEC,NS);
}

/**
 * Function which is not yet tested or documented
 **/ 
inline void VATime::setTimeSystemWithoutConvertingTime(TimeSystem s)
{
  fTimeSystem=s;
}

// ============================================================================
// Time Conversion
// ============================================================================

/**
 * Function which is not yet tested or documented
 **/ 
inline bool VATime::convertTimeToSystem(TimeSystem s)
{  
  int64_t offset = 
    int64_t(getOffsetToTAI(fMJD, fTimeSystem)-getOffsetToTAI(fMJD, s))*
    SECOND_NS_I;
  
  if(offset > 0){
     uint64_t day_len = getDayLength(fMJD, s);
     if(fDayNS+offset < day_len){
        fDayNS+=offset;
     }else{
        offset -= day_len-fDayNS;
        ++fMJD;
        fDayNS=offset;
     }
    }else if(offset < 0){
      offset = - offset;
      if(fDayNS >= (uint64_t)offset){
         fDayNS-=offset;
      }else{
         offset -= fDayNS;
         --fMJD;
         uint64_t day_len = getDayLength(fMJD, s);
         fDayNS=day_len-offset;
      }
    }

  fTimeSystem=s;
  return true;
}

// ==============
//  Constructors
// ==============

/**
 * Constructor to make time from VBF time components. See the setFromVBF()
 * function.
 **/ 
inline VATime::VATime(uint8_t gps_year, 
                      unsigned num_gps_words, 
                      const uint16_t* gps_words)
             :fDecodedOK(), 
              fGPSStatus(), 
              fTimeSystem(TS_UTC), 
              fMJD(), 
              fDayNS()
{
//warning What TimeSystem does the VERITAS GPS clock use ?
  setFromVBF(gps_year, num_gps_words, gps_words);
}

// ============================================================================
// Getter Functions
// ============================================================================

/**
 * Convert VATime to VBF time components. Inverse of the setFromVBF()
 * function.
 **/ 
inline bool VATime::
getForVBF(uint8_t& gps_year, unsigned num_gps_words, uint16_t* gps_words) const
{
  if(num_gps_words!=5)return false;

  uint32_t Y,M,D;
  getCalendarDate(Y,M,D);

  uint32_t DOY = calendarDateToDOY(Y,M,D);

  uint32_t h,m,s,ns;
  getTime(h,m,s,ns);

  // ==========================================================================

  gps_year = uint8_t(Y-2000);

  gps_words[0] = 
    (fGPSStatus     & 0x000F)<<4   |
    ((DOY/100)      % 10)<<0;
  
  gps_words[1] =
    ((DOY/10)       % 10)<<12 |
    (DOY            % 10)<<8  |
    ((h/10)         % 10)<<4  |
    (h              % 10)<<0;
  
  gps_words[2] =
    ((m/10)         % 10)<<12 |
    (m              % 10)<<8  |
    ((s/10)         % 10)<<4  |
    (s              % 10)<<0;

  gps_words[3] =
    ((ns/100000000) % 10)<<12 |
    ((ns/10000000)  % 10)<<8  |
    ((ns/1000000)   % 10)<<4  |
    ((ns/100000)    % 10)<<0;

  gps_words[4] =
    ((ns/10000)     % 10)<<12 |
    ((ns/1000)      % 10)<<8  |
    ((ns/100)       % 10)<<4;

  // ==========================================================================

  return true;
}

/**
 * Convert VATime to a double precision number giving the modified
 * Julian date.
 **/ 
inline double VATime::getMJDDbl() const
{
  uint64_t day_length = getDayLength(fMJD);
  double mjd = double(fMJD)+double(fDayNS)/double(day_length);
  return mjd;
}

/**
 * Convert VATime to a POSIX (UNIX) time, with second accuracy.
 **/ 
inline time_t VATime::getPOSIXTimeT(time_t* t) const
{
  // From the notes for "man 2 time" on Linux:
  //
  // "POSIX.1 defines seconds since the Epoch as a value to be
  // interpreted as the number of seconds between a specified time and
  // the Epoch, according to a formula for conversion from UTC
  // equivalent to conversion on the naïve basis that leap seconds are
  // ignored and all years divisible by 4 are leap years.  This value
  // is not the same as the actual number of seconds between the time
  // and the Epoch, because of leap seconds and because clocks are not
  // required to be synchronised to a standard reference.  The
  // intention is that the interpretation of seconds since the Epoch
  // values be consistent; see POSIX.1 Annex B 2.2.2 for further
  // rationale."
  //
  // The POSIX epoch is defined as: Thu Jan  1 00:00:00 1970 = MJD 40587
  //
  // Since every leap year between 1901 and 2099 is divisible by four
  // the following algorithm is good in that range

  time_t my_time;
  if(t==0)t=&my_time;
  *t = (time_t(fMJD)-40587)*JULIAN_S_I+time_t(fDayNS/SECOND_NS_U);
  return *t;
}

/**
 * Convert VATime to a POSIX (UNIX) time structure, with microsecond
 * accuracy.
 **/ 
inline void VATime::getPOSIXTimeVal(struct timeval* tv) const
{
#ifdef __APPLE__
//#warning Fix for POSIX non-complient Mac/BSD should be removed when no longer needed
  time_t tv_sec;
  getPOSIXTimeT(&tv_sec);
  tv->tv_sec = tv_sec;
#else
  getPOSIXTimeT(&tv->tv_sec);
#endif
  tv->tv_usec=suseconds_t((fDayNS%SECOND_NS_U)/1000);
}

/**
 * Convert VATime to a calendar date.
 **/ 
inline void VATime::getCalendarDate(uint32_t& y, 
                                    uint32_t& m, 
                                    uint32_t& d) const
{
  // Julian Date To Calendar Date Conversion Algorithm
  // Reference: John Meeus, Astronomical Algorithms, Second Edition, 
  // Willmann-Bell, 1998, Chapter 7
  uint32_t ALPHA = uint32_t(floor((double(fMJD)+532784.75)/36524.25));
  uint32_t A     = fMJD+2400002+ALPHA-ALPHA/4;
  uint32_t B     = A+1524;
  uint32_t C     = uint32_t(floor((double(B)-122.1)/365.25));
  uint32_t D     = uint32_t(floor(365.25*double(C)));
  uint32_t E     = uint32_t(floor(double(B-D)/30.6001));
  d = B-D-uint32_t(floor(30.6001*double(E)));
  m = (E<14)?E-1:E-13;
  y = (m>2)?C-4716:C-4715;
}

/**
 * Convert VATime to a UTC time.
 **/ 
inline void VATime::getTime(uint32_t& h, 
                            uint32_t& m, 
                            uint32_t& s, 
                            uint32_t& ns) const
{
  if(fDayNS>=JULIAN_NS_U){
      // We have a leap second situation
      h=23;
      m=59;
      s=60+(fDayNS-JULIAN_NS_U)/SECOND_NS_U;
      ns=fDayNS%SECOND_NS_U;
      return;
  }

  h=fDayNS/HOUR_NS_U;
  m=(fDayNS/MIN_NS_U)%60;
  s=(fDayNS/SECOND_NS_U)%60;
  ns=fDayNS%SECOND_NS_U;
  return;
}

/**
 * Convert VATime to a calendar year.
 **/ 
inline uint32_t VATime::getYear() const
{
  uint32_t Y,M,D;
  getCalendarDate(Y,M,D);
  return Y;
}

/**
 * Convert VATime to a calendar month.
 **/ 
inline uint32_t VATime::getMonth() const
{
  uint32_t Y,M,D;
  getCalendarDate(Y,M,D);
  return M;
}

/**
 * Convert VATime to a calendar day.
 **/ 
inline uint32_t VATime::getDay() const
{
  uint32_t Y,M,D;
  getCalendarDate(Y,M,D);
  return D;
}

/**
 * Convert VATime to a calendar day of year.
 **/ 
inline uint32_t VATime::getDayOfYear() const
{
  uint32_t Y,M,D;
  getCalendarDate(Y,M,D);
  return calendarDateToDOY(Y, M, D);
}

/**
 * Convert VATime to a calendar day of week.
 **/ 
inline uint32_t VATime::getDayOfWeek() const
{
  // Julian Date to Day Of Week Conversion Algorithm
  // Reference: John Meeus, Astronomical Algorithms, Second Edition, 
  // Willmann-Bell, 1998, Chapter 7
  // 
  // Note: (JD+1.5)%7 == (MJD+2400002)%7 == (MJD+3)%7
  return (fMJD+3)%7;
}

/**
 * Convert VATime to UTC hour.
 **/ 
inline uint32_t VATime::getHour() const
{
  uint32_t H,M,S,NS;
  getTime(H,M,S,NS);
  return H;
}

/**
 * Convert VATime to UTC minute.
 **/ 
inline uint32_t VATime::getMin() const
{
  uint32_t H,M,S,NS;
  getTime(H,M,S,NS);
  return M;
}

/**
 * Convert VATime to UTC second.
 **/ 
inline uint32_t VATime::getSec() const
{
  uint32_t H,M,S,NS;
  getTime(H,M,S,NS);
  return S;
}

/**
 * Convert VATime to UTC nanosecond.
 **/ 
inline uint32_t VATime::getNanoSec() const
{
  uint32_t H,M,S,NS;
  getTime(H,M,S,NS);
  return NS;
}

/**
 * Convert VATime to an integer time stamp with millisecond accuracy, in
 * the form YYYYMMDDhhmmssfff.
 **/ 
inline uint64_t VATime::getMSTimeStamp() const
{
  uint32_t y,m,d;
  getCalendarDate(y,m,d);
  uint32_t H,M,S,NS;
  getTime(H,M,S,NS);

  uint64_t timestamp =
    y   * 10000000000000ULL
    + m *   100000000000ULL
    + d *     1000000000ULL
    + H *       10000000ULL
    + M *         100000ULL
    + S *           1000ULL
    + NS / 1000000ULL;
  
  return timestamp;
}

/**
 * Convert VATime to an integer time stamp with second accuracy, in the
 * form YYYYMMDDhhmmss.
 **/ 
inline uint64_t VATime::getDBTimeStamp() const
{
  uint32_t y,m,d;
  getCalendarDate(y,m,d);
  uint32_t H,M,S,NS;
  getTime(H,M,S,NS);

  uint64_t timestamp =
    y   * 10000000000ULL
    + m *   100000000ULL
    + d *     1000000ULL
    + H *       10000ULL
    + M *         100ULL
    + S;
  
  return timestamp;
}

/**
 * Convert VATime to a stringified time stamp in the form
 * "YYYY-MM-DD hh:mm:ss.fffffffff"
 **/ 
inline void VATime::getString(std::string& str) const
{
  uint32_t Y,M,D,h,m,s,ns;
  getCalendarDate(Y,M,D);
  getTime(h,m,s,ns);
  std::ostringstream stream;
  stream << std::fixed << std::setfill('0')
	 << std::setw(4) << std::setprecision(4) << Y << '-' 
	 << std::setw(2) << std::setprecision(2) << M << '-' 
	 << std::setw(2) << std::setprecision(2) << D << ' '
	 << std::setw(2) << std::setprecision(2) << h << ':'
	 << std::setw(2) << std::setprecision(2) << m << ':'
	 << std::setw(2) << std::setprecision(2) << s << '.'
	 << std::setw(9) << std::setprecision(9) << ns;
  if((int(fTimeSystem) < (int)sTimeSystems.size())&&
     (!sTimeSystems[int(fTimeSystem)].fTimeSystemName.empty()))
    stream << " " << sTimeSystems[int(fTimeSystem)].fTimeSystemName;
  str = stream.str();
}

/**
 * Function which is not yet tested or documented
 **/ 
inline VATime::TimeSystem VATime::getTimeSystem() const
{
  return fTimeSystem;
}

// ============================================================================
// Mathematical Member Operators
// ============================================================================

/**
 * Add a nanosecond offset to the VATime.
 **/ 
inline VATime& VATime::operator += (const int64_t& ns)
{
  if(sIgnoreLeapSeconds)
    {
      // Faster code which ignores leapseconds
      int64_t offset_from_start_of_day = ns+int64_t(fDayNS);
      int64_t quot = offset_from_start_of_day/JULIAN_NS_I;
      int64_t rem = offset_from_start_of_day%JULIAN_NS_I;
      if(rem<0)quot--,rem+=JULIAN_NS_I;
      fMJD += quot;
      fDayNS = rem;
      return *this;
    }

  bool ts_good = (int(fTimeSystem) < (int)sTimeSystems.size());
  const std::map<uint32_t, int>* ts_tai = 
    ts_good?(&sTimeSystems[int(fTimeSystem)].fOffsetToTAI):0;
  const int ts_0 = 
    ts_good?(sTimeSystems[int(fTimeSystem)].fZeroOffsetToTAI):0;

  if(ns>0)
    {
      int64_t ns_remaining = ns;
      uint64_t day_len;

      std::map<uint32_t, int>::const_iterator i;
      if(ts_good)
	{
	  i = ts_tai->upper_bound(fMJD);
	  if((i==ts_tai->end())||(i->first!=fMJD+1))day_len=JULIAN_NS_U;
	  else 
	    {
	      std::map<uint32_t, int>::const_iterator j=i;
	      int t0 = (j==ts_tai->begin())?ts_0:(--j)->second;
	      day_len=uint64_t(JULIAN_S_I+i->second-t0)*SECOND_NS_U; 
	      i++;
	    }
	}
      else day_len=JULIAN_NS_U;

      if((uint64_t)ns_remaining <= (day_len-fDayNS))
	{
	  fDayNS += ns_remaining;
	  ns_remaining = 0;
	}
      else
	{
	  ns_remaining -= (day_len-fDayNS);
	  fDayNS=0;
	  fMJD++;
	}
      
      if(ts_good)
	while((ns_remaining)&&(i!=ts_tai->end()))
	  {
	    if(ns_remaining <= ((i->first-1 - fMJD)*JULIAN_NS_I))
	      {
		fMJD += ns_remaining/JULIAN_NS_I;
		fDayNS = ns_remaining%JULIAN_NS_I;
		ns_remaining=0;
		continue;
	      }
	    
	    ns_remaining -= (i->first-1 - fMJD)*JULIAN_NS_I;
	    fMJD = i->first-1;
	    
	    std::map<uint32_t, int>::const_iterator j=i;
	    int t0 = (j==ts_tai->begin())?ts_0:(--j)->second;
	    day_len=uint64_t(JULIAN_S_I+i->second-t0)*SECOND_NS_U;

	    if((uint64_t)ns_remaining <= day_len)
	      {
		fDayNS=ns_remaining;
		ns_remaining = 0;
		continue;
	      }
	    
	    ns_remaining -= day_len;
	    fMJD++;
	    i++;
	  }

      if(ns_remaining)
	{
	  fMJD += ns_remaining/JULIAN_NS_I;
	  fDayNS = ns_remaining%JULIAN_NS_I;
	}
    }
  else if(ns<0)
    {
      int64_t ns_remaining = -ns;

      if((uint64_t)ns_remaining <= fDayNS)
	{
	  fDayNS -= ns_remaining;
	  ns_remaining = 0;
	}
      else
	{
	  ns_remaining -= fDayNS;
	  fDayNS=0;
	}

      if(ts_good)
	{
	  uint64_t day_len;
	  std::map<uint32_t, int>::const_iterator i = 
	    ts_tai->upper_bound(fMJD);
	  
	  if(i!=ts_tai->begin())i--;
	  
	  while((ns_remaining)&&(i!=ts_tai->end())&&(i->first<=fMJD))
	    {
	      if(ns_remaining <= (fMJD-i->first)*JULIAN_NS_I)
		{
		  fMJD -= 1+(ns_remaining-1)/JULIAN_NS_I;
		  fDayNS = JULIAN_NS_I - ns_remaining%JULIAN_NS_I;
		  ns_remaining = 0;
		  continue;
		}
	      
	      ns_remaining -= (fMJD-i->first)*JULIAN_NS_I;
	      fMJD = i->first;
	      
	      std::map<uint32_t, int>::const_iterator j=i;
	      int t0 = (j==ts_tai->begin())?ts_0:(--j)->second;
	      day_len=uint64_t(JULIAN_S_I+i->second-t0)*SECOND_NS_U;
	      
	      if((uint64_t)ns_remaining <= day_len)
		{
		  fMJD--;
		  fDayNS = day_len - ns_remaining;
		  ns_remaining=0;
		  continue;
		}
	      
	      ns_remaining -= day_len;
	      fMJD--;
	      
	      if(i==ts_tai->begin())i=ts_tai->end();
	      else i--;
	    }
	}
      
      if(ns_remaining)
	{
	  fMJD -= 1+(ns_remaining-1)/JULIAN_NS_I;
	  fDayNS = JULIAN_NS_I - ns_remaining%JULIAN_NS_I;
	}
    }
  return *this;
}

/**
 * Subtract a nanosecond offset to the VATime.
 **/ 
inline VATime& VATime::operator -= (const int64_t& ns)
{
  int64_t nns = -ns;
  *this += nns;
  return *this;
}

/**
 * Subtract two VATimes, returning the number of nanoseconds between
 * them.
 **/ 
inline int64_t VATime::operator - (const VATime& o) const
{
  // Maximum time difference expressable as a uint64_t with nano second
  // resolution is 2^63 ns which is slightly more than 292.2 years

  int64_t answer=int64_t(fDayNS)-int64_t(o.fDayNS);

  if(fMJD != o.fMJD)
    {
      int leap_seconds = 0;
      if((!sIgnoreLeapSeconds)&&
	 (int(fTimeSystem) < (int)sTimeSystems.size()))
	{
	  const std::map<uint32_t, int>* ts_tai = 
	    &sTimeSystems[int(fTimeSystem)].fOffsetToTAI;
	  const int ts_0 = 
	    sTimeSystems[int(fTimeSystem)].fZeroOffsetToTAI;
	  
	  std::map<uint32_t, int>::const_iterator tai0 = 
	    ts_tai->upper_bound(((fMJD<o.fMJD)?fMJD:o.fMJD));
	  std::map<uint32_t, int>::const_iterator tai1 = 
	    ts_tai->upper_bound(((fMJD<o.fMJD)?o.fMJD:fMJD));

	  int t0 = ts_0;
	  if(tai0 != ts_tai->begin()) { tai0--; t0=tai0->second; }

	  int t1 = ts_0;
	  if(tai1 != ts_tai->begin()) { tai1--; t1=tai1->second; }

	  leap_seconds = t1-t0;

	  if(fMJD<o.fMJD)leap_seconds=-leap_seconds;
	}

      answer += 
	((int64_t(fMJD)-int64_t(o.fMJD))*JULIAN_S_I+leap_seconds)*SECOND_NS_I;
    }

  return answer;
}

// ============================================================================
// Boolean Operators
// ============================================================================

/**
 * Return true if time has been set correctly and if the GPS status was
 * good (if appropriate).
 **/ 
inline bool VATime::isGood() const
{
  if((!fDecodedOK)||(fGPSStatus!=0))return false;
  return true;
}

/**
 * Return true if time has been set correctly and if the GPS status was
 * good (if appropriate).
 **/ 
inline VATime::operator bool() const
{
  return isGood();
}

/**
 * Return true if two times are equal.
 **/ 
inline bool VATime::operator == (const VATime& o) const
{
  return((fMJD==o.fMJD)&&(fDayNS==o.fDayNS));
}

/**
 * Return true if first time is less than or equal to second.
 **/ 
inline bool VATime::operator <= (const VATime& o) const
{
  if(fMJD<o.fMJD)return true;
  else if(fMJD>o.fMJD)return false;
  return fDayNS<=o.fDayNS;
}

/**
 * Return true if first time is greater than or equal to second.
 **/ 
inline bool VATime::operator >= (const VATime& o) const
{
  if(fMJD>o.fMJD)return true;
  else if(fMJD<o.fMJD)return false;
  return fDayNS>=o.fDayNS;
}

/**
 * Return true if first time is strictly less than second.
 **/ 
inline bool VATime::operator < (const VATime& o) const
{
  if(fMJD<o.fMJD)return true;
  else if(fMJD>o.fMJD)return false;
  return fDayNS<o.fDayNS;
}

/**
 * Return true if first time is strictly greater than second.
 **/ 
inline bool VATime::operator > (const VATime& o) const
{
  if(fMJD>o.fMJD)return true;
  else if(fMJD<o.fMJD)return false;
  return fDayNS>o.fDayNS;
}

/**
 * Return true if two times are not equal.
 **/ 
inline bool VATime::operator != (const VATime& o) const
{
  return((fMJD!=o.fMJD)||(fDayNS!=o.fDayNS));
}

/**
 * Stream output insertion. Write stringified representation of time to
 * a stream.
 **/ 
inline std::ostream& operator << (std::ostream& stream, 
				  const VATime& gps)
{
  std::string str;
  gps.getString(str);
  stream << str;
  return stream;
}

/**
 * Static function to return new VATime set from the system clock.
 **/ 
inline VATime VATime::now()
{
  VATime n;
  n.setFromSystemTime();
  return n;
}

// ============================================================================
// Non-member Mathematical Operators
// ============================================================================

inline VATime operator + (const VATime& x, int64_t ns)
{
  VATime t(x);
  return t+=ns;
}

inline VATime operator + (int64_t ns, const VATime& x)
{
  VATime t(x);
  return t+=ns;
}

inline VATime operator - (const VATime& x, int64_t ns)
{
  VATime t(x);
  return t-=ns;
}

#endif // VATIME_H
