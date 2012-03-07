//-*-mode:c++; mode:font-lock;-*-

/**
 * \file VATime.cpp
 * \ingroup common
 *
 * Time conversion/manipulation class for VEGAS
 *
 * Original Author: Stephan Fegan
 * $Author: nepomuk $
 * $Date: 2010/05/25 23:20:20 $
 * $Revision: 1.1 $
 * $Tag$
 *
 **/

/**
 * \class VATime
 * \ingroup common
 * \brief Time conversion/manipulation class for VEGAS
 *
 * VATime is a general purpose time class for VEGAS. The class can
 * convert between many different time representations, such as: MJD,
 * calendar time, Symmetricon binary GPS format, UTC time strings and
 * many more. It provides arithmetic operations, such as subtraction
 * of two times to give the interval between them, addition of an
 * interval to a time and logical comparissions between two times.
 *
 * VATime aims to handle leap seconds correctly, although that
 * functionality is not yet tested.
 *
 **/

// These deinitions tell the makefile which library the cpp file
// should be included in
// VA_LIBRARY_TAG: libSP24common.a
// VA_LIBRARY_TAG: libSP24commonLite.a

// This is where the code starts.

#include "VATime.h"

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class VATime;
#endif

#ifndef NOROOT
ClassImp(VATime);
#endif

VATime::~VATime()
{
  // nothing to see here
}

bool                                      VATime::sIgnoreLeapSeconds=true;
std::vector<VATime::TimeSystemDefinition> VATime::sTimeSystems;

// Structure containing the difference in time between TAI and UTC. The
// "delta" parameter defines the difference between Coordinated
// Universal Time (UTC) and International Atomic Time (TAI) as
// specified in http://hpiers.obspm.fr/iers/bul/bulc/bulletinc.dat.
//
// delta=TAI-UTC
//
// The "mjd" parameter defines the Modified Julian Date (MJD) at which
// the value of "delta" arose, through the insertion of a positive or
// negative leap second in the previous day. The "date" parameter is
// for convenience only.

// For example on "1999-01-01" TAI-UTC=32. This change happened
// through an insertion of a leap second with UTC as follows:
//
// 1998-12-31 23:59:59 
// 1998-12-31 23:59:60 (*LEAP SECOND*)
// 1998-12-31 23:59:00

struct TAI_UTC_Difference
{
  int delta;
  uint32_t mjd;
  const char* date;
};

void VATime::configure(/* VSOptions* options, ConfigInfo* info */)
{
  // By definition the time difference between TAI and TAI is ALWAYS 0
  setZeroOffsetToTAI(0, "TAI", TS_TAI);

  // TAI-UTC adapted from ftp://tycho.usno.navy.mil/pub/series/ser14.txt
  TAI_UTC_Difference delta[] = 
    { { 10, 36204, "01-01-58" },
      { 11, 41499, "07-01-72" },
      { 12, 41683, "01-01-73" },
      { 13, 42048, "01-01-74" },
      { 14, 42413, "01-01-75" },
      { 15, 42778, "01-01-76" },
      { 16, 43144, "01-01-77" },
      { 17, 43509, "01-01-78" },
      { 18, 43874, "01-01-79" },
      { 19, 44239, "01-01-80" },
      { 20, 44786, "07-01-81" },
      { 21, 45151, "07-01-82" },
      { 22, 45516, "07-01-83" },
      { 23, 46247, "07-01-85" },
      { 24, 47161, "01-01-88" },
      { 25, 47892, "01-01-90" },
      { 26, 48257, "01-01-91" },
      { 27, 48804, "07-01-92" },
      { 28, 49169, "07-01-93" },
      { 29, 49534, "07-01-94" },
      { 30, 50083, "01-01-96" },
      { 31, 50630, "07-01-97" },
      { 32, 51179, "01-01-99" },
      { 33, 53736, "01-01-06" } 
    };

  setZeroOffsetToTAI(delta[0].delta, "UTC", TS_UTC);
  for(unsigned i=1; i<sizeof(delta)/sizeof(*delta); i++)
    setOffsetToTAI(delta[i].mjd, delta[i].delta, TS_UTC);

  // By definition the time difference between TAI and GPS is constant
  // at the value of the TAU-UTC offset on 1980-01-06 (MJD=44244)
  setZeroOffsetToTAI(getOffsetToTAI(44244, TS_UTC), "GPS", TS_GPS);
}

#ifdef TEST_MAIN

#include <iostream>
#include <set>
#include <cstdlib>

int main(int argc, char**argv)
{
#if 0
  for(unsigned round=0; round<1000; round++)
    {
      VATime::clearAllOffsets(VATime::TS_UTC);

      int N = rand()%1000;
      std::set<uint32_t> mjds;
      for(unsigned ls_num=0; ls_num<N; ls_num++)mjds.insert(50000+rand()%1000);
      int delta = rand()%20-10;
      VATime::setZeroOffsetToTAI(delta, "SJF", VATime::TS_UTC);
      for(std::set<uint32_t>::const_iterator i=mjds.begin();i!=mjds.end();i++)
	{
	  delta += rand()%20-10;
	  VATime::setOffsetToTAI(*i, delta, VATime::TS_UTC);
	}

      for(unsigned i=0;i<10000;i++)
	{
	  uint32_t mjd = 49800+rand()%1400;
	  uint64_t ns = 0;
	  for(unsigned j=0;j<5;j++)
	    {
	      ns = ns*uint64_t(RAND_MAX) + uint64_t(rand());
	      ns = ns%VATime::getDayLength(mjd);
	    }
	  VATime t1;
	  t1.setFromMJDIntAndNS(mjd, ns);


	  mjd = 49800+rand()%1400;
	  ns = 0;
	  for(unsigned j=0;j<5;j++)
	    {
	      ns = ns*uint64_t(RAND_MAX) + uint64_t(rand());
	      ns = ns%VATime::getDayLength(mjd);
	    }
	  VATime t2;
	  t2.setFromMJDIntAndNS(mjd, ns);

	  int64_t ns0 = t2-t1;

	  VATime t3 = t1;
	  t3 += ns0;

	  VATime t4 = t2;
	  t4 -= ns0;

	  int64_t ns1 = t3-t1;
	  int64_t ns2 = t2-t4;

	  if((t1!=t4)||(t2!=t3))
	    {
	      VATime::dumpTAIOffsets(std::cout);
	      std::cout << std::endl;
	      if(t1<=t2)
		{
		  std::cout << ns0 << ' ' << ns1 << ' ' << ns2 << ' ' << (t4-t1)/1000000000LL << ' ' << (t3-t2)/1000000000LL << std::endl
			    << t1 << std::setw(5) << "  " << t1.getMJDInt() << "   " << std::setw(5) << t1.getDayNS()/1000000000ULL << std::endl 
			    << t4 << std::setw(5) << "  " << t4.getMJDInt() << "   " << std::setw(5) << t4.getDayNS()/1000000000ULL <<  "   " << (t4==t1) <<  std::endl
			    << t2 << std::setw(5) << "  " << t2.getMJDInt() << "   " << std::setw(5) << t2.getDayNS()/1000000000ULL << std::endl 
			    << t3 << std::setw(5) << "  " << t3.getMJDInt() << "   " << std::setw(5) << t3.getDayNS()/1000000000ULL <<  "   " << (t3==t2) << std::endl 
			    <<  std::endl;
		}
	      else
		{
		  std::cout << ns0 << ' ' << ns1 << ' ' << ns2 << ' ' << (t4-t1)/1000000000LL << ' ' << (t3-t2)/1000000000LL << std::endl
			    << t2 << std::setw(5) << "  " << t2.getMJDInt() << "   " << std::setw(5) << t2.getDayNS()/1000000000ULL << std::endl 
			    << t3 << std::setw(5) << "  " << t3.getMJDInt() << "   " << std::setw(5) << t3.getDayNS()/1000000000ULL <<  "   " << (t3==t2) << std::endl 
			    << t1 << std::setw(5) << "  " << t1.getMJDInt() << "   " << std::setw(5) << t1.getDayNS()/1000000000ULL << std::endl 
			    << t4 << std::setw(5) << "  " << t4.getMJDInt() << "   " << std::setw(5) << t4.getDayNS()/1000000000ULL <<  "   " << (t4==t1) <<  std::endl
			    <<  std::endl;
		}
	    }
	}
    }
#else

#if 1

  VATime::configure();
  VATime t1;
  t1.setFromString("1999-01-01 00:00:00.000000000 GPS");
  std::cout << t1 << std::endl;
  t1.convertTimeToSystem(VATime::TS_UTC);
  std::cout << t1 << std::endl;
  t1.convertTimeToSystem(VATime::TS_TAI);
  std::cout << t1 << std::endl;
  t1.convertTimeToSystem(VATime::TS_GPS);
  std::cout << t1 << std::endl;
  t1.convertTimeToSystem(VATime::TS_UTC);
  std::cout << t1 << std::endl;
  t1.convertTimeToSystem(VATime::TS_GPS);
  std::cout << t1 << std::endl;
  t1.convertTimeToSystem(VATime::TS_TAI);
  std::cout << t1 << std::endl;
  t1.convertTimeToSystem(VATime::TS_UTC);
  std::cout << t1 << std::endl;


#else
  VATime t1;
  t1.setFromMJDIntAndNS(50000, 500000000000ULL);
  VATime t2;
  t2.setFromMJDIntAndNS(50005, 500000000001ULL);
  
  std::cout << t2-t1 << std::endl;

  VATime::setZeroOffsetToTAI(10, "SJF", VATime::TS_UTC);
  std::cout << t2-t1 << std::endl;

  VATime::setOffsetToTAI(50010, 10, VATime::TS_UTC);
  std::cout << t2-t1 << std::endl;

  VATime::setOffsetToTAI(40010, 10, VATime::TS_UTC);
  std::cout << t2-t1 << std::endl;


  VATime::setOffsetToTAI(50001, 10, VATime::TS_UTC);
  std::cout << t2-t1 << std::endl;

  VATime::setOffsetToTAI(50001, 12, VATime::TS_UTC);
  std::cout << t2-t1 << std::endl;

  VATime::setOffsetToTAI(50000, 8, VATime::TS_UTC);
  VATime::setOffsetToTAI(50003, 20, VATime::TS_UTC);
  VATime::setOffsetToTAI(50004, 5, VATime::TS_UTC);
  std::cout << t2-t1 << std::endl;

  VATime t3 = t1;
  t3 += (t2-t1);
  
  std::cout << t1 << "  " << std::setw(5) << t1.getMJDInt() << "   " << std::setw(5) << t1.getDayNS()/1000000000ULL << std::endl;
  std::cout << t2 << "  " << std::setw(5) << t2.getMJDInt() << "   " << std::setw(5) << t2.getDayNS()/1000000000ULL << std::endl;
  std::cout << t3 << "  " << std::setw(5) << t3.getMJDInt() << "   " << std::setw(5) << t3.getDayNS()/1000000000ULL << std::endl;
#endif

#endif

  return EXIT_SUCCESS;
}

#endif
