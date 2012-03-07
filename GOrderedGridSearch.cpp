/*  GOrderedGridSearch.cpp

    June 17, 2011
    Charlie Duke
    Grinnell College
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <deque>
#include <list>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include <limits>

using namespace std;

#include "TMath.h"

#include "GOrderedGrid.h"
#include "GOrderedGridSearch.h"

#define DEBUG(x) cout << #x << " = " << x << endl
#define DEBUGS(x) cout << "      " << #x << " = " << x << endl

GOrderedGridSearch::GOrderedGridSearch(const vector<double> &xe,
                                       const vector<double> &ye,
                                       const vector<double> &re,
                                       const vector<int> &se,
                                       const vector<double> &rote,
                                       const int &inbinsx, const int &inbinsy,
                                       const int &option,
                                       const string &filename) 
 {

   bool debug = false;
   if (debug) {
     cout << "  -- GOrderedGridSearch::GOrderedGridSearch " << endl;
   }

  // initialize pointers
  pGrid = 0;

  vxe = &xe;
  vye = &ye;
  vre = &re;
  vse = &se;
  vrote = &rote;

  nbinsx = inbinsx;
  nbinsy = inbinsy;
  iGridOption = option;
  fileNameGrid = filename;

  if (iGridOption > 0)
    bGridOption = true;
  else {
    bGridOption = false;
  }

  if (debug) {
    cout << "       input vector sizes " << endl;
    cout << "       vxe " << vxe->size() << endl;
    cout << "       vye " << vye->size() << endl;
    cout << "       vre " << vre->size() << endl;
    cout << "       vse " << vse->size() << endl;
    cout << "       vrote " << vrote->size() << endl;
    cout << "       inbinsx inbinsy " << inbinsx << "  " << inbinsy << endl;
    cout << "       iGridOption " << iGridOption << endl;
    cout << "       bGridOption " << bGridOption << endl;
    cout << "       fileNameGrid " << fileNameGrid << endl;

  }

  // make or read the element grid if necessary
  if (iGridOption != 0) {
    pGrid = new GOrderedGrid(*vxe,*vye,*vre,nbinsx,nbinsy,
                             iGridOption,fileNameGrid);
  }
};
//************************  end of GOrderedGridSearch ****************

GOrderedGridSearch::~GOrderedGridSearch() {
  delete pGrid;

};
//************************  end of ~GOrderedGridSearch ****************

int GOrderedGridSearch::getElemNumber(const double &x,const double &y) {

  bool debug = false;
  if (debug) {
    cout << "  -- GOrderedGridSearch::getElemNumber " << endl;
    cout << "         x  y  " << x << "  " << y << endl;
  }

  // pointer for gridList
  list<GridElem> *gridList = 0;
  list<GridElem>::iterator gridIter;

  int maxElem;  /* maximum number of elements for looping, 
                  either all elements or grid elements */
  bool gridFound;  // if true, grid found for x,y

  if (bGridOption) {

    gridFound = pGrid->getGridBinList(x,y,gridList);
    
    if (debug) {
      cout << "     gridFound " << gridFound << endl;
      if (gridFound) {
        cout << "     number of grid elements " 
             << gridList->size() << endl;
        cout << "         element num ";
        for (gridIter = gridList->begin();
             gridIter != gridList->end(); gridIter++) {
          cout << (*gridIter).elemNum << " ";
        }
        cout << endl;    
      }
    }
  }
  // we now have a group of grid elements to loop over, but build
  // in loop over all elements for testing.
  
  if (!bGridOption) {
    // will loop over all elements
    maxElem = vxe->size();
  }
  else {
    if (gridFound) {
      maxElem = gridList->size();  // number of elements in the gridbin
      gridIter = gridList->begin();
    }
    else {
      if (debug) {
        cout << "      points x,y outside the grid, return -1" << endl;       
      }
      return -1;
    }
  }
  
  // ready for looping over elements
  int pNum;
  if (debug) cout << "      starting loop over elements " << endl;
  for (int k1 = 0;k1<maxElem;k1++) {
    
    // set up element number for this value of k1
    if (!bGridOption) {
      pNum = k1;   // looping over all elements
    }
    else {
      GridElem gElem = *gridIter;
      pNum = gElem.elemNum;
      if (debug) {
        cout << "     pNum  " << pNum << "  "
             << gElem.dist << endl;  // dist is to bin center from x,y
      }
      gridIter++;
    }

    double xpp = x - (*vxe)[pNum];
    double ypp = y - (*vye)[pNum];
    double r = sqrt(xpp*xpp + ypp*ypp);
    if (debug) {
      cout << "     checking distance for elem " << pNum << endl;
      DEBUGS( (*vxe)[pNum] ); DEBUGS((*vye)[pNum]);
      DEBUGS(xpp); DEBUGS(ypp); DEBUGS(r);
    }
       
    if ( r < (*vre)[pNum] ) {
      // found an element within reach
      if (debug) {
        cout << "   FOUND ELEMENT WITHIN REACH: elem.num.  " 
             << pNum << endl;
        DEBUGS((*vxe)[pNum]); DEBUGS( (*vye)[pNum]);       
      }
      // need to test to see if it's inside the element polygon
      int sides = (*vse)[pNum];
      double alpha = (*vrote)[pNum];
      double radius = (*vre)[pNum];
      
      if (debug) {
        cout << "       entering polyInside with these parms" << endl;
        DEBUGS(sides);DEBUGS(alpha);DEBUGS(radius);
        DEBUGS(xpp);DEBUGS(ypp);
      }
      bool polyTest = polyInside(sides,alpha,radius,xpp,ypp);
      if (debug) cout << "      polyTest " << polyTest << endl;

      if (polyTest) {
        return pNum;
      }
      
    }
    
  }

return -1;
};
//************************  end of getElemNumber ****************

bool GOrderedGridSearch::polyInside(const int &sides, const double &alpha, 
                                    const double &radius, 
                                    const double &x, const double &y) {
  
  double PI = (TMath::Pi());   /*   Value of pi */

  double theta,dist;
  double phi;
  
  static double delta[50],beta[50],gamma[50];
  static int ifirst=0;
  int iside;
  double x1,y1;
  
  /* initialize constant parameters */
  if (ifirst==0) {
    ifirst = 1;
    delta[0]=delta[1]=delta[2]=0.0;
    beta[0] = beta[1]= beta[2]=0.0;
    gamma[0]=gamma[1]=gamma[2]=0.0;

    for (iside=3;iside<50;iside++) {
      delta[iside] = (PI)*( 0.5 - (1/(double)iside) );
      gamma[iside] = 2*PI/(double)iside;
      beta[iside]  = gamma[iside]/2.;
    }
  }
  /* number of valid sides determined by size of arrays 
     I only went to 49.............. Maybe that's a circle....
  */
  if (sides > 49) {
    cout << " number of sides greater than 49! " << endl;
    cout << " exiting code" << endl;
    exit(0);
  }

  /* if at least three sides, can proceed */
  /*-------------------------------------------*/
  if (sides > 2) {

    dist =  sqrt((x*x) + (y*y));

    /* first see if the distance isn't too large */
    if (dist > radius) return 0;

    /* rotate x and y through alpha, dist doesn't change */
    x1 = (cos(alpha)*x)- (sin(alpha)*y);
    y1 = (sin(alpha)*x)+ (cos(alpha)*y);

    /* if you use atan to get theta, phi will be incorrect */
    theta = acos(y1 / dist);
  
    phi = fabs( (fmod((theta + beta[sides]), gamma[sides])) - beta[sides]); 

    if (dist<=((sin(delta[sides])*radius)/sin(delta[sides] + phi))) 
      return true;
    
  }
  /*----------------------------------------------*/
  /* if sides = 1 or 2, then shape is a circle rather than a polygon */

  else if (sides > 0 ) {
    dist =  sqrt(x * x + y * y);
    if (dist <= radius) return 1;
  }

  /* if sides = 0, then always return 0 */
  return false;
};
/********************  end of polyInside ***************/
