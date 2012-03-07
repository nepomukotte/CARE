/*  GOrderedGridSearch class for determining 
  element number using an ordered 2D hash table produced by
  an instance of GOrderedGrid.  Increases search speed by about
  a factor of 7 over sequential loop over elements.
 */

#ifndef GORDEREDGRIDSEARCH
#define GORDEREDGRIDSEARCH

// forward declarations
class GOrderedGrid;

/*! \brief Determines element number using an ordered 2D hash table produced by
  an instance of GOrderedGrid.

  Increases search speed by about
  a factor of 7 over sequential loop over elements.

 */
class GOrderedGridSearch {


  GOrderedGrid *pGrid;   //!< GOrderedGrid instance composed within this class

  const vector<double> *vxe;   //!< x location of elements
  const vector<double> *vye;   //!< y location of elements
  const vector<double> *vre;   //!< radius of elements
  const vector<int> *vse;      //!< number of sides of elements
  const vector<double> *vrote; //!< rotation angle of elements (radians)

  int nbinsx;           //!< number of grid bins in x direction
  int nbinsy;           //!< number of grid bins in y direction
  int iGridOption;      //!< see documentation in GOrderedGrid class
  bool bGridOption;     //!< true if grid is active
  string fileNameGrid;  //!< grid filename 

  /*!  \brief Determines if x,y hit location relative to center
    of element falls within the polygon defining the element

    \param sides  number of sides of element
    \param alpha rotation angle of element (radians) 
    \param radius radius of element
    \param x x-location of hit, relative to element center
    \param y y-location of hit, relative to element center
    \return true if hit falls within the element polygon (or circle)
   */
  bool polyInside(const int &sides, const double &alpha, 
                  const double &radius, 
                  const double &x, const double &y);
 public:
  
  /*! \brief Constructor

    \param xe vector of x locations of elements
    \param ye vector of y locations of elements
    \param re vector of radii of elements
    \param se vector of number of sides of element (1 for circle)
    \param rote vector of rotation angles (radians) for elements
    \param inbinsx number of grid bins in x direction 
    \param inbinsy number of grid bins in y direction 
    \param option option as documented in GOrderedGrid.h
    \param filename grid filename (see GOrderedGrid.h)
   */
  GOrderedGridSearch(const vector<double> &xe,const vector<double> &ye,
                     const vector<double> &re,const vector<int> &se,
                     const vector<double> &rote,
                     const int &inbinsx, const int &inbinsy, 
                     const int &option = 1,
                     const string &filename = "");
  

  /*!  \brief GOrderedGridSearch destructor
   */
  ~GOrderedGridSearch();   
  
  /*! /brief eturns the element number containing x,y

    \param hit location in x direction (relative to grid center)
    \param hit location in y direction (relative to grid center)
    \return element number of element containing hit or -1.
   */
  int getElemNumber(const double &x,const double &y);
 
};

#endif

