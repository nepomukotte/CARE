/*  GOrderedGrid class for creating an ordered 2D
     hash table, e.g. used for quickly finding element number
     given hit location on the grid. Grid elements are ordered
     by maximum area overlap with the grid bin.
     Within the GrISU simulation code, this search code was at least
     seven times faster than a simple looping search over all elements.
 */

#ifndef GORDEREDGRID
#define GORDEREDGRID

#include <list>

/*! \brief structure for each grid element.

  Each grid bin holds an stl list of these GridElem's.

 */
struct GridElem { 
  int elemNum;  //!< element number, e.g. pixel or facet number
  double dist;  //!< distance of element from grid center

  /*!  \brief constructor, often used to set the elemNum and dist.
   */
  GridElem(const int &elemNum=0, const double &dist=0.0);
};

/*!  \brief makes ordered hash table

  Given a set of elements, e.g. pixels, with positions, radii, shape,
  and the number of sides as well as the number of bins in the x and y 
  directions, GOrderedGrid produces an ordered hash table where a given 
  input x,y hit location addresses the appropriate bin and its list of 
  elements that fully or partially overlap the grid bin. The elements in 
  the list of ordered by decreasing area of overlap. Thus, the x,y hit 
  location is most likely to fall within the range of the first element
  in the list, etc.  Element shapes may be circles, squares, or hexagons.
  Extensions to other regular polygons is very easy.
 
 */
class GOrderedGrid {
  
  vector<double> vmx;  //!< element x location
  vector<double> vmy;  //!< element y location
  vector<double> vmr;  //!< element radius
  int nbinsx;          //!< number of grid bins in x direction
  int nbinsy;          //!< number of grid bins in y direction
  int iGridOption;     //!< see details in constructor documentation

  bool bGridOption;    //!< true if grid is active

  string sFileName;    //!< grid filename

  // Grid parameters
  double fXmin;  //!< min. x that encloses all element areas
  double fXmax;  //!< max. x that encloses all element areas
  double fYmin;  //!< min. y that encloses all element areas
  double fYmax;  //!< max. y that encloses all element areas
  double fDelX;  //!< x grid spacing
  double fDelY;  //!< y grid spacing

  // grid, vector of lists of elements
  vector< list<GridElem> * > vGrid;  //!< ordered grid vector

  /*!  \brief initialize grid parameters
   */
  void initialize();

  /*! \brief make the parameters for the grid
   */
  bool makeGridParameters();

  /*! \brief construct the grid
   */
  bool makeGrid();

  /*! \brief read the grid from a text file
   */
  bool readGrid();

  /*! \brief print the grid to a text file
   */
  bool printGrid();

 public:

  /*!  \brief constructor

    \param xe vector of element x locations
    \param ye vector of element x locations
    \param re vector of element x locations
    \param nbinsx  number of bins in x direction
    \param nbinsy  number of bins in y direction
    \param option grid options 
    \param filename name of gridfile (see options)
    *
    *   option = 0 no grid is created <BR>
    *   option = 1      1, make grid        <BR>
    *   option = 2 make grid and write to file <BR>
    *   option = 3 read from file; if no file exists, go to option 2 <BR>

   */
  GOrderedGrid(const vector<double> &xe,const vector<double> &ye,
               const vector<double> &re,
               const int &nbinsx, const int &nbinsy, const int &option = 1,
               const string &filename = "");

  /*! \brief deconstructor
   */
  ~GOrderedGrid();

  /*!  \brief copy constructor
   */
  GOrderedGrid(const GOrderedGrid &orderedGrid);

  /*!   \brief getGridBinList returns the list of grid elements

    \param x  hit x location in grid coordinates
    \param y  hit y location in grid coordinates
    \gridList list of GridElem (see above) for grid bin addressed by x,y
    \return   true if GridList has at least one element
  
    Create "list<GridElem> *gridList" in the calling program, this method
    will assign the gridList pointer.

   */
  bool getGridBinList(const double &x, const double &y,
                 list<GridElem> *&gridList);


  /*! \brief  converts a string to a vector of tokens using "space"
    as a delimiter

   */
  void tokenizer(const string& str, 
                 vector<string>& tokens);

};

#endif
