#ifndef _COMMON_H_
#define _COMMON_H_
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <fstream>
#include <iostream>
#include <iomanip>

#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <sys/stat.h>

#define D2R 3.141592653589793238462643/180
#define R2D 180.0/3.141592653589793238462643

#define MIN(A,B) ((A) < (B) ? (A) : (B))

#define EXE "fmol"
#define SOFT "findMolecules"
#define VERS "0.0.5"
#define MAX_STRING_LENGTH 512


struct PeriodicTable{
    std::string element[110];
    double amu[110];
    double radius[110];
};

// structure to hold runtime parameters
typedef struct {

    double rcut;
    double minBond;
    double tolerance;


    int    verb;
    int    printNlist;
    int    printMolecules;
    int    save;
    int    nbonds;
    int    printClusterInfo;

    char*  prefix;
    char*  inFile;
    char*  minPair;

}param_t;


typedef struct {
    double x,y,z;
    double xs,ys,zs;

    int    gid;
    int    id;
    int    tid;
    int    clusterID;
    std::vector<double> rnn;
    std::vector<int> naborIds;
    std::vector<std::vector<int> > tvecs;
    std::string symb;

}point_t;



/*------------------------------------------------------------
    Convert a std::string to and number Type
--------------------------------------------------------------*/
template<typename retVal>
retVal 
str2num(const std::string& numstr)
{
  retVal number;
  std::stringstream stream(numstr);
  stream >> number;
  return number;
}


/*------------------------------------------------------------
    Convert any number Type to a std::string
--------------------------------------------------------------*/
template<typename Type>
std::string num2str(Type number)
{
    std::stringstream ss;
    ss << number;
    std::string str = ss.str();
    return str;
}

/*------------------------------------------------------------
    Convert a std::string to char*
--------------------------------------------------------------*/
char* 
str2char(std::string str);

/*------------------------------------------------------------
    Formulate a name based on the stoichiometry of the system
--------------------------------------------------------------*/
char* 
base_name( std::vector<std::string> species, 
           std::vector<int>         ntypes );

/*------------------------------------------------------------
    Scan periodic table to get the covalent radius of 
    "element". Return  covalent radius in Å
--------------------------------------------------------------*/
double 
elem_radii( std::string element );

/*------------------------------------------------------------
    Define a cutoff for a pair of atoms based on the
    sum of the covalent radius. dr represents a fudge factor
    in determining a bond.
--------------------------------------------------------------*/
double 
sum_covalent_radii( std::string atomi, 
                    std::string atomj, 
                    double dr);

/*------------------------------------------------------------
    Sort the list of atoms based on type ID, this is required
    by the format of POSCAR/CONTCAR files
--------------------------------------------------------------*/
void 
sort_ids( std::vector<point_t>& atoms );

/*------------------------------------------------------------
    Checks if a given character or string of characters is 
    a regular file which can be opened and read
------------------------------------------------------------*/
int 
is_regular_file(const char *ch);

/*------------------------------------------------------------
    Parse the command line arguments provided.
    Return 1 if successful, 0 if not.
--------------------------------------------------------------*/
int
get_arguments( int      nargs,
               char**   argv, 
               param_t* p);

void 
show_help();

void 
show_version();

void 
show_usage();


/*------------------------------------------------------------
    Sort a 2D vector based on the second column
--------------------------------------------------------------*/
struct col2_asc
{
    template<class T1>
    bool operator()(const std::vector<T1> &v1,const std::vector<T1> &v2) const { return v1[1] < v2[1]; }
};

inline 
int 
gcd(int a, int b) {return b == 0 ? a : gcd(b, a % b);}

char* 
stoichiometry( std::vector<std::string> species, 
               std::vector<int>         ntypes );


#endif
