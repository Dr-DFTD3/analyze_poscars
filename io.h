#ifndef _IO_H_
#define _IO_H_
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <fstream>
#include <iostream>
#include <iomanip>

#include <sstream>
#include <vector>
#include "common.h"

/*------------------------------------------------------------
    Write each atom and its nearest neighbors to a
    POSCAR file
--------------------------------------------------------------*/
void
print_atoms_neighbors( std::vector<point_t>&     atoms,
                       double**                  H);

/*------------------------------------------------------------
    Scan list of "atomClusters" and write the coordinates
    of each atom in a cluster to a single POSCAR file
--------------------------------------------------------------*/
void
print_atom_clusters( std::vector<point_t>&     atoms,
                     double**                  H,
                     std::vector<int>          nc,
                     param_t*                  param );

/*------------------------------------------------------------
    Print the neighbor list and bond lengths for each
    atom in the system	
--------------------------------------------------------------*/
void 
print_bonds_nnlist( std::vector<point_t>&     atoms, 
					std::vector<std::string>& species,
					param_t*                  p);

/*------------------------------------------------------------
    Scan list of "atomClusters" and identify what molecules 
    are present.
--------------------------------------------------------------*/
void 
print_molecules( std::vector<std::vector<int> > atomClusters,
                 std::vector<point_t>&          atoms,
                 param_t*                       p);


/*------------------------------------------------------------
    Write a single POSCAR/CONTCAR file. 
--------------------------------------------------------------*/
int 
write_poscar( char*                     name,
              int                       ix,
              std::vector<point_t>&     atoms,
              double**                  H, 
              std::vector<std::string>& species, 
              std::vector<int>&         ntypes);

/*------------------------------------------------------------
    Read a POSCAR/CONTCAR file. 
    Store lattice vectors in "H"
    Store atom positions in "atoms"
    Store unique species in "species"
    Store type counts in "ntypes"
    return 1 if successful, -1 if not
--------------------------------------------------------------*/
int
read_poscar( char*                     inFile, 
			 std::vector<point_t>&     atoms, 
			 double**                  H, 
			 std::vector<std::string>& species, 
			 std::vector<int>&         ntypes );


#endif