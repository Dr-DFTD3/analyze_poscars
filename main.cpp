#include "common.h"
#include "neighbors.h"
#include "lattice.h"
#include "io.h"
#include "mymath.h"


/*------------------------------------------------------------
 *
 *  fMol is a simple utility to identify known  
 *  and unknown molecules present in an atomic 
 *  structure file. Structure file must be in 
 *  the form of a standard POSCAR/CONTCAR file.
 *         
 *  List of known molecules:
 *         H2,O2,N2,C2,OH,CO,NO,H2O,NO2,CO2,CH4,CH2,Azide
 * 
 *  Functionality: - if known print molecule names, 
 *                 - uniform pairwise cutoff or covalent radii
 *                 - save each molecule/cluster to a POSCAR
 *                 - print bond lengths
 *                 - print atoms and their neighbor lists
 *
 *  Author:        Joseph M. Gonzalez
 *  Date:          June 6, 2017
 *  Contact:       jmgonza6@mail.usf.edu           
 *  Version:       0.5
 *
--------------------------------------------------------------*/


int verb;

int 
main( int nargs, 
      char** argv)
{

	param_t parm;

    // set default parameters
	parm.rcut              =0.;  // Å
	parm.verb              =0;    // silent mode
    parm.save              =0;
    parm.printNlist        =0;
    parm.printMolecules    =0;
    parm.tolerance         =0.;
    parm.inFile            =NULL;

    std::vector<std::string> species;
    std::vector<int> ntypes;

    // container for atomic data
    std::vector<point_t> atoms;

    // allocate matrix for cell vectors
    double** H = new double*[3];
    for (int i=0;i<3;i++) H[i] = new double[3];


    // try to parse the command line arguments passed
    // to the program.  Fails if no input file was found
    if (!get_arguments( nargs,argv,&parm ))
    {
        fprintf(stderr,"ERROR: Command line arguments not parsed properly...\n" );
        return EXIT_FAILURE;
    }
	
    verb = parm.verb;

    // try to rad the atomic data
    int err = read_poscar( parm.inFile,atoms,H,species,ntypes );

    // if something went wrong during input of atomic data then bail
	if (err) return EXIT_FAILURE;

    // construct the nearest neighbor list
    // build a 2x2x2 supercell to impose pbc when searching
    // use a fixed rcut or use sum of covalent radii
    build_nn_list( atoms,H,&parm );

    // get the reduced stoichiometry and form a 
    // filename from it i.e. P4_S8_H    
    parm.prefix =  stoichiometry( species,ntypes );

    std::vector<int> clusterTags;
    std::vector<std::vector<int> > atomClusters = find_clusters(parm.prefix, atoms, species,&parm,clusterTags); 

    if (parm.printNlist)
    {
        print_bonds_nnlist( atoms, species,&parm );
    }
    if (parm.save)
    {
        print_atom_clusters( atoms,H,clusterTags,&parm );
    }
    if (parm.printMolecules)
    {
        print_molecules( atomClusters,atoms,&parm );
    }

    if (4<verb)
        print_atoms_neighbors( atoms,H );

    fprintf(stdout,"Done....: Found %lu clusters\n",atomClusters.size() );
    

    return EXIT_SUCCESS;

}











