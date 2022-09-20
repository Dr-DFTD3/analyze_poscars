#include "common.h"

PeriodicTable ptab = 
    /* Elements list */ //std::string
    {"EM",
     "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
     "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", 
     "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
     "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", 
     "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
     "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
     "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
     "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
     "Ti", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
     "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
     "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    /* Atomic masses */                 
    {0.,
     1.0079, 4.0026, 6.941, 9.0122, 10.811, 12.011, 14.007, 15.999, 18.998, 20.180, 
     22.990, 24.305, 26.982, 28.086, 30.974, 32.065, 35.453, 39.948, 39.098, 40.078, 
     44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.390,
     69.723, 72.610, 74.992, 78.960, 79.904, 83.800, 85.468, 87.620, 88.906, 91.224,
     92.906, 95.940, 98.000, 101.070, 102.910, 106.420, 107.87, 112.41, 114.82, 118.71,
     121.76, 127.60, 126.90, 131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24,
     145.00, 150.36, 151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.04,
     174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59,
     204.38, 207.20, 208.98, 209.00, 210.00, 222.00, 223.00, 226.00, 227.00, 232.04,
     231.04, 238.03, 237.00, 244.00, 243.00, 247.00, 247.00, 251.00, 252.00, 257.00,
     258.00, 259.00, 262.00, 261.00, 262.00, 266.00, 264.00, 269.00, 268.00
    },
    /* Atomic radii */
    {0.,
     0.37, 0.93, 1.23, 0.90, 0.82, 0.77, 0.75, 0.73, 0.72, 0.71, 
     1.54, 1.36, 1.18, 1.11, 1.06, 1.02, 0.99, 0.98, 2.03, 1.74, 
     1.44, 1.32, 1.22, 1.18, 1.17, 1.17, 1.16, 1.15, 1.17, 1.25,
     1.26, 1.22, 1.20, 1.16, 1.14, 1.12, 2., 1.91, 1.62, 1.45,
     1.34, 1.30, 1.27, 1.25, 1.25, 1.28, 1.34, 1.48, 1.44, 1.41,
     1.41, 1.36, 1.33, 1.31, 2.35, 1.98, 1.69, 1.65, 1.65, 1.64,
     1.63, 1.62, 1.85, 1.61, 1.59, 1.59, 1.58, 1.57, 1.56, 1.74,
     1.56, 1.44, 1.34, 1.30, 1.28, 1.26, 1.27, 1.30, 1.34, 1.49,
     1.48, 1.47, 1.46, 1.46, 1.45, 1.65, 1.42, -1, -1, -1,
     -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     -1, -1, -1, -1, -1, -1, -1, -1, -1
    }
    };


/*------------------------------------------------------------
    Parse the command line arguments provided.
    Return 1 if successful, 0 if not.
--------------------------------------------------------------*/
int
get_arguments(int      nargs,
              char**   argv, 
              param_t* p)
{
    bool foundFile = false;
    if (nargs==1) show_help();
    for (int ii = 1; ii < nargs; ii++) 
    {
        if (argv[ii][0]=='-' && argv[ii][1]!='-' ) 
        {
            char* option =argv[ii];
            for (int l=1;l<strlen(option);l++) 
            {
                if (option[l]=='h') show_help( );
                else if (option[l]=='g') p->verb =5;
                else if (option[l]=='m') p->printMolecules =1;
                else if (option[l]=='n') p->printNlist =1;
                // else if (option[l]=='p') printAtomsNeighbors =1;
                else if (option[l]=='r') p->rcut = atof( argv[ii+1] );
                else if (option[l]=='s') p->save =1;
                else if (option[l]=='t') p->tolerance = atof( argv[ii+1] );
                else if (option[l]=='u') show_usage( ); 
                else if (option[l]=='v') p->verb++; 
            }
        } 
        else if (argv[ii][0]=='-' && argv[ii][1]=='-') 
        {
            argv[ii]++;argv[ii]++;
            if (strcmp(argv[ii],"help")==0) show_help( );
            else if (strcmp(argv[ii],"molecules")==0) p->printMolecules =1;
            else if (strcmp(argv[ii],"neighbors")==0) p->printNlist =1;
            else if (strcmp(argv[ii],"radius")==0) p->rcut = atof( argv[ii+1] );
            else if (strcmp(argv[ii],"save")==0) p->save =1;
            else if (strcmp(argv[ii],"tolerance")==0) p->tolerance = atof( argv[ii+1] );
            else if (strcmp(argv[ii],"debug")==0) p->verb =5;
            else if (strcmp(argv[ii],"usage")==0) show_usage( );
            else if (strcmp(argv[ii],"version")==0) show_version( );
        }
    }
    for (int ii = 1; ii < nargs; ii++) 
    {
        if ( is_regular_file( argv[ii] ) )
        {
            p->inFile = argv[ii];
            foundFile = true;
            break;
        }
    }

    if (!foundFile) return 0;
    else return 1;
}

/*------------------------------------------------------------
    Scan periodic table to get the covalent radius of 
    "element". Return  covalent radius in Å
--------------------------------------------------------------*/
double 
elem_radii(std::string element)
{
    int ii;
    for(ii=1;ii<110;++ii) 
    {
        if (element == ptab.element[ii]) return ptab.radius[ii];
    }

    return -1;
}

/*------------------------------------------------------------
    Define a cutoff for a pair of atoms based on the
    sum of the covalent radius. dr represents a fudge factor
    in determining a bond.
--------------------------------------------------------------*/
double 
sum_covalent_radii( std::string atomi, 
                    std::string atomj, 
                    double      dr)
{
    double cutoff = 0;
    double irad = elem_radii( atomi );
    double jrad = elem_radii( atomj );

    if (irad < 0.)
    {
        fprintf(stderr,"ERROR: No covalent radius information for %s\n",atomi.c_str() );
        return -1;
    }
    if (jrad < 0.)
    {
        fprintf(stderr,"ERROR: No covalent radius information for %s\n",atomj.c_str() );
        return -1;
    }
    else
    {
        cutoff = irad + jrad + dr;
        return cutoff; 
    }
}

/*------------------------------------------------------------
    Sort the list of atoms based on type ID, this is required
    by the format of POSCAR/CONTCAR files
--------------------------------------------------------------*/
void 
sort_ids( std::vector<point_t>& atoms )
{

    // need to pack by ids first, 1 1 1 2 2 2 ...
    std::vector<std::vector<double> > temp;
    std::vector<double> row(5);

    std::vector<std::vector<std::string> > elems;
    std::vector<std::string> row2(2);

    
    for (int i=0;i<atoms.size();i++) {
        row[0] = atoms[i].gid;
        row[1] = atoms[i].tid;
        row[2] = atoms[i].x;
        row[3] = atoms[i].y;
        row[4] = atoms[i].z;

        row2[0] = atoms[i].tid;
        row2[1] = atoms[i].symb;
        temp.push_back( row );
        elems.push_back( row2 );
    }
    

    // sort the temp array to get the atoms in 
    // the correct ordering
    std::sort( temp.begin(),temp.end(),col2_asc() );

    // sorrt the elems array to get the species
    // matched with the atoms
    std::sort( elems.begin(),elems.end(),col2_asc() );

    // if the output corrdinates are cartesian,
    // only copy those coordinates
    for (int i=0;i<atoms.size();i++) {
        atoms[i].gid  = temp[i][0]; 
        atoms[i].tid = temp[i][1];
        atoms[i].x = temp[i][2];
        atoms[i].y = temp[i][3];
        atoms[i].z = temp[i][4];
        atoms[i].symb = elems[i][1];
    }
}

/*------------------------------------------------------------
    Convert a std::string type to a char*
------------------------------------------------------------*/
char* 
str2char(std::string str)
{
    char *arg = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), arg);
    arg[str.size()] = '\0';
    return arg;
}

/*------------------------------------------------------------
    Formulate a name based on the stoichiometry of the system
--------------------------------------------------------------*/
char*
base_name( std::vector<std::string> species, 
           std::vector<int>         ntypes )
{
    std::string name;
    for(int i=0;i<species.size();i++) 
    {
        name += species[i];
        name += num2str( ntypes[i] );
    }
    return str2char( name );
}

char* 
stoichiometry( std::vector<std::string> species, 
               std::vector<int>         ntypes )
{

    int divisor = ntypes[0];
    for (int i=1;i<ntypes.size();i++) divisor = gcd( divisor,ntypes[i] );

    std::string unit;
    int last = species.size()-1;
    for(int i=0;i<species.size();i++) 
    {   
        unit += species[i];
        if ((ntypes[i]/divisor )!=1) unit += num2str( ntypes[i]/divisor );

        if (i!=last) unit += "_";
    }

    return str2char( unit );
}

/*------------------------------------------------------------
    Checks if a given character or string of characters is 
    a regular file which can be opened and read
------------------------------------------------------------*/
int 
is_regular_file(const char *ch)
{
    struct stat ch_status;
    stat(ch, &ch_status);
    return S_ISREG(ch_status.st_mode);
}

void 
show_help()
{
    fprintf(stderr,"   %s v %s\n",SOFT,VERS);
    fprintf(stderr,"  Utility to analyze molecules/clusters in a POSCAR/CONTCAR  \n");
    fprintf(stderr,"  structure file. Periodic boundaries are assumed.\n" );
    fprintf(stderr,"  Currently known molecules: \n");
    fprintf(stderr,"       H2,O2,N2,C2,OH,CO,NO,H2O,NO2,CO2,CH4,CH2,Azide\n\n");
    fprintf(stderr,"options...\n");
    fprintf(stderr,"     -m|--molecules,\n");
    fprintf(stderr,"             analyze clusters found to identify well known molecules, print \n");
    fprintf(stderr,"             them to STDOUT\n");
    fprintf(stderr,"     -n|--neighbors,\n");
    fprintf(stderr,"             print bonding/nearest neighbor info. Use a fixed cutoff \"RADIUS\"\n");
    fprintf(stderr,"             or use covalent radii for pair of atoms\n");
    fprintf(stderr,"     -r|--radius RADIUS,\n");
    fprintf(stderr,"             uniform search radius for assigning nearest neighbors, in units of Å.\n");
    fprintf(stderr,"             default == use sum of covalent radii of two neighboring atoms \n");
    fprintf(stderr,"     -s|--save,\n");
    fprintf(stderr,"             export clusters/molecules found to \"STOICH-mol-%%i.vasp\"\n");
    fprintf(stderr,"             STOICH = stoichiometry of the loaded structure, i.e. P4S8H1\n");
    fprintf(stderr,"             default == do not save\n");
    fprintf(stderr,"     -t|--tolerance TOL,\n");
    fprintf(stderr,"             tolerance for assigning neighbors when using the covalent radius\n");
    fprintf(stderr,"             as the cutoff, in units of Å.\n");
    fprintf(stderr,"             default == 0. \n");
    fprintf(stderr,"     -v|-vv|-vvv,\n");
    fprintf(stderr,"             amount of info printed to \"stderr\"\n");
    fprintf(stderr,"             default = silent\n");
    fprintf(stderr,"                 -v   == startup messages\n");    
    fprintf(stderr,"                 -vv  == + cluster and neighbor information\n"); 
    fprintf(stderr,"                 -vvv == + all intermediate data and subroutine calls\n");
    fprintf(stderr,"extras...\n");
    fprintf(stderr,"     -g|--debug, \n");
    fprintf(stderr,"            extensive info, equivalent to \"-vvvv\"\n");
    fprintf(stderr,"     -h|--help,\n");
    fprintf(stderr,"            this message\n" );
    fprintf(stderr,"     -u|--usage,\n");
    fprintf(stderr,"             brief list of parmaters and arguments\n" );
    fprintf(stderr,"       --version,\n");
    fprintf(stderr,"             display version and author info then exit\n\n");
    fprintf(stderr,"examples....\n");
    fprintf(stderr,"USAGE: %s [OPTIONS] ... [FILE] \n",EXE);
    fprintf(stderr,"  %s POSCAR -m > molecules.nfo \n",EXE);
    fprintf(stderr,"  %s -n CONTCAR -c > cluster-data.nfo \n",EXE);
    fprintf(stderr,"  %s sns2.vasp -vvv --radius 2 --tolerance 0.2 --clusters --save  \n",EXE);
    fprintf(stderr,"  %s -vvncs mos2.contcar \n",EXE);
    exit(EXIT_FAILURE);
}


void 
show_version()
{
    fprintf(stdout,"(MSL coreutils) %s v. %s\n",EXE,VERS );
    fprintf(stdout,"Copyright (C) 2016 Materials Simulation Laboratory\n" );
    fprintf(stdout,"This is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n\n" );
    fprintf(stdout,"Written by: Joseph M. Gonzalez\n");
    exit(EXIT_SUCCESS);
}


void 
show_usage()
{
    fprintf(stdout,"(MSL coreutils 2016) %s v. %s\n",EXE,VERS );
    fprintf(stdout,"USAGE: %s [OPTIONS] FILE \n",EXE );
    fprintf(stdout,"OPTIONS: [-m molecules info]     [-n neighbors info] [-r cutoff]\n" );
    fprintf(stdout,"         [-s save clusters]      [-t tolerance] \n" );
    fprintf(stdout,"         [-v|-vv|-vvv verbosity] [--version]\n" );
    exit(EXIT_SUCCESS);
}
