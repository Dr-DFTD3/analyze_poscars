#include "io.h"
#include "lattice.h"


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
             std::vector<int>&         ntypes )
{

    char message[256];
    const char *subroutine = "Vasp::read_poscar()";

    std::ifstream data(inFile);
    std::string line,field,comLine;
    std::istringstream iss;
    std::string arg1,arg2,arg3,arg4,arg5;


    if (!data.is_open()) 
    {
		fprintf(stderr,"could not open [ %s ] for reading!\n",inFile);
		return 1;
    }

    std::getline(data,line);
    iss.clear();iss.str(line);
    iss >> comLine;


    std::getline(data,line);
    iss.clear();iss.str(line);
    iss >> arg1;
    double scale = str2num<double>( arg1 );

    std::getline(data,line);
    iss.clear();iss.str(line);
    iss >> arg1 >> arg2 >> arg3;
    H[0][0] = str2num<double>( arg1 ) * scale;
    H[0][1] = str2num<double>( arg2 ) * scale;
    H[0][2] = str2num<double>( arg3 ) * scale;

    std::getline(data,line);
    iss.clear();iss.str(line);
    iss >> arg1 >> arg2 >> arg3;
    H[1][0] = str2num<double>( arg1 ) * scale;
    H[1][1] = str2num<double>( arg2 ) * scale;
    H[1][2] = str2num<double>( arg3 ) * scale;

    std::getline(data,line);
    iss.clear();iss.str(line);
    iss >> arg1 >> arg2 >> arg3;
    H[2][0] = str2num<double>( arg1 ) * scale;
    H[2][1] = str2num<double>( arg2 ) * scale;
    H[2][2] = str2num<double>( arg3 ) * scale;

    std::vector<std::string> tmp;
    std::string word;

    // sixth line
    std::getline(data,line);
    iss.clear();iss.str(line);
    iss >> arg1;

    char* tmpc = str2char(arg1);

    if (atoi(strtok (tmpc," "))) 
    {
        // this is VASP < 5.x.x
        char msg[MAX_STRING_LENGTH];
        
        fprintf(stderr,"Vasp::read_data() -> WARNING! no species information found...\n" );
        fprintf(stderr,"                     using \"C\" instead\n" );
        
        iss.clear();iss.str(line);
        // grab the count of each type
        while (iss >> word) ntypes.push_back( str2num<int>( word ) );

        for (int i=0;i<ntypes.size();i++) 
        {
            std::string s = "Type ";
            s += num2str( i+1 );
            species.push_back( s );
        }

    } 
    else 
    {
        // this is VASP > 5.x.x
        iss.clear();iss.str(line);
        // grab the element of each type
        while (iss >> word) species.push_back( word );

        word.clear();

        // grab the count of each type
        std::getline(data,line);
        iss.clear();iss.str(line);
        while (iss >> word) ntypes.push_back( str2num<int>( word ) );

    }

    int gtag =1;

    int idx =0;

    int scaled =0;

    int id;

    while (std::getline(data,line)) 
    {
    	iss.clear();iss.str(line);
	    iss >> arg1;

	    if ( strstr( arg1.c_str(),"art" ) ) 
        {
	    	for (int t=0;t<ntypes.size();t++) 
            {
                id =1;
	    		for (int i=0;i<ntypes[t];i++) 
                {
	    			std::getline(data,line);
				    iss.clear();iss.str(line);
				    iss >> arg1 >> arg2 >> arg3;
                    point_t pt;
                    pt.x = scale*str2num<double>( arg1 );
                    pt.y = scale*str2num<double>( arg2 );
                    pt.z = scale*str2num<double>( arg3 );
                    pt.id = id;
                    pt.gid = gtag;
                    pt.tid = t+1;
                    pt.symb = species[t];
                    atoms.push_back( pt );
                    idx++;
				    gtag++;
                    id++;
	    		}
	    	}
	    } 
        else if ( strstr( arg1.c_str(),"irec" ) ) 
        {
	    	scaled =1;
	    	for (int t=0;t<ntypes.size();t++) 
            {
                id =1;
	    		for (int i=0;i<ntypes[t];i++) 
                {
	    			std::getline(data,line);
				    iss.clear();iss.str(line);
				    iss >> arg1 >> arg2 >> arg3;
                    point_t pt;
                    pt.xs = str2num<double>( arg1 );
                    pt.ys = str2num<double>( arg2 );
                    pt.zs = str2num<double>( arg3 );
                    pt.id = id;
                    pt.gid = gtag;
                    pt.tid = t+1;
                    pt.symb = species[t];
                    atoms.push_back( pt );
				    gtag++;
                    idx++;
                    id++;
	    		}
	    	}
	    }
    }

    if (scaled) cart_coords( atoms,H );

    if (!scaled) reduce_coords( atoms,H );

    return 0;
}


/*------------------------------------------------------------
    Write a single POSCAR/CONTCAR file. 
--------------------------------------------------------------*/
int write_poscar( char*                     name,
                  int                       ix,
                  std::vector<point_t>&     atoms,
                  double**                  H, 
                  std::vector<std::string>& species, 
                  std::vector<int>&         ntypes)
{
    
    FILE* fptr;

    if (name==NULL) fptr = stdout;
    else fptr = fopen( name,"w" );

    fprintf(fptr,"cluster: %i\n",ix );   

    fprintf(fptr,"1.0000000\n" );    

    fprintf(fptr,"%12f %12f %12f\n",H[0][0],H[0][1],H[0][2] );
    fprintf(fptr,"%12f %12f %12f\n",H[1][0],H[1][1],H[1][2] );
    fprintf(fptr,"%12f %12f %12f\n",H[2][0],H[2][1],H[2][2] );
    for (int i=0;i<species.size();i++)
    {
        fprintf(fptr,"%s ",species[i].c_str() );
    }

    fprintf(fptr,"\n" );

    for (int n=0;n<ntypes.size();n++) {
        fprintf(fptr,"%i ",ntypes[n] );
    }
    fprintf(fptr,"\n" );


    // sort the atoms to be written by id as
    // required by VASP conventions
    sort_ids( atoms );

    fprintf(fptr,"Cartesian\n" );
    for (int n=0;n<atoms.size();n++) {
        fprintf(fptr,"%12f %12f %12f\n",atoms[n].x,atoms[n].y,atoms[n].z );
    }


    if (fptr!=stdout) fclose( fptr );

    return 1;

}

/*------------------------------------------------------------
    Scan list of "atomClusters" and write the coordinates
    of each atom in a cluster to a single POSCAR file
--------------------------------------------------------------*/
void
print_atom_clusters(std::vector<point_t>&     atoms,
                    double**                  H,
                    std::vector<int>          nc,
                    param_t*                  param )
{
    std::vector<int> temp = nc;
    std::vector<int> uniq;

    // get only unique cluster ids
    std::sort( temp.begin(),temp.end() );
    std::unique_copy(temp.begin(), temp.end(),std::back_inserter(uniq));


    bool found = false;
    bool empty = true;
    int  id    = 0;

    // loop over each unique id and check each atom
    for (int i=0;i<uniq.size();i++)
    {
        // container for cluster data
        std::vector<point_t> cluster;
        std::vector<std::string> species;
        std::vector<int> ids;
        int tag = uniq[i];

        empty = true;
        for(int n=0;n<atoms.size();n++)
        {   
            if(atoms[n].clusterID==tag)
            {
                empty = false;
                ids.push_back( n );
            }
        }

        if(!empty)
        {
            for (int p=0;p<ids.size();p++)
            {
                point_t atomi;
                id = ids[p];

                atomi.x = atoms[id].x;
                atomi.y = atoms[id].y;
                atomi.z = atoms[id].z;
                atomi.symb = atoms[id].symb;
                atomi.tid = atoms[id].tid;
                atomi.gid = atoms[id].gid;
                cluster.push_back( atomi );
             
                found = false;
                for (int s=0;s<species.size();s++)
                {
                    if (species[s]==atomi.symb) 
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    species.push_back( atomi.symb );
                }
            }

            // now the cluster is "full", identify
            // unique species and colunt them.
            std::vector<int> nt( species.size() );
            for(int t=0;t<species.size();t++)
            {
                nt[t] = 0;
                for(int c=0;c<cluster.size();c++)
                {
                    if (species[t]==cluster[c].symb) nt[t]++;
                }
            }
            char name[50];

            sprintf( name, "%s-mol-%i.vasp",param->prefix,i+1);

            if (param->save)
                write_poscar( name,i+1,cluster,H,species,nt );
            
        } // if cluster is not empty       
    } // loop over possible clusters
}

/*------------------------------------------------------------
    Write each atom and its nearest neighbors to a
    POSCAR file
--------------------------------------------------------------*/
void
print_atoms_neighbors(std::vector<point_t>&     atoms,
                      double**                  H)
{
    for (int i=0;i<atoms.size();i++)
    {
        // container for cluster data
        std::vector<point_t> cluster;
        std::vector<std::string> species;
        
        point_t atomi;
        atomi.x = atoms[i].x;
        atomi.y = atoms[i].y;
        atomi.z = atoms[i].z;
        atomi.symb = atoms[i].symb;
        atomi.tid = atoms[i].tid;
        atomi.gid = atoms[i].gid;
        cluster.push_back( atomi );
        species.push_back( atomi.symb );

        bool found = false;
        for (int j=0;j<atoms[i].naborIds.size();j++) 
        {
            point_t atomj;
            int nj = atoms[i].naborIds[j];
            atomj.x = atoms[nj].x;
            atomj.y = atoms[nj].y;
            atomj.z = atoms[nj].z;
            atomj.symb = atoms[nj].symb;
            atomj.tid = atoms[nj].tid;
            atomj.gid = atoms[nj].gid;
            cluster.push_back( atomj );

            for (int s=0;s<species.size();s++)
            {
                if (species[s]==atomj.symb) 
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                species.push_back( atomj.symb );
            }
        }

        std::vector<int> nt( species.size() );
        for(int n=0;n<species.size();n++)
        {
            nt[n] = 0;
            for(int i=0;i<cluster.size();i++)
            {
                if (species[n]==cluster[i].symb) nt[n]++;
            }
        }
        char name[50];

        sprintf( name, "atom-%i-neighs.vasp",i+1);

        write_poscar( name,atoms[i].gid,cluster,H,species,nt );
    }
}

/*------------------------------------------------------------
    Scan list of "atomClusters" and identify what molecules 
    are present.
--------------------------------------------------------------*/
void 
print_molecules(std::vector<std::vector<int> > atomClusters,
               std::vector<point_t>&           atoms,
               param_t*                        p)
{

    int v = p->verb;
    int no =0;
    int nh = 0;
    int nc = 0;
    int nn = 0;

    int n_azide = 0;
    int n_n2    = 0;
    int n_h2    = 0;
    int n_h2o   = 0;
    int n_oh    = 0;
    int n_no    = 0;
    int n_no2   = 0;
    int n_co    = 0;
    int n_co2   = 0;
    int n_ch2   = 0;
    int n_ch4   = 0;
    int n_unk   = 0;
    for (int i=0;i<atomClusters.size();i++)
    {
        nn = nh = nc = no = 0;
        int nmembers = atomClusters[i].size();
        for (int n=0;n<nmembers;n++)
        {
            int id = atomClusters[i][n];
            if (atoms[id].symb == "N") nn++;
            else if (atoms[id].symb == "H") nh++;
            else if (atoms[id].symb == "O") no++;
            else if (atoms[id].symb == "C") nc++;
        }

        if (nn==3 && nh==0 && nc==0 && no==0) 
        {
            n_azide++;
            if(2<v) fprintf(stderr,"FOUND AN AZIDE\n" );
        }
        else if (nn==2 && nh==0 && nc==0 && no==0) 
        {
            n_n2++;
            if(2<v) fprintf(stderr,"FOUND AN N2\n" );
        }
        else if (nn==0 && nh==2 && nc==0 && no==0) 
        {
            n_h2++;
            if(2<v) fprintf(stderr,"FOUND AN H2\n" );
        }
        else if (nn==0 && nh==1 && nc==0 && no==1) 
        {
            n_oh++;
            if(2<v) fprintf(stderr,"FOUND AN OH\n" );
        }
        else if (nn==0 && nh==2 && nc==0 && no==1) 
        {
            n_h2o++;
            if(2<v) fprintf(stderr,"FOUND AN H2O\n" );
        }
        else if (nn==1 && nh==0 && nc==0 && no==1) 
        {
            n_no++;
            if(2<v) fprintf(stderr,"FOUND AN NO\n" );
        }
        else if (nn==1 && nh==0 && nc==0 && no==2) 
        {
            n_no2++;
            if(2<v) fprintf(stderr,"FOUND AN NO2\n" );
        }
        else if (nn==0 && nh==0 && nc==1 && no==1) 
        {
            n_co++;
            if(2<v) fprintf(stderr,"FOUND AN CO\n" );
        }
        else if (nn==0 && nh==0 && nc==1 && no==2) 
        {
            n_co2++;
            if(2<v) fprintf(stderr,"FOUND AN CO2\n" );
        }
        else if (nn==0 && nh==2 && nc==1 && no==0) 
        {
            n_ch2++;
            if(2<v) fprintf(stderr,"FOUND AN CH2\n" );
        }
        else if (nn==0 && nh==4 && nc==1 && no==0) 
        {
            n_ch4++;
            if(2<v) fprintf(stderr,"FOUND AN CH4\n" );
        }
        else
        {
            n_unk++;
            if(2<v) fprintf(stderr,"FOUND AN UNKNOWN MOLECULE\n" );
        }
    }

    fprintf(stdout,"\n** Molecule Info: %s\n",p->prefix );
    if(n_azide) fprintf(stdout,"  %3i \"AZIDE\" molecules\n",n_azide );
    if (n_n2) fprintf(stdout,"  %3i \"N2\" molecules\n",n_n2 );
    if (n_h2) fprintf(stdout,"  %3i \"H2\" molecules\n",n_h2 );
    if (n_h2o) fprintf(stdout,"  %3i \"H2O\" molecules\n",n_h2o );
    if (n_oh) fprintf(stdout,"  %3i \"OH\" molecules\n",n_oh );
    if (n_no) fprintf(stdout,"  %3i \"NO\" molecules\n",n_no );
    if (n_no2) fprintf(stdout,"  %3i \"NO2\" molecules\n",n_no2 );
    if (n_co) fprintf(stdout,"  %3i \"CO\" molecules\n",n_co );
    if (n_co2) fprintf(stdout,"  %3i \"CO2\" molecules\n",n_co2 );
    if (n_ch2) fprintf(stdout,"  %3i \"CH2\" molecules\n",n_ch2 );
    if (n_ch4) fprintf(stdout,"  %3i \"CH4\" molecules\n",n_ch4 );
    if (n_unk) fprintf(stdout,"  %3i \"Unknown\" molecules\n",n_unk );

}


/*------------------------------------------------------------
    Print the neighbor list and bond lengths for each
    atom in the system  
--------------------------------------------------------------*/
void 
print_bonds_nnlist( std::vector<point_t>&     atoms, 
                    std::vector<std::string>& species,
                    param_t*                  p )
{
    char fname[20];

    sprintf( fname, "%s.bonds", p->prefix );

    FILE *fptr;

    if (p->save) fptr = fopen( fname,"a" );
    else fptr = stdout;

    fprintf(fptr,"\n** Nearest Neighbor Info: %s\n",p->prefix );

    for (int i=0;i<atoms.size();i++) {

        fprintf(fptr,"  %s%i atom(%i) has the following neighbors:\n",atoms[i].symb.c_str(),atoms[i].id,atoms[i].gid );

        for (int j=0;j<atoms[i].naborIds.size();j++) 
        {
            int nj = atoms[i].naborIds[j];
            std::string sj = atoms[nj].symb;
            int ij = atoms[nj].id;
            fprintf(fptr,"     %s%i - %s%i = %f\n",atoms[i].symb.c_str(),i+1,sj.c_str(),nj+1,atoms[i].rnn[j] );            
        }
    }
    if ( fptr != stdout ) fclose( fptr );
}


