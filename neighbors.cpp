#include "neighbors.h"

/*------------------------------------------------------------
    Construct a nearest neighbor list for each atom in the 
    system. Use an n^2 recursive search for nearest neighbors 
    within a fixed cutoff RADIUS, or use covalent radii for 
    pair of atoms. Returns nothing, update the neighborlist 
    of each atom and the translation vectors connecting each
    atom with each of its neighbors.
--------------------------------------------------------------*/
void 
build_nn_list( std::vector<point_t>& atoms,
			   double**              H,
			   param_t*              param )
{

    int found =0;
    double rc =0;
    bool useCovalent = true;
    if(param->rcut!=0) 
    {
    	useCovalent = false;
    	rc = param->rcut;
    }

    // loop over all periodic images
    // param->minBond =10;
    param->nbonds = 0;
    for(int iz=-1;iz<=1;iz++) {
        for(int iy=-1;iy<=1;iy++) {
            for(int ix=-1;ix<=1;ix++) {
                double shiftX = ix*H[0][0] + iy*H[1][0] + iz*H[2][0];
                double shiftY = ix*H[0][1] + iy*H[1][1] + iz*H[2][1];
                double shiftZ = ix*H[0][2] + iy*H[1][2] + iz*H[2][2];
                for (int i=0;i<atoms.size();i++) {
                	register double xi = atoms[i].x;
                	register double yi = atoms[i].y;
                	register double zi = atoms[i].z;
                	std::string     si = atoms[i].symb;
                    for (int j=0;j<atoms.size();j++) {
                        point_t atomj = atoms[j];
	                	
                        // compute the position of real and image atoms for
                        // neighbor consideration
                        register double xj = atoms[j].x + shiftX;
	                	register double yj = atoms[j].y + shiftY;
	                	register double zj = atoms[j].z + shiftZ;
	                	std::string     sj = atoms[j].symb;

                        // if the user did not provide
                        // a uniform cutoff, then use the
                        // sum of covalent radii: (i) + (j)
                        if (useCovalent)
                        {
                        	rc = sum_covalent_radii( si,sj,param->tolerance );
                        }
                        
                        double bond =  sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj) );
                        // do this calculation here, dont call func
                        // register double bond = distance( atoms[i],atomj );
                        if (bond!=0 && bond<=rc) {
                            // bond = sqrt(bond);

                            std::vector<int> tvec;
                            tvec.push_back(shiftX);
                            tvec.push_back(shiftY);
                            tvec.push_back(shiftZ);
                            atoms[i].tvecs.push_back( tvec );
                            atoms[i].rnn.push_back( bond );
                            found =1;
                            atoms[i].naborIds.push_back( j );
                            // min_pair( param,bond,atoms[i].id,si,atoms[j].id,sj );
                            param->nbonds++;
                            if (2<param->verb) {
                                fprintf(stderr," atomi %3s(%3i)  atomj %3s(%3i)  (%5.3f,%5.3f,%5.3f) d = %5.3f\n",
                                             si.c_str(),atoms[i].gid,sj.c_str(),atoms[j].gid,xj,yj,zj,bond );
                            } // print info
                        } // check distance  
                    } // j
                } // i
            } // ix
        } // iy
    }// iz

}


/*------------------------------------------------------------
    Cycle through each atoms' nearest neighbor list to 
    identify clusters of atoms, i.e. molecules.
    Assume each atom is in a cluster by itself, recursivley 
    update which cluster each atom  belongs to. Return 
    an array containing the member ids of each atom belonging
    to a cluster.
--------------------------------------------------------------*/
std::vector<std::vector<int> >
find_clusters( char*                     stdName, 
               std::vector<point_t>&     atoms, 
               std::vector<std::string>& species,
               param_t*                  param,
               std::vector<int>&         clusterTags )
{
    std::vector<std::vector<int> > atomClusters;
    for (int i=0;i<atoms.size();i++)
    {
        clusterTags.push_back(0);
        atoms[i].clusterID = atoms[i].gid;
    }

    bool finished;

    double xi,yi,zi;
    double xj,yj,zj;

    double rsq,rc=0;
    bool useCovalent = true;

    if(param->rcut!=0) 
    {
        useCovalent = false;
        rc = param->rcut;
    }
    else
    {
        // find the maximum allowed bond for all possible
        // species in the molecule. 
        for( int i=0;i<species.size();i++ )
        {
            for (int j=i;j<species.size();j++)
            {
                double rij = sum_covalent_radii( species[i],species[j],param->tolerance );
                rc = max(rc,rij);
                if(2<param->verb)
                    fprintf(stderr,"%s-%s ",species[i].c_str(),species[j].c_str() );
            }
        } 
    }
    if(2<param->verb)
        fprintf(stderr,"\n" );

    if(param->verb>2)
    {
        if(useCovalent)
        {
            fprintf(stderr,"CLUSTER: MAX COVALENT PAIR: %f\n",rc);
        }
        else
        {
            fprintf(stderr,"CLUSTER: MAX CUTOFF:         %f\n",rc);
        }
    }

    int tag = 100;
    while(1)
    {
        while(1)
        {
            finished = true;
            for (int i=0;i<atoms.size();i++)
            {
                xi = atoms[i].x;
                yi = atoms[i].y;
                zi = atoms[i].z;
                std::string     si = atoms[i].symb;
                register int itag = atoms[i].clusterID;
                for (int j=0;j<atoms[i].naborIds.size();j++) 
                {
                    int nj = atoms[i].naborIds[j];
                    register int jtag = atoms[nj].clusterID;
                    if (itag == jtag) continue;

                    double tx = atoms[i].tvecs[j][0];
                    double ty = atoms[i].tvecs[j][1];
                    double tz = atoms[i].tvecs[j][2];
                    std::string     sj = atoms[nj].symb;

                    xj = atoms[nj].x + tx;
                    yj = atoms[nj].y + ty;
                    zj = atoms[nj].z + tz;

                    double dx = xi - xj;
                    double dy = yi - yj;
                    double dz = zi - zj;

                    rsq = sqrt( dx*dx + dy*dy + dz*dz );
                    int pi = atoms[i].gid;
                    int pj = atoms[nj].gid;
                    
                    if (rsq <= rc)
                    {
                        if (2<param->verb)
                            fprintf(stderr,"%s%i-%s%i  rs=%f  rc=%f \n",si.c_str(),pi,sj.c_str(),pj,rsq,rc);
                        itag = jtag = min(itag,jtag);
                        atoms[i].clusterID = itag;
                        atoms[nj].clusterID = jtag;
                        finished = false;
                    }
                }
            }
            if (finished) break;
        }
        if (finished) break;
    }

    for (int i=0;i<atoms.size();i++) clusterTags[i] = atoms[i].clusterID;

    std::vector<int> temp = clusterTags;
    std::vector<int> uniq;
    bool empty;

    // get only unique cluster ids
    std::sort( temp.begin(),temp.end() );
    std::unique_copy(temp.begin(), temp.end(),std::back_inserter(uniq));

    // loop over each unique id and check each atom
    if (param->verb)
    {
        fprintf(stdout,"\n** Cluster Info: %s\n",param->prefix );
        fprintf(stdout,"    ID   Members\n");
    }
    for (int i=0;i<uniq.size();i++)
    {
        std::vector<int> members;
        int tag = uniq[i];

        empty = true;
        for(int n=0;n<atoms.size();n++)
        {   
            if(atoms[n].clusterID==tag)
            {
                empty = false;
                members.push_back( n );
            }
        }

        if(!empty)
        {
            atomClusters.push_back( members );
            int id;
            if (param->verb)
            {
                fprintf(stdout,"    %-3i: ",tag);
                for (int p=0;p<members.size();p++)
                {
                    point_t atomi;
                    id = members[p];
                    fprintf(stdout,"%s%i ",atoms[id].symb.c_str(),id );
                }
                fprintf(stdout,"\n");
            }
        }
    }
    return atomClusters;
}





