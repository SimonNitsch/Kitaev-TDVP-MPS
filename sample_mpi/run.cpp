#include <mpi.h>
#include "Kitaev-TDVP-MPS/main.h"

int main(int argc, char *argv[]){
    MPI_Init(&argc,&argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

    /*
    Getting the world_rank is important for saving the data
    as the data from each node is saved separately
    */

    TDVP_MPS::Hamiltonian H_Details;
    H_Details.set("K",1);
    H_Details.set("hz",0.1);

    int LX = 2;
    int LY = 6;

    int auxiliaries = 2;
    int secondary_auxiliaries = 0;

    int DoubleSpin = 1;
    std::string foldername = "Result";


    std::vector<int> timesteps = {50,80,90,100,300}; 
    std::vector<double> intervals = {3,5,12,15,75};

    int max_sites = 256;
    int min_sites = 32;

    int Evols = 4;

    /*
    Evolutions are parallelized with OpenMP, so each node runs 4 threads
    The total number of evolutions here is 4*nodes
    */


    auto Model = TDVP_MPS::Kitaev_Model(LX,LY,H_Details,DoubleSpin,auxiliaries,secondary_auxiliaries,"HoneycombPeriodic");
    
    Model.TPQ_MPS(timesteps,intervals,Evols,max_sites,min_sites);

    Model.Save(foldername,world_rank);
    /*
    It is important to give world_rank as second argument, because otherwise, the folders from each node will overwrite each other
    This way, each node attaches its world_rank number to the foldername
    Also, the unaveraged data is also saved in the _raw.txt files

    The _raw.txt files can be put together and averaged with the following command:
    python3 /path/to/Kitaev-TDVP-MPS/zip.py foldername nodes

    The final averages can be found in the folder foldername

    So if this program is run on 5 nodes, the data from the nodes will be saved in the folders Result0, Result1, Result2, Result3, Result4
    The command to put them together and calculate the average is:
    python3 /path/to/Kitaev-TDVP-MPS/zip.py Result 5

    This will create the folder Result and save the final averages into it
    */

    MPI_Finalize();
}




