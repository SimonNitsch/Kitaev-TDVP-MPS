#include "Kitaev-TDVP-MPS/main.h"

using namespace TDVP_MPS;

int main(){
    Hamiltonian H_Details;
    H_Details.set("K",1);
    H_Details.set("hz",0.1);

    /*
    The class Hamiltonian can be modified with .set(), as demonstrated here.
    Coefficients that can be changed, are:
    Kx ... Kitaev Interaction in x direction
    Ky ... Kitaev Interaction in y direction
    Kz ... Kitaev Interaction in z direction
    K ... Changes Kx, Ky and Kz at once
    J ... Heisenberg Interaction
    hx ... Magnetic Field in x direction
    hy ... Magnetic Field in y direction
    hz ... Magnetic Field in z direction
    Gamma ... Gamma Interaction
    GammaQ ... Gamma' Interaction
    */

    int LX = 2;
    int LY = 6;
    /*
    Size of the Lattice in x and y direction
    For the Honeycomb Lattice, LY has to be even
    */

    int auxiliaries = 2;
    int secondary_auxiliaries = 0;

    /*
    Number of Auxiliary Tensors used
    auxiliaries attaches tensors at each end of the MPS
    secondary_auxiliaries inserts tensors between first & second, and second-to-last & last column
    */

    int DoubleSpin = 1;
    /*
    Twice the Spin of the System
    For Spin1 System use 2, for Spin3/2 use 3
    */

    std::string foldername = "Result";
    /*
    Name of the folder where all the computed quantities will be saved to
    (The folder will be automatically generated)
    */

    int max_sites = 256; // Maximum Bond Dimension
    int min_sites = 32; // Initial Bond Dimension of the random MPS


    std::vector<int> timesteps = {50,80,90,100,300}; 
    std::vector<double> intervals = {3,5,12,15,75};
    /*
    timesteps and intervals tell the code about number and length of the time evolutions

    This modification means: 
    From beta = 0 to beta = 3, there are 50 steps
    From beta = 3 to beta = 8, there are 80 steps
    From beta = 8 to beta = 20, there are 90 steps
    From beta = 20 to beta = 35, there are 100 steps
    From beta = 35 to beta = 110, there are 300 steps

    Please note, that the imaginary time steps made have the length beta/2
    Also, the timesteps and intervals vectors have to have the same length

    There are probably some performance gains to be made by playing around with these numbers
    */


    int Evols = 5;
    /*
    Number of total Evolutions made

    One Evolution means
    One random MPS is generated and then evolved in imaginary time according to timesteps and intervals
    
    So here, 5 random MPS are generated and evolved in imaginary time
    */


    auto Model = Kitaev_Model(LX,LY,H_Details,DoubleSpin,auxiliaries,secondary_auxiliaries,"HoneycombPeriodic");
    /*
    Generates the Model
    
    For the Lattice Modification:
    Honeycomb ... Honeycomb Lattice with semi-open boundary conditions (efficient for large lattices)
    HoneycombPeriodic ... Honeycomb Lattice with periodic boundary conditions
    Triangular ... Triangular Lattice with semi-open boundary conditions (efficient for large lattices)
    TriangularPeriodic ... Triangular Lattice with periodic boundary conditions
    */
    
    Model.TPQ_MPS(timesteps,intervals,Evols,max_sites,min_sites);
    // Imaginary Time Evolutions with TPQ-MPS (Alternatively, the command tanTRG() can also be used)

    Model.Save(foldername);
    // Generates the folder foldername and saves the data obtained from TPQ-MPS() into it as .txt files
    // Additionally, the file x_data.txt is saved, which contains all temperature points (for plotting purposes)



}


