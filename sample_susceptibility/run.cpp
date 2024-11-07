#include "Kitaev-TDVP-MPS/main.h"

using namespace TDVP_MPS;

int main(){
    Hamiltonian H_Details;
    H_Details.set("K",1);
    H_Details.set("hz",0.1);

    int LX = 2;
    int LY = 6;

    int auxiliaries = 2;
    int secondary_auxiliaries = 0;

    int DoubleSpin = 1;

    std::string foldername = "Result";


    int max_sites = 256; // Maximum Bond Dimension
    int min_sites = 32; // Initial Bond Dimension of the random MPS


    std::vector<int> timesteps = {50,80,90,100,300}; 
    std::vector<double> intervals = {3,5,12,15,75};

    double susdiff = 0.005;

    /*
    The magnetic Susceptibility is calculated by computing the numerical derivation of the magnetization

    chi = [M(h + Δh) - M(h)] / Δh

    susdiff clarifies the Δh of this calculation, the direction(s) are specified later
    This also means, that an additional imaginary time evolution has to be performed on all random MPS

    susdiff defaults to 0. In this case, only one imaginary time evolution is done per random MPS
    */

    std::string directions = "xz";
    /*
    Specifies the direction(s) of the susceptibilities to be computed
    Here, the magnetic Susceptibilities in x-direction and in z-direction are calculated
    Defaults to "xyz" (Computing for all directions)

    The string has to contain the direction(s), so "x" would only compute the Susceptibility in x-direction
    For all calculations, susdiff is used as Δh
    */
    


    int Evols = 5;

    auto Model = Kitaev_Model(LX,LY,H_Details,DoubleSpin,auxiliaries,secondary_auxiliaries,"HoneycombPeriodic");
    
    Model.TPQ_MPS(timesteps,intervals,Evols,max_sites,min_sites,"TwoSite",susdiff,directions);
    /*
    "TwoSite" specifies the use of 2TDVP (recommended)
    This argument defaults to "TwoSite"

    Alternatively, will "OneTDVP" you can use 1TDVP
    */

    Model.Save(foldername);



}


