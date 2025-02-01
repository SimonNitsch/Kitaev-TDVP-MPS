#include <map>
#include <string>
#include <stdexcept>
#include <ostream>
#pragma once

namespace TDVP_MPS{

class Hamiltonian{

    private:
    std::map<std::string, double> H;

    public:
    Hamiltonian(){
        H["K"] = 0;
        H["Kx"] = 0;
        H["Ky"] = 0;
        H["Kz"] = 0;
        H["J"] = 0;
        H["hx"] = 0;
        H["hy"] = 0;
        H["hz"] = 0;
        H["Gamma"] = 0;
        H["GammaQ"] = 0;
    }

    double get(std::string x){
        if (H.count(x)){
            return H[x];
        }
        else{
            throw std::invalid_argument("Input Argument is not valid");
        }
    }

    void set(std::string x, double n){
        if (H.count(x)){
            H[x] = n;

            if (x == "K"){
                H[x+"x"] = n;
                H[x+"y"] = n;
                H[x+"z"] = n;
            }
        }
        else{
            throw std::invalid_argument("Input Argument is not valid");
        }
    }

    void print(){
        std::cout << "Kx: " << H["Kx"] << "\n";
        std::cout << "Ky: " << H["Ky"] << "\n";
        std::cout << "Kz: " << H["Kz"] << "\n";
        std::cout << "J: " << H["J"] << "\n";
        std::cout << "hx: " << H["hx"] << "\n";
        std::cout << "hy: " << H["hy"] << "\n";
        std::cout << "hz: " << H["hz"] << "\n";
        std::cout << "Gamma: " << H["Gamma"] << "\n";
        std::cout << "GammaQ: " << H["GammaQ"] << "\n" << std::flush;
    }

    std::vector<double> total_info(){
        std::vector<double> filev;
        filev = {H["Kx"], H["Ky"], H["Kz"], H["J"], H["hx"], H["hy"], H["hz"], H["Gamma"], H["GammaQ"]};
        return filev;
    }   


};


}


