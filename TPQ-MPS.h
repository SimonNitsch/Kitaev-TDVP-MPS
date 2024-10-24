#include <iostream>
#include <array>
#include <vector>
#include <chrono>
#include <ios>
#include <type_traits>
#include <cmath>
#include <cstdio>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <exception>
#include <filesystem>
#include <climits>
#include <cstdlib>
#include <future>
#include <mutex>
#include "itensor/all.h"
#include "TDVP/tdvp.h"
#include "TDVP/basisextension.h"
#include "customspin_nitsch.h"
#include "Hamiltonian.h"
#pragma once


std::vector<double> operator-(double x, std::vector<double>& v){
    std::vector<double> vx;
    vx.reserve(v.size());
    
    for (auto& i : v){
        vx.push_back(x-i);
    }
    return vx;
}

std::vector<double> operator-(std::vector<double>& v1, std::vector<double>& v2){
    std::vector<double> v;
    v.reserve(v1.size());
    
    for (int i = 0; i != v1.size(); i++){
        v.push_back(v1[i] - v2[i]);
    }
    return v;
}

std::vector<double> operator/(std::vector<double>& v, double x){
    std::vector<double> vx;
    vx.reserve(v.size());

    for (auto& i : v){
        vx.push_back(i/x);
    }
    return vx;
}



using namespace itensor;
namespace TDVP_MPS{



class Kitaev_Model{

    public:
    AutoMPO ampo;
    int Lattice_Type; // 1 for Honeycomb, 2 for Triangular, 3 for Periodic Honeycomb
    bool CalcDMRG = false;
    bool CalcTDVP = false;

    public:
    MPO H0, H2, H0x, H0y, H0z;
    std::array<MPO,3> M;
    std::array<MPO,3> M2;
    Hamiltonian H_Details;
    CustomSpinNitsch sites;
    CustomSpin sitestan;
    MPO H_flux;
    MPS Htan;
    int LX, LY, aux, sec_aux;

    private:
    double GSE;
    std::vector<std::array<double,2>> E;
    std::vector<std::array<double,2>> Cv;
    std::vector<std::array<double,2>> S;
    std::vector<std::array<double,2>> W;
    std::vector<std::array<double,2>> Mx;
    std::vector<std::array<double,2>> My;
    std::vector<std::array<double,2>> Mz;
    std::vector<std::array<double,2>> Mx2;
    std::vector<std::array<double,2>> My2;
    std::vector<std::array<double,2>> Mz2;
    std::vector<std::array<double,2>> Chix;
    std::vector<std::array<double,2>> Chiy;
    std::vector<std::array<double,2>> Chiz;

    std::vector<std::vector<double>> Energies;
    std::vector<std::vector<double>> Capacity;
    std::vector<std::vector<double>> Entropy;
    std::vector<std::vector<double>> Flux;
    std::array<std::vector<std::vector<double>>,3> Magnetization;
    std::array<std::vector<std::vector<double>>,3> Magnetization2;
    std::array<std::vector<std::vector<double>>,3> Susceptibility;

    std::vector<double> xdata;

    bool Calsusx = false;
    bool Calsusy = false;
    bool Calsusz = false;

    std::array<int,3> get_neighbour_data_hex(int pos);
    std::array<int,3> get_neighbour_data_hex_periodic(int pos);
    std::array<int,3> get_neighbour_data_hex_rev(int pos);
    std::array<int,3> get_neighbour_data_hex_rev2(int pos);
    std::array<int,3> get_neighbour_data_tri(int pos);
    std::array<int,3> get_neighbour_data_tri_periodic(int pos);

    std::array<int,3>(*neighfuncs[6])(int);
    
    std::array<int,3> get_neighbour_data(int pos){
        std::array<int,3> n;

        switch(Lattice_Type){
            case 1:
                n = get_neighbour_data_hex(pos);
                break;
            case 2:
                n = get_neighbour_data_tri(pos);
                break;
            case 3:
                n = get_neighbour_data_hex_periodic(pos);
                break;
            case 4:
                n = get_neighbour_data_hex_rev(pos);
                break;
            case 5:
                n = get_neighbour_data_hex_rev2(pos);
                break;
            case 6: 
                n = get_neighbour_data_tri_periodic(pos);
                break;
        }
        return n;
    }

    void add_kitaev_interaction(std::vector<int>& p_vec, int aux, int sec_aux);
    void add_magnetic_interaction(int aux, int sec_aux);
    void add_heisenberg_interaction(std::vector<int>& p_vec, int aux, int sec_aux);
    void add_gamma_interaction(std::vector<int>& p_vec, int aux, int sec_aux);
    void add_gammaq_interaction(std::vector<int>& p_vec, int aux, int sec_aux);
    MPO honeycomb_flux_operator(int aux, int sec_aux);
    MPO honeycomb_flux_operator_half(int aux, int sec_aux);
    std::array<std::array<MPO,3>,2> magnetization_operators(int aux, int sec_aux);
    int aux_num(int pos, int aux, int sec_aux);
    
    std::vector<std::array<double,2>> Mean(std::vector<std::vector<double>>& M);
    
    int tdvp_loop(std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& S_vec, std::vector<double>& W_vec,
    std::array<std::vector<double>,3>& M_vec, std::array<std::vector<double>,3>& M_vec2,
    MPO& H0, MPS& psi, Cplx& t, int TimeSteps, Args& args, itensor::Sweeps& sweeps, double& cb);

    void time_evolution(std::vector<std::vector<double>>& Energies, std::vector<std::vector<double>>& Capacity,
    std::vector<std::vector<double>>& Entropy, std::vector<std::vector<double>>& Flux,
    std::array<std::vector<std::vector<double>>,3>& Magnetization, std::array<std::vector<std::vector<double>>,3>& Magnetization2,
    std::array<std::vector<std::vector<double>>,3>& Susceptibility,
    std::vector<Cplx> T, std::vector<int> timesteps, int entries, double SusceptDiff, int init_rand_sites, int& max_bond, Args& args, itensor::Sweeps& sweeps);
        
    void tdvp_loop(std::array<std::vector<double>,3>& M_vec,
    MPO& H0, MPS& psi, Cplx& t, int TimeSteps, Args& args, itensor::Sweeps& sweeps);

    void chi_int(MPS& psi, double n, double t, std::array<std::vector<double>,3>& chi_vec, double step, itensor::Sweeps& sweeps, Args& args);

    template<typename T>
    void save_data(std::string filename, std::vector<std::vector<T>>& vec);

    template<std::size_t n, typename T>
    void save_data(std::string filename, std::vector<std::array<T,n>>& vec);

    template<typename T>
    void save_data(std::string filename, std::vector<T>& vec);
    
    template<typename T>
    void save_data(std::string filename, T v);

    std::vector<double> derivative(std::vector<double>& f, double dx);
    std::vector<double> integral(std::vector<double>& f, double dx, double c);
    std::vector<double> multiply(std::vector<double>& a, std::vector<double>& b);

    int dims;

    MPS mpo_to_tanmps(MPO& H);
    MPO tanmps_to_mpo(MPS& X);
    MPO mpo_to_tanmpo(MPO& H);
    MPS mps_to_tanmps(MPS& X);

    int tan_tdvp_loop(int steps, double dt, MPS& Hexptan, MPO& H0tan, MPO& H2tan,
    std::array<MPO,3>& Mtan, std::array<MPO,3>& M2tan, itensor::Sweeps& sweeps, Args& tdvp_args,
    std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& S_vec, std::vector<double>& W_vec,
    std::array<std::vector<double>,3>& M_vec, std::array<std::vector<double>,3>& M_vec2,
    std::array<std::vector<double>,3>& Chi_vec, MPO& Hfluxtan, double& cb, int KrylovExpansions);
    bool SusceptIntegral;

    double calculate_chi(MPS& Hexp, MPS& Hexpinv, std::string spin);
    std::vector<std::array<double,2>> mean_wrap(std::vector<double>& vec);


    public:
    Kitaev_Model(int LX, int LY, Hamiltonian H_Details, int DoubleSpin, int aux, int sec_aux, std::string shape){
        this -> sites = CustomSpinNitsch((LX*LY+2*aux+2*sec_aux),{"2S=",DoubleSpin,"ConserveQNs=",false});
        this -> ampo = AutoMPO(sites);
        this -> H_Details = H_Details;
        this -> dims = DoubleSpin + 1;
        this -> sitestan = CustomSpin((LX*LY+2*aux+2*sec_aux),{"2S=",(dims*dims-1),"ConserveQNs=",false});
        this -> LX = LX;
        this -> LY = LY;
        this -> aux = aux;
        this -> sec_aux = sec_aux;

        std::vector<int> full_points;
        if (shape == "Honeycomb" || shape == "HoneycombPeriodic" || shape == "HoneycombReverse" || shape == "HoneycombReverse2"){
            if (LY%2 != 0){
                std::invalid_argument("LY needs to be divisible by 2 when using the honeycomb lattice");
            }
            if (shape == "Honeycomb"){
                Lattice_Type = 1;
            } else if (shape == "HoneycombPeriodic") {
                Lattice_Type = 3;
            } else if (shape == "HoneycombReverse"){
                Lattice_Type = 4;
            } else if (shape == "HoneycombReverse2"){
                Lattice_Type = 5;
            }
            full_points.reserve(((LY+1)/2)*LX);
            for (int m = 0; m != LX; m++){
                for (int n = 1; n <= LY; n+=2){
                    full_points.push_back(n+LY*m);
                }   
            }
        } else if (shape == "Triangular" || shape == "TriangularPeriodic"){
            if (shape == "Triangular"){
                Lattice_Type = 2;
            } else if (shape == "TriangularPeriodic"){
                Lattice_Type = 6;
            }
            
            full_points.reserve(LX*LY);
            for (int i = 1; i != LX*LY+1; i++){
                full_points.push_back(i);
            }

        } else {
            std::invalid_argument("Invalid Lattice Shape");
        }

        SusceptIntegral = false;

        std::cout << "Lattice Type: " << Lattice_Type << "\n";
        PrintData(ampo);
        add_kitaev_interaction(full_points,aux,sec_aux);
        PrintData(ampo);
        add_magnetic_interaction(aux,sec_aux);
        PrintData(ampo);
        add_heisenberg_interaction(full_points,aux,sec_aux);
        PrintData(ampo);
        add_gamma_interaction(full_points,aux,sec_aux);
        PrintData(ampo);
        add_gammaq_interaction(full_points,aux,sec_aux);
        PrintData(ampo);

        if (DoubleSpin == 1){
            this -> H_flux = honeycomb_flux_operator_half(aux,sec_aux);
        } 
        else{        
            this -> H_flux = honeycomb_flux_operator(aux,sec_aux);
        }
        //PrintData(H_flux);
        
        this -> H0 = toMPO(this -> ampo);
        auto Ms = magnetization_operators(aux,sec_aux);
        M = Ms[0];
        M2 = Ms[1];

        std::cout << "Spin " << DoubleSpin << "/2 System" << "\n";
        std::cout << "Lattice Type: " << shape << "\n";
        std::cout << "LX = " << LX << ", LY = " << LY << "\n";
        std::cout << "Auxiliary Sites: " << aux << "\n";
        std::cout << "Secondary Auxiliary Sites: " << sec_aux << "\n\n";
        std::cout << "Hamiltonian Parameters \n";
        this->H_Details.print();
        std::cout << "\n\n";

    }

    Kitaev_Model(int LX, int LY, Hamiltonian H_Details, int DoubleSpin, std::string shape){
        Kitaev_Model(LX,LY,H_Details,DoubleSpin,0,0,shape);
    }


    void TPQ_MPS(std::vector<int> timesteps, std::vector<double> intervals, int Evols, int max_sites=512, int init_rand_sites=64, std::string TDVP_Type="TwoSite", double SusceptDiff=0, std::string Suscepts="xyz");
    void tanTRG(std::vector<int> timesteps, std::vector<double> intervals, int max_sites=256, int KrylovExpansions=0);

    Hamiltonian Get_Constants(){
        return H_Details;
    }

    void Print_Interactions(){
        PrintData(ampo);
    }

    void Save(std::string x, bool SaveRaw = false){
        std::filesystem::create_directory(x);
        if (CalcDMRG){
            std::string xGSE = x + "/" + "GSE";
            save_data(xGSE,GSE);
        } 
        if (CalcTDVP){
            std::string xd = x + "/" + "xdata";
            std::string xE = x + "/" + "E";
            std::string xC = x + "/" + "C";
            std::string xS = x + "/" + "S";
            std::string xW = x + "/" + "W";
            std::string xcx = x + "/" + "Mx";
            std::string xcy = x + "/" + "My";
            std::string xcz = x + "/" + "Mz";
            std::string xcx2 = x + "/" + "Mx2";
            std::string xcy2 = x + "/" + "My2";
            std::string xcz2 = x + "/" + "Mz2";

            save_data(xd,xdata);
            save_data(xE,E);
            save_data(xC,Cv);
            save_data(xS,S);
            save_data(xW,W);
            save_data(xcx,Mx);
            save_data(xcy,My);
            save_data(xcz,Mz);
            save_data(xcx2,Mx2);
            save_data(xcy2,My2);
            save_data(xcz2,Mz2);

            if (SaveRaw){
                std::string xEr = x + "/" + "E_raw";
                std::string xCr = x + "/" + "C_raw";
                std::string xSr = x + "/" + "S_raw";
                std::string xWr = x + "/" + "W_raw";
                std::string xcxr = x + "/" + "Mx_raw";
                std::string xcyr = x + "/" + "My_raw";
                std::string xczr = x + "/" + "Mz_raw";
                std::string xcx2r = x + "/" + "Mx2_raw";
                std::string xcy2r = x + "/" + "My2_raw";
                std::string xcz2r = x + "/" + "Mz2_raw";

                save_data(xEr,Energies);
                save_data(xCr,Capacity);
                save_data(xSr,Entropy);
                save_data(xWr,Flux);
                save_data(xcxr,Magnetization[0]);
                save_data(xcyr,Magnetization[1]);
                save_data(xczr,Magnetization[2]);
                save_data(xcx2r,Magnetization2[0]);
                save_data(xcy2r,Magnetization2[1]);
                save_data(xcz2r,Magnetization2[2]);
            }

            if (SusceptIntegral){
                std::string xchix = x + "/" + "Chix";
                std::string xchiy = x + "/" + "Chiy";
                std::string xchiz = x + "/" + "Chiz";

                if (Calsusx){
                    save_data(xchix,Chix);
                }
                if (Calsusy){
                    save_data(xchiy,Chiy);
                }
                if (Calsusz){
                    save_data(xchiz,Chiz);
                }

                if (SaveRaw){
                    std::string xchixr = x + "/" + "Chix_raw";
                    std::string xchiyr = x + "/" + "Chiy_raw";
                    std::string xchizr = x + "/" + "Chiz_raw";

                    if (Calsusx){
                        save_data(xchixr,Susceptibility[0]);
                    }
                    if (Calsusy){
                        save_data(xchiyr,Susceptibility[1]);
                    }
                    if (Calsusz){
                        save_data(xchizr,Susceptibility[2]);
                    }                    
                }
            }
        }
        
    }

    
    MPS DMRG(int Sweeps=10){
        CalcDMRG = true;
        auto psi0 = randomMPS(sites,16);
        auto [energy, psi] = dmrg(H0,psi0,Sweeps,{"Quiet",true});
        GSE = energy;
        return psi;
    }





};






}










