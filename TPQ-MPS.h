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
    double GSE;

    private:
    std::vector<std::array<double,2>> E;
    std::vector<std::array<double,2>> Cv;
    std::vector<std::array<double,2>> Cv_alt;
    std::vector<std::array<double,2>> Cv_O2;
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
    std::vector<std::vector<double>> Alternative_Capacity;
    std::vector<std::vector<double>> Alternative_Capacity_O2;
    std::vector<std::vector<double>> Entropy;
    std::vector<std::vector<double>> Flux;
    std::array<std::vector<std::vector<double>>,3> Magnetization;
    std::array<std::vector<std::vector<double>>,3> Magnetization2;
    std::array<std::vector<std::vector<double>>,3> Susceptibility;

    std::vector<double> xdata, xdatashift;

    bool Calsusx = false;
    bool Calsusy = false;
    bool Calsusz = false;

    std::array<int,3> get_neighbour_data_hex(int pos);
    std::array<int,3> get_neighbour_data_hex_periodic(int pos);
    std::array<int,3> get_neighbour_data_hex_rev(int pos);
    std::array<int,3> get_neighbour_data_hex_rev2(int pos);
    std::array<int,3> get_neighbour_data_tri(int pos);
    std::array<int,3> get_neighbour_data_tri_periodic(int pos);

    std::string backup_name;
    std::string Foldername;

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
    
    int tdvp_loop(std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& C_alt_vec, std::vector<double>& C_O2_vec, std::vector<double>& S_vec, std::vector<double>& W_vec,
    std::array<std::vector<double>,3>& M_vec, std::array<std::vector<double>,3>& M_vec2,
    MPO& H0, MPS& psi, Cplx& t, int TimeSteps, Args& args, itensor::Sweeps& sweeps, double& cb, int evnum, int secnum);

    void time_evolution(std::vector<Cplx> T, std::vector<int> timesteps, int entries, double SusceptDiff, int init_rand_sites, int& max_bond, Args& args, itensor::Sweeps& sweeps, int evnum);
    void time_evolution_cont(std::vector<Cplx> T, std::vector<int> timesteps, int entries, double SusceptDiff, int init_rand_sites, int& max_bond, Args& args, Sweeps sweeps, int evnum);

    void tdvp_loop(std::array<std::vector<double>,3>& M_vec, MPO& H0, MPS& psi,
    Cplx& t, int TimeSteps, Args& args, itensor::Sweeps& sweeps, std::string axis, int evnum, int secnum);
    void sus_func(std::array<std::vector<double>,3>& M_vec, MPO& H0, MPS& psi, 
    std::vector<itensor::Cplx>& T, std::vector<int> timesteps, Args& args, itensor::Sweeps& sweeps, std::string axis, int evnum);

    void chi_int(MPS& psi, double n, double t, std::array<std::vector<double>,3>& chi_vec, double step, itensor::Sweeps& sweeps, Args& args);

    public:
    template<typename T>
    void save_data(std::string filename, std::vector<std::vector<T>>& vec);

    template<std::size_t n, typename T>
    void save_data(std::string filename, std::vector<std::array<T,n>>& vec);

    template<typename T>
    void save_data(std::string filename, std::vector<T>& vec);
    
    template<typename T>
    void save_data(std::string filename, T v);

    private:
    std::vector<double> derivative(std::vector<double>& f, double dx);
    std::vector<double> integral(std::vector<double>& f, double dx, double c);
    std::vector<double> multiply(std::vector<double>& a, std::vector<double>& b);

    template<typename T>
    void load_vector(std::string filename, std::vector<T>& v);

    int dims;

    MPS mpo_to_tanmps(MPO& H);
    MPO tanmps_to_mpo(MPS& X);
    MPO mpo_to_tanmpo(MPO& H);
    MPS mps_to_tanmps(MPS& X);

    int tan_tdvp_loop(int steps, double dt, MPS& Hexptan, MPO& H0tan, MPO& H2tan,
    std::array<MPO,3>& Mtan, std::array<MPO,3>& M2tan, itensor::Sweeps& sweeps, Args& tdvp_args,
    std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& C_alt_vec,
    std::vector<double>& S_vec, std::vector<double>& W_vec,
    std::array<std::vector<double>,3>& M_vec, std::array<std::vector<double>,3>& M_vec2,
    MPO& Hfluxtan, double& cb, int KrylovExpansions);
    bool SusceptIntegral;

    double calculate_chi(MPS& Hexp, MPS& Hexpinv, std::string spin);
    std::vector<std::array<double,2>> mean_wrap(std::vector<double>& vec);


    public:
    Kitaev_Model(int LX, int LY, Hamiltonian H_Details, int DoubleSpin, int aux, int sec_aux, std::string shape, std::string backup = "NULL"){
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

        add_kitaev_interaction(full_points,aux,sec_aux);
        add_magnetic_interaction(aux,sec_aux);
        add_heisenberg_interaction(full_points,aux,sec_aux);
        add_gamma_interaction(full_points,aux,sec_aux);
        add_gammaq_interaction(full_points,aux,sec_aux);

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

        if (backup == "NULL"){
            srand(time(NULL));
            this -> backup_name = std::to_string(rand()) + "_backup";
        }
        else{
            this -> backup_name = backup + "_backup";
        }

        std::filesystem::create_directory(backup_name);
        save_data(backup_name+"/LX",LX);
        save_data(backup_name+"/LY",LY);
        save_data(backup_name+"/DS",DoubleSpin);
        save_data(backup_name+"/aux",aux);
        save_data(backup_name+"/sec_aux",sec_aux);
        save_data(backup_name+"/LT",Lattice_Type);

        std::vector<double> Hinf = this->H_Details.total_info();
        save_data(backup_name+"/H",Hinf);
        writeToFile(backup_name+"/sites",this->sites);

        std::cout << "Spin " << DoubleSpin << "/2 System" << "\n";
        std::cout << "Lattice Type: " << shape << "\n";
        std::cout << "LX = " << LX << ", LY = " << LY << "\n";
        std::cout << "Auxiliary Sites: " << aux << "\n";
        std::cout << "Secondary Auxiliary Sites: " << sec_aux << "\n\n";
        std::cout << "Hamiltonian Parameters \n";
        this->H_Details.print();
        std::cout << "\n\n";

    }


    void TPQ_MPS(std::vector<int> timesteps, std::vector<double> intervals, int Evols, int max_sites=512, int init_rand_sites=64, std::string TDVP_Type="TwoSite", double SusceptDiff=0, std::string Suscepts="xyz");
    void tanTRG(std::vector<int> timesteps, std::vector<double> intervals, int max_sites=256, int KrylovExpansions=0);

    Kitaev_Model(std::string backup_name);

    Hamiltonian Get_Constants(){
        return H_Details;
    }

    void Print_Interactions(){
        PrintData(ampo);
    }

    void Save(bool SaveRaw = false){
        std::filesystem::create_directory(Foldername);
        if (CalcDMRG){
            std::string xGSE = Foldername + "/" + "GSE";
            save_data(xGSE,GSE);
        } 
        if (CalcTDVP){
            std::string xd = Foldername + "/" + "xdata";
            std::string xds = Foldername + "/" + "xdatashift";
            std::string xE = Foldername + "/" + "E";
            std::string xC = Foldername + "/" + "C";
            std::string xCalt = Foldername + "/" + "C_alt";
            std::string xCo2 = Foldername + "/" + "C_alt_O2";
            std::string xS = Foldername + "/" + "S";
            std::string xW = Foldername + "/" + "W";
            std::string xcx = Foldername + "/" + "Mx";
            std::string xcy = Foldername + "/" + "My";
            std::string xcz = Foldername + "/" + "Mz";
            std::string xcx2 = Foldername + "/" + "Mx2";
            std::string xcy2 = Foldername + "/" + "My2";
            std::string xcz2 = Foldername + "/" + "Mz2";

            save_data(xd,xdata);
            save_data(xds,xdatashift);
            save_data(xE,E);
            save_data(xC,Cv);
            save_data(xCalt,Cv_alt);
            save_data(xCo2,Cv_O2);
            save_data(xS,S);
            save_data(xW,W);
            save_data(xcx,Mx);
            save_data(xcy,My);
            save_data(xcz,Mz);
            save_data(xcx2,Mx2);
            save_data(xcy2,My2);
            save_data(xcz2,Mz2);

            if (SaveRaw){
                std::string xEr = Foldername + "/" + "E_raw";
                std::string xCr = Foldername + "/" + "C_raw";
                std::string xCaltr = Foldername + "/" + "C_alt_raw";
                std::string xCo2 = Foldername + "/" + "C_alt_O2_raw";
                std::string xSr = Foldername + "/" + "S_raw";
                std::string xWr = Foldername + "/" + "W_raw";
                std::string xcxr = Foldername + "/" + "Mx_raw";
                std::string xcyr = Foldername + "/" + "My_raw";
                std::string xczr = Foldername + "/" + "Mz_raw";
                std::string xcx2r = Foldername + "/" + "Mx2_raw";
                std::string xcy2r = Foldername + "/" + "My2_raw";
                std::string xcz2r = Foldername + "/" + "Mz2_raw";

                save_data(xEr,Energies);
                save_data(xCr,Capacity);
                save_data(xCaltr,Alternative_Capacity);
                save_data(xCo2,Alternative_Capacity_O2);
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
                std::string xchix = Foldername + "/" + "Chix";
                std::string xchiy = Foldername + "/" + "Chiy";
                std::string xchiz = Foldername + "/" + "Chiz";

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
                    std::string xchixr = Foldername + "/" + "Chix_raw";
                    std::string xchiyr = Foldername + "/" + "Chiy_raw";
                    std::string xchizr = Foldername + "/" + "Chiz_raw";

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
        
        std::filesystem::remove_all(backup_name);
    }


    void Save(int world_rank){
        std::string xnew = Foldername + std::to_string(world_rank);
        Save(xnew,true);
    }

    void Save(std::string x, bool SaveRaw = false){
        Foldername = x;
        Save(SaveRaw);
    }

    void Save(std::string x, int world_rank){
        Foldername = x;
        Save(world_rank);
    }

    void Set_Foldername(std::string x){
        Foldername = x;
        save_data(backup_name+"/foldername",Foldername);
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










