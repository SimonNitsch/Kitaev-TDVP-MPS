#include "TPQ-MPS.h"
#pragma once


static std::mutex forloop_mutex;
static std::mutex suscept_mutex;

namespace TDVP_MPS{


std::array<int,3> Kitaev_Model::get_neighbour_data_hex(int pos){
    int Y = (pos-1) % LY;
    int X = (pos-1) / LY;
    std::array<int,3> neighbours{};

    if (Y != (LY-1)){
        neighbours[0] = pos+1;
    }
    if (X != 0){
        neighbours[1] = pos-LY+1;
    }
    if (pos != 1 && pos != LX*LY){
        neighbours[2] = pos-1;
    }
    return neighbours;

}


std::array<int,3> Kitaev_Model::get_neighbour_data_hex_periodic(int pos){
    int Y = (pos-1) % LY;
    int X = (pos-1) / LY;
    std::array<int,3> neighbours{};

    if (Y != (LY-1)){
        neighbours[0] = pos+1;
    }
    if (X != 0){
        neighbours[1] = pos-LY+1;
    }
    else{
        neighbours[1] = pos+(LX-1)*LY+1;
    }
    if (Y != 0){
        neighbours[2] = pos-1;
    } 
    else{
        neighbours[2] = pos+LY-1;
    }
    return neighbours;

} 


std::array<int,3> Kitaev_Model::get_neighbour_data_hex_rev(int pos){
    int Y = (pos-1) % LY;
    int X = (pos-1) / LY;
    std::array<int,3> neighbours{};

    if (X != LX-1){
        neighbours[0] = pos+LY+1;
    }
    else{
        neighbours[0] = pos-(LX-1)*LY+1;
    }

    neighbours[1] = pos+1;

    if (Y != 0){
        neighbours[2] = pos-1;
    }
    else{
        neighbours[2] = pos+LY-1;
    }
    return neighbours;

}


std::array<int,3> Kitaev_Model::get_neighbour_data_hex_rev2(int pos){
    int Y = LY - 1 - (pos-1) % LY;
    int X = (pos-1) / LY;
    std::array<int,3> neighbours{};

    if (X != 0){
        neighbours[0] = pos-LY+1;
    }
    else {
        neighbours[0] = pos+(LX-1)*LY+1;
    }
    
    neighbours[1] = pos+1;

    if (Y != LY-1){
        neighbours[2] = pos-1;
    }
    else {
        neighbours[2] = pos+LY-1;
    } 
    return neighbours;

}


std::array<int,3> Kitaev_Model::get_neighbour_data_tri(int pos){
    int X = (pos-1) / LY;
    int Y = (pos-1) % LY;
    std::array<int,3> neighbours{};

    if (X != (LX-1)){
        neighbours[0] = pos + LY;

        if (Y != 0){
            neighbours[2] = pos + LY - 1;
        }
    }

    if (Y != LY-1 || X != LX-1){
        neighbours[1] = pos + 1;
    }

    return neighbours;
}




std::array<int,3> Kitaev_Model::get_neighbour_data_tri_periodic(int pos){
    int X = (pos-1) / LY;
    int Y = (pos-1) % LY;
    std::array<int,3> neighbours{};

    if (X != (LX-1)){
        neighbours[0] = pos + LY;
        if (Y == 0){
            neighbours[2] = pos + 2*LY - 1;
        }
        else{
            neighbours[2] = pos + LY - 1;
        }
    }
    else{
        neighbours[0] = pos - (LX-1) * LY;
        if (Y == 0){
            neighbours[2] = pos - (LX-2) * LY - 1;
        }
        else{
            neighbours[2] = pos - (LX-1) * LY - 1;
        }
    }
    
    if (Y == (LY-1)){
        neighbours[1] = pos + 1 - LY;
    }
    else{
        neighbours[1] = pos + 1;
    }

    return neighbours;
}





void Kitaev_Model::add_kitaev_interaction(std::vector<int>& p_vec, int aux, int sec_aux){
    double Kx = H_Details.get("Kx");
    double Ky = H_Details.get("Ky");
    double Kz = H_Details.get("Kz");
    
    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(i);
        
        if (n[0] != 0){
            ampo += Kx,"Sx",aux_num(i,aux,sec_aux),"Sx",aux_num(n[0],aux,sec_aux);
        }
        if (n[1] != 0){
            ampo += Ky,"Sy",aux_num(i,aux,sec_aux),"Sy",aux_num(n[1],aux,sec_aux);
        }
        if (n[2] != 0){
            ampo += Kz,"Sz",aux_num(i,aux,sec_aux),"Sz",aux_num(n[2],aux,sec_aux);
        }
    }

}


void Kitaev_Model::add_magnetic_interaction(int aux, int sec_aux){
    double hx = H_Details.get("hx");
    double hy = H_Details.get("hy");
    double hz = H_Details.get("hz");

    for (int i = 1; i != LX*LY+1; i++){
        ampo += -1*hx,"Sx",aux_num(i,aux,sec_aux);
        ampo += -1*hy,"Sy",aux_num(i,aux,sec_aux);
        ampo += -1*hz,"Sz",aux_num(i,aux,sec_aux);
    }
}



void Kitaev_Model::add_heisenberg_interaction(std::vector<int>& p_vec, int aux, int sec_aux){
    double J = H_Details.get("J");

    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(i);
        if (n[0] != 0){
            ampo += J,"Sz",aux_num(i,aux,sec_aux),"Sz",aux_num(n[0],aux,sec_aux);
            ampo += J*0.5,"S+",aux_num(i,aux,sec_aux),"S-",aux_num(n[0],aux,sec_aux);
            ampo += J*0.5,"S-",aux_num(i,aux,sec_aux),"S+",aux_num(n[0],aux,sec_aux);
        }
        if (n[1] != 0){
            ampo += J,"Sz",aux_num(i,aux,sec_aux),"Sz",aux_num(n[1],aux,sec_aux);
            ampo += J*0.5,"S+",aux_num(i,aux,sec_aux),"S-",aux_num(n[1],aux,sec_aux);
            ampo += J*0.5,"S-",aux_num(i,aux,sec_aux),"S+",aux_num(n[1],aux,sec_aux);
        }
        if (n[2] != 0){
            ampo += J,"Sz",aux_num(i,aux,sec_aux),"Sz",aux_num(n[2],aux,sec_aux);
            ampo += J*0.5,"S+",aux_num(i,aux,sec_aux),"S-",aux_num(n[2],aux,sec_aux);
            ampo += J*0.5,"S-",aux_num(i,aux,sec_aux),"S+",aux_num(n[2],aux,sec_aux);
        }
    }

}


void Kitaev_Model::add_gamma_interaction(std::vector<int>& p_vec, int aux, int sec_aux){
    double G = H_Details.get("Gamma");

    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(i);
        if (n[0] != 0){
            ampo += G,"Sy",aux_num(i,aux,sec_aux),"Sz",aux_num(n[0],aux,sec_aux);
            ampo += G,"Sz",aux_num(i,aux,sec_aux),"Sy",aux_num(n[0],aux,sec_aux);
        }
        if (n[1] != 0){
            ampo += G,"Sx",aux_num(i,aux,sec_aux),"Sz",aux_num(n[1],aux,sec_aux);
            ampo += G,"Sz",aux_num(i,aux,sec_aux),"Sx",aux_num(n[1],aux,sec_aux);
        }
        if (n[2] != 0){
            ampo += G,"Sx",aux_num(i,aux,sec_aux),"Sy",aux_num(n[2],aux,sec_aux);
            ampo += G,"Sy",aux_num(i,aux,sec_aux),"Sx",aux_num(n[2],aux,sec_aux);
        }
    }

}



void Kitaev_Model::add_gammaq_interaction(std::vector<int>& p_vec, int aux, int sec_aux){
    double GQ = H_Details.get("GammaQ");

    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(i);
        if (n[0] != 0){
            ampo += GQ,"Sy",aux_num(i,aux,sec_aux),"Sx",aux_num(n[0],aux,sec_aux);
            ampo += GQ,"Sx",aux_num(i,aux,sec_aux),"Sy",aux_num(n[0],aux,sec_aux);
            ampo += GQ,"Sz",aux_num(i,aux,sec_aux),"Sx",aux_num(n[0],aux,sec_aux);
            ampo += GQ,"Sx",aux_num(i,aux,sec_aux),"Sz",aux_num(n[0],aux,sec_aux);
        }
        if (n[1] != 0){
            ampo += GQ,"Sx",aux_num(i,aux,sec_aux),"Sy",aux_num(n[1],aux,sec_aux);
            ampo += GQ,"Sy",aux_num(i,aux,sec_aux),"Sx",aux_num(n[1],aux,sec_aux);
            ampo += GQ,"Sz",aux_num(i,aux,sec_aux),"Sy",aux_num(n[1],aux,sec_aux);
            ampo += GQ,"Sy",aux_num(i,aux,sec_aux),"Sz",aux_num(n[1],aux,sec_aux);
        }
        if (n[2] != 0){
            ampo += GQ,"Sx",aux_num(i,aux,sec_aux),"Sz",aux_num(n[2],aux,sec_aux);
            ampo += GQ,"Sz",aux_num(i,aux,sec_aux),"Sx",aux_num(n[2],aux,sec_aux);
            ampo += GQ,"Sz",aux_num(i,aux,sec_aux),"Sy",aux_num(n[2],aux,sec_aux);
            ampo += GQ,"Sy",aux_num(i,aux,sec_aux),"Sz",aux_num(n[2],aux,sec_aux);
        }
    }

}


int Kitaev_Model::aux_num(int pos, int aux, int sec_aux){
    int totaux = pos + aux;

    if (pos > LY){
        totaux += sec_aux;
    }
    if (pos > (LX-1)*LY){
        totaux += sec_aux;
    }

    return totaux;
}



itensor::MPO Kitaev_Model::honeycomb_flux_operator_half(int aux, int sec_aux){
    int Py = std::max((LY/2)-1,0);
    int Px = LX-1;
    auto fluxop = itensor::AutoMPO(sites);
    double Wfac = 64. / static_cast<double>(Px) / static_cast<double>(Py);

    for (int i = 0; i != Px; i++){
        for (int j = 0; j != Py; j++){
            int f_prot = LY*i + 2*j + 2;

            int f1 = aux_num(f_prot,aux,sec_aux);
            int f2 = aux_num(f_prot+1,aux,sec_aux);
            int f3 = aux_num(f_prot+2,aux,sec_aux);
            int f4 = aux_num(f_prot+LY+1,aux,sec_aux);
            int f5 = aux_num(f_prot+LY,aux,sec_aux);
            int f6 = aux_num(f_prot+LY-1,aux,sec_aux);

            fluxop += Wfac,"Sx",f1,"Sy",f2,"Sz",f3,"Sx",f4,"Sy",f5,"Sz",f6;
        }
    }
    auto fluxH = itensor::toMPO(fluxop);
    return fluxH;
}




itensor::MPO Kitaev_Model::honeycomb_flux_operator(int aux, int sec_aux){
    int Py = std::max((LY/2)-1,0);
    int Px = LX-1;
    auto fluxop = itensor::AutoMPO(sites);
    double Wfac = 1. / static_cast<double>(Px) / static_cast<double>(Py);
    //std::cout << "Px: " << Px << " Py: " << Py << "\n" << std::flush;

    for (int i = 0; i != Px; i++){
        for (int j = 0; j != Py; j++){
            int f_prot = LY*i + 2*j + 2;

            int f1 = aux_num(f_prot,aux,sec_aux);
            int f2 = aux_num(f_prot+1,aux,sec_aux);
            int f3 = aux_num(f_prot+2,aux,sec_aux);
            int f4 = aux_num(f_prot+LY+1,aux,sec_aux);
            int f5 = aux_num(f_prot+LY,aux,sec_aux);
            int f6 = aux_num(f_prot+LY-1,aux,sec_aux);
            std::cout << Wfac << ", " << f1 << ", " << f2 << ", " << f3 << ", " << f4 << ", " << f5 << ", " << f6 << "\n" << std::flush;

            fluxop += Wfac,"Expz",f1,"Expz",f2,"Expx",f2,"Expx",f3,"Expy",f3,"Expy",f4,"Expz",f4,"Expz",f5,"Expx",f5,"Expx",f6,"Expy",f6,"Expy",f1;
        }
    }    
    itensor::PrintData(fluxop); 
    auto fluxH = itensor::toMPO(fluxop);
    return fluxH;
}





std::array<std::array<itensor::MPO,3>,2> Kitaev_Model::magnetization_operators(int aux, int sec_aux){
    auto mx = itensor::AutoMPO(sites);
    auto my = itensor::AutoMPO(sites);
    auto mz = itensor::AutoMPO(sites);
    auto mx2 = itensor::AutoMPO(sites);
    auto my2 = itensor::AutoMPO(sites);
    auto mz2 = itensor::AutoMPO(sites);

    for (int i = 1; i <= LX*LY; i++){
        mx += "Sx",aux_num(i,aux,sec_aux);
        my += "Sy",aux_num(i,aux,sec_aux);
        mz += "Sz",aux_num(i,aux,sec_aux);
        for (int j = 1; j < i; j++){
            mx2 += "Sx",aux_num(i,aux,sec_aux),"Sx",aux_num(j,aux,sec_aux);
            my2 += "Sy",aux_num(i,aux,sec_aux),"Sy",aux_num(j,aux,sec_aux);
            mz2 += "Sz",aux_num(i,aux,sec_aux),"Sz",aux_num(j,aux,sec_aux);
        }
    }
    std::array<itensor::MPO,3> m1 = {itensor::toMPO(mx),itensor::toMPO(my),itensor::toMPO(mz)};
    std::array<itensor::MPO,3> m2 = {itensor::toMPO(mx2),itensor::toMPO(my2),itensor::toMPO(mz2)};
    std::array<std::array<itensor::MPO,3>,2> m = {m1,m2};

    return m;
}







std::vector<std::array<double,2>> Kitaev_Model::Mean(std::vector<std::vector<double>>& M){
    int M0 = M[0].size();
    int M1 = M.size();
    std::vector<std::array<double,2>> vec;
    vec.reserve(M0-1);
    
    for (int i = 1; i != M0; i++){
        double v = 0;
        double v2 = 0;
        for (int j = 0; j != M1; j++){
            v += M[j][i];
            v2 += M[j][i] * M[j][i];
        }
        double mean = v/static_cast<double>(M1);
        double std = std::sqrt(v2/M1 - mean*mean);
        std::array<double,2> proto_vec = {mean,std};
        vec.emplace_back(proto_vec);
    }
    
    return vec;

}



int Kitaev_Model::tdvp_loop(std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& C_alt_vec, std::vector<double>& C_O2_vec, std::vector<double>& S_vec, std::vector<double>& W_vec,
std::array<std::vector<double>,3>& M_vec, std::array<std::vector<double>,3>& M_vec2,
itensor::MPO& H0, itensor::MPS& psi, itensor::Cplx& t, int TimeSteps, itensor::Args& args, itensor::Sweeps& sweeps, double& cb, int evnum, int secnum){
    std::complex<double> tcompl2 = t;
    double t_beta = std::real(tcompl2) * -2;
    int max_bond = 0;

    for (int j = 0; j != TimeSteps; j++){
        double cb_old = cb;
        if (cb_old == 0){
            cb_old = 1;
        }
        double cbhalf = cb + t_beta / 2.;
        cb += t_beta;

        double E = itensor::tdvp(psi,H0,t,sweeps,args);
        double E2 = std::real(itensor::innerC(psi,H0,H0,psi));
        max_bond = std::max(max_bond,itensor::maxLinkDim(psi));

        std::complex<double> w = itensor::innerC(psi,H_flux,psi);
        double calt = cb * cb * (E_vec.back() - E) / t_beta;
        double co2 = cbhalf * cbhalf * (E_vec.back() - E) / t_beta;
        double c = cb * cb * (E2 - E * E);
        double s = S_vec.back() + cb_old * E_vec.back() - cb * E + t_beta * 0.5 * (E + E_vec.back());

        E_vec.push_back(E);
        C_vec.push_back(c);
        C_alt_vec.push_back(calt);
        C_O2_vec.push_back(co2);
        S_vec.push_back(s);
        W_vec.push_back(std::real(w)); 
 
        std::cout << "Current Beta: " << cb << ", Energy: " << E << ", Heat Capacity: " << c << ", Alternative Heat Capacity: " << calt << "\n" << std::flush;


        for (int indi = 0; indi != 3; indi++){
            std::complex<double> m = itensor::innerC(psi,M[indi],psi);
            M_vec[indi].push_back(std::real(m));
            std::complex<double> m2 = itensor::innerC(psi,M2[indi],psi);
            M_vec2[indi].push_back(std::real(m2));
        }
    }
    
    std::string datf = backup_name + "/" + std::to_string(evnum);
    std::string datfbak = datf + "bak" + std::to_string(secnum);
    if (std::filesystem::exists(datfbak)){
        std::filesystem::remove_all(datfbak);
    }

    std::filesystem::rename(datf, datfbak);
    std::this_thread::sleep_for(std::chrono::seconds(2));

    std::filesystem::create_directory(datf);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));

    save_data(datf + "/Evec", E_vec);
    save_data(datf + "/Cvec", C_vec);
    save_data(datf + "/Caltvec", C_alt_vec);
    save_data(datf + "/CO2vec", C_O2_vec);
    save_data(datf + "/Svec", S_vec);
    save_data(datf + "/Wvec", W_vec);
    save_data(datf + "/Mvecx", M_vec[0]);
    save_data(datf + "/Mvecy", M_vec[1]);
    save_data(datf + "/Mvecz", M_vec[2]);
    save_data(datf + "/Mvecx2", M_vec2[0]);
    save_data(datf + "/Mvecy2", M_vec2[1]);
    save_data(datf + "/Mvecz2", M_vec2[2]);
    save_data(datf + "/section", secnum);

    itensor::writeToFile(datf + "/psi",psi);

    if (secnum != 0){
        std::filesystem::remove_all(datf + "bak" + std::to_string(secnum-1));
    }

    return max_bond;
}


void Kitaev_Model::tdvp_loop(std::array<std::vector<double>,3>& M_vec, itensor::MPO& H0, itensor::MPS& psi,
itensor::Cplx& t, int TimeSteps, itensor::Args& args, itensor::Sweeps& sweeps, std::string axis, int evnum, int secnum){
    std::complex<double> tcompl2 = t;
    double t_beta = std::real(tcompl2) * -2;

    for (int j = 0; j != TimeSteps; j++){
        double E = itensor::tdvp(psi,H0,t,sweeps,args);

        for (int indi = 0; indi != 3; indi++){
            std::complex<double> m = itensor::innerC(psi,M[indi],psi);
            M_vec[indi].push_back(std::real(m));
        }
    }

    std::string datf = backup_name + "/" + std::to_string(evnum) + axis;
    std::filesystem::rename(datf, datf + "bak" + std::to_string(secnum));
    std::filesystem::create_directory(datf);

    save_data(datf + "/Mvecx", M_vec[0]);
    save_data(datf + "/Mvecy", M_vec[1]);
    save_data(datf + "/Mvecz", M_vec[2]);

}



void Kitaev_Model::sus_func(std::array<std::vector<double>,3>& M_vec, itensor::MPO& H0, itensor::MPS& psi, 
std::vector<itensor::Cplx>& T, std::vector<int> timesteps, Args& args, itensor::Sweeps& sweeps, std::string axis, int evnum){
    for (int j = 0; j != timesteps.size(); j++){
        tdvp_loop(M_vec,H0,psi,T[j],timesteps[j],args,sweeps,axis,evnum,j);
    }
}


template<typename T>
void Kitaev_Model::save_data(std::string filename, std::vector<std::vector<T>>& vec){
    std::string ending = ".txt";
    std::string full_file = filename + ending;
    std::ofstream file(full_file);

    for (auto& varr : vec){
        auto v = varr.begin();
        file << *v;
        v++;

        for (; v != varr.end(); v++){
            file << ", " << *v;
        }
        file << "\n";

    }
    file.close();
    std::cout << "Data saved as: " << full_file << "\n" << std::flush;
}




template<std::size_t n, typename T>
void Kitaev_Model::save_data(std::string filename, std::vector<std::array<T,n>>& vec){
    std::string ending = ".txt";
    std::string full_file = filename + ending;
    std::ofstream file(full_file);

    for (auto& varr : vec){
        auto v = varr.begin();
        file << *v;
        v++;

        for (; v != varr.end(); v++){
            file << ", " << *v;
        }
        file << "\n";

    }
    file.close();
    std::cout << "Data saved as: " << full_file << "\n" << std::flush;
}


template<typename T>
void Kitaev_Model::save_data(std::string filename, std::vector<T>& vec){
    std::vector<std::array<T,1>> arrvec;
    arrvec.reserve(vec.size());

    for (auto& v : vec){
        std::array<T,1> a = {v};
        arrvec.push_back(a);
    }
    save_data(filename,arrvec);
}


template<typename T>
void Kitaev_Model::save_data(std::string filename, T v){
    std::vector<std::array<T,1>> vec;
    std::array<T,1> v_array = {v};
    vec.push_back(v_array);

    save_data(filename, vec);
}


void Kitaev_Model::chi_int(itensor::MPS& psi, double n, double t, std::array<std::vector<double>,3>& chi_vec, double step, itensor::Sweeps& sweeps, itensor::Args& args){
    int steps = std::ceil(t);
    itensor::Cplx dt = t / static_cast<double>(steps) * itensor::Cplx_1;
    auto appl_args = itensor::Args({"Cutoff=",1e-4,"MaxDim=",32});

    itensor::MPS sp;
    sp = itensor::applyMPO(M[0],psi,appl_args);
    itensor::println(sp);
    std::cout << 3333333 << "\n" << std::flush;
    std::complex<double> ss = itensor::innerC(sp,sp);
    std::cout << ss << "\n" << std::flush;
    sp.normalize();
    std::cout << steps << " " << dt << "\n" << std::flush;


    double E = itensor::tdvp(sp,H0,dt,sweeps,args);
    std::cout << 0 << " " << std::flush;
    
}



void Kitaev_Model::time_evolution(std::vector<Cplx> T, std::vector<int> timesteps, int entries, double SusceptDiff,
                                    int init_rand_sites,int& max_bond, itensor::Args& args, itensor::Sweeps& sweeps, int evnum){

    std::vector<double> E_vec;
    std::vector<double> C_vec;
     std::vector<double> C_O2_vec;
    std::vector<double> C_alt_vec;
    std::vector<double> S_vec;
    std::vector<double> W_vec;

    E_vec.reserve(entries);
    C_vec.reserve(entries);
    C_alt_vec.reserve(entries);
    C_O2_vec.reserve(entries);
    S_vec.reserve(entries);
    W_vec.reserve(entries);

    std::array<std::vector<double>,3> Mag_vec;
    for (auto& i : Mag_vec){
        i.reserve(entries);
    }
    std::array<std::vector<double>,3> Mag_vec2;
    for (auto& i : Mag_vec2){
        i.reserve(entries);
    }

    std::array<std::vector<double>,3> Mag_vec_nextx;
    std::array<std::vector<double>,3> Mag_vec_nexty;
    std::array<std::vector<double>,3> Mag_vec_nextz;


    double curr_beta = 0;
    auto t1 = std::chrono::system_clock::now();
    auto psi = itensor::randomMPS(sites,init_rand_sites);
    MPS psix = psi;
    MPS psiy = psi;
    MPS psiz = psi;

    std::complex<double> E = itensor::innerC(psi,H0,psi) / itensor::inner(psi,psi);
    E_vec.push_back(std::real(E));
    W_vec.push_back(0);
    C_vec.push_back(0);
    C_alt_vec.push_back(0);
    C_O2_vec.push_back(0);
    S_vec.push_back(0);
    for (auto& i : Mag_vec){
        i.push_back(0);
    }
    for (auto& i : Mag_vec2){
        i.push_back(0);
    }
        //double n = 1;

    std::vector<std::future<void>> SusFut;

    std::filesystem::create_directory(backup_name + "/" + std::to_string(evnum));

    if (SusceptIntegral){
        SusFut.reserve(3);
        for (auto& i : Mag_vec_nextx){
            i.reserve(entries);
        }
        for (auto& i : Mag_vec_nexty){
            i.reserve(entries);
        }
        for (auto& i : Mag_vec_nextz){
            i.reserve(entries);
        }

        for (auto& i : Mag_vec_nextx){
            i.push_back(0);
        }
        for (auto& i : Mag_vec_nextz){
            i.push_back(0);
        }
        for (auto& i : Mag_vec_nexty){
            i.push_back(0);
        }

        if (Calsusx){
            std::filesystem::create_directory(backup_name + "/" + std::to_string(evnum) + "x");
            SusFut.push_back(std::async(std::launch::async,[&](){sus_func(Mag_vec_nextx,H0x,psix,T,timesteps,args,sweeps,"x",evnum);}));
        }
        if (Calsusy){
            std::filesystem::create_directory(backup_name + "/" + std::to_string(evnum) + "y");
            SusFut.push_back(std::async(std::launch::async,[&](){sus_func(Mag_vec_nexty,H0y,psiy,T,timesteps,args,sweeps,"y",evnum);}));
        }
        if (Calsusz){
            std::filesystem::create_directory(backup_name + "/" + std::to_string(evnum) + "z");
            SusFut.push_back(std::async(std::launch::async,[&](){sus_func(Mag_vec_nextz,H0z,psiz,T,timesteps,args,sweeps,"z",evnum);}));
        }

    }


    for (int j = 0; j != timesteps.size(); j++){
        int curbond = tdvp_loop(E_vec,C_vec,C_alt_vec,C_O2_vec,S_vec,W_vec,Mag_vec,Mag_vec2,H0,psi,T[j],timesteps[j],args,sweeps,curr_beta,evnum,j);
        max_bond = std::max(max_bond,curbond);
    }
    S_vec = std::log(dims) * static_cast<double>(LX) * static_cast<double>(LY) - S_vec;

    {
        std::lock_guard<std::mutex> lock(forloop_mutex);
        Energies.push_back(E_vec);
        E_vec.clear();
        Capacity.push_back(C_vec);
        C_vec.clear();
        Alternative_Capacity.push_back(C_alt_vec);
        C_alt_vec.clear();
        Alternative_Capacity_O2.push_back(C_O2_vec);
        C_O2_vec.clear();
        Entropy.push_back(S_vec);
        S_vec.clear();
        Flux.push_back(W_vec);
        W_vec.clear();

        for (int j = 0; j != 3; j++){
            Magnetization[j].push_back(Mag_vec[j]);
            Magnetization2[j].push_back(Mag_vec2[j]);
        }
    }


    if (SusceptIntegral){
        for (auto& sf : SusFut){
            sf.wait();
        }

        std::lock_guard<std::mutex> lock2(suscept_mutex);
        if (Calsusx){
            std::vector<double> chix = Mag_vec_nextx[0] - Mag_vec[0];
            Susceptibility[0].push_back(chix/SusceptDiff);
        }
        if (Calsusy){
            std::vector<double> chiy = Mag_vec_nexty[1] - Mag_vec[1];
            Susceptibility[1].push_back(chiy/SusceptDiff);
        }
        if (Calsusz){
            std::vector<double> chiz = Mag_vec_nextz[2] - Mag_vec[2];
            Susceptibility[2].push_back(chiz/SusceptDiff);
        }
    }

    for (int j = 0; j != 3; j++){
        Mag_vec[j].clear();
        Mag_vec2[j].clear();
    }

        

    auto t2 = std::chrono::system_clock::now();
    auto time = std::chrono::duration<double>(t2-t1);
    std::cout << "Finished Evolution ID " << std::this_thread::get_id() << ", Time Needed: " << time.count() << " Seconds\n" << std::flush;



}




void Kitaev_Model::time_evolution_cont(std::vector<itensor::Cplx> T, std::vector<int> timesteps, int entries, double SusceptDiff, int init_rand_sites, int& max_bond, itensor::Args& args, itensor::Sweeps sweeps, int evnum){

    std::string cfname = backup_name + "/" + std::to_string(evnum);

    int sec;
    std::ifstream file1(cfname + "/section.txt");
    file1 >> sec;
    file1.close();

    std::vector<double> E_vec, C_vec, C_alt_vec, C_O2_vec, S_vec, W_vec;
    E_vec.reserve(entries);
    C_vec.reserve(entries);
    C_alt_vec.reserve(entries);
    C_O2_vec.reserve(entries);
    S_vec.reserve(entries);
    W_vec.reserve(entries);
    
    std::array<std::vector<double>,3> Mag_vec, Mag_vec2;
    for (int i = 0; i != 3; i++){
        Mag_vec[i].reserve(entries);
        Mag_vec2[i].reserve(entries);
    }

    load_vector(cfname + "/Evec.txt", E_vec);
    load_vector(cfname + "/Cvec.txt", C_vec);
    load_vector(cfname + "/Caltvec.txt", C_alt_vec);
    load_vector(cfname + "/CO2vec.txt", C_O2_vec);
    load_vector(cfname + "/Svec.txt", S_vec);
    load_vector(cfname + "/Wvec.txt", W_vec);

    load_vector(cfname + "/Mvecx.txt", Mag_vec[0]);
    load_vector(cfname + "/Mvecy.txt", Mag_vec[1]);
    load_vector(cfname + "/Mvecz.txt", Mag_vec[2]);
    load_vector(cfname + "/Mvecx2.txt", Mag_vec2[0]);
    load_vector(cfname + "/Mvecy2.txt", Mag_vec2[1]);
    load_vector(cfname + "/Mvecz2.txt", Mag_vec2[2]);

    auto psi = itensor::MPS(sites);
    itensor::readFromFile(cfname + "/psi", psi);

    double curr_beta = 0;
    for (int i = 0; i != sec+1; i++){
        curr_beta += std::real(T[i]) * -2 * timesteps[i];
    }

    
    auto t1 = std::chrono::system_clock::now();

    for (int i = sec+1; i != timesteps.size(); i++){
        int curbond = tdvp_loop(E_vec,C_vec,C_alt_vec,C_O2_vec,S_vec,W_vec,Mag_vec,Mag_vec2,
            H0,psi,T[i],timesteps[i],args,sweeps,curr_beta,evnum,i);
        max_bond = std::max(max_bond,curbond);
    }
    S_vec = std::log(dims) * LX * LY - S_vec;

    {
        std::lock_guard<std::mutex> lock(forloop_mutex);
        Energies.push_back(E_vec);
        E_vec.clear();
        Capacity.push_back(C_vec);
        C_vec.clear();
        Alternative_Capacity.push_back(C_alt_vec);
        C_alt_vec.clear();
        Alternative_Capacity_O2.push_back(C_O2_vec);
        C_O2_vec.clear();
        Entropy.push_back(S_vec);
        S_vec.clear();
        Flux.push_back(W_vec);
        W_vec.clear();

        for (int j = 0; j != 3; j++){
            Magnetization[j].push_back(Mag_vec[j]);
            Magnetization2[j].push_back(Mag_vec2[j]);
            Mag_vec[j].clear();
            Mag_vec2[j].clear();
        }
    }


    auto t2 = std::chrono::system_clock::now();
    auto time = std::chrono::duration<double>(t2-t1);
    std::cout << "Finished Evolution ID " << std::this_thread::get_id() << ", Time Needed: " << time.count() << " Seconds\n" << std::flush;




}


template<typename T>
void Kitaev_Model::load_vector(std::string filename, std::vector<T>& v){
    std::cout << "File: " << filename << "\n" << std::flush;
    T ent;
    std::ifstream file(filename);
    
    while (file >> ent){
        v.push_back(ent);
    }
    std::cout << filename << " finished\n" << std::flush;
}




/*
    std::array<itensor::MPS*,3> Spsi;

    for (int i = 0; i != 3; i++){
        *(Spsi[i]) = itensor::applyMPO(M[i],psi);

        for (int j = 0; j != steps; j++){
            double E = itensor::tdvp(*(Spsi[i]),H0,dt,sweeps,args);
            std::cout << j << " " << std::flush;
        }
        double chip = itensor::inner(*(Spsi[i]),*(Spsi[i])) * n;
        double chiold = chi_vec[i].back();
        chi_vec[i].push_back(step/2 * (chip + chiold));
    }
    */






void Kitaev_Model::TPQ_MPS(std::vector<int> timesteps, std::vector<double> intervals, int Evols, int max_sites, int init_rand_sites, std::string TDVP_Type, double SusceptDiff, std::string Suscepts){
    this->CalcTDVP = true;
    this->SusceptIntegral = SusceptDiff != 0;

    if (SusceptIntegral){
        double hx = H_Details.get("hx");
        double hy = H_Details.get("hy");
        double hz = H_Details.get("hz");

        H_Details.set("hx",0);
        H_Details.set("hy",0);
        H_Details.set("hz",0);

        H_Details.set("hx",SusceptDiff);
        add_magnetic_interaction(aux,sec_aux);
        H0x = toMPO(ampo);

        H_Details.set("hx",-1.*SusceptDiff);
        H_Details.set("hy",SusceptDiff);
        add_magnetic_interaction(aux,sec_aux);
        H0y = toMPO(ampo);

        H_Details.set("hx",0);
        H_Details.set("hy",-1.*SusceptDiff);
        H_Details.set("hz",SusceptDiff);
        add_magnetic_interaction(aux,sec_aux);
        H0z = toMPO(ampo);

        H_Details.set("hy",0);
        H_Details.set("hz",-1*SusceptDiff);
        add_magnetic_interaction(aux,sec_aux);

        H_Details.set("hx",hx);
        H_Details.set("hy",hy);
        H_Details.set("hz",hz);
    }

    if (Evols < 1){
        std::invalid_argument("Please type in a valid argument for Evols");
    }

    Calsusx = Suscepts.find("x") != std::string::npos;
    Calsusy = Suscepts.find("y") != std::string::npos;
    Calsusz = Suscepts.find("z") != std::string::npos;

    int casintx = static_cast<int>(Calsusx);
    int casinty = static_cast<int>(Calsusy);
    int casintz = static_cast<int>(Calsusz);

    save_data(backup_name+"/chix",casintx);
    save_data(backup_name+"/chiy",casinty);
    save_data(backup_name+"/chiz",casintz);
    
    if (timesteps.size() != intervals.size()){
        std::invalid_argument("Time Steps vector and Intervals vector have to have the same length");
    }
    int entries = 1;
    for (int& i : timesteps){
        entries += i;
    }
    xdata.reserve(entries-1);

    Energies.reserve(Evols);
    Capacity.reserve(Evols);
    Alternative_Capacity.reserve(Evols);
    Entropy.reserve(Evols);
    Flux.reserve(Evols);


    for (auto& i : Magnetization){
        i.reserve(Evols);
    }
    for (auto& i : Magnetization2){
        i.reserve(Evols);
    }
    for (auto& i : Susceptibility){
        i.reserve(Evols);
    }



    std::cout << "Intervals: ";
    for (auto& i : intervals){
        std::cout << i << " ";
    }
    std::cout << "\n" << "Time Steps: ";
    for (auto& i : timesteps){
        std::cout << i << " ";
    }
    std::cout << "\n\n";
    std::cout << "Initial Random State Bond Dimension: " << init_rand_sites  << "\n";
    std::cout << "Maximum Bond Dimension: " << max_sites << "\n" << std::flush;


    auto Sweeps = itensor::Sweeps(1);
    itensor::Args tdvp_args;
    if (TDVP_Type == "TwoSite") {
        tdvp_args = itensor::Args({"Silent",true,"ErrGoal",1E-7});
        int tsave = 2;
        save_data(backup_name+"/TDVP",tsave);
    } else {
        tdvp_args = itensor::Args({"Silent",true,"ErrGoal",1E-7,"NumCenter",1});
        int tsave = 1;
        save_data(backup_name+"/TDVP",tsave);
    }
    Sweeps.maxdim() = max_sites;
    Sweeps.cutoff() = 1e-10;
    Sweeps.mindim() = init_rand_sites;
    Sweeps.niter() = 30;
    std::cout << "TDVP Technique used: " << TDVP_Type << "\n";

    std::vector<itensor::Cplx> T;
    std::vector<double> Tsave;
    double xsum = 0.;
    for (int i = 0; i != timesteps.size(); i++){
        double t_doub = -0.5 * intervals[i] / static_cast<double>(timesteps[i]);
        itensor::Cplx t = t_doub * itensor::Cplx_1;
        Tsave.push_back(t_doub);
        T.emplace_back(t);

        for (int j = 0; j != timesteps[i]; j++){
            double delb = std::real(t) * -2.;
            double xsumshift = xsum + delb / 2.;
            xsum += delb;
            xdata.push_back(1./xsum);
            xdatashift.push_back(1./xsumshift);
        }
    }

    save_data(backup_name+"/Evols",Evols);
    save_data(backup_name+"/T",Tsave);
    save_data(backup_name+"/timesteps",timesteps);
    save_data(backup_name+"/minsites",init_rand_sites);
    save_data(backup_name+"/maxsites",max_sites);
    save_data(backup_name+"/Susdiff",SusceptDiff);


    auto t0 = std::chrono::system_clock::now();
    int max_bond = 0;
    
    std::vector<std::future<void>> ThreadVector;
    ThreadVector.reserve(Evols-1);
    std::cout << "\n\n";
    for (int i = 0; i != Evols-1; i++){
        int evnum = i+1;
        ThreadVector.push_back(std::async(std::launch::async,[&](){time_evolution(T,timesteps,entries,SusceptDiff,init_rand_sites,std::ref(max_bond),tdvp_args,Sweeps,evnum);}));
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }


    time_evolution(T,timesteps,entries,SusceptDiff,init_rand_sites,max_bond,tdvp_args,Sweeps,0);


    for (auto& i : ThreadVector){
        i.wait();
    }


    E = Mean(Energies);
    Cv = Mean(Capacity);
    Cv_alt = Mean(Alternative_Capacity);
    Cv_O2 = Mean(Alternative_Capacity_O2);
    S = Mean(Entropy);
    W = Mean(Flux);

    Mx = Mean(Magnetization[0]);
    My = Mean(Magnetization[1]);
    Mz = Mean(Magnetization[2]);

    Mx2 = Mean(Magnetization2[0]);
    My2 = Mean(Magnetization2[1]);
    Mz2 = Mean(Magnetization2[2]);

    if (SusceptIntegral){
        if (Calsusx){
            Chix = Mean(Susceptibility[0]);
        }
        if (Calsusy){
            Chiy = Mean(Susceptibility[1]);
        }
        if (Calsusz){
            Chiz = Mean(Susceptibility[2]);
        }
    }

    auto t3 = std::chrono::system_clock::now();
    auto time_total = std::chrono::duration<double>(t3-t0);

    auto hours = std::chrono::duration_cast<std::chrono::hours>(time_total);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(time_total-hours);
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(time_total-hours-minutes);
    std::cout << "\n\n";
    std::cout << "Finished Imaginary Time Evolution, Time Needed: " << hours.count() << " Hours, " << minutes.count() << " Minutes, " << seconds.count() << " Seconds\n" << std::flush;

    if (max_bond == max_sites){
        std::cout << "Warning: MPS bond dimension has reached the bond dimension maximum. For more accurate results increase max_sites.\n\n" << std::flush;
    }
    else {
        std::cout << "Maximum MPS bond dimension: " << max_bond << "\n\n" << std::flush;
    }
}



Kitaev_Model::Kitaev_Model(std::string backup_name){

    this->backup_name = backup_name;

    std::ifstream file1(backup_name + "/LX.txt");
    file1 >> this->LX;
    file1.close();

    std::ifstream file2(backup_name + "/LY.txt");
    file2 >> this->LY;
    file2.close();

    int DoubleSpin;
    std::ifstream file3(backup_name + "/DS.txt");
    file3 >> DoubleSpin;
    file3.close();

    std::ifstream file4(backup_name + "/aux.txt");
    file4 >> this->aux;
    file4.close();

    std::ifstream file5(backup_name + "/sec_aux.txt");
    file5 >> this->sec_aux;
    file5.close();

    std::ifstream file6(backup_name + "/LT.txt");
    file6 >> this->Lattice_Type;
    file6.close();

    int Evols;
    std::ifstream file7(backup_name + "/Evols.txt");
    file7 >> Evols;
    file7.close();


    std::ifstream filestr(backup_name + "/foldername.txt");
    if (filestr){
        filestr >> this->Foldername;
        filestr.close();
    }

    this -> sites = itensor::readFromFile<itensor::CustomSpinNitsch>(backup_name + "/sites");
    this -> ampo = AutoMPO(sites);
    this -> dims = DoubleSpin + 1;
    this -> CalcTDVP = true;

    std::vector<double> Hinf;
    Hinf.reserve(9);
    load_vector(backup_name + "/H.txt", Hinf);

    H_Details.set("Kx",Hinf[0]);
    H_Details.set("Ky",Hinf[1]);
    H_Details.set("Kz",Hinf[2]);
    H_Details.set("J",Hinf[3]);
    H_Details.set("hx",Hinf[4]);
    H_Details.set("hy",Hinf[5]);
    H_Details.set("hz",Hinf[6]);
    H_Details.set("Gamma",Hinf[7]);
    H_Details.set("GammaQ",Hinf[8]);
    file7.close();


    std::vector<int> full_points;
    switch(Lattice_Type){
        case 1:
        case 3:
        case 4:
        case 5:
            full_points.reserve(((LY+1)/2)*LX);
            for (int m = 0; m != LX; m++){
                for (int n = 1; n <= LY; n+=2){
                    full_points.push_back(n+LY*m);
                }   
            }
            break;
        case 2:
        case 6:
            full_points.reserve(LX*LY);
            for (int i = 1; i != LX*LY+1; i++){
                full_points.push_back(i);
            }
            break;
    }



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
        
    this -> H0 = toMPO(this -> ampo);
    auto Ms = magnetization_operators(aux,sec_aux);
    M = Ms[0];
    M2 = Ms[1];

    std::cout << "Spin " << DoubleSpin << "/2 System" << "\n";
    std::cout << "Lattice Type Reloaded: " << Lattice_Type << "\n";
    std::cout << "LX = " << LX << ", LY = " << LY << "\n";
    std::cout << "Auxiliary Sites: " << aux << "\n";
    std::cout << "Secondary Auxiliary Sites: " << sec_aux << "\n\n";
    std::cout << "Hamiltonian Parameters \n";
    this->H_Details.print();
    std::cout << "\n\n";

    int max_bond = 0;
    auto t0 = std::chrono::system_clock::now();

    std::vector<itensor::Cplx> T;
    std::vector<int> timesteps;
    load_vector(backup_name + "/timesteps.txt", timesteps);

    double Tsmall;
    std::ifstream file8(backup_name + "/T.txt");
    while (file8 >> Tsmall){
        T.push_back(Tsmall * itensor::Cplx_1);
    }
    file8.close();

    int entries = 0;
    for (auto& i : timesteps){
        entries += i;
    }


    int SusceptDiff = 0;

    int init_rand_sites;
    std::ifstream file9(backup_name + "/minsites.txt");
    file9 >> init_rand_sites;
    file9.close();

    int max_sites;
    std::ifstream file10(backup_name + "/maxsites.txt");
    file10 >> max_sites;
    file10.close();

    int TDVP;
    std::ifstream file11(backup_name + "/TDVP.txt");
    file11 >> TDVP;
    file11.close();

    auto Sweeps = itensor::Sweeps(1);
    auto tdvp_args = itensor::Args({"Silent",true,"ErrGoal",1E-7,"NumCenter",TDVP});
    Sweeps.maxdim() = max_sites;
    Sweeps.cutoff() = 1e-10;
    Sweeps.mindim() = init_rand_sites;
    Sweeps.niter() = 30;

    /*
    for (int i = 0; i != timesteps.size(); i++){
        for (int j = 0; j != timesteps[i]; j++){

        }
    }
*/



    std::vector<std::future<void>> ThreadVector;
    ThreadVector.reserve(Evols-1);
    for (int i = 0; i != Evols-1; i++){
        int evnum = i+1;
        ThreadVector.push_back(std::async(std::launch::async,[&](){time_evolution_cont(T,timesteps,entries,SusceptDiff,init_rand_sites,std::ref(max_bond),tdvp_args,Sweeps,evnum);}));
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    time_evolution_cont(T,timesteps,entries,SusceptDiff,init_rand_sites,max_bond,tdvp_args,Sweeps,0);


    E = Mean(Energies);
    Cv = Mean(Capacity);
    Cv_alt = Mean(Alternative_Capacity);
    Cv_O2 = Mean(Alternative_Capacity_O2);
    S = Mean(Entropy);
    W = Mean(Flux);

    Mx = Mean(Magnetization[0]);
    My = Mean(Magnetization[1]);
    Mz = Mean(Magnetization[2]);

    Mx2 = Mean(Magnetization2[0]);
    My2 = Mean(Magnetization2[1]);
    Mz2 = Mean(Magnetization2[2]);
    
    
    auto t3 = std::chrono::system_clock::now();
    auto time_total = std::chrono::duration<double>(t3-t0);

    auto hours = std::chrono::duration_cast<std::chrono::hours>(time_total);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(time_total-hours);
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(time_total-hours-minutes);
    std::cout << "\n\n";
    std::cout << "Finished Imaginary Time Evolution, Time Needed: " << hours.count() << " Hours, " << minutes.count() << " Minutes, " << seconds.count() << " Seconds\n" << std::flush;

    if (max_bond == max_sites){
        std::cout << "Warning: MPS bond dimension has reached the bond dimension maximum. For more accurate results increase max_sites.\n\n" << std::flush;
    }
    else {
        std::cout << "Maximum MPS bond dimension: " << max_bond << "\n\n" << std::flush;
    }

    

}






}






