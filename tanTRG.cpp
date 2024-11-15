#include "TPQ-MPS.h"
#pragma once



using namespace itensor;
namespace TDVP_MPS{


MPS Kitaev_Model::mpo_to_tanmps(MPO& H){
    auto Htan = MPS(sitestan);
    //PrintData(Htan);

    auto H1 = H(1);
    H1 = permute(H1,{rightLinkIndex(H,1),siteInds(H,1)(1),siteInds(H,1)(2)});
    auto ind11 = siteIndex(Htan,1);
    auto indlinkold = sim(rightLinkIndex(H,1));
    auto ten1 = ITensor(ind11,indlinkold);
    //std::cout << "aaaaaaaaaaaaaaaaaaa\n" << std::flush;
    //println(H1);
    //PrintData(ten1);

    for (auto I : iterInds(H1)){
        Cplx value = eltC(H1,I);
        int newind = val(I[1]) + (val(I[2]) - 1) * dims;
        ten1.set(newind,val(I[0]),value);
    }
    Htan.set(1,ten1);
    //std::cout << "aaaaaaaaaaaaaaaaaaa\n" << std::flush;


    for (int i = 2; i != length(Htan); i++){
        auto Hiten = H(i);
        Hiten = permute(Hiten,leftLinkIndex(H,i),rightLinkIndex(H,i),siteInds(H,i)(1),siteInds(H,i)(2));
        auto ind2 = siteInds(Htan,i)(2);
        auto indlinknew = sim(rightLinkIndex(H,i));
        auto ten = ITensor(indlinkold,ind2,indlinknew);
        //PrintData(ten);
        //println(Hiten);

        for (auto I : iterInds(Hiten)){ 
            Cplx value = eltC(Hiten,I);
            int newind = val(I[2]) + (val(I[3]) - 1) * dims;
            ten.set(val(I[0]),newind,val(I[1]),value);
        }
        Htan.set(i,ten);
        indlinkold = indlinknew;
    }

    int last = length(Htan);
    auto Hlast = H(last);
    Hlast = permute(Hlast,{leftLinkIndex(H,last),siteInds(H,last)(1),siteInds(H,last)(2)});

    auto indl2 = siteInds(Htan,last)(2);
    auto tenlast = ITensor(indlinkold,indl2);

    for (auto I : iterInds(Hlast)){
        Cplx value = eltC(Hlast,I);
        int newind = val(I[1]) + (val(I[2]) - 1) * dims;
        tenlast.set(val(I[0]),newind,value);
    }
    Htan.set(last,tenlast);

    return Htan;
}


int adjusted_modulo(int A, int dims){
    int y = A % dims;
    if (y == 0){
        y = dims;
    }
    return y;
}



MPO Kitaev_Model::tanmps_to_mpo(MPS& Xtan){
    auto H = MPO(sites);

    auto X1 = Xtan(1);
    X1 = permute(X1,{siteIndex(Xtan,1),rightLinkIndex(Xtan,1)});
    auto indlinkold = sim(rightLinkIndex(Xtan,1));
    auto ind11 = siteInds(H,1)(1);
    auto ind12 = siteInds(H,1)(2);
    auto ten1 = ITensor(indlinkold,ind11,ind12);

    for (auto I : iterInds(X1)){
        Cplx value = eltC(X1,I);
        auto newind1 = adjusted_modulo(val(I[0]),dims);
        auto newind2 = (val(I[0]) - 1) / dims + 1;
        ten1.set(val(I[1]),newind1,newind2,value);
    }
    H.set(1,ten1);

    for (int i = 2; i != length(Xtan); i++){
        auto Xiten = Xtan(i);
        Xiten = permute(Xiten,{leftLinkIndex(Xtan,i),siteIndex(Xtan,i),rightLinkIndex(Xtan,i)});
        auto indlinknew = sim(rightLinkIndex(Xtan,i));
        auto indi1 = siteInds(H,i)(1);
        auto indi2 = siteInds(H,i)(2);
        auto teni = ITensor(indlinkold,indlinknew,indi1,indi2);

        for (auto I : iterInds(Xiten)){
            Cplx value = eltC(Xiten,I);
            int newind1 = adjusted_modulo(val(I[1]),dims);
            int newind2 = (val(I[1]) - 1) / dims + 1;
            teni.set(val(I[0]),val(I[2]),newind1,newind2,value);
        }
        H.set(i,teni);
        indlinkold = indlinknew;
    }

    int last = length(Xtan);
    auto Xlast = Xtan(last);
    Xlast = permute(Xlast,{leftLinkIndex(Xtan,last),siteIndex(Xtan,last)});
    auto indl1 = siteInds(H,last)(1);
    auto indl2 = siteInds(H,last)(2);
    auto tenl = ITensor(indlinkold,indl1,indl2);

    for(auto I : iterInds(Xlast)){
        Cplx value = eltC(Xlast,I);
        auto newind1 = adjusted_modulo(val(I[1]),dims);
        auto newind2 = (val(I[1]) - 1) / dims + 1;
        tenl.set(val(I[0]),newind1,newind2,value);
    }
    H.set(last,tenl);
    

    return H;


}










MPO Kitaev_Model::mpo_to_tanmpo(MPO& H){
    auto Htan = MPO(sitestan);

    auto H1 = H(1);
    H1 = permute(H1,{rightLinkIndex(H,1),siteInds(H,1)(1),siteInds(H,1)(2)});
    auto ind11 = rightLinkIndex(H,1);
    auto ind12 = siteInds(Htan,1)(1);
    auto ind13 = siteInds(Htan,1)(2);
    auto ten1 = ITensor(ind11,ind12,ind13);

    for (auto I : iterInds(H1)){
        Cplx value = eltC(H1,I);        
        for (int j = 0; j != dims; j++){
            int newind1 = val(I[1]) + j * dims;
            int newind2 = val(I[2]) + j * dims;
            ten1.set(val(I[0]),newind1,newind2,value);
        }
    }
    Htan.set(1,ten1);


    for (int i = 2; i != length(Htan); i++){
        auto Hiten = H(i);
        Hiten = permute(Hiten,leftLinkIndex(H,i),rightLinkIndex(H,i),siteInds(H,i)(1),siteInds(H,i)(2));
        auto ind1 = leftLinkIndex(H,i);
        auto ind2 = rightLinkIndex(H,i);
        auto ind3 = siteInds(Htan,i)(1);
        auto ind4 = siteInds(Htan,i)(2);
        auto teni = ITensor(ind1,ind2,ind3,ind4);

        for (auto I : iterInds(Hiten)){
            Cplx value = eltC(Hiten,I);
            for (int j = 0; j != dims; j++){
                int newind1 = val(I[2]) + j * dims;
                int newind2 = val(I[3]) + j * dims;
                teni.set(val(I[0]),val(I[1]),newind1,newind2,value);
            }
        }
        Htan.set(i,teni);

    }

    int last = length(Htan);
    auto Hlast = H(last);
    Hlast = permute(Hlast,{leftLinkIndex(H,last),siteInds(H,last)(1),siteInds(H,last)(2)});

    auto indl1 = leftLinkIndex(H,last);
    auto indl2 = siteInds(Htan,last)(1);
    auto indl3 = siteInds(Htan,last)(2);
    auto tenl = ITensor(indl1,indl2,indl3);

    for (auto I : iterInds(Hlast)){
        Cplx value = eltC(Hlast,I);
        for (int j = 0; j != dims; j++){
            int newind1 = val(I[1]) + j * dims;
            int newind2 = val(I[2]) + j * dims;
            tenl.set(val(I[0]),newind1,newind2,value);
        }
    }
    Htan.set(last,tenl);


    return Htan;   
    
} 







MPS Kitaev_Model::mps_to_tanmps(MPS& X){
    auto Xtan = MPS(sitestan);

    auto X1 = X(1);
    auto indsite1 = siteIndex(Xtan,1);
    auto indold = sim(rightLinkIndex(X,1));
    auto ten1 = ITensor(indsite1,indold);

    for (auto I : iterInds(X1)){
        for (auto J : iterInds(X1)){
            if (val(I[1]) == val(J[1])){   
                auto newind = val(I[0]) + (val(J[0]) - 1) * dims;
                auto value = eltC(X1,I) * eltC(X1,J);
                ten1.set(newind,val(I[1]),value); 
            }
        }
    }
    Xtan.set(1,ten1);

    for (int i = 2; i != length(Xtan); i++){
        auto Xiter = X(i);
        auto indsite = siteInds(Xtan,i)(2);
        auto indlinknew = sim(rightLinkIndex(X,i));
        auto ten = ITensor(indold,indsite,indlinknew);
        for (auto I : iterInds(Xiter)){
            for (auto J : iterInds(Xiter)){
                if (val(I[0]) == val(J[0]) && val(I[2]) == val(I[2])){
                    auto newind = val(I[1]) + (val(J[1]) - 1) * dims;
                    auto value = eltC(Xiter,I) * eltC(Xiter,J);
                    ten.set(val(I[0]),newind,val(I[2]),value);
                }
            }
        }
        Xtan.set(i,ten);
        indold = indlinknew;
    }

    auto last = length(Xtan);
    auto Xlast = X(last);
    auto indsite = siteInds(Xtan,last)(2);
    auto tenlast = ITensor(indold,indsite);
    for (auto I : iterInds(Xlast)){
        for (auto J : iterInds(Xlast)){
            if (val(I[0]) == val(I[0])){
                auto newind = val(I[1]) + (val(J[1]) - 1) * dims;
                auto value = eltC(Xlast,I) * eltC(Xlast,J);
                tenlast.set(val(I[0]),newind,value);
            }
        }
    }
    Xtan.set(last,tenlast);
    return Xtan;


}


std::vector<std::array<double,2>> Kitaev_Model::mean_wrap(std::vector<double>& vec){
    std::vector<double> cutvec(vec.begin()+1,vec.end());
    std::vector<std::array<double,2>> wrap;
    wrap.reserve(cutvec.size());

    for (auto& i : cutvec){
            std::array<double,2> part{i,0.};
            wrap.push_back(part);
    }
    return wrap;
}





int Kitaev_Model::tan_tdvp_loop(int steps, double dt, MPS& Hexptan, MPO& H0tan, MPO& H2tan,
std::array<MPO,3>& Mtan, std::array<MPO,3>& M2tan, Sweeps& sweeps, Args& tdvp_args,
std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& C_alt_vec, 
std::vector<double>& S_vec, std::vector<double>& W_vec,
std::array<std::vector<double>,3>& M_vec, std::array<std::vector<double>,3>& M_vec2,
MPO& Hfluxtan, double& cb, int KrylovExpansions){

    int max_bond = 0;
    for (int i = 0; i != steps; i++){
        if (cb <= static_cast<double>(KrylovExpansions)*dt){
            std::vector<double> eps = {1e-10,1e-10};
            addBasis(Hexptan,H0tan,eps,{"Cutoff",1E-8,
                                        "Method","DensityMatrix",
                                        "KrylovOrd",3,
                                        "DoNormalize",true,
                                        "Quiet",true});
        }


        auto stupid_result = tdvp(Hexptan,H0tan,-0.5*dt*Cplx_1,sweeps,tdvp_args);

        cb += dt;
        double cure = std::real(innerC(Hexptan,H0tan,Hexptan));
        double cure2 = std::real(innerC(Hexptan,H2tan,Hexptan));
        double curc = cb * cb * std::real(cure2 - cure * cure);
        double curc_alt = cb * cb * (E_vec.back() - cure) / dt;
        double curs = S_vec.back() + dt * 0.5 * (curc/cb + C_vec.back()/(cb-dt));
        double curw = std::real(innerC(Hexptan,Hfluxtan,Hexptan));

        E_vec.push_back(cure);
        C_vec.push_back(curc);
        C_alt_vec.push_back(curc_alt);
        S_vec.push_back(curs);
        W_vec.push_back(curw);

        for (int indi = 0; indi != 3; indi++){
            std::complex<double> m = innerC(Hexptan,Mtan[indi],Hexptan);
            std::complex<double> m2 = innerC(Hexptan,M2tan[indi],Hexptan);
            M_vec[indi].push_back(std::real(m));
            M_vec2[indi].push_back(std::real(m2));
        }

        max_bond = std::max(max_bond,maxLinkDim(Hexptan));
    }
    
    return max_bond;
}



double Kitaev_Model::calculate_chi(MPS& Hexptan, MPS& Hexptaninv, std::string spin){
    MPO Hexpmpo = tanmps_to_mpo(Hexptan);
    double curchi = 0.;

    for (int i = 1; i <= length(Hexptan); i++){
        AutoMPO chiampo(sites);
        chiampo += spin,i;
        MPO Spin1 = toMPO(chiampo);

        MPO chimpo = nmultMPO(Spin1,prime(Hexpmpo));
        chimpo = nmultMPO(chimpo,prime(Hexpmpo));

        for (int j = 1; j <= length(Hexptan); j++){      
            AutoMPO chiampo(sites);
            chiampo += spin,j;
            MPO Spin2 = toMPO(chiampo);

            MPO chimpo_total = nmultMPO(chimpo,prime(Spin2));
            MPO chimpo_prod = mpo_to_tanmpo(chimpo_total);
            curchi += std::real(innerC(Hexptaninv,chimpo_prod,Hexptaninv));
        }
    }

    return curchi;
}







void Kitaev_Model::tanTRG(std::vector<int> timesteps, std::vector<double> intervals, int max_sites, int KrylovExpansions){
    auto t0 = std::chrono::system_clock::now();
    this -> CalcTDVP = true;

    MPO H2 = nmultMPO(H0,prime(H0));
    MPO H2tan = mpo_to_tanmpo(H2);

    auto Hexp1 = toExpH(ampo,0.25*intervals[0]/static_cast<double>(timesteps[0])*(Cplx_1+Cplx_i));
    auto Hexp2 = toExpH(ampo,0.25*intervals[0]/static_cast<double>(timesteps[0])*(Cplx_1-Cplx_i));
    auto Hexp = nmultMPO(Hexp2,prime(Hexp1));
    auto Hexptan = mpo_to_tanmps(Hexp);
    Hexptan.orthogonalize();
    Hexptan.normalize();

    auto H0tan = mpo_to_tanmpo(H0);
    auto Hfluxtan = mpo_to_tanmpo(H_flux);

    std::array<MPO,3> Mtan, M2tan;
    for (int i = 0; i != 3; i++){
        Mtan[i] = mpo_to_tanmpo(M[i]);
        M2tan[i] = mpo_to_tanmpo(M2[i]);
    }

    std::cout << "Intervals: ";
    for (auto& i : intervals){
        std::cout << i << " ";
    }
    std::cout << "\n\n";
    std::cout << "Time Steps: ";
    for (auto& i : timesteps){
        std::cout << i << " ";
    }
    std::cout << "\n\n";
    std::cout << "Maximum Sites: " << max_sites << "\n";
    std::cout << "Number of Krylov Expansions: " << KrylovExpansions << "\n" << std::flush;
    auto sweeps = Sweeps(1);
    sweeps.maxdim() = max_sites;
    sweeps.cutoff() = 1e-8;
    sweeps.niter() = 30;

    auto tdvp_args = Args({"Silent",true,"ErrGoal",1E-7});

    int entries = 1;
    for (int& i : timesteps){
        entries += i;
    }
    xdata.reserve(entries-1);
    
    std::vector<double> E_vec;
    E_vec.reserve(entries);

    std::vector<double> C_vec;
    C_vec.reserve(entries);

    std::vector<double> C_alt_vec;
    C_alt_vec.reserve(entries);

    std::vector<double> S_vec;
    S_vec.reserve(entries);

    std::vector<double> W_vec;
    W_vec.reserve(entries);

    std::array<std::vector<double>,3> Mag_vec;
    std::array<std::vector<double>,3> Mag_vec2;


    E_vec.push_back(0.);
    E_vec.push_back(std::real(innerC(Hexptan,H0tan,Hexptan)));
    W_vec.push_back(0.);
    W_vec.push_back(std::real(innerC(Hexptan,Hfluxtan,Hexptan)));

    std::vector<double> dT;
    double xsum = 0.;
    for (int i = 0; i != timesteps.size(); i++){
        double t = intervals[i] / static_cast<double>(timesteps[i]);
        dT.emplace_back(t);

        for (int j = 0; j != timesteps[i]; j++){
            double delb = t;
            xsum += delb;
            xdata.push_back(1./xsum);
        }
    }

    double cur_beta = dT[0];
    
    C_vec.push_back(0.);
    C_vec.push_back(cur_beta * (E_vec[1] - E_vec[0]));
    C_alt_vec.push_back(0.);
    C_alt_vec.push_back(cur_beta * (E_vec[1] - E_vec[0]));
    S_vec.push_back(0.);
    S_vec.push_back(dT[0] * 0.5 * C_vec[1] / cur_beta);

    for (int indi = 0; indi != 3; indi++){
        std::complex<double> m = innerC(Hexptan,Mtan[indi],Hexptan);
        std::complex<double> m2 = innerC(Hexptan,M2tan[indi],Hexptan);
        Mag_vec[indi].reserve(entries);
        Mag_vec2[indi].reserve(entries);
        Mag_vec[indi].push_back(std::real(m));
        Mag_vec2[indi].push_back(std::real(m2));
    }


    int max_bond = tan_tdvp_loop(timesteps[0]-1,dT[0],Hexptan,H0tan,H2tan,Mtan,M2tan,sweeps,tdvp_args,
                                E_vec,C_vec,C_alt_vec,S_vec,W_vec,Mag_vec,Mag_vec2,
                                Hfluxtan,cur_beta,KrylovExpansions);

    for (auto i = 1; i != intervals.size(); i++){
        int m2 = tan_tdvp_loop(timesteps[i],dT[i],Hexptan,H0tan,H2tan,Mtan,M2tan,sweeps,tdvp_args,
                                E_vec,C_vec,C_alt_vec,S_vec,W_vec,Mag_vec,Mag_vec2,
                                Hfluxtan,cur_beta,KrylovExpansions);
        max_bond = std::max(max_bond,m2);
    }

    E = mean_wrap(E_vec);
    Cv = mean_wrap(C_vec);
    Cv_alt = mean_wrap(C_alt_vec);
    S = mean_wrap(S_vec);
    W = mean_wrap(W_vec);

    Mx = mean_wrap(Mag_vec[0]);
    My = mean_wrap(Mag_vec[1]);
    Mz = mean_wrap(Mag_vec[2]);
    Mx2 = mean_wrap(Mag_vec2[0]);
    My2 = mean_wrap(Mag_vec2[1]);
    Mz2 = mean_wrap(Mag_vec2[2]);

    auto tfin = std::chrono::system_clock::now();
    auto time = std::chrono::duration<double>(tfin-t0);

    auto hours = std::chrono::duration_cast<std::chrono::hours>(time);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(time-hours);
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(time-hours-minutes);

    std::cout << "Finished Simulation, Time Needed: " << hours.count() << " Hours, " << minutes.count() << " Minutes, " << seconds.count() << " Seconds\n\n" << std::flush;
    if (max_bond == max_sites){
        std::cout << "Warning: MPS bond dimension has reached the bond dimension maximum. For more accurate results increase max_sites.\n\n" << std::flush;
    }
    else {
        std::cout << "Maximum MPS bond dimension: " << max_bond << "\n\n" << std::flush;
    }
}






}



