#pragma once

#include "itensor/mps/siteset.h"
#include <vector>
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"



namespace itensor {

    



class CustomSpinNitschSite;

using CustomSpinNitsch = BasicSiteSet<CustomSpinNitschSite>;

class CustomSpinNitschSite
    {
    Index s;
    static bool is_called;
    static std::array<std::vector<std::vector<Cplx>>,3> SpinExp;


    static std::array<std::vector<std::vector<Cplx>>,3> ExpMatrix(int DoubleSpin){
        int dims = DoubleSpin + 1;
        double DSdouble = static_cast<double>(DoubleSpin);
        Eigen::MatrixX<Cplx> Sx = Eigen::MatrixX<Cplx>::Zero(dims,dims);
        Eigen::MatrixX<Cplx> Sy = Eigen::MatrixX<Cplx>::Zero(dims,dims);
        Eigen::MatrixX<Cplx> Sz = Eigen::MatrixX<Cplx>::Zero(dims,dims);

        for (int i = 0; i != dims-1; i++){
            double sz = static_cast<double>(i) - DSdouble / 2.; 
            double f = std::sqrt(DSdouble/2.*(DSdouble/2.+1.) - sz*(sz+1.));
            
            Cplx X = f/2. * Cplx_1;
            Cplx Y = -f/2. * Cplx_i;

            Sx(i+1,i) = X * M_PI * Cplx_i;
            Sx(i,i+1) = X * M_PI * Cplx_i;
            Sy(i+1,i) = Y * M_PI * Cplx_i;
            Sy(i,i+1) = -Y * M_PI * Cplx_i;
        }
        for (int i = 0; i != dims; i++){
            double sz = static_cast<double>(i) - DSdouble / 2.;
            Sz(i,i) = sz * M_PI * Cplx_i;
        }
        std::array<Eigen::MatrixX<Cplx>,3> ExpEig;
        ExpEig[0] = Sx.exp();
        ExpEig[1] = Sy.exp();
        ExpEig[2] = Sz.exp();

        std::array<std::vector<std::vector<Cplx>>,3> ExpVec;
        for (int i = 0; i < 3; i++) {
            ExpVec[i].resize(dims);
            for (int j = 0; j < dims; j++) {
                ExpVec[i][j].resize(dims);
            }
        }
        for (int i = 0; i != 3; i++){
            for (int j = 0; j != dims; j++){
                for (int k = 0; k != dims; k++){
                    ExpVec[i][j][k] = ExpEig[i](j,k);
                }
            }
        }     

        return ExpVec;
    }

        

    public:

    CustomSpinNitschSite() { }

    CustomSpinNitschSite(Index I) : s(I) { }
    
    CustomSpinNitschSite(Args const& args = Args::global())
        {
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserveSz = args.getBool("ConserveSz",conserveQNs);

        auto DSmax = 3;
        if(args.defined("2S")) 
            {
            DSmax = args.getInt("2S");
            }
        else if(args.defined("S")) 
            {
            auto S = args.getReal("S");
            DSmax = std::round(2*S+0.01);
            }
        else
            {
            error("Must pass named args \"2S\" (integer) or \"S\" (real number) to CustomSpinNitsch");
            }

        if(DSmax < 1) error(tinyformat::format("Invalid spin value %d/2 in CustomSpinNitsch",DSmax));
        
        if (is_called == false){
            CustomSpinNitschSite::SpinExp = CustomSpinNitschSite::ExpMatrix(DSmax);
            CustomSpinNitschSite::is_called = true;
        }

        auto tags = TagSet(tinyformat::format("Site,S=%d/2",DSmax));
        if(args.defined("SiteNumber") )
            {
            auto n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
            }

        if(conserveQNs)
            {
            if(conserveSz)
                {
                auto qints = Index::qnstorage(DSmax+1);
                for(int n : range(DSmax+1)) 
                    {
                    qints[n] = QNInt(QN({"Sz",static_cast<QNum::qn_t>(2*n-DSmax)}),1);
                    }
                s = Index(std::move(qints),tags);
                }
            else
                {
                s = Index(QN(),DSmax+1,tags);
                }
            }
        else
            {
            if(conserveSz) throw ITError("ConserveSz cannot be true when ConserveQNs=false");
            s = Index(DSmax+1,tags);
            }
        }


    Index
    index() const { return s; }
    
    IndexVal
    state(std::string const& state)
        {
        auto DSmax = dim(s)-1;
        for(auto n : range(DSmax+1))
            {
            if(state == str(2*n-DSmax)) return s(1+n); //state name is "2Sz"
            }
        throw ITError("State " + state + " not recognized");
        return IndexVal{};
        }


	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        auto Op = ITensor(dag(s),sP);
        
        auto DSmax = dim(s)-1;

        if(opname == "Sz")
            {
            for(int i = 1; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i),sz);
                }
            //std::cout << "Sz: \n";
            }
        else
        if(opname == "S+")
            {
            for(int i = 1; i <= dim(s)-1; ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i+1),std::sqrt((DSmax/2.0-sz)*(DSmax/2.0+sz+1.0)));
                }
            }
        else
        if(opname == "S-")
            {
            for(int i = 2; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i-1),std::sqrt((DSmax/2.0+sz)*(DSmax/2.0-sz+1.0)));
                }
            }








	
	else
	if(opname == "Sx")
	    {
	    for (int i=1; i < dim(s); ++i){
		auto sz = -DSmax/2. + i - 1;
		Op.set(s(i),sP(i+1),0.5*std::sqrt((DSmax/2.0-sz)*(DSmax/2.0+sz+1.0)));
		}
	    for(int i = 2; i <= dim(s); ++i)
            {
            auto sz = -DSmax/2.0+i-1;
            Op.set(s(i),sP(i-1),0.5*std::sqrt((DSmax/2.0+sz)*(DSmax/2.0-sz+1.0)));
            }
        //std::cout << "Sx: \n";
	    }
	else
	if(opname == "Sy"){
	    for(int i = 1; i <= dim(s)-1; ++i)
            {
            auto sz = -DSmax/2.0+i-1;
            Op.set(s(i),sP(i+1),-0.5*std::sqrt((DSmax/2.0-sz)*(DSmax/2.0+sz+1.0))*Cplx_i);
            }
        for(int i = 2; i <= dim(s); ++i)
            {
            auto sz = -DSmax/2.0+i-1;
            Op.set(s(i),sP(i-1),0.5*std::sqrt((DSmax/2.0+sz)*(DSmax/2.0-sz+1.0))*Cplx_i);
            }
        //std::cout << "Sy: \n";
	    }
    else
    if(opname == "Expx"){
        auto ExMat = SpinExp[0];
        for (int i = 0; i != dim(s); i++){
            for (int j = 0; j != dim(s); j++){
                Cplx Ex = ExMat[i][j];
                Op.set(s(j+1),sP(i+1),Ex);
            }
        }
        //std::cout << "Expx: \n";
    }
    else
    if(opname == "Expy"){
        auto ExMat = SpinExp[1];
        for (int i = 0; i != dim(s); i++){
            for (int j = 0; j != dim(s); j++){
                Cplx Ex = ExMat[i][j];
                Op.set(s(j+1),sP(i+1),Ex);
            }
        }
        //std::cout << "Expy: \n";
    }
    else
    if(opname == "Expz"){
        auto ExMat = SpinExp[2];
        for (int i = 0; i != dim(s); i++){
            for (int j = 0; j != dim(s); j++){
                Cplx Ex = ExMat[i][j];
                Op.set(s(j+1),sP(i+1),Ex);        
            }
        }
        //std::cout << "Expz \n";
    }

	




	
	else
        if(opname == "U+")
            {
            for(int i = 1; i <= dim(s)-1; ++i)
                {
                Op.set(s(i),sP(i+1),1.0);
                }
            }
        else
        if(opname == "U-")
            {
            for(int i = 2; i <= dim(s); ++i)
                {
                Op.set(s(i),sP(i-1),1.0);
                }
            }
        else
        if(opname == "Sz2")
            {
            for(int i = 1; i <= dim(s); ++i)
                {
                auto sz = -DSmax/2.0+i-1;
                Op.set(s(i),sP(i),pow(sz, 2.0));
                }
            }
        else
            {
            throw ITError("Operator " + opname + " name not recognized");
            }
        //PrintData(Op);

        return Op;
        }
    };



std::array<std::vector<std::vector<Cplx>>,3> CustomSpinNitschSite::SpinExp;
bool CustomSpinNitschSite::is_called = false;





}
