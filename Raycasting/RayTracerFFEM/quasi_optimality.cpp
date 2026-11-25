#include<iostream>
#include<cmath>
#include<vector>

typedef double REAL ;

std::vector<REAL> k_indepenent_qo(int p, REAL k, REAL trappingorder, REAL visibleorder){
    std::vector<REAL> list(4,0);
    list[0]= std::pow(k, -(1+trappingorder/p)); 
    list[1]= std::pow(k,- (1.+(1+visibleorder)/(2.*p) )); 
    list[2]= std::pow(k,-1-1./p); 
    list[3]= 1./k; 
    
    return list;
}

std::vector<REAL> k_indepenent_qo_awayfromtrapping(int p, REAL k, REAL trappingorder, REAL visibleorder){
    std::vector<REAL> list(4,0);
    list[0]= std::pow(k, -(1+(trappingorder+1.)/(2.*p) ) ); 
    list[1]= std::pow(k,-(p+1.)/p); 
    list[2]= std::pow(k,-(p+1.)/p); 
    list[3]= 1./k; 
    
    return list;
}

std::vector<REAL> errorctrl(int p, REAL k, REAL trappingorder, REAL visibleorder){
    std::vector<REAL> list(4,0);
    list[0]= std::pow(k, -(1.+(trappingorder)/(2.*p) ) ); 
    list[1]= std::pow(k,-(1+(1+visibleorder)/(4.*p))); 
    list[2]= std::pow(k,-(1.+1./(2.*p))); 
    list[3]= 1./k; 
    
    return list;
}

std::vector<REAL> errorctrl_awayfromtrapping(int p, REAL k, REAL trappingorder, REAL visibleorder){
    std::vector<REAL> list(4,0);
    list[0]= std::pow(k, -(1.+(trappingorder)/(2.*p) ) ); 
    list[1]= std::pow(k,-(1.+1./(2.*p))); 
    list[2]= std::pow(k,-(1.+1./(2.*p))); 
    list[3]= 1./k; 
    
    return list;
}
int main(int argc,char ** argv){
    //REAL eps= 4e-2;
    //REAL k=std::atof(argv[1]);
    REAL k=std::atof(argv[1]);

    int p = std::atoi(argv[2]);
    std::cout<<" p= "<<p<<" wavenumber k= "<<k<<std::endl;
    REAL trappingorder=2,  visibleorder=2.;
    auto seuil = k_indepenent_qo(  p,  k,  trappingorder,  visibleorder);
    std::cout<<" k-independent QO:"<<std::endl;
    for(auto elm : seuil ) std::cout<<elm<<" ";
    std::cout<<std::endl;

    auto seuiawaytra = k_indepenent_qo_awayfromtrapping(  p,  k,  trappingorder,  visibleorder);
    std::cout<<" k-independent QO away from trapping:"<<std::endl;
    for(auto elm : seuiawaytra ) std::cout<<elm<<" ";
    std::cout<<std::endl;
    auto seuilerror = errorctrl(  p,  k,  trappingorder,  visibleorder);
    std::cout<<" error ctrl:"<<std::endl;
    for(auto elm : seuilerror ) std::cout<<elm<<" ";
    std::cout<<std::endl;

    auto seuilerrorawaytra = errorctrl_awayfromtrapping(  p,  k,  trappingorder,  visibleorder);
    std::cout<<" error ctrl away from trapping:"<<std::endl;
    for(auto elm : seuilerrorawaytra ) std::cout<<elm<<" ";
    std::cout<<std::endl;
    return 0;
}