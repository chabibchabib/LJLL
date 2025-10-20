/* Ceci est une implementation pour un element Pk */

/*#include "ff++.hpp"
#include "AddNewFE.h"
using namespace std; 

class TypeOfFE_P3Lagrange : public TypeOfFE {
public:
    static const int k = 3;
    static const int ndf = (k + 2) * (k + 1) / 2;
    static int Data[];
    static double Pi_h_coef[];
    static const int nn[ndf][k];
    static const int aa[ndf][k];
    static const int ff[ndf];
    static const int il[ndf];
    static const int jl[ndf];
    static const int kl[ndf];
    
    TypeOfFE_P3Lagrange( ) : TypeOfFE(3 + 2 * 3 + 1, 1, Data, 4, 1, 16, 10, 0) {
        static const R2 Pt[10] = {R2(0 / 3., 0 / 3.), R2(3 / 3., 0 / 3.), R2(0 / 3., 3 / 3.),
            R2(2 / 3., 1 / 3.), R2(1 / 3., 2 / 3.), R2(0 / 3., 2 / 3.),
            R2(0 / 3., 1 / 3.), R2(1 / 3., 0 / 3.), R2(2 / 3., 0 / 3.),
            R2(1 / 3., 1 / 3.)};
        // 3,4,5,6,7,8
        int other[10] = {-1, -1, -1, 4, 3, 6, 5, 8, 7, -1};
        int kk = 0;
        
        for (int i = 0; i < NbDoF; i++) {
            pij_alpha[kk++] = IPJ(i, i, 0);
            if (other[i] >= 0) {
                pij_alpha[kk++] = IPJ(i, other[i], 0);
            }
            
            P_Pi_h[i] = Pt[i];
        }
        
        assert(P_Pi_h.N( ) == NbDoF);
        assert(pij_alpha.N( ) == kk);
    }
};*/
#include <vector>
#include<iostream>
using namespace std; 
struct R2
{
    double i, j;
    R2(double ii, double jj){
        i=ii;
        j=jj;
    }
};

 vector<R2> PtConstruction(int k){
    int NptPerV = k-1;
    vector<R2> Pt;//(Ndof);
    Pt.push_back(R2(0,0));
    Pt.push_back(R2(1,0));
    Pt.push_back(R2(0,1));
    // arrete ou y est fixe(celui du bas)
    for(int i =0; i<NptPerV;i++) Pt.push_back(R2(double(i+1),0));
    // arrete ou lambda1 est fixe(celui de droite) 
    for(int i =0; i<NptPerV;i++) Pt.push_back(R2( double(k-i-1),double(i+1)));
    // arrete ou x est fixe(celui de gauche)
    for(int i =0; i<NptPerV;i++) Pt.push_back(R2(0.,double(k-i-1)));      
    return Pt;
}
int main(){
    vector<R2> Test= PtConstruction(4);
    for(auto &elem :Test ) cout<<"elem:"<<elem.i<<","<<elem.j<<endl;
    return 0; 
}