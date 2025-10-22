/* Ceci est une implementation pour un element Pk */

/*#include "ff++.hpp"
#include "AddNewFE.h"
using namespace std; 

class TypeOfFE_P4Lagrange : public TypeOfFE {
public:
    static const int k = 3;
    static const int ndf = (k + 2) * (k + 1) / 2;
    vector<int> Data;
    vector<double> Pi_h_coef;
    vector<vector<int>> nn(ndf,vector<int>(k,0))//[ndf][k];
    vector<vector<int>> aa(ndf,vector<int>(k,0))//[ndf][k];
    vector<int> ff(ndf,0);//[ndf];
    vector<int> il(ndf,0);//[ndf];
    vector<int> jl(ndf,0);//[ndf];
    vector<int> kl(ndf,0); //[ndf];
    //void FillDataLagrange( );
    TypeOfFE_P4Lagrange( ) : TypeOfFE(3 + 2 * 3 + 1, 1, Data, 4, 1, 16, 10, 0) {
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
#include<map>
using namespace std; 
struct R2
{
    double i, j;
    R2(double ii, double jj){
        i=ii;
        j=jj;
    }
};

// Constuction des coordonnees de points (TO DO: ajouter les points de l'interieur PB? Reste à savoir comment numeroter par convention)
vector<R2> PtConstruction(int k){
    int NptPerV = k-1;
    vector<R2> Pt;//(Ndof);
    // sommets
    Pt.push_back(R2(0,0));
    Pt.push_back(R2(1,0));
    Pt.push_back(R2(0,1));
    // arrete bas 
    for(int i =0; i<NptPerV;i++) Pt.push_back(R2( double(k-i-1),double(i+1)));
    // arrete droite
    for(int i =0; i<NptPerV;i++) Pt.push_back(R2(0,double(k-i-1)));
    // arrete gauche
    for(int i =0; i<NptPerV;i++) Pt.push_back(R2(double(i+1),0.)); 
    // Interieur de l'element

    return Pt;
}

// Remplir le tableau data
void FillDataLagrange(int k, vector<int> &Data, vector<int> &Pi_h_coef ){
    int ndof= (k+1)*(k+2)*0.5;
    Data.resize(5*ndof+3,0);
    // First row(the support number  of the node of the df)
    //Data[0]=0;
    Data[1]=1;
    Data[2]=2;
    // Second row(the number of the df on  the node)
    /*Data[ndof]=0;
    Data[ndof+1]=0;
    Data[ndof+2]=0;*/
    //3rd row
    Data[2*ndof]=0;
    Data[2*ndof+1]=1;
    Data[2*ndof+2]=2;
    //4th row
    /*Data[3*ndof]=0;
    Data[3*ndof+1]=0;
    Data[3*ndof+2]=0;*/
    //5th row
    //Data[4*ndof]=0;
    Data[4*ndof+1]=1;
    Data[4*ndof+2]=2;
    for (int i =3; i<ndof; i++){
        if(i<=k+1){
            Data[i]=3;
            Data[ndof+i]=i-3;
            
        }
        else{
            if(i<=2*k) {
                Data[i]=4;
                Data[ndof+i]=i-3-(k)+1; 
            }
            else if (i<=3*k-1) {
                Data[i]=5; 
                Data[ndof+i]=i-3-2*(k-1);
            }
            else {
                Data[i]=6;
                Data[ndof+i]=i-3-3*(k-1);
            } 
        }
        // 3row
        Data[2*ndof+i]=Data[i];
        // 4 th row
        //Data[3*ndof+i]=0;
        //5th row
        Data[4*ndof+i]=i;
    }
    Data[5*ndof+2]=ndof;
    /*Data[5*ndof+1]=0;
    Data[5*ndof]=0;*/
    Pi_h_coef.resize(ndof,1.);

}
// if orientation is <-1 it provides the indices of nodes to be exchanged  
vector<pair<int,int>> Exchangeidx(int k){
    vector<pair<int,int>> permut;
    permut.resize(3);

    permut[0].first=3;
    permut[0].second=3+(k-2);
    
    permut[1].first=3+k-1;
    permut[1].second=3+(2*k-3);
    
    permut[2].first=3+2*k-2;
    permut[2].second=3+(3*k-4);    
    return permut;
}

// Numerotation point interieur
void NumInternPoint(int k ){
    // l'indice i associe a la coord lambda3 et j a lambda2
    // donc celui de lambda1 va etre k-(i+j) on itere uniquement
    // entre 1 et k-2 pour eviter les noeuds des aretes  
    for (int i =1;i<=k-2;i++){
        for (int j =1;j<=k-2-i+1;j++){
            cout<<"("<<k-(i+j)<<","<<j<<")";
        }
    }
}
int main(int argv, char **argc){
    int PK= std::atoi(argc[1]);
    int Ndof=(PK+1)*(PK+2)*0.5;
    vector<R2> Test= PtConstruction(PK);
    //for(auto &elem :Test ) cout<<"elem:"<<elem.i<<","<<elem.j<<endl;
    vector<int> Data;
    vector<int> Pi_dataLagrange;

    FillDataLagrange(PK, Data,Pi_dataLagrange );
    /*cout<<"Row 1:"<<endl;
    for(int i=0; i<Ndof;i++ ) cout<<Data[i]<<", ";
    cout<<"\nRow 2:"<<endl;
    for(int i=0; i<Ndof;i++ ) cout<<Data[Ndof+i]<<", ";
    cout<<"\nRow 3:"<<endl;
    for(int i=0; i<Ndof;i++ ) cout<<Data[2*Ndof+i]<<", ";
    cout<<"\nRow 4:"<<endl;
    for(int i=0; i<Ndof;i++ ) cout<<Data[3*Ndof+i]<<", ";
    cout<<"\nRow 5:"<<endl;
    for(int i=0; i<Ndof;i++ ) cout<<Data[4*Ndof+i]<<", ";
    cout<<"\nRow 6:"<<endl;
     cout<<Data[5*Ndof];
    cout<<"\nRow 7:"<<endl;
     cout<<Data[5*Ndof+1]<<" "<<Data[5*Ndof+2]<<endl;*/
    //for (auto elem :Pi_dataLagrange ) cout<<" "<<elem<<endl;
    //vector<pair<int,int>> Idxs=Exchangeidx( PK);
    //for (auto elm : Idxs ) cout<<"("<<elm.first<<","<<elm.second<<")";
    NumInternPoint(PK );
    cout<<endl;
    return 0; 
}