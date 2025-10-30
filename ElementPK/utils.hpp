#include <vector>
#include<iostream>
#include<map>
#include<cmath>
#include "ff++.hpp"
#include "AddNewFE.h"
using namespace std; 
/*struct R2
{
    double i, j;
    R2(double ii, double jj){
        i=ii;
        j=jj;
    }
};*/

// Numerotation point interieur
void NumInternPoint(int k,vector<R2> &Pt ){
    // l'indice i associe a la coord lambda3 et j a lambda2
    // donc celui de lambda1 va etre k-(i+j) on itere uniquement
    // entre 1 et k-2 pour eviter les noeuds des aretes  
    for (int i =1;i<=k-2;i++){
        for (int j =1;j<=k-2-i+1;j++){
            //Pt.push_back(R2(double(k-(i+j)),double(j)));
            Pt.push_back(R2(double(j)/k,double(i)/k));

        }
    }
}

// Numerotation point interieur comme P4
void NumInternPoint2(int k,vector<R2> &Pt ){
 
    /*for (int i =1;i<=k-2;i++){
        for (int j =1;j<=k-2-i+1;j++){
            //Pt.push_back(R2(double(k-(i+j)),double(j)));
            Pt.push_back(R2(double(j)/k,double(k-2-i+1-(j-1))/k));
            //cout<<"("<<j<<","<<k-2-i+1-(j-1)<<")";
        }
    }*/

    for (int i =k-2;i>=1;i--){
            for (int j =1;j<=k-2-i+1;j++){
                Pt.push_back(R2(double(k-1-i-(j-1))/k,double(i)/k));
                cout<<"("<<j<<","<<k-1-i-(j-1)<<","<<i<<")";
            }
        }
}

// Constuction des coordonnees de points (TO DO: ajouter les points de l'interieur PB? Reste à savoir comment numeroter par convention)
vector<R2> PtConstruction(int k){
    int NptPerV = k-1;
    vector<R2> Pt;//(Ndof);
    // sommets
    Pt.push_back(R2(0,0));
    Pt.push_back(R2(k/k,0));
    Pt.push_back(R2(0,k/k));
    // arrete bas 
    for(int i =0; i<NptPerV;i++) Pt.push_back(R2( double(k-i-1)/k,double(i+1)/k));
    // arrete droite
    for(int i =0; i<NptPerV;i++) Pt.push_back(R2(0/k,double(k-i-1)/k));
    // arrete gauche
    for(int i =0; i<NptPerV;i++) Pt.push_back(R2(double(i+1)/k,0./k)); 
    // Interieur de l'element
    NumInternPoint2( k,Pt );
    return Pt;
}

// Remplir le tableau data
void FillDataLagrange(int k, vector<int> &Data, vector<double> &Pi_h_coef ){
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
// if orientation is <-1 it provides the indices of nodes to be exchanged  (Valable que pour K3 et K4)
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
// Generic de Exchangeidx
vector<pair<int,int>> ExchangeidxVector(int k){
    int nbrPerm=(k%2)? (k-1) : (k-2); // Number of permutation per axis
    vector<pair<int,int>> permut;
    for (int i=0;i<nbrPerm/2;i++){
        permut.push_back(make_pair(3+i,3+(k-2)-i));
    }
    for (int i=0;i<nbrPerm/2;i++){
        permut.push_back(make_pair(3+k-1+i,3+(2*k-3)-i));
    }
    for (int i=0;i<nbrPerm/2;i++){
        permut.push_back(make_pair(3+2*k-2+i,3+(3*k-4)-i));
    }
    return permut;
}
// Fill Other table
void FillOther(vector<int> &other, int PK, int Ndof){
    for(int i=0; i<Ndof;i++ ) other[i]=i;
    for(int i=0; i<3;i++ ) {
        for (int j=0; j<((PK-1)/2);j++ ){
            other[3+i*(PK-1)+j]=3+(i+1)*(PK-1)-j-1;
            other[3+(i+1)*(PK-1)-j-1]=3+i*(PK-1)+j;
        }
    }
}