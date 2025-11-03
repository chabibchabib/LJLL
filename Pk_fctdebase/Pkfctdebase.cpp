#include<iostream>
#include<vector>
#include <tuple>
#include <cmath>
#include <functional>
#include <algorithm>
#include <set>
#include<map>
#include <fstream>
#include <iomanip>      // std::setprecision

using namespace std;
typedef double REAL  ; 
vector<vector<int>> generate_partitions(int sum, int nparts) {
    vector<vector<int>> partitions;
    
    if (nparts == 1) {
        partitions.push_back({int(sum)});
        return partitions;
    }
    
    for (int first = 0; first <= sum; first++) {
        vector<vector<int>> rest_partitions = generate_partitions(sum - first, nparts - 1);
        for (auto &rest : rest_partitions) {
            vector<int> comp;
            comp.push_back(first);
            comp.insert(comp.end(), rest.begin(), rest.end());
            partitions.push_back(comp);
        }
    }
    return partitions;
}

// Fonction pour trier les partitions selon l'ordre geometrique P_K
vector<vector<int>> sort_partitions_geometric(vector<vector<int>> partitions, int K) {
    vector<vector<int>> sorted;
    
    // 1. SOMMETS (i+j+k=K avec deux coordonnees = 0)
    sorted.push_back({K,0,0});
    sorted.push_back({0,K,0});
    sorted.push_back({0,0,K});    
    // 2. ARETES (i+j+k=K avec une coordonnee = 0, les autres > 0)
    // Ordre: arete k=0 (sommet 0→1), puis j=0 (sommet 1→2), puis i=0 (sommet 2→0)
    
    // Arete k=0 (entre sommets i=K et j=K) : decroissant en i
    vector<vector<int>> edge1, edge2, edge3;
    for (auto &p : partitions) {
        if (p[2] == 0 && p[0] > 0 && p[1] > 0) edge1.push_back(p);
        if (p[1] == 0 && p[0] > 0 && p[2] > 0) edge2.push_back(p);
        if (p[0] == 0 && p[1] > 0 && p[2] > 0) edge3.push_back(p);
    }
    
    // Trier arete 1 (k=0): i decroissant
    sort(edge1.begin(), edge1.end(), [](const vector<int> &a, const vector<int> &b) {
        return a[0] > b[0];
    });
    
    // Trier arete 2 (j=0): i croissant  
    sort(edge2.begin(), edge2.end(), [](const vector<int> &a, const vector<int> &b) {
        return a[0] < b[0];
    });
    
    // Trier arete 3 (i=0): j decroissant
    sort(edge3.begin(), edge3.end(), [](const vector<int> &a, const vector<int> &b) {
        return a[1] > b[1];
    });
    
    // Ajouter les aretes
    sorted.insert(sorted.end(), edge3.begin(), edge3.end());
    sorted.insert(sorted.end(), edge2.begin(), edge2.end());
    sorted.insert(sorted.end(), edge1.begin(), edge1.end());
    
    // 3. INTERIEUR (i,j,k  > 0)
    vector<vector<int>> interior;
    for (auto &p : partitions) {
        if (p[0] > 0 && p[1] > 0 && p[2] > 0) {
            interior.push_back(p);
        }
    }
    
    // Trier interieur par k decroissant, puis j croissant
    sort(interior.begin(), interior.end(), [](const vector<int> &a, const vector<int> &b) {
        if (a[2] != b[2]) return a[2] > b[2];
        return a[1] < b[1];
    });
    
    sorted.insert(sorted.end(), interior.begin(), interior.end());
    
    return sorted;
}

REAL EvaluateBasisFct( int i, int j , int k, REAL lambda1, REAL lambda2, REAL lambda3, int K ){ 
    // Silverster formula ~1968 
    // Compute Phi_{i,j,k}(lambda1,lambda2,lambda3) ; while i+j+k=K and lambda are the basis coordinates  
    double result=1;
    double denom=tgamma(i+1)*tgamma(j+1)*tgamma(k+1);
    if( (i+j+k==K) && (lambda1 +lambda2 +lambda3 ==1)){
            if(i>0){
                for (int ii = 0; ii<=i-1;ii++) result *= (K*lambda1 - ii );
            }
            if(j>0){
                for (int jj = 0; jj<=j-1;jj++) result *= (K*lambda2 - jj );
            }
            if(k>0){
                for (int kk = 0; kk<=k-1;kk++) result *= (K*lambda3 - kk );
            }
        return result/denom;
        }
    else return -1;
}

void BasisFctPK(int K, vector<vector<int>> &coef,vector<vector<int>> &shift, vector<int> &ff, vector<int> &il,vector<int> &jl,vector<int> &kl){ 
    // Silverster formula ~1968 
    // Compute Phi_{i,j,k}(lambda1,lambda2,lambda3) basis functions; while i+j+k=K and lambda are the barycentric coordinates  
    REAL result=1;

    vector<vector<int>> unsortedpartitionK=generate_partitions(K, 3);
    vector<vector<int>> partitionK=sort_partitions_geometric(unsortedpartitionK,  K);
    int idx=0;
    for (auto &partition : partitionK ){
        /*vector<vector<REAL>> coef;//(3,vector<REAL>(K));
        vector<vector<REAL>> shift;//(3,vector<REAL>(K));
        coef.resize(3);
        shift.resize(3);

        for (int i=0; i<3;i++){
            coef[i].resize(K,-1);
            shift[i].resize(K,-1);

        }*/
        int i = partition[0];
        int j = partition[1];
        int k = partition[2];
        if (i+j+k==K){
            int ID=0;
            cout<<"Point -->"<< idx<<"( "<<i<<" "<<j<<" "<<k<<" )" << endl;
            REAL denom=tgamma(i+1)*tgamma(j+1)*tgamma(k+1);
            ff[idx]=denom;
            il[idx]=i;
            jl[idx]=j;
            kl[idx]=k;
            if(i>0){
                for (int ii = 0; ii<=i-1;ii++) {//result *= (K*lambda1 - ii );
                        /*coef[0][ii]=0;
                        shift[0][ii]=(REAL(ii));*/
                        coef[idx][ID]=0;
                        shift[idx][ID]=ii;
                        ID++;


                }
                //if(coef[0][0] !=0) {coef[0][0]/=denom; shift[0][0]/=denom;} 
            }
            if(j>0){
                for (int jj = 0; jj<=j-1;jj++) {//result *= (K*lambda2 - jj );
                        /*coef[1][jj]=1;
                        shift[1][jj]=(REAL(jj));*/
                        coef[idx][ID]=1;
                        shift[idx][ID]=jj;
                        ID++;
                }
                //if(coef[1][0] !=0 && coef[0][0] ==0 ) {coef[1][0]/=denom; shift[1][0]/=denom;}
            }
            if(k>0){
                for (int kk = 0; kk<=k-1;kk++){ //result *= (K*lambda3 - kk );
                        /*coef[2][kk]=2;
                        shift[2][kk]=(REAL(kk));*/
                        coef[idx][ID]=2;
                        shift[idx][ID]=kk;
                        ID++;
                }
                //if(coef[2][0] !=0 && coef[1][0] ==0 && coef[0][0] ==0 ) {coef[2][0]/=denom; shift[2][0]/=denom;}

            }
            //cout<<"ID="<<ID<<endl;
            idx++;
            /*for (int j =0; j<coef[0].size();j++){
                //if(coef[0][j] !=0 ||shift[0][j] !=0 )
                if(coef[0][j] !=-1)
                    cout<<"\t denom="<<denom<< "\t C[0]= "<<coef[0][j]<<" "<<"S[0]= "<<shift[0][j]<<endl;
                //if(coef[1][j] !=0 ||shift[1][j] !=0 )
                if(coef[1][j] !=-1)
                    cout<<"\t denom="<<denom<< "\t C[1]= "<< coef[1][j]<<" "<<"S[1]= "<<shift[1][j]<<endl;
                //if(coef[2][j] !=0 ||shift[2][j] !=0 )
                if(coef[2][j] !=-1)
                    cout<<"\t denom="<<denom<< "\t C[2]= "<< coef[2][j]<<" "<<"S[2]= "<<shift[2][j]<<endl;

            }*/
        }
    }
}

double fct_de_base(int idx, double x , double y, int k){
    if(idx==0) return k*(1.-x-y);
    else if(idx==1) return double(k*x);
    else return double(k*y);
}
int main(int argc, char **argv){
    int PK= 8;
    int i =1, j=1,k=0;
    int ndof=(PK+1)*(PK+2)*0.5;
    vector<vector<int>> coef;//(3,vector<REAL>(K));
    vector<vector<int>> shift;//(3,vector<REAL>(K));
    coef.resize(ndof);
    shift.resize(ndof);
    for (int i=0; i<ndof;i++){
        coef[i].resize(PK,-1);
        shift[i].resize(PK,-1);
    }
    vector<int> ff,il,jl,kl;
    ff.resize(ndof,-1);
    il.resize(ndof,-1);
    jl.resize(ndof,-1);
    kl.resize(ndof,-1);
    BasisFctPK( PK, coef,shift, ff, il,jl,kl);
    cout<<"coordinates:"<<endl;
    for (int i =0;i<ndof;i++){
        cout<<"{"<<il[i]<<","<<jl[i]<<","<<kl[i]<<"}\t";
    }
    cout<<"\n denom:"<<endl;
    for (int i =0;i<ndof;i++){
        cout<<"{"<<ff[i]<<","<<ff[i]<<","<<ff[i]<<"}\t";
    }
    cout<<"\nCoef:"<<endl;
    for (int i =0;i<ndof;i++){
        cout<<"{";
        for (int j=0; j<PK;j++) {cout<<coef[i][j]; if(j!=PK-1) cout<<",";}
        cout<<"}\t";    }
    cout<<"\n Shift:"<<endl;
    for (int i =0;i<ndof;i++){
        cout<<"{";
        for (int j=0; j<PK;j++) {cout<<shift[i][j]; if(j!=PK-1) cout<<",";}
        cout<<"}\t";

    }
    cout<<endl;
    /*// Verification
    for (int i=0; i<ndof;i++){
        long double prod=1./ff[i];
        for (int j=0; j<PK;j++) {prod*=fct_de_base(coef[i][j],jl[i],  kl[i],PK)-shift[i][j]; }
        cout<<"\nNdof= "<< i<<" prod= "<< prod<<endl;
    }*/
// Verification
vector< vector<double> > verif(ndof);
for (int i = 0; i < ndof; i++) verif[i].resize(ndof,-1);
for (int i = 0; i < ndof; i++) {
    for (int idx=0;idx<ndof;idx++){
        // Coordonnées barycentriques du point i
            double lambda1 = ( double)il[idx] / PK;
            double lambda2 = ( double)jl[idx] / PK;
            double lambda3 = ( double)kl[idx] / PK;
            
            double prod = 1.0 / ff[i];
            
            for (int j = 0; j < PK; j++) {
                if (coef[i][j] == -1) break; // Fin des termes valides
                
                // Choisir la bonne coordonnée barycentrique
                double lambda_val;
                if (coef[i][j] == 0) 
                    lambda_val = lambda1;
                else if (coef[i][j] == 1) 
                    lambda_val = lambda2;
                else // coef[i][j] == 2
                    lambda_val = lambda3;
                
                prod *= (PK * lambda_val - shift[i][j]);
            }
            verif[i][idx]=prod;
            //cout << "Ndof= " << i << " prod= " << prod << endl;
        }
    }
    double norm=0.;
    for (int i=0;i<ndof;i++){
        for (int j=0;j<ndof;j++){ //norm+=
            if( i==j  && verif[i][j]!=1) cout  << verif[i][j] << endl;

        }
    }
    
}