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


void BasisFctPK(int K, vector<vector<int>> &coef,vector<vector<int>> &shift, vector<int> &ff, vector<int> &il,vector<int> &jl,vector<int> &kl){ 
    // Silverster formula ~1968 
    // Compute Phi_{i,j,k}(lambda1,lambda2,lambda3) basis functions; while i+j+k=K and lambda are the barycentric coordinates  

    vector<vector<int>> unsortedpartitionK=generate_partitions(K, 3);
    vector<vector<int>> partitionK=sort_partitions_geometric(unsortedpartitionK,  K);
    int idx=0;
    for (auto &partition : partitionK ){
        int i = partition[0];
        int j = partition[1];
        int k = partition[2];
        if (i+j+k==K){
            int ID=0;
            cout<<"Point -->"<< idx<<"( "<<i<<" "<<j<<" "<<k<<" )" << endl;
            int denom=tgamma(i+1)*tgamma(j+1)*tgamma(k+1);
            ff[idx]=denom;
            il[idx]=i;
            jl[idx]=j;
            kl[idx]=k;
            if(i>0){
                for (int ii = 0; ii<=i-1;ii++) {

                        coef[idx][ID]=0;
                        shift[idx][ID]=ii;
                        ID++;
                }
            }
            if(j>0){
                for (int jj = 0; jj<=j-1;jj++) {
                        coef[idx][ID]=1;
                        shift[idx][ID]=jj;
                        ID++;
                }
            }
            if(k>0){
                for (int kk = 0; kk<=k-1;kk++){ 
                        coef[idx][ID]=2;
                        shift[idx][ID]=kk;
                        ID++;
                }
            }
            idx++;
        }
    }
}
