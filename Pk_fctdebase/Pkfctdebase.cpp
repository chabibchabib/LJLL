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

void BasisFctPK(int K ){ 
    // Silverster formula ~1968 
    // Compute Phi_{i,j,k}(lambda1,lambda2,lambda3) basis functions; while i+j+k=K and lambda are the barycentric coordinates  
    REAL result=1;

    vector<vector<int>> partitionK=generate_partitions(K, 3);
    int idx=0;
    for (auto &partition : partitionK ){
        vector<vector<REAL>> coef(3,vector<REAL>(K));
        vector<vector<REAL>> shift(3,vector<REAL>(K));
        int i = partition[0];
        int j = partition[1];
        int k = partition[2];
        if (i+j+k==K){
            cout<<"Point -->"<< idx<<"( "<<i<<" "<<j<<" "<<k<<" )" << endl;
            idx++;
            REAL denom=tgamma(i+1)*tgamma(j+1)*tgamma(k+1);
            if(i>0){
                for (int ii = 0; ii<=i-1;ii++) {//result *= (K*lambda1 - ii );
                        coef[0][ii]=1.;
                        shift[0][ii]=(REAL(ii));

                }
                //if(coef[0][0] !=0) {coef[0][0]/=denom; shift[0][0]/=denom;} 
            }
            if(j>0){
                for (int jj = 0; jj<=j-1;jj++) {//result *= (K*lambda2 - jj );
                        coef[1][jj]=1.;
                        shift[1][jj]=(REAL(jj));
                }
                //if(coef[1][0] !=0 && coef[0][0] ==0 ) {coef[1][0]/=denom; shift[1][0]/=denom;}
            }
            if(k>0){
                for (int kk = 0; kk<=k-1;kk++){ //result *= (K*lambda3 - kk );
                        coef[2][kk]=1.;
                        shift[2][kk]=(REAL(kk));
                }
                //if(coef[2][0] !=0 && coef[1][0] ==0 && coef[0][0] ==0 ) {coef[2][0]/=denom; shift[2][0]/=denom;}

            }

            for (int j =0; j<coef[0].size();j++){
                if(coef[0][j] !=0 ||shift[0][j] !=0 )
                    cout<<"\t denom="<<denom<< "\t C[0]= "<<coef[0][j]<<" "<<"S[0]= "<<shift[0][j]<<endl;
                if(coef[1][j] !=0 ||shift[1][j] !=0 )
                    cout<<"\t denom="<<denom<< "\t C[1]= "<< coef[1][j]<<" "<<"S[1]= "<<shift[1][j]<<endl;
                if(coef[2][j] !=0 ||shift[2][j] !=0 )                
                    cout<<"\t denom="<<denom<< "\t C[2]= "<< coef[2][j]<<" "<<"S[2]= "<<shift[2][j]<<endl;

            }
        }
    }
}

int main(int argc, char **argv){
    int PK= 3;
    int i =1, j=1,k=0;
    //REAL result= EvaluateBasisFct(  i,  j ,  k,REAL(i)/PK ,REAL(j)/PK,REAL(k)/PK,PK );   
    //cout<<"Result: "<<result<<endl;
    BasisFctPK(PK );

}