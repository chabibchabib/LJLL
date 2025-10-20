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

vector<vector<double>> generate_partitions(int sum, int nparts) {
    vector<vector<double>> partitions;
    
    if (nparts == 1) {
        partitions.push_back({double(sum)});
        return partitions;
    }
    
    for (int first = 0; first <= sum; first++) {
        vector<vector<double>> rest_partitions = generate_partitions(sum - first, nparts - 1);
        for (auto &rest : rest_partitions) {
            vector<double> comp;
            comp.push_back(first);
            comp.insert(comp.end(), rest.begin(), rest.end());
            partitions.push_back(comp);
        }
    }
    return partitions;
}

double integration_rule(int n, int s,function<double(const vector<double>&)> f){
    int d = 2*s+1;
    double result =0.;
    for (int i = 0 ; i<s+1; i++){
        double coef = pow(-1,i)*pow(2,-2*s)*pow(d+n-2*i,d)/(tgamma(i+1)*tgamma(d+n-i+1)); //C=(-1)^i x 2^(-2s) x (d+n-2i)^d/(i! x (d+n-i) !)

        vector<vector<double>> betas =  generate_partitions( s-i, n+1);

        double denom = d+n-2*i; 
        for (const auto & beta_orig : betas ){
            vector<double> beta;
            for ( auto x : beta_orig){
                //x = (2*x +1)/denom;
                beta.push_back((2*x + 1)/denom);
            }
        
        /*cout << "i=" << i << ", beta=(";
        for(auto xi : beta) cout << xi << " ";
        cout << "), coef=" << coef << endl;*/

    result += coef * f(beta);
        }
    }
    return result;
}



double fct (vector<double> x ){
    /*double y=.0;
    int idx=0;
    for(auto & xx : x) { y+=pow(xx,2); }//y[idx]=pow(xx,2); idx++;}*/
    double y = 0.0;
    /*for (auto it = x.begin(); it != x.end() -1; it++) {
        y += pow(*it, 2);
    }*/
    y =pow(x[0],5)+pow(x[1],2)+2*x[0]*x[1]+pow(x[2],2);
   return y;
}; 

double evaluate_integrale(const vector<pair<double, vector<double>>>& data, double (*fct)(vector<double> x)){
    double result = 0.; 
    for (const auto &pair : data){
        const auto & weight = pair.first;
        const auto & point = pair.second;
        double y = fct(point);
        result += weight * y;
    }
    return result;
}

vector<pair<double, vector<double>>> integration_weightspoints(int n, int s){
    int d = 2*s+1;
    vector<pair<double, vector<double>>> data;
    
    for (int i = 0; i < s+1; i++){
        double coef = pow(-1,i)*pow(2,-2*s)*pow(d+n-2*i,d)/(tgamma(i+1)*tgamma(d+n-i+1));
        vector<vector<double>> betas = generate_partitions(s-i, n+1);
        double denom = d+n-2*i;
        
        for (const auto & beta_orig : betas){
            vector<double> beta;
            for (auto x : beta_orig){
                beta.push_back((2*x + 1)/denom);
            }
            data.push_back({coef, beta});  //weight, points for each beta
        }
    }
    
    return data;
}

int nbrOfIntgPOints(int n,int s) {
    int nbr=0;
    for (int i=0; i<=s;i++)
        nbr+=tgamma(n+s-i+1)/(tgamma(s-i+1)*tgamma(n+1));
    return nbr;
}
int main(int argc,char * argv[] ){
    int n=atoi(argv[1]),s=atoi(argv[2]);

    auto f= [] (vector<double> x ){
    return pow(x[0],5)+pow(x[1],2)+2*x[0]*x[1]+pow(x[2],2);
    //return 1.;
    //return x[0];
    };

    double result = integration_rule( n, s,f);
    cout<<"result:"<<result<<endl;

    auto  data= integration_weightspoints( n,  s);
    double result2=evaluate_integrale( data,fct );
    cout<<"result2:"<<result2<<endl;

    /*for (auto & element : data){
        //cout<<element.first<<"\t"<<element.second[0]<<endl;
        cout<<endl;
        cout<<element.first<<"\t";
        for(auto &val : element.second) 
            cout<<val<<"\t";

        }*/
    cout<<endl;
    ofstream file("output.txt");
    if (!file) {
        cerr << "Error file" << endl;
        return 1;
    }

    file << fixed << setprecision(13);  
    for (auto & element : data){
        //file<<endl;
        file<<"PQP3(R3(";
        for (auto it = element.second.begin(); it != element.second.end()-1 ; it++) 
        
        {   if(it != element.second.end()-2)
                file<<*it<<",";
            else 
                file<<*it<<"),";
        }
        //file<<"),";
        file<<element.first<<"*6),"<<endl;
        
    }
    file.close();
    cout<<"N= "<<n<<" s= "<<s<<" Nbr Points= "<<nbrOfIntgPOints( n, s)<<endl;
    return 0;
}
