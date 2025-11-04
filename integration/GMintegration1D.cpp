#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;

#include "RNM.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "QuadratureFormular.hpp"
using namespace Fem2D;

#include<vector>

double factorial(int nbr){
    if (nbr<=1) return 1.;
    else{
        double prod=1.;
        for (int i=1;i<=nbr; i++){
            prod*=i;
        }
        return prod; 
    }
}

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

vector<pair<double, vector<double>>> integration_weightspoints(int n, int s){
    int d = 2*s+1;
    vector<pair<double, vector<double>>> data;
    
    for (int i = 0; i < s+1; i++){
        double coef = pow(-1,i)*pow(2,-2*s)*pow(d+n-2*i,d)/(factorial(i)*factorial(d+n-i));
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
/*const QuadratureFormular *GenerateQuadratureFormular(const int s, int n=1){
  vector<pair<double, vector<double>>> Vec = integration_weightspoints( n,s);
  int NbrP=Vec.size();
   QuadraturePoint *Tab = new QuadraturePoint [NbrP];
  for(int i=0; i<NbrP; i++ ){
    Tab[i]=QuadraturePoint(R1(Vec[i].first*factorial(n)),Vec[i].second[0]);
  }
  const QuadratureFormular * QuadratureFormular_GM_nD_Ss = new QuadratureFormular(s,NbrP, Tab);
  return QuadratureFormular_GM_nD_Ss;
}*/

static QuadratureFormular1d ** TabQuadrFormula = new QuadratureFormular1d* [1000];

const QuadratureFormular1d *GenerateQuadratureFormularForOperator(Stack stack, const long &s){
  if (TabQuadrFormula[s]!=nullptr) {return TabQuadrFormula[s];}
  else{
    int n=1;
    vector<pair<double, vector<double>>> Vec = integration_weightspoints( n,s);
    int NbrP=Vec.size();
     QuadratureFormular1dPoint *Tab = new QuadratureFormular1dPoint [NbrP];
    for(int i=0; i<NbrP; i++ ){
      Tab[i]=QuadratureFormular1dPoint((Vec[i].first*factorial(n)),(Vec[i].second[0]));
    }
    TabQuadrFormula[s] = new QuadratureFormular1d(s,NbrP, Tab);
    return TabQuadrFormula[s];
  }

}


#include "lex.hpp"
extern mylex *zzzfff;
static void Load_Init( ) {

  for (int i=0;i<10; i++) TabQuadrFormula[i]=nullptr;
  Global.Add(
      "GMQuadrature1D", "(",
      new OneOperator1s_<const QuadratureFormular1d*,long>(GenerateQuadratureFormularForOperator)
  );

}

LOADFUNC(Load_Init)
