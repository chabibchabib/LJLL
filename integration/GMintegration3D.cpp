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
typedef GQuadraturePoint< R3 > PQP3;
typedef GQuadratureFormular< R3 > PQF3;
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

static PQF3 ** TabQuadrFormula = new PQF3* [1000];

const PQF3 *GenerateQuadratureFormularForOperator(Stack stack, const long &s){
  if (TabQuadrFormula[s]!=nullptr) {return TabQuadrFormula[s];}
  else{
    int n=3;
    vector<pair<double, vector<double>>> Vec = integration_weightspoints( n,s);
    int NbrP=Vec.size();
     PQP3 *Tab = new PQP3 [NbrP];
    for(int i=0; i<NbrP; i++ ){
      Tab[i]=PQP3(R3(Vec[i].second[0],Vec[i].second[1],Vec[i].second[2]),Vec[i].first*tgamma(n+1));
    }
    TabQuadrFormula[s] = new PQF3(s,NbrP, Tab);
    //Add2StackOfPtr2Free(stack, TabQuadrFormula[s]);
    return TabQuadrFormula[s];
  }

}

#include "lex.hpp"
extern mylex *zzzfff;
static void Load_Init( ) {
  

  for (int i=0;i<1000; i++) TabQuadrFormula[i]=nullptr;
  Global.Add(
      "GMQuadrature2D", "(",
      new OneOperator1s_<const PQF3*,long>(GenerateQuadratureFormularForOperator)
  );

}

LOADFUNC(Load_Init)
//delete [] TabQuadrFormula;
