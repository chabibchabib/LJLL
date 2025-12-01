#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include "AFunction.hpp"
#include "error.hpp"

#include "RNM.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "QuadratureFormular.hpp"
using namespace Fem2D;
using namespace std;

// Compute the factorial
long double factorial(int nbr) {
    if (nbr <= 1)
        return 1.;
    else {
        long double prod = 1.;
        for (int i = 1; i <= nbr; i++) {
            prod *= i;
        }
        return prod;
    }
}

// 
vector<vector<long double>> generate_partitions(int sum, int nparts) {
    vector<vector<long double>> partitions;

    if (nparts == 1) {
        partitions.push_back({(long double)(sum)});
        return partitions;
    }

    for (int first = 0; first <= sum; first++) {
        vector<vector<long double>> rest_partitions = generate_partitions(sum - first, nparts - 1);
        for (auto &rest : rest_partitions) {
            vector<long double> comp;
            comp.push_back(first);
            comp.insert(comp.end(), rest.begin(), rest.end());
            partitions.push_back(comp);
        }
    }
    return partitions;
}

vector<pair<long double, vector<long double>>> integration_weightspoints(int n, int s) {
    int d = 2 * s + 1;
    vector<pair<long double, vector<long double>>> data;
    long double signe =1.;
    for (int i = 0; i < s + 1; i++) {
        long double dn_2i=d+n-2*i;
        long double coef = signe;
        int j_hi = max(n, max(d, d + n - i));
        // Loop to avoid the overflow for a better understanding check the formula 
        for(int j = 1; j <= j_hi; j++)
        {
        long double rj= j;
        if(j <= n) { coef *= rj; }
        if(j <= d) { coef *= dn_2i; }
        if(j <= 2 * s) { coef /= 2.0; }
        if(j <= i) { coef /= rj; }
        if(j <= d + n - i) { coef /= rj; }
        }
        signe = - signe;


        vector<vector<long double>> betas = generate_partitions(s - i, n + 1);
        long double denom = d + n - 2 * i;

        for (const auto &beta_orig : betas) {
            vector<long double> beta;
            for (auto x : beta_orig) {
                beta.push_back((2 * x + 1) / denom);
            }
            data.push_back({coef, beta}); // weight, points for each beta
        }
    }

    return data;
}

/*static QuadratureFormular1d **TabQuadrFormula1 = new QuadratureFormular1d *[1000];

const QuadratureFormular1d *GenerateQuadratureFormularForOperator(Stack stack, const long &s) {
    if (TabQuadrFormula1[s] != nullptr) {
        return TabQuadrFormula1[s];
    } else {
        int n = 1;
        vector<pair<long double, vector<long double>>> Vec = integration_weightspoints(n, s);
        int NbrP = Vec.size();
        QuadratureFormular1dPoint *Tab = new QuadratureFormular1dPoint[NbrP];
        for (int i = 0; i < NbrP; i++) {
            Tab[i] = QuadratureFormular1dPoint((Vec[i].first * factorial(n)), (Vec[i].second[0]));
        }
        TabQuadrFormula1[s] = new QuadratureFormular1d(s, NbrP, Tab);
        return TabQuadrFormula1[s];
    }
}*/

/*static QuadratureFormular **TabQuadrFormula2 = new QuadratureFormular *[1000];

const QuadratureFormular *GenerateQuadratureFormularForOperator2(Stack stack, const long &s) {
    if (TabQuadrFormula2[s] != nullptr) {
        return TabQuadrFormula2[s];
    } else {
        int n = 2;
        vector<pair<long double, vector<long double>>> Vec = integration_weightspoints(n, s);
        int NbrP = Vec.size();
        QuadraturePoint *Tab = new QuadraturePoint[NbrP];
        for (int i = 0; i < NbrP; i++) {
            Tab[i] = QuadraturePoint(Vec[i].first * factorial(n), Vec[i].second[0], Vec[i].second[1]);
        }
        TabQuadrFormula2[s] = new QuadratureFormular(s, NbrP, Tab);
        // Add2StackOfPtr2Free(stack, TabQuadrFormula[s]);
        return TabQuadrFormula2[s];
    }
}*/


typedef GQuadraturePoint<R1> PQP1;
typedef GQuadratureFormular<R1> PQF1;
static PQF1 **TabQuadrFormula1 = new PQF1 *[1000];

const PQF1 *GenerateQuadratureFormularForOperator1(Stack stack, const long &s) {
    if (TabQuadrFormula1[s] != nullptr) {
        return TabQuadrFormula1[s];
    } else {
        int n = 1;
        vector<pair<long double, vector<long double>>> Vec = integration_weightspoints(n, s);
        int NbrP = Vec.size();
        PQP1 *Tab = new PQP1[NbrP];
        for (int i = 0; i < NbrP; i++) {
            Tab[i] = PQP1(R1(Vec[i].second[0]), Vec[i].first * factorial(n ));
        }
        TabQuadrFormula1[s] = new PQF1(s, NbrP, Tab);
        // Add2StackOfPtr2Free(stack, TabQuadrFormula[s]);
        return TabQuadrFormula1[s];
    }
}

typedef GQuadraturePoint<R2> PQP2;
typedef GQuadratureFormular<R2> PQF2;
static PQF2 **TabQuadrFormula2 = new PQF2 *[1000];

const PQF2 *GenerateQuadratureFormularForOperator2(Stack stack, const long &s) {
    if (TabQuadrFormula2[s] != nullptr) {
        return TabQuadrFormula2[s];
    } else {
        int n = 2;
        vector<pair<long double, vector<long double>>> Vec = integration_weightspoints(n, s);
        int NbrP = Vec.size();
        PQP2 *Tab = new PQP2[NbrP];
        for (int i = 0; i < NbrP; i++) {
            Tab[i] = PQP2(R2(Vec[i].second[0], Vec[i].second[1]), Vec[i].first * factorial(n ));
        }
        TabQuadrFormula2[s] = new PQF2(s, NbrP, Tab);
        // Add2StackOfPtr2Free(stack, TabQuadrFormula[s]);
        return TabQuadrFormula2[s];
    }
}



typedef GQuadraturePoint<R3> PQP3;
typedef GQuadratureFormular<R3> PQF3;
static PQF3 **TabQuadrFormula3 = new PQF3 *[1000];

const PQF3 *GenerateQuadratureFormularForOperator3(Stack stack, const long &s) {
    if (TabQuadrFormula3[s] != nullptr) {
        return TabQuadrFormula3[s];
    } else {
        int n = 3;
        vector<pair<long double, vector<long double>>> Vec = integration_weightspoints(n, s);
        int NbrP = Vec.size();
        PQP3 *Tab = new PQP3[NbrP];
        for (int i = 0; i < NbrP; i++) {
            Tab[i] = PQP3(R3(Vec[i].second[0], Vec[i].second[1], Vec[i].second[2]), Vec[i].first * factorial(n ));
        }
        TabQuadrFormula3[s] = new PQF3(s, NbrP, Tab);
        // Add2StackOfPtr2Free(stack, TabQuadrFormula[s]);
        return TabQuadrFormula3[s];
    }
}

#include "lex.hpp"
extern mylex *zzzfff;
static void Load_Init() {

    for (int i = 0; i < 1000; i++){
        TabQuadrFormula1[i] = nullptr;
        TabQuadrFormula2[i] = nullptr;
        TabQuadrFormula3[i] = nullptr;

    }
    //Global.Add("GMQuadrature1D", "(", new OneOperator1s_<const QuadratureFormular1d *, long>(GenerateQuadratureFormularForOperator));
    //Global.Add("GMQuadrature2D", "(", new OneOperator1s_<const QuadratureFormular *, long>(GenerateQuadratureFormularForOperator2));
    Global.Add("GMQuadrature1D", "(", new OneOperator1s_<const PQF1 *, long>(GenerateQuadratureFormularForOperator1));
    Global.Add("GMQuadrature2D", "(", new OneOperator1s_<const PQF2 *, long>(GenerateQuadratureFormularForOperator2));
    Global.Add("GMQuadrature3D", "(", new OneOperator1s_<const PQF3 *, long>(GenerateQuadratureFormularForOperator3));

}

LOADFUNC(Load_Init)