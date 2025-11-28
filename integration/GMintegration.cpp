#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <vector>

#include "error.hpp"
#include "AFunction.hpp"

#include "RNM.hpp"
#include "rgraph.hpp"
#include "fem.hpp"

#include "QuadratureFormular.hpp"

using namespace Fem2D;
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

vector<pair<double, vector<double>>> integration_weightspoints(int n, int s) {
    int d = 2 * s + 1;
    vector<pair<double, vector<double>>> data;

    for (int i = 0; i < s + 1; i++) {
        double coef = pow(-1, i) * pow(2, -2 * s) * pow(d + n - 2 * i, d) / (tgamma(i + 1) * tgamma(d + n - i + 1));
        vector<vector<double>> betas = generate_partitions(s - i, n + 1);
        double denom = d + n - 2 * i;

        for (const auto &beta_orig : betas) {
            vector<double> beta;
            for (auto x : beta_orig) {
                beta.push_back((2 * x + 1) / denom);
            }
            data.push_back({coef, beta}); // weight, points for each beta
        }
    }

    return data;
}
const QuadratureFormular *GenerateQuadratureFormular(const int s, int n = 2) {
    vector<pair<double, vector<double>>> Vec = integration_weightspoints(n, s);
    int NbrP = Vec.size();
    QuadraturePoint *Tab = new QuadraturePoint[NbrP];
    for (int i = 0; i < NbrP; i++) {
        Tab[i] = QuadraturePoint(Vec[i].first * tgamma(n + 1), Vec[i].second[0], Vec[i].second[1]);
    }
    const QuadratureFormular *QuadratureFormular_GM_nD_Ss = new QuadratureFormular(s, NbrP, Tab);
    // Add2StackOfPtr2Free(stack, QuadratureFormular_GM_nD_Ss);
    return QuadratureFormular_GM_nD_Ss;
}

static QuadratureFormular **TabQuadrFormula = new QuadratureFormular *[1000];

const QuadratureFormular *GenerateQuadratureFormularForOperator(Stack stack, const long &s) {
    if (TabQuadrFormula[s] != nullptr) {
        return TabQuadrFormula[s];
    } else {
        int n = 2;
        vector<pair<double, vector<double>>> Vec = integration_weightspoints(n, s);
        int NbrP = Vec.size();
        QuadraturePoint *Tab = new QuadraturePoint[NbrP];
        for (int i = 0; i < NbrP; i++) {
            Tab[i] = QuadraturePoint(Vec[i].first * tgamma(n + 1), Vec[i].second[0], Vec[i].second[1]);
        }
        TabQuadrFormula[s] = new QuadratureFormular(s, NbrP, Tab);
        // Add2StackOfPtr2Free(stack, TabQuadrFormula[s]);
        return TabQuadrFormula[s];
    }
}
// 1D
/*static QuadraturePoint P_GM_1D_S0[] = {
QuadraturePoint(1.0000000000000	,0.5000000000000),
};

static QuadraturePoint P_GM_1D_S1[] = {
QuadraturePoint(0.6666666666667	,0.2500000000000),
QuadraturePoint(0.6666666666667	,0.7500000000000),
QuadraturePoint(-0.3333333333333,	0.5000000000000),
};

static QuadraturePoint P_GM_1D_S2[] = {
QuadraturePoint(0.6750000000000	,0.1666666666667),
QuadraturePoint(0.6750000000000	,0.5000000000000),
QuadraturePoint(0.6750000000000	,0.8333333333333),
QuadraturePoint(-0.5333333333333,	0.2500000000000),
QuadraturePoint(-0.5333333333333,	0.7500000000000),
QuadraturePoint(0.0416666666667,	0.5000000000000),
};

static QuadraturePoint P_GM_1D_S3[] = {
QuadraturePoint(0.8126984126984,0.1250000000000),
QuadraturePoint(0.8126984126984,0.3750000000000),
QuadraturePoint(0.8126984126984,0.6250000000000),
QuadraturePoint(0.8126984126984,0.8750000000000),
QuadraturePoint(-0.8678571428571,0.1666666666667),
QuadraturePoint(-0.8678571428571,0.5000000000000),
QuadraturePoint(-0.8678571428571,0.8333333333333),
QuadraturePoint(0.1777777777778,0.2500000000000),
QuadraturePoint(0.1777777777778,0.7500000000000),
QuadraturePoint(-0.0027777777778,0.5000000000000),
};*/
// 2D
static QuadraturePoint P_GM_2D_S0[] = {
    QuadraturePoint(0.5000000000000 * 2, 0.3333333333333, 0.3333333333333),
};
static QuadraturePoint P_GM_2D_S1[] = {
    QuadraturePoint(0.2604166666667 * 2, 0.2000000000000, 0.2000000000000),
    QuadraturePoint(0.2604166666667 * 2, 0.2000000000000, 0.6000000000000),
    QuadraturePoint(0.2604166666667 * 2, 0.6000000000000, 0.2000000000000),
    QuadraturePoint(-0.2812500000000 * 2, 0.3333333333333, 0.3333333333333),
};
static QuadraturePoint P_GM_2D_S2[] = {
    QuadraturePoint(0.2084201388889 * 2, 0.1428571428571, 0.1428571428571),
    QuadraturePoint(0.2084201388889 * 2, 0.1428571428571, 0.4285714285714),
    QuadraturePoint(0.2084201388889 * 2, 0.1428571428571, 0.7142857142857),
    QuadraturePoint(0.2084201388889 * 2, 0.4285714285714, 0.1428571428571),
    QuadraturePoint(0.2084201388889 * 2, 0.4285714285714, 0.4285714285714),
    QuadraturePoint(0.2084201388889 * 2, 0.7142857142857, 0.1428571428571),
    QuadraturePoint(-0.2712673611111 * 2, 0.2000000000000, 0.2000000000000),
    QuadraturePoint(-0.2712673611111 * 2, 0.2000000000000, 0.6000000000000),
    QuadraturePoint(-0.2712673611111 * 2, 0.6000000000000, 0.2000000000000),
    QuadraturePoint(0.0632812500000 * 2, 0.3333333333333, 0.3333333333333),
};
static QuadraturePoint P_GM_2D_S3[] = {
    QuadraturePoint(0.2059465680804 * 2, 0.1111111111111, 0.7777777777778),
    QuadraturePoint(0.2059465680804 * 2, 0.3333333333333, 0.5555555555556),
    QuadraturePoint(0.2059465680804 * 2, 0.5555555555556, 0.3333333333333),
    QuadraturePoint(0.2059465680804 * 2, 0.7777777777778, 0.1111111111111),
    QuadraturePoint(0.2059465680804 * 2, 0.1111111111111, 0.5555555555556),
    QuadraturePoint(0.2059465680804 * 2, 0.3333333333333, 0.3333333333333),
    QuadraturePoint(0.2059465680804 * 2, 0.5555555555556, 0.1111111111111),
    QuadraturePoint(0.2059465680804 * 2, 0.1111111111111, 0.3333333333333),
    QuadraturePoint(0.2059465680804 * 2, 0.3333333333333, 0.1111111111111),
    QuadraturePoint(0.2059465680804 * 2, 0.1111111111111, 0.1111111111111),
    QuadraturePoint(-0.3191433376736 * 2, 0.1428571428571, 0.7142857142857),
    QuadraturePoint(-0.3191433376736 * 2, 0.4285714285714, 0.4285714285714),
    QuadraturePoint(-0.3191433376736 * 2, 0.7142857142857, 0.1428571428571),
    QuadraturePoint(-0.3191433376736 * 2, 0.1428571428571, 0.4285714285714),
    QuadraturePoint(-0.3191433376736 * 2, 0.4285714285714, 0.1428571428571),
    QuadraturePoint(-0.3191433376736 * 2, 0.1428571428571, 0.1428571428571),
    QuadraturePoint(0.1211015004960 * 2, 0.2000000000000, 0.6000000000000),
    QuadraturePoint(0.1211015004960 * 2, 0.6000000000000, 0.2000000000000),
    QuadraturePoint(0.1211015004960 * 2, 0.2000000000000, 0.2000000000000),
    QuadraturePoint(-0.0079101562500 * 2, 0.3333333333333, 0.3333333333333),
};

const QuadratureFormular QuadratureFormular_GM_2D_S0(0, 1, P_GM_2D_S0);
const QuadratureFormular QuadratureFormular_GM_2D_S1(1, 4, P_GM_2D_S1);
const QuadratureFormular QuadratureFormular_GM_2D_S2(2, 10, P_GM_2D_S2);
const QuadratureFormular QuadratureFormular_GM_2D_S3(3, 20, P_GM_2D_S3);
// 3D
typedef GQuadraturePoint<R3> PQP3;
typedef GQuadratureFormular<R3> PQF3;

PQP3 P_GM_3D_S0[] = {
    PQP3(R3(0.2500000000000, 0.2500000000000, 0.2500000000000), 0.1666666666667 * 6),
};
PQF3 const QuadratureFormular_GM_3D_S0(0, 1, P_GM_3D_S0);

PQP3 P_GM_3D_S1[] = {
    PQP3(R3(0.1666666666667, 0.1666666666667, 0.1666666666667), 0.0750000000000 * 6),
    PQP3(R3(0.1666666666667, 0.1666666666667, 0.5000000000000), 0.0750000000000 * 6),
    PQP3(R3(0.1666666666667, 0.5000000000000, 0.1666666666667), 0.0750000000000 * 6),
    PQP3(R3(0.5000000000000, 0.1666666666667, 0.1666666666667), 0.0750000000000 * 6),
    PQP3(R3(0.2500000000000, 0.2500000000000, 0.2500000000000), -0.1333333333333 * 6),

};
PQF3 const QuadratureFormular_GM_3D_S1(1, 5, P_GM_3D_S1);

PQP3 P_GM_3D_S2[] = {
    PQP3(R3(0.1250000000000, 0.1250000000000, 0.1250000000000), 0.0507936507937 * 6),
    PQP3(R3(0.1250000000000, 0.1250000000000, 0.3750000000000), 0.0507936507937 * 6),
    PQP3(R3(0.1250000000000, 0.1250000000000, 0.6250000000000), 0.0507936507937 * 6),
    PQP3(R3(0.1250000000000, 0.3750000000000, 0.1250000000000), 0.0507936507937 * 6),
    PQP3(R3(0.1250000000000, 0.3750000000000, 0.3750000000000), 0.0507936507937 * 6),
    PQP3(R3(0.1250000000000, 0.6250000000000, 0.1250000000000), 0.0507936507937 * 6),
    PQP3(R3(0.3750000000000, 0.1250000000000, 0.1250000000000), 0.0507936507937 * 6),
    PQP3(R3(0.3750000000000, 0.1250000000000, 0.3750000000000), 0.0507936507937 * 6),
    PQP3(R3(0.3750000000000, 0.3750000000000, 0.1250000000000), 0.0507936507937 * 6),
    PQP3(R3(0.6250000000000, 0.1250000000000, 0.1250000000000), 0.0507936507937 * 6),
    PQP3(R3(0.1666666666667, 0.1666666666667, 0.1666666666667), -0.0964285714286 * 6),
    PQP3(R3(0.1666666666667, 0.1666666666667, 0.5000000000000), -0.0964285714286 * 6),
    PQP3(R3(0.1666666666667, 0.5000000000000, 0.1666666666667), -0.0964285714286 * 6),
    PQP3(R3(0.5000000000000, 0.1666666666667, 0.1666666666667), -0.0964285714286 * 6),
    PQP3(R3(0.2500000000000, 0.2500000000000, 0.2500000000000), 0.0444444444444 * 6),
};
PQF3 const QuadratureFormular_GM_3D_S2(2, 15, P_GM_3D_S2);

PQP3 P_GM_3D_S3[] = {
    PQP3(R3(0.1000000000000, 0.1000000000000, 0.1000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.1000000000000, 0.1000000000000, 0.3000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.1000000000000, 0.1000000000000, 0.5000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.1000000000000, 0.1000000000000, 0.7000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.1000000000000, 0.3000000000000, 0.1000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.1000000000000, 0.3000000000000, 0.3000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.1000000000000, 0.3000000000000, 0.5000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.1000000000000, 0.5000000000000, 0.1000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.1000000000000, 0.5000000000000, 0.3000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.1000000000000, 0.7000000000000, 0.1000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.3000000000000, 0.1000000000000, 0.1000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.3000000000000, 0.1000000000000, 0.3000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.3000000000000, 0.1000000000000, 0.5000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.3000000000000, 0.3000000000000, 0.1000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.3000000000000, 0.3000000000000, 0.3000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.3000000000000, 0.5000000000000, 0.1000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.5000000000000, 0.1000000000000, 0.1000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.5000000000000, 0.1000000000000, 0.3000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.5000000000000, 0.3000000000000, 0.1000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.7000000000000, 0.1000000000000, 0.1000000000000), 0.0430583112875 * 6),
    PQP3(R3(0.1250000000000, 0.1250000000000, 0.1250000000000), -0.0902998236332 * 6),
    PQP3(R3(0.1250000000000, 0.1250000000000, 0.3750000000000), -0.0902998236332 * 6),
    PQP3(R3(0.1250000000000, 0.1250000000000, 0.6250000000000), -0.0902998236332 * 6),
    PQP3(R3(0.1250000000000, 0.3750000000000, 0.1250000000000), -0.0902998236332 * 6),
    PQP3(R3(0.1250000000000, 0.3750000000000, 0.3750000000000), -0.0902998236332 * 6),
    PQP3(R3(0.1250000000000, 0.6250000000000, 0.1250000000000), -0.0902998236332 * 6),
    PQP3(R3(0.3750000000000, 0.1250000000000, 0.1250000000000), -0.0902998236332 * 6),
    PQP3(R3(0.3750000000000, 0.1250000000000, 0.3750000000000), -0.0902998236332 * 6),
    PQP3(R3(0.3750000000000, 0.3750000000000, 0.1250000000000), -0.0902998236332 * 6),
    PQP3(R3(0.6250000000000, 0.1250000000000, 0.1250000000000), -0.0902998236332 * 6),
    PQP3(R3(0.1666666666667, 0.1666666666667, 0.1666666666667), 0.0542410714286 * 6),
    PQP3(R3(0.1666666666667, 0.1666666666667, 0.5000000000000), 0.0542410714286 * 6),
    PQP3(R3(0.1666666666667, 0.5000000000000, 0.1666666666667), 0.0542410714286 * 6),
    PQP3(R3(0.5000000000000, 0.1666666666667, 0.1666666666667), 0.0542410714286 * 6),
    PQP3(R3(0.2500000000000, 0.2500000000000, 0.2500000000000), -0.0084656084656 * 6),
};
PQF3 const QuadratureFormular_GM_3D_S3(3, 35, P_GM_3D_S3);

#include "lex.hpp"
extern mylex *zzzfff;
static void Load_Init() {
    /*Global.New("gm2ds0", CConstant< const QuadratureFormular * >(&QuadratureFormular_GM_2D_S0));
    Global.New("gm2ds1", CConstant< const QuadratureFormular * >(&QuadratureFormular_GM_2D_S1));
    Global.New("gm2ds2", CConstant< const QuadratureFormular * >(&QuadratureFormular_GM_2D_S2));
    Global.New("gm2ds3", CConstant< const QuadratureFormular * >(&QuadratureFormular_GM_2D_S3));

    Global.New("gm3ds0", CConstant< const GQuadratureFormular< R3 > * >(&QuadratureFormular_GM_3D_S0));
    Global.New("gm3ds1", CConstant< const GQuadratureFormular< R3 > * >(&QuadratureFormular_GM_3D_S1));
    Global.New("gm3ds2", CConstant< const GQuadratureFormular< R3 > * >(&QuadratureFormular_GM_3D_S2));
    Global.New("gm3ds3", CConstant< const GQuadratureFormular< R3 > * >(&QuadratureFormular_GM_3D_S3));*/

    /*string tab[1000];
    for (int i=0; i<20;i++){
      tab[i]= ("gm2ds"+to_string(i));
      Global.New(tab[i].c_str(), CConstant< const QuadratureFormular * >(GenerateQuadratureFormular(i,2)));
    }*/
    for (int i = 0; i < 10; i++)
        TabQuadrFormula[i] = nullptr;
    Global.Add("GMQuadrature2D", "(",
               new OneOperator1s_<const QuadratureFormular *, long>(GenerateQuadratureFormularForOperator));
}

LOADFUNC(Load_Init)
// delete [] TabQuadrFormula;
