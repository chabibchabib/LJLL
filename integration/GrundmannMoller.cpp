#include "ff++.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

/* Generate the partitions of the number sum on nparts*/
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
/* Generate the weights and the integration points*/
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
            data.push_back({coef, beta}); // weight, point
        }
    }

    return data;
}

// Classe pour l'intégration 1D
class GrundmannMoller1D : public OneOperator {
  public:
    // Expression à intégrer
    class E_GMIntegral : public E_F0mps {
      public:
        Expression expr;  // Expression FreeFEM à intégrer
        Expression order; // Ordre s de la quadrature

        E_GMIntegral(const basicAC_F0 &args) {
            args.SetNameParam();
            expr = to<double>(args[0]);
            order = to<long>(args[1]);
        }

        AnyType operator()(Stack stack) const {
            long s = GetAny<long>((*order)(stack));
            int n = 1; // Dimension 1D

            // Récupérer les points de quadrature
            auto points = integration_weightspoints(n, s);

            double result = 0.0;

            // Créer un contexte pour évaluer l'expression
            MeshPoint &mp = *MeshPointStack(stack);
            R2 P_save(mp.P.x); // Sauvegarder les coordonnées

            for (const auto &point : points) {
                double weight = point.first;
                const auto &coords = point.second;

                // Transformer coordonnées barycentriques -> cartésiennes
                // Pour le triangle de référence: (0,0), (1,0), (0,1)
                double x = coords[0];
                // coords[2] est la troisième coordonnée barycentrique (implicite)

                mp.set(x);

                // Évaluer l'expression au point (x,y)
                double val = GetAny<double>((*expr)(stack));

                // Volume du simplexe de référence en 1D = 1/2
                result += weight * val;
            }

            mp.set(P_save.x); // Restaurer le point
            return result;
        }

        operator aType() const { return atype<double>(); }
    };

    E_F0 *code(const basicAC_F0 &args) const { return new E_GMIntegral(args); }

    GrundmannMoller1D()
        : OneOperator(atype<double>(),
                      atype<double>(), // Expression
                      atype<long>())   // Ordre
    {}
};

// Classe pour l'intégration 2D
class GrundmannMoller2D : public OneOperator {
  public:
    // Expression à intégrer
    class E_GMIntegral : public E_F0mps {
      public:
        Expression expr;  // Expression FreeFEM à intégrer
        Expression order; // Ordre s de la quadrature

        E_GMIntegral(const basicAC_F0 &args) {
            args.SetNameParam();
            expr = to<double>(args[0]);
            order = to<long>(args[1]);
        }

        AnyType operator()(Stack stack) const {
            long s = GetAny<long>((*order)(stack));
            int n = 2; // Dimension 2D

            // Récupérer les points de quadrature
            auto points = integration_weightspoints(n, s);

            double result = 0.0;

            // Créer un contexte pour évaluer l'expression
            MeshPoint &mp = *MeshPointStack(stack);
            R2 P_save(mp.P.x, mp.P.y); // Sauvegarder les coordonnées

            for (const auto &point : points) {
                double weight = point.first;
                const auto &coords = point.second;

                // Transformer coordonnées barycentriques -> cartésiennes
                // Pour le triangle de référence: (0,0), (1,0), (0,1)
                double x = coords[0];
                double y = coords[1];
                // coords[2] est la troisième coordonnée barycentrique (implicite)

                mp.set(x, y);

                // Évaluer l'expression au point (x,y)
                double val = GetAny<double>((*expr)(stack));

                // Volume du simplexe de référence en 2D = 1/2
                result += weight * val * 0.5;
            }

            mp.set(P_save.x, P_save.y); // Restaurer le point
            return result;
        }

        operator aType() const { return atype<double>(); }
    };

    E_F0 *code(const basicAC_F0 &args) const { return new E_GMIntegral(args); }

    GrundmannMoller2D()
        : OneOperator(atype<double>(),
                      atype<double>(), // Expression
                      atype<long>())   // Ordre
    {}
};

// Classe pour l'intégration 3D
class GrundmannMoller3D : public OneOperator {
  public:
    class E_GMIntegral3D : public E_F0mps {
      public:
        Expression expr;
        Expression order;

        E_GMIntegral3D(const basicAC_F0 &args) {
            args.SetNameParam();
            expr = to<double>(args[0]);
            order = to<long>(args[1]);
        }

        AnyType operator()(Stack stack) const {
            long s = GetAny<long>((*order)(stack));
            int n = 3; // Dimension 3D

            auto points = integration_weightspoints(n, s);

            double result = 0.0;
            MeshPoint &mp = *MeshPointStack(stack);
            R3 P_save(mp.P.x, mp.P.y, mp.P.z); // Sauvegarder les coordonnées

            for (const auto &point : points) {
                double weight = point.first;
                const auto &coords = point.second;

                // Coordonnées cartésiennes pour le tétraèdre de référence
                double x = coords[0];
                double y = coords[1];
                double z = coords[2];
                // coords[3] est implicite

                mp.set(x, y, z);
                double val = GetAny<double>((*expr)(stack));

                // Volume du tétraèdre de référence = 1/6
                result += weight * val / 6.0;
            }

            mp.set(P_save.x, P_save.y, P_save.z); // Restaurer le point
            return result;
        }

        operator aType() const { return atype<double>(); }
    };

    E_F0 *code(const basicAC_F0 &args) const { return new E_GMIntegral3D(args); }

    GrundmannMoller3D() : OneOperator(atype<double>(), atype<double>(), atype<long>()) {}
};

// Fonction d'initialisation du plugin
static void Init() {
    // Enregistrer les fonctions
    Global.Add("intGM1d", "(", new GrundmannMoller1D);
    Global.Add("intGM2d", "(", new GrundmannMoller2D);
    Global.Add("intGM3d", "(", new GrundmannMoller3D);
}

LOADFUNC(Init);