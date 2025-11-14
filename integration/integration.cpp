#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip> // std::setprecision
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

using namespace std;
double factorial(int nbr) {
    if (nbr <= 1)
        return 1.;
    else {
        double prod = 1.;
        for (int i = 1; i <= nbr; i++) {
            prod *= i;
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

double integration_rule(int n, int s, function<double(const vector<double> &)> f) {
    int d = 2 * s + 1;
    double result = 0.;
    for (int i = 0; i < s + 1; i++) {
        double coef = pow(-1, i) * pow(2, -2 * s) * pow(d + n - 2 * i, d) /
                      (factorial(i) * factorial(d + n - i)); // C=(-1)^i x 2^(-2s) x (d+n-2i)^d/(i! x (d+n-i) !)

        vector<vector<double>> betas = generate_partitions(s - i, n + 1);

        double denom = d + n - 2 * i;
        for (const auto &beta_orig : betas) {
            vector<double> beta;
            for (auto x : beta_orig) {
                // x = (2*x +1)/denom;
                beta.push_back((2 * x + 1) / denom);
            }

            /*cout << "i=" << i << ", beta=(";
            for(auto xi : beta) cout << xi << " ";
            cout << "), coef=" << coef << endl;*/

            result += coef * f(beta);
        }
    }
    return result;
}

double fct(vector<double> x) {
    /*double y=.0;
    int idx=0;
    for(auto & xx : x) { y+=pow(xx,2); }//y[idx]=pow(xx,2); idx++;}*/
    double y = 1.;
    /*for (auto it = x.begin(); it != x.end() -1; it++) {
        y += pow(*it, 2);
    }*/
    // y =21* pow(x[0], 1)+ pow(x[0], 2)-3;//+pow(x[1],2)+2*x[0]*x[1]+pow(x[2],2);
    return y;
};

double evaluate_integrale(const vector<pair<double, vector<double>>> &data, double (*fct)(vector<double> x)) {
    double result = 0.;
    double sum = 0;
    for (const auto &pair : data) {
        const auto &weight = pair.first;
        const auto &point = pair.second;
        double y = fct(point);
        if (abs(weight) > 1e-04) {
            result += weight * y;
            sum += weight;
        }
    }
    cout << scientific << setprecision(15);
    cout << "Poids sum= " << sum << endl;
    return result;
}

vector<pair<double, vector<double>>> integration_weightspoints(int n, int s) {
    int d = 2 * s + 1;
    vector<pair<double, vector<double>>> data;

    for (int i = 0; i < s + 1; i++) {
        double coef = pow(-1, i) * pow(2, -2 * s) * pow(d + n - 2 * i, d) / (factorial(i) * factorial(d + n - i));
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

int nbrOfIntgPOints(int n, int s) {
    int nbr = 0;
    for (int i = 0; i <= s; i++)
        nbr += factorial(n + s - i) / (factorial(s - i) * factorial(n));
    return nbr;
}
int main(int argc, char *argv[]) {
    int n = atoi(argv[1]), s = atoi(argv[2]);

    auto f = [](vector<double> x) {
        return pow(x[0], 0); //+pow(x[1],2)+2*x[0]*x[1]+pow(x[2],2);
        // return 1.;
        // return x[0];
    };
    /*double result = integration_rule( n, s,f);
    cout<<"result:"<<result<<endl;*/

    auto data = integration_weightspoints(n, s);
    double result2 = evaluate_integrale(data, fct);
    cout << scientific << setprecision(15);

    cout << "result2:" << 1 - result2 << endl;

    /*for (auto & element : data){
        //cout<<element.first<<"\t"<<element.second[0]<<endl;
        cout<<endl;
        cout<<element.first<<"\t";
        for(auto &val : element.second)
            cout<<val<<"\t";

        }*/
    cout << endl;
    ofstream file("output.txt");
    if (!file) {
        cerr << "Error file" << endl;
        return 1;
    }

    // file << fixed << setprecision(13);
    file << scientific << setprecision(15);

    for (auto &element : data) {
        // file<<endl;
        file << "PQP3(R3(";
        for (auto it = element.second.begin(); it != element.second.end() - 1; it++)

        {
            if (it != element.second.end() - 2)
                file << *it << ",";
            else
                file << *it << "),";
        }
        // file<<"),";
        file << element.first << "*6)," << endl;
    }
    file.close();
    cout << "N= " << n << " s= " << s << " Nbr Points= " << nbrOfIntgPOints(n, s) << endl;
    return 0;
}
