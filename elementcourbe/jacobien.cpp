#include<iostream>
#include<cmath>
#include<vector>

struct node{
    double x, y;
};

void shapefct( double zeta, double eta, std::vector<double> &N ){
    N[0]= 1.-zeta-eta;
    N[1]= zeta;
    N[2]= eta;
}

void dshapefct( std::vector<double> &dN_dzeta, std::vector<double> &dN_deta ){

    dN_dzeta[0]= -1;
    dN_dzeta[1]= 1.;
    dN_dzeta[2]= 0;
    
    dN_deta[0]= -1;
    dN_deta[1]= 0;
    dN_deta[2]= 1;

}

void jacobian (std::vector<node> Phycoordinayes ){
    std::vector<double> dN_dzeta(3,0), dN_deta(3,0);
    dshapefct( dN_dzeta,  dN_deta );
    std::vector<std::vector<double>> J(2,std::vector<double>(2,0));
    for (int i =0;i<Phycoordinayes.size(); i++){
        J[0][0]+=dN_dzeta[i]*Phycoordinayes[i].x;
        J[0][1]+=dN_deta[i]*Phycoordinayes[i].x;
        J[1][0]+=dN_dzeta[i]*Phycoordinayes[i].y;
        J[1][1]+=dN_deta[i]*Phycoordinayes[i].y;
    }
    double detJ=J[0][0]*J[1][1]-J[0][1]*J[1][0];
    std::cout << "J = \n";
    std::cout << J[0][0] << " " << J[0][1] << "\n";
    std::cout << J[1][0] << " " << J[1][1] << "\n";
    std::cout << "det(J) = " << detJ << std::endl;    
}

int main() {
    /*std::vector<node> nodes = {
        {0, 0},
        {2, 0},
        {1, 1}
    };*/

    std::vector<node> nodes = {
        {0, 0},
        {1, 0},
        {0, 1}
    };


    jacobian(nodes);
}