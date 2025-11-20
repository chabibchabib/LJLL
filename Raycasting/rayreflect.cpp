#include<iostream>
#include<vector>
#include<cassert>
#include<cmath>
using namespace std;
typedef double real; 

// Un rayon est definie par une origine
//et un vecteur directeur R=Origin+ t*Dir
template<int dim>
struct ray{
    vector<real> origin;
    vector<real> dir; 
    ray() : origin(dim), dir(dim) {}
}  ;

// Un cercle est defini
// par son Centre et
//son rayon 
template<int dim>
struct cercle{
    vector<real> center;
    real radius;
    cercle() :  center(dim) {}
}  ;

// Un rectangle est defini
// par ses points 
 
//
template<int dim>
struct rectangle{
    vector<real> P0;
    vector<real> P1;
    vector<real> P2;

    real width;
    real height;

    void computeDim(){
        width=sqrt(pow((P1[0]-P0[0]),2) + pow((P1[1]-P0[1]),2));
        height=sqrt(pow((P2[0]-P0[0]),2) + pow((P2[1]-P0[1]),2));
    }

    bool pointonsegment(const vector<real>& P,
                        const vector<real>& A,
                        const vector<real>& B,
                        real eps = 1e-12)
    {
        real dx = B[0] - A[0];
        real dy = B[1] - A[1];
        real dxp = P[0] - A[0];
        real dyp = P[1] - A[1];
        // colineaire ou non?
        real cross = dx*dyp - dy*dxp;
        if(fabs(cross) > eps) return false;

        real dot = dx*dxp + dy*dyp;
        if(dot < -eps) return false;
        // hors segment?
        real len2 = dx*dx + dy*dy;
        if(dot > len2 + eps) return false;

        return true;
    }

    bool pointonedge(const vector<real>& P,real eps = 1e-12)
    {
        vector<real> P3 = { P1[0] + (P2[0]-P0[0]),
                            P1[1] + (P2[1]-P0[1]) };

        return pointonsegment(P,P0,P1,eps) ||
            pointonsegment(P,P1,P3,eps) ||
            pointonsegment(P,P3,P2,eps) ||
            pointonsegment(P,P2,P0,eps);
    }


    vector<real> getnormal(vector<real> P){
        assert(P.size()==dim);
        assert(pointonedge(P));
        vector<real> A(dim); 
        vector<real> B(dim); 
        vector<real> P3 = { P1[0] + (P2[0]-P0[0]),
                            P1[1] + (P2[1]-P0[1]) };
        real eps = 1e-12;
        if(pointonsegment(P,P0,P1,eps)) {A=P0; B=P1;}
        if(pointonsegment(P,P1,P3,eps)) {A=P1; B=P3;}
        if(pointonsegment(P,P3,P2,eps)) {A=P3; B=P2;}
        if(pointonsegment(P,P2,P0,eps)) {A=P2; B=P0;}

        real dx = B[0] - A[0];
        real dy = B[1] - A[1];
        vector<real> normal = {-dy, dx};
        real norm = 1/sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
        normal[0]*=norm;
        normal[1]*=norm;

        return normal;
    }

    rectangle(vector<real> P00, vector<real> P11, vector<real> P22) :  P0(P00),P1(P11),P2(P22) {}
}  ;

int main(){
    ray<2> r;

    /*r.origin={0,0};
    r.dir={1,2};*/
    vector<real> P0={0,0};
    vector<real> P1={1,0};
    vector<real> P2={0,3};
    rectangle<2> Rec(P0,P1,P2);
    Rec.computeDim();
    cout<<"HxW= "<<Rec.height<<" , "<<Rec.width<<endl;
    vector<real> n=Rec.getnormal({1,2});
    cout<<"normal= "<<n[0]<<" , "<<n[1]<<endl;
    return 0;
}