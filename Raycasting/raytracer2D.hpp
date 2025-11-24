#include<iostream>
#include<vector>
#include<cassert>
#include<cmath>
using namespace std;
typedef double real; 
typedef double vec2[2] ;

inline real dot(const vec2& a, const vec2& b) {
    return a[0]*b[0] + a[1]*b[1];
}

/*inline vec2 operator*(real s, const vec2& v) {
    return static{s*v[0], s*v[1]};
}*/

// Structure de donnees Ray
struct ray{
    vec2 origin;
    vec2 dir;
    bool norm1=false; 
    ray(vec2 o, vec2 d){ //: origin({0,0}), dir({0,0}) {}
        origin[0] = o[0]; origin[1] = o[1];
        dir[0] = d[0]; dir[1] = d[1];
    }
    void normaliseD(){
        if(!norm1){
            real norm=sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
            norm=1/norm;
            dir[0]*=norm;
            dir[1]*=norm;
            norm1=true;
        }
        else return;

    };
}  ;

// Structure domain
struct domain{
    real xmin,xmax,ymin,ymax;
    domain(real xm,real xM, real ym,real yM) : xmin(xm), xmax(xM),ymin(ym), ymax(yM) {}
}  ;

// Structure edge
struct edge{
    vec2 A,B; 
    edge(vec2 AA, vec2 BB) {        
        A[0] = AA[0]; A[1] = AA[1];
        B[0] = BB[0]; B[1] = BB[1];}
}  ;

// Orienter correctement la normale
inline void orientNormalForRay(const ray r, vec2 &normal) {
    // si la normale pointe dans le mauvais sens, on l'inverse
    if (dot(r.dir, normal) > 0.0) {
        normal[0] = -normal[0];
        normal[1] = -normal[1];
    }
}



// Ray tracer class
template<int dim>
class raytracer {
    private:
    domain omega; //  domain 
    vector<edge> edges; // edges list
    public: 
    raytracer(real xm,real xM, real ym,real yM,vector<edge> egdList): omega(xm,xM,ym,yM), edges(egdList){};
    void getNormal(edge AB, vec2 normale);
    bool intessectionpoint(ray r, edge AB,vec2 point );
    void reflection(ray &r, edge AB );
    void trappedOrNot(ray &r); // Trapped ray ?

};

