

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include "RNM.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "ff++.hpp"
#include "AFunction_ext.hpp" // Extension of "AFunction.hpp" to deal with more than 3 parameters function


using namespace std;

struct ray {
    R2 origin={0,0};
    R2 dir={0,0};
    ray(R2 O, R2 D){
        origin[0]=O[0];
        origin[0]=O[1];
        dir[0]=D[0];
        dir[1]=D[1];
    }
}; 

R2 normalVector(int j,const Triangle &T){
    /* Compute the normal of the edge
     opposed to vertice j in a given triangle T */
    R2 normale={0,0};
    R2 A = T[(j + 1) % 3]; // Node edge 1
    R2 B = T[(j + 2) % 3]; // Node edge 1
    normale[0]=B.y - A.y; // Rotation 90° (dy,-dx)
    normale[1]=-(B.x - A.x);
    R norm = sqrt(normale[0] * normale[0] + normale[1] * normale[1]);
    normale[0]/=norm;
    normale[1]/=norm;
    return normale;
}

// Orienter correctement la normale
/*inline void orientNormalForRay(const ray r, vector<REAL> normal) {
    // si la normale pointe dans le mauvais sens, on l'inverse
    if ((r.dir[0]*normal[0] +r.dir[1]*normal[1]) > 0.0) {
        normal[0] = -normal[0];
        normal[1] = -normal[1];
    }
}*/

void reflection(ray &r, R2 &normale){
        R ddotn=(r.dir[0]*normale[0] +r.dir[1]*normale[1]);
        if ( ddotn> 0.0) {
            normale[0] = -normale[0];
            normale[1] = -normale[1];
            ddotn*=-1;
        }
        r.dir[0]-=2*ddotn*normale[0];
        r.dir[1]-=2*ddotn*normale[1];
        R norm = sqrt(r.dir[0] * r.dir[0] + r.dir[1] * r.dir[1]);
        r.dir[0]/=norm;
        r.dir[1]/=norm;
    }

long RayTracer(Fem2D::Mesh const *const &pTh, const long & itt,  KN<R> *const &RayRef, KN<long> * const &HitsRay)  {
    int it=itt;
    ray Ray(R2((*RayRef)[0],(*RayRef)[1]),R2((*RayRef)[2],(*RayRef)[3]));
    //HitsRay = 0; // Initialise tout à zéro
    static int count = 0;
    static R ddtp = 0;
    double dt= 1000;
    R ddt = dt;
    ddtp = ddt;
    R l[3];
    int k = 0;
    int j;
    while ((j = WalkInTriangle(pTh[it], it, l, Ray.dir[0], Ray.dir[1], dt)) >= 0) { // Dans la fct originale  u v sont maj grace à la ligne mpc.change
        cout<<"j="<<j<<" (j + 1) % 3="<<(j + 1) % 3<<" (j + 2) % 3="<<(j + 2) % 3<<" itt="<<itt<<endl;

        // MAJ origin ray
        Ray.origin[0]=l[1];
        Ray.origin[1]=l[2];
        ffassert(l[j] == 0);
        // int jj  = j;
        R a = l[(j + 1) % 3], b = l[(j + 2) % 3];
        int itt = pTh[it].ElementAdj(it, j);
        if (itt == it || itt < 0)
            break; // le bord
        // Reflection si bord: update Ray (To do ajout d'une condition de reflection)
        //vector<REAL> normale= normalVector( j,*(mpc.T));
        //reflection(Ray, normale);
        // MAJ mapping HitsRay
        HitsRay[it]+=1; 
        it = itt;
        l[j] = 0;
        l[(j + 1) % 3] = b;
        l[(j + 2) % 3] = a;

        /*if (k++ > 1000) {
            cerr << "Fatal  error  in RayTracer (R2) operator: loop  => velocity too high ???? or NaN "<< endl;
            ffassert(0);
        }*/
    }

    return 1;
}

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "raytracer"  a freefem++
  if (verbosity) {
    cout << " load: raytracer  " << endl;
  }

  Global.Add("RayTracer", "(", new OneOperator4_< long, pmesh, long,  KN<R>*, KN<long>* >(RayTracer));

    cout<<"END load: raytracer"<<endl;
}
LOADFUNC(Load_Init)
