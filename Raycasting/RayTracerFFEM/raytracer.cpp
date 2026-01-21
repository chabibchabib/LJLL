

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
#include<omp.h>

using namespace std;

struct ray {
    R2 origin={0,0};
    R2 dir={0,0};
    ray(R2 O, R2 D){
        origin[0]=O[0];
        origin[1]=O[1];
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
    //cout<<"inside reflection fct\n";
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
    //cout<<"values ray "<<r.dir[0]<<" "<<r.dir[1]<<" values normale "<<normale[0]<<" "<<normale[1]<<endl;
    ffassert(r.dir[0] || r.dir[1] );

    }

long RayTracer(Fem2D::Mesh const *const &pTh, const long & itt,  KN<R> *const &RayRef, KN<long> * const &HitsRay)  {
    int it=itt;
    R2 normale;
    ray Ray(R2((*RayRef)[0],(*RayRef)[1]),R2((*RayRef)[2],(*RayRef)[3]));
    //HitsRay = 0; // Initialise tout à zéro
    static int count = 0;
    double dt= 200;
    R l[3];
    l[1]=Ray.origin[0]  ;
    l[2]=Ray.origin[1]  ;
    l[0]=1-l[1]-l[2]  ;
    //if((l[1]< 1e-12)  || (l[1]>1.) || (l[2]< 1e-12) || (l[2]> 1.) || l[1] +l[2]>1. ) {cout<<"coord baryc erronee"<<endl; exit(0);}

    /*cout<<"Ray=("<<Ray.origin[0]<<","<<Ray.origin[1]<<" ),("<<Ray.dir[0]<<","<<Ray.dir[1]<<")"<<endl;
    cout<<"l=("<<l[0]<<","<<l[1]<<" ,"<<l[2]<<")"<<endl;*/
    int k = 0;
    int j=0;
    while (j >= 0) { 
        /*cout<<"\t triangle= "<<it<<" k="<<k<<endl;
        cout<<"l=("<<l[0]<<","<<l[1]<<" ,"<<l[2]<<")"<<endl;*/
        j = WalkInTriangle(*pTh, it, l, Ray.dir[0], Ray.dir[1], dt);

        const Triangle & T = (*pTh)[it];
        const R2 Q[3]={(const R2) T[0],(const R2) T[1],(const R2) T[2]};
        Ray.origin= l[0]*Q[0]  + l[1]*Q[1]  + l[2]*Q[2];
        /*Ray.origin[0] = l[1];
        Ray.origin[1] = l[2];*/

        // Reflect 
        //if((Ray.origin[0]>= 0.25 && Ray.origin[0]<= 0.75) || (Ray.origin[1]>= 0.25 && Ray.origin[1]<= 0.75)){
        if((abs(Ray.origin[1])<1e-8 || abs(Ray.origin[1]-1)<1e-8) && ( abs(Ray.origin[0])<0.5) ){
        //if(0){
            #pragma omp atomic update
            (*HitsRay)[it]=(*HitsRay)[it]+1;
            //out<<"Im on borders\n";
            normale=normalVector( j,T);
            /*ccout<<"Normale sur triangle="<<it<<" Edge= "<<j<<" Normale= "<<normale[0]<<" , "<<normale[1]<<endl;
            cout<<"Ray before=("<<Ray.origin[0]<<","<<Ray.origin[1]<<" ),("<<Ray.dir[0]<<","<<Ray.dir[1]<<")"<<endl;*/
            
            reflection(Ray, normale);
            //cout<<"Ray after=("<<Ray.origin[0]<<","<<Ray.origin[1]<<" ),("<<Ray.dir[0]<<","<<Ray.dir[1]<<")"<<endl;
            Ray.origin=Ray.origin+ Ray.dir * 1e-9;

            continue;

        }


        /*cout<<"Q=("<<Q[0][0]<<","<<Q[0][1]<<"),("<<Q[1][0]<<","<<Q[1][1]<<"),("<<Q[2][0]<<","<<Q[2][1]<<")"<<endl;
        cout<<"j="<<j<<" (j + 1) % 3="<<(j + 1) % 3<<" (j + 2) % 3="<<(j + 2) % 3<<" itt="<<itt<<endl;
        cout<<"l=("<<l[0]<<","<<l[1]<<" ,"<<l[2]<<")"<<endl;
        cout<<"Ray=("<<Ray.origin[0]<<","<<Ray.origin[1]<<" ),("<<Ray.dir[0]<<","<<Ray.dir[1]<<")"<<endl;
        cout<<"DT="<<dt<<endl;*/
        dt=100;
        ffassert(abs(l[j])<1e-10);
        // int jj  = j;
        R a = l[(j + 1) % 3], b = l[(j + 2) % 3];
        // MAJ mapping HitsRay


        int itt = (*pTh).ElementAdj(it, j);
        if (itt == it || itt < 0)
            break; // le bord
        it = itt;
        l[j] = 0;
        l[(j + 1) % 3] = b;
        l[(j + 2) % 3] = a;
        #pragma omp atomic update
        (*HitsRay)[it]=(*HitsRay)[it]+1;
        //k=k+1;
        if (++k > 10000) {
            cerr << "Trapping zone "<< endl;
            break;
        }

    }

    return 1;
}

inline R2 DiscretizeCercleAngle(R k, R n){
    ffassert(abs(k)<=n);
    R coef=k/(n)*Pi;
    return R2(cos(coef),sin(coef));
}

long MeshRayTracer(Fem2D::Mesh const *const &pTh,  const long & NbrRays, KN<long> * const &HitsRay)  {
    #pragma omp parallel for num_threads(8)
    for (int it =0; it<pTh->nt;it++) { // parallelisable
            cout<<"Num threads="<<omp_get_num_threads()<<endl;
        for (int angle= -(NbrRays-1);angle<=NbrRays;angle++){ // Parallelisable
            R2 Dir=DiscretizeCercleAngle( angle, NbrRays); 
            KN<R> Ray(4);
            Ray[0]=1./3;
            Ray[1]=1./3;
            Ray[2]=Dir[0];
            Ray[3]=Dir[1];
            RayTracer(pTh,it,  &Ray, HitsRay); // /!\ condition race
        }
    }
    return 1;
}
static void Load_Init( ) {    // le constructeur qui ajoute la fonction "raytracer"  a freefem++
  if (verbosity) {
    cout << " load: raytracer  " << endl;
  }

  Global.Add("RayTracer", "(", new OneOperator4_< long, pmesh, long,  KN<R>*, KN<long>* >(RayTracer));
  Global.Add("MeshRayTracer", "(", new OneOperator3_< long, pmesh, long, KN<long>* >(MeshRayTracer));

    cout<<"END load: raytracer"<<endl;
}
LOADFUNC(Load_Init)
