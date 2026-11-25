

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

// a ray is defined by its origin and direction
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


void reflection(ray &r, R2 &normale){
    /* This function computes the reflection of a 
    ray direction on a surface given its normal vector.
    */
    // the normal oposes the direction?
    R ddotn=(r.dir[0]*normale[0] +r.dir[1]*normale[1]);
    // if its not change the normal direction
    if ( ddotn> 0.0) {
        normale[0] = -normale[0];
        normale[1] = -normale[1];
        ddotn*=-1;
    }
    // Refclect then
    r.dir[0]-=2*ddotn*normale[0];
    r.dir[1]-=2*ddotn*normale[1];
    R norm = sqrt(r.dir[0] * r.dir[0] + r.dir[1] * r.dir[1]);
    r.dir[0]/=norm;
    r.dir[1]/=norm;
    ffassert(r.dir[0] || r.dir[1] );
    }
    
static KN<long> TabLabel; // This array will contain the labels of boundaries

// Principal fct 
long Raycaster(Fem2D::Mesh const *const &pTh, const long & itt,  KN<R> *const &RayRef, KN<double> * const &HitsRay)  {
    /* Pth:Mesh, itt: nbr triangle, RayRef: array containing the origin and the dir of the ray
        HitsRay: the output an array, containing the number of time each element hit by a ray  
    */

    int it=itt;
    R2 normale;
    ray Ray(R2((*RayRef)[0],(*RayRef)[1]),R2((*RayRef)[2],(*RayRef)[3]));
   
    static int count = 0;
    double dt= 200; // A big step used to walk in the mesh
    
    R l[3];
    l[1]=Ray.origin[0]  ;
    l[2]=Ray.origin[1]  ;
    l[0]=1-l[1]-l[2]  ;
    
    int k = 0;
    int j=0;
    int counter=0;

    while (j >= 0) { 
    
        j = WalkInTriangle(*pTh, it, l, Ray.dir[0], Ray.dir[1], dt);
        const Triangle & T = (*pTh)[it];
        const R2 Q[3]={(const R2) T[0],(const R2) T[1],(const R2) T[2]};
        Ray.origin= l[0]*Q[0]  + l[1]*Q[1]  + l[2]*Q[2]; // In the Cartisian global coordinates
        int oldj=j;
        int itt = (*pTh).ElementAdj(it, j); // Get adj neighboor, it changes the value of j, hence oldj 
        
        long labR=-1; // used to store the labels of boundaries
        if(( (itt==it) )  )  labR=TabLabel[3*it+j]; // RBC
            // Reflect 
        double eps= 1e-12;
        if(!(labR%376)){
            #pragma omp atomic update
            counter++;
            normale=normalVector( oldj,T);
            reflection(Ray, normale);
            Ray.origin=Ray.origin+ Ray.dir * 1e-9;
            continue;
        }
        dt=100;
        ffassert(abs(l[oldj])<1e-10); //Check whether the barycentric coordinate associated with node j on the element boundary is zero. 
        // int jj  = j;
        R a = l[(oldj + 1) % 3], b = l[(oldj + 2) % 3]; // Update for the neighboor elem
        if (itt == it || itt < 0)
            break; // Border
        it = itt;
        l[j] = 0;
        l[(j + 1) % 3] = b;
        l[(j + 2) % 3] = a;
        if (++k > 1e+05) {
            //cerr << "Trapping zone "<< endl;
            break;
        }
    }

    (*HitsRay)[itt]=(*HitsRay)[itt]+counter;
    return 1;
}

inline R2 GetCercleAngle(R k, R n){
    ffassert(abs(k)<=n);
    R coef=k/(n)*Pi;
    return R2(cos(coef),sin(coef));
}

long MeshRaycaster(Fem2D::Mesh const *const &pTh,  const long & NbrRays, KN<double> * const &HitsRay)  {

    TabLabel.resize((*pTh).nt *3); // Allocate 3*nbr of triangles
    // Get the labels
    for (int cpt =0; cpt <(*pTh).neb; cpt++){
        int ed; 
        int elm=(*pTh).BoundaryElement(cpt,ed);
        TabLabel[3*elm+ed]= (*pTh).be(cpt).lab;
    } 
    // Loop 
    int HalfNbrRays=(NbrRays>2)? NbrRays/2 : 1;
    #pragma omp parallel for //num_threads(8)
    for (int it =0; it<pTh->nt;it++) { // parallelisable
        for (int angle= -(NbrRays-1);angle<=NbrRays;angle++){ // Parallelisable
            R2 Dir=GetCercleAngle( angle, NbrRays); 
            // Launch from the elm center
            KN<R> Ray(4);
            Ray[0]=1./3;
            Ray[1]=1./3;
            Ray[2]=Dir[0];
            Ray[3]=Dir[1];
            // The raytacer/caster
            Raycaster(pTh,it,  &Ray, HitsRay); // /!\ condition race
        }
        #pragma omp atomic
        (*HitsRay)[it]/=double( ( (*pTh)[it]).area);
    }
    return 1;
}


long RaycasterVisibleArea(Fem2D::Mesh const *const &pTh, const long & itt,  KN<R> *const &RayRef,KN<double> * const &visibleelm)  {
    /*
    Function similar to the previous Raycaster fct, but used to detect 
    the visible/invisible zones.
    This time only the elements living in the Trapping zone have the right to launch the rays 
    Visible elem=1/ invisible elem=0
    */
    int it=itt;
    R2 normale;
    ray Ray(R2((*RayRef)[0],(*RayRef)[1]),R2((*RayRef)[2],(*RayRef)[3]));
    double dt= 200;
    R l[3];
    l[1]=Ray.origin[0]  ;
    l[2]=Ray.origin[1]  ;
    l[0]=1-l[1]-l[2]  ;
    int k = 0;
    int j=0;
    int counter=0;
    while (j >= 0) { 
        j = WalkInTriangle(*pTh, it, l, Ray.dir[0], Ray.dir[1], dt);
        int oldj=j;
        int itt = (*pTh).ElementAdj(it, j);
        const Triangle & T = (*pTh)[it];
        const R2 Q[3]={(const R2) T[0],(const R2) T[1],(const R2) T[2]};
        Ray.origin= l[0]*Q[0]  + l[1]*Q[1]  + l[2]*Q[2];

        double eps= 1e-12;
        long labR=-1; // label edge
        if(( (itt==it) )  )  labR=TabLabel[3*it+j]; //RBC
        if (!(labR%376)){
            normale=normalVector( oldj,T);
            reflection(Ray, normale);
            Ray.origin=Ray.origin+ Ray.dir * 1e-9;
            continue;

        }
        dt=100;
        ffassert(abs(l[oldj])<1e-10);
        R a = l[(oldj + 1) % 3], b = l[(oldj + 2) % 3];
        // MAJ mapping HitsRay
        if (itt == it || itt < 0)
            break; // le bord
        it = itt;
        int old_val;
        if (!(*visibleelm)[it]) {
            #pragma omp atomic write
            (*visibleelm)[it] = 1;
        }
        l[j] = 0;
        l[(j + 1) % 3] = b;
        l[(j + 2) % 3] = a;
        if (++k > 1e+05) {
            //cerr << "Trapping zone "<< endl;
            break;
        }

    }
    return 1;
}

long  MeshRaycasterVisibleArea(Fem2D::Mesh const *const &pTh,  const long & NbrRays, KN<double> * const &HitsRayKN, const long &threshold, KN<double> * const &visibleelm)  {
    #pragma omp parallel for //num_threads(8)
    for (int it =0; it<pTh->nt;it++) { // parallelisable
        if((*HitsRayKN)[it]>threshold){
            for (int angle= -(NbrRays-1);angle<=NbrRays;angle++){ // Parallelisable
            R2 Dir=GetCercleAngle( angle, NbrRays); 
            KN<R> Ray(4);
            Ray[0]=1./3;
            Ray[1]=1./3;
            Ray[2]=Dir[0];
            Ray[3]=Dir[1];
            RaycasterVisibleArea(pTh, it,  &Ray,visibleelm);

        }
        }

        else continue; 
    }
    return 1;
}
static void Load_Init( ) {    // le constructeur qui ajoute la fonction "Raycaster"  a freefem++
  if (verbosity) {
    cout << " load: Raycaster  " << endl;
  }

  Global.Add("MeshRaycaster", "(", new OneOperator3_< long, pmesh, long, KN<double>* >(MeshRaycaster));
  Global.Add("MeshRaycasterVisibleArea", "(", new OneOperator5_< long, pmesh, long, KN<double>*, long,KN<double>*  >(MeshRaycasterVisibleArea));

    cout<<"END load: Raycaster"<<endl;
}
LOADFUNC(Load_Init)
