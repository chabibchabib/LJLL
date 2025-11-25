#include"raytracer2D.hpp"
void coutvec(string MSG, vec2 v){
    cout<<MSG+" ="<<v[0]<<" "<<v[1]<<endl;
}
// Calculer la normale d'un edge
template<int dim>
void raytracer<dim>::getNormal(edge AB, vec2 normale){
        real eps = 1e-12;
        real dx = AB.B[0] - AB.A[0];
        real dy = AB.B[1] - AB.A[1];
        
        normale[0]= -dy;
        normale[1] = dx;
        real norm = 1/sqrt(normale[0]*normale[0] + normale[1]*normale[1]);
        normale[0]*=norm;
        normale[1]*=norm;
} 

template<int dim>
bool raytracer<dim>::intessectionpoint(ray r, edge AB,vec2 point ){
    //Crammer's rule
    real ABx = (AB.B[0]-AB.A[0]), ABy= (AB.B[1]-AB.A[1]), AOx=(r.origin[0]-AB.A[0]),AOy=(r.origin[1]-AB.A[1]);

    real det=-ABx*r.dir[1] + ABy*r.dir[0];
    assert(abs(det)> 1e-12);
    if(abs(det)< 1e-5) return false;
    det=1/det;
    real s= det*(r.dir[0]*AOy - r.dir[1]*AOx );
    cout<<"s="<<s<<"det="<< 1/det<<endl;
    //assert(s>=0 && s<=1); // appartient au segment?
    if(s< 1e-12 || s>1) return false;

    real t= det*(ABx*AOy - ABy*AOx );
    //assert(t>=1e-12); // appartient au faisceau?
    if(t< 1e-5 ) return false;

    point[0]= r.origin[0]+ t*r.dir[0];
    point[1]= r.origin[1]+ t*r.dir[1];
    cout<<"Test ="<<AB.A[0]+s*(AB.B[0]-AB.A[0])<<","<<AB.A[1]+s*(AB.B[1]-AB.A[1])<<endl;
    assert(abs(AB.A[0]+s*(AB.B[0]-AB.A[0])-point[0])<1e-12);
    assert(abs(AB.A[1]+s*(AB.B[1]-AB.A[1])-point[1])<1e-12);
    return true;
}

template<int dim>
void raytracer<dim>::reflection(ray &r, edge AB ){
        vec2 normale;
        getNormal(AB,normale);
        orientNormalForRay(  r,  normale);
        vec2 pt;
        intessectionpoint(r,  AB,pt );
        copy( pt, pt+dim, r.origin);
        real ddotn=dot(normale, r.dir);
        r.dir[0]-=2*ddotn*normale[0];
        r.dir[1]-=2*ddotn*normale[1];
        r.normaliseD();

    }

template<int dim>
void raytracer<dim>::trappedOrNot(ray &r, int & hit,int hitmax){
    hit=0.;
    while(hit<hitmax){
        cout<<"\nr="<<r.origin[0]<<","<<r.origin[1]<<endl;

        vector<real> params;
        vector<real> idx;
        int i=0;
        for (auto &E: edges){
            vec2 pt;
            if(intessectionpoint( r,  E, pt )){
                real d = (r.dir[0]!=0) ?  r.dir[0]: r.dir[1];
                real o = (d==r.dir[0]) ?  r.origin[0]: r.origin[1];
                real coord = (d==r.dir[0]) ?  pt[0]: pt[1];
                idx.push_back(i);
                params.push_back((coord-o)/d);
            };
            i++;
        }
        cout<<"params.empty()? "<<params.size()<<" "<<params.empty()<<endl;
        if(!params.empty()){

            sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2) {
                return params[i1] < params[i2];  // ordre croissant
            });

            reflection( r,  edges[idx[0]] );
            r.normaliseD();
            if(r.origin[0]<omega.xmin ||r.origin[0]>omega.xmax || r.origin[1]<omega.ymin ||r.origin[1]>omega.ymax  ) break;
            hit++;
        }
        else break;
        }

}
int main(){
    domain omega(0,6,0,6); // domaine carree
    vec2 A={2,1}, B={2,5};
    edge E(A,B);
    vector<edge> edgeList;
    edgeList.push_back(E);
    A[0]=4,A[1]=1;
    B[0]=4,B[1]=5;
    edgeList.push_back(edge (A,B));
    /*for (auto elm: edgeList){
        cout<<"elm="<<elm.B[0]<<","<<elm.B[1];
    }*/
    raytracer<2> Rt(0,6,0,6,edgeList);
    vec2 normal;
    Rt.getNormal(E,normal);
    vec2 O={3,2}, dir={0.5,0.25};
    ray r(O,dir);
    r.normaliseD();
    int hit;
    Rt.trappedOrNot(r,  hit,1000);
    cout<<"HIT=" <<hit<<endl;
    /*orientNormalForRay( r, normal);
    coutvec("normal", normal);
    
    // Intersection point
    vec2 point ;
    Rt.intessectionpoint( r,  E, point );
    coutvec("intessectionpoint", point);
    
    coutvec("Ancienne origine ", r.origin);
    coutvec("Ancienne direction ", r.dir);
    
    Rt.reflection(r,E);
    coutvec("nouvelle origine ", r.origin);
    coutvec("nouvelle direction ", r.dir);*/

    return 0;
}