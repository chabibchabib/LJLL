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
void raytracer<dim>::intessectionpoint(ray r, edge AB,vec2 point ){
    //Crammer's rule
    real ABx = (AB.B[0]-AB.A[0]), ABy= (AB.B[1]-AB.A[1]), AOx=(r.origin[0]-AB.A[0]),AOy=(r.origin[1]-AB.A[1]);

    real det=-ABx*r.dir[1] + ABy*r.dir[0];
    assert(abs(det)> 1e-12);
    det=1/det;
    real s= det*(r.dir[0]*AOy - r.dir[1]*AOx );
    //cout<<"s="<<s<<"det="<< 1/det<<endl;
    assert(s>=0 && s<=1); // appartient au segment?
    real t= det*(ABx*AOy - ABy*AOx );
    assert(t>=1e-12); // appartient au faisceau?

    point[0]= r.origin[0]+ t*r.dir[0];
    point[1]= r.origin[1]+ t*r.dir[1];
    //cout<<"Test ="<<AB.A[0]+s*(AB.B[0]-AB.A[0])<<","<<AB.A[1]+s*(AB.B[1]-AB.A[1])<<endl;
    assert(abs(AB.A[0]+s*(AB.B[0]-AB.A[0])-point[0])<1e-12);
    assert(abs(AB.A[1]+s*(AB.B[1]-AB.A[1])-point[1])<1e-12);

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


int main(){
    domain omega(0,6,0,6); // domaine carree
    vec2 A={0,0}, B={1,0};
    edge E(A,B);
    vector<edge> edgeList;
    edgeList.push_back(E);
    raytracer<2> Rt(0,6,0,6,edgeList);
    vec2 normal;
    Rt.getNormal(E,normal);
    vec2 O={0.5,-1}, dir={0,2};
    ray r(O,dir);
    r.normaliseD();

    orientNormalForRay( r, normal);
    coutvec("normal", normal);
    
    // Intersection point
    vec2 point ;
    Rt.intessectionpoint( r,  E, point );
    coutvec("intessectionpoint", point);
    
    coutvec("Ancienne origine ", r.origin);
    coutvec("Ancienne direction ", r.dir);
    
    Rt.reflection(r,E);
    coutvec("nouvelle origine ", r.origin);
    coutvec("nouvelle direction ", r.dir);

    return 0;
}