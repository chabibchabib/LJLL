#include<iostream>
#include<vector>
#include<cassert>
using namespace std;
template<int dim>
class polygon{
    private:
        vector<vector<double>> vertices; // Number of vertices
    public:
        polygon(int n){
            vertices.resize(n,vector<double>(dim,0.));
        };
        void affichage();
        void setter(vector<vector<double>> vectorofv);
        bool intersection(vector<double> point);

};

template<int dim>
void polygon<dim>::affichage(){
    cout<<"Vertices:\n";
    for (auto elms : vertices) {
        for (auto elm : elms){
        cout<<elm<<" ";
        }
            cout<<endl;
    }
}

template<int dim>
void polygon<dim>::setter(vector<vector<double>> vectorofv){
    vertices.assign(vectorofv.begin(),vectorofv.end());
}

template<int dim>
bool polygon<dim>::intersection(vector<double> point){
    assert(point.size()==dim);
    int num_intersections = 0;
    int num_vertices=vertices.size();
    for (int i=0; i<num_vertices;i++){
        vector<double> v1= vertices[i];
        vector<double> v2= vertices[(i+1)%num_vertices];
        if(((v1[1] > point[1]) != (v2[1] > point[1])) && (point[0] < (v2[0] - v1[0]) * (point[1] - v1[1]) / (v2[1] - v1[1]) + v1[0])){
          num_intersections += 1;  
        }

    }

    return num_intersections % 2 == 1;
}
 int main(){
    polygon<2> triangle(3);
    triangle.affichage();
    vector<vector<double> >v={{0,0},{1,0},{0,1}};
    triangle.setter(v);
    triangle.affichage();
    vector<double> point={0.5,0.5};
    cout<<"Point ("<<point[0]<<"," << point[1]<<") appartient? "<<triangle.intersection(point)<<endl;
    return 0;
 }