//
//  main.cpp
//  Pk-Rd-data
//
//  Created by Frédéric Hecht on 14/10/2025.
//

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include<vector>
using namespace std;

typedef int int4 [4] ;

// InvIntFunc

inline int NumSimplex2(int i) {return ((i)*(i+1))/2;}
inline int NumSimplex2_1(int l) { return sqrt(l*2)+3;}
inline int NumSimplex3(int i) {return ((i)*(i+1)*(i+2))/6;}
inline int NumSimplex3_1(int l) { return pow(l*6,1./3)+4;}

#define InvIntFunction  invNumSimplex2
#define F(i) NumSimplex2(i)
#define F_1(i) NumSimplex2_1(i)
#include "InvIntFunc.hpp"
#undef F
#undef F_1
#undef InvIntFunction

#define InvIntFunction invNumSimplex3
#define F(i) NumSimplex3(i)
#define F_1(i) NumSimplex3_1(i)
#include "InvIntFunc.hpp"
#undef F
#undef F_1
#undef InvIntFunction


// nb ndof dim 1, 2 3 pour PK
inline int ndofPk1(int N){return N+1;}
inline int ndofPk2(int N){return (N+1)*(N+2)/2;}
inline int ndofPk3(int N){return (N+1)*(N+2)*(N+3)/6;}
inline int ndofPk(int d, int N) {
    if(d==1) return N+1;
    else if(d==2) return (N+1)*(N+2)/2;
    else if (d==3)  return (N+1)*(N+2)*(N+3)/6;
    else assert(0);
    return 0;
}

// nomerotation des vertex(dof) dans le simplex
inline int NumSimplex1(int i) { return i;}
inline int NumSimplex2(int i,int j) { return j+NumSimplex2(i+j);}
inline int NumSimplex3(int i,int j,int k) { return NumSimplex3(i+j+k)+NumSimplex2(j+k)+k;}


 void  invNumSimplex2(int n,int &i,int &j)
{
  int k= invNumSimplex2(n); //( i+j)
  j=n-NumSimplex2(k);
  i= k-j;
  //  cout << n << " " << k << " -> " << i << " " << j << endl;
  assert( n == NumSimplex2(i,j));
}
 void  invNumSimplex3(int n,int &i,int &j,int &k)
{
  int l= invNumSimplex3(n); //=( i+j+k)
  invNumSimplex2(n-NumSimplex3(l),j,k);
  assert(j>=0 && k>=0);
  i=l-k-j;
  // cout << n << "   " << l << "-> " << i << " " << j << " " << k <<endl;
  assert( n == NumSimplex3(i,j,k)) ;
}
inline void  invNumSimplex(int d,int n,int i[])
{
    i[0]=i[1]=i[2]=0;
    if (d==1) i[0]=n;
    else if (d==2) invNumSimplex2( n, i[0],i[1]);
    else if (d==3) invNumSimplex3( n, i[0],i[1],i[2]);
    else assert(0); // bug !!!

}


class DataPk {
public:
    int d;
    int N; // degre
    int n; // numero Fb.
    vector<int> nn,aa;
    int ff;
    DataPk(int dd,int NN,int nnn) :
     d(dd),N(NN),n(nnn),nn(N),aa(N)
    {
        int N1 = N+1;
        int nbl =(d+1)*N1;//  (d+1)*N1 faiseau de droite
        int i[4];
        invNumSimplex(d,n,i);
        assert(d<=3);
        i[d]= N-i[0]-i[1]-i[2];// ok car i est fill by 0
        cout << n << " : " <<i[0] << " " << i[1] << " " << i[2] << " " <<i[3] << endl;
        
        // loop on les lignes  du faiseau
        
        int kbl =0;
        ff=1;
        for(int l=0; l<nbl;++l)
        {
            int kk = l/N1; // L_kk
            int ii = l%N1; // L_kk*N - ii
            assert(kk<d+1 && ii< N1);
           // cout << kk  << " " << ii << " " << i[ii] << endl;
            if(  ii < i[kk]) {
                
                ff*= (N-ii);
                nn[kbl]= kk;
                aa[kbl]= ii;
                 cout << kk  << " " << ii << " / " << i[ii] << endl;
                kbl++;
            }
        }
        

    }
}
;
// Overloading << operator 
ostream & operator<<(ostream & f, DataPk & Pk)
{
    f << Pk.d << " " << Pk.N << " " << Pk.n ;
    int i[4];
    invNumSimplex(Pk.d,Pk.n,i);
    i[Pk.d]= Pk.N-i[0]-i[1]-i[2];// ok car i est fill by 0
    for(int j=0; j < Pk.d+1;++j)
        f << " "<< i[j] ;
    f << " -- ";
    for (int l=0; l<Pk.N;++l)
        f << Pk.nn[l] << " " << Pk.aa[l] << "; ";
    f << Pk.ff << endl;
    return f;
}

int data(int d, int N, int* dd)
{
    // par point ii
   //  N+1 data : nn , aa , ff = (L_nn - aa )
 
    int N1 = N+1;
   
    int nbl =(d+1)*N1;//  (d+1)*N1 faiseau de droite
    int nv = ndofPk(d,N);
    int id=0;
    for (int ii=0; ii< nv;++ii) //
    {
        int i[4];
        invNumSimplex(d,ii,i);
        assert(d<=3);
        i[d]= N-i[0]-i[1]-i[2];// ok car i est fill by 0
        cout << ii << " : " <<i[0] << " " << i[1] << " " << i[2] << " " <<i[3] << endl;
        
        // loop on les lignes  du faiseau
        
        int kbl =0;
        int ff=1;
        int id0=id;
        for(int l=0; l<nbl;++l)
        {
            int kk = l/N1; // L_kk
            int ii = l%N1; // L_kk*N - ii
            assert(kk<d+1 && ii< N1);
           // cout << kk  << " " << ii << " " << i[ii] << endl;
            if(  ii < i[kk]) {
                kbl++;
                ff*= (N-ii);
                dd[id++]= kk;
                dd[id++]= ii;
                 cout << kk  << " " << ii << " / " << i[ii] << endl;
            }
        }
        

        dd[id++]=ff;
        for (int l=0; l<N;++l)
            cout << " " << l << "  : " << dd[id0+2*l]<< " "<< dd[id0+2*l+1] << endl;
        cout << " ff "<< ff << endl;
        assert(kbl == N);
    }
    return nv;
}

int main(int argc,const char **argv)
{
    int d =2, N= 2;
    int dd[100000];
    if(argc>2)
    {
         d = atoi(argv[1]);
         N = atoi(argv[2]);
    }
    std::cout << " d = " << d << " Pk :" << N << endl;
    int nPk =ndofPk(d,N);
    vector<DataPk> Pk;
    for(int i=0;i<nPk;++i)
        Pk.push_back(DataPk(d,N,i));
    cout<<"Results main "<<endl;
    for(int i=0;i<nPk;++i)
        std::cout <<i<< " " <<  Pk[i] << endl;
                 
    return 0;
}
