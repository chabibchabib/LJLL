//#include "ff++.hpp"
//#include "AddNewFE.h"
#include"utils.hpp"
#include"polynomial.hpp"

using namespace std; 
// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend
// de l'orientation des aretes
//
/// ---------------------------------------------------------------
namespace Fem2D {
// ------ P3  Hierarchical (just remove P1 node of the P2 finite element)  --------
class TypeOfFE_PkLagrange : public TypeOfFE {

    private:
    static vector<int> PrepareData(int PK) {
        cout<<"\nSTART\n";
        int ndof = (PK + 2) * (PK + 1) / 2;
        vector<int> data(5*ndof+3,0);
        vector<double> dummy_pi_coef;
        // Remplir data ici
        FillDataLagrange(PK, data, dummy_pi_coef);
        return data;
    }
    public:
    int k ;
    int ndf;
     vector<int> Data;
     vector<double> Pi_h_coef;
     vector<vector<int>> nn;//[ndf][k];
     vector<vector<int>>aa;//[ndf][k];
     vector<int> ff;//[ndf];
     vector<int> il;//[ndf];
     vector<int> jl;//[ndf];
     vector<int> kl;//[ndf];
    TypeOfFE_PkLagrange( int PK ) : TypeOfFE((PK+1)*(PK+2)/2, 1, PrepareData(PK).data(), PK, 1, (PK+1)*(PK+2)/2 +6, (PK+1)*(PK+2)/2, 0) {
        k=PK;
        ndf = (k + 2) * (k + 1) / 2;
        nn.resize(ndf);
        aa.resize(ndf);
        for (int i=0;i<ndf;i++){
            aa[i].resize(k,-1);
            nn[i].resize(k,-1);
        }
        ff.resize(ndf,-1);
        il.resize(ndf,-1);
        jl.resize(ndf,-1);
        kl.resize(ndf,-1);
        // Fill Data Lagrange
        FillDataLagrange( k, Data, Pi_h_coef );

        //Fill ff,il,jl,kl,nn,aa
        BasisFctPK(k, nn,aa, ff, il,jl,kl);
       vector<R2> Pt=PtConstruction( k);//[ndf] //= {R2(0 / 3., 0 / 3.), R2(3 / 3., 0 / 3.), R2(0 / 3., 3 / 3.),
            //R2(2 / 3., 1 / 3.), R2(1 / 3., 2 / 3.), R2(0 / 3., 2 / 3.),
            //R2(0 / 3., 1 / 3.), R2(1 / 3., 0 / 3.), R2(2 / 3., 0 / 3.),
            //R2(1 / 3., 1 / 3.)}; *
        cout<<endl;
        for (auto elm : il) cout<<elm<<" ";
        cout<<endl;
        for (auto elm : jl) cout<<elm<<" ";
        cout<<endl;        
        for (auto elm : kl) cout<<elm<<" ";
        cout<<endl;
        for (auto elm : nn) {for(auto elmelm : elm) cout<<elmelm<<" ";}
        cout<<endl;
        for (auto elm : aa) {for(auto elmelm : elm) cout<<elmelm<<" ";}
        cout<<endl;
        // 3,4,5,6,7,8
        vector<int> other(ndf,0);//[ndf] = {-1, -1, -1, 4, 3, 6, 5, 8, 7, -1}; // 
        FillOther(other, k,  ndf);
        int kk = 0;

        for (int i = 0; i < NbDoF; i++) {
            pij_alpha[kk++] = IPJ(i, i, 0);
            if (other[i] != i) {
                pij_alpha[kk++] = IPJ(i, other[i], 0);
            }
            
            P_Pi_h[i] = Pt[i];
        }        
        assert(P_Pi_h.N( ) == NbDoF);
        assert(pij_alpha.N( ) == kk);

    }
    
    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
        for (int i = 0; i < ndf+6; ++i) {
            v[i] = 1;
        }
        
        int e0 = K.EdgeOrientation(0);
        int e1 = K.EdgeOrientation(1);
        int e2 = K.EdgeOrientation(2);
        int ooo[6] = {e0, e0, e1, e1, e2, e2};
        int iii[6] = {};
        int jjj[6] = {};
        
        for (int i = 0; i < 6; ++i) {
            iii[i] = 3 + 2 * i;    // si  orient = 1
            jjj[i] = 4 + 2 * i;    // si orient = -1
        }
        
        for (int i = 0; i < 6; ++i) {
            if (ooo[i] == 1) {
                v[jjj[i]] = 0;
            } else {
                v[iii[i]] = 0;
            }
        }
    }
};

// on what     nu df on node node of df
/*int TypeOfFE_PkLagrange::Data[] = {
    0, 1, 2, 3, 3, 4, 4, 5, 5, 6,    // the support number  of the node of the df
    0, 0, 0, 0, 1, 0, 1, 0, 1, 0,    // the number of the df on  the node
    0, 1, 2, 3, 3, 4, 4, 5, 5, 6,    // the node of the df
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    // the df come from which FE (generaly 0)
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,    // which are de df on sub FE
    0,                               // for each compontant $j=0,N-1$ it give the sub FE associated
    0, 10};
double TypeOfFE_PkLagrange::Pi_h_coef[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};*/ // T
void TypeOfFE_PkLagrange::FB(const bool *whatd, const Mesh &, const Triangle &K,
                             const RdHat &PHat, RNMK_ &val) const {
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R L[3] = {l0 * k, l1 * k, l2 * k};
    
    throwassert(val.N( ) >= ndf);
    throwassert(val.M( ) == 1);
    // Attention il faut renumeroter les fonction de bases
    // car dans freefem++, il y a un node par sommet, arete or element
    // et la numerotation naturelle  mais 2 noud pas arete
    // donc p est la perumation
    // echange de numerotation si les arete sont dans le mauvais sens
    vector<int> p(ndf,0) ;
    
    for (int i = 0; i < ndf; ++i) {
        p[i] = i;
    }
    
    // Get index exchange 
    vector<pair<int,int>> ExId=Exchangeidx( k);

    if (K.EdgeOrientation(0) < 0) {
        Exchange(p[ExId[0].first], p[ExId[0].second]);   
    }
    
    if (K.EdgeOrientation(1) < 0) {
        Exchange(p[ExId[1].first], p[ExId[1].second]);    
    }
    
    if (K.EdgeOrientation(2) < 0) {
        Exchange(p[ExId[2].first], p[ExId[2].second]);    
    }
    
    val = 0;
    
    if (whatd[op_id]) {
        RN_ f0(val('.', 0, op_id));
        
        for (int df = 0; df < ndf; df++) {
            int pdf = p[df];
            R f = 1. / ff[df];
            
            for (int i = 0; i < k; ++i) {
                f *= L[nn[df][i]] - aa[df][i];
            }
            
            f0[pdf] = f;
        }
    }
    
    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
        R2 D[] = {K.H(0) * k, K.H(1) * k, K.H(2) * k};
        if (whatd[op_dx] || whatd[op_dy]) {
            for (int df = 0; df < ndf; df++) {
                int pdf = p[df];
                R fx = 0., fy = 0., f = 1. / ff[df];
                
                for (int i = 0; i < k; ++i) {
                    int n = nn[df][i];
                    R Ln = L[n] - aa[df][i];
                    fx = fx * Ln + f * D[n].x;
                    fy = fy * Ln + f * D[n].y;
                    f = f * Ln;
                }
                
                if (whatd[op_dx]) {
                    val(pdf, 0, op_dx) = fx;
                }
                
                if (whatd[op_dy]) {
                    val(pdf, 0, op_dy) = fy;
                }
            }
        }
        
        if (whatd[op_dyy] || whatd[op_dxy] || whatd[op_dxx]) {
            for (int df = 0; df < ndf; df++) {
                int pdf = p[df];
                R fx = 0., fy = 0., f = 1. / ff[df];
                R fxx = 0., fyy = 0., fxy = 0.;
                
                for (int i = 0; i < k; ++i) {
                    int n = nn[df][i];
                    R Ln = L[n] - aa[df][i];
                    fxx = fxx * Ln + 2. * fx * D[n].x;
                    fyy = fyy * Ln + 2. * fy * D[n].y;
                    fxy = fxy * Ln + fx * D[n].y + fy * D[n].x;
                    fx = fx * Ln + f * D[n].x;
                    fy = fy * Ln + f * D[n].y;
                    f = f * Ln;
                }
                
                if (whatd[op_dxx]) {
                    val(pdf, 0, op_dxx) = fxx;
                }
                
                if (whatd[op_dyy]) {
                    val(pdf, 0, op_dyy) = fyy;
                }
                
                if (whatd[op_dxy]) {
                    val(pdf, 0, op_dxy) = fxy;
                }
            }
        }
    }
}

static  TypeOfFE_PkLagrange PK3(3);
static  TypeOfFE_PkLagrange PK4(4);

static void init( ) {
/*for (int i=0;i<1000; i++) TablePK[i]=nullptr;

  Global.Add(
      "PK", "(",
      new OneOperator1s_< TypeOfFE_PkLagrange ,long>(GenerateTypeOfFE_PkLagrangeOperator)
  );
  //for (int i=1;i<5;i++)   AddNewFE("P_"+str(i), &PK(i));*/
AddNewFE("PK3", &PK3);
static ListOfTFE FE_P3("PK3", &PK3); // to add P3 in list of Common FE
AddNewFE("PK4", &PK4);
static ListOfTFE FE_P4("PK4", &PK4); // to add P4 in list of Common FE
}

}    // namespace Fem2D
LOADFUNC(Fem2D::init);



/*int main(int argv, char **argc){
    int PK= atoi(argc[1]);
    int Ndof=(PK+1)*(PK+2)*0.5;
    vector<R2> Test= PtConstruction(PK);
    for(auto &elem :Test ) cout<<"elem:"<<elem.i<<","<<elem.j<<endl;
    vector<int> Data;
    vector<int> Pi_dataLagrange;

    FillDataLagrange(PK, Data,Pi_dataLagrange );
    //cout<<"Row 1:"<<endl;
    //for(int i=0; i<Ndof;i++ ) cout<<Data[i]<<", ";
    //cout<<"\nRow 2:"<<endl;
    //for(int i=0; i<Ndof;i++ ) cout<<Data[Ndof+i]<<", ";
    //cout<<"\nRow 3:"<<endl;
    //for(int i=0; i<Ndof;i++ ) cout<<Data[2*Ndof+i]<<", ";
    //cout<<"\nRow 4:"<<endl;
    //for(int i=0; i<Ndof;i++ ) cout<<Data[3*Ndof+i]<<", ";
    //cout<<"\nRow 5:"<<endl;
    //for(int i=0; i<Ndof;i++ ) cout<<Data[4*Ndof+i]<<", ";
    //cout<<"\nRow 6:"<<endl;
    // cout<<Data[5*Ndof];
    //cout<<"\nRow 7:"<<endl;
    // cout<<Data[5*Ndof+1]<<" "<<Data[5*Ndof+2]<<endl;
    //for (auto elem :Pi_dataLagrange ) cout<<" "<<elem<<endl;
    //vector<pair<int,int>> Idxs=Exchangeidx( PK);
    //for (auto elm : Idxs ) cout<<"("<<elm.first<<","<<elm.second<<")";
    vector<int> other(Ndof,0);
    for(int i=0; i<Ndof;i++ ) other[i]=i;
    for(int i=0; i<3;i++ ) {
        for (int j=0; j<((PK-1)/2);j++ ){
            other[3+i*(PK-1)+j]=3+(i+1)*(PK-1)-j-1;
            other[3+(i+1)*(PK-1)-j-1]=3+i*(PK-1)+j;
        }
    }
     for(int i=0; i<Ndof;i++ ) cout<<other[i]<<" ";
    cout<<endl;
    return 0; 
}*/