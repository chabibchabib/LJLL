  class TypeOfFE_P4Lagrange : public TypeOfFE {
   public:
    static const int k = 4;
    static const int ndf = (k + 2) * (k + 1) / 2;
    static int Data[];
    static double Pi_h_coef[];
    static const int nn[15][4];
    static const int aa[15][4];
    static const int ff[15];
    static const int il[15];
    static const int jl[15];
    static const int kl[15];

    TypeOfFE_P4Lagrange( ) : TypeOfFE(15, 1, Data, 4, 1, 15 + 6, 15, 0) 
    {
      static const R2 Pt[15] = {R2(0 / 4., 0 / 4.), R2(4 / 4., 0 / 4.), 
                                R2(0 / 4., 4 / 4.)...}

      // 3,4,5, 6,7,8, 9,10,11,
      int other[15] = {0, 1, 2, 5, 4, 3, 8, 7, 6, 11, 10, 9, 12, 13, 14};
      ...
    }
}