void BasisFctPK(int , vector<vector<long>> &nn, 
        vector<vector<long>> &aa, 
        vector<long> &ff) {
    int idx = 0;
    for (auto &coordinate : coordinate_list) {
        int i = coordinate[0];
        int j = coordinate[1];
        int k = coordinate[2];
        if (i + j + k == p) {
            int ID = 0;
            ff[idx] = factorial(i)*factorial(j)
                    *factorial(k);
            if (i > 0) {
                for (int ii = 0; ii < i ; ii++) {
                    nn[idx][ID] = 0;
                    aa[idx][ID] = ii;
                    ID++;
                }
            }
            // same for j and k
            idx++;
            }
        }
    }