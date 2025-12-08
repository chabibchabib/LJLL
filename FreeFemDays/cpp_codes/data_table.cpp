  int TypeOfFE_P4Lagrange::Data[] = {
    // the support number  of the node of the dof
    0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5,  5,  6,  6,  6,
    // the number of the dof on the support
    0, 0, 0, 0, 1, 2, 0, 1, 2, 0, 1,  2,  0,  1,  2,
    // the node of the dof
    0, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5,  5,  6,  6,  6,
    // the dof come from which FE (generaly 0)
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,    
    // which are de dof on sub FE
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
    // for each compontant $j=0,N-1$ it give the sub FE associated    
    0,
    // First dof
    0, 
    // #dof
    15 
    };