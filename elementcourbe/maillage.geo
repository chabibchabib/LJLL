lc = 0.1;

// Carré
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};

// Cercle (trou)
Point(5) = {0.5, 0.5, 0, lc};
Point(6) = {0.7, 0.5, 0, lc};
Point(7) = {0.5, 0.7, 0, lc};
Point(8) = {0.3, 0.5, 0, lc};
Point(9) = {0.5, 0.3, 0, lc};

Circle(10) = {6,5,7};
Circle(11) = {7,5,8};
Circle(12) = {8,5,9};
Circle(13) = {9,5,6};

Line Loop(2) = {10,11,12,13};

Plane Surface(1) = {1,2};

//Mesh.ElementOrder = 2;

Mesh 2;