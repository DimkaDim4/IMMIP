cl = 0.007;
r = 0.08;
h = 0.09;
Point(1) = {0, 0, 0, cl};
Point(2) = {r, 0, 0, cl};
Point(3) = {r, h, 0, cl};
Point(4) = {0, h, 0, cl};
Point(5) = {r/2, h/2, 0, cl};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line(5) = {4, 5};
Line(6) = {5, 3};
Line(7) = {5, 2};
Line(8) = {5, 1};
Line Loop(9) = {7, -2, -6};
Plane Surface(10) = {9};
Line Loop(11) = {6, -1, 5};
Plane Surface(12) = {11};
Line Loop(13) = {4, 5, 8};
Plane Surface(14) = {13};
Line Loop(15) = {8, -3, -7};
Plane Surface(16) = {15};
Physical Surface(17) = {10};
Physical Surface(18) = {12};
Physical Surface(19) = {14};
Physical Surface(20) = {16};
Physical Line(21) = {2};
Physical Line(22) = {1};
Physical Line(23) = {4};
Physical Line(24) = {3};
Recombine Surface "*";
Mesh.SubdivisionAlgorithm = 1;
// Transfinite Surface "*";
// Mesh.Smoothing = 0;
