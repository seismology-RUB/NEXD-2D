//Create Points
Point(1) = {0, 0, 0, 35};
Point(2) = {4000, 0, 0, 35};
Point(3) = {4000, 0, 2705.31, 35};
Point(4) = {0, 0, 2000, 35};
//Create Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//Create Surface
Line Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};
//Set Material Properties
Physical Surface("elastic 1 3200 1847.5 2200 9999 9999 0") = {1};
//Set free and absorbing boundaries
Physical Line("free") = {3};
Physical Line("absorb") = {4, 1, 2};
