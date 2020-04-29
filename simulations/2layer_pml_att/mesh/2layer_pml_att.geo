//Add points
Point(1) = {0, 0, 0, 1};
Point(2) = {0, 0, 50, 1};
Point(3) = {50, 0, 50, 1};
Point(4) = {50, 0, 0, 1};
Point(5) = {0, 0, 25, 1};
Point(6) = {50, 0, 25, 1};
Point(7) = {5, 0, 0, 1};
Point(8) = {5, 0, 25, 1};
Point(9) = {5, 0, 50, 1};
Point(10) = {0, 0, 5, 1};
Point(11) = {5, 0, 5, 1};
Point(12) = {50, 0, 5, 1};
Point(13) = {45, 0, 0, 1};
Point(14) = {45, 0, 5, 1};
Point(15) = {45, 0, 25, 1};
Point(16) = {45, 0, 50, 1};
//Add Lines
Line(1) = {1, 7};
Line(2) = {7, 13};
Line(3) = {13, 4};
Line(4) = {4, 12};
Line(5) = {12, 6};
Line(6) = {6, 3};
Line(7) = {3, 16};
Line(8) = {16, 9};
Line(9) = {9, 2};
Line(10) = {2, 5};
Line(11) = {5, 10};
Line(12) = {10, 1};
Line(13) = {7, 11};
Line(14) = {11, 10};
Line(15) = {13, 14};
Line(16) = {14, 11};
Line(17) = {12, 14};
Line(18) = {11, 8};
Line(19) = {8, 5};
Line(20) = {14, 15};
Line(21) = {15, 8};
Line(22) = {6, 15};
Line(23) = {8, 9};
Line(24) = {15, 16};
//Create surfaces
Line Loop(1) = {1, 13, 14, 12};
Plane Surface(1) = {1};
Line Loop(2) = {2, 15, 16, -13};
Plane Surface(2) = {2};
Line Loop(3) = {3, 4, 17, -15};
Plane Surface(3) = {3};
Line Loop(4) = {-14, 18, 19, 11};
Plane Surface(4) = {4};
Line Loop(5) = {-16, 20, 21, -18};
Plane Surface(5) = {5};
Line Loop(6) = {-17, 5, 22, -20};
Plane Surface(6) = {6};
Line Loop(7) = {-19, 23, 9, 10};
Plane Surface(7) = {7};
Line Loop(8) = {-21, 24, 8, -23};
Plane Surface(8) = {8};
Line Loop(9) = {-22, 6, 7, -24};
Plane Surface(9) = {9};
//Set material properties
Physical Surface("elastic 1 3000 2200 2000 10 10 0") = {8};
Physical Surface("elastic 2 3000 2200 2000 10 10 1") = {7, 9};
Physical Surface("elastic 3 5000 3200 3000 9999 9999 0") = {5};
Physical Surface("elastic 4 5000 3200 3000 9999 9999 1") = {4, 1, 2, 3, 6};
//Define free and absorbing boundaries
Physical Line("free") = {9, 8, 7};
Physical Line("absorb") = {10, 11, 12, 1, 2, 3, 4, 5, 6};