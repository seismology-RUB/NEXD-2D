//Create Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {5, 0, 0, 1.0};
Point(3) = {45, 0, 0, 1.0};
Point(4) = {50, 0, 0, 1.0};
Point(5) = {0, 0, 5, 1.0};
Point(6) = {05, 0, 5, 1.0};
Point(7) = {45, 0, 5, 1.0};
Point(8) = {50, 0, 5, 1.0};
Point(9) = {0, 0, 25, 1.0};
Point(10) = {05, 0, 25, 1.0};
Point(11) = {45, 0, 25, 1.0};
Point(12) = {50, 0, 25, 1.0};
Point(13) = {0, 0, 50, 1.0};
Point(14) = {5, 0, 50, 1.0};
Point(15) = {45, 0, 50, 1.0};
Point(16) = {50, 0, 50, 1.0};
//Create Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {9, 10};
Line(8) = {10, 11};
Line(9) = {11, 12};
Line(10) = {13, 14};
Line(11) = {14, 15};
Line(12) = {15, 16};
Line(13) = {1, 5};
Line(14) = {5, 9};
Line(15) = {9, 13};
Line(16) = {2, 6};
Line(17) = {6, 10};
Line(18) = {10, 14};
Line(19) = {3, 7};
Line(20) = {7, 11};
Line(21) = {11, 15};
Line(22) = {4, 8};
Line(23) = {8, 12};
Line(24) = {12, 16};
//Create Surfaces
Line Loop(1) = {1, 16, -4, -13};
Plane Surface(1) = {1};
Line Loop(2) = {2, 19, -5, -16};
Plane Surface(2) = {2};
Line Loop(3) = {3, 22, -6, -19};
Plane Surface(3) = {3};
Line Loop(4) = {4, 17, -7, -14};
Plane Surface(4) = {4};
Line Loop(5) = {5, 20, -8, -17};
Plane Surface(5) = {5};
Line Loop(6) = {6, 23, -9, -20};
Plane Surface(6) = {6};
Line Loop(7) = {7, 18, -10, -15};
Plane Surface(7) = {7};
Line Loop(8) = {8, 21, -11, -18};
Plane Surface(8) = {8};
Line Loop(9) = {9, 24, -12, -21};
Plane Surface(9) = {9};

//Set material properties
Physical Surface("elastic 1 5800 3800 2600 9999 9999 0") = {5,8};             //normal
Physical Surface("elastic 2 5800 3800 2600 9999 9999 1") = {1, 2, 3, 4, 6, 7, 9}; //pml
//Set free and absorbing boundaries
Physical Line("free") = {10, 11, 12};
Physical Line("absorb") = {13, 14, 15, 1, 2, 3, 22, 23, 24};
