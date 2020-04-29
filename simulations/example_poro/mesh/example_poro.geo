//Create Points
Point(1) = {0, 0, 0, 2.5};
Point(2) = {100, 0, 0, 2.5};
Point(3) = {100, 0, 100, 2.5};
Point(4) = {0, 0, 100, 2.5};
//Create Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//Create Surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
//Set material properties
//poroelastic:     poro nr 1 rhos lambda^u my phi kappa b 1/T 1/N rho1 S1 K1 ny1 rho2 S2 K2 ny2 fitting_n fitting_chi Sr1 Sr2
//              OR poro nr 2 rhos lambda^u my phi kappa b 1/T 1/N rho1 S1 K1 ny1 rho2 S2 K2 ny2 p_b       lambda_BC   Sr1 Sr2
//              OR poro nr 3 rhos K^d      my phi kappa 0 1/T Ks  rho1 S1 K1 ny1 rho2 S2 K2 ny2 fitting_n fitting_chi Sr1 Sr2
//              OR poro nr 4 rhos K^d      my phi kappa 0 1/T Ks  rho1 S1 K1 ny1 rho2 S2 K2 ny2 p_b       lambda_BC   Sr1 Sr2
//              OR poro nr 5 rhos K^d      my phi kappa 0 1/T 1/N rho1 S1 K1 ny1 rho2 S2 K2 ny2 fitting_n fitting_chi Sr1 Sr2
//              OR poro nr 6 rhos K^d      my phi kappa 0 1/T 1/N rho1 S1 K1 ny1 rho2 S2 K2 ny2 p_b       lambda_BC   Sr1 Sr2
//              OR poro nr 7 rhos lambda^u my phi kappa b 1/T 1/N rho1 S1 K1 ny1 rho2 S2 K2 ny2 A         0           Sr1 Sr2
//              OR poro nr 8 rhos K^d      my phi kappa 0 1/T Ks  rho1 S1 K1 ny1 rho2 S2 K2 ny2 A         0           Sr1 Sr2
//              OR poro nr 9 rhos K^d      my phi kappa 0 1/T 1/N rho1 S1 K1 ny1 rho2 S2 K2 ny2 A         0           Sr1 Sr2
Physical Surface("poro 1 1 3000 13.56e9 30.720e9 0.3 1.0e-17 0.7 0.74 1.33e-7 0.099 0.5 20.8e5 0.0 0.055 0.5 4.0e1 0.0 2.86 0.79 0.00 0.00") = {1};
//Set absorbing boundary
Physical Line("absorb") = {1, 2, 3, 4};
