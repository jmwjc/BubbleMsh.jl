
n = 4;
Point(1) = {0.0, 0.0, 0.0};
Point(2) = {1.0, 0.0, 0.0};
Point(3) = {1.0, 1.0, 0.0};
Point(4) = {0.0, 1.0, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Transfinite Curve{1,2,3,4} = n+1;
Physical Curve("Γ") = {1,2,3,4};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh 2;
