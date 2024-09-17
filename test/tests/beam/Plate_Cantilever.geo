// Gmsh project created on Mon Mar 25 12:56:44 2024
SetFactory("OpenCASCADE");
//+
//+
//+
Point(1) = {0, 0, 0, 0.75};
//+
Point(2) = {10, 0, 0, 0.75};
//+
Point(3) = {10, 0, 1, 0.75};
//+
Point(4) = {0, 0, 1, 0.75};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Surface(1) = {1};
//+
Physical Curve("right", 5) = {2};
//+
Physical Curve("left", 6) = {4};
//+
Physical Curve("top", 7) = {3};
//+
Physical Curve("bottom", 8) = {1};
//+
Physical Surface("shell", 9) = {1};
//+
Transfinite Surface {1} Alternated;
//+
Transfinite Curve {4, 2} = 4 Using Progression 1;
