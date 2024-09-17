// Gmsh project created on Mon Mar 25 12:56:44 2024
SetFactory("OpenCASCADE");
//+
//+
//+
Point(1) = {0, 0, 0, 5};
//+
Point(2) = {300, 0, 0, 5};
//+
Point(3) = {0, 300, 0, 5};

//+
Circle(1) = {2, 1, 3};
//+
Extrude {0, 0, 300} {
  Curve{1}; 
}
//+
Physical Curve("left", 5) = {3};
//+
Physical Curve("right", 6) = {2};
//+
Physical Curve("back", 7) = {4};
//+
Physical Curve("front", 8) = {1};
//+
Physical Surface("shell", 9) = {1};
