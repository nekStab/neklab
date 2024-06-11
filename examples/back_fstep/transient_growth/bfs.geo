//////////////////////////////////////////////////////////////////////////////////
//////////                                                              //////////
//////////     MESH GENERATION FOR THE ROUNDED BACKWARD FACING-STEP     //////////
//////////                                                              //////////
//////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////
/////     PHYSICAL PARAMETERS OF THE COMPUTATIONAL DOMAIN     /////
///////////////////////////////////////////////////////////////////

h  =   1; // Step height.
Li = -20; // Inflow length.
Lo = 100; // Outflow length.
H  =  20; // Domain height.
Lb =  -2; // Upstream flat plate length.
R  = 5/2; // Radius of curvature of the rounded step.

///////////////////////////////////////////////////////////////
/////     BLOCK-PARTITION OF THE COMPUTATIONAL DOMAIN     /////
///////////////////////////////////////////////////////////////

// Parameters.
H_fs      =  3; // Starting height for the free-stream domain.
L_refined = 20; // Streamwise extent of the refined region to capture
                // the recirculation bubble.
/////
///// Free-stream domain.
/////

// Parameters of the mesh.
scale_fs   = 1.5; // Progression factor in the free stream direction.
scale_out  = 1.05; // Progression factor in the outflow direction.
scale_in   = 0.8; // Progression factor in the inflow part.

nx_inflow     = 12;
nx_flat_plate = 6;
nx_step       = 6;
nx_bubble     = 48;
nx_outflow    = 48;

ynpts_fs   = 9; // Number of elements in the vertical direction.

// List of points defining the different regions
// of the free-stream domain.
Point(1)  = {Li, H_fs, 0};
Point(2)  = {Lb, H_fs, 0};
Point(3)  = {0, H_fs, 0};
Point(4)  = {2*h, H_fs, 0};
Point(5)  = {L_refined, H_fs, 0};
Point(6)  = {Lo, H_fs, 0};
Point(7)  = {Lo, H, 0};
Point(8)  = {L_refined, H, 0};
Point(9)  = {2*h, H, 0};
Point(10) = {0, H, 0};
Point(11) = {Lb, H, 0};
Point(12) = {Li, H, 0};

// Upstream free-stream domain.
Line(1) = {1, 2};
Line(2) = {2, 11};
Line(3) = {11, 12};
Line(4) = {12, 1};

Transfinite Line {1, -3}  = nx_inflow Using Progression scale_in;
Transfinite Line {2, -4} = ynpts_fs Using Progression scale_fs;

Curve Loop(1) = {1, 2, 3, 4};
Surface(1) = {1}; Transfinite Surface{1}; Recombine Surface{1};

// Flat-plate free stream domain region.
Line(5) = {2, 3};
Line(6) = {3, 10};
Line(7) = {10, 11};

Transfinite Line {5, 7} = nx_flat_plate;
Transfinite Line {2, 6} = ynpts_fs Using Progression scale_fs;

Curve Loop(2) = {5, 6, 7, -2};
Surface(2) = {2} ; Transfinite Surface{2}; Recombine Surface{2};


// Step region free stream domain.
Line(8)  = {3, 4};
Line(9)  = {4, 9};
Line(10) = {9, 10};

Transfinite Line {8, 10} = nx_step;
Transfinite Line {9, 6} = ynpts_fs Using Progression scale_fs;

Curve Loop(3) = {8, 9, 10, -6};
Surface(3) = {3}; Transfinite Surface{3}; Recombine Surface{3};

// Reversed-flow region free-stream domain.
Line(11) = {4, 5};
Line(12) = {5, 8};
Line(13) = {8, 9};

Transfinite Line {11, 13} = nx_bubble;
Transfinite Line {12, 9} = ynpts_fs Using Progression scale_fs;

Curve Loop(4) = {11, 12, 13, -9};
Surface(4) = {4}; Transfinite Surface{4}; Recombine Surface{4};

// Output free-stream domain.
Line(14) = {5, 6};
Line(15) = {6, 7};
Line(16) = {7, 8};

Transfinite Line {14, -16} = nx_outflow Using Progression scale_out;
Transfinite Line {15, 12} = ynpts_fs Using Progression scale_fs;

Curve Loop(5) = {14, 15, 16, -12};
Surface(5) = {5}; Transfinite Surface{5}; Recombine Surface{5};

/////
/////     Boundary layer domain.
/////

// Parameters of the mesh.
scale_bl   = 1; // Progession factor.
ynpts_bl   = 17;    // Number of elements in the vertical direction.

// List of additional points defining the different
// regions of the boundary layer domain.
Point(13) = {Li, 1, 0};
Point(14) = {Lb, 1, 0};
Point(15) = {0, 1, 0};
Point(16) = {2*h, 0, 0};
Point(17) = {L_refined, 0, 0};
Point(18) = {Lo, 0, 0};

// Upstream boundary layer domain.
Line(17) = {13, 14};
Line(18) = {14, 2};
Line(19) = {1, 13};

Transfinite Line {17, 1} = nx_inflow Using Progression scale_in;
Transfinite Line {18, 19} = ynpts_bl Using Progression scale_bl;

Curve Loop(6) = {17, 18, -1, 19};
Surface(6) = {6}; Transfinite Surface{6}; Recombine Surface{6};

// Flat plate region
Line(20) = {14, 15};
Line(21) = {15, 3};

Transfinite Line {20, 5}  = nx_flat_plate;
Transfinite Line {21, 18} = ynpts_bl Using Progression scale_bl;

Curve Loop(7) = {20, 21, -5, -18};
Surface(7) = {7}; Transfinite Surface{7}; Recombine Surface{7};

// Step region domain.
Point(19) = {0, 1-R, 0}; // Center of the cylinder defining the step.
Circle(22) = {15, 19, 16};
Line(23) = {16, 4};

Transfinite Line {22, 8} = nx_step;
Transfinite Line {21, 23} = ynpts_bl Using Progression scale_bl;

Curve Loop(8) = {22, 23, -8, -21};
Surface(8) = {8}; Transfinite Surface{8}; Recombine Surface{8};

// Reversed flow region.
Line(24) = {16, 17};
Line(25) = {17, 5};

Transfinite Line {24, 11} = nx_bubble;
Transfinite Line {23, 25} = ynpts_bl Using Progression scale_bl;

Curve Loop(9) = {24, 25, -11, -23};
Surface(9) = {9}; Transfinite Surface{9}; Recombine Surface{9};

// Downstream boundary layer domain.
Line(26) = {17, 18};
Line(27) = {18, 6};

Transfinite Line {26, 14} = nx_outflow Using Progression scale_out;
Transfinite Line {25, 27} = ynpts_bl Using Progression scale_bl;

Curve Loop(10) = {26, 27, -14, -25};
Surface(10) = {10}; Transfinite Surface{10}; Recombine Surface{10};

Coherence;

///////////////////////////////////////
/////     BOUNDARY CONDITIONS     /////
///////////////////////////////////////

// Volume.
Physical Surface("Fluid", 1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

// Inflow.
Physical Curve("Inflow", 2) = {4, 19};

// Outflow.
Physical Curve("Outflow", 3) = {15, 27};

// Top.
Physical Curve("Sym", 4) = {3, 7, 10, 13, 16, 17};

// Wall.
Physical Curve("Wall", 5) = {20, 22, 24, 26};

////////////////////////////////////////
/////     MESHING THE GEOMETRY     /////
////////////////////////////////////////

// Meshing Process.
Mesh 1;
Mesh 2;
Mesh 3;
SetOrder = 2;
RenumberMeshElements;

// Export mesh.
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 0;
Mesh.Binary = 0;

Save "bfs.msh";

