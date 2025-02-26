SetFactory("OpenCASCADE");
  
Mesh.Algorithm = 6;
Mesh.CharacteristicLengthMin = 00.01;
Mesh.CharacteristicLengthMax = 20.00;
Mesh.PartitionOldStyleMsh2 = 1;
Mesh.PreserveNumberingMsh2 = 1;

//
// Author: Sohail Rathore
// Date  : 16/01/2025
//
// This is a parameterised
// model of the Oxford-Nano-
// systems setup for a single
// sample
//
//

//Tank dimensions
L_tank=1000.0;
W_tank=500.0;
D_tank=715.0;

//Anode dimensions
L_anode=840.0;
W_anode=710.0;
th_anode=5.0;

//Sample dimensions
L_sampl=700.0;
W_sampl=700.0;
th_sampl=0.50;

//Anode Spacing
S=150.0;

//Corner positions of
//anodes and samples
cx_1 =(L_tank-L_anode)/2.0;
cy_1 =(W_tank/2.0) - (th_anode/2.0) - S;
cz_1 =(D_tank-W_anode)/2.0;

cx_2 =(L_tank-L_anode)/2.0;
cy_2 =(W_tank/2.0) - (th_anode/2.0) + S;
cz_2 =(D_tank-W_anode)/2.0;

cx_3 =(L_tank - L_sampl)/2.0;
cy_3 =(W_tank - th_sampl)/2.0;
cz_3 =(D_tank - W_sampl)/2.0;

//The Plain geometry without Shielding
Box(1) = {0, 0, 0, L_tank, W_tank, D_tank};                                         //Liquid Tank
Box(2) = {cx_1, cy_1, cz_1, L_anode, th_anode, W_anode}; //Anode 1
Box(3) = {cx_2, cy_2, cz_2, L_anode, th_anode, W_anode}; //Anode 2
Box(4) = {cx_3, cy_3, cz_3, L_sampl, th_sampl, W_sampl}; //Sample

//Shielding dimensions
H_Shield=50.0;
DiaT_Shield=5.5;

L_Shield1 = L_sampl*(1.001); //Adding 0.001% to remove edge problems
th_Shield1= th_sampl+(2.0*H_Shield);
W_Shield1 = W_sampl*(1.001); //Adding 0.001% to remove edge problems

L_Shield2 = L_sampl-(2.0*DiaT_Shield);
th_Shield2= th_sampl+(2.2*H_Shield);
W_Shield2 = W_sampl-(2.0*DiaT_Shield);

cx_4 =(L_tank - L_Shield1)/2.0;
cy_4 =(W_tank - th_Shield1)/2.0;
cz_4 =(D_tank - W_Shield1)/2.0;

cx_5 =(L_tank - L_Shield2)/2.0;
cy_5 =(W_tank - th_Shield2)/2.0;
cz_5 =(D_tank - W_Shield2)/2.0;

Box(5) = {cx_4, cy_4, cz_4, L_Shield1, th_Shield1, W_Shield1}; //Sample Shielding Outside
Box(6) = {cx_5, cy_5, cz_5, L_Shield2, th_Shield2, W_Shield2}; //Sample Shielding Inside

BooleanDifference(7) = { Volume{5}; Delete; }{ Volume{6}; Delete; }; //Shielding
BooleanDifference(8) = { Volume{1}; Delete; }{ Volume{2,3,4,7}; Delete; };

Physical Volume(1)  = {8};
Physical Surface("sample", 1) = {13, 28};
Physical Surface("Anode1", 2) = { 7,  8,  9, 10, 11, 12};
Physical Surface("Anode2", 3) = {29, 30, 31, 32, 33, 34};
Physical Surface("Shield", 4) = {14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};
Physical Surface("TAndB",  5) = {3,5};
Physical Surface("inOut1", 6) = {1};
Physical Surface("inOut2", 7) = {6};
Physical Surface("inOut3", 8) = {2};
Physical Surface("inOut4", 9) = {4};
Mesh.MshFileVersion = 2.2;
