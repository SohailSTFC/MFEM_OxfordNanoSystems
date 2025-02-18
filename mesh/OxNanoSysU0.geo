SetFactory("OpenCASCADE");
  
Mesh.Algorithm = 6;
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.2;
Mesh.PartitionOldStyleMsh2 = 1;
Mesh.PreserveNumberingMsh2 = 1;

L=0.5;
L_anode=0.1;
L_sampl=0.05;
th_anode=0.005;
th_sampl=0.0005;
S=0.1;

Rectangle(1) = {0, 0, 0, L, L};                                         //Liquid Tank
Rectangle(2) = {L/2-L_anode/2, L/2-th_anode/2+S, 0, L_anode, th_anode}; //Anode 1
Rectangle(3) = {L/2-L_anode/2, L/2-th_anode/2-S, 0, L_anode, th_anode}; //Anode 2
Rectangle(4) = {L/2-L_sampl/2, L/2+th_sampl/2,   0, L_sampl, th_sampl}; //Sample

BooleanDifference(5) = { Surface{1}; Delete; }{ Surface{2,3,4}; Delete; };


Physical Surface(1) = {5};
Physical Curve(4) = {1, 2, 3, 4};
Physical Curve(3) = {13, 14, 15, 16};
Physical Curve(2) = {5, 6, 7, 8};
Physical Curve(1) = {9, 10, 11, 12};

Mesh.MshFileVersion = 2.2;
