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
L_Shield=0.03;
th_Shield=0.05;

Box(1) = {0, 0, 0, L, L, L};                                         //Liquid Tank
Box(2) = {L/2-L_anode/2, L/2-th_anode/2+S, 0, L_anode, th_anode, L}; //Anode 1
Box(3) = {L/2-L_anode/2, L/2-th_anode/2-S, 0, L_anode, th_anode, L}; //Anode 2
Box(4) = {L/2-L_sampl/2, L/2+th_sampl/2,   0, L_sampl, th_sampl, L}; //Sample

BooleanDifference(5) = { Volume{1}; Delete; }{ Volume{2,3,4}; Delete; };

Physical Volume(1)  = {5};
Physical Surface(1) = {11, 14, 12, 13};
Physical Surface(2) = {15, 18, 16, 17};
Physical Surface(3) = {7, 10, 8, 9};
Physical Surface(4) = {6, 2, 1, 4};
Physical Surface(5) = {3, 5};

Mesh.MshFileVersion = 2.2;