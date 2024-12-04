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

Recombine Surface {5};
Extrude{0,0,L}{
  Surface{5}; Layers{{8},{1}}; Recombine;
}

Physical Surface("OuterSurface", 4) = {8, 9, 7, 6, 5, 22};
Physical Surface("Anode1", 3) = {21, 20, 18, 19};
Physical Surface("Anode2", 2) = {12, 13, 10, 11};
Physical Surface("Sample", 1) = {14, 16, 17, 15};

Mesh.MshFileVersion = 2.2;
