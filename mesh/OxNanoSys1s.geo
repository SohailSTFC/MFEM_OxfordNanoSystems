SetFactory("OpenCASCADE");
  
Mesh.Algorithm = 6;
Mesh.CharacteristicLengthMin = 0.0005;
Mesh.CharacteristicLengthMax = 0.4;
Mesh.PartitionOldStyleMsh2 = 1;
Mesh.PreserveNumberingMsh2 = 1;

//Standard workipiece dimensions
ExtrL=0.1;
L=0.5;
L_anode=0.1;
L_sampl=0.05;
th_anode=0.005;
th_sampl=0.0005;
S=0.1;

//Shielding dimensions
H_Shield=2*0.03;
th_Shield=2*0.005;
T=th_Shield;

Rectangle(1) = {0, 0, 0, L, L};                                         //Liquid Tank
Rectangle(2) = {L/2-L_anode/2, L/2-th_anode/2+S, 0, L_anode, th_anode}; //Anode 1
Rectangle(3) = {L/2-L_anode/2, L/2-th_anode/2-S, 0, L_anode, th_anode}; //Anode 2
Rectangle(4) = {L/2-L_sampl/2, L/2+th_sampl/2,   0, L_sampl, th_sampl}; //Sample

//Shielding Geometry
Rectangle(5) = {L/2-L_sampl/2-T, L/2+th_sampl/2-H_Shield/2, 0, th_Shield, H_Shield}; //Shield Left
Rectangle(6) = {L/2+L_sampl/2,   L/2+th_sampl/2-H_Shield/2, 0, th_Shield, H_Shield}; //Shield Right

BooleanDifference(7) = { Surface{1}; Delete; }{ Surface{2,3,4,5,6}; Delete; };
Recombine Surface {7};

Extrude{0,0,ExtrL}{
  Surface{7}; Layers{{2},{1}}; Recombine;
}

Physical Volume(1)  = {1};
Physical Surface(10)= {23, 24, 25, 26, 27};//Shielding 1
Physical Surface(9) = {17, 18, 19, 20, 21};//Shielding 2
Physical Surface(8) = {7, 32};             //Top and bottom
Physical Surface(7) = {9};                 //Back
Physical Surface(6) = {8};                 //Left
Physical Surface(5) = {10};                //Right
Physical Surface(4) = {11};                //Front
Physical Surface(3) = {12, 13, 14, 15};    //Anode 1
Physical Surface(2) = {28, 29, 30, 31};    //Anode 2
Physical Surface(1) = {16, 22};            //Sample
Mesh.MshFileVersion = 2.2;
