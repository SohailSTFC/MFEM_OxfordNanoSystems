#ifndef DARCYVISUALISATION_HPP
#define DARCYVISUALISATION_HPP 

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

void ParaViewVisualise(std::vector<ParGridFunction*> Fields
                     , std::vector<std::string>      FieldNames
					 , int order, ParMesh *pmesh, double time){
  
  ParaViewDataCollection paraview_dc("Example5P", pmesh);
  paraview_dc.SetPrefixPath("ParaView");
  paraview_dc.SetLevelsOfDetail(order);
  paraview_dc.SetDataFormat(VTKFormat::BINARY);
  paraview_dc.SetHighOrderOutput(true);
  paraview_dc.SetCycle(0);
  paraview_dc.SetTime(time);
  for(int I=0; I< FieldNames.size(); I++) paraview_dc.RegisterField(FieldNames[I],Fields[I]);
  paraview_dc.Save();
};

#endif