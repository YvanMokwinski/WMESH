#pragma once

#include "MeshGeometry.hpp"
#include "MeshTopology.hpp"
#include "MeshTopology2D.hpp"
#include "MeshTopology3D.hpp"
#include "MeshProperty.hpp"

template <unsigned int _dimension,typename _real_t,DimensionType::enum_t _topologicalDimension> class Mesh
{
public: MeshTopology<_topologicalDimension> * m_topology;
private: MeshGeometry<_dimension,_real_t> * m_geometry;
  
public: Mesh(const char * filename_)
  {    
    Input::Medit inputMedit(filename_);    
    const auto numNodes = inputMedit.GetNumNodes();
    this->m_geometry = new MeshGeometry<_dimension,_real_t>(numNodes);    
    inputMedit.ReadGeometry(this->m_geometry->GetPoints());
    this->m_topology = new MeshTopology<_topologicalDimension>(inputMedit);
    MeshProperty<int,1,DimensionType::Volume> meshProperty(this->m_topology);    
  };
  
};
