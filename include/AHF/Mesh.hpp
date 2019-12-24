#pragma once
#include "CRTP_Mesh.hpp"

#include "AHF/MeshTopology3D.hpp"
#include "AHF/MeshGeometry.hpp"
#include "AHF/MeshProperty.hpp"

namespace AHF
{
  template <const DimensionType::enum_t _dimensionType> class Mesh;
};

template <const DimensionType::enum_t _dimensionType> 
struct Traits_CRTP_Mesh< AHF::Mesh<_dimensionType> >
{
public: using meshtopology_t = AHF::MeshTopology<_dimensionType>;
public: using meshgeometry_t = AHF::MeshGeometry;
};

namespace AHF
{
  
  template <const DimensionType::enum_t _dimensionType>
  class Mesh : public CRTP_Mesh< Mesh<_dimensionType> >
  {

  private: using this_t = Mesh;
  private: using traits_t = Traits_CRTP_Mesh<this_t>;
    
  public: using meshtopology_t = typename traits_t::meshtopology_t;
  public: using meshgeometry_t = typename traits_t::meshgeometry_t;
    
  private: meshtopology_t * m_topology;
  private: meshgeometry_t * m_geometry;
    
  public: Mesh(meshtopology_t * topology_,
	       meshgeometry_t * geometry_)
    : m_topology(topology_),
      m_geometry(geometry_)      
    {
      MeshProperty<int,1,DimensionType::Volume> meshProperty(topology_);
    };
    
    //
    // Constructor.
    //
  public: Mesh(const char * filename_)
    {
      
      {
	Input::Medit inputMedit(filename_);
	this->m_topology = meshtopology_t::Build(inputMedit);      
	this->m_geometry = meshgeometry_t::Build(this->m_topology,
						 inputMedit);
      }
      
    };
    
  public: const meshtopology_t * GetTopology() const noexcept
    {
      return this->m_topology;
    };
    
  public: const meshgeometry_t * GetGeometry() const noexcept
    {
      return this->m_geometry;
    };
    
  };

  using Mesh1D = Mesh<DimensionType::Edge>;  
  using Mesh2D = Mesh<DimensionType::Face>;  
  using Mesh3D = Mesh<DimensionType::Volume>;  

};
