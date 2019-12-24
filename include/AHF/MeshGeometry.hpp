#pragma once
#include "Point.hpp"
#include <valarray>
#include "CRTP_MeshGeometry.hpp"

namespace AHF
{
  class MeshGeometry;  
};

//!
//! @brief The traits.
//!    
template <> struct Traits_CRTP_MeshGeometry<AHF::MeshGeometry>
{
  
  //!
  //! @brief Define the type of the point.
  //!    
public: using pts_t = Point3d<double>;
  //!
  //! @brief Define the mesh topology.
  //!    
public: using meshtopology_t = AHF::MeshTopology3D;  
};

namespace AHF
{
  class MeshGeometry : public CRTP_MeshGeometry<MeshGeometry>
  {
    //!
    //! @brief This type.
    //!    
  private: using this_t = MeshGeometry;
    //!
    //! @brief Typedefs
    //!
  private: using traits_t = Traits_CRTP_MeshGeometry<this_t>;
    
    //!
    //! @brief Type of the mesh topology.
    //!
  public: using meshtopology_t 	= typename traits_t::meshtopology_t;
    
    //!
    //! @brief Type of the node entity of the topology.
    //!
  public: using node_t 	= meshtopology_t::entity_t<DimensionType::Node>;
    
    //!
    //! @brief Type of the point.
    //!
  public: using pts_t 	= typename traits_t::pts_t;

  public: static inline MeshGeometry*Build(const meshtopology_t*meshtopology_,
					   Input::Medit&inputMedit_) noexcept
    {
      MeshGeometry * meshGeometry = new MeshGeometry(meshtopology_);
      inputMedit_.ReadGeometry(meshGeometry->m_pts);
      return meshGeometry;
    };

    //!
    //! @brief Constructor.
    //! @param meshtopology_ The mesh topology the geometry is defined on.
    //!
  public: inline MeshGeometry(const meshtopology_t*meshtopology_) noexcept;
    
  public: inline unsigned int 		GetDimension() 			const noexcept;
  public: inline const meshtopology_t* 	GetMeshTopology() 		const noexcept;
  public: inline int_t 			GetNumPoints() 			const noexcept;
  public: inline pts_t 			GetPoint(const node_t node_) 	const noexcept;
    
  public: template<typename _int_t>
  inline pts_t& operator[] (const _int_t nodeIndex_) noexcept;
    
  public: template<typename _int_t>
  inline const pts_t& operator[] (const _int_t nodeIndex_) const noexcept;
    //!
    //! @brief The mesh topology the geometry is defined on.
    //!
  private: const meshtopology_t * m_meshtopology;
    //!
    //! @brief The array of points.
    //!
  private: std::valarray<pts_t> m_pts;    
  };
  
  //
  //
  //
  inline MeshGeometry::MeshGeometry(const MeshGeometry::meshtopology_t*meshtopology_) noexcept
    : m_meshtopology(meshtopology_),      
      m_pts(meshtopology_->GetNumEntities<DimensionType::Node>())
  {
  };
  
  //
  //
  //
  inline unsigned int MeshGeometry::GetDimension() const noexcept
  {
    return 3;
  };
  
  //
  //
  //
  inline const MeshGeometry::meshtopology_t* MeshGeometry::GetMeshTopology() const noexcept 
  {
    return this->m_meshtopology;
  };
  
  //
  //
  //
  inline int_t MeshGeometry::GetNumPoints() const noexcept 
  {
    return this->m_meshtopology->GetNumEntities<DimensionType::Node>();
  };

  //
  //
  //
  inline MeshGeometry::pts_t MeshGeometry::GetPoint(const node_t node_) const noexcept
  {
    return this->m_pts[this->m_meshtopology->GetEntityIndex<DimensionType::Node>(node_)];
  };

  //
  //
  //
  template<typename _int_t> inline MeshGeometry::pts_t& MeshGeometry::operator[] (const _int_t nodeIndex_) noexcept
  {
    return this->m_pts[nodeIndex_];
  };

  //
  //
  //
  template<typename _int_t> inline const MeshGeometry::pts_t& MeshGeometry::operator[] (const _int_t nodeIndex_) const noexcept
  {
    return this->m_pts[nodeIndex_];
  };
  
};
