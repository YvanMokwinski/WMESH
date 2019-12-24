#pragma once
#include "CRTP.hpp"
//!
//! @brief Traits of the class CRTP_MeshGeometry.
//!  
template <typename _impl_t> struct Traits_CRTP_MeshGeometry;

//!
//! @brief Define a Curiously Recursive Template Pattern for the mesh geometry.
//!
template <typename _impl_t> class CRTP_MeshGeometry : public CRTP<_impl_t>
{  
  //!
  //! @brief Type of the traits of the derived class.
  //!  
private: 	using traits_t	= Traits_CRTP_MeshGeometry<_impl_t>;
  
  //!
  //! @brief Type of the topology.
  //!  
public:   	using meshtopology_t 	= typename traits_t::meshtopology_t;

  //!
  //! @brief Type of the node entity of the topology.
  //!  
public:   	using node_t 	= typename traits_t::meshtopology_t::template entity_t<DimensionType::Node>;
  
  //!
  //! @brief Type of the point of the geometry.
  //!  
public:   	using pts_t 	= typename traits_t::pts_t;
  
  //!
  //! @brief Deleted copy constructor.
  //!  
public: CRTP_MeshGeometry(const CRTP_MeshGeometry<_impl_t>&) = delete;
  
  //!
  //! @brief Protected default constructor.
  //!  
protected: inline CRTP_MeshGeometry();
  
  //!
  //! @brief Get the topology.  
  //!
public: inline const meshtopology_t* 	GetMeshTopology() 		const noexcept;

  //!
  //! @brief Get the dimension of the geometry.
  //! @return The dimension of the geometry.
  //!
public: inline unsigned int 		GetDimension() 			const noexcept;

  //!
  //! @brief Get the number of points.
  //! @return The number of points.
  //!  
public: inline int_t 			GetNumPoints() 			const noexcept;

  //!
  //! @brief Get a point.
  //! @param node_ The node entity from the topology.
  //! @return The point.
  //!    
public: inline pts_t 			GetPoint(const node_t node_) 	const noexcept;
};

//
//
//
template <typename _impl_t> 
inline CRTP_MeshGeometry<_impl_t>::CRTP_MeshGeometry() 
{
};

//
//
//
template <typename _impl_t> 
inline const typename CRTP_MeshGeometry<_impl_t>::meshtopology_t* CRTP_MeshGeometry<_impl_t>::GetMeshTopology() const noexcept 
{
  return this->asImp().GetMeshTopology();
};

//
//
//
template <typename _impl_t> 
inline typename CRTP_MeshGeometry<_impl_t>::pts_t
CRTP_MeshGeometry<_impl_t>::GetPoint(const CRTP_MeshGeometry<_impl_t>::node_t node_) const noexcept
{
  return this->asImp().GetPoint(node_);
};

//
//
//
template <typename _impl_t> 
inline unsigned int CRTP_MeshGeometry<_impl_t>::GetDimension() const noexcept 
{
  return this->asImp().GetDimension();
};

//
//
//
template <typename _impl_t> 
inline int_t CRTP_MeshGeometry<_impl_t>::GetNumPoints() const noexcept 
{
  return this->asImp().GetNumPoints();
};

