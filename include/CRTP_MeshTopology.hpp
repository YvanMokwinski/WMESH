#pragma once

#include "CRTP_MeshEntity.hpp"
#include "entity_iterator.hpp"

#include "ReferenceShapeTraits.hpp"
template <typename _int_t> struct RangeMeshEntities
{

protected: _int_t m_start;
protected: _int_t m_bound;
public: using iterator = entity_iterator<_int_t>;

public: inline iterator begin()  noexcept
  {
    return iterator(m_start);
  };
  
public: inline iterator end() noexcept
  {
    return iterator(m_bound);     
  };
  
public: inline RangeMeshEntities(const _int_t start_,
				 const _int_t bound_)
  : m_start(start_),m_bound(bound_)
  {
  };  

};



template <class impl_t> struct CRTP_MeshTopology_traits
{
public: static constexpr bool IsMixed = false;
public: static constexpr DimensionType::enum_t Dimension = DimensionType::Node;
};

template <typename impl_t> class CRTP_MeshTopology
{
  //!
  //! @brief Constant static cast of the derived class.
  //!  
private: inline const impl_t & asImp() const { return static_cast<const impl_t&>(*this); }
  
  //!
  //! @brief Static cast of the derived class.
  //!  
private: inline impl_t & asImp()  { return static_cast<impl_t&>(*this); }    
  
  //
  // Typedefs
  //
public: using cellidx_t = int_t;

  //!
  //! @brief Associated traits of the class.
  //!  
private: using traits_t = CRTP_MeshTopology_traits<impl_t>;

  //
  // Type of the entity kind (tet,hex etc..)
  //
public: using entitykind_t 	= typename traits_t::entitykind_t;

  //!
  //! @brief Type of the mesh entity.
  //!  
public: template <DimensionType::enum_t _dimension>
using entity_t = typename traits_t::template entity_t<_dimension>;
  
  //!
  //! @brief Type of the range of entities.
  //!  
public: template <DimensionType::enum_t _dimension> using range_t = typename traits_t::template range_t<_dimension>;
  
  //!
  //! @brief The dimension as a constant expression.
  //!  
public: static constexpr DimensionType::enum_t Dimension = traits_t::Dimension;

public: using cell_t 	= entity_t<Dimension>;
public: using node_t 	= entity_t<DimensionType::Node>;
  
  //!
  //! @brief Deleted copy constructor.
  //!  
public: CRTP_MeshTopology(const CRTP_MeshTopology<impl_t>&) = delete;

  //!
  //! @brief Protected default constructor.
  //!  
protected: inline CRTP_MeshTopology(){};  
  
  //!
  //! @brief Get the dimension of the topology.
  //! @return The dimension of the topology.
  //!
public: inline DimensionType::enum_t GetDimension() const noexcept
  {
    return asImp().GetDimension();
  };

  //!
  //! @brief Get the number of entities for a specific dimension.
  //! @return The number of entities for a specific dimension.
  //!
public: template <DimensionType::enum_t _dimension> inline int_t GetNumEntities() const noexcept
  {
    return asImp().GetNumEntities<_dimension>();
  };
  
  //!
  //! @brief Get the number of entities for a specific dimension and a specific type.
  //! @return The number of entities for a specific dimension and a specific type.
  //!
public: template <DimensionType::enum_t _dimension>
inline int_t GetNumEntities(const typename ReferenceShapeTraits<_dimension>::type_t::enum_t cellType_) const noexcept
  {
    return asImp().GetNumEntities<_dimension>(cellType_);
  };

  
  //!
  //! @brief Get the index of an entity of a specific dimension.
  //! @param entity_ the mesh entity.
  //! @return The index of the mesh entity.
  //!
public: template <DimensionType::enum_t _dimension>
inline int_t GetEntityIndex(const entity_t<_dimension>& entity_) const noexcept
  {
    return asImp().GetEntityIndex<_dimension>(entity_);
  };
  
  //!
  //! @brief Get the nodes of a specific entity of a cell.
  //! @param target_ the dimension of the target entity. 
  //! @param cell_ the cell 
  //! @param localEntityIndex_ the local entity index.
  //! @param entityToNodes_ the nodes.
  //! @return the number of nodes of the entity.
  //!
public: template <DimensionType::enum_t _target>
inline unsigned int GetEntityToNodes(const cell_t& 	cell_,
				     const unsigned int localEntityIndex_,
				     node_t 		entityToNodes_[]) const noexcept
  {
    return asImp().GetEntityToNodes<_target>(cell_,
					     localEntityIndex_,
					     entityToNodes_);
  };


  //!
  //! @brief Get the set of entities of a given dimension.
  //! @param entityKind_ The kind of the entity.
  //! @return The range.
  //!
public: template <DimensionType::enum_t _dimension>
inline range_t<_dimension> GetEntities(const typename ReferenceShapeTraits<_dimension>::type_t::enum_t entityKind_) const noexcept
  {
    return asImp().GetEntities<_dimension>(entityKind_);
  };

public: template <DimensionType::enum_t _source,
		  DimensionType::enum_t _target>
inline unsigned int GetEntityToEntities(const entity_t<_source>& sourceEntity_,
					entity_t<_target> targetEntities_[]) const noexcept
  {
    return asImp().GetEntityToEntities<_source,_target>(sourceEntity_,
							targetEntities_);
  };
  

#if 0    
public: 
  inline entity_t<Dimension> GetNeighborHalfFacet(const entity_t<Dimension>& 		cell_,
						  const unsigned int 			localFacetIndex_,
						  entity_t<DimensionType::Node> 	entityToNodes_[]) const noexcept
  {
    return asImp().GetNeighborHalfFacet(cell_,
					localEntityIndex_,
					entityToNodes_);
  };
public: template <DimensionType::enum_t _dimension> inline int_t GetEntityTopologicalId(const entity_t<_dimension>& entity_) const noexcept
  {
    return asImp().GetEntityTopologicalId<_dimension>(entity_);
  };
#endif  
  
};


